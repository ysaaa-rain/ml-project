"""DNA foundation-model embedding interfaces and explicit smoke-test backend."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from typing import Protocol, Sequence

import numpy as np


class DNAEmbedder(Protocol):
    """Minimal interface shared by real and smoke-test embedding backends."""

    @property
    def metadata(self) -> dict:
        ...

    def embed(self, sequences: Sequence[str]) -> np.ndarray:
        ...


@dataclass
class HashKmerSmokeEmbedder:
    """Deterministic test backend; never use its vectors as scientific evidence."""

    dimension: int = 32
    kmer_size: int = 4
    seed: int = 20260913

    def __post_init__(self) -> None:
        if self.dimension < 2:
            raise ValueError("dimension must be at least 2")
        if self.kmer_size < 1:
            raise ValueError("kmer_size must be at least 1")

    @property
    def metadata(self) -> dict:
        return {
            "backend": "hash-kmer-smoke",
            "scientific_use": False,
            "warning": "Deterministic code-path test only; not a DNA foundation model.",
            "dimension": self.dimension,
            "kmer_size": self.kmer_size,
            "seed": self.seed,
        }

    def embed(self, sequences: Sequence[str]) -> np.ndarray:
        matrix = np.zeros((len(sequences), self.dimension), dtype=np.float32)
        for row_index, sequence in enumerate(sequences):
            sequence = sequence.upper()
            kmers = (
                [sequence]
                if len(sequence) < self.kmer_size
                else [
                    sequence[start : start + self.kmer_size]
                    for start in range(len(sequence) - self.kmer_size + 1)
                ]
            )
            for kmer in kmers:
                payload = f"{self.seed}:{kmer}".encode("ascii")
                digest = hashlib.sha256(payload).digest()
                bucket = int.from_bytes(digest[:4], "little") % self.dimension
                sign = 1.0 if digest[4] % 2 == 0 else -1.0
                matrix[row_index, bucket] += sign
            norm = float(np.linalg.norm(matrix[row_index]))
            if norm:
                matrix[row_index] /= norm
        return matrix


class HuggingFaceDNAEmbedder:
    """Frozen Hugging Face DNA model with attention-aware mean/CLS pooling."""

    def __init__(
        self,
        *,
        model_id: str,
        revision: str,
        model_class: str = "masked_lm",
        pooling: str = "mean",
        max_length: int = 512,
        batch_size: int = 8,
        device: str = "auto",
        trust_remote_code: bool = False,
        seed: int = 20260913,
    ) -> None:
        try:
            import torch
            from transformers import AutoConfig, AutoModel, AutoModelForMaskedLM, AutoTokenizer
            import transformers.modeling_utils as modeling_utils
            import transformers.pytorch_utils as pytorch_utils
            from transformers.modeling_utils import PreTrainedModel
        except ImportError as exc:
            raise RuntimeError(
                "The Hugging Face backend requires torch and transformers. "
                "Install the audited embedding environment before running L0."
            ) from exc

        # Some pinned DNA-model revisions still import these helpers from
        # modeling_utils, while Transformers 5 moved them to pytorch_utils.
        # Patch the public module namespace in memory; never modify Hub cache.
        if not hasattr(modeling_utils, "prune_linear_layer"):
            modeling_utils.prune_linear_layer = pytorch_utils.prune_linear_layer
        if not hasattr(modeling_utils, "find_pruneable_heads_and_indices"):
            def find_pruneable_heads_and_indices(heads, n_heads, head_size, already_pruned_heads):
                heads = set(heads) - already_pruned_heads
                mask = torch.ones(n_heads, head_size)
                for head in heads:
                    head -= sum(1 for pruned in already_pruned_heads if pruned < head)
                    mask[head] = 0
                mask = mask.view(-1).contiguous().eq(1)
                index = torch.arange(len(mask))[mask].long()
                return heads, index

            modeling_utils.find_pruneable_heads_and_indices = find_pruneable_heads_and_indices
        if not hasattr(PreTrainedModel, "all_tied_weights_keys"):
            @property
            def all_tied_weights_keys(self):
                cached = getattr(self, "_compat_all_tied_weights_keys", None)
                if cached is not None:
                    return cached
                keys = getattr(self, "_tied_weights_keys", None) or []
                return {key: key for key in keys}

            @all_tied_weights_keys.setter
            def all_tied_weights_keys(self, value):
                self._compat_all_tied_weights_keys = value

            PreTrainedModel.all_tied_weights_keys = all_tied_weights_keys
        if not hasattr(PreTrainedModel, "get_head_mask"):
            def get_head_mask(self, head_mask, num_hidden_layers, is_attention_chunked=False):
                if head_mask is None:
                    return [None] * num_hidden_layers
                if head_mask.dim() == 1:
                    head_mask = head_mask.unsqueeze(0).unsqueeze(0).unsqueeze(-1).unsqueeze(-1)
                    head_mask = head_mask.expand(num_hidden_layers, -1, -1, -1, -1)
                elif head_mask.dim() == 2:
                    head_mask = head_mask.unsqueeze(1).unsqueeze(-1).unsqueeze(-1)
                if head_mask.dim() != 5:
                    raise ValueError("head_mask must have dimension 1, 2 or 5")
                head_mask = head_mask.to(dtype=self.dtype)
                if is_attention_chunked:
                    head_mask = head_mask.unsqueeze(-1)
                return head_mask

            PreTrainedModel.get_head_mask = get_head_mask

        if not revision or revision == "main":
            raise ValueError("a fixed model revision is required for reproducible experiments")
        if pooling not in {"mean", "cls"}:
            raise ValueError("pooling must be 'mean' or 'cls'")
        if model_class not in {"auto", "masked_lm"}:
            raise ValueError("model_class must be 'auto' or 'masked_lm'")
        if max_length < 8 or batch_size < 1:
            raise ValueError("max_length must be >= 8 and batch_size must be >= 1")

        self._torch = torch
        self.model_id = model_id
        self.revision = revision
        self.model_class = model_class
        self.pooling = pooling
        self.max_length = max_length
        self.batch_size = batch_size
        self.trust_remote_code = trust_remote_code
        self.seed = seed
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
        torch.use_deterministic_algorithms(True, warn_only=True)
        self.device = (
            "cuda" if device == "auto" and torch.cuda.is_available() else "cpu" if device == "auto" else device
        )
        self.tokenizer = AutoTokenizer.from_pretrained(
            model_id,
            revision=revision,
            trust_remote_code=trust_remote_code,
        )
        model_config = AutoConfig.from_pretrained(
            model_id,
            revision=revision,
            trust_remote_code=trust_remote_code,
        )
        compatibility_defaults = {
            "is_decoder": False,
            "add_cross_attention": False,
            "chunk_size_feed_forward": 0,
            "output_attentions": False,
            "output_hidden_states": True,
            "use_return_dict": True,
        }
        for name, value in compatibility_defaults.items():
            if not hasattr(model_config, name):
                setattr(model_config, name, value)
        loader = AutoModelForMaskedLM if model_class == "masked_lm" else AutoModel
        self.model = loader.from_pretrained(
            model_id,
            revision=revision,
            trust_remote_code=trust_remote_code,
            config=model_config,
            use_safetensors=False,
        )
        self.model.to(self.device)
        self.model.eval()

    @property
    def metadata(self) -> dict:
        return {
            "backend": "huggingface",
            "scientific_use": True,
            "model_id": self.model_id,
            "revision": self.revision,
            "model_class": self.model_class,
            "pooling": self.pooling,
            "max_length": self.max_length,
            "batch_size": self.batch_size,
            "device": self.device,
            "trust_remote_code": self.trust_remote_code,
            "seed": self.seed,
        }

    @staticmethod
    def _hidden_states(outputs):
        hidden_states = getattr(outputs, "hidden_states", None)
        if hidden_states:
            return hidden_states[-1]
        last_hidden_state = getattr(outputs, "last_hidden_state", None)
        if last_hidden_state is not None:
            return last_hidden_state
        if isinstance(outputs, tuple) and outputs:
            return outputs[0]
        raise RuntimeError("model output does not expose hidden states")

    def embed(self, sequences: Sequence[str]) -> np.ndarray:
        if not sequences:
            return np.empty((0, 0), dtype=np.float32)
        batches: list[np.ndarray] = []
        torch = self._torch
        for start in range(0, len(sequences), self.batch_size):
            batch = list(sequences[start : start + self.batch_size])
            encoded = self.tokenizer(
                batch,
                padding=True,
                truncation=True,
                max_length=self.max_length,
                return_tensors="pt",
                return_special_tokens_mask=True,
            )
            special_tokens_mask = encoded.pop("special_tokens_mask", None)
            encoded = {name: tensor.to(self.device) for name, tensor in encoded.items()}
            with torch.no_grad():
                outputs = self.model(**encoded, output_hidden_states=True, return_dict=True)
            hidden = self._hidden_states(outputs)
            if self.pooling == "cls":
                pooled = hidden[:, 0, :]
            else:
                mask = encoded["attention_mask"].bool()
                if special_tokens_mask is not None and special_tokens_mask.shape == mask.shape:
                    mask = mask & ~special_tokens_mask.to(self.device).bool()
                if mask.shape[1] != hidden.shape[1]:
                    mask = mask[:, : hidden.shape[1]]
                mask_float = mask.unsqueeze(-1).to(hidden.dtype)
                denominator = mask_float.sum(dim=1).clamp(min=1.0)
                pooled = (hidden * mask_float).sum(dim=1) / denominator
            batches.append(pooled.detach().cpu().to(torch.float32).numpy())
        return np.concatenate(batches, axis=0)


def build_embedder(config: dict) -> DNAEmbedder:
    """Create an embedder while keeping the smoke backend explicitly gated."""

    backend = str(config.get("backend", "huggingface"))
    if backend == "hash-smoke":
        if not config.get("allow_test_backend", False):
            raise ValueError("hash-smoke requires allow_test_backend: true")
        return HashKmerSmokeEmbedder(
            dimension=int(config.get("dimension", 32)),
            kmer_size=int(config.get("kmer_size", 4)),
            seed=int(config.get("seed", 20260913)),
        )
    if backend != "huggingface":
        raise ValueError(f"unsupported embedding backend: {backend}")
    return HuggingFaceDNAEmbedder(
        model_id=str(config["model_id"]),
        revision=str(config["revision"]),
        model_class=str(config.get("model_class", "masked_lm")),
        pooling=str(config.get("pooling", "mean")),
        max_length=int(config.get("max_length", 512)),
        batch_size=int(config.get("batch_size", 8)),
        device=str(config.get("device", "auto")),
        trust_remote_code=bool(config.get("trust_remote_code", False)),
        seed=int(config.get("seed", 20260913)),
    )

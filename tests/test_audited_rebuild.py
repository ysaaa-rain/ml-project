import pandas as pd
from preprocessing.rebuild_audited_data import related, tss_window, dbtbs_records, rc
from preprocessing.windowing import apply_tss_window
from preprocessing.source_adapters import adapt_regulondb_gff3
from preprocessing.pipeline import normalize_metadata

def test_already_oriented_negative_strand_is_not_flipped_again(tmp_path):
    s='A'*60+'C'+'G'*20
    p=tmp_path/'x.gff3';p.write_text('NC_000913.3\tR\tpromoter\t100\t100\t.\t-\t.\tname=p;Sequence='+s[:60].lower()+s[60]+s[61:].lower()+';SigmaFactor=Sigma70;Confidence=Strong\n')
    df,_=adapt_regulondb_gff3(p,species='E.coli');out,_=apply_tss_window(df)
    assert out.iloc[0].sequence==s
    assert out.iloc[0].input_orientation=='transcription_forward'

def test_genomic_negative_window_uses_asymmetric_bounds_correctly():
    g='ACGT'*70;t=130
    df=pd.DataFrame([{'sequence_id':'x','sequence':g,'tss_position':t,'strand':'-','input_orientation':'genomic_forward'}])
    out,_=apply_tss_window(df)
    assert out.iloc[0].sequence==tss_window(g,t,'-')
    assert out.iloc[0].sequence[60]==rc(g[t-1])
    assert len(out.iloc[0].sequence)==81

def test_tss_at_sixty_cannot_supply_sixty_upstream_bases():
    df=pd.DataFrame([{'sequence_id':'x','sequence':'A'*100,'tss_position':60,'strand':'+'}])
    out,report=apply_tss_window(df)
    assert out.empty and report['dropped_by_reason']['window_out_of_range']==1

def test_sigma_dedup_retains_both_labels():
    df=pd.DataFrame([{'sequence_id':str(i),'species':'E','source_dataset':'s','sequence':'ACGT'*21,'sigma_factor_type':z} for i,z in enumerate(['Sigma70','Sigma38'])])
    out,_=normalize_metadata(df)
    assert len(out)==1 and set(__import__('json').loads(out.iloc[0].sigma_labels))=={'Sigma70','Sigma38'}
    assert out.iloc[0].sigma_label_status == 'multi_sigma'

def test_edit_and_reverse_complement_leakage():
    s='ACGT'*20+'A';changed=list(s)
    for i in range(9,81,10):changed[i]={'A':'C','C':'G','G':'T','T':'A'}[changed[i]]
    assert related(s,''.join(changed))
    assert related(s,rc(s))
    assert related(s,s[3:]+'AAA')
    assert not related('A'*81,'C'*81)

def test_dbtbs_nd_cannot_be_used_as_tss(tmp_path):
    p=tmp_path/'p.html';p.write_text('<table><tr><th>Genes</th><th>Synonyms</th><th>Direction</th><th>Genome position</th></tr><tr><td>x</td><td></td><td>+</td><td>100..200</td></tr></table><table><tr><th>Binding factor</th><th>Regulation</th><th>Location</th><th>Absolute position</th><th>Binding seq.</th><th>Evidence</th></tr><tr><td>SigA</td><td>Promoter</td><td>ND</td><td>100..103</td><td>AAAA</td><td>PE</td></tr></table>')
    rows=dbtbs_records(p,'A'*1000)
    assert rows[0]['reason']=='Location_ND_or_unresolved' and 'tss' not in rows[0]

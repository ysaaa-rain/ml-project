$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '../../..')).Path
$cacheRoot = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '../cache')).Path
$cleanupPlan = Get-Content -LiteralPath (Join-Path $PSScriptRoot 'cleanup_plan.json') -Raw -Encoding utf8 | ConvertFrom-Json
if ($cleanupPlan.file_count -ne 84) { throw 'Unexpected cleanup scope' }
# Validate every target before the first deletion. No recursive directory removal.
foreach ($entry in $cleanupPlan.files) {
    $targetPath = (Resolve-Path -LiteralPath $entry.absolute_path).Path
    if (-not $targetPath.StartsWith($cacheRoot + [IO.Path]::DirectorySeparatorChar, [StringComparison]::OrdinalIgnoreCase)) { throw 'Target escapes historical cache directory' }
    if ([IO.Path]::GetExtension($targetPath) -ne '.npy') { throw 'Unexpected file type' }
    if ((Get-Item -LiteralPath $targetPath).PSIsContainer) { throw 'Expected a file' }
    $ancestor = Get-Item -LiteralPath $targetPath
    while ($null -ne $ancestor -and $ancestor.FullName.Length -ge $cacheRoot.Length) {
        if ($ancestor.Attributes -band [IO.FileAttributes]::ReparsePoint) { throw 'Reparse-point targets are outside this cleanup scope' }
        $ancestor = if ($ancestor.PSIsContainer) { $ancestor.Parent } else { $ancestor.Directory }
    }
    $targetHash = (Get-FileHash -LiteralPath $targetPath -Algorithm SHA256).Hash.ToLowerInvariant()
    if ($targetHash -ne $entry.sha256) { throw 'Target changed since inventory' }
}
foreach ($entry in $cleanupPlan.files) {
    Remove-Item -LiteralPath $entry.absolute_path -Force
}
$absentCount = @($cleanupPlan.files | Where-Object { -not (Test-Path -LiteralPath $_.absolute_path) }).Count
if ($absentCount -ne $cleanupPlan.file_count) { throw 'Some inventoried caches remain' }
foreach ($entry in $cleanupPlan.files) {
    if (-not (Test-Path -LiteralPath (Join-Path $repoRoot $entry.manifest))) { throw 'Original manifest missing' }
}
$result = [ordered]@{
    status = 'COMPLETED'
    executed_utc = [DateTime]::UtcNow.ToString('o')
    plan_sha256 = (Get-FileHash -LiteralPath (Join-Path $PSScriptRoot 'cleanup_plan.json') -Algorithm SHA256).Hash.ToLowerInvariant()
    deleted_files = $absentCount
    reclaimed_bytes = $cleanupPlan.total_bytes
    retained = 'source DNA, weights, index.tsv, cache_manifest.json, all per-run evidence, negative results, logs and user work'
    availability = 'inventoried historical cache arrays intentionally evicted; regenerate in a new directory before reuse'
}
$result | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $PSScriptRoot 'cleanup_result.json') -Encoding utf8
$result | ConvertTo-Json -Compress

param(
    [ValidateSet('3d', 'seq')][string]$Variant = '3d',
    [ValidateRange(60, 7200)][int]$MaximumSeconds = 3600,
    [switch]$StartDownload
)

# Default is validation/preview only. A human/root-agent must explicitly choose
# -StartDownload; this file is not an automatic background downloader.
$ErrorActionPreference = 'Stop'
$weightDirectory = $PSScriptRoot
$expectedInfo = @{
    '3d' = @{ Length = [int64]4843266070; ETag = '"120ae6016-62f0c141fb800"'; Name = 'flexpert_3d_weights.bin' }
    'seq' = @{ Length = [int64]4842603885; ETag = '"120a4456d-62f0c13872180"'; Name = 'flexpert_seq_weights.verified.part.bin' }
}[$Variant]
$weightUrl = "https://data.ciirc.cvut.cz/public/projects/2025Flexpert/flexpert-weights/flexpert_${Variant}_weights.bin"
$headerFile = Join-Path $weightDirectory "resume_${Variant}_latest_head.headers.txt"
$headerArguments = @('--head', '--location', '--silent', '--show-error', '--fail', '--connect-timeout', '8', '--max-time', '20', '--dump-header', $headerFile, '--output', 'NUL', $weightUrl)
& curl.exe @headerArguments
if ($LASTEXITCODE -ne 0) { throw 'Official HEAD request failed; original partial files have not been modified.' }
$headerText = Get-Content -LiteralPath $headerFile -Raw
$lengthMatches = [regex]::Matches($headerText, '(?im)^Content-Length:\s*(\d+)\s*$')
$etagMatches = [regex]::Matches($headerText, '(?im)^ETag:\s*(.+?)\s*$')
if ($lengthMatches.Count -eq 0 -or $etagMatches.Count -eq 0) { throw 'Official HEAD lacks length/ETag; stop for audit.' }
$currentLength = [int64]$lengthMatches[$lengthMatches.Count - 1].Groups[1].Value
$currentETag = $etagMatches[$etagMatches.Count - 1].Groups[1].Value.Trim()
if ($currentLength -ne $expectedInfo.Length -or $currentETag -cne $expectedInfo.ETag) {
    throw 'Official remote fingerprint changed. Do not resume old pieces into a changed entity.'
}
$targetPath = Join-Path $weightDirectory $expectedInfo.Name
if ($Variant -eq '3d' -and -not (Test-Path -LiteralPath ($targetPath + '.aria2'))) {
    throw 'The original 3D aria2 control file is required. Never infer a resume offset from the logical .bin length.'
}
if ($Variant -eq 'seq' -and (Test-Path -LiteralPath $targetPath) -and -not (Test-Path -LiteralPath ($targetPath + '.aria2'))) {
    throw 'An unmanaged Seq target already exists at the distinct verified.part path. Stop rather than overwrite or assume a valid prefix.'
}
$ariaArguments = @(
    '--continue=true', '--always-resume=true', '--max-connection-per-server=1', '--split=1',
    '--file-allocation=none', '--auto-file-renaming=false', '--allow-overwrite=false',
    '--check-certificate=true', '--connect-timeout=10', '--timeout=30', '--max-tries=3',
    '--retry-wait=30', '--summary-interval=30', '--disk-cache=8M',
    "--stop=$MaximumSeconds", "--header=If-Match: $($expectedInfo.ETag)",
    "--dir=$weightDirectory", "--out=$($expectedInfo.Name)",
    "--log=$(Join-Path $weightDirectory ('resume_' + $Variant + '.log'))", $weightUrl
)
[pscustomobject]@{
    RemoteLengthBytes = $currentLength
    RemoteETag = $currentETag
    DownloadTarget = $targetPath
    MaximumRuntimeSeconds = $MaximumSeconds
    DiskFreeBytes = (Get-PSDrive -Name ([IO.Path]::GetPathRoot($weightDirectory).Substring(0,1))).Free
    Command = 'aria2c.exe'
    Arguments = $ariaArguments
    StartDownloadRequested = $StartDownload.IsPresent
} | Format-List
if (-not $StartDownload) {
    Write-Output 'Validation/command preview only; no full download was started.'
    return
}
& aria2c.exe @ariaArguments
if ($LASTEXITCODE -ne 0) { throw "aria2 stopped with exit $LASTEXITCODE; preserve both target and .aria2 map for the next bounded retry." }
if (-not (Test-Path -LiteralPath $targetPath) -or (Get-Item -LiteralPath $targetPath).Length -ne $expectedInfo.Length -or (Test-Path -LiteralPath ($targetPath + '.aria2'))) {
    throw 'Incomplete transfer remains. Do not load the model or rename the weight as complete.'
}
Get-FileHash -LiteralPath $targetPath -Algorithm SHA256
Write-Output 'Transfer length/control-file checks passed. Independently validate full ZIP CRC/inventory and retain the full SHA256 before any model load. ETag is not an official SHA256.'

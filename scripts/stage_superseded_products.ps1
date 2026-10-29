param([switch]$Apply)

# Move only. The source-relative layout makes every item recoverable.
$ErrorActionPreference = 'Stop'
$taskRoot = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$holdingRoot = [IO.Path]::GetFullPath((Join-Path $taskRoot '待删除'))
$manifestPath = Join-Path $holdingRoot '移动清单.json'
$previousItems = @{}
if (Test-Path -LiteralPath $manifestPath) {
    $previousManifest = Get-Content -LiteralPath $manifestPath -Raw -Encoding UTF8 | ConvertFrom-Json
    if ($previousManifest.project -ne $taskRoot -or $previousManifest.holding_folder -ne $holdingRoot) { throw 'Existing staging manifest belongs to another folder' }
    foreach ($previousItem in $previousManifest.items) { $previousItems[$previousItem.source] = $previousItem }
}
$plan = Get-Content -LiteralPath (Join-Path $taskRoot 'docs/cleanup_plan_20261023.json') -Raw -Encoding UTF8 | ConvertFrom-Json
if ($taskRoot -ne [IO.Path]::GetFullPath($plan.root)) { throw 'Plan belongs to another project' }
if (-not $holdingRoot.StartsWith($taskRoot + '\', [StringComparison]::OrdinalIgnoreCase)) { throw 'Unsafe holding folder' }

function Get-SharedSha256([string]$Path) {
    $stream = [IO.File]::Open($Path, [IO.FileMode]::Open, [IO.FileAccess]::Read,
        [IO.FileShare]::ReadWrite -bor [IO.FileShare]::Delete)
    $hasher = [Security.Cryptography.SHA256]::Create()
    try { return [Convert]::ToHexString($hasher.ComputeHash($stream)).ToLowerInvariant() }
    finally { $hasher.Dispose(); $stream.Dispose() }
}

$finalPath = [IO.Path]::GetFullPath($plan.final_document)
if ((Get-SharedSha256 $finalPath) -ne $plan.final_sha256) { throw 'Final Word has changed; review the plan again' }
$protected = @($finalPath, (Join-Path $taskRoot 'paper/figures_v58'),
    (Join-Path $taskRoot 'experiments/core_corrected_pipeline'),
    (Join-Path $taskRoot 'experiments/_case_validation_sources'),
    (Join-Path $taskRoot 'experiments/_legacy_paper_data'))
$items = [Collections.Generic.List[object]]::new()

# Complete the path and link checks before moving the first item.
foreach ($entry in $plan.delete) {
    $source = [IO.Path]::GetFullPath($entry.path)
    if (-not $source.StartsWith($taskRoot + '\', [StringComparison]::OrdinalIgnoreCase)) { throw ('Outside project: ' + $source) }
    $relative = $source.Substring($taskRoot.Length + 1)
    if (-not ($relative.StartsWith('paper\') -or $relative.StartsWith('experiments\_sources\'))) { throw ('Unexpected scope: ' + $relative) }
    foreach ($keep in $protected) {
        if ($source -eq $keep -or $source.StartsWith($keep + '\', [StringComparison]::OrdinalIgnoreCase) -or $keep.StartsWith($source + '\', [StringComparison]::OrdinalIgnoreCase)) {
            throw ('Protected source: ' + $source)
        }
    }
    $destination = [IO.Path]::GetFullPath((Join-Path $holdingRoot $relative))
    if (-not $destination.StartsWith($holdingRoot + '\', [StringComparison]::OrdinalIgnoreCase)) { throw 'Unsafe destination' }
    $exists = Test-Path -LiteralPath $source
    if (-not $exists -and (Test-Path -LiteralPath $destination)) {
        if (-not $previousItems.ContainsKey($source) -or $previousItems[$source].status -ne 'moved' -or $previousItems[$source].destination -ne $destination) {
            throw ('Existing destination has no matching successful move record: ' + $destination)
        }
        $items.Add($previousItems[$source])
        continue
    }
    if (-not $exists) { throw ('Missing planned source: ' + $source) }
    if (Test-Path -LiteralPath $destination) { throw ('Destination already exists; no overwrite: ' + $destination) }
    $cursor = Get-Item -LiteralPath $source -Force
    while ($null -ne $cursor -and $cursor.FullName -ne $taskRoot) {
        if ($cursor.Attributes -band [IO.FileAttributes]::ReparsePoint) { throw ('Linked path protected: ' + $source) }
        $cursor = if ($cursor.PSIsContainer) { $cursor.Parent } else { $cursor.Directory }
    }
    if ($entry.directory) {
        $descendants = @(Get-ChildItem -LiteralPath $source -Recurse -Force)
        if ($descendants | Where-Object { $_.Attributes -band [IO.FileAttributes]::ReparsePoint }) { throw ('Directory contains a link: ' + $source) }
        $files = @($descendants | Where-Object { -not $_.PSIsContainer })
        $bytes = ($files | Measure-Object -Property Length -Sum).Sum
    } else {
        $files = @(Get-Item -LiteralPath $source -Force)
        $bytes = $files[0].Length
    }
    $items.Add([pscustomobject]@{source=$source; destination=$destination; directory=[bool]$entry.directory;
        reason=$entry.reason; files=$files.Count; bytes=$bytes; status='pending'; error=$null})
}
if (-not $Apply) {
    [pscustomobject]@{destination=$holdingRoot; targets=$items.Count;
        files=($items | Measure-Object -Property files -Sum).Sum;
        bytes=($items | Measure-Object -Property bytes -Sum).Sum} | Format-List
    exit 0
}

[IO.Directory]::CreateDirectory($holdingRoot) | Out-Null
function Save-StagingManifest {
    [ordered]@{operation='move_only'; permanent_deletion=$false; project=$taskRoot; holding_folder=$holdingRoot;
        final_document=$finalPath; final_sha256=$plan.final_sha256; items=@($items.ToArray())} |
        ConvertTo-Json -Depth 7 | Set-Content -LiteralPath $manifestPath -Encoding UTF8
}
Save-StagingManifest
$counter = 0
foreach ($item in $items) {
    if ($item.status -eq 'moved') { continue }
    try {
        [IO.Directory]::CreateDirectory([IO.Path]::GetDirectoryName($item.destination)) | Out-Null
        if ($item.directory) { [IO.Directory]::Move($item.source, $item.destination) }
        else { [IO.File]::Move($item.source, $item.destination) }
        if ((Test-Path -LiteralPath $item.source) -or -not (Test-Path -LiteralPath $item.destination)) { throw 'Move verification failed' }
        $item.status = 'moved'
        $item.error = $null
    } catch {
        $item.status = 'failed'
        $item.error = $_.Exception.Message
    }
    $counter++
    if ($counter % 100 -eq 0) { Save-StagingManifest; Write-Output ('Moved/checkpointed ' + $counter + '/' + $items.Count + ' targets') }
}
Save-StagingManifest
$moved = @($items | Where-Object status -eq 'moved')
$failed = @($items | Where-Object status -eq 'failed')
$wordUnchanged = (Get-SharedSha256 $finalPath) -eq $plan.final_sha256
$statusPath = Join-Path $taskRoot 'docs/cleanup_status_20261023.json'
$status = Get-Content -LiteralPath $statusPath -Raw -Encoding UTF8 | ConvertFrom-Json
$status.deletion_status = 'not_requested_move_only'
$status.deleted_files = 0
$status.deleted_bytes = 0
$status.final_unchanged = $wordUnchanged
$status | Add-Member -Force -NotePropertyName staging -NotePropertyValue ([ordered]@{
    status=$(if ($failed.Count -eq 0) {'complete'} else {'incomplete'}); folder=$holdingRoot;
    manifest=$manifestPath; moved_targets=$moved.Count; moved_files=($moved | Measure-Object -Property files -Sum).Sum;
    moved_bytes=($moved | Measure-Object -Property bytes -Sum).Sum; failed_targets=$failed.Count})
$status | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $statusPath -Encoding UTF8
if (-not $wordUnchanged) { throw 'Final Word verification failed' }
Write-Output ('Moved ' + $moved.Count + ' targets; failed ' + $failed.Count + '; deleted 0 files.')
if ($failed.Count) { $failed | Select-Object source,error | Format-List; exit 1 }

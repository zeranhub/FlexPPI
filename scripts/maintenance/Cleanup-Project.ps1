# One-time, allowlisted housekeeping for the September 2026 v2 transition.
# Preview by default. Use -Apply only after reviewing the printed targets.
[CmdletBinding()]
param([switch]$Apply)

$ErrorActionPreference = 'Stop'
$projectRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '../..')).TrimEnd('\')
$sourcePrefix = $projectRoot + '\'
$archiveParent = Join-Path (Split-Path $projectRoot -Parent) 'FlexPPI-archive'
$archiveRoot = Join-Path $archiveParent ('cleanup-' + (Get-Date -Format 'yyyyMMdd-HHmmss-fff'))
$archivePrefix = [IO.Path]::GetFullPath($archiveRoot).TrimEnd('\') + '\'
if ($archivePrefix.StartsWith($sourcePrefix, [StringComparison]::OrdinalIgnoreCase)) {
    throw 'Archive must be outside the project.'
}
if (Test-Path -LiteralPath $archiveParent) {
    if ((Get-Item -LiteralPath $archiveParent -Force).Attributes -band [IO.FileAttributes]::ReparsePoint) {
        throw 'Archive parent must not be a junction or symbolic link.'
    }
}

$archiveTargets = @(
    'allinone.md', 'ideaspark_run', 'artifacts/mmseqs_work',
    'artifacts/dynarepo_interface_dynamics_smoke',
    'output/itrace_residue_biology_smoke', 'output/itrace_zero_shot_context_smoke',
    'tmp/itrace_benchmark_smoke_20260724', 'tmp/fig2_analysis', 'paper/tmp'
)
$deleteTargets = @(
    '.pytest_cache', '.ruff_cache',
    'paper/.render_check', 'paper/.render_check_final', 'paper/_qa_render_final',
    'paper/figures/qa_render', 'tmp/paper_figure_check', 'tmp/paper_render',
    'paper/compile_context_pass1.log', 'paper/compile_context_pass2.log',
    'paper/compile_pass1.log', 'paper/compile_pass2.log',
    'paper/paper.aux', 'paper/paper.blg', 'paper/paper.fdb_latexmk',
    'paper/paper.fls', 'paper/paper.log', 'paper/paper.out', 'paper/paper.xdv',
    'src/flexppi/__pycache__', 'src/flexppi/data/__pycache__',
    'src/flexppi/model/__pycache__', 'src/flexppi/models/__pycache__',
    'src/flexppi/training/__pycache__', 'src/flexppi/evaluation/__pycache__',
    'scripts/__pycache__', 'tests/__pycache__', 'tmp/__pycache__'
)
$plan = @()
foreach ($action in @('archive', 'delete', 'snapshot')) {
    $targets = switch ($action) {
        'archive' { $archiveTargets }
        'delete' { $deleteTargets }
        'snapshot' { @('README.md', 'BACKUP_MANIFEST.md') }
    }
    foreach ($relative in $targets) {
        $source = [IO.Path]::GetFullPath((Join-Path $projectRoot $relative))
        if (-not $source.StartsWith($sourcePrefix, [StringComparison]::OrdinalIgnoreCase)) {
            throw "Source escapes project: $source"
        }
        if (-not (Test-Path -LiteralPath $source)) { continue }
        # Validate every ancestor too, so nested targets cannot escape via a junction.
        $ancestor = $source
        while ($ancestor -ne $projectRoot) {
            if ((Get-Item -LiteralPath $ancestor -Force).Attributes -band [IO.FileAttributes]::ReparsePoint) {
                throw "Target or ancestor is a reparse point: $ancestor"
            }
            $ancestor = Split-Path $ancestor -Parent
        }
        $item = Get-Item -LiteralPath $source -Force
        $entries = if ($item.PSIsContainer) { @(Get-ChildItem -LiteralPath $source -Recurse -Force) } else { @($item) }
        $links = @($entries | Where-Object { $_.Attributes -band [IO.FileAttributes]::ReparsePoint })
        if ($action -eq 'delete' -and $links.Count -gt 0) {
            throw "Refusing recursive deletion containing links: $source"
        }
        $files = @($entries | Where-Object { -not $_.PSIsContainer -and -not ($_.Attributes -band [IO.FileAttributes]::ReparsePoint) })
        $destination = [IO.Path]::GetFullPath((Join-Path $archiveRoot $relative))
        if (-not $destination.StartsWith($archivePrefix, [StringComparison]::OrdinalIgnoreCase)) {
            throw "Destination escapes archive: $destination"
        }
        $plan += [PSCustomObject]@{
            action = $action; relative = $relative; source = $source; destination = $destination
            files = $files.Count; bytes = [long](($files | Measure-Object Length -Sum).Sum)
            links = $links.Count; status = 'planned'
        }
    }
}
$plan | Select-Object action, relative, files, bytes, links | Format-Table -AutoSize
if (-not $Apply) { Write-Output 'Preview only; nothing changed.'; return }

New-Item -ItemType Directory -Path $archiveRoot | Out-Null
$receiptPath = Join-Path $archiveRoot 'receipt.json'
$receipt = [ordered]@{
    project = $projectRoot; archive = $archiveRoot; started = (Get-Date -Format o)
    gitHead = (git -C $projectRoot rev-parse HEAD)
    gitStatusBefore = @(git -C $projectRoot status --short)
    note = 'Archive moves preserve regular files and link entries. Deleted caches must be regenerated. Restore archived relative paths only if their original locations are absent.'
    plan = $plan
}
$receipt | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $receiptPath -Encoding UTF8
foreach ($entry in $plan) {
    try {
        if ($entry.action -eq 'delete') {
            Remove-Item -LiteralPath $entry.source -Recurse -Force
        } else {
            $parent = Split-Path $entry.destination -Parent
            if (-not (Test-Path -LiteralPath $parent)) { New-Item -ItemType Directory -Path $parent -Force | Out-Null }
            if (Test-Path -LiteralPath $entry.destination) { throw "Destination exists: $($entry.destination)" }
            if ($entry.action -eq 'snapshot') {
                Copy-Item -LiteralPath $entry.source -Destination $entry.destination
            } else {
                # Same-volume directory move does not traverse or delete contained links.
                Move-Item -LiteralPath $entry.source -Destination $entry.destination
            }
            $archived = Get-Item -LiteralPath $entry.destination -Force
            $archivedFiles = if ($archived.PSIsContainer) {
                @(Get-ChildItem -LiteralPath $entry.destination -Recurse -File -Force | Where-Object { -not ($_.Attributes -band [IO.FileAttributes]::ReparsePoint) })
            } else { @($archived) }
            $archivedBytes = [long](($archivedFiles | Measure-Object Length -Sum).Sum)
            if ($archivedFiles.Count -ne $entry.files -or $archivedBytes -ne $entry.bytes) {
                throw "Archive count/size verification failed: $($entry.relative)"
            }
        }
        $entry.status = 'completed'
    } catch {
        $entry.status = 'failed: ' + $_.Exception.Message
        throw
    } finally {
        $receipt | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $receiptPath -Encoding UTF8
    }
}
Write-Output "Receipt: $receiptPath"

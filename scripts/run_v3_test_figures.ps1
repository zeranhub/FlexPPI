param([switch]$SkipInference)

$ErrorActionPreference = 'Stop'
$taskRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$taskPython = Join-Path $taskRoot '.venv-md\Scripts\python.exe'
Set-Location -LiteralPath $taskRoot

function Invoke-TaskPython {
    param([string[]]$Arguments)
    & $taskPython @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "Python task failed with exit code $LASTEXITCODE"
    }
}

if (-not $SkipInference) {
    Invoke-TaskPython -Arguments @(
        'scripts/export_v3_oof_predictions.py',
        '--output-dir', 'artifacts/evaluation_v3/oof',
        '--path-map', "/ranze/repos/FlexPPI/artifacts/data=$taskRoot/artifacts/data",
        '--path-map', '/ranze/data/ATLAS_extracted=E:/data/protein_md_datasets/ATLAS_extracted',
        '--path-map', '/ranze/data/dynarepo/dynarepo_non_nucleic_full=D:/data/dynarepo_non_nucleic_full',
        '--path-map', '/ranze/data/ITrace=E:/data/ITrace'
    )
}

Invoke-TaskPython -Arguments @('scripts/analyze_v3_rmsf_tests.py')
Invoke-TaskPython -Arguments @('scripts/render_v3_rmsf_structures.py', '--vmax', '5')
Invoke-TaskPython -Arguments @('scripts/plot_v3_rmsf_reference_figures.py', '--label-sd', '--vmax', '5')

$taskFigures = Join-Path $taskRoot 'output/flexppi_v3_test_figures'
foreach ($taskFile in @('TEST_REPORT.md', 'selected_cases.csv', 'selection_protocol.json', 'five_fold_mean_sd.csv')) {
    Copy-Item -LiteralPath (Join-Path $taskRoot "artifacts/evaluation_v3/analysis/$taskFile") -Destination $taskFigures
}
Write-Output "Finished. Figures and report: $taskFigures"

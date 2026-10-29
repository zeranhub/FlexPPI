param([int]$WaitForBackFlip2Pid = 0)
$ErrorActionPreference = 'Stop'
$taskDirectory = Split-Path -Parent $MyInvocation.MyCommand.Path
$pythonExecutable = Join-Path (Split-Path -Parent $taskDirectory) 'runtime_py312/Scripts/python.exe'
if ($WaitForBackFlip2Pid -gt 0) {
    Wait-Process -Id $WaitForBackFlip2Pid -ErrorAction SilentlyContinue
}
$backflip2Metadata = Get-Content -LiteralPath (Join-Path $taskDirectory 'results/backflip-2.1/run_metadata.json') -Raw | ConvertFrom-Json
if (-not $backflip2Metadata.completed_utc) {
    throw 'BackFlip-2 has not recorded completion. Inspect its log before starting further variants.'
}
foreach ($variantName in @('backflip-1.0', 'backflip-1.0-seq')) {
    $variantCheckpoint = "E:/Project/FlexPPI_baselines/model_cache/backflip/$variantName/$variantName.ckpt"
    $variantLog = Join-Path $taskDirectory ($variantName + '_full_local.log')
    & $pythonExecutable -u (Join-Path $taskDirectory 'run_benchmark.py') `
        --source 'E:/Project/FlexPPI_baselines/BackFlip' `
        --checkpoint $variantCheckpoint --tag $variantName --context both --threads 2 *> $variantLog
    if ($LASTEXITCODE -ne 0) {
        throw "Variant $variantName exited with $LASTEXITCODE. See $variantLog."
    }
    & $pythonExecutable -u (Join-Path $taskDirectory 'validate_outputs.py') --tag $variantName `
        *> (Join-Path $taskDirectory ($variantName + '_validation.log'))
    if ($LASTEXITCODE -ne 0) {
        throw "Variant $variantName failed its identity/hash validation."
    }
}
& $pythonExecutable -u (Join-Path $taskDirectory 'summarize_results.py') `
    *> (Join-Path $taskDirectory 'all_variants_summary.log')

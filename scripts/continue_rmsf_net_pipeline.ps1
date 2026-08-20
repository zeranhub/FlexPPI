param(
    [int]$PollSeconds = 30
)

$ErrorActionPreference = 'Stop'
$repoRoot = Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $repoRoot

$python = Join-Path $repoRoot '.venv-md\Scripts\python.exe'
$baselineRoot = Join-Path $repoRoot 'artifacts\evaluation\additional_baselines'
$rawRoot = Join-Path $baselineRoot 'raw'
$staticNames = 0..2 | ForEach-Object { "flexppi-rmsf-static-$_" }
$itraceNames = 0..2 | ForEach-Object { "flexppi-rmsf-itrace-$_" }

function Get-ContainerState {
    param([string]$Name)
    $json = docker inspect $Name --format '{{json .State}}'
    if ($LASTEXITCODE -ne 0) {
        throw "Cannot inspect container $Name"
    }
    return $json | ConvertFrom-Json
}

function Wait-ContainerGroup {
    param([string[]]$Names)
    while ($true) {
        $pending = $false
        foreach ($name in $Names) {
            $state = Get-ContainerState -Name $name
            if ($state.Status -in @('created', 'running', 'restarting')) {
                $pending = $true
                continue
            }
            if ($state.Status -ne 'exited' -or [int]$state.ExitCode -ne 0) {
                throw "$name ended with status=$($state.Status), exit=$($state.ExitCode)"
            }
        }
        if (-not $pending) {
            return
        }
        Start-Sleep -Seconds $PollSeconds
    }
}

function Assert-NativeSuccess {
    param([string]$Label)
    if ($LASTEXITCODE -ne 0) {
        throw "$Label failed with exit code $LASTEXITCODE"
    }
}

Write-Output "[$(Get-Date -Format o)] Waiting for static RMSF-net shards"
Wait-ContainerGroup -Names $staticNames

$staticInputs = 0..2 | ForEach-Object {
    Join-Path $rawRoot "rmsf_net_static_shard$_.csv"
}
& $python scripts\merge_baseline_shards.py `
    --inputs $staticInputs `
    --output (Join-Path $rawRoot 'rmsf_net_static.csv') `
    --expected-samples 276
Assert-NativeSuccess -Label 'Static shard merge'

foreach ($index in 0..2) {
    $name = $itraceNames[$index]
    $existing = docker ps -a --filter "name=^/$name`$" --format '{{.Names}}'
    Assert-NativeSuccess -Label "Container lookup for $name"
    if ($existing) {
        $state = Get-ContainerState -Name $name
        if ($state.Status -ne 'running') {
            docker start $name | Out-Null
            Assert-NativeSuccess -Label "Restart $name"
        }
        continue
    }
    docker run -d `
        --name $name `
        -v 'E:\data:/data_e:ro' `
        -v 'D:\data:/data_d:ro' `
        -v "${baselineRoot}:/baseline" `
        flexppi-baseline/rmsf-net:pdb-eb353b0 `
        --manifest /baseline/manifests/itrace.tsv `
        --output "/baseline/raw/rmsf_net_itrace_shard$index.csv" `
        --shard-index $index `
        --shard-count 3 `
        --threads 4 `
        --batch-size 2 | Out-Null
    Assert-NativeSuccess -Label "Start $name"
}

Write-Output "[$(Get-Date -Format o)] Waiting for ITrace RMSF-net shards"
Wait-ContainerGroup -Names $itraceNames

$itraceInputs = 0..2 | ForEach-Object {
    Join-Path $rawRoot "rmsf_net_itrace_shard$_.csv"
}
& $python scripts\merge_baseline_shards.py `
    --inputs $itraceInputs `
    --output (Join-Path $rawRoot 'rmsf_net_itrace.csv') `
    --expected-samples 244
Assert-NativeSuccess -Label 'ITrace shard merge'

& $python scripts\evaluate_additional_rmsf_baselines.py `
    --bootstrap-repetitions 2000 `
    --seed 2026
Assert-NativeSuccess -Label 'RMSF baseline evaluation'

& $python scripts\evaluate_alphaflow_ensembles.py `
    --ensemble-dir tmp\third_party\alphaflow\samples\alphaflow_12l_md_templates_distilled_202406 `
    --bootstrap-repetitions 2000 `
    --seed 2026
Assert-NativeSuccess -Label 'AlphaFlow overlap evaluation'
Write-Output "[$(Get-Date -Format o)] RMSF-net pipeline complete"

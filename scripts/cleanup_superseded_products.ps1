param([switch]$Apply)
# User requested reversible staging only. This compatibility entry never deletes.
$ErrorActionPreference = 'Stop'
& (Join-Path $PSScriptRoot 'stage_superseded_products.ps1') @PSBoundParameters
exit $LASTEXITCODE

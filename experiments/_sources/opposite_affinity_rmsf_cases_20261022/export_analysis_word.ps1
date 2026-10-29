$ErrorActionPreference = 'Stop'
$taskRoot = $PSScriptRoot
$documentPath = Join-Path $taskRoot '局部柔性下降与相反亲和结果_案例分析_三线表_20261022.docx'
$qaPath = Join-Path $taskRoot 'document_QA'
New-Item -ItemType Directory -Path $qaPath -Force | Out-Null
$pdfPath = Join-Path $qaPath 'case_analysis_render.pdf'
$wordInstance = $null
$reportDocument = $null
try {
    $wordInstance = New-Object -ComObject Word.Application
    $wordInstance.Visible = $false
    $wordInstance.DisplayAlerts = 0
    $reportDocument = $wordInstance.Documents.Open($documentPath, $false, $true)
    $reportDocument.Repaginate()
    $pageCount = $reportDocument.ComputeStatistics(2)
    $reportDocument.ExportAsFixedFormat($pdfPath, 17)
    Write-Output "Exported analysis report PDF; pages=$pageCount"
    Write-Output $pdfPath
}
finally {
    if ($null -ne $reportDocument) { $reportDocument.Close(0) }
    if ($null -ne $wordInstance) { $wordInstance.Quit(0) }
}

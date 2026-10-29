$ErrorActionPreference = 'Stop'
$taskDocx='E:\Project\FlexPPI\paper\revision_v14_20261010\paper_v53_20261022_中文修订_MEL5残基接触与双动态联合图_三线表.docx'
$taskQa=Join-Path $PSScriptRoot 'manuscript_QA'
New-Item -ItemType Directory -Path $taskQa -Force | Out-Null
$taskPdf=Join-Path $taskQa 'v53_word.pdf'
$taskWord=New-Object -ComObject Word.Application
$taskWord.Visible=$false
$taskWord.DisplayAlerts=0
$taskDocument=$null
try {
    $taskWord.Options.PrintComments=$false
    $taskDocument=$taskWord.Documents.Open($taskDocx,$false,$true)
    $taskDocument.Repaginate()
    $taskDocument.ExportAsFixedFormat($taskPdf,17)
    Write-Output ('Pages: '+$taskDocument.ComputeStatistics(2))
    Write-Output $taskPdf
} finally {
    if($null -ne $taskDocument){$taskDocument.Close(0)}
    $taskWord.Quit()
    [Runtime.InteropServices.Marshal]::ReleaseComObject($taskWord)|Out-Null
}

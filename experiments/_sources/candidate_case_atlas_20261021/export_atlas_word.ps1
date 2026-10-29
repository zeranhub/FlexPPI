$taskDocx='E:\Project\FlexPPI\paper\revision_v14_20261010\FlexPPI_32个潜在案例图册_背景证据与局限_20261021.docx'
$taskTmp='C:\Users\DELL\AppData\Local\Temp\flexppi_case_atlas_20261021'
New-Item -ItemType Directory -Path $taskTmp -Force | Out-Null
$taskPdf=Join-Path $taskTmp 'atlas_word.pdf'
$taskWord=New-Object -ComObject Word.Application
$taskWord.Visible=$false
$taskWord.DisplayAlerts=0
$taskDocument=$null
try {
    $taskWord.Options.PrintComments=$false
    $taskDocument=$taskWord.Documents.Open($taskDocx,$false,$true)
    $taskDocument.ExportAsFixedFormat($taskPdf,17)
    Write-Output ('Pages: '+$taskDocument.ComputeStatistics(2))
} finally {
    if($null -ne $taskDocument){$taskDocument.Close(0)}
    $taskWord.Quit()
    [Runtime.InteropServices.Marshal]::ReleaseComObject($taskWord)|Out-Null
}

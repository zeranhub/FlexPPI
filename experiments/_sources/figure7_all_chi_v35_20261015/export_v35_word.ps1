$taskDocx='E:\Project\FlexPPI\paper\revision_v14_20261010\paper_v35_20261015_中文修订_图7全残基逐χ占据关联_三线表.docx'
$taskPdf='C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010\v35_word.pdf'
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

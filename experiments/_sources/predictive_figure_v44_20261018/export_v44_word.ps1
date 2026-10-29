$taskDocx='E:\Project\FlexPPI\paper\revision_v14_20261010\paper_v44_20261018_中文修订_预测一致性与高动态识别_三线表.docx'
$taskPdf='C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010\v44_word.pdf'
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

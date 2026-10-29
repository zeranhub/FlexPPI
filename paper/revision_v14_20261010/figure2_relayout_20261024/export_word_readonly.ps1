param([Parameter(Mandatory=$true)][string]$DocumentPath,
      [Parameter(Mandatory=$true)][string]$PdfPath)
$ErrorActionPreference = 'Stop'
[Console]::OutputEncoding = [System.Text.UTF8Encoding]::new($false)
$taskWord = $null
$taskDocument = $null
try {
    $taskWord = New-Object -ComObject Word.Application
    $taskWord.Visible = $false
    $taskWord.DisplayAlerts = 0
    $taskWord.AutomationSecurity = 3
    $taskDocument = $taskWord.Documents.Open($DocumentPath, $false, $true)
    $taskDocument.ExportAsFixedFormat($PdfPath, 17)
    Write-Output ('Pages: ' + $taskDocument.ComputeStatistics(2))
    Write-Output ('PDF: ' + $PdfPath)
} finally {
    try {
        if ($null -ne $taskDocument) {
            $taskDoNotSave = 0
            $taskDocument.Close([ref]$taskDoNotSave)
            [Runtime.InteropServices.Marshal]::FinalReleaseComObject($taskDocument) | Out-Null
        }
    } finally {
        if ($null -ne $taskWord) {
            $taskWord.Quit()
            [Runtime.InteropServices.Marshal]::FinalReleaseComObject($taskWord) | Out-Null
        }
    }
}

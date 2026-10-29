param([Parameter(Mandatory=$true)][string]$DocumentPath)
$ErrorActionPreference = 'Stop'
$taskWordCleanup = $null
try {
    $taskWordCleanup = [Runtime.InteropServices.Marshal]::GetActiveObject('Word.Application')
    if ($taskWordCleanup.Visible) { throw 'Refusing cleanup: Word is visible.' }
    $taskDocsCleanup = @($taskWordCleanup.Documents)
    foreach ($taskDocCleanup in $taskDocsCleanup) {
        if (-not $taskDocCleanup.ReadOnly -or $taskDocCleanup.FullName -ne $DocumentPath) {
            throw 'Refusing cleanup: found a document other than the exact read-only export source.'
        }
    }
    foreach ($taskDocCleanup in $taskDocsCleanup) {
        $taskNoSaveCleanup = 0
        $taskDocCleanup.Close([ref]$taskNoSaveCleanup)
    }
    if ($taskWordCleanup.Documents.Count -eq 0) {
        $taskWordCleanup.Quit()
        Write-Output 'Closed only the verified read-only export document and empty hidden Word instance.'
    }
} catch {
    Write-Output $_.Exception.Message
} finally {
    if ($null -ne $taskWordCleanup) {
        [Runtime.InteropServices.Marshal]::FinalReleaseComObject($taskWordCleanup) | Out-Null
    }
}

$ErrorActionPreference='Stop'
$env:PYTHONIOENCODING='utf-8'
$taskRoot=$PSScriptRoot
$taskRuntime='C:\Users\DELL\.cache\flexppi-figure-tools\paperbanana-venv\Scripts\python.exe'
$taskRepo='C:\Users\DELL\.cache\flexppi-figure-tools\PaperBanana'
$taskKeyFile='C:\Users\DELL\.cache\flexppi-figure-tools\credentials\gemini-session.key'
$taskOutput=Join-Path $taskRoot ('generation_retry_'+(Get-Date -Format 'yyyyMMdd_HHmmss'))
$taskArgs=@((Join-Path $taskRoot 'run_from_scratch.py'),'--method',(Join-Path $taskRoot 'scientific_content.txt'),
 '--caption',(Join-Path $taskRoot 'visual_intent.txt'),'--output',$taskOutput,'--repo',$taskRepo,
 '--aspect-ratio','1:1','--resolution','4K','--critic-rounds','2','--candidates','2','--retrieval','auto')
if ([string]::IsNullOrWhiteSpace($env:GOOGLE_API_KEY)) {$taskArgs+=@('--key-file',$taskKeyFile)}
& $taskRuntime @taskArgs
if ($LASTEXITCODE -ne 0) {throw 'PaperBanana did not complete. See the safe failure.json in the output folder.'}

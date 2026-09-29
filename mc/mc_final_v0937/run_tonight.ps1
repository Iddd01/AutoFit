# run_tonight.ps1 -- final Monte Carlo of xtdpthresh 0.9.37 (PowerShell 5.1, Stata 17 MP).
# Runs, one after another, for each registry of
# final/ (point estimation, then inference) a smoke (execution and merge
# contracts only) followed, if the smoke merges, by the formal run and its
# verified merge. A failed step is logged and the next
# study continues. Keep every other Stata closed: the script waits for all
# StataMP-64 processes to finish between steps.
[CmdletBinding()]
param(
 [string]$Stata='C:\Program Files\Stata17\StataMP-64.exe',
 [ValidateRange(1,64)][int]$NShard=28,
 [string]$Tag='0937'
)
$ErrorActionPreference='Stop'
Set-StrictMode -Version 2.0
$root=$PSScriptRoot
$log=Join-Path $root ("run_tonight_"+(Get-Date -Format 'yyyyMMdd_HHmm')+".log")
function Note($m){$line=(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')+'  '+$m; Write-Host $line; Add-Content -LiteralPath $log -Value $line}
function Wait-Stata {
 # all shards of a step must finish; two empty polls a minute apart
 Start-Sleep -Seconds 90
 $empty=0
 $deadline=(Get-Date).AddHours(24)
 while($empty -lt 2){
  if((Get-Date) -gt $deadline){throw 'Stata still running after 24 hours'}
  $n=@(Get-Process StataMP-64 -ErrorAction SilentlyContinue).Count
  if($n -eq 0){$empty++}else{$empty=0}
  Start-Sleep -Seconds 60
 }
}
function Step($name,[scriptblock]$body){
 Note "START $name"
 # launcher output goes to the console, so Step returns only the boolean
 try {& $body | Out-Host; Note "OK    $name"; return $true}
 catch {
  Note "FAIL  $name : $($_.Exception.Message)"
  # a launch that failed part way can leave shards running: let them finish
  # before the next block starts (the incomplete run can be resumed later)
  try {Wait-Stata} catch {Note "      $($_.Exception.Message)"}
  return $false
 }
}
if(@(Get-Process StataMP-64 -ErrorAction SilentlyContinue).Count -gt 0){throw 'Close every Stata before starting.'}
Note "run_tonight: root=$root NShard=$NShard"

# preflight: xthenreg and moremata (XTH cells) must be installed
$pre=Join-Path $root 'preflight.log'
Remove-Item -LiteralPath $pre -ErrorAction SilentlyContinue
$null=Start-Process -FilePath $Stata -ArgumentList '/e do preflight.do' -WorkingDirectory $root -WindowStyle Hidden -PassThru -Wait
if(-not (Select-String -LiteralPath $pre -Pattern '^PREFLIGHT_PASS\s*$' -Quiet)){
 Note 'preflight failed: install xthenreg and moremata (ssc install xthenreg; ssc install moremata); see preflight.log'
 throw 'preflight failed'
}
Note 'preflight: xthenreg and moremata found'

# 1-2. final registries: point estimation, then inference
$F=Join-Path $root 'final'
function Block($label,$registry,$smoke,$formal){
 $ok=Step "$label smoke" {
  & (Join-Path $F 'run_final.ps1') -Action Fresh -Registry $registry -RunId $smoke -NShard 4 -RepCap 1 -B 19 -Grid 10 -GridCI 10 -Stata $Stata
  Wait-Stata
  & (Join-Path $F 'run_final.ps1') -Action Merge -RunId $smoke -Stata $Stata
 }
 if($ok){
  $null=Step "$label formal" {
   & (Join-Path $F 'run_final.ps1') -Action Fresh -Registry $registry -RunId $formal -NShard $NShard -Stata $Stata
   Wait-Stata
   & (Join-Path $F 'run_final.ps1') -Action Merge -RunId $formal -Stata $Stata
  }
 }
}
Block 'point estimation' 'final_point_cells.csv' "point_smoke_$Tag" "point_final_$Tag"
Block 'inference' 'final_inf_cells.csv' "inf_smoke_$Tag" "inf_final_$Tag"
Note 'run_tonight: finished'

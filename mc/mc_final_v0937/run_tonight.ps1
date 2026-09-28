# run_tonight.ps1 -- overnight sequence for xtdpthresh 0.9.37 (PowerShell 5.1, Stata 17 MP).
# Runs, one after another: the version check; then for each study a smoke
# (execution and merge contracts only) followed, if the smoke merges, by the
# formal run and its verified merge. A failed step is logged and the next
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
 while($empty -lt 2){
  $n=@(Get-Process StataMP-64 -ErrorAction SilentlyContinue).Count
  if($n -eq 0){$empty++}else{$empty=0}
  Start-Sleep -Seconds 60
 }
}
function Step($name,[scriptblock]$body){
 Note "START $name"
 # launcher output goes to the console, so Step returns only the boolean
 try {& $body | Out-Host; Note "OK    $name"; return $true}
 catch {Note "FAIL  $name : $($_.Exception.Message)"; return $false}
}
if(@(Get-Process StataMP-64 -ErrorAction SilentlyContinue).Count -gt 0){throw 'Close every Stata before starting.'}
Note "run_tonight: root=$root NShard=$NShard"

# 0. version check: 0.9.35/0.9.36 vs 0.9.37 on common samples
$null=Step 'version_check' {
 $dir=Join-Path $root 'version_check'
 $p=Start-Process -FilePath $Stata -ArgumentList '/e do version_check.do' -WorkingDirectory $dir -WindowStyle Hidden -PassThru -Wait
 $vlog=Join-Path $dir 'version_check.log'
 if(-not (Select-String -LiteralPath $vlog -Pattern 'VERSION_CHECK_PASS' -SimpleMatch -Quiet)){throw 'VERSION_CHECK_PASS not found in version_check.log'}
 (Select-String -LiteralPath $vlog -Pattern 'max reldif' -SimpleMatch) | ForEach-Object {Note ('      '+$_.Line.Trim())}
}

# 1. Study A core (point estimation, 116 cells)
$A=Join-Path $root 'study_a'
Push-Location $A
try {
 $smoke="smoke_a_$Tag"; $formal="final_a_$Tag"
 $ok=Step 'study_a smoke' {
  & (Join-Path $A 'run_study_a.ps1') -Fresh -RunId $smoke -NShard 6 -RepCap 1 -Stata $Stata
  Wait-Stata
  & (Join-Path $A "runs\$smoke\verify_and_merge_study_a.ps1") -RunId $smoke
 }
 if($ok){
  $null=Step 'study_a formal' {
   & (Join-Path $A 'run_study_a.ps1') -Fresh -RunId $formal -NShard $NShard -Stata $Stata
   Wait-Stata
   & (Join-Path $A "runs\$formal\verify_and_merge_study_a.ps1") -RunId $formal
  }
 }
} finally {Pop-Location}

# 2..5: blocks run through the run_supplement.ps1 / run_supp2.ps1 launchers
function Block($label,$dir,$launcher,$smoke,$formal,$extra,$smokeOpts){
 $ok=Step "$label smoke" {
  & (Join-Path $dir $launcher) -Action Fresh -RunId $smoke -NShard 4 -RepCap 1 @smokeOpts @extra -Stata $Stata
  Wait-Stata
  & (Join-Path $dir $launcher) -Action Merge -RunId $smoke -Stata $Stata
 }
 if($ok){
  $null=Step "$label formal" {
   & (Join-Path $dir $launcher) -Action Fresh -RunId $formal -NShard $NShard @extra -Stata $Stata
   Wait-Stata
   & (Join-Path $dir $launcher) -Action Merge -RunId $formal -Stata $Stata
  }
 }
}
$quick=@{B=19;Grid=10;GridCI=10}
$none=@{}
$B=Join-Path $root 'study_b_supp'
Block 'study_a supplement' $A 'run_supplement.ps1' "a_supp_smoke_$Tag" "a_supp_final_$Tag" $none @{Grid=10}
Block 'kink supplement' $B 'run_supplement.ps1' "kink_smoke_$Tag" "kink_final_$Tag" $none $quick
Block 'power (FD/FOD, maxlag 1 3)' $B 'run_supp2.ps1' "power_smoke_$Tag" "power_final_$Tag" $none $quick
Block 'power (Gong-Seo geometry)' $B 'run_supp2.ps1' "powergs_smoke_$Tag" "powergs_final_$Tag" @{Registry='supp2b_cells.csv'} $quick
Note 'run_tonight: finished'

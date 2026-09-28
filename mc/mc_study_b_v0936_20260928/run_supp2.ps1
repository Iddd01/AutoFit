# Supplement 2 (power of the threshold test at gamma0 + c, FD vs FOD); PowerShell 5.1 and Stata 17.
[CmdletBinding()]
param(
 [ValidateSet('Fresh','Resume','Status','Merge')][string]$Action='Status',
 [Parameter(Mandatory=$true)][ValidatePattern('^[A-Za-z][A-Za-z0-9_-]{0,70}$')][string]$RunId,
 [ValidateRange(1,64)][int]$NShard=28,
 [ValidateRange(0,500)][int]$RepCap=0,
 [ValidateRange(0,10000)][int]$B=0,
 [ValidateRange(0,10000)][int]$Grid=0,
 [ValidateRange(0,10000)][int]$GridCI=0,
 [string]$Stata='C:\Program Files\Stata17\StataMP-64.exe'
)
$ErrorActionPreference='Stop'
Set-StrictMode -Version 2.0
$source=$PSScriptRoot
$isStaged=Test-Path -LiteralPath (Join-Path $source 'manifest.json')
if($isStaged -and $Action -eq 'Fresh'){throw 'Launch a fresh run from the release folder, not a frozen run.'}
$run=if($isStaged){$source}else{Join-Path $source ('supp2_runs/'+$RunId)}
$members=@('supp2_worker.do','supp2_merge.do','supp2_cells.csv',
 'run_supp2.ps1','supp2_link_check.do','xtdpthresh.ado','xtdpthresh_p.ado','xtdpthresh.sthlp','SUPP2.md')
function Hashes($directory,$names) {
 $map=[ordered]@{}
 foreach($name in $names){$map[$name]=(Get-FileHash -LiteralPath (Join-Path $directory $name) -Algorithm SHA256).Hash}
 return $map
}
function Write-Utf8($path,$value) {
 [IO.File]::WriteAllText($path,$value,[Text.UTF8Encoding]::new($false))
}
# jobs.json stores the start time in ISO 8601 ('o'); ConvertFrom-Json may
# return it as a string or a DateTime, and culture-dependent Parse can fail.
function Started-Utc($v) {
 if($v -is [DateTime]){return $v.ToUniversalTime()}
 return [DateTime]::Parse([string]$v,[Globalization.CultureInfo]::InvariantCulture,
  [Globalization.DateTimeStyles]::RoundtripKind).ToUniversalTime()
}
function Live-Jobs {
 $jobsPath=Join-Path $run 'jobs.json'
 if(-not(Test-Path -LiteralPath $jobsPath)){return @()}
 $jobs=@(Get-Content -LiteralPath $jobsPath -Raw | ConvertFrom-Json)
 return @($jobs | Where-Object {
   $p=Get-Process -Id $_.pid -ErrorAction SilentlyContinue
   $null -ne $p -and $p.ProcessName -eq 'StataMP-64' -and
     [Math]::Abs(($p.StartTime.ToUniversalTime()-(Started-Utc $_.started)).TotalSeconds) -lt 2
 })
}
if($Action -eq 'Fresh') {
 if(Test-Path -LiteralPath $run){throw 'Run already exists; use a new RunId or Resume.'}
 if(($B -gt 0 -and $B -lt 10)-or($Grid -gt 0 -and $Grid -lt 10)-or($GridCI -gt 0 -and $GridCI -lt 10)){throw 'Positive B/Grid/GridCI overrides must be at least 10.'}
 $cells=@(Import-Csv -LiteralPath (Join-Path $source 'supp2_cells.csv'))
 $expected=0
 foreach($c in $cells){$expected+=if($RepCap){[Math]::Min($RepCap,[int]$c.R)}else{[int]$c.R}}
 if($NShard -gt $expected){throw 'NShard exceeds number of requested fits.'}
 New-Item -ItemType Directory -Path (Split-Path -Parent $run) -Force | Out-Null
 New-Item -ItemType Directory -Path $run | Out-Null
 foreach($f in $members){Copy-Item -LiteralPath (Join-Path $source $f) -Destination (Join-Path $run $f)}
 $config=@{
  RunId=$RunId;NShard=$NShard;RepCap=$RepCap;B=$B;Grid=$Grid;GridCI=$GridCI;
  Master=20260814;Expected=$expected;Cells=$cells.Count;Study='SUPP2';
  Formal=($RepCap -eq 0 -and $B -eq 0 -and $Grid -eq 0 -and $GridCI -eq 0)
 }
 $configDo=@"
global supp_run $RunId
global supp_nshard $NShard
global supp_repcap $RepCap
global supp_B $B
global supp_grid $Grid
global supp_gridci $GridCI
global supp_master 20260814
global supp_expected $expected
global supp_formal $([int]$config.Formal)
"@
 Write-Utf8 (Join-Path $run 'supp2_config.do') $configDo
 $members+= 'supp2_config.do'
 $manifest=[ordered]@{schema=1;config=$config;Stata=(Resolve-Path -LiteralPath $Stata).Path;
  StataSHA256=(Get-FileHash -LiteralPath $Stata -Algorithm SHA256).Hash;hashes=(Hashes $run $members)}
 Write-Utf8 (Join-Path $run 'manifest.json') ($manifest | ConvertTo-Json -Depth 8)
}
if(-not(Test-Path -LiteralPath (Join-Path $run 'manifest.json'))){throw 'Unknown run.'}
$manifest=Get-Content -LiteralPath (Join-Path $run 'manifest.json') -Raw | ConvertFrom-Json
if($manifest.config.RunId -cne $RunId){throw 'Manifest RunId mismatch.'}
if((Get-FileHash -LiteralPath $MyInvocation.MyCommand.Path).Hash -cne $manifest.hashes.'run_supp2.ps1'){throw "Launcher changed; invoke the frozen run_supp2.ps1 inside $run."}
$lock=[IO.File]::Open((Join-Path $run '.runlock'),[IO.FileMode]::OpenOrCreate,[IO.FileAccess]::ReadWrite,[IO.FileShare]::None)
try {
 if($Action -eq 'Merge') {
  foreach($name in @('attestation.json','supp2_merge.ok')) {
   $stale=Join-Path $run $name
   if(Test-Path -LiteralPath $stale -PathType Leaf){Remove-Item -LiteralPath $stale}
  }
 }
 foreach($entry in $manifest.hashes.PSObject.Properties){
  if((Get-FileHash -LiteralPath (Join-Path $run $entry.Name) -Algorithm SHA256).Hash -cne $entry.Value){throw "Frozen file changed: $($entry.Name)"}
 }
 $Stata=$manifest.Stata
 if((Get-FileHash -LiteralPath $Stata -Algorithm SHA256).Hash -cne $manifest.StataSHA256){throw 'Stata executable changed.'}
 $live=@(Live-Jobs)
 if($Action -eq 'Status') {
  $done=@(Get-ChildItem -LiteralPath $run -Filter 'supp2_done_SH*.txt')
  Write-Output "Run=$RunId; completed=$($done.Count)/$($manifest.config.NShard); live=$($live.Count); formal=$($manifest.config.Formal)"
  return
 }
 if($live.Count){throw 'Run still has live workers; wait before Resume or Merge.'}
 if($Action -eq 'Merge') {
  for($i=1;$i -le $manifest.config.NShard;$i++){
   $expectedShard=[Math]::Floor(($manifest.config.Expected-$i)/$manifest.config.NShard)+1
   $done=(Get-Content -LiteralPath (Join-Path $run "supp2_done_SH$i.txt") -Raw).Trim()
   if($done -cne "$RunId,$i,$expectedShard"){throw "Bad completion marker for shard $i"}
  }
  $rawNames=@(1..$manifest.config.NShard | ForEach-Object {"supp2_SH$_.csv"})
  $before=Hashes $run $rawNames
  $marker=Join-Path $run 'supp2_merge.ok'
  if(Test-Path -LiteralPath $marker){Remove-Item -LiteralPath $marker}
  $merge=Start-Process -FilePath $Stata -ArgumentList '-e do supp2_merge.do' -WorkingDirectory $run -WindowStyle Hidden -PassThru
  $merge.WaitForExit()
  if(-not(Test-Path -LiteralPath $marker)){throw 'Merge failed; inspect supp2_merge.log.'}
  if((Get-Content -LiteralPath $marker -Raw).Trim() -cne "$RunId,$($manifest.config.Expected)"){throw 'Bad merge marker.'}
  $after=Hashes $run $rawNames
  foreach($key in $before.Keys){if($before[$key] -cne $after[$key]){throw 'Raw input changed during merge.'}}
  foreach($entry in $manifest.hashes.PSObject.Properties) {
   if((Get-FileHash -LiteralPath (Join-Path $run $entry.Name)).Hash -cne $entry.Value){throw 'Source changed during merge.'}
  }
  $outputs=Hashes $run @('supp2_all.dta','supp2_all.csv','supp2_summary.csv','supp2_paired.csv')
  $att=[ordered]@{verifiedUTC=[DateTime]::UtcNow.ToString('o');config=$manifest.config;
    manifestSHA256=(Get-FileHash -LiteralPath (Join-Path $run 'manifest.json')).Hash;
    sourceHashes=$manifest.hashes;rawHashes=$after;outputHashes=$outputs}
  Write-Utf8 (Join-Path $run 'attestation.json') ($att | ConvertTo-Json -Depth 8)
  Write-Output "Verified supplement merge: $run"
  return
 }
 $jobs=@()
 for($i=1;$i -le $manifest.config.NShard;$i++) {
  $marker=Join-Path $run "supp2_done_SH$i.txt"
  if(Test-Path -LiteralPath $marker){continue}
  $deadline=[DateTime]::UtcNow.AddSeconds(60)
  do {
   $mem=Get-CimInstance Win32_PerfFormattedData_PerfOS_Memory
   $headroom=([double]$mem.CommitLimit-[double]$mem.CommittedBytes)/1MB
   if($mem.AvailableMBytes -ge 1536 -and $headroom -ge 1536){break}
   if([DateTime]::UtcNow -ge $deadline){throw 'Memory gate stopped launch; existing workers continue. Resume later.'}
   Start-Sleep -Seconds 2
  } while($true)
  $job="clear all`nset more off`ncd `"$($run.Replace('\','/'))`"`ndo supp2_worker.do $i`nexit, clear`n"
  Write-Utf8 (Join-Path $run "job_SH$i.do") $job
  $p=Start-Process -FilePath $Stata -ArgumentList "-e do job_SH$i.do" -WorkingDirectory $run -WindowStyle Hidden -PassThru
  $jobs+=@{pid=$p.Id;started=$p.StartTime.ToUniversalTime().ToString('o');shard=$i}
  Write-Utf8 (Join-Path $run 'jobs.json') (ConvertTo-Json -InputObject @($jobs) -Depth 5)
  Start-Sleep -Milliseconds 500
 }
 Write-Output "Launched $($jobs.Count) supplement shards: $RunId"
} finally {$lock.Dispose()}

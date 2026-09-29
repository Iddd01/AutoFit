# monitor_study_b.ps1 -- monitor exactly one isolated MC run (PowerShell 5.1).
[CmdletBinding()]
param([Parameter(Mandatory=$true)][ValidatePattern('^[A-Za-z0-9][A-Za-z0-9_.-]{0,79}$')][string]$RunId,
      [ValidateRange(0,86400)][int]$Watch=0)
$ErrorActionPreference='Stop';Set-StrictMode -Version 2.0
function Rows([string]$p){$r=[IO.File]::OpenText($p);try{$n=0L;while($null-ne$r.ReadLine()){$n++};return [Math]::Max(0L,$n-1L)}finally{$r.Dispose()}}
function Alive($s,[string]$field,[string]$kind){
 if($s.PSObject.Properties.Name -notcontains $field){return $false}
 $n=[int]$s.$field;if($n -le 0){return $false}
 $p=Get-Process -Id $n -ErrorAction SilentlyContinue;if($null -eq $p){return $false}
 if($kind -eq 'wrapper' -and $p.ProcessName -notlike '*powershell*'){return $false}
 if($kind -eq 'stata'){try{if($p.Path -ine $script:stata){return $false}}catch{return $false}}
 try{$u=[DateTime]::Parse([string]$s.updatedUtc).ToUniversalTime();if($p.StartTime.ToUniversalTime() -gt $u.AddSeconds(5)){return $false}}catch{return $false}
 return $true
}
$src=Split-Path -Parent $MyInvocation.MyCommand.Path
if(Test-Path -LiteralPath (Join-Path $src 'manifest.json') -PathType Leaf){$dir=$src}else{$dir=Join-Path (Join-Path $src 'runs') $RunId}
$mp=Join-Path $dir 'manifest.json'
if(-not(Test-Path -LiteralPath $mp -PathType Leaf)){throw "Manifest not found: $mp"};$m=Get-Content -LiteralPath $mp -Raw|ConvertFrom-Json
if([string]$m.config.RunId -cne $RunId){throw 'RunId/manifest mismatch'};$ns=[int]$m.config.NShard;$master=[long]$m.config.Master;$expected=[long]$m.config.ExpectedReplications;$cells=[int]$m.config.ExpectedCells;$cap=[int]$m.config.RepCap;$e5=if($cap-gt0){[Math]::Min(500,$cap)}else{500};$script:stata=(Resolve-Path -LiteralPath ([string]$m.stataExecutable)).Path
if([string]$m.config.HarnessVersion -cne 'study_b_type4_v0936'-or[string]$m.config.CodeVersionExpected-cne'0.9.36'-or$cells-ne52){throw 'Manifest is not the frozen Study B 0.9.36 registry (52 cells).'}
if($expected-ne[long]$cells*$e5){throw "Manifest replication contract is inconsistent: expected=$expected, cells=$cells, reps/cell=$e5."}
function ExpectedCells([int]$k){if($cap -eq 1 -and $ns -gt 1){if($k -gt $cells){return 0};return [int]([Math]::Floor(([double]($cells-$k))/$ns)+1)};return $cells}
function ExpectedRows([int]$k){if($cap -eq 1 -and $ns -gt 1){return [long](ExpectedCells $k)};return [long]$cells*([Math]::Floor($k*$e5/$ns)-[Math]::Floor(($k-1)*$e5/$ns))}
function Show([bool]$clear){if($clear){Clear-Host};$total=0L;$complete=0;$failed=0;$running=0;$pending=0;$out=@()
 for($k=1;$k-le$ns;$k++){$tag='{0:D3}'-f$k;$sp=Join-Path $dir "status_SH$tag.json";$csv=Join-Path $dir "study_b_results_SH$k.csv";$er=ExpectedRows $k;$ecells=ExpectedCells $k;$n=0L;$rowText='0';if(Test-Path -LiteralPath $csv){try{$n=Rows $csv;$rowText=[string]$n}catch{$rowText='locked'}};if($n-gt 0){$total+=$n};$state='pending';$ec=$null
  if(Test-Path -LiteralPath $sp){try{$s=Get-Content -LiteralPath $sp -Raw|ConvertFrom-Json;$state=[string]$s.state;$ec=$s.exitCode;if([int]$s.schemaVersion-ne1-or[string]$s.runId-cne$RunId-or[string]$s.harnessVersion-cne'study_b_type4_v0936'-or[string]$s.codeVersionExpected-cne'0.9.36'-or[long]$s.master-ne$master-or[int]$s.shard-ne$k-or[int]$s.nshard-ne$ns-or[long]$s.expectedRows-ne$er-or[int]$s.expectedCells-ne$ecells){$state='bad-status'};$alive=(Alive $s 'wrapperPid' 'wrapper')-or(Alive $s 'stataPid' 'stata');if(($state-eq'running'-or$state-eq'launching')-and-not$alive){$state='stale'};if($state-eq'completed'){$mk=Join-Path $dir "study_b_complete_SH$k.csv";if(-not(Test-Path -LiteralPath $mk -PathType Leaf)-or-not(Test-Path -LiteralPath $csv -PathType Leaf)){$state='bad-completion'}else{$cm=@(Import-Csv -LiteralPath $mk);if($cm.Count-ne1-or[string]$cm[0].run_id-cne$RunId-or[string]$cm[0].harness_version-cne'study_b_type4_v0936'-or[string]$cm[0].code_version_expected-cne'0.9.36'-or[long]$cm[0].master-ne$master-or[int]$cm[0].shard-ne$k-or[int]$cm[0].nshard-ne$ns-or[int]$cm[0].cells-ne$ecells-or[long]$cm[0].reps-ne$er-or[int]$cm[0].exit_code-ne0-or[long]$s.expectedRows-ne$er-or[long]$s.resultRows-ne$er-or$n-ne$er){$state='bad-completion'}}}}catch{$state='bad-status'}}
  switch($state){'completed'{$complete++};'running'{$running++};'launching'{$running++};'failed'{$failed++};'stale'{$failed++};'bad-status'{$failed++};'bad-completion'{$failed++};default{$pending++}};$et=if($null-ne$ec){" exit=$ec"}else{''};$out+=('  SH{0} rows={1,-7} state={2}{3}'-f$tag,$rowText,$state,$et)}
 Write-Host ('=== MC run {0} @ {1} ==='-f$RunId,(Get-Date -Format 'yyyy-MM-dd HH:mm:ss'));Write-Host "Directory: $dir";Write-Host ('Config: shards={0}, master={1}, B={2}, grid={3}, gridci={4}, repcap={5}'-f$ns,$m.config.Master,$m.config.B,$m.config.Grid,$m.config.GridCI,$m.config.RepCap);Write-Host ('Status: complete={0}, running={1}, failed/stale={2}, pending={3}'-f$complete,$running,$failed,$pending)
 if($expected-gt 0){$pct=[Math]::Min(100,[Math]::Round(100*$total/$expected,1));Write-Host "Replications: $total / $expected ($pct%)"}else{Write-Host "Replications written: $total (expected total not declared)"};Write-Host 'Shards:';$out|ForEach-Object{Write-Host $_};return [pscustomobject]@{Complete=$complete;Running=$running;Failed=$failed;Pending=$pending}}
if($Watch-eq 0){Show $false|Out-Null;exit 0}
while($true){$x=Show $true;if($x.Complete-eq$ns){Write-Host "`nAll shards completed for run $RunId.";break};if($x.Running-eq 0-and($x.Failed-gt 0-or$x.Pending-gt 0)){Write-Error 'No shard is running but the run is incomplete; inspect status/logs, then resume with identical configuration.';break};Start-Sleep -Seconds $Watch}

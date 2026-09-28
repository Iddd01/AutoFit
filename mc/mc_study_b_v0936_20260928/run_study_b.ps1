# run_study_b.ps1 -- isolated, resumable, auditable MC launcher (PowerShell 5.1).
[CmdletBinding()]
param(
    [switch]$Fresh, [switch]$Resume, [switch]$StageOnly,
    [ValidatePattern('^[A-Za-z0-9][A-Za-z0-9_.-]{0,79}$')][string]$RunId,
    [ValidateRange(1,256)][int]$NShard=28,
    [ValidateRange(1,2147483000)][long]$Master=20260814,
    [ValidateScript({$_ -eq '.' -or (($_ -match '^[0-9]+$') -and ([long]$_ -ge 10) -and ([long]$_ -le 100000))})][string]$B='.',
    [ValidateScript({$_ -eq '.' -or (($_ -match '^[0-9]+$') -and ([long]$_ -ge 10) -and ([long]$_ -le 10000))})][string]$Grid='.',
    [ValidateScript({$_ -eq '.' -or (($_ -match '^[0-9]+$') -and ([long]$_ -ge 10) -and ([long]$_ -le 10000))})][string]$GridCI='.',
    [ValidateRange(0,500)][int]$RepCap=0,
    [ValidateRange(0,1000000000)][long]$ExpectedReplications=0,
    [ValidatePattern('^[A-Za-z0-9_.-]+$')][string]$HarnessVersion='study_b_type4_v0936',
    [ValidatePattern('^[0-9]+\.[0-9]+\.[0-9]+$')][string]$CodeVersionExpected='0.9.36',
    [ValidateRange(30,3600)][int]$MemoryGateSeconds=600,
    [ValidateRange(256,32768)][int]$MinHeadroomMB=1536,
    [string]$Stata='C:\Program Files\Stata17\StataMP-64.exe'
)
$ErrorActionPreference='Stop'; Set-StrictMode -Version 2.0

function Move-AtomicWithRetry([string]$TemporaryPath,[string]$DestinationPath) {
    for($attempt=1;$attempt -le 10;$attempt++) {
        try {
            Move-Item -LiteralPath $TemporaryPath -Destination $DestinationPath -Force
            return
        } catch {
            if($attempt -eq 10){throw}
            Start-Sleep -Milliseconds (50*$attempt)
        }
    }
}
function Save-Json($Object,[string]$Path) {
    $tmp=$Path+'.tmp.'+[Guid]::NewGuid().ToString('N')
    $Object | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $tmp -Encoding UTF8
    Move-AtomicWithRetry $tmp $Path
}
function Is-Alive($Status) {
    foreach($p in @('wrapperPid','stataPid')) {
        if($Status.PSObject.Properties.Name -contains $p) {
            $id=[int]$Status.$p
            $proc=if($id -gt 0){Get-Process -Id $id -ErrorAction SilentlyContinue}else{$null}
            if($null -ne $proc) {
                if($p -eq 'wrapperPid' -and $proc.ProcessName -notlike '*powershell*'){continue}
                if($p -eq 'stataPid') {
                    try {if($proc.Path -ine $Stata){continue}} catch {continue}
                }
                try {$updated=[DateTime]::Parse([string]$Status.updatedUtc).ToUniversalTime();if($proc.StartTime.ToUniversalTime() -gt $updated.AddSeconds(5)){continue}} catch {continue}
                return $true
            }
        }
    }
    return $false
}
function Invoke-RegistryProbe([string]$Directory) {
    $probe=Join-Path $Directory '_registry_probe.do'
    $marker=Join-Path $Directory '_registry_probe.ok'
    if(-not(Test-Path -LiteralPath $probe -PathType Leaf)){throw "Registry probe is missing: $probe"}
    Remove-Item -LiteralPath $marker -Force -ErrorAction SilentlyContinue
    $a='/e do "'+$probe.Replace('"','\"')+'"'
    $p=Start-Process -FilePath $Stata -ArgumentList $a -WorkingDirectory $Directory -WindowStyle Hidden -Wait -PassThru
    $expected='cells=52;pairs=24;fdonly=4;reps=26000;R500=52;FULL=40;COVERAGE=12;B499=48;B500=4;G199=48;G46=4;GC100=40;GC46=4;GC10=8;RF4=48;RF0=4'
    if($p.ExitCode -ne 0 -or -not(Test-Path -LiteralPath $marker -PathType Leaf)){
        throw "Registry probe failed (Stata exit=$($p.ExitCode)). Inspect _registry_probe.log in $Directory."
    }
    $actual=(Get-Content -LiteralPath $marker -Raw).Trim()
    if($actual -cne $expected){throw "Registry probe marker mismatch: $actual"}
}
function Wait-MemoryHeadroom([int]$NextShard) {
    Start-Sleep -Seconds 3
    $deadline=[DateTime]::UtcNow.AddSeconds($MemoryGateSeconds)
    $probeWorked=$false;$freeMb=-1L;$commitFreeMb=-1L
    while([DateTime]::UtcNow -lt $deadline){
        try {
            $mem=Get-CimInstance Win32_PerfFormattedData_PerfOS_Memory
            $freeMb=[long]$mem.AvailableMBytes
            $commitFreeMb=[long][Math]::Floor((([double]$mem.CommitLimit)-([double]$mem.CommittedBytes))/1MB)
            $probeWorked=$true
        } catch {
            try {
                $os=Get-CimInstance Win32_OperatingSystem
                $freeMb=[long][Math]::Floor(([double]$os.FreePhysicalMemory)/1KB)
                $commitFreeMb=[long][Math]::Floor(([double]$os.FreeVirtualMemory)/1KB)
                $probeWorked=$true
            } catch {break}
        }
        if($freeMb -ge $MinHeadroomMB -and $commitFreeMb -ge $MinHeadroomMB){return}
        Start-Sleep -Seconds 2
    }
    if(-not$probeWorked){
        Write-Warning 'Memory telemetry unavailable; continuing with the three-second launch stagger.'
        return
    }
    throw "Hard memory gate stopped before shard $NextShard after $MemoryGateSeconds seconds: available RAM=${freeMb}MB, commit headroom=${commitFreeMb}MB, required=${MinHeadroomMB}MB. Existing shards continue; rerun -Resume later."
}
function Q([string]$s) { return '"'+$s.Replace('"','\"')+'"' }
function Get-ResultRows([string]$Path) {
    if(-not(Test-Path -LiteralPath $Path -PathType Leaf)){return 0L}
    $reader=[IO.File]::OpenText($Path)
    try{$n=0L;while($null-ne$reader.ReadLine()){$n++};return [Math]::Max(0L,$n-1L)}
    finally{$reader.Dispose()}
}

if(($Fresh -and $Resume)-or(-not $Fresh -and -not $Resume)){throw 'Choose exactly one: -Fresh or -Resume.'}
if($HarnessVersion -cne 'study_b_type4_v0936' -or $CodeVersionExpected -cne '0.9.36'){
    throw 'This frozen package supports only study_b_type4_v0936 with xtdpthresh 0.9.36.'
}
if($RepCap -gt 0 -and $RepCap -lt $NShard -and $RepCap -ne 1){throw 'A capped run requires RepCap >= NShard; RepCap=1 is the dedicated parallel-smoke exception.'}
if(-not [string]::IsNullOrWhiteSpace($RunId) -and $RunId -notmatch '[A-Za-z]'){throw 'RunId must contain at least one letter.'}
# Study B frozen registry: 52 cells, all at R=500.
$expectedCells=52
$eR500=if($RepCap -gt 0){[Math]::Min(500,$RepCap)}else{500}
$registryExpected=[long](52L*$eR500)
function Get-ExpectedShardRows([int]$Shard,[int]$ShardCount) {
    if($RepCap -eq 1 -and $ShardCount -gt 1){return [long](Get-ExpectedShardCells $Shard $ShardCount)}
    $u5=[long][Math]::Floor(([double]$Shard*$eR500)/$ShardCount)
    $l5=[long][Math]::Floor(([double]($Shard-1)*$eR500)/$ShardCount)
    return [long]52*($u5-$l5)
}
function Get-ExpectedShardCells([int]$Shard,[int]$ShardCount) {
    if($RepCap -eq 1 -and $ShardCount -gt 1){
        if($Shard -gt $expectedCells){return 0}
        return [int]([Math]::Floor(([double]($expectedCells-$Shard))/$ShardCount)+1)
    }
    return $expectedCells
}
if($RepCap -eq 1 -and $NShard -gt $expectedCells){throw 'Parallel smoke requires NShard <= 52 so every shard owns at least one cell.'}
if($ExpectedReplications -eq 0){$ExpectedReplications=$registryExpected}
elseif($ExpectedReplications -ne $registryExpected){throw "ExpectedReplications=$ExpectedReplications disagrees with the frozen registry ($registryExpected for RepCap=$RepCap)."}
$src=Split-Path -Parent $MyInvocation.MyCommand.Path
$isStagedResume=$Resume -and (Test-Path -LiteralPath (Join-Path $src 'manifest.json') -PathType Leaf)
if($isStagedResume){$runs=Split-Path -Parent $src}
else{$runs=Join-Path $src 'runs'}
if(-not(Test-Path -LiteralPath $Stata -PathType Leaf)){throw "Stata not found: $Stata"}
$Stata=(Resolve-Path -LiteralPath $Stata).Path
if($Fresh -and [string]::IsNullOrWhiteSpace($RunId)){$RunId='study_b_{0}_M{1}' -f (Get-Date -Format 'yyyyMMdd_HHmmss'),$Master}
if($Resume -and [string]::IsNullOrWhiteSpace($RunId)){throw '-Resume requires -RunId.'}

$files=@('study_b_shard.do','study_b_worker.do','study_b_cells.do','_registry_probe.do','merge_study_b.do',
 'xtdpthresh.ado','xtdpthresh_p.ado','run_study_b.ps1',
 'monitor_study_b.ps1','verify_and_merge_study_b.ps1','README.md','_DESIGN_B.md')
$dir=if($isStagedResume){$src}else{Join-Path $runs $RunId}
$manifestPath=Join-Path $dir 'manifest.json'
$config=[ordered]@{RunId=$RunId;NShard=$NShard;Master=$Master;B=$B;Grid=$Grid;
 GridCI=$GridCI;RepCap=$RepCap;ExpectedReplications=$ExpectedReplications;
 ExpectedCells=$expectedCells;
 HarnessVersion=$HarnessVersion;CodeVersionExpected=$CodeVersionExpected;RngKind='mt64';
 MemoryGateSeconds=$MemoryGateSeconds;MinHeadroomMB=$MinHeadroomMB}

if($Fresh) {
    foreach($f in $files){if(-not(Test-Path -LiteralPath (Join-Path $src $f) -PathType Leaf)){throw "Required file missing: $f"}}
    if(Test-Path -LiteralPath $dir){throw "Run already exists; nothing deleted: $dir"}
    New-Item -ItemType Directory -Path $dir | Out-Null
    $hash=[ordered]@{}
    foreach($f in $files){$to=Join-Path $dir $f; Copy-Item -LiteralPath (Join-Path $src $f) -Destination $to; $hash[$f]=(Get-FileHash -LiteralPath $to -Algorithm SHA256).Hash}
    $manifest=[ordered]@{schemaVersion=1;createdUtc=(Get-Date).ToUniversalTime().ToString('o');sourceDirectory=$src;stataExecutable=$Stata;stataExecutableSha256=(Get-FileHash -LiteralPath $Stata -Algorithm SHA256).Hash;config=$config;sha256=$hash}
    Save-Json $manifest $manifestPath
} else {
    if(-not(Test-Path -LiteralPath $manifestPath -PathType Leaf)){throw "Manifest not found: $manifestPath"}
    $manifest=Get-Content -LiteralPath $manifestPath -Raw | ConvertFrom-Json
    if($Stata -ine [string]$manifest.stataExecutable){throw 'Resume Stata executable differs from the manifest.'}
    $selfHash=(Get-FileHash -LiteralPath $MyInvocation.MyCommand.Path -Algorithm SHA256).Hash
    if($selfHash -cne [string]$manifest.sha256.'run_study_b.ps1'){throw "Current launcher differs from the immutable staged launcher. Resume with: $dir\run_study_b.ps1"}
    foreach($key in $config.Keys){$old=[string]$manifest.config.$key;$new=[string]$config[$key];if($old -cne $new){throw "Resume config mismatch: $key (manifest='$old', requested='$new')"}}
    foreach($p in $manifest.sha256.PSObject.Properties){$path=Join-Path $dir $p.Name;if(-not(Test-Path -LiteralPath $path)){throw "Staged file missing: $($p.Name)"};if((Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash -cne [string]$p.Value){throw "Staged file changed: $($p.Name)"}}
    if([string]$manifest.stataExecutableSha256 -notmatch '^[A-Fa-f0-9]{64}$' -or
       (Get-FileHash -LiteralPath $Stata -Algorithm SHA256).Hash -cne [string]$manifest.stataExecutableSha256){throw 'Resume Stata executable hash differs from the manifest.'}
}

$lockPath=Join-Path $dir '.launcher.lock'
try {
    $launcherLock=New-Object IO.FileStream($lockPath,[IO.FileMode]::OpenOrCreate,
        [IO.FileAccess]::ReadWrite,[IO.FileShare]::None)
} catch {
    throw "Another launcher owns this run (lock: $lockPath)."
}
try {
# Execute the frozen registry itself under the launcher lock before any shard
# can start.  Hashing the probe file alone does not certify its realized cells.
Invoke-RegistryProbe $dir
if($StageOnly){
    Write-Host "Stage/probe complete; no shard launched: $dir"
    return
}
$wrapper=Join-Path $dir '_run_shard_wrapper.ps1'
@'
[CmdletBinding()]
param([string]$Stata,[string]$DoFile,[string]$LogFile,[string]$StatusFile,
      [string]$ResultFile,[int]$Shard,[string]$RunId,
      [long]$ExpectedRows,[int]$ExpectedCells,[long]$Master,[int]$NShard,
      [string]$HarnessVersion,[string]$CodeVersionExpected)
$ErrorActionPreference='Stop'
function Rows([string]$Path){if(-not(Test-Path -LiteralPath $Path -PathType Leaf)){return 0L};$r=[IO.File]::OpenText($Path);try{$n=0L;while($null-ne$r.ReadLine()){$n++};return [Math]::Max(0L,$n-1L)}finally{$r.Dispose()}}
function M([string]$TemporaryPath,[string]$DestinationPath){for($attempt=1;$attempt-le10;$attempt++){try{Move-Item -LiteralPath $TemporaryPath -Destination $DestinationPath -Force;return}catch{if($attempt-eq10){throw};Start-Sleep -Milliseconds (50*$attempt)}}}
function S([string]$state,[int]$stataPid,$exitCode,[string]$message){$rows=0L;try{$rows=Rows $ResultFile}catch{};$o=[ordered]@{schemaVersion=1;runId=$RunId;harnessVersion=$HarnessVersion;codeVersionExpected=$CodeVersionExpected;master=$Master;shard=$Shard;nshard=$NShard;state=$state;wrapperPid=$PID;stataPid=$stataPid;updatedUtc=(Get-Date).ToUniversalTime().ToString('o');exitCode=$exitCode;expectedRows=$ExpectedRows;expectedCells=$ExpectedCells;resultRows=$rows;message=$message};$t=$StatusFile+'.tmp.'+[Guid]::NewGuid().ToString('N');$o|ConvertTo-Json|Set-Content -LiteralPath $t -Encoding UTF8;M $t $StatusFile}
$p=$null
try{S 'launching' 0 $null 'Starting Stata';$a='/e do "'+$DoFile.Replace('"','\"')+'"';$p=Start-Process -FilePath $Stata -ArgumentList $a -WorkingDirectory (Split-Path $DoFile) -WindowStyle Hidden -PassThru;S 'running' $p.Id $null 'Stata running';$p.WaitForExit();$marker=Join-Path (Split-Path $DoFile) ('study_b_complete_SH'+$Shard+'.csv');$ok=$false;if(Test-Path -LiteralPath $marker -PathType Leaf){$cm=@(Import-Csv -LiteralPath $marker);$rows=Rows $ResultFile;$ok=($cm.Count-eq1-and[string]$cm[0].run_id-ceq$RunId-and[string]$cm[0].harness_version-ceq$HarnessVersion-and[string]$cm[0].code_version_expected-ceq$CodeVersionExpected-and[long]$cm[0].master-eq$Master-and[int]$cm[0].shard-eq$Shard-and[int]$cm[0].nshard-eq$NShard-and[int]$cm[0].cells-eq$ExpectedCells-and[long]$cm[0].reps-eq$ExpectedRows-and[int]$cm[0].exit_code-eq0-and$rows-eq$ExpectedRows)};if($p.ExitCode -eq 0 -and $ok){S 'completed' 0 $p.ExitCode 'Completion file and row count verified'}else{S 'failed' 0 $p.ExitCode 'Nonzero exit or invalid completion/count contract';exit 1}}catch{$childPid=0;if($null-ne$p){$childPid=$p.Id;try{if(-not$p.HasExited){Stop-Process -Id $p.Id -Force;$p.WaitForExit()}}catch{}};try{S 'failed' 0 -1 $_.Exception.Message}catch{};exit 1}
'@ | Set-Content -LiteralPath $wrapper -Encoding UTF8

$ps=Join-Path $PSHome 'powershell.exe';$launched=0;$done=0;$active=0;$stataDir=$dir.Replace('"','""')
for($k=1;$k -le $NShard;$k++) {
    $tag='{0:D3}' -f $k;$status=Join-Path $dir "status_SH$tag.json";$log=Join-Path $dir "shard_SH$tag.log";$do=Join-Path $dir "launch_SH$tag.do"
    $expectedShardRows=Get-ExpectedShardRows $k $NShard
    $expectedShardCells=Get-ExpectedShardCells $k $NShard
    if(Test-Path -LiteralPath $status){
        $s=Get-Content -LiteralPath $status -Raw|ConvertFrom-Json
        if([int]$s.schemaVersion -ne 1 -or [string]$s.runId -cne $RunId -or
           [string]$s.harnessVersion -cne $HarnessVersion -or
           [string]$s.codeVersionExpected -cne $CodeVersionExpected -or
           [long]$s.master -ne $Master -or [int]$s.shard -ne $k -or
           [int]$s.nshard -ne $NShard -or [long]$s.expectedRows -ne $expectedShardRows -or
           [int]$s.expectedCells -ne $expectedShardCells){throw "Status identity/config mismatch: $status"}
        if($s.state -eq 'completed'){
            $marker=Join-Path $dir "study_b_complete_SH$k.csv"
            $result=Join-Path $dir "study_b_results_SH$k.csv"
            if(-not(Test-Path -LiteralPath $marker -PathType Leaf)-or-not(Test-Path -LiteralPath $result -PathType Leaf)){throw "Completed shard $k is missing its marker or result CSV."}
            $cm=@(Import-Csv -LiteralPath $marker)
            if($cm.Count -ne 1 -or [string]$cm[0].run_id -cne $RunId -or [string]$cm[0].harness_version -cne $HarnessVersion -or [string]$cm[0].code_version_expected -cne $CodeVersionExpected -or [long]$cm[0].master -ne $Master -or [int]$cm[0].shard -ne $k -or [int]$cm[0].nshard -ne $NShard -or [int]$cm[0].cells -ne $expectedShardCells -or [long]$cm[0].reps -ne $expectedShardRows -or [int]$cm[0].exit_code -ne 0 -or [long]$s.expectedRows -ne $expectedShardRows -or [long]$s.resultRows -ne $expectedShardRows -or (Get-ResultRows $result) -ne $expectedShardRows){throw "Invalid completion/count contract for shard $k."}
            $done++;continue
        }
        if(Is-Alive $s){$active++;continue}
    }
    if(($launched+$active) -gt 0){Wait-MemoryHeadroom $k}
    $stataLog=$log.Replace('"','""')
@"
version 15.0
clear all
set more off
cd "$stataDir"
capture log close _all
log using "$stataLog", text append name(mcbrun)
capture noisily do study_b_shard.do $k $NShard $Master $B $Grid $GridCI $RepCap "$RunId" "$HarnessVersion" "$CodeVersionExpected" "mt64"
local runrc = _rc
log close mcbrun
exit ``runrc', clear
"@ | Set-Content -LiteralPath $do -Encoding ASCII
    $result=Join-Path $dir "study_b_results_SH$k.csv"
    $marker=Join-Path $dir "study_b_complete_SH$k.csv"
    Remove-Item -LiteralPath $marker -Force -ErrorAction SilentlyContinue
    $q=[ordered]@{schemaVersion=1;runId=$RunId;harnessVersion=$HarnessVersion;codeVersionExpected=$CodeVersionExpected;master=$Master;shard=$k;nshard=$NShard;state='queued';wrapperPid=0;stataPid=0;updatedUtc=(Get-Date).ToUniversalTime().ToString('o');exitCode=$null;expectedRows=$expectedShardRows;expectedCells=$expectedShardCells;resultRows=(Get-ResultRows $result);message='queued'};Save-Json $q $status
    $args=@('-NoProfile','-ExecutionPolicy','Bypass','-File',(Q $wrapper),'-Stata',(Q $Stata),'-DoFile',(Q $do),'-LogFile',(Q $log),'-StatusFile',(Q $status),'-ResultFile',(Q $result),'-Shard',[string]$k,'-RunId',(Q $RunId),'-ExpectedRows',[string]$expectedShardRows,'-ExpectedCells',[string]$expectedShardCells,'-Master',[string]$Master,'-NShard',[string]$NShard,'-HarnessVersion',(Q $HarnessVersion),'-CodeVersionExpected',(Q $CodeVersionExpected))-join ' '
    $wp=Start-Process -FilePath $ps -ArgumentList $args -WorkingDirectory $dir -WindowStyle Hidden -PassThru
    $ready=$false;$deadline=[DateTime]::UtcNow.AddSeconds(10)
    do {
        Start-Sleep -Milliseconds 50
        try {$cur=Get-Content -LiteralPath $status -Raw|ConvertFrom-Json;if([int]$cur.wrapperPid-eq$wp.Id-and$cur.state-ne'queued'){$ready=$true}} catch {}
        if($wp.HasExited -and-not$ready){break}
    } while(-not$ready-and[DateTime]::UtcNow-lt$deadline)
    if(-not$ready){try{if(-not$wp.HasExited){Stop-Process -Id $wp.Id -Force}}catch{};$q.state='failed';$q.wrapperPid=0;$q.exitCode=-1;$q.updatedUtc=(Get-Date).ToUniversalTime().ToString('o');$q.message='Wrapper did not claim status within 10 seconds';Save-Json $q $status;throw "Shard wrapper $k did not claim its status within 10 seconds."}
    $launched++
}
Write-Host "Run ID: $RunId";Write-Host "Run directory: $dir";Write-Host "Shards launched: $launched; already active: $active; already complete: $done";Write-Host "Monitor: .\monitor_study_b.ps1 -RunId $RunId -Watch 60"
} finally {
    if($null -ne $launcherLock){$launcherLock.Dispose()}
}

# run_study_a.ps1 -- isolated, resumable Study A launcher (PowerShell 5.1).
[CmdletBinding()]
param(
    [switch]$Fresh,
    [switch]$Resume,
    [ValidatePattern('^[A-Za-z0-9][A-Za-z0-9_.-]{0,79}$')]
    [string]$RunId,
    [ValidateRange(1,256)]
    [int]$NShard = 28,
    [ValidateRange(1,2147483000)]
    [long]$Master = 20260814,
    [ValidateRange(0,500)]
    [int]$RepCap = 0,
    [ValidatePattern('^[A-Za-z0-9][A-Za-z0-9_.-]*$')]
    [string]$HarnessVersion = 'xtdpthresh_study_a_v0937_r400',
    [ValidatePattern('^[0-9]+\.[0-9]+\.[0-9]+$')]
    [string]$CodeVersionExpected = '0.9.37',
    [string]$Stata = 'C:\Program Files\Stata17\StataMP-64.exe'
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version 2.0

$expectedCells = 116
$rngKind = 'mt64'
$requiredFiles = @(
    'study_a_shard.do',
    'study_a_worker.do',
    'study_a_cells.do',
    '_registry_probe.do',
    'merge_study_a.do',
    'xtdpthresh.ado',
    'xtdpthresh_p.ado',
    'run_study_a.ps1',
    'monitor_study_a.ps1',
    'verify_and_merge_study_a.ps1',
    'README.md',
    '_DESIGN.md'
)

function Save-JsonAtomic {
    param($Object, [string]$Path)
    $temporaryPath = $Path + '.tmp.' + [Guid]::NewGuid().ToString('N')
    $Object | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $temporaryPath -Encoding UTF8
    for ($attempt = 1; $attempt -le 10; $attempt++) {
        try {
            Move-Item -LiteralPath $temporaryPath -Destination $Path -Force
            return
        }
        catch {
            if ($attempt -eq 10) {throw}
            Start-Sleep -Milliseconds (25 * $attempt)
        }
    }
}

function Quote-Argument {
    param([string]$Value)
    return '"' + $Value.Replace('"', '\"') + '"'
}

function Get-ResultRows {
    param([string]$Path)
    $reader = [IO.File]::OpenText($Path)
    try {
        $lineCount = 0L
        while ($null -ne $reader.ReadLine()) {
            $lineCount++
        }
        return [Math]::Max(0L, $lineCount - 1L)
    }
    finally {
        $reader.Dispose()
    }
}

function Get-ExpectedShardRows {
    param([int]$Shard, [int]$ShardCount)
    if ($RepCap -eq 1 -and $ShardCount -gt 1) {
        return [long](Get-ExpectedShardCells -Shard $Shard -ShardCount $ShardCount)
    }
    # Two replication groups: 24 stresses at R=500, 92 regular cells at R=400.
    # Each cell's shard split uses its own R, so a shard's row
    # count is the sum over the two groups.
    $u5 = [long][Math]::Floor(([double]$Shard * $script:eR500) / $ShardCount)
    $l5 = [long][Math]::Floor(([double]($Shard - 1) * $script:eR500) / $ShardCount)
    $u4 = [long][Math]::Floor(([double]$Shard * $script:eR400) / $ShardCount)
    $l4 = [long][Math]::Floor(([double]($Shard - 1) * $script:eR400) / $ShardCount)
    return [long]$script:cells500 * ($u5 - $l5) + `
        [long]$script:cells400 * ($u4 - $l4)
}

function Get-ExpectedShardCells {
    param([int]$Shard, [int]$ShardCount)
    if ($RepCap -eq 1 -and $ShardCount -gt 1) {
        if ($Shard -gt $expectedCells) {return 0}
        return [int]([Math]::Floor(([double]($expectedCells - $Shard)) / $ShardCount) + 1)
    }
    return $expectedCells
}

function Test-ProcessIdentity {
    param($Status)
    foreach ($pidField in @('wrapperPid', 'stataPid')) {
        if ($Status.PSObject.Properties.Name -notcontains $pidField) {
            continue
        }
        $processId = [int]$Status.$pidField
        if ($processId -le 0) {
            continue
        }
        $process = Get-Process -Id $processId -ErrorAction SilentlyContinue
        if ($null -eq $process) {
            continue
        }
        if ($pidField -eq 'wrapperPid' -and $process.ProcessName -notlike '*powershell*') {
            continue
        }
        if ($pidField -eq 'stataPid') {
            try {
                if ($process.Path -ine $Stata) {
                    continue
                }
            }
            catch {
                continue
            }
        }
        try {
            $statusTime = [DateTime]::Parse([string]$Status.updatedUtc).ToUniversalTime()
            if ($process.StartTime.ToUniversalTime() -gt $statusTime.AddSeconds(10)) {
                continue
            }
        }
        catch {
            continue
        }
        return $true
    }
    return $false
}

function Assert-Completion {
    param(
        [string]$RunDirectory,
        [int]$Shard,
        [long]$ExpectedRows,
        [int]$ExpectedCells
    )
    $markerPath = Join-Path $RunDirectory ('study_a_complete_SH{0}.csv' -f $Shard)
    $resultPath = Join-Path $RunDirectory ('study_a_results_SH{0}.csv' -f $Shard)
    $progressPath = Join-Path $RunDirectory ('study_a_progress_SH{0}.csv' -f $Shard)
    if (-not (Test-Path -LiteralPath $markerPath -PathType Leaf)) {
        throw "Completed shard $Shard is missing its completion marker."
    }
    if (-not (Test-Path -LiteralPath $resultPath -PathType Leaf)) {
        throw "Completed shard $Shard is missing its result CSV."
    }
    if (-not (Test-Path -LiteralPath $progressPath -PathType Leaf)) {
        throw "Completed shard $Shard is missing its progress marker."
    }
    $marker = @(Import-Csv -LiteralPath $markerPath)
    if ($marker.Count -ne 1) {
        throw "Completion marker for shard $Shard must contain exactly one record."
    }
    $record = $marker[0]
    if ([string]$record.run_id -cne $RunId -or
        [string]$record.harness_version -cne $HarnessVersion -or
        [string]$record.code_version_expected -cne $CodeVersionExpected -or
        [long]$record.master -ne $Master -or
        [int]$record.shard -ne $Shard -or
        [int]$record.nshard -ne $NShard -or
        [int]$record.cells -ne $ExpectedCells -or
        [long]$record.reps -ne $ExpectedRows -or
        [int]$record.exit_code -ne 0) {
        throw "Completion marker identity/count mismatch for shard $Shard."
    }
    $actualRows = Get-ResultRows -Path $resultPath
    if ($actualRows -ne $ExpectedRows) {
        throw "Result row count mismatch for completed shard ${Shard}: expected $ExpectedRows, found $actualRows."
    }
    $progress = @(Import-Csv -LiteralPath $progressPath)
    if ($progress.Count -ne 1) {
        throw "Progress marker for shard $Shard must contain exactly one record."
    }
    $p = $progress[0]
    if ([string]$p.run_id -cne $RunId -or
        [string]$p.harness_version -cne $HarnessVersion -or
        [string]$p.code_version_expected -cne $CodeVersionExpected -or
        [long]$p.master -ne $Master -or
        [int]$p.shard -ne $Shard -or
        [int]$p.nshard -ne $NShard -or
        [int]$p.cells_completed -ne $ExpectedCells -or
        [long]$p.reps_completed -ne $ExpectedRows) {
        throw "Progress marker identity/count mismatch for shard $Shard."
    }
}

if (($Fresh -and $Resume) -or (-not $Fresh -and -not $Resume)) {
    throw 'Choose exactly one launch mode: -Fresh or -Resume.'
}
if ($HarnessVersion -cne 'xtdpthresh_study_a_v0937_r400') {
    throw 'Study A is frozen at HarnessVersion=xtdpthresh_study_a_v0937_r400.'
}
if ($CodeVersionExpected -cne '0.9.37') {
    throw 'Study A is frozen at CodeVersionExpected=0.9.37.'
}
if (-not [string]::IsNullOrWhiteSpace($RunId) -and $RunId -notmatch '[A-Za-z]') {
    throw 'RunId must contain at least one letter.'
}
if ($RepCap -gt 0 -and $RepCap -lt $NShard -and $RepCap -ne 1) {
    throw 'A capped run requires RepCap >= NShard; RepCap=1 is the dedicated parallel-smoke exception.'
}
if ($RepCap -eq 1 -and $NShard -gt $expectedCells) {
    throw 'Parallel smoke requires NShard <= 116 so every shard owns at least one cell.'
}

$sourceDirectory = Split-Path -Parent $MyInvocation.MyCommand.Path
$sourceHasManifest = Test-Path -LiteralPath (Join-Path $sourceDirectory 'manifest.json') -PathType Leaf
if ($Fresh -and $sourceHasManifest) {
    throw 'A fresh run must be launched from the package root, not from an immutable run snapshot.'
}
$isStagedResume = $Resume -and $sourceHasManifest
$runsDirectory = if ($isStagedResume) {
    Split-Path -Parent $sourceDirectory
}
else {
    Join-Path $sourceDirectory 'runs'
}

if (-not (Test-Path -LiteralPath $Stata -PathType Leaf)) {
    throw "Stata executable not found: $Stata"
}
$Stata = (Resolve-Path -LiteralPath $Stata).Path

if ($Fresh -and [string]::IsNullOrWhiteSpace($RunId)) {
    $RunId = 'study_a_{0}_M{1}' -f (Get-Date -Format 'yyyyMMdd_HHmmss'), $Master
}
if ($Resume -and [string]::IsNullOrWhiteSpace($RunId)) {
    throw '-Resume requires -RunId.'
}

# Study A: 24 cells @ R=500, 92 cells @ R=400 (including N=1600).
# RepCap caps each group's effective R.
$script:cells500 = 24
$script:cells400 = 92
$script:eR500 = if ($RepCap -gt 0) { [Math]::Min(500, $RepCap) } else { 500 }
$script:eR400 = if ($RepCap -gt 0) { [Math]::Min(400, $RepCap) } else { 400 }
$expectedReplications = [long]$script:cells500 * $script:eR500 + `
    [long]$script:cells400 * $script:eR400
$runDirectory = if ($isStagedResume) {
    $sourceDirectory
}
else {
    Join-Path $runsDirectory $RunId
}
$manifestPath = Join-Path $runDirectory 'manifest.json'
$config = [ordered]@{
    RunId                 = $RunId
    NShard                = $NShard
    Master                = $Master
    RepCap                = $RepCap
    ExpectedCells         = $expectedCells
    ExpectedReplications  = $expectedReplications
    HarnessVersion        = $HarnessVersion
    CodeVersionExpected   = $CodeVersionExpected
    RngKind               = $rngKind
}

if ($Fresh) {
    foreach ($fileName in $requiredFiles) {
        $sourcePath = Join-Path $sourceDirectory $fileName
        if (-not (Test-Path -LiteralPath $sourcePath -PathType Leaf)) {
            throw "Required package file is missing: $fileName"
        }
    }
    if (-not (Test-Path -LiteralPath $runsDirectory -PathType Container)) {
        New-Item -ItemType Directory -Path $runsDirectory | Out-Null
    }
    if (Test-Path -LiteralPath $runDirectory) {
        throw "Run already exists; nothing was overwritten: $runDirectory"
    }
    New-Item -ItemType Directory -Path $runDirectory | Out-Null
    $hashes = [ordered]@{}
    foreach ($fileName in $requiredFiles) {
        $sourcePath = Join-Path $sourceDirectory $fileName
        $stagedPath = Join-Path $runDirectory $fileName
        Copy-Item -LiteralPath $sourcePath -Destination $stagedPath
        $hashes[$fileName] = (Get-FileHash -LiteralPath $stagedPath -Algorithm SHA256).Hash
    }
    $manifest = [ordered]@{
        schemaVersion   = 1
        createdUtc      = (Get-Date).ToUniversalTime().ToString('o')
        sourceDirectory = $sourceDirectory
        stataExecutable = $Stata
        config          = $config
        sha256          = $hashes
    }
    Save-JsonAtomic -Object $manifest -Path $manifestPath
    # Reload so Fresh and Resume use the same PSCustomObject representation.
    # PowerShell 5.1 exposes OrderedDictionary keys differently via PSObject.
    $manifest = Get-Content -LiteralPath $manifestPath -Raw | ConvertFrom-Json
}
else {
    if (-not (Test-Path -LiteralPath $manifestPath -PathType Leaf)) {
        throw "Manifest not found: $manifestPath"
    }
    $manifest = Get-Content -LiteralPath $manifestPath -Raw | ConvertFrom-Json
    if ([int]$manifest.schemaVersion -ne 1) {
        throw 'Unsupported or missing manifest schemaVersion.'
    }
    if ($Stata -ine [string]$manifest.stataExecutable) {
        throw 'Resume Stata executable differs from the immutable manifest.'
    }
    $stagedLauncherHash = (Get-FileHash -LiteralPath $MyInvocation.MyCommand.Path -Algorithm SHA256).Hash
    if ($stagedLauncherHash -cne [string]$manifest.sha256.'run_study_a.ps1') {
        throw "Launcher differs from the staged copy. Resume with: $runDirectory\run_study_a.ps1"
    }
    foreach ($key in $config.Keys) {
        $oldValue = [string]$manifest.config.$key
        $newValue = [string]$config[$key]
        if ($oldValue -cne $newValue) {
            throw "Resume config mismatch: $key (manifest='$oldValue', requested='$newValue')"
        }
    }
}

# Validate every frozen file before every launch, including the initial launch.
if ($null -eq $manifest.sha256 -or @($manifest.sha256.PSObject.Properties).Count -ne $requiredFiles.Count) {
    throw 'Manifest SHA256 inventory does not match the frozen package inventory.'
}
foreach ($fileName in $requiredFiles) {
    $hashProperty = $manifest.sha256.PSObject.Properties[$fileName]
    if ($null -eq $hashProperty -or [string]$hashProperty.Value -notmatch '^[A-Fa-f0-9]{64}$') {
        throw "Manifest has no valid SHA256 value for: $fileName"
    }
    $stagedPath = Join-Path $runDirectory $fileName
    if (-not (Test-Path -LiteralPath $stagedPath -PathType Leaf)) {
        throw "Frozen file is missing: $fileName"
    }
    $actualHash = (Get-FileHash -LiteralPath $stagedPath -Algorithm SHA256).Hash
    if ($actualHash -cne [string]$hashProperty.Value) {
        throw "Frozen file changed after staging: $fileName"
    }
}

$lockPath = Join-Path $runDirectory '.launcher.lock'
$launcherLock = $null
try {
    $launcherLock = New-Object IO.FileStream(
        $lockPath,
        [IO.FileMode]::OpenOrCreate,
        [IO.FileAccess]::ReadWrite,
        [IO.FileShare]::None
    )
}
catch {
    throw "Another launcher owns this run (lock: $lockPath)."
}

try {
    $activeStatuses = @()
    Get-ChildItem -LiteralPath $runDirectory -Filter 'status_SH*.json' -ErrorAction SilentlyContinue | ForEach-Object {
        try {
            $statusObject = Get-Content -LiteralPath $_.FullName -Raw | ConvertFrom-Json
            if ([string]$statusObject.runId -cne $RunId -or
                [int]$statusObject.shard -lt 1 -or
                [int]$statusObject.shard -gt $NShard) {
                throw 'status identity mismatch'
            }
            if (Test-ProcessIdentity -Status $statusObject) {
                $activeStatuses += $_.Name
            }
        }
        catch {
            throw "Invalid shard status file: $($_.FullName)"
        }
    }
    if ($activeStatuses.Count -gt 0) {
        throw ('This run still has active processes ({0}); nothing was launched.' -f ($activeStatuses -join ', '))
    }

    $wrapperPath = Join-Path $runDirectory '_run_study_a_shard_wrapper.ps1'
    @'
[CmdletBinding()]
param(
    [string]$Stata,
    [string]$DoFile,
    [string]$StatusFile,
    [int]$Shard,
    [int]$NShard,
    [long]$Master,
    [long]$ExpectedRows,
    [int]$ExpectedCells,
    [string]$RunId,
    [string]$HarnessVersion,
    [string]$CodeVersionExpected
)
$ErrorActionPreference = 'Stop'
function Save-Status([string]$State, [int]$StataPid, $ExitCode, [string]$Message, [long]$ResultRows) {
    $value = [ordered]@{
        schemaVersion = 1; runId = $RunId; shard = $Shard; state = $State
        wrapperPid = $PID; stataPid = $StataPid
        updatedUtc = (Get-Date).ToUniversalTime().ToString('o')
        exitCode = $ExitCode; expectedRows = $ExpectedRows; resultRows = $ResultRows
        message = $Message
    }
    $temporaryPath = $StatusFile + '.tmp.' + [Guid]::NewGuid().ToString('N')
    $value | ConvertTo-Json | Set-Content -LiteralPath $temporaryPath -Encoding UTF8
    for ($attempt = 1; $attempt -le 10; $attempt++) {
        try {
            Move-Item -LiteralPath $temporaryPath -Destination $StatusFile -Force
            return
        }
        catch {
            if ($attempt -eq 10) {throw}
            Start-Sleep -Milliseconds (25 * $attempt)
        }
    }
}
function Count-Rows([string]$Path) {
    $reader = [IO.File]::OpenText($Path)
    try {$count = 0L; while ($null -ne $reader.ReadLine()) {$count++}; return [Math]::Max(0L, $count - 1L)}
    finally {$reader.Dispose()}
}
function Quote-Argument([string]$Value) {return '"' + $Value.Replace('"', '\"') + '"'}
$stataProcess = $null
try {
    Save-Status 'launching' 0 $null 'Starting Stata' 0
    $stataArguments = '/e do ' + (Quote-Argument $DoFile)
    $stataProcess = Start-Process -FilePath $Stata -ArgumentList $stataArguments `
        -WorkingDirectory (Split-Path -Parent $DoFile) -WindowStyle Hidden -PassThru
    Save-Status 'running' $stataProcess.Id $null 'Stata running' 0
    $stataProcess.WaitForExit()
    $directory = Split-Path -Parent $DoFile
    $markerPath = Join-Path $directory ('study_a_complete_SH' + $Shard + '.csv')
    $resultPath = Join-Path $directory ('study_a_results_SH' + $Shard + '.csv')
    $progressPath = Join-Path $directory ('study_a_progress_SH' + $Shard + '.csv')
    $valid = $false
    $rowCount = -1L
    if ((Test-Path -LiteralPath $markerPath -PathType Leaf) -and
        (Test-Path -LiteralPath $resultPath -PathType Leaf) -and
        (Test-Path -LiteralPath $progressPath -PathType Leaf)) {
        $marker = @(Import-Csv -LiteralPath $markerPath)
        $progress = @(Import-Csv -LiteralPath $progressPath)
        $rowCount = Count-Rows $resultPath
        if ($marker.Count -eq 1 -and $progress.Count -eq 1) {
            $record = $marker[0]
            $p = $progress[0]
            $valid = ([string]$record.run_id -ceq $RunId -and
                [string]$record.harness_version -ceq $HarnessVersion -and
                [string]$record.code_version_expected -ceq $CodeVersionExpected -and
                [long]$record.master -eq $Master -and
                [int]$record.shard -eq $Shard -and
                [int]$record.nshard -eq $NShard -and
                [int]$record.cells -eq $ExpectedCells -and
                [long]$record.reps -eq $ExpectedRows -and
                [int]$record.exit_code -eq 0 -and
                [string]$p.run_id -ceq $RunId -and
                [string]$p.harness_version -ceq $HarnessVersion -and
                [string]$p.code_version_expected -ceq $CodeVersionExpected -and
                [long]$p.master -eq $Master -and
                [int]$p.shard -eq $Shard -and
                [int]$p.nshard -eq $NShard -and
                [int]$p.cells_completed -eq $ExpectedCells -and
                [long]$p.reps_completed -eq $ExpectedRows -and
                $rowCount -eq $ExpectedRows)
        }
    }
    if ($stataProcess.ExitCode -eq 0 -and $valid) {
        Save-Status 'completed' $stataProcess.Id $stataProcess.ExitCode 'Completion marker and row count verified' $rowCount
    }
    else {
        Save-Status 'failed' $stataProcess.Id $stataProcess.ExitCode 'Nonzero exit or invalid completion artifacts' $rowCount
        exit 1
    }
}
catch {
    $childPid = 0
    if ($null -ne $stataProcess) {
        $childPid = $stataProcess.Id
        try {
            if (-not $stataProcess.HasExited) {
                Stop-Process -Id $stataProcess.Id -Force
                $stataProcess.WaitForExit()
            }
        }
        catch {}
    }
    try {Save-Status 'failed' $childPid -1 $_.Exception.Message -1} catch {}
    exit 1
}
'@ | Set-Content -LiteralPath $wrapperPath -Encoding UTF8

    $powerShellExe = Join-Path $PSHome 'powershell.exe'
    $stataDirectory = $runDirectory.Replace('"', '""')
    $launchedCount = 0
    $completedCount = 0

    for ($shard = 1; $shard -le $NShard; $shard++) {
        $tag = '{0:D3}' -f $shard
        $statusPath = Join-Path $runDirectory "status_SH$tag.json"
        $logPath = Join-Path $runDirectory "study_a_shard_SH$tag.log"
        $doPath = Join-Path $runDirectory "launch_SH$tag.do"
        $expectedShardRows = Get-ExpectedShardRows -Shard $shard -ShardCount $NShard
        $expectedShardCells = Get-ExpectedShardCells -Shard $shard -ShardCount $NShard

        if (Test-Path -LiteralPath $statusPath -PathType Leaf) {
            $statusObject = Get-Content -LiteralPath $statusPath -Raw | ConvertFrom-Json
            if ([string]$statusObject.runId -cne $RunId -or [int]$statusObject.shard -ne $shard) {
                throw "Status identity mismatch: $statusPath"
            }
            if ([string]$statusObject.state -ceq 'completed') {
                Assert-Completion -RunDirectory $runDirectory -Shard $shard `
                    -ExpectedRows $expectedShardRows -ExpectedCells $expectedShardCells
                $completedCount++
                continue
            }
        }

        $stataLog = $logPath.Replace('"', '""')
        @"
version 15.0
clear all
set more off
set processors 1
cd "$stataDirectory"
capture log close _all
log using "$stataLog", text append name(studyarun)
capture noisily do study_a_shard.do $shard $NShard $Master $RepCap "$RunId" "$HarnessVersion" "$CodeVersionExpected" "$rngKind"
local runrc = _rc
log close studyarun
exit ``runrc', clear
"@ | Set-Content -LiteralPath $doPath -Encoding ASCII

        $queuedStatus = [ordered]@{
            schemaVersion = 1; runId = $RunId; shard = $shard; state = 'queued'
            wrapperPid = 0; stataPid = 0
            updatedUtc = (Get-Date).ToUniversalTime().ToString('o')
            exitCode = $null; expectedRows = $expectedShardRows; resultRows = 0
            message = 'Queued by launcher'
        }
        Save-JsonAtomic -Object $queuedStatus -Path $statusPath

        $wrapperArguments = @(
            '-NoLogo', '-NoProfile', '-ExecutionPolicy', 'Bypass',
            '-File', (Quote-Argument $wrapperPath),
            '-Stata', (Quote-Argument $Stata),
            '-DoFile', (Quote-Argument $doPath),
            '-StatusFile', (Quote-Argument $statusPath),
            '-Shard', [string]$shard,
            '-NShard', [string]$NShard,
            '-Master', [string]$Master,
            '-ExpectedRows', [string]$expectedShardRows,
            '-ExpectedCells', [string]$expectedShardCells,
            '-RunId', (Quote-Argument $RunId),
            '-HarnessVersion', (Quote-Argument $HarnessVersion),
            '-CodeVersionExpected', (Quote-Argument $CodeVersionExpected)
        ) -join ' '
        $wrapperProcess = Start-Process -FilePath $powerShellExe -ArgumentList $wrapperArguments `
            -WorkingDirectory $runDirectory -WindowStyle Hidden -PassThru

        $claimed = $false
        $deadline = [DateTime]::UtcNow.AddSeconds(15)
        do {
            Start-Sleep -Milliseconds 50
            try {
                $currentStatus = Get-Content -LiteralPath $statusPath -Raw | ConvertFrom-Json
                if ([int]$currentStatus.wrapperPid -eq $wrapperProcess.Id -and
                    [string]$currentStatus.state -cne 'queued') {
                    $claimed = $true
                }
            }
            catch {}
            if ($wrapperProcess.HasExited -and -not $claimed) {
                break
            }
        } while (-not $claimed -and [DateTime]::UtcNow -lt $deadline)

        if (-not $claimed) {
            try {
                if (-not $wrapperProcess.HasExited) {
                    Stop-Process -Id $wrapperProcess.Id -Force
                }
            }
            catch {}
            throw "Shard wrapper $shard did not claim its status within 15 seconds."
        }
        $launchedCount++

        # Launch stagger + adaptive RAM gate.  A shard "claims" its status as
        # soon as the wrapper flips it to 'running', which is BEFORE Stata has
        # finished loading the ado/Mata stack.  Launching all shards back to
        # back therefore piles every Stata initialization on top of the others;
        # each init transiently allocates far more than the ~175 MB steady state,
        # so on a small-RAM host the peak overflows physical memory and the OS
        # kills the late shards (access violation / failed module load).  Spacing
        # out the launches, and pausing while free memory is low, flattens that
        # peak so every requested core is actually used.  The added wall-clock is
        # a few minutes at most and is negligible against a multi-hour run.
        if ($shard -lt $NShard) {
            Start-Sleep -Seconds 3
            $ramDeadline = [DateTime]::UtcNow.AddSeconds(45)
            while ([DateTime]::UtcNow -lt $ramDeadline) {
                try {
                    $freeMb = [int]((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory / 1KB)
                }
                catch { break }
                if ($freeMb -ge 1536) { break }
                Start-Sleep -Seconds 2
            }
        }
    }

    Write-Host "Run ID: $RunId"
    Write-Host "Run directory: $runDirectory"
    Write-Host "Design: $expectedCells cells; $expectedReplications expected rows; point estimation only"
    Write-Host "Shards launched: $launchedCount; already complete: $completedCount"
    Write-Host "Monitor: .\monitor_study_a.ps1 -RunId $RunId -Watch 60"
}
finally {
    if ($null -ne $launcherLock) {
        $launcherLock.Dispose()
    }
}

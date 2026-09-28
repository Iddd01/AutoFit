# monitor_study_a.ps1 -- monitor exactly one isolated Study A run (PowerShell 5.1).
[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [ValidatePattern('^[A-Za-z0-9][A-Za-z0-9_.-]{0,79}$')]
    [string]$RunId,
    [ValidateRange(0,86400)]
    [int]$Watch = 0
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version 2.0

$expectedHarness = 'xtdpthresh_study_a_v0937_r400'
$expectedCodeVersion = '0.9.37'
$expectedRng = 'mt64'

$requiredFiles = @(
    'study_a_shard.do', 'study_a_worker.do', 'study_a_cells.do', '_registry_probe.do',
    'merge_study_a.do',
    'xtdpthresh.ado', 'xtdpthresh_p.ado',
    'run_study_a.ps1', 'monitor_study_a.ps1', 'verify_and_merge_study_a.ps1',
    'README.md', '_DESIGN.md'
)

function Get-ResultRows {
    param([string]$Path)
    $reader = [IO.File]::OpenText($Path)
    try {
        $lineCount = 0L
        while ($null -ne $reader.ReadLine()) {$lineCount++}
        return [Math]::Max(0L, $lineCount - 1L)
    }
    finally {$reader.Dispose()}
}

function Test-AliveProcess {
    param($Status, [string]$PidField, [string]$ExpectedStata)
    if ($Status.PSObject.Properties.Name -notcontains $PidField) {return $false}
    $processId = [int]$Status.$PidField
    if ($processId -le 0) {return $false}
    $process = Get-Process -Id $processId -ErrorAction SilentlyContinue
    if ($null -eq $process) {return $false}
    if ($PidField -eq 'wrapperPid' -and $process.ProcessName -notlike '*powershell*') {return $false}
    if ($PidField -eq 'stataPid') {
        try {if ($process.Path -ine $ExpectedStata) {return $false}}
        catch {return $false}
    }
    try {
        $statusTime = [DateTime]::Parse([string]$Status.updatedUtc).ToUniversalTime()
        if ($process.StartTime.ToUniversalTime() -gt $statusTime.AddSeconds(10)) {return $false}
    }
    catch {return $false}
    return $true
}

function Get-ExpectedShardRows {
    param([int]$Shard, [int]$ShardCount)
    if ($script:repCap -eq 1 -and $ShardCount -gt 1) {
        return [long](Get-ExpectedShardCells -Shard $Shard -ShardCount $ShardCount)
    }
    # Two replication groups: 24 stresses at R=500, 92 regular cells at R=400.
    $u5 = [long][Math]::Floor(([double]$Shard * $script:eR500) / $ShardCount)
    $l5 = [long][Math]::Floor(([double]($Shard - 1) * $script:eR500) / $ShardCount)
    $u4 = [long][Math]::Floor(([double]$Shard * $script:eR400) / $ShardCount)
    $l4 = [long][Math]::Floor(([double]($Shard - 1) * $script:eR400) / $ShardCount)
    return [long]$script:cells500 * ($u5 - $l5) + `
        [long]$script:cells400 * ($u4 - $l4)
}

function Get-ExpectedShardCells {
    param([int]$Shard, [int]$ShardCount)
    if ($script:repCap -eq 1 -and $ShardCount -gt 1) {
        if ($Shard -gt $script:expectedCells) {return 0}
        return [int]([Math]::Floor(([double]($script:expectedCells - $Shard)) / $ShardCount) + 1)
    }
    return $script:expectedCells
}

$scriptDirectory = Split-Path -Parent $MyInvocation.MyCommand.Path
if (Test-Path -LiteralPath (Join-Path $scriptDirectory 'manifest.json') -PathType Leaf) {
    $runDirectory = $scriptDirectory
}
else {
    $runDirectory = Join-Path (Join-Path $scriptDirectory 'runs') $RunId
}
$manifestPath = Join-Path $runDirectory 'manifest.json'
if (-not (Test-Path -LiteralPath $manifestPath -PathType Leaf)) {
    throw "Manifest not found: $manifestPath"
}
$manifest = Get-Content -LiteralPath $manifestPath -Raw | ConvertFrom-Json
if ([int]$manifest.schemaVersion -ne 1) {throw 'Unsupported or missing manifest schemaVersion.'}
if ([string]$manifest.config.RunId -cne $RunId) {throw 'RunId/manifest mismatch.'}

$shardCount = [int]$manifest.config.NShard
$master = [long]$manifest.config.Master
$repCap = [int]$manifest.config.RepCap
$expectedCells = [int]$manifest.config.ExpectedCells
$script:repCap = $repCap
$script:expectedCells = $expectedCells
$expectedReplications = [long]$manifest.config.ExpectedReplications
$harnessVersion = [string]$manifest.config.HarnessVersion
$codeVersionExpected = [string]$manifest.config.CodeVersionExpected
$rngKind = [string]$manifest.config.RngKind
$stataExecutable = [string]$manifest.stataExecutable
# Study A: 24 cells @ R=500, 92 cells @ R=400 (including N=1600).
$script:cells500 = 24
$script:cells400 = 92
$script:eR500 = if ($repCap -gt 0) {[Math]::Min(500, $repCap)} else {500}
$script:eR400 = if ($repCap -gt 0) {[Math]::Min(400, $repCap)} else {400}

if ($shardCount -lt 1 -or $expectedCells -ne 116 -or
    $repCap -lt 0 -or $repCap -gt 500 -or
    ($repCap -gt 0 -and $repCap -lt $shardCount -and $repCap -ne 1) -or
    ($repCap -eq 1 -and $shardCount -gt $expectedCells) -or
    $expectedReplications -ne ([long]$script:cells500 * $script:eR500 + `
        [long]$script:cells400 * $script:eR400) -or
    $harnessVersion -cne $expectedHarness -or
    $codeVersionExpected -cne $expectedCodeVersion -or
    $rngKind -cne $expectedRng) {
    throw 'Manifest Study A counts are internally inconsistent.'
}
if ($null -eq $manifest.sha256 -or @($manifest.sha256.PSObject.Properties).Count -ne $requiredFiles.Count) {
    throw 'Manifest SHA256 inventory does not match the frozen package inventory.'
}
foreach ($fileName in $requiredFiles) {
    $hashProperty = $manifest.sha256.PSObject.Properties[$fileName]
    if ($null -eq $hashProperty -or [string]$hashProperty.Value -notmatch '^[A-Fa-f0-9]{64}$') {
        throw "Manifest has no valid SHA256 value for: $fileName"
    }
    $path = Join-Path $runDirectory $fileName
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {throw "Frozen file is missing: $fileName"}
    if ((Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash -cne [string]$hashProperty.Value) {
        throw "Frozen file changed after staging: $fileName"
    }
}

function Test-CompletionArtifacts {
    param([int]$Shard, [long]$ExpectedRows, [int]$ExpectedCells)
    $markerPath = Join-Path $runDirectory ('study_a_complete_SH{0}.csv' -f $Shard)
    $resultPath = Join-Path $runDirectory ('study_a_results_SH{0}.csv' -f $Shard)
    $progressPath = Join-Path $runDirectory ('study_a_progress_SH{0}.csv' -f $Shard)
    if (-not (Test-Path -LiteralPath $markerPath -PathType Leaf) -or
        -not (Test-Path -LiteralPath $resultPath -PathType Leaf) -or
        -not (Test-Path -LiteralPath $progressPath -PathType Leaf)) {return $false}
    try {
        $marker = @(Import-Csv -LiteralPath $markerPath)
        if ($marker.Count -ne 1) {return $false}
        $record = $marker[0]
        if ([string]$record.run_id -cne $RunId -or
            [string]$record.harness_version -cne $harnessVersion -or
            [string]$record.code_version_expected -cne $codeVersionExpected -or
            [long]$record.master -ne $master -or
            [int]$record.shard -ne $Shard -or
            [int]$record.nshard -ne $shardCount -or
            [int]$record.cells -ne $ExpectedCells -or
            [long]$record.reps -ne $ExpectedRows -or
            [int]$record.exit_code -ne 0) {return $false}
        $progress = @(Import-Csv -LiteralPath $progressPath)
        if ($progress.Count -ne 1) {return $false}
        $p = $progress[0]
        if ([string]$p.run_id -cne $RunId -or
            [string]$p.harness_version -cne $harnessVersion -or
            [string]$p.code_version_expected -cne $codeVersionExpected -or
            [long]$p.master -ne $master -or
            [int]$p.shard -ne $Shard -or
            [int]$p.nshard -ne $shardCount -or
            [int]$p.cells_completed -ne $ExpectedCells -or
            [long]$p.reps_completed -ne $ExpectedRows) {return $false}
        return (Get-ResultRows -Path $resultPath) -eq $ExpectedRows
    }
    catch {return $false}
}

function Show-StudyAStatus {
    param([bool]$ClearScreen)
    if ($ClearScreen) {Clear-Host}
    $totalRows = 0L
    $completed = 0
    $running = 0
    $failed = 0
    $pending = 0
    $lines = @()

    for ($shard = 1; $shard -le $shardCount; $shard++) {
        $tag = '{0:D3}' -f $shard
        $statusPath = Join-Path $runDirectory "status_SH$tag.json"
        $resultPath = Join-Path $runDirectory ('study_a_results_SH{0}.csv' -f $shard)
        $progressPath = Join-Path $runDirectory ('study_a_progress_SH{0}.csv' -f $shard)
        $expectedRows = Get-ExpectedShardRows -Shard $shard -ShardCount $shardCount
        $expectedShardCells = Get-ExpectedShardCells -Shard $shard -ShardCount $shardCount
        $actualRows = 0L
        $rowText = '0'
        $lastProgress = ''
        $progressBad = $false
        $progressRows = 0L
        if (Test-Path -LiteralPath $progressPath -PathType Leaf) {
            try {
                $progress = $null
                for ($attempt = 1; $attempt -le 5; $attempt++) {
                    try {
                        $progress = @(Import-Csv -LiteralPath $progressPath)
                        break
                    }
                    catch {
                        if ($attempt -eq 5) {throw}
                        Start-Sleep -Milliseconds (25 * $attempt)
                    }
                }
                if ($progress.Count -ne 1) {throw 'progress record count mismatch'}
                $p = $progress[0]
                if ([string]$p.run_id -cne $RunId -or
                    [string]$p.harness_version -cne $harnessVersion -or
                    [string]$p.code_version_expected -cne $codeVersionExpected -or
                    [long]$p.master -ne $master -or
                    [int]$p.shard -ne $shard -or
                    [int]$p.nshard -ne $shardCount -or
                    [int]$p.cells_completed -lt 0 -or
                    [int]$p.cells_completed -gt $expectedShardCells -or
                    [long]$p.reps_completed -lt 0 -or
                    [long]$p.reps_completed -gt $expectedRows) {
                    throw 'progress identity/count mismatch'
                }
                $progressRows = [long]$p.reps_completed
                $lastProgress = ' cell={0}/{1} last={2}:{3}' -f `
                    [int]$p.cells_completed, $expectedShardCells,
                    [string]$p.last_cell_id, [string]$p.last_rep
            }
            catch {$progressBad = $true}
        }
        if (Test-Path -LiteralPath $resultPath -PathType Leaf) {
            try {
                $actualRows = Get-ResultRows -Path $resultPath
                $rowText = [string]$actualRows
                $totalRows += $actualRows
            }
            catch {
                $actualRows = $progressRows
                $rowText = if ($progressBad) {'locked'} else {"$progressRows (progress)"}
                $totalRows += $actualRows
            }
        }
        elseif (-not $progressBad -and $progressRows -gt 0) {
            $actualRows = $progressRows
            $rowText = "$progressRows (progress)"
            $totalRows += $actualRows
        }

        $state = 'pending'
        $exitCode = $null
        if (Test-Path -LiteralPath $statusPath -PathType Leaf) {
            try {
                $status = Get-Content -LiteralPath $statusPath -Raw | ConvertFrom-Json
                $state = [string]$status.state
                $exitCode = $status.exitCode
                if ([string]$status.runId -cne $RunId -or
                    [int]$status.shard -ne $shard -or
                    [long]$status.expectedRows -ne $expectedRows) {
                    $state = 'bad-status'
                }
                elseif ($state -notin @('queued', 'launching', 'running', 'completed', 'failed')) {
                    $state = 'bad-status'
                }
                else {
                    $wrapperAlive = Test-AliveProcess -Status $status -PidField 'wrapperPid' -ExpectedStata $stataExecutable
                    $stataAlive = Test-AliveProcess -Status $status -PidField 'stataPid' -ExpectedStata $stataExecutable
                    if ($state -in @('queued', 'launching', 'running') -and -not ($wrapperAlive -or $stataAlive)) {
                        $statusAgeSeconds = [double]::PositiveInfinity
                        try {
                            $statusAgeSeconds = ([DateTime]::UtcNow -
                                [DateTime]::Parse([string]$status.updatedUtc).ToUniversalTime()).TotalSeconds
                        }
                        catch {}
                        # Queued is a very short handoff state with PID=0.  Its
                        # grace period prevents a monitor poll racing the wrapper claim.
                        if ($state -cne 'queued' -or $statusAgeSeconds -gt 30) {
                            $state = 'stale'
                        }
                    }
                    if ($state -ceq 'completed' -and
                        -not (Test-CompletionArtifacts -Shard $shard `
                            -ExpectedRows $expectedRows -ExpectedCells $expectedShardCells)) {
                        $state = 'bad-completion'
                    }
                }
            }
            catch {$state = 'bad-status'}
        }
        if ($progressBad) {$state = 'bad-progress'}

        switch ($state) {
            'completed' {$completed++}
            'running' {$running++}
            'launching' {$running++}
            'queued' {$running++}
            'failed' {$failed++}
            'stale' {$failed++}
            'bad-status' {$failed++}
            'bad-completion' {$failed++}
            'bad-progress' {$failed++}
            default {$pending++}
        }
        $exitText = if ($null -ne $exitCode) {" exit=$exitCode"} else {''}
        $lines += ('  SH{0} rows={1,-18}/{2,-7} state={3}{4}{5}' -f `
            $tag, $rowText, $expectedRows, $state, $exitText, $lastProgress)
    }

    Write-Host ('=== Study A run {0} @ {1} ===' -f $RunId, (Get-Date -Format 'yyyy-MM-dd HH:mm:ss'))
    Write-Host "Directory: $runDirectory"
    Write-Host ('Config: shards={0}, master={1}, repcap={2}, cells={3}, mode=point/noboot' -f `
        $shardCount, $master, $repCap, $expectedCells)
    Write-Host ('Status: complete={0}, running/queued={1}, failed/stale={2}, pending={3}' -f `
        $completed, $running, $failed, $pending)
    $percent = [Math]::Min(100, [Math]::Round(100.0 * $totalRows / $expectedReplications, 1))
    Write-Host "Replications: $totalRows / $expectedReplications ($percent%)"
    Write-Host 'Shards:'
    $lines | ForEach-Object {Write-Host $_}
    return [pscustomobject]@{
        Complete = $completed; Running = $running; Failed = $failed; Pending = $pending
        Rows = $totalRows
    }
}

if ($Watch -eq 0) {
    $summary = Show-StudyAStatus -ClearScreen $false
    if ($summary.Failed -gt 0) {exit 1}
    exit 0
}
while ($true) {
    $summary = Show-StudyAStatus -ClearScreen $true
    if ($summary.Complete -eq $shardCount -and $summary.Rows -eq $expectedReplications) {
        Write-Host "`nAll $shardCount shards completed with exactly $expectedReplications Study A rows."
        break
    }
    if ($summary.Running -eq 0 -and ($summary.Failed -gt 0 -or $summary.Pending -gt 0)) {
        Write-Error 'No shard is running but the run is incomplete; inspect status/logs, then use -Resume with the identical snapshot/configuration.'
    }
    Start-Sleep -Seconds $Watch
}

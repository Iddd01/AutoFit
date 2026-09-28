# verify_and_merge_study_a.ps1 -- sole supported Study A merge entrypoint.
[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [ValidatePattern('^[A-Za-z0-9][A-Za-z0-9_.-]{0,79}$')]
    [string]$RunId
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version 2.0

$expectedCells = 116
# Two replication groups: 24 stresses at R=500, 92 regular cells at R=400.
$cells500 = 24
$cells400 = 92
$groupRMax = 500
$expectedHarness = 'xtdpthresh_study_a_v0937_r400'
$expectedCodeVersion = '0.9.37'
$expectedRng = 'mt64'
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
$publishedOutputs = @(
    'study_a_all.dta',
    'study_a_all.csv',
    'study_a_summary.dta',
    'study_a_summary.csv',
    'study_a_paired_fd_fod.dta',
    'study_a_paired_fd_fod.csv',
    'study_a_paired_summary.dta',
    'study_a_paired_summary.csv'
)

function Move-AtomicWithRetry {
    param([string]$TemporaryPath, [string]$DestinationPath)
    for ($attempt = 1; $attempt -le 10; $attempt++) {
        try {
            Move-Item -LiteralPath $TemporaryPath -Destination $DestinationPath -Force
            return
        }
        catch {
            if ($attempt -eq 10) {throw}
            Start-Sleep -Milliseconds (25 * $attempt)
        }
    }
}

function Save-JsonAtomic {
    param($Object, [string]$Path)
    $temporaryPath = $Path + '.tmp.' + [Guid]::NewGuid().ToString('N')
    $Object | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $temporaryPath -Encoding UTF8
    Move-AtomicWithRetry -TemporaryPath $temporaryPath -DestinationPath $Path
}

function Save-AsciiLinesAtomic {
    param([string[]]$Lines, [string]$Path)
    $temporaryPath = $Path + '.tmp.' + [Guid]::NewGuid().ToString('N')
    [IO.File]::WriteAllLines($temporaryPath, $Lines, [Text.Encoding]::ASCII)
    Move-AtomicWithRetry -TemporaryPath $temporaryPath -DestinationPath $Path
}

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

function Quote-Argument {
    param([string]$Value)
    return '"' + $Value.Replace('"', '\"') + '"'
}

function Get-FileHashMap {
    param([string]$Directory, [string[]]$Names)
    $hashes = [ordered]@{}
    foreach ($name in $Names) {
        $path = Join-Path $Directory $name
        if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
            throw "Required file is missing: $name"
        }
        $hashes[$name] = (Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash
    }
    return $hashes
}

function Assert-HashMapsEqual {
    param($Before, $After, [string]$Label)
    if ($Before.Count -ne $After.Count) {throw "$Label hash inventory changed."}
    foreach ($name in $Before.Keys) {
        if (-not $After.Contains($name) -or [string]$After[$name] -cne [string]$Before[$name]) {
            throw "$Label changed during verified merge: $name"
        }
    }
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
$manifestCells = [int]$manifest.config.ExpectedCells
$script:repCap = $repCap
$script:expectedCells = $expectedCells
$expectedReplications = [long]$manifest.config.ExpectedReplications
$harnessVersion = [string]$manifest.config.HarnessVersion
$codeVersion = [string]$manifest.config.CodeVersionExpected
$rngKind = [string]$manifest.config.RngKind
$script:cells500 = $cells500
$script:cells400 = $cells400
$script:eR500 = if ($repCap -gt 0) {[Math]::Min(500, $repCap)} else {500}
$script:eR400 = if ($repCap -gt 0) {[Math]::Min(400, $repCap)} else {400}

if ($shardCount -lt 1 -or $shardCount -gt 256 -or
    $master -lt 1 -or $master -gt 2147483000 -or
    $repCap -lt 0 -or $repCap -gt $groupRMax -or
    ($repCap -gt 0 -and $repCap -lt $shardCount -and $repCap -ne 1) -or
    ($repCap -eq 1 -and $shardCount -gt $expectedCells) -or
    $manifestCells -ne $expectedCells -or
    $expectedReplications -ne ([long]$script:cells500 * $script:eR500 + `
        [long]$script:cells400 * $script:eR400) -or
    $harnessVersion -cne $expectedHarness -or
    $codeVersion -cne $expectedCodeVersion -or
    $rngKind -cne $expectedRng) {
    throw 'Manifest Study A configuration is unsupported or internally inconsistent.'
}
if ($null -eq $manifest.sha256 -or @($manifest.sha256.PSObject.Properties).Count -ne $requiredFiles.Count) {
    throw 'Manifest SHA256 inventory does not match the frozen package inventory.'
}

$stataExecutable = [string]$manifest.stataExecutable
if (-not (Test-Path -LiteralPath $stataExecutable -PathType Leaf)) {
    throw "Manifest Stata executable is unavailable: $stataExecutable"
}
$stataExecutable = (Resolve-Path -LiteralPath $stataExecutable).Path

$lockPath = Join-Path $runDirectory '.launcher.lock'
$runLock = $null
$authorizationPath = Join-Path $runDirectory '_merge_authorized.txt'
$mergeMarkerPath = Join-Path $runDirectory 'study_a_merge_complete.csv'
$attestationPath = Join-Path $runDirectory 'study_a_merge_attestation.json'
$mergeSucceeded = $false
$ownsRunLock = $false
$mergeAttemptStarted = $false

try {
    try {
        $runLock = New-Object IO.FileStream(
            $lockPath,
            [IO.FileMode]::OpenOrCreate,
            [IO.FileAccess]::ReadWrite,
            [IO.FileShare]::None
        )
        $ownsRunLock = $true
    }
    catch {
        throw "The Study A run is active or another verifier owns it (lock: $lockPath)."
    }

    $sourceHashes = [ordered]@{}
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
            throw "Frozen file differs from the immutable manifest: $fileName"
        }
        $sourceHashes[$fileName] = $actualHash
    }

    $invokedVerifierHash = (Get-FileHash -LiteralPath $MyInvocation.MyCommand.Path -Algorithm SHA256).Hash
    if ($invokedVerifierHash -cne [string]$sourceHashes['verify_and_merge_study_a.ps1']) {
        throw "Verifier differs from the staged copy. Use: $runDirectory\verify_and_merge_study_a.ps1"
    }

    $statusFiles = @(Get-ChildItem -LiteralPath $runDirectory -Filter 'status_SH*.json' -File)
    $resultFiles = @(Get-ChildItem -LiteralPath $runDirectory -Filter 'study_a_results_SH*.csv' -File)
    $progressFiles = @(Get-ChildItem -LiteralPath $runDirectory -Filter 'study_a_progress_SH*.csv' -File)
    $completionFiles = @(Get-ChildItem -LiteralPath $runDirectory -Filter 'study_a_complete_SH*.csv' -File)
    if ($statusFiles.Count -ne $shardCount -or $resultFiles.Count -ne $shardCount -or
        $progressFiles.Count -ne $shardCount -or $completionFiles.Count -ne $shardCount) {
        throw 'Run artifact counts do not equal the immutable shard count.'
    }

    $rawNames = @('manifest.json')
    for ($shard = 1; $shard -le $shardCount; $shard++) {
        $tag = '{0:D3}' -f $shard
        $statusName = "status_SH$tag.json"
        $resultName = 'study_a_results_SH{0}.csv' -f $shard
        $progressName = 'study_a_progress_SH{0}.csv' -f $shard
        $completionName = 'study_a_complete_SH{0}.csv' -f $shard
        $statusPath = Join-Path $runDirectory $statusName
        $resultPath = Join-Path $runDirectory $resultName
        $progressPath = Join-Path $runDirectory $progressName
        $completionPath = Join-Path $runDirectory $completionName
        foreach ($path in @($statusPath, $resultPath, $progressPath, $completionPath)) {
            if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
                throw "Required shard artifact is missing: $path"
            }
        }

        $expectedRows = Get-ExpectedShardRows -Shard $shard -ShardCount $shardCount
        $expectedShardCells = Get-ExpectedShardCells -Shard $shard -ShardCount $shardCount
        $status = Get-Content -LiteralPath $statusPath -Raw | ConvertFrom-Json
        if ([int]$status.schemaVersion -ne 1 -or
            [string]$status.runId -cne $RunId -or
            [int]$status.shard -ne $shard -or
            [string]$status.state -cne 'completed' -or
            [int]$status.exitCode -ne 0 -or
            [long]$status.expectedRows -ne $expectedRows -or
            [long]$status.resultRows -ne $expectedRows) {
            throw "Shard $shard does not have an exact completed status."
        }

        $completion = @(Import-Csv -LiteralPath $completionPath)
        if ($completion.Count -ne 1) {throw "Shard $shard completion marker is malformed."}
        $c = $completion[0]
        if ([string]$c.run_id -cne $RunId -or
            [string]$c.harness_version -cne $harnessVersion -or
            [string]$c.code_version_expected -cne $codeVersion -or
            [long]$c.master -ne $master -or
            [int]$c.shard -ne $shard -or
            [int]$c.nshard -ne $shardCount -or
            [int]$c.cells -ne $expectedShardCells -or
            [long]$c.reps -ne $expectedRows -or
            [int]$c.exit_code -ne 0) {
            throw "Shard $shard completion marker identity/count mismatch."
        }

        $progress = @(Import-Csv -LiteralPath $progressPath)
        if ($progress.Count -ne 1) {throw "Shard $shard progress marker is malformed."}
        $p = $progress[0]
        if ([string]$p.run_id -cne $RunId -or
            [string]$p.harness_version -cne $harnessVersion -or
            [string]$p.code_version_expected -cne $codeVersion -or
            [long]$p.master -ne $master -or
            [int]$p.shard -ne $shard -or
            [int]$p.nshard -ne $shardCount -or
            [int]$p.cells_completed -ne $expectedShardCells -or
            [long]$p.reps_completed -ne $expectedRows) {
            throw "Shard $shard progress marker identity/count mismatch."
        }
        if ((Get-ResultRows -Path $resultPath) -ne $expectedRows) {
            throw "Shard $shard result row count differs from $expectedRows."
        }
        $rawNames += @($statusName, $resultName, $progressName, $completionName)
    }

    $rawHashesBefore = Get-FileHashMap -Directory $runDirectory -Names $rawNames
    $manifestHash = [string]$rawHashesBefore['manifest.json']
    if ($manifestHash -notmatch '^[A-Fa-f0-9]{64}$') {throw 'Could not establish manifest SHA256.'}

    foreach ($stalePath in @($authorizationPath, $mergeMarkerPath, $attestationPath)) {
        if (Test-Path -LiteralPath $stalePath -PathType Leaf) {
            Remove-Item -LiteralPath $stalePath -Force
        }
    }
    $mergeAttemptStarted = $true
    $nonce = [Guid]::NewGuid().ToString('N')
    Save-AsciiLinesAtomic -Path $authorizationPath -Lines @(
        "run_id=$RunId",
        "manifest_sha256=$manifestHash",
        "nonce=$nonce"
    )

    $mergePath = Join-Path $runDirectory 'merge_study_a.do'
    $stataArguments = '/e do ' + (Quote-Argument $mergePath) + ' ' +
        (Quote-Argument $RunId) + ' ' + (Quote-Argument $nonce)
    $stataProcess = Start-Process -FilePath $stataExecutable -ArgumentList $stataArguments `
        -WorkingDirectory $runDirectory -WindowStyle Hidden -PassThru
    $stataProcess.WaitForExit()
    if ($stataProcess.ExitCode -ne 0) {
        throw "Stata merger exited with code $($stataProcess.ExitCode)."
    }
    if (Test-Path -LiteralPath $authorizationPath -PathType Leaf) {
        throw 'Merger did not consume its single-use authorization.'
    }
    if (-not (Test-Path -LiteralPath $mergeMarkerPath -PathType Leaf)) {
        throw 'Merger did not publish study_a_merge_complete.csv.'
    }
    $expectedMarkerHeader = 'run_id,manifest_sha256,rows,cells,nshard,exit_code'
    $actualMarkerHeader = Get-Content -LiteralPath $mergeMarkerPath -TotalCount 1
    if ([string]$actualMarkerHeader -cne $expectedMarkerHeader) {
        throw 'Merge completion marker has an incompatible schema.'
    }
    $mergeMarker = @(Import-Csv -LiteralPath $mergeMarkerPath)
    if ($mergeMarker.Count -ne 1) {throw 'Merge completion marker must have exactly one row.'}
    $mm = $mergeMarker[0]
    if ([string]$mm.run_id -cne $RunId -or
        [string]$mm.manifest_sha256 -cne $manifestHash -or
        [long]$mm.rows -ne $expectedReplications -or
        [int]$mm.cells -ne $expectedCells -or
        [int]$mm.nshard -ne $shardCount -or
        [int]$mm.exit_code -ne 0) {
        throw 'Merge completion marker identity/count mismatch.'
    }

    $outputHashes = Get-FileHashMap -Directory $runDirectory -Names $publishedOutputs
    $sourceHashesAfter = Get-FileHashMap -Directory $runDirectory -Names $requiredFiles
    Assert-HashMapsEqual -Before $sourceHashes -After $sourceHashesAfter -Label 'Frozen source'
    $rawHashesAfter = Get-FileHashMap -Directory $runDirectory -Names $rawNames
    Assert-HashMapsEqual -Before $rawHashesBefore -After $rawHashesAfter -Label 'Raw input'

    $attestation = [ordered]@{
        schemaVersion = 1
        verifiedUtc = (Get-Date).ToUniversalTime().ToString('o')
        runId = $RunId
        manifestSha256 = $manifestHash
        nonce = $nonce
        config = [ordered]@{
            NShard = $shardCount
            Master = $master
            RepCap = $repCap
            ExpectedCells = $expectedCells
            ExpectedReplications = $expectedReplications
            HarnessVersion = $harnessVersion
            CodeVersionExpected = $codeVersion
            RngKind = $rngKind
        }
        stataExecutable = $stataExecutable
        sourceSha256 = $sourceHashesAfter
        rawSha256 = $rawHashesAfter
        mergeCompletionSha256 = (Get-FileHash -LiteralPath $mergeMarkerPath -Algorithm SHA256).Hash
        outputSha256 = $outputHashes
    }
    Save-JsonAtomic -Object $attestation -Path $attestationPath
    $mergeSucceeded = $true

    Write-Host "Study A verified merge complete: $RunId"
    Write-Host "Rows: $expectedReplications; cells: $expectedCells; shards: $shardCount"
    Write-Host "Attestation: $attestationPath"
}
finally {
    if ($ownsRunLock -and $mergeAttemptStarted -and -not $mergeSucceeded) {
        foreach ($cleanupPath in @($authorizationPath, $mergeMarkerPath, $attestationPath)) {
            try {
                if (Test-Path -LiteralPath $cleanupPath -PathType Leaf) {
                    Remove-Item -LiteralPath $cleanupPath -Force
                }
            }
            catch {}
        }
    }
    if ($null -ne $runLock) {$runLock.Dispose()}
}

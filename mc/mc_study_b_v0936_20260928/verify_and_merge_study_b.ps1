# verify_and_merge_study_b.ps1 -- sole supported Study B merge entrypoint.
[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [ValidatePattern('^[A-Za-z0-9][A-Za-z0-9_.-]{0,79}$')]
    [string]$RunId,
    [switch]$AllowNonFormal
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version 2.0

$expectedCells = 52
$cells500 = 52
$expectedHarness = 'study_b_type4_v0936'
$expectedCodeVersion = '0.9.36'
$expectedRng = 'mt64'
$canonicalMaster = 20260814L
$requiredFiles = @(
    'study_b_shard.do',
    'study_b_worker.do',
    'study_b_cells.do',
    '_registry_probe.do',
    'merge_study_b.do',
    'xtdpthresh.ado',
    'xtdpthresh_p.ado',
    'run_study_b.ps1',
    'monitor_study_b.ps1',
    'verify_and_merge_study_b.ps1',
    'README.md',
    '_DESIGN_B.md'
)
$publishedOutputs = @(
    'study_b_all.dta',
    'study_b_all.csv',
    'study_b_summary.dta',
    'study_b_summary.csv',
    'study_b_paired_fd_fod.dta',
    'study_b_paired_fd_fod.csv',
    'study_b_paired_summary.dta',
    'study_b_paired_summary.csv'
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
    $Object | ConvertTo-Json -Depth 8 |
        Set-Content -LiteralPath $temporaryPath -Encoding UTF8
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
    $u5 = [long][Math]::Floor(([double]$Shard * $script:eR500) / $ShardCount)
    $l5 = [long][Math]::Floor(([double]($Shard - 1) * $script:eR500) / $ShardCount)
    return [long]$script:cells500 * ($u5 - $l5)
}

function Get-ExpectedShardCells {
    param([int]$Shard, [int]$ShardCount)
    if ($script:repCap -eq 1 -and $ShardCount -gt 1) {
        if ($Shard -gt $script:expectedCells) {return 0}
        return [int]([Math]::Floor(
            ([double]($script:expectedCells - $Shard)) / $ShardCount) + 1)
    }
    return $script:expectedCells
}

function Quote-Argument {
    param([string]$Value)
    return '"' + $Value.Replace('"', '\"') + '"'
}

function Invoke-RegistryProbe {
    param([string]$Directory, [string]$StataExecutable)
    $probePath = Join-Path $Directory '_registry_probe.do'
    $markerPath = Join-Path $Directory '_registry_probe.ok'
    Remove-Item -LiteralPath $markerPath -Force -ErrorAction SilentlyContinue
    $arguments = '/e do ' + (Quote-Argument $probePath)
    $process = Start-Process -FilePath $StataExecutable -ArgumentList $arguments `
        -WorkingDirectory $Directory -WindowStyle Hidden -Wait -PassThru
    $expected = 'cells=52;pairs=24;fdonly=4;reps=26000;R500=52;FULL=40;COVERAGE=12;B499=48;B500=4;G199=48;G46=4;GC100=40;GC46=4;GC10=8;RF4=48;RF0=4'
    if ($process.ExitCode -ne 0 -or
        -not (Test-Path -LiteralPath $markerPath -PathType Leaf)) {
        throw "Registry probe failed (Stata exit=$($process.ExitCode))."
    }
    $actual = (Get-Content -LiteralPath $markerPath -Raw).Trim()
    if ($actual -cne $expected) {throw "Registry probe marker mismatch: $actual"}
    return (Get-FileHash -LiteralPath $markerPath -Algorithm SHA256).Hash
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
        if (-not $After.Contains($name) -or
            [string]$After[$name] -cne [string]$Before[$name]) {
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
if ([int]$manifest.schemaVersion -ne 1) {
    throw 'Unsupported or missing manifest schemaVersion.'
}
if ([string]$manifest.config.RunId -cne $RunId) {throw 'RunId/manifest mismatch.'}

$shardCount = [int]$manifest.config.NShard
$master = [long]$manifest.config.Master
$repCap = [int]$manifest.config.RepCap
$manifestCells = [int]$manifest.config.ExpectedCells
$expectedReplications = [long]$manifest.config.ExpectedReplications
$harnessVersion = [string]$manifest.config.HarnessVersion
$codeVersion = [string]$manifest.config.CodeVersionExpected
$rngKind = [string]$manifest.config.RngKind
$bOverride = [string]$manifest.config.B
$gridOverride = [string]$manifest.config.Grid
$gridCiOverride = [string]$manifest.config.GridCI
$memoryGateSeconds = [int]$manifest.config.MemoryGateSeconds
$minHeadroomMB = [int]$manifest.config.MinHeadroomMB
$script:cells500 = $cells500
$script:repCap = $repCap
$script:expectedCells = $expectedCells
$script:eR500 = if ($repCap -gt 0) {[Math]::Min(500, $repCap)} else {500}
$isFormalDesign = ($repCap -eq 0 -and $master -eq $canonicalMaster -and
    $bOverride -ceq '.' -and $gridOverride -ceq '.' -and
    $gridCiOverride -ceq '.')

foreach ($override in @($bOverride, $gridOverride, $gridCiOverride)) {
    if ($override -ne '.' -and $override -notmatch '^[0-9]+$') {
        throw 'Manifest contains an invalid B/grid override.'
    }
}
if (($bOverride -ne '.' -and
        ([long]$bOverride -lt 10 -or [long]$bOverride -gt 100000)) -or
    ($gridOverride -ne '.' -and
        ([long]$gridOverride -lt 10 -or [long]$gridOverride -gt 10000)) -or
    ($gridCiOverride -ne '.' -and
        ([long]$gridCiOverride -lt 10 -or [long]$gridCiOverride -gt 10000))) {
    throw 'Manifest B/grid override is outside the supported range.'
}
if ($shardCount -lt 1 -or $shardCount -gt 256 -or
    $master -lt 1 -or $master -gt 2147483000 -or
    $repCap -lt 0 -or $repCap -gt 500 -or
    ($repCap -gt 0 -and $repCap -lt $shardCount -and $repCap -ne 1) -or
    ($repCap -eq 1 -and $shardCount -gt $expectedCells) -or
    $manifestCells -ne $expectedCells -or
    $expectedReplications -ne
        ([long]$script:cells500 * $script:eR500) -or
    $harnessVersion -cne $expectedHarness -or
    $codeVersion -cne $expectedCodeVersion -or
    $rngKind -cne $expectedRng -or
    $memoryGateSeconds -lt 30 -or $memoryGateSeconds -gt 3600 -or
    $minHeadroomMB -lt 256 -or $minHeadroomMB -gt 32768) {
    throw 'Manifest Study B configuration is unsupported or internally inconsistent.'
}
if (-not $isFormalDesign -and -not $AllowNonFormal) {
    throw 'This run has smoke/sensitivity overrides. Re-run with -AllowNonFormal only if you intentionally want separately attested non-formal outputs.'
}
if ($null -eq $manifest.sha256 -or
    @($manifest.sha256.PSObject.Properties).Count -ne $requiredFiles.Count) {
    throw 'Manifest SHA256 inventory does not match the frozen package inventory.'
}

$stataExecutable = [string]$manifest.stataExecutable
if (-not (Test-Path -LiteralPath $stataExecutable -PathType Leaf)) {
    throw "Manifest Stata executable is unavailable: $stataExecutable"
}
$stataExecutable = (Resolve-Path -LiteralPath $stataExecutable).Path
$stataExecutableHash = (Get-FileHash -LiteralPath $stataExecutable -Algorithm SHA256).Hash
if ([string]$manifest.stataExecutableSha256 -notmatch '^[A-Fa-f0-9]{64}$' -or
    $stataExecutableHash -cne [string]$manifest.stataExecutableSha256) {
    throw 'Stata executable SHA256 differs from the immutable manifest.'
}

$lockPath = Join-Path $runDirectory '.launcher.lock'
$authorizationPath = Join-Path $runDirectory '_merge_authorized.txt'
$mergeMarkerPath = Join-Path $runDirectory 'study_b_merge_complete.csv'
$outputPointerPath = Join-Path $runDirectory 'study_b_outputs_current.json'
$attestationPath = $null
$stageDirectory = $null
$generationDirectory = $null
$runLock = $null
$ownsRunLock = $false
$mergeAttemptStarted = $false
$mergeSucceeded = $false

try {
    try {
        $runLock = New-Object IO.FileStream(
            $lockPath, [IO.FileMode]::OpenOrCreate,
            [IO.FileAccess]::ReadWrite, [IO.FileShare]::None
        )
        $ownsRunLock = $true
    }
    catch {
        throw "The Study B run is active or another verifier owns it (lock: $lockPath)."
    }

    $sourceHashes = [ordered]@{}
    foreach ($fileName in $requiredFiles) {
        $hashProperty = $manifest.sha256.PSObject.Properties[$fileName]
        if ($null -eq $hashProperty -or
            [string]$hashProperty.Value -notmatch '^[A-Fa-f0-9]{64}$') {
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

    $invokedVerifierHash =
        (Get-FileHash -LiteralPath $MyInvocation.MyCommand.Path -Algorithm SHA256).Hash
    if ($invokedVerifierHash -cne
        [string]$sourceHashes['verify_and_merge_study_b.ps1']) {
        throw "Verifier differs from staged copy. Use: $runDirectory\verify_and_merge_study_b.ps1"
    }
    $registryProbeHash = Invoke-RegistryProbe -Directory $runDirectory `
        -StataExecutable $stataExecutable

    $statusFiles = @(Get-ChildItem -LiteralPath $runDirectory -Filter 'status_SH*.json' -File)
    $resultFiles = @(Get-ChildItem -LiteralPath $runDirectory -Filter 'study_b_results_SH*.csv' -File)
    $completionFiles = @(Get-ChildItem -LiteralPath $runDirectory -Filter 'study_b_complete_SH*.csv' -File)
    if ($statusFiles.Count -ne $shardCount -or
        $resultFiles.Count -ne $shardCount -or
        $completionFiles.Count -ne $shardCount) {
        throw 'Run artifact counts do not equal the immutable shard count.'
    }

    $rawNames = @('manifest.json')
    for ($shard = 1; $shard -le $shardCount; $shard++) {
        $tag = '{0:D3}' -f $shard
        $statusName = "status_SH$tag.json"
        $resultName = 'study_b_results_SH{0}.csv' -f $shard
        $completionName = 'study_b_complete_SH{0}.csv' -f $shard
        $statusPath = Join-Path $runDirectory $statusName
        $resultPath = Join-Path $runDirectory $resultName
        $completionPath = Join-Path $runDirectory $completionName
        foreach ($path in @($statusPath, $resultPath, $completionPath)) {
            if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
                throw "Required shard artifact is missing: $path"
            }
        }

        $expectedRows = Get-ExpectedShardRows -Shard $shard -ShardCount $shardCount
        $expectedShardCells = Get-ExpectedShardCells -Shard $shard -ShardCount $shardCount
        $status = Get-Content -LiteralPath $statusPath -Raw | ConvertFrom-Json
        if ([int]$status.schemaVersion -ne 1 -or
            [string]$status.runId -cne $RunId -or
            [string]$status.harnessVersion -cne $harnessVersion -or
            [string]$status.codeVersionExpected -cne $codeVersion -or
            [long]$status.master -ne $master -or
            [int]$status.shard -ne $shard -or
            [int]$status.nshard -ne $shardCount -or
            [string]$status.state -cne 'completed' -or
            [int]$status.exitCode -ne 0 -or
            [long]$status.expectedRows -ne $expectedRows -or
            [int]$status.expectedCells -ne $expectedShardCells -or
            [long]$status.resultRows -ne $expectedRows) {
            throw "Shard $shard does not have an exact completed status."
        }

        $expectedCompletionHeader =
            'run_id,harness_version,code_version_expected,master,shard,nshard,cells,reps,exit_code'
        $completionHeader = Get-Content -LiteralPath $completionPath -TotalCount 1
        if ([string]$completionHeader -cne $expectedCompletionHeader) {
            throw "Shard $shard completion marker has an incompatible schema."
        }
        $completion = @(Import-Csv -LiteralPath $completionPath)
        if ($completion.Count -ne 1) {
            throw "Shard $shard completion marker is malformed."
        }
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
        if ((Get-ResultRows -Path $resultPath) -ne $expectedRows) {
            throw "Shard $shard result row count differs from $expectedRows."
        }
        $rawNames += @($statusName, $resultName, $completionName)
    }

    $rawHashesBefore = Get-FileHashMap -Directory $runDirectory -Names $rawNames
    $manifestHash = [string]$rawHashesBefore['manifest.json']
    foreach ($stalePath in @($authorizationPath, $mergeMarkerPath)) {
        if (Test-Path -LiteralPath $stalePath -PathType Leaf) {
            Remove-Item -LiteralPath $stalePath -Force
        }
    }
    $mergeAttemptStarted = $true
    $nonce = [Guid]::NewGuid().ToString('N')
    $stageName = "_merge_stage_$nonce"
    $generationName = "study_b_outputs_$nonce"
    $stageDirectory = Join-Path $runDirectory $stageName
    $generationDirectory = Join-Path $runDirectory $generationName
    $attestationPath = Join-Path $runDirectory "study_b_merge_attestation_$nonce.json"
    if ((Test-Path -LiteralPath $stageDirectory) -or
        (Test-Path -LiteralPath $generationDirectory) -or
        (Test-Path -LiteralPath $attestationPath)) {
        throw 'Nonce-specific merge artifacts already exist.'
    }
    Save-AsciiLinesAtomic -Path $authorizationPath -Lines @(
        "run_id=$RunId",
        "manifest_sha256=$manifestHash",
        "nonce=$nonce"
    )

    $mergePath = Join-Path $runDirectory 'merge_study_b.do'
    $stataArguments = '/e do ' + (Quote-Argument $mergePath) + ' ' +
        (Quote-Argument $RunId) + ' ' + (Quote-Argument $nonce)
    $stataProcess = Start-Process -FilePath $stataExecutable `
        -ArgumentList $stataArguments -WorkingDirectory $runDirectory `
        -WindowStyle Hidden -PassThru
    $stataProcess.WaitForExit()
    if ($stataProcess.ExitCode -ne 0) {
        throw "Stata merger exited with code $($stataProcess.ExitCode)."
    }
    if (Test-Path -LiteralPath $authorizationPath -PathType Leaf) {
        throw 'Merger did not consume its single-use authorization.'
    }
    if (-not (Test-Path -LiteralPath $mergeMarkerPath -PathType Leaf)) {
        throw 'Merger did not publish study_b_merge_complete.csv.'
    }
    $expectedMarkerHeader = 'run_id,manifest_sha256,output_generation,rows,cells,nshard,exit_code'
    $actualMarkerHeader = Get-Content -LiteralPath $mergeMarkerPath -TotalCount 1
    if ([string]$actualMarkerHeader -cne $expectedMarkerHeader) {
        throw 'Merge completion marker has an incompatible schema.'
    }
    $mergeMarker = @(Import-Csv -LiteralPath $mergeMarkerPath)
    if ($mergeMarker.Count -ne 1) {
        throw 'Merge completion marker must have exactly one row.'
    }
    $mm = $mergeMarker[0]
    if ([string]$mm.run_id -cne $RunId -or
        [string]$mm.manifest_sha256 -cne $manifestHash -or
        [string]$mm.output_generation -cne $stageName -or
        [long]$mm.rows -ne $expectedReplications -or
        [int]$mm.cells -ne $expectedCells -or
        [int]$mm.nshard -ne $shardCount -or
        [int]$mm.exit_code -ne 0) {
        throw 'Merge completion marker identity/count mismatch.'
    }

    if (-not (Test-Path -LiteralPath $stageDirectory -PathType Container)) {
        throw 'Merger did not create its nonce-specific output generation.'
    }
    $stageFiles = @(Get-ChildItem -LiteralPath $stageDirectory -File)
    $unexpectedStageFiles = @($stageFiles | Where-Object {$publishedOutputs -cnotcontains $_.Name})
    if ($stageFiles.Count -ne $publishedOutputs.Count -or $unexpectedStageFiles.Count -ne 0) {
        throw 'Output generation is not exactly the registered eight-file inventory.'
    }
    $stageHashes = Get-FileHashMap -Directory $stageDirectory -Names $publishedOutputs
    $sourceHashesAfter =
        Get-FileHashMap -Directory $runDirectory -Names $requiredFiles
    Assert-HashMapsEqual -Before $sourceHashes -After $sourceHashesAfter -Label 'Frozen source'
    $rawHashesAfter = Get-FileHashMap -Directory $runDirectory -Names $rawNames
    Assert-HashMapsEqual -Before $rawHashesBefore -After $rawHashesAfter -Label 'Raw input'

    # A directory rename on the same volume publishes all eight files as one
    # generation.  The atomic JSON pointer below is the sole current-output
    # selector; an interrupted attempt leaves the previous pointer untouched.
    Move-Item -LiteralPath $stageDirectory -Destination $generationDirectory
    $outputHashes = Get-FileHashMap -Directory $generationDirectory -Names $publishedOutputs
    Assert-HashMapsEqual -Before $stageHashes -After $outputHashes -Label 'Output generation'

    $attestation = [ordered]@{
        schemaVersion = 1
        verifiedUtc = (Get-Date).ToUniversalTime().ToString('o')
        runId = $RunId
        manifestSha256 = $manifestHash
        nonce = $nonce
        config = [ordered]@{
            FormalDesign = $isFormalDesign
            NShard = $shardCount
            Master = $master
            B = $bOverride
            Grid = $gridOverride
            GridCI = $gridCiOverride
            RepCap = $repCap
            ExpectedCells = $expectedCells
            ExpectedReplications = $expectedReplications
            HarnessVersion = $harnessVersion
            CodeVersionExpected = $codeVersion
            RngKind = $rngKind
            MemoryGateSeconds = $memoryGateSeconds
            MinHeadroomMB = $minHeadroomMB
        }
        stataExecutable = $stataExecutable
        stataExecutableSha256 = $stataExecutableHash
        registryProbeSha256 = $registryProbeHash
        outputGeneration = $generationName
        sourceSha256 = $sourceHashesAfter
        rawSha256 = $rawHashesAfter
        mergeCompletionSha256 =
            (Get-FileHash -LiteralPath $mergeMarkerPath -Algorithm SHA256).Hash
        outputSha256 = $outputHashes
    }
    Save-JsonAtomic -Object $attestation -Path $attestationPath
    $attestationHash = (Get-FileHash -LiteralPath $attestationPath -Algorithm SHA256).Hash
    $pointer = [ordered]@{
        schemaVersion = 1
        committedUtc = (Get-Date).ToUniversalTime().ToString('o')
        runId = $RunId
        nonce = $nonce
        outputGeneration = $generationName
        attestation = (Split-Path -Leaf $attestationPath)
        attestationSha256 = $attestationHash
        outputSha256 = $outputHashes
    }
    Save-JsonAtomic -Object $pointer -Path $outputPointerPath
    $mergeSucceeded = $true

    Write-Host "Study B verified merge complete: $RunId"
    Write-Host "Rows: $expectedReplications; cells: $expectedCells; shards: $shardCount"
    Write-Host "Outputs: $generationDirectory"
    Write-Host "Current pointer: $outputPointerPath"
    Write-Host "Attestation: $attestationPath"
}
finally {
    if ($ownsRunLock -and $mergeAttemptStarted -and -not $mergeSucceeded) {
        foreach ($cleanupPath in @($authorizationPath, $mergeMarkerPath, $attestationPath)) {
            try {
                if (-not [string]::IsNullOrWhiteSpace([string]$cleanupPath) -and
                    (Test-Path -LiteralPath $cleanupPath -PathType Leaf)) {
                    Remove-Item -LiteralPath $cleanupPath -Force
                }
            }
            catch {}
        }
        foreach ($cleanupDirectory in @($stageDirectory, $generationDirectory)) {
            try {
                if (-not [string]::IsNullOrWhiteSpace([string]$cleanupDirectory) -and
                    (Split-Path -Parent $cleanupDirectory) -ieq $runDirectory -and
                    (Split-Path -Leaf $cleanupDirectory) -match '^(_merge_stage_|study_b_outputs_)[0-9a-f]{32}$' -and
                    (Test-Path -LiteralPath $cleanupDirectory -PathType Container)) {
                    Remove-Item -LiteralPath $cleanupDirectory -Recurse -Force
                }
            }
            catch {}
        }
    }
    if ($null -ne $runLock) {$runLock.Dispose()}
}

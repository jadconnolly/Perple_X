param(
    [string]$TargetDir = (Join-Path $PSScriptRoot "..\optionfiles"),
    [string]$GitExePath = "C:\Users\jamie\AppData\Local\gitkraken\app-12.0.1\resources\app.asar.unpacked\git\mingw64\bin\git.exe",
    [switch]$Recurse,
    [string]$NamePattern = "*"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

function Get-RelativeGitPath {
    param(
        [Parameter(Mandatory = $true)][string]$BasePath,
        [Parameter(Mandatory = $true)][string]$TargetPath
    )
    $baseUri = New-Object System.Uri(($BasePath.TrimEnd('\') + '\'))
    $targetUri = New-Object System.Uri($TargetPath)
    return [System.Uri]::UnescapeDataString($baseUri.MakeRelativeUri($targetUri).ToString()).Replace('/', '\')
}

# Resolve git executable — try PATH first, then the explicit path.
$gitCommand = $null
if (Get-Command "git" -ErrorAction SilentlyContinue) {
    $gitCommand = "git"
}
elseif (Test-Path -LiteralPath $GitExePath) {
    $gitCommand = $GitExePath
}
else {
    throw "git not found on PATH or at '$GitExePath'. Cannot perform git-tracked renames."
}

$resolvedTargetDir = (Resolve-Path -LiteralPath $TargetDir).Path
$repoRoot = (& $gitCommand -C $resolvedTargetDir rev-parse --show-toplevel 2>&1)
if ($LASTEXITCODE -ne 0 -or -not $repoRoot) {
    throw "Target path is not inside a Git repository: $resolvedTargetDir"
}
$repoRoot = $repoRoot.Trim().Replace('/', '\')

$renamedCount = 0

$filesToRename = Get-ChildItem -LiteralPath $resolvedTargetDir -File -Recurse:$Recurse |
    Where-Object { ($_.Name -like $NamePattern) -and ($_.Name -cmatch '[A-Z]') }

foreach ($file in $filesToRename) {
    $lowerName = $file.Name.ToLowerInvariant()

    if ($file.Name -ceq $lowerName) { continue }

    $lowerFullPath = Join-Path $file.DirectoryName $lowerName
    if ((Test-Path -LiteralPath $lowerFullPath) -and ($file.FullName -ine $lowerFullPath)) {
        Write-Warning "Skipping '$($file.Name)' — lowercase target already exists as a distinct file."
        continue
    }

    # Only rename files tracked by git.
    $oldRel = Get-RelativeGitPath -BasePath $repoRoot -TargetPath $file.FullName
    & $gitCommand -C $repoRoot ls-files --error-unmatch -- "$oldRel" *> $null
    if ($LASTEXITCODE -ne 0) {
        Write-Warning "Skipping '$($file.Name)' — not tracked by git."
        continue
    }

    # Two-step git mv to force case-only rename on Windows.
    $tempName = "__git_case_tmp__{0}{1}" -f ([guid]::NewGuid().ToString("N")), $file.Extension
    $tempFullPath = Join-Path $file.DirectoryName $tempName
    $tmpRel  = Get-RelativeGitPath -BasePath $repoRoot -TargetPath $tempFullPath
    $newRel  = Get-RelativeGitPath -BasePath $repoRoot -TargetPath $lowerFullPath

    & $gitCommand -C $repoRoot mv -- "$oldRel" "$tmpRel"
    if ($LASTEXITCODE -ne 0) { throw "git mv failed: '$oldRel' -> '$tmpRel'" }

    & $gitCommand -C $repoRoot mv -- "$tmpRel" "$newRel"
    if ($LASTEXITCODE -ne 0) { throw "git mv failed: '$tmpRel' -> '$newRel'" }

    $renamedCount++
    Write-Host "Renamed '$($file.Name)' -> '$lowerName'"
}

Write-Host "Done. Renamed $renamedCount file(s)."
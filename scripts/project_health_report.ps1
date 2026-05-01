# project_health_report.ps1

$BuildDir = "build-debug"
$LogFile = "build_output.log"

Write-Host "--- Step 1: Initializing & Building ---" -ForegroundColor Cyan

# 1. Configuration Step (Creates CMakeCache.txt)
if (!(Test-Path "$BuildDir/CMakeCache.txt")) {
    Write-Host "[*] Configuring CMake for the first time..." -ForegroundColor Gray
    cmake -S . -B $BuildDir
}

# 2. Build Step
Write-Host "[*] Compiling project..." -ForegroundColor Gray
cmake --build $BuildDir 2>&1 | Tee-Object -FilePath $LogFile

# Parse the log for actual compiler errors
$Errors = Get-Content $LogFile | Select-String -Pattern "error [A-Z]\d+:", "fatal error", "LNK\d+:"

Write-Host "`n--- Build Summary ---" -ForegroundColor Cyan
if ($Errors) {
    Write-Host "[!] Found $($Errors.Count) compilation errors:" -ForegroundColor Red
    $Errors | ForEach-Object { Write-Host "  - $($_.Line.Trim())" }
}
else {
    Write-Host "[+] No compilation errors found." -ForegroundColor Green
}

Write-Host "`n--- Step 2: Scanning Source for Logical Red Flags ---" -ForegroundColor Cyan
# Using a Regex Boundary \b to avoid matching "toDouble" when looking for "TODO"
$RedFlags = Get-ChildItem -Recurse -Include *.cpp, *.hpp, *.h | 
Select-String -Pattern "\bTODO\b", "\bFIXME\b", "error_here"

if ($RedFlags) {
    Write-Host "[!] Found markers in code:" -ForegroundColor Yellow
    $RedFlags | ForEach-Object { Write-Host "  - $($_.FileName):$($_.LineNumber) -> $($_.Line.Trim())" }
}
else {
    Write-Host "[+] No logical red flags found." -ForegroundColor Green
}
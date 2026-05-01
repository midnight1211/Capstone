$CMakeFile = "CMakeLists.txt"
$Content = Get-Content $CMakeFile -Raw

# 1. Get all .cpp files
$SrcFiles = Get-ChildItem -Path "desktop/src/main/cpp", "src/main/cpp" -Recurse -Filter *.cpp -ErrorAction SilentlyContinue | 
Resolve-Path -Relative |
Where-Object { $_ -notlike "*CoreEngine.cpp" }

if ($SrcFiles.Count -eq 0) {
    Write-Error "No .cpp files found!"
    return
}

# --- THE FIX: Convert backslashes to forward slashes ---
$FormattedFiles = $SrcFiles | ForEach-Object { $_ -replace '\\', '/' }

# 2. Format the list for CMake
$NewSourceList = "`n    # [[SOURCE_START]]`n    " + ($FormattedFiles -join "`n    ") + "`n    # [[SOURCE_END]]"

# 3. Use Regex to replace the block
$Regex = "(?s)# \[\[SOURCE_START\]\].*?# \[\[SOURCE_END\]\]"
$NewContent = $Content -replace $Regex, $NewSourceList

$NewContent | Set-Content $CMakeFile
Write-Host "Success: Added $($FormattedFiles.Count) source files with forward slashes." -ForegroundColor Green
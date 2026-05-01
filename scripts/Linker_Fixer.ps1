# =============================================================================
# MathEngine Linker-Fixer Build Script
# =============================================================================

$CppDir = "C:\Users\marle\CLionProjects\mathengine\desktop\src\main\cpp"
$BuildDir = "$CppDir\build"
$DllPath = "$BuildDir\mathengine.dll"

# 1. Clean and Prep
if (Test-Path $BuildDir) { Remove-Item -Recurse -Force $BuildDir }
New-Item -ItemType Directory -Path $BuildDir | Out-Null

# 2. Gather ALL source files recursively
# This ensures Stat::dispatch and AA::ok are included from their respective folders
$Sources = Get-ChildItem -Path $CppDir -Include *.cpp -Recurse | Where-Object {
    $_.FullName -notlike "*MathBridgeJNI.cpp" # Exclude JNI if handled separately
}

$RequiredFiles = @("AA.cpp", "NA.cpp")
foreach ($req in $RequiredFiles) {
    if ($Sources.Name -notcontains $req) {
        Write-Host "WARNING: Critical file $req NOT found in source list!" -ForegroundColor Red
    }
}

Write-Host "Found $($Sources.Count) source files. Starting compilation..." -ForegroundColor Cyan

# 3. Compiler Flags
$CommonArgs = @(
    "/LD", "/std:c++20", "/utf-8", "/O2", "/EHsc", "/D_USE_MATH_DEFINES", "/FS",
    "/I`"$CppDir`"",
    "/I`"$env:JAVA_HOME\include`"",
    "/I`"$env:JAVA_HOME\include\win32`""
)

# 4. Compile Step-by-Step to prevent filename collisions (Fixes LNK4042)
$ObjFiles = @()

foreach ($file in $Sources) {
    $RelativePath = $file.FullName.Replace($CppDir, "").TrimStart("\").Replace("\", "_")
    $ObjName = $RelativePath.Replace(".cpp", ".obj")
    $ObjPath = Join-Path $BuildDir $ObjName

    Write-Host "Compiling: $($file.Name) -> $ObjName" -ForegroundColor Gray

    # Capture error output specifically
    $errorBuffer = & cl.exe @CommonArgs /c $file.FullName /Fo"$ObjPath" 2>&1

    if ($LASTEXITCODE -ne 0) {
        Write-Host "`n--- COMPILER ERROR IN $($file.Name) ---" -ForegroundColor Red
        $errorBuffer | Out-String | Write-Host -ForegroundColor Yellow
        Write-Error "Failed to compile $($file.Name). Check the output above."
        exit $LASTEXITCODE
    }
    $ObjFiles += $ObjPath
}

# 5. Link Step (Fixes LNK2019/LNK2001)
Write-Host "Linking all objects into mathengine.dll..." -ForegroundColor Green
& link.exe /DLL /OUT:"$DllPath" $ObjFiles

if ($LASTEXITCODE -eq 0) {
    Write-Host "SUCCESS: Build complete at $DllPath" -ForegroundColor Green
}
else {
    Write-Host "LINK ERROR: Check for missing function definitions (unresolved externals)." -ForegroundColor Red
}
# =============================================================================
# build_and_run.ps1
# Math Engine -- JNI Build + Run Script (Windows / PowerShell)
# =============================================================================

[CmdletBinding()]
param(
    [switch] $SkipJavaBuild,
    [switch] $SkipCppBuild,
    [switch] $RunOnly,
    [ValidateSet("msvc", "mingw", "auto")]
    [string] $Compiler = "auto",
    [switch] $VerboseOutput,
    [switch] $StartServer
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

# ---- Output helpers ----------------------------------------------------------

function Write-Header([string]$msg) {
    Write-Host ""
    Write-Host "== $msg ==" -ForegroundColor Cyan
}

function Write-Step([string]$msg) {
    Write-Host "   >> $msg" -ForegroundColor White
}

function Write-Ok([string]$msg) {
    Write-Host "   OK  $msg" -ForegroundColor Green
}

function Write-Warn([string]$msg) {
    Write-Host "   !!  $msg" -ForegroundColor Yellow
}

function Write-Fail([string]$msg) {
    Write-Host "   FAIL  $msg" -ForegroundColor Red
    exit 1
}

function Invoke-Cmd([string]$Cmd, [string[]]$Arguments) {
    if ($VerboseOutput) {
        Write-Host "   $ $Cmd $($Arguments -join ' ')" -ForegroundColor DarkGray
    }
    & $Cmd @Arguments
    if ($LASTEXITCODE -ne 0) {
        Write-Fail "Command failed (exit $LASTEXITCODE): $Cmd"
    }
}

function Test-MSVC { return [bool](Get-Command cl -ErrorAction SilentlyContinue) }
function Test-MinGW { return [bool](Get-Command g++ -ErrorAction SilentlyContinue) }

function Invoke-VcVarsAll {
    $vswhere = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
    if (-not (Test-Path $vswhere)) { Write-Warn "vswhere.exe not found"; return $false }

    $vsPath = & "$vswhere" -latest -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 `
        -property installationPath 2>$null
    if (-not $vsPath) { Write-Warn "No VS install with C++ tools found"; return $false }

    $vcvarsall = Join-Path $vsPath "VC\Auxiliary\Build\vcvarsall.bat"
    if (-not (Test-Path $vcvarsall)) { Write-Warn "vcvarsall.bat not found"; return $false }

    Write-Step "Running vcvarsall.bat x64..."
    $envDump = & cmd.exe /c "`"$vcvarsall`" x64 > nul 2>&1 && set" 2>&1
    foreach ($line in $envDump) {
        if ($line -match "^([^=]+)=(.*)$") {
            [System.Environment]::SetEnvironmentVariable($Matches[1], $Matches[2], "Process")
        }
    }
    return $true
}

# ---- Resolve paths -----------------------------------------------------------

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path

if (Test-Path (Join-Path $ScriptDir "pom.xml")) {
    $ProjectDir = $ScriptDir
}
else {
    $ProjectDir = Split-Path -Parent $ScriptDir
}

if (-not (Test-Path (Join-Path $ProjectDir "pom.xml"))) {
    Write-Fail "Cannot locate pom.xml. Run the script from the project root."
}

$CppDir = Join-Path $ProjectDir "desktop\src\main\cpp"
$BuildDir = Join-Path $CppDir "build"
$IncludeDir = Join-Path $CppDir "include"
$ParserDir = Join-Path $CppDir "parser"
$DllPath = Join-Path $BuildDir "mathengine.dll"

Write-Host ""
Write-Host "  Math Engine -- Build and Run" -ForegroundColor Magenta
Write-Host "  Project root: $ProjectDir"    -ForegroundColor DarkGray

# =============================================================================
# STEP 0 -- Prerequisite checks
# =============================================================================

Write-Header "Checking prerequisites"

# ---- JDK 21 -- use the known install location, fall back to JAVA_HOME / PATH -

$KnownJdk = "C:\Program Files\Java\jdk-21.0.11"

if (Test-Path $KnownJdk) {
    $env:JAVA_HOME = $KnownJdk
    Write-Ok "JDK 21 found at $KnownJdk"
}
elseif ($env:JAVA_HOME) {
    Write-Warn "JDK not found at $KnownJdk -- using JAVA_HOME: $($env:JAVA_HOME)"
}
else {
    $javaExe = Get-Command java -ErrorAction SilentlyContinue
    if ($javaExe) {
        $env:JAVA_HOME = Split-Path (Split-Path $javaExe.Source -Parent) -Parent
        Write-Warn "JDK not found at $KnownJdk and JAVA_HOME not set -- inferred: $($env:JAVA_HOME)"
    }
    else {
        Write-Fail "JDK not found. Expected JDK 21 at: $KnownJdk"
    }
}

# Define the dynamic include paths for the compiler (must come after JAVA_HOME is resolved)
$IncludeFlags = @(
    "/I`"$env:JAVA_HOME\include`"",
    "/I`"$env:JAVA_HOME\include\win32`"",
    "/I`"$IncludeDir`"",
    "/I`"$CppDir`"",
    "/I`"$ParserDir\include_week_four`"",
    "/I`"$CppDir\Linear_Algebra`"",
    "/I`"$CppDir\Statistics`"",
    "/I`"$CppDir\NumberTheory`""
)

$JniHeader = Join-Path $env:JAVA_HOME "include\jni.h"
if (-not (Test-Path $JniHeader)) {
    Write-Fail "jni.h not found at $JniHeader -- verify the JDK installation."
}
Write-Ok "JNI headers found."

if (-not (Get-Command mvn -ErrorAction SilentlyContinue)) {
    Write-Fail "mvn not found on PATH."
}
Write-Ok "Maven found."

# ---- Compiler detection ------------------------------------------------------

$UseMSVC = $false
$UseMingW = $false

if ((-not $SkipCppBuild) -and (-not $RunOnly)) {
    if ($Compiler -eq "msvc") {
        if (Test-MSVC) {
            $UseMSVC = $true
            Write-Ok "Compiler: MSVC (cl.exe)"
        }
        else {
            Write-Step "cl.exe not on PATH -- trying vcvarsall.bat..."
            if (Invoke-VcVarsAll) {
                if (Test-MSVC) { $UseMSVC = $true; Write-Ok "Compiler: MSVC (via vcvarsall)" }
                else { Write-Fail "vcvarsall ran but cl.exe not found" }
            }
            else { Write-Fail "Could not locate MSVC" }
        }
    }
    elseif ($Compiler -eq "mingw") {
        if (Test-MinGW) { $UseMingW = $true; Write-Ok "Compiler: MinGW (g++)" }
        else { Write-Fail "g++ not found on PATH" }
    }
    else {
        # auto
        if (Test-MSVC) {
            $UseMSVC = $true; Write-Ok "Compiler: MSVC (auto-detected)"
        }
        elseif (Test-MinGW) {
            $UseMingW = $true; Write-Ok "Compiler: MinGW (auto-detected)"
        }
        else {
            Write-Step "No compiler on PATH -- trying vcvarsall.bat..."
            if (Invoke-VcVarsAll) {
                if (Test-MSVC) { $UseMSVC = $true; Write-Ok "Compiler: MSVC (via vcvarsall)" }
                else { Write-Fail "vcvarsall ran but cl.exe not found" }
            }
            else { Write-Fail "No C++ compiler found. Install VS C++ workload or MinGW." }
        }
    }
}

# ---- week_4 header check -----------------------------------------------------

$HeadersDir = Join-Path $ParserDir "include_week_four"
$RequiredHeaders = @("evaluator.hpp", "lexer.hpp", "parser.hpp", "value.hpp")
$MissingHeaders = $RequiredHeaders | Where-Object { -not (Test-Path (Join-Path $HeadersDir $_)) }

if ($MissingHeaders) {
    Write-Warn "Missing week_4 headers: $($MissingHeaders -join ', ')"
    Write-Warn "Copy them to: $HeadersDir"
    if ((-not $SkipCppBuild) -and (-not $RunOnly)) {
        Write-Fail "Cannot build DLL without the parser headers."
    }
}
else {
    Write-Ok "All week_4 parser headers present"
}

# =============================================================================
# STEP 1 -- Maven compile
# =============================================================================

if ((-not $SkipJavaBuild) -and (-not $RunOnly)) {
    Write-Header "Compiling Java (mvn compile)"
    Push-Location $ProjectDir
    try { Invoke-Cmd mvn @("compile", "-q") } finally { Pop-Location }
    Write-Ok "Java compiled"
}

# =============================================================================
# STEP 2 -- Compile C++ -> mathengine.dll
# =============================================================================

if ((-not $SkipCppBuild) -and (-not $RunOnly)) {
    Write-Header "Compiling C++ -> mathengine.dll"

    if (-not (Test-Path $BuildDir)) { New-Item -ItemType Directory -Path $BuildDir | Out-Null }

    # Source collection
    $MainSources = @(
        (Join-Path $CppDir "MathBridgeJNI.cpp"),
        (Join-Path $CppDir "CoreEngine.cpp")
    )

    # Exclude stats.cpp and duplicate files if they cause LNK2005 errors
    $TopicSources = Get-ChildItem -Path $CppDir -Include *.cpp -Recurse |
    Where-Object {
        $_.FullName -notlike "*MathBridgeJNI.cpp" -and
        $_.FullName -notlike "*CoreEngine.cpp" -and
        $_.Name -notmatch "stats.cpp" # Removing monolithic stats file to avoid LNK2005
    } |
    Select-Object -ExpandProperty FullName

    $Sources = $MainSources + $TopicSources
    Write-Step "Source files: $($Sources.Count)"

    if ($UseMSVC) {
        Write-Step "Building with MSVC (cl.exe)"
        $ClArgs = @(
            "/LD",
            "/std:c++20",     # Updated for C++20 features
            "/utf-8",         # Added to fix C4566 (Unicode/Alpha characters)
            "/O2",
            "/EHsc",
            "/D_USE_MATH_DEFINES",
            "/Fe:`"$DllPath`"",
            "/Fo:`"$BuildDir\\`""
        )
        # Combine the fixed ClArgs, our new IncludeFlags, and the Sources
        Invoke-Cmd cl ($ClArgs + $IncludeFlags + $Sources)
    }

    if ($UseMingW) {
        Write-Step "Building with MinGW (g++)"
        $GppArgs = @(
            "-shared", "-std=c++17", "-O2",
            "-D_USE_MATH_DEFINES",
            "-I`"$env:JAVA_HOME\include`"",
            "-I`"$env:JAVA_HOME\include\win32`"",
            "-I`"$IncludeDir`"",
            "-I`"$CppDir`"",
            "-I`"$ParserDir\include_week_four`"",
            "-o", "`"$DllPath`""
        ) + $Sources
        Invoke-Cmd g++ $GppArgs
    }

    if (Test-Path $DllPath) {
        $size = [Math]::Round((Get-Item $DllPath).Length / 1KB, 1)
        Write-Ok "DLL built: $DllPath ($size KB)"
    }
    else {
        Write-Fail "DLL not found after build"
    }
}

# =============================================================================
# STEP 3 -- Launch
# =============================================================================

if ($StartServer) {
    Write-Header "Starting Math Engine Server"
    $ServerDir = Join-Path $ProjectDir "server"
    if (-not (Test-Path (Join-Path $ServerDir "pom.xml"))) {
        Write-Fail "server\pom.xml not found"
    }
    Push-Location $ServerDir
    try {
        Invoke-Cmd mvn @("spring-boot:run", "-Djava.library.path=`"$BuildDir`"")
    }
    finally { Pop-Location }
}
else {
    Write-Header "Launching Math Engine Desktop"
    Push-Location $ProjectDir
    try {
        Invoke-Cmd mvn @("javafx:run", "-Djava.library.path=`"$BuildDir`"")
    }
    finally { Pop-Location }
}

Write-Ok "Done."
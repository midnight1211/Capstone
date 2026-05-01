# build.ps1 -- project root shortcut
# Forwards all arguments to the real build script.
# Usage:  .\build.ps1  or  .\build.ps1 -StartServer  or  .\build.ps1 -SkipCppBuild
& "$PSScriptRoot\scripts\build_and_run.ps1" @args

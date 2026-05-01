# Math Engine

A cross-device math computation application built entirely in **C++ and Java**, connected via JNI.

- **Desktop** — JavaFX application with rendered LaTeX input/output and full ADA accessibility
- **Mobile/Tablet** — browser-based UI served over your local WiFi network
- **Engine** — your week_4 C++ Lexer/Parser/Evaluator, called from Java via a shared native library
- **Server** — Spring Boot REST API with SQLite user accounts, JWT auth, and history sync

---

## Table of Contents

- [Math Engine](#math-engine)
  - [Table of Contents](#table-of-contents)
  - [Architecture](#architecture)
  - [Project Structure](#project-structure)
  - [Prerequisites](#prerequisites)
  - [Quick Start](#quick-start)
    - [Step 1 -- Copy your week\_4 parser headers](#step-1----copy-your-week_4-parser-headers)
    - [Step 2 -- Build the C++ engine and launch the desktop app](#step-2----build-the-c-engine-and-launch-the-desktop-app)
    - [Step 3 -- Start the server for phone and tablet access (optional)](#step-3----start-the-server-for-phone-and-tablet-access-optional)
  - [Build Script Reference](#build-script-reference)
  - [REST API Reference](#rest-api-reference)
    - [Auth](#auth)
    - [Compute](#compute)
    - [History](#history)
    - [Preferences](#preferences)
    - [Export](#export)
  - [Supported Expressions](#supported-expressions)
    - [Precision modes](#precision-modes)
  - [Accounts and Sync](#accounts-and-sync)
  - [Accessibility Features](#accessibility-features)
  - [Configuration](#configuration)
  - [Troubleshooting](#troubleshooting)
    - [Phone cannot reach the server](#phone-cannot-reach-the-server)
    - [JNI handshake test](#jni-handshake-test)

---

## Architecture

```txt
Phone / Tablet (browser)
    |
    |  HTTP on local WiFi
    v
Spring Boot Server  (port 8080)
    |  calls via JNI
    v
mathengine.dll  <-- C++ Lexer -> Parser -> Evaluator (week_4)
    ^
    |  calls via JNI (or HTTP if server is running)
JavaFX Desktop App
```

The C++ engine is compiled once into a single `mathengine.dll`. The Spring Boot server loads
it via `System.loadLibrary` and exposes it over HTTP. Every device — desktop, phone, tablet —
reaches the same engine through the server's REST API. No C++ runs on any client device.

**Guest users** can compute freely without an account. Signing in unlocks history sync,
preference sync (theme + font size), and result export across all your devices.

---

## Project Structure

```txt
mathengine/
|
|
+---desktop
|   |   \---src
|   |       \---main
|   |           \---cpp
|   |               \---Calculus
|   |                   |   calculus.vcxproj
|   |                   |   calculus.vcxproj.filters
|   |                   |   cmake_install.cmake
|   |                   |   
|   |                   +---calculus.dir
|   |                   |   +---Debug
|   |                   |   +---MinSizeRel
|   |                   |   +---Release
|   |                   |   \---RelWithDebInfo
|   |                   \---CMakeFiles
|   |                           generate.stamp
|   |                           generate.stamp.depend
|   |                           
|   +---mathengine.dir
|   |   +---Debug
|   |   +---MinSizeRel
|   |   +---Release
|   |   \---RelWithDebInfo
|   \---ZERO_CHECK.dir
|       +---Debug
|       +---MinSizeRel
|       +---Release
|       \---RelWithDebInfo
+---build-debug
|   |   ALL_BUILD.vcxproj
|   |   ALL_BUILD.vcxproj.filters
|   |   CMakeCache.txt
|   |   cmake_install.cmake
|   |   mathengine.slnx
|   |   mathengine.vcxproj
|   |   mathengine.vcxproj.filters
|   |   ZERO_CHECK.vcxproj
|   |   ZERO_CHECK.vcxproj.filters
|   |   
|   +---ALL_BUILD.dir
|   |   +---Debug
|   |   +---MinSizeRel
|   |   +---Release
|   |   \---RelWithDebInfo
|   +---CMakeFiles
|   |   |   cmake.check_cache
|   |   |   CMakeConfigureLog.yaml
|   |   |   generate.stamp
|   |   |   generate.stamp.depend
|   |   |   generate.stamp.list
|   |   |   InstallScripts.json
|   |   |   TargetDirectories.txt
|   |   |   
|   |   +---1682133aee2f7af36bbc2ab94f3edfe7
|   |   |       generate.stamp.rule
|   |   |       
|   |   +---4.3.1
|   |   |   |   CMakeCXXCompiler.cmake
|   |   |   |   CMakeDetermineCompilerABI_CXX.bin
|   |   |   |   CMakeRCCompiler.cmake
|   |   |   |   CMakeSystem.cmake
|   |   |   |   VCTargetsPath.txt
|   |   |   |   VCTargetsPath.vcxproj
|   |   |   |   
|   |   |   +---CompilerIdCXX
|   |   |   |   |   CMakeCXXCompilerId.cpp
|   |   |   |   |   CompilerIdCXX.exe
|   |   |   |   |   CompilerIdCXX.vcxproj
|   |   |   |   |   
|   |   |   |   +---Debug
|   |   |   |   |   |   CMakeCXXCompilerId.obj
|   |   |   |   |   |   CompilerIdCXX.exe.recipe
|   |   |   |   |   |   
|   |   |   |   |   \---CompilerIdCXX.tlog
|   |   |   |   |           CL.command.1.tlog
|   |   |   |   |           Cl.items.tlog
|   |   |   |   |           CL.read.1.tlog
|   |   |   |   |           CL.write.1.tlog
|   |   |   |   |           CompilerIdCXX.lastbuildstate
|   |   |   |   |           link.command.1.tlog
|   |   |   |   |           link.read.1.tlog
|   |   |   |   |           link.secondary.1.tlog
|   |   |   |   |           link.write.1.tlog
|   |   |   |   |           
|   |   |   |   \---tmp
|   |   |   +---VCTargetsPath
|   |   |   |   \---x64
|   |   |   |       \---Debug
|   |   |   |           |   VCTargetsPath.recipe
|   |   |   |           |   
|   |   |   |           \---VCTargetsPath.tlog
|   |   |   |                   VCTargetsPath.lastbuildstate
|   |   |   |                   
|   |   |   \---x64
|   |   |       \---Debug
|   |   \---pkgRedirects
|   +---Debug
|   |       mathengine.exp
|   |       mathengine.lib
|   |       
|   +---desktop
|   |   \---src
|   |       \---main
|   |           \---cpp
|   |               \---Calculus
|   |                   |   calculus.vcxproj
|   |                   |   calculus.vcxproj.filters
|   |                   |   cmake_install.cmake
|   |                   |   
|   |                   +---calculus.dir
|   |                   |   +---Debug
|   |                   |   |   |   calculus.lib.recipe
|   |                   |   |   |   Calculus.obj
|   |                   |   |   |   Convergence.obj
|   |                   |   |   |   Derivative.obj
|   |                   |   |   |   Expression.obj
|   |                   |   |   |   Implicit.obj
|   |                   |   |   |   Limits.obj
|   |                   |   |   |   LineIntegral.obj
|   |                   |   |   |   Multivariable.obj
|   |                   |   |   |   Numerical.obj
|   |                   |   |   |   Optimize.obj
|   |                   |   |   |   Partial.obj
|   |                   |   |   |   Series.obj
|   |                   |   |   |   Simplify.obj
|   |                   |   |   |   SurfaceIntegral.obj
|   |                   |   |   |   Symbolic.obj
|   |                   |   |   |   Theorems.obj
|   |                   |   |   |   VectorOps.obj
|   |                   |   |   |   Week4Bridge.obj
|   |                   |   |   |   
|   |                   |   |   \---calculus.tlog
|   |                   |   |           calculus.lastbuildstate
|   |                   |   |           CL.command.1.tlog
|   |                   |   |           Cl.items.tlog
|   |                   |   |           CL.read.1.tlog
|   |                   |   |           CL.write.1.tlog
|   |                   |   |           CustomBuild.command.1.tlog
|   |                   |   |           CustomBuild.read.1.tlog
|   |                   |   |           CustomBuild.write.1.tlog
|   |                   |   |           Lib-link.read.1.tlog
|   |                   |   |           Lib-link.write.1.tlog
|   |                   |   |           Lib.command.1.tlog
|   |                   |   |           
|   |                   |   +---MinSizeRel
|   |                   |   +---Release
|   |                   |   \---RelWithDebInfo
|   |                   +---CMakeFiles
|   |                   |       generate.stamp
|   |                   |       generate.stamp.depend
|   |                   |       
|   |                   \---Debug
|   |                           calculus.lib
|   |                           calculus.pdb
|   |                           
|   +---mathengine.dir
|   |   +---Debug
|   |   |   |   AA.obj
|   |   |   |   Additional.obj
|   |   |   |   BasicOps.obj
|   |   |   |   CA.obj
|   |   |   |   CAA.obj
|   |   |   |   CAF.obj
|   |   |   |   Calculus.obj
|   |   |   |   CAR.obj
|   |   |   |   CAS.obj
|   |   |   |   CASing.obj
|   |   |   |   CASpecial.obj
|   |   |   |   Combinatorics.obj
|   |   |   |   CommonUtils.obj
|   |   |   |   Congruences.obj
|   |   |   |   Convergence.obj
|   |   |   |   CoreEngine.obj
|   |   |   |   Crypto.obj
|   |   |   |   DE.obj
|   |   |   |   DEAdv.obj
|   |   |   |   Derivative.obj
|   |   |   |   Determinant.obj
|   |   |   |   Diophantine.obj
|   |   |   |   Dispatch.obj
|   |   |   |   Divisibility.obj
|   |   |   |   EA.obj
|   |   |   |   Eigen.obj
|   |   |   |   EV.obj
|   |   |   |   evaluator.obj
|   |   |   |   Expression.obj
|   |   |   |   Fields.obj
|   |   |   |   GenFunc.obj
|   |   |   |   Groups.obj
|   |   |   |   GT.obj
|   |   |   |   Implicit.obj
|   |   |   |   Interp.obj
|   |   |   |   Inverse.obj
|   |   |   |   LA.obj
|   |   |   |   LASpecial.obj
|   |   |   |   lexer.obj
|   |   |   |   Limits.obj
|   |   |   |   LineIntegral.obj
|   |   |   |   Logic.obj
|   |   |   |   LS.obj
|   |   |   |   LUDecomp.obj
|   |   |   |   MathBridgeJNI.obj
|   |   |   |   Matrix.obj
|   |   |   |   Mod.obj
|   |   |   |   Multivariable.obj
|   |   |   |   NA.obj
|   |   |   |   ND.obj
|   |   |   |   NI.obj
|   |   |   |   Norms.obj
|   |   |   |   NTSpecial.obj
|   |   |   |   NumberTheory.obj
|   |   |   |   Numerical.obj
|   |   |   |   Optimize.obj
|   |   |   |   parser.obj
|   |   |   |   Partial.obj
|   |   |   |   Perms.obj
|   |   |   |   Primes.obj
|   |   |   |   PRings.obj
|   |   |   |   PT.obj
|   |   |   |   P_and_F.obj
|   |   |   |   QR.obj
|   |   |   |   RecRel.obj
|   |   |   |   RF.obj
|   |   |   |   Rings.obj
|   |   |   |   RowReduction.obj
|   |   |   |   Series.obj
|   |   |   |   Simplify.obj
|   |   |   |   Solving.obj
|   |   |   |   Strassen.obj
|   |   |   |   SurfaceIntegral.obj
|   |   |   |   SVD.obj
|   |   |   |   Symbolic.obj
|   |   |   |   Theorems.obj
|   |   |   |   vc145.pdb
|   |   |   |   Week4Bridge.obj
|   |   |   |   
|   |   |   +---desktop
|   |   |   |   \---src
|   |   |   |       \---main
|   |   |   |           \---cpp
|   |   |   |               +---AppliedMath
|   |   |   |               |   \---DM
|   |   |   |               +---Calculus
|   |   |   |               |   \---vectorcalc
|   |   |   |               +---DIscreteMath
|   |   |   |               \---Linear_Algebra
|   |   |   \---mathengine.tlog
|   |   |           CL.command.1.tlog
|   |   |           CL.read.1.tlog
|   |   |           CL.write.1.tlog
|   |   |           CustomBuild.command.1.tlog
|   |   |           CustomBuild.read.1.tlog
|   |   |           CustomBuild.write.1.tlog
|   |   |           link-cvtres.read.1.tlog
|   |   |           link-cvtres.write.1.tlog
|   |   |           link-rc.read.1.tlog
|   |   |           link-rc.write.1.tlog
|   |   |           link.command.1.tlog
|   |   |           link.read.1.tlog
|   |   |           link.read.2.tlog
|   |   |           link.read.3.tlog
|   |   |           link.read.4.tlog
|   |   |           link.write.1.tlog
|   |   |           mathengine.lastbuildstate
|   |   |           unsuccessfulbuild
|   |   |           
|   |   +---MinSizeRel
|   |   +---Release
|   |   \---RelWithDebInfo
|   +---x64
|   |   \---Debug
|   |       \---ZERO_CHECK
|   |           |   ZERO_CHECK.recipe
|   |           |   
|   |           \---ZERO_CHECK.tlog
|   |                   CustomBuild.command.1.tlog
|   |                   CustomBuild.read.1.tlog
|   |                   CustomBuild.write.1.tlog
|   |                   ZERO_CHECK.lastbuildstate
|   |                   
|   \---ZERO_CHECK.dir
|       +---Debug
|       +---MinSizeRel
|       +---Release
|       \---RelWithDebInfo
+---cmake-build-debug
|   |   .gitignore
|   |   .ninja_deps
|   |   .ninja_log
|   |   build.ninja
|   |   capstone_stuff.exe
|   |   CMakeCache.txt
|   |   cmake_install.cmake
|   |   compile_commands.json
|   |   
|   +---.cmake
|   |   \---api
|   |       \---v1
|   |           +---query
|   |           |       cache-v2
|   |           |       cmakeFiles-v1
|   |           |       codemodel-v2
|   |           |       toolchains-v1
|   |           |       
|   |           \---reply
|   |                   error-2026-04-15T23-16-00-0425.json
|   |                   error-2026-04-15T23-18-10-0524.json
|   |                   error-2026-04-15T23-30-02-0572.json
|   |                   error-2026-04-15T23-30-23-0549.json
|   |                   error-2026-04-15T23-30-48-0588.json
|   |                   
|   +---CMakeFiles
|   |   |   clion-Debug-log.txt
|   |   |   clion-environment.txt
|   |   |   cmake.check_cache
|   |   |   CMakeConfigureLog.yaml
|   |   |   
|   |   +---4.2.2
|   |   |   |   CMakeCXXCompiler.cmake
|   |   |   |   CMakeDetermineCompilerABI_CXX.bin
|   |   |   |   CMakeRCCompiler.cmake
|   |   |   |   CMakeSystem.cmake
|   |   |   |   
|   |   |   \---CompilerIdCXX
|   |   |       |   a.exe
|   |   |       |   CMakeCXXCompilerId.cpp
|   |   |       |   
|   |   |       \---tmp
|   |   \---pkgRedirects
|   +---desktop
|   |   \---src
|   |       \---main
|   |           \---cpp
|   |               \---calculus
|   |                   |   cmake_install.cmake
|   |                   |   
|   |                   \---CMakeFiles
|   |                       \---calculus.dir
|   |                           +---core
|   |                           +---differentiation
|   |                           +---integration
|   |                           +---limits
|   |                           +---optimization
|   |                           +---series
|   |                           \---vectorcalc
|   \---Testing
|       \---Temporary
|               LastTest.log
|               
+---desktop
|   |   capstone-workspace.code-workspace
|   |   
|   \---src
|       +---main
|       |   +---cpp
|       |   |   |   CommonUtils.cpp
|       |   |   |   CommonUtils.hpp
|       |   |   |   CoreEngine.cpp
|       |   |   |   CoreEngine.hpp
|       |   |   |   MathBridgeJNI.cpp
|       |   |   |   MathCore.hpp
|       |   |   |   Preprocessor.cpp
|       |   |   |   
|       |   |   +---AbstractAlgebra
|       |   |   |       AA.cpp
|       |   |   |       AA.hpp
|       |   |   |       
|       |   |   +---AppliedMath
|       |   |   |       AM.cpp
|       |   |   |       AM.hpp
|       |   |   |       
|       |   |   +---build
|       |   |   |       AA.obj
|       |   |   |       AbstractAlgebra_AA.obj
|       |   |   |       Advanced_Statistics.obj
|       |   |   |       AM.obj
|       |   |   |       AppliedMath_AM.obj
|       |   |   |       CA.obj
|       |   |   |       Calculus.obj
|       |   |   |       Calculus_Calculus.obj
|       |   |   |       Calculus_core_Expression.obj
|       |   |   |       Calculus_core_Simplify.obj
|       |   |   |       Calculus_core_Week4Bridge.obj
|       |   |   |       Calculus_differentiation_Derivative.obj
|       |   |   |       Calculus_differentiation_Implicit.obj
|       |   |   |       Calculus_differentiation_Partial.obj
|       |   |   |       Calculus_integration_Multivariable.obj
|       |   |   |       Calculus_integration_Numerical.obj
|       |   |   |       Calculus_integration_Symbolic.obj
|       |   |   |       Calculus_limits_Limits.obj
|       |   |   |       Calculus_optimization_Optimize.obj
|       |   |   |       Calculus_series_Convergence.obj
|       |   |   |       Calculus_series_Series.obj
|       |   |   |       Calculus_vectorcalc_LineIntegral.obj
|       |   |   |       Calculus_vectorcalc_SurfaceIntegral.obj
|       |   |   |       Calculus_vectorcalc_Theorems.obj
|       |   |   |       Calculus_vectorcalc_VectorOps.obj
|       |   |   |       CommonUtils.obj
|       |   |   |       ComplexAnalysis_CA.obj
|       |   |   |       Convergence.obj
|       |   |   |       CoreEngine.obj
|       |   |   |       DE.obj
|       |   |   |       DEAdv.obj
|       |   |   |       Derivative.obj
|       |   |   |       DiffEq_DE.obj
|       |   |   |       DiffEq_DEAdv.obj
|       |   |   |       DIscreteMath_DM.obj
|       |   |   |       DM.obj
|       |   |   |       evaluator.obj
|       |   |   |       Expression.obj
|       |   |   |       Geom.obj
|       |   |   |       Geometry_Geom.obj
|       |   |   |       Implicit.obj
|       |   |   |       LA.obj
|       |   |   |       lexer.obj
|       |   |   |       Limits.obj
|       |   |   |       Linear_Algebra_LA.obj
|       |   |   |       LineIntegral.obj
|       |   |   |       MathBridgeJNI.obj
|       |   |   |       mathengine.dll
|       |   |   |       mathengine.exp
|       |   |   |       mathengine.lib
|       |   |   |       Multivariable.obj
|       |   |   |       NA.obj
|       |   |   |       NumberTheory.obj
|       |   |   |       NumberTheory_NumberTheory.obj
|       |   |   |       Numerical.obj
|       |   |   |       NumericalAnalysis_NA.obj
|       |   |   |       Optimize.obj
|       |   |   |       parser.obj
|       |   |   |       parser_src_week_four_evaluator.obj
|       |   |   |       parser_src_week_four_lexer.obj
|       |   |   |       parser_src_week_four_parser.obj
|       |   |   |       Partial.obj
|       |   |   |       Preprocessor.obj
|       |   |   |       ProbabilityTheory_PT.obj
|       |   |   |       PT.obj
|       |   |   |       Series.obj
|       |   |   |       Simplify.obj
|       |   |   |       Statistics.obj
|       |   |   |       Statistics_Advanced_Statistics.obj
|       |   |   |       Statistics_Statistics.obj
|       |   |   |       SurfaceIntegral.obj
|       |   |   |       Symbolic.obj
|       |   |   |       Theorems.obj
|       |   |   |       VectorOps.obj
|       |   |   |       Week4Bridge.obj
|       |   |   |       
|       |   |   +---Calculus
|       |   |   |   |   Calculus.cpp
|       |   |   |   |   Calculus.hpp
|       |   |   |   |   CMakeLists.txt
|       |   |   |   |   
|       |   |   |   +---core
|       |   |   |   |       Expression.cpp
|       |   |   |   |       Expression.hpp
|       |   |   |   |       Simplify.cpp
|       |   |   |   |       Simplify.hpp
|       |   |   |   |       Week4Bridge.cpp
|       |   |   |   |       Week4Bridge.hpp
|       |   |   |   |       
|       |   |   |   +---differentiation
|       |   |   |   |       Derivative.cpp
|       |   |   |   |       Derivative.hpp
|       |   |   |   |       Implicit.cpp
|       |   |   |   |       Implicit.hpp
|       |   |   |   |       Partial.cpp
|       |   |   |   |       Partial.hpp
|       |   |   |   |       
|       |   |   |   +---integration
|       |   |   |   |       Multivariable.cpp
|       |   |   |   |       Multivariable.hpp
|       |   |   |   |       Numerical.cpp
|       |   |   |   |       Numerical.hpp
|       |   |   |   |       Symbolic.cpp
|       |   |   |   |       Symbolic.hpp
|       |   |   |   |       
|       |   |   |   +---limits
|       |   |   |   |       Limits.cpp
|       |   |   |   |       Limits.hpp
|       |   |   |   |       
|       |   |   |   +---optimization
|       |   |   |   |       Optimize.cpp
|       |   |   |   |       Optimize.hpp
|       |   |   |   |       
|       |   |   |   +---series
|       |   |   |   |       Convergence.cpp
|       |   |   |   |       Convergence.hpp
|       |   |   |   |       Series.cpp
|       |   |   |   |       Series.hpp
|       |   |   |   |       
|       |   |   |   \---vectorcalc
|       |   |   |           LineIntegral.cpp
|       |   |   |           LineIntegral.hpp
|       |   |   |           SurfaceIntegral.cpp
|       |   |   |           SurfaceIntegral.hpp
|       |   |   |           Theorems.cpp
|       |   |   |           Theorems.hpp
|       |   |   |           VectorOps.cpp
|       |   |   |           VectorOps.hpp
|       |   |   |           
|       |   |   +---ComplexAnalysis
|       |   |   |       CA.cpp
|       |   |   |       CA.hpp
|       |   |   |       
|       |   |   +---DiffEq
|       |   |   |       DE.cpp
|       |   |   |       DE.hpp
|       |   |   |       DEAdv.cpp
|       |   |   |       
|       |   |   +---DIscreteMath
|       |   |   |       DM.cpp
|       |   |   |       DM.hpp
|       |   |   |       
|       |   |   +---Geometry
|       |   |   |       Geom.cpp
|       |   |   |       Geom.hpp
|       |   |   |       
|       |   |   +---include
|       |   |   |       com_mathengine_jni_MathBridge.h
|       |   |   |       
|       |   |   +---Linear_Algebra
|       |   |   |       LA.cpp
|       |   |   |       LA.hpp
|       |   |   |       
|       |   |   +---NumberTheory
|       |   |   |       NumberTheory.cpp
|       |   |   |       NumberTheory.hpp
|       |   |   |       
|       |   |   +---NumericalAnalysis
|       |   |   |       NA.cpp
|       |   |   |       NA.hpp
|       |   |   |       
|       |   |   +---parser
|       |   |   |   +---include_week_four
|       |   |   |   |       evaluator.hpp
|       |   |   |   |       lexer.hpp
|       |   |   |   |       parser.hpp
|       |   |   |   |       value.hpp
|       |   |   |   |       
|       |   |   |   \---src_week_four
|       |   |   |           evaluator.cpp
|       |   |   |           lexer.cpp
|       |   |   |           parser.cpp
|       |   |   |           
|       |   |   +---ProbabilityTheory
|       |   |   |       PT.cpp
|       |   |   |       PT.hpp
|       |   |   |       
|       |   |   \---Statistics
|       |   |           Advanced_Statistics.cpp
|       |   |           Statistics.cpp
|       |   |           Statistics.hpp
|       |   |           Statistics_internal.hpp
|       |   |           
|       |   +---java
|       |   |   \---com
|       |   |       \---mathengine
|       |   |           +---jni
|       |   |           |       MathBridge.java
|       |   |           |       
|       |   |           +---model
|       |   |           |       PrecisionMode.java
|       |   |           |       
|       |   |           \---ui
|       |   |                   AccessibilityToolbar.java
|       |   |                   AuthService.java
|       |   |                   CalcExpressionBuilder.java
|       |   |                   DatasetImportPanel.java
|       |   |                   GraphPanel.java
|       |   |                   InputPanel.java
|       |   |                   LatexRenderer.java
|       |   |                   Launcher.java
|       |   |                   LoginScreen.java
|       |   |                   MainLayout.java
|       |   |                   MathEngineApp.java
|       |   |                   MathKeyboard.java
|       |   |                   OperationPanel.java
|       |   |                   OutputPanel.java
|       |   |                   QuickFunctionsPanel.java
|       |   |                   regUser.java
|       |   |                   SettingsPanel.java
|       |   |                   ThemeManager.java
|       |   |                   
|       |   \---resources
|       |       \---css
|       |               base.css
|       |               dark.css
|       |               font-large.css
|       |               font-normal.css
|       |               font-x-large.css
|       |               high-contrast.css
|       |               light.css
|       |               
|       \---target
|           +---classes
|           \---test-classes
+---logs
+---scripts
|       build_and_run.ps1
|       Linker_Fixer.ps1
|       project_health_report.ps1
|       update_cmake.ps1
|       
+---server
|   |   mathengine.db
|   |   mathengine.db-shm
|   |   mathengine.db-wal
|   |   pom.xml
|   |   
|   +---.vscode
|   |       settings.json
|   |       
|   +---src
|   |   \---main
|   |       +---java
|   |       |   \---com
|   |       |       \---mathengine
|   |       |           \---com
|   |       |               \---mathengine
|   |       |                   \---server
|   |       |                       |   MathEngineServer.java
|   |       |                       |   
|   |       |                       +---api
|   |       |                       |       AuthController.java
|   |       |                       |       ComputeController.java
|   |       |                       |       HistoryController.java
|   |       |                       |       
|   |       |                       +---auth
|   |       |                       |       JwtFilter.java
|   |       |                       |       JwtService.java
|   |       |                       |       SecurityConfig.java
|   |       |                       |       SecurityUtils.java
|   |       |                       |       
|   |       |                       +---db
|   |       |                       |       DatabaseManager.java
|   |       |                       |       
|   |       |                       +---engine
|   |       |                       |       ServerEngineService.java
|   |       |                       |       
|   |       |                       \---model
|   |       |                               Models.java
|   |       |                               
|   |       \---resources
|   |           |   application.properties
|   |           |   
|   |           \---static
|   |               |   index.html
|   |               |   mobile.css
|   |               |   mobile.js
|   |               |   
|   |               +---css
|   |               \---js
|   \---target
|       +---classes
|       |   |   application.properties
|       |   |   
|       |   +---com
|       |   |   \---mathengine
|       |   |       \---server
|       |   |           |   MathEngineServer.class
|       |   |           |   
|       |   |           +---api
|       |   |           |       AuthController.class
|       |   |           |       ComputeController.class
|       |   |           |       HistoryController.class
|       |   |           |       
|       |   |           +---auth
|       |   |           |       JwtFilter.class
|       |   |           |       JwtService.class
|       |   |           |       SecurityConfig.class
|       |   |           |       SecurityUtils.class
|       |   |           |       
|       |   |           +---db
|       |   |           |       DatabaseManager.class
|       |   |           |       
|       |   |           +---engine
|       |   |           |       ServerEngineService.class
|       |   |           |       
|       |   |           +---model
|       |   |           |       Models$AuthResponse.class
|       |   |           |       Models$ComputeRequest.class
|       |   |           |       Models$ComputeResponse.class
|       |   |           |       Models$ExportRequest.class
|       |   |           |       Models$ExportResponse.class
|       |   |           |       Models$HistoryEntry.class
|       |   |           |       Models$LoginRequest.class
|       |   |           |       Models$Preferences.class
|       |   |           |       Models$RegisterRequest.class
|       |   |           |       Models$User.class
|       |   |           |       Models.class
|       |   |           |       
|       |   |           \---service
|       |   \---static
|       |           index.html
|       |           mobile.css
|       |           mobile.js
|       |           
|       +---generated-sources
|       |   \---annotations
|       \---maven-status
|           \---maven-compiler-plugin
|               \---compile
|                   \---default-compile
|                           createdFiles.lst
|                           inputFiles.lst
|                           
+---src
|   \---main
|       \---cpp
|           +---build
|           |   \---Debug
|           |           mathengine.pdb
|           |           
|           \---parser
|               \---include_week_four
|                       evaluator.hpp
|                       lexer.hpp
|                       parser.hpp
|                       value.hpp
|                       
\---target
    +---classes
    |   +---com
    |   |   \---mathengine
    |   |       +---jni
    |   |       |       MathBridge$1.class
    |   |       |       MathBridge.class
    |   |       |       
    |   |       +---model
    |   |       |       PrecisionMode.class
    |   |       |       
    |   |       \---ui
    |   |               AccessibilityToolbar.class
    |   |               AuthService$AuthResult.class
    |   |               AuthService$ServerResp.class
    |   |               AuthService.class
    |   |               CalcExpressionBuilder$1.class
    |   |               CalcExpressionBuilder.class
    |   |               DatasetImportPanel$Spacer.class
    |   |               DatasetImportPanel.class
    |   |               GraphPanel$PlotEntry.class
    |   |               GraphPanel.class
    |   |               InputPanel$Spacer.class
    |   |               InputPanel.class
    |   |               LatexRenderer.class
    |   |               Launcher.class
    |   |               LoginScreen.class
    |   |               MainLayout$1.class
    |   |               MainLayout$Spacer.class
    |   |               MainLayout.class
    |   |               MathEngineApp.class
    |   |               MathKeyboard$Key.class
    |   |               MathKeyboard.class
    |   |               OperationPanel$HistoryEntry.class
    |   |               OperationPanel$Operation.class
    |   |               OperationPanel$Spacer.class
    |   |               OperationPanel.class
    |   |               OutputPanel$Spacer.class
    |   |               OutputPanel.class
    |   |               QuickFunctionsPanel.class
    |   |               regUser.class
    |   |               SettingsPanel.class
    |   |               ThemeManager$FontSize.class
    |   |               ThemeManager$Theme.class
    |   |               ThemeManager.class
    |   |               
    |   \---css
    |           base.css
    |           dark.css
    |           font-large.css
    |           font-normal.css
    |           font-x-large.css
    |           high-contrast.css
    |           light.css
    |           
    +---generated-sources
    |   \---annotations
    \---maven-status
        \---maven-compiler-plugin
            \---compile
                \---default-compile
                        createdFiles.lst
                        inputFiles.lst
```

---

## Prerequisites

| Tool | Version | Notes |
| ------ | --------- | ------- |
| JDK | 17 or higher | JAVA_HOME must be set, or java.exe on PATH |
| Maven | 3.8 or higher | mvn on PATH |
| Visual Studio | 2019 or higher | "Desktop development with C++" workload |
| Git | any | optional, for cloning |

The build script auto-detects and initializes the MSVC environment via `vcvarsall.bat` --
you do not need to run from a Developer PowerShell.

---

## Quick Start

### Step 1 -- Copy your week_4 parser headers

Copy these four files from your `week_4/include_week_four/` directory into
`src/main/cpp/parser/include_week_four/`:

```txt
evaluator.hpp
lexer.hpp
parser.hpp
value.hpp
```

The three `.cpp` source files are already present with their `#include` paths corrected.

### Step 2 -- Build the C++ engine and launch the desktop app

Open PowerShell in the project root and run:

```powershell
.\build.ps1
```

This compiles Java, generates the JNI header, compiles `mathengine.dll`, and launches
the JavaFX desktop application.

### Step 3 -- Start the server for phone and tablet access (optional)

In a second PowerShell window:

```powershell
.\build.ps1 -StartServer
```

The server prints something like:

```txt
========================================
  Math Engine Server is running!
----------------------------------------
  Desktop / localhost:
    http://localhost:8080

  Other devices on the same WiFi:
    http://192.168.1.45:8080
========================================
```

Open that address in a browser on any phone or tablet connected to the same WiFi network.

---

## Build Script Reference

All flags can be combined.

| Command | What it does |
| --------- | ------------- |
| `.\build.ps1` | Full build (Java + C++) then launch desktop app |
| `.\build.ps1 -StartServer` | Full build then start Spring Boot server |
| `.\build.ps1 -SkipJavaBuild` | Skip `mvn compile`, rebuild DLL only, then run |
| `.\build.ps1 -SkipCppBuild` | Skip DLL build, recompile Java only, then run |
| `.\build.ps1 -RunOnly` | Skip all builds, launch desktop app immediately |
| `.\build.ps1 -RunOnly -StartServer` | Skip all builds, start server immediately |
| `.\build.ps1 -Compiler mingw` | Force MinGW g++ instead of MSVC |
| `.\build.ps1 -VerboseOutput` | Print every command before running it |

The script validates all prerequisites before touching any build step, and prints
a clear error message if something is missing (e.g. headers not copied, DLL not found).

---

## REST API Reference

All endpoints are served on port `8080`. Authenticated endpoints require the header:

```txt
Authorization: Bearer <token>
```

Tokens are returned by `/api/auth/register` and `/api/auth/login`, and are valid for 72 hours.

### Auth

| Method | Endpoint | Auth | Description |
| -------- | ---------- | ------ | ------------- |
| POST | `/api/auth/register` | None | Create a new account |
| POST | `/api/auth/login` | None | Sign in, receive token |
| GET | `/api/auth/me` | Required | Get current user profile |

**Register / Login request body:**

```json
{ "username": "marle", "email": "marle@example.com", "password": "mypassword" }
```

**Register / Login response:**

```json
{ "token": "eyJ...", "username": "marle", "theme": "dark", "fontSize": "normal" }
```

### Compute

| Method | Endpoint | Auth | Description |
| -------- | ---------- | ------ | ------------- |
| POST | `/api/compute` | Optional | Evaluate an expression. History saved if signed in. |
| GET | `/api/engine/status` | None | Check if the C++ native library is loaded |

**Request body:**

```json
{ "expression": "\\int_0^8 x^2 dx", "precisionFlag": 0, "operation": "integral" }
```

`precisionFlag`: `0` = symbolic, `1` = numerical

**Response:**

```json
{
  "ok": true,
  "expression": "\\int_0^8 x^2 dx",
  "result": "512/3  ~  170.6666666667",
  "operation": "integral",
  "precisionMode": "symbolic"
}
```

Results containing `~` include both a symbolic and a numeric form.

### History

| Method | Endpoint | Auth | Description |
| -------- | ---------- | ------ | ------------- |
| GET | `/api/history?limit=50` | Required | Fetch recent history (max 200) |
| DELETE | `/api/history` | Required | Clear all history |
| DELETE | `/api/history/{id}` | Required | Delete one entry |

### Preferences

| Method | Endpoint | Auth | Description |
| -------- | ---------- | ------ | ------------- |
| GET | `/api/preferences` | Required | Load theme and font size |
| PUT | `/api/preferences` | Required | Save theme and font size |

**Request / response body:**

```json
{ "theme": "dark", "fontSize": "normal" }
```

Valid theme values: `light`, `dark`, `system`, `high-contrast`
Valid fontSize values: `normal`, `large`, `x-large`

### Export

| Method | Endpoint | Auth | Description |
| -------- | ---------- | ------ | ------------- |
| POST | `/api/export` | None | Convert result to an export format |

**Request body:**

```json
{ "expression": "x^2", "result": "512/3", "format": "latex" }
```

**Supported formats:**

| Format | Output example |
| -------- | ---------------- |
| `latex` | `\[ x^2 = 512/3 \]` |
| `unicode` | `x² = 512/3` |
| `plaintext` | `x^2 = 512/3` |
| `mathml` | `<math xmlns="...">...</math>` |

---

## Supported Expressions

Everything the week_4 parser handles, entered as LaTeX or plain text:

```tex
-- Arithmetic
2 + 3 * 4
sin(pi/2)
cos(pi/3) + sin(pi/6)
sqrt(2^2 + 3^2)
ln(exp(5))
2 ^ 8
-(2 + 3)
pi                    -> symbolic: "pi  ~  3.1415926536"
e                     -> symbolic: "e   ~  2.7182818285"

-- Calculus
\int_0^8 x^2\,dx
\frac{d}{dx}[x^3 + sin(x)]
\lim_{x \to 0} \frac{sin(x)}{x}

-- Matrix
matrix:[[1,2],[3,4]]
```

### Precision modes

| Mode | Toggle | C++ exactMode | Example output for "pi" |
| ------ | -------- | -------------- | ------------------------ |
| Symbolic | left pill | `true` | `pi  ~  3.1415926536` |
| Numerical | right pill | `false` | `3.1415926536` |

---

## Accounts and Sync

Accounts are entirely optional. Guest users can compute, export, and use all five
operation types without signing in.

**Signing in unlocks:**

- **History sync** -- every computation is saved to the server database and visible on all
  devices signed in to the same account
- **Preference sync** -- theme and font size choices are saved server-side and restored
  automatically on next login from any device
- **Cross-device continuity** -- start a calculation on your phone, pick it up from history
  on your desktop

Passwords are hashed with BCrypt (strength 12). Tokens are HMAC-SHA256 signed JWTs
with a 72-hour expiry. The database is a single SQLite file (`mathengine.db`) in the
project root, created automatically on first server start.

---

## Accessibility Features

All features apply to both the desktop app and the mobile browser UI.

| Feature | Desktop (JavaFX) | Mobile (browser) |
| --------- | ----------------- | ----------------- |
| Light theme | Yes | Yes |
| Dark theme | Yes | Yes |
| System theme (follows OS) | Yes | Yes |
| High Contrast (WCAG AAA) | Yes | Yes |
| Font size: Normal (13/15px) | Yes | Yes |
| Font size: Large (16/18px) | Yes | Yes |
| Font size: Extra-large (20/22px) | Yes | Yes |
| Full keyboard navigation | Yes | Yes |
| Screen reader labels (ARIA) | Yes | Yes |
| Theme persisted across restarts | Yes (java.util.prefs) | Yes (localStorage + server sync) |
| Reduced motion support | -- | Yes (prefers-reduced-motion) |
| iOS safe-area insets | -- | Yes (env(safe-area-inset-*)) |

---

## Configuration

All server configuration is in `server/src/main/resources/application.properties`.

```properties
# Port the server listens on
server.port=8080

# Bind address -- 0.0.0.0 lets all devices on the network connect
# Change to 127.0.0.1 to restrict to localhost only
server.address=0.0.0.0

# SQLite database file path (relative to project root)
mathengine.db.path=mathengine.db

# JWT secret -- change this before any real deployment
# Must be at least 32 characters
mathengine.jwt.secret=MathEngine-ChangeThis-Secret-Key-Min32Chars!!

# How long tokens stay valid
mathengine.jwt.expiry-hours=72

# CORS -- use * for local network development
# Set to specific origin(s) before deploying publicly
mathengine.cors.allowed-origins=*
```

---

## Troubleshooting

**`FAIL  No C++ compiler found`**
The build script searches for MSVC via `vswhere.exe`. Make sure Visual Studio is installed
with the "Desktop development with C++" workload. MinGW (`g++` on PATH) also works as a fallback.

**`FAIL  Cannot build DLL without the parser headers`**
Copy `evaluator.hpp`, `lexer.hpp`, `parser.hpp`, and `value.hpp` from your `week_4/include_week_four/`
folder into `src/main/cpp/parser/include_week_four/`.

**`WARNING: Native library not found -- running in Java stub mode`**
The server started but cannot find `mathengine.dll`. Run `.\build.ps1 -SkipJavaBuild` to
rebuild the DLL, or run `.\build.ps1 -StartServer` which builds everything first.

### Phone cannot reach the server

- Confirm the phone and PC are on the same WiFi network
- Check Windows Firewall: allow inbound connections on port 8080
- Use the IP address printed by the server on startup, not `localhost`

**`java.library.path` errors at runtime**
The build script sets this automatically. If running Maven directly, add:

```powershell
-Djava.library.path=src/main/cpp/build
```

### JNI handshake test

```powershell
cd mathengine
mvn test -Dtest=MathBridgeTest
```

Sends `PING` from Java to C++, expects `PING::PONG` back. If this passes, the full
engine pipeline is wired correctly.

## Known Errors and Bugs

| Error ID     | Description |
|--------------|-------------|
| ER00001      | In the Derivative quick function `curl[y,-x,0,x,y,z]`, users receive Computation Error: Unexpected token after expression |
| ER00002      | In the Derivative quick function `jacobian[x^2+y,x*y,x,y]`, users receive Computation Error: Unexpected token after expression |
| ER00003      | In the Integral quick function `romberg[sin(x)|x|0|3.14159]`, users receive Computation Error: Unexpected character: $|$ |
| ER00004      | In the Integral quick function `numerical_int[exp(-x^2)|x|-5|5]`, users receive Computation Error: Unexpected character: "_" |
| ER00005      | In the Vector quick function `line[x^2+y^2,x,y,t^2,t,0,1]`, users receive Computation Error: Undefined variable: y |
| ER00006      | In the Vector quick function `surface[x+y+z,x,y,z]`, users receive Computation Error: Undefined variable: z |
| ER00007      | In the Vector quick function `greens[y,-x,0,1,0,1]`, users receive Computation Error: Undefined variable: z |
| ER00008      | In the Vector quick function `stokes[y,-x,x,y,z]`, users receive Computation Error: Undefined variable: z |
| ER00009      | In the Optimize quick function `saddle:x^2-y^2`, users receive Computation Error: Undefined variable: y |
| ER00010      | In the Optimize quick function `gradient_descent: x^2+y^2,-5,5` users receive Computation Error: Undefined variable: y |
| ER00011      | In the Probability quick function `markov[0.7,0.3;0.4,0.5|0.6,0.4|5]`, users receive Computation Error: Unknown probability theory operation: markov |
| ER00012      | In the Numerical Analysis quick function `spline[0,1,2,3|0,1,4,9|1.5]`, users receive Computation Error: cubic spline requires at least 3 data points |
| ER00013      | In the Linear Algebra quick function `la:null|[[1,2,3],[4,5,6],[7,8,9]]`, users receive Computation Error: Unknown linear algebra operation: null |
| ER00014      | In the Linear Algebra quick function `la:qr|[[1,1],[1,-1],[0,1]]`, users receive Computation Error: Unknown linear algebra operation: qr |
| ER00015      | In the Linear Algebra quick function `la:lu|[[2,1,-1],[4,3,-3],[2,-1,2]]`, users receive Computation Error: Unknown linear algebra operation: lu |
| ER00016      | In the Linear Algebra quick function `la:svd|[[1,2],[3,4],[5,6]]`, users receive Computation Error: Unknown linear algebra operation: svd |
| ER00017      | In the Linear Algebra quick function `la:eigen|[[4,1],[2,3]]`, users receive Computation Error: Unknown linear algebra operation: eigen |
| ER00018      | In the Linear Algebra quick function `la:inv|[[2,1],[`, users receive Computation Error: Unknown linear algebra operation: inv |

| Bug ID | Description |
|--------|-------------|
| B00001 | Dot product and Convex Hull are known to utilize the `circle[]` operation instead of their respective operations |
| B00002 | Outputted fractions do not seem to be reduced to their simplest form |
| B00003 | Text output is in 'math mode' but shouldn't be |

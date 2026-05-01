cmake_minimum_required(VERSION 3.16)
project(mathengine CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ── Find JDK (required for jni.h) ─────────────────────────────────────────────
find_package(JNI REQUIRED)

# ── Calculus module (static lib) ──────────────────────────────────────────────
add_subdirectory(desktop/src/main/cpp/Calculus)

# ── Source File Definitions ───────────────────────────────────────────────────
set(SOURCES
    desktop/src/main/cpp/CoreEngine.cpp
    desktop/src/main/cpp/MathBridgeJNI.cpp
    desktop/src/main/cpp/CommonUtils.cpp
    # Check this filename carefully! Changed Preprocssor -> Preprocessor
    desktop/src/main/cpp/LaTeX_Parser/LaTeXPreprocessor.cpp
    desktop/src/main/cpp/AbstractAlgebra/AA.cpp
    desktop/src/main/cpp/AppliedMath/AM.cpp
    desktop/src/main/cpp/Calculus/Calculus.cpp
    desktop/src/main/cpp/Calculus/core/Expression.cpp
    desktop/src/main/cpp/Calculus/core/Simplify.cpp
    desktop/src/main/cpp/Calculus/core/Week4Bridge.cpp
    desktop/src/main/cpp/Calculus/differentiation/Derivative.cpp
    desktop/src/main/cpp/Calculus/differentiation/Implicit.cpp
    desktop/src/main/cpp/Calculus/differentiation/Partial.cpp
    desktop/src/main/cpp/Calculus/integration/Multivariable.cpp
    desktop/src/main/cpp/Calculus/integration/Numerical.cpp
    desktop/src/main/cpp/Calculus/integration/Symbolic.cpp
    desktop/src/main/cpp/Calculus/limits/Limits.cpp
    desktop/src/main/cpp/Calculus/optimization/Optimize.cpp
    desktop/src/main/cpp/Calculus/series/Convergence.cpp
    desktop/src/main/cpp/Calculus/series/Series.cpp
    desktop/src/main/cpp/Calculus/vectorcalc/LineIntegral.cpp
    desktop/src/main/cpp/Calculus/vectorcalc/SurfaceIntegral.cpp
    desktop/src/main/cpp/Calculus/vectorcalc/Theorems.cpp
    desktop/src/main/cpp/Calculus/vectorcalc/VectorOps.cpp
    desktop/src/main/cpp/ComplexAnalysis/CA.cpp
    desktop/src/main/cpp/DiffEq/DE.cpp
    desktop/src/main/cpp/DiffEq/DEAdv.cpp
    desktop/src/main/cpp/DiscreteMath/DM.cpp  # Fixed typo: DIscrete -> Discrete
    desktop/src/main/cpp/Geometry/Geom.cpp
    desktop/src/main/cpp/Linear_Algebra/LA.cpp
    desktop/src/main/cpp/NumberTheory/NumberTheory.cpp
    desktop/src/main/cpp/NumericalAnalysis/NA.cpp
    desktop/src/main/cpp/parser/src_week_four/evaluator.cpp
    desktop/src/main/cpp/parser/src_week_four/lexer.cpp
    desktop/src/main/cpp/parser/src_week_four/parser.cpp
    desktop/src/main/cpp/ProbabilityTheory/PT.cpp
    desktop/src/main/cpp/Statistics/Advanced_Statistics.cpp
    desktop/src/main/cpp/Statistics/Statistics.cpp
)

# ── Build the shared library ──────────────────────────────────────────────────
add_library(mathengine SHARED ${SOURCES})

target_include_directories(mathengine PRIVATE
    ${JNI_INCLUDE_DIRS}
    desktop/src/main/cpp/include
    desktop/src/main/cpp/parser/include_week_four
)

# ── Link the calculus static lib ──────────────────────────────────────────────
target_link_libraries(mathengine PRIVATE calculus)

# ── Output Configuration ──────────────────────────────────────────────────────
set_target_properties(mathengine PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/src/main/cpp/build"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/src/main/cpp/build"
)
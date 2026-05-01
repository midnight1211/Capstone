package com.mathengine.server.engine;

import org.springframework.beans.factory.annotation.Value;
import org.springframework.stereotype.Service;

import jakarta.annotation.PostConstruct;

/**
 * ServerEngineService
 * ────────────────────
 * Server-side wrapper around the C++ math engine, loaded via JNI.
 *
 * This is the single place in the entire server where the native
 * library is loaded. All REST controllers route through this service.
 *
 * The native library (mathengine.dll / libmathengine.so) must be
 * on java.library.path, which the build_and_run.ps1 script sets
 * automatically via -Djava.library.path=src/main/cpp/build.
 *
 * If the library fails to load (e.g. on a fresh clone before the
 * C++ has been compiled), the service falls back to a Java stub
 * so the server still starts and the REST API still responds.
 */
@Service
public class ServerEngineService {

    private boolean nativeLoaded = false;

    // ── JNI declarations (mirror com.mathengine.jni.MathBridge) ──────────────

    private native String nativeCompute(String expression, int precisionFlag);
    private native String nativeHandshake(String ping);
    private native String nativeVersion();

    // ── Initialisation ────────────────────────────────────────────────────────

    @PostConstruct
    public void init() {
        try {
            System.loadLibrary("mathengine");
            nativeLoaded = true;
            System.out.println("[Engine] Native library loaded successfully");
            System.out.println("[Engine] Version: " + nativeVersion());
        } catch (UnsatisfiedLinkError e) {
            System.err.println("[Engine] WARNING: Native library not found — " +
                               "running in Java stub mode");
            System.err.println("[Engine] Run build_and_run.ps1 to compile the C++ engine");
        }
    }

    // ── Public API ────────────────────────────────────────────────────────────

    /**
     * Evaluate a math expression.
     *
     * @param expression   LaTeX or plain-text expression from the user
     * @param precisionFlag  0 = symbolic, 1 = numerical
     * @return  Result string (may include "  ~  " separator for sym + approx)
     * @throws  RuntimeException if the C++ engine throws
     */
    public String compute(String expression, int precisionFlag) {
        if (nativeLoaded) {
            return nativeCompute(expression, precisionFlag);
        }
        return stubCompute(expression, precisionFlag);
    }

    public boolean isNativeLoaded() { return nativeLoaded; }

    public String getVersion() {
        return nativeLoaded ? nativeVersion() : "Java stub mode";
    }

    // ── Java stub (used when native lib is absent) ────────────────────────────

    private String stubCompute(String expression, int precisionFlag) {
        boolean symbolic = (precisionFlag == 0);
        String e = expression.toLowerCase().trim();

        if (e.equals("pi") || e.equals("\\pi"))
            return symbolic ? "pi  ~  3.1415926536" : "3.1415926536";
        if (e.equals("e"))
            return symbolic ? "e  ~  2.7182818285" : "2.7182818285";
        if (e.contains("sqrt(2)"))
            return symbolic ? "sqrt(2)  ~  1.4142135624" : "1.4142135624";
        if (e.contains("sin(pi/2)"))
            return symbolic ? "sin(pi/2) = 1" : "1";
        if (e.contains("2^8") || e.contains("2 ^ 8"))
            return "256";

        return "[stub] " + expression + " = ?  (compile C++ engine for real results)";
    }
}
package com.mathengine.server.api;

import com.mathengine.server.auth.SecurityUtils;
import com.mathengine.server.db.DatabaseManager;
import com.mathengine.server.engine.ServerEngineService;
import com.mathengine.server.model.Models;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;

import java.sql.SQLException;
import java.util.Optional;

/**
 * ComputeController
 * ──────────────────
 * POST /api/compute
 *
 * Accepts a math expression and precision flag, passes it to the
 * C++ engine via ServerEngineService, and returns the result.
 *
 * If the request carries a valid JWT, the result is automatically
 * saved to the user's history in SQLite.
 *
 * Guest requests (no JWT) work fine — history just isn't persisted.
 */
@RestController
@RequestMapping("/api")
public class ComputeController {

    private final ServerEngineService engine;
    private final DatabaseManager     db;

    public ComputeController(ServerEngineService engine, DatabaseManager db) {
        this.engine = engine;
        this.db     = db;
    }

    @PostMapping("/compute")
    public ResponseEntity<?> compute(@RequestBody Models.ComputeRequest req) {

        // Validate
        if (req.expression == null || req.expression.isBlank())
            return ResponseEntity.badRequest().body(
                Models.ComputeResponse.failure("", "expression is required"));

        if (req.precisionFlag != 0 && req.precisionFlag != 1)
            return ResponseEntity.badRequest().body(
                Models.ComputeResponse.failure(req.expression,
                    "precisionFlag must be 0 (symbolic) or 1 (numerical)"));

        String operation = req.operation != null ? req.operation : "arithmetic";
        String modeName  = req.precisionFlag == 0 ? "symbolic" : "numerical";

        try {
            // Call C++ engine
            String result = engine.compute(req.expression, req.precisionFlag);

            // If user is authenticated, save to history
            Optional<Long> uid = SecurityUtils.currentUserId();
            if (uid.isPresent()) {
                try {
                    db.saveHistory(uid.get(), req.expression, result,
                                   operation, modeName);
                } catch (SQLException e) {
                    // History save failure is non-fatal — still return the result
                    System.err.println("[Compute] Failed to save history: " + e.getMessage());
                }
            }

            return ResponseEntity.ok(
                Models.ComputeResponse.success(req.expression, result,
                                               operation, modeName));

        } catch (Exception e) {
            return ResponseEntity.status(500).body(
                Models.ComputeResponse.failure(req.expression, e.getMessage()));
        }
    }

    /** GET /api/engine/status — health check for the desktop client */
    @GetMapping("/engine/status")
    public ResponseEntity<?> engineStatus() {
        return ResponseEntity.ok(java.util.Map.of(
            "nativeLoaded", engine.isNativeLoaded(),
            "version",      engine.getVersion(),
            "status",       "ok"
        ));
    }
}
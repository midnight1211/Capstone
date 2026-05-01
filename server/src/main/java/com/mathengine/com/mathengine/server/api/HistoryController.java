package com.mathengine.server.api;

import com.mathengine.server.auth.SecurityUtils;
import com.mathengine.server.db.DatabaseManager;
import com.mathengine.server.model.Models;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;

import java.sql.SQLException;
import java.util.List;
import java.util.Map;
import java.util.Optional;

/**
 * HistoryController
 * ──────────────────
 * All endpoints require a valid JWT (enforced in SecurityConfig).
 *
 * GET    /api/history            — fetch recent history (default last 50)
 * DELETE /api/history            — clear all history for the user
 * DELETE /api/history/{id}       — delete a single entry
 *
 * GET    /api/preferences        — load theme + font size
 * PUT    /api/preferences        — save theme + font size
 */
@RestController
@RequestMapping("/api")
public class HistoryController {

    private final DatabaseManager db;

    public HistoryController(DatabaseManager db) {
        this.db = db;
    }

    // ── History ───────────────────────────────────────────────────────────────

    @GetMapping("/history")
    public ResponseEntity<?> getHistory(
            @RequestParam(defaultValue = "50") int limit) {

        Optional<Long> uid = SecurityUtils.currentUserId();
        if (uid.isEmpty()) return unauthorized();

        try {
            List<Models.HistoryEntry> entries = db.getHistory(uid.get(),
                Math.min(limit, 200)); // cap at 200
            return ResponseEntity.ok(entries);
        } catch (SQLException e) {
            return serverError(e.getMessage());
        }
    }

    @DeleteMapping("/history")
    public ResponseEntity<?> clearHistory() {
        Optional<Long> uid = SecurityUtils.currentUserId();
        if (uid.isEmpty()) return unauthorized();

        try {
            db.clearHistory(uid.get());
            return ResponseEntity.ok(Map.of("cleared", true));
        } catch (SQLException e) {
            return serverError(e.getMessage());
        }
    }

    @DeleteMapping("/history/{id}")
    public ResponseEntity<?> deleteEntry(@PathVariable long id) {
        Optional<Long> uid = SecurityUtils.currentUserId();
        if (uid.isEmpty()) return unauthorized();

        try {
            db.deleteHistoryEntry(id, uid.get());
            return ResponseEntity.ok(Map.of("deleted", true));
        } catch (SQLException e) {
            return serverError(e.getMessage());
        }
    }

    // ── Preferences ───────────────────────────────────────────────────────────

    @GetMapping("/preferences")
    public ResponseEntity<?> getPreferences() {
        Optional<Long> uid = SecurityUtils.currentUserId();
        if (uid.isEmpty()) return unauthorized();

        try {
            var user = db.findById(uid.get());
            if (user.isEmpty()) return unauthorized();
            return ResponseEntity.ok(
                new Models.Preferences(user.get().theme, user.get().fontSize));
        } catch (SQLException e) {
            return serverError(e.getMessage());
        }
    }

    @PutMapping("/preferences")
    public ResponseEntity<?> savePreferences(
            @RequestBody Models.Preferences prefs) {

        Optional<Long> uid = SecurityUtils.currentUserId();
        if (uid.isEmpty()) return unauthorized();

        // Validate values
        if (!List.of("light","dark","system","high-contrast")
                 .contains(prefs.theme))
            return ResponseEntity.badRequest()
                .body(Map.of("error", "invalid theme value"));
        if (!List.of("normal","large","x-large")
                 .contains(prefs.fontSize))
            return ResponseEntity.badRequest()
                .body(Map.of("error", "invalid fontSize value"));

        try {
            db.updatePreferences(uid.get(), prefs.theme, prefs.fontSize);
            return ResponseEntity.ok(prefs);
        } catch (SQLException e) {
            return serverError(e.getMessage());
        }
    }

    // ── Helpers ───────────────────────────────────────────────────────────────

    private ResponseEntity<Map<String,String>> unauthorized() {
        return ResponseEntity.status(401).body(Map.of("error", "Authentication required"));
    }

    private ResponseEntity<Map<String,String>> serverError(String msg) {
        return ResponseEntity.status(500).body(Map.of("error", msg));
    }
}
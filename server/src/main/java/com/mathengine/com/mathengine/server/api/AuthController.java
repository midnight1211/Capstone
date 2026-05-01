package com.mathengine.server.api;

import com.mathengine.server.auth.JwtService;
import com.mathengine.server.auth.SecurityUtils;
import com.mathengine.server.db.DatabaseManager;
import com.mathengine.server.model.Models;
import org.springframework.http.ResponseEntity;
import org.springframework.security.crypto.password.PasswordEncoder;
import org.springframework.web.bind.annotation.*;

import java.sql.SQLException;
import java.util.Map;
import java.util.Optional;

/**
 * AuthController
 * ───────────────
 * Handles user registration and login.
 *
 * POST /api/auth/register   — create account, returns JWT
 * POST /api/auth/login      — authenticate, returns JWT
 * GET  /api/auth/me         — returns current user info (requires JWT)
 */
@RestController
@RequestMapping("/api/auth")
public class AuthController {

    private final DatabaseManager db;
    private final JwtService      jwtService;
    private final PasswordEncoder passwordEncoder;

    public AuthController(DatabaseManager db, JwtService jwtService,
                          PasswordEncoder passwordEncoder) {
        this.db              = db;
        this.jwtService      = jwtService;
        this.passwordEncoder = passwordEncoder;
    }

    // ── POST /api/auth/register ───────────────────────────────────────────────

    @PostMapping("/register")
    public ResponseEntity<?> register(@RequestBody Models.RegisterRequest req) {
        // Validate input
        if (req.username == null || req.username.isBlank())
            return badRequest("username is required");
        if (req.email == null || req.email.isBlank())
            return badRequest("email is required");
        if (req.password == null || req.password.length() < 8)
            return badRequest("password must be at least 8 characters");
        if (req.username.length() > 32)
            return badRequest("username must be 32 characters or fewer");
        if (!req.email.contains("@"))
            return badRequest("email address is not valid");

        try {
            // Check for duplicates
            if (db.findByUsername(req.username).isPresent())
                return badRequest("Username already taken");
            if (db.findByEmail(req.email).isPresent())
                return badRequest("Email already registered");

            // Hash password and create user
            String hash   = passwordEncoder.encode(req.password);
            long   userId = db.createUser(req.username, req.email, hash);

            String token = jwtService.generateToken(userId, req.username);
            return ResponseEntity.ok(
                new Models.AuthResponse(token, req.username, "light", "normal"));

        } catch (SQLException e) {
            return serverError("Registration failed: " + e.getMessage());
        }
    }

    // ── POST /api/auth/login ──────────────────────────────────────────────────

    @PostMapping("/login")
    public ResponseEntity<?> login(@RequestBody Models.LoginRequest req) {
        if (req.username == null || req.password == null)
            return badRequest("username and password are required");

        try {
            Optional<Models.User> userOpt = db.findByUsername(req.username);
            if (userOpt.isEmpty())
                return unauthorized("Invalid username or password");

            Models.User user = userOpt.get();
            if (!passwordEncoder.matches(req.password, user.passwordHash))
                return unauthorized("Invalid username or password");

            String token = jwtService.generateToken(user.id, user.username);
            return ResponseEntity.ok(
                new Models.AuthResponse(token, user.username,
                                        user.theme, user.fontSize));

        } catch (SQLException e) {
            return serverError("Login failed: " + e.getMessage());
        }
    }

    // ── GET /api/auth/me ──────────────────────────────────────────────────────

    @GetMapping("/me")
    public ResponseEntity<?> me() {
        Optional<Long> uid = SecurityUtils.currentUserId();
        if (uid.isEmpty()) return unauthorized("Not authenticated");

        try {
            Optional<Models.User> user = db.findById(uid.get());
            if (user.isEmpty()) return unauthorized("User not found");

            Models.User u = user.get();
            // Return profile without the password hash
            return ResponseEntity.ok(Map.of(
                "id",       u.id,
                "username", u.username,
                "email",    u.email,
                "theme",    u.theme,
                "fontSize", u.fontSize,
                "since",    u.createdAt.toString()
            ));
        } catch (SQLException e) {
            return serverError(e.getMessage());
        }
    }

    // ── Helpers ───────────────────────────────────────────────────────────────

    private ResponseEntity<Map<String,String>> badRequest(String msg) {
        return ResponseEntity.badRequest().body(Map.of("error", msg));
    }

    private ResponseEntity<Map<String,String>> unauthorized(String msg) {
        return ResponseEntity.status(401).body(Map.of("error", msg));
    }

    private ResponseEntity<Map<String,String>> serverError(String msg) {
        return ResponseEntity.status(500).body(Map.of("error", msg));
    }
}
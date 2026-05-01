package com.mathengine.server.db;

import com.mathengine.server.model.Models;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.stereotype.Component;

import jakarta.annotation.PostConstruct;
import java.sql.*;
import java.time.Instant;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

/**
 * DatabaseManager
 * ────────────────
 * Manages the SQLite connection and all database operations.
 *
 * Uses raw JDBC — no ORM. This keeps the dependency count low
 * and makes the SQL explicit and easy to follow for a capstone project.
 *
 * Schema:
 *   users        — accounts with hashed passwords and preferences
 *   history      — calculation history entries per user
 *
 * The database file is created automatically at first run.
 * Path is configured via mathengine.db.path in application.properties.
 *
 * Thread safety: SQLite in WAL mode handles concurrent reads fine.
 * We use a single connection with synchronized methods for writes.
 */
@Component
public class DatabaseManager {

    @Value("${mathengine.db.path}")
    private String dbPath;

    private Connection connection;

    // ── Initialisation ────────────────────────────────────────────────────────

    @PostConstruct
    public void init() throws SQLException {
        String url = "jdbc:sqlite:" + dbPath;
        connection = DriverManager.getConnection(url);

        // WAL mode: better concurrent read performance
        try (Statement s = connection.createStatement()) {
            s.execute("PRAGMA journal_mode=WAL");
            s.execute("PRAGMA foreign_keys=ON");
        }

        createSchema();
        System.out.println("[DB] SQLite database ready at: " + dbPath);
    }

    private void createSchema() throws SQLException {
        try (Statement s = connection.createStatement()) {

            // Users table
            s.execute("""
                CREATE TABLE IF NOT EXISTS users (
                    id            INTEGER PRIMARY KEY AUTOINCREMENT,
                    username      TEXT    NOT NULL UNIQUE COLLATE NOCASE,
                    email         TEXT    NOT NULL UNIQUE COLLATE NOCASE,
                    password_hash TEXT    NOT NULL,
                    theme         TEXT    NOT NULL DEFAULT 'light',
                    font_size     TEXT    NOT NULL DEFAULT 'normal',
                    created_at    INTEGER NOT NULL DEFAULT (strftime('%s','now'))
                )
            """);

            // History table
            s.execute("""
                CREATE TABLE IF NOT EXISTS history (
                    id             INTEGER PRIMARY KEY AUTOINCREMENT,
                    user_id        INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
                    expression     TEXT    NOT NULL,
                    result         TEXT    NOT NULL,
                    operation      TEXT    NOT NULL DEFAULT 'arithmetic',
                    precision_mode TEXT    NOT NULL DEFAULT 'symbolic',
                    created_at     INTEGER NOT NULL DEFAULT (strftime('%s','now'))
                )
            """);

            // Index for fast per-user history lookups
            s.execute("""
                CREATE INDEX IF NOT EXISTS idx_history_user
                ON history(user_id, created_at DESC)
            """);
        }
    }

    // ── User operations ───────────────────────────────────────────────────────

    /**
     * Insert a new user. Returns the new user's id.
     * Throws SQLException if username or email already exists.
     */
    public synchronized long createUser(String username, String email,
                                        String passwordHash) throws SQLException {
        String sql = """
            INSERT INTO users (username, email, password_hash)
            VALUES (?, ?, ?)
        """;
        try (PreparedStatement ps = connection.prepareStatement(
                sql, Statement.RETURN_GENERATED_KEYS)) {
            ps.setString(1, username);
            ps.setString(2, email);
            ps.setString(3, passwordHash);
            ps.executeUpdate();
            try (ResultSet keys = ps.getGeneratedKeys()) {
                if (keys.next()) return keys.getLong(1);
            }
        }
        throw new SQLException("Failed to retrieve generated user id");
    }

    public Optional<Models.User> findByUsername(String username) throws SQLException {
        String sql = "SELECT * FROM users WHERE username = ? COLLATE NOCASE";
        try (PreparedStatement ps = connection.prepareStatement(sql)) {
            ps.setString(1, username);
            try (ResultSet rs = ps.executeQuery()) {
                if (rs.next()) return Optional.of(mapUser(rs));
            }
        }
        return Optional.empty();
    }

    public Optional<Models.User> findByEmail(String email) throws SQLException {
        String sql = "SELECT * FROM users WHERE email = ? COLLATE NOCASE";
        try (PreparedStatement ps = connection.prepareStatement(sql)) {
            ps.setString(1, email);
            try (ResultSet rs = ps.executeQuery()) {
                if (rs.next()) return Optional.of(mapUser(rs));
            }
        }
        return Optional.empty();
    }

    public Optional<Models.User> findById(long id) throws SQLException {
        try (PreparedStatement ps = connection.prepareStatement(
                "SELECT * FROM users WHERE id = ?")) {
            ps.setLong(1, id);
            try (ResultSet rs = ps.executeQuery()) {
                if (rs.next()) return Optional.of(mapUser(rs));
            }
        }
        return Optional.empty();
    }

    /** Update theme + font size preferences for a user. */
    public synchronized void updatePreferences(long userId, String theme,
                                               String fontSize) throws SQLException {
        try (PreparedStatement ps = connection.prepareStatement(
                "UPDATE users SET theme = ?, font_size = ? WHERE id = ?")) {
            ps.setString(1, theme);
            ps.setString(2, fontSize);
            ps.setLong(3, userId);
            ps.executeUpdate();
        }
    }

    // ── History operations ────────────────────────────────────────────────────

    /**
     * Save a history entry for a user.
     * Returns the new entry's id.
     */
    public synchronized long saveHistory(long userId, String expression,
                                         String result, String operation,
                                         String precisionMode) throws SQLException {
        String sql = """
            INSERT INTO history (user_id, expression, result, operation, precision_mode)
            VALUES (?, ?, ?, ?, ?)
        """;
        try (PreparedStatement ps = connection.prepareStatement(
                sql, Statement.RETURN_GENERATED_KEYS)) {
            ps.setLong(1, userId);
            ps.setString(2, expression);
            ps.setString(3, result);
            ps.setString(4, operation);
            ps.setString(5, precisionMode);
            ps.executeUpdate();
            try (ResultSet keys = ps.getGeneratedKeys()) {
                if (keys.next()) return keys.getLong(1);
            }
        }
        throw new SQLException("Failed to retrieve generated history id");
    }

    /**
     * Return the most recent N history entries for a user,
     * newest first.
     */
    public List<Models.HistoryEntry> getHistory(long userId, int limit) throws SQLException {
        String sql = """
            SELECT * FROM history
            WHERE user_id = ?
            ORDER BY created_at DESC
            LIMIT ?
        """;
        List<Models.HistoryEntry> entries = new ArrayList<>();
        try (PreparedStatement ps = connection.prepareStatement(sql)) {
            ps.setLong(1, userId);
            ps.setInt(2, limit);
            try (ResultSet rs = ps.executeQuery()) {
                while (rs.next()) entries.add(mapHistory(rs));
            }
        }
        return entries;
    }

    /** Delete all history for a user. */
    public synchronized void clearHistory(long userId) throws SQLException {
        try (PreparedStatement ps = connection.prepareStatement(
                "DELETE FROM history WHERE user_id = ?")) {
            ps.setLong(1, userId);
            ps.executeUpdate();
        }
    }

    /** Delete a single history entry (only if it belongs to the user). */
    public synchronized void deleteHistoryEntry(long entryId,
                                                long userId) throws SQLException {
        try (PreparedStatement ps = connection.prepareStatement(
                "DELETE FROM history WHERE id = ? AND user_id = ?")) {
            ps.setLong(1, entryId);
            ps.setLong(2, userId);
            ps.executeUpdate();
        }
    }

    // ── Mapping helpers ───────────────────────────────────────────────────────

    private Models.User mapUser(ResultSet rs) throws SQLException {
        return new Models.User(
            rs.getLong("id"),
            rs.getString("username"),
            rs.getString("email"),
            rs.getString("password_hash"),
            rs.getString("theme"),
            rs.getString("font_size"),
            Instant.ofEpochSecond(rs.getLong("created_at"))
        );
    }

    private Models.HistoryEntry mapHistory(ResultSet rs) throws SQLException {
        return new Models.HistoryEntry(
            rs.getLong("id"),
            rs.getLong("user_id"),
            rs.getString("expression"),
            rs.getString("result"),
            rs.getString("operation"),
            rs.getString("precision_mode"),
            Instant.ofEpochSecond(rs.getLong("created_at"))
        );
    }
}
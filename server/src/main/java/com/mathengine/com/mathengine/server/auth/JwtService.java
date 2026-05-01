package com.mathengine.server.auth;

import io.jsonwebtoken.*;
import io.jsonwebtoken.security.Keys;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.stereotype.Service;

import javax.crypto.SecretKey;
import java.nio.charset.StandardCharsets;
import java.util.Date;

/**
 * JwtService
 * ───────────
 * Issues and validates JSON Web Tokens for stateless authentication.
 *
 * Tokens are signed with HMAC-SHA256 using the secret from
 * application.properties. They carry the user's id and username
 * as claims, and expire after the configured number of hours.
 *
 * Both the JavaFX desktop client and mobile browser pass the token
 * in the Authorization header:
 *   Authorization: Bearer <token>
 */
@Service
public class JwtService {

    @Value("${mathengine.jwt.secret}")
    private String secret;

    @Value("${mathengine.jwt.expiry-hours}")
    private int expiryHours;

    // ── Token generation ──────────────────────────────────────────────────────

    public String generateToken(long userId, String username) {
        Date now    = new Date();
        Date expiry = new Date(now.getTime() + expiryHours * 3_600_000L);

        return Jwts.builder()
            .subject(String.valueOf(userId))
            .claim("username", username)
            .issuedAt(now)
            .expiration(expiry)
            .signWith(getKey())
            .compact();
    }

    // ── Token validation ──────────────────────────────────────────────────────

    /**
     * Validates and parses a token.
     * Returns null if the token is invalid or expired.
     */
    public Claims validateToken(String token) {
        try {
            return Jwts.parser()
                .verifyWith(getKey())
                .build()
                .parseSignedClaims(token)
                .getPayload();
        } catch (JwtException | IllegalArgumentException e) {
            return null;
        }
    }

    public long getUserId(Claims claims) {
        return Long.parseLong(claims.getSubject());
    }

    public String getUsername(Claims claims) {
        return claims.get("username", String.class);
    }

    // ── Helpers ───────────────────────────────────────────────────────────────

    private SecretKey getKey() {
        byte[] keyBytes = secret.getBytes(StandardCharsets.UTF_8);
        return Keys.hmacShaKeyFor(keyBytes);
    }
}
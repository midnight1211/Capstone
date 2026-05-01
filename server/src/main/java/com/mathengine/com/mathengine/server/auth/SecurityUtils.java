package com.mathengine.server.auth;

import org.springframework.security.core.Authentication;
import org.springframework.security.core.context.SecurityContextHolder;
import org.springframework.security.authentication.UsernamePasswordAuthenticationToken;

import java.util.Optional;

/**
 * SecurityUtils
 * ──────────────
 * Static helpers for controllers to retrieve the authenticated user's
 * id from the Spring Security context without boilerplate.
 *
 * Usage in a controller:
 *   Optional<Long> uid = SecurityUtils.currentUserId();
 *   if (uid.isEmpty()) return ResponseEntity.status(401).build();
 */
public final class SecurityUtils {

    private SecurityUtils() {}

    /**
     * Returns the authenticated user's id, or empty if the request
     * is from a guest (no valid JWT).
     */
    public static Optional<Long> currentUserId() {
        Authentication auth = SecurityContextHolder.getContext().getAuthentication();
        if (auth instanceof UsernamePasswordAuthenticationToken token) {
            Object details = token.getDetails();
            if (details instanceof Long id) return Optional.of(id);
        }
        return Optional.empty();
    }

    /** Returns true if the current request is from an authenticated user. */
    public static boolean isAuthenticated() {
        return currentUserId().isPresent();
    }
}
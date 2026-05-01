package com.mathengine.server;

import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;

import java.net.InetAddress;
import java.net.NetworkInterface;
import java.util.Collections;

/**
 * MathEngineServer — Spring Boot entry point.
 *
 * On startup, prints the local IP address so you can point phones
 * and tablets at the right address on your WiFi network.
 *
 * Run via Maven:
 *   cd server
 *   mvn spring-boot:run
 *
 * Or via the build_and_run.ps1 script:
 *   .\build.ps1 -StartServer
 */
@SpringBootApplication
public class MathEngineServer {

    public static void main(String[] args) throws Exception {
        SpringApplication.run(MathEngineServer.class, args);
        printNetworkInfo();
    }

    /**
     * Prints all non-loopback IPv4 addresses so the user knows
     * which address to type into a phone browser.
     */
    private static void printNetworkInfo() {
        System.out.println("\n========================================");
        System.out.println("  Math Engine Server is running!");
        System.out.println("----------------------------------------");
        System.out.println("  Desktop / localhost:");
        System.out.println("    http://localhost:8080");
        System.out.println();
        System.out.println("  Other devices on the same WiFi:");

        try {
            for (NetworkInterface iface :
                    Collections.list(NetworkInterface.getNetworkInterfaces())) {
                if (iface.isLoopback() || !iface.isUp()) continue;
                for (InetAddress addr :
                        Collections.list(iface.getInetAddresses())) {
                    String ip = addr.getHostAddress();
                    // Skip IPv6 (contains colon) and loopback
                    if (!ip.contains(":") && !ip.startsWith("127.")) {
                        System.out.printf("    http://%s:8080%n", ip);
                    }
                }
            }
        } catch (Exception e) {
            System.out.println("    (Could not determine IP address)");
        }

        System.out.println("========================================\n");
    }
}
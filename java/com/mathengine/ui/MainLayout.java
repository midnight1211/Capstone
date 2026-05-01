package com.mathengine.ui;

import com.mathengine.jni.MathBridge;
import com.mathengine.model.PrecisionMode;
import javafx.application.Platform;
import javafx.concurrent.Task;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.Label;
import javafx.scene.control.ProgressIndicator;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
import javafx.scene.layout.*;
import javafx.scene.control.Button;

/**
 * MainLayout
 * ───────────
 * Top-level BorderPane that assembles all panels and wires their
 * callbacks together. This is the root node passed to the Scene.
 *
 * Layout:
 *
 *   ┌─────────────────────────────────────────────────────────┐
 *   │  HEADER: logo  ·  AccessibilityToolbar                  │  TOP
 *   ├───────────────┬─────────────────────────────────────────┤
 *   │               │  InputPanel                             │
 *   │ OperationPanel│  ─────────────────────────────────────  │  CENTER
 *   │  (sidebar)    │  OutputPanel                            │
 *   │               │                                         │
 *   ├───────────────┴─────────────────────────────────────────┤
 *   │  STATUS BAR: dot · message · version                    │  BOTTOM
 *   └─────────────────────────────────────────────────────────┘
 *
 * Solve flow:
 *   InputPanel.onSolve → Task<String> (background thread)
 *       → MathBridge.compute(expression, PrecisionMode)
 *       → OutputPanel.showResult / showError
 *       → OperationPanel.addHistory
 */
public class MainLayout extends BorderPane {

    // ── Panels ────────────────────────────────────────────────────────────────
    private final OperationPanel       opPanel      = new OperationPanel();
    private final InputPanel           inputPanel   = new InputPanel();
    private final OutputPanel          outputPanel  = new OutputPanel();
    private final GraphPanel           graphPanel   = new GraphPanel();
    private final QuickFunctionsPanel  quickPanel      = new QuickFunctionsPanel();
    private final DatasetImportPanel   datasetPanel    = new DatasetImportPanel();

    // ── Auth ──────────────────────────────────────────────────────────────────
    private final String          currentUser;
    private final AuthService     auth = AuthService.getInstance();

    // ── Status bar widgets ────────────────────────────────────────────────────
    private final Label             statusDot  = new Label("●");
    private final Label             statusText = new Label("Ready");
    private final ProgressIndicator spinner    = new ProgressIndicator();

    // ── Bridge ────────────────────────────────────────────────────────────────
    private final MathBridge    bridge       = MathBridge.getInstance();
    private final SettingsPanel settingsPanel = new SettingsPanel();

    // ── Constructor ───────────────────────────────────────────────────────────

    public MainLayout(String username) {
        this.currentUser = username;
        getStyleClass().add("main-layout");

        setTop(buildHeader());
        setLeft(opPanel);
        setCenter(buildCenter());
        setBottom(buildStatusBar());

        wireCallbacks();
        inputPanel.focusInput();
        // Owner window for DatasetImportPanel file chooser (set after stage shows)
        javafx.application.Platform.runLater(() -> {
            if (getScene() != null && getScene().getWindow() != null)
                datasetPanel.setOwnerWindow(getScene().getWindow());
        });
    }

    /** Convenience no-arg constructor (guest mode). */
    public MainLayout() { this(null); }

    // ── Build sections ────────────────────────────────────────────────────────

    private HBox buildHeader() {
        // Logo mark
        Label mark = new Label("M");
        mark.getStyleClass().addAll("logo-mark", "logo-mark-letter");
        mark.setFocusTraversable(false);

        Label name = new Label("Math");
        name.getStyleClass().add("logo-name");
        Label accent = new Label("Engine");
        accent.getStyleClass().add("logo-accent");

        HBox logo = new HBox(8, mark, name, accent);
        logo.setAlignment(Pos.CENTER_LEFT);
        logo.getStyleClass().add("logo");

        Button settingsBtn = new Button("⚙");
        settingsBtn.getStyleClass().add("settings-btn");
        settingsBtn.setTooltip(new javafx.scene.control.Tooltip("Settings"));
        settingsBtn.setOnAction(e -> settingsPanel.toggle(settingsBtn));

        // User identity / logout
        HBox userRow;
        if (currentUser != null) {
            Label userLbl = new Label("👤 " + currentUser);
            userLbl.getStyleClass().add("header-user-label");
            Button logoutBtn = new Button("Log Out");
            logoutBtn.getStyleClass().add("header-logout-btn");
            logoutBtn.setOnAction(e -> doLogout());
            userRow = new HBox(8, userLbl, logoutBtn);
        } else {
            Label guestLbl = new Label("Guest");
            guestLbl.getStyleClass().add("header-user-label");
            userRow = new HBox(guestLbl);
        }
        userRow.setAlignment(javafx.geometry.Pos.CENTER_RIGHT);

        HBox header = new HBox(8, logo, new Spacer(), userRow, settingsBtn);
        header.setAlignment(Pos.CENTER);
        header.setPadding(new Insets(12, 20, 12, 20));
        header.getStyleClass().add("app-header");
        return header;
    }

    private VBox buildCenter() {
        // Tabs: Compute | Graph
        TabPane tabs = new TabPane();
        tabs.setTabClosingPolicy(TabPane.TabClosingPolicy.UNAVAILABLE);
        tabs.getStyleClass().add("main-tabs");

        Tab computeTab = new Tab("Compute");
        datasetPanel.setVisible(false);
        datasetPanel.setManaged(false);
        VBox computeContent = new VBox(0, inputPanel, quickPanel, datasetPanel, outputPanel);
        computeContent.setPadding(new Insets(16, 20, 20, 0));
        VBox.setVgrow(outputPanel, Priority.ALWAYS);
        computeTab.setContent(computeContent);
        VBox.setMargin(inputPanel,   new Insets(0, 0, 0, 0));
        VBox.setMargin(quickPanel,   new Insets(6, 0, 4, 0));
        VBox.setMargin(datasetPanel, new Insets(0, 0, 6, 0));

        Tab graphTab = new Tab("Graph");
        VBox graphContent = new VBox(graphPanel);
        VBox.setVgrow(graphPanel, Priority.ALWAYS);
        graphContent.setPadding(new Insets(8, 12, 12, 0));
        graphTab.setContent(graphContent);

        tabs.getTabs().addAll(computeTab, graphTab);
        VBox.setVgrow(tabs, Priority.ALWAYS);

        VBox center = new VBox(tabs);
        VBox.setVgrow(tabs, Priority.ALWAYS);
        return center;
    }

    private HBox buildStatusBar() {
        statusDot.getStyleClass().add("status-dot");
        statusDot.setFocusTraversable(false);
        statusText.getStyleClass().add("status-text");

        spinner.setMaxSize(14, 14);
        spinner.setVisible(false);
        spinner.setFocusTraversable(false);

        // ADA: status updates are announced to screen readers via aria-live equivalent
        statusText.setAccessibleText("Status: " + statusText.getText());

        Label version = new Label("week_4 parser");
        version.getStyleClass().add("status-version");

        String nativeStatus = bridge.isNativeLoaded()
            ? "C++ engine loaded"
            : "Stub mode (native lib not found)";
        Label nativeLabel = new Label(nativeStatus);
        nativeLabel.getStyleClass().add("status-version");
        nativeLabel.setAccessibleText(nativeStatus);

        HBox bar = new HBox(8, statusDot, statusText, spinner,
                            new Spacer(), nativeLabel, version);
        bar.setAlignment(Pos.CENTER_LEFT);
        bar.setPadding(new Insets(6, 20, 6, 20));
        bar.getStyleClass().add("status-bar");
        return bar;
    }

    // ── Callback wiring ───────────────────────────────────────────────────────

    private void wireCallbacks() {
        // Operation selection → update badge, placeholder, quick functions, dataset panel
        opPanel.setOnOperationChanged(op -> {
            String badgeStyle = "badge-" + op.name().toLowerCase();
            inputPanel.setOperationBadge(op.label, badgeStyle);
            inputPanel.setPlaceholder(op.placeholder);
            quickPanel.setOperation(op);
            boolean isStats = op == OperationPanel.Operation.STATISTICS;
            datasetPanel.setVisible(isStats);
            datasetPanel.setManaged(isStats);
        });

        // Dataset import → load expression into input
        datasetPanel.setOnExpressionRequested(expr -> {
            opPanel.selectOperation(OperationPanel.Operation.STATISTICS);
            inputPanel.loadExpression(expr);
            inputPanel.focusInput();
        });
        datasetPanel.setOwnerWindow(null); // set after stage is shown

        // Quick function selected → switch to the right operation, then load expression
        quickPanel.setOnExpressionSelected((expr, targetOp) -> {
            if (targetOp != null && targetOp != opPanel.getCurrentOperation()) {
                // Tell the sidebar to switch — this fires onOperationChanged which
                // updates the badge, placeholder, and quickPanel itself
                opPanel.selectOperation(targetOp);
            }
            inputPanel.loadExpression(expr);
            inputPanel.focusInput();
        });

        // Init quick panel with default operation
        quickPanel.setOperation(opPanel.getCurrentOperation());

        // History item click → re-load expression into input
        opPanel.setOnHistoryItemLoaded(entry -> {
            inputPanel.loadExpression(entry.expression());
            // Restore precision toggle
            // (precision is in entry.mode() — "symbolic" or "numerical")
        });

        // Solve → dispatch to JNI on a background thread
        inputPanel.setOnSolve((expression, precisionFlag) -> {
            dispatchSolve(expression, precisionFlag);
        });
    }

    // ── Solve dispatch ────────────────────────────────────────────────────────

    /**
     * Runs MathBridge.compute on a daemon thread so the UI never freezes.
     * Updates OutputPanel and OperationPanel on the FX thread when done.
     */
    private void dispatchSolve(String expression, int precisionFlag) {
        PrecisionMode mode = (precisionFlag == 0)
            ? PrecisionMode.SYMBOLIC
            : PrecisionMode.NUMERICAL;

        String opName    = opPanel.getCurrentOperation().label;
        String precName  = mode == PrecisionMode.SYMBOLIC ? "symbolic" : "numerical";

        // ── UI: loading state ─────────────────────────────────────────────────
        setStatus("Computing…", "neutral");
        spinner.setVisible(true);
        statusDot.getStyleClass().removeAll("dot-ok", "dot-error", "dot-warn");
        statusDot.getStyleClass().add("dot-neutral");

        // ── Translate expression + operation into the correct protocol string ──
        // Plain arithmetic stays as-is; derivative/integral/limit/etc. become
        // "calc:<op>|{json}" strings that CoreEngine routes to the C++ calculus module.
        String engineInput = CalcExpressionBuilder.build(
            expression, opPanel.getCurrentOperation());

        // ── Background task ───────────────────────────────────────────────────
        Task<String> task = new Task<>() {
            @Override protected String call() {
                return bridge.compute(engineInput, mode);
            }
        };

        task.setOnSucceeded(e -> {
            spinner.setVisible(false);
            String result = task.getValue();
            boolean isArith = opPanel.getCurrentOperation() == OperationPanel.Operation.ARITHMETIC;
            outputPanel.showResult(expression, result, null, opName, precName, isArith);
            // If result contains ODE solution points, offer to plot them
            if (result.contains("Points: [")) {
                int idx = result.indexOf("Points: [");
                graphPanel.plotPoints(result.substring(idx + 9), engineInput, null);
            }
            opPanel.addHistory(new OperationPanel.HistoryEntry(
                expression, result, opPanel.getCurrentOperation(), precName));
            // Push to server if logged in (fire-and-forget)
            auth.pushHistoryAsync(expression, result,
                opPanel.getCurrentOperation().label);
            setStatus("Done", "ok");
        });

        task.setOnFailed(e -> {
            spinner.setVisible(false);
            Throwable ex = task.getException();
            String msg = ex != null ? ex.getMessage() : "Unknown error";
            outputPanel.showError(expression, msg);
            setStatus("Error — " + msg, "error");
        });

        Thread t = new Thread(task);
        t.setDaemon(true);
        t.start();
    }

    // ── Status bar helpers ────────────────────────────────────────────────────

    private void setStatus(String message, String type) {
        Platform.runLater(() -> {
            statusText.setText(message);
            statusText.setAccessibleText("Status: " + message);
            statusDot.getStyleClass().removeAll("dot-ok", "dot-error", "dot-warn", "dot-neutral");
            switch (type) {
                case "ok"      -> statusDot.getStyleClass().add("dot-ok");
                case "error"   -> statusDot.getStyleClass().add("dot-error");
                case "warn"    -> statusDot.getStyleClass().add("dot-warn");
                default        -> statusDot.getStyleClass().add("dot-neutral");
            }
        });
    }

    // ── Auth ──────────────────────────────────────────────────────────────────

    private void doLogout() {
        auth.logout();
        // Transition back to login screen
        javafx.application.Platform.runLater(() -> {
            LoginScreen ls = new LoginScreen();
            javafx.scene.Scene scene = getScene();
            if (scene != null) {
                ls.setOnSuccess(u -> {
                    MainLayout main = new MainLayout(u);
                    scene.setRoot(main);
                });
                ls.setOnGuest(() -> {
                    MainLayout main = new MainLayout(null);
                    scene.setRoot(main);
                });
                scene.setRoot(ls);
            }
        });
    }

    // ── Helper ────────────────────────────────────────────────────────────────
    private static class Spacer extends Region {
        Spacer() { HBox.setHgrow(this, Priority.ALWAYS); }
    }
}
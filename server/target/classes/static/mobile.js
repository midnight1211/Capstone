'use strict';

// ── State ─────────────────────────────────────────────────────────────────────

const state = {
    token:     localStorage.getItem('me_token')    || null,
    username:  localStorage.getItem('me_username') || null,
    theme:     localStorage.getItem('me_theme')    || 'dark',
    fontSize:  localStorage.getItem('me_fontsize') || 'normal',
    operation: 'arithmetic',
    mode:      'symbolic',
    lastExpr:  null,
    lastResult: null,
};

const THEMES     = ['dark', 'light', 'system', 'high-contrast'];
const FONTSIZES  = ['normal', 'large', 'x-large'];

// ── Init ──────────────────────────────────────────────────────────────────────

document.addEventListener('DOMContentLoaded', () => {
    applyTheme(state.theme);
    applyFontSize(state.fontSize);
    updateAuthUI();
    checkServerStatus();

    // Ctrl+Enter / Enter (mobile keyboard "Go") solves
    document.getElementById('expr-input').addEventListener('keydown', e => {
        if (e.key === 'Enter' && (e.ctrlKey || e.metaKey)) {
            e.preventDefault();
            solve();
        }
    });

    // If signed in, load history
    if (state.token) loadHistory();
});

// ── Theme ─────────────────────────────────────────────────────────────────────

function cycleTheme() {
    const idx = THEMES.indexOf(state.theme);
    state.theme = THEMES[(idx + 1) % THEMES.length];
    localStorage.setItem('me_theme', state.theme);
    applyTheme(state.theme);

    // Sync to server if signed in
    if (state.token) savePreferences();
}

function applyTheme(preference) {
    let resolved = preference;
    if (preference === 'system') {
        resolved = window.matchMedia('(prefers-color-scheme: dark)').matches
            ? 'dark' : 'light';
    }
    document.body.setAttribute('data-theme', resolved);
}

// ── Font size ─────────────────────────────────────────────────────────────────

function cycleFontSize() {
    const idx = FONTSIZES.indexOf(state.fontSize);
    state.fontSize = FONTSIZES[(idx + 1) % FONTSIZES.length];
    localStorage.setItem('me_fontsize', state.fontSize);
    applyFontSize(state.fontSize);

    const btn = document.getElementById('btn-font');
    btn.textContent = { normal: 'A', large: 'A+', 'x-large': 'A++' }[state.fontSize];

    if (state.token) savePreferences();
}

function applyFontSize(size) {
    document.body.setAttribute('data-font-size', size);
}

// ── Operation selector ────────────────────────────────────────────────────────

function setOp(btn, op) {
    state.operation = op;
    document.querySelectorAll('.op-tab').forEach(b => {
        const active = b === btn;
        b.classList.toggle('active', active);
        b.setAttribute('aria-selected', String(active));
    });

    const placeholders = {
        arithmetic: 'sin(pi/2) + 5   or   sqrt(2^2 + 3^2)',
        integral:   '\\int_0^8 x^2\\,dx   or   integrate(x^2, 0, 8)',
        derivative: 'd/dx[x^3 + sin(x)]',
        limit:      '\\lim_{x \\to 0} \\frac{sin(x)}{x}',
        matrix:     'matrix:[[1,2],[3,4]]',
    };
    document.getElementById('expr-input').placeholder =
        placeholders[op] || 'Enter expression...';
}

// ── Precision mode ────────────────────────────────────────────────────────────

function setMode(mode) {
    state.mode = mode;
    document.getElementById('pill-sym').classList.toggle('active', mode === 'symbolic');
    document.getElementById('pill-sym').setAttribute('aria-pressed', String(mode === 'symbolic'));
    document.getElementById('pill-num').classList.toggle('active', mode === 'numerical');
    document.getElementById('pill-num').setAttribute('aria-pressed', String(mode === 'numerical'));
}

// ── Solve ─────────────────────────────────────────────────────────────────────

async function solve() {
    const expression = document.getElementById('expr-input').value.trim();
    if (!expression) return;

    const btn = document.getElementById('btn-solve');
    btn.disabled = true;
    btn.textContent = 'Working...';
    setStatus('Computing...', 'warn');
    hideResult();

    try {
        const body = {
            expression,
            precisionFlag: state.mode === 'symbolic' ? 0 : 1,
            operation: state.operation,
        };

        const res = await apiFetch('/api/compute', {
            method:  'POST',
            headers: { 'Content-Type': 'application/json' },
            body:    JSON.stringify(body),
        });

        if (res.ok) {
            state.lastExpr   = expression;
            state.lastResult = res.result;
            showResult(expression, res.result);
            setStatus('Ready', 'ok');

            // Refresh history panel if signed in
            if (state.token) loadHistory();
        } else {
            showError(res.error || 'Computation failed');
            setStatus('Error', 'error');
        }
    } catch (err) {
        showError('Could not reach server. Is it running?');
        setStatus('Server unreachable', 'error');
    } finally {
        btn.disabled = false;
        btn.innerHTML = 'Solve &#9654;';
    }
}

// ── Display result ────────────────────────────────────────────────────────────

function showResult(expression, result) {
    document.getElementById('result-empty').style.display   = 'none';
    document.getElementById('result-content').classList.remove('hidden');
    document.getElementById('result-error').classList.add('hidden');

    document.getElementById('result-input-expr').textContent = expression;

    // Split "sym  ~  num" if present
    if (result && result.includes('  ~  ')) {
        const [sym, num] = result.split('  ~  ');
        document.getElementById('result-value').textContent  = sym.trim();
        document.getElementById('result-approx').textContent = '\u2248 ' + num.trim();
    } else {
        document.getElementById('result-value').textContent  = result || '';
        document.getElementById('result-approx').textContent = '';
    }
    document.getElementById('copy-toast').textContent = '';
}

function hideResult() {
    document.getElementById('result-content').classList.add('hidden');
    document.getElementById('result-error').classList.add('hidden');
    document.getElementById('result-empty').style.display = '';
}

function showError(msg) {
    document.getElementById('result-empty').style.display   = 'none';
    document.getElementById('result-content').classList.add('hidden');
    document.getElementById('result-error').classList.remove('hidden');
    document.getElementById('error-message').textContent = msg;
}

function clearAll() {
    document.getElementById('expr-input').value = '';
    hideResult();
    document.getElementById('expr-input').focus();
}

// ── Export ────────────────────────────────────────────────────────────────────

async function exportResult(format) {
    if (!state.lastExpr || !state.lastResult) return;

    try {
        const res = await apiFetch('/api/export', {
            method:  'POST',
            headers: { 'Content-Type': 'application/json' },
            body:    JSON.stringify({
                expression: state.lastExpr,
                result:     state.lastResult,
                format,
            }),
        });

        if (res.content) {
            await navigator.clipboard.writeText(res.content);
            showToast('Copied as ' + format + '!');
        }
    } catch (err) {
        showToast('Copy failed: ' + err.message);
    }
}

function showToast(msg) {
    const toast = document.getElementById('copy-toast');
    toast.textContent = msg;
    setTimeout(() => { toast.textContent = ''; }, 2500);
}

// ── Auth panel ────────────────────────────────────────────────────────────────

function toggleAuthPanel() {
    const panel = document.getElementById('auth-panel');
    const open  = panel.classList.toggle('open');
    panel.setAttribute('aria-hidden', String(!open));
}

function setAuthTab(tab) {
    document.getElementById('tab-login').classList.toggle('active', tab === 'login');
    document.getElementById('tab-register').classList.toggle('active', tab === 'register');
    document.getElementById('tab-login').setAttribute('aria-selected', String(tab === 'login'));
    document.getElementById('tab-register').setAttribute('aria-selected', String(tab === 'register'));
    document.getElementById('form-login').classList.toggle('hidden', tab !== 'login');
    document.getElementById('form-register').classList.toggle('hidden', tab !== 'register');
}

// ── Auth actions ──────────────────────────────────────────────────────────────

async function doLogin() {
    const username = document.getElementById('login-user').value.trim();
    const password = document.getElementById('login-pass').value;
    const errEl    = document.getElementById('login-error');
    errEl.textContent = '';

    if (!username || !password) {
        errEl.textContent = 'Username and password are required.';
        return;
    }

    try {
        const res = await apiFetch('/api/auth/login', {
            method:  'POST',
            headers: { 'Content-Type': 'application/json' },
            body:    JSON.stringify({ username, password }),
        });

        if (res.token) {
            storeAuth(res);
            loadHistory();
            document.getElementById('auth-panel').classList.remove('open');
        } else {
            errEl.textContent = res.error || 'Login failed.';
        }
    } catch (err) {
        errEl.textContent = 'Server error: ' + err.message;
    }
}

async function doRegister() {
    const username = document.getElementById('reg-user').value.trim();
    const email    = document.getElementById('reg-email').value.trim();
    const password = document.getElementById('reg-pass').value;
    const errEl    = document.getElementById('register-error');
    errEl.textContent = '';

    if (!username || !email || !password) {
        errEl.textContent = 'All fields are required.';
        return;
    }

    try {
        const res = await apiFetch('/api/auth/register', {
            method:  'POST',
            headers: { 'Content-Type': 'application/json' },
            body:    JSON.stringify({ username, email, password }),
        });

        if (res.token) {
            storeAuth(res);
            loadHistory();
            document.getElementById('auth-panel').classList.remove('open');
        } else {
            errEl.textContent = res.error || 'Registration failed.';
        }
    } catch (err) {
        errEl.textContent = 'Server error: ' + err.message;
    }
}

function doLogout() {
    localStorage.removeItem('me_token');
    localStorage.removeItem('me_username');
    state.token    = null;
    state.username = null;
    updateAuthUI();
    document.getElementById('history-section').classList.add('hidden');
    document.getElementById('auth-panel').classList.remove('open');
}

function storeAuth(res) {
    state.token    = res.token;
    state.username = res.username;
    localStorage.setItem('me_token',    res.token);
    localStorage.setItem('me_username', res.username);

    // Apply server-side preferences
    if (res.theme)    { state.theme    = res.theme;    applyTheme(res.theme);    }
    if (res.fontSize) { state.fontSize = res.fontSize; applyFontSize(res.fontSize); }
    updateAuthUI();
}

function updateAuthUI() {
    const authBtn = document.getElementById('btn-auth');
    if (state.token) {
        authBtn.classList.add('signed-in');
        authBtn.setAttribute('aria-label', 'Account: ' + state.username);
        document.getElementById('auth-username').textContent = state.username;
        document.getElementById('form-loggedin').classList.remove('hidden');
        document.getElementById('form-login').classList.add('hidden');
        document.getElementById('form-register').classList.add('hidden');
    } else {
        authBtn.classList.remove('signed-in');
        authBtn.setAttribute('aria-label', 'Sign in or create account');
        document.getElementById('form-loggedin').classList.add('hidden');
        document.getElementById('form-login').classList.remove('hidden');
    }
}

// ── History ───────────────────────────────────────────────────────────────────

async function loadHistory() {
    if (!state.token) return;

    try {
        const entries = await apiFetch('/api/history?limit=30');
        const section = document.getElementById('history-section');
        const list    = document.getElementById('history-list');

        if (!Array.isArray(entries) || entries.length === 0) {
            section.classList.add('hidden');
            return;
        }

        section.classList.remove('hidden');
        list.innerHTML = entries.map(e => `
            <li class="history-item"
                role="button" tabindex="0"
                aria-label="Load: ${escHtml(e.expression)}"
                onclick="loadHistoryItem('${escAttr(e.expression)}')"
                onkeydown="if(event.key==='Enter')loadHistoryItem('${escAttr(e.expression)}')">
                <div class="history-expr">${escHtml(e.expression)}</div>
                <div class="history-meta">${escHtml(e.operation)} &middot; ${escHtml(e.precisionMode)}</div>
            </li>
        `).join('');
    } catch (err) {
        console.warn('[History] Load failed:', err);
    }
}

function loadHistoryItem(expression) {
    document.getElementById('expr-input').value = expression;
    document.getElementById('expr-input').focus();
    window.scrollTo({ top: 0, behavior: 'smooth' });
}

async function clearHistory() {
    if (!state.token) return;
    try {
        await apiFetch('/api/history', { method: 'DELETE' });
        document.getElementById('history-list').innerHTML = '';
        document.getElementById('history-section').classList.add('hidden');
    } catch (err) {
        console.warn('[History] Clear failed:', err);
    }
}

// ── Preferences sync ──────────────────────────────────────────────────────────

async function savePreferences() {
    if (!state.token) return;
    try {
        await apiFetch('/api/preferences', {
            method:  'PUT',
            headers: { 'Content-Type': 'application/json' },
            body:    JSON.stringify({ theme: state.theme, fontSize: state.fontSize }),
        });
    } catch (err) {
        console.warn('[Prefs] Save failed:', err);
    }
}

// ── Server status ─────────────────────────────────────────────────────────────

async function checkServerStatus() {
    try {
        const res = await fetch('/api/engine/status');
        if (res.ok) {
            const data = await res.json();
            setStatus(data.nativeLoaded ? 'C++ engine ready' : 'Stub mode (no native lib)', 'ok');
        } else {
            setStatus('Server error', 'error');
        }
    } catch {
        setStatus('Server unreachable', 'error');
    }
}

// ── Status bar ────────────────────────────────────────────────────────────────

function setStatus(msg, type) {
    document.getElementById('status-text').textContent = msg;
    const dot = document.getElementById('status-dot');
    dot.className = 'status-dot' + (type === 'ok' ? '' : ' ' + type);
}

// ── API fetch helper ──────────────────────────────────────────────────────────

async function apiFetch(url, options = {}) {
    if (state.token) {
        options.headers = options.headers || {};
        options.headers['Authorization'] = 'Bearer ' + state.token;
    }
    const res  = await fetch(url, options);
    const data = await res.json();
    return data;
}

// ── Utilities ─────────────────────────────────────────────────────────────────

function escHtml(str) {
    return String(str)
        .replace(/&/g,'&amp;').replace(/</g,'&lt;')
        .replace(/>/g,'&gt;').replace(/"/g,'&quot;');
}

function escAttr(str) {
    return String(str).replace(/'/g, "\\'").replace(/"/g, '\\"');
}
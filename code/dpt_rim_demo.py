#!/usr/bin/env python3
"""
DPT Residual Integrity Monitor (RIM) demo
- Null-calibrated thresholds (Pi, CB)
- Tail veto to suppress outlier-dominated artifacts
- Multi-GNSS ISB trap + Common-Mode Oscillator (shared receiver clock)
- Autonomous triage: clock vs ISB/frame leak vs mixed vs nominal

Run:
  python code/dpt_rim_demo.py
"""

import numpy as np

EPS = 1e-12


# -----------------------------
# DPT core
# -----------------------------
def second_diff(x: np.ndarray) -> np.ndarray:
    return x[2:] - 2 * x[1:-1] + x[:-2]


def penalized_smoother_d2(w: np.ndarray, lam: float) -> np.ndarray:
    """
    Solve: u = argmin ||w-u||^2 + lam ||D^2 u||^2
    via linear system: (I + lam * D2^T D2) u = w
    """
    n = len(w)
    if n < 5:
        return w.copy()

    A = np.eye(n)
    main = np.r_[1, 5, np.full(n - 4, 6), 5, 1]
    off1 = np.r_[-2, np.full(n - 3, -4), -2]
    off2 = np.full(n - 2, 1)

    A += lam * (
        np.diag(main)
        + np.diag(off1, 1) + np.diag(off1, -1)
        + np.diag(off2, 2) + np.diag(off2, -2)
    )
    return np.linalg.solve(A, w)


def spectral_entropy(r: np.ndarray) -> float:
    """Normalized spectral entropy in [0,1]."""
    n = len(r)
    if n < 8:
        return 1.0
    psd = np.abs(np.fft.rfft(r * np.hanning(n))) ** 2 + EPS
    p = psd / (psd.sum() + EPS)
    H = -np.sum(p * np.log(p + EPS))
    Hmax = np.log(len(p) + EPS)
    return float(H / (Hmax + EPS))


def tail_fraction(r: np.ndarray, T: float = 3.0) -> float:
    return float(np.mean(np.abs(r) > T))


def autocorr_peak(r: np.ndarray, lag_min: int = 3, lag_max: int = 60) -> float:
    n = len(r)
    lag_max = min(lag_max, max(lag_min, n // 4))
    if lag_max <= lag_min:
        return 0.0
    x0 = r - np.mean(r)
    denom = float(np.dot(x0, x0) + EPS)
    best = 0.0
    for l in range(lag_min, lag_max + 1):
        val = float(np.dot(x0[l:], x0[:-l]) / denom)
        best = max(best, abs(val))
    return best


def compute_metrics(w: np.ndarray, lam: float, tail_T: float = 3.0):
    """
    Returns:
      Pi, CB, extras dict
    """
    u = penalized_smoother_d2(w, lam)
    r = w - u

    E = float(np.mean(r**2))
    H = spectral_entropy(r)
    tail = tail_fraction(r, tail_T)
    A = autocorr_peak(r)

    # Complexity Budget: pure curvature effort / residual energy
    kappa = float(np.sum(second_diff(u)**2)) if len(u) >= 3 else 0.0
    res_energy = float(np.sum(r**2))
    CB = kappa / (res_energy + EPS)

    # Persistence score (tail veto included)
    Pi = E * (1.0 - H) * (1.0 - tail) * A

    extras = {"E": E, "Hspec": H, "tail": tail, "A": A, "kappa": kappa, "res_energy": res_energy}
    return Pi, CB, extras


# -----------------------------
# Null calibration
# -----------------------------
def calibrate_null(n: int, lambdas, n_trials: int = 300, alpha: float = 1e-3, seed: int = 11):
    """
    Empirical thresholds for Pi and CB under w ~ N(0,1).
    """
    rng = np.random.default_rng(seed)
    pi_samp = {lam: [] for lam in lambdas}
    cb_samp = {lam: [] for lam in lambdas}

    for _ in range(n_trials):
        w = rng.normal(0.0, 1.0, n)
        for lam in lambdas:
            Pi, CB, _ = compute_metrics(w, lam)
            pi_samp[lam].append(Pi)
            cb_samp[lam].append(CB)

    thr_pi = {lam: float(np.quantile(pi_samp[lam], 1 - alpha)) for lam in lambdas}
    thr_cb = {lam: float(np.quantile(cb_samp[lam], 1 - alpha)) for lam in lambdas}
    return thr_pi, thr_cb


def stable_hit_across_lambdas(Pis, CBs, thr_pi, thr_cb, lambdas, need: int = 2) -> tuple[bool, int]:
    """
    Stable alarm if >= need lambdas exceed either Pi or CB threshold.
    """
    hits = 0
    for lam, Pi, CB in zip(lambdas, Pis, CBs):
        hits += int((Pi > thr_pi[lam]) or (CB > thr_cb[lam]))
    return hits >= need, hits


# -----------------------------
# Signal models: CMO + ISB mismatch
# -----------------------------
def make_common_mode_oscillator(t: np.ndarray, amp: float = 0.15, period: float = 400.0) -> np.ndarray:
    """Shared receiver clock wobble (simple sinusoid)."""
    return amp * np.sin(2 * np.pi * t / period)


def make_isb_mismatch(t: np.ndarray, n: int, ramp_amp: float = 0.12, sin_amp: float = 0.08, sin_period: float = 200.0) -> np.ndarray:
    """Galileo-only mismatch: drift + periodic component."""
    ramp = ramp_amp * (t / max(1, n))
    sinus = sin_amp * np.sin(2 * np.pi * t / sin_period)
    return ramp + sinus


# -----------------------------
# Triage logic
# -----------------------------
def perform_triage(gps_hit: bool, gal_hit: bool, diff_hit: bool) -> str:
    if diff_hit and not (gps_hit and gal_hit):
        return "DIAGNOSIS: ISB TRAP / RELATIVISTIC FRAME LEAK (structure isolated in difference stream)."
    if gps_hit and gal_hit and not diff_hit:
        return "DIAGNOSIS: RECEIVER CLOCK INSTABILITY (shared persistence cancels in difference)."
    if gps_hit and gal_hit and diff_hit:
        return "DIAGNOSIS: MIXED FAULT (receiver clock instability + ISB/frame leak)."
    if gps_hit or gal_hit:
        return "DIAGNOSIS: CONSTELLATION-SPECIFIC ANOMALY (one stream persistent)."
    return "STATUS: NOMINAL (no persistence beyond null-calibrated thresholds)."


# -----------------------------
# Demo run
# -----------------------------
def run_demo():
    # Config
    n = 1000
    t = np.arange(n)
    lambdas = [25.0, 100.0]

    # Reproducibility
    np.random.seed(7)

    # 1) Common-mode receiver clock drift
    common_drift = make_common_mode_oscillator(t, amp=0.15, period=400.0)

    # 2) ISB mismatch (Galileo-only)
    isb_mismatch = make_isb_mismatch(t, n=n, ramp_amp=0.12, sin_amp=0.08, sin_period=200.0)

    # 3) Refined noise injection: shared + independent
    shared_noise = np.random.normal(0.0, 0.7, n)   # shared front-end noise
    gps_ind = np.random.normal(0.0, 0.3, n)        # independent per constellation
    gal_ind = np.random.normal(0.0, 0.3, n)

    w_gps = shared_noise + gps_ind + common_drift
    w_gal = shared_noise + gal_ind + common_drift + isb_mismatch
    w_diff = w_gps - w_gal  # shared noise + common drift cancel; mismatch remains

    # 4) Null calibration (fast)
    thr_pi, thr_cb = calibrate_null(n=n, lambdas=lambdas, n_trials=250, alpha=1e-3, seed=11)

    print("\n--- DPT-RIM Multi-GNSS Triage Demo (Calibrated) ---")
    print(f"n={n}, lambdas={lambdas}, alpha=1e-3, stable_need=2\n")

    results = {}
    for name, stream in [("GPS", w_gps), ("GAL", w_gal), ("DIFF", w_diff)]:
        Pis, CBs, extras = [], [], []
        for lam in lambdas:
            Pi, CB

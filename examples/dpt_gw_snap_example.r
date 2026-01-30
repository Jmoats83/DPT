# DPT Example: Detecting a Persistent Snap in a Smooth Chirp Signal

This file provides a concrete, runnable example of how Distinction-Persistence Theory (DPT) is used in practice.

We construct a smooth, predictable signal, inject a brief but structured disturbance ("snap"), subtract the best smooth model, and then compute a DPT score that isolates the persistent distinction as a localized peak.

This example is synthetic but falsifiable. Anyone can run it, modify it, or attempt to break it.

---

## Conceptual Setup (DPT framing)

We define:

Observed signal:
x(t)

Best smooth model:
m(t)

Residual:
r(t) = x(t) - m(t)

The DPT coherence-retention score is defined as:

Lambda_hat(t) =
Q(t) * (1 - H_spec(t)) * (1 - Tail_frac(t)) * (1 - E_str(t))

where Lambda_hat(t) is bounded in [0, 1].

A persistent distinction appears as a localized peak in Lambda_hat(t).

---

## What This Example Tests

1. A smooth chirp signal is generated (frequency and amplitude increase smoothly).
2. A brief, high-frequency snap is injected near the peak.
3. Noise is added.
4. A smooth template is reconstructed from the noisy signal.
5. The residual is evaluated using the DPT score.

The injected snap is designed to be:
- strong (high Q)
- non-noise-like (low spectral entropy)
- temporally localized (low late-energy fraction)
- locally coherent (low structure error)

---

## How to Run

Copy the Python code below into a local file named:

dpt_gw_snap_example.py

Then run:

python3 dpt_gw_snap_example.py

Dependencies:
- numpy
- matplotlib
- scipy (optional; only used for smoother template fitting)

---

## Runnable Python Code

```python
import numpy as np
import matplotlib.pyplot as plt

try:
    from scipy.signal import savgol_filter
except Exception:
    savgol_filter = None


def clip01(x):
    return np.clip(x, 0.0, 1.0)


def spectral_entropy(window, eps=1e-12):
    w = window - np.mean(window)
    spec = np.abs(np.fft.rfft(w)) ** 2 + eps
    p = spec / np.sum(spec)
    H = -np.sum(p * np.log(p))
    return clip01(H / np.log(len(p)))


def late_energy_fraction(resid, idx, half_win):
    start = max(0, idx - half_win)
    end = min(len(resid), idx + half_win + 1)
    w = resid[start:end]
    center = idx - start
    right = w[center + 1 :]
    total = np.sum(w ** 2) + 1e-12
    late = np.sum(right ** 2)
    return clip01(late / total)


def structure_error(resid, idx, seg_len=128, max_lag=24, eps=1e-12):
    half = seg_len // 2
    start = max(0, idx - half)
    end = min(len(resid), idx + half)
    seg = resid[start:end].copy()

    if len(seg) < 32:
        return 1.0

    seg = seg - np.mean(seg)
    best = 0.0

    for lag in range(1, max_lag + 1):
        if lag >= len(seg):
            break
        a = seg[:-lag]
        b = seg[lag:]
        denom = (np.sqrt(np.sum(a ** 2)) + eps) * (np.sqrt(np.sum(b ** 2)) + eps)
        corr = np.dot(a, b) / denom
        best = max(best, corr)

    return clip01(1.0 - best)


def dpt_score(x, m, fs, win_seconds=0.22):
    r = x - m
    half_win = int(max(8, (win_seconds * fs) // 2))
    ref = np.median(np.abs(r - np.median(r))) + 1e-12

    Q = np.zeros_like(r)
    H = np.zeros_like(r)
    Tail = np.zeros_like(r)
    Estr = np.zeros_like(r)

    for i in range(len(r)):
        Q[i] = clip01(abs(r[i]) / (6.0 * ref))

        a = max(0, i - half_win)
        b = min(len(r), i + half_win + 1)
        w = r[a:b]

        H[i] = spectral_entropy(w)
        Tail[i] = late_energy_fraction(r, i, half_win)
        Estr[i] = structure_error(r, i)

    Lambda = Q * (1.0 - H) * (1.0 - Tail) * (1.0 - Estr)
    return r, clip01(Lambda)


def make_chirp(t, f0=20.0, f1=220.0):
    T = t[-1] - t[0]
    k = (f1 - f0) / T
    phase = 2.0 * np.pi * (f0 * t + 0.5 * k * t ** 2)
    amp = (t / T) ** 1.5
    return amp * np.sin(phase)


def inject_snap(x, t, t0, strength=0.35, width=0.012, freq=420.0):
    env = np.exp(-0.5 * ((t - t0) / width) ** 2)
    snap = strength * env * np.sin(2 * np.pi * freq * (t - t0))
    return x + snap


def smooth_template(x, fs):
    if savgol_filter is not None:
        w = int(max(31, (0.12 * fs) // 2 * 2 + 1))
        return savgol_filter(x, w, 3)
    else:
        w = int(max(31, 0.12 * fs))
        kernel = np.ones(w) / w
        return np.convolve(x, kernel, mode="same")


def main():
    fs = 2048
    duration = 2.0
    t = np.arange(0, duration, 1.0 / fs)

    base = make_chirp(t)
    t0 = 1.55
    x_clean = inject_snap(base, t, t0)

    rng = np.random.default_rng(7)
    noise = 0.10 * rng.standard_normal(len(t))
    x = x_clean + noise

    m = smooth_template(x, fs)
    r, Lambda = dpt_score(x, m, fs)

    plt.figure(figsize=(10, 6))
    plt.subplot(3, 1, 1)
    plt.plot(t, x)
    plt.plot(t, m)
    plt.title("Observed signal and smooth model")

    plt.subplot(3, 1, 2)
    plt.plot(t, r)
    plt.axvline(t0, linestyle="--")
    plt.title("Residual")

    plt.subplot(3, 1, 3)
    plt.plot(t, Lambda)
    plt.axvline(t0, linestyle="--")
    plt.title("DPT score (Lambda_hat)")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()

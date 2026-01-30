# DPT Example: Rejecting Tail-Smeared Energy as a Persistent Distinction

This example demonstrates why Distinction-Persistence Theory (DPT) explicitly penalizes
late-arriving or smeared energy.

Two signals can have similar total energy, but only one represents a coherent,
time-localized distinction. The other spreads its energy into the tail and should
not be classified as a real event.

This example isolates the role of the Tail_frac term in the DPT score.

---

## Conceptual Setup (DPT framing)

We define:

Observed signal:
x(t)

Best smooth model:
m(t)

Residual:
r(t) = x(t) - m(t)

The DPT score is:

Lambda_hat(t) =
Q(t) * (1 - H_spec(t)) * (1 - Tail_frac(t)) * (1 - E_str(t))

In this example, Q and spectral properties are intentionally similar.
The only meaningful difference is how much residual energy arrives late.

---

## What This Example Tests

We compare two residual events:

Event A: Compact burst
- Energy arrives quickly
- Minimal late-time leakage
- Low Tail_frac

Event B: Tail-smeared burst
- Same total energy
- Energy leaks slowly into the future
- High Tail_frac

DPT should strongly favor Event A and suppress Event B.

---

## How to Run

Copy the Python code below into a local file named:

dpt_tail_energy_example.py

Then run:

python3 dpt_tail_energy_example.py

Dependencies:
- numpy
- matplotlib

---

## Runnable Python Code

```python
import numpy as np
import matplotlib.pyplot as plt


def clip01(x):
    return np.clip(x, 0.0, 1.0)


def late_energy_fraction(resid, idx, half_win):
    start = max(0, idx - half_win)
    end = min(len(resid), idx + half_win + 1)
    w = resid[start:end]
    center = idx - start
    right = w[center + 1 :]
    total = np.sum(w ** 2) + 1e-12
    late = np.sum(right ** 2)
    return clip01(late / total)


def generate_compact_event(t, t0, width=0.02):
    return np.exp(-0.5 * ((t - t0) / width) ** 2)


def generate_tail_event(t, t0, decay=0.15):
    env = np.zeros_like(t)
    idx = t >= t0
    env[idx] = np.exp(-(t[idx] - t0) / decay)
    return env


def main():
    fs = 2000
    duration = 1.0
    t = np.arange(0, duration, 1.0 / fs)

    t0 = 0.4
    noise = 0.02 * np.random.randn(len(t))

    # Event A: compact
    r_compact = generate_compact_event(t, t0)
    r_compact = r_compact / np.linalg.norm(r_compact)

    # Event B: tail-smeared
    r_tail = generate_tail_event(t, t0)
    r_tail = r_tail / np.linalg.norm(r_tail)

    r_compact += noise
    r_tail += noise

    half_win = int(0.15 * fs)

    Tail_A = np.array([
        late_energy_fraction(r_compact, i, half_win)
        for i in range(len(t))
    ])

    Tail_B = np.array([
        late_energy_fraction(r_tail, i, half_win)
        for i in range(len(t))
    ])

    plt.figure(figsize=(10, 7))

    plt.subplot(3, 1, 1)
    plt.plot(t, r_compact)
    plt.axvline(t0, linestyle="--")
    plt.title("Residual A: Compact Event")

    plt.subplot(3, 1, 2)
    plt.plot(t, r_tail)
    plt.axvline(t0, linestyle="--")
    plt.title("Residual B: Tail-Smeared Event")

    plt.subplot(3, 1, 3)
    plt.plot(t, Tail_A, label="Compact event")
    plt.plot(t, Tail_B, label="Tail-smeared event")
    plt.axvline(t0, linestyle="--")
    plt.title("Tail_frac comparison")
    plt.legend()

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()

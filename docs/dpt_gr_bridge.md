# DPT and General Relativity: Residual Analysis as a Shared Interface

This file explains how Distinction-Persistence Theory (DPT) relates to General Relativity (GR)
in a practical and testable way.

DPT is not presented as a replacement for GR. The clean framing is that GR supplies the best
smooth physical model, while DPT evaluates what remains after that model is subtracted.

In short:
- GR explains the smooth, expected structure.
- DPT evaluates the residual for persistent, ordered structure.

This mirrors how modern data analysis is already performed in practice.

---

## Shared Setup

Observed data:
x(t)

GR template family:
m(t; theta)

Best-fit parameters:
theta_hat

Residual:
r(t) = x(t) - m(t; theta_hat)

DPT score:
Lambda_hat(t) =
Q(t) * (1 - H_spec(t)) * (1 - Tail_frac(t)) * (1 - E_str(t))

Lambda_hat(t) is bounded between 0 and 1.

Interpretation:
- GR answers: what smooth physical waveform best fits the data?
- DPT answers: after subtracting that waveform, does structured residue persist?

---

## Cooperative Use of DPT and GR

In standard gravitational-wave pipelines, GR waveforms are fit to the data using matched
filtering or Bayesian inference. Analysts then examine residuals to assess data quality and
model adequacy.

DPT can be used as a disciplined residual diagnostic:

1. Fit a GR waveform to obtain m(t; theta_hat).
2. Compute residual r(t).
3. Compute Lambda_hat(t) over time.
4. Flag time intervals where Lambda_hat(t) is elevated.

Elevated DPT scores indicate residual structure that is:
- stronger than expected
- non-noise-like
- temporally localized
- internally coherent

This does not imply new physics by itself. It identifies where further scrutiny is warranted.

---

## Where DPT Can Challenge GR (Narrow and Falsifiable)

The head-to-head claim is not that GR is incorrect. The claim is narrower and testable:

Claim:
After subtracting the best-fit GR waveform, the residual sometimes contains localized,
persistent structure more often than expected under the noise model and known systematics.

If such structure exists, possible explanations include:
1. waveform modeling omissions (higher modes, precession, eccentricity, memory)
2. calibration or instrumental artifacts
3. environmental coupling
4. genuinely new physics (last and most conservative hypothesis)

DPT does not identify the cause. It identifies the anomaly.

---

## A Concrete Residual Test

Goal:
Test whether persistent GR residuals occur more frequently than expected.

Inputs:
- x(t): observed strain data
- m(t; theta_hat): best-fit GR waveform
- control data: off-source segments or noise-only injections

Procedure:
1. Fit GR waveform to obtain theta_hat.
2. Compute residual r(t).
3. Compute Lambda_hat(t) using a sliding window.
4. Identify peaks where Lambda_hat(t) exceeds a threshold.
5. Compare peak statistics against controls.

Metrics to report:
- peak heights
- timing relative to merger or peak amplitude
- false-positive rate under controls

Pass condition:
Residual peaks cluster near a consistent physical phase and exceed control expectations.

Fail condition:
Residual peaks match noise statistics or correlate with known artifacts.

Both outcomes are informative.

---

## Practical Definitions (Implementation-Oriented)

Q(t):
Normalized residual amplitude relative to a robust scale estimate.

H_spec(t):
Spectral entropy of a short residual window.
Low values indicate ordered frequency content.

Tail_frac(t):
Fraction of residual energy arriving after the center of a symmetric window.
This penalizes smeared, late-arriving energy.

E_str(t):
Structure error proxy based on short-lag autocorrelation.
This penalizes fragmented or incoherent residue.

Lambda_hat(t):
Product of the above terms, emphasizing strong, ordered, localized residue.

---

## Recommended Public Framing

To remain credible and collaborative, use language such as:
- DPT complements GR by analyzing residual structure.
- High DPT scores flag regions for deeper modeling or data-quality checks.
- New-physics interpretations are treated as a last hypothesis.

Avoid claims that DPT disproves GR. The correct framing is residual analysis.

---

## Suggested Repository Placement

Place this file at:
docs/dpt_gr_bridge.md

Optionally reference it in the README with a single line pointing readers to the GR connection.

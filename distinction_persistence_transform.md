\documentclass[11pt]{article}

\usepackage{amsmath,amssymb}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{geometry}
\usepackage{hyperref}

\geometry{margin=1in}

\title{\textbf{The Distinction--Persistence Transform:}\\
Signal Validation via Manifold Perturbation}

\author{
Anonymous Author(s)\\
\small{Affiliation withheld for review}
}

\date{}

\begin{document}
\maketitle

\begin{abstract}
We present the Distinction--Persistence Transform (DPT), a non-parametric framework for identifying periodic anomalies in noisy time-series data. Unlike matched-filtering techniques that maximize cross-correlation with fixed templates, DPT evaluates the information-theoretic stability of a signal under manifold perturbation. By requiring signals to simultaneously satisfy low-entropy and high-persistence criteria across varying detrending scales, DPT effectively demotes common astrophysical and instrumental false positives (e.g., stellar rotation, radio frequency interference, and red noise). We demonstrate that DPT reliably recovers low signal-to-noise transits where traditional power-threshold methods fail, providing a robust path for technosignature and exoplanet validation in high-cadence surveys.
\end{abstract}

\section{Introduction}

The current generation of transit, radial velocity, and SETI pipelines has achieved unprecedented sensitivity, largely through the optimization of power-threshold methods such as Box Least Squares and matched filtering. However, as surveys push into lower signal-to-noise ratio (SNR) regimes, the bottleneck is no longer detection, but belief calibration. Systematic artifacts, stochastic red noise, and Earth-based radio frequency interference (RFI) frequently manifest with sufficient periodic power to clear traditional significance thresholds, leading to a high rate of false positives that require expensive follow-up resources.

In this paper, we introduce the Distinction--Persistence Transform (DPT), a non-parametric diagnostic lens designed to be applied post-candidate generation. Rather than competing with discovery engines, DPT serves as a jurisdictional filter. We shift the fundamental question of signal validation: we ask not whether a signal is strong, but whether it is difficult to destroy. By subjecting a candidate to coordinate stress, DPT reframes detection as a measure of survival under abuse rather than optimization under fit.

\section{Conceptual Framework: Distinction vs. Persistence}

The DPT framework is built upon two information-theoretic pillars: distinction, defined as deviation from a smooth, entropic null model, and persistence, defined as the stability of that deviation across a perturbed parameter manifold.

A candidate signal $\mathcal{S}$ is evaluated through a multiplicative ``AND-gate'' logic, where the final DPT score $\hat{\Lambda}$ is the product of independent survival metrics. This structure ensures that any single failure mode---such as high morphological entropy or instability under detrending---collapses the total belief in the signal.

\begin{equation}
\hat{\Lambda}(t) = Q(t)\,(1-H(t))\,(1-\mathcal{T}(t))\,(1-\mathcal{E}(t)),
\end{equation}

where $Q$ represents raw distinction, $H$ the spectral entropy, $\mathcal{T}$ the tail or smear fraction, and $\mathcal{E}$ the fragmentation or echo coefficient. In this regime, a genuine signal must maintain structural integrity across all dimensions simultaneously.

\section{The Stress Test Protocol: Manifold Perturbation}

To distinguish physically persistent structures from transient or instrumental artifacts, we define a stability manifold $\mathcal{M}$ around the candidate parameters (period $P$, phase $\phi$, and duration $\tau$).

\subsection{Parameter Jitter (The ``Jiggle Test'')}

We subject the folding engine to stochastic Gaussian perturbations of the period, $P \pm \Delta P$. Brittle signals rely on narrow, accidental alignments and rapidly decohere. Genuine astrophysical transits exhibit a broad basin of attraction, maintaining elevated DPT scores under coordinate stress.

\subsection{Detrending Sensitivity (The ``Morphology Test'')}

Systematics are frequently induced by the detrending process itself. We therefore recompute the DPT score across a range of filter kernels (0.5--4.0 days).

\textit{Collapse Rule:} If the DPT score fluctuates by more than 50\% across detrending scales, the signal is flagged as instrumental. Genuine physical distinctions must remain invariant to baseline normalization.

\subsection{Scrambled Phase Control}

As a safeguard against pareidolia, the original time series is phase-scrambled. DPT must return a null result under this control, confirming sensitivity to time-ordered structural recurrence rather than global spectral shape.

\section{Validation on Exoplanet Transit Data}

\subsection{Validation Strategy}

DPT is evaluated strictly in a post-detection regime. The goal is not rediscovery, but differential reward of physical persistence and penalization of known false-positive classes.

\subsection{Case Study: Kepler-10 b}

Kepler-10 b is a low-depth ($\sim$150 ppm), short-period (0.837 d) confirmed exoplanet embedded in significant stellar and instrumental variability. Individual transit events are not visually compelling after detrending.

Applying DPT reveals localized persistence spikes coincident with transit epochs. Folding on trial periods yields a dominant persistence peak at the known orbital period, while harmonic aliases present in Lomb--Scargle periodograms are strongly suppressed.

This demonstrates that DPT amplifies repeatability rather than depth.

\subsection{False-Positive Controls}

Signals dominated by stellar rotation, red noise, or instrumental systematics collapse under manifold perturbation. Elevated entropy, phase smear, and detrending instability systematically demote such candidates without heuristic intervention.

\subsection{Sparse and Long-Period Regime}

In low-multiplicity regimes (2--3 observed transits), DPT promotes candidates that maintain structural recurrence across perturbations while demoting accidental alignments.

\section{Extension to Technosignature Vetting}

The mathematical structure of DPT is domain-agnostic. In time--frequency surveys, candidate signals are evaluated across alignment parameters such as Doppler drift rate, polarization state, and observation frame.

Terrestrial RFI is typically brittle under these perturbations, while astrophysical emitters are constrained by known physical processes. DPT does not classify intent; it identifies non-entropic persistence inconsistent with known interference classes.

\section{Discussion and Conclusion}

\subsection{Epistemic Boundaries}

The primary output of DPT is not a detection claim, but a refusal to discard. By shifting emphasis from power optimization to manifold stability, DPT provides a principled mechanism for ranking candidates in the low-SNR ``grey zone.''

\subsection{Implementation in Future Surveys}

The non-parametric nature of DPT makes it well suited to forthcoming high-volume surveys such as the \textit{Vera C. Rubin Observatory (LSST)} and the \textit{Next-Generation Transit Survey (NGTS)}, where prioritization of follow-up resources on assets like the \textit{James Webb Space Telescope (JWST)} will be essential.

\subsection{Conclusion}

In domains where false positives are cheap and belief is expensive, survival under stress is the appropriate currency. The Distinction--Persistence Transform provides a disciplined methodology for evaluating the structural integrity of anomalies across astrophysical and artificial domains, advancing signal validation beyond optimization toward principled skepticism.

\section*{Figures}

\textbf{Figure 1: The Stability Manifold $\mathcal{M}$.} Conceptual visualization of DPT score distribution across perturbed parameter space.

\textbf{Figure 2: Persistence vs. Power.} Lomb--Scargle periodogram compared with DPT persistence map for identical low-SNR data.

\textbf{Figure 3: Technosignature Diagnostic Ranking.} DPT response to synthetic narrowband signal versus terrestrial RFI.

\end{document}

"""
Project 8: Clinical Trial Analyzer
Module 4: Theory vs. Experiment — Power Analysis & Monte Carlo Validation

Compares exact theoretical statistical power (non-central t-distribution)
against empirical power estimated from 2,000 simulated two-arm trials per
sample size across multiple effect sizes (d = 0.3, 0.5, 0.8), plus:
  - Empirical Type I error verification under the null (d = 0.0, alpha = 0.05)
  - Assumption-violation experiment: heavy-tailed / skewed outcomes showing
    where empirical t-test power drops below normal-theory predictions.
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

SEED = 42


def theoretical_power_two_sample_t(n_per_arm: int, d: float, alpha: float = 0.05) -> float:
    """Exact two-sided power for independent two-sample t-test via non-central t."""
    df = 2 * n_per_arm - 2
    t_crit = stats.t.ppf(1.0 - alpha / 2.0, df)
    nc = d * np.sqrt(n_per_arm / 2.0)
    # Power = P(|T| > t_crit | nc)
    power = stats.nct.sf(t_crit, df, nc) + stats.nct.cdf(-t_crit, df, nc)
    return float(power)


def run_power_validation(project_dir: Path, n_sims: int = 2000, seed: int = SEED) -> dict:
    rng = np.random.default_rng(seed)
    fig_dir = project_dir / "reports" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    alpha = 0.05
    sample_sizes = [10, 20, 30, 40, 50, 65, 80, 100, 120, 150]
    effect_sizes = [0.3, 0.5, 0.8]

    records = []
    for d in effect_sizes:
        for n in sample_sizes:
            theo = theoretical_power_two_sample_t(n, d, alpha=alpha)
            # Vectorized simulation of n_sims two-arm trials under Normal(0,1) vs Normal(d,1)
            x_ctrl = rng.normal(0.0, 1.0, size=(n_sims, n))
            x_trt = rng.normal(d, 1.0, size=(n_sims, n))
            _, p_vals = stats.ttest_ind(x_trt, x_ctrl, axis=1, equal_var=False)
            emp = float(np.mean(p_vals < alpha))
            mc_se = float(np.sqrt(emp * (1.0 - emp) / n_sims))
            records.append(
                {
                    "effect_size_d": d,
                    "n_per_arm": n,
                    "theoretical_power": round(theo, 4),
                    "empirical_power": round(emp, 4),
                    "abs_gap": round(abs(emp - theo), 4),
                    "mc_se": round(mc_se, 4),
                }
            )

    df_pow = pd.DataFrame(records)

    # 1. Empirical Type I Error Check (Null d = 0.0, n = 80, 5,000 trials)
    n_null_sims = 5000
    x0 = rng.normal(0.0, 1.0, size=(n_null_sims, 80))
    x1 = rng.normal(0.0, 1.0, size=(n_null_sims, 80))
    _, p_null = stats.ttest_ind(x1, x0, axis=1, equal_var=False)
    emp_type1 = float(np.mean(p_null < alpha))

    # 2. Normality-Violation Experiment (d = 0.5, Log-Normal skewed data with same mean/SD)
    # LogNormal(0, 1) has mean = exp(0.5), var = (e - 1)*e -> standardize to mean 0, sd 1
    ln_mean = np.exp(0.5)
    ln_sd = np.sqrt((np.e - 1.0) * np.e)
    skew_records = []
    for n in [15, 30, 50, 80, 100]:
        theo = theoretical_power_two_sample_t(n, 0.5, alpha=alpha)
        raw_c = (rng.lognormal(0.0, 1.0, size=(n_sims, n)) - ln_mean) / ln_sd
        raw_t = (rng.lognormal(0.0, 1.0, size=(n_sims, n)) - ln_mean) / ln_sd + 0.5
        _, p_t = stats.ttest_ind(raw_t, raw_c, axis=1, equal_var=False)
        # Mann-Whitney U on first 800 sims for speed
        mw_hits = sum(
            stats.mannwhitneyu(raw_t[i], raw_c[i], alternative="two-sided").pvalue < alpha
            for i in range(800)
        )
        skew_records.append(
            {
                "n_per_arm": n,
                "theoretical_normal_power": round(theo, 4),
                "empirical_welch_skewed": round(float(np.mean(p_t < alpha)), 4),
                "empirical_mw_skewed": round(mw_hits / 800.0, 4),
            }
        )
    df_skew = pd.DataFrame(skew_records)

    # Required sample size for 80% power at d = 0.5 and d = 0.8
    req_n = {}
    for d_target in [0.3, 0.5, 0.8]:
        for n_cand in range(5, 400):
            if theoretical_power_two_sample_t(n_cand, d_target, alpha=alpha) >= 0.80:
                req_n[f"d_{d_target}"] = n_cand
                break

    # ------------------------------------------------------------------
    # Plot Figure 5: Theory vs. Experiment Power Curves + Skewness Stress Test
    # ------------------------------------------------------------------
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": "#CBD5E1", "axes.linewidth": 0.9})
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0), dpi=220)

    colors = {0.3: "#D97706", 0.5: "#2A9D8F", 0.8: "#0F4C5C"}
    n_smooth = np.linspace(10, 150, 141)
    for d in effect_sizes:
        theo_curve = [theoretical_power_two_sample_t(int(round(ns)), d, alpha=alpha) for ns in n_smooth]
        sub = df_pow[df_pow["effect_size_d"] == d]
        axes[0].plot(n_smooth, theo_curve, color=colors[d], lw=2.0, label=f"Theory (d = {d})")
        axes[0].errorbar(
            sub["n_per_arm"],
            sub["empirical_power"],
            yerr=1.96 * sub["mc_se"],
            fmt="o",
            color=colors[d],
            markersize=5,
            capsize=2.5,
            alpha=0.9,
            label=f"Simulated (d = {d})",
        )
    axes[0].axhline(0.80, color="#DC2626", linestyle="--", lw=1.2, label="80% Power Target")
    axes[0].set_title(
        "Theory vs. Experiment: Normal Outcomes\n(2,000 Simulated Trials per Point, α = 0.05)",
        fontsize=10.5,
        fontweight="bold",
        color="#0F172A",
    )
    axes[0].set_xlabel("Sample Size per Arm (n)")
    axes[0].set_ylabel("Statistical Power (1 − β)")
    axes[0].set_ylim(0.0, 1.05)
    axes[0].legend(loc="lower right", fontsize=8.0, ncol=2, frameon=True, facecolor="white", edgecolor="#CBD5E1")
    axes[0].grid(linestyle="--", alpha=0.4)
    sns.despine(ax=axes[0])

    # Panel 2: Skewed Data Stress Test (d = 0.5)
    axes[1].plot(
        df_skew["n_per_arm"],
        df_skew["theoretical_normal_power"],
        color="#64748B",
        linestyle="--",
        lw=2.0,
        marker="s",
        label="Normal Theory (d = 0.5)",
    )
    axes[1].plot(
        df_skew["n_per_arm"],
        df_skew["empirical_welch_skewed"],
        color="#D97706",
        lw=2.0,
        marker="o",
        label="Welch t-test (Skewed Data)",
    )
    axes[1].plot(
        df_skew["n_per_arm"],
        df_skew["empirical_mw_skewed"],
        color="#0F4C5C",
        lw=2.2,
        marker="^",
        label="Mann-Whitney U (Skewed Data)",
    )
    axes[1].axhline(0.80, color="#DC2626", linestyle=":", lw=1.1)
    axes[1].set_title(
        "When Assumptions Fail: Skewed Biomarker (d = 0.5)\nMann-Whitney U Outperforms Welch's t-test",
        fontsize=10.5,
        fontweight="bold",
        color="#0F172A",
    )
    axes[1].set_xlabel("Sample Size per Arm (n)")
    axes[1].set_ylabel("Empirical Power (1 − β)")
    axes[1].set_ylim(0.0, 1.05)
    axes[1].legend(loc="lower right", fontsize=8.5, frameon=True, facecolor="white", edgecolor="#CBD5E1")
    axes[1].grid(linestyle="--", alpha=0.4)
    sns.despine(ax=axes[1])

    plt.tight_layout()
    fig.savefig(fig_dir / "05_theory_vs_experiment_power.png", dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    return {
        "n_sims_per_point": n_sims,
        "alpha": alpha,
        "empirical_type1_error_null": round(emp_type1, 4),
        "max_abs_gap_normal": round(float(df_pow["abs_gap"].max()), 4),
        "mean_abs_gap_normal": round(float(df_pow["abs_gap"].mean()), 4),
        "required_n_per_arm_80pct": req_n,
        "power_table": records,
        "skewed_stress_test": skew_records,
    }

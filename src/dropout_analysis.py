"""
Project 8: Clinical Trial Analyzer
Module 5 (Stretch Goal): Missing Follow-Up Data & Dropout Bias Analysis

Compares four ways of estimating the High-Dose vs. Placebo treatment effect
on 12-week SBP change when dropout is informative (patients with adverse events
or poor Week-6 response drop out disproportionately in the active arm):
  1. True Latent Full Cohort (Benchmark — unobserved in real life)
  2. Complete-Case Analysis (CCA — drops missing Week-12 patients; biased optimistic!)
  3. Last Observation Carried Forward (LOCF — uses Week-6 SBP for dropouts)
  4. Baseline Observation Carried Forward (BOCF — assumes zero change for dropouts; conservative)
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from stats_analysis import welch_ci, cohens_d_and_ci


def run_dropout_bias_analysis(project_dir: Path) -> dict:
    data_dir = project_dir / "data"
    fig_dir = project_dir / "reports" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(data_dir / "clean_trial_cohort.csv")
    truth = pd.read_csv(data_dir / "latent_ground_truth.csv")

    # Dropout rates by arm
    dropout_summary = {}
    for arm in ["Placebo", "Low_Dose_50mg", "High_Dose_100mg"]:
        sub = df[df["arm"] == arm]
        n_tot = len(sub)
        n_drop = int((sub["completed_w12"] == 0).sum())
        dropout_summary[arm] = {
            "enrolled": n_tot,
            "dropped": n_drop,
            "dropout_pct": round(100.0 * n_drop / n_tot, 1),
        }

    # Build 4 analytical datasets for High-Dose vs. Placebo
    two_arm = df[df["arm"].isin(["Placebo", "High_Dose_100mg"])].copy()
    two_truth = truth[truth["arm"].isin(["Placebo", "High_Dose_100mg"])].copy()

    # 1. Ground Truth (Full Latent Cohort)
    gt_h = two_truth.loc[two_truth["arm"] == "High_Dose_100mg", "sbp_delta_true"].values
    gt_p = two_truth.loc[two_truth["arm"] == "Placebo", "sbp_delta_true"].values

    # 2. Complete-Case Analysis (CCA)
    cca = two_arm.dropna(subset=["sbp_delta"])
    cca_h = cca.loc[cca["arm"] == "High_Dose_100mg", "sbp_delta"].values
    cca_p = cca.loc[cca["arm"] == "Placebo", "sbp_delta"].values

    # 3. LOCF (If Week-12 missing, carry forward Week-6 SBP; if Week-6 also missing/invalid, use Week-0)
    locf_sbp12 = two_arm["sbp_w12"].fillna(two_arm["sbp_w6"]).fillna(two_arm["sbp_w0"])
    two_arm["sbp_delta_locf"] = locf_sbp12 - two_arm["sbp_w0"]
    locf_h = two_arm.loc[two_arm["arm"] == "High_Dose_100mg", "sbp_delta_locf"].values
    locf_p = two_arm.loc[two_arm["arm"] == "Placebo", "sbp_delta_locf"].values

    # 4. BOCF (If Week-12 missing, assume return to baseline -> delta = 0)
    two_arm["sbp_delta_bocf"] = two_arm["sbp_delta"].fillna(0.0)
    bocf_h = two_arm.loc[two_arm["arm"] == "High_Dose_100mg", "sbp_delta_bocf"].values
    bocf_p = two_arm.loc[two_arm["arm"] == "Placebo", "sbp_delta_bocf"].values

    methods = [
        ("True Full Cohort (Benchmark)", gt_h, gt_p),
        ("Complete-Case (CCA)", cca_h, cca_p),
        ("LOCF (Week-6 Carried)", locf_h, locf_p),
        ("BOCF (Zero Change Imputed)", bocf_h, bocf_p),
    ]

    gt_diff = float(np.mean(gt_h) - np.mean(gt_p))
    comparison = []
    for name, h_vals, p_vals in methods:
        w = welch_ci(h_vals, p_vals)
        d = cohens_d_and_ci(h_vals, p_vals)
        bias = w["mean_diff"] - gt_diff
        comparison.append(
            {
                "method": name,
                "n_high": int(len(h_vals)),
                "n_placebo": int(len(p_vals)),
                "mean_diff_mmhg": round(w["mean_diff"], 2),
                "ci_low": round(w["ci_low"], 2),
                "ci_high": round(w["ci_high"], 2),
                "cohens_d": round(d["d"], 3),
                "bias_vs_truth_mmhg": round(bias, 2),
                "p_value": w["p_value"],
            }
        )

    # Plot Figure 6: Forest Plot of Treatment Effect Under Dropout Handling Methods
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": "#CBD5E1", "axes.linewidth": 0.9})
    fig, ax = plt.subplots(figsize=(8.2, 3.8), dpi=220)

    y_pos = np.arange(len(comparison))[::-1]
    colors = ["#0F172A", "#DC2626", "#2A9D8F", "#1E6091"]

    ax.axvline(gt_diff, color="#0F172A", linestyle="--", lw=1.2, alpha=0.7, label=f"True Effect ({gt_diff:.2f} mmHg)")
    for yp, row, col in zip(y_pos, comparison, colors):
        est = row["mean_diff_mmhg"]
        lo, hi = row["ci_low"], row["ci_high"]
        ax.errorbar(
            est,
            yp,
            xerr=[[est - lo], [hi - est]],
            fmt="o",
            color=col,
            markersize=8,
            capsize=5,
            lw=2.2,
        )
        bias_txt = "Benchmark" if "Benchmark" in row["method"] else f"Bias: {row['bias_vs_truth_mmhg']:+.2f} mmHg"
        ax.text(
            hi + 0.45,
            yp,
            f"{est:.2f} mmHg [{lo:.1f}, {hi:.1f}]  ({bias_txt})",
            va="center",
            ha="left",
            fontsize=8.8,
            color="#1E293B",
            fontweight="bold" if "CCA" in row["method"] else "normal",
        )

    ax.set_yticks(y_pos)
    ax.set_yticklabels([r["method"] for r in comparison], fontsize=9.8)
    ax.set_xlabel("Estimated Treatment Difference: High-Dose − Placebo (mmHg, 95% CI)")
    ax.set_xlim(min(r["ci_low"] for r in comparison) - 1.5, max(r["ci_high"] for r in comparison) + 12.5)
    ax.set_title(
        "Stretch Goal: Impact of Informative Dropout on Estimated Drug Efficacy",
        fontsize=11,
        fontweight="bold",
        color="#0F172A",
    )
    ax.grid(axis="x", linestyle="--", alpha=0.4)
    sns.despine(ax=ax)
    plt.tight_layout()
    fig.savefig(fig_dir / "06_dropout_bias_comparison.png", dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    return {"dropout_summary": dropout_summary, "method_comparison": comparison}

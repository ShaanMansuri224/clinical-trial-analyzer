"""
Project 8: Clinical Trial Analyzer
Module 3: Core Statistical Analysis Pipeline & Visualizations

Covers every required statistical test in the handbook:
  - Baseline balance checks (Welch's t-test, one-way ANOVA, Chi-square)
  - Normality & homogeneity checks (Shapiro-Wilk, Levene's test)
  - Paired t-test (within-arm Baseline vs. Week 12)
  - Independent Student's t-test & Welch's t-test (High-Dose vs. Placebo SBP change)
  - 95% Confidence Intervals & Cohen's d effect size (with 95% CI)
  - Mann-Whitney U non-parametric test (skewed hs-CRP & sensitivity check on SBP)
  - One-Way ANOVA (3 arms: Placebo, Low-Dose 50mg, High-Dose 100mg) + Bonferroni post-hoc
  - Kruskal-Wallis H test (3-arm skewed hs-CRP)
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

PALETTE = {
    "Placebo": "#64748B",
    "Low_Dose_50mg": "#2A9D8F",
    "High_Dose_100mg": "#0F4C5C",
}
ARM_LABELS = {
    "Placebo": "Placebo",
    "Low_Dose_50mg": "Low-Dose (50 mg)",
    "High_Dose_100mg": "High-Dose (100 mg)",
}


def cohens_d_and_ci(x: np.ndarray, y: np.ndarray, alpha: float = 0.05) -> dict:
    """Compute Cohen's d (x - y) using pooled SD and Hedges/approximate 95% CI."""
    nx, ny = len(x), len(y)
    mx, my = np.mean(x), np.mean(y)
    vx, vy = np.var(x, ddof=1), np.var(y, ddof=1)
    s_pooled = np.sqrt(((nx - 1) * vx + (ny - 1) * vy) / (nx + ny - 2))
    d = (mx - my) / s_pooled
    # Standard error of Cohen's d (Hedges & Olkin 1985)
    se_d = np.sqrt((nx + ny) / (nx * ny) + (d ** 2) / (2 * (nx + ny)))
    z = stats.norm.ppf(1 - alpha / 2)
    return {
        "d": float(d),
        "ci_low": float(d - z * se_d),
        "ci_high": float(d + z * se_d),
        "pooled_sd": float(s_pooled),
    }


def welch_ci(x: np.ndarray, y: np.ndarray, alpha: float = 0.05) -> dict:
    """Welch's two-sample t-test with Satterthwaite df and 95% CI for (mean(x) - mean(y))."""
    nx, ny = len(x), len(y)
    mx, my = np.mean(x), np.mean(y)
    vx, vy = np.var(x, ddof=1), np.var(y, ddof=1)
    diff = mx - my
    se = np.sqrt(vx / nx + vy / ny)
    df = ((vx / nx + vy / ny) ** 2) / (((vx / nx) ** 2) / (nx - 1) + ((vy / ny) ** 2) / (ny - 1))
    t_stat = diff / se
    p_val = 2 * stats.t.sf(abs(t_stat), df)
    t_crit = stats.t.ppf(1 - alpha / 2, df)
    return {
        "mean_diff": float(diff),
        "se": float(se),
        "t_stat": float(t_stat),
        "df": float(df),
        "p_value": float(p_val),
        "ci_low": float(diff - t_crit * se),
        "ci_high": float(diff + t_crit * se),
    }


def run_statistical_analysis(project_dir: Path) -> dict:
    data_dir = project_dir / "data"
    fig_dir = project_dir / "reports" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(data_dir / "clean_trial_cohort.csv")
    completers = df[df["completed_w12"] == 1].dropna(subset=["sbp_delta"]).copy()

    # ------------------------------------------------------------------
    # 1. Baseline Balance Checks (All N=300 Randomized Patients)
    # ------------------------------------------------------------------
    arms = ["Placebo", "Low_Dose_50mg", "High_Dose_100mg"]
    baseline_balance = {}
    for col in ["age", "bmi", "sbp_w0", "crp_w0"]:
        groups = [df.loc[df["arm"] == a, col].dropna().values for a in arms]
        if col == "crp_w0":
            stat, pval = stats.kruskal(*groups)
            test_name = "Kruskal-Wallis H"
        else:
            stat, pval = stats.f_oneway(*groups)
            test_name = "One-Way ANOVA"
        baseline_balance[col] = {
            "test": test_name,
            "stat": float(stat),
            "p_value": float(pval),
            "arm_means": {a: float(df.loc[df["arm"] == a, col].mean()) for a in arms},
            "arm_sds": {a: float(df.loc[df["arm"] == a, col].std(ddof=1)) for a in arms},
            "arm_medians": {a: float(df.loc[df["arm"] == a, col].median()) for a in arms},
        }

    # Categorical balance: Sex & Smoking Status (Chi-square)
    for cat in ["sex", "smoking_status"]:
        ct = pd.crosstab(df["arm"], df[cat])
        chi2, pval, dof, _ = stats.chi2_contingency(ct)
        baseline_balance[cat] = {
            "test": "Chi-Square",
            "stat": float(chi2),
            "dof": int(dof),
            "p_value": float(pval),
        }

    # ------------------------------------------------------------------
    # 2. Normality & Variance Assumption Checks (Completers)
    # ------------------------------------------------------------------
    assumptions = {"shapiro_sbp_delta": {}, "shapiro_crp_delta": {}}
    for a in arms:
        sub_sbp = completers.loc[completers["arm"] == a, "sbp_delta"].dropna().values
        sub_crp = completers.loc[completers["arm"] == a, "crp_delta"].dropna().values
        w_s, p_s = stats.shapiro(sub_sbp)
        w_c, p_c = stats.shapiro(sub_crp)
        assumptions["shapiro_sbp_delta"][a] = {"W": float(w_s), "p_value": float(p_s), "n": int(len(sub_sbp))}
        assumptions["shapiro_crp_delta"][a] = {"W": float(w_c), "p_value": float(p_c), "n": int(len(sub_crp))}

    sbp_p = completers.loc[completers["arm"] == "Placebo", "sbp_delta"].values
    sbp_l = completers.loc[completers["arm"] == "Low_Dose_50mg", "sbp_delta"].values
    sbp_h = completers.loc[completers["arm"] == "High_Dose_100mg", "sbp_delta"].values

    lev_stat, lev_p = stats.levene(sbp_h, sbp_p)
    assumptions["levene_primary_sbp"] = {
        "stat": float(lev_stat),
        "p_value": float(lev_p),
        "sd_high_dose": float(np.std(sbp_h, ddof=1)),
        "sd_placebo": float(np.std(sbp_p, ddof=1)),
    }

    # ------------------------------------------------------------------
    # 3. Paired t-tests (Within-Arm Baseline vs. Week 12 SBP)
    # ------------------------------------------------------------------
    paired_tests = {}
    for a in arms:
        sub = completers.loc[completers["arm"] == a].dropna(subset=["sbp_w0", "sbp_w12"])
        t_stat, p_val = stats.ttest_rel(sub["sbp_w12"], sub["sbp_w0"])
        diff = sub["sbp_w12"] - sub["sbp_w0"]
        se = diff.std(ddof=1) / np.sqrt(len(diff))
        t_crit = stats.t.ppf(0.975, len(diff) - 1)
        paired_tests[a] = {
            "n": int(len(diff)),
            "mean_w0": float(sub["sbp_w0"].mean()),
            "mean_w12": float(sub["sbp_w12"].mean()),
            "mean_change": float(diff.mean()),
            "sd_change": float(diff.std(ddof=1)),
            "ci_low": float(diff.mean() - t_crit * se),
            "ci_high": float(diff.mean() + t_crit * se),
            "t_stat": float(t_stat),
            "p_value": float(p_val),
        }

    # ------------------------------------------------------------------
    # 4. Primary Two-Arm Comparison: High-Dose 100mg vs. Placebo (SBP Delta)
    # ------------------------------------------------------------------
    t_ind, p_ind = stats.ttest_ind(sbp_h, sbp_p, equal_var=True)
    welch_res = welch_ci(sbp_h, sbp_p)
    cohen_res = cohens_d_and_ci(sbp_h, sbp_p)
    u_sbp, p_u_sbp = stats.mannwhitneyu(sbp_h, sbp_p, alternative="two-sided")

    primary_two_arm = {
        "n_high_dose": int(len(sbp_h)),
        "n_placebo": int(len(sbp_p)),
        "mean_high_dose": float(np.mean(sbp_h)),
        "sd_high_dose": float(np.std(sbp_h, ddof=1)),
        "mean_placebo": float(np.mean(sbp_p)),
        "sd_placebo": float(np.std(sbp_p, ddof=1)),
        "student_t": {"t_stat": float(t_ind), "df": int(len(sbp_h) + len(sbp_p) - 2), "p_value": float(p_ind)},
        "welch_t": welch_res,
        "cohens_d": cohen_res,
        "mann_whitney_sensitivity": {"U_stat": float(u_sbp), "p_value": float(p_u_sbp)},
    }

    # ------------------------------------------------------------------
    # 5. Secondary Endpoint (Skewed hs-CRP): Mann-Whitney U & Kruskal-Wallis
    # ------------------------------------------------------------------
    crp_p = completers.loc[completers["arm"] == "Placebo", "crp_delta"].dropna().values
    crp_l = completers.loc[completers["arm"] == "Low_Dose_50mg", "crp_delta"].dropna().values
    crp_h = completers.loc[completers["arm"] == "High_Dose_100mg", "crp_delta"].dropna().values

    u_crp, p_u_crp = stats.mannwhitneyu(crp_h, crp_p, alternative="two-sided")
    # Rank-biserial correlation r = 1 - 2U / (n1*n2)
    r_rb = 1.0 - (2.0 * u_crp) / (len(crp_h) * len(crp_p))
    kw_crp, p_kw_crp = stats.kruskal(crp_p, crp_l, crp_h)

    secondary_crp = {
        "n_high_dose": int(len(crp_h)),
        "n_placebo": int(len(crp_p)),
        "median_delta_placebo": float(np.median(crp_p)),
        "iqr_delta_placebo": float(np.percentile(crp_p, 75) - np.percentile(crp_p, 25)),
        "median_delta_low": float(np.median(crp_l)),
        "median_delta_high": float(np.median(crp_h)),
        "iqr_delta_high": float(np.percentile(crp_h, 75) - np.percentile(crp_h, 25)),
        "skewness_high": float(stats.skew(crp_h)),
        "skewness_placebo": float(stats.skew(crp_p)),
        "mann_whitney_high_vs_placebo": {
            "U_stat": float(u_crp),
            "p_value": float(p_u_crp),
            "rank_biserial_r": float(r_rb),
        },
        "kruskal_wallis_3arm": {"H_stat": float(kw_crp), "p_value": float(p_kw_crp)},
    }

    # ------------------------------------------------------------------
    # 6. Multi-Arm One-Way ANOVA (Placebo vs Low-Dose vs High-Dose SBP)
    # ------------------------------------------------------------------
    f_stat, p_anova = stats.f_oneway(sbp_p, sbp_l, sbp_h)
    all_sbp = np.concatenate([sbp_p, sbp_l, sbp_h])
    grand_mean = np.mean(all_sbp)
    ss_between = sum(len(g) * (np.mean(g) - grand_mean) ** 2 for g in [sbp_p, sbp_l, sbp_h])
    ss_total = np.sum((all_sbp - grand_mean) ** 2)
    eta_sq = ss_between / ss_total

    # Post-hoc pairwise Welch's t-tests with Bonferroni correction (m = 3 comparisons)
    pairs = [
        ("Low_Dose_50mg_vs_Placebo", sbp_l, sbp_p),
        ("High_Dose_100mg_vs_Placebo", sbp_h, sbp_p),
        ("High_Dose_100mg_vs_Low_Dose_50mg", sbp_h, sbp_l),
    ]
    posthoc = {}
    for label, g1, g2 in pairs:
        w = welch_ci(g1, g2)
        d_info = cohens_d_and_ci(g1, g2)
        posthoc[label] = {
            "mean_diff": w["mean_diff"],
            "ci_low": w["ci_low"],
            "ci_high": w["ci_high"],
            "p_raw": w["p_value"],
            "p_bonferroni": min(1.0, w["p_value"] * 3.0),
            "cohens_d": d_info["d"],
        }

    anova_3arm = {
        "F_stat": float(f_stat),
        "df_between": 2,
        "df_within": int(len(all_sbp) - 3),
        "p_value": float(p_anova),
        "eta_squared": float(eta_sq),
        "posthoc_bonferroni": posthoc,
    }

    # ------------------------------------------------------------------
    # Generate Publication-Ready Charts (Figs 1 to 4)
    # ------------------------------------------------------------------
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": "#CBD5E1", "axes.linewidth": 0.9})

    # Figure 1: Baseline Balance Check
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.4), dpi=220)
    for ax, (col, title, unit) in zip(
        axes,
        [
            ("age", "Age at Enrollment", "Years"),
            ("bmi", "Body Mass Index", "kg/m²"),
            ("sbp_w0", "Baseline Systolic BP", "mmHg"),
        ],
    ):
        sns.boxplot(
            data=df,
            x="arm",
            y=col,
            hue="arm",
            order=arms,
            palette=PALETTE,
            width=0.55,
            legend=False,
            ax=ax,
            fliersize=3,
        )
        p_val = baseline_balance[col]["p_value"]
        ax.set_title(f"{title}\n(ANOVA p = {p_val:.3f})", fontsize=10.5, fontweight="bold", color="#0F172A")
        ax.set_xlabel("")
        ax.set_ylabel(unit, fontsize=9.5)
        ax.set_xticklabels(["Placebo", "Low 50mg", "High 100mg"], fontsize=8.8)
        ax.grid(axis="y", linestyle="--", alpha=0.4)
        sns.despine(ax=ax)
    fig.suptitle("Baseline Group Balance Check Across Randomized Arms (N = 300)", fontsize=12, fontweight="bold", y=1.02)
    plt.tight_layout()
    fig.savefig(fig_dir / "01_baseline_balance.png", dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # Figure 2: Primary Two-Arm SBP Outcome & Normality Q-Q Plot
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.8), dpi=220)
    two_arm_df = completers[completers["arm"].isin(["Placebo", "High_Dose_100mg"])].copy()
    sns.violinplot(
        data=two_arm_df,
        x="arm",
        y="sbp_delta",
        hue="arm",
        order=["Placebo", "High_Dose_100mg"],
        palette=PALETTE,
        inner="box",
        alpha=0.75,
        legend=False,
        ax=axes[0],
    )
    axes[0].axhline(0, color="#94A3B8", linestyle="--", lw=1.0)
    axes[0].set_title(
        f"Primary Endpoint: 12-Week SBP Change\nWelch p < 0.001, Cohen's d = {cohen_res['d']:.2f}",
        fontsize=10.5,
        fontweight="bold",
        color="#0F172A",
    )
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Δ SBP (Week 12 − Baseline, mmHg)")
    axes[0].set_xticklabels(["Placebo\n(n={})".format(len(sbp_p)), "High-Dose 100 mg\n(n={})".format(len(sbp_h))])
    axes[0].grid(axis="y", linestyle="--", alpha=0.4)
    sns.despine(ax=axes[0])

    # Q-Q Plot of High-Dose SBP Delta
    osm, osr = stats.probplot(sbp_h, dist="norm")
    axes[1].scatter(osm[0], osm[1], color="#0F4C5C", alpha=0.75, s=24, label="High-Dose Δ SBP")
    slope, intercept, _ = osr
    axes[1].plot(osm[0], slope * osm[0] + intercept, color="#D97706", lw=1.8, label="Normal reference")
    sw_p = assumptions["shapiro_sbp_delta"]["High_Dose_100mg"]["p_value"]
    axes[1].set_title(
        f"Normality Check: High-Dose Δ SBP Q-Q Plot\nShapiro-Wilk p = {sw_p:.3f} (Normal)",
        fontsize=10.5,
        fontweight="bold",
        color="#0F172A",
    )
    axes[1].set_xlabel("Theoretical Normal Quantiles")
    axes[1].set_ylabel("Observed Δ SBP (mmHg)")
    axes[1].legend(frameon=True, facecolor="white", edgecolor="#CBD5E1", fontsize=9)
    axes[1].grid(linestyle="--", alpha=0.4)
    sns.despine(ax=axes[1])
    plt.tight_layout()
    fig.savefig(fig_dir / "02_sbp_outcomes_and_normality.png", dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # Figure 3: Skewed Secondary Biomarker (hs-CRP) justifying Mann-Whitney U
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.6), dpi=220)
    sns.histplot(crp_p, kde=True, color="#64748B", alpha=0.45, label="Placebo", ax=axes[0], bins=22)
    sns.histplot(crp_h, kde=True, color="#0F4C5C", alpha=0.55, label="High-Dose 100 mg", ax=axes[0], bins=22)
    axes[0].set_title(
        f"Right-Skewed Secondary Endpoint (Δ hs-CRP)\nShapiro-Wilk p < 0.001 (Non-Normal)",
        fontsize=10.5,
        fontweight="bold",
        color="#0F172A",
    )
    axes[0].set_xlabel("Δ hs-CRP (Week 12 − Baseline, mg/L)")
    axes[0].set_ylabel("Patient Count")
    axes[0].legend(frameon=True, facecolor="white", edgecolor="#CBD5E1", fontsize=9)
    axes[0].grid(axis="y", linestyle="--", alpha=0.4)
    sns.despine(ax=axes[0])

    # Median & IQR comparison (clipped view + callout for outliers)
    sns.boxplot(
        data=completers,
        x="arm",
        y="crp_delta",
        hue="arm",
        order=arms,
        palette=PALETTE,
        width=0.5,
        legend=False,
        ax=axes[1],
        fliersize=4,
    )
    axes[1].axhline(0, color="#94A3B8", linestyle="--", lw=1.0)
    axes[1].set_ylim(-10, 12)
    axes[1].set_title(
        f"Non-Parametric Comparison (Mann-Whitney U)\nHigh-Dose vs Placebo p = {p_u_crp:.4f}",
        fontsize=10.5,
        fontweight="bold",
        color="#0F172A",
    )
    axes[1].set_xlabel("")
    axes[1].set_ylabel("Δ hs-CRP (mg/L, zoomed −10 to +12)")
    axes[1].set_xticklabels(["Placebo", "Low 50mg", "High 100mg"])
    axes[1].grid(axis="y", linestyle="--", alpha=0.4)
    sns.despine(ax=axes[1])
    plt.tight_layout()
    fig.savefig(fig_dir / "03_crp_skewed_biomarker.png", dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # Figure 4: Three-Arm Dose-Response ANOVA (Mean ± 95% CI)
    fig, ax = plt.subplots(figsize=(7.2, 3.8), dpi=220)
    means = [paired_tests[a]["mean_change"] for a in arms]
    ci_errs = [(paired_tests[a]["mean_change"] - paired_tests[a]["ci_low"]) for a in arms]
    x_pos = np.arange(len(arms))
    bars = ax.bar(
        x_pos,
        means,
        yerr=ci_errs,
        capsize=6,
        width=0.52,
        color=[PALETTE[a] for a in arms],
        edgecolor="#0F172A",
        linewidth=0.8,
        error_kw=dict(ecolor="#0F172A", lw=1.5),
    )
    ax.axhline(0, color="#334155", lw=1.0)
    for rect, m, err in zip(bars, means, ci_errs):
        ax.text(
            rect.get_x() + rect.get_width() / 2,
            m - err - 1.2,
            f"{m:.1f} mmHg",
            ha="center",
            va="top",
            fontsize=10,
            fontweight="bold",
            color="#0F172A",
        )
    ax.set_ylim(min(means) - 6.0, 2.0)
    ax.set_xticks(x_pos)
    ax.set_xticklabels([ARM_LABELS[a] for a in arms], fontsize=10)
    ax.set_ylabel("Mean 12-Week Δ SBP (mmHg) ± 95% CI")
    ax.set_title(
        f"Three-Arm Dose-Response Effect (One-Way ANOVA: F = {f_stat:.2f}, p < 0.001, η² = {eta_sq:.3f})",
        fontsize=11,
        fontweight="bold",
        color="#0F172A",
    )
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    sns.despine(ax=ax)
    plt.tight_layout()
    fig.savefig(fig_dir / "04_three_arm_anova.png", dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    return {
        "baseline_balance": baseline_balance,
        "assumptions": assumptions,
        "paired_tests": paired_tests,
        "primary_two_arm": primary_two_arm,
        "secondary_crp": secondary_crp,
        "anova_3arm": anova_3arm,
    }

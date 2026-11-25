# Project 8: Clinical Trial Analyzer (CARDIO-BETA Trial)

**Course:** Probability & Statistics | School of Engineering, JG University (BTech CSE)  
**Project Difficulty:** Intermediate  
**Reproducibility Seed:** `SEED = 42`

## Student & Submission Details
- **Name:** MANSURI MAHAMADSHAAN SAJIDBHAI
- **Division:** D2D A
- **Enrolment No.:** 2601442038
- **Roll No.:** 267138
- **Group No.:** 8

---

## Project Overview
This repository implements an end-to-end **Clinical Trial Analyzer** for the simulated 12-week, three-arm **CARDIO-BETA** antihypertensive randomized controlled trial (`N = 300`: Placebo, Low-Dose 50 mg, High-Dose 100 mg). It integrates:
1. **Relational Database (SQLite, 3NF):** Normalized `patients`, `visits`, and `measurements` tables (`sql/schema.sql`), visual ER diagram (`sql/er_diagram.png`), and 6 analytical SQL queries (`sql/analytical_queries.sql`).
2. **Data-Quality Audit:** Automated detection and logging of duplicate rows, impossible physiological entry errors, missing Week-6 labs, acute hs-CRP outliers, and Week-12 patient dropouts (`data/data_quality_log.md`).
3. **Statistical Analysis Pipeline:** Baseline balance checks, Shapiro-Wilk & Levene's assumption checks, paired t-tests, independent Student's & Welch's t-tests, 95% CIs, Cohen's *d*, Mann-Whitney U (for skewed hs-CRP), and 3-arm One-Way ANOVA with Bonferroni post-hoc tests.
4. **Theory vs. Experiment Validation:** 2,000-trial Monte Carlo simulation per sample size comparing empirical power against exact non-central *t* theoretical power (plus Type I error verification and a skewed-distribution stress test).
5. **Stretch Goal (Dropout Bias):** Quantitative comparison of Complete-Case Analysis (CCA), Last Observation Carried Forward (LOCF), and Baseline Observation Carried Forward (BOCF) against latent ground truth under informative attrition.

---

## Repository Structure
```text
project8_clinical_trial_analyzer/
├── README.md
├── requirements.txt
├── data/
│   ├── raw_patients.csv
│   ├── raw_visits.csv
│   ├── raw_measurements.csv
│   ├── clean_trial_cohort.csv
│   ├── clinical_trial.db
│   └── data_quality_log.md
├── sql/
│   ├── schema.sql
│   ├── analytical_queries.sql
│   └── er_diagram.png
├── src/
│   ├── generate_data.py
│   ├── db_pipeline.py
│   ├── stats_analysis.py
│   ├── power_simulation.py
│   ├── dropout_analysis.py
│   └── run_all.py
├── notebooks/
│   └── clinical_trial_analysis.ipynb
└── reports/
    ├── figures/                  # 6 publication-ready PNG charts
    ├── proposal.pdf              # Project Proposal
    ├── final_report.pdf          # 10-page Final Report (9-section template)
    ├── clinical_trial_presentation.pptx  # 8-minute slide deck with speaker notes
    ├── viva_cheat_sheet.md       # Quick 1-minute viva Q&A prep guide
    └── results_summary.json      # Full reproducible numerical outputs
```

---

## Quick Start (Reproduce Everything)
```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Run the full end-to-end pipeline (fixed SEED = 42)
python3 src/run_all.py
```
Running `python3 src/run_all.py` regenerates the raw CSVs, builds `data/clinical_trial.db`, executes all 6 SQL queries, runs all statistical tests and Monte Carlo power simulations, creates all 6 figures + ER diagram, writes `reports/results_summary.json`, and updates `notebooks/clinical_trial_analysis.ipynb`.

---

## Key Verified Results (`SEED = 42`)
| Metric | Result | Interpretation |
| :--- | :--- | :--- |
| **Baseline Balance** | All 6 covariates *p* > 0.05 | Randomization balanced age, BMI, sex, smoking, SBP, and hs-CRP. |
| **Primary ΔSBP (High-Dose vs. Placebo)** | **−9.55 mmHg** (95% CI [−12.65, −6.45]), Welch *p* = 1.09×10⁻⁸ | Statistically decisive and clinically meaningful blood-pressure drop. |
| **Effect Size (Cohen's *d*)** | **−0.971** (95% CI [−1.293, −0.649]) | Large standardized effect size. |
| **Skewed Δhs-CRP (Mann-Whitney U)** | *U* = 851.0, *p* = 5.45×10⁻¹⁷ (*r* = 0.754) | Significant inflammation reduction where parametric t-test fails normality. |
| **3-Arm One-Way ANOVA (ΔSBP)** | *F*(2, 244) = 19.26, *p* = 1.71×10⁻⁸, η² = 0.136 | Clear monotonic dose-response; High-Dose beats both Low-Dose & Placebo. |
| **Theory vs. Experiment Power** | Mean gap **0.61 pp**; Null Type I error = **0.0514** | Simulated power closely matches exact non-central *t* theory. |
| **Informative Dropout Bias (Stretch)** | CCA **−9.55** vs. Truth **−8.65** (Bias **−0.90 mmHg**) | Complete-case analysis exaggerates efficacy; BOCF gives strict **−6.51 mmHg** lower bound. |

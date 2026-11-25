# 1-Minute Viva Q&A Cheat Sheet — Project 8: Clinical Trial Analyzer

Use these crisp, plain-language answers to ace your 1-minute viva and 8-minute presentation.

---

### Q1: Why did you use Welch's t-test instead of Student's t-test for systolic blood pressure?
**Answer:** After confirming 12-week ΔSBP is normally distributed (Shapiro-Wilk *p* ≥ 0.152 in all arms), Levene's test showed borderline unequal variances (*p* = 0.0628; High-Dose SD = 11.30 mmHg vs. Placebo SD = 8.42 mmHg) and unequal completer sizes (76 vs. 91). Welch's t-test does not assume equal variances, so it preserves the true 5% Type I error rate.

### Q2: Why did you use Mann-Whitney U for hs-CRP?
**Answer:** High-sensitivity CRP is strongly right-skewed (skewness = +5.12 in Placebo) with acute-inflammation outliers above 18 mg/L, failing Shapiro-Wilk normality (*p* < 0.0001). Mann-Whitney U compares ranks rather than means, so it isn't distorted by extreme outliers—and our skewed power simulation proved it achieves **96% power** at *n* = 30 vs. only **58.7%** for Welch's t-test.

### Q3: Why wasn't a before-and-after paired t-test enough?
**Answer:** Even our **Placebo** arm dropped significantly by **−4.99 mmHg** (*p* = 1.84×10⁻⁷) due to the placebo effect and regression to the mean. Comparing High-Dose (−14.54 mmHg) against Placebo isolates the true drug effect of **−9.55 mmHg** (95% CI [−12.65, −6.45]).

### Q4: What does Cohen's *d* = −0.971 mean in plain English?
**Answer:** A *p*-value only tells us the difference is unlikely to be chance; Cohen's *d* measures *how big* the effect is. A *d* of −0.971 means the average High-Dose patient improved by nearly one full standard deviation more than the average Placebo patient—a **large**, clinically meaningful benefit (~10 mmHg drop lowers long-term stroke and heart-attack risk).

### Q5: How did your Theory vs. Experiment power check turn out?
**Answer:** We compared exact non-central *t* theoretical power against 2,000 simulated trials per sample size across *d* = 0.3, 0.5, and 0.8. Under normality, empirical power matched theory within **0.61 percentage points** on average (and null Type I error was **5.14%** vs. 5%). The small discrepancies (~1–2 pp) come from finite Monte Carlo sampling noise, whereas under log-normal skewness, normality failure causes Welch's t-test to diverge from normal theory while Mann-Whitney U gains power.

### Q6: What did your Stretch Goal (dropout bias) show?
**Answer:** High-Dose had 24% Week-12 dropout vs. 9% in Placebo, and High-Dose patients with Week-6 adverse events dropped out at **47.6%**. Because non-responders and side-effect patients left disproportionately, Complete-Case Analysis (−9.55 mmHg) overstated the true full-cohort effect (−8.65 mmHg) by **0.90 mmHg**. Imputing dropouts with LOCF (−7.70 mmHg) or strict BOCF (−6.51 mmHg) provides conservative regulatory bounds.

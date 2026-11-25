# Data-Quality Log — CARDIO-BETA Clinical Trial

| Issue Type | Count | Detection Rule | Action Taken |
| :--- | ---: | :--- | :--- |
| Duplicate measurement rows | 4 | Duplicate (measurement_id, visit_id) keys | Retained first record; removed exact duplicates before DB load. |
| Physiologically impossible SBP values | 2 | sbp_mmhg < 70 or sbp_mmhg > 260 mmHg (found -142.0, 385.0) | Flagged is_valid=0 (ERR_SBP_RANGE) and excluded from primary SBP analysis. |
| Negative hs-CRP lab measurement | 1 | crp_mg_l < 0 mg/L (found -3.40) | Flagged is_valid=0 (ERR_CRP_NEGATIVE) and excluded from biomarker analysis. |
| Missing Week-6 hs-CRP lab draws | 5 | crp_mg_l IS NULL on attended visit | Retained row with NULL crp_mg_l; kept valid SBP/DBP readings. |
| High hs-CRP clinical outliers (>18 mg/L) | 8 | crp_mg_l > 18.0 mg/L (consistent with acute inflammation) | Retained (is_valid=1) and analyzed using non-parametric Mann-Whitney U / Kruskal-Wallis. |
| Missing Week-12 follow-up visits (dropouts) | 53 | visit_week = 12 AND attended = 0 | Evaluated via Complete-Case, LOCF, and BOCF sensitivity analyses (Stretch Goal). |

**Summary:** Loaded 300 patients, 900 scheduled visits, and 847 deduplicated measurement rows into `data/clinical_trial.db` (844 valid rows; 3 excluded range errors).
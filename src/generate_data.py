"""
Project 8: Clinical Trial Analyzer
Module 1: Synthetic Messy Clinical Trial Data Generator (CARDIO-BETA Trial)

Generates realistic 3-arm RCT data (Placebo, Low_Dose_50mg, High_Dose_100mg)
with known ground-truth treatment effects, right-skewed secondary biomarker
(hs-CRP), and realistic clinical data messiness:
  - Duplicate measurement rows
  - Physiologically impossible entry errors
  - Genuine clinical outliers
  - Informative Week-12 dropout (MNAR/MAR) linked to adverse events and poor Week-6 response
"""

from pathlib import Path
import numpy as np
import pandas as pd

SEED = 42


def generate_clinical_trial_data(data_dir: Path, seed: int = SEED) -> dict:
    rng = np.random.default_rng(seed)
    data_dir.mkdir(parents=True, exist_ok=True)

    n_per_arm = 100
    arms = ["Placebo", "Low_Dose_50mg", "High_Dose_100mg"]
    sites = ["SITE-AHM-01", "SITE-MUM-02", "SITE-BLR-03", "SITE-DEL-04"]

    # True population parameters (for Theory vs Experiment & Dropout Bias validation)
    # Primary endpoint: 12-week SBP change (mmHg), negative = improvement
    true_sbp_delta_mean = {
        "Placebo": -3.0,
        "Low_Dose_50mg": -8.2,
        "High_Dose_100mg": -12.5,
    }
    # Slightly unequal SDs across arms to motivate Welch's t-test
    true_sbp_delta_sd = {
        "Placebo": 9.2,
        "Low_Dose_50mg": 10.1,
        "High_Dose_100mg": 11.4,
    }
    # Secondary endpoint: multiplicative reduction in skewed hs-CRP (mg/L)
    true_crp_mult = {
        "Placebo": 0.96,
        "Low_Dose_50mg": 0.78,
        "High_Dose_100mg": 0.62,
    }

    patients_rows = []
    visits_rows = []
    measurements_rows = []
    latent_truth_rows = []

    patient_id = 1001
    visit_id = 5001
    meas_id = 9001

    base_dates = pd.date_range("2025-01-10", periods=90, freq="D")

    for arm in arms:
        for _ in range(n_per_arm):
            age = int(np.clip(rng.normal(54, 10), 28, 80))
            sex = rng.choice(["M", "F"], p=[0.52, 0.48])
            bmi = round(float(np.clip(rng.normal(27.8, 4.2), 18.5, 42.0)), 1)
            smoking = rng.choice(["Never", "Former", "Current"], p=[0.55, 0.27, 0.18])
            site = rng.choice(sites)
            enroll_dt = pd.Timestamp(rng.choice(base_dates))

            patients_rows.append(
                {
                    "patient_id": patient_id,
                    "age": age,
                    "sex": sex,
                    "bmi": bmi,
                    "smoking_status": smoking,
                    "arm": arm,
                    "site_id": site,
                    "enrollment_date": enroll_dt.strftime("%Y-%m-%d"),
                }
            )

            # Baseline measurements (balanced across arms via randomization)
            sbp_w0 = float(np.clip(rng.normal(152.0, 9.8), 132.0, 182.0))
            dbp_w0 = float(np.clip(0.55 * sbp_w0 + rng.normal(10.0, 4.5), 75.0, 115.0))
            # Right-skewed baseline hs-CRP (lognormal)
            crp_w0 = float(np.clip(rng.lognormal(mean=1.15, sigma=0.65), 0.3, 25.0))

            # True latent 12-week SBP change & Week-6 intermediate trajectory
            sbp_delta_12 = float(rng.normal(true_sbp_delta_mean[arm], true_sbp_delta_sd[arm]))
            sbp_w12_true = sbp_w0 + sbp_delta_12
            # Week 6 reflects ~60% of final response plus measurement noise
            sbp_w6_true = sbp_w0 + 0.60 * sbp_delta_12 + float(rng.normal(0.0, 4.0))

            dbp_w6_true = dbp_w0 + 0.45 * (sbp_w6_true - sbp_w0) + float(rng.normal(0.0, 2.5))
            dbp_w12_true = dbp_w0 + 0.45 * sbp_delta_12 + float(rng.normal(0.0, 2.5))

            crp_w6_true = float(np.clip(crp_w0 * np.sqrt(true_crp_mult[arm]) * rng.lognormal(0, 0.18), 0.2, 30.0))
            crp_w12_true = float(np.clip(crp_w0 * true_crp_mult[arm] * rng.lognormal(0, 0.22), 0.2, 30.0))

            # Adverse events higher in High-Dose arm
            ae_prob_w6 = {"Placebo": 0.08, "Low_Dose_50mg": 0.14, "High_Dose_100mg": 0.24}[arm]
            ae_w6 = int(rng.random() < ae_prob_w6)
            ae_w12 = int(rng.random() < (ae_prob_w6 * 0.6))

            # Informative Week-12 dropout (Stretch Goal):
            # Patients who had an AE OR poor Week-6 response (delta > -3 mmHg) drop out more often,
            # especially in High_Dose_100mg -> Complete-Case Analysis overestimates efficacy!
            poor_wk6_response = (sbp_w6_true - sbp_w0) > -3.0
            base_drop = {"Placebo": 0.07, "Low_Dose_50mg": 0.10, "High_Dose_100mg": 0.11}[arm]
            if ae_w6 == 1:
                base_drop += 0.20 if arm == "High_Dose_100mg" else 0.10
            if poor_wk6_response and arm != "Placebo":
                base_drop += 0.22 if arm == "High_Dose_100mg" else 0.12
            dropped_w12 = int(rng.random() < min(base_drop, 0.65))

            latent_truth_rows.append(
                {
                    "patient_id": patient_id,
                    "arm": arm,
                    "sbp_w0": round(sbp_w0, 2),
                    "sbp_w6": round(sbp_w6_true, 2),
                    "sbp_w12_true": round(sbp_w12_true, 2),
                    "sbp_delta_true": round(sbp_delta_12, 2),
                    "crp_w0": round(crp_w0, 2),
                    "crp_w12_true": round(crp_w12_true, 2),
                    "crp_delta_true": round(crp_w12_true - crp_w0, 2),
                    "ae_w6": ae_w6,
                    "dropped_w12": dropped_w12,
                }
            )

            # Create 3 scheduled visits (Week 0, Week 6, Week 12)
            visit_specs = [
                (0, 1, 0, sbp_w0, dbp_w0, crp_w0),
                (6, 1, ae_w6, sbp_w6_true, dbp_w6_true, crp_w6_true),
                (12, 1 - dropped_w12, ae_w12 if not dropped_w12 else 0, sbp_w12_true, dbp_w12_true, crp_w12_true),
            ]
            for wk, attended, ae, sbp_val, dbp_val, crp_val in visit_specs:
                v_date = (enroll_dt + pd.Timedelta(days=wk * 7 + int(rng.integers(-2, 3)))).strftime("%Y-%m-%d")
                visits_rows.append(
                    {
                        "visit_id": visit_id,
                        "patient_id": patient_id,
                        "visit_week": wk,
                        "visit_date": v_date if attended else None,
                        "attended": attended,
                        "adverse_event": ae,
                    }
                )
                if attended:
                    measurements_rows.append(
                        {
                            "measurement_id": meas_id,
                            "visit_id": visit_id,
                            "sbp_mmhg": round(sbp_val, 1),
                            "dbp_mmhg": round(dbp_val, 1),
                            "crp_mg_l": round(crp_val, 2),
                        }
                    )
                    meas_id += 1
                visit_id += 1

            patient_id += 1

    df_patients = pd.DataFrame(patients_rows)
    df_visits = pd.DataFrame(visits_rows)
    df_meas = pd.DataFrame(measurements_rows)
    df_truth = pd.DataFrame(latent_truth_rows)

    # ------------------------------------------------------------------
    # Inject realistic clinical data messiness into raw_measurements.csv
    # ------------------------------------------------------------------
    # 1. Impossible physiological entry errors (3 records)
    w6_visit_ids = df_visits.loc[df_visits["visit_week"] == 6, "visit_id"].values
    err_indices = df_meas.index[df_meas["visit_id"].isin(w6_visit_ids[:3])].tolist()
    df_meas.loc[err_indices[0], "sbp_mmhg"] = -142.0   # Negative SBP typo
    df_meas.loc[err_indices[1], "sbp_mmhg"] = 385.0    # Impossible SBP typo
    df_meas.loc[err_indices[2], "crp_mg_l"] = -3.40    # Negative biomarker error

    # 2. Genuine acute-inflammation hs-CRP outliers (4 records, kept after review)
    w12_indices = df_meas.index[
        df_meas["visit_id"].isin(df_visits.loc[df_visits["visit_week"] == 12, "visit_id"])
    ].tolist()
    for idx, spike in zip(w12_indices[10:14], [19.4, 22.8, 18.6, 24.1]):
        df_meas.loc[idx, "crp_mg_l"] = spike

    # 3. Sporadic missing lab values (5 missing hs-CRP lab draws at Week 6)
    for idx in df_meas.index[df_meas["visit_id"].isin(w6_visit_ids[15:20])]:
        df_meas.loc[idx, "crp_mg_l"] = np.nan

    # 4. Exact duplicate measurement rows (4 duplicated rows)
    dup_rows = df_meas.iloc[[12, 45, 110, 205]].copy()
    df_meas = pd.concat([df_meas, dup_rows], ignore_index=True)

    # Save raw CSVs and ground-truth reference
    df_patients.to_csv(data_dir / "raw_patients.csv", index=False)
    df_visits.to_csv(data_dir / "raw_visits.csv", index=False)
    df_meas.to_csv(data_dir / "raw_measurements.csv", index=False)
    df_truth.to_csv(data_dir / "latent_ground_truth.csv", index=False)

    return {
        "n_patients": len(df_patients),
        "n_visits": len(df_visits),
        "n_raw_measurements": len(df_meas),
        "n_duplicates_injected": len(dup_rows),
        "n_range_errors_injected": 3,
        "n_missing_crp_injected": 5,
        "n_crp_outliers_injected": 4,
    }


if __name__ == "__main__":
    out = Path(__file__).resolve().parents[1] / "data"
    info = generate_clinical_trial_data(out)
    print("Generated clinical trial raw data:", info)

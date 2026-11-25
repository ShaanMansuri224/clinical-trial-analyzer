"""
Project 8: Clinical Trial Analyzer
Module 2: Database ETL, Data-Quality Logging, Analytical SQL Execution & ER Diagram
"""

from pathlib import Path
import sqlite3
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches


def build_database_and_clean(project_dir: Path) -> dict:
    data_dir = project_dir / "data"
    sql_dir = project_dir / "sql"
    db_path = data_dir / "clinical_trial.db"

    df_patients = pd.read_csv(data_dir / "raw_patients.csv")
    df_visits = pd.read_csv(data_dir / "raw_visits.csv")
    df_meas_raw = pd.read_csv(data_dir / "raw_measurements.csv")

    qc_events = []

    # 1. Detect and drop exact duplicate measurement rows
    n_dups = int(df_meas_raw.duplicated(subset=["measurement_id", "visit_id"]).sum())
    df_meas = df_meas_raw.drop_duplicates(subset=["measurement_id", "visit_id"], keep="first").copy()
    qc_events.append(
        {
            "issue": "Duplicate measurement rows",
            "count": n_dups,
            "detection_rule": "Duplicate (measurement_id, visit_id) keys",
            "action_taken": "Retained first record; removed exact duplicates before DB load.",
        }
    )

    # 2. Range & physiological validity checks
    df_meas["is_valid"] = 1
    df_meas["qc_flag"] = "OK"

    bad_sbp = df_meas["sbp_mmhg"].notna() & ((df_meas["sbp_mmhg"] < 70) | (df_meas["sbp_mmhg"] > 260))
    n_bad_sbp = int(bad_sbp.sum())
    df_meas.loc[bad_sbp, "is_valid"] = 0
    df_meas.loc[bad_sbp, "qc_flag"] = "ERR_SBP_RANGE"
    qc_events.append(
        {
            "issue": "Physiologically impossible SBP values",
            "count": n_bad_sbp,
            "detection_rule": "sbp_mmhg < 70 or sbp_mmhg > 260 mmHg (found -142.0, 385.0)",
            "action_taken": "Flagged is_valid=0 (ERR_SBP_RANGE) and excluded from primary SBP analysis.",
        }
    )

    bad_crp = df_meas["crp_mg_l"].notna() & (df_meas["crp_mg_l"] < 0)
    n_bad_crp = int(bad_crp.sum())
    df_meas.loc[bad_crp, "is_valid"] = 0
    df_meas.loc[bad_crp, "qc_flag"] = "ERR_CRP_NEGATIVE"
    qc_events.append(
        {
            "issue": "Negative hs-CRP lab measurement",
            "count": n_bad_crp,
            "detection_rule": "crp_mg_l < 0 mg/L (found -3.40)",
            "action_taken": "Flagged is_valid=0 (ERR_CRP_NEGATIVE) and excluded from biomarker analysis.",
        }
    )

    miss_crp = df_meas["crp_mg_l"].isna()
    n_miss_crp = int(miss_crp.sum())
    df_meas.loc[miss_crp, "qc_flag"] = "MISSING_CRP_LAB"
    qc_events.append(
        {
            "issue": "Missing Week-6 hs-CRP lab draws",
            "count": n_miss_crp,
            "detection_rule": "crp_mg_l IS NULL on attended visit",
            "action_taken": "Retained row with NULL crp_mg_l; kept valid SBP/DBP readings.",
        }
    )

    # 3. Genuine high hs-CRP outliers (>18 mg/L) — flagged for awareness, NOT deleted
    high_crp = df_meas["crp_mg_l"].notna() & (df_meas["crp_mg_l"] > 18.0) & (df_meas["is_valid"] == 1)
    n_high_crp = int(high_crp.sum())
    df_meas.loc[high_crp, "qc_flag"] = "OUTLIER_CRP_HIGH"
    qc_events.append(
        {
            "issue": "High hs-CRP clinical outliers (>18 mg/L)",
            "count": n_high_crp,
            "detection_rule": "crp_mg_l > 18.0 mg/L (consistent with acute inflammation)",
            "action_taken": "Retained (is_valid=1) and analyzed using non-parametric Mann-Whitney U / Kruskal-Wallis.",
        }
    )

    # 4. Missing Week-12 follow-up visits (patient dropouts)
    n_dropouts = int(((df_visits["visit_week"] == 12) & (df_visits["attended"] == 0)).sum())
    qc_events.append(
        {
            "issue": "Missing Week-12 follow-up visits (dropouts)",
            "count": n_dropouts,
            "detection_rule": "visit_week = 12 AND attended = 0",
            "action_taken": "Evaluated via Complete-Case, LOCF, and BOCF sensitivity analyses (Stretch Goal).",
        }
    )

    # Build SQLite database from schema.sql
    if db_path.exists():
        db_path.unlink()

    conn = sqlite3.connect(db_path)
    schema_sql = (sql_dir / "schema.sql").read_text(encoding="utf-8")
    conn.executescript(schema_sql)

    df_patients.to_sql("patients", conn, if_exists="append", index=False)
    df_visits.to_sql("visits", conn, if_exists="append", index=False)
    df_meas.to_sql("measurements", conn, if_exists="append", index=False)
    conn.commit()

    # Execute all 6 analytical SQL queries and store outputs
    sql_text = (sql_dir / "analytical_queries.sql").read_text(encoding="utf-8")
    raw_blocks = [b.strip() for b in sql_text.split(";") if b.strip()]
    sql_results = {}
    for idx, stmt in enumerate(raw_blocks, start=1):
        # Strip leading comments to verify it's a SELECT/WITH
        lines = [ln for ln in stmt.splitlines() if not ln.strip().startswith("--") and ln.strip()]
        if not lines:
            continue
        clean_stmt = "\n".join(lines)
        df_q = pd.read_sql_query(clean_stmt, conn)
        sql_results[f"query_{idx}"] = df_q.to_dict(orient="records")

    # Export clean wide cohort table via SQL join for Python statistical analysis
    cohort_sql = """
    WITH w0 AS (
        SELECT v.patient_id, m.sbp_mmhg AS sbp_w0, m.dbp_mmhg AS dbp_w0, m.crp_mg_l AS crp_w0
        FROM visits v JOIN measurements m ON v.visit_id = m.visit_id
        WHERE v.visit_week = 0 AND v.attended = 1 AND m.is_valid = 1
    ),
    w6 AS (
        SELECT v.patient_id, v.adverse_event AS ae_w6, m.sbp_mmhg AS sbp_w6, m.crp_mg_l AS crp_w6
        FROM visits v LEFT JOIN measurements m ON v.visit_id = m.visit_id AND m.is_valid = 1
        WHERE v.visit_week = 6
    ),
    w12 AS (
        SELECT v.patient_id, v.attended AS completed_w12, v.adverse_event AS ae_w12,
               m.sbp_mmhg AS sbp_w12, m.dbp_mmhg AS dbp_w12, m.crp_mg_l AS crp_w12
        FROM visits v LEFT JOIN measurements m ON v.visit_id = m.visit_id AND m.is_valid = 1
        WHERE v.visit_week = 12
    )
    SELECT
        p.*,
        w0.sbp_w0, w0.dbp_w0, w0.crp_w0,
        w6.ae_w6, w6.sbp_w6, w6.crp_w6,
        w12.completed_w12, w12.ae_w12, w12.sbp_w12, w12.dbp_w12, w12.crp_w12,
        (w12.sbp_w12 - w0.sbp_w0) AS sbp_delta,
        (w12.crp_w12 - w0.crp_w0) AS crp_delta
    FROM patients p
    LEFT JOIN w0 ON p.patient_id = w0.patient_id
    LEFT JOIN w6 ON p.patient_id = w6.patient_id
    LEFT JOIN w12 ON p.patient_id = w12.patient_id
    ORDER BY p.patient_id;
    """
    df_cohort = pd.read_sql_query(cohort_sql, conn)
    df_cohort.to_csv(data_dir / "clean_trial_cohort.csv", index=False)
    conn.close()

    # Write data_quality_log.md
    log_lines = [
        "# Data-Quality Log — CARDIO-BETA Clinical Trial",
        "",
        "| Issue Type | Count | Detection Rule | Action Taken |",
        "| :--- | ---: | :--- | :--- |",
    ]
    for ev in qc_events:
        log_lines.append(
            f"| {ev['issue']} | {ev['count']} | {ev['detection_rule']} | {ev['action_taken']} |"
        )
    log_lines.extend(
        [
            "",
            f"**Summary:** Loaded {len(df_patients)} patients, {len(df_visits)} scheduled visits, and "
            f"{len(df_meas)} deduplicated measurement rows into `data/clinical_trial.db` "
            f"({int((df_meas['is_valid'] == 1).sum())} valid rows; {int((df_meas['is_valid'] == 0).sum())} excluded range errors).",
        ]
    )
    (data_dir / "data_quality_log.md").write_text("\n".join(log_lines), encoding="utf-8")

    # Draw clean visual ER diagram
    draw_er_diagram(sql_dir / "er_diagram.png")

    return {"qc_events": qc_events, "sql_results": sql_results, "cohort_rows": len(df_cohort)}


def draw_er_diagram(out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10.5, 4.4), dpi=220)
    ax.set_xlim(0, 10.5)
    ax.set_ylim(0, 4.4)
    ax.axis("off")

    tables = [
        (
            0.3,
            0.55,
            2.85,
            3.35,
            "PATIENTS",
            "#0F4C5C",
            [
                "PK  patient_id : INTEGER",
                "    age : INTEGER",
                "    sex : TEXT (M/F)",
                "    bmi : REAL",
                "    smoking_status : TEXT",
                "    arm : TEXT (3 arms)",
                "    site_id : TEXT",
                "    enrollment_date : TEXT",
            ],
        ),
        (
            3.85,
            0.95,
            2.8,
            2.65,
            "VISITS",
            "#1E6091",
            [
                "PK  visit_id : INTEGER",
                "FK  patient_id : INTEGER",
                "    visit_week : INT (0,6,12)",
                "    visit_date : TEXT",
                "    attended : INT (0/1)",
                "    adverse_event : INT (0/1)",
            ],
        ),
        (
            7.35,
            0.75,
            2.85,
            2.95,
            "MEASUREMENTS",
            "#2A7F62",
            [
                "PK  measurement_id : INTEGER",
                "FK  visit_id : INTEGER",
                "    sbp_mmhg : REAL",
                "    dbp_mmhg : REAL",
                "    crp_mg_l : REAL",
                "    is_valid : INT (0/1)",
                "    qc_flag : TEXT",
            ],
        ),
    ]

    for x, y, w, h, title, color, fields in tables:
        box = patches.FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.04,rounding_size=0.12",
            facecolor="#F8FAFC",
            edgecolor=color,
            linewidth=2.0,
        )
        ax.add_patch(box)
        header = patches.FancyBboxPatch(
            (x, y + h - 0.52),
            w,
            0.52,
            boxstyle="round,pad=0.04,rounding_size=0.12",
            facecolor=color,
            edgecolor=color,
            linewidth=1.2,
        )
        ax.add_patch(header)
        ax.text(
            x + w / 2,
            y + h - 0.26,
            title,
            ha="center",
            va="center",
            color="white",
            fontsize=11,
            fontweight="bold",
            family="DejaVu Sans",
        )
        for i, fld in enumerate(fields):
            fy = y + h - 0.82 - i * 0.34
            is_key = fld.startswith("PK") or fld.startswith("FK")
            ax.text(
                x + 0.14,
                fy,
                fld,
                ha="left",
                va="center",
                color="#0F172A" if is_key else "#334155",
                fontsize=9.0,
                fontweight="bold" if is_key else "normal",
                family="DejaVu Sans Mono",
            )

    # Relationships (1:N arrows with cardinality labels)
    ax.annotate(
        "",
        xy=(3.82, 2.35),
        xytext=(3.18, 2.35),
        arrowprops=dict(arrowstyle="-|>", color="#334155", lw=2.0, mutation_scale=14),
    )
    ax.text(3.50, 2.60, "1 : N", ha="center", va="bottom", fontsize=9.5, fontweight="bold", color="#0F4C5C")
    ax.text(3.50, 2.08, "attends", ha="center", va="top", fontsize=8.5, color="#475569")

    ax.annotate(
        "",
        xy=(7.32, 2.35),
        xytext=(6.68, 2.35),
        arrowprops=dict(arrowstyle="-|>", color="#334155", lw=2.0, mutation_scale=14),
    )
    ax.text(7.00, 2.60, "1 : 0..1", ha="center", va="bottom", fontsize=9.5, fontweight="bold", color="#1E6091")
    ax.text(7.00, 2.08, "records", ha="center", va="top", fontsize=8.5, color="#475569")

    ax.text(
        5.25,
        4.18,
        "CARDIO-BETA Clinical Trial — Normalized 3NF Relational Schema",
        ha="center",
        va="center",
        fontsize=12.5,
        fontweight="bold",
        color="#0F172A",
    )
    plt.tight_layout()
    fig.savefig(out_path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)

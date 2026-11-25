-- Project 8: Clinical Trial Analyzer
-- Six Meaningful Analytical SQL Queries (Joins, Aggregations, Window Functions, Cohort Filtering)

-- ============================================================================
-- QUERY 1: Baseline Demographic & Clinical Balance Summary by Treatment Arm
-- ============================================================================
SELECT
    p.arm,
    COUNT(DISTINCT p.patient_id) AS n_patients,
    ROUND(AVG(p.age), 2) AS mean_age,
    ROUND(AVG(p.bmi), 2) AS mean_bmi,
    ROUND(100.0 * SUM(CASE WHEN p.sex = 'F' THEN 1 ELSE 0 END) / COUNT(DISTINCT p.patient_id), 1) AS pct_female,
    ROUND(AVG(m.sbp_mmhg), 2) AS mean_baseline_sbp,
    ROUND(AVG(m.crp_mg_l), 2) AS mean_baseline_crp
FROM patients p
JOIN visits v ON p.patient_id = v.patient_id AND v.visit_week = 0
JOIN measurements m ON v.visit_id = m.visit_id AND m.is_valid = 1
GROUP BY p.arm
ORDER BY CASE p.arm
    WHEN 'Placebo' THEN 1
    WHEN 'Low_Dose_50mg' THEN 2
    WHEN 'High_Dose_100mg' THEN 3
END;

-- ============================================================================
-- QUERY 2: Cohort Retention & Dropout Rates by Arm and Clinical Site
-- ============================================================================
SELECT
    p.arm,
    p.site_id,
    COUNT(DISTINCT p.patient_id) AS enrolled_n,
    SUM(CASE WHEN v.visit_week = 12 AND v.attended = 1 THEN 1 ELSE 0 END) AS completed_wk12,
    SUM(CASE WHEN v.visit_week = 12 AND v.attended = 0 THEN 1 ELSE 0 END) AS dropped_out,
    ROUND(100.0 * SUM(CASE WHEN v.visit_week = 12 AND v.attended = 0 THEN 1 ELSE 0 END) / COUNT(DISTINCT p.patient_id), 1) AS dropout_pct
FROM patients p
JOIN visits v ON p.patient_id = v.patient_id
WHERE v.visit_week = 12
GROUP BY p.arm, p.site_id
HAVING COUNT(DISTINCT p.patient_id) >= 5
ORDER BY p.arm, dropout_pct DESC;

-- ============================================================================
-- QUERY 3: Paired Baseline-to-Week-12 Outcome Change (Complete-Case Cohort)
-- ============================================================================
WITH baseline AS (
    SELECT v.patient_id, m.sbp_mmhg AS sbp_w0, m.crp_mg_l AS crp_w0
    FROM visits v
    JOIN measurements m ON v.visit_id = m.visit_id
    WHERE v.visit_week = 0 AND v.attended = 1 AND m.is_valid = 1
),
week12 AS (
    SELECT v.patient_id, m.sbp_mmhg AS sbp_w12, m.crp_mg_l AS crp_w12
    FROM visits v
    JOIN measurements m ON v.visit_id = m.visit_id
    WHERE v.visit_week = 12 AND v.attended = 1 AND m.is_valid = 1
)
SELECT
    p.arm,
    COUNT(*) AS n_completers,
    ROUND(AVG(b.sbp_w0), 2) AS mean_sbp_w0,
    ROUND(AVG(w.sbp_w12), 2) AS mean_sbp_w12,
    ROUND(AVG(w.sbp_w12 - b.sbp_w0), 2) AS mean_sbp_change,
    ROUND(AVG(w.crp_w12 - b.crp_w0), 2) AS mean_crp_change
FROM patients p
JOIN baseline b ON p.patient_id = b.patient_id
JOIN week12 w ON p.patient_id = w.patient_id
GROUP BY p.arm
ORDER BY mean_sbp_change ASC;

-- ============================================================================
-- QUERY 4: Window Functions — Patient Visit-to-Visit SBP Delta & Arm Rank
-- ============================================================================
WITH valid_series AS (
    SELECT
        p.patient_id,
        p.arm,
        v.visit_week,
        m.sbp_mmhg,
        LAG(m.sbp_mmhg) OVER (PARTITION BY p.patient_id ORDER BY v.visit_week) AS prev_sbp,
        FIRST_VALUE(m.sbp_mmhg) OVER (PARTITION BY p.patient_id ORDER BY v.visit_week) AS baseline_sbp
    FROM patients p
    JOIN visits v ON p.patient_id = v.patient_id
    JOIN measurements m ON v.visit_id = m.visit_id
    WHERE v.attended = 1 AND m.is_valid = 1 AND m.sbp_mmhg IS NOT NULL
)
SELECT
    patient_id,
    arm,
    visit_week,
    ROUND(sbp_mmhg, 1) AS sbp_mmhg,
    ROUND(sbp_mmhg - prev_sbp, 1) AS step_delta,
    ROUND(sbp_mmhg - baseline_sbp, 1) AS cumulative_delta,
    ROUND(100.0 * PERCENT_RANK() OVER (PARTITION BY arm, visit_week ORDER BY (sbp_mmhg - baseline_sbp) ASC), 1) AS response_percentile
FROM valid_series
WHERE visit_week = 12
ORDER BY arm, cumulative_delta ASC
LIMIT 15;

-- ============================================================================
-- QUERY 5: Adverse Events & Relationship to Week-12 Dropout by Arm
-- ============================================================================
WITH wk6_ae AS (
    SELECT patient_id, MAX(adverse_event) AS had_ae
    FROM visits
    WHERE visit_week IN (0, 6)
    GROUP BY patient_id
)
SELECT
    p.arm,
    a.had_ae,
    COUNT(*) AS n_patients,
    SUM(CASE WHEN v12.attended = 0 THEN 1 ELSE 0 END) AS n_dropped_wk12,
    ROUND(100.0 * SUM(CASE WHEN v12.attended = 0 THEN 1 ELSE 0 END) / COUNT(*), 1) AS dropout_rate_pct
FROM patients p
JOIN wk6_ae a ON p.patient_id = a.patient_id
JOIN visits v12 ON p.patient_id = v12.patient_id AND v12.visit_week = 12
GROUP BY p.arm, a.had_ae
ORDER BY p.arm, a.had_ae DESC;

-- ============================================================================
-- QUERY 6: Data-Quality Audit — Flagged & Missing Measurements Summary
-- ============================================================================
SELECT
    v.visit_week,
    m.qc_flag,
    COUNT(*) AS record_count,
    SUM(CASE WHEN m.is_valid = 0 THEN 1 ELSE 0 END) AS excluded_count
FROM measurements m
JOIN visits v ON m.visit_id = v.visit_id
GROUP BY v.visit_week, m.qc_flag
ORDER BY v.visit_week, record_count DESC;

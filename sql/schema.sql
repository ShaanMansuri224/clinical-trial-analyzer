-- Project 8: Clinical Trial Analyzer
-- Relational Database Schema (SQLite / Standard SQL)

PRAGMA foreign_keys = ON;

DROP TABLE IF EXISTS measurements;
DROP TABLE IF EXISTS visits;
DROP TABLE IF EXISTS patients;

CREATE TABLE patients (
    patient_id      INTEGER PRIMARY KEY,
    age             INTEGER NOT NULL CHECK (age BETWEEN 18 AND 90),
    sex             TEXT NOT NULL CHECK (sex IN ('M', 'F')),
    bmi             REAL NOT NULL CHECK (bmi BETWEEN 15.0 AND 55.0),
    smoking_status  TEXT NOT NULL CHECK (smoking_status IN ('Never', 'Former', 'Current')),
    arm             TEXT NOT NULL CHECK (arm IN ('Placebo', 'Low_Dose_50mg', 'High_Dose_100mg')),
    site_id         TEXT NOT NULL,
    enrollment_date TEXT NOT NULL
);

CREATE TABLE visits (
    visit_id        INTEGER PRIMARY KEY,
    patient_id      INTEGER NOT NULL,
    visit_week      INTEGER NOT NULL CHECK (visit_week IN (0, 6, 12)),
    visit_date      TEXT,
    attended        INTEGER NOT NULL CHECK (attended IN (0, 1)),
    adverse_event   INTEGER NOT NULL DEFAULT 0 CHECK (adverse_event IN (0, 1)),
    FOREIGN KEY (patient_id) REFERENCES patients(patient_id),
    UNIQUE (patient_id, visit_week)
);

CREATE TABLE measurements (
    measurement_id  INTEGER PRIMARY KEY,
    visit_id        INTEGER NOT NULL,
    sbp_mmhg        REAL,
    dbp_mmhg        REAL,
    crp_mg_l        REAL,
    is_valid        INTEGER NOT NULL DEFAULT 1 CHECK (is_valid IN (0, 1)),
    qc_flag         TEXT DEFAULT 'OK',
    FOREIGN KEY (visit_id) REFERENCES visits(visit_id)
);

CREATE INDEX idx_patients_arm ON patients(arm);
CREATE INDEX idx_visits_patient_week ON visits(patient_id, visit_week);
CREATE INDEX idx_measurements_visit ON measurements(visit_id);

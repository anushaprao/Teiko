"""
load_data.py
Part 1: Relational database schema design and data loading.

Schema rationale:
- projects: top-level entity; one row per project, scales to hundreds of projects
- subjects: one row per patient/subject; foreign key to project; avoids repeating
  demographic data across every sample
- samples: one row per biological sample; foreign key to subject; captures sample
  metadata (type, time point)
- cell_counts: one row per sample containing all five cell-population counts;
  kept in a single wide row (rather than EAV) so aggregate queries remain simple
  and indexed; for hundreds of analyte types a separate analytes/measurements
  EAV pair would be more appropriate

This design avoids redundancy, keeps PK/FK relationships clear, and lets you
add new projects, subjects, or samples without schema changes.
"""

import sqlite3
import csv
import os

DB_PATH = os.path.join(os.path.dirname(__file__), "immune_trial.db")
CSV_PATH = os.path.join(os.path.dirname(__file__), "cell-count.csv")


def init_db(conn: sqlite3.Connection) -> None:
    cur = conn.cursor()
    cur.executescript("""
        PRAGMA journal_mode=WAL;

        CREATE TABLE IF NOT EXISTS projects (
            project_id   TEXT PRIMARY KEY
        );

        CREATE TABLE IF NOT EXISTS subjects (
            subject_id   TEXT PRIMARY KEY,
            project_id   TEXT NOT NULL REFERENCES projects(project_id),
            condition    TEXT,
            age          INTEGER,
            sex          TEXT,
            treatment    TEXT,
            response     TEXT
        );

        CREATE TABLE IF NOT EXISTS samples (
            sample_id               TEXT PRIMARY KEY,
            subject_id              TEXT NOT NULL REFERENCES subjects(subject_id),
            project_id              TEXT NOT NULL REFERENCES projects(project_id),
            sample_type             TEXT,
            time_from_treatment_start INTEGER
        );

        CREATE TABLE IF NOT EXISTS cell_counts (
            sample_id   TEXT PRIMARY KEY REFERENCES samples(sample_id),
            b_cell      INTEGER,
            cd8_t_cell  INTEGER,
            cd4_t_cell  INTEGER,
            nk_cell     INTEGER,
            monocyte    INTEGER
        );

        CREATE INDEX IF NOT EXISTS idx_subjects_project  ON subjects(project_id);
        CREATE INDEX IF NOT EXISTS idx_samples_subject   ON samples(subject_id);
        CREATE INDEX IF NOT EXISTS idx_samples_project   ON samples(project_id);
    """)
    conn.commit()


def load_csv(conn: sqlite3.Connection, csv_path: str) -> None:
    cur = conn.cursor()
    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        projects, subjects, samples, counts = set(), {}, {}, []
        for row in reader:
            proj = row["project"]
            subj = row["subject"]
            samp = row["sample"]

            projects.add(proj)

            if subj not in subjects:
                subjects[subj] = (
                    subj,
                    proj,
                    row["condition"],
                    int(row["age"]) if row["age"] else None,
                    row["sex"],
                    row["treatment"],
                    row["response"] if row["response"] else None,
                )

            if samp not in samples:
                samples[samp] = (
                    samp,
                    subj,
                    proj,
                    row["sample_type"],
                    int(row["time_from_treatment_start"])
                    if row["time_from_treatment_start"] != ""
                    else None,
                )

            counts.append(
                (
                    samp,
                    int(row["b_cell"]),
                    int(row["cd8_t_cell"]),
                    int(row["cd4_t_cell"]),
                    int(row["nk_cell"]),
                    int(row["monocyte"]),
                )
            )

    cur.executemany(
        "INSERT OR IGNORE INTO projects(project_id) VALUES (?)",
        [(p,) for p in sorted(projects)],
    )
    cur.executemany(
        """INSERT OR IGNORE INTO subjects
           (subject_id, project_id, condition, age, sex, treatment, response)
           VALUES (?,?,?,?,?,?,?)""",
        subjects.values(),
    )
    cur.executemany(
        """INSERT OR IGNORE INTO samples
           (sample_id, subject_id, project_id, sample_type, time_from_treatment_start)
           VALUES (?,?,?,?,?)""",
        samples.values(),
    )
    cur.executemany(
        """INSERT OR IGNORE INTO cell_counts
           (sample_id, b_cell, cd8_t_cell, cd4_t_cell, nk_cell, monocyte)
           VALUES (?,?,?,?,?,?)""",
        counts,
    )
    conn.commit()


def main() -> None:
    if os.path.exists(DB_PATH):
        os.remove(DB_PATH)
    conn = sqlite3.connect(DB_PATH)
    init_db(conn)
    load_csv(conn, CSV_PATH)
    cur = conn.cursor()
    cur.execute("SELECT COUNT(*) FROM cell_counts")
    n = cur.fetchone()[0]
    conn.close()
    print(f"Database created at {DB_PATH} ({n} cell-count rows loaded).")


if __name__ == "__main__":
    main()

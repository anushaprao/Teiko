"""
analysis.py
Parts 2-4: frequency table, statistical analysis, subset analysis.
Writes CSV outputs and PNG plots to ./outputs/.
"""

import os
import sqlite3
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

DB_PATH  = os.path.join(os.path.dirname(__file__), "immune_trial.db")
OUT_DIR  = os.path.join(os.path.dirname(__file__), "outputs")
os.makedirs(OUT_DIR, exist_ok=True)

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


# ── helpers ──────────────────────────────────────────────────────────────────

def get_conn() -> sqlite3.Connection:
    return sqlite3.connect(DB_PATH)


# ── Part 2: frequency table ───────────────────────────────────────────────────

def build_frequency_table(conn: sqlite3.Connection) -> pd.DataFrame:
    df = pd.read_sql_query(
        """SELECT s.sample_id AS sample,
                  cc.b_cell, cc.cd8_t_cell, cc.cd4_t_cell, cc.nk_cell, cc.monocyte
           FROM cell_counts cc
           JOIN samples s ON cc.sample_id = s.sample_id""",
        conn,
    )
    df["total_count"] = df[POPULATIONS].sum(axis=1)
    rows = []
    for pop in POPULATIONS:
        tmp = df[["sample", "total_count", pop]].copy()
        tmp["population"] = pop
        tmp["count"] = tmp[pop]
        tmp["percentage"] = (tmp["count"] / tmp["total_count"] * 100).round(4)
        rows.append(tmp[["sample", "total_count", "population", "count", "percentage"]])
    result = pd.concat(rows, ignore_index=True).sort_values(["sample", "population"])
    path = os.path.join(OUT_DIR, "frequency_table.csv")
    result.to_csv(path, index=False)
    print(f"[Part 2] Frequency table saved → {path}  ({len(result)} rows)")
    return result


# ── Part 3: statistical analysis ─────────────────────────────────────────────

def run_statistical_analysis(conn: sqlite3.Connection, freq: pd.DataFrame) -> None:
    # Filter: melanoma, miraclib, PBMC only
    meta = pd.read_sql_query(
        """SELECT s.sample_id, sub.condition, sub.treatment, sub.response,
                  s.sample_type
           FROM samples s
           JOIN subjects sub ON s.subject_id = sub.subject_id
           WHERE sub.condition = 'melanoma'
             AND sub.treatment = 'miraclib'
             AND s.sample_type = 'PBMC'
             AND sub.response IN ('yes','no')""",
        conn,
    )
    merged = freq.merge(meta, left_on="sample", right_on="sample_id")

    # ── boxplot ──
    fig, axes = plt.subplots(1, 5, figsize=(20, 6), sharey=False)
    palette = {"yes": "#2ecc71", "no": "#e74c3c"}
    label_map = {"yes": "Responder", "no": "Non-responder"}
    merged["response_label"] = merged["response"].map(label_map)

    sig_pops = []
    stats_rows = []
    for ax, pop in zip(axes, POPULATIONS):
        sub = merged[merged["population"] == pop]
        yes_vals = sub[sub["response"] == "yes"]["percentage"]
        no_vals  = sub[sub["response"] == "no"]["percentage"]
        stat, pval = stats.mannwhitneyu(yes_vals, no_vals, alternative="two-sided")

        sns.boxplot(
            data=sub, x="response_label", y="percentage",
            hue="response_label",
            palette={"Responder": "#2ecc71", "Non-responder": "#e74c3c"},
            order=["Responder", "Non-responder"],
            legend=False,
            ax=ax, width=0.5, linewidth=1.5,
        )
        ax.set_title(pop.replace("_", " ").title(), fontsize=11, fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel("Relative frequency (%)" if ax == axes[0] else "")
        sig_label = f"p = {pval:.4f}" + (" *" if pval < 0.05 else "")
        ax.text(0.5, 0.97, sig_label, transform=ax.transAxes,
                ha="center", va="top", fontsize=9,
                color="navy" if pval < 0.05 else "gray")
        if pval < 0.05:
            sig_pops.append(pop)
        stats_rows.append({"population": pop, "U_statistic": stat,
                            "p_value": pval, "significant": pval < 0.05})

    fig.suptitle(
        "Melanoma PBMC – Miraclib: Responders vs Non-responders\n"
        "(Mann-Whitney U, two-sided)",
        fontsize=13, fontweight="bold", y=1.02,
    )
    plt.tight_layout()
    plot_path = os.path.join(OUT_DIR, "responder_boxplots.png")
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[Part 3] Boxplot saved → {plot_path}")

    stats_df = pd.DataFrame(stats_rows)
    stats_path = os.path.join(OUT_DIR, "statistical_results.csv")
    stats_df.to_csv(stats_path, index=False)
    print(f"[Part 3] Stats table saved → {stats_path}")
    if sig_pops:
        print(f"[Part 3] Significant populations (p<0.05): {sig_pops}")
    else:
        print("[Part 3] No populations reached p<0.05 significance.")


# ── Part 4: subset analysis ───────────────────────────────────────────────────

def run_subset_analysis(conn: sqlite3.Connection) -> None:
    # Melanoma PBMC baseline miraclib samples
    baseline = pd.read_sql_query(
        """SELECT s.sample_id, s.project_id,
                  sub.subject_id, sub.response, sub.sex
           FROM samples s
           JOIN subjects sub ON s.subject_id = sub.subject_id
           WHERE sub.condition  = 'melanoma'
             AND s.sample_type  = 'PBMC'
             AND s.time_from_treatment_start = 0
             AND sub.treatment  = 'miraclib'""",
        conn,
    )

    # How many samples per project
    samples_per_project = (
        baseline.groupby("project_id")["sample_id"].count()
        .reset_index().rename(columns={"sample_id": "sample_count"})
    )

    # Subject-level (deduplicate)
    subjects = baseline.drop_duplicates("subject_id")
    response_counts = (
        subjects.groupby("response")["subject_id"].count()
        .reset_index().rename(columns={"subject_id": "subject_count"})
    )
    sex_counts = (
        subjects.groupby("sex")["subject_id"].count()
        .reset_index().rename(columns={"subject_id": "subject_count"})
    )

    print("\n[Part 4] Melanoma PBMC baseline miraclib samples")
    print("  Samples per project:\n", samples_per_project.to_string(index=False))
    print("  Responders/Non-responders:\n", response_counts.to_string(index=False))
    print("  Males/Females:\n", sex_counts.to_string(index=False))

    # Save
    samples_per_project.to_csv(os.path.join(OUT_DIR, "p4_samples_per_project.csv"), index=False)
    response_counts.to_csv(os.path.join(OUT_DIR, "p4_response_counts.csv"), index=False)
    sex_counts.to_csv(os.path.join(OUT_DIR, "p4_sex_counts.csv"), index=False)
    print("[Part 4] Subset results saved to outputs/")


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    conn = get_conn()
    freq = build_frequency_table(conn)
    run_statistical_analysis(conn, freq)
    run_subset_analysis(conn)
    conn.close()


if __name__ == "__main__":
    main()

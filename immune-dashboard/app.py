"""
app.py – Interactive dashboard (Flask + Plotly)
"""
import os, sqlite3, json
import pandas as pd
from scipy import stats
from flask import Flask, render_template, jsonify

BASE   = os.path.dirname(__file__)
DB     = os.path.join(BASE, "immune_trial.db")
OUT    = os.path.join(BASE, "outputs")
POPS   = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

app = Flask(__name__, template_folder="templates", static_folder="static")


def conn():
    c = sqlite3.connect(DB)
    c.row_factory = sqlite3.Row
    return c


def load_freq() -> pd.DataFrame:
    return pd.read_csv(os.path.join(OUT, "frequency_table.csv"))


# ── API endpoints ─────────────────────────────────────────────────────────────

@app.route("/api/overview")
def api_overview():
    with conn() as c:
        projects = c.execute("SELECT COUNT(*) n FROM projects").fetchone()["n"]
        subjects = c.execute("SELECT COUNT(*) n FROM subjects").fetchone()["n"]
        samples  = c.execute("SELECT COUNT(*) n FROM samples").fetchone()["n"]
    return jsonify({"projects": projects, "subjects": subjects, "samples": samples})


@app.route("/api/frequency")
def api_frequency():
    freq = load_freq()
    # Return mean pct per population for a quick bar chart
    agg = freq.groupby("population")["percentage"].mean().reset_index()
    agg.columns = ["population", "mean_pct"]
    return jsonify(agg.to_dict(orient="records"))


@app.route("/api/boxplot")
def api_boxplot():
    freq = load_freq()
    with conn() as c:
        meta = pd.read_sql_query(
            """SELECT s.sample_id, sub.condition, sub.treatment,
                      sub.response, s.sample_type
               FROM samples s JOIN subjects sub ON s.subject_id=sub.subject_id
               WHERE sub.condition='melanoma' AND sub.treatment='miraclib'
                 AND s.sample_type='PBMC' AND sub.response IN ('yes','no')""", c)
    merged = freq.merge(meta, left_on="sample", right_on="sample_id")
    out = {}
    for pop in POPS:
        sub = merged[merged["population"]==pop]
        yes = sub[sub["response"]=="yes"]["percentage"].tolist()
        no  = sub[sub["response"]=="no"]["percentage"].tolist()
        _, pval = stats.mannwhitneyu(yes, no, alternative="two-sided")
        out[pop] = {"responder": yes, "non_responder": no, "pval": round(pval,4)}
    return jsonify(out)


@app.route("/api/subset")
def api_subset():
    with conn() as c:
        rows = pd.read_sql_query(
            """SELECT s.project_id, sub.subject_id, sub.response, sub.sex
               FROM samples s JOIN subjects sub ON s.subject_id=sub.subject_id
               WHERE sub.condition='melanoma' AND s.sample_type='PBMC'
                 AND s.time_from_treatment_start=0 AND sub.treatment='miraclib'""", c)
    spp = rows.groupby("project_id")["project_id"].count().reset_index(name="samples")
    subj = rows.drop_duplicates("subject_id")
    resp = subj.groupby("response")["subject_id"].count().reset_index(name="count")
    sex  = subj.groupby("sex")["subject_id"].count().reset_index(name="count")
    return jsonify({
        "samples_per_project": spp.to_dict(orient="records"),
        "response": resp.to_dict(orient="records"),
        "sex": sex.to_dict(orient="records"),
    })


@app.route("/api/stats")
def api_stats():
    df = pd.read_csv(os.path.join(OUT, "statistical_results.csv"))
    return jsonify(df.to_dict(orient="records"))


@app.route("/")
def index():
    return render_template("index.html")


if __name__ == "__main__":
    app.run(debug=False, port=5050, host="0.0.0.0")

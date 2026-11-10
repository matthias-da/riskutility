"""Experiment A, step 2: MST (differentially private) releases + Anonymeter inference risk.

Continuous variables are coarsened with fixed, data-independent bins before MST (the domain is
treated as public); synthetic values are mapped back to bin midpoints. Run with:
  ~/.venvs/rapid-epjds/bin/python -I paper/rapid/A2_mst_anonymeter.py <step>
step = mst | anonymeter
"""
import sys, glob, os, json
import numpy as np, pandas as pd

D = "paper/rapid/data/A"
R = "paper/rapid/results"
BINS = {
    "age": [16, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 91],
    "hours.per.week": [0, 20, 30, 39, 40, 45, 50, 60, 100],
    "capital.gain": [-1, 0, 5000, 10000, 100000],
    "capital.loss": [-1, 0, 1500, 2000, 5000],
}
REP = {"capital.gain": [0, 2500, 7500, 15000], "capital.loss": [0, 750, 1750, 2200]}


def coarsen(df):
    out = df.copy()
    for c, b in BINS.items():
        out[c] = pd.cut(df[c], b, labels=False).astype(str)
    out["education.num"] = df["education.num"].astype(str)
    return out.astype(str)


def decode(df):
    out = df.copy()
    for c, b in BINS.items():
        k = out[c].astype(float).astype(int)
        mids = REP.get(c, [(b[i] + b[i + 1]) / 2 for i in range(len(b) - 1)])
        out[c] = [mids[i] for i in k]
    out["education.num"] = out["education.num"].astype(int)
    return out


def run_mst():
    from snsynth import Synthesizer
    train = pd.read_csv(f"{D}/train.csv")
    tc = coarsen(train)
    for eps in [0.1, 1.0, 10.0]:
        for i in range(1, 6):
            np.random.seed(4000 + i)
            s = Synthesizer.create("mst", epsilon=eps, verbose=False)
            s.fit(tc, preprocessor_eps=0.0, categorical_columns=list(tc.columns))
            syn = decode(s.sample(len(train))[tc.columns])
            syn.to_csv(f"{D}/syn_mst{eps:g}_{i}.csv", index=False)
            print("mst", eps, i, flush=True)


def run_anonymeter():
    from anonymeter.evaluators import InferenceEvaluator
    train = pd.read_csv(f"{D}/train.csv"); hold = pd.read_csv(f"{D}/holdout.csv")
    aux = [c for c in train.columns if c != "income"]
    rows = []
    for f in sorted(glob.glob(f"{D}/syn_*.csv")):
        syn = pd.read_csv(f)[train.columns]
        ev = InferenceEvaluator(ori=train, syn=syn, control=hold, aux_cols=aux,
                                secret="income", n_attacks=2000)
        ev.evaluate(n_jobs=-1)
        r = ev.risk(confidence_level=0.95)
        name = os.path.basename(f)[4:-4]
        rows.append({"release": name, "anon_risk": r.value, "anon_lo": r.ci[0], "anon_hi": r.ci[1]})
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(f"{R}/A_anonymeter.csv", index=False)


if __name__ == "__main__":
    {"mst": run_mst, "anonymeter": run_anonymeter}[sys.argv[1]]()

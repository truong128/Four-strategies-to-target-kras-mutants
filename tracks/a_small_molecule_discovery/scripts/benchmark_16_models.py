#!/usr/bin/env python
"""Benchmark the 16 BRFSS classifiers plus four tabular neural nets.

The 16 names and hyperparameters follow
https://github.com/truong128/brfss-scd-ml-benchmark
(shared_scripts/run_ml_pipeline.py get_classifiers). Gradient boosting is the
histogram booster used by that pipeline on CPU, which is also fast enough to
score a multi-billion library. The four nets are ft_transformer,
ft_transformer_lite, realmlp, and tabular_resnet.
"""
import argparse, json, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd, joblib
from sklearn.base import clone
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import (RandomForestClassifier, ExtraTreesClassifier, AdaBoostClassifier,
                              HistGradientBoostingClassifier, BaggingClassifier)
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis, QuadraticDiscriminantAnalysis
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.isotonic import IsotonicRegression

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import load_yaml, fingerprint_and_desc, scaffold_smiles, enrichment_factor, recall_at_fraction, model_score
from tabular_nets import TabularNetClassifier

warnings.filterwarnings("ignore")

CLASSICAL = [
    "logistic_regression", "multinom_reg", "svm_rbf", "svm_linear", "decision_tree",
    "random_forest", "gradient_boosting", "extra_trees", "adaboost", "bag_tree",
    "bagging", "knn", "lda", "qda", "naive_bayes", "mlp",
]
NEURAL = ["ft_transformer", "ft_transformer_lite", "realmlp", "tabular_resnet"]


def build_classifiers(seed, ml):
    """Exact 16 BRFSS names, plus the four neural nets."""
    dt = DecisionTreeClassifier(max_depth=12, min_samples_leaf=5, class_weight="balanced", random_state=seed)
    dl = dict(epochs=int(ml.get("dl_epochs", 8)), batch_size=int(ml.get("dl_batch_size", 2048)),
              patience=int(ml.get("dl_patience", 3)), learning_rate=float(ml.get("dl_learning_rate", 1e-3)),
              max_rows=int(ml.get("dl_max_rows", 60000)), seed=seed)
    models = {
        "logistic_regression": Pipeline([("s", StandardScaler()), ("m", LogisticRegression(max_iter=2000, class_weight="balanced", random_state=seed))]),
        "multinom_reg": Pipeline([("s", StandardScaler()), ("m", LogisticRegression(solver="saga", l1_ratio=0, max_iter=400, class_weight="balanced", random_state=seed))]),
        "svm_rbf": Pipeline([("s", StandardScaler()), ("m", SVC(kernel="rbf", probability=True, class_weight="balanced", random_state=seed))]),
        "svm_linear": Pipeline([("s", StandardScaler()), ("m", SVC(kernel="linear", probability=True, class_weight="balanced", random_state=seed))]),
        "decision_tree": clone(dt),
        "random_forest": RandomForestClassifier(n_estimators=300, max_depth=12, min_samples_leaf=5, class_weight="balanced", random_state=seed, n_jobs=-1),
        "gradient_boosting": HistGradientBoostingClassifier(max_iter=200, max_depth=4, learning_rate=0.05, random_state=seed),
        "extra_trees": ExtraTreesClassifier(n_estimators=300, max_depth=12, class_weight="balanced", random_state=seed, n_jobs=-1),
        "adaboost": AdaBoostClassifier(n_estimators=100, random_state=seed),
        "bag_tree": BaggingClassifier(
            estimator=DecisionTreeClassifier(max_depth=8, min_samples_leaf=10, class_weight="balanced", random_state=seed),
            n_estimators=25, max_samples=0.8, random_state=seed, n_jobs=-1),
        "bagging": BaggingClassifier(estimator=clone(dt), n_estimators=100, random_state=seed, n_jobs=-1),
        "knn": KNeighborsClassifier(n_neighbors=15, weights="distance"),
        "lda": Pipeline([("s", StandardScaler()), ("m", LinearDiscriminantAnalysis())]),
        "qda": Pipeline([("s", StandardScaler()), ("m", QuadraticDiscriminantAnalysis(solver="eigen", shrinkage="auto", reg_param=0.5))]),
        "naive_bayes": GaussianNB(),
        "mlp": Pipeline([("s", StandardScaler()), ("m", MLPClassifier(hidden_layer_sizes=(64, 32), max_iter=400, early_stopping=True, validation_fraction=0.15, random_state=seed))]),
    }
    for name in NEURAL:
        models[name] = TabularNetClassifier(arch=name, **dl)
    missing = (set(CLASSICAL) | set(NEURAL)) - set(models)
    if missing or len(models) != 20:
        raise RuntimeError(f"Expected 20 classifiers, missing {sorted(missing)}")
    return models


def _subsample(X, d, n, rng):
    if len(d) <= n:
        return X, d
    pos = np.where(d.y.values == 1)[0]
    neg = np.where(d.y.values == 0)[0]
    p = pos if len(pos) <= n // 2 else rng.choice(pos, n // 2, replace=False)
    q = rng.choice(neg, min(n - len(p), len(neg)), replace=False)
    ix = np.r_[p, q]
    rng.shuffle(ix)
    return X[ix], d.iloc[ix]


def run(a):
    C = load_yaml(a.config)
    M = C["ml"]
    seed = int(C["project"].get("seed", 42))
    rng = np.random.default_rng(seed)
    df = pd.read_csv(a.input).dropna(subset=["smiles", "id", "docking_score"]).copy()
    if len(df) < 200:
        raise SystemExit("Need >=200 successful docking rows.")
    pct = float(C["docking"].get("active_percent", 1))
    thr = float(np.percentile(df.docking_score, pct))
    df["y"] = (df.docking_score <= thr).astype(int)
    if df.y.sum() < 20:
        raise SystemExit("Too few virtual actives; increase training set or active_percent.")
    df["scaffold"] = [scaffold_smiles(s) for s in df.smiles]
    fr = M["split"]
    targets = {k: fr[k] * len(df) for k in ["train", "validation", "calibration", "test"]}
    counts = {k: 0 for k in targets}
    assign = {k: [] for k in targets}
    for sc, n in df.groupby("scaffold").size().sort_values(ascending=False).items():
        k = max(targets, key=lambda z: targets[z] - counts[z])
        assign[k].append(sc)
        counts[k] += int(n)
    mp = {s: k for k, L in assign.items() for s in L}
    df["split"] = df.scaffold.map(mp)
    df[["id", "smiles", "docking_score", "y", "scaffold", "split"]].to_csv(a.split_manifest, index=False)

    def fx(sub):
        X, ix = [], []
        for i, s in zip(sub.index, sub.smiles):
            x = fingerprint_and_desc(s, M["radius"], M["nbits"], M["include_descriptors"])
            if x is not None:
                X.append(x)
                ix.append(i)
        return np.asarray(X, dtype=np.float32), sub.loc[ix].copy()

    P = {}
    for sp in ["train", "validation", "calibration", "test"]:
        X, d = fx(df[df.split == sp])
        if len(d) == 0 or d.y.nunique() < 2:
            raise SystemExit(f"Split {sp} lacks both classes.")
        P[sp] = (X, d)

    allmods = build_classifiers(seed, M)
    caps = {name: int(M.get("benchmark_max_rows", 250000)) for name in allmods}
    for name in ("svm_rbf", "svm_linear"):
        caps[name] = int(M.get("svm_max_rows", 5000))
    caps["knn"] = int(M.get("knn_max_rows", 20000))
    for name in ("qda", "mlp"):
        caps[name] = int(M.get("slow_model_max_rows", 20000))
    caps["multinom_reg"] = int(M.get("multinom_max_rows", 40000))
    for name in NEURAL:
        caps[name] = int(M.get("dl_max_rows", 60000))

    Xtr, dtr = P["train"]
    Xv, dv = P["validation"]
    rows = []
    eligible = set(M["production_eligible"])
    for name, mo in allmods.items():
        family = "neural_net" if name in NEURAL else "brfss16"
        Xt, dt = _subsample(Xtr, dtr, caps[name], rng)
        try:
            mo.fit(Xt, dt.y)
            p = model_score(mo, Xv)
            rows.append(dict(
                model=name, family=family, n_train=len(dt),
                roc_auc=roc_auc_score(dv.y, p), average_precision=average_precision_score(dv.y, p),
                ef_1pct=enrichment_factor(dv.y, p, .01), ef_0_1pct=enrichment_factor(dv.y, p, .001),
                recall_top_1pct=recall_at_fraction(dv.y, p, .01),
                production_eligible=name in eligible))
        except Exception as e:
            rows.append(dict(model=name, family=family, error=repr(e), production_eligible=False))
    R = pd.DataFrame(rows)
    R.to_csv(a.comparison, index=False)
    if R.model.nunique() != 20:
        raise SystemExit(f"Comparison has {R.model.nunique()} models, expected 20.")
    metric = M.get("primary_metric", "average_precision")
    E = R[(R.production_eligible == True) & R[metric].notna()].sort_values(metric, ascending=False)
    if E.empty:
        raise SystemExit("No production-eligible model succeeded.")
    winner_name = E.iloc[0].model

    def refit(name):
        Xt, dt = _subsample(Xtr, dtr, caps[name], rng)
        return clone(allmods[name]).fit(Xt, dt.y)

    winner = refit(winner_name)
    top3 = E.head(3).model.tolist()
    mods = [refit(n) for n in top3]
    ps = model_score(winner, Xv)
    pv = np.mean([model_score(m, Xv) for m in mods], axis=0)
    use_vote = average_precision_score(dv.y, pv) > average_precision_score(dv.y, ps)
    prod = mods if use_vote else winner
    prod_name = "soft_vote_top3" if use_vote else winner_name
    Xc, dc = P["calibration"]
    rawc = np.mean([model_score(m, Xc) for m in prod], axis=0) if isinstance(prod, list) else model_score(prod, Xc)
    iso = IsotonicRegression(out_of_bounds="clip").fit(rawc, dc.y)
    pc = iso.transform(rawc)
    alpha0 = np.sort(pc[dc.y.values == 0])
    alpha1 = np.sort(1 - pc[dc.y.values == 1])

    def cp(prob):
        p = np.asarray(prob)
        out = np.zeros((len(p), 2))
        for c, x, ref in [(0, p, alpha0), (1, 1 - p, alpha1)]:
            i = np.searchsorted(ref, x, side="left")
            out[:, c] = (len(ref) - i + 1) / (len(ref) + 1)
        return out

    Xt, dt = P["test"]
    raw = np.mean([model_score(m, Xt) for m in prod], axis=0) if isinstance(prod, list) else model_score(prod, Xt)
    pt = iso.transform(raw)
    pvls = cp(pt)
    eps = float(M["conformal_epsilon"])
    sel = (pvls[:, 1] > eps) & (pvls[:, 0] <= eps)
    met = {
        "n_classifiers": 20,
        "classifiers": CLASSICAL + NEURAL,
        "label_active_percent": pct,
        "docking_threshold": thr,
        "winner_single": winner_name,
        "top3": top3,
        "production_model": prod_name,
        "locked_test": {
            "n": len(dt),
            "prevalence": float(dt.y.mean()),
            "roc_auc": float(roc_auc_score(dt.y, pt)),
            "average_precision": float(average_precision_score(dt.y, pt)),
            "brier": float(brier_score_loss(dt.y, pt)),
            "ef_1pct": enrichment_factor(dt.y, pt, .01),
            "ef_0_1pct": enrichment_factor(dt.y, pt, .001),
            "recall_top_1pct": recall_at_fraction(dt.y, pt, .01),
            "conformal_active_fraction": float(sel.mean()),
            "conformal_active_recall": float(dt.y.values[sel].sum() / max(1, dt.y.sum())),
        },
    }
    json.dump(met, open(a.metrics, "w"), indent=2)
    joblib.dump({
        "production": prod, "production_name": prod_name, "isotonic": iso,
        "alpha0": alpha0, "alpha1": alpha1, "epsilon": eps,
        "feature_config": {"radius": M["radius"], "nbits": M["nbits"], "include_descriptors": M["include_descriptors"]},
        "metrics": met,
    }, a.bundle, compress=3)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    for x in ["input", "config", "comparison", "bundle", "metrics", "split_manifest"]:
        ap.add_argument("--" + x.replace("_", "-"), dest=x, required=True)
    run(ap.parse_args())

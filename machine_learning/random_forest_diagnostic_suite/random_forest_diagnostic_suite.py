"""Random-forest binary-classification diagnostics.

Run ``python random_forest_diagnostic_suite.py --demo`` to create 18 figures
on a synthetic credit-risk-like data set.  For a CSV, pass ``--data`` and
``--target``.  The program deliberately chooses tuning thresholds using only
training out-of-fold predictions; the test set is reserved for reporting.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.calibration import calibration_curve
from sklearn.datasets import make_classification
from sklearn.decomposition import PCA
from sklearn.ensemble import (
    BaggingClassifier,
    GradientBoostingClassifier,
    HistGradientBoostingClassifier,
    RandomForestClassifier,
)
from sklearn.inspection import PartialDependenceDisplay, permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    cohen_kappa_score,
    confusion_matrix,
    f1_score,
    precision_recall_curve,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import (
    StratifiedKFold,
    cross_val_predict,
    cross_val_score,
    learning_curve,
    train_test_split,
)
from sklearn.tree import DecisionTreeClassifier

COLORS = {
    "teal": "#009391", "orange": "#EFB47A", "purple": "#B3A9ED",
    "blue": "#4194C6", "red": "#F6847A", "gray": "#666666",
}
plt.rcParams.update({"figure.dpi": 120, "savefig.dpi": 160, "axes.unicode_minus": False})


def save(fig: plt.Figure, output: Path, name: str) -> None:
    fig.tight_layout()
    fig.savefig(output / name, bbox_inches="tight")
    plt.close(fig)


def demo_data() -> tuple[pd.DataFrame, pd.Series]:
    """Return a reproducible binary classification data set with credit-style names."""
    names = ["age", "education", "workyears", "resideyears", "income", "debtratio", "creditdebt", "otherdebt"]
    x, y = make_classification(
        n_samples=1200, n_features=len(names), n_informative=6, n_redundant=1,
        weights=[0.72, 0.28], class_sep=1.0, random_state=10,
    )
    frame = pd.DataFrame(x, columns=names)
    # Make selected columns resemble common credit variables while retaining signal.
    frame["age"] = 40 + 9 * frame["age"]
    frame["workyears"] = np.clip(12 + 6 * frame["workyears"], 0, None)
    frame["income"] = np.exp(10.4 + 0.35 * frame["income"])
    frame["debtratio"] = np.clip(0.32 + 0.12 * frame["debtratio"], 0, 1)
    return frame, pd.Series(y, name="credit")


def load_data(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.Series]:
    if args.demo:
        return demo_data()
    if not args.data or not args.target:
        raise ValueError("Provide --data CSV --target COLUMN, or use --demo.")
    data = pd.read_csv(args.data)
    if args.target not in data:
        raise ValueError(f"Target column not found: {args.target}")
    y = data[args.target]
    if y.nunique() != 2:
        raise ValueError("This suite currently requires a binary target column.")
    if set(y.dropna().unique()) != {0, 1}:
        values = sorted(y.dropna().unique())
        y = y.map({values[0]: 0, values[1]: 1})
    x = data.drop(columns=[args.target]).select_dtypes(include=np.number)
    if x.empty:
        raise ValueError("No numeric features found after excluding the target.")
    return x, y.astype(int)


def cost_per_sample(y: np.ndarray, prob: np.ndarray, thresholds: np.ndarray, fn_cost: float, fp_cost: float) -> np.ndarray:
    return np.array([
        (fn_cost * ((prob < threshold) & (y == 1)).sum() + fp_cost * ((prob >= threshold) & (y == 0)).sum()) / len(y)
        for threshold in thresholds
    ])


def ks_score(y: np.ndarray, prob: np.ndarray) -> tuple[float, float]:
    points = np.sort(np.unique(prob))
    good = np.array([(prob[y == 0] <= point).mean() for point in points])
    bad = np.array([(prob[y == 1] <= point).mean() for point in points])
    index = int(np.argmax(good - bad))
    return float(good[index] - bad[index]), float(points[index])


def plot_suite(x: pd.DataFrame, y: pd.Series, output: Path, skip_shap: bool) -> None:
    output.mkdir(parents=True, exist_ok=True)
    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.30, stratify=y, random_state=10)
    cv5 = StratifiedKFold(5, shuffle=True, random_state=10)
    rf = RandomForestClassifier(n_estimators=300, max_features="sqrt", class_weight="balanced", oob_score=True, random_state=10, n_jobs=-1)
    rf.fit(x_train, y_train)
    prob = rf.predict_proba(x_test)[:, 1]
    pred = rf.predict(x_test)
    labels = list(x.columns)

    # 01: confusion matrix.
    fig, ax = plt.subplots(figsize=(5, 4))
    sns.heatmap(confusion_matrix(y_test, pred), annot=True, fmt="d", cmap="Blues", ax=ax)
    ax.set(xlabel="Predicted class", ylabel="True class", title="01 Confusion matrix")
    save(fig, output, "01_confusion_matrix.png")

    # 02: OOB and CV tuning curves.
    n_values = list(range(50, 401, 50))
    oob = [RandomForestClassifier(n_estimators=n, max_features="sqrt", oob_score=True, random_state=10, n_jobs=-1).fit(x_train, y_train).oob_score_ for n in n_values]
    feature_values = list(range(1, min(x.shape[1], 10) + 1))
    cv_auc = [cross_val_score(RandomForestClassifier(n_estimators=250, max_features=value, random_state=10, n_jobs=-1), x_train, y_train, cv=cv5, scoring="roc_auc").mean() for value in feature_values]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].plot(n_values, oob, "o-", color=COLORS["teal"]); axes[0].set(title="OOB vs number of trees", xlabel="n_estimators", ylabel="OOB accuracy")
    axes[1].plot(feature_values, cv_auc, "s-", color=COLORS["red"]); axes[1].set(title="CV AUC vs max_features", xlabel="max_features", ylabel="CV AUC")
    save(fig, output, "02_tuning_curves.png")

    # 03: impurity and permutation importance.
    perm = permutation_importance(rf, x_test, y_test, n_repeats=15, scoring="roc_auc", random_state=10, n_jobs=-1)
    order = np.argsort(rf.feature_importances_)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    axes[0].barh(np.arange(len(order)), rf.feature_importances_[order], color=COLORS["teal"]); axes[0].set_yticks(np.arange(len(order)), np.array(labels)[order]); axes[0].set_title("Impurity importance")
    axes[1].barh(np.arange(len(order)), perm.importances_mean[order], xerr=perm.importances_std[order], color=COLORS["blue"]); axes[1].set_title("Permutation importance (test AUC)")
    save(fig, output, "03_feature_importance_comparison.png")

    # 04: PDP + ICE; only use a manageable number of features.
    features_for_pdp = labels[:min(3, len(labels))]
    fig, axes = plt.subplots(1, len(features_for_pdp), figsize=(5 * len(features_for_pdp), 4))
    PartialDependenceDisplay.from_estimator(rf, x_train, features_for_pdp, kind="both", ax=np.atleast_1d(axes), n_jobs=-1)
    fig.suptitle("04 Partial dependence and ICE", y=1.02)
    save(fig, output, "04_pdp_ice.png")

    # 05: learning curve.
    sizes, train_scores, valid_scores = learning_curve(rf, x_train, y_train, cv=cv5, scoring="roc_auc", train_sizes=np.linspace(.15, 1, 7), n_jobs=-1)
    fig, ax = plt.subplots(figsize=(7, 4.5)); ax.plot(sizes, train_scores.mean(1), "o-", label="train", color=COLORS["teal"]); ax.plot(sizes, valid_scores.mean(1), "s-", label="CV", color=COLORS["red"]); ax.fill_between(sizes, valid_scores.mean(1)-valid_scores.std(1), valid_scores.mean(1)+valid_scores.std(1), alpha=.15, color=COLORS["red"]); ax.set(title="05 Learning curve", xlabel="Training samples", ylabel="ROC AUC"); ax.legend()
    save(fig, output, "05_learning_curve.png")

    models = {
        "Logistic": LogisticRegression(max_iter=2000), "Decision tree": DecisionTreeClassifier(max_depth=6, random_state=10),
        "Bagging": BaggingClassifier(n_estimators=200, random_state=10, n_jobs=-1), "Random forest": rf,
        "GBDT": GradientBoostingClassifier(random_state=10), "HistGBDT": HistGradientBoostingClassifier(random_state=10),
    }
    for name, model in models.items():
        if name != "Random forest": model.fit(x_train, y_train)
    probs = {name: model.predict_proba(x_test)[:, 1] for name, model in models.items()}

    # 06: ROC comparison.
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    for name, values in probs.items():
        fpr, tpr, _ = roc_curve(y_test, values); ax.plot(fpr, tpr, label=f"{name} ({roc_auc_score(y_test, values):.3f})")
    ax.plot([0, 1], [0, 1], "k--", lw=1); ax.set(title="06 ROC comparison", xlabel="False positive rate", ylabel="True positive rate"); ax.legend(fontsize=8)
    save(fig, output, "06_multi_model_roc.png")

    # 07: calibration.
    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.plot([0, 1], [0, 1], "k--", label="perfect calibration")
    for name in ("Random forest", "GBDT", "Logistic"):
        observed, predicted = calibration_curve(y_test, probs[name], n_bins=10, strategy="quantile"); ax.plot(predicted, observed, "o-", label=name)
    ax.set(title="07 Calibration curves", xlabel="Predicted probability", ylabel="Observed frequency"); ax.legend()
    save(fig, output, "07_calibration.png")

    shap_values = None
    if not skip_shap:
        try:
            import shap
            xs = x_test.sample(min(150, len(x_test)), random_state=10)
            raw = shap.TreeExplainer(rf).shap_values(xs)
            shap_values = raw[1] if isinstance(raw, list) else (raw[..., 1] if np.ndim(raw) == 3 else raw)
            plt.figure(); shap.summary_plot(shap_values, xs, show=False); plt.title("08 SHAP summary"); plt.savefig(output / "08_shap_summary.png", dpi=160, bbox_inches="tight"); plt.close()
        except ImportError:
            (output / "08_shap_summary.SKIPPED.txt").write_text("Install shap to generate SHAP figures.\n", encoding="utf-8")

    # 09: compact dashboard.
    threshold = .5
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    fpr, tpr, _ = roc_curve(y_test, prob); axes[0, 0].plot(fpr, tpr, color=COLORS["red"]); axes[0, 0].plot([0,1],[0,1],"k--"); axes[0,0].set_title("ROC")
    precision, recall, _ = precision_recall_curve(y_test, prob); axes[0, 1].plot(recall, precision, color=COLORS["blue"]); axes[0,1].set_title("Precision-recall")
    observed, predicted = calibration_curve(y_test, prob, n_bins=10); axes[0,2].plot(predicted, observed, "o-"); axes[0,2].plot([0,1],[0,1],"k--"); axes[0,2].set_title("Calibration")
    sns.heatmap(confusion_matrix(y_test, prob >= threshold), annot=True, fmt="d", cmap="Blues", cbar=False, ax=axes[1,0]); axes[1,0].set_title("Confusion matrix")
    thresholds = np.linspace(.05, .95, 91); axes[1,1].plot(thresholds, [precision_score(y_test, prob >= t, zero_division=0) for t in thresholds], label="precision"); axes[1,1].plot(thresholds, [recall_score(y_test, prob >= t) for t in thresholds], label="recall"); axes[1,1].plot(thresholds, [f1_score(y_test, prob >= t) for t in thresholds], label="F1"); axes[1,1].legend(fontsize=8); axes[1,1].set_title("Threshold metrics")
    axes[1,2].hist(prob[y_test.to_numpy()==0], bins=25, alpha=.6, label="negative"); axes[1,2].hist(prob[y_test.to_numpy()==1], bins=25, alpha=.6, label="positive"); axes[1,2].legend(); axes[1,2].set_title("Probability distribution")
    save(fig, output, "09_evaluation_dashboard.png")

    # 10: three importance measures.
    if shap_values is not None:
        fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True); shap_mean = np.abs(shap_values).mean(0)
        for axis, values, title, color in zip(axes, (rf.feature_importances_, perm.importances_mean, shap_mean), ("Impurity", "Permutation", "mean |SHAP|"), (COLORS["teal"], COLORS["blue"], COLORS["purple"])):
            axis.barh(np.arange(len(order)), values[order], color=color); axis.set_title(title)
        axes[0].set_yticks(np.arange(len(order)), np.array(labels)[order]); save(fig, output, "10_three_importance_views.png")

    # 11: joint tuning heatmap.
    small_n = [50, 100, 200, 300]; grid = np.empty((len(feature_values), len(small_n)))
    for i, mf in enumerate(feature_values):
        for j, n in enumerate(small_n): grid[i, j] = cross_val_score(RandomForestClassifier(n_estimators=n, max_features=mf, random_state=10, n_jobs=-1), x_train, y_train, cv=cv5, scoring="roc_auc").mean()
    fig, ax = plt.subplots(figsize=(8, 5)); sns.heatmap(grid, annot=True, fmt=".3f", xticklabels=small_n, yticklabels=feature_values, cmap="YlGnBu", ax=ax); ax.set(title="11 Joint tuning", xlabel="n_estimators", ylabel="max_features")
    save(fig, output, "11_joint_tuning_heatmap.png")

    # 12: two-dimensional PCA decision surface.
    pca = PCA(n_components=2, random_state=10); tr2, te2 = pca.fit_transform(x_train), pca.transform(x_test); rf2 = RandomForestClassifier(n_estimators=250, class_weight="balanced", random_state=10).fit(tr2, y_train)
    gx, gy = np.meshgrid(np.linspace(te2[:,0].min()-1, te2[:,0].max()+1, 180), np.linspace(te2[:,1].min()-1, te2[:,1].max()+1, 180)); zz = rf2.predict_proba(np.c_[gx.ravel(), gy.ravel()])[:,1].reshape(gx.shape)
    fig, ax = plt.subplots(figsize=(7, 5)); contour=ax.contourf(gx, gy, zz, levels=20, cmap="RdBu_r", alpha=.65); ax.contour(gx, gy, zz, levels=[.5], colors="black"); ax.scatter(te2[:,0], te2[:,1], c=y_test, cmap="bwr", edgecolor="k", s=18); fig.colorbar(contour, ax=ax, label="P(positive)"); ax.set(title="12 PCA decision surface", xlabel="PC1", ylabel="PC2")
    save(fig, output, "12_pca_decision_surface.png")

    # 13 and 17 need OOF probabilities to avoid using test data for threshold selection.
    oof = cross_val_predict(rf, x_train, y_train, cv=cv5, method="predict_proba", n_jobs=-1)[:, 1]
    costs_train, costs_test = cost_per_sample(y_train.to_numpy(), oof, thresholds, 5, 1), cost_per_sample(y_test.to_numpy(), prob, thresholds, 5, 1)
    best_cost_threshold = float(thresholds[np.argmin(costs_train)])
    fig, ax = plt.subplots(figsize=(7, 4.5)); ax.plot(thresholds, costs_train, label="train OOF"); ax.plot(thresholds, costs_test, "--", label="test (report only)"); ax.axvline(best_cost_threshold, color="k", ls="--", label=f"cost threshold={best_cost_threshold:.2f}"); ax.set(title="13 Cost curve", xlabel="Threshold", ylabel="Cost per sample"); ax.legend()
    save(fig, output, "13_cost_curve.png")

    # 14: decision curve.
    net_benefit = [(np.sum((prob >= t) & (y_test.to_numpy()==1)) / len(y_test)) - (np.sum(prob >= t) / len(y_test)) * t/(1-t) for t in thresholds]
    fig, ax = plt.subplots(figsize=(7, 4.5)); ax.plot(thresholds, net_benefit, color=COLORS["teal"], label="random forest"); ax.axhline(0, color="k", lw=1, label="treat none"); ax.set(title="14 Decision curve", xlabel="Threshold probability", ylabel="Net benefit"); ax.legend()
    save(fig, output, "14_decision_curve.png")

    # 15: six-view model comparison.
    fig, axes = plt.subplots(2, 3, figsize=(16, 9)); names = list(models); palette = plt.cm.tab10(np.linspace(0, 1, len(names)))
    for name, color in zip(names, palette):
        values=probs[name]; fpr,tpr,_=roc_curve(y_test,values); axes[0,0].plot(fpr,tpr,label=name,color=color); pr,re,_=precision_recall_curve(y_test,values); axes[0,1].plot(re,pr,label=name,color=color); observed,predicted=calibration_curve(y_test,values,n_bins=8); axes[0,2].plot(predicted,observed,"o-",label=name,color=color)
    axes[0,0].plot([0,1],[0,1],"k--"); axes[0,0].set_title("ROC"); axes[0,1].set_title("PR"); axes[0,2].plot([0,1],[0,1],"k--"); axes[0,2].set_title("Calibration")
    rows=[]
    for name, values in probs.items():
        hard=(values>=.5).astype(int); rows.append([name, roc_auc_score(y_test,values), average_precision_score(y_test,values), f1_score(y_test,hard), cohen_kappa_score(y_test,hard)])
    metrics=pd.DataFrame(rows,columns=["model","AUC","AP","F1","Kappa"]).set_index("model"); metrics.plot(kind="bar", ax=axes[1,0]); axes[1,0].set_title("Metric comparison")
    for name, color in zip(names, palette): axes[1,1].hist(probs[name], bins=20, histtype="step", density=True, label=name, color=color)
    axes[1,1].set_title("Predicted-probability distributions"); axes[1,1].legend(fontsize=7)
    for name, color in zip(names, palette):
        ranking=np.argsort(-probs[name]); captured=np.cumsum(y_test.to_numpy()[ranking])/y_test.sum(); axes[1,2].plot(np.arange(1,len(y_test)+1)/len(y_test),captured,label=name,color=color)
    axes[1,2].plot([0,1],[0,1],"k--"); axes[1,2].set_title("Cumulative capture / lift"); axes[1,2].legend(fontsize=7)
    for axis in axes.ravel(): axis.grid(alpha=.2)
    save(fig, output, "15_six_view_model_comparison.png")

    # 16: global SHAP magnitude plus representative local waterfall-like bars.
    if shap_values is not None:
        xs = x_test.sample(min(150, len(x_test)), random_state=10); mean_abs=np.abs(shap_values).mean(0); high,low=int(np.argmax(shap_values.sum(1))),int(np.argmin(shap_values.sum(1)))
        fig, axes = plt.subplots(1,3,figsize=(16,5)); axes[0].barh(labels,np.sort(mean_abs)); axes[0].set_title("Global mean |SHAP|")
        for axis, index, title in ((axes[1],high,"Highest-risk sample"),(axes[2],low,"Lowest-risk sample")):
            order_local=np.argsort(np.abs(shap_values[index]))[-8:]; axis.barh(np.array(labels)[order_local],shap_values[index][order_local],color=[COLORS["red"] if v>0 else COLORS["teal"] for v in shap_values[index][order_local]]); axis.axvline(0,color="k",lw=1); axis.set_title(title)
        save(fig, output, "16_shap_global_and_local.png")

    # 17: business decision canvas (cost, KS, deciles, approval).
    ks, ks_threshold = ks_score(y_test.to_numpy(), prob); decile=pd.qcut(pd.Series(prob).rank(method="first"),10,labels=False); tab=pd.DataFrame({"decile":9-decile,"target":y_test.to_numpy(),"prob":prob}).groupby("decile").agg(count=("target","size"), event_rate=("target","mean"), mean_probability=("prob","mean")).sort_index(); tab["lift"]=tab.event_rate/y_test.mean()
    fig, axes = plt.subplots(2,2,figsize=(13,9)); axes[0,0].plot(thresholds,costs_train,label="OOF cost"); axes[0,0].axvline(best_cost_threshold,ls="--",color="k"); axes[0,0].set_title("Cost-optimal threshold")
    axes[0,1].bar(np.arange(10),tab.event_rate,color=COLORS["red"]); axes[0,1].plot(np.arange(10),tab.mean_probability,color=COLORS["blue"],marker="o"); axes[0,1].set_title("Risk deciles (D1 = highest risk)")
    sorted_p=np.sort(np.unique(prob)); good=np.array([(prob[y_test.to_numpy()==0]<=v).mean() for v in sorted_p]); bad=np.array([(prob[y_test.to_numpy()==1]<=v).mean() for v in sorted_p]); axes[1,0].plot(sorted_p,good,label="negative CDF"); axes[1,0].plot(sorted_p,bad,label="positive CDF"); axes[1,0].axvline(ks_threshold,ls="--",color="k",label=f"KS={ks:.3f}"); axes[1,0].legend(); axes[1,0].set_title("KS curve")
    approvals=np.array([(prob<t).mean() for t in thresholds]); bad_rates=np.array([y_test.to_numpy()[prob<t].mean() if (prob<t).any() else np.nan for t in thresholds]); scatter=axes[1,1].scatter(approvals,bad_rates,c=thresholds,cmap="RdYlGn_r"); fig.colorbar(scatter,ax=axes[1,1],label="threshold"); axes[1,1].set(title="Approval strategy",xlabel="Approval rate",ylabel="Bad-event rate")
    save(fig, output, "17_business_decision_canvas.png")

    # 18: misclassification diagnosis.
    hard=(prob>=best_cost_threshold).astype(int); groups=np.select([(y_test.to_numpy()==1)&(hard==1),(y_test.to_numpy()==1)&(hard==0),(y_test.to_numpy()==0)&(hard==1)],["TP","FN","FP"],default="TN"); z=(x_test-x_test.mean())/x_test.std(); profile=pd.DataFrame({group:z.loc[groups==group].mean() for group in ("TP","FN","FP","TN")}).T
    fig, axes = plt.subplots(1,3,figsize=(16,5)); pd.Series(groups).value_counts().reindex(["TP","FN","FP","TN"]).plot(kind="bar",ax=axes[0],color=[COLORS["red"],COLORS["orange"],COLORS["purple"],COLORS["teal"]]); axes[0].set_title("Misclassification counts")
    sns.heatmap(profile,center=0,cmap="RdBu_r",ax=axes[1]); axes[1].set_title("Standardized feature profiles")
    for group,color in zip(("TP","FN","FP","TN"),(COLORS["red"],COLORS["orange"],COLORS["purple"],COLORS["teal"])): axes[2].hist(prob[groups==group],bins=20,alpha=.5,label=group,color=color)
    axes[2].axvline(best_cost_threshold,color="k",ls="--"); axes[2].legend(); axes[2].set_title("Probability by outcome group")
    save(fig, output, "18_misclassification_diagnosis.png")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate a random-forest diagnostic figure suite.")
    parser.add_argument("--data", help="CSV data path")
    parser.add_argument("--target", help="binary target column (0/1 or two categories)")
    parser.add_argument("--output", default="rf_diagnostics", help="directory for generated figures")
    parser.add_argument("--demo", action="store_true", help="use a synthetic credit-risk-like data set")
    parser.add_argument("--skip-shap", action="store_true", help="skip optional SHAP figures")
    args = parser.parse_args()
    x, y = load_data(args)
    plot_suite(x, y, Path(args.output), args.skip_shap)
    print(f"Saved diagnostics to {Path(args.output).resolve()}")


if __name__ == "__main__":
    main()

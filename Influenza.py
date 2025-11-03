# =============================================================
# Gene Expression Classification: Influenza A effect on pDCs
# With ANOVA Feature Selection + Detailed Model Reports
# =============================================================
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_validate
from sklearn.metrics import (
    accuracy_score, f1_score, roc_auc_score, classification_report, confusion_matrix
)
from sklearn.decomposition import PCA
from sklearn.feature_selection import SelectKBest, f_classif

# Models
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.neighbors import KNeighborsClassifier
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier



import warnings
warnings.filterwarnings("ignore")


try:
    from gseapy import enrichr
except ImportError:
    enrichr = None
    print("⚠️ gseapy not installed — skipping enrichment analysis.")

# ------------------------------------------------------------
# 1. Load Local Expression Data
# ------------------------------------------------------------
file_path = "GSE30550_series_matrix.txt"
print(f"📂 Loading data from: {file_path}")

with open(file_path, "r") as f:
    lines = f.readlines()

start_idx = None
for i, line in enumerate(lines):
    if line.startswith("ID_REF") or line.count("\t") > 5:
        start_idx = i
        break

if start_idx is None:
    raise ValueError("❌ Could not find a tabular data section in the file.")
print(f"✅ Detected table header at line {start_idx + 1}")

data = pd.read_csv(file_path, sep="\t", header=0, skiprows=start_idx)
data = data.dropna(axis=1, how="all")
print(f"✅ Expression matrix loaded: {data.shape[0]} rows × {data.shape[1]} columns")

# ------------------------------------------------------------
# 2. Preprocessing
# ------------------------------------------------------------
id_col = data.columns[0]
data = data.set_index(id_col)
data = data.apply(pd.to_numeric, errors="coerce")
data = data.dropna(axis=0, thresh=data.shape[1] * 0.5)
data = data.dropna(axis=1, thresh=data.shape[0] * 0.5)

print(f"✅ Clean numeric matrix: {data.shape[0]} genes × {data.shape[1]} samples")

X = data.T.values
sample_names = data.columns
labels = [
    "Infected" if any(x in name.lower() for x in ["hour", "post", "flu"]) else "Control"
    for name in sample_names
]
y = LabelEncoder().fit_transform(labels)

# ------------------------------------------------------------
# 3. Normalization + ANOVA Feature Selection
# ------------------------------------------------------------
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

selector = SelectKBest(score_func=f_classif, k=min(200, X_scaled.shape[1]))
X_selected = selector.fit_transform(X_scaled, y)
selected_genes = data.index[selector.get_support()]
print(f"🧬 Selected top {len(selected_genes)} genes using ANOVA F-test")

X_train, X_test, y_train, y_test = train_test_split(
    X_selected, y, test_size=0.2, random_state=42, stratify=y
)

# ------------------------------------------------------------
# 4. Define Models
# ------------------------------------------------------------
models = {
    "Logistic Regression": LogisticRegression(max_iter=1000),
    "SVM (RBF)": SVC(kernel="rbf", probability=True),
    "Random Forest": RandomForestClassifier(n_estimators=300, random_state=42),
    "Gradient Boosting": GradientBoostingClassifier(random_state=42),
    "KNN": KNeighborsClassifier(n_neighbors=5),
    "XGBoost": XGBClassifier(use_label_encoder=False, eval_metric='logloss', random_state=42),
    "LightGBM": LGBMClassifier(random_state=42, verbose=-1)
}

# ------------------------------------------------------------
# 5. Cross-Validation + Full Evaluation for Each Model
# ------------------------------------------------------------
scoring = {"accuracy": "accuracy", "f1": "f1_weighted", "roc_auc": "roc_auc_ovr"}
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
cv_results = []

for name, model in models.items():
    print(f"\n🔍 Evaluating {name}...")

    # Limit threads for stability on macOS
    if "n_jobs" in model.get_params().keys():
        model.set_params(n_jobs=1)

    try:
        # safer to avoid parallelism here
        scores = cross_validate(model, X_selected, y, cv=cv, scoring=scoring, n_jobs=1)

        mean_acc, mean_f1, mean_auc = (
            scores["test_accuracy"].mean(),
            scores["test_f1"].mean(),
            scores["test_roc_auc"].mean(),
        )
        cv_results.append({
            "Model": name,
            "Accuracy": mean_acc,
            "F1": mean_f1,
            "ROC-AUC": mean_auc
        })
        print(f"✅ CV Accuracy: {mean_acc:.3f} | F1: {mean_f1:.3f} | ROC-AUC: {mean_auc:.3f}")

        # Train & evaluate on holdout
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        y_proba = model.predict_proba(X_test)[:, 1] if hasattr(model, "predict_proba") else None

        acc = accuracy_score(y_test, y_pred)
        f1 = f1_score(y_test, y_pred, average="weighted")
        auc = roc_auc_score(y_test, y_proba) if y_proba is not None else np.nan

        print(f"🧾 Test Accuracy: {acc:.3f} | F1: {f1:.3f} | ROC-AUC: {auc:.3f}")
        print("📋 Classification Report:")
        print(classification_report(y_test, y_pred, target_names=["Control", "Infected"]))

        # Confusion matrix
        cm = confusion_matrix(y_test, y_pred)
        plt.figure(figsize=(5, 4))
        sns.heatmap(cm, annot=True, fmt="d", cmap="coolwarm",
                    xticklabels=["Control", "Infected"],
                    yticklabels=["Control", "Infected"])
        plt.title(f"Confusion Matrix — {name}")
        plt.xlabel("Predicted")
        plt.ylabel("True")
        plt.tight_layout()
        plt.show()

    except Exception as e:
        print(f"⚠️ Skipping {name} due to error: {e}")
        continue

# ------------------------------------------------------------
# 6. Cross-Validation Summary Table
# ------------------------------------------------------------
cv_summary = pd.DataFrame(cv_results).sort_values(by="ROC-AUC", ascending=False)
print("\n📊 Cross-Validation Summary:")
print(cv_summary)

# ------------------------------------------------------------
# 7. PCA Visualization (Top 200 ANOVA Genes)
# ------------------------------------------------------------
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_selected)
plt.figure(figsize=(6, 5))
sns.scatterplot(x=X_pca[:, 0], y=X_pca[:, 1],
                hue=labels, palette="Set2", s=60, alpha=0.8)
plt.title("PCA (Top 200 ANOVA-Selected Genes)")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.tight_layout()
plt.show()

# ------------------------------------------------------------
# 8. Heatmap of Top ANOVA Genes
# ------------------------------------------------------------
plt.figure(figsize=(10, 7))
sns.heatmap(data.loc[selected_genes], cmap="vlag", xticklabels=False)
plt.title("Top 200 ANOVA-Selected Genes — Heatmap")
plt.tight_layout()
plt.show()


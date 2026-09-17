# 02_exploratory_analysis/02_umap.py

from pathlib import Path

import pandas as pd
from sklearn.preprocessing import StandardScaler
import umap
import matplotlib.pyplot as plt

project_root = Path(
    "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"
)

analysis_dataset = project_root / "01_build_analysis_dataset" / "outputs"
out_dir = project_root / "02_exploratory_analysis" / "outputs" / "umap"
out_dir.mkdir(parents=True, exist_ok=True)

input_file = analysis_dataset / "umap_input_clean.csv"
if not input_file.exists():
    raise FileNotFoundError(f"Missing UMAP input file: {input_file}")

df = pd.read_csv(input_file)

if "sample_id" not in df.columns:
    raise ValueError("umap_input_clean.csv must contain a 'sample_id' column")

sample_id = df["sample_id"].copy()
X = df.drop(columns=["sample_id"]).copy()

# Convert everything to numeric and fail fast if anything is missing.
X = X.apply(pd.to_numeric, errors="coerce")
if X.isna().any().any():
    missing_rows = X.index[X.isna().any(axis=1)].tolist()
    raise ValueError(
        f"UMAP input contains missing values. First problematic row index: {missing_rows[0]}"
    )

if X.shape[0] == 0:
    raise ValueError("No rows left in UMAP input.")

X_scaled = StandardScaler().fit_transform(X)

embedding = umap.UMAP(
    n_neighbors=15,
    min_dist=0.1,
    metric="euclidean",
    random_state=42,
).fit_transform(X_scaled)

out = pd.DataFrame(
    {
        "sample_id": sample_id.values,
        "UMAP1": embedding[:, 0],
        "UMAP2": embedding[:, 1],
    }
)

out.to_csv(out_dir / "umap_coords_clean.csv", index=False)

plt.figure(figsize=(8, 6))
plt.scatter(out["UMAP1"], out["UMAP2"], s=10)
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.title("UMAP of metabolomics data")
plt.tight_layout()
plt.savefig(out_dir / "umap_scatter.png", dpi=300)
plt.close()

print(f"Saved UMAP coordinates to: {out_dir / 'umap_coords_clean.csv'}")
print(f"Saved UMAP plot to: {out_dir / 'umap_scatter.png'}")
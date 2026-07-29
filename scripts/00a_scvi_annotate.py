import argparse
import json
import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scvi


def ensure_parent(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)


def to_history_dict(history_obj):
    out = {}
    if history_obj is None:
        return out
    if isinstance(history_obj, dict):
        items = history_obj.items()
    else:
        items = []
    for k, v in items:
        if hasattr(v, "tolist"):
            out[k] = v.tolist()
        elif hasattr(v, "to_list"):
            out[k] = v.to_list()
        elif isinstance(v, (list, tuple)):
            out[k] = list(v)
        else:
            try:
                out[k] = list(v)
            except Exception:
                out[k] = [str(v)]
    return out


def main():
    parser = argparse.ArgumentParser(description="Train scVI/SCANVI on raw data and annotate query cells")
    parser.add_argument("--query-h5ad", required=True)
    parser.add_argument("--reference-h5ad", required=True)
    parser.add_argument("--output-h5ad", required=True)
    parser.add_argument("--history-json", required=True)
    parser.add_argument("--label-key", default="final_label_transfer")
    parser.add_argument("--batch-key", default="sample_id")
    parser.add_argument("--unlabeled-category", default="Unknown")
    parser.add_argument("--n-latent", type=int, default=30)
    parser.add_argument("--max-epochs-scvi", type=int, default=200)
    parser.add_argument("--max-epochs-scanvi", type=int, default=40)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    scvi.settings.seed = int(args.seed)

    query = sc.read_h5ad(args.query_h5ad)
    ref = sc.read_h5ad(args.reference_h5ad)

    if args.label_key not in ref.obs.columns:
        raise ValueError(f"Reference label key not found: {args.label_key}")

    common = ref.var_names.intersection(query.var_names)
    if len(common) < 200:
        raise ValueError(f"Too few shared genes for scVI transfer: {len(common)}")

    ref2 = ref[:, common].copy()
    query2 = query[:, common].copy()

    ref2.obs_names = [f"ref::{x}" for x in ref2.obs_names]
    query2.obs_names = [f"qry::{x}" for x in query2.obs_names]

    ref2.obs[args.label_key] = ref2.obs[args.label_key].astype(str)
    query2.obs[args.label_key] = args.unlabeled_category

    merged = ad.concat([ref2, query2], label="_dataset", keys=["reference", "query"], join="inner")

    effective_batch_key = args.batch_key if args.batch_key in merged.obs.columns else "_dataset"

    scvi.model.SCVI.setup_anndata(merged, batch_key=effective_batch_key)
    model = scvi.model.SCVI(merged, n_latent=int(args.n_latent))
    model.train(max_epochs=int(args.max_epochs_scvi))

    scanvi = scvi.model.SCANVI.from_scvi_model(
        model,
        labels_key=args.label_key,
        unlabeled_category=args.unlabeled_category,
    )
    scanvi.train(max_epochs=int(args.max_epochs_scanvi))

    pred = scanvi.predict(merged).astype(str)
    latent = scanvi.get_latent_representation(merged)

    query_idx = merged.obs["_dataset"].astype(str) == "query"
    merged_query = merged[query_idx].copy()
    pred_q = np.asarray(pred[query_idx])
    lat_q = np.asarray(latent[query_idx])

    out = query.copy()
    out.obs["cluster_pred"] = pred_q
    out.obs[args.label_key] = pred_q
    out.obsm["X_scVI"] = lat_q

    hist = to_history_dict(getattr(scanvi, "history", None))
    out.uns["scvi_history"] = hist

    ensure_parent(args.output_h5ad)
    out.write_h5ad(args.output_h5ad)

    ensure_parent(args.history_json)
    with open(args.history_json, "w", encoding="utf-8") as f:
        json.dump({"history": hist, "n_shared_genes": int(len(common))}, f, indent=2)


if __name__ == "__main__":
    main()

"""Probe 4 -- text-branch benchmark: hash-64 (shipped) vs hash-512 vs
word TF-IDF, all feeding the same Lloyd k-means (k=3) and Cramer's V.

Tasks:
 (a) the shipped B6 sentiment-phrase task (20 seeds);
 (b) an overlapping-vocabulary task: both classes draw most tokens from
     a shared pool, class signal carried by a minority of tokens;
 (c) REAL text: Olist review comments (Portuguese) vs review band
     (score >= 4 positive, <= 2 negative), n = 2000 sampled orders.
"""
import json
import numpy as np
import pandas as pd
import facade_v2 as F

import os
FIELD = os.environ.get("BUKPI_FIELD_DATA", "data/field")  # DEPOSIT COPY: archival absolute path replaced with a configurable one (see MANIFEST.md, "Deposit edits")

def v_of_clusters(clusters, label):
    return F.cramers_v_from_table(F.contingency(pd.Series(clusters),
                                               pd.Series(label)))[0]

def bench(texts, label, seed):
    out = {}
    X64 = np.vstack([F.hash_embed(t, 64) for t in texts])
    out["hash64"] = v_of_clusters(F.kmeans(X64, 3, seed=seed), label)
    X512 = np.vstack([F.hash_embed(t, 512) for t in texts])
    out["hash512"] = v_of_clusters(F.kmeans(X512, 3, seed=seed), label)
    Xt = F.tfidf_embed(texts)
    out["tfidf"] = v_of_clusters(F.kmeans(Xt, 3, seed=seed), label)
    return out

res = {}
# (a) shipped B6 task
pos = ["love this product", "great support team", "fast checkout",
       "very satisfied"]
neg = ["terrible bug", "slow as molasses", "support never replied",
       "lost my data"]
rows = []
for s in range(20):
    rng = np.random.RandomState(1000 + s)
    texts, label = [], []
    for _ in range(500):
        if rng.rand() < 0.5:
            texts.append(rng.choice(pos)); label.append(1)
        else:
            texts.append(rng.choice(neg)); label.append(0)
    rows.append(bench(texts, label, seed=s))
res["task_a_shipped_phrases"] = {
    k: {"mean": float(np.mean([r[k] for r in rows])),
        "sd": float(np.std([r[k] for r in rows]))} for k in rows[0]}

# (b) overlapping vocabulary
shared = ["the order arrived", "package was", "item is", "delivery was",
          "box contained", "service was", "product is", "seller was"]
posmark = ["excellent", "perfect", "delighted"]
negmark = ["broken", "awful", "furious"]
rows = []
for s in range(20):
    rng = np.random.RandomState(2000 + s)
    texts, label = [], []
    for _ in range(500):
        base = rng.choice(shared)
        if rng.rand() < 0.5:
            texts.append(base + " " + rng.choice(posmark)); label.append(1)
        else:
            texts.append(base + " " + rng.choice(negmark)); label.append(0)
    rows.append(bench(texts, label, seed=s))
res["task_b_overlapping_vocab"] = {
    k: {"mean": float(np.mean([r[k] for r in rows])),
        "sd": float(np.std([r[k] for r in rows]))} for k in rows[0]}

# (c) real Olist review text
rev = pd.read_csv(f"{FIELD}/olist_order_reviews_dataset.csv")
rev = rev.dropna(subset=["review_comment_message"])
rev = rev[rev["review_comment_message"].astype(str).str.strip().str.len() >= 10]
rev["label"] = np.where(rev["review_score"] >= 4, 1,
                        np.where(rev["review_score"] <= 2, 0, -1))
rev = rev[rev["label"] >= 0]
sample = rev.sample(n=2000, random_state=42)
texts = sample["review_comment_message"].astype(str).tolist()
label = sample["label"].tolist()
res["task_c_olist_real_text"] = {"n": len(texts),
                                 **bench(texts, label, seed=42)}
with open("results/probe4_text.json", "w") as fh:
    json.dump(res, fh, indent=2)
print(json.dumps(res, indent=2))

#!/usr/bin/env python3
"""validate.py -- deposit validator for the Paper 8 replication package.

Recomputes the manuscript's headline numbers from the deposited artifacts
themselves (run JSONs and full pair tables) and checks them against the
values printed in the submitted manuscript:

  * field families and promotions: Olist 107 pairs -> 13 promoted,
    Telco 231 -> 88, Bank Marketing 253 -> 33; no negative-control pair
    promoted on any of the three;
  * calibration audit: strength-only rule R0 averages 97.9 false
    discoveries per null world at n = 60; the conjoint rules R3/R4
    average 0.0;
  * family error certificate (Theorem 3.6) evaluated on the three field
    families: at most 4.2e-6 (largest deposited value 4.14e-6, Bank).

Stages:
  A  presence and parse of every deposited artifact class;
  B  headline recomputation from the pair tables and run JSONs;
  C  near-miss negatives (the split-half condition is load-bearing;
     planted control columns are never promoted; the Olist table's
     109 rows are 105 non-text pairs + both reductions of the 2 text
     target pairs the 107-pair family counts once each).

Exits 0 and prints OVERALL PASS only if every stage passes.
Standard library only.
"""
import csv
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent
FAILURES = []


def check(stage, name, ok, detail=""):
    print(f"[{stage}] {'PASS' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILURES.append(f"{stage}:{name}")


def load_json(rel):
    with open(ROOT / rel, encoding="utf-8") as fh:
        return json.load(fh)


def load_table(rel):
    with open(ROOT / rel, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def truthy(v):
    return str(v).strip().lower() == "true"


# ---------------------------------------------------------------- Stage A
REQUIRED = [
    "pipeline/facade_v2.py",
    "probes/probe1_calibration_vstar.py", "probes/probe1b_bias_smalln.py",
    "probes/probe2_fdr.py", "probes/probe2b_fdr_smalln.py",
    "probes/probe3_persistence.py", "probes/probe3b_splithalf.py",
    "probes/probe4_text.py", "probes/probe5_olist_run.py",
    "probes/probe5b_promoted_signed.py", "probes/probe6_phik_crosscheck.py",
    "probes/probe7_telco_run.py",
    "probes/results/probe5_full_pair_table.csv",
    "probes/results/probe7_full_pair_table.csv",
    "upgrade/theorem/probe8_certificate.py",
    "upgrade/cluster/probe10_cluster_certificate.py",
    "upgrade/fieldcert/probe11_field_certificates.py",
    "upgrade/thirddata/probe9_bank_run.py",
    "upgrade/thirddata/probe9_full_pair_table.csv",
    "upgrade/baseline/cpss_engine.py",
    "preregistrations/olist_preregistration.md",
    "preregistrations/telco_preregistration.md",
    "preregistrations/bank_preregistration.md",
    "preregistrations/cluster_preregistration.md",
]
missing = [r for r in REQUIRED if not (ROOT / r).is_file()]
check("A", "all promised artifact classes present", not missing,
      ("missing: " + ", ".join(missing)) if missing else f"{len(REQUIRED)} required paths")

json_paths = sorted(ROOT.glob("probes/results/*.json")) + \
    sorted(ROOT.glob("upgrade/*/*.json"))
bad_json = []
for p in json_paths:
    try:
        json.loads(p.read_text(encoding="utf-8"))
    except Exception as exc:  # noqa: BLE001
        bad_json.append(f"{p.relative_to(ROOT)}: {exc}")
check("A", "every deposited JSON parses", not bad_json,
      "; ".join(bad_json) if bad_json else f"{len(json_paths)} JSON files")

# ------------------------------------------------------------- Stage B
# Olist
p5 = load_json("probes/results/probe5_olist_run.json")
t5 = load_table("probes/results/probe5_full_pair_table.csv")
check("B", "Olist family size 107", p5["family_pairs"] == 107, f"json={p5['family_pairs']}")
check("B", "Olist BH rejections 69", p5["bh_rejections"] == 69, f"json={p5['bh_rejections']}")
check("B", "Olist promoted 13", len(p5["promoted_final"]) == 13,
      f"json list={len(p5['promoted_final'])}")
check("B", "Olist no control pair promoted", p5["control_pairs_promoted"] == 0)
t5_pairs = {(r["a"], r["b"]) for r in t5}
prom5 = {(d["a"], d["b"]) for d in p5["promoted_final"]}
check("B", "Olist promoted pairs all in deposited pair table", prom5 <= t5_pairs,
      f"{len(prom5)} promoted, table {len(t5_pairs)} unique pairs")

# Telco (pair table carries bh/floor columns -> recompute from the table)
p7 = load_json("probes/results/probe7_telco_run.json")
t7 = load_table("probes/results/probe7_full_pair_table.csv")
bh_floor7 = {(r["a"], r["b"]) for r in t7 if truthy(r["bh"]) and truthy(r["floor"])}
check("B", "Telco family size 231", p7["family_pairs"] == 231 and len(t7) == 231,
      f"json={p7['family_pairs']}, table rows={len(t7)}")
check("B", "Telco BH-and-floor recomputed from table = 89", len(bh_floor7) == 89,
      f"recomputed={len(bh_floor7)}, json bh_and_floor={p7['bh_and_floor']}")
failed7 = {(d["a"], d["b"]) for d in p7["promoted_but_splithalf_failed"]}
prom7 = {(d["a"], d["b"]) for d in p7["promoted_final"]}
check("B", "Telco promoted 88 = BH&floor minus split-half failures",
      prom7 == bh_floor7 - failed7 and len(prom7) == 88,
      f"promoted={len(prom7)}, failed={len(failed7)}")
check("B", "Telco no control pair promoted", p7["control_pairs_promoted"] == 0)

# Bank Marketing
p9 = load_json("upgrade/thirddata/probe9_bank_run.json")
t9 = load_table("upgrade/thirddata/probe9_full_pair_table.csv")
bh_floor9 = {(r["a"], r["b"]) for r in t9 if truthy(r["bh"]) and truthy(r["floor"])}
check("B", "Bank family size 253", p9["family_pairs"] == 253 and len(t9) == 253,
      f"json={p9['family_pairs']}, table rows={len(t9)}")
check("B", "Bank BH-and-floor recomputed from table = 34", len(bh_floor9) == 34,
      f"recomputed={len(bh_floor9)}, json bh_and_floor={p9['bh_and_floor']}")
check("B", "Bank promoted 33", p9["promoted_final_count"] == 33 and len(p9["promoted_final"]) == 33)
check("B", "Bank no control pair promoted (of 43 controls)",
      p9["control_pairs_promoted"] == 0 and p9["control_pairs_total"] == 43
      and p9["control_pairs_bh_and_floor"] == 0)

# Calibration audit
p2b = load_json("probes/results/probe2b_fdr_smalln.json")
n60 = p2b["null_n=60"]
check("B", "strength-only rule R0 mean false discoveries at n=60 is 97.9",
      n60["R0"]["mean"] == 97.9, f"json={n60['R0']['mean']}")
check("B", "conjoint rules R3/R4 mean false discoveries at n=60 are 0.0",
      n60["R3"]["mean"] == 0.0 and n60["R4"]["mean"] == 0.0)

# Family error certificate on the field families (manuscript: at most 4.2e-6)
fc = load_json("upgrade/fieldcert/field_certificates.json")
bounds = {k: v["C_U_full_family"] for k, v in fc["datasets"].items()}
check("B", "certificate bound at most 4.2e-6 on every field family",
      all(v <= 4.2e-6 for v in bounds.values()),
      ", ".join(f"{k}={v:.3g}" for k, v in bounds.items()))
check("B", "certificate file's observed promotions match 13/88/33, controls 0",
      [fc["datasets"][k]["promoted_pairs_observed"] for k in ("olist", "telco", "bank")] == [13, 88, 33]
      and all(fc["datasets"][k]["control_pairs_promoted_observed"] == 0 for k in bounds))
p8v = load_json("upgrade/theorem/probe8_validation.json")
check("B", "certificate engine agrees with the facade to <1e-6",
      p8v["max_abs_diff_strength_vs_facade"] < 1e-6 and p8v["engine_max_abs_diff_strength"] < 1e-6)
ev = load_json("upgrade/baseline/engine_validation.json")
check("B", "CPSS comparison engine validation gate passed", ev["gate_passed"] is True)

# ------------------------------------------------------------- Stage C
check("C", "split-half condition is load-bearing (Telco 89 -> 88, Bank 34 -> 33)",
      len(bh_floor7) - len(prom7) == 1 and p9["bh_and_floor"] - p9["promoted_final_count"] == 1)
rand_promoted = [o for o in p9["control_pair_outcomes"]
                 if (o["a"] == "rand_num" or o["b"] == "rand_num") and o.get("bh") and o.get("floor")]
check("C", "planted random control column never passes BH+floor on Bank",
      not rand_promoted and p9["control_pairs_bh_and_floor"] == 0)
text_rows = [r for r in t5 if str(r["text_pair"]).strip().lower() == "true"]
check("C", "Olist table = 105 non-text pairs + 4 text rows (2 text targets x 2 reductions); "
           "the 107 family counts each text target once",
      len(t5) == 109 and len(text_rows) == 4
      and fc["datasets"]["olist"]["text_pairs_in_family"] == 2,
      f"rows={len(t5)}, text rows={len(text_rows)}, family text pairs=2")

# ---------------------------------------------------------------- done
print()
if FAILURES:
    print(f"OVERALL FAIL  ({len(FAILURES)} check(s): {', '.join(FAILURES)})")
    sys.exit(1)
print("OVERALL PASS")

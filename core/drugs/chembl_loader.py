import json
import os
import random
import time
import requests
from core.drugs.drug_model import DrugModel

CH_EMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
CACHE_DIR = "data/drugs"
SEED_DIR = "data/seeds"


# ============================
# SEED LOADER
# ============================

def load_seed_drugs(disease):
    """
    Load disease-specific seed drugs from:
    data/seeds/<disease>_drugs.json
    """
    path = os.path.join(SEED_DIR, f"{disease}_drugs.json")

    try:
        with open(path, "r") as f:
            data = json.load(f)
        return data.get("drugs", [])
    except Exception:
        print(f"⚠️ No seed file found at {path}")
        return []


# ============================
# ChEMBL HELPERS
# ============================

def _get_target_id(gene_name):
    url = f"{CH_EMBL_API}/target/search.json?q={gene_name}"
    try:
        res = requests.get(url, timeout=10).json()
        targets = res.get("targets", [])
        if not targets:
            return None
        return targets[0]["target_chembl_id"]
    except Exception:
        return None


def _activity_to_confidence(activity):
    """
    Convert IC50/Ki/Kd values into a 0–1 confidence score
    Lower value = stronger binding
    """
    try:
        value = float(activity.get("standard_value", 0))
        if value <= 0:
            return None

        score = 1.0 / (1.0 + (value / 1000.0))
        return round(min(max(score, 0.05), 1.0), 3)
    except Exception:
        return None


def _fetch_activity_page(target_id, limit=50, offset=0, retries=3):
    url = (
        f"{CH_EMBL_API}/activity.json"
        f"?target_chembl_id={target_id}"
        f"&limit={limit}"
        f"&offset={offset}"
    )

    for attempt in range(retries):
        try:
            res = requests.get(url, timeout=10)
            if res.status_code == 200:
                return res.json().get("activities", [])
        except Exception:
            time.sleep(1.5 * (attempt + 1))

    return []


# ============================
# MAIN DISCOVERY FUNCTION
# ============================

def fetch_drugs_for_targets(
    targets,
    disease="ms",
    max_drugs=50,
    use_cache=True,
    mode="therapeutic"  # "binding" or "therapeutic"
):
    """
    Fetch drug panel for target proteins

    Cache:
      data/drugs/<disease>_panel.json

    Modes:
      binding     → only drugs with direct IC50/Ki/Kd data
      therapeutic → seed with known therapies + expand via ChEMBL
    """
    os.makedirs(CACHE_DIR, exist_ok=True)
    cache_file = os.path.join(CACHE_DIR, f"{disease}_panel.json")

    # ----------------------------
    # Load cache
    # ----------------------------
    if use_cache and os.path.exists(cache_file):
        print(f"📦 Loading drug panel from cache ({disease})")
        with open(cache_file, "r") as f:
            data = json.load(f)
        return [DrugModel(**d) for d in data]

    print(f"🌐 Fetching drug panel from ChEMBL ({mode} mode, disease={disease})...")
    drugs = {}

    # ----------------------------
    # Seed with known therapies
    # ----------------------------
    if mode == "therapeutic":
        seeds = load_seed_drugs(disease)
        print(f"🧬 Loaded {len(seeds)} seed drugs for {disease}")
        for name in seeds:
            drugs[name] = {}

    # ----------------------------
    # Expand using bioactivity
    # ----------------------------
    for gene in targets:
        target_id = _get_target_id(gene)
        if not target_id:
            print(f"⚠️ No ChEMBL target for {gene}")
            continue

        print(f"🔍 Target {gene} → {target_id}")

        offset = 0
        page_size = 50

        while len(drugs) < max_drugs:
            activities = _fetch_activity_page(
                target_id,
                limit=page_size,
                offset=offset
            )

            if not activities:
                break

            for act in activities:
                name = act.get("molecule_pref_name")
                if not name:
                    continue

                if mode == "therapeutic":
                    # Use available activity data if present, otherwise assign weak baseline
                    confidence = _activity_to_confidence(act) or round(0.1 + 0.2 * random.random(), 3)
                else:
                    confidence = _activity_to_confidence(act)
                    if confidence is None:
                        continue

                if name not in drugs:
                    drugs[name] = {}

                drugs[name][gene] = confidence

                if len(drugs) >= max_drugs:
                    break

            offset += page_size
            time.sleep(0.2)

        print(f"  → Collected {len(drugs)} drugs so far")

        if len(drugs) >= max_drugs:
            break

    # ----------------------------
    # Build Drug Models
    # ----------------------------
    drug_models = [
        DrugModel(name=name, targets=targets_map)
        for name, targets_map in drugs.items()
    ]

    # ----------------------------
    # Cache Results
    # ----------------------------
    with open(cache_file, "w") as f:
        json.dump(
            [{"name": d.name, "targets": d.targets} for d in drug_models],
            f,
            indent=2
        )

    print(f"✅ Cached {len(drug_models)} drugs to {cache_file}")
    return drug_models

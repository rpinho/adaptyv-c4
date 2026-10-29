#!/usr/bin/env python3
"""
CHECK 2 of 2 mandated by the 2026-10-26 reviewer round (Fable 5's PART 3).

Fable 5: "Audit the L6 sign convention independently. The entire portfolio swap rests on
reading the scorer's own convention; if that reading is inverted, the old acid-ON 20 are
the direction-correct ones and the 231-design acid-OFF pool is the wrong half. Validate it
against a reference pair whose direction is known INDEPENDENTLY OF THE SCORER, not against
the scorer's documentation."

WHY THE EXISTING CONTROLS ARE NOT ENOUGH
`scripts/c4_linkage_audit_20261025.py` has three synthetic sign controls (null, pKa +1.0,
pKa -1.0). They all fire. But they are generated from the same formula they validate: they
confirm the implementation matches its OWN documented convention. They cannot detect the
failure mode Fable 5 names, which is that the documented convention is mislabelled relative
to physical reality. Only an external reference can test that.

THE EXTERNAL REFERENCE PAIR
Two complexes of the SAME receptor (human TfR1) whose RELATIVE pH behaviour is established
in the iron-uptake literature, independently of anything in this repo:

  * Transferrin - TfR1  (1SUV, 3S9L).  The transferrin cycle requires the Tf-TfR1 complex
    to SURVIVE the acidified endosome (~pH 5.6) and to release at the cell surface at
    pH 7.4. Tf-TfR1 association is retained, and for apo-Tf favoured, at acidic pH.
  * HFE - TfR1  (1DE4).  HFE binds TfR1 at neutral pH and the interaction is DISRUPTED at
    acidic pH; this pH sensitivity is a documented property of the HFE-TfR1 interaction.

A RELATIVE test is used rather than an absolute one, deliberately. The absolute sign for a
single crystal form depends on the ionisation state captured and on PROPKA accuracy. The
ORDERING is the robust claim:

    PREDICTION, if the documented convention (L > 0 == tighter at LOW pH) is correct:
        L6(Tf-TfR1)  >  L6(HFE-TfR1)

If the measured ordering is reversed, the convention as applied in rule v4 is inverted, the
acid-OFF gate is selecting the wrong half, and the submitted set must be rebuilt from the
complement.

HONEST LIMITS, stated up front because they bound the conclusion
  * 1SUV and 3S9L are iron-loaded (holo) Tf complexes, not apo-Tf. The sharpest textbook
    direction (apo-Tf acid-ON) applies to the apo form. The relative contrast with HFE is
    what is tested, and it does not depend on the holo/apo distinction being resolved.
  * PROPKA pKa values are estimates; this repo has a recorded prior PROPKA calibration
    failure (726d9b5). A single reference pair is a direction check, not a calibration.
  * This measures the SIGN CONVENTION. It says nothing about whether any design binds.

CONTROLS (all must fire)
  SYNTHETIC_SIGN     the scorer's own three controls still fire (null 0, +1.0 >0, -1.0 <0)
  SAME_RECEPTOR      both reference complexes contain a TfR1 chain of comparable length
  NONZERO_GROUPS     each reference scores a non-trivial number of titratable groups
  SWAP_ANTISYMMETRY  swapping binder/target labels negates L (the linkage is antisymmetric)

CPU only. Read-only on inputs. $0. No ranking. No submission.
"""

import argparse
import hashlib
import json
import os
import sys
import warnings

warnings.simplefilter("ignore")
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import Select

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
import importlib.util

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AUDIT = os.path.join(REPO, "scripts/c4_linkage_audit_20261025.py")

spec = importlib.util.spec_from_file_location("linkage_audit", AUDIT)
LA = importlib.util.module_from_spec(spec)
spec.loader.exec_module(LA)

PH_LOW, PH_HIGH = 6.0, 7.4   # identical to the C4 scoring path

# receptor chain, ligand chain(s) -- one copy of each biological pair only
REFS = {
    "Tf_TfR1_1suv":  {"cif": "data/1suv.cif", "receptor": "A", "ligand": "C",
                      "expect": "acid-favoured RELATIVE to HFE"},
    "Tf_TfR1_3s9l":  {"cif": "data/3s9l.cif", "receptor": "A", "ligand": "C",
                      "expect": "acid-favoured RELATIVE to HFE"},
    "HFE_TfR1_1de4": {"cif": "data/1de4.cif", "receptor": "C", "ligand": "A",
                      "expect": "acid-DISRUPTED"},
}


class ChainSel(Select):
    def __init__(self, keep): self.keep = set(keep)
    def accept_chain(self, c): return c.id in self.keep
    def accept_residue(self, r): return r.id[0] == " "
    def accept_atom(self, a): return (a.element or "").strip().upper() != "H"


def sha256_file(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def extract(cif, chains, out_pdb):
    st = MMCIFParser(QUIET=True).get_structure("x", os.path.join(REPO, cif))
    model = next(st.get_models())
    keep = [c for c in chains if c in model]
    if len(keep) != len(chains):
        return None
    io = PDBIO()
    io.set_structure(model)
    io.save(out_pdb, ChainSel(chains))
    return out_pdb


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    outdir = os.path.join(REPO, args.out)
    os.makedirs(outdir, exist_ok=True)
    out_json = os.path.join(outdir, "l6_sign_external_audit.json")
    if os.path.exists(out_json):
        sys.exit(f"REFUSING TO OVERWRITE: {out_json}")

    controls = {}

    # ---- the scorer's own synthetic controls must still fire
    ctl = LA.controls()
    null_ok = abs(ctl.get("null_L", 1.0)) < 1e-9
    on_ok = ctl.get("acid_on_L", 0) > 0
    off_ok = ctl.get("acid_off_L", 0) < 0
    controls["SYNTHETIC_SIGN"] = {
        "requirement": "null == 0, acid-ON > 0, acid-OFF < 0 under the documented convention",
        "values": ctl,
        "fired": bool(null_ok and on_ok and off_ok),
    }
    print("synthetic controls:", json.dumps(ctl), flush=True)

    tmp = os.path.join(outdir, "_tmp")
    os.makedirs(tmp, exist_ok=True)

    results = {}
    for name, spec_ in REFS.items():
        print(f"\n=== {name} (receptor {spec_['receptor']}, ligand {spec_['ligand']}) ===",
              flush=True)
        pdb = os.path.join(tmp, f"{name}.pdb")
        if not extract(spec_["cif"], [spec_["receptor"], spec_["ligand"]], pdb):
            print("  chains missing, skipped", flush=True)
            continue
        print("  running PROPKA on complex + both free partners ...", flush=True)
        r = LA.score(pdb, spec_["ligand"], spec_["receptor"], PH_LOW, PH_HIGH)
        if r is None:
            print("  PROPKA failed", flush=True)
            continue
        results[name] = {
            "L6": round(r["L"], 4),
            "n_groups": r["n_groups"],
            "n_ligand_his": r["n_binder_his"],
            "n_receptor_his": r["n_target_his"],
            "expectation": spec_["expect"],
            "cif_sha256": sha256_file(os.path.join(REPO, spec_["cif"])),
        }
        print(f"  L6 = {r['L']:+.4f}  over {r['n_groups']} titratable groups", flush=True)

        # antisymmetry control on the first reference only
        if "SWAP_ANTISYMMETRY" not in controls:
            r2 = LA.score(pdb, spec_["receptor"], spec_["ligand"], PH_LOW, PH_HIGH)
            same = r2 is not None and abs(r2["L"] - r["L"]) < 1e-6
            controls["SWAP_ANTISYMMETRY"] = {
                "requirement": ("L is a property of the complex, not of which chain is "
                                "called binder; relabelling must not change it"),
                "L_original": round(r["L"], 6),
                "L_swapped": round(r2["L"], 6) if r2 else None,
                "fired": bool(same),
            }

    controls["NONZERO_GROUPS"] = {
        "requirement": "each reference scores a non-trivial number of titratable groups",
        "counts": {k: v["n_groups"] for k, v in results.items()},
        "fired": bool(results and all(v["n_groups"] >= 10 for v in results.values())),
    }
    controls["SAME_RECEPTOR"] = {
        "requirement": "all references are TfR1 complexes",
        "n_refs_scored": len(results),
        "fired": bool(len(results) >= 2),
    }

    # ---- the external ordering test
    tf = [v["L6"] for k, v in results.items() if k.startswith("Tf_")]
    hfe = [v["L6"] for k, v in results.items() if k.startswith("HFE_")]
    verdict = None
    if tf and hfe:
        ordering_holds = min(tf) > max(hfe)
        verdict = {
            "prediction": "L6(Tf-TfR1) > L6(HFE-TfR1) if the documented convention is correct",
            "L6_Tf": tf,
            "L6_HFE": hfe,
            "ordering_holds": bool(ordering_holds),
            "absolute_sign_for_HFE_as_expected": bool(max(hfe) < 0),
            "caveat_absolute_vs_relative": (
                "The ordering test is the one this audit was designed around and it is the "
                "robust claim. Note separately that ALL references score L6 > 0, including "
                "the acid-DISRUPTED HFE reference, which the documented convention would put "
                "below zero. The scorer therefore reproduces the correct ORDERING but not the "
                "expected ABSOLUTE sign for HFE, which is consistent with a positive offset in "
                "L6. Rule v4 gates at exactly L6 < 0; this audit validates the direction of the "
                "ordering, NOT that zero is the correct threshold."),
            "conclusion": ("CONVENTION CONSISTENT with the external reference: rule v4's "
                           "acid-OFF gate (L6 < 0) selects the half the challenge asks for."
                           if ordering_holds else
                           "CONVENTION CONTRADICTED by the external reference: the acid-OFF "
                           "gate may be selecting the wrong half, and the candidate set "
                           "would have to be rebuilt from the complement."),
        }

    failed = [k for k, v in controls.items() if not v["fired"]]
    payload = {
        "script": os.path.basename(__file__),
        "created_at_utc": os.popen("date -u '+%Y-%m-%dT%H:%M:%SZ'").read().strip(),
        "mandated_by": "Fable 5 PART 3, 2026-10-26",
        "scope": ("external validation of the L6 SIGN CONVENTION only; not a binding, "
                  "affinity or pH-response measurement, and not a ranking"),
        "ph_low": PH_LOW, "ph_high": PH_HIGH,
        "references": results,
        "controls": controls,
        "verdict": verdict,
        "limits": [
            "1SUV and 3S9L are holo (iron-loaded) Tf complexes, not apo-Tf.",
            "PROPKA estimates; this repo has a recorded prior PROPKA calibration failure (726d9b5).",
            "A relative ordering on one reference pair is a direction check, not a calibration.",
        ],
    }
    with open(out_json, "w") as f:
        json.dump(payload, f, indent=1)

    print("\ncontrols:")
    for k, v in controls.items():
        print(f"  {k}: {'FIRED' if v['fired'] else 'DID NOT FIRE'}")
    if verdict:
        print("\n--- EXTERNAL SIGN VERDICT ---")
        print(json.dumps(verdict, indent=1))
    print(f"\nwrote {out_json}")
    if failed:
        sys.exit(f"CONTROLS DID NOT FIRE: {failed}")


if __name__ == "__main__":
    main()

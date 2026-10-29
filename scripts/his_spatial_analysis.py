"""3D spatial analysis of His residues in C4 v8 ESMFold structures.

Checks:
- His-His Cα distances (actimo-labs criterion: >8Å apart for independent switching)
- His-D/E Cα distances (intramolecular pairing risk: <6Å means likely salt bridge)
- His solvent accessibility proxy (Cα distance to centroid — crude but fast)

Usage:
    python3 -I scripts/his_spatial_analysis.py
"""
import csv
import math
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / "results"
CSV = RESULTS / "submission_v8_20261023.csv"
PDB_DIR = RESULTS / "c4_v8_pdbs"


def parse_ca_coords(pdb_path):
    coords = {}
    for line in pdb_path.read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            res_num = int(line[22:26].strip())
            res_name = line[17:20].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords[res_num] = {"name": res_name, "xyz": (x, y, z)}
    return coords


def dist(a, b):
    return math.sqrt(sum((ai - bi)**2 for ai, bi in zip(a, b)))


def centroid(coords):
    xs = [c["xyz"][0] for c in coords.values()]
    ys = [c["xyz"][1] for c in coords.values()]
    zs = [c["xyz"][2] for c in coords.values()]
    n = len(coords)
    return (sum(xs)/n, sum(ys)/n, sum(zs)/n)


def analyze_pdb(name, pdb_path):
    coords = parse_ca_coords(pdb_path)
    if not coords:
        return {"name": name, "error": "no CA atoms"}

    cent = centroid(coords)

    his_res = {k: v for k, v in coords.items() if v["name"] == "HIS"}
    de_res = {k: v for k, v in coords.items() if v["name"] in ("ASP", "GLU")}

    if not his_res:
        return {"name": name, "his_count": 0}

    his_nums = sorted(his_res.keys())

    # His-His distances
    hh_pairs = []
    close_hh = 0
    for i, h1 in enumerate(his_nums):
        for h2 in his_nums[i+1:]:
            d = dist(his_res[h1]["xyz"], his_res[h2]["xyz"])
            hh_pairs.append((h1, h2, round(d, 1)))
            if d < 8.0:
                close_hh += 1

    # His-D/E distances
    hde_close = []
    for h in his_nums:
        for d_num, d_info in de_res.items():
            d_val = dist(his_res[h]["xyz"], d_info["xyz"])
            if d_val < 6.0:
                hde_close.append((h, d_num, d_info["name"], round(d_val, 1)))

    # His distance from centroid (surface proxy)
    max_dist = max(dist(c["xyz"], cent) for c in coords.values())
    his_surface = []
    for h in his_nums:
        d_cent = dist(his_res[h]["xyz"], cent)
        rel = d_cent / max_dist if max_dist > 0 else 0
        his_surface.append((h, round(d_cent, 1), round(rel, 2)))

    surface_his = sum(1 for _, _, r in his_surface if r > 0.6)
    core_his = sum(1 for _, _, r in his_surface if r < 0.4)

    # Min His-His distance
    min_hh = min(d for _, _, d in hh_pairs) if hh_pairs else 999

    # His with D/E within 6Å
    his_with_de = len(set(h for h, _, _, _ in hde_close))

    return {
        "name": name,
        "his_count": len(his_res),
        "de_count": len(de_res),
        "close_hh_pairs": close_hh,
        "min_hh_dist": min_hh,
        "his_with_de_6A": his_with_de,
        "surface_his": surface_his,
        "core_his": core_his,
        "hh_pairs": hh_pairs,
        "hde_close": hde_close,
        "his_surface": his_surface,
    }


def main():
    designs = {}
    with open(CSV) as f:
        for r in csv.DictReader(f):
            designs[r["name"]] = r["sequence"]

    pdb_files = sorted(PDB_DIR.glob("*.pdb"))
    print(f"3D His spatial analysis for {len(pdb_files)} C4 v8 PDB files")
    print("=" * 120)

    all_results = []
    for pdb_path in pdb_files:
        name = pdb_path.stem
        r = analyze_pdb(name, pdb_path)
        all_results.append(r)

    # Summary table
    print(f"\n{'Design':<40s} {'His':>4s} {'D/E':>4s} {'HH<8Å':>6s} {'minHH':>6s} "
          f"{'H-DE<6Å':>8s} {'Surf':>5s} {'Core':>5s}")
    print("-" * 120)

    for r in sorted(all_results, key=lambda x: x.get("close_hh_pairs", 0)):
        if "error" in r or r.get("his_count", 0) == 0:
            continue
        print(f"{r['name']:<40s} {r['his_count']:>4d} {r['de_count']:>4d} "
              f"{r['close_hh_pairs']:>6d} {r['min_hh_dist']:>6.1f} "
              f"{r['his_with_de_6A']:>8d} {r['surface_his']:>5d} {r['core_his']:>5d}")

    # Detailed report for concerning designs
    print("\n\nHIS-HIS CLOSE PAIRS (<8Å) — actimo-labs criterion: prefer >8Å separation")
    print("=" * 120)
    for r in sorted(all_results, key=lambda x: -x.get("close_hh_pairs", 0)):
        if r.get("close_hh_pairs", 0) == 0:
            continue
        print(f"\n{r['name']} ({r['close_hh_pairs']} close pairs):")
        for h1, h2, d in r["hh_pairs"]:
            if d < 8.0:
                label = "CLOSE" if d < 5.0 else "NEAR"
                print(f"  H{h1}-H{h2}: {d:.1f}Å [{label}]")

    print("\n\nHIS-D/E CONTACTS (<6Å) — intramolecular salt bridge risk")
    print("=" * 120)
    for r in sorted(all_results, key=lambda x: -x.get("his_with_de_6A", 0)):
        if r.get("his_with_de_6A", 0) == 0:
            continue
        print(f"\n{r['name']} ({r['his_with_de_6A']}/{r['his_count']} His near D/E):")
        for h, d_num, d_name, d_val in r["hde_close"]:
            print(f"  H{h}-{d_name}{d_num}: {d_val:.1f}Å")

    # Quality score
    print("\n\nCOMPOSITE QUALITY RANKING")
    print("=" * 120)
    print("Score = switching_his (from PROPKA) - close_hh_pairs - his_with_de*0.5")
    print("Higher = better pH-switch potential with less intramolecular interference")


if __name__ == "__main__":
    main()

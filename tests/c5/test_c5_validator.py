#!/usr/bin/env python3
"""Fail-closed tests for scripts/c5_validate_submission.py. No GPU, no network, no sequences."""
import csv, json, subprocess, sys, tempfile, pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]
V = ROOT / "scripts" / "c5_validate_submission.py"
PY3 = sys.executable
AA = "ACDEFGHIKLMNPQRSTVWY"
results = []


def run(rows, *extra, cols=("name", "sequence", "molecule_class")):
    with tempfile.TemporaryDirectory() as d:
        p = pathlib.Path(d) / "sub.csv"
        with open(p, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(cols)
            w.writerows(rows)
        r = subprocess.run([PY3, str(V), str(p), *extra], capture_output=True, text=True)
        return r.returncode, r.stdout + r.stderr


def seq(n, ch="A"):
    return ch * n


def uniq(n, i, ch="A"):
    """n residues, unique per row i (so uniqueness checks test the intended rule)."""
    a, b = AA[i % 20], AA[(i // 20) % 20]
    return (a + b + AA[(i + 7) % 20]) * 2 + ch * (n - 6)


def case(name, rows, expect_rc, *extra, cols=("name", "sequence", "molecule_class"), needle=None):
    rc, out = run(rows, *extra, cols=cols)
    ok = rc == expect_rc and (needle is None or needle in out)
    results.append(ok)
    print(("PASS " if ok else "FAIL ") + name + ("" if ok else f"  rc={rc} want={expect_rc}\n{out[:400]}"))


N = 20
# --- happy path: minimum legal length 10, official default class ---
case("official single_chain class + 10aa minimum passes",
     [[f"d{i}", uniq(10, i), "single_chain"] for i in range(N)], 0)

# --- official enum: single_chain is the demonstrated upload value; protein is NOT default ---
case("single_chain is ACCEPTED by default (official portal/template value)",
     [[f"d{i}", uniq(60, i), "single_chain"] for i in range(N)], 0)
case("protein is REJECTED by default (public-FAQ value contradicted by the form)",
     [[f"d{i}", uniq(60, i), "protein"] for i in range(N)], 1, needle="molecule_class")
case("protein remains available only via explicit --classes override",
     [[f"d{i}", uniq(60, i), "protein"] for i in range(N)], 0,
     "--classes", "protein")
case("9aa is below the official 10aa floor",
     [[f"d{i}", seq(9), "single_chain"] for i in range(N)], 1, needle="outside")
case("250aa is allowed (inclusive upper bound)",
     [[f"d{i}", uniq(250, i), "single_chain"] for i in range(N)], 0)
case("251aa is rejected",
     [[f"d{i}", uniq(251, i), "single_chain"] for i in range(N)], 1, needle="outside")

# --- slot counts ---
case("19 rows fails the default 20 slots",
     [[f"d{i}", uniq(50, i), "single_chain"] for i in range(19)], 1, needle="row count")
case("40 rows fails 20 slots but passes --slots 40",
     [[f"d{i}", uniq(50, i), "single_chain"] for i in range(40)], 1, needle="row count")
case("40 rows passes --slots 40",
     [[f"d{i}", uniq(50, i), "single_chain"] for i in range(40)], 0, "--slots", "40")

# --- uniqueness / alphabet ---
case("duplicate names rejected",
     [[f"d{i}", seq(50), "single_chain"] for i in range(N - 1)] + [["d0", seq(51), "single_chain"]],
     1, needle="duplicate names")
case("duplicate sequences rejected",
     [["a", uniq(50, 0), "single_chain"], ["b", uniq(50, 0), "single_chain"]] +
     [[f"d{i}", uniq(50 + i, i), "single_chain"] for i in range(N - 2)], 1, needle="duplicate sequences")
case("non-canonical residue X rejected",
     [["d0", "X" + seq(49), "single_chain"]] + [[f"d{i}", seq(50 + i), "single_chain"] for i in range(N - 1)],
     1, needle="non-canonical")
case("lowercase sequence rejected",
     [["d0", uniq(50, 0).lower(), "single_chain"]] + [[f"d{i}", uniq(50 + i, i), "single_chain"] for i in range(N - 1)],
     1, needle="non-canonical")
case("empty sequence rejected",
     [["d0", "", "single_chain"]] + [[f"d{i}", uniq(50 + i, i), "single_chain"] for i in range(N - 1)],
     1, needle="empty sequence")
case("name with whitespace rejected",
     [["d 0", uniq(50, 0), "single_chain"]] + [[f"d{i}", uniq(50 + i, i), "single_chain"] for i in range(N - 1)],
     1, needle="invalid name")

# --- columns ---
case("missing molecule_class column rejected",
     [[f"d{i}", seq(50)] for i in range(N)], 1,
     cols=("name", "sequence"), needle="missing required column")
case("custom --columns honoured",
     [[f"d{i}", uniq(50, i)] for i in range(N)], 0,
     "--columns", "name,sequence", "--classes", "", cols=("name", "sequence"))
case("omitting molecule_class while its class filter is active fails once, not per row",
     [[f"d{i}", uniq(50, i)] for i in range(N)], 1,
     cols=("name", "sequence"), needle="missing required column 'molecule_class'")

# --- fabs ---
vh, vl = seq(120), seq(110)
case("fab_kappa with VH:VL passes without --vh-vl",
     [[f"d{i}", f"{vh}:{vl}{AA[i % 20]}", "fab_kappa"] for i in range(N)], 0)
case("fab_lambda with | separator passes",
     [[f"d{i}", f"{vh}|{vl}{AA[i % 20]}", "fab_lambda"] for i in range(N)], 0)
case("fab_kappa with no separator rejected",
     [[f"d{i}", vh + AA[i % 20], "fab_kappa"] for i in range(N)], 1, needle="separator")
case("fab with a too-short VH rejected",
     [[f"d{i}", f"{seq(50)}:{vl}", "fab_kappa"] for i in range(N)], 1, needle="VH length")
case("plain single_chain with a colon is not split by default",
     [[f"d{i}", f"{seq(25)}:{seq(25)}", "single_chain"] for i in range(N)], 1, needle="non-canonical")

# --- known-binder guard ---
with tempfile.TemporaryDirectory() as d:
    fa = pathlib.Path(d) / "known.fa"
    fa.write_text(">known1\n" + seq(50, "W") + "\n")
    rows = [["d0", seq(50, "W"), "single_chain"]] + [[f"d{i}", seq(50 + i), "single_chain"] for i in range(N - 1)]
    rc, out = run(rows, "--forbid-fasta", str(fa))
    ok = rc == 1 and "known" in out
    results.append(ok)
    print(("PASS " if ok else "FAIL ") + "verbatim known-binder sequence rejected")
    rows2 = [[f"d{i}", seq(50 + i), "single_chain"] for i in range(N)]
    rc, out = run(rows2, "--forbid-fasta", str(fa))
    ok = rc == 0
    results.append(ok)
    print(("PASS " if ok else "FAIL ") + "disjoint set passes the known-binder guard")

# --- manifest ---
with tempfile.TemporaryDirectory() as d:
    p = pathlib.Path(d) / "sub.csv"
    man = pathlib.Path(d) / "man.json"
    with open(p, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["name", "sequence", "molecule_class"])
        w.writerows([[f"d{i}", seq(50 + i), "single_chain"] for i in range(N)])
    r = subprocess.run([PY3, str(V), str(p), "--manifest", str(man)], capture_output=True, text=True)
    m = json.loads(man.read_text())
    ok = r.returncode == 0 and m["pass"] is True and m["reviewed"] is False and len(m["sha256"]) == 64
    results.append(ok)
    print(("PASS " if ok else "FAIL ") + "manifest records pass=false review state + sha256")
    r = subprocess.run([PY3, str(V), str(p), "--manifest", str(man), "--reviewer", "Ricardo"],
                       capture_output=True, text=True)
    m = json.loads(man.read_text())
    ok = r.returncode == 0 and m["reviewed"] is True and m["reviewer"] == "Ricardo"
    results.append(ok)
    print(("PASS " if ok else "FAIL ") + "reviewer sign-off recorded in manifest")

print(f"\n{sum(results)}/{len(results)} passed")
sys.exit(0 if all(results) else 1)

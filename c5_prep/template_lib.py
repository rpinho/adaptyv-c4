#!/usr/bin/env python3
"""C5 target-independent template tooling (Kimi plan §6, 2026-10-25).

CPU-only. No target-specific work: nothing here knows what C5's target is.
Everything operates on experimental PDB/mmCIF files supplied by the caller.

The rules encoded here are the ones Kimi listed:
  * parse experimental structures and extract chains with metadata
  * build a template manifest (PDB ID, chain, method, resolution, release date,
    source checksum, sequence coverage, missing-residue annotations)
  * reject predicted structures used as "templates" (fail closed)
  * CPU alignment/superposition and template ranking
  * explicit output schema + fail-closed behaviour
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import json
import math
import pathlib
import re

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}
AA = set("ACDEFGHIKLMNPQRSTVWYUO")

# File-name / header markers that identify a *predicted* model rather than an
# experimentally determined structure. Any hit is a hard failure.
PREDICTION_MARKERS = (
    "alphafold", "af-", "af_", "esmfold", "esm2", "boltz", "rfd", "rfdiffusion",
    "colabfold", "model_", "_predicted", ".pred.", "unrelaxed", "rank_",
)
# _exptl.method values that denote a computational model.
PREDICTED_METHODS = {"PREDICTED", "COMPUTATIONAL MODEL", "THEORETICAL MODEL",
                     "AB INITIO", "DEEP LEARNING"}

MANIFEST_SCHEMA = {
    "type": "object",
    "required": ["generated_at", "tool", "templates", "rejected", "snapshot"],
    "properties": {
        "generated_at": {"type": "string"},
        "tool": {"type": "object"},
        "snapshot": {"type": "object"},
        "templates": {"type": "array"},
        "rejected": {"type": "array"},
    },
}


class TemplateError(Exception):
    """Any condition that must fail closed."""


def sha256_file(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def guess_format(path: pathlib.Path) -> str:
    suffix = path.suffix.lower()
    if suffix in (".cif", ".mmcif"):
        return "cif"
    if suffix == ".pdb":
        return "pdb"
    if suffix in (".ent", ".pdb1"):
        return "pdb"
    raise TemplateError(f"{path.name}: unsupported extension (want .cif/.mmcif/.pdb)")


def _parse_pdb(path: pathlib.Path):
    atoms = []
    with open(path, errors="replace") as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            alt = line[16].strip()
            if alt not in ("", "A"):
                continue
            name = line[12:16].strip()
            resn = line[17:20].strip()
            chain = line[21]
            try:
                resseq = int(line[22:26])
            except ValueError:
                continue
            icode = line[26].strip()
            try:
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            b = line[60:66].strip()
            atoms.append({"name": name, "resn": resn, "chain": chain, "resseq": resseq,
                          "icode": icode, "xyz": xyz,
                          "bfac": float(b) if b else 0.0,
                          "group": "HETATM" if line.startswith("HETATM") else "ATOM"})
    return atoms


def _parse_cif(path: pathlib.Path):
    try:
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    except ImportError as exc:  # pragma: no cover - environment guard
        raise TemplateError("Biopython is required to parse mmCIF") from exc
    d = MMCIF2Dict(str(path))
    need = ["_atom_site.group_PDB", "_atom_site.label_atom_id",
            "_atom_site.label_comp_id", "_atom_site.Cartn_x", "_atom_site.Cartn_y",
            "_atom_site.Cartn_z"]
    missing = [k for k in need if k not in d]
    if missing:
        raise TemplateError(f"{path.name}: mmCIF lacks {missing}")
    # Grouping MUST use label_asym_id: author chain IDs may be reused across
    # entities (6X18 does exactly this), which would silently merge chains.
    asym = d.get("_atom_site.label_asym_id")
    if asym is None:
        raise TemplateError(f"{path.name}: mmCIF has no label_asym_id")
    auth = d.get("_atom_site.auth_asym_id", asym)
    seqid = d.get("_atom_site.auth_seq_id") or d.get("_atom_site.label_seq_id")
    alt = d.get("_atom_site.label_alt_id", ["."] * len(d["_atom_site.group_PDB"]))
    bfac = d.get("_atom_site.B_iso_or_equiv", ["0"] * len(d["_atom_site.group_PDB"]))
    atoms = []
    for i, group in enumerate(d["_atom_site.group_PDB"]):
        a = (alt[i] or ".").strip()
        if a not in (".", "?", "A"):
            continue
        raw_seq = str(seqid[i]).strip()
        m = re.match(r"-?\d+", raw_seq)
        resseq = int(m.group(0)) if m else None
        try:
            xyz = (float(d["_atom_site.Cartn_x"][i]), float(d["_atom_site.Cartn_y"][i]),
                   float(d["_atom_site.Cartn_z"][i]))
        except ValueError:
            continue
        try:
            b = float(bfac[i])
        except (TypeError, ValueError):
            b = 0.0
        atoms.append({"name": d["_atom_site.label_atom_id"][i], "resn": d["_atom_site.label_comp_id"][i],
                      "chain": asym[i], "auth_chain": auth[i], "resseq": resseq, "icode": "",
                      "xyz": xyz, "bfac": b, "group": group, "raw_seq": raw_seq})
    meta = {}
    for key in ("_exptl.method", "_struct.title", "_entry.id",
                "_em_3d_reconstruction.resolution", "_refine.ls_d_res_high",
                "_reflns.d_resolution_high", "_database_PDB_rev.date_original",
                "_pdbx_database_status.recvd_initial_deposition_date",
                "_struct_keywords.pdbx_keywords"):
        if key in d:
            v = d[key]
            meta[key] = v[0] if isinstance(v, list) else v
    return atoms, meta


def parse_structure(path: pathlib.Path):
    """Return (atoms, meta). Raises TemplateError on unreadable input."""
    path = pathlib.Path(path)
    if not path.is_file():
        raise TemplateError(f"{path}: not a file")
    fmt = guess_format(path)
    if fmt == "pdb":
        return _parse_pdb(path), {}
    return _parse_cif(path)


def chains(atoms):
    """Group POLYMER atoms by chain, preserving residue order.

    Non-polymer records (waters, glycans, ions, ligands) are excluded from the
    chain sequence and counted separately by `hetero_summary`. Returns
    {chain: [ {resseq, icode, resn, atoms:{name:xyz}, bfac:{name:float}} ]}
    """
    out, order = {}, {}
    for a in atoms:
        if a.get("group", "ATOM") != "ATOM":
            continue
        ch = a["chain"]
        key = (a["resseq"], a["icode"])
        bucket = out.setdefault(ch, {})
        if key not in bucket:
            bucket[key] = {"resseq": a["resseq"], "icode": a["icode"], "resn": a["resn"],
                           "atoms": {}, "bfac": {}}
            order.setdefault(ch, []).append(key)
        bucket[key]["atoms"].setdefault(a["name"], a["xyz"])
        bucket[key]["bfac"].setdefault(a["name"], a["bfac"])
    return {ch: [out[ch][k] for k in order[ch]] for ch in out}


def sequence(residues):
    return "".join(THREE_TO_ONE.get(r["resn"], "X") for r in residues)


def chain_sequence(residues):
    seq = sequence(residues)
    if set(seq) - AA:
        bad = sorted(set(seq) - AA)
        raise TemplateError(f"non-amino-acid residues in chain: {bad}")
    return seq


def coverage(residues, seq_one):
    """(fraction of polymer residues with complete N/CA/C/O, missing-residue runs).

    `seq_one` is the deposited one-letter sequence (entity level) when available,
    otherwise the observed chain sequence.
    """
    complete = sum(1 for r in residues if {"N", "CA", "C", "O"} <= set(r["atoms"]))
    observed_numbers = sorted(r["resseq"] for r in residues if r["resseq"] is not None)
    gaps = []
    if observed_numbers:
        for a, b in zip(observed_numbers, observed_numbers[1:]):
            if b - a > 1:
                gaps.append([a + 1, b - 1])
        lo, hi = observed_numbers[0], observed_numbers[-1]
        if len(observed_numbers) < (hi - lo + 1):
            missing = (hi - lo + 1) - len(observed_numbers)
        else:
            missing = 0
    else:
        missing = 0
    return {"residues_observed": len(residues), "residues_complete_backbone": complete,
            "complete_backbone_fraction": round(complete / len(residues), 4) if residues else 0.0,
            "missing_residue_runs": gaps, "missing_residue_count": missing,
            "observed_range": [observed_numbers[0], observed_numbers[-1]] if observed_numbers else []}


def reject_if_predicted(path: pathlib.Path, meta, atoms):
    """Fail closed when the file is a computational model, not an experimental one."""
    name = path.name.lower()
    hits = [m for m in PREDICTION_MARKERS if m in name]
    if hits:
        raise TemplateError(f"{path.name}: filename looks like a predicted model ({hits})")
    method = str(meta.get("_exptl.method", "")).strip().upper()
    if method in PREDICTED_METHODS:
        raise TemplateError(f"{path.name}: _exptl.method={method!r} is a computational model")
    if not method:
        # No method means we cannot prove it is experimental -> refuse.
        raise TemplateError(f"{path.name}: no _exptl.method; cannot prove experimental origin")
    # AlphaFold/ESMFold style metadata blocks.
    plddt_like = [a for a in atoms if a["name"] in ("pLDDT", "PLDDT")]
    if plddt_like:
        raise TemplateError(f"{path.name}: contains pLDDT atoms (predicted model)")
    return method


def _resolution(meta):
    for key in ("_em_3d_reconstruction.resolution", "_refine.ls_d_res_high",
                "_reflns.d_resolution_high"):
        raw = meta.get(key)
        if raw in (None, "", "?", "."):
            continue
        try:
            return float(str(raw).strip())
        except ValueError:
            continue
    return None


def _release_date(meta):
    for key in ("_database_PDB_rev.date_original",
                "_pdbx_database_status.recvd_initial_deposition_date"):
        raw = meta.get(key)
        if raw in (None, "", "?", "."):
            continue
        m = re.match(r"(\d{4})-(\d{2})-(\d{2})", str(raw).strip())
        if m:
            try:
                return _dt.date(*(int(x) for x in m.groups())).isoformat()
            except ValueError:
                continue
    return None


def describe_chain(path, chain_id, residues, meta, source_seq=None):
    try:
        seq = chain_sequence(residues)
    except TemplateError as exc:
        raise TemplateError(f"{pathlib.Path(path).name}:{chain_id}: {exc}") from exc
    entry = {
        "pdb_id": str(meta.get("_entry.id", pathlib.Path(path).stem)).strip().upper(),
        "chain": chain_id,
        "source_file": str(pathlib.Path(path).name),
        "source_sha256": sha256_file(pathlib.Path(path)),
        "method": str(meta.get("_exptl.method", "")).strip().upper(),
        "resolution_a": _resolution(meta),
        "release_date": _release_date(meta),
        "title": str(meta.get("_struct.title", "")).strip(),
        "length": len(seq),
        "sequence": seq,
        "sequence_sha256": hashlib.sha256(seq.encode()).hexdigest(),
        "coverage": coverage(residues, source_seq or seq),
    }
    return entry


def hetero_summary(atoms):
    """Count non-polymer records per chain (waters/ions/ligands are not templates)."""
    counts = {}
    for a in atoms:
        if a.get("group", "ATOM") == "ATOM":
            continue
        key = (a["chain"], a["resn"])
        counts[key] = counts.get(key, 0) + 1
    return {f"{ch}:{resn}": n for (ch, resn), n in sorted(counts.items())}


def build_manifest(paths, snapshot=None, tool_version="c5_prep/1.0.0",
                   fail_closed=True):
    """Parse each file, extract every polymer chain, reject predictions, and record.

    With ``fail_closed`` (the default, and the only mode template selection may
    use) ANY unreadable file, predicted input or chain that cannot be typed as a
    polymer raises TemplateError immediately instead of being collected. Pass
    ``fail_closed=False`` only when building an inventory of a heterogeneous
    directory that is allowed to contain unusable files.
    """
    templates, rejected = [], []
    for raw in paths:
        p = pathlib.Path(raw)
        try:
            atoms, meta = parse_structure(p)
            method = reject_if_predicted(p, meta, atoms)
            ch = chains(atoms)
            if not ch:
                raise TemplateError(f"{p.name}: no polymer atoms")
            hetero = hetero_summary(atoms)
            for cid, residues in sorted(ch.items()):
                if len(residues) < 3:
                    raise TemplateError(f"{p.name}:{cid}: only {len(residues)} polymer residues")
                entry = describe_chain(p, cid, residues, meta)
                entry["method_verified"] = method
                entry["hetero_records"] = {k: v for k, v in hetero.items()
                                           if k.startswith(f"{cid}:")}
                templates.append(entry)
        except TemplateError as exc:
            if fail_closed:
                raise
            rejected.append({"source_file": p.name, "chain": None, "reason": str(exc)})
    manifest = {
        "generated_at": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "tool": {"name": "c5_prep/template_lib.py", "version": tool_version, "cpu_only": True},
        "snapshot": snapshot or {"note": "no template-database snapshot recorded"},
        "templates": templates,
        "rejected": rejected,
    }
    return manifest


# ---------------------------------------------------------------- alignment --

BLOSUM62 = {
    "A": " 4 -1 -2 -2  0 -2 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0",
    "R": "-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3",
    "N": "-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3",
    "D": "-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3",
    "C": " 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1",
    "Q": "-2  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2",
    "E": "-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2",
    "G": " 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3",
    "H": "-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3",
    "I": "-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3",
    "L": "-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1",
    "K": "-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2",
    "M": "-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1",
    "F": "-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1",
    "P": "-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2",
    "S": " 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2",
    "T": " 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0",
    "W": "-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3",
    "Y": "-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1",
    "V": " 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4",
}
_AAS = "ARNDCQEGHILKMFPSTWYV"
BLOSUM = {a: {b: int(v) for b, v in zip(_AAS, BLOSUM62[a].split())}
          for a in BLOSUM62}
GAP_OPEN, GAP_EXTEND = -11, -1


def needlet(seq_a, seq_b):
    """Global Needleman-Wunsch on BLOSUM62. Returns (score, aligned_a, aligned_b).

    The cost column is only ever read through ``alignment_metrics``; the returned
    aligned strings keep gaps as '-'.
    """
    n, m = len(seq_a), len(seq_b)
    if n == 0 or m == 0:
        raise TemplateError("cannot align an empty sequence")
    neg = float("-inf")
    M = [[neg] * (m + 1) for _ in range(n + 1)]
    X = [[neg] * (m + 1) for _ in range(n + 1)]
    Y = [[neg] * (m + 1) for _ in range(n + 1)]
    M[0][0] = 0
    for i in range(1, n + 1):
        X[i][0] = GAP_OPEN + GAP_EXTEND * (i - 1)
    for j in range(1, m + 1):
        Y[0][j] = GAP_OPEN + GAP_EXTEND * (j - 1)
    for i in range(1, n + 1):
        ai = seq_a[i - 1]
        for j in range(1, m + 1):
            s = BLOSUM[ai][seq_b[j - 1]]
            M[i][j] = max(M[i - 1][j - 1], X[i - 1][j - 1], Y[i - 1][j - 1]) + s
            X[i][j] = max(M[i - 1][j] + GAP_OPEN, X[i - 1][j] + GAP_EXTEND,
                          Y[i - 1][j] + GAP_OPEN)
            Y[i][j] = max(M[i][j - 1] + GAP_OPEN, Y[i][j - 1] + GAP_EXTEND,
                          X[i][j - 1] + GAP_OPEN)
    best = max(M[n][m], X[n][m], Y[n][m])
    state = max(((M[n][m], "M"), (X[n][m], "X"), (Y[n][m], "Y")))[1]
    out_a, out_b = [], []
    i, j = n, m
    while i > 0 or j > 0:
        if state == "M":
            out_a.append(seq_a[i - 1]); out_b.append(seq_b[j - 1])
            prev = max(((M[i - 1][j - 1], "M"), (X[i - 1][j - 1], "X"),
                        (Y[i - 1][j - 1], "Y")))[1]
            i, j, state = i - 1, j - 1, prev
        elif state == "X":
            out_a.append(seq_a[i - 1]); out_b.append("-")
            prev = max(((M[i - 1][j] + GAP_OPEN, "M"), (X[i - 1][j] + GAP_EXTEND, "X"),
                        (Y[i - 1][j] + GAP_OPEN, "Y")))[1]
            i, state = i - 1, prev
        else:
            out_a.append("-"); out_b.append(seq_b[j - 1])
            prev = max(((M[i][j - 1] + GAP_OPEN, "M"), (X[i][j - 1] + GAP_OPEN, "X"),
                        (Y[i][j - 1] + GAP_EXTEND, "Y")))[1]
            j, state = j - 1, prev
    return best, "".join(reversed(out_a)), "".join(reversed(out_b))


def alignment_metrics(aligned_a, aligned_b):
    """identity/coverage/similarity plus the residue mapping for superposition."""
    pairs, identical, similar = [], 0, 0
    for ga, gb in zip(aligned_a, aligned_b):
        if ga != "-" and gb != "-":
            pairs.append((ga, gb))
            if ga == gb:
                identical += 1
            if BLOSUM[ga][gb] > 0:
                similar += 1
    n = len(pairs)
    return {
        "aligned_pairs": n,
        "identical": identical,
        "similar": similar,
        "identity": round(identical / n, 4) if n else 0.0,
        "similarity": round(similar / n, 4) if n else 0.0,
        "query_coverage": round(n / len(aligned_a.replace("-", "")), 4) if aligned_a.replace("-", "") else 0.0,
        "template_coverage": round(n / len(aligned_b.replace("-", "")), 4) if aligned_b.replace("-", "") else 0.0,
    }


def align_sequences(query, template_seq):
    if not query or not template_seq:
        raise TemplateError("empty sequence in alignment")
    _, aa, ab = needlet(query, template_seq)
    out = alignment_metrics(aa, ab)
    out["aligned_query"] = aa
    out["aligned_template"] = ab
    return out


def _kabsch_rmsd(P, Q):
    """RMSD of P onto Q after optimal rigid superposition (no pruning)."""
    import numpy as np
    P = np.asarray(P, dtype=float); Q = np.asarray(Q, dtype=float)
    pc, qc = P - P.mean(0), Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(pc.T @ qc)
    d = 1.0 if np.linalg.det(Vt.T @ U.T) >= 0 else -1.0
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return float(np.sqrt((((R @ pc.T).T - qc) ** 2).sum(1).mean()))


def superimpose_ca(query_ca, template_ca, max_rmsd=None):
    """Superpose template CA onto query CA (same length, matched pairs)."""
    if len(query_ca) != len(template_ca):
        raise TemplateError("superposition requires equal-length CA lists")
    if len(query_ca) < 3:
        raise TemplateError("need >= 3 C-alpha pairs to superimpose")
    rmsd = _kabsch_rmsd(query_ca, template_ca)
    out = {"n_ca": len(query_ca), "ca_rmsd": round(rmsd, 3)}
    if max_rmsd is not None and rmsd > max_rmsd:
        raise TemplateError(f"CA RMSD {rmsd:.2f} exceeds limit {max_rmsd}")
    return out


def rank_templates(entries, weights=None):
    """Rank template entries for one aligned query.

    Each entry needs alignment_metrics keys (identity, similarity, query_coverage,
    template_coverage) and resolution_a/method. Missing resolution sorts last.
    No biological threshold is baked in: this is a deterministic ordering.
    """
    w = {"identity": 0.45, "similarity": 0.15, "query_coverage": 0.25,
         "resolution": 0.15}
    if weights:
        unknown = set(weights) - set(w)
        if unknown:
            raise TemplateError(f"unknown ranking weights: {sorted(unknown)}")
        w.update(weights)
    rows = []
    for e in entries:
        res = e.get("resolution_a")
        res_score = 0.0 if res is None else max(0.0, 1.0 - (res / 5.0))
        score = (w["identity"] * e.get("identity", 0.0)
                 + w["similarity"] * e.get("similarity", 0.0)
                 + w["query_coverage"] * e.get("query_coverage", 0.0)
                 + w["resolution"] * res_score)
        row = dict(e)
        row["rank_score"] = round(score, 4)
        row["resolution_score"] = round(res_score, 4)
        rows.append(row)
    rows.sort(key=lambda r: (-r["rank_score"], -(r.get("identity") or 0),
                             r.get("resolution_a") if r.get("resolution_a") is not None else 1e9,
                             r["pdb_id"], r["chain"]))
    for i, r in enumerate(rows, 1):
        r["rank"] = i
    return rows


def validate_target_input(record):
    """Validate a future C5 target spec. Fail closed; no target data exists yet.

    `record` is a dict with at least: name, sequence, and optionally chain,
    numbering_start, construct_notes. Returns a normalized dict.
    """
    if not isinstance(record, dict):
        raise TemplateError("target record must be a mapping")
    for key in ("name", "sequence"):
        if key not in record:
            raise TemplateError(f"target record missing {key!r}")
    name = str(record["name"]).strip()
    if not name:
        raise TemplateError("target name is empty")
    seq = re.sub(r"\s+", "", str(record["sequence"])).upper()
    if not seq:
        raise TemplateError("target sequence is empty")
    bad = sorted(set(seq) - AA)
    if bad:
        raise TemplateError(f"target sequence contains non-AA characters {bad}")
    if not (10 <= len(seq) <= 250):
        raise TemplateError(f"target length {len(seq)} outside the expected 10-250 window")
    start = record.get("numbering_start", 1)
    if not isinstance(start, int) or start < 0:
        raise TemplateError("numbering_start must be a non-negative integer")
    return {
        "name": name,
        "sequence": seq,
        "length": len(seq),
        "numbering_start": start,
        "numbering_end": start + len(seq) - 1,
        "sequence_sha256": hashlib.sha256(seq.encode()).hexdigest(),
        "chain": record.get("chain"),
        "construct_notes": record.get("construct_notes", ""),
    }


def write_manifest(manifest, path):
    """Write the manifest, validating its shape first (fail closed)."""
    missing = [k for k in MANIFEST_SCHEMA["required"] if k not in manifest]
    if missing:
        raise TemplateError(f"manifest missing keys: {missing}")
    for t in manifest["templates"]:
        for k in ("pdb_id", "chain", "method", "length", "sequence_sha256", "coverage"):
            if k not in t:
                raise TemplateError(f"template entry missing {k!r}")
        if t["sequence_sha256"] != hashlib.sha256(t["sequence"].encode()).hexdigest():
            raise TemplateError(f"template {t['pdb_id']}:{t['chain']} sequence hash mismatch")
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    return path

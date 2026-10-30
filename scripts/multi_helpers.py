# =============================================================================
# Helper functions — nothing here needs editing to run a design.
# Collapse this cell (click the gutter arrow) and work from Section 1 down.
#
#   PDB parsing & validation ... read_pdb_ca, parse_hotspots, validate_templates, summarize
#   3D viewer ................. show_templates, show_template, save_templates_html
#   Diagnostics ............... viewer_selftest
#   Export .................... build_target_json, write_target_json
#   Runtime setup ............. check_runtime, init_pyrosetta
# =============================================================================

# --- PDB parsing & validation ------------------------------------------------

VALID_STRATEGIES = ("target", "detarget")


def read_pdb_ca(pdb_path):
    '''Parse first model of a PDB. Returns (chain_order, residues, atom_lines).

    residues: OrderedDict[(chain, resnum)] -> dict(resname, x, y, z)
    '''
    chain_order, residues, atom_lines = [], OrderedDict(), []
    with open(pdb_path) as fh:
        for line in fh:
            rec = line[:6]
            if rec == "ENDMDL":
                break
            if rec != "ATOM  ":
                continue
            atom_lines.append(line.rstrip("\n"))
            if line[12:16].strip() != "CA":
                continue
            chain = line[21]
            try:
                resnum = int(line[22:26])
            except ValueError:
                continue
            if chain not in chain_order:
                chain_order.append(chain)
            residues[(chain, resnum)] = {
                "resname": line[17:20].strip(),
                "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
            }
    return chain_order, residues, atom_lines


def parse_hotspots(spec, default_chain):
    '''Expand a ColabDesign hotspot string into [(chain, resnum), ...].'''
    out = []
    spec = (spec or "").strip()
    if not spec:
        return out
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        lo, hi = token.split("-") if "-" in token else (token, None)
        if lo[0].isalpha():
            chain, start = lo[0], int(lo[1:])
        else:
            chain, start = default_chain, int(lo)
        if hi is None:
            end = start
        else:
            end = int(hi[1:] if hi[0].isalpha() else hi)
        out.extend((chain, r) for r in range(start, end + 1))
    return out


def validate_templates(templates, general, verbose=True):
    '''Normalize + validate. Returns (validated_templates, errors, warnings).'''
    target_chains = [c.strip() for c in general["chains"].split(",") if c.strip()]
    default_chain = target_chains[0]
    errors, warnings, validated = [], [], []
    seen_names = set()

    for i, raw in enumerate(templates, start=1):
        tpl = dict(raw)
        tag = tpl.get("template_name") or f"<template #{i}>"

        if not tpl.get("template_name"):
            errors.append(f"[{tag}] template_name is empty")
        elif tpl["template_name"] in seen_names:
            errors.append(f"[{tag}] duplicate template_name — outputs would overwrite each other")
        else:
            seen_names.add(tpl["template_name"])

        # strategy: '' -> 'target' (same default the pipeline applies)
        strategy = (tpl.get("template_targeting_strategy") or "").strip() or "target"
        if strategy not in VALID_STRATEGIES:
            errors.append(f"[{tag}] template_targeting_strategy={strategy!r} not in {VALID_STRATEGIES}")
        tpl["template_targeting_strategy"] = strategy

        # loss_weighting must survive int()
        try:
            int(str(tpl.get("loss_weighting", "1")))
        except (TypeError, ValueError):
            errors.append(
                f"[{tag}] loss_weighting={tpl.get('loss_weighting')!r} is not integer-valued; "
                "step_multi casts it with int() and would raise"
            )

        # filters file
        fpath = tpl.get("template_mpnn_filters_path", "")
        if not fpath or not Path(fpath).is_file():
            errors.append(f"[{tag}] template_mpnn_filters_path not found: {fpath}")

        # structure
        pdb_path = Path(tpl.get("template_pdb", ""))
        if not pdb_path.is_file():
            errors.append(f"[{tag}] template_pdb not found: {pdb_path}")
            validated.append(tpl)
            continue

        chain_order, residues, _ = read_pdb_ca(pdb_path)
        tpl["_chains"] = chain_order
        tpl["_residues"] = residues

        missing_chains = [c for c in target_chains if c not in chain_order]
        if missing_chains:
            errors.append(
                f"[{tag}] chain(s) {missing_chains} from GENERAL['chains'] absent; PDB has {chain_order}"
            )

        target_res = [k for k in residues if k[0] in target_chains]
        tpl["_target_length"] = len(target_res)
        if target_res:
            nums = [r for _, r in target_res]
            tpl["_resnum_range"] = (min(nums), max(nums))
        else:
            tpl["_resnum_range"] = (None, None)

        extra = [c for c in chain_order if c not in target_chains]
        if extra:
            warnings.append(
                f"[{tag}] chain(s) {extra} present in the PDB but not in GENERAL['chains'] — "
                "they are ignored during design"
            )

        # hotspots
        try:
            hs = parse_hotspots(tpl.get("template_hostspot_residues", ""), default_chain)
        except (ValueError, IndexError) as exc:
            errors.append(f"[{tag}] could not parse template_hostspot_residues: {exc}")
            hs = []
        tpl["_hotspots"] = hs

        unresolved = [f"{c}{r}" for c, r in hs if (c, r) not in residues]
        if unresolved:
            errors.append(
                f"[{tag}] hotspot residue(s) absent from the PDB: {', '.join(unresolved)} "
                "(prep_inputs would raise AssertionError)"
            )
        if not hs:
            warnings.append(f"[{tag}] no hotspots set — ColabDesign will choose the interface freely")

        validated.append(tpl)

    if verbose:
        for w in warnings:
            print(f"WARNING  {w}")
        for e in errors:
            print(f"ERROR    {e}")
        if not errors:
            print(f"OK       {len(validated)} template(s) validated, no blocking problems")
    return validated, errors, warnings


def summarize(templates):
    rows = []
    for t in templates:
        lo, hi = t.get("_resnum_range", (None, None))
        rows.append({
            "name": t.get("template_name"),
            "strategy": t.get("template_targeting_strategy"),
            "weight": t.get("loss_weighting"),
            "chains": ",".join(t.get("_chains", [])),
            "n_res": t.get("_target_length"),
            "resnum_range": f"{lo}-{hi}" if lo is not None else "-",
            "n_hotspots": len(t.get("_hotspots", [])),
            "pdb": Path(t.get("template_pdb", "")).name,
        })
    if pd is not None:
        return pd.DataFrame(rows)
    for r in rows:
        print(r)
    return rows


# --- 3D viewer ---------------------------------------------------------------

# Distinct cartoon colour per template; hotspots are always red.
PALETTE = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860"]
HOTSPOT_COLOR = "#E24A33"
OFF_TARGET_COLOR = "#D9D9D9"


def _hotspots_by_chain(tpl):
    by_chain = {}
    for ch, rn in (tpl.get("_hotspots") or []):
        by_chain.setdefault(ch, []).append(rn)
    return by_chain


def _draw_panel(view, tpl, target_chains, color, cell=None,
                show_surface=False, stick_radius=0.35, show_resi_labels=False):
    """Draw one template into `view`.

    `cell` is the (row, col) of a viewergrid, or None for a single-panel view.
    py3Dmol applies a call to *every* panel unless `viewer=` names one, so the kwarg
    must be omitted entirely rather than passed as None — and it must be a keyword:
    passed positionally it leaks into the 3Dmol call as a stray argument.
    """
    kw = {} if cell is None else {"viewer": cell}

    view.addModel(Path(tpl["template_pdb"]).read_text(), "pdb", **kw)

    # base: everything pale, target chain(s) in the template colour
    view.setStyle({}, {"cartoon": {"color": OFF_TARGET_COLOR, "opacity": 0.55}}, **kw)
    view.setStyle({"chain": target_chains},
                  {"cartoon": {"color": color, "opacity": 0.85}}, **kw)

    # hotspots: side-chain sticks that visibly protrude from the cartoon
    for ch, nums in _hotspots_by_chain(tpl).items():
        sel = {"chain": ch, "resi": nums}
        view.addStyle(sel, {"cartoon": {"color": HOTSPOT_COLOR}}, **kw)
        view.addStyle(sel, {"stick": {"radius": stick_radius,
                                      "colorscheme": "redCarbon"}}, **kw)
        view.addStyle({"chain": ch, "resi": nums, "atom": "CA"},
                      {"sphere": {"radius": stick_radius * 1.7,
                                  "color": HOTSPOT_COLOR}}, **kw)
        if show_surface:
            view.addSurface("VDW", {"opacity": 0.30, "color": HOTSPOT_COLOR}, sel, **kw)
        if show_resi_labels:
            view.addResLabels(sel, {"fontSize": 10, "fontColor": "black",
                                    "backgroundColor": "white",
                                    "backgroundOpacity": 0.7,
                                    "showBackground": True}, **kw)

    view.zoomTo({"chain": target_chains}, **kw)


def _header_html(tpl, color, width):
    """Script-free header: PDB filename on top, template details underneath."""
    pdb_name = Path(tpl.get("template_pdb", "")).name or "(no pdb)"
    n_hs = len(tpl.get("_hotspots") or [])
    lo, hi = tpl.get("_resnum_range", (None, None))
    bits = [
        f"<b>{tpl.get('template_name', '?')}</b>",
        f"chain {','.join(tpl.get('_chains', [])) or '?'}",
        f"{tpl.get('_target_length', '?')} res",
        (f"numbering {lo}–{hi}" if lo is not None else "-"),
        (f"<span style='color:{HOTSPOT_COLOR}'><b>{n_hs} hotspot"
         f"{'s' if n_hs != 1 else ''}</b></span>" if n_hs else
         "<span style='color:#999'>no hotspots</span>"),
        tpl.get("template_targeting_strategy", "target"),
    ]
    return (
        f"<div style='flex:0 0 {width}px;box-sizing:border-box;padding:8px 10px;"
        "border:1px solid #ddd;border-bottom:none;border-radius:6px 6px 0 0;"
        "background:#fafafa;font-family:system-ui,sans-serif;'>"
        f"<div style='font-size:15px;font-weight:600;color:#222;"
        f"border-left:5px solid {color};padding-left:8px;line-height:1.3;'>"
        f"{pdb_name}</div>"
        "<div style='font-size:11px;color:#555;margin-top:3px;padding-left:13px;'>"
        + " &nbsp;·&nbsp; ".join(bits) +
        "</div></div>"
    )


def _header_row(templates, panel_w, ncols):
    """One flex row of headers, each exactly one panel wide so it sits over its panel.

    ponytail: only the first grid row lines up when templates wrap onto several rows —
    a single canvas cannot have HTML interleaved between its rows. Irrelevant at the
    2-3 templates this is built for; for more than that, drop the headers and
    read the in-canvas labels instead.
    """
    chips = "".join(_header_html(t, PALETTE[i % len(PALETTE)], panel_w)
                    for i, t in enumerate(templates))
    return (f"<div style='display:flex;flex-wrap:wrap;margin-top:10px;"
            f"max-width:{panel_w * ncols}px;'>{chips}</div>")


def show_templates(templates, general=None, panel=(420, 380), ncols=None, linked=False,
                   show_surface=False, stick_radius=0.35, show_resi_labels=False):
    """All templates side by side, as one viewergrid in a single output.

    One viewer object means one <script>, which is the render path VSCode handles most
    reliably. The headers go out as a separate HTML output above it — wrapping py3Dmol's
    own HTML inside a document of ours stops its bootstrap script from running.

    panel  -- (width, height) of each panel in px
    ncols  -- columns; defaults to one row holding every template
    linked -- True ties the cameras together, so rotating one panel rotates all of them.
              Off by default: linked panels share a zoom, which crops templates that are
              not already in a common frame. Turn it on for orthologs of equal length.
    """
    general = general or GENERAL
    templates = [t for t in templates if Path(t.get("template_pdb", "")).is_file()]
    if not templates:
        print("No readable template PDBs to display.")
        return

    if not any(t.get("_hotspots") for t in templates):
        print("Note: no hotspots set on any template — the red stick styling "
              "appears once template_hostspot_residues is filled in (Section 2).")

    target_chains = [c.strip() for c in general["chains"].split(",") if c.strip()]
    n = len(templates)
    ncols = max(1, min(ncols or n, n))
    nrows = math.ceil(n / ncols)

    view = py3Dmol.view(viewergrid=(nrows, ncols), linked=linked,
                        width=panel[0] * ncols, height=panel[1] * nrows)
    for i, tpl in enumerate(templates):
        _draw_panel(view, tpl, target_chains, PALETTE[i % len(PALETTE)],
                    cell=(i // ncols, i % ncols), show_surface=show_surface,
                    stick_radius=stick_radius, show_resi_labels=show_resi_labels)
        # in-canvas name, so an exported PNG or a wrapped header row stays readable
        view.addLabel(tpl.get("template_name", "?"),
                      {"useScreen": True, "screenOffset": {"x": 6, "y": 6},
                       "fontSize": 12, "fontColor": "black", "backgroundColor": "white",
                       "backgroundOpacity": 0.85, "inFront": True},
                      {}, viewer=(i // ncols, i % ncols))
    view.setBackgroundColor("white")

    display(HTML(_header_row(templates, panel[0], ncols)))
    view.show()
    return view


def show_template(tpl, general=None, size=(760, 560), show_surface=False,
                  spin=False, stick_radius=0.38, show_resi_labels=True):
    """Single template, larger, residue labels on by default."""
    general = general or GENERAL
    target_chains = [c.strip() for c in general["chains"].split(",") if c.strip()]
    view = py3Dmol.view(width=size[0], height=size[1])
    _draw_panel(view, tpl, target_chains, PALETTE[0], cell=None,
                show_surface=show_surface, stick_radius=stick_radius,
                show_resi_labels=show_resi_labels)
    view.setBackgroundColor("white")
    if spin:
        view.spin(True)
    display(HTML(_header_row([tpl], size[0], 1)))
    view.show()
    return view


def save_templates_html(templates, path, general=None, panel=(420, 380), ncols=None,
                        linked=False, show_surface=False, stick_radius=0.35,
                        show_resi_labels=False):
    """Write the same grid to a standalone HTML file (3Dmol.js from CDN).

    Scripts run normally in a real browser page, so this always works even when the
    notebook renderer does not.
    """
    general = general or GENERAL
    templates = [t for t in templates if Path(t.get("template_pdb", "")).is_file()]
    if not templates:
        print("No readable template PDBs to write.")
        return

    target_chains = [c.strip() for c in general["chains"].split(",") if c.strip()]
    n = len(templates)
    ncols = max(1, min(ncols or n, n))
    nrows = math.ceil(n / ncols)

    view = py3Dmol.view(viewergrid=(nrows, ncols), linked=linked,
                        width=panel[0] * ncols, height=panel[1] * nrows)
    for i, tpl in enumerate(templates):
        _draw_panel(view, tpl, target_chains, PALETTE[i % len(PALETTE)],
                    cell=(i // ncols, i % ncols), show_surface=show_surface,
                    stick_radius=stick_radius, show_resi_labels=show_resi_labels)
        view.addLabel(tpl.get("template_name", "?"),
                      {"useScreen": True, "screenOffset": {"x": 6, "y": 6},
                       "fontSize": 12, "fontColor": "black", "backgroundColor": "white",
                       "backgroundOpacity": 0.85, "inFront": True},
                      {}, viewer=(i // ncols, i % ncols))
    view.setBackgroundColor("white")

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("<!doctype html><html><head><meta charset='utf-8'>"
                    "<title>BindCraft templates</title></head>"
                    "<body style='margin:16px;background:#fff;'>"
                    + _header_row(templates, panel[0], ncols)
                    + view._make_html() + "</body></html>")
    print(f"wrote {path}")
    return path


# --- Diagnostics -------------------------------------------------------------

def viewer_selftest():
    """Minimal py3Dmol render - isolates 'my code is wrong' from 'the webview can't
    run 3Dmol at all'.

    Emits the simplest possible viewer: 3 atoms, one style call, nothing else. If this
    panel is BLANK, the problem is environmental (notebook not trusted, webview CSP
    blocking the CDN, kernel/extension state) and no change to the viewer code will fix
    it. If a small molecule appears, 3Dmol works and the fault is in the styling code.
    """
    tiny = (
        "ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  GLY A   1       1.458   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   GLY A   1       2.009   1.420   0.000  1.00  0.00           C\n"
        "END\n"
    )
    v = py3Dmol.view(width=320, height=240)
    v.addModel(tiny, "pdb")
    v.setStyle({}, {"sphere": {"radius": 0.5}})
    v.zoomTo()
    print("Three spheres below means 3Dmol.js is working in this notebook.")
    print("If the area below is blank/empty, it is an environment problem - see docstring.")
    v.show()


# --- Export ------------------------------------------------------------------

def build_target_json(general, templates):
    payload = {"general_information": dict(general)}
    for i, tpl in enumerate(templates, start=1):
        payload[f"template_{i}_information"] = {
            "template_name":               tpl["template_name"],
            "template_pdb":                tpl["template_pdb"],
            "template_hostspot_residues":  tpl.get("template_hostspot_residues", ""),
            "template_targeting_strategy": tpl.get("template_targeting_strategy", "target"),
            "template_mpnn_filters_path":  tpl["template_mpnn_filters_path"],
            "loss_weighting":              str(tpl.get("loss_weighting", "1")),
        }
    return payload


def write_target_json(general, templates, filename=None, overwrite=False):
    payload = build_target_json(general, templates)
    filename = filename or f"{general['binder_name']}.json"
    out = SETTINGS_TARGET / filename
    if out.exists() and not overwrite:
        raise FileExistsError(f"{out} exists — pass overwrite=True to replace it")
    Path(general["design_path"]).mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=4))
    print(f"wrote {out}")
    return out


# --- Runtime setup -----------------------------------------------------------
# Everything below needs a GPU node. JAX / ColabDesign / PyRosetta are imported
# lazily *inside* these functions so Sections 0-4 stay runnable on a login node.

AF_MODEL_NAMES = {True:  [f"model_{k}_multimer_v3" for k in range(1, 6)],
                  False: [f"model_{k}_ptm" for k in range(1, 6)]}


GPU_HELP = """No GPU available to JAX — binder design needs one.

  * On a cluster: request a GPU node (`srun --gres=gpu:1 --pty bash`) and restart the
    kernel there. `nvidia-smi` showing a card is not enough — the kernel must be on the
    node that holds the allocation.
  * Seen on this cluster: cuInit intermittently fails with CUDA_ERROR_NOT_INITIALIZED
    for several minutes at a time while `nvidia-smi` stays healthy, then recovers on its
    own. If the setup is unchanged, wait and re-run this cell before changing
    anything. Switching conda environments does not help — the fault is below Python."""


def check_gpu(verbose=True):
    """The pipeline's check_jax_gpu(), raising something actionable instead of exit().

    jax.devices() itself raises when no backend initialises at all, which is the most
    common failure — so that case has to be caught here or the advice below never prints.
    """
    import jax
    try:
        devices = jax.devices()
    except RuntimeError as exc:
        raise RuntimeError(f"{GPU_HELP}\n\n  JAX reported: {exc}") from exc
    if verbose:
        for d in devices:
            print(f"    {d.device_kind} ({d.platform})")
    if not any(d.platform == "gpu" for d in devices):
        raise RuntimeError(f"{GPU_HELP}\n\n  JAX sees only: "
                           + ", ".join(f"{d.device_kind} ({d.platform})" for d in devices))
    return devices


def check_af_params(af_params_dir, use_multimer=True):
    """ColabDesign loads <af_params_dir>/params/params_<model>.npz — check they are there.

    This is the one dependency too large to keep next to the notebook (~5 GB). Point
    AF_PARAMS_DIR at any existing AlphaFold2 / BindCraft / ColabDesign weights directory.
    """
    root = Path(af_params_dir)
    names = AF_MODEL_NAMES[bool(use_multimer)]
    missing = [n for n in names if not (root / "params" / f"params_{n}.npz").is_file()]
    if missing:
        raise FileNotFoundError(
            f"AlphaFold2 weights missing under {root / 'params'}: "
            + ", ".join(f"params_{n}.npz" for n in missing)
            + "\nSet AF_PARAMS_DIR to a directory that contains a 'params/' folder with these, "
              "or fetch them:\n"
              "  mkdir -p params && cd params && "
              "wget -qO- https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar | tar x"
        )
    return [str(root / "params" / f"params_{n}.npz") for n in names]


def check_executables(**paths):
    """Each path must exist and be executable — same check the pipeline makes, named."""
    for label, path in paths.items():
        if not Path(path).is_file():
            raise FileNotFoundError(f"{label} not found: {path}")
        if not os.access(path, os.X_OK):
            raise PermissionError(f"{label} is not executable: {path}\n  fix: chmod +x {path}")
    return True


def check_runtime(af_params_dir, advanced_json, dssp_path, dalphaball_path,
                  use_multimer=None, verbose=True):
    """Run every startup check `multi_binder_design.main()` makes, before committing a GPU.

    Returns (advanced_settings, design_models, prediction_models, multimer_validation).
    Raises on anything that would otherwise fail deep inside a trajectory.
    """
    say = print if verbose else (lambda *a, **k: None)

    say("GPU")
    devices = check_gpu(verbose=verbose)

    say("advanced settings")
    with open(advanced_json) as fh:
        advanced_settings = json.load(fh)
    say(f"    {Path(advanced_json).resolve()}  ({len(advanced_settings)} keys)")

    if use_multimer is None:
        use_multimer = advanced_settings["use_multimer_design"]

    say("AlphaFold2 weights")
    found = check_af_params(af_params_dir, use_multimer)
    say(f"    {len(found)} x {'multimer_v3' if use_multimer else 'ptm'} under {Path(af_params_dir) / 'params'}")

    say("executables")
    check_executables(dssp=dssp_path, DAlphaBall=dalphaball_path)
    say(f"    dssp        {dssp_path}")
    say(f"    DAlphaBall  {dalphaball_path}")

    # the pipeline overwrites these three from its own folder layout; set them explicitly
    advanced_settings["af_params_dir"] = str(af_params_dir)
    advanced_settings["dssp_path"] = str(dssp_path)
    advanced_settings["dalphaball_path"] = str(dalphaball_path)

    # normalize omit_AAs exactly as perform_advanced_settings_check() does
    if advanced_settings.get("omit_AAs") in (None, False, ""):
        advanced_settings["omit_AAs"] = None
    elif isinstance(advanced_settings["omit_AAs"], str):
        advanced_settings["omit_AAs"] = advanced_settings["omit_AAs"].strip()

    # load_af2_models(): design and prediction must not share models
    if use_multimer:
        design_models, prediction_models, multimer_validation = [0, 1, 2, 3, 4], [0, 1], False
    else:
        design_models, prediction_models, multimer_validation = [0, 1], [0, 1, 2, 3, 4], True

    say("model split")
    say(f"    design {design_models}  prediction {prediction_models}  "
        f"multimer_validation={multimer_validation}")
    say(f"\nready — {len(devices)} device(s), omit_AAs={advanced_settings['omit_AAs']!r}")
    return advanced_settings, design_models, prediction_models, multimer_validation


def init_pyrosetta(advanced_settings, extra_flags=""):
    """Initialise PyRosetta with the pipeline's flags. Safe to re-run; it re-inits."""
    import pyrosetta as pr
    flags = (
        "-ignore_unrecognized_res -ignore_zero_occupancy -mute all "
        f"-holes:dalphaball {advanced_settings['dalphaball_path']} "
        "-corrections::beta_nov16 true -relax:default_repeats 1 "
        + extra_flags
    ).strip()
    pr.init(flags)
    print(f"PyRosetta initialised\n  {flags}")
    return pr

"""What a learn records as its input (the marker finder's input and --in): learn from a file, the file's path; learn from
an in-memory AnnData, 'in-memory AnnData' with its shape; each learn of transfer (cells in memory), the input file's
path and sha256 and how the cells were used: thinned in draw i with the depth-matching seed s, or not thinned."""
from __future__ import annotations

import hashlib
import json
import zipfile

import numpy as np
from scipy import sparse

from conftest import make_synthetic_adata


def _sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def _recorded(path):
    with zipfile.ZipFile(path) as z:
        mp = json.loads(z.read("marker_finder_params.json"))
        prov = json.loads(z.read("provenance.json"))
    argv = mp["context"]["marker_finder_argv"]
    assert argv[argv.index("--in") + 1] == mp["marker_finder"]["input"]
    return mp["marker_finder"]["input"], prov


def test_learn_records_its_input(tmp_path):
    from scarlet import api
    a = make_synthetic_adata(seed=5)
    p = tmp_path / "ref.h5ad"
    a.write_h5ad(p)
    api.learn(str(p), "cell_type", str(tmp_path / "file.scarlet"), blocklist=None)
    api.learn(a, "cell_type", str(tmp_path / "mem.scarlet"), blocklist=None)
    assert _recorded(tmp_path / "file.scarlet")[0] == str(p)
    assert _recorded(tmp_path / "mem.scarlet")[0] == f"in-memory AnnData ({a.n_obs} cells x {a.n_vars} genes)"


def test_transfer_learns_record_their_source(tmp_path, monkeypatch):
    import anndata as ad
    from scarlet import api
    a = make_synthetic_adata(seed=3)
    idx = np.random.default_rng(0).permutation(a.n_obs)
    ref, q = a[np.sort(idx[: a.n_obs // 2])].copy(), a[np.sort(idx[a.n_obs // 2:])].copy()
    X = q.X.toarray()
    o = np.where(np.arange(q.n_obs) % 2 == 0, "shallow", "deep").astype(object)
    X[o == "shallow"] = X[o == "shallow"] * (np.random.default_rng(1).random(X[o == "shallow"].shape) < 0.5)
    q.X = sparse.csr_matrix(X.astype(np.float32))
    q.obs["origin"] = o
    rp, qp = tmp_path / "ref.h5ad", tmp_path / "query.h5ad"
    ref.write_h5ad(rp)
    q.write_h5ad(qp)
    src = {"reference": f"{rp} (sha256 {_sha(rp)})", "query": f"{qp} (sha256 {_sha(qp)})"}
    seen = []
    orig = api._learn_into

    def spy(model_path, data_path, labels, kw, resolved, command, params, inputs, **k):
        seen.append((params.get("step"), data_path))
        return orig(model_path, data_path, labels, kw, resolved, command, params, inputs, **k)
    monkeypatch.setattr(api, "_learn_into", spy)
    out = tmp_path / "t"
    api.transfer(str(rp), "cell_type", str(qp), str(out), origin_col="origin", blocklist=None, draws=2)
    rec = json.loads((out / "transfer_record.json").read_text())
    plans = {d["draw"]: {p["origin"]: p for p in d["origins"]} for d in rec["depth_match"]}
    assert any(p["thinned"] == "reference" for p in plans[1].values())

    def expected(step, origin, draw):
        if step.startswith("refine"):
            return f"{src['query']}: in memory, not thinned, draw {draw}, refinement (origin {origin})"
        if step.startswith("qc"):
            return f"{src['query']}: in memory, not thinned, transfer check (origin {origin})"
        p = plans.get(draw, {}).get(origin)
        if p is not None and p["thinned"] == "reference":
            return f"{src['reference']}: in memory, thinned, draw {draw}, seed {p['seed']} (origin {origin})"
        return f"{src['reference']}: in memory, not thinned, draw {draw}"

    saved = sorted(out.glob("*.scarlet"))
    assert {f.name.split("_")[0].split(".")[0] for f in saved} == {"reference", "refine"}
    kinds = set()
    for f in saved:   # the written models (draw 1)
        inp, prov = _recorded(f)
        step, origin = prov["parameters"]["step"], prov["inputs"].get("origin")
        assert inp == expected(step, origin, 1), (f.name, inp)
        kinds.add(": in memory, thinned, " in inp)
    assert kinds == {True, False}   # both a thinned and an unthinned learn written
    steps = {s.split(" ")[0] for s, _ in seen}
    assert {"reference", "refine", "qc"} <= steps, steps
    for step, name in seen:   # every learn of the run, draw 2 and the transfer check included
        assert name.startswith(src["reference"] + ": in memory, ") or name.startswith(src["query"] + ": in memory, ")
        assert "<AnnData>" not in name and ".h5ad (sha256 " in name
    assert sum(name.endswith("transfer check (origin deep)") or name.endswith("transfer check (origin shallow)")
               for _, name in seen) >= 1
    assert ad.read_h5ad(rp).n_obs == ref.n_obs   # the inputs are not changed

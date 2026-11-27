"""transfer on synthetic data: depth matching and the conserved-marker QC (the refinement's anchors: see
test_transfer_draws.py)."""
from __future__ import annotations

import numpy as np
from scipy import sparse
import pytest

from conftest import make_synthetic_adata


def test_transfer_depth_match_and_qc(tmp_path):
    """transfer with depth matching thins the deeper side per origin and records it; the QC flags identities
    whose reference markers are not conserved in the query programme (pre-registered cut-off 0.2)."""
    from scarlet import api, transfer_steps
    ref = make_synthetic_adata(seed=2)
    q = make_synthetic_adata(seed=3)
    rng = np.random.default_rng(0)                                 # 0.6.0: depth = genes detected per cell
    Xq = q.X.toarray()
    q.X = sparse.csr_matrix((Xq * (rng.random(Xq.shape) < 0.5)).astype(np.float32))   # the query detects ~half
    q.obs["batch"] = np.where(np.arange(q.n_obs) % 2, "b1", "b2")
    res = api.transfer(ref, "cell_type", q, str(tmp_path / "t"), origin_col="batch", depth_match="capture",
                       n_iterations=10, draws=7)
    assert [d["draw"] for d in res.record["depth_match"]] == list(range(1, 8))   # one plan per draw
    plans = {p["origin"]: p for p in res.record["depth_match"][0]["origins"]}   # 0.6.0: capture (d2)
    assert set(plans) == {"b1", "b2"} and all(p["thinned"] == "reference" and p["method"] == "capture"
                                               for p in plans.values())
    assert set(res.reference_models) == {"b1", "b2"}
    assert {"origin", "identity", "share_conserved", "doubtful_transfer"} <= set(res.qc.columns)
    ok = res.qc.dropna(subset=["share_conserved"])
    assert (ok.doubtful_transfer == (ok.share_conserved < transfer_steps.QC_FLAG)).all()
    assert res.record["qc_performance"]["precision_for_true_f1_below_0.5"] == 0.61
    from scarlet.input_checks import NonCountInputError
    with pytest.raises(NonCountInputError):                        # depth matching needs non-negative values
        q2 = q.copy()
        q2.X = q2.X * -0.5
        api.transfer(ref, "cell_type", q2, str(tmp_path / "t2"), depth_match="capture", n_iterations=5, draws=7)

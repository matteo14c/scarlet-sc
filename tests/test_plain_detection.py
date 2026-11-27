"""Both technology presets read detection (count > 0): there is no binarisation. assign with a dense model reads the
query as given, so a query of counts and the same query as 0 / 1 detection get the same calls; the annotated output
keeps the input's values; the settings record names the detection rule and no binarisation."""
import anndata as ad
import numpy as np
import pandas as pd
import pytest

from scarlet import api
from synthetic_dense import _data

KW = dict(n_iterations=6, n_cells=20)


@pytest.fixture(scope="module")
def dense(tmp_path_factory):
    d = tmp_path_factory.mktemp("r063")
    A = _data(d / "counts.h5ad")
    assert (np.diff(A.X.indptr) > 5000).sum() > 20   # deep cells, as on plate data
    api.learn(str(d / "counts.h5ad"), "label", str(d / "hd.scarlet"), technology="dense", seed=42, **KW)
    return d, A


def _q(A):
    Q = A.copy()
    Q.obs = pd.DataFrame(index=Q.obs_names)
    return Q


def _detection(A):
    B = _q(A)
    B.X = (B.X > 0).astype(np.float32)
    return B


def test_dense_query_read_as_detection(dense):
    d, A = dense
    pd.testing.assert_frame_equal(api.assign(str(d / "hd.scarlet"), _q(A), None),
                                  api.assign(str(d / "hd.scarlet"), _detection(A), None))


def test_dense_learn_reads_detection(dense, tmp_path):
    d, A = dense
    B = _detection(A)
    B.obs["label"] = A.obs["label"].to_numpy()
    B.write_h5ad(tmp_path / "det.h5ad")
    m1 = api.load_model(str(d / "hd.scarlet"))
    m2 = api.learn(str(tmp_path / "det.h5ad"), "label", None, technology="dense", seed=42, **KW)
    assert m1.markers.equals(m2.markers)


def test_settings_record_detection(dense):
    d, _ = dense
    t = api.load_model(str(d / "hd.scarlet")).settings["resolved"]["technology"]
    assert t == {"name": "dense", "detection": "count > 0", "marker_gamma": 0.0}


@pytest.mark.parametrize("from_file", [False, True])
def test_the_annotated_output_keeps_the_input_values(dense, tmp_path, from_file):
    d, A = dense
    Q = _q(A)
    src = tmp_path / "q.h5ad"
    Q.write_h5ad(src)
    lab = api.assign(str(d / "hd.scarlet"), str(src) if from_file else Q, str(tmp_path / "out.h5ad"))
    O = ad.read_h5ad(tmp_path / "out.h5ad")
    assert (O.X != A.X).nnz == 0
    assert "query" not in O.uns["scarlet_assign"]
    assert (O.obs["scarlet_best"].astype(str).to_numpy() == lab["scarlet_best"].astype(str).to_numpy()).all()
    if not from_file:
        assert (Q.X != A.X).nnz == 0   # the caller's object is not modified

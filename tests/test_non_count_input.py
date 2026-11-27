"""Non-count input stops SCARLET: a matrix with negative values (e.g. batch-corrected) is not read as detection;
the marker finder, the classifier, every command (CLI and API) stop with a clear message."""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from conftest import make_synthetic_adata
from scarlet.input_checks import NON_COUNT_MESSAGE, NonCountInputError, check_non_negative

REPO = Path(__file__).resolve().parents[1]
CLI = ["-c", "import sys; sys.argv[0] = 'scarlet'; from scarlet.cli import main; main()"]   # the scarlet command
MSG = "non-count input, e.g. batch-corrected, is outside SCARLET's scope"


def _corrected(seed=0):
    """The synthetic counts with a few entries made negative, as after a batch correction."""
    a = make_synthetic_adata(seed)
    X = a.X.tocsr().copy()
    X.data[::7] = -0.25
    a.X = X
    return a


def test_message_is_the_decided_one():
    assert NON_COUNT_MESSAGE == MSG


@pytest.mark.parametrize("fmt", ["csr", "csc", "dense"])
def test_check_raises_on_negative_values(fmt):
    X = _corrected().X
    X = X.toarray() if fmt == "dense" else getattr(X, f"to{fmt}")()
    with pytest.raises(NonCountInputError, match="outside SCARLET's scope"):
        check_non_negative(X)


def test_check_accepts_counts_and_normalised_values():
    a = make_synthetic_adata()
    check_non_negative(a.X)
    check_non_negative(np.log1p(a.X.toarray() / 3.0))
    check_non_negative(sparse.csr_matrix((5, 4)))


def test_learn_api_stops(tmp_path):
    from scarlet import api
    with pytest.raises(NonCountInputError):
        api.learn(_corrected(), "cell_type", str(tmp_path / "m.scarlet"), n_iterations=5, n_cells=20)
    assert not (tmp_path / "m.scarlet").exists()


def test_classifier_stops():
    from scarlet.classifier import classify
    a = _corrected()
    genes = list(a.var_names)
    universe = pd.DataFrame({"A": 0.5, "B": 0.4, "bg": 0.1}, index=genes)
    with pytest.raises(NonCountInputError):
        classify(a, universe)


def test_assign_transfer_api_stop(tmp_path):
    from scarlet import api
    good = make_synthetic_adata()
    api.learn(good, "cell_type", str(tmp_path / "m.scarlet"), n_iterations=5, n_cells=20)
    with pytest.raises(NonCountInputError):
        api.assign(str(tmp_path / "m.scarlet"), _corrected(), str(tmp_path / "a.h5ad"))
    with pytest.raises(NonCountInputError):
        api.transfer(good, "cell_type", _corrected(), str(tmp_path / "tr"), n_iterations=5, n_cells=20)


def _cli(args):
    env = dict(os.environ, PYTHONPATH=str(REPO))
    return subprocess.run([sys.executable, *map(str, args)], capture_output=True, text=True, env=env, timeout=600)


def test_cli_learn_assign_stop(tmp_path):
    h5 = tmp_path / "corrected.h5ad"
    _corrected().write_h5ad(h5)
    r = _cli([*CLI, "learn", "--reference", h5, "--labels", "cell_type", "--n-iterations", "5",
              "--model", tmp_path / "m.scarlet"])
    assert r.returncode != 0 and MSG in r.stderr and "Traceback" not in r.stderr
    assert not (tmp_path / "m.scarlet").exists()
    # the dense preset reads detection as the sparse one does: stops at the marker search
    r = _cli([*CLI, "learn", "--reference", h5, "--labels", "cell_type", "--n-iterations", "5",
              "--technology", "dense", "--model", tmp_path / "m2.scarlet"])
    assert r.returncode != 0 and MSG in r.stderr and "Traceback" not in r.stderr
    good = tmp_path / "good.h5ad"
    make_synthetic_adata().write_h5ad(good)
    r = _cli([*CLI, "learn", "--reference", good, "--labels", "cell_type", "--n-iterations", "5",
              "--model", tmp_path / "g.scarlet"])
    assert r.returncode == 0, r.stderr[-2000:]
    r = _cli([*CLI, "assign", "--model", tmp_path / "g.scarlet", "--data", h5, "--output",
              tmp_path / "a.h5ad"])
    assert r.returncode != 0 and MSG in r.stderr and "Traceback" not in r.stderr and not (tmp_path / "a.h5ad").exists()

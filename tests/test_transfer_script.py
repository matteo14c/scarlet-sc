"""transfer with ncpu above 1 from a script file. The draws run in new Python processes, which import the calling
script; a script that calls transfer without the `if __name__ == "__main__":` guard would run it again in every
worker. The workers then stop before reading any data, and the caller gets one clear error instead of a broken
process pool. With the guard, ncpu=2 gives the same files as ncpu=1."""
import multiprocessing as mp
import os
import subprocess
import sys
import textwrap

import pytest

from conftest import make_synthetic_adata

CALL = 'scarlet.transfer(r"{ref}", "cell_type", r"{qry}", r"{out}", ncpu={ncpu})'


def _run(tmp_path, name, body):
    script = tmp_path / f"{name}.py"
    script.write_text(textwrap.dedent(body))
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(p for p in sys.path if p))
    return subprocess.run([sys.executable, str(script)], capture_output=True, text=True, env=env, timeout=900,
                          cwd=tmp_path)


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    d = tmp_path_factory.mktemp("script")
    make_synthetic_adata(seed=5).write_h5ad(d / "ref.h5ad")
    make_synthetic_adata(seed=6).write_h5ad(d / "qry.h5ad")
    return d


def test_unguarded_script_gets_one_clear_error(data):
    out = data / "out_unguarded"
    call = CALL.format(ref=data / "ref.h5ad", qry=data / "qry.h5ad", out=out, ncpu=2)
    p = _run(data, "unguarded", f"import scarlet\n{call}\n")
    assert p.returncode != 0
    last = p.stderr.strip().splitlines()[-1]
    assert last.startswith("RuntimeError: transfer with ncpu above 1 runs its draws in new Python processes"), last
    assert 'if __name__ == "__main__":' in last and "ncpu=1" in last
    assert "BrokenProcessPool" not in p.stderr
    assert not (out / "labels.tsv").exists()


def test_unguarded_script_with_learn_and_assign_first_exits(data):
    """A script that calls learn and assign before transfer, without the guard (as the getting-started example would).
    The draw workers must stop at the script's first SCARLET call (no learn or assign again in the workers) and the
    caller must exit non-zero with the message, not wait (a hang would hit the subprocess timeout)."""
    out = data / "out_unguarded_all"
    model = data / "unguarded_all.scarlet"
    body = (f'import scarlet\n'
            f'scarlet.learn(r"{data / "ref.h5ad"}", "cell_type", r"{model}", ncpu=2)\n'
            f'scarlet.assign(r"{model}", r"{data / "qry.h5ad"}", r"{data / "unguarded_all.h5ad"}", ncpu=2)\n'
            + CALL.format(ref=data / "ref.h5ad", qry=data / "qry.h5ad", out=out, ncpu=2) + "\n")
    p = _run(data, "unguarded_all", body)
    assert p.returncode != 0
    last = p.stderr.strip().splitlines()[-1]
    assert last.startswith("RuntimeError: transfer with ncpu above 1 runs its draws in new Python processes"), last
    assert "BrokenProcessPool" not in p.stderr
    assert p.stdout.count("learn finished") == 1, p.stdout[-3000:]   # the caller's learn only; no worker re-ran it
    assert not (out / "labels.tsv").exists()


@pytest.mark.skipif("fork" not in mp.get_all_start_methods(), reason="learn and assign spawn without fork")
def test_unguarded_learn_and_assign_ncpu2_run(data):
    """learn and assign at ncpu 2 use fork pools where the platform has fork: a script without the guard just runs."""
    model = data / "unguarded_la.scarlet"
    body = (f'import scarlet\n'
            f'scarlet.learn(r"{data / "ref.h5ad"}", "cell_type", r"{model}", ncpu=2)\n'
            f'scarlet.assign(r"{model}", r"{data / "qry.h5ad"}", r"{data / "unguarded_la.h5ad"}", ncpu=2)\n')
    p = _run(data, "unguarded_la", body)
    assert p.returncode == 0, p.stderr[-2000:]
    assert (data / "unguarded_la.h5ad").exists()


def test_guarded_script_ncpu2_equals_ncpu1(data):
    outs = {}
    for ncpu in (1, 2):
        out = data / f"out_guarded{ncpu}"
        call = CALL.format(ref=data / "ref.h5ad", qry=data / "qry.h5ad", out=out, ncpu=ncpu)
        p = _run(data, f"guarded{ncpu}", f'import scarlet\nif __name__ == "__main__":\n    {call}\n')
        assert p.returncode == 0, p.stderr[-2000:]
        outs[ncpu] = out
    for f in ("labels.tsv", "first_pass.tsv", "qc.tsv"):
        assert (outs[1] / f).read_bytes() == (outs[2] / f).read_bytes(), f

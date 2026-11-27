"""Release 0.6.3: each marker search frees its data when it ends. With a blocklist that removes genes, the search works
on a column view of its input; the view and its parent form reference cycles that Python's cyclic collector frees only
when it runs, so a transfer's twenty or so searches used to pile them up (mouse brain: 7.2 GB instead of 3.1 GB with
one CPU). The engine now collects once per search; the results do not change."""
import gc
import weakref

from scarlet import api, marker_finder
from synthetic_dense import _data

KW = dict(n_iterations=6, n_cells=20)


def test_marker_search_data_freed_when_learn_returns(tmp_path, monkeypatch):
    A = _data(tmp_path / "human.h5ad")
    bl = tmp_path / "list.tsv"
    bl.write_text("Symbol\tCategory\n" + "".join(f"{g}\tcore_stress\n" for g in list(A.var_names[:3])))
    refs = []
    orig = marker_finder.apply_universal_blocklist

    def keep(adata, args):
        out = orig(adata, args)
        refs.append((weakref.ref(adata), weakref.ref(out), out is not adata))
        return out
    monkeypatch.setattr(marker_finder, "apply_universal_blocklist", keep)
    gc.disable()   # only an explicit collection may free the cycles
    try:
        api.learn(str(tmp_path / "human.h5ad"), "label", str(tmp_path / "m.scarlet"), blocklist=str(bl), **KW)
        assert refs and all(made for _, _, made in refs)   # the blocklist removed genes: a view was made
        assert all(i() is None and v() is None for i, v, _ in refs)
    finally:
        gc.enable()

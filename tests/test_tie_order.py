"""The stable genes of a bootstrap round when scores tie at the cut-off: a fixed order (the more specific gene first,
then the gene name as a string), not argpartition's CPU-dependent choice."""
import numpy as np

from scarlet.marker_finder import _gene_name_rank, _top_positive_tied


def test_ties_take_the_more_specific_gene_then_the_name():
    names = ["GENE_C", "GENE_A", "GENE_B", "GENE_D", "GENE_E"]
    rank = _gene_name_rank(names, 5)
    assert list(rank) == [2, 0, 1, 3, 4]
    # class 0: gene 3 above the cut-off; genes 0, 1, 2 tied at 1.0 for k - 1 = 1 slot
    s0 = np.array([1.0, 1.0, 1.0, 2.0, 0.0])
    s1 = np.array([0.5, 0.2, 0.2, 0.0, 3.0])   # specificity of the tied genes: 0.5, 0.2, 0.2
    sel = _top_positive_tied(np.vstack([s0, s1]), 0, 2, rank)
    assert sorted(sel.tolist()) == [1, 3]       # genes 1 and 2 equally specific: the name decides (GENE_A < GENE_B)
    s1b = np.array([0.1, 0.2, 0.2, 0.0, 3.0])   # now gene 0 is the most specific
    assert sorted(_top_positive_tied(np.vstack([s0, s1b]), 0, 2, rank).tolist()) == [0, 3]


def test_without_a_straddling_tie_nothing_changes():
    s = np.array([3.0, 2.0, 1.0, 0.0])
    assert sorted(_top_positive_tied(np.vstack([s, s]), 0, 2, np.arange(4)).tolist()) == [0, 1]
    assert sorted(_top_positive_tied(np.vstack([s, s]), 0, 5, np.arange(4)).tolist()) == [0, 1, 2]
    assert _top_positive_tied(np.vstack([np.zeros(4), s]), 0, 2, np.arange(4)) is None

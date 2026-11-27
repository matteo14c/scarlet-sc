"""The value array of detection patterns (boolean CSR matrices whose stored values are all True).

A CSR matrix needs one value per stored entry, and a detection pattern's are all True: patterns take a slice of a
shared read-only all-True buffer instead of allocating their own. SciPy copies a value array that is a view of less
than half of its base (scipy._lib._util._prune_array, called by every CSR constructor), so there is one buffer per
size class [2^(k-1), 2^k), as long as the longest pattern of its class seen so far: every slice is more than half of
its buffer and is kept as it is, and no buffer is longer than a pattern that would have had its own array. Read-only,
so any in-place write to a pattern's values fails loudly."""
from __future__ import annotations

import numpy as np

_BUFFERS = {}   # size class (n.bit_length()) -> read-only all-True buffer


def all_true(n):
    """A read-only boolean array of n True values (a slice of the shared buffer of n's size class)."""
    n = int(n)
    k = n.bit_length()
    buf = _BUFFERS.get(k)
    if buf is None or len(buf) < n:
        buf = np.ones(n, dtype=bool)
        buf.flags.writeable = False
        _BUFFERS[k] = buf
    return buf[:n]


def shares_buffer(a):
    """Whether the array a is a slice of one of the shared buffers (for the tests)."""
    return any(np.shares_memory(a, b) for b in _BUFFERS.values())

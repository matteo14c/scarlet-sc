---
name: A call looks wrong
about: SCARLET called (or refused) cells in a way you did not expect
title: "[call] "
labels: calls
---

**What you ran**
The command or Python call (learn / assign / transfer), the SCARLET version (`scarlet --version` or
`scarlet.__version__`) and the profile, if you set one.

**The data**
Technology (e.g. 10x, Smart-seq2), number of cells and genes, the label column, and whether the data are raw counts.
If you can share a small file that shows the problem, please say so (do not attach data you are not allowed to share).

**The cells that look wrong**
For a few of them, please paste:
- the call: `scarlet_best`, `scarlet_unique` (and, after transfer, the first-pass calls);
- the **reason** and the **flags** from the report's bundle (`bundle/tables/cells.tsv` after assign,
  `bundle/tables/query_cells.tsv` after transfer): `reason_best`, `reason_unique`, flip, low positive evidence,
  not refined;
- what you expected, and why (e.g. known markers, an independent annotation).

**The pre-flight summary**
The block SCARLET prints before it runs (technology profile and how it was resolved, marker sizes, blocklist).

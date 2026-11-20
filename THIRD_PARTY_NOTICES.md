# Third-party method attribution

`R/gbs_tassel_native.R` implements a native R, biallelic adaptation based on
TASSEL 5.2.96's `LDKNNiImputationPlugin`, `LinkageDisequilibrium`, and `Tuple`
sources. The algorithm is attributed upstream to Daniel Money and the TASSEL
conversion to Ed Buckler and contributors. This project's R implementation,
reference/query interface, memory guards and deterministic LD tie policy were
added in 2026. The implementation is marked `LGPL-2.1-only`, consistent with
the version declared in the upstream release metadata.

`R/gbs_tassel_search.R` is an LGPL-2.1-only extension of that adaptation. It
reuses the same donor distances and voting rules across neighbour counts,
and adds protected masked-call calibration and bounded coarse/refined search.

`R/gbs_preprocessing.R` wraps those LGPL components with separate predictor and
analysis marker pools, training-partition calibration/scaling and explicit
residual handling. It is also marked `LGPL-2.1-only`.

`R/gbs_parallel.R` is an LGPL-2.1-only target-block execution of the same native
LD, distance, tie and voting rules, with global CPU/RAM budgeting and checkpoints.
It retains the complete predictor/reference context for each target block.

Upstream release and source provenance:

- https://repo.maven.apache.org/maven2/net/maizegenetics/tassel/5.2.96/
- https://www.maizegenetics.net/tassel
- LGPL 2.1 text: `licenses/LGPL-2.1.txt`

The implementation is not affiliated with or endorsed by TASSEL. Numerical
equivalence has not yet been verified. The original project's other files do
not acquire a new license declaration through this notice.

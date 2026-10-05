# Cutoff-free positive clock compression

A constructive replacement for the eta^-2 two-edge quadrature factor. Status: independently audited PASS for the fixed-family coefficient/node theorem.

Main theorem: `CUTOFF-FREE-SECTOR-GAUSS.md`.
Implementation: `cutoff_free_sector_gauss.py`; exact rational Legendre core: `gauss_rational_core.py`.
Full execution contract: `COMPLETE-EXECUTION-BILL.md`.
Author checks: 1,356 passed in `checks.json` and `check-run.txt`. Independent audit: 3,175 checks passed; see `independent-audit/INDEPENDENT-COMPLEX-OU-CLOCK-AUDIT.md` and `REVIEW-PIN.json`. A full production scalar preparation and C0-C4 original-VALUE callback smoke test also passed; no native sampler was instantiated.

The actual three-copy integrand factors as a product of three Hilbert-valued complex OU transforms. An elementary L3 energy estimate gives local analytic radii proportional to the sector gaps. Positive dyadic Gauss then uses

    Q_tree=2(128 J m)^2,
    J=ceil(log2(32M/epsilon)), m=ceil(log4(32M/epsilon)),

where M is a complete original-target envelope with polynomial D and eta^-1 dependence. Both enter the node count only through logarithms. Main executed primitive order stays C3 for triple cubic and C4 for full/diagonal 3/3/6.

This removes the new quadrature cutoff power, including in arbitrary finite physical dimension. It preserves the original bank-mean target, literal positive weights, full cross-history census, exact original heat, proper-cut accuracy, all actual callers and the separate native replay bill. It does not remove mixed caller losses, old cutoff debts, native readout factors, or inherited work exponents.

The displayed constants are conservative and the certified rule can be very large; this is an asymptotic constructive improvement, not a claim that these raw node counts are practical. A caller may select the smaller of this fully certified rule and the sealed prior rule after evaluating both valid cost certificates.

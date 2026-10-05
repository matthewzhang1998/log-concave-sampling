# Nonlinear ancestry bridge repair

2026-10-05. A finite positive original-VALUE correction to the matrix-free resummed canonical-m3 graph.

## What is new

The original graph applies g to the averaged Gaussian history v. For nonlinear g this misses the average of g along the history and creates an order-A^2 bias. The new graph retains the original shift

    S=x-g(v)+g(g(w))

and adds centered VALUE finite differences about S. A known four-root Gaussian bridge retains the correct joint law of each nonlinear ancestor with x,v,w. Positive O(log^2(1/delta)) time quadrature matches the first nonlinear ancestry term without a path discretization, source-dependent matrix root, Hessian producer, or inverse-A sample bank.

## Results and limits

- Quantitative bias O(A^3 sqrt(D)) for analytically blockwise sources with bounded block dimension and Lip(Dg_block)<=B A. The algorithm need not know the block basis.
- This includes the sine source that provably leaves order-A^2 bias in the old graph, plus genuinely noncommuting nonlinear two-dimensional blocks.
- Exact canonical m3 mean for every anisotropic matrix quadratic.
- For any fixed dimension and fixed C2 potential shape g_A=A h, qualitative bias o(A^2) without the block/curvature-Lipschitz qualification. No dimension-uniform rate is asserted in this general-C2 corollary.
- Four private D-Gaussian roots; residual cost (5+3J)N_out original VALUES per fully charged replay; normalized private first <=4A and curl <=12A^2; actual caller origins, normalized callers and absolute floors exported.
- The same guarded finite own-mean consumer contributes only O(A^4 sqrt(D)). Its native guards and every replay/fill/clock/precision cost remain imposed.
- Not full order four, not growing-order closure, and not an unrestricted high-dimensional C2 rate. The new stencil does not automatically inherit the old five-residual certificate. Fixed structural certificate selection can retain the old graph on its original grade-four class.

## Files

- NONLINEAR-ANCESTRY-BRIDGE-REPAIR.md: complete construction, bias theorem, source ports, caller and executable ledger.
- bridge-quadrature-proof.md: independent Gaussian bridge, exact complex-disk certificate, and positive logarithmic-node operator quadrature proof.
- INDEPENDENT-REPAIRED-SOURCE-AUDIT.md: independent cancellation, constants, first/curl, energy/origin and numerical-floor audit.
- verify_complex_disk.py: exact-rational complex-contraction certificates.
- check_bridge_source.py / .json: 4,787 author diagnostics of the literal source and envelopes.
- independent_repaired_source_checks.py: independent exact radius/curl constants and 100 noncommuting shared-root derivative checks.

The numerical programs do not implement or re-certify the imported completed mean compilers. Real-arithmetic quadrature existence and operator error are proved analytically; finite-precision coefficient errors are separate absolute floors.

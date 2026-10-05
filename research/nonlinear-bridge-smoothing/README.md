# C2 extension of the stage-two nonlinear bridge

Read C2-STAGE-TWO-SMOOTHING-TRANSFER.md first.

## Established result

The unchanged finite original-VALUE stage-two graph has canonical-m3 mean-law error

    <= 110 A^(8/3) b^(2/3) sqrt(D)
       + original quadrature/completion/floor allowances,

under only the anchored C2-potential Hessian interval on orthogonal blocks of size at most b, A<=1/36, and all actual native compiler guards. The block basis is not queried. There is no bound on a Hessian modulus and no need to fix the source shape as A changes.

Smoothing is used only for an analytical comparison. Both the genuine nonlinear target remainder and the executed residual stencil have an extra-A weak source stability bound, proved using near-identity Gaussian entropy and their exact baseline cancellation. No smoothed source, additional Gaussian root, or convolution node is executed. The original 5D private roots and complete VALUE/HVP/caller/replay bill are unchanged.

For unrestricted dimension use b=D, giving 110 A^(8/3)D^(7/6), capped by the separate 5A^2sqrt(D) exact-outer bound. This is not a dimension-uniform improved power at natural sqrt(D) scale. It is one canonical-m3 component, not a full endpoint join or all-order recurrence.

## Files

- C2-STAGE-TWO-SMOOTHING-TRANSFER.md: integrated theorem, full cost, outer/completion interfaces, and exact scope.
- BASELINE-SMOOTHING-TRANSFER.md: detailed Gaussian entropy, true-history, residual-stencil, and analytical-mollification proofs.
- FINITE-SMOOTHING-IMPLEMENTATION.md: optional finite positive convolution implementation and its fully charged costs. Not needed by the main theorem.
- INDEPENDENT-SMOOTHING-AUDIT.md: independent review.
- check_smoothing_constants.py / .json: author elementary certificates and parameter checks.
- independent_smoothing_checks.py / .json: independent noncommuting/entropy/stencil/constant diagnostics.
- MANIFEST.json and SHA256SUMS: source pins and package verification.

## Reproduce checks

    python check_smoothing_constants.py
    python independent_smoothing_checks.py
    sha256sum -c SHA256SUMS

Diagnostics support the written arguments; they do not implement or re-prove the imported native mean compilers. The main new assertions are the analytical transfer estimates. Sealed stage-one and stage-two inputs are unchanged.

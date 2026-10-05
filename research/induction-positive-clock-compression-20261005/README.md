# Positive clock compression: a concrete improvement with a retained cutoff cost

The main result is in `POSITIVE-SECTOR-GAUSS.md`. It gives a positive, deterministic bridge rule for the complete triple-cubic and 3/3/6 families, including their literal old histories and caller ledgers.

For a two-edge BKAR tree,

    Q <= C [eta^-1 + log(B_*/epsilon)]^2 log^2(B_*/epsilon).

The sealed refined midpoint rule for three cubics had

    Q <= C [epsilon^-1 L^17 eta^-1 + log(1/eta)]^2.

Thus the new rule removes polynomial inverse-accuracy work. It does not remove the eta^-2 cutoff work. If eta>=c alpha^q and epsilon=alpha^p at fixed q, the node exponent improves from 2p+2q to 2q, up to the separately stated old-clock factors. Consequently this specific quadrature contribution is sublinear in requested accuracy grade p. Any q=q(p), readout dilution, native gap loss, and complete old/replay work must still be inserted into the true induction ledger.

Mechanism:

1. Explicit local Gaussian all-cut/HS bounds have factorial growth, giving a certified eta-scale analytic radius under Price differentiation.
2. Two positive Duffy sector maps remove the covariance min-path walls.
3. Positive high-order Gauss rules on subdivided dyadic gap panels have exponential degree accuracy.
4. Their literal weights satisfy the same bridge-gap mass majorants, so the sealed main, complete first, and root-square budgets survive with no node-count penalty.
5. Original private heat and the exact retained-bank node realization stay unchanged. Every numerical floor and replay remains charged.

Files:

* `POSITIVE-SECTOR-GAUSS.md`: theorem, proof, complete caller/private-heat qualification, rational precision construction, and work bill.
* `positive_sector_gauss.py`: executable exact-rational scalar certificates and lazy sector-node builder. It does not execute a native sampler or query original tensors.
* `checks.json`, `check-run.txt`: deterministic checks and conservative sample node censuses.
* `independent-audit/`: separately written audit and independent diagnostics.
* `INPUT-PINS.json`, `MANIFEST.json`, `SHA256SUMS`: exact input/output integrity records.

There is no negative-weight cubature, Monte Carlo N^-1/2 accuracy claim, target reheating, conditional-law replacement from marginal LAW, or assertion that the whole any-order sampler is finished.

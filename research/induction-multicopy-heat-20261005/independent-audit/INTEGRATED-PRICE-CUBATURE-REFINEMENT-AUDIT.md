# Independent audit: integrated-majorant cubature

2026-10-05. **PASS for the stated normalized-coefficient cubature certificate.** This separate addendum does not alter or unpin the four previously frozen proof snapshots.

Reviewed file: `INTEGRATED-MAJORANT-CUBATURE.md`.
SHA256: `3bdb992544dbba0af0581d0ae4cfc023e2528277fe144fb1b617fd43a8de5e56`.
A byte-identical `REVIEWED-INTEGRATED-MAJORANT-CUBATURE.md` snapshot and `integrated_majorant_cubature_review_pin.json` are stored beside this audit.

## Checks

1. **The correct positive envelope is used.** B is formed after summing the positive original-clock weights. Its local factors are (delta_i+eta)^(-a_i), with total exponent n/2 for the two Price hits. A raw old-node inverse-width product can have a larger variation exponent, so this order of operations is necessary. Including ALL potentially active Price pairs, with no sector-dependent zero indicators, supplies one envelope valid across every ordering wall.

2. **Price adds an ordinary graph edge, not a self-trace.** Its two hits lie in distinct replicas and therefore at distinct force occurrences. The augmented graph retains the original marked spanning tree and has no self-loop. Applying the graph theorem directly gives the one-Hilbert bound. One must not instead contract a general externally differentiated tensor using a trace estimate. The printed proof makes the correct choice.

3. **The q=2 integrated estimate has exponent n-2.** Each Price term has total local hit count 2(n-1)+2, and each replica degree is positive. Thus sum a_i=n/2. The audited all-n theorem gives e(n,2)=n-2. Its scalar proof applies to the direct augmented graph; the same local width exponents and ordered-gap majorants are available. The added finite Price-pair multiplicity is explicitly charged if not already in the q-hit convention.

4. **Dyadic comparison includes the last panel.** In a product panel, every s_e+eta changes by at most a factor two, including s_e in [0,eta]. Since delta_i+eta is the minimum of the incident regularized gaps, it has the same comparison. Each positive term and therefore their sum B varies by at most 2^(n/2). Edge ordering may change inside the cell without affecting this comparison.

5. **The straight-segment bound avoids an extra edge-direction factor.** Along t(lambda)=t0+lambda(t-t0), each min-path covariance is the minimum of finitely many affine functions. At almost every lambda its slope has absolute value at most ||t-t0||_infinity. Price differentiation therefore bounds the segment derivative by ||t-t0||_infinity times the ALL-pair envelope B, already summing all pair terms. No additional sum over edge directions is needed. Private heat and continuity justify sector-wall passage and singular covariance regularization under the imported contracts.

6. **The mesh and count follow.** On a midpoint cell of side length at most h, the maximum coordinate displacement is h/2. Integrating the segment estimate, applying the factor 2^(n/2) comparison, and summing cells gives exactly

       error <= (h/2) 2^(n/2) C(n,2) L^(6n-1)
                         eta^(-(n-2)) sqrt(D).

   The printed h is therefore sufficient. The existing bound 2+ceil(1/h)+ceil(log2(1/eta)) on each edge gives the stated (n-1)-fold tensor node count. All midpoint weights remain positive, satisfy u<=2s, and preserve panel masses. No tensor entry or expectation is evaluated to choose this rule.

7. **Normalization and remaining costs are correctly separated.** F and epsilon are normalized by alpha^(3n); the actual coefficient error carries that factor. Actual source, input-row, readout and fixed-rank constants must remain in the scalar coefficient and work bill when present. The count is still polynomial in inverse tolerance and cutoff. It does not establish a public-log bound or admission of the actual native sources at arbitrary grades, and does not supply a new retained/common-carrier ownership theorem.

## Conclusion

The refinement is valid with the two important envelope choices built into the final text: original-clock sums precede the majorant, and all possible Price pairs are included uniformly across sectors. The straight-segment argument justifies the sharper constant without an extra edge-direction factor. No further correction is required.

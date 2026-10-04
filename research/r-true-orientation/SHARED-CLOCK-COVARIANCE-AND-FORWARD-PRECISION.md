# Shared-clock covariance and forward precision

Publication copy: nonmathematical context and/or local paths were sanitized. Original and public SHA-256 values are recorded in `INVENTORY.json`; source/audit pins below identify their historical versions, not these edited bytes.

New supplement to the positive orientation clock67ca0e5c and affine fork consumer3b6fd06a. It specifies the stronger numerical tolerance needed when the old forward core is regenerated on the SAME clock.

Let J_j(X)=E Df(c_jX+v_jZ), with the common positive weights w_j of67ca. Define

    K_clock=sum_j w_j E[J_j J_j*],
    B_clock=sum_j w_j Sym E[J_j^2],
    O_clock=K_clock-B_clock.

Their infinite-clock targets are Cov(f), Sym B_f and O_f, respectively. All three use the SAME source, nodes, weights and original-input genealogy. They are analytical coefficients; no unknown matrix is queried.

## Spectral error

Let q_k=2 sum_j omega_j exp(-2k t_j), and suppose the clock obeys

    sup_(k>=1) |1-kq_k|/sqrt(k) <= delta.

For covariance, the error is E[(f-Ef) tensor sum_k (1-kq_k)f_k]. For every output row u,

    sum_k |1-kq_k|^2 E|u*f_k|^2
      <= delta^2 sum_k k E|u*f_k|^2
      = delta^2 E|u*Df|^2 <= delta^2 A^2 |u|^2.

One HS factor e gives

    ||Cov(f)-K_clock||HS <= delta A e.

The corresponding operator bound is delta A^2. Since B_clock=K_clock-O_clock, the historical orientation estimate gives

    ||Sym B_f-B_clock||HS <= delta(A+kappa)e,
    ||Sym B_f-B_clock||op <= delta A(A+kappa).

Thus the orientation-specific delta kappa e error is not also the generic covariance/forward error. No hidden extra curl factor is claimed for those targets.

## Required shared tolerance and versioning

To fit the fractional orientation budget b kappa e, choose the shared clock tolerance

    delta_clock <= c b kappa/(A+kappa),

with a fixed numerical allocation c. At kappa=A^(1+g), b=A^z this is O(A^(g+z)), rather than merely O(A^z). The positive clock still has only logarithmically many panels and degree. Its actual weighted source/path sums remain the same; no inverse-floor first loss is introduced by this smaller scalar tolerance.

Regenerate the complete old forward source on these nodes if it did not already use a compatible clock. Its main/covariance actions, first-response filters, coherent square, H/inner modules, numerical floors, root arrays, source versions and old-root mixture comparison must agree on that clock. The old covariance-stage hypothesis of3b6 is supplied only after those actual comparisons and with its mean independent of the owned roots.

This supplement establishes the common constant coefficient and precision identities. It is not permission to attach a new random-root fork to a previously averaged marginal old-core theorem. Any source-level old provider that is re-instantiated must retain its actual finite call/width/row guards and source-return proof.

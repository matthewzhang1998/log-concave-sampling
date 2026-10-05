# Independent audit: weighted heated queue and simultaneous family join

2026-10-05. Reviewed without editing the two primary notes. Exact input hashes and independent rational checks are in `independent-checks.json`. The pinned shared-variables manifest hash matches the advertised `875eca039776397c9f711ea71811e4236bcebac6b7411f430de617124592da3c`.

## Verdict

**The weighted potential identities and the Riesz energy-times-old-first improvement are valid under their stated hypotheses.** The rank-eight grades are exactly heat/old mixing **57/7**, native own **93/7**, new-bank self **99/7**, and root-caller first **13/2**. The improvement over 65/8 is **1/56**.

One displayed quantitative constant needs correction or explicit interpretation: JOIN §3 must charge **cross-group native return**, as the imported rank-eight theorem already does. The finite heated queue is a conditional grammar/termination result; it does not prove that an arbitrary actual mixed endpoint admits that grammar. Physical rank reductions are valid only for admissible contractions. These qualifications do not invalidate the one-stage 57/7 gain.

## 1. Weighted recurrence and termination

For `Psi=a-beta*N-gamma*(K+2)`, a bank tree with r arguments adds exactly `2(r-1)` local derivative slots and `r-1` edges. Therefore

    Psi' = sum_i Psi_i.

This is occurrence-level bookkeeping: repeated coefficient labels remain distinct vertices, repeated hits at a vertex increase its local k, and actual injection losses must be charged. The union of the input trees and the replica-tree bridges preserves a marked spanning tree. The positive Gaussian replica formula represents the coefficient only after its complete bank is averaged, as the note correctly states.

The heat-balanced allocation gives each nonroot an effective alpha exponent beta, while a genuine root gets beta+d. A conditional cluster with h genuine roots, including any subset of admissible side spines, has

    a' = beta*N' + gamma*K' + h*d,
    Psi' = h*Psi + 2*gamma*(h-1).

The argument does not require N or physical rank to increase. A captured-center path has exponent `d+ell*beta-gamma`; its worst one-force case is `beta+Psi+gamma`. Thus the stated sufficient first guard is algebraically correct. Its preservation concerns exponents, not automatically the numerical path sum or a previously executed version.

With every hittable source port having `Psi>=delta>0`, each genuine source-repair edge raises Psi by at least delta. For retained `G<P`, both Psi and N are bounded. Bank arity and root multiplicity are bounded by `P/delta`; tracing an occurrence through at most that many operations gives the proposed coarse quadratic bound on k. With finite jets, fixed q-lowering base branches, and certified finite quadrature at each node, this proves finite branching and finite depth.

Dropping a `G>=P` current is justified only by its actual one-Hilbert boundary-current/keep estimate. It is not justified by hypothetical descendants being above cutoff: G can decrease after side omissions. The note's boundary-error wording makes the correct choice.

### Structural restriction that should be explicit

WEIGHTED §6's physical Wick rule is an algebraic identity, **not an unrestricted certificate-preserving operation**. Only contractions in the admitted grammar may be used: no self-loop and no loss of the protected selected-endpoint marks on the spanning-tree leaves. Pairing arbitrary physical legs can fail both requirements. For example, tracing both physical slots of a one-vertex C0 tensor creates a self-loop/scalar trace; for the identity tensor this is D, not a dimension-safe one-Hilbert current. The imported conditional-grammar proof avoids precisely this operation.

Likewise, §7(2), requiring the *complete* source-current extraction to have the audited grammar, is a substantive hypothesis. The static tree theorem and coefficient cumulant formula do not alone establish it for every same-endpoint residual, observer or bank-Riesz branch. The result should remain named and used as a **conditional finite heated-queue lemma**, not an unconditional closure theorem. Keeping a bank-tree replacement's conditional coefficient in the ledger is necessary; future branches still need that source-current certificate.

## 2. Exact simultaneous Riesz estimate

For centered F, Gaussian integration by parts in the full B bank gives the displayed identity with

    A_t = (R_B F)(D_B W + t*D_B F)^T.

It is exact at the positive path endpoint, with no cumulant truncation and no extra coefficient derivative in the old term. Transfer the physical Hessian index through the independent keep. Conditional on A, the squared norm of `A^T G` is exactly `||A||_HS^2`, so no second Hilbert/dimension factor is introduced. Conditional expectation onto X_t supplies a continuity velocity with the same bound.

At fixed B,Y, degree-m hypercontractivity in H controls the Hilbert-valued Riesz field's L4 norm by its L2 norm. The assumed uniform-in-(B,Y) operator-Jacobian moments then allow conditional Hölder, followed by Riesz L2 contraction in B. Hence

    W2(X_1,X_0) <= 3^(m/2)/sqrt(kappa) * E*(L*alpha^s + J/2)*sqrt(D).

For the old term, independence of H and U even removes the hypercontractive factor; keeping it is conservative. There is no missing extra private inverse-width factor in this calculation.

The moment/ownership assumptions are essential and correctly visible: the whole stacked observer Jacobian must satisfy the old certificate, the bound must be uniform in B,Y rather than only jointly averaged, H must have B-independent geometry and be independent of U, and the keep must be independent of every coefficient/readout field. An exposed unattenuated B copy fails the hypotheses. Ordinary finite-second-moment/weak-differentiability conditions for the positive path and W2 are implicit and should be enforced in applications.

## 3. Native own-return constant needs the joint bound

JOIN §3 currently writes

    C_own = sum_h Lambda_native,h * [tau^K K_h]^2.

If `Lambda_native,h` means the **individual frozen-graph** constant, this does not follow from those ports: at a common endpoint the current differentiated in group h can hit group g's nonlinear reference/readout defect. The imported `RANK8-SMOOTHED-NATIVE-RETURN.md` §4 explicitly charges these cross-group branches by

    Gamma_mix * (sum_h Gamma_h*|rho_h|)^2
      <= Gamma_mix * n_groups * sum_h Gamma_h^2*rho_h^2.

Use that bound explicitly, or define each displayed `Lambda_native,h` to include the jointly verified factor `Gamma_mix*n_groups*Gamma_h^2`. The reference must remain the common-H transport of the **summed** tensor. Separate unrecoupled H_h comparators do not make this correction unnecessary.

This is a quantitative repair, not a change to 93/7, provided all added group/readout/moment factors are the declared fixed-rank/public-log quantities. An algebraic dependence on alpha, heat or dimension cannot be hidden in them.

## 4. Rank-eight arithmetic, floors and scope

At `(a,N,K,s,beta,gamma)=(8,8,6,1,1/14,1/7)`:

- Root exponent: `15/2-6/7=93/14`.
- Heat: `8+1/7=57/7`.
- Old mixing: `9-6/7=57/7`.
- Self mixing: `16-13/7=99/7`.
- Joint native own: twice the root exponent, `93/7`.
- Captured-root caller: `15/2-7/7=13/2`.

The changed clock and response tolerance is indeed proportional to alpha; the old alpha^(7/8) tolerance no longer meets this grade. With minimum width alpha^(1/7), the corresponding VALUE floor is proportional to `A*alpha^(8/7)`. Numerical first/Sobolev, filter, mean, calibration, covariance-gap and replay floors remain separate; a response tolerance cannot certify curl. The stated conservative starting native order and actual-prior check are compatible with beta=1/14.

The target heat estimate must apply to the complete summed rank-eight family. The proof correctly avoids assigning aggregate heat cancellation to individual histories and does not replace original shared banks by independent history banks. A fixed-rank/public-log **grade** statement must retain its logarithms; a strict log-free O(alpha^P) claim would additionally need bounded leading constants or a margin.

For the weighted floor-only restoration strategy, the ceiling `8+(8-8/14)/8=125/14` is correct. The C2 high-frequency example genuinely excludes obtaining a universal improved VALUE heat power from regularity alone. Neither a heat norm nor a finished LAW error is a source-qualified residual port. No all-order re-entry is proved or needed for the one-stage claim.

Finally, the expanded-DAG max-plus query recurrence is a legitimate accounting rule once each replication exponent is supplied. A finite quadrature's existence is not a public-log cost estimate. The actual `M_native` high-accuracy/order exponents and the complete restoration/version/replay costs remain missing; no `c(P)=o(P)` conclusion follows.

## Checks performed

- Verified the pinned manifest and the listed graph/frozen-port file hashes.
- Independently checked 3,510 exact-rational recurrence instances and every displayed rank-eight balance above.
- Read the marked-spanning-tree addendum and uniform same-endpoint current lemma to distinguish structural closure, scalar cumulants, and actual positive-current consumers.

No primary file was changed. Recommendation: accept the one-stage energy-times-old-first gain and conditional weighted termination lemma after the joint-own-constant clarification and explicit admissibility wording.

Subsequent audit updates: the author's revised notes addressed those clarifications. `COMMON-CARRIER-AUDIT-ADDENDUM.md` checks the separate raw-carrier proposal and identifies its missing carrier-compatible native hybrid. `FIXED-NATIVE-COST-SUPPLEMENT.md` locates an explicit fixed-order polylogarithmic native serial recurrence in the pinned LOW30 source; this refines the earlier statement that the local native exponent was unspecified. Full replay, cubature and growing-order closure remain separate obligations.

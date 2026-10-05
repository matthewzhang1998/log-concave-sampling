# Independent audit: finite-dimensional shared-carrier replacement

2026-10-05. The primary notes were not edited by this auditor. This report concerns the revised two-energy/operator-Jacobian taming lemma, not its earlier, looser coefficient-only Lipschitz estimate.

## Verdict

**The actual shared-carrier native-to-polynomial hybrid is proved under the explicitly imported fixed-order native inputs and original product-bank/read-set conditions.** The proof-only cap and radial projection give the missing global carrier-Lipschitz envelope at every mixed comparison stage. The marginal native error is used before the group's final independent keep; no Gaussian deconvolution or carrier-retaining native prior is assumed. The resulting admission window is genuinely finite-dimensional: the new factor grows approximately as `D^((m-1)/2)` at fixed polynomial degree. It is not a public-log-dimensionality theorem.

**The independent variance-partition alternative is also valid under its stated full native/frame inputs and share/census guards.** Its rank-nine inverse-share powers `a^(-9)`, `a^(-17)`, and `a^(-8)` are correct. Its conversion to public-log-cost, public-log-guard admission remains conditional on a complete admissible packet census, which has not been supplied by this replacement argument.

This is an audit of the new replacement implications against the pinned source interfaces, not a fresh verification of every theorem in LOW30, the original ideal common-carrier current theorem, or the full triple-cubic heat/cubature construction. Their numerical constants, exact finite compiler degree, original-source eligibility, complete observer ledger, and all old target debts remain real imported obligations.

## 1. Product-bank and source interfaces

Condition on retained external labels `E`, including the original coefficient bank and retained callers. The actual construction uses a Gaussian `H` independent of a bank of mutually independent packet tapes `Z_h`, conditionally on `E`. This is stronger than merely saying that the packets are conditionally independent given `H`: the fixed-law base tapes can be drawn independently of a coupling of two different carrier values. It follows from the displayed Gaussian coisometry and independent null/private tapes.

The Gaussian geometry is exact:

- `Q_h=L_h^T L_h/v` is an orthogonal projection and `Pi_h=I-Q_h`.
- `P_h=L_h^T H/v+Pi_h V_h` is standard Gaussian, with `L_h P_h=H` pointwise.
- The map from `(H/sqrt(v),V_h)` to `P_h` has operator norm one.
- The old component and final keep read no new packet tape. An old random component can be coupled identically along with `E`; it need not be a deterministic constant.

The pinned source is LOW30, SHA256 `7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8`.

Relevant interfaces checked directly:

1. The finite compiler theorem around lines 5703–5740 permits fixed-order constants and declared public logarithms, with original external labels frozen.
2. The prior interface at lines 5835–5855 retains incoming input and original labels, while discarding private tapes. This is an ordinary conditional law interface, not the missing carrier-retaining one.
3. The complete-spine tree replacement at lines 6029–6061 prices actual prior errors before adding the final independent readout buffer. Its proof keeps declared publics, collapses internal spines with no other observer, and uses internal terminal innovations.
4. Lines 6200–6206 put all primitives with nonzero source-zero linear readout in `P`; the remaining zero-readout roots are independent `U`. The polynomial comparison at lines 6210–6249 and its moment proof through approximately line 6300 give finite `P` degree, every required fixed value moment, and operator-first moments. The lack of a finite degree requirement in original-query `U` is harmless: those coordinates are coefficient labels in the taming lemma.
5. The raw-path proposition at lines 6488–6506 and its complete path proof through approximately line 6578 give global actual residual Lipschitz bounds. Varying original labels need their own genuine uniform source first; a matrix moment alone is not enough. The active radial-pullback source supplies the stated private/probe/caller bounds, with each caller injection and ancestor path charged.

The actual and polynomial residual budgets must be root-linear for the concrete completed packet. Merely knowing a selected formal tensor or an ordinary output law would not establish this. The primary note imports the complete path and polynomial budgets rather than attempting to derive them from a marginal law comparison.

### Why the local marginal error is available before the keep

This point is essential. A Wasserstein estimate only after convolution by a Gaussian cannot simply be stripped of that Gaussian.

Here the source supplies a different route: first couple the actual and ideal completed spines using the tree-replacement proof, before the final independent readout/buffer. Add the retained public readout within that coupling. Next compare the ideal forest and finite polynomial on their common Gaussian roots. In the latter strong comparison the same final `eta Z` appears on both sides and cancels identically from the vector difference. Triangle inequality therefore gives the required ordinary marginal error `e_h` between the two pre-keep packets. This uses neither a deconvolution inequality nor an actual-prior coupling with the carrier held fixed.

The later intrinsic current comparison does use the final keep. It is a separate step and is not included in `e_h` by removing that keep.

## 2. Hermite evaluation and matrix-Jacobian bound

Let `A=1+D+R^2`, `rho=1/A`, and `D>=1`. Mehler's positive series gives

    M(rho,x)=(1-rho^2)^(-D/2) exp(rho |x|^2/(1+rho)).

For `|x|<=R`, its logarithm is at most

    D/[2(A^2-1)] + R^2/(A+1) < 7/6.

Indeed `D<=A-1`, `A>=2`, and `R^2<A+1`. Thus `M<exp(7/6)<4`; the note's constant 4 is safe even without using `R>=sqrt(D)`. Positivity gives

    sum_{|alpha|<=j} phi_alpha(x)^2 <= 4 A^j.

The exact derivative sum is

    sum_{|alpha|<=m} |grad phi_alpha(x)|^2
      =sum_{|beta|<=m-1}(D+|beta|)phi_beta(x)^2,

so the stated upper bound `4(D+m) A^(m-1)` is valid. The revised proof correctly avoids spending this additional dimension factor in the Lipschitz constant.

For the matrix polynomial `J=D_x u^p=sum C_beta phi_beta`, matrix Parseval is an identity of positive semidefinite matrices:

    sum C_beta^T C_beta=E_x[J^T J].

For each unit vector `w`, ordinary Cauchy–Schwarz in the coefficient index gives

    |J(x)w|^2 <= (sum phi_beta(x)^2) E_x|Jw|^2
               <= 4 A^(m-1) E_x||J||op^2.

No Hilbert–Schmidt norm of `J` is introduced. This proves the sharpened pointwise operator estimate in the primary note.

## 3. The combined cap and strong taming bias

Write

    A0^2=E_x|u^p|^2,
    A1^2=E_x||D_x u^p||op^2,
    Astar=max(A0/sqrt(D),A1).

The combined twelfth-moment input is sufficient, exactly as written:

    E Astar^12
      <= D^(-6)E A0^12 + E A1^12
      <= D^(-6)E|u^p|^12 + E||D_x u^p||op^12
      <= s_h^12.

Enlarging two individually valid envelope constants by at most `2^(1/12)` supplies this combined norm. With `q=min(1,K s_h/Astar)`, both `q A0<=K s_h sqrt(D)` and `q A1<=K s_h` hold. Therefore projection onto the radius-`R` ball gives a globally carrier-Lipschitz comparator with

    Lip_H u^t <= 2K A^((m-1)/2) s_h/sqrt(v).

The projection is 1-Lipschitz and the polynomial Jacobian is bounded throughout the convex ball; differentiability at its boundary is unnecessary. For `m=0`, the residual is carrier-independent and the second term is omitted.

The exact paid bias calculations are:

- Cap: `||u^p-q u^p||2^2 <= D E[Astar^2 1{Astar>K s_h}] <= s_h^2 D K^(-10)`.
- Original tail: Holder with exponent 12 gives the probability power `5/12`; Gaussian concentration then gives `s_h sqrt(D) exp(-5 ell^2/24)`.
- Projected tail: the capped value is bounded everywhere on the ball by `2K s_h sqrt(D) A^(m/2)`, giving its tail factor `exp(-ell^2/4)`.

Thus `d_h<=delta s_h sqrt(D)` with exactly the note's three terms. The coupling uses the identical carrier and identical tapes. The displayed bias is an integrated `L2` bound over that bank, not a uniform conditional error for every carrier value. This distinction is preserved in the revised note.

## 4. Full mixed-environment disintegration

For a fixed retained label value, choose an optimal, or arbitrarily near-optimal, coupling of the completed pre-keep outputs `X_h^a` and `X_h^t`, at cost `epsilon_h<=e_h+d_h`. Disintegrate each packet's own joint law of carrier and output conditional on its output. Sample its carrier from that regular conditional law. This recovers the correct `(H_a,X_h^a)` and `(H_t,X_h^t)` marginals. Finite Euclidean Gaussian programs are standard Borel, so the needed kernels exist.

The pointwise relation `X=H+u` gives

    ||H_a-H_t||2
      <=epsilon_h+||u_h^a||2+||u_h^t||2
      <=epsilon_h+2s_h sqrt(D)+d_h.

No equality of the carriers is asserted. Their Gaussian covariance alone is not the reason for this bound; their correct marginal joint laws and the displayed identity are.

Now draw all other base tapes independently of that coupling, using exactly the same other tape on both sides. At any stage each unchanged packet is actual or tamed, so the entire unchanged residual has deterministic global carrier-Lipschitz constant at most `L S_-h`. Each side has the correct full joint hybrid law, including dependence through its own carrier. Consequently

    W2(stage_before,stage_after)
      <= (1+L S_-h)(e_h+d_h)
           + L S_-h(2s_h sqrt(D)+d_h).

Summing produces `(1+LS)E0+(1+2LS)sum d_h+2LS^2 sqrt(D)`. Restoring all original polynomial packets on a common original bank adds `sum d_h`. This is precisely the claimed equation (3). At no stage is an untamed polynomial residual used as a globally Lipschitz environment.

All statements may be conditioned on `E`, then integrated by conditional coupling and Minkowski. If an application replaces the deterministic uniform envelopes by random envelopes, its product-integrability bounds must be supplied; plain unweighted marginal averages are not a substitute.

## 5. Exact finite-dimensional guard and grade

The amplitude arithmetic is exact:

    17/2 + 8(1/16)=9,
    17/2 + 5/2=11,
    17 - 1/2=33/2,
    33/2 - 10=13/2.

For `0<alpha<=1`, take `K=alpha^(-1/2)` and the least positive integer `ell` satisfying the stated two-tail inequality. Existence follows because the exponential in `ell^2` dominates the fixed polynomial radius factor. At fixed `m`, `ell^2=O_m(log(2+D)+log(1/alpha))`; the ceiling changes only the constant. Hence `delta<=2alpha^(5/2)`.

For `m>=1`, the exact new coefficient is

    C_dim=2(1+D+(sqrt(D)+ell)^2)^((m-1)/2)/sqrt(v).

With `S<=B alpha^(17/2)`, `LS<=1`, and `E0<=c_e alpha^11 sqrt(D)`, the entire replacement error is at most

    [2 c_e alpha + 8 B alpha + 2 L B^2 alpha^7]
         alpha^10 sqrt(D).

Thus the guard `2 L B^2 alpha^7<=c_hybrid` is valid. Because `L` is a maximum, this means checking both the base contribution and the alpha-dependent contribution `2 C_dim B^2 alpha^(13/2)<=c_hybrid`. The latter alone is not a replacement for the full maximum in general. Native-radius, ideal-current, caller, and old-target/error-budget guards remain separate.

At each fixed finite `D`, fixed `m`, and fixed `v>0`, with `B` of the declared public-log size, sufficiently small alpha satisfies the new guards. Allowing `m` to grow with `log(1/alpha)` would require a separate analysis and is not covered by that statement. The proof's polynomial dimension factor is real, except for low-degree special cases: for `m=1` the sharpened Lipschitz coefficient has no polynomial `D` factor; for `m=0` no carrier taming is needed. Neither special case removes the general fixed-degree obstruction relevant to the higher-rank packets.

## 6. Independent variance partition

Take deterministic positive shares with `sum a_h^2=1` and independent complete packet banks conditional on the original retained labels. Then the summed physical carriers and summed keeps have covariances `vI` and `kappa I` and are independent. The redesign does not expose any packet tape externally.

Rank-nine scaling is checked by conditioning after the local native comparisons. For standard Gaussian carriers with correlation `a_h` to the aggregate carrier, homogeneous chaos of degree eight projects with factor `a_h^8`. Therefore multiplying the packet itself by `a_h` transports its rank-nine current by `a_h^9`. This verifies the root allocation `b_h a_h^(-9) alpha^(17/2)`.

The scaled native intrinsic square has the factor

    a_h (b_h a_h^(-9) alpha^(17/2))^2
       =b_h^2 a_h^(-17) alpha^17.

The local drift and complete first have factor `a_h*a_h^(-9)=a_h^(-8)`. Equal shares consequently give `J^(17/2)` in the native-square sum, `J^4` in the l1 drift/first budget, and `J^8` in its quadratic null-current row. The per-packet root guard has `J^(9/2)`.

### One-Hilbert null-bank reduction

After the orthogonal Gaussian rewrite, let `F=E_loc-E_sum`, centered in the whole null bank at fixed original labels and `H`. Then `D_Z F=D_Z E_loc`. Gaussian Riesz gives

    E[F dot grad phi(X_t)]
      =t E[(R_Z F)(D_Z F)^T : Hess phi(X_t)].

The fixed-degree/frame inputs give one vector-energy bound `||F||4<=C alpha^9 B_a sqrt(D)` and an operator-first bound `||D_Z F||4,op<=C alpha^9 B_a`. The Hilbert-valued Riesz estimate gives `||R_Z F||4,HS<=C||F||4`. Hence

    ||(R_Z F)(D_Z F)^T||2,HS
      <= C alpha^18 B_a^2 sqrt(D).

The matrix coefficient reads no keep. Integration by parts through the independent `K~N(0,kappa I)` turns its Hessian current into a velocity whose `L2` norm is at most its HS norm divided by `sqrt(kappa)`. This follows by conditioning on the matrix and using `E|M K|^2=kappa ||M||HS^2`; it does not use a second vector energy. Integrating `t` along the positive path proves the claimed one-Hilbert Wasserstein row. Gamma must include the actual matrix derivative frames, as explicitly required in the note.

The projected target tensor is exactly the original sum; its original-bank first has no inverse-share factor. The actual executed packet derivatives still have the conservative l1 share inflation. The note correctly distinguishes these two statements and retains the complete terminal chain rule, numerical Sobolev floors, and old residual first.

### Census limitation remains substantive

The full list must include ordered old histories, endpoint panels, bridge nodes, force hits, permutations, and every required native version. Public-log total coefficient mass is not by itself a public-log bound on this list. If the count is algebraic in alpha, its powers must be substituted into the share guards and error rows. The midpoint fallback does not establish the missing public-log census. No all-order or growing-native-order conclusion follows.

The final integration addendum cites `POSITIVE-SECTOR-GAUSS.md`. I checked its displayed node bound and independent audit. Its new two-edge count is `2 m_cub^2[ceil(32P/eta)+J_cub]^2`, with the degree/panel factors logarithmic in the requested accuracy and coefficient envelope. At fixed original `eta>0` and otherwise public-log parameters, that new bridge count is public-log. At `eta>=c alpha^q_eta`, it has safe new exponent `2q_eta`. Its positive-rule majorant argument preserves the summed cubic coefficient and first budgets without an additional node-count factor.

Thus, with all remaining census/frame constants public-log, substitution in the equal-share rows gives exactly

    native grade = 17 - (17/2)(2q_eta) = 17 - 17q_eta,
    null-bank grade = 18 - 8(2q_eta) = 18 - 16q_eta.

Strictly exceeding grade 10 requires `q_eta<7/17` and `q_eta<1/2`, respectively; the former is stronger. Strictness leaves room for the declared logarithmic factors. A nonzero old-census exponent `nu_old` instead gives grades `17-17q_eta-(17/2)nu_old` and `18-16q_eta-8nu_old`. Root, caller, native gap, source-version, and replay conditions still require their actual ledger values. This arithmetic validates the added sufficient window; it does not prove that the inherited old cutoff or complete old census satisfies it.


## 7. Independent reproducible checks and exact scope

`check_independent_carrier_audit.py` and `independent-check-results.json` accompany this report. They perform:

- 480 direct multivariate Hermite kernel/gradient checks and the exact derivative-index identity;
- matrix-valued Parseval by exact-degree Gaussian quadrature and 100 matrix evaluation checks;
- exact symbolic conditional-Hermite identities at degrees 0 through 9;
- an exact nonlinear conditional null-bank Riesz identity;
- exact rational alpha/share-exponent checks;
- enumeration of all 243 labelled force-hit histories, with maximum 44 raw VALUE calls;
- twelve illustrative literal finite-D radius/guard computations, including `D=10^6` and `m=16`.

All checks pass. Numerical examples supplement the proof and certify neither actual application constants nor absent compiler/cubature inputs. The JSON records SHA256 hashes of the audited primary notes; rerunning the checker refreshes those hashes after any primary-note edit.

No new sampler operation, source tensor, conditional expectation, Hermite coefficient, matrix Jacobian, or projection is executed by the shared-carrier repair. It changes only the coupling proof. It does not remove source/calibration/mean/arithmetic/Sobolev floors, original target-heat debts, complete replay costs, retained-observer requirements, or the distinction between private source curl and the actual completed terminal curl.

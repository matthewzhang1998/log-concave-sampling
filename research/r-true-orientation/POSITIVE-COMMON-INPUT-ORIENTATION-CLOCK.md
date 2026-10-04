# A finite positive common-input orientation clock

Publication copy: nonmathematical context and/or local paths were sanitized. Original and public SHA-256 values are recorded in `INVENTORY.json`; source/audit pins below identify their historical versions, not these edited bytes.

## 0. New proof and exact scope

This is a new proof, not a reconstruction. It supplies a finite positive clock for the constant true-orientation coefficient of the WHOLE same source. Its Gaussian input genealogy must be preserved at every node. It does not independently smooth primitive Hessians or by itself execute a coefficient source.

Let f:R^n->R^n be the same finite source, with centered Gaussian energy e, full first A and curl C=Df-Df* satisfying ||C||op<=kappa. Let P_t be whole-input OU heat and Aop=I-E-T the output-swap defect from the exact historical whole-source heat lemma (4b2b72a4).

Define

    Gamma_t = -Sym E[(P_t Df)(P_t C)],
    O_f = 2 integral_0^infinity exp(-2t) Gamma_t dt.         (1)

The expectation in Gamma uses one common Gaussian root X; the two OU innovations may be independent conditional X. Within each evaluation of f, all original aliases, finite twins and query records remain on one substituted input. Set c_t=exp(-t), v_t=sqrt(1-exp(-2t)).

For any 0<delta<=1/4, the construction below gives positive weights w_j and nodes t_j with

    || O_f-sum_j w_j Gamma_(t_j) ||HS <= delta kappa e,    (2)
    || O_f-sum_j w_j Gamma_(t_j) ||op <= delta A kappa.    (3)

It has O(log(1/delta)[log(1/delta)+log log(1/delta)]) nodes. Scalar precision is logarithmic; complete source Gaussian-row encoding has its separate assigned numerical budget.

## 1. The spectral error criterion

Write f=sum_(k>=0) f_k by vector Gaussian chaos. The exact identity is

    k Aop f_k = Gaussian-divergence(Curl f_k),
    O_k=Sym E[f_k tensor Aop f_k]
       =-(1/k) Sym E[(Df_k)(Curl f_k)].

Hence (1) follows from integral_0^infinity 2exp(-2kt)dt=1/k. For a quadrature q_k approximating this last scalar integral, the coefficient error is

    Sym E[(f-Ef) tensor sum_(k>=1)(1-k q_k) Aop f_k].

The divergence row frame and orthogonality imply the dimension-free criterion

    error_HS <= kappa e sup_(k>=1) |1-k q_k|/sqrt(k),
    error_op <= A kappa sup_(k>=1) |1-k q_k|/sqrt(k).       (4)

This is one Hilbert energy times one ordinary row frame. It does not estimate a sampled outer product by two Hilbert energies.

## 2. Truncation and positive Gaussian quadrature

Take

    t0=delta^2/128,
    T=.5 log(16/delta),
    J=ceil(log_2(T/t0)),    a_l=2^l t0,
    Tprime=2^J t0 >= T.

Use the J dyadic panels [a_l,2a_l], l=0,...,J-1. On every panel use the positive M-node Gauss-Legendre rule, where

    M >= ceil(log_4(384/delta)).

Let its nodes and ordinary integration weights be t_j, omega_j, and set

    w_j=2 omega_j exp(-2t_j)>0.

The corresponding scalar multiplier is q_k=2 sum_j omega_j exp(-2k t_j).

The omitted initial interval contributes

    sup_k (1-exp(-2kt0))/sqrt(k) <= sqrt(2t0)=delta/8.

The omitted final interval contributes at most exp(-2Tprime)<=delta/16.

For an individual panel [a,2a], map the Bernstein ellipse of parameter 2 from [-1,1]. Its real part is at least 7a/8. Therefore the analytic function 2exp(-2kt) is bounded there by 2exp(-7ka/4). Truncating its Chebyshev series at degree 2M-1 has uniform error at most 4B 4^-M. Exactness of the M-point rule for that polynomial and positivity/mass a bound the quadrature error by 8Ba4^-M. Thus

    |panel scalar error| <= 16 a exp(-7ka/4) 4^-M.

Put x_l=k a_l. The dyadic sum sum_l x_l exp(-7x_l/4) is at most 3: below one use the geometric sum, and above one use the rapidly decreasing dyadic exponential tail. Multiplying the total scalar error by sqrt(k) gives at most 48 4^-M/sqrt(k)<=delta/8. Together the ideal truncation and quadrature errors are less than delta/2 in (4), leaving room for finite scalar encoding.

No analyticity of a finite VALUE-error field is invoked. Only the exact scalar OU spectral multipliers are approximated.

## 3. Positive masses and actual complete path sums

All nodes stay inside their panels. With v_j^2=1-exp(-2t_j), positivity gives fixed constants such that

    sum_j w_j <= C,
    w_j <= 3 v_j^2,
    sum_j w_j/v_j <= C,
    sum_j w_j/v_j^2 <= J.                              (5)

For the last row, each panel contributes at most 2a/(exp(2a)-1)<=1. For the third row use a geometric sqrt(a) sum below a=1 and exponential decay above one. The individual-weight row follows from omega_j<=a and the same panel comparison.

Also

    sum_j sqrt(w_j)/v_j <= J sqrt(M),                   (6)

by Cauchy-Schwarz and (5). Thus independent per-clock banks combined with readouts sqrt(w_j) have complete root first O(sqrt(J)) times their unit-width first, while shared callers/absolute outgoing masses cost at most Jsqrt(M). Those are actual logarithmic path bills. No per-node inverse t0 power is needed after this positive weighted aggregation.

When a positive fork is run at each node, its target covariance coefficient is multiplied by w_j and its complete output by sqrt(w_j). Its variance share is therefore also multiplied by w_j. A numerical common variance per node and one final known fill can be chosen because sum w_j is bounded. The target coefficient and every readout must be scaled together; a law to a saved unknown covariance matrix is not an executed reserve.

## 4. Scalar rounding and complete Gaussian-row encoding

Use the ideal construction above for an error budget delta/2, or retain its displayed slack. Isolate every node within its original panel with |Delta t_j|<=eta a and positive weights renormalized to exact panel mass a, with per-panel total weight error <=eta a. For sufficiently small fixed eta=C0 delta, the extra scalar spectral error is at most delta/4.

Indeed node motion, after multiplying by sqrt(k), is bounded by a constant times eta/sqrt(k) sum_l (k a_l)^2 exp(-c k a_l). Weight error is bounded by a constant times eta/sqrt(k) sum_l (k a_l) exp(-c k a_l). Both dyadic sums are uniformly bounded. Positive node/weight isolation has logarithmic bit precision in delta and t0, and does not change the node count. Exponential weight encoding can be included in this positive relative-weight budget.

The actual source at the rounded node must use its matching c=exp(-t), v=sqrt(1-c^2). If these rows are themselves encoded numerically, preserve the exact known covariance normalization or pay the complete source VALUE comparison under the actual row perturbation. Its budget uses the actual full input dimension, source first, readouts, caller profiles and finite path amplification. Such a source comparison is numerical; it is not obtained by differentiating the quadrature error. No uncharged precision or discarded replay is claimed.

## 5. Source and observer boundary

For the marker at node j, the conditional source

    f_(X,j)(z)=[f(c_j X+v_j z)-f(c_j X)]/v_j

has first A and the same curl bound. Its centered conditional L2 energy e_j(X) satisfies ||e_j||L2(X)<=e/v_j by total variance. The arbitrary constant subtraction does not change the selected derivative or the centered energy. Actual finite pair contrast/translation identities must be used when executing this anchored source.

The common conditional Jacobian is E_Z Df(c_j X+v_jZ). Any skew-channel source at that node must match its curl at the SAME X and the SAME whole-input covariance. The affine-state pullback construction supplies such a match for its declared known-map class; independent primitive additive heat does not.

Equations (2)--(3) concern the averaged coefficient. If node roots X_j appear in the output covariance before they are integrated, either cancel the same-root fields at the same endpoint or supply their actual root/Price currents. Positive scalar weights and (5) alone do not control a random covariance mixture. The finite source work at every node, and all shared/owned labels in that comparison, remain explicit.

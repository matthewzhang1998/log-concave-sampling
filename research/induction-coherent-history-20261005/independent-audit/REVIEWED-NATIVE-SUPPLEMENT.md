# A retained-carrier native port for the first marked pair-covariance correction

2026-10-05. This supplements CURRENT-CARRYING-RANK-TWO-RESUMMATION.md. Unlike merely writing a likelihood identity, it gives a bounded finite original-VALUE source, its complete/caller first and curl, explicit clipping/truncation errors, a literal zero-baseline native-mean substitution, and the correct retained-carrier current contract. It does not supply the genuine history base/singleton coefficients.

## 1. Retune the existing actual pair source before differentiating it

In TRUE-PAIR-HISTORY-RESUMMATION.md, keep all levels, geometry, smoothing roots, shared midpoint roots, caps and positive weights. Change only its secant scale from

    delta_k=e_k A^2 sqrt(tanh(d_k/2))

to

    delta_k=e_k sqrt(tanh(d_k/2)).                       (1.1)

The complete source is still exactly fourteen original VALUES for two nuisance replicas. Its uniform conditional-inner first remains ell<=C A^2, since that proof is independent of delta. Every branch still uses the full original nonlinear shifted VALUES and its saved same-bank anchors.

The normalized scale now satisfies a_k delta_k/e_k<=C, instead of O(A^2). The original smoothed secant proof therefore gives mean-field bias C A^3 n and covariance bias

    C[A^5 n^(3/2)+A^6 n^2]                             (1.2)

in one physical block of dimension n. In a supplied block decomposition n<=b, this is C[A^5 b+A^6 b^(3/2)]sqrt(D). It is still smaller than the original fixed-b A^(38/9) actual-history covariance error. No unsmoothed Taylor expansion is made.

The frozen-label COMPLETE source first improves by A^2. The original bound C A^(-3/2)d_K^(-5/4) becomes

    L_F<=C A^(1/2)d_K^(-5/4)<=C A^(-8/9),             (1.3)

for d_K comparable to A^(10/9), with the same declared coarse-scale constants. Radial cap projection does not increase this first. The source first is not inferred from its A^2 energy or conditional-inner radius. This parameter change must be recorded in every source version and replay.

## 2. The finite bounded current source

For one physical block retain Y and a NEW standard Gaussian G in R^n as CALLER LABELS. They are not private roots of the correction compiler. Let S contain all sampled level/cell selectors and local-time/geometry clock labels. Retain S as a frozen mode through each native call. Let Omega be the complete GAUSSIAN source bank generating the retuned capped F,F', and a coherent finite original-VALUE mark M in R^n at that S. M and H may have arbitrary common ancestry. The entire (S,Omega) bank is independent of G conditional on Y. The uniform first bounds below differentiate Gaussian/private and declared continuous caller coordinates with S frozen; no global first across inverse-CDF selector boundaries is claimed.

Set H,B,K,E,h,e as in the companion, so E=H/v+H^2/(4v^2), h=R^2/(2v)<=1/8, e=2h+h^2. Choose deterministic caps R_M,B_G>0 and execute

    Mbar=projection_R_M(M), Gbar=projection_B_G(G),
    Q=(Gbar^T E Gbar-tr E)/2,
    N(Y,G;Omega)=Mbar Q.                                (2.1)

No matrix is materialized. Compute HG, H^2G, tr H=F dot F', and

    tr H^2=[|F|^2 |F'|^2+(F dot F')^2]/2.

Thus Q uses O(n) arithmetic. The source is a finite original-VALUE graph with every H/M ancestor charged. It is distinct from a derivative/HVP producer. The projections are Lipschitz; where a C1 native interface requires a smooth cap, use a fixed smooth radial cap with the same support/radius and fixed-factor Lipschitz bounds. The constants below enlarge by those fixed factors, and the same tail proof applies.

Suppose the uniform complete/caller firsts of the mark and replica vectors are at most L_M,L_F, respectively. Then valid matrix-source first bounds are

    L_H<=2R L_F,
    L_E<=(1+h)L_H/v.                                    (2.2)

For a private/caller direction, the image of DE is contained in span{F,F',DF,DF'}, hence rank(DE)<=4 and |tr DE|<=4||DE||op. This avoids inserting a full ambient-dimension trace factor.

The complete private first, and a safe combined private/retained-(Y,G) first, are bounded by

    L_N <= L_M e(B_G^2+2)/2
           +R_M L_E(B_G^2+4)/2
           +R_M e B_G.                                 (2.3)

The last term is the retained-G path; it is absent from a purely private Omega derivative. All caller versions of L_M,L_F enter the same formula. v and the geometry labels are frozen here; a changing v has its extra explicit paths and is not hidden by (2.3).

A deterministic source-energy bound is

    |N|<=E_N:=R_M e(B_G^2+2)/2.                         (2.4)

Any centered L2 energy is at most E_N. Embed the output with a fixed known coisometry P from the complete private bank (pad with an unused n-Gaussian row if needed). Then the entire square lift P* N has

    first<=L_N, curl<=2L_N.                             (2.5)

This is a conservative generic curl, not a claimed gradient property. Its whole-bank lift and all off-row derivatives are included.

## 3. Paid clipping and current truncation

Let Q0=(G^T E G-tr E)/2. Conditional on Omega,Y, Gaussian isometry gives E_G Q0^2=tr(E^2)/2<=e^2. Therefore

    ||M Q0-Mbar Q||_2
       <=e||M-Mbar||_2
          +(e R_M/2)(E[|G|^4 1_{|G|>B_G}])^(1/2).     (3.1)

All norms use the actual retained Y law. For a local source with |M|<=c A^3|x0| and Gaussian affine x0 of covariance at most c0 I under standard Y, choose

    R_M=C A^3 sqrt(n+L), B_G=sqrt(8(n+L)), L>=1.

Elementary Gaussian tails make (3.1) at most C A^3 e (n+L)^(3/2) exp(-c L), with constants depending only on the declared affine rows. This is an integrated-standard-Y restoration, not an arbitrary-caller claim. The executed capped source and its firsts are uniform in Y regardless.

The exact likelihood transfer and its p=1 Hermite tail give, for every bounded measurable vector callback a with ||a||infinity<=1,

    |E[M dot a(Y,Z/sqrt(v))]
        -E[(M+M Q0) dot a(Y,G)]|
       <=2sqrt(2) e^2 ||M||_2.                         (3.2)

The M in the second line is independent of G conditional on Y, and M Q0 retains the coherent correlation of M with H. It is the missing marked covariance insertion. After clipping, add (3.1).

## 4. Genuine native admission, with the carrier kept visible

The native source variables are the Gaussian Omega coordinates only; (Y,G,S) are retained callers/mode labels. The discrete selectors and clock coordinates S are frozen during the native call and are averaged only afterward. This distinction is essential: integrating G as a private root would make the untruncated Q0 correction have exactly zero mean and destroy the current coefficient.

Choose an auxiliary marker buffer u>0. This buffer belongs to the marker output, not to the endpoint X=sqrt(v)G. Normalize N by sqrt(u), and let

    l=L_N/sqrt(u), r=sqrt(l),
    a_curl=2r, mu=r,
    delta_native=E_N/(sqrt(u) r sqrt(d_private)).         (4.1)

If the imported interface requires zero origin, execute and save n0=N(Y,G,S;0), use N0=N−n0, and restore that same n0 after the physical readout. Then |N0|<=2E_N, its private first/curl are unchanged, and its continuous caller first is at most doubled. Use 2E_N in delta_native and the corresponding enlarged caller majorant; all constants below absorb this factor. The saved origin costs a full at-most-nineteen-VALUE replay for each distinct (Y,G,S,source-version) key unless exact-key reuse is proved. Use deterministic majorants with fixed harmless factors, so the complete raw-split source meets every declared first/energy/curl bound. If l<=1, its true first l is at most the padded r, and normalized curl<=2l=r a_curl. The genuine-gradient baseline is identically zero under the valid coisometry. Its Gaussian part can be emitted exactly; there is no gradient-prior error to invent.

Run the imported zero-baseline raw-split self-reserve construction at these ACTUAL parameters and complete private dimension. Require its literal radius, relative-curl, self-reserve, clipping/filter, dimension and precision guards. In particular r must be within its numerical small-radius window and a_curl within its allowed range. The possible energy parameter is a deterministic upper bound, never a measured realized energy. If necessary choose a larger deterministic E_N majorant; all errors change with it.

The five imported physical remainder rows are, before simplification,

    E_N [r^2+r mu+r a_curl+r^3/sqrt(mu)+r^3].            (4.2)

With (4.1), r<=1, they are bounded by

    epsilon_native<=C Lambda E_N L_N/sqrt(u)+floors.    (4.3)

The resulting positive mean LAW R_N targets

    N(mu_N(Y,G,S),u I), mu_N(Y,G,S)=E_Omega[N(Y,G;Omega,S)|Y,G,S].    (4.4)

Every comparator conditions on the same retained (Y,G,S). After S is averaged, this is a mixture of these Gaussian references, not a claimed N(E_S mu_N,uI). Its whole private bank, including all native filters, reserves and the auxiliary keep, is consumed. One must use the source-path theorem, not differentiate (4.3), for any later caller-first request. The precise paths depend on the literal imported native program; no new derivative of a LAW error is asserted.

## 5. Actual positive joint output and current contract

Draw X=sqrt(v)G and sample the clocks/selectors S. Draw a baseline mark M0 from the original finite mark graph at S, independent of the WHOLE native correction bank conditional on (Y,G,S); it may reuse an already available unchanged-key M when that bank remains unread by the callback. Execute R_N conditional on the retained (Y,G,S), and return

    (X,J)=(sqrt(v)G, M0+R_N).                            (5.1)

This is an ordinary positive pushforward law. For any callback a(Y,G) of norm at most one that reads neither S nor a private mark/native bank, the independent auxiliary Gaussian in reference (4.4) has mean zero. Thus (3.1)-(3.2) and the conditional mean consequence of the native W2 comparison imply

    |E[J dot a(Y,G)]-E[M dot a(Y,Z/sqrt(v))]|
       <=2sqrt(2)e^2||M||_2
          +clipping_error+epsilon_native.              (5.2)

This includes every admissible bounded gradient callback defining a current. It is a retained-carrier weak-current certificate on the exact executing joint law (5.1), not tensor matching.

The actual output J has its auxiliary Gaussian marker variance; its full L2 energy is not E_N or the original small mark energy. It must not be inserted nonlinearly into a force, declared a RAW history, or used by a callback that reads native private roots. No Gaussian carrier is subtracted. The only cancellation is the legitimate zero conditional mean of an independent Gaussian in a LINEAR mark observable after its bank is consumed.

## 6. Concrete grades and genuine limits

For the existing pair graph and a finite local five-VALUE mark with energy A^3, at fixed supplied physical block size and frozen public logarithms, equations (1.3),(2.2)-(2.4) give

    E_N<=Lambda_b A^7/v,
    L_N<=Lambda_b A^(37/9)/v.                           (6.1)

All n,B_G,R_M and cap-log factors are in the explicitly finite polynomial Lambda_b; this is not a dimension-uniform unrestricted-block assertion. The three main first paths have A powers 5,37/9,7, divided by v, so 37/9 is the controlling one for A<=1.

The substantive error rows in (5.2) are

    current truncation: C_b A^11/v^2 sqrt(n),
    native correction: Lambda_b A^(100/9)/(v^2 sqrt(u)) sqrt(n),
    clipping: the explicit exponential-tail term (3.1). (6.2)

When u is comparable to v and v>=c A^2 K^2, the padded normalized radius satisfies

    r<=Lambda_b A^(37/18) v^(-3/4)
      <=Lambda_b A^(5/9) K^(-3/2).                     (6.3)

Thus it has a strict positive A margin at the existing endpoint cutoff. At v comparable to A^2, the native row has exponent 55/9, and the truncation row exponent 7. All literal native numerical guards still must pass; fixed b and fixed public-log exponents admit a nonempty sufficiently-small-A window. This statement does not transfer to arbitrary b=D without retaining the b-dependent guard and error factors.

The rows (6.2) are single-block current rows. For supplied block-separable sources, block-local current vector fields may be combined in Euclidean L2, yielding the root-sum-square of their errors. No claim covers an arbitrary globally coupled callback by multiplying block likelihoods; that would increase the rank and lose the dimension-free energy constant.

No original-gradient count depends on the selected likelihood-current rank p: the p=1 result is explicit above, while higher p changes scalar polynomial work and source first constants. This does not eliminate all depth/order costs in the unresolved true-history induction.

## 7. Fully expanded costs and floors

One N evaluation uses the full retuned pair graph Q_H=14 plus every added coherent mark VALUE. For the concrete local five-call mark, Q_N<=19 before exact-key reuse. Its complete private bank excludes the retained G, hence the pair part has at most 16n roots; add every mark root not already among those rows, and any coisometry padding. Vector/scalar work is O(n) for p=1. The baseline M0 costs five additional VALUES if not already available. Every native replay retains the sampled S; it does not resample or differentiate its discrete selectors. S is integrated only after the conditional linear-current comparison. Scalar clock randomness and any separately requested inactive-slot Gaussian buffer retain their own costs; they are not included in the coordinate-root abbreviation.

Let N_raw(r,a_curl,mu,L,epsilon) and d_raw be the fully expanded original-source occurrence and Gaussian-dimension recurrences of the ACTUAL imported zero-baseline self-reserve program. Then a safe uncached bill is

    Q_total <=5+19 N_raw+Q_capture+Q_restore+Q_numeric,
    roots <=n + d_M0 + d_raw,                            (7.1)

where the first n is the retained carrier G, d_raw counts every full fresh private source bank plus native roots, and d_M0 is the complete baseline mark bank. If the interface anchors source zero, each origin costs another full nineteen-call replay unless exact-key reuse is established. Changing Y,G,v,cap scales, frozen clock labels or source version replays all affected sites and the corresponding saved origin. No source or conditional coefficient is a unit-cost leaf.

At fixed native orders, the imported recurrence is a fixed public-log polynomial multiplying these nineteen-call sources, with the actual padded radius (4.1), complete dimension, buffers and floors. The program is not numerically instantiated against an unspecified g oracle, so N_raw is a declared imported recurrence rather than an invented numeric count.

Freeze the retuned delta, geometry, caps, u,v and complete source graph before assigning numerical floors. Use actual downstream multipliers including 1/v, 1/v^2, source secants, scalar Hermite contractions, cap derivatives, native normalization, self-reserve/filter weights and output readouts. The finite precision ledger must sum those weighted absolute floors to its allocated tolerance. The target truncation, retuned secant bias, clipping and native substantive rows above are never absorbed into numerical floors.

## 8. Exact contribution to the open task

This constructs the previously absent finite retained-carrier coefficient for the one-mark/pair-covariance current, including the E[M H] correlation that a marginal pair LAW forgets. It is admitted through an actual clipped original-VALUE source and a retuned source-first calculation. Its native comparison retains the correct visible Gaussian carrier.

It does not provide an original-VALUE graph for the true Delta3 mark, the base conditional covariance, or exact singleton history fields. It does not certify higher marked proper cuts or permit nonlinear reuse of an auxiliary noisy marker. Consequently it is a useful concrete local port for a future coherent history join, not completion of the general-D m4 or arbitrary-order induction.

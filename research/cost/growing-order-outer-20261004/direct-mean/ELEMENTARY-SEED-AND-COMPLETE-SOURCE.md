# Elementary mode-Gaussian seed and the complete same-heat mean source

2026-10-04. Companion to ANY-ORDER-DIRECT-MEAN-OUTER.md. This source construction supplies only the VALUE, actual first/caller, complete pair and numerical ports needed there. No retained graph, genuine-gradient reference, proxy, covariance compiler or marked-domain theorem is used.

## Theorem

For every fixed maximum integer J>=1, fixed hidden depth d, and finite moment list, the selected original-gradient Picard history admits finite sources X_(d,j), 1<=j<=J, with

    state law order j+1/2,
    complete adjacent canonical-force Lp gap O(sqrt(D) A^j),
    complete original-query cost and Gaussian dimension/Dim <= Lambda A^-2(j-1),
    actual caller first O(1), actual canonical-force caller first O(sqrt(A)),
    complete canonical-force private first O(A).

All law and strong estimates include separately restored original numerical/caller profiles. Conditional source moments are integrated only at their actual original caller; internal ancestors remain private. The physical history coefficient on original gradients is O(A), with public-log factors, and every source call uses A or A/2, recursively only fixed-depth constant ratios. The selected outer histories have this property since their original-gradient predecessor is h^2A_il and h^2<=A/L^C.

The source has no algebraic condition of the form D<=A^-b. Dimension appears in sqrt(D) moments and in public logarithms/precision only. The effective small-A guard may depend on J,d and public-log row masses. The theorem is about query count, not dimension-free arithmetic, storage or Gaussian generation.

## 1. Finite seed, law, firsts and moment profiles

Fix 0<A<=1/8 and a physical caller y. The exact mode b solves b=y-A grad V(b). Execute the fixed count

    b_0=y, b_(m+1)=y-A grad V(b_m),
    S(A,y)=b_M+sqrt(A) Z.

Here M is selected before caller differentiation, from the chosen numerical tolerance and physical-query profile. Since grad V is 1-Lipschitz,

    |b_M-b| <= A^(M+1)|grad V(y)|/(1-A).

The exact b is an analytical reference, never a queried proximal oracle. The derivative recurrence D b_(m+1)=I-A H_m D b_m gives

    ||D b_m||<=1/(1-A),
    ||D b_m-I||<=A/(1-A).

Each requested first/adjoint is an original HVP at the already recorded gradient point. The private first of S is sqrt(A)I. Its centered Gaussian moment is C_p sqrt(D A). The absolute seed/mode error retains the displayed caller profile instead of being bounded uniformly over unbounded y.

For the law bound, standardize X=b+sqrt(A)U. The posterior potential is

    W(U)=|U|^2/2+V(b+sqrt(A)U)-V(b)-sqrt(A) grad V(b)·U.

Its Hessian lies between I and (1+A)I. Couple the stationary diffusions dU=-grad W(U)dt+sqrt(2)dB and dG=-Gdt+sqrt(2)dB synchronously. The added drift r(U)=sqrt(A)[grad V(b+sqrt(A)U)-grad V(b)] has |r(U)|<=A|U|. Strong convexity and grad W(0)=0 give E|U|^2<=D. The stationary contraction inequality yields

    ||U-G||_2 <= A sqrt(D),
    W2(Q_(A,y),N(b,AI)) <= sqrt(D) A^(3/2).

For the finite seed add |b_M-b|. This is the familiar mode-Gaussian estimate used in LOW30's terminal serial approximation, but the above proof explicitly supplies it here with no D-domain guard.

## 2. Literal known-center quarter map

No short refresh is needed for this mean-only source. Let T=pi/2, t_i=iT/N, and

    w_iq=cos(t_i-t_(q+1))-cos(t_i-t_q), q<i.

The weights are nonnegative and sum to 1-cos(t_i). Given a previous state x0, draw fresh Z after its caller is exposed. Set

    x_i^[0]=y+cos(t_i)(x0-y)+sqrt(A)sin(t_i)Z,
    x_i^[m+1]=x_i^[0]-A sum_(q<i) w_iq grad V(x_q^[m]).

Return x_N^[M]. All gradient locations and deterministic rows are literal. Cos(T)=0 is used exactly (or any encoding residual gets its actual numerical budget). The Hessian bound and positive row sums imply the finite endpoint's x0-Lipschitz constant <=A/(1-A); this follows from the actual recurrence, not posterior invariance. Internal path firsts are bounded by (1-A)^-1.

The exact Hamiltonian system with potential |x-y|^2/(2A)+V(x), time scaled by sqrt(A), preserves Q_(A,y) after refreshing its velocity sqrt(A)Z. Its same-Z position endpoint at T contracts an initial position difference by A/(1-A), by the same Volterra inequality.

For discretization, use the ANALYTICAL exact mode b to bound the exact Volterra solution and the discrete FIXED POINT. The identities b-y=-A grad V(b) and sum w_iq=1-cos(t_i) cancel constant forces in those two equations. They do NOT cancel in the finite layer-zero predictor: its centered path contains A(1-cos(t_i))grad V(b). The exact and discrete-fixed-point centered paths have moment bound C_p(||x0-b||_p+sqrt(D A)), while the finite-M iterate has the additional residual C A^(M+1)|grad V(b)|. Compare the exact solution with the discrete fixed point, then terminate the Picard contraction. The actual bound is

    Lp error <= C_p sqrt(D)[A^(M+3/2)+A^(3/2)/N]
                 +C A^(M+1)|grad V(b)|
                 +propagated finite-seed/mode numerical residual.

The last two terms are caller-dependent NUMERICAL allowances, not intrinsic sqrt(D) errors. Since |grad V(b)|<=|grad V(y)|/(1-A), they have an explicit original caller profile. When the input centered path has its own numerical displacement profile, propagate it through the actual A/(1-A) endpoint contraction. This is precisely the quarter recurrence proved in the pinned protected-short source, Section2; the proof uses only these path moments, Lipschitz original gradient and positive kernel, so its earlier CW7/retained assumptions are irrelevant here. Finite mode residuals add their explicitly retained numerical profile. The same bound can be verified by the Volterra Picard contraction A and the O(sqrt(D A)) path time-Lipschitz estimate; the rectangle quadrature error is A times its O(1/N) modulus.

Freeze level-j N_j=ceil(C_J A^-(j-1)), and M_j=max(j,M_num), with M_num selected so A^(M_num+1)/(1-A) meets the fixed original numerical/caller budget after every downstream amplification. M_num is O_J of public logarithms and does not add an inverse-heat cost exponent. Define K_j by j such maps from one complete S seed, all maps at the j-grid. Define K_0=K_1 literally. Exact invariance, contraction and the seed law bound give

    W2(K_j,Q_(A,y)) <= C_J sqrt(D) A^(j+1/2)+e_num.

A complete K_j uses O_J(N_j) original-gradient values plus its finite mode; its heat exponent is j-1. It is not run at a finest common grid for free.

The actual finite caller first is I+O_J(A): differentiate all rows with respect to y, subtract I, and use D b_M-I=O(A), cos(T)=0 and bounded Hessians. The complete private first is sqrt(A)P_j+O_J(sqrt(A)A); only the last quarter's Z enters its affine endpoint, and every older private path gains the terminal factor A. The whole private first remains O_J(sqrt(A)). These are actual finite derivatives.

## 3. Complete known adjacent pairs

For j>=2, draw one complete seed and j complete quarter momentum blocks. Fine uses all j on its j-grid; coarse uses the last j-1 on its own (j-1)-grid. Each is its exact named marginal. At a shared packet, finite same-grid input contraction is C A; changing grids costs C sqrt(D) A^(j-1/2), by comparing both with the same exact quarter map. The fine first-step versus seed displacement is O_p(sqrt(D A)) plus the same numerical profile. Thus the literal recurrence over the common suffix gives

    ||K_j-K_(j-1)||_p <= C_J sqrt(D) A^(j-1/2),
    ||sqrt(A)grad V(K_j)-sqrt(A)grad V(K_(j-1))||_p
           <= C_J sqrt(D) A^j.

For j=1 the pair K_1/K_0 is identical at a common caller. At two callers add the actual K_j center first, then draw the fresh pair after both entering callers are known. No coupling is inferred from a W2 estimate.

Deterministic all-private zeros are treated separately. The seed equals b_M; its distance to the analytical mode is the original finite residual above. The unanchored layer-zero predictor does not fix b, so each terminated quarter also has the deterministic A^(M+1)|grad V(b)| residual just identified. Its actual endpoint contraction transports previous zero errors, and these new finite-Picard/mode residuals are included separately. One may assign them the same target numerical allowances at their actual caller. In particular there is no use of a Gaussian Lp estimate evaluated at zero. Strong mean identity below requires only matched zeros, but any requested zero-gap bound can also be paid this way.

## 4. Empirical hidden recursion and target chronology

The selected history at depth d has physical form

    H=h_base+sqrt(A) F_history,
    F_history=sum_s B_s F_(d-1,j,s),
    F_(d-1,j,s)=sqrt(A)grad V(X_(d-1,j,s)),
    sum_s||B_s||<=Lambda.

For the outer Picard history B_s=-(h^2/A)A_is. The complete right-hand blocks are independent full replays conditional on the original phase caller, as literally prescribed in LOW30. Every alias within each block is retained. The same source hierarchy supplies coherent adjacent history pairs by executing the constituent full pairs in their original chronology. More general correlated histories require the actual topology/first/row contract; no such generalization is needed here.

For each history level j form M_j from the complete banks of Section2 of the outer theorem:

    N0=ceil(C A^-2(j-1)), N_l=ceil(C A^-2(j-l)),
    M_j=average F_history,0+sum_l average Delta_history,l.

The conditional mean is exactly E F_history,j, and the strong error is C sqrt(D) A^j. For level0 use one complete history sample. Then execute the sufficient statistic and observation decoder

    T=h_base+sqrt(A)M_j+sqrt(A/(k+1)) Z_T,
    R_i=T+sqrt(A)(E_i-Ebar), i=0,...,k,
    X_0=finite mode at R_0,
    X_i=K_j(A/2,(X_(i-1)+R_i)/2), i=1,...,k.

All K_j children are fresh after their entering caller is exposed. One fixed k=k_J=O_J(log(1/A)+L) is used for every compared order level. With exact T replaced by E T plus its displayed independent Gaussian, the R_i are iid N(q,A I), q=E T. The exact posterior Gibbs transition contracts by 1/2, and the decoder law proof gives

    state error <= C[ sqrt(A)||M_j-EF_history,j||_2
                     +known-child law error+2^-k sqrt(D A)].

The finite history center differs from its exact-score ideal center by at most C A delta_(d-1,j), because its original-gradient row is O(A) and grad V is 1-Lipschitz. Therefore

    delta_(d,j)<=C sqrt(D)[A^(j+1/2)+A^(J+1/2)]
                        +C A delta_(d-1,j)

with the actual separate numerical profiles. This proves state order j+1/2 at every fixed d.

The actual adjacent mean banks have difference O_p(sqrt(D)A^(j-1)) for j>=2. Share the untouched statistic Gaussian, the observation array, and known-child suffix couplings after both callers are exposed. The decoder's center contraction is (1+C A)/2<3/4, so its state gap is O_p(sqrt(D)A^(j-1/2)) and terminal canonical-force gap O_p(sqrt(D)A^j). At j=1 the coarse statistical errors are O(A), adequate for the same weaker declared gap. This is the complete paired-marginal construction, not separate marginal matching.

The bank's private first is <=[sum_l Lip(F_l)^2/N_l]^1/2=O_J(A); its coherent caller first is O_J(sqrt(A)). The statistic multiplies by sqrt(A); the actual known children have I+O(A) caller first; decoder differentiation once gives the same complete O(sqrt(A)) state-private first and bounded state caller first. The affine observation/noise carrier is coisometric by the covariance identity, independently of any proxy projection. These actual firsts supply the next depth's hypothesis.

## 5. Cost fixed point, numerical profiles, and absence of dimension gates

Write C_(d,j) for a complete source or pair and its requested original first sweep. Then

    C_(d,j)<=Lambda_J[k Q(K_j)
       +N0(j) C_(d-1,0)
       +sum_(l=1)^j N_l(j)(C_(d-1,l)+C_(d-1,l-1))
       +actual finite mode/numerical work].

Induction yields C_(d,j)<=Lambda_(J,d) A^-2(j-1). The actual Gaussian dimension divided by Dim obeys the same recurrence. No ambient-dimensional tail estimate is used: all moment estimates use complete physical output norms, Hilbert-space averaging, finite positive-kernel propagation and exact Gaussian decoder covariance. Thus the large tape introduces no D<=A^-b premise. The coarsest sample level has exponent zero at every fixed d.

The mode recurrence gives an explicit initial numerical envelope A^(M+1)|grad V(y)|/(1-A). Every downstream path has a known fixed-power amplification after J,d are fixed. Source centers are finite affine combinations of original caller coordinates, Gaussian physical vectors and Lipschitz original gradient values. Their DISPLACEMENT FROM THE ORIGINAL CALLER and their ORIGINAL-GRADIENT profiles, rather than their absolute coordinate norms, are bounded by an order/public-log factor times |grad V(y)|+sqrt(D A) and any captured phase velocity. Indeed |grad V(z)|<=|grad V(y)|+|z-y| and the physical history row contributes O(A) times earlier gradient profiles; observations contribute physical Gaussian moments and the decoder has a bounded geometric resolvent. This gives a fixed-depth profile recurrence B_d<=Lambda_J(B_(d-1)+1), with B_0 finite, independent of the number of empirical samples because complete averages have their literal absolute weights. No maximum over the whole ambient tape is used. Enumerate the complete graph and reduce each primitive/mode tolerance to pay its known path amplification. Because A^-b has logarithmic bit length O_(J,d)(log(1/A)), fixed finite iterations require only public-log overhead. Source banks sum common numerical bias with their literal absolute weights; it is never divided by sqrt(N).

At the outer caller theta=(x,v^-), these profiles become Lambda[|grad V(x)|+|v^-|+sqrt(D)], and the fresh G0 moment is integrated afterward. This is exactly the numerical profile used in the chronological outer moment recurrence. The required mode counts remain caller-independent; no adaptive unrecorded stopping test based on grad V at a random center is used.

The only inequalities requiring a guard are fixed-J,d resolvents, A times public-log row masses, interpolation/precision requirements, and comparison contraction. None contains a positive power of D except inside public logarithms or the explicitly factored sqrt(D) allowance. This proves the domain-free-in-physical-dimension source specialization needed for the direct outer rate, subject to independent verification of the finite quarter bound and the stated original Gaussian decoder identity.

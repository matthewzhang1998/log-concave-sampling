# The genuine next history: m4 and its missing marked-current port

2026-10-05. A source-qualified target-isolation result for the anchored general C2 convex-gradient class. This is not a new accuracy-raising producer. It identifies the exact next algebraic object, proves its difference bounds, and specifies why the sealed canonical-m3 family does not supply it.

## 1. Result

Let X be stationary standard OU and g=grad U, g(0)=0, 0<=Dg<=A I, A<=1/2. Retain the original history definitions

    F_0,t=0,
    F_k,t=integral_0^infinity e^-s g(X_(t+s)-F_(k-1,t+s)) ds,
    m_k(y)=E[F_k,0 | X_0=y].

The actual NEXT history for a mean reentry is F3, and the next mean target is m4. Set, on ONE coherent future conditioned on Y,

    V=F2, U=F3, Delta=U-V,
    U_theta=V+theta Delta, 0<=theta<=1.

The exact new conditional targets, including every cross term, are

    delta_m=E[Delta|Y],
    delta_Sigma=Cov(U|Y)-Cov(V|Y),
    delta_K3=kappa3(U|Y)-kappa3(V|Y),
    delta_K4=kappa4(U|Y)-kappa4(V|Y).

Their integrated one-energy bounds are respectively

    ||delta_m||_2 <= A^3 sqrt(D),
    ||delta_Sigma||_(2;HS) <= C A^4 sqrt(D),
    ||delta_K3||_(2;HS) <= C A^5 sqrt(D),
    ||delta_K4||_(2;HS) <= C A^6 sqrt(D).                 (1)

The target mean change itself satisfies ||m4-m3||_2<=A^4 sqrt(D). It is nonzero at exactly that order even for g(x)=A x:

    m4(x)-m3(x)=-A^4 x/16,
    Cov(F3|x)-Cov(F2|x)
        =[3A^4/8-9A^5/16+19A^6/64] I.                  (2)

Thus an unshifted higher-cumulant correction alone cannot perform the next history step. All cumulants of ranks three and higher vanish in (2), while the next mean and covariance obligations remain.

The new algebraic port is a finite positive original-VALUE implementation of the HISTORY-MARKED current of the coherent pair (F2,F3-F2), with the full shifted nested ancestry. A positive joint block-Gram is the precise covariance target; separate independent F2 and Delta banks lose its off-diagonal block. The original-gradient path/star packets have no such supplied shifted-history vertex. Sections 6–8 give the exact port, native/readset obligations and honest cost placeholders. No small first for Delta, RAW promotion of the sealed m3 graph, signed probability law, or execution of an analytical field is used.

## 2. True histories, differences and their regularity

Write Delta_k=F_k-F_(k-1). The A-Lipschitz contraction on the SAME stationary path gives

    ||F_k||_2 <= A/(1-A) sqrt(D),
    ||Delta_k||_2 <= A^k sqrt(D),
    Lip_complete(F_k) <= A/(1-A)=:L.                    (3)

Conditional centering is an L2 contraction, so the centered error

    R=Delta-E[Delta|Y]

has integrated energy e=||R||_2<=A^3 sqrt(D). The available full private first of Delta is at most 2L. It is not A^3 merely because its energy is A^3 sqrt(D).

The true histories and identities below are analytical references. One may prove them first on common finite OU-time Gaussian-bank approximants, with every nested history on the same finite future, and pass to the strong limit. This is a justification of analytical laws only. It is not a proposed executing path-grid or cold nested-history algorithm, and no query-complexity claim is based on such a limit.

The common-bank requirement is substantive. Fresh banks are required for independent COMPLETE services in the executing join, but F2 and Delta INSIDE a history-marked reference or coefficient supplier must have their coherent shared ancestry. Independence between services never licenses independence between slots of one covariance/cumulant target.

## 3. Exact mixed targets and one-energy grades

Suppress Y, and center V, U and Delta. Denote the centered variables by v,u,r. For a mixed tensor, physical slot permutations are retained; Sym below means the AVERAGE over the indicated positions, so the displayed binomial multipliers are literal.

### 3.1 Covariance

    delta_Sigma=E[v r^T+r v^T+r r^T]
               =integral_0^1 K2(theta) dtheta,
    K2(theta)=Cov(Delta,U_theta)+Cov(U_theta,Delta).     (4)

The covariance identity with the derivative/Riesz operation on the Lipschitz slot gives

    ||Cov(R,W)||HS <= Lip(W) ||R||_2.

Using the telescope Cov(U)-Cov(V)=Cov(R,U)+Cov(V,R), not a falsely graded quadratic expansion, proves

    ||delta_Sigma||HS <=2 L e_R(Y).                    (5)

The r r^T term has only

    Tr Cov(Delta) <= e_R(Y)^2,
    ||Cov(Delta)||HS <=2 L e_R(Y).                      (6)

After integration these are A^6 D and C A^4 sqrt(D), respectively. They are different norms. The first does NOT imply a dimension-uniform A^6 sqrt(D) Hilbert-Schmidt bound.

### 3.2 Third cumulant

Multilinearity gives the exact seven-term expansion

    delta_K3=3 Sym kappa(Delta,V,V)
              +3 Sym kappa(Delta,Delta,V)
              +kappa3(Delta)
             =integral_0^1 K3(theta) dtheta,
    K3(theta)=3 Sym kappa(Delta,U_theta,U_theta).       (7)

For centered Lipschitz P,Q with firsts at most L, Gaussian Poincare implies

    Var(P^T M Q)<=4 L^4 ||M||HS^2.

Consequently ||E[R tensor P tensor Q]||HS<=2 L^2 ||R||_2. The three-slot telescope with remaining slots U or V gives

    ||delta_K3||HS<=6 L^2 e_R(Y).                      (8)

If (7) is expanded instead, any slot occupied by a second or third Delta has only the known first 2L. All seven nonempty marked terms therefore have the common certified grade A^5 sqrt(D), up to constants. The mark count does not give extra powers for free.

### 3.3 Fourth cumulant

All mixed cumulants include all three pairwise covariance subtractions. Exact multilinearity gives the fifteen nonempty marked terms

    delta_K4=4 Sym kappa(Delta,V,V,V)
              +6 Sym kappa(Delta,Delta,V,V)
              +4 Sym kappa(Delta,Delta,Delta,V)
              +kappa4(Delta)
             =integral_0^1 K4(theta) dtheta,
    K4(theta)=4 Sym kappa(Delta,U_theta,U_theta,U_theta). (9)

The sealed dimension-safe fourth comparison proves, with one arbitrary centered error slot,

    ||kappa(R,P,Q,S)||HS<=6 L^3 ||R||_2.

Its covariance-subtracted mixed Wick cubic is essential; an unsubtracted third product would introduce dimension-sized trace pieces. Telescoping four slots with remaining U/V firsts at most L yields

    ||delta_K4||HS<=24 L^3 e_R(Y).                     (10)

Again every nonempty marked term in the expanded form has only grade A^6 sqrt(D) under the available ports. In particular multiple Delta slots cannot be deleted as allegedly higher order.

Equations (5), (8), (10) hold conditionally on the actual Y. Integrating their squared bounds gives (1), and applies only with the same actual retained-Y law. It is not a uniform substitution theorem for an arbitrary non-Gaussian caller.

### 3.4 Why the missing first/cut grade is a genuine logical gap

The following is a counterexample to an inference from PORTS, not an example of actual force histories. Take A=1/n, D=n^4, and Delta=A G_1 e_1. Then

    Lip(Delta)=A,
    ||Delta||_2=A=A^3 sqrt(D),
    ||Cov(Delta)||HS=A^2,
    A^6 sqrt(D)=A^4.

Thus energy and first alone cannot certify the proposed stronger grade; the ratio is n^2. An additional actual-history directional-energy or marked proper-cut theorem could improve it, but is a genuinely new theorem. This packet makes no impossibility claim for such a theorem.

## 4. The full new buffered reference and paid target differences

At an old outer LAW node keep its actual retained tuple (z,Y,N,G_h), conditional variance v, and q. The correct updated analytical reference is

    -q m3(Y)+S3
      +[-q^3 kappa3(F3|Y)/6]:H2^(C3)(S3)
      +[+q^4 kappa4(F3|Y)/24]:H3^(C3)(S3)
      +sqrt(v/2) G_keep,
    C3=(v/2)I+q^2 Cov(F3|Y),  S3~N(0,C3).              (11)

H2 and H3 have exactly the same covariance-normalized Wick conventions as the sealed quartic join. This is one positive pushforward law. No Gaussian square root, conditional moment, Hermite coefficient or exact history is an executing leaf.

The changes from the F2 reference are not only its mean. They are delta_Sigma, delta_K3, delta_K4, and the induced CHANGE OF WICK COVARIANCE from C2 to C3. Coefficient maps at the old covariance cannot simply be carried over after adding a purported covariance-increment sample. Their carrier regression, trace corrections and same-law cross feedback must be paid on the same visible new carrier. A fresh standalone tensor packet does not establish that comparison.

For orientation, if each exact target restoration is paid using the sealed stable buffered comparison at v, the NEW history-only rows after the outer A-Lipschitz consumer have scales

    mean supplier:       A e3,
    covariance delta:    C A^5/sqrt(v) sqrt(D),
    third delta:         C A^6/v sqrt(D),
    fourth delta:        C A^7/v^(3/2) sqrt(D).          (12)

The coefficient rows in (12) presuppose the same-carrier comparison and proper cuts required by that sealed comparison; (12) is not a new executing update program. A safe implementation must either rebuild the full new references and prove their joint comparison or prove the incremental same-carrier comparison explicitly.

For the old actual positive outer rule, the corresponding scalar upper ledgers are

    A^5(w^-1/2+1),
    A^6(w^-1+log(1/eta)),
    A^7(w^-3/2+eta^-1/2),                              (13)

times sqrt(D) and the declared factors. These are history differences ONLY. The larger existing mismatches between F2 and its finite covariance target, between F2 and H in the cubic target, and between F2 and H_tau in the quartic target remain separate. Neither the prefix, endpoint, native mean, corrected-Gram, nor final completion debt disappears.

At w=A, the bulk leading powers in (13) are 9/2,5,11/2. This is a diagnostic threshold, not permission to use eta=A^2 or violate any literal native guard. In particular m3 mean error O(A^(4-epsilon)) becomes O(A^(5-epsilon)) through its actual mean service, but unchanged forcing prevents this from being an accuracy-raising theorem.

## 5. Exact original-g nested response for the mean increment

There is an exact four-force RESPONSE CHAIN behind the next mean, distinct from an unshifted cumulant tree. For k=2,3,4 define the analytical matrix

    J_(k,t)=integral_0^1
       Dg(X_t-F_(k-2,t)-theta Delta_(k-1,t)) dtheta.

Then

    Delta_(k,t)=-integral_0^infinity e^-s
                  J_(k,t+s) Delta_(k-1,t+s) ds.        (14)

In particular, with T_j=s1+...+sj,

    Delta_(4,0)=-integral_(R_+^4) e^-(s1+s2+s3+s4)
       J_(4,T1) J_(3,T2) J_(2,T3) g(X_T4) ds1...ds4.   (15)

Each J contains its own theta integral. Every F1,F2 and Delta in these shifted sites is on the SAME future used by its descendants. Matrix order in (15) is the displayed order; in multiple dimensions it may not be permuted. Norms of the J's are at most A, so (15) recovers ||Delta4||_2<=A^4 sqrt(D). Conditional expectation given X0 yields m4-m3.

The derivative matrices in (14)–(15) are analytical. They are NOT new producer leaves. A VALUE realization must implement their response using actual finite original-g queries, preserve all shifted ancestor values and their replays, and prove its own admission. The already sealed C0/C1/C2 four-force paths are built for original-gradient coefficient targets of H; (15) contains a history-shifted nonlinear composite source. The original g gradient qualification alone is not a qualification of that composite source.

General C2 does not permit replacing all shifted J's by Dg(X_t) with a uniform new A-power: no uniform Hessian-Lipschitz or higher-derivative modulus is supplied. Gaussian smoothing is allowed, but the exact coherent smoothing bias, derivative/cut losses, VALUE implementation, and all affected history ancestors must enter a new proof. Choosing a smoothing scale is not a certificate for these obligations.

## 6. Exact positive block-Gram target: the minimal covariance algebra

Here is a precise new target rather than the phrase 'improve covariance.' On a common finite Gaussian-bank reference conditioned on Y, center the pair

    W=(V-E[V|Y], Delta-E[Delta|Y]) in R^(2D).

Let P_s denote OU on its COMPLETE private bank, not on the physical Y. Its covariance has the exact positive Gram representation

    B(Y)=Cov(W|Y)
        =2 integral_0^infinity E[(D P_s W)(D P_s W)^T|Y] ds
        =[[Sigma2, C],[C^T, SigmaDelta]],              (16)
    C=Cov(F2,Delta|Y).

Every integrand in (16) is positive semidefinite. The actual new covariance is

    Sigma3=[I I] B [I I]^T
           =Sigma2+C+C^T+SigmaDelta.                  (17)

One can equivalently use the path [I theta I] B [I theta I]^T. Its derivative is exactly K2(theta) in (4), including the 2theta SigmaDelta term. The interpolation stays positive because it is a Gram contraction, irrespective of whether the covariance difference itself has a sign.

The standard energy identity proving (16) is

    -d/ds E[(P_s W)(P_s W)^T]
       =2 E[(D P_s W)(D P_s W)^T].

The terminal limit is zero because W is conditionally centered. Limits from finite references justify the identity for actual histories when needed; (16) remains an analytical target.

A positive increment-only covariance service independent of the old service can add only a PSD matrix; the available target bounds do not certify that delta_Sigma is PSD. If its positivity were separately proved, a service explicitly targeting the FULL delta_Sigma could in principle encode the cross terms. By contrast, an independent Delta-only bank necessarily adds only SigmaDelta, omits C+C^T, and is wrong even when delta_Sigma happens to be PSD. In the linear test, the missing cross block is

    C+C^T=3A^4/8-9A^5/16,
    SigmaDelta=19A^6/64.

This is a failure already at grade four. No sign claim for the true nonlinear-history delta_Sigma is needed for this conclusion. Generic admissible port data do not force a PSD difference, but this packet does not claim a nonlinear-history counterexample establishing indefiniteness.

A sufficient covariance producer would return the COMPLETE positive block law N(0,(u/2)I_(2D)+B(Y)), with a finite original-VALUE realization of (16), then apply [I I]. Its readout has norm sqrt(2), the output Gaussian buffer is exactly u I_D, and its target covariance is u I_D+Sigma3. The normalized pair-source first and doubled physical dimension must be charged; a normalized error epsilon_block becomes at most sqrt(2) epsilon_block after this readout. Alternatively it could directly realize the full (17) through an equally complete positive law. Neither method needs a signed Gaussian law. Both need a new supplier for the marked shifted-history coefficient; an exact block Gram formula is not that supplier.

### 6.1 A positive current formulation that keeps the mark

An alternative exact formulation exposes the joint law that an incremental implementation must preserve. Conditionally on Y let G be a fresh D-Gaussian independent of the WHOLE coherent history and put

    Z_theta=sqrt(v) G-q U_theta,
    nu_theta=Law(Z_theta|Y),
    b_theta(x)=-q E[Delta | Z_theta=x,Y].               (17a)

For every compactly supported smooth scalar test phi, differentiation under the finite-bank expectation and then its strong-limit extension give

    d/dtheta E[phi(Z_theta)|Y]
        =-q E[<Dphi(Z_theta),Delta>|Y]
        =E[<Dphi(Z_theta),b_theta(Z_theta)>|Y].         (17b)

Thus nu_theta is a positive conditional law path satisfying its continuity equation, and its current is a joint (Delta,Z_theta) observable. Its mean, covariance and rank-three/four cumulant derivatives are exactly -q delta_m, q^2 K2(theta), -q^3 K3(theta), q^4 K4(theta). A marginal approximation to nu_theta alone does not specify this marked current.

Conditional Jensen gives E[|b_theta(Z_theta)|^2|Y]<=q^2 E[|Delta|^2|Y]. Equivalently, the original SAME-G and SAME-history coupling directly gives

    W2(nu_1,nu_0|Y)<=q ||Delta||_(2|Y).                (17c)

After an A-Lipschitz terminal this is the honest grade-four change, not a smaller error. To exploit the mark in a higher-order remainder, a finite positive implementation needs a relative weak/current certificate that preserves this exact joint ancestry; two unrelated absolute LAW guarantees do not yield an error proportional to the small Delta energy. The native current/response supplier for (17a)-(17b) is another precise formulation of the same missing port. Conditional expectations in (17a) are analysis only and never a score or current oracle.

## 7. Concrete missing port and contracts that would constitute progress

The following bounded port is the exact next target for construction. Its input is the actual original g, the retained Y and unread labels, the coherent-history specification in (14)–(16), an allocated buffer u>0 and an absolute floor budget. It may not take F3, Delta, conditional tensors, Dg, a path integral, or a previous m3 LAW as an executing oracle.

It must produce:

1. A complete positive buffered covariance service whose target is Sigma3, preferably by the block-Gram (16)–(17), with an explicit finite original-VALUE program. It must improve the existing unresolved A^4 sqrt(D) covariance target error at the actual buffers if it is to advance beyond its threshold in (13).
2. Compatible complete positive skew and quartic services for kappa3(F3|Y), kappa4(F3|Y), or explicit history-marked updates K3(theta), K4(theta), retaining every term of (7),(9). New targets must be restored under the actual Y law, with a simultaneous carrier/Wick regression proof.
3. A qualified original-VALUE response vertex for the shifted nested chain (15), or an equally explicit finite source specification replacing it. Its full private and caller firsts, energy, curl/split and each marked proper-cut norm must be proved from that graph. A mere smaller mean error or smaller Hilbert energy does not satisfy a first/cut port.
4. Native guards using the actual normalized radii of those vertices, positive allocations, every shield/readout factor, a bound on the full dimension, and the exact finite-clock approximation error. If a difference slot is assigned a small radius, its small FIRST or relevant proper-cut norm must actually be proved.
5. Whole-law error and same-law cross feedback. The old mean, old corrected-Gram, new history update, cubic, quartic and keep cannot be compared separately and then reassembled by cumulant matching. All consumed banks stay consumed.

Within this five-service mean/covariance/skew/quartic architecture, updating the suppressed F3 ancestry requires these target obligations or a proved replacement for them. This proposed port is not asserted necessary for every possible next-mean architecture, and it is not sufficient by itself for a grade-above-four mean theorem. Such a theorem also needs an improved prefix/endpoint and a service/completion error ledger whose every active row advances. For a desired terminal mean grade p, simple per-node sufficient target tolerances are

    mean error <= A^(p-1) sqrt(D),
    covariance HS error <= A^(p-1) sqrt(v) sqrt(D),
    cubic HS error <= A^(p-1) v sqrt(D),
    quartic HS error <= A^(p-1) v^(3/2) sqrt(D),         (18)

subject to the already stated joint comparison. Integrated allocations can be weaker than these pointwise sufficient tolerances and must use the actual finite positive weights.

No existing theorem is silently relabeled as satisfying items 1–5. The exact algebraic port remains open, and that is the principal result of this target-isolation packet.

## 8. Readsets, positivity, native guards and costs

The sealed m3 RAW graph S3 can be reused ONLY as the mean supplier. Rerun its native compiler on -q S3 with its COMPLETE fresh tape and original split at the actual allocation. Do not reuse or subtract the known carrier of its completed terminal LAW. Its possible order-A^2 RAW covariance defect is irrelevant to the mean service but fatal to a RAW force promotion.

For any proposed new port, use the following exact ownership rules:

- Retained external labels are (z,Y,N,G_h), and the declared node/origin/mode labels. N and G_h are unread by every future-history service.
- All F2/Delta slots within a single marked coefficient share its coherent ancestral bank. An independent copy may be introduced only at the location required by its actual Riesz/coarse construction, not by altering the target's slots.
- Distinct COMPLETE mean/covariance/skew/quartic/native service banks are fresh conditionally on the same Y. Their internal coarse, coefficient, shield, filter, path, side, complement and keep roots are integrated by their whole-law comparison.
- A caller or test matrix may depend on retained unread labels. It may not read a consumed bank. In particular it may not adapt to a private marked history or a covariance-field realization after that bank has been forgotten.
- The top keep variance is allocated once; every internal keep remains inside its own assigned share. Block readout [I I] must be included in the variance and radius normalization.
- Changed Y, shifted source argument, theta, history version, coefficient heat or source mode replays every affected F1/F2/Delta/J ancestor. Only an explicitly same-key unchanged VALUE site can be shared. Discarded primals and all first/adjoint ancestor sweeps remain charged.

For illustration, REPLACING rather than appending old coefficient services leaves the top pattern

    u_M=u_H3=u_K3=u_3,3=u_4,3=v/10, u_keep=v/2.

This is merely a valid TOP variance budget, not admission of the unspecified services. Appending a separate new service requires an explicit reallocation and recompilation of all changed shares; it cannot borrow the untouched v/2 without changing the consumer proof.

Let Q_H3,d_H3 etc. denote the fully expanded ORIGINAL-VALUE/root counts of actual future implementations, not unit-cost tensors. The honest prospective node bill is

    Q_next <= sum_LAW[Q_M[S3;u_M]+Q_H3+Q_K3+Q_3,3+Q_4,3
                       +J_bridge+1+Q_capture]
              +sum_RAW[Q_R,new+J_bridge+1+Q_capture]
              +Q_external.                            (19)

Here

    Q_M[S3;u] <= Q_cap+N_B(u) Q_B3
                  +N_E(u)(Q_S3+Q_B3)+Q_ext,

with native counts evaluated at the actual S3 dimensions, public logs, allocations and floors. The root bill before known-row alignment is the sum of every complete service dimension plus G_Y,G_keep,N,G_h; common-carrier alignment removes only duplicated known D-rows, retaining all perpendicular/service roots. A new joint-history compiler must supply its own literal Q and d recurrence before (19) has numerical or polylogarithmic content.

The sealed b=25, mu=full-radius setting and its A^(1/20), A^(3/40) margins certify its old actual sources only. They do not certify a history-marked composite vertex. Normalizing a putative Delta source by A^3 when its only known first is O(A) makes its normalized first O(A^-2). Thus a false difference grade can destroy native admission before any error ledger is considered.

For every actual proposed primitive occurrence l enumerate its downstream amplification P_l and choose absolute floors with sum P_l epsilon_l<=e_abs. Include new block readouts, inverse u, new shields, theta/clock gaps, coefficient heat and all changed ancestors. Substantive target differences (1), prefix errors, conditional law remainders and native priors cannot be moved into this floor budget.

No new cost theorem follows until Q_H3,Q_K3,Q_3,3,Q_4,3,Q_R,new and their native radii are constructed. In particular this packet does not claim fixed-log, any-order, sublinear-in-order, or endpoint-distribution complexity for an unspecified port.

## 9. Linear witness in full

For g(x)=A x let

    Z_j=integral_0^infinity e^-s s^(j-1)/(j-1)! X_s ds.

Then

    F_k=sum_(j=1)^k (-1)^(j-1) A^j Z_j,
    E[Z_j|X0=x]=x/2^j,
    Cov(Z_n,Z_m|X0)=[(binomial(n+m,n)-1)/2^(n+m)] I.    (20)

To verify the last equality split the double integral at s=t. The unconditional covariance is binomial(n+m,n)/2^(n+m); subtract the product 2^-n 2^-m of conditional means. For j<=3 the scalar conditional covariance matrix is

    [ 1/4    1/4     3/16  ]
    [ 1/4    5/16    9/32  ]
    [ 3/16   9/32   19/64  ].                          (21)

Thus

    Cov(F2|x)=[A^2/4-A^3/2+5A^4/16] I,
    Cov(F3|x)=[A^2/4-A^3/2+11A^4/16-9A^5/16+19A^6/64] I.

The canonical mean formulas give (2). These identities concern the TRUE force genealogy; substituting a common-carrier mean graph would change the covariance already at A^2 and would fail this test.

### 9.1 The next increment reverses sign, and the algebra repeats honestly

For every fixed history depth k>=2 the same proof, with Delta_k=F_k-F_(k-1), gives the rank-one through rank-four analytical envelope

    ||E[Delta_k|Y]||_2 <= A^k sqrt(D),
    ||Cov(F_k|Y)-Cov(F_(k-1)|Y)||_(2;HS)
        <=2 L A^k sqrt(D),
    ||kappa3(F_k|Y)-kappa3(F_(k-1)|Y)||_(2;HS)
        <=6 L^2 A^k sqrt(D),
    ||kappa4(F_k|Y)-kappa4(F_(k-1)|Y)||_(2;HS)
        <=24 L^3 A^k sqrt(D).                          (22)

This is a reusable target-restoration identity and difference bound, not a reusable accuracy-raising VALUE program. A marked slot has energy grade k and only the known first O(A); multiple marks still do not multiply energy grades.

There is also an exact obstruction to indefinitely adding independent POSITIVE covariance-increment banks. It uses the SAME original linear g(x)=A x and the next genuine force F4:

    Cov(F4|x)-Cov(F3|x)
       =[-A^5/4+7A^6/16-17A^7/32+69A^8/256] I.        (23)

For 0<A<=1/2, the bracket divided by A^5 is at most -1/32: the last two terms combine as A^2(-17/32+69A/256)<=0, and -1/4+7A/16<=-1/32. Thus this covariance increment is STRICTLY NEGATIVE. A fresh independent Gaussian or other independent positive-law bank can only add a PSD covariance and cannot implement (23).

The reference law for the FULL new covariance remains positive. For example the block-Gram (16)-(17) now uses (F3,F4-F3), and its [I I] contraction is exactly Cov(F4). The sign is handled by the coherent cross block, not by a negative probability measure or a negative Gaussian variance. This gives a same-g actual-history reason for preserving cross ancestry in any repeatable version of this architecture, beyond the generic port nonimplication in Section 3.4.

More generally, the leading linear adjacent covariance increment at depth k is

    (-1)^(k-1) k A^(k+1)/2^k times I,

as follows from 2 Cov(Z1,Z_k)=k/2^k in (20). Its sign alternates with depth. The associated next-mean increment m_(k+1)-m_k is (-1)^k A^(k+1)x/2^(k+1). These exact identities continue to hold when all higher cumulants are zero.

## 10. Conclusion

The next target is fully determined and contains a nontrivial four-force shifted mean chain, a grade-four cross covariance, and complete grade-five/six marked cumulant differences. Their exact algebra, positive reference, standard-Y scope and full source ownership are now explicit. The missing theorem is an executing, admitted history-marked original-VALUE supplier, with same-carrier feedback and cut grades. More unshifted force cumulant ranks or more native order alone do not supply that theorem. No independent accuracy-raising step has been proved here.

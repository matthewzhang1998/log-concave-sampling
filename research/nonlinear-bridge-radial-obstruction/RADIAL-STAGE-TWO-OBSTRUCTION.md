# A genuine C2 radial obstruction for the unchanged stage-two graph

2026-10-05. This packet adds a lower bound. It does not alter a sealed graph, execute a smoothing oracle, or claim an impossibility theorem for other algorithms.

## 1. Result

Throughout let 0<A<=1/36. Let mu_g=Q_out psi_Q,g be the own mean of the ORIGINAL finite stage-two VALUE graph. Its sources, resummed S, genuine conditional Gaussian rows, positive clock rules, and all recorded descendants are exactly those of STAGE-TWO-POSITIVE-VALUE-THEOREM.md. Assume the single and outer rules have their stated exact first exponential moments, and the quadratic-pair rule has its certified error deltaQ with the actual block parameter b=D. Linear-pair weights need only be positive and mass one for this lower bound.

There are explicitly specified anchored convex gradients g=grad U with U in C2, 0<=Dg<=A I, on the full D-dimensional space, for which

    ||mu_g-m3(g)||_(L2 gamma_D)
       >= [c_star-8A-deltaQ/4] A^2 sqrt(D)-7A,       (1.1)

    c_star=(17/sqrt(2)-11)/192
           =0.0053167462508922... .

The same first-chaos witness works for exact outer R1. No arbitrary numerical error is credited with improving the lower bound; an outer moment defect or any implementation floor must be subtracted explicitly.

In particular choose, for each integer n>=10000,

    A=1/n, D=(3000 n)^2, deltaQ<=A^2.

Then the source below is computable and belongs to the stated C2 class, and

    ||mu_g-m3(g)||2 >= (1/500) A^2 sqrt(D).           (1.2)

Thus this graph cannot have a dimension-uniform o(A^2)sqrt(D) certificate over the unrestricted original class. For a uniform monomial bound C A^p D^(1/2+gamma), with p>2 and C independent of A,D, this example requires

    gamma >= (p-2)/2.                               (1.3)

At p=8/3, at least an extra D^(1/3) is necessary. The current smoothing certificate pays D^(2/3); this result DOES NOT establish that its full exponent is necessary, nor construct a matching upper bound. It identifies a real graph obstruction beyond the use of a loose abstract tensor norm.

A second explicit sequence A=1/n,D=n^6, deltaQ<=A^2 has a first-chaos witness converging to c_star after division by A^2 sqrt(D). This is a witness asymptotic; no equality for the entire L2 norm is asserted.

## 2. The actual smooth radial source

Set

    alpha=beta=1/2,
    cv=alpha(1-1/sqrt(2)), c1=alpha/2=1/4,
    tau=cv/2, r1=1-A tau, eta=A^2.

For eta>0 define the continuously differentiable hinge

    H_eta(t)=0                         if t<=-eta,
             (t+eta)^2/(4eta)           if -eta<=t<=eta,
             t                         if t>=eta.

It satisfies 0<=H_eta'<=1 and 0<=H_eta(t)-max(t,0)<=eta/4. Put

    phi_A,eta(q)=alpha H_eta(q-1/2)+beta H_eta(q-r1),
    g_A,D(x)=A sqrt(D) phi_A,eta(|x|/sqrt(D)) x/|x|,
    g_A,D(0)=0,
    U_A,D(x)=A D integral_0^(|x|/sqrt(D)) phi_A,eta(s) ds.

This is an original source, not an analytical convolution oracle. Its value uses a norm, scalar arithmetic, and the displayed piecewise polynomial. It is constant zero in a neighborhood of the origin. Since phi is C1, U is C2. The radial and tangential Hessian eigenvalues are respectively

    A phi'(q), A phi(q)/q.

Both lie in [0,A]: the derivative bound follows from alpha+beta=1, and H_eta(q-r)<=q when r>=eta. Here both thresholds exceed eta for A<=1/36. The construction is full-dimensional and rotationally equivariant. No block decomposition or special basis is supplied to the graph.

For analysis only, replace H_eta by t_+ and call the resulting Lipschitz convex gradient g^0. Its potential need not be C2 at the shells, but

    sup_x |g_A,D(x)-g^0_A,D(x)| <= A^3 sqrt(D)/4.     (2.1)

The final source in (1.1) is always the displayed C2 source. Section 5 transfers the calculation back to it by a strong coherent-source estimate, not by freezing descendants.

## 3. Dimension-uniform reduction to scalar Gaussian rows

Write phi(q)=alpha(q-1/2)_++beta(q-r1)_+ for this section. If Y is a centered scalar Gaussian row with covariance q^2 I_D, compare g^0(Y) with

    T(Y)=A [phi(q)/q] Y,

with the zero row interpreted as zero. Since 0<=phi'<=1 almost everywhere and 0<=phi(q)/q<=1, the scalar function

    r -> sqrt(D)phi(r/sqrt(D))-[phi(q)/q]r

vanishes at q sqrt(D) and is 1-Lipschitz. Consequently

    ||g^0(Y)-T(Y)||2
       <= A || |Y|-q sqrt(D) ||2 <= sqrt(2) A q.     (3.1)

The last estimate follows by | |G|-sqrt(D) |<=||G|^2-D|/sqrt(D) and E(|G|^2-D)^2=2D. This calculation is valid even though the source itself varies with A,D.

Apply the deterministic radial map T to coefficient rows in their Gaussian Hilbert space, through EVERY original VALUE site. That radial map is A-Lipschitz in any Hilbert dimension. All reference input row norms are below 2 at A<=1/36: v has norm 1/sqrt(2), w has norm k=sqrt(3/8), all U,V have norm 1, and a largest terminal norm is bounded by

    1+(2+3/sqrt(2))A+k A^2 < 2.

The coherent source-error propagation from the sealed graph is

    K_graph=7/2+14A+(11/2)A^2 <4.

It applies just as well to L2 local VALUE defects at these Gaussian reference inputs. It includes g(v), g(w), g(g(w)), every g(U), g(V), g(U-g(V)), and all signed terminal sites. Positive averages prevent any factor equal to a clock node count. Equations (3.1) and this propagation give a graph/reference error below 12A in L2. This remains true when the quadratic-pair average is an exact probability integral.

For the genuine target retain the ENTIRE original OU future. Let H1,a and H2 denote its original exponentially weighted linear histories. Every X_a is marginally standard. Write

    cU=phi(1)=c1+beta A tau,
    q2=sqrt(1-A cU+(A^2/2)cU^2), c2=phi(q2)/q2.

Minkowski and (3.1), with the A-Lipschitz propagation at the next source, give

    F1,a = A cU H1,a + L2 error <=sqrt(2) A,
    F2   = A c2 (v-A cU w)
                 + L2 error <=sqrt(2)A(1+2A).

The exact identity integral exp(-a)H1,a da=H2=w is used here; no independent-time surrogate replaces the genealogy. One more original source application gives

    g^0(x-F2) = T(x-A c2 v+A^2 c2 cU w)
                   + L2 error <2A.

Thus the full graph-minus-genuine-target first-chaos calculation below has L2 reduction error less than 14A, uniformly in D and all finite single/linear node lists. It is not an asymptotic appeal to concentration with an uncontrolled growing number of nodes.

## 4. Leading defect on the actual rows

For a Gaussian reference row

    Y=x+A y+A^2 z, ||y||_row,||z||_row<=1,
    yx=Cov(x,y) in any coordinate,

put q=||Y||_row and n=Cov(x,Y). For A<=1/36,

    |q-1-A yx|<=2A^2,  0<=1-n/q<=A^2.

Indeed ||A y+A^2 z||<=A(1+A), and q-n is the squared norm of the perpendicular component divided by q+n. Here q+n>1.9 and q>0.97. Hence the normalized first-chaos coefficient of T(Y) is

    E<x,T(Y)>/D
       =A c1+A^2[alpha yx+beta(yx+tau)_+]+r,
    |r|<=3A^3.                                    (4.1)

The lower-shell amplitude is alpha(q-1/2); the upper-shell amplitude is beta(q-1+A tau)_+. Apply their Lipschitz bounds and the two preceding inequalities. The factors multiplying the O(A^2) errors are bounded by 3 after summing alpha+beta=1.

The original retained rows are v=x/2+N/2 and w=x/4+N/2+M/4. In the reference graph, g(v)=A cv v. The nested g(g(w)) is exactly zero, because g(w)'s row norm is below 1/2. Therefore

    S=x-A cv v.

At a one-time row U_s, put r=exp(-s). It has its ORIGINAL covariance with x,v,w, not a replacement history. The single displacement is

    d(U)=A cv v-A(c1+beta A tau)U.

Its first-order radial projection is

    a(r)=cv/2-c1 r.

The signed leading terminal projections relative to the upper shell are +a(r),-a(r) for S+/-d. The reference S itself is exactly on that shell to first order, because tau=cv/2. The finite single rule has mean r=1/2, so its upper-shell centered contribution is EXACTLY

    beta A^2 E[a(r)]/2,

regardless of its node count or nonlinear quadrature properties. No nonlinear one-time operator estimate is being assumed.

For Q, use the genuine ordered clock earlier a~Exp(2), gap h~Exp(1). Its projections are a(exp(-a)),a(exp(-(a+h))). Q is symmetric in its two displacement slots. Splitting two iid Exp(1) times by order therefore gives the identical expectation with two independent r,s~Uniform[0,1]. Its upper-shell leading coefficient is

    (beta A^2/8) E[|a(r)+a(s)|-|a(r)-a(s)|].

This uses the true pair law; the retained Gaussian row correlations are all still present. The leading radial projection happens to depend only on the x coefficients.

Set u=cv/(2c1)=1-1/sqrt(2), so a(r)=c1(u-r), with 0<u<1/2. Elementary integration gives

    E a =c1(u-1/2),
    E|a(r)-a(s)|=c1/3,
    E|a(r)+a(s)|=c1[1-2u+(8/3)u^3].

The total upper-shell coefficient is consequently

    beta c1[u/4-1/6+u^3/3] A^2
       =(11-17/sqrt(2)) A^2/96
       =-2 c_star A^2.                             (4.2)

The lower-shell contributions cancel to order A^2 against the genuine target. In detail g(S) contributes A c1-alpha cv A^2/2; the single correction contributes alpha(cv-c1)A^2/2. Their sum is A c1-alpha c1 A^2/2, exactly the target lower-shell coefficient. Q annihilates the affine lower-shell coefficient.

For the genuine target, its row is x-A c2 v+A^2 c2 cU w. Here cU<0.252, q2<r1, |c2-c1|<0.064A, and c2<=1/4, so its y,z satisfy the bounds in (4.1). Its leading radial projection is -c1/2; since tau<c1/2 the target has zero upper-shell leading coefficient. This target calculation uses the actual two-level genuine force from Section 3.

All g(S), single, and Q terminal rows also satisfy ||y||,||z||<=1: the largest y norm is at most 3cv/sqrt(2)+2c1<0.819 and the largest z norm is at most 2 beta tau<0.075. Their total absolute terminal coefficient mass is 1+1+1/2=5/2. Their (4.1) errors contribute at most (15/2)A^3. The target contributes at most 3A^3. The linear-pair centered term has row norm at most (1+k)A^3<2A^3, including its entire nested source and any positive mass-one linear clock rule. The resulting reference first-chaos defect is therefore

    -2 c_star A^2 + error of absolute value <13A^3.  (4.3)

## 5. Finite pair rule, C2 transfer, and outer witness

First replace ONLY the finite quadratic-pair sum on the smooth source by its exact ordered probability integral. The sealed general-b=D pair theorem bounds the inner mean difference by

    (deltaQ/2) A^2 sqrt(D).

The proof uses its actual C_D and finite node count; neither is silently dimension independent. All other finite nodes are left in place.

Next use (2.1) to pass from the smooth source to the hinge source throughout that hybrid graph AND throughout the genuine target. For coherent sources of Lipschitz constant A, separated by sup norm delta, the graph costs K_graph delta, while the genuine target costs (1+A+A^2)delta. Their sum is

    [9/2+15A+(13/2)A^2]delta <5delta.

This explicitly includes the changed descendants. With delta=A^3 sqrt(D)/4 the combined transfer is below (5/4)A^3 sqrt(D).

For every vector f in L2(gamma_D) and standard Gaussian Z, stationarity and Gaussian regression give

    E<Z,R1 f(Z)>=(1/2) E<X,f(X)>.

The finite outer rule has EXACTLY the same identity because sum_i omega_i t_i=1/2. Consequently there is no outer-rule weakening of this first-chaos witness, and outer endpoint singularities are irrelevant. Cauchy-Schwarz gives

    ||mu_g-m3(g)||2 >= |E<Z,mu_g(Z)-m3(g)(Z)>|/sqrt(D).

Combine (4.3), the <14A radial reduction, the (5/4)A^3 sqrt(D) coherent smoothing transfer, the quadratic-pair allowance, and the outer factor 1/2. The A^3 coefficient is (13+5/4)/2=57/8<8. This proves (1.1).

For n>=10000 and D=(3000n)^2 the lower ratio is at least

    c_star-8/n-1/(4n^2)-7/3000 > 1/500.

For a purely rational check, sqrt(2)<99/70 implies c_star>53/10000; this already proves the strict 1/500 inequality above. This proves (1.2). Substituting this dimension sequence into any proposed uniform monomial upper bound proves (1.3).

For D=n^6, the errors in the first-chaos witness divided by A^2 sqrt(D) vanish: the expansion/smoothing cost is O(1/n), the Gaussian reduction O(1/n^2), and the pair cost O(1/n^2). Thus the absolute first-chaos witness tends to c_star, as stated.

## 6. What remains unchanged, and what this does not prove

This is a counterexample for the existing graph's own mean relative to its genuine-history canonical m3. There is no new VALUE site, smoothing grid, root, root dimension, HVP action, or source-query multiplier. Every retained shift and descendant in the original graph is preserved. On the displayed smooth source, the original Hessian interval holds, so the already-audited algebraic first/curl/source-energy bounds apply with their original constants and all original live caller/anchor/scale dependencies. This packet does not independently re-audit imported native compilers, or assert that their caller/radius/clock/filter/precision/mode guards hold without checking them.

The original private dimension is 5D and its full original-VALUE bill remains

    Q_rule/certificate_setup + Q_captured + N_out N_B
       +(5+3J1+5JL+6JQ) N_out N_E
       +Q_known/numerical/replay.

The actual general-class pair setup uses b=D and its huge dimension-dependent C_D. The lower bound assumes that certified rule, rather than replacing it by a cheaper unproved one. Numerical coefficients, moments, source values, row computations, caller captures, original HVPs, and full replay keep every original absolute floor. The explicit source formula is computable; this does not make the existing finite rule or compiler practically cheap.

If an actual guarded completed service is within epsilon_comp, in integrated W2, of N(mu_g(Z),I), the triangle inequality gives its distance to N(m3(g)(Z),I) at least (1.1) minus epsilon_comp. In particular a genuinely smaller certified completion allowance cannot repair this bias. No lower bound for a completed program is claimed when its guard or numerical/completion allowance has not been verified.

The example does not show that the Pinsker/transport sqrt(D) loss alone is sharp for coherent gradient differences. It does show that improving that estimate and/or the shifted remainder cannot lead to a dimension-uniform improved power of A for this same graph over the entire C2 class. A change of graph, additional history cancellation, a source qualification, or an exposed dimension regime would be needed for that claim.

The D^(2/3) extra factor in the current A^(8/3) certificate remains open between the necessary D^(1/3) and that existing upper bound. This packet does not complete an any-order recurrence, an endpoint join, an algorithm-wide lower bound, or any old Picard/cold-start route.

## 7. Pinned inputs

Smoothing manifest SHA256:
27eb73bef8e9e3fa7461dda875ecdf4c9d488b1a04ddbcc0765348a16c944417

Smoothing C2 transfer SHA256:
245bc253a619452ab3112b8e40efcbcb298089756c7b170a47a4a77a57526498

Stage-two manifest SHA256:
cd1111e28dc399643f27bcdd0f97dc735ada196f8a43f220211500bcb45be3ae

Stage-two theorem SHA256:
c7275de9261bce2eeaee1ba054351055b1e236590340ee7b1554ff064ffa2f13

Stage-two shifted derivation SHA256:
9650c96e96a24e64e4afa4da183b652bb5dbc410040b69e40d20783058acd682

Stage-one manifest SHA256:
5d030acf6f366539239ddfea9dca3e6b5c9b05adc323f2432cf49ecc266c161d

Stage-one theorem SHA256:
301ace0ae1a52f1cc62193c716702e0d908affab149ffd147082b23ad4d5637e

No sealed input file is modified. The separate INDEPENDENT-RADIAL-AUDIT.md accepts the proof, explicit constants, both dimension sequences, and the monomial dimension-exponent floor. Author and independent scripts provide diagnostics only; neither executes the imported completed compiler.

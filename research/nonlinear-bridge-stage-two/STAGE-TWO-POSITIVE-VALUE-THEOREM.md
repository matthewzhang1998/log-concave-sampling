# One finite nonlinear stage-two VALUE repair of canonical m3

2026-10-05. A bounded order-four result on an explicitly smoother nonlinear class. This is one additional accuracy step, not an induction or an unrestricted-C2 theorem.

## 0. The result and its precise qualification

The original source is anchored and convex-gradient:

    g(0)=0, 0<=Dg<=A I, 0<A<=1/2.

For the order-four target certificate, assume that in a fixed orthogonal decomposition into blocks of dimension at most b, g acts blockwise and, globally on each block,

    Lip(Dg_block)<=B A,
    D²g_block is continuous and Lip(D²g_block)<=C A.

All tensor norms are induced Euclidean multilinear norms. The constants b,B,C are declared source qualifications; b is fixed independently of D and A. The executing graph does not read the block basis. These assumptions include the sine source and bounded-size genuinely noncommuting nonlinear blocks. They are stronger than the stage-one curvature qualification and much stronger than unrestricted C2. Their verification is not inferred from finitely many black-box VALUES.

The construction below preserves the original resummed shift S=x-g(v)+g(g(w)), adds two concrete original-VALUE corrections, and has

    integrated W2 error to N(m3(Z),I)
      <= [K_(b,B,C,A) A^4 + delta1 A^2
           +sqrt(11/8)deltaL A^3+(deltaQ/2)A^2
           +delta_out A C_F2+Lambda_comp A^4]sqrt(D)
           +restored absolute floors.

Here m3 is the genuine-history canonical target. The pair rules use the explicit normalization constant C_b in Section 3, so their displayed deltaL,deltaQ errors already include that factor. Lambda_comp is the collected guarded completion constant, including the numerical factors of the normalized ports below; it is not silently identified with an old sharp native prefactor. Set

    k=sqrt(3/8), d0=1+1/sqrt(2), e0=1+k,
    u_A=d0+A e0,
    K_(b,B,C,A)
      =B[d0 e0+(A/2)e0²+1/2]sqrt(b+2)
        +(C/6)[u_A³+4d0³+A³e0³]sqrt((b+2)(b+4)),
    C_F2=1+A[3/2+5/(2sqrt(2))+A(1+2k)].

Choosing delta1<=A², deltaL<=A, deltaQ<=A², and delta_out<=A³ gives order-four target accuracy. The pair quadrature estimate itself needs no B or C; these enter the shifted target-bias theorem. All these choices have polylogarithmic node counts and zero inverse-A original-VALUE exponent at this fixed grade.

The literal residual has 5D private Gaussian coordinates and costs

    (5+3J1+5JL+6JQ)N_out

original VALUES per fully replayed occurrence. J1=O(log²(1/delta1)), JL,JQ=O(log⁴(1/delta)), and N_out is the pinned outer positive Hermite-rule count. The safe half-variance normalized private first is at most 9A and curl is below 25A². Native compiler guards, caller captures, and all replays remain mandatory.

## 1. The actual original-VALUE graph

Draw fresh independent standard D-Gaussians N,M,L1,L2, and retain

    v=x/2+N/2, w=x/4+N/2+M/4,
    b0=g(g(w)), S=x-g(v)+b0.

The unchanged centered primitive and the new polarized primitive are

    C0(S,a)=[g(S+a)-g(S-a)]/2,
    Q0(S,a,b)=[g(S+a+b)+g(S-a-b)
                 -g(S+a-b)-g(S-a+b)]/8.

Use three finite positive probability rules:

1. The sealed one-time conditional bridge rule (p_j,s_j), with operator error delta1, mass one and sum p_j exp(-s_j)=1/2.
2. An ordered-pair tensor rule for Exp(1)(a) times Exp(1)(h), with error deltaL and exact coordinate exponential moments 1/2,1/2.
3. An ordered-pair tensor rule for Exp(2)(a) times Exp(1)(h), with error deltaQ and exact coordinate exponential moments 2/3,1/2.

For each single node, U_s has the sealed exact joint law with x,v,w and may use L1. For each pair node, generate U,V with the TRUE conditional law of X_a,X_(a+h) given x,v,w, using Section 2. Define

    d(U)=g(v)-g(U),
    e(U,V)=g(U)-g(U-g(V))-b0,
    F2_node=g(S)
        +sum_single p_j C0(S,d(U_j))
        +sum_linear_pair q_j C0(S,e(U_j,V_j))
        +sum_quadratic_pair r_j Q0(S,d(U_j),d(V_j)).

All nodes share N,M,L1,L2 within one occurrence. Sharing these roots across nodes does not assert that their collection is a genuine multi-time OU path. Each node's required one-time or pair law is exact; expectation linearity is all that is used across nodes.

For the pinned outer rule (omega_i,t_i), let c_i=sqrt(1-t_i²), draw another independent D-Gaussian G, and use the SAME G,N,M,L1,L2 at all outer nodes:

    x_i=t_i Z+c_i G,
    F2=sum_i omega_i F2_node(x_i,N,M,L1,L2),
    B_raw=sum_i omega_i g(x_i), E2=F2-B_raw.

The exact own mean is Q_out psi_stage2(Z), where psi_stage2 is the mean of the literal inner finite graph. There are no expectation, derivative-action, or continuous-history producer leaves.

## 2. Exact two-time conditional Gaussian rows

For a,h>0 put

    sigma_s=sqrt(1-exp(-2s)), f(s)=2s/sqrt(exp(2s)-1),
    R(a,h)=[[f(a),f(a)(a-1)],
            [exp(-a)f(h),exp(-a)f(h)(2a+h-1)]],
    K=I2-RR^T.

The two standardized OU Markov innovations have cross-covariance R with (N,M). A fully rational certificate plus analytic tails proves

    I2/16 <= K <= I2

for every positive a,h. Execute only the known scalar 2-by-2 root

    Kroot=(K+sqrt(det K)I2)/sqrt(tr K+2sqrt(det K)),
    (W1,W2)=R(N,M)+Kroot(L1,L2),
    U=exp(-a)x+sigma_a W1,
    V=exp(-h)U+sigma_h W2.

This is the genuine ordered pair conditional on x,N,M. Each full U or V scalar row in (x,N,M,L1,L2) has Euclidean norm one. The uniform K gap avoids a hidden inverse small-time conditional variance in the root calculation. Only known scalar functions and a fixed-size source-independent matrix are used; no D-by-D matrix or source-dependent covariance root is reconstructed.

The exact rational certificate is certify_pair_covariance.py: 2,792 accepted boxes on [0,8]² prove K>=I/10 there; explicit exponential tails imply the global I/16 bound. The complete argument is TWO-TIME-POSITIVE-QUADRATURE.md.

## 3. Why the pair rules are finite, positive, and polylogarithmic

For each declared block, the conditional Gaussian density of U,V is holomorphic in both positive time variables on uniform relative complex disks. Only the density is complexified. The original g, retained S(x,N,M), and all VALUE arguments under integration remain real; no source analyticity is assumed.

The pair proof supplies a certified radius eta=2^-24 after normalizing by the real Markov innovations. The covariance stays uniformly close to a real positive matrix, and the absolute complex density is bounded by an integrable Gaussian times exp(2^-23|Y|²), where Y=(x,N,M) on that block. Its square is integrable. Polynomial-growth real integrands therefore produce uniformly bounded L2(Y)-valued holomorphic functions.

On each time axis, use positive Gaussian quadrature for the exact exponential measure on dyadic panels, plus positive small head/tail atoms. There are O(log(1/delta)) panels and O(log(1/delta)) nodes per panel. Tensoring the two axes yields O(log⁴(1/delta)) pair nodes. Positive one-atom moment repairs preserve exact mass and the stated coordinate exponential moments. The construction uses product weights, so its required joint moments are exactly

    linear pair: E exp(-(a+h))=1/4,
    quadratic pair: E[exp(-a)+exp(-(a+h))]=1.

These identities would not follow merely from unrelated marginal moment claims for a non-product rule.

For R5=sqrt(|x|²+|N|²+|M|²+|U|²+|V|²), the actual integrands have GLOBAL LINEAR bounds

    |C0(S,e(U,V))|<=sqrt(11/8)A³ R5,
    |Q0(S,d(U),d(V))|<=(A²/2)R5.

For a block of dimension d<=b, the uniform complex L2 envelope is at most C_b H sqrt(d) whenever the real integrand is bounded by H R5, with

    C_b=12*8^b*(1-2^-21)^(-(3b+2)/4).

Run the rule at internal tolerance delta/C_b to obtain error <=delta H sqrt(d), giving the precise pair errors in Section 0 after summing squared block norms. It is not a claim that the sealed one-time operator rule already controls these nonlinear stencils.

For clarity, actual finite node counts can be bounded explicitly. For either pair tolerance delta, put epsilon=delta/(32C_b), ell=epsilon/8, Htime=log(8/epsilon), rho=1+2^-24, and

    n=ceil(log(128/[epsilon(rho-1)])/[2log rho]),
    J_marginal<=n ceil(log2(Htime/ell))+3,
    J_pair<=J_marginal².

The last three marginal nodes cover the two head/tail atoms and the optional moment-repair atom. These conservative constants are VERY large; the theorem makes a finite fixed-grade and zero inverse-A-exponent claim, not a practical runtime claim. Every actual node is still included in the occurrence bill.

## 4. The specific order-four cancellation

On the genuine history write

    F1=int exp(-t)g(X_t)dt,
    F2=int exp(-a)g(X_a-F1,a)da,
    d=g(v)-F1,
    e=F1-F2-b0.

Then x-F2=S+d+e exactly. The first repair matches Dg(S)d. The new C0(S,e(U,V)) matches Dg(S)e up to order four: expand the inner g(X_a-g(X_(a+h))) at the actual X_a-F1,a and use the pathwise variance identity. Its mean error before the terminal derivative is at most (B/2)A³ sqrt(d_block(d_block+2)). No unshifted covariance proxy is used.

The new Q0 term matches (1/2)D²g(S)[d,d]. The exact rectangle identity

    Q0(S,a,b)=(1/8)int_[-1,1]^2 D²g(S+s a+t b)[a,b] ds dt

gives its leading term and a quantitatively controlled remainder. Splitting two independent Exp(1) history times by their order gives earlier time Exp(2) and gap Exp(1), exactly explaining the second pair rule. Independent one-time innovations would give the wrong cross moment.

All Taylor expansions are at S or at X_a-F1,a. The source assumptions bound the remaining cross terms and third-order remainders, yielding the explicit K_(b,B,C,A) in Section 0. The complete derivation is SHIFTED-BIAS-AND-PORTS-DERIVATION.md, Sections 1–2. No third derivative is queried; D²g and its Lipschitz bound occur only in analysis.

For g(x)=Kx, Q0 is identically zero and each new linear-pair correction is K³(V_j-w). Its WEIGHTED mean is zero because sum_j q_j E[V_j|x]=x/4 and E[w|x]=x/4. A fixed pair node instead has E[V_j|x]=exp(-(a_j+h_j))x; it is not individually x/4. Consequently all anisotropic quadratics retain exact canonical m3 means, including the exact outer first moment. The actual raw covariance is allowed to change; the consumer uses the literal source's own mean.

The formerly obstructing sine family g_A(x)=A(x+sin x)/2 satisfies b=1,B=C=1/2, so this packet raises its certified target grade from the original A² obstruction through the first repair's A³ bound to A⁴. A genuinely noncommuting two-dimensional example is

    g(z)=A[z/2+(sin(a.z)a+sin(b.z)b)/8],

with nonparallel nonorthogonal unit vectors a,b. It has A I/4<=Dg<=3A I/4 and B=C=1/4. Arbitrary fixed orthogonal sums and rotations of such blocks remain covered without disclosing their basis to the executor.

## 5. Complete first, curl, caller and energy ports

The exact noncommuting derivative identities are recorded in SHIFTED-BIAS-AND-PORTS-DERIVATION.md, Section 4. They only require the original PSD Hessian interval, not B,C or the block qualification. In particular the nested e derivative is

    e_y=[Dg(U)-Dg(U-g(V))]U_y
          +Dg(U-g(V))Dg(V)V_y-Dg(g(w))Dg(w)w_y,

with the factor order retained. The aggregate coefficient of S_y is symmetric with norm, and difference from Dg(x), at most 9A/4. Put

    eta=A/2+A²/4, q=A/2+A²/2, r0=A²/4,
    Cx=9(1+eta)/4+A(13/4+5A/4),
    Cn=9q/4+A(13/4+3A/2),
    Cm=9r0/4+A(5/2+5A/4),
    Cl=A(5/2+A), beta=sum omega_i c_i<=sqrt(3)/2.

Then

    raw private first <=A sqrt(beta² Cx²+Cn²+Cm²+2Cl²),
    raw retained-Z caller first <=A Cx/2,
    raw lifted curl
       <=2beta[(9/4)A eta+A²(13/4+5A/4)]
           +A sqrt(Cn²+Cm²+2Cl²).

The square lift is now E2 -> (E2,0,0,0,0). Its leading G-block difference is symmetric; every other skew contribution is an explicitly small product. No Hessians are commuted. Under half-variance normalization these give first<=9A and curl<25A²<(9A)(3A).

The conditional energy is

    ||E2||_(Lp|Z=z)
      <=A²[(1+3A/8)|z|
         +kappa_p(3/2+5/(2sqrt(2))+A(1+2k))sqrt(D)].

In particular Q0 is bounded by A(|d1|+|d2|)/4, so this is a linear-growth source-energy estimate, not a hidden quadratic caller envelope.

Capture the ACTUAL caller origin at G=N=M=L1=L2=0, including every nonzero U,V descendant. Its magnitude is <=A²(1+3A/8)|z| and its cost is the full residual occurrence bill before exact-key sharing. Subtract and restore it. Anchoring leaves private first and curl unchanged and gives raw retained-Z caller bound A Cx, normalized bound sqrt(2)A Cx. The baseline's normalized anchored retained-Z bound is sqrt(2)A. Native caller guards use those normalized values. Any additional exterior source/anchor/scale parameters retain their actual recorded original-source chains and their own guards; these Z bounds are not substituted for arbitrary exterior-parameter bounds.

At z and all private roots zero every VALUE descendant is literally zero. At a nonzero caller, its captured origin is not erased. At g=0 the raw source vanishes identically and the imported completed mean service retains its literal Gaussian zero-source row. No normalization divides by residual energy, B,C, or a measured quantity that might vanish.

## 6. Positive own-mean completion and full original-VALUE bill

Use the same pinned gradient and near-gradient finite own-mean compilers and coisometry adapter as the sealed first-repair packet. The baseline and residual use independent COMPLETE banks conditional on Z and every actual exterior label, with half variance shares. Declare

    gradient order >=4,
    ell_E=9A, a_seed=3A, padding mu=A,
    residual private dimension=5D.

A<=1/36 suffices ONLY for ell_E<=1/4. Every native radius, curl, caller, active-dimension, clock/filter, mode, tolerance and precision guard remains imposed at the actual parameters. With energy O(A²)(|z|+sqrt(D)), the imported near-gradient bracket is

    36A²+729A³+729A^(5/2)=O(A²),

and therefore gives O(Lambda_comp A⁴ sqrt(D)) after integration. The gradient branch retains its separate order-four allowance. Product coupling is applied only to completed independent outputs. They target N(Q_out psi_stage2(Z),I); no raw private root, source record or trained statistic is appended as an observer.

The safe complete original-VALUE bill is

    Q_rule/certificate_setup+Q_captured+N_out N_B
       +(5+3J1+5JL+6JQ)N_out N_E
       +Q_known/numerical/replay.

N_B,N_E include every original occurrence, later replay, bank, filter, clock, pair, mark and finite-mode anchor in the fully expanded native programs, at the NEW normalized radii and dimension 5D. Caller captures and repeated captures are charged. Each derivative/adjoint sweep uses one original HVP at every recorded original VALUE site; changed inputs or coefficient versions replay the whole graph, and discarded primals pay complete replay. No HVP is differentiated. Exterior source parameters and original caller anchors remain live.

Every output is an ordinary finite Gaussian pushforward. Positive quadrature weights alone would not certify the completed law: positivity additionally uses the native positive variance shares, fills and buffers through this exact guarded reentry. Signed arithmetic in C0,Q0 is not a signed measure. Vector arithmetic costs O(D) per original leaf apart from the source's own cost; there is no dimension-independent arithmetic claim.

## 7. Absolute numerical floors

For uniform absolute original-VALUE error nu at the actually requested sites, a conservative raw F2 floor is

    (7/2+14A+(11/2)A²)nu.

E2 adds one baseline nu, and separately executing the captured origin doubles the residual floor before exact-key reuse. This accounts for every nested descendant, including the g(U-g(V)) action. It does not divide by energy or A-dependent measured norms.

More detailed propagation uses Delta_b<=nu_b+A nu_w and Delta_S<=nu_v+Delta_b. The total terminal sensitivity to Delta_S is at most 7A/2. Besides terminal original-VALUE errors, the remaining terms are bounded by 2A nu_v+A Delta_b, the averaged single-node A nu_U, the linear-pair A(nu_U+nu_shift+A nu_V), and the quadratic-pair (A/2)(nu_U+nu_V). Together these yield the displayed uniform floor.

For a pair node, Euclidean absolute scalar-row errors epsilon_U,epsilon_V contribute integrated floors at most

    linear-pair C0: A²(epsilon_U+A epsilon_V)sqrt(D),
    quadratic-pair Q0: (A²/4)(epsilon_U+epsilon_V)sqrt(D),

weighted by the node probabilities and then summed. Single-row error gives A² epsilon sqrt(D). Absolute weight perturbations have floors

    single: A²d0 sqrt(D)||Delta p||_1,
    linear pair: A³e0 sqrt(D)||Delta q||_1,
    quadratic pair: (A²d0/2)sqrt(D)||Delta r||_1.

Check positivity, normalized mass and native source ports on the actual numerical coefficient version. The stable K root has a certified fixed gap; do not silently clip an erroneous covariance. Any first-moment defects, original-root errors, weight/node/arithmetic errors, finite-mode/clock/filter/precision residuals and replays retain their absolute allowances. In particular a single-rule exponential-moment defect epsilon1 costs up to A²epsilon1 sqrt(D) in a quadratic raw mean; a linear-pair product-moment defect epsilonL costs up to A³epsilonL sqrt(D). Mathematical exactness is not declared for a finite-precision rule with unpriced moment defects.

## 8. Conditional normalization and the remaining boundary

The global block qualification has a useful, explicit normalization rule. For s>=0 and a live exterior anchor a, let

    f(y)=s[g(a+s y)-g(a)], alpha=A s².

The same block decomposition is retained and its valid declared constants are

    Lip(Df_block)<= (B s)alpha,
    Lip(D²f_block)<= (C s²)alpha.

For s<=1 these do not enlarge B,C. This only applies because the present bounds are GLOBAL; no localized certificate is being translated without proof. Every f VALUE must execute and charge its actual original-g sites, including g(a) as an actual reusable caller capture, with all a,s dependencies live. This exports only the structural qualification: applying the theorem to f still requires the actual alpha, every new caller/anchor/scale and native compiler guard, and the full original-g implementation bill to be checked. The rule/certificate/numerical versions are fixed for a requested sweep, while the actual a,s source dependencies remain live. This normalization fact is not an induction over the new history stencil.

The theorem proves one actual second repair of this canonical-m3 mean-law component, on its explicit smoother block class. It does not complete the full fourth-order endpoint join, show that repeating the correction gives further powers of A, remove growing smoothness requirements, or supply an unrestricted high-dimensional C2 rate. The sealed first-stage fixed-(D,h) little-o corollary remains qualitative and cannot be advertised as a uniform computable any-order mechanism. The old five-residual certificate is not automatically inherited at these new shifted VALUE sites; use a fixed structural certificate branch if retaining that separate old class.

The next unresolved recurrence question begins AFTER this finite order-four step: identify and cancel the next genuine conditional ancestry/cumulant debt, prove its actual finite positive joint-history rule, and export its own new source ports and smoothness costs. No such induction is claimed here.

## 9. Provenance and checks

Sealed stage-one theorem SHA256:
301ace0ae1a52f1cc62193c716702e0d908affab149ffd147082b23ad4d5637e.

Sealed stage-one manifest SHA256:
5d030acf6f366539239ddfea9dca3e6b5c9b05adc323f2432cf49ecc266c161d.

Author raw-graph diagnostic check_stage_two_source.py checks 8,170 assertions, including genuine pair covariances, noncommuting finite-difference firsts, normalized first/curl bounds, exact VALUE counts, coherent zeros, anisotropic quadratic means and absolute numerical floors. These diagnostics do not implement or re-certify the native completed mean compilers. The new independent audit is maintained separately in this packet.

# Stage-two shifted cancellation and actual source ports

2026-10-05. Independent derivation for the specific VALUE stencil. The positive pair-quadrature theorem is an external input to this document, not a claim proved here. No sealed first-repair file is modified.

## 1. Definitions and hypotheses

Let g act in a fixed orthogonal decomposition into blocks of dimension d<=b, with g(0)=0, symmetric 0<=Dg<=A I, 0<A<=1/2, Lip(Dg_block)<=B A, and continuous D²g_block with Lip(D²g_block)<=C A. All tensor norms below are Euclidean induced multilinear norms. In particular ||D²g_block||<=B A.

Retain exactly

    v=x/2+N/2, w=x/4+N/2+M/4,
    b0=g(g(w)), S=x-g(v)+b0.

The VALUE primitives are

    C0(S,a)=[g(S+a)-g(S-a)]/2,
    Q0(S,a,b)=[g(S+a+b)+g(S-a-b)
                 -g(S+a-b)-g(S-a+b)]/8.

On the genuine stationary OU history conditioned on x,v,w, write

    F1=int exp(-t) g(X_t) dt,
    F1,a=int exp(-h) g(X_(a+h)) dh,
    F2=int exp(-a) g(X_a-F1,a) da,
    d=g(v)-F1,
    e=F1-F2-b0.

Then x-F2=S+d+e exactly. The canonical inner target is psi2(x)=E[g(x-F2)|x].

At a genuine ordered pair U=X_a, V=X_(a+h), define

    d1=g(v)-g(U), d2=g(v)-g(V),
    e12=g(U)-g(U-g(V))-b0.

The exact analytical correction is

    g(S)+E_(a~Exp1) C0(S,g(v)-g(X_a))
        +E_(a~Exp1,h~Exp1) C0(S,e12)
        +E_(a~Exp2,h~Exp1) Q0(S,d1,d2).

These are analytical expressions only. An executable source replaces every integral with its finite positive rule and each time pair by the actual jointly conditioned Gaussian pair. A pair is not constructed from independent one-time innovations. Sharing residual roots across different nodes does not assert a joint multi-time history law.

## 2. The exact order-four bias constant

Put

    k=sqrt(3/8), d0=1+1/sqrt(2), e0=1+k,
    u_A=d0+A e0,
    K_(b,B,C,A)
      = B[d0 e0+(A/2)e0^2+1/2] sqrt(b+2)
        +(C/6)[u_A^3+4 d0^3+A^3 e0^3]
                         sqrt((b+2)(b+4)).

For exact analytical pair integrals and either the exact single integral or a finite positive single rule, the own inner mean obeys

    ||psi_stage2-psi2||_(L2 gamma_D)
      <= [K_(b,B,C,A) A^4 + delta1 A^2] sqrt(D),

where delta1 is the sealed one-time Gaussian-operator rule error. With finite pair rules, add their certified errors for the actual C0(S,e12) and Q0(S,d1,d2) integrands. This display does not replace those pair errors by an unproved rate.

Here is a proof using only Taylor expansions at the actual shifted S and at the actual nested mean X_a-F1,a.

For a block of dimension d, write

    sigma4=[d(d+2)]^(1/4),
    sigma6=[d(d+2)(d+4)]^(1/6).

Gaussian marginal laws and Minkowski imply, for every finite p>=2,

    ||d||p <= A d0 sigma_p,
    ||e||p <= A^2 e0 sigma_p,
    ||d+e||p <= A u_A sigma_p,
    ||d1||p,||d2||p <= A d0 sigma_p,
    ||e12||p <= A^2 e0 sigma_p.

The same d1 bound holds at every actual positive single-rule node. These bounds do not use a false joint law across quadrature nodes.

### 2.1 Actual shifted terminal expansion

With J=Dg(S), H=D²g(S),

    g(S+d+e)=g(S)+Jd+Je+(1/2)H[d,d]+R_T,

and

    ||R_T||2
       <= B A^4[d0 e0+(A/2)e0^2] sigma4^2
          +(C/6)A^4 u_A^3 sigma6^3.

Indeed the omitted second-order terms are H[d,e]+H[e,e]/2, and the third-order remainder is bounded pointwise by C A |d+e|^3/6. This is an expansion around S, not around x or an unshifted covariance proxy.

### 2.2 The nested one-sample substitution

Let ehat be the pathwise Exp1 x Exp1 average of e12. For each fixed complete history and each a, set Y_h=g(X_(a+h)) and F=F1,a. Taylor-expand g(X_a-Y_h) around X_a-F. The first-order term vanishes after integration over h because int exp(-h)(Y_h-F)dh=0. Therefore

    |int exp(-h)g(X_a-Y_h)dh-g(X_a-F)|
      <= (B A/2) int exp(-h)|Y_h-F|^2 dh
      = (B A/2)[int exp(-h)|Y_h|^2dh-|F|^2].

Dropping the last nonnegative term and applying Minkowski gives

    ||ehat-e||2 <= (B/2) A^3 sigma4^2,
    ||J(ehat-e)||2 <= (B/2) A^4 sigma4^2.

This variance identity avoids the looser factor four from bounding Y_h-F separately.

### 2.3 First centered stencil

For every vector a,

    C0(S,a)=Ja+R_C(S,a),
    |R_C(S,a)| <= (C A/6)|a|^3.

The quadratic terms at S cancel exactly. Consequently the single stencil remainder is at most (C/6)A^4 d0^3 sigma6^3. The nested e12 stencil remainder is at most (C/6)A^7 e0^3 sigma6^3.

For a finite positive single rule, apply the remainder bound directly at its actual nodes. The conditional linear difference is J times the sealed single Gaussian-operator quadrature error applied to g, which contributes delta1 A^2 sqrt(D). Thus no nonlinear operator-norm claim about the complete C0 stencil is being made.

### 2.4 Polarized stencil and the correct ordered-clock identity

The exact integral identity is

    Q0(S,a,b)=(1/8) int_[-1,1]^2
                    D²g(S+s a+t b)[a,b] ds dt.

It yields

    Q0(S,a,b)=(1/2)H[a,b]+R_Q(S,a,b),
    |R_Q(S,a,b)|
       <= (C A/4)|a||b|(|a|+|b|).

Hence ||R_Q(S,d1,d2)||2 <= (C/2)A^4 d0^3 sigma6^3.

For fixed complete history, H[d_t,d_u] is symmetric under t,u exchange. Splitting the iid Exp1 square into t<u and u<t, and setting a=min(t,u), h=|t-u|, gives density 2 exp(-2a-h). Therefore

    int int exp(-t-u) H[d_t,d_u] dt du
      = E_(a~Exp2,h~Exp1) H[d_a,d_(a+h)]
      = H[d,d].

The same ordered identity holds for Q0 because Q0 is symmetric in its last two arguments. Only the true ordered pair is needed; no third history time is required for this redesigned pair of separately weighted stencils.

Summing 2.1--2.4, using conditional Jensen, then sigma4^2<=sqrt(b+2)sqrt(d) and sigma6^3<=sqrt((b+2)(b+4))sqrt(d), and summing squared block norms proves K above. Neither Hessian commutativity nor A sqrt(D) smallness is used.

## 3. Actual finite source and its VALUE bill

Let the three positive, normalized rule sizes be J1, JL, JQ. The linear pair rule approximates Exp1(a) Exp1(h); the quadratic pair rule approximates Exp2(a) Exp1(h). Let U,V at every pair node have their genuine conditional joint law given x,N,M, realized by two standard D-roots L1,L2. Each U and V is a scalar row in (x,N,M,L1,L2), of unit Euclidean row norm. The single bridge may use L1.

All finite nodes share the same roots. With the pinned outer rule and x_i=t_i Z+c_i G, all outer nodes also share the same G,N,M,L1,L2. Define F_node by the displayed finite sums and E_node=F_node-g(x). Then:

    F: (4+3J1+5JL+6JQ) N_out original VALUES,
    E: (5+3J1+5JL+6JQ) N_out original VALUES,
    private root dimension: 5D.

The four shared base sites are g(v),g(w),g(g(w)),g(S). A single node adds g(U),g(S+d),g(S-d). A linear pair node adds g(U),g(V),g(U-g(V)),g(S+e12),g(S-e12). A quadratic pair node adds g(U),g(V) and its four terminal sites. If the pair rules use exactly the same node pairs, g(U),g(V) may be shared and the pair bill becomes 9 J_pair. Do not assume that sharing if their actual supports differ.

Every first or adjoint sweep requires one original HVP at every recorded original VALUE site, with complete replay when the primal was discarded. All captured origins and native compiler repetitions retain their costs.

## 4. Exact derivative products and conservative ports

This section only needs symmetric 0<=Dg<=A I, not B,C or the block assumption. Write

    V0=Dg(v), W=Dg(w), J0=Dg(g(w)),
    P=I-V0/2+J0 W/4,
    Qn=-V0/2+J0 W/2, Rm=J0 W/4.

Then (S_x,S_N,S_M,S_L1,S_L2)=(P,Qn,Rm,0,0). For a descendant U with scalar row coefficients u_y,

    d(U)_y = V0 v_y - Dg(U) u_y,

where v_x=v_N=I/2 and all other v_y=0. For e(U,V)=g(U)-g(U-g(V))-g(g(w)),

    e_y = [Dg(U)-Dg(U-g(V))] u_y
          +Dg(U-g(V))Dg(V) vpair_y -J0 W w_y,

where (w_x,w_N,w_M,w_L1,w_L2)=(I/4,I/2,I/4,0,0). None of these generally noncommuting factors has been reordered.

A first centered stencil has derivative

    D_y C0(S,a)=T_C S_y+K_C a_y,
    T_C=[Dg(S+a)-Dg(S-a)]/2,
    K_C=[Dg(S+a)+Dg(S-a)]/2,

with symmetric ||T_C||<=A/2 and ||K_C||<=A. For H_st=Dg(S+s a+t b), s,t in {+1,-1}, the polarized stencil has

    D_y Q0(S,a,b)=T_Q S_y+U_Q a_y+V_Q b_y,
    T_Q=sum_(s,t) s t H_st/8,
    U_Q=sum_(s,t) t H_st/8,
    V_Q=sum_(s,t) s H_st/8.

Each of these three symmetric matrices has norm <=A/4, using the PSD range [0,A I], rather than merely a four-term triangle inequality.

The aggregate coefficient Hbar of S_y is Dg(S) plus both averaged T_C terms and the averaged T_Q term. It is symmetric and satisfies

    ||Hbar||<=9A/4, ||Hbar-Dg(x)||<=9A/4.

Put

    eta=A/2+A^2/4, q=A/2+A^2/2, r0=A^2/4,
    Cx=9(1+eta)/4+A(13/4+5A/4),
    Cn=9q/4+A(13/4+3A/2),
    Cm=9r0/4+A(5/2+5A/4),
    Cl=A(5/2+A).

Every scalar row coefficient has modulus <=1, so these are conservative simultaneous bounds even without exploiting the first pair row's absent L2 coefficient. With beta=sum omega_i c_i<=sqrt(3)/2,

    raw private first
       <= A sqrt(beta^2 Cx^2+Cn^2+Cm^2+2Cl^2),
    raw caller-Z first <= A Cx/2.

For the square lift E -> (E,0,0,0,0), the G-block's leading Hbar-Dg(x) term is symmetric. Therefore

    raw curl <= 2 beta[(9/4)A eta+A^2(13/4+5A/4)]
                 +A sqrt(Cn^2+Cm^2+2Cl^2).

For A<=1/2 these imply, after half-variance normalization,

    private first < 9A,
    curl < 25 A^2 < (9A)(3A).

All coefficients after division by the displayed A or A² powers are increasing in A, so endpoint evaluation suffices. At A=1/2, the squared normalized first divided by A² is 547655/8192<81. The normalized curl divided by A² is at most (169 sqrt(6)+sqrt(126874))/32<25 (use sqrt(6)<49/20 and sqrt(126874)<357). Safe compiler declarations are ell_E=9A, a_seed=3A, mu=A. The ell_E<=1/4 guard follows from A<=1/36; all other native guards still need their actual checks at dimension 5D. This does not transfer any old sharper normalization constant.

## 5. Energy, caller origins, and literal zeros

Pointwise,

    |g(S)-g(x)| <= A^2|v|+A^3|w|,
    |C0(S,d(U))| <= A^2(|v|+|U|),
    |C0(S,e(U,V))| <= A^3(|V|+|w|),
    |Q0(S,d(U),d(V))|
       <= (A^2/4)(2|v|+|U|+|V|).

The last bound follows from |Q0(S,a,b)|<=A min(|a|,|b|)/2 <= A(|a|+|b|)/4. In particular the residual has linear Gaussian growth O(A²), even without curvature bounds.

Assume the actual numerical rules also have these exact exponential moments:

    single E exp(-s)=1/2,
    linear pair E exp(-(a+h))=1/4,
    quadratic pair E[exp(-a)+exp(-(a+h))]=1.

A product rule with the appropriate separately repaired coordinate moments supplies these identities. For an arbitrary pair rule, separate marginal moments alone do not imply the product moment; it must be checked explicitly or priced as an error.

For every fixed finite p>=2 with ||N(0,I_D)||p<=kappa_p sqrt(D),

    ||E||_(Lp|Z=z)
      <= A^2[(1+3A/8)|z|
          +kappa_p(3/2+5/(2sqrt(2))+A(1+2k))sqrt(D)].

At the actual captured caller origin G=N=M=L1=L2=0, every descendant is evaluated normally with x_i=t_i z. The origin costs the full residual occurrence bill before exact-key reuse and obeys

    |E_origin(z)|<=A^2(1+3A/8)|z|.

Without the additional moments a safe substitute is

    ||E||_(Lp|Z=z)
      <= A^2[(11/8+3A/4)|z|
          +kappa_p(3/2+5/(2sqrt(2))+A(1+2k))sqrt(D)],

since all exponentials are in [0,1]. The same coefficient bounds the origin. Thus the moment identities improve constants; the O(A²) source energy and O(A^4) completion grade do not depend on them.

Subtracting and restoring the actual caller origin leaves private first/curl unchanged and doubles the raw caller first bound to A Cx; its normalized bound is sqrt(2)A Cx. The gradient baseline has its existing normalized caller bound sqrt(2)A. At all caller and root coordinates zero every descendant is literally zero; at nonzero z the caller origin is not zero. At g=0 the entire source is identically zero for every root.

## 6. Source growth for the external pair quadrature theorem

These inequalities identify the actual integrands to which the external positive pair-rule theorem must be applied. They do not assert that the sealed one-time operator rule already controls them.

For arbitrary deterministic y=(x,n,m), U,V, define

    R5=(|x|^2+|n|^2+|m|^2+|U|^2+|V|^2)^(1/2).

The two pair integrands have global linear growth

    |C0(S,e(U,V))| <= A^3 sqrt(11/8) R5,
    |Q0(S,d(U),d(V))| <= (A^2/2) R5.

Indeed ||w||<=sqrt(3/8)||y|| and ||v||<=||y||/sqrt(2). These are source-uniform polynomial-growth bounds, valid without bounded blocks. The additional bound

    |Q0(S,d(U),d(V))| <= (B A/2)|d(U)||d(V)|

is available blockwise, but is not needed for logarithmic-node target accuracy if the linear-growth pair theorem is available. A certified pair error proportional to delta_L times the first growth size and delta_Q times the second yields O(A^4 sqrt(D)) with delta_L=O(A), delta_Q=O(A²), with all theorem-specific dimension-independent constants shown. Taking both tolerances O(A²) is simpler and still polylogarithmic. The sealed single rule needs delta1=O(A²); the outer rule needs delta_out=O(A³).

The full source has integrated norm <=A C_F2 sqrt(D), where one convenient moment-independent bound is

    C_F2=1+A[3/2+5/(2sqrt(2))+A(1+2k)].

This follows directly under the stationary outer coupling from the pointwise residual bound. The outer L² operator rule then contributes delta_out A C_F2 sqrt(D).

## 7. Completion scope

Given the externally certified pair rules and their stated errors, the original guarded positive gradient/near-gradient own-mean compilers can be re-entered using private dimension 5D and the actual normalized ports above. Source energy is O(A²)(|z|+sqrt(D)); the imported near-gradient bracket

    ell_E(a_seed+mu)+ell_E^3(1+mu^(-1/2))

is O(A²), so its allowance remains O(A^4 sqrt(D)) plus restored absolute floors. The gradient branch must separately run at order at least four. Independent complete banks, actual caller-origin restoration, full replay charges, positivity rows, coisometry adapter, and every native caller/radius/curl/dimension/clock/filter/mode guard remain mandatory. This derivation is not an independent audit of those imported compilers or of the new pair-quadrature lemma.

## 8. Quadratic exactness check

For g(x)=Kx with 0<=K<=A I, Q0 vanishes, C0(S,d1)=K²(v-U), and C0(S,e12)=K³(V-w). The stage-two finite source reduces algebraically to

    F_node=Kx-K² sum_single p U +K³ sum_linear-pair q V.

If the actual single and linear-pair rules have E exp(-a)=1/2 and E exp(-(a+h))=1/4, its inner mean is (K-K²/2+K³/4)x. The pinned outer first exponential moment gives the exact canonical m3. This is an own-mean identity; shared-root raw covariance is not identified with a true multi-time covariance. No eigenbasis is read.

## 9. Absolute VALUE, row, and weight floors

Let nu_v, nu_w, nu_b, nu_S be the original VALUE errors at g(v), g(w), g(g(w)), g(S). Use nu_U, nu_V, nu_inner for the actual original VALUE errors at g(U),g(V),g(U-g(V)); terminal errors are those at each actual S plus/minus descendant site. Direct Lipschitz propagation gives the following absolute F floor:

    nu_S
    +sum_single p(nu_plus+nu_minus)/2
    +sum_linear-pair q(nu_plus+nu_minus)/2
    +sum_quadratic-pair r(sum_four_terminal nu)/8
    +(11/2)A nu_v+(9/2)A nu_b+(9/2)A² nu_w
    +A sum_single p nu_U
    +A sum_linear-pair q(nu_U+nu_inner+A nu_V)
    +(A/2)sum_quadratic-pair r(nu_U+nu_V).

For uniform original VALUE error nu this is

    [7/2+14A+(11/2)A²] nu.

E adds the baseline VALUE error; origin subtraction doubles the resulting worst-case floor. These bounds require no division by source energy or by A. They are deliberately conservative; coherent exact-key sharing may reduce the executed bill but must not erase the retained errors.

For an absolute Euclidean coefficient error epsilon in the scalar row of U, the single-stencil integrated floor is at most A² epsilon sqrt(D). For pair rows with coefficient errors epsilon_U,epsilon_V, the corresponding integrated floors are

    linear-pair stencil: A²(epsilon_U+A epsilon_V) sqrt(D),
    quadratic-pair stencil: (A²/4)(epsilon_U+epsilon_V) sqrt(D).

These follow from the exact source derivative products, without division by a small conditional variance or a covariance eigenvalue. They compare each pair row with its intended exact row using the same Gaussian roots. Actual numerical rows and coefficients must still be checked against the row/first/curl declarations used by the compiler.

Absolute probability-weight changes contribute respectively

    single: A² d0 sqrt(D) ||Delta p||1,
    linear pair: A³ e0 sqrt(D) ||Delta q||1,
    quadratic pair: (A² d0/2) sqrt(D) ||Delta r||1.

The actual numerical rules must have nonnegative weights and the declared mass. Moment errors must be priced rather than silently treated as identities. Arithmetic, roots, outer nodes, coefficients, compiler floors, finite-mode errors, and complete replay retain their separate absolute allowances.

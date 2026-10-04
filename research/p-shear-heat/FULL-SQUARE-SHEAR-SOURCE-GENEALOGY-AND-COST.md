# Full-square shear source, exact genealogy, and literal work

New supplement to the held nonlinear-shear candidate
`81c01db6986b0f86fe7315b07928fbb7a3ad3d8c42487a3d371141bf24de34fd`.
This specifies the full-square source and exact fork constants, and certifies
the marked/exposed-ancestor independence from the executed callback graph.
No primitive ancestor query is independently refreshed within an occurrence.

## 1. Why the exposed scalar is not an input of the marked callback

At one positive whole-input clock let c^2+s^2=1. The true orientation integrand
is, by the diagonal OU identity,

    -Sym E_X [(E_U Df(cX+sU))(E_V Curl f(cX+sV))].       (1)

Here U,V are independent complete standard vectors conditional on X. This is
an exact representation of a product of two conditional means, not a claimed
independence between two factors at one realized f query.

For the shear source f(W)=r[g(W+a u g0(v*W))-g(W)], expose ONLY the curl-side
innovation V as

    V=V_perp+Z_v v,   z=c v*X+s Z_v.

The marked-side conditional Jacobian is still

    J_f(X)=E_U Df(cX+sU).

Every original alias inside the marked f(cX+sU) stays on cX+sU. In particular
its ancestor is g0(c v*X+s v*U), not g0(z). The scalar Z_v is absent from this
marked source's read set. After U is averaged, J_f is a deterministic function
of X and therefore independent of Z_v at fixed X, exactly as required in the
weak scalar-heat proof.

The finite implementation has the same ownership:

1. Sample X and the new scalar Z_v as coefficient records; compute z.
2. Define the marked VALUE callback f_X(w)=[f(cX+s w)-f(cX)]/s. Its code reads
   X, w, and its actual old caller only. It does not read Z_v or z.
3. Run the admitted marked pair Q_E(p) on f_X with its own COMPLETE private
   records. Every scalar inside a marked f evaluation comes from that same
   evaluation's w. No ancestor is independently resampled inside it.
4. Define the unmarked side callback below using the captured(X,z). Its two
   complete pair banks are private-independent of the marked bank conditional
   on(X,z,p). They deliberately share the same fresh passive p.

The marked pair's selected coefficient is c_sel J_f(X). This statement is
calibrated from its finite conditional mean; no sampled Jacobian is read or
frozen. All marker numerical/source errors retain their actual conditional
(X,p) envelopes and are averaged over the additional unused Z_v afterward.
Thus independence is certified by the input graph and the exact identity(1).
It is not added as a free assumption to a shared-record derivative product.

## 2. Full n-square genuine-gradient side source

Let Pi=I-vv*, eta=b s^2, and q0=c Pi X+z v+a u g0(z). Define, on R^n,

    h_(X,z)(x)
      =(rho/s) Pi[g(q0+s Pi x)-g(q0)]
        +(rho/eta) v[g0(z+eta v*x)-g0(z)].               (2)

This is a genuine full gradient: its two terms are gradients on the orthogonal
v-perp and v subspaces. Its active first has norm at most rho up to declared
fixed constants. Its selected mean matrix is

    E Dh =rho D,  D=Hbar_X(z)+m_eta(z) vv*,
    Hbar v=0,  m_eta(z)=E g0'(z+eta Z).

No perpendicular basis is needed. All projections use the fixed known u,v.
Put

    J=u v*-v u*,  Pi0=I-uu*-vv*.

Then ||J||=1, J J*=uu*+vv*, and Pi0 is the exact known square root of
I-JJ*. Let P1,P2 be independent COMPLETE genuine-gradient pair banks on the
SAME source(2), with the same captured X,z and finite versions. Execute

    Q_C(p)=P2(J P1(p)+Pi0 Z_mid).

Both pairs and Z_mid are private-independent conditional on(X,z,p). The
middle fill gives standard input to the outer pair after the first pair's
reference input is integrated. The exact selected coefficient is

    K=c_sel^2 rho^2 D J D
     =c_sel^2 rho^2 m_eta(z)
          [Hbar_X(z)u v*-v u*Hbar_X(z)].                (3)

It is skew. Conditional on p, its reference covariance is I-KK*. There are
no input/output endpoint pads and no factor1/2 in(3).

For M=c_sel J_f(X), tau=kappa0=ra and
O_(X,z),eta=-tau Sym[J_f(X)m_eta(z)(Hbar u v*-v u*Hbar)],

    M K*+K M*=2 c_sel^3 rho^2 O_(X,z),eta/tau.

Hence the positive fork at fixed numerical variance v0 uses

    b Q_E+c_side Q_C+sqrt(v0-b^2-c_side^2)Z_keep,
    b c_side=-q tau/(2 c_sel^3 rho^2).                 (4)

Known signed tau and q determine the actual readout. If tau=0, use the direct
Gaussian branch; never divide by tau. All records are attached before the
same-endpoint comparison. The positive reference is conditional on the same
X,z, and the new p is independent of the old core when invoking its marginal.

## 3. Complete paths and numerical widths

The actual old-label source firsts of(2) are O(rho/s) from terminal centers and
O(rho/eta) from scalar centers, including BOTH anchors. The path through the
stored g0(z) inside q0 has the actual factor a g0'(z), bounded by a; it is not
an uncharged constant. The two-pair chain gives the conservative whole-root
row O(rho/eta), because an outer source-root path need not gain inner rho.
After(4) its side row is O(kappa0/(b rho eta)).

Multiply each COMPLETE fork by its actual sqrt(w_j). With eta_j=b s_j^2 and
orientation-only s_min comparable to b, the full absolute source/caller sum is
Lambda kappa0/(b^3 rho), as in the candidate. The copied-root first has its
corresponding square sum. This is a conservative complete-path bound. No
sharper block-projection cancellation is needed for the stated candidate.

eta is an EXECUTED scalar Gaussian width. Making it small does not cause an
inverse-eta sample multiplicity in a fixed-order gradient pair, but it does:
- enlarge the actual captured-root/caller row above;
- set a real minimum query width and inverse readout;
- tighten finite original-source/row/filter VALUE tolerances;
- add its log(1/eta), known query sizes, and all path amplification to the
  inherited fixed-order serial cost and small-radius guard.
All those are declared before the matched pair banks are built. A bound on
the averaged matrix does not replace any of these executed costs.

## 4. Original-gradient calls, cache, and FIRST sweeps

At fixed X,z compute g0(z), q0 and g(q0) once and keep their actual VALUE and
caller records. A source call(2) at a new x needs one terminal g query at
q0+s Pi x and one scalar g0 query at z+eta v*x; the anchors can be reused only
because their exact inputs, finite versions and captured records agree.
If g0 is the single-potential chart v*[g(x0+t v)-g(x0)]-lambda t, each new
scalar query is one original g call with a known readout and linear term;
g(x0) is its separately stored common anchor. Every such original call is
counted. If caches are discarded, replay their true costs.

Let Q_P,B be the fully expanded call count of one genuine-gradient pair on
its complete source, and Q_f the actual complete whole-f VALUE count. One
clock costs a complete marked-pair program on f_X, plus

    2 Q_P,B Q_h + Q_anchor + Q_known,

where Q_h includes both differences in(2), and Q_known includes the middle
projector, vector arithmetic and Gaussian records. This does not evaluate an
unknown Hbar or m_eta. Sum over the positive clock list; no N^2 or
inverse-eta replication is introduced. There is also no uniform polynomial
bound in a growing fixed-order index B.

One requested forward/adjoint sweep differentiates this literal VALUE graph
using original HVPs at saved queries. It accumulates all contributions into
the shared g0(z) and g(q0) anchors before their own sweeps. Derivatives of HVP
outputs are not used. The complete actual source version and all aliases are
held across law, current and retention comparisons.

## 5. Remaining join boundary

This supplement certifies the source typing/genealogy and exact two-pair
normalizer. The scalar weak coefficient proof remains candidate Sections2-3.
The owned-root Price/independent constant-reference and retained-host returns
must be supplied at the stated complete root frames. Reusing an old mean law
that already integrated a still-observed public is not permitted. An
orientation-only clock differs from a regenerated forward/Cov clock; its
legality rests on that independent negative-kernel consumer.

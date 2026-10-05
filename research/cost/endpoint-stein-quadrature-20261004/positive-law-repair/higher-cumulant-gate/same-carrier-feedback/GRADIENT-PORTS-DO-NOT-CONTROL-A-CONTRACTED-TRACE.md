# Gradient, first, and one-energy ports do not control a contracted trace

2026-10-05. Port-insufficiency counterexample; independent review requested. This is NOT a counterexample to the actual same-g m3 covariance expansion. The field and the preserved Hessian coefficient below are chosen independently. Their missing common genealogy is exactly why the example cannot settle that expansion.

## Result

There are smooth, anchored, genuine-gradient vector fields E_D with

    E_D(0)=0, curl(E_D)=0,
    Lip(E_D)<=C A,
    ||E_D||_(L2(gamma_D))<=C A2 sqrt(D), A=D^(−1/2),

and a constant preserved coefficient W_D=A2 I, for which

    ||R1[W_D:D2 E_D]||_(L2(gamma_D))>=c A4 D.        (1)

Here [W:D2E]_i=sum_(a,b)W_ab partial_a partial_b E_i. The proposed O(A4 sqrt(D)) bound therefore fails by a sqrt(D)=1/A factor, including any fixed public-log loss. Even genuine-gradient output and exact zero do not remove this trace amplification.

The matrix W is itself the product of two constant symmetric Hessians A I. However E is NOT the chord/current produced by that same quadratic gradient. Any successful m3 bound must use additional literal same-source genealogy or a stronger current identity; the listed ports alone cannot imply it.

## 1. A globally smooth anchored construction

Choose fixed smooth cutoff functions on [0,infinity):

    chi_out(s)=1 for s<=3/2, chi_out(s)=0 for s>=2,
    chi_in(s)=1 for s<=1/4, chi_in(s)=0 for s>=1/2.

They may take values in [0,1], with all needed derivatives bounded independently of D. Define

    phi_D(x)=A2 x1 [(|x|2−D−2) chi_out(|x|2/D)
                                  +(D+2) chi_in(|x|2/D)],
    E_D=grad phi_D.                                  (2)

The inner cutoff cancels the gradient at the origin exactly. The outer cutoff makes the field globally controlled. E is C-infinity, a genuine gradient, and E(0)=0.

Set z=x/sqrt(D). Then

    phi_D(x)=A2 D^(3/2) f_D(z),
    f_D(z)=z1[(|z|2−1−2/D)chi_out(|z|2)
                           +(1+2/D)chi_in(|z|2)].

The Hessians D2 f_D are uniformly bounded for D>=1. Hence

    ||D E_D||op<=C A2 sqrt(D)=C A.                   (3)

No dimension/heat restriction beyond the specified counterfamily is hidden in (3).

## 2. Gaussian energy and the exact trace

On the annulus D/2<=|x|2<=3D/2, the cutoffs give the cubic Hermite potential

    phi_0(x)=A2 x1(|x|2−D−2).

Its gradient is

    (E_0)_1=A2[3(x1²−1)+sum_(j>1)(xj²−1)],
    (E_0)_i=2 A2 x1 xi, i>1.

Gaussian orthogonality gives exactly

    E E_0=0,
    ||E_0||2²=A4(6D+12),
    Delta E_0=A2(2D+4)e1.                           (4)

The Gaussian probability of leaving the annulus is exponentially small in D. E_D and E_0 have polynomial growth of size at most C A2(D+|x|2), so standard chi-square Gaussian tails imply

    ||E_D−E_0||2<=C A2 D^c exp(−c' D)               (5)

for fixed constants c,c'>0. In particular ||E_D||2<=C A2 sqrt(D) for large D (and constants can absorb the finitely many smaller D).

Gaussian integration by parts twice, componentwise, yields

    E Delta E_D=E[(|X|2−D)E_D(X)].

Use (5) and finite Gaussian fourth moments to compare this with (4). The result is

    |E Delta E_D−A2(2D+4)e1|
                        <=C A2 D^c exp(−c' D).     (6)

Thus |E Delta E_D|>=c A2 D for all large D. Since R1 preserves Gaussian expectations and W=A2 I,

    ||R1[W:D2E_D]||2
      >=|E R1[A2 Delta E_D]|
       =A2 |E Delta E_D|>=c A4 D.                  (7)

At A=D^(−1/2), the left lower scale is A2 while A4 sqrt(D)=A3. This proves (1).

## 3. Exact scope for the m3 investigation

The example falsifies a proposed PORT-ONLY implication using:

- a preserved Hessian product with operator size A2;
- a correction with actual first O(A), energy O(A2 sqrt(D)), exact zero and even zero curl;
- two Gaussian scores or an unqualified Hodge contraction.

It does not use a joint original-g construction. For the true m3 source, E and W have shared ancestors and constrained algebraic relationships. The same-g conditional Gaussian interpolation current must remain intact; this note gives no license to declare it false or to replace it by an unrelated field.

The role of this separator is to prevent an invalid generic lemma from closing that gate. A valid proof needs an additional trace cancellation, a stronger typed source/current condition, or an exact same-g regrouping that the present independently chosen fields do not satisfy.

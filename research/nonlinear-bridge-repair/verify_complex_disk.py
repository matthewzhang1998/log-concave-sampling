"""Exact rational audit of the complex-disk certificate; requires SymPy."""
import sympy as S
s = S.symbols('s')
a,b = S.Rational(65,64), S.Rational(63,64)
E = sum((2*b*s)**j / S.factorial(j) for j in range(9))
Pm = 1+4*a*a*s*s*(1+(1-b*s)**2)
Pp = 1+4*a*a*s*s*(1+(a*s-1)**2)
R = S.Poly(S.cancel((E-Pm)/s),s)
B = [sum(R.nth(j)*S.binomial(k,j)/S.binomial(7,j) for j in range(k+1)) for k in range(8)]
Bden = 24629060462182400
Bnum = [48488462784921600,26273173943091200,15076247889510400,12524793794396160,16419361684717568,24828792912732160,36265100707349760,50076245130694685]
assert B == [S.Rational(n,Bden) for n in Bnum]
assert all(x>S.Rational(1,2) for x in B)
Cp = S.Poly((E-Pp).subs(s,s+1),s)
C = [Cp.nth(k) for k in range(9)]
Cden = 703687441776640
Cnum = [1430749860876991,4011516377289976,3642743757137380,228474934592712,2444527082810,1071611807151816,279243854976996,47517861678840,3938980639167]
assert C == [S.Rational(n,Cden) for n in Cnum]
assert all(x>0 for x in C)
assert S.expand(sum(B[k]*S.binomial(7,k)*s**k*(1-s)**(7-k) for k in range(8))-R.as_expr()) == 0
print('PASS: exact Bernstein and shifted-power certificates; relative disk radius 1/64.')

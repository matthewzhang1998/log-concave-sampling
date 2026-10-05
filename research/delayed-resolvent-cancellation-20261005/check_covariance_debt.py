import json, math
from fractions import Fraction
import numpy as np
from numpy.polynomial.legendre import leggauss
eta=.5

def M(a):
    return .5*((1+a)*np.exp(-.5*(1+a)**2)+(1-a)*np.exp(-.5*(1-a)**2))

def integrand(r,s):
    vr=1-r*r; vs=1-s*s
    dC=np.sqrt(vr*vs)
    dT=np.minimum(r,s)/np.maximum(r,s)-r*s
    diff=dC-dT
    bracket=diff*(math.exp(-.5)+eta*(np.exp(-vs/2)*M(s)+np.exp(-vr/2)*M(r)))
    bracket+=eta*eta/2*np.exp(-(vr+vs)/2)*((np.exp(dC)-np.exp(dT))*M(r-s)-(np.exp(-dC)-np.exp(-dT))*M(r+s))
    return eta/(4*(1+eta)**3)*bracket,diff

rows=[]
for n in [4,8,16,32,64,128,256]:
    r,p=leggauss(n); r=(r+1)/2;p=p/2
    c,I=integrand(r[:,None],r[None,:])
    cn=float(p@c@p); In=float(p@I@p)
    rows.append(dict(kind='Gauss-Legendre square',nodes=n,c_eta=cn,covariance_mass=In,first_moment=float(p@r),A001_remainder_adjusted_coefficient=cn-.01/6))
# Smooth-triangle integration: r=s*u and symmetry, with Jacobian 2s.
for n in [32,64,128,256]:
    z,p=leggauss(n); z=(z+1)/2;p=p/2
    s=z[:,None];u=z[None,:]; r=s*u
    c,I=integrand(r,s)
    rows.append(dict(kind='Gauss-Legendre triangle',nodes_each=n,c_eta=float(p@(2*s*c)@p),covariance_mass=float(p@(2*s*I)@p)))
# Four actual finite positive clock nodes, no limiting integral needed.
r=np.array([1,3,5,7],float)/8;p=np.ones(4)/4
c,I=integrand(r[:,None],r[None,:])
finite=dict(nodes=list(r),weights=list(p),c_eta=float(p@c@p),covariance_mass=float(p@I@p),rational_covariance_mass_lower=str(Fraction(253,320)**2-Fraction(127,420)),rational_excess_over_5_16=str(Fraction(253,320)**2-Fraction(127,420)-Fraction(5,16)))
# The exact finite proof uses these four squared rational root lower bounds.
assert Fraction(79,10)**2<63 and Fraction(74,10)**2<55 and Fraction(62,10)**2<39 and Fraction(38,10)**2<15
assert Fraction(253,320)**2-Fraction(127,420)>Fraction(5,16)
assert Fraction(7,1728)-Fraction(1,600)>Fraction(1,500)
checks=dict(eta=eta,coefficient_lower_bound=str(Fraction(7,1728)),continuous_covariance_mass=math.pi**2/16-.25,continuous_covariance_lower_coefficient=(math.exp(-.5)-.25)*(math.pi**2/16-.25)/27,finite_four_node=finite,quadrature_diagnostics=rows,quadratic_defects=[dict(A=A,coefficient=-3/8+math.sqrt(A)*(1-math.sqrt(A))**2/4) for A in [.01,.001,.0001]])
print(json.dumps(checks,indent=2))

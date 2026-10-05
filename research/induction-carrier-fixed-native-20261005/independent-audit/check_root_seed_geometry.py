#!/usr/bin/env python3
"""Independent finite-matrix checks of the pinned root-seed port.

These are algebra checks, not substitutes for the analytic proof in the audit.
No primary source files are edited.
"""
import hashlib, json
from pathlib import Path
from fractions import Fraction
import numpy as np

OUT = Path(__file__).resolve().parent
rng = np.random.default_rng(531805)
D, J, m = 6, 4, 6
v = .7
# P_h consists of three root primitive Gaussians, two marked publics,
# and one numerical readout buffer. The root innovation is a normalized
# affine row of the first three coordinates.
a = np.array([2., -1., 3.]); a /= np.linalg.norm(a)
c = .31
b = np.array([.22, -.27])
g = np.sqrt(v-c*c-b@b)
ell = np.r_[c*a, b, g]
L = np.kron(ell[None, :], np.eye(D))
G = np.kron(np.r_[a, np.zeros(3)][None, :], np.eye(D))
Q = np.kron(np.array([[0,0,0,1,0,0],[0,0,0,0,1,0]]), np.eye(D))
Pi = np.eye(m*D)-L.T@L/v
N = D+J*m*D+J*D
H = np.zeros((D,N)); H[:,:D] = np.sqrt(v)*np.eye(D)
max_cov = max_carrier = max_gq = max_root_h = 0.
max_endpoint = 0.
max_covroot_ratio = 0.
for h in range(J):
    A = np.zeros((m*D,N))
    A[:,:D] = L.T/np.sqrt(v)
    off = D+h*m*D
    A[:,off:off+m*D] = Pi
    max_cov = max(max_cov, float(np.max(np.abs(A@A.T-np.eye(m*D)))))
    max_carrier = max(max_carrier, float(np.max(np.abs(L@A-H))))
    Gr, Qr = G@A, Q@A
    max_gq = max(max_gq, float(np.max(np.abs(Gr@Qr.T))))
    max_root_h = max(max_root_h, float(np.max(np.abs(Gr@H.T-c*np.eye(D)))))
    M = rng.normal(size=(D,D)); M *= .18/np.linalg.norm(M,2)
    ev,U = np.linalg.eigh(np.eye(D)-M@M.T)
    C = (U*np.sqrt(ev))@U.T
    K = np.zeros((D,2*D)); K[:,:D] = .21*np.eye(D)
    child = K@Qr
    child[:, D+J*m*D+h*D:D+J*m*D+(h+1)*D] += np.sqrt(1-.21**2)*np.eye(D)
    # Explicit ideal root and marked affine readout, in the SAME geometry.
    root = M@child+C@Gr
    markbuf = (L-c*G)@A
    direct = c*root+markbuf
    restored = H+c*M@child+c*(C-np.eye(D))@Gr
    max_endpoint = max(max_endpoint, float(np.max(np.abs(direct-restored))))
    max_covroot_ratio = max(max_covroot_ratio, float(np.linalg.norm(C-np.eye(D),'fro')/(np.sqrt(D)*np.linalg.norm(M,2)**2)))

# The one-Hilbert/two-frame inequality used in the current bound.
ratios=[]
for _ in range(50):
    RF = rng.normal(size=(D,9*D))
    DF = rng.normal(size=(D,9*D))
    ratios.append(float(np.linalg.norm(RF@DF.T,'fro')/(np.linalg.norm(RF,'fro')*np.linalg.norm(DF,2))))

# Wrong completion test: in scalar dimension H=c G+b Q (var H=v),
# replacing (C-1)G with an independent G' changes covariance with H.
r=.15; C=np.sqrt(1-r*r)
correct_cross = v+c*c*(C-1)
wrong_cross = v

results={
    'D':D, 'packet_count':J, 'primitive_blocks_per_packet':m, 'variance_v':v,
    'individual_packet_covariance_identity_max_error':max_cov,
    'literal_common_carrier_identity_max_error':max_carrier,
    'root_innovation_marked_public_covariance_max_error':max_gq,
    'root_innovation_H_covariance_minus_cI_max_error':max_root_h,
    'exact_restored_endpoint_identity_max_error':max_endpoint,
    'covariance_root_HS_over_sqrtD_opM_squared_max_ratio':max_covroot_ratio,
    'one_Hilbert_matrix_inequality_max_ratio':max(ratios),
    'wrong_independent_completion_cross_covariance_error':wrong_cross-correct_cross,
    'triple_cubic_root_grade':str(Fraction(17,2)),
    'root_quadratic_grade':str(2*Fraction(17,2)),
    'root_times_child_grade_b24':str(Fraction(17,2)+Fraction(24,16)),
    'root_times_child_grade_b25':str(Fraction(17,2)+Fraction(25,16)),
}
assert max(max_cov,max_carrier,max_gq,max_root_h,max_endpoint)<1e-12
assert max_covroot_ratio<1
assert max(ratios)<=1+1e-12
assert wrong_cross-correct_cross>0
(OUT/'geometry-check-results.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))

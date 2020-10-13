import Solutions as sol
import Potentials as pt
from scipy.linalg import eig
import numpy as np

## define the background matrix K0

def K0(sol_backg,params):
    r0 , l0 = sol_backg[0] , sol_backg[2]
    r0inv = 1/r0
    r0_m2 = r0inv*r0inv
    r0_m3 = r0_m2*r0inv
    omega0_inv = 1/omega(r0)
    omega0_m2 = omega0_inv*omega0_inv
    return np.array([[0.,0.,omega0_inv,0.],[2*l0*r0_m3*omega0_m2,0.,0.,-r0_m2*omega0_m2],[-r0_m3*omega0_inv,0.,0.,-2*l0*r0_m3*omega_0inv],[0.,0.,0.,0.]])

## define the modewise perturbed potential matrix

def K1(sol_backg,params):
    m = params[2]
    psi = pt.modewise_Phi_grav(params)
    dpsi = pt.modewise_dPhi_grav(params)
    ddpsi = pt.modewise_ddPhi_grav(params)
    r0 , l0 = sol_backg[0] , sol_backg[2]
    omega0_inv = 1/sol.omega(r0)
    omega0_m2 = omega0_inv*omega0_inv
    r0_inv = 1/r0
    r0_m2 = r0_inv*r0_inv
    r0_m3 = r0_inv*r0_m2
    K1_pr = omega0_inv*(ddpsi + 2*l0*ro_m3*omega0_inv*dpsi)
    K1_lr = -m*omega0_inv*(dpsi + 2*l0*ro_m3*omega0_inv*psi)
    K1_pl = -dpsi*r0_m2*omega0_m2
    K1_ll = m*psi*r0_m2*omega0_m2
    return np.array([[0.,0.,0.,0.],[0.,0.,0.,0.],[K1_pr,0.,0.,K1_pl],[K1_lr,0.,0.,K1_ll]])

## define the modewise perturbed momentum matrix

def K2(sol_pert,sol_backg,params):
    r0 , l0 = sol_backg[0] , sol_backg[2]
    r0_inv = 1/r_0
    r0_m2 = r0_inv*r0_inv
    r0_m3 = r0_m2*r0_inv
    p1 = sol_pert[1]
    omega0_inv = 1/sol.omega(r0)
    omega0_m2 = omega0_inv*omega0_inv
    term1 = 2*p1*l0*r0_m3*omega_0_m2
    term2 = -p1*r0_m2*omega0_m2
    return np.array([[term1,0.,0.,term2],[0.,0.,0.,0.],[0.,0.,0.,0.],[0.,0.,0.,0.]])

## define the modewise perturbed solution matrix

def K3(sol_pert,sol_backg,params):
    r0 , l0 = sol_backg[0] , sol_backg[2]
    r0_inv = 1/r_0
    r0_m2 = r0_inv*r0_inv
    r0_m3 = r0_m2*r0_inv
    r0_m5 = r0_m2*r0_m3
    r1 , l1 = sol_pert[0],sol_pert[2]
    omega0_inv = 1/sol.omega(r0)
    omega0_m2 = omega0_inv*omega0_inv
    r1_norm = r1*r0_inv
    l1_norm = l1/l0
    omega1_norm = (l1*r0_m2 -2*l0*r1*r0_m3)*omega0_inv\
    K3_rp = -omega1_norm*omega0_inv
    K3_tr = 2*l0*r0_m3*omega0_m2*(l1_norm - 3*r1_norm - 2*omega1_norm)
    K3_tl = r0_m2*omega0_inv*(2*r1_norm + 2*omega1_norm)
    K3_pr = -r0_m3*omega0_inv*(-omega1_norm + 6*l1_norm - 6*r1_norm) + 2*l0*r0_m5*omega_0_m2*(2*l1_norm -r1_norm)
    K3_pl = 2*l0*r0_m3*omega0_inv*(-omega1_norm + l1_norm - 3*r1_norm) - r0_m2*r0_m2*omega0_m2*(2*l1_norm - r1_norm)
    return np.array([[0.,0.,K3_rp,0.],[K3_tr,0.,0.,K3_tl],[K3_pr,0.,0.,K3_pl],[0.,0.,0.,0.]])

## define the averaged stability matrix

def K(params):
    sol_backg = sol.backg_sol(0,params)
    r1 , l1 = sol.r_m(params) , sol.l_m(params)
    sol_pert = np.array([r1,-m*r1,l1])
    mat_0 = K0(sol_backg,params)
    mat_1 = K1(sol_backg,params)
    mat_2 = K2(sol_pert,sol_backg,params)
    mat_3 = K3(sol_pert,sol_backg,params)
    return mat_1 + mat_2 + mat_3 + mat_4
    

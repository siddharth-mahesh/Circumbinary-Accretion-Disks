from numpy import sin , cos
import Potentials as pt

## solution outputted as r , pr , l , t

def omega(r):
    return r**(-1.5) - 1

def omega1(backg_sol,pert_sol,params):
    r0 = params[0]
    r0inv = 1/r0
    r0_m2 = r0inv*r0inv
    r0_m3 = r0inv*r0_m2
    return pert_sol[2]*r0_m2 - 2*backg_sol[2]*pert_sol[1]*r0_m3

## define the background solutions

def backg_sol(phi,params):
    r_0 = params[0]
    omega_0 = omega(r_0)
    return [r_0 , 0 , (omega_0 + 1)**(-1/3) , phi/omega_0]

## define the modewise perturbed solutions

def r_m(params):
    r_0 , q , m = params[0] , params[1] , params[2]
    omega_0 = omega(r_0)
    omega_bar = 1 + 1/omega_0
    psi_m = pt.modewise_Phi_grav(params)
    dpsi_m = pt.modewise_dPhi_grav(params)
    delta_m = m*m - omega_bar*omega_bar
    return (2*omega_bar*psi_m + r_0*dpsi_m)/(r_0*omega_0*omega_0*delta_m)

def l_m(params):
    r_0 , q , m = params[0] , params[1] , params[2]
    omega_0 = omega(r_0)
    psi_m = pt.modewise_Phi_grav(params)
    return -psi_m/omega_0

def modewise_pert_sol(phi,params):
    r_0 , m = params[0] , params[2]
    s_phi = sin(m*phi)
    c_phi = cos(m*phi)
    R_m = r_m(params)
    L_m = l_m(params)
    omega_0 = omega(r_0)
    T_m = 2*(1+1/omega_0)*R_m/r_0 - L_m/(r_0*r_0*omega_0*omega_0)
    return [R_m*c_phi,-m*R_m*omega_0*s_phi,L_m*c_phi,T_m*s_phi/m]

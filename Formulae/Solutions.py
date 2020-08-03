from numpy import sin , cos
import Potentials as pt

## solution outputted as r , pr , l , t

def omega(r):
    return r**(-1.5) - 1

## define the background solutions

def backg_sol(phi,params):
    params[0] = r_0
    omega_0 = omega(r_0)
    return [r_0 , 0 , (omega_0 + 1)**(-1/3) , phi/omega_0]

## define the modewise perturbed solutions

def r_m(params):
    r_0 , q , m = params[0] , params[1] , params[2]
    omega_0 = omega(r_0)
    omega_bar = 1 + 1/omega_0
    psi_m = pt.phi_grav(r_0,q)
    dpsi_m = pt.dphi_grav(r_0,q)
    delta_m = m*m - omega_bar*omega_bar
    return (2*omega_bar*psi_m + r_0*dpsi_m)/(r_0*omega_0*omega_0*delta_m)

def l_m(params):
    r_0 , q , m = params[0] , params[1] , params[2]
    omega_0 = omega(r_0)
    psi_m = pt.phi_grav(r_0,q)
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

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import factorial,lpmv

## Define the background potential

def backg_Phi(params):
    return -1/params[0]

## Define the background l-polar potential

def backg_multipole_Phi(l,params):
    r , q = params[0] , params[1]
    Q_l = ((1-q)*(-q)**l + q*((1-q)**l))
    W_lm = W(l,0)
    return -2*Q_l*W_lm*W_lm/(r**(l+1))

## define the W coefficient from integrating the Legendre Polynomial

def W(l,m):
  return ( ( factorial(l-m)/factorial(l+m) ) **0.5) * lpmv(m,l,0)

## define the modewise perturbed potential coefficient and its radial derivative

def modewise_Phi_grav(params):
  r , q , m = params[0] , params[1] , params[2]
  lmin = max(2,m)
  lmax = lmin + 6
  #lmax = 2
  pmn = 0
  for l in range(lmin,lmax+1):
    Q_l = ((1-q)*(-q)**l + q*((1-q)**l))
    W_lm = W(l,m)
    pmn += -2*Q_l*W_lm*W_lm/(r**(l+1))
  return pmn

def modewise_dPhi_grav(params):
  r , q , m = params[0] , params[1] , params[2]
  lmin = max(2,m)
  lmax = lmin + 6
  #lmax = 2
  dpmn = 0
  for l in range(lmin,lmax+1):
    Q_l = ((1-q)*(-q)**l + q*((1-q)**l))
    W_lm = W(l,m)
    dpmn += -2*(-l-1)*Q_l*W_lm*W_lm/(r**(l+2))
  return dpmn

def modewise_ddPhi_grav(params):
  r , q , m = params[0] , params[1] , params[2]
  lmin = max(2,m)
  lmax = lmin + 6
  #lmax = 2
  ddpmn = 0
  for l in range(lmin,lmax+1):
    Q_l = ((1-q)*(-q)**l + q*((1-q)**l))
    W_lm = W(l,m)
    ddpmn += -2*(-l-2)*(-l-1)*Q_l*W_lm*W_lm/(r**(l+3))
  return ddpmn    

## define the perturbed potential

def pert_Phi(phi,params):
    m = params[2]
    phi_m = modewise_Phi_grav(params)
    return phi_m*np.cos(m*phi)

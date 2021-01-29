import numpy as np
import matplotlib.pyplot as plt
from scipy.special import factorial,lpmv

## Define the Wigner little-d matrix with inputs:
## j - the angular momentum quantum number
## m - the magnetic quantum number in the frame of choice 
## mp - the magnetic quantum number in the rotated frame 
## beta - the inclination between the x-y planes of both frames

def wignerd(j,m,mp,beta):
  prefactor = np.sqrt(factorial(j+mp)*factorial(j-mp)*factorial(j+m)*factorial(j-m))
  temp = 0
  s = 0
  run = 1
  th = 0.5*beta
  while run != 0:
    f1 , f2 , f3 = j + m - s , m - m + s , j - mp - s
    if f1 < 0 or f2 < 0 or f3 < 0:
      run = 0
    else:
      temp += ((-1)**(mp - m + s) * (np.cos(th))**(2*j + m - mp - 2*s) * (np.sin(th))**(mp - m + 2*s))/(factorial(f1)*factorial(f2)*factorial(f3)*factorial(s))
      s += 1
  return prefactor * temp

## define the W coefficient from integrating the Legendre Polynomial

def W(l,m):
  return ( ( factorial(l-m)/factorial(l+m) ) **0.5) * lpmv(m,l,0)


## Define the background potential

def backg_Phi(params):
    return -1/params[0]

## Define the background l-polar potential

def backg_multipole_Phi(l,params):
    r , q = params[0] , params[1]
    Q_l = ((1-q)*(-q)**l + q*((1-q)**l))
    W_lm = W(l,0)
    return -2*Q_l*W_lm*W_lm/(r**(l+1))

## define the modewise perturbed potential coefficient and its radial derivative

def modewise_Phi_grav(params):
  r , q , i , m , n = params[0] , params[1] , params[2] , params[3] , params[4]
  lmin = max(2,m)
  lmax = lmin + 6
  #lmax = 2
  pmn = 0
  for l in range(lmin,lmax+1):
    Q_l = ((1-q)*(-q)**l + q*((1-q)**l))
    wigner = wignerd(l,n,m,i)
    pmn -= 2*Ml*W(l,m)*W(l,n)*wigner/(r**(l+1))
  return pmn

def modewise_dPhi_grav(params):
  r , q , i , m , n = params[0] , params[1] , params[2] , params[3] , params[4]
  lmin = max(2,m)
  lmax = lmin + 6
  #lmax = 2
  dpmn = 0
  for l in range(lmin,lmax+1):
    Q_l = ((1-q)*(-q)**l + q*((1-q)**l))
    wigner = wignerd(l,n,m,i)
    dpmn += 2*(l+1)*Ml*W(l,m)*W(l,n)*wigner/(r**(l+2))
  return dpmn

def modewise_ddPhi_grav(params):
  r , q , i , m , n = params[0] , params[1] , params[2] , params[3] , params[4]
  lmin = max(2,m)
  lmax = lmin + 6
  #lmax = 2
  ddpmn = 0
  for l in range(lmin,lmax+1):
    Q_l = ((1-q)*(-q)**l + q*((1-q)**l))
    wigner = wignerd(l,n,m,i)
    dpmn -= 2*(l+1)*(l+2)*Ml*W(l,m)*W(l,n)*wigner/(r**(l+3))
  return ddpmn    

## define the perturbed potential

def pert_Phi(phi,params):
    m = params[4]
    phi_m = modewise_Phi_grav(params)
    return phi_m*np.cos(m*phi)
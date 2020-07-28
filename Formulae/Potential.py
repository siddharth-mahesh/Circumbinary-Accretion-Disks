import numpy as np
import matplotlib.pyplot as plt
from scipy.special import factorial,lpmv

## define the W coefficient from integrating the Legendre Polynomial

def W(l,m):
  return ( ( factorial(l-m)/factorial(l+m) ) **0.5) * lpmv(m,l,0)

## define the potential and its radial derivative

def Phi_grav(m,r,q=0.5):
  lmin = max(2,m)
  lmax = lmin + 6
  pmn = 0
  for l in range(lmin,lmax+1):
    Q_l = ((1-q)*(-q)**l + q*((1-q)**l))
    W_lm = W(l,m)
    pmn -= 2*Q_l*W_lm*W_lm/(r**(l+1))
  return pmn

def dPhi_grav(m,r,q=0.5):
  lmin = max(2,m)
  lmax = lmin + 6
  dpmn = 0
  for l in range(lmin,lmax+1):
    Q_l = ((1-q)*(-q)**l + q*((1-q)**l))
    W_lm = W(l,m)
    dpmn += 2*(l+1)*Q_l*W_lm*W_lm/(r**(l+2))
  return dpmn

# -*- coding: utf-8 -*-
# http://pymiescatt.readthedocs.io/en/latest/forwardCS.html
import numpy as np
from scipy.special import jv, yv
from PyMieScatt.Mie import MieQ, MiePiTau

def MieQCoreShell(mCore,mShell,wavelength,dCore,dShell):
#  http://pymiescatt.readthedocs.io/en/latest/forwardCS.html#MieQCoreShell
  xCore = np.pi*dCore/wavelength
  xShell = np.pi*dShell/wavelength
  if xCore==xShell:
    return MieQ(mCore,wavelength,dShell)
  elif xCore==0:
    return MieQ(mShell,wavelength,dShell)
  elif mCore==mShell:
    return MieQ(mCore,wavelength,dShell)
  elif xCore>0:
    nmax = np.round(2+xShell+4*(xShell**(1/3)))
    n = np.arange(1,nmax+1)
    n1 = 2*n+1
    n2 = n*(n+2)/(n+1)
    n3 = n1/(n*(n+1))
    xShell2 = xShell**2

    an, bn = CoreShell_ab(mCore,mShell,xCore,xShell)

    qext = (2/xShell2)*np.sum(n1*(an.real+bn.real))
    qsca = (2/xShell2)*np.sum(n1*(an.real**2+an.imag**2+bn.real**2+bn.imag**2))
    qabs = qext-qsca

    g1 = [an.real[1:int(nmax)],an.imag[1:int(nmax)],bn.real[1:int(nmax)],bn.imag[1:int(nmax)]]
    g1 = [np.append(x, 0.0) for x in g1]
    g = (4/(qsca*xShell2))*np.sum((n2*(an.real*g1[0]+an.imag*g1[1]+bn.real*g1[2]+bn.imag*g1[3]))+(n3*(an.real*bn.real+an.imag*bn.imag)))

    qpr = qext-qsca*g
    qback = (1/xShell2)*(np.abs(np.sum(n1*((-1)**n)*(an-bn)))**2)
    qratio = qback/qsca

    return qext, qsca, qabs, g, qpr, qback, qratio

def CoreShell_ab(mCore,mShell,xCore,xShell):
#  http://pymiescatt.readthedocs.io/en/latest/forwardCS.html#CoreShell_ab
  m = mShell/mCore
  u = mCore*xCore
  v = mShell*xCore
  w = mShell*xShell

  mx = max(np.abs(mCore*xShell),np.abs(mShell*xShell))
  nmax = np.round(2+xShell+4*(xShell**(1/3)))
  nmx = np.round(max(nmax,mx)+16)
  n = np.arange(1,nmax+1)
  nu = n+0.5

  sv = np.sqrt(0.5*np.pi*v)
  sw = np.sqrt(0.5*np.pi*w)
  sy = np.sqrt(0.5*np.pi*xShell)

  pv = sv*jv(nu,v)
  pw = sw*jv(nu,w)
  py = sy*jv(nu,xShell)

  chv = -sv*yv(nu,v)
  chw = -sw*yv(nu,w)
  chy = -sy*yv(nu,xShell)

  p1y = np.append([np.sin(xShell)], [py[0:int(nmax)-1]])
  ch1y = np.append([np.cos(xShell)], [chy[0:int(nmax)-1]])
  gsy = py-(0+1.0j)*chy
  gs1y = p1y-(0+1.0j)*ch1y

  # B&H Equation 4.89
  Dnu = np.zeros((int(nmx)),dtype=complex)
  Dnv = np.zeros((int(nmx)),dtype=complex)
  Dnw = np.zeros((int(nmx)),dtype=complex)
  for i in range(int(nmx)-1,1,-1):
    Dnu[i-1] = i/u-1/(Dnu[i]+i/u)
    Dnv[i-1] = i/v-1/(Dnv[i]+i/v)
    Dnw[i-1] = i/w-1/(Dnw[i]+i/w)

  Du = Dnu[1:int(nmax)+1]
  Dv = Dnv[1:int(nmax)+1]
  Dw = Dnw[1:int(nmax)+1]

  uu = m*Du-Dv
  vv = Du/m-Dv
  fv = pv/chv

  dns = ((uu*fv/pw)/(uu*(pw-chw*fv)+(pw/pv)/chv))+Dw
  gns = ((vv*fv/pw)/(vv*(pw-chw*fv)+(pw/pv)/chv))+Dw
  a1 = dns/mShell+n/xShell
  b1 = mShell*gns+n/xShell

  an = (py*a1-p1y)/(gsy*a1-gs1y)
  bn = (py*b1-p1y)/(gsy*b1-gs1y)

  return an, bn

def CoreShellScatteringFunction(mCore,mShell,wavelength,dCore,dShell,minAngle=0, maxAngle=180, angularResolution=0.5, normed=False):
#  http://pymiescatt.readthedocs.io/en/latest/forwardCS.html#CoreShellScatteringFunction
  xCore = np.pi*dCore/wavelength
  xShell = np.pi*dShell/wavelength
  theta = np.linspace(minAngle,maxAngle,int((maxAngle-minAngle)/angularResolution))*np.pi/180
  thetaSteps = len(theta)
  SL = np.zeros(thetaSteps)
  SR = np.zeros(thetaSteps)
  SU = np.zeros(thetaSteps)
  for j in range(thetaSteps):
    u = np.cos(theta[j])
    S1,S2 = CoreShellS1S2(mCore,mShell,xCore,xShell,u)
    SL[j] = (np.sum((np.conjugate(S1)*S1))).real
    SR[j] = (np.sum((np.conjugate(S2)*S2))).real
    SU[j] = (SR[j]+SL[j])/2
  if normed:
    SL /= np.max(SL)
    SR /= np.max(SR)
    SU /= np.max(SU)
  return theta,SL,SR,SU

def CoreShellS1S2(mCore,mShell,xCore,xShell,mu):
#  http://pymiescatt.readthedocs.io/en/latest/forwardCS.html#CoreShellS1S2
  nmax = np.round(2+xShell+4*(xShell**(1/3)))
  an,bn = CoreShell_ab(mCore,mShell,xCore,xShell)
  pin,taun = MiePiTau(mu,nmax)
  n = np.arange(1,int(nmax)+1)
  n2 = (2*n+1)/(n*(n+1))
  pin *= n2
  taun *= n2
  S1=np.sum(an*np.conjugate(pin))+np.sum(bn*np.conjugate(taun))
  S2=np.sum(an*np.conjugate(taun))+np.sum(bn*np.conjugate(pin))
  return S1,S2

def CoreShellMatrixElements(mCore,mShell,xCore,xShell,mu):
#  http://pymiescatt.readthedocs.io/en/latest/forwardCS.html#CoreShellMatrixElements
  S1,S2 = CoreShellS1S2(mCore,mShell,xCore,xShell,mu)
  S11 = 0.5*(np.abs(S2)**2+np.abs(S1)**2)
  S12 = 0.5*(np.abs(S2)**2-np.abs(S1)**2)
  S33 = 0.5*(np.conjugate(S2)*S1+S2*np.conjugate(S1))
  S34 = 0.5j*(S1*np.conjugate(S2)-S2*np.conjugate(S1))
  return S11, S12, S33, S34
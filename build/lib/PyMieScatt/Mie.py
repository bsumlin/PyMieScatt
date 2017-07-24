# -*- coding: utf-8 -*-

import numpy as np
from scipy.special import jv, yv
from scipy.integrate import trapz

def coerceDType(d):
  if type(d) is not np.ndarray:
    return np.array(d)
  else:
    return d

def MieQ(m, wavelength, diameter, asDict=False):
  # Calculates extinction, scattering, absorption and backscatter efficiencies,
  # as well as asymmetry parameter and qratio for a given m = n+ik, wavelength
  # (in nm) and particle diameter (in nm).
  x = np.pi*diameter/wavelength
  if x==0:
    return 0, 0, 0, 1.5, 0, 0, 0
  elif x>0:
    nmax = np.round(2+x+4*(x**(1/3)))
    n = np.arange(1,nmax+1)
    n1 = 2*n+1
    n2 = n*(n+2)/(n+1)
    n3 = n1/(n*(n+1))
    x2 = x**2

    an,bn = Mie_ab(m,x)

    qext = (2/x2)*np.sum(n1*(an.real+bn.real))
    qsca = (2/x2)*np.sum(n1*(an.real**2+an.imag**2+bn.real**2+bn.imag**2))
    qabs = qext-qsca

    g1 = [an.real[1:int(nmax)],an.imag[1:int(nmax)],bn.real[1:int(nmax)],bn.imag[1:int(nmax)]]
    g1 = [np.append(x, 0.0) for x in g1]
    g = (4/(qsca*x2))*np.sum((n2*(an.real*g1[0]+an.imag*g1[1]+bn.real*g1[2]+bn.imag*g1[3]))+(n3*(an.real*bn.real+an.imag*bn.imag)))

    qpr = qext-qsca*g
    qback = (1/x2)*(np.abs(np.sum(n1*((-1)**n)*(an-bn)))**2)
    qratio = qback/qsca
    if asDict:
      return dict(Qext=qext,Qsca=qsca,Qabs=qabs,g=g,Qpr=qpr,Qback=qback,Qratio=qratio)
    else:
      return qext, qsca, qabs, g, qpr, qback, qratio

def Mie_ab(m,x):
  mx = m*x
  nmax = np.round(2+x+4*(x**(1/3)))
  nmx = np.round(max(nmax,np.abs(mx))+16)
  n = np.arange(1,nmax+1)
  nu = n + 0.5

  sx = np.sqrt(0.5*np.pi*x)
  px = sx*jv(nu,x)

  p1x = np.append(np.sin(x), px[0:int(nmax)-1])
  chx = -sx*yv(nu,x)

  ch1x = np.append(np.cos(x), chx[0:int(nmax)-1])
  gsx = px-(0+1j)*chx
  gs1x = p1x-(0+1j)*ch1x

  # B&H Equation 4.89
  Dn = np.zeros(int(nmx),dtype=complex)
  for i in range(int(nmx)-1,1,-1):
    Dn[i-1] = (i/mx)-(1/(Dn[i]+i/mx))

  D = Dn[1:int(nmax)+1] # Dn(mx), drop terms beyond nMax
  da = D/m+n/x
  db = m*D+n/x

  an = (da*px-p1x)/(da*gsx-gs1x)
  bn = (db*px-p1x)/(db*gsx-gs1x)

  return an, bn

def Mie_cd(m,x):
  mx = m*x
  nmax = np.round(2+x+4*(x**(1/3)))
  nmx = np.round(max(nmax,np.abs(mx))+16)
  n = np.arange(1,int(nmax)+1)
  nu = n+0.5

  cnx = np.zeros(int(nmx),dtype=complex)

  for j in np.arange(nmx,1,-1):
    cnx[int(j)-2] = j-mx*mx/(cnx[int(j)-1]+j)

  cnn = np.array([cnx[b] for b in range(0,len(n))])

  jnx = np.sqrt(np.pi/(2*x))*jv(nu, x)
  jnmx = np.sqrt((2*mx)/np.pi)/jv(nu, mx)

  yx = np.sqrt(np.pi/(2*x))*yv(nu, x)
  hx = jnx+(1.0j)*yx

  b1x = np.append(np.sin(x)/x, jnx[0:int(nmax)-1])
  y1x = np.append(-np.cos(x)/x, yx[0:int(nmax)-1])

  hn1x =  b1x+(1.0j)*y1x
  ax = x*b1x-n*jnx
  ahx =  x*hn1x-n*hx

  numerator = jnx*ahx-hx*ax
  c_denominator = ahx-hx*cnn
  d_denominator = m*m*ahx-hx*cnn

  cn = jnmx*numerator/c_denominator
  dn = jnmx*m*numerator/d_denominator

  return cn, dn

def RayleighMieQ(m, wavelength, diameter, asDict=False):
  x = np.pi*diameter/wavelength
  if x==0:
    return 0, 0, 0, 1.5, 0, 0, 0
  elif x>0:
    LL = (m**2-1)/(m**2+2) # Lorentz-Lorenz term
    LLabsSq = np.abs(LL)**2
    qsca = 8*LLabsSq*(x**4)/3 # B&H eq 5.8
    qabs=4*x*LL.imag # B&H eq. 5.11
    qext=qsca+qabs
    qback = 1.5*qsca # B&H eq. 5.9
    qratio = 1.5
    g = 0
    qpr = qext
    if asDict:
      return dict(Qext=qext,Qsca=qsca,Qabs=qabs,g=g,Qpr=qpr,Qback=qback,Qratio=qratio)
    else:
      return qext, qsca, qabs, g, qpr, qback, qratio

def LowFrequencyMieQ(m, wavelength, diameter, asDict=False):
  # A computationally-cheap way to do full Mie calculations in the Rayleigh limit.
  # Only retains two terms in an and bn per pg 131 of B&H
  x = np.pi*diameter/wavelength
  if x==0:
    return 0, 0, 0, 1.5, 0, 0, 0
  elif x>0:
    n = np.arange(1,3)
    n1 = 2*n+1
    n2 = n*(n+2)/(n+1)
    n3 = n1/(n*(n+1))
    x2 = x**2

    an,bn = LowFrequencyMie_ab(m,x)

    qext = (2/x2)*np.sum(n1*(an.real+bn.real))
    qsca = (2/x2)*np.sum(n1*(an.real**2+an.imag**2+bn.real**2+bn.imag**2))
    qabs = qext-qsca

    g1 = [an.real[1:2],an.imag[1:2],bn.real[1:2],bn.imag[1:2]]
    g1 = [np.append(x, 0.0) for x in g1]
    g = (4/(qsca*x2))*np.sum((n2*(an.real*g1[0]+an.imag*g1[1]+bn.real*g1[2]+bn.imag*g1[3]))+(n3*(an.real*bn.real+an.imag*bn.imag)))

    qpr = qext-qsca*g
    qback = (1/x2)*(np.abs(np.sum(n1*((-1)**n)*(an-bn)))**2)
    qratio = qback/qsca

    if asDict:
      return dict(Qext=qext,Qsca=qsca,Qabs=qabs,g=g,Qpr=qpr,Qback=qback,Qratio=qratio)
    else:
      return qext, qsca, qabs, g, qpr, qback, qratio

def LowFrequencyMie_ab(m,x):
  # B&H page 131
  m2 = m**2
  LL = (m**2-1)/(m**2+2)
  x3 = x**3
  x5 = x**5
  x6 = x**6

  a1 = (-2j*x3/3)*LL-(2j*x5/5)*LL*(m2-2)/(m2+2)+(4*x6/9)*(LL**2)
  a2 = (-1j*x5/15)*(m2-1)/(2*m2+3)
  b1 = (-1j*x5/45)*(m2-1)
  b2 = 0+0j
  an = np.append(a1,a2)
  bn = np.append(b1,b2)
  return an,bn

def MieQ_withSizeDistribution(m, wavelength, sizeDistributionDiameterBins, sizeDistribution, asDict=False):
  # Calculates Mie efficencies integrated across a size distribution, with bins
  # given in nm and the distribution given in particles/cm^3. When using a TSI
  # SMPS, be sure that the units of the output are dW. This function also takes
  # a given m=n+ik, wavelength (in nm) and particle diameter (in nm).
  sizeDistributionDiameterBins = coerceDType(sizeDistributionDiameterBins)
  sizeDistribution = coerceDType(sizeDistribution)
  _length = np.size(sizeDistributionDiameterBins)
  Q_ext = np.zeros(_length)
  Q_sca = np.zeros(_length)
  Q_abs = np.zeros(_length)
  Q_pr = np.zeros(_length)
  Q_back = np.zeros(_length)
  Q_ratio = np.zeros(_length)
  g = np.zeros(_length)

  # scaling of 1e-6 to cast in units of inverse megameters
  aSDn = np.pi*((sizeDistributionDiameterBins/2)**2)*sizeDistribution*(1e-6)

  for i in range(_length):
    Q_ext[i], Q_sca[i], Q_abs[i], g[i], Q_pr[i], Q_back[i], Q_ratio[i] = MieQ(m,wavelength,sizeDistributionDiameterBins[i])

  Bext = trapz(Q_ext*aSDn)
  Bsca = trapz(Q_sca*aSDn)
  Babs = Bext-Bsca
  Bback = trapz(Q_back*aSDn)
  Bratio = trapz(Q_ratio*aSDn)
  bigG = trapz(g*Q_sca*aSDn)/trapz(Q_sca*aSDn)
  Bpr = Bext - bigG*Bsca

  if asDict:
    return dict(Bext=Bext, Bsca=Bsca, Babs=Babs, G=bigG, Bpr=Bpr, Bback=Bback, Bratio=Bratio)
  else:
    return Bext, Bsca, Babs, bigG, Bpr, Bback, Bratio

def ScatteringFunction(m, wavelength, diameter, minAngle=0, maxAngle=180, angularResolution=0.5, normed=False):
  x = np.pi*diameter/wavelength
  theta = np.linspace(minAngle,maxAngle,int((maxAngle-minAngle)/angularResolution))*np.pi/180
  thetaSteps = len(theta)
  SL = np.zeros(thetaSteps)
  SR = np.zeros(thetaSteps)
  SU = np.zeros(thetaSteps)
  for j in range(thetaSteps):
    u = np.cos(theta[j])
    S1, S2 = MieS1S2(m,x,u)
    SL[j] = (np.sum(np.conjugate(S1)*S1)).real
    SR[j] = (np.sum(np.conjugate(S2)*S2)).real
    SU[j] = (SR[j]+SL[j])/2
  if normed:
    SL /= np.max(SL)
    SR /= np.max(SR)
    SU /= np.max(SU)
  return theta,SL,SR,SU

def SF_withSizeDistribution(m, wavelength, diameters, bins, maxAngle=180, angularResolution=0.5, space='theta', normed=False):
  theta = np.arange(0,np.pi,angularResolution)
  thetaSteps = len(theta)
  diameters = coerceDType(diameters)
  bins = coerceDType(bins)
  SL = np.zeros(thetaSteps)
  SR = np.zeros(thetaSteps)
  SU = np.zeros(thetaSteps)
  for num,size in zip(diameters,bins):
    if space=='qspace':
      measure,l,r,u = qSpaceScatteringFunction(m,wavelength,size)
    else:
      measure,l,r,u = ScatteringFunction(m,wavelength,size,thetaSteps)
    SL += l*num
    SR += r*num
    SU += u*num
  SL /= np.sum(diameters)
  SR /= np.sum(diameters)
  SU /= np.sum(diameters)
  if normed:
    SL /= np.max(SL)
    SR /= np.max(SR)
    SU /= np.max(SU)
  return measure,SL,SR,SU

def qSpaceScatteringFunction(m,wavelength,diameter,normed=False):
  x = np.pi*diameter/wavelength
  theta = np.linspace(1,3600,num=3600,endpoint=False)*np.pi/3600
  qR = (4*np.pi/wavelength)*np.sin(theta/2)*(diameter/2)
  SL = np.zeros(3600)
  SR = np.zeros(3600)
  SU = np.zeros(3600)
  for j in range(3600):
    u = np.cos(theta[j])
    S1, S2 = MieS1S2(m,x,u)
    SL[j] = (np.sum(np.conjugate(S1)*S1)).real
    SR[j] = (np.sum(np.conjugate(S2)*S2)).real
    SU[j] = (SR[j]+SL[j])/2
  if normed:
    SL /= np.max(SL)
    SR /= np.max(SR)
    SU /= np.max(SU)
  return qR,SL,SR,SU

def MieS1S2(m,x,mu):
  # B&H eqs. 4.74, mu=cos(theta)
  nmax = np.round(2+x+4*np.power(x,1/3))
  an, bn = Mie_ab(m,x)
  pin, taun = MiePiTau(mu,nmax)
  n = np.arange(1,int(nmax)+1)
  n2 = (2*n+1)/(n*(n+1))
  S1 = np.sum(n2*(an*pin+bn*taun))
  S2 = np.sum(n2*(an*taun+bn*pin))
  return S1, S2

def MiePiTau(mu,nmax):
  # B&H eqs. 4.47 modified for proper list indexing
  # mu = cosine of scattering angle
  p = np.zeros(int(nmax))
  t = np.zeros(int(nmax))
  p[0] = 1
  p[1] = 3*mu
  t[0] = mu
  t[1] = 3.0*np.cos(2*np.arccos(mu))
  for n in range(2,int(nmax)):
    p[n] = ((2*n+1)*(mu*p[n-1])-(n+1)*p[n-2])/n
    t[n] = (n+1)*mu*p[n]-(n+2)*p[n-1]
  return p, t

def MatrixElements(m,wavelength,diameter,mu):
  x = np.pi*diameter/wavelength
  # B&H eqs. 4.77, where mu=cos(theta)
  S1,S2 = MieS1S2(m,x,mu)
  S11 = 0.5*(np.abs(S2)**2+np.abs(S1)**2)
  S12 = 0.5*(np.abs(S2)**2-np.abs(S1)**2)
  S33 = 0.5*(np.conjugate(S2)*S1+S2*np.conjugate(S1))
  S34 = 0.5j*(S1*np.conjugate(S2)-S2*np.conjugate(S1))
  return S11, S12, S33, S34

def MieQ_withDiameterRange(m, wavelength, diameterRange=[10,1000], nd=1000, logD=False):
  if logD:
    diameters = np.logspace(np.log10(diameterRange[0]),np.log10(diameterRange[1]),nd)
  else:
    diameters = np.linspace(diameterRange[0],diameterRange[1],nd)
  _qD = [MieQ(m,wavelength,d) for d in diameters]
  qext = [q[0] for q in _qD]
  qsca = [q[1] for q in _qD]
  qabs = [q[2] for q in _qD]
  g = [q[3] for q in _qD]
  qpr = [q[4] for q in _qD]
  qback = [q[5] for q in _qD]
  qratio = [q[6] for q in _qD]
  return diameters, qext, qsca, qabs, g, qpr, qback, qratio

def MieQ_withWavelengthRange(m, diameter, wavelengthRange=[100,1600], nw=1000, logW=False):
  if logW:
    wavelengths = np.logspace(np.log10(wavelengthRange[0]),np.log10(wavelengthRange[1]),nw)
  else:
    wavelengths = np.linspace(wavelengthRange[0],wavelengthRange[1],nw)
  _qD = [MieQ(m,wavelength,diameter) for wavelength in wavelengths]
  qext = [q[0] for q in _qD]
  qsca = [q[1] for q in _qD]
  qabs = [q[2] for q in _qD]
  g = [q[3] for q in _qD]
  qpr = [q[4] for q in _qD]
  qback = [q[5] for q in _qD]
  qratio = [q[6] for q in _qD]
  return wavelengths, qext, qsca, qabs, g, qpr, qback, qratio

def MieQ_withSizeParameterRange(m, xRange=[1,10], nx=1000, logX=False):
  if logX:
    xValues = list(np.logspace(np.log10(xRange[0]),np.log10(xRange[1]),nx))
  else:
    xValues = list(np.linspace(xRange[0],xRange[1], nx))
  dValues = [1000*x/np.pi for x in xValues]
  _qD = [MieQ(m,1000,d) for d in dValues]
  qext = [q[0] for q in _qD]
  qsca = [q[1] for q in _qD]
  qabs = [q[2] for q in _qD]
  g = [q[3] for q in _qD]
  qpr = [q[4] for q in _qD]
  qback = [q[5] for q in _qD]
  qratio = [q[6] for q in _qD]
  return xValues, qext, qsca, qabs, g, qpr, qback, qratio

def MieQ_withLognormalDistribution(m,wavelength,geoStdDev,geoMean,numberOfParticles,numberOfBins=1000,lower=1,upper=1000,returnDistribution=False,asDict=False):
  def Lognormal(GSD,numberOfParticles,geometricMean,diameters,numberOfBins):
    # Lognormal distribution per Friedlander eq. 1.27
    n_r = (numberOfParticles/np.sqrt(2*np.pi))*(1/np.log(GSD))*(1/diameters)*(np.exp(-1.0*np.square(np.log(diameters)-np.log(geometricMean))/(2*np.square(np.log(GSD)))))
    return n_r
  diameters = np.linspace(lower,upper,numberOfBins)
  nd = Lognormal(geoStdDev,numberOfParticles,geoMean,diameters,numberOfBins)
  Bext, Bsca, Babs, bigG, Bpr, Bback, Bratio = MieQ_withSizeDistribution(m,wavelength,diameters,nd)
  if returnDistribution:
    if asDict==True:
      return dict(Bext=Bext, Bsca=Bsca, Babs=Babs, bigG=bigG, Bpr=Bpr, Bback=Bback, Bratio=Bratio), diameters, nd
    else:
      return Bext, Bsca, Babs, bigG, Bpr, Bback, Bratio, diameters, nd
  else:
    if asDict==True:
      return dict(Bext=Bext, Bsca=Bsca, Babs=Babs, bigG=bigG, Bpr=Bpr, Bback=Bback, Bratio=Bratio)
    else:
      return Bext, Bsca, Babs, bigG, Bpr, Bback, Bratio


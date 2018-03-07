# -*- coding: utf-8 -*-
# http://pymiescatt.readthedocs.io/en/latest/forward.html
import numpy as np
from scipy.special import jv, yv
from scipy.integrate import trapz
import warnings

def coerceDType(d):
  if type(d) is not np.ndarray:
    return np.array(d)
  else:
    return d

def MieQ(m, wavelength, diameter, asDict=False, asCrossSection=False):
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#MieQ
  x = np.pi*diameter/wavelength
  if x==0:
    return 0, 0, 0, 1.5, 0, 0, 0
  elif x<=0.05:
    return RayleighMieQ(m, wavelength, diameter, asDict)
  elif x>0.05:
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

    g1 = [an.real[1:int(nmax)],
          an.imag[1:int(nmax)],
          bn.real[1:int(nmax)],
          bn.imag[1:int(nmax)]]
    g1 = [np.append(x, 0.0) for x in g1]
    g = (4/(qsca*x2))*np.sum((n2*(an.real*g1[0]+an.imag*g1[1]+bn.real*g1[2]+bn.imag*g1[3]))+(n3*(an.real*bn.real+an.imag*bn.imag)))

    qpr = qext-qsca*g
    qback = (1/x2)*(np.abs(np.sum(n1*((-1)**n)*(an-bn)))**2)
    qratio = qback/qsca
    if asCrossSection:
      css = np.pi*(diameter/2)**2
      cext = css*qext
      csca = css*qsca
      cabs = css*qabs
      cpr = css*qpr
      cback = css*qback
      cratio = css*qratio
      if asDict:
        return dict(Cext=cext,Csca=csca,Cabs=cabs,g=g,Cpr=cpr,Cback=cback,Cratio=cratio)
      else:
        return cext, csca, cabs, g, cpr, cback, cratio
    else:
      if asDict:
        return dict(Qext=qext,Qsca=qsca,Qabs=qabs,g=g,Qpr=qpr,Qback=qback,Qratio=qratio)
      else:
        return qext, qsca, qabs, g, qpr, qback, qratio

def Mie_ab(m,x):
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#Mie_ab
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
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#Mie_cd
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

def RayleighMieQ(m, wavelength, diameter, asDict=False, asCrossSection=False):
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#RayleighMieQ
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
    if asCrossSection:
      css = np.pi*(diameter/2)**2
      cext = css*qext
      csca = css*qsca
      cabs = css*qabs
      cpr = css*qpr
      cback = css*qback
      cratio = css*qratio
      if asDict:
        return dict(Cext=cext,Csca=csca,Cabs=cabs,g=g,Cpr=cpr,Cback=cback,Cratio=cratio)
      else:
        return cext, csca, cabs, g, cpr, cback, cratio
    else:
      if asDict:
        return dict(Qext=qext,Qsca=qsca,Qabs=qabs,g=g,Qpr=qpr,Qback=qback,Qratio=qratio)
      else:
        return qext, qsca, qabs, g, qpr, qback, qratio

def AutoMieQ(m, wavelength, diameter, crossover=0.01, asDict=False, asCrossSection=False):
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#AutoMieQ
  x = np.pi*diameter/wavelength
  if x==0:
    return 0, 0, 0, 1.5, 0, 0, 0
  elif x<crossover:
    return RayleighMieQ(m, wavelength, diameter, asDict=asDict, asCrossSection=asCrossSection)
  else:
    return MieQ(m, wavelength, diameter, asDict=asDict, asCrossSection=asCrossSection)

def LowFrequencyMieQ(m, wavelength, diameter, asDict=False, asCrossSection=False):
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#LowFrequencyMieQ
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

    if asCrossSection:
      css = np.pi*(diameter/2)**2
      cext = css*qext
      csca = css*qsca
      cabs = css*qabs
      cpr = css*qpr
      cback = css*qback
      cratio = css*qratio
      if asDict:
        return dict(Cext=cext,Csca=csca,Cabs=cabs,g=g,Cpr=cpr,Cback=cback,Cratio=cratio)
      else:
        return cext, csca, cabs, g, cpr, cback, cratio
    else:
      if asDict:
        return dict(Qext=qext,Qsca=qsca,Qabs=qabs,g=g,Qpr=qpr,Qback=qback,Qratio=qratio)
      else:
        return qext, qsca, qabs, g, qpr, qback, qratio

def LowFrequencyMie_ab(m,x):
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#LowFrequencyMie_ab
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

def AutoMie_ab(m,x):
  if x<0.5:
    return LowFrequencyMie_ab(m,x)
  else:
    return Mie_ab(m,x)

def Mie_SD(m, wavelength, dp, ndp, interpolate=False, asDict=False):
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#Mie_SD
  dp = coerceDType(dp)
  ndp = coerceDType(ndp)
  _length = np.size(dp)
  Q_ext = np.zeros(_length)
  Q_sca = np.zeros(_length)
  Q_abs = np.zeros(_length)
  Q_pr = np.zeros(_length)
  Q_back = np.zeros(_length)
  Q_ratio = np.zeros(_length)
  g = np.zeros(_length)
  
  # scaling of 1e-6 to cast in units of inverse megameters - see docs
  aSDn = np.pi*((dp/2)**2)*ndp*(1e-6)
#  _logdp = np.log10(dp)

  for i in range(_length):
    Q_ext[i], Q_sca[i], Q_abs[i], g[i], Q_pr[i], Q_back[i], Q_ratio[i] = AutoMieQ(m,wavelength,dp[i])

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

def ScatteringFunction(m, wavelength, diameter, minAngle=0, maxAngle=180, angularResolution=0.5, space='theta', angleMeasure='radians', normalization=None):
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#ScatteringFunction
  x = np.pi*diameter/wavelength

  _steps = int(1+(maxAngle-minAngle)/angularResolution) # default 361

  if angleMeasure in ['radians','RADIANS','rad','RAD']:
    adjust = np.pi/180
  elif angleMeasure in ['gradians','GRADIANS','grad','GRAD']:
    adjust = 1/200
  else:
    adjust = 1

  if space in ['q','qspace','QSPACE','qSpace']:
    _steps *= 10
    if minAngle==0:
      minAngle = 1e-5
    measure = np.logspace(np.log10(minAngle),np.log10(maxAngle),_steps)*np.pi/180
    _q = True
  else:
    measure = np.linspace(minAngle,maxAngle,_steps)*adjust
    _q = False
  if x == 0:
    return measure,0,0,0
  _measure = np.linspace(minAngle,maxAngle,_steps)*np.pi/180
  SL = np.zeros(_steps)
  SR = np.zeros(_steps)
  SU = np.zeros(_steps)
  for j in range(_steps):
    u = np.cos(_measure[j])
    S1, S2 = MieS1S2(m,x,u)
    SL[j] = (np.sum(np.conjugate(S1)*S1)).real
    SR[j] = (np.sum(np.conjugate(S2)*S2)).real
    SU[j] = (SR[j]+SL[j])/2
  if normalization in ['m','M','max','MAX']:
    SL /= np.max(SL)
    SR /= np.max(SR)
    SU /= np.max(SU)
  elif normalization in ['t','T','total','TOTAL']:
    SL /= trapz(SL,measure)
    SR /= trapz(SR,measure)
    SU /= trapz(SU,measure)
  if _q:
    measure = (4*np.pi/wavelength)*np.sin(measure/2)*(diameter/2)
  return measure,SL,SR,SU

def SF_SD(m, wavelength, dp, ndp, minAngle=0, maxAngle=180, angularResolution=0.5, space='theta', angleMeasure='radians', normalization=None):
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#SF_SD
  _steps = int(1+(maxAngle-minAngle)/angularResolution)
  ndp = coerceDType(ndp)
  dp = coerceDType(dp)
  SL = np.zeros(_steps)
  SR = np.zeros(_steps)
  SU = np.zeros(_steps)
  kwargs = {'minAngle':minAngle,
            'maxAngle':maxAngle,
            'angularResolution':angularResolution,
            'space':space,
            'normalization':None}
  for n,d in zip(ndp,dp):
    measure,l,r,u = ScatteringFunction(m,wavelength,d,**kwargs)
    SL += l*n
    SR += r*n
    SU += u*n
  if normalization in ['n','N','number','particles']:
    _n = trapz(ndp,dp)
    SL /= _n
    SR /= _n
    SU /= _n
  elif normalization in ['m','M','max','MAX']:
    SL /= np.max(SL)
    SR /= np.max(SR)
    SU /= np.max(SU)
  elif normalization in ['t','T','total','TOTAL']:
    SL /= trapz(SL,measure)
    SR /= trapz(SR,measure)
    SU /= trapz(SU,measure)
  return measure,SL,SR,SU

def MieS1S2(m,x,mu):
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#MieS1S2
  nmax = np.round(2+x+4*np.power(x,1/3))
  an, bn = AutoMie_ab(m,x)
  pin, taun = MiePiTau(mu,nmax)
  n = np.arange(1,int(nmax)+1)
  n2 = (2*n+1)/(n*(n+1))
  S1 = np.sum(n2[0:len(an)]*(an*pin[0:len(an)]+bn*taun[0:len(bn)]))
  S2 = np.sum(n2[0:len(an)]*(an*taun[0:len(an)]+bn*pin[0:len(bn)]))
  return S1, S2

def MiePiTau(mu,nmax):
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#MiePiTau
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
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#MatrixElements
  x = np.pi*diameter/wavelength
  # B&H eqs. 4.77, where mu=cos(theta)
  S1,S2 = MieS1S2(m,x,mu)
  S11 = 0.5*(np.abs(S2)**2+np.abs(S1)**2)
  S12 = 0.5*(np.abs(S2)**2-np.abs(S1)**2)
  S33 = 0.5*(np.conjugate(S2)*S1+S2*np.conjugate(S1))
  S34 = 0.5j*(S1*np.conjugate(S2)-S2*np.conjugate(S1))
  return S11, S12, S33, S34

def MieQ_withDiameterRange(m, wavelength, diameterRange=(10,1000), nd=1000, logD=False):
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#MieQ_withDiameterRange
  if logD:
    diameters = np.logspace(np.log10(diameterRange[0]),np.log10(diameterRange[1]),nd)
  else:
    diameters = np.linspace(diameterRange[0],diameterRange[1],nd)
  _qD = [AutoMieQ(m,wavelength,diameter) for diameter in diameters]
  qext = np.array([q[0] for q in _qD])
  qsca = np.array([q[1] for q in _qD])
  qabs = np.array([q[2] for q in _qD])
  g = np.array([q[3] for q in _qD])
  qpr = np.array([q[4] for q in _qD])
  qback = np.array([q[5] for q in _qD])
  qratio = np.array([q[6] for q in _qD])
  return diameters, qext, qsca, qabs, g, qpr, qback, qratio

def MieQ_withWavelengthRange(m, diameter, wavelengthRange=(100,1600), nw=1000, logW=False):
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#MieQ_withWavelengthRange
  if type(m) == complex and len(wavelengthRange)==2:
    if logW:
      wavelengths = np.logspace(np.log10(wavelengthRange[0]),np.log10(wavelengthRange[1]),nw)
    else:
      wavelengths = np.linspace(wavelengthRange[0],wavelengthRange[1],nw)
    _qD = [AutoMieQ(m,wavelength,diameter) for wavelength in wavelengths]
  elif type(m) in [np.ndarray,list,tuple] and len(wavelengthRange)==len(m):
    wavelengths=wavelengthRange
    _qD = [MieQ(emm,wavelength,diameter) for emm,wavelength in zip(m,wavelengths)]
  else:
    warnings.warn("Error: the size of the input data is minmatched. Please examine your inputs and try again.")
    return

  qext = np.array([q[0] for q in _qD])
  qsca = np.array([q[1] for q in _qD])
  qabs = np.array([q[2] for q in _qD])
  g = np.array([q[3] for q in _qD])
  qpr = np.array([q[4] for q in _qD])
  qback = np.array([q[5] for q in _qD])
  qratio = np.array([q[6] for q in _qD])
  return wavelengths, qext, qsca, qabs, g, qpr, qback, qratio

def MieQ_withSizeParameterRange(m, xRange=(1,10), nx=1000, logX=False):
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#MieQ_withSizeParameterRange
  if logX:
    xValues = list(np.logspace(np.log10(xRange[0]),np.log10(xRange[1]),nx))
  else:
    xValues = list(np.linspace(xRange[0],xRange[1], nx))
  dValues = [1000*x/np.pi for x in xValues]
  _qD = [AutoMieQ(m,1000,d) for d in dValues]
  qext = np.array([q[0] for q in _qD])
  qsca = np.array([q[1] for q in _qD])
  qabs = np.array([q[2] for q in _qD])
  g = np.array([q[3] for q in _qD])
  qpr = np.array([q[4] for q in _qD])
  qback = np.array([q[5] for q in _qD])
  qratio = np.array([q[6] for q in _qD])
  return xValues, qext, qsca, qabs, g, qpr, qback, qratio

def Mie_Lognormal(m,wavelength,geoStdDev,geoMean,numberOfParticles,numberOfBins=10000,lower=1,upper=1000,gamma=[1],returnDistribution=False,decomposeMultimodal=False,asDict=False):
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#Mie_Lognormal
  ithPart = lambda gammai, dp, dpgi, sigmagi: (gammai/(np.sqrt(2*np.pi)*np.log(sigmagi)*dp))*np.exp(-(np.log(dp)-np.log(dpgi))**2/(2*np.log(sigmagi)**2))
  dp = np.logspace(np.log10(lower),np.log10(upper),numberOfBins)
  if all([type(x) in [list, tuple, np.ndarray] for x in [geoStdDev, geoMean]]):
    # multimodal
    if len(gamma)==1 and (len(geoStdDev)==len(geoMean)>1):
      # gamma is distributed equally among modes
      gamma = [1 for x in geoStdDev]
      gamma = [float(x/np.sum(gamma)) for x in gamma]
      ndpi = [numberOfParticles*ithPart(g,dp,dpg,sg) for g,dpg,sg in zip(gamma,geoMean,geoStdDev)]
      ndp = np.sum(ndpi,axis=0)
    elif len(gamma)==len(geoStdDev)==len(geoMean):
      # gamma is fully specified for each mode
      gamma = [float(x/np.sum(gamma)) for x in gamma]
      ndpi = [numberOfParticles*ithPart(g,dp,dpg,sg) for g,dpg,sg in zip(gamma,geoMean,geoStdDev)]
      ndp = np.sum(ndpi,axis=0)
    else:
      # user fucked up
      warnings.warn("Not enough parameters to fully specify each mode.")
      return None
  else:
    # unimodal
    decomposeMultimodal = False
    ndp = numberOfParticles*ithPart(1,dp,geoMean,geoStdDev)
  if ndp[-1]>np.max(ndp)/100 or ndp[0]>np.max(ndp)/100:
    warnings.warn("Warning: distribution may not be compact on the specified interval. Consider using a higher upper bound.")
  Bext, Bsca, Babs, bigG, Bpr, Bback, Bratio = Mie_SD(m,wavelength,dp,ndp)
  if returnDistribution:
    if decomposeMultimodal:
      if asDict==True:
        return dict(Bext=Bext, Bsca=Bsca, Babs=Babs, bigG=bigG, Bpr=Bpr, Bback=Bback, Bratio=Bratio), dp, ndp, ndpi
      else:
        return Bext, Bsca, Babs, bigG, Bpr, Bback, Bratio, dp, ndp, ndpi
    else:
      if asDict==True:
        return dict(Bext=Bext, Bsca=Bsca, Babs=Babs, bigG=bigG, Bpr=Bpr, Bback=Bback, Bratio=Bratio), dp, ndp
      else:
        return Bext, Bsca, Babs, bigG, Bpr, Bback, Bratio, dp, ndp
  else:
    if asDict==True:
      return dict(Bext=Bext, Bsca=Bsca, Babs=Babs, bigG=bigG, Bpr=Bpr, Bback=Bback, Bratio=Bratio)
    else:
      return Bext, Bsca, Babs, bigG, Bpr, Bback, Bratio
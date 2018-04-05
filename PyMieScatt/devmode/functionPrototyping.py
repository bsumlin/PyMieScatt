import numpy as np
from PyMieScatt import RayleighMieQ
from scipy.special import jv, yv

def MieQ(m, wavelength, diameter, nMedium=1, asDict=False, asCrossSection=False):
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#MieQ
  wavelength /= nMedium
  m /= nMedium
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

m = 1.5+0.5j
nMedium = 1.3

q1 = MieQ(m, 530, 500, asDict=True)
q2 = MieQ(m, 530, 500, nMedium=nMedium,asDict=True)

m /= nMedium
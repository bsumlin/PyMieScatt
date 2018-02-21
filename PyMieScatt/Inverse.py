# -*- coding: utf-8 -*-
# http://pymiescatt.readthedocs.io/en/latest/inverse.html
from PyMieScatt.Mie import Mie_ab
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from scipy.ndimage import zoom
from scipy.integrate import trapz
from shapely import geometry
import warnings

def coerceDType(d):
  if type(d) is not np.ndarray:
    return np.array(d)
  else:
    return d

def Inversion(Qsca,Qabs,wavelength,diameter,nMin=1,nMax=3,kMin=0.001,kMax=1,scatteringPrecision=0.010,absorptionPrecision=0.010,spaceSize=120,interp=2):
  error = lambda measured,calculated: np.abs((calculated-measured)/measured)
  
  nRange = np.linspace(nMin,nMax,spaceSize)
  kRange = np.logspace(np.log10(kMin),np.log10(kMax),spaceSize)
  scaSpace = np.zeros((spaceSize,spaceSize))
  absSpace = np.zeros((spaceSize,spaceSize))

  for ni,n in enumerate(nRange):
    for ki,k in enumerate(kRange):
      _derp = fastMieQ(n+(1j*k),wavelength,diameter)
      scaSpace[ni][ki] = _derp[0]
      absSpace[ni][ki] = _derp[1]
  if interp is not None:
    nRange = zoom(nRange,interp)
    kRange = zoom(kRange,interp)
    scaSpace = zoom(scaSpace,interp)
    absSpace = zoom(absSpace,interp)
    
  scaSolutions = np.where(np.logical_and(Qsca*(1-scatteringPrecision)<scaSpace, scaSpace<Qsca*(1+scatteringPrecision)))
  absSolutions = np.where(np.logical_and(Qabs*(1-absorptionPrecision)<absSpace, absSpace<Qabs*(1+absorptionPrecision)))

  validScattering = nRange[scaSolutions[0]]+1j*kRange[scaSolutions[1]]
  validAbsorption = nRange[absSolutions[0]]+1j*kRange[absSolutions[1]]
  
  solution = np.intersect1d(validScattering,validAbsorption)
#  errors = [error()]

  return solution

def Inversion_SD(Bsca,Babs,wavelength,dp,ndp,nMin=1,nMax=3,kMin=0,kMax=1,scatteringPrecision=0.001,absorptionPrecision=0.001,spaceSize=40,interp=2):
  dp = coerceDType(dp)
  ndp = coerceDType(ndp)

  nRange = np.linspace(nMin,nMax,spaceSize)
  kRange = np.linspace(kMin,kMax,spaceSize)
  scaSpace = np.zeros((spaceSize,spaceSize))
  absSpace = np.zeros((spaceSize,spaceSize))

  for ni,n in enumerate(nRange):
    for ki,k in enumerate(kRange):
      _derp = fastMie_SD(n+(1j*k),wavelength,dp,ndp)
      scaSpace[ni][ki] = _derp[0]
      absSpace[ni][ki] = _derp[1]
  if interp is not None:
    nRange = zoom(nRange,interp)
    kRange = zoom(kRange,interp)
    scaSpace = zoom(scaSpace,interp)
    absSpace = zoom(absSpace,interp)

  scaSolutions = np.where(np.logical_and(Bsca*(1-scatteringPrecision)<scaSpace, scaSpace<Bsca*(1+scatteringPrecision)))
  absSolutions = np.where(np.logical_and(Babs*(1-absorptionPrecision)<absSpace, absSpace<Babs*(1+absorptionPrecision)))

  validScattering = nRange[scaSolutions[0]]+1j*kRange[scaSolutions[1]]
  validAbsorption = nRange[absSolutions[0]]+1j*kRange[absSolutions[1]]

  return np.intersect1d(validScattering,validAbsorption)

def ContourIntersection(Qsca,Qabs,wavelength,diameter,Qback=None,nMin=1,nMax=3,kMin=0.00001,kMax=1,gridPoints=100,interpolationFactor=2,maxError=0.005,fig=None,ax=None,axisOption=0):
#  http://pymiescatt.readthedocs.io/en/latest/inverse.html#ContourIntersection
  error = lambda measured,calculated: np.abs((calculated-measured)/measured)
  if Qback is not None:
    if gridPoints*interpolationFactor<400:
      gridPoints = 2*gridPoints
  labels = []
  incErrors = False
  if type(Qsca) in [list, tuple, np.ndarray]:
    incErrors = True
    scaError = Qsca[1]
    Qsca = Qsca[0]
    labels.append("Qsca = {b:1.3f}±{e:1.3f}".format(b=Qsca,e=scaError))
  else:
    scaError = None
    labels.append("Qsca = {b:1.3f}".format(b=Qsca))
  if type(Qabs) in [list, tuple, np.ndarray]:
    absError = Qabs[1]
    Qabs = Qabs[0]
    labels.append("Qabs = {b:1.3f}±{e:1.3f}".format(b=Qabs,e=absError))
  else:
    absError = None
    labels.append("Qabs = {b:1.3f}".format(b=Qabs))
  if type(Qback) in [list, tuple, np.ndarray]:
    backError = Qback[1]
    Qback = Qback[0]
    labels.append("Qback = {b:1.3f}±{e:1.3f}".format(b=Qback,e=backError))
  elif Qback is not None:
    backError = None
    labels.append("Qback - {b:1.3f}".format(b=Qback))
  else:
    backError = None

  nRange = np.linspace(nMin,nMax,gridPoints)
  kRange = np.logspace(np.log10(kMin),np.log10(kMax),gridPoints)
  QscaList, QabsList, QbackList = [], [], []
  for n in nRange:
    s, a, b = [], [], []
    for k in kRange:
      m = n+k*1.0j
      _Qsca,_Qabs,_Qback = fastMieQ(m,wavelength,diameter)
      s.append(_Qsca)
      a.append(_Qabs)
      b.append(_Qback)
    QscaList.append(s)
    QabsList.append(a)
    QbackList.append(b)
  QscaList = zoom(np.transpose(np.array(QscaList)),interpolationFactor)
  QabsList = zoom(np.transpose(np.array(QabsList)),interpolationFactor)
  QbackList = zoom(np.transpose(np.array(QbackList)),interpolationFactor)

  n = zoom(nRange,interpolationFactor)
  k = zoom(kRange,interpolationFactor)

  if fig is None and ax is None:
    fig, ax = plt.subplots()
  elif fig is None:
    fig = ax.get_figure()
  elif ax is None:
    ax = fig.gca()

  scaLevels = np.array([Qsca])
  absLevels = np.array([Qabs])

  if Qback is not None:
    backLevels = np.array([Qback])
    if backError is not None:
      backErrorLevels = np.array([Qback+x for x in [-backError,backError]])

  scaChart = ax.contour(n,k,QscaList,scaLevels,origin='lower',linestyles='dashdot',linewidths=1.5,colors=('red'))
  absChart = ax.contour(n,k,QabsList,absLevels,origin='lower',linewidths=1.5,colors=('blue'))
  if Qback is not None:
    backChart = ax.contour(n,k,QbackList,backLevels,origin='lower',linestyles='dotted',linewidths=1.5,colors=('green'))

  if scaError is not None:
    scaErrorLevels = np.array([Qsca+x for x in [-scaError, scaError]])
    ax.contourf(n,k,QscaList,scaErrorLevels,origin='lower',colors=('red'),alpha=0.15)
    ax.contour(n,k,QscaList,scaErrorLevels,origin='lower',linewidths=0.5,colors=('red'),alpha=0.5)

  if absError is not None:
    absErrorLevels = np.array([Qabs+x for x in [-absError, absError]])
    ax.contourf(n,k,QabsList,absErrorLevels,origin='lower',colors=('blue'),alpha=0.15)
    ax.contour(n,k,QabsList,absErrorLevels,origin='lower',linewidths=0.5,colors=('blue'),alpha=0.5)

  if backError is not None:
    backErrorLevels = np.array([Qback+x for x in [-backError, backError]])
    ax.contourf(n,k,QbackList,backErrorLevels,origin='lower',colors=('green'),alpha=0.15)
    ax.contour(n,k,QbackList,backErrorLevels,origin='lower',linewidths=0.5,colors=('green'),alpha=0.5)

  m1 = find_intersections(scaChart,absChart)
  print(m1)
  if Qback is not None:
    m2 = find_intersections(scaChart,backChart)
    r1 = [np.round(x+y*1j,2) for x,y in zip(m1[0],m1[1])]
    r2 = [np.round(x+y*1j,2) for x,y in zip(m2[0],m2[1])]
    m_sol = list(set(r1).intersection(r2))
    nSolution,kSolution = [xx.real for xx in m_sol],[xx.imag for xx in m_sol]
  else:
    nSolution,kSolution = m1[0],m1[1]

  if type(nSolution)==np.float64:
    solutionSet = [nSolution + (0+1j)*kSolution]
  else:
    solutionSet = [(x+y*1j) for x,y in zip(nSolution,kSolution)]

  forwardCalculations = []
  for s in solutionSet:
    _s,_a,_ = fastMieQ(s,wavelength,diameter)
    forwardCalculations.append([_s,_a])
  solutionErrors = []
  for f in forwardCalculations:
    solutionErrors.append([error(f[0],Qsca),error(f[1],Qabs)])

  solutionSet = np.array(solutionSet)
  forwardCalculations = np.array(forwardCalculations)
  solutionErrors = np.array(solutionErrors)

  proper = solutionErrors <= maxError
  solution = []
  for r,c in proper:
    if r and c:
      solution.append(True)
    else:
      solution.append(False)

  solutionSet = solutionSet[solution]
  forwardCalculations = forwardCalculations[solution]
  solutionErrors = solutionErrors[solution]
  nSolutionsToPlot,kSolutionsToPlot = [x.real for x in solutionSet],[x.imag for x in solutionSet]

  ax.scatter(nSolutionsToPlot,kSolutionsToPlot,marker='o',s=128,linewidth=1.5,edgecolor='k',facecolor='none',zorder=3)
  ax.scatter(nSolutionsToPlot,kSolutionsToPlot,marker='o',s=128,linewidth=0,edgecolor='none',facecolor='c',zorder=1,alpha=0.25)
  
  for x,y,s in zip(nSolutionsToPlot,kSolutionsToPlot,solutionErrors):
    ax.axhline(y,linewidth=0.5,alpha=0.5,zorder=0)
    ax.axvline(x,linewidth=0.5,alpha=0.5,zorder=0)

  ax.set_xlabel('n',fontsize=16)
  ax.set_ylabel('k',fontsize=16)

  ax.set_xlim((np.min(nRange),np.max(nRange)))
  ax.set_ylim((np.min(kRange),np.max(kRange)))
  ax.tick_params(which='both',direction='in')

  if axisOption == 0:
    if max(kSolutionsToPlot) <= 0.5 or kMax <= 1:
      ax.set_yscale('log')
    else:
      ax.set_yscale('linear')
  elif axisOption == 1:
    ax.set_xscale('linear')
    ax.set_yscale('linear')
  elif axisOption == 2:
    ax.set_yscale('log')
  elif axisOption == 3:
    ax.set_xscale('log')
  elif axisOption == 4:
    ax.set_xscale('log')
    ax.set_yscale('log')
  else:
    pass

  _c = ax.get_children()
  if Qback is None:
    if incErrors:
      # no Qback, with error bounds
      graphElements = {'Qsca':_c[0],'Qabs':_c[1], # contours
                       'QscaErrFill':_c[2],'QscaErrOutline1':_c[3],'QscaErrOutline2':_c[4],
                       'QabsErrFill':_c[5],'QabsErrOutline1':_c[6],'QabsErrOutline2':_c[7],
                       'SolMark':_c[8],'SolFill':_c[9], # the circly thingies at each solutions
                       'CrosshairsH':_c[10:-10:2],'CrosshairsV':_c[11:-10:2], # solution crosshairs
                       'LeftSpine':_c[-10],'RightSpine':_c[-9],'BottomSpine':_c[-8],'TopSpine':_c[-7], # spines
                       'XAxis':_c[-6],'YAxis':_c[-5]} # the axes
    else:
      # no Qback, no error bounds
      graphElements = {'Qsca':_c[0],'Qabs':_c[1], # contours
                       'SolFill':_c[2],'SolMark':_c[3], # the circly thingies at each solutions
                       'CrosshairsH':_c[4:-10:2],'CrosshairsV':_c[5:-10:2], # solution crosshairs
                       'LeftSpine':_c[-10],'RightSpine':_c[-9],'BottomSpine':_c[-8],'TopSpine':_c[-7], # spines
                       'XAxis':_c[-6],'YAxis':_c[-5]} # the axes
                       
  else:
    if incErrors:
      # with Qback and error bounds
      graphElements = {'Qsca':_c[0],'Qabs':_c[1],'Qback':_c[2], # contours
                       'QscaErrFill':_c[4],'QscaErrOutline1':_c[5],'QscaErrOutline2':_c[6],
                       'QabsErrFill':_c[7],'QabsErrOutline1':_c[8],'QabsErrOutline2':_c[9],
                       'SolMark':_c[10],'SolFill':_c[11], # the circly thingies at each solutions
                       'CrosshairsH':_c[12:-10:2],'CrosshairsV':_c[13:-10:2], # solution crosshairs
                       'LeftSpine':_c[-10],'RightSpine':_c[-9],'BottomSpine':_c[-8],'TopSpine':_c[-7], # spines
                       'XAxis':_c[-6],'YAxis':_c[-5]} # the axes
    else:
      # with Qback, no error bounds
      graphElements = {'Qsca':_c[0],'Qabs':_c[1],'Qback':_c[2], # contours
                       'SolFill':_c[3],'SolMark':_c[4], # the circly thingies at each solution
                       'CrosshairsH':_c[5:-10:2],'CrosshairsV':_c[6:-10:2], # solution crosshairs
                       'LeftSpine':_c[-10],'RightSpine':_c[-9],'BottomSpine':_c[-8],'TopSpine':_c[-7], # spines
                       'XAxis':_c[-6],'YAxis':_c[-5]} # the axes

  return solutionSet,forwardCalculations,solutionErrors, fig, ax, graphElements

def ContourIntersection_SD(Bsca,Babs,wavelength,dp,ndp,nMin=1,nMax=3,kMin=0.00001,kMax=1,Bback=None,gridPoints=60,interpolationFactor=2,maxError=0.005,fig=None,ax=None,axisOption=0):
#  http://pymiescatt.readthedocs.io/en/latest/inverse.html#ContourIntersection_SD
  error = lambda measured,calculated: np.abs((calculated-measured)/measured)
  if Bback is not None:
    if gridPoints*interpolationFactor<120:
      gridPoints = 2*gridPoints
  labels = []
  incErrors = False
  if type(Bsca) in [list, tuple, np.ndarray]:
    incErrors = True
    scaError = Bsca[1]
    Bsca = Bsca[0]
    labels.append("Bsca = {b:1.1f}±{e:1.1f}".format(b=Bsca,e=scaError))
  else:
    scaError = None
    labels.append("Bsca = {b:1.1f}".format(b=Bsca))
  if type(Babs) in [list, tuple, np.ndarray]:
    absError = Babs[1]
    Babs = Babs[0]
    labels.append("Babs = {b:1.1f}±{e:1.1f}".format(b=Babs,e=absError))
  else:
    absError = None
    labels.append("Babs = {b:1.1f}".format(b=Babs))
  if type(Bback) in [list, tuple, np.ndarray]:
    backError = Bback[1]
    Bback = Bback[0]
    labels.append("Bback = {b:1.1f}±{e:1.1f}".format(b=Bback,e=backError))
  elif Bback is not None:
    backError = None
    labels.append("Bback - {b:1.1f}".format(b=Bback))
  else:
    backError = None

  dp = coerceDType(dp)
  ndp = coerceDType(ndp)
  nRange = np.linspace(nMin,nMax,gridPoints)
  kRange = np.logspace(np.log10(kMin),np.log10(kMax),gridPoints)
  BscaList, BabsList, BbackList = [], [], []
  for n in nRange:
    s, a, b = [], [], []
    for k in kRange:
      m = n+k*1.0j
      _Bsca,_Babs,_Bback = fastMie_SD(m,wavelength,dp,ndp)
      s.append(_Bsca)
      a.append(_Babs)
      b.append(_Bback)
    BscaList.append(s)
    BabsList.append(a)
    BbackList.append(b)
  BscaList = zoom(np.transpose(np.array(BscaList)),interpolationFactor)
  BabsList = zoom(np.transpose(np.array(BabsList)),interpolationFactor)
  BbackList = zoom(np.transpose(np.array(BbackList)),interpolationFactor)

  n = zoom(nRange,interpolationFactor)
  k = zoom(kRange,interpolationFactor)

  if fig is None and ax is None:
    fig, ax = plt.subplots()
  elif fig is None:
    fig = ax.get_figure()
  elif ax is None:
    ax = fig.gca()

  scaLevels = np.array([Bsca])
  absLevels = np.array([Babs])

  if Bback is not None:
    backLevels = np.array([Bback])
    if backError is not None:
      backErrorLevels = np.array([Bback+x for x in [-backError,backError]])

  scaChart = ax.contour(n,k,BscaList,scaLevels,origin='lower',linestyles='dashdot',linewidths=1.5,colors=('red'))
  absChart = ax.contour(n,k,BabsList,absLevels,origin='lower',linewidths=1.5,colors=('blue'))
  if Bback is not None:
    backChart = ax.contour(n,k,BbackList,backLevels,origin='lower',linestyles='dotted',linewidths=1.5,colors=('green'))

  if scaError is not None:
    scaErrorLevels = np.array([Bsca+x for x in [-scaError, scaError]])
    ax.contourf(n,k,BscaList,scaErrorLevels,origin='lower',colors=('red'),alpha=0.15)
    ax.contour(n,k,BscaList,scaErrorLevels,origin='lower',linewidths=0.5,colors=('red'),alpha=0.5)

  if absError is not None:
    absErrorLevels = np.array([Babs+x for x in [-absError, absError]])
    ax.contourf(n,k,BabsList,absErrorLevels,origin='lower',colors=('blue'),alpha=0.15)
    ax.contour(n,k,BabsList,absErrorLevels,origin='lower',linewidths=0.5,colors=('blue'),alpha=0.5)

  if backError is not None:
    backErrorLevels = np.array([Bback+x for x in [-backError, backError]])
    ax.contourf(n,k,BbackList,backErrorLevels,origin='lower',colors=('green'),alpha=0.15)
    ax.contour(n,k,BbackList,backErrorLevels,origin='lower',linewidths=0.5,colors=('green'),alpha=0.5)

  m1 = find_intersections(scaChart,absChart)
  print(m1)
  if Bback is not None:
    m2 = find_intersections(scaChart,backChart)
    r1 = [np.round(x+y*1j,2) for x,y in zip(m1[0],m1[1])]
    r2 = [np.round(x+y*1j,2) for x,y in zip(m2[0],m2[1])]
    m_sol = list(set(r1).intersection(r2))
    nSolution,kSolution = [xx.real for xx in m_sol],[xx.imag for xx in m_sol]
  else:
    nSolution,kSolution = m1[0],m1[1]

  if type(nSolution)==np.float64:
    solutionSet = [nSolution + kSolution*1j]
  else:
    solutionSet = [(x+y*1j) for x,y in zip(nSolution,kSolution)]

  forwardCalculations = []
  for s in solutionSet:
    _s,_a,_ = fastMie_SD(s,wavelength,dp,ndp)
    forwardCalculations.append([_s,_a])
  solutionErrors = []
  for f in forwardCalculations:
    solutionErrors.append([error(f[0],Bsca),error(f[1],Babs)])

  solutionSet = np.array(solutionSet)
  forwardCalculations = np.array(forwardCalculations)
  solutionErrors = np.array(solutionErrors)

  proper = solutionErrors <= maxError
  solution = []
  for r,c in proper:
    if r and c:
      solution.append(True)
    else:
      solution.append(False)

  solutionSet = solutionSet[solution]
  forwardCalculations = forwardCalculations[solution]
  solutionErrors = solutionErrors[solution]
  nSolutionsToPlot,kSolutionsToPlot = [x.real for x in solutionSet],[x.imag for x in solutionSet]

  ax.scatter(nSolutionsToPlot,kSolutionsToPlot,marker='o',s=128,linewidth=1.5,edgecolor='k',facecolor='none',zorder=3)
  ax.scatter(nSolutionsToPlot,kSolutionsToPlot,marker='o',s=128,linewidth=0,edgecolor='none',facecolor='c',zorder=1,alpha=0.25)
  
  for x,y,s in zip(nSolutionsToPlot,kSolutionsToPlot,solutionErrors):
    ax.axhline(y,linewidth=0.5,alpha=0.5,zorder=0)
    ax.axvline(x,linewidth=0.5,alpha=0.5,zorder=0)

  ax.set_xlabel('n',fontsize=16)
  ax.set_ylabel('k',fontsize=16)

  ax.set_xlim((np.min(nRange),np.max(nRange)))
  ax.set_ylim((np.min(kRange),np.max(kRange)))
  ax.tick_params(which='both',direction='in')

  if axisOption == 0:
    if max(kSolutionsToPlot) <= 0.5 or kMax <= 1:
      ax.set_yscale('log')
    else:
      ax.set_yscale('linear')
  elif axisOption == 1:
    ax.set_xscale('linear')
    ax.set_yscale('linear')
  elif axisOption == 2:
    ax.set_yscale('log')
  elif axisOption == 3:
    ax.set_xscale('log')
  elif axisOption == 4:
    ax.set_xscale('log')
    ax.set_yscale('log')
  else:
    pass

  _c = ax.get_children()
  if Bback is None:
    if incErrors:
      # no Bback, with error bounds
      graphElements = {'Bsca':_c[0],'Babs':_c[1], # contours
                       'Labels':labels,
                       'BscaErrFill':_c[2],'BscaErrOutline1':_c[3],'BscaErrOutline2':_c[4],
                       'BabsErrFill':_c[5],'BabsErrOutline1':_c[6],'BabsErrOutline2':_c[7],
                       'SolMark':_c[8],'SolFill':_c[9], # the circly thingies at each solutions
                       'CrosshairsH':_c[10:-10:2],'CrosshairsV':_c[11:-10:2], # solution crosshairs
                       'LeftSpine':_c[-10],'RightSpine':_c[-9],'BottomSpine':_c[-8],'TopSpine':_c[-7], # spines
                       'XAxis':_c[-6],'YAxis':_c[-5]} # the axes
    else:
      # no Bback, no error bounds
      graphElements = {'Bsca':_c[0],'Babs':_c[1], # contours
                       'Labels':labels,
                       'SolFill':_c[2],'SolMark':_c[3], # the circly thingies at each solutions
                       'CrosshairsH':_c[4:-10:2],'CrosshairsV':_c[5:-10:2], # solution crosshairs
                       'LeftSpine':_c[-10],'RightSpine':_c[-9],'BottomSpine':_c[-8],'TopSpine':_c[-7], # spines
                       'XAxis':_c[-6],'YAxis':_c[-5]} # the axes
                       
  else:
    if incErrors:
      # with Bback and error bounds
      graphElements = {'Bsca':_c[0],'Babs':_c[1],'Bback':_c[2], # contours
                       'Labels':labels,
                       'BscaErrFill':_c[4],'BscaErrOutline1':_c[5],'BscaErrOutline2':_c[6],
                       'BabsErrFill':_c[7],'BabsErrOutline1':_c[8],'BabsErrOutline2':_c[9],
                       'SolMark':_c[10],'SolFill':_c[11], # the circly thingies at each solutions
                       'CrosshairsH':_c[12:-10:2],'CrosshairsV':_c[13:-10:2], # solution crosshairs
                       'LeftSpine':_c[-10],'RightSpine':_c[-9],'BottomSpine':_c[-8],'TopSpine':_c[-7], # spines
                       'XAxis':_c[-6],'YAxis':_c[-5]} # the axes
    else:
      # with Bback, no error bounds
      graphElements = {'Bsca':_c[0],'Babs':_c[1],'Bback':_c[2], # contours
                       'Labels':labels,
                       'SolFill':_c[3],'SolMark':_c[4], # the circly thingies at each solution
                       'CrosshairsH':_c[5:-10:2],'CrosshairsV':_c[6:-10:2], # solution crosshairs
                       'LeftSpine':_c[-10],'RightSpine':_c[-9],'BottomSpine':_c[-8],'TopSpine':_c[-7], # spines
                       'XAxis':_c[-6],'YAxis':_c[-5]} # the axes

  return solutionSet,forwardCalculations,solutionErrors, fig, ax, graphElements

def find_intersections(A,B):
  X = []
  Y = []
  for p1 in A.collections[0].get_paths():
    for p2 in B.collections[0].get_paths():
      v1 = p1.vertices
      v2 = p2.vertices
      
      try:
        poly1 = geometry.LineString(v1)
        poly2 = geometry.LineString(v2)
        sol = poly1.intersection(poly2)
        for s in sol:
          X.append(s.x)
          Y.append(s.y)
      except:
        pass

  return X,Y

def SurveyIteration(Qsca,Qabs,wavelength,diameter,tolerance=0.0005):
#  http://pymiescatt.readthedocs.io/en/latest/inverse.html#IterativeInversion
  error = lambda measured,calculated: np.abs((calculated-measured)/measured)
  initial_m = Inversion(Qsca,Qabs,wavelength,diameter,scatteringPrecision=0.015,absorptionPrecision=0.015,spaceSize=85,interp=2)
  print(len(initial_m))
  factors = [2.5, 5.0, 10.0, 25.0, 50.0]
  errors = [tolerance*x for x in [25,15,5,3,1]]
  resultM = []
  resultScaErr = []
  resultAbsErr = []
  for m in list(initial_m):
    nResolution = 100
    kResolution = 1000
    trialQsca,trialQabs,trialQback = fastMieQ(m, wavelength, diameter)

    scaError = np.abs(trialQsca-Qsca)/Qsca
    absError = np.abs(trialQabs-Qabs)/Qabs

    for f,e in zip(factors,errors):
      nStep = 10/(f*nResolution)
      flippyFloppy = 0
      while scaError > e:
        if flippyFloppy == 4:
          break
        scaError_prev = scaError
        m = m+nStep
        trialQsca,trialQabs,trialQback = fastMieQ(m,wavelength,diameter)
        scaError = error(Qsca,trialQsca)#error(Qsca,trialQsca)
        if scaError_prev - scaError < 0:
          nStep *= -1
          flippyFloppy += 1

    for f,e in zip(factors,errors):
      kStep = 1.0j/(f*kResolution)
      flippyFloppy = 0
      while (absError > e):
        if flippyFloppy == 4:
          break
        absError_prev = absError
        m = m + kStep
        trialQsca,trialQabs,trialQback = fastMieQ(m,wavelength,diameter)
        absError = error(Qabs,trialQabs)
        if absError_prev-absError < 0:
          kStep *= -1
          flippyFloppy += 1

    scaError = error(Qsca,trialQsca)
    nStep = 0.001
    flippyFloppy = 0
    while (scaError > 0.01):
      if flippyFloppy == 4:
        break
      scaError_prev = scaError
      m += nStep
      trialQsca,trialQabs,trialQback = fastMieQ(m,wavelength,diameter)
      scaError = error(Qsca,trialQsca)
      if scaError_prev - scaError < 0:
        nStep *= -1
        flippyFloppy += 1

    absError = error(Qabs,trialQabs)
    kStep = 0.0+0.00001j
    flippyFloppy = 0
    while absError > 0.005:
      if flippyFloppy == 4:
        break
      absError_prev = absError
      m += kStep
      trialQsca,trialQabs,trialQback = fastMieQ(m,wavelength,diameter)
      absError = error(Qabs,trialQabs)
      if absError_prev - absError < 0:
        kStep *= -1
        flippyFloppy += 1

    scaError = error(Qsca,trialQsca)
    nStep = 0.0001
    flippyFloppy = 0
    while (scaError > 0.005):
      if flippyFloppy == 4:
        break
      scaError_prev = scaError
      m += nStep
      trialQsca,trialQabs,trialQback = fastMieQ(m,wavelength,diameter)
      scaError = error(Qsca,trialQsca)
      if scaError_prev - scaError < 0:
        nStep *= -1
        flippyFloppy += 1

    trialQsca,trialQabs,trialQback = fastMieQ(m,wavelength,diameter)
    scaError = error(Qsca,trialQsca)
    absError = error(Qabs,trialQabs)
    resultM.append(m)
    resultScaErr.append(scaError)
    resultAbsErr.append(absError)

  return resultM, resultScaErr, resultAbsErr

def SurveyIteration_SD(Bsca,Babs,wavelength,dp,ndp,tolerance=0.0005):
#  http://pymiescatt.readthedocs.io/en/latest/inverse.html#IterativeInversion_SD
  dp = coerceDType(dp)
  ndp = coerceDType(ndp)
  error = lambda measured,calculated: np.abs((calculated-measured)/measured)
  initial_m = Inversion_SD(Bsca,Babs,wavelength,dp,ndp,scatteringPrecision=0.08,absorptionPrecision=0.08,spaceSize=40,interp=2)
  print(len(initial_m),flush=True)
  factors = [2.5, 5.0, 10.0, 25.0, 50.0]
  errors = [tolerance*x for x in [25,15,5,3,1]]
  resultM = []
  resultScaErr = []
  resultAbsErr = []

  for m in initial_m:
    nResolution = 100
    kResolution = 1000
    trialBsca,trialBabs,_ = fastMie_SD(m,wavelength,dp,ndp)

    scaError = error(Bsca,trialBsca)
    absError = error(Babs,trialBabs)

    for f,e in zip(factors,errors):
      nStep = 10/(f*nResolution)
      flippyFloppy = 0
      while scaError > e:
        if flippyFloppy == 4:
          break
        scaError_prev = scaError
        m = m+nStep
        trialBsca,trialBabs,_ = fastMie_SD(m,wavelength,dp,ndp)
        scaError = error(Bsca,trialBsca)
        if scaError_prev - scaError < 0:
          nStep *= -1
          flippyFloppy += 1

    for f,e in zip(factors,errors):
      kStep = 1.0j/(f*kResolution)
      flippyFloppy = 0
      while (absError > e):
        if flippyFloppy == 4:
          break
        absError_prev = absError
        m = m + kStep
        trialBsca,trialBabs,_ = fastMie_SD(m,wavelength,dp,ndp)
        absError = error(Babs,trialBabs)
        if absError_prev-absError < 0:
          kStep *= -1
          flippyFloppy += 1

    scaError = error(Bsca,trialBsca)
    nStep = 0.001
    flippyFloppy = 0
    while (scaError > 0.01):
      if flippyFloppy == 4:
        break
      scaError_prev = scaError
      m += nStep
      trialBsca,trialBabs,_ = fastMie_SD(m,wavelength,dp,ndp)
      scaError = error(Bsca,trialBsca)
      if scaError_prev - scaError < 0:
        nStep *= -1
        flippyFloppy += 1

    absError = error(Babs,trialBabs)
    kStep = 0.0+0.00001j
    flippyFloppy = 0
    while absError > 0.005:
      if flippyFloppy == 4:
        break
      absError_prev = absError
      m += kStep
      trialBsca,trialBabs,_ = fastMie_SD(m,wavelength,dp,ndp)
      absError = error(Babs,trialBabs)
      if absError_prev - absError < 0:
        kStep *= -1
        flippyFloppy += 1

    scaError = error(Bsca,trialBsca)
    nStep = 0.0001
    flippyFloppy = 0
    while (scaError > 0.005):
      if flippyFloppy == 4:
        break
      scaError_prev = scaError
      m += nStep
      trialBsca,trialBabs,_ = fastMie_SD(m,wavelength,dp,ndp)
      scaError = error(Bsca,trialBsca)
      if scaError_prev - scaError < 0:
        nStep *= -1
        flippyFloppy += 1

    trialBsca,trialBabs,_ = fastMie_SD(m,wavelength,dp,ndp)
    scaError = error(Bsca,trialBsca)
    absError = error(Babs,trialBabs)
    resultM.append(m)
    resultScaErr.append(scaError)
    resultAbsErr.append(absError)

  if len(resultM)==1:
    return resultM[0], resultScaErr[0], resultAbsErr[0]
  else:
    return resultM, resultScaErr, resultAbsErr

def fastMie_SD(m, wavelength, dp, ndp):
#  http://pymiescatt.readthedocs.io/en/latest/inverse.html#fastMie_SD
  dp = coerceDType(dp)
  ndp = coerceDType(ndp)
  _length = np.size(dp)
  Q_sca = np.zeros(_length)
  Q_abs = np.zeros(_length)
  Q_back = np.zeros(_length)

  aSDn = np.pi*((dp/2)**2)*ndp*(1e-6)

  for i in range(_length):
    Q_sca[i],Q_abs[i],Q_back[i] = fastMieQ(m,wavelength,dp[i])

  Bsca = trapz(Q_sca*aSDn)
  Babs = trapz(Q_abs*aSDn)
  Bback = trapz(Q_back*aSDn)

  return Bsca, Babs, Bback

def fastMieQ(m, wavelength, diameter):
#  http://pymiescatt.readthedocs.io/en/latest/inverse.html#fastMieQ
  x = np.pi*diameter/wavelength
  if x==0:
    return 0, 0, 0
  elif x>0:
    nmax = np.round(2+x+4*(x**(1/3)))
    n = np.arange(1,nmax+1)
    n1 = 2*n+1
    x2 = x**2
    an,bn = Mie_ab(m,x)
    qext = (2/x2)*np.sum(n1*(an.real+bn.real))
    qsca = (2/x2)*np.sum(n1*(an.real**2+an.imag**2+bn.real**2+bn.imag**2))
    qback = (1/x2)*(np.abs(np.sum(n1*((-1)**n)*(an-bn)))**2)
    qabs = qext-qsca
    return qsca, qabs, qback
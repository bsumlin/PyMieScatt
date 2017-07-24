# -*- coding: utf-8 -*-

from PyMieScatt.Mie import Mie_ab
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from scipy.ndimage import zoom
from scipy.integrate import trapz
import warnings

def coerceDType(d):
  if type(d) is not np.ndarray:
    return np.array(d)
  else:
    return d

def GraphicalInversion(QscaMeasured,QabsMeasured,wavelength,diameter,nMin=1,nMax=3,kMin=0.00001,kMax=1,QbackMeasured=None,gridPoints=200,interpolationFactor=2,maxError=0.005,makePlot=True,axisOption=0,returnGraphElements=False,annotation=True):
  error = lambda measured,calculated: np.abs((calculated-measured)/measured)
  if (gridPoints*interpolationFactor) < 400:
    if QbackMeasured is not None:
      warnings.warn("If the product of gridPoints and interpolationFactor is less than 400, a unique solution may not be found!")
    else:
      warnings.warn("If the product of gridPoints and interpolationFactor is less than 400, not all solutions may be found, and errors may be high. Increase gridPoints to improve performance.")
  nRange = np.linspace(nMin,nMax,gridPoints)
  kRange = np.logspace(np.log10(kMin),np.log10(kMax),gridPoints)
  QscaList = []
  QabsList = []
  QbackList = []
  for n in nRange:
    s = []
    a = []
    b = []
    for k in kRange:
      m = n+k*1.0j
      Qsca,Qabs,Qback = fastMieQ(m,wavelength,diameter)
      s.append(Qsca)
      a.append(Qabs)
      b.append(Qback)
    QscaList.append(s)
    QabsList.append(a)
    QbackList.append(b)
  QscaList = zoom(np.transpose(np.array(QscaList)),interpolationFactor)
  QabsList = zoom(np.transpose(np.array(QabsList)),interpolationFactor)
  QbackList = zoom(np.transpose(np.array(QbackList)),interpolationFactor)

  n = zoom(nRange,interpolationFactor)
  k = zoom(kRange,interpolationFactor)
  if makePlot:
    plt.ion()
  else:
    plt.ioff()

  fig = plt.figure(figsize=(8.5,8.5))
  ax = fig.gca()

  scaLevels = np.array([QscaMeasured])
  absLevels = np.array([QabsMeasured])
  if QbackMeasured is not None:
    backLevels = np.array([QbackMeasured])

  scaChart = ax.contour(n,k,QscaList,scaLevels,origin='lower',linestyles='dashdot',linewidths=1.5,colors=('red'))
  absChart = ax.contour(n,k,QabsList,absLevels,origin='lower',linewidths=1.5,colors=('blue'))
  if QbackMeasured is not None:
    backChart = ax.contour(n,k,QbackList,backLevels,origin='lower',linestyles='dotted',linewidths=1.5,colors=('green'))

  distinctScaPaths = len(scaChart.collections[0].get_paths())
  distinctAbsPaths = len(absChart.collections[0].get_paths())
  if QbackMeasured is not None:
    distinctBackPaths = len(backChart.collections[0].get_paths())
  scaVertices = []
  absVertices = []
  if QbackMeasured is not None:
    backVertices = []
  for i in range(distinctScaPaths):
    scaVertices.append(scaChart.collections[0].get_paths()[i].vertices)
  for i in range(distinctAbsPaths):
    absVertices.append(absChart.collections[0].get_paths()[i].vertices)
  if QbackMeasured is not None:
    for i in range(distinctBackPaths):
      backVertices.append(backChart.collections[0].get_paths()[i].vertices)

  scaVertices = np.concatenate(scaVertices)
  absVertices = np.concatenate(absVertices)
  if QbackMeasured is not None:
    backVertices = np.concatenate(backVertices)

  nSolution,kSolution = find_intersections(scaVertices,absVertices)
  trials = [(x+(0+1j)*y) for x,y in zip(nSolution,kSolution)]
  if QbackMeasured is not None:
    oneTrueN,oneTrueK = find_intersections(backVertices,absVertices)
    _these = [np.round(x+(0+1j)*y,3) for x,y in zip(oneTrueN,oneTrueK)]
    _match = list(set(np.round(trials,3)).intersection(np.round(_these,3)))
    if len(_match) > 0:
      _matchIndex = list(np.round(trials,2)).index(np.round(_match[0],2))
      nSolution = nSolution[_matchIndex]
      kSolution = kSolution[_matchIndex]

  if type(nSolution)==np.float64:
    solutionSet = [nSolution + (0+1j)*kSolution]
  else:
    solutionSet = [(x+(0+1j)*y) for x,y in zip(nSolution,kSolution)]

  forwardCalculations = []
  for s in solutionSet:
    _s,_a,_ = fastMieQ(s,wavelength,diameter)
    forwardCalculations.append([_s,_a])
  solutionErrors = []
  for f in forwardCalculations:
    solutionErrors.append([error(f[0],QscaMeasured),error(f[1],QabsMeasured)])

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

  solutions = ax.scatter(nSolutionsToPlot,kSolutionsToPlot,marker='o',s=128,linewidth=1.5,edgecolor='k',facecolor='none',zorder=3)
  solutions = ax.scatter(nSolutionsToPlot,kSolutionsToPlot,marker='o',s=128,linewidth=0,edgecolor='none',facecolor='c',zorder=1,alpha=0.25)
  ax.set_xlabel('n',fontsize=16)
  ax.set_ylabel('k',fontsize=16)

  plt.title("Graphical Inversion for λ={w} nm, d={d} nm".format(w=wavelength,d=diameter),fontsize=18)

  ax.set_xlim((np.min(nRange),np.max(nRange)))
  ax.set_ylim((np.min(kRange),np.max(kRange)))
  ax.tick_params(which='both',direction='in')
  if annotation:
    solutionText = ''
    if len(solutionSet)==0:
      solutionText = "No solutions found!"
    else:
      for i,(s,f,e) in enumerate(zip(solutionSet,forwardCalculations,solutionErrors)):
        solutionText += "m{num}={n:1.3f}+{k:1.3f}i Esca={se:1.3f}% Eabs={ae:1.3f}%\n".format(num=i+1,n=s.real,k=s.imag,se=100*e[0],ae=100*e[1])
    st = AnchoredText(solutionText[:-1],loc=2)
    st.prop.set_size(14)
    st.patch.set_edgecolor('k')
    st.patch.set_alpha(0.85)
    ax.add_artist(st)

  if axisOption ==0:
    if max(kSolutionsToPlot) <= 0.5 or kMax <= 1:
      ax.set_yscale('log')
    else:
      ax.set_yscale('linear')
  elif axisOption == 1:
    ax.set_yscale('log')
  elif axisOption == 2:
    ax.set_xscale('log')
  elif axisOption == 3:
    ax.set_xscale('log')
    ax.set_yscale('log')
  else:
    pass

  lines = [scaChart.collections[0], absChart.collections[0]]
  labels = ["Qsca = {q:1.3f}".format(q=QscaMeasured),"Qabs = {q:1.3f}".format(q=QabsMeasured)]
  graphelements = [fig,ax,st,solutions,scaChart,absChart]
  if QbackMeasured is not None:
    lines.append(backChart.collections[0])
    labels.append("Qback = {q:1.3f}".format(q=QbackMeasured))
    graphelements.append(backChart)

  plt.legend(lines,labels,fontsize=16)
  if makePlot:
    plt.show(fig)
  else:
    plt.close(fig)

  if returnGraphElements:
    return solutionSet,forwardCalculations,solutionErrors, graphelements
  else:
    return solutionSet,forwardCalculations,solutionErrors

def GraphicalInversion_withSizeDistribution(BscaMeasured,BabsMeasured,wavelength,diameterBins,concentration,nMin=1,nMax=3,kMin=0.00001,kMax=1,gridPoints=60,interpolationFactor=2,maxError=0.005,axisOption=0,returnGraphElements=False,annotation=True,SDinset=False):
  error = lambda measured,calculated: np.abs((calculated-measured)/measured)
  diameterBins = coerceDType(diameterBins)
  concentration = coerceDType(concentration)
  nRange = np.linspace(nMin,nMax,gridPoints)
  kRange = np.logspace(np.log10(kMin),np.log10(kMax),gridPoints)
  solutionSet = []
  BscaList = []
  BabsList = []
  for n in nRange:
    s = []
    a = []
    for k in kRange:
      m = n+k*1.0j
      Bsca,Babs = fastMie_withSizeDistribution(m,wavelength,diameterBins,concentration)
      s.append(Bsca)
      a.append(Babs)
    BscaList.append(s)
    BabsList.append(a)
  BscaList = zoom(np.transpose(np.array(BscaList)),interpolationFactor)
  BabsList = zoom(np.transpose(np.array(BabsList)),interpolationFactor)

  n = zoom(nRange,interpolationFactor)
  k = zoom(kRange,interpolationFactor)

  fig = plt.figure(figsize=(8.5,8.5))
  ax = fig.gca()

  scaLevels = np.array([BscaMeasured])
  absLevels = np.array([BabsMeasured])

  scaChart = ax.contour(n,k,BscaList,scaLevels,origin='lower',linestyles='dashdot',linewidths=1.5,colors=('red'))
  absChart = ax.contour(n,k,BabsList,absLevels,origin='lower',linewidths=1.5,colors=('blue'))

  distinctScaPaths = len(scaChart.collections[0].get_paths())
  distinctAbsPaths = len(absChart.collections[0].get_paths())
  scaVertices = []
  absVertices = []
  for i in range(distinctScaPaths):
    scaVertices.append(scaChart.collections[0].get_paths()[i].vertices)
  for i in range(distinctAbsPaths):
    absVertices.append(absChart.collections[0].get_paths()[i].vertices)

  scaVertices = np.concatenate(scaVertices)
  absVertices = np.concatenate(absVertices)

  nSolution,kSolution = find_intersections(scaVertices,absVertices)
  solutionSet = [(x+(0+1j)*y) for x,y in zip(nSolution,kSolution)]

  forwardCalculations = []
  for s in solutionSet:
    _s,_a = fastMie_withSizeDistribution(s,wavelength,diameterBins,concentration)
    forwardCalculations.append([_s,_a])
  solutionErrors = []
  for f in forwardCalculations:
    solutionErrors.append([error(f[0],BscaMeasured),error(f[1],BabsMeasured)])

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

  solutions = ax.scatter(nSolution,kSolution,marker='o',s=128,linewidth=1.5,edgecolor='k',facecolor='none',zorder=3)
  solutions = ax.scatter(nSolution,kSolution,marker='o',s=128,linewidth=0,edgecolor='none',facecolor='c',zorder=1,alpha=0.25)

  ax.set_xlabel('n',fontsize=16)
  ax.set_ylabel('k',fontsize=16)

  plt.title("Graphical Solution for λ={w} nm".format(w=wavelength),fontsize=18)

  ax.set_xlim((np.min(nRange),np.max(nRange)))
  ax.set_ylim((np.min(kRange),np.max(kRange)))

  ax.tick_params(which='both',direction='in')

  if annotation:
    solutionText = ''
    if len(solutionSet)==0:
      solutionText = "No solutions found!"
    else:
      for i,(s,f,e) in enumerate(zip(solutionSet,forwardCalculations,solutionErrors)):
        solutionText += "m={n:1.3f}+{k:1.3f}i Esca={se:1.3f}% Eabs={ae:1.3f}%\n".format(n=s.real,k=s.imag,se=100*e[0],ae=100*e[1])
    st = AnchoredText(solutionText[:-1],loc=4)
    st.prop.set_size(16)
    st.patch.set_edgecolor('k')
    st.patch.set_alpha(0.85)
    ax.add_artist(st)

  if axisOption ==0:
    if max(kSolution) <= 0.5 or kMax <= 1:
      ax.set_yscale('log')
    else:
      ax.set_yscale('linear')
  elif axisOption == 1:
    ax.set_yscale('log')
  elif axisOption == 2:
    ax.set_xscale('log')
  elif axisOption == 3:
    ax.set_xscale('log')
    ax.set_yscale('log')
  else:
    pass

  lines = [scaChart.collections[0], absChart.collections[0]]
  labels = ["Bsca = {b:1.3f}".format(b=BscaMeasured),"Babs = {b:1.3f}".format(b=BabsMeasured)]
  plt.legend(lines,labels,fontsize=16)

  if(SDinset):
    left,bottom,width,height = [0.55,0.5,0.3,0.25]
    ax2=fig.add_axes([left,bottom,width,height],zorder=5)
    ax2.semilogx(diameterBins,concentration,color='k')
    ax2.set_xlabel("Diameter (nm)")
    ax2.set_ylabel("Concentration $\mathregular{(cm^{-3})}$")
    ax2.patch.set_alpha(0.85)
    ax2.tick_params(which='both',direction='in')
  if returnGraphElements:
    return solutionSet,forwardCalculations,solutionErrors,(fig,ax,st,scaChart,absChart,solutions)
  else:
    return solutionSet,forwardCalculations,solutionErrors

def find_intersections(A,B):
  arrayMinimum = lambda x1, x2: np.where(x1<x2, x1, x2)
  arrayMaximum = lambda x1, x2: np.where(x1>x2, x1, x2)
  arrayAll = lambda abools: np.dstack(abools).all(axis=2)
  slope = lambda line: (lambda d: d[:,1]/d[:,0])(np.diff(line, axis=0))

  x11, x21 = np.meshgrid(A[:-1, 0], B[:-1, 0])
  x12, x22 = np.meshgrid(A[1:, 0], B[1:, 0])
  y11, y21 = np.meshgrid(A[:-1, 1], B[:-1, 1])
  y12, y22 = np.meshgrid(A[1:, 1], B[1:, 1])

  m1, m2 = np.meshgrid(slope(A), slope(B))
  # Here we use masked arrays to properly treat the rare case where a line segment is perfectly vertical
  _m1 = np.ma.masked_array(m1,m1==-np.inf)
  _m2 = np.ma.masked_array(m2,m2==-np.inf)
  yi = (_m1*(x21-x11-y21/_m2)+y11)/(1-_m1/_m2)
  xi = (yi-y21)/_m2+x21

  xconds = (arrayMinimum(x11, x12) < xi, xi <= arrayMaximum(x11, x12),
            arrayMinimum(x21, x22) < xi, xi <= arrayMaximum(x21, x22) )
  yconds = (arrayMinimum(y11, y12) < yi, yi <= arrayMaximum(y11, y12),
            arrayMinimum(y21, y22) < yi, yi <= arrayMaximum(y21, y22) )

  return xi[arrayAll(xconds)], yi[arrayAll(yconds)]

def IterativeInversion(Qsca,Qabs,wavelength,diameter,tolerance=0.0005,initial_m=1.3+0.0j):
  error = lambda measured,calculated: np.abs((calculated-measured)/measured)
  if type(initial_m) not in [list,tuple]:
    initial_m = [initial_m]
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

def IterativeInversion_withSizeDistribution(Bsca,Babs,wavelength,diameterBins,sizeDistribution,tolerance=0.0005,initial_m=1.3+0.0j):
  error = lambda measured,calculated: np.abs((calculated-measured)/measured)
  if type(initial_m) not in [list,tuple]:
    initial_m = [initial_m]
  factors = [2.5, 5.0, 10.0, 25.0, 50.0]
  errors = [tolerance*x for x in [25,15,5,3,1]]
  resultM = []
  resultScaErr = []
  resultAbsErr = []

  for m in initial_m:
    nResolution = 100
    kResolution = 1000
    trialBsca,trialBabs = fastMie_withSizeDistribution(m,wavelength,diameterBins,sizeDistribution)

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
        trialBsca,trialBabs = fastMie_withSizeDistribution(m,wavelength,diameterBins,sizeDistribution)
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
        trialBsca,trialBabs = fastMie_withSizeDistribution(m,wavelength,diameterBins,sizeDistribution)
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
      trialBsca,trialBabs = fastMie_withSizeDistribution(m,wavelength,diameterBins,sizeDistribution)
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
      trialBsca,trialBabs = fastMie_withSizeDistribution(m,wavelength,diameterBins,sizeDistribution)
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
      trialBsca,trialBabs = fastMie_withSizeDistribution(m,wavelength,diameterBins,sizeDistribution)
      scaError = error(Bsca,trialBsca)
      if scaError_prev - scaError < 0:
        nStep *= -1
        flippyFloppy += 1

    trialBsca,trialBabs = fastMie_withSizeDistribution(m,wavelength,diameterBins,sizeDistribution)
    scaError = error(Bsca,trialBsca)
    absError = error(Babs,trialBabs)
    resultM.append(m)
    resultScaErr.append(scaError)
    resultAbsErr.append(absError)
    
  if len(resultM)==1:
    return resultM[0], resultScaErr[0], resultAbsErr[0]
  else:
    return resultM, resultScaErr, resultAbsErr

def fastMie_withSizeDistribution(m, wavelength, sizeDistributionDiameterBins, sizeDistribution):
  sizeDistributionDiameterBins = coerceDType(sizeDistributionDiameterBins)
  sizeDistribution = coerceDType(sizeDistribution)
  _length = np.size(sizeDistributionDiameterBins)
  Q_sca = np.zeros(_length)
  Q_abs = np.zeros(_length)

  aSDn = np.pi*((sizeDistributionDiameterBins/2)**2)*sizeDistribution*(1e-6)

  for i in range(_length):
    Q_sca[i],Q_abs[i],_ = fastMieQ(m,wavelength,sizeDistributionDiameterBins[i])

  Bsca = trapz(Q_sca*aSDn)
  Babs = trapz(Q_abs*aSDn)

  return Bsca, Babs

def fastMieQ(m, wavelength, diameter):
  x = np.pi*diameter/wavelength
  if x==0:
    return 0, 0
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


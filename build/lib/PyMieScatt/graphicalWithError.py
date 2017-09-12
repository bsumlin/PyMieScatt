import PyMieScatt as ps
q=ps.MieQ(1.77+0.66j,375,300)
ps.GraphicalInversion(q[1],q[2],375,300,axisOption=6)
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from PyMieScatt.Mie import Mie_ab
from scipy.ndimage import zoom
import warnings

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

def fastMieQ(m, wavelength, diameter):
#  http://pymiescatt.readthedocs.io/en/latest/inverse.html#fastMieQ
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

#def GraphicalInversion(QscaMeasured,QabsMeasured,wavelength,diameter,nMin=1,nMax=3,kMin=0.00001,kMax=1,QbackMeasured=None,gridPoints=200,interpolationFactor=2,maxError=0.005,makePlot=True,axisOption=0,returnGraphElements=False,annotation=True):
#  http://pymiescatt.readthedocs.io/en/latest/inverse.html#GraphicalInversion
plt.close('all')

#QscaMeasured = (1.315,0.02)
#QabsMeasured = (1.544, 0.008)
#QbackMeasured = (0.201,0.005)
q = ps.MieQ(1.77+0.63j,375,300)
QscaMeasured = (q[1],q[1]-1.3)
QabsMeasured = (q[2],q[2]-1.5)
QbackMeasured = (q[5],q[5]-0.2)

wavelength = 375
diameter = 300

# Default kwargs
nMin=1
nMax=3
kMin=0.00001
kMax=1
#QbackMeasured=None
gridPoints=200
interpolationFactor=2
maxError=0.001
makePlot=True
axisOption=4
returnGraphElements=False
annotation=True

error = lambda measured,calculated: np.abs((calculated-measured)/measured)

if type(QscaMeasured) in [list, tuple, np.ndarray]:
  scaError = QscaMeasured[1]
  QscaMeasured = QscaMeasured[0]
else:
  scaError = None
if type(QabsMeasured) in [list, tuple, np.ndarray]:
  absError = QabsMeasured[1]
  QabsMeasured = QabsMeasured[0]
else:
  absError = None
if type(QbackMeasured) in [list, tuple, np.ndarray]:
  backError = QbackMeasured[1]
  QbackMeasured = QbackMeasured[0]
else:
  backError = None
  
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
#%%
fig = plt.figure(figsize=(8.5,8.5))
ax = fig.gca()

scaLevels = np.array([QscaMeasured])
absLevels = np.array([QabsMeasured])

if QbackMeasured is not None:
  backLevels = np.array([QbackMeasured])
  if backError is not None:
    backErrorLevels = np.array([QbackMeasured+x for x in [-backError,backError]])

scaChart = ax.contour(n,k,QscaList,scaLevels,origin='lower',linestyles='dashdot',linewidths=1.5,colors=('red'))
absChart = ax.contour(n,k,QabsList,absLevels,origin='lower',linewidths=1.5,colors=('blue'))
if QbackMeasured is not None:
  backChart = ax.contour(n,k,QbackList,backLevels,origin='lower',linestyles='dotted',linewidths=1.5,colors=('green'))

if scaError is not None:
  scaErrorLevels = np.array([QscaMeasured+x for x in [-scaError, scaError]])
  absErrorLevels = np.array([QabsMeasured+x for x in [-absError, absError]])
  scaErrorChart = ax.contourf(n,k,QscaList,scaErrorLevels,origin='lower',colors=('red'),alpha=0.15)
  absErrorChart = ax.contourf(n,k,QabsList,absErrorLevels,origin='lower',colors=('blue'),alpha=0.15)
  ax.contour(n,k,QscaList,scaErrorLevels,origin='lower',linewidths=0.5,colors=('red'),alpha=0.5)
  ax.contour(n,k,QabsList,absErrorLevels,origin='lower',linewidths=0.5,colors=('blue'),alpha=0.5)
  if QbackMeasured is not None:
    backErrorLevels = np.array([QbackMeasured+x for x in [-backError, backError]])
    backErrorChart = ax.contourf(n,k,QbackList,backErrorLevels,origin='lower',colors=('green'),alpha=0.15)
    ax.contour(n,k,QbackList,backErrorLevels,origin='lower',linewidths=0.5,colors=('green'),alpha=0.5)

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
solutions = ax.scatter(nSolutionsToPlot,kSolutionsToPlot,marker='o',s=128,linewidth=0,edgecolor='none',facecolor='w',zorder=1,alpha=0.25)
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
#graphelements = [fig,ax,st,solutions,scaChart,absChart]
if QbackMeasured is not None:
  lines.append(backChart.collections[0])
  labels.append("Qback = {q:1.3f}".format(q=QbackMeasured))
#  graphelements.append(backChart)

plt.legend(lines,labels,fontsize=16)
if makePlot:
  plt.show(fig)
else:
  plt.close(fig)

#if returnGraphElements:
#  return solutionSet,forwardCalculations,solutionErrors, graphelements
#else:
#  return solutionSet,forwardCalculations,solutionErrors

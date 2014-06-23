import math
import numpy
import pylab
import random
import time
import steps.model as smodel
import steps.solver as solvmod
import steps.geom as stetmesh
import steps.rng as srng
import pylab


def gen_model():
	mdl = smodel.Model()
	A = smodel.Spec('A', mdl)
	vsys = smodel.Volsys('cytosolv', mdl)
	diff_A = smodel.Diff('diff_A', vsys, A, DCST)
	return mdl

## Initializing Model
mdl = smodel.Model() #Initialize model
molA = smodel.Spec('molA', mdl)
molB = smodel.Spec('molB', mdl)
molC = smodel.Spec('molC', mdl)
vsys = smodel.Volsys('vsys', mdl) #Volume
kreac_f = smodel.Reac('kreac_f', vsys, lhs=[molA, molB], rhs=[molC], kcst=0.3e6)
kreac_b = smodel.Reac('kreac_b', vsys, lhs=[molC], rhs=[molA, molB])
kreac_b.kcst = 0.7

## Initializing Geometry
wmgeom = swm.Geom()
comp = swm.Comp('comp', wmgeom)
comp.addVolsys('vsys')
comp.setVol(1.6667e-21)

##Random Number Generator
r = srng.create('mt19937', 256)
r.initialize(23412)

##Generating Solver (Well mixed)
sim = ssolver.Wmdirect(mdl, wmgeom, r)
sim.reset()

## Set Initial States
sim.setCompConc('comp', 'molA', 31.4e-6)
sim.setCompConc('comp', 'molB', 22.3e-6)

## Storing Variables
import numpy
tpnt = numpy.arange(0.0, 2.001, 0.001)
res = numpy.zeros([2001,3])

## Run Simulation
for t in range(0,2001):
	sim.run(tpnt[t])
	res[t,0] = sim.getCompCount('comp', 'molA')
	res[t,1] = sim.getCompCount('comp', 'molB')
	res[t,2] = sim.getCompCount('comp', 'molC')

## Plot Simulation
# Plot number of molecules of 'molA' over the time range:
pylab.plot(tpnt, res[:,0], label = 'A')
# Plot number of molecules of 'molB' over the time range:
pylab.plot(tpnt, res[:,1], label = 'B')
# Plot number of molecules of 'molC' over the time range:
pylab.plot(tpnt, res[:,2], label = 'C')
pylab.xlabel('Time (sec)')
pylab.ylabel('#molecules')
pylab.legend()
pylab.show()

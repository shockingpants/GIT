import math
import numpy
import pylab
import random
import time

import steps.model as smodel
import steps.solver as solvmod
import steps.geom as stetmesh
import steps.rng as srng

##################################
####	   Initial State
##################################
# The number of iterations to run
NITER = 10

# The data collection time increment (s)
DT = 0.001

# The simulation endtime (s)
INT = 0.101

# The number of molecules to be injected into the centre
NINJECT = 10000

# The number of tetrahedral elements to sample data from.
SAMPLE = 2000

# The diffusion constant for our diffusing species (m^2/s)
DCST = 20.0e-12


# Array to hold tetrahedron indices (integers)
tetidxs = numpy.zeros(SAMPLE, dtype = 'int')

# Array to hold tetrahedron radial distances (floats)
tetrads = numpy.zeros(SAMPLE)

##################################
####	  Methods
##################################


def gen_model():
	mdl = smodel.Model()
	A = smodel.Spec('A', mdl)
	vsys = smodel.Volsys('cytosolv', mdl)
	diff_A = smodel.Diff('diff_A', vsys, A, DCST)
	return mdl

def gen_geom():
	import steps.utilities.meshio as meshio
	mesh = smeshio.loadMesh('meshes/sphere_rad10_11Ktets')[0]

	# Find the total number of tetrahedrons in the mesh
	ntets = mesh.countTets()
	# Create a compartment object containing all tetrahedrons
	comp = stetmesh.TmComp('cyto', mesh, range(ntets))
	comp.addVolsys('cytosolv')

	# Fetch the central tetrahedron index and store:
	ctetidx = mesh.findTetByPoint([0.0, 0.0, 0.0])
	tetidxs[0] = ctetidx

	# Find the central tetrahedron's four neighbours:
	neighbidcs = mesh.getTetTetNeighb(ctetidx)
	tetidxs[1], tetidxs[2], tetidxs[3], tetidxs[4] = neighbidcs

	# Keep track how many tetrahedron we have stored so far
	stored = 5

	# Find the maximum and minimum coordinates of the mesh
	max = mesh.getBoundMax()
	min = mesh.getBoundMin()

	# Run a loop until we have stored as many indices as we require
	while (stored < SAMPLE):

		# Fetch 3 random numbers between 0 and 1:
		rnx = random.random()
		rny = random.random()
		rnz = random.random()

		# Find the related coordinates in the mesh:
		xcrd = min[0] + (max[0]-min[0])*rnx
		ycrd = min[1] + (max[1]-min[1])*rny
		zcrd = min[2] + (max[2]-min[2])*rnz

		# Find the tetrahedron that encompasses this point:
		tidx = mesh.findTetByPoint([xcrd, ycrd, zcrd])

		# -1 was returned if point is outside the mesh:
		if (tidx == -1): continue
		if (tidx not in tetidxs):
			tetidxs[stored] = tidx
			stored += 1

	# Find the barycenter of the central tetrahedron
	cbaryc = mesh.getTetBarycenter(ctetidx)

	for i in range(SAMPLE):
		# Fetch the barycenter of the tetrahedron:
		baryc = mesh.getTetBarycenter(tetidxs[i])

		# Find the radial distance of this tetrahedron:
		r = math.sqrt(math.pow((baryc[0]-cbaryc[0]),2) \
			+ math.pow((baryc[1]-cbaryc[1]),2) \
				+ math.pow((baryc[2]-cbaryc[2]),2))

		# Store the radial distance (in microns):
		tetrads[i] = r*1.0e6

	return mesh

model  = gen_model()
tmgeom = gen_geom()
rng = srng.create('mt19937', 512)
rng.initialize(2903)

sim = solvmod.Tetexact(model, tmgeom, rng)

tpnts = numpy.arange(0.0, INT, DT)
# Find how many "time points" we have
ntpnts = tpnts.shape[0]
# Create the data structure: iterations x time points x tet samples
res = numpy.zeros((NITER, ntpnts, SAMPLE))

# Fetch the index of the tetrahedron at the centre of the mesh
ctetidx = tmgeom.findTetByPoint([0.0, 0.0, 0.0])
# Run NITER number of iterations:
for i in range(NITER):
	sim.reset()
	# Inject all molecules into the central tet:
	sim.setTetCount(ctetidx, 'A', NINJECT)
	for j in range(ntpnts):
		sim.run(tpnts[j])
		# Loop over the tetrahedrons we are saving data for
		for k in range(SAMPLE):
			# Save the concentration in the tetrahedron, in uM
			res[i,j,k] = sim.getTetConc(tetidxs[k], 'A')*1.0e6

def plotres(tidx):
	# Check the time index is in range
	if (tidx >= INT/DT):
		print 'Time index is out of range.'
		return

	pylab.scatter(tetrads, res_mean[tidx], s=2)
	pylab.xlabel('Radial distance of tetrahedron ($\mu$m)')
	pylab.ylabel('Concentration in tetrahedron ($\mu$M)')
	t = tpnts[tidx]
	pylab.title('Unbounded diffusion. Time: ' + str(t) + 's')
	# Plot the analytical solution
	plotanlyt(t)
	pylab.xlim(0,0, 10.0)
	pylab.ylim(0.0)
	pylab.show()

def plotanlyt(t):
	segs = 100
	anlytconc = numpy.zeros((segs))
	radialds = numpy.zeros((segs))
	maxrad = 0.0
	for i in tetrads:
		if (i > maxrad): maxrad = i
	maxrad *= 1e-6
	intervals = maxrad/segs
	rad = 0.0
	for i in range((segs)):
		# Find the analytical concentration, and convert to mol/L
		anlytconc[i]=1.0e3*(1/6.022e23)*\
			((NINJECT/(math.pow((4*math.pi*DCST*t),1.5)))*\
			(math.exp((-1.0*(rad*rad))/(4*DCST*t))))
		radialds[i] = rad*1e6
		rad += intervals
	pylab.plot(radialds, anlytconc, color = 'red')

res_mean = numpy.mean(res, axis = 0)

import fipy as fp
import numpy as np
nx = 200
ny = 200
dx = 1.
dy = dx
L = dx * nx
mesh = fp.Grid2D(dx=dx, dy=dy, nx=nx, ny=ny)
phi = fp.CellVariable(name = "phi", mesh = mesh, value = 0.)
psi = fp.CellVariable(name = "psi", mesh = mesh, value = 0.)
D = 1.
eq = fp.TransientTerm() == fp.DiffusionTerm(coeff=D)
eq2 = fp.TransientTerm() == fp.DiffusionTerm(coeff=D*10)
#We apply Dirichlet boundary conditions
valueTopLeft = 0
valueBottomRight = 1
#to the top-left and bottom-right corners. Neumann boundary conditions are automatically applied to the top-right and bottom-left corners.

X, Y = mesh.faceCenters
facesTopLeft = ((mesh.facesLeft & (Y > L / 2)) | (mesh.facesTop & (X < L / 2)))
facesBottomRight = ((mesh.facesRight & (Y < L / 2)) | (mesh.facesBottom & (X > L / 2)))
x=mesh.cellCenters[0]
y=mesh.cellCenters[1]
# Setting value and getting value from a particular cell
#phi.setValue(1., where=((x < 55) & (x > 45 ) & (y < 55) & (y > 45 )))
#ind=fipy.tools.numerix.nearest(mesh.cellCenters.globalValue,np.array([[1.2],[1.2]]))
#pos=((50,51,52,53,54,54,56,57),(50,51,52,53,54,54,56,57))
#ID=mesh._getNearestCellID(pos)
ID=array([20100, 20066])
cellsnum=mesh.getNumberOfCells()
whre=np.zeros(cellsnum)
whre[ID]=1
conc=20.
phi.setValue(conc,where=whre)
psi.setValue(conc,where=whre)
#phi.constrain(valueTopLeft, facesTopLeft)
#phi.constrain(valueBottomRight, facesBottomRight)
#We create a viewer to see the results
if __name__ == '__main__':
	viewer = fp.Viewer(vars=phi, datamin=0., datamax=20)
	viewer2 = fp.Viewer(vars=psi, datamin=0., datamax=20)
	viewer.plot()
	viewer2.plot()
#and solve the equation by repeatedly looping in time:
timeStepDuration = 0.1 * dx**2 / (2 * D)
steps = 1000
for step in range(steps):
	eq.solve(var=phi, dt=timeStepDuration)
	eq2.solve(var=psi, dt=timeStepDuration)
	phi.setValue(phi+1, where=whre)
	if __name__ == '__main__':
		viewer.plot()
	if np.mod(step,10)==0:
		raw_input("Press key to continue")


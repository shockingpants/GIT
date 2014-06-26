import fipy as fp
nx = 100
ny = 100
dx = 1.
dy = dx
L = dx * nx
mesh = fp.Grid2D(dx=dx, dy=dy, nx=nx, ny=ny)
phi = fp.CellVariable(name = "solution variable", mesh = mesh, value = 0.)
D = 1.
eq = fp.TransientTerm() == fp.DiffusionTerm(coeff=D)
#We apply Dirichlet boundary conditions
valueTopLeft = 0
valueBottomRight = 1
#to the top-left and bottom-right corners. Neumann boundary conditions are automatically applied to the top-right and bottom-left corners.

X, Y = mesh.faceCenters
facesTopLeft = ((mesh.facesLeft & (Y > L / 2)) | (mesh.facesTop & (X < L / 2)))
facesBottomRight = ((mesh.facesRight & (Y < L / 2)) | (mesh.facesBottom & (X > L / 2)))
x=mesh.cellCenters[0]
y=mesh.cellCenters[1]
phi.setValue(1., where=((x < 55) & (x > 45 ) & (y < 55) & (y > 45 )))
#phi.constrain(valueTopLeft, facesTopLeft)
#phi.constrain(valueBottomRight, facesBottomRight)
#We create a viewer to see the results
if __name__ == '__main__':
	viewer = fp.Viewer(vars=phi, datamin=0., datamax=1.)
	viewer.plot()
#and solve the equation by repeatedly looping in time:
timeStepDuration = 10 * 0.9 * dx**2 / (2 * D)
steps = 1000
for step in range(steps):
	phi.setValue(phi.value*0.95, where=((x < 55) & (x > 45 ) & (y < 55) & (y > 45 )))
	eq.solve(var=phi, dt=timeStepDuration)
	if __name__ == '__main__':
		viewer.plot()
	if mod(step,10)==0:
		raw_input("Press key to continue")


from fipy import *
nx = 5
ny = 5
dx = 1.
dy = dx
L = dx * nx
mesh = Grid2D(dx=dx, dy=dy, nx=nx, ny=ny)
phi = CellVariable(name = "solution variable", mesh = mesh, value = 0.)
D = 1.
eq = TransientTerm() == DiffusionTerm(coeff=D)
#We apply Dirichlet boundary conditions
valueTopLeft = 0
valueBottomRight = 1
#to the top-left and bottom-right corners. Neumann boundary conditions are automatically applied to the top-right and bottom-left corners.

X, Y = mesh.faceCenters
facesTopLeft = ((mesh.facesLeft & (Y > L / 2)) | (mesh.facesTop & (X < L / 2)))
facesBottomRight = ((mesh.facesRight & (Y < L / 2)) | (mesh.facesBottom & (X > L / 2)))
#phi.constrain(valueTopLeft, facesTopLeft)
#phi.constrain(valueBottomRight, facesBottomRight)
#We create a viewer to see the results
if __name__ == '__main__':
	viewer = Viewer(vars=phi, datamin=0., datamax=1.)
	viewer.plot()
#and solve the equation by repeatedly looping in time:

timeStepDuration = 10 * 0.9 * dx**2 / (2 * D)
steps = 10
for step in range(steps):
	eq.solve(var=phi, dt=timeStepDuration)
	if __name__ == '__main__':
		viewer.plot()

print numerix.allclose(phi(((L,), (0,))), valueBottomRight, atol = 1e-2)
if __name__ == '__main__':
	raw_input("Implicit transient diffusion. Press <return> to proceed...")

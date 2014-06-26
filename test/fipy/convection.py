from fipy import *
L = 10.
nx = 5000
dx =  L / nx
mesh = Grid1D(dx=dx, nx=nx)
phi0 = 1.0
alpha = 1.0
phi = CellVariable(name=r"$\phi$", mesh=mesh, value=phi0)
solution = CellVariable(name=r"solution", mesh=mesh, value=phi0 * numerix.exp(-alpha * mesh.cellCenters[0]))

if __name__ == "__main__":
	viewer = Viewer(vars=(phi, solution))
	viewer.plot()
	raw_input("press key to continue")
raise Exception('Cuz i wanna end it')
phi.constrain(phi0, mesh.facesLeft)
## fake outflow condition
phi.faceGrad.constrain([0], mesh.facesRight)
eq = PowerLawConvectionTerm((1,)) + ImplicitSourceTerm(alpha)
eq.solve(phi)
print numerix.allclose(phi, phi0 * numerix.exp(-alpha * mesh.cellCenters[0]), atol=1e-3)

if __name__ == "__main__":
	viewer = Viewer(vars=(phi, solution))
	viewer.plot()
	raw_input("finished")    

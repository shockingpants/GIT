import scipy
import scipy.integrate
od=scipy.integrate.ode
import numpy as np
import os
import datafitting as df
import sys
import progressbar as pb
import __main__
import math

def test_func(t,y,param):
	##{{{
	"""
	example function to be input ODE
	The first arg must be t
	The second argument must be the variables
	The third argumanent are optional parameters
	"""
	y1,y2,y3=y
	p=param

	dy1dt=-p['a']*y1
	dy2dt=-p['b']*y2
	dy3dt=-p['c']*y2 + y1

	return np.array([dy1dt,dy2dt,dy3dt])
	##}}}

def test_osc(t,y,param=None):
	##{{{
	"""
	"""
	c,i=y
	
	p=param
	beta=p['beta']
	kappa=p['kappa']

	dcdt=c**2/(beta*i**2+c**2)-c
	didt=kappa*(c-i)
	return np.array([dcdt,didt])

	##}}}

def ODE(func,y0,t0=0,t1=50,dt=0.1,param=None,integrator='vode',int_method='Adams',**kwargs):
	##{{{
	"""
	Full ODE without Jacobian
	####------ARGUMENTS-----------
	func --> see test_func for an example
	y0 is the initial value for the various variables
	t0 is initial time point (only important if func depends on time)
	t1 is last time point
	dt is the size of time step
	param is the parameter used in func
	integrator is the integration type
	int_method is the int_method
	**kwargs are the optional arguments that scipy.integrate.ode uses
	####-------RETURNS------------
	time is a list of time points
	values contain
	(Use 'for i in values' to extract values for each individual variable)
	####---------------------------
	>>> import cellsim as ode
	>>> ode.ODE(ode.test_func,[1,1,3],param=[1,1,1])
	"""
	# Set up
	de=od(func).set_integrator(integrator,method=int_method,**kwargs) 
	de.set_initial_value(y0,t0).set_f_params(param) #Adds parameter
	
	# Run ODE
	time=[t0]
	values=np.array([y0])
	while de.successful() and de.t < t1:
		pbar = pb.ProgressBar(widgets=[pb.Percentage(), pb.Bar()], maxval=t1).start()
		pbar.update(de.t)
		de.integrate(de.t+dt)
		time.append(de.t) 
		values=np.concatenate((values,np.array([de.y])),axis=0)

	pbar.finish()
  	#values=zip(*values) #Note this may be a slow step	
	return time,values
	##}}}

class SDE_integrator(object):
	##{{{
	def __init__(self,func,vol=None,method="euler"):
		"""
		Applies Chemical Langevin Equation
		func is a list of functions for integration
		func should require 2 or 3 parameters --> func(t,y,param) where kwargs are the param
		method is the type of integration
		vol is the cell volume, used to modulate the noise
		"""
		self.func=func
		self.method=method
		self.success=True #Checks the success of integration. Based on rel and abs err
		self.param=None
		if vol is None: # Vol in L
			self.vol=3.14*1e-15  #1/(Na*V)  Vol of bacteria 3.14fL
		else:
			self.vol=vol
		self.volnorm=np.sqrt(1./6e23/self.vol) # Normalization of noise based on volume

	def set_initial_value(self,y0,t0):
		"""
		Sets the initial reactant amount
		Number of reactant
		"""
		self.y=y0
		self.t=t0
	
	def set_parameters(self,param):
		"""
		Sets param needed for the ODE
		param can be in any form that the function codes for
		"""
		self.param=param 
	
	def successful(self):
		"""
		Should return an error message if something wrong happened.
		"""
		return self.success
	
	def integrate(self,new_t):
		"""
		t is the current time. Important for time dependent SDE
		dt is the timestep
		"""
		dt=new_t-self.t
		if self.method=="rk":
			self.y = self.rk(self.t,dt)
			self.t = new_t
		elif self.method=="euler":
			self.y = self.euler(self.t,dt)
			self.t = new_t
		else:
			raise Exception("Unknown Method")
		## Check and negate negative values in an ad-hoc way. Need some ways to quantify it
		check = self.y < 0
		if check.any():
			self.y[check] *= -1

		return self.y #
	
	def evolve(self,t,y,param):
		"""
		Changes func into something that takes either param or not
		"""
		if param is None:
			return self.func(t,y)
		else:
			return self.func(t,y,param)
		

	def rk(self,t,dt):
		"""
		Applies Runge Kutta 
		http://arc.aiaa.org/doi/pdf/10.2514/3.56665
		current is t
		timestep is dt
		"""
		a21 =   2.71644396264860
		a31 = - 6.95653259006152
		a32 =   0.78313689457981
		a41 =   0.0
		a42 =   0.48257353309214
		a43 =   0.26171080165848
		a51 =   0.47012396888046
		a52 =   0.36597075368373
		a53 =   0.08906615686702
		a54 =   0.07483912056879

		q1 =   2.12709852335625
		q2 =   2.73245878238737
		q3 =  11.22760917474960
		q4 =  13.36199560336697
		
		assert "y" in vars(self)
		rv_n=np.random.randn(len(self.y))
		x1 = self.y
		k1 = dt * self.evolve(t,x1,self.param) + np.sqrt(dt) * self.volnorm * np.dot(np.sqrt(x1), rv_n)
		#Make sure output of func is np.array

		x2 = x1 + a21 * k1
		k2 = dt * self.evolve(t,x2,self.param) + np.sqrt(dt) * self.volnorm * np.dot(np.sqrt(x1), rv_n)

		x3 = x1 + a31 * k1 + a32 * k2
		k3 = dt * self.evolve(t,x3,self.param) + np.sqrt(dt) * self.volnorm * np.dot(np.sqrt(x1), rv_n)

		x4 = x1 + a41 * k1 + a42 * k2
		k4 = dt * self.evolve(t,x4,self.param) + np.sqrt(dt) * self.volnorm * np.dot(np.sqrt(x1), rv_n)

		return x1 + a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4
	
	def euler(self,t,dt):
		"""
		Applies Euler
		http://arc.aiaa.org/doi/pdf/10.2514/3.56665
		current is t
		timestep is dt
		"""
		x = self.y
		assert "y" in vars(self)
		rv_n=np.random.randn(len(self.y))
		return x + dt * self.evolve(t, x, self.param) + np.sqrt(dt) * self.volnorm * np.dot(np.sqrt(x), rv_n)
	##}}}
	
def SDE(func,y0,t0=0,t1=50,dt=0.1,param=None,vol=None,method='rk',**kwargs):
	##{{{
	"""
	####------ARGUMENTS-----------
	func --> see test_func for an example
	y0 is the initial value for the various variables
	t0 is initial time point (only important if func depends on time)
	t1 is last time point
	dt is the size of time step
	param is the parameter used in func
	integrator is the integration type
	int_method is the int_method
	**kwargs are the optional arguments that scipy.integrate.ode uses
	####-------RETURNS------------
	time is a list of time points
	values contain
	(Use 'for i in values' to extract values for each individual variable)
	"""
	# Set up
	de=SDE_integrator(func,method=method,vol=vol)
	de.set_initial_value(y0,t0) #Adds parameter
	de.set_parameters(param)
		
	# Run ODE
	time=[t0]
	values=np.array([y0])

	while de.successful() and de.t < t1:
		pbar = pb.ProgressBar(widgets=[pb.Percentage(), pb.Bar()], maxval=t1).start()
		pbar.update(de.t)
		de.integrate(de.t+dt)
		time.append(de.t)
		values=np.concatenate((values,np.array([de.y])),axis=0)

	pbar.finish()
  	#values=zip(*values) #Note this may be a slow step
	#values=[np.array(i) for i in values]
	return time,values
	##}}}

def example():
	##{{{
	"""
import cellsim
import cellsim.integrator as ode
param=dict([('a',1),('b',1),('c',2)])
time,values=ode.SDE(ode.test_func,[0.5,0.5,0.5],param=param,vol=1e-21,t0=0,t1=10,dt=0.01,method="euler")
plt.plot(time,values[:,0],label='euler')
time,values=ode.SDE(ode.test_func,[0.5,0.5,0.5],param=param,vol=1e-21,t0=0,t1=10,dt=0.01,method="rk")
plt.plot(time,values[:,0],label='rk')
time,values=ode.ODE(ode.test_func,[0.5,0.5,0.5],param=param,t0=0,t1=10,dt=0.01)
plt.plot(time,values[:,0],label='ODE')
plt.legend()
param=dict([("beta",4),("kappa",0.5)])
time,values=ode.SDE(ode.test_osc,[0.4,0.4],param=param,vol=1e-21,t0=0,t1=100,dt=0.01,method="euler")
plt.figure(2)
plt.plot(time,values[:,0],label='euler')
plt.figure(3)
plt.plot(values[:,0],values[:,1],label="euler")
time,values=ode.SDE(ode.test_osc,[0.4,0.4],param=param,vol=1e-21,t0=0,t1=100,dt=0.01,method="rk")
plt.figure(2)
plt.plot(time,values[:,0],label='rk')
plt.figure(3)
plt.plot(values[:,0],values[:,1],label="rk")
time,values=ode.ODE(ode.test_osc,[0.4,0.4],param=param,t0=0,t1=100,dt=0.01)
plt.figure(2)
plt.plot(time,values[:,0],label='ODE')
plt.legend()
plt.figure(3)
plt.plot(values[:,0],values[:,1],label="ODE")
plt.legend()
plt.show()
	"""
	pass
	##}}}

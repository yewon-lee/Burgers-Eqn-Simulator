import numpy as np
import matplotlib.pyplot as plt

nx = 41
ny = 41

dx = 2/float(nx-1)
dy = 2/float(ny-1)

u = numpy.ones((ny,nx))
v = numpy.ones((ny,nx))

sigma = .001
nu = 0.01
dt = sigma * dx * dy / u

def equation_of_motion(u,v,dx,dy,dt,nu):
	#generate the next state as a function of the old state
	un = u.copy()
	vn = v.copy()

	u[1:-1,1:-1] = (un[1:-1,1:-1] - dt/dx * un[1:-1,1:-1] * (un[1:-1,1:-1] - un[1:-1,0:-2]) - dt/dy * vn[1:-1,1:-1] * (un[1:-1,1:-1] - un[0:-2, 1:-1]) + nu*dt / dx**2 * (un[1:-1,2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) + nu * dt / dy**2 * (un[2:,1:-1] -2 * un[1:-1,1:-1] + un[0:-2,1:-1]))
	
	v[1:-1, 1:-1] = (vn[1:-1, 1:-1] - dt / dx * un[1:-1, 1:-1] * (vn[1:-1,1:-1] - vn[1:-1,0:-2]) - dt/dy * vn[1:-1, 1:-1] * (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) + nu * dt /dx**2 * (vn[1:-1,2:] - 2 * vn[1:-1,1:-1] + vn[1:-1,0:-2]) + nu * dt / dy**2 * (vn[2:,1:-1] - 2 * vn[1:-1,1:-1] + vn[0:-2, 1:-1])) 
	
	return (u,v)

def boundary(u, v, nozzle_u, nozzle_v, nx, ny, t_step):
	u[0,:] = 0
	u[-1,:] = 0
	u[:,0] = 0
	u[:,-1] = 0

	v[0,:] = 0
	v[-1,:] = 0
	v[:,0] = 0
	v[:,-1] = 0
	
	#special nozzle BC
	u[ny//2-2:ny//2+2,0] = nozzle_u[t_step]
	v[ny//2-2:ny//2+2,0] = nozzle_v[t_step]
	
	return (u,v)

def simulate(f, u, v, dx, dy, dt, nu, nozzle_u, nozzle_v, nx, ny,steps):
	for i in range(steps):
		(u,v) = equation_of_motion(u,v,dx,dy,dt,nu)
		(u,v) = boundary(u,v,nozzle_u,nozzle_v,nx,ny,i)
	return (u,v)

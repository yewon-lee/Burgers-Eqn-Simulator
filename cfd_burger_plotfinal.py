import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot, cm
from matplotlib.colors import Normalize
from cfd_burger import * 
import numpy as numpy

nx = 41
ny = 41
 
dx = 2 / float(nx - 1)
dy = 2 / float(ny - 1)

u = numpy.ones((ny, nx))
v = numpy.ones((ny, nx))

sigma = .001
nu = .01
dt = sigma * dx * dy / nu

nt = 2510 # the number of steps we're simulating

### assign initial conditions
final_u = numpy.zeros((nx, ny))
final_v = numpy.zeros((nx, ny))

### special BC for nozzle
#located at (0,1)
nozzle_u = numpy.append(10*numpy.ones(1000), numpy.zeros(nt))
nozzle_v = numpy.append(10*numpy.ones(1000), numpy.zeros(nt))

(final_u, final_v) = evolve(final_u, final_v, dt, dx, dy, nu, nozzle_u, nozzle_v, nx, ny, nt)

#plot
i = nt
ax = pyplot.figure()
norm = Normalize()
magnitude = numpy.sqrt(final_u[::2]**2 + final_v[::2]**2)
pyplot.quiver(final_u[::2],final_v[::2], norm(magnitude), scale = 60, cmap = pyplot.cm.jet)
ax.savefig('frame' + str(i).zfill(5)+'.png',dpi=300)
ax.clear()


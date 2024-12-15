import matplotlib.pyplot as plt
from scipy import integrate
from numpy import sin, sqrt, pi, cos
import numpy as np

ax = plt.figure().add_subplot(projection='3d')

# Make the grid
x, y, z = np.meshgrid(np.arange(-0.8, 5, 1),
                      np.arange(-0.8, 5, 1),
                      np.arange(-0.8, 5, 1))

d = float(input("The radius of the charged ring? : "))

# Define the integrand for X-component
def integrand_x(u, a, b, c, d):
    return (a - d * cos(u)) / np.sqrt((a**2 + b**2 + c**2 - 2 * a * d * cos(u) - 2 * b * d * sin(u) + d**2)**3)

# Define the integrand for Y-component
def integrand_y(v, a, b, c, d):
    return (b - d * cos(v)) / np.sqrt((a**2 + b**2 + c**2 - 2 * a * d * cos(v) - 2 * b * d * sin(v) + d**2)**3)

# Define the integrand for Z-component
def integrand_z(w, a, b, c, d):
    return c / np.sqrt((a**2 + b**2 + c**2 - 2 * a * d * cos(w) - 2 * b * d * sin(w) + d**2)**3)

# Constants
C = (8.85418782 * 10**-12) * 8 * (pi**2)
Q = float(input("Uniform charge of the ring = ", ))

# Arrays to store the electric field components
u = np.zeros(x.shape)
v = np.zeros(y.shape)
w = np.zeros(z.shape)

# Loop over each grid point and compute the integrals
for i in range(x.shape[0]):
    for j in range(x.shape[1]):
        for k in range(x.shape[2]):
            a, b, c = x[i, j, k], y[i, j, k], z[i, j, k]
            
            # Perform the integration for each component
            Ix = integrate.quad(integrand_x, 0, pi * 2, args=(a, b, c, d))[0]
            Iy = integrate.quad(integrand_y, 0, pi * 2, args=(a, b, c, d))[0]
            Iz = integrate.quad(integrand_z, 0, pi * 2, args=(a, b, c, d))[0]
            
            # Compute the field components at this point
            u[i, j, k] = Q * Ix / C
            v[i, j, k] = Q * Iy / C
            w[i, j, k] = Q * Iz / C

# Plot the vector field
ax.quiver(x, y, z, u, v, w, length=0.1, normalize=True)
plt.show()

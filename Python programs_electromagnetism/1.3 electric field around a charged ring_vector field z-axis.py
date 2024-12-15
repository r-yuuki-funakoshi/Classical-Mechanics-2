import matplotlib.pyplot as plt
from scipy import integrate
from numpy import sin, pi, cos
import numpy as np

ax = plt.figure().add_subplot(projection='3d')

# Define the z-axis points (x = 0, y = 0, z varies)
x = np.zeros(5)  # x = 0 along the z-axis
y = np.zeros(5)  # y = 0 along the z-axis
z = np.linspace(-0.8, 1, 5)  # Vary z coordinates

d = float(input("The radius of the charged ring? : "))

# Define the integrand for X-component
def integrand_x(x):
    return cos(x)

# Define the integrand for Y-component
def integrand_y(y):
    return sin(y)

# Define the integrand for Z-component
def integrand_z(z):
    return 1

# Constants
C = (8.85418782 * 10**-12) * 8 * (pi**2)
Q = float(input("Uniform charge of the ring = "))

Ix = integrate.quad(integrand_x, 0, pi * 2)
Iy = integrate.quad(integrand_y, 0, pi * 2)
Iz = integrate.quad(integrand_z, 0, pi * 2)

print("Integral answers")
print(Ix)
print(Iy)
print(Iz)

# Arrays to store the electric field components, but only for z-axis
u = np.zeros(5)  # No electric field along x-axis
v = np.zeros(5)  # No electric field along y-axis
w = Q * Iz[0] / C * np.ones(5)  # Field along z-axis (uniform)

# Plot the vector field on z-axis
ax.quiver(x, y, z, u, v, w, length=0.1, normalize=True)

plt.show()

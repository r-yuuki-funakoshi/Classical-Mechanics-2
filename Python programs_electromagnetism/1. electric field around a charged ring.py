from vpython import *
from scipy import integrate
from numpy import sin, sqrt, pi, cos
import numpy as np

op=0.3
visibility=True
indicator1=arrow(pos=vec(0,0,0), axis=vec(100,0,0), shaftwidth=0.1, visible=visibility, opacity=op, color=color.red)
indicator2=arrow(pos=vec(0,0,0), axis=vec(0,100,0), shaftwidth=0.1, visible=visibility, opacity=op, color=color.cyan)
indicator3=arrow(pos=vec(0,0,0), axis=vec(0,0,100), shaftwidth=0.1, visible=visibility, opacity=op, color=color.green)

# X COMPONENT----------------------------------------------------------------------------------

# Define the integrand function
def integrand(x, a, b, c, d):
    return (a - d * cos(x)) / np.sqrt((a**2 + b**2 + c**2 - 2 * a * d * cos(x) - 2 * b * d * sin(x) + d**2)**3)

# Convert inputs to float
a = float(input("x component of the position subject to the electric field caused by a ring = "))
b = float(input("y component of the position subject to the electric field caused by a ring = "))
c = float(input("z component of the position subject to the electric field caused by a ring = "))
d = float(input("The radius of the ring = "))

# Perform the integration
Ix = integrate.quad(integrand, 0, pi * 2, args=(a, b, c, d))

# Y COMPONENT----------------------------------------------------------------------------------

# Define the integrand function
def integrand(x, a, b, c, d):
    return (b - d * cos(x)) / np.sqrt((a**2 + b**2 + c**2 - 2 * a * d * cos(x) - 2 * b * d * sin(x) + d**2)**3)

# Perform the integration
Iy = integrate.quad(integrand, 0, pi * 2, args=(a, b, c, d))

# Z COMPONENT----------------------------------------------------------------------------------

# Define the integrand function
def integrand(x, a, b, c, d):
    return c / np.sqrt((a**2 + b**2 + c**2 - 2 * a * d * cos(x) - 2 * b * d * sin(x) + d**2)**3)

# Perform the integration
Iz = integrate.quad(integrand, 0, pi * 2, args=(a, b, c, d))

C = (8.85418782 * 10**-12) * 8 * (pi**2)
Q = float(input("Uniform charge of the ring (ex. 1E-12) = ",))

if a == 0 and b == 0:
    E = vec(0, 0, Iz[0]*Q / C)
else:
    E = vec(Ix[0]*Q / C, Iy[0]*Q / C, Iz[0]*Q / C)
print("Initial elecric field on the object", E)
print()
print("Integral answers: ")
print(Ix)
print(Iy)
print(Iz)

Ring = ring(pos = vec(0, 0, 0), radius = d, axis = vec(0, 0, 1), thickness = 0.05)
Object = sphere(pos = vec(a, b, c) , EF = E, radius = 0.2, mass = float(input("Mass of the object = ",)),
                charge = float(input("The charge of the object = ",)), v = vec(0, 0, 0), color = color.green, make_trail = True)
ElecField = arrow(pos = Object.pos, round = True, color = color.blue, axis = E / mag(E)) #NOT TO SCALE

for t in range (0, 400):
    dt = 1E-6
    t += dt
    rate(8)

# X COMPONENT----------------------------------------------------------------------------------
    def integrand(x, a, b, c, d):
        return (a - d * cos(x)) / np.sqrt((a**2 + b**2 + c**2 - 2 * a * d * cos(x) - 2 * b * d * sin(x) + d**2)**3)

# Convert inputs to float
    a = Object.pos.x
    b = Object.pos.y
    c = Object.pos.z
    d = d

    Ix = integrate.quad(integrand, 0, pi * 2, args=(a, b, c, d))

# Y COMPONENT----------------------------------------------------------------------------------

    def integrand(x, a, b, c, d):
        return (b - d * cos(x)) / np.sqrt((a**2 + b**2 + c**2 - 2 * a * d * cos(x) - 2 * b * d * sin(x) + d**2)**3)

    Iy = integrate.quad(integrand, 0, pi * 2, args=(a, b, c, d))

# Z COMPONENT----------------------------------------------------------------------------------

    def integrand(x, a, b, c, d):
        return c / np.sqrt((a**2 + b**2 + c**2 - 2 * a * d * cos(x) - 2 * b * d * sin(x) + d**2)**3)

    Iz = integrate.quad(integrand, 0, pi * 2, args=(a, b, c, d))

    if a == 0 and b == 0:
        E = Object.EF = vec(0, 0, Iz[0]*Q / C)
    else:
        E = Object.EF = vec(Ix[0]*Q / C, Iy[0]*Q / C, Iz[0]*Q / C)
    
    F = Object.charge * E
    Object.v += (F / Object.mass) * dt
    Object.pos += Object.v * dt
    ElecField.pos = Object.pos
    ElecField.axis = E / mag(E)

print("Final electric field on the object", E)
print("FInal force on the object", F)
print("Final velocity of the object", Object.v)

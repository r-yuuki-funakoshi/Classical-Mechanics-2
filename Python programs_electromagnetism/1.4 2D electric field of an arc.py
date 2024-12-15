from vpython import *
import matplotlib.pyplot as plt
from scipy import integrate
from numpy import sin, sqrt, pi, cos
import numpy as np

op=0.2
visibility=True
indicator1=arrow(pos=vec(0,0,0), axis=vec(100,0,0), shaftwidth=0.05, visible=visibility, opacity=op, color=color.red)
indicator2=arrow(pos=vec(0,0,0), axis=vec(0,100,0), shaftwidth=0.05, visible=visibility, opacity=op, color=color.cyan)
indicator3=arrow(pos=vec(0,0,0), axis=vec(0,0,100), shaftwidth=0.05, visible=visibility, opacity=op, color=color.green)

R = float(input("The radius of the charged arc? : "))
a = float(input("The angle of the arc in degrees with respect to the origin: ", ))

# Constants
C = 1 / (8.85418782 * 10**-12) * 4 * pi
Q = float(input("Uniform charge of the arc = ", ))

# Arrays to store the electric field components
E_x = C * Q * np.sin(pi * a / 180) / (a * R)
E_y = C * Q * (1 - np.cos(pi * a / 180)) / (a * R)
E = vec(E_x, E_y, 0)

arc = ring(radius = R, pos = vec(0, 0, 0), axis = vec(0, 0, 1))
ef = arrow(pos = vec(0, 0, 0), axis = E / mag(E), color = color.red)

print("Electric field at the origin: ", E)



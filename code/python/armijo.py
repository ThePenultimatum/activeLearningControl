import math
import numpy as np
import scipy.optimization as opt

# initial conditions
piOver2 = math.pi / 2
twoPi = 2 * math.pi

N = 100

T = 4

dt = float(T) / N

xd = 4
yd = 0
thetad = piOver2

x0 = 0
y0 = 0
theta0 = piOver2

x = np.array([[x0], [y0], [theta0]])

u1 = 1
u2 = -0.5

u = np.array([u1, u2])

# dynamics
def xDot(theta, u1):
	return math.cos(theta) * u1

def yDot(theta, u1):
	return math.sin(theta) * u1

def thetaDot(u2):
	return u2

def step(xVector, dt, uVector):
	u1val, u2val = uVector[0], uVector[1]
	thetaval = xVector[2][0]
	return (xVector + dt * np.array([[xDot(thetaval, u1)], [yDot(thetaval, u1)], [thetaDot(u2)]]))

# simulation
def simulate()
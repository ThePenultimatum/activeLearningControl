import math
import numpy as np
import scipy.optimize as opt

# initial conditions
piOver2 = math.pi / 2
twoPi = 2 * math.pi

N = 100

T = twoPi

dt = float(T) / N

def xdesired(t):
	return 4 * t / twoPi

def ydesired(t):
	return 0

def thetadesired(t):
	return piOver2

xdesiredMatrix = np.array([[xdesired(i * dt), ydesired(i * dt), thetadesired(i * dt)] for i in xrange(N+1)])

x0 = 0
y0 = 0
theta0 = piOver2

xinit = np.array([[x0], [y0], [theta0]])

u1init = 1
u2init = -0.5

uinit = np.array([u1init, u2init])#[[u1, u2] for i in xrange(N)])

# dynamics
def xDot(theta, u1val):
	return math.cos(theta) * u1val

def yDot(theta, u1val):
	return math.sin(theta) * u1val

def thetaDot(u2val):
	return u2val

def step(xVector, dt, uVector):
	u1val, u2val = uVector[0], uVector[1]
	thetaval = xVector[2][0]
	return (xVector + dt * np.array([[xDot(thetaval, u1)], [yDot(thetaval, u1)], [thetaDot(u2)]]))

# cost setup

Q = np.array([[1,0,0],
	          [0,1,0],
	          [0,0,1]]) # where each of these columns isolates weight for one of the x variables

R = np.array([[1,0],
	          [0,1]]) # where these are weights of the two control variables

P1 = np.array(([[10,  0],
	            [ 0, 10]]))

def getThreeMult(l, m, r):
	return np.matmul(np.matmul(l.transpose(), m), r)

def minJi(xlast, u, i):
    xdesiredcurr = xdesiredMatrix[i] # row vector
    xToUse = xcurr - xdesiredcurr
    

def j(xs, us):
	sumSoFar = 0
	xToUse = 0
	for i in xrange(N+1):
		xcurr = xs[i]
		xdesiredcurr = xdesiredMatrix[i] # row vector
		xToUse = xcurr - xdesiredcurr
		controlscurr = us[i]
		sumSoFar += (getThreeMult(xToUse, Q, xToUse)) + (getThreeMult(controlscurr, R, controlscurr))
	sumSoFar += np.matmul(np.matmul(xToUse.transpose(), P1), xToUse)
	return sumSoFar

# simulation
def simulate():
	print("Entering simulation")
	xs = np.array([[x0], [y0], [theta0]])
	us = np.array([[u1, u2]])
	# without optimization
	controlsAllNoOpt = np.array([[u1, u2] for i in xrange(N)])
	for i in xrange(N):

	return 1
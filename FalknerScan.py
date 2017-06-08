# Greg Hayrapetyan
# This program contains code for solving the Falkner-Skan equation using
# Runge-Kutta solver and uses the code to estimate the aerodynamic performance
# of a flat plate at an angle of attack, alpha, operating at a chord Reynolds
# number of 50,000.  The performance is also evaluated with addition of suction to
# extend the operating range beyond the stall angle.

import numpy as np
import matplotlib.pyplot as plt  # the tidy way
import scipy.optimize as opt

Re = 50000
N = 2000
betaArr = [1, 0, -0.1, -0.19884, -0.3]
deltastArr = [0.648, 1.217, 1.443, 2.3314, 2.3314]
#initGuess = [1.23259, 0.4696, 0.32, 0.0, 1]
tmax = 10.0
t = np.linspace(0, tmax, N)


# Falkner Scan
# f''' + f f'' + beta*(1-(f')^2) = 0
# y0 = f, y1 = f', y2 = f'', y2' = f''' = -y0 y2 + beta(1-y1^2)
# [y0, y1, y2, y3]' = [y1, y2, -y0*y2-beta(1-y1^2) ]
def F(t, y):
    rv = [y[1], y[2],  -y[0]*y[2]-beta*(1-y[1]*y[1])]
    return np.array(rv)

# Runge-Kutta Solver
def rkg( F, y0, t ):
    a = (np.sqrt(2)-1)/2
    b = (2 - np.sqrt(2))/2
    c = -np.sqrt(2)/2
    d = (2+np.sqrt(2))/2

    n = len( t )
    y = np.array( [ y0 ] * n )
    for i in xrange( n - 1 ):
        h = t[i+1] - t[i]
        k1 = h * F( t[i], y[i] )
        k2 = h * F( t[i] + 0.5 * h, y[i] + 0.5 * k1,  )
        k3 = h * F( t[i] + 0.5 * h, y[i] + a*k1 + b * k2 )
        k4 = h * F( t[i+1], y[i] + c*k2 + d*k3 )
        y[i+1] = y[i] + ( k1 + 2 * b * k2 + 2 * d * k3 + k4 ) / 6.0

    return y

def sfs(F,t,f0):
    numIter = 20
    initGuess = np.linspace(0, 1, numIter)
    initGuess[0] = 1.5
    initGuess[1] = 1.4
    curSoln = rkg(F, [f0,0,initGuess[0]], t)
    fpCur = curSoln[len(t)-1][1]
    #print fpCur
    for j in xrange( numIter - 2 ):
        j = j+1
        fpPrev = fpCur
        curSoln = rkg(F, [f0,0,initGuess[j]], t)
        fpCur = curSoln[len(t)-1][1]

        initGuess[j+1] = initGuess[j] - (fpCur-1) * (initGuess[j] - initGuess[j-1]) / (fpCur-fpPrev)

        if (np.absolute(fpCur-1) < 0.0001):
            print "f'' init: ", initGuess[j]
            return curSoln
    return curSoln

# Solve for various values of beta and plot the results
for i in xrange (len(betaArr)):
    beta = betaArr[i]
    deltast= deltastArr[i]
    tmax = 20.0*deltast
    t = np.linspace(0, tmax, N)
    soln = sfs(F, t, 0)

    plt.figure()
    plt.xlabel("$\eta$")
    #plt.xlim(0, 15)
    plt.ylim(0, 2)
    plt.ylabel("f")
    plt.plot(t, soln, linewidth=2.5, linestyle="-")

    #plt.plot(t, soln[1], linewidth=2.5, linestyle="-", label="f'")
    #plt.plot(t, soln[2], linewidth=2.5, linestyle="-", label="f''")
    plt.legend(loc='upper right')
    plt.show
    H = (tmax-soln[N-1][0])/soln[0][2]
    th = 0
    for j in xrange(N-2):
        th = th + soln[j+1][1]*(1-soln[j+1][1])*(t[j+1]-t[j])
    H = (tmax-soln[N-1][0])/th
    print "beta: ", beta, "delta*: ", (tmax-soln[N-1][0]), "theta*: ", th, "H:", H

import sys;
sys.exit(0)

# Calculate and plot Lift and Drag coefficients at various angles of attack

alphaSize = 30
alphaSt = 0.3123
alphaArr = np.linspace(0, alphaSt, alphaSize)
fps = np.linspace(0, alphaSt, alphaSize)
fss = np.linspace(0, alphaSt, alphaSize)
Fa = np.linspace(0, alphaSt, alphaSize)
Fn = np.linspace(0, alphaSt, alphaSize)
Lift = np.linspace(0, alphaSt, alphaSize)
Drag = np.linspace(0, alphaSt, alphaSize)
LiftSh = np.linspace(0, alphaSt, alphaSize)
DragSh = np.linspace(0, alphaSt, alphaSize)

for i in xrange (len(alphaArr)):
    beta = -2*alphaArr[i]/np.pi
    fss[i] = sfs(F, t,0)[0][2]
    beta = 2*alphaArr[i]/np.pi
    fps[i] = sfs(F,t,0)[0][2]


    Fa[i] = (4/np.sqrt(2))/np.sqrt(Re)*((np.sqrt(np.pi/(np.pi+alphaArr[i]))
    *((np.pi+alphaArr[i])/(np.pi-2*alphaArr[i]))*fss[i]) + (np.sqrt(np.pi/(np.pi-alphaArr[i]))
    *((np.pi-alphaArr[i])/(np.pi+2*alphaArr[i]))*fps[i]))

    Fn[i] = -(np.pi-alphaArr[i])/(np.pi+alphaArr[i]) + (np.pi+alphaArr[i])/(np.pi-alphaArr[i])

    Lift[i] = Fn[i]*np.cos(alphaArr[i]) - Fa[i]*np.sin(alphaArr[i])
    Drag[i] = Fn[i]*np.sin(alphaArr[i]) + Fa[i]*np.cos(alphaArr[i])

    LiftSh[i] = Fa[i]*np.sin(alphaArr[i])
    DragSh[i] = Fa[i]*np.cos(alphaArr[i])
    print "alpha: ", alphaArr[i], "C_L: ", Lift[i], "C_D: ", Drag[i]

    #if (np.isnan(fps[i]) or np.isnan(fss[i])):
    #    break

plt.figure()
plt.xlim(0,alphaSt)
plt.plot(alphaArr, Lift, linewidth=2.5, linestyle="-", label="$C_l$")
plt.plot(alphaArr, Drag, linewidth=2.5, linestyle="-", label="$C_d$")
plt.plot(alphaArr, DragSh, linewidth=2.5, linestyle="-", label="$C_d$ only shear")
plt.legend(loc='upper left')

plt.show()

plt.figure()
plt.xlim(0,alphaSt)
plt.plot(alphaArr, Lift/Drag, linewidth=2.5, linestyle="-", label="$C_l$/$C_d$")
plt.plot(alphaArr, LiftSh/DragSh, linewidth=2.5, linestyle="-", label="$C_l$/$C_d$ only shear")
plt.legend(loc='upper right')
plt.show()


alphaSize = 30
alphaSt = 0.45
alphaArr = np.linspace(0, alphaSt, alphaSize)
fps = np.linspace(0, alphaSt, alphaSize)
fss = np.linspace(0, alphaSt, alphaSize)

Fa = np.linspace(0, alphaSt, alphaSize)
Fn = np.linspace(0, alphaSt, alphaSize)
Lift = np.linspace(0, alphaSt, alphaSize)
Drag = np.linspace(0, alphaSt, alphaSize)
LiftSh = np.linspace(0, alphaSt, alphaSize)
DragSh = np.linspace(0, alphaSt, alphaSize)

for i in xrange (len(alphaArr)):
    beta = -2*alphaArr[i]/np.pi
    # nonzero in case of suction
    f0 = 0.001*np.sqrt(Re)*((2.0+beta)/(2.0-beta))/np.sqrt(2)
    fss[i] = sfs(F,t,f0)[0][2]

    beta = 2*alphaArr[i]/np.pi
    fps[i] = sfs(F,t,0)[0][2]


    Fa[i] = (4/np.sqrt(2))/np.sqrt(Re)*((np.sqrt(np.pi/(np.pi+alphaArr[i]))
    *((np.pi+alphaArr[i])/(np.pi-2*alphaArr[i]))*fss[i]) + (np.sqrt(np.pi/(np.pi-alphaArr[i]))
    *((np.pi-alphaArr[i])/(np.pi+2*alphaArr[i]))*fps[i]))

    Fn[i] = -(np.pi-alphaArr[i])/(np.pi+alphaArr[i]) + (np.pi+alphaArr[i])/(np.pi-alphaArr[i])

    Lift[i] = Fn[i]*np.cos(alphaArr[i]) - Fa[i]*np.sin(alphaArr[i])
    Drag[i] = Fn[i]*np.sin(alphaArr[i]) + Fa[i]*np.cos(alphaArr[i])

    LiftSh[i] = Fa[i]*np.sin(alphaArr[i])
    DragSh[i] = Fa[i]*np.cos(alphaArr[i])
    print "alpha: ", alphaArr[i], "C_L: ", Lift[i], "C_D: ", Drag[i]

    #if (np.isnan(fps[i]) or np.isnan(fss[i])):
    #    break

plt.figure()
plt.xlim(0,alphaSt)
plt.plot(alphaArr, Lift, linewidth=2.5, linestyle="-", label="$C_l$")
plt.plot(alphaArr, Drag, linewidth=2.5, linestyle="-", label="$C_d$")
plt.plot(alphaArr, DragSh, linewidth=2.5, linestyle="-", label="$C_d$ only shear")
plt.legend(loc='upper left')
#plt.plot(alphaArr, 6.28*alphaArr)
plt.show()

plt.figure()
plt.xlim(0,alphaSt)
plt.plot(alphaArr, Lift/Drag, linewidth=2.5, linestyle="-", label="$C_l$/$C_d$")
plt.plot(alphaArr, LiftSh/DragSh, linewidth=2.5, linestyle="-", label="$C_l$/$C_d$ only shear")
plt.legend(loc='upper right')
plt.show()

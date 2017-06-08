# Greg Hayrapetyan
# This program evaluates the NACA0010 airfoil performance for low Reynolds
# numbers (Re_c = 50,000) and incompressible flow (M=0).
# Thwaites method is used for cf calculation and the viscous cp data from XFLR is used to compute
# the normal and axial force coefficients on the airfoil and ultimately the cl and cd.
# Then cl and cd are predicted using both the inviscid and viscous cp distributions from XFLR from
# alpha=0 up to the stall angle of attack.
# Thwaites predictions are compared with the XFLR predictions.
 
import numpy as np
import matplotlib.pyplot as plt  # the tidy way

alpha, cl, cd = np.loadtxt("xflr.txt", usecols=(0,1,2), unpack=True)

clth = np.array(cl)
cdth = np.array(cd)
clth_i = np.array(cl)
cdth_i = np.array(cd)

ximp = []
cpi = []
cpv = []

xss = []
theta_ss = []
cf_ss = []
H_ss = []

xps = []
theta_ps = []
cf_ps = []
H_ps = []

for i in [0,1,2,3,4,5,6,7,8,9]:
    xt, cpit, cpvt = np.loadtxt("cp"+str(i)+"deg.txt", skiprows=6, usecols=(0,1,2), unpack=True)
    ximp.append(xt)
    cpi.append(cpit)
    cpv.append(cpvt)

for i in [0,3,6]:
    xss0, theta_ss0, cf_ss0, H_ss0 = np.loadtxt("cf"+str(i)+"degBD.txt", skiprows=7, usecols=(0, 7,3,1), unpack=True)
    xss.append(xss0)
    theta_ss.append(theta_ss0)
    cf_ss.append(cf_ss0)
    H_ss.append(H_ss0)

for i in [0,3,6]:
    xps0, theta_ps0, cf_ps0, H_ps0 = np.loadtxt("cf"+str(i)+"degBDps.txt", skiprows=7, usecols=(0, 7,3,1), unpack=True)
    xps.append(xps0)
    theta_ps.append(theta_ps0)
    cf_ps.append(cf_ps0)
    H_ps.append(H_ps0)

for i in [0,3,6]:
    plt.figure()
    plt.xlabel("x")
    plt.ylabel("angle: "+str(i))
    plt.plot(ximp[i],cpi[i], linewidth=2.5, linestyle="-", label="cpi")
    plt.plot(ximp[i],cpv[i], linewidth=2.5, linestyle="-", label="cpv")
    plt.legend(loc='upper right')
    plt.show

Re = 50000.0

def integ(a, b, n):
    s = 0.0
    s += np.power(1 - cp[0], 2.5)/2.0
    for i in range(1, n):
        s += np.power(1 - cp[i], 2.5)
    s += np.power(1 - cp[n], 2.5)/2.0
    return s * (b-a) / n

def integ_inviscid(a, b, n):
    s = 0.0
    s += np.power(1 - cp_i[0], 2.5)/2.0
    for i in range(1, n):
        s += np.power(1 - cp_i[i], 2.5)
    s += np.power(1 - cp_i[n], 2.5)/2.0
    return s * (b-a) / n

for j in [0,1,2,3,4,5,6,7,8,9]:

    cp = list(reversed(cpv[j][0:100]))
    cp_i = list(reversed(cpi[j][0:100]))

    x = list(reversed(ximp[j][0:100]))

    thetar = np.array(cp)
    thetasqss = np.array(cp)
    thetasqps = np.array(cp)
    l =  np.array(cp)
    cfss =  np.array(cp)
    cfps =  np.array(cp)
    Hss =  np.array(cp)
    Hps =  np.array(cp)

    if (j==0): st = 0
    if (j==1): st = 4
    if (j>1): st = 8

    # Viscous, SS
    fn = 0;
    ft = 0;
    for i in range(st,100):
        thetar[i] = np.sqrt(1-cp[i])
        thetasqss[i] = 0.45*np.power(thetar[i],-6)*integ(x[st], x[i], i )
        l[i] = 0.45*np.power(thetar[i],-6) * integ(x[st], x[i], i) * (thetar[i] - thetar[i-1]) / (x[i]-x[i-1])
        cfss[i] = 2*np.power(l[i] + 0.09, 0.62) / (thetar[i]*np.power(Re,0.5)*np.sqrt(thetasqss[i]))
        z = 0.25 - l[i]
        Hss[i] = 2.0 + 4.14*z -83.5*z*z+854*z*z*z-3337*z*z*z*z+4576*z*z*z*z*z
        if ((i>st) and (i < 98) and (l[i] > -0.09)):
            fn = fn - cp[i]*(x[i+1]-x[i])
            ft = ft + cfss[i]*(x[i+1]-x[i])
    # Inviscid, SS
    fn_i = 0;
    ft_i = 0;
    for i in range(st,100):
        thetar[i] = np.sqrt(1-cp_i[i])
        thetasqss[i] = 0.45*np.power(thetar[i],-6)*integ_inviscid(x[st], x[i], i )
        l[i] = 0.45*np.power(thetar[i],-6) * integ_inviscid(x[st], x[i], i) * (thetar[i] - thetar[i-1]) / (x[i]-x[i-1])
        cfss[i] = 2*np.power(l[i] + 0.09, 0.62) / (thetar[i]*np.power(Re,0.5)*np.sqrt(thetasqss[i]))
        z = 0.25 - l[i]
        Hss[i] = 2.0 + 4.14*z -83.5*z*z+854*z*z*z-3337*z*z*z*z+4576*z*z*z*z*z
        if ((i>st) and (i < 98) and (l[i] > -0.09)):
            fn_i = fn_i - cp_i[i]*(x[i+1]-x[i])
            ft_i = ft_i + cfss[i]*(x[i+1]-x[i])

    cp = cpv[j][99:199]
    cp_i = cpi[j][99:199]

    # Viscous, PS
    for i in range(st,100):
        thetar[i] = np.sqrt(1-cp[i])
        thetasqps[i] = 0.45*np.power(thetar[i],-6)*integ(x[st], x[i], i )
        l[i] = 0.45*np.power(thetar[i],-6) * integ(x[st], x[i], i) * (thetar[i] - thetar[i-1]) / (x[i]-x[i-1])
        cfps[i] = 2*np.power(l[i] + 0.09, 0.62) / (thetar[i]*np.power(Re,0.5)*np.sqrt(thetasqps[i]))
        z = 0.25 - l[i]
        Hps[i] = 2.0 + 4.14*z -83.5*z*z+854*z*z*z-3337*z*z*z*z+4576*z*z*z*z*z
        if ((i>st) and (i < 98) and (l[i] > -0.09)):
            fn = fn + cp[i]*(x[i+1]-x[i])
            ft = ft + cfps[i]*(x[i+1]-x[i])

    # Inviscid, PS
    for i in range(st,100):
        thetar[i] = np.sqrt(1-cp_i[i])
        thetasqps[i] = 0.45*np.power(thetar[i],-6)*integ_inviscid(x[st], x[i], i )
        l[i] = 0.45*np.power(thetar[i],-6) * integ_inviscid(x[st], x[i], i) * (thetar[i] - thetar[i-1]) / (x[i]-x[i-1])
        cfps[i] = 2*np.power(l[i] + 0.09, 0.62) / (thetar[i]*np.power(Re,0.5)*np.sqrt(thetasqps[i]))
        z = 0.25 - l[i]
        Hps[i] = 2.0 + 4.14*z -83.5*z*z+854*z*z*z-3337*z*z*z*z+4576*z*z*z*z*z
        if ((i>st) and (i < 98) and (l[i] > -0.09)):
            fn_i = fn_i + cp_i[i]*(x[i+1]-x[i])
            ft_i = ft_i + cfps[i]*(x[i+1]-x[i])

    clth[j] = fn*np.cos(j*np.pi/180.0) - ft*np.sin(j*np.pi/180.0)
    cdth[j] = ft*np.cos(j*np.pi/180.0) + fn*np.sin(j*np.pi/180.0)

    clth_i[j] = fn_i*np.cos(j*np.pi/180.0) - ft_i*np.sin(j*np.pi/180.0)
    cdth_i[j] = ft_i*np.cos(j*np.pi/180.0) + fn_i*np.sin(j*np.pi/180.0)

    if (j==0) or (j==3) or (j==6):
        plt.figure()
        plt.xlabel("x")
        plt.ylabel("angle: "+str(j))
        plt.xlim(x[st], 1)
        plt.ylim(0, 0.040)
        plt.plot(x,np.sqrt(thetasqss/Re), linewidth=2.5, linestyle="-", label="theta, ss, thwaites")
        plt.plot(x,np.sqrt(thetasqps/Re), linewidth=2.5, linestyle="-", label="theta, ps, thwaites")
        plt.plot(xss[j/3],theta_ss[j/3], linewidth=2.5, linestyle="-", label="theta, ss, xflr")
        plt.plot(xps[j/3],theta_ps[j/3], linewidth=2.5, linestyle="-", label="theta, ps, xflr")
        plt.legend(loc='upper right')

        plt.figure()
        plt.xlabel("x")
        plt.ylabel("angle: "+str(j))
        plt.xlim(x[st], 1)
        plt.ylim(0, 0.040)
        plt.plot(x,cfss, linewidth=2.5, linestyle="-", label="cf, ss, thwaites")
        plt.plot(x,cfps, linewidth=2.5, linestyle="-", label="cf, ps, thwaites")
        plt.plot(xss[j/3],cf_ss[j/3], linewidth=2.5, linestyle="-", label="cf, ss, xflr")
        plt.plot(xps[j/3],cf_ps[j/3], linewidth=2.5, linestyle="-", label="cf, ps, xflr")
        plt.legend(loc='upper right')

        # for problem 2
        if (j==3):
            cfps3 = cfps

        plt.figure()
        plt.xlabel("x")
        plt.ylabel("angle: "+str(j))
        plt.ylim(0, 5)
        plt.xlim(x[st], 1)
        plt.plot(x,Hss, linewidth=2.5, linestyle="-", label="H, ss, thwaites")
        plt.plot(x,Hps, linewidth=2.5, linestyle="-", label="H, ps, thwaites")
        plt.plot(xss[j/3],H_ss[j/3], linewidth=2.5, linestyle="-", label="H, ss, xflr")
        plt.plot(xps[j/3],H_ps[j/3], linewidth=2.5, linestyle="-", label="H, ps, xflr")
        plt.legend(loc='upper right')
        plt.show

plt.figure()
plt.xlabel("alpha")
plt.ylabel("cl")
plt.plot(alpha,cl, linewidth=2.5, linestyle="-", label="xflr")
plt.plot(alpha,2*np.pi*alpha*np.pi/180.0, linewidth=2.5, linestyle="-", label="thin airfoil")
plt.plot(alpha,clth, linewidth=2.5, linestyle="-", label="thwaites, viscous")
plt.plot(alpha,clth_i, linewidth=2.5, linestyle="-", label="thwaites, inviscid")
plt.legend(loc='upper left')
plt.show

plt.figure()
plt.xlabel("alpha")
plt.ylabel("cd")
plt.plot(alpha,cd, linewidth=2.5, linestyle="-", label="xflr")
plt.plot(alpha,cdth, linewidth=2.5, linestyle="-", label="thwaites, viscous")
plt.plot(alpha,cdth_i, linewidth=2.5, linestyle="-", label="thwaites, inviscid")
plt.legend(loc='upper left')
plt.show


# Problem 2
Pr = 0.713
rho = 1.284
mu = 1.725e-5
#thermal conductivity
k = 0.024
c_p = 1000
q = np.array(cp)
q2 = np.array(cp)
for i in range(0,99):
    Ch = cfps3[i]/(2*np.power(Pr, 2.0/3.0))
    if (x[i]>0.5) and (x[i]<=0.7):
        q[i] = Ch*c_p*rho*5*20
    if (x[i]>0.7) and (x[i]<=0.9):
        q[i] = Ch*c_p*rho*5*40
    if (x[i]>0.9):
        q[i] = 0

for i in range(0,99):
    Re = 5.0*rho*x[i]/mu
    if (x[i]>0.5) and (x[i]<=0.7):
        q2[i] = (0.332*k*np.power(Pr,1.0/3)*np.sqrt(Re))/(x[i]*np.power((1-np.power(0.5/x[i],3.0/4)),1.0/3))*20
    if (x[i]>0.7):
        q2[i] = (0.332*k*np.power(Pr,1.0/3)*np.sqrt(Re))/(x[i]*np.power((1-np.power(0.5/x[i],3.0/4)),1.0/3))*20
        q2[i] = q2[i] + (0.332*k*np.power(Pr,1.0/3)*np.sqrt(Re))/(x[i]*np.power((1-np.power(0.7/x[i],3.0/4)),1.0/3))*40


plt.figure()
plt.xlabel("x")
plt.ylabel("q")
plt.xlim(0.5,1)
plt.plot(x,q, linewidth=2.5, linestyle="-", label="Reynolds Analogy")
plt.plot(x,q2, linewidth=2.5, linestyle="-", label="Indicial Solns")
plt.legend(loc='upper left')
plt.show

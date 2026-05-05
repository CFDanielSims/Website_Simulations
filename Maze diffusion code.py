import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import math

# Prescribed parameters

untilsteady = True              #If False, runs until maxiter!!!
maxiter = 100000                #CHANGE FOR NUMBER OF ITERATIONS
Tthr = 1e-5                     #THRESHOLD VALUE, USED 1e-5 FOR QUESTION 2.3

alpha = 1.0e-3              
U0 = 1.0
dt = 0.005                      #TIME STEP, 0.01 IS STILL STABLE FOR N = 20, Q2.3
Lx = Ly = 1.0
Nlist = [20, 40, 60, 80]

T_left = 1.0
T_right = 0.0

Llist = len(Nlist)
Thc = [None] * Llist
Tvc = [None] * Llist
x = [None] * Llist
y = [None] * Llist
niter = [None] * Llist
ulist = [None] * Llist
vlist = [None] * Llist
Tlist = [None] * Llist
Xlist = [None] * Llist
Ylist = [None] * Llist
timelist = [None] * Llist
Tgradmaxlist = [None] * Llist

for res in np.arange(len(Nlist)):
    Tmaxgradplot = np.zeros(maxiter)

    Nx = Ny = Nlist[res]

    dx = Lx / Nx
    dy = Ly / Ny

    x[res] = np.linspace(0.5 * dx, Lx - 0.5 * dx, Nx)
    y[res] = np.linspace(0.5 * dy, Ly - 0.5 * dy, Ny)

    X, Y = np.meshgrid(x[res], y[res])

    u = -U0 * np.sin(np.pi * X) * np.cos(np.pi * Y)
    v = U0 * np.cos(np.pi * X) * np.sin(np.pi * Y)

    ubig = np.pad(u, pad_width=1, mode="constant", constant_values=0)
    vbig = np.pad(v, pad_width=1, mode="constant", constant_values=0)

    ubig[:, 0] = -ubig[:, 1]
    ubig[:, Nx+1] = -ubig[:, Nx]
    vbig[0, :] = -vbig[1, :]
    vbig[Ny+1, :] = -vbig[Ny, :]

    T = np.full((Ny + 2, Nx + 2), T_right)

    iter = np.arange(maxiter)

    Fox = alpha * dt / (dx**2)
    Foy = alpha * dt / (dy**2)
    print("Current resolution:", Nlist[res])
    print("Fox:", Fox, "Foy:", Foy)

    #Check for diffusion stability
    dtcriticaldiff = 1/(alpha*((1/dx**2)+(1/dy**2)))
    if dt > dtcriticaldiff:
        print("CODE STOPPED, Diffusion instability", "Fox:", Fox, "Foy:", Foy)
        sys.exit()

    Tmaxgrad = 0

    # Discretized temperature

    for it in iter:
        if it % 200 == 0:
            # print("Maximum gradient: ", Tmaxgrad)
            print("Number of iter's:", it)
        T_old = T.copy()

        for i in range(1, Nx+1):
            for j in range(1, Ny+1):

                uw = (ubig[j, i] + ubig[j, i-1]) / 2
                ue = (ubig[j, i+1] + ubig[j, i]) / 2
                vs = (vbig[j, i] + vbig[j-1, i]) / 2
                vn = (vbig[j+1, i] + vbig[j, i]) / 2

                Tw = T_old[j, i-1] if uw >= 0 else T_old[j, i]
                Te = T_old[j, i] if ue >= 0 else T_old[j, i+1]
                Ts = T_old[j-1, i] if vs >= 0 else T_old[j, i]
                Tn = T_old[j, i] if vn >= 0 else T_old[j+1, i]

                Cow = uw * dt / dx
                Coe = ue * dt / dx
                Cos = vs * dt / dy
                Con = vn * dt / dy
                
                #Check for advection stability
                dtcriticaladv = 1/((abs(ue)/dx)+(abs(uw)/dx)+(abs(vn)/dy)+(abs(vs)/dy))
                if dt > dtcriticaladv:
                    print("CODE STOPPED, CLS instability", "Cow:", Cow,"Coe:", Coe,"Cos:", Cos,"Con:", Con)
                    sys.exit()

                T[j, i] = (
                    (Cow * Tw - Coe * Te) +
                    (Cos * Ts - Con * Tn) +
                    Fox * (T_old[j, i-1] - 2*T_old[j, i] + T_old[j, i+1]) +
                    Foy * (T_old[j-1, i] - 2*T_old[j, i] + T_old[j+1, i]) +
                    T_old[j, i]
                )

        #Dirichlet BC (Left and right)
        T[:, 0] = 2*T_left - T[:, 1]
        T[:, Nx+1] = 2*T_right - T[:, Nx]

        #Neumann BC (Bottom and top)
        T[0, :] = T[1, :]
        T[Ny+1, :] = T[Ny, :]

        #STEADY STATE CRITERIA, IF NOT MET, GOES TO MAX ITER (VERY HIGH)

        Tdiff  = T - T_old
        Tmaxgrad = np.max(np.abs(Tdiff))
        Tmaxgradplot[it] = Tmaxgrad
        
        if untilsteady:
            if Tmaxgrad < Tthr:
                print("Steady stage reached at it", it)
                break

    Tlist[res] = T.copy()
    vlist[res] = vbig
    ulist[res] = ubig
    Xlist[res] = X
    Ylist[res] = Y
    Tgradmaxlist[res] = Tmaxgradplot.copy()

    timelist[res] = dt * it
    niter[res] = it
    umax = np.max([np.max(np.abs(ubig)), np.max(np.abs(vbig))])
    Comax = umax * dt / dx
    print("Fox:", Fox, "Foy: ", Foy,  "Comax: ", Comax)
    print("Finished resolution:", Nlist[res], "Number of iter's:", it)

    #Extraction of T(x, y = 0.5) and T(x = 0.5, y) Centerlines)
    Thc[res] = 0.5 * (T[int(Ny/2), 1:-1] + T[int((Ny/2) + 1), 1:-1]) #T(x, y = 0.5)
    Tvc[res] = 0.5 * (T[1:-1, int(Nx/2)] + T[1:-1, int((Nx/2)+1)]) #T(x = 0.5, y)

#END OF LOOP


def rootmeansquared(y_coarse, y_fine):

    N = len(y_coarse)
    M = len(y_fine)

    #x gives the real "physical" coordinates
    x_coarse = np.linspace(0, 1, N)
    x_fine = np.linspace(0, 1, M)

    #This interpolates the values at the "physical coordinates" of the fine data set to match the size of the coarse data set
    y_interpfine = np.interp(x_coarse, x_fine, y_fine)

    #Take the root mean squared error between two temperature fields based on temperature difference at maxiter
    rmse = np.sqrt(np.mean((y_coarse - y_interpfine)**2))

    return rmse

#Grid convergence calculation
rmsgridhor = np.zeros(Llist-1)
rmsgridver = np.zeros(Llist-1)

for i in np.arange(Llist-1):
    rmsgridhor[i] = rootmeansquared(Thc[i], Thc[i+1])
    rmsgridver[i] = rootmeansquared(Tvc[i], Tvc[i+1])

rmsgridhor_ratio = np.array(rmsgridhor)/T_left
rmsgridver_ratio = np.array(rmsgridver)/T_left

print("rmshor:", rmsgridhor, "rmsver:", rmsgridver)

n_plots = Llist + 5
n_cols = 4
n_rows = math.ceil(n_plots / n_cols)

fig, axs = plt.subplots(n_rows, n_cols, figsize=(15, 4*n_rows))
axs = axs.flatten()

# 2D temperature field

for i in range(Llist):
    im = axs[i].imshow(
        Tlist[i][1:-1, 1:-1],
        origin="lower",
        extent=[0.0, Lx, 0.0, Ly],
        cmap="inferno",
        aspect="equal",
    )

for i in range(Llist):
    axs[i].quiver(Xlist[i], Ylist[i], ulist[i][1:-1, 1:-1], vlist[i][1:-1, 1:-1], color="white")
    axs[i].set_title(f"Temperature field of N= {Nlist[i]}\n"
                     f"iter= {niter[i]}, time= {timelist[i]:.1f}s")
    axs[i].set_xlabel("x")
    axs[i].set_ylabel("y")

# colorbar linked to first plot
fig.colorbar(im, ax=axs[0])

#Horizontal centerline
for i in range(Llist):
    axs[Llist].plot(x[i], Thc[i], label = f'N={Nlist[i]}')
axs[Llist].set_title("Horizontal centerline")
axs[Llist].set_xlabel("x")
axs[Llist].set_ylabel("T")
axs[Llist].grid()
axs[Llist].legend()

#Vertical centerline
for i in range(Llist):
    axs[Llist + 1].plot(y[i], Tvc[i], label = f'N={Nlist[i]}')
axs[Llist + 1].set_title("Vertical centerline")
axs[Llist + 1].set_xlabel("y")
axs[Llist + 1].set_ylabel("T")
axs[Llist + 1].grid()
axs[Llist + 1].legend()

#RMS difference between coarse and finer grid
axs[Llist + 2].bar(Nlist[0:Llist-1], rmsgridhor_ratio, label = f'Hor rms(TNn+1 - TNn)/Thot')
axs[Llist + 2].set_xlabel("N")
axs[Llist + 2].set_ylabel("rms/Thot")
axs[Llist + 2].grid()
axs[Llist + 2].legend()

axs[Llist + 3].bar(Nlist[0:Llist-1], rmsgridver_ratio, label = f'Ver rms(TNn+1 - TNn)/Thot')
axs[Llist + 3].set_xlabel("N")
axs[Llist + 3].set_ylabel("rms/Thot")
axs[Llist + 3].grid()
axs[Llist + 3].legend()

for i in range(Llist):
    axs[Llist + 4].plot(np.arange(niter[i]), Tgradmaxlist[i][0:niter[i]], label = f'dTmax-It, N={Nlist[i]}')
axs[Llist + 4].set_xlabel("It")
axs[Llist + 4].set_ylabel("Maximum gradient per iteration")
axs[Llist + 4].grid()
axs[Llist + 4].legend()
axs[Llist + 4].set_yscale("log")
axs[Llist + 4].legend()

plt.tight_layout()
plt.show()

df1 = pd.DataFrame(T)
df2 = pd.DataFrame(ubig)

df1.to_excel(r"C:\Users\danie\OneDrive - Delft University of Technology\AA MSc EFPT\Q3\Computational Fluid Dynamics\Assignment 1\T.xlsx")
df2.to_excel(r"C:\Users\danie\OneDrive - Delft University of Technology\AA MSc EFPT\Q3\Computational Fluid Dynamics\Assignment 1\u.xlsx")
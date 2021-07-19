#*********************EXPERIMENT 5 PROJECT : Resistor Problem*******************
#*********************Done by Abhishek Sekar(EE18B067)**************************



#importing the necessary libraries
from pylab import *
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3


Niter=1500
Nx=int(input('What do you want Nx value to be?')) #taking Nx value from the user
Ny=Nx
Radius=Nx*8/25
phi=np.zeros((Ny,Nx))  #initialising phi matrix  
y=linspace(-0.5,0.5,Ny)
x=linspace(-0.5,0.5,Nx)
X,Y=np.meshgrid(x,y)# creates a coordinate system
ii = where(square(X) + square(Y) <= pow(0.35, 2))
phi[ii] = 1.0  # assigning V=1 for the circular region

#plotting the potential
plt.contourf(phi, cmap=cm.jet)#cmap gives the colour(cmap.jet=default) contourf fills contour lines with colour
plt.title("Figure 1 : Contour Plot of $\phi$")
ax = gca()
plt.colorbar( ax=ax, orientation='vertical')
plt.xlabel("x->")
plt.ylabel("y->")
plt.show()

# function to create Matrix for finding the Best fit using lstsq
# with no_of rows, columns by default  and vector x as arguments

def createAmatrix(nrow, x):
    A = zeros((nrow, 2))  # allocate space for A
    A[:, 0] = 1
    A[:, 1] = x
    return A


# function to find best fit errors using lstsq
def fitForError(errors, x):
    A = createAmatrix(len(errors), x)
    return A, np.linalg.lstsq(A, np.log(errors))[0]

# Function to compute function back from Matrix and Coefficients A and B
def computeErrorFit(M, c):
    return np.exp(dot(M,c))

errors = zeros(Niter)  # initialise error array to zeros
iterations = []  # array from 0 to Niter used for findind lstsq

for k in range(Niter):
    # copy the old phi
    oldphi = phi.copy()
    # Updating the potential
    phi[1:-1, 1:-1] = 0.25 *(phi[1:-1, 0:-2]+phi[1:-1, 2:]+phi[0:-2, 1:-1]+phi[2:, 1:-1])
    # applying boundary conditions: constant in the normal direction
    phi[:, 0] = phi[:, 1]  # Left edge
    phi[:, -1] = phi[:, -2]  # right edge
    phi[0, :] = phi[1, :]  # Top edge
    # Bottom edge is grounded so no boundary conditions
    # Assign 1 V to electrode region
    phi[ii] = 1.0
    # Appending errors for each iterations
    errors[k] = (abs(phi-oldphi)).max()
    iterations.append(k)

plt.semilogy(iterations[0::50], errors[0::50],'go', markersize=8, label="Original Error")
plt.semilogy(iterations, errors, markersize=8, label="Original Error cont",color='black')
plt.legend()
plt.title("Figure 2a : Error Vs No of iterations (Semilog)")
plt.xlabel("Niter ->")
plt.ylabel("Error ->")
plt.grid()
plt.show()

plt.loglog(iterations[0::50], errors[0::50],'go', markersize=8, label="Original Error")
plt.loglog(iterations, errors, markersize=8, label="Original Error",color='black')
plt.legend()
plt.title("Figure 2b : Error Vs No of iterations (Loglog)")
plt.xlabel("Niter ->")
plt.ylabel("Error ->")
plt.grid()
plt.show()

# to find the coefficients of fit1 and fit2
# M1 and M2 are matrices and c1 and c2 are coefficients
M1, c1 = fitForError(errors, iterations)  # fit1
M2, c2 = fitForError(errors[500:], iterations[500:])  # fit2

print("Fit1 : A = %g , B = %g" % ((np.around(np.exp(c1[0]),4), np.around(c1[1],4))))#%g is replaced by the corresponding terms
print("Fit2 : A = %g , B = %g" % ((np.around(np.exp(c2[0]),4), np.around(c2[1]))))#np.around rounds off the value

print("The time Constant (1/B) all iterations considered: %g" % (np.around(abs(1/c1[1]),4)))
print("The time Constant (1/B) for higher iterations (from 500) : %g" %(np.around(abs(1/c2[1]),4)))


# Calculating the fit using Matrix M and Coefficents C obained
error_fit1 = computeErrorFit(M1, c1)  # fit1
M2new = createAmatrix(len(errors), iterations)

# Error calculated for all iterations using coefficients found using lstsq
error_fit2 = computeErrorFit(M2new, c2)  # fit2 calculated

# Plotting the estimated error_fits using lstsq
# plotted for every 200 points for fit1 and fit2

plt.semilogy(iterations[0::200], error_fit1[0::200],'ro', markersize=8, label="Fit1")
plt.semilogy(iterations[0::200], error_fit2[0::200],'bo', markersize=6, label="Fit2")
plt.semilogy(iterations, errors, 'k', markersize=6, label="Actual Error")
plt.legend()
plt.title(r"Figure 3 : Error Vs No of iterations (Semilog)")
plt.xlabel("Niter  ->")
plt.ylabel("Error ->")
plt.grid()
plt.show()

#function calculating cumulative error
def cumerror(error, N, A, B):
    return -(A/B)*exp(B*(N+0.5))

#finds the stopping condition
def findStopCondn(errors, Niter, error_tol):
    cum_error = []
    for n in range(1, Niter):
        cum_error.append(cumerror(errors[n], n, exp(c1[0]), c1[1]))
        if(cum_error[n-1] <= error_tol):
            print("last per-iteration change in the error is %g"
                  % (np.around(cum_error[-1]-cum_error[-2],4)))
            return cum_error[n-1], n   #reduces the number of iterations based on tolerance



    print("last per-iteration change in the error is %g"
          % (np.abs(cum_error[-1]-cum_error[-2])))
    return cum_error[-1], Niter

#error tolerance
error_tol = pow(10, -8)
cum_error, Nstop = findStopCondn(errors, Niter, error_tol)
print("Stopping Condition N: %g and Error is %g" % (Nstop, cum_error))


fig1 = plt.figure()  # open a new figure
ax1= p3.Axes3D(fig1)  # Axes3D is the means to do a surface plot
plt.title("Figure 4: 3-D surface plot of the potential $\phi$")
surf = ax1.plot_surface(X, Y, phi, rstride=1, cstride=1, cmap=cm.jet)
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel('z')
cax = fig1.add_axes([1, 0, 0.1, 1])
fig1.colorbar(surf, cax=cax, orientation='vertical')
plt.show()


###### Contour Plot of the Potential:
plt.contourf(X, Y, phi, cmap=cm.jet)
plt.title("Figure 6 : Contour plot of Updated potential $\phi$")
ax = gca()
plt.colorbar(plt6, ax=ax, orientation='vertical')
plt.xlabel("x ->")
plt.ylabel("y ->")
plt.grid()
plt.show()
      
      
Jx = zeros((Ny, Nx))
Jy = zeros((Ny, Nx))
#from differential equation for current
Jx[1:-1, 1:-1] = 0.5*(phi[1:-1, 0:-2] - phi[1:-1, 2:])
Jy[1:-1, 1:-1] = 0.5*(phi[2:, 1:-1] - phi[0:-2, 1:-1])

      
# #### To Plot the current density using quiver, and mark the electrode via red dots :
plt.scatter(x[ii[0]], y[ii[1]], color='r', s=10, label="V = 1 region")#s is for size of the scatter 
plt.quiver(X, Y, Jx[::-1,:], Jy[::-1,:])
plt.xlabel('x ->')
plt.ylabel('y ->')
plt.legend()
plt.title("The Vector plot of the current flow")
plt.show()

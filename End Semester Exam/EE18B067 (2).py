#******************************END SEMESTER CODING ASSIGNMENT*******************
#******************************Done by Abhishek Sekar(EE18B067)*****************

#importing the necessary libraries
from pylab import *
import numpy as np
from numpy import *
import math
from matplotlib import cm #for  color map
import matplotlib.pyplot as plt


#Q: Create a function that solves the laplace equation and obtains phi(m,n) on the grid

#Defining all the arguments of the function
M=20 #no of nodes along x axis
N=40  #no of nodes along y axis, taken in this manner to ensure there's a node at the interface
delta=pow(10,-9) #error tolerance/desired accuracy
Niter=2500  # max no of iterations to complete
k=20 #value of height of the fluid
Er=2 #The relative permittivity or the dielectric constant of the fluid


# function to create Matrix for finding the Best fit using lstsq
# with no_of rows, 2 columns by default  and vector x as arguments

def A_matrix(nrow, x):
    A = np.zeros((nrow, 2))  # initializing the A matrix
    A[:, 0] = 1
    A[:, 1] = x
    return A

# function to find best fit for errors using lstsq
def Error_fitter_lstsq(errors, x):
    A = A_matrix(len(errors), x)
    return A, np.linalg.lstsq(A, np.log(errors))[0]

#function calculating cumulative error
def cum_error(N, A, B):
    return -(A/B)*exp(B*(N+0.5)) #the error is modeled in this form

#finds the stopping condition of the iterations based on tolerance value
def find_Stop_Condn(errors, Niter, error_tol,c1):
    cumerror = []
    for n in range(1, Niter):
        cumerror.append(cum_error(n, exp(c1[0]), c1[1])) #creates the cumulative error array based on error values plus directly takes values for the exp
        if(cumerror[n-1] <= error_tol):  #compares the cumulative error and the tolerance value
            return n   #deduces the number of iterations based on tolerance

def Laplace_solver(M,N,k,delta,Niter):
    y=linspace(0,20,N)
    x=linspace(0,10,M)
    phi=np.zeros((N,M))
    X,Y=np.meshgrid(x,y)# creates a coordinate system
    phi[39,:]=1#making the top of the tank have 1V   

    errors = zeros(Niter)  # initialise error array to zeros
    iterations = []  # array from 0 to Niter used for findind lstsq

    for i in range(Niter):
        #copy the old phi
        oldphi=phi.copy()
        #updating the full phi matrix
        phi[1:-1,1:-1]=0.25 *(phi[1:-1, 0:-2]+phi[1:-1, 2:]+phi[0:-2, 1:-1]+phi[2:, 1:-1])
        #selectively doing it for n=k alone so it overrides the previous data
        phi[k,1:-1]=(Er*phi[k-1,1:-1]+phi[k+1,1:-1])/(1+Er)
        #Appending errors for each iteration
        errors[i] = (abs(phi-oldphi)).max()
        iterations.append(i)
        
    #The corresponding least square fitting of the errors
    c = Error_fitter_lstsq(errors, iterations)[1]
        
    #error tolerance
    error_tol = delta
    Nstop = find_Stop_Condn(errors, Niter, error_tol,c)
    return Nstop,phi,errors

iterations=np.linspace(1,Niter,Niter) #an array containing the number of iterations
#initialization
errors=np.zeros(Niter)
phi=np.zeros((N,M))
Nstop=Laplace_solver(M,N,k,delta,Niter)[0]
y=linspace(0,20,N)
x=linspace(0,10,M)
X,Y=np.meshgrid(x,y)# creates a coordinate system or a grid of sorts
phi=Laplace_solver(M,N,k,delta,Niter)[1]
errors=Laplace_solver(M,N,k,delta,Niter)[2]
print("The number of iterations executed =",Nstop)
    
###### Contour Plot of the Potential:
plt.contourf(X, Y, phi, cmap=cm.jet)#cmap gives the colour(cmap.jet=default) contourf fills contour lines with colour
plt.title("Figure 1: Contour plot of Updated Potential $\phi$")
plt.colorbar()
plt.xlabel("x ->")
plt.ylabel("y ->")
plt.grid()
plt.show()

#plots of the errors: This will explain why we modelled the errors as an exponential
plt.semilogy(iterations[0::50], errors[0::50],'go', markersize=8, label="Original Error")#plots for every 50 values
plt.semilogy(iterations[0::50], errors[0::50], markersize=8, label="Original Error cont",color='black')
plt.legend()
plt.title("Figure 2a : Error Vs No of iterations (Semilog)")
plt.xlabel("Niter ->")
plt.ylabel("Error ->")
plt.grid()
plt.show()

plt.loglog(iterations[0::50], errors[0::50],'go', markersize=8, label="Original Error")
plt.loglog(iterations[0::50], errors[0::50], markersize=8, label="Original Error",color='black')
plt.legend()
plt.title("Figure 2b : Error Vs No of iterations (Loglog)")
plt.xlabel("Niter ->")
plt.ylabel("Error ->")
plt.grid()
plt.show()
#Question F
#For h/Ly=0.5 compute Ex and Ey for m+0.5 and n+0.5  and prove that Dn is continuous at m=k
#Now for a conductor like the tank,the volume charge density will be zero inside the tank.
#We could use this to our advantage to calculate the Ex and Ey at the centre of the mesh
def field_finder(phi,N,M):
    Ex = np.zeros((N,M))
    Ey = np.zeros((N,M))
    #finding the E field components using approximation of partial derivatives
    Ex[1:-1,1:-1]= 0.5*(phi[1:-1, 0:-2] - phi[1:-1, 2:])
    Ey[1:-1,1:-1]= 0.5*(phi[0:-2, 1:-1] - phi[2:, 1:-1])
    return Ex,Ey
 
#finding the electric fields as an average of its surrounding electric field components     
def centred_field_finder(Ex,Ey):
    Exmn=np.zeros((N,M))
    Eymn=np.zeros((N,M))
    Exmn=0.25*(Ex[:-1,:-1]+Ex[:-1,1:]+Ex[1:,:-1]+Ex[1:,1:])
    Eymn=0.25*(Ey[:-1,:-1]+Ey[:-1,1:]+Ey[1:,:-1]+Ey[1:,1:])
    return Exmn,Eymn

#Quiver Plot of the Electric Field
Ex,Ey=field_finder(phi,N,M)
Exmn,Eymn=centred_field_finder(Ex,Ey)
cs=plt.contour(2*X,2*Y,phi)  #to create the contour graph
plt.clabel(cs, cs.levels[:5], inline=1, fontsize=10) #this is to write values of the potential in the contour
plt.quiver(Exmn,Eymn,label="electric field")  #quiver plot with the Electric Field
plt.xlabel('x ->')
plt.ylabel('y ->')
plt.legend(loc='lower center')
plt.grid()
plt.title("Figure 4:The Vector plot of the Electric Field")
plt.show()


# Question E
#Run code for different values of h and plot Qtop and Qfluid vs h
# Use the expression Q=CV to find Qtop and Qfluid
# We can get the value of V from phi so all that's required is C
#using for loop to test this across multiple values of k
Qtop=zeros(9)
Qfluid=zeros(9)
iterations_1=[]
for j in range(2,20,2):#k starts from 2cm and advances to 18cm in steps of 2cm
    phi_1=Laplace_solver(M,N,2*j,delta,Niter)[1] #2j since everything is scaled by 2 as per the values of M and N
    Ex,Ey=field_finder(phi_1,N,M)
    Ex_1,Ey_1=centred_field_finder(Ex,Ey)
    Qtop[int(j/2-1)]=-(Ey_1[-2]).mean()
    Qfluid[int(j/2-1)]= Er*(-(Ex_1[:2*j,2]).mean()-(Ex_1[:2*j,-2]).mean()-(Ey_1[2:]).mean())
    iterations_1.append(j)
     
#Plotting Q vs h
plt.plot(iterations_1,Qtop,'go-',label="Qtop")
plt.plot(iterations_1,Qfluid,'bo-',label='Qfluid')
plt.title("Figure 3: The relationship of Q vs h")
plt.legend()
plt.xlabel("h->")
plt.ylabel("Q->")
plt.grid()
plt.show()

#now lets consider Dn right above and below the interface. For an ideal case the ratio between the two should be 1.
print("The ratio of Dn on either side is",np.around(Eymn[k].mean()/(Er*Eymn[k-1].mean()),4))
theta1=np.zeros(9)
theta2=np.zeros(9)
for i in range (1,9):   #because the angle is approx symmetric with respect to x, it is sufficient to do this
    theta1[i]=math.atan(Ey[k+1,i]/Ex[k+1,i])
    theta2[i]=math.atan(Ey[k-1,i]/Ex[k-1,i])
print("The value of theta 1 is",np.around(theta1.mean(),5)) #angle incidence
print("The value of theta 2 is",np.around(theta2.mean(),5)) #angle transmission


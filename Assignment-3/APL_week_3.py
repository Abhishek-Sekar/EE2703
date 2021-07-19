'''
*****************PROGRAM BY ABHISHEK SEKAR************************************
*******This Program fits the given data into a model and finds the error present in the model***********
BRIEF DESCRIPTION OF THE CODE:
The code performs the following functions:
    -> Plotting: This program makes use of matplotlib.pyplot for the various types of plots involved.
                 Here,four different plots are done in this program. A basic plot, a contour plot,an errorbar plot,
                 and a loglog plot. All these plots are labelled and annotated appropriately using legend() and annotate()
                 respectively.
    ->Calculative part:This program performs three major calculations. Firstly,the inner product of the matrices holding the
                        function value,time values and the coefficient values is verified if its equivalent to the function as
                        in theory.Secondly,the function plots mean square error to try approaching an estimate of the
                        coefficients.Lastly,the same thing is done using the least squares function and it is further used
                        in calculating the mean square error in coefficients tested.
                 

'''


#importing the necessary libraries

import numpy as np
import pandas as pd
import  scipy.special as sp
import scipy
import matplotlib.pyplot as plt
from numpy import *
from sys import exit

#defining the constants

Ao= 1.05
Bo = -0.105
sigma=logspace(-1, -3, 9)


# Function to generate the data(bessel function)

def g(t, A=1.05, B=-0.105):

    return(A*sp.jv(2, t)+B*t)


try:
        f=np.loadtxt('fitting.dat')
except Exception:

        print("fitting.dat not found!") #if the file is not found: display error message and exit
        exit()
       
       
df=pd.DataFrame(f,columns=['0','1','2','3','4','5','6','7','8','9']) #assigning a dataframe to deal with the given data



#Graph 1: Plotting all the columns of the given data

for i in range(1,10):
   
    plt.plot(df['0'],df[str(i)], label=r'$\sigma_'+str(i)+'$ = '+ str(around(sigma[i-1],3)))
    plt.legend(loc='best',fontsize=6.5)
    plt.grid()
plt.plot(df['0'],g(df['0']), label='True value',color='black')
plt.title('Data to be fitted to theory')
plt.legend(loc='best',fontsize=6.5)
plt.show()


#Graph 2: Plotting errorbars to compare the first column of the function with the true data

plt.plot(df['0'],g(df['0']),label='f(t)',color='black')
plt.errorbar(df.iloc[::5, 0], df.iloc[::5, 1], sigma[0], fmt="ro",label="Errorbar")
plt.xlabel("t ->")
plt.ylabel('f(t) with error ->')
plt.legend()


plt.annotate("True Data Curve",  (df.iloc[50, 0], g(df.iloc[50, 0])-0.1), xytext=(20, 48), textcoords="offset points", arrowprops={"arrowstyle": "->"})
plt.annotate("Noise present in data", (df.iloc[5, 0], df.iloc[5, 1]-0.06), xytext=(-20, -60), textcoords="offset points", arrowprops={"arrowstyle": "->"})
plt.title(r"Data points for $\sigma$ = 0.10 along with exact function")
plt.grid()
plt.show()



# part 3: Check if the inner product of these matrices are equal to the required function


M_matrix =c_[sp.jv(2,df['0']),df['0']]
p_matrix =c_[[Ao,Bo]]
val=dot(M_matrix,p_matrix)
flag = 0

for i in range(size(df.iloc[:,0])):
        flag += val[i]**3-g(df.iloc[i, 0])**3
if(flag == 0):
    print('The equality check has been completed and is successful')
else:
    print('The equality check has been completed and is unsuccessful')
   
   
A = linspace(0, 2, 21)     #using linspace function to create a series of 21 numbers for A and B in the given range
B = linspace(-0.2, 0, 21)  
e = zeros((len(A), len(B))) #initializing the error matrix


# Mean square error calculation to plot contours

for i in range(len(A)):
    for j in range(len(B)):

        e[i, j] = mean(square(g(df['0'], A[i], B[j])-g(df['0'])))
        #print(e[i,j])
       

# Plotting contours

cs = plt.contour(A, B, e, levels=20)
plt.xlabel("A ->")
plt.ylabel("B ->")
plt.title(r"Contour of $\epsilon_{ij}$")
plt.clabel(cs, cs.levels[:5], inline=1, fontsize=10) #to indicate values on the contour plots
plt.plot([1.05], [-0.105], 'ro')
#plt.grid()
plt.annotate("Exact Location", (1.05, -0.105),xytext=(-50, -30) ,textcoords="offset points")
plt.show()

# Solving for A,B using lstsq function

AB, *rest = scipy.linalg.lstsq(M_matrix, df.iloc[:, 1:10])
Aerr = np.array([square(AB[0, i]-Ao)for i in range(9)])
Berr = np.array([square(AB[1, i]-Bo)for i in range(9)])

#Error vs noise relationship

# Plotting the error in A and B vs stddev

plt.plot(sigma, Aerr, 'o--', linewidth=0.3, label='Aerr', dashes=(10, 10))
plt.plot(sigma, Berr, 'o--', linewidth=0.3, label='Berr', dashes=(10, 10))
plt.xlabel("Noise Standard Deviation ->")
plt.title("Variation of error with noise")
plt.ylabel("MSerror ->")
plt.legend()
plt.grid()
plt.show()

# Plotting the log log plot of error in A and B vs stddev

plt.loglog(sigma, Aerr, 'ro', label="Aerr")
plt.loglog(sigma, Berr, 'bo', label="Berr")
plt.legend()
plt.errorbar(sigma, Aerr, std(Aerr), fmt='ro')
plt.errorbar(sigma, Berr, std(Berr), fmt='bo')
plt.title("Variation of error with noise :loglog plot")
plt.xlabel(r'$\sigma_n$ ->')
plt.ylabel("MSerror ->")
plt.grid()
plt.show()

  
  
  

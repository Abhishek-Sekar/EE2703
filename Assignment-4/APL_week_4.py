'''
*****************PROGRAM BY ABHISHEK SEKAR************************************
*************This program finds the fourier coefficients********************** 
*************by two different methods and compares them **********************
BRIEF DESCRIPTION OF THE CODE:
    The code can be broken down into 7 modules,namely:
     1. Creating functions that will generate exp(x) and cos(cos(x)) and plotting the curves 
        and the expected curves which will be generated by fourier series.
        
     2. Creating functions that will generate the fourier series coefficients(by integration),thereby generating the 
        parent function itself and calculating the coefficients appropriately.
        
     3.Plotting magnitude spectrum of the fourier series coefficients and plotting the function generated
     
     4.Using the least squares method to calculate the fourier series and obtaining the parent function using those
       coefficients.
       
     5.Plotting magnitude spectrum of these coefficients obtained using least squares.
     
     6.Taking the absolute difference between the coefficients obtained by either methods
       and plotting graphs of the deviation to compare them. Also finds the maximum deviation present.
    
     7.Plotting the functions itself ,obtained using the least square fourier coefficients
'''

#importing the necessary libraries 

from numpy import *
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.integrate import quad
from sys import exit


#part1

def expo(x):         #function designed to calculate exp(x)

    return exp(x)

def coscos(x):      #function designed to calculate cos(cos(x))

    return np.cos(np.cos(x))

x = linspace(-2*np.pi, 4*np.pi, 400)   #linspace gives a linear distribution of numbers in the designated range

period=2*np.pi   #period of the function as per fourier series

#plot 1: semilog plot of exp(x) and the expected function from fourier series

plt.semilogy(x,expo(x),linewidth=3,label='$e^{x}$',color='orange')
plt.semilogy(x,expo(x%period),'--',label="expected function from fourier series $e^{x}$",color='black')
plt.xlabel("x ->")
plt.ylabel('$e^{x}$ ->')
plt.title('Expected function to be generated by fourier series ')
plt.grid()
plt.legend()
plt.show()

#plot 2: plot of cos(cos(x)) and the expected function from fourier series

plt.plot(x,coscos(x),linewidth=3,label='cos(cos(x))',color='yellow')
plt.plot(x,coscos(x%period),'--',label="expected function from fourier series",color='black')
plt.xlabel('x ->')
plt.ylabel('coscos(x) ->')
plt.title('Expected function to be generated by fourier series')
plt.grid()
plt.legend()
plt.axis([-5, 5, 0, 2])
plt.show()


#part2

def a_coeff(x,k,f):   
    return f(x)*np.cos(k*x)

def b_coeff(x,k,f):
    return f(x)*np.sin(k*x)

def coeff_calc(f):     #this function calculates the fourier coefficients and packages them in a single array
    coeff_val=[]
    
    coeff_val.append(scipy.integrate.quad(f, 0, 2*np.pi)[0]/(2*np.pi))
    
    for i in range(1,26):
        coeff_val.append((scipy.integrate.quad(a_coeff, 0, 2*np.pi, args=(i, f))[0])/np.pi) #giving function arguments inside quad
        coeff_val.append((scipy.integrate.quad(b_coeff, 0, 2*np.pi, args=(i, f))[0])/np.pi)
        
    return coeff_val

def A_matrix(x):     
    A = zeros((400,51))  # allocate space for A
    A[:,0]=1  # col 1 is all ones

    for k in range(1,26):
        A[:, 2*k-1] = np.cos(k*x)  # cos(kx) column
        A[:, 2*k] = np.sin(k*x)  # sin(kx) column

    return A

def fourier_function(c,num):  #this function obtains the parent function using the fourier series coefficients with them as the argument
    A=A_matrix(x)if(num==0)else A_matrix(x1)#num=1 for lst sq  function generation, done using conditional operator 
    b=dot(A,c)
    return b

#lists(arrays) store the various values present
exp_coeff=coeff_calc(expo)
coscos_coeff=coeff_calc(coscos)
exp_coeff_abs=np.abs(exp_coeff)   #np.abs provides the magnitude of the argument
coscos_coeff_abs=np.abs(coscos_coeff)
exp_fourier=fourier_function(exp_coeff,0)
coscos_fourier=fourier_function(coscos_coeff,0)

#part 3

#plot 3:  semilog plot of exp(x) obtained using fourier series coefficients
plt.semilogy(x, exp_fourier,'ro', label="Function obtained using Fourier Coefficients")
plt.ylim([pow(10, -1), pow(10, 4)])
plt.xlabel('x ->')
plt.ylabel('fourier $e^{x}$ ->')
plt.title('Fourier series function: $e^{x}$')
plt.legend()
plt.grid()
plt.show()

#plot4: plot of cos(cos(x)) obtained using fourier series coefficients
plt.plot(x, coscos_fourier,'ro',label="Function obtained using Fourier Coefficients")
plt.legend(loc='upper right')
plt.xlabel('x ->')
plt.ylabel('fourier cos(cos(x))->')
plt.title('Fourier series function: cos(cos(x))')
plt.axis([-5, 5, -0.5, 2])
plt.grid()
plt.show()


#Note: an and bn are plotted in the same graph for simplicity.

#plot5: semilog plot of the magnitude spectrum of exp(x) coefficients
plt.semilogy((exp_coeff_abs[1::2]),'ro', label="$a_{k}$s using Integration")   # By using array indexing methods we separate all odd indexes starting from 1 (ak)
plt.semilogy((exp_coeff_abs[2::2]),'bo', label="$b_{k}$s using Integration")   # and all even indexes starting from 2 (bk)
plt.legend()
plt.title("Fourier coefficients of $e^{x}$ (semi-log graph)")
plt.xlabel("k ->")
plt.ylabel("Magnitude of coeffients ->")
plt.grid()
plt.show()


#plot6: loglog plot of the magnitude spectrum of exp(x) coefficients
plt.loglog((exp_coeff_abs[1::2]),'ro',label="$a_{k}$s using Integration")
plt.loglog((exp_coeff_abs[2::2]),'bo',label="$b_{k}$s using Integration")
plt.legend(loc='upper right')
plt.title("Fourier coefficients of $e^{x}$ (Log-Log graph)")
plt.xlabel("k ->")
plt.grid()
plt.ylabel("Magnitude of coeffients ->")
plt.show()


#plot7: semilog plot of the magnitude spectrum of cos(cos(x)) coefficients
plt.semilogy((coscos_coeff_abs[1::2]),'ro',label="$a_{k}$s using Integration")
plt.semilogy((coscos_coeff_abs[2::2]),'bo',label="$b_{k}$s using Integration")
plt.legend(loc='upper right')
plt.title("Fourier coefficients of cos(cos(x)) (semi-log graph)")
plt.xlabel("k ->")
plt.grid()
plt.ylabel("Magnitude of coeffients ->")
plt.show()


#plot8: loglog plot of the magnitude spectrum of cos(cos(x)) coefficients
plt.loglog((coscos_coeff_abs[1::2]),'ro',label="$a_{k}$s using Integration")
plt.loglog((coscos_coeff_abs[2::2]),'bo',label="$b_{k}$s using Integration")
plt.legend(loc='upper right')
plt.title("Fourier coefficients of cos(cos(x)) (log-log graph)")
plt.xlabel("k ->")
plt.grid()
plt.ylabel("Magnitude of coeffients ->")
plt.show()


#part 4

def lstsq_coeff(f): #function to calculate the coefficients using the least squares method

    x1 = linspace(0,2*np.pi,401)
    x1 = x1[:-1]            # drop last term to have a proper periodic integral
    b = f(x1)
    A = A_matrix(x1)       #uses this function declared earlier to create the A_matrix
    c = scipy.linalg.lstsq(A, b)[0]  # the ’[0]’ is to pull out the  best fit vector coz lstsq returns a list with error as well.
    return c


# storing coefficients in respective vectors.
coeff_exp = lstsq_coeff(expo)
coeff_coscos = lstsq_coeff(coscos)

coeff_exp_abs = np.abs(coeff_exp)
coeff_coscos_abs= np.abs(coeff_coscos)

#part5

#plot9:semilog plot of the magnitude spectrum of exp(x) coefficients obtained using lstsq
plt.semilogy((coeff_exp_abs[1::2]),'go',label="$a_{k}$s using Least Squares")
plt.semilogy((coeff_exp_abs[2::2]),'bo',label="$b_{k}$s using Least Squares")
plt.grid()
plt.legend(loc='upper right')
plt.title("Fourier coefficients of $e^{x}$ (semilog graph)")
plt.xlabel("k ->")
plt.ylabel("Magnitude of coeffients ->")
plt.show()

#plot10:loglog plot of the magnitude spectrum of exp(x) coefficients obtained using lstsq
plt.loglog((coeff_exp_abs[1::2]),'go',label="$a_{k}$s using Least Squares ")
plt.loglog((coeff_exp_abs[2::2]),'bo',label="$b_{k}$s using Least Squares")
plt.grid()
plt.legend(loc='lower left')
plt.title("Fourier coefficients of $e^{x}$ (log-log graph)")
plt.xlabel("k ->")
plt.ylabel("Magnitude of coeffients ->")
plt.show()

#plot11:semilog plot of the magnitude spectrum of cos(cos(x)) coefficients obtained using lstsq
plt.semilogy((coeff_coscos_abs[1::2]),'go',label="$a_{k}$s using Least Squares")
plt.semilogy((coeff_coscos_abs[2::2]),'bo',label="$b_{k}$s using Least Squares")
plt.grid()
plt.legend(loc='upper right')
plt.title("Fourier coefficients of cos(cos(x)) (semilog graph)")
plt.xlabel("k ->")
plt.ylabel("Magnitude of coeffients ->")
plt.show()

#plot12:loglog plot of the magnitude spectrum of cos(cos(x)) coefficients obtained using lstsq
plt.loglog((coeff_coscos_abs[1::2]),'go',label="$a_{k}$s using Least Squares")
plt.loglog((coeff_coscos_abs[2::2]),'bo',label="$b_{k}$s using Least Squares")
plt.grid()
plt.legend(loc='best')
plt.title("Fourier coefficients of cos(cos(x)) (log-log graph)")
plt.xlabel("k ->")
plt.ylabel("Magnitude of coeffients ->")
plt.show()


#part 6

def coeff_compare(f): # Function to compare the coefficients got by integration and least squares and find largest deviation

    max_dev = 0
    deviations = np.abs(exp_coeff - coeff_exp) if(f==1) else np.abs(coscos_coeff - coeff_coscos) #conditional operator 1 : exp(x),2 : cos(cos(x)) 
    max_dev = np.amax(deviations) #finds maximum value of deviations

    return deviations, max_dev

dev1, maxdev1 = coeff_compare(1)
dev2, maxdev2 = coeff_compare(2)

print("Maximum deviation obtained in comparing coefficients of $e^{x}$ : ", around(maxdev1,4))
print("Maximum deviation obtained in comparing coefficients of cos(cos(x)) : ", around(maxdev2,4))


#plot13:Plotting the deviation vs n for exp(x)

plt.plot(dev1, 'g')
plt.title("Deviation between Coefficients for $e^{x}$")
plt.grid()
plt.xlabel("n ->")
plt.ylabel("Magnitude of Deviations ->")
plt.ylim([-1,3])
plt.show()


#plot 14: Plotting the deviation vs n for cos(cos(x))

plt.plot(dev2, 'g')
plt.title("Figure 8 : Deviation between coefficients for $\cos(\cos(x))$")
plt.grid()
plt.xlabel("n ->")
plt.ylabel("Magnitude of Deviations ->")
plt.show()


#part 7

# Define x1 from 0 to 2pi

x1 = linspace(0, 2*np.pi, 400)

exp_fn_lstsq = fourier_function(coeff_exp,1) #uses fourier_function defined earlier
coscos_fn_lstsq = fourier_function(coeff_coscos,1)



# comparing plots between lst square and the actual function
#plot 15 for exp(x)
plt.semilogy(x1, exp_fn_lstsq, 'go',label="Obtaining the function $e^{x}$ from Least Squares")
plt.semilogy(x1,expo(x1),label='$e^{x}$',color='black')
plt.legend()
plt.grid()
plt.title('$e^{x}$ obtained from lst sq coefficients')
plt.axis([0,2*np.pi,pow(10, -2), pow(10, 3)])
plt.show()

#plot 16 for cos(cos(x))
plt.plot(x1,coscos(x1),label='cos(cos(x))',color='black')
plt.legend()
plt.grid()
plt.title('cos(cos(x)) obtained from lst sq coefficients')
plt.axis([0, 2*np.pi,0.5, 1.3])
plt.show()
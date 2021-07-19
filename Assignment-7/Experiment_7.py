#*********Experiment7:Sympy*****************
#*********Abhishek Sekar EE18B067***********


import scipy 
import scipy.signal as sp
import matplotlib.pyplot as plt
import numpy as np
import sympy
from sympy import symbol,expand,simplify
from sympy.utilities.lambdify import lambdify

#lowpass filter function definition
def lowpass(R1=10000,R2=10000,C1=10**(-9),C2=10**(-9),G=1.586,Vi=1):#function to solve for the output voltage given the parameters
    s=sympy.symbols('s') #treats s as a symbol in the below equations
    A=sympy.Matrix([[0,0,1,-1/G],[-1/(1+s*R2*C2),1,0,0],[0,-G,G,1],[(-1/R1)-(1/R2)-s*C1,1/R2,0,s*C1]]) #matrix defines a matrix
    b=sympy.Matrix([0,0,0,-Vi/R1])
    V=A.inv()*b #solution is in s-space
    return A,b,V

#high pass filter function definition
def highpass(R1,R3,C1,C2,G,Vi):
    s=  sympy.symbols("s")
    A = Matrix([[0,-1,0,1/G],[s*C2*R3/(s*C2*R3+1),0,-1,0],[0,G,-G,1],[-s*C2-1/R1-s*C1,0,s*C2,1/R1]])
    b = Matrix([0,0,0,-Vi*s*C1])
    V = A.inv()*b
    return A,b,V
    
def laplace_converter(V):
    Vout=sympy.expand(sympy.simplify(V))#simplify simplifies the expression ie x^2+x becomes x(x+1)
    #expand :expands the simplified expression to a polynomial notation like x^2+x
    n,d=sympy.fraction(Vout) #writes separately as num and denom 
    n,d=sympy.Poly(n,s),sympy.Poly(d,s)# expresses n and d as polynomials of s
    num,denom=n.all_coeffs(),d.all_coeffs() #extracts the coefficients
    num,denom= [float(f) for f in num], [float(f) for f in denom]
    return num,denom


def inputfn(t):   #function for the input
    return (np.sin(2000*np.pi*t)+np.cos(2e6*np.pi*t))

def damped(t,freq=2e6,decay=500):  #default is for high freq function 
    return np.cos(freq*np.pi*t)*np.exp(-decay*t) * (t>0)
    

# The below code snippet will calculate and plot the transfer function for the lowpass filter.
s = symbols('s')
A,b,V = lowpass(10000,10000,1e-9,1e-9,1.586,1)
Vimp = V[3]#output is in third element of V (its in terms of s)
H = laplace_converter(Vimp)
ww = np.logspace(0,8,801)
ss = 1j*ww
hf = lambdify(s,Vo,'numpy')#like a short function definition for replacing s
vimp = hf(ss)

print('The transfer function of the lowpass filter is:',Vimp)

# The plot for Magnitude of transfer function of the lowpass filter.
V_bode=20*np.log10(np.abs(vimp))   #definition of bode
plt.semilogx(ww,V_bode,label='Impulse response')
plt.title('Fig1:Magnitude response:Impulse response of the lowpass filter')
plt.xlabel('freq ->')
plt.ylabel('|V_impulse|->')
plt.grid()
plt.legend(loc='best')
plt.show()

# The below code snippet will calculate and plot the step response for the lowpass filter.
A1,b1,V1 = lowpass(10000,10000,1e-9,1e-9,1.586,1/s)
Vstep = V1[3]
Hstep = laplace_converter(Vstep)
t,vstep = sp.impulse(Hstep,None,np.linspace(0,5e-3,10000))

# The plot for step response of the lowpass filter.
plt.plot(t,vstep,label='Step response')
plt.title("Fig2:Step Response of the lowpass filter")
plt.xlabel('t ->')
plt.ylabel('Vstep ->')
plt.grid()
plt.legend(loc='best')
plt.show()


# The response is also calculated and plotted for sum of sinusoids.
t,Vsum,svec = sp.lsim(H,inputfn(t),t)

# The plot for output response for sum of sinusoids of the lowpass filter.
plt.plot(t,Vsum,label='output response')
plt.title("Fig3:Output response for sum of solonoids of the lowpass filter")
plt.xlabel('t ->')
plt.ylabel('Vout ->')
plt.grid()
plt.legend(loc='best')
plt.show()

# The below snippet will calculate and the transfer function for the highpass filter.
A2,b2,V2 = highpass(10000,10000,1e-9,1e-9,1.586,1)
V_imp = V2[3]
H2 = laplace_converter(V_imp)
hf2 = lambdify(s,V_imp,'numpy')
v_imp = hf2(ss)

print('The transfer function of the highpass filter is:',V_imp)

# The plot for Magnitude of transfer function of the highpass filter.
V_bode2=20*np.log10(np.abs(v_imp))   #definition of bode
plt.semilogx(ww,V_bode2,label='Impulse response')
plt.title('Fig4:Magnitude response:Impulse response of the highpass filter')
plt.xlabel('freq ->')
plt.ylabel('|V_impulse|->')
plt.grid()
plt.legend(loc='best')
plt.show()

# The output response is calculated and plotted when the input is a damped sinusoid for both high and low frequency.
# High frequency.
t2 = np.linspace(0,1e-2,1e5) 
t2,Vhf,svec = sp.lsim(H2,damped(t2),t2)

# The plot for high frequency damped sinusoids.
plt.plot(t2,damped(t2),label='high frequency damped sinusoid')
plt.title("Fig5:plot of high frequency damped sinusoid")
plt.xlabel('t ->')
plt.ylabel('V ->')
plt.grid()
plt.legend(loc='best')
plt.show()

# The plot for high frequency damped sinusoid response from highpass filter.
plt.plot(t2,Vhf,label='high frequency damped sinusoid response')
plt.title("Fig7:plot of high frequency damped sinusoid response from highpass filter")
plt.xlabel('t ->')
plt.ylabel('Vout ->')
plt.grid()
plt.legend(loc='best')
plt.show()


# Low frequency.
t3 = np.linspace(0,1e-2,1e5)
t3,Vlf,svec = sp.lsim(H2,damped(t3,2e3),t3)

# The plot for low frequency damped sinusoids.
plt.plot(t3,damped(t3,2e3),label='low frequency damped sinusoid')
plt.title("Fig6:plot of low frequency damped sinusoid")
plt.xlabel('t ->')
plt.ylabel('V ->')
plt.grid()
plt.legend(loc='best')
plt.show()


# The plot for low frequency damped sinusoid response from  highpass filter.
plt.plot(t2,Vlf,label='low frequency damped sinusoid response')
plt.title("Fig8:plot of low frequency damped sinusoid response from highpass filter")
plt.xlabel('t ->')
plt.ylabel('Vout ->')
plt.grid()
plt.legend(loc='best')
plt.show()

# The step response is calculated for a highpass filter.
A3,b3,V3 = highpass(10000,10000,1e-9,1e-9,1.586,1/s)
V_step = V3[3]
H_step = laplace_converter(V_step)
t1,v_step = sp.impulse(H_step,None,np.linspace(0,5e-3,10000))

# The plot for step response of the highpass filter.
plt.plot(t1,v_step,label='Step response')
plt.title("Fig9:Step Response of the highpass filter")
plt.xlabel('t ->')
plt.ylabel('Vstep ->')
plt.grid()
plt.legend(loc='best')
plt.show()
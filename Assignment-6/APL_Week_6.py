#*****************EXPERIMENT 6: LAPLACE EQUATIONS************
#*****************Done by Abhishek Sekar(EE18B067)***********

#importing necessary header files
import scipy
import scipy.signal as sp
import matplotlib.pyplot as plt
import numpy as np


#obtaining the x signals using a function
def x_solver(decay,freq):
    num = np.poly1d([1, decay])
    denom = np.polymul([1, 0, 2.25], [1, 2*decay, (pow(freq, 2)+pow(decay, 2))])
    # Computes the impulse response of the transfer function
    X = sp.lti(num, denom)
    t, x = sp.impulse(X, None, np.linspace(0, 50, 100))
    return t,x

#obtaining the x signals for Q1 and Q2
t1,x1= x_solver(0.5,1.5)
t2,x2= x_solver(0.05,1.5)

#plotting the x signals
#plot for Q1
plt.plot(t1,x1,label='decay=0.5',color='black')
plt.title('Fig 1:Decay=0.5')
plt.xlabel('time ->')
plt.ylabel('x1(t)->')
plt.ylim((-1,1))
plt.grid()
plt.show()


#plot for Q2
plt.plot(t2,x2,label='decay=0.05',color='black')
plt.title('Fig 2:Decay=0.05')
plt.xlabel('time ->')
plt.ylabel('x2(t)->')
plt.grid()
plt.show()

#Q3
terms=linspace(1.4,1.6,4)
for i in (terms):
    t,x=x_solver(0.05,i)
    plt.plot(t,x)
    
    
plt.title('Fig3: Variation with frequency')
plt.xlabel('time ->')
plt.ylabel('$x_{i}$(t) ->')#subscript
plt.grid()
plt.show()


#Q4
#solving the coupled differential equations we get fourth order equations for x or y
#which reduce to the following transfer function
X=sp.lti([1,0,2],[1,0,3,0])
Y=sp.lti([2],[1,0,3,0])
t,x=sp.impulse(X,None,np.linspace(0,20,200))
t,y=sp.impulse(Y,None,np.linspace(0,20,200))

#plotting x and y

plt.plot(t,x,label='x(t)')
plt.legend()
plt.title('Fig4a:Solution to X in the coupled differential equations')
plt.xlabel('time ->')
plt.ylabel('x(t)->')
plt.grid()
plt.show()

plt.plot(t,y,label='y(t)')
plt.legend()
plt.title('Fig4b:Solution to Y in the coupled differential equations')
plt.xlabel('time ->')
plt.ylabel('y(t)->')
plt.grid()
plt.show()

#phase space y vs x plot
plt.plot(x,y)
plt.title('Fig4c:Phase Space plot of y vs x')
plt.xlabel('x ->')
plt.ylabel('y ->')
plt.grid()
plt.show()

#Q5 RLC circuit
H=sp.lti([1],[10**(-12),10**(-4),1])
w,S,phi=H.bode()  #gives bode parameters
plt.subplot(2,1,1)  
plt.semilogx(w,S)
plt.title('Fig5a: Bode Magnitude Plot')
plt.ylabel('|h(t)|->')
plt.xlabel('w')
plt.show()
plt.subplot(2,1,2) 
plt.semilogx(w,phi)
plt.title('Fig5b: Bode Phase Plot')
plt.ylabel('phase(h(t))->')
plt.xlabel('w')
plt.show()

#Q6
#part6
#define functioning returning input
def V_in(t):
    V_in=(np.cos(10**3*t)-np.cos(10**6*t))*(np.heaviside(t,1)) #heaviside(1) gives u(t)
    return V_in

#Short term
t1=linspace(0,30e-6,1000)
v1=V_in(t1)
t1,y,svec=sp.lsim(H,v1,t1)#convolves y=h*v1
plt.plot(t1,y)
plt.title('Fig6a:Response till 30 us')
plt.xlabel('time ->')
plt.ylabel('$V_{out}$(t) ->')
plt.grid()
plt.show()

#Long term
t2=linspace(0,1e-2,1000)
v2=V_in(t2)
t2,y,svec=sp.lsim(H,v2,t2)#convolves y=h*v2
plt.plot(t2,y)
plt.title('Fig6b:Response till 10ms')
plt.xlabel('time ->')
plt.ylabel('$V_{out}$(t) ->')
plt.grid()
plt.show()
'''
*****************PROGRAM BY ABHISHEK SEKAR************************************
*******This Program Solves R L C circuits with compatability for AC***********
BRIEF DESCRIPTION OF THE CODE:
    The code can be roughly split into 4 sections, namely:
        -> Class and function declarations:
            There is one user function called value_parser which has been defined inorder to 
            convert the alphanumeric value in the file to a float value.
            Classes replete with constructors have been defined for resistors , capacitors
            inductors , voltage source and current sources.
            
        -> Reading the file and taking values:
            Firstly, the given file is tested to see if its a netlist file or not,since
            I'd like the code to run only on valid netlist files.
            After this, testing for the ac nature of the circuit is done and values
            are then updated into objects of the classes appropriately.
        
        -> Updating the values in the matrices:
            This is by far the most important section in the code. Here the values are input to the matrix
            after performing modified nodal analysis. For a more efficient program ,I've reduced the equations
            involving the ground by performing a trick complying modified nodal analysis(remove the ground's row or column
            signatures of such elements where the ground is a node)
         
        -> Solving the equation and printing the result:    
            This part is self explanatory 

'''


#importing the necessary libraries
import numpy as np
from numpy import *
from sys import argv,exit
import math



#a function which converts the string for the value to the actual value
def value_parser(word):
    if('e' in word):
        val=float(word.split('e')[0])*pow(10,int(word.split('e')[1]))
        return val
    else:
        return int(word)   


#Declaring different classes for each circuit element for ease of handling
class resistor:
    def __init__(self,name,node1,node2,val):
        self.name=name
        self.val=value_parser(val)
        if(node1 =='GND'):
            self.node1=0
        else:
            self.node1=int(node1)    
        if(node2 =='GND'):
            self.node2=0
        else:
            self.node2=int(node2)


class inductor:
    def __init__(self,name,node1,node2,val,freq=0):
        self.name=name
        if(freq!=0):
            self.val=complex(0,value_parser(val)*freq*22/7)
        else:
            self.val=value_parser(val)
        if(node1 =='GND'):
            self.node1=0
        else:
            self.node1=int(node1)    
        if(node2 =='GND'):
            self.node2=0
        else:
            self.node2=int(node2)


class capacitor:
     def __init__(self,name,node1,node2,val,freq=0):
        self.name=name
        if(freq!=0):
            self.val=complex(0,(-(7/22)*(value_parser(val)*freq)**(-1)))
        else:
            self.val=value_parser(val)
        if(node1 =='GND'):
            self.node1=0
        else:
            self.node1=int(node1)    
        if(node2 =='GND'):
            self.node2=0
        else:
            self.node2=int(node2)


class voltage_inde:
    def __init__(self,name,node1,node2,char,val,phase=0):
        self.name=name
        if(char=='dc'):
            self.val=value_parser(val)
        elif(char=='ac'):
            #assuming Vpp is half the amplitude of the source
             #self.val=complex((value_parser(val)/2)*(math.cos(int(phase))),(value_parser(val)/2)*(math.sin(int(phase))))
             self.val=complex((value_parser(val))*(math.cos(int(phase))),(value_parser(val))*(math.sin(int(phase))))

        if(node1 =='GND'):
            self.node1=0
        else:
            self.node1=int(node1)    
        if(node2 =='GND'):
            self.node2=0
        else:
            self.node2=int(node2)


class current_inde:
    def __init__(self,name,node1,node2,char,val,phase=0):
        self.name=name
        if(char=='dc'):
            self.val=value_parser(val)
        elif(char=='ac'):
           #self.val=complex((value_parser(val)/2)*(math.cos(int(phase))),((value_parser(val)/2)*math.sin(int(phase))))
           self.val=complex((value_parser(val))*(math.cos(int(phase))),((value_parser(val))*math.sin(int(phase))))
        if(node1 =='GND'):
            self.node1=0
        else:
            self.node1=int(node1)    
        if(node2 =='GND'):
            self.node2=0
        else:
            self.node2=int(node2)
            
            
            
#declaring a dictionary for storing the lists of class objects
components={'resistor': [ ], 'capacitor':[ ], 'inductor':[ ],'voltage_source':[ ],'current_source':[ ]}
nodes=[]
     

file_start='.circuit'             #since a valid netlist file starts with .circuit. can be altered as per desire.
file_end='.end'
if(len(argv)==2):                 #We'd like to deal with only two arguments, one being the code's location and the other being the file location
    
    try:
        with open(argv[1]) as f:  
            lines=f.readlines()  
        start=-1
        end=-2
        flag=0
        for line in lines:
            if('.ac' in line): #to see if the circuit is in ac or dc mode
                flag=1             
                freq=value_parser(list(line.split("#")[0].split())[2])  #obtaining the freq of the circuit if it is ac mode
                                                  
            if(file_start==line[:len(file_start)]):
                start=lines.index(line)
            elif (file_end==line[:len(file_end)]):
                end=lines.index(line)
        if (start>=end):
            print("not a valid netlist")
            exit(0)
        for line in lines[start+1:end]:    
            vals=list(line.split("#")[0].split())
        
            if vals[1] not in nodes:              #creating a list which stores all the nodes 

                    nodes.append(vals[1])

            if(vals[2] not in nodes):
                nodes.append(vals[2])
                
              #obtaining data for all the circuit elements from the netlist file  
            if (line[0]=='R'):
                components['resistor'].append(resistor( vals[0], vals[1],vals[2],vals[3]))
            elif (line[0]=='L'):
                components['inductor'].append(inductor( vals[0], vals[1],vals[2],vals[3],freq))
            elif (line[0]=='C'):
                components['capacitor'].append(capacitor( vals[0], vals[1],vals[2],vals[3],freq))
            elif (line[0]=='V'):
                if(flag==1):
                    components['voltage_source'].append(voltage_inde( vals[0],vals[1],vals[2],vals[3],vals[4],vals[5]))
                else:
                    components['voltage_source'].append(voltage_inde( vals[0],vals[1],vals[2],vals[3],vals[4]))

            elif(line[0] =='I'):
                if(flag==1):
                    components['current_source'].append(current_inde( vals[0],vals[1],vals[2],vals[3],vals[4],vals[5]))
                else:
                    components['current_source'].append(current_inde( vals[0],vals[1],vals[2],vals[3],vals[4]))
   
        if(flag==0):                 #incorporating the dc behaviour of capacitors and inductors
           if(len(components['capacitor'])!=0):
               for C in components['capacitor']:
                   components['current_source'].append(current_inde( C.name,C.node1,C.node2,'dc',0))

           if(len(components['inductor'])!=0):
               for L in components['inductor']:
                   components['voltage_source'].append(voltage_inde( L.name,L.node1,L.node2,'dc',0))        
           
            #converting the nodes list to an integer list for ease of calculation
        for i in range(len(nodes)):
                        nodes[i]=i   
        
        #now we'd like to initialise and create the requisite matrices to solve the given circuit
        M_matrix=np.zeros(((len(nodes)+len(components['voltage_source'])-1),(len(nodes)+len(components['voltage_source'])-1)),dtype=complex)
        b_matrix=np.zeros(((len(nodes)+len(components['voltage_source']))-1),dtype=complex)
        x_matrix=np.zeros(((len(nodes)+len(components['voltage_source']))-1),dtype=complex)
        #writing the matrix or the modified nodal analysis equation for all the components
        
        
          #incorporates the ac characteristics of capacitor and inductor if the circuit is an ac one     
        if(flag==1):
            if(len(components['capacitor'])!=0):
                for C in components['capacitor']:
                    if (C.node1==0 or C.node2==0):
                        M_matrix[max(C.node1,C.node2)-1][max(C.node1,C.node2)-1] +=1/C.val
                    
                    else:
                        M_matrix[C.node1-1][C.node1-1] += 1/C.val
                        M_matrix[C.node1-1][C.node2-1] -= 1/C.val
                        M_matrix[C.node2-1][C.node1-1] -= 1/C.val   
                        M_matrix[C.node2-1][C.node2-1] += 1/C.val 
            if(len(components['inductor'])!=0):
                for L in components['inductor']:
                    if (L.node1==0 or L.node2==0):
                        M_matrix[max(L.node1,L.node2)-1][max(L.node1,L.node2)-1] +=1/L.val
                    
                    else:
                        M_matrix[L.node1-1][L.node1-1] += 1/L.val
                        M_matrix[L.node1-1][L.node2-1] -= 1/L.val
                        M_matrix[L.node2-1][L.node1-1] -= 1/L.val   
                        M_matrix[L.node2-1][L.node2-1] += 1/L.val
                        
        
        if(len(components['resistor'])!=0):
            for R in components['resistor']:
                
                if (R.node1==0 or R.node2==0):
                    M_matrix[max(R.node1,R.node2)-1][max(R.node1,R.node2)-1] +=1/R.val
                    
                else:
                    M_matrix[R.node1-1][R.node1-1] += 1/R.val
                    M_matrix[R.node1-1][R.node2-1] -= 1/R.val
                    M_matrix[R.node2-1][R.node1-1] -= 1/R.val   
                    M_matrix[R.node2-1][R.node2-1] += 1/R.val 
        
        
        
        if(len(components['voltage_source'])!=0):            
            for V in components['voltage_source']:
                if (V.node2!=0):
                    M_matrix[V.node2-1][len(nodes)+components['voltage_source'].index(V)-1] -=1
                    M_matrix[len(nodes)+components['voltage_source'].index(V)-1][V.node2-1] -=1
                if(V.node1!=0):
                    M_matrix[V.node1-1][len(nodes)+components['voltage_source'].index(V)-1] +=1
                    M_matrix[len(nodes)+components['voltage_source'].index(V)-1][V.node1-1] +=1
                b_matrix[len(nodes)+components['voltage_source'].index(V)-1] += V.val
            
        
        if(len(components['current_source'])!=0):   
            for I in components['current_source']:
                if (I.node1!=0):
                    b_matrix[I.node1-1] +=I.val
                if(I.node2!= 0):
                    b_matrix[I.node2-1] -=I.val
                    
             
                    
        
        #This command solves for x in the Matrix equation Mx=b
        x_matrix=np.linalg.solve(M_matrix , b_matrix)
        
        #print all the values of the voltages and currents
    
        print('The node voltages are as follows')
        print('The node voltage at node GND is taken to be 0V')
        for node in nodes:
            if(node!=0):
                    print('The node voltage at node {} is {} V'.format(node,x_matrix[node-1]))
                
        if(len(components['voltage_source'])!=0):
                    for V in components['voltage_source']:
                        print('The current in {} is {} A'.format(V.name,x_matrix[len(nodes)+components['voltage_source'].index(V)-1]))
                        
                        
                        
    except IOError:                                    #for exception case ,exit 
         print('Invalid file')
         exit()

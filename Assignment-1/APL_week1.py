# -*- coding: utf-8
"""
APL WEEK 1 Project by Abhishek Sekar
This code essentially takes a .netlist file and gives an output which 
has the characteristics of each line reversed besides reversing the order
in which these lines occur.
"""

from sys import argv,exit
file_start='.circuit'             #since a valid netlist file starts with .circuit. can be altered as per desire.
file_end='.end'
if(len(argv)==2):                 #We'd like to deal with only two arguments, one being the code's location and the other being the file location
    
    try:
        with open(argv[1]) as f:  
            lines=f.readlines()  
        start=-1
        end=-2
        sentence=list(reversed(lines))   #reversing the lines here, so we just have to reverse the characteristics in a line

        for line in sentence:
            if(file_end==line[:len(file_end)]):
                start=sentence.index(line)            #search for start and end expressions
            elif(file_start==line[:len(file_start)]): #since we're dealing with the reversed versions here start becomes end!
                end=sentence.index(line)
            
        if (start>=end):
            print("not a valid netlist")
            exit(0)   
        for line in sentence[start+1:end]:
            print(' '.join(list(reversed(line.split("#")[0].split())))) #seperates comments from line,reverses characteristics and prints
    except IOError: 
        print('Invalid file')
        exit()        
              
              
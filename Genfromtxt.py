'''
Created on May 18, 2014

@author: Keisuke
'''
from pylab import * 

def Genfromtxt(name,folder):
    First_part = "C:\\Users\\Keiin\\OneDrive\\Desktop\\figures\\"
    Second_part= folder + "\\" + name
    return genfromtxt(First_part + Second_part, delimiter=",")
'''
Created on May 18, 2014

@author: Keisuke
'''
from pylab import * 

def Savetxt2(parameter,savefolder,name):
    First_part="C:\\Users\\Keiin\\OneDrive\\Desktop\\figures\\"
    Second_part=savefolder+"\\"+name
    Last_part=".csv"
    savetxt(First_part+Second_part+Last_part, parameter, delimiter=",",fmt='%s')
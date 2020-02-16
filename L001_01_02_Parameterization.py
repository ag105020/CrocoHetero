'''
Created on Apr 26, 2019
Light limited growth rate of Croco
@author: keiin
'''

from pylab import *
from FigSetting2 import *
from Savefig3 import *
from Savetxt2 import *
from Genfromtxt import *

rcParams.update({'ytick.right': 'True'})
rcParams.update({'ytick.direction': 'in'})
rcParams.update({'xtick.top': 'True'})
rcParams.update({'xtick.direction': 'in'})

#==================================
# Main Funcation
#==================================
def MuZ(a0,b0,c0):
    z = arange(0,150+1,1) #(m) depth
    z0 = 30         #(m) Depth of 1/e light of surface
    I0 = 1000       #(umol m-2 s-1) Surface light intensity
    I = I0*exp(-z/z0) #(umol m-2 s-1)
    Mu = zeros(size(z)) 
    Qc = 2.4*10**(-13)*10**15 #(fmol C cell-1) from Masuda-san's info + Azoto stoichiometry
    
    a = a0 * 24/Qc #a*Mu + b (Kei 215-15~16) 24 is to convert from h-1 to d-1
    b = b0 * 24/Qc
    c = c0 * 24/Qc #c*MU (Kei 215-15~17)
    
    Pmax = 7. #(fmol C cell-1 h-1)
    #I0 = 100 #(umol m-2 s-1)
    I0P = 100
    
    P = Pmax*(1-exp(-I/I0P))
    
    Mu1 = (P-b)/(a+c)
    Mu2 = P/(2*c)
    
    Mu[Mu1<Mu2] = Mu1[Mu1<Mu2]
    Mu[Mu1>=Mu2] = Mu2[Mu1>=Mu2]
    
    Mu[Mu<0] = 0

    if a0 == 35.373:
        figure(1)
        plot(Mu,z,color='b')

        figure(2)
        plot(I,z,color='orange')
    
    elif a0 == 17.708:
        figure(1)
        plot(Mu,z,color='r',dashes=(5,2))
    
    zNan = copy(z).astype(float)
    zNan[Mu==0]=nan
    zMax = max(zNan)
    
    return zMax
        

a00 = 35.373  #From C:\Users\keiin\OneDrive\Desktop\figures\01\HeteroC00\HighRes100\04 BoxSum0.xlsx run:617 04 23 00
b00 = 11.184
c00 = 78.847

a01 = 31.623  #From C:\Users\keiin\OneDrive\Desktop\figures\01\HeteroC00\HighRes101\BoxSum.xlsx run:617 04 23 01
b01 = 11.184
c01 = 62.192

a02 = 17.708  #From C:\Users\keiin\OneDrive\Desktop\figures\01\HeteroC00\HighRes103\BoxSum.xlsx run:617 04 23 03
b02 = 22.369
c02 = 78.847

zMax0 = MuZ(a00,b00,c00)
#MuZ(a01,b01,c01)
zMax2 = MuZ(a02,b02,c02)
#===================================
# Supporting funcation
#===================================
def ud(FigNumber):
    figure(FigNumber)
    gca().invert_yaxis()

#========================
# Plot preparation
#========================

Ylabel = 'Depth (m)'

Savefolder = '02\\05 LightLimitationCroco'
Ylim = (0,150)
Yticks = arange(0,150+25,25)

DPI = 600
#=======================
# Plotting
#=======================

figure(1)
#---For legend---
rcParams.update({'legend.fancybox': False})
plot([],[],color='r',dashes=(5,2),label='Homo.')
plot([],[],color='b',label='Hetero.')
legend(loc=4)
#----------------
ylabel(Ylabel)
xlabel('$\mathit{\mu}$ (d$^{-1}$)')
ylim(Ylim)
yticks(Yticks)
xlim(left=0,right=0.5)
ud(1)
Savefig3(Savefolder,'Growth rate',DPI)

figure(2)
ylabel(Ylabel)
xlabel('Light ($\mu$mol m$^{-2}$ s$^{-1}$)')
ylim(Ylim)
yticks(Yticks)
xlim(left=0)
ud(2)
Savefig3(Savefolder,'Light',DPI)

pHOT = genfromtxt('..\\Data\\HOT_PO4.csv',delimiter=',')[:17].T #[:17] is for up to 150 (m)
figure(3)
errorbar(pHOT[1],pHOT[0],xerr=(pHOT[2],pHOT[2]),fmt='o-',color='r',markeredgecolor='k',elinewidth=1,capthick=2,capsize=7)

ylabel(Ylabel)
ylim(Ylim)
xlabel('PO$_4^{3-}$ ($\mu$mol kg$^{-1}$)')


PO4range = arange(0,0.2+0.01,0.01)
y1 = zMax0*ones(size(PO4range))
y2 = zMax2*ones(size(PO4range))
fill_between(PO4range, y1, y2, where=y1 >= y2, facecolor='#FDD2D0',edgecolor = "none")

xlim(0.03,0.19)
ud(3)
Savefig3(Savefolder,'PO4hot',DPI)


show()


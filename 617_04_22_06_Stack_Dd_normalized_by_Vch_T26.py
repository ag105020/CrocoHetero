'''
Created on 2018/05/16

#----------------------------------------------------------------------------------------------
about 617 04 01 ~
Here I develop a model for Masuda-san
loop can be learned from 317 06 23
stack plot can be learned from 605 06 04 09
Combining ammonium uptake cells can be learned from 317 03 38 01
stackplot from 605_06_04_09

#----------------------------------------------------------------------------------------------
#note
Start with various vitality Rv (heterogenaity)
and its stackplot

Here we got same Vch as in 617 04 00 at 100% nitrogen fixer

#----------------------------------------------------------------------------------------------
@author: ag105020
'''
from pylab import *
from Kuhla_88_lab_data_01 import *
from Savefig2 import *
from Savetxt import *
rc('text', usetex=True)
rc('font', family='serif')
rc('font', serif='Times New Roman')

##########################################################################################################################################
def Croco():
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #1.Parameters----------------------------
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #defining basic parameters
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #============================================================================
    DdA=arange(0.0006,0.6+0.0006,0.0006)
    Ef=0.2
    Rv=0.5
    
    CreatingNormal = 0
    if CreatingNormal == 1:
        Rv=1
    Urv=arange(size(DdA))
    o=zeros(size(DdA))
    
    #oooooooooooooooooooooooooooooooooo
    # Output parameters
    #oooooooooooooooooooooooooooooooooo
    Vch = copy(o)
    
    CCHbiomass = copy(o)
    CCHrespbiomass = copy(o)
    CCHN2fix = copy(o)
    CCHrespN2fix = copy(o)
    CCHextraresp = copy(o)
    
    #============================================================================
    for i in Urv:
        Dd=DdA[i]
        mu25=8.903*10**(-4)           #(kg/m/s) Dynamic viscosity (dounenseikeisuu) of water at T=25C degrees Kestin 1978 (15 Computation of diffusion coefficient for different temperature.xlsx)
        mu30=7.975*10**(-4)          #(kg/m/s) Dynamic viscosity (dounenseikeisuu) of water at T=30C degrees Kestin 1978 (15 Computation of diffusion coefficient for different temperature.xlsx)
        mu26 = mu25 + (mu30-mu25)/(30-25)*(26-25)
        Kelvin=273.15
        TeffectonD25to26=(26+Kelvin)/(25+Kelvin)/mu26*mu25  #(dimensionless) converting diffusion coefficient from 25C to 30C
        O2=0.208   #(mol/m-3) O2 saturation at Salinity = 35 and T=26C (Benson 1984)
        D=Dd/86400          #(s-1) growth rate  
        RLg=2.5e-6    #(m) Radius of the cell
        h=0.0000151
        ep3=h      #diffusivity of glycolipid layer compared to alginate layer (132-23)    
        Rg=1/29            #ratio of glycolipid layer to the cell radius
        R=RLg/(1+Rg)        #Radius of the cell without glycolipid layer(m)    
        Lg=R*Rg             #(m) thickness of glycolipid layer
        Ra=0
        La=R*Ra         #(m) thickness of alginate layer (132-23)
        x0=1/h*(1/R-1/(R+Lg))+ep3/h*1/(R+Lg)+(1/(R+Lg+La))*(1-ep3/h)  #effect of glycolipid layer (Kei 133-11, 134-17)
        r5=1/(x0*R)         #diffusivity efficiency of cell membrane
        V=4/3*pi*RLg**3   #Volume of the cell (m3/cell)
        C=6             #number of C in one carbohydrate (ex. C(glucose)=6)
        Kch=1             #half saturation of protein CH for biomass production
        Do2_25=2.12*10**(-9)  #Diffusion coefficient of O2 in the water at 25C (m2/s)
        Do2_30=Do2_25*TeffectonD25to26  #(m2/2) diffusion coefficient of O2 in the water at 30C
        Do2=Do2_30*r5         #Diffusion coefficient of O2 in the water (m2/s)
        Dch_25=6.728*10**(-10)  #Diffusion coefficient of glucose in the water (m2/s)
        Dch_30=Dch_25*TeffectonD25to26  #Diffusion coefficient of glucoe in water at 30C (m2/s)
        Dch=Dch_30*r5          #Diffusion coefficient of glucose in the water (m2/s)
        Dnh4_25=1.98*10**(-9)  #Diffusion coefficient of ammonium in the water at 25C (m2/s)
        Dnh4_30=Dnh4_25*TeffectonD25to26     #Diffusion coefficient of ammonium in water at 30C (m2/s)
        Dnh4=Dnh4_30*r5        #Diffusion coefficient of ammonium in the water (m2/s)
        O2cri=0.00000000000001
        lmax00=1000         #(molCs-1cell-1m-3) maximum biomass production rate per volume
        lmax=lmax00*V      #(molCs-1cell-1) maximum biomass production rate
        rho=0.22*10**6/12                #(molCbiomass/m3/cell) biomass density in the cell (105-10)

        CHin=100    #(mol/m3): CH concentration in the poring water
        n=6             #number of C in one bacterial biomass (BB)
        a=10.8             #number of H in one BB
        b=2.9             #number of O in one BB
        c=1.5             #number of N in one BB
        d=4*n+a-2*b-3*c #inverse number in the coefficient of BB in the half reaction for one e-
        BB=12*n+1*a+16*b+14*c #(g/mol): mass of BB (bacterial biomass)
        y=c/d           #coefficient of NH4+ in the half reaction of BB production (Kei 97) (Kei 210-87)
        z=1/4           #coefficient of NH4+ in the half reaction of nitrogen-fixation (Kei 97) (Kei 210-87)
        pr=1/1.32      #(dimensionless): protein ration in biomass
        Q1=rho*V             #(molC/cell) biomass per cell
        Q = 60e-15/c*n

        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        #Computation of fe0, fpr and fn considering material, redox and energy balance
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        #=====================================
        # Computing M Kei 210-90,91
        #=====================================
        M=1+(1-Rv)/Rv/Ef
        
        #=====================================
        fs00=(1/y)/(1/y+M/z)        #The ratio of electron used for protein synthesis  (Kei 97) (Kei 210-87) (Kei 210-90)
        fn00=(M/z)/(1/y+M/z)        #The ration of electron used for nitrogen fixation (Kei 97) (Kei 210-87) (Kei 210-90)
        #=====================================
        ep=0.22                   #energy efficiency for the production of energy and the consumption of energy
        dgc0=41.35      #The free energy necessary for the half reaction of glucose production (kJ/e-mol)
        dgATP=50      #(kJ/ATP): energy produced by the reaction of ATP -> ADP (147-19)
        dgn=2*dgATP-dgc0*ep   #The free energy necessary (dg) for the half reaction of nitrogen fixation (kJ/e-mol)
        dgp=35.09-dgc0  #dg for production of pyruvate from glucose (kJ/e-mol))
        dgpc=3.33*1/d*(12*n+1*a+16*b+14*c)  #dg for the production of BB (bacterial biomass) from pyruvate) (147-17)
        dgr=-120.07     #-dg for the energy production pathway (kJ/e-mol)   
    
        if dgp<0:       
            ep1=1/ep    #change ep1 depending on the sign of dgp    
        else:
            ep1=ep
        A=(fn00*dgn+fs00*(dgp/ep1+dgpc/ep))/(-ep*dgr) #A is related to fs0 and fe0
        fe0=A/(1+A)     #the ratio of electron used for energy production
        fs0=1-fe0       #the ratio of electron used for biomass synthesis+nitrogen fixation
        fpr=fs0*fs00    #the ratio of electron used for biomass synthesis
        fn=fs0*fn00     #the ratio of electron used for nitrogen fixation        
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #2.--Stoichiometry (to get E1(E for the case O2cri>O2in))------------------
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        S1=array([["","CH","H2O","CO2","O2","HCO3-","NH4+","N2","H2","BB","H+","e-"],
                  ["'-Rd",1/24,0.25,-0.25,0.,0.,0.,0.,0.,0.,-1.,-1.],
                  ["Ra",0.,-0.5,0.,0.25,0.,0.,0.,0.,0.,1.,1.],
                  ["Rpr",0.,-(2*n-b+c)/d,(n-c)/d,0.,c/d,c/d,0.,0.,-1/d,1.,1.],
                  ["Rn",0.,0.,0.,0.,0.,-0.25,0.125,-0.125,0.,1.25,1.]])
        
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        #for creating S2 (f*R for electron acceptance)
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        Mu=array([[1],[fe0],[fpr],[fn]])  #column of f
        S2=copy(S1[1:5,1:12])           #use copy so S2 does not respond to the change in S1
        S2=S2.astype(float64)
        S21=arange(1,12,1)
        S22=arange(0,5,1)
        S22=S22.reshape(5,1)
        S2=vstack((S21,S2))
        S2=hstack((S22,S2))
        S2[1:5,1:12]=Mu*S2[1:5,1:12]
        
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        #for creating S3 (f*R for electron donation)
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        S3=vstack((S2[1,1:12],S2[1,1:12],S2[1,1:12]))
        S3=Mu[1:]*S3
        
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        #for creating S4 (f*R for "electron donation + electron acceptance")
        # and RR, which is the entire reaction
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        S4=S2[2:,1:]+S3     #S4 is f*R for "electron donation + electron acceptance"
        RR=S4[0]+S4[1]+S4[2]
        RR1=copy(RR)
        RR1=vstack((S1[0,1:],RR1))
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        #output of each array into CSV files
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        savetxt("S1.csv", S1, delimiter=",",fmt='%s')
        savetxt("S2.csv", S2, delimiter=",",fmt='%2.3f')
        savetxt("S3.csv", S3, delimiter=",",fmt='%2.8f')
        savetxt("S4.csv", S4, delimiter=",",fmt='%2.8f')
        savetxt("RR.csv", RR, delimiter=",",fmt='%2.8f')
        savetxt("RR1.csv", RR1, delimiter=",",fmt='%s')
     
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        #Getting yield (Y) and the ratio of CO2 production rate to CH consumption (E)
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        Y=-(RR[8]*n)/(C*RR[0])              #Yield
        E1=1/Y-1                            #the ratio of CO2 production rate to CH consumption
        
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #3.Calculation to obtain E2 (E when [O2]in>[O2]cri
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        r1=C*(S4[1,0]+S4[2,0])/(-S4[1,8]*n) #(CH consumption except for energy production)/(biomass synthesis)
        r2=A
        r3=S4[0,3]/(S4[0,0]*C)              #(O2 consumption)/(CH consumption for energy production))
        p=-S4[1,8]*n                        #(molC/e-mol) CH consumption for protein production
        h1=(S4[1,0]+S4[2,0])*C               #(molC/e-mol) CH consumption for other than energy production
        alp=p/h1                             #(CO2 from (Ra-Rd))/(CH for protein production) (see 72-5)
        beta=-(S4[1,2]+S4[1,4]+S4[2,2])/h1   #(CO2 from (Rpr-Rd) + CO2 from (RN-Rd))/(CH for protein production) (see 72-5)
      
        #AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        #=====================================================================================
        #NH4+ absorption case
        #=====================================================================================
        epNH4=0.54
        epNH4=ep
        
        if dgp<0:       
            epNH4_1=1/epNH4    #change ep1 depending on the sign of dgp    
        else:
            epNH4_1=epNH4
        
        A1=(dgp/epNH4_1+dgpc/epNH4)/(-epNH4*dgr) #A is related to fs0 and fe0
        fe01=A1/(1+A1)     #the ratio of electron used for energy production
        fs01=1-fe01       #the ratio of electron used for biomass synthesis+nitrogen fixation
        fpr1=fs01         #the ratio of electron used for biomass synthesis
        
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        #for creating S2 (f*R for electron acceptance)
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
                
        Mua=array([[1],[fe01],[fpr1],[0.0]])  #column of f
        S2a=copy(S1[1:5,1:12])               #use copy so S2 does not respond to the change in S1
        S2a=S2a.astype(float64) 
        S21a=arange(1,12,1)
        S22a=arange(0,5,1)
        S22a=S22a.reshape(5,1)
        S2a=vstack((S21a,S2a))
        S2a=hstack((S22a,S2a))
        S2a[1:5,1:12]=Mua*S2a[1:5,1:12]
        
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        #for creating S3 (f*R for electron donation)
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        S3a=vstack((S2a[1,1:12],S2a[1,1:12],S2a[1,1:12]))
        S3a=Mua[1:]*S3a
        
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        #for creating S4 (f*R for "electron donation + electron acceptance")
        # and RR, which is the entire reaction
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        S4a=S2a[2:,1:]+S3a     #S4 is f*R for "electron donation + electron acceptance"
        RRa=S4a[0]+S4a[1]+S4a[2]
        RR1a=copy(RRa)
        RR1a=vstack((S1[0,1:],RR1a))    
        pa=-S4a[1,8]*n                        #(molC/e-mol) CH consumption for protein production
        ha=(S4a[1,0]+S4a[2,0])*C               #(molC/e-mol) CH consumption for other than energy production
        alpa=pa/ha                             #(CO2 from (Ra-Rd))/(CH for protein production) (see 72-5)
        betaa=-(S4a[1,2]+S4a[1,4]+S4a[2,2])/ha   #(CO2 from (Rpr-Rd) + CO2 from (RN-Rd))/(CH for protein production) (see 72-5)
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        #Getting yield (Y) and the ratio of CO2 production rate to CH consumption (E)
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    
        Y1=-(RRa[8]*n)/(C*RRa[0])              #Yield
        E3=1/Y1-1 
        #AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAaa
      
        CHc=(D*Kch*Q)/(lmax-D*Q)            #(mol/m3) CH concentration inside the cell
        ls=lmax*CHc/(CHc+Kch)                   #(mol/s/cell): biomass synthesis rate   
        lsgrow=ls#/Rv                #(mol/s/cell): biomass synthesis rate for the growing cells
        lmaxgrow=lmax#/Rv          #(mol s-1 cell-1):maximum biomass synthesis rate for the growing cells (181-11)
        dd=(4*pi*R*Do2*(O2-O2cri))/(alp*lmaxgrow*r1*r3)   
        E2=dd*((CHc+Kch)/(CHc))+beta/alp    #E when [O2]in>[O2]cri
        O2in=-(lsgrow*r1*r2*r3)/(4*pi*Do2*R)+O2     #(mol/m3): O2 concentration in the cell given there is not Ocri limit (72-6)
        
        #important part#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        if O2in<O2cri:
            E=E1
        else:
            E=E2
    
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        VchF=ls*(1+E)                        #(molC cell-1 s-1)
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #Nitrogen fixing cells (F: Nitrogen Fixing) :For stack plot of CH consumption (178-1~2)(181-13) Based on 601_00_03
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        CCHbiomassF=ls/(-S4[1,8]*n)*S4[1,0]*C                #(molC cell-1 s-1) consumption of CH for biomsas production
        CCHN2fixF=ls/(-S4[1,8]*n)*S4[2,0]*C                  #(molC cell-1 s-1) consumption of CH for N2 fixation
        CCHbalancedrespF=ls/(-S4[1,8]*n)*S4[0,0]*C           #(molC cell-1 s-1) consumption of CH for balanced respiration
    
        CCHrespbiomassF=CCHbalancedrespF/(fn00*dgn+fs00*(dgp/ep1+dgpc/ep))*(fs00*(dgp/ep1+dgpc/ep))     #(molC cell-1 s-1)
                                                        #Consumption of CH for respiration for biomass production
        CCHrespN2fixF=CCHbalancedrespF/(fn00*dgn+fs00*(dgp/ep1+dgpc/ep))*(fn00*dgn)   #(molC cell-1 s-1)
                                                        #Consumption of CH for respiration for N2 fixation
        CCHextrarespF=VchF-CCHbiomassF-CCHN2fixF-CCHbalancedrespF    #(molC cell-1 s-1) Carbohydrate consumption for extra respiration
        
        N2fixFperSec = ls/(-S4[1,8]*n) * S4[2,6]*2    #(molN cell-1 s-1) per second
        N2fixF = N2fixFperSec*3600   #(molN cell-1 s-1) per hour
        if N2fixF > 6.1e-15:
            DdA[i] = nan

        #AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        # NH4 limited case  (NF: Non-N-Fixing)
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        #AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        VchNF=ls*(1+E3)                               #(molC cell-1 s-1)
        
        CCHbiomassNF=ls/(-S4a[1,8])*S4a[1,0]                #(molC cell-1 s-1) consumption of CH for biomsas production
        CCHrespbiomassNF=ls/(-S4a[1,8])*S4a[0,0]                 #(molC cell-1 s-1) consumption of CH for balanced respiration
        
        #CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        # Combining
        #CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        Vch[i]=VchF*Rv + VchNF*(1-Rv)     #(molC cell-1 s-1)
        CCHbiomass[i] = CCHbiomassF*Rv + CCHbiomassNF*(1-Rv)
        CCHrespbiomass[i] = CCHrespbiomassF*Rv + CCHrespbiomassNF*(1-Rv)
        CCHN2fix[i] = CCHN2fixF*Rv
        CCHrespN2fix[i] = CCHrespN2fixF*Rv
        CCHextraresp[i] = CCHextrarespF*Rv
    
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #5.Plot
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    rcParams.update({'font.size': 30,
                     'lines.markersize':10,
                     'lines.markeredgewidth':0.5})
    rcParams.update({'xtick.major.pad': 15})
    rcParams.update({'xtick.major.pad': 15})
    rcParams.update({'font.serif': 'Times New Roman'})
    rcParams.update({'figure.autolayout': True})
    rcParams['figure.figsize']=8,6.5
    rcParams.update({'figure.facecolor':'W'})
    rcParams.update({'text.usetex': True})   #to call real latex
    rcParams.update({'text.latex.preamble': ['\\usepackage[greek,english]{babel}']})
    rcParams.update({'lines.linewidth':2.5})
    rcParams.update({'patch.edgecolor':'none'})
    
    def gft(FileName):
        return genfromtxt('C:\\Users\\Keiin\\OneDrive\\Desktop\\figures\\01\\HeteroC00\\HighRes13\\'+FileName+'.csv',delimiter=',')
    
    if CreatingNormal == 1:
        Savetxt(Vch,'HeteroC00\\HighRes13','Vch0-06')
    else:
        N0=gft('Vch0-06')
        
        Xticks=arange(0,0.5,0.1)
        Ym=1.20000000001
        Xlab='$\mu$ (d$^{-1}$)'
        
        figure(1)
        plot(arange(0,0.45+0.001,0.001),ones(size(arange(0,0.45+0.001,0.001))),dashes=(2.5,2.5),color='k',linewidth=2.5)
        plot(DdA,Vch/N0)
        xlabel(Xlab)
        xticks(Xticks)
        ylabel('$C_S$/$C_S^0$')
        ylim(0,Ym)
        Savefig2('HeteroC00\\HighRes14\\Revision','D0',600)
    
        StackPolotColors=('#008000','#64FF32','#FF9900','#FFFF00','cyan')
#         figure(2)
#         stackplot(DdA,CCHbiomass/N0,CCHN2fix/N0,CCHrespbiomass/N0,CCHrespN2fix/N0,CCHextraresp/N0,colors=StackPolotColors)
#         plot(DdA,Vch/N0,'k')
#         xticks(Xticks)
#       #  ylim(0,Ym)
#         xlabel(Xlab)
#         ylabel('$C_S$/$C_S^0$')
#         Savefig2('HeteroC00\\HighRes14\\Revision','D1',300)
#         
        figure(3)  #Plot in fmol C cell-1 d-1
        Conversion = 1e15*3600 #from molC cell-1 s-1 to fmol cell-1 day-1
        stackplot(DdA,CCHbiomass*Conversion,CCHrespbiomass*Conversion,CCHN2fix*Conversion,CCHrespN2fix*Conversion,CCHextraresp*Conversion,colors=StackPolotColors)
        plot(DdA,Vch*Conversion,'k')
        #ylim(0,1.20000001)
        xlabel(Xlab)
        ylabel('$C_S$ (fmol C cell$^{-1}$ h$^{-1}$)')
        xticks(Xticks)
        #xlim(0,1)
        Savefig2('HeteroC00\\HighRes14\\Revision','D2',600)

        
    return 

Croco()

show()

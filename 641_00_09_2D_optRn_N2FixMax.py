'''
Created on 2018/05/16

#----------------------------------------------------------------------------------------------
about 640 00 00 ~
This is for plotting

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

#GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
# Reading file and prepare
#GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG

Step1=0.001
Min1=0.001
Max1=1
Array1=arange(Min1,Max1+Step1,Step1)

Step2=0.001
Min2=0.001
Max2=1
Array2=arange(Min2,Max2+Step2,Step2)

def gft(FileName):
    return genfromtxt('C:\\Users\\Keisuke\\Desktop\\figures\\HeteroC00\\HighRes13\\'+FileName+'.csv',delimiter=',')

Box = gft('01_Rn-Ef')/gft('01_Rn-Ef_Normal')

Conversion = 1e15*3600  #conversion from molC cell-1 s-1 to fmolC cell-1 h-1)
Box2 = gft('01_Rn-Ef') * Conversion

B1A = gft('01_B1A')
Box_N2fixF = gft('01_N2FixF')
B1A = zeros(size(B1A))
B1A[B1A==0] = nan

U1=arange(size(Array1))
U2=arange(size(Array2))

for i in U2:            #i is y axis
    for j in U1:        #j is x axis
        if j>0:
            if Box[i,j]<=1.00000001 and Box[i,j-1]>1:
                B1A[i]=Array2[j]
    if i>0:
        if isnan(B1A[i])==True and (B1A[i-1]!=Min1 and isnan(B1A[i-1])==False):# and (B1A[i-1]!=Min1):
            B1A[i]=Min1

B1B=Max1-argmin(Box[:,::-1],axis=1)/(Max1/Step1)     #here Argmin gives first min value but at Ef=1, there are multiple and I want it to point to the last min, so I change left and right of Box_Vch
                                                          #that is the "1-" and ",[::-1,:]" part for 

Box_mask = ma.masked_where(Box_N2fixF > 6.1e-15,Box)

Box2_mask = ma.masked_where(Box_N2fixF > 6.1e-15,Box2)
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

rcParams['image.cmap'] = 'hot'

BackColor = '#7B7B7B'

#Color inversion
cmap = cm.RdYlBu
cmaplist = [cmap(i) for i in range(cmap.N)]
cmaplist = flipud(cmaplist)
cmap = cmap.from_list('Accent', cmaplist, cmap.N)

cmap = cm.gnuplot2

###############

figure(1)
gca().set_axis_bgcolor(BackColor)
pcolormesh(Array1,Array2,Box_mask,cmap=cmap)
xlabel('$f_N$')
ylabel('$E_{N}$')
title('$C_S$/$C_S^0$')
xlim((0.075,Max1))
ylim((0,Max2))
xticks(arange(0.0,Max1+0.0001,0.2))
yticks(arange(0.0,Max2+0.0001,0.2))
#colorbar()

cbar = colorbar(ticks=[0.5,0.6,0.7,0.8,0.9,1.0])
clim(0.5,1)
ax = gca()
cbar.ax.set_yticklabels(['$<$0.5','0.6','0.7','0.8','0.9','$>$1.0'])
#plot(B1A,Array2,color='white')
plot(B1B,Array2,'--',dashes=(7, 3),linewidth=2.5,color='#00FFFF')
Savefig2('HeteroC00\\HighRes15\\gnuplot2',4,300)
#Savetxt(Box,'HeteroC00\\HighRes10','01_Box')
#Savetxt(B1A,'HeteroC00\\HighRes09','01_B1A')

figure(2)
plot(B1A,Array2)
plot(B1B,Array2)

figure(3)
gca().set_axis_bgcolor(BackColor)
pcolormesh(Array1,Array2,Box2_mask,cmap=cmap)
xlabel('$f_N$')
ylabel('$E_{N}$')
title('$C_S$ (fmol C cell$^{-1}$ h$^{-1}$)')
xlim((0.075,Max1))
ylim((0,Max2))
xticks(arange(0.0,Max1+0.0001,0.2))
yticks(arange(0.0,Max2+0.0001,0.2))
#colorbar()
cbar = colorbar(ticks=[10,15,20,25,30])
clim(10,30)
ax = gca()
cbar.ax.set_yticklabels(['10','15','20','25','$>$30'])
Savefig2('HeteroC00\\HighRes15\\gnuplot2',41,300)

print(Box_mask.max())
print(Box_mask.min())

show()

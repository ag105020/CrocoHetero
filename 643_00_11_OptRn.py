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

Step1=0.0006
Min1=0.0006
Max1=0.6
Array1=arange(Min1,Max1+Step1,Step1)

Step2=0.001
Min2=0.001
Max2=1
Array2=arange(Min2,Max2+Step2,Step2)

def gft(FileName):
    return genfromtxt('C:\\Users\\Keisuke\\Desktop\\figures\\HeteroC00\\HighRes13\\'+FileName+'.csv',delimiter=',')

Box=gft('03_Dd-Rn')/gft('03_NormalD')
B1A=gft('03_B1A')
B1B=argmin(Box,axis=0)[0:Max1/Step1]/(Max2/Step2)+Step2
Box_N2fixF = gft('03_N2FixF')

Conversion = 1e15*3600 #from molC cell-1 s-1 to fmol cell-1 day-1
Box2=gft('03_Dd-Rn')*Conversion

Box_mask = ma.masked_where(Box_N2fixF > 6.1e-15,Box)
#B1B=argmin(Box_mask,axis=0)[0:Max1/Step1]/(Max2/Step2)+Step2
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
#BackColor = '#C9C9C9'
#BackColor = 'k'

#Color inversion
cmap = cm.RdYlBu
cmap = cm.gnuplot2
# cmaplist = [cmap(i) for i in range(cmap.N)]
# cmaplist = flipud(cmaplist)
# cmap = cmap.from_list('Accent', cmaplist, cmap.N)

###############

figure(1)
gca().set_axis_bgcolor(BackColor)
pcolormesh(Array1,Array2,Box_mask,cmap=cmap)
xlabel('$\mu$ (d$^{-1}$)')
ylabel('$f_N$')
title('$C_S$/$C_S^0$')
xlim((Min1,Max1))
ylim((Min2,Max2))
xticks(arange(0.,Max1+0.0001,0.2))
yticks(arange(0.,Max2+0.0001,0.2))
#=====================
cbar = colorbar(ticks=arange(0,1.40001,0.2))
clim(0.,1.4)
ax = gca()
cbar.ax.set_yticklabels(arange(0,1.40001,0.2))
#=====================
#plot(B1A,Array2,color='w')
plot(Array1,B1B,'--',dashes=(7, 3),color='#00FFFF')
Savefig2('HeteroC00\\HighRes15\\gnuplot2',6,300)
Savetxt(Box,'HeteroC00\\HighRes15\\gnuplot2','03_Box')

figure(3)

gca().set_axis_bgcolor(BackColor)
pcolormesh(Array1,Array2,Box2_mask,cmap=cmap)
xlabel('$\mu$ (d$^{-1}$)')
ylabel('$f_N$')
title('$C_S$ (fmol C cell$^{-1}$ h$^{-1}$)')
xlim((Min1,Max1))
ylim((Min2,Max2))
xticks(arange(0.,Max1+0.0001,0.2))
yticks(arange(0.,Max2+0.0001,0.2))
#=====================
cbar = colorbar(ticks=arange(0,50.00001,10))
clim(0.,50.00001)
ax = gca()
cbar.ax.set_yticklabels(arange(0,51,10))
#=====================
Savefig2('HeteroC00\\HighRes15\\gnuplot2',61,300)


figure(2)
plot(B1A,Array2)
plot(Array1,B1B)

show()

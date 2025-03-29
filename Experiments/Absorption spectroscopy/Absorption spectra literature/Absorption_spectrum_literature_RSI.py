# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#%% Imports
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter,MultipleLocator

#%% Catalog of defect peaks
# Nameofdefect, peak number, OA peak eV , FWHM eV, sources (list)

#GR91 = Griscom 1991
#GI19 = Girard 2019
#KA22 = Kashaykin 2022
#MO20 = Morana, 2020
#DI89 = Dianov Chernov 1989  

#{"peak":,"FWHM":,"sources":["GI19"]}
#FWHM=1e-6 when unknown 

DICOPEAKS = {
    "ODC-II":[
        {"peak":5.05,"FWHM":0.32,"sources":["GI19"]},
        {"peak":3.15,"FWHM":0.30,"sources":["GI19"]},
        {"peak":6.9,"FWHM":0.4,"sources":["GI19"]}
        ],
    "ODC-I":[
        {"peak":7.6,"FWHM":(0.5+0.6)/2,"sources":["GI19"]}
        ],
    "NBOHC":[
        {"peak":1.97,"FWHM":0.17,"sources":["GI19"]},
        {"peak":4.8,"FWHM":1.0,"sources":["GI19"]},
        {"peak":6.4,"FWHM":1.7,"sources":["GI19"]}
        ],
    "E-prime":[
        {"peak":5.8,"FWHM":0.7,"sources":["GI19"]}
        ],
    "STH1":[
        {"peak":2.61,"FWHM":1.2,"sources":["GI19"]},
        {"peak":1.88,"FWHM":(0.2+0.5)/2,"sources":["GI19"]}
        ],
    "STH2":[
        {"peak":2.16,"FWHM":(0.3+0.6)/2,"sources":["GI19"]},
        {"peak":1.63,"FWHM":(0.3+0.7)/2,"sources":["GI19"]}
        ],
    "POL":[
        {"peak":3.8,"FWHM":0.2,"sources":["GI19"]},
        {"peak":4.2,"FWHM":0.6,"sources":["GI19"]},
        {"peak":7.3,"FWHM":0.2,"sources":["GI19"]},
        {"peak":7.5,"FWHM":0.1,"sources":["GI19"]}
        ],
    "POR":[
        {"peak":2.02,"FWHM":1e-6,"sources":["GI19"]},
        {"peak":4.08,"FWHM":1e-6,"sources":["GI19"]},
        {"peak":5.02,"FWHM":1e-6,"sources":["GI19"]},
        {"peak":2.0,"FWHM":1e-6,"sources":["GI19"]},
        {"peak":4.8,"FWHM":1e-6,"sources":["GI19"]}
        ],
    "STE":[
        {"peak":3.7,"FWHM":1e-6,"sources":["GI19"]},
        {"peak":4.6,"FWHM":1e-6,"sources":["GI19"]},
        {"peak":6.4,"FWHM":1e-6,"sources":["GI19"]}
        ],
    "STEX":[
        {"peak":4.2,"FWHM":1.16,"sources":["GI19"]},
        {"peak":5.3,"FWHM":0.78,"sources":["GI19"]}
        ],
    "Composite 1-eV band":[
        {"peak":1.2,"FWHM":0.56,"sources":["KA22","MO20"]},
        {"peak":0.93,"FWHM":0.42,"sources":["KA22","MO20"]}
        ]
    }

DICOIMPUR = {
    "O3":[
        {"peak":4.8,"FWHM":(0.8+0.86)/2,"sources":["GI19"]}
        ],
    "O2":[
        {"peak":0.97,"FWHM":0.013,"sources":["GI19"]},
        {"peak":1.62,"FWHM":0.013,"sources":["GI19"]}
        ],
    "Cl0":[
        {"peak":3.26,"FWHM":1e-6,"sources":["GI19"]},
        {"peak":3.65,"FWHM":1e-6,"sources":["GI19"]}
        ],
    "Cl2":[
        {"peak":3.78,"FWHM":0.6,"sources":["GI19"]},
        {"peak":2.3,"FWHM":1e-6,"sources":["GI19"]}
        ],
    "H(I)":[], #not observed
    "LTIRA":[] #non gaussian, see Girard 2019, Dianov Chernov 1989  
    }

DICOLTIRA = {"D":0.39, #eV (even if they say eV^-1?) #From Ka22 and Dianov-Chernov 89
             "S":0.11, #eV (even if they say eV^-1?)
             "R":0.43, #eV (even if they say eV^-1?)
             "Ep":0.65 #eV
             }

#%% Plotting functions
COLORS_v1=["tab:blue","tab:orange","tab:green","tab:red","tab:purple","tab:brown","tab:pink","tab:gray","tab:olive","tab:cyan","b","g","r","c","m","y","k"]
COLORS_v2 = ["b", "g", "r", "c", "m", "y", "k"]

my_cmap=plt.get_cmap("tab20")  # https://matplotlib.org/3.1.1/gallery/color/colormap_reference.html
COLORS_v3 = [my_cmap(i) for i in range(my_cmap.N)]
# COLORS = [COLORS[2*i] for i in range(10)]+[COLORS[2*i+1] for i in range(10)]

#Let's try to be more colorblind-friendly
COLORS = ["#D81B60",
          "#1E88E5",
          "#FFC107",
          "#004D40",
          "#79A6CD",
          "#96CAC4",
          "#AA2968",
          "#5D8759",
          "#6E378B",
          "#C6742C",
          "#36BB0D",
          "#FA93A4"]
MARKERS = ["o",
           "v",
           "^",
           "s",
           "p",
           "P",
           "*",
           "X",
           "D",
           "h",
           "<",
           ">"]

def Sigma(FWHM):
    return FWHM/(2*np.sqrt(2*np.log(2)))

def plotlists_LTIRA(D,S,R,Ep,magnitude=1e3, depth=1e-3,model="Kubo-Greenwood",toleV=1e-5):
    assert depth<1
    def LTIRAfunc(E):
        if model=="Kubo-Greenwood":
            Q=D+S/(D-S)
            B0=magnitude
            return B0*E*np.exp(-E/S)*(np.exp(E/Q)-np.exp(Ep/Q))
        elif model=="exp":
            C=magnitude
            return C*np.exp(-E/R)
        else:
            raise NameError("You fucked up")
    if model=="Kube-Greenwood":
        x=np.linspace(0.5,2)
        indmax=np.argmax(LTIRAfunc(x))
        maxLTIRA,EmaxLTIRA=LTIRAfunc(x[indmax]),x[indmax]
        #find left side
        left=[EmaxLTIRA/2,EmaxLTIRA]
        while left[1]-left[0]>toleV:
            mid=(left[1]+left[0])/2
            if LTIRAfunc(mid)>maxLTIRA*depth:
                left[1]=mid
            else:
                left[0]=mid
        left=(left[0]+left[1])/2
        #find right side
        right=[EmaxLTIRA,EmaxLTIRA+3.0]
        while right[1]-right[0]>toleV:
            mid=(right[1]+right[0])/2
            if LTIRAfunc(mid)>maxLTIRA*depth:
                right[0]=mid
            else:
                right[1]=mid
        right=(right[0]+right[1])/2
        E_ARRAY=np.linspace(left, right)
    else:
        E_ARRAY=np.linspace(0,-R*np.log(depth))
    return E_ARRAY,LTIRAfunc(E_ARRAY)
    

def plotlists_peak(E_eV,FWHM_eV,magnitude=1,depth=1e-3,normalizemax=False):
    assert depth<1
    if FWHM_eV==1e-6:
        return [E_eV,E_eV],[depth,10]
    sigma = Sigma(FWHM_eV)
    Delta_E = np.sqrt(-2*np.log(depth))*sigma
    E_ARRAY = np.linspace(E_eV-Delta_E, E_eV+Delta_E)
    coef = 1
    if not(normalizemax):
        coef = magnitude/(sigma*np.sqrt(2*np.pi))
    PEAK = coef*np.exp(-(E_ARRAY-E_eV)**2/(2*sigma**2))
    return E_ARRAY,PEAK

eV = 1.60218e-19 #J
h_Planck = 6.62607015e-34 #m^2.kg/s
c_light = 299792458 #m/s

def eV2WLnm(E):
    return h_Planck*c_light/(E*eV*1e-9) #in nm
def WLnm2eV(lmbda):
    return h_Planck*c_light/(lmbda*1e-9*eV)

print(eV2WLnm(6.2),WLnm2eV(200))


# fig, ax = plt.subplots(constrained_layout=True)
# ax.semilogy(*plotlists_peak(2.61,1.2),color="red")
# ax.set_xlabel("E (eV)")
# ax.set_ylabel("Absorption (arbitrary unit)")
# ax.set_title('Absorption peaks from literature')
# secax = ax.secondary_xaxis('top', functions=(eV2WLnm, WLnm2eV))
# secax.set_xlabel('$\lambda$ (nm)')
# secax.set_xscale('log')
# secax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
# plt.legend()
# plt.show()


    
#%% PLOTS

#%% Plot whole catalog

for DICO in DICOPEAKS,DICOIMPUR:
    fig, ax = plt.subplots(constrained_layout=True,figsize=(7.5,3))
    fig.set_dpi(1000)
    #ax.set_aspect(0.6)
    axisfontsize=14
    mksize=4
    
    LEGENDCOLORS,LEGENDLABELS=[plt.Line2D([0], [0], color='black', linestyle='dashdot', linewidth=1),plt.Line2D([0], [0], linestyle="--",color="k")],["1550 nm","Undefined FWHM"]
    for i,defect in enumerate(DICO):
        for j,dicpeak in enumerate(DICO[defect]):
            lstyle="-"
            if dicpeak["FWHM"]==1e-6:
                lstyle="--"
            if j==0:
                LEGENDCOLORS += [plt.Line2D([0], [0], color=COLORS[i],marker=MARKERS[i],markersize=mksize)]
                sources=""
                for s in dicpeak["sources"]:
                    if len(sources)>0:
                        sources+=","
                    sources+=s
                if defect in ["Composite 1-eV band"]:
                    LEGENDLABELS += [defect+f"\n[{sources}]"]
                else:
                    LEGENDLABELS += [defect+f" [{sources}]"]
            xdata,ydata=plotlists_peak(dicpeak["peak"],dicpeak["FWHM"],normalizemax=False)
            ax.semilogy(xdata,ydata,color=COLORS[i],linestyle=lstyle)
            ax.plot(xdata[np.argmax(ydata)],ydata[np.argmax(ydata)],marker=MARKERS[i],color=COLORS[i],markersize=mksize)
    
    if DICO==DICOPEAKS:
        ax.semilogy(*plotlists_LTIRA(DICOLTIRA["D"],DICOLTIRA["S"],DICOLTIRA["R"],DICOLTIRA["Ep"],
                                     magnitude=6e3, depth=1e-2,model="Kubo-Greenwood",toleV=1e-5),
                    color=COLORS[i+1],marker="o",markersize=2)
        LEGENDCOLORS += [plt.Line2D([0], [0], color=COLORS[i+1],marker="o",markersize=2)]
        LEGENDLABELS += [f"LTIRA (K-G model)\n[DI89,GI19]"]
    
    line_position_nm = 1550
    ax.axvline(WLnm2eV(line_position_nm), color='black', linestyle='dashdot', linewidth=1, label='1550 nm')
    # ax.set_xlim(0.3,11)
    ax.set_xlim(0.3,9.4)
    ax.set_xlabel("E (eV)",fontsize=axisfontsize)
    ax.set_ylabel("Absorption (arbitrary unit)",fontsize=axisfontsize)
    #ax.set_title('Absorption peaks from literature')
    ax.set_xticks(range(1, 10))
    
    secax = ax.secondary_xaxis('top', functions=(eV2WLnm, WLnm2eV))
    secax.set_xlabel('$\lambda$ (nm)',fontsize=axisfontsize)
    # Set all the ticks on the secondary x-axis
    secax.set_xticks([150, 200, 300, 400, 500, 600, 700, 800, 900, 1000,2000])
    
    # Define a custom function to format tick labels
    def format_func(value, tick_number):
        if value in [150, 200, 300, 400, 500, 600, 1000, 2000]:
            return int(value)
        else:
            return ""
    
    # Apply the custom tick formatter
    secax.get_xaxis().set_major_formatter(plt.FuncFormatter(format_func))
    ax.get_yaxis().set_major_formatter(plt.FuncFormatter((lambda v,ticknb : "")))
    #plt.legend(LEGENDCOLORS,LEGENDLABELS,loc="best",fontsize=7)
    plt.legend(LEGENDCOLORS, LEGENDLABELS, loc='center left', bbox_to_anchor=(0.79, 0.5), fontsize=6, borderaxespad=1.5)
    
    save = False
    if save:
        base=r"C:/Users/aurel/Downloads/"
        if DICO==DICOPEAKS:
            plt.savefig(base+"Absorption_peaks_defects.pdf")
        elif DICO==DICOIMPUR:
            plt.savefig(base+"Absorption_peaks_impurities.pdf")
    else:
        plt.show()



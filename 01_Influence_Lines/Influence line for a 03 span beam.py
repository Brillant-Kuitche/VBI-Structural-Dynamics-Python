# -*- coding: utf-8 -*-
"""
Created on Sun May  3 15:45:03 2026

@author: Brillant
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

input("\n Ce Programme permet de tracer la ligne d'influence pour un point quelconque sur une poutre à 3 travées. \n Les Reaction d'appuis sont R0; R1; R2; R3, Respectivement aux noeuds 0, 1, 2, et 3. ")
L = float(input(" Entrer la longueur de la travée mediane L = "))
psi = float(input(" Entrer la valeur de Psi: Psi = "))
#s = float(input(" Entrer la valeur s_tild ( compris entre [0; 1] et indiquand la position relative de la force uniitaire sur l'une des travées). s_tild = "))
nbr_point = int(input(" Vous souhaitez que la force passe par combien de point sur chaque travée: "))
x = float(input(" Vous souhaitez la ligne d'influence de quelle point sur la poutre ? Entrer sa posion en (m) par rapport au Noeud 0 (X compris entre 0 et "+ str(f"{(L+2*psi*L):.3f}")+ "m) : X = "))
n = int(input(" Vous souhaitez: \n      1. La ligne d'influence du moment fléchissant en x = " +str(f"{x:.3f}")+ "m \n      2. La ligne d'influence de l'effort tranchant en x = " +str(f"{x:.3f}")+ "m \n      3. La ligne d'influence des reactions R0, R1, R2, et R3 \n Choix N°: "))

s = np.linspace(0, 1, nbr_point)
l = psi*L
reaction_de_contrainte = []
reaction_de_contrainte_du_TS = []
moment_en_X = []
moment_TS = []
eff_tranchant_en_X = []
eff_tranchant_TS = []
position_X = []

#  la force est sur la travée 1
for i in range(len(s)-1):
    a = float(s[i])*l
    a_transl = (a+1.2)*(a+1.2 <= l) + (a+1.2-l)*(a+1.2 >= l)
    position_X.append(a)
    #position_X_TS.append(a)
    b = l-a
    b_transl = (l - a_transl)*(a+1.2 <= l) + (L - a_transl)*((a+1.2 >= l))
    X2 = -(a*b*L*(l+a))/(l*(L**2-4*(l+L)**2))
    X2_transl = (-(a_transl*b_transl*L*(l+a_transl))/(l*(L**2-4*(l+L)**2)))*(a+1.2 <= l) + (((a_transl*b_transl*(L+b_transl))/L-(2*a_transl*b_transl*(L+a_transl)*(L+l))/(L**2))*L/(-(L**2)+4*(l+L)**2))*((a+1.2 >= l))
    X1 = -2*(L+l)*X2/L
    X1_transl = (-2*(L+l)*X2_transl/L)*(a+1.2 <= l) + (-a_transl*b_transl*(L+a_transl)/(L**2)-2*(l+L)*X2_transl/L)*((a+1.2 >= l))
    R0 = (X1+b)/l
    R0_transl = ((X1_transl+b_transl)/l)*(a+1.2 <= l) + (X1_transl/l)*(a+1.2 >= l)
    R3 = X2/l
    R3_transl = (X2_transl/l)*(a+1.2 <= l) + (X2_transl/l)*(a+1.2 >= l)
    R2 = (X1-R3*(l+L))/L
    R2_transl = ((X1_transl-R3_transl*(l+L))/L)*(a+1.2 <= l) + ((X1_transl-R3_transl*(l+L)+a_transl)/L)*(a+1.2 >= l)
    R1 = 1-(R0+R2+R3)
    R1_transl = (1-(R0_transl+R2_transl+R3_transl))
    R = np.array([R0, R1, R2, R3])
    R_TS = np.array([R0+R0_transl, R1+R1_transl, R2+R2_transl, R3+R3_transl])
    reaction_de_contrainte.append(R)
    reaction_de_contrainte_du_TS.append(R_TS)
    
    if 0<=x<=a:
        Mp = R0*x
        Tp = R0
    elif a<=x<=l:
        Mp = R0*x-(x-a)
        Tp = R0-1
    elif l<=x<=l+L:
        Mp = R0*x+R1*(x-l)-(b+(x-l))
        Tp = R0+R1-1
    elif l+L<=x<=l+L+l:
        Mp = R0*x+R1*(x-l)+R2*(x-l-L)-(b+L+(x-l-L))
        Tp = -R3
    else:
        input(" Le point que vous avez specifiez n'appartient pas à la poutre")
     
    if a+1.2 <= l:    
        Mp_transl = (R0_transl*x)*(0<=x<=a_transl) + (R0_transl*x-(x-a_transl))*(a_transl<=x<=l) + (R0_transl*x+R1_transl*(x-l)-(b_transl+(x-l)))*(l<=x<=l+L) + (R0_transl*x+R1_transl*(x-l)+R2_transl*(x-l-L)-(b_transl+L+(x-l-L)))*(l+L<=x<=l+L+l)   
        Tp_transl = (R0_transl)*(0<=x<=a_transl) + (R0_transl-1)*(a_transl<=x<=l) + (R0_transl+R1_transl-1)*(l<=x<=l+L) + (-R3_transl)*(l+L<=x<=l+L+l)
    else:
        Mp_transl = (R0_transl*x)*(0<=x<=l) + (R0_transl*x+R1_transl*(x-l))*(l<=x<=l+a_transl) + (R0_transl*x+R1_transl*(x-l)-(x-l-a_transl))*(l+a_transl<=x<=l+L) + (R0_transl*x+R1_transl*(x-l)+R2_transl*(x-l-L)-(b_transl+x-l-L))*(l+L<=x<=l+L+l)
        Tp_transl = (R0_transl)*(0<=x<=l) + (R0_transl+R1_transl)*(l<=x<=l+a_transl) + (R0_transl+R1_transl-1)*(l+a_transl<=x<=l+L) + (-R3_transl)*((l+L<=x<=l+L+l))
        
    moment_en_X.append(Mp)
    moment_TS.append(Mp + Mp_transl)
    eff_tranchant_en_X.append(Tp)
    eff_tranchant_TS.append(Tp + Tp_transl)
    
#  la force est sur la travée 2
for i in range(np.shape(s)[0]-1):
    a = float(s[i])*L
    a_transl = (a+1.2)*(a+1.2 <= L) + (a+1.2-L)*(a+1.2 >= L)
    position_X.append(l+a)
    b = L-a
    b_transl = (L - a_transl)*(a+1.2 <= L) + (l - a_transl)*((a+1.2 >= L))
    X2= ((a*b*(L+b))/L-(2*a*b*(L+a)*(L+l))/(L**2))*L/(-(L**2)+4*(l+L)**2)
    X2_transl= (((a_transl*b_transl*(L+b_transl))/L-(2*a_transl*b_transl*(L+a_transl)*(L+l))/(L**2))*L/(-(L**2)+4*(l+L)**2))*(a+1.2 <= L) + ((a_transl*b_transl*(l+b_transl)*2*(l+L))/(l*(L**2-4*(l+L)**2)))*(a+1.2 >= L)
    X1= -a*b*(L+a)/(L**2)-2*(l+L)*X2/L
    X1_transl= (-a_transl*b_transl*(L+a_transl)/(L**2)-2*(l+L)*X2_transl/L)*(a+1.2 <= L) + (-L*X2_transl/(2*(l+L)))*(a+1.2 >= L)
    R0 = X1/l
    R0_transl = X1_transl/l
    R1 = (X2+b-R0*(l+L))/L
    R1_transl = ((X2_transl+b_transl-R0_transl*(l+L))/L)*(a+1.2 <= L) + ((X2_transl-R0_transl*(l+L))/L)*(a+1.2 >= L)
    R3 = X2/l
    R3_transl = (X2_transl/l)*(a+1.2 <= L) + ((X2_transl+a_transl)/l)*(a+1.2 >= L)
    R2 = 1-(R0+R1+R3)
    R2_transl = 1-(R0_transl+R1_transl+R3_transl)
    R = np.array([R0, R1, R2, R3])
    R_TS = np.array([R0+R0_transl, R1+R1_transl, R2+R2_transl, R3+R3_transl])
    reaction_de_contrainte.append(R)
    reaction_de_contrainte_du_TS.append(R_TS)
    
    if 0<=x<=l:
        Mp = R0*x
        Tp = R0
    elif l<=x<=l+a:
        Mp = R0*x+R1*(x-l)
        Tp = R0+R1
    elif l+a<=x<=l+L:
        Mp = R0*x+R1*(x-l)-(x-l-a)
        Tp = R0+R1-1
    elif l+L<=x<=l+L+l:
        Mp = R0*x+R1*(x-l)+R2*(x-l-L)-(b+x-l-L)
        Tp = -R3
    else:
        input(" Le point que vous avez specifiez n'appartient pas à la poutre")
    
    if a+1.2 <= L:
        Mp_transl = (R0_transl*x)*(0<=x<=l) + (R0_transl*x+R1_transl*(x-l))*(l<=x<=l+a_transl) + (R0_transl*x+R1_transl*(x-l)-(x-l-a_transl))*(l+a_transl<=x<=l+L) + (R0_transl*x+R1_transl*(x-l)+R2_transl*(x-l-L)-(b_transl+x-l-L))*(l+L<=x<=l+L+l)
        Tp_transl = (R0_transl)*(0<=x<=l) + (R0_transl+R1_transl)*(l<=x<=l+a_transl) + (R0_transl+R1_transl-1)*(l+a_transl<=x<=l+L) + (-R3_transl)*((l+L<=x<=l+L+l))
    else:
        Mp_transl = (R0_transl*x)*(0<=x<=l) + (R0_transl*x+R1_transl*(x-l))*(l<=x<=L+l) + ( R3_transl*(l-(x-l-L))-(a_transl-(x-l-L)))*(L+l<=x<=l+L+a_transl) + (R3_transl*(b_transl-(x-L-l-a_transl)))*(l+L+a_transl<=x<=l+L+l)
        Tp_transl = (R0_transl)*(0<=x<=l) + (R0_transl+R1_transl)*(l<=x<=L+l) + (-R3_transl+1)*((L+l<=x<=l+L+a_transl)) + (-R3_transl)*(l+L+a_transl<=x<=l+L+l)    

    moment_en_X.append(Mp)
    moment_TS.append(Mp + Mp_transl)
    eff_tranchant_en_X.append(Tp)
    eff_tranchant_TS.append(Tp + Tp_transl)

# la force est sur la travée 3
for i in range(np.shape(s)[0]):
    a = float(s[i])*l
    a_transl = (a+1.2)*(a+1.2 <= l)
    position_X.append(l+L+a)
    b = l-a
    b_transl = (l - a_transl)*(a+1.2 <= l)
    X1 = (-a*b*L*(l+b))/(l*(L**2-4*(l+L)**2))
    X1_transl = (-a_transl*b_transl*L*(l+b_transl))/(l*(L**2-4*(l+L)**2))
    X2 = -2*(l+L)*X1/L
    X2_transl = -2*(l+L)*X1_transl/L
    R0 = X1/l
    R0_transl = X1_transl/l
    R1 = (X2-R0*(l+L))/L
    R1_transl = (X2_transl-R0_transl*(l+L))/L
    R3 = (X2+a)/l
    R3_transl = (X2_transl+a_transl)/l
    R2 = (1-(R0+R1+R3))
    R2_transl = (1-(R0_transl+R1_transl+R3_transl))
    R = np.array([R0, R1, R2, R3])
    R_TS = np.array([R0+R0_transl, R1+R1_transl, R2+R2_transl, R3+R3_transl])
    reaction_de_contrainte.append(R)
    reaction_de_contrainte_du_TS.append(R_TS)

    if 0<=x<=l:
        Mp = R0*x
        Tp = R0
    elif l<=x<=L+l:
        Mp = R0*x+R1*(x-l)
        Tp = R0+R1
    elif L+l<=x<=l+L+a:
        Mp = R3*(l-(x-l-L))-(a-(x-l-L))
        Tp = -R3+1
    elif l+L+a<=x<=l+L+l:
        Mp = R3*(b-(x-L-l-a))
        Tp = -R3
    else:
        input(" Le point que vous avez specifiez n'appartient pas à la poutre")
        
    Mp_transl = (R0_transl*x)*(0<=x<=l) + (R0_transl*x+R1_transl*(x-l))*(l<=x<=L+l) + ( R3_transl*(l-(x-l-L))-(a_transl-(x-l-L)))*(L+l<=x<=l+L+a_transl) + (R3_transl*(b_transl-(x-L-l-a_transl)))*(l+L+a_transl<=x<=l+L+l)
    Tp_transl = (R0_transl)*(0<=x<=l) + (R0_transl+R1_transl)*(l<=x<=L+l) + (-R3_transl+1)*((L+l<=x<=l+L+a_transl)) + (-R3_transl)*(l+L+a_transl<=x<=l+L+l)    
    
    moment_en_X.append(Mp)
    moment_TS.append(Mp + Mp_transl)
    eff_tranchant_en_X.append(Tp)
    eff_tranchant_TS.append(Tp + Tp_transl)

if n==1:
    plt.figure(figsize=(10,6))
    plt.plot(position_X, moment_en_X, color='blue', linewidth=2, linestyle='-', label='M')
    plt.title("Ligne d'influence du moment en X = "+str(x)+" m)", fontsize=14)
    plt.xlabel("x(m)", fontsize=12)
    plt.ylabel("M(kNm)", fontsize=12) 
    plt.legend()
 #  plt.gca().invert_xaxis()
    # affiche les valeur de l'axe y avec 3 decimales et impose un pas 1 sur l'axe x
    ax = plt.gca() # Recupere les axes
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(2))
    
    plt.fill_between(position_X, moment_en_X, color='blue',alpha=0.2)
    plt.margins(0)
    plt.grid()
    plt.tight_layout()
    plt.savefig("Graph des LI/LI_Moment.png", dpi=800)
    plt.show()
    
    #LI du TS
    plt.figure(figsize=(10,6))
    plt.plot(position_X, moment_TS, color='blue', linewidth=2, linestyle='-', label='M')
    plt.title("Ligne d'influence du moment de TS en X = "+str(x)+" m)", fontsize=14)
    plt.xlabel("x(m)", fontsize=12)
    plt.ylabel("M(kNm)", fontsize=12) 
    plt.legend()
    ax = plt.gca() # Recupere les axes
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(2)) 
    plt.fill_between(position_X, moment_TS, color='blue',alpha=0.2)
    plt.margins(0)
    plt.grid()
    plt.tight_layout()
    plt.savefig("Graph des LI/LI_Moment_TS.png", dpi=800)
    plt.show()
    
    
elif n==2:
    plt.figure(figsize=(10,6)) 
    plt.plot(position_X, eff_tranchant_en_X, color='blue', linewidth=2, linestyle='-', label='T')
    plt.title("Ligne d'influence de l' effort tranchant en X = "+str(x)+" m)", fontsize=14)
    plt.xlabel("x(m)", fontsize=12)
    plt.ylabel("T(kN)", fontsize=12) 
    plt.legend()
    # affiche les valeur de l'axe y avec 3 decimales
    ax = plt.gca()
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(2))
    plt.fill_between(position_X, eff_tranchant_en_X, color='blue',alpha=0.2)
    plt.margins(0)
    plt.grid()
    plt.tight_layout()
    plt.savefig("Graph des LI/LI_effort tranchant.png", dpi=800) 
    plt.show()
    
    #LI du TS
    plt.figure(figsize=(10,6)) 
    plt.plot(position_X, eff_tranchant_TS, color='blue', linewidth=2, linestyle='-', label='T')
    plt.title("Ligne d'influence de l' effort tranchant de TS en X = "+str(x)+" m)", fontsize=14)
    plt.xlabel("x(m)", fontsize=12)
    plt.ylabel("T(kN)", fontsize=12) 
    plt.legend()
    # affiche les valeur de l'axe y avec 3 decimales
    ax = plt.gca()
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(2))
    plt.fill_between(position_X, eff_tranchant_TS, color='blue',alpha=0.2)
    plt.margins(0)
    plt.grid()
    plt.tight_layout()
    plt.savefig("Graph des LI/LI_effort tranchant_TS.png", dpi=800) 
    plt.show()
    
elif n==3:
    for i in range(4):
        Reaction_i = []
        for j in range(len(reaction_de_contrainte)):
            Reaction_i.append(float(reaction_de_contrainte[j][i]))
            
        plt.figure(figsize=(10,6)) 
        plt.plot(position_X, Reaction_i, color='blue', linewidth=2, linestyle='-', label='R')
        plt.title("Ligne d'influence de la Reaction R"+str(i), fontsize=14)
        plt.xlabel("x(m)", fontsize=12)
        plt.ylabel("R"+str(i)+"(kN)", fontsize=12) 
        plt.legend()
        ax = plt.gca()
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
        ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
        plt.fill_between(position_X, Reaction_i, color='blue',alpha=0.2)
        plt.margins(0)
        plt.grid()
        plt.tight_layout()
        plt.savefig("Reactions/LI_Reaction_"+str(i)+".png", dpi=800) 
        plt.show()
else:
    input("Veuillez choisir 1, 2 ou 3 ") 


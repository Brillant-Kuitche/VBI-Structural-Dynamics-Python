# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 09:33:58 2026

@author: Brillant

Ce programme calcul les lignes d'influences longitudinales du projet PSI DP
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

input("\n\n            Ligne d'influence d'une poutre à 2 travée. \n \n Les Reaction d'appuis sont R0; R1; R2, Respectivement aux noeuds 0, 1 et 2.  ")
l1 = float(input(" Entrer la portée de la 1ère travée (en m): L1 = "))
l2 = float(input(" Entrer la portée de la 2ème travée (en m): L2 = "))
nbr_point = int(input(" Vous souhaitez que la force passe par combien de point sur chaque travée: "))
s = np.zeros([1,nbr_point]) # position relative de la charge roulante
s_barre = np.linspace(0, 1, nbr_point)
s[0, :] = s_barre
x = np.zeros([1,2*nbr_point-1]) # position relative de la section
x[0, 0:nbr_point] = s*l1 # Travée1
x[0, nbr_point:] = s[0, 1:]*l2 + l1 # Travée2
#x_trans = (x+1.2)*(x+1.2<=l1+l2)

# Intialisation des matrice
M = np.zeros([2*nbr_point-1,2*nbr_point-1]) #initialisation
M_trans = np.zeros([2*nbr_point-1,2*nbr_point-1]) #initialisation de la 2eme matrice de l'effort translater de 1.2m
M_TS = np.zeros([2*nbr_point-1,2*nbr_point-1])

T = np.zeros([2*nbr_point-1,2*nbr_point-1]) #initialisation
T_trans = np.zeros([2*nbr_point-1,2*nbr_point-1]) #initialisation de la 2eme matrice de l'effort translater de 1.2m 
T_TS = np.zeros([2*nbr_point-1,2*nbr_point-1]) #initialisation de la matrice des effort de TS

Reactions = np.zeros([3,2*nbr_point-1])  #initialisation
Reactions_trans = np.zeros([3,2*nbr_point-1])  #initialisation de la 2eme matrice de l'effort translater de 1.2m
Reactions_TS = np.zeros([3,2*nbr_point-1])  #initialisation

# Charge roulante dans la travée 1
# inconnues hyperstatiques
X1 = np.zeros([1,nbr_point])  # Initialisation des inconnues hyperstatiques
X1_trans = np.zeros([1,nbr_point])

X1 = (-s*l1*(l1-s*l1)*(l1+s*l1))/(2*l1*(l1+l2))

X1_trans = ((-(s*l1+1.2)*(l1-(s*l1+1.2))*(l1+s*l1+1.2))/(2*l1*(l1+l2)))*(s*l1+1.2<=l1) + ((-(s*l1+1.2-l1)*(l2-(s*l1+1.2-l1))*(2*l2-(s*l1+1.2-l1)))/(2*l2*(l1+l2)))*(s*l1+1.2>l1)

#Réactions d'appui
Reactions[0,0:nbr_point] = (X1+l1-s*l1)/l1
Reactions_trans[0,0:nbr_point] = ((X1_trans+l1-(s*l1+1.2))/l1)*(s*l1+1.2<=l1) + (X1_trans/l1)*((s*l1+1.2>l1))

Reactions[2,0:nbr_point] = X1/l2
Reactions_trans[2,0:nbr_point] = (X1_trans/l2)*(s*l1+1.2<=l1) + ((X1_trans+(s*l1+1.2-l1))/l2)*(s*l1+1.2>l1)

Reactions[1,0:nbr_point] = 1-(Reactions[0,0:nbr_point]+Reactions[2,0:nbr_point])
Reactions_trans[1,0:nbr_point] = 1-(Reactions_trans[0,0:nbr_point]+Reactions_trans[2,0:nbr_point])

# Defintion des foctions de calculs des moments
def M1_T1(y,R0):  # si 0<=y<=s*l1
    return R0*y
def M2_T1(y,R0,x_force_T1):  # si s*l1<=y<=l1
    return R0*y-(y-x_force_T1)
def M3_T1(y, R2):
    return R2*(l1+l2-y)
#Moment fléchissants et efforts tranchants
for x_force in range(s.shape[1]):
    M[x_force, :] = (Reactions[0,x_force]*x)*(0<=x)*(x<=s[0, x_force]*l1)+(Reactions[0,x_force]*x-(x-s[0, x_force]*l1))*(s[0, x_force]*l1<x)*(x<=l1)+(Reactions[2,x_force]*(l1+l2-x))*(l1<x)*(x<=l1+l2)
    T[x_force, :] = Reactions[0,x_force]*(0<=x)*(x<=s[0, x_force]*l1)+(Reactions[0,x_force]-1)*(s[0, x_force]*l1<x)*(x<=l1)+(-Reactions[2,x_force])*(l1<x)*(x<=l1+l2)

    if s[0, x_force]*l1+1.2 <= l1:
        M_trans[x_force, :] = (Reactions_trans[0,x_force]*x)*(0<=x)*(x<s[0, x_force]*l1+1.2)+(Reactions_trans[0,x_force]*x-(x-(s[0, x_force]*l1+1.2)))*(s[0, x_force]*l1+1.2<=x)*(x<l1)+(Reactions_trans[2,x_force]*(l1+l2-x))*(l1<=x)*(x<=l1+l2)
        T_trans[x_force, :] = Reactions_trans[0,x_force]*(0<=x)*(x<=s[0, x_force]*l1+1.2)+(Reactions_trans[0,x_force]-1)*(s[0, x_force]*l1+1.2<x)*(x<=l1)+(-Reactions_trans[2,x_force])*(l1<x)*(x<=l1+l2)
    else:
        M_trans[x_force, :] = (Reactions_trans[0,x_force]*x)*(0<=x)*(x<l1)+(Reactions_trans[0,x_force]*x + Reactions_trans[1,x_force]*(x-l1))*(l1<=x)*(x<l1+s[0, x_force]*l1+1.2-l1)+(Reactions_trans[2,x_force]*(l1+l2-x))*(l1+s[0, x_force]*l1+1.2-l1<=x)*(x<=l1+l2)
        T_trans[x_force, :] = Reactions_trans[0,x_force]*(0<=x)*(x<l1)+(Reactions_trans[0,x_force]+Reactions_trans[1,x_force])*(l1<=x)*(x<l1+s[0, x_force]*l1+1.2-l1)+(-Reactions_trans[2,x_force])*(l1+s[0, x_force]*l1+1.2-l1<=x)*(x<=l1+l2)


# Charge roulante dans la travée 2
# inconnues hyperstatiques
X1 = np.zeros([1,nbr_point])  # Initialisation des inconnues hyperstatiques
X1_trans = np.zeros([1,nbr_point])

X1 = (-s*l2*(l2-s*l2)*(2*l2-s*l2))/(2*l2*(l1+l2))
X1_trans = ((-(s*l2+1.2)*(l2-(s*l2+1.2))*(2*l2-(s*l2+1.2)))/(2*l2*(l1+l2)))*(s*l2+1.2 <= l2)
#Réactions d'appui
Reactions[0,nbr_point-1:] = X1/l1
Reactions_trans[0,nbr_point-1:] = X1_trans/l1

Reactions[2,nbr_point-1:] = (X1+s*l2)/l2
Reactions_trans[2,nbr_point-1:] = (X1_trans+(s*l2+1.2))/l2

Reactions[1,nbr_point-1:] = 1-(Reactions[0,nbr_point-1:]+Reactions[2,nbr_point-1:])
Reactions_trans[1,nbr_point-1:] = (1-(Reactions_trans[0,nbr_point-1:]+Reactions_trans[2,nbr_point-1:]))*(s*l2+1.2 <= l2)

# Defintion des foctions de calculs des moments
def M1_T2(y,R0):  # si 0<=y<=l1
    return R0*y
def M2_T2(y, R0, R1):  # si l1<=y<=l1+s*l2
    return R0*y + R1*(y-l1)
def M3_T2(y, R2):
    return R2*(l1+l2-y)
    
#Moment fléchissants et efforts tranchants
for x_force1 in range(s.shape[1]):
    x_force = x_force1 + nbr_point-1
    M[x_force, :] = (Reactions[0,x_force]*x)*(0<=x)*(x<l1)+(Reactions[0,x_force]*x + Reactions[1,x_force]*(x-l1))*(l1<=x)*(x<l1+s[0, x_force1]*l2)+(Reactions[2,x_force]*(l1+l2-x))*(l1+s[0, x_force1]*l2<=x)*(x<=l1+l2)
    T[x_force, :] = Reactions[0,x_force]*(0<=x)*(x<l1)+(Reactions[0,x_force]+Reactions[1,x_force])*(l1<=x)*(x<l1+s[0, x_force1]*l2)+(-Reactions[2,x_force])*(l1+s[0, x_force1]*l2<=x)*(x<=l1+l2)

    M_trans[x_force, :] = (Reactions_trans[0,x_force]*x)*(0<=x)*(x<l1)+(Reactions_trans[0,x_force]*x + Reactions_trans[1,x_force]*(x-l1))*(l1<=x)*(x<l1+s[0, x_force1]*l2+1.2)+(Reactions_trans[2,x_force]*(l1+l2-x))*(l1+s[0, x_force1]*l2+1.2<=x)*(x<=l1+l2)
    T_trans[x_force, :] = Reactions_trans[0,x_force]*(0<=x)*(x<l1)+(Reactions_trans[0,x_force]+Reactions_trans[1,x_force])*(l1<=x)*(x<l1+s[0, x_force1]*l2+1.2)+(-Reactions_trans[2,x_force])*(l1+s[0, x_force1]*l2+1.2<=x)*(x<=l1+l2)

M_TS = M + M_trans
T_TS = T + T_trans
Reactions_TS = Reactions + Reactions_trans

# Exploitation des lignes d'influence
# Moments [x, Omega+, Omega-, ymax, ymin, ymax1, ymax2, ymin1, ymin2, ]
data_moments = np.zeros([2*nbr_point,7])
#entete = np.array([["x (m)", "Omega +", "Omega -", "ymax (LM2)", "ymin (LM2)", "ymax1+ymax2 (TS)", "ymin1+ymin2 (TS)"]])
x_Matrix = np.repeat(x, 2*nbr_point-1, 0).T
x_T1 = (x_Matrix<=l1)
x_T2 = (x_Matrix>l1)*(x_Matrix<=l1+l2)

data_moments[1:, 0:1] = x.T
data_moments[1:, 1] = np.sum(M*(M>=0)*x_T1,0)*(1/(nbr_point-1))*l1 + np.sum(M*(M>=0)*x_T2,0)*(1/(nbr_point-1))*l2  #omega+
data_moments[1:, 2] = np.sum(M*(M<=0)*x_T1,0)*(1/(nbr_point-1))*l1 + np.sum(M*(M<=0)*x_T2,0)*(1/(nbr_point-1))*l2  #omega-
data_moments[1:, 3] = np.max(M*(M>=0),0) #ymax
data_moments[1:, 4] = np.min(M*(M<=0),0) #ymin
data_moments[1:, 5] = np.max(M_TS*(M>=0),0) #ymax de TS
data_moments[1:, 6] = np.min(M_TS*(M<=0),0) #ymin de TS

# Effort tranchat [x, Omega+, Omega-, ymax, ymin, ymax1, ymax2, ymin1, ymin2, ]
data_shear = np.zeros([2*nbr_point,7])
data_shear[1:, 0:1] = x.T
data_shear[1:, 1] = np.sum(T*(T>=0)*x_T1,0)*(1/(nbr_point-1))*l1 + np.sum(T*(T>=0)*x_T2,0)*(1/(nbr_point-1))*l2  #omega+
data_shear[1:, 2] = np.sum(T*(T<=0)*x_T1,0)*(1/(nbr_point-1))*l1 + np.sum(T*(T<=0)*x_T2,0)*(1/(nbr_point-1))*l2  #omega-
data_shear[1:, 3] = np.max(T*(T>=0),0) #ymax
data_shear[1:, 4] = np.min(T*(T<=0),0) #ymin
data_shear[1:, 5] = np.max(T_TS*(T>=0),0) #ymax de TS
data_shear[1:, 6] = np.min(T_TS*(T<=0),0) #ymin de TS

# Graphiques
#Moment
plt.figure(figsize=(10,6))
plt.plot(x[0, :], T_TS[:, 500], color='blue', linewidth=2, linestyle='-', label='T TS')
plt.title("Ligne d'influence de l'effort tranchant TS en X = "+str(float(x[0, 500]))+" m)", fontsize=14)
plt.xlabel("x(m)", fontsize=12)
plt.ylabel("T TS(kNm)", fontsize=12) 
plt.legend()
# affiche les valeur de l'axe y avec 3 decimales et impose un pas 1 sur l'axe x
ax = plt.gca() # Recupere les axes
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
plt.fill_between(x[0, :], T_TS[:, 500], color='blue',alpha=0.2)
plt.margins(0)
plt.grid()
plt.tight_layout()
plt.savefig("TS effort.png", dpi=800)
plt.show()

#  Mts
plt.figure(figsize=(10,6)) 
plt.plot(x[0, :], T_TS[:, 1000], color='blue', linewidth=2, linestyle='-', label='T TS')
plt.title("Ligne d'influence de l'effort tranchant TS en X = "+str(float(x[0, 1000]))+" m)", fontsize=14)
plt.xlabel("x(m)", fontsize=12)
plt.ylabel("T TS(kNm)", fontsize=12) 
plt.legend()
# affiche les valeur de l'axe y avec 3 decimales
ax = plt.gca()
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
plt.fill_between(x[0, :], T_TS[:, 1000], color='blue',alpha=0.2)
plt.margins(0)
plt.grid()
plt.tight_layout()
plt.savefig("LI_effortTS2001.png", dpi=800) 
plt.show()







































"""
# Reaction 
plt.figure(figsize=(10,6))
plt.plot(x[0, :], Reactions[0, :], color='blue', linewidth=2, linestyle='-', label='R0')
plt.plot(x[0, :], Reactions[1, :], color='red', linewidth=2, linestyle='-', label='R1')
plt.plot(x[0, :], Reactions[2, :], color='green', linewidth=2, linestyle='-', label='R2')

plt.title("Ligne d'influence des reactions", fontsize=14)
plt.xlabel("x(m)", fontsize=12)
plt.ylabel("R(kN)", fontsize=12) 
plt.legend()
# affiche les valeur de l'axe y avec 3 decimales et impose un pas 1 sur l'axe x
ax = plt.gca() # Recupere les axes
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
plt.fill_between(x[0, :], Reactions[0, :], color='blue',alpha=0.1)
plt.fill_between(x[0, :], Reactions[1, :], color='red',alpha=0.1)
plt.fill_between(x[0, :], Reactions[2, :], color='green',alpha=0.1)
plt.margins(0)
plt.grid()
plt.tight_layout()
plt.savefig("LI_Reactions.png", dpi=800)
plt.show()

# Reaction du TS
plt.figure(figsize=(10,6))
plt.plot(x[0, :], Reactions_TS[0, :], color='blue', linewidth=2, linestyle='-', label='R0_TS')
plt.plot(x[0, :], Reactions_TS[1, :], color='red', linewidth=2, linestyle='-', label='R1_TS')
plt.plot(x[0, :], Reactions_TS[2, :], color='green', linewidth=2, linestyle='-', label='R2_TS')

plt.title("Ligne d'influence des reactions du TS", fontsize=14)
plt.xlabel("x(m)", fontsize=12)
plt.ylabel("R_TS(kN)", fontsize=12) 
plt.legend()
# affiche les valeur de l'axe y avec 3 decimales et impose un pas 1 sur l'axe x
ax = plt.gca() # Recupere les axes
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
plt.fill_between(x[0, :], Reactions_TS[0, :], color='blue',alpha=0.1)
plt.fill_between(x[0, :], Reactions_TS[1, :], color='red',alpha=0.1)
plt.fill_between(x[0, :], Reactions_TS[2, :], color='green',alpha=0.1)
plt.margins(0)
plt.grid()
plt.tight_layout()
plt.savefig("LI_Reactions_TS.png", dpi=800)
plt.show()
"""
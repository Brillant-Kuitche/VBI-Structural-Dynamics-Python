# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 23:55:42 2026

@author: BRILLANT
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# définition de la table de connectivité
#elt: [nd_dep,nd_arr,E,L,A_dep,A_arr]
E=36000000 #kN/m²
gamma= 25 #kN/m^3
L=10 #m

n = int(input("\n Entrer le nombre d'élémént de la structure: "))
phi_1 = float(input(" Entrer le diamètre inferieur (en m) de la structure: D_inf = "))
phi_2 = float(input(" Entrer le diamètre superieur (en m) de la structure: D_sup = "))
nb_noeuds = n+1
x = np.linspace(0, L, nb_noeuds)
z = np.zeros(len(x))

# Discretisation de l'element
plt.figure(figsize=(5,7))
plt.plot(z, x, color='blue', marker='o', markersize=5, linewidth=2, linestyle='-', label='Discrétisation')
plt.title(" Discrétisation en "+str(n)+" élèments", fontsize=14)
plt.ylabel("x(m)", fontsize=12)
plt.legend()
plt.margins(0)
plt.grid()
plt.tight_layout()
plt.savefig("Discrétisation.png", dpi=800)
plt.show()

#Code pour ressortir automatiquement la table de connectivité de l'élément poutre 2D
Tab_connect = np.zeros([n,6]) 
r=L/n
L1=L/n
S1=(np.pi)*(phi_1**2)/4
S2=(np.pi)*(phi_2**2)/4
A_arr=S1
for i in range (n):
    L_elt=L1 + i*r
    S_elt=(S2-S1)*L_elt/L + S1
    A_dep=A_arr
    A_arr= S_elt
    nd_dep=i
    nd_arr=i+1
    Caract_elt = np.array([nd_dep,nd_arr,E,L/n,A_dep,A_arr])
    Tab_connect[i] = Tab_connect[i] + Caract_elt
    
#Création de la matrice de rigidité élémentaire d'un element poutre 2D à section linéairement variable
def Matrice_rig_el_poutre_2D(E,L,A_dep,A_arr):
    a=(A_arr - A_dep)/L
    b=A_dep
    E1=E*(((a*L)/2)+b)/L
    E2=E/(4*(np.pi))
    P1=((12*b**2)/L**3)+((12*b*a)/L**2)+((24*a**2)/(5*L))
    P2=((6*b**2)/L**2)+((4*b*a)/L)+(7/5)*a**2
    P3=((6*b**2)/L**2)+((8*b*a)/L)+(17/5)*a**2
    P4=((4*b**2)/L)+((2*b*a))+(8/15)*L*a**2
    P5=((2*b**2)/L)+((2*b*a))+(13/15)*L*a**2
    P6=((4*b**2)/L)+((6*b*a))+(38/15)*L*a**2
    return np.array([[E1,0,0,-E1,0,0],
                     [0,E2*P1,E2*P2,0,-E2*P1,E2*P3],
                     [0,E2*P2,E2*P4,0,-E2*P2,E2*P5],
                     [-E1,0,0,E1,0,0],
                     [0,-E2*P1,-E2*P2,0,E2*P1,-E2*P3],
                     [0,E2*P3,E2*P5,0,-E2*P3,E2*P6]])

#Matrice de rigidité élémentaire des éléments  discrétisés 
Mat_elem=[] # liste devrant stocker les matrices de rig. élementaire
Force_elm=[]
for elt in range(np.shape(Tab_connect)[0]):
    k_elem = Matrice_rig_el_poutre_2D(Tab_connect[elt,2],
                                     Tab_connect[elt,3],
                                     Tab_connect[elt,4],
                                     Tab_connect[elt,5])
    Mat_elem.append(k_elem)
    
    # Vecteur force energetiquement equivalent
    l = Tab_connect[elt,3]
    ai = (Tab_connect[elt,5]-Tab_connect[elt,4])/l
    bi = Tab_connect[elt,4]
    F1 = -gamma*(bi*l/2+ai*(l**2)/6)
    F2 = -20*(elt*l)*l/2-(20*(elt+1)*l-20*elt*l)*3*l/20
    F3 = -20*elt*l*l**2/12-(20*(elt+1)*l-20*elt*l)*l**2/30
    F4 = -gamma*(bi*l/2+ai*(l**2)/3)
    F5 = -20*elt*l*l/2-(20*(elt+1)*l-20*elt*l)*7*l/20
    F6 = 20*elt*l*l**2/12+(20*(elt+1)*l-20*elt*l)*l**2/20
    F = np.array([F1,F2,F3,F4,F5,F6])
    Force_elm.append(F)
Force_elm[n-1][3]= Force_elm[n-1][3]-2000

# definition des matrices de rigidité et des forces 
k_global = np.zeros([nb_noeuds*3,nb_noeuds*3])  
F_global = np.zeros(nb_noeuds*3)
for i in range(np.shape(Tab_connect)[0]):
    # Assemblage des matrices de rigidité
    k_global[3*(int(Tab_connect[i,0])):3*(int(Tab_connect[i,0]))+3,3*(int(Tab_connect[i,0])):3*(int(Tab_connect[i,0]))+3]=k_global[3*(int(Tab_connect[i,0])):3*(int(Tab_connect[i,0]))+3,3*(int(Tab_connect[i,0])):3*(int(Tab_connect[i,0]))+3]+Mat_elem[i][0:3,0:3]
    k_global[3*(int(Tab_connect[i,1])):3*(int(Tab_connect[i,1]))+3,3*(int(Tab_connect[i,1])):3*(int(Tab_connect[i,1]))+3]=k_global[3*(int(Tab_connect[i,1])):3*(int(Tab_connect[i,1]))+3,3*(int(Tab_connect[i,1])):3*(int(Tab_connect[i,1]))+3]+Mat_elem[i][3:6,3:6]
    k_global[3*(int(Tab_connect[i,0])):3*(int(Tab_connect[i,0]))+3,3*(int(Tab_connect[i,1])):3*(int(Tab_connect[i,1]))+3]=k_global[3*(int(Tab_connect[i,0])):3*(int(Tab_connect[i,0]))+3,3*(int(Tab_connect[i,1])):3*(int(Tab_connect[i,1]))+3]+Mat_elem[i][0:3,3:6]
    k_global[3*(int(Tab_connect[i,1])):3*(int(Tab_connect[i,1]))+3,3*(int(Tab_connect[i,0])):3*(int(Tab_connect[i,0]))+3]=k_global[3*(int(Tab_connect[i,1])):3*(int(Tab_connect[i,1]))+3,3*(int(Tab_connect[i,0])):3*(int(Tab_connect[i,0]))+3]+Mat_elem[i][3:6,0:3]
    
    # Assemblage du vecteur force
    F_global[3*(int(Tab_connect[i,0])):3*(int(Tab_connect[i,0]))+3]=F_global[3*(int(Tab_connect[i,0])):3*(int(Tab_connect[i,0]))+3]+Force_elm[i][0:3]
    F_global[3*(int(Tab_connect[i,1])):3*(int(Tab_connect[i,1]))+3]=F_global[3*(int(Tab_connect[i,1])):3*(int(Tab_connect[i,1]))+3]+Force_elm[i][3:6]
    
k_reduit = k_global[3:,3:]                                                                                                                                  
F_reduit = F_global[3:]

sol_deplacement = np.linalg.solve(k_reduit, F_reduit)

# Determination des efforts internes par element
sol_par_noeud=[]
sol_par_noeud.append(np.array([0,0,0]))
moment=[]
eff_tranchant=[]
eff_normal=[]
for i in range(n):
    sol_par_noeud.append(np.array([float(sol_deplacement[3*i]),float(sol_deplacement[3*i+1]),float(sol_deplacement[3*i+2])]))
    x_var =np.linspace(i*l, l*(i+1), 100)
    u1 = float(sol_par_noeud[i][0])
    v1 = float(sol_par_noeud[i][1])
    o1 = float(sol_par_noeud[i][2])
    u2 = float(sol_par_noeud[i+1][0])
    v2 = float(sol_par_noeud[i+1][1])
    o2 = float(sol_par_noeud[i+1][2])
    l = Tab_connect[i,3]
    a = (Tab_connect[i,5]-Tab_connect[i,4])/l
    b = Tab_connect[i,4]
    M = (E/(4*(np.pi)))*((a**2)*(x_var-i*l)**2 + 2*a*b*(x_var-i*l)+b**2)*(v1*(-6/l**2 + 12*(x_var-i*l)/l**3)+ o1*(-4/l+6*(x_var-i*l)/l**2)+v2*(6/l**2 - 12*(x_var-i*l)/l**3)+o2*(-2/l+6*(x_var-i*l)/l**2))
    moment.append(M)
    T= (-E/(4*(np.pi)))*((2*(a**2)*(x_var-i*l)+2*a*b)*(v1*(-6/l**2 + 12*(x_var-i*l)/l**3)+ o1*(-4/l+6*(x_var-i*l)/l**2)+v2*(6/l**2 - 12*(x_var-i*l)/l**3)+o2*(-2/l+6*(x_var-i*l)/l**2)) + ((a**2)*(x_var-i*l)**2 + 2*a*b*(x_var-i*l)+b**2)*(v1*12/l**3+o1*6/l**2-v2*12/l**3+o2*6/l**2))
    eff_tranchant.append(T)
    N = E*(a*(x_var-i*l)+b)*(-u1/l+u2/l)
    eff_normal.append(N)

# Graphique de convergence de U, V, Theta en fonction du nombre d'element (Section variable, D_inf=0.5m et D_sup=0.3m)
u1 = np.array([-0.00419, -0.00443, -0.00449, -0.00451, -0.00452, -0.00453, -0.00453, -0.00453, -0.00454, -0.00454, -0.00454])    
v1 = np.array([-2.34084, -2.35920, -2.35945, -2.35933, -2.35926, -2.35922, -2.35921, -2.35920, -2.35919, -2.35919, -2.35919])
theta = np.array([-0.36859, -0.36090, -0.36074, -0.36076, -0.36078, -0.36078, -0.36079, -0.36079, -0.36079, -0.36079, -0.36079])
x1 = np.array([1,2,3,4,5,6,7,8,9,10,11])
# Graphique de U a section variable
plt.figure(figsize=(8,4))
plt.plot(x1, u1, color='blue', linewidth=2, linestyle='-', label='Deplacement longitudinal')
xmin,xmax = plt.xlim()
ymin,ymax = plt.ylim()
plt.xticks(np.arange(0, xmax, 1))
plt.ylim(bottom=ymin)
plt.title("Courbe des Deplacements: D_inf=0.5m et D_sup=0.3m ", fontsize=14)
plt.xlabel("Nombre d'elements", fontsize=12)
plt.ylabel("U(m)", fontsize=12)
plt.legend()
plt.margins(0)
plt.grid()
plt.tight_layout()
plt.savefig("conv solution U.png", dpi=800)
plt.show()
# Graphique de V a section variable
plt.figure(figsize=(8,4))
plt.plot(x1, v1, color='red', linewidth=2, linestyle='-', label='Deplacement transversal')
xmin,xmax = plt.xlim()
ymin,ymax = plt.ylim()
plt.xticks(np.arange(0, xmax, 1))
plt.ylim(bottom=ymin)
plt.title("Courbe des Deplacements: D_inf=0.5m et D_sup=0.3m ", fontsize=14)
plt.xlabel("Nombre d'elements", fontsize=12)
plt.ylabel("V(m)", fontsize=12)
plt.legend()
plt.margins(0)
plt.grid()
plt.tight_layout()
plt.savefig("conv solution V.png", dpi=800)
plt.show()
# Graphique de theta a section variable
plt.figure(figsize=(8,4))
plt.plot(x1, theta, color='green', linewidth=2, linestyle='-', label='Rotation')
xmin,xmax = plt.xlim()
ymin,ymax = plt.ylim()
plt.xticks(np.arange(0, xmax, 1))
plt.ylim(bottom=ymin)
plt.title("Courbe des Deplacements: D_inf=0.5m et D_sup=0.3m ", fontsize=14)
plt.xlabel("Nombre d'elements", fontsize=12)
plt.ylabel("theta(m/m)", fontsize=12)
plt.legend()
plt.margins(0)
plt.grid()
plt.tight_layout()
plt.savefig("conv solution O.png", dpi=800)
plt.show()

# Graphique de convergence de U, V, Theta en fonction du nombre d'element (Section constante, diametre D=0.4m)
uf = np.array([-0.00446, -0.00446, -0.00446, -0.00446, -0.00446, -0.00446, -0.00446, -0.00446, -0.00446, -0.00446, -0.00446])    
vf = np.array([-4.05256, -4.05256, -4.05256, -4.05256, -4.05256, -4.05256, -4.05256, -4.05256, -4.05256, -4.05256, -4.05256])
theta_f = np.array([-0.55262, -0.55262, -0.55262, -0.55262, -0.55262, -0.55262, -0.55262, -0.55262, -0.55262, -0.55262, -0.55262])
# Graphique de U a section fixe
plt.figure(figsize=(8,4))
plt.plot(x1, uf, color='blue', linewidth=2, linestyle='-', label='Deplacement longitudinal')
xmin,xmax = plt.xlim()
ymin,ymax = plt.ylim()
plt.xticks(np.arange(0, xmax, 1))
plt.ylim(bottom=ymin)
plt.title("Courbe des Deplacements, section fixe: D =0.4m ", fontsize=14)
plt.xlabel("Nombre d'elements", fontsize=12)
plt.ylabel("U(m)", fontsize=12)
plt.legend()
plt.margins(0)
plt.grid()
plt.tight_layout()
plt.savefig("conv solution Uf.png", dpi=800)
plt.show()
# graphique de V a section fixe
plt.figure(figsize=(8,4))
plt.plot(x1, v1, color='red', linewidth=2, linestyle='-', label='Deplacement transversal')
xmin,xmax = plt.xlim()
ymin,ymax = plt.ylim()
plt.xticks(np.arange(0, xmax, 1))
plt.ylim(bottom=ymin)
plt.title("Courbe des Deplacements, section fixe: D =0.4m ", fontsize=14)
plt.xlabel("Nombre d'elements", fontsize=12)
plt.ylabel("V(m)", fontsize=12)
plt.legend()
plt.margins(0)
plt.grid()
plt.tight_layout()
plt.savefig("conv solution Vf.png", dpi=800)
plt.show()
# Graphique de theta a section fixe
plt.figure(figsize=(8,4))
plt.plot(x1, theta, color='green', linewidth=2, linestyle='-', label='Rotation')
xmin,xmax = plt.xlim()
ymin,ymax = plt.ylim()
plt.xticks(np.arange(0, xmax, 1))
plt.ylim(bottom=ymin)
plt.title("Courbe des Deplacements, section fixe: D =0.4m ", fontsize=14)
plt.xlabel("Nombre d'elements", fontsize=12)
plt.ylabel("theta(m/m)", fontsize=12)
plt.legend()
plt.margins(0)
plt.grid()
plt.tight_layout()
plt.savefig("conv solution Of.png", dpi=800)
plt.show()
 
# Solution élement fini     
M = np.concatenate(moment)    
T = np.concatenate(eff_tranchant) 
N = np.concatenate(eff_normal) 
print("\n    **************************************************")
print("    **** Deplacement du noeud superieur du poteau **** ")
print("    **************************************************") 
print("\n    - Deplacement longitudinal: U = "+str(f"{sol_par_noeud[n][0]:.5f}")+" m")
print("    - Deplacement Transversal: V = "+str(f"{sol_par_noeud[n][1]:.5f}")+" m")
print("    - Rotation: Theta = "+str(f"{sol_par_noeud[n][2]:.5f}")+" m/m")
print("\n\n            Sous la supervision de Dr-Ing.SEGNING TCHIOTSOP\n")

# Solution RDM 
x= np.linspace(0, L, n*100)
M_RDM = -6666.6667+1000*x-10*x**3/3
V_RDM = -1000+10*x**2
if a==0:
    N_RDM = -2033.38+3.338*x
else:
    N_RDM = -2033.38+4.90875*x-0.157075*(x**2)

# Tracé des diagrammes de comparaison
# Moment
plt.figure(figsize=(10,6))
plt.plot(M, x, color='blue', linewidth=2, linestyle='-', label='M_EF')
plt.plot(M_RDM, x, color='red', linewidth=2, linestyle='--', label='M_RDM')
plt.title("Courbe des Moments (Discrétisation en "+str(n)+" élèments)", fontsize=14)
plt.xlabel("M(kNm)", fontsize=12)
plt.ylabel("x(m)", fontsize=12)
plt.legend()
plt.gca().invert_xaxis()
plt.fill_between(M, x, color='blue',alpha=0.2)
plt.fill_between(M_RDM, x, color='red',alpha=0.1)
plt.margins(0)
plt.grid()
plt.tight_layout()
plt.savefig("Moment.png", dpi=800)
plt.show()

#Effort tranchant
plt.figure(figsize=(10,6))
plt.plot(T, x, color='blue', linewidth=2, linestyle='-', label='T_EF')
plt.plot(V_RDM, x, color='red', linewidth=2, linestyle='--', label='T_RDM')
plt.title("Courbe des efforts tranchants (Discrétisation en "+str(n)+" élèments)", fontsize=14)
plt.xlabel("T(kN)", fontsize=12)
plt.ylabel("x(m)", fontsize=12)
plt.legend()
plt.gca().invert_xaxis()
plt.margins(0)
plt.fill_between(T, x, color='blue',alpha=0.2)
plt.fill_between(V_RDM, x, color='red',alpha=0.1)
plt.grid()
plt.tight_layout()
plt.savefig("Effort tranchant.png", dpi=800)
plt.show()

#Effort Normal
plt.figure(figsize=(10,6))
plt.plot(N, x, color='blue', linewidth=2, linestyle='-', label='N_EF')
plt.plot(N_RDM, x, color='red', linewidth=2, linestyle='--', label='N_RDM')
plt.title("Courbe des efforts Normales (Discrétisation en "+str(n)+" élèments)", fontsize=14)
plt.xlabel("N(kN)", fontsize=12)
plt.ylabel("x(m)", fontsize=12)
plt.legend()
plt.gca().invert_xaxis()
plt.margins(0)
#plt.fill_between(N, x, color='blue',alpha=0.2)
#plt.fill_between(N_RDM, x, color='red',alpha=0.1)
plt.grid()
plt.tight_layout()
plt.savefig("Effort Normal.png", dpi=800)
plt.show()


# Representation 3D du poteau
r1=phi_1/2
r2=phi_2/2
z=np.linspace(0, L, 50)
#plt.figure(figsize=(15,7))
theta = np.linspace(0, 2*np.pi, 50)
Z, theta = np.meshgrid(z,theta)
a=(S2-S1)/L
b= S1
aire = a*z + b 
R=(4*aire/np.pi)**0.5/2
X=R*np.cos(theta)
Y=R*np.sin(theta)
fig = plt.figure()
ax=fig.add_subplot(111, projection='3d')
ax.plot_surface(X,Y,Z)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.margins(0)
plt.title("Structure en 3D", fontsize=14)
plt.savefig("Structure.png", dpi=800)
plt.show()


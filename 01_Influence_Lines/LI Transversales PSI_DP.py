# -*- coding: utf-8 -*-
"""
Tracé des lignes d'influence de Guyon-Massonet-Barès
On modifie les variable directement à l'interieure du code

@author: Brillant
"""
import numpy as np
import matplotlib.pyplot as plt
from numpy import sinh as sh
from numpy import cosh as ch
import matplotlib.ticker as ticker

b = 4.13898                                                 # mi-largeur du tablier (m)
L_equi = 13.5 
mu = 0              # ELU
nbre_point = 1001 #9      # Subdivision de la largeur en nbre_point                                                        # Longueur de la travée (m)
largeur_ep_fixe = 6.64
demie_largeur_total = 5

e = np.linspace(-b, b, nbre_point)            # generation des points e, l' excentricité de l'effort
B = np.linspace(-demie_largeur_total, demie_largeur_total, nbre_point)            # Subdivision de la largeur du tablier
h_pas = B[1]-B[0]
B_extremité = B*(B>=largeur_ep_fixe/2)+ B*(B<=-largeur_ep_fixe/2)

y_s = np.zeros([1,nbre_point]) # position relative de la section
y_s[0, :] = np.linspace(0, demie_largeur_total, nbre_point)   

alpha = 1

#y = 4.13898


# Donnée sur les forces
x0 = L_equi # position suivant x (longitudinale de la charge ponctuelle)
Q = 0   # Intensité de la charge ponctuelle
q = 14.592   # Intensité de la charge repartie
x = 5.625   # position de calcul du moment Mxx 

#on calculera la surface de k sur l'interval de definition de la charge: debut_charge, fin_charge
debut_charge = -3.8       # point de debut de la charge repartie
fin_charge = 3.8          # point de fin de la charge repartie

position_force_ponctuelle = -0.58        #(demie_largeur_total/(nbre_point-1))  # toujour s'assurer que la postion de la charge fait partit avec l'une des subdivision de B

"""

# Debut du block de code qui recupere l'enveloppe de K 
#y = 3
harmonique = np.array([[1,3,5,7]])                                                       # N° harmonique   
Matrice_Km = []
for j in range(np.shape(harmonique)[1]):
    m = float(harmonique[0, j])
    theta = m*b/L_equi #m*(b/L)*(rho_p/rho_e)**(1/4) 
    psi = np.pi*e/b
    
    K = np.zeros([nbre_point, nbre_point])    # Matrice des coefficients de repartion transversal
    for i in range(np.shape(y_s)[1]):
        y = float(y_s[0, i])                                            # position de la mesure
    
        beta = np.pi*y/b
        sig = np.pi*theta
        ksi = np.pi - abs(beta-psi)
        lbda = np.pi*theta/(np.sqrt(2)*b)
        eta = 2*lbda*b
        
        # définition des fonctions
        D1 = (1-mu)*((3+mu)*np.sinh(sig)*np.cosh(sig) - (1-mu)*sig)
        D2 = (1-mu)*((3+mu)*np.sinh(sig)*np.cosh(sig) + (1-mu)*sig)
        def F(u):
            return ((1-mu)*sig*np.cosh(sig) - (1+mu)*np.sinh(sig))*np.cosh(theta*u) - (1-mu)*theta*u*np.sinh(sig)*np.sinh(theta*u)
        def G(u):
            return (2*np.sinh(sig) + (1-mu)*sig*np.cosh(sig))*np.sinh(theta*u) - (1-mu)*theta*u*np.sinh(sig)*np.cosh(theta*u)
        
        
        km_1 = (sig/(2*sh(sig)**2))*((sig*ch(sig)+sh(sig))*ch(theta*ksi)  -  theta*ksi*sh(sig)*sh(theta*ksi) + F(beta)*F(psi)/D1 + G(beta)*G(psi)/D2)
        # km_alph = km_1 car le tablier est une dalle
    
        K[:, i] = km_1
    Matrice_Km.append(K)
    
    
"""
# Position y; ymax; ymin
#data_K = np.zeros([nbre_point, 3])
#data_K[:, 0:1] = y_s.T
#data_K[:, 1] = np.max(K*(K>=0),0) #ymax
#data_K[:, 2] = np.min(K*(K<=0),0) #ymin
"""
# Harmonique m = 1
data_K1 = np.zeros([nbre_point, 3])
data_K1[:, 0:1] = y_s.T
data_K1[:, 1] = np.max(Matrice_Km[0]*(Matrice_Km[0]>=0),0) #ymax
data_K1[:, 2] = np.min(Matrice_Km[0]*(Matrice_Km[0]<=0),0) #ymin

# Harmonique m = 3
data_K3 = np.zeros([nbre_point, 3])
data_K3[:, 0:1] = y_s.T
data_K3[:, 1] = np.max(Matrice_Km[1]*(Matrice_Km[1]>=0),0) #ymax
data_K3[:, 2] = np.min(Matrice_Km[1]*(Matrice_Km[1]<=0),0) #ymin

# Harmonique m = 5
data_K5 = np.zeros([nbre_point, 3])
data_K5[:, 0:1] = y_s.T
data_K5[:, 1] = np.max(Matrice_Km[2]*(Matrice_Km[2]>=0),0) #ymax
data_K5[:, 2] = np.min(Matrice_Km[2]*(Matrice_Km[2]<=0),0) #ymin

# Harmonique m = 7
data_K7 = np.zeros([nbre_point, 3])
data_K7[:, 0:1] = y_s.T
data_K7[:, 1] = np.max(Matrice_Km[3]*(Matrice_Km[3]>=0),0) #ymax
data_K7[:, 2] = np.min(Matrice_Km[3]*(Matrice_Km[3]<=0),0) #ymin

# Fin du block de code qui recupere l'enveloppe de K

"""


# Debut de l' exploitation de K
 
data_coef_longi_K_pour_harmonique_m = [] # le premier harmonique correspont a l'indice zero
harmonique = np.arange(50)+1                                                       # N° harmonique   

for j in range(len(harmonique)):
    m = float(harmonique[j])

    theta = m*b/L_equi #m*(b/L)*(rho_p/rho_e)**(1/4) 
    psi = np.pi*e/b
    
    #K = np.zeros([nbre_point, nbre_point])    # Matrice des coefficients de repartion transversal
    data_coef_longi_K = np.zeros([nbre_point, 4]) # remplacer 1 par la taille de la matrice y
    for i in range(np.shape(y_s)[1]):
        y = float(y_s[0, i]) 
        #y = float(y_s[0, i])                                            # position de la mesure
        
        beta = np.pi*y/b
        sig = np.pi*theta
        ksi = np.pi - abs(beta-psi)
        lbda = np.pi*theta/(np.sqrt(2)*b)
        eta = 2*lbda*b
        
        # définition des fonctions
        D1 = (1-mu)*((3+mu)*np.sinh(sig)*np.cosh(sig) - (1-mu)*sig)
        D2 = (1-mu)*((3+mu)*np.sinh(sig)*np.cosh(sig) + (1-mu)*sig)
        def F(u):
            return ((1-mu)*sig*np.cosh(sig) - (1+mu)*np.sinh(sig))*np.cosh(theta*u) - (1-mu)*theta*u*np.sinh(sig)*np.sinh(theta*u)
        def G(u):
            return (2*np.sinh(sig) + (1-mu)*sig*np.cosh(sig))*np.sinh(theta*u) - (1-mu)*theta*u*np.sinh(sig)*np.cosh(theta*u)
        
        km_1 = (sig/(2*sh(sig)**2))*((sig*ch(sig)+sh(sig))*ch(theta*ksi)  -  theta*ksi*sh(sig)*sh(theta*ksi) + F(beta)*F(psi)/D1 + G(beta)*G(psi)/D2)
        # km_alph = km_1 car le tablier est une dalle
        
        kcorige = km_1*(B>=-largeur_ep_fixe/2)*(B<=largeur_ep_fixe/2)
        indice_non_nul = np.nonzero(kcorige)[0]
        # coté des b negatif
        k1n = float(kcorige[indice_non_nul[0]])
        k2n = float(kcorige[indice_non_nul[1]])
        B1n = float(B[indice_non_nul[0]])
        B2n = float(B[indice_non_nul[1]])
        an = (k1n-k2n)/(B1n-B2n)
        bn = k1n-an*B1n
        
        # coté des b positif
        k1p = float(kcorige[indice_non_nul[-2]])
        k2p = float(kcorige[indice_non_nul[-1]])
        B1p = float(B[indice_non_nul[-2]])
        B2p = float(B[indice_non_nul[-1]])
        ap = (k1p-k2p)/(B1p-B2p)
        bp = k1p-ap*B1p
        K_extremité = (an*B_extremité+bn)*(B<-largeur_ep_fixe/2) + (ap*B_extremité+bp)*(B>largeur_ep_fixe/2)
          
        km_alph = kcorige + K_extremité   # Valeur de K interpolé sur toute la largeur du tablier  
        #on plot B (largeur du tablier, identique à l'excentricité) vs km_alph
        
        K_charge = km_alph*(B>=debut_charge)*(B<=fin_charge)
        # decomposition de K en valeur positive et negative
        K_charge_positif = K_charge*(K_charge>=0)
        K_charge_negatif = K_charge*(K_charge<=0)
        B_charge = B*(B>=debut_charge)*(B<=fin_charge)
        
        if np.any(K_charge_positif>0) :
            indice_non_nul_charge_positive = np.nonzero(K_charge_positif)[0]
            K_charge_extremité_gauche_positif = float(K_charge_positif[indice_non_nul_charge_positive[0]])
            K_charge_extremité_droit_positif = float(K_charge_positif[indice_non_nul_charge_positive[-1]])
        
        if np.any(K_charge_negatif<0) :
            indice_non_nul_charge_negatif = np.nonzero(K_charge_negatif)[0]
            K_charge_extremité_gauche_negatif = float(K_charge_negatif[indice_non_nul_charge_negatif[0]])
            K_charge_extremité_droit_negatif = float(K_charge_negatif[indice_non_nul_charge_negatif[-1]])
        
        # Moments [y, Omega+, Omega-, k_pour une force ponctuelle]
        #indice_de_stockage =  i
        data_coef_longi_K[i, 0] = y
        if np.any(K_charge_positif>0) :
            data_coef_longi_K[i, 1] = (np.sum(K_charge_positif*(K_charge_positif>=0),0)-(K_charge_extremité_gauche_positif+K_charge_extremité_droit_positif)/2)*h_pas
        if np.any(K_charge_negatif<0) :
            data_coef_longi_K[i, 2] = (np.sum(K_charge_negatif*(K_charge_negatif<=0),0)-(K_charge_extremité_gauche_negatif+K_charge_extremité_droit_negatif)/2)*h_pas
        
        indice_force_ponctuelle = np.where(np.isclose(B, position_force_ponctuelle))[0][0]
        data_coef_longi_K[i, 3] =  km_alph[indice_force_ponctuelle]
    data_coef_longi_K_pour_harmonique_m.append(data_coef_longi_K)  # pour chaque harmonique, on a [y, Omega+, Omega-, k_pour une force ponctuelle] pour toute les positions y positif de la dalle

# Debut du calcul des moments transversauts

"""
x0 = L_equi # position suivant x (longitudinale de la charge ponctuelle)
Q = 0   # Intensité de la charge ponctuelle
q = 0   # Intensité de la charge repartie
x = 0   # position de calcul du moment Mxx  
"""

Mxx_divB = np.zeros([nbre_point,4])
Mxx = np.zeros([nbre_point,4])
Mxx[:, 0] = data_coef_longi_K_pour_harmonique_m[0][:, 0]
for i in range(len(harmonique)):
    m = float(harmonique[i])
    pm_Q = (2*Q*np.sin(np.radians(m*180*x0/L_equi)))/L_equi  
    pm_q = (4*q/(m*np.pi))*(m % 2 != 0)                 # 0 si m pair
    
    km_Q = data_coef_longi_K_pour_harmonique_m[i][:, 3]
    km_q_positif = data_coef_longi_K_pour_harmonique_m[i][:, 1]
    km_q_negatif = data_coef_longi_K_pour_harmonique_m[i][:, 2]
    
    Mxm_Q =  pm_Q*km_Q*np.sin(np.radians(m*180*x/L_equi))/m**2   
    Mxm_q_positif =  pm_q*km_q_positif*np.sin(np.radians(m*180*x/L_equi))/m**2   
    Mxm_q_negatif =  pm_q*km_q_negatif*np.sin(np.radians(m*180*x/L_equi))/m**2   
    
    Mxx_divB[:, 1] = Mxx_divB[:, 1] + Mxm_q_positif
    Mxx_divB[:, 2] = Mxx_divB[:, 2] + Mxm_q_negatif
    Mxx_divB[:, 3] = Mxx_divB[:, 3] + Mxm_Q
# position y, Mxx_positif du a la charge uniforme q, Mxx_Negatif due a q, Mxx du a la charge poctuelle Q
Mxx[:, 1] = Mxx_divB[:, 1]*L_equi**2/(2*demie_largeur_total*np.pi**2) 
Mxx[:, 2] = Mxx_divB[:, 2]*L_equi**2/(2*demie_largeur_total*np.pi**2)  
Mxx[:, 3] = Mxx_divB[:, 3]*L_equi**2/(2*demie_largeur_total*np.pi**2)  
    
    

# Fin de l' exploitation de K



"""
eps = 1*(y<=e) -1*(y>e)
F1 = sh(eta)*np.cos(lbda*(b+eps*e))*ch(lbda*(b-eps*e)) - np.sin(eta)*ch(lbda*(b+eps*e))*np.cos(lbda*(b-eps*e))
F2 = ch(lbda*(b+eps*y))*np.sin(lbda*(b+eps*y)) + sh(lbda*(b+eps*y))*np.cos(lbda*(b+eps*y))
F3 = np.sin(lbda*(b+eps*e))*ch(lbda*(b-eps*e)) - np.cos(lbda*(b+eps*e))*sh(lbda*(b-eps*e))
F4 = sh(lbda*(b+eps*e))*np.cos(lbda*(b-eps*e)) - ch(lbda*(b+eps*e))*np.sin(lbda*(b-eps*e))

km_0 = (eta/(sh(eta)**2 - (np.sin(eta))**2))*(2*F1*ch(lbda*(eps*y+b))*np.cos(lbda*(eps*y+b)) + F2*(F3*sh(eta) + F4*np.sin(eta)))

f_theta = 0.05*(theta>0)*(theta<=0.1) + (1-np.exp((0.065-theta)/.663))*(theta>0.1)*(theta<=1) + 0.5*(theta>1)

km_alph = km_0 + (km_1 - km_0)*alpha**(f_theta)
"""

# Coefficient de repartition tranversale mu

data_coef_longi_mu_pour_harmonique_m = [] # le premier harmonique correspont a l'indice zero

for j in range(len(harmonique)):
    m = float(harmonique[j])

    theta = m*b/L_equi #m*(b/L)*(rho_p/rho_e)**(1/4) 
    psi = np.pi*e/b
    
    #K = np.zeros([nbre_point, nbre_point])    # Matrice des coefficients de repartion transversal
    data_coef_longi_mu = np.zeros([nbre_point, 4]) # remplacer 1 par la taille de la matrice y
    for i in range(np.shape(y_s)[1]):
        y = float(y_s[0, i])                # position de la mesure
        
        beta = np.pi*y/b
        sig = np.pi*theta
        ksi = np.pi - abs(beta-psi)
        lbda = np.pi*theta/(np.sqrt(2)*b)
        eta = 2*lbda*b
        
        # définition des fonctions
        D1 = (1-mu)*((3+mu)*np.sinh(sig)*np.cosh(sig) - (1-mu)*sig)
        D2 = (1-mu)*((3+mu)*np.sinh(sig)*np.cosh(sig) + (1-mu)*sig)
        def F(u):
            return ((1-mu)*sig*np.cosh(sig) - (1+mu)*np.sinh(sig))*np.cosh(theta*u) - (1-mu)*theta*u*np.sinh(sig)*np.sinh(theta*u)
        def G(u):
            return (2*np.sinh(sig) + (1-mu)*sig*np.cosh(sig))*np.sinh(theta*u) - (1-mu)*theta*u*np.sinh(sig)*np.cosh(theta*u)
        def H(u):
            return F(u) - 2*sh(sig)*ch(theta*u)

        def L(u):
            return G(u) - 2*sh(sig)*sh(theta*u)

        M1 = D1/(1-mu)
        M2 = D2/(1-mu)

        mu_1 = (-1/(4*sig*sh(sig)**2))*(((1-mu)*sig*ch(sig) - (1+mu)*sh(sig))*ch(theta*ksi) - (1-mu)*theta*ksi*sh(sig)*sh(theta*ksi) + H(beta)*F(psi)/M1 + L(beta)*G(psi)/M2 )
        # mu_m_alph = mu_1 car le tablier est une dalle


        # km_alph = km_1 car le tablier est une dalle
        
        mu_corige = mu_1*(B>=-largeur_ep_fixe/2)*(B<=largeur_ep_fixe/2)
        indice_non_nul = np.nonzero(mu_corige)[0]
        # coté des b negatif
        mu1n = float(mu_corige[indice_non_nul[0]])
        mu2n = float(mu_corige[indice_non_nul[1]])
        B1n = float(B[indice_non_nul[0]])
        B2n = float(B[indice_non_nul[1]])
        an = (mu1n-mu2n)/(B1n-B2n)
        bn = mu1n-an*B1n
        
        # coté des b positif
        mu1p = float(mu_corige[indice_non_nul[-2]])
        mu2p = float(mu_corige[indice_non_nul[-1]])
        B1p = float(B[indice_non_nul[-2]])
        B2p = float(B[indice_non_nul[-1]])
        ap = (mu1p-mu2p)/(B1p-B2p)
        bp = mu1p-ap*B1p
        mu_extremité = (an*B_extremité+bn)*(B<-largeur_ep_fixe/2) + (ap*B_extremité+bp)*(B>largeur_ep_fixe/2)
          
        mu_m_alph = mu_corige + mu_extremité   # Valeur de mu interpolé sur toute la largeur du tablier  
        #on plot B (largeur du tablier, identique à l'excentricité) vs mu_m_alph
        
        mu_charge = mu_m_alph*(B>=debut_charge)*(B<=fin_charge)
        # decomposition de K en valeur positive et negative
        mu_charge_positif = mu_charge*(mu_charge>=0)
        mu_charge_negatif = mu_charge*(mu_charge<=0)
        #B_charge = B*(B>=debut_charge)*(B<=fin_charge)
        
        if np.any(mu_charge_positif>0) :
            indice_non_nul_charge_positive = np.nonzero(mu_charge_positif)[0]
            mu_charge_extremité_gauche_positif = float(mu_charge_positif[indice_non_nul_charge_positive[0]])
            mu_charge_extremité_droit_positif = float(mu_charge_positif[indice_non_nul_charge_positive[-1]])
        
        if np.any(mu_charge_negatif<0) :
            indice_non_nul_charge_negatif = np.nonzero(mu_charge_negatif)[0]
            mu_charge_extremité_gauche_negatif = float(mu_charge_negatif[indice_non_nul_charge_negatif[0]])
            mu_charge_extremité_droit_negatif = float(mu_charge_negatif[indice_non_nul_charge_negatif[-1]])
        
        # Moments [y, Omega+, Omega-, k_pour une force ponctuelle]
        #indice_de_stockage =  i
        data_coef_longi_mu[i, 0] = y
        if np.any(mu_charge_positif>0) :
            data_coef_longi_mu[i, 1] = (np.sum(mu_charge_positif*(mu_charge_positif>=0),0)-(mu_charge_extremité_gauche_positif+mu_charge_extremité_droit_positif)/2)*h_pas
        if np.any(mu_charge_negatif<0) :
            data_coef_longi_mu[i, 2] = (np.sum(mu_charge_negatif*(mu_charge_negatif<=0),0)-(mu_charge_extremité_gauche_negatif+mu_charge_extremité_droit_negatif)/2)*h_pas
        
        indice_force_ponctuelle = np.where(np.isclose(B, position_force_ponctuelle))[0][0]
        data_coef_longi_mu[i, 3] =  mu_m_alph[indice_force_ponctuelle]
    data_coef_longi_mu_pour_harmonique_m.append(data_coef_longi_mu)  # pour chaque harmonique, on a [y, Omega+, Omega-, k_pour une force ponctuelle] pour toute les positions y positif de la dalle

# Debut du calcul des moments transversauts

Myy_divB = np.zeros([nbre_point,4])
Myy = np.zeros([nbre_point,4])
Myy[:, 0] = data_coef_longi_mu_pour_harmonique_m[0][:, 0]
for i in range(len(harmonique)):
    m = float(harmonique[i])
    pm_Q = (2*Q*np.sin(np.radians(m*180*x0/L_equi)))/L_equi  
    pm_q = 4*q/(m*np.pi)*(m % 2 != 0)               # 0 si m pair
    
    mu_m_Q = data_coef_longi_mu_pour_harmonique_m[i][:, 3]
    mu_m_q_positif = data_coef_longi_mu_pour_harmonique_m[i][:, 1]
    mu_m_q_negatif = data_coef_longi_mu_pour_harmonique_m[i][:, 2]
    
    Mym_Q =  pm_Q*mu_m_Q*np.sin(np.radians(m*180*x/L_equi))   
    Mym_q_positif =  pm_q*mu_m_q_positif*np.sin(np.radians(m*180*x/L_equi))   
    Mym_q_negatif =  pm_q*mu_m_q_negatif*np.sin(np.radians(m*180*x/L_equi))   
    
    Myy_divB[:, 1] = Myy_divB[:, 1] + Mym_q_positif
    Myy_divB[:, 2] = Myy_divB[:, 2] + Mym_q_negatif
    Myy_divB[:, 3] = Myy_divB[:, 3] + Mym_Q
# position y, Mxx_positif du a la charge uniforme q, Mxx_Negatif due a q, Myy du a la charge poctuelle Q
Myy[:, 1] = Myy_divB[:, 1]*demie_largeur_total 
Myy[:, 2] = Myy_divB[:, 2]*demie_largeur_total 
Myy[:, 3] = Myy_divB[:, 3]*demie_largeur_total 
    
    
# Fin de l' exploitation de mu













"""
F2p = ch(lbda*(b+eps*y))*np.sin(lbda*(b+eps*y)) - sh(lbda*(b+eps*y))*np.cos(lbda*(b+eps*y))
    
mu_0 = (1/(eta*(sh(eta)**2 - np.sin(eta)**2)))*(2*F1*sh(lbda*(eps*y+b))*np.sin(lbda*(eps*y+b)) + F2p*(F3*sh(eta) + F4*np.sin(eta)))

mu_m_alph = mu_0 + (mu_1-mu_0)*alpha**f_theta
"""


#plt.plot(e,km_1)
#plt.plot(e,mu_1)
"""
print(mu_1[8])
print(km_1[8])

"""


"""
# Graphiques
#Moment
plt.figure(figsize=(10,6))
plt.plot(data_K1[:, 0], data_K1[:, 1], color='orange', linewidth=2, linestyle='-', label='K pour m=1')
plt.plot(data_K3[:, 0], data_K3[:, 1], color='red', linewidth=2, linestyle='-', label='K pour m=3')
plt.plot(data_K5[:, 0], data_K5[:, 1], color='green', linewidth=2, linestyle='-', label='K pour m=5')
plt.plot(data_K7[:, 0], data_K7[:, 1], color='blue', linewidth=2, linestyle='-', label='K pour m=7')

plt.title("Enveloppe de K pour les Harmoniques 1,3,5, et 7" , fontsize=14)
plt.xlabel("y(m)", fontsize=12)
plt.ylabel("K", fontsize=12) 
plt.legend()
# affiche les valeur de l'axe y avec 3 decimales et impose un pas 1 sur l'axe x
ax = plt.gca() # Recupere les axes
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
plt.fill_between(data_K1[:, 0], data_K1[:, 1], color='orange',alpha=0.2)
plt.fill_between(data_K3[:, 0], data_K3[:, 1], color='red',alpha=0.2)
plt.fill_between(data_K5[:, 0], data_K5[:, 1], color='green',alpha=0.2)
plt.fill_between(data_K7[:, 0], data_K7[:, 1], color='blue',alpha=0.2)

plt.margins(0)
plt.grid()
plt.tight_layout()
#plt.savefig("Moment/LI_Moment.png", dpi=800)
plt.show()

"""

plt.figure(figsize=(10,6))
plt.plot(B, km_alph, color='blue', linewidth=2, linestyle='-', label='K pour m')

plt.title(" de K pour les " , fontsize=14)
plt.xlabel("y(m)", fontsize=12)
plt.ylabel("K", fontsize=12) 
plt.legend()
# affiche les valeur de l'axe y avec 3 decimales et impose un pas 1 sur l'axe x
ax = plt.gca() # Recupere les axes
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
plt.fill_between(B, km_alph, color='blue',alpha=0.2)

plt.margins(0)
plt.grid()
plt.tight_layout()
#plt.savefig("Moment/LI_Moment.png", dpi=800)
plt.show()



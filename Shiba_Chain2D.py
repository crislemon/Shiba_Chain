#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 12:37:27 2018

@author: cristina
"""
#everything in atomic units
import numpy as np

def Shiba_Chain2(nstep, N_atoms, state, alpha, borde, ancho, k_f, U, j, DOS, s, delta, N_omega, range_omega):
    
    # d = nstep*a distance between sites
    # N_atoms in the chain
    #state = AF or FM
    #alpha in eV*A for SOC

    "SC matrix" 
    #N_x = number of sites along the chain
    #N_y = number of sites perpendicular to the chain
    
    pi=np.pi
    N_x = N_atoms + 2*borde
    N_y = ancho
    
    #for d = 2a
    #N_x = 2*N_atoms + 2*borde
    
    "Magnetic impurities parameters and spin state"
    S = 5.0/2.0 #Cr
    S = s

    if (state == 'FM'):
        thetaS = np.zeros(N_atoms, dtype = 'float')
    elif (state == 'AF'):
        thetaS = np.zeros(N_atoms, dtype = 'float')
        for i in range(int(len(thetaS)/2)):
            thetaS[2*i+1] = pi
        if N_atoms % 2 == 0:
            thetaS[-1] = pi
    elif (state == 'inplane'):
        thetaS = np.full(N_atoms, - pi/2.0, dtype = 'float')
        

    phi = np.zeros(N_atoms)
    
    #spin spiral
    #thetaS =np.linspace(0, 2.0 * pi, N_atoms)
    
    "Material data Bi2Pd"
    Damping = 0.02/27211.6 #Dynes damping
    Delta=0.75/27211.6 #SC gap
    Delta = delta
    DOS_o = DOS #Normal phase DOS
    Fermi_k = k_f
    mass_eff=1 #SC Band effective mass
    a_interatomic=nstep*3.36/0.529

    "spin-orbit coupling"
    lamda = (alpha/(2 * 3.36))/27.2116
    #lamda = (alpha/(3.36))/27.2116

    "we define the omega vector"
    "from -N_delta/2 to N_delta/2"
    #N_omega = 2001
    N_delta = range_omega
    
    Romega = np.zeros([N_omega])
    Romega=np.array(Romega, np.longdouble)
    step_omega=N_delta*Delta/(N_omega-1)

    for i_omega in range(N_omega):
        Romega[i_omega] = (-N_delta/2.*Delta+(i_omega)*step_omega)
     
        
    Romega = np.array(Romega)
    vv=Romega*27211.6

    "Kondo hamiltonian"
    J = 1800./27211.6
    J = j
    

    "We calculate the Green's functions and solve Dyson eq"
    
    #impurity Hamiltonian
    #    import Self_Energy2D as SE
    #    Self = SE.Self_Energy(J, S, thetaS, phi, U, N_atoms, N_x, N_y, borde, lamda)
    
    import Self_Energy_loop as SL
    Self2 = SL.Self_Energy(J, S, thetaS, phi, U, N_atoms, N_x, N_y, borde, lamda)
    
    GG = np.zeros([4 * N_y * N_x , 4 * N_y * N_x, N_omega], dtype=complex)
    
    
    for i_omega in range(N_omega):
        
        omega = Romega[i_omega]
    
    
        #BCS Green's function
        import Free_Green_new as FG
        Go = FG.Free_Green(N_x, N_y, omega, Damping, Fermi_k, mass_eff, DOS_o, Delta, a_interatomic)
        
        #import Free_Gree_loop as FL
        #(Go2, go) = FL.Free_Green(N_x, N_y, omega, Damping, Fermi_k, mass_eff, DOS_o, Delta, a_interatomic)
        
        #Solve Dyson's equation
        import Dyson as Dy
        gg = Dy.Dyson_eq(Go , Self2 , N_x, N_y)
        
        GG[:,:, i_omega] = gg
        
        
        
    return(GG , N_x, N_y, N_omega , vv, Go, Self2)


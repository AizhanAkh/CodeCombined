#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 15:37:17 2018

@author: aizhan.akh
"""

def pk_nonlinear(k, P0):
    
    import numpy as np
    from numpy.fft import fft, ifft , rfft, irfft , fftfreq
    from numpy import exp, log, log10, cos, sin, pi, cosh, sinh , sqrt
    from scipy.special import gamma
    from scipy.special import hyp2f1
    from scipy import interpolate
    import sys
    from time import time
    from scipy.integrate import quad
    import scipy.integrate as integrate
    from scipy import special
    from time import time
    from scipy.optimize import fsolve
    from scipy.special import factorial
    
    #name = str(file)
    #d=np.loadtxt(file)
    #d=np.loadtxt('/Users/michalychforever/Dropbox/Docs/science/class_public-2.5.0/output/cmassDR10_z2_pk.dat', skiprows = 4)
    #k=d[:,0]; P0=d[:,1]
    Plin = interpolate.InterpolatedUnivariateSpline(k,P0)
    
    Nmax = 240
    #Nmax = len(k)
    b = -0.3
    k0 = k[0]
    kmax = k[-1]
    kcutoff = 5.
    
    ### here I prepare kbins and sample the power spectrum
    
    Delta = log(kmax/k0) / (Nmax - 1)
    jsNm = np.arange(-Nmax/2,Nmax/2+1,1)
    etam = b + 2*1j*pi*(jsNm)/Nmax/Delta
    Pdiscrin = np.zeros(Nmax);
    Plin0 = np.zeros(Nmax);
    kbins3 = np.zeros(Nmax);
    for i in range(Nmax):
    	Pdiscrin[i] = Plin(k0 * exp(Delta * i)) * exp(-1.*b*i*Delta) 
    	kbins3[i] = k0 * exp(Delta * i)
    	Plin0[i] = Plin(kbins3[i])
    
    ### here I compute FFT coeffs, throw out high-freq harmonics and symmetrize the rest 
    
    cm = np.fft.fft(Pdiscrin)/ Nmax
    cmsym = [0 for i in range(Nmax+1)]
    for i in range(0, Nmax+1):
        if (i+2 - Nmax/2) < 1:
            cmsym[i] =  k0**(-etam[i])*np.conjugate(cm[int(-i + Nmax/2)])
        else:
    	       cmsym[i] = k0**(-etam[i])* cm[int(i - Nmax/2)]
    	
    cmsym[-1] = cmsym[-1] / 2
    cmsym[0] = cmsym[0] / 2
    
    ### here I define the functions which evaluate the non-linear matrices
    
    def Gamma(x):
        return special.gamma(x)
    
    
    def J(nu1,nu2):
        return (Gamma(1.5 - nu1) * Gamma(1.5 - nu2) * Gamma(nu1 + nu2 - 1.5) / (Gamma(nu1) * Gamma(nu2) * Gamma(3. - nu2 - nu1))) / (8. * pi**(1.5))
    
    
    def M13(nu1):
        return  (1 + 9 * nu1) * np.tan(nu1*pi)/(4. * 28. * pi * (nu1 + 1.) * nu1 * (nu1 - 1.) * (nu1 - 2.) * (nu1 - 3.))
    
    
    def M22(nu1,nu2):
        return (1.5 - nu1 - nu2) * (0.5 - nu1 - nu2) * (nu1*nu2*(98.*(nu1+nu2)**2.-14.*(nu1+nu2)+36.) - 91.*(nu1+nu2)**2. + 3.*(nu1 + nu2) + 58.)*J(nu1,nu2)/(196. * nu1 * (1. + nu1) * (0.5 - nu1) * nu2 * (1. + nu2) * (0.5 - nu2))
    
    
    ### here I compute the UV part of P13
    
    sigmaV = integrate.quad(lambda x: (4*pi)*Plin(x)/(3*(2*pi)**3), k0, k[-1])[0]
    def P13UV(x):
        return -61 * 1 * Plin(x) * x**2 * sigmaV / 105
    
    
    ### here I compute 1-loop PS = multiply the non-linear matrices and Fourier harmonics
    
    P13 = np.zeros(Nmax)
    m13 = np.zeros((Nmax+1),dtype=complex)
    b13 = np.zeros((Nmax+1),dtype=complex)
    for i in range(Nmax+1):
        m13[i] = M13(-0.5 * etam[i])
    
    t03 = time()
    
    
    for j in range(Nmax):
        for i in range(Nmax+1):
            b13[i] = cmsym[i]* kbins3[j]**(etam[i])
        f13 = np.dot(b13,m13);
        P13[j] = np.real((kbins3[j]**3 * Plin(kbins3[j]) * f13 + P13UV(kbins3[j]))* exp(-(kbins3[j]/kcutoff)**6))
    
    t04 = time()
    
    #print('time', t04-t03)
    
    P22 = np.zeros(Nmax)
    P1loop = np.zeros(Nmax)
    b1 = np.zeros((Nmax+1),dtype=complex)
    m22mat = np.zeros((Nmax+1,Nmax+1),dtype=complex)
    
    for j1 in range(Nmax+1):
        for j2 in range(Nmax+1):
            if j1 - Nmax/2 < 1 :
                m22mat[j1][j2] = M22(-0.5*etam[j1],-0.5*etam[j2])		
            else:
                m22mat[j1][j2] = np.conjugate(m22mat[Nmax - j1][Nmax - j2])
    
    t1 = time()
    
    for j in range (Nmax):
        for i in range(Nmax+1):
            b1[i] = cmsym[i] * kbins3[j]**(etam[i])
        f22in = np.dot(m22mat,b1)
        f22out = np.dot(b1,f22in)
        P22[j] = np.real(kbins3[j]**3 * f22out * exp(-(kbins3[j]/kcutoff)**6))
        P1loop[j] = P13[j] + P22[j]
    
    t2=time()
    #print('time', t2-t1)
    
    
    #np.savetxt('wmap_P_1loop.dat', np.c_[kbins3,P1loop])
    #np.savetxt('wmap_P13.dat', np.c_[kbins3,P13])
    #np.savetxt('wmap_P22.dat', np.c_[kbins3,P22])
    return kbins3, P1loop

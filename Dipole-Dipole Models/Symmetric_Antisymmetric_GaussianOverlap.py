#!/bin/bash
#SBATCH --partition=long        ### Partition (like a queue in PBS)
#SBATCH --job-name=Abs_Spec_Fit      ### Job Name
#SBATCH --output=Abs.out         ### File in which to store job output
#SBATCH --error=Abs.err          ### File in which to store job error messages
#SBATCH --time=0-00:05:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=1     ### Number of tasks to be launched per Node

# You can ignore what's above this, relevant for talapas only. 
################################################################
"""
<><> SPECTRA FITTING <><>

Created on Thu Jul  5 09:37:29 2018

Fits 1D spectra for the Duplex Cy-3 Dimer.

All relevant data files for these experiments are stored in a
subdirectory under this one called "Data/" and are labeled 
in the form 

"Data/" + construct + "_Abs_visible_" + str(temp) + "C.txt"

@author: Dylan J Heussman
"""
############################################
############################################
##Andy Variables############################

phi_sm = 70              # phi, theta, r, sigma used for the single molecule
theta_sm = 20
r_sm = 7.7
sigma_sm = 316
saveimag = False


############# IMPORT STATEMENTS ##############
# Tools
#JWM: Standard import statements. with some renaming going on to reduce the name length 
#for various classes and their methods
import scipy as sp
from scipy.optimize import minimize, minimize_scalar

#JWM can import from the TQ model instead to see if there is a difference
import ext_dip_New as ed  # This file contains the definitions for 
import atomistic_tq_jan23 as tc
  # the dipole orientation angles and the
                        # extended dipole J calculation
import numpy as np
import math as ma
import matplotlib.pyplot as plt
import scipy.linalg as la

from numpy import loadtxt
from scipy.stats import norm

# Time/date stamp runs:
from datetime import datetime as dt

# Stop Warning Messages
import warnings
warnings.filterwarnings('ignore')

###############################################

########### VARIABLE DEFINITIONS ##############
Pi = np.pi

D = 3.33e-30  #Debye
mumonomer = 12.8 * D  #EDTM for Cy-3

gamma = 93 # In this version of the code 2gamma=FWHM*)
low = 16700
high = 22222 # Data ranges
c0 = 2.9979e8
H = 6.626e-34
Hbar = 1.054e-34
nubar2nu = 100*c0
permFree = 8.8542e-12

#convert J to wavenumbers value 
J2nubar = 1/(100*c0*H)


##############################################################
##############################################################
#######Variables for Single Molecule Simulation###############
#model_sm = 0             # 0 is point dipole, 1 is extended dipole 
model = 1
useTQ = True
useED = False


laser_center = 18796.99  #laser center for single molecule table
laser_sigma = 20 #/ (10**7) #width of the sm laser (arb) 

AbsPlusEnergy = []
AbsPlusAmp = []
AbsMinusEnergy = []
AbsMinusAmp = [] 


################################################

################## OPERATORS ###################
c = sp.zeros((2,2))
c[0,1] = 1
cD = c.T    # electronic raising and lowering operators

muOp = cD + c

# Vibrational Modes                                #***#   6
nVib = 6

b = sp.zeros((nVib,nVib))
for i in range(nVib-1):
    b[i,i+1] = sp.sqrt(i+1)  # vibrational raising and lowering ops
bD = b.T

Iel = sp.eye(2)
Ivib = sp.eye(nVib)  # identity ops
cDc = sp.dot(cD,c)
bDb = sp.dot(bD,b)

#################################################

################# LINEAR ALGEBRA ################
def kr(a,b): return sp.kron(a,b)
def kr4(a,b,c,d): return kr(kr(kr(a,b),c),d)
def dot(a,b): return np.dot(a,b)    

#JWM in this version of a function defintion, the if statement causses the returned
#output to default to the original vector if its normalization factor is zero, otherwise
#proceed with the divison as usual. The indent indicates which contents belong with
#which case.
def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0: 
       return v
    return v / norm

#################################################

########## BROADENING DISTRIBUTIONS #############

#Pseudo Voigt Distribution
    
# (to be placed on top of stick spectra)
#JWM: going to cross verify that the definition used in the calculator script match the 
#definitions used in the "full on" TQ model script

#JWM Same def for NormDist in both codes
def NormalDist(x, mu, sig):
    return ma.exp(-(x - mu)**2/(2 * sig**2))/(ma.sqrt(2*Pi) * sig)

#JWM Same def for CauchyDist in both Codes 
def CauchyDist(x, mu, gam):
    return 1/(Pi * gam * (((x - mu)/gam)**2 + 1))
### HACKED FOR GAUSSIAN
#JWM Q; Where are these constants all coming from? Why?
#JWM Same defition for pseduo voigt in both codes 
def PseudoVoigtDistribution(x, gamma, sigma, epsilon):
    g = ((gamma**5) + (sigma**5) + (2.69296)*(sigma**4)*(gamma) + \
        2.42843*(sigma**3)*(gamma**2) + 4.47163*(sigma**2)*(gamma**3) + \
                0.07842*(sigma)*(gamma**4))**(1/5)
    eta = gamma/g
    eta = eta * (1.36603 - 0.47719*eta + 0.11116*eta**2)
    return eta * CauchyDist(x, epsilon, g) + (1 - eta) * NormalDist(0, epsilon, g)

#DON'T USE 
#def PseudoVoigtDistribution(x, gamma, sigma, epsilon):
#    g = ((gamma**5) + (sigma**5) + (2.69296)*(sigma**4)*(gamma) + \
#        2.42843*(sigma**3)*(gamma**2) + 4.47163*(sigma**2)*(gamma**3) + \
#                0.07842*(sigma)*(gamma**4))**(1/5)
#    eta = gamma/g
#    eta = eta * (1.36603 - 0.47719*eta + 0.11116*eta**2)    
#    return NormalDist(x,epsilon,sigma)
##    return eta * CauchyDist(x, epsilon, sigma)
    

def SimData(stickdata, data, gamma, sigma, norm):
    #JWM might be important to consider the parameters in the psedoVoigt
    #broadening profile for the SM data, as we would ideally take the 
    #limit of the inhomogenous down to zero for a "single" conformer -
    #also true that we assume these macrostates have some width, just less
    # than the broadening of the full ensemble
    '''simulates a broadened spectra from a stick spectra
    array, array, num, num, num -> array'''
    output = []
    for point in data:
        simval = 0
        for stick in stickdata:
            simval += stick[1]/norm * \
                PseudoVoigtDistribution(point[0], gamma, sigma, stick[0])
                #CauchyDist(point[0],stick[0],gamma)
                #NormalDist(point[0],stick[0], sigma)
        output.append([point[0], simval])
    return np.array(output)

###################################################

################## CALCULATION ####################

def Coupling(mu1, mu2, R12):
    '''calculate coupling strength for dimer system
    #JWM This looks like a simple defition of the point dipole coupling for two 
    vectors with some seperation and direction. Used anymore?
    
    vec, vec, sep -> J (joules)'''
    return (1/(4*Pi*permFree*(la.norm(R12))**3))*(np.dot(mu1, mu2)- \
                                            3*(np.dot(R12,mu1)*(np.dot(mu2,R12)/(la.norm(R12))**2)))

# def TargetAbsCDSpectra(phiN_thetaN_R12ang_sigma_chiAbs_chiCD):
#     '''determines difference between a calculated spectra 
#     and a measured one.
    
#     6 element array-like -> num (chi^2 for spectra)'''
    
#     # first we unpack the input parameter array
    
#     phiN = phiN_thetaN_R12ang_sigma_chiAbs_chiCD[0]
#             # Twist
#     thetaN = phiN_thetaN_R12ang_sigma_chiAbs_chiCD[1]
#             # Tilt
#     R12ang = phiN_thetaN_R12ang_sigma_chiAbs_chiCD[2]
#             # Separation
#     sigma = phiN_thetaN_R12ang_sigma_chiAbs_chiCD[3]
#             # broadening
#     chiAbs = phiN_thetaN_R12ang_sigma_chiAbs_chiCD[4]
#     chiCD = phiN_thetaN_R12ang_sigma_chiAbs_chiCD[5]
#             # scale factors

    
#     phi = phiN * (Pi/180)       #
#     theta = thetaN * (Pi/180)   # some quick conversions
#     lam = ma.sqrt(lambdaSq)     #Huang Rhys param
#     R12 = R12ang*10**(-10)      #
    
#     # vector definitions in ext_dip.py
    
#     muA, muB = ed.dip_vecs(phi, theta)
    
#     # here we generate the vector dependent operators:
    
#     muTot = np.array([muA[i]*kr4(muOp, Iel, Ivib, Ivib) + \
#              muB[i]*kr4(Iel, muOp, Ivib, Ivib) for i in range(3)])
    
#     unitR = np.array([0, 0, 1])
#     Rvec = R12*unitR
    
#     magVecA = np.cross(unitR, muA)
#     magVecB = np.cross(-unitR, muB) 
    
#     magA = [magVecA[i]*muOp for i in (0,1,2)]
#     magB = [magVecB[i]*muOp for i in (0,1,2)]
    
#     op1 = kr4(magA[0], Iel, Ivib, Ivib) + kr4(Iel, magB[0], Ivib, Ivib)
#     op2 = kr4(magA[1], Iel, Ivib, Ivib) + kr4(Iel, magB[1], Ivib, Ivib)
#     op3 = kr4(magA[2], Iel, Ivib, Ivib) + kr4(Iel, magB[2], Ivib, Ivib)
    
      
#     magOps = [op1, op2, op3]
    
#     mVecA = muA*mumonomer 
#     mVecB = muB*mumonomer
    
#     # model, a global parameter, 0 for point dipole
#     #                            1 for exten dipole

#     if model == 0:
        
#         #JWM confirmed that the above "Coupling" model for simpel point dipole, appears
#         #to be a different definition than the one used in the TQ code - 
#         #might warrant a comparsion of the plots when calling one set of function defs
#         #versus another. 
        
#         # Point Coupling
#         J = J2nubar*Coupling(mVecA, mVecB, Rvec)
#         if make_plots:
#             print("Point: ", J)           # debugging
    
#     if model == 1:
#         # Extended
#         #l and c are constants for the vector length and charge respectively
#         j = ed.J_exten(phi, theta, R12, ed.l,ed.c)
#         J = J2nubar * j
#         if make_plots:
#             print("Extended: ", J)           # debugging
#             k = 8.99e9
    
# #    if model == 2:
# #        #Transition
# #        J = TransitionCoupling(phi, theta, Rvec)

#     RS = (H*nubar2nu*epsilon0/(4*Hbar)) * dot(np.cross(mVecA, mVecB), Rvec)
    
#     #JWM where does this constant come from? 7.659e-54??
#     Area = RS*epsilon0*nubar2nu/(7.659e-54)
#     Height = Area/(sigma * ma.sqrt(2 * Pi) * nubar2nu) / 2
    
    
#     # generate the full system hamiltonian
    
#     h1 = epsilon0*kr4(cDc, Iel, Ivib, Ivib)
#     h2 = epsilon0*kr4(Iel, cDc, Ivib, Ivib)
#     h3 = J*kr(kr(kr(cD, c) + kr(c, cD), Ivib), Ivib)
#     h4 = omega0*kr4(Iel, Iel, bDb, Ivib)
#     h5 = omega0*kr4(Iel, Iel, Ivib, bDb)
#     h6 = omega0*kr4(cDc, Iel, lam * (bD + b) + (lam**2) * Ivib, Ivib)
#     h7 = omega0*kr4(Iel, cDc, Ivib, lam * (bD + b) + (lam**2 * Ivib))
#     ham = h1 + h2 + h3 + h4 + h5 + h6 + h7

#     # Diagonalize Hamiltonian
#     eps, vecs = la.eig(ham)
        
#     idx = eps.argsort()[::-1]   
#     eps = eps[idx]
#     vecs = vecs[:,idx]
# #    print(vecs)           # debugging
#     eps = np.flip(eps, 0)
#     vecs = np.fliplr(vecs)

#     # absorbtion intensities
#     Ix = dot(muTot[0], vecs)[0]
#     Iy = dot(muTot[1], vecs)[0]
#     Iz = dot(muTot[2], vecs)[0]
#     SimI = (Ix**2 + Iy**2 + Iz**2)*(2/3)
    
#     AbsData = np.transpose([eps, SimI])
#     #print(SimI)
#     # CD intensities
#     cdk1 = - dot(muTot[0], vecs)[0] * dot(magOps[0], vecs)[0]
#     cdk2 = - dot(muTot[1], vecs)[0] * dot(magOps[1], vecs)[0]
#     cdk3 = - dot(muTot[2], vecs)[0] * dot(magOps[2], vecs)[0]
    
#     cdTot = Height*(cdk1 + cdk2 + cdk3)
#     np.set_printoptions(threshold=1000)
    
#     CDdata = np.transpose([eps, cdTot])
#     #print(cdTot)
    
#     closeabsdata = Closeup(AbsData, low, high)
#     closecddata = Closeup(CDdata, low, high)
    
#     normAbs = PseudoVoigtDistribution(18282, gamma, sigma, 18282)
    
#     # The variables cAbsSpectrum and cCDSpectrum (see below)
#     # are external to the function definition (global) and are 
#     # created from data files via the function "get_c_spectra"
#     # which is defined below.
    
#     simAbsSpectra = SimData(closeabsdata, cAbsSpectrum, gamma, sigma, normAbs)
#     simCDSpectra = SimData(closecddata, cCDSpectrum, gamma, sigma, normAbs)
    
#     # Implement internal chi fitting, to allow coarse grid searching for
#     # start values:

#     if chi_int:
#         # Abs #
#         #JWM what is this structure doing?
#         fit_abs = lambda chi: sum([(cAbsSpectrum[i][1]-(chi*simAbsSpectra[i][1]))**2 
#                                for i in range(len(simAbsSpectra))])
#         res = minimize_scalar(fit_abs, bounds = [0, 1])
#         chiAbs = res.x
#         # CD #
#         fit_cd = lambda chi: sum([(cCDSpectrum[i][1] - (chi * simCDSpectra[i][1]) )**2 
#                                for i in range(len(simCDSpectra))])
#         res = minimize_scalar(fit_cd, bounds = [0, 1])
#         chiCD = res.x 
#     #    print(chiAbs, chiCD)
        
#     simAbsSpectra[:,1] *= chiAbs        
#     simCDSpectra[:,1] *= chiCD
    
    
#     # determine the chi^2 for both Abs and CD
#     AbsRes = sum([(cAbsSpectrum[i][1] - simAbsSpectra[i][1])**2 for i in range(len(simAbsSpectra))])
#     CDRes = sum([(cCDSpectrum[i][1] - simCDSpectra[i][1])**2 for i in range(len(simCDSpectra))])

#     # For plotting: make_plots is a global Boolean flag which 
#     # dictates whether a target function calculation generates
#     # a plot.
    
#     if make_plots:
        
#         fig1 = plt.figure(1, figsize = (15,5))
        
#         # sort data for plots #
#         CDplus = []
#         Absplus = []
#         CDminus = []
#         Absminus = []
#         for i in range(len(CDdata)):
#             if CDdata[i][1]>0:
#                 CDplus.append(CDdata[i])
#                 Absplus.append(AbsData[i])
#             else:
#                 CDminus.append(CDdata[i])
#                 Absminus.append(AbsData[i])
        
#         CDplus = np.array(CDplus)
#         Absplus = np.array(Absplus)
#         CDminus = np.array(CDminus)
#         Absminus = np.array(Absminus)
        
#         ### ABS ###
#         ax1 = fig1.add_subplot(121)
        
#         # bars #
#         simAbsplus = SimData(Absplus, cAbsSpectrum, 2*gamma, sigma, normAbs)
#         simAbsminus = SimData(Absminus, cAbsSpectrum, 2*gamma, sigma, normAbs)
        
#         ax1.stem(Absplus[:,0], Absplus[:,1]*chiAbs, linefmt='r-', markerfmt='ro',basefmt=' ')
#         ax1.plot(simAbsplus[:,0], simAbsplus[:,1]*chiAbs, color = "r",linestyle="dashed")
        
#         ax1.stem(Absminus[:,0], Absminus[:,1]*chiAbs, linefmt='b-', markerfmt='bo',basefmt=' ')
#         ax1.plot(simAbsminus[:,0], simAbsminus[:,1]*chiAbs, color = "b",linestyle = "dashed")

#      #   somearr = np.linspace(simAbsminus[0,0],simAbsminus[296,0],len(simAbsminus[:,0]))

#         ax1.fill_between(simAbsplus[:,0], 8*norm.pdf(simAbsminus[:,0], laser_center, laser_sigma), 0, color = 'g')
# #Retrieval of the parameters from a given simulation         
#       #  np.savetxt('absplus.txt',(Absplus[:3,0]))
#      #   print(Absplus[:,1]*chiAbs)
#      #   np.savetxt('absplusheights.txt',(Absplus[:3,1]*chiAbs))
     
#     #    print((Absplus[0,1]*chiAbs)/(Absminus[36,1]*chiAbs))
#         print("The Antisymmetric Eigenenergies are:")
#         print(Absplus[:3,0]) 
#         for i in Absplus[:3,0]:
#             AbsPlusEnergy.append(i)
#         print("The Symmetric Eigenenergies are:")
#         print(Absminus[36:39,0])

#         for i in Absminus[36:39,0]:
#             AbsMinusEnergy.append(i)
#         print("The Antisymmetric Amplitudes are:")       
#         print(Absplus[:3,1])
#         for i in Absplus[:3,1]:
#             AbsPlusAmp.append(i)         
#         print("The Symmetric Amplitudes are:")
#         print(Absminus[36:39,1])
#         for i in Absminus[36:39,1]:
#             AbsMinusAmp.append(i)
        
        
#         #np.savetxt('absminus.txt', (Absminus[:3,0]))    
#        # print(Absminus[:,1]*chiAbs)
#     #    np.savetxt('absminusheights.txt',(Absminus[:,1]*chiAbs))
    
        
        
        
#         # spectra #
#    #     ax1.plot(cAbsSpectrum[:,0], cAbsSpectrum[:,1], color = "g",linewidth=3)
#         ax1.plot(simAbsSpectra[:,0], simAbsSpectra[:,1], color = "k",linewidth=3)

#         ax1.set_xlim((low, high))
#         ax1.set_title("Absorption")
#         ax1.set_ylabel("Intensity (arb. units)")
#         ax1.set_xlabel(r'$\overline{\nu}\ (cm^{-1})$')
#         ax1.axhline(0, color = "k")
#         plt.annotate('\u03A6 = %s\N{DEGREE SIGN} \n \u03B8 = %s\N{DEGREE SIGN} \n R = %s $\AA$' % (phi_sm, theta_sm, r_sm),
#                      xy=(0.02, 0.8), xycoords='axes fraction', fontsize=12) 

        
#         ### CD ###
#         ax2 = fig1.add_subplot(122)
        
#         # bars #
#         simCDplus = SimData(CDplus, cCDSpectrum, gamma, sigma, normAbs)
#         simCDminus = SimData(CDminus, cCDSpectrum, gamma, sigma, normAbs)
        
#         ax2.stem(CDplus[:,0], CDplus[:,1]*chiCD, linefmt='r-', markerfmt='ro',basefmt=' ')
#         ax2.plot(simCDplus[:,0], simCDplus[:,1]*chiCD, color = "r",linestyle="dashed")
        
#         ax2.stem(CDminus[:,0], CDminus[:,1]*chiCD, linefmt='b-', markerfmt='bo',basefmt=' ')
#         ax2.plot(simCDminus[:,0], simCDminus[:,1]*chiCD, color = "b",linestyle = "dashed")
        
#         #print(CDminus[:,0])

#         #print(CDplus[:,1]*chiCD)
#         #print((CDplus[0,1]*chiCD)/(CDminus[36,1]*chiCD))

#        # print(CDminus[:,1]*chiCD)
#         # spectra #
#    #     ax2.plot(cCDSpectrum[:,0], cCDSpectrum[:,1], color = "g",linewidth=3)
#         ax2.plot(simCDSpectra[:,0], simCDSpectra[:,1], color = "k",linewidth=3)
        
#         ax2.set_xlim((low, high))
#         ax2.set_title("CD")
#         ax2.set_ylabel("Intensity (arb. units)")
#         ax2.set_xlabel(r'$\overline{\nu}\ (cm^{-1})$')
#         ax2.axhline(0, color = "k")
   
#         if saveimag:
#             fig1.savefig("homogeneous_broadened_spectra/lorentzianweight_abscd.pdf_phi"+str(phi_sm)+"_theta"+str(theta_sm)+"_r"+str(r_sm)+".pdf", bbox_inches = "tight")
#         plt.show()
    
#     # For minimization:
#     #print("*", end = '')  # this debugging print statement just 
#                             # makes it so you know each function
#                             # evaluation happens.
    
    
#     # return the residues with appropriate weights:
#     #print(AbsRes, CDRes)
#    # return np.real(AbsRes * 1000 + CDRes * 0.001) # Target value!!
def TargetAbsCDSpectra(phiN_thetaN_rollN_shiftN_shearN_R12ang_sigma_chiAbs_chiCD):
    '''determines difference between a calculated spectra 
    and a measured one.
    
    6 element array-like -> num (chi^2 for spectra)'''
    
    # first we unpack the input parameter array
    
    phiN = phiN_thetaN_rollN_shiftN_shearN_R12ang_sigma_chiAbs_chiCD[0]
            # Twist
    thetaN = phiN_thetaN_rollN_shiftN_shearN_R12ang_sigma_chiAbs_chiCD[1]
            # Tilt
    rollN = phiN_thetaN_rollN_shiftN_shearN_R12ang_sigma_chiAbs_chiCD[2]
            #Shift
    shiftN = phiN_thetaN_rollN_shiftN_shearN_R12ang_sigma_chiAbs_chiCD[3]
            # Roll 
    shearN = phiN_thetaN_rollN_shiftN_shearN_R12ang_sigma_chiAbs_chiCD[4]
            #Shear
    R12ang = phiN_thetaN_rollN_shiftN_shearN_R12ang_sigma_chiAbs_chiCD[5]
            # Separation
    sigma = phiN_thetaN_rollN_shiftN_shearN_R12ang_sigma_chiAbs_chiCD[6]
            # broadening
    chiAbs = phiN_thetaN_rollN_shiftN_shearN_R12ang_sigma_chiAbs_chiCD[7]
    chiCD = phiN_thetaN_rollN_shiftN_shearN_R12ang_sigma_chiAbs_chiCD[8]
            # scale factors
    
#def TargetAbsCDSpectra(phiN_thetaN_R12ang_sigma_chiAbs_chiCD):
#    '''determines difference between a calculated spectra 
#    and a measured one.
#    
#    6 element array-like -> num (chi^2 for spectra)'''
#    
#    # first we unpack the input parameter array
#    
#    phiN = phiN_thetaN_R12ang_sigma_chiAbs_chiCD[0]
#            # Twist
#    thetaN = phiN_thetaN_R12ang_sigma_chiAbs_chiCD[1]
#            # Tilt
#    R12ang = phiN_thetaN_R12ang_sigma_chiAbs_chiCD[2]
#            # Separation
#    sigma = phiN_thetaN_R12ang_sigma_chiAbs_chiCD[3]
#            # broadening
#    chiAbs = phiN_thetaN_R12ang_sigma_chiAbs_chiCD[4]
#    chiCD = phiN_thetaN_R12ang_sigma_chiAbs_chiCD[5]
#            # scale factors


    
    phi = phiN * (Pi/180)       #
    theta = thetaN * (Pi/180)   # some quick conversions
    roll = rollN * (Pi/180)
    shift = shiftN *10**-(10)
    shear = shearN *10**-(10)
    R12 = R12ang  *10**-(10)
    lam = ma.sqrt(lambdaSq)
    # vector definitions in ext_dip.py
    
    #JWM First potential boolean toggle for swtiching between TQ and ED
    if useED:
        muA, muB = ed.dip_vecs(phi, theta)
    elif useTQ:
        muA, muB = tc.dip_vecs(phi,theta,roll,shift,shear,R12)
    # here we generate the vector dependent operators:
    
    muTot = np.array([muA[i]*kr4(muOp, Iel, Ivib, Ivib) + \
             muB[i]*kr4(Iel, muOp, Ivib, Ivib) for i in range(3)])
    
    unitR = np.array([0, 0, 1])
    Rvec = R12*unitR
    
    magVecA = np.cross(unitR, muA)
    magVecB = np.cross(-unitR, muB) 
    
    magA = [magVecA[i]*muOp for i in (0,1,2)]
    magB = [magVecB[i]*muOp for i in (0,1,2)]
    
    op1 = kr4(magA[0], Iel, Ivib, Ivib) + kr4(Iel, magB[0], Ivib, Ivib)
    op2 = kr4(magA[1], Iel, Ivib, Ivib) + kr4(Iel, magB[1], Ivib, Ivib)
    op3 = kr4(magA[2], Iel, Ivib, Ivib) + kr4(Iel, magB[2], Ivib, Ivib)
    
      
    magOps = [op1, op2, op3]
    
    mVecA = muA*mumonomer 
    mVecB = muB*mumonomer
    
    # model, a global parameter, 0 for point dipole
    #                            1 for exten dipole

    if model == 0:
        # Point Coupling
        J = ed.J_point(phi, theta, roll, shift, shear, R12)
       # J = J2nubar*Coupling(mVecA, mVecB, Rvec)
        if make_plots:
            print("Point: ", J)           # debugging
    
    # if model == 1:
    #     # Extended
    #     j = ed.J_exten(phi, theta, R12, ed.l,ed.c)
    #     J = J2nubar * j
    #     if make_plots:
    #         print("Extended: ", J)           # debugging
            
    if model == 1:
        # Extended
        # JWM here in the defintion of J, two possible calls exist. The first is
        #using the standard Ext Dip type calc, the other useing the full atomistic
        #coord map to determine the mag of J FOR THE EXT DIP (unphysical)
       
        #JWM second potential toggle between TQ and ED
        if useED:
            j = ed.J_exten(phi, theta, R12, ed.l,ed.c)
            J = J2nubar * j
        elif useTQ:
            J = tc.J_exten(phi, theta, roll, shift, shear, R12)
            
        if make_plots:
            print("Extended: ", J)         # debugging
            
    if model == 2:
        J = tc.J_trans(phi, theta, roll, shift, shear, R12)
        if make_plots:
            print("Transition: ", J)
        
    RS = (H*nubar2nu*epsilon0/(4*Hbar)) * dot(np.cross(mVecA, mVecB), Rvec)
    
    Area = RS*epsilon0*nubar2nu/(7.659e-54)
    Height = Area/(sigma * ma.sqrt(2 * Pi) * nubar2nu) / 2
    
    
    # generate the full system hamiltonian
    
    h1 = epsilon0*kr4(cDc, Iel, Ivib, Ivib)
    h2 = epsilon0*kr4(Iel, cDc, Ivib, Ivib)
    h3 = J*kr(kr(kr(cD, c) + kr(c, cD), Ivib), Ivib)
    h4 = omega0*kr4(Iel, Iel, bDb, Ivib)
    h5 = omega0*kr4(Iel, Iel, Ivib, bDb)
    h6 = omega0*kr4(cDc, Iel, lam * (bD + b) + (lam**2) * Ivib, Ivib)
    h7 = omega0*kr4(Iel, cDc, Ivib, lam * (bD + b) + (lam**2 * Ivib))
    ham = h1 + h2 + h3 + h4 + h5 + h6 + h7

    # Diagonalize Hamiltonian
    eps, vecs = la.eig(ham)
        
    idx = eps.argsort()[::-1]   
    eps = eps[idx]
    vecs = vecs[:,idx]
    #print(vecs)           # debugging
    eps = np.flip(eps, 0)
    vecs = np.fliplr(vecs)

    # absorbtion intensities
    Ix = dot(muTot[0], vecs)[0]
    Iy = dot(muTot[1], vecs)[0]
    Iz = dot(muTot[2], vecs)[0]
    SimI = (Ix**2 + Iy**2 + Iz**2)*(2/3)
    
    AbsData = np.transpose([eps, SimI])
    
    # CD intensities
    cdk1 = - dot(muTot[0], vecs)[0] * dot(magOps[0], vecs)[0]
    cdk2 = - dot(muTot[1], vecs)[0] * dot(magOps[1], vecs)[0]
    cdk3 = - dot(muTot[2], vecs)[0] * dot(magOps[2], vecs)[0]
    
    cdTot = Height*(cdk1 + cdk2 + cdk3)
    np.set_printoptions(threshold=1000)
    
    CDdata = np.transpose([eps, cdTot])

    closeabsdata = Closeup(AbsData, low, high)
    closecddata = Closeup(CDdata, low, high)
    
    normAbs = PseudoVoigtDistribution(18282, gamma, sigma, 18282)
    
    # The variables cAbsSpectrum and cCDSpectrum (see below)
    # are external to the function definition (global) and are 
    # created from data files via the function "get_c_spectra"
    # which is defined below.
    
    simAbsSpectra = SimData(closeabsdata, cAbsSpectrum, gamma, sigma, normAbs)
    simCDSpectra = SimData(closecddata, cCDSpectrum, gamma, sigma, normAbs)
    
    # Implement internal chi fitting, to allow coarse grid searching for
    # start values:

    if chi_int:
        # Abs #
        fit_abs = lambda chi: sum([(cAbsSpectrum[i][1]-(chi*simAbsSpectra[i][1]))**2 
                              for i in range(len(simAbsSpectra))])
        res = minimize_scalar(fit_abs, bounds = [0, 1])
        chiAbs = res.x
        # CD #
        fit_cd = lambda chi: sum([(cCDSpectrum[i][1] - (chi * simCDSpectra[i][1]) )**2 
                              for i in range(len(simCDSpectra))])
        res = minimize_scalar(fit_cd, bounds = [0, 1])
        chiCD = res.x 
 #       print(chiAbs, chiCD)
        
    simAbsSpectra[:,1] *= chiAbs        
    simCDSpectra[:,1] *= chiCD
    
    
    # determine the chi^2 for both Abs and CD
    AbsRes = sum([(cAbsSpectrum[i][1] - simAbsSpectra[i][1])**2 for i in range(len(simAbsSpectra))])
    CDRes = sum([(cCDSpectrum[i][1] - simCDSpectra[i][1])**2 for i in range(len(simCDSpectra))])

    # For plotting: make_plots is a global Boolean flag which 
    # dictates whether a target function calculation generates
    # a plot.
    
    if make_plots:
        
        fig1 = plt.figure(1, figsize = (15,5))
        
        # sort data for plots #
        CDplus = []
        Absplus = []
        CDminus = []
        Absminus = []
        for i in range(len(CDdata)):
            if CDdata[i][1]>0:
                CDplus.append(CDdata[i])
                Absplus.append(AbsData[i])
            else:
                CDminus.append(CDdata[i])
                Absminus.append(AbsData[i])
        
        CDplus = np.array(CDplus)
        Absplus = np.array(Absplus)
        CDminus = np.array(CDminus)
        Absminus = np.array(Absminus)
        
        ### ABS ###
        ax1 = fig1.add_subplot(121)
        
        # bars #
        simAbsplus = SimData(Absplus, cAbsSpectrum, gamma, sigma, normAbs)
        simAbsminus = SimData(Absminus, cAbsSpectrum, gamma, sigma, normAbs)

        ax1.stem(Absplus[:,0], Absplus[:,1]*chiAbs, linefmt='r-', markerfmt='ro',basefmt=' ')
        ax1.plot(simAbsplus[:,0], simAbsplus[:,1]*chiAbs, color = "r",linestyle="dashed")
        
        ax1.stem(Absminus[:,0], Absminus[:,1]*chiAbs, linefmt='b-', markerfmt='bo',basefmt=' ')
        ax1.plot(simAbsminus[:,0], simAbsminus[:,1]*chiAbs, color = "b",linestyle = "dashed")
        
        print("The Antisymmetric Eigenenergies are:")
        print(Absplus[:3,0]) 
        for i in Absplus[:3,0]:
            AbsPlusEnergy.append(i)
        print("The Symmetric Eigenenergies are:")
        print(Absminus[36:39,0])

        for i in Absminus[36:39,0]:
            AbsMinusEnergy.append(i)
        print("The Antisymmetric Amplitudes are:")       
        print(Absplus[:3,1])
        for i in Absplus[:3,1]:
            AbsPlusAmp.append(i)         
        print("The Symmetric Amplitudes are:")
        print(Absminus[36:39,1])
        for i in Absminus[36:39,1]:
            AbsMinusAmp.append(i)
        # spectra #
        ax1.plot(cAbsSpectrum[:,0], cAbsSpectrum[:,1], color = "g",linewidth=3)
        ax1.plot(simAbsSpectra[:,0], simAbsSpectra[:,1], color = "k",linewidth=3)

        ax1.set_xlim((low, high))
        ax1.set_title("Absorption")
        ax1.set_ylabel("Intensity (arb. units)")
        ax1.set_xlabel(r'$\overline{\nu}\ (cm^{-1})$')
        ax1.axhline(0, color = "k")
        
        ### CD ###
        ax2 = fig1.add_subplot(122)
        
        # bars #
        simCDplus = SimData(CDplus, cCDSpectrum, gamma, sigma, normAbs)
        simCDminus = SimData(CDminus, cCDSpectrum, gamma, sigma, normAbs)
        
        ax2.stem(CDplus[:,0], CDplus[:,1]*chiCD, linefmt='r-', markerfmt='ro',basefmt=' ')
        ax2.plot(simCDplus[:,0], simCDplus[:,1]*chiCD, color = "r",linestyle="dashed")
        
        ax2.stem(CDminus[:,0], CDminus[:,1]*chiCD, linefmt='b-', markerfmt='bo',basefmt=' ')
        ax2.plot(simCDminus[:,0], simCDminus[:,1]*chiCD, color = "b",linestyle = "dashed")
        
        # spectra #
        ax2.plot(cCDSpectrum[:,0], cCDSpectrum[:,1], color = "g",linewidth=3)
        ax2.plot(simCDSpectra[:,0], simCDSpectra[:,1], color = "k",linewidth=3)
        
        ax2.set_xlim((low, high))
        ax2.set_title("CD")
        ax2.set_ylabel("Intensity (arb. units)")
        ax2.set_xlabel(r'$\overline{\nu}\ (cm^{-1})$')
        ax2.axhline(0, color = "k")
        fig1.savefig("p15m1sym_vdw2A_test.pdf", bbox_inches = "tight")
        plt.show()
    
 #   For minimization:
    #print("*", end = '')  # this debugging print statement just 
                           #  makes it so you know each function
                            # evaluation happens.
    
    
  #   return the residues with appropriate weights:
 #   print(AbsRes, CDRes)
    return np.real(AbsRes * 10000 + CDRes * .005) # Target value!!




def Closeup(data, low, high):
    '''trims a data set down to fit in x range between low and high
    
    array, num, num -> array'''
    
    temp = []
    for dat in data:
        if (dat[0] > low) and (dat[0] < high):
            temp.append(dat)
    return np.array(temp)



def get_c_spectra(temp):
    '''generates Abs and CD spectra from a text file
    
    str (e.g. "J89"), str (e.g. "15") -> list of arrays'''
    #fileAbs = open("Duplex_Temp_Data/D12_"+str(temp)+"C_Abs_visible.txt", "r")
    #fileCD = open("Duplex_Temp_Data/D12_"+str(temp)+"C_CD_visible.txt", "r") 
    fileAbs = open("D12_"+str(temp)+"C_Abs_visible.txt", "r")
    fileCD = open("D12_"+str(temp)+"C_CD_visible.txt", "r") 
    AbsSpec = []
    CDSpec = []
    
    for line in fileAbs:
        AbsSpec.append([float(line.split()[i]) for i in range(len(line.split()))])
    fileAbs.close()
    for line in fileCD:
        CDSpec.append([float(line.split()[i]) for i in range(len(line.split()))])
    fileCD.close()
    
    return [Closeup(AbsSpec, low, high), Closeup(CDSpec, low, high)]
    
    CDSpec = CDSpec
    
    
    
    
    
################################################################
##################### SUBROUTINES #########################
#### Single Run ##
'''produces a target function value between a simulated spectra and taken data
   for a single construct and temperature'''
#JWM what is the role of Chi_int? It's turned on in the TQ code but off
#here. 
chi_int = False

model = 1

# Temp Depen Params #

sampleNumber = 1

temperatures = np.array([15, 25, 35, 45, 55, 
                        65, 75, 85])
epsilon0s = np.array([18285, 18277, 18266, 18262, 18280, 
                      18301,18323, 18309])
omega0s = np.array([1116, 1109, 1119, 1113, 1124, 
                    1103, 1072, 1091])
lambdaSqs = np.array([0.54, 0.56, 0.56, 0.56, 0.55, 
                      0.54, 0.54, 0.56])
    

temp = temperatures[sampleNumber]
epsilon0 = epsilon0s[sampleNumber]
omega0 = omega0s[sampleNumber]
lambdaSq = lambdaSqs[sampleNumber]

cAbsSpectrum, cCDSpectrum = get_c_spectra(temp)

make_plots = True
graphic_out = "duplex_e.svg"

#JWM above parameters for the "sample" arrays are same in the TQ and this
#code, but the start values below here are different betweeen the two for 
#seemingly the same single run use case. The Calc Code aims for SM, might
#explain the discrepancy. The final two params are the ChiAbs ChiCD

#JWM old start vector given the alt defintions
# if useED:
#     start = [phi_sm, theta_sm, r_sm, sigma_sm,
#        7.32614032e-01, 1.03765240e-01]
# elif useTQ:
#     #starting points from some opt fit in TQ
start = [ phi_sm,  theta_sm, 104, 0, 0, r_sm, sigma_sm, 0.077028, 0.37574886]
    #old roll value =104.27 (3rd entry in TQ case)

print(TargetAbsCDSpectra(start))

##################################################################

#######Single Molecule Visibility Calculation######################


xspace = np.linspace(17000,21000,8000)

laser_gaussian = norm.pdf(xspace, laser_center, laser_sigma)

#JWM as written here, it seems this builds up the array "l" by 
#scaling a cauchy distribution by magnitude j "AbsPlusAmp", then centering
# that dist at the energy i "AbsPlusEnergy". Then overlapped with the laser. 
#The broadeing is fixed at 
#93 here, need to understand why that parameter is fixed at this value

l = []
for i,j in zip(AbsPlusEnergy, AbsPlusAmp):
         l.append(j*CauchyDist(xspace, i , 93) * laser_gaussian)
    
m = []
for i,j in zip(AbsMinusEnergy, AbsMinusAmp):
         m.append(j*CauchyDist(xspace, i , 93) * laser_gaussian)
         
#pluses = (np.max(l[0]) + np.max(l[1]) + np.max(l[2]))
#minuses = (np.max(m[0]) + np.max(m[1]) + np.max(m[2]))

#JWM if you were to include a relative brightness contribution to the
# various excited states of the dimer, those scalars (quantum yield) 
#would be best implemented here, since this is the raw magntiude sum over
#the various electronic-vibrational states.
pluses = (np.sum(l[0]) + np.sum(l[1]) + np.sum(l[2]))
minuses = (np.sum(m[0]) + np.sum(m[1]) + np.sum(m[2]))

#JWM define the visiblity here. 
A = round(abs((pluses - minuses)/(minuses + pluses)),4)

numpy_arrays = [l, m]
result = np.max(numpy_arrays)



#ABS Plus
plt.figure(figsize=(10,8))
plt.plot(xspace,l[0]/result, color = 'r', linewidth= 1.75, linestyle = '--', label = '$1^{st}$ Anti-symmetric')
plt.plot(xspace,l[1]/result, 'r', linewidth= 1.75, linestyle = ':', label = '$2^{nd}$ Anti-symmetric') 
plt.plot(xspace,l[2]/result, 'r', linewidth= 1.75, linestyle = '-.',label = '$3^{rd}$ Anti-symmetric') 
plt.plot(xspace, l[0]/result + l[1]/result + l[2]/result, 'r', linewidth = 4, label = "Anti-symmetric Sum")
#ABS Minus
plt.plot(xspace,m[0]/result, color = 'b', linewidth= 1.75, linestyle = '--', label = '$1^{st}$ Symmetric')
plt.plot(xspace,m[1]/result, 'b', linestyle = ':', linewidth= 1.75, label = '$2^{nd}$ Symmetric') 
plt.plot(xspace,m[2]/result, 'b', linewidth= 1.75, linestyle = '-.', label = '$3^{rd}$ Symmetric')
plt.plot(xspace, m[0]/result + m[1]/result + m[2]/result, 'b', linewidth = 4, label = "Symmetric Sum")

#JWM might consider generalizing the bounds for these plots as passed arguments from
#the start 
plt.xlim(18700,18900)  
plt.legend()
plt.xlabel('Energy ($cm^{-1}$)', fontsize = 16)
plt.ylabel('Peak Overlap (au)', fontsize = 16)
plt.title('CW Laser Weighted Lorentzian Transition', fontsize = 18)

#plt.annotate('Something', xy=(0.05, 0.95), xycoords='axes fraction')
#plt.text(18000,.4,'yo')
plt.annotate(' Visibility = %s\n \u03A6 = %s\N{DEGREE SIGN} \n \u03B8 = %s\N{DEGREE SIGN} \n R = %s $\AA$' % (A, phi_sm, theta_sm, r_sm),
             xy=(0.02, 0.8), xycoords='axes fraction', fontsize=14) 

if saveimag:
    plt.savefig("cw_weighted_lorentzian/sim_anti_phi"+str(phi_sm)+"_theta"+str(theta_sm)+"_r"+str(r_sm)+".pdf")
plt.show()



#plt.text(0,0, 'Visibility ='A  '\n' 'Phi = 'phi_sm  
#         '\n' 'Theta= 'theta_sm   '\n'   'r = 'r_sm)



print("The visbility is:")

print(abs((pluses - minuses)/(minuses + pluses)))   


####################################################################
#############Make an ellipse for gif##############
####################################################################

import numpy as np
from matplotlib import pyplot as plt
from math import pi, cos, sin

u=0.       #x-position of the center
v=0.      #y-position of the center
a= pluses/result       #radius on the x-axis
b= minuses/result      #radius on the y-axis
t_rot=pi/4 #rotation angle

t = np.linspace(0, 2*pi, 100)
Ell = np.array([a*np.cos(t) , b*np.sin(t)])  
     #u,v removed to keep the same center location
R_rot = np.array([[cos(t_rot) , -sin(t_rot)],[sin(t_rot) , cos(t_rot)]])  
     #2-D rotation matrix

Ell_rot = np.zeros((2,Ell.shape[1]))
for i in range(Ell.shape[1]):
    Ell_rot[:,i] = np.dot(R_rot,Ell[:,i])


plt.figure(figsize=(10,8))
plt.plot( u+Ell[0,:] , v+Ell[1,:] )     #initial ellipse
#plt.plot( u+Ell_rot[0,:] , v+Ell_rot[1,:],'darkorange' )    #rotated ellipse
plt.grid(color='lightgray',linestyle='--')

#add some vectors 

X1 = [0, 0]
Y1 = [0, minuses/result]
plt.plot(X1,Y1, 'b', linewidth = 4)

X2 = [0, pluses/result]
Y2 = [0, 0]
plt.plot(X2,Y2, 'r', linewidth = 4)
plt.axis('equal')
plt.show()


































##################################################################

#### FITTING LOOP  ###
''' Loops through all the temps for a particular construct and 
    generates a list of optimum values to fit that spectra'''

#model = model_sm
#
#sampleNumbers = range(1)
#
##Temp dependent params:
#
#temperatures = np.array([15,25,35,45,55,65,75,85])
#
#epsilon0s = np.array([18285, 18277, 18266, 18262, 18280, 
#                      18301, 18323, 18309])
#
#omega0s = np.array([1116, 1109, 1119, 1113, 1124, 
#                    1103, 1072, 1091])
#
#lambdaSqs = np.array([0.54, 0.56, 0.56, 0.56, 0.55, 
#                      0.54, 0.54, 0.56])
#
#      
#now = str(dt.now())
#dat_run = str(now[0:10])
#
#if model == 0:
#    x = "Point"
#else:
#    x = "Exten"
#print(x)
#
#vals = open("D12_" + x +"_"+dat_run+"_gauss.txt", "w")
#
#
## Grid start points:
#
#phis = [70, 80, 90, 100]

#phi_sm_list = np.linspace(60,120,60)
#
#chi_int = False
#
#starts = []
#for phii in phi_sm_list:
#            starts.append([phii, theta_sm, r_sm, 1, .2, .05])
#for starts in starts:
#    TargetAbsCDSpectra(start)
#chi_int = False
#
## Loop through temps:
#
#for sampleNumber in sampleNumbers:
#    
#    ranges = [(0, 180),(0, 90),(2, 20), (150, 700), (0.005, 0.9), (0.001, 0.9)]
#    
#    # get params #
#    temp = temperatures[sampleNumber]
#    epsilon0 = epsilon0s[sampleNumber]
#    omega0 = omega0s[sampleNumber]
#    lambdaSq = lambdaSqs[sampleNumber]
#    
#    # get spectra #
#    cAbsSpectrum, cCDSpectrum = get_c_spectra(temp)
#    
#    # find keep starts #
#    make_plots = False
#    
#    chi_int = True
#    
#    start_funs = {}
#    for start in starts:
#        val = TargetAbsCDSpectra(start)
#        print("*", end = "")
#        start_funs[val] = start
#        
#    
#    # choose top 3 minimum starts:
#    keep_starts = []
#    for i in range(3):
#         min_key = min(start_funs)
#         keep_starts.append(start_funs[min_key])
#         del start_funs[min_key]
#    
#    chi_int = False     
#    print(keep_starts)
#    
#    # minimize #
#    targs = {}
#    
#    for start in keep_starts:
#        
#        make_plots = False
#
#        output = minimize(TargetAbsCDSpectra, start, \
#                              bounds = ranges,
#                              options = {"maxiter":55, "ftol":1.11e-9, "maxls":15})
#        print("\nTemp: " + str(temp))
#        print("Func: " + str(output.fun))
#        print("Xmin: " + str(list(output.x)))
#        targs[output.fun] = list(output.x)
#        
#    
#    min_x = targs[min(targs)]
#    vals.write(str(temp) +", "+ str(min(targs)) + ", " + str(min_x)[1:-1] + "\n")
#    make_plots = True
#    graphic_out = str(temp) + "_plot.svg"
#    
#    TargetAbsCDSpectra(min_x)
#    
#vals.close()


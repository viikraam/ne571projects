# -*- coding: utf-8 -*-
"""
Created on Sat Dec  2 20:01:26 2017

@author: vsingh10
"""

import numpy as np

nrich  = 5.0

TF = 813.6 # avg. temp of the fuel
TM = 560.0 # avg. temp of the moderator

def mat_val(void, bu):
    if void == 0 and bu == 0:
        
        diffC = np.array([1.4679E+00, 3.4556E-01])
        sigA = np.array([7.7991E-03, 7.1618E-02])
        nSigF = np.array([6.8063E-03, 1.2327E-01])
        gtrans = np.array(([0, 0], [1.7401E-02-sigA[0], 0]))
    
    if void == 0  and bu == 20:
        
        diffC = np.array([1.4647E+00, 3.4056E-01])
        sigA = np.array([8.3491E-03, 7.0372E-02])
        nSigF = np.array([5.3935E-03, 1.0485E-01])
        gtrans = np.array(([0, 0], [1.7649E-02-sigA[0], 0]))
        
    if void == 0 and bu == 40:
        
        diffC = np.array([1.4632E+00, 3.3976E-01])
        sigA = np.array([8.6973E-03, 6.0780E-02])
        nSigF = np.array([4.1358E-03, 7.8613E-02])
        gtrans = np.array(([0, 0], [1.7989E-02-sigA[0], 0]))
        
    if void == 40 and bu == 0:
        
        diffC = np.array([1.7036E+00, 4.3268E-01])
        sigA = np.array([7.3507E-03, 6.8018E-02])
        nSigF = np.array([6.4986E-03, 1.1913E-01])
        gtrans = np.array(([0, 0], [1.2232E-02-sigA[0], 0]))
    
    if void == 40  and bu == 20:
        
        diffC = np.array([1.7002E+00, 4.2242E-01])
        sigA = np.array([7.9075E-03, 7.0200E-02])
        nSigF = np.array([5.2651E-03, 1.0696E-01])
        gtrans = np.array(([0, 0], [1.2322E-02-sigA[0], 0]))
        
    if void == 40 and bu == 40:
        
        diffC = np.array([1.7001E+00, 4.2129E-01])
        sigA = np.array([8.2395E-03, 6.3531E-02])
        nSigF = np.array([4.2060E-03,  8.6471E-02])
        gtrans = np.array(([0, 0], [1.2464E-02-sigA[0], 0]))

        
    if void == 80 and bu == 0:
        
        diffC = np.array([2.0287E+00, 5.7035E-01])
        sigA = np.array([6.5707E-03, 6.3557E-02])
        nSigF = np.array([5.9675E-03, 1.1315E-01])
        gtrans = np.array(([0, 0], [7.2985E-03-sigA[0], 0]))
    
    if void == 80  and bu == 20:
        
        diffC = np.array([2.0246E+00, 5.4287E-01])
        sigA = np.array([7.0784E-03, 7.1050E-02])
        nSigF = np.array([4.9971E-03, 1.1044E-01])
        gtrans = np.array(([0, 0], [7.2782E-03-sigA[0], 0]))
        
    if void == 80 and bu == 40:
        
        diffC = np.array([2.0265E+00, 5.3788E-01])
        sigA = np.array([7.3691E-03, 6.8782E-02])
        nSigF = np.array([4.2214E-03, 9.8281E-02])
        gtrans = np.array(([0, 0], [7.2868E-03-sigA[0], 0]))
    
    return diffC, sigA, nSigF, gtrans
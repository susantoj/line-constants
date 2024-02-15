#!python3
#
# Copyright (C) 2017 Julius Susanto. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

"""
Overhead Line Constants Calculation Library

Author: Julius Susanto
Last update: January 2017

Reference: 
    [1] Dommel, H. W., "Electromagnetic Transients Program Reference Manual (EMTP Theory Book)", Chapter 4, "Overhead Transmission Lines".
    [2] Arrillaga, J., and Watson, N. R., "Computer Modelling of Electrical Power Systems", 2nd Edition, Wiley, 2005, Chapter 2.6

Functions:
    calc_L_int          Calculates internal inductance of solid or tubular conductor
    calc_GMR            Calculates geometric mean radius (GMR) of solid or tubular conductor
    carsons             Calculates Carson's earth return correction factors Rp and Xp for self or mutual terms
    calc_self_Z         Calculates self impedance term (in Ohm/km)
    calc_mutual_Z       Calculates mutual impedance term (in Ohm/km)
    calc_Dubanton_Z     Calculates Dubanton approximation for self or mutual impedance (in Ohm/km)  
    calc_Z_matrix       Calculates primitive impedance matrix
    calc_Y_matrix       Calculates primitive admittance matrix
    calc_kron_Z         Calculates Kron reduced matrix
    
Features currently not supported:
    - Skin effect calculation for internal self impedance
    - Input data validation, e.g. check that all phase conductor vectors are the same size
"""

import numpy as np

"""
Calculates internal inductance of solid or tubular conductor
Note that calculations assume uniform current distribution in the conductor, thus conductor stranding is not taken into account.

Usage:
    L_int = calc_L_int(type, r, q)
    
where   type is 'solid' or 'tube'
        r is the radius of the conductor (mm)
        q is the radius of the inner tube (mm)

Returns: 
        L_int the internal inductance of the conductor (mH/m)
"""
def calc_L_int(type, r, q):
    # Permeability of free space
    mu_0 = 4 * np.pi * 1e-7
    
    if type == 'solid':
        # Solid conductor
        L_int = mu_0 / 8 / np.pi
    else:
        # Tubular conductor
        L_int =  mu_0 / 2 / np.pi * (q**4 / (r**2 - q**2)**2 * np.log(r/q) - (3 * q **2 - r **2)/(4 * (r**2 - q**2)))
    
    return L_int

"""
Calculates geometric mean radius (GMR) of solid or tubular conductor
Note that calculations assume uniform current distribution in the conductor, thus conductor stranding is not taken into account.

Usage:
    GMR = calc_GMR(type, r, q)
    
where   type is 'solid' or 'tube'
        r is the radius of the conductor (mm)
        q is the radius of the inner tube (mm)

Returns: 
        GMR the geometric mean radius (mm)
"""
def calc_GMR(type, r, q):
    if type == 'solid':
        # Solid conductor
        GMR = r * np.exp(-0.25)
    else:
        # Tubular conductor
        GMR = r * np.exp((3 * q **2 - r **2)/(4 * (r**2 - q**2)) - q**4 / (r**2 - q**2)**2 * np.log(r/q))
    
    return GMR

"""
Calculates Carson's earth return correction factors Rp and Xp for both self and mutual terms.
The number of terms evaluated in the infinite loop is based on convergence to the desired error tolerance.

Usage:
    Rp, Xp = carsons(type, h_i, h_k, x_ik, f, rho, err_tol)
    
where   type is 'self' or 'mutual'
        h_i is the height of conductor i above ground (m)
        h_k is the height of conductor k above ground (m)
        x_ik is the horizontal distance between conductors i and k (m)
        f is the frequency (Hz)
        rho is the earth resistivity (Ohm.m)
        err_tol is the error tolerance for the calculation (default = 1e-6)

Returns: 
        Rp, Xp the Carson earth return correction factors (in Ohm/km)
"""
def carsons(type, h_i, h_k, x_ik, f, rho, err_tol=1e-6):  
    
    # Geometrical calculations
    if type == 'self':
        D = 2 * h_i
        cos_phi = 1
        sin_phi = 0
        phi = 0
    else:
        D = np.sqrt((h_i + h_k)**2 + x_ik ** 2)     # Distance between conductor i and image of conductor k
        cos_phi = (h_i + h_k) / D
        sin_phi = (x_ik) / D
        phi = np.arccos(cos_phi)
    
    # Initialise parameters
    i = 1
    err = 1
    sgn = 1
    
    # Initial values and constants for calculation
    omega = 2 * np.pi * f
    a = 4 * np.pi * np.sqrt(5) * 1e-4 * D * np.sqrt(f / rho)
    acosphi = a * cos_phi
    asinphi = a * sin_phi
    b = np.array([np.sqrt(2)/6, 1/16])
    c = np.array([0, 1.3659315])
    d = np.pi / 4 * b
    
    # First two terms of carson correction factor
    Rp = np.pi / 8 - b[0] * acosphi
    Xp = 0.5 * (0.6159315 - np.log(a)) + b[0] * acosphi
    
    # Loop through carson coefficient terms starting with i = 2
    while(err > err_tol):
        term = np.mod(i, 4)
        # Check sign for b term
        if term == 0:
            sgn = -1 * sgn
        
        # Calculate coefficients
        bi = b[i-1] * sgn / ((i+1) * (i+3))
        ci = c[i-1] + 1 / (i+1) + 1 / (i+3)
        di = np.pi / 4 * bi
        b = np.append(b, bi)
        c = np.append(c, ci)
        d = np.append(d, di)
        
        # Recursively calculate powers of acosphi and asinphi
        acosphi_prev = acosphi
        asinphi_prev = asinphi
        acosphi = (acosphi_prev * cos_phi - asinphi_prev * sin_phi) * a
        asinphi = (acosphi_prev * sin_phi + asinphi_prev * cos_phi) * a
        
        Rp_prev = Rp
        Xp_prev = Xp
        
        # First term
        if term == 0:
            Rp = Rp - bi * acosphi
            Xp = Xp + bi * acosphi
            
        # Second term
        elif term == 1:
            Rp = Rp + bi * ((ci - np.log(a)) * acosphi + phi * asinphi)
            Xp = Xp - di * acosphi
            
        # Third term
        elif term == 1:
            Rp = Rp + bi * acosphi
            Xp = Xp + bi * acosphi
            
        # Fourth term
        else:
            Rp = Rp - di * acosphi
            Xp = Xp - bi * ((ci - np.log(a)) * acosphi + phi * asinphi)
        
        i = i = 1
        err = np.sqrt((Rp - Rp_prev) **2 + (Xp - Xp_prev)**2)
        
    Rp = 4 * omega * 1e-04 * Rp
    Xp = 4 * omega * 1e-04 * Xp
    return Rp, Xp

"""
Calculates self impedance term (in Ohm/km)
NOTE: No allowance has been made for skin effects

Usage:
    self_Z = calc_self_Z(R_int, cond_type, r, q, h_i, f, rho, err_tol=1e-6)

where   R_int is the AC conductor resistance (Ohm/km)
        cond_type is the conductor type ('solid' or 'tube')
        r is the radius of the conductor (mm)
        q is the radius of the inner tube (mm)
        h_i is the height of conductor i above ground (m)
        f is the frequency (Hz)
        rho is the earth resistivity (Ohm.m)
        err_tol is the error tolerance for the calculation (default = 1e-6)
        
Returns: 
        self_Z the self impedance term of line impedance matrix (Ohm/km)
"""
def calc_self_Z(R_int, cond_type, r, q, h_i, f, rho, err_tol=1e-6):
    # Constants
    omega = 2 * np.pi * f                       # Nominal angular frequency
    mu_0 = 4 * np.pi * 1e-7                     # Permeability of free space
    
    # Calculate internal conductor reactance (in Ohm/km)
    X_int = 1000 * omega * calc_L_int(cond_type, r, q)
    
    # Calculate geometrical reactance (in Ohm/km)
    X_geo = 1000 * omega * mu_0 / 2 / np.pi * np.log(2 * h_i / r * 1000)
        
    # Calculate Carson's correction factors (in Ohm/km)
    Rp, Xp = carsons('self', h_i, 0, 0, f, rho, err_tol)
    
    self_Z = complex(R_int + Rp, X_int + X_geo + Xp)
    
    return self_Z

"""
Calculates mutual impedance term (in Ohm/km)

Usage:
    mutual_Z = calc_mutual_Z(cond_type, r, q, h_i, h_k, x_ik, f, rho, err_tol=1e-6)

where   cond_type is the conductor type ('solid' or 'tube')
        r is the radius of the conductor (mm)
        q is the radius of the inner tube (mm)
        h_i is the height of conductor i above ground (m)
        h_k is the height of conductor k above ground (m)
        x_ik is the horizontal distance between conductors i and k (m)
        f is the frequency (Hz)
        rho is the earth resistivity (Ohm.m)
        err_tol is the error tolerance for the calculation (default = 1e-6)
        
Returns: 
        mutual_Z the self impedance term of line impedance matrix (Ohm/km)
"""
def calc_mutual_Z(cond_type, r, q, h_i, h_k, x_ik, f, rho, err_tol=1e-6):
    # Constants
    omega = 2 * np.pi * f                       # Nominal angular frequency
    mu_0 = 4 * np.pi * 1e-7                     # Permeability of free space
    D = np.sqrt((h_i + h_k)**2 + x_ik ** 2)     # Distance between conductor i and image of conductor k
    d = np.sqrt((h_i - h_k)**2 + x_ik ** 2)     # Distance between conductors i and k
    
    # Calculate geometrical mutual reactance (in Ohm/km)
    X_geo = 1000 * omega * mu_0 / 2 / np.pi * np.log(D/d)
        
    # Calculate Carson's correction factors (in Ohm/km)
    Rp, Xp = carsons('mutual', h_i, h_k, x_ik, f, rho, err_tol)
    
    self_Z = complex(Rp, X_geo + Xp)
    
    return self_Z

"""
Calculates Dubanton approximation for self or mutual impedance (in Ohm/km)

Usage:
    Dubanton_Z = calc_Dubanton_Z(type, R_int, cond_type, r, q, h_i, h_k, x_ik, f, rho)

where   type is 'self' or 'mutual'
        cond_type is the conductor type ('solid' or 'tube')
        r is the radius of the conductor (mm)
        q is the radius of the inner tube (mm)
        h_i is the height of conductor i above ground (m)
        h_k is the height of conductor k above ground (m)
        x_ik is the horizontal distance between conductors i and k (m)
        f is the frequency (Hz)
        rho is the earth resistivity (Ohm.m)
        
Returns: 
        Dubanton_Z the self or mutual impedance term of line impedance matrix (Ohm/km)        
"""    
def calc_Dubanton_Z(type, R_int, cond_type, r, q, h_i, h_k, x_ik, f, rho):
    # Constants
    omega = 2 * np.pi * f           # Nominal angular frequency
    mu_0 = 4 * np.pi * 1e-7         # Permeability of free space
    p = np.sqrt(rho/omega/mu_0)     # Complex depth below earth
    
    if type == 'self':
        # Self impedance        
        # Calculate internal conductor reactance (in Ohm/km)
        X_int = 1000 * omega * calc_L_int(cond_type, r, q)
        
        # Calculate geometrical reactance (in Ohm/km)
        X_geo = 1000 * omega * mu_0 / 2 / np.pi * np.log(2000 * (h_i + p) / r)
        
        Dubanton_Z = complex(R_int, X_int + X_geo)
        
    else:
        # Mutual impedance
        d = np.sqrt((h_i - h_k)**2 + x_ik ** 2)     # Distance between conductors i and k
        X_geo = 1000 * omega * mu_0 / 2 / np.pi * np.log(np.sqrt((h_i + h_k + 2 * p)**2 + x_ik ** 2) / d)
        
        Dubanton_Z = complex(0, X_geo)
    
    return Dubanton_Z
    
"""
Calculates primitive impedance matrix
NOTE: all phase conductor vectors must be the same size. No checks are made to enforce this. Same goes for earth conductor vectors.

Usage:
    Z = calc_Z_matrix(line_dict)
    
where   line_dict is a dictionary of overhead line parameters:
            'mode' is the calculate mode ('carson' or 'dubanton')
            'f' is the nominal frequency (Hz)
            'rho' is the earth resistivity (Ohm.m)
            'err_tol' is the error tolerance for the calculation (default = 1e-6)
            'phase_h' is a vector of phase conductor heights above ground (m)
            'phase_x' is a vector of phase conductor horizontal spacings with arbitrary reference point (m)
            'phase_cond' is a vector of phase conductor types ('solid' or 'tube')
            'phase_R' is a vector of phase conductor AC resistances (Ohm/km)
            'phase_r' is a vector of phase conductor radii (mm)
            'phase_q' is a vector of phase conductor inner tube radii (mm) - use 0 for solid conductors
            'earth_h' is a vector of earth conductor heights above ground (m)
            'earth_x' is a vector of earth conductor horizontal spacings with arbitrary reference point (m)
            'earth_cond' is a vector of earth conductor types ('solid' or 'tube')
            'earth_R' is a vector of earth conductor AC resistances (Ohm/km)
            'earth_r' is a vector of earth conductor radii (mm)
            'earth_q' is a vector of earth conductor inner tube radii (mm) - use 0 for solid conductors
    
Returns:
        Z is the primitive impedance matrix (with earth conductors shown first)
        n_p is the number of phase conductors
        n_e is the number of earth conductors
"""    
def calc_Z_matrix(line_dict):
    # Unpack line dictionary
    mode = line_dict['mode']
    f = line_dict['f']
    rho = line_dict['rho']
    cond_h = line_dict['earth_h'] + line_dict['phase_h']
    cond_x = line_dict['earth_x'] + line_dict['phase_x']
    cond_type = line_dict['earth_cond'] + line_dict['phase_cond']
    cond_R = line_dict['earth_R'] + line_dict['phase_R']
    cond_r = line_dict['earth_r'] + line_dict['phase_r']
    cond_q = line_dict['earth_q'] + line_dict['phase_q'] 
    
    # Set error tolerance for carsons equations
    if 'err_tol' in line_dict:
        err_tol = line_dict['err_tol']
    else:
        err_tol = 1e-6
    
    # Number of phase and earth conductors
    n_p = len(line_dict['phase_h'])
    n_e = len(line_dict['earth_h'])
    n_c = n_p + n_e
        
    # Set up primitive Z matrix
    Z = np.mat(np.zeros((n_c, n_c)), dtype='complex')
    for i in range(n_c):
        for j in range(n_c):
            if i == j:
                Z[i,j] = calc_self_Z(cond_R[i],cond_type[i],cond_r[i],cond_q[i],cond_h[i],f,rho,err_tol)
            else:
                Z[i,j] = calc_mutual_Z(cond_type[i],cond_r[i],cond_q[i],cond_h[i],cond_h[j],cond_x[i]-cond_x[j],f,rho,err_tol)
    
    return Z, n_p, n_e

"""
Calculates primitive admittance matrix
Assumes that conductance of air is zero.
NOTE: all phase conductor vectors must be the same size. No checks are made to enforce this. Same goes for earth conductor vectors.

Usage:
    Y = calc_Y_matrix(line_dict)
    
where   line_dict is a dictionary of overhead line parameters:
            'f' is the nominal frequency (Hz)
            'rho' is the earth resistivity (Ohm.m)
            'err_tol' is the error tolerance for the calculation (default = 1e-6)
            'phase_h' is a vector of phase conductor heights above ground (m)
            'phase_x' is a vector of phase conductor horizontal spacings with arbitrary reference point (m)
            'phase_r' is a vector of phase conductor radii (mm)
            'earth_h' is a vector of earth conductor heights above ground (m)
            'earth_x' is a vector of earth conductor horizontal spacings with arbitrary reference point (m)
            'earth_r' is a vector of earth conductor radii (mm)
    
Returns:
        Y is the primitive admittance matrix (with earth conductors shown first)
        n_p is the number of phase conductors
        n_e is the number of earth conductors
"""    
def calc_Y_matrix(line_dict):
    # Unpack line dictionary
    f = line_dict['f']
    rho = line_dict['rho']
    cond_h = line_dict['earth_h'] + line_dict['phase_h']
    cond_x = line_dict['earth_x'] + line_dict['phase_x']
    cond_r = line_dict['earth_r'] + line_dict['phase_r'] 
    
    # Number of phase and earth conductors
    n_p = len(line_dict['phase_h'])
    n_e = len(line_dict['earth_h'])
    n_c = n_p + n_e
    
    # Constants
    omega = 2 * np.pi * f                       # Nominal angular frequency
    e_0 = 8.85418782 * 1e-12                    # Permittivity of free space
    
    # Set up primitive Y matrix
    Y = np.mat(np.zeros((n_c, n_c)), dtype='complex')
    # Build up potential coefficients
    for i in range(n_c):
        for j in range(n_c):
            if i == j:
                # Self potential coefficient
                Y[i,j] = (1 / 2 / np.pi / e_0) * np.log(2 * cond_h[i] / cond_r[i] * 1000)
            else:
                # Mutual potential coefficient
                D = np.sqrt((cond_h[i] + cond_h[j])**2 + (cond_x[i] - cond_x[j]) ** 2)     # Distance between conductor i and image of conductor k
                d = np.sqrt((cond_h[i] - cond_h[j])**2 + (cond_x[i] - cond_x[j]) ** 2)     # Distance between conductors i and k
                Y[i,j] = (1 / 2 / np.pi / e_0) * np.log(D/d)

    Y = 1000j * omega * Y.I
    
    return Y, n_p, n_e
    
"""
Calculates Kron reduced matrix
Reduction can be used for impedance or admittance matrix reductions

Usage:
    kron_Z = calc_kron_Z(Z, n_e)
    
where   Z is the primitive matrix (with earth conductors shown first)
        n_e is the number of earth conductors

Returns:
        kron_Z is the kron reduced matrix              
"""    
def calc_kron_Z(Z, n_e):
    # Slice up primtive matrix into constituent parts
    Zee = Z[0:n_e,0:n_e]
    Zep = Z[0:n_e,n_e:]
    Zpe = Z[n_e:,0:n_e]
    Zpp = Z[n_e:,n_e:]
    
    # Calculate Kron reduced matrix
    kron_Z = Zpp - Zpe * Zee.I * Zep
    
    return kron_Z

# Example   
if __name__ == '__main__':
    # Overhead line parameters (Single circuit tower with an overhead earth wire)
    line_dict = {
        'mode' : 'carson',
        'f' :           50,                         # Nominal frequency (Hz)
        'rho' :         100,                        # Earth resistivity (Ohm.m)
        'phase_h' :     [16, 20, 24],               # Phase conductor heights (m)
        'phase_x' :     [3, 2, 1],                  # Phase conductor x-axis coordinates (m)
        'phase_cond' :  ['tube','tube','tube'],     # Phase conductor types ('tube' or 'solid')
        'phase_R' :     [0.1199, 0.1199, 0.1199],   # Phase conductor AC resistances (Ohm/km)
        'phase_r' :     [10.895, 10.895, 10.895],   # Phase conductor radi (mm)
        'phase_q' :     [3, 3, 3],                  # Phase conductor inner tube radii (mm)
        'earth_h' :     [30],                       # Earth conductor heights (m)
        'earth_x' :     [0],                        # Earth conductor x-axis coordinates (m)
        'earth_cond' :  ['solid'],                  # Earth conductor types ('tube' or 'solid')
        'earth_R' :     [0.1199],                   # Earth conductor AC resistances (Ohm/km)
        'earth_r' :     [10.895],                   # Earth conductor radi (mm)
        'earth_q' :     [3]                         # Earth conductor inner tube radii (mm)
    }

    # Impedance matrix
    Zp, n_p, n_e = calc_Z_matrix(line_dict)
    Z = calc_kron_Z(Zp,n_e)
    print(Z)
    
    # Admittance matrix
    Yp, n_p, n_e = calc_Y_matrix(line_dict)
    Y = calc_kron_Z(Yp,n_e)
    print(Y)
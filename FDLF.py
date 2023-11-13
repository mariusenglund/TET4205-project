from functions import *
from time import time

import pandas as pd
import numpy as np
import sympy as sp
import cmath
import math
import numpy as np


##################### FAST DECOUPLED LOAD FLOW ########################


def get_B_matrices(file,line_data,bus_data):
    '''Obtaining B matrices for FDLF, B1 by deleting row and column of SLACK bus, B2 by deleting row and column for slack + pv buses '''
    # Get the number of buses
    num_buses = len(bus_data)

    df_ld = pd.read_excel(file, line_data)
    df_bd = pd.read_excel(file, bus_data)

    # Get relevant columns from line data
    f_bus = df_ld["From bus"].dropna().tolist() 
    t_bus = df_ld["To bus"].dropna().tolist() 
    X_line = df_ld["X[pu]"].tolist()

    Transformer = df_ld["Transformer(1 = yes, 0 = no)"].tolist()
    S_base_global = df_ld["S_base_global [MVA]"].tolist()
    S_base_global = float(S_base_global[0])

    V_base_global = df_ld["V_base_global [kV]"].tolist()
    V_base_global = float(V_base_global[0]) 

    for i in range(len(f_bus)):
        if Transformer[i] == 1:
            X_trafo_pu_local = df_ld["X [pu local]"].tolist()
            Trafo_p_voltage = df_ld["Primary voltage"].tolist()
            Trafo_s_voltage = df_ld["Secondary voltage"].tolist()
            Trafo_rating = df_ld["MVA Rating"].tolist()
            result = X_trafo_pu_local[i] * (Trafo_p_voltage[i] / V_base_global)**2 * (S_base_global / Trafo_rating[i])
            rounded_result = round(result, 6)
            X_line[i] = rounded_result

    # Convert bus numbers to integers
    f_bus = [int(x) for x in f_bus]
    t_bus = [int(x) for x in t_bus]

    # Extract bus data
    n_bus = df_bd["Bus"].tolist()
    voltage = df_bd["Voltage [pu]"].tolist()

    # Get the bus types and determine maximum bus number
    bus_type = df_bd["Bus type (slack=0, PV=1, PQ=2)"].tolist()   
    n = max(n_bus)

    # Initialize the Ybus matrix with zeros
    Ybus = np.zeros((n, n))


    # Define a 2x2 network matrix
    Net_matrix = np.array([[1, -1], [-1, 1]])


    for i in range(len(f_bus)):
        # Admittance for the line
        f = f_bus[i] - 1
        t = t_bus[i] - 1
        
        # Neglecting the line-reactance (and imag part) R_line[i] + 1j *
        Y_line = 1 / (X_line[i])

        # Update Ybus matrix based on line data, Y_line * Net_matrix
        Ybus[f][t] += Y_line * Net_matrix[0][1]
        Ybus[t][f] += Y_line * Net_matrix[0][1]
        Ybus[f][f] += Y_line * Net_matrix[0][0] 
        Ybus[t][t] += Y_line * Net_matrix[1][1]

    #Obtaining B matrices.

    # B1 by deleting row and column of SLACK bus
    for i in range(len(bus_type)):
        if bus_type[i] == 0:
            slack_bus_index = i                          
    B1 = np.delete(Ybus, slack_bus_index, axis=0)  # Delete the row
    B1 = np.delete(B1, slack_bus_index, axis=1)  # Delete the column

   # Find and store indices of PV buses
    PV_bus_indices = [i for i, bus_type in enumerate(bus_type) if bus_type == 1]

    # Initialize B2 as a copy of B1
    B2 = B1.copy()

    # Delete rows and columns for each PV bus
    for PV_bus_index in PV_bus_indices:
        B2 = np.delete(B2, PV_bus_index, axis=0)
        B2 = np.delete(B2, PV_bus_index, axis=1)

    return B1,B2,bus_type

def PQ_spec_FDLF(P_spec,Q_spec):
    '''Obtaining P and Q specified for only known P and Q '''
    P_spec_FDL = [P_spec[i] for i in range(len(P_spec)) if P_spec[i] != "Unknown"]
    Q_spec_FDL = [Q_spec[i] for i in range(len(P_spec)) if Q_spec[i] != "Unknown"]
    return P_spec_FDL,Q_spec_FDL
def PQ_calc_FDLF(Ybus,voltage, delta, P_spec, Q_spec):
    'Generates lists of mathematical expressions for active Power (P) and reactive power (Q) at bus locations where they are known. Unknown is being used to signify P and Q at buses where their values are yet to be determined'
    # Initiate lists
    n = len(P_spec)
    P_calc_FDL = []
    Q_calc_FDL = []

    for i in range(n):
        P_i = 0
        Q_i = 0

        for j in range(n):
            P_i += voltage[i] * voltage[j] * (sp.re(Ybus[i][j]) * sp.cos(delta[i] - delta[j]) + sp.im(Ybus[i][j]) * sp.sin(delta[i] - delta[j]))

            Q_i += voltage[i] * voltage[j] * (sp.re(Ybus[i][j]) * sp.sin(delta[i] - delta[j]) - sp.im(Ybus[i][j]) * sp.cos(delta[i] - delta[j]))

        if P_spec[i] != "Unknown":
            P_calc_FDL.append(P_i)
        if Q_spec[i] != "Unknown":
            Q_calc_FDL.append(Q_i)

    return P_calc_FDL, Q_calc_FDL

def voltages_delta_PQ_FDLF(voltage, Q_spec, P_spec,delta):
    '''
    Obtaining Voltages values excluding slackbus, also getting delta values for known P

    '''

    Voltage_P = [voltage[i] for i in range(len(P_spec)) if P_spec[i] != 'Unknown']
    Voltage_Q = [voltage[i] for i in range(len(Q_spec)) if Q_spec[i] != 'Unknown']
    
    delta_for_known_P = []
    for i in range(len(P_spec)):
        if P_spec[i] != 'Unknown':
            delta_for_known_P.append(delta[i])
    return Voltage_P, Voltage_Q,delta_for_known_P


def flat_start_FDLF(Voltage_P, Voltage_Q,delta_for_known_P):
    """
    Provides a flat start for Voltage_P, Voltage_Q, and delta, including known P values.
    
    """
    Voltage_P_fs = [1.0 if str(symbol).startswith('V') else symbol for symbol in Voltage_P]
    Voltage_Q_fs = [1.0 if str(symbol).startswith('V') else symbol for symbol in Voltage_Q]
    delta_fs = [0.0 if symbol.name[0] == 'd' else symbol for symbol in delta_for_known_P]
    
    return Voltage_P_fs, Voltage_Q_fs, delta_fs
def fs_voltage_P_unknowns_only(Voltage_P):
    """
    Provides a flat start for all unknown voltage magnitudes, excluding known values.
    
    """

    Voltage_P = [symbol for symbol in Voltage_P if isinstance(symbol, sp.Symbol)]
    fs_voltage_P_unknowns = []
    for symbol in Voltage_P:
        if symbol.name[0] == 'V':
            fs_voltage_P_unknowns.append(1.0)

    return fs_voltage_P_unknowns

def update_vd_i(delta_fs, fs_voltage_P):
    '''Updates values from mismatches and combines them into a single vd_i vector.'''

    vd_i = np.concatenate((delta_fs,fs_voltage_P))

    return vd_i
def update_Voltage_P_fs(voltage, Voltage_P, Voltage_P_fs, voltage_mismatch):
    '''Updating the list of Voltage_P_fs values corresponding to the known P-values.'''
    
    n = len(voltage)
    
    # Create a copy of Voltage_P_fs to avoid modifying the original list
    Voltage_P_fs_updated = Voltage_P_fs.copy()
    
    # Iterate through the specified voltages and update Voltage_P_fs_updated
    for i in range(n):
        if isinstance(voltage[i], (int, float)) and voltage[i] in Voltage_P:
            Voltage_P_fs_updated[Voltage_P.index(voltage[i])] += voltage_mismatch[i]
    
    return Voltage_P_fs_updated

def Fast_Decoupled_load_flow(file,line_data,bus_data,tolerance, max_iterations,printLaTeX):
    '''
    This function performs a Fast Decoupled Load Flow analysis for a power system

    '''
    startTime = time() # Acquire the time reference at the beginning
    convergenceTime = 0 # Initialize a time reference for the point of convergence
    i = 0   # Initialize the iteration counter


    # Extract relevant data from the power system
    Ybus, f_bus, t_bus, number_of_buses,S_base,V_base = Y_bus(file,bus_data, line_data)
    voltage,delta,P_spec,Q_spec, vd_unknown = get_Voltage_deltas_PQ_spec(file, bus_data,line_data)
    
    # Separate specified P and Q values for the FDLF
    P_spec_FDLF,Q_spec_FDLF = PQ_spec_FDLF(P_spec,Q_spec)

    # Obtain B matrices and bus types
    B1,B2,bus_type = get_B_matrices(file,line_data,bus_data)
    
    # Get voltages and deltas, excluding slack bus.
    Voltage_P, Voltage_Q,delta_for_known_P = voltages_delta_PQ_FDLF(voltage, Q_spec, P_spec,delta)

    # Initialize flat_start for voltage and delta. Including already known voltages and deltas
    Voltage_P_fs, Voltage_Q_fs, delta_fs = flat_start_FDLF(Voltage_P,Voltage_Q,delta_for_known_P)
    
    # Initialize flat_start for voltages for only unknown P
    fs_voltage_P_unknowns = fs_voltage_P_unknowns_only(Voltage_P)

    # P and Q calculated equations
    P_calc_FDL,Q_calc_FDL = PQ_calc_FDLF(Ybus,voltage,delta,P_spec,Q_spec)


    # Flat start with combined delta and voltage values for unknowns
    vd_i = flat_start(vd_unknown)

    converged = False
    while not converged: 
        # Calculate numeric P and Q, and mismatches 

        P_calc_numeric = convert_symbolic_list_to_numeric(P_calc_FDL, vd_unknown, vd_i)
        P_mismatch = P_spec_FDLF - P_calc_numeric
        Q_calc_numeric = convert_symbolic_list_to_numeric(Q_calc_FDL,vd_unknown,vd_i)
        Q_mismatch = Q_spec_FDLF - Q_calc_numeric

        # Calculates delta mismatch with B1 matrix, and updates delta_fs with mismatch values
        delta_mismatch = np.linalg.inv(B1) @ (P_mismatch/Voltage_P_fs)
        delta_fs += delta_mismatch
        
        # Calculates voltage_mismatch, and updates initial values for only unknown voltages with mismatch values
        voltage_mismatch = np.linalg.inv(B2) @ (Q_mismatch/Voltage_Q_fs)
        fs_voltage_P_unknowns +=voltage_mismatch
        Voltage_Q_fs += voltage_mismatch
        
        # Update combined flat start for P and Q calc numeric
        vd_i = update_vd_i(delta_fs, fs_voltage_P_unknowns)
        Voltage_P_fs = update_Voltage_P_fs(voltage, Voltage_P, Voltage_P_fs, voltage_mismatch)
        
        # Check for convergence
        if np.amax(np.abs(P_mismatch)) < tolerance and np.amax(np.abs(Q_mismatch)) < tolerance:
            converged = True  # Set the convergence indicator to True
            convergenceTime = time() # Acquire the time reference at the point of convergence
            print(f'FDLF converged in {i} iterations and took {convergenceTime-startTime:.{3}} seconds.') # Printing a convergence message
            break  # # Exit the loop since convergence is achieved

        i += 1

        if i >= max_iterations:
            print('FDLF did not converge in', i, 'iterations')
            break


    if converged: 
        # Update voltage and angles with converged values from vd_i
        voltage_updated = voltage.copy()
        delta_updated = delta.copy()
        for i, symbol in enumerate(vd_unknown):
            if symbol.name[0] == 'd':
                delta_updated[delta.index(symbol)] = vd_i[i]
            elif symbol.name[0] == 'V':
                voltage_updated[voltage.index(symbol)] = vd_i[i]
        
        

        P_values_updated,Q_values_updated = numeric_PQ_values(Ybus,voltage_updated,delta_updated,P_spec,Q_spec)        
        I_ij,P_ij,P_ji, Q_ij, Q_ji, P_loss, Q_loss, S_ij, S_ji = line_flows(file,bus_data, line_data,np.array(voltage_updated),np.array(delta_updated))

        bus_data_pu_df = print_bus_data_pu(voltage_updated, delta_updated, P_values_updated, Q_values_updated, printLaTeX, 3)
        bus_data_df = print_bus_data(np.array(voltage_updated), delta_updated, P_values_updated, Q_values_updated, S_base, V_base, printLaTeX, 3)
        branch_data_pu_df = print_branch_data_pu(f_bus, t_bus, P_ij, Q_ij, P_ji, Q_ji, P_loss, Q_loss, printLaTeX, 3)
        branch_data_df = print_branch_data(f_bus, t_bus, P_ij, Q_ij, P_ji, Q_ji, P_loss, Q_loss, S_base, printLaTeX, 3)

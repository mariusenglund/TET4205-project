from functions import *
from time import time

import pandas as pd
import numpy as np
import sympy as sp
import cmath
import math




##################### DECOUPLED LOAD FLOW ########################

def split_variables(vd_unknown):
    """Splits a list of variables vd_unknown into two separate lists delta_unknown and V_unknown"""
    delta_unknown = []
    V_unknown = []

    # Iterate through each variable in the input list, 
    for a in vd_unknown:
        if a.name.startswith('delta_'):
            delta_unknown.append(a)
        elif a.name.startswith('V_'):
            V_unknown.append(a)

    return delta_unknown, V_unknown

## Construct submatrix J1
def symbolic_J1(P_spec, P_calc,vd_unknown):
    """J1 Jacobian Matrix for active power"""

    delta_unknown, V_unknown = split_variables(vd_unknown)
    
    n = len(P_calc)
    P_known_calc_exp = []

    # P_known_calc_exp returns the calculated expression for P and excludes Unknowns.
    for i in range(n):
        if P_calc[i] != "Unknown":
            P_known_calc_exp.append(P_calc[i])

    # Create an empty Jacobian matrix
    J = sp.zeros(n, n)

    # Calculate the partial derivative of the known active power expression with respect to delta_unknown[b]
    for a in range(n):
        for b in range(n):
            J[a, b] = sp.diff(P_known_calc_exp[a], delta_unknown[b])

    return J,P_known_calc_exp


## Construct submatrix J4
def symbolic_J4(Q_spec,Q_calc,vd_unknown):
    """J4 Jacobian Matrix for reactive power"""
    delta_unknown, V_unknown = split_variables(vd_unknown)
    
    n = len(Q_calc)
    Q_known_calc_exp = []

    # Q_known_calc_exp returns the calculated expression for P and excludes Unknowns.
    for i in range(n):
        if Q_calc[i] != "Unknown":
            Q_known_calc_exp.append(Q_calc[i])

    # Create an empty Jacobian matrix
    J = sp.zeros(n, n)

    # Calculate the partial derivative of the known reactive power expression with respect to V_unknown[b]
    for a in range(n):
        for b in range(n):
            J[a, b] = sp.diff(Q_known_calc_exp[a], V_unknown[b])

    return J,Q_known_calc_exp


def Decoupled_load_flow(file,line_data,bus_data,tolerance, max_iterations,printLaTeX=False):
    '''
    This function performs a Decoupled Load Flow analysis for a power system

    '''
    
    
    startTime = time() # Acquire the time reference at the beginning
    convergenceTime = 0 # Initialize a time reference for the point of convergence
    i = 0   # Initialize the iteration counter

    # Extract relevant data from the power system
    Ybus,f_bus, t_bus,number_of_buses, S_base, V_base= Y_bus(file,bus_data, line_data)
    bus_type, number_of_buses,p_gen, p_load, q_gen, q_load = analyze_bus_data(file, bus_data,line_data)
    voltage,delta,P_spec,Q_spec, vd_unknown = get_Voltage_deltas_PQ_spec(file, bus_data,line_data)
    P_calc, Q_calc,PQ_calc, PQ_known_calc_exp,PQ_spec = PQ_calculation(Ybus, voltage, delta, P_spec, Q_spec)
    delta_unknown, V_unknown = split_variables(vd_unknown)

    # Initialize voltage and delta values using the flat start method
    vd_i = flat_start(vd_unknown)
    
    # From vdi, split into separate flat start initial values
    s = 0
    delta_i = np.array([vd_i[s] for s in range(len(vd_i)) if vd_unknown[s] in delta_unknown])
    voltage_i = np.array([vd_i[s] for s in range(len(vd_i)) if vd_unknown[s] in V_unknown])


    n = len(P_spec)
    # Symbolic Jacobians for active and reactive power
    J1,P_known_calc_exp = symbolic_J1(P_spec, P_calc,vd_unknown)
    J4,Q_known_calc_exp = symbolic_J4(Q_spec, Q_calc,vd_unknown)

    # Extract the known P and Q values
    known_P_values = [P_spec[i] for i in range(n) if P_spec[i] != "Unknown"]
    known_Q_values = [Q_spec[i] for i in range(n) if Q_spec[i] != "Unknown"]



    converged = False  # Add an indicator for whether convergence is achieved
    while not converged:  # Continue iterations until convergence or maximum iterations reached
        P_calc_numeric = convert_symbolic_list_to_numeric(P_known_calc_exp, vd_unknown, vd_i)    # Calculate numeric P and Q
        Q_calc_numeric = convert_symbolic_list_to_numeric(Q_known_calc_exp, vd_unknown, vd_i) 
        
        # Calculate mismatch in P and Q
        P_mismatch_num = known_P_values - P_calc_numeric
        Q_mismatch_num = known_Q_values - Q_calc_numeric      # Calculate mismatch in P and Q


        if np.amax(np.abs(P_mismatch_num)) < tolerance and np.amax(np.abs(Q_mismatch_num)) < tolerance:
            converged = True  # Set the convergence indicator to True
            convergenceTime = time() # Acquire the time reference at the point of convergence
            print(f'DLF converged in {i} iterations and took {convergenceTime-startTime:.{3}} seconds.') # Printing a convergence message            
            break  # # Exit the loop since convergence is achieved

        
        # Numerical jacobian
        J1_filled = numerical_jacobian(J1, vd_unknown, vd_i)
        J4_filled = numerical_jacobian(J4, vd_unknown, vd_i)

        # Calculate voltage and angle mismatche
        delta_mismatch = np.linalg.inv(J1_filled) @ P_mismatch_num
        voltage_mismatch = np.linalg.inv(J4_filled) @ Q_mismatch_num    # Calculate voltage mismatch
        
        # Update initial values
        i += 1  # Increment the iteration counter
        voltage_i += voltage_mismatch
        delta_i += delta_mismatch
        vd_i = np.concatenate((delta_i,voltage_i))

        if i >= max_iterations:
                    print('Decoupled did not converge in', i, 'iterations')
                    break
        

    # Print the final voltage and angles only once after the loop is finished
    if converged:

        # Update voltage and v_angle with converged values from voltage_i
        voltage_updated = voltage.copy()
        delta_updated = delta.copy()
        
        for i, symbol in enumerate(vd_unknown):
            if symbol.name[0] == 'd':
                delta_updated[delta.index(symbol)] = vd_i[i]
            elif symbol.name[0] == 'V':
                voltage_updated[voltage.index(symbol)] = vd_i[i]


        # For printing of values
        P_values_updated,Q_values_updated  = numeric_PQ_values(Ybus,voltage_updated,delta_updated,P_spec,Q_spec)
        
        I_ij,P_ij,P_ji, Q_ij, Q_ji, P_loss, Q_loss, S_ij, S_ji = line_flows(file,bus_data, line_data,np.array(voltage_updated),np.array(delta_updated))
        
        
        bus_data_pu_df = print_bus_data_pu(voltage_updated, delta_updated, P_values_updated, Q_values_updated, printLaTeX, 3)
        bus_data_df = print_bus_data(np.array(voltage_updated), delta_updated, P_values_updated, Q_values_updated, S_base, V_base, printLaTeX, 3)
        branch_data_pu_df = print_branch_data_pu(f_bus, t_bus, P_ij, Q_ij, P_ji, Q_ji, P_loss, Q_loss, printLaTeX, 3)
        branch_data_df = print_branch_data(f_bus, t_bus, P_ij, Q_ij, P_ji, Q_ji, P_loss, Q_loss, S_base, printLaTeX, 3)
    
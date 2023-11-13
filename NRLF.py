from functions import *
from time import time

import pandas as pd
import numpy as np
import sympy as sp
import cmath
import math

##################### NEWTON RHAPSON LOAD FLOW ########################

def NR_load_flow(file,line_data,bus_data,tolerance, max_iterations,printLaTeX):
    '''
    This function utilizes the Newton-Raphson method to iteratively solve the power flow equations
    and converge to a solution that satisfies the specified convergence criteria. It considers both
    active (P) and reactive (Q) power mismatches to determine the convergence of the solution.
    
    '''
    
    startTime = time() # Acquire the time reference at the beginning
    convergenceTime = 0 # Initialize a time reference for the point of convergence
    i = 0   # Initialize the iteration counter
    
     # Get the Y-bus matrix and other relevant data from the power system
    Ybus,f_bus, t_bus,number_of_buses, S_base, V_base= Y_bus(file,bus_data, line_data)
    bus_type, number_of_buses,p_gen, p_load, q_gen, q_load = analyze_bus_data(file, bus_data,line_data)
    voltage,delta,P_spec,Q_spec, vd_unknown = get_Voltage_deltas_PQ_spec(file, bus_data,line_data)
    P_calc, Q_calc,PQ_calc, PQ_known_calc_exp,PQ_spec = PQ_calculation(Ybus, voltage, delta, P_spec, Q_spec)
    J = symbolic_jacobian(PQ_known_calc_exp, vd_unknown)
    vd_i = flat_start(vd_unknown)  # Initilize voltage and values uing the flat start method

    

    converged = False  # Convergence indicator for the loop


    # Continue iterations until convergence or maximum iterations reached
    while not converged:  
        PQ_calc_numeric_i = convert_symbolic_list_to_numeric(PQ_known_calc_exp, vd_unknown, vd_i)    # Calculate numeric P and Q
        PQ_mismatch_num = PQ_spec - PQ_calc_numeric_i       # Calculate mismatch in P and Q

        
        if np.amax(np.abs(PQ_mismatch_num)) < tolerance:
            converged = True  # Set the convergence indicator to True
            convergenceTime = time() # Acquire the time reference at the point of convergence
            print(f'NRLF converged in {i} iterations and took {convergenceTime-startTime:.{3}} seconds.') # Printing a convergence message            break  # # Exit the loop since convergence is achieved

        
        J_filled = numerical_jacobian(J, vd_unknown, vd_i) # Get numerical values in  the jacobian amtrix
        voltage_mismatch = np.linalg.inv(J_filled) @ PQ_mismatch_num    # Calculate voltage mismatch
        i += 1  # Increment the iteration counter
        vd_i += voltage_mismatch   # Update voltage and angles
        
        if i >= max_iterations:
                    print('NR did not converge in', i, 'iterations')
                    break
       

    # Print the final voltage and angles only once after the loop is finished
    if converged:

        # Update voltage and v_angle with converged values from vd_i
        voltage_updated = voltage.copy()
        delta_updated = delta.copy()
        for i, symbol in enumerate(vd_unknown):
            if symbol.name[0] == 'd':
                delta_updated[delta.index(symbol)] = vd_i[i]
            elif symbol.name[0] == 'V':
                voltage_updated[voltage.index(symbol)] = vd_i[i]
        
        
    
        # Calculate numeric P and Q values based on the converged solution
        P_values_updated,Q_values_updated = numeric_PQ_values(Ybus,voltage_updated,delta_updated,P_spec,Q_spec)   

        # Calculate line flows, losses, and apparent power based on the converged solution     
        I_ij,P_ij,P_ji, Q_ij, Q_ji, P_loss, Q_loss, S_ij, S_ji = line_flows(file,bus_data, line_data,np.array(voltage_updated),np.array(delta_updated))



        # Print bus data and branch data
        bus_data_pu_df = print_bus_data_pu(voltage_updated, delta_updated, P_values_updated, Q_values_updated, printLaTeX, 3)
        bus_data_df = print_bus_data(np.array(voltage_updated), delta_updated, P_values_updated, Q_values_updated, S_base, V_base, printLaTeX, 3)
        branch_data_pu_df = print_branch_data_pu(f_bus, t_bus, P_ij, Q_ij, P_ji, Q_ji, P_loss, Q_loss, printLaTeX, 3)
        branch_data_df = print_branch_data(f_bus, t_bus, P_ij, Q_ij, P_ji, Q_ji, P_loss, Q_loss, S_base, printLaTeX, 3)


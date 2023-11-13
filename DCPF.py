from functions import *

import pandas as pd
import numpy as np
import sympy as sp
import cmath
import math

##################### DC POWER FLOW ########################

def Y_bus_DC(file,bus_data, line_data):
    '''
    Obtain Ybus matrix for DC load flow, by neglecting all real parts of elements in Ybus, and removing row and column of slack bus. 
    
    '''


    # Read line data and bus data from Excel file
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

    
    # Loop through each line to calculate transformer reactance and update line reactance
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
    
    # Set all voltage values to 1, DC power flow
    voltage = [1 for _ in voltage] 
    
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

    #Deleting rows and columns of SLACK-bus
    for i in range(len(bus_type)):
        if bus_type[i] == 0: # SLACK-bus
            Ybus = np.delete(Ybus, i, axis=0)  # Delete the row
            Ybus = np.delete(Ybus, i, axis=1)  # Delete the column
    
    return Ybus, voltage,f_bus,t_bus, S_base_global,V_base_global

def generate_P_spec_DC(file,line_data,bus_data):
    '''
    Calculates specified active and reactive power from P and Q specified for DC Load flow
    
    '''
    # Calculate the Ybus matrix and recieve bus data
    Ybus_DC, voltage_DC,f_bus,t_bus,S_base,V_base  = Y_bus_DC(file,bus_data, line_data)
    bus_type, number_of_buses,p_gen, p_load, q_gen, q_load = analyze_bus_data(file,bus_data, line_data)
    voltage,delta,P_spec,Q_spec, vd_unknown = get_Voltage_deltas_PQ_spec(file, bus_data,line_data)
    # Number of buses

    n = len(voltage_DC)

    # Initialize a list to store specified power generation values
    P_spec_DC = []
    # Iterate through the buses to calculate spedified power at buses, and give 'Unknown' at slack buses
    for i in range(n):
        if bus_type[i]== 0:
            P_spec_DC.append('Unknown') # Slack bus, power generation is "Unknown"
        else: 
            P_spec_DC.append(p_gen[i]-p_load[i])
            delta[i] = sp.symbols("delta_" + str(i+1)) # Define a symbolic angle variable for the delta values
    
    # Create a list of spedificed power values with 'Unknowns' excluded
    P_spec_known = [P_spec_DC[i] for i in range(len(P_spec_DC)) if P_spec_DC[i] != "Unknown"]
    
    return P_spec_DC,P_spec_known,delta



def DC_line_flows(file, line_data,bus_data,delta_updated):

    '''
    Calculates line flows and losses based on updated voltage and delta values for DC load flow.

    '''
    # Read line and bus data from Excel files
    df_ld = pd.read_excel(file, line_data)
    df_bd = pd.read_excel(file, bus_data)  

    # List of columns in the excel file Line data
    f_bus = df_ld["From bus"].dropna().tolist() 
    t_bus = df_ld["To bus"].dropna().tolist() 
    R_line = df_ld["R[pu]"].tolist()
    X_line = df_ld["X[pu]"].tolist()
    Transformer = df_ld["Transformer(1 = yes, 0 = no)"].tolist()
    S_base_global = df_ld["S_base_global [MVA]"].tolist()
    S_base_global = float(S_base_global[0])

    V_base_global = df_ld["V_base_global [kV]"].tolist()
    V_base_global = float(V_base_global[0])

    
    # Loop through each line to calculate transformer reactance and update line reactance
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

    # Finding the number of buses and number of lines
    n_bus = df_bd["Bus"].tolist()
    n = len(f_bus)

    # Create an empty array for line impedances, which is only X_line values for DC power flow
    Z_line = np.empty(n)
    for i in range(n):
        Z_line[i] = X_line[i]

    # Initialize arrays for various line flow parameters
    P_ij = np.array([])
    P_ji = np.array([])
    Q_ij = np.array([])
    Q_ji = np.array([])
    P_loss = np.array([])
    Q_loss = np.array([])

    S_ij = np.array([])
    S_ji = np.array([])

    # Calculate line power flows and losses
    for i in range(n):
        Pij = (1 / Z_line[i]) * (delta_updated[f_bus[i]-1] - delta_updated[t_bus[i]-1])
        P_ij = np.append(P_ij, Pij)
        
        Pji = (1 / Z_line[i]) * (delta_updated[t_bus[i]-1] - delta_updated[f_bus[i]-1])
        P_ji = np.append(P_ji,Pji)

        # Initialize arrays for reactive power (Q) with zeros
        Q_ij = np.zeros(n)
        Q_ji = np.zeros(n)

        # Calculate total power losses
        P_loss = np.append(P_loss, P_ij[i] + P_ji[i])
        Q_loss = np.append(Q_loss, Q_ij[i] + Q_ji[i])

         # Define complex power (S) values
        S_ij = P_ij
        S_ji = P_ji
    # Return the calculated line power flows and losses
    return P_ij,P_ji, Q_ij, Q_ji, P_loss, Q_loss, S_ij, S_ji



def DC_power_flow(file,bus_data,line_data,printLaTeX=False):
    ''' 
    This function performs the DC load flow analysis for a power system
    By assuming all voltages =1 , and calculating angles (delta) from modified Ybus

    '''

    # Obtain values
    P_spec_DC,P_spec_known_DC,delta_unknown = generate_P_spec_DC(file,line_data,bus_data)
    Ybus_DC, voltage_DC,f_bus,t_bus,S_base,V_base  = Y_bus_DC(file,bus_data, line_data)
    
    # Calculate the unknown delta values in radians and degrees using matrix inversion  
    delta_radians = np.linalg.inv(Ybus_DC) @ P_spec_known_DC
    delta_radians = list(delta_radians)
    delta_degrees = np.degrees(delta_radians)
    delta_degrees = list(delta_degrees)
    
    # Assign calculated delta values to the unknown angles
    for i in range(len(delta_unknown)):
        if delta_unknown[i] != 0:
            delta_unknown[i] = delta_radians.pop(0)
    delta_updated=delta_unknown


   # Calculate total power at the slack bus
    P_slack = 0
    P_values_updated_DC = []
    n = len(P_spec_known_DC)
    n2 = len(P_spec_DC)

    for i in range (n):
        P_slack -= P_spec_known_DC[i]

    # Assign calculated power for slack bus to the other calculated specified power values
    for s in range (n2):
        if P_spec_DC[s] == 'Unknown':
            P_spec_DC[s] = P_slack
            P_values_updated_DC = np.array(P_spec_DC.copy())

    # Initialize reactive power values as zeros
    Q_values_updated_DC = np.zeros(n2)

    # Calculate line flows and losses
    P_ij,P_ji, Q_ij, Q_ji, P_loss, Q_loss, S_ij, S_ji = DC_line_flows(file, line_data,bus_data,delta_updated)

    # Print bus and branch data in per unit or standard units (if requested)
    bus_data_pu_df = print_bus_data_pu(voltage_DC, delta_updated, P_values_updated_DC, Q_values_updated_DC, printLaTeX, 3)
    bus_data_df = print_bus_data(np.array(voltage_DC), delta_updated, P_values_updated_DC, Q_values_updated_DC, S_base, V_base, printLaTeX, 3)
    branch_data_pu_df = print_branch_data_pu(f_bus, t_bus, P_ij, Q_ij, P_ji, Q_ji, P_loss, Q_loss, printLaTeX, 2)
    branch_data_df = print_branch_data(f_bus, t_bus, P_ij, Q_ij, P_ji, Q_ji, P_loss, Q_loss, S_base, printLaTeX, 2)
    
    
















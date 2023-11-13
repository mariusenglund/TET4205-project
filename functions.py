import pandas as pd
import numpy as np
import sympy as sp
import cmath
import math

        
def Y_bus(file,bus_data, line_data):
    '''
    This function calculates the Y-bus matrix for a power system based on line data and bus data from Excel files.

    Parameters:
        file (str): The name of the Excel file containing line and bus data.
        line_data (str): The sheet name in the Excel file with line data.

    Returns:
        np.ndarray: The Y-bus matrix of the power system.
    '''

    # Read data from Excel files
    df_ld = pd.read_excel(file, line_data)
    df_bd = pd.read_excel(file, bus_data)

    # Extract relevant data from line data sheet
    f_bus = df_ld["From bus"].dropna().tolist()
    t_bus = df_ld["To bus"].dropna().tolist()
    R_line = df_ld["R[pu]"].tolist()
    X_line = df_ld["X[pu]"].tolist()
    adm_line = df_ld["Line Charging Admittance [pu]"]
    
    Transformer = df_ld["Transformer(1 = yes, 0 = no)"].tolist()
    S_base_global = df_ld["S_base_global [MVA]"].tolist()
    S_base_global = float(S_base_global[0])

    V_base_global = df_ld["V_base_global [kV]"].tolist()
    V_base_global = float(V_base_global[0])

    
    # Convert to integers
    f_bus = [int(x) for x in f_bus]
    t_bus = [int(x) for x in t_bus]

    # Extract relevant data from bus data sheet
    n_bus = df_bd["Bus"].tolist()
    number_of_buses = max(n_bus)


    
    # Loop through each line to calculate transformer reactance and update line reactance
    for i in range(len(f_bus)):
        if Transformer[i] == 1:
            # Extract relevant data for transformer
            X_trafo_pu_local = df_ld["X [pu local]"].tolist()
            Trafo_p_voltage = df_ld["Primary voltage"].tolist()
            Trafo_s_voltage = df_ld["Secondary voltage"].tolist()
            Trafo_rating = df_ld["MVA Rating"].tolist()

            # Calculate transformer reactance in per unit
            result = X_trafo_pu_local[i] * (Trafo_p_voltage[i] / V_base_global)**2 * (S_base_global / Trafo_rating[i])
            rounded_result = round(result, 6)
            X_line[i] = rounded_result

    # Initialize Y-bus matrix with zeros
    Ybus = np.zeros((number_of_buses, number_of_buses), dtype='complex_')

    # Defining network matrix 
    Net_matrix = np.array([[1, -1], [-1, 1]])
    
    # Iterate through each branch in the system
    for i in range(len(f_bus)):
        # Get the from and to bus number, subtracting for the first name rows.
        f = f_bus[i] - 1
        t = t_bus[i] - 1

        # Calculate Y_line for the branch
        Y_line = 1 / (complex(R_line[i], X_line[i]))
        
        # Update the Y-bus matrix by adding contributions from each branch
        Ybus[f][t] += Y_line * Net_matrix[0][1]
        Ybus[t][f] += Y_line * Net_matrix[0][1]
        # Adding line charging admittance at diagonal of Ybus
        Ybus[f][f] += Y_line * Net_matrix[0][0] + complex(0, adm_line[i] / 2)
        Ybus[t][t] += Y_line * Net_matrix[1][1] + complex(0, adm_line[i] / 2)

    return Ybus, f_bus, t_bus, number_of_buses,S_base_global,V_base_global


def analyze_bus_data(file,bus_data, line_data):
    '''This function analyzes bus data from an excel file and returns relevant data'''
    # Read bus and line data from Excel sheets
    df_bd = pd.read_excel(file, bus_data) #Bus data from Excel-sheet
    df_ld = pd.read_excel(file,line_data)
    
    # Extract relevant bus data
    bus = df_bd['Bus'].tolist()
    bus_type = df_bd["Bus type (slack=0, PV=1, PQ=2)"].tolist()
    
    # Check the unit of the data (SI or Pu)
    unit = df_bd['Unit(SI or Pu)'].tolist()
    unit = str(unit[0])
    
    # If the unit is in SI, convert relevant columns to per unit values
    if unit == 'SI':

        S_base_global = df_ld["S_base_global [MVA]"].tolist()
        S_base_global = float(S_base_global[0])

        V_base_global = df_ld["V_base_global [kV]"].tolist()
        V_base_global = float(V_base_global[0])

        relevant_columns = ['P_gen [MW]', 'P_load [MW]', 'Q_load [MVar]', 'Q_gen [MVar]']
        for col in relevant_columns:
            df_bd[col] = df_bd[col] / S_base_global
            df_bd[col].fillna(0, inplace=True)

        p_gen = df_bd["P_gen [MW]"].tolist()
        p_load = df_bd["P_load [MW]"].tolist()
        q_load = df_bd["Q_load [MVar]"].tolist()
        q_gen = df_bd["Q_gen [MVar]"].tolist()

    # If the unit is in per unit (Pu), directly extract the values
    elif unit =='Pu':

        p_gen = df_bd["P_gen [pu]"].tolist()
        p_load = df_bd["P_load [pu]"].tolist()
        q_load = df_bd["Q_load [pu]"].tolist()
        q_gen = df_bd["Q_gen [pu]"].tolist()
    
    # If the unit is not recognized, print an error message
    else:
        print('Error: Unit is not recognized. Please use either "SI" or "Pu".')
    
    number_of_buses = max(bus)

   
    return bus_type, number_of_buses,p_gen, p_load, q_gen, q_load


def get_Voltage_deltas_PQ_spec(file, bus_data,line_data):
    '''This function analyzes bus data, calculates specified power values, and returns relevant data'''
    
    # Analyze bus data to get relevant information
    bus_type, number_of_buses,p_gen, p_load, q_gen, q_load = analyze_bus_data(file, bus_data,line_data)
    
    # Read bus and line data from Excel sheets
    df_bd = pd.read_excel(file, bus_data) #Bus data from Excel-sheet
    df_ld = pd.read_excel(file,line_data)
    
    # Extract voltage and delta from bus data
    voltage = df_bd["Voltage [pu]"].tolist()
    delta = df_bd["delta[rad]"].tolist() 
    
    # Initialize lists for P_spec and Q_spec values
    P_spec = []
    Q_spec = []


    # Loop through each bus to check bus type, update specified power values and get voltage and delta from bus type
    for i in range(number_of_buses):
        if bus_type[i] == 2:
            Q_spec.append(q_gen[i] - q_load[i])
            voltage[i] = sp.symbols("V_" + str(i + 1))
        else:
            Q_spec.append('Unknown')

        if bus_type[i] != 0:
            P_spec.append(p_gen[i] - p_load[i])
            delta[i] = sp.symbols("delta_" + str(i + 1))
        else:
            P_spec.append('Unknown')
    
    # Returns a list of sympy symbols for all unknown voltages and angles combined.
    vd_unknown = [symbol for symbol in delta if isinstance(symbol, sp.Symbol)]
    vd_unknown += [symbol for symbol in voltage if isinstance(symbol, sp.Symbol)]



    return voltage,delta,P_spec,Q_spec, vd_unknown

def PQ_calculation(Ybus, voltage, delta, P_spec, Q_spec):
    '''
    Generates lists of mathematical expressions for active Power (P) and reactive power (Q)
    at bus locations where they are known. "Unknown" is being used to signify P and Q
    at buses where their values are yet to be determined.
    
    '''
    n = len(P_spec)
    P_calc = []
    Q_calc = []
    PQ_calc = []  # Initialize a list to store the combined P and Q calculations
    PQ_spec = []  # Initialize a list to store the specified and known P and Q values

    # Separate Unknown P and Q values from P_spec and Q_spec
    known_P_values = [P_spec[i] for i in range(n) if P_spec[i] != "Unknown"]
    known_Q_values = [Q_spec[i] for i in range(n) if Q_spec[i] != "Unknown"]

    for i in range(n):
        P_i = 0
        Q_i = 0

        # P_calc equation if it is known
        if P_spec[i] != "Unknown":
            for j in range(n):
                P_i += voltage[i] * voltage[j] * (sp.re(Ybus[i][j]) * sp.cos(delta[i] - delta[j]) + sp.im(Ybus[i][j]) * sp.sin(delta[i] - delta[j]))
            P_calc.append(P_i)
            PQ_calc.append(P_i)  # Add P to PQ_calc
        # Q_calc equation if it is known
        if Q_spec[i] != "Unknown":
            for j in range(n):
                Q_i += voltage[i] * voltage[j] * (sp.re(Ybus[i][j]) * sp.sin(delta[i] - delta[j]) - sp.im(Ybus[i][j]) * sp.cos(delta[i] - delta[j]))
            Q_calc.append(Q_i)
            PQ_calc.append(Q_i)

    PQ_spec = known_P_values + known_Q_values     
    
    # PQ_spec combines the specified values and given values for P and Q, only knowns

    # PQ_known_calc_exp combines the calculated expressions for P and Q calculated, and excludes Unknowns. 
    PQ_known_calc_exp = []
    
    for i in range(len(P_calc)):
        if P_calc[i] != "Unknown":
            PQ_known_calc_exp.append(P_calc[i])
            
    for i in range(len(Q_calc)):
        if Q_calc[i] != "Unknown":
            PQ_known_calc_exp.append(Q_calc[i])

   
    
    return P_calc, Q_calc, np.array(PQ_calc), np.array(PQ_known_calc_exp),PQ_spec

def symbolic_jacobian(PQ_known_calc_exp, vd_unknown):
    '''Calculates the symbolic Jacobian matrix for a given set of symbolic equations for P and Q and unknown voltage and delta.'''
    
    Matrix_PQ_known_calc = sp.Matrix(PQ_known_calc_exp)
    matrix_vd_unknown = sp.Matrix(vd_unknown)
    # Returns symbolic Jacobian by taking symbolic equations and partial derivate with respect to unknown symbols for voltage and deltas
    J = Matrix_PQ_known_calc.jacobian(matrix_vd_unknown)
    return J

def flat_start(vd_unknown):
    '''
    Initializes a flat start (voltages set to 1, and angles 0) for the unknown voltages and angles
    
    '''
    flat_start = [0.0 if symbol.name[0] == 'd' else 1.0 for symbol in vd_unknown]
    return np.array(flat_start)

def convert_symbolic_to_numeric(expression, variables_to_replace, numerical_values):
    """Replace symbolic variables with their numerical values and convert the expression to a numeric value."""
    numeric_expression = expression
    for variable, numerical_value in zip(variables_to_replace, numerical_values):
        numeric_expression = numeric_expression.subs(variable, numerical_value)
    return float(numeric_expression)  # Converts to float for Jacobian inversion compatibility

def convert_symbolic_list_to_numeric(expression_list, variables_to_replace, numerical_values):
    """Replace symbolic variables with their numerical values in a list of expressions and convert them to numeric values."""
    numeric_list = []
    for expression in expression_list:
        numeric_list.append(convert_symbolic_to_numeric(expression, variables_to_replace, numerical_values))
    return np.array(numeric_list)
def numerical_jacobian(J, voltage_unk, voltage_estimate):
    """Converts the Jacobian matrix to its numerical form by substituting symbolic variables."""
    # Get the shape of the Jacobian matrix
    shape_of_J = J.shape
    rows, columns = shape_of_J

    # Use vectorize to apply the conversion function to each element of the Jacobian matrix
    J_filled = np.vectorize(lambda expr: convert_symbolic_to_numeric(expr, voltage_unk, voltage_estimate))(J)

    return J_filled

def numeric_PQ_values(Ybus,voltage,delta,P_spec,Q_spec):
    '''Calculates numerical values for the P_calc and Q_calc equations,
       based on Ybus, voltages, angles and which P and Q are known'''
    

    P_values_updated = np.array([])
    Q_values_updated = np.array([])

    # Calculate numerical values for active power (P)
    n = len(P_spec)
    for i in range(n):
        P_i = 0
        if P_spec[i] == "Unknown":
                for j in range(n):
                    P_i += voltage[i] * voltage[j] * (sp.re(Ybus[i][j]) * sp.cos(delta[i] - delta[j]) + sp.im(Ybus[i][j]) * sp.sin(delta[i] - delta[j]))
                P_values_updated= np.append(P_values_updated,float(P_i))
        else:
            P_values_updated = np.append(P_values_updated, float(P_spec[i]))
    # Calculate numerical values for reactive power (Q)
    for i in range(n):
        Q_i = 0
        if Q_spec[i] == "Unknown":
                for j in range(n):
                    Q_i += voltage[i] * voltage[j] * (sp.re(Ybus[i][j]) * sp.sin(delta[i] - delta[j]) - sp.im(Ybus[i][j]) * sp.cos(delta[i] - delta[j]))
                Q_values_updated = np.append(Q_values_updated,float(Q_i))

        else:
            Q_values_updated = np.append(Q_values_updated,float(Q_spec[i]))
    return P_values_updated,Q_values_updated

def polar_to_rectangular(voltage,delta):
    # Convert polar coordinates to rectangular form using Euler's formula
    # voltage is an array of amplitudes, and delta is an array of angles in radians

    voltage_rect = voltage * np.exp(1j * delta)
    return voltage_rect

def line_flows(file,bus_data, line_data,updated_voltage,updated_delta):
    '''
    Calculates line flows and losses based on updated voltage and delta values.

    '''
    # Read bus and line data from Excel sheets
    df_ld = pd.read_excel(file, line_data)
    df_bd = pd.read_excel(file, bus_data)  # Assuming bus_data is defined somewhere

    # Extract line data from Excel sheet
    f_bus = df_ld["From bus"].dropna().tolist() 
    t_bus = df_ld["To bus"].dropna().tolist() 
    R_line = df_ld["R[pu]"].tolist()
    X_line = df_ld["X[pu]"].tolist()
    adm_line = df_ld["Line Charging Admittance [pu]"]
    

    
    Transformer = df_ld["Transformer(1 = yes, 0 = no)"].tolist()
    S_base_global = df_ld["S_base_global [MVA]"].tolist()
    S_base_global = float(S_base_global[0])

    V_base_global = df_ld["V_base_global [kV]"].tolist()
    V_base_global = float(V_base_global[0])

    

   
    # Update transformer reactance
    for i in range(len(f_bus)):
        if Transformer[i] == 1:
            X_trafo_pu_local = df_ld["X [pu local]"].tolist()
            Trafo_p_voltage = df_ld["Primary voltage"].tolist()
            Trafo_s_voltage = df_ld["Secondary voltage"].tolist()
            Trafo_rating = df_ld["MVA Rating"].tolist()


            result = X_trafo_pu_local[i] * (Trafo_p_voltage[i] / V_base_global)**2 * (S_base_global / Trafo_rating[i])
            rounded_result = round(result, 6)
            X_line[i] = rounded_result

     # Converts list to list with integers
    f_bus = [int(x) for x in f_bus]
    t_bus = [int(x) for x in t_bus]

    # Finding the number of buses
    n_bus = df_bd["Bus"].tolist()
    n = len(f_bus)

    # Calculate line impedance
    Z_line = np.empty(n, dtype = complex)
    for i in range(n):
        Z_line[i] = complex(R_line[i], X_line[i])


    # Convert polar coordinates to rectangular form
    voltage_rect = polar_to_rectangular(updated_voltage,updated_delta)

    P_ij = np.array([])
    P_ji = np.array([])
    Q_ij = np.array([])
    Q_ji = np.array([])
    P_loss = np.array([])
    Q_loss = np.array([])

    S_ij = np.array([])
    S_ji = np.array([])


    # Calculate line flows and losses
    
    for i in range(n):
        I_ij = (voltage_rect[f_bus[i] - 1] - voltage_rect[t_bus[i] - 1]) / Z_line[i]
        P_ij = np.append(P_ij, (voltage_rect[f_bus[i] - 1] * np.conj(I_ij)).real)
        P_ji = np.append(P_ji, (voltage_rect[t_bus[i] - 1] * np.conj(I_ij)).real)
        Q_ij = np.append(Q_ij, (voltage_rect[f_bus[i] - 1] * np.conj(I_ij)).imag)
        Q_ji = np.append(Q_ji, (voltage_rect[t_bus[i] - 1] * np.conj(I_ij)).imag)

        P_loss = np.append(P_loss, P_ij[i] - P_ji[i])
        Q_loss = np.append(Q_loss, Q_ij[i] - Q_ji[i])

        S_ij = np.append(S_ij, np.abs(voltage_rect[f_bus[i] - 1]) * np.abs(np.conj(I_ij)))
        S_ji = np.append(S_ji, np.abs(voltage_rect[t_bus[i] - 1]) * np.abs(np.conj(I_ij)))

    return I_ij,P_ij,P_ji, Q_ij, Q_ji, P_loss, Q_loss, S_ij, S_ji

def print_dataframe_as_latex(dataframe, includeTotalRow=False):
    """
    Generates LaTeX code to present tabular data using the tabularx package.
    
    If includeTotalRow is set to True, the last row in the table is treated as a summary total row.
    """
    
    # Generate LaTeX code
    column_format_tabular = "c" * len(dataframe.columns)  # Create a "c" for each column using the tabular package
    latex_code = dataframe.to_latex(index=False,
                                    header=False,
                                    column_format=column_format_tabular,
                                    position="H",
                                    label="tab:change-me",
                                    caption="\color{red}Change me...")
    
    # Modify the LaTeX package
    latex_code = latex_code.replace("tabular", "tabularx") # Uses the tabularx package instead of the tabular one
    latex_code = latex_code.replace("\\caption", "\\centering\n\setstretch{1.0}\n\caption") # Adds the centering parameter to center the table

    # Modify the column format
    column_format_tabularx = "\n>{\\centering\\footnotesize\\arraybackslash}c" # First column uses the "c" parameter
    for i in range(1, len(dataframe.columns)):
        column_format_tabularx += "\n>{\\centering\\footnotesize\\arraybackslash}X" # All the other columns uses the "X" parameter
    latex_code = latex_code.replace(column_format_tabular, "0.9\\textwidth}{" + column_format_tabularx) # Replaces the column format

    # Use bold text for the header
    header_row = " & ".join(["\\footnotesize{\\textbf{" + col + "}}" for col in dataframe.columns])
    latex_code = latex_code.replace("\\toprule", "\\toprule\n" + header_row + " \\\\")

    # Add total row
    lines = latex_code.split('\n') # Split LaTeX code into lines
    if includeTotalRow:
        total_row_index = 0
        for i, line in enumerate(lines):
            if line.strip() == "\\bottomrule":
                lines.insert(i-1, '\\midrule') # Insert midrule before total row
                total_row_index = i
                break
        # Use bold text for the total row
        total_row = dataframe.iloc[-1].values # Get values in the total row
        total_row_formatted = " & ".join(["\\textbf{" + str(val) + "}" for val in total_row])
        lines[total_row_index] = total_row_formatted + " \\\\"

    # Define indentation rules
    indented_lines = [lines[0]]
    for line in lines[1:-1]:
        if any(line.strip().startswith(prefix) for prefix in ["\\begin{table}", "\\end{table}"]):
            indented_lines.append(line) # No indentation
        elif any(line.strip().startswith(prefix) for prefix in ["\\begin{tabularx}", "\\end{tabularx}", "\\setstretch", "\\centering", "\\caption", "\\label"]):
            indented_lines.append("    " + line) # Single indentation
        else:
            indented_lines.append("        " + line) # Double indentation
    indented_lines.append(lines[-1])
    latex_code = '\n'.join(indented_lines) # Rejoin the lines
    
    print("\n---------------- LaTeX Code -----------------")
    print(latex_code)


def print_bus_data_pu(voltage, delta, P, Q, printLaTeX=False, decimals=3):
    "Prints bus data (in per-unit) to the terminal in a pretty format. If `printLaTeX` = `True`, the function also prints corresponding LaTeX code. Lastly the function returns the dataframe."

    # Store the bus data in a data frame
    bus_data_df = pd.DataFrame()
    bus_data_df["Bus No."] = [str(i) for i in range(1, len(voltage) + 1)]
    bus_data_df["voltage [pu]"] = voltage
    bus_data_df[u"\u03B4" + " [deg]"] = np.degrees(delta)
    bus_data_df["P [pu]"] = P
    bus_data_df["Q [pu]"] = Q

    # Add total row
    bus_data_df = bus_data_df._append({
        "Bus No.": "Totals:",
        "voltage [pu]": "-",
        u"\u03B4" + " [deg]": "-",
        "P [pu]": sum(P),
        "Q [pu]": sum(Q)
    }, ignore_index=True)

    # Round the numeric values to the preferred number of decimals
    bus_data_df = bus_data_df.map(lambda x: f"{x:.{decimals}f}" if isinstance(x, (int, float)) else x)

    print("\n-------------- Bus Data (PU) ----------------")
    print(bus_data_df.to_string(index=False)) # Print bus data without the index column

    if printLaTeX:
        print_dataframe_as_latex(bus_data_df, True) # Print LaTeX code
    
    return bus_data_df


def print_bus_data(voltage, delta, P, Q, S_base, V_base, printLaTeX=False, decimals=3):
    "Prints bus data (in SI) to the terminal in a pretty format. If `printLaTeX` = `True`, the function also prints corresponding LaTeX code. Lastly the function returns the dataframe."

    # Store the bus data in a data frame
    bus_data_df = pd.DataFrame()
    bus_data_df["Bus No."] = [str(i) for i in range(1, len(voltage) + 1)]
    bus_data_df["voltage [kV]"] = voltage*V_base
    bus_data_df[u"\u03B4" + " [deg]"] = np.degrees(delta)
    bus_data_df["P [MW]"] = P*S_base
    bus_data_df["Q [MVAr]"] = Q*S_base

    # Add total row
    bus_data_df = bus_data_df._append({
        "Bus No.": "Totals:",
        "voltage [kV]": "-",
        u"\u03B4" + " [deg]": "-",
        "P [MW]": sum(P)*S_base,
        "Q [MVAr]": sum(Q)*S_base
    }, ignore_index=True)

    # Round the numeric values to the preferred number of decimals
    bus_data_df = bus_data_df.map(lambda x: f"{x:.{decimals}f}" if isinstance(x, (int, float)) else x)

    print("\n-------------- Bus Data (SI) ----------------")
    print(bus_data_df.to_string(index=False)) # Print bus data without the index column

    if printLaTeX:
        print_dataframe_as_latex(bus_data_df, True) # Print LaTeX code
    
    return bus_data_df


def print_branch_data_pu(fromBus, toBus, P_ij, Q_ij, P_ji, Q_ji, P_loss, Q_loss, printLaTeX=False, decimals=3):
    "Prints branch data (in per-unit) to the terminal in a pretty format. If `printLaTeX` = `True`, the function also prints corresponding LaTeX code. Lastly the function returns the dataframe."

    # Store the branch data in a data frame
    branch_data_df = pd.DataFrame()
    branch_data_df["Branch No."] = [str(i) for i in range(1, len(fromBus) + 1)]
    branch_data_df["From Bus"] = [str(int(bus)) for bus in fromBus]
    branch_data_df["To Bus"] = [str(int(bus)) for bus in toBus]
    branch_data_df["P-ij [pu]"] = P_ij
    branch_data_df["Q-ij [pu]"] = Q_ij
    branch_data_df["P-ji [pu]"] = P_ji
    branch_data_df["Q-ji [pu]"] = Q_ji
    branch_data_df["P-loss [pu]"] = P_loss
    branch_data_df["Q-loss [pu]"] = Q_loss

    # Add total row
    branch_data_df = branch_data_df._append({
        "Branch No.": "Totals:",
        "From Bus": "-",
        "To Bus": "-",
        "P-ij [pu]": "-",
        "Q-ij [pu]": "-",
        "P-ji [pu]": "-",
        "Q-ji [pu]": "-",
        "P-loss [pu]": sum(P_loss),
        "Q-loss [pu]": sum(Q_loss)
    }, ignore_index=True)

    # Round the numeric values to the preferred number of decimals
    branch_data_df = branch_data_df.map(lambda x: f"{x:.{decimals}f}" if isinstance(x, (int, float)) else x)

    print("\n------------ Branch Data (PU) ---------------")
    print(branch_data_df.to_string(index=False)) # Print bus data without the index column

    if printLaTeX:
        print_dataframe_as_latex(branch_data_df, True) # Print LaTeX code
    
    return branch_data_df


def print_branch_data(fromBus, toBus, P_ij, Q_ij, P_ji, Q_ji, P_loss, Q_loss, S_base, printLaTeX=False, decimals=3):
    "Prints branch data (in SI) to the terminal in a pretty format. If `printLaTeX` = `True`, the function also prints corresponding LaTeX code. Lastly the function returns the dataframe."

    # Store the branch data in a data frame
    branch_data_df = pd.DataFrame()
    branch_data_df["Branch No."] = [str(i) for i in range(1, len(fromBus) + 1)]
    branch_data_df["From Bus"] = [str(int(bus)) for bus in fromBus]
    branch_data_df["To Bus"] = [str(int(bus)) for bus in toBus]
    branch_data_df["P-ij [MW]"] = P_ij*S_base
    branch_data_df["Q-ij [MVAr]"] = Q_ij*S_base
    branch_data_df["P-ji [MW]"] = P_ji*S_base
    branch_data_df["Q-ji [MVAr]"] = Q_ji*S_base
    branch_data_df["P-loss [MW]"] = P_loss*S_base
    branch_data_df["Q-loss [MVAr]"] = Q_loss*S_base

    # Add total row
    branch_data_df = branch_data_df._append({
        "Branch No.": "Totals:",
        "From Bus": "-",
        "To Bus": "-",
        "P-ij [MW]": "-",
        "Q-ij [MVAr]": "-",
        "P-ji [MW]": "-",
        "Q-ji [MVAr]": "-",
        "P-loss [MW]": sum(P_loss)*S_base,
        "Q-loss [MVAr]": sum(Q_loss)*S_base
    }, ignore_index=True)

    # Round the numeric values to the preferred number of decimals
    branch_data_df = branch_data_df.map(lambda x: f"{x:.{decimals}f}" if isinstance(x, (int, float)) else x)

    print("\n------------ Branch Data (SI) ---------------")
    print(branch_data_df.to_string(index=False)) # Print bus data without the index column

    if printLaTeX:
        print_dataframe_as_latex(branch_data_df, True) # Print LaTeX code
    
    return branch_data_df





        



        
    
        



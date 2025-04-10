import numpy as np

def Climate_Sys(year, x, u, c):
    """
    The ODE's of the EnRoads climate module, converted from MATLAB.
    
    Parameters:
        year : float
            The current year.
        x : numpy.ndarray
            The state vector. x[0:5] represents Q (heat in surface and ocean layers),
            x[5] is C (atmospheric carbon), x[6] is N (N2O), and x[7] is M (CH4).
        u : any
            The control vector (unused in this function).
        c : object
            A constants object (or dict) with required parameters as attributes.
    
    Returns:
        dxdt : numpy.ndarray
            The derivative of the state vector.
        y : numpy.ndarray
            The output vector (temperatures in atmosphere and ocean).
    """
    # Extract states from x (MATLAB indices 1:5 become Python indices 0:5)
    Q = x[:5]        # Heat in surface and ocean layers [Watt year / sq.m]
    C = x[5]         # Carbon in atmosphere [Gton]
    N = x[6]         # N2O in atmosphere [Mton]
    M = x[7]         # CH4 in atmosphere [Mton]
    # Note: Other gases (PFCs, SF6, HFCs, ODS) are not used in this code.

    # Algebraic equations -----------------------------------------------------
    # Calculate parts of the total GHG radiative forcing using the provided constants.
    co2_ppm = C * c.ppm_co2_per_GtonC / c.unit_ppm
    ch4_ppb = M * c.ppb_ch4_per_mton_ch4 / c.unit_ppb
    n2o_ppb = N * c.ppb_n2o_per_mton_n2o / c.unit_ppb
    
    # Determine alpha_co2 based on co2_ppm
    if co2_ppm < c.co2_ref:
        alpha_co2 = c.d1
    elif co2_ppm > c.co2_alpha_max:
        alpha_co2 = c.d1 - c.b1**2 / (4 * c.a1)
    else:
        alpha_co2 = c.d1 + c.a1 * (co2_ppm - c.co2_ref)**2 + c.b1 * (co2_ppm - c.co2_ref)
    
    alpha_n2o = c.c1 * np.sqrt(n2o_ppb)
    
    rf_co2 = (alpha_co2 + alpha_n2o) * np.log(C / c.c_preindustrial)
    
    rf_ch4 = (c.a3 * np.sqrt(ch4_ppb) + c.b3 * np.sqrt(n2o_ppb) + c.d3) * \
             (np.sqrt(ch4_ppb) - np.sqrt(c.ch4_ref))
    
    rf_n2o = (c.a2 * np.sqrt(co2_ppm) + c.b2 * np.sqrt(n2o_ppb) + c.c2 * np.sqrt(ch4_ppb) + c.d2) * \
             (np.sqrt(n2o_ppb) - np.sqrt(c.n2o_ref))
    
    # For now, these are set to zero:
    rf_f_gasses = 0
    rf_ods      = 0
    
    well_mixed_ghg_forcing = rf_co2 + rf_ch4 + rf_n2o + rf_f_gasses + rf_ods
    
    other_forcing = -1.0 + (year - 1990) * ((-0.3) - (-1)) / (2100 - 1990)
    
    hydrogen_released_per_year_in_Tg = 0.0  # Constant, for now.
    h2_effect_on_rf = hydrogen_released_per_year_in_Tg * (
        c.ch4_erf_per_h2_flux + c.o3_erf_per_h2_flux + 
        c.strat_h2o_erf_per_h2_flux + c.aerosol_rf_per_h2_flux)
    
    total_radiative_forcing = well_mixed_ghg_forcing + other_forcing + h2_effect_on_rf
    
    # Temperature at all layers -----------------------------------------------
    # Element-wise division (assuming c.heat_capacity is a NumPy array of appropriate shape)
    temperature = Q / c.heat_capacity
    
    # Feedback cooling --------------------------------------------------------
    x2_co2_forcing = (alpha_co2 + alpha_n2o) * np.log(2)
    climate_feedback_param = x2_co2_forcing / c.climate_sensitivity_to_2x_co2
    feedback_cooling = temperature[0] * climate_feedback_param
    
    # Heat transfer between ocean layers:
    # Assuming c.heat_transfer_coeff is an array with 4 elements and 
    # c.layer_depth is an array with 5 elements.
    heat_transfer_ocean_layers = (temperature[:4] - temperature[1:5]) * c.heat_transfer_coeff / \
                                 ((c.layer_depth[:4] + c.layer_depth[1:5]) / 2)
    
    # Differential Equations --------------------------------------------------
    # dQdt for the surface and ocean layers
    dQdt_surface = total_radiative_forcing - feedback_cooling - heat_transfer_ocean_layers[0]
    
    # For ocean layers (indices 1 to 4 in MATLAB correspond to 0 to 3 in Python)
    # The MATLAB code subtracts the next element, with the last element subtracting 0.
    dQdt_ocean = np.empty(4)
    dQdt_ocean[:-1] = heat_transfer_ocean_layers[:-1] - heat_transfer_ocean_layers[1:]
    dQdt_ocean[-1] = heat_transfer_ocean_layers[-1]
    
    # Derivatives for gases (noise terms are effectively zero)
    dCdt = 0.005 * C  # + 0.000 * np.random.randn() if noise were needed
    dNdt = 0.004 * N  # + 0.000 * np.random.randn()
    dMdt = 0.003 * M  # + 0.000 * np.random.randn()
    
    # Combine derivatives into one state derivative vector
    # dQdt has 5 elements: surface (1) and 4 ocean layers
    dxdt = np.concatenate((
        np.array([dQdt_surface]),
        dQdt_ocean,
        np.array([dCdt, dNdt, dMdt])
    ))
    
    # Outputs: temperatures in atmosphere and ocean
    y = temperature
    
    return dxdt, y
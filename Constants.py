import numpy as np
class EnRoadsConstants:
    # Universal constants
    c100_percent = 100             # [percent]
    one_year = 1.0                 # [year]
    unit_ppb = 1.0                 # [ppb]
    unit_ppm = 1.0                 # [ppm]
    hours_per_year = 24 * 365.25   # [hours/year]
    seconds_per_year = hours_per_year * 60 * 60  # [sec/year]
    ppt_per_ppb = 1000             # [ppt/ppb]
    ppt_per_mol = 5.68e-09         # [ppt/mole]
    g_per_ton = 1e6                # [g/ton]
    kg_per_tg = 1e9                # [kg/tg]
    ton_per_mton = 1e6             # [ton/Mton]
    mton_c_per_gton_c = 1000       # [MtonsC/GtonsC]
    watt_s_per_J = 1.0             # [Watt-sec/Joule]
    MJ_per_EJ = 1e12               # [MJ/EJ]
    GJ_per_kWh = 0.0036            # [GJ/kWh]
    density = 1000                 # [kg/m^3] density of water
    mass_heat_cap = 4186           # [J/kg/°C] heat capacity of water

    present_year = 2025            # used to shift from data-driven to policy-driven structure

    # Earth model constants
    earth_area = 5.1e14            # [m^2]
    land_area_fraction = 0.292     # fraction of earth that is land
    land_thickness = 8.4           # [m]
    layer_depth = np.array([100, 300, 300, 1300, 1800])  # [m] (MATLAB column vector converted to Python list)
    initial_temp_change_from_1750 = np.array([0.7, 0.34, 0.17, 0.0156, 5.9e-04])  # [°C]
    mean_temperature_change_from_v_1750_to_v_1850 = 0.05  # [°C]
    land_layers = np.r_[1, np.zeros(len(layer_depth)-1)]

    # Compute layer_volume:
    # First element: effective surface layer (land + ocean)
    # Next elements: ocean layers only.
    layer_volume = earth_area* ((land_area_fraction*land_thickness)*land_layers # Land area depth
                                +(1-land_area_fraction)*layer_depth)
 
    # Volumetric heat capacity [Watt-year/(°C·m^3)]
    volumetric_heat_capacity = mass_heat_cap * watt_s_per_J / seconds_per_year * density

    # Heat capacity of each layer (per unit area) [Watt-year/(°C·m^2)]
    heat_capacity = layer_volume * volumetric_heat_capacity / earth_area

    # Initial heat in each layer (per unit area) [Watt-year/m^2]
    heat_initial = np.array([dt * hc for dt, hc in zip(initial_temp_change_from_1750, heat_capacity)])

    heat_transfer_rate = 1.23      # [W/(m^2·°C)]
    heat_diffusion_covar = 1       # (dimensionless)
    eddy_diff_mean = 4400          # [m^2/year] mean eddy diffusion coefficient
    ocean_heat_mixing_time = 1.0   # [year]
    mixing_time = 1.0              # [year]

    eddy_diff_coeff_index = 1      # (dimensionless)
    eddy_diff_coeff = eddy_diff_coeff_index * eddy_diff_mean  # [m^2/year]

    # Compute heat_transfer_coeff for each interface between layers.
    # In MATLAB: (layer_depth(1:4) + layer_depth(2:5))/2 gives the mean depth between layers.
    
    heat_transfer_coeff = heat_transfer_rate*(layer_depth[:-1] + layer_depth[1:])/2*heat_diffusion_covar*(eddy_diff_coeff/eddy_diff_mean+(1-heat_diffusion_covar))

    # Climate model constant
    climate_sensitivity_to_2x_co2 = 3  # [°C]

    # GHG constants
    co2_per_c = 3.66         # [GtonsCO2/GtonsC]
    ch4_per_c = 1.33         # [GtonsCO2/GtonsC]

    # Radiant Heat Forcing Coefficients for CO2 (IPCC AR6)
    a1 = -2.4785e-7          # [W/m^2/ppm^2]
    b1 = 7.5906e-4           # [W/m^2/ppm]
    c1 = -2.1492e-3          # [W/m^2/ppb^(0.5)]
    d1 = 5.2488              # [W/m^2]
    C0 = 353.98              # [ppm] initial carbon concentration (1990 Mauna Loa, 2020)
    co2_ref = 277.15         # [ppm] CO2 reference concentration
    ppm_co2_per_GtonC = 0.4695  # [ppm/GtonC]
    c_preindustrial = 590    # [ppm] preindustrial Carbon concentration (not CO2)
    co2_alpha_max = C0 - b1 / (2 * a1)  # CO2 threshold in EnRoads

    # N2O coefficients
    a2 = -3.4197e-4          # [W/m^2/ppm]
    b2 = 2.5455e-4           # [W/m^2/ppb]
    c2 = -2.4357e-4          # [W/m^2/ppb]
    d2 = 0.12173             # [W/m^2/ppb^(0.5)]
    N0 = 308.725             # [ppb] initial N2O concentration (1990 levels)
    n2o_ref = 273.87         # [ppb] N2O reference concentration
    natural_n2o_emissions = 11.2  # [Mton/year]
    n2o_N_molar_mass = 28    # [g/mole]
    ppb_n2o_per_mton_n2o = (ppt_per_mol / n2o_N_molar_mass *
                             g_per_ton * ton_per_mton) / ppt_per_ppb
    time_const_for_n2o = 121   # [years]
    initial_n2o_in_atm = N0 / ppb_n2o_per_mton_n2o  # [Mtons]

    # CH4 coefficients
    a3 = -8.9603e-5          # [W/m^2/ppb]
    b3 = -1.2462e-4          # [W/m^2/ppb]
    d3 = 0.045194            # [W/m^2/ppb^(0.5)]
    M0 = 1714.0              # [ppm] initial CH4 concentration (1990 levels)
    ch4_ref = 731.41         # [ppb] CH4 reference concentration
    ch4_molar_mass = 16      # [g/mole]
    ppb_ch4_per_mton_ch4 = (ppt_per_mol / ch4_molar_mass *
                            g_per_ton * ton_per_mton) / ppt_per_ppb
    time_const_for_ch4 = 10.3  # [years]
    initial_ch4_in_atm = M0 / ppb_ch4_per_mton_ch4  # [Mtons]

    # ODS (Ozone Depleting Substances)
    ods_molar_mass = 120     # [g/mole]

    # PFC (Perfluorocarbons)
    preindustrial_pfc_conc = 40         # [ppt]
    pfc_radiative_efficiency = 0.099    # [W/(ppb·m^2)]

    # SF6 (Sulfur Hexafluoride)
    sf6_molar_mass = 146     # [g/mole]
    preindustrial_sf6_conc = 0           # [ppt]
    sf6_radiative_efficiency = 0.567     # [W/(ppb·m^2)]

    # HFC (Hydrofluorocarbons)
    hfc_molar_mass = [102, 70, 52, 120, 84, 66, 170, 134, 252]  # [g/mole]
    preindustrial_hfc_conc = 0           # [ppt]
    hfc_radiative_efficiency = [0.167, 0.191, 0.111, 0.234,
                                0.168, 0.102, 0.237, 0.24, 0.357]  # [W/(ppb·m^2)]

    # H2 (Hydrogen) related constants
    h2_storage_leakage_rate = 2          # [percent/year]
    average_utilization_of_storage_capacity = 0.5  # (dimensionless)
    energy_intensity_of_hydrogen = 120   # [MJ/kg H2]
    ch4_erf_per_h2_flux = 0.00046        # [W/(m^2)/(tg H2/year)]
    o3_erf_per_h2_flux = 0.00040         # [W/(m^2)/(tg H2/year)]
    strat_h2o_erf_per_h2_flux = 0.00019   # [W/(m^2)/(tg H2/year)]
    aerosol_rf_per_h2_flux = 0           # [W/(m^2)/(tg H2/year)]

    # Carbon concentrations
    C_init = np.array([C0, 1044.4, 3114.46, 3099.9, 13361.1, 18483])  # initial C in each layer
    preind_ocean_c_per_meter = 10.2373   # [GtonsC/meter]
    preind_c_in_mixed_layer = preind_ocean_c_per_meter * layer_depth[0]  # [GtonsC]

    sensitivity_of_c_uptake_to_temperature = 1  # (dimensionless)
    sensitivity_of_pco2_dic_to_temperature_mean = 0.3  # [percent/°C]
    sensitivity_of_pco2_dic_to_temperature = (
        sensitivity_of_c_uptake_to_temperature * sensitivity_of_pco2_dic_to_temperature_mean / c100_percent
    )
    flux_ocean_c_to_ch4_as_ch4 = 30  # [Mtons/Year] (as CH4)
    flux_ocean_c_to_ch4 = flux_ocean_c_to_ch4_as_ch4 / ch4_per_c / mton_c_per_gton_c
    buff_c_coeff = 3.92
    ref_buffer_factor = 9.7
    buffer_factor_init = ref_buffer_factor

    baseline_change= np.array([-0.015,-0.005,-0.01,-0.007,-0,-0.007]) #Every baseline chage except CO2
    Tj_per_Ej = 1e6
    bioenergy_C_of_forests = 0.027

    nonenergy_annual_improvement = -0.015
    lifetime_products = 300
    
    tons_per_Mton = 1e6
    historic_change_in_F_gas_for_solvents = -0.05
    expected_change_in_F_gas_for_solvents = -0.08
    N2O_per_N = 1.57143
    initial_fraction_of_GHG_emissions_from_demand_capital = np.array([0, 0.0005, 0.092, 1, 1, 1])
    initial_fraction_of_demand_emissions_from_solvents = np.array([0, 0, 0, 0.7, 0.5, 0.05])
    base_F_gas_leak_rate = 2
    base_F_gas_recycling_rate = 25
    solvent_use_lifetime = 2
    baseline_change_in_HFCs_for_cooling = -0.011
    initial_cumulative_emissions = 1225
    historic_ODS_stockpile = 1e6
    lower_limit_to_waste_production_ratio_and_intensity = 0.4
    years_to_achieve_mature_forest_degradation_policy = 10
    switch_use_land_detailed_settings = 0
    base_ODS_smooth_time = 10
    initial_ODS_stock_GtCO2 = 20
    GWP_of_ODS = 6000
    tonsCO2_per_GtonCO2 = 1e9
    tonsC_per_GtonC = 1e9
    target_reduction_in_mature_forest_degradation = 0
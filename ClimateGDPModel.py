import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class ClimateGDPModel:
    def __init__(self, start_year=2023, end_year=2100, damage_model='burke_sr'):
        # Basic settings
        self.time = np.arange(start_year, end_year + 1)
        self.damage_model = damage_model

        # Constants
        self.last_gdp_year = 2023
        self.present_year = 2025
        self.T2021_PPP = 1e12
        self.social_discount_rate = 4.25
        self.one_year = 1
        self.people_per_thousand = 1000
        self.regional_near_term_gdp_rate = 2.5
        self.long_term_gdp_rate = 1.5
        self.transition_time = 75
        self.factor_exp_change = 3
        self.temperature_change_1750 = 2.0
        self.mean_temp_change_1750_to_1850 = 0.05
        self.unit_degreeC = 1
        self.population = np.linspace(8e9, 10.2e9, len(self.time))
        self.sectors = ['Residential', 'Industry', 'Transport']
        self.initial_sector_share = np.array([0.3, 0.5, 0.2])
        self.normal_capacity_utilization = 0.9
        self.last_historical_gdp_per_capita = 10000

        # Coefficients for damage functions
        self.damage_coeffs = {
            'burke_sr': {'coeff1': 0.226613, 'coeff2': -0.0459135, 'coeff3': 0.0045574, 'exp1': 2, 'exp2': 3, 'base': 0.88},
            'burke_lr': {'coeff1': 0.324059, 'coeff2': 0.0821253, 'coeff3': 0.0116978, 'exp1': 2, 'exp2': 3, 'base': 0.8},
            'dietz_stern': {'denom1': 18.8, 'exp1': 2, 'denom2': 4, 'exp2': 6.754},
            'howard_sterner': {'coeff': 1.145, 'power': 2},
        }

        # Outputs
        self.gdp_rate_ratio = None
        self.gdp_proj = None
        self.population_series = None
        self.gdp_by_sector = None
        self.gdp_growth_rate = None

        # Run model steps
        self._simulate_gdp_rate_ratio()
        self._simulate_gdp_per_capita()
        self._simulate_damage_and_loss()

    def _zidz(self, n, d):
        return n / d if d != 0 else 0

    def _simulate_gdp_rate_ratio(self):
        long_term_rate_ratio = self._zidz(self.long_term_gdp_rate, self.regional_near_term_gdp_rate)
        rate_ratio = [1.0]

        for t in range(1, len(self.time)):
            current = rate_ratio[-1]
            if self.time[t] < self.last_gdp_year:
                delta = 0
            else:
                gap = self._zidz(long_term_rate_ratio - current, current)
                adjustment = gap / (self.transition_time / self.factor_exp_change)
                delta = adjustment * current
            rate_ratio.append(current + delta)

        self.gdp_rate_ratio = np.array(rate_ratio)

    def _simulate_gdp_per_capita(self):
        gdp_proj = pd.DataFrame(index=self.time, columns=self.sectors, dtype=float)
        for sector in self.sectors:
            gdp_proj.at[self.last_gdp_year, sector] = self.last_historical_gdp_per_capita * self.initial_sector_share[self.sectors.index(sector)]

        for t in range(1, len(self.time)):
            year = self.time[t]
            prev = gdp_proj.iloc[t - 1]
            growth = 0 if year < self.last_gdp_year else (self.regional_near_term_gdp_rate * self.gdp_rate_ratio[t]) / 100
            gdp_proj.iloc[t] = prev + prev * growth

        self.gdp_proj = gdp_proj
        self.population_series = self.population
        self.gdp_by_sector = gdp_proj.multiply(self.population[:, np.newaxis]) / self.T2021_PPP
        self.gdp_growth_rate = self.gdp_by_sector.pct_change() * 100

    def _compute_damage(self, delta_T):
        model = self.damage_model
        if model in ['burke_sr', 'burke_lr']:
            p = self.damage_coeffs[model]
            x = (delta_T - p['base']) / self.unit_degreeC
            return 100 * (1 - 1 / (1 + p['coeff1'] * x + p['coeff2'] * x ** p['exp1'] + p['coeff3'] * x ** p['exp2']))
        elif model == 'dietz_stern':
            p = self.damage_coeffs[model]
            x = delta_T / self.unit_degreeC
            return 100 * (1 - 1 / (1 + (x / p['denom1']) ** p['exp1'] + (x / p['denom2']) ** p['exp2']))
        elif model == 'howard_sterner':
            p = self.damage_coeffs[model]
            x = delta_T / self.unit_degreeC
            return p['coeff'] * x ** p['power']
        return 0.0

    def _simulate_damage_and_loss(self):
        delta_temp = self.temperature_change_1750 - self.mean_temp_change_1750_to_1850
        damage = self._compute_damage(delta_temp)
        damage_series = np.full_like(self.time, damage)

        gdp_weighted = self.gdp_proj.mul(self.initial_sector_share * self.normal_capacity_utilization, axis=1).sum(axis=1)
        gdp_per_capita_loss = gdp_weighted * (damage_series / 100)
        gdp_loss = gdp_per_capita_loss.values * self.population / self.T2021_PPP
        discount_factor = 1 / ((1 + self.social_discount_rate / 100) ** np.maximum(self.time - self.present_year, 0))
        self.present_value_gdp_loss = gdp_loss * discount_factor
        self.cumulative_pv_loss = np.cumsum(self.present_value_gdp_loss)

# Run the model
model = ClimateGDPModel()

# Output GDP values and growth rates
import ace_tools as tools; tools.display_dataframe_to_user(
    name="GDP Values and Growth Rates by Sector",
    dataframe=pd.concat([model.gdp_by_sector, model.gdp_growth_rate.add_suffix(" Growth (%)")], axis=1)
)

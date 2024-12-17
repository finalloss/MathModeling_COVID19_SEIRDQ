import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import differential_evolution

def loss_func(params, model, t, observed):
    predicted = model(t, *params)
    return np.sum((predicted - observed) ** 2)


class EpidemicModel:
    def __init__(self, data_file, country, model='SEIR'):
        self.data_file = data_file
        self.country = country
        self.model = model
        self.new_cases_or_deaths = None
        self.fitted_data = None
        self.t = None
        self.popt = None

    @staticmethod
    def seir_model(t, N, beta, gamma, sigma):
        def dydt(y, t, N, beta, gamma, sigma):
            S, E, I, R = y
            dSdt = -beta * S * (I+E) / N
            dEdt = beta * S * (I+E) / N - sigma * E
            dIdt = sigma * E - gamma * I
            dRdt = gamma * I
            return dSdt, dEdt, dIdt, dRdt

        S0, E0, I0, R0 = N - 1, 1, 0, 0
        y0 = S0, E0, I0, R0
        ret = odeint(dydt, y0, t, args=(N, beta, gamma, sigma))
        _, _, I, _ = ret.T
        return np.diff(np.concatenate(([0], I)))

    @staticmethod
    def seirdq_model(t, L, t0, N, beta, gamma, sigma, mu, alpha ):
        def logistic(t, L, t0):
            return L * np.exp((t - t0)) / (1 + np.exp((t - t0))) if t > t0 else 0

        def dydt(y, t, L, t0, N, beta, gamma, sigma, mu, alpha):
            S, E, I, R, D, Q = y
            rho_t = logistic(t, L, t0)
            dSdt = -beta * S * (I+E) / (N-D) -  rho_t * S + alpha * R
            dEdt = beta * S * (I+E) / (N-D) - sigma * E -  rho_t * E
            dIdt = sigma * E - (gamma + mu) * I -  rho_t * I
            dRdt = gamma * (I + Q) - alpha * R
            dDdt = mu * I
            dQdt =  rho_t * (S + E + I)
            return dSdt, dEdt, dIdt, dRdt, dDdt, dQdt

        S0, E0, I0, R0, D0, Q0 = N - 100, 100, 0, 0, 0, 0
        y0 = S0, E0, I0, R0, D0, Q0
        ret = odeint(dydt, y0, t, args=(L, t0, N, beta, gamma, sigma, mu, alpha))
        _, _, I, _, _, _ = ret.T
        return np.diff(np.concatenate(([0], I)))

    def fit(self, start_date, end_date, bounds):
        all_data = pd.read_csv(self.data_file)
        self.start_date = start_date
        self.end_date = end_date
        country_data = all_data[all_data['Country'] == self.country]
        filtered_data = country_data[(country_data['Date_reported'] >= start_date) & (country_data['Date_reported'] <= end_date)]
        self.t = np.linspace(0, len(filtered_data), len(filtered_data))
        
        if self.model == 'SEIR':
            self.new_cases_or_deaths = filtered_data['New_cases'].values
            result = differential_evolution(loss_func, bounds, args=(self.seir_model, self.t, self.new_cases_or_deaths))
            self.popt = result.x
            self.fitted_data = self.seir_model(self.t, *self.popt)
        elif self.model == 'SEIRDQ':
            self.new_cases = filtered_data['New_cases'].values
            result = differential_evolution(loss_func, bounds, args=(self.seirdq_model, self.t, self.new_cases))
            self.popt_cases = result.x
            self.fitted_cases = self.seirdq_model(self.t, *self.popt_cases)
            print("Fitted parameters: ", self.popt_cases)
            print("Fitted cases: ", self.fitted_cases)
        else:
            raise ValueError("Invalid model. Choose 'SEIR' or 'SEIRDQ'.")


    def plot_fitting(self):
        sns.set(style='whitegrid', font_scale=1.2)
        if self.model == 'SEIR':
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(self.t, self.new_cases_or_deaths, 'bo', label='Observed data', markersize=4)
            ax.plot(self.t, self.fitted_data, 'r-', label='Fitted data', linewidth=2)
            ax.set_xlabel('Days since start date', fontsize=14)
            ax.set_ylabel('Number of cases', fontsize=14)
            ax.set_title(f'{self.model} Model Fitting for {self.country} From {self.start_date} to {self.end_date}' , fontsize=16)
            ax.legend(loc='upper right', fontsize=12)
            ax.set_xlim(min(self.t), max(self.t))
            ax.set_ylim(0, max(max(self.new_cases_or_deaths), max(self.fitted_data)) * 1.1)
            ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
            ax.grid(True, linestyle='--', linewidth=0.5)
            N, beta, gamma, sigma = self.popt
            ax.text(max(self.t) * 0.6, max(max(self.new_cases_or_deaths), max(self.fitted_data)) * 0.9, f"N = {N:.1f}")
            ax.text(max(self.t) * 0.6, max(max(self.new_cases_or_deaths), max(self.fitted_data)) * 0.85, f"beta = {beta:.4f}")
            ax.text(max(self.t) * 0.6, max(max(self.new_cases_or_deaths), max(self.fitted_data)) * 0.8, f"gamma = {gamma:.4f}")
            ax.text(max(self.t) * 0.6, max(max(self.new_cases_or_deaths), max(self.fitted_data)) * 0.75, f"sigma = {sigma:.4f}")
            plt.savefig(f"./img/{self.model} Model Fitting for {self.country} From {self.start_date} to {self.end_date}.png")
            plt.show()
        elif self.model == 'SEIRDQ':
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(self.t, self.new_cases, 'bo', label='Observed cases', markersize=4)
            ax.plot(self.t, self.fitted_cases, 'r-', label='Fitted cases', linewidth=2)
            ax.text(max(self.t) * 0.6, max(max(self.new_cases), max(self.fitted_cases)) * 0.9, 
                    f"L: {self.popt_cases[0]:.4f}, t0: {self.popt_cases[1]:.1f}, N: {self.popt_cases[2]:.1f}", fontsize=10)
            ax.text(max(self.t) * 0.6, max(max(self.new_cases), max(self.fitted_cases)) * 0.85, 
                    f"beta: {self.popt_cases[3]:.4f}, gamma: {self.popt_cases[4]:.4f}, sigma: {self.popt_cases[5]:.4f}", fontsize=10)
            ax.text(max(self.t) * 0.6, max(max(self.new_cases), max(self.fitted_cases)) * 0.8, 
                    f"mu: {self.popt_cases[6]:.4f}, alpha: {self.popt_cases[7]:.4f}", fontsize=10)

            ax.set_xlabel('Days since start date', fontsize=14)
            ax.set_ylabel('Number of cases', fontsize=14)
            ax.set_title(f'{self.model} Model Fitting for {self.country} (Cases) From {self.start_date} to {self.end_date}', fontsize=16)
            ax.legend(loc='upper right', fontsize=12)
            ax.set_xlim(min(self.t), max(self.t))
            ax.set_ylim(0, max(max(self.new_cases), max(self.fitted_cases)) * 1.1)
            ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
            ax.grid(True, linestyle='--', linewidth=0.5)
            plt.savefig(f"./img/{self.model} Model Fitting cases for {self.country} (Cases) From {self.start_date} to {self.end_date}.png")
            plt.show()

        else:
            raise ValueError("Invalid model. Choose 'SEIR' or 'SEIRDQ'.")


if __name__ == '__main__':
    data_file = 'data/WHO-COVID-19-global-data.csv'
    country = 'China'
    start_date = '2022-11-25'
    end_date = '2023-01-15'

    """ seir_model = EpidemicModel(data_file, model='SEIR', country = country)
    seir_model.fit(start_date, end_date, [(1000000,1400000000), (0.1,1), (0.001,1),(0.1,1)])
    seir_model.plot_fitting() """
    seirdq_model = EpidemicModel(data_file, model='SEIRDQ', country = country)
    seirdq_model.fit(start_date, end_date,  [(0.0001,0.05), (20,40), (1000000,1400000000), (0.2,1), (0.001,0.4), (0.1,1), (0.005,0.1),(0.00005,0.01)])
    seirdq_model.plot_fitting()

import pandas as pd
import numpy as np
from scipy.integrate import simps
import warnings
import matplotlib.pyplot as plt

def get_massless_data(s, l, m):
    path = "./Z_s" + str(s) + "_l" + str(l) + "_m" + str(m) + ".dat"
    heada0 = ["eta", "omega", "Z"]
    head = ["a", "eta", "omega", "Z"]
    a0_data = pd.read_csv("./Z_a0_s" + str(s) + "_l" + str(l) + ".dat",
                            sep = "\t", names = heada0, dtype = float)
    data = a0_data.append(pd.read_csv(path, sep = "\t", names = head, dtype = float),
                            sort = False, ignore_index = True)
    data = data.fillna(0)[head]
    return data

def get_massive_data():
    head = ["a", "eta", "mu", "omega", "Z"]
    path = "./Z_s0_l1_m1_massive.dat"
    return pd.read_csv(path, sep = "\t", names = head, dtype = float)


# Read and manipulate data
class AmplificationFactors:
    # Called at instantiation
    def __init__(self, s, l, m, massive : bool = False):
        self.massive = massive
        self.s, self.l, self.m = s, l, m
        self.check_slm()
        self.data = self.read_data()

    # Function that checks eventual errors in user's input
    def check_slm(self):
        # Check type of s, l and m
        if type(self.s) is not int:
            raise ValueError("s = " + str(self.s) + " is not of type int")
        if type(self.l) is not int:
            raise ValueError("l = " + str(self.l) + " is not of type int")
        if type(self.m) is not int:
            raise ValueError("m = " + str(self.m) + " is not of type int")
        # If user's selecting massive data, s = 0, l = 1 and m = 1 or raise error
        if self.massive:
            if self.s != 0:
                raise ValueError("s = " + str(self.s) + " is not valid for massive data")
            if self.l != 1:
                raise ValueError("l = " + str(self.l) + " is not valid for massive data")
            if self.m != 1:
                raise ValueError("m = " + str(self.m) + " is not valid for massive data")
        # If user's selecting massless data, s <= 1, s <= l, l <= 2 and abs(m) <= l
        else:
            if self.s > 1 or self.s < 0:
                raise ValueError("There are no data with s = " + str(self.s))
            if self.l > 2:
                raise ValueError("There are no data with l = " + str(self.l))
            if self.l < self.s:
                raise ValueError("l must be greater or equal to s")
            if self.l < np.abs(self.m):
                raise ValueError("The absolute value of m must be smaller or equal to l")

    # Data reading function
    def read_data(self):
        if self.massive:
            return get_massive_data()
        return get_massless_data(self.s, self.l, self.m)

    # Selects data by the spin a of the object and returns pd.DataFrame
    def select_data_spin(self, a : float, data : pd.DataFrame = pd.DataFrame()):
        if data.empty:
            data = self.data
        df = data[data["a"] == a]
        if df.empty:
            warnings.warn("There are no data for a = " + str(a), RuntimeWarning)
        return df

    def select_data_eta(self, eta : float, data : pd.DataFrame = pd.DataFrame()):
        if data.empty:
            data = self.data
        df = data[data["eta"] == eta]
        if df.empty:
            warnings.warn("There are no data for eta = " + str(eta), RuntimeWarning)
        return df

    def select_data_mu(self, mu : float, data : pd.DataFrame = pd.DataFrame()):
        if not self.massive:
            raise ValueError("Loaded data are not massive data")
        if data.empty:
            data = self.data
        df = data[data["mu"] == mu]
        if df.empty:
            warnings.warn("There are no data for mu = " + str(mu), RuntimeWarning)
        return df

    def get_Zmax(self, data : pd.DataFrame):
        return data[data["Z"] == data["Z"].max()]

    def get_Zmax_fixed_spin(self, a : float, data : pd.DataFrame = pd.DataFrame()):
        if data.empty:
            data = self.data
        if not (a in data["a"].unique()):
            warnings.warn("No data for a = " + str(a), RuntimeWarning)
        a_data = data[data["a"] == a]
        Zmax = pd.DataFrame(columns = a_data.columns)
        if not self.massive:
            for eta in a_data["eta"].unique():
                eta_data = a_data[a_data["eta"] == eta]
                Zmax = Zmax.append(eta_data[eta_data["Z"] == eta_data["Z"].max()],
                        ignore_index = True, sort = False)
            return Zmax
        else:
            for mu in a_data["mu"].unique():
                mu_data = a_data[a_data["mu"] == mu]
                for eta in mu_data["eta"].unique():
                    eta_data = mu_data[mu_data["eta"] == eta]
                    Zmax = Zmax.append(eta_data[eta_data["Z"] == eta_data["Z"].max()],
                            ignore_index = True, sort = False)
            return Zmax

    def get_Zmax_fixed_pars(self, a : float, eta : float, mu : float = None,
                            data : pd.DataFrame = pd.DataFrame()):
        if data.empty:
            data = self.data
        if not (a in data["a"].unique()):
            warnings.warn("No data for a = " + str(a), RuntimeWarning)
        if not (eta in data["eta"].unique()):
            warnings.warn("No data for eta = " + str(eta))
        if mu != None and (not (mu in data["mu"].unique())):
            warnings.warn("No data for mu = " + str(mu), RuntimeWarning)
        if data.empty:
            data = self.data
        if not self.massive:
            df = data[(data["a"] == a) & (data["eta"] == eta)]
            return df[df["Z"] == df["Z"].max()]
        else:
            df = data[(data["a"] == a) & (data["eta"] == eta) & (data["mu"] == mu)]
            return df[df["Z"] == df["Z"].max()]

    def plot_spectra(self, ax : plt.Axes, a : float, eta : float, mu : float = None, label : str = None):
        if not self.massive:
            if label == None:
                label = "$\eta/M^3 = %g$" % eta
            condition = (self.data["a"] == a) & (self.data["eta"] == eta)
            df = self.data[condition]
            if df.empty:
                warnings.warn("No data for a = %g and eta = %g" % (a, eta), RuntimeWarning)
            ax.plot(df["omega"], df["Z"], label = label)
            ax.set_xlabel("$\omega M$", fontsize = 15)
            ax.set_ylabel("$Z_{%g, %g, %g}$" % (self.s, self.l, self.m), fontsize = 15)
            return ax
        else:
            condition = (self.data["a"] == a) & (self.data["eta"] == eta) & (self.data["mu"] == mu)
            df = self.data[condition]
            if df.empty:
                warnings.warn("No data for a = %g, eta = %g and mu = %g" % (a, eta, mu), RuntimeWarning)
            ax.plot(df["omega"], df["Z"], label = "$M\mu = %g$" % mu)
            ax.set_xlabel("$M\omega$", fontsize = 15)
            ax.set_ylabel("$Z_{%g, %g, %g}$" % (self.s, self.l, self.m), fontsize = 15)
            return ax

    def plot_Zmax(self, ax : plt.Axes, a : float,  mu : float = None):
        if not self.massive:
            df = self.select_data_spin(a)
            df_Zmax = pd.DataFrame()
            for eta in df["eta"].unique():
                df_Zmax = df_Zmax.append(self.get_Zmax_fixed_pars(a, eta),
                        ignore_index = True, sort = False)
            ax.plot(df_Zmax["eta"], df_Zmax["Z"],
                    label = "$a/M = %g$" % a)
            ax.set_xlabel("$\eta/M^3$", fontsize = 15)
            ax.set_ylabel("$Z_{%g,%g,%g}^{Max}$" \
                            % (self.s, self.l, self.m), fontsize = 15)
            return ax
        else:
            df = self.select_data_spin(a)
            df = self.select_data_mu(mu, df)
            df_Zmax = pd.DataFrame()
            for eta in df["eta"].unique():
                df_Zmax = df_Zmax.append(self.get_Zmax_fixed_pars(a, eta, mu),
                        ignore_index = True, sort = False)
            ax.plot(df_Zmax["eta"], df_Zmax["Z"],
                    label = "$M\mu = %g$" % (mu))
            ax.set_xlabel("$\eta/M^3$", fontsize = 15)
            ax.set_ylabel("$Z_{%g,%g,%g}^{Max}$" \
                            % (self.s, self.l, self.m), fontsize = 15)
            return ax

    def plot_Zmax_compared(self, ax : plt.Axes, a : float,  mu : float = None):
        if not self.massive:
            df = self.select_data_spin(a)
            df_Zmax = pd.DataFrame()
            for eta in df["eta"].unique():
                df_Zmax = df_Zmax.append(self.get_Zmax_fixed_pars(a, eta),
                        ignore_index = True, sort = False)
            ax.plot(df_Zmax["eta"], df_Zmax["Z"]/df_Zmax[df_Zmax["eta"] == 0]["Z"].values[0] - 1,
                    label = "$a/M = %g$" % a)
            ax.set_xlabel("$\eta/M^3$", fontsize = 15)
            ax.set_ylabel("$Z_{%g,%g,%g}^{Max}/Z_{%g,%g,%g}^{Max, Kerr} - 1$" \
                            % (self.s, self.l, self.m, self.s, self.l, self.m), fontsize = 15)
            return ax
        else:
            df = self.select_data_spin(a)
            df = self.select_data_mu(mu, df)
            df_Zmax = pd.DataFrame()
            for eta in df["eta"].unique():
                df_Zmax = df_Zmax.append(self.get_Zmax_fixed_pars(a, eta, mu),
                        ignore_index = True, sort = False)
            ax.plot(df_Zmax["eta"], df_Zmax["Z"]/df_Zmax[df_Zmax["eta"] == 0]["Z"].values[0] - 1,
                    label = "$M\mu = %g$" % (mu))
            ax.set_xlabel("$\eta/M^3$", fontsize = 15)
            ax.set_ylabel("$Z_{%g,%g,%g}^{Max}/Z_{%g,%g,%g}^{Max, Kerr} - 1$" \
                            % (self.s, self.l, self.m, self.s, self.l, self.m), fontsize = 15)
            return ax

    def integrate_superradiant_spectra(self, a : float, eta : float, mu : float = None):
        if mu == None:
            condition = (self.data["a"] == a) & (self.data["eta"] == eta)
        else:
            condition = (self.data["a"] == a) & (self.data["eta"] == eta) & (self.data["mu"] == mu)
        df = self.data[condition]
        Z = df["Z"].values[0:49]
        omega = df["omega"].values[0:49]
        return simps(Z, omega)

    def plot_integral(self, ax, a : float, mu : float = None):
        if mu == None:
            condition = (self.data["a"] == a)
            label = "$a/M = %g$" % a
        else:
            condition = (self.data["a"] == a) & (self.data["mu"] == mu)
            label = "$M\mu = %g$" % mu
        I = []
        eta_list = self.data[condition]["eta"].unique()
        for eta in eta_list:
            I.append(self.integrate_superradiant_spectra(a, eta, mu))
        ax.plot(eta_list, I, label = label)
        ax.set_xlabel("$\eta/M^3$", fontsize = 15)
        ax.set_ylabel("$M\,I_{%g,%g,%g}$" % (self.s, self.l, self.m), fontsize = 15)
        return ax

# The `AmplificationFactors` class

Here we give some information about the AmplificationFactors class, used in the Jupyter notebook to manipulate and plot data.

- `AmplificationFactors` is instantiated as:
    `af = AmplificationFactors(s, l, m, massive = False)`\
    `s`, `l`, `m` are `int` representing the spin-weight of the field, the harmonic number and the azimuthal number, respectively, while `massive` is a `bool`, which is `True` when chosen data have &mu;&ne;0.
    Each `AmplificationFactors` object has attributes:
    - `massive`: `bool` which is `True` when data have &mu;&ne;0;
    - `s` is an `int` representing the spin-weight of the field;
    - `l` is an `int` representing the harmonic number;
    - `m` is an `int` representing the azimuthal number;
    - `data` is a `pandas.DataFrame` containing all data.
    
- `check_slm` method checks, at instantiation, if `s`, `l`, `m` and `massive` values, inserted by the user, are correct.

- `read_data` method reads data from file.

- `select_data_spin(a, data = pandas.DataFrame())` method selects data with spin parameter `a`.
    - `a` is a `float` representing the spin parameter of the black hole;
    - `data` is a `pandas.DataFrame` representing data to be selected; if empty, data is set to be equal to `af.data`. This is an optional parameter.
    
- `select_data_eta(eta, data = pandas.DataFrame())` method selects data with spin parameter `a`.
    - `eta` is a `float` representing the deformation parameter;
    - `data` is a `pandas.DataFrame` representing data to be selected; if empty, data is set to be equal to `af.data`. This is an optional parameter.
    
- `select_data_mu(mu, data = pandas.DataFrame())` method selects data with mass parameter `mu`.
    - `mu` is a `float` representing the mass parameter of the perturbing field;
    - `data` is a `pandas.DataFrame` representing data to be selected; if empty, data is set to be equal to `af.data`. This is an optional parameter.
    
- `get_Zmax(data)` method selects data with Z<sub>s, l, m</sub>=Z<sub>s, l, m</sub><sup>Max</sup>.
    - `data` is a `pandas.DataFrame` representing selected data.
    
- `get_Zmax_fixed_spin(a, data)` method returns a `pandas.DataFrame` containing Z<sub>s, l, m</sub><sup>Max</sup> for each value of the deformation parameter and for the selected spin parameter `a`.
    - `a` is a `float` representing the spin parameter;
    - `data`is a `pandas.DataFrame` representing data to be selected; if empty, data is set to be equal to `af.data`. This is an optional parameter.
    
- `get_Zmax_fixed_pars(a, eta, mu = None, data = pandas.DataFrame())` method returns a `pandas.DataFrame` containing Z<sub>s, l, m</sub><sup>Max</sup> for the selected parameters.
    - `a` is a `float` representing the spin parameter;
    - `eta` is a `float` representing the deformation parameter;
    - `mu` is a `float` representing the mass parameter of the perturbing field. This is an optional parameter to be used only if `af.massive == True`;
    - `data`is a `pandas.DataFrame` representing data to be selected; if empty, data is set to be equal to `af.data`. This is an optional parameter.
    
- `plot_spectra(ax, a, eta, mu = None)` method create plot of spectra with selected parameters in `ax` and returns `ax` with desired plot.
    - `ax` is a `matplotlib.pyplot.Axes` object where plots are put;
    - `a` is a `float` representing the spin parameter;
    - `eta` is a `float` representing the deformation parameter;
    - `mu` is a `float` representing the mass parameter of the perturbing field. This is an optional parameter to be used only if `af.massive == True`.
    
- `plot_Zmax(ax, a, mu = None)` method create plot of Z<sub>s, l, m</sub><sup>Max</sup> as a function of `eta` for selected spin and mass parameters and returns `ax` with desired plot.
    - `ax` is a `matplotlib.pyplot.Axes` object where plots are put;
    - `a` is a `float` representing the spin parameter;
    - `mu` is a `float` representing the mass parameter of the perturbing field. This is an optional parameter to be used only if `af.massive == True`.
    
- `plot_Zmax_compared(ax, a, mu = None)` method is the same as `plot_Zmax` but plots Z<sub>s, l, m</sub><sup>Max</sup>/Z<sub>s, l, m</sub><sup>Max, Kerr</sup>-1.
    - `ax` is a `matplotlib.pyplot.Axes` object where plots are put;
    - `a` is a `float` representing the spin parameter;
    - `mu` is a `float` representing the mass parameter of the perturbing field. This is an optional parameter to be used only if `af.massive == True`.
    
- `integrate_superradiant_spectra(a, eta, mu = None)` method returns the integral I<sub>s, l, m</sub> of the amplification fator over the superradiant frequency range.
    - `a` is a `float` representing the spin parameter;
    - `eta` is a `float` representing the deformation parameter;
    - `mu` is a `float` representing the mass parameter of the perturbing field. This is an optional parameter to be used only if `af.massive == True`.
    
- `plot_integral(ax, a, mu = None)` method plots the integral I<sub>s, l, m</sub> of the amplification factor over the superradiant frequency range as a function of the deformation parameter for the selected parameters.
    - `ax` is a `matplotlib.pyplot.Axes` object where plots are put;
    - `a` is a `float` representing the spin parameter;
    - `mu` is a `float` representing the mass parameter of the perturbing field. This is an optional parameter to be used only if `af.massive == True`.
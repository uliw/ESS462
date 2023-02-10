"""
A simple P-cycle model, based on chpater 8 of Modeling Methods for the marine
sciences. P-export depends only on P-concentration.
"""
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

F_r = 1.5e3  # River (weathering) flux of PO4 mol/s
V = 1.33e18  # Volume of ocean [m^3]

# some misc defintions
C_0 = [0]  # initial P concentration in ocean
time = 4e6  # model run time in years
year_to_seconds = 60 * 60 * 24 * 365.2425  # yr to s
t_span = (0, time * year_to_seconds)  #
p = (V, F_r)  # function arguments


def dCdt(t, C_0, V, F_r):
    """Calculate the change in concentration as function of time

    :param t: array, not used but must be present
    :param C: initial conditions (concentrations)
    :params V: Volume of ocean [m^3]
    :param F_r: River (weathering) flux of PO4 mol/s

    :returns: Change in Concentration
    """

    C = C_0[0]  # concentration
    F_b = C * 7e4  # determine burial flux
    dCdt = (F_r - F_b) / V  # dC/dt ocean

    return dCdt


# integrate dC/dt
result = solve_ivp(dCdt, t_span, C_0, args=p, max_step=t_span[1] / 50)
t = result.t / year_to_seconds / 1e6  # time in Myr
C = result.y.T  # concentration in mol
print(f"Final PO4 concentration {C[-1]} [mol/m^3]")
# rt = C * V / Fr mol/m^3 / mol/s -> s
rt = C[-1][0] * V / F_r
print(f"Rt = {rt/year_to_seconds:.2e}")

fig, ax = plt.subplots()
ax.plot(t, C)
ax.set_xlabel("Time [Myr]")
ax.set_ylabel("P concentration [mol/m^3]")
fig.tight_layout()
ax.grid()
plt.show()

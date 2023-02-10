"""
A simple P-cycle model, based on chapter 8 of Modeling Methods for the marine
sciences.

This model uses a 2 box ocean, and the export flux depends on the biological
productivity in the surface ocean
"""
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np

F_r = 1.5e3 * 1  # River (weathering) flux of PO4 mol/s
tau = 1  # P residency time in surface ocean in years
r = 0.99  # remineralization efficiency
V = [30e15, 1.33e18]  # Volume of surface and deep box m^3
THC = 30e6  # m^3/s
C_0 = [1.578e-04, 5.11e-03]  # initial P in surface and deep box
time = 10  # years
show_deep_ocean_data = False  # True or False

# some misc defintions
year_to_seconds = 60 * 60 * 24 * 365.2425  #
tau = tau * year_to_seconds  # convert to seconds
t_span = (0, time * year_to_seconds)  # model run time
p = (V, F_r, THC, tau, r)  # function arguments
max_step = t_span[1] / 100


def dCdt(t, C_0, V, F_r, THC, tau, r):
    """Calculate the change in concentration as function of time

    :param t: array, not used but must be present
    :param C_0: array of initial conditions (concentrations)

    :params V_s: Volume of surface ocean [m^3]

    :param V_d: Volume of deep ocean [m^3]
    :param F_r: River (weathering) flux of PO4 mol/s
    :param THC: Thermo haline circulation flux
    :param r: remineralization efficiency from 0 to 1

    :returns: array of dC/dt values

        The surface to deep ocean flux has two components:

            1. The dissolved PO4 that is removed via the thermohaline
               circulation, i.e., THC * C_s

            2. The particulate organic matter that sinks into the deep ocean and
               carries PO4 as as port of the cell structure.  This flux depends
               on the marine export productivity which is limited by the PO4
               concentratio, i.e., C_s* V / tau wher tau is a scaling factor
    """

    C_s = C_0[0]  # surface ocean concentration
    C_d = C_0[1]  # deep ocean concentration
    V_s = V[0]  # surface ocean volume
    V_d = V[1]  # deep ocean volume

    P = C_s * V_s / tau  # export productivity
    F_b = P * (1 - r)  # P burial
    F_s2d = THC * C_s  # Surface to deep ocean flux
    F_d2s = THC * C_d  # Deep to surface ocean flux

    # dC/dt surface ocean
    dCdt_s = (F_r + F_d2s - F_s2d - P) / V_s
    # dC/dt deep ocean
    dCdt_d = (F_s2d + P - F_d2s - F_b) / V_d

    return [dCdt_s, dCdt_d]


# integrate dC/dt
result = solve_ivp(
    dCdt,
    t_span,
    C_0,
    args=p,
    method="BDF",
    max_step=max_step,
)
t = result.t / year_to_seconds / 1e6  # time in Myr
C_s = result.y[0]  # surface box concentration in mol
C_d = result.y[1]  # deep box concentration in mol


rt = C_d[-1] * V[1] / F_r
D_C = C_d[-1] - C_s[-1]
Rt = rt / year_to_seconds
F_POP = C_s[-1] * V[0] / tau
print("C_s mol/m^3, C_d mol/m^3, D_C, F_POP, Rt")
print(f"{C_s[-1]:.2e}, {C_d[-1]:.2e}, {D_C:.2e}, {F_POP:.2e}, {Rt:.2e}\n")

fig, ax = plt.subplots()
ax.plot(t, C_s, color="C0", label="Surface Ocean")
if show_deep_ocean_data:
    ax.plot(t, C_d, color="C1", label="Deep Ocean")
ax.set_xlabel("Time [Myr]")
ax.set_ylabel("P concentration [mol/m^3]")
ax.legend()
ax.grid()
fig.tight_layout()
plt.show(block=False)

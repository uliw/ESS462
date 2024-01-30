def check_convergence(C, c):
    """check if first derivative is smaller than c
    C: Array with numbers
    c: convergence criterion
    """
    import numpy as np

    lc = len(C)
    lc2 = int(lc / 2)
    diff = np.diff(C)
    for i, e in enumerate(diff):
        if e < lc and i > lc2:
            break
    return i


def dCdt1(t, C_0, V, F_r):
    """Calculate the change in concentration as function of time
    No export productivity

    :param t: array, not used but must be present
    :param C: initial conditions (concentrations)
    :params V: Volume of ocean [m^3]
    :param F_r: River (weathering) flux of PO4 mol/s
    :returns: Change in Concentration
    """

    C = C_0[0]
    F_b = C * 7e4  # determine burial flux
    dCdt = (F_r - F_b) / V  # dC/dt ocean

    return dCdt


def dCdt2(t, C_0, V, F_r, THC, tau, r):
    """Calculate the change in concentration as function of time
    including particulate p transport by export productivity

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

    C_s = C_0[0] / 1000  # surface ocean concentration
    C_d = C_0[1] / 1000  # deep ocean concentration
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


def po41(time, C_0, V, F_r, dpi, fn, plot=True):
    """Calculate P concentration over time. Iput flux is
    fixed, output flux depends on concentration
    """
    from scipy.integrate import solve_ivp
    import matplotlib.pyplot as plt
    import numpy as np

    year_to_seconds = 60 * 60 * 24 * 365.2425  # yr to s
    t_span = (0, time * year_to_seconds)  #
    p = (V, F_r)  # function arguments
    C_0 = np.array(C_0) / 1000  # convert to mol/m^3

    # integrate dC/dt
    result = solve_ivp(dCdt1, t_span, C_0, args=p, max_step=t_span[1] / 50)
    t = result.t / year_to_seconds / 1e6  # time in Myr
    C = result.y.T * 1e3  # concentration in umol/l
    print(f"Final PO4 concentration {C[-1][0]:.2f} [umol/l]")
    i = check_convergence(C[:, 0], 1e-3)
    if i >= len(C[:, 0]):
        print("Failed to reach steady state. Double the run time")
    else:
        print(f"Steady state reached after {t[i]:.1e} Myrs")

    if plot:
        # plot data
        fig, ax = plt.subplots()
        ax.plot(t, C)
        ax.set_xlabel("Time [Myr]")
        ax.set_ylabel(r"P concentration [$\mu$mol/l]")
        fig.tight_layout()
        ax.grid()
        plt.show()
        fig.set_dpi(dpi)
        fig.savefig(fn)

    return C[:, 0]


def po42(time, C_0, V, F_r, dpi, fn, tau, r, THC):
    """Calculate P concentration over time. Iput flux is
    fixed, output flux depends on concentration. Consider
    transfer of particular P by export productivity.
    """
    from scipy.integrate import solve_ivp
    import matplotlib.pyplot as plt

    year_to_seconds = 60 * 60 * 24 * 365.2425  #
    tau = tau * year_to_seconds  # convert to seconds
    t_span = (0, time * year_to_seconds)  # model run time
    p = (V, F_r, THC, tau, r)  # function arguments
    max_step = t_span[1] / 100
    show_deep_ocean_data = True  # True or False

    result = solve_ivp(  # integrate dC/dt
        dCdt2,
        t_span,
        C_0,
        args=p,
        method="BDF",
        max_step=max_step,
    )
    t = result.t / year_to_seconds / 1e6  # time in Myr
    C_s = result.y[0]  # surface box concentration in mol/m^3
    C_d = result.y[1]  # deep box concentration in mol/m^3

    print(f"tau = C_s = {C_s[-1]:.2e} mol/m^2, C_d = {C_d[-1]:.2e} mol/m^3")
    fig, ax = plt.subplots()
    ax.plot(t, C_s, color="C0", label="Surface Ocean")
    if show_deep_ocean_data:
        ax.plot(t, C_d, color="C1", label="Deep Ocean")
        ax.set_xlabel("Time [Myr]")
        ax.set_ylabel("P concentration [\mu mol/l]")
        ax.legend()
        ax.grid()
        fig.tight_layout()
        plt.show()
        fig.set_dpi(dpi)
        fig.savefig(fn)


if __name__ == "__main__":
    import po4

    F_r = 1e3  # River (weathering) flux of PO4 mol/s
    V = 1.33e18  # Volume of ocean [m^3]
    C_0 = [6]  # initial P concentration in ocean umol/l
    time = 4e6  # model run time in years
    plot_dpi = 300  #
    figure_name = "po41.png"

    # ---------- no user serviceable parts below ---------- #
    print(f"Weathering Flux = {F_r:.2e} [mol/s])")
    print(f"Volume = {V:.2e} m^3"),
    print(f"C_0 = {C_0[0]:.2f} umol/l)")
    C = po4.po41(time, C_0, V, F_r, plot_dpi, figure_name, plot=False)

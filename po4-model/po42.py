def check_convergence(C, c):
    """check if first derivative is smaller than c, and return
    the index position where this condition is achieved


    C: Array with numbers
    c: convergence criterion
    """
    import numpy as np

    lc = len(C)  # length of vector
    lc2 = int(lc / 2)  # 1/5 length
    diff = np.diff(C)
    index = lc
    for i, e in enumerate(diff):
        if i > lc2:  # ignore the first half
            if e < c:
                index = i
                break

    return index


def check_halfway(C, t):
    """check if first derivative is smaller than c
    C: Array with numbers
    c: convergence criterion
    """
    import numpy as np

    c_min = np.min(C)
    c_max = np.max(C)
    h = (c_max - c_min) / 2

    for i, e in enumerate(C):
        if e >= h + c_min:
            break

    return t[i]


def dCdt2(t, C_0, V, F_r, THC, tau, eff, burial):
    """Calculate the change in concentration as function of time
    including particulate p transport by export productivity

    :param t: array, not used but must be present
    :param C_0: array of initial conditions (concentrations)
    :params V_s: Volume of surface ocean [m^3]
    :param V_d: Volume of deep ocean [m^3]
    :param F_r: River (weathering) flux of PO4 mol/s
    :param THC: Thermo haline circulation flux
    :param tau: Residence time in the surface water
    :param eff: uptake efficiency
    :param burial: remineralization efficiency from 0 to 1
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

    F_s2d = C_s * THC  # Surface to deep ocean flux
    F_d2s = C_d * THC  # Deep to surface ocean flux

    if eff is None:
        P = C_s * V_s / tau  # export productivity
    elif tau is None:
        P = F_d2s * eff / 100
    else:
        raise ValueError("This should never happen")

    F_b = P * burial / 100  # P burial flux
    # dC/dt surface ocean
    dCdt_s = (F_r + F_d2s - F_s2d - P) / V_s
    # dC/dt deep ocean
    dCdt_d = (F_s2d + P - F_d2s - F_b) / V_d

    return [dCdt_s, dCdt_d]


def po42(time, C_0, V, F_r, dpi, fn, tau, eff, burial, THC, export=False, ss=True):
    """Calculate P concentration over time. Iput flux is
    fixed, output flux depends on concentration. Consider
    transfer of particular P by export productivity.
    """
    from scipy.integrate import solve_ivp
    import matplotlib.pyplot as plt
    import numpy as np

    year_to_seconds = 60 * 60 * 24 * 365
    THC = THC * 1e6  # Sv to m^3/s
    if tau is None:
        pass
    else:
        tau = tau * year_to_seconds  # convert to seconds
    t_span = (0, time * year_to_seconds)  # model run time
    p = (V, F_r, THC, tau, eff, burial)  # function arguments
    max_step = t_span[1] / 100
    tpoints = np.linspace(0, t_span[1], 100)
    show_deep_ocean_data = True  # True or False

    result = solve_ivp(  # integrate dC/dt
        dCdt2,
        t_span,
        C_0,
        args=p,
        t_eval=tpoints,
        method="BDF",
        max_step=max_step,
    )
    t = result.t / year_to_seconds / 1e6  # time in Myr
    C_s = result.y[0] * 1000  # surface box concentration in mmol/m^3
    C_d = result.y[1] * 1000  # deep box concentration in mmol/m^3

    if ss:
        i = check_convergence(C_s, 1e-3)
        if i >= len(C_s):
            raise ValueError("\nFailed to reach steady state, increase your run time\n")
        else:
            print(f"Steady state reached after {t[i]:.1e} Myrs")

    print(f"C_s = {C_s[-1]:.1f} mmol/m^2, C_d = {C_d[-1]:.1f} mmol/m^3")
    fig, ax = plt.subplots()
    if export:
        if eff is None:
            P = C_s * V[0] / tau  # export productivity
        elif tau is None:
            P = C_d * THC * eff / 100
        else:
            raise ValueError("This should never happen")
        ax.plot(t*1e6, P, color="C0", label="Export Production")
        ax.set_xlabel("Time [yr]")
        ax.set_ylabel("Export Production [mol/yr]")
        ax.ticklabel_format(useOffset=False)
    else:
        ax.plot(t, C_s, color="C0", label="Surface Ocean")
        ax.plot(t, C_d, color="C1", label="Deep Ocean")
        ax.set_xlabel("Time [Myr]")
        ax.set_ylabel("P concentration [$\mu$mol/l]")
        # ax.plt.legend()
        
    ax.legend()
    ax.grid()
    fig.tight_layout()
    plt.show()
    fig.set_dpi(dpi)
    fig.savefig(fn)

    return [C_s, C_d]


if __name__ == "__main__":
    import po42

    F_r = 1.43e3  # River (weathering) flux of PO4 mol/s
    V = [30e15, 1.33e18]  # Volume of surface and deep box m^3 [1, 3]
    THC = 19  # Sverdrup
    tau = 1  # P residency time in surface ocean in years
    uptake_efficiency = None
    burial = 1  # P-burial in %
    C_0 = [0, 0]  # initial P in surface and deep box
    time = 8e6  # run time in years
    plot_dpi = 120  #
    figure_name = "po42.png"

    # ---------- no user serviceable parts below ---------- #
    fig = po42.po42(
        time,
        C_0,
        V,
        F_r,
        plot_dpi,
        figure_name,
        tau,
        uptake_efficiency,
        burial,
        THC,
    )

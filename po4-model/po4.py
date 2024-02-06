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
    result = solve_ivp(
        dCdt1,
        t_span,
        C_0,
        args=p,
        max_step=t_span[1] / 50,
    )
    t = result.t / year_to_seconds / 1e6  # time in Myr
    C = result.y.T * 1e3  # concentration in umol/l
    print(f"Final PO4 concentration {C[-1][0]:.2f} [umol/l]")
    t2 = check_halfway(C, t)
    print(f"Halfway mark at t = {t2:.2f} [Myr]")

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


if __name__ == "__main__":
    import po4

    F_r = 1e3  # River (weathering) flux of PO4 mol/s
    V = 1.33e18  # Volume of ocean [m^3]
    C_0 = [6]  # initial P concentration in ocean umol/l
    time = 2e6  # model run time in years
    plot_dpi = 120  #
    figure_name = "po41.png"

    # ---------- no user serviceable parts below ---------- #
    print(f"Weathering Flux = {F_r:.2e} [mol/s])")
    print(f"Volume = {V:.2e} m^3"),
    print(f"C_0 = {C_0[0]:.2f} umol/l)")
    C = po42.po42(time, C_0, V, F_r, plot_dpi, figure_name, plot=False)
    # convert flux into mol/kyear/liter
    f = F_r * 60 * 60 * 24 * 365 * 1e6
    Rt = V * C[-1] / f
    print(f"Rt = {Rt:.0f} [kyr]")

from model_calibration import *
from scipy.optimize import curve_fit

def Pmodel(t, *pars):
    """ Return reservoir pressure at a set of times for a given set of parameters.

        Parameters
        ----------
        t : array-like
            Vector of time (year).
        pars : list
            List of parameters to be fitted [a, b1, c].

        Returns
        -------
        P : array-like
            Vector of reservoir pressure (bar).
    """
    q = interpolate_mass_extraction_rate(t)
    dqdt = dq_dt(t, q)

    t0 = t[0]
    t1 = t[-1]
    dt = t[1]-t[0]

    (a, b, c) = pars
    P0 = 56.26
    _, P = solve_reservoir_ode(reservoir_ode, t0, t1, dt, P0, q, dqdt, [a,b,c,P0])
    return P

def Umodel(t, *pars):
    """ Return subsidence at a set of times for a given set of parameters.

        Parameters
        ----------
        t : array-like
            Vector of time (year).
        pars : list
            List of parameters to be fitted [b2, mv].

        Returns
        -------
        U : Vector of subsidence (m).

    """
    q = interpolate_mass_extraction_rate(t)
    dqdt = dq_dt(t, q)

    t0 = t[0]
    t1 = t[-1]
    dt = t[1] - t[0]

    # fixed parameters
    a = 0.0015377621974240604
    b1 = 0.06624924440742241
    c = 0.007883715689885993
    P0 = 56.26
    Pm0 = 56.26
    L = 100

    a = 0.0016753161225277887
    b1 = 0.07412774351741645
    c = 0.007630305388155743

    # parameters to be fitted
    (b2, mv) = pars

    # solve for subsidence
    _, P = solve_reservoir_ode(reservoir_ode, t0, t1, dt, P0, q, dqdt, [a, b1, c, P0])
    _, Pm = solve_mudstone_ode(mudstone_ode, t0, t1, dt, P, Pm0, [b2])
    U = subsidence_eqn(Pm, Pm0, mv, L)

    return U

def calibrate_reservoir_pressure():
    """ Calibrate parameters to best fit the reservoir pressure data.

        Returns
        -------
        a : float
            Fitted lumped parameter for forced term.
        b : float
            Fitted lumped parameter for recharge term.
        c : float
            Fitted lumped parameter for slow drainage term.
    """
    to,Po = load_pressure_data()

    # interpolate mass extraction rate data and find derivative
    q = interpolate_mass_extraction_rate(to)
    dqdt = dq_dt(to, q)

    p0 = [1.62e-3,7.15e-2,0.007]  # original guess
    constants = curve_fit(Pmodel, to[:-10], Po[:-10], p0)
    (a, b, c) = constants[0]
    P0 = 56.26
    ti,Pi = solve_reservoir_ode(reservoir_ode, to[0], to[-1], to[1]-to[0], P0, q, dqdt, [a,b,c,P0])

    fig, ax = plt.subplots()

    # plot the observations and the numerical data
    ax.plot(to, Po, '.', ti, Pi)
    ax.set_ylabel('reservoir pressure [bars]')
    ax.set_xlabel('time [years]')
    ax.legend(("observations", "model"))
    plt.show()

    # best fit
    # a = 0.0015377621974240604
    # b1 = 0.06624924440742241
    # c = 0.007883715689885993

    return a, b, c

def calibrate_subsidence():
    """ Calibrate parameters to best fit to the subsidence data.

        Returns
        -------
        b2 : float
            Fitted lumped parameter for recharge term.
        mv : float
            Fitted compressibility of mudstone.
    """
    to, Uo = load_subsidence_data()

    # interpolate subsidence at defined time intervals
    t = time_range(to[0], to[-1], 1)
    U = np.interp(t, to, Uo)

    p0 = [0.0355,0.0076]  # initial guess
    constants = curve_fit(Umodel, t[:-10], U[:-10], p0)
    (b2, mv) = constants[0]  # extract fitted parameters

    q = interpolate_mass_extraction_rate(t)
    dqdt = dq_dt(t, q)

    t0 = t[0]
    t1 = t[-1]
    dt = t[1] - t[0]

    # fixed parameters
    a = 0.0015377621974240604
    b1 = 0.06624924440742241
    c = 0.007883715689885993
    P0 = 56.26
    Pm0 = 56.26
    L = 100

    a = 0.0016753161225277887
    b1 = 0.07412774351741645
    c = 0.007630305388155743

    # solve for subsidence
    t, P = solve_reservoir_ode(reservoir_ode, t0, t1, dt, P0, q, dqdt, [a, b1, c, P0])
    t, Pm = solve_mudstone_ode(mudstone_ode, t0, t1, dt, P, Pm0, [b2])
    subsidence = subsidence_eqn(Pm, Pm0, mv, L)

    fig, ax = plt.subplots()

    # plot the observations and the numerical data
    ax.plot(to, Uo, '.', t, subsidence)
    ax.set_ylabel('subsidence [m]')
    ax.set_xlabel('time, [yr]')
    ax.legend(("observations", "model"))
    plt.show()

    # best fit
    # b2 = 0.031776570956506885
    # mv = 0.007905952772346621

    return b2, mv


if __name__ == "__main__":
    # pars = calibrate_reservoir_pressure()
    # print(pars)
    pars = calibrate_subsidence()
    print(pars)
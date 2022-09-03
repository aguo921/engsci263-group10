import numpy as np


def reservoir_ode(t, P, q, dqdt, a, b, c, P0):
    """ Reservoir pressure ODE.

        Parameters
        ----------
        t : float
            Time (year).
        P : float
            Fluid pressure in reservoir (bar).
        q : float
            Mass extraction rate (kg/s).
        dqdt : float
            Rate of change of mass extraction rate (kg/s^2).
        a : float
            Lumped parameter for forcing term.
        b : float
            Lumped parameter for recharge term.
        c : float
            Lumped parameter for slow drainage term.
        P0 : float
            Ambient reservoir pressure (bar).

        Returns
        -------
        dPdt : float
            Rate of change of reservoir pressure (bar/year).
    """
    return -a*q - b*(P-P0) - c*dqdt


def mudstone_ode(t, Pm, P, b):
    """ Mudstone pressure ODE.

        Parameters
        ----------
        t : float
            Time (year).
        Pm : float
            Fluid pressure in mudstone (bar).
        P : float
            Fluid pressure in reservoir (bar).
        b : float
            Lumped parameter for recharge term.

        Returns
        -------
        dPdt : float
            Rate of change of the mudstone pressure (bar/year).
    """
    return -b*(Pm-P)


def subsidence_eqn(P, P0, mv, L):
    """ Subsidence equation.

        Parameters
        ----------
        P : np.array
            Vector of fluid pressure in mudstone (bar).
        P0 : float
            Ambient mudstone pressure (bar).
        mv : float
            Bulk compressibility of mudstone (bar^-1).
        L : float
            Thickness of mudstone (m).

        Returns
        -------
        U : np.array
            Vector of subsidence at the centre of the bowl (m).
    """
    return mv*L*(P0-P)


def time_range(t0, t1, dt):
    """
        Return vector of times with starting time t0 and ending time t1 (inclusive), with time step dt.
    """
    t = np.arange(t0, t1, dt)
    if t[-1] == t1 - dt:
        t = np.append(t, t1)
    return np.array(t)


def solve_reservoir_ode(f, t0, t1, dt, P0, q, dqdt, pars):
    """ Solve the reservoir pressure ODE numerically.

        Parameters
        ----------
        f : callable
            Reservoir pressure ODE.
        t0 : float
            Initial time (year).
        t1 : float
            Final time (year).
        dt : float
            Time step (year).
        P0 : float
            Initial reservoir pressure (bar).
        q : np.array
            Vector of interpolated mass extraction values (kg/s).
        dqdt : np.array
            Vector of derivative of mass extraction interpolation function (kg/s^2).
        pars : list
            List of reservoir pressure ODE parameters.

        Returns
        -------
        t : np.array
            Vector of time (year).
        P : np.array
            Vector of reservoir pressure (bar).

        Notes
        -----
        ODE should be solved using the Improved Euler Method.

        Assume that ODE function f takes the following inputs, in order:
            1. independent variable
            2. dependent variable
            3. forcing term, q
            4. rate of change of forcing term, dqdt
            5. all other parameters

        Assume q and dqdt vectors correspond to a time vector with start t0, end t1 and step dt.
    """
    # initialise time and pressure vectors
    t = time_range(t0, t1, dt)
    P = [P0]

    # iteratively perform Improved Euler's method
    for i in range(len(t)-1):
        f0 = f(t[i], P[i], q[i], dqdt[i], *pars)
        f1 = f(t[i+1], P[i]+dt*f0, q[i+1], dqdt[i+1], *pars)
        P.append(P[i] + dt / 2 * (f0 + f1))

    return t, np.array(P)


def solve_mudstone_ode(f, t0, t1, dt, P, Pm0, pars):
    """ Solve the reservoir pressure ODE numerically.

        Parameters
        ----------
        f : callable
            Mudstone pressure ODE.
        t0 : float
            Initial time (year).
        t1 : float
            Final time (year).
        dt : float
            Time step (year).
        P : np.array
            Vector of reservoir pressure values (bar).
        Pm0 : float
            Initial mudstone pressure (bar).
        pars : list
            List of mudstone pressure ODE parameters.

        Returns
        -------
        t : np.array
            Vector of time (year).
        Pm : np.array
            Vector of mudstone pressure (bar).

        Notes
        -----
        ODE should be solved using the Improved Euler Method.

        Assume that ODE function f takes the following inputs, in order:
            1. independent variable
            2. dependent variable
            3. forcing term, P
            4. all other parameters

        Assume P vector corresponds to a time vector with start t0, end t1 and step dt.
    """
    # initialise time and pressure vectors
    Pm = [Pm0]
    t = time_range(t0, t1, dt)

    # iteratively perform Improved Euler's method
    for i in range(len(t) - 1):
        f0 = f(t[i], Pm[i], P[i], *pars)
        f1 = f(t[i+1], Pm[i]+dt*f0, P[i+1], *pars)
        Pm.append(Pm[i] + dt / 2 * (f0 + f1))

    return t, np.array(Pm)

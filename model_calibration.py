from model_functions import *
from matplotlib import pyplot as plt

def load_pressure_data():
    ''' Returns time and pressure from sb_pres.txt

        Parameters:
        -----------
        none

        Returns:
        --------
        t : array-like
            Vector of times (years) at which measurements were taken.
        p : array-like
            Vector of Fluid pressure in reservoir (bar).
    '''

    t, p = np.genfromtxt("sb_pres.txt", delimiter=",", skip_header=1).T
    return t, p


def load_subsidence_data():
    ''' Returns time and subsidenc from sb_pres.txt

        Parameters:
        -----------
        none

        Returns:
        --------
        t : array-like
            Vector of times (years) at which measurements were taken.
        s : array-like
            Vector of subsidence in reservoir (bar).
    '''

    t, s = np.genfromtxt("sb_disp.txt", delimiter=",", skip_header=1).T
    return t, s


def interpolate_mass_extraction_rate(t):
    ''' Return mass extraction rate q.

        Parameters:
        -----------
        t : array-like
            Vector of times at which to interpolate the mass extraction rate.

        Returns:
        --------
        q : array-like
            Mass extraction rate interpolated at t.
    '''
    t0, qt = np.genfromtxt("sb_mass.txt", delimiter=",", skip_header=1).T

    q = np.interp(t, t0, qt)
    return q


def dq_dt(t, q):
    ''' Return the rate of change of mass extraction rate.

        Parameters:
        -----------
        t : array-like
            Vector of times at which to interpolate the mass extraction rate.
        q : array-like
            Mass extraction rate interpolated at t.

        Returns:
        --------
        dqdt : array-like
            Vector of the rate of change of the mass extraction rate.

    '''
    dqdt = np.divide(np.diff(q), np.diff(t))
    dqdt = np.append(dqdt, dqdt[-1])
    return dqdt

def plot_reservoir_pressure_model():
    ''' Plot the Reservoir Pressure ODE over top of the data.

        Parameters:
        -----------
        none

        Returns:
        --------
        none

        Notes:
        ------
        This function called within if __name__ == "__main__":
    '''
    time, pres = load_pressure_data()
    t0 = time[0]
    t1 = time[-1]
    dt = 0.1
    t = time_range(t0, t1, dt)
    q = interpolate_mass_extraction_rate(t)
    dqdt = dq_dt(t, q)

    a = 0.0015377621974240604  # lumped parameter for forced extraction term
    b1 = 0.06624924440742241  # lumped parameter for recharge term (into reservoir)
    c = 0.007883715689885993  # lumped parameter for slow drainage term
    P0 = 56.26  # initial/ambient reservoir pressure

    # find numerical solution for reservoir pressure
    t, p = solve_reservoir_ode(reservoir_ode, t0, t1, dt, P0, q, dqdt, [a, b1, c, P0])

    fig, ax = plt.subplots()

    # plot the observations and the numerical data
    ax.plot(time, pres, '.', t, p)
    ax.set_ylabel('reservoir pressure [bars]')
    ax.set_xlabel('time [years]')
    ax.legend(("observations", "model"))
    plt.show()


def plot_model():
    ''' Plot the reservoir pressure and subsidence over top of the data.

        Parameters:
        -----------
        none

        Returns:
        --------
        none

        Notes:
        ------
        This function called within if __name__ == "__main__":
    '''
    a = 0.0015377621974240604  # lumped parameter for forced extraction term
    b1 = 0.06624924440742241  # lumped parameter for recharge term (into reservoir)
    c = 0.007883715689885993  # lumped parameter for slow drainage term
    P0 = 56.26  # initial/ambient reservoir pressure
    b2 = 0.031776570956506885  # lumped parameter for recharge term (into mudstone)
    mv = 0.007905952772346621  # compressibility of mudstone
    L = 100  # thickness of mudstone
    Pm0 = 56.26  # initial/ambient mudstone pressure

    a = 0.0016753161225277887
    b1 = 0.07412774351741645
    c = 0.007630305388155743
    b2 = 0.024184021965282205
    mv = 0.009539371878306748

    time_pres, pres = load_pressure_data()
    time_sub, sub = load_subsidence_data()
    t0 = time_pres[0]
    t1 = time_sub[-1]
    dt = 1

    t = time_range(t0, t1, dt)
    q = interpolate_mass_extraction_rate(t)
    dqdt = dq_dt(t, q)

    _, P = solve_reservoir_ode(reservoir_ode, t0, t1, dt, P0, q, dqdt, [a, b1, c, P0])
    _, Pm = solve_mudstone_ode(mudstone_ode, t0, t1, dt, P, Pm0, [b2])

    U = subsidence_eqn(Pm, Pm0, mv, L)

    fig, (ax1, ax2) = plt.subplots(ncols=2)

    # plot pressure observations and numerical solution
    ax1.plot(time_pres, pres, '.', t, P)
    ax1.set_ylabel('reservoir pressure [bars]')
    ax1.set_xlabel('time [years]')
    ax1.legend(("observations", "model"))

    # plot subsidence observations and numerical solution
    ax2.plot(time_sub, sub, '.', t, U)
    ax2.set_ylabel('subsidence [m]')
    ax2.set_xlabel('time [years]')
    ax2.legend(("observations", "model"))

    fig.set_size_inches(10, 6)
    plt.show()


def plot_misfit(slow_drainage=False):
    if slow_drainage:
        a = 0.0015377621974240604  # lumped parameter for forced extraction term
        b1 = 0.06624924440742241  # lumped parameter for recharge term (into reservoir)
        c = 0.007883715689885993  # lumped parameter for slow drainage term
    else:
        a = 0.0029363380044929903
        b1 = 0.13757916629903563
        c = 0
    P0 = 56.26  # initial/ambient reservoir pressure

    time, pres = load_pressure_data()
    t0 = time[0]
    t1 = time[-1]
    dt = time[1] - time[0]

    t = time_range(t0, t1, dt)
    q = interpolate_mass_extraction_rate(t)
    dqdt = dq_dt(t, q)

    _, P = solve_reservoir_ode(reservoir_ode, t0, t1, dt, P0, q, dqdt, [a, b1, c, P0])

    misfit = P - pres

    fig, (ax1, ax2) = plt.subplots(ncols=2)

    # plot pressure observations and numerical solution
    ax1.plot(time, pres, '.', t, P)
    ax1.set_ylabel('reservoir pressure [bars]')
    ax1.set_xlabel('time [years]')
    ax1.legend(("observations", "model"))

    # plot subsidence observations and numerical solution
    ax2.plot(time, misfit, 'x', time, np.zeros(len(time)), '--')
    ax2.set_ylabel('pressure misfit [bars]')
    ax2.set_xlabel('time [years]')

    fig.set_size_inches(10, 6)
    plt.show()


if __name__ == "__main__":
    plot_model()
    # plot_misfit(True)
    # plot_misfit(False)


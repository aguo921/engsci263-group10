from model_calibration import *
from matplotlib import pyplot as plt

def forecast(q1):
    ''' Forecast the subsidence.

        Parameters:
        -----------
        q1: list
            List of future mass extraction rates.

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

    time_sub, sub = load_subsidence_data()
    time_pres, pres = load_pressure_data()

    t0 = time_pres[0]
    t1 = 2012
    dt = 1

    ti = time_range(t0, t1, dt)
    tP = time_range(t1+dt,2042,dt)
    t = np.concatenate((ti,tP))

    q0 = interpolate_mass_extraction_rate(ti)

    subsidence = []
    pressure = []

    for i in range(len(q1)):
        qp = q1[i]*np.ones(len(tP))
        q=np.concatenate((q0,qp))
        dqdt = dq_dt(t, q)

        _, P = solve_reservoir_ode(reservoir_ode, t0, 2042, dt, P0, q, dqdt, [a, b1, c, P0])
        _, Pm = solve_mudstone_ode(mudstone_ode, t0, 2042, dt, P, Pm0, [b2])

        pressure.append(P)
        subsidence.append(subsidence_eqn(Pm, Pm0, mv, L))


    fig, (ax1, ax2) = plt.subplots(ncols=2)

    # plot forecasts
    for i in range(len(q1)):
        ax1.plot(t, pressure[i], '-', label=f"q={q1[i]}kg/s")
        ax2.plot(t, subsidence[i], '-', label=f"q={q1[i]}kg/s")

    # plot observations
    ax1.plot(time_pres, pres, '.', label="observations")
    ax2.plot(time_sub, sub, '.', label="observations")

    # plot best-fit solution
    ax1.plot(ti, pressure[0][:len(ti)], label="best-fit model")
    ax2.plot(ti, subsidence[0][:len(ti)], label="best-fit model")

    ax1.set_ylabel('pressure [bar]')
    ax2.set_ylabel('subsidence [m]')
    ax1.set_xlabel('time [yrs]')
    ax2.set_xlabel('time [yrs]')
    ax1.legend()
    ax2.legend()

    fig.set_size_inches(12, 6)
    plt.show()




if __name__ == "__main__":
    forecast([1250, 900, 600, 0])
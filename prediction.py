from model_calibration import *
from matplotlib import pyplot as plt

def forecast(q1):
    ''' Plot the forecast of the subsidence and the reservoir pressure for the next 30 years when given the list of
        proposed future mass extraction rates.

        Parameters:
        -----------
        q1: list
            List of proposed future mass extraction rate. The mass extraction rates given are held constant.

        Returns:
        --------
        none
    '''
    a = 0.0015377621974240604  # lumped parameter for forced extraction term
    b1 = 0.06624924440742241  # lumped parameter for recharge term (into reservoir)
    b2 = 0.031776570956506885  # lumped parameter for recharge term (into mudstone)
    c = 0.007883715689885993  # lumped parameter for slow drainage term

    mv = 0.007905952772346621  # compressibility of mudstone
    L = 100  # thickness of mudstone

    P0 = 56.26  # initial/ambient reservoir pressure
    Pm0 = 56.26  # initial/ambient mudstone pressure

    # load subsidence and time data from sb_disp.txt
    time_sub, sub = load_subsidence_data()

    # load pressure and time data from sb_pres.txt
    time_pres, pres = load_pressure_data()

    t0 = time_pres[0]  # first time data point
    t1 = time_pres[-1]  # last time data point
    dt = 1  # step size for time data

    ti = time_range(t0, t1, dt)  # time range for given period
    tP = time_range(t1+dt, 2042, dt)  # time range for the prediction period
    t = np.concatenate((ti, tP))  # total time range

    # interpolate the mass extraction rate within the given period
    q0 = interpolate_mass_extraction_rate(ti)

    # create empty lists for the subsidence and pressure values
    subsidence = []
    pressure = []

    # use a for loop to predict the reservoir pressure, mudstone pressure and subsidence under four different conditions
    for i in range(len(q1)):
        # create a vector of constant mass extraction rate for the prediction period
        qp = q1[i]*np.ones(len(tP))

        # get the vector of the rate of change in mass extraction rate for the full time range
        q = np.concatenate((q0, qp))

        # get the vector of the rate of the rate of change in mass extraction rate
        dqdt = dq_dt(t, q)

        # compute the reservoir pressure
        _, P = solve_reservoir_ode(reservoir_ode, t0, 2042, dt, P0, q, dqdt, [a, b1, c, P0])

        # compute the mudstone pressure
        _, Pm = solve_mudstone_ode(mudstone_ode, t0, 2042, dt, P, Pm0, [b2])

        # compute the subsidence
        U = subsidence_eqn(Pm, Pm0, mv, L)

        # append computed reservoir pressure and subsidence to lists
        pressure.append(P)
        subsidence.append(U)

    # create figure of two graphs (one for pressure, one for subsidence)
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

    # label title and axes
    ax1.set_title('Reservoir pressure predictions')
    ax2.set_title('Subsidence predictions')
    ax1.set_ylabel('pressure [bar]')
    ax2.set_ylabel('subsidence [m]')
    ax1.set_xlabel('time [yrs]')
    ax2.set_xlabel('time [yrs]')

    # create legend
    ax1.legend()
    ax2.legend()

    # display figure
    fig.tight_layout()
    fig.set_size_inches(12, 6)
    plt.show()


if __name__ == "__main__":
    forecast([1250, 900, 600, 0])
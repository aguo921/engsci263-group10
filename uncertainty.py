from model_calibration import *
from scipy.optimize import curve_fit

def Pmodel(t, *pars):
    """ Return reservoir pressure at a set of times for a given set of parameters.

        Parameters
        ----------
        t : array-like
            Vector of time (year).
        pars : list
            List of parameters to be fitted [a, b1].

        Returns
        -------
        P : array-like
            Vector of reservoir pressure (bar).
    """
    q = interpolate_mass_extraction_rate(t)
    dqdt = dq_dt(t, q)

    t0 = t[0]
    t1 = t[-1]
    dt = t[1] - t[0]

    # fixed parameters
    c = 0.007883715689885993
    P0 = 56.26

    # parameters to be fitted
    (a, b1) = pars

    _, P = solve_reservoir_ode(reservoir_ode, t0, t1, dt, P0, q, dqdt, [a, b1, c, P0])
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

    # parameters to be fitted
    (b2, mv) = pars

    _, P = solve_reservoir_ode(reservoir_ode, t0, t1, dt, P0, q, dqdt, [a, b1, c, P0])
    _, Pm = solve_mudstone_ode(mudstone_ode, t0, t1, dt, P, Pm0, [b2])
    U = subsidence_eqn(Pm, Pm0, mv, L)
    return U


def best_fit_pressure(v):
    """ Calculate best fit parameters and covariance to reservoir pressure data.

        Parameters
        ----------
        v : float
            Observational error of the reservoir pressure data.

        Returns
        -------
        pars : tuple of float
            Tuple of the best fit parameters a and b1 respectively.
        cov : array-like
            Covariance matrix.
    """
    to, Po = load_pressure_data()

    # interpolate subsidence at defined time intervals
    t = time_range(to[0], to[-1], 1)
    U = np.interp(t, to, Po)

    p0 = [1.62e-3,7.15e-2]  # initial guess of a and b1

    sigma = [v] * len(t)

    pars, cov = curve_fit(Pmodel, to, Po, p0, sigma=sigma)

    return pars, cov


def best_fit_subsidence(v):
    """ Calculate best fit parameters and covariance to subsidence data.

        Parameters
        ----------
        v : float
            Observational error of the subsidence data.

        Returns
        -------
        pars : tuple of float
            Tuple of the best fit parameters b2 and mv respectively.
        cov : array-like
            Covariance matrix.
    """
    # obtain subsidence data to fit model to
    to, Uo = load_subsidence_data()

    # interpolate subsidence at equal time intervals
    t = time_range(to[0], to[-1], 1)
    U = np.interp(t, to, Uo)

    p0 = [1.62e-3,7.15e-2]  # initial guess of b2 and mv

    sigma = [v] * len(t)  # create vector of observational error

    # find best fit parameters and covariance matrix
    pars, cov = curve_fit(Umodel, t, U, p0, sigma=sigma)

    return pars, cov


def model_ensemble():
    """ Plot the ensemble of models based on uncertainty of the data. """
    # obtain subsidence data
    to, Uo = load_subsidence_data()

    # calculate time range
    t = time_range(to[0], to[-1], 1)

    v_pres = 1.  # uncertainty on pressure observations
    p_pres, cov_pres = best_fit_pressure(v_pres)  # best fit parameters and covariance
    ps_pres = np.random.multivariate_normal(p_pres, cov_pres, 25)  # samples from posterior

    v_sub = 0.5  # uncertainty on subsidence observations
    p_sub, cov_sub = best_fit_subsidence(v_sub)  # best fit parameters and covariance
    ps_sub = np.random.multivariate_normal(p_sub, cov_sub, 25)  # samples from posterior

    # calculate mass extraction and derivative
    q = interpolate_mass_extraction_rate(t)
    dqdt = dq_dt(t, q)

    # get start and end times and time step
    t0 = t[0]
    t1 = t[-1]
    dt = t[1] - t[0]

    # fixed parameters
    c = 0.007883715689885993
    P0 = 56.26
    Pm0 = 56.26
    L = 100

    # create figure
    fig, (ax1, ax2) = plt.subplots(ncols=2)

    # loop through posterior samples
    for a, b1 in ps_pres:
        for b2, mv in ps_sub:
            # solve model for set of sampled parameters
            t, P = solve_reservoir_ode(reservoir_ode, t0, t1, dt, P0, q, dqdt, [a, b1, c, P0])
            t, Pm = solve_mudstone_ode(mudstone_ode, t0, t1, dt, P, Pm0, [b2])
            U = subsidence_eqn(Pm, Pm0, mv, L)

            # plot model for set of sampled parameters
            ax1.plot(t, P, 'k-', lw=0.25, alpha=0.2)
            ax2.plot(t, U, 'k-', lw=0.25, alpha=0.2)

    # labels
    ax1.plot([], [], 'k-', lw=0.5, label='posterior samples')
    ax2.plot([], [], 'k-', lw=0.5, label='posterior samples')

    # best fit parameters
    (a, b1) = p_pres
    (b2, mv) = p_sub

    # plot best fit model
    t, P = solve_reservoir_ode(reservoir_ode, t0, t1, dt, P0, q, dqdt, [a, b1, c, P0])
    t, Pm = solve_mudstone_ode(mudstone_ode, t0, t1, dt, P, Pm0, [b2])
    U = subsidence_eqn(Pm, Pm0, mv, L)
    ax1.plot(t, P, 'b-', label='best-fit')
    ax2.plot(t, U, 'b-', label='best-fit')

    # plot pressure observations with error bars
    ax1.errorbar(*load_pressure_data(), yerr=v_pres, fmt='ro', label='data')
    ax1.set_title('Reservoir pressure')
    ax1.set_xlabel('time [years]')
    ax1.set_ylabel('reservoir pressure [bars]')
    ax1.legend()

    # plot subsidence observations with error bars
    ax2.errorbar(*load_subsidence_data(), yerr=v_sub, fmt='ro', label='data')
    ax2.set_title('Subsidence')
    ax2.set_xlabel('time [years]')
    ax2.set_ylabel('subsidence [m]')
    ax2.legend()

    # display the figure
    fig.tight_layout()
    fig.set_size_inches(10,5)
    fig.show()

def forecast_ensemble():
    """ Forecast the reservoir pressure and subsidence with uncertainty. """
    # fixed parameters
    c = 0.007883715689885993  # lumped parameter for slow drainage term
    P0 = 56.26  # initial/ambient reservoir pressure
    L = 100  # thickness of mudstone
    Pm0 = 56.26  # initial/ambient mudstone pressure

    # obtain pressure observations
    to, Po = load_pressure_data()

    # obtain start and end times and time step
    t0 = to[0]
    t1 = to[-1]
    dt = 1

    # get observed time range plus prediction time range
    ti = time_range(t0, t1, dt)
    tP = time_range(t1+dt,2042,dt)
    t = np.concatenate((ti,tP))

    # get mass extraction rate and derivative
    q0 = interpolate_mass_extraction_rate(ti)

    v_pres = 1.  # uncertainty on pressure observations
    p_pres, cov_pres = best_fit_pressure(v_pres)  # best fit parameters and covariance
    ps_pres = np.random.multivariate_normal(p_pres, cov_pres, 25)  # samples from posterior

    v_sub = 0.5  # uncertainty on subsidence observations
    p_sub, cov_sub = best_fit_subsidence(v_sub)  # best fit parameters and covariance
    ps_sub = np.random.multivariate_normal(p_sub, cov_sub, 25)  # samples from posterior

    # create figure
    fig, (ax1, ax2) = plt.subplots(ncols=2)

    # forecast for q = 1250 kg/s
    qp = 1250*np.ones(len(tP))
    q=np.concatenate((q0,qp))
    dqdt = dq_dt(t, q)

    # loop through posterior samples
    for a, b1 in ps_pres:
        for b2, mv in ps_sub:
            # solve model for set of sampled parameters
            _, P = solve_reservoir_ode(reservoir_ode, t[0], t[-1], dt, P0, q, dqdt, [a, b1, c, P0])
            _, Pm = solve_mudstone_ode(mudstone_ode, t[0], t[-1], dt, P, Pm0, [b2])
            U = subsidence_eqn(Pm, Pm0, mv, L)

            # plot model for set of sampled parameters
            ax1.plot(t, P, 'b-', lw=0.25, alpha=0.05)
            ax2.plot(t, U, 'b-', lw=0.25, alpha=0.05)

    # labels
    ax1.plot([], [], 'b-', lw=0.5, label='q=1250kg/s')
    ax2.plot([], [], 'b-', lw=0.5, label='q=1250kg/s')

    # forecast for q = 900 kg/s
    qp = 900*np.ones(len(tP))
    q=np.concatenate((q0,qp))
    dqdt = dq_dt(t, q)

    # loop through posterior samples
    for a, b1 in ps_pres:
        for b2, mv in ps_sub:
            # solve model for set of sampled parameters
            _, P = solve_reservoir_ode(reservoir_ode, t[0], t[-1], dt, P0, q, dqdt, [a, b1, c, P0])
            _, Pm = solve_mudstone_ode(mudstone_ode, t[0], t[-1], dt, P, Pm0, [b2])
            U = subsidence_eqn(Pm, Pm0, mv, L)

            # plot model for set of sampled parameters
            ax1.plot(t, P, 'g-', lw=0.25, alpha=0.05)
            ax2.plot(t, U, 'g-', lw=0.25, alpha=0.05)

    # labels
    ax1.plot([], [], 'g-', lw=0.5, label='q=900kg/s')
    ax2.plot([], [], 'g-', lw=0.5, label='q=900kg/s')

    # forecast for q = 600 kg/s
    qp = 600*np.ones(len(tP))
    q=np.concatenate((q0,qp))
    dqdt = dq_dt(t, q)

    # loop through posterior samples
    for a, b1 in ps_pres:
        for b2, mv in ps_sub:
            # solve model for set of sampled parameters
            _, P = solve_reservoir_ode(reservoir_ode, t[0], t[-1], dt, P0, q, dqdt, [a, b1, c, P0])
            _, Pm = solve_mudstone_ode(mudstone_ode, t[0], t[-1], dt, P, Pm0, [b2])

            # plot model for set of sampled parameters
            U = subsidence_eqn(Pm, Pm0, mv, L)
            ax1.plot(t, P, 'r-', lw=0.25, alpha=0.05)
            ax2.plot(t, U, 'r-', lw=0.25, alpha=0.05)

    # labels
    ax1.plot([], [], 'r-', lw=0.5, label='q=600kg/s')
    ax2.plot([], [], 'r-', lw=0.5, label='q=600kg/s')

    # forecast for q = 0 kg/s
    qp = np.zeros(len(tP))
    q=np.concatenate((q0,qp))
    dqdt = dq_dt(t, q)

    # loop through posterior samples
    for a, b1 in ps_pres:
        for b2, mv in ps_sub:
            # solve model for set of sampled parameters
            _, P = solve_reservoir_ode(reservoir_ode, t[0], t[-1], dt, P0, q, dqdt, [a, b1, c, P0])
            _, Pm = solve_mudstone_ode(mudstone_ode, t[0], t[-1], dt, P, Pm0, [b2])
            U = subsidence_eqn(Pm, Pm0, mv, L)

            # plot model for set of sampled parameters
            ax1.plot(t, P, 'm-', lw=0.25, alpha=0.05)
            ax2.plot(t, U, 'm-', lw=0.25, alpha=0.05)

    # labels
    ax1.plot([], [], 'm-', lw=0.5, label='q=0kg/s')
    ax2.plot([], [], 'm-', lw=0.5, label='q=0kg/s')

    # best fit parameters
    a, b1 = p_pres
    b2, mv = p_sub

    # solve best fit model
    qi = interpolate_mass_extraction_rate(ti)
    _, P = solve_reservoir_ode(reservoir_ode, t0, t1, dt, P0, qi, dq_dt(ti, qi), [a, b1, c, P0])
    _, Pm = solve_mudstone_ode(mudstone_ode, t0, t1, dt, P, Pm0, [b2])
    U = subsidence_eqn(Pm, Pm0, mv, L)

    # plot best-fit solution
    ax1.plot(ti, P, label="best-fit model")
    ax2.plot(ti, U, label="best-fit model")

    # plot observations
    ax1.plot(*load_pressure_data(), '.', label="observations")
    ax2.plot(*load_subsidence_data(), '.', label="observations")

    # set title, axis labels and legend
    ax1.set_title('Forecast of pressure')
    ax2.set_title('Forecast of subsidence')
    ax1.set_ylabel('pressure [bar]')
    ax2.set_ylabel('subsidence [m]')
    ax1.set_xlabel('time [yrs]')
    ax2.set_xlabel('time [yrs]')
    ax1.legend()
    ax2.legend()

    # display figure
    fig.tight_layout()
    fig.set_size_inches(10,5)
    plt.show()


def parameter_histogram():
    """ Plot probability distribution histograms of uncertain parameters. """

    v_pres = 1.  # uncertainty on pressure observations
    p_pres, cov_pres = best_fit_pressure(v_pres)  # best fit parameters and covariance matrix
    ps_pres = np.random.multivariate_normal(p_pres, cov_pres, 1000)  # samples from posterior

    v_sub = 0.5  # uncertainty on subsidence observations
    p_sub, cov_sub = best_fit_subsidence(v_sub)  # best fit parameters and covariance matrix
    ps_sub = np.random.multivariate_normal(p_sub, cov_sub, 1000)  # samples from posterior

    # initialise parameter lists
    a = []
    b1 = []
    b2 = []
    mv = []

    # add samples to parameter lists
    for pi in ps_pres:
        a.append(pi[0])
        b1.append(pi[1])

    for pi in ps_sub:
        b2.append(pi[0])
        mv.append(pi[1])

    # plot distribution of parameter a
    fig, (ax1, ax2) = plt.subplots(ncols=2)
    ax1.hist(a,bins=20,density=True)
    ax1.axvline(np.percentile(a,5), color='r')
    ax1.axvline(np.percentile(a,95), color='r')
    ax1.set_title('Posterior distribution of a')
    ax1.set_xlabel('a')
    ax1.set_ylabel('probability density')

    # plot distribution of parameter b1
    ax2.hist(b1,bins=20,density=True)
    ax2.axvline(np.percentile(b1,5), color='r')
    ax2.axvline(np.percentile(b1,95), color='r')
    ax2.set_title('Posterior distribution of b1')
    ax2.set_xlabel('b1')
    ax2.set_ylabel('probability density')

    # display figure
    fig.tight_layout()
    fig.set_size_inches(10, 5)
    fig.show()

    fig, (ax1, ax2) = plt.subplots(ncols=2)
    # plot distribution of parameter b2
    ax1.hist(b2,bins=20,density=True)
    ax1.axvline(np.percentile(b2,5), color='r')
    ax1.axvline(np.percentile(b2,95), color='r')
    ax1.set_title('Posterior distribution of b2')
    ax1.set_xlabel('b2')
    ax1.set_ylabel('probability density')

    # plot distribution of parameter mv
    ax2.hist(mv,bins=20,density=True)
    ax2.axvline(np.percentile(mv,5), color='r')
    ax2.axvline(np.percentile(mv,95), color='r')
    ax2.set_title('Posterior distribution of mv')
    ax2.set_xlabel('mv')
    ax2.set_ylabel('probability density')

    # display figure
    fig.tight_layout()
    fig.set_size_inches(10, 5)
    fig.show()


if __name__ == "__main__":
    model_ensemble()
    forecast_ensemble()
    parameter_histogram()
from model_functions import *
import matplotlib.pyplot as plt


def solve_analytical(t0, t1, q0, P0, params):
    """ Solve for reservoir pressure, mudstone pressure and subsidence analytically.

        Parameters
        ----------
        t0 : float
            Start time (s).
        t1 : float
            End time (s).
        q0 : float
            Constant mass extraction rate (kg/s).
        P0 : float
            Ambient pressure (Pa).
        params: list
            List of parameters in the order [a, b1, b2, c, mv, L].

        Returns
        -------
        (t, P) : tuple of np.array
            Tuple of time vector (s) and analytical reservoir pressure solution vector (Pa).
        (t, Pm) : tuple of np.array
            Tuple of time vector (s) and analytical mudstone pressure solution vector (Pa).
        (t, P) : tuple of np.array
            Tuple of time vector (s) and analytical subsidence solution vector (m).

        Notes
        -----
        Assume constant mass extraction rate `q0`.
        Assume ambient pressure `P0` is equal between the reservoir and mudstone.
        Assume initial pressure is equal to the ambient pressure.
    """
    # unpack parameters
    [a, b1, b2, c, mv, L] = params
    
    # create range of times
    t = np.linspace(t0, t1)
    
    # solve reservoir pressure analytically
    P = P0 + a/b1*q0*(np.exp(-b1*t) - 1)

    # solve mudstone pressure analytically
    Pm = P0 + a/b1*q0*(b2/(b2-b1)*np.exp(-b1*t) - a*b1/(b2-b1)*np.exp(-b2*t) - 1)

    # calculate subsidence analytically
    U = mv*L*(P0 - Pm)
    
    return (t, P), (t, Pm), (t, U)


def solve_numerical(t0, t1, h, q0, P0, params):
    """ Solve for reservoir pressure, mudstone pressure and subsidence numerically.

        Parameters
        ----------
        t0 : float
            Start time (s).
        t1 : float
            End time (s).
        h: float
            Step size (s).
        q0 : float
            Constant mass extraction rate (kg/s).
        P0 : float
            Ambient pressure (Pa).
        params: list
            List of parameters in the order [a, b1, b2, c, mv, L].

        Returns
        -------
        (t, P) : tuple of np.array
            Tuple of time vector (s) and numerical reservoir pressure solution vector (Pa).
        (t, Pm) : tuple of np.array
            Tuple of time vector (s) and numerical mudstone pressure solution vector (Pa).
        (t, P) : tuple of np.array
            Tuple of time vector (s) and numerical subsidence solution vector (m).

        Notes
        -----
        Assume constant mass extraction rate `q0`.
        Assume ambient pressure `P0` is equal between the reservoir and mudstone.
        Assume initial pressure is equal to the ambient pressure.
    """
    # unpack parameters
    [a, b1, b2, c, mv, L] = params
    
    Pm0 = P0  # assume equal ambient pressure between mudstone and reservoir

    # create range of times at time step apart
    t = time_range(t0, t1, h)
    
    # calculate mass extraction and derivative
    q = q0*np.ones(len(t),)
    dqdt = np.divide(np.diff(q), np.diff(t))
    dqdt = np.append(dqdt, dqdt[-1])

    # solve reservoir pressure numerically
    _, P = solve_reservoir_ode(reservoir_ode, t0, t1, h, P0, q, dqdt, [a, b1, c, P0])

    # solve mudstone pressure numerically
    _, Pm = solve_mudstone_ode(mudstone_ode, t0, t1, h, P, Pm0, [b2])

    # calculate subsidence numerically
    U = subsidence_eqn(Pm, P0, mv, L)
    
    return (t, P), (t, Pm), (t, U)


def error_analysis(analytical, numerical):
    """  Calculate relative error between analytical and numerical solutions.

        Parameters
        ----------
        analytical : tuple of np.array
            Tuple containing the analytical time vector and analytical solution vector in order.
        numerical : tuple of np.array
            Tuple containing the numerical time vector and numerical solution vector in order.

        Returns
        -------
        t : np.array
            Numerical time vector (excluding bad points).
        error : np.array
            Relative error between analytical and numerical solutions (excluding bad points).

        Notes
        -----
        Size of np.array within each tuple should be the same.
        Error with denominator zero and corresponding time removed.
    """
    # find numerical times
    t = numerical[0]
    y_num = numerical[1]

    # interpolate analytical solution at numerical times
    y_interp = np.interp(numerical[0], *analytical)

    error = np.zeros(len(t))  # initialise error vector
    delete_indices = []  # initialise list of indices to delete

    # calculate relative error between numerical and absolute solution
    for i in range(len(t)):
        # check denominator is not zero
        if y_interp[i] == 0:
            delete_indices.append(i)
        else:
            error[i] = abs(y_interp[i] - y_num[i]) / y_interp[i]

    # remove points with zero denominator
    t = np.delete(t, delete_indices)
    error = np.delete(error, delete_indices)

    return t, error


def time_step_convergence(t0, t1, q0, P0, params):
    """ Solve numerical solutions on varying time steps to test for convergence.

        Parameters
        ----------
        t0 : float
            Start time.
        t1 : float
            End time.
        q0 : float
            Constant mass extraction rate.
        P0 : float
            Ambient pressure.
        params : list
            List of parameters in the order [a, b1, b2, c, mv, L].

        Returns
        -------
        (time_step_recip, P_final) : tuple of np.array
            Tuple containing reciprocals of time step and final reservoir pressure solutions.
        (time_step_recip, Pm_final) : tuple of np.array
            Tuple containing reciprocals of time step and final mudstone pressure solutions.
        (time_step_recip, U_final) : tuple of np.array
            Tuple containing reciprocals of time step and final subsidence solutions.

        Notes
        -----
        Assume constant mass extraction rate `q0`.
        Assume ambient pressure `P0` is equal between the reservoir and mudstone.
        Assume initial pressure is equal to the ambient pressure.
    """
    # unpack parameters
    [a, b1, b2, c, mv, L] = params

    Pm0 = P0  # assume equal ambient pressure between mudstone and reservoir

    # range of time steps
    time_steps = np.logspace(-2, 0, 20)

    # final values for reservoir pressure, mudstone pressure and subsidence
    final = {
        'P': [],
        'Pm': [],
        'U': []
    }

    # step through each time step
    for h in time_steps:
        # solve reservoir pressure
        t = time_range(t0, t1, h)

        # calculate mass extraction and derivative
        q = q0 * np.ones(len(t), )
        dqdt = np.divide(np.diff(q), np.diff(t))
        dqdt = np.append(dqdt, dqdt[-1])

        # solve reservoir pressure
        _, P = solve_reservoir_ode(reservoir_ode, t0, t1, h, P0, q, dqdt, [a, b1, c, P0])
        final['P'].append(P[-1])

        # solve mudstone pressure
        _, Pm = solve_mudstone_ode(mudstone_ode, t0, t1, h, P, Pm0, [b2])
        final['Pm'].append(Pm[-1])

        # calculate subsidence
        final['U'].append(subsidence_eqn(Pm[1], Pm0, mv, L))

    P_final = np.array(final['P'])
    Pm_final = np.array(final['Pm'])
    U_final = np.array(final['U'])

    time_steps_recip = 1 / time_steps

    return (time_steps_recip, P_final), (time_steps_recip, Pm_final), (time_steps_recip, U_final)


def benchmark():
    """ Perform benchmarking, error analysis and time convergence testing for reservoir pressure, mudstone pressure
        and subsidence models, and plot the results.

        Parameters
        ----------
        None

        Returns
        -------
        None
    """
    a = 1  # lumped parameter for forced extraction term
    b1 = 2  # lumped parameter for recharge term (into reservoir)
    b2 = 0.1  # lumped parameter for recharge term (into mudstone)
    c = 0  # lumped parameter for slow drainage term
    mv = 3  # compressibility of mudstone
    L = 2  # thickness of mudstone

    params = [a, b1, b2, c, mv, L]

    q0 = 3  # mass extraction rate (constant)
    P0 = 4  # initial/ambient reservoir pressure

    t0 = 0  # initial time
    t1 = 10  # final time
    h = 0.5  # time step

    # calculate analytical solutions
    P_analytical, Pm_analytical, U_analytical = solve_analytical(t0, t1, q0, P0, params)

    # calculate numerical solutions
    P_numerical, Pm_numerical, U_numerical = solve_numerical(t0, t1, h, q0, P0, params)

    # calculate relative errors
    P_error = error_analysis(P_analytical, P_numerical)
    Pm_error = error_analysis(Pm_analytical, Pm_numerical)
    U_error = error_analysis(U_analytical, U_numerical)

    # perform timestep convergence
    P_final, Pm_final, U_final = time_step_convergence(t0, t1, q0, P0, params)

    # plot benchmarking results for reservoir pressure
    plot_benchmark(P_analytical, P_numerical, P_error, P_final, "reservoir pressure [Pa]", "P(t=10)")

    # plot benchmarking results for mudstone pressure
    plot_benchmark(Pm_analytical, Pm_numerical, Pm_error, Pm_final, "mudstone pressure [Pa]", "Pmud(t=10)")

    # plot benchmarking results for subsidence
    plot_benchmark(U_analytical, U_numerical, U_error, U_final, "subsidence [m]", "U(t=10)")


def plot_benchmark(analytical, numerical, error, timestep, label, final_label):
    """ Plot benchmarking, error analysis and timestep convergence results.

        Parameters
        ----------
        analytical : tuple of np.array
            Tuple containing the analytical time vector and analytical solution vector in order.
        numerical : tuple of np.array
            Tuple containing the numerical time vector and numerical solution vector in order.
        error : tuple of np.array
            Tuple containing the numerical time vector and relative error vector in order.
        timestep : tuple of np.array
            Tuple containing the time step reciprocal vector and final solution vector in order.
        label : string
            y-axis label of value to be plotted.
        final_label : string
            y-axis label of final value to be plotted.

        Returns
        -------
        None
    """
    # create figure and row of axes
    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3)

    # plot numerical solution and analytical solution
    ax1.plot(*analytical, label="analytical solution")
    ax1.plot(*numerical, "x", label="numerical solution")
    ax1.legend()
    ax1.set_title("benchmark")
    ax1.set_xlabel("time [s]")
    ax1.set_ylabel(label)

    # plot error of numerical solution against analytical solution
    ax2.plot(*error, "x")
    ax2.set_title("error analysis")
    ax2.set_xlabel("time [s]")
    ax2.set_ylabel("error against benchmark")

    # plot final solution at t = 10 for varying timesteps
    ax3.plot(*timestep, "x")
    ax3.set_title("timestep convergence")
    ax3.set_xlabel("1/Δt")
    ax3.set_ylabel(final_label)
    ax3.set_xscale("log")

    # display plots
    fig.set_size_inches(15, 5)
    plt.show()


if __name__ == "__main__":
    benchmark()

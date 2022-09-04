from model_calibration import *
import matplotlib.pyplot as plt

oranges = plt.get_cmap('Oranges')

# TODO: select uncertain parameters (2-3)
# TODO: future predictions with uncertainty
# TODO: histogram of final value
# TODO: document code
# TODO: refactor code

def grid_search():
    a_best = 0.0015377621974240604
    b_best = 0.06624924440742241

    N = 51

    a = np.linspace(a_best*0.8, a_best*1.2, N)
    b = np.linspace(b_best*0.8, b_best*1.2, N)

    A, B = np.meshgrid(a, b, indexing='ij')

    S = np.zeros(A.shape)

    to, Po = load_pressure_data()

    t0 = to[0]
    t1 = to[-1]
    dt = to[1] - to[0]

    t = time_range(t0, t1, dt)
    q = interpolate_mass_extraction_rate(t)
    dqdt = dq_dt(t, q)

    P0 = 56.26
    c = 0.007883715689885993

    v = 2.  # error variance - 2 bar

    for i in range(len(a)):
        for j in range(len(b)):
            _, P = solve_reservoir_ode(reservoir_ode, t0, t1, dt, P0, q, dqdt, [a[i], b[j], c, P0])
            S[i,j] = np.sum((Po - P)**2)/v

    P = np.exp(-S/2.)

    Pint = np.sum(P) * (a[1] - a[0]) * (b[1] - b[0])

    P = P/Pint

    plot_posterior(a, b, P=P)

    return a, b, P

def plot_posterior(a, b, c=None, P=None):
    if c is None:
        plot_posterior2D(a, b, P)
    else:
        plot_poseterior3D(a, b, c, P)

def plot_posterior2D(a, b, P):
    A, B = np.meshgrid(a, b)

    fig = plt.figure(figsize=[10., 7.])
    ax1 = fig.add_subplot(111, projection='3d')
    ax1.plot_surface(A, B, P, rstride=1, cstride=1, cmap=oranges, lw=0.5, edgecolor='k')

    ax1.set_xlabel('a')
    ax1.set_ylabel('b')
    ax1.set_xlabel('P')

    ax1.set_xlim(a[0], a[-1])
    ax1.set_ylim(b[0], b[-1])
    ax1.set_zlim(0., )

    ax1.view_init(40, 100.)

    plt.show()

def plot_poseterior3D(a, b, c, P):
    azim = 15.

    Ab, Ba = np.meshgrid(a, b, indexing='ij')
    Pab = np.zeros(Ab.shape)
    for i in range(len(a)):
        for j in range(len(b)):
            Pab[i][j] = sum([P[i][j][k] for k in range(len(c))])

    Ac, Ca = np.meshgrid(a, c, indexing='ij')
    Pac = np.zeros(Ac.shape)
    for i in range(len(a)):
        for k in range(len(c)):
            Pac[i][k] = sum([P[i][j][k] for j in range(len(b))])

    Bc, Cb = np.meshgrid(b, c, indexing='ij')
    Pbc = np.zeros(Bc.shape)
    for j in range(len(b)):
        for k in range(len(c)):
            Pbc[j][k] = sum([P[i][j][k] for i in range(len(a))])

    fig = plt.figure(figsize=[20.0, 15.])
    ax1 = fig.subplot(221, projection='3d')
    ax1.plot_surface(Ab, Ba, Pab, rstride=1, cstride=1, cmap=oranges, lw=0.5)

    ax1.set_xlabel('a')
    ax1.set_ylabel('b')
    ax1.set_zlabel('P')
    ax1.set_xlim([a[0], a[-1]])
    ax1.set_ylim([b[0], b[-1]])
    ax1.set_zlim(0., )
    ax1.view_init(40, azim)

    ax2 = fig.add_subplot(222, projection='3d')
    ax2.plot_surface(Ac, Ca, Pac, rstride=1, cstride=1, cmap=oranges, lw=0.5)
    ax2.set_xlabel('a')
    ax2.set_ylabel('c')
    ax2.set_zlabel('P')
    ax2.set_xlim([a[0], a[-1]])
    ax2.set_ylim([c[0], c[-1]])
    ax2.set_zlim(0., )
    ax2.view_init(40, azim)

    ax3 = fig.add_subplot(223, projection='3d')
    ax3.plot_surface(Bc, Cb, Pbc, rstride=1, cstride=1, cmap=oranges, lw=0.5)
    ax3.set_xlabel('b')
    ax3.set_ylabel('c')
    ax3.set_zlabel('P')
    ax3.set_xlim([b[0], b[-1]])
    ax3.set_ylim([c[0], c[-1]])
    ax3.set_zlim(0., )
    ax3.view_init(40, azim)

    # save and show
    plt.show()

def fit_mvn(parspace, dist):
    N = len(parspace)

    parspace = [p.flatten() for p in parspace]
    dist = dist.flatten()

    mean = [np.sum(dist*par)/np.sum(dist) for par in parspace]

    cov = np.zeros((N,N))
    for i in range(N):
        for j in range(i,N):
            cov[i,j] = np.sum(dist*(parspace[i] - mean[i]) * (parspace[j] - mean[j])) / np.sum(dist)
            if i != j: cov[j,i] = cov[i,j]

    return np.array(mean), np.array(cov)

def construct_samples(a, b, P, N_samples):
    A, B = np.meshgrid(a, b, indexing='ij')
    mean, covariance = fit_mvn([A, B], P)

    samples = np.random.multivariate_normal(mean, covariance, N_samples)

    plot_samples(a, b, P=P, samples=samples)

    return samples

def plot_samples(a, b, c=None, P=None, samples=None):
    if c is None:
        plot_samples2D(a, b, P, samples)
    else:
        plot_samples3D(a, b, c, P, samples)

def plot_samples2D(a, b, P, samples):
    fig = plt.figure(figsize=[10., 7.])
    ax1 = fig.add_subplot(111, projection='3d')
    A, B = np.meshgrid(a, b, indexing='ij')
    ax1.plot_surface(A, B, P, rstride=1, cstride=1, cmap=oranges, lw=0.5)

    to, Po = load_pressure_data()
    v = 2

    t0 = to[0]
    t1 = to[-1]
    dt = to[1] - to[0]

    t = time_range(t0, t1, dt)
    q = interpolate_mass_extraction_rate(t)
    dqdt = dq_dt(t, q)

    P0 = 56.26
    c = 0.007883715689885993

    s = []
    for a, b in samples:
        _, pres = solve_reservoir_ode(reservoir_ode, t0, t1, dt, P0, q, dqdt, [a, b, c, P0])
        s.append(np.sum((pres - Po)**2) / v)
    s = np.array(s)

    p = np.exp(-s/2.)
    p = p/np.max(p) * np.max(P)

    ax1.plot(*samples.T, p, 'k.')

    ax1.set_xlabel('a')
    ax1.set_ylabel('b')
    ax1.set_zlabel('P')
    ax1.set_zlim(0., )
    ax1.view_init(40, 100.)

    # save and show
    plt.show()

def plot_samples3D(a, b, c, P, samples):
    azim = 15.

    # a and b combination
    Ab, Ba = np.meshgrid(a, b, indexing='ij')
    Pab = np.zeros(Ab.shape)
    for i in range(len(a)):
        for j in range(len(b)):
            Pab[i][j] = sum([P[i][j][k] for k in range(len(c))])

    # a and c combination
    Ac, Ca = np.meshgrid(a, c, indexing='ij')
    Pac = np.zeros(Ac.shape)
    for i in range(len(a)):
        for k in range(len(c)):
            Pac[i][k] = sum([P[i][j][k] for j in range(len(b))])

    # b and c combination
    Bc, Cb = np.meshgrid(b, c, indexing='ij')
    Pbc = np.zeros(Bc.shape)
    for j in range(len(b)):
        for k in range(len(c)):
            Pbc[j][k] = sum([P[i][j][k] for i in range(len(a))])

    to, Po = load_pressure_data()
    v = 2

    t0 = to[0]
    t1 = to[-1]
    dt = to[1] - to[0]

    t = time_range(t0, t1, dt)
    q = interpolate_mass_extraction_rate(t)
    dqdt = dq_dt(t, q)

    P0 = 56.26

    s = []
    for a, b, c in samples:
        _, pres = solve_reservoir_ode(reservoir_ode, t0, t1, dt, P0, q, dqdt, [a, b, c, P0])
        s.append(np.sum((pres - Po)**2) / v)
    s = np.array(s)

    p = np.exp(-s / 2.)
    p = p / np.max(p) * np.max(P) * 1.2

    # plotting
    fig = plt.figure(figsize=[20.0, 15.])
    ax1 = fig.add_subplot(221, projection='3d')
    ax1.plot_surface(Ab, Ba, Pab, rstride=1, cstride=1, cmap=oranges, lw=0.5)
    ax1.set_xlabel('a')
    ax1.set_ylabel('b')
    ax1.set_zlabel('P')
    ax1.set_xlim([a[0], a[-1]])
    ax1.set_ylim([b[0], b[-1]])
    ax1.set_zlim(0., )
    ax1.view_init(40, azim)
    ax1.plot(samples[:, 0], samples[:, 1], p, 'k.')

    ax1 = fig.add_subplot(222, projection='3d')
    ax1.plot_surface(Ac, Ca, Pac, rstride=1, cstride=1, cmap=oranges, lw=0.5)
    ax1.set_xlabel('a')
    ax1.set_ylabel('c')
    ax1.set_zlabel('P')
    ax1.set_xlim([a[0], a[-1]])
    ax1.set_ylim([c[0], c[-1]])
    ax1.set_zlim(0., )
    ax1.view_init(40, azim)
    ax1.plot(samples[:, 0], samples[:, -1], p, 'k.')

    ax1 = fig.add_subplot(223, projection='3d')
    ax1.plot_surface(Bc, Cb, Pbc, rstride=1, cstride=1, cmap=oranges, lw=0.5)
    ax1.set_xlabel('b')
    ax1.set_ylabel('c')
    ax1.set_zlabel('P')
    ax1.set_xlim([b[0], b[-1]])
    ax1.set_ylim([c[0], c[-1]])
    ax1.set_zlim(0., )
    ax1.view_init(40, azim)
    ax1.plot(samples[:, 1], samples[:, -1], p, 'k.')

    # save and show
    plt.show()

def model_ensemble(samples):
    to, Po = load_pressure_data()

    f, ax = plt.subplots(1, 1)

    t0 = to[0]
    t1 = to[-1]
    dt = to[1] - to[0]

    t = time_range(t0, t1, dt)
    q = interpolate_mass_extraction_rate(t)
    dqdt = dq_dt(t, q)

    P0 = 56.26
    c = 0.007883715689885993

    for a, b in samples:
        _, pres = solve_reservoir_ode(reservoir_ode, t0, t1, dt, P0, q, dqdt, [a, b, c, P0])
        ax.plot(t, pres, 'k-', lw=0.25, alpha=0.2)
    ax.plot([], [], 'k-', lw=0.5, alpha=0.4, label='model ensemble')

    v = 2.
    ax.errorbar(to, Po, yerr=v, fmt='ro', label='data')
    ax.set_xlabel('time [years]')
    ax.set_ylabel('pressure [bars]')
    ax.legend()
    plt.show()

if __name__ == "__main__":
    a, b, posterior = grid_search()

    N = 100
    samples = construct_samples(a, b, posterior, N)

    model_ensemble(samples)

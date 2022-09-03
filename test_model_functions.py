from model_functions import *


def test_reservoir_ode():
    got = reservoir_ode(0, 1, 2, 3, 4, 5, 6, 7)
    exact = 4

    assert got == exact


def test_mudstone_ode():
    got = mudstone_ode(0, 1, 2, 3)
    exact = 3

    assert got == exact


def test_subsidence_eqn():
    got = subsidence_eqn(0, 1, 2, 3)
    exact = 6

    assert got == exact


def test_solve_reservoir_ode():
    # 1 step
    q = np.array([3, 3])
    dqdt = np.array([0, 0])
    _, soln = solve_reservoir_ode(reservoir_ode, 0, 0.1, 0.1, 5, q, dqdt, [1, 2, 0, 4])
    got = soln[-1]
    exact = 4.55

    assert got == exact

    # 2 steps
    q = np.array([3, 3, 3])
    dqdt = np.array([0, 0, 0])
    _, soln = solve_reservoir_ode(reservoir_ode, 0, 0.2, 0.1, 5, q, dqdt, [1, 2, 0, 4])
    got = soln[-1]
    exact = 4.181

    assert got == exact


def test_solve_mudstone_ode():
    # 1 step
    P = np.array([4, 4])
    _, soln = solve_mudstone_ode(mudstone_ode, 0, 0.1, 0.1, P, 5, [0.5])
    got = soln[-1]
    exact = 4.95125

    assert got == exact

    # 2 steps
    P = np.array([4, 4, 4])
    _, soln = solve_mudstone_ode(mudstone_ode, 0, 0.2, 0.1, P, 5, [0.5])
    got = soln[-1]
    exact = 4.9048765625

    assert got == exact

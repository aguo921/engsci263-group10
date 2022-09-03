from benchmark import *
from curve_fitting import *
from prediction import *

if __name__ == "__main":
    benchmark()

    calibrate_reservoir_pressure()
    calibrate_subsidence()

    plot_model()

    plot_misfit(True)
    plot_misfit(False)
    
    forecast([1250, 900, 600, 0])

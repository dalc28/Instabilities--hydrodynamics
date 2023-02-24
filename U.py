import numpy as np
def U(x):
    vel = (0.5)*(1+np.tanh(x))  # Velocity profile
    return vel
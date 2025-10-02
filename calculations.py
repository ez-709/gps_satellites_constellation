import numpy as np

def calculate_samples_from_hours(end_time_hours, step_seconds=120): 
    total_time_seconds = end_time_hours * 3600
    samples = int(total_time_seconds / step_seconds)
    return samples, step_seconds

def mean_anomaly(M_0, a, dt_seconds):
    mu = 13.986004418e14
    n = np.sqrt(mu / a**3)
    M = M_0 + dt_seconds * n
    return M % (2 * np.pi) 

def kepler_equation(M, e, method_steps=5):
    E = M
    for i in range(method_steps):
        numerator = E - e * np.sin(E) - M
        denominator = 1 - e * np.cos(E)
        E = E - numerator / denominator
    return E

def true_anomaly(e, E):
    tan_half = np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2)
    v = 2 * np.arctan(tan_half)
    return v

def calculate_r(a, e, E):
    return a * (1 - e * np.cos(E))

def calculate_r_ISK(r, v, i, w, omega):
    x = r * (np.cos(v) * np.cos(omega) - np.sin(v) * np.cos(i) * np.sin(omega))
    y = r * (np.cos(v) * np.sin(omega) - np.sin(v) * np.cos(i) * np.cos(omega))
    z = r * (np.sin(v) * np.sin(i))
    return np.array([[x], [y], [z]])
        
def calculate_ISK_coordinates(a, e, M_0, i, w, omega, time_agp_start, end_time_hours, step_seconds=600):
    samples, step_seconds = calculate_samples_from_hours(end_time_hours, step_seconds)
    total_time_seconds = end_time_hours * 3600
    times = np.linspace(time_agp_start, time_agp_start + total_time_seconds, samples) 
    
    r_isk_vectors = []
    for time in times:
        dt_seconds = time - time_agp_start 
        M = mean_anomaly(M_0, a, dt_seconds) 
        E = kepler_equation(M, e) 
        v = true_anomaly(e, E)  
        r = calculate_r(a, e, E)
        r_isk = calculate_r_ISK(r, v, i, w, omega)
        r_isk_vectors.append({'r_isk': r_isk, 'time_point': time})
    return r_isk_vectors

def from_ECI_to_DGCS(r_ECI):
    matrix = 
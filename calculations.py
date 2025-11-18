import numpy as np

def matrix_of_direction_cos(psi=0, gamma=0, theta=0):
    matrix = np.eye(3) 

    Rx = np.array([
        [1, 0, 0],
        [0, np.cos(psi), -np.sin(psi)],
        [0, np.sin(psi), np.cos(psi)]
    ])
    matrix = np.dot(matrix, Rx)

    Ry = np.array([
        [np.cos(gamma), 0, np.sin(gamma)],
        [0, 1, 0],
        [-np.sin(gamma), 0, np.cos(gamma)]
    ])
    matrix = np.dot(matrix, Ry)

    Rz = np.array([
        [np.cos(theta), -np.sin(theta), 0],
        [np.sin(theta), np.cos(theta), 0],
        [0, 0, 1]
    ])
    matrix = np.dot(matrix, Rz)
    
    return matrix


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
    sin_v = np.sqrt(1 - e**2) * np.sin(E) / (1 - e * np.cos(E))
    cos_v = (np.cos(E) - e) / (1 - e * np.cos(E))
    return np.arctan2(sin_v, cos_v)

def calculate_r(a, e, E):
    return a * (1 - e * np.cos(E))

def calculate_r_ECI(r, v, w, i, Omega):
    u = v + w
    x = r * (np.cos(u) * np.cos(Omega) - np.sin(u) * np.sin(Omega) * np.cos(i))
    y = r * (np.cos(u) * np.sin(Omega) + np.sin(u) * np.cos(Omega) * np.cos(i))
    z = r * (np.sin(u) * np.sin(i))
    return np.array([[x], [y], [z]])
        
def calculate_ECI_coordinates(a, e, M_0,w, i, omega, time_agp_start, end_time_hours, step_seconds=600):
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
        r_isk = calculate_r_ECI(r, v, w, i, omega)
        r_isk_vectors.append({'r_ECI': r_isk, 'time_point': time})
    return r_isk_vectors

def from_ECI_to_DGCS(r_ECI_dicts):
    DGSK_coors = []
    u = 7.2921158554e-5
    for dict in r_ECI_dicts:
        r_ECI = dict['r_ECI']
        time = dict['time_point']
        angle = u * time
        A = np.array([
            [np.cos(angle),  np.sin(angle), 0],
            [-np.sin(angle), np.cos(angle), 0],
            [0, 0, 1]
        ])
    
        r_DGSK = A @ r_ECI
        DGSK_coors.append(r_DGSK)
    return DGSK_coors

def geodetic_to_DGCS(lat, lon, alt):
    a = 6378137
    f = 1 / 298.257223563  
    e2 = 2 * f - f**2       

    N = a / np.sqrt(1 - e2 * np.sin(lat)**2)

    x = (N + alt) * np.cos(lat) * np.cos(lon)
    y = (N + alt) * np.cos(lat) * np.sin(lon)
    z = ((1 - e2) * N + alt) * np.sin(lat)

    return np.array([[x], [y], [z]])

def calculate_EL_AZ(DGSK_dicts, r_o_dgsk):
    x_o, y_o, z_o = float(r_o_dgsk[0]), float(r_o_dgsk[1]), float(r_o_dgsk[2])
    R = np.sqrt(x_o**2 + y_o**2 + z_o**2)
    R_g = np.sqrt(x_o**2 + y_o**2)
    
    res_dict = []
    for dict in DGSK_dicts:
        r_s = dict['r_s_dgsk']
        time = dict['time_point']
        x_s, y_s, z_s = float(r_s[0]), float(r_s[1]), float(r_s[2])
        r_dif = np.array([x_s - x_o, y_s - y_o, z_s - z_o])

        if R_g == 0:
            res_dict.append({'el': None, 
                             'az' : +None, 
                             'time' : time})
            continue

        matrix = np.array([
            [-y_o/R_g, x_o/R_g, 0.0],
            [-x_o*z_o/(R*R_g), -y_o*z_o/(R*R_g), R_g/R],
            [ x_o/R, y_o/R, z_o/R]
        ])
        
        r_x_sr = matrix @ r_dif
        abs_r = np.linalg.norm(r_x_sr)

        sin_el = r_x_sr[2] / abs_r
        sin_el = np.clip(sin_el, -1.0, 1.0)
        el = np.arcsin(sin_el)
        az = np.arctan2(r_x_sr[0], r_x_sr[1])
        if az < 0:
            az += 2 * np.pi
        res_dict.append({'el' : el, 'az':az, 'time':time})

    return res_dict

def count_sats_in_the_sky(all_sats_angels, r_dgsk_dict, gamma_mask):
    first_sat = next(iter(all_sats_angels), None)

    num_time_steps = len(all_sats_angels[first_sat])
    result = []

    for i in range(num_time_steps):
        time = all_sats_angels[first_sat][i]['time']
        count = 0

        for sat_name in all_sats_angels:
            el = all_sats_angels[sat_name][i]['el']
            r_dgsk_entry = r_dgsk_dict[sat_name][i]

            if el >= np.radians(gamma_mask):
                count += 1
            else:
                r_dgsk_entry['r_s_dgsk'] = None

        result.append({'time': time, 'visible_count': count})

    return result, r_dgsk_dict

'''
def calculate_ionospheric_delay(el, az, phi_o, lamda_o, alpha, beta, true_distance):
    if el <= 0:
        return 0.0

    Az = az / np.pi
    El = el / np.pi
    phi_o = phi_o / np.pi
    lamda_o = lamda_o / np.pi

    phi = 0.0137 * El + 0.11

    phi_i = phi_o + phi * np.cos(Az * np.pi)

    if phi_i > 0.416:
        phi_i = 0.416
    elif phi_i < -0.416:
        phi_i = -0.416

    lamda_i = lamda_o + np.sin(az * np.pi)/np.cos(phi_i * np.pi)
    phi_m = phi_i + 0.064 * np.cos((lamda_i - 1.617) * np.pi)
    T_em = true_distance / 299792458.0 
    T = 4.32 * 1e4 * lamda_i + T_em

    F = 1 + 16 * (0.53 - el)^3
'''


def psevdo_r(r_dgsk_dict, r_o_dgsk, noise_level=2.0):
    r_o = np.array(r_o_dgsk)
    psevdo_r_dict = {}

    for sat_name, sat_list in r_dgsk_dict.items():
        pseudorange_list = []
        for entry in sat_list:
            r_s_raw = entry['r_s_dgsk']
            time_point = entry['time_point']

            if r_s_raw is None:
                pseudorange = None
            else:
                r_s = np.array(r_s_raw)
                true_range = np.linalg.norm(r_s - r_o)
                noise = np.random.uniform(-noise_level, noise_level)
                pseudorange = true_range + noise

            pseudorange_list.append({
                'psevdo r': pseudorange,
                'time_point': time_point
            })

        psevdo_r_dict[sat_name] = pseudorange_list

    return psevdo_r_dict

def calculate_possition(time_point_data):
    if len(time_point_data) < 4:
        return None
    
    satellite_list = list(time_point_data.values())

    pos = np.array([0.0, 0.0, 0.0])

    for _ in range(100):
        residuals = []
        H = []

        for sat_data in satellite_list:
            pseudorange = sat_data['pseudo_r']
            sat_coords_dgsk = np.array(sat_data['r_dgsk']).flatten()

            dist_vec = sat_coords_dgsk - pos
            dist = np.linalg.norm(dist_vec)

            if dist < 1e-8:
                continue

            residuals.append(pseudorange - dist)

            ex, ey, ez = dist_vec / dist

            H.append([ex, ey, ez, 1])

        H = np.array(H)
        residuals = np.array(residuals)

        if len(residuals) < 4:
            return None  

        try:
            a1 = H.T @ H
            a2 = np.linalg.inv(a1)
            a3 = a2 @ H.T
            delta = a3 @ residuals

            pos = pos - delta[:3]

        except np.linalg.LinAlgError:
            return None

    return pos.tolist()

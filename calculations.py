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


def psevdo_r(r_dgsk_dict, r_o_dgsk, noise_level=3.0):
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

def least_squares_position(visible_sat_coords, pseudoranges, max_iter=5, max_error=0.01):
    r_est = np.array([[0.0], [0.0], [0.0]])

    N = len(visible_sat_coords)
    if N < 4:
        raise ValueError("Для решения МНК необходимо минимум 4 видимых НКА.")

    for _ in range(max_iter):
        H = []  
        z = [] 

        for i in range(N):
            r_sat = visible_sat_coords[i]  
            rho_meas = pseudoranges[i]  

            rho_calc = np.linalg.norm(r_sat - r_est)

            z_i = rho_meas - rho_calc
            z.append(z_i)

            e_x = (r_sat[0, 0] - r_est[0, 0]) / rho_calc
            e_y = (r_sat[1, 0] - r_est[1, 0]) / rho_calc
            e_z = (r_sat[2, 0] - r_est[2, 0]) / rho_calc

            H.append([e_x, e_y, e_z, 1.0])

        H = np.array(H)            
        z = np.array(z).reshape(-1, 1) 

        x_ = np.linalg.inv(H.T @ H)@ (H.T @ z)  

        delta_r = x_[:3].reshape(3, 1)
        r_new = r_est + delta_r

        if np.linalg.norm(delta_r) < max_error:
            r_est = r_new
            break

        r_est = r_new

    return r_est

def calculate_possiotion_errors_for_end_time_hours(all_sats_DGSK, psevdo_r_dict, r_o_DGSK, time_agp, step, end_hours = 3):
    start_time = time_agp
    end_time = start_time + end_hours * 3600  
    step_lr3 = 60                      
    num_steps = int((end_time - start_time) / step_lr3)

    errors_X = []
    errors_Y = []
    errors_Z = []
    time_hours = []

    for k in range(num_steps):
        t_current = start_time + k * step_lr3

        visible_sats = []
        pseudoranges = []

        for sat_name in all_sats_DGSK:
            sat_list = all_sats_DGSK[sat_name]
            pr_list = psevdo_r_dict[sat_name]

            for i, entry in enumerate(sat_list):
                if abs(entry['time_point'] - t_current) < step/ 2:
                    r_s = entry['r_s_dgsk']
                    rho = pr_list[i]['psevdo r']
                    if r_s is not None and rho is not None:
                        visible_sats.append(r_s)
                        pseudoranges.append(rho)
                    break

        if len(visible_sats) < 4:
            continue

        try:
            r_mnk = least_squares_position(visible_sats, pseudoranges)
        except ValueError:
            continue

        err_x = abs(r_o_DGSK[0, 0] - r_mnk[0, 0])
        err_y = abs(r_o_DGSK[1, 0] - r_mnk[1, 0])
        err_z = abs(r_o_DGSK[2, 0] - r_mnk[2, 0])

        errors_X.append(err_x)
        errors_Y.append(err_y)
        errors_Z.append(err_z)
        time_hours.append(k * step_lr3 / 3600.0)

    return errors_X, errors_Y, errors_Z, time_hours
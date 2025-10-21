import os
import matplotlib.pyplot as plt
import numpy as np

from utils import read_config
from calculations import calculate_ECI_coordinates, from_ECI_to_DGCS, count_sats_in_the_sky, geodetic_to_DGCS, calculate_EL_AZ, psevdo_r
from visualization import vizualize_orbits_3d, plot_visible_satellites, plot_psevdo_r
from parser import parse_orbit_file

cd = os.getcwd()
cd_agp = os.path.join(cd, 'data', '2025', 'mcct_250908.agp')

if __name__ == "__main__":
    config_path = os.path.join(os.path.dirname(__file__), 'config.json')
    observer_longitude, observer_latitude, observer_altitude, step_seconds, gamma_max, noise_level = read_config(config_path)

    sats_coords = parse_orbit_file(cd_agp)
    all_orbits = []
    all_sats_ECI = []
    
    for idx, sat in enumerate(sats_coords):        
        orbit_data = calculate_ECI_coordinates(
            a=sat['a'],  
            e=sat['e'],
            M_0=sat['M_0'],
            i=sat['i'],
            w = sat['w'],
            omega=sat['omega'],
            time_agp_start=sat['unix_time'],
            end_time_hours=12,
            step_seconds=step_seconds
        )

        time_agp = sat['unix_time']
        x_coords = []
        y_coords = []
        z_coords = []
        
        for point in orbit_data:
            r_isk = point['r_ECI']
            x_coords.append(float(r_isk[0][0]) / 1000)
            y_coords.append(float(r_isk[1][0]) / 1000)
            z_coords.append(float(r_isk[2][0]) / 1000)
        
        all_sats_ECI.append(orbit_data)
            
        all_orbits.append({
            'x': x_coords,
            'y': y_coords,
            'z': z_coords
        })
    
    fig1, ax1 = vizualize_orbits_3d(all_orbits, time_agp, observer_altitude = observer_altitude, 
                     observer_latitude=observer_latitude, observer_longitude=observer_longitude)
    
    r_o_dgsk = geodetic_to_DGCS(np.radians(observer_latitude), np.radians(observer_longitude), observer_altitude)
    all_sats_DGSK = {}
    all_sats_angels = {}
    for num, ECI_data in enumerate(all_sats_ECI):
        dgsk_coords = from_ECI_to_DGCS(ECI_data) 
        sat_dgsk = []
        for i in range(len(dgsk_coords)):
            sat_dgsk.append({
                'r_s_dgsk': dgsk_coords[i],
                'time_point': ECI_data[i]['time_point']
            })

        all_sats_DGSK[f'sat {num}'] = sat_dgsk
        all_sats_angels[f'sat {num}'] = calculate_EL_AZ(sat_dgsk, r_o_dgsk) 
        
    sats_in_the_sky_by_time, all_sats_DGSK_vis = count_sats_in_the_sky(all_sats_angels, all_sats_DGSK, gamma_max)
    psevdo_r_dict = psevdo_r(all_sats_DGSK_vis, r_o_dgsk, noise_level)
    fig2, ax2 = plot_visible_satellites(sats_in_the_sky_by_time, sats_coords[0]['unix_time'])

    fig3, ax3 = plot_psevdo_r(psevdo_r_dict)

    plt.show() 

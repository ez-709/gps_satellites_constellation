import os

from utils import read_config
from calculations import calculate_ISK_coordinates
from visualization import vizualize_orbits
from parser import parse_orbit_file

cd = os.getcwd()
cd_agp = os.path.join(cd, 'data', '2025', 'mcct_250908.agp')

if __name__ == "__main__":
    config_path = os.path.join(os.path.dirname(__file__), 'config.json')
    observer_longitude, observer_latitude, observer_altitude, step_seconds = read_config(config_path)

    sats_coords = parse_orbit_file(cd_agp)
    all_orbits = []
    
    for idx, sat in enumerate(sats_coords):        
        orbit_data = calculate_ISK_coordinates(
            a=sat['a'],  
            e=sat['e'],
            M_0=sat['M_0'],
            omega=sat['omega'],
            i=sat['i'],
            Omega=sat['Omega'],
            time_agp_start=sat['unix_time'],
            end_time_hours=12,
            step_seconds=step_seconds
        )

        x_coords = []
        y_coords = []
        z_coords = []
        
        for point in orbit_data:
            r_isk = point['r_isk']
            x_coords.append(float(r_isk[0][0]) / 1000)
            y_coords.append(float(r_isk[1][0]) / 1000)
            z_coords.append(float(r_isk[2][0]) / 1000)
            
        all_orbits.append({
            'x': x_coords,
            'y': y_coords,
            'z': z_coords
        })
    
    vizualize_orbits(all_orbits, observer_altitude = observer_altitude, 
                     observer_latitude=observer_latitude, observer_longitude=observer_longitude)
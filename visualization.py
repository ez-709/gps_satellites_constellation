import numpy as np
import matplotlib.pyplot as plt
import datetime

from utils import unix_to_utc, proections_from_coordinates
from calculations import count_sats_in_the_sky

def vizualize_orbits_3d(orbit_data, time_agp, observer_longitude=None, observer_latitude=None, observer_altitude=None):
    
    colors = [
        "#D40A0A", "#E25C04", "#F1B501", "#C1B202", "#6E9E03",
        "#007034", "#055F54", "#0A4E75", "#2800BA", "#5A00B0",
        "#8C00A5", "#A23B72", "#896284", "#6F8996", "#00DCBF",
        "#00BFB3", "#00A2A6", "#008599", "#00688C", "#004B7F",
        "#002E72", "#000000", "#201203", "#402406", "#603609",
        "#803A0F", "#864A1C", "#8C5A29", "#926A36", "#987A43",
        "#9E8A50", "#A49A5D"
    ]
    
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d') 

    R_earth = 6371  
    
    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 25)
    x_earth = R_earth * np.outer(np.cos(u), np.sin(v))
    y_earth = R_earth * np.outer(np.sin(u), np.sin(v))
    z_earth = R_earth * np.outer(np.ones(np.size(u)), np.cos(v))
    
    ax.plot_surface(x_earth, y_earth, z_earth, 
                   color='lightblue', alpha=0.3, linewidth=0, antialiased=True)

    lon_g = 0 
    lat_g = np.linspace(-np.pi/2, np.pi/2, 100)
    x_g, y_g, z_g = proections_from_coordinates(lat_g, lon_g, 0)
    ax.plot(x_g, y_g, z_g, color='black', linewidth=1.5, linestyle='--', label='Prime Meridian')

    if observer_longitude is not None and observer_latitude is not None:
        lon_rad = np.radians(observer_longitude)
        lat_rad = np.radians(observer_latitude)
        alt = observer_altitude if observer_altitude is not None else 0
        x_obs, y_obs, z_obs = proections_from_coordinates(lat_rad, lon_rad, alt)
        ax.scatter([x_obs], [y_obs], [z_obs], color='red', s=50, label='Observer')

    for idx, orbit in enumerate(orbit_data):
        x = orbit['x']
        y = orbit['y']
        z = orbit['z']
        times = orbit.get('times', None)

        ax.plot(x, y, z, 
                color=colors[idx % len(colors)],
                linewidth=1.5)

        if len(x) > 0:
            ax.scatter([x[0]], [y[0]], [z[0]], color=colors[idx % len(colors)], s=30)
            if times is not None and len(times) > 0:
                utc_label = unix_to_utc(times[0])
                ax.text(x[0], y[0], z[0], f'  {utc_label}', fontsize=8)

    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')  
    ax.set_zlabel('Z (km)')
    ax.set_title(f'Satellite GPS Orbits at {datetime.datetime.utcfromtimestamp(time_agp)}')

    max_coord = R_earth * 1.1
    for orbit in orbit_data:
        orbit_max = max(max(abs(np.array(orbit['x']))), 
                       max(abs(np.array(orbit['y']))), 
                       max(abs(np.array(orbit['z']))))
        max_coord = max(max_coord, orbit_max * 1.1)
    
    ax.set_xlim([-max_coord, max_coord])
    ax.set_ylim([-max_coord, max_coord])
    ax.set_zlim([-max_coord, max_coord])
    
    plt.tight_layout()
    plt.legend()
    return fig, ax

def plot_visible_satellites(sats_in_the_sky_by_time, start_time_unix=None):
    times_unix = np.array([item['time'] for item in sats_in_the_sky_by_time])
    counts = [item['visible_count'] for item in sats_in_the_sky_by_time]

    if start_time_unix is None:
        start_time_unix = times_unix[0]

    times_hours = (times_unix - start_time_unix) / 3600.0

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(times_hours, counts)
    ax.set_xlabel('Time')
    ax.set_ylabel('Visivle sattellites')
    ax.set_title('Number of satellites in the sky by time')
    ax.grid(True)
    ax.set_xlim(left=0)
    ax.set_xticks(np.arange(0, max(times_hours) + 1, 1))

    return fig, ax

def plot_psevdo_r(psevdo_r_dict):
    sat_names = [f'sat {i}' for i in range(5)]
    colors = [
        "#D40A0A", "#E25C04", "#F1B501", "#C1B202", "#6E9E03",
        "#007034", "#055F54", "#0A4E75", "#2800BA", "#5A00B0",
        "#8C00A5", "#A23B72", "#896284", "#6F8996", "#00DCBF",
        "#00BFB3", "#00A2A6", "#008599", "#00688C", "#004B7F",
        "#002E72", "#000000", "#201203", "#402406", "#603609",
        "#803A0F", "#864A1C", "#8C5A29", "#926A36", "#987A43",
        "#9E8A50", "#A49A5D"
    ]
    color_dict = {sat_name: color for sat_name, color in zip(sat_names, colors)}

    start_time_unix = psevdo_r_dict['sat 0'][0]['time_point']

    fig, ax = plt.subplots(figsize=(12, 6))

    for sat_name in sat_names:
        data = psevdo_r_dict[sat_name]
        times = np.array([entry['time_point'] for entry in data])
        prs_raw = [entry['psevdo r'] for entry in data]
        
        segments = []
        current_segment_times = []
        current_segment_prs = []
        
        for i, pr in enumerate(prs_raw):
            if pr is not None:
                current_segment_times.append(times[i])
                current_segment_prs.append(pr)
            else:
                if current_segment_prs:
                    segments.append({
                        'times': current_segment_times,
                        'prs': current_segment_prs
                    })
                    current_segment_times = []
                    current_segment_prs = []

        if current_segment_prs:
            segments.append({
                'times': current_segment_times,
                'prs': current_segment_prs
            })
        
        color = color_dict[sat_name]
        for j, segment in enumerate(segments):
            times_rel = (np.array(segment['times']) - start_time_unix) / 3600.0
            prs_km = np.array(segment['prs']) / 1000.0

            label = sat_name if j == 0 else ""
            ax.plot(times_rel, prs_km, color=color, label=label, linewidth=1.5)

    ax.set_xlabel('Time [hours from start]')
    ax.set_ylabel('Pseudorange [km]')
    ax.set_title('Pseudoranges')
    ax.grid(True)
    ax.legend()
    ax.set_xlim(left=0)

    return fig, ax
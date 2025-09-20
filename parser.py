import numpy as np
from datetime import datetime, timezone

def parse_orbit_file(filename):    
    def line_to_floats(line):
        parts = line.strip().split()
        return [float(part.replace('E+', 'E').replace('E-', 'E-')) for part in parts]
    
    def date_to_unix(day, month, year, seconds_from_year_start):
        dt = datetime(year, month, day, tzinfo=timezone.utc)
        unix_time = dt.timestamp() + seconds_from_year_start
        return unix_time
    
    with open(filename, 'r') as f:
        lines = [line.rstrip() for line in f if line.strip()]
    
    result = []
    i = 0
    
    while i + 2 < len(lines):  
    
        line1_parts = lines[i].strip().split()
        day_recv = int(line1_parts[0])
        month_recv = int(line1_parts[1]) 
        year_recv = int(line1_parts[2])
        
        params2 = line_to_floats(lines[i + 1])
        
        params3 = line_to_floats(lines[i + 2])
        
        seconds_from_year_start = params2[7]
        health_flag = int(params2[8])   

        if health_flag != 0:
            i += 3
            continue

        Omega = params3[0]      
        i_orbit = params3[1]    
        omega = params3[2]     
        e = params3[3]         
        sqrt_a = params3[4]    
        M_0 = params3[5]      
        
        a = sqrt_a * sqrt_a
        
        Omega_rad = Omega * np.pi
        i_rad = i_orbit * np.pi
        omega_rad = omega * np.pi
        M_0_rad = M_0 * np.pi
        
        unix_time = date_to_unix(day_recv, month_recv, year_recv, seconds_from_year_start)
        
        data = {
            'unix_time': unix_time,
            'a': a,               
            'e': e,                
            'M_0': M_0_rad,         
            'omega': omega_rad,    
            'i': i_rad,         
            'Omega': Omega_rad   
        }
        
        result.append(data)
            
        i += 3
    return result
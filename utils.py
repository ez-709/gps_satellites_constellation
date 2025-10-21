from datetime import datetime, timezone, timedelta
import json
import numpy as np

def proections_from_coordinates(lat, lon, alt):
    R_e = 6371
    R_x = (R_e + alt) * np.cos(lat) * np.cos(lon)
    R_y = (R_e + alt) * np.cos(lat) * np.sin(lon)
    R_z = (R_e + alt) * np.sin(lat) 
    return R_x, R_y, R_z

def unix_to_utc(time_unix, time_zone = 3):
    time_str = str(datetime.fromtimestamp(time_unix, tz=timezone(timedelta(hours=time_zone))))
    datetime_part = time_str[:19]
    tz_part = time_str[26:29]
    tz_sign = tz_part[0]   
    tz_hours = tz_part[1:].lstrip('0') 

    utc_time = f"{datetime_part} {tz_sign}{tz_hours} UTC"
    return utc_time

def json_to_py(cd_json):
    '''Читает json с путем  'cd_json' и возвращает список словарей'''
    with open(cd_json, 'r', encoding='utf-8') as file:
        data = json.load(file)
    return data

def read_config(cd_config, observer_longitude=True, observer_latitude=True, 
                observer_altitude=True, step_seconds=True, gamma_max = True, noise_level = True):
    '''
    функция читает конфиг и возвращает список с нужными параметрами упорядоченными так же как и сам конфиг
    '''
    config = json_to_py(cd_config)
    out = []
    
    if observer_longitude == True:
        out.append(config.get('observer longitude'))
    
    if observer_latitude == True:
        out.append(config.get('observer latitude'))
    
    if observer_altitude == True:
        out.append(config.get('observer altitude'))
    if step_seconds == True:
        out.append(config.get('step seconds'))
    if gamma_max == True:
        out.append(gamma_max)
    if noise_level == True:
        out.append(config.get('noise level'))
    return out
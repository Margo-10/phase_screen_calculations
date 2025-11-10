import numpy as np
from pyproj import Transformer
import PyIRI
import PyIRI.main_library as ml
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import interp1d


transformer = Transformer.from_crs("EPSG:4326", "EPSG:4978", always_xy=True)
inv_transformer = Transformer.from_crs("EPSG:4978", "EPSG:4326", always_xy=True)

def geodetic_to_ecef(lat, lon, alt_m):
    return np.array(transformer.transform(lon, lat, alt_m))

def calculate_path_length(src_lat, src_lon, src_alt_m, rcv_lat, rcv_lon, rcv_alt_m):
    p1 = geodetic_to_ecef(src_lat, src_lon, src_alt_m)
    p2 = geodetic_to_ecef(rcv_lat, rcv_lon, rcv_alt_m)
    return np.linalg.norm(p2 - p1)

def sample_path_ecef(src_lat, src_lon, src_alt_m, rcv_lat, rcv_lon, rcv_alt_m, n=300):
    p1 = geodetic_to_ecef(src_lat, src_lon, src_alt_m)
    p2 = geodetic_to_ecef(rcv_lat, rcv_lon, rcv_alt_m)
    t = np.linspace(0, 1, n)
    return p1 + t[:, None] * (p2 - p1)

def ecef_to_geodetic(path_ecef):
    lon, lat, alt = inv_transformer.transform(path_ecef[:, 0], path_ecef[:, 1], path_ecef[:, 2])
    return lat, lon, alt

def get_ne_profile(path_ecef, year, month, day, hour, f107, ccir_or_ursi=0):
    lat, lon, alt = ecef_to_geodetic(path_ecef)
    alon = np.array(lon)
    alat = np.array(lat)
    aalt = np.array(alt / 1000)  # km
    ahr = np.full(len(lon), hour)
    
    _, _, _, _, _, _, edp = ml.IRI_density_1day(
        year, month, day, ahr, alon, alat, aalt, f107, PyIRI.coeff_dir, ccir_or_ursi
    )
    return edp[0, :, 0]  # [m^-3]

def compute_iri_path(src_lat, src_lon, src_alt_m, rcv_lat, rcv_lon, rcv_alt_m,
                     year, month, day, hour, num_screens, dN_N, f107, ccir_or_ursi=0):
    z_len = calculate_path_length(src_lat, src_lon, src_alt_m, rcv_lat, rcv_lon, rcv_alt_m)
    path_ecef = sample_path_ecef(src_lat, src_lon, src_alt_m, rcv_lat, rcv_lon, rcv_alt_m, n=300)
    ne_profile = get_ne_profile(path_ecef, year, month, day, hour, f107, ccir_or_ursi)
    z_profile = np.linspace(0, z_len, len(ne_profile))

   
    weights = ne_profile 
    cum_weight = cumulative_trapezoid(weights, z_profile, initial=0)
    total_weight = cum_weight[-1]
    sum_norm = cum_weight / total_weight
    z_func = interp1d(sum_norm, z_profile, kind='linear')

    u = np.linspace(0, 1, num_screens + 1)
    z_bounds = z_func(u)
    dz_list = np.diff(z_bounds)
  
    ne_local_list = []
    for i in range(num_screens):
        z1, z2 = z_bounds[i], z_bounds[i+1]
        mask = (z_profile >= z1) & (z_profile <= z2)
        ne_seg = ne_profile[mask]
        z_seg = z_profile[mask]
        ne_local = np.trapz(ne_seg, z_seg) / (z2 - z1)
        ne_local_list.append(ne_local)

    
    dne_list = [dN_N * ne for ne in ne_local_list]

    return z_len, dz_list, ne_local_list, dne_list
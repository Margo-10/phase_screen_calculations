import numpy as np
from scipy import special
from scipy.fft import fft, ifft, fftfreq, fftshift, ifftshift
import matplotlib.pyplot as plt

C_LIGHT = 299792458.0   # [m/s]
ELECTRON_RADIUS = 2.818e-15   # [m]

# === Default global parameters ===
PSD_TYPE = 2       # 0 - Gaussian; 1 - von Karman; 2 - Shkarofsky 
OUTER_SCALE = 5000.0  # [m]
INNER_SCALE = 10.0    # [m]
SPECTRAL_INDEX = 1.5  # spectral index


# === Grid initialization ===
def init_grid(lx: float, nx: int):
    N = nx
    LN = lx
    dx = LN / N
    x = dx * (np.arange(N) - N // 2)  
    kx = fftshift(fftfreq(N, d=dx) * 2.0 * np.pi)
    dkx = np.abs(kx[1] - kx[0])

    return {
        'N': N, 'LN': LN, 'dx': dx, 'x': x,
        'kx': kx, 'dkx': dkx
    }

# === PSD functions ===
def von_karman_dphi_1d_psd(dne, screen_width, radio_frequency, kx, outer_scale):
    k0 = 2.0 * np.pi / outer_scale
    gamma_ratio = special.gamma(SPECTRAL_INDEX / 2.0) / special.gamma((SPECTRAL_INDEX - 1) / 2.0)
    s_dne = (dne**2 / np.sqrt(np.pi)) * gamma_ratio * (k0**(SPECTRAL_INDEX - 1)) / ((kx**2 + k0**2)**(SPECTRAL_INDEX / 2.0))
    s_dphi = s_dne * (ELECTRON_RADIUS * C_LIGHT / radio_frequency)**2 * screen_width
    return s_dphi


def shkarovsky_dphi_1d_psd(dne, screen_width, radio_frequency, kx, outer_scale, inner_scale):
    k0 = 2.0 * np.pi / outer_scale
    ki = 2.0 * np.pi / inner_scale

    kv_ratio = special.kv(SPECTRAL_INDEX / 2.0, np.sqrt(kx**2 + k0**2) / ki) / special.kv((SPECTRAL_INDEX - 1.0) / 2.0, k0 / ki)
    denom = (np.sqrt(kx**2 + k0**2) / ki)**(SPECTRAL_INDEX / 2.0)

    s_dne = (dne**2 / np.sqrt(2 * np.pi) / ki) * (k0 / ki)**((SPECTRAL_INDEX - 1) / 2.0) * kv_ratio / denom

    arg1 = inner_scale / outer_scale
    l_eff = np.sqrt(2.0 * np.pi * outer_scale * inner_scale) * special.kv((SPECTRAL_INDEX - 1.0) / 2.0, arg1) / special.kv((SPECTRAL_INDEX - 2.0) / 2.0, arg1)

    s_dphi = s_dne * (ELECTRON_RADIUS * C_LIGHT / radio_frequency)**2 * screen_width * l_eff
    return s_dphi


def gauss_dphi_1d_psd(dne, screen_width, radio_frequency, kx, outer_scale):
    s_dne = dne**2 * outer_scale / 2.0 / np.sqrt(np.pi) * np.exp(-(kx * outer_scale / 2.0)**2)
    s_dphi = s_dne * (ELECTRON_RADIUS * C_LIGHT / radio_frequency)**2 * screen_width * outer_scale * np.sqrt(np.pi)
    return s_dphi

# === Phase screen generation ===
def generate_phase_screen(ne: float, screen_width: float, radio_frequency: float, grid: dict):
    kx = grid['kx']
    LN = grid['LN']
    N = grid['N']

   
    if PSD_TYPE == 0:
        f = gauss_dphi_1d_psd(ne, screen_width, radio_frequency, kx, OUTER_SCALE)
    elif PSD_TYPE == 1:
        f = von_karman_dphi_1d_psd(ne, screen_width, radio_frequency, kx, OUTER_SCALE)
    elif PSD_TYPE == 2:
        f = shkarovsky_dphi_1d_psd(ne, screen_width, radio_frequency, kx, OUTER_SCALE, INNER_SCALE)
    else:
        raise ValueError("Unknown PSD type")
        exit(1)

    f[0] = 1e-6  
    g_real = np.random.normal(0, 1, N)
    g_imag = np.random.normal(0, 1, N)
    herm = (g_real + 1j * g_imag) / np.sqrt(2.0)

    phi_k = np.sqrt(ifftshift(f)) * herm
    phase_screen_complex = np.sqrt(2.0 * np.pi / LN) * ifft(phi_k, norm='forward')
    return phase_screen_complex.real


# === Wave propagation ===
def propagate_wave(amp_field, radio_frequency, screen_width, phase_screen, grid):
    k_wave = 2.0 * np.pi * radio_frequency / C_LIGHT
    kx = ifftshift(grid['kx'])

    # Adding noize
    noize = np.random.uniform(-1e-6, 1e-6, len(amp_field))
    amp_field = amp_field * np.exp(1j * (phase_screen + noize))

 
    fft_field = fft(amp_field, norm='forward')
    propagator = np.exp(1j * (kx**2) / (2.0 * k_wave) * screen_width)
    fft_field_prop = fft_field * propagator
    amp_field_out = ifft(fft_field_prop, norm='forward')

    return amp_field_out
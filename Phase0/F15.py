cd ~/lrd_workspace/analysis
python3 << 'EOF'
import numpy as np
import os

# --- Generate synthetic dataset of spectra with known z ---
# Input: Hβ + [OIII]4959 + [OIII]5007 + continuum + noise
# Output: true z

np.random.seed(42)
N = 5000  # number of spectra

# Observed wavelength grid
wave_obs = np.linspace(29000, 53000, 1500)  # Å (G395M)

# Parameters to sample
z_true = np.random.uniform(5.0, 7.5, N)
A_hb   = np.random.uniform(0.5, 5.0, N)
A_4959 = A_hb * 2.5 * np.random.uniform(0.8, 1.2, N)
A_5007 = A_4959 * np.random.uniform(2.7, 3.0, N)
sigma  = np.random.uniform(2.0, 4.0, N)
cont   = np.random.uniform(0.0, 0.1, N)

# Rest-frame lines
LAM_HB, LAM_4959, LAM_5007 = 4862.68, 4960.30, 5008.24

spectra = np.zeros((N, len(wave_obs)), dtype=np.float32)

for i in range(N):
    # Convert to observed
    lam_hb_obs   = LAM_HB   * (1 + z_true[i])
    lam_4959_obs = LAM_4959 * (1 + z_true[i])
    lam_5007_obs = LAM_5007 * (1 + z_true[i])
    
    # Model
    s = cont[i] * np.ones_like(wave_obs)
    s += A_hb[i]   * np.exp(-0.5 * ((wave_obs - lam_hb_obs)   / sigma[i])**2)
    s += A_4959[i] * np.exp(-0.5 * ((wave_obs - lam_4959_obs) / sigma[i])**2)
    s += A_5007[i] * np.exp(-0.5 * ((wave_obs - lam_5007_obs) / sigma[i])**2)
    
    # Noise
    noise_level = 0.01
    s += np.random.normal(0, noise_level, len(wave_obs))
    
    spectra[i] = s.astype(np.float32)

# Normalize
spectra = (spectra - spectra.mean(axis=1, keepdims=True)) / spectra.std(axis=1, keepdims=True)

# Save
out_dir = os.path.expanduser('~/lrd_workspace/analysis/ml_data')
os.makedirs(out_dir, exist_ok=True)
np.save(os.path.join(out_dir, 'spectra.npy'), spectra)
np.save(os.path.join(out_dir, 'z_true.npy'), z_true.astype(np.float32))
print(f"✅ Dataset: {spectra.shape} spectra")
print(f"   z range: {z_true.min():.2f} - {z_true.max():.2f}")
print(f"   Saved to: {out_dir}")
EOF
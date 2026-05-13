#!/usr/bin/env python
import numpy as np
import xarray as xr
from scipy.special import gamma
from scipy.stats import gamma as gamma_dist

# ---------- Functions ----------
# Initially written with GPT assistancei by SB
def lambda_from_nr_qr(nr, qr, rho_water, rho_air, mu):
    """Compute slope parameter lambda for CASIM gamma PSD."""
    return ((np.pi * nr * rho_water * (mu + 3) * (mu + 2) * (mu + 1)) /
            (6.0 * qr * rho_air )) ** (1.0/3.0)

def sample_mass_weighted(nr, qr, rho_water, rho_air, mu, size=1):
    """Draw diameters weighted by mass from CASIM gamma PSD."""
    lam = lambda_from_nr_qr(nr, qr, rho_water, rho_air, mu)
    shape_mass = mu + 4.0
    scale = 1.0 / lam
    return gamma_dist.rvs(a=shape_mass, scale=scale, size=size)

def equivalent_number_concentration(D, qr, rho_water, rho_air):
    """Equivalent N_r if all mass q were in droplets of diameter D."""
    mass_per_drop = (np.pi / 6.0) * rho_water * D**3
    return rho_air * qr / mass_per_drop

r_plume = 500.0
centre = [6395, 3000]
dx_parcel = 10.0 # Distance between parcels

# shape factor for spheroids; >1 is prolate, <1 is oblate
shape_factor = 1

# spheroid axes: make `a` the x-radius and `b` the z-radius.
# For shape_factor > 1 we elongate in x (prolate), < 1 elongates in z.
a = r_plume / shape_factor
b = r_plume * shape_factor

# Ensure the grid of parcel shifts covers the full ellipse (use max radius)
n_steps = int(np.ceil(max(a, b) / dx_parcel))
parcel_shifts = np.arange(-n_steps, n_steps + 1) * dx_parcel
n_shifts=len(parcel_shifts)

qr_parcels=0.001
Nr_parcels=60000
mu=2.5
rho_water = 1000.0
rho_air = 1.2256 
multiplicity=10
parcel_volume=dx_parcel*dx_parcel

i_parcel=0
# oversized_initialisation
x_array=np.zeros((1,n_shifts*n_shifts*multiplicity))
z_array=np.zeros((1,n_shifts*n_shifts*multiplicity))

for ii in range(n_shifts):
    parcel_shift_x=parcel_shifts[ii]
    for jj in range(n_shifts):
        parcel_shift_z=parcel_shifts[jj]
        if ((parcel_shift_x / a) ** 2 + (parcel_shift_z / b) ** 2) < 1.0:
             for mm in range(multiplicity):
                 x_array[0,i_parcel] = centre[0] + parcel_shift_x
                 z_array[0,i_parcel] = centre[1] + parcel_shift_z
                 i_parcel=i_parcel+1

x_array=x_array[:,:i_parcel]
z_array=z_array[:,:i_parcel]

len_parcels=np.shape(x_array)[1]

volume_array=np.ones((1,len_parcels))*parcel_volume
qr_array=np.ones((1,len_parcels))*qr_parcels/multiplicity
Nr_array=np.ones((1,len_parcels))

D_samples = sample_mass_weighted(Nr_parcels, qr_parcels, rho_water, rho_air, mu, len_parcels)
Nr_array[0,:] = equivalent_number_concentration(D_samples, qr_parcels, rho_water, rho_air)/multiplicity
n_parcels = np.arange(1, len_parcels+1,dtype=np.int32)

time = np.array([0.0])

# Coordinates
coords = {
    "time": ("time", time, {
        "units": "seconds since 1970-01-01 00:00:00",
        "calendar": "proleptic_gregorian"
    }),
    "n_parcels": ("n_parcels", n_parcels)
}

# Create the actual dataset
ds = xr.Dataset(
    {
        "x_position": xr.DataArray(x_array, dims=["time", "n_parcels"], coords=coords,
                              attrs={"units": "m", "long_name": "x position component"}),
        "z_position": xr.DataArray(z_array, dims=["time", "n_parcels"], coords=coords,
                              attrs={"units": "m", "long_name": "z position component"}),
        "volume": xr.DataArray(volume_array, dims=["time", "n_parcels"], coords=coords,
                              attrs={"units": "m^2", "long_name": "parcel volume"}),
        "qr": xr.DataArray(qr_array, dims=["time", "n_parcels"], coords=coords,
                              attrs={"units": "kg/kg", "long_name": "rain mixing ratio"}),
        "Nr": xr.DataArray(Nr_array, dims=["time", "n_parcels"], coords=coords,
                              attrs={"units": "1/kg", "long_name": "rain number concentration"})
    }
)

print(sum(Nr_array[0,:]))

# Save with unlimited time dimension
ds.to_netcdf("1g60000_stochastic_rain_input.nc", unlimited_dims=["time"])

# Quick x-z scatter (small, optional):
try:
    import matplotlib.pyplot as plt

    x_vals = x_array[0, :].ravel()
    z_vals = z_array[0, :].ravel()

    plt.figure(figsize=(6, 6))
    plt.scatter(x_vals, z_vals, s=6, alpha=0.8)
    plt.gca().set_aspect("equal", adjustable="box")
    plt.xlabel("x position (m)")
    plt.ylabel("z position (m)")
    plt.title("Parcel positions (x vs z)")
    plt.tight_layout()
    plt.savefig("xz_scatter.png", dpi=150)
    plt.close()
    print("Saved xz_scatter.png")
except Exception as _err:
    # plotting is optional; continue if matplotlib isn't available
    pass

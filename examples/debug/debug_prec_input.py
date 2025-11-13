#!/usr/bin/env python
import numpy as np
import xarray as xr

centre = [3000, 3000]
dx_parcel = 5  # Distance between parcels

# Only one parcel at the centre
x_array = np.array([[centre[0]]])
z_array = np.array([[centre[1]]])

qr_parcels = 0.002
Nr_parcels = 10000
theta_parcels = 300.0

#Some constants for caluclting qv
surf_press = 1.0e5
pressure_scale_height = 7000.0
ref_press = 1.0e5
r_d = 287.04
c_p = 1004.0
RH = 0.7  # relative humidity

#For calculating qv
press = surf_press * np.exp(- (z_array) / pressure_scale_height)
exn = (press / ref_press)**(r_d / c_p)
temp = theta_parcels * exn
ws = 3.8 / ((0.01 * press) * np.exp(-17.2693882 * (temp - 273.15) / (temp - 35.86)) - 6.109)
qv_parcels = float(RH * ws[0, 0])



parcel_volume = dx_parcel * dx_parcel

volume_array = np.array([[parcel_volume]])
qr_array = np.array([[qr_parcels]])
Nr_array = np.array([[Nr_parcels]])
theta_array = np.array([[theta_parcels]])
qv_array = np.array([[qv_parcels]])

n_parcels = np.array([1], dtype=np.int32)
time = np.array([0.0])

# Coordinates
dims = ["time", "n_parcels"]
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
        "x_position": xr.DataArray(x_array, dims=dims, coords=coords,
                              attrs={"units": "m", "long_name": "x position component"}),
        "z_position": xr.DataArray(z_array, dims=dims, coords=coords,
                              attrs={"units": "m", "long_name": "z position component"}),
        "volume": xr.DataArray(volume_array, dims=dims, coords=coords,
                              attrs={"units": "m^2", "long_name": "parcel volume"}),
        "qr": xr.DataArray(qr_array, dims=dims, coords=coords,
                              attrs={"units": "kg/kg", "long_name": "rain mixing ratio"}),
        "Nr": xr.DataArray(Nr_array, dims=dims, coords=coords,
                              attrs={"units": "1/kg", "long_name": "rain number concentration"}),
        "theta": xr.DataArray(theta_array, dims=dims, coords=coords,
                              attrs={"units": "K", "long_name": "potential temperature"}),
        "qv": xr.DataArray(qv_array, dims=dims, coords=coords,
                              attrs={"units": "kg/kg", "long_name": "water vapor mixing ratio"})
    }
)
print(ds)

# Save with unlimited time dimension
ds.to_netcdf("debug_prec_input.nc", unlimited_dims=["time"])
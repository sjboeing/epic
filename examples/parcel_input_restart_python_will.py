from tools.nc_parcels import nc_parcels
import numpy as np
try:
    ncp = nc_parcels()
    ncp.open('negative_theta_parcels_restart.nc')
    # Set the grid size
    nx = 32
    nz = 32
    # Set the origin
    origin = (0.0, 0.0)
    # Set the extent
    extent = (6000, 6000)
    # Set the domain centre
    centre = (3000, 3000)
    # Set the grid spacing
    dx = extent[0] / nx
    dz = extent[1] / nz
    # Set the radius of the bubble
    r_bubble = 1000.0
    # parcels per grid box in each dimension
    n_par_res = 2
    parcels_per_x = nx * n_par_res
    parcels_per_z = nz * n_par_res
    num_parcel = parcels_per_x * parcels_per_z
    dx_parcel = extent[0] / parcels_per_x
    dz_parcel = extent[1] / parcels_per_z
    d_parcel = np.array((dx_parcel, dz_parcel))
    # Initialise arrays
    position = np.zeros((num_parcel, 2))
    theta = np.zeros(num_parcel)
    volume = np.ones(num_parcel) * d_parcel[0] * d_parcel[1]
    b_diag = (0.75 * d_parcel[0] * d_parcel[1] / np.pi) ** (2.0 / 2.0)
    B = np.zeros((num_parcel, 3))
    B[:, 0] = b_diag
    B[:, 2] = b_diag
    vorticity = np.zeros((num_parcel, 2))
    theta_env = 300.0
    theta_pert = -2.0
    iparcel = 0
    for ix in range(parcels_per_x):
        for iz in range(parcels_per_z):
            pos = origin + d_parcel * np.array([ix + 0.5, iz + 0.5])
            position[iparcel, :] = pos
            r = np.sqrt((pos[0] - centre[0]) ** 2 + (pos[1] - centre[1]) ** 2)
            if r < r_bubble:
                theta[iparcel] = theta_env + theta_pert
            else:
                theta[iparcel] = theta_env
            iparcel += 1
    # Write datasets
    ncp.add_dataset('x_position', position[:, 0], unit='m')
    ncp.add_dataset('z_position', position[:, 1], unit='m')
    ncp.add_dataset('theta', theta, unit='K')
    ncp.add_dataset('volume', volume, unit='m^2')
    ncp.add_dataset('x_vorticity', vorticity[:, 0], unit='1/s')
    ncp.add_dataset('z_vorticity', vorticity[:, 1], unit='1/s')
    ncp.add_dataset('B11', B[:, 0], unit='m^2')
    ncp.add_dataset('B12', B[:, 1], unit='m^2')
    ncp.add_dataset('B22', B[:, 2], unit='m^2')
    ncp.add_box(origin, extent, [nx, nz])
    ncp.close()
except Exception as err:
    print(err)

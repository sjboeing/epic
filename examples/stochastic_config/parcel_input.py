try:
    from tools.nc_parcels import nc_parcels
    import numpy as np
    import math 

    nx = 1280
    nz = 500
    # Set the origin
    origin = np.array((0.0, 0.0))
    # Set the extent
    extent = np.array((12800,5000))
    #Set the domin centre
    centre = np.array((6400,3000))
    
    #Set the radius of the bubble
    r_bubble = 500
    #initialise the potential temperature
    
    theta_env = 300.0
    theta_pert = -2.0
    
    
    #Set other constants
    surf_press = 1.0e5
    pressure_scale_height = 7000.0
    ref_press = 1.0e5
    r_d = 287.04
    c_p = 1004.0
    RH = 0.70  # relative humidity
    RH_bubble = 0.95  # saturated inside the bubble


    n_par_res = 3
    parcels_per_dim = np.array([n_par_res*nx,n_par_res*nz])
    dx_parcel = extent/parcels_per_dim
    tuple_ncells = (np.int32(nx),np.int32(nz))

    ncp = nc_parcels()
    ncp.open('EV10.nc')
    ncp.add_box(origin=origin, ncells=tuple_ncells, extent=extent)



    num_parcel = parcels_per_dim[0]*parcels_per_dim[1]
    position = np.zeros((num_parcel, 2))
    theta = np.zeros(num_parcel)
    qv = np.zeros(num_parcel)
    volume = np.ones(num_parcel)*dx_parcel[0]*dx_parcel[1]
    vorticity = np.zeros(num_parcel)
    b_diag = dx_parcel[0] * dx_parcel[1] / np.pi
   
    B = np.zeros((num_parcel, 3))
    B[:, 0] = b_diag   # xx
    B[:, 2] = b_diag   # yy
        # Choose a stretch factor along x (s > 1 stretches, s < 1 squashes)
    s = 1.0  # example

    # Compute semi-axes to keep area same
    a = r_bubble / s      # semi-axis along x
    b = r_bubble * s        # semi-axis along z

    # --- Hydrostatic pressure profile from reference potential-temperature ---
    # prepare vertical levels for parcels (one level per parcel row)
    Pz = int(parcels_per_dim[1])
    z_levels = origin[1] + (np.arange(Pz) + 0.5) * dx_parcel[1]

    # reference potential temperature profile at those levels
    # (use uniform theta_env here; replace with an array if you have vertical profile)
    theta_ref_levels = np.full(Pz, theta_env)

    # constants for hydrostatic integration
    g = 9.80665
    alpha = r_d / c_p

    def compute_pressure_profile(z_levels, theta_ref, p_ref, p_surface, R, cp, g):
        # Transform variable Y = (p/p_ref)^alpha, integrate downward from top
        K = len(z_levels)
        # initial guess for top pressure (small)
        ptop = max(1e-6, p_surface * 1e-3)
        for ipass in range(2):
            Y = np.zeros(K)
            Y[-1] = (ptop / p_ref) ** alpha
            # integrate downward (from top index K-1 to 0)
            for k in range(K-2, -1, -1):
                theta_bar = 0.5 * (theta_ref[k] + theta_ref[k+1])
                dz = z_levels[k+1] - z_levels[k]
                delta_y = (g * dz) / (cp * theta_bar)
                Y[k] = Y[k+1] + delta_y
            p = p_ref * (Y ** (1.0 / alpha))
            # estimate near-surface pressure p0 from first two levels
            if K > 1:
                p0 = 0.5 * (p[0] + p[1])
            else:
                p0 = p[0]
            if ipass == 0:
                Y_top = (ptop / p_ref) ** alpha
                Y0 = (p0 / p_ref) ** alpha
                Y_target = (p_surface / p_ref) ** alpha
                Y_top_new = Y_target + Y_top - Y0
                if Y_top_new <= 0:
                    Y_top_new = 1e-12
                ptop = p_ref * (Y_top_new ** (1.0 / alpha))
        return p

    # compute pressure at each parcel-row level and saturation mixing ratio (kg/kg)
    p_profile = compute_pressure_profile(z_levels, theta_ref_levels, ref_press, surf_press, r_d, c_p, g)
    # temperature corresponding to theta_ref at each level: T = theta_ref * (p/p_ref)^alpha
    T_profile = theta_ref_levels * (p_profile / ref_press) ** alpha
    # Magnus formula for saturation vapour pressure (Pa) and saturation mixing ratio q_s (kg/kg)
    eps = 0.622
    es_hPa = 6.1094 * np.exp(17.2693882 * (T_profile - 273.15) / (T_profile - 35.86))
    es_Pa = es_hPa * 100.0
    p_hPa = p_profile * 0.01
    # avoid division by tiny number
    denom = np.maximum(p_profile - es_Pa, 1e-6)
    qsat_profile = eps * es_Pa / denom

    iparcel = 0
    for i in range(parcels_per_dim[0]):
        for j in range(parcels_per_dim[1]):

            pos = origin + dx_parcel * np.array([i+0.5, j+0.5])
            position[iparcel, :] = pos

            # Elliptical normalized coordinates
            delx = (pos[0] - centre[0]) / a
            delz = (pos[1] - centre[1]) / b

            # Elliptical radius
            rad = np.sqrt(delx*delx + delz*delz)

            # Cosine taper following ellipse with an inner core + taper region
            # inner_core is the fraction of the ellipse radius that is fully perturbed,
            # the region between inner_core and 1.0 is tapered smoothly to zero.
            inner_core = 0.8  # 1 - taper_width (previously attempted 0.2)

            # use hydrostatic profile computed per parcel-row
            level_idx = int(j)
            press = p_profile[level_idx]
            temp = T_profile[level_idx]
            qsat = qsat_profile[level_idx]

            # Ellipse mask with core + taper
            if rad < inner_core:
                # fully inside core: full perturbation* (1 - 0.4*delx*delz)
                theta[iparcel] = theta_env + theta_pert * (1 - 0.4*delx*delz)
                qv[iparcel] = RH_bubble * qsat
            elif rad < 1.0:
                # taper region: map rad from [inner_core, 1.0] -> [0, 1]
                t = (rad - inner_core) / (1.0 - inner_core)
                xi = math.cos(0.5 * math.pi * t)
                xi = max(min(xi, 1.0), 0.0)
                theta[iparcel] = theta_env + theta_pert *  (1 - 0.4*delx*delz) * xi * xi
                # smoothly blend humidity between bubble and environment
                qv[iparcel] = (RH_bubble * xi + RH * (1.0 - xi)) * qsat
            else:
                # outside ellipse: environment
                theta[iparcel] = theta_env
                qv[iparcel] = RH * qsat

            iparcel += 1


    # write all provided datasets
    ncp.add_dataset('theta', theta, unit='K')
    ncp.add_dataset('vorticity', vorticity, unit='s^-1')
    ncp.add_dataset('qv', qv, unit='kg kg^-1')
    ncp.add_dataset("volume", volume)
    ncp.add_dataset('x_position', position[:,0])
    ncp.add_dataset('z_position', position[:,1])
    ncp.add_dataset('B11', B[:, 0], unit='m^2')
    ncp.add_dataset('B12', B[:, 1], unit='m^2')
    ncp.add_dataset('B22', B[:, 2], unit='m^2')

    import matplotlib.pyplot as plt
    # --- Plot vertical profiles: RH, potential temperature, virtual temperature, pressure, density ---
    RH_profile = np.full_like(p_profile, RH)
    qv_profile_env = RH_profile * qsat_profile
    Tv_profile = T_profile * (1.0 + 0.61 * qv_profile_env)
    rho_profile = p_profile / (r_d * Tv_profile)

    # smaller multi-panel: 2 rows x 3 columns
    # print qv at 3150 m (nearest profile level)
    target_z = 3150.0
    idx_near = int(np.argmin(np.abs(z_levels - target_z)))
    print(f"qv at {z_levels[idx_near]:.1f} m = {qv_profile_env[idx_near]:.6e} kg/kg")

    fig, axes = plt.subplots(2, 3, figsize=(12, 6), sharey=True)
    axes = axes.flatten()

    axes[0].plot(RH_profile, z_levels)
    axes[0].set_xlabel('RH')
    axes[0].grid(True)

    axes[1].plot(theta_ref_levels, z_levels)
    axes[1].set_xlabel('theta (K)')
    axes[1].grid(True)

    # show qv profile (mixing ratio)
    axes[2].plot(qv_profile_env, z_levels)
    axes[2].set_xlabel('qv (kg/kg)')
    axes[2].grid(True)

    axes[3].plot(Tv_profile, z_levels)
    axes[3].set_xlabel('Tv (K)')
    axes[3].grid(True)

    axes[4].plot(p_profile / 100.0, z_levels)
    axes[4].set_xlabel('p (hPa)')
    axes[4].grid(True)

    axes[5].plot(rho_profile, z_levels)
    axes[5].set_xlabel('rho (kg/m3)')
    axes[5].grid(True)

    axes[0].set_ylabel('z (m)')
    axes[0].invert_yaxis()
    fig.suptitle('Profiles vs height')
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()

    # # Also show parcel scatter colored by theta
    # plt.figure(figsize=(8, 4))
    # sc = plt.scatter(position[:, 0], position[:, 1], c=theta, s=1, cmap='viridis', marker='o')
    
    # plt.colorbar(sc, label='theta (K)')
    # plt.xlabel('x')
    # plt.ylabel('z')
    # plt.title('Parcel positions colored by theta')
    # plt.gca().set_aspect('equal', adjustable='box')
    # plt.tight_layout()
    # plt.show()

    # # Parcel scatter colored by qv
    # plt.figure(figsize=(8, 4))
    # sc = plt.scatter(position[:, 0], position[:, 1], c=qv, s=1, cmap='viridis', marker='o')
    # plt.colorbar(sc, label='qv (kg/kg)')
    # plt.xlabel('x')
    # plt.ylabel('z')
    # plt.title('Parcel positions colored by qv')
    # plt.gca().set_aspect('equal', adjustable='box')
    # plt.tight_layout()
    # plt.show()
    

    ncp.close()

except Exception as err:
    print(err)

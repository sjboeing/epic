try:
    from tools.nc_parcels import nc_parcels
    import numpy as np
    

    nx = 256
    nz = 100
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
    theta_pert = -2
    
    
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
    ncp.open('P_stochastic_parcel_input.nc')
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
    s = 2.0  # example

    # Compute semi-axes to keep area same
    a = r_bubble / s      # semi-axis along x
    b = r_bubble * s        # semi-axis along z

    iparcel = 0
    for i in range(parcels_per_dim[0]):
        for j in range(parcels_per_dim[1]):

            pos = origin + dx_parcel*np.array([i+0.5,j+0.5])
            position[iparcel,:] = pos
            

            press = surf_press * np.exp(-pos[1] / pressure_scale_height)
            exn = (press / ref_press)**(r_d / c_p)
            temp = theta_env * exn

            ws = 3.8 / (
                (0.01 * press)
                * np.exp(-17.2693882 * (temp - 273.15) / (temp - 35.86))
                - 6.109
            )

            # Normalized distance for ellipse (unchanged)
            r_ellipse = np.sqrt(
                ((pos[0] - centre[0]) / a)**2
                + ((pos[1] - centre[1]) / b)**2
            )

            if r_ellipse < 1.0:
                theta[iparcel] = theta_env + theta_pert
                qv[iparcel] = RH_bubble * ws
            else:
                theta[iparcel] = theta_env
                qv[iparcel] = RH * ws

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

    plt.figure(figsize=(8, 4))
    sc = plt.scatter(position[:, 0], position[:, 1], c=theta, s=1, cmap='viridis', marker='o')
    plt.colorbar(sc, label='theta (K)')
    plt.xlabel('x')
    plt.ylabel('z')
    plt.title('Parcel positions colored by theta')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.tight_layout()
    plt.savefig('parcels_theta.png', dpi=300)
    plt.show()
    

    ncp.close()

except Exception as err:
    print(err)

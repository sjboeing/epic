try:
    from tools.nc_fields import nc_fields
    import numpy as np
    #Set the path to the netcdf file
    ncf = nc_fields()
    ncf.open('191125_field_input.nc')
    #Set the grid size

    nx = 160
    nz = 160
    # Set the origin
    origin = (0.0,0.0)
    # Set the extent
    extent = (6000,6000)
    #Set the domin centre
    centre = (3000,3000)
    #Set the grid spacing
    dx = extent[0]/nx
    dz = extent[1]/nz
    ngrid = [nx,nz]
    #Set the radius of the bubble
    r_bubble = 500
    #initialise the potential temperature
    theta = np.zeros((nz+1,nx))
    theta_env = 300.0
    theta_pert = 0
    #initialise vorticity
    vorticity = np.zeros((nz+1,nx))
    #initialise specific humidity
    qv = np.zeros((nz+1,nx))
    #Set other constants
    surf_press = 1.0e5
    pressure_scale_height = 7000.0
    ref_press = 1.0e5
    r_d = 287.04
    c_p = 1004.0
    RH = 0.7  # relative humidity
    RH_bubble = 0.95  # saturated inside the bubble

        # Choose a stretch factor along x (s > 1 stretches, s < 1 squashes)
    s = 1  # example

    # Compute semi-axes to keep area same
    a = r_bubble * s        # semi-axis along x
    b = r_bubble / s        # semi-axis along z

    for i in range(nx):
        for j in range(nz+1):
            x = origin[0] + (i + 0.5) * dx
            z = origin[1] + j * dz
            press = surf_press * np.exp(- (z) / pressure_scale_height)
            exn = (press / ref_press)**(r_d / c_p)
            temp = theta_env * exn
            ws = 3.8 / ((0.01 * press) * np.exp(-17.2693882 * (temp - 273.15) / (temp - 35.86)) - 6.109)
            
            # Normalized distance for ellipse
            r_ellipse = np.sqrt(((x - centre[0]) / a)**2 + ((z - centre[1]) / b)**2)

            if r_ellipse < 1.0:
                theta[j, i] = theta_env + theta_pert
                qv[j, i] = RH_bubble*ws  # saturated inside the bubble
            else:
                theta[j, i] = theta_env
                qv[j, i] = RH*ws

    # write all provided fields
    ncf.add_field('theta', theta, unit='K')
    ncf.add_field('vorticity', vorticity, unit='s^-1')
    ncf.add_field('qv', qv, unit='kg kg^-1')
    x_array = origin[0] + (np.arange(nx) + 0.5) * dx
    z_array = origin[1] + np.arange(nz+1) * dz
    ncf.add_axis('x', x_array)
    ncf.add_axis('z', z_array)
    ncf.add_axis('t', [0.0])
    
    ncf.add_box(origin, extent, ngrid)

    

    ncf.close()

except Exception as err:
    print(err)

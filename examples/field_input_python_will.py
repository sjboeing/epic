try:
    from tools.nc_fields import nc_fields
    import numpy as np
    #Set the path to the netcdf file
    ncf = nc_fields()
    ncf.open('negative_theta_bubble_field.nc')
    #Set the grid size

    nx = 32
    nz = 32
    # Set the origin
    origin = (0.0,0.0)
    # Set the extent
    extent = (6000,6000)
    #Set the domin centre
    centre = (3000,3000)
    #Set the grid spacing
    dx = extent[0]/nx
    dz = extent[1]/nz
    ngrid = (np.int32(nx),np.int32(nz))
    #Set the radius of the bubble
    r_bubble = 1000.0
    #initialise the potential temperature
    theta = np.zeros((nz+1,nx))
    theta_env = 300.0
    theta_pert = -2.0 

    for i in range(nx):
        for j in range(nz+1):
            x = origin[0] + (i+0.5)*dx
            z = origin[1] + j*dz
            r = np.sqrt((x-centre[0])**2+(z-centre[1])**2)
            if(r<r_bubble):
                theta[j,i] = theta_env + theta_pert
            else:
                theta[j,i] = theta_env

    # write all provided fields
    ncf.add_field('theta', theta, unit='K')

    ncf.add_box(origin, extent, ngrid)

    ncf.close()

except Exception as err:
    print(err)
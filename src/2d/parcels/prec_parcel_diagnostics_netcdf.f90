! =============================================================================
!                      Write prec parcel diagnostics to NetCDF
! =============================================================================
module prec_parcel_diagnostics_netcdf
    use prec_parcel_diagnostics
    use netcdf_utils
    use netcdf_writer
    use netcdf_reader
    use precipitation_parcels, only : prec_parcels, n_prec_parcels
    use prec_parcel_diagnostics
    use parameters, only : lower, extent, nx, nz
    use config, only : package_version, cf_version
    use omp_lib
    use timer, only : start_timer, stop_timer
    use options, only : write_netcdf_options
    use physics, only : write_physical_quantities
    use parcel_types, only : summed_precipitation, summed_deletion
    use constants, only : one
    implicit none

    private

    integer :: prec_parcel_stats_io_timer

    character(len=512) :: ncfname
    integer            :: ncid
    integer            :: t_axis_id, t_dim_id, n_writes,            &
                          n_prec_par_id,                            &
                          int_qr_id,                                &
                          zqr_id,                                   &
                          int_Nr_id,                                &
                          summed_precipitation_id,                  &
                          summed_deletion_id
    double precision   :: restart_time

    public :: create_netcdf_prec_parcel_stats_file,  &
              write_netcdf_prec_parcel_stats,        &
              prec_parcel_stats_io_timer

    contains

        ! Create the parcel diagnostic file.
        ! @param[in] basename of the file
        ! @param[in] overwrite the file
        subroutine create_netcdf_prec_parcel_stats_file(basename, overwrite, l_restart)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical,      intent(in)  :: l_restart
            logical                   :: l_exist

            ncfname =  basename // '_prec_parcel_stats.nc'

            call exist_netcdf_file(ncfname, l_exist)

            restart_time = -one
            n_writes = 1

            if (l_restart .and. l_exist) then
                call open_netcdf_file(ncfname, NF90_NOWRITE, ncid)
                call get_num_steps(ncid, n_writes)
                call get_time(ncid, restart_time)
                call read_netcdf_prec_parcel_stats_content
                call close_netcdf_file(ncid)
                n_writes = n_writes + 1
                return
            endif

            call create_netcdf_file(ncfname, overwrite, ncid)

            ! define global attributes
            call write_netcdf_info(ncid=ncid,                    &
                                   version_tag=package_version,  &
                                   file_type='prec_parcel_stats',     &
                                   cf_version=cf_version)

            call write_netcdf_box(ncid, lower, extent, (/nx, nz/))

            call write_physical_quantities(ncid)

            call write_netcdf_options(ncid)

            call define_netcdf_temporal_dimension(ncid, t_dim_id, t_axis_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='n_prec_parcels',                                      &
                long_name='number of prec parcels',                         &
                std_name='',                                                &
                unit='1',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=n_prec_par_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='int_qr',                                              &
                long_name='qr*volume parcels',                              &
                std_name='',                                                &
                unit='m^3',                                                 &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=int_qr_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='z_mean_qr',                                           &
                long_name='mean height of precip',                          &
                std_name='',                                                &
                unit='m',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=zqr_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='int_Nr',                                              &
                long_name='Nr*volume parcels',                              &
                std_name='',                                                &
                unit='1',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=int_Nr_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='summed_precipitation',                                &
                long_name='qr*vol precipitation',                           &
                std_name='',                                                &
                unit='m^3',                                                 &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=summed_precipitation_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='summed_deletion',                                     &
                long_name='qr*vol deleted',                                 &
                std_name='',                                                &
                unit='m^3',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=summed_deletion_id)

            call close_definition(ncid)

        end subroutine create_netcdf_prec_parcel_stats_file

        ! Pre-condition: Assumes an open file
        subroutine read_netcdf_prec_parcel_stats_content

            call get_dim_id(ncid, 't', t_dim_id)

            call get_var_id(ncid, 't', t_axis_id)

            call get_var_id(ncid, 'n_prec_parcels', n_prec_par_id)

            call get_var_id(ncid, 'int_qr', int_qr_id)

            call get_var_id(ncid, 'z_mean_qr', zqr_id)

            call get_var_id(ncid, 'int_Nr', int_Nr_id)

            call get_var_id(ncid, 'summed_precipition', summed_precipitation_id)

            call get_var_id(ncid, 'summed_deletion', summed_deletion_id)

        end subroutine read_netcdf_prec_parcel_stats_content

        ! Write a step in the parcel diagnostic file.
        ! @param[in] t is the time
        subroutine write_netcdf_prec_parcel_stats(t)
            double precision, intent(in)    :: t

            call start_timer(prec_parcel_stats_io_timer)

            if (t <= restart_time) then
                call stop_timer(prec_parcel_stats_io_timer)
                return
            endif

            call open_netcdf_file(ncfname, NF90_WRITE, ncid)

            ! write time
            call write_netcdf_scalar(ncid, t_axis_id, t, n_writes)

            !
            ! write diagnostics
            !
            call write_netcdf_scalar(ncid, n_prec_par_id, n_prec_parcels, n_writes)
            call write_netcdf_scalar(ncid, int_qr_id, int_qr, n_writes)
            call write_netcdf_scalar(ncid, zqr_id, int_zqr/int_qr, n_writes)
            call write_netcdf_scalar(ncid, int_Nr_id, int_Nr, n_writes)
            call write_netcdf_scalar(ncid, summed_precipitation_id, summed_precipitation, n_writes)
            call write_netcdf_scalar(ncid, summed_deletion_id, summed_deletion, n_writes)

            ! increment counter
            n_writes = n_writes + 1

            call close_netcdf_file(ncid)

            call stop_timer(prec_parcel_stats_io_timer)

        end subroutine write_netcdf_prec_parcel_stats

end module prec_parcel_diagnostics_netcdf

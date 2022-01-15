module parcel_hdf5
    use parcel_container, only : parcels, n_parcels
    use options, only : verbose, write_h5_options
    use parameters, only : nx, nz, extent, lower
    use hdf5
    use h5_utils
    use h5_writer
    use h5_reader, only : get_num_steps
    use timer, only : start_timer, stop_timer
    implicit none

    integer :: hdf5_parcel_timer

    ! h5 file handle
    integer(hid_t)     :: h5file_id
    character(len=512) :: h5fname

    private :: h5file_id, h5fname

    contains

        ! Create the parcel file.
        ! @param[in] basename of the file
        ! @param[in] overwrite the file
        subroutine create_h5_parcel_file(basename, overwrite, l_restart, step)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical,      intent(in)  :: l_restart
            integer,      intent(out) :: step
            logical                   :: l_exist

            h5fname =  basename // '_parcels.hdf5'

            call exist_h5_file(h5fname, l_exist)

            if (l_restart .and. l_exist) then
                call open_h5_file(h5fname, H5F_ACC_RDWR_F, h5file_id)
                call get_num_steps(h5file_id, step)
                call close_h5_file(h5file_id)
                return
            endif

            step = 0

            call create_h5_file(h5fname, overwrite, h5file_id)

            call write_h5_scalar_attrib(h5file_id, 'output_type', 'parcels')

            call write_h5_timestamp(h5file_id)
            call write_h5_options(h5file_id)
            call write_h5_box(h5file_id, lower, extent, (/nx, nz/))

            call close_h5_file(h5file_id)

        end subroutine create_h5_parcel_file

        ! Write diagnostics for the current time step in the parcel file.
        ! @param[inout] nw counts the number of writes
        ! @param[in] t is the time
        ! @param[in] dt is the time step
        subroutine write_h5_parcel_step(nw, t, dt)
            integer,          intent(inout) :: nw
            double precision, intent(in)    :: t
            double precision, intent(in)    :: dt

            call start_timer(hdf5_parcel_timer)

#ifdef ENABLE_VERBOSE
            if (verbose) then
                print "(a19)", "write parcels to h5"
            endif
#endif

            call open_h5_file(h5fname, H5F_ACC_RDWR_F, h5file_id)

            call write_h5_scalar_step_attrib(h5file_id, nw, "t", t)

            call write_h5_scalar_step_attrib(h5file_id, nw, "dt", dt)

            call write_h5_scalar_step_attrib(h5file_id, nw, "num parcel", n_parcels)

            call write_h5_parcels(nw)

            ! increment counter
            nw = nw + 1

            ! update number of iterations to h5 file
            call write_h5_num_steps(h5file_id, nw)

            call close_h5_file(h5file_id)

            call stop_timer(hdf5_parcel_timer)

        end subroutine write_h5_parcel_step


        ! Write parcel datasets (called from write_h5_parcel_step).
        ! @param[in] iter is the number of the write
        subroutine write_h5_parcels(iter)
            integer, intent(in)           :: iter ! iteration
            integer(hid_t)                :: group
            character(:), allocatable     :: name
            logical                       :: created

            name = trim(get_step_group_name(iter))

            call create_h5_group(h5file_id, name, group, created)

            if (.not. created) then
                call open_h5_group(h5file_id, name, group)
            endif

            !
            ! write parcel data
            !

            call write_h5_dataset(h5file_id, name, "position", &
                                  parcels%position(1:n_parcels, :))

            call write_h5_dataset(h5file_id, name, "B", &
                                  parcels%B(1:n_parcels, :))

            call write_h5_dataset(h5file_id, name, "volume", &
                                  parcels%volume(1:n_parcels))

            call write_h5_dataset(h5file_id, name, "vorticity", &
                                  parcels%vorticity(1:n_parcels))

            call write_h5_dataset(h5file_id, name, "buoyancy", &
                                  parcels%buoyancy(1:n_parcels))
#ifndef ENABLE_DRY_MODE
            call write_h5_dataset(h5file_id, name, "humidity", &
                                  parcels%humidity(1:n_parcels))
#endif
            call close_h5_group(group)

        end subroutine write_h5_parcels
end module parcel_hdf5

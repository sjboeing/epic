! =============================================================================
!                           Module to split ellipsoids
! =============================================================================
module parcel_split_mod
    use options, only : verbose, parcel
    use constants, only : pi, three, five, f12, f34
    use parameters, only : amax
    use parcel_container, only : parcels                &
                               , n_parcels              &
                               , parcel_resize
    use parcel_bc, only : apply_reflective_bc
    use parcel_ellipsoid, only : diagonalise, get_aspect_ratio, get_eigenvalues
    use timer, only : start_timer, stop_timer
    use omp_lib
    implicit none

    double precision, parameter :: dh = f12 * dsqrt(three / five)

    private :: dh

    integer :: split_timer

    ! number of parcel splits (is reset in every write step)
    integer :: n_parcel_splits = 0


    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Split elongated parcels (semi-major axis larger than amax) or
        ! parcels with aspect ratios larger than parcel%lambda_max.
        subroutine parcel_split
            double precision     :: B(5)
            double precision     :: vol, lam
            double precision     :: D(3), V(3, 3)
            integer              :: last_index, n_indices
            integer              :: grown_size, shrunk_size, n_required
            integer              :: i, n, n_thread_loc
            integer, allocatable :: indices(:)

            call start_timer(split_timer)

            !------------------------------------------------------------------
            ! Check which parcels split and store the indices in *pid*:
            !$omp parallel default(shared)
            !$omp do private(n, B, vol, lam, D)
            do n = 1, n_parcels
                B = parcels%B(:, n)
                vol = parcels%volume(n)

                D = get_eigenvalues(B, vol)

                ! evaluate maximum aspect ratio (a2 >= b2 >= c2)
                lam = get_aspect_ratio(D)

                if (lam < threshold .and. D(1) < amax ** 2) then
                    cycle
                endif

                pid(n) = n

            enddo
            !$omp end do
            !$omp end parallel

            ! contains all indices of parcels that split
            indices = pack(pid(1:n_parcels), pid(1:n_parcels) /= 0)

            n_indices = size(indices)

            !------------------------------------------------------------------
            ! Adapt container size if needed:

            ! we get additional "n_indices" parcels
            n_required = n_parcels + n_indices

            shrunk_size = max(n_required, max_num_parcels)
            shrunk_size = nint(parcel%shrink_factor * shrunk_size)

            call stop_timer(split_timer)

            if (n_required > max_num_parcels) then
                grown_size = nint(parcel%grow_factor * n_required)
                call parcel_resize(grown_size)
            else if (n_required < int(f34 * shrunk_size)) then
                call parcel_resize(shrunk_size)
            endif

            call start_timer(split_timer)

            !------------------------------------------------------------------
            ! Loop over all parcels that really split:

            last_index = n_parcels

            !$omp parallel default(shared)
            !$omp do private(i, n, B, vol, lam, D, V, n_thread_loc)
            do i = 1, n_indices

                ! get parcel index
                n = indices(i)

                B = parcels%B(:, n)
                vol = parcels%volume(n)

                call diagonalise(B, vol, D, V)

                pid(n) = 0

                !
                ! this ellipsoid is split, i.e., add a new parcel
                !
                parcels%B(1, n) = B(1) - f34 * D(1) * V(1, 1) ** 2
                parcels%B(2, n) = B(2) - f34 * D(1) * V(1, 1) * V(2, 1)
                parcels%B(3, n) = B(3) - f34 * D(1) * V(1, 1) * V(3, 1)
                parcels%B(4, n) = B(4) - f34 * D(1) * V(2, 1) ** 2
                parcels%B(5, n) = B(5) - f34 * D(1) * V(2, 1) * V(3, 1)


                parcels%volume(n) = f12 * vol

                !$omp critical
                n_thread_loc = n_parcels + 1

                ! we only need to add one new parcel
                n_parcels = n_parcels + 1
                !$omp end critical


                parcels%B(:, n_thread_loc) = parcels%B(:, n)

                parcels%vorticity(:, n_thread_loc) = parcels%vorticity(:, n)
                parcels%volume(n_thread_loc) = parcels%volume(n)
                parcels%buoyancy(n_thread_loc) = parcels%buoyancy(n)
#ifndef ENABLE_DRY_MODE
                parcels%humidity(n_thread_loc) = parcels%humidity(n)
#endif
                V(:, 1) = V(:, 1) * dh * dsqrt(D(1))
                parcels%position(:, n_thread_loc) = parcels%position(:, n) - V(:, 1)
                parcels%position(:, n) = parcels%position(:, n) + V(:, 1)

                ! child parcels need to be reflected into domain, if their center
                ! is inside the halo region
                call apply_reflective_bc(parcels%position(:, n_thread_loc), &
                                         parcels%B(:, n_thread_loc))

                call apply_reflective_bc(parcels%position(:, n), parcels%B(:, n))

            enddo
            !$omp end do
            !$omp end parallel

            n_parcel_splits = n_parcel_splits + n_parcels - last_index

#ifdef ENABLE_VERBOSE
            if (verbose) then
                print "(a36, i0, a3, i0)", &
                      "no. parcels before and after split: ", last_index, "...", n_parcels
            endif
#endif

            call stop_timer(split_timer)

            ! subtract one call as we start and stop the timer twice here:
            timings(split_timer)%n_calls = timings(split_timer)%n_calls - 1

        end subroutine parcel_split

end module parcel_split_mod

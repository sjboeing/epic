! =============================================================================
!                       Module to determine minimal parcel volume
! =============================================================================
module parcel_vmin
    use options, only : parcel, time
    use parcel_interpl, only : par2grid
    use fields, only : velgradg, velog, vortg, vtend, tbuoyg
    use tri_inversion, only : vor2vel, vorticity_tendency
    use constants, only : zero, f12
    use parameters, only : nx, nz, dxi, vkeep, vcutoff, vmin_dt_factor
    use parcel_interpl, only : trilinear, ngp
    use parcel_container
    use timer, only : start_timer, stop_timer

    implicit none
    integer:: vmin_timer

    ! interpolation indices
    ! (first dimension x, y; second dimension k-th index)
    integer :: is(ngp), js(ngp), n, l

    ! interpolation weights
    double precision :: weights(ngp)

    contains

        ! Merge small parcels into neighbouring equal-sized parcels or bigger
        ! parcels which are close by.
        ! @param[inout] parcels is the parcel container
        subroutine get_parcel_vmin(parcels)
            type(parcel_container_type), intent(inout) :: parcels
            double precision                :: dt

            call par2grid

            ! need to be called in order to set initial time step;
            ! this is also needed for the first ls-rk4 substep
            call vor2vel(vortg, velog, velgradg)

            call vorticity_tendency(tbuoyg, vtend)

            call start_timer(vmin_timer)

            ! update the time step
            dt = get_parcel_time_step(parcels)

            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, n_parcels
              if(parcels%vmin(n)<parcel%vmin_dt_factor*dt) then
                parcels%vmin(n)=vcutoff
              elseif(vcutoff*(parcels%vmin(n)/(parcel%vmin_dt_factor*dt))<vkeep) then
                parcels%vmin(n)=vcutoff*(parcels%vmin(n)/(parcel%vmin_dt_factor*dt))
              else
                parcels%vmin(n)=vkeep
              end if
            end do
            !$omp end do
            !$omp end parallel

            call stop_timer(vmin_timer)

        end subroutine get_parcel_vmin

        ! Estimate a suitable time step based on the velocity strain
        ! and buoyancy gradient.
        ! @returns the parcel time step
        function get_parcel_time_step(parcels) result(dt)
            type(parcel_container_type), intent(inout) :: parcels

            double precision             :: dt
            double precision             :: gmax, bmax
            double precision             :: gfield(0:nz, 0:nx-1),dbdz(0:nz, 0:nx-1),bfield(0:nz, 0:nx-1)
            integer                      :: kx,iz

            gfield = f12 * dsqrt(((velgradg(0:nz, :, 1) - velgradg(0:nz, :, 4)) ** 2 + &
                                        (velgradg(0:nz, :, 2) + velgradg(0:nz, :, 3)) ** 2))

            dbdz(0:nz, :) = f12 * dxi(2) * (tbuoyg(1:nz+1, :) - tbuoyg(-1:nz-1, :))

            bfield = dsqrt(dsqrt(vtend(0:nz, :) ** 2 + dbdz ** 2))

            do kx=0,nx-1
              do iz=0,nz
                gfield(iz,kx)=max(epsilon(gfield(iz,kx)),gfield(iz,kx))
                bfield(iz,kx)=max(epsilon(bfield(iz,kx)),bfield(iz,kx))
              enddo
            enddo

            parcels%vmin(1:n_parcels)=zero

            !$omp parallel default(shared)
            !$omp do private(n, l, is, js, weights)
            do n = 1, n_parcels
              call trilinear(parcels%position(n, :), is, js, weights)
              do l = 1, ngp
                parcels%vmin(n) = parcels%vmin(n) + weights(l)*&
                min(time%alpha / gfield(js(l), is(l)), time%alpha / bfield(js(l), is(l)))
              enddo
            enddo
            !$omp end do
            !$omp end parallel

            dt = min(time%alpha / maxval(gfield), time%alpha / maxval(bfield))

        end function get_parcel_time_step

end module parcel_vmin


! =============================================================================
!                       Module to determine minimal parcel volume
! =============================================================================
module parcel_vmin
    use options, only : parcel, time
    use parcel_interpl, only : par2grid
    use fields, only : velgradg, velog, vortg, vtend, tbuoyg
    use tri_inversion, only : vor2vel, vorticity_tendency, xtrig, xfactors, hrkx
    use constants, only : zero, f12
    use parameters, only : nx, nz, dxi, vkeep, vcutoff, vmin_dt_factor
    use parcel_interpl, only : trilinear, ngp
    use parcel_container
    use stafft
    use deriv1d
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
            double precision             :: dt, vortgrad2max,vortgrad2mag
            double precision             :: vortgrad2magg(-1:nz+1, 0:nx-1)

            call par2grid

            ! need to be called in order to set initial time step;
            ! this is also needed for the first ls-rk4 substep
            call vor2vel(vortg, velog, velgradg)

            call vorticity_tendency(tbuoyg, vtend)

            call start_timer(vmin_timer)

            ! calculate the magnitude of the vorticity gradient
            call get_vortgrad2mag(vortg,vortgrad2magg,vortgrad2max)

            ! update the time step
            call get_parcel_time_step(parcels,dt)

            !$omp parallel default(shared)
            !$omp do private(n,is,js,weights,vort2gradmag)
            do n = 1, n_parcels
              parcels%vmin(n)=vcutoff*(parcels%vmin(n)/(parcel%vmin_dt_factor*dt))
              call trilinear(parcels%position(n, :), is, js, weights)
              vortgrad2mag=zero
              do l = 1, ngp
                vortgrad2mag= vortgrad2mag+weights(l)*vortgrad2magg(js(l), is(l))
              enddo
              parcels%vmin(n)=min(parcels%vmin(n),vcutoff*(vortgrad2max)/(parcel%vmin_dt_factor*vortgrad2mag))
              if(parcels%vmin(n)<vcutoff) then
                parcels%vmin(n)=vcutoff
              elseif(parcels%vmin(n)>vkeep) then
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
        subroutine get_parcel_time_step(parcels,dt)
            type(parcel_container_type), intent(inout) :: parcels
            double precision,intent(out) :: dt
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

        end subroutine get_parcel_time_step

        subroutine get_vortgrad2mag(vortg, vortgrad2magg,vortgrad2max)
            double precision, intent(in)  :: vortg(-1:nz+1, 0:nx-1)
            double precision, intent(out) :: vortgrad2magg(-1:nz+1, 0:nx-1)
            double precision, intent(out) :: vortgrad2max
            double precision              :: dvortdz(-1:nz+1, 0:nx-1)
            double precision              :: psig(0:nz, 0:nx-1)

            psig = vortg(0:nz, 0:nx-1)

            ! Forward x FFT:
            call forfft(nz+1, nx, psig, xtrig, xfactors)

            ! Compute x derivative spectrally of psig:
            call deriv(nz+1, nx, hrkx, psig, vortgrad2magg(0:nz, :))

            ! Reverse x FFT
            call revfft(nz+1, nx, vortgrad2magg(0:nz, :), xtrig, xfactors)

            ! Fill z grid lines outside domain:
            vortgrad2magg(-1,   :) = two * vortgrad2magg(0,  :) - vortgrad2magg(1,    :)
            vortgrad2magg(nz+1, :) = two * vortgrad2magg(nz, :) - vortgrad2magg(nz-1, :)

            dvortdz(0:nz, :)=(f12 * dxi(2) * (vortg(1:nz+1, :) - vortg(-1:nz-1, :)))
            dvortdz(-1,   :) = two * dvortdz(0,  :) - dvortdz(1,    :)
            dvortdz(nz+1, :) = two * dvortdz(nz, :) - dvortdz(nz-1, :)

            vortgrad2magg = vortgrad2magg**2 +dvortdz**2

            vortgrad2max = maxval(vortgrad2magg)

        end subroutine get_vortgrad2mag

end module parcel_vmin


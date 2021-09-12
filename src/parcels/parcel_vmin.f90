! =============================================================================
!                       Module to determine minimal parcel volume
! =============================================================================
module parcel_vmin
    use options, only : parcel, time
    use parcel_interpl, only : par2grid
    use fields, only : velgradg, velog, vortg, vtend, tbuoyg
    use tri_inversion, only : vor2vel, vorticity_tendency, xtrig, xfactors, hrkx
    use constants, only : zero, f12, f32, f23
    use parameters, only : nx, nz, dxi, vkeep, vcutoff, vmin_nu_factor
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

            call par2grid

            ! need to be called in order to set initial time step;
            ! this is also needed for the first ls-rk4 substep
            call vor2vel(vortg, velog, velgradg)

            call vorticity_tendency(tbuoyg, vtend)

            call start_timer(vmin_timer)

            call get_parcel_vmin_part_2(parcels)

            call stop_timer(vmin_timer)

        end subroutine get_parcel_vmin

        ! Estimate a suitable time step based on the velocity strain
        ! and buoyancy gradient.
        ! @returns the parcel time step
        subroutine get_parcel_vmin_part_2(parcels)
            type(parcel_container_type), intent(inout) :: parcels
            double precision             :: gmax, bmax, vortgradmax, vmintemp
            double precision             :: glocal, blocal, vortgradlocal
            double precision             :: gfield(0:nz, 0:nx-1),dbdz(0:nz, 0:nx-1),bfield(0:nz, 0:nx-1)
            double precision             :: vortgradmagg(-1:nz+1, 0:nx-1)
            double precision             :: dvortdz(-1:nz+1, 0:nx-1)
            double precision             :: psig(0:nz, 0:nx-1)
            integer                      :: kx,iz

            psig = vortg(0:nz, 0:nx-1)

            ! Forward x FFT:
            call forfft(nz+1, nx, psig, xtrig, xfactors)

            ! Compute x derivative spectrally of psig:
            call deriv(nz+1, nx, hrkx, psig, vortgradmagg(0:nz, :))

            ! Reverse x FFT
            call revfft(nz+1, nx, vortgradmagg(0:nz, :), xtrig, xfactors)

            ! Fill z grid lines outside domain:
            vortgradmagg(-1,   :) = two * vortgradmagg(0,  :) - vortgradmagg(1,    :)
            vortgradmagg(nz+1, :) = two * vortgradmagg(nz, :) - vortgradmagg(nz-1, :)

            dvortdz(0:nz, :)=(f12 * dxi(2) * (vortg(1:nz+1, :) - vortg(-1:nz-1, :)))
            dvortdz(-1,   :) = two * dvortdz(0,  :) - dvortdz(1,    :)
            dvortdz(nz+1, :) = two * dvortdz(nz, :) - dvortdz(nz-1, :)

            vortgradmagg = dsqrt(vortgradmagg**2 +dvortdz**2)

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

            do kx=0,nx-1
              do iz=0,nz
                vortgradmagg(iz,kx)=max(epsilon(vortgradmagg(iz,kx)),vortgradmagg(iz,kx))
              enddo
            enddo

            gmax=maxval(gfield)
            bmax=maxval(bfield)
            vortgradmax = maxval(vortgradmagg)

            parcels%vmin(1:n_parcels)=zero

            !$omp parallel default(shared)
            !$omp do private(n, l, is, js, weights,glocal,blocal,vmintemp,vortgradlocal)
            do n = 1, n_parcels
              call trilinear(parcels%position(n, :), is, js, weights)
              glocal=zero
              do l = 1, ngp
                glocal = glocal + weights(l)*gfield(js(l), is(l))
              enddo
              vmintemp=vmin_nu_factor*vcutoff*gmax/glocal
              if(bmax>2.0*epsilon(bmax)) then
                blocal=zero
                do l = 1, ngp
                  blocal = blocal + weights(l)*bfield(js(l), is(l))
                enddo
                vmintemp=min(vmintemp,vmin_nu_factor*vcutoff*bmax/blocal)
              endif
              vortgradlocal=zero
              do l = 1, ngp
                vortgradlocal = vortgradlocal + weights(l)*vortgradmagg(js(l), is(l))
              enddo
              vmintemp=min(vmintemp,vcutoff*(vmin_nu_factor*vortgradmax/vortgradlocal)**f23)
              if(vmintemp>vkeep) then
                vmintemp=vkeep
              elseif(vmintemp<vcutoff) then
                vmintemp=vcutoff
              endif
              parcels%vmin(n)=vmintemp
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine get_parcel_vmin_part_2

end module parcel_vmin


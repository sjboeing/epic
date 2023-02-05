! =============================================================================
!                             Field diagnostics
! =============================================================================
module field_diagnostics
    use parameters, only : vcell, vcelli, nx, nz, ngridi, ncelli, vdomaini
    use constants, only : f12, f14
    use fields
    use timer, only : start_timer, stop_timer
    use physics, only : ape_calculation
    use ape_density, only : ape_den
    implicit none

    integer :: field_stats_timer

    double precision :: rms_v,      &       ! rms volume error
                        abserr_v,   &       ! max absolute normalised volume error
                        max_npar,   &       ! max num parcels per cell
                        min_npar,   &       ! min num parcels per cell
                        avg_npar,   &       ! average num parcels per cell
                        avg_nspar,  &       ! average num small parcels per cell
                        keg,        &       ! domain-averaged kinetic energy calculated on the grid
                        apeg,       &       ! domain-average available potential energy on the grid
                        eng,        &       ! domain-averaged enstrophy calculated on the grid
                        min_buoyg,  &       ! minimum gridded buoyancy value
                        max_buoyg           ! maximum gridded buoyancy value
    contains

        subroutine calculate_field_diagnostics
            double precision :: sqerrsum, z(0:nz)
            integer          :: ix, iy, iz

            call start_timer(field_stats_timer)

            ! do not take halo cells into account
            sqerrsum = sum((volg(0:nz, :, :) - vcell) ** 2)
            rms_v = dsqrt(sqerrsum * ngridi) * vcelli

            abserr_v = maxval(abs(volg(0:nz, :, :)  - vcell)) * vcelli

            max_npar = maxval(nparg(0:nz-1, :, :))

            min_npar = minval(nparg(0:nz-1, :, :))

            avg_npar = sum(nparg(0:nz-1, :, :)) * ncelli

            avg_nspar = sum(nsparg(0:nz-1, :, :)) * ncelli

            min_buoyg = minval(tbuoyg(0:nz, :, :))
            max_buoyg = maxval(tbuoyg(0:nz, :, :))

            ! use half weights for boundary grid points
            keg = f12 * sum(volg(1:nz-1, :, :) * ( velog(1:nz-1, :, :, 1) ** 2   &
                                                 + velog(1:nz-1, :, :, 2) ** 2   &
                                                 + velog(1:nz-1, :, :, 3) ** 2)) &
                + f14 * sum(volg(0,  :, :) * ( velog(0,  :, :, 1) ** 2           &
                                             + velog(0,  :, :, 2) ** 2           &
                                             + velog(0,  :, :, 3) ** 2))         &
                + f14 * sum(volg(nz, :, :) * ( velog(nz, :, :, 1) ** 2           &
                                             + velog(nz, :, :, 2) ** 2           &
                                             + velog(nz, :, :, 3) ** 2))

            ! divide by domain volume to get domain-averaged kinetic energy
            keg = keg * vdomaini

            eng = f12 * sum(volg(1:nz-1, :, :) * ( vortg(1:nz-1, :, :, 1) ** 2   &
                                                 + vortg(1:nz-1, :, :, 2) ** 2   &
                                                 + vortg(1:nz-1, :, :, 3) ** 2)) &
                + f14 * sum(volg(0,  :, :) * ( vortg(0,  :, :, 1) ** 2           &
                                             + vortg(0,  :, :, 2) ** 2           &
                                             + vortg(0,  :, :, 3) ** 2))         &
                + f14 * sum(volg(nz, :, :) * ( vortg(nz, :, :, 1) ** 2           &
                                             + vortg(nz, :, :, 2) ** 2           &
                                             + vortg(nz, :, :, 3) ** 2))

            ! divide by domain volume to get domain-averaged enstrophy
            eng = eng * vdomaini


            if (ape_calculation == 'ape density') then
                do iz = 0, nz
                    z(iz) = lower(3) + dble(iz) * dx(3)
                enddo

                apeg = zero
                do ix = 0, nx-1
                    do iy = 0, ny-1
                        apeg = apeg + sum(volg(1:nz-1, iy, ix) * ape_den(tbuoyg(1:nz-1, iy, ix), z(1:nz-1))) &
                             + f12 *      volg(0,      iy, ix) * ape_den(tbuoyg(0,      iy, ix), z(0))       &
                             + f12 *      volg(nz,     iy, ix) * ape_den(tbuoyg(nz,     iy, ix), z(nz))
                    enddo
                enddo

                apeg = apeg * vdomaini
            endif

            call stop_timer(field_stats_timer)

        end subroutine calculate_field_diagnostics

end module field_diagnostics

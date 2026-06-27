! =============================================================================
! This module contains the subroutines to do damping of parcel properties to
! gridded fields in a conservative manner. It does this by nudging all the
! parcels associated with a grid point to the grid point value, with the strength
! of damping proportional to the grid point contribution to the gridded value and
! the strain rate at the grid point.
! =============================================================================
module parcel_damping
    use constants, only :  f14, zero, one
    use timer, only : start_timer, stop_timer
    use parameters, only : nx, nz, vmin, upper, lower
    use dynamic_parcels, only : parcels, n_parcels
    use parcel_ellipse
    use parcel_interpl
    use fields
    use omp_lib
    use options, only : damping
    use parcel_interpl, only : bilinear, par2grid
    use tri_inversion, only : vor2vel
    use rk4_utils, only : get_strain_magnitude_field
    implicit none


    private

    ! interpolation indices
    ! (first dimension x, z; second dimension l-th index)
    integer :: is, js
    integer :: damping_timer

    ! interpolation weights
    double precision :: weights(0:1,0:1)
    double precision :: weight(0:1,0:1)
    double precision :: time_fact(0:1,0:1)

    public :: parcel_damp, damping_timer

    contains

        subroutine parcel_damp(dt)
            double precision, intent(in)  :: dt

            if (damping%l_vorticity .and. damping%l_surface_vorticity) then
                print *, "damping%l_vorticity and  damping%l_surface_vorticity both activated, only one allowed"
                stop
            elseif (damping%l_scalars .and. damping%l_surface_scalars) then
                print *, "damping%l_scalars and  damping%l_surface_scalars both activated, only one allowed"
                stop
            endif
            if (damping%l_vorticity .or. damping%l_scalars .or. &
                damping%l_surface_vorticity .or. damping%l_surface_scalars) then
                ! ensure gridded fields are up to date
                call par2grid
                call vor2vel(vortg, velog, velgradg)

                select type (parcels)
                type is (idealised_parcel_type)
                    call reflect_idealised(parcels)
                type is (realistic_parcel_type)
                    call reflect_realistic(parcels)
                end select

                call get_strain_magnitude_field
                
                select type (parcels)
                type is (idealised_parcel_type)
                    call perturbation_damping_idealised(parcels, dt, .true.)
                type is (realistic_parcel_type)
                    call perturbation_damping_realistic(parcels, dt, .true.)
                end select
            end if

        end subroutine parcel_damp

        subroutine reflect_idealised(parcels)
                class(idealised_parcel_type), intent(in) :: parcels
                ! Reflect beyond boundaries to ensure damping is conservative
                ! This is because the points below the surface contribute to the level above
                !$omp parallel workshare
                vortg(-1,   :) = vortg(1, :)
                vortg(nz+1, :) = vortg(nz-1, :)

#ifndef ENABLE_DRY_MODE
                humg(-1,   :) = humg(1, :)
                humg(nz+1, :) = humg(nz-1, :)
                dbuoyg(-1,   :) = dbuoyg(1, :)
                dbuoyg(nz+1, :) = dbuoyg(nz-1, :)
#else
                tbuoyg(-1,   :) = tbuoyg(1, :)
                tbuoyg(nz+1, :) = tbuoyg(nz-1, :)
#endif
                !$omp end parallel workshare
        end subroutine reflect_idealised

        subroutine reflect_realistic(parcels)
                 class(realistic_parcel_type), intent(in) :: parcels
                ! Reflect beyond boundaries to ensure damping is conservative
                ! This is because the points below the surface contribute to the level above
                !$omp parallel workshare
                vortg(-1,   :) = vortg(1, :)
                vortg(nz+1, :) = vortg(nz-1, :)
                !$omp end parallel workshare

                if(parcels%is_moist) then
                  !$omp parallel workshare
                  qvg(-1,   :) = qvg(1, :)
                  qvg(nz+1, :) = qvg(nz-1, :)
                  qlg(-1,   :) = qlg(1, :)
                  qlg(nz+1, :) = qlg(nz-1, :)
                  !$omp end parallel workshare
                endif
                
                if(parcels%has_droplets) then
                  !$omp parallel workshare
                  Nlg(-1,   :) = Nlg(1, :)
                  Nlg(nz+1, :) = Nlg(nz-1, :)
                  !$omp end parallel workshare
                endif

                !$omp parallel workshare
                thetag(-1,   :) = thetag(1, :)
                thetag(nz+1, :) = thetag(nz-1, :)
                !$omp end parallel workshare
        end subroutine reflect_realistic

        !
        ! @pre: the strain must be calculated and the gridded fields updated
        subroutine perturbation_damping_idealised(parcels, dt, l_reuse)
            class(idealised_parcel_type), intent(inout) :: parcels
            double precision, intent(in)  :: dt
            logical, intent(in)           :: l_reuse
            integer                       :: n, p, l, surface_index
            double precision              :: points(2, 2)
            double precision              :: pvol
            ! tendencies need to be summed up between associated 4 points
            ! before modifying the parcel attribute
            double precision              :: vortend
#ifndef ENABLE_DRY_MODE
            double precision              :: dbuoytend
            double precision              :: humtend
#else
            double precision              :: tbuoytend
#endif
            call start_timer(damping_timer)

            
            !$omp parallel default(shared)
            !$omp do private(n, p, l, points, pvol, weight, surface_index) &
#ifndef ENABLE_DRY_MODE
            !$omp& private(is, js, weights, vortend, humtend, dbuoytend, time_fact)
#else
            !$omp& private(is, js, weights, vortend, tbuoytend, time_fact)
#endif
            do n = 1, n_parcels
                ! check if only surface damping applies and we are far from surfaces
                ! put in a buffer here as parcels can get stretched in integration
                if(.not.(damping%l_vorticity .or. damping%l_scalars)) then
                    if(parcels%position(2, n) > lower(2) + 2 * dx(2)) then
                    if(parcels%position(2, n) < upper(2) - 2 * dx(2)) then
                        cycle
                    end if
                    end if
                endif

                pvol = parcels%volume(n)

                points = get_ellipse_points(parcels%position(:, n), &
                                            pvol, parcels%B(:, n))

                vortend = zero
#ifndef ENABLE_DRY_MODE
                humtend = zero
                dbuoytend = zero
#else
                tbuoytend = zero
#endif

                ! we have 2 points per ellipse
                do p = 1, 2

                    ! ensure point is within the domain
                    call apply_periodic_bc(points(:, p))

                    ! get interpolation weights and mesh indices
                    call bilinear(points(:, p), is, js, weights)

                    ! loop over grid points which are part of the interpolation
                    ! the weight is halved due to 2 points per ellipse
                    weight = f12 * weights * pvol

                    if (damping%l_vorticity) then
                        ! Note this exponential factor can be different for vorticity/scalars
                        time_fact = one - exp(-damping%vorticity_prefactor * strain_mag(js:js+1, is:is+1) * dt)
                        vortend = vortend + sum(weight * time_fact * (vortg(js:js+1, is:is+1) &
                                       - parcels%vorticity(1, n)))
                    endif

                    if (damping%l_scalars) then
                        time_fact = one - exp(-damping%scalars_prefactor * strain_mag(js:js+1, is:is+1) * dt)
#ifndef ENABLE_DRY_MODE
                        humtend = humtend + sum(weight * time_fact * (humg(js:js+1, is:is+1) - parcels%humidity(n)))
                        dbuoytend = dbuoytend + sum(weight * time_fact * (dbuoyg(js:js+1, is:is+1) - parcels%buoyancy(n)))
#else
                        tbuoytend = tbuoytend + sum(weight * time_fact * (tbuoyg(js:js+1, is:is+1) - parcels%buoyancy(n)))
#endif
                    endif

                    if (damping%l_surface_vorticity .or. damping%l_surface_scalars) then
                        ! Index to keep track of grid cells right above/below boundary
                        ! This is because the damping only happens at the boundary level
                        ! Consistent with reflection used in parcel_damp
                        if ((js == -1) .or. (js == nz-1)) then
                            surface_index = 1 ! below lower or below upper boundary
                        elseif ((js == 0) .or. (js == nz)) then
                            surface_index = 0 ! above lower or above upper boundary
                        else
                            cycle ! continue loop if not near a surface
                        endif
                    else
                        cycle ! continue loop if no surface damping
                    endif

                    if (damping%l_surface_vorticity) then
                        ! Note this exponential factor can be different for vorticity/scalars
                        time_fact = one - exp(-damping%vorticity_prefactor * strain_mag(js:js+1, is:is+1) * dt)
                        vortend = vortend+sum(weight(surface_index, :) * time_fact(surface_index, :) * &
                                         (vortg(js+surface_index, is:is+1)  - parcels%vorticity(1, n)))
                    endif

                    if (damping%l_surface_scalars) then
                        ! Note this exponential factor can be different for vorticity/scalars
                        time_fact = one - exp(-damping%scalars_prefactor * strain_mag(js:js+1, is:is+1) * dt)
#ifndef ENABLE_DRY_MODE
                        humtend = humtend + sum(weight(surface_index, :) * time_fact(surface_index, :) * &
                                                  (humg(js+surface_index, is:is+1) - parcels%humidity(n)))
                        dbuoytend = dbuoytend + sum(weight(surface_index, :)  * time_fact(surface_index, :) * &
                                                  (dbuoyg(js+surface_index, is:is+1) - parcels%buoyancy(n)))
#else
                        tbuoytend = tbuoytend + sum(weight(surface_index, :) * time_fact(surface_index, :) * &
                                                        (tbuoyg(js+surface_index, is:is+1) - parcels%buoyancy(n)))
#endif
                    endif
                enddo
                ! Add all the tendencies only at the end
                if (damping%l_vorticity .or. damping%l_surface_vorticity) then
                    parcels%vorticity(1, n) = parcels%vorticity(1, n) + vortend
                endif
                if (damping%l_scalars .or. damping%l_surface_scalars) then
#ifndef ENABLE_DRY_MODE
                    parcels%humidity(n) = parcels%humidity(n) + humtend
                    parcels%buoyancy(n) = parcels%buoyancy(n) + dbuoytend
#else
                    parcels%buoyancy(n) = parcels%buoyancy(n) + tbuoytend
#endif
                endif
            enddo
            !$omp end do
            !$omp end parallel

            call stop_timer(damping_timer)

        end subroutine perturbation_damping_idealised

        !
        ! @pre: the strain must be calculated and the gridded fields updated
        subroutine perturbation_damping_realistic(parcels, dt, l_reuse)
            class(realistic_parcel_type), intent(inout) :: parcels
            double precision, intent(in)  :: dt
            logical, intent(in)           :: l_reuse
            integer                       :: n, p, l, surface_index
            double precision              :: points(2, 2)
            double precision              :: pvol
            ! tendencies need to be summed up between associated 4 points
            ! before modifying the parcel attribute
            double precision              :: vortend
            double precision              :: qvtend
            double precision              :: qltend
            double precision              :: Nltend
            double precision              :: thetatend

            call start_timer(damping_timer)

            
            !$omp parallel default(shared)
            !$omp do private(n, p, l, points, pvol, weight, surface_index) &
            !$omp& private(is, js, weights, vortend, qvtend, qltend, Nltend) &
            !$omp& private(thetatend, time_fact)
            do n = 1, n_parcels
                ! check if only surface damping applies and we are far from surfaces
                ! put in a buffer here as parcels can get stretched in integration
                if(.not.(damping%l_vorticity .or. damping%l_scalars)) then
                    if(parcels%position(2, n) > lower(2) + 2 * dx(2)) then
                    if(parcels%position(2, n) < upper(2) - 2 * dx(2)) then
                        cycle
                    end if
                    end if
                endif

   
                pvol = parcels%volume(n)

                points = get_ellipse_points(parcels%position(:, n), &
                                            pvol, parcels%B(:, n))

                vortend = zero
                thetatend = zero
                if(parcels%is_moist) then
                    qvtend = zero
                    qltend = zero
                endif
                if(parcels%has_droplets) then
                    Nltend = zero
                endif

                ! we have 2 points per ellipse
                do p = 1, 2

                    ! ensure point is within the domain
                    call apply_periodic_bc(points(:, p))

                    ! get interpolation weights and mesh indices
                    call bilinear(points(:, p), is, js, weights)

                    ! loop over grid points which are part of the interpolation
                    ! the weight is halved due to 2 points per ellipse
                    weight = f12 * weights * pvol

                    if (damping%l_vorticity) then
                        ! Note this exponential factor can be different for vorticity/scalars
                        time_fact = one - exp(-damping%vorticity_prefactor * strain_mag(js:js+1, is:is+1) * dt)
                        vortend = vortend + sum(weight * time_fact * (vortg(js:js+1, is:is+1) &
                                       - parcels%vorticity(1, n)))
                    endif

                    if (damping%l_scalars) then
                        time_fact = one - exp(-damping%scalars_prefactor * strain_mag(js:js+1, is:is+1) * dt)
                        thetatend = thetatend + sum(weight * time_fact * (thetag(js:js+1, is:is+1) - parcels%theta(n)))
                        if(parcels%is_moist) then
                            qltend = qltend + sum(weight * time_fact * (qlg(js:js+1, is:is+1) - parcels%ql(n)))
                            qvtend = qvtend + sum(weight * time_fact * (qvg(js:js+1, is:is+1) - parcels%qv(n)))
                         endif
                        if(parcels%has_droplets) then
                            Nltend = Nltend + sum(weight * time_fact * (Nlg(js:js+1, is:is+1) - parcels%Nl(n)))
                         endif
                    endif

                    if (damping%l_surface_vorticity .or. damping%l_surface_scalars) then
                        ! Index to keep track of grid cells right above/below boundary
                        ! This is because the damping only happens at the boundary level
                        ! Consistent with reflection used in parcel_damp
                        if ((js == -1) .or. (js == nz-1)) then
                            surface_index = 1 ! below lower or below upper boundary
                        elseif ((js == 0) .or. (js == nz)) then
                            surface_index = 0 ! above lower or above upper boundary
                        else
                            cycle ! continue loop if not near a surface
                        endif
                    else
                        cycle ! continue loop if no surface damping
                    endif

                    if (damping%l_surface_vorticity) then
                        ! Note this exponential factor can be different for vorticity/scalars
                        time_fact = one - exp(-damping%vorticity_prefactor * strain_mag(js:js+1, is:is+1) * dt)
                        vortend = vortend+sum(weight(surface_index, :) * time_fact(surface_index, :) * &
                                         (vortg(js+surface_index, is:is+1)  - parcels%vorticity(1, n)))
                    endif

                    if (damping%l_surface_scalars) then
                        ! Note this exponential factor can be different for vorticity/scalars
                        time_fact = one - exp(-damping%scalars_prefactor * strain_mag(js:js+1, is:is+1) * dt)
                        thetatend = thetatend + sum(weight(surface_index, :) * time_fact(surface_index, :) * &
                                                  (thetag(js+surface_index, is:is+1) - parcels%theta(n)))
                        if(parcels%is_moist) then
                            qvtend = qvtend + sum(weight(surface_index, :)  * time_fact(surface_index, :) * &
                                                  (qvg(js+surface_index, is:is+1) - parcels%qv(n)))
                            qltend = qltend + sum(weight(surface_index, :)  * time_fact(surface_index, :) * &
                                                  (qlg(js+surface_index, is:is+1) - parcels%ql(n)))
                        endif
                        if(parcels%has_droplets) then
                            Nltend = Nltend + sum(weight(surface_index, :)  * time_fact(surface_index, :) * &
                                                  (Nlg(js+surface_index, is:is+1) - parcels%Nl(n)))
                        endif                        
                    endif
                enddo
                ! Add all the tendencies only at the end
                if (damping%l_vorticity .or. damping%l_surface_vorticity) then
                    parcels%vorticity(1, n) = parcels%vorticity(1, n) + vortend
                endif
                if (damping%l_scalars .or. damping%l_surface_scalars) then
                    parcels%theta(n) = parcels%theta(n) + thetatend
                    if(parcels%is_moist) then
                        parcels%qv(n) = parcels%qv(n) + qvtend
                        parcels%ql(n) = parcels%ql(n) + qltend
                    endif 
                    if(parcels%has_droplets) then
                        parcels%Nl(n) = parcels%Nl(n) + Nltend
                    endif 
                endif
            enddo
            !$omp end do
            !$omp end parallel

            call stop_timer(damping_timer)

        end subroutine perturbation_damping_realistic

end module parcel_damping

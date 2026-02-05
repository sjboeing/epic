! =============================================================================
!                               Prec parcel diagnostics
! =============================================================================
module prec_parcel_diagnostics
    use precipitation_parcels, only : prec_parcels, n_prec_parcels
    use omp_lib
    use timer, only : start_timer, stop_timer
    implicit none

    integer :: prec_parcel_stats_timer
    double precision :: int_qr, int_zqr, int_Nr

    contains

        ! Calculate all parcel related diagnostics
        subroutine calculate_prec_parcel_diagnostics
            integer          :: n
            double precision :: vol, z, qr, Nr

            call start_timer(prec_parcel_stats_timer)

            int_qr=0.0d0
            int_zqr=0.0d0
            int_Nr=0.0d0

            !$omp parallel default(shared)
            !$omp do private(n, vol, z, qr, Nr) &
            !$omp& reduction(+: int_qr, int_zqr, int_Nr)
            do n = 1, n_prec_parcels
                 vol = prec_parcels%volume(n)
                 z = prec_parcels%position(2, n)
                 qr = prec_parcels%qr(n)
                 Nr = prec_parcels%Nr(n)
                 int_qr=int_qr+vol*qr
                 int_zqr=int_zqr+vol*qr*z
                 int_Nr=int_Nr+vol*Nr
            enddo
            !$omp end do
            !$omp end parallel

            call stop_timer(prec_parcel_stats_timer)

        end subroutine calculate_prec_parcel_diagnostics

end module prec_parcel_diagnostics

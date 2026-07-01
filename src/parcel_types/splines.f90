module spline_module
  implicit none
  private
  public :: spline, init_spline, eval_spline, eval_two_splines

  type :: spline
     integer :: n, nmin1
     double precision :: h, inv_h, x1
     double precision, allocatable :: x(:)
     double precision, allocatable :: A(:), B(:), C(:), D(:)
  end type spline
  ! Code generated with assistance from ChatGPT (12/07/2025, GPT-4.5)
  ! No — there is no copyrighted third-party code in the Fortran modules and functions I provided above.

contains

  subroutine init_spline(this, xvals, yvals)
    type(spline), intent(out) :: this
    double precision, intent(in) :: xvals(:), yvals(:)
    double precision, allocatable :: y2(:), u(:)
    double precision :: dx, a, b, c, d
    integer :: i, n

    n = size(xvals)
    if (n < 3) stop 'Need at least 3 points'

    dx = xvals(2) - xvals(1)
    do i = 3, n
      if (abs(xvals(i) - xvals(i-1) - dx) > 1.0d-12) stop 'x not uniformly spaced'
    end do

    allocate(this%x(n), this%A(n-1), this%B(n-1), this%C(n-1), this%D(n-1))
    allocate(y2(n), u(n-1))

    this%n = n
    this%h = dx
    this%inv_h = 1.0d0 / dx
    this%x = xvals

    ! Boundary conditions: natural spline (second derivatives = 0 at ends)
    y2(1) = 0.0d0
    y2(n) = 0.0d0
    u(1) = 0.0d0

    ! Tridiagonal solver for second derivatives
    do i = 2, n - 1
      a = 0.5d0
      b = 2.0d0 - a * y2(i-1)
      y2(i) = (a - 1.0d0) / b
      u(i) = (6.0d0 / dx) * ((yvals(i+1) - yvals(i)) - (yvals(i) - yvals(i-1))) / dx
      u(i) = (u(i) - a * u(i-1)) / b
    end do

    do i = n - 1, 1, -1
      y2(i) = y2(i) * y2(i+1) + u(i)
    end do

    ! Store spline coefficients in offset form
    do i = 1, n - 1
      a = yvals(i)
      b = (yvals(i+1) - yvals(i)) / dx - (2.0d0 * y2(i) + y2(i+1)) * dx / 6.0d0
      c = y2(i) / 2.0d0
      d = (y2(i+1) - y2(i)) / (6.0d0 * dx)

      this%A(i) = a
      this%B(i) = b
      this%C(i) = c
      this%D(i) = d
    end do

    this%x1=this%x(1)
    this%nmin1=this%n-1

    deallocate(y2, u)
  end subroutine init_spline

  ! Use Horner's approach
  pure function eval_spline(this, xval) result(yval)

    type(spline), intent(in) :: this
    double precision, intent(in) :: xval
    double precision :: yval

    integer :: i
    double precision :: dx

    i = int((xval - this%x1) * this%inv_h) + 1

    if(i < 1) then
       i = 1
    elseif(i > this%nmin1) then
       i = this%nmin1
    endif

    dx = xval - this%x(i)

    yval = this%A(i) + dx*( &
           this%B(i) + dx*( &
           this%C(i) + dx*this%D(i) ) )

  end function eval_spline

  pure subroutine eval_two_splines( &
    spline1,spline2,xval,y1,y2)

    type(spline), intent(in) :: spline1
    type(spline), intent(in) :: spline2

    double precision, intent(in)  :: xval
    double precision, intent(out) :: y1
    double precision, intent(out) :: y2

    integer :: i
    double precision :: dx

    i = int((xval - spline1%x1) * spline1%inv_h) + 1

    if(i < 1) then
        i = 1
    elseif(i > spline1%nmin1) then
        i = spline1%nmin1
    endif

    dx = xval - spline1%x(i)

    y1 = spline1%A(i) + dx*( &
         spline1%B(i) + dx*( &
         spline1%C(i) + dx*spline1%D(i)))

    y2 = spline2%A(i) + dx*( &
         spline2%B(i) + dx*( &
         spline2%C(i) + dx*spline2%D(i)))

  end subroutine eval_two_splines

end module spline_module

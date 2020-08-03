! quasirandom.f90 --
!     Generate points in a unit hypercube that are distributed
!     quasi-randomly. The algorithm is based on:
!     "The Unreasonable Effectiveness of Quasirandom Sequences" by Martin Roberts,
!     http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/
!
!     TODO:
!     Support for points on a circle, disk, sphere or ball
!     Integration
!
module quasirandom
    implicit none

    private
    public :: quasirandom_generator

    integer, parameter :: dp = kind(1.0d0)

    ! Implementation still to be done
    integer, parameter, public :: quasi_circle = -1
    integer, parameter, public :: quasi_disk   = -2
    integer, parameter, public :: quasi_sphere = -3
    integer, parameter, public :: quasi_ball   = -4

    type quasirandom_generator
        integer :: dimin    = -1
        integer :: dim      = -1
        integer :: step     = 1
        integer :: stepsize = 1
        real(kind=dp), dimension(:), allocatable :: factor
    contains
        procedure :: init => init_quasi
        procedure :: restart => restart_quasi
        procedure :: single_next => single_next_quasi
        procedure :: double_next => double_next_quasi
        generic   :: next => single_next, double_next
    end type quasirandom_generator

contains

! init_quasi --
!     Set the dimension and the start point for the quasi-random generator
!
! Arguments:
!     dim                 Dimension of the space for which to generate the points
!     stepsize            (Optional) use a different step than 1 to select the next point
!
subroutine init_quasi( this, dimin, stepsize )
    class(quasirandom_generator)  :: this
    integer, intent(in)           :: dimin
    integer, intent(in), optional :: stepsize

    integer                       :: i, n, dim
    real(kind=dp)                 :: f

    if ( present(stepsize) ) then
        this%stepsize = stepsize
    else
        this%stepsize = 1
    endif

    if ( dimin == 0 .or. dimin < quasi_ball ) then
        write(*,*) 'Quasirandom numbers: Dimension must be positive!'
        return
    endif

    select case ( dimin )
        case( quasi_circle )
            dim = 1
        case( quasi_disk )
            dim = 2
        case( quasi_sphere )
            dim = 2
        case( quasi_ball )
            dim = 3
        case default
            dim = dimin
    end select

    this%dimin = dimin
    this%dim   = dim

    if ( allocated(this%factor) ) then
        deallocate( this%factor )
    endif
    allocate( this%factor(dim) )

    !
    ! Solve the equation x**n - x + 1 = 0
    !
    n = dim + 1
    f = 1.0_dp
    do i = 1,10
        f = f - (f**n - f - 1.0_dp) / (n*f**(n-1)-1.0_dp)
    enddo

    this%factor(1) = 1.0_dp / f

    do i = 2,dim
        this%factor(i) = this%factor(i-1) / f
    enddo
end subroutine init_quasi

! restart_quasi --
!     Restart the generator
!
! Arguments:
!     start               (Optional) Starting point (simply the n'th point that would be generated)
!     stepsize            (Optional) use a different step than 1 to select the next point
!
subroutine restart_quasi( this, start, stepsize )
    class(quasirandom_generator)  :: this
    integer, intent(in), optional :: start
    integer, intent(in), optional :: stepsize

    if ( present(start) ) then
        this%step = start
    else
        this%step = 1
    endif

    if ( present(stepsize) ) then
        this%stepsize = stepsize
    else
        this%stepsize = 1
    endif
end subroutine restart_quasi

! single_next_quasi --
!     Generate the next point (single-precision reals)
!
! Arguments:
!     coords              Array for holding the new coordinates
!
subroutine single_next_quasi( this, coords )
    class(quasirandom_generator)    :: this
    real, dimension(:), intent(out) :: coords

    coords(1:this%dim) = modulo( this%step * this%factor, 1.0_dp )

    if ( this%dimin < 0 ) then
        call single_transform( this%dimin, coords )
    endif

    this%step = this%step + this%stepsize

end subroutine single_next_quasi

! double_next_quasi --
!     Generate the next point (double-precision reals)
!
! Arguments:
!     coords              Array for holding the new coordinates
!
subroutine double_next_quasi( this, coords )
    class(quasirandom_generator)             :: this
    real(kind=dp), dimension(:), intent(out) :: coords

    coords(1:this%dim) = modulo( this%step * this%factor, 1.0_dp )

    if ( this%dimin < 0 ) then
        call double_transform( this%dimin, coords )
    endif

    this%step = this%step + this%stepsize

end subroutine double_next_quasi

! single_transform --
!     Transform the coordinates to the required shape
!
subroutine single_transform( type, coords )
    integer, intent(in)               :: type
    real, dimension(:), intent(inout) :: coords

    real, parameter                   :: twopi = 2.0 * acos(-1.0)
    real                              :: x, y, z, r, theta, phi

    select case( type )
        case( quasi_circle )
            x = cos( twopi * coords(1) )
            y = sin( twopi * coords(1) )
            coords = [x, y]

        case( quasi_disk )
            x = sqrt(coords(1)) * cos( twopi * coords(2) )
            y = sqrt(coords(1)) * sin( twopi * coords(2) )
            coords = [x, y]

        case( quasi_sphere )
            theta = asin( 2.0 * coords(1) - 1.0 )
            phi   = twopi * coords(2)
            x     = cos(phi) * cos(theta)
            y     = sin(phi) * cos(theta)
            z     = 2.0 * coords(1) - 1.0 ! sin(theta)
            coords = [ x, y, z ]

        case( quasi_ball )
            theta = asin( 2.0 * coords(1) - 1.0 )
            phi   = twopi * coords(2)
            r     = coords(3) ** (1.0/3.0)
            x     = cos(phi) * cos(theta)
            y     = sin(phi) * cos(theta)
            z     = 2.0 * coords(1) - 1.0 ! sin(theta)
            coords = r * [ x, y, z ]
     end select
end subroutine single_transform

! double_transform --
!     Transform the coordinates to the required shape (double precision)
!
subroutine double_transform( type, coords )
    integer, intent(in)                   :: type
    real(dp), dimension(:), intent(inout) :: coords

    real(dp), parameter                   :: twopi = 2.0_dp * acos(-1.0_dp)
    real(dp)                              :: x, y, z, r, theta, phi

    select case( type )
        case( quasi_circle )
            x = cos( twopi * coords(1) )
            y = sin( twopi * coords(1) )
            coords = [x, y]

        case( quasi_disk )
            x = sqrt(coords(1)) * cos( twopi * coords(2) )
            y = sqrt(coords(1)) * sin( twopi * coords(2) )
            coords = [x, y]

        case( quasi_sphere )
            theta = asin( 2.0 * coords(1) - 1.0 )
            phi   = twopi * coords(2)
            x     = cos(phi) * cos(theta)
            y     = sin(phi) * cos(theta)
            z     = 2.0 * coords(1) - 1.0 ! sin(theta)
            coords = [ x, y, z ]

        case( quasi_ball )
            theta = asin( 2.0 * coords(1) - 1.0 )
            phi   = twopi * coords(2)
            r     = coords(3) ** (1.0/3.0)
            x     = cos(phi) * cos(theta)
            y     = sin(phi) * cos(theta)
            z     = 2.0 * coords(1) - 1.0 ! sin(theta)
            coords = r * [ x, y, z ]
     end select
end subroutine double_transform

end module quasirandom

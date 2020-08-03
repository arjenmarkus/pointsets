! test_quasirandom.f90 --
!     Simple test of the quasi-random points generator
!
!     The algorithm is based on:
!     "The Unreasonable Effectiveness of Quasirandom Sequences" by Martin Roberts,
!     http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/
!
!     TODO:
!     Support for points on a circle, disk, sphere or ball
!     Integration
!

! test_quasirandom
!     Simple test program
!
program test_quasirandom
    use quasirandom

    implicit none

    integer, parameter :: dp = kind(1.0d0)

    type(quasirandom_generator) :: q, qcircle, qdisk, qsphere, qball

    real, dimension(1)          :: coords
    real, dimension(2)          :: coords2
    real, dimension(3)          :: coords3
    real                        :: radius2

    real(kind=dp), dimension(1) :: dcoords
    real(kind=dp), dimension(2) :: dcoords2
    real(kind=dp), dimension(3) :: dcoords3
    real(kind=dp)               :: dradius2

    integer                     :: i, points_in, points_off, points_x, points_y, points_z

    write( *, * ) 'SINGLE-PRECISION VERSION'
    write( *, * ) '------------------------'

    call q%init( 1 )

    write( *, * ) 'Quasirandom points on the interval [0,1]:'
    do i = 1,10
        call q%next( coords )
        write( *, * ) i, coords
    enddo

    call q%init( 3 )

    write( *, * ) 'Quasirandom points on the unit cube [0,1]^3:'
    do i = 1,10
        call q%next( coords3 )
        write( *, * ) i, coords3
    enddo

    !
    ! Test the special forms:
    ! - For circle and sphere:
    !   - the length of the vector should be one
    !   - the number of points with x > 0.7 should be roughly equal to
    !     the number of points with y > 0.7
    ! - For disk and ball:
    !   - the number of points with distance to the origin less than 0.5
    !     should roughly be 1/4 and 1/8 of the total
    !   - the number of points with x > 0.7 should be roughly equal to
    !     the number of points with y > 0.7 and to the number of points
    !     with z > 0.7
    !
    write( *, * ) 'Quasirandom points on a circle:'

    points_off = 0
    points_x   = 0
    points_y   = 0
    call qcircle%init( quasi_circle )
    do i = 1,1000
        call qcircle%next( coords2 )
        radius2 = sum( coords2 ** 2 )
        if ( abs(radius2 - 1.0) > 1.0e-6 ) then
            points_off = points_off + 1
        endif
        if ( coords2(1) > 0.7 ) then
            points_x = points_x  + 1
        endif
        if ( coords2(2) > 0.7 ) then
            points_y = points_y  + 1
        endif
    enddo

    write( *, * ) 'Number of points off the circle:  ', points_off
    write( *, * ) 'Number of points right of x = 0.7:', points_x
    write( *, * ) 'Number of points above of y = 0.7:', points_y

    write( *, * ) 'Quasirandom points on a disk:'
    points_in  = 0
    points_x   = 0
    points_y   = 0
    call qdisk%init( quasi_disk )
    do i = 1,1000
        call qdisk%next( coords2 )
        radius2 = sum( coords2 ** 2 )
        if ( radius2 < 0.25 ) then
            points_in = points_in + 1
        endif
        if ( coords2(1) > 0.7 ) then
            points_x = points_x  + 1
        endif
        if ( coords2(2) > 0.7 ) then
            points_y = points_y  + 1
        endif
    enddo

    write( *, * ) 'Number of points inner part of disk:  ', points_in, ' - expected: ', 1000/4
    write( *, * ) 'Number of points right of x = 0.7:', points_x
    write( *, * ) 'Number of points above of y = 0.7:', points_y

    !
    ! Sphere and ball
    !
    write( *, * ) 'Quasirandom points on a sphere:'

    points_off = 0
    points_x   = 0
    points_y   = 0
    points_z   = 0
    call qsphere%init( quasi_sphere )
    do i = 1,1000
        call qsphere%next( coords3 )
        radius2 = sum( coords3 ** 2 )
        if ( abs(radius2 - 1.0) > 1.0e-6 ) then
            points_off = points_off + 1
        endif
        if ( coords3(1) > 0.7 ) then
            points_x = points_x  + 1
        endif
        if ( coords3(2) > 0.7 ) then
            points_y = points_y  + 1
        endif
        if ( coords3(3) > 0.7 ) then
            points_z = points_z  + 1
        endif
    enddo

    write( *, * ) 'Number of points off the sphere:  ', points_off
    write( *, * ) 'Number of points right of x = 0.7:', points_x
    write( *, * ) 'Number of points above of y = 0.7:', points_y
    write( *, * ) 'Number of points above of z = 0.7:', points_z

    write( *, * ) 'Quasirandom points on a ball:'
    points_in  = 0
    points_x   = 0
    points_y   = 0
    points_z   = 0
    call qball%init( quasi_ball )
    do i = 1,1000
        call qball%next( coords3 )
        radius2 = sum( coords3 ** 2 )
        if ( radius2 < 0.25 ) then
            points_in = points_in + 1
        endif
        if ( coords3(1) > 0.7 ) then
            points_x = points_x  + 1
        endif
        if ( coords3(2) > 0.7 ) then
            points_y = points_y  + 1
        endif
        if ( coords3(3) > 0.7 ) then
            points_z = points_z  + 1
        endif
    enddo

    write( *, * ) 'Number of points inner part of ball:  ', points_in, ' - expected: ', 1000/8
    write( *, * ) 'Number of points right of x = 0.7:', points_x
    write( *, * ) 'Number of points above of y = 0.7:', points_y
    write( *, * ) 'Number of points above of z = 0.7:', points_z

    !
    ! "Just" a copy of the previous code
    !
    write( *, * ) ''
    write( *, * ) 'DOUBLE-PRECISION VERSION'
    write( *, * ) '------------------------'

    call q%init( 1 )

    write( *, * ) 'Quasirandom points on the interval [0,1]:'
    do i = 1,10
        call q%next( dcoords )
        write( *, * ) i, dcoords
    enddo

    call q%init( 3 )

    write( *, * ) 'Quasirandom points on the unit cube [0,1]^3:'
    do i = 1,10
        call q%next( dcoords3 )
        write( *, * ) i, dcoords3
    enddo

    !
    ! Test the special forms:
    ! - For circle and sphere:
    !   - the length of the vector should be one
    !   - the number of points with x > 0.7 should be roughly equal to
    !     the number of points with y > 0.7
    ! - For disk and ball:
    !   - the number of points with distance to the origin less than 0.5
    !     should roughly be 1/4 and 1/8 of the total
    !   - the number of points with x > 0.7 should be roughly equal to
    !     the number of points with y > 0.7 and to the number of points
    !     with z > 0.7
    !
    write( *, * ) 'Quasirandom points on a circle:'

    points_off = 0
    points_x   = 0
    points_y   = 0
    call qcircle%init( quasi_circle )
    do i = 1,1000
        call qcircle%next( dcoords2 )
        dradius2 = sum( dcoords2 ** 2 )
        if ( abs(dradius2 - 1.0_dp) > 1.0e-12_dp ) then
            points_off = points_off + 1
        endif
        if ( dcoords2(1) > 0.7_dp ) then
            points_x = points_x  + 1
        endif
        if ( dcoords2(2) > 0.7_dp ) then
            points_y = points_y  + 1
        endif
    enddo

    write( *, * ) 'Number of points off the circle:  ', points_off
    write( *, * ) 'Number of points right of x = 0.7:', points_x
    write( *, * ) 'Number of points above of y = 0.7:', points_y

    write( *, * ) 'Quasirandom points on a disk:'
    points_in  = 0
    points_x   = 0
    points_y   = 0
    call qdisk%init( quasi_disk )
    do i = 1,1000
        call qdisk%next( dcoords2 )
        dradius2 = sum( dcoords2 ** 2 )
        if ( dradius2 < 0.25_dp ) then
            points_in = points_in + 1
        endif
        if ( dcoords2(1) > 0.7_dp ) then
            points_x = points_x  + 1
        endif
        if ( dcoords2(2) > 0.7_dp ) then
            points_y = points_y  + 1
        endif
    enddo

    write( *, * ) 'Number of points inner part of disk:  ', points_in, ' - expected: ', 1000/4
    write( *, * ) 'Number of points right of x = 0.7:', points_x
    write( *, * ) 'Number of points above of y = 0.7:', points_y

    !
    ! Sphere and ball
    !
    write( *, * ) 'Quasirandom points on a sphere:'

    points_off = 0
    points_x   = 0
    points_y   = 0
    points_z   = 0
    call qsphere%init( quasi_sphere )
    do i = 1,1000
        call qsphere%next( dcoords3 )
        dradius2 = sum( dcoords3 ** 2 )
        if ( abs(dradius2 - 1.0_dp) > 1.0e-12_dp ) then
            points_off = points_off + 1
        endif
        if ( dcoords3(1) > 0.7_dp ) then
            points_x = points_x  + 1
        endif
        if ( dcoords3(2) > 0.7_dp ) then
            points_y = points_y  + 1
        endif
        if ( dcoords3(3) > 0.7_dp ) then
            points_z = points_z  + 1
        endif
    enddo

    write( *, * ) 'Number of points off the sphere:  ', points_off
    write( *, * ) 'Number of points right of x = 0.7:', points_x
    write( *, * ) 'Number of points above of y = 0.7:', points_y
    write( *, * ) 'Number of points above of z = 0.7:', points_z

    write( *, * ) 'Quasirandom points on a ball:'
    points_in  = 0
    points_x   = 0
    points_y   = 0
    points_z   = 0
    call qball%init( quasi_ball )
    do i = 1,1000
        call qball%next( dcoords3 )
        dradius2 = sum( dcoords3 ** 2 )
        if ( dradius2 < 0.25_dp ) then
            points_in = points_in + 1
        endif
        if ( dcoords3(1) > 0.7_dp ) then
            points_x = points_x  + 1
        endif
        if ( dcoords3(2) > 0.7_dp ) then
            points_y = points_y  + 1
        endif
        if ( dcoords3(3) > 0.7_dp ) then
            points_z = points_z  + 1
        endif
    enddo

    write( *, * ) 'Number of points inner part of ball:  ', points_in, ' - expected: ', 1000/8
    write( *, * ) 'Number of points right of x = 0.7:', points_x
    write( *, * ) 'Number of points above of y = 0.7:', points_y
    write( *, * ) 'Number of points above of z = 0.7:', points_z
end program test_quasirandom

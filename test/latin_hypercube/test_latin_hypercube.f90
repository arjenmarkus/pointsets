! test_latin_hypercube.f90 --
!     Test program for Latin hypercube sampling
!
program test_latin_hypercube
    use latin_hypercube

    implicit none

    real(wp) :: coords(5,10)

    call generate_points( coords )

    write(*,'(a)')       'Sample points:'
    write(*,'(5g12.4)' ) coords

    write(*,'(a)')       'Sample points - check:'
    write(*,'(5g12.0)' ) int(coords*size(coords,2))
    write(*,'(a)')       'Note: the coordinates per column should be unique'


end program test_latin_hypercube

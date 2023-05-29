! latin_hypercube.f90 --
!     Use the Latin hybercube method to generate points in a hypercube
!     that are uniformly distributed.
!
!     More information:
!     https://en.wikipedia.org/wiki/Latin_hypercube_sampling
!
!     Partly inspired by code from John Burkardt.
!
module latin_hypercube
    use random_permutations
    use select_precision

    implicit none

contains

! generate_points --
!     Generate n points in an m-dimensional space according to
!     the rules of Latin hypercube sampling.
!
! Arguments:
!     coords              Coordinates of the selected points
!
! Note:
!     The dimensions of the two-dimensional array are used to
!     define the number of space dimensions and the number
!     of points in that space.
!     The first dimension is the number of space dimensions,
!     the second the number of points.
!
!
subroutine generate_points( coords )
    real(wp), intent(out) :: coords(:,:)

    integer               :: perm(size(coords,2))
    integer               :: i, j

    call random_number( coords )

    perm = [(i, i = 1,size(perm))]

    do i = 1,size(coords,1)
        call randomly_permute( perm )

        do j = 1,size(coords,2)
            coords(i,j) = ( perm(j) - 1 + coords(i,j) ) / size(coords,2)
        enddo
    enddo
end subroutine generate_points

end module latin_hypercube

module initialize

    use precision
    use read_parameters

    implicit none

    type, public :: drop_type
        real(pr), dimension(:), allocatable :: position  ! 3D position vector
        real(pr), dimension(:), allocatable :: velocity  ! 3D velocity vector
        real(pr) :: radius         ! Particle radius
        real(pr) :: mass           ! Particle mass
        real(pr) :: temperature    ! Particle temperature
    end type drop_type
    
contains

    subroutine init_drops(param, drops)

        !!This subroutine allow to initialize the drops field

        !In
        type(parameters), intent(in) :: param
        !Out
        type(drop_type), dimension(param%nb_drops), intent(inout) :: drops

        !Local
        integer :: p,i
        real(pr) :: radius
        real(pr) :: pos_rand(param%dim)  ! Array to store random positions for each particle
        real(pr) :: radius_rand


        do p=1,param%nb_drops 

            allocate(drops(p)%position(param%dim))
            allocate(drops(p)%velocity(param%dim))

            ! Generate a random position for each component of the position vector
            call random_number(pos_rand)
            drops(p)%position = pos_rand  ! Random positions in the range [0, 1)

            ! Set the same velocity for each particle
            do i=1,param%dim
                drops(p)%velocity(i) = 0.0
            end do

            ! Generate a random radius within the range [1.0e-3, 1.0e-2]
            call random_number(radius_rand) !Tirer selon N(r)
            drops(p)%radius = 1.0e-3 + radius_rand * (1.0e-3 - 1.0e-6)

            ! Compute the particle mass based on its radius (assuming spherical volume)
            radius = drops(p)%radius
            drops(p)%mass = 4./3_pr * pi * radius**3 * param%rho_w

            ! Set the same temperature for each particle
            drops(p)%temperature = param%T_air

        end do

    end subroutine init_drops


    
end module initialize
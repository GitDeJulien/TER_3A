module initialize

    use precision

    implicit none

    type, public :: drop_type
        real(pr) :: position(dim)    ! 3D position vector
        real(pr) :: velocity(dim)    ! 3D velocity vector
        real(pr) :: radius         ! Particle radius
        real(pr) :: mass           ! Particle mass
        real(pr) :: temperature    ! Particle temperature
    end type drop_type

    type, public :: NM_type
        integer                :: sch
        real(pr)               :: tmax, cfl, time_step
    end type NM_type
    
contains

    subroutine init_drops(drops)

        !!This subroutine allow to initialize the probleme

        !Out
        type(drop_type), intent(out) :: drops(nb_drops)

        !Local
        integer :: p,i
        real(pr) :: radius, volume
        real(pr) :: pos_rand(dim)  ! Array to store random positions for each particle
        real(pr) :: radius_rand


        do p=1,nb_drops 

            ! Generate a random position for each component of the position vector
            call random_number(pos_rand)
            drops(p)%position = pos_rand  ! Random positions in the range [0, 1)

            ! Set the same velocity for each particle
            do i=1,dim
                drops(p)%velocity(i) = 5.0_pr
            end do

            ! Generate a random radius within the range [1.0e-3, 1.0e-2]
            call random_number(radius_rand)
            drops(p)%radius = 1.0e-3 + radius_rand * (1.0e-2 - 1.0e-3)

            ! Compute the particle mass based on its radius (assuming spherical volume)
            radius = drops(p)%radius
            volume = (4.0_pr / 3.0_pr) * pi * radius**3
            drops(p)%mass = volume * rho_eau

            ! Set the same temperature for each particle
            drops(p)%temperature = 20.0_pr + 273.0_pr

        end do

    end subroutine init_drops


    subroutine init_model(NM)

        !Out
        type(NM_type), intent(out) :: NM

        NM%sch = 0
        NM%tmax = 1.0
        NM%cfl = 0.5

    end subroutine init_model


    
end module initialize
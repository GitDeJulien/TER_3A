module precision

    implicit none

    integer, parameter  :: pr = 8
    real(pr), parameter :: pi = acos(-1._pr)
    real(pr), parameter :: gravity = 9.80665
    real(pr), parameter :: rho_eau = 1020.0
    real(pr), parameter :: dyn_viscosity = 8.90E10-4
    real(pr), parameter :: air_velocity = 5.0
    integer, parameter  :: nb_drops = 10
    integer, parameter  :: dim = 2
    
end module precision
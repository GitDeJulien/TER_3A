module read_parameters
    
    use precision

    implicit none

    type, public :: parameters

        ! In the data file
        real(pr) :: U_air
        real(pr) :: T_air_celcius
        real(pr) :: T_p_0_celcius
        real(pr) :: Q_RH
        real(pr) :: gravity
        real(pr) :: M_w
        real(pr) :: M_s
        real(pr) :: M_a
        real(pr) :: R_g
        real(pr) :: p_0
        real(pr) :: Salinity_w
        real(pr) :: Salinity_p
        real(pr) :: r_p_0
        real(pr) :: Delta_v
        real(pr) :: alpha_c
        real(pr) :: Delta_T
        real(pr) :: alpha_T
        real(pr) :: I

        ! Simulation parameters
        integer :: cas
        integer :: dim
        integer :: nb_drops
        integer :: Nmax_iter
        integer :: sch
        real(pr) :: cfl
        real(pr) :: time_step

        ! Computed from parameters in data file
        real(pr) :: k_a
        real(pr) :: T_air
        real(pr) :: T_p_0
        real(pr) :: qs0
        real(pr) :: rs0
        real(pr) :: r0
        real(pr) :: q0
        real(pr) :: Tv0
        real(pr) :: rho_air
        real(pr) :: nu_air
        real(pr) :: mu_air
        real(pr) :: p_air
        real(pr) :: pv_sat_T_air
        real(pr) :: L_v
        real(pr) :: D_v
        real(pr) :: Gamma_p
        real(pr) :: Phi_s
        real(pr) :: rho_w
        real(pr) :: rho_p
        real(pr) :: ms
        real(pr) :: c_p_s
        real(pr) :: c_p_air_sec
        real(pr) :: rho_air_sec

    end type parameters


contains

    subroutine read_param(param, filename)
        
        !! This subroutine read the parameters of the parameters.dat file in the
        !! input directory and attribut the value with varaibles with the same name

        type(parameters), intent(inout) :: param
        character(len=*), intent(in) :: filename

        character(len=100) :: line, str_value, key
        integer :: ios, equal_pos


        ! Open the parameters.dat file
        open(unit=10, file="../input/parameters.dat", status="old", action="read", iostat=ios)

        if (ios /= 0) then
            print *, "Error: Could not open file", filename
            return
        end if

        ! Read each parameter from the file
        do
            read(10, '(A)', iostat=ios) line
            if (ios /= 0) exit

            ! Skip comment and header lines
            if (trim(line) == "" .or. line(1:1) == "#" .or. line(1:1) == "[") cycle

            ! Find the position of '=' in the line
            equal_pos = index(line, "=")

            ! Check if '=' was found
            if (equal_pos > 0) then

                str_value = trim(adjustl(line(equal_pos+1:)))
                key = trim(adjustl(line(1:equal_pos-2)))


                if (key == "U_air") then 
                    read(str_value, *, iostat=ios) param%U_air
                else if (key == "T_air_celcius")then
                    read(str_value, *, iostat=ios) param%T_air_celcius
                else if (key == "T_p_0_celcius") then
                    read(str_value, *, iostat=ios) param%T_p_0_celcius
                else if (key == "Q_RH") then
                    read(str_value, *, iostat=ios) param%Q_RH
                else if (key == "g") then
                    read(str_value, *, iostat=ios) param%gravity
                else if (key == "M_w") then
                    read(str_value, *, iostat=ios) param%M_w
                else if (key == "M_s") then
                    read(str_value, *, iostat=ios) param%M_s
                else if (key == "M_a") then
                    read(str_value, *, iostat=ios) param%M_a
                else if (key == "R_g") then
                    read(str_value, *, iostat=ios) param%R_g
                else if (key == "p_0") then
                    read(str_value, *, iostat=ios) param%p_0
                else if (key == "Salinity_w") then
                    read(str_value, *, iostat=ios) param%Salinity_w
                else if (key == "Salinity_p") then
                    read(str_value, *, iostat=ios) param%Salinity_p
                else if (key == "r_p_0") then
                    read(str_value, *, iostat=ios) param%r_p_0
                else if (key == "Delta_v") then
                    read(str_value, *, iostat=ios) param%Delta_V
                else if (key == "alpha_c")  then
                    read(str_value, *, iostat=ios) param%alpha_c
                else if (key == "Delta_T") then
                    read(str_value, *, iostat=ios) param%Delta_T
                else if (key == "alpha_T") then
                    read(str_value, *, iostat=ios) param%alpha_T
                else if (key == "I") then
                    read(str_value, *, iostat=ios) param%I
                else if (key == "cas") then
                    read(str_value, *, iostat=ios) param%cas
                else if (key == "dim") then
                    read(str_value, *, iostat=ios) param%dim
                else if (key == "nb_drops") then
                    read(str_value, *, iostat=ios) param%nb_drops
                else if (key == "Nmax_iter") then
                    read(str_value, *, iostat=ios) param%Nmax_iter
                else if (key == "key_scheme") then
                    read(str_value, *, iostat=ios) param%sch
                else if (key == "cfl") then
                    read(str_value, *, iostat=ios) param%cfl
                else if (key == "time_step") then
                    read(str_value, *, iostat=ios) param%time_step
                
                else 
                    print*, "parameter not found in line:", line

                end if

            end if

        end do

        ! Close the file
        close(10)

        ! Computed parameters
        param%k_a = 0.02411*(1+3.309e-3*(param%T_p_0_celcius) -1.441e-6*(param%T_p_0_celcius)**2)
        param%T_air = param%T_air_celcius + 273.15
        param%T_p_0 = param%T_p_0_celcius + 273.15
        param%qs0 = 0.62197*611.2/param%p_0*EXP((17.67*(param%T_p_0_celcius))/(param%T_p_0-29.66))
        param%rs0 = param%qs0/(1-param%qs0)
        param%r0 = 98.0/100.0*param%rs0
        param%q0 = param%r0/(1+param%r0)
        param%Tv0 = param%T_p_0*(1+0.6078*param%q0)
        param%rho_air = 1.2929*273.15/param%Tv0
        param%nu_air = 1.33e-5 +(0.0084*(param%Tv0-273.15))*1e-5
        param%mu_air = param%nu_air*param%rho_air
        param%p_air = param%p_0-param%rho_air*param%gravity*10.0
        param%pv_sat_T_air = 2.33e3
        param%L_v = (25.00895-0.02274*(param%T_air-273.15))*1e5
        param%D_v = 2.11e-5*(param%T_air-273.15)**1.94*(1013.25/(param%p_air/100.0))
        param%Gamma_p = (75.63 - 0.144*param%T_p_0_celcius + 0.221*param%Salinity_w)*1e-3
        param%Phi_s = 0.91154614+1.7317496707e-4*param%Salinity_p+4.7616058412e-6*param%Salinity_p**2-&
        & 9.2541509027e-9*param%Salinity_p**3+7.3475024678e-12*param%Salinity_p**4
        param%rho_w = (999.842594+6.793952e-2*param%T_p_0_celcius-9.095290e-3*param%T_p_0_celcius**2+&
        & 1.001685e-4*param%T_p_0_celcius**3-1.120083e-6*param%T_p_0_celcius**4+6.536332e-9*param%T_p_0_celcius**5)
        param%rho_p = param%rho_w + param%Salinity_p*(0.824493-4.0899e-3*param%T_p_0_celcius + &
        & 7.6438e-5*param%T_p_0_celcius**2-8.2467e-7*param%T_p_0_celcius**3+5.3875e-9*param%T_p_0_celcius**4)+&
        & param%Salinity_p**(3./2.)*(-5.72466e-3+1.0227e-4*param%T_p_0_celcius&
        & -1.6546e-6*param%T_p_0_celcius**2) +4.8314e-4*param%Salinity_p**2
        param%ms = (4.0/3.0)*3.14159265358979323846264338327950288*param%r_p_0**3*&
        & param%rho_p*param%Salinity_p/1000.0
        param%c_p_s = 4217.4 -3.720283*(param%T_air-273.15)+0.1412855*(param%T_air-273.15)**2&
        & -2.654387e-3*(param%T_air-273.15)**3 +2.093236e-5*(param%T_air-273.15)**4 +&
        & param%Salinity_p*(-7.6444+0.107276*(param%T_air-273.15)-1.3839e-3*param%T_p_0_celcius**2)+&
        & param%Salinity_p**(3./2)*(0.17709-4.0772e-3*(param%T_air-273.15)+5.3539e-5*param%T_p_0_celcius**2)
        param%c_p_air_sec = 1.9327e-10*param%T_p_0**4-7.9999e-7*param%T_p_0**3+&
        & 1.1407e-3*param%T_p_0**2-4.4890e-1*param%T_p_0+1.0575e+3
        param%rho_air_sec = 1.2929*273.13/param%T_air


    end subroutine read_param

    subroutine display_parameters(param)

        !This subroutine allow to print all the parameter of the simulation

        type(parameters) :: param

        print*, "Phyisical parameters: "
        print*, "------------------------------------"
        print*, "U_air        :", param%U_air
        print*, "T_air_celcius:", param%T_air_celcius
        print*, "T_p_0_celcius:", param%T_p_0_celcius
        print*, "Q_RH         :", param%Q_RH
        print*, "g            :", param%gravity
        print*, "M_w          :", param%M_w
        print*, "M_s          :", param%M_s
        print*, "M_a          :", param%M_a
        print*, "R_g          :", param%R_g
        print*, "p_0          :", param%p_0
        print*, "Salinity_w   :", param%Salinity_w
        print*, "Salinity_p   :", param%Salinity_p
        print*, "r_p_0        :", param%r_p_0
        print*, "Delta_v      :", param%Delta_v
        print*, "alpha_c      :", param%alpha_c
        print*, "Delta_T      :", param%Delta_T
        print*, "alpha_T      :", param%alpha_T
        print*, "I            :", param%I
        print*, "k_a          :", param%k_a
        print*, "T_air        :", param%T_air
        print*, "T_p_0        :", param%T_p_0
        print*, "qs0          :", param%qs0
        print*, "rs0          :", param%rs0
        print*, "r0           :", param%r0
        print*, "q0           :", param%q0
        print*, "Tv0          :", param%Tv0
        print*, "rho_air      :", param%rho_air
        print*, "nu_air       :", param%nu_air
        print*, "mu_air       :", param%mu_air
        print*, "p_air        :", param%p_air
        print*, "pv_sat(T_air):", param%pv_sat_T_air
        print*, "L_v          :", param%L_v
        print*, "D_v          :", param%D_v
        print*, "Gamma_p      :", param%Gamma_p
        print*, "Phi_s        :", param%Phi_s
        print*, "rhow_w       :", param%rho_w
        print*, "rho_p        :", param%rho_p
        print*, "ms           :", param%ms
        print*, "c_p_s        :", param%c_p_s
        print*, "c_p_air_sec  :", param%c_p_air_sec
        print*, "rho_air_sec  :", param%rho_air_sec
        print*, "..."
        print*, "Parameters of the simulation: "
        print*, "------------------------------------"
        print*, "cas       :", param%cas
        print*, "dim       :", param%dim
        print*, "nb_drops  :", param%nb_drops
        print*, "Nmax_iter :", param%Nmax_iter
        print*, "key_scheme:", param%sch
        print*, "cfl       :", param%cfl
        print*, "time_step :", param%time_step

    end subroutine display_parameters

end module read_parameters

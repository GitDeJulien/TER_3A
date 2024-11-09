module schemes

    use precision
    use initialize

    implicit none
    
contains

    function advance(NM, drops_n) result(drops_np1)

        !! This function computing Xp_(n+1) at the time t_(n+1)
        !! taking the vector Xp_n=(xp_n, vp_n, rp_n, mp_n, Tp_n)
        !! at the time t_n.
        !! It can use different time scheme : 
        !! sch = 0 : Euler explicite


        !In
        type(drop_type), intent(in)  :: drops_n(nb_drops)

        !Inout
        type(NM_type), intent(inout)  :: NM
        
        !Out
        type(drop_type)              :: drops_np1(nb_drops)


        !Local
        integer  :: p
        real(pr) :: dt

        Select case(NM%sch)

        case(0)
            

            do p=1,nb_drops  
                
                !Time step comutation for each particles at each time step
                !! prenre taup = taué
                NM%time_step = NM%cfl * drops_n(p)%mass / (6*pi*drops_n(p)%radius*dyn_viscosity)
                dt = NM%time_step
                print*, 'dt:', dt
                !NM%time_step = 0.01

                !EE Schemes for the 2 first equations and nothing on the remainings
                drops_np1(p)%position = drops_n(p)%position + dt*drops_n(p)%velocity
                drops_np1(p)%velocity = drops_n(p)%velocity + dt*(6*pi*drops_n(p)%radius*dyn_viscosity/drops_n(p)%mass&
                *(air_velocity-drops_n(p)%velocity) + gravity)
                drops_np1(p)%radius = drops_n(p)%radius
                drops_np1(p)%mass = drops_n(p)%mass
                drops_np1(p)%temperature = drops_n(p)%temperature 

            end do

        end Select


    end function advance


    function adimentionalization(drops_n, drops_init_n) result(drops_n_ad)

        !In
        type(drop_type), intent(in)  :: drops_n(nb_drops)
        type(drop_type), intent(in)  :: drops_init_n(nb_drops)

        !Out
        type(drop_type)              :: drops_n_ad(nb_drops)

        !Local
        integer :: p


        do p=1,nb_drops

            drops_n_ad(p)%position = drops_n(p)%position / drops_init_n(p)%position
            drops_n_ad(p)%velocity = drops_n(p)%velocity / drops_init_n(p)%velocity
            drops_n_ad(p)%radius = drops_n(p)%radius / drops_init_n(p)%radius
            drops_n_ad(p)%mass = drops_n(p)%mass / drops_init_n(p)%mass
            drops_n_ad(p)%temperature = drops_n(p)%temperature / drops_init_n(p)%temperature

        end do

    end function



    subroutine out(NM,n,tn,drops)

        type(NM_type), intent(in)    :: NM
        type(drop_type), intent(in)  :: drops(nb_drops)
        integer, intent(in)          :: n
        real(pr), intent(in)         :: tn


        character(len=40) :: cn
        
        integer :: p

        write(cn,'(I4.1)') n

        if (NM%sch == 0) then
            
            open(unit=10,file='../output/model_0/sol_'//trim(adjustl(cn))//'.dat', action='write')

        else 

            print*, 'Error: No such numerical model code exist'

        end if

        write(10,*) '## tn = ', tn
        write(10,*) '## x(m),  y(m),  z(m),  vx(m/s),  vy(m/s),  vz(m/s),  r(m),  m(kg),  T(K)'
        do p=1,nb_drops
            write(10,*) drops(p)%position, drops(p)%velocity, drops(p)%radius, drops(p)%mass, drops(p)%temperature
        end do


        close(10)

    end subroutine out


    
end module schemes
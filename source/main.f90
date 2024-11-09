program embruns

    use precision
    use read_parameters
    use initialize
    use schemes


    implicit none
    
    type(parameters)                           :: param
    type(drop_type), dimension(:), allocatable :: drops
    type(drop_type), dimension(:), allocatable :: drops_init

    character(len=100) :: filename

    filename = '../input/parameters.dat'

    ! integer :: n
    ! real(pr) :: tn

    
    call read_param(param, filename)
    call display_parameters(param)

    allocate(drops(param%nb_drops))
    allocate(drops_init(param%nb_drops))

    call init_drops(param, drops)
    drops_init = drops

    ! open(unit=11,file='../output/model_0/vitesse.dat', action='write')


    ! do n=1,Nmax

    !     !Adimentionalization of the variable
    !     ! drops = adimentionalization(drops, drops_init)

    !     write(11,*)  tn, drops(2)%position, drops(2)%velocity, drops(2)%radius, drops(2)%mass, drops(2)%temperature

    !     !write the initial drop field
    !     call out(NM, n, tn, drops)

    !     !Advance of on step of time
    !     drops = advance(NM, drops)


    !     tn = n*NM%time_step

    !     ! print*, 'tn:', tn

    ! end do

    ! close(11)

    deallocate(drops, drops_init)



    
end program embruns


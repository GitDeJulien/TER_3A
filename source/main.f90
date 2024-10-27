program embruns

    use precision
    use schemes


    implicit none
    
    type(drop_type) :: drops(nb_drops)
    type(drop_type) :: drops_init(nb_drops)
    type(NM_type)   :: NM

    integer :: n, Nmax
    real(pr) :: tn


    call init_drops(drops)
    drops_init = drops
    call init_model(NM)
    tn = 0.0
    Nmax = 100

    do n=1,Nmax

        !Adimentionalization of the variable
        ! drops = adimentionalization(drops, drops_init)

        !write the initial drop field
        call out(NM, n, tn, drops)

        !Advance of on step of time
        drops = advance(NM, drops)


        tn = n*NM%time_step

        ! print*, 'tn:', tn

    end do


    
end program embruns


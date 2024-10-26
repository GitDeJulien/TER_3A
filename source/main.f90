program embruns

    use precision
    use schemes


    implicit none
    
    type(drop_type) :: drops(nb_drops)
    type(NM_type)   :: NM

    real(pr) :: tn


    call init_drops(drops)
    call init_model(NM)
    tn = 0.0


    do while (tn <= NM%tmax) 

        !write the initial drop field
        call out(NM, tn, drops)

        !Advance of on step of time
        drops = advance(NM, drops)

        tn = tn + NM%time_step

    end do


    
end program embruns


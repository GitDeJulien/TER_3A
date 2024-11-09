module elementray_functions

    use precision
    use read_parameters

    implicit none
    
contains


    function Re_p(param, r_p, v_p) result(res)

        !In
        type(parameters), intent(in) :: param
        real(pr), intent(in)         :: r_p
        real(pr), dimension(param%dim), intent(in) :: v_p 

        !Out
        real(pr), dimension(param%dim) :: res

        !Local
        real(pr), dimension(param%dim) :: v_s

        v_s = v_p - param%U_air
        res = 2.0*r_p*abs(v_s)/param%nu_air

    end function Re_p

    function tau_p(param, r_p, m_p) result(res)

        !In
        type(parameters), intent(in) :: param
        real(pr), intent(in)         :: r_p
        real(pr), intent(in)         :: m_p 

        !Out
        real(pr)     :: res

        res = (2.0*r_p**2*(m_p/(4./3.*acos(-1.0)*r_p**3)))/(9.0*param%mu_air);

    end function tau_p

    function F(param, r_p, v_p, m_p) result(res)

        !In
        type(parameters), intent(in) :: param
        real(pr), intent(in)         :: r_p, m_p
        real(pr), dimension(param%dim), intent(in) :: v_p 

        !Out
        real(pr), dimension(param%dim) :: res

        !Local
        real(pr) :: tau_D
        real(pr), dimension(param%dim) :: fct

        tau_D = tau_p(param, r_p, m_p)
        fct = (m_p/tau_D)*(param%U_air - v_p)
        if (param%dim == 1) then
            res = fct
        else if (param%dim == 2) then
            res(1) = fct(1)
            res(2) = m_p*param%gravity + fct(2)
        else if (param%dim == 3) then
            res(1) = fct(1)
            res(2) = m_p*param%gravity + fct(2)
            res(3) = fct(3)
        else
            print*, "Dimension greater than 3 not allowed!"
        end if


    end function F

    function f_v(param, r_p, v_p) result(res)

        !In
        type(parameters), intent(in) :: param
        real(pr), intent(in)         :: r_p
        real(pr), dimension(param%dim), intent(in) :: v_p 

        !Out
        real(pr), dimension(param%dim) :: res

        res = 1.0 + sqrt(Re_p(param, r_p,v_p))/4.0

    end function f_v



    function Dv_star_function(param, r_p, v_p) result(Dv_star)

        !! Function computing the Dv* expression taking 
        !! the vector Xp_n=(xp_n, vp_n, rp_n, mp_n, Tp_n)
        !! Dv* is the modified diffusivity for water vapour

        implicit none

        !In
        type(parameters), intent(in) :: param
        real(pr), intent(in)         :: r_p
        real(pr), dimension(param%dim), intent(in) :: v_p 

        !Out
        real(pr), dimension(param%dim) :: Dv_star

        !Local
        real(pr) :: b, c
        real(pr), dimension(param%dim) :: a


        a = f_v(param, r_p,v_p)*param%D_v
        b = r_p/(r_p + param%Delta_v)
        c = (param%D_v/(r_p*param%alpha_c))*sqrt((2.0*pi*param%M_w)/(param%R_g*param%T_air))

        Dv_star = a/(b+c)

    end function Dv_star_function


    ! function ka_star_function(Xp_n) result(ka_star)

    !     !! Function computing the ka* expression taking 
    !     !! the vector Xp_n=(xp_n, vp_n, rp_n, mp_n, Tp_n)
    !     !! ka* is the modified thermal conductivity of air

    !     implicit none

    !     !In
    !     real(pr), dimension(:), intent(in) :: Xp_n

    !     !Out
    !     real(pr)            :: ka_star

    !     ka_star = fh*ka / ((Xp_n(3)/(Xp_n(3)+Delta_T)) + ka/(Xp_n(3)*alpha_T*rho_a*c_pa)*sqrt(2*pi*Ma/(Rg*T_air)))


    ! end function ka_star_function

    ! function M_function(Xp_n) result(Mp_n)

    !     !! Function computing the Mp_n expression taking 
    !     !! the vector X_n=(xp_n, vp_n, rp_n, mp_n, Tp_n)

    !     implicit none

    !     !In
    !     real(pr), dimension(:), intent(in) :: Xp_n

    !     !Out
    !     real(pr)            :: Mp_n

    !     Mp_n = (4*pi*Dv_star_function(Xp_n)*Mw*pv_sat)/(Rg*T_air) * &
    !     (Q_rh - T_air/Xp_n(5)*exp(Lv*Mw/Rg*(1/T_air - 1/Xp_n(5)) + &
    !     2*Mw*Gamma_p/(Rg*rhop_w*Xp_n(3)*Xp_n(5)) - &
    !     Ions*Phi_s*m_s*(Mw/Ms)/(Xp_n(4) - m_s)))


    ! end function M_function
    
end module elementray_functions
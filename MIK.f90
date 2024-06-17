    !  MIK.f90

    !****************************************************************************
    !
    !  PROGRAM: MIK model - Malama & Izmodenov & Korolkov model
    !
    !****************************************************************************
    
    ! ceiling(a) возвращает наименьшее целое число, большее или равное a. Тип - integer по умолчанию
	
	
	include "Storage_modul.f90"

    module My_func                     ! Модуль интерфейсов для внешних функций
        real(8), parameter :: MF_par_pi = acos(-1.0_8) 
        real(8), parameter :: MF_meDmp = (1.0_8/1836.15_8)   ! Отношение массы электрона к массе протона

        interface
        
            !@cuf attributes(host, device) & 
            real(8) pure function calc_square(n) ! Вычисляет площадь грани под номером n в общем списке граней
                use STORAGE, only: gl_all_Gran, gl_x, gl_y, gl_z
                implicit none
                integer, intent (in) :: n
            end function

            !@cuf attributes(host, device) & 
            real(8) pure function  tetrahedron(t1, t2, t3, t4)
                implicit none
                real(8), intent(in) :: t1(3), t2(3), t3(3), t4(3)
            end function

        end interface
        
        contains 

        !@cuf attributes(host, device) & 
        subroutine Calc_sourse_MF(plasma, fluid, sourse, zone, n_He_, MAS_PUI_)  ! Считаются мультифлюидные источники
            ! Variables
            use GEO_PARAM, only: par_a_2, par_n_p_LISM, par_Kn, par_pi_8, par_n_H_LISM_, par_n_sort
            implicit none
            real(8), intent(in) :: plasma(9)
            real(8), intent(in) :: fluid(5, par_n_sort)
            real(8), intent(in), OPTIONAL :: n_He_
            real(8), intent(in), optional :: MAS_PUI_(2)
            real(8) :: rho_Th, p_Th, p_Pui, rho_Pui
            real(8), intent(out) :: sourse(5, par_n_sort + 1)  ! (масса, три импульса и энергия)
            integer(4), intent(in) :: zone
            integer(4) :: al
            
            integer(4) :: i, kk(par_n_sort)
            real(8) :: U_M_H(par_n_sort), U_H(par_n_sort), sigma(par_n_sort), nu(par_n_sort), S1, S2, n_He
            real(8) :: MAS_PUI(2)
            
            MAS_PUI = 0.0
            n_He = 0.0
            if(PRESENT(n_He_)) n_He = n_He_
            !if(PRESENT(MAS_PUI_)) MAS_PUI = MAS_PUI_

            al = 1
            if(zone <= 2) al = 2
            rho_Pui = MAS_PUI(1)
            call Sootnosheniya(plasma(1), plasma(5), n_He, MAS_PUI(1), MAS_PUI(2), al, rho_Th = rho_Th, p_Th = p_Th, &
                               p_Pui = p_Pui )

            if(rho_Pui < 0.000001) then
                p_Pui = 0.0_8
                rho_Pui = 0.000001_8
            end if

            if(rho_Th < 0.000001) then
                p_Th = 0.0_8
                rho_Th = 0.000001_8
            end if

            if(p_Th < 0.0) p_Th = 0.000001

            !!   -----------------
            ! p_Pui = 0.0_8
            ! rho_Pui = 0.000001_8
            ! rho_Th = plasma(1)
            ! p_Th = plasma(5)/2.0
            !! 
            
            sourse = 0.0
            kk = 0
            kk(zone) = 1
            S1 = 0.0
            S2 = 0.0
            
            ! Body of Calc_sourse_MF
            !if(zone > 2) then
            if(.True.) then
                do i = 1, par_n_sort
                    U_M_H(i) = sqrt( (plasma(2) - fluid(2, i))**2 + (plasma(3) - fluid(3, i))**2 + (plasma(4) - fluid(4, i))**2 + &
                    (64.0 / (9.0 * par_pi_8)) * (plasma(5) / plasma(1) + 2.0 * fluid(5, i) / fluid(1, i)) )
                    U_H(i) = sqrt( (plasma(2) - fluid(2, i))**2 + (plasma(3) - fluid(3, i))**2 + (plasma(4) - fluid(4, i))**2 + &
                    (4.0 / par_pi_8) * (plasma(5) / plasma(1) + 2.0 * fluid(5, i) / fluid(1, i)) )
                    sigma(i) = (1.0 - par_a_2 * log(U_M_H(i)))**2
                    nu(i) = plasma(1) * fluid(1, i) * U_M_H(i) * sigma(i)
                end do

                do i = 1, par_n_sort
                    sourse(2, 1) =  sourse(2, 1) + nu(i) * (fluid(2, i) - plasma(2))
                    sourse(3, 1) =  sourse(3, 1) + nu(i) * (fluid(3, i) - plasma(3))
                    sourse(4, 1) =  sourse(4, 1) + nu(i) * (fluid(4, i) - plasma(4))
                    sourse(5, 1) = sourse(5, 1) + nu(i) * ( (fluid(2, i)**2 + fluid(3, i)**2 + fluid(4, i)**2 - &
                        plasma(2)**2 - plasma(3)**2 - plasma(4)**2)/2.0 + (U_H(i)/U_M_H(i)) * ( 2.0 * fluid(5, i)/fluid(1, i) - plasma(5)/plasma(1) ) )
                end do

                do i = 1, par_n_sort
                    S1 = S1 + nu(i)
                    S2 = S2 + nu(i) * ( (plasma(2)**2 + plasma(3)**2 + plasma(4)**2)/2.0 + (U_H(i)/U_M_H(i)) * (plasma(5)/plasma(1)) )
                end do

            else
                do i = 1, par_n_sort
                    U_M_H(i) = sqrt( (plasma(2) - fluid(2, i))**2 + (plasma(3) - fluid(3, i))**2 + (plasma(4) - fluid(4, i))**2 + &
                    (64.0 / (9.0 * par_pi_8)) * (2.0 * p_Th / rho_Th + 2.0 * fluid(5, i) / fluid(1, i)) )
                    U_H(i) = sqrt( (plasma(2) - fluid(2, i))**2 + (plasma(3) - fluid(3, i))**2 + (plasma(4) - fluid(4, i))**2 + &
                    (4.0 / par_pi_8) * (2.0 * p_Th / rho_Th + 2.0 * fluid(5, i) / fluid(1, i)) )
                    sigma(i) = (1.0 - par_a_2 * log(U_M_H(i)))**2
                    nu(i) = rho_Th * fluid(1, i) * U_M_H(i) * sigma(i)
                end do
                
                do i = 1, par_n_sort
                    sourse(2, 1) =  sourse(2, 1) + nu(i) * (fluid(2, i) - plasma(2))
                    sourse(3, 1) =  sourse(3, 1) + nu(i) * (fluid(3, i) - plasma(3))
                    sourse(4, 1) =  sourse(4, 1) + nu(i) * (fluid(4, i) - plasma(4))
                    sourse(5, 1) = sourse(5, 1) + nu(i) * ( (fluid(2, i)**2 + fluid(3, i)**2 + fluid(4, i)**2 - &
                        plasma(2)**2 - plasma(3)**2 - plasma(4)**2)/2.0 + (U_H(i)/U_M_H(i)) * ( 2.0 * fluid(5, i)/fluid(1, i) - 2.0 * p_Th/rho_Th ) )
                end do

                do i = 1, par_n_sort
                    S1 = S1 + nu(i)
                    S2 = S2 + nu(i) * ( (plasma(2)**2 + plasma(3)**2 + plasma(4)**2)/2.0 + (U_H(i)/U_M_H(i)) * (2.0 * p_Th/rho_Th) )
                end do

                do i = 1, par_n_sort
                    U_M_H(i) = sqrt( (plasma(2) - fluid(2, i))**2 + (plasma(3) - fluid(3, i))**2 + (plasma(4) - fluid(4, i))**2 + &
                    (64.0 / (9.0 * par_pi_8)) * (2.0 * p_Pui / rho_Pui + 2.0 * fluid(5, i) / fluid(1, i)) )
                    U_H(i) = sqrt( (plasma(2) - fluid(2, i))**2 + (plasma(3) - fluid(3, i))**2 + (plasma(4) - fluid(4, i))**2 + &
                    (4.0 / par_pi_8) * (2.0 * p_Pui / rho_Pui + 2.0 * fluid(5, i) / fluid(1, i)) )
                    sigma(i) = (1.0 - par_a_2 * log(U_M_H(i)))**2
                    nu(i) = rho_Pui * fluid(1, i) * U_M_H(i) * sigma(i)
                end do

                do i = 1, par_n_sort
                    sourse(2, 1) =  sourse(2, 1) + nu(i) * (fluid(2, i) - plasma(2))
                    sourse(3, 1) =  sourse(3, 1) + nu(i) * (fluid(3, i) - plasma(3))
                    sourse(4, 1) =  sourse(4, 1) + nu(i) * (fluid(4, i) - plasma(4))
                    sourse(5, 1) = sourse(5, 1) + nu(i) * ( (fluid(2, i)**2 + fluid(3, i)**2 + fluid(4, i)**2 - &
                        plasma(2)**2 - plasma(3)**2 - plasma(4)**2)/2.0 + (U_H(i)/U_M_H(i)) * ( 2.0 * fluid(5, i)/fluid(1, i) - 2.0 * p_Pui/rho_Pui ) )
                end do


                do i = 1, par_n_sort
                    S1 = S1 + nu(i)
                    S2 = S2 + nu(i) * ( (plasma(2)**2 + plasma(3)**2 + plasma(4)**2)/2.0 + (U_H(i)/U_M_H(i)) * (2.0 * p_Pui/rho_Pui) )
                end do
            end if

            ! if( sourse(2, 1) < 0.0 .or.  sourse(2, 1) > 100000.0) then
            !             print*,  "NUN u3  potok	"
            !             print*, "___"
            !             print*, fluid(1, 1), fluid(1, 2), fluid(1, 3), fluid(1, 4), fluid(1, 5), fluid(1, 6)
            !             print*, "___"
            !             STOP
            !     end if
            
            sourse(2:5, 1) =  sourse(2:5, 1) * (par_n_p_LISM/par_Kn)
            
            
            ! do i = 1, par_n_sort
            !     sourse(1, i + 1) = (par_n_H_LISM_/par_Kn) * (kk(i) * S1 - nu(i))
            !     sourse(2, i + 1) = (par_n_H_LISM_/par_Kn) * (kk(i) * S1 * plasma(2) - nu(i) * fluid(2, i))
            !     sourse(3, i + 1) = (par_n_H_LISM_/par_Kn) * (kk(i) * S1 * plasma(3) - nu(i) * fluid(3, i))
            !     sourse(4, i + 1) = (par_n_H_LISM_/par_Kn) * (kk(i) * S1 * plasma(4) - nu(i) * fluid(4, i))
            !     sourse(5, i + 1) = (par_n_H_LISM_/par_Kn) * (kk(i) * S2 - nu(i) * ( (fluid(2, i)**2 + fluid(3, i)**2 + fluid(4, i)**2)/2.0 + &
            !         (U_H(i)/U_M_H(i)) * 2.0 * (fluid(5, i) / fluid(1, i)) ) )
            ! end do
            
        
        end subroutine Calc_sourse_MF

        subroutine get_bazis(ex, ey, ez)
            ! По заданному вектору ex подбирает ey, ez ему перпендикулярные, образующие правую тройку
            !! ex - должен быть единичным
            real(8), intent(in) :: ex(3)
            real(8), intent(out) :: ey(3)
            real(8), intent(out) :: ez(3)
            real(8) :: norm

            ez(1) = 1.0
            ez(2) = 0.0
            ez(3) = 0.0

            if( 1.0 - dabs(DOT_PRODUCT(ex, ez)) < 0.001) then
                ez(1) = 0.0
                ez(2) = 1.0
                ez(3) = 0.0
            end if

            ey(1) = ex(2) * ez(3) - ex(3) * ez(2)
            ey(2) = ex(3) * ez(1) - ex(1) * ez(3)
            ey(3) = ex(1) * ez(2) - ex(2) * ez(1)
            norm = sqrt(ey(1)**2 + ey(2)**2 + ey(3)**2)
            ey = ey/norm

            ez(1) = ex(2) * ey(3) - ex(3) * ey(2)
            ez(2) = ex(3) * ey(1) - ex(1) * ey(3)
            ez(3) = ex(1) * ey(2) - ex(2) * ey(1)
            norm = sqrt(ez(1)**2 + ez(2)**2 + ez(3)**2)
            ez = ez/norm
        end subroutine get_bazis
        
        real(8) pure function kvv(x, y, z)
            real(8), intent (in) :: x, y, z
            kvv = x**2 + y**2 + z**2
        end function kvv
        
        !@cuf attributes(host, device) & 
        real(8) pure function polar_angle(x, y)
            use GEO_PARAM
            implicit none
            real(8), intent(in) :: x, y

            if (dabs(x) + dabs(y) < 0.00001 / par_R_character) then
                polar_angle = 0.0_8
                return
            end if


            if (x < 0) then
                polar_angle = atan(y / x) + 1.0 * par_pi_8
            elseif (x > 0 .and. y >= 0) then
                polar_angle = atan(y / x)
            elseif (x > 0 .and. y < 0) then
                polar_angle = atan(y / x) + 2.0 * par_pi_8
            elseif (y > 0 .and. x >= 0 .and. x <= 0) then
                polar_angle = par_pi_8 / 2.0
            elseif (y < 0 .and. x >= 0 .and. x <= 0) then
                polar_angle =  3.0 * par_pi_8 / 2.0
            end if

            return
        end function polar_angle
        
        !@cuf attributes(host, device) & 
        real(8) pure function angle_cilindr(x, par_al1)
            ! Функция для неравномерного распределения сетки по углу
            ! x от 0 до 1 и возвращает функция от 0 до 1
            implicit none
            real(8), intent(in) :: x, par_al1

            if (x < 0.5) then
                angle_cilindr = x * (par_al1 - 6 * (-1 + par_al1) * x + 8 * (-1 + par_al1) * x**2)
            else
                angle_cilindr = 3 + par_al1 * (-1 + x) * (-1 + 2 * x) * (-3 + 4 *x) - 2 * x * (6 + x *(-9 + 4 * x))
            end if
            
            return
        end function angle_cilindr
        
        !@cuf attributes(host, device) & 
        real(8) pure function sgushenie_1(x, par_al1)
            ! x от 0 до 1 и возвращает функция от 0 до 1
            ! Сгущение точек к обоим концам отрезка
            implicit none
            real(8), intent(in) :: x, par_al1

            sgushenie_1 = x * (par_al1 - 3.0 * (par_al1 - 1.0) * x + 2.0 * (par_al1 - 1.0)*x*x)
            
            return
        end function sgushenie_1
        
        !@cuf attributes(host, device) & 
        real(8) pure function sgushenie_2(x, all)
            ! x от 0 до 1 и возвращает функция от 0 до 1
            ! Сгущение точек к обоим концам отрезка (сильнее, чем предыдущая функция)
            implicit none
            real(8), intent(in) :: x, all

            sgushenie_2 = all * x - 10 * (-1 + all) * x**3 + 15 * (-1 + all) * x**4 - 6 * (-1 + all) * x**5
            
            return
        end function sgushenie_2
        
        !@cuf attributes(host, device) & 
        real(8) pure function sgushenie_3(x, all)
            ! x от 0 до 1 и возвращает функция от 0 до 1
            ! Сгущение точек к 1
            implicit none
            real(8), intent(in) :: x, all

            sgushenie_3 = -(-x + 1)**all + 1.0
            
            return
        end function sgushenie_3
        
        !@cuf attributes(host, device) & 
        real(8) pure function sgushenie_4(x, x0, y0)
            ! x от 0 до 1 и возвращает функция от 0 до 1
            ! Сгущение точек за BS
            implicit none
            real(8), intent(in) :: x, x0, y0

            if(x < x0) then
                sgushenie_4 = (x - x0) * y0/x0 + y0
            else
                sgushenie_4 = (1.0 - y0) * (x - x0)/(1.0 - x0) + y0
            end if
            
            return
        end function sgushenie_4
        
        !@cuf attributes(host, device) & 
        integer(4) pure function signum(x)
            implicit none
            real(8), intent(in) :: x
            
            if (x > 0) then
                signum = 1
                return
            else if (x < 0) then
                signum = -1
                return
            else 
                signum = 0
                return
            end if
        end function signum
        
        !@cuf attributes(host, device) & 
        real(8) pure function minmod(x, y)
            implicit none
            real(8), intent(in) :: x, y
            
            if (signum(x) + signum(y) == 0) then
                minmod = 0.0_8
                return
            else
                minmod = ((signum(x) + signum(y)) / 2.0) * min(dabs(x), dabs(y)) ! minmod
                return   
            end if
        end function minmod
        
        
        !@cuf attributes(host, device) & 
        real(8) pure function  linear(x1, t1, x2, t2, x3, t3, y)
            ! Главное значение с параметрами 2
            ! Строим линии между 1 и 2,  2 и 3, потом находим минмодом значение в y
            implicit none
            real(8), intent(in) :: x1, x2, x3, y, t1, t2, t3
            real(8) :: d

            d = minmod((t1 - t2) / (x1 - x2), (t2 - t3) / (x2 - x3))
            linear =  (d * (y - x2) + t2)
            return
            
        end function linear
        
        !@cuf attributes(host, device) & 
        real(8) pure function  linear2(x1, t1, x2, t2, y)
            ! Просто строим линию по параметрам 1 и 2
            implicit none
            real(8), intent(in) :: x1, x2, y, t1, t2
            real(8) :: d

            d = (t1 - t2) / (x1 - x2)
            linear2 =  (d * (y - x2) + t2)
            return
            
        end function linear2
        
        !@cuf attributes(host, device) & 
        subroutine polyar_skorost(phi, Vy, Vz, Vr, Vphi)
            ! Variables
            implicit none
            real(8), intent(in) :: phi, Vy, Vz
            real(8), intent(out) :: Vr, Vphi

            Vr = Vy * cos(phi) + Vz * sin(phi)
            Vphi = Vz * cos(phi) - Vy * sin(phi)

        end subroutine polyar_skorost
        
        !@cuf attributes(host, device) & 
        subroutine dekard_polyar_skorost(phi, Vr, Vphi, Vy, Vz)
            ! Variables
            implicit none
            real(8), intent(in) :: phi, Vr, Vphi 
            real(8), intent(out) :: Vy, Vz

            Vy = Vr * cos(phi) - Vphi * sin(phi)
            Vz = Vr * sin(phi) + Vphi * cos(phi)

        end subroutine dekard_polyar_skorost
        
        real(8) pure function MK_sigma(x)
            USE GEO_PARAM
            real(8), intent (in) :: x
            MK_sigma = (1.0 - par_a_2 * log(x))**2
        end function MK_sigma

        real(8) pure function MK_sigma2(x, y)
            USE GEO_PARAM
            real(8), intent (in) :: x, y
            MK_sigma2 = (1.0 - par_a_2 * log(x * y))**2
        end function MK_sigma2

        !@cuf attributes(host, device) &
        subroutine Sootnosheniya(rho, p, rho_He, rho_Pui, T_Pui, al, rho_Th, rho_E, p_Th, p_Pui, T_Th, T_E)
            ! Функция, определяющая температуры и концентрации гелия, пикапов и т.д.
            ! al - это заряд гелия
            ! если al = 1 то вне гелиопаузы
            ! если al = 2, то внутри гелиопаузы
            ! Th - термальные протоны, He - гелий, Pui - пикапы, E - электроны
            ! без параметров, это общие (те, что считаются в МГД)
            implicit none
            real(8), intent(in) :: rho, p, rho_He, rho_Pui, T_Pui
            integer(4), intent(in) :: al
            real(8), intent(out), optional :: rho_Th, rho_E, p_Th, p_Pui, T_Th, T_E

            if(PRESENT(rho_Th)) then
                rho_Th = -(MF_meDmp * al * rho_He + 4.0 * (-rho + rho_He))/(4.0 * (1.0 + MF_meDmp)) - rho_Pui
            end if

            if(PRESENT(rho_E)) then
                rho_E = MF_meDmp * (4.0 * rho + (-4.0 + al) * rho_He)/(4.0 * (1.0 + MF_meDmp))
            end if

            if(PRESENT(p_Th)) then
                p_Th = (p - rho_Pui * T_Pui) * (-4.0 * rho + (4.0 + al * MF_meDmp) * rho_He + 4.0 * (1.0 + MF_meDmp) * rho_Pui)/&
                (-8.0 * rho + (7.0 - al + (-1.0 + al) * MF_meDmp) * rho_He + 4.0 * (1.0 + MF_meDmp) * rho_Pui)
            end if

            if(PRESENT(p_Pui)) then
                p_Pui = T_pui * rho_Pui
            end if

            if(PRESENT(T_Th)) then
                T_Th = -4.0 * (1.0 + MF_meDmp) * (p - rho_Pui * T_Pui)/&
                (-8.0 * rho + (7.0 - al + (-1.0 + al)*MF_meDmp) * rho_He + 4.0 * (1.0 - MF_meDmp) * rho_Pui)
            end if

            if(PRESENT(T_E)) then
                if(PRESENT(T_Th)) then
                    T_E = T_Th
                else
                    T_E = -4.0 * (1.0 + MF_meDmp) * (p - rho_Pui * T_Pui)/&
                    (-8.0 * rho + (7.0 - al + (-1.0 + al)*MF_meDmp) * rho_He + 4.0 * (1.0 - MF_meDmp) * rho_Pui)
                end if
            end if


            return
        end subroutine Sootnosheniya

	end module My_func
	
	! Описание модулей
	
    include "Solvers.f90"
	include "cgod_3D.f90"
    include "Interpolation.f90"
    include "Help_func.f90"
    include "Move_func.f90"
	include "TVD.f90"
	include "Surface_setting.f90"
	include "PUI.f90"
	include "Monte-Karlo.f90"
	!@cuf include "cuf_kernel.cuf"

    ! Вспомогательные функции

	!@cuf attributes(host, device) & 
    subroutine spherical_skorost(z, x, y, Vz, Vx, Vy, Vr, Vphi, Vtheta)
        ! Variables
        use My_func
        implicit none
        real(8), intent(in) :: x, y, z, Vx, Vy, Vz
        real(8), intent(out) :: Vr, Vphi, Vtheta
        real(8) :: r_1, the_1, phi_1

        r_1 = sqrt(x * x + y * y + z * z)
        the_1 = acos(z / r_1)
        phi_1 = polar_angle(x, y)

        Vr = Vx * sin(the_1) * cos(phi_1) + Vy * sin(the_1) * sin(phi_1) + Vz * cos(the_1);
        Vtheta = Vx * cos(the_1) * cos(phi_1) + Vy * cos(the_1) * sin(phi_1) - Vz * sin(the_1);
        Vphi = -Vx * sin(phi_1) + Vy * cos(phi_1);
	end subroutine spherical_skorost

	!@cuf attributes(host, device) & 
    subroutine dekard_skorost(z, x, y, Vr, Vphi, Vtheta, Vz, Vx, Vy)
        use My_func
        implicit none
        real(8), intent(in) :: x, y, z,  Vr, Vphi, Vtheta
        real(8), intent(out) :: Vx, Vy, Vz
        real(8) :: r_2, the_2, phi_2

        r_2 = sqrt(x * x + y * y + z * z);
        the_2 = acos(z / r_2);
        phi_2 = polar_angle(x, y);
        
        !print*, r_2, the_2, phi_2, Vr, Vphi, Vtheta
        !print*, sin(the_2), cos(phi_2), cos(the_2), sin(phi_2)

        Vx = Vr * sin(the_2) * cos(phi_2) + Vtheta * cos(the_2) * cos(phi_2) - Vphi * sin(phi_2);
        Vy = Vr * sin(the_2) * sin(phi_2) + Vtheta * cos(the_2) * sin(phi_2) + Vphi * cos(phi_2);
        Vz = Vr * cos(the_2) - Vtheta * sin(the_2);
	end subroutine dekard_skorost
	
	!@cuf attributes(host, device) & 
    real(8) pure function  tetrahedron(t1, t2, t3, t4)
        implicit none
        
        real(8), intent(in) :: t1(3), t2(3), t3(3), t4(3)
        real(8) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4
        real(8) :: a1, b1, c1, a2, b2, c2, a3, b3, c3
        real(8) :: det
        
        x1 = t1(1)
        y1 = t1(2)
        z1 = t1(3)
        
        x2 = t2(1)
        y2 = t2(2)
        z2 = t2(3)
        
        x3 = t3(1)
        y3 = t3(2)
        z3 = t3(3)
        
        x4 = t4(1)
        y4 = t4(2)
        z4 = t4(3)

        a1 = x2 - x1
        b1 = y2 - y1
        c1 = z2 - z1
        a2 = x3 - x1
        b2 = y3 - y1
        c2 = z3 - z1
        a3 = x4 - x1
        b3 = y4 - y1
        c3 = z4 - z1

        ! вычисление определителя матрицы
        det = a1 * (b2 * c3 - b3 * c2) - b1 * (a2 * c3 - a3 * c2) + c1 * (a2 * b3 - a3 * b2)

        ! вычисление объема тетраэдра
        tetrahedron = abs(det) / 6.0
    end function tetrahedron

    ! ****************************************************************************************************************************************************
    ! Блок функций для начального построения сетки

    subroutine Build_Mesh_start()        ! Начальное построение сетки
        use STORAGE
        use GEO_PARAM
        use My_func
        implicit none

        integer(4) :: i, j, k, N1, N2, N3, i1, kk, node, kk2, ni, num1
        real(8) :: r, phi, the, xx, x, y, z, rr, x2, y2, z2

        ! Сетка строится в декартовой системе координат, хотя сама сетка цилиндрически симметричная
        ! Нужно уметь преобразовывать x,y,z <---> x, r, phi  (x - ось симметрии), а также сферические координаты

        node = 1
        ! Задаём первый узел
        gl_x(node) = 0.0
        gl_y(node) = 0.0
        gl_z(node) = 0.0
        node = node + 1
        !gl_x = [gl_x, 0.0]
        !gl_y = [gl_y, 0.0]
        !gl_z = [gl_z, 0.0]

        N3 = size(gl_RAY_A(1, 1, :))
        N2 = size(gl_RAY_A(1, :, 1))
        N1 = size(gl_RAY_A(:, 1, 1))

        ! **************************************************************** Начнём с узлов сетки

        ! Цикл генерации точек на лучах А и их связывание с этими лучами ************************************************************
        do k = 1, N3
            do j = 1, N2
                do i = 1, N1
                    if (i == 1) then
                        gl_RAY_A(i, j, k) = 1
                        CYCLE
                    end if

                    if (k /= 1 .and. j == 1) then
                        gl_RAY_A(i, j, k) = i
                        CYCLE
                    end if

                    ! Вычисляем координаты текущего луча в пространстве
                    the = (j - 1) * par_pi_8/2.0/(N2 - 1)
                    phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1) ! (k - 1) * 2.0_8 * par_pi_8/(N3)
                    ! Вычисляем координаты точки на луче

                    ! до TS
                    if (i <= par_n_TS) then  ! До расстояния = par_R_character
                        if(i == 2) then
                            r =  par_R0 - (par_R_character - par_R0) * (DBLE(3 - 2)/(par_n_TS - 2))**par_kk1
                            if(r < 0.0) then
                                print*, "Error iouihjgfdcydygy  ", r
                                STOP
                            end if
                        else
                            r =  par_R0 + (par_R_character - par_R0) * (DBLE(i - 2)/(par_n_TS - 2))**par_kk1
                        end if
                    else if (i <= par_n_HP) then  ! До расстояния = par_R_character * 1.3
                        r = par_R_character + (i - par_n_TS) * 0.3 * par_R_character/(par_n_HP - par_n_TS)
                    else if (i <= par_n_BS) then  ! До расстояния = par_R_character * 2
                        r = 1.3 * par_R_character + (i - par_n_HP) * 0.7 * par_R_character/(par_n_BS - par_n_HP)
                    else  ! До расстояния = par_R_character * 1.3_8
                        !r = 2.0 * par_R_character + (i - par_n_BS) * (par_R_END - 2.0 * par_R_character)/(par_n_END - par_n_BS)
                        r = 2.0 * par_R_character + (par_R_END - 2.0 * par_R_character) * (DBLE(i- par_n_BS)/(par_n_END - par_n_BS))**par_kk2
                    end if

                    ! Создаём точку
                    !gl_x = [gl_x, r * cos(the)]
                    !gl_y = [gl_y, r * sin(the) * cos(phi)]
                    !gl_z = [gl_z, r * sin(the) * sin(phi)]
                    gl_x(node) = r * cos(the)
                    gl_y(node) = r * sin(the) * cos(phi)
                    gl_z(node) = r * sin(the) * sin(phi)
                    node = node + 1

                    gl_RAY_A(i, j, k) = node - 1
                end do
            end do
        end do

        N3 = size(gl_RAY_B(1, 1, :))
        N2 = size(gl_RAY_B(1, :, 1))
        N1 = size(gl_RAY_B(:, 1, 1))

        ! Цикл генерации точек на лучах B и их связывание с этими лучами ************************************************************
        do k = 1, N3
            do j = 1, N2
                do i = 1, N1

                    if (i == 1) then
                        gl_RAY_B(i, j, k) = 1
                        CYCLE
                    end if

                    ! Вычисляем координаты текущего луча в пространстве
                    the = par_pi_8/2.0 + (j) * par_triple_point/(N2)
                    phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1) !(k - 1) * 2.0_8 * par_pi_8/(N3)
                    ! Вычисляем координаты точки на луче

                    ! до TS
                    if (i <= par_n_TS) then  ! До расстояния = par_R_character
                        if(i == 2) then
                            r =  par_R0 - (par_R_character - par_R0) * (DBLE(3 - 2)/(par_n_TS - 2))**par_kk1
                            if(r < 0.0) then
                                print*, "Error iouihjgfdcydygy  ", r
                                STOP
                            end if
                        else
                            r =  par_R0 + (par_R_character - par_R0) * (DBLE(i - 2)/(par_n_TS - 2))**par_kk1
                        end if
                        !r =  par_R0 + (par_R_character - par_R0) * (REAL(i - 2, KIND = 4)/(par_n_TS - 2))**par_kk1
                        !print *, r
                        !pause
                    else if (i <= par_n_HP) then  ! До расстояния = par_R_character * 1.3
                        xx = 1.3 * par_R_character * (1.0_8 - cos(the - par_pi_8/2.0_8))/cos(the - par_pi_8/2.0_8)
                        r = par_R_character + (i - par_n_TS) * (0.3 * par_R_character + xx) /(par_n_HP - par_n_TS)
                    end if

                    ! Создаём точку
                    !gl_x = [gl_x, r * cos(the)]
                    !gl_y = [gl_y, r * sin(the) * cos(phi)]
                    !gl_z = [gl_z, r * sin(the) * sin(phi)]

                    gl_x(node) = r * cos(the)
                    gl_y(node) = r * sin(the) * cos(phi)
                    gl_z(node) = r * sin(the) * sin(phi)
                    node = node + 1

                    gl_RAY_B(i, j, k) = node - 1
                end do
            end do
        end do


        N3 = size(gl_RAY_C(1, 1, :))
        N2 = size(gl_RAY_C(1, :, 1))
        N1 = size(gl_RAY_C(:, 1, 1))

        ! Цикл генерации точек на лучах C и их связывание с этими лучами ************************************************************
        do k = 1, N3
            do j = 1, N2
                do i = 1, N1

                    if(i == 1) then
                        gl_RAY_C(i, j, k) = gl_RAY_B(par_n_HP, j, k)
                        CYCLE
                    end if

                    ! Вычисляем координаты текущего луча в пространстве
                    phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1) !(k - 1) * 2.0_8 * par_pi_8/(N3)
                    ! Вычисляем координаты точки на луче

                    x = gl_x(gl_RAY_B(par_n_HP, j, k))
                    y = gl_y(gl_RAY_B(par_n_HP, j, k))
                    z = gl_z(gl_RAY_B(par_n_HP, j, k))
                    rr = (y**2 + z**2)**(0.5)

                    if (i <= par_n_BS - par_n_HP + 1) then
                        r = rr + (i - 1) * (2.0 * par_R_character - rr)/(par_n_BS - par_n_HP)
                    else
                        !r = 2.0 * par_R_character + (i - (par_n_BS - par_n_HP + 1)) * (par_R_END - 2.0 * par_R_character) /(N1 - (par_n_BS - par_n_HP + 1) )
                        r = 2.0 * par_R_character + (DBLE(i - (par_n_BS - par_n_HP + 1))/(N1 - (par_n_BS - par_n_HP + 1) ))**par_kk2 * (par_R_END - 2.0 * par_R_character)
                    end if

                    ! Создаём точку
                    !gl_x = [gl_x, x]
                    !gl_y = [gl_y, r * cos(phi)]
                    !gl_z = [gl_z, r * sin(phi)]

                    gl_x(node) = x
                    gl_y(node) = r * cos(phi)
                    gl_z(node) = r * sin(phi)
                    node = node + 1

                    gl_RAY_C(i, j, k) = node - 1
                end do
            end do
        end do

        N3 = size(gl_RAY_O(1, 1, :))
        N2 = size(gl_RAY_O(1, :, 1))
        N1 = size(gl_RAY_O(:, 1, 1))

        xx = gl_x(gl_RAY_B(par_n_HP, par_m_BC, 1))

        ! Цикл генерации точек на лучах O и их связывание с этими лучами ************************************************************
        do k = 1, N3
            do j = 1, N2
                !x = xx - j * (xx - par_R_LEFT)/N2
                x = xx - (DBLE(j)/N2)**par_kk3 * (xx - par_R_LEFT)
                do i = 1, N1

                    ! Вычисляем координаты текущего луча в пространстве
                    phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1) !(k - 1) * 2.0_8 * par_pi_8/(N3)
                    ! Вычисляем координаты точки на луче


                    if (i <= par_n_BS - par_n_HP + 1) then
                        r = 1.3 * par_R_character + (i - 1) * par_R_character * (0.7)/(par_n_BS - par_n_HP)
                    else
                        !r = 2.0 * par_R_character + (i - (par_n_BS - par_n_HP + 1)) * (par_R_END - 2.0 * par_R_character) /(N1 - (par_n_BS - par_n_HP + 1) )
                        r = 2.0 * par_R_character + (DBLE(i - (par_n_BS - par_n_HP + 1))/(N1 - (par_n_BS - par_n_HP + 1) ))**par_kk2 * (par_R_END - 2.0 * par_R_character)
                    end if

                    ! Создаём точку
                    !gl_x = [gl_x, x]
                    !gl_y = [gl_y, r * cos(phi)]
                    !gl_z = [gl_z, r * sin(phi)]

                    gl_x(node) = x
                    gl_y(node) = r * cos(phi)
                    gl_z(node) = r * sin(phi)
                    node = node + 1

                    gl_RAY_O(i, j, k) = node - 1
                end do
            end do
        end do

        N3 = size(gl_RAY_K(1, 1, :))
        N2 = size(gl_RAY_K(1, :, 1))
        N1 = size(gl_RAY_K(:, 1, 1))

        ! Цикл генерации точек на лучах K и их связывание с этими лучами ************************************************************
        do k = 1, N3
            do j = 1, N2
                do i = 1, N1

                    if (i == 1) then
                        gl_RAY_K(i, j, k) = 1
                        CYCLE
                    end if

                    if (k /= 1 .and. j == 1) then
                        gl_RAY_K(i, j, k) = gl_RAY_K(i, 1, 1)
                        CYCLE
                    end if

                    ! Вычисляем координаты текущего луча в пространстве
                    the = par_pi_8/2.0 + par_triple_point + (N2 - j + 1) * (par_pi_8/2.0 - par_triple_point)/(N2)
                    phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1) !(k - 1) * 2.0_8 * par_pi_8/(N3)
                    ! Вычисляем координаты точки на луче
                    if(i == 2) then
                            r =  par_R0 - (par_R_character - par_R0) * (DBLE(3 - 2)/(par_n_TS - 2))**par_kk1
                            if(r < 0.0) then
                                print*, "Error iouihjgfdcydygy  ", r
                                STOP
                            end if
                    else
                        r =  par_R0 + (par_R_character - par_R0) * (DBLE(i - 2)/(par_n_TS - 2))**par_kk1
                    end if
                        
                    !r =  par_R0 + (par_R_character - par_R0) * (REAL(i - 2, KIND = 4)/(par_n_TS - 2))**par_kk1


                    ! Создаём точку
                    !gl_x = [gl_x, r * cos(the)]
                    !gl_y = [gl_y, r * sin(the) * cos(phi)]
                    !gl_z = [gl_z, r * sin(the) * sin(phi)]

                    gl_x(node) = r * cos(the)
                    gl_y(node) = r * sin(the) * cos(phi)
                    gl_z(node) = r * sin(the) * sin(phi)
                    node = node + 1

                    gl_RAY_K(i, j, k) = node - 1
                end do
            end do
        end do

        N3 = size(gl_RAY_D(1, 1, :))
        N2 = size(gl_RAY_D(1, :, 1))
        N1 = size(gl_RAY_D(:, 1, 1))

        ! Цикл генерации точек на лучах D и их связывание с этими лучами ************************************************************
        do k = 1, N3
            do j = 1, N2
                do i = 1, N1

                    if (i == 1) then
                        if (j < N2) then
                            gl_RAY_D(i, j, k) = gl_RAY_K(par_n_TS, j, k)
                        else
                            gl_RAY_D(i, j, k) = gl_RAY_B(par_n_TS, par_m_BC, k)
                        end if
                        CYCLE
                    end if

                    if (k /= 1 .and. j == 1) then
                        gl_RAY_D(i, j, k) = gl_RAY_D(i, 1, 1)
                        CYCLE
                    end if

                    ! Вычисляем координаты текущего луча в пространстве
                    phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1) !(k - 1) * 2.0_8 * par_pi_8/(N3)
                    ! Вычисляем координаты точки на луче

                    if (j < N2) then
                        xx = gl_x(gl_RAY_K(par_n_TS, j, k))
                        y = gl_y(gl_RAY_K(par_n_TS, j, k))
                        z = gl_z(gl_RAY_K(par_n_TS, j, k))
                    else
                        xx = gl_x(gl_RAY_B(par_n_TS, par_m_BC, k))
                        y = gl_y(gl_RAY_B(par_n_TS, par_m_BC, k))
                        z = gl_z(gl_RAY_B(par_n_TS, par_m_BC, k))
                    end if

                    r = sqrt(y**2 + z**2)
                    !x = xx + (i - 1) * (par_R_LEFT - xx)/(N1 - 1)
                    x = xx + (DBLE(i - 1)/(N1 - 1))**par_kk3 * (par_R_LEFT - xx)


                    ! Создаём точку
                    !gl_x = [gl_x, x]
                    !gl_y = [gl_y, r * cos(phi)]
                    !gl_z = [gl_z, r * sin(phi)]

                    gl_x(node) = x
                    gl_y(node) = r * cos(phi)
                    gl_z(node) = r * sin(phi)
                    node = node + 1

                    gl_RAY_D(i, j, k) = node - 1
                end do
            end do
        end do

        N3 = size(gl_RAY_E(1, 1, :))
        N2 = size(gl_RAY_E(1, :, 1))
        N1 = size(gl_RAY_E(:, 1, 1))

        ! Цикл генерации точек на лучах E и их связывание с этими лучами ************************************************************
        do k = 1, N3
            do j = 1, N2
                do i = 1, N1

                    if (i == 1) then
                        gl_RAY_E(i, j, k) = gl_RAY_D(j + 1, size(gl_RAY_D(i + 1, :, k)), k)
                        CYCLE
                    end if

                    if (i == N1) then
                        gl_RAY_E(i, j, k) = gl_RAY_O(1, j, k)
                        CYCLE
                    end if

                    x = gl_x(gl_RAY_E(1, j, k))
                    y = gl_y(gl_RAY_E(1, j, k))
                    z = gl_z(gl_RAY_E(1, j, k))

                    x2 = gl_x(gl_RAY_O(1, j, k))
                    y2 = gl_y(gl_RAY_O(1, j, k))
                    z2 = gl_z(gl_RAY_O(1, j, k))


                    ! Создаём точку
                    !gl_x = [gl_x, x + (x2 - x) * (i - 1)/(N1 - 1)]
                    !gl_y = [gl_y, y + (y2 - y) * (i - 1)/(N1 - 1)]
                    !gl_z = [gl_z, z + (z2 - z) * (i - 1)/(N1 - 1)]

                    gl_x(node) = x + (x2 - x) * (i - 1)/(N1 - 1)
                    gl_y(node) = y + (y2 - y) * (i - 1)/(N1 - 1)
                    gl_z(node) = z + (z2 - z) * (i - 1)/(N1 - 1)
                    node = node + 1

                    gl_RAY_E(i, j, k) = node - 1
                end do
            end do
        end do

        ! Строим сами ячейки и связываем их с точками

        ! А - группа ячеек ************************************************************************************************************************
        N3 = size(gl_Cell_A(1, 1, :))
        N2 = size(gl_Cell_A(1, :, 1))
        N1 = size(gl_Cell_A(:, 1, 1))

        i1 = 1

        do k = 1, N3
            do j = 1, N2
                do i = 1, N1

                    kk = k + 1
                    if (kk > par_l_phi) kk = 1

                    if (j < par_m_A) then
                        gl_all_Cell(1, i1) = gl_RAY_A(i, j, k)
                        gl_all_Cell(2, i1) = gl_RAY_A(i + 1, j, k)
                        gl_all_Cell(3, i1) = gl_RAY_A(i + 1, j + 1, k)
                        gl_all_Cell(4, i1) = gl_RAY_A(i, j + 1, k)
                        gl_all_Cell(5, i1) = gl_RAY_A(i, j, kk)
                        gl_all_Cell(6, i1) = gl_RAY_A(i + 1, j, kk)
                        gl_all_Cell(7, i1) = gl_RAY_A(i + 1, j + 1, kk)
                        gl_all_Cell(8, i1) = gl_RAY_A(i, j + 1, kk)
                    else if (j > par_m_A) then
                        if (i < par_n_HP) then
                            gl_all_Cell(1, i1) = gl_RAY_B(i, j - par_m_A, k)
                            gl_all_Cell(2, i1) = gl_RAY_B(i + 1, j - par_m_A, k)
                            gl_all_Cell(3, i1) = gl_RAY_B(i + 1, j - par_m_A + 1, k)
                            gl_all_Cell(4, i1) = gl_RAY_B(i, j - par_m_A + 1, k)
                            gl_all_Cell(5, i1) = gl_RAY_B(i, j - par_m_A, kk)
                            gl_all_Cell(6, i1) = gl_RAY_B(i + 1, j - par_m_A, kk)
                            gl_all_Cell(7, i1) = gl_RAY_B(i + 1, j - par_m_A + 1, kk)
                            gl_all_Cell(8, i1) = gl_RAY_B(i, j - par_m_A + 1, kk)
                        else if (i == par_n_HP) then
                            gl_all_Cell(1, i1) = gl_RAY_B(i, j - par_m_A, k)
                            gl_all_Cell(2, i1) = gl_RAY_C(2, j - par_m_A, k)
                            gl_all_Cell(3, i1) = gl_RAY_C(2, j - par_m_A + 1, k)
                            gl_all_Cell(4, i1) = gl_RAY_B(i, j - par_m_A + 1, k)
                            gl_all_Cell(5, i1) = gl_RAY_B(i, j - par_m_A, kk)
                            gl_all_Cell(6, i1) = gl_RAY_C(2, j - par_m_A, kk)
                            gl_all_Cell(7, i1) = gl_RAY_C(2, j - par_m_A + 1, kk)
                            gl_all_Cell(8, i1) = gl_RAY_B(i, j - par_m_A + 1, kk)
                        else ! i > par_n_HP
                            gl_all_Cell(1, i1) = gl_RAY_C(i - par_n_HP + 1, j - par_m_A, k)
                            gl_all_Cell(2, i1) = gl_RAY_C(i - par_n_HP + 2, j - par_m_A, k)
                            gl_all_Cell(3, i1) = gl_RAY_C(i - par_n_HP + 2, j - par_m_A + 1, k)
                            gl_all_Cell(4, i1) = gl_RAY_C(i - par_n_HP + 1, j - par_m_A + 1, k)
                            gl_all_Cell(5, i1) = gl_RAY_C(i - par_n_HP + 1, j - par_m_A, kk)
                            gl_all_Cell(6, i1) = gl_RAY_C(i - par_n_HP + 2, j - par_m_A, kk)
                            gl_all_Cell(7, i1) = gl_RAY_C(i - par_n_HP + 2, j - par_m_A + 1, kk)
                            gl_all_Cell(8, i1) = gl_RAY_C(i - par_n_HP + 1, j - par_m_A + 1, kk)
                        end if
                    else ! j == par_m_A
                        if (i < par_n_HP) then
                            gl_all_Cell(1, i1) = gl_RAY_A(i, j, k)
                            gl_all_Cell(2, i1) = gl_RAY_A(i + 1, j, k)
                            gl_all_Cell(3, i1) = gl_RAY_B(i + 1, 1, k)
                            gl_all_Cell(4, i1) = gl_RAY_B(i, 1, k)
                            gl_all_Cell(5, i1) = gl_RAY_A(i, j, kk)
                            gl_all_Cell(6, i1) = gl_RAY_A(i + 1, j, kk)
                            gl_all_Cell(7, i1) = gl_RAY_B(i + 1, 1, kk)
                            gl_all_Cell(8, i1) = gl_RAY_B(i, 1, kk)
                        else if (i == par_n_HP) then
                            gl_all_Cell(1, i1) = gl_RAY_A(i, j, k)
                            gl_all_Cell(2, i1) = gl_RAY_A(i + 1, j, k)
                            gl_all_Cell(3, i1) = gl_RAY_C(2, 1, k)
                            gl_all_Cell(4, i1) = gl_RAY_B(i, 1, k)
                            gl_all_Cell(5, i1) = gl_RAY_A(i, j, kk)
                            gl_all_Cell(6, i1) = gl_RAY_A(i + 1, j, kk)
                            gl_all_Cell(7, i1) = gl_RAY_C(2, 1, kk)
                            gl_all_Cell(8, i1) = gl_RAY_B(i, 1, kk)
                        else ! i > par_n_HP
                            gl_all_Cell(1, i1) = gl_RAY_A(i, j, k)
                            gl_all_Cell(2, i1) = gl_RAY_A(i + 1, j, k)
                            gl_all_Cell(3, i1) = gl_RAY_C(i - par_n_HP + 2, 1, k)
                            gl_all_Cell(4, i1) = gl_RAY_C(i - par_n_HP + 1, 1, k)
                            gl_all_Cell(5, i1) = gl_RAY_A(i, j, kk)
                            gl_all_Cell(6, i1) = gl_RAY_A(i + 1, j, kk)
                            gl_all_Cell(7, i1) = gl_RAY_C(i - par_n_HP + 2, 1, kk)
                            gl_all_Cell(8, i1) = gl_RAY_C(i - par_n_HP + 1, 1, kk)
                        end if
                    end if

                    gl_Cell_A(i, j, k) = i1
                    i1 = i1 + 1

                end do
            end do
        end do

        ! B - группа ячеек ************************************************************************************************************************
        N3 = size(gl_Cell_B(1, 1, :))
        N2 = size(gl_Cell_B(1, :, 1))
        N1 = size(gl_Cell_B(:, 1, 1))

        do k = 1, N3
            do j = 1, N2
                do i = 1, N1

                    kk = k + 1
                    if (kk > par_l_phi) kk = 1

                    if (j < par_m_K) then
                        if (i < par_n_TS) then
                            gl_all_Cell(1, i1) = gl_RAY_K(i, j, k)
                            gl_all_Cell(2, i1) = gl_RAY_K(i, j + 1, k)
                            gl_all_Cell(3, i1) = gl_RAY_K(i + 1, j + 1, k)
                            gl_all_Cell(4, i1) = gl_RAY_K(i + 1, j, k)
                            gl_all_Cell(5, i1) = gl_RAY_K(i, j, kk)
                            gl_all_Cell(6, i1) = gl_RAY_K(i, j + 1, kk)
                            gl_all_Cell(7, i1) = gl_RAY_K(i + 1, j + 1, kk)
                            gl_all_Cell(8, i1) = gl_RAY_K(i + 1, j, kk)
                        else if (i == par_n_TS) then
                            gl_all_Cell(1, i1) = gl_RAY_K(i, j, k)
                            gl_all_Cell(2, i1) = gl_RAY_K(i, j + 1, k)
                            gl_all_Cell(3, i1) = gl_RAY_D(2, j + 1, k)
                            gl_all_Cell(4, i1) = gl_RAY_D(2, j, k)
                            gl_all_Cell(5, i1) = gl_RAY_K(i, j, kk)
                            gl_all_Cell(6, i1) = gl_RAY_K(i, j + 1, kk)
                            gl_all_Cell(7, i1) = gl_RAY_D(2, j + 1, kk)
                            gl_all_Cell(8, i1) = gl_RAY_D(2, j, kk)
                        else  ! i > par_n_TS
                            gl_all_Cell(1, i1) = gl_RAY_D(i - par_n_TS + 1, j, k)
                            gl_all_Cell(2, i1) = gl_RAY_D(i - par_n_TS + 1, j + 1, k)
                            gl_all_Cell(3, i1) = gl_RAY_D(i - par_n_TS + 1 + 1, j + 1, k)
                            gl_all_Cell(4, i1) = gl_RAY_D(i - par_n_TS + 1 + 1, j, k)
                            gl_all_Cell(5, i1) = gl_RAY_D(i - par_n_TS + 1, j, kk)
                            gl_all_Cell(6, i1) = gl_RAY_D(i - par_n_TS + 1, j + 1, kk)
                            gl_all_Cell(7, i1) = gl_RAY_D(i - par_n_TS + 1 + 1, j + 1, kk)
                            gl_all_Cell(8, i1) = gl_RAY_D(i - par_n_TS + 1 + 1, j, kk)
                        end if
                    else  ! j == par_m_K
                        if (i < par_n_TS) then
                            gl_all_Cell(1, i1) = gl_RAY_K(i, j, k)
                            gl_all_Cell(2, i1) = gl_RAY_B(i, size(gl_RAY_B(i, :, k)), k)
                            gl_all_Cell(3, i1) = gl_RAY_B(i + 1, size(gl_RAY_B(i, :, k)), k)
                            gl_all_Cell(4, i1) = gl_RAY_K(i + 1, j, k)
                            gl_all_Cell(5, i1) = gl_RAY_K(i, j, kk)
                            gl_all_Cell(6, i1) = gl_RAY_B(i, size(gl_RAY_B(i, :, k)), kk)
                            gl_all_Cell(7, i1) = gl_RAY_B(i + 1, size(gl_RAY_B(i, :, k)), kk)
                            gl_all_Cell(8, i1) = gl_RAY_K(i + 1, j, kk)
                        else if (i == par_n_TS) then
                            gl_all_Cell(1, i1) = gl_RAY_K(i, j, k)
                            gl_all_Cell(2, i1) = gl_RAY_B(i, size(gl_RAY_B(i, :, k)), k)
                            gl_all_Cell(3, i1) = gl_RAY_D(2, size(gl_RAY_D(2, :, k)), k)
                            gl_all_Cell(4, i1) = gl_RAY_D(2, size(gl_RAY_D(2, :, k)) - 1, k)
                            gl_all_Cell(5, i1) = gl_RAY_K(i, j, kk)
                            gl_all_Cell(6, i1) = gl_RAY_B(i, size(gl_RAY_B(i, :, k)), kk)
                            gl_all_Cell(7, i1) = gl_RAY_D(2, size(gl_RAY_D(2, :, k)), kk)
                            gl_all_Cell(8, i1) = gl_RAY_D(2, size(gl_RAY_D(2, :, k)) - 1, kk)
                        else   ! i > par_n_TS
                            gl_all_Cell(1, i1) = gl_RAY_D(i - par_n_TS + 1, j, k)
                            gl_all_Cell(2, i1) = gl_RAY_D(i - par_n_TS + 1, j + 1, k)
                            gl_all_Cell(3, i1) = gl_RAY_D(i - par_n_TS + 1 + 1, j + 1, k)
                            gl_all_Cell(4, i1) = gl_RAY_D(i - par_n_TS + 1 + 1, j, k)
                            gl_all_Cell(5, i1) = gl_RAY_D(i - par_n_TS + 1, j, kk)
                            gl_all_Cell(6, i1) = gl_RAY_D(i - par_n_TS + 1, j + 1, kk)
                            gl_all_Cell(7, i1) = gl_RAY_D(i - par_n_TS + 1 + 1, j + 1, kk)
                            gl_all_Cell(8, i1) = gl_RAY_D(i - par_n_TS + 1 + 1, j, kk)
                        end if
                    end if


                    gl_Cell_B(i, j, k) = i1
                    i1 = i1 + 1

                end do
            end do
        end do

        ! C - группа ячеек ************************************************************************************************************************
        N3 = size(gl_Cell_C(1, 1, :))
        N2 = size(gl_Cell_C(1, :, 1))
        N1 = size(gl_Cell_C(:, 1, 1))

        do k = 1, N3
            do j = 1, N2
                do i = 1, N1

                    kk = k + 1
                    if (kk > par_l_phi) kk = 1

                    if (j > 1) then
                        if (i < par_n_HP - par_n_TS + 1) then
                            gl_all_Cell(1, i1) = gl_RAY_E(i, j - 1, k)
                            gl_all_Cell(2, i1) = gl_RAY_E(i + 1, j - 1, k)
                            gl_all_Cell(3, i1) = gl_RAY_E(i + 1, j - 1 + 1, k)
                            gl_all_Cell(4, i1) = gl_RAY_E(i, j - 1 + 1, k)
                            gl_all_Cell(5, i1) = gl_RAY_E(i, j - 1, kk)
                            gl_all_Cell(6, i1) = gl_RAY_E(i + 1, j - 1, kk)
                            gl_all_Cell(7, i1) = gl_RAY_E(i + 1, j - 1 + 1, kk)
                            gl_all_Cell(8, i1) = gl_RAY_E(i, j - 1 + 1, kk)
                        else if (i == par_n_HP - par_n_TS + 1) then
                            gl_all_Cell(1, i1) = gl_RAY_E(i, j - 1, k)
                            gl_all_Cell(2, i1) = gl_RAY_O(2, j - 1, k)
                            gl_all_Cell(3, i1) = gl_RAY_O(2, j - 1 + 1, k)
                            gl_all_Cell(4, i1) = gl_RAY_E(i, j - 1 + 1, k)
                            gl_all_Cell(5, i1) = gl_RAY_E(i, j - 1, kk)
                            gl_all_Cell(6, i1) = gl_RAY_O(2, j - 1, kk)
                            gl_all_Cell(7, i1) = gl_RAY_O(2, j - 1 + 1, kk)
                            gl_all_Cell(8, i1) = gl_RAY_E(i, j - 1 + 1, kk)
                        else  ! i > par_n_HP - par_n_TS + 1
                            gl_all_Cell(1, i1) = gl_RAY_O(i - (par_n_HP - par_n_TS + 1) + 1, j - 1, k)
                            gl_all_Cell(2, i1) = gl_RAY_O(i - (par_n_HP - par_n_TS + 1) + 2, j - 1, k)
                            gl_all_Cell(3, i1) = gl_RAY_O(i - (par_n_HP - par_n_TS + 1) + 2, j - 1 + 1, k)
                            gl_all_Cell(4, i1) = gl_RAY_O(i - (par_n_HP - par_n_TS + 1) + 1, j - 1 + 1, k)
                            gl_all_Cell(5, i1) = gl_RAY_O(i - (par_n_HP - par_n_TS + 1) + 1, j - 1, kk)
                            gl_all_Cell(6, i1) = gl_RAY_O(i - (par_n_HP - par_n_TS + 1) + 2, j - 1, kk)
                            gl_all_Cell(7, i1) = gl_RAY_O(i - (par_n_HP - par_n_TS + 1) + 2, j - 1 + 1, kk)
                            gl_all_Cell(8, i1) = gl_RAY_O(i - (par_n_HP - par_n_TS + 1) + 1, j - 1 + 1, kk)
                        end if
                    else  ! j == 1
                        if (i < par_n_HP - par_n_TS + 1) then
                            gl_all_Cell(1, i1) = gl_RAY_B(par_n_TS + i - 1, size(gl_RAY_B(par_n_TS + i - 1, :, k)), k)
                            gl_all_Cell(2, i1) = gl_RAY_B(par_n_TS + i, size(gl_RAY_B(par_n_TS + i - 1, :, k)), k)
                            gl_all_Cell(3, i1) = gl_RAY_E(i + 1, j, k)
                            gl_all_Cell(4, i1) = gl_RAY_E(i, j, k)
                            gl_all_Cell(5, i1) = gl_RAY_B(par_n_TS + i - 1, size(gl_RAY_B(par_n_TS + i - 1, :, k)), kk)
                            gl_all_Cell(6, i1) = gl_RAY_B(par_n_TS + i, size(gl_RAY_B(par_n_TS + i - 1, :, k)), kk)
                            gl_all_Cell(7, i1) = gl_RAY_E(i + 1, j, kk)
                            gl_all_Cell(8, i1) = gl_RAY_E(i, j, kk)
                        else if (i == par_n_HP - par_n_TS + 1) then
                            gl_all_Cell(1, i1) = gl_RAY_B(par_n_TS + i - 1, size(gl_RAY_B(par_n_TS + i - 1, :, k)), k)
                            gl_all_Cell(2, i1) = gl_RAY_C(2, size(gl_RAY_C(2, :, k)), k)
                            gl_all_Cell(3, i1) = gl_RAY_O(2, j, k)
                            gl_all_Cell(4, i1) = gl_RAY_E(i, j, k)
                            gl_all_Cell(5, i1) = gl_RAY_B(par_n_TS + i - 1, size(gl_RAY_B(par_n_TS + i - 1, :, k)), kk)
                            gl_all_Cell(6, i1) = gl_RAY_C(2, size(gl_RAY_C(2, :, k)), kk)
                            gl_all_Cell(7, i1) = gl_RAY_O(2, j, kk)
                            gl_all_Cell(8, i1) = gl_RAY_E(i, j, kk)
                        else  ! i > par_n_HP - par_n_TS + 1
                            gl_all_Cell(1, i1) = gl_RAY_C(i - (par_n_HP - par_n_TS + 1) + 1, size(gl_RAY_C(1, :, k)), k)
                            gl_all_Cell(2, i1) = gl_RAY_C(i - (par_n_HP - par_n_TS + 1) + 2, size(gl_RAY_C(1, :, k)), k)
                            gl_all_Cell(3, i1) = gl_RAY_O(i - (par_n_HP - par_n_TS + 1) + 2, j, k)
                            gl_all_Cell(4, i1) = gl_RAY_O(i - (par_n_HP - par_n_TS + 1) + 1, j, k)
                            gl_all_Cell(5, i1) = gl_RAY_C(i - (par_n_HP - par_n_TS + 1) + 1, size(gl_RAY_C(1, :, k)), kk)
                            gl_all_Cell(6, i1) = gl_RAY_C(i - (par_n_HP - par_n_TS + 1) + 2, size(gl_RAY_C(1, :, k)), kk)
                            gl_all_Cell(7, i1) = gl_RAY_O(i - (par_n_HP - par_n_TS + 1) + 2, j, kk)
                            gl_all_Cell(8, i1) = gl_RAY_O(i - (par_n_HP - par_n_TS + 1) + 1, j, kk)
                        end if
                    end if


                    gl_Cell_C(i, j, k) = i1
                    i1 = i1 + 1

                end do
            end do
        end do

        ! Свяжем ячейки с их соседями (имеется строгий порядок)

        ! Начнём с ячеек группы А
        N3 = size(gl_Cell_A(1, 1, :))
        N2 = size(gl_Cell_A(1, :, 1))
        N1 = size(gl_Cell_A(:, 1, 1))

        do k = 1, N3
            do j = 1, N2
                do i = 1, N1

                    kk = k + 1                      ! Слой ячеек выше
                    if (kk > par_l_phi) kk = 1

                    kk2 = k - 1                     ! Слой ячеек ниже
                    if (kk2 < 1) kk2 = par_l_phi


                    ! Заполняем первого соседа (по лучу вверх)
                    if (i == N1) then
                        if (j < par_m_A) then
                            gl_Cell_neighbour(1, gl_Cell_A(i, j, k)) = -1   ! Граница (набегающий поток)
                        else
                            gl_Cell_neighbour(1, gl_Cell_A(i, j, k)) = -3   ! Граница верхний цилиндр
                        end if
                    else
                        gl_Cell_neighbour(1, gl_Cell_A(i, j, k)) = gl_Cell_A(i + 1, j, k)
                    end if

                    ! Заполняем второго соседа (по лучу вниз)
                    if (i == 1) then
                        continue
                    else
                        gl_Cell_neighbour(2, gl_Cell_A(i, j, k)) = gl_Cell_A(i - 1, j, k)
                    end if

                    ! Заполняем третьего соседа (по углу в плоскости, в сторону увеличения)
                    if (j == N2) then
                        if (i < par_n_TS) then
                            gl_Cell_neighbour(3, gl_Cell_A(i, j, k)) = gl_Cell_B(i, par_m_K, k)
                        else
                            gl_Cell_neighbour(3, gl_Cell_A(i, j, k)) = gl_Cell_C(i - par_n_TS + 1, 1, k)
                        end if
                    else
                        gl_Cell_neighbour(3, gl_Cell_A(i, j, k)) = gl_Cell_A(i, j + 1, k)
                    end if

                    ! Заполняем четвёртого соседа (по углу в плоскости, в сторону уменьшения)
                    if (j == 1) then
                        continue
                    else
                        gl_Cell_neighbour(4, gl_Cell_A(i, j, k)) = gl_Cell_A(i, j - 1, k)
                    end if

                    ! Заполняем пятого и шестого соседа (вверх и вниз)
                    gl_Cell_neighbour(5, gl_Cell_A(i, j, k)) = gl_Cell_A(i, j, kk)
                    gl_Cell_neighbour(6, gl_Cell_A(i, j, k)) = gl_Cell_A(i, j, kk2)
                end do
            end do
        end do

        ! Ячееки группы B
        N3 = size(gl_Cell_B(1, 1, :))
        N2 = size(gl_Cell_B(1, :, 1))
        N1 = size(gl_Cell_B(:, 1, 1))

        do k = 1, N3
            do j = 1, N2
                do i = 1, N1

                    kk = k + 1                      ! Слой ячеек выше
                    if (kk > par_l_phi) kk = 1

                    kk2 = k - 1                     ! Слой ячеек ниже
                    if (kk2 < 1) kk2 = par_l_phi


                    ! Заполняем первого соседа (по лучу вверх)
                    if (i == N1) then
                        gl_Cell_neighbour(1, gl_Cell_B(i, j, k)) = -2   ! Выходная граница
                    else
                        gl_Cell_neighbour(1, gl_Cell_B(i, j, k)) = gl_Cell_B(i + 1, j, k)
                    end if

                    ! Заполняем второго соседа (по лучу вниз)
                    if (i == 1) then
                        continue
                    else
                        gl_Cell_neighbour(2, gl_Cell_B(i, j, k)) = gl_Cell_B(i - 1, j, k)
                    end if

                    ! Заполняем третьего соседа (по углу в плоскости наверх на схеме)
                    if (j == N2) then
                        if(i < par_n_TS) then
                            gl_Cell_neighbour(3, gl_Cell_B(i, j, k)) = gl_Cell_A(i, size(gl_Cell_A(i,:, k)), k)
                        else
                            gl_Cell_neighbour(3, gl_Cell_B(i, j, k)) = gl_Cell_C(1, i - par_n_TS + 1, k)
                        end if
                    else
                        gl_Cell_neighbour(3, gl_Cell_B(i, j, k)) = gl_Cell_B(i, j + 1, k)
                    end if

                    ! Заполняем четвёртого соседа (по углу в плоскости)
                    if (j == 1) then
                        continue
                    else
                        gl_Cell_neighbour(4, gl_Cell_B(i, j, k)) = gl_Cell_B(i, j - 1, k)
                    end if

                    ! Заполняем пятого и шестого соседа (вверх и вниз)
                    gl_Cell_neighbour(5, gl_Cell_B(i, j, k)) = gl_Cell_B(i, j, kk)
                    gl_Cell_neighbour(6, gl_Cell_B(i, j, k)) = gl_Cell_B(i, j, kk2)
                end do
            end do
        end do

        ! Ячееки группы C
        N3 = size(gl_Cell_C(1, 1, :))
        N2 = size(gl_Cell_C(1, :, 1))
        N1 = size(gl_Cell_C(:, 1, 1))

        do k = 1, N3
            do j = 1, N2
                do i = 1, N1

                    kk = k + 1                      ! Слой ячеек выше
                    if (kk > par_l_phi) kk = 1

                    kk2 = k - 1                     ! Слой ячеек ниже
                    if (kk2 < 1) kk2 = par_l_phi


                    ! Заполняем первого соседа (по лучу вверх)
                    if (i == N1) then
                        gl_Cell_neighbour(1, gl_Cell_C(i, j, k)) = -3   ! Верхний цилиндр
                    else
                        gl_Cell_neighbour(1, gl_Cell_C(i, j, k)) = gl_Cell_C(i + 1, j, k)
                    end if

                    ! Заполняем второго соседа (по лучу вниз)
                    if (i == 1) then
                        gl_Cell_neighbour(2, gl_Cell_C(i, j, k)) = gl_Cell_B(j + par_n_TS - 1, size(gl_Cell_B(j + par_n_TS - 1, :, k)), k)
                    else
                        gl_Cell_neighbour(2, gl_Cell_C(i, j, k)) = gl_Cell_C(i - 1, j, k)
                    end if

                    ! Заполняем третьего соседа (по углу в плоскости, в сторону увеличения угла!)
                    if (j == N2) then
                        gl_Cell_neighbour(3, gl_Cell_C(i, j, k)) = -2  ! Выходная граница
                    else
                        gl_Cell_neighbour(3, gl_Cell_C(i, j, k)) = gl_Cell_C(i, j + 1, k)
                    end if

                    ! Заполняем четвёртого соседа (по углу в плоскости, в сторону уменьшения)
                    if (j == 1) then
                        gl_Cell_neighbour(4, gl_Cell_C(i, j, k)) = gl_Cell_A(i + par_n_TS - 1, size(gl_Cell_A(i + par_n_TS - 1, : , k)), k)
                    else
                        gl_Cell_neighbour(4, gl_Cell_C(i, j, k)) = gl_Cell_C(i, j - 1, k)
                    end if

                    ! Заполняем пятого и шестого соседа (вверх и вниз)
                    gl_Cell_neighbour(5, gl_Cell_C(i, j, k)) = gl_Cell_C(i, j, kk)
                    gl_Cell_neighbour(6, gl_Cell_C(i, j, k)) = gl_Cell_C(i, j, kk2)
                end do
            end do
        end do

        if (par_developer_info) print*, "Build_Mesh_start: Number all points: ", size(gl_x), " == ", node - 1

        ! Создадим грани

        node = 1

        ! Начнём с ячеек группы А
        N3 = size(gl_Cell_A(1, 1, :))
        N2 = size(gl_Cell_A(1, :, 1))
        N1 = size(gl_Cell_A(:, 1, 1))

        do k = 1, N3
            do j = 1, N2
                do i = 1, N1

                    kk = k + 1                      ! Слой ячеек выше
                    if (kk > par_l_phi) kk = 1

                    ni = gl_Cell_A(i, j, k)

                    ! Заполняем 1-ую грань для этой ячейки
                    if (i < N1) then
                        gl_all_Gran(1, node) = gl_all_Cell(2, ni)
                        gl_all_Gran(2, node) = gl_all_Cell(3, ni)
                        gl_all_Gran(3, node) = gl_all_Cell(7, ni)
                        gl_all_Gran(4, node) = gl_all_Cell(6, ni)

                        gl_Gran_neighbour(1, node) = ni
                        gl_Gran_neighbour(2, node) = gl_Cell_A(i + 1, j, k)

                        gl_Cell_gran(1, ni) = node
                        gl_Cell_gran(2, gl_Cell_A(i + 1, j, k)) = node

                        node = node + 1
                    else
                        gl_all_Gran(1, node) = gl_all_Cell(2, ni)
                        gl_all_Gran(2, node) = gl_all_Cell(3, ni)
                        gl_all_Gran(3, node) = gl_all_Cell(7, ni)
                        gl_all_Gran(4, node) = gl_all_Cell(6, ni)

                        gl_Gran_neighbour(1, node) = ni

                        if (j < par_m_A) then
                            gl_Gran_neighbour(2, node) = -1
                        else
                            gl_Gran_neighbour(2, node) = -3
                        end if

                        gl_Cell_gran(1, ni) = node
                        node = node + 1
                    end if

                    ! Заполняем 4-ую грань для этой ячейки
                    if (j /= 1) then
                        gl_all_Gran(1, node) = gl_all_Cell(1, ni)
                        gl_all_Gran(2, node) = gl_all_Cell(2, ni)
                        gl_all_Gran(3, node) = gl_all_Cell(6, ni)
                        gl_all_Gran(4, node) = gl_all_Cell(5, ni)

                        gl_Gran_neighbour(1, node) = ni
                        gl_Gran_neighbour(2, node) = gl_Cell_A(i, j - 1, k)

                        gl_Cell_gran(4, ni) = node
                        gl_Cell_gran(3, gl_Cell_A(i, j - 1, k)) = node

                        node = node + 1
                    end if

                    ! Заполняем 5-ую грань для этой ячейки
                    gl_all_Gran(1, node) = gl_all_Cell(5, ni)
                    gl_all_Gran(2, node) = gl_all_Cell(6, ni)
                    gl_all_Gran(3, node) = gl_all_Cell(7, ni)
                    gl_all_Gran(4, node) = gl_all_Cell(8, ni)

                    gl_Gran_neighbour(1, node) = ni
                    gl_Gran_neighbour(2, node) = gl_Cell_A(i, j, kk)

                    gl_Cell_gran(5, ni) = node
                    gl_Cell_gran(6, gl_Cell_A(i, j, kk)) = node

                    node = node + 1
                end do
            end do
        end do

        ! Ячееки группы B
        N3 = size(gl_Cell_B(1, 1, :))
        N2 = size(gl_Cell_B(1, :, 1))
        N1 = size(gl_Cell_B(:, 1, 1))

        do k = 1, N3
            do j = 1, N2
                do i = 1, N1

                    kk = k + 1                      ! Слой ячеек выше
                    if (kk > par_l_phi) kk = 1

                    ni = gl_Cell_B(i, j, k)

                    ! Заполняем 1-ую грань для этой ячейки
                    if (i < N1) then
                        gl_all_Gran(1, node) = gl_all_Cell(3, ni)
                        gl_all_Gran(2, node) = gl_all_Cell(4, ni)
                        gl_all_Gran(3, node) = gl_all_Cell(8, ni)
                        gl_all_Gran(4, node) = gl_all_Cell(7, ni)

                        gl_Gran_neighbour(1, node) = ni
                        gl_Gran_neighbour(2, node) = gl_Cell_B(i + 1, j, k)

                        gl_Cell_gran(1, ni) = node
                        gl_Cell_gran(2, gl_Cell_B(i + 1, j, k)) = node

                        node = node + 1
                    else
                        gl_all_Gran(1, node) = gl_all_Cell(3, ni)
                        gl_all_Gran(2, node) = gl_all_Cell(4, ni)
                        gl_all_Gran(3, node) = gl_all_Cell(8, ni)
                        gl_all_Gran(4, node) = gl_all_Cell(7, ni)

                        gl_Gran_neighbour(1, node) = ni
                        gl_Gran_neighbour(2, node) = -2

                        gl_Cell_gran(1, ni) = node

                        node = node + 1
                    end if

                    ! Заполняем 4-ую грань для этой ячейки
                    if (j /= 1) then
                        gl_all_Gran(1, node) = gl_all_Cell(4, ni)
                        gl_all_Gran(2, node) = gl_all_Cell(1, ni)
                        gl_all_Gran(3, node) = gl_all_Cell(5, ni)
                        gl_all_Gran(4, node) = gl_all_Cell(8, ni)

                        gl_Gran_neighbour(1, node) = ni
                        gl_Gran_neighbour(2, node) = gl_Cell_B(i, j - 1, k)

                        gl_Cell_gran(4, ni) = node
                        gl_Cell_gran(3, gl_Cell_B(i, j - 1, k)) = node

                        node = node + 1
                    end if

                    ! Заполняем 5-ую грань для этой ячейки
                    gl_all_Gran(1, node) = gl_all_Cell(5, ni)
                    gl_all_Gran(2, node) = gl_all_Cell(6, ni)
                    gl_all_Gran(3, node) = gl_all_Cell(7, ni)
                    gl_all_Gran(4, node) = gl_all_Cell(8, ni)

                    gl_Gran_neighbour(1, node) = ni
                    gl_Gran_neighbour(2, node) = gl_Cell_B(i, j, kk)

                    gl_Cell_gran(5, ni) = node
                    gl_Cell_gran(6, gl_Cell_B(i, j, kk)) = node

                    node = node + 1

                    if (j == N2) then  ! 3-яя грань для границы
                        gl_all_Gran(1, node) = gl_all_Cell(2, ni)
                        gl_all_Gran(2, node) = gl_all_Cell(3, ni)
                        gl_all_Gran(3, node) = gl_all_Cell(7, ni)
                        gl_all_Gran(4, node) = gl_all_Cell(6, ni)

                        gl_Gran_neighbour(1, node) = ni

                        if(i < par_n_TS) then
                            gl_Gran_neighbour(2, node) = gl_Cell_A(i, size(gl_Cell_A(i,:, k)), k)
                            gl_Cell_gran(3, gl_Cell_A(i, size(gl_Cell_A(i,:, k)), k)) = node
                        else
                            gl_Gran_neighbour(2, node) = gl_Cell_C(1, i - par_n_TS + 1, k)
                            gl_Cell_gran(2, gl_Cell_C(1, i - par_n_TS + 1, k)) = node
                        end if

                        gl_Cell_gran(3, ni) = node

                        node = node + 1
                    end if


                end do
            end do
        end do

        ! Ячейки группы C
        N3 = size(gl_Cell_C(1, 1, :))
        N2 = size(gl_Cell_C(1, :, 1))
        N1 = size(gl_Cell_C(:, 1, 1))

        do k = 1, N3
            do j = 1, N2
                do i = 1, N1

                    kk = k + 1                      ! Слой ячеек выше
                    if (kk > par_l_phi) kk = 1

                    ni = gl_Cell_C(i, j, k)

                    ! Заполняем 1-ую грань для этой ячейки
                    if (i < N1) then
                        gl_all_Gran(1, node) = gl_all_Cell(2, ni)
                        gl_all_Gran(2, node) = gl_all_Cell(3, ni)
                        gl_all_Gran(3, node) = gl_all_Cell(7, ni)
                        gl_all_Gran(4, node) = gl_all_Cell(6, ni)

                        gl_Gran_neighbour(1, node) = ni
                        gl_Gran_neighbour(2, node) = gl_Cell_C(i + 1, j, k)

                        gl_Cell_gran(1, ni) = node
                        gl_Cell_gran(2, gl_Cell_C(i + 1, j, k)) = node

                        node = node + 1
                    else
                        gl_all_Gran(1, node) = gl_all_Cell(2, ni)
                        gl_all_Gran(2, node) = gl_all_Cell(3, ni)
                        gl_all_Gran(3, node) = gl_all_Cell(7, ni)
                        gl_all_Gran(4, node) = gl_all_Cell(6, ni)

                        gl_Gran_neighbour(1, node) = ni
                        gl_Gran_neighbour(2, node) = -3    ! Верхний цилиндр

                        gl_Cell_gran(1, ni) = node

                        node = node + 1
                    end if

                    ! Заполняем 4-ую грань для этой ячейки
                    if (j /= 1) then
                        gl_all_Gran(1, node) = gl_all_Cell(1, ni)
                        gl_all_Gran(2, node) = gl_all_Cell(2, ni)
                        gl_all_Gran(3, node) = gl_all_Cell(6, ni)
                        gl_all_Gran(4, node) = gl_all_Cell(5, ni)

                        gl_Gran_neighbour(1, node) = ni
                        gl_Gran_neighbour(2, node) = gl_Cell_C(i, j - 1, k)

                        gl_Cell_gran(4, ni) = node
                        gl_Cell_gran(3, gl_Cell_C(i, j - 1, k)) = node

                        node = node + 1
                    end if

                    if (j == N2) then  ! Третья грань на границе
                        gl_all_Gran(1, node) = gl_all_Cell(3, ni)
                        gl_all_Gran(2, node) = gl_all_Cell(4, ni)
                        gl_all_Gran(3, node) = gl_all_Cell(8, ni)
                        gl_all_Gran(4, node) = gl_all_Cell(7, ni)

                        gl_Gran_neighbour(1, node) = ni
                        gl_Gran_neighbour(2, node) = -2

                        gl_Cell_gran(3, ni) = node
                        node = node + 1
                    end if


                    ! Заполняем 5-ую грань для этой ячейки
                    gl_all_Gran(1, node) = gl_all_Cell(5, ni)
                    gl_all_Gran(2, node) = gl_all_Cell(6, ni)
                    gl_all_Gran(3, node) = gl_all_Cell(7, ni)
                    gl_all_Gran(4, node) = gl_all_Cell(8, ni)

                    gl_Gran_neighbour(1, node) = ni
                    gl_Gran_neighbour(2, node) = gl_Cell_C(i, j, kk)

                    gl_Cell_gran(5, ni) = node
                    gl_Cell_gran(6, gl_Cell_C(i, j, kk)) = node

                    node = node + 1

                    if (j == 1) then  ! 4-яя грань для границы
                        gl_all_Gran(1, node) = gl_all_Cell(1, ni)
                        gl_all_Gran(2, node) = gl_all_Cell(2, ni)
                        gl_all_Gran(3, node) = gl_all_Cell(6, ni)
                        gl_all_Gran(4, node) = gl_all_Cell(5, ni)

                        gl_Gran_neighbour(1, node) = ni

                        gl_Gran_neighbour(2, node) = gl_Cell_A(i + par_n_TS - 1, size(gl_Cell_A(i + par_n_TS - 1, : , k)), k)
                        gl_Cell_gran(3, gl_Cell_A(i + par_n_TS - 1, size(gl_Cell_A(i + par_n_TS - 1, : , k)), k)) = node

                        gl_Cell_gran(4, ni) = node

                        node = node + 1
                    end if

                end do
            end do
        end do

        
        ! Теперь для каждой ячейки добавим ей информацию какого она типа и какой у неё номер (тройка) в этом типе
        ! Начнём с группы А
        N3 = size(gl_Cell_A(1, 1, :))
        N2 = size(gl_Cell_A(1, :, 1))
        N1 = size(gl_Cell_A(:, 1, 1))

        do k = 1, N3
            do j = 1, N2
                do i = 1, N1
                    num1 = gl_Cell_A(i, j, k)
                    gl_Cell_type(num1) = "A"
                    gl_Cell_number(:, num1) = (/ i, j, k /)
                end do
            end do
        end do
        
        ! Начнём с группы B
        N3 = size(gl_Cell_B(1, 1, :))
        N2 = size(gl_Cell_B(1, :, 1))
        N1 = size(gl_Cell_B(:, 1, 1))

        do k = 1, N3
            do j = 1, N2
                do i = 1, N1
                    num1 = gl_Cell_B(i, j, k)
                    gl_Cell_type(num1) = "B"
                    gl_Cell_number(:, num1) = (/ i, j, k /)
                end do
            end do
        end do
        
        ! Начнём с группы C
        N3 = size(gl_Cell_C(1, 1, :))
        N2 = size(gl_Cell_C(1, :, 1))
        N1 = size(gl_Cell_C(:, 1, 1))

        do k = 1, N3
            do j = 1, N2
                do i = 1, N1
                    num1 = gl_Cell_C(i, j, k)
                    gl_Cell_type(num1) = "C"
                    gl_Cell_number(:, num1) = (/ i, j, k /)
                end do
            end do
        end do
        
        

        if (par_developer_info) print *, "Build_Mesh_start: All Cell = ", i1 - 1

        if (par_developer_info) print*, "Build_Mesh_start: Number all gran: ", size(gl_all_Gran(1,:)), " == ", node - 1

        !print*, par_n_END * par_l_phi * (par_m_A + par_m_BC) + par_m_K * (par_n_TS + par_m_O) * par_l_phi + par_l_phi * (par_n_END - par_n_TS + 1) * par_m_O - &
        !        (par_m_A + par_m_BC + par_m_K - 1) * par_l_phi - par_n_END * (par_l_phi - 1) - (par_n_TS + par_m_O - 1) * (par_l_phi - 1)
    end subroutine Build_Mesh_start


    subroutine Find_inner()  ! Находит внутренние ячейки и грани (действует после calc_all_Gran)
        use STORAGE
        use GEO_PARAM
        implicit none

        integer :: n, Ngran, iter
        
        if(allocated(gl_all_Gran_inner)) deallocate(gl_all_Gran_inner)
        if(allocated(gl_all_Cell_inner)) deallocate(gl_all_Cell_inner)

        n = 0
        Ngran = size(gl_all_Gran(1,:))

        do  iter = 1, Ngran
            if (gl_Gran_info(iter) /= 0) n = n + 1
        end do

        allocate(gl_all_Gran_inner(n))
        n = 1

        do  iter = 1, Ngran
            if (gl_Gran_info(iter) /= 0) then
                gl_all_Gran_inner(n) = iter
                n = n + 1
            end if
        end do


        n = 0
        Ngran = size(gl_all_Cell(1,:))

        do  iter = 1, Ngran
            if (gl_Cell_info(iter) /= 2) n = n + 1
        end do

        allocate(gl_all_Cell_inner(n))
        n = 1

        do  iter = 1, Ngran
            if (gl_Cell_info(iter) /= 2) then
                gl_all_Cell_inner(n) = iter
                n = n + 1
            end if
        end do
	end subroutine Find_inner
    ! ************************************************************************************************************************************************
    !! Блок физики
	
	subroutine Test_raspadnik()
        ! Variables
        USE Solvers
        real(8) :: dsl, dsp, dsc
        real(8) :: al, be, ge, w
        integer(4) :: n_state
        logical :: null_bn1, p_correct_
        integer :: n_disc
        real(8) :: konvect_(3, 1)
        real(8) :: qqq1(9), qqq2(9), qqq(9)
        
        null_bn1 = .False.; p_correct_ = .True.; n_disc = 1
        al = 1.0; be = 0.0; ge = 0; w = 0.0; n_state = 3
        
        qqq1 = (/10.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 10.5, 100.0/)
        qqq2 = (/1.0, 0.0, 0.0, 0.0, 6.0, 0.0, 0.0, 10.5, 100.0/)
        konvect_(1, 1) = 2.0
        konvect_(2, 1) = 100.0
        
        call chlld_Q(n_state, al, be, ge, &
                                    w, qqq1, qqq2, &
                                    dsl, dsp, dsc, &
                                    qqq, null_bn1 = null_bn1, n_disc = n_disc, p_correct_ = p_correct_, konvect_ = konvect_)
        
        print*, "qqq = ", qqq
        print*, "konvect_ = ", konvect_
        print*, "____"
        w = dsc
        call chlld_Q(n_state, al, be, ge, &
                                    w, qqq1, qqq2, &
                                    dsl, dsp, dsc, &
                                    qqq, null_bn1 = null_bn1, n_disc = n_disc, p_correct_ = p_correct_, konvect_ = konvect_)
        
        print*, "qqq = ", qqq
        print*, "konvect_ = ", konvect_
        pause
	end subroutine Test_raspadnik
	
	subroutine Test_koordinate(x, y, z)
	    use STORAGE
        use GEO_PARAM
        implicit none
		
		real(8), intent(in) :: x, y, z
		
        real(8) :: c(3), Matr(3, 3), Matr2(3, 3), cc(3), vv(3), r, the, tt, ro
	    integer :: ijk
		
		Matr(1,1) = -0.995868
        Matr(1,2) = 0.0177307
        Matr(1,3) = 0.0890681
        Matr(2,1) = 0.0730412
        Matr(2,2) = 0.739193
        Matr(2,3) = 0.669521
        Matr(3,1) = -0.0539675
        Matr(3,2) = 0.67326
        Matr(3,3) = -0.737434
	
	    Matr2(1,1) = -0.995868
        Matr2(1,2) = 0.0730412
        Matr2(1,3) = -0.0539675
        Matr2(2,1) = 0.0177307
        Matr2(2,2) = 0.739193
        Matr2(2,3) = 0.67326
        Matr2(3,1) = 0.0890681
        Matr2(3,2) = 0.669521
        Matr2(3,3) = -0.737434
		
		!print*, "pr = ", MATMUL(Matr, (/0, 0, 1/))
		
		c = (/x, y, z/)
		
		r = norm2(c)
	
	    cc = MATMUL(Matr2, c)
		print*, "cc = ", cc
	    the = acos(cc(3)/r)
		the = -the + par_pi_8/2.0   ! Т.к. в данных по СВ на 1 а.е. угол от -90 до 90 у Алексашова
		
		print*, "the = ", the * 180/par_pi_8
		
		tt = ((the + par_pi_8/2.0)/(par_pi_8/18.0) - FLOOR((the + par_pi_8/2.0)/(par_pi_8/18.0)))
	    ijk = INT( (the + par_pi_8/2.0)/(par_pi_8/18.0) )
	    ro = (par_density_in(ijk + 1) * (1.0 - tt) + tt * par_density_in(ijk + 2))/ 0.04
		
		c = DBLE(c) * (par_velocity_in(ijk + 1) * (1.0 - tt) + tt * par_velocity_in(ijk + 2))/DBLE(r) * 0.0963179
		
		print*, ro * 0.04
		print*, c(1) * 10.3823, c(2) * 10.3823, c(3) * 10.3823
	
	
	end subroutine Test_koordinate
    
    subroutine Inner_conditions(ncell)  ! Задаём начальные условия в ячейке (функция принимает номер ячейки, в которой надо задать эти условия)
        use STORAGE
        use GEO_PARAM
        implicit none
        
        integer(4), intent(in) :: ncell
        real(8) :: r, tt
        real(8) :: ro, P_E, the, BPHI, BR, V1, V2, V3
        real(8) :: c(3), Matr(3, 3), Matr2(3, 3), cc(3), vv(3)
        integer :: ijk

        
        Matr(1,1) = -0.9958639688067077
        Matr(1,2) = 0.01776569097515556
        Matr(1,3) = 0.08910295088675518
        Matr(2,1) = 0.07561695085992419
        Matr(2,2) = 0.7057402284561812
        Matr(2,3) = 0.7044237408557894
        Matr(3,1) = -0.05036896241933166
        Matr(3,2) = 0.7082479157489926
        Matr(3,3) = -0.7041646522383864
        
        Matr2(1,1) = -0.9958639688067080
        Matr2(1,2) = 0.0756169508599243
        Matr2(1,3) = -0.0503689624193315
        Matr2(2,1) = 0.0177656909751554
        Matr2(2,2) = 0.7057402284561816
        Matr2(2,3) = 0.7082479157489927
        Matr2(3,1) = 0.0891029508867553
        Matr2(3,2) = 0.7044237408557898
        Matr2(3,3) = -0.7041646522383865
        
        !   Matr(1,1) = -0.9958639688067080
        !   Matr(1,2) = 0.0177656909751554
        !   Matr(1,3) = 0.0891029508867553
        !   Matr(2,1) = 0.0756169508599243
        !   Matr(2,2) = 0.7057402284561816
        !   Matr(2,3) = 0.7044237408557898
        !   Matr(3,1) = -0.0503689624193315
        !   Matr(3,2) = 0.7082479157489927
        !   Matr(3,3) = -0.7041646522383864
        !
        !   Matr2(1,1) = -0.99586397
        !   Matr2(1,2) = 0.07561695
        !   Matr2(1,3) = -0.05036896
        !   Matr2(2,1) = 0.01776569
        !   Matr2(2,2) = 0.70574023
        !   Matr2(2,3) = 0.70824792
        !   Matr2(3,1) = 0.08910295
        !   Matr2(3,2) = 0.70442374
        !   Matr2(3,3) = -0.70416465
        
        
        c = gl_Cell_center(:, ncell)
        r = norm2(c)
        
        
        cc = MATMUL(Matr2, c)
        
        the = acos(cc(3)/r)
        
        !cc(1) = 1.0
        !cc(2) = -1.0
        !cc(3) = 1.0
        !r = norm2(cc)
        
        !print*, acos(1.0), acos(0.5), acos(0.0), acos(-1.0)
        !print*, cc
        !Pause
        
        !! НАДО ТУТ ПЕРЕДЕЛЫВАТЬ ОПРЕДЕЛЕНИЕ МАГНИТНОГО ПОЛЯ
        ! BE = sqrt(cpi4 * par_kk)/(par_Mach_alf * r)

        BR = -par_B_0 * (par_R0/r)**2
        BPHI = -BR * sin(the) * (r/par_R0)
        
        
        !print*, "Be = ", BE, sin(the), the
        
        call dekard_skorost(cc(3), cc(1), cc(2), BR, BPHI, 0.0_8, V3, V1, V2)
        
        
        vv(1) = V1
        vv(2) = V2
        vv(3) = V3
        
        ! print*, vv
        !PAUSE
        
        the = -the + par_pi_8/2.0   ! Т.к. в данных по СВ на 1 а.е. угол от -90 до 90 у Алексашова
        
        cc = MATMUL(Matr, vv)
        
        tt = ((the + par_pi_8/2.0)/(par_pi_8/18.0) - FLOOR((the + par_pi_8/2.0)/(par_pi_8/18.0)))
        ijk = INT( (the + par_pi_8/2.0)/(par_pi_8/18.0) )
        ro = (par_density_in(ijk + 1) * (1.0 - tt) + tt * par_density_in(ijk + 2))/ par_n_p_inf
        
        !print*, the + par_pi_8/2.0, (the + par_pi_8/2.0)/(par_pi_8/18.0), tt, INT( (the + par_pi_8/2.0)/(par_pi_8/18.0) )
        
        !ro =  par_kk/(par_chi**2 * r**2)
        !P_E = (par_kk/(par_chi**2 * par_R0**2)) * par_chi**2 / (ggg * par_Mach_0**2)
        
        !c = DBLE(c) * par_chi/DBLE(r)
        
        c = DBLE(c) * (par_velocity_in(ijk + 1) * (1.0 - tt) + tt * par_velocity_in(ijk + 2))/DBLE(r) / par_V_character
        P_E = par_p_0 ! (ro * norm2(c)**2) / (ggg * par_Mach_0**2)
        
        ! Задаём плазму  (ro, u, v, w, p, bx, by, bz, Q)
        !gl_Cell_par(:, ncell) = (/ro, DBLE(c(1)), DBLE(c(2)), DBLE(c(3)), P_E * (par_R0/r) ** (2.0 * ggg), cc(1), cc(2), cc(3), ro/)
        gl_Cell_par(:, ncell) = (/ro * (par_R0/r)**2, DBLE(c(1)), DBLE(c(2)), DBLE(c(3)), & 
            P_E * (par_R0/r) ** (2.0 * ggg), -cc(1), -cc(2), -cc(3), ro * (par_R0/r)**2/)
        
        ! Если включаем гелий
        if(par_helium) then	
            gl_Cell_par2(1, ncell) = gl_Cell_par(1, ncell) * par_mHe_0   ! He++
            
            gl_Cell_par(1, ncell) = gl_Cell_par(1, ncell) * (1.0_8 + par_mHe_0)
            gl_Cell_par(9, ncell) = gl_Cell_par(9, ncell) * (1.0_8 + par_mHe_0)
            !gl_Cell_par(5, ncell) = gl_Cell_par(5, ncell) * 1.0525
        end if
                
        ! Задаём мультифлюид (только первый сорт)
        gl_Cell_par_MF(:, 1, ncell) = (/0.0000001_8, c(1), c(2), c(3), 0.0000001_8 * P_E/ro /)
        
        gl_Cell_par_MF(1, 2, ncell) = 0.0000001_8
        gl_Cell_par_MF(5, 2, ncell) = 0.00001_8
    end subroutine Inner_conditions

    subroutine Initial_conditions()  ! Задаём начальные условия
        use STORAGE
        use GEO_PARAM
        implicit none

        integer(4) :: ncell, gr, N1, N2, N3, j, k, ncell2
        real(8) :: qqq(9), dist, dist2
        
        
        ! Перенормировка параметров, если это необходимо
        ncell = size(gl_all_Cell(1, :))
        
        do gr = 1, ncell
            
            !if(gl_zone_Cell(gr) == 3 .or. gl_zone_Cell(gr) == 4) then
            !if(gl_Cell_center(1, gr) > 0.0 .and. norm2(gl_Cell_center(:, gr)) > 200.0) then
            !	gl_Cell_par(1, gr) = 1.3
            !	gl_Cell_par2(1, gr) = 0.3
            !	gl_Cell_par(9, gr) = 1.3
            !end if
            
            !else
            !	gl_Cell_par2(1, gr) = gl_Cell_par2(1, gr) /0.07 * 0.14
            !	gl_Cell_par(1, gr) = gl_Cell_par(1, gr) /1.07 * 1.14
            !	gl_Cell_par(9, gr) = gl_Cell_par(9, gr) /1.07 * 1.14
            !end if
            
        end do
        

        
        ! Задаём граничные условия (параметры в первых ячейках на внутренней сфере)
        
        N3 = size(gl_Cell_A(1, 1, :))
        N2 = size(gl_Cell_A(1, :, 1))
        N1 = size(gl_Cell_A(:, 1, 1))
        
 
        do k = 1, N3
            do j = 1, N2
                ncell = gl_Cell_A(1, j, k)
                call Inner_conditions(ncell)
                
                ncell = gl_Cell_A(2, j, k)
                call Inner_conditions(ncell)
                
            end do
        end do
        
        
        N3 = size(gl_Cell_B(1, 1, :))
        N2 = size(gl_Cell_B(1, :, 1))
        N1 = size(gl_Cell_B(:, 1, 1))
        
        do k = 1, N3
            do j = 1, N2
                ncell = gl_Cell_B(1, j, k)
                call Inner_conditions(ncell)
                
                ncell = gl_Cell_B(2, j, k)
                call Inner_conditions(ncell)
                
            end do
        end do
        
        
        ! Пробежимся по всем ячейкам и зададим какие-то условия, если нужно

        ncell = size(gl_all_Cell(1, :))
        do gr = 1, ncell
            
            if(gl_zone_Cell(gr) == 1 .or. gl_zone_Cell(gr) == 2) then
            ! gl_Cell_par(6:8, gr) = -gl_Cell_par(6:8, gr)
            end if
            !if(gl_Cell_center(1, gr) > 0.0 .and. norm2(gl_Cell_center(:, gr)) > 200.0) then
            !	gl_Cell_par(1, gr) = 1.3
            !	gl_Cell_par2(1, gr) = 0.3
            !	gl_Cell_par(9, gr) = 1.3
            !end if
            
            !else
            !	gl_Cell_par2(1, gr) = gl_Cell_par2(1, gr) /0.07 * 0.14
            !	gl_Cell_par(1, gr) = gl_Cell_par(1, gr) /1.07 * 1.14
            !	gl_Cell_par(9, gr) = gl_Cell_par(9, gr) /1.07 * 1.14
            !end if
            
        end do

    end subroutine Initial_conditions

    subroutine calc_all_Gran()   ! Программа расчёта площадей и нормалей граней и Объёмов ячеек
        ! Эта функция нужна только в неподвижной программе (т.к. в подвижном случае быстрее будет это всё считать "налету"
        ! Либо её можно использовать в проверке геометрии сетки
        use STORAGE
        use GEO_PARAM
        use My_func
        implicit none

        integer :: Ngran, iter, Ncell
        real(8) :: p(3, 4), Vol, D, m(3, 4)
        real(8) :: a(3), b(3), c(3), S, node1(3), node2(3)
        real(8) :: dist, di, gr_center(3), ger1, ger2, ger3, ger
        integer :: i, j, k, ll, grc, N1, N2, N3
        real(8) :: pp(3, 8)

        ! Цикл по граням - считаем площадь грани, её нормаль
        Ngran = size(gl_all_Gran(1,:))

        do  iter = 1, Ngran
            grc = 0
            p = 0.0
            a = 0.0
            b = 0.0
            c = 0.0
            gr_center = 0.0
            node1 = 0.0
            node2 = 0.0
            ! Считываем из глобальной памяти координаты точек грани

            !print*, gl_all_Gran(:, iter)
            !pause

            if (gl_all_Gran(1, iter) == gl_all_Gran(3, iter)) then
                print*, "EROROR nierhfue 1", iter
                pause
            end if

            if (gl_all_Gran(2, iter) == gl_all_Gran(4, iter)) then
                print*, "EROROR nierhfue 2", iter
                pause
            end if


            do j = 1, 4
                i = gl_all_Gran(j, iter)
                p(:,j) = (/gl_x(i), gl_y(i), gl_z(i)/)
            end do
            a = p(:,3) - p(:,1)
            b = p(:,4) - p(:,2)

            ! Нужно сохранить центр грани
            gr_center = p(:, 1) + p(:, 2) + p(:, 3) + p(:, 4)
            grc = 4
            !gr_center = p(:, 1)
            !grc = 1
            !if (gl_all_Gran(2, iter) /= gl_all_Gran(1, iter)) then
            !    grc = grc + 1
            !    gr_center = gr_center + p(:,2)
            !end if
            !
            !if (gl_all_Gran(3, iter) /= gl_all_Gran(2, iter) .and. gl_all_Gran(3, iter) /= gl_all_Gran(1, iter)) then
            !    grc = grc + 1
            !    gr_center = gr_center + p(:,3)
            !end if
            !
            !if (gl_all_Gran(4, iter) /= gl_all_Gran(1, iter) .and. gl_all_Gran(4, iter) /= gl_all_Gran(2, iter) &
            !    .and. gl_all_Gran(4, iter) /= gl_all_Gran(3, iter)) then
            !    grc = grc + 1
            !    gr_center = gr_center + p(:,4)
            !end if

            gr_center = gr_center/grc

            gl_Gran_center(:, iter) = gr_center

            
            c(1) = a(2) * b(3) - a(3) * b(2)
            c(2) = a(3) * b(1) - a(1) * b(3)
            c(3) = a(1) * b(2) - a(2) * b(1)

            S = norm2(c)  ! S = S/2
            c = c/S
            S = S/2.0
            gl_Gran_square(iter) = S
            
            ! Альтернативное вычисление площади грани:
            if (.False.) then
            S = 0.0
            ger1 = norm2(p(:,1) - p(:,2))
            ger2 = norm2(p(:,3) - p(:,2))
            ger3 = norm2(p(:,1) - p(:,3))
            ger = (ger1 + ger2 + ger3) / 2.0
            ! вычисление площади по формуле Герона
            S = sqrt(ger * (ger - ger1) * (ger - ger2) * (ger - ger3))
            
            ger1 = norm2(p(:,1) - p(:,4))
            ger2 = norm2(p(:,3) - p(:,4))
            ger = (ger1 + ger2 + ger3) / 2.0
            ! вычисление площади по формуле Герона
            S = S + sqrt(ger * (ger - ger1) * (ger - ger2) * (ger - ger3))
            
            if((S - gl_Gran_square(iter))/S * 100 > 1E-3) then
                print*, "1516 GGFGYTFYG  ",  S, gl_Gran_square(iter), (S - gl_Gran_square(iter))/S * 100
                print*, gr_center
                PAUSE
            end if
            end if
            
            ! if (grc == 4) then
            ! 	m(:, 1) = (p(:, 1) + p(:, 2))/2.0
            ! 	m(:, 2) = (p(:, 2) + p(:, 3))/2.0
            ! 	m(:, 3) = (p(:, 3) + p(:, 4))/2.0
            ! 	m(:, 4) = (p(:, 4) + p(:, 1))/2.0
            ! 	a = m(:,3) - m(:,1)
            !  b = m(:,4) - m(:,2)
            ! 	c(1) = a(2) * b(3) - a(3) * b(2)
            !  c(2) = a(3) * b(1) - a(1) * b(3)
            !  c(3) = a(1) * b(2) - a(2) * b(1)
            ! 	S = norm2(c)  ! S = S/2
            !  c = c/S
            ! end if
            
            gl_Gran_normal(:, iter) = c
            
            ! Можно один раз проверить, правильно ли ориентирована нормаль!\
            
            

            ! call Get_center(gl_Gran_neighbour(1, iter), node1)
            ! node2 = (p(:,1) + p(:,2) + p(:,3) + p(:,4))/4.0
            ! node1 = node2 - node1
            
            
            ! if(DOT_PRODUCT(node1, c) < 0.0) then
            !         call Get_center(gl_Gran_neighbour(1, iter), node1)
            !     print*, "ERROR 1333 ", node1
            !     pause
            ! end if

            !     Нужно записать площадь грани и нормаль в общий массив!
                

            !     if (S < 0.000001) then
            !        print *, "ERROR 134443   gran = 0 ", c, a, b
            !        print *, p(:,1), p(:,2), p(:,3), p(:,4)
            !        pause
            !     end if


            !     print *, gl_Gran_square(iter)
            !     pause
        end do




        ! Теперь посчитаем объёмы ячеек
        Ngran = size(gl_all_Cell(1,:))

        do  iter = 1, Ngran   ! Пробегаемся по всем ячейкам
            Vol = 0.0
            ll = 0
            dist = 10000.0 * par_R_character
            c = 0.0
            ! Для вычисления правильного центра ячейки, нужно понять какая она (некоторые узлы пропущены)
            if (gl_all_Cell(1, iter) /= gl_all_Cell(5, iter)) then
                do j = 1, 8
                    i = gl_all_Cell(j, iter)
                    c = c + (/gl_x(i), gl_y(i), gl_z(i)/)
                end do
                c = c/8.0
            else
                if (gl_all_Cell(1, iter) == gl_all_Cell(2, iter)) then
                    if (gl_all_Cell(4, iter) == gl_all_Cell(8, iter)) then ! Класс ячейки = 4
                        do j = 1,8
                            ! if (j == 2 .or. j == 5 .or. j == 6 .or. j == 8) CYCLE 
                            i = gl_all_Cell(j, iter)
                            c = c + (/gl_x(i), gl_y(i), gl_z(i)/)
                        end do
                        c = c/8.0
                    else  ! Класс ячейки = 6
                        do j = 1,8
                            ! if (j == 2 .or. j == 5 .or. j == 6) CYCLE
                            i = gl_all_Cell(j, iter)
                            c = c + (/gl_x(i), gl_y(i), gl_z(i)/)
                        end do
                        c = c/8.0
                    end if
                else
                    if (gl_all_Cell(2, iter) == gl_all_Cell(6, iter)) then
                        if (gl_all_Cell(4, iter) == gl_all_Cell(5, iter)) then  ! Класс ячейки = 1
                            do j = 1,8
                                ! if (j == 4 .or. j == 5 .or. j == 6 .or. j == 8) CYCLE
                                i = gl_all_Cell(j, iter)
                                c = c + (/gl_x(i), gl_y(i), gl_z(i)/)
                            end do
                            c = c/8.0
                        else  ! Класс ячейки = 2
                            do j = 1,8
                                ! if (j == 5 .or. j == 6) CYCLE
                                i = gl_all_Cell(j, iter)
                                c = c + (/gl_x(i), gl_y(i), gl_z(i)/)
                            end do
                            c = c/8.0   ! 6  ----------------------------------------------------------------------------
                        end if
                    else
                        if (gl_all_Cell(4, iter) == gl_all_Cell(5, iter)) then ! Класс ячейки = 3
                            do j = 1,8
                                !if (j == 4 .or. j == 5 .or. j == 8) CYCLE
                                i = gl_all_Cell(j, iter)
                                c = c + (/gl_x(i), gl_y(i), gl_z(i)/)
                            end do
                            c = c/8.0
                        else  ! Класс ячейки = 5
                            do j = 1,8
                                ! if (j == 5 .or. j == 8) CYCLE        ! ----------------------------------------------------------------------------
                                i = gl_all_Cell(j, iter)
                                c = c + (/gl_x(i), gl_y(i), gl_z(i)/)
                            end do
                            c = c/8.0                            ! 6  ----------------------------------------------------------------------------
                        end if
                    end if
                end if
            end if

            gl_Cell_center(:, iter) = c

            !print *, c
            !pause

            ! Вычислили центр ячейки теперь считаем объём пирамиды на каждую грань
            do j = 1, 6
                i = gl_Cell_gran(j, iter)   ! Берём по очереди все грани ячейки

                !if (iter == gl_Cell_C(1, size(gl_Cell_C(1,:,1)) ,1)) then
                !    print*, gl_Cell_gran(:, iter)
                !    pause
                !end if

                if (i == 0) CYCLE
                k = gl_all_Gran(1, i)  ! Номер первого узла грани
                b = (/gl_x(k), gl_y(k), gl_z(k)/)
                a = gl_Gran_normal(:,i)
                di = dabs(DOT_PRODUCT(c,a) - DOT_PRODUCT(a,b))  ! Расстояние от точки до плоскости
                Vol = Vol + gl_Gran_square(i) * di/3.0
                dist = MIN(dist, di)
                ll = ll + 1
            end do



            !if (iter == gl_Cell_A(5, 1, 1)) then
            !    write(*, *) c
            !    print*, Vol
            !    print*, gl_Cell_gran(:, iter)
            !    pause
            !end if


            gl_Cell_Volume(iter) = Vol
            gl_Cell_dist(iter) = dist
            
            ! Альтернативное вычисление объёма пирамиды

            !        print*, Vol, ll
            !        pause
            
            ! do j = 1, 8
            !      i = gl_all_Cell(j, iter)
            !      pp(:, j) = (/gl_x(i), gl_y(i), gl_z(i)/)
            ! end do
            !
            !Vol = tetrahedron(pp(:, 1), pp(:, 2), pp(:, 4), pp(:, 5)) + tetrahedron(pp(:, 2), pp(:, 3), pp(:, 4), pp(:, 7)) + &
            !	tetrahedron(pp(:, 4), pp(:, 5), pp(:, 7), pp(:, 8)) + tetrahedron(pp(:, 2), pp(:, 4), pp(:, 5), pp(:, 7))
            !
            !if (gl_Cell_number(2, iter) /= 1) then
            !    gl_Cell_Volume(iter) = Vol
            !end if
            
            !print*, gl_Cell_type(i), gl_Cell_number(1, iter), gl_Cell_number(2, iter), gl_Cell_number(3, iter), gl_Cell_Volume(iter), Vol, &
            !	(gl_Cell_Volume(iter) - Vol)/Vol * 100
            !print*, c
            !print*, "_______________________________________"
            !PAUSE
            

        end do


        ! Вычислим, что за грань и ячейки       для разделения сетки на внутреннюю и внешнюю области, которые считаются отдельно
        
        ! Цикл для ячеек
        N1 = size(gl_Cell_A(:, 1, 1))
        N2 = size(gl_Cell_A(1, :, 1))
        N3 = size(gl_Cell_A(1, 1, :))
        
        gl_Gran_info = 0
        
        do k = 1, N3
            do j = 1, N2
                do i = 1, N1
                    if(i < par_n_IA) then
                        gl_Cell_info(gl_Cell_A(i, j, k)) = 0
                    else if (i < par_n_IB) then
                        gl_Cell_info(gl_Cell_A(i, j, k)) = 1
                    else
                        gl_Cell_info(gl_Cell_A(i, j, k)) = 2
                    end if
                end do
            end do
        end do
        
        ! Цикл для ячеек
        N1 = size(gl_Cell_B(:, 1, 1))
        N2 = size(gl_Cell_B(1, :, 1))
        N3 = size(gl_Cell_B(1, 1, :))
        
        do k = 1, N3
            do j = 1, N2
                do i = 1, N1
                    if(i < par_n_IA) then
                        gl_Cell_info(gl_Cell_B(i, j, k)) = 0
                    else if (i < par_n_IB) then
                        gl_Cell_info(gl_Cell_B(i, j, k)) = 1
                    else
                        gl_Cell_info(gl_Cell_B(i, j, k)) = 2
                    end if
                end do
            end do
        end do
        
        Ncell = size(gl_all_Cell(1, :))
        
        do  iter = 1, Ncell
            if(gl_Cell_info(iter) /= 0) CYCLE
            do j = 1, 6
                if(gl_Cell_gran(j, iter) > 0) then
                    gl_Gran_info(gl_Cell_gran(j,iter)) = 2
                end if
            end do
        end do
        
        do  iter = 1, Ncell
            if(gl_Cell_info(iter) /= 1) CYCLE
            do j = 1, 6
                if(gl_Cell_gran(j,iter) > 0) then
                    gl_Gran_info(gl_Cell_gran(j,iter)) = 1
                end if
            end do
        end do
        
        
        
        !Ngran = size(gl_all_Gran(1,:))
        !gl_Gran_info = 0
        !
        !do  iter = 1, Ngran
        !    i = gl_Gran_neighbour(1, iter)
        !    if (norm2(gl_Cell_center(:, i)) <= par_R0 * 45) then
        !        gl_Cell_info(i) = 0                                        ! Внутренняя сфера
        !        gl_Gran_info(iter) = gl_Gran_info(iter) + 1
        !        if (norm2(gl_Cell_center(:, i)) >= par_R0 * 40) then
        !            gl_Cell_info(i) = 1                                        ! Пограничный слой
        !        end if
        !    else
        !        gl_Cell_info(i) = 2
        !    end if
        !
        !    i = gl_Gran_neighbour(2, iter)
        !    if (i < 1) CYCLE
        !    if (norm2(gl_Cell_center(:, i)) <= par_R0 * 45) then
        !        gl_Cell_info(i) = 0                                        ! Внутренняя сфера
        !        gl_Gran_info(iter) = gl_Gran_info(iter) + 1
        !        if (norm2(gl_Cell_center(:, i)) >= par_R0 * 40) then
        !            gl_Cell_info(i) = 1                                        ! Пограничный слой
        !        end if
        !    else
        !        gl_Cell_info(i) = 2
        !    end if
        !
        !end do


    end subroutine calc_all_Gran

    subroutine Start_GD_1(steps)
        ! Считается просто газовая динамика (без мультифлюида, без магнитных полей)
        ! Без атомов и т.д. Без ТВД
        ! Функция паспараллелена на OPEN_MP
        use STORAGE
        use GEO_PARAM
        USE OMP_LIB
        use Solvers
        implicit none

        integer(4), intent(in) :: steps
        integer(4) :: st, gr, ngran, ncell, s1, s2, i, j, k
        real(8) :: qqq1(8), qqq2(8), qqq(8)  ! Переменные в ячейке
        real(8) :: dist, dsl, dsc, dsp
        real(8) :: POTOK(8)
        real(8) :: time, Volume, TT, U8, rad1, rad2, aa, bb, cc

        real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3

        ngran = size(gl_all_Gran(1, :))
        ncell = size(gl_all_Cell(1, :))

        call calc_all_Gran()   ! Посчитаем все объёмы, центры, площади и т.д.
        call Initial_conditions()  ! Задаём начальные условия

        time = 0.00000001

        !k = gl_Cell_C(1, size(gl_Cell_C(1,:,1)) ,1)
        !print*, gl_Cell_Volume(k)
        !
        !k = gl_Cell_C(1, size(gl_Cell_C(1,:,1)) - 1 ,1)
        !print*, gl_Cell_Volume(k)
        !
        !pause


        do st = 1, steps
            ! Делаем цикл по граням и считаем потоки через них
            TT = time
            time = 100000.0
            if (mod(st, 50) == 0)  print*, st, TT

            !$omp parallel

            !$omp do private(POTOK, s1, s2, qqq1, qqq2, dist, dsl, dsp, dsc, rad1, rad2, aa, bb, cc) &
            !$omp & reduction(min:time)
            do gr = 1, ngran
                POTOK = 0.0
                s1 = gl_Gran_neighbour(1, gr)
                if(s1 < 1) print*, "ERROR 1465gdfdyhe"
                s2 = gl_Gran_neighbour(2, gr)
                qqq1 = gl_Cell_par(:, s1)

                ! Попробуем снести плотность пропорционально квадрату
                if(norm2(qqq1(2:4))/sqrt(ggg*qqq1(5)/qqq1(1)) > 5.0) then
                    rad1 = norm2(gl_Cell_center(:, s1))
                    rad2 = norm2(gl_Gran_center(:, gr))
                    qqq1(1) = qqq1(1) * rad1**2 / rad2**2
                    qqq1(5) = qqq1(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
                    ! Скорости сносим в сферической С.К.
                    call spherical_skorost(gl_Cell_center(1, s1), gl_Cell_center(2, s1), gl_Cell_center(3, s1), &
                        qqq1(2), qqq1(3), qqq1(4), aa, bb, cc)
                    call dekard_skorost(gl_Gran_center(1, gr), gl_Gran_center(2, gr), gl_Gran_center(3, gr), &
                        aa, bb, cc, qqq1(2), qqq1(3), qqq1(4))

                end if

                if (s2 >= 1) then
                    if ( norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int .and. norm2(gl_Cell_center(:, s2)) <= par_R0 * par_R_int) CYCLE
                    qqq2 = gl_Cell_par(:, s2 )
                    dist = min(gl_Cell_dist(s1), gl_Cell_dist(s2))

                    ! Попробуем снести плотность пропорционально квадрату
                    if(norm2(qqq2(2:4))/sqrt(ggg*qqq2(5)/qqq2(1)) > 5.0) then
                        rad1 = norm2(gl_Cell_center(:, s2))
                        rad2 = norm2(gl_Gran_center(:, gr))
                        qqq2(1) = qqq2(1) * rad1**2 / rad2**2
                        qqq2(5) = qqq2(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
                        call spherical_skorost(gl_Cell_center(1, s2), gl_Cell_center(2, s2), gl_Cell_center(3, s2), &
                            qqq2(2), qqq2(3), qqq2(4), aa, bb, cc)
                        call dekard_skorost(gl_Gran_center(1, gr), gl_Gran_center(2, gr), gl_Gran_center(3, gr), &
                            aa, bb, cc, qqq2(2), qqq2(3), qqq2(4))
                    end if

                else  ! В случае граничных ячеек - граничные условия
                    if (norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int) CYCLE
                    if(s2 == -1) then  ! Набегающий поток
                        dist = gl_Cell_dist(s1)
                        qqq2 = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 1.0_8, 0.0_8, 0.0_8, 0.0_8/)
                    else  ! Здесь нужны мягкие условия
                        dist = gl_Cell_dist(s1)
                        qqq2 = qqq1
                        qqq2(5) = 1.0_8
                        if(qqq2(2) > par_Velosity_inf) then
                            qqq2(2) = par_Velosity_inf ! Отсос жидкости
                        end if
                    end if
                end if

                !print*, POTOK
                call chlld(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                    0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)
                !print*, POTOK
                !pause

                time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )   ! REDUCTION
                gl_Gran_POTOK(:, gr) = POTOK * gl_Gran_square(gr)
            end do
            !$omp end do

            !$omp barrier ! Синхронизация нитей
            ! Теперь цикл по ячейкам
            !$omp do private(POTOK, Volume, qqq, i, j, ro3, u3, v3, w3, p3)
            do gr = 1, ncell
                if (norm2(gl_Cell_center(:, gr)) <= par_R0 * par_R_int) CYCLE    ! Не считаем внутри сферы
                POTOK = 0.0
                Volume = gl_Cell_Volume(gr)
                qqq = gl_Cell_par(:, gr)
                ! Просуммируем потоки через грани
                do i = 1, 6
                    j = gl_Cell_gran(i, gr)
                    if (j == 0) CYCLE
                    if (gl_Gran_neighbour(1, j) == gr) then
                        POTOK = POTOK + gl_Gran_POTOK(1:8, j)
                    else
                        POTOK = POTOK - gl_Gran_POTOK(1:8, j)
                    end if
                end do

                ro3 = qqq(1) - TT * POTOK(1) / Volume
                if (ro3 <= 0.0_8) then
                    print*, "Ro < 0  1490 ", ro3
                    pause
                end if
                u3 = (qqq(1) * qqq(2) - TT * POTOK(2) / Volume) / ro3
                v3 = (qqq(1) * qqq(3) - TT * POTOK(3) / Volume) / ro3
                w3 = (qqq(1) * qqq(4) - TT * POTOK(4) / Volume) / ro3
                !bx3 = (bx * Volume_do / Volume - T[now1] * (K->Potok[4] + qqq(2) * K->Potok[8]) / Volume)
                !by3 = (by * Volume_do / Volume - T[now1] * (K->Potok[5] + qqq(3) * K->Potok[8]) / Volume)
                !bz3 = (bz * Volume_do / Volume - T[now1] * (K->Potok[6] + qqq(4) * K->Potok[8]) / Volume)
                !U8 = ( qqq(5) / (ggg - 1.0) + 0.5 * qqq(1) * norm2(qqq(2:4))**2 )
                p3 = ((  ( qqq(5) / (ggg - 1.0) + 0.5 * qqq(1) * norm2(qqq(2:4))**2 ) &
                    - TT * POTOK(5)/ Volume) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) ) * (ggg - 1.0)

                if (p3 <= 0.0_8) then
                    print*, "p < 0  1490 ", p3 , gl_Cell_center(:, gr)
                    p3 = 0.000001
                    !pause
                end if

                gl_Cell_par(:, gr) = (/ro3, u3, v3, w3, p3, 0.0_8, 0.0_8, 0.0_8/)

            end do
            !$omp end do

            !$omp end parallel

        end do

    end subroutine Start_GD_1

    subroutine Start_GD_2(steps)
        ! Считается просто газовая динамика + параметр Q (без мультифлюида, без магнитных полей)
        ! Без атомов и т.д. Без ТВД
        ! Функция паспараллелена на OPEN_MP
        use STORAGE
        use GEO_PARAM
        USE OMP_LIB
        use My_func
        use Solvers
        implicit none

        integer(4), intent(in) :: steps
        integer(4) :: st, gr, ngran, ncell, s1, s2, i, j, k
        real(8) :: qqq1(9), qqq2(9), qqq(9)  ! Переменные в ячейке
        real(8) :: dist, dsl, dsc, dsp
        real(8) :: POTOK(9)
        real(8) :: time, Volume, TT, U8, rad1, rad2, aa, bb, cc

        real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3

        ngran = size(gl_all_Gran(1, :))
        ncell = size(gl_all_Cell(1, :))

        call calc_all_Gran()   ! Посчитаем все объёмы, центры, площади и т.д.
        call Initial_conditions()  ! Задаём начальные условия

        time = 0.00000001

        !k = gl_Cell_C(1, size(gl_Cell_C(1,:,1)) ,1)
        !print*, gl_Cell_Volume(k)
        !
        !k = gl_Cell_C(1, size(gl_Cell_C(1,:,1)) - 1 ,1)
        !print*, gl_Cell_Volume(k)
        !
        !pause


        do st = 1, steps
            ! Делаем цикл по граням и считаем потоки через них
            TT = time
            time = 100000.0
            if (mod(st, 100) == 0)  print*, st, TT

            !$omp parallel

            !$omp do private(POTOK, s1, s2, qqq1, qqq2, dist, dsl, dsp, dsc, rad1, rad2, aa, bb, cc) &
            !$omp & reduction(min:time)
            do gr = 1, ngran
                POTOK = 0.0
                s1 = gl_Gran_neighbour(1, gr)
                if(s1 < 1) print*, "ERROR 1465gdfdyhe"
                s2 = gl_Gran_neighbour(2, gr)
                qqq1 = gl_Cell_par(:, s1)

                ! Попробуем снести плотность пропорционально квадрату
                if(norm2(qqq1(2:4))/sqrt(ggg*qqq1(5)/qqq1(1)) > 5.0) then
                    rad1 = norm2(gl_Cell_center(:, s1))
                    rad2 = norm2(gl_Gran_center(:, gr))
                    qqq1(1) = qqq1(1) * rad1**2 / rad2**2
                    qqq1(9) = qqq1(9) * rad1**2 / rad2**2
                    qqq1(5) = qqq1(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
                    ! Скорости сносим в сферической С.К.
                    call spherical_skorost(gl_Cell_center(1, s1), gl_Cell_center(2, s1), gl_Cell_center(3, s1), &
                        qqq1(2), qqq1(3), qqq1(4), aa, bb, cc)
                    call dekard_skorost(gl_Gran_center(1, gr), gl_Gran_center(2, gr), gl_Gran_center(3, gr), &
                        aa, bb, cc, qqq1(2), qqq1(3), qqq1(4))

                end if

                if (s2 >= 1) then
                    if ( norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int .and. norm2(gl_Cell_center(:, s2)) <= par_R0 * par_R_int) CYCLE
                    qqq2 = gl_Cell_par(:, s2 )
                    dist = min(gl_Cell_dist(s1), gl_Cell_dist(s2))

                    ! Попробуем снести плотность пропорционально квадрату
                    if(norm2(qqq2(2:4))/sqrt(ggg*qqq2(5)/qqq2(1)) > 5.0) then
                        rad1 = norm2(gl_Cell_center(:, s2))
                        rad2 = norm2(gl_Gran_center(:, gr))
                        qqq2(1) = qqq2(1) * rad1**2 / rad2**2
                        qqq2(9) = qqq2(9) * rad1**2 / rad2**2
                        qqq2(5) = qqq2(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
                        call spherical_skorost(gl_Cell_center(1, s2), gl_Cell_center(2, s2), gl_Cell_center(3, s2), &
                            qqq2(2), qqq2(3), qqq2(4), aa, bb, cc)
                        call dekard_skorost(gl_Gran_center(1, gr), gl_Gran_center(2, gr), gl_Gran_center(3, gr), &
                            aa, bb, cc, qqq2(2), qqq2(3), qqq2(4))
                    end if

                else  ! В случае граничных ячеек - граничные условия
                    if (norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int) CYCLE
                    if(s2 == -1) then  ! Набегающий поток
                        dist = gl_Cell_dist(s1)
                        qqq2 = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 1.0_8, 0.0_8, 0.0_8, 0.0_8, 100.0_8/)
                    else  ! Здесь нужны мягкие условия
                        dist = gl_Cell_dist(s1)
                        qqq2 = qqq1
                        qqq2(5) = 1.0_8
                        if(qqq2(2) > par_Velosity_inf) then
                            qqq2(2) = par_Velosity_inf ! Отсос жидкости
                        end if
                    end if
                end if

                !print*, POTOK
                !print*, qqq1, qqq2
                call chlld_Q(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                    0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK, .False.)
                !print*, POTOK
                !pause
                !
                !call chlld(2, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                !               0.0_8, qqq1(1:8), qqq2(1:8), dsl, dsp, dsc, POTOK(1:8))
                !print*, POTOK
                !pause

                time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )   ! REDUCTION
                gl_Gran_POTOK(1:9, gr) = POTOK * gl_Gran_square(gr)
            end do
            !$omp end do

            !$omp barrier ! Синхронизация нитей
            ! Теперь цикл по ячейкам
            !$omp do private(POTOK, Volume, qqq, i, j, ro3, u3, v3, w3, p3, Q3)
            do gr = 1, ncell
                if (norm2(gl_Cell_center(:, gr)) <= par_R0 * par_R_int) CYCLE    ! Не считаем внутри сферы
                POTOK = 0.0
                Volume = gl_Cell_Volume(gr)
                qqq = gl_Cell_par(:, gr)
                ! Просуммируем потоки через грани
                do i = 1, 6
                    j = gl_Cell_gran(i, gr)
                    if (j == 0) CYCLE
                    if (gl_Gran_neighbour(1, j) == gr) then
                        POTOK = POTOK + gl_Gran_POTOK(:, j)
                    else
                        POTOK = POTOK - gl_Gran_POTOK(:, j)
                    end if
                end do

                ro3 = qqq(1) - time * POTOK(1) / Volume
                Q3 = qqq(9) - time * POTOK(9) / Volume
                if (ro3 <= 0.0_8) then
                    print*, "Ro < 0  1490 ", ro3
                    pause
                end if
                u3 = (qqq(1) * qqq(2) - time * POTOK(2) / Volume) / ro3
                v3 = (qqq(1) * qqq(3) - time * POTOK(3) / Volume) / ro3
                w3 = (qqq(1) * qqq(4) - time * POTOK(4) / Volume) / ro3
                !bx3 = (bx * Volume_do / Volume - T[now1] * (K->Potok[4] + qqq(2) * K->Potok[8]) / Volume)
                !by3 = (by * Volume_do / Volume - T[now1] * (K->Potok[5] + qqq(3) * K->Potok[8]) / Volume)
                !bz3 = (bz * Volume_do / Volume - T[now1] * (K->Potok[6] + qqq(4) * K->Potok[8]) / Volume)
                !U8 = ( qqq(5) / (ggg - 1.0) + 0.5 * qqq(1) * norm2(qqq(2:4))**2 )
                p3 = ((  ( qqq(5) / (ggg - 1.0) + 0.5 * qqq(1) * norm2(qqq(2:4))**2 ) &
                    - time * POTOK(5)/ Volume) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) ) * (ggg - 1.0)

                if (p3 <= 0.0_8) then
                    print*, "p < 0  1490 ", p3 , gl_Cell_center(:, gr)
                    p3 = 0.000001
                    !pause
                end if

                gl_Cell_par(:, gr) = (/ro3, u3, v3, w3, p3, 0.0_8, 0.0_8, 0.0_8, Q3/)

            end do
            !$omp end do

            !$omp end parallel

        end do

    end subroutine Start_GD_2

    subroutine Start_GD_3(steps)
        ! Считается газовая динамика + параметр Q + мультифлюид (без магнитных полей)
        ! Без атомов и т.д. Без ТВД
        ! Функция паспараллелена на OPEN_MP
        use STORAGE
        use GEO_PARAM
        USE OMP_LIB
        use My_func
        use Solvers
        implicit none

        integer(4), intent(in) :: steps
        integer(4) :: st, gr, ngran, ncell, s1, s2, i, j, k, zone
        real(8) :: qqq1(9), qqq2(9), qqq(9)  ! Переменные в ячейке
        real(8) :: fluid1(5, par_n_sort), fluid2(5, par_n_sort)
        real(8) :: dist, dsl, dsc, dsp
        real(8) :: POTOK(9)
        real(8) :: POTOK_MF(5)
        real(8) :: POTOK_MF_all(5, 4)
        real(8) :: time, Volume, TT, U8, rad1, rad2, aa, bb, cc
        real(8) :: SOURSE(5,par_n_sort + 1)  ! Источники массы, импульса и энергии для плазмы и каждого сорта мультифлюида

        real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3

        logical :: l_1

        ngran = size(gl_all_Gran(1, :))
        ncell = size(gl_all_Cell(1, :))

        !call calc_all_Gran()   ! Посчитаем все объёмы, центры, площади и т.д.
        !call Initial_conditions()  ! Задаём начальные условия

        time = 0.000000001

        do st = 1, steps
            ! Делаем цикл по граням и считаем потоки через них
            TT = time
            time = 100000.0

            !$omp parallel

            !$omp do private(POTOK, s1, s2, qqq1, qqq2, dist, dsl, dsp, dsc, rad1, rad2, aa, bb, cc, fluid1, fluid2, POTOK_MF) &
            !$omp & reduction(min:time)
            do gr = 1, ngran
                if(gl_Gran_info(gr) == 2) CYCLE
                POTOK = 0.0
                s1 = gl_Gran_neighbour(1, gr)
                s2 = gl_Gran_neighbour(2, gr)
                qqq1 = gl_Cell_par(:, s1)
                fluid1 = gl_Cell_par_MF(:, :, s1)   ! Загрузили параметры жидкостей для мультифлюида

                ! Попробуем снести плотность пропорционально квадрату
                if(norm2(qqq1(2:4))/sqrt(ggg*qqq1(5)/qqq1(1)) > 5.0) then
                    rad1 = norm2(gl_Cell_center(:, s1))
                    rad2 = norm2(gl_Gran_center(:, gr))
                    qqq1(1) = qqq1(1) * rad1**2 / rad2**2
                    qqq1(9) = qqq1(9) * rad1**2 / rad2**2
                    qqq1(5) = qqq1(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
                    ! Скорости сносим в сферической С.К.
                    call spherical_skorost(gl_Cell_center(1, s1), gl_Cell_center(2, s1), gl_Cell_center(3, s1), &
                        qqq1(2), qqq1(3), qqq1(4), aa, bb, cc)
                    call dekard_skorost(gl_Gran_center(1, gr), gl_Gran_center(2, gr), gl_Gran_center(3, gr), &
                        aa, bb, cc, qqq1(2), qqq1(3), qqq1(4))

                    call spherical_skorost(gl_Cell_center(1, s1), gl_Cell_center(2, s1), gl_Cell_center(3, s1), &
                        fluid1(2, 1), fluid1(3, 1), fluid1(4, 1), aa, bb, cc)
                    call dekard_skorost(gl_Gran_center(1, gr), gl_Gran_center(2, gr), gl_Gran_center(3, gr), &
                        aa, bb, cc, fluid1(2, 1), fluid1(3, 1), fluid1(4, 1))
                end if

                if (s2 >= 1) then
                    !if ( norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int .and. norm2(gl_Cell_center(:, s2)) <= par_R0 * par_R_int) CYCLE
                    qqq2 = gl_Cell_par(:, s2)
                    fluid2 = gl_Cell_par_MF(:, :, s2)   ! Загрузили параметры жидкостей для мультифлюида
                    dist = min(gl_Cell_dist(s1), gl_Cell_dist(s2))

                    ! Попробуем снести плотность пропорционально квадрату
                    if(norm2(qqq2(2:4))/sqrt(ggg*qqq2(5)/qqq2(1)) > 5.0) then
                        rad1 = norm2(gl_Cell_center(:, s2))
                        rad2 = norm2(gl_Gran_center(:, gr))
                        qqq2(1) = qqq2(1) * rad1**2 / rad2**2
                        qqq2(9) = qqq2(9) * rad1**2 / rad2**2
                        qqq2(5) = qqq2(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
                        call spherical_skorost(gl_Cell_center(1, s2), gl_Cell_center(2, s2), gl_Cell_center(3, s2), &
                            qqq2(2), qqq2(3), qqq2(4), aa, bb, cc)
                        call dekard_skorost(gl_Gran_center(1, gr), gl_Gran_center(2, gr), gl_Gran_center(3, gr), &
                            aa, bb, cc, qqq2(2), qqq2(3), qqq2(4))

                        call spherical_skorost(gl_Cell_center(1, s2), gl_Cell_center(2, s2), gl_Cell_center(3, s2), &
                            fluid2(2, 1), fluid2(3, 1), fluid2(4, 1), aa, bb, cc)
                        call dekard_skorost(gl_Gran_center(1, gr), gl_Gran_center(2, gr), gl_Gran_center(3, gr), &
                            aa, bb, cc, fluid2(2, 1), fluid2(3, 1), fluid2(4, 1))
                    end if

                else  ! В случае граничных ячеек - граничные условия
                    !if (norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int) CYCLE
                    if(s2 == -1) then  ! Набегающий поток
                        dist = gl_Cell_dist(s1)
                        qqq2 = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 1.0_8, 0.0_8, 0.0_8, 0.0_8, 100.0_8/)
                        fluid2(:, 1) = fluid1(:, 1)
                        fluid2(:, 2) = fluid1(:, 2)
                        fluid2(:, 3) = fluid1(:, 3)
                        fluid2(:, 4) = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 0.5_8/)
                    else  ! Здесь нужны мягкие условия
                        dist = gl_Cell_dist(s1)
                        qqq2 = qqq1
                        fluid2 = fluid1
                        qqq2(5) = 1.0_8
                        if(qqq2(2) > par_Velosity_inf) then
                            qqq2(2) = par_Velosity_inf ! Отсос жидкости
                        end if

                        !if(fluid2(2,1) > par_Velosity_inf) then
                        !    fluid2(2,1) = par_Velosity_inf ! Отсос жидкости
                        !end if
                        !if(fluid2(2,2) > par_Velosity_inf) then
                        !    fluid2(2,2) = par_Velosity_inf ! Отсос жидкости
                        !end if
                        !if(fluid2(2,3) > par_Velosity_inf) then
                        !    fluid2(2,3) = par_Velosity_inf ! Отсос жидкости
                        !end if
                        !if(fluid2(2,4) > par_Velosity_inf) then
                        !    fluid2(2,4) = par_Velosity_inf ! Отсос жидкости
                        !end if
                    end if
                end if


                call chlld_Q(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                    0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)
                time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )   ! REDUCTION
                gl_Gran_POTOK(1:9, gr) = POTOK * gl_Gran_square(gr)

                !print*, time
                !pause

                call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                    0.0_8, fluid1(:, 1), fluid2(:, 1), dsl, dsp, dsc, POTOK_MF)
                time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )
                gl_Gran_POTOK_MF(:, 1, gr) = POTOK_MF * gl_Gran_square(gr)

                !print*, time
                !print*, fluid1(:, 1)
                !print*, fluid2(:, 1)
                !print*, dsl, dsp, dsc
                !pause

                call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                    0.0_8, fluid1(:, 2), fluid2(:, 2), dsl, dsp, dsc, POTOK_MF)
                time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )
                gl_Gran_POTOK_MF(:, 2, gr) = POTOK_MF * gl_Gran_square(gr)

                !print*, time
                !pause

                call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                    0.0_8, fluid1(:, 3), fluid2(:, 3), dsl, dsp, dsc, POTOK_MF)
                time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )
                gl_Gran_POTOK_MF(:, 3, gr) = POTOK_MF * gl_Gran_square(gr)

                !print*, time
                !pause

                call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                    0.0_8, fluid1(:, 4), fluid2(:, 4), dsl, dsp, dsc, POTOK_MF)
                time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )
                gl_Gran_POTOK_MF(:, 4, gr) = POTOK_MF * gl_Gran_square(gr)

                !if (dabs(POTOK_MF(1)) > 10.0) then
                !    print*, "POTOK ", POTOK_MF
                !    print*, "fluid1 ",fluid1(:, 4)
                !    print*, "fluid2 ",fluid2(:, 4)
                !    pause
                !end if

            end do
            !$omp end do



            !$omp barrier ! Синхронизация нитей

            !!$omp single
            !    if (mod(st, 5) == 0)  print*, st, time
            !!$omp end single

            ! Теперь цикл по ячейкам
            !$omp do private(POTOK, Volume, qqq, i, j, ro3, u3, v3, w3, p3, Q3, POTOK_MF_all, zone, l_1, fluid1, SOURSE)
            do gr = 1, ncell
                if(gl_Cell_info(gr) == 0) CYCLE
                l_1 = .TRUE.
                if (norm2(gl_Cell_center(:, gr)) <= 1.0 * par_R0) l_1 = .FALSE.    ! Не считаем внутри сферы
                POTOK = 0.0
                SOURSE = 0.0
                POTOK_MF_all = 0.0
                Volume = gl_Cell_Volume(gr)
                qqq = gl_Cell_par(:, gr)
                fluid1 = gl_Cell_par_MF(:, :, gr)
                ! Просуммируем потоки через грани
                do i = 1, 6
                    j = gl_Cell_gran(i, gr)
                    if (j == 0) CYCLE
                    if (gl_Gran_neighbour(1, j) == gr) then
                        POTOK = POTOK + gl_Gran_POTOK(1:9, j)
                        POTOK_MF_all(:, 1) = POTOK_MF_all(:, 1) +  gl_Gran_POTOK_MF(:, 1, j)
                        POTOK_MF_all(:, 2) = POTOK_MF_all(:, 2) +  gl_Gran_POTOK_MF(:, 2, j)
                        POTOK_MF_all(:, 3) = POTOK_MF_all(:, 3) +  gl_Gran_POTOK_MF(:, 3, j)
                        POTOK_MF_all(:, 4) = POTOK_MF_all(:, 4) +  gl_Gran_POTOK_MF(:, 4, j)
                    else
                        POTOK = POTOK - gl_Gran_POTOK(1:9, j)
                        POTOK_MF_all(:, 1) = POTOK_MF_all(:, 1) -  gl_Gran_POTOK_MF(:, 1, j)
                        POTOK_MF_all(:, 2) = POTOK_MF_all(:, 2) -  gl_Gran_POTOK_MF(:, 2, j)
                        POTOK_MF_all(:, 3) = POTOK_MF_all(:, 3) -  gl_Gran_POTOK_MF(:, 3, j)
                        POTOK_MF_all(:, 4) = POTOK_MF_all(:, 4) -  gl_Gran_POTOK_MF(:, 4, j)
                    end if
                end do


                ! Определяем зону в которой находимся
                if(qqq(9)/qqq(1) < 50.0) then

                    ! Перенормируем параметры
                    !qqq(2:4) = qqq(2:4) * (par_chi/par_chi_real)
                    !qqq(1) = qqq(1) / (par_chi/par_chi_real)**2

                    if(norm2(qqq(2:4))/sqrt(ggg*qqq(5)/qqq(1)) > 3.0) then
                        zone = 1
                    else
                        zone = 2
                    end if
                else
                    if(norm2(qqq(2:4))/sqrt(ggg*qqq(5)/qqq(1)) > 1.1) then
                        zone = 4
                    else
                        zone = 3
                    end if
                end if

                ! Перенормируем первую жидкость
                !fluid1(2:4, 1) = fluid1(2:4, 1) * (par_chi/par_chi_real)
                !fluid1(1, 1) = fluid1(1, 1) / (par_chi/par_chi_real)**2

                !fluid1(2:4, 2) = fluid1(2:4, 2) * (3.0_8)
                !fluid1(1, 2) = fluid1(1, 2) / (3.0_8)**2

                call Calc_sourse_MF(qqq, fluid1, SOURSE, zone)  ! Вычисляем источники

                ! Перенормируем первую жидкость обратно
                !fluid1(2:4, 1) = fluid1(2:4, 1) / (par_chi/par_chi_real)
                !fluid1(1, 1) = fluid1(1, 1) * (par_chi/par_chi_real)**2
                !SOURSE(5, 2) = SOURSE(5, 2)/ (par_chi/par_chi_real)
                !SOURSE(1, 2) = SOURSE(1, 2)* (par_chi/par_chi_real)

                !fluid1(2:4, 2) = fluid1(2:4, 2) / (3.0_8)
                !fluid1(1, 2) = fluid1(1, 2) * (3.0_8)**2
                !SOURSE(5, 3) = SOURSE(5, 3)/ (3.0_8)
                !SOURSE(1, 3) = SOURSE(1, 3)* (3.0_8)

                ! Перенормируем обратно
                !if(zone == 1 .or. zone == 2) then
                !    qqq(2:4) = qqq(2:4) / (par_chi/par_chi_real)
                !    qqq(1) = qqq(1) * (par_chi/par_chi_real)**2
                !    SOURSE(5, 1) = SOURSE(5, 1)/ (par_chi/par_chi_real)
                !end if


                if (l_1 == .TRUE.) then
                    ro3 = qqq(1) - time * POTOK(1) / Volume
                    Q3 = qqq(9) - time * POTOK(9) / Volume
                    if (ro3 <= 0.0_8) then
                        print*, "Ro < 0  1490 ", ro3, gl_Cell_center(:, gr)
                        pause
                    end if
                    u3 = (qqq(1) * qqq(2) - time * POTOK(2) / Volume + time * SOURSE(2, 1)) / ro3
                    v3 = (qqq(1) * qqq(3) - time * POTOK(3) / Volume + time * SOURSE(3, 1)) / ro3
                    w3 = (qqq(1) * qqq(4) - time * POTOK(4) / Volume + time * SOURSE(4, 1)) / ro3
                    p3 = ((  ( qqq(5) / (ggg - 1.0) + 0.5 * qqq(1) * norm2(qqq(2:4))**2 ) &
                        - time * POTOK(5)/ Volume + time * SOURSE(5, 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) ) * (ggg - 1.0)

                    if (p3 <= 0.0_8) then
                        !print*, "p < 0  plasma 2028 ", p3 , gl_Cell_center(:, gr)
                        p3 = 0.000001
                        !pause
                    end if

                    gl_Cell_par(:, gr) = (/ro3, u3, v3, w3, p3, 0.0_8, 0.0_8, 0.0_8, Q3/)

                end if

                ! Теперь посчитаем законы сохранения для остальных жидкостей

                do i = 1, 4
                    if (i == 1 .and. l_1 == .FALSE.) CYCLE       ! Пропускаем внутреннюю сферу для сорта 1
                    if (l_1 == .FALSE.) SOURSE(:, i + 1) = 0.0       ! Пропускаем внутреннюю сферу для сорта 1
                    ro3 = fluid1(1, i) - time * POTOK_MF_all(1, i) / Volume + time * SOURSE(1, i + 1)
                    if (ro3 <= 0.0_8) then
                        print*, "Ro < 0  in ", i,  ro3
                        print*, "gl_Cell_center(:, gr)", gl_Cell_center(:, gr)
                        print*, "gl_Cell_info(gr)", gl_Cell_info(gr)
                        print*, fluid1(:, i)
                        print*, l_1
                        print*, SOURSE(:, i + 1)
                        print*, time
                        print*, "POTOK ", POTOK_MF_all(:, i)
                        pause
                    end if
                    u3 = (fluid1(1, i) * fluid1(2, i) - time * POTOK_MF_all(2, i) / Volume + time * SOURSE(2, i + 1)) / ro3
                    v3 = (fluid1(1, i) * fluid1(3, i) - time * POTOK_MF_all(3, i) / Volume + time * SOURSE(3, i + 1)) / ro3
                    w3 = (fluid1(1, i) * fluid1(4, i) - time * POTOK_MF_all(4, i) / Volume + time * SOURSE(4, i + 1)) / ro3
                    p3 = ((  ( fluid1(5, i) / (ggg - 1.0) + 0.5 * fluid1(1, i) * norm2(fluid1(2:4, i))**2 ) &
                        - time * POTOK_MF_all(5, i)/ Volume + time * SOURSE(5, i + 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) ) * (ggg - 1.0)
                    if (p3 <= 0.0_8) then
                        !print*, "p < 0  in ", i, p3, gl_Cell_center(:, gr)
                        !print*, fluid1(:, i)
                        !print*, l_1
                        !print*, "Sourse", SOURSE(:, i + 1)
                        !print*, time
                        !print*, "zone", zone, norm2(qqq(2:4))/sqrt(ggg*qqq(5)/qqq(1))
                        !print*, "POTOK ", POTOK_MF_all(:, i)
                        !print*, "qqq", qqq
                        !pause
                        p3 = 0.000001
                        !pause
                    end if

                    gl_Cell_par_MF(:, i, gr) = (/ro3, u3, v3, w3, p3/)
                end do

            end do
            !$omp end do

            !$omp end parallel

        end do

    end subroutine Start_GD_3

    subroutine Start_GD_3_inner(steps)
        ! Считается газовая динамика + параметр Q + мультифлюид (без магнитных полей)
        ! Без атомов и т.д. Без ТВД
        ! Функция паспараллелена на OPEN_MP
        use STORAGE
        use GEO_PARAM
        USE OMP_LIB
        use My_func
        use Solvers
        implicit none

        integer(4), intent(in) :: steps
        integer(4) :: st, gr, ngran, ncell, s1, s2, i, j, k, zone, iter
        real(8) :: qqq1(9), qqq2(9), qqq(9)  ! Переменные в ячейке
        real(8) :: fluid1(5, par_n_sort), fluid2(5, par_n_sort)
        real(8) :: dist, dsl, dsc, dsp
        real(8) :: POTOK(9)
        real(8) :: POTOK_MF(5)
        real(8) :: POTOK_MF_all(5, 4)
        real(8) :: time, Volume, TT, U8, rad1, rad2, aa, bb, cc
        real(8) :: SOURSE(5,par_n_sort + 1)  ! Источники массы, импульса и энергии для плазмы и каждого сорта мультифлюида

        real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3

        logical :: l_1

        ngran = size(gl_all_Gran_inner(:))
        ncell = size(gl_all_Cell_inner(:))

        !call calc_all_Gran()   ! Посчитаем все объёмы, центры, площади и т.д.
        !call Initial_conditions()  ! Задаём начальные условия

        time = 0.000000001

        do st = 1, steps
            ! Делаем цикл по граням и считаем потоки через них
            TT = time
            time = 100000.0
            !if (mod(st, 10) == 0)  print*, "inner",st, TT

            !$omp parallel

            !$omp do private(POTOK, s1, s2, qqq1, qqq2, dist, dsl, dsp, dsc, rad1, rad2, aa, bb, cc, fluid1, fluid2, POTOK_MF, gr) &
            !$omp & reduction(min:time)
            do iter = 1, ngran
                gr = gl_all_Gran_inner(iter)
                !if(gl_Gran_info(gr) == 0) CYCLE
                POTOK = 0.0
                s1 = gl_Gran_neighbour(1, gr)
                s2 = gl_Gran_neighbour(2, gr)
                qqq1 = gl_Cell_par(:, s1)
                fluid1 = gl_Cell_par_MF(:, :, s1)   ! Загрузили параметры жидкостей для мультифлюида

                ! Попробуем снести плотность пропорционально квадрату
                if(norm2(qqq1(2:4))/sqrt(ggg*qqq1(5)/qqq1(1)) > 2.5) then
                    rad1 = norm2(gl_Cell_center(:, s1))
                    rad2 = norm2(gl_Gran_center(:, gr))
                    qqq1(1) = qqq1(1) * rad1**2 / rad2**2
                    qqq1(9) = qqq1(9) * rad1**2 / rad2**2
                    qqq1(5) = qqq1(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
                    ! Скорости сносим в сферической С.К.
                    call spherical_skorost(gl_Cell_center(1, s1), gl_Cell_center(2, s1), gl_Cell_center(3, s1), &
                        qqq1(2), qqq1(3), qqq1(4), aa, bb, cc)
                    call dekard_skorost(gl_Gran_center(1, gr), gl_Gran_center(2, gr), gl_Gran_center(3, gr), &
                        aa, bb, cc, qqq1(2), qqq1(3), qqq1(4))

                    call spherical_skorost(gl_Cell_center(1, s1), gl_Cell_center(2, s1), gl_Cell_center(3, s1), &
                        fluid1(2, 1), fluid1(3, 1), fluid1(4, 1), aa, bb, cc)
                    call dekard_skorost(gl_Gran_center(1, gr), gl_Gran_center(2, gr), gl_Gran_center(3, gr), &
                        aa, bb, cc, fluid1(2, 1), fluid1(3, 1), fluid1(4, 1))
                end if

                if (s2 >= 1) then
                    !if ( norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int .and. norm2(gl_Cell_center(:, s2)) <= par_R0 * par_R_int) CYCLE
                    qqq2 = gl_Cell_par(:, s2)
                    fluid2 = gl_Cell_par_MF(:, :, s2)   ! Загрузили параметры жидкостей для мультифлюида
                    dist = min(gl_Cell_dist(s1), gl_Cell_dist(s2))

                    ! Попробуем снести плотность пропорционально квадрату
                    if(norm2(qqq2(2:4))/sqrt(ggg*qqq2(5)/qqq2(1)) > 2.5) then
                        rad1 = norm2(gl_Cell_center(:, s2))
                        rad2 = norm2(gl_Gran_center(:, gr))
                        qqq2(1) = qqq2(1) * rad1**2 / rad2**2
                        qqq2(9) = qqq2(9) * rad1**2 / rad2**2
                        qqq2(5) = qqq2(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
                        call spherical_skorost(gl_Cell_center(1, s2), gl_Cell_center(2, s2), gl_Cell_center(3, s2), &
                            qqq2(2), qqq2(3), qqq2(4), aa, bb, cc)
                        call dekard_skorost(gl_Gran_center(1, gr), gl_Gran_center(2, gr), gl_Gran_center(3, gr), &
                            aa, bb, cc, qqq2(2), qqq2(3), qqq2(4))

                        call spherical_skorost(gl_Cell_center(1, s2), gl_Cell_center(2, s2), gl_Cell_center(3, s2), &
                            fluid2(2, 1), fluid2(3, 1), fluid2(4, 1), aa, bb, cc)
                        call dekard_skorost(gl_Gran_center(1, gr), gl_Gran_center(2, gr), gl_Gran_center(3, gr), &
                            aa, bb, cc, fluid2(2, 1), fluid2(3, 1), fluid2(4, 1))
                    end if

                else  ! В случае граничных ячеек - граничные условия
                    !if (norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int) CYCLE
                    if(s2 == -1) then  ! Набегающий поток
                        dist = gl_Cell_dist(s1)
                        qqq2 = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 1.0_8, 0.0_8, 0.0_8, 0.0_8, 100.0_8/)
                        fluid2(:, 1) = (/0.000001_8, 0.0_8, 0.0_8, 0.0_8, 0.000001_8/)
                        fluid2(:, 2) = (/0.000001_8, 0.0_8, 0.0_8, 0.0_8, 0.000001_8/)
                        fluid2(:, 3) = (/0.000001_8, 0.0_8, 0.0_8, 0.0_8, 0.000001_8/)
                        fluid2(:, 4) = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 0.5_8/)
                    else  ! Здесь нужны мягкие условия (это задняя стенка)
                        dist = gl_Cell_dist(s1)
                        qqq2 = qqq1
                        fluid2 = fluid1
                        qqq2(5) = 1.0_8
                        if(qqq2(2) > par_Velosity_inf) then
                            qqq2(2) = par_Velosity_inf ! Отсос жидкости
                        end if

                        !if(fluid2(2,1) > par_Velosity_inf) then
                        !    fluid2(2,1) = par_Velosity_inf ! Отсос жидкости
                        !end if
                        !if(fluid2(2,2) > par_Velosity_inf) then
                        !    fluid2(2,2) = par_Velosity_inf ! Отсос жидкости
                        !end if
                        !if(fluid2(2,3) > par_Velosity_inf) then
                        !    fluid2(2,3) = par_Velosity_inf ! Отсос жидкости
                        !end if
                        !if(fluid2(2,4) > par_Velosity_inf) then
                        !    fluid2(2,4) = par_Velosity_inf ! Отсос жидкости
                        !end if
                    end if
                end if


                call chlld_Q(2, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                    0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)
                time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )   ! REDUCTION
                gl_Gran_POTOK(1:9, gr) = POTOK * gl_Gran_square(gr)

                !print*, time
                !pause

                call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                    0.0_8, fluid1(:, 1), fluid2(:, 1), dsl, dsp, dsc, POTOK_MF)
                time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )
                gl_Gran_POTOK_MF(:, 1, gr) = POTOK_MF * gl_Gran_square(gr)

                !print*, time
                !print*, fluid1(:, 1)
                !print*, fluid2(:, 1)
                !print*, dsl, dsp, dsc
                !pause

                call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                    0.0_8, fluid1(:, 2), fluid2(:, 2), dsl, dsp, dsc, POTOK_MF)
                time = min(time, 0.8 * dist/max(dabs(dsl), dabs(dsp)) )
                gl_Gran_POTOK_MF(:, 2, gr) = POTOK_MF * gl_Gran_square(gr)

                !print*, time
                !pause

                call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                    0.0_8, fluid1(:, 3), fluid2(:, 3), dsl, dsp, dsc, POTOK_MF)
                time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )
                gl_Gran_POTOK_MF(:, 3, gr) = POTOK_MF * gl_Gran_square(gr)

                !print*, time
                !pause

                call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                    0.0_8, fluid1(:, 4), fluid2(:, 4), dsl, dsp, dsc, POTOK_MF)
                time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )
                gl_Gran_POTOK_MF(:, 4, gr) = POTOK_MF * gl_Gran_square(gr)

                !if (dabs(POTOK_MF(1)) > 10.0) then
                !    print*, "POTOK ", POTOK_MF
                !    print*, "fluid1 ",fluid1(:, 4)
                !    print*, "fluid2 ",fluid2(:, 4)
                !    pause
                !end if

            end do
            !$omp end do

            !$omp barrier ! Синхронизация нитей

            !!$omp single
            !    if (mod(st, 50) == 0)  print*, st, time
            !!$omp end single

            ! Теперь цикл по ячейкам
            !$omp do private(POTOK, Volume, qqq, i, j, ro3, u3, v3, w3, p3, Q3, POTOK_MF_all, zone, l_1, fluid1, SOURSE, gr)
            do iter = 1, ncell
                gr = gl_all_Cell_inner(iter)
                !if(gl_Cell_info(gr) == 2) CYCLE
                l_1 = .TRUE.
                if ((gl_Cell_type(gr) == "A" .or. gl_Cell_type(gr) == "B").and.(gl_Cell_number(1, gr) <= 2) ) l_1 = .FALSE.    ! Не считаем в первых двух ячейках
                POTOK = 0.0
                SOURSE = 0.0
                POTOK_MF_all = 0.0
                Volume = gl_Cell_Volume(gr)
                qqq = gl_Cell_par(:, gr)
                fluid1 = gl_Cell_par_MF(:, :, gr)
                ! Просуммируем потоки через грани
                do i = 1, 6
                    j = gl_Cell_gran(i, gr)
                    if (j == 0) CYCLE
                    if (gl_Gran_neighbour(1, j) == gr) then
                        POTOK = POTOK + gl_Gran_POTOK(:, j)
                        POTOK_MF_all(:, 1) = POTOK_MF_all(:, 1) +  gl_Gran_POTOK_MF(:, 1, j)
                        POTOK_MF_all(:, 2) = POTOK_MF_all(:, 2) +  gl_Gran_POTOK_MF(:, 2, j)
                        POTOK_MF_all(:, 3) = POTOK_MF_all(:, 3) +  gl_Gran_POTOK_MF(:, 3, j)
                        POTOK_MF_all(:, 4) = POTOK_MF_all(:, 4) +  gl_Gran_POTOK_MF(:, 4, j)
                    else
                        POTOK = POTOK - gl_Gran_POTOK(:, j)
                        POTOK_MF_all(:, 1) = POTOK_MF_all(:, 1) -  gl_Gran_POTOK_MF(:, 1, j)
                        POTOK_MF_all(:, 2) = POTOK_MF_all(:, 2) -  gl_Gran_POTOK_MF(:, 2, j)
                        POTOK_MF_all(:, 3) = POTOK_MF_all(:, 3) -  gl_Gran_POTOK_MF(:, 3, j)
                        POTOK_MF_all(:, 4) = POTOK_MF_all(:, 4) -  gl_Gran_POTOK_MF(:, 4, j)
                    end if
                end do


                ! Определяем зону в которой находимся

                zone = 1


                call Calc_sourse_MF(qqq, fluid1, SOURSE, zone)  ! Вычисляем источники

                if (l_1 == .TRUE.) then
                    ro3 = qqq(1) - time * POTOK(1) / Volume
                    Q3 = qqq(9) - time * POTOK(9) / Volume
                    if (ro3 <= 0.0_8) then
                        print*, "Ro < 0  1490 ", ro3, gl_Cell_center(:, gr)
                        pause
                    end if
                    u3 = (qqq(1) * qqq(2) - time * POTOK(2) / Volume + time * SOURSE(2, 1)) / ro3
                    v3 = (qqq(1) * qqq(3) - time * POTOK(3) / Volume + time * SOURSE(3, 1)) / ro3
                    w3 = (qqq(1) * qqq(4) - time * POTOK(4) / Volume + time * SOURSE(4, 1)) / ro3
                    p3 = ((  ( qqq(5) / (ggg - 1.0) + 0.5 * qqq(1) * norm2(qqq(2:4))**2 ) &
                        - time * POTOK(5)/ Volume + time * SOURSE(5, 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) ) * (ggg - 1.0)

                    if (p3 <= 0.0_8) then
                        !print*, "p < 0  plasma 2028 ", p3 , gl_Cell_center(:, gr)
                        p3 = 0.000001
                        !pause
                    end if

                    gl_Cell_par(:, gr) = (/ro3, u3, v3, w3, p3, 0.0_8, 0.0_8, 0.0_8, Q3/)

                end if

                ! Теперь посчитаем законы сохранения для остальных жидкостей

                do i = 1, 4
                    if (i == 1 .and. l_1 == .FALSE.) CYCLE
                    
                    if (l_1 == .FALSE.) SOURSE(:, i + 1) = 0.0       ! Не перезаряжаем жидкости внутри мини сферы
                    ro3 = fluid1(1, i) - time * POTOK_MF_all(1, i) / Volume + time * SOURSE(1, i + 1)
                    if (ro3 <= 0.0_8) then
                        print*, "Ro < 0  in ", i,  ro3, gl_Cell_center(:, gr)
                        print*, fluid1(:, i)
                        print*, l_1
                        print*, SOURSE(:, i + 1)
                        print*, time
                        print*, "POTOK ", POTOK_MF_all(:, i)
                        pause
                    end if
                    u3 = (fluid1(1, i) * fluid1(2, i) - time * POTOK_MF_all(2, i) / Volume + time * SOURSE(2, i + 1)) / ro3
                    v3 = (fluid1(1, i) * fluid1(3, i) - time * POTOK_MF_all(3, i) / Volume + time * SOURSE(3, i + 1)) / ro3
                    w3 = (fluid1(1, i) * fluid1(4, i) - time * POTOK_MF_all(4, i) / Volume + time * SOURSE(4, i + 1)) / ro3
                    p3 = ((  ( fluid1(5, i) / (ggg - 1.0) + 0.5 * fluid1(1, i) * norm2(fluid1(2:4, i))**2 ) &
                        - time * POTOK_MF_all(5, i)/ Volume + time * SOURSE(5, i + 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) ) * (ggg - 1.0)
                    if (p3 <= 0.0_8) then
                        !print*, "p < 0  in ", i, p3, gl_Cell_center(:, gr)
                        !print*, fluid1(:, i)
                        !print*, l_1
                        !print*, "Sourse", SOURSE(:, i + 1)
                        !print*, time
                        !print*, "zone", zone, norm2(qqq(2:4))/sqrt(ggg*qqq(5)/qqq(1))
                        !print*, "POTOK ", POTOK_MF_all(:, i)
                        !print*, "qqq", qqq
                        !pause
                        p3 = 0.000001
                        !pause
                    end if

                    gl_Cell_par_MF(:, i, gr) = (/ro3, u3, v3, w3, p3/)
                end do

            end do
            !$omp end do

            !$omp end parallel

        end do

	end subroutine Start_GD_3_inner
	
	subroutine Start_MGD_3_inner_MK(steps)
        ! Считается газовая динамика + параметр Q + мультифлюид (без магнитных полей)
        ! Без атомов и т.д. Без ТВД
        ! Функция паспараллелена на OPEN_MP
        use STORAGE
        use GEO_PARAM
        USE OMP_LIB
        use My_func
        use Solvers
        implicit none

        integer(4), intent(in) :: steps
        integer(4) :: st, gr, ngran, ncell, s1, s2, i, j, k, zone, iter
        integer(4) :: ss1, ss2
        real(8) :: qqq1(9), qqq2(9), qqq(9), n_He  ! Переменные в ячейке
        real(8) :: fluid1(5, par_n_sort), fluid2(5, par_n_sort), MK_kk(5), POTOK2
        real(8) :: dist, dsl, dsc, dsp, distant(3), konvect_(3, 1)
        real(8) :: df1, df2, dff1, dff2, rast(3), rad3, rad4, rad5
        real(8) :: POTOK(9)
        real(8) :: POTOK_MF(5)
        real(8) :: POTOK_MF_all(5, 4)
        real(8) :: time, Volume, TT, U8, rad1, rad2, aa, bb, cc
        real(8) :: SOURSE(5,par_n_sort + 1)  ! Источники массы, импульса и энергии для плазмы и каждого сорта мультифлюида
        real(8) :: Matr(3, 3), Matr2(3, 3), vv(3), kord(3), the1, the2, the3, the4, the5
        real(8) :: qqq11(9), qqq22(9), qqq1_TVD(9), qqq2_TVD(9), qqq11_2(1), qqq22_2(1), qqq1_TVD_2(1), qqq2_TVD_2(1)

        real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3, sks

        logical :: l_1

        ngran = size(gl_all_Gran_inner(:))
        ncell = size(gl_all_Cell_inner(:))
        
        Matr(1,1) = -0.995868
        Matr(1,2) = 0.0177307
        Matr(1,3) = 0.0890681
        Matr(2,1) = 0.0730412
        Matr(2,2) = 0.739193
        Matr(2,3) = 0.669521
        Matr(3,1) = -0.0539675
        Matr(3,2) = 0.67326
        Matr(3,3) = -0.737434
        
        Matr2(1,1) = -0.995868
        Matr2(1,2) = 0.0730412
        Matr2(1,3) = -0.0539675
        Matr2(2,1) = 0.0177307
        Matr2(2,2) = 0.739193
        Matr2(2,3) = 0.67326
        Matr2(3,1) = 0.0890681
        Matr2(3,2) = 0.669521
        Matr2(3,3) = -0.737434

        !call calc_all_Gran()   ! Посчитаем все объёмы, центры, площади и т.д.
        !call Initial_conditions()  ! Задаём начальные условия

        time = 0.000000001

        do st = 1, steps
            ! Делаем цикл по граням и считаем потоки через них
            TT = time
            time = 100000.0
            !if (mod(st, 10) == 0)  print*, "inner",st, TT

            !$omp parallel

            !$omp do private(rad3, rad4, rad5, qqq11, qqq22, qqq1_TVD, qqq2_TVD, qqq11_2, qqq22_2, qqq1_TVD_2, qqq2_TVD_2, df1, df2, dff1, dff2, rast, ss1, ss2, konvect_, distant, vv, kord, the1, the2, the3, the4, the5, POTOK, s1, s2, qqq1, qqq2, dist, dsl, dsp, dsc, rad1, rad2, aa, bb, cc, fluid1, fluid2, POTOK_MF, gr) &
            !$omp & reduction(min:time)
            do iter = 1, ngran
                gr = gl_all_Gran_inner(iter)
                !if(gl_Gran_info(gr) == 0) CYCLE
                POTOK = 0.0
                konvect_(3, 1) = 0.0
                s1 = gl_Gran_neighbour(1, gr)
            s2 = gl_Gran_neighbour(2, gr)
            qqq1 = gl_Cell_par(:, s1)
            konvect_(1, 1) = gl_Cell_par2(1, s1)
            
            distant = gl_Gran_center(:, gr) - gl_Cell_center(:, s1)
            dist = norm2(distant)
    
            !! Попробуем снести плотность пропорционально квадрату
            !rad1 = norm2(gl_Cell_center(:, s1))
            !rad2 = norm2(gl_Gran_center(:, gr))
            !qqq1(1) = qqq1(1) * rad1**2 / rad2**2
            !qqq1(9) = qqq1(9) * rad1**2 / rad2**2
            !qqq1(5) = qqq1(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
            !! Скорости сносим в сферической С.К.
            !call spherical_skorost(gl_Cell_center(3, s1), gl_Cell_center(1, s1), gl_Cell_center(2, s1), &
            !    qqq1(4), qqq1(2), qqq1(3), aa, bb, cc)
            !call dekard_skorost(gl_Gran_center(3, gr), gl_Gran_center(1, gr), gl_Gran_center(2, gr), &
            !    aa, bb, cc, qqq1(4), qqq1(2), qqq1(3))
    
    
            qqq2 = gl_Cell_par(:, s2)
            konvect_(2, 1) = gl_Cell_par2(1, s2)
                
            distant = gl_Gran_center(:, gr) - gl_Cell_center(:, s2)
            dist = min( dist, norm2(distant))
        
            !  ! Попробуем снести плотность пропорционально квадрату
            !  rad1 = norm2(gl_Cell_center(:, s2))
            !  rad2 = norm2(gl_Gran_center(:, gr))
            !  qqq2(1) = qqq2(1) * rad1**2 / rad2**2
            !  qqq2(9) = qqq2(9) * rad1**2 / rad2**2
            !  qqq2(5) = qqq2(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
            !  call spherical_skorost(gl_Cell_center(3, s2), gl_Cell_center(1, s2), gl_Cell_center(2, s2), &
            !      qqq2(4), qqq2(2), qqq2(3), aa, bb, cc)
            !  call dekard_skorost(gl_Gran_center(3, gr), gl_Gran_center(1, gr), gl_Gran_center(2, gr), &
            !      aa, bb, cc, qqq2(4), qqq2(2), qqq2(3))
            
            ss1 = gl_Gran_neighbour_TVD(1, gr)
            ss2 = gl_Gran_neighbour_TVD(2, gr)
            
            if (gl_Cell_number(1, s2) == 1 .and. gl_Cell_number(1, s1) /= 1) write(*, *) "Error cuf_solvers uytuyruyt94847667"
            
            if (gl_Cell_number(1, s1) == 1 .and. gl_Cell_number(1, s2) /= 1) then

                rast = gl_Gran_center(:, gr) - gl_Cell_center(:, s1)
                df1 = norm2(rast)
                rast = gl_Gran_center(:, gr) - gl_Cell_center(:, s2)
                df2 = norm2(rast)
                rast = gl_Gran_center(:, gr) - gl_Cell_center(:, ss2)
                dff2 = norm2(rast)
                
                qqq22 = gl_Cell_par(:, ss2)
                
                
                qqq22_2 = gl_Cell_par2(:, ss2)
                    
                    
                rad1 = norm2(gl_Cell_center(:, s1))                              
                rad2 = norm2(gl_Cell_center(:, s2))     
                rad4 = norm2(gl_Cell_center(:, ss2))                              
                rad5 = norm2(gl_Gran_center(:, gr))
                        
                qqq2_TVD(1) = linear(-dff2, qqq22(1) * rad4**2, -df2, qqq2(1) * rad2**2, df1, qqq1(1) * rad1**2, 0.0_8)/ rad5**2
                qqq2_TVD_2(1) = linear(-dff2, qqq22_2(1) * rad4**2, -df2, konvect_(2, 1) * rad2**2, df1, konvect_(1, 1) * rad1**2, 0.0_8)/ rad5**2
                qqq2_TVD(9) = linear(-dff2, qqq22(9) * rad4**2, -df2, qqq2(9) * rad2**2, df1, qqq1(9) * rad1**2, 0.0_8)/ rad5**2
                qqq2_TVD(5) = linear(-dff2, qqq22(5) * rad4**(2 * ggg), -df2, qqq2(5) * rad2**(2 * ggg), df1,&
                    qqq1(5) * rad1**(2 * ggg), 0.0_8)/ rad5**(2 * ggg)
                        
                ! Переводим скорости в сферическую с.к.
            
                        
                call spherical_skorost(gl_Cell_center(3, s2), gl_Cell_center(1, s2), gl_Cell_center(2, s2), &
                    qqq2(4), qqq2(2), qqq2(3), aa, bb, cc)
                qqq2(4) = aa
                qqq2(2) = bb
                qqq2(3) = cc
                
                        
                call spherical_skorost(gl_Cell_center(3, ss2), gl_Cell_center(1, ss2), gl_Cell_center(2, ss2), &
                    qqq22(4), qqq22(2), qqq22(3), aa, bb, cc)
                qqq22(4) = aa
                qqq22(2) = bb
                qqq22(3) = cc
                
                vv = MATMUL(Matr2, qqq22(6:8))
                kord = MATMUL(Matr2, gl_Cell_center(1:3, ss2))
                the4 = acos(kord(3)/rad4)
                call spherical_skorost(kord(3), kord(1), kord(2), &
                    vv(3), vv(1), vv(2), aa, bb, cc)
                qqq22(8) = aa
                qqq22(6) = bb
                qqq22(7) = cc
                
                vv = MATMUL(Matr2, qqq2(6:8))
                kord = MATMUL(Matr2, gl_Cell_center(1:3, s2))
                the2 = acos(kord(3)/rad2)
                call spherical_skorost(kord(3), kord(1), kord(2), &
                    vv(3), vv(1), vv(2), aa, bb, cc)
                qqq2(8) = aa
                qqq2(6) = bb
                qqq2(7) = cc
                
                vv = MATMUL(Matr2, qqq1(6:8))
                kord = MATMUL(Matr2, gl_Cell_center(1:3, s1))
                the1 = acos(kord(3)/rad1)
                call spherical_skorost(kord(3), kord(1), kord(2), &
                    vv(3), vv(1), vv(2), aa, bb, cc)
                qqq1(8) = aa
                qqq1(6) = bb
                qqq1(7) = cc
                
                kord = MATMUL(Matr2, gl_Gran_center(1:3, gr))
                the5 = acos(kord(3)/rad5)
                        
                do i = 2, 4
                    qqq2_TVD(i) = linear(-dff2, qqq22(i), -df2, qqq2(i), df1, qqq1(i), 0.0_8)
                end do
                        
                do i = 6, 6
                    qqq2_TVD(i) = linear(-dff2, qqq22(i) * rad4**2, -df2, qqq2(i) * rad2**2, df1, qqq1(i) * rad1**2, 0.0_8)/rad5**2
                end do
                
                do i = 7, 7
                    qqq2_TVD(i) = linear(-dff2, qqq22(i) * rad4 / sin(the4), -df2, qqq2(i) * rad2 / sin(the2), df1, qqq1(i) * rad1 / sin(the1), 0.0_8)/rad5 * sin(the5)
                end do
                
                do i = 8, 8
                    qqq2_TVD(i) = linear(-dff2, qqq22(i), -df2, qqq2(i), df1, qqq1(i), 0.0_8)
                end do
                        
                call dekard_skorost(gl_Gran_center(3, gr), gl_Gran_center(1, gr), gl_Gran_center(2, gr), &
                    qqq2_TVD(4), qqq2_TVD(2), qqq2_TVD(3), aa, bb, cc)
                qqq2_TVD(4) = aa
                qqq2_TVD(2) = bb
                qqq2_TVD(3) = cc
                
                
                qqq1(6) = qqq1(6) * rad1**2 / rad5**2
                qqq1(7) = qqq1(7) * rad1 / sin(the1) / rad5 * sin(the5)
                call dekard_skorost(kord(3), kord(1), kord(2), &
                    qqq1(8), qqq1(6), qqq1(7), aa, bb, cc)
                qqq1(8) = aa
                qqq1(6) = bb
                qqq1(7) = cc
                qqq1(6:8) = MATMUL(Matr, qqq1(6:8))
                
                
                call dekard_skorost(kord(3), kord(1), kord(2), &
                    qqq2_TVD(8), qqq2_TVD(6), qqq2_TVD(7), aa, bb, cc)
                qqq2_TVD(8) = aa
                qqq2_TVD(6) = bb
                qqq2_TVD(7) = cc
                qqq2_TVD(6:8) = MATMUL(Matr, qqq2_TVD(6:8))
                
                        
                
                qqq2 = qqq2_TVD
                
                konvect_(2, 1) = qqq2_TVD_2(1)
                
                rad2 = norm2(gl_Gran_center(:, gr))
                qqq1(1) = qqq1(1) * rad1**2 / rad2**2
                konvect_(1, 1) = konvect_(1, 1) * rad1**2 / rad2**2
                qqq1(9) = qqq1(9) * rad1**2 / rad2**2
                qqq1(5) = qqq1(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
                ! Скорости сносим в сферической С.К.
                call spherical_skorost(gl_Cell_center(3, s1), gl_Cell_center(1, s1), gl_Cell_center(2, s1), &
                    qqq1(4), qqq1(2), qqq1(3), aa, bb, cc)
                call dekard_skorost(gl_Gran_center(3, gr), gl_Gran_center(1, gr), gl_Gran_center(2, gr), &
                    aa, bb, cc, qqq1(4), qqq1(2), qqq1(3))
            
            else if (ss1 /= 0 .and. ss2 /= 0) then
                rast = gl_Gran_center(:, gr) - gl_Cell_center(:, s1)
                df1 = norm2(rast)
                rast = gl_Gran_center(:, gr) - gl_Cell_center(:, s2)
                df2 = norm2(rast)
                rast = gl_Gran_center(:, gr) - gl_Cell_center(:, ss1)
                dff1 = norm2(rast)
                rast = gl_Gran_center(:, gr) - gl_Cell_center(:, ss2)
                dff2 = norm2(rast)
                qqq11 = gl_Cell_par(:, ss1)
                qqq22 = gl_Cell_par(:, ss2)
                
                qqq11_2 = gl_Cell_par2(:, ss1)
                qqq22_2 = gl_Cell_par2(:, ss2)
                    
                    
                rad1 = norm2(gl_Cell_center(:, s1))                              
                rad2 = norm2(gl_Cell_center(:, s2))                              
                rad3 = norm2(gl_Cell_center(:, ss1))                              
                rad4 = norm2(gl_Cell_center(:, ss2))                              
                rad5 = norm2(gl_Gran_center(:, gr))
                        
                qqq1_TVD(1) = linear(-dff1, qqq11(1) * rad3**2, -df1, qqq1(1) * rad1**2, df2, qqq2(1) * rad2**2, 0.0_8)/ rad5**2
                qqq1_TVD_2(1) = linear(-dff1, qqq11_2(1) * rad3**2, -df1, konvect_(1, 1) * rad1**2, df2, konvect_(2, 1) * rad2**2, 0.0_8)/ rad5**2
                qqq1_TVD(9) = linear(-dff1, qqq11(9) * rad3**2, -df1, qqq1(9) * rad1**2, df2, qqq2(9) * rad2**2, 0.0_8)/ rad5**2
                qqq1_TVD(5) = linear(-dff1, qqq11(5) * rad3**(2 * ggg), -df1, qqq1(5) * rad1**(2 * ggg), df2,&
                    qqq2(5) * rad2**(2 * ggg), 0.0_8)/ rad5**(2 * ggg)
                        
                qqq2_TVD(1) = linear(-dff2, qqq22(1) * rad4**2, -df2, qqq2(1) * rad2**2, df1, qqq1(1) * rad1**2, 0.0_8)/ rad5**2
                qqq2_TVD_2(1) = linear(-dff2, qqq22_2(1) * rad4**2, -df2, konvect_(2, 1) * rad2**2, df1, konvect_(1, 1) * rad1**2, 0.0_8)/ rad5**2
                qqq2_TVD(9) = linear(-dff2, qqq22(9) * rad4**2, -df2, qqq2(9) * rad2**2, df1, qqq1(9) * rad1**2, 0.0_8)/ rad5**2
                qqq2_TVD(5) = linear(-dff2, qqq22(5) * rad4**(2 * ggg), -df2, qqq2(5) * rad2**(2 * ggg), df1,&
                    qqq1(5) * rad1**(2 * ggg), 0.0_8)/ rad5**(2 * ggg)
                        
                ! Переводим скорости в сферическую с.к.
                call spherical_skorost(gl_Cell_center(3, s1), gl_Cell_center(1, s1), gl_Cell_center(2, s1), &
                    qqq1(4), qqq1(2), qqq1(3), aa, bb, cc)
                qqq1(4) = aa
                qqq1(2) = bb
                qqq1(3) = cc
                
                vv = MATMUL(Matr2, qqq1(6:8))
                kord = MATMUL(Matr2, gl_Cell_center(1:3, s1))
                the1 = acos(kord(3)/rad1)
                call spherical_skorost(kord(3), kord(1), kord(2), &
                    vv(3), vv(1), vv(2), aa, bb, cc)
                qqq1(8) = aa
                qqq1(6) = bb
                qqq1(7) = cc
                        
                call spherical_skorost(gl_Cell_center(3, s2), gl_Cell_center(1, s2), gl_Cell_center(2, s2), &
                    qqq2(4), qqq2(2), qqq2(3), aa, bb, cc)
                qqq2(4) = aa
                qqq2(2) = bb
                qqq2(3) = cc
                
                vv = MATMUL(Matr2, qqq2(6:8))
                kord = MATMUL(Matr2, gl_Cell_center(1:3, s2))
                the2 = acos(kord(3)/rad2)
                call spherical_skorost(kord(3), kord(1), kord(2), &
                    vv(3), vv(1), vv(2), aa, bb, cc)
                qqq2(8) = aa
                qqq2(6) = bb
                qqq2(7) = cc
                        
                call spherical_skorost(gl_Cell_center(3, ss1), gl_Cell_center(1, ss1), gl_Cell_center(2, ss1), &
                    qqq11(4), qqq11(2), qqq11(3), aa, bb, cc)
                qqq11(4) = aa
                qqq11(2) = bb
                qqq11(3) = cc
                
                vv = MATMUL(Matr2, qqq11(6:8))
                kord = MATMUL(Matr2, gl_Cell_center(1:3, ss1))
                the3 = acos(kord(3)/rad3)
                call spherical_skorost(kord(3), kord(1), kord(2), &
                    vv(3), vv(1), vv(2), aa, bb, cc)
                qqq11(8) = aa
                qqq11(6) = bb
                qqq11(7) = cc
                        
                call spherical_skorost(gl_Cell_center(3, ss2), gl_Cell_center(1, ss2), gl_Cell_center(2, ss2), &
                    qqq22(4), qqq22(2), qqq22(3), aa, bb, cc)
                qqq22(4) = aa
                qqq22(2) = bb
                qqq22(3) = cc
                
                vv = MATMUL(Matr2, qqq22(6:8))
                kord = MATMUL(Matr2, gl_Cell_center(1:3, ss2))
                the4 = acos(kord(3)/rad4)
                call spherical_skorost(kord(3), kord(1), kord(2), &
                    vv(3), vv(1), vv(2), aa, bb, cc)
                qqq22(8) = aa
                qqq22(6) = bb
                qqq22(7) = cc
                
                kord = MATMUL(Matr2, gl_Gran_center(1:3, gr))
                the5 = acos(kord(3)/rad5)
                        
                do i = 2, 4
                    qqq1_TVD(i) = linear(-dff1, qqq11(i), -df1, qqq1(i), df2, qqq2(i), 0.0_8)
                    qqq2_TVD(i) = linear(-dff2, qqq22(i), -df2, qqq2(i), df1, qqq1(i), 0.0_8)
                end do
                        
                do i = 6, 6
                    qqq1_TVD(i) = linear(-dff1, qqq11(i) * rad3**2, -df1, qqq1(i) * rad1**2, df2, qqq2(i) * rad2**2, 0.0_8)/rad5**2
                    qqq2_TVD(i) = linear(-dff2, qqq22(i) * rad4**2, -df2, qqq2(i) * rad2**2, df1, qqq1(i) * rad1**2, 0.0_8)/rad5**2
                end do
                
                do i = 7, 7
                    qqq1_TVD(i) = linear(-dff1, qqq11(i) * rad3 / sin(the3), -df1, qqq1(i) * rad1 / sin(the1), df2, qqq2(i) * rad2 / sin(the2), 0.0_8)/rad5 * sin(the5)
                    qqq2_TVD(i) = linear(-dff2, qqq22(i) * rad4 / sin(the4), -df2, qqq2(i) * rad2 / sin(the2), df1, qqq1(i) * rad1 / sin(the1), 0.0_8)/rad5 * sin(the5)
                end do
                
                do i = 8, 8
                    qqq1_TVD(i) = linear(-dff1, qqq11(i), -df1, qqq1(i), df2, qqq2(i), 0.0_8)
                    qqq2_TVD(i) = linear(-dff2, qqq22(i), -df2, qqq2(i), df1, qqq1(i), 0.0_8)
                end do
                        
                call dekard_skorost(gl_Gran_center(3, gr), gl_Gran_center(1, gr), gl_Gran_center(2, gr), &
                    qqq1_TVD(4), qqq1_TVD(2), qqq1_TVD(3), aa, bb, cc)
                qqq1_TVD(4) = aa
                qqq1_TVD(2) = bb
                qqq1_TVD(3) = cc
                call dekard_skorost(gl_Gran_center(3, gr), gl_Gran_center(1, gr), gl_Gran_center(2, gr), &
                    qqq2_TVD(4), qqq2_TVD(2), qqq2_TVD(3), aa, bb, cc)
                qqq2_TVD(4) = aa
                qqq2_TVD(2) = bb
                qqq2_TVD(3) = cc
                
                
                call dekard_skorost(kord(3), kord(1), kord(2), &
                    qqq1_TVD(8), qqq1_TVD(6), qqq1_TVD(7), aa, bb, cc)
                qqq1_TVD(8) = aa
                qqq1_TVD(6) = bb
                qqq1_TVD(7) = cc
                qqq1_TVD(6:8) = MATMUL(Matr, qqq1_TVD(6:8))
                
                call dekard_skorost(kord(3), kord(1), kord(2), &
                    qqq2_TVD(8), qqq2_TVD(6), qqq2_TVD(7), aa, bb, cc)
                qqq2_TVD(8) = aa
                qqq2_TVD(6) = bb
                qqq2_TVD(7) = cc
                qqq2_TVD(6:8) = MATMUL(Matr, qqq2_TVD(6:8))
                        
                
                    
                qqq1 = qqq1_TVD
                qqq2 = qqq2_TVD
                
                konvect_(1, 1) = qqq1_TVD_2(1)
                konvect_(2, 1) = qqq2_TVD_2(1)
            !end if
        else
            rad1 = norm2(gl_Cell_center(:, s1))
            rad2 = norm2(gl_Gran_center(:, gr))
            qqq1(1) = qqq1(1) * rad1**2 / rad2**2
            konvect_(1, 1) = konvect_(1, 1) * rad1**2 / rad2**2
            qqq1(9) = qqq1(9) * rad1**2 / rad2**2
            qqq1(5) = qqq1(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
            ! Скорости сносим в сферической С.К.
            call spherical_skorost(gl_Cell_center(3, s1), gl_Cell_center(1, s1), gl_Cell_center(2, s1), &
                qqq1(4), qqq1(2), qqq1(3), aa, bb, cc)
            call dekard_skorost(gl_Gran_center(3, gr), gl_Gran_center(1, gr), gl_Gran_center(2, gr), &
                aa, bb, cc, qqq1(4), qqq1(2), qqq1(3))
            
            rad1 = norm2(gl_Cell_center(:, s2))
            qqq2(1) = qqq2(1) * rad1**2 / rad2**2
            konvect_(2, 1) = konvect_(2, 1) * rad1**2 / rad2**2
            qqq2(9) = qqq2(9) * rad1**2 / rad2**2
            qqq2(5) = qqq2(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
            call spherical_skorost(gl_Cell_center(3, s2), gl_Cell_center(1, s2), gl_Cell_center(2, s2), &
                qqq2(4), qqq2(2), qqq2(3), aa, bb, cc)
            call dekard_skorost(gl_Gran_center(3, gr), gl_Gran_center(1, gr), gl_Gran_center(2, gr), &
                aa, bb, cc, qqq2(4), qqq2(2), qqq2(3))
        end if
    
    
    
    
            call chlld_Q(3, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK, p_correct_ = .True., konvect_ = konvect_)
            time = min(time, 0.99 * dist/max(dabs(dsl), dabs(dsp)) )   ! REDUCTION
            gl_Gran_POTOK(1:9, gr) = POTOK * gl_Gran_square(gr)
            gl_Gran_POTOK2(1, gr) = konvect_(3, 1) * gl_Gran_square(gr)
            gl_Gran_POTOK(10, gr) = 0.5 * DOT_PRODUCT(gl_Gran_normal(:, gr), qqq1(6:8) + qqq2(6:8)) * gl_Gran_square(gr)
            

            end do
            !$omp end do

            !$omp barrier ! Синхронизация нитей

            !!$omp single
            !    if (mod(st, 50) == 0)  print*, st, time
            !!$omp end single

            ! Теперь цикл по ячейкам
            !$omp do private(n_He, POTOK2, POTOK, Volume, qqq, i, j, ro3, u3, v3, w3, p3, Q3, bx3, by3, bz3, POTOK_MF_all, zone, l_1, fluid1, SOURSE, gr, sks)
            do iter = 1, ncell
                gr = gl_all_Cell_inner(iter)
                !if(gl_Cell_info(gr) == 2) CYCLE
                l_1 = .TRUE.
                if ((gl_Cell_type(gr) == "A" .or. gl_Cell_type(gr) == "B").and.(gl_Cell_number(1, gr) <= 2) ) l_1 = .FALSE. ! 2   ! Не считаем в первых двух ячейках
                POTOK = 0.0
                POTOK2 = 0.0
                sks = 0.0
                SOURSE = 0.0
                Volume = gl_Cell_Volume(gr)
                qqq = gl_Cell_par(:, gr)
                n_He = gl_Cell_par2(1, gr)
                fluid1 = gl_Cell_par_MK(1:5, :, gr)
                MK_kk = gl_Cell_par_MK(6:10, 1, gr)
                !MK_kk = 1.0
                ! Просуммируем потоки через грани
                do i = 1, 6
                    j = gl_Cell_gran(i, gr)
                    if (j == 0) CYCLE
                    if (gl_Gran_neighbour(1, j) == gr) then
                        POTOK = POTOK + gl_Gran_POTOK(1:9, j)
                        POTOK2 = POTOK2 + gl_Gran_POTOK2(1, j)
                        sks = sks + gl_Gran_POTOK(10, j)
                    else
                        POTOK = POTOK - gl_Gran_POTOK(1:9, j)
                        POTOK2 = POTOK2 - gl_Gran_POTOK2(1, j)
                        sks = sks - gl_Gran_POTOK(10, j)
                    end if
                end do
    
    
                ! Определяем зону в которой находимся
    
                zone = 1
    
                qqq(2:4) = qqq(2:4) * (par_chi/par_chi_real)
                qqq(1) = qqq(1) / (par_chi/par_chi_real)**2
                
                ! He
                qqq(1) = qqq(1) - n_He
                qqq(5) = qqq(5) / (8.0 * qqq(1) + 3.0 * n_He) * 8.0 * qqq(1)
                
                
                call Calc_sourse_MF(qqq, fluid1, SOURSE, zone)  ! Вычисляем источники
                
                qqq(5) = qqq(5) * (8.0 * qqq(1) + 3.0 * n_He) / (8.0 * qqq(1))
                qqq(1) = qqq(1) + n_He
                
                qqq(2:4) = qqq(2:4) / (par_chi/par_chi_real)
                qqq(1) = qqq(1) * (par_chi/par_chi_real)**2
                SOURSE(5, 1) = SOURSE(5, 1)/ (par_chi/par_chi_real)
    
                if (l_1 == .TRUE.) then
                    ro3 = qqq(1) - time * POTOK(1) / Volume + time * MK_kk(1)
                    gl_Cell_par2(1, gr) = n_He - time * POTOK2 / Volume
                    Q3 = qqq(9) - time * POTOK(9) / Volume + (qqq(9)/qqq(1)) *  time * MK_kk(1)
                    if (ro3 <= 0.0_8) then
                        write(*, *) "Ro < 0  3242341234 "
                        stop
                    end if
                    
                    if (gl_Cell_par2(1, gr) <= 0.0_8) then
                        write(*, *) "gl_Cell_par2(1, gr) < 0  3242341234 ", gl_Cell_par2(1, gr), n_He
                        stop
                    end if
                    
                    u3 = (qqq(1) * qqq(2) - time * (POTOK(2) + (qqq(6)/cpi4) * sks) / Volume + time * SOURSE(2, 1)) / ro3
                    v3 = (qqq(1) * qqq(3) - time * (POTOK(3) + (qqq(7)/cpi4) * sks) / Volume + time * SOURSE(3, 1)) / ro3
                    w3 = (qqq(1) * qqq(4) - time * (POTOK(4) + (qqq(8)/cpi4) * sks) / Volume + time * SOURSE(4, 1)) / ro3
                    
                    bx3 = qqq(6) - time * (POTOK(6) + qqq(2) * sks) / Volume
                    by3 = qqq(7) - time * (POTOK(7) + qqq(3) * sks) / Volume
                    bz3 = qqq(8) - time * (POTOK(8) + qqq(4) * sks) / Volume
                    
                    p3 = ((  ( qqq(5) / (ggg - 1.0) + 0.5 * qqq(1) * norm2(qqq(2:4))**2 + (qqq(6)**2 + qqq(7)**2 + qqq(8)**2) / 25.13274122871834590768 ) &
                        - time * ( POTOK(5) + (DOT_PRODUCT(qqq(2:4), qqq(6:8))/cpi4) * sks)/ Volume + time * SOURSE(5, 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) - (bx3**2 + by3**2 + bz3**2) / 25.13274122871834590768 ) * (ggg - 1.0)
                    
                    
                    
    
                    if (p3 <= 0.0_8) then
                        !print*, "p < 0  plasma 2028 ", p3 , gl_Cell_center(:, gr)
                        p3 = 0.000001
                        !pause
                    end if
    
                    gl_Cell_par(:, gr) = (/ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3/)
                    
                    !if(gr == 5) then
                    !write(*,*) ro3, u3, v3, w3
                    !write(*,*) POTOK(1), Volume, qqq(1), time
                    !write(*,*) "____________"
                    !	
                    !end if
                    

                end if


            end do
            !$omp end do

            !$omp end parallel

        end do

	end subroutine Start_MGD_3_inner_MK
	
	subroutine Start_MGD_3_inner(steps)
        ! Считается газовая динамика + параметр Q + мультифлюид (без магнитных полей)
        ! Без атомов и т.д. Без ТВД
        ! Функция паспараллелена на OPEN_MP
        use STORAGE
        use GEO_PARAM
        USE OMP_LIB
        use My_func
        use Solvers
        implicit none

        integer(4), intent(in) :: steps
        integer(4) :: st, gr, ngran, ncell, s1, s2, i, j, k, zone, iter
        real(8) :: qqq1(9), qqq2(9), qqq(9)  ! Переменные в ячейке
        real(8) :: fluid1(5, par_n_sort), fluid2(5, par_n_sort)
        real(8) :: dist, dsl, dsc, dsp
        real(8) :: POTOK(9)
        real(8) :: POTOK_MF(5)
        real(8) :: POTOK_MF_all(5, 4)
        real(8) :: time, Volume, TT, U8, rad1, rad2, aa, bb, cc
        real(8) :: SOURSE(5,par_n_sort + 1)  ! Источники массы, импульса и энергии для плазмы и каждого сорта мультифлюида

        real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3, sks

        logical :: l_1

        ngran = size(gl_all_Gran_inner(:))
        ncell = size(gl_all_Cell_inner(:))

        !call calc_all_Gran()   ! Посчитаем все объёмы, центры, площади и т.д.
        !call Initial_conditions()  ! Задаём начальные условия

        time = 0.000000001

        do st = 1, steps
            ! Делаем цикл по граням и считаем потоки через них
            TT = time
            time = 100000.0
            !if (mod(st, 10) == 0)  print*, "inner",st, TT

            !$omp parallel

            !$omp do private(POTOK, s1, s2, qqq1, qqq2, dist, dsl, dsp, dsc, rad1, rad2, aa, bb, cc, fluid1, fluid2, POTOK_MF, gr) &
            !$omp & reduction(min:time)
            do iter = 1, ngran
                gr = gl_all_Gran_inner(iter)
                !if(gl_Gran_info(gr) == 0) CYCLE
                POTOK = 0.0
                s1 = gl_Gran_neighbour(1, gr)
                s2 = gl_Gran_neighbour(2, gr)
                qqq1 = gl_Cell_par(:, s1)
                fluid1 = gl_Cell_par_MF(:, :, s1)   ! Загрузили параметры жидкостей для мультифлюида
                dist = norm2(gl_Gran_center(:, gr) - gl_Cell_center(:, s1))

                ! Попробуем снести плотность пропорционально квадрату
                if(norm2(qqq1(2:4))/sqrt(ggg*qqq1(5)/qqq1(1)) > 2.5) then
                    rad1 = norm2(gl_Cell_center(:, s1))
                    rad2 = norm2(gl_Gran_center(:, gr))
                    qqq1(1) = qqq1(1) * rad1**2 / rad2**2
                    qqq1(9) = qqq1(9) * rad1**2 / rad2**2
                    qqq1(5) = qqq1(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
                    ! Скорости сносим в сферической С.К.
                    call spherical_skorost(gl_Cell_center(3, s1), gl_Cell_center(1, s1), gl_Cell_center(2, s1), &
                        qqq1(4), qqq1(2), qqq1(3), aa, bb, cc)
                    call dekard_skorost(gl_Gran_center(3, gr), gl_Gran_center(1, gr), gl_Gran_center(2, gr), &
                        aa, bb, cc, qqq1(4), qqq1(2), qqq1(3))

                    call spherical_skorost(gl_Cell_center(3, s1), gl_Cell_center(1, s1), gl_Cell_center(2, s1), &
                        fluid1(4, 1), fluid1(2, 1), fluid1(3, 1), aa, bb, cc)
                    call dekard_skorost(gl_Gran_center(3, gr), gl_Gran_center(1, gr), gl_Gran_center(2, gr), &
                        aa, bb, cc, fluid1(4, 1), fluid1(2, 1), fluid1(3, 1))
                end if

                if (s2 >= 1) then
                    !if ( norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int .and. norm2(gl_Cell_center(:, s2)) <= par_R0 * par_R_int) CYCLE
                    qqq2 = gl_Cell_par(:, s2)
                    fluid2 = gl_Cell_par_MF(:, :, s2)   ! Загрузили параметры жидкостей для мультифлюида
                    !dist = min(gl_Cell_dist(s1), gl_Cell_dist(s2))
                    dist = min( dist, norm2(gl_Gran_center(:, gr) - gl_Cell_center(:, s2)))

                    ! Попробуем снести плотность пропорционально квадрату
                    if(norm2(qqq2(2:4))/sqrt(ggg*qqq2(5)/qqq2(1)) > 2.5) then
                        rad1 = norm2(gl_Cell_center(:, s2))
                        rad2 = norm2(gl_Gran_center(:, gr))
                        qqq2(1) = qqq2(1) * rad1**2 / rad2**2
                        qqq2(9) = qqq2(9) * rad1**2 / rad2**2
                        qqq2(5) = qqq2(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
                        call spherical_skorost(gl_Cell_center(3, s2), gl_Cell_center(1, s2), gl_Cell_center(2, s2), &
                            qqq2(4), qqq2(2), qqq2(3), aa, bb, cc)
                        call dekard_skorost(gl_Gran_center(3, gr), gl_Gran_center(1, gr), gl_Gran_center(2, gr), &
                            aa, bb, cc, qqq2(4), qqq2(2), qqq2(3))

                        call spherical_skorost(gl_Cell_center(3, s2), gl_Cell_center(1, s2), gl_Cell_center(2, s2), &
                            fluid2(4, 1), fluid2(2, 1), fluid2(3, 1), aa, bb, cc)
                        call dekard_skorost(gl_Gran_center(3, gr), gl_Gran_center(1, gr), gl_Gran_center(2, gr), &
                            aa, bb, cc, fluid2(4, 1), fluid2(2, 1), fluid2(3, 1))
                    end if

                else  ! В случае граничных ячеек - граничные условия
                    print*, " effefefeIJHGFTYUIJHGUJHGYUIEEE "
                end if


                call chlld_Q(3, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                    0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)
                time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )   ! REDUCTION
                gl_Gran_POTOK(1:9, gr) = POTOK * gl_Gran_square(gr)
                gl_Gran_POTOK(10, gr) = 0.5 * DOT_PRODUCT(gl_Gran_normal(:, gr), qqq1(6:8) + qqq2(6:8)) * gl_Gran_square(gr)

                !print*, time
                !pause

                call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                    0.0_8, fluid1(:, 1), fluid2(:, 1), dsl, dsp, dsc, POTOK_MF)
                time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )
                gl_Gran_POTOK_MF(:, 1, gr) = POTOK_MF * gl_Gran_square(gr)

                !print*, time
                !print*, fluid1(:, 1)
                !print*, fluid2(:, 1)
                !print*, dsl, dsp, dsc
                !pause

                call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                    0.0_8, fluid1(:, 2), fluid2(:, 2), dsl, dsp, dsc, POTOK_MF)
                time = min(time, 0.8 * dist/max(dabs(dsl), dabs(dsp)) )
                gl_Gran_POTOK_MF(:, 2, gr) = POTOK_MF * gl_Gran_square(gr)

                !print*, time
                !pause

                call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                    0.0_8, fluid1(:, 3), fluid2(:, 3), dsl, dsp, dsc, POTOK_MF)
                time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )
                gl_Gran_POTOK_MF(:, 3, gr) = POTOK_MF * gl_Gran_square(gr)

                !print*, time
                !pause

                call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                    0.0_8, fluid1(:, 4), fluid2(:, 4), dsl, dsp, dsc, POTOK_MF)
                time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )
                gl_Gran_POTOK_MF(:, 4, gr) = POTOK_MF * gl_Gran_square(gr)

                !if (dabs(POTOK_MF(1)) > 10.0) then
                !    print*, "POTOK ", POTOK_MF
                !    print*, "fluid1 ",fluid1(:, 4)
                !    print*, "fluid2 ",fluid2(:, 4)
                !    pause
                !end if

            end do
            !$omp end do

            !$omp barrier ! Синхронизация нитей

            !!$omp single
            !    if (mod(st, 50) == 0)  print*, st, time
            !!$omp end single

            ! Теперь цикл по ячейкам
            !$omp do private(POTOK, Volume, qqq, i, j, ro3, u3, v3, w3, p3, Q3, bx3, by3, bz3, POTOK_MF_all, zone, l_1, fluid1, SOURSE, gr, sks)
            do iter = 1, ncell
                gr = gl_all_Cell_inner(iter)
                !if(gl_Cell_info(gr) == 2) CYCLE
                l_1 = .TRUE.
                if ((gl_Cell_type(gr) == "A" .or. gl_Cell_type(gr) == "B").and.(gl_Cell_number(1, gr) <= 2) ) l_1 = .FALSE.    ! Не считаем в первых двух ячейках
                POTOK = 0.0
                sks = 0.0
                SOURSE = 0.0
                POTOK_MF_all = 0.0
                Volume = gl_Cell_Volume(gr)
                qqq = gl_Cell_par(:, gr)
                fluid1 = gl_Cell_par_MF(:, :, gr)
                ! Просуммируем потоки через грани
                do i = 1, 6
                    j = gl_Cell_gran(i, gr)
                    if (j == 0) CYCLE
                    if (gl_Gran_neighbour(1, j) == gr) then
                        POTOK = POTOK + gl_Gran_POTOK(1:9, j)
                        sks = sks + gl_Gran_POTOK(10, j)
                        POTOK_MF_all(:, 1) = POTOK_MF_all(:, 1) +  gl_Gran_POTOK_MF(:, 1, j)
                        POTOK_MF_all(:, 2) = POTOK_MF_all(:, 2) +  gl_Gran_POTOK_MF(:, 2, j)
                        POTOK_MF_all(:, 3) = POTOK_MF_all(:, 3) +  gl_Gran_POTOK_MF(:, 3, j)
                        POTOK_MF_all(:, 4) = POTOK_MF_all(:, 4) +  gl_Gran_POTOK_MF(:, 4, j)
                    else
                        POTOK = POTOK - gl_Gran_POTOK(1:9, j)
                        sks = sks - gl_Gran_POTOK(10, j)
                        POTOK_MF_all(:, 1) = POTOK_MF_all(:, 1) -  gl_Gran_POTOK_MF(:, 1, j)
                        POTOK_MF_all(:, 2) = POTOK_MF_all(:, 2) -  gl_Gran_POTOK_MF(:, 2, j)
                        POTOK_MF_all(:, 3) = POTOK_MF_all(:, 3) -  gl_Gran_POTOK_MF(:, 3, j)
                        POTOK_MF_all(:, 4) = POTOK_MF_all(:, 4) -  gl_Gran_POTOK_MF(:, 4, j)
                    end if
                end do


                ! Определяем зону в которой находимся

                zone = 1
                


                call Calc_sourse_MF(qqq, fluid1, SOURSE, zone)  ! Вычисляем источники

                if (l_1 == .TRUE.) then
                    ro3 = qqq(1) - time * POTOK(1) / Volume
                    Q3 = qqq(9) - time * POTOK(9) / Volume
                    if (ro3 <= 0.0_8) then
                        print*, "Ro < 0  1490 ", ro3, gl_Cell_center(:, gr)
                        pause
                    end if
                    u3 = (qqq(1) * qqq(2) - time * (POTOK(2) + (qqq(6)/cpi4) * sks) / Volume + time * SOURSE(2, 1)) / ro3
                    v3 = (qqq(1) * qqq(3) - time * (POTOK(3) + (qqq(7)/cpi4) * sks) / Volume + time * SOURSE(3, 1)) / ro3
                    w3 = (qqq(1) * qqq(4) - time * (POTOK(4) + (qqq(8)/cpi4) * sks) / Volume + time * SOURSE(4, 1)) / ro3
                    
                    bx3 = qqq(6) - time * (POTOK(6) + qqq(2) * sks) / Volume
                    by3 = qqq(7) - time * (POTOK(7) + qqq(3) * sks) / Volume
                    bz3 = qqq(8) - time * (POTOK(8) + qqq(4) * sks) / Volume
                    
                    p3 = ((  ( qqq(5) / (ggg - 1.0) + 0.5 * qqq(1) * norm2(qqq(2:4))**2 + (qqq(6)**2 + qqq(7)**2 + qqq(8)**2) / 25.13274122871834590768 ) &
                        - time * ( POTOK(5) + (DOT_PRODUCT(qqq(2:4), qqq(6:8))/cpi4) * sks)/ Volume + time * SOURSE(5, 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) - (bx3**2 + by3**2 + bz3**2) / 25.13274122871834590768 ) * (ggg - 1.0)
                    
                    
                    !p3 = ((  ( qqq(5) / (ggg - 1.0) + 0.5 * qqq(1) * norm2(qqq(2:4))**2 ) &
                    !    - time * ( POTOK(5) + (DOT_PRODUCT(qqq(2:4), qqq(6:8))/cpi4) * sks)/ Volume + time * SOURSE(5, 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) ) * (ggg - 1.0)
                    

                    if (p3 <= 0.0_8) then
                        !print*, "p < 0  plasma 2028 ", p3 , gl_Cell_center(:, gr)
                        p3 = 0.000001
                        !pause
                    end if

                    gl_Cell_par(:, gr) = (/ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3/)
                    
                    !if(gr == 5) then
                    !write(*,*) ro3, u3, v3, w3
                    !write(*,*) POTOK(1), Volume, qqq(1), time
                    !write(*,*) "____________"
                    !	
                    !end if
                    

                end if

                ! Теперь посчитаем законы сохранения для остальных жидкостей

                do i = 1, 4
                    if (i == 1 .and. l_1 == .FALSE.) CYCLE
                    
                    if (l_1 == .FALSE.) SOURSE(:, i + 1) = 0.0       ! Не перезаряжаем жидкости внутри мини сферы
                    ro3 = fluid1(1, i) - time * POTOK_MF_all(1, i) / Volume + time * SOURSE(1, i + 1)
                    if (ro3 <= 0.0_8) then
                        print*, "Ro < 0  in ", i,  ro3, gl_Cell_center(:, gr)
                        print*, fluid1(:, i)
                        print*, l_1
                        print*, SOURSE(:, i + 1)
                        print*, time
                        print*, "POTOK ", POTOK_MF_all(:, i)
                        pause
                    end if
                    u3 = (fluid1(1, i) * fluid1(2, i) - time * POTOK_MF_all(2, i) / Volume + time * SOURSE(2, i + 1)) / ro3
                    v3 = (fluid1(1, i) * fluid1(3, i) - time * POTOK_MF_all(3, i) / Volume + time * SOURSE(3, i + 1)) / ro3
                    w3 = (fluid1(1, i) * fluid1(4, i) - time * POTOK_MF_all(4, i) / Volume + time * SOURSE(4, i + 1)) / ro3
                    p3 = ((  ( fluid1(5, i) / (ggg - 1.0) + 0.5 * fluid1(1, i) * norm2(fluid1(2:4, i))**2 ) &
                        - time * POTOK_MF_all(5, i)/ Volume + time * SOURSE(5, i + 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) ) * (ggg - 1.0)
                    if (p3 <= 0.0_8) then
                        !print*, "p < 0  in ", i, p3, gl_Cell_center(:, gr)
                        !print*, fluid1(:, i)
                        !print*, l_1
                        !print*, "Sourse", SOURSE(:, i + 1)
                        !print*, time
                        !print*, "zone", zone, norm2(qqq(2:4))/sqrt(ggg*qqq(5)/qqq(1))
                        !print*, "POTOK ", POTOK_MF_all(:, i)
                        !print*, "qqq", qqq
                        !pause
                        p3 = 0.000001
                        !pause
                    end if

                    gl_Cell_par_MF(:, i, gr) = (/ro3, u3, v3, w3, p3/)
                end do

            end do
            !$omp end do

            !$omp end parallel

        end do

	end subroutine Start_MGD_3_inner
	
    subroutine Start_GD_move()
        ! Функция газовой динамики + мультифлюид + подвижная сетка
        ! Это управляющая функция, в отличие от предыдущих версий, здесь показано как правильно вызывать процедуры для движения сетки
        use GEO_PARAM
        use STORAGE
        use My_func
        use Solvers
        implicit none
        integer :: step, now, now2, step2
        integer(4) :: st, gr, ngran, ncell, s1, s2, i, j, k, zone, iter
        real(8) :: qqq1(9), qqq2(9), qqq(9)  ! Переменные в ячейке
        real(8) :: fluid1(5, par_n_sort), fluid2(5, par_n_sort)
        real(8) :: dist, dsl, dsc, dsp
        real(8) :: POTOK(9)
        real(8) :: POTOK_MF(5)
        real(8) :: POTOK_MF_all(5, 4)
        real(8) :: time, Volume, Volume2, TT, U8, rad1, rad2, aa, bb, cc, wc
        real(8) :: SOURSE(5,par_n_sort + 1)  ! Источники массы, импульса и энергии для плазмы и каждого сорта мультифлюида

        real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3

        logical :: l_1


        ! Сначала подготовим все массивы для правильного движения
        if (allocated(gl_Vx) == .False.) then   ! Если движение запускается впервый раз, то нужно выделить память под используемые массивы и задать значения
            allocate(gl_Vx(size(gl_x)))
            allocate(gl_Vy(size(gl_x)))
            allocate(gl_Vz(size(gl_x)))
            allocate(gl_Point_num(size(gl_x)))
            allocate(gl_x2(size(gl_x), 2))
            allocate(gl_y2(size(gl_x), 2))
            allocate(gl_z2(size(gl_x), 2))
            allocate(gl_Cell_Volume2(size(gl_Cell_Volume), 2))
            allocate(gl_Gran_normal2(  size(gl_Gran_normal(:, 1)), size(gl_Gran_normal(1, :)), 2  ))
            allocate(gl_Gran_center2(  size(gl_Gran_center(:, 1)), size(gl_Gran_center(1, :)), 2  ))
            allocate(gl_Cell_center2(size(gl_Cell_center(:, 1)), size(gl_Cell_center(1, :)), 2))
            allocate(gl_Gran_square2(size(gl_Gran_square), 2))
            
            ! Начальная инициализвация
            gl_Vx = 0.0
            gl_Vy = 0.0
            gl_Vz = 0.0
            gl_Point_num = 0
            gl_x2(:, 1) = gl_x
            gl_x2(:, 2) = gl_x
            gl_y2(:, 1) = gl_y
            gl_y2(:, 2) = gl_y
            gl_z2(:, 1) = gl_z
            gl_z2(:, 2) = gl_z
            gl_Cell_Volume2(:, 1) = gl_Cell_Volume
            gl_Cell_Volume2(:, 2) = gl_Cell_Volume
            gl_Gran_normal2(:, :, 1) = gl_Gran_normal
            gl_Gran_normal2(:, :, 2) = gl_Gran_normal
            gl_Gran_center2(:, :, 1) = gl_Gran_center
            gl_Gran_center2(:, :, 2) = gl_Gran_center
            gl_Cell_center2(:, :, 1) = gl_Cell_center
            gl_Cell_center2(:, :, 2) = gl_Cell_center
            gl_Gran_square2(:, 1) = gl_Gran_square
            gl_Gran_square2(:, 2) = gl_Gran_square
        end if
        
        
        
        ! Нужно посчитать новые нормали граней, их площади, центры, объёмы ячеек, центры ячеек (зачем?)

        ! Запускаем глобальный цикл
        now = 2                           ! Какие параметры сейчас будут считаться (1 или 2). Они меняются по очереди
        time = 0.00000000001               ! Начальная инициализация шага по времени (в данной программе это не нужно, так как шаг вычисляется налету)
        do step = 1, 20000 * 2 * 4 * 4   !    ! Нужно чтобы это число было чётным!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            if (mod(step, 100) == 0) print*, "Step = ", step , "  step_time = ", time
            now2 = now
            now = mod(now, 2) + 1
            
            TT = time
            time = 100000.0
            
            ! Вычисляем новые скорости управляющих узлов
            call Calc_move(now)   ! Записали скорости каждого узла в этот узел для последующего движения
            
            
            ! Двигаем все узлы сетки в соответствии с расчитанными скоростями в предыдущей функции
            call Move_all(now, TT)   
            
            call calc_all_Gran_move(now2)   ! Расчитываются новые объёмы\площади\нормали и т.д.
            
            ngran = size(gl_all_Gran_inner(:))
            ncell = size(gl_all_Cell_inner(:))
            
            !$omp parallel

            ! Теперь по основным ячейкам
            ngran = size(gl_all_Gran(1, :))
            ncell = size(gl_all_Cell(1, :))
            
            !$omp do private(POTOK, s1, s2, qqq1, qqq2, dist, dsl, dsp, dsc, rad1, rad2, aa, bb, cc, fluid1, fluid2, POTOK_MF, wc) &
            !$omp & reduction(min:time)
            do gr = 1, ngran
                if(gl_Gran_info(gr) == 2) CYCLE
                POTOK = 0.0
                s1 = gl_Gran_neighbour(1, gr)
                s2 = gl_Gran_neighbour(2, gr)
                qqq1 = gl_Cell_par(:, s1)
                fluid1 = gl_Cell_par_MF(:, :, s1)   ! Загрузили параметры жидкостей для мультифлюида

                ! Попробуем снести плотность пропорционально квадрату
                if(norm2(qqq1(2:4))/sqrt(ggg*qqq1(5)/qqq1(1)) > 2.2) then
                    rad1 = norm2(gl_Cell_center2(:, s1, now))
                    rad2 = norm2(gl_Gran_center2(:, gr, now))
                    qqq1(1) = qqq1(1) * rad1**2 / rad2**2
                    qqq1(9) = qqq1(9) * rad1**2 / rad2**2
                    qqq1(5) = qqq1(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
                    ! Скорости сносим в сферической С.К.
                    call spherical_skorost(gl_Cell_center2(1, s1, now), gl_Cell_center2(2, s1, now), gl_Cell_center2(3, s1, now), &  
                        qqq1(2), qqq1(3), qqq1(4), aa, bb, cc)
                    call dekard_skorost(gl_Gran_center2(1, gr, now), gl_Gran_center2(2, gr, now), gl_Gran_center2(3, gr, now), &  
                        aa, bb, cc, qqq1(2), qqq1(3), qqq1(4))

                    call spherical_skorost(gl_Cell_center2(1, s1, now), gl_Cell_center2(2, s1, now), gl_Cell_center2(3, s1, now), & 
                        fluid1(2, 1), fluid1(3, 1), fluid1(4, 1), aa, bb, cc)
                    call dekard_skorost(gl_Gran_center2(1, gr, now), gl_Gran_center2(2, gr, now), gl_Gran_center2(3, gr, now), & 
                        aa, bb, cc, fluid1(2, 1), fluid1(3, 1), fluid1(4, 1))
                end if

                if (s2 >= 1) then
                    !if ( norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int .and. norm2(gl_Cell_center(:, s2)) <= par_R0 * par_R_int) CYCLE
                    qqq2 = gl_Cell_par(:, s2)
                    fluid2 = gl_Cell_par_MF(:, :, s2)   ! Загрузили параметры жидкостей для мультифлюида
                    dist = min(gl_Cell_dist(s1), gl_Cell_dist(s2))                                              

                    ! Попробуем снести плотность пропорционально квадрату
                    if(norm2(qqq2(2:4))/sqrt(ggg*qqq2(5)/qqq2(1)) > 2.2) then 
                        rad1 = norm2(gl_Cell_center2(:, s2, now))                              
                        rad2 = norm2(gl_Gran_center2(:, gr, now))
                        qqq2(1) = qqq2(1) * rad1**2 / rad2**2
                        qqq2(9) = qqq2(9) * rad1**2 / rad2**2
                        qqq2(5) = qqq2(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
                        call spherical_skorost(gl_Cell_center2(1, s2, now), gl_Cell_center2(2, s2, now), gl_Cell_center2(3, s2, now), &
                            qqq2(2), qqq2(3), qqq2(4), aa, bb, cc)
                        call dekard_skorost(gl_Gran_center2(1, gr, now), gl_Gran_center2(2, gr, now), gl_Gran_center2(3, gr, now), &
                            aa, bb, cc, qqq2(2), qqq2(3), qqq2(4))

                        call spherical_skorost(gl_Cell_center2(1, s2, now), gl_Cell_center2(2, s2, now), gl_Cell_center2(3, s2, now), &
                            fluid2(2, 1), fluid2(3, 1), fluid2(4, 1), aa, bb, cc)
                        call dekard_skorost(gl_Gran_center2(1, gr, now), gl_Gran_center2(2, gr, now), gl_Gran_center2(3, gr, now), &
                            aa, bb, cc, fluid2(2, 1), fluid2(3, 1), fluid2(4, 1))
                    end if

                else  ! В случае граничных ячеек - граничные условия
                    !if (norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int) CYCLE
                    if(s2 == -1) then  ! Набегающий поток
                        dist = gl_Cell_dist(s1)
                        qqq2 = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 1.0_8, 0.0_8, 0.0_8, 0.0_8, 100.0_8/)
                        fluid2(:, 1) = fluid1(:, 1)
                        fluid2(:, 2) = fluid1(:, 2)
                        fluid2(:, 3) = fluid1(:, 3)
                        fluid2(:, 4) = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 0.5_8/)
                    else  ! Здесь нужны мягкие условия
                        dist = gl_Cell_dist(s1)
                        qqq2 = qqq1
                        fluid2 = fluid1
                        qqq2(5) = 1.0_8
                        if(qqq2(2) > par_Velosity_inf) then
                            qqq2(2) = par_Velosity_inf ! Отсос жидкости
                        end if

                    end if
                end if
                
                ! Нужно вычислить скорость движения грани
                wc = DOT_PRODUCT((gl_Gran_center2(:, gr, now2) -  gl_Gran_center2(:, gr, now))/TT, gl_Gran_normal2(:, gr, now))
                
                call chlld_Q(2, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
                    wc, qqq1, qqq2, dsl, dsp, dsc, POTOK)
                time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )   ! REDUCTION
                gl_Gran_POTOK(1:9, gr) = POTOK * gl_Gran_square2(gr, now)

                call chlld_gd(1, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
                    wc, fluid1(:, 1), fluid2(:, 1), dsl, dsp, dsc, POTOK_MF)
                time = min(time, 0.9 * dist/(max(dabs(dsl), dabs(dsp)) + dabs(wc)) )
                gl_Gran_POTOK_MF(:, 1, gr) = POTOK_MF * gl_Gran_square2(gr, now)

                call chlld_gd(1, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
                    wc, fluid1(:, 2), fluid2(:, 2), dsl, dsp, dsc, POTOK_MF)
                time = min(time, 0.9 * dist/(max(dabs(dsl), dabs(dsp)) + dabs(wc)) )
                gl_Gran_POTOK_MF(:, 2, gr) = POTOK_MF * gl_Gran_square2(gr, now)
                
                call chlld_gd(1, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
                    wc, fluid1(:, 3), fluid2(:, 3), dsl, dsp, dsc, POTOK_MF)
                time = min(time, 0.9 * dist/(max(dabs(dsl), dabs(dsp)) + dabs(wc)) )
                gl_Gran_POTOK_MF(:, 3, gr) = POTOK_MF * gl_Gran_square2(gr, now)

                call chlld_gd(1, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
                    wc, fluid1(:, 4), fluid2(:, 4), dsl, dsp, dsc, POTOK_MF)
                time = min(time, 0.9 * dist/(max(dabs(dsl), dabs(dsp)) + dabs(wc)) )
                gl_Gran_POTOK_MF(:, 4, gr) = POTOK_MF * gl_Gran_square2(gr, now)


            end do
            !$omp end do



            !$omp barrier ! Синхронизация нитей

            !!$omp single
            !    if (mod(st, 5) == 0)  print*, st, time
            !!$omp end single

            ! Теперь цикл по ячейкам
            !$omp do private(POTOK, Volume, Volume2, qqq, i, j, ro3, u3, v3, w3, p3, Q3, POTOK_MF_all, zone, l_1, fluid1, SOURSE)
            do gr = 1, ncell
                if(gl_Cell_info(gr) == 0) CYCLE
                l_1 = .TRUE.
                if (norm2(gl_Cell_center2(:, gr, now)) <= par_R0 + (par_R_character - par_R0) * (3.0_8/par_n_TS)**par_kk1) l_1 = .FALSE.    ! Не считаем внутри сферы
                POTOK = 0.0
                SOURSE = 0.0
                POTOK_MF_all = 0.0
                Volume = gl_Cell_Volume2(gr, now)
                Volume2 = gl_Cell_Volume2(gr, now2)
                qqq = gl_Cell_par(:, gr)
                fluid1 = gl_Cell_par_MF(:, :, gr)
                ! Просуммируем потоки через грани
                do i = 1, 6
                    j = gl_Cell_gran(i, gr)
                    if (j == 0) CYCLE
                    if (gl_Gran_neighbour(1, j) == gr) then
                        POTOK = POTOK + gl_Gran_POTOK(1:9, j)
                        POTOK_MF_all(:, 1) = POTOK_MF_all(:, 1) +  gl_Gran_POTOK_MF(:, 1, j)
                        POTOK_MF_all(:, 2) = POTOK_MF_all(:, 2) +  gl_Gran_POTOK_MF(:, 2, j)
                        POTOK_MF_all(:, 3) = POTOK_MF_all(:, 3) +  gl_Gran_POTOK_MF(:, 3, j)
                        POTOK_MF_all(:, 4) = POTOK_MF_all(:, 4) +  gl_Gran_POTOK_MF(:, 4, j)
                    else
                        POTOK = POTOK - gl_Gran_POTOK(1:9, j)
                        POTOK_MF_all(:, 1) = POTOK_MF_all(:, 1) -  gl_Gran_POTOK_MF(:, 1, j)
                        POTOK_MF_all(:, 2) = POTOK_MF_all(:, 2) -  gl_Gran_POTOK_MF(:, 2, j)
                        POTOK_MF_all(:, 3) = POTOK_MF_all(:, 3) -  gl_Gran_POTOK_MF(:, 3, j)
                        POTOK_MF_all(:, 4) = POTOK_MF_all(:, 4) -  gl_Gran_POTOK_MF(:, 4, j)
                    end if
                end do


                ! Определяем зону в которой находимся
                if(qqq(9)/qqq(1) < 50.0) then

                    ! Перенормируем параметры
                    !qqq(2:4) = qqq(2:4) * (par_chi/par_chi_real)
                    !qqq(1) = qqq(1) / (par_chi/par_chi_real)**2

                    if(norm2(qqq(2:4))/sqrt(ggg*qqq(5)/qqq(1)) > 2.3) then
                        zone = 1
                    else
                        zone = 2
                    end if
                else
                    if(norm2(qqq(2:4))/sqrt(ggg*qqq(5)/qqq(1)) > 1.1) then
                        zone = 4
                    else
                        zone = 3
                    end if
                end if

                ! Перенормируем первую жидкость
                !fluid1(2:4, 1) = fluid1(2:4, 1) * (par_chi/par_chi_real)
                !fluid1(1, 1) = fluid1(1, 1) / (par_chi/par_chi_real)**2
                !
                !fluid1(2:4, 2) = fluid1(2:4, 2) * (3.0_8)
                !fluid1(1, 2) = fluid1(1, 2) / (3.0_8)**2

                call Calc_sourse_MF(qqq, fluid1, SOURSE, zone)  ! Вычисляем источники

                ! Перенормируем первую жидкость обратно
                !fluid1(2:4, 1) = fluid1(2:4, 1) / (par_chi/par_chi_real)
                !fluid1(1, 1) = fluid1(1, 1) * (par_chi/par_chi_real)**2
                !SOURSE(5, 2) = SOURSE(5, 2)/ (par_chi/par_chi_real)
                !SOURSE(1, 2) = SOURSE(1, 2)* (par_chi/par_chi_real)
                !
                !fluid1(2:4, 2) = fluid1(2:4, 2) / (3.0_8)
                !fluid1(1, 2) = fluid1(1, 2) * (3.0_8)**2
                !SOURSE(5, 3) = SOURSE(5, 3)/ (3.0_8)
                !SOURSE(1, 3) = SOURSE(1, 3)* (3.0_8)

                ! Перенормируем обратно
                !if(zone == 1 .or. zone == 2) then
                !    qqq(2:4) = qqq(2:4) / (par_chi/par_chi_real)
                !    qqq(1) = qqq(1) * (par_chi/par_chi_real)**2
                !    SOURSE(5, 1) = SOURSE(5, 1)/ (par_chi/par_chi_real)
                !end if


                if (l_1 == .TRUE.) then
                    ro3 = qqq(1)* Volume / Volume2 - time * POTOK(1) / Volume2
                    Q3 = qqq(9)* Volume / Volume2 - time * POTOK(9) / Volume2
                    if (ro3 <= 0.0_8) then
                        print*, "Ro < 0  1490 ", ro3, gl_Cell_center(:, gr)
                        pause
                    end if
                    u3 = (qqq(1) * qqq(2)* Volume / Volume2 - time * POTOK(2) / Volume2 + time * SOURSE(2, 1)) / ro3
                    v3 = (qqq(1) * qqq(3)* Volume / Volume2 - time * POTOK(3) / Volume2 + time * SOURSE(3, 1)) / ro3
                    w3 = (qqq(1) * qqq(4)* Volume / Volume2 - time * POTOK(4) / Volume2 + time * SOURSE(4, 1)) / ro3
                    p3 = ((  ( qqq(5) / (ggg - 1.0) + 0.5 * qqq(1) * norm2(qqq(2:4))**2 )* Volume / Volume2 &
                        - time * POTOK(5)/ Volume2 + time * SOURSE(5, 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) ) * (ggg - 1.0)

                    if (p3 <= 0.0_8) then
                        !print*, "p < 0  plasma 2028 ", p3 , gl_Cell_center(:, gr)
                        p3 = 0.000001
                        !pause
                    end if

                    gl_Cell_par(:, gr) = (/ro3, u3, v3, w3, p3, 0.0_8, 0.0_8, 0.0_8, Q3/)

                end if

                ! Теперь посчитаем законы сохранения для остальных жидкостей

                do i = 1, 4
                    if (i == 1 .and. l_1 == .FALSE.) CYCLE       ! Пропускаем внутреннюю сферу для сорта 1
                    if (l_1 == .FALSE.) SOURSE(:, i + 1) = 0.0       ! Пропускаем внутреннюю сферу для сорта 1
                    ro3 = fluid1(1, i)* Volume / Volume2 - time * POTOK_MF_all(1, i) / Volume2 + time * SOURSE(1, i + 1)
                    if (ro3 <= 0.0_8) then
                        print*, "Ro < 0  in ", i,  ro3
                        print*, "gl_Cell_center(:, gr)", gl_Cell_center(:, gr)
                        print*, "gl_Cell_info(gr)", gl_Cell_info(gr)
                        print*, fluid1(:, i)
                        print*, l_1
                        print*, SOURSE(:, i + 1)
                        print*, time
                        print*, "POTOK ", POTOK_MF_all(:, i)
                        pause
                    end if
                    u3 = (fluid1(1, i) * fluid1(2, i)* Volume / Volume2 - time * POTOK_MF_all(2, i) / Volume2 + time * SOURSE(2, i + 1)) / ro3
                    v3 = (fluid1(1, i) * fluid1(3, i)* Volume / Volume2 - time * POTOK_MF_all(3, i) / Volume2 + time * SOURSE(3, i + 1)) / ro3
                    w3 = (fluid1(1, i) * fluid1(4, i)* Volume / Volume2 - time * POTOK_MF_all(4, i) / Volume2 + time * SOURSE(4, i + 1)) / ro3
                    p3 = ((  ( fluid1(5, i) / (ggg - 1.0) + 0.5 * fluid1(1, i) * norm2(fluid1(2:4, i))**2 )* Volume / Volume2 &
                        - time * POTOK_MF_all(5, i)/ Volume2 + time * SOURSE(5, i + 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) ) * (ggg - 1.0)
                    if (p3 <= 0.0_8) then
                        !print*, "p < 0  in ", i, p3, gl_Cell_center(:, gr)
                        !print*, fluid1(:, i)
                        !print*, l_1
                        !print*, "Sourse", SOURSE(:, i + 1)
                        !print*, time
                        !print*, "zone", zone, norm2(qqq(2:4))/sqrt(ggg*qqq(5)/qqq(1))
                        !print*, "POTOK ", POTOK_MF_all(:, i)
                        !print*, "qqq", qqq
                        !pause
                        p3 = 0.000001
                        !pause
                    end if

                    gl_Cell_par_MF(:, i, gr) = (/ro3, u3, v3, w3, p3/)
                end do

            end do
            !$omp end do
            
            !$omp end parallel
            
            gl_x = gl_x2(:, now2)
            gl_y = gl_y2(:, now2)
            gl_z = gl_z2(:, now2)
            gl_Cell_Volume = gl_Cell_Volume2(:, now2)
            gl_Gran_normal = gl_Gran_normal2(:, :, now2)
            gl_Gran_center = gl_Gran_center2(:, :, now2)
            gl_Cell_center = gl_Cell_center2(:, :, now2)
            gl_Gran_square = gl_Gran_square2(:, now2)
            
            if (mod(step, 100) == 0) then
                call Initial_conditions()
            end if
            
            call Start_GD_3_inner(3)
            
            if (mod(step, 5000) == 0) then
                call Print_surface_2D()
                call Print_Setka_2D()
                call Print_par_2D()
            end if
            
            
            
            
        end do


            gl_x = gl_x2(:, 2)
            gl_y = gl_y2(:, 2)
            gl_z = gl_z2(:, 2)
            gl_Cell_Volume = gl_Cell_Volume2(:, 2)
            gl_Gran_normal = gl_Gran_normal2(:, :, 2)
            gl_Gran_center = gl_Gran_center2(:, :, 2)
            gl_Cell_center = gl_Cell_center2(:, :, 2)
            gl_Gran_square = gl_Gran_square2(:, 2)



	end subroutine Start_GD_move
	
	subroutine Start_MGD_move_MK()
        ! Функция газовой динамики + мультифлюид + подвижная сетка
        ! Это управляющая функция, в отличие от предыдущих версий, здесь показано как правильно вызывать процедуры для движения сетки
        use GEO_PARAM
        use STORAGE
        USE ieee_arithmetic 
        use My_func
        USE OMP_LIB
        use Solvers
        use cgod
        implicit none
        integer :: step, now, now2, step2, min_sort, ijk, ijk2, min_sort2, N2, N3
        integer(4) :: st, gr, ngran, ncell, s1, s2, i, j, k, zone, iter, metod, mincell, ss1, ss2
        real(8) :: qqq1(9), qqq2(9), qqq(9), distant(3), n_He  ! Переменные в ячейке
        real(8) :: fluid1(5, par_n_sort), fluid2(5, par_n_sort), MK_kk(5)
        real(8) :: dist, dsl, dsc, dsp, start_time, end_time, rast(3), df1, df2, dff1, dff2, qqq11(9), qqq22(9)
        real(8) :: POTOK(9), qqq1_TVD(9), qqq2_TVD(9), POTOK2
        real(8) :: POTOK_MF(5), konvect_(3, 1)
        real(8) :: POTOK_MF_all(5, 4)
        logical :: tvd1, tvd2, tvd3, tvd4  ! Нужно ли делать особый снос в гиперзвуковом источнике
        real(8) :: time, Volume, Volume2, TT, U8, rad1, rad2, aa, bb, cc, wc, sks, loc_time, loc_time2, rr, y, z
        real(8) :: SOURSE(5, par_n_sort + 1)  ! Источники массы, импульса и энергии для плазмы и каждого сорта мультифлюида
        real(8) :: qqq11_2(1), qqq22_2(1), qqq1_TVD_2(1), qqq2_TVD_2(1)
        real(8) :: rad3, rad4, rad5
        integer(4) :: kdir, KOBL, idgod

        real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3

        logical :: l_1, null_bn

        min_sort = 0

        ! Сначала подготовим все массивы для правильного движения
        if (allocated(gl_Vx) == .False.) then   ! Если движение запускается впервый раз, то нужно выделить память под используемые массивы и задать значения
            allocate(gl_Vx(size(gl_x)))
            allocate(gl_Vy(size(gl_x)))
            allocate(gl_Vz(size(gl_x)))
            allocate(gl_Point_num(size(gl_x)))
            allocate(gl_x2(size(gl_x), 2))
            allocate(gl_y2(size(gl_x), 2))
            allocate(gl_z2(size(gl_x), 2))
            allocate(gl_Cell_Volume2(size(gl_Cell_Volume), 2))
            allocate(gl_Gran_normal2(  size(gl_Gran_normal(:, 1)), size(gl_Gran_normal(1, :)), 2  ))
            allocate(gl_Gran_center2(  size(gl_Gran_center(:, 1)), size(gl_Gran_center(1, :)), 2  ))
            allocate(gl_Cell_center2(size(gl_Cell_center(:, 1)), size(gl_Cell_center(1, :)), 2))
            allocate(gl_Gran_square2(size(gl_Gran_square), 2))
            
            ! Начальная инициализация
            gl_Vx = 0.0
            gl_Vy = 0.0
            gl_Vz = 0.0
            gl_Point_num = 0
            gl_x2(:, 1) = gl_x
            gl_x2(:, 2) = gl_x
            gl_y2(:, 1) = gl_y
            gl_y2(:, 2) = gl_y
            gl_z2(:, 1) = gl_z
            gl_z2(:, 2) = gl_z
            gl_Cell_Volume2(:, 1) = gl_Cell_Volume
            gl_Cell_Volume2(:, 2) = gl_Cell_Volume
            gl_Gran_normal2(:, :, 1) = gl_Gran_normal
            gl_Gran_normal2(:, :, 2) = gl_Gran_normal
            gl_Gran_center2(:, :, 1) = gl_Gran_center
            gl_Gran_center2(:, :, 2) = gl_Gran_center
            gl_Cell_center2(:, :, 1) = gl_Cell_center
            gl_Cell_center2(:, :, 2) = gl_Cell_center
            gl_Gran_square2(:, 1) = gl_Gran_square
            gl_Gran_square2(:, 2) = gl_Gran_square
        end if
        
        
        
        start_time = omp_get_wtime()
        
        mincell = 1
        
        !call Start_MGD_3_inner(10000)
        
        
        ! Запускаем глобальный цикл
        now = 2                           ! Какие параметры сейчас будут считаться (1 или 2). Они меняются по очереди
        time = 0.00002_8               ! Начальная инициализация шага по времени 
        do step = 1, 100  !    ! Нужно чтобы это число было чётным!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            if (mod(step, 1) == 0) then
                print*, "Step = ", step , "  step_time = ", time, "  mingran = ", mincell, & 
                "  gran Center = ", gl_Gran_center2(:, mincell, now), "  minsort = ", min_sort, "  cell = ", gl_Gran_neighbour(1, mincell), gl_Gran_neighbour(2, mincell) 
                !print*, gl_Gran_normal2(:, mincell, now)
                print*, "  ******************************************************************************** "
                print*, "  ******************************************************************************** "
                print*, "  ******************************************************************************** "
                print*, "  ******************************************************************************** "
                print*, "  ******************************************************************************** "
                print*, "  ******************************************************************************** "
            end if
                
            now2 = now
            now = mod(now, 2) + 1
            
            TT = time
            time = 100000.0
            !TT = 0.0
            ! Вычисляем новые скорости управляющих узлов
            
            call Calc_move(now)   ! Записали скорости каждого узла в этот узел для последующего движения
            
            
            ! Двигаем все узлы сетки в соответствии с расчитанными скоростями в предыдущей функции
            call Move_all(now, TT) 
            
            call calc_all_Gran_move(now2)   ! Расчитываются новые объёмы\площади\нормали и т.д.
            
            
            ! Теперь по основным граням
            ngran = size(gl_all_Gran(1, :))
            ncell = size(gl_all_Cell(1, :))
            
            !$omp parallel
            
        !$omp do private(n_He, kdir, KOBL, idgod, rad3, rad4, rad5, qqq11_2, qqq22_2, qqq1_TVD_2, qqq2_TVD_2, tvd1, tvd2, tvd3, tvd4, distant, konvect_, metod, POTOK, ss1, ss2, qqq1_TVD, qqq2_TVD, rast, df1, df2, dff1, dff2, qqq11, qqq22, s1, s2, qqq1, qqq2, dist, dsl, dsp, dsc, rad1, rad2, aa, bb, cc, fluid1, fluid2, POTOK_MF, wc, null_bn, loc_time, loc_time2, min_sort2 ) 
        !   !$omp & reduction(min:time, min_sort, mincell)
            do gr = 1, ngran
                metod = 2
                if(gl_Gran_info(gr) == 2) CYCLE
                POTOK = 0.0
                konvect_(3, 1) = 0.0
                KOBL = 0
                kdir = 0
                idgod = 0
                s1 = gl_Gran_neighbour(1, gr)
                s2 = gl_Gran_neighbour(2, gr)
                qqq1 = gl_Cell_par(:, s1)
                konvect_(1, 1) = gl_Cell_par2(1, s1)
                null_bn = .False.
        
                distant = gl_Gran_center2(:, gr, now) - gl_Cell_center2(:, s1, now)
                dist = norm2(distant)
        
                !tvd1 = (gl_zone_Cell(s1) == 1) 
                tvd1 = (norm2(qqq1(2:4))/sqrt(ggg*qqq1(5)/qqq1(1)) > 2.2)
        
                ! Попробуем снести плотность пропорционально квадрату
                        if (.False.) then !if(norm2(qqq1(2:4))/sqrt(ggg*qqq1(5)/qqq1(1)) > 2.2) then
                        !if(gl_zone_Cell(s1) == 1) then
                            rad1 = norm2(gl_Cell_center2(:, s1, now))
                            rad2 = norm2(gl_Gran_center2(:, gr, now))
                            qqq1(1) = qqq1(1) * rad1**2 / rad2**2
                            qqq1(9) = qqq1(9) * rad1**2 / rad2**2
                            qqq1(5) = qqq1(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
                            ! Скорости сносим в сферической С.К.
                            call spherical_skorost(gl_Cell_center2(1, s1, now), gl_Cell_center2(2, s1, now), gl_Cell_center2(3, s1, now), &  
                                qqq1(2), qqq1(3), qqq1(4), aa, bb, cc)
                            call dekard_skorost(gl_Gran_center2(1, gr, now), gl_Gran_center2(2, gr, now), gl_Gran_center2(3, gr, now), &  
                                aa, bb, cc, qqq1(2), qqq1(3), qqq1(4))
    
                        end if
    
                        if (s2 >= 1) then
                            !if ( norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int .and. norm2(gl_Cell_center(:, s2)) <= par_R0 * par_R_int) CYCLE
                            qqq2 = gl_Cell_par(:, s2)
                            konvect_(2, 1) = gl_Cell_par2(1, s2)
                            !dist = min(gl_Cell_dist(s1), gl_Cell_dist(s2))   
                    
                            distant = gl_Gran_center2(:, gr, now) - gl_Cell_center2(:, s2, now)
                            dist = min( dist, norm2(distant))
    
                            ! Попробуем снести плотность пропорционально квадрату
                            if (.False.) then !(norm2(qqq2(2:4))/sqrt(ggg*qqq2(5)/qqq2(1)) > 2.2) then 
                            !if(gl_zone_Cell(s1) == 1) then
                                rad1 = norm2(gl_Cell_center2(:, s2, now))                              
                                rad2 = norm2(gl_Gran_center2(:, gr, now))
                                qqq2(1) = qqq2(1) * rad1**2 / rad2**2
                                qqq2(9) = qqq2(9) * rad1**2 / rad2**2
                                qqq2(5) = qqq2(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
                                call spherical_skorost(gl_Cell_center2(1, s2, now), gl_Cell_center2(2, s2, now), gl_Cell_center2(3, s2, now), &
                                    qqq2(2), qqq2(3), qqq2(4), aa, bb, cc)
                                call dekard_skorost(gl_Gran_center2(1, gr, now), gl_Gran_center2(2, gr, now), gl_Gran_center2(3, gr, now), &
                                    aa, bb, cc, qqq2(2), qqq2(3), qqq2(4))
    
                            end if
                    
    
                        else  ! В случае граничных ячеек - граничные условия
                            !if (norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int) CYCLE
                            if(s2 == -1) then  ! Набегающий поток
                                !dist = gl_Cell_dist(s1)
                        
                                qqq2 = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 1.0_8, -par_B_inf * cos(par_alphaB_inf), -par_B_inf * sin(par_alphaB_inf), 0.0_8, 100.0_8/)
                    
                                if(par_helium) then	
                                    konvect_(2, 1) = qqq2(1) * 0.3
            
                                    qqq2(1) = qqq2(1) * 1.3
                                    qqq2(9) = qqq2(9) * 1.3
                                    qqq2(5) = qqq2(5) * 1.075
                                end if
                    
                            else if(s2 == -3) then
                                !dist = gl_Cell_dist(s1)
                        
                                !if (gl_Cell_center2(2, s1, now) > 0.0_8) then
                !                   qqq2 = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 1.0_8, -par_B_inf * cos(par_alphaB_inf), -par_B_inf * sin(par_alphaB_inf), 0.0_8, 100.0_8/)
                                !else
                                !	qqq2 = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 1.0_8, qqq1(6), qqq1(7), qqq1(8), 100.0_8/)
                                !end if
                        
                                qqq2 = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 1.0_8, -par_B_inf * cos(par_alphaB_inf), -par_B_inf * sin(par_alphaB_inf), 0.0_8, 100.0_8/)
                        
                                !qqq2 = (/1.0_8, qqq1(2), qqq1(3), qqq1(4), 1.0_8, -par_B_inf * cos(par_alphaB_inf), -par_B_inf * sin(par_alphaB_inf), 0.0_8, 100.0_8/)
                                !if(qqq2(3) < 0.0) qqq2(3) = 0.0
                        
                                if(par_helium) then	
                                    konvect_(2, 1) = qqq2(1)  * 0.3
            
                                    qqq2(1) = qqq2(1) * 1.3
                                    qqq2(9) = qqq2(9) * 1.3
                                    qqq2(5) = qqq2(5) * 1.075
                                end if
                        
                            else  ! Здесь нужны мягкие условия
                                !dist = gl_Cell_dist(s1)
                                qqq2 = qqq1
                                konvect_(2, 1) = konvect_(1, 1)
                                !qqq2(5) = 1.0_8
                                if(qqq2(2) > -0.01) then
                                    qqq2(2) = 0.5 * par_Velosity_inf ! Отсос жидкости
                                end if
                        
                                !qqq2(5) = 0.00001   ! Маленькое противодавление
                        
                                !if(qqq2(6) > 0.0) then
                !                   qqq2(6) = -0.1 ! Отсос магнитного поля
                !               end if
    
                            end if
                end if
        
                !tvd2 = (gl_zone_Cell(s2) == 1) !
                tvd2 = (norm2(qqq2(2:4))/sqrt(ggg*qqq2(5)/qqq2(1)) > 2.2)
        
                ! Делаем ТВД
                if (s2 >= 1 .and. par_TVD == .True. .and. gl_Gran_type(gr) /= 2) then
                    !if(norm2(qqq1(2:4))/sqrt(ggg*qqq1(5)/qqq1(1)) < 2.2 .and. norm2(qqq2(2:4))/sqrt(ggg*qqq2(5)/qqq2(1)) < 2.2) then
                        ss1 = gl_Gran_neighbour_TVD(1, gr)
                        ss2 = gl_Gran_neighbour_TVD(2, gr)
                        if (ss1 /= 0 .and. ss2 /= 0) then
                            rast = gl_Gran_center2(:, gr, now) - gl_Cell_center2(:, s1, now)
                            df1 = norm2(rast)
                            rast = gl_Gran_center2(:, gr, now) - gl_Cell_center2(:, s2, now)
                            df2 = norm2(rast)
                            rast = gl_Gran_center2(:, gr, now) - gl_Cell_center2(:, ss1, now)
                            dff1 = norm2(rast)
                            rast = gl_Gran_center2(:, gr, now) - gl_Cell_center2(:, ss2, now)
                            dff2 = norm2(rast)
                            qqq11 = gl_Cell_par(:, ss1)
                            qqq22 = gl_Cell_par(:, ss2)
                            qqq11_2 = gl_Cell_par2(:, ss1)
                            qqq22_2 = gl_Cell_par2(:, ss2)
                    
                            !tvd3 = (gl_zone_Cell(ss1) == 1) 
                            tvd3 = (norm2(qqq11(2:4))/sqrt(ggg*qqq11(5)/qqq11(1)) > 2.2)
                            !tvd4 = (gl_zone_Cell(ss2) == 1) !
                            tvd4 = (norm2(qqq22(2:4))/sqrt(ggg*qqq22(5)/qqq22(1)) > 2.2)
                    
                            if(tvd1 == .True. .and. tvd2 == .False.) then
                                rad1 = norm2(gl_Cell_center2(:, s1, now))                              
                                rad5 = norm2(gl_Gran_center2(:, gr, now))
                        
                                qqq1_TVD = qqq1
                                qqq1_TVD_2 = konvect_(1, :)
                        
                                qqq1_TVD(1) = qqq1_TVD(1) * rad1**2 / rad5**2
                                qqq1_TVD_2(1) = qqq1_TVD_2(1) * rad1**2 / rad5**2
                                qqq1_TVD(9) = qqq1_TVD(9) * rad1**2 / rad5**2
                                qqq1_TVD(5) = qqq1_TVD(5) * rad1**(2 * ggg) / rad5**(2 * ggg)
                        
                                call spherical_skorost(gl_Cell_center2(1, s1, now), gl_Cell_center2(2, s1, now), gl_Cell_center2(3, s1, now), &
                                    qqq1_TVD(2), qqq1_TVD(3), qqq1_TVD(4), aa, bb, cc)
                        
                                call dekard_skorost(gl_Gran_center2(1, gr, now), gl_Gran_center2(2, gr, now), gl_Gran_center2(3, gr, now), &
                                    aa, bb, cc, qqq1_TVD(2), qqq1_TVD(3), qqq1_TVD(4))
                        
                        
                                do i = 1, 9
                                    !qqq2_TVD(i) = qqq2(i) ! linear(-dff2, qqq22(i), -df2, qqq2(i), df1, qqq1(i), 0.0_8)
                                    !qqq2_TVD(i) = linear2(-dff2, qqq22(i), -df2, qqq2(i), 0.0_8)
                                end do
                                
                                qqq2_TVD = qqq2
                                qqq2_TVD_2(1) = konvect_(2, 1)
                        
                                !qqq2_TVD_2(1) = konvect_(2, 1)! linear(-dff2, qqq22_2(1), -df2, konvect_(2, 1), df1, konvect_(1, 1), 0.0_8)	
                                !qqq2_TVD_2(1) = linear2(-dff2, qqq22_2(1), -df2, konvect_(2, 1), 0.0_8)	
                        
                            else if(tvd1 == .False. .and. tvd2 == .True.) then                              
                                rad2 = norm2(gl_Cell_center2(:, s2, now))                               
                                rad5 = norm2(gl_Gran_center2(:, gr, now))
                        
                                qqq2_TVD = qqq2
                                qqq2_TVD_2 = konvect_(2, :)
                        
                                qqq2_TVD(1) = qqq2_TVD(1) * rad2**2 / rad5**2
                                qqq2_TVD_2(1) = qqq2_TVD_2(1) * rad2**2 / rad5**2
                                qqq2_TVD(9) = qqq2_TVD(9) * rad2**2 / rad5**2
                                qqq2_TVD(5) = qqq2_TVD(5) * rad2**(2 * ggg) / rad5**(2 * ggg)
                        
                                call spherical_skorost(gl_Cell_center2(1, s2, now), gl_Cell_center2(2, s2, now), gl_Cell_center2(3, s2, now), &
                                    qqq2_TVD(2), qqq2_TVD(3), qqq2_TVD(4), aa, bb, cc)
                        
                                call dekard_skorost(gl_Gran_center2(1, gr, now), gl_Gran_center2(2, gr, now), gl_Gran_center2(3, gr, now), &
                                    aa, bb, cc, qqq2_TVD(2), qqq2_TVD(3), qqq2_TVD(4))
                        
                                do i = 1, 9
                                    qqq1_TVD(i) = linear(-dff1, qqq11(i), -df1, qqq1(i), df2, qqq2(i), 0.0_8)
                                end do
                        
                                qqq1_TVD_2(1) = linear(-dff1, qqq11_2(1), -df1, konvect_(1, 1), df2, konvect_(2, 1), 0.0_8)
                        
                            else if(tvd1 == .False. .and. tvd2 == .False.) then
                                do i = 1, 9
                                    qqq1_TVD(i) = linear(-dff1, qqq11(i), -df1, qqq1(i), df2, qqq2(i), 0.0_8)
                                    qqq2_TVD(i) = linear(-dff2, qqq22(i), -df2, qqq2(i), df1, qqq1(i), 0.0_8)
                                end do
                        
                                qqq1_TVD_2(1) = linear(-dff1, qqq11_2(1), -df1, konvect_(1, 1), df2, konvect_(2, 1), 0.0_8)
                                qqq2_TVD_2(1) = linear(-dff2, qqq22_2(1), -df2, konvect_(2, 1), df1, konvect_(1, 1), 0.0_8)	
                            else
                                rad1 = norm2(gl_Cell_center2(:, s1, now))                              
                                rad2 = norm2(gl_Cell_center2(:, s2, now))                              
                                rad3 = norm2(gl_Cell_center2(:, ss1, now))                              
                                rad4 = norm2(gl_Cell_center2(:, ss2, now))                              
                                rad5 = norm2(gl_Gran_center2(:, gr, now))
                        
                                qqq1_TVD(1) = linear(-dff1, qqq11(1) * rad3**2, -df1, qqq1(1) * rad1**2, df2, qqq2(1) * rad2**2, 0.0_8)/ rad5**2
                                qqq1_TVD_2(1) = linear(-dff1, qqq11_2(1) * rad3**2, -df1, konvect_(1, 1) * rad1**2, df2, konvect_(2, 1) * rad2**2, 0.0_8)/ rad5**2
                                qqq1_TVD(9) = linear(-dff1, qqq11(9) * rad3**2, -df1, qqq1(9) * rad1**2, df2, qqq2(9) * rad2**2, 0.0_8)/ rad5**2
                                qqq1_TVD(5) = linear(-dff1, qqq11(5) * rad3**(2 * ggg), -df1, qqq1(5) * rad1**(2 * ggg), df2,&
                                    qqq2(5) * rad2**(2 * ggg), 0.0_8)/ rad5**(2 * ggg)
                        
                                qqq2_TVD(1) = linear(-dff2, qqq22(1) * rad4**2, -df2, qqq2(1) * rad2**2, df1, qqq1(1) * rad1**2, 0.0_8)/ rad5**2
                                qqq2_TVD_2(1) = linear(-dff2, qqq22_2(1) * rad4**2, -df2, konvect_(2, 1) * rad2**2, df1, konvect_(1, 1) * rad1**2, 0.0_8)/ rad5**2
                                qqq2_TVD(9) = linear(-dff2, qqq22(9) * rad4**2, -df2, qqq2(9) * rad2**2, df1, qqq1(9) * rad1**2, 0.0_8)/ rad5**2
                                qqq2_TVD(5) = linear(-dff2, qqq22(5) * rad4**(2 * ggg), -df2, qqq2(5) * rad2**(2 * ggg), df1,&
                                    qqq1(5) * rad1**(2 * ggg), 0.0_8)/ rad5**(2 * ggg)
                        
                                ! Переводим скорости в сферическую с.к.
                                call spherical_skorost(gl_Cell_center2(3, s1, now), gl_Cell_center2(1, s1, now), gl_Cell_center2(2, s1, now), &
                                    qqq1(4), qqq1(2), qqq1(3), aa, bb, cc)
                                qqq1(4) = aa
                                qqq1(2) = bb
                                qqq1(3) = cc
                        
                                call spherical_skorost(gl_Cell_center2(3, s2, now), gl_Cell_center2(1, s2, now), gl_Cell_center2(2, s2, now), &
                                    qqq2(4), qqq2(2), qqq2(3), aa, bb, cc)
                                qqq2(4) = aa
                                qqq2(2) = bb
                                qqq2(3) = cc
                        
                                call spherical_skorost(gl_Cell_center2(3, ss1, now), gl_Cell_center2(1, ss1, now), gl_Cell_center2(2, ss1, now), &
                                    qqq11(4), qqq11(2), qqq11(3), aa, bb, cc)
                                qqq11(4) = aa
                                qqq11(2) = bb
                                qqq11(3) = cc
                        
                                call spherical_skorost(gl_Cell_center2(3, ss2, now), gl_Cell_center2(1, ss2, now), gl_Cell_center2(2, ss2, now), &
                                    qqq22(4), qqq22(2), qqq22(3), aa, bb, cc)
                                qqq22(4) = aa
                                qqq22(2) = bb
                                qqq22(3) = cc
                        
                                do i = 2, 4
                                    qqq1_TVD(i) = linear(-dff1, qqq11(i), -df1, qqq1(i), df2, qqq2(i), 0.0_8)
                                    qqq2_TVD(i) = linear(-dff2, qqq22(i), -df2, qqq2(i), df1, qqq1(i), 0.0_8)
                                end do
                        
                                do i = 6, 8
                                    qqq1_TVD(i) = linear(-dff1, qqq11(i), -df1, qqq1(i), df2, qqq2(i), 0.0_8)
                                    qqq2_TVD(i) = linear(-dff2, qqq22(i), -df2, qqq2(i), df1, qqq1(i), 0.0_8)
                                end do
                        
                                call dekard_skorost(gl_Gran_center2(3, gr, now), gl_Gran_center2(1, gr, now), gl_Gran_center2(2, gr, now), &
                                    qqq1_TVD(4), qqq1_TVD(2), qqq1_TVD(3), aa, bb, cc)
                                qqq1_TVD(4) = aa
                                qqq1_TVD(2) = bb
                                qqq1_TVD(3) = cc
                                call dekard_skorost(gl_Gran_center2(3, gr, now), gl_Gran_center2(1, gr, now), gl_Gran_center2(2, gr, now), &
                                    qqq2_TVD(4), qqq2_TVD(2), qqq2_TVD(3), aa, bb, cc)
                                qqq2_TVD(4) = aa
                                qqq2_TVD(2) = bb
                                qqq2_TVD(3) = cc
                            end if
                    
                            qqq1 = qqq1_TVD
                            qqq2 = qqq2_TVD
                    
                            konvect_(1, 1) = qqq1_TVD_2(1)
                            konvect_(2, 1) = qqq2_TVD_2(1)
                        !end if
                    end if
                end if
        
                ! На ударной волне не делаем снос справа
                if(gl_Gran_type(gr) == 1) then ! .and. gl_Gran_center2(1, gr, now) <= -5.0) then
                    !qqq1 = gl_Cell_par(:, s1)
                    qqq2 = gl_Cell_par(:, s2)
                end if
        
        
                ! Вычитаем для снесённых значений нормальною компоненту магнитного поля
                !if (gl_Gran_type(gr) == 2 .and. sqrt(gl_Gran_center2(2, gr, now)**2 + gl_Gran_center2(3, gr, now)**2) <= 15.0 .and. par_null_bn == .True.) then
                if (gl_Gran_type(gr) == 2 .and. gl_Gran_center2(1, gr, now) >= par_null_bn_x .and. par_null_bn == .True.) then
                    !write(*, *) " bn1 = ", DOT_PRODUCT(gl_Gran_normal2(:, gr, now), qqq1(6:8))
                    !write(*, *) " bn2 = ", DOT_PRODUCT(gl_Gran_normal2(:, gr, now), qqq2(6:8))
            
                    qqq1(6:8) = qqq1(6:8) - DOT_PRODUCT(gl_Gran_normal2(:, gr, now), qqq1(6:8)) * gl_Gran_normal2(:, gr, now)
                    qqq2(6:8) = qqq2(6:8) - DOT_PRODUCT(gl_Gran_normal2(:, gr, now), qqq2(6:8)) * gl_Gran_normal2(:, gr, now)
                end if
        
                
                ! Нужно вычислить скорость движения грани
                wc = DOT_PRODUCT((gl_Gran_center2(:, gr, now2) -  gl_Gran_center2(:, gr, now))/TT, gl_Gran_normal2(:, gr, now))
            
                metod = gl_Gran_scheme(gr)
                
                if(gl_Gran_type(gr) == 2 .or. gl_Gran_type(gr) == 1) metod = 3
                
                
                
                if (gl_Gran_type(gr) == 2) then ! Для гелиопаузы особая процедура
                    qqq1(5) = qqq1(5) + (norm2(qqq1(6:8))**2)/(8.0 * par_pi_8)
                    qqq2(5) = qqq2(5) + (norm2(qqq2(6:8))**2)/(8.0 * par_pi_8)
                    
                    call cgod3d(KOBL, 0, 0, 0, kdir, idgod, &
                                            gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), 1.0_8, &
                                            wc, qqq1(1:8), qqq2(1:8), &
                                            dsl, dsp, dsc, 1.0_8, 1.66666666666666_8, &
                                            POTOK)
                    
                    !print*, "igod", wc, dsc, gr
                    !PAUSE
                    
                    
                    POTOK(6:8) = 0.0
                    if (idgod == 2) then ! Если не удалось посчитать Годуновым
                        write(*, *) "Ne ydalos godunov 098767890987678923874224243234"
                        qqq1(5) = qqq1(5) - (norm2(qqq1(6:8))**2)/(8.0 * par_pi_8)
                        qqq2(5) = qqq2(5) - (norm2(qqq2(6:8))**2)/(8.0 * par_pi_8)
                        call chlld_Q(metod, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
                            wc, qqq1, qqq2, dsl, dsp, dsc, POTOK, null_bn, p_correct_ = .True., konvect_ = konvect_)
                    else
                        ! Нужно правильно посчитать конвективный снос
                        if(dsc >= wc) then
                            konvect_(3, 1) = konvect_(1, 1) * POTOK(1) / qqq1(1)
                        else
                            konvect_(3, 1) = konvect_(2, 1) * POTOK(1) / qqq2(1)
                        end if
                    end if
                else
                    if(qqq1(5) < 0.0 .or. qqq2(5) < 0.0) then
                        print*, "p < 0  678uygbnjytfgweyfgwyg3r23" 
                        STOP
                    end if
                    
                    call chlld_Q(metod, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
                    wc, qqq1, qqq2, dsl, dsp, dsc, POTOK, null_bn, p_correct_ = .True., konvect_ = konvect_)
                end if
                
                if(ieee_is_nan(POTOK(2)) .or. ieee_is_nan(POTOK(6))) then
                        print*,  "NUN u3  potok	"
                        print*, "___"
                        print*, POTOK
                        print*, "___"
                        print*, qqq1
                        print*, "___"
                        print*, qqq2
                        print*, "___"
                        print*, wc
                        print*, "___"
                        print*, gl_Gran_normal2(:, gr, now)
                        STOP
                end if
                
                loc_time = 0.99 * dist/( max(dabs(dsl), dabs(dsp)) + dabs(wc))
                
                if(loc_time < 0.000001) then
                    print*, dsl, dsp, dsc, wc, dist
                    print*, gl_Gran_normal2(:, gr, now)
                    print*, gl_Gran_center2(:, gr, now2)
                    print*, gl_Gran_center2(:, gr, now)
                    STOP
                end if
                
                
                gl_Gran_POTOK(1:9, gr) = POTOK * gl_Gran_square2(gr, now)
                gl_Gran_POTOK2(1, gr) = konvect_(3, 1) * gl_Gran_square2(gr, now)
                
                
                gl_Gran_POTOK(10, gr) = 0.5 * DOT_PRODUCT(gl_Gran_normal2(:, gr, now), qqq1(6:8) + qqq2(6:8)) * gl_Gran_square2(gr, now)

                
                ! Теперь сравниваем с глобальным временем
                if (gl_Cell_center2(1, s1, now) > -100000.0) then
                    if(loc_time < time) then
                        !$omp critical
                        if(loc_time < time) then
                            time = loc_time
                            ! time = min(time, loc_time )   ! REDUCTION
                            mincell = gr
                            min_sort = min_sort2
                        end if
                        !$omp end critical
                    end if
                end if
                


            end do
            !$omp end do

            
            !$omp barrier ! Синхронизация нитей
            

            !!$omp single
            !    if (mod(st, 5) == 0)  print*, st, time
            !!$omp end single

            ! Теперь цикл по ячейкам
            !$omp do private(n_He, POTOK2, POTOK, Volume, Volume2, qqq, i, j, ro3, u3, v3, w3, p3, Q3, bx3, by3, bz3, zone, l_1, fluid1, SOURSE, sks, ijk)
            do gr = 1, ncell
                if(gl_Cell_info(gr) == 0) CYCLE
                l_1 = .TRUE.
                if (norm2(gl_Cell_center2(:, gr, now)) <= par_R0 + (par_R_character - par_R0) * (3.0_8/par_n_TS)**par_kk1) l_1 = .FALSE. 
                !if ((gl_Cell_type(gr) == "A" .or. gl_Cell_type(gr) == "B").and.(gl_Cell_number(1, gr) <= 2) ) l_1 = .FALSE.    ! Не считаем внутри сферы
                POTOK = 0.0
                POTOK2 = 0.0
                sks = 0.0
                SOURSE = 0.0
                Volume = gl_Cell_Volume2(gr, now)
                Volume2 = gl_Cell_Volume2(gr, now2)
                
                
        
                !Volume2 = Volume
        
                qqq = gl_Cell_par(:, gr)
                n_He = gl_Cell_par2(1, gr)
                
                fluid1 = gl_Cell_par_MK(1:5, :, gr)
                MK_kk = gl_Cell_par_MK(6:10, 1, gr)
                !if(gl_Cell_center2(1, gr, now) < -150.0) then
                !	MK_kk = 1.0
                !end if
        
                !MK_kk = 1.0
                ! Просуммируем потоки через грани
                do i = 1, 6
                    j = gl_Cell_gran(i, gr)
                    if (j == 0) CYCLE
                    if (j < 0) write(*, *) "ERROR 3876tfghjuyghejk"
                    if (gl_Gran_neighbour(1, j) == gr) then
                        POTOK = POTOK + gl_Gran_POTOK(1:9, j)
                        POTOK2 = POTOK2 + gl_Gran_POTOK2(1, j)
                        sks = sks + gl_Gran_POTOK(10, j)
                    else
                        POTOK = POTOK - gl_Gran_POTOK(1:9, j)
                        POTOK2 = POTOK2 - gl_Gran_POTOK2(1, j)
                        sks = sks - gl_Gran_POTOK(10, j)
                    end if
                end do
    
    
                ! Определяем зону в которой находимся
                if(qqq(9)/qqq(1) < 50.0) then
                    !if(norm2(qqq(2:4))/sqrt(ggg*qqq(5)/qqq(1)) > 1.3) then
                    if(gl_zone_Cell(gr) == 1) then
                        zone = 1
                    else
                        zone = 2
                    end if
                else
                    if(norm2(qqq(2:4))/sqrt(ggg*qqq(5)/qqq(1)) > 1.1) then
                        zone = 4
                    else
                        zone = 3
                    end if
                end if
        
                if(zone == 1 .or. zone == 2) then ! Перенормировка
                    qqq(2:4) = qqq(2:4) * (par_chi/par_chi_real)
                    qqq(1) = qqq(1) / (par_chi/par_chi_real)**2
                end if
        
                ! He
                qqq(1) = qqq(1) - n_He
                
                if(zone == 1 .or. zone == 2) then
                    qqq(5) = qqq(5) / (8.0 * qqq(1) + 3.0 * n_He) * 8.0 * qqq(1)
                else
                    qqq(5) = qqq(5) / (4.0 * qqq(1) + n_He) * 4.0 * qqq(1)
                end if
                
                call Calc_sourse_MF(qqq, fluid1, SOURSE, zone)  ! Вычисляем источники
                
                if(zone == 1 .or. zone == 2) then
                    qqq(5) = qqq(5) * (8.0 * qqq(1) + 3.0 * n_He) / (8.0 * qqq(1))
                else
                    qqq(5) = qqq(5) * (4.0 * qqq(1) + n_He) / (4.0 * qqq(1))
                end if
        
                qqq(1) = qqq(1) + n_He
                
                if(zone == 1 .or. zone == 2) then ! Перенормировка обратно
                    qqq(2:4) = qqq(2:4) / (par_chi/par_chi_real)
                    qqq(1) = qqq(1) * (par_chi/par_chi_real)**2
                    SOURSE(5, 1) = SOURSE(5, 1)/ (par_chi/par_chi_real)
                end if
            
                
                ! Домножаем на коэффицент 
                do i = 2, 5
                    SOURSE(i, 1) = SOURSE(i, 1) * MK_kk(i)
                end do
        
    
                if (l_1 == .TRUE.) then
                    ro3 = qqq(1)* Volume / Volume2 - time * POTOK(1) / Volume2 + time * MK_kk(1)
                    gl_Cell_par2(1, gr) = n_He * Volume / Volume2 - time * POTOK2 / Volume2
                    Q3 = qqq(9)* Volume / Volume2 - time * POTOK(9) / Volume2 + (qqq(9)/qqq(1)) *  time * MK_kk(1)
                    if (ro3 <= 0.0_8) then
                        write(*, *) "Ro < 0  3688"
                        write(*, *) " ---  ", ro3, qqq(1), gl_Cell_center2(1, gr, now), gl_Cell_center2(2, gr, now), gl_Cell_center2(3, gr, now)! , MK_kk(1), &
                            ! qqq(1), time * POTOK(1) / Volume2
                        !write(*, *) qqq(1), Q3
                        !write(*, *) Volume , Volume2
                        ro3 = qqq(1) * 0.8
                        Q3 = (qqq(9)/qqq(1)) * qqq(1) * 0.8
                    end if
                    
                    if(gl_Cell_par2(1, gr) <= 0.0) then
                        gl_Cell_par2(1, gr) = 0.000001
                        !write(*, *) "gl_Cell_par2(1, gr) < 0  789uhgtyefef"
                        !write(*, *) gl_Cell_par2(1, gr), qqq2, POTOK2, gl_Cell_center2(1, gr, now), gl_Cell_center2(2, gr, now), gl_Cell_center2(3, gr, now)
                        !stop
                    end if
                    
                    
                    
                    
                    u3 = (qqq(1) * qqq(2)* Volume / Volume2 - time * (POTOK(2) + (qqq(6)/cpi4) * sks) / Volume2 + time * SOURSE(2, 1)) / ro3
                    v3 = (qqq(1) * qqq(3)* Volume / Volume2 - time * (POTOK(3) + (qqq(7)/cpi4) * sks) / Volume2 + time * SOURSE(3, 1)) / ro3
                    w3 = (qqq(1) * qqq(4)* Volume / Volume2 - time * (POTOK(4) + (qqq(8)/cpi4) * sks) / Volume2 + time * SOURSE(4, 1)) / ro3
                    
                    bx3 = qqq(6) * Volume / Volume2 - time * (POTOK(6) + qqq(2) * sks) / Volume2
                    by3 = qqq(7) * Volume / Volume2 - time * (POTOK(7) + qqq(3) * sks) / Volume2
                    bz3 = qqq(8) * Volume / Volume2 - time * (POTOK(8) + qqq(4) * sks) / Volume2
                    
                    p3 = ((  ( qqq(5) / (ggg - 1.0) + 0.5 * qqq(1) * norm2(qqq(2:4))**2 + (qqq(6)**2 + qqq(7)**2 + qqq(8)**2) / 25.13274122871834590768 )* Volume / Volume2 &
                        - time * ( POTOK(5) + (DOT_PRODUCT(qqq(2:4), qqq(6:8))/cpi4) * sks)/ Volume2 + time * SOURSE(5, 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) - (bx3**2 + by3**2 + bz3**2) / 25.13274122871834590768 ) * (ggg - 1.0)
                    
                    
                    if (p3 <= 0.0_8) then
                        !print*, "p < 0  plasma 2028 ", p3 , gl_Cell_center(:, gr)
                        !write(*, *) "p < 0  3688"
                        !write(*, *) "p ----", p3, gl_Cell_center2(1, gr, now), gl_Cell_center2(2, gr, now), gl_Cell_center2(3, gr, now)
                        p3 = 0.000001
                        !pause
                    end if
    
                    gl_Cell_par(:, gr) = (/ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3/)

                end if

                
            end do
            !$omp end do
            
            !$omp barrier ! Синхронизация нитей
            
            !$omp end parallel
            
            
            gl_Vx = 0.0
            gl_Vy = 0.0
            gl_Vz = 0.0
            gl_Point_num = 0
            
            
            gl_x = gl_x2(:, now2)
            gl_y = gl_y2(:, now2)
            gl_z = gl_z2(:, now2)
            gl_Cell_Volume = gl_Cell_Volume2(:, now2)
            gl_Gran_normal = gl_Gran_normal2(:, :, now2)
            gl_Gran_center = gl_Gran_center2(:, :, now2)
            gl_Cell_center = gl_Cell_center2(:, :, now2)
            gl_Gran_square = gl_Gran_square2(:, now2)
            
            if (mod(step, 5000) == 0) then
                call Get_MK_to_MHD()
                call Find_TVD_sosed()
            end if
            
            if (.False. .and. mod(step, 2000) == 0) then
                call Start_MGD_3_inner_MK(2000)
            end if
            
            if (mod(step, 500) == 0 .or. step == 1000 .or. step == 3000) then
                print*, "PECHAT"
                call Print_surface_2D()
                call Print_Setka_2D()
                call Print_par_2D()
                call Print_Contact_3D()
                call Print_TS_3D()
            end if
            
            
            
            
        end do


        gl_x = gl_x2(:, 2)
        gl_y = gl_y2(:, 2)
        gl_z = gl_z2(:, 2)
        gl_Cell_Volume = gl_Cell_Volume2(:, 2)
        gl_Gran_normal = gl_Gran_normal2(:, :, 2)
        gl_Gran_center = gl_Gran_center2(:, :, 2)
        gl_Cell_center = gl_Cell_center2(:, :, 2)
        gl_Gran_square = gl_Gran_square2(:, 2)

        end_time = omp_get_wtime()
        print *, "Time work: ", (end_time-start_time)/60.0, "   in minutes"
		


	end subroutine Start_MGD_move_MK
	
	subroutine Start_MGD_move()
        ! Функция газовой динамики + мультифлюид + подвижная сетка
        ! Это управляющая функция, в отличие от предыдущих версий, здесь показано как правильно вызывать процедуры для движения сетки
        use GEO_PARAM
        use STORAGE
        USE ieee_arithmetic 
        use My_func
        USE OMP_LIB
        use Solvers
        implicit none
        integer :: step, now, now2, step2, min_sort, ijk, ijk2, min_sort2, N2, N3
        integer(4) :: st, gr, ngran, ncell, s1, s2, i, j, k, zone, iter, metod, mincell, ss1, ss2
        real(8) :: qqq1(9), qqq2(9), qqq(9)  ! Переменные в ячейке
        real(8) :: fluid1(5, par_n_sort), fluid2(5, par_n_sort)
        real(8) :: dist, dsl, dsc, dsp, start_time, end_time, rast(3), df1, df2, dff1, dff2, qqq11(9), qqq22(9)
        real(8) :: POTOK(9), qqq1_TVD(9), qqq2_TVD(9)
        real(8) :: POTOK_MF(5)
        real(8) :: POTOK_MF_all(5, 4)
        real(8) :: time, Volume, Volume2, TT, U8, rad1, rad2, aa, bb, cc, wc, sks, loc_time, loc_time2, rr, y, z
        real(8) :: SOURSE(5,par_n_sort + 1)  ! Источники массы, импульса и энергии для плазмы и каждого сорта мультифлюида

        real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3

        logical :: l_1, null_bn

        min_sort = 0

        ! Сначала подготовим все массивы для правильного движения
        if (allocated(gl_Vx) == .False.) then   ! Если движение запускается впервый раз, то нужно выделить память под используемые массивы и задать значения
            allocate(gl_Vx(size(gl_x)))
            allocate(gl_Vy(size(gl_x)))
            allocate(gl_Vz(size(gl_x)))
            allocate(gl_Point_num(size(gl_x)))
            allocate(gl_x2(size(gl_x), 2))
            allocate(gl_y2(size(gl_x), 2))
            allocate(gl_z2(size(gl_x), 2))
            allocate(gl_Cell_Volume2(size(gl_Cell_Volume), 2))
            allocate(gl_Gran_normal2(  size(gl_Gran_normal(:, 1)), size(gl_Gran_normal(1, :)), 2  ))
            allocate(gl_Gran_center2(  size(gl_Gran_center(:, 1)), size(gl_Gran_center(1, :)), 2  ))
            allocate(gl_Cell_center2(size(gl_Cell_center(:, 1)), size(gl_Cell_center(1, :)), 2))
            allocate(gl_Gran_square2(size(gl_Gran_square), 2))
            
            ! Начальная инициализация
            gl_Vx = 0.0
            gl_Vy = 0.0
            gl_Vz = 0.0
            gl_Point_num = 0
            gl_x2(:, 1) = gl_x
            gl_x2(:, 2) = gl_x
            gl_y2(:, 1) = gl_y
            gl_y2(:, 2) = gl_y
            gl_z2(:, 1) = gl_z
            gl_z2(:, 2) = gl_z
            gl_Cell_Volume2(:, 1) = gl_Cell_Volume
            gl_Cell_Volume2(:, 2) = gl_Cell_Volume
            gl_Gran_normal2(:, :, 1) = gl_Gran_normal
            gl_Gran_normal2(:, :, 2) = gl_Gran_normal
            gl_Gran_center2(:, :, 1) = gl_Gran_center
            gl_Gran_center2(:, :, 2) = gl_Gran_center
            gl_Cell_center2(:, :, 1) = gl_Cell_center
            gl_Cell_center2(:, :, 2) = gl_Cell_center
            gl_Gran_square2(:, 1) = gl_Gran_square
            gl_Gran_square2(:, 2) = gl_Gran_square
        end if
        
        
        
        start_time = omp_get_wtime()
        
        mincell = 1
        
        !call Start_MGD_3_inner(10000)
        
        
        ! Запускаем глобальный цикл
        now = 2                           ! Какие параметры сейчас будут считаться (1 или 2). Они меняются по очереди
        time = 0.00002_8               ! Начальная инициализация шага по времени 
        do step = 1, 1000  !    ! Нужно чтобы это число было чётным!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            if (mod(step, 10) == 0) then
                print*, "Step = ", step , "  step_time = ", time, "  mingran = ", mincell, & 
                "  gran Center = ", gl_Gran_center2(:, mincell, now), "  minsort = ", min_sort, "  cell = ", gl_Gran_neighbour(1, mincell), gl_Gran_neighbour(2, mincell) 
                !print*, gl_Gran_normal2(:, mincell, now)
                print*, "   "
                !print*, "____1____"
                !print*,  gl_Cell_par(:, gl_Gran_neighbour(1, mincell))
                !print*, "____2____"
                !print*, gl_Cell_par(:, gl_Gran_neighbour(2, mincell))
                !print*, "____3____"
                !print*, gl_Cell_dist(gl_Gran_neighbour(1, mincell)), gl_Cell_dist(gl_Gran_neighbour(2, mincell))
                !print*, "____4____"
                !print*, gl_Cell_Volume2(gl_Gran_neighbour(1, mincell), now), gl_Cell_Volume2(gl_Gran_neighbour(2, mincell), now)
                !print*, "____5____"
                !print*, gl_Cell_Volume2(gl_Gran_neighbour(1, mincell), now2), gl_Cell_Volume2(gl_Gran_neighbour(2, mincell), now2)
                !
                !
                !N3 = size(gl_RAY_O(1, 1, :))
                !N2 = size(gl_RAY_O(1, :, 1))
                !rr = 10000.0
                !do k = 1, N3
                !	do j = 1, N2
                !        y = gl_y2(gl_RAY_O(1, j , k), now2)
                !        z = gl_z2(gl_RAY_O(1, j, k), now2)
                !        print*, sqrt(y**2 + z**2)
                !	end do
                !end do
                
            end if
                
            now2 = now
            now = mod(now, 2) + 1
            
            TT = time
            time = 100000.0
            !TT = 0.0
            ! Вычисляем новые скорости управляющих узлов
            
            call Calc_move(now)   ! Записали скорости каждого узла в этот узел для последующего движения
            
            
            ! Двигаем все узлы сетки в соответствии с расчитанными скоростями в предыдущей функции
            call Move_all(now, TT) 
            
            call calc_all_Gran_move(now2)   ! Расчитываются новые объёмы\площади\нормали и т.д.
            !CYCLE
            
            ! Теперь по основным граням
            ngran = size(gl_all_Gran(1, :))
            ncell = size(gl_all_Cell(1, :))
            
            !$omp parallel
            
        !$omp do private(metod, POTOK, ss1, ss2, qqq1_TVD, qqq2_TVD, rast, df1, df2, dff1, dff2, qqq11, qqq22, s1, s2, qqq1, qqq2, dist, dsl, dsp, dsc, rad1, rad2, aa, bb, cc, fluid1, fluid2, POTOK_MF, wc, null_bn, loc_time, loc_time2, min_sort2 ) 
        !   !$omp & reduction(min:time, min_sort, mincell)
            do gr = 1, ngran
                metod = 2
                if(gl_Gran_info(gr) == 2) CYCLE
                POTOK = 0.0_8
                loc_time = 10000000000.0
                s1 = gl_Gran_neighbour(1, gr)
                s2 = gl_Gran_neighbour(2, gr)
                qqq1 = gl_Cell_par(:, s1)
                fluid1 = gl_Cell_par_MF(:, :, s1)   ! Загрузили параметры жидкостей для мультифлюида
                null_bn = .False.
                
                metod = gl_Gran_scheme(gr)
                
                dist = norm2(gl_Gran_center2(:, gr, now) - gl_Cell_center2(:, s1, now))

                ! Попробуем снести плотность пропорционально квадрату
                if(norm2(qqq1(2:4))/sqrt(ggg*qqq1(5)/qqq1(1)) > 2.2) then
                    !metod = 1
                    rad1 = norm2(gl_Cell_center2(:, s1, now))
                    rad2 = norm2(gl_Gran_center2(:, gr, now))
                    qqq1(1) = qqq1(1) * rad1**2 / rad2**2
                    qqq1(9) = qqq1(9) * rad1**2 / rad2**2
                    qqq1(5) = qqq1(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
                    ! Скорости сносим в сферической С.К.
                    call spherical_skorost(gl_Cell_center2(1, s1, now), gl_Cell_center2(2, s1, now), gl_Cell_center2(3, s1, now), &  
                        qqq1(2), qqq1(3), qqq1(4), aa, bb, cc)
                    call dekard_skorost(gl_Gran_center2(1, gr, now), gl_Gran_center2(2, gr, now), gl_Gran_center2(3, gr, now), &  
                        aa, bb, cc, qqq1(2), qqq1(3), qqq1(4))

                    call spherical_skorost(gl_Cell_center2(1, s1, now), gl_Cell_center2(2, s1, now), gl_Cell_center2(3, s1, now), & 
                        fluid1(2, 1), fluid1(3, 1), fluid1(4, 1), aa, bb, cc)
                    call dekard_skorost(gl_Gran_center2(1, gr, now), gl_Gran_center2(2, gr, now), gl_Gran_center2(3, gr, now), & 
                        aa, bb, cc, fluid1(2, 1), fluid1(3, 1), fluid1(4, 1))
                end if

                if (s2 >= 1) then
                    !if ( norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int .and. norm2(gl_Cell_center(:, s2)) <= par_R0 * par_R_int) CYCLE
                    qqq2 = gl_Cell_par(:, s2)
                    fluid2 = gl_Cell_par_MF(:, :, s2)   ! Загрузили параметры жидкостей для мультифлюида
                    !dist = min(gl_Cell_dist(s1), gl_Cell_dist(s2))   
                    
                    dist = min( dist, norm2(gl_Gran_center2(:, gr, now) - gl_Cell_center2(:, s2, now)))

                    ! Попробуем снести плотность пропорционально квадрату
                    if(norm2(qqq2(2:4))/sqrt(ggg*qqq2(5)/qqq2(1)) > 2.2) then 
                        !metod = 1
                        rad1 = norm2(gl_Cell_center2(:, s2, now))                              
                        rad2 = norm2(gl_Gran_center2(:, gr, now))
                        qqq2(1) = qqq2(1) * rad1**2 / rad2**2
                        qqq2(9) = qqq2(9) * rad1**2 / rad2**2
                        qqq2(5) = qqq2(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
                        call spherical_skorost(gl_Cell_center2(1, s2, now), gl_Cell_center2(2, s2, now), gl_Cell_center2(3, s2, now), &
                            qqq2(2), qqq2(3), qqq2(4), aa, bb, cc)
                        call dekard_skorost(gl_Gran_center2(1, gr, now), gl_Gran_center2(2, gr, now), gl_Gran_center2(3, gr, now), &
                            aa, bb, cc, qqq2(2), qqq2(3), qqq2(4))

                        call spherical_skorost(gl_Cell_center2(1, s2, now), gl_Cell_center2(2, s2, now), gl_Cell_center2(3, s2, now), &
                            fluid2(2, 1), fluid2(3, 1), fluid2(4, 1), aa, bb, cc)
                        call dekard_skorost(gl_Gran_center2(1, gr, now), gl_Gran_center2(2, gr, now), gl_Gran_center2(3, gr, now), &
                            aa, bb, cc, fluid2(2, 1), fluid2(3, 1), fluid2(4, 1))
                    end if

                else  ! В случае граничных ячеек - граничные условия
                    !if (norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int) CYCLE
                    if(s2 == -1) then  ! Набегающий поток
                        !dist = gl_Cell_dist(s1)
                        
                        qqq2 = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 1.0_8, -par_B_inf * cos(par_alphaB_inf), -par_B_inf * sin(par_alphaB_inf), 0.0_8, 100.0_8/)
                        fluid2(:, 1) = fluid1(:, 1)
                        fluid2(:, 2) = fluid1(:, 2)
                        fluid2(:, 3) = fluid1(:, 3)
                        fluid2(:, 4) = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 0.5_8/)
                    else if(s2 == -3) then
                        !dist = gl_Cell_dist(s1)
                        
                        if (gl_Cell_center2(2, s1, now) > 0.0_8) then
                            qqq2 = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 1.0_8, -par_B_inf * cos(par_alphaB_inf), -par_B_inf * sin(par_alphaB_inf), 0.0_8, 100.0_8/)
                        else
                            qqq2 = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 1.0_8, qqq1(6), qqq1(7), qqq1(8), 100.0_8/)
                        end if
                        
                        !qqq2 = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 1.0_8, -par_B_inf * cos(par_alphaB_inf), -par_B_inf * sin(par_alphaB_inf), 0.0_8, 100.0_8/)
                        qqq2 = (/1.0_8, qqq1(2), qqq1(3), qqq1(4), 1.0_8, -par_B_inf * cos(par_alphaB_inf), -par_B_inf * sin(par_alphaB_inf), 0.0_8, 100.0_8/)
                        if(qqq2(3) < 0.0) qqq2(3) = 0.0
                        
                        fluid2(:, 1) = fluid1(:, 1)
                        fluid2(:, 2) = fluid1(:, 2)
                        fluid2(:, 3) = fluid1(:, 3)
                        fluid2(:, 4) = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 0.5_8/)
                    else  ! Здесь нужны мягкие условия
                        !dist = gl_Cell_dist(s1)
                        qqq2 = qqq1
                        fluid2 = fluid1
                        !qqq2(5) = 1.0_8
                        if(qqq2(2) > 0.2 * par_Velosity_inf) then
                            qqq2(2) = 0.2 * par_Velosity_inf ! Отсос жидкости
                        end if

                    end if
                end if
                
                if (s2 >= 1 .and. par_TVD == .True. .and. gl_Gran_type(gr) /= 2) then
                    if(norm2(qqq1(2:4))/sqrt(ggg*qqq1(5)/qqq1(1)) < 2.2 .and. norm2(qqq2(2:4))/sqrt(ggg*qqq2(5)/qqq2(1)) < 2.2) then
                        ss1 = gl_Gran_neighbour_TVD(1, gr)
                        ss2 = gl_Gran_neighbour_TVD(2, gr)
                        if (ss1 /= 0 .and. ss2 /= 0) then
                            rast = gl_Gran_center2(:, gr, now) - gl_Cell_center2(:, s1, now)
                            df1 = norm2(rast)
                            rast = gl_Gran_center2(:, gr, now) - gl_Cell_center2(:, s2, now)
                            df2 = norm2(rast)
                            rast = gl_Gran_center2(:, gr, now) - gl_Cell_center2(:, ss1, now)
                            dff1 = norm2(rast)
                            rast = gl_Gran_center2(:, gr, now) - gl_Cell_center2(:, ss1, now)
                            dff2 = norm2(rast)
                            qqq11 = gl_Cell_par(:, ss1)
                            qqq22 = gl_Cell_par(:, ss2)
                    
                            do i = 1, 9
                                qqq1_TVD(i) = linear(-dff1, qqq11(i), -df1, qqq1(i), df2, qqq2(i), 0.0_8)
                                qqq2_TVD(i) = linear(-dff2, qqq22(i), -df2, qqq2(i), df1, qqq1(i), 0.0_8)
                            end do
                    
                            qqq1 = qqq1_TVD
                            qqq2 = qqq2_TVD
                        end if
                    end if
                end if
                
                ! Вычитаем для снесённых значений нормальною компоненту магнитного поля
                !if (gl_Gran_type(gr) == 2 .and. sqrt(gl_Gran_center2(2, gr, now)**2 + gl_Gran_center2(3, gr, now)**2) <= 15.0 .and. par_null_bn == .True.) then
                if (gl_Gran_type(gr) == 2 .and. gl_Gran_center2(1, gr, now) >= par_null_bn_x .and. par_null_bn == .True.) then
                    qqq1(6:8) = qqq1(6:8) - DOT_PRODUCT(gl_Gran_normal2(:, gr, now), qqq1(6:8)) * gl_Gran_normal2(:, gr, now)
                    qqq2(6:8) = qqq2(6:8) - DOT_PRODUCT(gl_Gran_normal2(:, gr, now), qqq2(6:8)) * gl_Gran_normal2(:, gr, now)

        end if
                
                ! Нужно вычислить скорость движения грани
                wc = DOT_PRODUCT((gl_Gran_center2(:, gr, now2) -  gl_Gran_center2(:, gr, now))/TT, gl_Gran_normal2(:, gr, now))
            
                
                
                metod = gl_Gran_scheme(gr)
                
                if(gl_Gran_type(gr) == 2 .or. gl_Gran_type(gr) == 1) metod = 3
                
                !if(gl_Gran_type(gr) == 2) null_bn = .True.
                
                !if (gr == 236608) then
                !	write(*, *) qqq1(1), qqq1(2), qqq1(3), qqq1(4), qqq1(5), qqq1(6), qqq1(7), qqq1(8), qqq1(9)
                !	write(*, *) "___"
                !	write(*, *) qqq2(1), qqq2(2), qqq2(3), qqq2(4), qqq2(5), qqq2(6), qqq2(7), qqq2(8), qqq2(9)
                !	write(*, *) "___"
                !	write(*, *) metod, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now)
                !	write(*, *) "___"
                !	write(*, *) wc
                !	write(*, *) "___"
                !	write(*, *) null_bn
                !	write(*, *) "___"
                !	continue
                !end if
                
                
                if (.True.) then !(gl_Gran_type(gr) == 1) then
                    call chlld_Q(metod, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
                    wc, qqq1, qqq2, dsl, dsp, dsc, POTOK, null_bn)
                else
                    call chlld_Q(metod, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
                    wc, qqq1, qqq2, dsl, dsp, dsc, POTOK, null_bn, 0)
                end if
                
                if(ieee_is_nan(POTOK(2)) .or. ieee_is_nan(POTOK(6))) then
                        print*,  "NUN u3  potok	"
                        print*, "___"
                        print*, POTOK
                        print*, "___"
                        print*, qqq1
                        print*, "___"
                        print*, qqq2
                        print*, "___"
                        print*, wc
                        print*, "___"
                        print*, gl_Gran_normal2(:, gr, now)
                        STOP
                    end if
                
                if (.False.) then
                !if (gr == 77) then
                    write(*,*) "__***__ 614302"
                    write(*,*) gl_Gran_info(gr), s1, s2
                    write(*,*) "__A__"
                    write(*,*)  POTOK(1), POTOK(2), POTOK(3), POTOK(4), POTOK(5), POTOK(6), POTOK(7), POTOK(8), POTOK(9)
                    write(*,*) "__B__"
                    write(*,*) gl_Gran_square2(gr, now), gl_Gran_square2(gr, now2)
                    write(*,*) "__C__"
                    write(*,*) dist, wc, dsl, dsp, dsc, TT
                    write(*,*) "__D__"
                    write(*,*) gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now)
                    write(*,*) "__D2__"
                    write(*,*) gl_Gran_normal2(1, gr, now2), gl_Gran_normal2(2, gr, now2), gl_Gran_normal2(3, gr, now2)
                    write(*,*) "__E__"
                    write(*,*) qqq1(1), qqq1(2), qqq1(3), qqq1(4), qqq1(5), qqq1(6), qqq1(7), qqq1(8), qqq1(9)
                    write(*,*) "__F__"
                    write(*,*) qqq2(1), qqq2(2), qqq2(3), qqq2(4), qqq2(5), qqq2(6), qqq2(7), qqq2(8), qqq2(9)
                    write(*,*) "__G__"
                    write(*,*) gl_Cell_Volume2(s1, now), gl_Cell_Volume2(s1, now2)
                    write(*,*) "__H__"
                    write(*,*) gl_Cell_Volume2(s2, now), gl_Cell_Volume2(s2, now2)
                    write(*,*) "************"
                    write(*,*) "************"
                end if
                
                !call chlld_Q(metod, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
                !    wc, qqq1, qqq2, dsl, dsp, dsc, POTOK, null_bn)
                loc_time2 = 0.99 * dist/( max(dabs(dsl), dabs(dsp)) + dabs(wc))
                
                if(loc_time2 < loc_time) then
                    loc_time = loc_time2
                    !loc_time = min(loc_time2, loc_time )   ! REDUCTION
                    min_sort2 = 0
                end if
                
                gl_Gran_POTOK(1:9, gr) = POTOK * gl_Gran_square2(gr, now)
                
                
                
                gl_Gran_POTOK(10, gr) = 0.5 * DOT_PRODUCT(gl_Gran_normal2(:, gr, now), qqq1(6:8) + qqq2(6:8)) * gl_Gran_square2(gr, now)

                
                
                call chlld_gd(0, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
                    wc, fluid1(:, 1), fluid2(:, 1), dsl, dsp, dsc, POTOK_MF)
                
                loc_time2 = 0.99 * dist/( max(dabs(dsl), dabs(dsp)) + dabs(wc))
                
                if(loc_time2 < loc_time) then
                    loc_time = loc_time2
                    !loc_time = min(loc_time2, loc_time )   ! REDUCTION
                    min_sort2 = 1
                end if
                
                
                gl_Gran_POTOK_MF(:, 1, gr) = POTOK_MF * gl_Gran_square2(gr, now)

                call chlld_gd(1, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
                    wc, fluid1(:, 2), fluid2(:, 2), dsl, dsp, dsc, POTOK_MF)
                
                
                loc_time2 = 0.99 * dist/( max(dabs(dsl), dabs(dsp)) + dabs(wc))
                
                if(loc_time2 < loc_time) then
                    loc_time = loc_time2
                    !loc_time = min(loc_time2, loc_time )   ! REDUCTION
                    min_sort2 = 2
                end if
                
                gl_Gran_POTOK_MF(:, 2, gr) = POTOK_MF * gl_Gran_square2(gr, now)
                
                call chlld_gd(1, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
                    wc, fluid1(:, 3), fluid2(:, 3), dsl, dsp, dsc, POTOK_MF)
                
                loc_time2 = 0.99 * dist/( max(dabs(dsl), dabs(dsp)) + dabs(wc))
                
                if(loc_time2 < loc_time) then
                    loc_time = loc_time2
                    !loc_time = min(loc_time2, loc_time )   ! REDUCTION
                    min_sort2 = 3
                end if
                
                
                gl_Gran_POTOK_MF(:, 3, gr) = POTOK_MF * gl_Gran_square2(gr, now)

                call chlld_gd(1, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
                    wc, fluid1(:, 4), fluid2(:, 4), dsl, dsp, dsc, POTOK_MF)
                
                loc_time2 = 0.99 * dist/( max(dabs(dsl), dabs(dsp)) + dabs(wc))
                
                if(loc_time2 < loc_time) then
                    loc_time = loc_time2
                    !loc_time = min(loc_time2, loc_time )   ! REDUCTION
                    min_sort2 = 4
                end if
                
                
                if( .False.) then
                !if (gr == 562093 .and. mod(step, 100) == 0 ) then
                    write(*,*) "__***__ 562093"
                    write(*,*) gl_Gran_info(gr), s1, s2
                    write(*,*) "__A__"
                    write(*,*)  gl_Gran_POTOK(1, gr), gl_Gran_POTOK(2, gr), gl_Gran_POTOK(3, gr), gl_Gran_POTOK(4, gr), &
                            gl_Gran_POTOK(5, gr), gl_Gran_POTOK(6, gr), gl_Gran_POTOK(7, gr), gl_Gran_POTOK(8, gr), gl_Gran_POTOK(9, gr), &
                            gl_Gran_POTOK(10, gr)
                    write(*,*) "__B__"
                    write(*,*) gl_Gran_square2(gr, now)
                    write(*,*) "__C__"
                    write(*,*) dist, wc, dsl, dsp, loc_time, min_sort2 
                    write(*,*) "__D__"
                    write(*,*) gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now)
                    write(*,*) "__E__"
                    write(*,*) qqq1(1), qqq1(2), qqq1(3), qqq1(4), qqq1(5), qqq1(6), qqq1(7), qqq1(8), qqq1(9)
                    write(*,*) "__F__"
                    write(*,*) qqq2(1), qqq2(2), qqq2(3), qqq2(4), qqq2(5), qqq2(6), qqq2(7), qqq2(8), qqq2(9)
                end if
                
                
                ! Теперь сравниваем с глобальным временем
                if (gl_Cell_center2(1, s1, now) > -100000.0) then
                    if(loc_time < time) then
                        !$omp critical
                        if(loc_time < time) then
                            time = loc_time
                            ! time = min(time, loc_time )   ! REDUCTION
                            mincell = gr
                            min_sort = min_sort2
                        end if
                        !$omp end critical
                    end if
                end if
                
                gl_Gran_POTOK_MF(:, 4, gr) = POTOK_MF * gl_Gran_square2(gr, now)


            end do
            !$omp end do

            
            !$omp barrier ! Синхронизация нитей
            

            !!$omp single
            !    if (mod(st, 5) == 0)  print*, st, time
            !!$omp end single

            ! Теперь цикл по ячейкам
            !$omp do private(POTOK, Volume, Volume2, qqq, i, j, ro3, u3, v3, w3, p3, Q3, bx3, by3, bz3, POTOK_MF_all, zone, l_1, fluid1, SOURSE, sks, ijk)
            do gr = 1, ncell
                if(gl_Cell_info(gr) == 0) CYCLE
                l_1 = .TRUE.
                if (norm2(gl_Cell_center2(:, gr, now)) <= par_R0 + (par_R_character - par_R0) * (3.0_8/par_n_TS)**par_kk1) l_1 = .FALSE.    ! Не считаем внутри сферы
                POTOK = 0.0
                sks = 0.0
                SOURSE = 0.0
                POTOK_MF_all = 0.0
                Volume = gl_Cell_Volume2(gr, now)
                Volume2 = gl_Cell_Volume2(gr, now2)
                qqq = gl_Cell_par(:, gr)
                fluid1 = gl_Cell_par_MF(:, :, gr)
                
                ! Просуммируем потоки через грани
                do i = 1, 6
                    j = gl_Cell_gran(i, gr)
                    if (j == 0) CYCLE
                    if (gl_Gran_neighbour(1, j) == gr) then
                        POTOK = POTOK + gl_Gran_POTOK(1:9, j)
                        sks = sks + gl_Gran_POTOK(10, j)
                        POTOK_MF_all(:, 1) = POTOK_MF_all(:, 1) +  gl_Gran_POTOK_MF(:, 1, j)
                        POTOK_MF_all(:, 2) = POTOK_MF_all(:, 2) +  gl_Gran_POTOK_MF(:, 2, j)
                        POTOK_MF_all(:, 3) = POTOK_MF_all(:, 3) +  gl_Gran_POTOK_MF(:, 3, j)
                        POTOK_MF_all(:, 4) = POTOK_MF_all(:, 4) +  gl_Gran_POTOK_MF(:, 4, j)
                    else
                        POTOK = POTOK - gl_Gran_POTOK(1:9, j)
                        sks = sks - gl_Gran_POTOK(10, j)
                        POTOK_MF_all(:, 1) = POTOK_MF_all(:, 1) -  gl_Gran_POTOK_MF(:, 1, j)
                        POTOK_MF_all(:, 2) = POTOK_MF_all(:, 2) -  gl_Gran_POTOK_MF(:, 2, j)
                        POTOK_MF_all(:, 3) = POTOK_MF_all(:, 3) -  gl_Gran_POTOK_MF(:, 3, j)
                        POTOK_MF_all(:, 4) = POTOK_MF_all(:, 4) -  gl_Gran_POTOK_MF(:, 4, j)
                    end if
                end do


                ! Определяем зону в которой находимся
                if(qqq(9)/qqq(1) < 50.0) then

                    ! Перенормируем параметры
                    !qqq(2:4) = qqq(2:4) * (par_chi/par_chi_real)
                    !qqq(1) = qqq(1) / (par_chi/par_chi_real)**2

                    if(norm2(qqq(2:4))/sqrt(ggg*qqq(5)/qqq(1)) > 1.3) then
                        zone = 1
                    else
                        zone = 2
                    end if
                else
                    if(norm2(qqq(2:4))/sqrt(ggg*qqq(5)/qqq(1)) > 1.1) then
                        zone = 4
                    else
                        zone = 3
                    end if
                end if

                ! Перенормируем первую жидкость
                !fluid1(2:4, 1) = fluid1(2:4, 1) * (par_chi/par_chi_real)
                !fluid1(1, 1) = fluid1(1, 1) / (par_chi/par_chi_real)**2
                !
                !fluid1(2:4, 2) = fluid1(2:4, 2) * (3.0_8)
                !fluid1(1, 2) = fluid1(1, 2) / (3.0_8)**2

                call Calc_sourse_MF(qqq, fluid1, SOURSE, zone)  ! Вычисляем источники

                ! Перенормируем первую жидкость обратно
                !fluid1(2:4, 1) = fluid1(2:4, 1) / (par_chi/par_chi_real)
                !fluid1(1, 1) = fluid1(1, 1) * (par_chi/par_chi_real)**2
                !SOURSE(5, 2) = SOURSE(5, 2)/ (par_chi/par_chi_real)
                !SOURSE(1, 2) = SOURSE(1, 2)* (par_chi/par_chi_real)
                !
                !fluid1(2:4, 2) = fluid1(2:4, 2) / (3.0_8)
                !fluid1(1, 2) = fluid1(1, 2) * (3.0_8)**2
                !SOURSE(5, 3) = SOURSE(5, 3)/ (3.0_8)
                !SOURSE(1, 3) = SOURSE(1, 3)* (3.0_8)

                ! Перенормируем обратно
                !if(zone == 1 .or. zone == 2) then
                !    qqq(2:4) = qqq(2:4) / (par_chi/par_chi_real)
                !    qqq(1) = qqq(1) * (par_chi/par_chi_real)**2
                !    SOURSE(5, 1) = SOURSE(5, 1)/ (par_chi/par_chi_real)
                !end if


                if (l_1 == .TRUE.) then
                    ro3 = qqq(1)* Volume / Volume2 - TT * POTOK(1) / Volume2
                    Q3 = qqq(9)* Volume / Volume2 - TT * POTOK(9) / Volume2
                    if (ro3 <= 0.0_8) then
                        print*, "Ro < 0  1490 to\from	", ro3, qqq(1)
                        print*, "Center = ", gl_Cell_center(1, gr), gl_Cell_center(2, gr), gl_Cell_center(3, gr)
                        print*, "Plotnost' opustilas' nizhe nulya v yachejke pod nomerom = ", gr
                        !print*, "gl_Cell_info(gr)", gl_Cell_info(gr)
                        print*, "time = ", TT
                        print*, "Ob'yomy =  ", Volume, Volume2
                        print*, "Potok = ", POTOK(1), POTOK(2), POTOK(3), POTOK(4)
                        print*, "Potok = ", POTOK(5), POTOK(6), POTOK(7), POTOK(8), POTOK(9)
                        ro3 = 0.1
                        !pause
                    end if
                    u3 = (qqq(1) * qqq(2)* Volume / Volume2 - TT * (POTOK(2) + (qqq(6)/cpi4) * sks) / Volume2 + TT * SOURSE(2, 1)) / ro3
                    v3 = (qqq(1) * qqq(3)* Volume / Volume2 - TT * (POTOK(3) + (qqq(7)/cpi4) * sks) / Volume2 + TT * SOURSE(3, 1)) / ro3
                    w3 = (qqq(1) * qqq(4)* Volume / Volume2 - TT * (POTOK(4) + (qqq(8)/cpi4) * sks) / Volume2 + TT * SOURSE(4, 1)) / ro3
                    
                    bx3 = qqq(6) * Volume / Volume2 - TT * (POTOK(6) + qqq(2) * sks) / Volume2
                    by3 = qqq(7) * Volume / Volume2 - TT * (POTOK(7) + qqq(3) * sks) / Volume2
                    bz3 = qqq(8) * Volume / Volume2 - TT * (POTOK(8) + qqq(4) * sks) / Volume2
                    
                    p3 = ((  ( qqq(5) / (ggg - 1.0) + 0.5 * qqq(1) * norm2(qqq(2:4))**2 + (norm2(qqq(6:8))**2) / 25.13274122871834590768 )* Volume / Volume2 &
                        - TT * ( POTOK(5) + (DOT_PRODUCT(qqq(2:4), qqq(6:8))/cpi4) * sks)/ Volume2 + TT * SOURSE(5, 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) - (bx3**2 + by3**2 + bz3**2) / 25.13274122871834590768 ) * (ggg - 1.0)
                    
                    if(ieee_is_nan(u3)) then
                        print*, "NUN   1490 u3", u3, qqq(2)
                        print*, "Center = ", gl_Cell_center(1, gr), gl_Cell_center(2, gr), gl_Cell_center(3, gr)
                        print*, "Potok = ", POTOK(1), POTOK(2), POTOK(3), POTOK(4)
                        print*, "Potok = ", POTOK(5), POTOK(6), POTOK(7), POTOK(8), POTOK(9)
                        print*, "Ob'yomy =  ", Volume, Volume2
                        print*, "qqq = ", qqq
                        STOP
                    end if
                    
                    if(ieee_is_nan(bx3)) then
                        print*, "NUN  1490 bx3	", bx3, TT, sks
                        print*, "Center = ", gl_Cell_center(1, gr), gl_Cell_center(2, gr), gl_Cell_center(3, gr)
                        print*, "Potok = ", POTOK(1), POTOK(2), POTOK(3), POTOK(4)
                        print*, "Potok = ", POTOK(5), POTOK(6), POTOK(7), POTOK(8), POTOK(9)
                        print*, "Ob'yomy =  ", Volume, Volume2
                        print*, "qqq = ", qqq
                        STOP
                    end if
                    
                    
                    !if(gr == 50) then
                    !	write(*,*) ro3, qqq(1), Volume, Volume2, TT, POTOK(1)
                    !end if

                    if (p3 <= 0.0_8) then
                        !print*, "p < 0  plasma 2028 ", p3 , gl_Cell_center(:, gr)
                        p3 = 0.000001
                        !pause
                    end if
                    
                    !if(gr == 50) then
                    !	write(*,*) ro3, Q3, u3, v3, w3, p3, bx3, by3, bz3
                    !end if

                    gl_Cell_par(:, gr) = (/ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3/)

                end if

                ! Теперь посчитаем законы сохранения для остальных жидкостей
                

                do i = 1, 4
                    if (i == 1 .and. l_1 == .FALSE.) CYCLE       ! Пропускаем внутреннюю сферу для сорта 1
                    if (l_1 == .FALSE.) SOURSE(:, i + 1) = 0.0       ! Пропускаем внутреннюю сферу для сорта 1
                    ro3 = fluid1(1, i)* Volume / Volume2 - TT * POTOK_MF_all(1, i) / Volume2 + TT * SOURSE(1, i + 1)
                    if (ro3 <= 0.0_8) then
                        print*, "Ro < 0  in sort ", i,  ro3, fluid1(1, i), "|||", gl_Cell_center(:, gr)
                        ro3 = 0.01
                        !print*, "gl_Cell_center(:, gr)", gl_Cell_center(:, gr)
                        !print*, "gl_Cell_info(gr)", gl_Cell_info(gr)
                        !print*, fluid1(:, i)
                        !print*, l_1
                        !print*, SOURSE(:, i + 1)
                        !print*, TT
                        !print*, "POTOK ", POTOK_MF_all(:, i)
                        !pause
                    end if
                    u3 = (fluid1(1, i) * fluid1(2, i)* Volume / Volume2 - TT * POTOK_MF_all(2, i) / Volume2 + TT * SOURSE(2, i + 1)) / ro3
                    v3 = (fluid1(1, i) * fluid1(3, i)* Volume / Volume2 - TT * POTOK_MF_all(3, i) / Volume2 + TT * SOURSE(3, i + 1)) / ro3
                    w3 = (fluid1(1, i) * fluid1(4, i)* Volume / Volume2 - TT * POTOK_MF_all(4, i) / Volume2 + TT * SOURSE(4, i + 1)) / ro3
                    p3 = ((  ( fluid1(5, i) / (ggg - 1.0) + 0.5 * fluid1(1, i) * norm2(fluid1(2:4, i))**2 )* Volume / Volume2 &
                        - TT * POTOK_MF_all(5, i)/ Volume2 + TT * SOURSE(5, i + 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) ) * (ggg - 1.0)
                    if (p3 <= 0.0_8) then
                        !print*, "p < 0  in ", i, p3, gl_Cell_center(:, gr)
                        !print*, fluid1(:, i)
                        !print*, l_1
                        !print*, "Sourse", SOURSE(:, i + 1)
                        !print*, TT
                        !print*, "zone", zone, norm2(qqq(2:4))/sqrt(ggg*qqq(5)/qqq(1))
                        !print*, "POTOK ", POTOK_MF_all(:, i)
                        !print*, "qqq", qqq
                        !pause
                        p3 = 0.000001
                        !pause
                    end if
                    
                    if(ieee_is_nan(u3)) then
                        print*,  "NUN u3  sort	", i
                        print*, "___"
                        print*, fluid1(:, i)
                        print*, "___"
                        print*, POTOK_MF_all(:, i)
                        print*, "___"
                        print*, fluid1(:, i)
                        print*, "___"
                        print*, gl_Cell_center2(:, gr, now)
                        STOP
                    end if
                    
                    
                    if(.False.) then 
                    !if(gr == 206135 .and. i == 1) then
                        write(*,*) "__1__"
                        write(*,*) ro3, u3, v3, w3, p3
                        write(*,*) "__2__"
                        write(*,*) fluid1(1, i), fluid1(2, i), fluid1(3, i), fluid1(4, i), fluid1(5, i)
                        write(*,*) "__3__"
                        write(*,*) POTOK_MF_all(1, i), POTOK_MF_all(2, i), POTOK_MF_all(3, i), POTOK_MF_all(4, i), POTOK_MF_all(5, i)
                        write(*,*) "__4__"
                        write(*,*) Volume, Volume2
                        write(*,*) "__5__"
                        write(*,*)  SOURSE(1, i + 1), SOURSE(2, i + 1), SOURSE(3, i + 1), SOURSE(4, i + 1), SOURSE(5, i + 1)
                        write(*,*) "__6__"
                        write(*,*) TT
                        write(*,*) "_____"
                        
                        do ijk = 1, 6
                    j = gl_Cell_gran(ijk, gr)
                    write(*,*) "*********************"
                    write(*,*) "gran = ", j, "type = ", gl_Gran_type(j), "Info = ",gl_Gran_info(j)
                    if (j == 0) CYCLE
                        
                    if (gl_Gran_neighbour(1, j) == gr) then 
                            write(*,*) "__+__"
                        else
                            write(*,*) "__-__"
                        end if
                        write(*,*) "__A__"
                        write(*,*)  gl_Gran_POTOK(1, j), gl_Gran_POTOK(2, j), gl_Gran_POTOK(3, j), gl_Gran_POTOK(4, j), &
                            gl_Gran_POTOK(5, j), gl_Gran_POTOK(6, j), gl_Gran_POTOK(7, j), gl_Gran_POTOK(8, j), gl_Gran_POTOK(9, j), &
                            gl_Gran_POTOK(10, j)
                        write(*,*) "__B__"
                        write(*,*) gl_Gran_POTOK_MF(1, 1, j), gl_Gran_POTOK_MF(2, 1, j), gl_Gran_POTOK_MF(3, 1, j), &
                            gl_Gran_POTOK_MF(4, 1, j), gl_Gran_POTOK_MF(5, 1, j)
                        write(*,*) "__C__"
                        write(*,*) gl_Gran_POTOK_MF(1, 2, j), gl_Gran_POTOK_MF(2, 2, j), gl_Gran_POTOK_MF(3, 2, j), &
                            gl_Gran_POTOK_MF(4, 2, j), gl_Gran_POTOK_MF(5, 2, j)
                        write(*,*) "__D__"
                        write(*,*) gl_Gran_POTOK_MF(1, 3, j), gl_Gran_POTOK_MF(2, 3, j), gl_Gran_POTOK_MF(3, 3, j), &
                            gl_Gran_POTOK_MF(4, 3, j), gl_Gran_POTOK_MF(5, 3, j)
                        write(*,*) "__F__"
                        write(*,*) gl_Gran_POTOK_MF(1, 4, j), gl_Gran_POTOK_MF(2, 4, j), gl_Gran_POTOK_MF(3, 4, j), &
                            gl_Gran_POTOK_MF(4, 4, j), gl_Gran_POTOK_MF(5, 4, j)
                    end do
                        
                    end if

                    gl_Cell_par_MF(:, i, gr) = (/ro3, u3, v3, w3, p3/)
                end do
                
                !if(gr == 50) then
                !		write(*,*) "____________________"
                !	end if

            end do
            !$omp end do
            
            !$omp barrier ! Синхронизация нитей
            
            !$omp end parallel
            
            
            gl_Vx = 0.0
            gl_Vy = 0.0
            gl_Vz = 0.0
            gl_Point_num = 0
            
        !	print*, "____"
        !print*, gl_Cell_par(:, 10)
        !print*, "____"
        !print*, gl_Cell_par(:, 1000)
        !print*, "____"
        !print*, "____"
            
            gl_x = gl_x2(:, now2)
            gl_y = gl_y2(:, now2)
            gl_z = gl_z2(:, now2)
            gl_Cell_Volume = gl_Cell_Volume2(:, now2)
            gl_Gran_normal = gl_Gran_normal2(:, :, now2)
            gl_Gran_center = gl_Gran_center2(:, :, now2)
            gl_Cell_center = gl_Cell_center2(:, :, now2)
            gl_Gran_square = gl_Gran_square2(:, now2)
            
            !if (mod(step, 100) == 0) then
            !	call Initial_conditions()
            !end if
            
            call Start_MGD_3_inner(5)
            
            if (mod(step, 50000) == 0 .or. step == 1000 .or. step == 3000) then
                print*, "PECHAT"
                call Print_surface_2D()
                call Print_Setka_2D()
                call Print_par_2D()
                call Print_Contact_3D()
                call Print_TS_3D()
            end if
            
            
            
            
        end do


            gl_x = gl_x2(:, 2)
            gl_y = gl_y2(:, 2)
            gl_z = gl_z2(:, 2)
            gl_Cell_Volume = gl_Cell_Volume2(:, 2)
            gl_Gran_normal = gl_Gran_normal2(:, :, 2)
            gl_Gran_center = gl_Gran_center2(:, :, 2)
            gl_Cell_center = gl_Cell_center2(:, :, 2)
            gl_Gran_square = gl_Gran_square2(:, 2)

            end_time = omp_get_wtime()
            print *, "Time work: ", (end_time-start_time)/60.0, "   in minutes"
            


    end subroutine Start_MGD_move
                             
      
    ! ************************************************************************************************************************************************
    !! Блок геометрии, необходимой для физики - Считаем площади, объёмы

    real(8) pure function calc_square(n) ! Вычисляет площадь грани под номером n в общем списке граней
        use STORAGE, only: gl_all_Gran, gl_x, gl_y, gl_z
        implicit none

        integer, intent (in) :: n
        real(8) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4
        real(8) :: v1(3), v2(3), d1, d2

        x1 = gl_x(gl_all_Gran(1, n))
        y1 = gl_y(gl_all_Gran(1, n))
        z1 = gl_z(gl_all_Gran(1, n))

        x2 = gl_x(gl_all_Gran(2, n))
        y2 = gl_y(gl_all_Gran(2, n))
        z2 = gl_z(gl_all_Gran(2, n))

        x3 = gl_x(gl_all_Gran(3, n))
        y3 = gl_y(gl_all_Gran(3, n))
        z3 = gl_z(gl_all_Gran(3, n))

        x4 = gl_x(gl_all_Gran(4, n))
        y4 = gl_y(gl_all_Gran(4, n))
        z4 = gl_z(gl_all_Gran(4, n))

        v1 = (/x3-x1, y3-y1, z3-z1/)
        v2 = (/x4-x2, y4-y2, z4-z2/)

        d1 = norm2(v1)
        d2 = norm2(v2)

        calc_square = 0.5 * d1 * d2 * sqrt(1.0 - (DOT_PRODUCT(v1, v2)/d1/d2)**2 )
    end function calc_square

    pure subroutine Get_center(n, CENTER)  ! Получить массив - центр ячейки
        use STORAGE, only: gl_all_Cell, gl_x, gl_y, gl_z
        implicit none

        integer, intent(in) :: n
        real(8), intent(out) :: CENTER(3)
        real(8) :: c(3)
        integer :: i, j

        CENTER = 0.0

        do i = 1, 8
            j = gl_all_Cell(i, n)
            c = (/gl_x(j), gl_y(j), gl_z(j)/)
            CENTER = CENTER + c
        end do

        CENTER = CENTER/8.0

    end subroutine Get_center

    ! ****************************************************************************************************************************************************
    ! Блок функций для печати и проверки геометрии сетки

    subroutine Print_all_surface(str)  ! Печать поверхностей, которые выделяются (TS, HP, BS)
        
        use GEO_PARAM
        use STORAGE
        use My_func
        implicit none

        !interface
        !    real(8) pure function calc_square(n)
        !    integer, intent(in) :: n
        !    end function
        !end interface

        character, intent(in) :: str   ! Задаёт поверхность, которую надо напечатать буквой (см. ниже)
        integer :: m, n, j, i

        if (str == "C") then ! Contact
            m = size(gl_Contact)
            open(1, file = 'Contact_geometry.txt')
            write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z', 'M'  ZONE T= 'HP', N= ", 4 * (m), ", E =  ", 4 * (m) , ", F=FEPOINT, ET=LINESEG "
            do j = 1, m
                do i = 1, 4
                    n = gl_Contact(j)
                    write(1,*) gl_x(gl_all_Gran(i, n)), gl_y(gl_all_Gran(i, n)), gl_z(gl_all_Gran(i, n)), calc_square(n)
                end do
            end do

            do j = 0, m
                write(1,*) 4 * j + 1, 4 * j + 2
                write(1,*) 4 * j + 2, 4 * j + 3
                write(1,*) 4 * j + 3, 4 * j + 4
                write(1,*) 4 * j + 4, 4 * j + 1
            end do
            close(1)
        elseif(str == "B") then
            m = size(gl_BS)
            open(1, file = 'BS_geometry.txt')
            write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z', 'M'  ZONE T= 'HP', N= ", 4 * (m), ", E =  ", 4 * (m) , ", F=FEPOINT, ET=LINESEG "
            do j = 1, m
                do i = 1, 4
                    n = gl_BS(j)
                    write(1,*) gl_x(gl_all_Gran(i, n)), gl_y(gl_all_Gran(i, n)), gl_z(gl_all_Gran(i, n)), calc_square(n)
                end do
            end do

            do j = 0, m
                write(1,*) 4 * j + 1, 4 * j + 2
                write(1,*) 4 * j + 2, 4 * j + 3
                write(1,*) 4 * j + 3, 4 * j + 4
                write(1,*) 4 * j + 4, 4 * j + 1
            end do
            close(1)
        elseif(str == "T") then
            m = size(gl_TS)
            open(1, file = 'TS_geometry.txt')
            write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z', 'M'  ZONE T= 'HP', N= ", 4 * (m), ", E =  ", 4 * (m) , ", F=FEPOINT, ET=LINESEG "
            do j = 1, m
                do i = 1, 4
                    n = gl_TS(j)
                    write(1,*) gl_x(gl_all_Gran(i, n)), gl_y(gl_all_Gran(i, n)), gl_z(gl_all_Gran(i, n)), calc_square(n)
                end do
            end do

            do j = 0, m
                write(1,*) 4 * j + 1, 4 * j + 2
                write(1,*) 4 * j + 2, 4 * j + 3
                write(1,*) 4 * j + 3, 4 * j + 4
                write(1,*) 4 * j + 4, 4 * j + 1
            end do
            close(1)
        else
            print*, "Print_all_surface: ERROR 1304"
        end if

    end subroutine Print_all_surface

    subroutine Print_Setka_2D()  ! Печатает 2Д сетку с линиями в Техплот
        use GEO_PARAM
        use STORAGE
        implicit none

        integer :: N1, N2, kk, k, i, j, N

        print*, "Print_Setka_2D()"
        N = size(gl_Cell_A(1, :, 1)) * size(gl_Cell_A(:, 1, 1)) + &
            size(gl_Cell_B(1, :, 1)) * size(gl_Cell_B(:, 1, 1)) + &
            size(gl_Cell_C(1, :, 1)) * size(gl_Cell_C(:, 1, 1))
        N = N * 2

        open(1, file = 'print_setka_2D.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y'  ZONE T= 'HP', N= ", 4 * N, ", E =  ", 4 * N , ", F=FEPOINT, ET=LINESEG "


        kk = 1
        N2 = size(gl_Cell_A(1, :, 1))
        N1 = size(gl_Cell_A(:, 1, 1))
        do j = 1, N2
            do i = 1, N1
                do k = 1, 4
                    write(1,*) gl_x(gl_all_Cell(k, gl_Cell_A(i, j, kk) )), gl_y(gl_all_Cell(k, gl_Cell_A(i, j, kk) ))
                end do
            end do
        end do
        N2 = size(gl_Cell_B(1, :, 1))
        N1 = size(gl_Cell_B(:, 1, 1))
        do j = 1, N2
            do i = 1, N1
                do k = 1, 4
                    write(1,*) gl_x(gl_all_Cell(k, gl_Cell_B(i, j, kk) )), gl_y(gl_all_Cell(k, gl_Cell_B(i, j, kk) ))
                end do
            end do
        end do
        N2 = size(gl_Cell_C(1, :, 1))
        N1 = size(gl_Cell_C(:, 1, 1))
        do j = 1, N2
            do i = 1, N1
                do k = 1, 4
                    write(1,*) gl_x(gl_all_Cell(k, gl_Cell_C(i, j, kk) )), gl_y(gl_all_Cell(k, gl_Cell_C(i, j, kk) ))
                end do
            end do
        end do

        ! Нижняя часть сетки

        kk = par_l_phi/2 + 1
        N2 = size(gl_Cell_A(1, :, 1))
        N1 = size(gl_Cell_A(:, 1, 1))
        do j = 1, N2
            do i = 1, N1
                do k = 1, 4
                    write(1,*) gl_x(gl_all_Cell(k, gl_Cell_A(i, j, kk) )), gl_y(gl_all_Cell(k, gl_Cell_A(i, j, kk) ))
                end do
            end do
        end do
        N2 = size(gl_Cell_B(1, :, 1))
        N1 = size(gl_Cell_B(:, 1, 1))
        do j = 1, N2
            do i = 1, N1
                do k = 1, 4
                    write(1,*) gl_x(gl_all_Cell(k, gl_Cell_B(i, j, kk) )), gl_y(gl_all_Cell(k, gl_Cell_B(i, j, kk) ))
                end do
            end do
        end do
        N2 = size(gl_Cell_C(1, :, 1))
        N1 = size(gl_Cell_C(:, 1, 1))
        do j = 1, N2
            do i = 1, N1
                do k = 1, 4
                    write(1,*) gl_x(gl_all_Cell(k, gl_Cell_C(i, j, kk) )), gl_y(gl_all_Cell(k, gl_Cell_C(i, j, kk) ))
                end do
            end do
        end do


        ! Connectivity list
        do j = 0, N
            write(1,*) 4 * j + 1, 4 * j + 2
            write(1,*) 4 * j + 2, 4 * j + 3
            write(1,*) 4 * j + 3, 4 * j + 4
            write(1,*) 4 * j + 4, 4 * j + 1
        end do

        close(1)
	end subroutine Print_Setka_2D

	subroutine Print_Setka_3D_part()  ! Печатает 2Д сетку с линиями в Техплот
        use GEO_PARAM
        use STORAGE
        implicit none

        integer :: N1, N2, kk, k, i, j, N, N3, c


        N = (size( gl_Cell_B(1, 1, :) ) * size( gl_Cell_B(1, :, 1) ) + size( gl_Cell_C(1, 1, :) ) * size( gl_Cell_C(:, 1, 1) ))

        open(1, file = 'Print_Setka_3D_part.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z'  ZONE T= 'HP', N= ", N * 4, ", E =  ", 4 * N , ", F=FEPOINT, ET=LINESEG "

        N3 = size( gl_Cell_B(1, 1, :) )
        N2 = size( gl_Cell_B(1, :, 1) )
        
        do k = 1, N3
            do j = 1, N2
                c = 1
                write(1,*) gl_x(gl_all_Cell(c, gl_Cell_B(13 + par_n_TS, j, k) )), gl_y(gl_all_Cell(c, gl_Cell_B(13 +  + par_n_TS, j, k) )), gl_z(gl_all_Cell(c, gl_Cell_B(13 +  + par_n_TS, j, k) ))
                c = 2
                write(1,*) gl_x(gl_all_Cell(c, gl_Cell_B(13 +  + par_n_TS, j, k) )), gl_y(gl_all_Cell(c, gl_Cell_B(13 +  + par_n_TS, j, k) )), gl_z(gl_all_Cell(c, gl_Cell_B(13 +  + par_n_TS, j, k) ))
                c = 6
                write(1,*) gl_x(gl_all_Cell(c, gl_Cell_B(13 +  + par_n_TS, j, k) )), gl_y(gl_all_Cell(c, gl_Cell_B(13 +  + par_n_TS, j, k) )), gl_z(gl_all_Cell(c, gl_Cell_B(13 +  + par_n_TS, j, k) ))
                c = 5
                write(1,*) gl_x(gl_all_Cell(c, gl_Cell_B(13 +  + par_n_TS, j, k) )), gl_y(gl_all_Cell(c, gl_Cell_B(13 +  + par_n_TS, j, k) )), gl_z(gl_all_Cell(c, gl_Cell_B(13 +  + par_n_TS, j, k) ))
                
            end do
        end do
        
        N3 = size( gl_Cell_C(1, 1, :) )
        N1 = size( gl_Cell_C(:, 1, 1) )
        
        do k = 1, N3
            do i = 1, N1
                c = 3
                write(1,*) gl_x(gl_all_Cell(c, gl_Cell_C(i, 13, k) )), gl_y(gl_all_Cell(c, gl_Cell_C(i, 13, k) )), gl_z(gl_all_Cell(c, gl_Cell_C(i, 13, k) ))
                c = 7
                write(1,*) gl_x(gl_all_Cell(c, gl_Cell_C(i, 13, k) )), gl_y(gl_all_Cell(c, gl_Cell_C(i, 13, k) )), gl_z(gl_all_Cell(c, gl_Cell_C(i, 13, k) ))
                c = 8
                write(1,*) gl_x(gl_all_Cell(c, gl_Cell_C(i, 13, k) )), gl_y(gl_all_Cell(c, gl_Cell_C(i, 13, k) )), gl_z(gl_all_Cell(c, gl_Cell_C(i, 13, k) ))
                c = 4
                write(1,*) gl_x(gl_all_Cell(c, gl_Cell_C(i, 13, k) )), gl_y(gl_all_Cell(c, gl_Cell_C(i, 13, k) )), gl_z(gl_all_Cell(c, gl_Cell_C(i, 13, k) ))
            end do
        end do
        
        
        


        ! Connectivity list
        do j = 0, N - 1
            write(1,*) 4 * j + 1, 4 * j + 2
            write(1,*) 4 * j + 2, 4 * j + 3
            write(1,*) 4 * j + 3, 4 * j + 4
            write(1,*) 4 * j + 4, 4 * j + 1
        end do

        close(1)


	end subroutine Print_Setka_3D_part
	
	subroutine Print_Setka_y_2D()  ! Печатает 2Д сетку с линиями в Техплот
        use GEO_PARAM
        use STORAGE
        implicit none

        integer :: N1, N2, kk, k, i, j, N


        N = size(gl_Cell_A(1, :, 1)) * size(gl_Cell_A(:, 1, 1)) + &
            size(gl_Cell_B(1, :, 1)) * size(gl_Cell_B(:, 1, 1)) + &
            size(gl_Cell_C(1, :, 1)) * size(gl_Cell_C(:, 1, 1))
        N = N * 2

        open(1, file = 'print_setka_y_2D.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Z'  ZONE T= 'HP', N= ", 4 * N, ", E =  ", 4 * N , ", F=FEPOINT, ET=LINESEG "


        kk = par_l_phi/4
        N2 = size(gl_Cell_A(1, :, 1))
        N1 = size(gl_Cell_A(:, 1, 1))
        do j = 1, N2
            do i = 1, N1
                do k = 1, 4
                    write(1,*) gl_x(gl_all_Cell(k, gl_Cell_A(i, j, kk) )), gl_z(gl_all_Cell(k, gl_Cell_A(i, j, kk) ))
                end do
            end do
        end do
        N2 = size(gl_Cell_B(1, :, 1))
        N1 = size(gl_Cell_B(:, 1, 1))
        do j = 1, N2
            do i = 1, N1
                do k = 1, 4
                    write(1,*) gl_x(gl_all_Cell(k, gl_Cell_B(i, j, kk) )), gl_z(gl_all_Cell(k, gl_Cell_B(i, j, kk) ))
                end do
            end do
        end do
        N2 = size(gl_Cell_C(1, :, 1))
        N1 = size(gl_Cell_C(:, 1, 1))
        do j = 1, N2
            do i = 1, N1
                do k = 1, 4
                    write(1,*) gl_x(gl_all_Cell(k, gl_Cell_C(i, j, kk) )), gl_z(gl_all_Cell(k, gl_Cell_C(i, j, kk) ))
                end do
            end do
        end do

        ! Нижняя часть сетки

        kk = 3 * par_l_phi/4 + 1
        N2 = size(gl_Cell_A(1, :, 1))
        N1 = size(gl_Cell_A(:, 1, 1))
        do j = 1, N2
            do i = 1, N1
                do k = 1, 4
                    write(1,*) gl_x(gl_all_Cell(k, gl_Cell_A(i, j, kk) )), gl_z(gl_all_Cell(k, gl_Cell_A(i, j, kk) ))
                end do
            end do
        end do
        N2 = size(gl_Cell_B(1, :, 1))
        N1 = size(gl_Cell_B(:, 1, 1))
        do j = 1, N2
            do i = 1, N1
                do k = 1, 4
                    write(1,*) gl_x(gl_all_Cell(k, gl_Cell_B(i, j, kk) )), gl_z(gl_all_Cell(k, gl_Cell_B(i, j, kk) ))
                end do
            end do
        end do
        N2 = size(gl_Cell_C(1, :, 1))
        N1 = size(gl_Cell_C(:, 1, 1))
        do j = 1, N2
            do i = 1, N1
                do k = 1, 4
                    write(1,*) gl_x(gl_all_Cell(k, gl_Cell_C(i, j, kk) )), gl_z(gl_all_Cell(k, gl_Cell_C(i, j, kk) ))
                end do
            end do
        end do


        ! Connectivity list
        do j = 0, N
            write(1,*) 4 * j + 1, 4 * j + 2
            write(1,*) 4 * j + 2, 4 * j + 3
            write(1,*) 4 * j + 3, 4 * j + 4
            write(1,*) 4 * j + 4, 4 * j + 1
        end do

        close(1)
	end subroutine Print_Setka_y_2D
	
	subroutine Print_Contact_3D()
        use GEO_PARAM
        use STORAGE
        implicit none
        integer(4) :: N1, N2, N3, j, k, cel, gr, kk, M1, M2, M3
        real(8) :: c(3)
        
        N3 = size(gl_Cell_A(1, 1, :))
        N2 = size(gl_Cell_A(1, :, 1))
        N1 = size(gl_Cell_A(:, 1, 1))
        
        M3 = size(gl_Cell_C(1, 1, :))
        M2 = size(gl_Cell_C(1, :, 1))
        M1 = size(gl_Cell_C(:, 1, 1))
        
        open(1, file = 'Print_Contact_3D.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z', 'Vn_1'  ZONE T= 'HP', N= ",  4 * N3 * (N2 - 1) + 4 * M3 * (M2 - 1) + 4 * N3 , ", E =  ", N3 * (N2 - 1) + M3 * (M2 - 1) + N3, ", F=FEPOINT, ET=quadrilateral "
        
        
        do k = 1, N3
            do j = 1, N2 - 1
                kk = k + 1
                if (kk > N3) kk = 1
                
                cel = gl_Cell_A(par_n_HP - 1, j, k)
                gr = gl_Cell_gran(1, cel)
                c = gl_Gran_center(:, gr)
                write(1,*) c(1), c(2), c(3), 1.0_8
                
                cel = gl_Cell_A(par_n_HP - 1, j + 1, k)
                gr = gl_Cell_gran(1, cel)
                c = gl_Gran_center(:, gr)
                write(1,*) c(1), c(2), c(3), 1.0_8
                
                cel = gl_Cell_A(par_n_HP - 1, j + 1, kk)
                gr = gl_Cell_gran(1, cel)
                c = gl_Gran_center(:, gr)
                write(1,*) c(1), c(2), c(3), 1.0_8
                
                cel = gl_Cell_A(par_n_HP - 1, j, kk)
                gr = gl_Cell_gran(1, cel)
                c = gl_Gran_center(:, gr)
                write(1,*) c(1), c(2), c(3), 1.0_8
                
            end do	
        end do
        
        
        do k = 1, M3
            do j = 1, M2 - 1
                kk = k + 1
                if (kk > N3) kk = 1
                
                cel = gl_Cell_C(par_n_HP - par_n_TS, j, k)
                gr = gl_Cell_gran(1, cel)
                c = gl_Gran_center(:, gr)
                write(1,*) c(1), c(2), c(3), 1.0_8
                
                cel = gl_Cell_C(par_n_HP - par_n_TS, j + 1, k)
                gr = gl_Cell_gran(1, cel)
                c = gl_Gran_center(:, gr)
                write(1,*) c(1), c(2), c(3), 1.0_8
                
                cel = gl_Cell_C(par_n_HP - par_n_TS, j + 1, kk)
                gr = gl_Cell_gran(1, cel)
                c = gl_Gran_center(:, gr)
                write(1,*) c(1), c(2), c(3), 1.0_8
                
                cel = gl_Cell_C(par_n_HP - par_n_TS, j, kk)
                gr = gl_Cell_gran(1, cel)
                c = gl_Gran_center(:, gr)
                write(1,*) c(1), c(2), c(3), 1.0_8
                
            end do	
        end do
        
        do k = 1, N3
                kk = k + 1
                if (kk > N3) kk = 1
                
                cel = gl_Cell_A(par_n_HP - 1, N2, k)
                gr = gl_Cell_gran(1, cel)
                c = gl_Gran_center(:, gr)
                write(1,*) c(1), c(2), c(3), 1.0_8
                
                cel = gl_Cell_C(par_n_HP - par_n_TS, 1, k)
                gr = gl_Cell_gran(1, cel)
                c = gl_Gran_center(:, gr)
                write(1,*) c(1), c(2), c(3), 1.0_8
                
                cel = gl_Cell_C(par_n_HP - par_n_TS, 1, kk)
                gr = gl_Cell_gran(1, cel)
                c = gl_Gran_center(:, gr)
                write(1,*) c(1), c(2), c(3), 1.0_8
                
                cel = gl_Cell_A(par_n_HP - 1, N2, kk)
                gr = gl_Cell_gran(1, cel)
                c = gl_Gran_center(:, gr)
                write(1,*) c(1), c(2), c(3), 1.0_8
        end do
        
        do k = 1, N3 * (N2 - 1) + M3 * (M2 - 1) + N3
            write(1,*) 4 * k - 3, 4 * k - 2, 4 * k - 1, 4 * k
        end do
        
        
        close(1)
	
	end subroutine Print_Contact_3D

    subroutine Print_BS_3D()
        use GEO_PARAM
        use STORAGE
        implicit none
        integer(4) :: N1, N2, N3, j, k, cel, gr, kk
        real(8) :: c(3)
        
        N3 = size(gl_Cell_A(1, 1, :))
        N2 = size(gl_Cell_A(1, :, 1))
        N1 = size(gl_Cell_A(:, 1, 1))
        
        
        open(1, file = 'Print_BS_3D.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z', 'Vn_1'  ZONE T= 'HP', N= ",  4 * N3 * (N2 - 1), ", E =  ", N3 * (N2 - 1), ", F=FEPOINT, ET=quadrilateral "
        
        
        do k = 1, N3
            do j = 1, N2 - 1
                kk = k + 1
                if (kk > N3) kk = 1
                
                cel = gl_Cell_A(par_n_BS - 1, j, k)
                gr = gl_Cell_gran(1, cel)
                c = gl_Gran_center(:, gr)
                write(1,*) c(1), c(2), c(3), 1.0_8
                
                cel = gl_Cell_A(par_n_BS - 1, j + 1, k)
                gr = gl_Cell_gran(1, cel)
                c = gl_Gran_center(:, gr)
                write(1,*) c(1), c(2), c(3), 1.0_8
                
                cel = gl_Cell_A(par_n_BS - 1, j + 1, kk)
                gr = gl_Cell_gran(1, cel)
                c = gl_Gran_center(:, gr)
                write(1,*) c(1), c(2), c(3), 1.0_8
                
                cel = gl_Cell_A(par_n_BS - 1, j, kk)
                gr = gl_Cell_gran(1, cel)
                c = gl_Gran_center(:, gr)
                write(1,*) c(1), c(2), c(3), 1.0_8
                
            end do	
        end do
        
        
        do k = 1, N3 * (N2 - 1)
            write(1,*) 4 * k - 3, 4 * k - 2, 4 * k - 1, 4 * k
        end do
        
        
        close(1)
	
	end subroutine Print_BS_3D
	
	subroutine Print_TS_3D()
        use GEO_PARAM
        use STORAGE
        implicit none
        integer(4) :: N1, N2, N3, j, k, cel, gr, kk, N, E
        real(8) :: c(3)
        
        N3 = size(gl_Cell_A(1, 1, :))
        N2 = size(gl_Cell_A(1, :, 1))
        N1 = size(gl_Cell_A(:, 1, 1))
        
        N = 4 * N3 * (N2 - 1) + 4 * size(gl_Cell_B(1, 1, :)) * (size(gl_Cell_B(1, :, 1)) - 1) + size(gl_Cell_B(1, 1, :)) * 4
        
        E = N3 * (N2 - 1) + size(gl_Cell_B(1, 1, :)) * (size(gl_Cell_B(1, :, 1)) - 1) + size(gl_Cell_B(1, 1, :))
        
        open(1, file = 'Print_TS_3D.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z', 'I'  ZONE T= 'HP', N= ",  N , ", E =  ", E, ", F=FEPOINT, ET=quadrilateral "
        
        
        do k = 1, N3
            do j = 1, N2 - 1
                kk = k + 1
                if (kk > N3) kk = 1
                
                cel = gl_Cell_A(par_n_TS - 1, j, k)
                gr = gl_Cell_gran(1, cel)
                c = gl_Gran_center(:, gr)
                write(1,*) c(1), c(2), c(3), 1.0_8
                
                cel = gl_Cell_A(par_n_TS - 1, j + 1, k)
                gr = gl_Cell_gran(1, cel)
                c = gl_Gran_center(:, gr)
                write(1,*) c(1), c(2), c(3), 1.0_8
                
                cel = gl_Cell_A(par_n_TS - 1, j + 1, kk)
                gr = gl_Cell_gran(1, cel)
                c = gl_Gran_center(:, gr)
                write(1,*) c(1), c(2), c(3), 1.0_8
                
                cel = gl_Cell_A(par_n_TS - 1, j, kk)
                gr = gl_Cell_gran(1, cel)
                c = gl_Gran_center(:, gr)
                write(1,*) c(1), c(2), c(3), 1.0_8
                
            end do	
        end do
        
        N3 = size(gl_Cell_B(1, 1, :))
        N2 = size(gl_Cell_B(1, :, 1))
        N1 = size(gl_Cell_B(:, 1, 1))
        
        do k = 1, N3
            do j = 1, N2 - 1
                kk = k + 1
                if (kk > N3) kk = 1
                
                cel = gl_Cell_B(par_n_TS - 1, j, k)
                gr = gl_Cell_gran(1, cel)
                c = gl_Gran_center(:, gr)
                write(1,*) c(1), c(2), c(3), 1.0_8
                
                cel = gl_Cell_B(par_n_TS - 1, j + 1, k)
                gr = gl_Cell_gran(1, cel)
                c = gl_Gran_center(:, gr)
                write(1,*) c(1), c(2), c(3), 1.0_8
                
                cel = gl_Cell_B(par_n_TS - 1, j + 1, kk)
                gr = gl_Cell_gran(1, cel)
                c = gl_Gran_center(:, gr)
                write(1,*) c(1), c(2), c(3), 1.0_8
                
                cel = gl_Cell_B(par_n_TS - 1, j, kk)
                gr = gl_Cell_gran(1, cel)
                c = gl_Gran_center(:, gr)
                write(1,*) c(1), c(2), c(3), 1.0_8
                
            end do	
        end do
        
        N2 = size(gl_Cell_A(1, :, 1))
        
        do k = 1, N3
            kk = k + 1
            if (kk > N3) kk = 1
                
            cel = gl_Cell_A(par_n_TS - 1, N2, k)
            gr = gl_Cell_gran(1, cel)
            c = gl_Gran_center(:, gr)
            write(1,*) c(1), c(2), c(3), 1.0_8
                
            cel = gl_Cell_B(par_n_TS - 1, size(gl_Cell_B(1, :, 1)), k)
            gr = gl_Cell_gran(1, cel)
            c = gl_Gran_center(:, gr)
            write(1,*) c(1), c(2), c(3), 1.0_8
                
            cel = gl_Cell_B(par_n_TS - 1, size(gl_Cell_B(1, :, 1)), kk)
            gr = gl_Cell_gran(1, cel)
            c = gl_Gran_center(:, gr)
            write(1,*) c(1), c(2), c(3), 1.0_8
                
            cel = gl_Cell_A(par_n_TS - 1, N2, kk)
            gr = gl_Cell_gran(1, cel)
            c = gl_Gran_center(:, gr)
            write(1,*) c(1), c(2), c(3), 1.0_8
        end do
        
        
        do k = 1, E
            write(1,*) 4 * k - 3, 4 * k - 2, 4 * k - 1, 4 * k
        end do
        
        
        close(1)
	
	end subroutine Print_TS_3D

    subroutine Print_surface_2D()
        use GEO_PARAM
        use STORAGE
        implicit none
        integer :: kk, j, yzel, N, M, L
        ! Body of Print_surface
        kk = 1
        N = size(gl_RAY_A(1, :, kk)) + size(gl_RAY_B(1, :, kk)) + size(gl_RAY_K(1, :, kk))
        M = size(gl_RAY_A(1, :, kk)) + size(gl_RAY_C(1, :, kk)) + size(gl_RAY_O(1, :, kk))
        L = M
        print*, "Print_surface_2D()"
        open(1, file = 'print_surface_2D.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y'  ZONE T= 'HP', N= ", 2 * N + 2 * M + 2 * L, ", E =  ", 2 * (N - 1) + 2 * (M - 1) + &
            + 2 * (L - 1), ", F=FEPOINT, ET=LINESEG "

        kk = 1
        do j = 1, size(gl_RAY_A(1, :, kk))
            yzel = gl_RAY_A(par_n_TS, j, kk)
            write(1,*) gl_x(yzel), gl_y(yzel)
        end do

        do j = 1, size(gl_RAY_B(1, :, kk))
            yzel = gl_RAY_B(par_n_TS, j, kk)
            write(1,*) gl_x(yzel), gl_y(yzel)
        end do

        do j = size(gl_RAY_K(1, :, kk)), 1, -1
            yzel = gl_RAY_K(par_n_TS, j, kk)
            write(1,*) gl_x(yzel), gl_y(yzel)
        end do

        kk = par_l_phi/2 + 1
        do j = 1, size(gl_RAY_A(1, :, kk))
            yzel = gl_RAY_A(par_n_TS, j, kk)
            write(1,*) gl_x(yzel), gl_y(yzel)
        end do

        do j = 1, size(gl_RAY_B(1, :, kk))
            yzel = gl_RAY_B(par_n_TS, j, kk)
            write(1,*) gl_x(yzel), gl_y(yzel)
        end do

        do j = size(gl_RAY_K(1, :, kk)), 1, -1
            yzel = gl_RAY_K(par_n_TS, j, kk)
            write(1,*) gl_x(yzel), gl_y(yzel)
        end do

        ! Контакт
        kk = 1
        do j = 1, size(gl_RAY_A(1, :, kk))
            yzel = gl_RAY_A(par_n_HP, j, kk)
            write(1,*) gl_x(yzel), gl_y(yzel)
        end do

        do j = 1, size(gl_RAY_C(1, :, kk))
            yzel = gl_RAY_C(1, j, kk)
            write(1,*) gl_x(yzel), gl_y(yzel)
        end do

        do j =  1, size(gl_RAY_O(1, :, kk))
            yzel = gl_RAY_O(1, j, kk)
            write(1,*) gl_x(yzel), gl_y(yzel)
        end do

        kk = par_l_phi/2 + 1
        do j = 1, size(gl_RAY_A(1, :, kk))
            yzel = gl_RAY_A(par_n_HP, j, kk)
            write(1,*) gl_x(yzel), gl_y(yzel)
        end do

        do j = 1, size(gl_RAY_C(1, :, kk))
            yzel = gl_RAY_C(1, j, kk)
            write(1,*) gl_x(yzel), gl_y(yzel)
        end do

        do j =  1, size(gl_RAY_O(1, :, kk))
            yzel = gl_RAY_O(1, j, kk)
            write(1,*) gl_x(yzel), gl_y(yzel)
        end do

        ! BS
        kk = 1
        do j = 1, size(gl_RAY_A(1, :, kk))
            yzel = gl_RAY_A(par_n_BS, j, kk)
            write(1,*) gl_x(yzel), gl_y(yzel)
        end do

        do j = 1, size(gl_RAY_C(1, :, kk))
            yzel = gl_RAY_C(par_n_BS - par_n_HP + 1, j, kk)
            write(1,*) gl_x(yzel), gl_y(yzel)
        end do

        do j =  1, size(gl_RAY_O(1, :, kk))
            yzel = gl_RAY_O(par_n_BS - par_n_HP + 1, j, kk)
            write(1,*) gl_x(yzel), gl_y(yzel)
        end do

        kk = par_l_phi/2 + 1
        do j = 1, size(gl_RAY_A(1, :, kk))
            yzel = gl_RAY_A(par_n_BS, j, kk)
            write(1,*) gl_x(yzel), gl_y(yzel)
        end do

        do j = 1, size(gl_RAY_C(1, :, kk))
            yzel = gl_RAY_C(par_n_BS - par_n_HP + 1, j, kk)
            write(1,*) gl_x(yzel), gl_y(yzel)
        end do

        do j =  1, size(gl_RAY_O(1, :, kk))
            yzel = gl_RAY_O(par_n_BS - par_n_HP + 1, j, kk)
            write(1,*) gl_x(yzel), gl_y(yzel)
        end do


        ! Connectivity list
        do j = 1, N - 1
            write(1,*) j, j + 1
        end do

        do j = N + 1, 2 * N - 1
            write(1,*) j, j + 1
        end do


        do j = 2 * N + 1, 2 * N + 1 + M - 2
            write(1,*) j, j + 1
        end do

        do j = 2 * N + 1 + M, 2 * N + 1 + 2 * M - 2
            write(1,*) j, j + 1
        end do

        do j = 2 * N + 1 + 2 * M, 2 * N + 1 + 2 * M + L - 2
            write(1,*) j, j + 1
        end do

        do j = 2 * N + 1 + 2 * M + L, 2 * N + 1 + 2 * M + L + L - 2
            write(1,*) j, j + 1
        end do

        close(1)
	end subroutine Print_surface_2D
	
	subroutine Print_surface_y_2D()
        use GEO_PARAM
        use STORAGE
        implicit none
        integer :: kk, j, yzel, N, M, L
        ! Body of Print_surface
        kk = par_l_phi/4
        N = size(gl_RAY_A(1, :, kk)) + size(gl_RAY_B(1, :, kk)) + size(gl_RAY_K(1, :, kk))
        M = size(gl_RAY_A(1, :, kk)) + size(gl_RAY_C(1, :, kk)) + size(gl_RAY_O(1, :, kk))
        L = M

        open(1, file = 'print_surface_y_2D.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Z'  ZONE T= 'HP', N= ", 2 * N + 2 * M + 2 * L, ", E =  ", 2 * (N - 1) + 2 * (M - 1) + &
            + 2 * (L - 1), ", F=FEPOINT, ET=LINESEG "

        kk = par_l_phi/4
        do j = 1, size(gl_RAY_A(1, :, kk))
            yzel = gl_RAY_A(par_n_TS, j, kk)
            write(1,*) gl_x(yzel), gl_z(yzel)
        end do

        do j = 1, size(gl_RAY_B(1, :, kk))
            yzel = gl_RAY_B(par_n_TS, j, kk)
            write(1,*) gl_x(yzel), gl_z(yzel)
        end do

        do j = size(gl_RAY_K(1, :, kk)), 1, -1
            yzel = gl_RAY_K(par_n_TS, j, kk)
            write(1,*) gl_x(yzel), gl_z(yzel)
        end do

        kk = 3 * par_l_phi/4 + 1
        do j = 1, size(gl_RAY_A(1, :, kk))
            yzel = gl_RAY_A(par_n_TS, j, kk)
            write(1,*) gl_x(yzel), gl_z(yzel)
        end do

        do j = 1, size(gl_RAY_B(1, :, kk))
            yzel = gl_RAY_B(par_n_TS, j, kk)
            write(1,*) gl_x(yzel), gl_z(yzel)
        end do

        do j = size(gl_RAY_K(1, :, kk)), 1, -1
            yzel = gl_RAY_K(par_n_TS, j, kk)
            write(1,*) gl_x(yzel), gl_z(yzel)
        end do

        ! Контакт
        kk = par_l_phi/4
        do j = 1, size(gl_RAY_A(1, :, kk))
            yzel = gl_RAY_A(par_n_HP, j, kk)
            write(1,*) gl_x(yzel), gl_z(yzel)
        end do

        do j = 1, size(gl_RAY_C(1, :, kk))
            yzel = gl_RAY_C(1, j, kk)
            write(1,*) gl_x(yzel), gl_z(yzel)
        end do

        do j =  1, size(gl_RAY_O(1, :, kk))
            yzel = gl_RAY_O(1, j, kk)
            write(1,*) gl_x(yzel), gl_z(yzel)
        end do

        kk = 3 * par_l_phi/4 + 1
        do j = 1, size(gl_RAY_A(1, :, kk))
            yzel = gl_RAY_A(par_n_HP, j, kk)
            write(1,*) gl_x(yzel), gl_z(yzel)
        end do

        do j = 1, size(gl_RAY_C(1, :, kk))
            yzel = gl_RAY_C(1, j, kk)
            write(1,*) gl_x(yzel), gl_z(yzel)
        end do

        do j =  1, size(gl_RAY_O(1, :, kk))
            yzel = gl_RAY_O(1, j, kk)
            write(1,*) gl_x(yzel), gl_z(yzel)
        end do

        ! BS
        kk = par_l_phi/4
        do j = 1, size(gl_RAY_A(1, :, kk))
            yzel = gl_RAY_A(par_n_BS, j, kk)
            write(1,*) gl_x(yzel), gl_z(yzel)
        end do

        do j = 1, size(gl_RAY_C(1, :, kk))
            yzel = gl_RAY_C(par_n_BS - par_n_HP + 1, j, kk)
            write(1,*) gl_x(yzel), gl_z(yzel)
        end do

        do j =  1, size(gl_RAY_O(1, :, kk))
            yzel = gl_RAY_O(par_n_BS - par_n_HP + 1, j, kk)
            write(1,*) gl_x(yzel), gl_z(yzel)
        end do

        kk = 3 * par_l_phi/4 + 1
        do j = 1, size(gl_RAY_A(1, :, kk))
            yzel = gl_RAY_A(par_n_BS, j, kk)
            write(1,*) gl_x(yzel), gl_z(yzel)
        end do

        do j = 1, size(gl_RAY_C(1, :, kk))
            yzel = gl_RAY_C(par_n_BS - par_n_HP + 1, j, kk)
            write(1,*) gl_x(yzel), gl_z(yzel)
        end do

        do j =  1, size(gl_RAY_O(1, :, kk))
            yzel = gl_RAY_O(par_n_BS - par_n_HP + 1, j, kk)
            write(1,*) gl_x(yzel), gl_z(yzel)
        end do


        ! Connectivity list
        do j = 1, N - 1
            write(1,*) j, j + 1
        end do

        do j = N + 1, 2 * N - 1
            write(1,*) j, j + 1
        end do


        do j = 2 * N + 1, 2 * N + 1 + M - 2
            write(1,*) j, j + 1
        end do

        do j = 2 * N + 1 + M, 2 * N + 1 + 2 * M - 2
            write(1,*) j, j + 1
        end do

        do j = 2 * N + 1 + 2 * M, 2 * N + 1 + 2 * M + L - 2
            write(1,*) j, j + 1
        end do

        do j = 2 * N + 1 + 2 * M + L, 2 * N + 1 + 2 * M + L + L - 2
            write(1,*) j, j + 1
        end do

        close(1)

	end subroutine Print_surface_y_2D
	
	subroutine Print_par_1D()
        use GEO_PARAM
        use STORAGE
        use Interpolate2
        implicit none

        integer :: N1, N2, kk, k, i, j, N, m, num
        real(8) :: c(3), Mach
        real(8) :: PAR(9)     ! Выходные параметры
        real(8) :: PAR_MOMENT(par_n_moment, par_n_sort)
        real(8) :: Dr, DI1, DI2, DI3
        
        print*, "Print_par_1D()"
        Dr = 1.0! 1.0/par_R0
        DI1 = 0.0264
        DI2 = 0.007622
        DI3 = 0.00292921

        num = 3
        N = size(gl_Cell_A(:, 1, 1)) 

        open(1, file = 'print_par_1D.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = r, rho, Ur, u, v, w, p, p+BB, BB, Bx, By, Bz, n_He,"
        write(1,*) "Max, In, Iu, IT"


        do j = 1, N
            c = gl_Cell_center(:, gl_Cell_A(j, 1, 1))
            m = gl_Cell_A(j, 1, 1)
            
            if(allocated(int2_Moment)) then
                call Int2_Get_par_fast(c(1), c(2), c(3), num, PAR, PAR_MOMENT = PAR_MOMENT)
            end if
            
            write(1,*) norm2(c) * Dr, gl_Cell_par(1, m ), DOT_PRODUCT(c, gl_Cell_par(2:4, m ))/norm2(c),&
                gl_Cell_par(2, m ), gl_Cell_par(3, m ), gl_Cell_par(4, m ), &
                gl_Cell_par(5, m ), gl_Cell_par(5, m ) + norm2(gl_Cell_par(6:8, m ))**2/(8.0 * par_pi_8), norm2(gl_Cell_par(6:8, m ))**2/(8.0 * par_pi_8),&
                gl_Cell_par(6, m ), gl_Cell_par(7, m ), gl_Cell_par(8, m ), &
                gl_Cell_par2(1, m ), &
                norm2(gl_Cell_par(2:4, m ))/sqrt(ggg * gl_Cell_par(5, m )/gl_Cell_par(1, m )),&
                sum(PAR_MOMENT(19, :)) * DI1, sum(PAR_MOMENT(6, :)) * DI2, sum(PAR_MOMENT(9, :)) * DI3
        end do
        
        close(1)
        
        print*, "****   TS = ", (norm2(gl_Cell_center(:, gl_Cell_A(par_n_TS - 1, 1, 1))) + norm2(gl_Cell_center(:, gl_Cell_A(par_n_TS, 1, 1))))/2.0
        print*, "****   HP = ", (norm2(gl_Cell_center(:, gl_Cell_A(par_n_HP - 1, 1, 1))) + norm2(gl_Cell_center(:, gl_Cell_A(par_n_HP, 1, 1))))/2.0
	
	end subroutine Print_par_1D

    subroutine Print_par_1D_PUI()
        use GEO_PARAM
        use STORAGE
        use Interpolate2
        use PUI
        implicit none

        integer :: N1, N2, kk, k, i, j, N, m, num
        real(8) :: c(3), Mach
        real(8) :: PAR(9)     ! Выходные параметры
        real(8) :: PAR_MOMENT(par_n_moment, par_n_sort)
        real(8) :: MAS_PUI(2)
        real(8) :: Dr, DI1, DI2, DI3
        real(8) :: Tpui, rho, rhoHe, rhopui, Tp, pp, p_th, rho_th, T_th
        
        print*, "Start Print_par_1D_PUI"
        Dr = 1.0/par_R0
        DI1 = 0.0264
        DI2 = 0.007622
        DI3 = 0.00292921

        num = 3
        N = size(gl_Cell_A(:, 1, 1)) 

        open(1, file = 'print_par_1D_PUI.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = r, rho, Ur, u, v, w, p, p+BB, BB, Bx, By, Bz, n_He,"
        write(1,*) "Max, In, Iu, IT, T, n_pui, T_pui, p_pui, T_th, rho_th, p_th"


        do j = 1, N
            c = gl_Cell_center(:, gl_Cell_A(j, 1, 1))
            m = gl_Cell_A(j, 1, 1)
            
            MAS_PUI = 0.0
            if(allocated(int2_Moment)) then
                call Int2_Get_par_fast(c(1), c(2), c(3), num, PAR, PAR_MOMENT = PAR_MOMENT, MAS_PUI = MAS_PUI)
            end if

            rho = gl_Cell_par(1, m)
            rhopui = MAS_PUI(1)
            rhoHe = gl_Cell_par2(1, m )
            Tpui = MAS_PUI(2)
            pp = gl_Cell_par(5, m )

            rho_th = rho - rhopui - rhoHe

            p_th = 4.0 * rho_th * (pp - rhopui * Tpui)/(8.0 * rho - 5.0 * rhoHe - 4.0 * rhopui)
            T_th = p_th/rho_th

            
            write(1,*) norm2(c) * Dr, gl_Cell_par(1, m ), DOT_PRODUCT(c, gl_Cell_par(2:4, m ))/norm2(c),&
                gl_Cell_par(2, m ), gl_Cell_par(3, m ), gl_Cell_par(4, m ), &
                gl_Cell_par(5, m ), gl_Cell_par(5, m ) + norm2(gl_Cell_par(6:8, m ))**2/(8.0 * par_pi_8), norm2(gl_Cell_par(6:8, m ))**2/(8.0 * par_pi_8),&
                gl_Cell_par(6, m ), gl_Cell_par(7, m ), gl_Cell_par(8, m ), &
                gl_Cell_par2(1, m ), &
                norm2(gl_Cell_par(2:4, m ))/sqrt(ggg * gl_Cell_par(5, m )/gl_Cell_par(1, m )),&
                sum(PAR_MOMENT(19, :)) * DI1, sum(PAR_MOMENT(6, :)) * DI2, sum(PAR_MOMENT(9, :)) * DI3, &
                gl_Cell_par(5, m )/(2.0 * gl_Cell_par(1, m )), MAS_PUI(1), MAS_PUI(2), MAS_PUI(1) * MAS_PUI(2), T_th, rho_th, p_th
        end do
        
        close(1)
        
	end subroutine Print_par_1D_PUI
	
	subroutine Print_tok_layer()
	    use GEO_PARAM
        use STORAGE
		USE Interpolate2
        implicit none
		
		real(8) :: n1, n2, n3, phi0, the0, t, phi, the
		real(8) :: x, y, z
		real(8) :: Matr(3, 3), Matr2(3, 3), cc(3), a(3), b(3)
		integer :: num, i, j, ii, k, NN1, NN2, jj, m
		real(8) :: PAR(9), dt     ! Выходные параметры
		real(8), allocatable :: SLOY(:, :, :)     ! Выходные параметры
		
		
		NN1 = 180
		NN2 = 20000
		allocate(SLOY(3, NN1, NN2))

		SLOY = 0.0
	    num = 3
		
	    Matr(1,1) = -0.995868
        Matr(1,2) = 0.0177307
        Matr(1,3) = 0.0890681
        Matr(2,1) = 0.0730412
        Matr(2,2) = 0.739193
        Matr(2,3) = 0.669521
        Matr(3,1) = -0.0539675
        Matr(3,2) = 0.67326
        Matr(3,3) = -0.737434
	
		
		
		t = 0.0
		dt = 0.0001
		
		do i = 1, NN2
		    if (mod(i, 100) == 0) print*, "i = ", i, " from", NN2
			t =	t + dt
			
			
			
			phi0 = cos(t * 2.0 * par_pi_8/0.013) * 2.0 * par_pi_8
		    the0 = cos(t * 2.0 * par_pi_8/9.4837) * par_pi_8/30.0
		
		    n1 = par_R0 * sin(the0) * cos(phi0)
		    n2 = par_R0 * sin(the0) * sin(phi0)
		    n3 = par_R0 * cos(the0)
			
			if (i <= NN2) then
			    do j = 1, NN1
		            phi = (j) * 2.0 * par_pi_8/(NN1 - 1)
		            the = par_pi_8/2.0 - atan((n1 * cos(phi) + n2 * sin(phi))/(-n3))
					!print*, i, j, phi/par_pi_8 * 180, the/par_pi_8 * 180
		
		            x = par_R0 * sin(the) * cos(phi)
		            y = par_R0 * sin(the) * sin(phi)
		            z = par_R0 * cos(the)
	
		            cc = MATMUL(Matr, (/ x, y, z /) )
				    SLOY(:, j, i) = cc
			    end do
			end if
			
			do ii = 1, min(NN2, i)
				do j = 1, NN1
					call Int2_Get_par_fast(SLOY(1, j, ii), SLOY(2, j, ii), SLOY(3, j, ii), num, PAR)
					SLOY(:, j, ii) = SLOY(:, j, ii) + dt * PAR(2:4)
				end do
			end do
		
		end do
		
		
		open(1, file = 'Print_Tokoviy_sloy.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z'  ZONE T= 'HP', N= ", (NN2 - 1) * (NN1) * 4  , ", E =  ", (NN2 - 1) * (NN1) , ", F=FEPOINT, ET=quadrilateral "
	    
		do i = 1, NN2 - 1
			do j = 1, NN1
				jj = j + 1
				if(jj > NN1) jj = 1
				write(1,*) SLOY(:, j, i)
				write(1,*) SLOY(:, j, i + 1)
				write(1,*) SLOY(:, jj, i + 1)
				write(1,*) SLOY(:, jj, i)
			end do
		end do
		
		do k = 1, (NN2 - 1) * (NN1)
		    write(1,*) 4 * k - 3, 4 * k - 2, 4 * k - 1, 4 * k
		end do
		
		close(1)
		
		open(2, file = 'Print_Tokoviy_sloy_mini.txt')
        write(2,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z'  ZONE T= 'HP', N= ", (NN2 - 1) * (NN1 - 1) * 4  , ", E =  ", (NN2 - 1) * (NN1 - 1) , ", F=FEPOINT, ET=quadrilateral "
	    
		do i = 1, NN2 - 1
			do j = 1, NN1 - 1
				jj = j + 1
				write(2,*) SLOY(:, j, i)
				write(2,*) SLOY(:, j, i + 1)
				write(2,*) SLOY(:, jj, i + 1)
				write(2,*) SLOY(:, jj, i)
			end do
		end do
		
		do k = 1, (NN2 - 1) * (NN1 - 1)
		    write(2,*) 4 * k - 3, 4 * k - 2, 4 * k - 1, 4 * k
		end do
		close(2)
		
		open(3, file = '2D_Print_Tokoviy_sloy.txt')
        write(3,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y'  ZONE T= 'HP', N= ", 0 , ", E =  ", 0 , ", F=FEPOINT, ET=LINESEG "
	    m = 0
		do i = 1, NN2 - 1
			do j = 1, NN1
				jj = j + 1
				if(jj > NN1) jj = 1
				
				a = SLOY(:, j, i)
				b = SLOY(:, j, i + 1)
				t = -a(3)/(b(3) - a(3))
				if(t > 0.0 .and. t < 1.0) then
					write(3,*) a(1) + (b(1) - a(1)) * t, a(2) + (b(2) - a(2)) * t
					m = m + 1
				end if
				
				a = SLOY(:, j, i + 1)
				b = SLOY(:, jj, i + 1)
				t = -a(3)/(b(3) - a(3))
				if(t > 0.0 .and. t < 1.0) then
					write(3,*) a(1) + (b(1) - a(1)) * t, a(2) + (b(2) - a(2)) * t
					m = m + 1
				end if
				
				a = SLOY(:, jj, i + 1)
				b = SLOY(:, jj, i)
				t = -a(3)/(b(3) - a(3))
				if(t > 0.0 .and. t < 1.0) then
					write(3,*) a(1) + (b(1) - a(1)) * t, a(2) + (b(2) - a(2)) * t
					m = m + 1
				end if
				
				a = SLOY(:, jj, i)
				b = SLOY(:, j, i)
				t = -a(3)/(b(3) - a(3))
				if(t > 0.0 .and. t < 1.0) then
					write(3,*) a(1) + (b(1) - a(1)) * t, a(2) + (b(2) - a(2)) * t
					m = m + 1
				end if
				
			end do
		end do
		
		print*, "N = ", m
		print*, "E = ", m/2
		
		do k = 1, m/2
		    write(3,*) 2 * k - 1, 2 * k
		end do
		
		close(3)
		
		deallocate(SLOY)
	
	end subroutine Print_tok_layer

    subroutine Cucl_div_V()
        USE STORAGE
        ! Посчитаем дивергенцию скорости для ИГОРЯ
        integer(4) :: i, N, j, gr, sosed
        real(8) :: V1(3), V2(3), DIV, d1, d2, SS, normal(3)

        N = size(gl_all_Cell(1, :))

        do i = 1, N
            V = 0.0
            V1 = gl_Cell_par(2:4, i)
            do j = 1, 6
                gr = gl_Cell_gran(j, i)
                if(gr == 0) CYCLE

                d1 = norm2(gl_Gran_center(:, gr) - gl_Cell_center(:, i))

                sosed = gl_Gran_neighbour(1, gr)
                if(sosed == i) sosed = gl_Gran_neighbour(2, gr)

                SS = gl_Gran_square(gr)

                normal = gl_Gran_normal(:, gr)

                if(gl_Gran_neighbour(1, gr) /= i) normal = -normal

                if(sosed > 0) then
                    if(gl_zone_Cell(sosed) == gl_zone_Cell(i)) then
                        V2 = gl_Cell_par(2:4, sosed)
                        d2 = norm2(gl_Gran_center(:, gr) - gl_Cell_center(:, sosed))

                        V = V + DOT_PRODUCT(normal, (d2 * V1 + d1 * V2)/(d1 + d2)) * SS
                    else
                        V = V + DOT_PRODUCT(normal, V1) * SS
                    end if
                else
                    V = V + DOT_PRODUCT(normal, V1) * SS
                end if

            end do

            DIV = V/gl_Cell_Volume(i)
            gl_Cell_par_div(i) = DIV
            ! print*, "DIV = ", DIV
            ! pause
        end do
    end subroutine Cucl_div_V

    subroutine Send_request()
        ! Функция для передачи параметров в двумерную сетку

        use STORAGE
        use GEO_PARAM
        use Interpolate2
        implicit none
        integer :: n, i, num
        real(8) :: x, y
        logical :: exists
        real(8) :: PAR(9), PAR_MOMENT(par_n_moment, par_n_sort)
        
        
        num = 3
        inquire(file="request.bin", exist=exists)
        
        if (exists == .False.) then
            pause "net faila!!! f fwf234r435234r3rw3r4"
            STOP "net faila!!!"
        end if
        
        
        open(1, file = "request.bin", FORM = 'BINARY', ACTION = "READ")
        open(2, file = "answer_request.bin", FORM = 'BINARY')

        read(1) n
        write(2) n

        do i = 1, n
            read(1) x, y
            call Int2_Get_par_fast(x, y, 0.0001_8, num, PAR, PAR_MOMENT)
            write(2) PAR(1), PAR(2), PAR(3), PAR(5)
            write(2) PAR_MOMENT(1, 1), PAR_MOMENT(2, 1), PAR_MOMENT(3, 1), PAR_MOMENT(5, 1)
            write(2) PAR_MOMENT(1, 2), PAR_MOMENT(2, 2), PAR_MOMENT(3, 2), PAR_MOMENT(5, 2)
            write(2) PAR_MOMENT(1, 3), PAR_MOMENT(2, 3), PAR_MOMENT(3, 3), PAR_MOMENT(5, 3)
            write(2) PAR_MOMENT(1, 4), PAR_MOMENT(2, 4), PAR_MOMENT(3, 4), PAR_MOMENT(5, 4)
        end do

        close(1)
        close(2)

    end subroutine Send_request


    subroutine Print_par_2D()  ! Печатает 2Д сетку с линиями в Техплот
        use GEO_PARAM
        use STORAGE
        implicit none

        integer :: N1, N2, kk, k, i, j, N, m
        real(8) :: c(3), Mach

        print*, "Print_par_2D()"
        N = size(gl_Cell_A(1, :, 1)) * size(gl_Cell_A(:, 1, 1)) + &
            size(gl_Cell_B(1, :, 1)) * size(gl_Cell_B(:, 1, 1)) + &
            size(gl_Cell_C(1, :, 1)) * size(gl_Cell_C(:, 1, 1))
        N = N * 2

        open(1, file = 'print_par_2D.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z', 'rho', 'u', 'v', 'w', 'p',"
        write(1,*) "'bx', 'by', 'bz', 'bb', 'Volume', 'Mach', 'Q', 'He',"
        write(1,*) "'Zone','T','rho1', 'u1', 'v1', 'w1', 'p1', 'rho2',"
        write(1,*)" 'u2', 'v2', 'w2', 'p2', 'rho3', 'u3', 'v3', 'w3', 'p3', "
        write(1,*) "'rho4', 'u4', 'v4', 'w4', 'p4', ZONE T= 'HP'"


        kk = 1
        N2 = size(gl_Cell_A(1, :, 1))
        N1 = size(gl_Cell_A(:, 1, 1))
        do j = 1, N2
            do i = 1, N1
                c = gl_Cell_center(:, gl_Cell_A(i, j, kk))
                m = gl_Cell_A(i, j, kk)
                Mach = norm2(gl_Cell_par(2:4, m ))/sqrt(ggg*gl_Cell_par(5, m )/gl_Cell_par(1, m ))
                write(1,*) c, gl_Cell_par(1:8, m ), norm2(gl_Cell_par(6:8, m ))/(8*par_pi_8), gl_Cell_Volume(m), Mach, &
                    gl_Cell_par(9, m )/gl_Cell_par(1, m ), gl_Cell_par2(1, m), gl_zone_Cell(m), gl_Cell_par(5, m )/gl_Cell_par(1, m ),&
                    gl_Cell_par_MF(:, :, m)
            end do
        end do
        N2 = size(gl_Cell_B(1, :, 1))
        N1 = size(gl_Cell_B(:, 1, 1))
        do j = 1, N2
            do i = 1, N1
                m = gl_Cell_B(i, j, kk)
                c = gl_Cell_center(:, m)
                Mach = norm2(gl_Cell_par(2:4, m ))/sqrt(ggg*gl_Cell_par(5, m )/gl_Cell_par(1, m ))
                write(1,*)  c, gl_Cell_par(1:8, m ), norm2(gl_Cell_par(6:8, m ))/(8*par_pi_8), gl_Cell_Volume(m), Mach, &
                    gl_Cell_par(9, m )/gl_Cell_par(1, m ), gl_Cell_par2(1, m), gl_zone_Cell(m), gl_Cell_par(5, m )/gl_Cell_par(1, m ), &
                    gl_Cell_par_MF(:, :, m)
            end do
        end do
        N2 = size(gl_Cell_C(1, :, 1))
        N1 = size(gl_Cell_C(:, 1, 1))
        do j = 1, N2
            do i = 1, N1
                m = gl_Cell_C(i, j, kk)
                c = gl_Cell_center(:, m)
                Mach = norm2(gl_Cell_par(2:4, m ))/sqrt(ggg*gl_Cell_par(5, m )/gl_Cell_par(1, m ))
                write(1,*)  c, gl_Cell_par(1:8, m ), norm2(gl_Cell_par(6:8, m ))/(8*par_pi_8), gl_Cell_Volume(m), Mach, &
                    gl_Cell_par(9, m )/gl_Cell_par(1, m ), gl_Cell_par2(1, m), gl_zone_Cell(m), gl_Cell_par(5, m )/gl_Cell_par(1, m ),&
                    gl_Cell_par_MF(:, :, m)
            end do
        end do

        ! Нижняя часть сетки

        kk = par_l_phi/2 + 1
        N2 = size(gl_Cell_A(1, :, 1))
        N1 = size(gl_Cell_A(:, 1, 1))
        do j = 1, N2
            do i = 1, N1
                m = gl_Cell_A(i, j, kk)
                c = gl_Cell_center(:, m)
                Mach = norm2(gl_Cell_par(2:4, m ))/sqrt(ggg*gl_Cell_par(5, m )/gl_Cell_par(1, m ))
                write(1,*)  c, gl_Cell_par(1:8, m ), norm2(gl_Cell_par(6:8, m ))/(8*par_pi_8), gl_Cell_Volume(m), Mach, &
                    gl_Cell_par(9, m )/gl_Cell_par(1, m ), gl_Cell_par2(1, m), gl_zone_Cell(m), gl_Cell_par(5, m )/gl_Cell_par(1, m ), &
                    gl_Cell_par_MF(:, :, m)
            end do
        end do
        N2 = size(gl_Cell_B(1, :, 1))
        N1 = size(gl_Cell_B(:, 1, 1))
        do j = 1, N2
            do i = 1, N1
                m = gl_Cell_B(i, j, kk)
                c = gl_Cell_center(:, m)
                Mach = norm2(gl_Cell_par(2:4, m ))/sqrt(ggg*gl_Cell_par(5, m )/gl_Cell_par(1, m ))
                write(1,*)  c, gl_Cell_par(1:8, m ), norm2(gl_Cell_par(6:8, m ))/(8*par_pi_8), gl_Cell_Volume(m), Mach, &
                    gl_Cell_par(9, m )/gl_Cell_par(1, m ), gl_Cell_par2(1, m), gl_zone_Cell(m), gl_Cell_par(5, m )/gl_Cell_par(1, m ), &
                    gl_Cell_par_MF(:, :, m)
            end do
        end do
        N2 = size(gl_Cell_C(1, :, 1))
        N1 = size(gl_Cell_C(:, 1, 1))
        do j = 1, N2
            do i = 1, N1
                m = gl_Cell_C(i, j, kk)
                c = gl_Cell_center(:, m)
                Mach = norm2(gl_Cell_par(2:4, m ))/sqrt(ggg*gl_Cell_par(5, m )/gl_Cell_par(1, m ))
                write(1,*)  c, gl_Cell_par(1:8, m ), norm2(gl_Cell_par(6:8, m ))/(8*par_pi_8), gl_Cell_Volume(m), Mach, &
                    gl_Cell_par(9, m )/gl_Cell_par(1, m ), gl_Cell_par2(1, m), gl_zone_Cell(m), gl_Cell_par(5, m )/gl_Cell_par(1, m ), &
                    gl_Cell_par_MF(:, :, m)
            end do
        end do


        close(1)
	end subroutine Print_par_2D
	
	subroutine Print_par_y_2D()  ! Печатает 2Д сетку с линиями в Техплот
        use GEO_PARAM
        use STORAGE
        implicit none

        integer :: N1, N2, kk, k, i, j, N, m
        real(8) :: c(3), Mach


        N = size(gl_Cell_A(1, :, 1)) * size(gl_Cell_A(:, 1, 1)) + &
            size(gl_Cell_B(1, :, 1)) * size(gl_Cell_B(:, 1, 1)) + &
            size(gl_Cell_C(1, :, 1)) * size(gl_Cell_C(:, 1, 1))
        N = N * 2

        open(1, file = 'print_par_y_2D.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Z', 'Y', 'rho', 'u', 'v', 'w', 'p',"
        write(1,*) "'bx', 'by', 'bz', 'bb', 'Volume', 'Mach', 'Q',"
        write(1,*) "'Zone', 'T','rho1', 'u1', 'v1', 'w1', 'p1', 'rho2',"
        write(1,*)" 'u2', 'v2', 'w2', 'p2', 'rho3', 'u3', 'v3', 'w3', 'p3', "
        write(1,*) "'rho4', 'u4', 'v4', 'w4', 'p4', ZONE T= 'HP'"


        kk = par_l_phi/4
        N2 = size(gl_Cell_A(1, :, 1))
        N1 = size(gl_Cell_A(:, 1, 1))
        do j = 1, N2
            do i = 1, N1
                c = gl_Cell_center(:, gl_Cell_A(i, j, kk))
                m = gl_Cell_A(i, j, kk)
                Mach = norm2(gl_Cell_par(2:4, m ))/sqrt(ggg*gl_Cell_par(5, m )/gl_Cell_par(1, m ))
                write(1,*) c(1), c(3), c(2), gl_Cell_par(1:8, m ), norm2(gl_Cell_par(6:8, m )), gl_Cell_Volume(m), Mach, gl_Cell_par(9, m )/gl_Cell_par(1, m ), gl_zone_Cell(m), gl_Cell_par(5, m )/gl_Cell_par(1, m ), gl_Cell_par_MF(:, :, m)
            end do
        end do
        N2 = size(gl_Cell_B(1, :, 1))
        N1 = size(gl_Cell_B(:, 1, 1))
        do j = 1, N2
            do i = 1, N1
                m = gl_Cell_B(i, j, kk)
                c = gl_Cell_center(:, m)
                Mach = norm2(gl_Cell_par(2:4, m ))/sqrt(ggg*gl_Cell_par(5, m )/gl_Cell_par(1, m ))
                write(1,*)  c(1), c(3), c(2), gl_Cell_par(1:8, m ), norm2(gl_Cell_par(6:8, m )), gl_Cell_Volume(m), Mach, gl_Cell_par(9, m )/gl_Cell_par(1, m ), gl_zone_Cell(m), gl_Cell_par(5, m )/gl_Cell_par(1, m ), gl_Cell_par_MF(:, :, m)
            end do
        end do
        N2 = size(gl_Cell_C(1, :, 1))
        N1 = size(gl_Cell_C(:, 1, 1))
        do j = 1, N2
            do i = 1, N1
                m = gl_Cell_C(i, j, kk)
                c = gl_Cell_center(:, m)
                Mach = norm2(gl_Cell_par(2:4, m ))/sqrt(ggg*gl_Cell_par(5, m )/gl_Cell_par(1, m ))
                write(1,*)  c(1), c(3), c(2), gl_Cell_par(1:8, m ), norm2(gl_Cell_par(6:8, m )), gl_Cell_Volume(m), Mach, gl_Cell_par(9, m )/gl_Cell_par(1, m ), gl_zone_Cell(m), gl_Cell_par(5, m )/gl_Cell_par(1, m ), gl_Cell_par_MF(:, :, m)
            end do
        end do

        ! Нижняя часть сетки

        kk = 3 * par_l_phi/4 + 1
        N2 = size(gl_Cell_A(1, :, 1))
        N1 = size(gl_Cell_A(:, 1, 1))
        do j = 1, N2
            do i = 1, N1
                m = gl_Cell_A(i, j, kk)
                c = gl_Cell_center(:, m)
                Mach = norm2(gl_Cell_par(2:4, m ))/sqrt(ggg*gl_Cell_par(5, m )/gl_Cell_par(1, m ))
                write(1,*)  c(1), c(3), c(2), gl_Cell_par(1:8, m ), norm2(gl_Cell_par(6:8, m )), gl_Cell_Volume(m), Mach, gl_Cell_par(9, m )/gl_Cell_par(1, m ), gl_zone_Cell(m), gl_Cell_par(5, m )/gl_Cell_par(1, m ), gl_Cell_par_MF(:, :, m)
            end do
        end do
        N2 = size(gl_Cell_B(1, :, 1))
        N1 = size(gl_Cell_B(:, 1, 1))
        do j = 1, N2
            do i = 1, N1
                m = gl_Cell_B(i, j, kk)
                c = gl_Cell_center(:, m)
                Mach = norm2(gl_Cell_par(2:4, m ))/sqrt(ggg*gl_Cell_par(5, m )/gl_Cell_par(1, m ))
                write(1,*)  c(1), c(3), c(2), gl_Cell_par(1:8, m ), norm2(gl_Cell_par(6:8, m )), gl_Cell_Volume(m), Mach, gl_Cell_par(9, m )/gl_Cell_par(1, m ), gl_zone_Cell(m), gl_Cell_par(5, m )/gl_Cell_par(1, m ), gl_Cell_par_MF(:, :, m)
            end do
        end do
        N2 = size(gl_Cell_C(1, :, 1))
        N1 = size(gl_Cell_C(:, 1, 1))
        do j = 1, N2
            do i = 1, N1
                m = gl_Cell_C(i, j, kk)
                c = gl_Cell_center(:, m)
                Mach = norm2(gl_Cell_par(2:4, m ))/sqrt(ggg*gl_Cell_par(5, m )/gl_Cell_par(1, m ))
                write(1,*)  c(1), c(3), c(2), gl_Cell_par(1:8, m ), norm2(gl_Cell_par(6:8, m )), gl_Cell_Volume(m), Mach, gl_Cell_par(9, m )/gl_Cell_par(1, m ), gl_zone_Cell(m), gl_Cell_par(5, m )/gl_Cell_par(1, m ), gl_Cell_par_MF(:, :, m)
            end do
        end do


        close(1)
    end subroutine Print_par_y_2D

    subroutine Print_Point_Plane()
        ! Точки в плоскости XOY - имеют две координаты
        ! Функция для проверки геометрии
        use GEO_PARAM
        use STORAGE
        implicit none

        integer :: N1, N2, k, l, i, j

        open(1, file = 'all_point_plane.txt')

        N2 = size(gl_RAY_A(1, :, 1))
        N1 = size(gl_RAY_A(:, 1, 1))

        k = 1
        do j = 1, N2
            do i = 1, N1
                l = gl_RAY_A(i, j, k)
                write(1,*) gl_x(l), gl_y(l)
            end do
        end do

        k = par_l_phi/2 + 1
        do j = 1, N2
            do i = 1, N1
                l = gl_RAY_A(i, j, k)
                write(1,*) gl_x(l), gl_y(l)
            end do
        end do

        ! ----------------------------------------------------- B
        N2 = size(gl_RAY_B(1, :, 1))
        N1 = size(gl_RAY_B(:, 1, 1))

        k = 1
        do j = 1, N2
            do i = 1, N1
                l = gl_RAY_B(i, j, k)
                write(1,*) gl_x(l), gl_y(l)
            end do
        end do

        k = par_l_phi/2 + 1
        do j = 1, N2
            do i = 1, N1
                l = gl_RAY_B(i, j, k)
                write(1,*) gl_x(l), gl_y(l)
            end do
        end do

        ! ----------------------------------------------------- C
        N2 = size(gl_RAY_C(1, :, 1))
        N1 = size(gl_RAY_C(:, 1, 1))

        k = 1
        do j = 1, N2
            do i = 1, N1
                l = gl_RAY_C(i, j, k)
                write(1,*) gl_x(l), gl_y(l)
            end do
        end do

        k = par_l_phi/2 + 1
        do j = 1, N2
            do i = 1, N1
                l = gl_RAY_C(i, j, k)
                write(1,*) gl_x(l), gl_y(l)
            end do
        end do

        ! ----------------------------------------------------- O
        N2 = size(gl_RAY_O(1, :, 1))
        N1 = size(gl_RAY_O(:, 1, 1))

        k = 1
        do j = 1, N2
            do i = 1, N1
                l = gl_RAY_O(i, j, k)
                write(1,*) gl_x(l), gl_y(l)
            end do
        end do

        k = par_l_phi/2 + 1
        do j = 1, N2
            do i = 1, N1
                l = gl_RAY_O(i, j, k)
                write(1,*) gl_x(l), gl_y(l)
            end do
        end do

        ! ----------------------------------------------------- K
        N2 = size(gl_RAY_K(1, :, 1))
        N1 = size(gl_RAY_K(:, 1, 1))

        k = 1
        do j = 1, N2
            do i = 1, N1
                l = gl_RAY_K(i, j, k)
                write(1,*) gl_x(l), gl_y(l)
            end do
        end do

        k = par_l_phi/2 + 1
        do j = 1, N2
            do i = 1, N1
                l = gl_RAY_K(i, j, k)
                write(1,*) gl_x(l), gl_y(l)
            end do
        end do

        ! ----------------------------------------------------- D
        N2 = size(gl_RAY_D(1, :, 1))
        N1 = size(gl_RAY_D(:, 1, 1))

        k = 1
        do j = 1, N2
            do i = 1, N1
                l = gl_RAY_D(i, j, k)
                write(1,*) gl_x(l), gl_y(l)
            end do
        end do

        k = par_l_phi/2 + 1
        do j = 1, N2
            do i = 1, N1
                l = gl_RAY_D(i, j, k)
                write(1,*) gl_x(l), gl_y(l)
            end do
        end do

        ! ----------------------------------------------------- E
        N2 = size(gl_RAY_E(1, :, 1))
        N1 = size(gl_RAY_E(:, 1, 1))

        k = 1
        do j = 1, N2
            do i = 1, N1
                l = gl_RAY_E(i, j, k)
                write(1,*) gl_x(l), gl_y(l)
            end do
        end do

        k = par_l_phi/2 + 1
        do j = 1, N2
            do i = 1, N1
                l = gl_RAY_E(i, j, k)
                write(1,*) gl_x(l), gl_y(l)
            end do
        end do

        close(1)

    end subroutine Print_Point_Plane

    subroutine Print_All_Points()   ! Функция просто печатает все точки сетки (далее можно просмотреть их в текплоте)
        ! Variables
        use STORAGE

        open(1, file = 'all_point.txt')

        do i=1,size(gl_x)
            write(1,*) gl_x(i), gl_y(i), gl_z(i)
        end do

        close(1)

    end subroutine Print_All_Points

    subroutine Print_A_Ray(n1, n2)   ! Функция просто печатает один луч (в плоскости, в пространстве) + красит ключевые точки

        use GEO_PARAM
        use STORAGE

        implicit none
        integer, intent(in) :: n1, n2
        integer(4) :: i, k

        open(1, file = 'Ray.txt')

        do i = 1, size(gl_RAY_A(:, n1, n2))
            k = 0
            if (i == par_n_TS .or. i == par_n_HP .or. i == par_n_BS .or. i == par_n_END) k = 1
            write(1,*) gl_x(gl_RAY_A(i, n1, n2)), gl_y(gl_RAY_A(i, n1, n2)), gl_z(gl_RAY_A(i, n1, n2)),  k
        end do

        close(1)

    end subroutine Print_A_Ray

    subroutine Print_Cell(n1, n2, n3, str)  ! Печатает ячейку с номерами n1, n2, n3 из множества ячеек str
        use GEO_PARAM
        use STORAGE

        implicit none
        integer, intent(in) :: n1, n2, n3
        character, intent(in) :: str
        integer(4) :: i, k

        open(1, file = 'one_Cell.txt')


        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z'  ZONE T= 'HP', N= ", 8, ", E =  ", 12 , ", F=FEPOINT, ET=LINESEG "

        do i = 1, 8
            if (str == "A") then
                if (n1 > size(gl_Cell_A(:, 1, 1))) PAUSE "n1 too mach! ERROR 822"
                if (n2 > size(gl_Cell_A(1, :, 1))) PAUSE "n2 too mach! ERROR 822"
                if (n3 > size(gl_Cell_A(1, 1, :))) PAUSE "n3 too mach! ERROR 822"
                write(1,*) gl_x(gl_all_Cell(i, gl_Cell_A(n1, n2, n3))), gl_y(gl_all_Cell(i, gl_Cell_A(n1, n2, n3))), gl_z(gl_all_Cell(i, gl_Cell_A(n1, n2, n3)))
            else if (str == "B") then
                if (n1 > size(gl_Cell_B(:, 1, 1))) PAUSE "n1 too mach! ERROR 827"
                if (n2 > size(gl_Cell_B(1, :, 1))) PAUSE "n2 too mach! ERROR 827"
                if (n3 > size(gl_Cell_B(1, 1, :))) PAUSE "n3 too mach! ERROR 827"
                write(1,*) gl_x(gl_all_Cell(i, gl_Cell_B(n1, n2, n3))), gl_y(gl_all_Cell(i, gl_Cell_B(n1, n2, n3))), gl_z(gl_all_Cell(i, gl_Cell_B(n1, n2, n3)))
            else
                if (n1 > size(gl_Cell_C(:, 1, 1))) PAUSE "n1 too mach! ERROR 829"
                if (n2 > size(gl_Cell_C(1, :, 1))) PAUSE "n2 too mach! ERROR 829"
                if (n3 > size(gl_Cell_C(1, 1, :))) PAUSE "n3 too mach! ERROR 829"
                write(1,*) gl_x(gl_all_Cell(i, gl_Cell_C(n1, n2, n3))), gl_y(gl_all_Cell(i, gl_Cell_C(n1, n2, n3))), gl_z(gl_all_Cell(i, gl_Cell_C(n1, n2, n3)))
            end if
        end do

        write(1,*) "1 2"
        write(1,*) "2 3"
        write(1,*) "3 4"
        write(1,*) "4 1"
        write(1,*) "1 5"
        write(1,*) "2 6"
        write(1,*) "3 7"
        write(1,*) "4 8"
        write(1,*) "5 6"
        write(1,*) "6 7"
        write(1,*) "7 8"
        write(1,*) "8 5"

        close(1)

    end subroutine Print_Cell

    subroutine Print_all_Cell()  ! Печатает все ячейки сетки (с линями) - Либо ячейки с каким-то условием
        use GEO_PARAM
        use STORAGE
        implicit none

        integer(4) :: N1, i, j

        N1 = size(gl_Cell_A(:, :, :)) + size(gl_Cell_B(:, :, :)) + size(gl_Cell_C(:, :, :))

        open(1, file = 'all_Cell.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z'  ZONE T= 'HP', N= ", 8 * N1, ", E =  ", 12 * N1 , ", F=FEPOINT, ET=LINESEG "



        do j = 1, N1
            do i = 1, 8
                write(1,*) gl_x(gl_all_Cell(i,j)), gl_y(gl_all_Cell(i,j)), gl_z(gl_all_Cell(i,j))
            end do
        end do

        do j = 0, N1-1
            write(1,*) 8 * j + 1, 8 * j + 2
            write(1,*) 8 * j + 2, 8 * j + 3
            write(1,*) 8 * j + 3, 8 * j + 4
            write(1,*) 8 * j + 4, 8 * j + 1
            write(1,*) 8 * j + 1, 8 * j + 5
            write(1,*) 8 * j + 2, 8 * j + 6
            write(1,*) 8 * j + 3, 8 * j + 7
            write(1,*) 8 * j + 4, 8 * j + 8
            write(1,*) 8 * j + 5, 8 * j + 6
            write(1,*) 8 * j + 6, 8 * j + 7
            write(1,*) 8 * j + 7, 8 * j + 8
            write(1,*) 8 * j + 8, 8 * j + 5
        end do


        close(1)


    end subroutine Print_all_Cell

    subroutine Print_cell_and_neighbour(n1, n2, n3)  ! Печатает всех соседей ячейки и саму ячейку
        use GEO_PARAM
        use STORAGE
        implicit none

        integer, intent(in) :: n1, n2, n3
        integer :: n, m, i, j, k

        n = gl_Cell_C(n1, n2, n3)  ! Получили номер ячейки в общем списке

        m = COUNT(gl_Cell_neighbour(:, n) > 0)

        if (par_developer_info) print *, "Kolichestvo sosedey = ", m

        open(1, file = 'Cell+neighbour.txt')


        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z', 'M'  ZONE T= 'HP', N= ", 8 * (m + 1), ", E =  ", 12 * (m + 1) , ", F=FEPOINT, ET=LINESEG "

        ! Печатаем саму ячейку
        do i = 1, 8
            write(1,*) gl_x(gl_all_Cell(i, n)), gl_y(gl_all_Cell(i, n)), gl_z(gl_all_Cell(i, n)), 0
        end do

        ! Печатаем соседей
        do j = 1, 6
            if (gl_Cell_neighbour(j, n) <= 0) CYCLE
            k = gl_Cell_neighbour(j, n)
            do i = 1, 8
                write(1,*) gl_x(gl_all_Cell(i, k)), gl_y(gl_all_Cell(i, k)), gl_z(gl_all_Cell(i, k)), j
            end do
        end do


        do j = 0, m
            write(1,*) 8 * j + 1, 8 * j + 2
            write(1,*) 8 * j + 2, 8 * j + 3
            write(1,*) 8 * j + 3, 8 * j + 4
            write(1,*) 8 * j + 4, 8 * j + 1
            write(1,*) 8 * j + 1, 8 * j + 5
            write(1,*) 8 * j + 2, 8 * j + 6
            write(1,*) 8 * j + 3, 8 * j + 7
            write(1,*) 8 * j + 4, 8 * j + 8
            write(1,*) 8 * j + 5, 8 * j + 6
            write(1,*) 8 * j + 6, 8 * j + 7
            write(1,*) 8 * j + 7, 8 * j + 8
            write(1,*) 8 * j + 8, 8 * j + 5
        end do

        close(1)

    end subroutine Print_cell_and_neighbour

    subroutine Print_gran()  ! Печатает грани сетки (можно с каким-то условием)
        use GEO_PARAM
        use STORAGE
        implicit none

        integer :: n, m, i, j, k

        m = 0

        !do i = 1, size(gl_all_Gran(1,:))
        !    if ( ANY( gl_Gran_neighbour(:,i) == -2 ) ) m = m + 1   ! Это условие нужно задавать и ниже (в самом цикле печати)
        !end do
        m = 1

        if (par_developer_info) print *, "Kolichestvo graney = ", m

        open(1, file = 'Grans.txt')


        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z', 'M'  ZONE T= 'HP', N= ", 4 * (m), ", E =  ", 4 * (m) , ", F=FEPOINT, ET=LINESEG "

        ! Печатаем
        do n = 3953, 3953 !1, size(gl_all_Gran(1,:))
            !if ( ANY( gl_Gran_neighbour(:,n) == -2 ) == .FALSE. ) CYCLE
            do i = 1, 4
                write(1,*) gl_x(gl_all_Gran(i, n)), gl_y(gl_all_Gran(i, n)), gl_z(gl_all_Gran(i, n)), n
            end do
        end do


        do j = 0, m-1
            write(1,*) 4 * j + 1, 4 * j + 2
            write(1,*) 4 * j + 2, 4 * j + 3
            write(1,*) 4 * j + 3, 4 * j + 4
            write(1,*) 4 * j + 4, 4 * j + 1
        end do

        close(1)

    end subroutine Print_gran

    subroutine Geometry_check()
        use GEO_PARAM
        use STORAGE

        implicit none
        integer(4) :: N1, N2, N3, i, j, k, j1
        integer(4) :: a1, a2, a3
        real(8) :: c(3), SS

        if (par_developer_info) print *, "Geometry_check: Proverka 1"
        ! Проверяем каждому ли лучу задали все точки (не осталось ли неинициализированных)
        ! ----------------------------------------------------- A
        N2 = size(gl_RAY_A(1, :, 1))
        N1 = size(gl_RAY_A(:, 1, 1))
        N3 = size(gl_RAY_A(1, 1, :))
        do k = 1, N3
            do j = 1, N2
                do i = 1, N1
                    if (gl_RAY_A(i, j, k) < 1) STOP "Error gl_RAY_A 618"
                end do
            end do
        end do

        ! ----------------------------------------------------- B
        N2 = size(gl_RAY_B(1, :, 1))
        N1 = size(gl_RAY_B(:, 1, 1))
        N3 = size(gl_RAY_B(1, 1, :))
        do k = 1, N3
            do j = 1, N2
                do i = 1, N1
                    if (gl_RAY_B(i, j, k) < 1) STOP "Error gl_RAY_B 618"
                end do
            end do
        end do

        ! ----------------------------------------------------- C
        N2 = size(gl_RAY_C(1, :, 1))
        N1 = size(gl_RAY_C(:, 1, 1))
        N3 = size(gl_RAY_C(1, 1, :))
        do k = 1, N3
            do j = 1, N2
                do i = 1, N1
                    if (gl_RAY_C(i, j, k) < 1) STOP "Error gl_RAY_C 618"
                end do
            end do
        end do

        ! ----------------------------------------------------- O
        N2 = size(gl_RAY_O(1, :, 1))
        N1 = size(gl_RAY_O(:, 1, 1))
        N3 = size(gl_RAY_O(1, 1, :))
        do k = 1, N3
            do j = 1, N2
                do i = 1, N1
                    if (gl_RAY_O(i, j, k) < 1) STOP "Error gl_RAY_O 618"
                end do
            end do
        end do

        ! ----------------------------------------------------- K
        N2 = size(gl_RAY_K(1, :, 1))
        N1 = size(gl_RAY_K(:, 1, 1))
        N3 = size(gl_RAY_K(1, 1, :))
        do k = 1, N3
            do j = 1, N2
                do i = 1, N1
                    if (gl_RAY_K(i, j, k) < 1) STOP "Error gl_RAY_K 618"
                end do
            end do
        end do

        ! ----------------------------------------------------- D
        N2 = size(gl_RAY_D(1, :, 1))
        N1 = size(gl_RAY_D(:, 1, 1))
        N3 = size(gl_RAY_D(1, 1, :))
        do k = 1, N3
            do j = 1, N2
                do i = 1, N1
                    if (gl_RAY_D(i, j, k) < 1) STOP "Error gl_RAY_D 618"
                end do
            end do
        end do

        ! ----------------------------------------------------- E
        N2 = size(gl_RAY_E(1, :, 1))
        N1 = size(gl_RAY_E(:, 1, 1))
        N3 = size(gl_RAY_E(1, 1, :))
        do k = 1, N3
            do j = 1, N2
                do i = 1, N1
                    if (gl_RAY_E(i, j, k) < 1) STOP "Error gl_RAY_E 618"
                end do
            end do
        end do

        if (par_developer_info) print *, "Geometry_check: Proverka 2"
        ! Проверяем правильность соседей (чтобы была обоюдность)
        do i = 1, size(gl_all_Cell(1, :))
            do j = 1, 6
                if (gl_Cell_neighbour(j, i) > 0) then
                    j1 = gl_Cell_neighbour(j, i)
                    if( ANY(gl_Cell_neighbour(:, j1) == i) == .FALSE.) then
                        pause "ERROR 1763"
                    end if
                end if
            end do
        end do

        if (par_developer_info) print *, "Geometry_check: Proverka 3"
        ! Проверяем правильность соседей граней (чтобы они реально связывали соседей, а не абы кого)
        do i = 1, size(gl_Gran_neighbour(1, :))
            if (gl_Gran_neighbour(1, i) < 1 .or. gl_Gran_neighbour(2, i) < 1) CYCLE
            j1 = gl_Gran_neighbour(1, i)
            if( ANY(gl_Cell_neighbour(:, j1) == gl_Gran_neighbour(2, i)) == .FALSE.) then
                pause "ERROR 1775"
            end if
        end do
        
        if (par_developer_info) print *, "Geometry_check: Proverka 4"
        ! Проверяем сумму нормалей на площадь грани в каждой ячейке
        N1 = size(gl_all_Cell(1, :))
        do i = 1, N1
            c = 0.0
            do j = 1, 6
                j1 = gl_Cell_gran(j, i)
                a1 = 1
                if (j1 == 0) CYCLE
                if (gl_Gran_neighbour(2, j1) == i) a1 = -1
                c = c + gl_Gran_square(j1) * gl_Gran_normal(:, j1) * a1
                if( (norm2(gl_Gran_normal(:, j1)) - 1.0) > 0.0001 ) print*, "ERROR 5678ijhgftydwed3"
            end do
            !print*, "________________"
            !print*, gl_Cell_type(i), gl_Cell_number(1, i), gl_Cell_number(2, i), gl_Cell_number(3, i), norm2(c), gl_Cell_Volume(i)
            !print*, gl_x(gl_all_Cell(1, i)), gl_y(gl_all_Cell(1, i)), gl_z(gl_all_Cell(1, i))
            !print*, gl_x(gl_all_Cell(2, i)), gl_y(gl_all_Cell(2, i)), gl_z(gl_all_Cell(2, i))
            !print*, gl_x(gl_all_Cell(3, i)), gl_y(gl_all_Cell(3, i)), gl_z(gl_all_Cell(3, i))
            !print*, gl_x(gl_all_Cell(4, i)), gl_y(gl_all_Cell(4, i)), gl_z(gl_all_Cell(4, i))
            !print*, gl_x(gl_all_Cell(5, i)), gl_y(gl_all_Cell(5, i)), gl_z(gl_all_Cell(5, i))
            !print*, gl_x(gl_all_Cell(6, i)), gl_y(gl_all_Cell(6, i)), gl_z(gl_all_Cell(6, i))
            !print*, gl_x(gl_all_Cell(7, i)), gl_y(gl_all_Cell(7, i)), gl_z(gl_all_Cell(7, i))
            !print*, gl_x(gl_all_Cell(8, i)), gl_y(gl_all_Cell(8, i)), gl_z(gl_all_Cell(8, i))
            !print*, "________________"
            !PAUSE
            if (norm2(c) > 1E-10) print*, "Error 2966765gvbfv ", gl_Cell_type(i), gl_Cell_number(1, i), gl_Cell_number(2, i), gl_Cell_number(3, i), norm2(c)
        end do
        
        
        

        if (par_developer_info) print *, "Geometry_check: OK"

        ! Какие-то временные проверки можно сюда вставлять
        !print *, "B111:"
        !a1 = gl_Cell_B(1,1,1)
        !do a2 = 1, 8
        !	print*, gl_all_Cell(a2, a1)
        !end do
        
        ! Проверяем параметры среды
        N1 = size(gl_all_Cell(1, :))
        do i = 1, N1
            if (gl_Cell_par(1, i) <= 0.0) gl_Cell_par(1, i) = 0.000001
            if (gl_Cell_par_MF(1, 1, i) <= 0.0) then
                !print*, gl_Cell_par_MF(:, 1, i)
                !print*, "Center = ", gl_Cell_center(:, i)
                !PAUSE "ERROR ro sfergh453"
                gl_Cell_par_MF(1, 1, i) = 0.000001
            end if
            if (gl_Cell_par_MF(1, 2, i) <= 0.0) gl_Cell_par_MF(1, 2, i) = 0.000001
            if (gl_Cell_par_MF(1, 3, i) <= 0.0) gl_Cell_par_MF(1, 3, i) = 0.000001
            if (gl_Cell_par_MF(1, 4, i) <= 0.0) gl_Cell_par_MF(1, 4, i) = 0.000001
            
            if (gl_Cell_par(5, i) <= 0.0) gl_Cell_par(5, i) = 0.000001
            if (gl_Cell_par_MF(5, 1, i) <= 0.0) gl_Cell_par_MF(5, 1, i) = 0.000001
            if (gl_Cell_par_MF(5, 2, i) <= 0.0) gl_Cell_par_MF(5, 2, i) = 0.000001
            if (gl_Cell_par_MF(5, 3, i) <= 0.0) gl_Cell_par_MF(5, 3, i) = 0.000001
            if (gl_Cell_par_MF(5, 4, i) <= 0.0) gl_Cell_par_MF(5, 4, i) = 0.000001
        end do
    end subroutine Geometry_check

    ! ****************************************************************************************************************************************************
    !! Основная программа

	subroutine Download_setka(num)
        ! Правильная загрузка сетки из файла (при этом все нужные функции вызываются)
        use STORAGE
        use GEO_PARAM
        integer(4), intent(in) :: num
        
        print*, " Download_setka 1"
        call Read_setka_bin(num)            ! Считываем сетку с файла 
        print*, " Download_setka 2"
        call Find_Surface()                ! Ищем поверхности, которые будем выделять (вручную)
        print*, " Download_setka 3"
        call calc_all_Gran()               ! Программа расчёта объёмов ячеек, площадей и нормалей граней (обязательна здесь)
        print*, " Download_setka 4"
        call Find_inner()                  ! Находит ячейки внутри небольшой сферы, в которых счёт будет происходить отдельно (обязательно после 
                                        ! предыдущей функции)
        print*, " Download_setka 5"
        call Geometry_check()              ! Проверка геометрии сетки, чтобы не было ошибок в построении
        print*, " Download_setka 6"
        call Initial_conditions()  ! Задаём граничные условия. Нужно проверить с каким chi задаётся (если проводится перенормировка на каждом шаге)
        print*, " Download_setka 7"
        call Find_TVD_sosed()
        print*, " Download_setka 8"
	
	end subroutine Download_setka
	
	subroutine Get_MK_to_MHD(koeff)
        ! Получим результаты работы Монте-Карло для использования их в МГД расчётах
        use Interpolate2
        LOGICAL, intent(in), OPTIONAL :: koeff
        integer(4) :: i, num, s1, j
        real(8) :: x, y, z, dd
        real(8) :: PAR(9)     ! Выходные параметры
        real(8) :: PAR_MOMENT(par_n_moment, par_n_sort)
        real(8) :: PAR_k(5)
        real(8) :: MAS_PUI(2)
        LOGICAL :: koeff_local

        koeff_local = .True.
        if(present(koeff)) then
            koeff_local = koeff
        end if
        
        num = 3
        gl_Cell_par_MK(1:5, :, :) = 0.0
        gl_Cell_par_MK(6:10, :, :) = 1.0
        
        print*, "subroutine Get_MK_to_MHD()"
        
        !call Int2_Read_bin(2)  ! Загрузка файла интерполяции

        
        dd = 1.0
        
        do i = 1, size(gl_Cell_center(1, :))
            x = gl_Cell_center(1, i)
            y = gl_Cell_center(2, i)
            z = gl_Cell_center(3, i)
            call Int2_Get_par_fast(x, y, z, num, PAR, PAR_MOMENT, PAR_k, MAS_PUI = MAS_PUI)
            !if(par_n_sort /= 4) STOP "ERROR 7890okjhyuio98765rtyuikgyui"
            gl_Cell_par_MK(1:5, :, i) = PAR_MOMENT(1:5, :) * dd  ! Если сортов - 4
            gl_Cell_par_pui(:, i) = MAS_PUI

            if(koeff_local) gl_Cell_par_MK(6:10, 1, i) = PAR_k(:) * dd  !TODO Нужно ли интерполировать коэффициенты? Или их лучше оставить на сетке?
            

            if(num < 1) then
                !print*, "mini-problem with interpolation", x, y, z
                num = 3
                call Int2_Get_tetraedron_inner(x, y, z, num)
                if(num < 1) STOP "ERROR  89uyfvbnm[;.xsw4567u"
                s1 = int2_all_tetraendron_point(1, num)
                gl_Cell_par_MK(1:5, 1:4, i) = int2_Moment(1:5, :, s1) * dd
                if(koeff_local) gl_Cell_par_MK(6:10, 1, i) = int2_Moment_k(:, s1) * dd
            end if
            
            
            
            do j = 1, par_n_sort
                if(gl_Cell_par_MK(1, j, i) < 0.0000001) then
                    gl_Cell_par_MK(1, j, i) = 0.0000001
                end if
            end do
            
        end do
        
        print*, "end subroutine Get_MK_to_MHD()"

        return 
        
        call Int2_Read_bin(3)  ! Загрузка файла интерполяции
        print*, "_____ menyaem fayl ____"
        
        dd = 1.0 - dd
        
        do i = 1, size(gl_Cell_center(1, :))
            x = gl_Cell_center(1, i)
            y = gl_Cell_center(2, i)
            z = gl_Cell_center(3, i)
            call Int2_Get_par_fast(x, y, z, num, PAR, PAR_MOMENT, PAR_k)
            !if(par_n_sort /= 4) STOP "ERROR 7890okjhyuio98765rtyuikgyui"
            gl_Cell_par_MK(1:5, :, i) = gl_Cell_par_MK(1:5, :, i) + PAR_MOMENT(1:5, :) * dd  ! Если сортов - 4
            gl_Cell_par_MK(6:10, 1, i) = gl_Cell_par_MK(6:10, 1, i) + PAR_k(:) * dd
            
            if(num < 1) then
                !print*, "mini-problem with interpolation", x, y, z
                num = 3
                call Int2_Get_tetraedron_inner(x, y, z, num)
                if(num < 1) STOP "ERROR  89uyfvbnm[;.xsw4567u"
                s1 = int2_all_tetraendron_point(1, num)
                gl_Cell_par_MK(1:5, :, i) = gl_Cell_par_MK(1:5, :, i) + int2_Moment(1:5, :, s1) * dd
                gl_Cell_par_MK(6:10, 1, i) = gl_Cell_par_MK(6:10, 1, i) + int2_Moment_k(:, s1) * dd
            end if
            
            
            
            do j = 1, 4
                if(gl_Cell_par_MK(1, j, i) < 0.000001) then
                    gl_Cell_par_MK(1, j, i) = 0.000001
                end if
            end do
            
        end do
        
        call Int2_Dell_interpolate()
	end subroutine Get_MK_to_MHD
    
	subroutine PRINT_ALL()

        print*, "= PRINT_ALL()"
        ! Печатает текущую сетку (разные файлы)
        print*, "C1"
        call Print_par_2D()
        print*, "C1"
        call Print_surface_2D()
        print*, "C2"
        call Print_Setka_2D()
        print*, "C1"
        call Print_Setka_y_2D()
        print*, "C1"
        call Print_TS_3D()
        print*, "C3"
        call Print_Contact_3D()
        call Print_BS_3D()
        print*, "C1"
        call Print_surface_y_2D()
        print*, "C1"
        call Print_par_y_2D()
        print*, "C4"
        call Print_par_1D()
        print*, "C1"
        
        ! Body of PRINT_ALL
	
	end subroutine PRINT_ALL
    

    program MIK
        use STORAGE
        use GEO_PARAM
        use Interpolate
        use Interpolate2
        use Surface_setting
        !@cuf use MY_CUDA
        use Monte_Karlo
        use MY_CUDA_smooth
        implicit none

        integer(4) :: i, NGRAN, j, nn, name_do, name_posle
        integer :: istat, s1,s2,s3, n1,n2, n3, k
        real(8) :: local1, F(9), b, xx, yy, zz
        integer :: ierrSync, ierrAsync
        
        print*, "START PROGRAM"
        
        
        
        name_do = 186
        name_posle = 999

        
        !call Test_koordinate(1.0_8, 0.0_8, 0.0_8)
        !call Test_raspadnik()
        
        !call Smooth_kvadr3(0.0_8, 0.0_8, 0.0_8,  3.0_8, 3.0_8, 0.0_8,  4.0_8, 4.5_8, 0.0_8,   5.0_8, 3.0_8, 0.0_8, 7.0_8, 0.0_8, 0.0_8,   xx, yy, zz)
        !call Smooth_kvadr4(0.0_8, 0.0_8, 0.0_8,  3.0_8, 2.0_8, 0.0_8,  4.0_8, 0.3_8, 0.0_8,  &
        !	5.0_8, 0.5_8, 0.0_8,   7.0_8, 0.0_8, 0.0_8,   8.0_8, 2.0_8, 0.0_8,   9.0_8, 0.0_8, 0.0_8,   xx, yy, zz)
        !print*, xx, yy, zz
        call Play_iter_algoritm()
        
        
        STOP
        
        ! НИЖЕ СТАРЫЕ ЭЛЕМЕНТЫ УПРАВЛЕНИЯ (потом можно будет удалить)
        
        !name_do = 224
        !name_posle = 225
        
        !@cuf call CUDA_info()
        !call EXIT()

        ! Процесс построения сетки (не менять, все шаги необходимы для корректной работы)
        !call Set_STORAGE()                 ! Выделяем память под все массимы рограммы
        !call Build_Mesh_start()            ! Запускаем начальное построение сетки (все ячейки связываются, но поверхности не выделены)
        
        call Read_setka_bin(name_do)            ! Либо считываем сетку с файла (при этом всё равно вызывается предыдущие функции под капотом)
        !
        !
        call Find_Surface()                ! Ищем поверхности, которые будем выделять (вручную)
        call calc_all_Gran()               ! Программа расчёта объёмов ячеек, площадей и нормалей граней (обязательна здесь)
        call Find_inner()                  ! Находит ячейки внутри небольшой сферы, в которых счёт будет происходить отдельно (обязательно после 
                                        ! предыдущей функции)
        
        
        !call Print_All_Points()
        !call Print_A_Ray(1,3)
        !call Print_Point_Plane()
        call Geometry_check()              ! Проверка геометрии сетки, чтобы не было ошибок в построении
        !call Print_Cell(1, 2, 1, "C")
        !call Print_all_Cell()
        
        
        ! call Print_par_2D()
        ! call Print_surface_2D()
        ! call Print_Setka_2D()
        ! call Print_Setka_y_2D()
        ! call Print_TS_3D()
        ! call Print_Contact_3D()
        ! call Print_surface_y_2D()
        
        
        !call Print_cell_and_neighbour(1,2,1)
        !call Print_gran()
        !call Print_all_surface("C")
        !print *, "Start_GD_2"
        !call Start_GD_2(10000)

        !call Start_GD_2(100000)
        !call Initial_conditions()

        !print*, "Size inner gran and cell = " , size(gl_all_Gran_inner), size(gl_all_Cell_inner)
        
        call Initial_conditions()  ! Задаём граничные условия. Нужно проверить с каким chi задаётся (если проводится перенормировка на каждом шаге)
        call Find_TVD_sosed()

        ! Перенормировка сортов для более быстрого счёта. Перенормируем первый и второй сорта. Делаем это ВЕЗДЕ
        !gl_Cell_par_MF(2:4, 1, :) = gl_Cell_par_MF(2:4, 1, :) / (par_chi/par_chi_real)
        !gl_Cell_par_MF(1, 1, :) =  gl_Cell_par_MF(1, 1, :) * (par_chi/par_chi_real)**2

        !gl_Cell_par_MF(2:4, 2, :) = gl_Cell_par_MF(2:4, 2, :) / (3.0)
        !gl_Cell_par_MF(1, 2, :) =  gl_Cell_par_MF(1, 2, :) * (3.0)**2

        !call Print_par_2D()
        
        !call Set_Interpolate_main()
        
        
        !call Read_interpolate_bin(163)
        !
        !call Re_interpolate()
        !
        !call Dell_Interpolate()
        
        !par_kk12 = 1.0_8
        !par_kk31 = 1.3_8
        !par_R_LEFT = -460.0
        
        !print*, par_n_IA, par_n_IB, par_R_inner
        
        !par_n_IA = 16
        !par_n_IB = 18
        !par_R_inner = 9.0_8
        
        !par_kk13 = 0.3
        
        !call Start_MGD_move()
        !pause
        
        ! call Set_Interpolate_main()
        ! call Surface_setup()
        !   i = 1
        ! call Get_Cell_Interpolate( -226.83946549657_8, 97.1182680029345_8, -3.5212386335546_8, i) 
        ! call Dell_Interpolate()
        
        
        !call Find_tetraedr_Interpolate(0.9_8, 0.04_8, 0.001_8, i)
        
        !PAUSE
        
        !call CUDA_START_MGD_move()
        
        !call Set_Interpolate_main()       ! Проверим интерполяцию
        
        !call Set_Interpolate_main()
        !call Streem_line(50.0_8, 0.001_8, 0.001_8, 1)
        !call Streem_line(50.0_8, 20.001_8, 20.001_8, 2)
        
        !call CUDA_START_GD_3()
        
        !i = 211242
        !call Print_Cell(gl_Cell_number(1, i), gl_Cell_number(2, i), gl_Cell_number(3, i), gl_Cell_type(i))
        
        
        !pause
        
        
        
        ! do i = 1, 1
        !    if(mod(i, 1000) == 0) print*, i
        !    call Start_GD_3(1)
        !    call Start_GD_3_inner(1)
        !    if(mod(i, 10000) == 0) then
            !   call Save_setka_bin(14)
        !       call Print_surface_2D()
        !       call Print_Setka_2D()
        !       call Print_par_2D()
        !    end if
        ! end do
        
        !   call Find_tetraedr_Interpolate(0.9_8, 0.04_8, 0.001_8, 4)
        
        !call Interpolate_point(4.0_8, 2.0_8, 0.003_8, F, istat, nn)
        
        
        !call Set_Interpolate_main()
        !call Save_interpolate_bin(163)
        

        ! Перенормировка сортов обратно
        !gl_Cell_par_MF(2:4, 1, :) = gl_Cell_par_MF(2:4, 1, :) * (par_chi/par_chi_real)
        !gl_Cell_par_MF(1, 1, :) =  gl_Cell_par_MF(1, 1, :) / (par_chi/par_chi_real)**2

        !gl_Cell_par_MF(2:4, 2, :) = gl_Cell_par_MF(2:4, 2, :) * (3.0)
        !gl_Cell_par_MF(1, 2, :) =  gl_Cell_par_MF(1, 2, :) / (3.0)**2

        !call Re_interpolate()
        
        !call Dell_Interpolate()
        
        
        ! ИНТЕРПОЛЯЦИОННЫЙ БЛОК _______________________________________________________________________________________
        !pause
        
        
        call Int2_Set_Interpolate()
        call Int2_Initial()
        call Int2_Set_interpol_matrix()
        call Int2_Print_setka_2()
        call Int2_Print_sosed()
        
        !n1 = size(int2_Cell_B(:, 1, 1))
        !n2 = size(int2_Cell_B(1, :, 1))
        !n3 = size(int2_Cell_B(1, 1, :))
        !
        !do k = 1, n3
        !	do j = 1, n2
        !		do i = 1, n1
        !			if(int2_Cell_B(i, j, k) == 167674) print*, "167674", i, j, k
        !		end do
        !	end do
        !end do
        !
        !print*, n1, n2, n3
        !
        !call Int2_Print_Cell(167674)
        
        
        !print*, "DDD", int2_all_tetraendron(:, 1195701)
        !do i = 1, 4
        !print*, int2_gran_sosed(int2_all_tetraendron(i, 1195701))
        !end do
        
        !pause
        
        !i = 3
        !call Int2_Get_tetraedron(10.0_8, 10.0_8, 10.0_8, i)
        !print*, "PAR"
        !call Int2_Get_par(10.0_8, 10.0_8, 10.0_8, i, F)
        
        !print*, "M_K"
        !call M_K_start()
        !call Int2_Print_tetraedron(109)
        call Int_2_Print_par_2D(0.0_8, 0.0_8, 1.0_8, -0.000001_8, 1)
        call Int_2_Print_par_2D(0.0_8, 1.0_8, 0.0_8, -0.000001_8, 2)
        
        !pause
        !call Set_Interpolate_main()
        !call Save_interpolate_bin(186)
        
        !call Surf_Save_bin(186)
        !call Surf_Read_setka_bin(186)
        !
        !do i = 1, 100
        !	if(mod(i, 20) == 0) print*, i
        !	call Surf_Set_surf(20.0_8)
        !end do
        !
        !do i = 1, 25
        !	if(mod(i, 20) == 0) print*, i
        !	call Surf_Set_surf(1.0_8)
        !end do
        !
        !do i = 1, 5
        !	if(mod(i, 20) == 0) print*, i
        !	call Surf_Set_surf(0.2_8)
        !end do
        !
        !do i = 1, 10
        !	if(mod(i, 20) == 0) print*, i
        !	call Surf_Set_surf(0.03_8)
        !end do
        !
        !call calc_all_Gran()               ! Программа расчёта объёмов ячеек, площадей и нормалей граней (обязательна здесь)
        !call Read_interpolate_bin(186)
        !print*, "Read interpolate"
        !pause
        !call Re_interpolate()
        !print*, "RE interpolate"
        !pause
    
        ! call Print_surface_2D()
        ! call Print_Setka_2D()
        ! call Print_Setka_y_2D()
        ! Конец интерполяционного блока _____________________________________________________________________________
        
        !call Print_all_surface("C")
        !call Print_all_surface("B")
        !call Print_all_surface("T")
        
        !    call Print_par_2D()
        ! 	call Print_par_y_2D()
        ! 	call Print_surface_y_2D()
        !    call Save_setka_bin(name_posle)
        !    ! Variables
        !    call Print_Contact_3D()
        ! 	call Print_TS_3D()
        ! 	call Print_Setka_3D_part()
        ! 	call Save_param(name_posle)
        ! 	pause
	
	end program MIK
	

	subroutine Play_iter_algoritm()
	    !! Основная программа запуска задачи
        !* Есть различные сценарий запуска, в зависимости от того, что именно считается (МГД, Монте-Карло и т.д.)
	    use STORAGE
        use GEO_PARAM
	    use Interpolate2
	    use Surface_setting
	    use Monte_Karlo
        use PUI
	    integer(4) :: step, name, name2, i, j, k, ij1, ij2, ij3, ii
		real(8) :: PAR(9)     ! Выходные параметры
	    integer(4) :: num  ! Тетраэдр, в котором предположительно находится точка (num по умолчанию должен быть равен 3)
	    real(8) :: PAR_MOMENT(18, par_n_sort), uz, u, cp, aa

        aa = 10E30

        print*, "test = ", aa/(10E28)
		
        
		name = 597 !548 544    536 !? 534 Номер файла основной сетки   533 - до PUI
        ! 551 - до изменения chi
        ! 574 до изменения сечения на стебегенса
        !? 535 - до того, как поменять определение давления в PUI

        ! 334! 307  ! С 237 надо перестроить сетку ! Имя основной сетки  начало с 224
		! 132 до экспериментов  134  138
        ! 152 (до того, как поменять сгущение после HP)		
		! 509 до изменения разрешения на Димино		
		! 519 до включения газовой динамики на контакте
		! 249 до фотоионизации
        ! 258 с гелием (только ввёл) до того, как поменять схему	
		! 259 вторая и третья область HLLC
		! 269 до того, как убрал осреднение скоростей в распаднике
		! 277 закончил мгд (278 аналогичное)
		! 284 до изменения граничных условий в источнике
		! 288 до того, как увеличили число точек в нуле
		! 298 до того, как изменили схему на гелиопаузе
		! 304 до изменения знака поля внутри
		! 315 перед тем, как перестроить сетку
		name2 = 34 !?19   9 8 Номер интерполяционного файла сетки с источниками    8 - до PUI
        ! 26 до изменения числа атомов водорода
		!name3 = 237  ! Имя сетки интерполяции для М-К
		step = 3  !? 3 Номер алгоритма

		!PAR_MOMENT = 0.0
		!call Int2_Read_bin(name2)
		
		!u = 0.5
		!cp = 11.0
		!uz = MK_Velosity_1(u, cp)
		!print*, "MK = ", u/cp, MK_int_1(u, cp), uz * MK_sigma(uz)
		!pause
		
		!Делаем файл интерполяции для ИГОРЯ
		!call Download_setka(name)  ! Загрузка основной сетки (со всеми нужными функциями)
		
		!do i = 1, size(gl_Cell_par2(1, :))
		!	if( ieee_is_nan(gl_Cell_par2(1, i)) .or. ieee_is_normal(gl_Cell_par2(1, i)) == .False. ) then
		!		print*, "No = ", gl_Cell_par2(1, i)
		!		print*, gl_Cell_center(:, i)
		!	end if
		!	
		!end do
		
		! He-лий
		!do i = 1, size(gl_Cell_par(1, :))
		!	gl_Cell_par(1, i) = gl_Cell_par(1, i) - gl_Cell_par2(1, i)
		!	
		!	if(gl_zone_Cell(i) == 1 .or. gl_zone_Cell(i) == 2) then
		!		gl_Cell_par(5, i) = gl_Cell_par(5, i) / (2.0 * gl_Cell_par(1, i) + 1.5 * gl_Cell_par2(1, i)) * 2.0 * gl_Cell_par(1, i)
		!	else
		!		gl_Cell_par(5, i) = gl_Cell_par(5, i) / (2.0 * gl_Cell_par(1, i) + gl_Cell_par2(1, i)) * 2.0 * gl_Cell_par(1, i)
		!	end if
		!end do
		
		
		!print*, "Interpol"
		!call Int2_Set_Interpolate()      ! Выделение памяти под	сетку интерполяции
		!call Int2_Initial()			     ! Создание сетки интерполяции
		!call Int2_Set_interpol_matrix()	 ! Заполнение интерполяционной матрицы в каждом тетраэдре с помощью Lapack
		!call Int2_Save_interpol_for_all_MHD(259)
		
		! 243 конец первой итерации
		
		
	
	    ! Описание алгоритма
        ! 1) Считаем МГД на подробной сетке + создаём файл интерполяции №1 + создаём файл поверхностей
        ! 2) Строим грубую сетку, двигаем поверхности разрывов на основе подробной сетки
		!	 интерполируем значения в центрах ячеек + создаём файл интерполяции №2 для этой сетки
		! + Считаем Монте-Карло на грубой сетке (файл интерполяции №2) + обновляем файл интерполяции №2
	    !	 далее переходим к шагу 1) и используем файл интерполяции №2 в МГД 
		
		print*, "********************************************************************************************"
        print*, "Vypolnyaetsya shag nomer ", step
		print*, "********************************************************************************************"
		
		if(step == -2) then  ! Меняем структуру сетки
			! СОЗДАЁМ СЕТКУ
			! задаём параметры мини-сетки
			par_l_phi = 60
			par_m_A = 30! 30      ! Количество лучей A в плоскости
            par_m_BC = 13! 18      ! Количество лучей B/C в плоскости
            par_m_O = 15! 17      ! Количество лучей O в плоскости
            par_m_K = 8! 7      ! Количество лучей K в плоскости
            ! Количество точек по лучам A
            par_n_TS =  33! 57! 52! 26                    ! Количество точек до TS (TS включается)
            par_n_HP =  63! 87! 82! 40                 ! Количество точек HP (HP включается)  всё от 0 считается
            par_n_BS =  103! 89! 127! 152! 60! 5                 ! Количество точек BS (BS включается)
            par_n_END = 121! 98! 142! 167! 72! 6                ! Количество точек до конца сетки (конец включается)
            par_n_IA =  20! 34! 12                   ! Количество точек, которые входят во внутреннюю область
	        par_n_IB =  22! 36! 14                   ! Количество точек, которые входят во внутреннюю область (с зазором)
			par_kk1 = 2.0_8
			par_kk14 = 1.0_8 !0.5_8
			par_kk12 = 1.0_8!1.3_8
			!par_kk13 = 1.7_8
			par_kk2 = 1.83_8
			par_al1 = 1.0_8
            par_R_END = 550.0
			
			call Set_STORAGE()                 ! Выделяем память под все массивы программы
            call Build_Mesh_start()            ! Запускаем начальное построение сетки (все ячейки связываются, но поверхности не выделены)
    
            call Find_Surface()                ! Ищем поверхности, которые будем выделять (вручную)
            call calc_all_Gran()               ! Программа расчёта объёмов ячеек, площадей и нормалей граней (обязательна здесь)
            call Find_inner()                  ! Находит ячейки внутри небольшой сферы, в которых счёт будет происходить отдельно (обязательно после 
                                               ! предыдущей функции)
            call Geometry_check()              ! Проверка геометрии сетки, чтобы не было ошибок в построении
			
			print*, "CENTR = ", norm2(gl_Cell_center(:, gl_Cell_A(1, 1, 1)))
			print*, "CENTR 2 = ", norm2(gl_Cell_center(:, gl_Cell_A(2, 1, 1)))
			
			! Считываем файл поверхностей разрыва и двигаем сетку
			
			call Surf_Read_setka_bin(name)
			
	        print*, "Nachinaem dvizhenie setki"
	        do i = 1, 270
	        	call Surf_Set_surf(1.0_8)
	        end do
	        
	        do i = 1, 55
	        	call Surf_Set_surf(1.0_8)
	        end do
	        
	        do i = 1, 50
	        	call Surf_Set_surf(0.2_8)
	        end do
	        
	        do i = 1, 55
	        	call Surf_Set_surf(0.03_8)
			end do
			
			do i = 1, 65
	        	call Surf_Set_surf(0.003_8)
			end do
			
			print*, "CENTR = ", norm2(gl_Cell_center(:, gl_Cell_A(1, 1, 1)))
			print*, "CENTR 2 = ", norm2(gl_Cell_center(:, gl_Cell_A(2, 1, 1)))
			
            call Find_Surface()                ! Ищем поверхности, которые будем выделять (вручную)
            call calc_all_Gran()               ! Программа расчёта объёмов ячеек, площадей и нормалей граней (обязательна здесь)
            call Find_inner()                  ! Находит ячейки внутри небольшой сферы, в которых счёт будет происходить отдельно (обязательно после 
                                               ! предыдущей функции)
            call Geometry_check()              ! Проверка геометрии сетки, чтобы не было ошибок в построении
			
			print*, "Dvizhenie setki zaversheno" 
			! Считываем файл интерполяции и интерполируем переменные на мини-сетку
			print*, "A"
			call Int2_Read_bin(name)
			print*, "B"
			call Int2_Re_interpol()
			print*, "C"
			call Int2_Dell_interpolate()
			print*, "D"
			
			call PRINT_ALL()
			print*, "Save"
			call Save_param(name + 1)
			call Surf_Save_bin(name + 1)   ! Сохранение поверхностей разрыва
			call Save_setka_bin(name + 1)
			print*, "Save 2"
			
		else if(step == -1) then  ! Перестройка сетки (когда поменяли структуру)
			call Download_setka(name)  ! Загрузка основной сетки
			call Surf_Read_setka_bin(name)
			
			par_triple_point_2 = 7.0 * par_pi_8/40.0
			
	        print*, "Nachinaem dvizhenie setki"
	        do i = 1, 120
	        	call Surf_Set_surf(20.0_8)
	        end do
	        
	        do i = 1, 35
	        	call Surf_Set_surf(1.0_8)
	        end do
	        
	        do i = 1, 10
	        	call Surf_Set_surf(0.2_8)
	        end do
	        
	        do i = 1, 15
	        	call Surf_Set_surf(0.03_8)
			end do
			
			do i = 1, 20
	        	call Surf_Set_surf(0.003_8)
			end do
			
            call Find_Surface()                ! Ищем поверхности, которые будем выделять (вручную)
            call calc_all_Gran()               ! Программа расчёта объёмов ячеек, площадей и нормалей граней (обязательна здесь)
            call Find_inner()                  ! Находит ячейки внутри небольшой сферы, в которых счёт будет происходить отдельно (обязательно после 
                                               ! предыдущей функции)
            call Geometry_check()              ! Проверка геометрии сетки, чтобы не было ошибок в построении

			call Int2_Read_bin(name)  ! Загрузка файла интерполяции
			call Int2_Re_interpol()
			call Int2_Dell_interpolate()
			
			
			! Печатаем сетку (для просмотра)
			call PRINT_ALL()
			call Save_setka_bin(name + 1)
			
			
			print*, "Dvizhenie setki zaversheno" 
			
			
		else if(step == 0) then  ! Интерполяция
			print*, "Download"
			call Download_setka(name)  ! Загрузка основной сетки
			! СОХРАНЕНИЕ
			print*, "Save"
			call Surf_Save_bin(name)   ! Сохранение поверхностей разрыва
			
			call Int2_Set_Interpolate()      ! Выделение памяти под	сетку интерполяции
	        call Int2_Initial()			     ! Создание сетки интерполяции
	        call Int2_Set_interpol_matrix()	 ! Заполнение интерполяционной матрицы в каждом тетраэдре с помощью Lapack
			call Int2_Save_bin(name)			 ! Сохранение полной сетки интерполяции
		else if(step == 1) then ! CUDA MHD ----------------------------------------------------------------------------
			! ЗАГРУЗКА СЕТКИ
			print*, "-----   Download"
			call Download_setka(name)  ! Загрузка основной сетки (со всеми нужными функциями)

			call Initial_conditions()  ! Задаём граничные условия на внутренней сфере

			print*, "-----   Download int"
			call Int2_Read_bin(name2)  ! Загрузка файла интерполяции
			!call Int2_Save_interpol_for_all_MK(name2)
			!call Int2_Dell_interpolate()

            if(par_PUI) then
                call PUI_Set()                !! PUI
                call PUI_f_Set()              !! PUI
                call PUI_f_Set2()             !! PUI
                call PUI_Read_bin(name2)
                call PUI_Read_f_bin(name2)
                call PUI_F_integr_Set()       !! PUI
                call PUI_Read_for_MK_bin(name2)   !! PUI
            end if
			
			call Get_MK_to_MHD() ! Заполняем центры ячеек параметрами водорода и коэффициентами интерполяции
			
			!call Int_2_Print_par_2D_set()
			call Int_2_Print_par_1D()
			
			!do ii = 1, 20
			!    do i = par_n_HP - 6, par_n_HP + 6
			!	    do j = 1, 20
			!		    do k = 1, 40
			!	            ij = gl_Cell_A(i, j, k)
			!	            ij2 = gl_Cell_A(i-1, j, k)
			!	            ij3 = gl_Cell_A(i+1, j, k)
			!			    gl_Cell_par(:, ij) = (gl_Cell_par(:, ij2) + gl_Cell_par(:, ij3))/2.0
			!		    end do
			!	    end do
			!    end do
			!end do
			
			
			
			! Перенормируем параметры плазмы в гелиосфере
			do i = 1, size(gl_Cell_par(1, :))
				if(gl_zone_Cell(i) <= 2) then
					gl_Cell_par(2:4, i) = gl_Cell_par(2:4, i) / (par_chi/par_chi_real)
				    gl_Cell_par(1, i) = gl_Cell_par(1, i) * (par_chi/par_chi_real)**2
				    gl_Cell_par2(1, i) = gl_Cell_par2(1, i) * (par_chi/par_chi_real)**2
				    gl_Cell_par(9, i) = gl_Cell_par(9, i) * (par_chi/par_chi_real)**2
				end if
			end do
			
			call PRINT_ALL()
			!return
			!go to 101
			
			! call Print_tok_layer()
			! return
			
		    !@cuf call CUDA_info()
			print*, "------------------------------ START ----------------------------"
			!@cuf call CUDA_START_MGD_move_MK() ! РАСЧЁТЫ
			!!call Start_MGD_move_MK()
			
			! Перенормируем параметры плазмы обратно
			do i = 1, size(gl_Cell_par(1, :))
				if(gl_zone_Cell(i) <= 2) then
					gl_Cell_par(2:4, i) = gl_Cell_par(2:4, i) * (par_chi/par_chi_real)
				    gl_Cell_par(1, i) = gl_Cell_par(1, i) / (par_chi/par_chi_real)**2
				    gl_Cell_par2(1, i) = gl_Cell_par2(1, i) / (par_chi/par_chi_real)**2
				    gl_Cell_par(9, i) = gl_Cell_par(9, i) / (par_chi/par_chi_real)**2
				end if
			end do
			
			
			! Печатаем сетку (для просмотра)
			call PRINT_ALL()
			! СОХРАНЕНИЕ
			print*, "Save"
			call Save_param(name + 1) ! + 1
			call Surf_Save_bin(name + 1)   ! Сохранение поверхностей разрыва  + 1
			call Save_setka_bin(name + 1) !  + 1
			print*, "Save 2"
			!
			call Int2_Dell_interpolate()
			print*, "Save 3"
			call Int2_Set_Interpolate()      ! Выделение памяти под	сетку интерполяции
			print*, "Save 4"
	        call Int2_Initial()			     ! Создание сетки интерполяции
			print*, "Save 5"
			call Int2_Set_interpol_matrix()	 ! Заполнение интерполяционной матрицы в каждом тетраэдре с помощью Lapack
			print*, "Save 6"
			!call Int2_Save_bin(name + 1)			 ! Сохранение полной сетки интерполяции  + 1
			print*, "Save 7"
			!! Сохранение сетки для общего пользования
		    !! call Int2_Save_interpol_for_all_MHD(name) !  + 1
			!
		
		else if(step == 2) then  !  МОНТЕ-КАРЛО ----------------------------------------------------------------------------------------
			! РАБОТА С МОНТЕ-КАРЛО
            ! Создаём все необходимые файлы из файла основной сетки
			call Download_setka(name)  ! Загрузка основной сетки (со всеми нужными функциями)
            print*, "** A **"
			call Surf_Save_bin(name)   ! Сохранение поверхностей разрыва
			call Int2_Set_Interpolate()      ! Выделение памяти под	сетку интерполяции
	        call Int2_Initial()			     ! Создание сетки интерполяции
			call Int2_Set_interpol_matrix()	 ! Заполнение интерполяционной матрицы в каждом тетраэдре с помощью Lapack
			call Int2_Save_bin(name)			 ! Сохранение полной сетки интерполяции
            print*, "** B **"
			call Int2_Dell_interpolate()
			call Dell_STORAGE()
			
            print*, "** Create Mesh **"
			! СОЗДАЁМ СЕТКУ
			! задаём параметры мини-сетки
			par_l_phi = 40
			par_m_A = 20! 30      ! Количество лучей A в плоскости
            par_m_BC = 10! 18      ! Количество лучей B/C в плоскости
            par_m_O = 10! 17      ! Количество лучей O в плоскости
            par_m_K = 8! 7      ! Количество лучей K в плоскости
            ! Количество точек по лучам A
            par_n_TS =  27! 26                    ! Количество точек до TS (TS включается)
            par_n_HP =  37! 40                 ! Количество точек HP (HP включается)  всё от 0 считается
            par_n_BS =  57! 60! 5                 ! Количество точек BS (BS включается)
            par_n_END = 65! 72! 6                ! Количество точек до конца сетки (конец включается)
            par_n_IA =  12! 12                   ! Количество точек, которые входят во внутреннюю область
	        par_n_IB =  14! 14                   ! Количество точек, которые входят во внутреннюю область (с зазором)
			
            print*, "** C **"
			call Set_STORAGE()                 ! Выделяем память под все массивы программы
            print*, "** D **"
            call Build_Mesh_start()            ! Запускаем начальное построение сетки (все ячейки связываются, но поверхности не выделены)
            print*, "** F **"
            call Find_Surface()                ! Ищем поверхности, которые будем выделять (вручную)
            call calc_all_Gran()               ! Программа расчёта объёмов ячеек, площадей и нормалей граней (обязательна здесь)
            call Find_inner()                  ! Находит ячейки внутри небольшой сферы, в которых счёт будет происходить отдельно (обязательно после 
                                               ! предыдущей функции)
            call Geometry_check()              ! Проверка геометрии сетки, чтобы не было ошибок в построении
			print*, "** G **"
			! Считываем файл поверхностей разрыва и двигаем сетку
			
			call Surf_Read_setka_bin(name)
			
	        print*, "Nachinaem dvizhenie setki"
	        do i = 1, 120
	        	call Surf_Set_surf(20.0_8)
	        end do
	        
	        do i = 1, 55
	        	call Surf_Set_surf(1.0_8)
	        end do
	        
	        do i = 1, 35
	        	call Surf_Set_surf(0.2_8)
	        end do
	        
	        do i = 1, 35
	        	call Surf_Set_surf(0.03_8)
			end do
			
			do i = 1, 35
	        	call Surf_Set_surf(0.003_8)
			end do

            do i = 1, 10
	        	call Surf_Set_surf(0.0001_8)
			end do
			
            call Find_Surface()                ! Ищем поверхности, которые будем выделять (вручную)
            call calc_all_Gran()               ! Программа расчёта объёмов ячеек, площадей и нормалей граней (обязательна здесь)
            call Find_inner()                  ! Находит ячейки внутри небольшой сферы, в которых счёт будет происходить отдельно (обязательно после 
                                               ! предыдущей функции)
            call Geometry_check()              ! Проверка геометрии сетки, чтобы не было ошибок в построении
			
			print*, "Dvizhenie setki zaversheno" 
			! Считываем файл интерполяции и интерполируем переменные на мини-сетку
			call Int2_Read_bin(name)
			call Int2_Re_interpol()
			call Int2_Dell_interpolate()
			par_n_moment = 19
			print*, "End Interpolatiya" 
			! Печатаем сетку (для просмотра)
			call PRINT_ALL()
			
			! Делаем файл интерполяции № 2 из мини-сетки
			call Int2_Set_Interpolate()      ! Выделение памяти под	сетку интерполяции
		    call Int2_Initial()			     ! Создание сетки интерполяции
		    call Int2_Set_interpol_matrix()	 ! Заполнение интерполяционной матрицы в каждом тетраэдре с помощью Lapack
			
			call Int2_Print_setka_2()
		    call Int2_Print_sosed()
			
			
			call Dell_STORAGE()  ! Удаляем основную сетку (т.к. считается только монте-карло - для экономии памяти)
			call Int_2_Print_par_2D(0.0_8, 0.0_8, 1.0_8, -0.000001_8, 1)
			call Int2_Print_center()
			call Int2_Print_point_plane()
			call Int2_Print_setka_2()
            
            ! Считываем нормальную сетку

            call Download_setka(name)  ! Загрузка основной сетки (со всеми нужными функциями)

            !! PUI 
            if(par_PUI) then
                call PUI_Set()                !! PUI
                call PUI_f_Set()              !! PUI
                call PUI_f_Set2()             !! PUI
                call PUI_Read_bin(name2)
                call PUI_Read_f_bin(name2)
                call PUI_Culc_h0()            !! PUI
                call PUI_F_integr_Set()       !! PUI
                call PUI_Read_for_MK_bin(name2)   !! PUI
                call Print_par_1D_PUI()
            end if
			
			! СЧИТАЕМ Монте-Карло на мини-сетке
			print*, "START MK"
			!!call Helium_off()
			call M_K_start()
            ! call PUI_proverka(20.0_8, 0.0_8, 0.0_8)

			!call M_K_sum()
			call Int2_culc_k()
			!!call Helium_on()
			!call Int_2_Print_par_2D(0.0_8, 0.0_8, 1.0_8, -0.000001_8, 1)
			!call Int_2_Print_par_2D(0.0_8, 1.0_8, 0.0_8, -0.000001_8, 2)
			call Int_2_Print_par_1D()


			
			! Сохраняем интерполяционный файл - мини - сетки
			call Int2_Save_bin(name2 + 1)			 ! Сохранение полной сетки интерполяции
			call Int2_Save_interpol_for_all_MK(name2 + 1)
			call Int_2_Print_par_2D_set()	
            ! call PUI_Save_bin(name2 + 1)
            if(par_PUI) then
                call PUI_calc_Sm()
                call PUI_Save_bin(name2 + 1)
            end if
			
            call PRINT_ALL()
            
		else if(step == 3) then  !----------------------------------------------------------------------------------------
			call Download_setka(name)  ! Загрузка основной сетки (со всеми нужными функциями)
            call Surf_Save_bin(name)   ! Сохранение поверхностей разрыва
            call Surf_Read_setka_bin(name)

			call Int2_Read_bin(name2)
			call PUI_Set()
			call PUI_f_Set()
			call PUI_f_Set2()
			call PUI_Read_bin(name2)

            call Int_2_Print_par_1D()
            
            ! open(11, file = "pui_num_tetr.bin", FORM = 'BINARY')
            ! write(11) size(pui_num_tetr)
            ! write(11) pui_num_tetr
            ! close(11)

            !return
			! call PUI_Read_bin(9)

			! call PUI_calc_Sm()
            ! call PUI_Save_bin(name2 + 1)

			call Culc_f_pui()
            call Cut_f_pui()
            call PUI_Save_f_bin(name2)
            ! call PUI_Read_f_bin(name2 + 1)

            call PUI_Culc_h0()
            call PUI_F_integr_Set()

            call PUI_F_integr_Culc()
            call PUI_n_T_culc()
            call PUI_Save_for_MK_bin(name2)
            ! call PUI_Read_for_MK_bin(name2 + 1)
            
            
			call PUI_print(1, 13.0_8, 0.00001_8, 0.00001_8)
			call PUI_print(2, 20.0_8, 0.00001_8, 0.00001_8)
			call PUI_print(3, 17.0_8, 10.00001_8, 0.00001_8)
			call PUI_print(4, 5.0_8, 0.00001_8, 0.00001_8)
			call PUI_print(5, 1.0_8, 0.00001_8, 0.00001_8)
			call PUI_print(6, 10.0_8, 0.00001_8, 0.00001_8)
			call PUI_print(7, 0.0_8, 10.00001_8, 0.00001_8)
			call PUI_print(8, -17.0_8, 0.00001_8, 0.00001_8)
			call PUI_print(9, -50.0_8, 0.00001_8, 0.00001_8)
			call PUI_print(10, -40.0_8, 30.00001_8, 0.00001_8)
			call PUI_print(11, 24.0_8, 0.00001_8, 0.00001_8)
			call PUI_print(12, 24.5_8, 0.00001_8, 0.00001_8)
            
            ! call Int_2_Print_par_1D()
            call Print_par_1D_PUI()

			!call PUI_print(2, 20.0_8, 0.00001_8, 0.00001_8)
			!call PUI_print(3, 17.0_8, 10.00001_8, 0.00001_8)
        else if(step == 4) then
            !? Посчитали МГД, теперь ходим посчитать PUI, используя новые поля плазмы и старые S+, S-
			! СОЗДАЁМ СЕТКУ
			! задаём параметры мини-сетки
			par_l_phi = 40
			par_m_A = 20! 30      ! Количество лучей A в плоскости
            par_m_BC = 10! 18      ! Количество лучей B/C в плоскости
            par_m_O = 10! 17      ! Количество лучей O в плоскости
            par_m_K = 8! 7      ! Количество лучей K в плоскости
            ! Количество точек по лучам A
            par_n_TS =  27! 26                    ! Количество точек до TS (TS включается)
            par_n_HP =  37! 40                 ! Количество точек HP (HP включается)  всё от 0 считается
            par_n_BS =  57! 60! 5                 ! Количество точек BS (BS включается)
            par_n_END = 65! 72! 6                ! Количество точек до конца сетки (конец включается)
            par_n_IA =  12! 12                   ! Количество точек, которые входят во внутреннюю область
	        par_n_IB =  14! 14                   ! Количество точек, которые входят во внутреннюю область (с зазором)
			
			call Set_STORAGE()                 ! Выделяем память под все массивы программы
            call Build_Mesh_start()            ! Запускаем начальное построение сетки (все ячейки связываются, но поверхности не выделены)
    
            call Find_Surface()                ! Ищем поверхности, которые будем выделять (вручную)
            call calc_all_Gran()               ! Программа расчёта объёмов ячеек, площадей и нормалей граней (обязательна здесь)
            call Find_inner()                  ! Находит ячейки внутри небольшой сферы, в которых счёт будет происходить отдельно (обязательно после 
                                               ! предыдущей функции)
            call Geometry_check()              ! Проверка геометрии сетки, чтобы не было ошибок в построении
			
			! Считываем файл поверхностей разрыва и двигаем сетку
			
			call Surf_Read_setka_bin(name)
			
	        print*, "Nachinaem dvizhenie setki"
	        do i = 1, 120
	        	call Surf_Set_surf(20.0_8)
	        end do
	        
	        do i = 1, 55
	        	call Surf_Set_surf(1.0_8)
	        end do
	        
	        do i = 1, 35
	        	call Surf_Set_surf(0.2_8)
	        end do
	        
	        do i = 1, 35
	        	call Surf_Set_surf(0.03_8)
			end do
			
			do i = 1, 35
	        	call Surf_Set_surf(0.003_8)
			end do
			
            call Find_Surface()                ! Ищем поверхности, которые будем выделять (вручную)
            call calc_all_Gran()               ! Программа расчёта объёмов ячеек, площадей и нормалей граней (обязательна здесь)
            call Find_inner()                  ! Находит ячейки внутри небольшой сферы, в которых счёт будет происходить отдельно (обязательно после 
                                               ! предыдущей функции)
            call Geometry_check()              ! Проверка геометрии сетки, чтобы не было ошибок в построении
			
			print*, "Dvizhenie setki zaversheno" 
			! Считываем файл интерполяции и интерполируем переменные на мини-сетку
			call Int2_Read_bin(name2)
			call Int2_Re_interpol()
			call Int2_Dell_interpolate()
			par_n_moment = 19
			print*, "End Interpolatiya" 
			! Печатаем сетку (для просмотра)
			call PRINT_ALL()
			
			! Делаем файл интерполяции № 2 из мини-сетки
			call Int2_Set_Interpolate()      ! Выделение памяти под	сетку интерполяции
		    call Int2_Initial()			     ! Создание сетки интерполяции
		    call Int2_Set_interpol_matrix()	 ! Заполнение интерполяционной матрицы в каждом тетраэдре с помощью Lapack
			
			call Int2_Print_setka_2()
		    call Int2_Print_sosed()
			
			
			call Dell_STORAGE()  ! Удаляем основную сетку (т.к. считается только монте-карло - для экономии памяти)
            
            ! Считываем нормальную сетку

            call Download_setka(name)  ! Загрузка основной сетки (со всеми нужными функциями)

            !! PUI 
            call PUI_Set()                !! PUI
            call PUI_f_Set()              !! PUI
            call PUI_f_Set2()             !! PUI
            call PUI_Read_bin(name2)

            call Culc_f_pui()
            call Cut_f_pui()
            call PUI_Save_f_bin(name2)

            call PUI_Culc_h0()
            call PUI_F_integr_Set()

            call PUI_F_integr_Culc()
            call PUI_n_T_culc()
            call PUI_Save_for_MK_bin(name2)
            
            
			! call PUI_print(40, 8.0_8, 0.00001_8, 0.00001_8)
			! call PUI_print(50, 10.0_8, 0.00001_8, 0.00001_8)
			! call PUI_print(60, 12.0_8, 0.00001_8, 0.00001_8)
			! call PUI_print(70, 14.0_8, 0.00001_8, 0.00001_8)
            ! call PUI_print(80, 16.0_8, 0.00001_8, 0.00001_8)
            ! call PUI_print(90, 18.0_8, 0.00001_8, 0.00001_8)
            ! call PUI_print(100, 20.0_8, 0.00001_8, 0.00001_8)
            ! call PUI_print(110, 22.0_8, 0.00001_8, 0.00001_8)
            ! call PUI_print(1, 14.2_8, 0.00001_8, 0.00001_8)
            ! call PUI_print(2, 14.5_8, 0.00001_8, 0.00001_8)
            ! call PUI_print(3, 15.1_8, 0.00001_8, 0.00001_8)
            ! call PUI_print(4, 15.2_8, 0.00001_8, 0.00001_8)
            ! call PUI_print(5, 15.3_8, 0.00001_8, 0.00001_8)
            ! call PUI_print(6, 15.5_8, 0.00001_8, 0.00001_8)
            
            call Print_par_1D_PUI()
        else if(step == 5) then
            ! Считаем дивергенцию скорости для игоря
            call Download_setka(name)  ! Загрузка основной сетки (со всеми нужными функциями)
            call Cucl_div_V()
            !print*, gl_Cell_par_div
            !call Int2_Read_bin(name2)  ! Загрузка файла интерполяции
			call Int2_Set_Interpolate()      ! Выделение памяти под	сетку интерполяции
	        call Int2_Initial()			     ! Создание сетки интерполяции
			call Int2_Set_interpol_matrix()	 ! Заполнение интерполяционной матрицы в каждом тетраэдре с помощью Lapack
            call Int_2_Print_par_1D()
            !print*, int2_Cell_par_div

            call Int2_Save_interpol_for_all_MHD(name + 1)
            call Send_request()

        else if(step == 6) then
            ! Забыл добавить в источники ещё один сорт водорода
            call Download_setka(name)  ! Загрузка основной сетки (со всеми нужными функциями)
            print*, "A------------------"
            call Int2_Read_bin(name2)
            print*, "B------------------"
            call PUI_Set()
            print*, "C------------------"
			call PUI_f_Set()
            print*, "D------------------"
			call PUI_f_Set2()
            print*, "E------------------"
            call PUI_Read_bin(name2)
            print*, "F------------------"
            ! call Helium_off()
            print*, "F1------------------"
			call Int2_culc_k(.False.)
            print*, "F2------------------"
			! call Helium_on()
            print*, "G------------------"
            call Int2_Save_bin(name2)
            call Int_2_Print_par_1D()
        end if !----------------------------------------------------------------------------------------
		
        101     continue
		print*, "********************************************************************************************"
        print*, "Programma uspeshno zavershena!"
		print*, "********************************************************************************************"
	
	end subroutine Play_iter_algoritm
	
	
	! Сохраненные сетки
	! 145 - до введения ТВД, но вариант ещё не установлен в хвосте
	! 148 - до введения неизотропного с.в.
	! 163 - до переделки сетки (чтобы вернуться придётся возвращать старое движение точек)
	! 177 - до ручной поправки
	! 184 - до включения null_bn везде на контакте
	! 186 - до изменения движения контакта в хвосте (до 15 )   ! Считаем последней рабочей версией
	
	! 200 - Новая сетка   203 до актвного движение контакта

    !! Виды комментариев
    !? Привет
    !TODO Привет
    !* param Привет
    !/// Привет
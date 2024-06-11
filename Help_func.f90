

    subroutine Save_setka_bin(num)  ! Сохранение сетки в бинарном файле
    ! Variables
    use STORAGE
    use GEO_PARAM
    implicit none
    integer, intent(in) :: num
    character(len=3) :: name
	integer :: i
    
    write(unit=name,fmt='(i3.3)') num
    
    open(1, file = "save_" // name // ".bin", FORM = 'BINARY')
    
    write(1)  par_l_phi, par_m_A, par_m_BC, par_m_O, par_m_K, par_triple_point, par_n_TS, par_n_HP, &
        par_n_BS, par_n_END
    
    
    write(1) size(gl_x(:))
    write(1) gl_x, gl_y, gl_z
    
    ! Печатаем параметры сетки
    
    write(1) size(gl_Cell_par(1, :))     ! 9 штук!  ro, u, v, w, p, bx, by, bz
    write(1) gl_Cell_par(:, :)
    
    
    ! Для  будующих дополнений
    write(1) 1
    write(1) gl_Cell_par_MF
    
    write(1) 1
	write(1) par_n_IA, par_n_IB, par_R_inner
	
    write(1) 1
	write(1) par_R_END, par_R_LEFT, par_kk1, par_kk12, par_kk2, par_kk3, par_kk31
	
    write(1) 1
	write(1) par_al1
	
    write(1) 1
	write(1) par_kk13
	
    write(1) 1
    write(1) par_triple_point_2
	
	
    write(1) 1
	do i = 1, size(gl_Cell_par2(1, :))
		    write(1) gl_Cell_par2(1, i)
	end do
	
    write(1) 1
    write(1) par_kk14
	
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    
    
    close(1)
    end subroutine Save_setka_bin
    
    
    subroutine Read_setka_bin(num)  ! Сохранение сетки в бинарном файле
    ! Variables
    use STORAGE
    use GEO_PARAM
    implicit none
    integer, intent(in) :: num
    character(len=3) :: name
    integer :: n, i
    logical :: exists
    
    write(unit=name,fmt='(i3.3)') num
    
    inquire(file="save_" // name // ".bin", exist=exists)
    
    if (exists == .False.) then
		pause "net faila!!!"
        STOP "net faila!!!"
    end if
    
    
    open(1, file = "save_" // name // ".bin", FORM = 'BINARY', ACTION = "READ")
    
    read(1)  par_l_phi, par_m_A, par_m_BC, par_m_O, par_m_K, par_triple_point, par_n_TS, par_n_HP, &
        par_n_BS, par_n_END
    
    call Set_STORAGE()
    call Build_Mesh_start()
    
    read(1) n                            ! Число узлов
    if (size(gl_x) /= n) print*, "ERROR 75 tfevh 765678 23412313"
    read(1) gl_x, gl_y, gl_z
    
    ! Печатаем параметры сетки
    
    read(1) n
    read(1) gl_Cell_par(:, :)
    
    
    ! Для  будующих дополнений
    read(1) n
    if (n == 1) then
        read(1) gl_Cell_par_MF
	end if
    
    read(1) n
	if (n == 1) then
        read(1) par_n_IA, par_n_IB, par_R_inner
	end if
	
	
    read(1) n
	if (n == 1) then
        read(1) par_R_END, par_R_LEFT, par_kk1, par_kk12, par_kk2, par_kk3, par_kk31
	end if
	
	
    read(1) n
	if (n == 1) then
        read(1) par_al1
	end if
	
	
    read(1) n
	if (n == 1) then
        read(1) par_kk13
	end if
	
	
    read(1) n
	if (n == 1) then
		read(1) par_triple_point_2
	else
		par_triple_point_2 = par_triple_point
	end if
	
    read(1) n
	if (n == 1) then
		do i = 1, size(gl_Cell_par2(1, :))
		    read(1) gl_Cell_par2(1, i)
		end do
	end if
	
	
    read(1) n
	if (n == 1) then
		read(1) par_kk14
	end if
	
	
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    
    
    close(1)
	end subroutine Read_setka_bin
	
	subroutine Save_param(num)
	! Сохраняет параметры программы (перечисленные в namelist в модуле GEO_PARAM) 
	USE GEO_PARAM
	implicit none
    integer, intent(in) :: num
    character(len=5) :: name
    
    write(unit=name,fmt='(i5.5)') num
    
    open(1, file = "Setka_param_" // name // ".txt")
	
	WRITE (1, SETKA_PARAM) 
	
	close(1)
	
	end subroutine Save_param

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
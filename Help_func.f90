

    subroutine Save_setka_bin(num)  ! Сохранение сетки в бинарном файле
    ! Variables
    use STORAGE
    use GEO_PARAM
    implicit none
    integer, intent(in) :: num
    character(len=3) :: name
	integer :: i
    
    write(unit=name,fmt='(i3.3)') num
    
    open(1, file = par_NAME // "save_" // name // ".bin", FORM = 'BINARY')
    
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
    
    inquire(file= par_NAME//"save_" // name // ".bin", exist=exists)
    
    if (exists == .False.) then
		pause "net faila!!!"
        STOP "net faila!!!"
    end if
    
    
    open(1, file = par_NAME//"save_" // name // ".bin", FORM = 'BINARY', ACTION = "READ")
    
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

    if(ALLOCATED(gl_Cell_par_pui) == .False.) then
        allocate( gl_Cell_par_pui(2, size(gl_Cell_A(:,:,:)) + size(gl_Cell_B(:,:,:)) + size(gl_Cell_C(:,:,:)) ) )
        gl_Cell_par_pui = 0.0
    end if
    
    close(1)
	end subroutine Read_setka_bin
	
	subroutine Save_param(num)
	! Сохраняет параметры программы (перечисленные в namelist в модуле GEO_PARAM) 
	USE GEO_PARAM
	implicit none
    integer, intent(in) :: num
    character(len=5) :: name
    
    write(unit=name,fmt='(i5.5)') num
    
    open(1, file = par_NAME // "Setka_param_" // name // ".txt")
	
	WRITE (1, SETKA_PARAM) 
	
	close(1)
	
	end subroutine Save_param



	
module Surface_setting  
	! Модуль для сохранения положения поверхностей на одной сетке
    ! для последующей переинтерполяции их на другую
	! Модуль работает и завершён
	USE GEO_PARAM
	USE STORAGE
	USE My_func
	
	integer(4) :: Surf_l_phi   ! Количество вращений плоскости 
	integer(4) :: Surf_m_A      ! Количество лучей A в плоскости
    integer(4) :: Surf_m_BC      ! Количество лучей B/C в плоскости
    integer(4) :: Surf_m_O      ! Количество лучей O в плоскости
    integer(4) :: Surf_m_K      ! Количество лучей K в плоскости
    real(8) :: Surf_triple_point
	real(8) :: Surf_al1
	
	
	real(8), allocatable :: Surf_phi(:)   ! Набо углов плоскостей вращения
	
	real(8), allocatable :: Surf_TS_the(:, :)   ! (:, :)   (j, k)  -  угол в этой плоскости
	real(8), allocatable :: Surf_TS_r(:, :)     ! (:, :)   (j, k)  -  Расстояние до центра
	
	real(8), allocatable :: Surf_BS_the(:, :)   ! (:, :)   (j, k)  -  угол в этой плоскости
	real(8), allocatable :: Surf_BS_r(:, :)     ! (:, :)   (j, k)  -  Расстояние до центра
	
	! Для A - Лучей
	real(8), allocatable :: Surf_HP_the(:, :)   ! (:, :)   (j, k)  -  угол в этой плоскости
	real(8), allocatable :: Surf_HP_r_A(:, :)     ! (:, :)   (j, k)  -  Расстояние до центра
	
	! Для остальных лучей
	real(8), allocatable :: Surf_HP_x(:, :)   ! (:, :)   (j, k)  -  х координата в этой плоскости
	real(8), allocatable :: Surf_HP_r(:, :)     ! (:, :)   (j, k)  -  Расстояние до оси симметрии
	
	contains
	
	subroutine Surf_Save_bin(num)
	! Сохраняет положения поверхностей основной сетки (которая находится в модуле Storage)
    use STORAGE
    use GEO_PARAM
    implicit none
    integer, intent(in) :: num
    character(len=5) :: name
	integer :: i, j, k, N1, N2, N3, yzel
	real(8):: kord(3)
    
    write(unit=name, fmt='(i5.5)') num
    
    open(1, file = "Surf_Save" // name // ".bin", FORM = 'BINARY')
	
    
    write(1) par_l_phi, par_m_A, par_m_BC, par_m_O, par_m_K, par_triple_point, par_al1
	
	do k = 1, par_l_phi
		
		yzel = gl_RAY_A(5, 3, k)
		kord = (/gl_x(yzel), gl_y(yzel), gl_z(yzel)/)
			
		write(1) polar_angle(kord(2), kord(3))
		
		! A - лучи
		do j = 1, par_m_A
			
			! TS
			yzel = gl_RAY_A(par_n_TS, j, k)
			kord = (/gl_x(yzel), gl_y(yzel), gl_z(yzel)/)
			
			write(1) polar_angle(kord(1), sqrt(kord(2)**2 + kord(3)**2)), norm2(kord)
			
			
			! Гелиопауза
			yzel = gl_RAY_A(par_n_HP, j, k)
			kord = (/gl_x(yzel), gl_y(yzel), gl_z(yzel)/)
			write(1) norm2(kord)
			
			! BS
			yzel = gl_RAY_A(par_n_BS, j, k)
			kord = (/gl_x(yzel), gl_y(yzel), gl_z(yzel)/)
			write(1) norm2(kord)
			
		end do
		
		! B - лучи
		do j = 1, par_m_BC
			
			! the = (j - 1) * par_pi_8/2.0/(par_m_A - 1)
			
			! TS
			yzel = gl_RAY_B(par_n_TS, j, k)
			kord = (/gl_x(yzel), gl_y(yzel), gl_z(yzel)/)
			write(1) polar_angle(kord(1), sqrt(kord(2)**2 + kord(3)**2)),  norm2(kord)
			
			
			! Гелиопауза
			yzel = gl_RAY_B(par_n_HP, j, k)
			kord = (/ gl_x(yzel), gl_y(yzel), gl_z(yzel)/)
			write(1) kord(1), sqrt(kord(2)**2 + kord(3)**2)
			
		end do
		
		
		! O - лучи
		do j = 1, par_m_O
			
			! Гелиопауза
			yzel = gl_RAY_O(1, j, k)
			kord = (/ gl_x(yzel), gl_y(yzel), gl_z(yzel)/)
			write(1) kord(1), sqrt(kord(2)**2 + kord(3)**2)
			
		end do
		
		
		! K - лучи
		do j = par_m_K, 1, -1
			! TS
			yzel = gl_RAY_K(par_n_TS, j, k)
			kord = (/gl_x(yzel), gl_y(yzel), gl_z(yzel)/)
			write(1) polar_angle(kord(1), sqrt(kord(2)**2 + kord(3)**2)), norm2(kord)
		end do
		
		
	end do
	
	
	close(1)
	
	
	end subroutine Surf_Save_bin
	
	
	subroutine Surf_Read_setka_bin(num)
	! Variables
    use STORAGE
    use GEO_PARAM
    implicit none
    integer, intent(in) :: num
    character(len=5) :: name
    integer :: n, i, j, k, N1, N2, N3
    logical :: exists
    
    write(unit=name,fmt='(i5.5)') num
    
    inquire(file="Surf_Save" // name // ".bin", exist=exists)
    
    if (exists == .False.) then
		pause "net faila!!!  8idjgeye0-0987fhjdkeye "
        STOP "net faila!!!"
    end if
    
    
    
    
    open(1, file = "Surf_Save" // name // ".bin", FORM = 'BINARY', ACTION = "READ")
    
    read(1)  Surf_l_phi, Surf_m_A, Surf_m_BC, Surf_m_O, Surf_m_K, Surf_triple_point, Surf_al1
    
	allocate( Surf_phi(Surf_l_phi) )
	
	allocate( Surf_TS_the(Surf_m_A + Surf_m_BC + Surf_m_K, Surf_l_phi) )
	allocate( Surf_TS_r(Surf_m_A + Surf_m_BC + Surf_m_K, Surf_l_phi) )
	
	allocate( Surf_BS_the(Surf_m_A, Surf_l_phi) )
	allocate( Surf_BS_r(Surf_m_A, Surf_l_phi) )
	
	allocate( Surf_HP_the(Surf_m_A, Surf_l_phi) )
	allocate( Surf_HP_r_A(Surf_m_A, Surf_l_phi) )
	
	allocate( Surf_HP_x(Surf_m_BC + Surf_m_O, Surf_l_phi) )
	allocate( Surf_HP_r(Surf_m_BC + Surf_m_O, Surf_l_phi) )
	
	do k = 1, Surf_l_phi
		read(1) Surf_phi(k)
		
		! A
		do j = 1, Surf_m_A
			read(1) Surf_TS_the(j, k), Surf_TS_r(j, k)
			Surf_HP_the(j, k) = Surf_TS_the(j, k)
			Surf_BS_the(j, k) = Surf_TS_the(j, k)
			read(1) Surf_HP_r_A(j, k)
			read(1) Surf_BS_r(j, k)
		end do
		
		! B
		do j = 1, Surf_m_BC
			read(1) Surf_TS_the(j + Surf_m_A, k), Surf_TS_r(j + Surf_m_A, k)
			read(1) Surf_HP_x(j, k), Surf_HP_r(j, k)
		end do
		
		! O
		do j = 1, Surf_m_O
			read(1) Surf_HP_x(j + Surf_m_BC, k), Surf_HP_r(j + Surf_m_BC, k)
		end do
		
		! K - лучи
		do j = 1, Surf_m_K
			read(1) Surf_TS_the(j + Surf_m_A + Surf_m_BC, k), Surf_TS_r(j + Surf_m_A + Surf_m_BC, k)
		end do
		
	end do
	
	
	close(1)
	
	end subroutine Surf_Read_setka_bin
	
	
	real(8) pure function Surf_Get_TS(x, y, z)
	! Variables
    implicit none
    real(8), intent(in) :: x, y, z
    real(8) :: the, phi, t1, t2, a, b, c, d
	integer :: k1, k, j, j1
	
	 phi = polar_angle(y, z)
	 the = polar_angle(x, sqrt(y**2 + z**2) )
	 
	 k1 = Surf_l_phi
	 do k = 1, Surf_l_phi
		 if(Surf_phi(k) > phi) then
			 k1 = k
			 EXIT
		 end if
	 end do
	 
	 t1 = (phi - Surf_phi(k1 - 1))/( Surf_phi(k1) - Surf_phi(k1 - 1) )
	 
	 j1 = size(Surf_TS_the(:, 1)) 
	 do j = 1, size(Surf_TS_the(:, 1)) 
		 if(Surf_TS_the(j, k1) >= the) then
			 j1 = j
			 EXIT
		 end if
	 end do
	 
	 if (j1 == 1) then
		 j1 = 2
		 t2 = 0.0
	 else
		t2 = ( the - Surf_TS_the(j1 - 1, k1) )/( Surf_TS_the(j1, k1) - Surf_TS_the(j1 - 1, k1) )
	 end if
	 
	 a = Surf_TS_r(j1 - 1, k1 - 1)
	 b = Surf_TS_r(j1, k1 - 1)
	 c = Surf_TS_r(j1 - 1, k1)
	 d = Surf_TS_r(j1, k1)
	
	 
	 Surf_Get_TS = ((1.0 - t1) * a + t1 * c) * (1.0 - t2) + ((1.0 - t1) * b + t1 * d) * t2
	 
	end function Surf_Get_TS
	
	
	real(8) pure function Surf_Get_HP(x, y, z)
	! Variables
    implicit none
    real(8), intent(in) :: x, y, z
    real(8) :: the, phi, t1, t2, a, b, c, d
	integer :: k1, k, j, j1
	
	phi = polar_angle(y, z)
	
	k1 = Surf_l_phi
	 do k = 1, Surf_l_phi
		 if(Surf_phi(k) > phi) then
			 k1 = k
			 EXIT
		 end if
	 end do
	 
	 t1 = (phi - Surf_phi(k1 - 1))/( Surf_phi(k1) - Surf_phi(k1 - 1) )
	
	if (x >= 0) then
	 
	 the = polar_angle(x, sqrt(y**2 + z**2) )
	 
	 j1 = size(Surf_HP_the(:, 1)) 
	 do j = 1, size(Surf_HP_the(:, 1)) 
		 if(Surf_HP_the(j, k1) > the) then
			 j1 = j
			 EXIT
		 end if
	 end do
	 
	 t2 = ( the - Surf_HP_the(j1 - 1, k1) )/( Surf_HP_the(j1, k1) - Surf_HP_the(j1 - 1, k1) )
	 
	 a = Surf_HP_r_A(j1 - 1, k1 - 1)
	 b = Surf_HP_r_A(j1, k1 - 1)
	 c = Surf_HP_r_A(j1 - 1, k1)
	 d = Surf_HP_r_A(j1, k1)
	
	 
	 Surf_Get_HP = ((1.0 - t1) * a + t1 * c) * (1.0 - t2) + ((1.0 - t1) * b + t1 * d) * t2
	else
		
	 j1 = size(Surf_HP_x(:, k1)) 
	 do j = 1, size(Surf_HP_x(:, k1)) 
		 if(Surf_HP_x(j, k1) < x) then
			 j1 = j
			 EXIT
		 end if
	 end do
	 
	 if(j1 == 1) then
		 t2 = (x)/( Surf_HP_x(j1, k1))
		 a = Surf_HP_r_A(size(Surf_HP_the(:, 1)) , k1 - 1)
		 b = Surf_HP_r(j1, k1 - 1)
		 c = Surf_HP_r_A(size(Surf_HP_the(:, 1)) , k1)
		 d = Surf_HP_r(j1, k1)
		 Surf_Get_HP = ((1.0 - t1) * a + t1 * c) * (1.0 - t2) + ((1.0 - t1) * b + t1 * d) * t2
	 else
		 
		t2 = ( x - Surf_HP_x(j1 - 1, k1) )/( Surf_HP_x(j1, k1) - Surf_HP_x(j1 - 1, k1) )
		 a = Surf_HP_r(j1 - 1, k1 - 1)
		 b = Surf_HP_r(j1, k1 - 1)
		 c = Surf_HP_r(j1 - 1, k1)
		 d = Surf_HP_r(j1, k1)
		 Surf_Get_HP = ((1.0 - t1) * a + t1 * c) * (1.0 - t2) + ((1.0 - t1) * b + t1 * d) * t2
	 
	 end if
		
	end if
	
	 
	end function Surf_Get_HP
	
	real(8) pure function Surf_Get_BS(x, y, z)
	! Variables
    implicit none
    real(8), intent(in) :: x, y, z
    real(8) :: the, phi, t1, t2, a, b, c, d
	integer :: k1, k, j, j1
	
	 phi = polar_angle(y, z)
	 the = polar_angle(x, sqrt(y**2 + z**2) )
	 
	 k1 = Surf_l_phi
	 do k = 1, Surf_l_phi
		 if(Surf_phi(k) > phi) then
			 k1 = k
			 EXIT
		 end if
	 end do
	 
	 t1 = (phi - Surf_phi(k1 - 1))/( Surf_phi(k1) - Surf_phi(k1 - 1) )
	 
	 j1 = size(Surf_BS_the(:, 1)) 
	 do j = 1, size(Surf_BS_the(:, 1)) 
		 if(Surf_BS_the(j, k1) > the) then
			 j1 = j
			 EXIT
		 end if
	 end do
	 
	 t2 = ( the - Surf_BS_the(j1 - 1, k1) )/( Surf_BS_the(j1, k1) - Surf_BS_the(j1 - 1, k1) )
	 
	 a = Surf_BS_r(j1 - 1, k1 - 1)
	 b = Surf_BS_r(j1, k1 - 1)
	 c = Surf_BS_r(j1 - 1, k1)
	 d = Surf_BS_r(j1, k1)
	 
	 Surf_Get_BS = ((1.0 - t1) * a + t1 * c) * (1.0 - t2) + ((1.0 - t1) * b + t1 * d) * t2
	 
	end function Surf_Get_BS
	
	subroutine Surf_Set_surf(dr)
	! Передвигаем поверхности в соответствии с запомненными положениями
	use STORAGE
    use GEO_PARAM
	use My_func
	! Variables
	real(8), intent(in) :: dr  ! Шаг движения
	real(8) :: R_TS, proect, vel(3), R_HP, R_BS, KORD(3), dist, ddt, ER(3)
    integer :: yzel, N1, N2, N3, i, j, k, yzel2
    real(8) :: the, phi, r, x, y, z, rr, xx, x2, y2, z2, rrr, r1, r2, r3, r4, rd, kk13
	
	! Body of Set_surf
	
	N3 = size(gl_RAY_A(1, 1, :))
    N2 = size(gl_RAY_A(1, :, 1))
    N1 = size(gl_RAY_A(:, 1, 1))
	
	! Цикл движения точек на лучах А  ************************************************************
    do k = 1, N3
        do j = 1, N2
            
            if (k /= 1 .and. j == 1) then
                    CYCLE
            end if
            
            ! Вычисляем координаты текущего луча в пространстве
            the = (j - 1) * par_pi_8/2.0/(N2 - 1)
            phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
            
            ! TS
            yzel = gl_RAY_A(par_n_TS, j, k)
            KORD(1) = gl_x(yzel); KORD(2) = gl_y(yzel); KORD(3) = gl_z(yzel) 
			R_TS = norm2(KORD)  ! Новое расстояние до TS
			
			if (R_TS > Surf_Get_TS(KORD(1), KORD(2), KORD(3))) then
				R_TS = R_TS - dr
			else
				R_TS = R_TS + dr
			end if
			
			
			
            ! HP
            yzel = gl_RAY_A(par_n_HP, j, k)
			KORD(1) = gl_x(yzel); KORD(2) = gl_y(yzel); KORD(3) = gl_z(yzel) 
            R_HP = norm2(KORD)  ! Новое расстояние до HP
			
			if (R_HP > Surf_Get_HP(KORD(1), KORD(2), KORD(3))) then
				R_HP = R_HP - dr
			else
				R_HP = R_HP + dr
			end if
			
			
            ! BS
            yzel = gl_RAY_A(par_n_BS, j, k)
            
			KORD(1) = gl_x(yzel); KORD(2) = gl_y(yzel); KORD(3) = gl_z(yzel) 
            R_BS = norm2(KORD)  ! Новое расстояние до BS
			
			if (R_BS > Surf_Get_BS(KORD(1), KORD(2), KORD(3))) then
				R_BS = R_BS - dr
			else
				R_BS = R_BS + dr
			end if
			
            
            ! Далее обычный цикл нахождения координат точек, такой же, как и при построении сетки
            do i = 1, N1
                
                if (i == 1) then
                    CYCLE
                end if

                yzel = gl_RAY_A(i, j, k)
                ! Вычисляем координаты точки на луче
				
				kk13 = par_kk13 * (par_pi_8/2.0 - the)/(par_pi_8/2.0)  +  1.0 * (the/(par_pi_8/2.0))**2

                ! до TS
				if (i <= par_n_IB) then  ! NEW
						r =  par_R0 + (par_R_inner - par_R0) * (DBLE(i)/(par_n_IB))**par_kk1
				else if (i <= par_n_TS) then  
						r =  par_R_inner + (R_TS - par_R_inner) * (DBLE(i - par_n_IB)/(par_n_TS - par_n_IB))**par_kk12
				!if (i <= par_n_TS) then  ! До расстояния = R_TS
				!    r =  par_R0 + (R_TS - par_R0) * (DBLE(i)/par_n_TS)**par_kk1
				else if (i <= par_n_HP) then  
					r = R_TS + (i - par_n_TS) * (R_HP - R_TS)/(par_n_HP - par_n_TS)
				else if (i == par_n_HP + 1) then 
					rd = R_TS + (par_n_HP - 1 - par_n_TS) * (R_HP - R_TS)/(par_n_HP - par_n_TS)
					r = R_HP + (R_HP - rd)
				else if (i <= par_n_BS) then 
					rd = R_TS + (par_n_HP - 1 - par_n_TS) * (R_HP - R_TS)/(par_n_HP - par_n_TS)
					rd = R_HP + (R_HP - rd)
					r = rd + (R_BS - rd) * sgushenie_1( 1.0_8 * (i - par_n_HP - 1)/(par_n_BS - par_n_HP - 1), kk13 )
					!r = R_HP + (R_BS - R_HP) * angle_cilindr( 1.0_8 * (i - par_n_HP)/(par_n_BS - par_n_HP), kk13 )
					!r = R_HP + (i - par_n_HP) * (R_BS - R_HP)/(par_n_BS - par_n_HP)
				else
					r = R_BS + (par_R_END - R_BS) * (DBLE(i- par_n_BS)/(par_n_END - par_n_BS))**(par_kk2 * (0.55 + 0.45 * cos(the)) )
				end if
				
				if (i == par_n_TS - 1) then
					r = R_TS - 0.5               
				end if

                ! Записываем новые координаты
                gl_x(yzel) = r * cos(the)
                gl_y(yzel) = r * sin(the) * cos(phi)
                gl_z(yzel) = r * sin(the) * sin(phi)
				
            end do
        end do
	end do
    
		
    N3 = size(gl_RAY_B(1, 1, :))
    N2 = size(gl_RAY_B(1, :, 1))
    N1 = size(gl_RAY_B(:, 1, 1))

    ! Цикл движения точек на лучах B  ************************************************************
    do k = 1, N3
        do j = 1, N2
            
            ! Вычисляем координаты текущего луча в пространстве
            the = par_pi_8/2.0 + (j) * par_triple_point/(N2)
			the2 = par_pi_8/2.0 + (j) * par_triple_point_2/(N2)
            phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
            
            ! TS
            yzel = gl_RAY_B(par_n_TS, j, k)
            
			ER(1) = cos(the); ER(2) = sin(the) * cos(phi); ER(3) = sin(the) * sin(phi)
			KORD(1) = gl_x(yzel); KORD(2) = gl_y(yzel); KORD(3) = gl_z(yzel) 
            R_TS = norm2(KORD)  ! Новое расстояние до TS
			
			if (R_TS > Surf_Get_TS(KORD(1), KORD(2), KORD(3))) then
				R_TS = R_TS - dr
			else
				R_TS = R_TS + dr
			end if
            
            ! HP
            yzel = gl_RAY_B(par_n_HP, j, k)
			KORD(1) = gl_x(yzel); KORD(2) = gl_y(yzel); KORD(3) = gl_z(yzel) 
            !R_HP = norm2(KORD)  ! Новое расстояние до HP
			R_HP = norm2(KORD - (R_TS * ER))  ! Новое расстояние до HP от края TS
			
			if ( sqrt(KORD(2)**2 + KORD(3)**2) > Surf_Get_HP(KORD(1), KORD(2), KORD(3))) then
				R_HP = R_HP - dr
			else
				R_HP = R_HP + dr
			end if
			
			!if( k == 1 ) then
			!print*, KORD(1), Surf_Get_HP(KORD(1), KORD(2), KORD(3))
			!pause
			!end if
                
            do i = 1, N1

                if (i == 1) CYCLE
                
                yzel = gl_RAY_B(i, j, k)
                ! Вычисляем координаты точки на луче

                ! до TS
				if (i <= par_n_IB) then  ! NEW
                    r =  par_R0 + (par_R_inner - par_R0) * (DBLE(i)/(par_n_IB))**par_kk1
                else if (i <= par_n_TS) then  ! До расстояния = R_TS
                    r =  par_R_inner + (R_TS - par_R_inner) * (DBLE(i - par_n_IB)/(par_n_TS - par_n_IB))**par_kk12
                !if (i <= par_n_TS) then  ! До расстояния = R_TS
                !    r =  par_R0 + (R_TS - par_R0) * (REAL(i, KIND = 4)/par_n_TS)**par_kk1
				else if (i <= par_n_HP) then  ! До расстояния = par_R_character * 1.3
                    !r = R_TS + (i - par_n_TS) * (R_HP - R_TS) /(par_n_HP - par_n_TS)
					
					r1 = par_R_inner + (R_TS - par_R_inner)
                    r =  (i - par_n_TS) * (R_HP) /(par_n_HP - par_n_TS)
					gl_x(yzel) = r1 * cos(the) + r * cos(the2)
					gl_y(yzel) = r1 * sin(the) * cos(phi) + r * sin(the2) * cos(phi)
					gl_z(yzel) = r1 * sin(the) * sin(phi) + r * sin(the2) * sin(phi)
					CYCLE
				end if
				
				if (i == par_n_TS - 1) then
					r = R_TS - 0.5               
				end if

                ! Записываем новые координаты
                gl_x(yzel) = r * cos(the)
                gl_y(yzel) = r * sin(the) * cos(phi)
                gl_z(yzel) = r * sin(the) * sin(phi)
                
               
            end do
        end do
	end do
    
	
    N3 = size(gl_RAY_C(1, 1, :))
    N2 = size(gl_RAY_C(1, :, 1))
    N1 = size(gl_RAY_C(:, 1, 1))

    ! Цикл движения точек на лучах C  ************************************************************
    do k = 1, N3
        do j = 1, N2
            
            ! Вычисляем координаты текущего луча в пространстве
                phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
                
                ! Вычисляем координаты точки на луче
                x = gl_x(gl_RAY_B(par_n_HP, j, k))
                y = gl_y(gl_RAY_B(par_n_HP, j, k))
                z = gl_z(gl_RAY_B(par_n_HP, j, k))
                rr = (y**2 + z**2)**(0.5)
				
				
				y = gl_y(gl_RAY_B(par_n_HP - 1, j, k))
                z = gl_z(gl_RAY_B(par_n_HP - 1, j, k))
                rd = (y**2 + z**2)**(0.5)
				rr = rr + (rr - rd)
                
                ! BS     Нужно взять положение BS из её положения на крайнем луче A
                yzel = gl_RAY_A(par_n_BS, size(gl_RAY_A(1, :, 1)), k)
				ER(1) = 0.0_8; ER(2) = gl_y(yzel); ER(3) = gl_z(yzel)
                R_BS = norm2(ER)  ! Новое расстояние до BS
				
            
            do i = 1, N1
                
                if(i == 1) CYCLE
                
                yzel = gl_RAY_C(i, j, k)

                if(i == 2) then
					r = rr
				else if (i <= par_n_BS - par_n_HP + 1) then
                    r = rr + (i - 2) * (R_BS - rr)/(par_n_BS - par_n_HP - 1)
                else
                    r = R_BS + (DBLE(i - (par_n_BS - par_n_HP + 1))/(N1 - (par_n_BS - par_n_HP + 1) ))**(0.55 * par_kk2) * (par_R_END - R_BS)
                end if
                
                gl_x(yzel) = x
                gl_y(yzel) = r * cos(phi)
                gl_z(yzel) = r * sin(phi)
                

            end do
        end do
    end do

    N3 = size(gl_RAY_O(1, 1, :))
    N2 = size(gl_RAY_O(1, :, 1))
    N1 = size(gl_RAY_O(:, 1, 1))


    ! Цикл движения точек на лучах O   ************************************************************
    do k = 1, N3
        phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
        
        do j = 1, N2
            
            yzel = gl_RAY_O(1, j, k)
			KORD(1) = 0.0_8; KORD(2) = gl_y(yzel); KORD(3) = gl_z(yzel)
            R_HP = norm2(KORD)  ! Новое расстояние до HP
		
			if (R_HP > Surf_Get_HP(gl_x(yzel), KORD(2), KORD(3))) then
				R_HP = R_HP - dr
			else
				R_HP = R_HP + dr
			end if
			
			
			
			! Блокируем схлопывание контакта к оси
			if(R_HP < 20.0_8) then
				R_HP = 20.0_8
			end if
			
            
            xx = gl_x(gl_RAY_B(par_n_HP, par_m_BC, k))              ! Отталкиваемся от x - координаты крайней точки B на гелиопаузе в этой плоскости (k)
            x = xx - (DBLE(j)/N2)**par_kk31 * (xx - par_R_LEFT)
            
            ! BS     Нужно взять положение BS из её положения на крайнем луче A
            yzel = gl_RAY_A(par_n_BS, size(gl_RAY_A(1, :, 1)), k)
			KORD(1) = 0.0_8; KORD(2) = gl_y(yzel); KORD(3) = gl_z(yzel)
            R_BS = norm2(KORD)  ! Новое расстояние до BS
			
			
            
            do i = 1, N1
                yzel = gl_RAY_O(i, j, k)
                
                if (i <= par_n_BS - par_n_HP + 1) then
                    r = R_HP + (i - 1) * (R_BS - R_HP)/(par_n_BS - par_n_HP)
                else
                    r = R_BS + (DBLE(i - (par_n_BS - par_n_HP + 1))/(N1 - (par_n_BS - par_n_HP + 1) ))**(0.55 * par_kk2) * (par_R_END - R_BS)
                end if


                gl_x(yzel) = x
                gl_y(yzel) = r * cos(phi)
                gl_z(yzel) = r * sin(phi)
                
            end do
        end do
    end do
    
    N3 = size(gl_RAY_K(1, 1, :))
    N2 = size(gl_RAY_K(1, :, 1))
    N1 = size(gl_RAY_K(:, 1, 1))

    ! Цикл движения точек на лучах K  ************************************************************
    do k = 1, N3
        do j = 1, N2
            
             ! Вычисляем координаты текущего луча в пространстве
            the = par_pi_8/2.0 + par_triple_point + (N2 - j + 1) * (par_pi_8/2.0 - par_triple_point)/(N2)
            phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
            
            if (k /= 1 .and. j == 1) CYCLE
            
            yzel = gl_RAY_K(par_n_TS, j, k)
            
            KORD(1) = gl_x(yzel); KORD(2) = gl_y(yzel); KORD(3) = gl_z(yzel) 
			R_TS = norm2(KORD)  ! Новое расстояние до TS
			
			if (R_TS > Surf_Get_TS(KORD(1), KORD(2), KORD(3))) then
				R_TS = R_TS - dr
			else
				R_TS = R_TS + dr
			end if
            
            do i = 1, N1

                if (i == 1) CYCLE
                
                ! Вычисляем координаты точки на луче
                yzel = gl_RAY_K(i, j, k)
				
				if (i <= par_n_IB) then  ! NEW
                    r =  par_R0 + (par_R_inner - par_R0) * (DBLE(i)/(par_n_IB))**par_kk1
                else 
                    r =  par_R_inner + (R_TS - par_R_inner) * (DBLE(i - par_n_IB)/(par_n_TS - par_n_IB))**par_kk12
				end if
				
				if (i == par_n_TS - 1) then
					r = R_TS - 0.5              
				end if
					
                !r =  par_R0 + (R_TS - par_R0) * (REAL(i, KIND = 4)/par_n_TS)**par_kk1


                ! Записываем новые координаты
                gl_x(yzel) = r * cos(the)
                gl_y(yzel) = r * sin(the) * cos(phi)
                gl_z(yzel) = r * sin(the) * sin(phi)
                
                
            end do
        end do
    end do
    
    
     N3 = size(gl_RAY_D(1, 1, :))
    N2 = size(gl_RAY_D(1, :, 1))
    N1 = size(gl_RAY_D(:, 1, 1))

    ! Цикл движения точек на лучах D ************************************************************
	
	
	
    do k = 1, N3
        do j = 1, N2
            
            ! Вычисляем координаты текущего луча в пространстве
                phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
				
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
				
				
            do i = 1, N1

                if (i == 1) CYCLE

                if (k /= 1 .and. j == 1) CYCLE
				

                yzel = gl_RAY_D(i, j, k)
                
                
                x = xx + (DBLE(i - 1)/(N1 - 1))**par_kk3 * (par_R_LEFT - xx)
				
				
                ! Записываем новые координаты
                gl_x(yzel) = x
                gl_y(yzel) = r * cos(phi)
                gl_z(yzel) = r * sin(phi)
                
            end do
        end do
    end do
    
    N3 = size(gl_RAY_E(1, 1, :))
    N2 = size(gl_RAY_E(1, :, 1))
    N1 = size(gl_RAY_E(:, 1, 1))
    
    ! Цикл движения точек на лучах E  ************************************************************
    do k = 1, N3
        do j = 1, N2
            do i = 1, N1

                if (i == 1) CYCLE

                if (i == N1) CYCLE

                x = gl_x(gl_RAY_E(1, j, k))
                y = gl_y(gl_RAY_E(1, j, k))
                z = gl_z(gl_RAY_E(1, j, k))

                x2 = gl_x(gl_RAY_O(1, j, k))
                y2 = gl_y(gl_RAY_O(1, j, k))
                z2 = gl_z(gl_RAY_O(1, j, k))

                yzel = gl_RAY_E(i, j, k)

                gl_x(yzel) = x + (x2 - x) * (i - 1)/(N1 - 1)
                gl_y(yzel) = y + (y2 - y) * (i - 1)/(N1 - 1)
                gl_z(yzel) = z + (z2 - z) * (i - 1)/(N1 - 1)
            end do
        end do
	end do
	
	end subroutine Surf_Set_surf
	
end module
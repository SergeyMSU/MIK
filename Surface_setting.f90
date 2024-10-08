
	
module Surface_setting  
	! ������ ��� ���������� ��������� ������������ �� ����� �����
    ! ��� ����������� ���������������� �� �� ������
	! ������ �������� � ��������
	USE GEO_PARAM
	USE STORAGE
	USE My_func
	USE MY_CUDA_smooth
	
	integer(4) :: Surf_l_phi   ! ���������� �������� ��������� 
	integer(4) :: Surf_m_A      ! ���������� ����� A � ���������
    integer(4) :: Surf_m_BC      ! ���������� ����� B/C � ���������
    integer(4) :: Surf_m_O      ! ���������� ����� O � ���������
    integer(4) :: Surf_m_K      ! ���������� ����� K � ���������
    real(8) :: Surf_triple_point
	real(8) :: Surf_al1
	
	
	real(8), allocatable :: Surf_phi(:)   ! ���� ����� ���������� ��������
	
	real(8), allocatable :: Surf_TS_the(:, :)   ! (:, :)   (j, k)  -  ���� � ���� ���������
	real(8), allocatable :: Surf_TS_r(:, :)     ! (:, :)   (j, k)  -  ���������� �� ������
	
	real(8), allocatable :: Surf_BS_the(:, :)   ! (:, :)   (j, k)  -  ���� � ���� ���������
	real(8), allocatable :: Surf_BS_r(:, :)     ! (:, :)   (j, k)  -  ���������� �� ������
	
	! ��� A - �����
	real(8), allocatable :: Surf_HP_the(:, :)   ! (:, :)   (j, k)  -  ���� � ���� ���������
	real(8), allocatable :: Surf_HP_r_A(:, :)     ! (:, :)   (j, k)  -  ���������� �� ������
	
	! ��� ��������� �����
	real(8), allocatable :: Surf_HP_x(:, :)   ! (:, :)   (j, k)  -  � ���������� � ���� ���������
	real(8), allocatable :: Surf_HP_r(:, :)     ! (:, :)   (j, k)  -  ���������� �� ��� ���������
	
	contains
	
	subroutine Surf_Save_bin(num)
		! ��������� ��������� ������������ �������� ����� (������� ��������� � ������ Storage)
		use STORAGE
		use GEO_PARAM
		implicit none
		integer, intent(in) :: num
		character(len=5) :: name
		integer :: i, j, k, N1, N2, N3, yzel
		real(8):: kord(3)
		
		write(unit=name, fmt='(i5.5)') num
		
		open(1, file = par_NAME // "Surf_Save" // name // ".bin", FORM = 'BINARY')
		
		
		write(1) par_l_phi, par_m_A, par_m_BC, par_m_O, par_m_K, par_triple_point, par_al1
		
		do k = 1, par_l_phi
			
			yzel = gl_RAY_A(5, 3, k)
			kord = (/gl_x(yzel), gl_y(yzel), gl_z(yzel)/)
				
			write(1) polar_angle(kord(2), kord(3))
			
			! A - ����
			do j = 1, par_m_A
				
				! TS
				yzel = gl_RAY_A(par_n_TS, j, k)
				kord = (/gl_x(yzel), gl_y(yzel), gl_z(yzel)/)
				
				write(1) polar_angle(kord(1), sqrt(kord(2)**2 + kord(3)**2)), norm2(kord)
				
				
				! ����������
				yzel = gl_RAY_A(par_n_HP, j, k)
				kord = (/gl_x(yzel), gl_y(yzel), gl_z(yzel)/)
				write(1) norm2(kord)
				
				! BS
				yzel = gl_RAY_A(par_n_BS, j, k)
				kord = (/gl_x(yzel), gl_y(yzel), gl_z(yzel)/)
				write(1) norm2(kord)
				
			end do
			
			! B - ����
			do j = 1, par_m_BC
				
				! the = (j - 1) * par_pi_8/2.0/(par_m_A - 1)
				
				! TS
				yzel = gl_RAY_B(par_n_TS, j, k)
				kord = (/gl_x(yzel), gl_y(yzel), gl_z(yzel)/)
				write(1) polar_angle(kord(1), sqrt(kord(2)**2 + kord(3)**2)),  norm2(kord)
				
				
				! ����������
				yzel = gl_RAY_B(par_n_HP, j, k)
				kord = (/ gl_x(yzel), gl_y(yzel), gl_z(yzel)/)
				write(1) kord(1), sqrt(kord(2)**2 + kord(3)**2)
				
			end do
			
			
			! O - ����
			do j = 1, par_m_O
				
				! ����������
				yzel = gl_RAY_O(1, j, k)
				kord = (/ gl_x(yzel), gl_y(yzel), gl_z(yzel)/)
				write(1) kord(1), sqrt(kord(2)**2 + kord(3)**2)
				
			end do
			
			
			! K - ����
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
		
		inquire(file=par_NAME//"Surf_Save" // name // ".bin", exist=exists)
		
		if (exists == .False.) then
			pause "net faila!!!  8idjgeye0-0987fhjdkeye "
			STOP "net faila!!!"
		end if
		
		
		
		
		open(1, file = par_NAME//"Surf_Save" // name // ".bin", FORM = 'BINARY', ACTION = "READ")
		
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
			
			! K - ����
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
	! ����������� ����������� � ������������ � ������������ �����������
	use STORAGE
    use GEO_PARAM
	use My_func
	! Variables
	real(8), intent(in) :: dr  ! ��� ��������
	real(8) :: R_TS, proect, vel(3), R_HP, R_BS, KORD(3), dist, ddt, ER(3)
    integer :: yzel, N1, N2, N3, i, j, k, yzel2, ij
    real(8) :: the, phi, r, x, y, z, rr, xx, x2, y2, z2, rrr, r1, r2, r3, r4, rd, kk13, kk14, kk15, dk13, ddr, rr0, dddr
    real(8) :: yy1, yy2, yy3
    real(8) :: y3, z3
	
	! Body of Set_surf
	
	N3 = size(gl_RAY_A(1, 1, :))
    N2 = size(gl_RAY_A(1, :, 1))
    N1 = size(gl_RAY_A(:, 1, 1))
	
	! ���� �������� ����� �� ����� �  ************************************************************
    do k = 1, N3
        do j = 1, N2
            
            if (k /= 1 .and. j == 1) then
                    CYCLE
            end if
            
            ! ��������� ���������� �������� ���� � ������������
            the = par_pi_8/2.0 * (DBLE(j - 1.0)/(N2 - 1.0))
            phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
            
            ! TS
            yzel = gl_RAY_A(par_n_TS, j, k)
            KORD(1) = gl_x(yzel); KORD(2) = gl_y(yzel); KORD(3) = gl_z(yzel) 
			R_TS = norm2(KORD)  ! ����� ���������� �� TS
			
			if (R_TS > Surf_Get_TS(KORD(1), KORD(2), KORD(3))) then
				R_TS = R_TS - dr
			else
				R_TS = R_TS + dr
			end if
			
			
			
            ! HP
            yzel = gl_RAY_A(par_n_HP, j, k)
			KORD(1) = gl_x(yzel); KORD(2) = gl_y(yzel); KORD(3) = gl_z(yzel) 
            R_HP = norm2(KORD)  ! ����� ���������� �� HP
			
			if (R_HP > Surf_Get_HP(KORD(1), KORD(2), KORD(3))) then
				R_HP = R_HP - dr
			else
				R_HP = R_HP + dr
			end if
			
			
            ! BS
            yzel = gl_RAY_A(par_n_BS, j, k)
            
			KORD(1) = gl_x(yzel); KORD(2) = gl_y(yzel); KORD(3) = gl_z(yzel) 
            R_BS = norm2(KORD)  ! ����� ���������� �� BS
			
			if (R_BS > Surf_Get_BS(KORD(1), KORD(2), KORD(3))) then
				R_BS = R_BS - dr
			else
				R_BS = R_BS + dr
			end if
			
            
            ! ����� ������� ���� ���������� ��������� �����, ����� ��, ��� � ��� ���������� �����
            do i = 1, N1
                
               if (i == 1) CYCLE

                yzel = gl_RAY_A(i, j, k)
                ! ��������� ���������� ����� �� ����
				
				!kk13 = par_kk13 * (1.0 + dabs(the)/18.0)! * (par_pi_8/2.0 - dabs(the))/(par_pi_8/2.0)  +  (par_kk13 - 0.2) * (dabs(the))/(par_pi_8/2.0)
				!dk13 = 0.1 + (dabs(the)/(par_pi_8/2.0)/2.1)**2
				
				dk13 = (R_TS + (R_HP - R_TS) * sgushenie_2(DBLE(par_n_HP - par_n_TS)/(par_n_HP - par_n_TS), par_kk14)) - &
						(R_TS + (R_HP - R_TS) * sgushenie_2(DBLE(par_n_HP - 1 - par_n_TS)/(par_n_HP - par_n_TS), par_kk14))

                ! �� TS
				
				r = Setka_A(i, R_TS, R_HP, R_BS, the, dk13, par_R0, par_R_inner, par_n_IB, par_kk1, par_kk12, &
		par_n_TS, par_n_HP, par_n_BS, par_kk14, par_R_END, par_kk2, par_n_END)

                ! ���������� ����� ����������
                gl_x(yzel) = r * cos(the)
                gl_y(yzel) = r * sin(the) * cos(phi)
                gl_z(yzel) = r * sin(the) * sin(phi)
				
            end do
        end do
	end do
    
		
    N3 = size(gl_RAY_B(1, 1, :))
    N2 = size(gl_RAY_B(1, :, 1))
    N1 = size(gl_RAY_B(:, 1, 1))

    ! ���� �������� ����� �� ����� B  ************************************************************
    do k = 1, N3
        do j = 1, N2
            
            ! ��������� ���������� �������� ���� � ������������
            the = par_pi_8/2.0 + (j) * par_triple_point/(N2)
			the2 = par_pi_8/2.0 + (j) * par_triple_point_2/(N2)
            phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
            
            ! TS
            yzel = gl_RAY_B(par_n_TS, j, k)
            
			ER(1) = cos(the); ER(2) = sin(the) * cos(phi); ER(3) = sin(the) * sin(phi)
			KORD(1) = gl_x(yzel); KORD(2) = gl_y(yzel); KORD(3) = gl_z(yzel) 
            R_TS = norm2(KORD)  ! ����� ���������� �� TS
			
			if (R_TS > Surf_Get_TS(KORD(1), KORD(2), KORD(3))) then
				R_TS = R_TS - dr
			else
				R_TS = R_TS + dr
			end if
            
            ! HP
            yzel = gl_RAY_B(par_n_HP, j, k)
			KORD(1) = gl_x(yzel); KORD(2) = gl_y(yzel); KORD(3) = gl_z(yzel) 
            !R_HP = norm2(KORD)  ! ����� ���������� �� HP
			R_HP = norm2(KORD - (R_TS * ER))  ! ����� ���������� �� HP �� ���� TS
			
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
                ! ��������� ���������� ����� �� ����

                ! �� TS
				if (i <= par_n_IB) then  ! NEW
					if(i == 2) then
						r =  par_R0 - (par_R_inner - par_R0) * (DBLE(3 - 2)/(par_n_IB - 2))**par_kk1
						if(r < 0.0) then
							print*, "Error iouihjgfdcydygy  ", r
							STOP
						end if
					else
                        r =  par_R0 + (par_R_inner - par_R0) * (DBLE(i - 2)/(par_n_IB - 2))**par_kk1
					end if
                else if (i <= par_n_TS) then  ! �� ���������� = R_TS
                    r =  par_R_inner + (R_TS - par_R_inner) * sgushenie_3( (DBLE(i - par_n_IB)/(par_n_TS - par_n_IB)) , par_kk12)
				else if (i <= par_n_HP) then  ! �� ���������� = par_R_character * 1.3
                    !r = R_TS + (i - par_n_TS) * (R_HP - R_TS) /(par_n_HP - par_n_TS)
					
					r1 = R_TS
                    !r =  (i - par_n_TS) * (R_HP) /(par_n_HP - par_n_TS)
                    r =  R_HP * sgushenie_2(DBLE(i - par_n_TS)/(par_n_HP - par_n_TS), par_kk14)
					gl_x(yzel) = r1 * cos(the) + r * cos(the2)
					gl_y(yzel) = r1 * sin(the) * cos(phi) + r * sin(the2) * cos(phi)
					gl_z(yzel) = r1 * sin(the) * sin(phi) + r * sin(the2) * sin(phi)
					CYCLE
				end if
				
				!if (i == par_n_TS - 1) then
				!	r = R_TS - 0.5               
				!end if

                ! ���������� ����� ����������
                gl_x(yzel) = r * cos(the)
                gl_y(yzel) = r * sin(the) * cos(phi)
                gl_z(yzel) = r * sin(the) * sin(phi)
                
               
            end do
        end do
	end do
    
	
    N3 = size(gl_RAY_C(1, 1, :))
    N2 = size(gl_RAY_C(1, :, 1))
    N1 = size(gl_RAY_C(:, 1, 1))

    ! ���� �������� ����� �� ����� C  ************************************************************
    do k = 1, N3
		yzel = gl_RAY_A(par_n_HP, size(gl_RAY_A(1, :, 1)), k)
		ER(1) = 0.0_8; ER(2) = gl_y(yzel); ER(3) = gl_z(yzel)
		rr0 = norm2(ER)  ! ����� ���������� �� BS

        do j = 1, N2
            
            ! ��������� ���������� �������� ���� � ������������
                phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
                
                ! ��������� ���������� ����� �� ����
                x = gl_x(gl_RAY_B(par_n_HP, j, k))
                y = gl_y(gl_RAY_B(par_n_HP, j, k))
                z = gl_z(gl_RAY_B(par_n_HP, j, k))
                rr = (y**2 + z**2)**(0.5)
				
				xx = gl_x(gl_RAY_B(par_n_HP, size(gl_RAY_B(1, :, 1)), k))
				x2 = gl_x(gl_RAY_B(par_n_HP, 1, k))
				
				kk13 = (par_kk13 - 0.4) * dabs(xx - x)/dabs(xx - x2)  +  (1.0) * dabs(x - x2)/dabs(xx - x2)
				
				
				y = gl_y(gl_RAY_B(par_n_HP - 1, j, k))
                z = gl_z(gl_RAY_B(par_n_HP - 1, j, k))
                rd = (y**2 + z**2)**(0.5)
				
				!rr = rr + (rr - rd)
                
                ! BS     ����� ����� ��������� BS �� � ��������� �� ������� ���� A
                yzel = gl_RAY_A(par_n_BS, size(gl_RAY_A(1, :, 1)), k)
				ER(1) = 0.0_8; ER(2) = gl_y(yzel); ER(3) = gl_z(yzel)
                R_BS = norm2(ER)  ! ����� ���������� �� BS


				yzel = gl_RAY_A(par_n_BS - 1, size(gl_RAY_A(1, :, 1)), k)
				ER(1) = 0.0_8; ER(2) = gl_y(yzel); ER(3) = gl_z(yzel)
                ddr = R_BS - norm2(ER)  ! ����� ���������� �� BS
				
            
            do i = 1, N1
                
                if(i == 1) CYCLE
                
                yzel = gl_RAY_C(i, j, k)
				
				dk13 = sqrt((gl_x(gl_RAY_B(size(gl_RAY_B(:,1,1)), j, k)) - gl_x(gl_RAY_B(size(gl_RAY_B(:,1,1)) - 1, j, k)))**2 + &
					(gl_y(gl_RAY_B(size(gl_RAY_B(:,1,1)), j, k)) - gl_y(gl_RAY_B(size(gl_RAY_B(:,1,1)) - 1, j, k)))**2 + &
					(gl_z(gl_RAY_B(size(gl_RAY_B(:,1,1)), j, k)) - gl_z(gl_RAY_B(size(gl_RAY_B(:,1,1)) - 1, j, k)))**2)

    !            if(i == 2) then
				!	r = rr
				!else 
				
				!if (i <= 11) then
    !                r = rr + dk13 * (i - 1)
				!else if (i <= par_n_BS - par_n_HP + 1) then
				!	rrr = rr + dk13 * (10)
				!	r1 = log((1.5 * dk13)/(R_BS - rr))/log(DBLE(1.0)/(par_n_BS - par_n_HP - 10))
    !                r = rrr + (R_BS - rrr) * (DBLE(i - 11)/(par_n_BS - par_n_HP - 10))**r1
    !            else
    !                r = R_BS + (DBLE(i - (par_n_BS - par_n_HP + 1))/(N1 - (par_n_BS - par_n_HP + 1) ))**(0.55 * par_kk2) * (par_R_END - R_BS)
				!end if
				
				r = Setka_C(i, R_BS, dk13, par_n_HP, par_n_BS, par_kk2, par_R_END, N1, rr, dabs(rr - rd), ddr, rr0)
                
                gl_x(yzel) = x
                gl_y(yzel) = r * cos(phi)
                gl_z(yzel) = r * sin(phi)
                

            end do
        end do
    end do

    N3 = size(gl_RAY_O(1, 1, :))
    N2 = size(gl_RAY_O(1, :, 1))
    N1 = size(gl_RAY_O(:, 1, 1))


    ! ���� �������� ����� �� ����� O   ************************************************************
    do k = 1, N3
        phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
        
		ij = size(gl_RAY_C(1, :, k))
		
		yzel = gl_RAY_C(1, ij, k)
		yy1 = sqrt(gl_y(yzel)**2 + gl_z(yzel)**2)
		
		yzel = gl_RAY_C(2, ij, k)
		yy2 = sqrt(gl_y(yzel)**2 + gl_z(yzel)**2)
		
		yzel = gl_RAY_C(size(gl_RAY_C(:, ij, k)), ij, k)
		yy3 = sqrt(gl_y(yzel)**2 + gl_z(yzel)**2)
		
		dk13 = 2.0 * (yy2 - yy1)/(yy3 - yy1)
		r1 = log(dk13)/log(1.0/(par_n_BS - par_n_HP - 9))

		yzel = gl_RAY_A(par_n_HP, size(gl_RAY_A(1, :, 1)), k)
		ER(1) = 0.0_8; ER(2) = gl_y(yzel); ER(3) = gl_z(yzel)
		rr0 = norm2(ER)  ! ����� ���������� �� BS
		
        do j = 1, N2
            
            yzel = gl_RAY_O(1, j, k)
			KORD(1) = 0.0_8; KORD(2) = gl_y(yzel); KORD(3) = gl_z(yzel)
            R_HP = norm2(KORD)  ! ����� ���������� �� HP
		
			if (R_HP > Surf_Get_HP(gl_x(yzel), KORD(2), KORD(3))) then
				R_HP = R_HP - dr
			else
				R_HP = R_HP + dr
			end if

			
			
			
			
			! ��������� ����������� �������� � ���
			if(R_HP < 20.0_8) then
				R_HP = 20.0_8
			end if
			
            
            xx = gl_x(gl_RAY_B(par_n_HP, par_m_BC, k))              ! ������������� �� x - ���������� ������� ����� B �� ���������� � ���� ��������� (k)
            x = xx - (DBLE(j)/N2)**par_kk31 * (xx - par_R_LEFT)
            
            ! BS     ����� ����� ��������� BS �� � ��������� �� ������� ���� A
            yzel = gl_RAY_A(par_n_BS, size(gl_RAY_A(1, :, 1)), k)
			KORD(1) = 0.0_8; KORD(2) = gl_y(yzel); KORD(3) = gl_z(yzel)
            R_BS = norm2(KORD)  ! ����� ���������� �� BS

			yzel = gl_RAY_A(par_n_BS - 1, size(gl_RAY_A(1, :, 1)), k)
			ER(1) = 0.0_8; ER(2) = gl_y(yzel); ER(3) = gl_z(yzel)
			ddr = R_BS - norm2(ER)  ! ����� ���������� �� BS

			
			yzel = gl_RAY_E(size(gl_RAY_E(:, 1, 1)) - 1, j, k)
			ER(1) = 0.0_8; ER(2) = gl_y(yzel); ER(3) = gl_z(yzel)
			dddr = max(min(dabs(R_HP - norm2(ER)), 10.0_8), 0.01_8)  ! ����� ���������� �� BS
			!dddr = 1.0

			
            do i = 1, N1
                yzel = gl_RAY_O(i, j, k)
                
				!if (i <= 10) then
    !                r = R_HP + (i - 1) * dk13 * (R_BS - R_HP)
				!else if (i <= par_n_BS - par_n_HP + 1) then
				!	rr = R_HP + 9 * dk13 * (R_BS - R_HP)
    !                r = rr + (R_BS - rr) * (DBLE(i - 10)/(par_n_BS - par_n_HP - 9))**r1
    !            else
    !                r = R_BS + (DBLE(i - (par_n_BS - par_n_HP + 1))/(N1 - (par_n_BS - par_n_HP + 1) ))**(0.55 * par_kk2) * (par_R_END - R_BS)
    !            end if
				
				r = Setka_O(i, r1, R_HP, R_BS, dk13, par_n_HP, par_n_BS, par_kk2, par_R_END, N1, ddr = ddr, rr0 = rr0, dddr = dddr)

                gl_x(yzel) = x
                gl_y(yzel) = r * cos(phi)
                gl_z(yzel) = r * sin(phi)
                
            end do
        end do
    end do
    
    N3 = size(gl_RAY_K(1, 1, :))
    N2 = size(gl_RAY_K(1, :, 1))
    N1 = size(gl_RAY_K(:, 1, 1))

    ! ���� �������� ����� �� ����� K  ************************************************************
    do k = 1, N3
        do j = 1, N2
            
             ! ��������� ���������� �������� ���� � ������������
            the = par_pi_8/2.0 + par_triple_point + (N2 - j + 1) * (par_pi_8/2.0 - par_triple_point)/(N2)
            phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
            
            if (k /= 1 .and. j == 1) CYCLE
            
            yzel = gl_RAY_K(par_n_TS, j, k)
            
            KORD(1) = gl_x(yzel); KORD(2) = gl_y(yzel); KORD(3) = gl_z(yzel) 
			R_TS = norm2(KORD)  ! ����� ���������� �� TS
			
			if (R_TS > Surf_Get_TS(KORD(1), KORD(2), KORD(3))) then
				R_TS = R_TS - dr
			else
				R_TS = R_TS + dr
			end if
            
            do i = 1, N1

                if (i == 1) CYCLE
                
                ! ��������� ���������� ����� �� ����
                yzel = gl_RAY_K(i, j, k)
				
				if (i <= par_n_IB) then  ! NEW
                    if(i == 2) then
						r =  par_R0 - (par_R_inner - par_R0) * (DBLE(3 - 2)/(par_n_IB - 2))**par_kk1
						if(r < 0.0) then
							print*, "Error iouihjgfdcydygy  ", r
							STOP
						end if
					else
                        r =  par_R0 + (par_R_inner - par_R0) * (DBLE(i - 2)/(par_n_IB - 2))**par_kk1
					end if
					!r =  par_R0 + (par_R_inner - par_R0) * (DBLE(i - 2)/(par_n_IB - 2))**par_kk1
                else 
                    r =  par_R_inner + (R_TS - par_R_inner) * sgushenie_3(DBLE(i - par_n_IB)/(par_n_TS - par_n_IB), par_kk12)
				end if
				
				!if (i == par_n_TS - 1) then
				!	r = R_TS - 0.5              
				!end if
					
                !r =  par_R0 + (R_TS - par_R0) * (REAL(i, KIND = 4)/par_n_TS)**par_kk1


                ! ���������� ����� ����������
                gl_x(yzel) = r * cos(the)
                gl_y(yzel) = r * sin(the) * cos(phi)
                gl_z(yzel) = r * sin(the) * sin(phi)
                
                
            end do
        end do
    end do
    
    
     N3 = size(gl_RAY_D(1, 1, :))
    N2 = size(gl_RAY_D(1, :, 1))
    N1 = size(gl_RAY_D(:, 1, 1))

    ! ���� �������� ����� �� ����� D ************************************************************
	
	
	
    do k = 1, N3

		y2 = gl_y(gl_RAY_O(1, size(gl_RAY_O(1, :, 1)), k))
		z2 = gl_z(gl_RAY_O(1, size(gl_RAY_O(1, :, 1)), k))

		y3 = gl_y(gl_RAY_C(1, 1, k))
		z3 = gl_z(gl_RAY_C(1, 1, k))

		r1 = sqrt(y3**2 + z3**2)/sqrt(y2**2 + z2**2)

        do j = 1, N2
            
            ! ��������� ���������� �������� ���� � ������������
                phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
				
				 ! ��������� ���������� ����� �� ����

                if (j < N2) then
                    xx = gl_x(gl_RAY_K(par_n_TS, j, k))
                    y = gl_y(gl_RAY_K(par_n_TS, j, k))
                    z = gl_z(gl_RAY_K(par_n_TS, j, k))
                else
                    xx = gl_x(gl_RAY_B(par_n_TS, par_m_BC, k))
                    y = gl_y(gl_RAY_B(par_n_TS, par_m_BC, k))
                    z = gl_z(gl_RAY_B(par_n_TS, par_m_BC, k))
				end if

				
            do i = 1, N1

				r = (i - N1)/(1.0 - N1) * sqrt(y**2 + z**2) + (i - 1.0)/(N1 - 1.0) * sqrt(y**2 + z**2)/r1

                if (i == 1) CYCLE

                if (k /= 1 .and. j == 1) CYCLE
				

                yzel = gl_RAY_D(i, j, k)
                
                
                x = xx + (DBLE(i - 1)/(N1 - 1))**par_kk3 * (par_R_LEFT - xx)
				
				
                ! ���������� ����� ����������
                gl_x(yzel) = x
                gl_y(yzel) = r * cos(phi)
                gl_z(yzel) = r * sin(phi)
                
            end do
        end do
    end do
    
    N3 = size(gl_RAY_E(1, 1, :))
    N2 = size(gl_RAY_E(1, :, 1))
    N1 = size(gl_RAY_E(:, 1, 1))
    
    ! ���� �������� ����� �� ����� E  ************************************************************
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

                gl_x(yzel) = x + (x2 - x) * sgushenie_2(DBLE(i - 1)/(N1 - 1), par_kk14)
                gl_y(yzel) = y + (y2 - y) * sgushenie_2(DBLE(i - 1)/(N1 - 1), par_kk14)
                gl_z(yzel) = z + (z2 - z) * sgushenie_2(DBLE(i - 1)/(N1 - 1), par_kk14)
            end do
        end do
	end do
	
	end subroutine Surf_Set_surf
	
end module
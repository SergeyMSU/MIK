!******7**10**********************************************************70
!**                                                                  **
!**    chlld   for                                                   **
!**                                                                  **
!******7**10**********************************************************70
!**    shema HLLD(modification) for 3D MHD                           **
!******7**10**********************************************************70

	
	attributes(global) subroutine Cuda_Move_all_1(now)  ! ������������� ��������� �� � - �����
	use MY_CUDA, gl_TS => dev_gl_TS, gl_Gran_neighbour => dev_gl_Gran_neighbour, gl_Gran_normal2 => dev_gl_Gran_normal2, &
		gl_Cell_par => dev_gl_Cell_par, gl_all_Gran => dev_gl_all_Gran, gl_Point_num => dev_gl_Point_num, gl_x2 => dev_gl_x2, &
		gl_y2 => dev_gl_y2, gl_z2 => dev_gl_z2, gl_Gran_center2 => dev_gl_Gran_center2, gl_Vx => dev_gl_Vx, &
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS, &
		gl_RAY_A => dev_gl_RAY_A, gl_RAY_B => dev_gl_RAY_B, gl_RAY_C => dev_gl_RAY_C, gl_RAY_O => dev_gl_RAY_O, &
		gl_RAY_K => dev_gl_RAY_K, gl_RAY_D => dev_gl_RAY_D, gl_RAY_E => dev_gl_RAY_E, norm2 => dev_norm2, par_n_TS => dev_par_n_TS, &
		par_n_HP => dev_par_n_HP, par_n_BS => dev_par_n_BS
	use GEO_PARAM
	use cudafor
	
	implicit none
	integer, intent(in) :: now
	
	real(8) :: Time, dist, r, rr, r1, r2, r3, r4, ddt
    real(8) :: vel(3), Ak(3), Bk(3), Ck(3), Dk(3), Ek(3)
    integer :: i, j, k
	
	integer(4):: N2, N3
	integer(4) :: yzel, yzel2
	
	  
	Time = time_step
	
	k = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! ����� ������
	j = blockDim%y * (blockIdx%y - 1) + threadIdx%y   ! ����� ������
	
	N3 = size(gl_RAY_A(1, 1, :))
    N2 = size(gl_RAY_A(1, :, 1))
	
	if(k > N3 .or. j > N2) return
	
	
	if (j == 1) then
                    return
	end if
	
	ddt = Time/0.000127
			
	! ������� �����
	yzel = gl_RAY_A(par_n_TS, j, k)
		Ak(1) = gl_x2(yzel, now); Ak(2) = gl_y2(yzel, now); Ak(3) = gl_z2(yzel, now)
			
	if (j < N2) then
		yzel2 = gl_RAY_A(par_n_TS, j + 1, k)
	else
		yzel2 = gl_RAY_B(par_n_TS, 1, k)
	end if
	Bk(1) = gl_x2(yzel2, now); Bk(2) = gl_y2(yzel2, now); Bk(3) = gl_z2(yzel2, now)
			
	if (k > 1) then
		yzel2 = gl_RAY_A(par_n_TS, j, k - 1)
	else
		yzel2 = gl_RAY_A(par_n_TS, j, N3)
	end if
	Ck(1) = gl_x2(yzel2, now); Ck(2) = gl_y2(yzel2, now); Ck(3) = gl_z2(yzel2, now)
			
	if (j > 1) then
		yzel2 = gl_RAY_A(par_n_TS, j - 1, k)
	else
		yzel2 = gl_RAY_A(par_n_TS, j, k + ceiling(1.0 * N3/2))
	end if
	Dk(1) = gl_x2(yzel2, now); Dk(2) = gl_y2(yzel2, now); Dk(3) = gl_z2(yzel2, now)
			
	if(k < N3) then
		yzel2 = gl_RAY_A(par_n_TS, j, k + 1)
	else
		yzel2 = gl_RAY_A(par_n_TS, j, 1)
	end if
	Ek(1) = gl_x2(yzel2, now); Ek(2) = gl_y2(yzel2, now); Ek(3) = gl_z2(yzel2, now)
			
	if (gl_Point_num(yzel) > 0) then
		vel = gl_Point_num(yzel) * par_nat_TS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time * ddt
	else
		vel = par_nat_TS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time * ddt
	end if
		
			
	gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
	gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
	gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
	
	
	 ! �������
			
	yzel = gl_RAY_A(par_n_HP, j, k)
	Ak(1) = gl_x2(yzel, now); Ak(2) = gl_y2(yzel, now); Ak(3) = gl_z2(yzel, now)
	r = norm2(Ak)
			
	if (j < N2) then
		yzel2 = gl_RAY_A(par_n_HP, j + 1, k)
	else
		yzel2 = gl_RAY_B(par_n_HP, 1, k)
	end if
	Bk(1) = gl_x2(yzel2, now); Bk(2) = gl_y2(yzel2, now); Bk(3) = gl_z2(yzel2, now)
	r1 = norm2(Bk)
			
	if (k > 1) then
		yzel2 = gl_RAY_A(par_n_HP, j, k - 1)
	else
		yzel2 = gl_RAY_A(par_n_HP, j, N3)
	end if
			
	Ck(1) = gl_x2(yzel2, now); Ck(2) = gl_y2(yzel2, now); Ck(3) = gl_z2(yzel2, now)
	r2 = norm2(Ck)
			
	if (j > 1) then
		yzel2 = gl_RAY_A(par_n_HP, j - 1, k)
	else
		yzel2 = gl_RAY_A(par_n_HP, j, k + ceiling(1.0 * N3/2))
	end if
	Dk(1) = gl_x2(yzel2, now); Dk(2) = gl_y2(yzel2, now); Dk(3) = gl_z2(yzel2, now)
	r3 = norm2(Dk)
			
	if(k < N3) then
		yzel2 = gl_RAY_A(par_n_HP, j, k + 1)
	else
		yzel2 = gl_RAY_A(par_n_HP, j, 1)
	end if
	Ek(1) = gl_x2(yzel2, now); Ek(2) = gl_y2(yzel2, now); Ek(3) = gl_z2(yzel2, now)
	r4 = norm2(Ek)
	
	rr = (r1 + r2 + r3 + r4)/4.0
	!dist = sqrt( (Dk(1) - Ak(1))**2 + (Dk(2) - Ak(2))**2 + (Dk(3) - Ak(3))**2 )
	
	!dist = max(dist, 1.0_8)
			
	if (gl_Point_num(yzel) > 0) then
		!vel = gl_Point_num(yzel) * par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time/dist
		vel = par_nat_HP * 0.035 * gl_Point_num(yzel) * ((Ak/r) * (rr - r)) * ddt
	else
		!vel = par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time/dist
		vel = par_nat_HP * 0.035 * ((Ak/r) * (rr - r)) * ddt
	end if
			
			
	gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
	gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
	gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
			
			
	! ������� ������� �����
			
	yzel = gl_RAY_A(par_n_BS, j, k)
	Ak(1) = gl_x2(yzel, now); Ak(2) = gl_y2(yzel, now); Ak(3) = gl_z2(yzel, now)
			
	if (j < N2) then
		yzel2 = gl_RAY_A(par_n_BS, j + 1, k)
	else
		return !yzel2 = yzel                                                    ! ��� ��� ���� ��������������� ��� ������ � ������� �����
	end if
	Bk(1) = gl_x2(yzel2, now); Bk(2) = gl_y2(yzel2, now); Bk(3) = gl_z2(yzel2, now)
			
	if (k > 1) then
		yzel2 = gl_RAY_A(par_n_BS, j, k - 1)
	else
		yzel2 = gl_RAY_A(par_n_BS, j, N3)
	end if
	Ck(1) = gl_x2(yzel2, now); Ck(2) = gl_y2(yzel2, now); Ck(3) = gl_z2(yzel2, now)
			
	if (j > 1) then
		yzel2 = gl_RAY_A(par_n_BS, j - 1, k)
	else
		yzel2 = gl_RAY_A(par_n_BS, j, k + ceiling(1.0 * N3/2))
	end if
	Dk(1) = gl_x2(yzel2, now); Dk(2) = gl_y2(yzel2, now); Dk(3) = gl_z2(yzel2, now)
			
	if(k < N3) then
		yzel2 = gl_RAY_A(par_n_BS, j, k + 1)
	else
		yzel2 = gl_RAY_A(par_n_BS, j, 1)
	end if
	Ek(1) = gl_x2(yzel2, now); Ek(2) = gl_y2(yzel2, now); Ek(3) = gl_z2(yzel2, now)
	
	dist = (sqrt( (Dk(1) - Ak(1))**2 + (Dk(2) - Ak(2))**2 + (Dk(3) - Ak(3))**2 ) + &
				sqrt( (Ek(1) - Ak(1))**2 + (Ek(2) - Ak(2))**2 + (Ek(3) - Ak(3))**2 ))/2.0
	
	dist = max(dist, 1.0_8)
			
	if (gl_Point_num(yzel) > 0) then
		vel = gl_Point_num(yzel) * par_nat_BS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time/dist * ddt
	else
		vel = par_nat_BS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time/dist * ddt
	end if
			
			
	gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
	gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
	gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
	
	end subroutine Cuda_Move_all_1
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	attributes(global) subroutine Cuda_Move_all_2(now)   ! ������������� ��������� �� 
	use MY_CUDA, gl_TS => dev_gl_TS, gl_Gran_neighbour => dev_gl_Gran_neighbour, gl_Gran_normal2 => dev_gl_Gran_normal2, &
		gl_Cell_par => dev_gl_Cell_par, gl_all_Gran => dev_gl_all_Gran, gl_Point_num => dev_gl_Point_num, gl_x2 => dev_gl_x2, &
		gl_y2 => dev_gl_y2, gl_z2 => dev_gl_z2, gl_Gran_center2 => dev_gl_Gran_center2, gl_Vx => dev_gl_Vx, &
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS, &
		gl_RAY_A => dev_gl_RAY_A, gl_RAY_B => dev_gl_RAY_B, gl_RAY_C => dev_gl_RAY_C, gl_RAY_O => dev_gl_RAY_O, &
		gl_RAY_K => dev_gl_RAY_K, gl_RAY_D => dev_gl_RAY_D, gl_RAY_E => dev_gl_RAY_E, norm2 => dev_norm2, par_n_TS => dev_par_n_TS, &
		par_n_HP => dev_par_n_HP, par_n_BS => dev_par_n_BS
	use GEO_PARAM
	use cudafor
	
	implicit none
	integer, intent(in) :: now
	
	real(8) :: Time, dist, r, rr, r1, r2, r3, r4, ddt
    real(8) :: vel(3), Ak(3), Bk(3), Ck(3), Dk(3), Ek(3)
    integer :: i, j, k
	
	integer(4):: N2, N3
	integer(4) :: yzel, yzel2
	
	N3 = size(gl_RAY_B(1, 1, :))
    N2 = size(gl_RAY_B(1, :, 1))
	  
	Time = time_step
	
	k = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! ����� ������
	j = blockDim%y * (blockIdx%y - 1) + threadIdx%y   ! ����� ������
	
	if(k > N3 .or. j > N2) return
	
	ddt = Time/0.000127
	! ������� �����
			
			yzel = gl_RAY_B(par_n_TS, j, k)
			Ak(1) = gl_x2(yzel, now); Ak(2) = gl_y2(yzel, now); Ak(3) = gl_z2(yzel, now)
			! Ak = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/)
			
			if (j < N2) then
			    yzel2 = gl_RAY_B(par_n_TS, j + 1, k)
			else
				yzel2 = gl_RAY_K(par_n_TS, size(gl_RAY_K(par_n_TS, :, k)), k)
			end if
			Bk(1) = gl_x2(yzel2, now); Bk(2) = gl_y2(yzel2, now); Bk(3) = gl_z2(yzel2, now)
			! Bk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if (k > 1) then
			    yzel2 = gl_RAY_B(par_n_TS, j, k - 1)
			else
			    yzel2 = gl_RAY_B(par_n_TS, j, N3)
			end if
			Ck(1) = gl_x2(yzel2, now); Ck(2) = gl_y2(yzel2, now); Ck(3) = gl_z2(yzel2, now)
			! Ck = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if (j > 1) then
			    yzel2 = gl_RAY_B(par_n_TS, j - 1, k)
			else
				yzel2 = gl_RAY_A(par_n_TS, size( gl_RAY_A(par_n_TS, :, k)), k)
			end if
			Dk(1) = gl_x2(yzel2, now); Dk(2) = gl_y2(yzel2, now); Dk(3) = gl_z2(yzel2, now)
			! Dk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if(k < N3) then
			    yzel2 = gl_RAY_B(par_n_TS, j, k + 1)
			else
				yzel2 = gl_RAY_B(par_n_TS, j, 1)
			end if
			Ek(1) = gl_x2(yzel2, now); Ek(2) = gl_y2(yzel2, now); Ek(3) = gl_z2(yzel2, now)
			! Ek = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			
			if (gl_Point_num(yzel) > 0) then
			    vel = gl_Point_num(yzel) * par_nat_TS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time * ddt
			else
				vel = par_nat_TS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time * ddt
			end if
			
			
			gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
			gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
			gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
			
			!return
			! �������
			
			yzel = gl_RAY_B(par_n_HP, j, k)
			Ak(1) = gl_x2(yzel, now); Ak(2) = gl_y2(yzel, now); Ak(3) = gl_z2(yzel, now)
			! Ak = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/)
			r = sqrt(Ak(2)**2 + Ak(3)**2)
			
			if (j < N2) then
			    yzel2 = gl_RAY_B(par_n_HP, j + 1, k)
			else
				!yzel2 = gl_RAY_O(1, 1, k)
				return
			end if
			Bk(1) = gl_x2(yzel2, now); Bk(2) = gl_y2(yzel2, now); Bk(3) = gl_z2(yzel2, now)
			! Bk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			r1 = sqrt(Bk(2)**2 + Bk(3)**2)
			
			if (k > 1) then
			    yzel2 = gl_RAY_B(par_n_HP, j, k - 1)
			else
			    yzel2 = gl_RAY_B(par_n_HP, j, N3)
			end if
			Ck(1) = gl_x2(yzel2, now); Ck(2) = gl_y2(yzel2, now); Ck(3) = gl_z2(yzel2, now)
			! Ck = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			r2 = sqrt(Ck(2)**2 + Ck(3)**2)
			
			if (j > 1) then
			    yzel2 = gl_RAY_B(par_n_HP, j - 1, k)
			else
				yzel2 = gl_RAY_A(par_n_HP, size( gl_RAY_A(par_n_HP, :, k)), k)
			end if
			Dk(1) = gl_x2(yzel2, now); Dk(2) = gl_y2(yzel2, now); Dk(3) = gl_z2(yzel2, now)
			! Dk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			r3 = sqrt(Dk(2)**2 + Dk(3)**2)
			
			if(k < N3) then
			    yzel2 = gl_RAY_B(par_n_HP, j, k + 1)
			else
				yzel2 = gl_RAY_B(par_n_HP, j, 1)
			end if
			Ek(1) = gl_x2(yzel2, now); Ek(2) = gl_y2(yzel2, now); Ek(3) = gl_z2(yzel2, now)
			! Ek = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			r4 = sqrt(Ek(2)**2 + Ek(3)**2)
			
			rr = (r1 + r2 + r3 + r4)/4.0
			
			!dist = sqrt( (Dk(1) - Ak(1))**2 + (Dk(2) - Ak(2))**2 + (Dk(3) - Ak(3))**2 )
			!dist = max(dist, 1.0_8)
			
			if (gl_Point_num(yzel) > 0) then
			    !vel = gl_Point_num(yzel) * par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time/(dist)
				vel = par_nat_HP * 0.03 * gl_Point_num(yzel) * ((Ak/r) * (rr - r)) * ddt  ! 0.001
				vel(1) = 0.0
			else
				!vel = par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time/(dist)
				vel = par_nat_HP * 0.03 * ((Ak/r) * (rr - r)) * ddt
				vel(1) = 0.0
			end if
			
			
			gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
			gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
			gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
	
	end subroutine Cuda_Move_all_2
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	attributes(global) subroutine Cuda_Move_all_3(now)   ! ������������� ��������� �� ����� K
	use MY_CUDA, gl_TS => dev_gl_TS, gl_Gran_neighbour => dev_gl_Gran_neighbour, gl_Gran_normal2 => dev_gl_Gran_normal2, &
		gl_Cell_par => dev_gl_Cell_par, gl_all_Gran => dev_gl_all_Gran, gl_Point_num => dev_gl_Point_num, gl_x2 => dev_gl_x2, &
		gl_y2 => dev_gl_y2, gl_z2 => dev_gl_z2, gl_Gran_center2 => dev_gl_Gran_center2, gl_Vx => dev_gl_Vx, &
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS, &
		gl_RAY_A => dev_gl_RAY_A, gl_RAY_B => dev_gl_RAY_B, gl_RAY_C => dev_gl_RAY_C, gl_RAY_O => dev_gl_RAY_O, &
		gl_RAY_K => dev_gl_RAY_K, gl_RAY_D => dev_gl_RAY_D, gl_RAY_E => dev_gl_RAY_E, norm2 => dev_norm2, par_n_TS => dev_par_n_TS, &
		par_n_HP => dev_par_n_HP, par_n_BS => dev_par_n_BS
	use GEO_PARAM
	use cudafor
	
	implicit none
	integer, intent(in) :: now
	
	real(8) :: Time,  r, rr, r1, r2, r3, r4, ddt
    real(8) :: vel(3), Ak(3), Bk(3), Ck(3), Dk(3), Ek(3)
    integer :: i, j, k
	
	integer(4):: N2, N3
	integer(4) :: yzel, yzel2
	
	
	N3 = size(gl_RAY_K(1, 1, :))
    N2 = size(gl_RAY_K(1, :, 1))
	  
	Time = time_step
	
	k = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! ����� ������
	j = blockDim%y * (blockIdx%y - 1) + threadIdx%y   ! ����� ������
	
	if(k > N3 .or. j > N2) return
	
	if (j == 1) then
            return
	end if
	
	ddt = Time/0.000127
			
	yzel = gl_RAY_K(par_n_TS, j, k)
	Ak(1) = gl_x2(yzel, now); Ak(2) = gl_y2(yzel, now); Ak(3) = gl_z2(yzel, now)
	
	
			
	if (j < N2) then
		yzel2 = gl_RAY_K(par_n_TS, j + 1, k)
	else
		yzel2 = gl_RAY_B(par_n_TS, size(gl_RAY_B(par_n_TS, :, k)), k)
	end if
	Bk(1) = gl_x2(yzel2, now); Bk(2) = gl_y2(yzel2, now); Bk(3) = gl_z2(yzel2, now)
	
			
	if (k > 1) then
		yzel2 = gl_RAY_K(par_n_TS, j, k - 1)
	else
		yzel2 = gl_RAY_K(par_n_TS, j, N3)
	end if
	Ck(1) = gl_x2(yzel2, now); Ck(2) = gl_y2(yzel2, now); Ck(3) = gl_z2(yzel2, now)
			
	if (j > 1) then
		yzel2 = gl_RAY_K(par_n_TS, j - 1, k)
	else
		yzel2 = gl_RAY_K(par_n_TS, j, k + ceiling(1.0 * N3/2))
	end if
	Dk(1) = gl_x2(yzel2, now); Dk(2) = gl_y2(yzel2, now); Dk(3) = gl_z2(yzel2, now)
			
	if(k < N3) then
		yzel2 = gl_RAY_K(par_n_TS, j, k + 1)
	else
		yzel2 = gl_RAY_K(par_n_TS, j, 1)
	end if
	Ek(1) = gl_x2(yzel2, now); Ek(2) = gl_y2(yzel2, now); Ek(3) = gl_z2(yzel2, now)
			
	if (gl_Point_num(yzel) > 0) then
		vel = gl_Point_num(yzel) * par_nat_TS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time * ddt
	else
		vel = par_nat_TS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time * ddt
	end if
		
			
	gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
	gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
	gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
		
	
	end subroutine Cuda_Move_all_3
	
	
	
	
	
	
	attributes(global) subroutine Cuda_Move_all_4(now)   ! ������������� ��������� 
	use MY_CUDA, gl_TS => dev_gl_TS, gl_Gran_neighbour => dev_gl_Gran_neighbour, gl_Gran_normal2 => dev_gl_Gran_normal2, &
		gl_Cell_par => dev_gl_Cell_par, gl_all_Gran => dev_gl_all_Gran, gl_Point_num => dev_gl_Point_num, gl_x2 => dev_gl_x2, &
		gl_y2 => dev_gl_y2, gl_z2 => dev_gl_z2, gl_Gran_center2 => dev_gl_Gran_center2, gl_Vx => dev_gl_Vx, &
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS, &
		gl_RAY_A => dev_gl_RAY_A, gl_RAY_B => dev_gl_RAY_B, gl_RAY_C => dev_gl_RAY_C, gl_RAY_O => dev_gl_RAY_O, &
		gl_RAY_K => dev_gl_RAY_K, gl_RAY_D => dev_gl_RAY_D, gl_RAY_E => dev_gl_RAY_E, norm2 => dev_norm2, par_n_TS => dev_par_n_TS, &
		par_n_HP => dev_par_n_HP, par_n_BS => dev_par_n_BS
	use GEO_PARAM
	use cudafor
	
	implicit none
	integer, intent(in) :: now
	
	real(8) :: Time, dist,  r, rr, r1, r2, r3, r4, ddt
    real(8) :: vel(3), Ak(3), Bk(3), Ck(3), Dk(3), Ek(3)
    integer :: i, j, k
	
	integer(4):: N2, N3
	integer(4) :: yzel, yzel2
	
	
	N3 = size(gl_RAY_O(1, 1, :))
    N2 = size(gl_RAY_O(1, :, 1))
	  
	Time = time_step
	
	k = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! ����� ������
	j = blockDim%y * (blockIdx%y - 1) + threadIdx%y   ! ����� ������
	
	if(k > N3 .or. j > N2) return
	
	if (k /= 1 .and. j == 1) then
            return
	end if
	
	ddt = Time/0.000127
			
	! �������
	yzel = gl_RAY_O(1, j, k)
	Ak(1) = gl_x2(yzel, now); Ak(2) = gl_y2(yzel, now); Ak(3) = gl_z2(yzel, now)
	! Ak = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/)
	r = sqrt(Ak(2)**2 + Ak(3)**2)
			
	if (j < N2) then
		yzel2 = gl_RAY_O(1, j + 1, k)
	else
		yzel2 = yzel
	end if
	Bk(1) = gl_x2(yzel2, now); Bk(2) = gl_y2(yzel2, now); Bk(3) = gl_z2(yzel2, now)
	! Bk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
	r1 = sqrt(Bk(2)**2 + Bk(3)**2)
			
	if (k > 1) then
		yzel2 = gl_RAY_O(1, j, k - 1)
	else
		yzel2 = gl_RAY_O(1, j, N3)
	end if
	Ck(1) = gl_x2(yzel2, now); Ck(2) = gl_y2(yzel2, now); Ck(3) = gl_z2(yzel2, now)
	! Ck = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
	r2 = sqrt(Ck(2)**2 + Ck(3)**2)
			
	if (j > 1) then
		yzel2 = gl_RAY_O(1, j - 1, k)
	else
		return
		!yzel2 = gl_RAY_C(1, size(gl_RAY_C(1, :, k)), k)
	end if
	Dk(1) = gl_x2(yzel2, now); Dk(2) = gl_y2(yzel2, now); Dk(3) = gl_z2(yzel2, now)
	! Dk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
	r3 = sqrt(Dk(2)**2 + Dk(3)**2)
			
	if(k < N3) then
		yzel2 = gl_RAY_O(1, j, k + 1)
	else
		yzel2 = gl_RAY_O(1, j, 1)
	end if
	Ek(1) = gl_x2(yzel2, now); Ek(2) = gl_y2(yzel2, now); Ek(3) = gl_z2(yzel2, now)
	! Ek = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
	r4 = sqrt(Ek(2)**2 + Ek(3)**2)
	
	!dist = sqrt( (Dk(1) - Ak(1))**2 + (Dk(2) - Ak(2))**2 + (Dk(3) - Ak(3))**2 )
	!dist = max(dist, 1.0_8)
	
	rr = (r1 + r2 + r3 + r4)/4.0
			
	if (gl_Point_num(yzel) > 0) then
		!vel = gl_Point_num(yzel) * par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time/dist
		vel = par_nat_HP * 0.0003 * gl_Point_num(yzel) * ((Ak/r) * (rr - r)) * ddt
		vel(1) = 0.0
	else
		!vel = par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time/dist
		vel = par_nat_HP * 0.0003 * ((Ak/r) * (rr - r)) * ddt
		vel(1) = 0.0
	end if
			
			
	gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
	gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
	gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
		
	
	end subroutine Cuda_Move_all_4
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	attributes(global) subroutine Cuda_Move_all_A(now)   
	! �������� ����� �� � - �����
	use MY_CUDA, gl_TS => dev_gl_TS, gl_Gran_neighbour => dev_gl_Gran_neighbour, gl_Gran_normal2 => dev_gl_Gran_normal2, &
		gl_Cell_par => dev_gl_Cell_par, gl_all_Gran => dev_gl_all_Gran, gl_Point_num => dev_gl_Point_num, gl_x2 => dev_gl_x2, &
		gl_y2 => dev_gl_y2, gl_z2 => dev_gl_z2, gl_Gran_center2 => dev_gl_Gran_center2, gl_Vx => dev_gl_Vx, &
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS, &
		gl_RAY_A => dev_gl_RAY_A, gl_RAY_B => dev_gl_RAY_B, gl_RAY_C => dev_gl_RAY_C, gl_RAY_O => dev_gl_RAY_O, &
		gl_RAY_K => dev_gl_RAY_K, gl_RAY_D => dev_gl_RAY_D, gl_RAY_E => dev_gl_RAY_E, norm2 => dev_norm2, par_n_TS => dev_par_n_TS, &
		par_n_HP => dev_par_n_HP, par_n_BS => dev_par_n_BS, par_n_END => dev_par_n_END, &
		par_n_IA => dev_par_n_IA, par_n_IB => dev_par_n_IB, par_R_inner => dev_par_R_inner, &
		par_kk1 => dev_par_kk1, par_kk2 => dev_par_kk2, par_R_END => dev_par_R_END, par_kk12 => dev_par_kk12
	use GEO_PARAM
	use cudafor
	
	implicit none
    integer, intent(in) :: now
	
	real(8), device :: Time
    real(8), device :: vel(3), ER(3), KORD(3), ER2(3)
	real(8) :: R_TS, proect, R_HP, R_BS
    integer:: i, j, k
    real(8), device :: the, phi, r
    integer :: now2             ! ��� ��������� �� ������ ������ �� ������ now
	
	integer(4):: N1, N2, N3
	integer(4) :: yzel
	
    now2 = mod(now, 2) + 1  
	Time = time_step
	
	N3 = size(gl_RAY_A(1, 1, :))
    N2 = size(gl_RAY_A(1, :, 1))
	N1 = size(gl_RAY_A(:, 1, 1))
	  
	
	k = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! ����� ������
	j = blockDim%y * (blockIdx%y - 1) + threadIdx%y   ! ����� ������
	
	if(k > N3 .or. j > N2) return
	
	if (k /= 1 .and. j == 1) then
                    return
            end if
            
        ! ��������� ���������� �������� ���� � ������������
        the = (j - 1) * par_pi_8/2.0/(N2 - 1)
        phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
            
        ! TS
        yzel = gl_RAY_A(par_n_TS, j, k)
        if(gl_Point_num(yzel) == 0) then
            vel = 0.0
        else
            ! vel = (/gl_Vx(yzel), gl_Vy(yzel), gl_Vz(yzel)/)
			vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
            vel = vel/gl_Point_num(yzel)                       ! ����� �������� �������� ����� ����
        end if
            
        ! ������� ��� ���������� �������������
        gl_Point_num(yzel) = 0
        gl_Vx(yzel) = 0.0
        gl_Vy(yzel) = 0.0
        gl_Vz(yzel) = 0.0
            
		ER(1) = cos(the); ER(2) = sin(the) * cos(phi); ER(3) = sin(the) * sin(phi)
		proect = DOT_PRODUCT(vel * Time, ER)
        ! proect = DOT_PRODUCT(vel * Time, (/cos(the), sin(the) * cos(phi), sin(the) * sin(phi)/))  !  ������� �������� ����������� �� ������ ������ ����
        KORD(1) = gl_x2(yzel, now); KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now) 
		ER2 = proect * ER
		ER2  = KORD + ER2
		R_TS = norm2(ER2)  ! ����� ���������� �� TS
		
		!if(j == 1 .and. k==1)  write(*,*) "=    ", gl_x2(yzel, now)
            
        ! HP
        yzel = gl_RAY_A(par_n_HP, j, k)
        if(gl_Point_num(yzel) == 0) then
            vel = 0.0
        else
            ! vel = (/gl_Vx(yzel), gl_Vy(yzel), gl_Vz(yzel)/)
			vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
            vel = vel/gl_Point_num(yzel)                       ! ����� �������� �������� ����� ����
        end if
            
        ! ������� ��� ���������� �������������
        gl_Point_num(yzel) = 0
        gl_Vx(yzel) = 0.0
        gl_Vy(yzel) = 0.0
        gl_Vz(yzel) = 0.0
            
        proect = DOT_PRODUCT(vel * Time, ER)
		KORD(1) = gl_x2(yzel, now); KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now) 
		ER2 = proect * ER
		ER2  = KORD + ER2
        R_HP = norm2(ER2)  ! ����� ���������� �� HP
            
        ! BS
        yzel = gl_RAY_A(par_n_BS, j, k)
        if(gl_Point_num(yzel) == 0) then
            vel = 0.0
        else
            ! vel = (/gl_Vx(yzel), gl_Vy(yzel), gl_Vz(yzel)/)
			vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
            vel = vel/gl_Point_num(yzel)                       ! ����� �������� �������� ����� ����
		end if
	
            
        ! ������� ��� ���������� �������������
        gl_Point_num(yzel) = 0
        gl_Vx(yzel) = 0.0
        gl_Vy(yzel) = 0.0
        gl_Vz(yzel) = 0.0
            
        proect = DOT_PRODUCT(vel * Time, ER)
		KORD(1) = gl_x2(yzel, now); KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now) 
		ER2 = proect * ER
		ER2  = KORD + ER2
        R_BS = norm2(ER2)  ! ����� ���������� �� BS
		
            
        ! ����� ������� ���� ���������� ��������� �����, ����� ��, ��� � ��� ���������� �����
        do i = 1, N1
                
            if (i == 1) then
                CYCLE
            end if

            yzel = gl_RAY_A(i, j, k)
            ! ��������� ���������� ����� �� ����

            ! �� TS
			if (i <= par_n_IB) then  ! NEW
                    r =  par_R0 + (par_R_inner - par_R0) * (DBLE(i)/(par_n_IB))**par_kk1
			else if (i <= par_n_TS) then  ! �� ���������� = R_TS
                    r =  par_R_inner + (R_TS - par_R_inner) * (DBLE(i - par_n_IB)/(par_n_TS - par_n_IB))**par_kk12
            !if (i <= par_n_TS) then  ! �� ���������� = R_TS
            !    r =  par_R0 + (R_TS - par_R0) * (DBLE(i)/par_n_TS)**par_kk1
            else if (i <= par_n_HP) then  
                r = R_TS + (i - par_n_TS) * (R_HP - R_TS)/(par_n_HP - par_n_TS)
            else if (i <= par_n_BS) then 
                r = R_HP + (i - par_n_HP) * (R_BS - R_HP)/(par_n_BS - par_n_HP)
            else
                r = R_BS + (par_R_END - R_BS) * (DBLE(i- par_n_BS)/(par_n_END - par_n_BS))**(par_kk2 * (0.55 + 0.45 * cos(the)) )
            end if

            ! ���������� ����� ����������
            gl_x2(yzel, now2) = r * cos(the)
            gl_y2(yzel, now2) = r * sin(the) * cos(phi)
            gl_z2(yzel, now2) = r * sin(the) * sin(phi)
			
			!if(j == 1 .and. k==1 .and. i==20)  write(*,*) "=    ", gl_x2(yzel, now2)
                
        end do
	
	
	end subroutine Cuda_Move_all_A
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	attributes(global) subroutine Cuda_Move_all_B(now)   
	! �������� ����� �� B - �����
	use MY_CUDA, gl_TS => dev_gl_TS, gl_Gran_neighbour => dev_gl_Gran_neighbour, gl_Gran_normal2 => dev_gl_Gran_normal2, &
		gl_Cell_par => dev_gl_Cell_par, gl_all_Gran => dev_gl_all_Gran, gl_Point_num => dev_gl_Point_num, gl_x2 => dev_gl_x2, &
		gl_y2 => dev_gl_y2, gl_z2 => dev_gl_z2, gl_Gran_center2 => dev_gl_Gran_center2, gl_Vx => dev_gl_Vx, &
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS, &
		gl_RAY_A => dev_gl_RAY_A, gl_RAY_B => dev_gl_RAY_B, gl_RAY_C => dev_gl_RAY_C, gl_RAY_O => dev_gl_RAY_O, &
		gl_RAY_K => dev_gl_RAY_K, gl_RAY_D => dev_gl_RAY_D, gl_RAY_E => dev_gl_RAY_E, norm2 => dev_norm2, par_n_TS => dev_par_n_TS, &
		par_n_HP => dev_par_n_HP, par_n_BS => dev_par_n_BS, par_n_END => dev_par_n_END, par_triple_point => dev_par_triple_point, &
		par_n_IA => dev_par_n_IA, par_n_IB => dev_par_n_IB, par_R_inner => dev_par_R_inner, &
		par_kk1 => dev_par_kk1, par_kk12 => dev_par_kk12
	use GEO_PARAM
	use cudafor
	
	implicit none
    integer, intent(in) :: now
	
	real(8), device :: Time
    real(8), device :: vel(3), ER(3), KORD(3), ER2(3)
	real(8) :: R_TS, proect, R_HP, R_BS
    integer:: i, j, k
    real(8), device :: the, phi, r
    integer :: now2             ! ��� ��������� �� ������ ������ �� ������ now
	
	integer(4):: N1, N2, N3
	integer(4) :: yzel
	
    now2 = mod(now, 2) + 1  
	Time = time_step
	
	N3 = size(gl_RAY_B(1, 1, :))
    N2 = size(gl_RAY_B(1, :, 1))
	N1 = size(gl_RAY_B(:, 1, 1))
	  
	
	k = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! ����� ������
	j = blockDim%y * (blockIdx%y - 1) + threadIdx%y   ! ����� ������
	
	if(k > N3 .or. j > N2) return
	
	! ��������� ���������� �������� ���� � ������������
            the = par_pi_8/2.0 + (j) * par_triple_point/(N2)
            phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
            
            ! TS
            yzel = gl_RAY_B(par_n_TS, j, k)
            if(gl_Point_num(yzel) == 0) then
                vel = 0.0
			else
				vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
                vel = vel/gl_Point_num(yzel)                       ! ����� �������� �������� ����� ����
            end if
            
            ! ������� ��� ���������� �������������
            gl_Point_num(yzel) = 0
            gl_Vx(yzel) = 0.0
            gl_Vy(yzel) = 0.0
            gl_Vz(yzel) = 0.0
            
			ER(1) = cos(the); ER(2) = sin(the) * cos(phi); ER(3) = sin(the) * sin(phi)
            proect = DOT_PRODUCT(vel * Time, ER)  !  ������� �������� ����������� �� ������ ������ ����
			KORD(1) = gl_x2(yzel, now); KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now) 
			ER2 = KORD + proect * ER
            R_TS = norm2(ER2)  ! ����� ���������� �� TS
            
            ! HP
            yzel = gl_RAY_B(par_n_HP, j, k)
            if(gl_Point_num(yzel) == 0) then
                vel = 0.0
			else
				vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
                vel = vel/gl_Point_num(yzel)                       ! ����� �������� �������� ����� ����
			end if

            
            ! ������� ��� ���������� �������������
            gl_Point_num(yzel) = 0
            gl_Vx(yzel) = 0.0
            gl_Vy(yzel) = 0.0
            gl_Vz(yzel) = 0.0
            
            proect = DOT_PRODUCT(vel * Time, ER)
			KORD(1) = gl_x2(yzel, now); KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now) 
			ER2 = KORD + proect * ER
            R_HP = norm2(ER2)  ! ����� ���������� �� HP
			
                
            do i = 1, N1

                if (i == 1) CYCLE
                
                yzel = gl_RAY_B(i, j, k)
                ! ��������� ���������� ����� �� ����

                ! �� TS
				if (i <= par_n_IB) then  ! NEW
                    r =  par_R0 + (par_R_inner - par_R0) * (DBLE(i)/(par_n_IB))**par_kk1
                else if (i <= par_n_TS) then  ! �� ���������� = R_TS
                    r =  par_R_inner + (R_TS - par_R_inner) * (DBLE(i - par_n_IB)/(par_n_TS - par_n_IB))**par_kk12
                !if (i <= par_n_TS) then  ! �� ���������� = R_TS
                !    r =  par_R0 + (R_TS - par_R0) * (REAL(i, KIND = 4)/par_n_TS)**par_kk1
                else if (i <= par_n_HP) then  ! �� ���������� = par_R_character * 1.3
                    r = R_TS + (i - par_n_TS) * (R_HP - R_TS) /(par_n_HP - par_n_TS)
				end if
				

                ! ���������� ����� ����������
                gl_x2(yzel, now2) = r * cos(the)
                gl_y2(yzel, now2) = r * sin(the) * cos(phi)
                gl_z2(yzel, now2) = r * sin(the) * sin(phi)
                
               
            end do
	
	
	end subroutine Cuda_Move_all_B
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	attributes(global) subroutine Cuda_Move_all_C(now)   
	! �������� ����� �� C - �����
	use MY_CUDA, gl_TS => dev_gl_TS, gl_Gran_neighbour => dev_gl_Gran_neighbour, gl_Gran_normal2 => dev_gl_Gran_normal2, &
		gl_Cell_par => dev_gl_Cell_par, gl_all_Gran => dev_gl_all_Gran, gl_Point_num => dev_gl_Point_num, gl_x2 => dev_gl_x2, &
		gl_y2 => dev_gl_y2, gl_z2 => dev_gl_z2, gl_Gran_center2 => dev_gl_Gran_center2, gl_Vx => dev_gl_Vx, &
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS, &
		gl_RAY_A => dev_gl_RAY_A, gl_RAY_B => dev_gl_RAY_B, gl_RAY_C => dev_gl_RAY_C, gl_RAY_O => dev_gl_RAY_O, &
		gl_RAY_K => dev_gl_RAY_K, gl_RAY_D => dev_gl_RAY_D, gl_RAY_E => dev_gl_RAY_E, norm2 => dev_norm2, par_n_TS => dev_par_n_TS, &
		par_n_HP => dev_par_n_HP, par_n_BS => dev_par_n_BS, par_n_END => dev_par_n_END, par_kk2 => dev_par_kk2, &
		par_R_END => dev_par_R_END
	use GEO_PARAM
	use cudafor
	
	implicit none
    integer, intent(in) :: now
	
	real(8), device :: Time
    real(8), device :: vel(3), ER(3), KORD(3), ER2(3)
	real(8) :: R_TS, proect, R_HP, R_BS
    integer:: i, j, k
    real(8), device :: the, phi, r,  x, y, z, rr
    integer :: now2             ! ��� ��������� �� ������ ������ �� ������ now
	
	integer(4):: N1, N2, N3
	integer(4) :: yzel
	
    now2 = mod(now, 2) + 1  
	Time = time_step
	
	N3 = size(gl_RAY_C(1, 1, :))
    N2 = size(gl_RAY_C(1, :, 1))
	N1 = size(gl_RAY_C(:, 1, 1))
	  
	
	k = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! ����� ������
	j = blockDim%y * (blockIdx%y - 1) + threadIdx%y   ! ����� ������
	
	if(k > N3 .or. j > N2) return
	
	
	! ��������� ���������� �������� ���� � ������������
                phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
                
                ! ��������� ���������� ����� �� ����
                x = gl_x2(gl_RAY_B(par_n_HP, j, k), now2)
                y = gl_y2(gl_RAY_B(par_n_HP, j, k), now2)
                z = gl_z2(gl_RAY_B(par_n_HP, j, k), now2)
                rr = (y**2 + z**2)**(0.5)
                
                ! BS     ����� ����� ��������� BS �� � ��������� �� ������� ���� A
                yzel = gl_RAY_A(par_n_BS, size(gl_RAY_A(1, :, 1)), k)
				ER(1) = 0.0_8; ER(2) = gl_y2(yzel, now2); ER(3) = gl_z2(yzel, now2)
                R_BS = norm2(ER)  ! ����� ���������� �� BS
				
            
            do i = 1, N1
                
                if(i == 1) CYCLE
                
                yzel = gl_RAY_C(i, j, k)

                if (i <= par_n_BS - par_n_HP + 1) then
                    r = rr + (i - 1) * (R_BS - rr)/(par_n_BS - par_n_HP)
                else
                    r = R_BS + (DBLE(i - (par_n_BS - par_n_HP + 1))/(N1 - (par_n_BS - par_n_HP + 1) ))**par_kk2 * (par_R_END - R_BS)
                end if
                
                gl_x2(yzel, now2) = x
                gl_y2(yzel, now2) = r * cos(phi)
                gl_z2(yzel, now2) = r * sin(phi)
                

            end do
	
	
	
	end subroutine Cuda_Move_all_C
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	attributes(global) subroutine Cuda_Move_all_O(now)   
	! �������� ����� �� O - �����
	use MY_CUDA, gl_TS => dev_gl_TS, gl_Gran_neighbour => dev_gl_Gran_neighbour, gl_Gran_normal2 => dev_gl_Gran_normal2, &
		gl_Cell_par => dev_gl_Cell_par, gl_all_Gran => dev_gl_all_Gran, gl_Point_num => dev_gl_Point_num, gl_x2 => dev_gl_x2, &
		gl_y2 => dev_gl_y2, gl_z2 => dev_gl_z2, gl_Gran_center2 => dev_gl_Gran_center2, gl_Vx => dev_gl_Vx, &
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS, &
		gl_RAY_A => dev_gl_RAY_A, gl_RAY_B => dev_gl_RAY_B, gl_RAY_C => dev_gl_RAY_C, gl_RAY_O => dev_gl_RAY_O, &
		gl_RAY_K => dev_gl_RAY_K, gl_RAY_D => dev_gl_RAY_D, gl_RAY_E => dev_gl_RAY_E, norm2 => dev_norm2, par_n_TS => dev_par_n_TS, &
		par_n_HP => dev_par_n_HP, par_n_BS => dev_par_n_BS, par_n_END => dev_par_n_END, par_m_BC => dev_par_m_BC, &
		par_kk2 => dev_par_kk2, par_kk3 => dev_par_kk3, par_R_END => dev_par_R_END, par_R_LEFT => dev_par_R_LEFT, &
		par_kk31 => dev_par_kk31
	use GEO_PARAM
	use cudafor
	
	implicit none
    integer, intent(in) :: now
	
	real(8), device :: Time
    real(8), device :: vel(3), ER(3), KORD(3), ER2(3)
	real(8) :: R_TS, proect, R_HP, R_BS
    integer:: i, j, k
    real(8), device :: the, phi, r, x, xx
    integer :: now2             ! ��� ��������� �� ������ ������ �� ������ now
	
	integer(4):: N1, N2, N3
	integer(4) :: yzel
	
    now2 = mod(now, 2) + 1  
	Time = time_step
	
	N3 = size(gl_RAY_O(1, 1, :))
    N2 = size(gl_RAY_O(1, :, 1))
	N1 = size(gl_RAY_O(:, 1, 1))
	  
	
	k = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! ����� ������
	j = blockDim%y * (blockIdx%y - 1) + threadIdx%y   ! ����� ������
	
	if(k > N3 .or. j > N2) return
	
	phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
            yzel = gl_RAY_O(1, j, k)
            if(gl_Point_num(yzel) == 0) then
                vel = 0.0
			else
				vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
                vel = vel/gl_Point_num(yzel)                       ! ����� �������� �������� ����� ����
            end if
            
            ! ������� ��� ���������� �������������
            gl_Point_num(yzel) = 0
            gl_Vx(yzel) = 0.0_8
            gl_Vy(yzel) = 0.0_8
            gl_Vz(yzel) = 0.0_8
            
			ER(1) = 0.0_8; ER(2) = cos(phi); ER(3) = sin(phi)
            proect = DOT_PRODUCT(vel * Time, ER)  !  ������� �������� ����������� �� ������ ������ ����
			KORD(1) = 0.0_8; KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now)
			ER2 = proect * ER
			ER2  = KORD + ER2
            R_HP = norm2(ER2)  ! ����� ���������� �� HP
			
			! ��������� ����������� �������� � ���
			if(R_HP < 20.0_8) then
				R_HP = 20.0_8
			end if
            
            xx = gl_x2(gl_RAY_B(par_n_HP, par_m_BC, k), now2)              ! ������������� �� x - ���������� ������� ����� B �� ���������� � ���� ��������� (k)
            x = xx - (DBLE(j)/N2)**par_kk31 * (xx - par_R_LEFT)
            
            ! BS     ����� ����� ��������� BS �� � ��������� �� ������� ���� A
            yzel = gl_RAY_A(par_n_BS, size(gl_RAY_A(1, :, 1)), k)
			KORD(1) = 0.0_8; KORD(2) = gl_y2(yzel, now2); KORD(3) = gl_z2(yzel, now2)
            R_BS = norm2(KORD)  ! ����� ���������� �� BS
			
			if (R_BS < R_HP + 10.0) then
				R_BS = R_HP + 10.0
			end if
            
            do i = 1, N1
                yzel = gl_RAY_O(i, j, k)
                
                if (i <= par_n_BS - par_n_HP + 1) then
                    r = R_HP + (i - 1) * (R_BS - R_HP)/(par_n_BS - par_n_HP)
                else
                    r = R_BS + (DBLE(i - (par_n_BS - par_n_HP + 1))/(N1 - (par_n_BS - par_n_HP + 1) ))**par_kk2 * (par_R_END - R_BS)
                end if


                gl_x2(yzel, now2) = x
                gl_y2(yzel, now2) = r * cos(phi)
                gl_z2(yzel, now2) = r * sin(phi)
                
            end do
	
	
	end subroutine Cuda_Move_all_O
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	attributes(global) subroutine Cuda_Move_all_K(now)   
	! �������� ����� �� B - �����
	use MY_CUDA, gl_TS => dev_gl_TS, gl_Gran_neighbour => dev_gl_Gran_neighbour, gl_Gran_normal2 => dev_gl_Gran_normal2, &
		gl_Cell_par => dev_gl_Cell_par, gl_all_Gran => dev_gl_all_Gran, gl_Point_num => dev_gl_Point_num, gl_x2 => dev_gl_x2, &
		gl_y2 => dev_gl_y2, gl_z2 => dev_gl_z2, gl_Gran_center2 => dev_gl_Gran_center2, gl_Vx => dev_gl_Vx, &
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS, &
		gl_RAY_A => dev_gl_RAY_A, gl_RAY_B => dev_gl_RAY_B, gl_RAY_C => dev_gl_RAY_C, gl_RAY_O => dev_gl_RAY_O, &
		gl_RAY_K => dev_gl_RAY_K, gl_RAY_D => dev_gl_RAY_D, gl_RAY_E => dev_gl_RAY_E, norm2 => dev_norm2, par_n_TS => dev_par_n_TS, &
		par_n_HP => dev_par_n_HP, par_n_BS => dev_par_n_BS, par_n_END => dev_par_n_END, par_triple_point => dev_par_triple_point, &
		par_n_IA => dev_par_n_IA, par_n_IB => dev_par_n_IB, par_R_inner => dev_par_R_inner, &
		par_kk1 => dev_par_kk1, par_kk12 => dev_par_kk12
	use GEO_PARAM
	use cudafor
	
	implicit none
    integer, intent(in) :: now
	
	real(8), device :: Time
    real(8), device :: vel(3), ER(3), KORD(3), ER2(3)
	real(8) :: R_TS, proect, R_HP, R_BS
    integer:: i, j, k
    real(8), device :: the, phi, r
    integer :: now2             ! ��� ��������� �� ������ ������ �� ������ now
	
	integer(4):: N1, N2, N3
	integer(4) :: yzel
	
    now2 = mod(now, 2) + 1  
	Time = time_step
	
	N3 = size(gl_RAY_K(1, 1, :))
    N2 = size(gl_RAY_K(1, :, 1))
	N1 = size(gl_RAY_K(:, 1, 1))
	  
	
	k = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! ����� ������
	j = blockDim%y * (blockIdx%y - 1) + threadIdx%y   ! ����� ������
	
	if(k > N3 .or. j > N2) return
	if (k /= 1 .and. j == 1) return
	
	! ��������� ���������� �������� ���� � ������������
    the = par_pi_8/2.0 + par_triple_point + (N2 - j + 1) * (par_pi_8/2.0 - par_triple_point)/(N2)
    phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
            
            
    yzel = gl_RAY_K(par_n_TS, j, k)
    if(gl_Point_num(yzel) == 0) then
        vel = 0.0
	else
		vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
        vel = vel/gl_Point_num(yzel)                       ! ����� �������� �������� ����� ����
    end if
            
    ! ������� ��� ���������� �������������
    gl_Point_num(yzel) = 0
    gl_Vx(yzel) = 0.0
    gl_Vy(yzel) = 0.0
    gl_Vz(yzel) = 0.0
            
	ER(1) = cos(the); ER(2) = sin(the) * cos(phi); ER(3) = sin(the) * sin(phi)
    proect = DOT_PRODUCT(vel * Time, ER)  !  ������� �������� ����������� �� ������ ������ ����
    KORD(1) = gl_x2(yzel, now); KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now)
	ER2 = proect * ER
	ER2  = KORD + ER2
	R_TS = norm2(ER2)  ! ����� ���������� �� TS
            
    do i = 1, N1

        if (i == 1) CYCLE
                
        ! ��������� ���������� ����� �� ����
        yzel = gl_RAY_K(i, j, k)
		
		if (i <= par_n_IB) then  ! NEW
            r =  par_R0 + (par_R_inner - par_R0) * (DBLE(i)/(par_n_IB))**par_kk1
        else 
            r =  par_R_inner + (R_TS - par_R_inner) * (DBLE(i - par_n_IB)/(par_n_TS - par_n_IB))**par_kk12
		end if
				
        !r =  par_R0 + (R_TS - par_R0) * (REAL(i, KIND = 4)/par_n_TS)**par_kk1


        ! ���������� ����� ����������
        gl_x2(yzel, now2) = r * cos(the)
        gl_y2(yzel, now2) = r * sin(the) * cos(phi)
        gl_z2(yzel, now2) = r * sin(the) * sin(phi)
                
                
    end do
	
	
	
	end subroutine Cuda_Move_all_K
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	attributes(global) subroutine Cuda_Move_all_D(now)   
	! �������� ����� �� D - �����
	use MY_CUDA, gl_TS => dev_gl_TS, gl_Gran_neighbour => dev_gl_Gran_neighbour, gl_Gran_normal2 => dev_gl_Gran_normal2, &
		gl_Cell_par => dev_gl_Cell_par, gl_all_Gran => dev_gl_all_Gran, gl_Point_num => dev_gl_Point_num, gl_x2 => dev_gl_x2, &
		gl_y2 => dev_gl_y2, gl_z2 => dev_gl_z2, gl_Gran_center2 => dev_gl_Gran_center2, gl_Vx => dev_gl_Vx, &
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS, &
		gl_RAY_A => dev_gl_RAY_A, gl_RAY_B => dev_gl_RAY_B, gl_RAY_C => dev_gl_RAY_C, gl_RAY_O => dev_gl_RAY_O, &
		gl_RAY_K => dev_gl_RAY_K, gl_RAY_D => dev_gl_RAY_D, gl_RAY_E => dev_gl_RAY_E, norm2 => dev_norm2, par_n_TS => dev_par_n_TS, &
		par_n_HP => dev_par_n_HP, par_n_BS => dev_par_n_BS, par_n_END => dev_par_n_END, par_m_BC => dev_par_m_BC, &
		par_kk3 => dev_par_kk3, par_R_LEFT => dev_par_R_LEFT
	use GEO_PARAM
	use cudafor
	
	implicit none
    integer, intent(in) :: now
	
	real(8), device :: Time
    real(8), device :: vel(3), ER(3), KORD(3), ER2(3)
	real(8) :: R_TS, proect, R_HP, R_BS
    integer:: i, j, k
    real(8), device :: the, phi, r, x, y, z, xx, rr, rrr
    integer :: now2             ! ��� ��������� �� ������ ������ �� ������ now
	
	integer(4):: N1, N2, N3
	integer(4) :: yzel
	
    now2 = mod(now, 2) + 1  
	Time = time_step
	
	N3 = size(gl_RAY_D(1, 1, :))
    N2 = size(gl_RAY_D(1, :, 1))
	N1 = size(gl_RAY_D(:, 1, 1))
	  
	
	k = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! ����� ������
	j = blockDim%y * (blockIdx%y - 1) + threadIdx%y   ! ����� ������
	
	if(k > N3 .or. j > N2) return
	
	if (k /= 1 .and. j == 1) return
	
	! ��������� ���������� �������� ���� � ������������
                phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
				
				! ��������� ���������� ����� �� ����

				if (j < N2) then
                    xx = gl_x2(gl_RAY_K(par_n_TS, j, k), now2)
                    y = gl_y2(gl_RAY_K(par_n_TS, j, k), now2)
                    z = gl_z2(gl_RAY_K(par_n_TS, j, k), now2)
                else
                    xx = gl_x2(gl_RAY_B(par_n_TS, par_m_BC, k), now2)
                    y = gl_y2(gl_RAY_B(par_n_TS, par_m_BC, k), now2)
                    z = gl_z2(gl_RAY_B(par_n_TS, par_m_BC, k), now2)
				end if
				
				!y = gl_y2(gl_RAY_B(par_n_TS, par_m_BC, k), now2)
				!z = gl_z2(gl_RAY_B(par_n_TS, par_m_BC, k), now2)
				
				r = sqrt(y**2 + z**2)
                
            do i = 1, N1

                if (i == 1) CYCLE

                ! ��������� ���������� ����� �� ����

    !            y = gl_y2(gl_RAY_O(1, i - 1, k), now2)
				!z = gl_z2(gl_RAY_O(1, i - 1, k), now2)
				!
				!rr = sqrt(y**2 + z**2) * 0.4

                yzel = gl_RAY_D(i, j, k)
                
                
                x = xx + (DBLE(i - 1)/(N1 - 1))**par_kk3 * (par_R_LEFT - xx)
				
				!rrr = min(r, rr) * (DBLE(j - 1)/(N2 - 1))
				
				! rrr = r * (1.0 - (DBLE(i - 1)**2/(N1 - 1)**2)) + (rr * (DBLE(j - 1)/(N2 - 1))) * ( (DBLE(i - 1)**2) /(N1 - 1)**2)


                ! ���������� ����� ����������
                gl_x2(yzel, now2) = x
                gl_y2(yzel, now2) = r * cos(phi)
                gl_z2(yzel, now2) = r * sin(phi)
                
            end do
	
	
	end subroutine Cuda_Move_all_D
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	attributes(global) subroutine Cuda_Move_all_E(now)   
	! �������� ����� �� E - �����
	use MY_CUDA, gl_TS => dev_gl_TS, gl_Gran_neighbour => dev_gl_Gran_neighbour, gl_Gran_normal2 => dev_gl_Gran_normal2, &
		gl_Cell_par => dev_gl_Cell_par, gl_all_Gran => dev_gl_all_Gran, gl_Point_num => dev_gl_Point_num, gl_x2 => dev_gl_x2, &
		gl_y2 => dev_gl_y2, gl_z2 => dev_gl_z2, gl_Gran_center2 => dev_gl_Gran_center2, gl_Vx => dev_gl_Vx, &
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS, &
		gl_RAY_A => dev_gl_RAY_A, gl_RAY_B => dev_gl_RAY_B, gl_RAY_C => dev_gl_RAY_C, gl_RAY_O => dev_gl_RAY_O, &
		gl_RAY_K => dev_gl_RAY_K, gl_RAY_D => dev_gl_RAY_D, gl_RAY_E => dev_gl_RAY_E, norm2 => dev_norm2, par_n_TS => dev_par_n_TS, &
		par_n_HP => dev_par_n_HP, par_n_BS => dev_par_n_BS, par_n_END => dev_par_n_END
	use GEO_PARAM
	use cudafor
	
	implicit none
    integer, intent(in) :: now
	
	real(8), device :: Time
    real(8), device :: vel(3), ER(3), KORD(3), ER2(3)
	real(8) :: R_TS, proect, R_HP, R_BS
    integer:: i, j, k
    real(8), device :: the, phi, r, x, y, z, x2, y2, z2
    integer :: now2             ! ��� ��������� �� ������ ������ �� ������ now
	
	integer(4):: N1, N2, N3
	integer(4) :: yzel
	
    now2 = mod(now, 2) + 1  
	Time = time_step
	
	N3 = size(gl_RAY_E(1, 1, :))
    N2 = size(gl_RAY_E(1, :, 1))
	N1 = size(gl_RAY_E(:, 1, 1))
	  
	
	k = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! ����� ������
	j = blockDim%y * (blockIdx%y - 1) + threadIdx%y   ! ����� ������
	
	if(k > N3 .or. j > N2) return
	
	do i = 1, N1

        if (i == 1) CYCLE

        if (i == N1) CYCLE

        x = gl_x2(gl_RAY_E(1, j, k), now2)
        y = gl_y2(gl_RAY_E(1, j, k), now2)
        z = gl_z2(gl_RAY_E(1, j, k), now2)

        x2 = gl_x2(gl_RAY_O(1, j, k), now2)
        y2 = gl_y2(gl_RAY_O(1, j, k), now2)
        z2 = gl_z2(gl_RAY_O(1, j, k), now2)

        yzel = gl_RAY_E(i, j, k)

        gl_x2(yzel, now2) = x + (x2 - x) * (i - 1)/(N1 - 1)
        gl_y2(yzel, now2) = y + (y2 - y) * (i - 1)/(N1 - 1)
        gl_z2(yzel, now2) = z + (z2 - z) * (i - 1)/(N1 - 1)
                
    end do
	
	
	end subroutine Cuda_Move_all_E
	
	
	
	
	
	
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	attributes(global) subroutine Cuda_calc_all_Gran_move_1(now)   
	! �������� ����� �� E - �����
	use MY_CUDA, gl_TS => dev_gl_TS, gl_Gran_neighbour => dev_gl_Gran_neighbour, gl_Gran_normal2 => dev_gl_Gran_normal2, &
		gl_Cell_par => dev_gl_Cell_par, gl_all_Gran => dev_gl_all_Gran, gl_Point_num => dev_gl_Point_num, gl_x2 => dev_gl_x2, &
		gl_y2 => dev_gl_y2, gl_z2 => dev_gl_z2, gl_Gran_center2 => dev_gl_Gran_center2, gl_Vx => dev_gl_Vx, &
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS, &
		gl_RAY_A => dev_gl_RAY_A, gl_RAY_B => dev_gl_RAY_B, gl_RAY_C => dev_gl_RAY_C, gl_RAY_O => dev_gl_RAY_O, &
		gl_RAY_K => dev_gl_RAY_K, gl_RAY_D => dev_gl_RAY_D, gl_RAY_E => dev_gl_RAY_E, norm2 => dev_norm2, par_n_TS => dev_par_n_TS, &
		par_n_HP => dev_par_n_HP, par_n_BS => dev_par_n_BS, par_n_END => dev_par_n_END, gl_Gran_square2 => dev_gl_Gran_square2
	use GEO_PARAM
	use cudafor
	
	implicit none
    integer, intent(in) :: now
	
	integer :: Ngran, iter
    real(8) :: p(3, 4), Vol, D
    real(8) :: a(3), b(3), c(3), S
    real(8) :: dist, di, gr_center(3)
    integer :: i, j, k, ll, grc, ii

    !print*, "calc_all_Gran_move"
    ! ���� �� ������ - ������� ������� �����, � �������
    Ngran = size(gl_all_Gran(1,:))
	
	iter = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! ����� ������
	
	if(iter > Ngran) return
	
	grc = 0
    p = 0.0
    a = 0.0
    b = 0.0
    c = 0.0
    gr_center = 0.0


    do j = 1, 4
        i = gl_all_Gran(j, iter)
		p(1,j) = gl_x2(i, now); p(2,j) = gl_y2(i, now); p(3,j) = gl_z2(i, now)
    end do
    a = p(:,3) - p(:,1)
    b = p(:,4) - p(:,2)

    ! ����� ��������� ����� �����
    gr_center = p(:, 1)
    grc = 1
    if (gl_all_Gran(2, iter) /= gl_all_Gran(1, iter)) then
        grc = grc + 1
        gr_center = gr_center + p(:,2)
    end if

    if (gl_all_Gran(3, iter) /= gl_all_Gran(2, iter) .and. gl_all_Gran(3, iter) /= gl_all_Gran(1, iter)) then
        grc = grc + 1
        gr_center = gr_center + p(:,3)
    end if

    if (gl_all_Gran(4, iter) /= gl_all_Gran(1, iter) .and. gl_all_Gran(4, iter) /= gl_all_Gran(2, iter) &
        .and. gl_all_Gran(4, iter) /= gl_all_Gran(3, iter)) then
        grc = grc + 1
        gr_center = gr_center + p(:,4)
    end if

    gr_center = gr_center/grc

    gl_Gran_center2(:, iter, now) = gr_center

    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)

    S = norm2(c)  ! S = S/2
    c = c/S
    S = S/2.0


    ! ����� �������� ������� ����� � ������� � ����� ������!
    gl_Gran_normal2(:, iter, now) = c
    gl_Gran_square2(iter, now) = S
	
	end subroutine Cuda_calc_all_Gran_move_1
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	attributes(global) subroutine Cuda_calc_all_Gran_move_2(now)   
	use MY_CUDA, gl_TS => dev_gl_TS, gl_Gran_neighbour => dev_gl_Gran_neighbour, gl_Gran_normal2 => dev_gl_Gran_normal2, &
		gl_Cell_par => dev_gl_Cell_par, gl_all_Gran => dev_gl_all_Gran, gl_Point_num => dev_gl_Point_num, gl_x2 => dev_gl_x2, &
		gl_y2 => dev_gl_y2, gl_z2 => dev_gl_z2, gl_Gran_center2 => dev_gl_Gran_center2, gl_Vx => dev_gl_Vx, &
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS, &
		gl_RAY_A => dev_gl_RAY_A, gl_RAY_B => dev_gl_RAY_B, gl_RAY_C => dev_gl_RAY_C, gl_RAY_O => dev_gl_RAY_O, &
		gl_RAY_K => dev_gl_RAY_K, gl_RAY_D => dev_gl_RAY_D, gl_RAY_E => dev_gl_RAY_E, norm2 => dev_norm2, par_n_TS => dev_par_n_TS, &
		par_n_HP => dev_par_n_HP, par_n_BS => dev_par_n_BS, par_n_END => dev_par_n_END, gl_Gran_square2 => dev_gl_Gran_square2, &
		gl_Cell_Volume2 => dev_gl_Cell_Volume2, gl_Cell_dist => dev_gl_Cell_dist, gl_all_Cell => dev_gl_all_Cell, gl_Cell_center2 => dev_gl_Cell_center2, &
		gl_Cell_gran => dev_gl_Cell_gran
	use GEO_PARAM
	use cudafor
	
	implicit none
    integer, intent(in) :: now
	
	integer :: Ngran, iter
    real(8) :: p(3, 4), Vol, D
    real(8) :: a(3), b(3), c(3), S
    real(8) :: dist, di, gr_center(3)
    integer :: i, j, k, ll, grc, ii, now2

    !print*, "calc_all_Gran_move"
    now2 = mod(now, 2) + 1 
	
    Ngran = size(gl_all_Cell(1,:))
	
	iter = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! ����� ������
	
	if(iter > Ngran) return
	
	Vol = 0.0
        ll = 0
        dist = 10000.0 * par_R_character
        c = 0.0
		
		do j = 1, 8
                i = gl_all_Cell(j, iter)
				c(1) = c(1) + gl_x2(i, now)
				c(2) = c(2) + gl_y2(i, now)
				c(3) = c(3) + gl_z2(i, now)
        end do
        c = c/8.0

        gl_Cell_center2(:, iter, now) = c


        ! ��������� ����� ������ ������ ������� ����� �������� �� ������ �����
        do j = 1, 6
            i = gl_Cell_gran(j, iter)   ! ���� �� ������� ��� ����� ������

            if (i == 0) CYCLE
            k = gl_all_Gran(1, i)  ! ����� ������� ���� �����
			b(1) = gl_x2(k, now)
			b(2) = gl_y2(k, now)
			b(3) = gl_z2(k, now)
            a = gl_Gran_normal2(:,i, now)
            di = dabs(DOT_PRODUCT(c,a) - DOT_PRODUCT(a,b))  ! ���������� �� ����� �� ���������
            Vol = Vol + gl_Gran_square2(i, now) * di/3.0
            dist = MIN(dist, di)
            ll = ll + 1
        end do



        gl_Cell_Volume2(iter, now) = Vol
        gl_Cell_dist(iter) = dist
	
	
	end subroutine Cuda_calc_all_Gran_move_2
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	attributes(global) subroutine Cuda_Calc_move_TS(now)
	use MY_CUDA, gl_TS => dev_gl_TS, gl_Gran_neighbour => dev_gl_Gran_neighbour, gl_Gran_normal2 => dev_gl_Gran_normal2, &
		gl_Cell_par => dev_gl_Cell_par, gl_all_Gran => dev_gl_all_Gran, gl_Point_num => dev_gl_Point_num, gl_x2 => dev_gl_x2, &
		gl_y2 => dev_gl_y2, gl_z2 => dev_gl_z2, gl_Gran_center2 => dev_gl_Gran_center2, norm2 => dev_norm2, gl_Vx => dev_gl_Vx, &
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS
	use GEO_PARAM
	use cudafor
	use My_func
	
	implicit none
	integer, intent(in) :: now
	
	integer i, Num, gr, s1, s2, j, yzel
    real(8) :: normal(3), qqq1(8), qqq2(8), POTOK(8), ray(3)
    real(8) :: a1, a2, v1, v2, norm, b1, b2, c1, c2
	real(8) :: www, dsl, dsp, dsc, the1, the2, center(3), center2(3)
	integer(4) :: metod
	
	i = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! ����� ������
	Num = size(gl_TS)
	
	
	
	if (i > Num) return
	
	metod = 3 
	www = 0.0
	
		
    gr = gl_TS(i)
    s1 = gl_Gran_neighbour(1, gr)
    s2 = gl_Gran_neighbour(2, gr)
    normal = gl_Gran_normal2(:, gr, now)
    qqq1 = gl_Cell_par(1:8, s1)
    qqq2 = gl_Cell_par(1:8, s2)
		
        
    call chlld(metod, normal(1), normal(2), normal(3), &
            www, qqq1, qqq2, dsl, dsp, dsc, POTOK)
	
		
    dsl = dsl * koef1
	
	center = gl_Gran_center2(:, gr, now)
		
	the1 = polar_angle(center(1), sqrt(center(2)**2 + center(3)**2))
		
	do j = 1, 4
        yzel = gl_all_Gran(j, gr)
			
		center2(1) = gl_x2(yzel, now)
		center2(2) = gl_y2(yzel, now)
		center2(3) = gl_z2(yzel, now)
		the2 = polar_angle(center2(1), sqrt(center2(2)**2 + center2(3)**2))
			
		if (the2 < the1 .and. the2 > 0.18) CYCLE
			
			
        gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsl
        gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsl
        gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsl
        gl_Point_num(yzel) = gl_Point_num(yzel) + 1
	end do
        
    return  ! ����������� � ���� ������, ��������� � ���������
	
	
	
    a1 = sqrt(ggg * qqq1(5)/qqq1(1))  ! �������� �����
    a2 = sqrt(ggg * qqq2(5)/qqq2(1))
	  
        
    if (norm2(qqq1(2:4)) <= a1 .and. norm2(qqq2(2:4)) <= a2) then ! ���������� �������� �� ��� ����
        do j = 1, 4
            yzel = gl_all_Gran(j, gr)
			
			! ���� ����������� �������
			select case (mod(yzel, 4))
			case(0)
				do while ( atomiccas(dev_mutex_1, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case(1)
				do while ( atomiccas(dev_mutex_2, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case(2)
				do while ( atomiccas(dev_mutex_3, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case(3)
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case default
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			end select
			
				gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsl
				gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsl
				gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsl
				gl_Point_num(yzel) = gl_Point_num(yzel) + 1
				
			call threadfence()                               ! ���������� ������ � ������
			select case (mod(yzel, 4))
			case(0)
				dev_mutex_1 = 0 
			case(1)
				dev_mutex_2 = 0 
			case(2)
				dev_mutex_3 = 0 
			case(3)
				dev_mutex_4 = 0 
			case default
				dev_mutex_4 = 0 
			end select
				
        end do
        return  ! ����������� � ���� ������, ��������� � ���������
    end if
        
    do j = 1, 4
            yzel = gl_all_Gran(j, gr)
                
            !ray = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/) - gl_Gran_center2(:, gr, now)
			ray(1) = gl_x2(yzel, now) - gl_Gran_center2(1, gr, now)
			ray(2) = gl_y2(yzel, now) - gl_Gran_center2(2, gr, now)
			ray(3) = gl_z2(yzel, now) - gl_Gran_center2(3, gr, now)
				
            norm = norm2(ray)
            ray = ray/norm
            v1 = DOT_PRODUCT(qqq1(2:4), ray)
            v2 = DOT_PRODUCT(qqq2(2:4), ray)
			
                
            if (v1 >= 0 .or. v2 >= 0) then
				
				select case (mod(yzel, 4))
			case(0)
				do while ( atomiccas(dev_mutex_1, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case(1)
				do while ( atomiccas(dev_mutex_2, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case(2)
				do while ( atomiccas(dev_mutex_3, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case(3)
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case default
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			end select
			
                gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsl
                gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsl
                gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsl
                gl_Point_num(yzel) = gl_Point_num(yzel) + 1
				call threadfence()                               ! ���������� ������ � ������
				select case (mod(yzel, 4))
				case(0)
					dev_mutex_1 = 0 
				case(1)
					dev_mutex_2 = 0 
				case(2)
					dev_mutex_3 = 0 
				case(3)
					dev_mutex_4 = 0 
				case default
					dev_mutex_4 = 0 
				end select
				
				
                CYCLE   ! ��������� � ���������� ����
            end if
                
            b1 = DOT_PRODUCT(qqq1(6:8), ray)
            b2 = DOT_PRODUCT(qqq2(6:8), ray)
                
            c1 = dabs(b1)/dsqrt(4.0 * par_pi_8 * qqq1(1))
            c2 = dabs(b2)/dsqrt(4.0 * par_pi_8 * qqq2(1))
                
            b1 = norm2(qqq1(6:8))
            b2 = norm2(qqq2(6:8))
			
                
            ! ���� ������ �������� ������ ��� �������� �� �������� ����� � ������������� �������� ����� �� ����� �����������
            if( dabs(v1) <= max(a1, 0.5 * (sqrt(b1 * b1/(4.0 * par_pi_8 * qqq1(1)) + a1 * a1+ 2.0 * a1 * c1) &
                + sqrt(b1 * b1/(4.0 * par_pi_8 * qqq1(1)) + a1 * a1 - 2.0 * a1 * c1)  ) ) &
                .and. dabs(v2) <= max(a2, 0.5 * (sqrt(b2 * b2/(4.0 * par_pi_8 * qqq2(1)) + a2 * a2+ 2.0 * a2 * c2) &
                + sqrt(b2 * b2/(4.0 * par_pi_8 * qqq2(1)) + a2 * a2 - 2.0 * a2 * c2)  ) ) ) then
				
				select case (mod(yzel, 4))
				case(0)
					do while ( atomiccas(dev_mutex_1, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
					end do                                        ! ���������� ������ � ������
				case(1)
					do while ( atomiccas(dev_mutex_2, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
					end do                                        ! ���������� ������ � ������
				case(2)
					do while ( atomiccas(dev_mutex_3, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
					end do                                        ! ���������� ������ � ������
				case(3)
					do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
					end do                                        ! ���������� ������ � ������
				case default
					do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
					end do                                        ! ���������� ������ � ������
				end select
			
                gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsl
                gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsl
                gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsl
                gl_Point_num(yzel) = gl_Point_num(yzel) + 1
				call threadfence()                               ! ���������� ������ � ������
				select case (mod(yzel, 4))
				case(0)
					dev_mutex_1 = 0 
				case(1)
					dev_mutex_2 = 0 
				case(2)
					dev_mutex_3 = 0 
				case(3)
					dev_mutex_4 = 0 
				case default
					dev_mutex_4 = 0 
				end select
				
                CYCLE
            end if
                
    end do
	
	end subroutine Cuda_Calc_move_TS
	
	
	attributes(global) subroutine Cuda_Calc_move_HP(now)
	use MY_CUDA, gl_TS => dev_gl_TS, gl_Gran_neighbour => dev_gl_Gran_neighbour, gl_Gran_normal2 => dev_gl_Gran_normal2, &
		gl_Cell_par => dev_gl_Cell_par, gl_all_Gran => dev_gl_all_Gran, gl_Point_num => dev_gl_Point_num, gl_x2 => dev_gl_x2, &
		gl_y2 => dev_gl_y2, gl_z2 => dev_gl_z2, gl_Gran_center2 => dev_gl_Gran_center2, norm2 => dev_norm2, gl_Vx => dev_gl_Vx, &
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS
	use GEO_PARAM
	use cudafor
	
	implicit none
	integer, intent(in) :: now
	
	integer i, Num, gr, s1, s2, j, yzel
    real(8) :: normal(3), qqq1(8), qqq2(8), POTOK(8), ray(3), vec(3), center(3), center2(3)
    real(8) :: a1, a2, v1, v2, norm, b1, b2, c1, c2
	real(8) :: www, dsl, dsp, dsc
	integer(4) :: metod
	
	i = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! ����� ������
	Num = size(gl_Contact)
	
	
	if (i > Num) return
	
	metod = 3
	www = 0.0
	

    gr = gl_Contact(i)
    s1 = gl_Gran_neighbour(1, gr)
    s2 = gl_Gran_neighbour(2, gr)
    normal = gl_Gran_normal2(:, gr, now)
    qqq1 = gl_Cell_par(1:8, s1)
    qqq2 = gl_Cell_par(1:8, s2)
	
	
	! ������ ���������� ���������� ���������� ����
	if (gl_Gran_center2(1, gr, now) > -100.0) then
		qqq1(6:8) = qqq1(6:8) - DOT_PRODUCT(normal, qqq1(6:8)) * normal
		qqq2(6:8) = qqq2(6:8) - DOT_PRODUCT(normal, qqq2(6:8)) * normal
		gl_Cell_par(6:8, s1) = qqq1(6:8)
		gl_Cell_par(6:8, s2) = qqq2(6:8)
	end if
	
	
	
	
	vec = 0.5 * (qqq1(2:4) + qqq2(2:4))
	vec = vec - DOT_PRODUCT(vec, normal) * normal
	center = gl_Gran_center2(:, gr, now)
        
    call chlld(metod, normal(1), normal(2), normal(3), &
            www, qqq1, qqq2, dsl, dsp, dsc, POTOK)
        
    dsc = (dsc + 0.0 * DOT_PRODUCT(0.5 * (qqq1(2:4) + qqq2(2:4)), normal)) * koef2
	
	do j = 1, 4
            yzel = gl_all_Gran(j, gr)
			center2(1) = gl_x2(yzel, now)
			center2(2) = gl_y2(yzel, now)
			center2(3) = gl_z2(yzel, now)
			
			if ( sqrt(center2(2)**2 + center2(3)**2) > 15.0) then
			    if( DOT_PRODUCT(vec, center2 - center) < 0.0) CYCLE
			end if
			
			! ���� ����������� �������
			select case (mod(yzel, 4))
			case(0)
				do while ( atomiccas(dev_mutex_1, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case(1)
				do while ( atomiccas(dev_mutex_2, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case(2)
				do while ( atomiccas(dev_mutex_3, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case(3)
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case default
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			end select
			
            gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsc
            gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsc
            gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsc
            gl_Point_num(yzel) = gl_Point_num(yzel) + 1
			
			call threadfence()                               ! ���������� ������ � ������
			select case (mod(yzel, 4))
			case(0)
				dev_mutex_1 = 0 
			case(1)
				dev_mutex_2 = 0 
			case(2)
				dev_mutex_3 = 0 
			case(3)
				dev_mutex_4 = 0 
			case default
				dev_mutex_4 = 0 
			end select
				
        end do
	return  ! ����������� � ���� ������, ��������� � ���������
	
        
    a1 = sqrt(ggg * qqq1(5)/qqq1(1))  ! �������� �����
    a2 = sqrt(ggg * qqq2(5)/qqq2(1))
	
        
    if (norm2(qqq1(2:4)) <= a1 .and. norm2(qqq2(2:4)) <= a2) then ! ���������� �������� �� ��� ����
        do j = 1, 4
            yzel = gl_all_Gran(j, gr)
			
			! ���� ����������� �������
			select case (mod(yzel, 4))
			case(0)
				do while ( atomiccas(dev_mutex_1, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case(1)
				do while ( atomiccas(dev_mutex_2, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case(2)
				do while ( atomiccas(dev_mutex_3, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case(3)
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case default
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			end select
			
            gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsc
            gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsc
            gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsc
            gl_Point_num(yzel) = gl_Point_num(yzel) + 1
			
			call threadfence()                               ! ���������� ������ � ������
				select case (mod(yzel, 4))
				case(0)
					dev_mutex_1 = 0 
				case(1)
					dev_mutex_2 = 0 
				case(2)
					dev_mutex_3 = 0 
				case(3)
					dev_mutex_4 = 0 
				case default
					dev_mutex_4 = 0 
				end select
				
        end do
        return
    end if
        
    do j = 1, 4
            yzel = gl_all_Gran(j, gr)
                
            !ray = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/) - gl_Gran_center2(:, gr, now)
			ray(1) = gl_x2(yzel, now) - gl_Gran_center2(1, gr, now)
			ray(2) = gl_y2(yzel, now) - gl_Gran_center2(2, gr, now)
			ray(3) = gl_z2(yzel, now) - gl_Gran_center2(3, gr, now)
				
            norm = norm2(ray)
            ray = ray/norm
            v1 = DOT_PRODUCT(qqq1(2:4), ray)
            v2 = DOT_PRODUCT(qqq2(2:4), ray)
                
            if (v1 >= 0 .or. v2 >= 0) then
				
				! ���� ����������� �������
			select case (mod(yzel, 4))
			case(0)
				do while ( atomiccas(dev_mutex_1, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case(1)
				do while ( atomiccas(dev_mutex_2, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case(2)
				do while ( atomiccas(dev_mutex_3, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case(3)
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case default
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			end select
			
                gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsc
                gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsc
                gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsc
                gl_Point_num(yzel) = gl_Point_num(yzel) + 1
				
				call threadfence()                               ! ���������� ������ � ������
				select case (mod(yzel, 4))
				case(0)
					dev_mutex_1 = 0 
				case(1)
					dev_mutex_2 = 0 
				case(2)
					dev_mutex_3 = 0 
				case(3)
					dev_mutex_4 = 0 
				case default
					dev_mutex_4 = 0 
				end select
				
                CYCLE   ! ��������� � ���������� ����
            end if
                
            b1 = DOT_PRODUCT(qqq1(6:8), ray)
            b2 = DOT_PRODUCT(qqq2(6:8), ray)
                
            c1 = dabs(b1)/dsqrt(4.0 * par_pi_8 * qqq1(1))
            c2 = dabs(b2)/dsqrt(4.0 * par_pi_8 * qqq2(1))
                
            b1 = norm2(qqq1(6:8))
            b2 = norm2(qqq2(6:8))
                
            ! ���� ������ �������� ������ ��� �������� �� �������� ����� � ������������� �������� ����� �� ����� �����������
            if( dabs(v1) <= max(a1, 0.5 * (sqrt(b1 * b1/(4.0 * par_pi_8 * qqq1(1)) + a1 * a1+ 2.0 * a1 * c1) &
                + sqrt(b1 * b1/(4.0 * par_pi_8 * qqq1(1)) + a1 * a1 - 2.0 * a1 * c1)  ) ) &
                .and. dabs(v2) <= max(a2, 0.5 * (sqrt(b2 * b2/(4.0 * par_pi_8 * qqq2(1)) + a2 * a2+ 2.0 * a2 * c2) &
                + sqrt(b2 * b2/(4.0 * par_pi_8 * qqq2(1)) + a2 * a2 - 2.0 * a2 * c2)  ) ) ) then
				
				! ���� ����������� �������
			select case (mod(yzel, 4))
			case(0)
				do while ( atomiccas(dev_mutex_1, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case(1)
				do while ( atomiccas(dev_mutex_2, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case(2)
				do while ( atomiccas(dev_mutex_3, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case(3)
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case default
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			end select
			
                gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsc
                gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsc
                gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsc
                gl_Point_num(yzel) = gl_Point_num(yzel) + 1
				
				call threadfence()                               ! ���������� ������ � ������
				select case (mod(yzel, 4))
				case(0)
					dev_mutex_1 = 0 
				case(1)
					dev_mutex_2 = 0 
				case(2)
					dev_mutex_3 = 0 
				case(3)
					dev_mutex_4 = 0 
				case default
					dev_mutex_4 = 0 
				end select
				
                CYCLE
            end if
                
    end do
        
	
	
	end subroutine Cuda_Calc_move_HP
	
	
	attributes(global) subroutine Cuda_Calc_move_BS(now)
	use MY_CUDA, gl_TS => dev_gl_TS, gl_Gran_neighbour => dev_gl_Gran_neighbour, gl_Gran_normal2 => dev_gl_Gran_normal2, &
		gl_Cell_par => dev_gl_Cell_par, gl_all_Gran => dev_gl_all_Gran, gl_Point_num => dev_gl_Point_num, gl_x2 => dev_gl_x2, &
		gl_y2 => dev_gl_y2, gl_z2 => dev_gl_z2, gl_Gran_center2 => dev_gl_Gran_center2, norm2 => dev_norm2, gl_Vx => dev_gl_Vx, &
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS
	use GEO_PARAM
	use cudafor
	
	implicit none
	integer, intent(in) :: now
	
	integer i, Num, gr, s1, s2, j, yzel
    real(8) :: normal(3), qqq1(8), qqq2(8), POTOK(8), ray(3)
    real(8) :: a1, a2, v1, v2, norm, b1, b2, c1, c2
	real(8) :: www, dsl, dsp, dsc
	integer(4) :: metod
	
	i = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! ����� ������
	Num = size(gl_BS)
	
	
	if (i > Num) return
	
	metod = 2
	www = 0.0
	
    
    gr = gl_BS(i)
    s1 = gl_Gran_neighbour(1, gr)
    s2 = gl_Gran_neighbour(2, gr)
    normal = gl_Gran_normal2(:, gr, now)
    qqq1 = gl_Cell_par(1:8, s1)
    qqq2 = gl_Cell_par(1:8, s2)
        
    call chlld(metod, normal(1), normal(2), normal(3), &
            www, qqq1, qqq2, dsl, dsp, dsc, POTOK)
	
	
        
    dsp = dsp * koef3
	
	
        
    a1 = sqrt(ggg * qqq1(5)/qqq1(1))  ! �������� �����
    a2 = sqrt(ggg * qqq2(5)/qqq2(1))
        
    if (norm2(qqq1(2:4)) <= a1 .and. norm2(qqq2(2:4)) <= a2) then ! ���������� �������� �� ��� ����
        do j = 1, 4
            yzel = gl_all_Gran(j, gr)
			
			! ���� ����������� �������
			select case (mod(yzel, 4))
			case(0)
				do while ( atomiccas(dev_mutex_1, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case(1)
				do while ( atomiccas(dev_mutex_2, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case(2)
				do while ( atomiccas(dev_mutex_3, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case(3)
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case default
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			end select
			
            gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsp
            gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsp
            gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsp
            gl_Point_num(yzel) = gl_Point_num(yzel) + 1
			
			call threadfence()                               ! ���������� ������ � ������
				select case (mod(yzel, 4))
				case(0)
					dev_mutex_1 = 0 
				case(1)
					dev_mutex_2 = 0 
				case(2)
					dev_mutex_3 = 0 
				case(3)
					dev_mutex_4 = 0 
				case default
					dev_mutex_4 = 0 
				end select
				
        end do
        return
    end if
        
    do j = 1, 4
            yzel = gl_all_Gran(j, gr)
                
            !ray = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/) - gl_Gran_center2(:, gr, now)
			ray(1) = gl_x2(yzel, now) - gl_Gran_center2(1, gr, now)
			ray(2) = gl_y2(yzel, now) - gl_Gran_center2(2, gr, now)
			ray(3) = gl_z2(yzel, now) - gl_Gran_center2(3, gr, now)
				
            norm = norm2(ray)
            ray = ray/norm
            v1 = DOT_PRODUCT(qqq1(2:4), ray)
            v2 = DOT_PRODUCT(qqq2(2:4), ray)
                
            if (v1 >= 0 .or. v2 >= 0) then
				
				! ���� ����������� �������
			select case (mod(yzel, 4))
			case(0)
				do while ( atomiccas(dev_mutex_1, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case(1)
				do while ( atomiccas(dev_mutex_2, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case(2)
				do while ( atomiccas(dev_mutex_3, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case(3)
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case default
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			end select
			
                gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsp
                gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsp
                gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsp
                gl_Point_num(yzel) = gl_Point_num(yzel) + 1
				
				call threadfence()                               ! ���������� ������ � ������
				select case (mod(yzel, 4))
				case(0)
					dev_mutex_1 = 0 
				case(1)
					dev_mutex_2 = 0 
				case(2)
					dev_mutex_3 = 0 
				case(3)
					dev_mutex_4 = 0 
				case default
					dev_mutex_4 = 0 
				end select
				
                CYCLE   ! ��������� � ���������� ����
            end if
                
            b1 = DOT_PRODUCT(qqq1(6:8), ray)
            b2 = DOT_PRODUCT(qqq2(6:8), ray)
                
            c1 = dabs(b1)/dsqrt(4.0 * par_pi_8 * qqq1(1))
            c2 = dabs(b2)/dsqrt(4.0 * par_pi_8 * qqq2(1))
                
            b1 = norm2(qqq1(6:8))
            b2 = norm2(qqq2(6:8))
                
            ! ���� ������ �������� ������ ��� �������� �� �������� ����� � ������������� �������� ����� �� ����� �����������
            if( dabs(v1) <= max(a1, 0.5 * (sqrt(b1 * b1/(4.0 * par_pi_8 * qqq1(1)) + a1 * a1+ 2.0 * a1 * c1) &
                + sqrt(b1 * b1/(4.0 * par_pi_8 * qqq1(1)) + a1 * a1 - 2.0 * a1 * c1)  ) ) &
                .and. dabs(v2) <= max(a2, 0.5 * (sqrt(b2 * b2/(4.0 * par_pi_8 * qqq2(1)) + a2 * a2+ 2.0 * a2 * c2) &
                + sqrt(b2 * b2/(4.0 * par_pi_8 * qqq2(1)) + a2 * a2 - 2.0 * a2 * c2)  ) ) ) then
				
				! ���� ����������� �������
			select case (mod(yzel, 4))
			case(0)
				do while ( atomiccas(dev_mutex_1, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case(1)
				do while ( atomiccas(dev_mutex_2, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case(2)
				do while ( atomiccas(dev_mutex_3, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case(3)
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			case default
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! ���������� ������ � ������ dev_mutex_1
				end do                                        ! ���������� ������ � ������
			end select
			
                gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsp
                gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsp
                gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsp
                gl_Point_num(yzel) = gl_Point_num(yzel) + 1
				
				call threadfence()                               ! ���������� ������ � ������
				select case (mod(yzel, 4))
				case(0)
					dev_mutex_1 = 0 
				case(1)
					dev_mutex_2 = 0 
				case(2)
					dev_mutex_3 = 0 
				case(3)
					dev_mutex_4 = 0 
				case default
					dev_mutex_4 = 0 
				end select
				
                CYCLE
            end if
                
    end do
        
	
	end subroutine Cuda_Calc_move_BS
	
	
	attributes(global) subroutine CUF_MGD_grans_MF(now)
	! ��������� �����, ��� + �����������
	use MY_CUDA, gl_TS => dev_gl_TS, gl_Gran_neighbour => dev_gl_Gran_neighbour, gl_Gran_normal2 => dev_gl_Gran_normal2, &
		gl_Cell_par => dev_gl_Cell_par, gl_all_Gran => dev_gl_all_Gran, gl_Point_num => dev_gl_Point_num, gl_x2 => dev_gl_x2, &
		gl_y2 => dev_gl_y2, gl_z2 => dev_gl_z2, gl_Gran_center2 => dev_gl_Gran_center2, gl_Vx => dev_gl_Vx, &
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS, &
		gl_RAY_A => dev_gl_RAY_A, gl_RAY_B => dev_gl_RAY_B, gl_RAY_C => dev_gl_RAY_C, gl_RAY_O => dev_gl_RAY_O, &
		gl_RAY_K => dev_gl_RAY_K, gl_RAY_D => dev_gl_RAY_D, gl_RAY_E => dev_gl_RAY_E, norm2 => dev_norm2, par_n_TS => dev_par_n_TS, &
		par_n_HP => dev_par_n_HP, par_n_BS => dev_par_n_BS, par_n_END => dev_par_n_END, gl_Gran_square2 => dev_gl_Gran_square2, &
		gl_Cell_Volume2 => dev_gl_Cell_Volume2, gl_Cell_dist => dev_gl_Cell_dist, gl_all_Cell => dev_gl_all_Cell, gl_Cell_center2 => dev_gl_Cell_center2, &
		gl_Cell_gran => dev_gl_Cell_gran, gl_Cell_par_MF => dev_gl_Cell_par_MF, gl_Gran_type => dev_gl_Gran_type, &
		gl_Gran_POTOK => dev_gl_Gran_POTOK, gl_Gran_POTOK_MF => dev_gl_Gran_POTOK_MF, gl_Gran_info => dev_gl_Gran_info, &
		gl_Gran_scheme => dev_gl_Gran_scheme
	use GEO_PARAM
	implicit none
	integer, intent(in) :: now
	
	integer(4) :: gr  ! ���������� ����� ������� �����
	
    integer(4) :: st, s1, s2, i, j, k, zone, metod, now2
    real(8) :: qqq1(9), qqq2(9), qqq(9)  ! ���������� � ������
    real(8) :: fluid1(5, 4), fluid2(5, 4)
    real(8) :: dist, dsl, dsc, dsp, distant(3)
    real(8) :: POTOK(9), ttest(3)
    real(8) :: POTOK_MF(5)
    real(8) :: POTOK_MF_all(5, 4)
    real(8) :: time, Volume, TT, U8, rad1, rad2, aa, bb, cc, wc
    real(8) :: SOURSE(5,5)  ! ��������� �����, �������� � ������� ��� ������ � ������� ����� ������������
    real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3
	
	logical :: null_bn
	
	now2 = mod(now, 2) + 1
	
	time = 100000.0
	gr = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! ����� ������
	
	!if(gr == 1) then
	!	print*, "Hellow from potok (1, 1) ", dev_Ngran
	!end if
	
	TT = time_step
	
	
	! ----------------------------------------------------------- ����� ������ ����������� ��� �� ����� ----------------------
	if (gr > dev_Ngran) return
	
	metod = 2
	if(gl_Gran_info(gr) == 2) return
	
	POTOK = 0.0
	s1 = gl_Gran_neighbour(1, gr)
	s2 = gl_Gran_neighbour(2, gr)
	qqq1 = gl_Cell_par(:, s1)
	fluid1 = gl_Cell_par_MF(:, :, s1)   ! ��������� ��������� ��������� ��� ������������
	null_bn = .False.
	
	distant = gl_Gran_center2(:, gr, now) - gl_Cell_center2(:, s1, now)
	dist = norm2(distant)
	
	! ��������� ������ ��������� ��������������� ��������
            if(norm2(qqq1(2:4))/sqrt(ggg*qqq1(5)/qqq1(1)) > 2.2) then
                rad1 = norm2(gl_Cell_center2(:, s1, now))
                rad2 = norm2(gl_Gran_center2(:, gr, now))
                qqq1(1) = qqq1(1) * rad1**2 / rad2**2
                qqq1(9) = qqq1(9) * rad1**2 / rad2**2
                qqq1(5) = qqq1(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
                ! �������� ������ � ����������� �.�.
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
                fluid2 = gl_Cell_par_MF(:, :, s2)   ! ��������� ��������� ��������� ��� ������������
                !dist = min(gl_Cell_dist(s1), gl_Cell_dist(s2))   
				
				distant = gl_Gran_center2(:, gr, now) - gl_Cell_center2(:, s2, now)
				dist = min( dist, norm2(distant))

                ! ��������� ������ ��������� ��������������� ��������
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

            else  ! � ������ ��������� ����� - ��������� �������
                !if (norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int) CYCLE
                if(s2 == -1) then  ! ���������� �����
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
					
					qqq2 = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 1.0_8, -par_B_inf * cos(par_alphaB_inf), -par_B_inf * sin(par_alphaB_inf), 0.0_8, 100.0_8/)
                    fluid2(:, 1) = fluid1(:, 1)
                    fluid2(:, 2) = fluid1(:, 2)
                    fluid2(:, 3) = fluid1(:, 3)
                    fluid2(:, 4) = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 0.5_8/)
                else  ! ����� ����� ������ �������
                    !dist = gl_Cell_dist(s1)
                    qqq2 = qqq1
                    fluid2 = fluid1
                    !qqq2(5) = 1.0_8
                    if(qqq2(2) > 0.2 * par_Velosity_inf) then
                        qqq2(2) = 0.2 * par_Velosity_inf ! ����� ��������
                    end if

                end if
            end if
            
            ! ����� ��������� �������� �������� �����
            wc = DOT_PRODUCT((gl_Gran_center2(:, gr, now2) -  gl_Gran_center2(:, gr, now))/TT, gl_Gran_normal2(:, gr, now))
			
	
			
			metod = gl_Gran_scheme(gr)
			
			if(gl_Gran_type(gr) == 2 .or. gl_Gran_type(gr) == 1) metod = 3
			
			!if(gl_Gran_type(gr) == 2) null_bn = .True.
			
            if (.False.) then !(gl_Gran_type(gr) == 1) then
				call chlld_Q(metod, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
                wc, qqq1, qqq2, dsl, dsp, dsc, POTOK, null_bn, 0)
			else
				call chlld_Q(metod, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
                wc, qqq1, qqq2, dsl, dsp, dsc, POTOK, null_bn)
			end if
	
	
	
            time = min(time, 0.9 * dist/(max(dabs(dsl), dabs(dsp))+ dabs(wc)) )   ! REDUCTION
            gl_Gran_POTOK(1:9, gr) = POTOK * gl_Gran_square2(gr, now)
			
			gl_Gran_POTOK(10, gr) = 0.5 * DOT_PRODUCT(gl_Gran_normal2(:, gr, now), qqq1(6:8) + qqq2(6:8)) * gl_Gran_square2(gr, now)

			if( .False.) then
				!if (gr == 614302) then
				write(*,*) "__***__ 614302"
				write(*,*) gl_Gran_info(gr), s1, s2
				write(*,*) "__A__"
				write(*,*)  gl_Gran_POTOK(1, gr), gl_Gran_POTOK(2, gr), gl_Gran_POTOK(3, gr), gl_Gran_POTOK(4, gr), &
						gl_Gran_POTOK(5, gr), gl_Gran_POTOK(6, gr), gl_Gran_POTOK(7, gr), gl_Gran_POTOK(8, gr), gl_Gran_POTOK(9, gr), &
						gl_Gran_POTOK(10, gr)
				write(*,*) "__B__"
				write(*,*) gl_Gran_square2(gr, now)
				write(*,*) "__C__"
				write(*,*) dist, wc
				write(*,*) "__D__"
				write(*,*) gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now)
				write(*,*) "__E__"
				write(*,*) qqq1(1), qqq1(2), qqq1(3), qqq1(4), qqq1(5), qqq1(6), qqq1(7), qqq1(8), qqq1(9)
				write(*,*) "__F__"
				write(*,*) qqq2(1), qqq2(2), qqq2(3), qqq2(4), qqq2(5), qqq2(6), qqq2(7), qqq2(8), qqq2(9)
				
				
			end if
	
			
			
            call chlld_gd(0, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
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
	
	! ----------------------------------------------------------- �������� �� ����� ������� ----------------------
	
	!if (time_step2 > time) then 
		! �� ������ ��� �� ������� � ������ ��-�� ������ �����
		if (gl_Cell_center2(1, s1, now) > -150.0) then
			time =  atomicmin(time_step2, time)   ! ��������� �������� ������ ������������ ��������
		end if
	!end if
	
	
	end subroutine  CUF_MGD_grans_MF
	
	
	attributes(global) subroutine CUF_MGD_cells_MF(now)
	! ��������� �����, ��� + �����������
	use MY_CUDA, gl_TS => dev_gl_TS, gl_Gran_neighbour => dev_gl_Gran_neighbour, gl_Gran_normal2 => dev_gl_Gran_normal2, &
		gl_Cell_par => dev_gl_Cell_par, gl_all_Gran => dev_gl_all_Gran, gl_Point_num => dev_gl_Point_num, gl_x2 => dev_gl_x2, &
		gl_y2 => dev_gl_y2, gl_z2 => dev_gl_z2, gl_Gran_center2 => dev_gl_Gran_center2, gl_Vx => dev_gl_Vx, &
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS, &
		gl_RAY_A => dev_gl_RAY_A, gl_RAY_B => dev_gl_RAY_B, gl_RAY_C => dev_gl_RAY_C, gl_RAY_O => dev_gl_RAY_O, &
		gl_RAY_K => dev_gl_RAY_K, gl_RAY_D => dev_gl_RAY_D, gl_RAY_E => dev_gl_RAY_E, norm2 => dev_norm2, par_n_TS => dev_par_n_TS, &
		par_n_HP => dev_par_n_HP, par_n_BS => dev_par_n_BS, par_n_END => dev_par_n_END, gl_Gran_square2 => dev_gl_Gran_square2, &
		gl_Cell_Volume2 => dev_gl_Cell_Volume2, gl_Cell_dist => dev_gl_Cell_dist, gl_all_Cell => dev_gl_all_Cell, gl_Cell_center2 => dev_gl_Cell_center2, &
		gl_Cell_gran => dev_gl_Cell_gran, gl_Cell_par_MF => dev_gl_Cell_par_MF, gl_Gran_type => dev_gl_Gran_type, &
		gl_Gran_POTOK => dev_gl_Gran_POTOK, gl_Gran_POTOK_MF => dev_gl_Gran_POTOK_MF, gl_Gran_info => dev_gl_Gran_info, gl_Cell_info => dev_gl_Cell_info, &
		Calc_sourse_MF => dev_Calc_sourse_MF, par_kk1 => dev_par_kk1
	use GEO_PARAM
	implicit none
	integer, intent(in) :: now
	
	integer(4) :: gr
	
    integer(4) :: st, s1, s2, i, j, k, zone, now2, ijk
    real(8) :: qqq1(9), qqq2(9), qqq(9)  ! ���������� � ������
    real(8) :: fluid1(5, 4), fluid2(5, 4)
    real(8) :: dist, dsl, dsc, dsp
    real(8) :: POTOK(9), ttest(3)
    real(8) :: POTOK_MF(5)
    real(8) :: POTOK_MF_all(5, 4)
    real(8) :: time, Volume, U8, rad1, rad2, aa, bb, cc, Volume2, sks
    real(8) :: SOURSE(5,5)  ! ��������� �����, �������� � ������� ��� ������ � ������� ����� ������������
    real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3
	logical :: l_1
	
	now2 = mod(now, 2) + 1
	gr = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! ����� ������
	time = time_step
	
	if (gr > dev_Ncell) return
	
	! ----------------------------------------------------------- ����� ������ ����������� ��� �� ����� ----------------------
	
	if(gl_Cell_info(gr) == 0) return
            l_1 = .TRUE.
            if (norm2(gl_Cell_center2(:, gr, now)) <= par_R0 + (par_R_character - par_R0) * (3.0_8/par_n_TS)**par_kk1) l_1 = .FALSE.    ! �� ������� ������ �����
            POTOK = 0.0
			sks = 0.0
            SOURSE = 0.0
            POTOK_MF_all = 0.0
            Volume = gl_Cell_Volume2(gr, now)
            Volume2 = gl_Cell_Volume2(gr, now2)
			
			
	
			!Volume2 = Volume
	
            qqq = gl_Cell_par(:, gr)
            fluid1 = gl_Cell_par_MF(:, :, gr)
            ! ������������ ������ ����� �����
            do i = 1, 6
                j = gl_Cell_gran(i, gr)
                if (j == 0) CYCLE
				if (j < 0) write(*, *) "ERROR 3876tfghjuyghejk"
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


            ! ���������� ���� � ������� ���������
            if(qqq(9)/qqq(1) < 50.0) then
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

            call Calc_sourse_MF(qqq, fluid1, SOURSE, zone)  ! ��������� ���������


            if (l_1 == .TRUE.) then
                ro3 = qqq(1)* Volume / Volume2 - time * POTOK(1) / Volume2
                Q3 = qqq(9)* Volume / Volume2 - time * POTOK(9) / Volume2
                if (ro3 <= 0.0_8) then
                    write(*, *) "Ro < 0  1490 ", ro3, gl_Cell_center2(1, gr, now), gl_Cell_center2(2, gr, now), gl_Cell_center2(3, gr, now)
					write(*, *) qqq(1), Q3
					write(*, *) Volume , Volume2
					ro3 = 0.01
				end if
				
				
				
                u3 = (qqq(1) * qqq(2)* Volume / Volume2 - time * (POTOK(2) + (qqq(6)/cpi4) * sks) / Volume2 + time * SOURSE(2, 1)) / ro3
                v3 = (qqq(1) * qqq(3)* Volume / Volume2 - time * (POTOK(3) + (qqq(7)/cpi4) * sks) / Volume2 + time * SOURSE(3, 1)) / ro3
                w3 = (qqq(1) * qqq(4)* Volume / Volume2 - time * (POTOK(4) + (qqq(8)/cpi4) * sks) / Volume2 + time * SOURSE(4, 1)) / ro3
                
				bx3 = qqq(6) * Volume / Volume2 - time * (POTOK(6) + qqq(2) * sks) / Volume2
				by3 = qqq(7) * Volume / Volume2 - time * (POTOK(7) + qqq(3) * sks) / Volume2
				bz3 = qqq(8) * Volume / Volume2 - time * (POTOK(8) + qqq(4) * sks) / Volume2
				
                p3 = ((  ( qqq(5) / (ggg - 1.0) + 0.5 * qqq(1) * norm2(qqq(2:4))**2 + (qqq(6)**2 + qqq(7)**2 + qqq(8)**2) / 25.13274122871834590768 )* Volume / Volume2 &
                    - time * ( POTOK(5) + (DOT_PRODUCT(qqq(2:4), qqq(6:8))/cpi4) * sks)/ Volume2 + time * SOURSE(5, 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) - (bx3**2 + by3**2 + bz3**2) / 25.13274122871834590768 ) * (ggg - 1.0)
				
				!p3 = ((  ( qqq(5) / (ggg - 1.0) + 0.5 * qqq(1) * norm2(qqq(2:4))**2 )* Volume / Volume2 &
    !                - time * ( POTOK(5) + (DOT_PRODUCT(qqq(2:4), qqq(6:8))/cpi4) * sks)/ Volume2 + time * SOURSE(5, 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) ) * (ggg - 1.0)
				!
				!bx3 = qqq(6) * Volume / Volume2 - time * (POTOK(6) + qqq(2) * sks) / Volume2
				!by3 = qqq(7) * Volume / Volume2 - time * (POTOK(7) + qqq(3) * sks) / Volume2
				!bz3 = qqq(8) * Volume / Volume2 - time * (POTOK(8) + qqq(4) * sks) / Volume2
				
				!if(gr == 50) then
				!	write(*,*) ro3, Q3, u3, v3, w3, p3, bx3, by3, bz3
				!end if
				

                if (p3 <= 0.0_8) then
                    !print*, "p < 0  plasma 2028 ", p3 , gl_Cell_center(:, gr)
                    p3 = 0.000001
                    !pause
                end if

                gl_Cell_par(:, gr) = (/ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3/)

            end if

            ! ������ ��������� ������ ���������� ��� ��������� ���������

            do i = 1, 4
                if (i == 1 .and. l_1 == .FALSE.) CYCLE       ! ���������� ���������� ����� ��� ����� 1
                if (l_1 == .FALSE.) SOURSE(:, i + 1) = 0.0       ! ���������� ���������� ����� ��� ����� 1
                ro3 = fluid1(1, i)* Volume / Volume2 - time * POTOK_MF_all(1, i) / Volume2 + time * SOURSE(1, i + 1)
                if (ro3 <= 0.0_8) then
					ro3 = 0.01
                    write(*, *) "Ro < 0  in ", i,  ro3
                end if
                u3 = (fluid1(1, i) * fluid1(2, i)* Volume / Volume2 - time * POTOK_MF_all(2, i) / Volume2 + time * SOURSE(2, i + 1)) / ro3
                v3 = (fluid1(1, i) * fluid1(3, i)* Volume / Volume2 - time * POTOK_MF_all(3, i) / Volume2 + time * SOURSE(3, i + 1)) / ro3
                w3 = (fluid1(1, i) * fluid1(4, i)* Volume / Volume2 - time * POTOK_MF_all(4, i) / Volume2 + time * SOURSE(4, i + 1)) / ro3
                p3 = ((  ( fluid1(5, i) / (ggg - 1.0) + 0.5 * fluid1(1, i) * norm2(fluid1(2:4, i))**2 )* Volume / Volume2 &
                    - time * POTOK_MF_all(5, i)/ Volume2 + time * SOURSE(5, i + 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) ) * (ggg - 1.0)
                if (p3 <= 0.0_8) then
                    p3 = 0.000001
				end if
				
				!if(gr == 50) then
				!	write(*,*) ro3, Q3, u3, v3, w3, p3
				!end if
				
				if(.False.) then !
				!if (gr == 206135 .and. i == 1) then
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
					write(*,*) time
					write(*,*) "_____"
					
					do ijk = 1, 6
                j = gl_Cell_gran(ijk, gr)
				write(*,*) "*********************"
				write(*,*) "gran = ", j, "type = ", gl_Gran_type(j)
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
	!				write(*,*) "____________________"
	!			end if
	
	end subroutine CUF_MGD_cells_MF
	
	
	attributes(global) subroutine CUF_MGD_grans_MF_inner()
	! ��������� �����, ��� + �����������
	use MY_CUDA, gl_TS => dev_gl_TS, gl_Gran_neighbour => dev_gl_Gran_neighbour, gl_Gran_normal2 => dev_gl_Gran_normal2, &
		gl_Cell_par => dev_gl_Cell_par, gl_all_Gran => dev_gl_all_Gran, gl_Point_num => dev_gl_Point_num, gl_x2 => dev_gl_x2, &
		gl_y2 => dev_gl_y2, gl_z2 => dev_gl_z2, gl_Gran_center2 => dev_gl_Gran_center2, gl_Vx => dev_gl_Vx, &
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS, &
		gl_RAY_A => dev_gl_RAY_A, gl_RAY_B => dev_gl_RAY_B, gl_RAY_C => dev_gl_RAY_C, gl_RAY_O => dev_gl_RAY_O, &
		gl_RAY_K => dev_gl_RAY_K, gl_RAY_D => dev_gl_RAY_D, gl_RAY_E => dev_gl_RAY_E, norm2 => dev_norm2, par_n_TS => dev_par_n_TS, &
		par_n_HP => dev_par_n_HP, par_n_BS => dev_par_n_BS, par_n_END => dev_par_n_END, gl_Gran_square2 => dev_gl_Gran_square2, &
		gl_Cell_Volume2 => dev_gl_Cell_Volume2, gl_Cell_dist => dev_gl_Cell_dist, gl_all_Cell => dev_gl_all_Cell, gl_Cell_center2 => dev_gl_Cell_center2, &
		gl_Cell_gran => dev_gl_Cell_gran, gl_Cell_par_MF => dev_gl_Cell_par_MF, gl_Gran_type => dev_gl_Gran_type, &
		gl_Gran_POTOK => dev_gl_Gran_POTOK, gl_Gran_POTOK_MF => dev_gl_Gran_POTOK_MF, gl_Gran_info => dev_gl_Gran_info, &
		gl_all_Gran_inner => dev_gl_all_Gran_inner, gl_Gran_center => dev_gl_Gran_center, gl_Cell_center => dev_gl_Cell_center, &
		gl_Gran_normal => dev_gl_Gran_normal, gl_Gran_square => dev_gl_Gran_square, gl_Cell_type => dev_gl_Cell_type, &
		gl_Cell_number => dev_gl_Cell_number, gl_Cell_Volume => dev_gl_Cell_Volume
	use GEO_PARAM
	implicit none
	
	integer(4) :: gr  ! ���������� ����� ������� �����
	
    integer(4) :: st, s1, s2, i, j, k, zone, metod, now2, iter
    real(8) :: qqq1(9), qqq2(9), qqq(9)  ! ���������� � ������
    real(8) :: fluid1(5, 4), fluid2(5, 4)
    real(8) :: dist, dsl, dsc, dsp, distant(3)
    real(8) :: POTOK(9), ttest(3)
    real(8) :: POTOK_MF(5)
    real(8) :: POTOK_MF_all(5, 4)
    real(8) :: time, Volume, TT, U8, rad1, rad2, aa, bb, cc, wc
    real(8) :: SOURSE(5,5)  ! ��������� �����, �������� � ������� ��� ������ � ������� ����� ������������
    real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3
	
	logical :: null_bn
	
	
	time = 100000.0
	iter = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! ����� ������
	
	! ----------------------------------------------------------- ����� ������ ����������� ��� �� ����� ----------------------
	if (iter > size(gl_all_Gran_inner)) return
	
	gr = gl_all_Gran_inner(iter)
        !if(gl_Gran_info(gr) == 0) CYCLE
        POTOK = 0.0
        s1 = gl_Gran_neighbour(1, gr)
        s2 = gl_Gran_neighbour(2, gr)
        qqq1 = gl_Cell_par(:, s1)
        fluid1 = gl_Cell_par_MF(:, :, s1)   ! ��������� ��������� ��������� ��� ������������
		
		distant = gl_Gran_center(:, gr) - gl_Cell_center(:, s1)
		dist = norm2(distant)

        ! ��������� ������ ��������� ��������������� ��������
        if(norm2(qqq1(2:4))/sqrt(ggg*qqq1(5)/qqq1(1)) > 2.5) then
            rad1 = norm2(gl_Cell_center(:, s1))
            rad2 = norm2(gl_Gran_center(:, gr))
            qqq1(1) = qqq1(1) * rad1**2 / rad2**2
            qqq1(9) = qqq1(9) * rad1**2 / rad2**2
            qqq1(5) = qqq1(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
            ! �������� ������ � ����������� �.�.
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
            fluid2 = gl_Cell_par_MF(:, :, s2)   ! ��������� ��������� ��������� ��� ������������
            !dist = min(gl_Cell_dist(s1), gl_Cell_dist(s2))
			
			distant = gl_Gran_center(:, gr) - gl_Cell_center(:, s2)
			dist = min( dist, norm2(distant))

            ! ��������� ������ ��������� ��������������� ��������
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
			else
				write(*,*), "YTGVBNYTGBJYGBNJUYHBNJUHGBNMIUHKIJ"
			end if
	end if



        call chlld_Q(3, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
            0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)
        time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )   ! REDUCTION
        gl_Gran_POTOK(1:9, gr) = POTOK * gl_Gran_square(gr)
		gl_Gran_POTOK(10, gr) = 0.5 * DOT_PRODUCT(gl_Gran_normal(:, gr), qqq1(6:8) + qqq2(6:8)) * gl_Gran_square(gr)


        call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
            0.0_8, fluid1(:, 1), fluid2(:, 1), dsl, dsp, dsc, POTOK_MF)
        time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )
        gl_Gran_POTOK_MF(:, 1, gr) = POTOK_MF * gl_Gran_square(gr)

        call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
            0.0_8, fluid1(:, 2), fluid2(:, 2), dsl, dsp, dsc, POTOK_MF)
        time = min(time, 0.8 * dist/max(dabs(dsl), dabs(dsp)) )
        gl_Gran_POTOK_MF(:, 2, gr) = POTOK_MF * gl_Gran_square(gr)


        call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
            0.0_8, fluid1(:, 3), fluid2(:, 3), dsl, dsp, dsc, POTOK_MF)
        time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )
        gl_Gran_POTOK_MF(:, 3, gr) = POTOK_MF * gl_Gran_square(gr)

        call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
            0.0_8, fluid1(:, 4), fluid2(:, 4), dsl, dsp, dsc, POTOK_MF)
        time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )
        gl_Gran_POTOK_MF(:, 4, gr) = POTOK_MF * gl_Gran_square(gr)

	
	! ----------------------------------------------------------- �������� �� ����� ������� ----------------------
	
	if (dev_time_step_inner > time) then 
		time =  atomicmin(dev_time_step_inner, time)   ! ��������� �������� ������ ������������ ��������
	end if
	
	
	end subroutine  CUF_MGD_grans_MF_inner
	
	
	attributes(global) subroutine CUF_MGD_cells_MF_inner()
	! ��������� �����, ��� + �����������
	use MY_CUDA, gl_TS => dev_gl_TS, gl_Gran_neighbour => dev_gl_Gran_neighbour, gl_Gran_normal2 => dev_gl_Gran_normal2, &
		gl_Cell_par => dev_gl_Cell_par, gl_all_Gran => dev_gl_all_Gran, gl_Point_num => dev_gl_Point_num, gl_x2 => dev_gl_x2, &
		gl_y2 => dev_gl_y2, gl_z2 => dev_gl_z2, gl_Gran_center2 => dev_gl_Gran_center2, gl_Vx => dev_gl_Vx, &
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS, &
		gl_RAY_A => dev_gl_RAY_A, gl_RAY_B => dev_gl_RAY_B, gl_RAY_C => dev_gl_RAY_C, gl_RAY_O => dev_gl_RAY_O, &
		gl_RAY_K => dev_gl_RAY_K, gl_RAY_D => dev_gl_RAY_D, gl_RAY_E => dev_gl_RAY_E, norm2 => dev_norm2, par_n_TS => dev_par_n_TS, &
		par_n_HP => dev_par_n_HP, par_n_BS => dev_par_n_BS, par_n_END => dev_par_n_END, gl_Gran_square2 => dev_gl_Gran_square2, &
		gl_Cell_Volume2 => dev_gl_Cell_Volume2, gl_Cell_dist => dev_gl_Cell_dist, gl_all_Cell => dev_gl_all_Cell, gl_Cell_center2 => dev_gl_Cell_center2, &
		gl_Cell_gran => dev_gl_Cell_gran, gl_Cell_par_MF => dev_gl_Cell_par_MF, gl_Gran_type => dev_gl_Gran_type, &
		gl_Gran_POTOK => dev_gl_Gran_POTOK, gl_Gran_POTOK_MF => dev_gl_Gran_POTOK_MF, gl_Gran_info => dev_gl_Gran_info, gl_Cell_info => dev_gl_Cell_info, &
		Calc_sourse_MF => dev_Calc_sourse_MF, gl_all_Cell_inner => dev_gl_all_Cell_inner, gl_Gran_center => dev_gl_Gran_center, gl_Cell_center => dev_gl_Cell_center , &
		gl_Gran_normal => dev_gl_Gran_normal, gl_Gran_square => dev_gl_Gran_square, gl_Cell_type => dev_gl_Cell_type, &
		gl_Cell_number => dev_gl_Cell_number, gl_Cell_Volume => dev_gl_Cell_Volume
	use GEO_PARAM
	implicit none
	
	integer(4) :: gr
	
    integer(4) :: st, s1, s2, i, j, k, zone, iter
    real(8) :: qqq1(9), qqq2(9), qqq(9)  ! ���������� � ������
    real(8) :: fluid1(5, 4), fluid2(5, 4)
    real(8) :: dist, dsl, dsc, dsp
    real(8) :: POTOK(9), ttest(3)
    real(8) :: POTOK_MF(5)
    real(8) :: POTOK_MF_all(5, 4)
    real(8) :: Volume, TT, U8, rad1, rad2, aa, bb, cc, Volume2, sks
    real(8) :: SOURSE(5,5)  ! ��������� �����, �������� � ������� ��� ������ � ������� ����� ������������
    real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3
	logical :: l_1
	
	iter = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! ����� ������
	
	if (iter > size(gl_all_Cell_inner)) return
	
	! ----------------------------------------------------------- ����� ������ ����������� ��� �� ����� ----------------------
	
	gr = gl_all_Cell_inner(iter)
            !if(gl_Cell_info(gr) == 2) CYCLE
            l_1 = .TRUE.
            if ((gl_Cell_type(gr) == "A" .or. gl_Cell_type(gr) == "B").and.(gl_Cell_number(1, gr) <= 2) ) l_1 = .FALSE.    ! �� ������� � ������ ���� �������
            POTOK = 0.0
			sks = 0.0
            SOURSE = 0.0
            POTOK_MF_all = 0.0
            Volume = gl_Cell_Volume(gr)
            qqq = gl_Cell_par(:, gr)
            fluid1 = gl_Cell_par_MF(:, :, gr)
            ! ������������ ������ ����� �����
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


            ! ���������� ���� � ������� ���������

            zone = 1


            call Calc_sourse_MF(qqq, fluid1, SOURSE, zone)  ! ��������� ���������

            if (l_1 == .TRUE.) then
                ro3 = qqq(1) - dev_time_step_inner * POTOK(1) / Volume
                Q3 = qqq(9) - dev_time_step_inner * POTOK(9) / Volume
                if (ro3 <= 0.0_8) then
                    write(*, *) "Ro < 0  1490 "
                end if
                u3 = (qqq(1) * qqq(2) - dev_time_step_inner * (POTOK(2) + (qqq(6)/cpi4) * sks) / Volume + dev_time_step_inner * SOURSE(2, 1)) / ro3
                v3 = (qqq(1) * qqq(3) - dev_time_step_inner * (POTOK(3) + (qqq(7)/cpi4) * sks) / Volume + dev_time_step_inner * SOURSE(3, 1)) / ro3
                w3 = (qqq(1) * qqq(4) - dev_time_step_inner * (POTOK(4) + (qqq(8)/cpi4) * sks) / Volume + dev_time_step_inner * SOURSE(4, 1)) / ro3
                
				bx3 = qqq(6) - dev_time_step_inner * (POTOK(6) + qqq(2) * sks) / Volume
				by3 = qqq(7) - dev_time_step_inner * (POTOK(7) + qqq(3) * sks) / Volume
				bz3 = qqq(8) - dev_time_step_inner * (POTOK(8) + qqq(4) * sks) / Volume
				
				p3 = ((  ( qqq(5) / (ggg - 1.0) + 0.5 * qqq(1) * norm2(qqq(2:4))**2 + (qqq(6)**2 + qqq(7)**2 + qqq(8)**2) / 25.13274122871834590768 ) &
                    - dev_time_step_inner * ( POTOK(5) + (DOT_PRODUCT(qqq(2:4), qqq(6:8))/cpi4) * sks)/ Volume + dev_time_step_inner * SOURSE(5, 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) - (bx3**2 + by3**2 + bz3**2) / 25.13274122871834590768 ) * (ggg - 1.0)
				
				
				!p3 = ((  ( qqq(5) / (ggg - 1.0) + 0.5 * qqq(1) * norm2(qqq(2:4))**2 ) &
                !    - dev_time_step_inner * ( POTOK(5) + (DOT_PRODUCT(qqq(2:4), qqq(6:8))/cpi4) * sks)/ Volume + dev_time_step_inner * SOURSE(5, 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) ) * (ggg - 1.0)
				
				

                if (p3 <= 0.0_8) then
                    !print*, "p < 0  plasma 2028 ", p3 , gl_Cell_center(:, gr)
                    p3 = 0.000001
                    !pause
                end if

                gl_Cell_par(:, gr) = (/ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3/)

            end if

            ! ������ ��������� ������ ���������� ��� ��������� ���������

            do i = 1, 4
                if (i == 1 .and. l_1 == .FALSE.) CYCLE
                
                if (l_1 == .FALSE.) SOURSE(:, i + 1) = 0.0       ! �� ������������ �������� ������ ���� �����
                ro3 = fluid1(1, i) - dev_time_step_inner * POTOK_MF_all(1, i) / Volume + dev_time_step_inner * SOURSE(1, i + 1)
                if (ro3 <= 0.0_8) then
                    write(*,*) "Ro < 0  in sort ", i
                end if
                u3 = (fluid1(1, i) * fluid1(2, i) - dev_time_step_inner * POTOK_MF_all(2, i) / Volume + dev_time_step_inner * SOURSE(2, i + 1)) / ro3
                v3 = (fluid1(1, i) * fluid1(3, i) - dev_time_step_inner * POTOK_MF_all(3, i) / Volume + dev_time_step_inner * SOURSE(3, i + 1)) / ro3
                w3 = (fluid1(1, i) * fluid1(4, i) - dev_time_step_inner * POTOK_MF_all(4, i) / Volume + dev_time_step_inner * SOURSE(4, i + 1)) / ro3
                p3 = ((  ( fluid1(5, i) / (ggg - 1.0) + 0.5 * fluid1(1, i) * norm2(fluid1(2:4, i))**2 ) &
                    - dev_time_step_inner * POTOK_MF_all(5, i)/ Volume + dev_time_step_inner * SOURSE(5, i + 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) ) * (ggg - 1.0)
                if (p3 <= 0.0_8) then
                    p3 = 0.000001
                end if

                gl_Cell_par_MF(:, i, gr) = (/ro3, u3, v3, w3, p3/)
            end do
	
	
	end subroutine CUF_MGD_cells_MF_inner
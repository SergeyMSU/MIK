!******7**10**********************************************************70
!**                                                                  **
!**    chlld   for                                                   **
!**                                                                  **
!******7**10**********************************************************70
!**    shema HLLD(modification) for 3D MHD                           **
!******7**10**********************************************************70

	
	attributes(global) subroutine Cuda_Move_all_1(now)  ! Поверхностное натяжение на А - лучах
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
	
	real(8) :: Time, dist, r, rr, r1, r2, r3, r4, ddt, k1
    real(8) :: vel(3), Ak(3), Bk(3), Ck(3), Dk(3), Ek(3)
    integer :: i, j, k
	
	integer(4):: N2, N3
	integer(4) :: yzel, yzel2, yzel3, yzel22, yzel33, yzel222, yzel333
	
	  
	Time = time_step
	
	k = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	j = blockDim%y * (blockIdx%y - 1) + threadIdx%y   ! Номер потока
	
	N3 = size(gl_RAY_A(1, 1, :))
    N2 = size(gl_RAY_A(1, :, 1))
	
	if(k > N3 .or. j > N2) return
	
	k1 = 0.002! 0.001
	if(j <= 18) k1 = 0.03   ! 0.03
	
	
	if (j == 1) then
                    return
	end if
	
	ddt = Time/0.000127
			
	! Ударная волна
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
		yzel2 = gl_RAY_A(par_n_TS, j, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
	end if
	Dk(1) = gl_x2(yzel2, now); Dk(2) = gl_y2(yzel2, now); Dk(3) = gl_z2(yzel2, now)
			
	if(k < N3) then
		yzel2 = gl_RAY_A(par_n_TS, j, k + 1)
	else
		yzel2 = gl_RAY_A(par_n_TS, j, 1)
	end if
	Ek(1) = gl_x2(yzel2, now); Ek(2) = gl_y2(yzel2, now); Ek(3) = gl_z2(yzel2, now)
			
	if (gl_Point_num(yzel) > 0) then
		vel = gl_Point_num(yzel) * 0.25 * par_nat_TS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time * ddt
	else
		vel = 0.25 * par_nat_TS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time * ddt
	end if
		
			
	gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
	gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
	gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
	
	
	 ! Контакт
		yzel = gl_RAY_A(par_n_HP, j, k)
		Ak(1) = gl_x2(yzel, now); Ak(2) = gl_y2(yzel, now); Ak(3) = gl_z2(yzel, now)
		
		if (j < N2) then
			yzel2 = gl_RAY_A(par_n_HP, j + 1, k)
		else
			yzel2 = gl_RAY_B(par_n_HP, 1, k)
		end if
	
		if (j < N2 - 1) then
			yzel22 = gl_RAY_A(par_n_HP, j + 2, k)
		else if (j == N2 - 1) then
			yzel22 = gl_RAY_B(par_n_HP, 1, k)
		else
			yzel22 = gl_RAY_B(par_n_HP, 2, k)
		end if
	
		if (j < N2 - 2) then
			yzel222 = gl_RAY_A(par_n_HP, j + 3, k)
		else if (j == N2 - 2) then
			yzel222 = gl_RAY_B(par_n_HP, 1, k)
		else if (j == N2 - 1) then
			yzel222 = gl_RAY_B(par_n_HP, 2, k)
		else
			yzel222 = gl_RAY_B(par_n_HP, 3, k)
		end if
		
		if (j > 1) then
			yzel3 = gl_RAY_A(par_n_HP, j - 1, k)
		else
			yzel3 = gl_RAY_A(par_n_HP, 1, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
		end if
		
		if (j > 2) then
			yzel33 = gl_RAY_A(par_n_HP, j - 2, k)
		else if(j == 2) then
			yzel33 = gl_RAY_A(par_n_HP, 1, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
		else
			yzel33 = gl_RAY_A(par_n_HP, 2, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
		end if
	
		if (j > 3) then
			yzel333 = gl_RAY_A(par_n_HP, j - 3, k)
		else if(j == 3) then
			yzel333 = gl_RAY_A(par_n_HP, 1, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
		else if(j == 2) then
			yzel333 = gl_RAY_A(par_n_HP, 2, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
		else
			yzel333 = gl_RAY_A(par_n_HP, 3, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
		end if
	
		call Smooth_kvadr4(gl_x2(yzel333, now), gl_y2(yzel333, now), gl_z2(yzel333, now), &
			gl_x2(yzel33, now), gl_y2(yzel33, now), gl_z2(yzel33, now), &
			gl_x2(yzel3, now), gl_y2(yzel3, now), gl_z2(yzel3, now), &
			Ak(1), Ak(2), Ak(3), &
			gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now), &
			gl_x2(yzel22, now), gl_y2(yzel22, now), gl_z2(yzel22, now), &
			gl_x2(yzel222, now), gl_y2(yzel222, now), gl_z2(yzel222, now), &
			Ek(1), Ek(2), Ek(3))
	
	!call Smooth_kvadr3(gl_x2(yzel33, now), gl_y2(yzel33, now), gl_z2(yzel33, now), &
	!		gl_x2(yzel3, now), gl_y2(yzel3, now), gl_z2(yzel3, now), &
	!		Ak(1), Ak(2), Ak(3), &
	!		gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now), &
	!		gl_x2(yzel22, now), gl_y2(yzel22, now), gl_z2(yzel22, now), &
	!		Ek(1), Ek(2), Ek(3))
		
		!call Smooth_kvadr(gl_x2(yzel3, now), gl_y2(yzel3, now), gl_z2(yzel3, now), &
		!	gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now), &
		!	gl_x2(yzel22, now), gl_y2(yzel22, now), gl_z2(yzel22, now), &
		!	Ak(1), Ak(2), Ak(3), &
		!	Ek(1), Ek(2), Ek(3))
		!
		!call Smooth_kvadr(gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now), &
		!	gl_x2(yzel3, now), gl_y2(yzel3, now), gl_z2(yzel3, now), &
		!	gl_x2(yzel33, now), gl_y2(yzel33, now), gl_z2(yzel33, now), &
		!	Ak(1), Ak(2), Ak(3), &
		!	Bk(1), Bk(2), Bk(3))
		
		if(k < N3) then
			yzel2 = gl_RAY_A(par_n_HP, j, k + 1)
		else
			yzel2 = gl_RAY_A(par_n_HP, j, 1)
		end if
		
		if(k < N3 - 1) then
			yzel22 = gl_RAY_A(par_n_HP, j, k + 2)
		else if(k == N3 - 1) then
			yzel22 = gl_RAY_A(par_n_HP, j, 1)
		else
			yzel22 = gl_RAY_A(par_n_HP, j, 2)
		end if
	
		if(k < N3 - 2) then
			yzel222 = gl_RAY_A(par_n_HP, j, k + 2)
		else if(k == N3 - 2) then
			yzel222 = gl_RAY_A(par_n_HP, j, 1)
		else if(k == N3 - 1) then
			yzel222 = gl_RAY_A(par_n_HP, j, 2)
		else
			yzel222 = gl_RAY_A(par_n_HP, j, 3)
		end if
	
		if (k > 1) then
			yzel3 = gl_RAY_A(par_n_HP, j, k - 1)
		else
			yzel3 = gl_RAY_A(par_n_HP, j, N3)
		end if
	
		if (k > 2) then
			yzel33 = gl_RAY_A(par_n_HP, j, k - 2)
		else if(k == 2) then
			yzel33 = gl_RAY_A(par_n_HP, j, N3)
		else
			yzel33 = gl_RAY_A(par_n_HP, j, N3 - 1)
		end if
	
		if (k > 3) then
			yzel333 = gl_RAY_A(par_n_HP, j, k - 3)
		else if(k == 3) then
			yzel333 = gl_RAY_A(par_n_HP, j, N3)
		else if(k == 2) then
			yzel333 = gl_RAY_A(par_n_HP, j, N3 - 1)
		else
			yzel333 = gl_RAY_A(par_n_HP, j, N3 - 2)
		end if
	
		call Smooth_kvadr4(gl_x2(yzel333, now), gl_y2(yzel333, now), gl_z2(yzel333, now), &
			gl_x2(yzel33, now), gl_y2(yzel33, now), gl_z2(yzel33, now), &
			gl_x2(yzel3, now), gl_y2(yzel3, now), gl_z2(yzel3, now), &
			Ak(1), Ak(2), Ak(3), &
			gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now), &
			gl_x2(yzel22, now), gl_y2(yzel22, now), gl_z2(yzel22, now), &
			gl_x2(yzel222, now), gl_y2(yzel222, now), gl_z2(yzel222, now), &
			Ck(1), Ck(2), Ck(3))
	
	!call Smooth_kvadr3(gl_x2(yzel33, now), gl_y2(yzel33, now), gl_z2(yzel33, now), &
	!		gl_x2(yzel3, now), gl_y2(yzel3, now), gl_z2(yzel3, now), &
	!		Ak(1), Ak(2), Ak(3), &
	!		gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now), &
	!		gl_x2(yzel22, now), gl_y2(yzel22, now), gl_z2(yzel22, now), &
	!		Ck(1), Ck(2), Ck(3))
	
		!call Smooth_kvadr(gl_x2(yzel3, now), gl_y2(yzel3, now), gl_z2(yzel3, now), &
		!	gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now), &
		!	gl_x2(yzel22, now), gl_y2(yzel22, now), gl_z2(yzel22, now), &
		!	Ak(1), Ak(2), Ak(3), &
		!	Ck(1), Ck(2), Ck(3))
		!
		!call Smooth_kvadr(gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now), &
		!	gl_x2(yzel3, now), gl_y2(yzel3, now), gl_z2(yzel3, now), &
		!	gl_x2(yzel33, now), gl_y2(yzel33, now), gl_z2(yzel33, now), &
		!	Ak(1), Ak(2), Ak(3), &
		!	Dk(1), Dk(2), Dk(3))
		
		if (gl_Point_num(yzel) > 0) then
			!vel = par_nat_HP * 0.006 * (Bk/4.0 + Ck/4.0 + Dk/4.0 + Ek/4.0 - Ak) * gl_Point_num(yzel)/Time  ! 0.00006
			vel = par_nat_HP * 0.001 * (Ck/2.0 + Ek/2.0 - Ak) * gl_Point_num(yzel)/Time  ! 0.00006
		else
			!vel = par_nat_HP * 0.006 * (Bk/4.0 + Ck/4.0 + Dk/4.0 + Ek/4.0 - Ak)/Time   !  0.00006
			vel = par_nat_HP * 0.001 * (Ck/4.0 + Ek/2.0 - Ak)/Time   !  0.00006
		end if
			
	if(.False.) then ! Старая версия
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
			yzel2 = gl_RAY_A(par_n_HP, j, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
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
			!vel = gl_Point_num(yzel) * par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time/dist * ddt
		
			!vel = par_nat_HP * k1 * gl_Point_num(yzel) * ((Ak/r) * (rr - r)) * ddt  !0.035   0.0001
		
			vel = par_nat_HP * 0.001 * (Bk/4.0 + Ck/4.0 + Dk/4.0 + Ek/4.0 - Ak) * gl_Point_num(yzel)/Time  ! 0.00006
		else
			!vel = par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time/dist * ddt
		
			!vel = par_nat_HP * k1 * ((Ak/r) * (rr - r)) * ddt
		
			vel = par_nat_HP * 0.001 * (Bk/4.0 + Ck/4.0 + Dk/4.0 + Ek/4.0 - Ak)/Time   !  0.00006
		end if
	
		!if(j >= N2 - 8) then ! Натяжение на стыке А и Б лучей
		!	vel = vel + par_nat_HP * 0.006 * (Bk/2.0 + Dk/2.0 - Ak)/Time
		!end if
	end if
			
	gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
	gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
	gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
			
			
	! Внешняя ударная волна
			
	yzel = gl_RAY_A(par_n_BS, j, k)
	Ak(1) = gl_x2(yzel, now); Ak(2) = gl_y2(yzel, now); Ak(3) = gl_z2(yzel, now)
			
	if (j < N2) then
		yzel2 = gl_RAY_A(par_n_BS, j + 1, k)
	else
		return !yzel2 = yzel                                                    ! Вот тут есть неопределённость что делать с крайним узлом
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
		yzel2 = gl_RAY_A(par_n_BS, j, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
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
	
	attributes(global) subroutine Cuda_Move_all_2(now)   ! Поверхностное натяжение на лучах B
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
	
	real(8) :: Time, dist, r, rr, r1, r2, r3, r4, ddt, kkk
    real(8) :: vel(3), Ak(3), Bk(3), Ck(3), Dk(3), Ek(3)
    integer :: i, j, k
	
	integer(4):: N2, N3
	integer(4) :: yzel, yzel2, yzel3, yzel22, yzel33, yzel222, yzel333, yzel2222, yzel3333
	
	N3 = size(gl_RAY_B(1, 1, :))
    N2 = size(gl_RAY_B(1, :, 1))
	  
	Time = time_step
	
	k = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	j = blockDim%y * (blockIdx%y - 1) + threadIdx%y   ! Номер потока
	
	if(k > N3 .or. j > N2) return
	
	ddt = Time/0.000127
	! Ударная волна
			
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
			kkk = 0.2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!if(j < 15) kkk = 6.0
			! Контакт
			if (.True.) then !(j < N2 - 3) then
			yzel = gl_RAY_B(par_n_HP, j, k)
			Ak(1) = gl_x2(yzel, now); Ak(2) = gl_y2(yzel, now); Ak(3) = gl_z2(yzel, now)
			
			!if(Ak(2) < 0 .and. Ak(3) > 0) kkk = 6.0
		
			if (j < N2) then
				yzel2 = gl_RAY_B(par_n_HP, j + 1, k)
			else
				yzel2 = gl_RAY_O(1, 1, k)
			end if
	
			if (j < N2 - 1) then
				yzel22 = gl_RAY_B(par_n_HP, j + 2, k)
			else if (j == N2 - 1) then
				yzel22 = gl_RAY_O(1, 1, k)
			else
				yzel22 = gl_RAY_O(1, 2, k)
			end if
	
			if (j < N2 - 2) then
				yzel222 = gl_RAY_B(par_n_HP, j + 3, k)
			else if (j == N2 - 2) then
				yzel222 = gl_RAY_O(1, 1, k)
			else if (j == N2 - 1) then
				yzel222 = gl_RAY_O(1, 2, k)
			else
				yzel222 = gl_RAY_O(1, 3, k)
			end if
			
			if (j < N2 - 3) then
				yzel2222 = gl_RAY_B(par_n_HP, j + 4, k)
			else if (j == N2 - 3) then
				yzel2222 = gl_RAY_O(1, 1, k)
			else if (j == N2 - 2) then
				yzel2222 = gl_RAY_O(1, 2, k)
			else if (j == N2 - 1) then
				yzel2222 = gl_RAY_O(1, 3, k)
			else
				yzel2222 = gl_RAY_O(1, 4, k)
			end if
		
			if (j > 1) then
				yzel3 = gl_RAY_B(par_n_HP, j - 1, k)
			else
				yzel3 = gl_RAY_A(par_n_HP, size( gl_RAY_A(par_n_HP, :, k)), k)
			end if
		
			if (j > 2) then
				yzel33 = gl_RAY_B(par_n_HP, j - 2, k)
			else if(j == 2) then
				yzel33 = gl_RAY_A(par_n_HP, size( gl_RAY_A(par_n_HP, :, k)), k)
			else
				yzel33 = gl_RAY_A(par_n_HP, size( gl_RAY_A(par_n_HP, :, k)) - 1, k)
			end if
	
			if (j > 3) then
				yzel333 = gl_RAY_B(par_n_HP, j - 3, k)
			else if(j == 3) then
				yzel333 = gl_RAY_A(par_n_HP, size( gl_RAY_A(par_n_HP, :, k)), k)
			else if(j == 2) then
				yzel333 = gl_RAY_A(par_n_HP, size( gl_RAY_A(par_n_HP, :, k)) - 1, k)
			else
				yzel333 = gl_RAY_A(par_n_HP, size( gl_RAY_A(par_n_HP, :, k)) - 2, k)
			end if
			
			if (j > 4) then
				yzel3333 = gl_RAY_B(par_n_HP, j - 4, k)
			else if(j == 4) then
				yzel3333 = gl_RAY_A(par_n_HP, size( gl_RAY_A(par_n_HP, :, k)), k)
			else if(j == 3) then
				yzel3333 = gl_RAY_A(par_n_HP, size( gl_RAY_A(par_n_HP, :, k)) - 1, k)
			else if(j == 2) then
				yzel3333 = gl_RAY_A(par_n_HP, size( gl_RAY_A(par_n_HP, :, k)) - 2, k)
			else
				yzel3333 = gl_RAY_A(par_n_HP, size( gl_RAY_A(par_n_HP, :, k)) - 3, k)
			end if
			
			call Smooth_kvadr5(gl_x2(yzel3333, now), gl_y2(yzel3333, now), gl_z2(yzel3333, now), &
				gl_x2(yzel333, now), gl_y2(yzel333, now), gl_z2(yzel333, now), &
			gl_x2(yzel33, now), gl_y2(yzel33, now), gl_z2(yzel33, now), &
			gl_x2(yzel3, now), gl_y2(yzel3, now), gl_z2(yzel3, now), &
			Ak(1), Ak(2), Ak(3), &
			gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now), &
			gl_x2(yzel22, now), gl_y2(yzel22, now), gl_z2(yzel22, now), &
			gl_x2(yzel222, now), gl_y2(yzel222, now), gl_z2(yzel222, now), &
				gl_x2(yzel2222, now), gl_y2(yzel2222, now), gl_z2(yzel2222, now), &
			Ek(1), Ek(2), Ek(3))
	
		!call Smooth_kvadr4(gl_x2(yzel333, now), gl_y2(yzel333, now), gl_z2(yzel333, now), &
		!	gl_x2(yzel33, now), gl_y2(yzel33, now), gl_z2(yzel33, now), &
		!	gl_x2(yzel3, now), gl_y2(yzel3, now), gl_z2(yzel3, now), &
		!	Ak(1), Ak(2), Ak(3), &
		!	gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now), &
		!	gl_x2(yzel22, now), gl_y2(yzel22, now), gl_z2(yzel22, now), &
		!	gl_x2(yzel222, now), gl_y2(yzel222, now), gl_z2(yzel222, now), &
		!	Ek(1), Ek(2), Ek(3))
	
			!call Smooth_kvadr3(gl_x2(yzel33, now), gl_y2(yzel33, now), gl_z2(yzel33, now), &
			!gl_x2(yzel3, now), gl_y2(yzel3, now), gl_z2(yzel3, now), &
			!Ak(1), Ak(2), Ak(3), &
			!gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now), &
			!gl_x2(yzel22, now), gl_y2(yzel22, now), gl_z2(yzel22, now), &
			!Ek(1), Ek(2), Ek(3))
		
			!call Smooth_kvadr(gl_x2(yzel3, now), gl_y2(yzel3, now), gl_z2(yzel3, now), &
			!	gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now), &
			!	gl_x2(yzel22, now), gl_y2(yzel22, now), gl_z2(yzel22, now), &
			!	Ak(1), Ak(2), Ak(3), &
			!	Ek(1), Ek(2), Ek(3))
		 !
			!call Smooth_kvadr(gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now), &
			!	gl_x2(yzel3, now), gl_y2(yzel3, now), gl_z2(yzel3, now), &
			!	gl_x2(yzel33, now), gl_y2(yzel33, now), gl_z2(yzel33, now), &
			!	Ak(1), Ak(2), Ak(3), &
			!	Bk(1), Bk(2), Bk(3))
		
			if(k < N3) then
				yzel2 = gl_RAY_B(par_n_HP, j, k + 1)
			else
				yzel2 = gl_RAY_B(par_n_HP, j, 1)
			end if
		
			if(k < N3 - 1) then
				yzel22 = gl_RAY_B(par_n_HP, j, k + 2)
			else if(k == N3 - 1) then
				yzel22 = gl_RAY_B(par_n_HP, j, 1)
			else
				yzel22 = gl_RAY_B(par_n_HP, j, 2)
			end if
	
			if(k < N3 - 2) then
				yzel222 = gl_RAY_B(par_n_HP, j, k + 2)
			else if(k == N3 - 2) then
				yzel222 = gl_RAY_B(par_n_HP, j, 1)
			else if(k == N3 - 1) then
				yzel222 = gl_RAY_B(par_n_HP, j, 2)
			else
				yzel222 = gl_RAY_B(par_n_HP, j, 3)
			end if
	
			if (k > 1) then
				yzel3 = gl_RAY_B(par_n_HP, j, k - 1)
			else
				yzel3 = gl_RAY_B(par_n_HP, j, N3)
			end if
	
			if (k > 2) then
				yzel33 = gl_RAY_B(par_n_HP, j, k - 2)
			else if(k == 2) then
				yzel33 = gl_RAY_B(par_n_HP, j, N3)
			else
				yzel33 = gl_RAY_B(par_n_HP, j, N3 - 1)
		end if
	
		if (k > 3) then
			yzel333 = gl_RAY_B(par_n_HP, j, k - 3)
		else if(k == 3) then
			yzel333 = gl_RAY_B(par_n_HP, j, N3)
		else if(k == 2) then
			yzel333 = gl_RAY_B(par_n_HP, j, N3 - 1)
		else
			yzel333 = gl_RAY_B(par_n_HP, j, N3 - 2)
		end if
	
		call Smooth_kvadr4(gl_x2(yzel333, now), gl_y2(yzel333, now), gl_z2(yzel333, now), &
			gl_x2(yzel33, now), gl_y2(yzel33, now), gl_z2(yzel33, now), &
			gl_x2(yzel3, now), gl_y2(yzel3, now), gl_z2(yzel3, now), &
			Ak(1), Ak(2), Ak(3), &
			gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now), &
			gl_x2(yzel22, now), gl_y2(yzel22, now), gl_z2(yzel22, now), &
			gl_x2(yzel222, now), gl_y2(yzel222, now), gl_z2(yzel222, now), &
			Ck(1), Ck(2), Ck(3))
	
			!call Smooth_kvadr3(gl_x2(yzel33, now), gl_y2(yzel33, now), gl_z2(yzel33, now), &
			!		gl_x2(yzel3, now), gl_y2(yzel3, now), gl_z2(yzel3, now), &
			!		Ak(1), Ak(2), Ak(3), &
			!		gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now), &
			!		gl_x2(yzel22, now), gl_y2(yzel22, now), gl_z2(yzel22, now), &
			!		Ck(1), Ck(2), Ck(3))
	
			!call Smooth_kvadr(gl_x2(yzel3, now), gl_y2(yzel3, now), gl_z2(yzel3, now), &
			!	gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now), &
			!	gl_x2(yzel22, now), gl_y2(yzel22, now), gl_z2(yzel22, now), &
			!	Ak(1), Ak(2), Ak(3), &
			!	Ck(1), Ck(2), Ck(3))
		 !
			!call Smooth_kvadr(gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now), &
			!	gl_x2(yzel3, now), gl_y2(yzel3, now), gl_z2(yzel3, now), &
			!	gl_x2(yzel33, now), gl_y2(yzel33, now), gl_z2(yzel33, now), &
			!	Ak(1), Ak(2), Ak(3), &
			!	Dk(1), Dk(2), Dk(3))
		
			if (gl_Point_num(yzel) > 0) then
				!vel = par_nat_HP * 0.006 * (Bk/4.0 + Ck/4.0 + Dk/4.0 + Ek/4.0 - Ak) * gl_Point_num(yzel)/Time  ! 0.00006
				vel = kkk * par_nat_HP * 0.0002 * (Ck/2.0 + Ek/2.0 - Ak) * gl_Point_num(yzel)/Time  ! 0.00006
			else
				!vel = par_nat_HP * 0.006 * (Bk/4.0 + Ck/4.0 + Dk/4.0 + Ek/4.0 - Ak)/Time   !  0.00006
				vel = kkk * par_nat_HP * 0.0002 * (Ck/2.0 + Ek/2.0 - Ak)/Time   !  0.00006
			end if
			
			else  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				yzel = gl_RAY_B(par_n_HP, j, k)
				Ak(1) = gl_x2(yzel, now); Ak(2) = gl_y2(yzel, now); Ak(3) = gl_z2(yzel, now)
				! Ak = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/)
				r = sqrt(Ak(2)**2 + Ak(3)**2)
			
				if (j < N2) then
					yzel2 = gl_RAY_B(par_n_HP, j + 1, k)
				else
					yzel2 = gl_RAY_O(1, 1, k)
					!return
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
				!rr = (r1 + r3)/2.0
			
				!dist = sqrt( (Dk(1) - Ak(1))**2 + (Dk(2) - Ak(2))**2 + (Dk(3) - Ak(3))**2 )
				!dist = max(dist, 1.0_8)
			
				if (gl_Point_num(yzel) > 0) then
					!vel = gl_Point_num(yzel) * par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time/dist
				
					!vel = par_nat_HP * 0.0001 * gl_Point_num(yzel) * ((Ak/r) * (rr - r)) * ddt  !0.001
					!vel(1) = 0.0
				
					vel = par_nat_HP * 0.00001 * (Bk/4.0 + Ck/4.0 + Dk/4.0 + Ek/4.0 - Ak) * gl_Point_num(yzel)/Time  !0.000003 
				else
					!vel = par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time/dist
				
					!vel = par_nat_HP * 0.0001 * ((Ak/r) * (rr - r)) * ddt
					!vel(1) = 0.0
		
					vel = par_nat_HP * 0.00001 * (Bk/4.0 + Ck/4.0 + Dk/4.0 + Ek/4.0 - Ak)/Time  !0.000003
				end if
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
	
	
	attributes(global) subroutine Cuda_Move_all_3(now)   ! Поверхностное натяжение на лучах K
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
	
	k = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	j = blockDim%y * (blockIdx%y - 1) + threadIdx%y   ! Номер потока
	
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
		yzel2 = gl_RAY_K(par_n_TS, j, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
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
	
	
	
	
	
	
	attributes(global) subroutine Cuda_Move_all_4(now)   ! Поверхностное натяжение на лучах O 
	use MY_CUDA, gl_TS => dev_gl_TS, gl_Gran_neighbour => dev_gl_Gran_neighbour, gl_Gran_normal2 => dev_gl_Gran_normal2, &
		gl_Cell_par => dev_gl_Cell_par, gl_all_Gran => dev_gl_all_Gran, gl_Point_num => dev_gl_Point_num, gl_x2 => dev_gl_x2, &
		gl_y2 => dev_gl_y2, gl_z2 => dev_gl_z2, gl_Gran_center2 => dev_gl_Gran_center2, gl_Vx => dev_gl_Vx, &
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS, &
		gl_RAY_A => dev_gl_RAY_A, gl_RAY_B => dev_gl_RAY_B, gl_RAY_C => dev_gl_RAY_C, gl_RAY_O => dev_gl_RAY_O, &
		gl_RAY_K => dev_gl_RAY_K, gl_RAY_D => dev_gl_RAY_D, gl_RAY_E => dev_gl_RAY_E, norm2 => dev_norm2, par_n_TS => dev_par_n_TS, &
		par_n_HP => dev_par_n_HP, par_n_BS => dev_par_n_BS, gl_Vn => dev_gl_Vn
	use GEO_PARAM
	use cudafor
	
	implicit none
	integer, intent(in) :: now
	
	real(8) :: Time, dist,  r, rr, r1, r2, r3, r4, ddt, nd, nd2, kk
    real(8) :: vel(3), Ak(3), Bk(3), Ck(3), Dk(3), Ek(3)
    integer :: i, j, k
	
	integer(4):: N2, N3
	integer(4) :: yzel, yzel2
	
	
	N3 = size(gl_RAY_O(1, 1, :))
    N2 = size(gl_RAY_O(1, :, 1))
	  
	Time = time_step
	
	k = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	j = blockDim%y * (blockIdx%y - 1) + threadIdx%y   ! Номер потока
	
	if(k > N3 .or. j > N2) return
	
	!if (k /= 1 .and. j == 1) then
 !           return
	!end if
	
	ddt = Time/0.000127
			
	! Контакт
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
		yzel2 = gl_RAY_C(1, size(gl_RAY_C(1, :, k)), k)
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
	
	!rr = (r1 + r3)/2.0
	
	
	
			
	if (gl_Point_num(yzel) > 0) then
		!vel = gl_Point_num(yzel) * par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time
		
		!vel = par_nat_HP * 0.0001 * gl_Point_num(yzel) * ((Ak/r) * (rr - r)) * ddt   ! 0.0003
		!vel(1) = 0.0
		
		vel = par_nat_HP * 0.00006 * (Bk/4.0 + Ck/4.0 + Dk/4.0 + Ek/4.0 - Ak) * gl_Point_num(yzel)/Time  ! надо ещё уменьшать
		
	else
		!vel = par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time
		
		!vel = par_nat_HP * 0.0001 * ((Ak/r) * (rr - r)) * ddt
		!vel(1) = 0.0
		
		vel = par_nat_HP * 0.00006 * (Bk/4.0 + Ck/4.0 + Dk/4.0 + Ek/4.0 - Ak)/Time  !0.000002
	end if
	
	!Ak = Bk/4.0 + Ck/4.0 + Dk/4.0 + Ek/4.0 - Ak
	!nd = norm2(Ak)
	!vel = 10.0 * gl_Vn * (Ak)/nd
			
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
	
	attributes(global) subroutine Cuda_Move_all_5(now)   ! Отдельное поверхностное натяжение для точки на оси симметрии
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
	
	real(8) :: Time, r, r2
    real(8) :: Ak(3), Bk(3)
    integer :: i, j, k
	
	integer(4) :: yzel, yzel2
	
	
	  
	Time = time_step
	
	i = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	
	if (i > 1) return
	
	r2 = 0.0
	k = size(gl_RAY_A(1, 1, :))
	yzel = gl_RAY_A(par_n_TS, 1, 1)
	Bk(1) = gl_x2(yzel, now)
	Bk(2) = gl_y2(yzel, now)
	Bk(3) = gl_z2(yzel, now)
	
	r = norm2(Bk)
	
	do j = 1, k
		yzel2 = gl_RAY_A(par_n_TS, 2, k)
		Ak(1) = gl_x2(yzel2, now)
		Ak(2) = gl_y2(yzel2, now)
		Ak(3) = gl_z2(yzel2, now)
		r2 = r2 + norm2(Ak)
	end do
	
	r2 = r2/k
	
	gl_Vx(yzel) = (r2 - r)/Time
	gl_Vy(yzel) = 0.0
	gl_Vz(yzel) = 0.0
	
	gl_Point_num(yzel) = 1
	
	end subroutine Cuda_Move_all_5
	
	attributes(global) subroutine Cuda_Move_all_A(now)   
	! Движение точек на А - лучах
	use MY_CUDA, gl_TS => dev_gl_TS, gl_Gran_neighbour => dev_gl_Gran_neighbour, gl_Gran_normal2 => dev_gl_Gran_normal2, &
		gl_Cell_par => dev_gl_Cell_par, gl_all_Gran => dev_gl_all_Gran, gl_Point_num => dev_gl_Point_num, gl_x2 => dev_gl_x2, &
		gl_y2 => dev_gl_y2, gl_z2 => dev_gl_z2, gl_Gran_center2 => dev_gl_Gran_center2, gl_Vx => dev_gl_Vx, &
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS, &
		gl_RAY_A => dev_gl_RAY_A, gl_RAY_B => dev_gl_RAY_B, gl_RAY_C => dev_gl_RAY_C, gl_RAY_O => dev_gl_RAY_O, &
		gl_RAY_K => dev_gl_RAY_K, gl_RAY_D => dev_gl_RAY_D, gl_RAY_E => dev_gl_RAY_E, norm2 => dev_norm2, par_n_TS => dev_par_n_TS, &
		par_n_HP => dev_par_n_HP, par_n_BS => dev_par_n_BS, par_n_END => dev_par_n_END, &
		par_n_IA => dev_par_n_IA, par_n_IB => dev_par_n_IB, par_R_inner => dev_par_R_inner, &
		par_kk1 => dev_par_kk1, par_kk2 => dev_par_kk2, par_R_END => dev_par_R_END, par_kk12 => dev_par_kk12, &
		par_kk13 => dev_par_kk13
	use GEO_PARAM
	use cudafor
	
	implicit none
    integer, intent(in) :: now
	
	real(8), device :: Time
    real(8), device :: vel(3), ER(3), KORD(3), ER2(3)
	real(8) :: R_TS, proect, R_HP, R_BS
    integer:: i, j, k
    real(8), device :: the, phi, r, kk13, rd
    integer :: now2             ! Эти параметры мы сейчас меняем на основе now
	
	integer(4):: N1, N2, N3
	integer(4) :: yzel
	
    now2 = mod(now, 2) + 1  
	Time = time_step
	
	N3 = size(gl_RAY_A(1, 1, :))
    N2 = size(gl_RAY_A(1, :, 1))
	N1 = size(gl_RAY_A(:, 1, 1))
	  
	
	k = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	j = blockDim%y * (blockIdx%y - 1) + threadIdx%y   ! Номер потока
	
	if(k > N3 .or. j > N2) return
	
	if (k /= 1 .and. j == 1) then
                    return
            end if
            
        ! Вычисляем координаты текущего луча в пространстве
        the = (j - 1) * par_pi_8/2.0/(N2 - 1)
        phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), dev_par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
            
        ! TS
        yzel = gl_RAY_A(par_n_TS, j, k)
        if(gl_Point_num(yzel) == 0) then
            vel = 0.0
        else
            ! vel = (/gl_Vx(yzel), gl_Vy(yzel), gl_Vz(yzel)/)
			vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
            vel = vel/gl_Point_num(yzel)                       ! Нашли скорость движения этого узла
        end if
            
        ! Обнулим для следующего использования
        gl_Point_num(yzel) = 0
        gl_Vx(yzel) = 0.0
        gl_Vy(yzel) = 0.0
        gl_Vz(yzel) = 0.0
            
		ER(1) = cos(the); ER(2) = sin(the) * cos(phi); ER(3) = sin(the) * sin(phi)
		proect = DOT_PRODUCT(vel * Time, ER)
        ! proect = DOT_PRODUCT(vel * Time, (/cos(the), sin(the) * cos(phi), sin(the) * sin(phi)/))  !  Находим проекцию перемещения на радиус вектор луча
        KORD(1) = gl_x2(yzel, now); KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now) 
		ER2 = proect * ER
		ER2  = KORD + ER2
		R_TS = norm2(ER2)  ! Новое расстояние до TS
		
		!if(j == 1 .and. k==1)  write(*,*) "=    ", gl_x2(yzel, now)
            
        ! HP
        yzel = gl_RAY_A(par_n_HP, j, k)
        if(gl_Point_num(yzel) == 0) then
            vel = 0.0
        else
            ! vel = (/gl_Vx(yzel), gl_Vy(yzel), gl_Vz(yzel)/)
			vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
            vel = vel/gl_Point_num(yzel)                       ! Нашли скорость движения этого узла
	end if
	
	
            
        ! Обнулим для следующего использования
        gl_Point_num(yzel) = 0
        gl_Vx(yzel) = 0.0
        gl_Vy(yzel) = 0.0
        gl_Vz(yzel) = 0.0
            
        proect = DOT_PRODUCT(vel * Time, ER)
		KORD(1) = gl_x2(yzel, now); KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now) 
		ER2 = proect * ER
		ER2  = KORD + ER2
        R_HP = norm2(ER2)  ! Новое расстояние до HP
		
		!if(j == 1 .and. k == 1) then
		!		write(*,*) yzel, now
		!		write(*,*) R_HP
		!		write(*,*) vel(1), vel(2), vel(3)
		!		write(*,*) Time, proect
		!		write(*,*) ER(1), ER(2), ER(3)
		!		write(*, *) KORD(1), KORD(2), KORD(3)
		!		write(*,*) "______________________"
		!end if
            
        ! BS
        yzel = gl_RAY_A(par_n_BS, j, k)
        if(gl_Point_num(yzel) == 0) then
            vel = 0.0
        else
            ! vel = (/gl_Vx(yzel), gl_Vy(yzel), gl_Vz(yzel)/)
			vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
            vel = vel/gl_Point_num(yzel)                       ! Нашли скорость движения этого узла
		end if
	
            
        ! Обнулим для следующего успользования
        gl_Point_num(yzel) = 0
        gl_Vx(yzel) = 0.0
        gl_Vy(yzel) = 0.0
        gl_Vz(yzel) = 0.0
            
        proect = DOT_PRODUCT(vel * Time, ER)
		KORD(1) = gl_x2(yzel, now); KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now) 
		ER2 = proect * ER
		ER2  = KORD + ER2
        R_BS = norm2(ER2)  ! Новое расстояние до BS
		
            
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
	! Движение точек на B - лучах
	use MY_CUDA, gl_TS => dev_gl_TS, gl_Gran_neighbour => dev_gl_Gran_neighbour, gl_Gran_normal2 => dev_gl_Gran_normal2, &
		gl_Cell_par => dev_gl_Cell_par, gl_all_Gran => dev_gl_all_Gran, gl_Point_num => dev_gl_Point_num, gl_x2 => dev_gl_x2, &
		gl_y2 => dev_gl_y2, gl_z2 => dev_gl_z2, gl_Gran_center2 => dev_gl_Gran_center2, gl_Vx => dev_gl_Vx, &
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS, &
		gl_RAY_A => dev_gl_RAY_A, gl_RAY_B => dev_gl_RAY_B, gl_RAY_C => dev_gl_RAY_C, gl_RAY_O => dev_gl_RAY_O, &
		gl_RAY_K => dev_gl_RAY_K, gl_RAY_D => dev_gl_RAY_D, gl_RAY_E => dev_gl_RAY_E, norm2 => dev_norm2, par_n_TS => dev_par_n_TS, &
		par_n_HP => dev_par_n_HP, par_n_BS => dev_par_n_BS, par_n_END => dev_par_n_END, par_triple_point => dev_par_triple_point, &
		par_n_IA => dev_par_n_IA, par_n_IB => dev_par_n_IB, par_R_inner => dev_par_R_inner, &
		par_kk1 => dev_par_kk1, par_kk12 => dev_par_kk12, par_triple_point_2 => dev_par_triple_point_2
	use GEO_PARAM
	use cudafor
	
	implicit none
    integer, intent(in) :: now
	
	real(8), device :: Time
    real(8), device :: vel(3), ER(3), ER2(3), KORD(3), ER3(3)
	real(8) :: R_TS, proect, R_HP, R_BS
    integer:: i, j, k
    real(8), device :: the, phi, r, r1, the2
    integer :: now2             ! Эти параметры мы сейчас меняем на основе now
	
	integer(4):: N1, N2, N3
	integer(4) :: yzel
	
    now2 = mod(now, 2) + 1  
	Time = time_step
	
	N3 = size(gl_RAY_B(1, 1, :))
    N2 = size(gl_RAY_B(1, :, 1))
	N1 = size(gl_RAY_B(:, 1, 1))
	  
	
	k = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	j = blockDim%y * (blockIdx%y - 1) + threadIdx%y   ! Номер потока
	
	if(k > N3 .or. j > N2) return
	
	! Вычисляем координаты текущего луча в пространстве
            the = par_pi_8/2.0 + (j) * par_triple_point/(N2)
			the2 = par_pi_8/2.0 + (j) * par_triple_point_2/(N2)
            phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), dev_par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
            
            ! TS
            yzel = gl_RAY_B(par_n_TS, j, k)
            if(gl_Point_num(yzel) == 0) then
                vel = 0.0
			else
				vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
                vel = vel/gl_Point_num(yzel)                       ! Нашли скорость движения этого узла
            end if
            
            ! Обнулим для следующего использования
            gl_Point_num(yzel) = 0
            gl_Vx(yzel) = 0.0
            gl_Vy(yzel) = 0.0
            gl_Vz(yzel) = 0.0
            
			ER(1) = cos(the); ER(2) = sin(the) * cos(phi); ER(3) = sin(the) * sin(phi)
            proect = DOT_PRODUCT(vel * Time, ER)  !  Находим проекцию перемещения на радиус вектор луча
			KORD(1) = gl_x2(yzel, now); KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now) 
			ER3 = KORD + proect * ER
            R_TS = norm2(ER3)  ! Новое расстояние до TS
            
            ! HP
            yzel = gl_RAY_B(par_n_HP, j, k)
            if(gl_Point_num(yzel) == 0) then
                vel = 0.0
			else
				vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
                vel = vel/gl_Point_num(yzel)                       ! Нашли скорость движения этого узла
			end if

            
            ! Обнулим для следующего использования
            gl_Point_num(yzel) = 0
            gl_Vx(yzel) = 0.0
            gl_Vy(yzel) = 0.0
            gl_Vz(yzel) = 0.0
            
			ER2(1) = cos(the2); ER2(2) = sin(the2) * cos(phi); ER2(3) = sin(the2) * sin(phi)
            proect = DOT_PRODUCT(vel * Time, ER2)
			KORD(1) = gl_x2(yzel, now); KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now)
			ER3 = KORD + proect * ER2 - R_TS * ER
            R_HP = norm2(ER3)  ! Новое расстояние до HP
			
                
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
					
                    r =  (i - par_n_TS) * (R_HP) /(par_n_HP - par_n_TS)
					gl_x2(yzel, now2) = R_TS * cos(the) + r * cos(the2)
					gl_y2(yzel, now2) = R_TS * sin(the) * cos(phi) + r * sin(the2) * cos(phi)
					gl_z2(yzel, now2) = R_TS * sin(the) * sin(phi) + r * sin(the2) * sin(phi)
					CYCLE
				end if
				
				if (i == par_n_TS - 1) then
					r = R_TS - 0.5               
				end if

                ! Записываем новые координаты
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
	! Движение точек на C - лучах
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
    real(8), device :: the, phi, r,  x, y, z, rr, rd
    integer :: now2             ! Эти параметры мы сейчас меняем на основе now
	
	integer(4):: N1, N2, N3
	integer(4) :: yzel
	
    now2 = mod(now, 2) + 1  
	Time = time_step
	
	N3 = size(gl_RAY_C(1, 1, :))
    N2 = size(gl_RAY_C(1, :, 1))
	N1 = size(gl_RAY_C(:, 1, 1))
	  
	
	k = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	j = blockDim%y * (blockIdx%y - 1) + threadIdx%y   ! Номер потока
	
	if(k > N3 .or. j > N2) return
	
	
	! Вычисляем координаты текущего луча в пространстве
                phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), dev_par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
                
                ! Вычисляем координаты точки на луче
                x = gl_x2(gl_RAY_B(par_n_HP, j, k), now2)
                y = gl_y2(gl_RAY_B(par_n_HP, j, k), now2)
                z = gl_z2(gl_RAY_B(par_n_HP, j, k), now2)
                rr = (y**2 + z**2)**(0.5)
				
				y = gl_y2(gl_RAY_B(par_n_HP - 1, j, k), now2)
                z = gl_z2(gl_RAY_B(par_n_HP - 1, j, k), now2)
                rd = (y**2 + z**2)**(0.5)
				rr = rr + (rr - rd)
                
                ! BS     Нужно взять положение BS из её положения на крайнем луче A
                yzel = gl_RAY_A(par_n_BS, size(gl_RAY_A(1, :, 1)), k)
				ER(1) = 0.0_8; ER(2) = gl_y2(yzel, now2); ER(3) = gl_z2(yzel, now2)
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
	! Движение точек на O - лучах
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
    integer :: now2             ! Эти параметры мы сейчас меняем на основе now
	
	integer(4):: N1, N2, N3
	integer(4) :: yzel
	
    now2 = mod(now, 2) + 1  
	Time = time_step
	
	N3 = size(gl_RAY_O(1, 1, :))
    N2 = size(gl_RAY_O(1, :, 1))
	N1 = size(gl_RAY_O(:, 1, 1))
	  
	
	k = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	j = blockDim%y * (blockIdx%y - 1) + threadIdx%y   ! Номер потока
	
	if(k > N3 .or. j > N2) return
	
	phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), dev_par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
            yzel = gl_RAY_O(1, j, k)
            if(gl_Point_num(yzel) == 0) then
                vel = 0.0
			else
				vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
                vel = vel/gl_Point_num(yzel)                       ! Нашли скорость движения этого узла
            end if
            
            ! Обнулим для следующего использования
            gl_Point_num(yzel) = 0
            gl_Vx(yzel) = 0.0_8
            gl_Vy(yzel) = 0.0_8
            gl_Vz(yzel) = 0.0_8
            
			ER(1) = 0.0_8; ER(2) = cos(phi); ER(3) = sin(phi)
            proect = DOT_PRODUCT(vel * Time, ER)  !  Находим проекцию перемещения на радиус вектор луча
			KORD(1) = 0.0_8; KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now)
			ER2 = proect * ER
			ER2  = KORD + ER2
            R_HP = norm2(ER2)  ! Новое расстояние до HP
			
			! Блокируем схлопывание контакта к оси  20.0
			if(R_HP < 17.0_8) then
				R_HP = 17.0_8
			end if
            
            xx = gl_x2(gl_RAY_B(par_n_HP, par_m_BC, k), now2)              ! Отталкиваемся от x - координаты крайней точки B на гелиопаузе в этой плоскости (k)
            x = xx - (DBLE(j)/N2)**par_kk31 * (xx - par_R_LEFT)
            
            ! BS     Нужно взять положение BS из её положения на крайнем луче A
            yzel = gl_RAY_A(par_n_BS, size(gl_RAY_A(1, :, 1)), k)
			KORD(1) = 0.0_8; KORD(2) = gl_y2(yzel, now2); KORD(3) = gl_z2(yzel, now2)
            R_BS = norm2(KORD)  ! Новое расстояние до BS
			
			if (R_BS < R_HP + 10.0) then
				R_BS = R_HP + 10.0
			end if
            
            do i = 1, N1
                yzel = gl_RAY_O(i, j, k)
                
                if (i <= par_n_BS - par_n_HP + 1) then
                    r = R_HP + (i - 1) * (R_BS - R_HP)/(par_n_BS - par_n_HP)
                else
                    r = R_BS + (DBLE(i - (par_n_BS - par_n_HP + 1))/(N1 - (par_n_BS - par_n_HP + 1) ))**(0.55 * par_kk2) * (par_R_END - R_BS)
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
	! Движение точек на B - лучах
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
    integer :: now2             ! Эти параметры мы сейчас меняем на основе now
	
	integer(4):: N1, N2, N3
	integer(4) :: yzel
	
    now2 = mod(now, 2) + 1  
	Time = time_step
	
	N3 = size(gl_RAY_K(1, 1, :))
    N2 = size(gl_RAY_K(1, :, 1))
	N1 = size(gl_RAY_K(:, 1, 1))
	  
	
	k = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	j = blockDim%y * (blockIdx%y - 1) + threadIdx%y   ! Номер потока
	
	if(k > N3 .or. j > N2) return
	if (k /= 1 .and. j == 1) return
	
	! Вычисляем координаты текущего луча в пространстве
    the = par_pi_8/2.0 + par_triple_point + (N2 - j + 1) * (par_pi_8/2.0 - par_triple_point)/(N2)
    phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), dev_par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
            
            
    yzel = gl_RAY_K(par_n_TS, j, k)
    if(gl_Point_num(yzel) == 0) then
        vel = 0.0
	else
		vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
        vel = vel/gl_Point_num(yzel)                       ! Нашли скорость движения этого узла
    end if
            
    ! Обнулим для следующего использования
    gl_Point_num(yzel) = 0
    gl_Vx(yzel) = 0.0
    gl_Vy(yzel) = 0.0
    gl_Vz(yzel) = 0.0
            
	ER(1) = cos(the); ER(2) = sin(the) * cos(phi); ER(3) = sin(the) * sin(phi)
    proect = DOT_PRODUCT(vel * Time, ER)  !  Находим проекцию перемещения на радиус вектор луча
    KORD(1) = gl_x2(yzel, now); KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now)
	ER2 = proect * ER
	ER2  = KORD + ER2
	R_TS = norm2(ER2)  ! Новое расстояние до TS
            
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
					r = R_TS - 0.5               ! UBRAT
		end if
				
        !r =  par_R0 + (R_TS - par_R0) * (REAL(i, KIND = 4)/par_n_TS)**par_kk1


        ! Записываем новые координаты
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
	! Движение точек на D - лучах
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
    integer :: now2             ! Эти параметры мы сейчас меняем на основе now
	
	integer(4):: N1, N2, N3
	integer(4) :: yzel
	
    now2 = mod(now, 2) + 1  
	Time = time_step
	
	N3 = size(gl_RAY_D(1, 1, :))
    N2 = size(gl_RAY_D(1, :, 1))
	N1 = size(gl_RAY_D(:, 1, 1))
	  
	
	k = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	j = blockDim%y * (blockIdx%y - 1) + threadIdx%y   ! Номер потока
	
	if(k > N3 .or. j > N2) return
	
	if (k /= 1 .and. j == 1) return
	
	
	!rrr = 100.0
	
	! Вычисляем координаты текущего луча в пространстве
                phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), dev_par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
				
				! Вычисляем координаты точки на луче

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

                ! Вычисляем координаты точки на луче

				!y = gl_y2(gl_RAY_O(1, i - 1, k), now2)
				!z = gl_z2(gl_RAY_O(1, i - 1, k), now2)
				!!
				!rr = sqrt(y**2 + z**2)
				!
				!
				!rrr = min(rrr, (rr/20.0))
				!if (rr < 20.0) then
				!	rr = r * rrr**2              ! Сдвигаем немного ниже
				!else
				!	rr = r
				!end if
				
				

                yzel = gl_RAY_D(i, j, k)
                
                
                x = xx + (DBLE(i - 1)/(N1 - 1))**par_kk3 * (par_R_LEFT - xx)
				
				!rrr = min(r, rr) * (DBLE(j - 1)/(N2 - 1))
				
				! rrr = r * (1.0 - (DBLE(i - 1)**2/(N1 - 1)**2)) + (rr * (DBLE(j - 1)/(N2 - 1))) * ( (DBLE(i - 1)**2) /(N1 - 1)**2)


                ! Записываем новые координаты
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
	! Движение точек на E - лучах
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
    real(8), device :: the, phi, r, x, y, z, x2, y2, z2, rd
    integer :: now2             ! Эти параметры мы сейчас меняем на основе now
	
	integer(4):: N1, N2, N3
	integer(4) :: yzel
	
    now2 = mod(now, 2) + 1  
	Time = time_step
	
	N3 = size(gl_RAY_E(1, 1, :))
    N2 = size(gl_RAY_E(1, :, 1))
	N1 = size(gl_RAY_E(:, 1, 1))
	  
	
	k = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	j = blockDim%y * (blockIdx%y - 1) + threadIdx%y   ! Номер потока
	
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
	! Движение точек на E - лучах
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
    ! Цикл по граням - считаем площадь грани, её нормаль
    Ngran = size(gl_all_Gran(1,:))
	
	iter = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	
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

    ! Нужно сохранить центр грани
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


    ! Нужно записать площадь грани и нормаль в общий массив!
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
	
	iter = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	
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


        ! Вычислили центр ячейки теперь считаем объём пирамиды на каждую грань
        do j = 1, 6
            i = gl_Cell_gran(j, iter)   ! Берём по очереди все грани ячейки

            if (i == 0) CYCLE
            k = gl_all_Gran(1, i)  ! Номер первого узла грани
			b(1) = gl_x2(k, now)
			b(2) = gl_y2(k, now)
			b(3) = gl_z2(k, now)
            a = gl_Gran_normal2(:,i, now)
            di = dabs(DOT_PRODUCT(c,a) - DOT_PRODUCT(a,b))  ! Расстояние от точки до плоскости
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
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS, gl_Gran_neighbour_TVD => dev_gl_Gran_neighbour_TVD, &
		gl_Cell_center2 => dev_gl_Cell_center2, gl_Cell_par2 => dev_gl_Cell_par2
	use GEO_PARAM
	use cudafor
	use My_func
	
	implicit none
	integer, intent(in) :: now
	
	integer i, Num, gr, s1, s2, j, yzel, ss1, ss2
    real(8) :: normal(3), qqq1(8), qqq2(8), POTOK(8), ray(3)
    real(8) :: a1, a2, v1, v2, norm, b1, b2, c1, c2, rad1, rad5
	real(8) :: www, dsl, dsp, dsc, the1, the2, center(3), center2(3)
	integer(4) :: metod
	real(8) :: qqq11(8), qqq22(8), qqq1_TVD(8), qqq2_TVD(8)
	real(8) :: df1, df2, dff1, dff2, rast(3), aa, bb, cc
	
	i = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	Num = size(gl_TS)
	
	
	
	if (i > Num) return
	
	metod = 3 
	www = 0.0
	
		
    gr = gl_TS(i)
    s1 = gl_Gran_neighbour(1, gr)
    s2 = gl_Gran_neighbour(2, gr)
	ss1 = gl_Gran_neighbour_TVD(1, gr)
	ss2 = gl_Gran_neighbour_TVD(2, gr)
    normal = gl_Gran_normal2(:, gr, now)
    qqq1 = gl_Cell_par(1:8, s1)
    qqq2 = gl_Cell_par(1:8, s2)
	
	rast = gl_Gran_center2(:, gr, now) - gl_Cell_center2(:, s1, now)
	df1 = norm2(rast)
	rast = gl_Gran_center2(:, gr, now) - gl_Cell_center2(:, s2, now)
	df2 = norm2(rast)
	rast = gl_Gran_center2(:, gr, now) - gl_Cell_center2(:, ss1, now)
	dff2 = norm2(rast)
	qqq22 = gl_Cell_par(1:8, ss2)
	
	
	rad1 = norm2(gl_Cell_center2(:, s1, now))                                
    rad5 = norm2(gl_Gran_center2(:, gr, now))
	
	qqq1_TVD = qqq1
					
	qqq1_TVD(1) = qqq1_TVD(1) * rad1**2 / rad5**2
    qqq1_TVD(5) = qqq1_TVD(5) * rad1**(2 * ggg) / rad5**(2 * ggg)
					
    call spherical_skorost(gl_Cell_center2(1, s1, now), gl_Cell_center2(2, s1, now), gl_Cell_center2(3, s1, now), &
        qqq1_TVD(2), qqq1_TVD(3), qqq1_TVD(4), aa, bb, cc)
					
    call dekard_skorost(gl_Gran_center2(1, gr, now), gl_Gran_center2(2, gr, now), gl_Gran_center2(3, gr, now), &
        aa, bb, cc, qqq1_TVD(2), qqq1_TVD(3), qqq1_TVD(4))
					
	do i = 1, 8
		qqq2_TVD(i) = qqq2(i)! linear(-dff2, qqq22(i), -df2, qqq2(i), df1, qqq1(i), 0.0_8)
	end do
					
	
	qqq1 = qqq1_TVD
	qqq2 = qqq2_TVD
        
    call chlld(metod, normal(1), normal(2), normal(3), &
            www, qqq1, qqq2, dsl, dsp, dsc, POTOK)
	
		
    dsl = dsl * koef1
	
	if (.True.) then   ! Новый вариант (плохо работает в носике)
	!center = gl_Gran_center2(:, gr, now)
		
	!the1 = polar_angle(center(1), sqrt(center(2)**2 + center(3)**2))
		
	do j = 1, 4
        yzel = gl_all_Gran(j, gr)
			
		!center2(1) = gl_x2(yzel, now)
		!center2(2) = gl_y2(yzel, now)
		!center2(3) = gl_z2(yzel, now)
		!the2 = polar_angle(center2(1), sqrt(center2(2)**2 + center2(3)**2))
		!	
		!if (the2 < the1) CYCLE
			
		! Блок безопасного доступа
			select case (mod(yzel, 4))
			case(0)
				do while ( atomiccas(dev_mutex_1, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(1)
				do while ( atomiccas(dev_mutex_2, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(2)
				do while ( atomiccas(dev_mutex_3, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(3)
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case default
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			end select
			
        gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsl
        gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsl
        gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsl
        gl_Point_num(yzel) = gl_Point_num(yzel) + 1
		
		call threadfence()                               ! Безопасный доступ к памяти
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
        
    return  ! Заканчиваем с этой гранью, переходим к следующей
	end if
	
	
    a1 = sqrt(ggg * qqq1(5)/qqq1(1))  ! Скорости звука
    a2 = sqrt(ggg * qqq2(5)/qqq2(1))
	  
        
    if (norm2(qqq1(2:4)) <= a1 .and. norm2(qqq2(2:4)) <= a2) then ! Записываем скорость во все узлы
        do j = 1, 4
            yzel = gl_all_Gran(j, gr)
			
			! Блок безопасного доступа
			select case (mod(yzel, 4))
			case(0)
				do while ( atomiccas(dev_mutex_1, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(1)
				do while ( atomiccas(dev_mutex_2, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(2)
				do while ( atomiccas(dev_mutex_3, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(3)
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case default
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			end select
			
				gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsl
				gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsl
				gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsl
				gl_Point_num(yzel) = gl_Point_num(yzel) + 1
				
			call threadfence()                               ! Безопасный доступ к памяти
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
        return  ! Заканчиваем с этой гранью, переходим к следующей
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
				do while ( atomiccas(dev_mutex_1, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(1)
				do while ( atomiccas(dev_mutex_2, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(2)
				do while ( atomiccas(dev_mutex_3, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(3)
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case default
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			end select
			
                gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsl
                gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsl
                gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsl
                gl_Point_num(yzel) = gl_Point_num(yzel) + 1
				call threadfence()                               ! Безопасный доступ к памяти
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
				
				
                CYCLE   ! Переходим к следующему узлу
            end if
                
            b1 = DOT_PRODUCT(qqq1(6:8), ray)
            b2 = DOT_PRODUCT(qqq2(6:8), ray)
                
            c1 = dabs(b1)/dsqrt(4.0 * par_pi_8 * qqq1(1))
            c2 = dabs(b2)/dsqrt(4.0 * par_pi_8 * qqq2(1))
                
            b1 = norm2(qqq1(6:8))
            b2 = norm2(qqq2(6:8))
			
                
            ! Если модуль скорости меньше чем максимум из скорости звука и альфвеновской скорости звука по этому направлению
            if( dabs(v1) <= max(a1, 0.5 * (sqrt(b1 * b1/(4.0 * par_pi_8 * qqq1(1)) + a1 * a1+ 2.0 * a1 * c1) &
                + sqrt(b1 * b1/(4.0 * par_pi_8 * qqq1(1)) + a1 * a1 - 2.0 * a1 * c1)  ) ) &
                .and. dabs(v2) <= max(a2, 0.5 * (sqrt(b2 * b2/(4.0 * par_pi_8 * qqq2(1)) + a2 * a2+ 2.0 * a2 * c2) &
                + sqrt(b2 * b2/(4.0 * par_pi_8 * qqq2(1)) + a2 * a2 - 2.0 * a2 * c2)  ) ) ) then
				
				select case (mod(yzel, 4))
				case(0)
					do while ( atomiccas(dev_mutex_1, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
					end do                                        ! Безопасный доступ к памяти
				case(1)
					do while ( atomiccas(dev_mutex_2, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
					end do                                        ! Безопасный доступ к памяти
				case(2)
					do while ( atomiccas(dev_mutex_3, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
					end do                                        ! Безопасный доступ к памяти
				case(3)
					do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
					end do                                        ! Безопасный доступ к памяти
				case default
					do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
					end do                                        ! Безопасный доступ к памяти
				end select
			
                gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsl
                gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsl
                gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsl
                gl_Point_num(yzel) = gl_Point_num(yzel) + 1
				call threadfence()                               ! Безопасный доступ к памяти
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
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS, &
		gl_Gran_neighbour_TVD => dev_gl_Gran_neighbour_TVD, gl_Cell_center2 => dev_gl_Cell_center2, gl_Vn => dev_gl_Vn
	use GEO_PARAM
	use cudafor
	use cgod
	
	implicit none
	integer, intent(in) :: now
	
	integer i, Num, gr, s1, s2, j, yzel
    real(8) :: normal(3), qqq1(8), qqq2(8), POTOK(8), ray(3), vec(3), center(3), center2(3)
    real(8) :: a1, a2, v1, v2, norm, b1, b2, c1, c2
	real(8) :: www, dsl, dsp, dsc
	integer(4) :: metod
	real(8) :: df1, df2, dff1, dff2, rast(3)
	real(8) :: qqq11(8), qqq22(8), qq, qqq1_TVD(8), qqq2_TVD(8), k1
	integer(4) :: ss1, ss2, kdir, KOBL, idgod
	
	KOBL = 0
	kdir = 0
	idgod = 0
	
	i = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
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
	
	k1 = 1.0
	!if (gl_Gran_center2(1, gr, now) >= 2.0) k1 = 0.7 !0.14
	
	! Вычтем нормальную компоненту магнитного поля для значений в самой ячейке
	if (gl_Gran_center2(1, gr, now) >= par_null_bn_x .and. par_null_bn == .True.) then
		qqq1(6:8) = qqq1(6:8) - DOT_PRODUCT(normal, qqq1(6:8)) * normal
		qqq2(6:8) = qqq2(6:8) - DOT_PRODUCT(normal, qqq2(6:8)) * normal
		gl_Cell_par(6:8, s1) = qqq1(6:8)
		gl_Cell_par(6:8, s2) = qqq2(6:8)
		
		
	end if
	
	
	if (.False. .and. par_TVD == .True.) then
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
				
				do i = 1, 8
					qqq1_TVD(i) = linear(-dff1, qqq11(i), -df1, qqq1(i), df2, qqq2(i), 0.0_8)
					qqq2_TVD(i) = linear(-dff2, qqq22(i), -df2, qqq2(i), df1, qqq1(i), 0.0_8)
				end do
				
				
				! Вычтем нормальную компоненту магнитного поля для снесённых значений
			!if ( sqrt(gl_Gran_center2(2, gr, now)**2 + gl_Gran_center2(3, gr, now)**2) <= 15.0 .and. par_null_bn == .True.) then
			if ( gl_Gran_center2(1, gr, now) >= par_null_bn_x .and. par_null_bn == .True.) then
				qqq1_TVD(6:8) = qqq1_TVD(6:8) - DOT_PRODUCT(normal, qqq1_TVD(6:8)) * normal
				qqq2_TVD(6:8) = qqq2_TVD(6:8) - DOT_PRODUCT(normal, qqq2_TVD(6:8)) * normal
				!if (par_TVD == .False.) then
				!	gl_Cell_par(6:8, s1) = qqq1(6:8)
				!	gl_Cell_par(6:8, s2) = qqq2(6:8)
				!end if
			end if
				
				
				qqq1 = qqq1_TVD
				qqq2 = qqq2_TVD
			end if
	end if
	
	
	! Вычтем нормальную компоненту магнитного поля
	!if (gl_Gran_center2(1, gr, now) >= -100.0 .and. par_null_bn == .True.) then
	!	qqq1(6:8) = qqq1(6:8) - DOT_PRODUCT(normal, qqq1(6:8)) * normal
	!	qqq2(6:8) = qqq2(6:8) - DOT_PRODUCT(normal, qqq2(6:8)) * normal
	!	!if (par_TVD == .False.) then
	!	!	gl_Cell_par(6:8, s1) = qqq1(6:8)
	!	!	gl_Cell_par(6:8, s2) = qqq2(6:8)
	!	!end if
	!end if
	
	
	
	
	vec = 0.5 * (qqq1(2:4) + qqq2(2:4))
	vec = vec - DOT_PRODUCT(vec, normal) * normal
	center = gl_Gran_center2(:, gr, now)
        
	dsc = 0.0_8
	
	if (gl_Gran_center2(1, gr, now) >= par_null_bn_x .and. par_null_bn == .True.) then
		! Теперь сделаем для газовой динамики
		qqq1(5) = qqq1(5) + (norm2(qqq1(6:8))**2)/(8.0 * par_pi_8)
		qqq2(5) = qqq2(5) + (norm2(qqq2(6:8))**2)/(8.0 * par_pi_8)
		call cgod3d(KOBL, 0, 0, 0, kdir, idgod, &
                                 normal(1), normal(2), normal(3), 1.0_8, &
                                 www, qqq1(1:8), qqq2(1:8), &
                                 dsl, dsp, dsc, 1.0_8, 1.66666666666666_8, &
                                 POTOK)
		if (idgod == 2) then
		call chlld_gd(3, normal(1), normal(2), normal(3), &
                www, qqq1(1:5), qqq2(1:5), dsl, dsp, dsc, POTOK)
		end if
		
	else
		call chlld(metod, normal(1), normal(2), normal(3), &
				www, qqq1, qqq2, dsl, dsp, dsc, POTOK)
	end if
        
    dsc = k1 * (dsc + 0.0 * DOT_PRODUCT(0.5 * (qqq1(2:4) + qqq2(2:4)), normal)) * koef2
	
	!if (gl_Gran_center2(1, gr, now) < -40.0 .and. gl_Cell_par(9, s1) > 50.0) then
	!		!print*, gl_Gran_center2(:, gr, now) 
	!		dsc = dsc - 20.0
	!end if
	
	!
	!
	!if (gl_Gran_center2(1, gr, now) < -170.0 .and. normal(2) > 0 .and. &
	!		dabs(gl_Gran_center2(3, gr, now)) < 5.0) then
	!		!print*, gl_Gran_center2(:, gr, now) 
	!		dsc = dsc + 10.0
	!end if
	
	
	
	!if(gr == 77) write(*, *) dsc, qqq1(1), qqq2(1), qqq1(5), qqq2(5)
	
	do j = 1, 4
            yzel = gl_all_Gran(j, gr)
			center2(1) = gl_x2(yzel, now)
			center2(2) = gl_y2(yzel, now)
			center2(3) = gl_z2(yzel, now)
			
			if ( sqrt(center2(2)**2 + center2(3)**2) > 5.0) then
			!if ( center2(1) < -40.0) then
			    if( DOT_PRODUCT(vec, center2 - center) < 0.0) CYCLE
			end if
			
			! Блок безопасного доступа
			select case (mod(yzel, 4))
			case(0)
				do while ( atomiccas(dev_mutex_1, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(1)
				do while ( atomiccas(dev_mutex_2, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(2)
				do while ( atomiccas(dev_mutex_3, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(3)
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case default
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			end select
			
            gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsc
            gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsc
            gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsc
			!qqq1(2:4) = qqq1(2:4) + qqq2(2:4)
			!gl_Vn(yzel) = gl_Vn(yzel) + dabs(0.5 * norm2(qqq1(2:4)))
            gl_Point_num(yzel) = gl_Point_num(yzel) + 1
			
			call threadfence()                               ! Безопасный доступ к памяти
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
	return  ! Заканчиваем с этой гранью, переходим к следующей
	
        
    a1 = sqrt(ggg * qqq1(5)/qqq1(1))  ! Скорости звука
    a2 = sqrt(ggg * qqq2(5)/qqq2(1))
	
        
    if (norm2(qqq1(2:4)) <= a1 .and. norm2(qqq2(2:4)) <= a2) then ! Записываем скорость во все узлы
        do j = 1, 4
            yzel = gl_all_Gran(j, gr)
			
			! Блок безопасного доступа
			select case (mod(yzel, 4))
			case(0)
				do while ( atomiccas(dev_mutex_1, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(1)
				do while ( atomiccas(dev_mutex_2, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(2)
				do while ( atomiccas(dev_mutex_3, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(3)
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case default
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			end select
			
            gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsc
            gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsc
            gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsc
            gl_Point_num(yzel) = gl_Point_num(yzel) + 1
			
			call threadfence()                               ! Безопасный доступ к памяти
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
				
				! Блок безопасного доступа
			select case (mod(yzel, 4))
			case(0)
				do while ( atomiccas(dev_mutex_1, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(1)
				do while ( atomiccas(dev_mutex_2, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(2)
				do while ( atomiccas(dev_mutex_3, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(3)
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case default
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			end select
			
                gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsc
                gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsc
                gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsc
                gl_Point_num(yzel) = gl_Point_num(yzel) + 1
				
				call threadfence()                               ! Безопасный доступ к памяти
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
				
                CYCLE   ! Переходим к следующему узлу
            end if
                
            b1 = DOT_PRODUCT(qqq1(6:8), ray)
            b2 = DOT_PRODUCT(qqq2(6:8), ray)
                
            c1 = dabs(b1)/dsqrt(4.0 * par_pi_8 * qqq1(1))
            c2 = dabs(b2)/dsqrt(4.0 * par_pi_8 * qqq2(1))
                
            b1 = norm2(qqq1(6:8))
            b2 = norm2(qqq2(6:8))
                
            ! Если модуль скорости меньше чем максимум из скорости звука и альфвеновской скорости звука по этому направлению
            if( dabs(v1) <= max(a1, 0.5 * (sqrt(b1 * b1/(4.0 * par_pi_8 * qqq1(1)) + a1 * a1+ 2.0 * a1 * c1) &
                + sqrt(b1 * b1/(4.0 * par_pi_8 * qqq1(1)) + a1 * a1 - 2.0 * a1 * c1)  ) ) &
                .and. dabs(v2) <= max(a2, 0.5 * (sqrt(b2 * b2/(4.0 * par_pi_8 * qqq2(1)) + a2 * a2+ 2.0 * a2 * c2) &
                + sqrt(b2 * b2/(4.0 * par_pi_8 * qqq2(1)) + a2 * a2 - 2.0 * a2 * c2)  ) ) ) then
				
				! Блок безопасного доступа
			select case (mod(yzel, 4))
			case(0)
				do while ( atomiccas(dev_mutex_1, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(1)
				do while ( atomiccas(dev_mutex_2, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(2)
				do while ( atomiccas(dev_mutex_3, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(3)
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case default
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			end select
			
                gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsc
                gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsc
                gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsc
                gl_Point_num(yzel) = gl_Point_num(yzel) + 1
				
				call threadfence()                               ! Безопасный доступ к памяти
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
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS, &
		gl_Gran_neighbour_TVD => dev_gl_Gran_neighbour_TVD, gl_Cell_center2 => dev_gl_Cell_center2
	use GEO_PARAM
	use cudafor
	
	implicit none
	integer, intent(in) :: now
	
	integer i, Num, gr, s1, s2, j, yzel
    real(8) :: normal(3), qqq1(8), qqq2(8), POTOK(8), ray(3)
    real(8) :: a1, a2, v1, v2, norm, b1, b2, c1, c2
	real(8) :: www, dsl, dsp, dsc
	integer(4) :: metod
	real(8) :: df1, df2, dff1, dff2, rast(3)
	real(8) :: qqq11(8), qqq22(8), qq, qqq1_TVD(8), qqq2_TVD(8)
	integer(4) :: ss1, ss2
	
	i = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
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
	
	if (par_TVD == .True.) then
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
				
				do i = 1, 8
					qqq1_TVD(i) = linear(-dff1, qqq11(i), -df1, qqq1(i), df2, qqq2(i), 0.0_8)
					qqq2_TVD(i) = linear(-dff2, qqq22(i), -df2, qqq2(i), df1, qqq1(i), 0.0_8)
				end do
				
				qqq1 = qqq1_TVD
				qqq2 = qqq2_TVD
			end if
	end if
        
    call chlld(metod, normal(1), normal(2), normal(3), &
            www, qqq1, qqq2, dsl, dsp, dsc, POTOK)
	
	
        
    dsp = dsp * koef3
	
	
        
    a1 = sqrt(ggg * qqq1(5)/qqq1(1))  ! Скорости звука
    a2 = sqrt(ggg * qqq2(5)/qqq2(1))
        
    if (norm2(qqq1(2:4)) <= a1 .and. norm2(qqq2(2:4)) <= a2) then ! Записываем скорость во все узлы
        do j = 1, 4
            yzel = gl_all_Gran(j, gr)
			
			! Блок безопасного доступа
			select case (mod(yzel, 4))
			case(0)
				do while ( atomiccas(dev_mutex_1, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(1)
				do while ( atomiccas(dev_mutex_2, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(2)
				do while ( atomiccas(dev_mutex_3, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(3)
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case default
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			end select
			
            gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsp
            gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsp
            gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsp
            gl_Point_num(yzel) = gl_Point_num(yzel) + 1
			
			call threadfence()                               ! Безопасный доступ к памяти
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
				
				! Блок безопасного доступа
			select case (mod(yzel, 4))
			case(0)
				do while ( atomiccas(dev_mutex_1, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(1)
				do while ( atomiccas(dev_mutex_2, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(2)
				do while ( atomiccas(dev_mutex_3, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(3)
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case default
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			end select
			
                gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsp
                gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsp
                gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsp
                gl_Point_num(yzel) = gl_Point_num(yzel) + 1
				
				call threadfence()                               ! Безопасный доступ к памяти
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
				
                CYCLE   ! Переходим к следующему узлу
            end if
                
            b1 = DOT_PRODUCT(qqq1(6:8), ray)
            b2 = DOT_PRODUCT(qqq2(6:8), ray)
                
            c1 = dabs(b1)/dsqrt(4.0 * par_pi_8 * qqq1(1))
            c2 = dabs(b2)/dsqrt(4.0 * par_pi_8 * qqq2(1))
                
            b1 = norm2(qqq1(6:8))
            b2 = norm2(qqq2(6:8))
                
            ! Если модуль скорости меньше чем максимум из скорости звука и альфвеновской скорости звука по этому направлению
            if( dabs(v1) <= max(a1, 0.5 * (sqrt(b1 * b1/(4.0 * par_pi_8 * qqq1(1)) + a1 * a1+ 2.0 * a1 * c1) &
                + sqrt(b1 * b1/(4.0 * par_pi_8 * qqq1(1)) + a1 * a1 - 2.0 * a1 * c1)  ) ) &
                .and. dabs(v2) <= max(a2, 0.5 * (sqrt(b2 * b2/(4.0 * par_pi_8 * qqq2(1)) + a2 * a2+ 2.0 * a2 * c2) &
                + sqrt(b2 * b2/(4.0 * par_pi_8 * qqq2(1)) + a2 * a2 - 2.0 * a2 * c2)  ) ) ) then
				
				! Блок безопасного доступа
			select case (mod(yzel, 4))
			case(0)
				do while ( atomiccas(dev_mutex_1, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(1)
				do while ( atomiccas(dev_mutex_2, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(2)
				do while ( atomiccas(dev_mutex_3, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case(3)
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			case default
				do while ( atomiccas(dev_mutex_4, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                        ! Безопасный доступ к памяти
			end select
			
                gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsp
                gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsp
                gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsp
                gl_Point_num(yzel) = gl_Point_num(yzel) + 1
				
				call threadfence()                               ! Безопасный доступ к памяти
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
	! Подвижная сетка, МГД + мультифлюид
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
		gl_Gran_scheme => dev_gl_Gran_scheme, gl_zone_Cell => dev_gl_zone_Cell, gl_Gran_neighbour_TVD => dev_gl_Gran_neighbour_TVD
	use GEO_PARAM
	implicit none
	integer, intent(in) :: now
	
	integer(4) :: gr  ! Глобальный номер текущей грани
	
    integer(4) :: st, s1, s2, i, j, k, zone, metod, now2, ss1, ss2
    real(8) :: qqq1(9), qqq2(9), qqq(9)  ! Переменные в ячейке
    real(8) :: fluid1(5, 4), fluid2(5, 4)
    real(8) :: dist, dsl, dsc, dsp, distant(3)
    real(8) :: POTOK(9), ttest(3)
    real(8) :: POTOK_MF(5)
    real(8) :: POTOK_MF_all(5, 4)
    real(8) :: time, Volume, TT, U8, rad1, rad2, aa, bb, cc, wc
    real(8) :: SOURSE(5,5)  ! Источники массы, импульса и энергии для плазмы и каждого сорта мультифлюида
    real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3
	real(8) :: df1, df2, dff1, dff2, rast(3)
	real(8) :: qqq11(9), qqq22(9), qq, qqq1_TVD(9), qqq2_TVD(9)
	
	logical :: null_bn
	
	now2 = mod(now, 2) + 1
	
	time = 100000.0
	gr = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	
	!if(gr == 1) then
	!	print*, "Hellow from potok (1, 1) ", dev_Ngran
	!end if
	
	TT = time_step
	
	
	! ----------------------------------------------------------- можно просто скопировать код из хоста ----------------------
	if (gr > dev_Ngran) return
	
	metod = 2
	if(gl_Gran_info(gr) == 2) return
	
	POTOK = 0.0
	s1 = gl_Gran_neighbour(1, gr)
	s2 = gl_Gran_neighbour(2, gr)
	qqq1 = gl_Cell_par(1:9, s1)
	fluid1 = gl_Cell_par_MF(:, :, s1)   ! Загрузили параметры жидкостей для мультифлюида
	null_bn = .False.
	
	distant = gl_Gran_center2(:, gr, now) - gl_Cell_center2(:, s1, now)
	dist = norm2(distant)
	
	! Попробуем снести плотность пропорционально квадрату
            if(norm2(qqq1(2:4))/sqrt(ggg*qqq1(5)/qqq1(1)) > 2.2) then
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

                call spherical_skorost(gl_Cell_center2(1, s1, now), gl_Cell_center2(2, s1, now), gl_Cell_center2(3, s1, now), & 
                    fluid1(2, 1), fluid1(3, 1), fluid1(4, 1), aa, bb, cc)
                call dekard_skorost(gl_Gran_center2(1, gr, now), gl_Gran_center2(2, gr, now), gl_Gran_center2(3, gr, now), & 
                    aa, bb, cc, fluid1(2, 1), fluid1(3, 1), fluid1(4, 1))
            end if

            if (s2 >= 1) then
                !if ( norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int .and. norm2(gl_Cell_center(:, s2)) <= par_R0 * par_R_int) CYCLE
                qqq2 = gl_Cell_par(1:9, s2)
                fluid2 = gl_Cell_par_MF(:, :, s2)   ! Загрузили параметры жидкостей для мультифлюида
                !dist = min(gl_Cell_dist(s1), gl_Cell_dist(s2))   
				
				distant = gl_Gran_center2(:, gr, now) - gl_Cell_center2(:, s2, now)
				dist = min( dist, norm2(distant))

                ! Попробуем снести плотность пропорционально квадрату
                if(norm2(qqq2(2:4))/sqrt(ggg*qqq2(5)/qqq2(1)) > 2.2) then 
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
					
					!if(qqq2(6) > 0.0) then
     !                   qqq2(6) = -0.1 ! Отсос магнитного поля
     !               end if

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
			
			if(gl_Gran_type(gr) == 2 .or. gl_Gran_type(gr) == 1) metod = 3 !2
			
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
	
			
			
            call chlld_gd(2, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
                wc, fluid1(:, 1), fluid2(:, 1), dsl, dsp, dsc, POTOK_MF)
            time = min(time, 0.99 * dist/(max(dabs(dsl), dabs(dsp)) + dabs(wc)) )
            gl_Gran_POTOK_MF(:, 1, gr) = POTOK_MF * gl_Gran_square2(gr, now)
            
            call chlld_gd(2, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
                wc, fluid1(:, 2), fluid2(:, 2), dsl, dsp, dsc, POTOK_MF)
            time = min(time, 0.99 * dist/(max(dabs(dsl), dabs(dsp)) + dabs(wc)) )
            gl_Gran_POTOK_MF(:, 2, gr) = POTOK_MF * gl_Gran_square2(gr, now)
            
            call chlld_gd(2, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
                wc, fluid1(:, 3), fluid2(:, 3), dsl, dsp, dsc, POTOK_MF)
            time = min(time, 0.99 * dist/(max(dabs(dsl), dabs(dsp)) + dabs(wc)) )
            gl_Gran_POTOK_MF(:, 3, gr) = POTOK_MF * gl_Gran_square2(gr, now)
            
            call chlld_gd(2, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
                wc, fluid1(:, 4), fluid2(:, 4), dsl, dsp, dsc, POTOK_MF)
            time = min(time, 0.99 * dist/(max(dabs(dsl), dabs(dsp)) + dabs(wc)) )
            gl_Gran_POTOK_MF(:, 4, gr) = POTOK_MF * gl_Gran_square2(gr, now)
	
	! ----------------------------------------------------------- копируем до этого момента ----------------------
	
	!if (time_step2 > time) then 
		! Не смотри шаг по времени в хвосте из-за плохих ячеек
		if (.True.) then !(gl_Cell_center2(1, s1, now) > -150.0) then
			time =  atomicmin(time_step2, time)   ! Атомарная операция взятия минимального значения
		end if
	!end if
	
	
	end subroutine  CUF_MGD_grans_MF
	
	
	attributes(global) subroutine CUF_MGD_cells_MF(now)
	! Подвижная сетка, МГД + мультифлюид
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
		Calc_sourse_MF => dev_Calc_sourse_MF, par_kk1 => dev_par_kk1, gl_Cell_type => dev_gl_Cell_type, gl_Cell_number => dev_gl_Cell_number, &
		gl_zone_Cell => dev_gl_zone_Cell
	use GEO_PARAM
	implicit none
	integer, intent(in) :: now
	
	integer(4) :: gr
	
    integer(4) :: st, s1, s2, i, j, k, zone, now2, ijk
    real(8) :: qqq1(9), qqq2(9), qqq(9)  ! Переменные в ячейке
    real(8) :: fluid1(5, 4), fluid2(5, 4)
    real(8) :: dist, dsl, dsc, dsp
    real(8) :: POTOK(9), ttest(3)
    real(8) :: POTOK_MF(5)
    real(8) :: POTOK_MF_all(5, 4)
    real(8) :: time, Volume, U8, rad1, rad2, aa, bb, cc, Volume2, sks
    real(8) :: SOURSE(5,5)  ! Источники массы, импульса и энергии для плазмы и каждого сорта мультифлюида
    real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3
	logical :: l_1
	
	now2 = mod(now, 2) + 1
	gr = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	time = time_step
	
	if (gr > dev_Ncell) return
	
	! ----------------------------------------------------------- можно просто скопировать код из хоста ----------------------
	
	if(gl_Cell_info(gr) == 0) return
            l_1 = .TRUE.
            if (norm2(gl_Cell_center2(:, gr, now)) <= par_R0 + (par_R_character - par_R0) * (3.0_8/par_n_TS)**par_kk1) l_1 = .FALSE. 
			!if ((gl_Cell_type(gr) == "A" .or. gl_Cell_type(gr) == "B").and.(gl_Cell_number(1, gr) <= 2) ) l_1 = .FALSE.    ! Не считаем внутри сферы
            POTOK = 0.0
			sks = 0.0
            SOURSE = 0.0
            POTOK_MF_all = 0.0
            Volume = gl_Cell_Volume2(gr, now)
            Volume2 = gl_Cell_Volume2(gr, now2)
			
			
	
			!Volume2 = Volume
	
            qqq = gl_Cell_par(:, gr)
            fluid1 = gl_Cell_par_MF(:, :, gr)
            ! Просуммируем потоки через грани
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


            ! Определяем зону в которой находимся
            if(qqq(9)/qqq(1) < 50.0) then
                if(norm2(qqq(2:4))/sqrt(ggg*qqq(5)/qqq(1)) > 1.3) then
				!if(gl_zone_Cell(gr) == 1) then
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

            call Calc_sourse_MF(qqq, fluid1, SOURSE, zone)  ! Вычисляем источники
		


            if (l_1 == .TRUE.) then
                ro3 = qqq(1)* Volume / Volume2 - time * POTOK(1) / Volume2
                Q3 = qqq(9)* Volume / Volume2 - time * POTOK(9) / Volume2
                if (ro3 <= 0.0_8) then
                    write(*, *) "Ro < 0  1490 ", ro3, gl_Cell_center2(1, gr, now), gl_Cell_center2(2, gr, now), gl_Cell_center2(3, gr, now)
					write(*, *) qqq(1), Q3
					!write(*, *) Volume , Volume2
					ro3 = 0.15
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

            ! Теперь посчитаем законы сохранения для остальных жидкостей
			
			
            do i = 1, 4
                if (i == 1 .and. l_1 == .FALSE.) CYCLE       ! Пропускаем внутреннюю сферу для сорта 1
                if (l_1 == .FALSE.) SOURSE(:, i + 1) = 0.0       ! Пропускаем источники на внутреннюю сферу для всех сортов
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
	! Подвижная сетка, МГД + мультифлюид
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
	
	integer(4) :: gr  ! Глобальный номер текущей грани
	
    integer(4) :: st, s1, s2, i, j, k, zone, metod, now2, iter
    real(8) :: qqq1(9), qqq2(9), qqq(9)  ! Переменные в ячейке
    real(8) :: fluid1(5, 4), fluid2(5, 4)
    real(8) :: dist, dsl, dsc, dsp, distant(3)
    real(8) :: POTOK(9), ttest(3)
    real(8) :: POTOK_MF(5)
    real(8) :: POTOK_MF_all(5, 4)
    real(8) :: time, Volume, TT, U8, rad1, rad2, aa, bb, cc, wc
    real(8) :: SOURSE(5,5)  ! Источники массы, импульса и энергии для плазмы и каждого сорта мультифлюида
    real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3
	
	logical :: null_bn
	
	
	time = 100000.0
	iter = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	
	! ----------------------------------------------------------- можно просто скопировать код из хоста ----------------------
	if (iter > size(gl_all_Gran_inner)) return
	
	gr = gl_all_Gran_inner(iter)
        !if(gl_Gran_info(gr) == 0) CYCLE
        POTOK = 0.0
        s1 = gl_Gran_neighbour(1, gr)
        s2 = gl_Gran_neighbour(2, gr)
        qqq1 = gl_Cell_par(1:9, s1)
        fluid1 = gl_Cell_par_MF(:, :, s1)   ! Загрузили параметры жидкостей для мультифлюида
		
		distant = gl_Gran_center(:, gr) - gl_Cell_center(:, s1)
		dist = norm2(distant)

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
			
			distant = gl_Gran_center(:, gr) - gl_Cell_center(:, s2)
			dist = min( dist, norm2(distant))

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
			else
				write(*,*), "YTGVBNYTGBJYGBNJUYHBNJUHGBNMIUHKIJ"
			end if
	end if



        call chlld_Q(3, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
            0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)
        time = min(time, 0.99 * dist/max(dabs(dsl), dabs(dsp)) )   ! REDUCTION
        gl_Gran_POTOK(1:9, gr) = POTOK * gl_Gran_square(gr)
		gl_Gran_POTOK(10, gr) = 0.5 * DOT_PRODUCT(gl_Gran_normal(:, gr), qqq1(6:8) + qqq2(6:8)) * gl_Gran_square(gr)


        call chlld_gd(2, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
            0.0_8, fluid1(:, 1), fluid2(:, 1), dsl, dsp, dsc, POTOK_MF)
        time = min(time, 0.99 * dist/max(dabs(dsl), dabs(dsp)) )
        gl_Gran_POTOK_MF(:, 1, gr) = POTOK_MF * gl_Gran_square(gr)
        
        call chlld_gd(2, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
            0.0_8, fluid1(:, 2), fluid2(:, 2), dsl, dsp, dsc, POTOK_MF)
        time = min(time, 0.99 * dist/max(dabs(dsl), dabs(dsp)) )
        gl_Gran_POTOK_MF(:, 2, gr) = POTOK_MF * gl_Gran_square(gr)
        
        
        call chlld_gd(2, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
            0.0_8, fluid1(:, 3), fluid2(:, 3), dsl, dsp, dsc, POTOK_MF)
        time = min(time, 0.99 * dist/max(dabs(dsl), dabs(dsp)) )
        gl_Gran_POTOK_MF(:, 3, gr) = POTOK_MF * gl_Gran_square(gr)
        
        call chlld_gd(2, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
            0.0_8, fluid1(:, 4), fluid2(:, 4), dsl, dsp, dsc, POTOK_MF)
        time = min(time, 0.99 * dist/max(dabs(dsl), dabs(dsp)) )
        gl_Gran_POTOK_MF(:, 4, gr) = POTOK_MF * gl_Gran_square(gr)

	
	! ----------------------------------------------------------- копируем до этого момента ----------------------
	
	if (dev_time_step_inner > time) then 
		time =  atomicmin(dev_time_step_inner, time)   ! Атомарная операция взятия минимального значения
	end if
	
	
	end subroutine  CUF_MGD_grans_MF_inner
	
	
	attributes(global) subroutine CUF_MGD_cells_MF_inner()
	! Подвижная сетка, МГД + мультифлюид
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
    real(8) :: qqq1(9), qqq2(9), qqq(9)  ! Переменные в ячейке
    real(8) :: fluid1(5, 4), fluid2(5, 4)
    real(8) :: dist, dsl, dsc, dsp
    real(8) :: POTOK(9), ttest(3)
    real(8) :: POTOK_MF(5)
    real(8) :: POTOK_MF_all(5, 4)
    real(8) :: Volume, TT, U8, rad1, rad2, aa, bb, cc, Volume2, sks
    real(8) :: SOURSE(5,5)  ! Источники массы, импульса и энергии для плазмы и каждого сорта мультифлюида
    real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3
	logical :: l_1
	
	iter = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	
	if (iter > size(gl_all_Cell_inner)) return
	
	! ----------------------------------------------------------- можно просто скопировать код из хоста ----------------------
	
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
				
				!if(gr == 5) then
				!write(*,*) ro3, u3, v3, w3
				!write(*,*) POTOK(1), Volume, qqq(1), dev_time_step_inner
				!write(*,*) "____________"
				!	
				!end if

            end if

            ! Теперь посчитаем законы сохранения для остальных жидкостей

            do i = 1, 4
                if (i == 1 .and. l_1 == .FALSE.) CYCLE
				
				if (i == 2 .and. l_1 == .FALSE.) CYCLE
                
                if (l_1 == .FALSE.) SOURSE(:, i + 1) = 0.0       ! Не перезаряжаем жидкости внутри мини сферы
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
	
	attributes(global) subroutine CUF_MGD_grans_MK(now)
	! Подвижная сетка, МГД + мультифлюид
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
		gl_Gran_scheme => dev_gl_Gran_scheme, gl_zone_Cell => dev_gl_zone_Cell, gl_Gran_neighbour_TVD => dev_gl_Gran_neighbour_TVD, &
		gl_Cell_par2 => dev_gl_Cell_par2, gl_Gran_POTOK2 => dev_gl_Gran_POTOK2
	use GEO_PARAM
	implicit none
	integer, intent(in) :: now
	
	integer(4) :: gr  ! Глобальный номер текущей грани
	
    integer(4) :: st, s1, s2, i, j, k, zone, metod, now2, ss1, ss2
    real(8) :: qqq1(9), qqq2(9), qqq(9)  ! Переменные в ячейке
    real(8) :: dist, dsl, dsc, dsp, distant(3), konvect_(3, 1)
    real(8) :: POTOK(9), ttest(3)
    real(8) :: time, Volume, TT, U8, rad1, rad2, aa, bb, cc, wc, rad3, rad4, rad5
    real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3
	real(8) :: df1, df2, dff1, dff2, rast(3)
	real(8) :: qqq11(9), qqq22(9), qqq1_TVD(9), qqq2_TVD(9), qqq11_2(1), qqq22_2(1), qqq1_TVD_2(1), qqq2_TVD_2(1)
	logical :: tvd1, tvd2, tvd3, tvd4  ! Нужно ли делать особый снос в гиперзвуковом источнике
	
	logical :: null_bn
	
	now2 = mod(now, 2) + 1

	konvect_(3, 1) = 0.0
	time = 100000.0
	gr = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	
	!if(gr == 1) then
	!	print*, "Hellow from potok (1, 1) ", dev_Ngran
	!end if
	
	TT = time_step
	
	
	! ----------------------------------------------------------- можно просто скопировать код из хоста ----------------------
	if (gr > dev_Ngran) return
	
	metod = 2
	if(gl_Gran_info(gr) == 2) return
	
	POTOK = 0.0
	s1 = gl_Gran_neighbour(1, gr)
	s2 = gl_Gran_neighbour(2, gr)
	qqq1 = gl_Cell_par(:, s1)
	konvect_(1, 1) = gl_Cell_par2(1, s1)
	null_bn = .False.
	
	distant = gl_Gran_center2(:, gr, now) - gl_Cell_center2(:, s1, now)
	dist = norm2(distant)
	
	tvd1 = (gl_zone_Cell(s1) == 1) !(norm2(qqq1(2:4))/sqrt(ggg*qqq1(5)/qqq1(1)) > 2.2)
	
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
						konvect_(2, 1) = qqq2(1) * (1.0 - 1.0/1.15)
		
						qqq2(1) = qqq2(1) * 1.15
						qqq2(9) = qqq2(9) * 1.15
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
						konvect_(2, 1) = qqq2(1)  * (1.0 - 1.0/1.15)
		
						qqq2(1) = qqq2(1) * 1.15
						qqq2(9) = qqq2(9) * 1.15
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
	
	tvd2 = (gl_zone_Cell(s2) == 1) !(norm2(qqq2(2:4))/sqrt(ggg*qqq2(5)/qqq2(1)) > 2.2)
	
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
				rast = gl_Gran_center2(:, gr, now) - gl_Cell_center2(:, ss1, now)
				dff2 = norm2(rast)
				qqq11 = gl_Cell_par(:, ss1)
				qqq22 = gl_Cell_par(:, ss2)
				qqq11_2 = gl_Cell_par2(:, ss1)
				qqq22_2 = gl_Cell_par2(:, ss2)
				
				tvd3 = (gl_zone_Cell(ss1) == 1) !(norm2(qqq11(2:4))/sqrt(ggg*qqq11(5)/qqq11(1)) > 2.2)
				tvd4 = (gl_zone_Cell(ss2) == 1) !(norm2(qqq22(2:4))/sqrt(ggg*qqq22(5)/qqq22(1)) > 2.2)
				
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
						qqq2_TVD(i) = qqq2(i) ! linear(-dff2, qqq22(i), -df2, qqq2(i), df1, qqq1(i), 0.0_8)
					end do
					
					qqq2_TVD_2(1) = konvect_(2, 1)! linear(-dff2, qqq22_2(1), -df2, konvect_(2, 1), df1, konvect_(1, 1), 0.0_8)	
					
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
	
	! Вычитаем для снесённых значений нормальною компоненту магнитного поля
	!if (gl_Gran_type(gr) == 2 .and. sqrt(gl_Gran_center2(2, gr, now)**2 + gl_Gran_center2(3, gr, now)**2) <= 15.0 .and. par_null_bn == .True.) then
	if (gl_Gran_type(gr) == 2 .and. gl_Gran_center2(1, gr, now) >= par_null_bn_x .and. par_null_bn == .True.) then
		qqq1(6:8) = qqq1(6:8) - DOT_PRODUCT(gl_Gran_normal2(:, gr, now), qqq1(6:8)) * gl_Gran_normal2(:, gr, now)
		qqq2(6:8) = qqq2(6:8) - DOT_PRODUCT(gl_Gran_normal2(:, gr, now), qqq2(6:8)) * gl_Gran_normal2(:, gr, now)
	end if
	
            
            ! Нужно вычислить скорость движения грани
            wc = DOT_PRODUCT((gl_Gran_center2(:, gr, now2) -  gl_Gran_center2(:, gr, now))/TT, gl_Gran_normal2(:, gr, now))
			
	
			
			metod = gl_Gran_scheme(gr)
			
			
			if(gl_Gran_type(gr) == 2 .or. gl_Gran_type(gr) == 1) metod = 3 !2
			
            if (.False.) then !(gl_Gran_type(gr) == 1) then
				call chlld_Q(metod, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
                wc, qqq1, qqq2, dsl, dsp, dsc, POTOK, null_bn, 0, p_correct_ = .True., konvect_ = konvect_)
			else
				call chlld_Q(metod, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
                wc, qqq1, qqq2, dsl, dsp, dsc, POTOK, null_bn, p_correct_ = .True., konvect_ = konvect_)
			end if
	
	
	
            time = min(time, 0.99 * dist/(max(dabs(dsl), dabs(dsp))+ dabs(wc)) )   ! REDUCTION
            gl_Gran_POTOK(1:9, gr) = POTOK * gl_Gran_square2(gr, now)
            gl_Gran_POTOK2(1, gr) = konvect_(3, 1) * gl_Gran_square2(gr, now)
			
			gl_Gran_POTOK(10, gr) = 0.5 * DOT_PRODUCT(gl_Gran_normal2(:, gr, now), qqq1(6:8) + qqq2(6:8)) * gl_Gran_square2(gr, now)
			
	
	
	! ----------------------------------------------------------- копируем до этого момента ----------------------
	
	!if (time_step2 > time) then 
		! Не смотри шаг по времени в хвосте из-за плохих ячеек
		if (.True.) then !(gl_Cell_center2(1, s1, now) > -150.0) then
			time =  atomicmin(time_step2, time)   ! Атомарная операция взятия минимального значения
		end if
	!end if
	
	
	end subroutine  CUF_MGD_grans_MK
	
	
	
	attributes(global) subroutine CUF_MGD_cells_MK(now)
	! Подвижная сетка, МГД + источники Монте-Карло
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
		Calc_sourse_MF => dev_Calc_sourse_MF, par_kk1 => dev_par_kk1, gl_Cell_type => dev_gl_Cell_type, gl_Cell_number => dev_gl_Cell_number, &
		gl_zone_Cell => dev_gl_zone_Cell, gl_Cell_par_MK => dev_gl_Cell_par_MK, gl_Gran_POTOK2 => dev_gl_Gran_POTOK2, &
		gl_Cell_par2 => dev_gl_Cell_par2
	use GEO_PARAM
	implicit none
	integer, intent(in) :: now
	
	integer(4) :: gr
	
    integer(4) :: st, s1, s2, i, j, k, zone, now2, ijk
    real(8) :: qqq(9), qqq2  ! Переменные в ячейке
    real(8) :: fluid1(5, 4), MK_kk(5)
    real(8) :: dist, dsl, dsc, dsp
    real(8) :: POTOK(9), ttest(3), POTOK2
    real(8) :: time, Volume, U8, rad1, rad2, aa, bb, cc, Volume2, sks
    real(8) :: SOURSE(5,5)  ! Источники массы, импульса и энергии для плазмы и каждого сорта мультифлюида
    real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3
	logical :: l_1
	
	now2 = mod(now, 2) + 1
	gr = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	time = time_step
	
	if (gr > dev_Ncell) return
	
	! ----------------------------------------------------------- можно просто скопировать код из хоста ----------------------
	
	if(gl_Cell_info(gr) == 0) return
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
            qqq2 = gl_Cell_par2(1, gr)
			
            fluid1 = gl_Cell_par_MK(1:5, 1:4, gr)
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
			qqq(1) = qqq(1) - qqq2
			
			if(zone == 1 .or. zone == 2) then
				qqq(5) = qqq(5) / (2.0 * qqq(1) + 1.5 * qqq2) * 2.0 * qqq(1)
			else
				qqq(5) = qqq(5) / (2.0 * qqq(1) + qqq2) * 2.0 * qqq(1)
			end if
			
            call Calc_sourse_MF(qqq, fluid1, SOURSE, zone)  ! Вычисляем источники
			
			
			if(zone == 1 .or. zone == 2) then
				qqq(5) = qqq(5) * (2.0 * qqq(1) + 1.5 * qqq2) / (2.0 * qqq(1))
			else
				qqq(5) = qqq(5) * (2.0 * qqq(1) + qqq2) / (2.0 * qqq(1))
			end if
			qqq(1) = qqq(1) + qqq2
			
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
				gl_Cell_par2(1, gr) = qqq2 * Volume / Volume2 - time * POTOK2 / Volume2
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
     !               write(*, *) gl_Cell_par2(1, gr), qqq2, POTOK2, gl_Cell_center2(1, gr, now), gl_Cell_center2(2, gr, now), gl_Cell_center2(3, gr, now)
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
 
	end subroutine CUF_MGD_cells_MK
	
	
	attributes(global) subroutine CUF_MGD_grans_MK_inner()
	! Подвижная сетка, МГД + мультифлюид
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
		gl_Cell_number => dev_gl_Cell_number, gl_Cell_Volume => dev_gl_Cell_Volume, gl_Gran_neighbour_TVD => dev_gl_Gran_neighbour_TVD, &
		gl_Cell_par2 => dev_gl_Cell_par2, gl_Gran_POTOK2 => dev_gl_Gran_POTOK2
	use GEO_PARAM
	implicit none
	
	integer(4) :: gr  ! Глобальный номер текущей грани
	
    integer(4) :: st, s1, s2, i, j, k, zone, metod, now2, iter, ss1, ss2
    real(8) :: qqq1(9), qqq2(9), qqq(9)  ! Переменные в ячейке
    real(8) :: dist, dsl, dsc, dsp, distant(3), konvect_(3, 1)
    real(8) :: POTOK(9), ttest(3)
    real(8) :: POTOK_MF(5)
    real(8) :: time, Volume, TT, U8, rad1, rad2, aa, bb, cc, wc, rad3, rad4, rad5
	real(8) :: df1, df2, dff1, dff2, rast(3)
    real(8) :: SOURSE(5,5)  ! Источники массы, импульса и энергии для плазмы и каждого сорта мультифлюида
    real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3
	real(8) :: qqq11(9), qqq22(9), qqq1_TVD(9), qqq2_TVD(9), qqq11_2(1), qqq22_2(1), qqq1_TVD_2(1), qqq2_TVD_2(1)
	
	logical :: null_bn
	
	konvect_(3, 1) = 0.0
	time = 100000.0
	iter = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	
	! ----------------------------------------------------------- можно просто скопировать код из хоста ----------------------
	if (iter > size(gl_all_Gran_inner)) return
	
	gr = gl_all_Gran_inner(iter)
        !if(gl_Gran_info(gr) == 0) CYCLE
        POTOK = 0.0
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
  !
  !      ! Попробуем снести плотность пропорционально квадрату
  !      rad1 = norm2(gl_Cell_center(:, s2))
  !      rad2 = norm2(gl_Gran_center(:, gr))
  !      qqq2(1) = qqq2(1) * rad1**2 / rad2**2
  !      qqq2(9) = qqq2(9) * rad1**2 / rad2**2
  !      qqq2(5) = qqq2(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
  !      call spherical_skorost(gl_Cell_center(3, s2), gl_Cell_center(1, s2), gl_Cell_center(2, s2), &
  !          qqq2(4), qqq2(2), qqq2(3), aa, bb, cc)
  !      call dekard_skorost(gl_Gran_center(3, gr), gl_Gran_center(1, gr), gl_Gran_center(2, gr), &
  !          aa, bb, cc, qqq2(4), qqq2(2), qqq2(3))
		
		ss1 = gl_Gran_neighbour_TVD(1, gr)
		ss2 = gl_Gran_neighbour_TVD(2, gr)
		if (ss1 /= 0 .and. ss2 /= 0) then
			rast = gl_Gran_center(:, gr) - gl_Cell_center(:, s1)
			df1 = norm2(rast)
			rast = gl_Gran_center(:, gr) - gl_Cell_center(:, s2)
			df2 = norm2(rast)
			rast = gl_Gran_center(:, gr) - gl_Cell_center(:, ss1)
			dff1 = norm2(rast)
			rast = gl_Gran_center(:, gr) - gl_Cell_center(:, ss1)
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
					
			call spherical_skorost(gl_Cell_center(3, s2), gl_Cell_center(1, s2), gl_Cell_center(2, s2), &
                qqq2(4), qqq2(2), qqq2(3), aa, bb, cc)
			qqq2(4) = aa
			qqq2(2) = bb
			qqq2(3) = cc
					
			call spherical_skorost(gl_Cell_center(3, ss1), gl_Cell_center(1, ss1), gl_Cell_center(2, ss1), &
                qqq11(4), qqq11(2), qqq11(3), aa, bb, cc)
			qqq11(4) = aa
			qqq11(2) = bb
			qqq11(3) = cc
					
			call spherical_skorost(gl_Cell_center(3, ss2), gl_Cell_center(1, ss2), gl_Cell_center(2, ss2), &
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
 
 
	
	! ----------------------------------------------------------- копируем до этого момента ----------------------
	
	if (dev_time_step_inner > time) then 
		time =  atomicmin(dev_time_step_inner, time)   ! Атомарная операция взятия минимального значения
	end if
	
	
	end subroutine  CUF_MGD_grans_MK_inner
	
	
	attributes(global) subroutine CUF_MGD_cells_MK_inner()
	! Подвижная сетка, МГД + мультифлюид
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
		gl_Cell_number => dev_gl_Cell_number, gl_Cell_Volume => dev_gl_Cell_Volume, gl_Cell_par_MK => dev_gl_Cell_par_MK, gl_Gran_POTOK2 => dev_gl_Gran_POTOK2, &
		gl_Cell_par2 => dev_gl_Cell_par2
	use GEO_PARAM
	implicit none
	
	integer(4) :: gr
	
    integer(4) :: st, s1, s2, i, j, k, zone, iter
    real(8) :: qqq(9), qqq2  ! Переменные в ячейке
    real(8) :: fluid1(5, 4), MK_kk(5)
    real(8) :: dist, dsl, dsc, dsp
    real(8) :: POTOK(9), ttest(3), POTOK2
    real(8) :: POTOK_MF(5)
    real(8) :: Volume, TT, U8, rad1, rad2, aa, bb, Volume2, sks
    real(8) :: SOURSE(5,5)  ! Источники массы, импульса и энергии для плазмы и каждого сорта мультифлюида
    real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3
	logical :: l_1
	
	iter = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	
	if (iter > size(gl_all_Cell_inner)) return
	
	! ----------------------------------------------------------- можно просто скопировать код из хоста ----------------------
	
	gr = gl_all_Cell_inner(iter)
            !if(gl_Cell_info(gr) == 2) CYCLE
            l_1 = .TRUE.
            if ((gl_Cell_type(gr) == "A" .or. gl_Cell_type(gr) == "B").and.(gl_Cell_number(1, gr) <= 2) ) l_1 = .FALSE.    ! Не считаем в первых двух ячейках
            POTOK = 0.0
            POTOK2 = 0.0
			sks = 0.0
            SOURSE = 0.0
            Volume = gl_Cell_Volume(gr)
            qqq = gl_Cell_par(:, gr)
			qqq2 = gl_Cell_par2(1, gr)
            fluid1 = gl_Cell_par_MK(1:5, 1:4, gr)
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
			qqq(1) = qqq(1) - qqq2
			qqq(5) = qqq(5) / (2.0 * qqq(1) + 1.5 * qqq2) * 2.0 * qqq(1)
            
			call Calc_sourse_MF(qqq, fluid1, SOURSE, zone)  ! Вычисляем источники
			
			qqq(5) = qqq(5) * (2.0 * qqq(1) + 1.5 * qqq2) / (2.0 * qqq(1))
			qqq(1) = qqq(1) + qqq2
			
			qqq(2:4) = qqq(2:4) / (par_chi/par_chi_real)
			qqq(1) = qqq(1) * (par_chi/par_chi_real)**2
			SOURSE(5, 1) = SOURSE(5, 1)/ (par_chi/par_chi_real)
 
            if (l_1 == .TRUE.) then
                ro3 = qqq(1) - dev_time_step_inner * POTOK(1) / Volume + dev_time_step_inner * MK_kk(1)
				gl_Cell_par2(1, gr) = qqq2 - dev_time_step_inner * POTOK2 / Volume
                Q3 = qqq(9) - dev_time_step_inner * POTOK(9) / Volume + (qqq(9)/qqq(1)) *  dev_time_step_inner * MK_kk(1)
                if (ro3 <= 0.0_8) then
                    write(*, *) "Ro < 0  3242341234 "
					stop
				end if
				
				if (gl_Cell_par2(1, gr) <= 0.0_8) then
                    write(*, *) "gl_Cell_par2(1, gr) < 0  3242341234 ", gl_Cell_par2(1, gr), qqq2
					stop
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
				
				!if(gr == 5) then
				!write(*,*) ro3, u3, v3, w3
				!write(*,*) POTOK(1), Volume, qqq(1), dev_time_step_inner
				!write(*,*) "____________"
				!	
				!end if
 
            end if
	
	end subroutine CUF_MGD_cells_MK_inner
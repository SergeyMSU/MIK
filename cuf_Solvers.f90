!******7**10**********************************************************70
!**                                                                  **
!**    chlld   for                                                   **
!**                                                                  **
!******7**10**********************************************************70
!**    shema HLLD(modification) for 3D MHD                           **
!******7**10**********************************************************70

	subroutine Cuda_Move_all(now)
	use MY_CUDA, gl_TS => dev_gl_TS, gl_Gran_neighbour => dev_gl_Gran_neighbour, gl_Gran_normal2 => dev_gl_Gran_normal2, &
		gl_Cell_par => dev_gl_Cell_par, gl_all_Gran => dev_gl_all_Gran, gl_Point_num => dev_gl_Point_num, gl_x2 => dev_gl_x2, &
		gl_y2 => dev_gl_y2, gl_z2 => dev_gl_z2, gl_Gran_center2 => dev_gl_Gran_center2, gl_Vx => dev_gl_Vx, &
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS, &
		gl_RAY_A => dev_gl_RAY_A, gl_RAY_B => dev_gl_RAY_B, gl_RAY_C => dev_gl_RAY_C, gl_RAY_O => dev_gl_RAY_O, &
		gl_RAY_K => dev_gl_RAY_K, gl_RAY_D => dev_gl_RAY_D, gl_RAY_E => dev_gl_RAY_E, norm2 => dev_norm2
	use GEO_PARAM
	implicit none
    integer, intent(in) :: now
	
	real(8), device :: Time
    real(8), device :: vel(3), Ak(3), Bk(3), Ck(3), Dk(3), Ek(3), ER(3), KORD(3), ER2(3)
	real(8) :: R_TS, proect, R_HP, R_BS
    integer:: i, j, k, yzel2
    real(8), device :: the, phi, r, x, y, z, rr, xx, x2, y2, z2
    integer :: now2             ! Эти параметры мы сейчас меняем на основе now
	
	integer(4):: ierrSync, ierrAsync, a1, a2, a3, N1, N2, N3
	integer(4) :: yzel
	
	print*, "START A "
	
    now2 = mod(now, 2) + 1  
	Time = time_step
	
	print*, "START B "
	! Поверхностное натяжение ------------------------------------------------------------------------------
	
	N3 = size(gl_RAY_A(1, 1, :))
    N2 = size(gl_RAY_A(1, :, 1))
    N1 = size(gl_RAY_A(:, 1, 1))
	
	print*, "START C "
	
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc start: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync start: ', cudaGetErrorString(cudaGetLastError())
	
	print*, "START D "

	a1 = N1
	a2 = N2
	a3 = N3
	print*, "START A1 ", a1, a2, a3, par_n_TS
	
    ! Цикл движения точек на лучах А 
	!$cuf kernel do(2) <<<*, *>>>
    do k = 1, N3
        do j = 1, N2
			
			if (k /= 1 .and. j == 1) then
                    CYCLE
			end if
			
			! Ударная волна
			yzel = yzel + gl_RAY_A(par_n_TS, j, k)
			 Ak(1) = gl_x2(yzel, now); Ak(2) = gl_y2(yzel, now); Ak(3) = gl_z2(yzel, now)
			 ! Ak = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/)
			
			if (j < N2) then
			    yzel2 = gl_RAY_A(par_n_TS, j + 1, k)
			else
				yzel2 = gl_RAY_B(par_n_TS, 1, k)
			end if
			Bk(1) = gl_x2(yzel2, now); Bk(2) = gl_y2(yzel2, now); Bk(3) = gl_z2(yzel2, now)
			! Bk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if (k > 1) then
			    yzel2 = gl_RAY_A(par_n_TS, j, k - 1)
			else
			    yzel2 = gl_RAY_A(par_n_TS, j, N3)
			end if
			Ck(1) = gl_x2(yzel2, now); Ck(2) = gl_y2(yzel2, now); Ck(3) = gl_z2(yzel2, now)
			!Ck = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if (j > 1) then
			    yzel2 = gl_RAY_A(par_n_TS, j - 1, k)
			else
				yzel2 = gl_RAY_A(par_n_TS, j, k + ceiling(1.0 * N3/2))
			end if
			Dk(1) = gl_x2(yzel2, now); Dk(2) = gl_y2(yzel2, now); Dk(3) = gl_z2(yzel2, now)
			! Dk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if(k < N3) then
			    yzel2 = gl_RAY_A(par_n_TS, j, k + 1)
			else
				yzel2 = gl_RAY_A(par_n_TS, j, 1)
			end if
			Ek(1) = gl_x2(yzel2, now); Ek(2) = gl_y2(yzel2, now); Ek(3) = gl_z2(yzel2, now)
			! Ek = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if (gl_Point_num(yzel) > 0) then
			    vel = gl_Point_num(yzel) * par_nat_TS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time
			else
				vel = par_nat_TS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time
			end if
			
			
			gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
			gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
			gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
			
			 ! Контакт
			
			yzel = gl_RAY_A(par_n_HP, j, k)
			Ak(1) = gl_x2(yzel, now); Ak(2) = gl_y2(yzel, now); Ak(3) = gl_z2(yzel, now)
			! Ak = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/)
			
			if (j < N2) then
			    yzel2 = gl_RAY_A(par_n_HP, j + 1, k)
			else
				yzel2 = gl_RAY_B(par_n_HP, 1, k)
			end if
			Bk(1) = gl_x2(yzel2, now); Bk(2) = gl_y2(yzel2, now); Bk(3) = gl_z2(yzel2, now)
			! Bk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if (k > 1) then
			    yzel2 = gl_RAY_A(par_n_HP, j, k - 1)
			else
			    yzel2 = gl_RAY_A(par_n_HP, j, N3)
			end if
			
			Ck(1) = gl_x2(yzel2, now); Ck(2) = gl_y2(yzel2, now); Ck(3) = gl_z2(yzel2, now)
			!Ck = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if (j > 1) then
			    yzel2 = gl_RAY_A(par_n_HP, j - 1, k)
			else
				yzel2 = gl_RAY_A(par_n_HP, j, k + ceiling(1.0 * N3/2))
			end if
			Dk(1) = gl_x2(yzel2, now); Dk(2) = gl_y2(yzel2, now); Dk(3) = gl_z2(yzel2, now)
			! Dk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if(k < N3) then
			    yzel2 = gl_RAY_A(par_n_HP, j, k + 1)
			else
				yzel2 = gl_RAY_A(par_n_HP, j, 1)
			end if
			Ek(1) = gl_x2(yzel2, now); Ek(2) = gl_y2(yzel2, now); Ek(3) = gl_z2(yzel2, now)
			! Ek = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if (gl_Point_num(yzel) > 0) then
			    vel = gl_Point_num(yzel) * par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time
			else
				vel = par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time
			end if
			
			
			gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
			gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
			gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
			
			
			! Внешняя ударная волна
			
			yzel = gl_RAY_A(par_n_BS, j, k)
			Ak(1) = gl_x2(yzel, now); Ak(2) = gl_y2(yzel, now); Ak(3) = gl_z2(yzel, now)
			! Ak = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/)
			
			if (j < N2) then
			    yzel2 = gl_RAY_A(par_n_BS, j + 1, k)
			else
				CYCLE !yzel2 = yzel                                                    ! Вот тут есть неопределённость что делать с крайним узлом
			end if
			Bk(1) = gl_x2(yzel2, now); Bk(2) = gl_y2(yzel2, now); Bk(3) = gl_z2(yzel2, now)
			! Bk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if (k > 1) then
			    yzel2 = gl_RAY_A(par_n_BS, j, k - 1)
			else
			    yzel2 = gl_RAY_A(par_n_BS, j, N3)
			end if
			Ck(1) = gl_x2(yzel2, now); Ck(2) = gl_y2(yzel2, now); Ck(3) = gl_z2(yzel2, now)
			! Ck = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if (j > 1) then
			    yzel2 = gl_RAY_A(par_n_BS, j - 1, k)
			else
				yzel2 = gl_RAY_A(par_n_BS, j, k + ceiling(1.0 * N3/2))
			end if
			Dk(1) = gl_x2(yzel2, now); Dk(2) = gl_y2(yzel2, now); Dk(3) = gl_z2(yzel2, now)
			!Dk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if(k < N3) then
			    yzel2 = gl_RAY_A(par_n_BS, j, k + 1)
			else
				yzel2 = gl_RAY_A(par_n_BS, j, 1)
			end if
			Ek(1) = gl_x2(yzel2, now); Ek(2) = gl_y2(yzel2, now); Ek(3) = gl_z2(yzel2, now)
			! Ek = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if (gl_Point_num(yzel) > 0) then
			    vel = gl_Point_num(yzel) * par_nat_BS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time
			else
				vel = par_nat_BS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time
			end if
			
			
			gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
			gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
			gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
			
		end do
	end do
	
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc start: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync start: ', cudaGetErrorString(cudaGetLastError())
	
	N3 = size(gl_RAY_B(1, 1, :))
    N2 = size(gl_RAY_B(1, :, 1))
    N1 = size(gl_RAY_B(:, 1, 1))

	a1 = N1
	a2 = N2
	a3 = N3
	print*, "START A2 ", a1, a2, a3
    ! Цикл движения точек на лучах B  
	!$cuf kernel do(2) <<<*, *>>>
    do k = 1, N3
        do j = 1, N2
			
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
			    vel = gl_Point_num(yzel) * par_nat_TS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time
			else
				vel = par_nat_TS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time
			end if
			
			
			gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
			gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
			gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
			
			
			! Контакт
			
			yzel = gl_RAY_B(par_n_HP, j, k)
			Ak(1) = gl_x2(yzel, now); Ak(2) = gl_y2(yzel, now); Ak(3) = gl_z2(yzel, now)
			! Ak = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/)
			
			if (j < N2) then
			    yzel2 = gl_RAY_B(par_n_HP, j + 1, k)
			else
				!yzel2 = gl_RAY_O(1, 1, k)
				CYCLE
			end if
			Bk(1) = gl_x2(yzel2, now); Bk(2) = gl_y2(yzel2, now); Bk(3) = gl_z2(yzel2, now)
			! Bk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if (k > 1) then
			    yzel2 = gl_RAY_B(par_n_HP, j, k - 1)
			else
			    yzel2 = gl_RAY_B(par_n_HP, j, N3)
			end if
			Ck(1) = gl_x2(yzel2, now); Ck(2) = gl_y2(yzel2, now); Ck(3) = gl_z2(yzel2, now)
			! Ck = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if (j > 1) then
			    yzel2 = gl_RAY_B(par_n_HP, j - 1, k)
			else
				yzel2 = gl_RAY_A(par_n_HP, size( gl_RAY_A(par_n_HP, :, k)), k)
			end if
			Dk(1) = gl_x2(yzel2, now); Dk(2) = gl_y2(yzel2, now); Dk(3) = gl_z2(yzel2, now)
			! Dk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if(k < N3) then
			    yzel2 = gl_RAY_B(par_n_HP, j, k + 1)
			else
				yzel2 = gl_RAY_B(par_n_HP, j, 1)
			end if
			Ek(1) = gl_x2(yzel2, now); Ek(2) = gl_y2(yzel2, now); Ek(3) = gl_z2(yzel2, now)
			! Ek = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if (gl_Point_num(yzel) > 0) then
			    vel = gl_Point_num(yzel) * par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time
			else
				vel = par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time
			end if
			
			
			gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
			gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
			gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
			
		end do
	end do
	
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc start: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync start: ', cudaGetErrorString(cudaGetLastError())
	
	N3 = size(gl_RAY_K(1, 1, :))
    N2 = size(gl_RAY_K(1, :, 1))
    N1 = size(gl_RAY_K(:, 1, 1))

	a1 = N1
	a2 = N2
	a3 = N3
	print*, "START A3 ", a1, a2, a3
    ! Цикл движения точек на лучах K  
	!$cuf kernel do(2) <<<*, *>>>
    do k = 1, N3
        do j = 1, N2
			
			if (k /= 1 .and. j == 1) then
                    CYCLE
			end if
			
			yzel = gl_RAY_K(par_n_TS, j, k)
			Ak(1) = gl_x2(yzel, now); Ak(2) = gl_y2(yzel, now); Ak(3) = gl_z2(yzel, now)
			! Ak = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/)
			
			if (j < N2) then
			    yzel2 = gl_RAY_K(par_n_TS, j + 1, k)
			else
				yzel2 = gl_RAY_B(par_n_TS, size(gl_RAY_B(par_n_TS, :, k)), k)
			end if
			Bk(1) = gl_x2(yzel2, now); Bk(2) = gl_y2(yzel2, now); Bk(3) = gl_z2(yzel2, now)
			! Bk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if (k > 1) then
			    yzel2 = gl_RAY_K(par_n_TS, j, k - 1)
			else
			    yzel2 = gl_RAY_K(par_n_TS, j, N3)
			end if
			Ck(1) = gl_x2(yzel2, now); Ck(2) = gl_y2(yzel2, now); Ck(3) = gl_z2(yzel2, now)
			! Ck = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if (j > 1) then
			    yzel2 = gl_RAY_K(par_n_TS, j - 1, k)
			else
				yzel2 = gl_RAY_K(par_n_TS, j, k + ceiling(1.0 * N3/2))
			end if
			Dk(1) = gl_x2(yzel2, now); Dk(2) = gl_y2(yzel2, now); Dk(3) = gl_z2(yzel2, now)
			! Dk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if(k < N3) then
			    yzel2 = gl_RAY_K(par_n_TS, j, k + 1)
			else
				yzel2 = gl_RAY_K(par_n_TS, j, 1)
			end if
			Ek(1) = gl_x2(yzel2, now); Ek(2) = gl_y2(yzel2, now); Ek(3) = gl_z2(yzel2, now)
			! Ek = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if (gl_Point_num(yzel) > 0) then
			    vel = gl_Point_num(yzel) * par_nat_TS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time
			else
				vel = par_nat_TS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time
			end if
			
			
			gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
			gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
			gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
			
		end do
	end do
	
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc start: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync start: ', cudaGetErrorString(cudaGetLastError())
	
	N3 = size(gl_RAY_O(1, 1, :))
    N2 = size(gl_RAY_O(1, :, 1))
    N1 = size(gl_RAY_O(:, 1, 1))

    ! Цикл движения точек на лучах O  
	if (.False.) then   ! Можно отключить этот цикл
		
		!$cuf kernel do(2) <<<*, *>>>
        do k = 1, N3
            do j = 1, N2
			
			    ! Контакт
			    yzel = gl_RAY_O(1, j, k)
				Ak(1) = gl_x2(yzel, now); Ak(2) = gl_y2(yzel, now); Ak(3) = gl_z2(yzel, now)
			    ! Ak = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/)
			
			    if (j < N2) then
			        yzel2 = gl_RAY_O(1, j + 1, k)
			    else
				    yzel2 = yzel
				end if
				Bk(1) = gl_x2(yzel2, now); Bk(2) = gl_y2(yzel2, now); Bk(3) = gl_z2(yzel2, now)
			    ! Bk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			    if (k > 1) then
			        yzel2 = gl_RAY_O(1, j, k - 1)
			    else
			        yzel2 = gl_RAY_O(1, j, N3)
				end if
				Ck(1) = gl_x2(yzel2, now); Ck(2) = gl_y2(yzel2, now); Ck(3) = gl_z2(yzel2, now)
			    ! Ck = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			    if (j > 1) then
			        yzel2 = gl_RAY_O(1, j - 1, k)
			    else
				    CYCLE
				    !yzel2 = gl_RAY_C(1, size(gl_RAY_C(1, :, k)), k)
				end if
				Dk(1) = gl_x2(yzel2, now); Dk(2) = gl_y2(yzel2, now); Dk(3) = gl_z2(yzel2, now)
			    ! Dk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			    if(k < N3) then
			        yzel2 = gl_RAY_O(1, j, k + 1)
			    else
				    yzel2 = gl_RAY_O(1, j, 1)
				end if
				Ek(1) = gl_x2(yzel2, now); Ek(2) = gl_y2(yzel2, now); Ek(3) = gl_z2(yzel2, now)
			    ! Ek = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			    if (gl_Point_num(yzel) > 0) then
			        vel = gl_Point_num(yzel) * par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time
			    else
				    vel = par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time
			    end if
			
			
			    gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
			    gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
			    gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
			
		    end do
	    end do
	end if
	
	
	
	! Конец поверхностного натяжения --------------------------------------------------------------------------------------------------------------------
	
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc start: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync start: ', cudaGetErrorString(cudaGetLastError())
    
    N3 = size(gl_RAY_A(1, 1, :))
    N2 = size(gl_RAY_A(1, :, 1))
    N1 = size(gl_RAY_A(:, 1, 1))

	a1 = N1
	a2 = N2
	a3 = N3
	print*, "START A4 ", a1, a2, a3
	
    ! Цикл движения точек на лучах А  ************************************************************
	!$cuf kernel do(2) <<<*, *>>>
    do k = 1, N3
        do j = 1, N2
            
            if (k /= 1 .and. j == 1) then
                    CYCLE
            end if
            
            ! Вычисляем координаты текущего луча в пространстве
            the = (j - 1) * par_pi_8/2.0/(N2 - 1)
            phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
            
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
            
            ! HP
            yzel = gl_RAY_A(par_n_HP, j, k)
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
            R_HP = norm2(ER2)  ! Новое расстояние до HP
            
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

                ! до TS
                if (i <= par_n_TS) then  ! До расстояния = R_TS
                    r =  par_R0 + (R_TS - par_R0) * (DBLE(i)/par_n_TS)**par_kk1
                    !print *, r
                    !pause
                else if (i <= par_n_HP) then  
                    r = R_TS + (i - par_n_TS) * (R_HP - R_TS)/(par_n_HP - par_n_TS)
                else if (i <= par_n_BS) then 
                    r = R_HP + (i - par_n_HP) * (R_BS - R_HP)/(par_n_BS - par_n_HP)
                else
                    r = R_BS + (par_R_END - R_BS) * (DBLE(i- par_n_BS)/(par_n_END - par_n_BS))**par_kk2
                end if

                ! Записываем новые координаты
                gl_x2(yzel, now2) = r * cos(the)
                gl_y2(yzel, now2) = r * sin(the) * cos(phi)
                gl_z2(yzel, now2) = r * sin(the) * sin(phi)
                
            end do
        end do
	end do
    
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc start: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync start: ', cudaGetErrorString(cudaGetLastError())
	
    N3 = size(gl_RAY_B(1, 1, :))
    N2 = size(gl_RAY_B(1, :, 1))
    N1 = size(gl_RAY_B(:, 1, 1))

	a1 = N1
	a2 = N2
	a3 = N3
	print*, "START A5 ", a1, a2, a3
    ! Цикл движения точек на лучах B  ************************************************************
	!$cuf kernel do(2) <<<*, *>>>
    do k = 1, N3
        do j = 1, N2
            
            ! Вычисляем координаты текущего луча в пространстве
            the = par_pi_8/2.0 + (j) * par_triple_point/(N2)
            phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
            
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
			ER2 = proect * ER
			ER2  = KORD + ER2
            R_TS = norm2(ER2)  ! Новое расстояние до TS
            
            ! HP
            yzel = gl_RAY_B(par_n_HP, j, k)
            if(gl_Point_num(yzel) == 0) then
                vel = 0.0
			else
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
            R_HP = norm2(ER2)  ! Новое расстояние до HP
                
            do i = 1, N1

                if (i == 1) CYCLE
                
                yzel = gl_RAY_B(i, j, k)
                ! Вычисляем координаты точки на луче

                ! до TS
                if (i <= par_n_TS) then  ! До расстояния = R_TS
                    r =  par_R0 + (R_TS - par_R0) * (REAL(i, KIND = 4)/par_n_TS)**par_kk1
                else if (i <= par_n_HP) then  ! До расстояния = par_R_character * 1.3
                    r = R_TS + (i - par_n_TS) * (R_HP - R_TS) /(par_n_HP - par_n_TS)
                end if

                ! Записываем новые координаты
                gl_x2(yzel, now2) = r * cos(the)
                gl_y2(yzel, now2) = r * sin(the) * cos(phi)
                gl_z2(yzel, now2) = r * sin(the) * sin(phi)
                
               
            end do
        end do
	end do
	
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc start: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync start: ', cudaGetErrorString(cudaGetLastError())
    
    N3 = size(gl_RAY_C(1, 1, :))
    N2 = size(gl_RAY_C(1, :, 1))
    N1 = size(gl_RAY_C(:, 1, 1))

	a1 = N1
	a2 = N2
	a3 = N3
	print*, "START A6 ", a1, a2, a3
    ! Цикл движения точек на лучах C  ************************************************************
	!$cuf kernel do(2) <<<*, *>>>
    do k = 1, N3
        do j = 1, N2
            
            ! Вычисляем координаты текущего луча в пространстве
                phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
                
                ! Вычисляем координаты точки на луче
                x = gl_x2(gl_RAY_B(par_n_HP, j, k), now2)
                y = gl_y2(gl_RAY_B(par_n_HP, j, k), now2)
                z = gl_z2(gl_RAY_B(par_n_HP, j, k), now2)
                rr = (y**2 + z**2)**(0.5)
                
                ! BS     Нужно взять положение BS из её положения на крайнем луче A
                yzel = gl_RAY_A(par_n_BS, size(gl_RAY_A(1, :, 1)), k)
				ER(1) = 0.0_8; ER(2) = gl_y2(yzel, now2); ER(3) = gl_z2(yzel, now2)
                R_BS = norm2(ER)  ! Новое расстояние до BS
            
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
        end do
	end do
	
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc start: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync start: ', cudaGetErrorString(cudaGetLastError())

    N3 = size(gl_RAY_O(1, 1, :))
    N2 = size(gl_RAY_O(1, :, 1))
    N1 = size(gl_RAY_O(:, 1, 1))

	a1 = N1
	a2 = N2
	a3 = N3
	print*, "START A7 ", a1, a2, a3

    ! Цикл движения точек на лучах O   ************************************************************
	!$cuf kernel do(2) <<<*, *>>>
    do k = 1, N3
        do j = 1, N2
            phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
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
            
            xx = gl_x2(gl_RAY_B(par_n_HP, par_m_BC, k), now2)              ! Отталкиваемся от x - координаты крайней точки B на гелиопаузе в этой плоскости (k)
            x = xx - (DBLE(j)/N2)**par_kk3 * (xx - par_R_LEFT)
            
            ! BS     Нужно взять положение BS из её положения на крайнем луче A
            yzel = gl_RAY_A(par_n_BS, size(gl_RAY_A(1, :, 1)), k)
			KORD(1) = 0.0_8; KORD(2) = gl_y2(yzel, now2); KORD(3) = gl_z2(yzel, now2)
            R_BS = norm2(KORD)  ! Новое расстояние до BS
            
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
        end do
	end do
    
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc start: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync start: ', cudaGetErrorString(cudaGetLastError())
	
    N3 = size(gl_RAY_K(1, 1, :))
    N2 = size(gl_RAY_K(1, :, 1))
    N1 = size(gl_RAY_K(:, 1, 1))

	a1 = N1
	a2 = N2
	a3 = N3
	print*, "START A8 ", a1, a2, a3
	
    ! Цикл движения точек на лучах K  ************************************************************
	!$cuf kernel do(2) <<<*, *>>>
    do k = 1, N3
        do j = 1, N2
            
             ! Вычисляем координаты текущего луча в пространстве
            the = par_pi_8/2.0 + par_triple_point + (N2 - j + 1) * (par_pi_8/2.0 - par_triple_point)/(N2)
            phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
            
            if (k /= 1 .and. j == 1) CYCLE
            
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
                r =  par_R0 + (R_TS - par_R0) * (REAL(i, KIND = 4)/par_n_TS)**par_kk1


                ! Записываем новые координаты
                gl_x2(yzel, now2) = r * cos(the)
                gl_y2(yzel, now2) = r * sin(the) * cos(phi)
                gl_z2(yzel, now2) = r * sin(the) * sin(phi)
                
                
            end do
        end do
	end do
    
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc start: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync start: ', cudaGetErrorString(cudaGetLastError())
    
     N3 = size(gl_RAY_D(1, 1, :))
    N2 = size(gl_RAY_D(1, :, 1))
    N1 = size(gl_RAY_D(:, 1, 1))

	a1 = N1
	a2 = N2
	a3 = N3
	print*, "START A9 ", a1, a2, a3
	
    ! Цикл движения точек на лучах D ************************************************************
	!$cuf kernel do(2) <<<*, *>>>
    do k = 1, N3
        do j = 1, N2
            
            ! Вычисляем координаты текущего луча в пространстве
                phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
                
            do i = 1, N1

                if (i == 1) CYCLE

                if (k /= 1 .and. j == 1) CYCLE

                
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

                yzel = gl_RAY_D(i, j, k)
                r = sqrt(y**2 + z**2)
                
                x = xx + (DBLE(i - 1)/(N1 - 1))**par_kk3 * (par_R_LEFT - xx)

                ! Записываем новые координаты
                gl_x2(yzel, now2) = x
                gl_y2(yzel, now2) = r * cos(phi)
                gl_z2(yzel, now2) = r * sin(phi)
                
            end do
        end do
	end do
	
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc start: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync start: ', cudaGetErrorString(cudaGetLastError())
    
    N3 = size(gl_RAY_E(1, 1, :))
    N2 = size(gl_RAY_E(1, :, 1))
    N1 = size(gl_RAY_E(:, 1, 1))
    
	a1 = N1
	a2 = N2
	a3 = N3
	print*, "START A10 ", a1, a2, a3
	
    ! Цикл движения точек на лучах E  ************************************************************
	!$cuf kernel do(2) <<<*, *>>>
    do k = 1, N3
        do j = 1, N2
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
        end do
	end do
	
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc start: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync start: ', cudaGetErrorString(cudaGetLastError())
	
	end subroutine Cuda_Move_all
	
	
	
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
	
	subroutine Cuda_Calc_move(now)
	! Аналог функции на хосте Calc_move
	use MY_CUDA, gl_TS => dev_gl_TS, gl_Gran_neighbour => dev_gl_Gran_neighbour, gl_Gran_normal2 => dev_gl_Gran_normal2, &
		gl_Cell_par => dev_gl_Cell_par, gl_all_Gran => dev_gl_all_Gran, gl_Point_num => dev_gl_Point_num, gl_x2 => dev_gl_x2, &
		gl_y2 => dev_gl_y2, gl_z2 => dev_gl_z2, gl_Gran_center2 => dev_gl_Gran_center2, norm2 => dev_norm2, gl_Vx => dev_gl_Vx, &
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS
	use GEO_PARAM
	implicit none
	integer, intent(in) :: now
	
	integer i, Num, gr, s1, s2, j, yzel
    real(8), device :: normal(3), qqq1(8), qqq2(8), POTOK(8), ray(3)
    real(8) :: a1, a2, v1, v2, norm, b1, b2, c1, c2, koef1, koef2, koef3
	real(8) :: www, dsl, dsp, dsc
	integer(4), device :: metod
	
    koef1 = 0.3
    koef2 = 0.7
    koef3 = 0.3
	
	metod = 2
	www = 0.0
    
    ! Пробегаемся по всем граням и вычисляем скорости их движения
    
    ! TS
    Num = size(gl_TS)
	
	write (*,*) "START dw234"
    
	!$cuf kernel do <<<*, *>>>
    do i = 1, Num
        dsl = 0.0
		dsp = 0.0
		dsc = 0.0
		!gr = 0
		!s1 = 0
		!s2 = 0
		!normal = 0.0
		!qqq1 = 0.0
		!qqq2 = 0.0
		!POTOK = 0.0
		!a1 = 0.0
		!a2 = 0.0
		!yzel = 0
		!metod = 2
		!www = 0.0
		
        gr = gl_TS(i)
    	s1 = gl_Gran_neighbour(1, gr)
        s2 = gl_Gran_neighbour(2, gr)
        normal = gl_Gran_normal2(:, gr, now)
        qqq1 = gl_Cell_par(1:8, s1)
        qqq2 = gl_Cell_par(1:8, s2)
		
        
        !call chlld(metod, normal(1), normal(2), normal(3), &
        !        www, qqq1, qqq2, dsl, dsp, dsc, POTOK)
		
        dsl = dsl * koef1
        
        a1 = sqrt(ggg * qqq1(5)/qqq1(1))  ! Скорости звука
        a2 = sqrt(ggg * qqq2(5)/qqq2(1))
        
        if (norm2(qqq1(2:4)) <= a1 .and. norm2(qqq2(2:4)) <= a2) then ! Записываем скорость во все узлы
            do j = 1, 4
            	yzel = gl_all_Gran(j, gr)
                gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsl
                gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsl
                gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsl
                gl_Point_num(yzel) = gl_Point_num(yzel) + 1
            end do
            CYCLE  ! Заканчиваем с этой гранью, переходим к следующей
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
                    gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsl
                    gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsl
                    gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsl
                    gl_Point_num(yzel) = gl_Point_num(yzel) + 1
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
                    gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsl
                    gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsl
                    gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsl
                    gl_Point_num(yzel) = gl_Point_num(yzel) + 1
                    CYCLE
                end if
                
        end do
        
	end do
        
	write (*,*) "STOP dw234"

    ! Контакт
    Num = size(gl_Contact)
    
	!$cuf kernel do <<<*, *>>>
    do i = 1, Num
        
        gr = gl_Contact(i)
    	s1 = gl_Gran_neighbour(1, gr)
        s2 = gl_Gran_neighbour(2, gr)
        normal = gl_Gran_normal2(:, gr, now)
        qqq1 = gl_Cell_par(1:8, s1)
        qqq2 = gl_Cell_par(1:8, s2)
        
        !call chlld(metod, normal(1), normal(2), normal(3), &
        !        www, qqq1, qqq2, dsl, dsp, dsc, POTOK)
        
        dsc = dsc * koef2
        
        a1 = sqrt(ggg * qqq1(5)/qqq1(1))  ! Скорости звука
        a2 = sqrt(ggg * qqq2(5)/qqq2(1))
        
        if (norm2(qqq1(2:4)) <= a1 .and. norm2(qqq2(2:4)) <= a2) then ! Записываем скорость во все узлы
            do j = 1, 4
            	yzel = gl_all_Gran(j, gr)
                gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsc
                gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsc
                gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsc
                gl_Point_num(yzel) = gl_Point_num(yzel) + 1
            end do
            CYCLE  ! Заканчиваем с этой гранью, переходим к следующей
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
                    gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsc
                    gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsc
                    gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsc
                    gl_Point_num(yzel) = gl_Point_num(yzel) + 1
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
                    gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsc
                    gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsc
                    gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsc
                    gl_Point_num(yzel) = gl_Point_num(yzel) + 1
                    CYCLE
                end if
                
        end do
        
    end do
    
    
    ! BS
    Num = size(gl_BS)
    
	!$cuf kernel do <<<*, *>>>
    do i = 1, Num
        
        gr = gl_BS(i)
    	s1 = gl_Gran_neighbour(1, gr)
        s2 = gl_Gran_neighbour(2, gr)
        normal = gl_Gran_normal2(:, gr, now)
        qqq1 = gl_Cell_par(1:8, s1)
        qqq2 = gl_Cell_par(1:8, s2)
        
        !call chlld(metod, normal(1), normal(2), normal(3), &
        !        www, qqq1, qqq2, dsl, dsp, dsc, POTOK)
        
        dsp = dsp * koef3
        
        a1 = sqrt(ggg * qqq1(5)/qqq1(1))  ! Скорости звука
        a2 = sqrt(ggg * qqq2(5)/qqq2(1))
        
        if (norm2(qqq1(2:4)) <= a1 .and. norm2(qqq2(2:4)) <= a2) then ! Записываем скорость во все узлы
            do j = 1, 4
            	yzel = gl_all_Gran(j, gr)
                gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsp
                gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsp
                gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsp
                gl_Point_num(yzel) = gl_Point_num(yzel) + 1
            end do
            CYCLE  ! Заканчиваем с этой гранью, переходим к следующей
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
                    gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsp
                    gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsp
                    gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsp
                    gl_Point_num(yzel) = gl_Point_num(yzel) + 1
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
                    gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsp
                    gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsp
                    gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsp
                    gl_Point_num(yzel) = gl_Point_num(yzel) + 1
                    CYCLE
                end if
                
        end do
        
    end do
	
	
	end subroutine Cuda_Calc_move
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	attributes(global) subroutine Cuda_Calc_move_TS(now)
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
    real(8) :: a1, a2, v1, v2, norm, b1, b2, c1, c2, koef1, koef2, koef3
	real(8) :: www, dsl, dsp, dsc
	integer(4) :: metod
	
	i = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	Num = size(gl_TS)
	
	
	
	if (i > Num) return
	
    koef1 = 0.3
    koef2 = 0.7
    koef3 = 0.3
	metod = 2
	www = 0.0
	
	!dsl = 0.0
	!dsp = 0.0
	!dsc = 0.0
	!gr = 0
	!s1 = 0
	!s2 = 0
	!normal = 0.0
	!qqq1 = 0.0
	!qqq2 = 0.0
	!POTOK = 0.0
	!a1 = 0.0
	!a2 = 0.0
	!yzel = 0
		
    gr = gl_TS(i)
    s1 = gl_Gran_neighbour(1, gr)
    s2 = gl_Gran_neighbour(2, gr)
    normal = gl_Gran_normal2(:, gr, now)
    qqq1 = gl_Cell_par(1:8, s1)
    qqq2 = gl_Cell_par(1:8, s2)
		
        
    call chlld(metod, normal(1), normal(2), normal(3), &
            www, qqq1, qqq2, dsl, dsp, dsc, POTOK)
	
		
    dsl = dsl * koef1
        
    a1 = sqrt(ggg * qqq1(5)/qqq1(1))  ! Скорости звука
    a2 = sqrt(ggg * qqq2(5)/qqq2(1))
	
        
    if (norm2(qqq1(2:4)) <= a1 .and. norm2(qqq2(2:4)) <= a2) then ! Записываем скорость во все узлы
        do j = 1, 4
            yzel = gl_all_Gran(j, gr)
            gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsl
            gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsl
            gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsl
            gl_Point_num(yzel) = gl_Point_num(yzel) + 1
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
				
				do while ( atomiccas(dev_mutex_1, 0, 1) == 1)  ! Безопасный доступ к памяти dev_mutex_1
				end do                                          ! Безопасный доступ к памяти
                gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsl
                gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsl
                gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsl
                gl_Point_num(yzel) = gl_Point_num(yzel) + 1
				call threadfence()                               ! Безопасный доступ к памяти
				dev_mutex_1 = 0                                  ! Безопасный доступ к памяти
				
				
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
                gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsl
                gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsl
                gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsl
                gl_Point_num(yzel) = gl_Point_num(yzel) + 1
				
                CYCLE
            end if
                
    end do
	
	end subroutine Cuda_Calc_move_TS
	
	
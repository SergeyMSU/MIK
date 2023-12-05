module PUI
	! Модуль интерполяции по тетраэдрам
    ! Строится двойственная сетка к основной (с дополнительными узлами на разрывах)
	! Строиться сетка тетраэдров на двойственной сетке
	USE GEO_PARAM
	USE STORAGE
	USE My_func
	USE Solvers
	USE Interpolate2
	USE OMP_LIB
	implicit none
	
	integer :: pui_nW = 40      ! 50
	real(8) :: pui_wR = 100.0    ! 150.0
	
	integer, allocatable :: pui_num_tetr(:)        ! По номеру интерполяционного тетраэдра определяем номер S
	integer, allocatable :: pui_num_tetr_2(:)	   ! Обратно по номеру S определяем номер интерполяционного тетраэдра
	real(8), allocatable :: pui_Sm(:, :)           ! (pui_nW, :) 
	real(8), allocatable :: pui_Sp(:, :)           ! (pui_nW, :) 
	real(8), allocatable :: f_pui(:, :)           ! (pui_nW, :) 
	real(8), allocatable :: f_pui_num(:)           ! По номеру в массиве пуи, определяем номер узла в интерполяционной сетке
	real(8), allocatable :: f_pui_num2(:)		   ! По номеру узла в интерполяционной сетке, определяем номер в массиве PUI 
	integer (kind=omp_lock_kind), allocatable :: pui_lock(:)  ! Для openMP
	
	contains
	
	subroutine PUI_Set()
		! Выделяем память под все массивы
		integer :: N, i, k, j, N2
		! Сначала нужно понять, сколько у нас внутренних точек, в которых мы считаем PUI
		! Для этого пробегаемся по массивам тетраэдров
		N = 0
		k = 1
		allocate(pui_num_tetr(size(int2_all_tetraendron_point(1, :))))
		pui_num_tetr = 0
		
		do i = 1, size(int2_all_tetraendron_point(1, :))
			j = int2_all_tetraendron_point(1, i)
			if(int2_Cell_par2(1, j) <= 2) then
				N = N + 1
				pui_num_tetr(i) = k
				k = k + 1
			else
				pui_num_tetr(i) = 0
			end if
		end do

		
		
		allocate(pui_Sm(pui_nW, N))
		allocate(pui_Sp(pui_nW, N))
		allocate(pui_num_tetr_2(N))
		allocate(pui_lock(N))
		
		
		do i = 1, N
			call omp_init_lock(pui_lock(i))
		end do
		
		
		pui_Sm = 0.0
		pui_Sp = 0.0
		f_pui = 0.0
		pui_num_tetr_2 = 0
		
		N = 0
		do i = 1, size(int2_all_tetraendron(1, :))
			j = int2_all_tetraendron_point(1, i)
			if(int2_Cell_par2(1, j) <= 2) then
				N = N + 1
				pui_num_tetr_2(N) = i
			end if
		end do
		
		
	end subroutine PUI_Set
	
	subroutine PUI_f_Set()
		! Создаём массивы для функции распределения
		integer :: N, i, k, j, N2
		
		allocate(f_pui_num2(size(int2_Cell_par2(1, :))))
		
		N2 = 0
		k = 1
		do i = 1, size(int2_Cell_par2(1, :))
			f_pui_num2(i) = 0
			if(int2_Cell_par2(1, i) <= 2) then
				N2 = N2 + 1
				f_pui_num2(i) = k
				k = k + 1
			end if
		end do
		
		allocate(f_pui(pui_nW, N2))
		allocate(f_pui_num(N2))
		
		
		k = 1
		do i = 1, size(int2_Cell_par2(1, :))
			if(int2_Cell_par2(1, i) <= 2) then
				f_pui_num(k) = i
				k = k + 1
			end if
		end do
	
	end subroutine PUI_f_Set
	
	subroutine Culc_f_pui()
		integer :: N, i, k, num, yzel, tetraedron, iw, numw, num_do, i1, num2
		real(8) :: r(3), PAR(9), dt, Sm, w0, q1, rho, rho0, w, qInt, Sp, rho_do
		real(8) :: s, cospsi, C, A, B, Sm2, Sp2, qInt2
		real(8) :: PAR_k(5), normal(3)
		real(8) :: f0_pui(pui_nW)
		real(8) :: mas_w(pui_nW)
		real(8) :: mas_w0(pui_nW)
		real(8) :: mas_Sm(pui_nW)
		logical :: find_n

		pui_Sm = pui_Sm * par_n_p_LISM
		pui_Sp = pui_Sp * par_n_p_LISM

		
		print*, "START Culc_f_pui"

		dt = 0.0005
		do iw = 1, pui_nW
				mas_w0(iw) = ((iw - 0.5) * pui_wR / pui_nW)
		end do

		! Находим функцию распределения для ячеек перед ударной волной
		do i = 1, size(f_pui_num)
			if(mod(i, 5000) == 0) print*, "i = ", i, " from", size(f_pui_num)
			k = f_pui_num(i)      ! Номер узла сетки интерполяции, в которой считаем PUI
			if(int2_Cell_par2(1, k) == 2) CYCLE
			rho0 = int2_Cell_par(1, k)
			f0_pui = 0.0
			mas_Sm = 0.0

			r = int2_coord(:, k)  ! Координаты этого узла
			mas_w = mas_w0

			num = 3
			qInt = 0.0            ! Интеграл от источника массы при ионизации
				
			! Бежим до Солнца
			do while (.TRUE.)
				call Int2_Get_par_fast(r(1), r(2), r(3), num, PAR, PAR_k = PAR_k)
				q1 = PAR_k(1)
				rho = PAR(1)
				qInt = qInt + dt * q1/rho                                    
				r = r - PAR(2:4) * dt
				tetraedron = pui_num_tetr(num) ! Номер тетраэдра в массиве источников
				
				do iw = 1, pui_nW
					!print*, iw, mas_w(iw)
					numw = min(INT(mas_w(iw)/pui_wR * pui_nW) + 1, pui_nW)
					!print*, iw, f0_pui(iw), pui_Sp(numw, tetraedron), dt
					if(mas_w(iw) < pui_wR .and. mas_w(iw) > 0) then
						mas_Sm(iw) = mas_Sm(iw) + pui_Sm(numw, tetraedron) * dt                         
						f0_pui(iw) = f0_pui(iw) + pui_Sp(numw, tetraedron) * dt * exp(-mas_Sm(iw))  ! Это S+, просто сразу накапливаем в функцию распределения
					end if
					!print*, iw, f0_pui(iw)
					mas_w(iw) = mas_w0(iw) / ( (rho0/rho)**(1.0/3.0) * exp(-1.0/3.0 * qInt) )
					!print*, iw, mas_w(iw)
					!pause
				end do
				
				if(norm2(r) <= par_R0) EXIT
			end do

			f_pui(:, i) = f0_pui(:)

			! do iw = 1, pui_nW
			! 	r = int2_coord(:, k)  ! Координаты этого узла
			! 	w0 = ((iw - 0.5) * pui_wR / pui_nW)
			! 	w = w0				
			! 	num = 3
			! 	dt = 0.0005
			! 	Sm = 0.0
			! 	Sp = 0.0
			! 	qInt = 0.0            ! Интеграл от источника массы при ионизации
			! 	! Бежим до Солнца
			! 	do while (.TRUE.)
			! 		num_do = num
			! 		!print*, r, num
			! 		!pause
			! 		call Int2_Get_par_fast(r(1), r(2), r(3), num, PAR, PAR_k = PAR_k)
			! 		q1 = PAR_k(1)
			! 		rho_do = rho
			! 		rho = PAR(1)
			! 		qInt = qInt + dt * q1/rho                                              !TODO
			! 		r = r - PAR(2:4) * dt
			! 		tetraedron = pui_num_tetr(num) ! Номер тетраэдра в массиве источников
			! 		numw = min(INT(w/pui_wR * pui_nW) + 1, pui_nW)
			! 		if(w < pui_wR .and. w > 0) then
			! 			Sm = Sm + pui_Sm(numw, tetraedron) * dt                               !TODO
			! 			Sp = Sp + pui_Sp(numw, tetraedron) * dt * exp(-Sm)
			! 		end if
					
			! 		w = w0 / ( (rho0/rho)**(1.0/3.0) * exp(-1.0/3.0 * qInt) )
			! 		!w = w + dt * w/3.0 * 2.0 * norm2(PAR(2:4))/norm2(r)
					
			! 		if(norm2(r) <= par_R0) EXIT
			! 	end do
				
			! 	f_pui(iw, i) = Sp
				
			! end do
		end do
		
		return

		! Находим функцию распределения для ячеек перед ударной волной
		do i = 1, size(f_pui_num)
			if(mod(i, 5000) == 0) print*, "i = ", i, " from", size(f_pui_num)
			k = f_pui_num(i)      ! Номер узла сетки интерполяции, в которой считаем PUI
			if(int2_Cell_par2(1, k) == 2) CYCLE
			rho0 = int2_Cell_par(1, k)
		
				
			do iw = 1, pui_nW
				r = int2_coord(:, k)  ! Координаты этого узла
				w0 = ((iw - 0.5) * pui_wR / pui_nW)
				w = w0				
				num = 3
				dt = 0.0005
				Sm = 0.0
				Sp = 0.0
				qInt = 0.0            ! Интеграл от источника массы при ионизации
				! Бежим до Солнца
				do while (.TRUE.)
					num_do = num
					!print*, r, num
					!pause
					call Int2_Get_par_fast(r(1), r(2), r(3), num, PAR, PAR_k = PAR_k)
					q1 = PAR_k(1)
					rho_do = rho
					rho = PAR(1)
					qInt = qInt + dt * q1/rho                                              !TODO
					r = r - PAR(2:4) * dt
					tetraedron = pui_num_tetr(num) ! Номер тетраэдра в массиве источников
					numw = min(INT(w/pui_wR * pui_nW) + 1, pui_nW)
					if(w < pui_wR .and. w > 0) then
						Sm = Sm + pui_Sm(numw, tetraedron) * dt                               !TODO
						Sp = Sp + pui_Sp(numw, tetraedron) * dt * exp(-Sm)
					end if
					
					w = w0 / ( (rho0/rho)**(1.0/3.0) * exp(-1.0/3.0 * qInt) )
					!w = w + dt * w/3.0 * 2.0 * norm2(PAR(2:4))/norm2(r)
					
					if(norm2(r) <= par_R0) EXIT
				end do
				
				f_pui(iw, i) = Sp
				
			end do
		end do
		
		return
		
		! Находим функцию распределения для ячеек за ударной волной
		do i = 1, size(f_pui_num)
			
			k = f_pui_num(i)      ! Номер узла сетки интерполяции, в которой считаем PUI
			if(int2_Cell_par2(1, k) == 1) CYCLE
			r = int2_coord(:, k)  ! Координаты этого узла
			rho0 = int2_Cell_par(1, k)
		
				
			do iw = 1, pui_nW
				w0 = ((iw - 0.5) * pui_wR / pui_nW)
				w = w0				
				num = 3
				dt = 0.001
				Sm = 0.0
				Sp = 0.0
				qInt = 0.0            ! Интеграл от источника массы при ионизации
				! Бежим до ударной волны
				do while (.TRUE.)
					num_do = num
					call Int2_Get_par_fast(r(1), r(2), r(3), num, PAR, PAR_k = PAR_k)
					q1 = PAR_k(1)
					rho_do = rho
					rho = PAR(1)
					qInt = qInt + dt * q1/rho
					r = r - PAR(2:4) * dt
					tetraedron = pui_num_tetr(num) ! Номер тетраэдра в массиве источников
					numw = min(INT(w/pui_wR * pui_nW) + 1, pui_nW)
					Sm = Sm + pui_Sm(numw, tetraedron)
					Sp = Sp + pui_Sp(numw, tetraedron) * exp(-Sm)
					
					w = w0 * (rho/rho0)**(1.0/3.0) * exp(-1.0/3.0 * qInt)
					
					yzel = int2_all_tetraendron_point(1, num)
					if(int2_Cell_par2(1, yzel) == 1) EXIT
				end do
				
				! Теперь надо получить нормаль к поверхности TS
				find_n = .FALSE.
				do i1 = 1, 4
					num2 = int2_gran_sosed(int2_all_tetraendron(i1, num))
					if(num2 == num_do) then
						normal = int2_plane_tetraendron(1:3, i1, num)
						find_n = .TRUE.
						EXIT
					end if
				end do
				
				if(find_n /= .TRUE.) then
					print*, "Ne nashol normal"
					normal = r
				end if
				
				s = rho_do/rho
				cospsi = DOT_PRODUCT(PAR(6:8), normal)/(norm2(normal) * norm2(PAR(6:8)))
				A = sqrt(cospsi**2 + s**2 * (1.0 - cospsi**2))
				B = s**2 / A**2
				C = (2.0 * A + B)/3.0
				
				! Бежим до бежим после ударной волны
				do while (.TRUE.)
					num_do = num
					call Int2_Get_par_fast(r(1), r(2), r(3), num, PAR, PAR_k = PAR_k)
					q1 = PAR_k(1)
					rho_do = rho
					rho = PAR(1)
					qInt = qInt + dt * q1/rho
					r = r - PAR(2:4) * dt
					tetraedron = pui_num_tetr(num) ! Номер тетраэдра в массиве источников
					numw = min(INT(w/pui_wR * pui_nW) + 1, pui_nW)
					Sm = Sm + pui_Sm(numw, tetraedron)
					Sp = Sp + pui_Sp(numw, tetraedron) * exp(-Sm)
					
					w = w0 * (rho/rho0)**(1.0/3.0) * exp(-1.0/3.0 * qInt)
					
					yzel = int2_all_tetraendron_point(1, num)
					if(int2_Cell_par2(1, yzel) == 1) EXIT
				end do
				
			end do
		end do
	
	end subroutine Culc_f_pui
	
	subroutine PUI_Add(cell, wr, nu_ex, nu_ph, mu, time)
		! wr - скорость атома в СК, связанной со средней скоростью плазмы (модуль этой скорости)
		integer, intent(in) :: cell
		real(8), intent(in) :: wr, nu_ex, nu_ph, mu, time
		integer i, j
		
		j = pui_num_tetr(cell)
		if(j > 0) then
			i = min(INT(wr/pui_wR * pui_nW) + 1, pui_nW)
			
			call omp_set_lock(pui_lock(j))
			pui_Sm(i, j) = pui_Sm(i, j) + mu * time
			pui_Sp(i, j) = pui_Sp(i, j) + mu * time * nu_ex
		
			!pui_Sm(i, j) = pui_Sm(i, j) + mu * time
			pui_Sp(i, j) = pui_Sp(i, j) + mu * time * nu_ph
			call omp_unset_lock(pui_lock(j))
			
		end if
		
	
	end subroutine PUI_Add
	
	subroutine PUI_calc_Sm()
		! Расчёт Sm - он не считается "на лету" в Монте-карло и нужна постобработка
		USE OMP_LIB
		real(8) :: pui_Sm2(pui_nW)           ! (pui_nW, :, potok) 
		integer :: ij, i, j, k, num_all
		real(8) :: dthe, Vh, ff, the, d, w
		
		dthe = par_pi_8/40.0
		num_all = 0
		
	 	!$omp parallel
        
	 	!$omp do private(pui_Sm2, ij, j, k, ff, Vh, the, d, w)
		do i = 1, size(pui_Sm(1, :))

			!$omp critical
			num_all = num_all + 1
			if(mod(num_all, 5000) == 0) then
				print*, num_all, "from", size(pui_Sm(1, :))
			end if
			!$omp end critical
			
			pui_Sm2 = 0.0
			do ij = 1, pui_nW
				w = ((ij - 0.5) * pui_wR / pui_nW);
				do j = 1, pui_nW
					Vh = ((j - 0.5) * pui_wR / pui_nW);
					ff = pui_Sm(j, i)
					if (ff <= 0.0) CYCLE
					do k = 1, 40
						the = dthe * k
						d = sqrt(Vh**2 + w**2 - 2.0 * w * Vh * cos(the))
						if (d > 0.000000001) then
							pui_Sm2(ij) = pui_Sm2(ij) + ff * d * MK_sigma(d) * sin(the) * dthe * 2.0 * par_pi_8
						end if
					end do
				end do
				pui_Sm2(ij) = pui_Sm2(ij) / (4.0 * par_pi_8)
			end do
			
			pui_Sm(:, i) = pui_Sm2/par_Kn/2.0  !TODO НУЖНО БУДЕТ УБРАТЬ ДЕЛЕНИЕ НА ДВА В СЛЕДУЮЩИЙ РАЗ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		end do	
		!$omp end do
		!$omp end parallel
		
	end subroutine PUI_calc_Sm
	
	subroutine PUI_print(num, x, y, z)
		implicit none
		integer, intent(in) :: num
		real(8), intent(in) :: x, y, z
		real(8) :: w, S, SS
		character(len=5) :: name
		integer(4) :: n1, i, n2
    
		n1 = 3
		call Int2_Get_tetraedron_inner(x, y, z, n1)
		write(unit=name,fmt='(i5.5)') num
    
		open(1, file = "S+_" // name // ".txt")
		n2 = pui_num_tetr(n1)
		if(n2 > 0) then
			do i = 1, pui_nW
				w = (i-0.5) * pui_wR/pui_nW
				write(1,*) w, pui_Sp(i, n2)
			end do
		end if
		
		close(1)
		
		open(2, file = "S-_" // name // ".txt")
		if(n2 > 0) then
			do i = 1, pui_nW
				w = (i-0.5) * pui_wR/pui_nW
				write(2,*) w, pui_Sm(i, n2)
			end do
		end if
		
		close(2)
		
		n2 = f_pui_num2(int2_all_tetraendron_point(1, n1))
		print*, n1, n2
		
		open(3, file = "f_PUI_" // name // ".txt")
		if(n2 > 0) then
			do i = 1, pui_nW
				w = (i-0.5) * pui_wR/pui_nW
				write(3,*) w, f_pui(i, n2)
			end do
		end if
		
		close(3)

		S = 0.0
		!? Печатаем концентрацию пикапов (в виде графика от скорости, чтобы посмотреть, на каких скоростях накапливается максимум)
		open(4, file = "n_pui_" // name // ".txt")
		if(n2 > 0) then
			do i = 1, pui_nW
				w = (i-0.5) * pui_wR/pui_nW
				S = S + f_pui(i, n2) * 4 * par_pi_8 * w**2
				write(4,*) w, S
			end do
		end if
		
		close(4)

		SS = 0.0
		open(5, file = "T_pui_" // name // ".txt")
		if(n2 > 0) then
			do i = 1, pui_nW
				w = (i-0.5) * pui_wR/pui_nW
				SS = SS + f_pui(i, n2) * 4 * par_pi_8 * w**4 /(3.0 * S)
				write(5,*) w, SS
			end do
		end if
		
		close(5)
	
	end subroutine PUI_print
	
	subroutine PUI_Save_bin(num)
		implicit none
		integer, intent(in) :: num
		character(len=5) :: name
		
		write(unit=name,fmt='(i5.5)') num
		
		open(1, file = "pui_save_" // name // ".bin", FORM = 'BINARY')
		
		write(1) size(pui_num_tetr_2)
		write(1) pui_nW
		write(1) pui_wR
		
		write(1) size(pui_Sm(:, 1)), size(pui_Sm(1, :))
		write(1) pui_Sm(:, :)
		
		write(1) size(pui_Sp(:, 1)), size(pui_Sp(1, :))
		write(1) pui_Sp(:, :)
		
		
		
		write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
		write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
		write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
		write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
		write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
		write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
		write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
		write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
		write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
		write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
		
		close(1)
	
	end subroutine PUI_Save_bin
	
	subroutine PUI_Read_bin(num)
		implicit none
		integer, intent(in) :: num
		character(len=5) :: name
		integer(4) :: n1, n2, n3
		real(8) :: aa
		logical :: exists
    
		write(unit=name,fmt='(i5.5)') num
	
		inquire(file="pui_save_" // name // ".bin", exist=exists)
    
		if (exists == .False.) then
			pause "ERROR net faila 171 PUI!!!"
			STOP "ERROR net faila 171 PUI!!!"
		end if
    
		open(1, file = "pui_save_" // name // ".bin", FORM = 'BINARY', ACTION = "READ")
	
		read(1) n1
		if(size(pui_num_tetr_2) /= n1) then
			STOP "ERROR PUI size(pui_num_tetr_2) /= n1"
		end if
	
		read(1) n2
		read(1) aa
		if(pui_nW /= n2) then
			STOP "ERROR PUI pui_nW /= n2"
		end if
		if(pui_wR /= aa) then
			print*, pui_wR, aa
			STOP "ERROR pui_wR /= aa"
		end if
	
		read(1) n1, n2
		if(size(pui_Sm(:, 1)) /= n1) then
			STOP "ERROR PUI size(pui_Sm(:, 1)) /= n1"
		end if
	
		if(size(pui_Sm(1, :)) /= n2) then
			STOP "ERROR PUI size(pui_Sm(1, :)) /= n2"
		end if
	
		read(1) pui_Sm
		read(1) n1, n2
		read(1) pui_Sp
	
		close(1)
	
	end subroutine PUI_Read_bin
	
end module PUI
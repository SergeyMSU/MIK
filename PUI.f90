module PUI
	USE GEO_PARAM
	USE STORAGE
	USE My_func
	USE Solvers
	USE Interpolate2
	USE OMP_LIB
	implicit none

	!! Важная информация
	!! Источники S+ и S- считаются в тетраэдрах (я не стал их суммировать по тетраэрам)
	!! Функция f_pui считается в центрах ячеек исходной сетки! 
	!! Т.е. она считается в тех же точках, что и основная параметры (+ ещё дополнительные точки на разрывах посчитаны)
	!! Плюс ещё в том, что значения можно интерполировать по тетраэдрам 
	!TODO (нужно написать свою программу интерполяции)


    !! Ntetr - число тетраэдров во внутри HP (области 1 и 2)
	integer, allocatable :: pui_num_tetr(:)        ! По номеру интерполяционного тетраэдра определяем номер S
	integer, allocatable :: pui_num_tetr_2(:)	   ! Обратно по номеру S определяем номер интерполяционного тетраэдра
	real(8), allocatable :: pui_Sm(:, :)           ! (pui_nW, Ntetr)
	real(8), allocatable :: pui_Sp(:, :)           ! (pui_nW, Ntetr)
	integer (kind=omp_lock_kind), allocatable :: pui_lock(:)  ! Для openMP

	!! Переменные для обработки функции распределения PUI для последующего розыгрыша
	


	!! функция h0(U_H) - см. документацию PUI
	integer :: pui_h0_n = 1000      
	real(8), PARAMETER :: pui_h0_wc = 100.0    
	real(8), allocatable :: h0_pui(:)   ! Для функции отказов при розыгрыше - её максимум для нормировки (см. документацию)

	!! Розыгрышь пикапов (заранее считаем функцию розыгрыша)
	integer :: pui_F_n = 100      ! На сколько частей мы разбиваем первообразную для розыгрыша PUI 
	! т.е. мы будем разыгрывать ksi от 0 до 1 и брать значения скорости при этой ksi. 
	! интеграл посчитан для dksi = 1.0/pui_F_n а в значениях между придётся линейно интерполировать
	real(8), allocatable :: F_integr_pui(:, :)           ! (pui_F_n, :) первообразная для розыгрыша
	real(8), allocatable :: nu_integr_pui(:, :)           ! (pui_F_n, :)   частота перезарядки
	real(8), allocatable :: Mz_integr_pui(:, :)           ! (pui_F_n, :)   источник импульса
	real(8), allocatable :: E_integr_pui(:, :)           ! (pui_F_n, :)	   источник энергии


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
			if(j == 0) then
				pui_num_tetr(i) = 0
				CYCLE
			end if

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
		pui_num_tetr_2 = 0

		N = 0
		do i = 1, size(int2_all_tetraendron(1, :))
			j = int2_all_tetraendron_point(1, i)
			if(j == 0) then
				pui_num_tetr(i) = 0
				CYCLE
			end if

			if(int2_Cell_par2(1, j) <= 2) then
				N = N + 1
				pui_num_tetr_2(N) = i
			end if
		end do

		allocate(h0_pui(pui_h0_n))

	end subroutine PUI_Set

	subroutine PUI_f_Set()
		! Создаём вспомогательные массивы для функции распределения
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

		allocate(n_pui(N2))
		allocate(T_pui(N2))
		allocate(f_pui_num(N2))
		allocate(f_pui_cut(N2))
		
		f_pui_cut = pui_nW
		n_pui = 0.0
		T_pui = 0.0

		k = 1
		do i = 1, size(int2_Cell_par2(1, :))
			if(int2_Cell_par2(1, i) <= 2) then
				f_pui_num(k) = i
				k = k + 1
			end if
		end do

	end subroutine PUI_f_Set

	subroutine PUI_f_Set2()
		! Создаём массивы для самой функции распределения
		integer :: N, i, k, j, N2

		N2 = 0
		k = 1
		do i = 1, size(int2_Cell_par2(1, :))
			if(int2_Cell_par2(1, i) <= 2) then
				N2 = N2 + 1
			end if
		end do

		allocate(f_pui(pui_nW, N2))
		allocate(f_pui2(pui_nW, N2))
		f_pui = 0.0
		f_pui2 = 0.0
	end subroutine PUI_f_Set2

	subroutine PUI_n_T_culc()
		integer :: n2, N, i, n3
		real(8) :: w, S, S2, rho

		N = size(f_pui(1, :))

		! !$omp parallel
	 	! !$omp do private(i, w, S, S2, n3, rho)
		do n2 = 1, N
	        111  continue
			S = 0.0
			S2 = 0.0
			do i = 1, pui_nW
				w = (i-0.5) * pui_wR/pui_nW
				S = S + f_pui(i, n2) * 4 * par_pi_8 * w**2 * (pui_wR/pui_nW)
				S2 = S2 + f_pui(i, n2) * 4 * par_pi_8 * w**4 * (pui_wR/pui_nW)
			end do
			S2 = S2/(S * 3.0)
			n_pui(n2) = S
			T_pui(n2) = S2

			n3 = f_pui_num(n2) 
			rho = int2_Cell_par(1, n3) - int2_Cell_par_2(1, n3)
			if(rho < 0.0)  then
				print*, int2_Cell_par(1, n3), int2_Cell_par_2(1, n3) 
				print*, S
				print*, "coord", int2_coord(:, n3)
				print*, "------------"
				STOP 
			end if

			if(S > 0.99 * rho) then
				f_pui(:, n2) = f_pui(:, n2) * (0.9899 * rho/S)
				GOTO 111
			end if

		end do
		! !$omp end do
		! !$omp end parallel

	end subroutine PUI_n_T_culc

	subroutine PUI_F_integr_Set()
		! Создаём массивы для интегрирования функции распределения для последующего розыгрыша
		! и сразу вычисляет их
		! Также считаем концентрацию PUI
		integer :: N, n2, i, step, j, k
		real(8) :: S, SS, w, S1, S2, ff, UH, u, the

		print*, "PUI_F_integr_Set  START"
		N = size(f_pui(1, :))
		allocate(F_integr_pui(pui_F_n, N))
		allocate(nu_integr_pui(pui_F_n, N))
		allocate(Mz_integr_pui(pui_F_n, N))
		allocate(E_integr_pui(pui_F_n, N))
		print*, "PUI_F_integr_Set  END"
	end subroutine PUI_F_integr_Set

	subroutine PUI_F_integr_Culc()
		! Создаём массивы для интегрирования функции распределения для последующего розыгрыша
		! и сразу вычисляет их
		! Также считаем концентрацию PUI
		integer :: N, n2, i, step, j, k
		real(8) :: S, SS, w, S1, S2, ff, UH, u, the

		print*, "PUI_F_integr_Culc  START"
		N = size(f_pui(1, :))

		step = 0
		! Посчитаем этот интеграл (см. документацию PUI \rho_w)
		!$omp parallel
	 	!$omp do private(SS, S, w, S1, S2, i, ff, UH, j, u, the, k)
		do n2 = 1, N
			!$omp critical
				step = step + 1
				if(mod(step, 500) == 0) then
					print*, step, "  from = ", N
				end if
			!$omp end critical
			SS = 0.0

			do i = 1, 10000
				w = (i-0.5) * pui_wR/10000
				ff = PUI_get_f(n2, w)
				SS = SS + ff * w**2 * (w + pui_h0_wc) * MK_sigma(w + pui_h0_wc) * (pui_wR/10000)
			end do

			S1 = 0.0
			w = 0.0
			!! Далее алгоритм посложнее, там шаг по w должен быть как можно меньше
			do i = 1, pui_F_n
				S2 = (i - 0.5) * 1.0/pui_F_n
				!print*, "S2 = ", S2
				do while (S1 < S2)
					w = w + (pui_wR/10000)
					ff = PUI_get_f(n2, w)
					S1 = S1 + ff * w**2 * (w + pui_h0_wc) * MK_sigma(w + pui_h0_wc) * (pui_wR/10000)/SS
				end do
				F_integr_pui(i, n2) = w
			end do


			! Считаем частоту перезарядки и источники импульса и энергии
			do i = 1, pui_F_n
				S = 0.0
				SS = 0.0
				S1 = 0.0
				UH = (i - 0.5) * pui_wR/pui_F_n
				do j = 1, 1000
					w = (j - 0.5) * pui_wR/1000
					ff = PUI_get_f(n2, w)
					do k = 1, 180
						the = par_pi_8 * k/180	
						u = sqrt(w**2 * sin(the)**2 + (w * cos(the) - UH)**2)
						u = max(u, 0.000001_8)
						S = S + 2.0 * par_pi_8 * ff * w**2 * u * MK_sigma(u) * sin(the) * (par_pi_8/180) * (pui_wR/1000)
						SS = SS + 2.0 * par_pi_8 * ff * w**2 * u * MK_sigma(u) * sin(the) * w * cos(the) * (par_pi_8/180) * (pui_wR/1000)
						S1 = S1 + par_pi_8 * ff * w**2 * u * MK_sigma(u) * sin(the) * w**2 * (par_pi_8/180) * (pui_wR/1000)   !TODO 
					end do
				end do
				nu_integr_pui(i, n2) = S
				Mz_integr_pui(i, n2) = SS
				E_integr_pui(i, n2) = S1
				!print*, UH, S, SS, SS/S, UH - SS/S
				!pause
			end do


		end do
		!$omp end do
		!$omp end parallel
		print*, "PUI_F_integr_Culc  END"
	end subroutine PUI_F_integr_Culc

	subroutine PUI_Culc_h0()
		implicit none
		integer :: i, j, k
		real(8) :: the, UH, w, u, S

		!$omp parallel
	 	!$omp do private(S, UH, u, the, w, i, j)
		do k = 1, pui_h0_n
			S = 0.0
			UH = (k - 0.5) * pui_wR/pui_h0_n
			do i = 1, 1000
				w = (i - 0.5) * pui_wR / 1000
				do j = 1, 180
					the = j * par_pi_8 / 180.0
					u = sqrt((w * sin(the))**2 + (w * cos(the) - UH)**2)
					S = max(S, u * MK_sigma(u)/((w + pui_h0_wc) * MK_sigma(w + pui_h0_wc)) )
				end do
			end do
			h0_pui(k) = S
		end do
		!$omp end do
		!$omp end parallel

		open(1, file = "h0_PUI.txt")
		do k = 1, pui_h0_n
			UH = (k - 0.5) * pui_wR/pui_h0_n
			write(1, *) UH, h0_pui(k)
		end do
		close(1)

		open(1, file = "n0_PUI_2d_UH=40.txt")
		UH = 30.0
		k = min(INT(UH/pui_wR * pui_h0_n) + 1, pui_h0_n)
		do i = 1, 1000
			w = (i - 0.5) * pui_wR / 1000
			do j = 1, 180
				the = j * par_pi_8 / 180.0
				u = sqrt((w * sin(the))**2 + (w * cos(the) - UH)**2)
				S = u * MK_sigma(u)/((w + pui_h0_wc) * MK_sigma(w + pui_h0_wc))/h0_pui(k)
				write(1, *) w, the, S
			end do
		end do
		close(1)

	end subroutine PUI_Culc_h0

	subroutine Culc_f_pui()
		USE Surface_setting
		integer :: i, k, num, yzel, tetraedron, iw, numw, num_do, i1, num2, cell_n, step
		real(8) :: r(3), PAR(9), dt, Sm, w0, q1, rho, rho0, w, qInt, Sp, rho_do
		real(8) :: aa(3), bb(3), cc(3)
		real(8) :: s, cospsi, C, A, B
		real(8) :: PAR_k(5), normal(3)
		real(8) :: f0_pui(pui_nW)
		real(8) :: mas_w(pui_nW)
		real(8) :: mas_w0(pui_nW)
		real(8) :: mas_Sm(pui_nW)
		real(8) :: mas_Sm2(pui_nW)
		logical :: find_n

		pui_Sm = pui_Sm * par_n_p_LISM
		pui_Sp = pui_Sp * par_n_p_LISM


		print*, "START Culc_f_pui"

		dt = 0.0001 !0.0005

		! Находим функцию распределения для ячеек перед ударной волной
		print*, "Do TS"
		!$omp parallel
	 	!$omp do private(step, k, num, yzel, tetraedron, iw, numw, num_do, i1, num2, r, PAR, Sm, w0, q1, rho, rho0, w, qInt, Sp, rho_do, s, cospsi, C, A, B, PAR_k, normal, f0_pui, mas_w, mas_w0, mas_Sm, mas_Sm2, find_n)
		do i = 1, size(f_pui_num)
			if(mod(i, 5000) == 0) print*, "i = ", i, " from", size(f_pui_num)
			k = f_pui_num(i)      ! Номер узла сетки интерполяции, в которой считаем PUI
			if(int2_Cell_par2(1, k) == 2) CYCLE

			do iw = 1, pui_nW
				mas_w0(iw) = ((iw - 0.5) * pui_wR / pui_nW)
			end do

			rho0 = int2_Cell_par(1, k)
			f0_pui = 0.0
			mas_Sm = 0.0

			r = int2_coord(:, k)  ! Координаты этого узла
			mas_w = mas_w0

			num = 3
			qInt = 0.0            ! Интеграл от источника массы при ионизации

			step = 0
			! Бежим до Солнца
			do while (.TRUE.)
				step = step + 1
				call Int2_Get_par_fast(r(1), r(2), r(3), num, PAR, PAR_k = PAR_k)
				q1 = PAR_k(1)
				rho = PAR(1)
				qInt = qInt + dt * q1/rho
				r = r - PAR(2:4) * dt
				tetraedron = pui_num_tetr(num) ! Номер тетраэдра в массиве источников

				do iw = 1, pui_nW
					numw = min(INT(mas_w(iw)/pui_wR * pui_nW) + 1, pui_nW)
					if(mas_w(iw) < pui_wR .and. mas_w(iw) > 0) then
						mas_Sm(iw) = mas_Sm(iw) + pui_Sm(numw, tetraedron) * dt
						f0_pui(iw) = f0_pui(iw) + pui_Sp(numw, tetraedron) * dt * exp(-mas_Sm(iw))  ! Это S+, просто сразу накапливаем в функцию распределения
					end if
					mas_w(iw) = mas_w0(iw) / ( (rho0/rho)**(1.0/3.0) * exp(-1.0/3.0 * qInt) )
				end do

				if(norm2(r) <= par_R0) EXIT

				if(step > 100000) then
					print*, "Step = ", step, " r = ", r
				end if
			end do

			f_pui(:, i) = f0_pui(:)
		end do
		!$omp end do
		!$omp end parallel


		! Находим функцию распределения для ячеек за ударной волной
		print*, "Posle TS"
		!$omp parallel
	 	!$omp do private(aa, bb, cc, k, num, yzel, tetraedron, iw, numw, num_do, i1, num2, r, PAR, Sm, w0, q1, rho, rho0, w, qInt, Sp, rho_do, s, cospsi, C, A, B, PAR_k, normal, f0_pui, mas_w, mas_w0, mas_Sm, mas_Sm2, find_n)
		do i = 1, size(f_pui_num)

			do iw = 1, pui_nW
				mas_w0(iw) = ((iw - 0.5) * pui_wR / pui_nW)
			end do

			if(mod(i, 500) == 0) print*, "i = ", i, " from", size(f_pui_num)
			k = f_pui_num(i)      ! Номер узла сетки интерполяции, в которой считаем PUI
			if(int2_Cell_par2(1, k) == 1) CYCLE   ! Пропускаем ячейки до TS
			rho0 = int2_Cell_par(1, k)
			rho = rho0
			f0_pui = 0.0
			mas_Sm = 0.0

			r = int2_coord(:, k)  ! Координаты этого узла
			mas_w = mas_w0

			num = 3
			qInt = 0.0            ! Интеграл от источника массы при ионизации

			!print*, r
			!pause

			! Бежим до TS
			do while (.TRUE.)
				num_do = num
				call Int2_Get_par_fast(r(1), r(2), r(3), num, PAR, PAR_k = PAR_k)
				q1 = PAR_k(1)
				rho_do = rho
				rho = PAR(1)
				yzel = int2_all_tetraendron_point(1, num)
				if(int2_Cell_par2(1, yzel) == 1) EXIT   ! Если попали в ячейку из области 1 - область до TS

				if(int2_Cell_par2(1, yzel) == 3) then   ! Если попали в ячейку из области 3 - область за HP
					!print*, "Popal v zonu 3"
					r = r * 0.999
					CYCLE
				end if

				qInt = qInt + dt * q1/rho
				r = r - PAR(2:4) * dt
				!print*, r
				!pause
				tetraedron = pui_num_tetr(num) ! Номер тетраэдра в массиве источников

				do iw = 1, pui_nW
					numw = min(INT(mas_w(iw)/pui_wR * pui_nW) + 1, pui_nW)
					if(mas_w(iw) < pui_wR .and. mas_w(iw) > 0) then
						mas_Sm(iw) = mas_Sm(iw) + pui_Sm(numw, tetraedron) * dt
						f0_pui(iw) = f0_pui(iw) + pui_Sp(numw, tetraedron) * dt * exp(-mas_Sm(iw))  ! Это S+, просто сразу накапливаем в функцию распределения
					end if
					mas_w(iw) = mas_w0(iw) / ( (rho0/rho)**(1.0/3.0) * exp(-1.0/3.0 * qInt) )
				end do

			end do

			!f_pui(:, i) = f0_pui(:)
			! Также нужно сохранить mas_Sm, для умножения посчитанной дальше функции распределения на него
			!f0_pui = 0.0
			!mas_Sm2 = 0.0

			! Теперь надо получить нормаль к поверхности TS
			! num - это номер тетраэдра в котором сейчас находимся 

			! cell_n = INT(num/6) + 1  !! Номер ячейки двойственной сетки, в которой сейчас находимся   Проверить это!
			! aa = int2_coord(:, int2_all_Cell(7, cell_n)) - int2_coord(:, int2_all_Cell(2, cell_n))
        	! bb = int2_coord(:, int2_all_Cell(6, cell_n)) - int2_coord(:, int2_all_Cell(3, cell_n))
			! Проблема в том, для ячеек из группы B по другому располагаются узлы, и такой набор задают не ту грань
			! normal(1) = aa(2) * bb(3) - aa(3) * bb(2)
			! normal(2) = aa(3) * bb(1) - aa(1) * bb(3)
			! normal(3) = aa(1) * bb(2) - aa(2) * bb(1)

			! aa(1) = norm2(normal)  ! S = S/2
			! normal = normal/aa(1)

			! if(DOT_PRODUCT(r, normal) < 0) then
			! 	print*, r
			! 	print*, normal
			! 	print*, int2_coord(:, int2_all_Cell(2, cell_n))
			! 	print*, int2_coord(:, int2_all_Cell(3, cell_n))
			! 	print*, int2_coord(:, int2_all_Cell(7, cell_n))
			! 	print*, int2_coord(:, int2_all_Cell(6, cell_n))
			! 	print*, "ERROR normal 424"
			! 	pause
			! end if

			find_n = .FALSE.
			do i1 = 1, 4
				num2 = int2_gran_sosed(int2_all_tetraendron(i1, num))
				if(num2 == num_do) then
					normal = int2_plane_tetraendron(1:3, i1, num)
					find_n = .TRUE.
					EXIT
				end if
			end do

			if(find_n == .FALSE.) then
				aa = r
				bb = r
				cc = r
				aa(1) = aa(1) + 0.1 * par_R0
				aa(3) = aa(3) + 0.12 * par_R0
				bb(2) = bb(2) + 0.1 * par_R0
				bb(3) = bb(3) + 0.1 * par_R0
				aa = aa / (norm2(aa) / Surf_Get_TS(aa(1), aa(2), aa(3)))
				bb = bb / (norm2(bb) / Surf_Get_TS(bb(1), bb(2), bb(3)))
				cc = cc / (norm2(cc) / Surf_Get_TS(cc(1), cc(2), cc(3)))
				bb = cc - bb
				aa = cc - aa
				if(norm2(aa) > 0.0001 .and. norm2(bb) > 0.0001) then
					normal(1) = aa(2) * bb(3) - aa(3) * bb(2)
					normal(2) = aa(3) * bb(1) - aa(1) * bb(3)
					normal(3) = aa(1) * bb(2) - aa(2) * bb(1)
					aa(1) = norm2(normal)  ! S = S/2
					normal = normal/aa(1)
					if(DOT_PRODUCT(r, normal) < 0) normal = -normal
					find_n = .TRUE.
				end if
			end if

			if(find_n /= .TRUE.) then
				print*, "Ne nashol normal ", r
				normal = r
				normal = normal/norm2(r)
			end if

			s = rho_do/rho
			cospsi = DOT_PRODUCT(PAR(6:8), normal)/(norm2(normal) * norm2(PAR(6:8)))
			A = sqrt(cospsi**2 + s**2 * (1.0 - cospsi**2))
			B = s**2 / A**2
			C = (2.0 * A + B)/3.0

			rho0 = rho
			qInt = 0.0
			mas_w0 = mas_w/sqrt(C)
			mas_w = mas_w0

			! Бежим до Солнца
			do while (.TRUE.)
				call Int2_Get_par_fast(r(1), r(2), r(3), num, PAR, PAR_k = PAR_k)
				q1 = PAR_k(1)
				rho = PAR(1)
				qInt = qInt + dt * q1/rho
				r = r - PAR(2:4) * dt
				!print*, r
				!pause
				tetraedron = pui_num_tetr(num) ! Номер тетраэдра в массиве источников

				do iw = 1, pui_nW
					numw = min(INT(mas_w(iw)/pui_wR * pui_nW) + 1, pui_nW)
					if(mas_w(iw) < pui_wR .and. mas_w(iw) > 0) then
						mas_Sm(iw) = mas_Sm(iw) + pui_Sm(numw, tetraedron) * dt
						f0_pui(iw) = f0_pui(iw) + pui_Sp(numw, tetraedron) * dt * exp(-mas_Sm(iw)) * s/C**(1.5)  ! Это S+, просто сразу накапливаем в функцию распределения
					end if
					mas_w(iw) = mas_w0(iw) / ( (rho0/rho)**(1.0/3.0) * exp(-1.0/3.0 * qInt) )
				end do

				if(norm2(r) <= par_R0) EXIT
			end do

			f_pui(:, i) = f0_pui(:)

		end do
		!$omp end do
		!$omp end parallel

		print*, "END Culc_f_pui"
	end subroutine Culc_f_pui

	real(8) pure function PUI_get_f(n, w)
		!real(8) function PUI_get_f(n, w)
		! Получить значение функции распределения PUI с номером n, со скоростью w (используется линейная интерполяция)
		real(8), intent (in) :: w
		integer, intent (in) :: n
		integer :: k1, k2
		real(8) :: x1, x2,  f1, f2
		k1 = max(min(INT(w/pui_wR * pui_nW + 0.5), pui_nW), 1)
		k2 = min(k1 + 1, pui_nW)
		x1 = (k1 - 0.5) * pui_wR/pui_nW
		x2 = (k2 - 0.5) * pui_wR/pui_nW
		f1 = f_pui(k1, n)
		f2 = f_pui(k2, n)
		!PUI_get_f = (f1 - f2)/(x1 - x2) * x + (f1 * x2 - f2 * x1)/(x2 - x1)

		! if(w > 5.0 .and.(w < x1 .or. w > x2)) then
		! 	print*, "ERROR PUI_get_f 390", w, x1, x2
		! end if

		if(w < x1 .or. dabs(x1 - x2) <= 0.0000001) then
			PUI_get_f = f1
		else if (w > x2) then
			PUI_get_f = f2
		else
			PUI_get_f = f1 * (w - x2)/(x1 - x2) + f2 * (w - x1)/(x2 - x1)
		end if
    end function PUI_get_f

	real(8) pure function PUI_get_h0(w)
		! Получить значение функции h0 
		real(8), intent (in) :: w

		integer :: k1, k2
		real(8) :: x1, x2,  f1, f2
		k1 = max(min(INT(w/pui_wR * pui_h0_n + 0.5), pui_h0_n), 1)
		k2 = min(k1 + 1, pui_h0_n)
		x1 = (k1 - 0.5) * 1.0/pui_h0_n
		x2 = (k2 - 0.5) * 1.0/pui_h0_n
		f1 = h0_pui(k1)
		f2 = h0_pui(k2)

		if(dabs(x1 - x2) <= 0.0000001) then
			PUI_get_h0 = f1
		else if(w < x1) then
			PUI_get_h0 = f1   ! * (ksi)/(x1)   
			!TODO
		else if (w > x2) then
			PUI_get_h0 = f2
		else
			PUI_get_h0 = f1 * (w - x2)/(x1 - x2) + f2 * (w - x1)/(x2 - x1)
		end if
    end function PUI_get_h0

	real(8) pure function PUI_get_F_integer(ksi, n)
		! Получить значение функции F_integr_pui  
		real(8), intent (in) :: ksi
		integer, intent (in) :: n

		integer :: k1, k2
		real(8) :: x1, x2,  f1, f2
		k1 = max(min(INT(ksi * pui_F_n + 0.5), pui_F_n), 1)
		k2 = min(k1 + 1, pui_F_n)
		x1 = (k1 - 0.5) * 1.0/pui_F_n
		x2 = (k2 - 0.5) * 1.0/pui_F_n
		f1 = F_integr_pui(k1, n)
		f2 = F_integr_pui(k2, n)

		if(dabs(x1 - x2) <= 0.0000001) then
			PUI_get_F_integer = f1
		else if(ksi < x1) then
			PUI_get_F_integer = f1 * (ksi)/(x1)
		else if (ksi > x2) then
			PUI_get_F_integer = f2
		else
			PUI_get_F_integer = f1 * (ksi - x2)/(x1 - x2) + f2 * (ksi - x1)/(x2 - x1)
		end if
    end function PUI_get_F_integer

	real(8) pure function PUI_get_nu_integr(n, w)
		!real(8) function PUI_get_f(n, w)
		! Получить значение функции распределения PUI с номером n, со скоростью w (используется линейная интерполяция)
		real(8), intent (in) :: w
		integer, intent (in) :: n
		integer :: k1, k2
		real(8) :: x1, x2,  f1, f2
		k1 = max(min(INT(w/pui_wR * pui_F_n + 0.5), pui_F_n), 1)
		k2 = min(k1 + 1, pui_F_n)
		x1 = (k1 - 0.5) * pui_wR/pui_F_n
		x2 = (k2 - 0.5) * pui_wR/pui_F_n
		f1 = nu_integr_pui(k1, n)
		f2 = nu_integr_pui(k2, n)
		!PUI_get_f = (f1 - f2)/(x1 - x2) * x + (f1 * x2 - f2 * x1)/(x2 - x1)

		! if(w > 5.0 .and.(w < x1 .or. w > x2)) then
		! 	print*, "ERROR PUI_get_f 390", w, x1, x2
		! end if

		if(w < x1 .or. dabs(x1 - x2) <= 0.0000001) then
			PUI_get_nu_integr = f1
		else if (w > x2) then
			PUI_get_nu_integr = f2
		else
			PUI_get_nu_integr = f1 * (w - x2)/(x1 - x2) + f2 * (w - x1)/(x2 - x1)
		end if
    end function PUI_get_nu_integr

	real(8) pure function PUI_get_Mz_integr(n, w)
		!real(8) function PUI_get_f(n, w)
		! Получить значение функции распределения PUI с номером n, со скоростью w (используется линейная интерполяция)
		real(8), intent (in) :: w
		integer, intent (in) :: n
		integer :: k1, k2
		real(8) :: x1, x2,  f1, f2
		k1 = max(min(INT(w/pui_wR * pui_F_n + 0.5), pui_F_n), 1)
		k2 = min(k1 + 1, pui_F_n)
		x1 = (k1 - 0.5) * pui_wR/pui_F_n
		x2 = (k2 - 0.5) * pui_wR/pui_F_n
		f1 = Mz_integr_pui(k1, n)
		f2 = Mz_integr_pui(k2, n)

		if(w < x1 .or. dabs(x1 - x2) <= 0.0000001) then
			PUI_get_Mz_integr = f1
		else if (w > x2) then
			PUI_get_Mz_integr = f2
		else
			PUI_get_Mz_integr = f1 * (w - x2)/(x1 - x2) + f2 * (w - x1)/(x2 - x1)
		end if
    end function PUI_get_Mz_integr

	real(8) pure function PUI_get_E_integr(n, w)
		!real(8) function PUI_get_f(n, w)
		! Получить значение функции распределения PUI с номером n, со скоростью w (используется линейная интерполяция)
		real(8), intent (in) :: w
		integer, intent (in) :: n
		integer :: k1, k2
		real(8) :: x1, x2,  f1, f2
		k1 = max(min(INT(w/pui_wR * pui_F_n + 0.5), pui_F_n), 1)
		k2 = min(k1 + 1, pui_F_n)
		x1 = (k1 - 0.5) * pui_wR/pui_F_n
		x2 = (k2 - 0.5) * pui_wR/pui_F_n
		f1 = E_integr_pui(k1, n)
		f2 = E_integr_pui(k2, n)
		!PUI_get_f = (f1 - f2)/(x1 - x2) * x + (f1 * x2 - f2 * x1)/(x2 - x1)

		! if(w > 5.0 .and.(w < x1 .or. w > x2)) then
		! 	print*, "ERROR PUI_get_f 390", w, x1, x2
		! end if

		if(w < x1 .or. dabs(x1 - x2) <= 0.0000001) then
			PUI_get_E_integr = f1
		else if (w > x2) then
			PUI_get_E_integr = f2
		else
			PUI_get_E_integr = f1 * (w - x2)/(x1 - x2) + f2 * (w - x1)/(x2 - x1)
		end if
    end function PUI_get_E_integr

	subroutine PUI_Add(cell, wr, nu_ex, mu, time)
		! wr - скорость атома в СК, связанной со средней скоростью плазмы (модуль этой скорости)
		integer, intent(in) :: cell
		real(8), intent(in) :: wr, nu_ex, mu, time
		integer i, j

		j = pui_num_tetr(cell)
		if(j > 0) then
			i = min(INT(wr/pui_wR * pui_nW) + 1, pui_nW)

			call omp_set_lock(pui_lock(j))
			pui_Sm(i, j) = pui_Sm(i, j) + mu * time
			pui_Sp(i, j) = pui_Sp(i, j) + mu * time * nu_ex
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
		print*, "Start PUI_calc_Sm"
	 	!$omp parallel

	 	!$omp do private(pui_Sm2, ij, j, k, ff, Vh, the, d, w)
		do i = 1, size(pui_Sm(1, :))

			!$omp critical
			num_all = num_all + 1
			if(mod(num_all, 50000) == 0) then
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

			pui_Sm(:, i) = pui_Sm2/par_Kn !/2.0  !TODO НУЖНО БУДЕТ УБРАТЬ ДЕЛЕНИЕ НА ДВА В СЛЕДУЮЩИЙ РАЗ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		end do
		!$omp end do
		!$omp end parallel

		print*, "End PUI_calc_Sm"

	end subroutine PUI_calc_Sm

	subroutine Cut_f_pui()
		implicit none
		integer(4) :: i, j, n_yzel, zone
		real(8) :: SS, w, SS2, S, k_norm
		real(8) :: rho_Th, p_Th, n_he_, cp, n_sw, p_sw, n_pui_, T_pui_

		do i = 1, size(f_pui(1, :)) !! Нормируем функцию распределения Пикапов, если она слишком большая
			S = 0.0 ! Концентрация пикапов
			SS = 0.0 ! Температура пикапов пикапов

			do j = 1, pui_nW
				w = (j-0.5) * pui_wR/pui_nW
				S = S + f_pui(j, i) * 4 * par_pi_8 * w**2
			end do

			do j = 1, pui_nW
				w = (j-0.5) * pui_wR/pui_nW
				SS = SS + f_pui(j, i) * 4 * par_pi_8 * w**4 /(3.0 * S)
			end do

			n_sw = int2_Cell_par(1, n_yzel)
			n_he_ = int2_Cell_par_2(1, n_yzel)

			if(S > n_sw - n_he_) k_norm = 0.999 * (n_sw - n_he_)/S

			f_pui(:, i) = f_pui(:, i) * k_norm
		end do

		f_pui2 = f_pui


		do i = 1, size(f_pui(1, :))

			S = 0.0 ! Концентрация пикапов
			SS = 0.0 ! Температура пикапов пикапов

			do j = 1, pui_nW
				w = (j-0.5) * pui_wR/pui_nW
				S = S + f_pui(j, i) * 4 * par_pi_8 * w**2
			end do

			do j = 1, pui_nW
				w = (j-0.5) * pui_wR/pui_nW
				SS = SS + f_pui(j, i) * 4 * par_pi_8 * w**4 /(3.0 * S)
			end do

			n_yzel = f_pui_num(i)  ! Номер узла в интерполяционной сетке

			n_sw = int2_Cell_par(1, n_yzel)
			p_sw = int2_Cell_par(5, n_yzel)
			n_he_ = int2_Cell_par_2(1, n_yzel)
			zone = int2_Cell_par2(1, n_yzel)
			call Sootnosheniya(n_sw, p_sw, n_he_, S, SS, zone, rho_Th = rho_Th, p_Th = p_Th)
			cp = p_Th/rho_Th


			! do j = 1, pui_nW
			! 	w = (j-0.5) * pui_wR/pui_nW
			! 	SS2 = SS2 + f_pui(j, i) * 4 * par_pi_8 * w**4 /(3.0 * S)
			! 	if(dabs(SS2 - SS)/SS * 100.0 < 0.01) then
			! 		f_pui_cut(i) = j
			! 		EXIT
			! 	end if
			! end do

		end do

	end subroutine Cut_f_pui

	subroutine PUI_print(num, x, y, z)
		implicit none
		integer, intent(in) :: num
		real(8), intent(in) :: x, y, z
		real(8) :: w, S, SS, ff, UH
		character(len=5) :: name
		integer(4) :: n1, i, n2, number, zone
		real(8) :: MAS_PUI(2)
		real(8) :: PAR_MOMENT(par_n_moment, par_n_sort)
		real(8) :: PAR(9)     ! Выходные параметры
		real(8) :: cp, n_pui_, T_pui_, ro_He, rho_Th, p_Th

		n1 = 3
		number = 3
		

		call Int2_Get_tetraedron_inner(x, y, z, n1)
		call Int2_Get_par_fast(x, y, z, number, PAR, PAR_MOMENT = PAR_MOMENT, MAS_PUI = MAS_PUI, rho_He = ro_He)
		



		write(unit=name,fmt='(i5.5)') num

		print*, "Nomer = ", num, " tetraedr = ", n1, " cell number = ", int2_all_tetraendron_point(1, n1), &
			"   nomer f_pui = ", f_pui_num2(int2_all_tetraendron_point(1, n1)), "area = ", int2_Cell_par2(1, f_pui_num(f_pui_num2(int2_all_tetraendron_point(1, n1))))
		zone = int2_Cell_par2(1, f_pui_num(f_pui_num2(int2_all_tetraendron_point(1, n1))))
		! open(1, file = "S+_" // name // ".txt")
		! n2 = pui_num_tetr(n1)
		! if(n2 > 0) then
		! 	do i = 1, pui_nW
		! 		w = (i-0.5) * pui_wR/pui_nW
		! 		write(1,*) w, pui_Sp(i, n2)
		! 	end do
		! end if

		! close(1)

		! open(2, file = "S-_" // name // ".txt")
		! if(n2 > 0) then
		! 	do i = 1, pui_nW
		! 		w = (i-0.5) * pui_wR/pui_nW
		! 		write(2,*) w, pui_Sm(i, n2)
		! 	end do
		! end if

		! close(2)

		n2 = f_pui_num2(int2_all_tetraendron_point(1, n1))
		n_pui_ = n_pui(n2)
		T_pui_ = T_pui(n2)
		print*, n1, n2

		call Sootnosheniya(PAR(1), PAR(5), ro_He, n_pui_, T_pui_, zone, rho_Th = rho_Th, p_Th = p_Th)
		if(rho_Th > 0.0) then
			cp = sqrt(p_Th/rho_Th)
		else
			cp = 1
		end if

		open(3, file = "f_PUI_" // name // ".txt")

		write(3,*) "TITLE = 'HP'  VARIABLES = 'w', 'pui', 'th', 'pui + th'"
		write(3,*) ", ZONE T= 'HP'"

		if(n2 > 0) then
			do i = 1, pui_nW
				w = (i-0.5) * pui_wR/pui_nW
				ff = exp(-w**2/cp**2)/(cp**3 * par_pi_8**(1.5))
				write(3,*) w, f_pui(i, n2), ff, ff + f_pui(i, n2)
			end do

		end if

		close(3)

		S = 0.0
		!? Печатаем концентрацию пикапов (в виде графика от скорости, чтобы посмотреть, на каких скоростях накапливается максимум)
		! open(4, file = "n_pui_" // name // ".txt")
		! if(n2 > 0) then
		! 	do i = 1, pui_nW
		! 		w = (i-0.5) * pui_wR/pui_nW
		! 		S = S + f_pui(i, n2) * 4 * par_pi_8 * w**2
		! 		write(4,*) w, S
		! 	end do
		! end if
		! close(4)

		! SS = 0.0
		! !? Также печаетаем температуру пикапов
		! open(5, file = "T_pui_" // name // ".txt")
		! if(n2 > 0) then
		! 	do i = 1, pui_nW
		! 		w = (i-0.5) * pui_wR/pui_nW
		! 		SS = SS + f_pui(i, n2) * 4 * par_pi_8 * w**4 /(3.0 * S)
		! 		write(5,*) w, SS
		! 	end do
		! end if
		! close(5)

		!? Печатаем функцию розыгрыша (см. документацию "PUI")
		! open(5, file = "F(w)_" // name // ".txt")
		! SS = 0.0
		! if(n2 > 0) then
		! 	do i = 1, pui_F_n
		! 		w = (i-0.5) * pui_wR/pui_F_n
		! 		ff = PUI_get_f(n2, w)
		! 		SS = SS + ff * w**2 * (w + pui_h0_wc) * MK_sigma(w + pui_h0_wc)
		! 		!write(5,*) w, SS
		! 	end do

		! 	S = 0.0
		! 	do i = 1, pui_F_n
		! 		w = (i-0.5) * pui_wR/pui_F_n
		! 		ff = PUI_get_f(n2, w)
		! 		S = S + ff * w**2 * (w + pui_h0_wc) * MK_sigma(w + pui_h0_wc)/SS
		! 		if(ff < 0) then
		! 			print*, "ERROR  PUI_print  594 uikjhgbn"	
		! 			print*, n2, w, ff
		! 			print*, f_pui(:, n2)
		! 			pause
		! 		end if
		! 		write(5,*) S, w
		! 	end do
		! end if
		! close(5)

		!? Печатаем интеграллы (частоту и источники импульса и энергии, см. документацию)
		! open(5, file = "nu_pui_" // name // ".txt")
		! open(6, file = "Mz_pui_" // name // ".txt")
		! open(7, file = "E_pui_" // name // ".txt")
		! do i = 1, pui_F_n
		! 	UH = (i - 0.5) * pui_wR/pui_F_n
		! 	write(5,*) UH, nu_integr_pui(i, n2)
		! 	write(6,*) UH, Mz_integr_pui(i, n2)/nu_integr_pui(i, n2)
		! 	write(7,*) UH, E_integr_pui(i, n2)/nu_integr_pui(i, n2)
		! end do
		! close(5)
		! close(6)
		! close(7)

	end subroutine PUI_print

	subroutine PUI_print_in_cell(num, cell)
		implicit none
		integer, intent(in) :: num, cell
		
		real(8) :: w, S, SS, ff, UH
		character(len=5) :: name
		integer(4) :: n1, i, n2

		write(unit=name,fmt='(i5.5)') num


		open(3, file = "f_PUI_" // name // ".txt")
		do i = 1, pui_nW
			w = (i-0.5) * pui_wR/pui_nW
			write(3,*) w, f_pui(i, cell)
		end do
		close(3)

	end subroutine PUI_print_in_cell

	subroutine PUI_Save_bin(num)
		implicit none
		integer, intent(in) :: num
		character(len=5) :: name

		write(unit=name,fmt='(i5.5)') num

		open(1, file = par_NAME//"pui_save_" // name // ".bin", FORM = 'BINARY')

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

	subroutine PUI_Save_for_MK_bin(num)
		implicit none
		integer, intent(in) :: num
		character(len=5) :: name

		write(unit=name,fmt='(i5.5)') num

		open(1, file = par_NAME//"pui_save_for_MK_" // name // ".bin", FORM = 'BINARY')


		write(1) pui_F_n, pui_nW, pui_wR
		write(1) F_integr_pui
		write(1) nu_integr_pui
		write(1) Mz_integr_pui
		write(1) E_integr_pui
		write(1) n_pui
		write(1) T_pui


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

	end subroutine PUI_Save_for_MK_bin

	subroutine PUI_Read_for_MK_bin(num)
		implicit none
		integer, intent(in) :: num
		character(len=5) :: name
		logical :: exists
		integer :: n1, n2
		real(8) :: a

		write(unit=name,fmt='(i5.5)') num

		inquire(file= par_NAME//"pui_save_for_MK_" // name // ".bin", exist=exists)

		if (exists == .False.) then
			pause "ERROR net faila 171 PUI dfwefcwcw!!!"
			STOP "ERROR net faila 171 PUI wfecwvwv!!!"
		end if

		print*, "PUI_Read_for_MK_bin"
		open(1, file = par_NAME//"pui_save_for_MK_" // name // ".bin", FORM = 'BINARY', ACTION = "READ")

		read(1) n1, n2, a
		read(1) F_integr_pui
		read(1) nu_integr_pui
		read(1) Mz_integr_pui
		read(1) E_integr_pui
		read(1) n_pui
		read(1) T_pui

		close(1)
		print*, "END PUI_Read_for_MK_bin"

	end subroutine PUI_Read_for_MK_bin

	subroutine PUI_Save_f_bin(num)
		!? Сохранение функций распределения атомов (если они уже были посчитаны, чтобы второй раз не считать)
		!? Также сохраняет диапазон для урезанной функции распределения атомов
		implicit none
		integer, intent(in) :: num
		character(len=5) :: name

		write(unit=name,fmt='(i5.5)') num

		open(1, file = par_NAME//"pui_save_f_" // name // ".bin", FORM = 'BINARY')

		print*, pui_nW, pui_wR

		write(1) pui_nW
		write(1) pui_wR

		write(1) size(f_pui(1, :))
		write(1) f_pui

		write(1) f_pui_cut


		write(1) 1; 
		write(1) f_pui2
		
		write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
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
	end subroutine PUI_Save_f_bin

	subroutine PUI_Read_bin(num)
		!! Считывает не саму функцию PUI, а только S+ и S-
		implicit none
		integer, intent(in) :: num
		character(len=5) :: name
		integer(4) :: n1, n2, n3
		real(8) :: aa
		logical :: exists

		write(unit=name,fmt='(i5.5)') num

		inquire(file= par_NAME//"pui_save_" // name // ".bin", exist=exists)

		if (exists == .False.) then
			pause "ERROR net faila 174444 PUI!!!"
			STOP "ERROR net faila 174444 PUI!!!"
		end if

		open(1, file = par_NAME//"pui_save_" // name // ".bin", FORM = 'BINARY', ACTION = "READ")

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

		print*, "1 Sp = ", pui_Sp(10, 1), pui_Sp(10, 10), pui_Sp(10, 100), pui_Sp(10, 1000)
		print*, "2 Sp = ", pui_Sp(10, 20000), pui_Sp(10, 30000), pui_Sp(10, 40000), pui_Sp(10, 50000)

		close(1)

	end subroutine PUI_Read_bin

	subroutine PUI_Read_f_bin(num)
		implicit none
		integer, intent(in) :: num
		character(len=5) :: name
		integer(4) :: n1, n2, n3
		real(8) :: aa
		logical :: exists

		write(unit=name,fmt='(i5.5)') num

		inquire(file= par_NAME//"pui_save_f_" // name // ".bin", exist=exists)

		if (exists == .False.) then
			pause "ERROR net faila 171 PUI!!!"
			STOP "ERROR net faila 171 PUI!!!"
		end if

		open(1, file = par_NAME//"pui_save_f_" // name // ".bin", FORM = 'BINARY', ACTION = "READ")

		read(1) n2
		read(1) aa
		print*, "schitano = ", n2, aa

		if(pui_nW /= n2) then
			print*, n2, pui_nW
			STOP "PUI_Read_f_bin ERROR PUI pui_nW /= n2"
		end if
		if(pui_wR /= aa) then
			print*, pui_wR, aa
			STOP "PUI_Read_f_bin ERROR pui_wR /= aa"
		end if

		read(1) n1
		if(n1 /= size(f_pui(1, :))) then
			pause "ERROR n1 /= size(f_pui(1, :))"
			STOP "ERROR n1 /= size(f_pui(1, :))"
		end if
		read(1) f_pui

		read(1) f_pui_cut

		read(1) n1
		if(n1 == 1) read(1) f_pui2

		close(1)

	end subroutine PUI_Read_f_bin

	
end module PUI
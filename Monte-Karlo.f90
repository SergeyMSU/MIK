
	
module Monte_Karlo  
	! Модуль Монте-Карло, использует интерполяционную сетку тетраэдров
	USE GEO_PARAM
	USE STORAGE
	USE Interpolate2
	USE My_func
	USE ieee_arithmetic 
	USE PUI
	
	implicit none
	
	integer(4), parameter :: par_stek = 2000  ! Глубина стека (заранее выделяется память под него)
	logical, parameter :: MK_is_NaN = .False.    ! Нужны ли проверки на nan
	logical, parameter :: MK_Mu_stat = .False.    ! Нужно ли накапливать веса для статистики и весовых каэффициентов
	logical, parameter :: MK_photoionization = .True.    ! Нужна ли фотоионизация
	
	real(8), parameter :: MK_Mu_mult = 100.0_8  ! На что домножаем веса для избежания потери точности
	
	real(8) :: sqv_1, sqv_2, sqv_3, sqv_4, sqv   ! Потоки частиц через 
	real(8) :: MK_mu1, MK_mu2, MK_mu3, MK_mu4
	integer(4) :: MK_N                  ! Сколько всего частиц запущено (сумма по всем потокам)
	
	
	real(8) :: MK_R_zone(par_n_zone)   ! Радиусы зон
	real(8) :: MK_al_zone(par_m_zone)   ! Лучи зон
	real(8) :: MK_SINKR(par_m_zone + 1)   ! Критические синусы для каждой зоны по углу
	real(8), allocatable :: MK_Mu(:, :, :)   ! Веса зон (par_n_zone + 1, par_m_zone + 1, сортов)
	real(8), allocatable :: MK_Mu_statistic(:, :, :)   ! Веса зон (par_n_zone + 1, par_m_zone + 1, сортов)
	! Для накапливания весов зон
	
	real(8) :: MK_gam_zone(par_n_zone)   ! Параметр гамма для зон
	real(8) :: MK_A0_, MK_A1_   ! Параметры для начального запуска
	
	
	
	
	real(8), allocatable :: M_K_particle(:, :, :)   ! Частицы (8, par_stek, число потоков)
	! (три координаты, три скорости, вес, радиус перегелия)
	integer(4), allocatable :: M_K_particle_2(:, :, :)  ! Частицы 
	! (в какой ячейке частица, сорт, зона назначения по r, зона назначения по углу)
	logical(4), allocatable :: M_K_particle_3(:, :, :, :)  ! Частицы (par_n_zone + 1, par_m_zone + 1, par_stek, число потоков)
	
	integer(4), allocatable :: sensor(:, :, :)  !(3, 2, : par_n_potok число потоков)  ! датчики случайных чисел 
	! Каждому потоку по два датчика
	
	integer(4), allocatable :: stek(:)   ! (: число потоков) Переменная чтения и записи в стек
	! Где стоит переменная, там что-то лежит, чтобы записать, нужно увеличить значение на 1
	
	real(8), allocatable :: M_K_Moment(:, :, :, :)  ! (19, par_n_sort, :, par_n_potok) То, что накапливаем в ячейках (по каждому сорту отдельно)
	!(rho, u, v, w, T, Iu, Iv, Iw, IT, Huu, Huv, Huw, Hvv, Hvw, Hww, Huuu, Hvvv, Hwww,  In)
	!(1  , 2, 3, 4, 5,  6,  7,  8,  9,  10,  11,  12,  13,  14,  15,   16,   17,   18,  19)
	
	
	contains
	
	
	subroutine M_K_start()
		USE OMP_LIB
		!$MPI include 'mpif.h'
		! Variables
		integer(4) :: potok, num, i, cell, to_i, to_j, j, pp, iter, step, k
		real(8) :: mu_(par_n_zone + 1), Wt_(par_n_zone + 1), Wp_(par_n_zone + 1), Wr_(par_n_zone + 1), X_(par_n_zone + 1)
		logical :: bb
		real(8) :: sin_, x, phi, y, z, ksi, Vx, Vy, Vz, r_peregel, no, ksi1, ksi2, ksi3, ksi4, ksi5
		real(8) :: ll, rr, Vphi, Vr, pui_w2, pui_w1
		real(8), allocatable :: vol_sr(:)                                    ! Для осреднения в узлах
		real(8) :: start_time, end_time
		integer mpi_process_Rank, mpi_size_Of_Cluster, mpi_ierror, mpi_rank
		real(8), allocatable :: buff(:, :, :)  ! ДЛЯ MPI
		real(8), allocatable :: M_K_Moment_print(:, :, :)
		!$MPI integer :: mpi_status(MPI_STATUS_SIZE)
		real(8) :: PAR(9), MAS_PUI(2)
		integer :: nk
		
		call M_K_Set()    ! Создали массивы
		call M_K_init()   ! Инициализируем веса и т.д.

		!! PUI
		if(par_PUI == .True.) then
			pui_Sm = 0.0   ! Нужно обнулить массивы, в которые будем накапливать всё.
			pui_Sp = 0.0
		end if
		
		! call PUI_f_Set()
		! call PUI_f_Set2()
		! call PUI_Set()                !! PUI
		! call PUI_F_integr_Set()       !! PUI
		! call PUI_Read_for_MK_bin(9)   !! PUI
		! call PUI_Culc_h0()            !! PUI

		! nk = 206852
		! call Int2_Get_par_fast(-6.133_8, -40.0_8, 9.688_8, nk, PAR, MAS_PUI = MAS_PUI)
		! print*, PAR
		! print*, "____________________"
		! print*, MAS_PUI

		! pause
		
		print*, "Vesa = ", MK_mu1, MK_mu2, MK_mu3, MK_mu4
		print*, "MK_k_multiply = ", MK_k_multiply
		end_time = 0.0
		start_time = 0.0
		
		! Запускаем каждый поток в параллельном цикле

		call omp_set_num_threads(par_n_potok)
		start_time = omp_get_wtime()
		step = 1
		
		mpi_size_Of_Cluster = par_n_claster
		mpi_rank = 0
		!$MPI call MPI_INIT(mpi_ierror)
		!$MPI call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size_Of_Cluster, mpi_ierror)
		!$MPI call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, mpi_ierror)
		if(mpi_size_Of_Cluster /= par_n_claster) then
			print*, "MPI error par_n_claster  789tyuighjbnmp;edzxcvrertgd  ", mpi_size_Of_Cluster, par_n_claster
			!$MPI call MPI_Abort(MPI_COMM_WORLD, 102, mpi_ierror)
			STOP
		end if
		
		!call Get_sensor(mpi_rank) ! Считали датчики случайных чисел
		call Get_sensor_sdvig(1)
		
		!$omp parallel
		
		!$omp do private(potok, num, mu_, Wt_, Wp_, Wr_, X_, bb, i, Vx, Vy, Vz, cell, sin_, x, phi, y, z, ksi, r_peregel, no, to_i, to_j, ksi1, ksi2, ksi3, ksi4, ksi5, ll, rr, Vphi, Vr)
		do iter = 1, par_n_potok * par_n_parallel
		
			potok = (omp_get_thread_num() + 1) 
			
			!$omp critical
			print*, "start potok = ", potok, " iter = ", iter, "   step = ", step, "from = ", par_n_potok * par_n_parallel, &
			"  computer № ", mpi_rank
			step = step + 1
			!$omp end critical
			
			cell = 3
			! Запускаем частицы первого типа (с полусферы)
			do num = 1, MK_N1
				!if( mod(num, 100000) == 0) then
				!	print*, "num = ", num, "  from ", MK_N1, "  potok = ", potok, "  computer № ", mpi_rank
				!	!print*, sensor(:, 1, potok), sensor(:, 2, potok)
				!end if
				
				! sensor(:, 1, potok) = (/  730   ,    18493      ,    61 /)
				! sensor(:, 2, potok) = (/  21338   ,    11299    ,     833/)

				
				call MK_Init_Parametrs(potok, mu_, Wt_, Wp_, Wr_, X_, bb)
				
				
				do i = 1, par_n_zone + 1
					sin_ = sqrt(1.0 - (X_(i)**2))
					x = (par_Rmax) * X_(i)
					call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi)
					phi = 2.0 * par_pi_8 * ksi
					y = (par_Rmax) * sin_ * cos(phi)
					z = (par_Rmax) * sin_ * sin(phi)
					
					call Int2_Get_tetraedron( x, y, z, cell)
					call dekard_skorost(x, y, z, Wr_(i), Wp_(i), Wt_(i), Vx, Vy, Vz)
				
					if(cell < 1) then
						!$MPI call MPI_Abort(MPI_COMM_WORLD, 103, mpi_ierror)
						STOP "Error 0lhy976yihko  "
					end if
					
					
					if(i /= par_n_zone + 1 .or. bb == .True.) then
						! Добавляем частицу в стек
						stek(potok) = stek(potok) + 1
						M_K_particle(1:7, stek(potok), potok) = (/ x, y, z, Vx, Vy, Vz, mu_(i) * MK_mu1 * MK_Mu_mult /)
						M_K_particle_2(1, stek(potok), potok) = cell       ! В какой ячейке находится
						M_K_particle_2(2, stek(potok), potok) = int2_Cell_par2(1, int2_all_tetraendron_point(1, cell)) ! Сорт
						call MK_Distination( M_K_particle(1:3, stek(potok), potok), M_K_particle(4:6, stek(potok), potok),&
							to_i, to_j, r_peregel)
						M_K_particle(8, stek(potok), potok) = r_peregel
						M_K_particle_2(3, stek(potok), potok) = to_i  ! Зона назначения
						M_K_particle_2(4, stek(potok), potok) = to_j  ! Зона назначения
				
					end if
				end do

				if(par_PUI == .False.) then
				!if(.True.) then
					call M_K_Fly(potok)
				else
					call M_K_Fly_PUI(potok)    !! PUI
				end if
				
			end do
			
			! Запускаем частицы второго типа (вылет сверху)
			do num = 1, MK_N2
				call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi1)
				call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi2)
				call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi3)
				call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi4)
				call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi5)
				
				ll = par_Rleft
				rr = -0.001;
				x = ll + ksi1 * (rr - ll)
				phi = ksi2 * 2.0 * par_pi_8
				Vphi = cos(2.0 * par_pi_8 * ksi3) * sqrt(-log(1.0 - ksi4))
				Vx = par_Velosity_inf + sin(2.0 * par_pi_8 * ksi3) * sqrt(-log(1.0 - ksi4))
				Vr = -sqrt(-log(ksi5))
				y = par_Rup * cos(phi)
				z = par_Rup * sin(phi)
				
				call Int2_Get_tetraedron( x, y, z, cell)
				
				if(cell < 1) then
					print*, x, y, z, cell
					!$MPI call MPI_Abort(MPI_COMM_WORLD, 104, mpi_ierror)
					STOP "Error 0lhy976yihkoqwewqdqwd  "
				end if
				
				
				stek(potok) = stek(potok) + 1
				M_K_particle(1:7, stek(potok), potok) = (/ x, y, z, Vx, cos(phi) * Vr - sin(phi) * Vphi,&
					sin(phi) * Vr + cos(phi) * Vphi,  MK_mu2 * MK_Mu_mult /)
				M_K_particle_2(1, stek(potok), potok) = cell       ! В какой ячейке находится
				M_K_particle_2(2, stek(potok), potok) = int2_Cell_par2(1, int2_all_tetraendron_point(1, cell)) ! Сорт
				call MK_Distination( M_K_particle(1:3, stek(potok), potok), M_K_particle(4:6, stek(potok), potok),&
					to_i, to_j, r_peregel)
				M_K_particle(8, stek(potok), potok) = r_peregel
				M_K_particle_2(3, stek(potok), potok) = to_i  ! Зона назначения
				M_K_particle_2(4, stek(potok), potok) = to_j  ! Зона назначения

				if(par_PUI == .False.) then
				!if(.True.) then
					call M_K_Fly(potok)
				else
					call M_K_Fly_PUI(potok)    !! PUI
				end if

			end do
			
			! Запускаем частицы третьего типа (вылет сзади)
			do num = 1, MK_N3
				call MK_Velosity_initial2(potok, Vx, Vy, Vz)
				call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi1)
				call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi2)
				rr = sqrt(ksi1 * (par_Rup) * (par_Rup))
				phi = ksi2 * 2.0 * par_pi_8
				y = rr * cos(phi)
				z = rr * sin(phi)
				
				call Int2_Get_tetraedron(par_Rleft, y, z, cell)
				if(cell < 1) then
					!$MPI call MPI_Abort(MPI_COMM_WORLD, 105, mpi_ierror)
					STOP "Error 0lhy976yihkodfresdfre  "
				end if
				
				stek(potok) = stek(potok) + 1
				M_K_particle(1:7, stek(potok), potok) = (/ par_Rleft, y, z, Vx, Vy, Vz,  MK_mu3 * MK_Mu_mult /)
				M_K_particle_2(1, stek(potok), potok) = cell       ! В какой ячейке находится
				M_K_particle_2(2, stek(potok), potok) = int2_Cell_par2(1, int2_all_tetraendron_point(1, cell)) ! Сорт
				call MK_Distination( M_K_particle(1:3, stek(potok), potok), M_K_particle(4:6, stek(potok), potok),&
					to_i, to_j, r_peregel)
				M_K_particle(8, stek(potok), potok) = r_peregel
				M_K_particle_2(3, stek(potok), potok) = to_i  ! Зона назначения
				M_K_particle_2(4, stek(potok), potok) = to_j  ! Зона назначения

				if(par_PUI == .False.) then
				!if(.True.) then
					call M_K_Fly(potok)
				else
					call M_K_Fly_PUI(potok)    !! PUI
				end if
			end do
			
			! Запускаем частицы четвёрного типа (вылет спереди с части плоскости)
			do num = 1, MK_N4
				call MK_Velosity_initial(potok, Vx, Vy, Vz)
				call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi1)
				call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi2)
				rr = sqrt(ksi1 * ((par_Rup**2) - (par_Rmax**2)) + (par_Rmax**2))
				phi = ksi2 * 2.0 * par_pi_8
				y = rr * cos(phi)
				z = rr * sin(phi)
				
				call Int2_Get_tetraedron(-0.001_8, y, z, cell)
				if(cell < 1) then
					!$MPI call MPI_Abort(MPI_COMM_WORLD, 106, mpi_ierror)
					STOP "Error 0lhy976yihko133131312  "
				end if
				
				stek(potok) = stek(potok) + 1
				M_K_particle(1:7, stek(potok), potok) = (/ -0.001_8, y, z, Vx, Vy, Vz,  MK_mu4 * MK_Mu_mult /)
				M_K_particle_2(1, stek(potok), potok) = cell       ! В какой ячейке находится
				M_K_particle_2(2, stek(potok), potok) = int2_Cell_par2(1, int2_all_tetraendron_point(1, cell)) ! Сорт
				call MK_Distination( M_K_particle(1:3, stek(potok), potok), M_K_particle(4:6, stek(potok), potok),&
					to_i, to_j, r_peregel)
				M_K_particle(8, stek(potok), potok) = r_peregel
				M_K_particle_2(3, stek(potok), potok) = to_i  ! Зона назначения
				M_K_particle_2(4, stek(potok), potok) = to_j  ! Зона назначения

				if(par_PUI == .False.) then
				!if(.True.) then
					call M_K_Fly(potok)
				else
					call M_K_Fly_PUI(potok)    !! PUI
				end if
			end do
		end do
		!$omp end do

		!$omp end parallel

		print*, "END ALL particle"
		
		end_time = omp_get_wtime()
		
		if(mpi_rank == 0) print *, "Time work: ", (end_time-start_time)/60.0, "   in minutes"
		
		no = MK_Mu_mult * MK_N * par_n_claster
		M_K_Moment(:, :, :, :) = M_K_Moment(:, :, :, :) / no  ! Вынес сюда для избежания потери точности при сложении
		if(par_pogloshenie == .True.) pogloshenie(:, :, :) = pogloshenie(:, :, :) / no

		if(par_PUI == .True.) then
			pui_Sm(:, :) = sqv * pui_Sm(:, :) / no
			do i = 1, pui_nW
				pui_w1 = (i - 1) * pui_wR/pui_nW 
				pui_w2 = i * pui_wR/pui_nW 
				pui_Sp(i, :) = sqv * pui_Sp(i, :) / (no * 4.0 * par_pi_8 * (1.0/3.0) * (pui_w2**3 - pui_w1**3))
			end do
		end if
		
		do i = 2, par_n_potok
			M_K_Moment(:, :, :, 1) = M_K_Moment(:, :, :, 1) + M_K_Moment(:, :, :, i)
		end do
		
		
		! Печатаем результаты в файл для последующего суммирования
		open(3, file = "M-K_param_007.bin", FORM = 'BINARY')
		no = 1.0_8 * MK_N * par_n_claster
		allocate(M_K_Moment_print, mold = M_K_Moment(:, :, :, 1))
		M_K_Moment_print = M_K_Moment(:, :, :, 1)
		WRITE(3) no
		WRITE(3) M_K_Moment_print
		close(3)
		deallocate(M_K_Moment_print)
		no = MK_Mu_mult * MK_N * par_n_claster
		
		
		! Сложим все MPI потоки
		!$MPI if(mpi_rank  == 0) allocate(buff(par_n_moment, par_n_sort, size(int2_all_tetraendron(1, :))))
		!$MPI call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierror)
		
		!$MPI do i = 1, par_n_claster - 1
		!$MPI 	if(mpi_rank == i) call MPI_SEND(M_K_Moment(:, :, :, 1), size(buff), MPI_DOUBLE_PRECISION, 0, &
		!$MPI 		100, MPI_COMM_WORLD, mpi_ierror)
		!$MPI 	if(mpi_rank == 0) then
		!$MPI 		call MPI_RECV(buff, size(buff), MPI_DOUBLE_PRECISION, i, 100, MPI_COMM_WORLD, mpi_status, mpi_ierror)
		!$MPI 		M_K_Moment(:, :, :, 1) = M_K_Moment(:, :, :, 1) + buff
		!$MPI 	end if
		!$MPI end do
		
		!$MPI if(mpi_rank  == 0) deallocate(buff)
		!$MPI call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierror)
		!$MPI call MPI_FINALIZE(mpi_ierror)
		
		!$MPI if(mpi_rank  /= 0) then
		!$MPI print*, "Process  ", mpi_rank, "  zavershon"
		!$MPI STOP 
		!$MPI end if
		
		! Бежим по всем тетраэдрам и нормируем моменты
		do i = 1, size(M_K_Moment(1, 1, :, 1))
			
			if(int2_all_tetraendron_point(1, i) == 0) CYCLE
			
			if( int2_all_Volume(i) <= 0.00000001) then
				!print*, "Error  dfgdfgdg89346767809098742577"
				!print*, M_K_Moment(:, 4, i, 1)
				print*, int2_all_Volume(i)
				continue
			end if
			
			!no = MK_Mu_mult * MK_N * int2_all_Volume(i)
			no = int2_all_Volume(i)
			
			if(MK_is_NaN == .True. .and. ieee_is_nan(no)) then
					print*, "NaN lj098inbh5dgfdfghed"
					pause
			end if
			
			M_K_Moment(:, :, i, 1) = sqv * M_K_Moment(:, :, i, 1) / no
			
			if(par_PUI == .True.) then
				j = pui_num_tetr(i)
				if(j > 0) then
					pui_Sm(:, j) = pui_Sm(:, j) / no
					pui_Sp(:, j) = pui_Sp(:, j) / no
				end if
			end if

			if(par_pogloshenie == .True.) pogloshenie(:, :, i) = pogloshenie(:, :, i) * sqv / no
			
			if(MK_is_NaN == .True. .and. ieee_is_nan(M_K_Moment(1, 1, i, 1))) then
					print*, "NaN 098uiknhuuyhjh"
					pause
			end if

			do j = 1, par_n_sort
				if(M_K_Moment(1, j, i, 1) > 0.0000001) then
					M_K_Moment(2:4, j, i, 1) = M_K_Moment(2:4, j, i, 1)/M_K_Moment(1, j, i, 1)  ! Скорости
					M_K_Moment(5, j, i, 1) = (2.0/3.0) * ( M_K_Moment(5, j, i, 1)/M_K_Moment(1, j, i, 1) - &
						kvv(M_K_Moment(2, j, i, 1), M_K_Moment(3, j, i, 1), M_K_Moment(4, j, i, 1)) )  ! Temp
					
					if(par_n_moment > 9) then
						M_K_Moment(10, j, i, 1) = M_K_Moment(10, j, i, 1) / M_K_Moment(1, j, i, 1) - &
							M_K_Moment(2, j, i, 1)**2
						M_K_Moment(11, j, i, 1) = M_K_Moment(11, j, i, 1) / M_K_Moment(1, j, i, 1) - &
							M_K_Moment(2, j, i, 1)*M_K_Moment(3, j, i, 1)
						M_K_Moment(12, j, i, 1) = M_K_Moment(12, j, i, 1) / M_K_Moment(1, j, i, 1) - &
							M_K_Moment(2, j, i, 1)*M_K_Moment(4, j, i, 1)
						M_K_Moment(13, j, i, 1) = M_K_Moment(13, j, i, 1) / M_K_Moment(1, j, i, 1) - &
							M_K_Moment(3, j, i, 1)**2
						M_K_Moment(14, j, i, 1) = M_K_Moment(14, j, i, 1) / M_K_Moment(1, j, i, 1) - &
							M_K_Moment(3, j, i, 1)*M_K_Moment(4, j, i, 1)
						M_K_Moment(15, j, i, 1) = M_K_Moment(15, j, i, 1) / M_K_Moment(1, j, i, 1) - &
							M_K_Moment(4, j, i, 1)**2
						M_K_Moment(16, j, i, 1) = 2.0 * M_K_Moment(2, j, i, 1)**3 + 3.0 * M_K_Moment(2, j, i, 1) * M_K_Moment(10, j, i, 1) - &
							M_K_Moment(16, j, i, 1) / M_K_Moment(1, j, i, 1) 
						M_K_Moment(17, j, i, 1) = 2.0 * M_K_Moment(3, j, i, 1)**3 + 3.0 * M_K_Moment(3, j, i, 1) * M_K_Moment(13, j, i, 1) - &
							M_K_Moment(17, j, i, 1) / M_K_Moment(1, j, i, 1) 
						M_K_Moment(18, j, i, 1) = 2.0 * M_K_Moment(4, j, i, 1)**3 + 3.0 * M_K_Moment(4, j, i, 1) * M_K_Moment(15, j, i, 1) - &
							M_K_Moment(18, j, i, 1) / M_K_Moment(1, j, i, 1) 
					end if
					
				end if
			end do
			
			M_K_Moment(6:9, :, i, 1) = M_K_Moment(6:9, :, i, 1) * par_n_p_LISM
			M_K_Moment(19, :, i, 1) = M_K_Moment(19, :, i, 1) * par_n_p_LISM
		end do
		
		! Вне расчётной области нужно заполнить значения в тетраэдрах
		loop2: do i = 1, size(M_K_Moment(1, 1, :, 1))
			do j = 1, 4
				pp = int2_all_tetraendron_point(j, i)
				if(pp == 0) CYCLE loop2
				
				if((int2_coord(1, pp) >= 0.0 .and. norm2(int2_coord(:, pp)) >= par_Rmax)) then
					M_K_Moment(:, :, i, 1) = 0.0
					M_K_Moment(1, 4, i, 1) = 1.0
					M_K_Moment(5, 4, i, 1) = 1.0
					M_K_Moment(10, 4, i, 1) = 0.5
					M_K_Moment(13, 4, i, 1) = 0.5
					M_K_Moment(15, 4, i, 1) = 0.5
					M_K_Moment(2, 4, i, 1) = par_Velosity_inf
					CYCLE loop2
				end if
			end do
		end do loop2
		
		
		! Теперь нужно сохранить моменты не в тетраэдрах, а в их вершинах (с осреднением в вершине) ************************
		int2_Moment = 0.0
		allocate(vol_sr(size(int2_Moment(1, 1, :))))
		vol_sr = 0.0
		
		do i = 1, size(M_K_Moment(1, 1, :, 1))        ! Бежим по тетраэрам
			if(int2_all_tetraendron_point(1, i) == 0) CYCLE
			do j = 1, 4
				pp = int2_all_tetraendron_point(j, i)
				vol_sr(pp) = vol_sr(pp) + int2_all_Volume(i)
				int2_Moment(:, :, pp) = int2_Moment(:, :, pp) + M_K_Moment(:, :, i, 1) * int2_all_Volume(i)
			end do
		end do
		
		do i = 1, size(int2_Moment(1, 1, :))
			if(vol_sr(i) <= 0.0000001) then
				int2_Moment(:, :, i) = 0.0
				int2_Moment(1, 4, i) = 1.0
				int2_Moment(5, 4, i) = 1.0
				int2_Moment(2, 4, i) = par_Velosity_inf
				CYCLE
			end if
			int2_Moment(:, :, i) = int2_Moment(:, :, i) / vol_sr(i)
		end do
		
		deallocate(vol_sr)
		! ******************************************************************************************************************
		! Собираем статистику по весам
		if(MK_Mu_stat) then
			open(1, file = "MK_Mu_statistic_.txt")
			
			do k = 1, par_n_sort
				do j = 1, par_m_zone + 1
					do i = 1, par_n_zone + 1
						write(1, *) i, j, k, MK_Mu_statistic(i, j, k) * &
							(par_Rmax/MK_R_zone( min(i, par_n_zone) ))**2 * (par_m_zone + 1)/ MK_N
					end do
				end do
			end do
			
			close(1)
		end if
	
	end subroutine M_K_start
	
	subroutine M_K_sum()
		! Сбор источником из разных файлов
		USE OMP_LIB
		!$MPI include 'mpif.h'
		! Variables
		integer(4) :: num, i, j, pp
		real(8), allocatable :: vol_sr(:)                                    ! Для осреднения в узлах
		character(len=3) :: name
		real(8) :: no, N, N1
		logical :: exists
		real(8), allocatable :: M_K_Moment_print(:, :, :)
		
		call M_K_Set()    ! Создали массивы
		call M_K_init()
		
		
		! Собираем результаты со всех файлов
		allocate(M_K_Moment_print, mold = M_K_Moment(:, :, :, 1))
		
		M_K_Moment(:, :, :, 1) = 0.0
		N = 0
		
		do num = 1, 7
			write(unit = name, fmt='(i3.3)') num
			inquire(file="M-K_param_" // name // ".bin", exist=exists)
			if (exists == .False.) then
				pause "net faila!!!"
				STOP "net faila!!!"
			end if
			open(2, file = "M-K_param_" // name // ".bin", FORM = 'BINARY', ACTION = "READ")
			read(2)  N1
			N = N + N1
			close(2)
			
			print*, "N1 = ", N1
		end do
		
		do num = 1, 7
			write(unit = name, fmt='(i3.3)') num
		
			inquire(file="M-K_param_" // name // ".bin", exist=exists)
		
			if (exists == .False.) then
				pause "net faila!!!"
				STOP "net faila!!!"
			end if
		
			open(2, file = "M-K_param_" // name // ".bin", FORM = 'BINARY', ACTION = "READ")
		
			read(2)  N1
			read(2)  M_K_Moment_print
			M_K_Moment(:, :, :, 1) = M_K_Moment(:, :, :, 1) + M_K_Moment_print * (N1/N)
			
			print*, "(N1/N) = ", (N1/N), M_K_Moment_print(1, 1, 100), M_K_Moment_print(1, 1, 1000)
			close(2)
		end do
		
		! M_K_Moment(:, :, :, 1) = M_K_Moment(:, :, :, 1)/N
		deallocate(M_K_Moment_print)
		
		! Бежим по всем тетраэдрам и нормируем моменты
		do i = 1, size(M_K_Moment(1, 1, :, 1))
			
			if(int2_all_tetraendron_point(1, i) == 0) CYCLE
			
			if( int2_all_Volume(i) <= 0.00000001) then
				!print*, "Error  dfgdfgdg89346767809098742577"
				M_K_Moment(:, 4, i, 1) = 0.0
				continue
			end if
			
			!no = MK_Mu_mult * MK_N * int2_all_Volume(i)
			no = int2_all_Volume(i)
			
			if(MK_is_NaN == .True. .and. ieee_is_nan(no)) then
					print*, "NaN lj098inbh5dgfdfghed"
					pause
			end if
			
			M_K_Moment(:, :, i, 1) = sqv * M_K_Moment(:, :, i, 1) / no
			
			if(MK_is_NaN == .True. .and. ieee_is_nan(M_K_Moment(1, 1, i, 1))) then
					print*, "NaN 098uiknhuuyhjh"
					pause
			end if

			do j = 1, par_n_sort
				if(M_K_Moment(1, j, i, 1) > 0.000001) then
					M_K_Moment(2:4, j, i, 1) = M_K_Moment(2:4, j, i, 1)/M_K_Moment(1, j, i, 1)  ! Скорости
					M_K_Moment(5, j, i, 1) = (2.0/3.0) * ( M_K_Moment(5, j, i, 1)/M_K_Moment(1, j, i, 1) - &
						kvv(M_K_Moment(2, j, i, 1), M_K_Moment(3, j, i, 1), M_K_Moment(4, j, i, 1)) )  ! Temp
					
					if(par_n_moment > 9) then
						M_K_Moment(10, j, i, 1) = M_K_Moment(10, j, i, 1) / M_K_Moment(1, j, i, 1) - &
							M_K_Moment(2, j, i, 1)**2
						M_K_Moment(11, j, i, 1) = M_K_Moment(11, j, i, 1) / M_K_Moment(1, j, i, 1) - &
							M_K_Moment(2, j, i, 1)*M_K_Moment(3, j, i, 1)
						M_K_Moment(12, j, i, 1) = M_K_Moment(12, j, i, 1) / M_K_Moment(1, j, i, 1) - &
							M_K_Moment(2, j, i, 1)*M_K_Moment(4, j, i, 1)
						M_K_Moment(13, j, i, 1) = M_K_Moment(13, j, i, 1) / M_K_Moment(1, j, i, 1) - &
							M_K_Moment(3, j, i, 1)**2
						M_K_Moment(14, j, i, 1) = M_K_Moment(14, j, i, 1) / M_K_Moment(1, j, i, 1) - &
							M_K_Moment(3, j, i, 1)*M_K_Moment(4, j, i, 1)
						M_K_Moment(15, j, i, 1) = M_K_Moment(15, j, i, 1) / M_K_Moment(1, j, i, 1) - &
							M_K_Moment(4, j, i, 1)**2
						M_K_Moment(16, j, i, 1) = 2.0 * M_K_Moment(2, j, i, 1)**3 + 3.0 * M_K_Moment(2, j, i, 1) * M_K_Moment(10, j, i, 1) - &
							M_K_Moment(16, j, i, 1) / M_K_Moment(1, j, i, 1) 
						M_K_Moment(17, j, i, 1) = 2.0 * M_K_Moment(3, j, i, 1)**3 + 3.0 * M_K_Moment(3, j, i, 1) * M_K_Moment(13, j, i, 1) - &
							M_K_Moment(17, j, i, 1) / M_K_Moment(1, j, i, 1) 
						M_K_Moment(18, j, i, 1) = 2.0 * M_K_Moment(4, j, i, 1)**3 + 3.0 * M_K_Moment(4, j, i, 1) * M_K_Moment(15, j, i, 1) - &
							M_K_Moment(18, j, i, 1) / M_K_Moment(1, j, i, 1) 
					end if
					
				end if
			end do
			
			M_K_Moment(6:9, :, i, 1) = M_K_Moment(6:9, :, i, 1) * par_n_p_LISM
			M_K_Moment(19, :, i, 1) = M_K_Moment(19, :, i, 1) * par_n_p_LISM
		end do
		
		! Вне расчётной области нужно заполнить значения в тетраэдрах
		loop2: do i = 1, size(M_K_Moment(1, 1, :, 1))
			do j = 1, 4
				pp = int2_all_tetraendron_point(j, i)
				if(pp == 0) CYCLE loop2
				
				if((int2_coord(1, pp) >= 0.0 .and. norm2(int2_coord(:, pp)) >= par_Rmax)) then
					M_K_Moment(:, :, i, 1) = 0.0
					M_K_Moment(1, 4, i, 1) = 1.0
					M_K_Moment(5, 4, i, 1) = 1.0
					M_K_Moment(10, 4, i, 1) = 0.5
					M_K_Moment(13, 4, i, 1) = 0.5
					M_K_Moment(15, 4, i, 1) = 0.5
					M_K_Moment(2, 4, i, 1) = par_Velosity_inf
					CYCLE loop2
				end if
			end do
		end do loop2
		
		
		! Теперь нужно сохранить моменты не в тетраэдрах, а в их вершинах (с осреднением в вершине) ************************
		int2_Moment = 0.0
		allocate(vol_sr(size(int2_Moment(1, 1, :))))
		vol_sr = 0.0
		
		do i = 1, size(M_K_Moment(1, 1, :, 1))        ! Бежим по тетраэрам
			if(int2_all_tetraendron_point(1, i) == 0) CYCLE
			do j = 1, 4
				pp = int2_all_tetraendron_point(j, i)
				vol_sr(pp) = vol_sr(pp) + int2_all_Volume(i)
				int2_Moment(:, :, pp) = int2_Moment(:, :, pp) + M_K_Moment(:, :, i, 1) * int2_all_Volume(i)
			end do
		end do
		
		do i = 1, size(int2_Moment(1, 1, :))
			if(vol_sr(i) <= 0.0000001) then
				int2_Moment(:, :, i) = 0.0
				int2_Moment(1, 4, i) = 1.0
				int2_Moment(5, 4, i) = 1.0
				int2_Moment(2, 4, i) = par_Velosity_inf
				CYCLE
			end if
			int2_Moment(:, :, i) = int2_Moment(:, :, i) / vol_sr(i)
		end do
		
		deallocate(vol_sr)
	
	end subroutine M_K_sum
	
	subroutine M_K_init()
		! Variables
		real(8) :: Y, betta
		real(8) :: PAR(9) 
		integer :: cell

		print*, "M_K_init()"
		
		! Инициализация некоторых параметров
		par_n_moment = 19
		betta = 2.0 * par_pi_8/(par_l_phi - 2)
		par_Rleft = par_R_LEFT + 0.0001
		par_Rup = par_R_END - 2.0_8! par_R_END * sqrt( 1.0 - 0.5 * sin(betta)**2 / (1.0 + cos(betta)) ) - 0.0001  !  Верхняя стенка
		print*, "par_Rup = ", par_Rup
		
		Y = dabs(par_Velosity_inf)
		sqv_1 = (par_Rmax) * (0.5 * (par_Rmax) * par_pi_8 * Y * &
			(erf(Y) * (1.0 + 1.0 / (2.0 * (Y**2))) + 1.0 + exp(-(Y**2)) / (Y * par_sqrtpi)))  ! Запуск с полусферы
		sqv_2 = par_sqrtpi * (par_Rup) * dabs(par_Rleft) ! Верхняя стенка
		sqv_3 = par_pi_8 * (par_Rup**2) * exp(-(par_Velosity_inf**2)) * &
		(1.0 + exp(par_Velosity_inf**2) * par_sqrtpi * par_Velosity_inf * (1.0 + erf(par_Velosity_inf))) / (2.0 * par_sqrtpi)
		sqv_4 = (par_sqrtpi / 2.0) * ((par_Rmax**2) - (par_Rup**2)) * &
			exp(-(par_Velosity_inf**2)) * (par_sqrtpi * par_Velosity_inf * erfc(par_Velosity_inf) * exp(par_Velosity_inf**2) - 1.0) ! Передняя стенка (выше полусферы)
		sqv = sqv_1 + sqv_2 + sqv_3 + sqv_4
		
		
		
		MK_N = MK_N1 + MK_N2 + MK_N3 + MK_N4  
		MK_mu1 = (sqv_1/ sqv) * (1.0 * MK_N / MK_N1)
		MK_mu2 = (sqv_2/ sqv) * (1.0 * MK_N / MK_N2)
		MK_mu3 = (sqv_3/ sqv) * (1.0 * MK_N / MK_N3)
		MK_mu4 = (sqv_4/ sqv) * (1.0 * MK_N / MK_N4)
		! Body of M_K_init
		MK_N = MK_N * par_n_potok * par_n_parallel
		
		
		! Проверка параметров
		cell = 3
		call Int2_Get_par_fast(260.0_8, 1.0_8, 1.0_8, cell, PAR)
		print*, "Proverka MK, 745465635", "must be 1 = ", sqrt(PAR(5)/PAR(1)), PAR(1)
	end subroutine M_K_init
	
	subroutine M_K_Fly(n_potok)
		! Функция запускает все частицы в стеке потока + все дочерние частицы
		
		integer(4), intent(in) :: n_potok  ! Номер потока 
		
		real(8) :: particle(8)
		integer(4):: particle_2(4), i, ijk
		logical :: particle_3(par_n_zone + 1, par_m_zone + 1)
		
		integer(4) :: num  ! Номер частицы, верхняя в стеке
		integer(4) :: cell ! Номер ячейки, в которой находится частица
		integer(4) :: next ! Номер ячейки, в которую попадёт частица в следующий раз
		integer(4) :: area2  ! Зона, в которой сейчас находится ячейка
		integer(4) :: II  ! На сколько атомов расщепляется атом при перезарядке
		integer(4) :: to_i, to_j, from_i, from_j
		logical :: bb, bb2
		
		real(8) :: time ! Оценочное время до вылета частицы из ячейки
		
		real(8) :: cp, vx, vy, vz, ro, p, PAR(9)  ! Параметры плазмы в ячейке
		real(8) :: uz, nu_ex, kappa, ksi, t_ex, t2, mu_ex, mu2, r_ex(3), r, mu, u, V(3), mu3
		real(8) :: uz_M, uz_E, k1, k2, k3, u1, u2, u3, skalar
		real(8) :: Ur, Uthe, Uphi, Vr, Vthe, Vphi
		real(8) :: v1, v2, v3, r_peregel, ddt, rho_He
		
		real(8) :: nu_ph, kappa_ph, kappa_all, mu_ph, mu_perez
		
		real(8) :: Wr(par_n_zone + 1), Wthe(par_n_zone + 1), Wphi(par_n_zone + 1), mu_(par_n_zone + 1)
		
		integer(4) :: step
	
		step = 0
		
		do while (stek(n_potok) >= 1)
			
			step = step + 1
			
			!if (mod(step, 1) == 0) print*, "step = ", step
			
			!pause
			
			if(stek(n_potok) > par_stek * 0.95) then
				print*, "1234543fj976r  Perepolnen stek", stek(n_potok) 
				!pause
				STOP
			end if
			
				
			num = stek(n_potok)
			stek(n_potok) = stek(n_potok) - 1
			! Берём все параметры частицы
			particle = M_K_particle(:, num, n_potok)
			particle_2 = M_K_particle_2(:, num, n_potok)
			if(MK_Mu_stat) particle_3 = M_K_particle_3(:, :, num, n_potok)
			!print*, "stek(n_potok) = ", stek(n_potok)
			!print*, particle
			!print*, "______"
			!print*, particle_2
			!pause
			
			
			loop1: do  ! пока частица не вылетит из области
				cell = particle_2(1)
				mu = particle(7)                                ! Вес частицы
			
				call Int2_Time_fly(particle(1:3), particle(4:6), time, cell, next)  ! Находим время time до вылета из ячейки
				
				time = max(0.00000001_8, time * 1.0001) ! Увеличим время, чтобы частица точно вышла из ячейки
				
				area2 = int2_Cell_par2(1, int2_all_tetraendron_point(1, cell)) ! Зона рождения
				r = norm2(particle(1:3) + time/2.0 * particle(4:6))

				kappa = 0.0
				do ijk = 1, 3
					
					select case (ijk)
						case(1)
							ddt = 1.0/6.0                                      
						case(2)
							ddt = 5.0/6.0                                    
						case(3)
							ddt = 3.0/6.0                                  
						case default
							print*, "Error uijkhjgfbnbnn hbuhuefw"
							STOP
					end select
					
					call Int2_Get_par_fast2(particle(1) + time * ddt * particle(4), particle(2)+ time * ddt * particle(5),&
						particle(3) + time * ddt * particle(6), cell, PAR = PAR, rho_He = rho_He)
				
					call Sootnosheniya(PAR(1), PAR(5), rho_He, 0.0_8, 0.0_8, area2, rho_Th = ro, p_Th = p)

					if(p < 0.0) p = par_p_0/((r/par_1ae)**(2.0 * ggg))

					cp = sqrt(2.0 * p/ro)
					vx = PAR(2)
					vy = PAR(3)
					vz = PAR(4)
					!ro = PAR(1)
				
					if(ro <= 0.0 .or. ro > 1000.0) then
						print*, PAR
						print*, "___"
						print*, cell, int2_all_tetraendron_point(:, cell)
						pause "ERROR ro MK 157 6787yutr4dfghhghjuhj0089"
					end if
				
					! Найдём время до перезарядки и веса частиц  ****************************************************************************************
					u = sqrt(kvv(particle(4) - vx, particle(5) - vy, particle(6) - vz))
					u1 =  vx - particle(4)
					u2 =  vy - particle(5)
					u3 =  vz - particle(6)
					skalar = particle(4) * u1 + particle(5) * u2 + particle(6) * u3
				
					if (u / cp > 7.0) then
						uz = MK_Velosity_1(u, cp);
						nu_ex = (ro * uz * MK_sigma(uz)) / par_Kn
					else
						nu_ex = (ro * MK_int_1(u, cp)) / par_Kn        ! Пробуем вычислять интеграллы численно
					end if
			
					kappa = kappa + (nu_ex * time/3.0)  ! по перезарядке
				end do
				
				
				
				
				if(MK_photoionization) then
					nu_ph = par_nu_ph * (par_1ae/r)**2
					kappa_ph = (nu_ph * time)     ! по фотоионизации
				end if
				
				kappa_all = kappa
				if(MK_photoionization) kappa_all = kappa_all + kappa_ph
				
				call M_K_rand(sensor(1, 2, n_potok), sensor(2, 2, n_potok), sensor(3, 2, n_potok), ksi)
				
				t_ex = -(time / kappa_all) * log(1.0 - ksi * (1.0 - exp(-kappa_all)))  ! Время до перезарядки
				t2 = time - t_ex  ! Время сколько лететь после того, как атом перезарядился
				mu_perez = mu * (1.0 - exp(-kappa_all)) ! вес перезаряженного атома по всем процессам
				mu2 = mu * exp(-kappa_all)  ! вес оставшегося неперезаряженного атома
				mu_ph = 0.0
				if(MK_photoionization) mu_ph = (kappa_ph / kappa_all) * mu_perez  ! Вес ионизированного атома
				mu_ex = mu_perez - mu_ph      ! Вес перезаряженного на протонах атома
				
				
				if(mu2 < 0.0) then
					print*, mu2, mu, mu_ex, kappa, kappa_all
					print*, particle(7)
					print*, "p = ", particle 
					STOP "Eror oiuyyuiojhu987uio9i1313143"
				end if
				

				r_ex = particle(1:3) + t_ex * particle(4:6)   ! Координаты перезарядки
				
				! проверка, если перезарядка произошла за пределами ячейки (нужно немного сдвигать её в этом случае)
				if (Belong_tetraedron(r_ex(1), r_ex(2), r_ex(3), cell) == .False.) then
					t_ex = time * 0.998
					t2 = time - t_ex
					r_ex = particle(1:3) + t_ex * particle(4:6)  
				end if
				
				r = norm2(r_ex)  ! Расстояние от точки перезарядки до Солнца
				
				from_i = MK_geo_zones(r, 1.0_8)     ! Зона по r в точке перезарядки
				from_j = MK_alpha_zones( polar_angle( r_ex(1), sqrt(r_ex(2)**2 + r_ex(3)**2) ) ) ! Зона по углу в точке перезарядки
				
				!print*, from_i, from_j, r_ex, polar_angle( r_ex(1), sqrt(r_ex(2)**2 + r_ex(3)**2) )
				!pause
				
				if(MK_Mu_stat .and. particle_3(from_i, from_j) == .False.) then
					particle_3(from_i, from_j) = .True.
					
					!$omp critical
					MK_Mu_statistic(from_i, from_j, particle_2(2)) = MK_Mu_statistic(from_i, from_j, particle_2(2)) + &
						mu/max(0.3 * MK_SINKR(from_j), sin(polar_angle( r_ex(1), sqrt(r_ex(2)**2 + r_ex(3)**2) )))
					!$omp end critical
				end if
				
				
				
				! Накапливаем моменты и т.д. ______________________________________________________________________________________________________________________________
				if(MK_is_NaN == .True. .and. ieee_is_nan(t_ex * mu + t2 * mu2)) then
					print*, "NaN 789olhgyuimnhyuiolkuytgf"
					pause
				end if
				
				M_K_Moment(1, particle_2(2), cell, n_potok) = M_K_Moment(1, particle_2(2), cell, n_potok) + t_ex * mu + t2 * mu2
				M_K_Moment(2:4, particle_2(2), cell, n_potok) = M_K_Moment(2:4, particle_2(2), cell, n_potok) + (t_ex * mu + t2 * mu2) * particle(4:6)
				M_K_Moment(5, particle_2(2), cell, n_potok) = M_K_Moment(5, particle_2(2), cell, n_potok) + &
					(t_ex * mu + t2 * mu2) * kvv(particle(4), particle(5), particle(6))
				
				if(par_n_moment > 9) then
					M_K_Moment(10, particle_2(2), cell, n_potok) = M_K_Moment(10, particle_2(2), cell, n_potok) + &
						(t_ex * mu + t2 * mu2) * particle(4)**2
					M_K_Moment(11, particle_2(2), cell, n_potok) = M_K_Moment(11, particle_2(2), cell, n_potok) + &
						(t_ex * mu + t2 * mu2) * particle(4)*particle(5)
					M_K_Moment(12, particle_2(2), cell, n_potok) = M_K_Moment(12, particle_2(2), cell, n_potok) + &
						(t_ex * mu + t2 * mu2) * particle(4)*particle(6)
					M_K_Moment(13, particle_2(2), cell, n_potok) = M_K_Moment(13, particle_2(2), cell, n_potok) + &
						(t_ex * mu + t2 * mu2) * particle(5)*particle(5)
					M_K_Moment(14, particle_2(2), cell, n_potok) = M_K_Moment(14, particle_2(2), cell, n_potok) + &
						(t_ex * mu + t2 * mu2) * particle(5)*particle(6)
					M_K_Moment(15, particle_2(2), cell, n_potok) = M_K_Moment(15, particle_2(2), cell, n_potok) + &
						(t_ex * mu + t2 * mu2) * particle(6)*particle(6)
					M_K_Moment(16, particle_2(2), cell, n_potok) = M_K_Moment(16, particle_2(2), cell, n_potok) + &
						(t_ex * mu + t2 * mu2) * particle(4)**3
					M_K_Moment(17, particle_2(2), cell, n_potok) = M_K_Moment(17, particle_2(2), cell, n_potok) + &
						(t_ex * mu + t2 * mu2) * particle(5)**3
					M_K_Moment(18, particle_2(2), cell, n_potok) = M_K_Moment(18, particle_2(2), cell, n_potok) + &
						(t_ex * mu + t2 * mu2) * particle(6)**3
				end if
				
				
				if (u / cp > 7.0) then
					uz_M = MK_Velosity_2(u, cp)/ (uz * (cp**2) * cp * par_pi_8 * par_sqrtpi)
					uz_E = MK_Velosity_3(u, cp)
					
					if(MK_is_NaN == .True. .and. ieee_is_nan(uz_E)) then
						print*, "NaN 34rtyu765rdfgh"
						pause
					end if
					
					M_K_Moment(6, particle_2(2), cell, n_potok) = M_K_Moment(6, particle_2(2), cell, n_potok) - mu_ex * uz_M * u1 / u
					M_K_Moment(7, particle_2(2), cell, n_potok) = M_K_Moment(7, particle_2(2), cell, n_potok) - mu_ex * uz_M * u2 / u
					M_K_Moment(8, particle_2(2), cell, n_potok) = M_K_Moment(8, particle_2(2), cell, n_potok) - mu_ex * uz_M * u3 / u
					M_K_Moment(9, particle_2(2), cell, n_potok) = M_K_Moment(9, particle_2(2), cell, n_potok) + &
						mu_ex * (-0.25 * (3.0 * cp**2 + 2.0 * u**2) * (uz_E / uz) - uz_M * skalar / u)
				else
					k1 = MK_int_1(u, cp)
					k2 = MK_int_2(u, cp)
					k3 = MK_int_3(u, cp)
					
					if( MK_is_NaN == .True. .and. (ieee_is_nan(k1) .or. ieee_is_nan(k2) .or. ieee_is_nan(k3))) then
						print*, "NaN klkjhgfddfujyukjhjhgfghozpem"
						pause
					end if
					
					M_K_Moment(6, particle_2(2), cell, n_potok) = M_K_Moment(6, particle_2(2), cell, n_potok) + mu_ex * (k2/k1) * u1 / u
					M_K_Moment(7, particle_2(2), cell, n_potok) = M_K_Moment(7, particle_2(2), cell, n_potok) + mu_ex * (k2/k1) * u2 / u
					M_K_Moment(8, particle_2(2), cell, n_potok) = M_K_Moment(8, particle_2(2), cell, n_potok) + mu_ex * (k2/k1) * u3 / u
					M_K_Moment(9, particle_2(2), cell, n_potok) = M_K_Moment(9, particle_2(2), cell, n_potok) + &
						mu_ex * (-0.5 * k3/k1 + k2/k1 * skalar / u)
					
				end if
				
				! Добавим интеграллы по фотоионизации  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				M_K_Moment(19, particle_2(2), cell, n_potok) = M_K_Moment(19, particle_2(2), cell, n_potok) + mu_ph 
				M_K_Moment(6:8, particle_2(2), cell, n_potok) = M_K_Moment(6:8, particle_2(2), cell, n_potok) + mu_ph * particle(4:6)
				M_K_Moment(9, particle_2(2), cell, n_potok) = M_K_Moment(9, particle_2(2), cell, n_potok) + &
					mu_ph * (0.5 * norm2(particle(4:6)) + par_E_ph)
				
				! Добавляем для расчёта PUI
				if(area2 <= 2 .and. par_PUI == .True.) then ! 
					!? call PUI_Add(cell, u, kappa/time, nu_ph, mu_ex, mu_ph, t_ex, time)
					!call PUI_Add(cell, u, kappa/time, nu_ph, mu, t_ex)
					!call PUI_Add(cell, u, kappa/time, nu_ph, mu2, t2)
					call PUI_Add(cell, u, kappa_all/time, mu, t_ex)
					call PUI_Add(cell, u, kappa_all/time, mu2, t2)
				end if
				!!(t_ex * mu + t2 * mu2)
				!_________________________________________________________________________________________
				
				!! Добавляем расчёт поглощения  ---------------------------------------------------
				if(par_pogloshenie == .True.) then
					uz_M = DOT_PRODUCT(particle(4:6), r_ex)/r   ! Проекция скорости на ось зрения (от звёзды)
					if(uz_M > pogl_v_min .and. uz_M < pogl_v_max) then
						ijk = MK_poglosh_nomer(uz_M)
						call omp_set_lock(pogl_lock(cell))
						pogloshenie(particle_2(2), ijk, cell) = pogloshenie(particle_2(2), ijk, cell) + t_ex * mu + t2 * mu2
						call omp_unset_lock(pogl_lock(cell))
					end if
				end if
				!! --------------------------------------------------------------------------------
				
				
				
				call spherical_skorost(r_ex(1), r_ex(2), r_ex(3), vx, vy, vz, Ur, Uphi, Uthe)
				call spherical_skorost(r_ex(1), r_ex(2), r_ex(3), particle(4), particle(5), particle(6), Vr, Vphi, Vthe)
			
				! Перезаряжаем
				if (area2 == 1 .or. Ur / cp > 1.8) then   ! Без геометрического расщепления
					II = 0  ! II - это сколько дополнительных атомов запускается (помимо основного)
				else
					II = MK_geo_zones(r, 1.2_8) - 1
				end if
				
				call M_K_Change_Velosity4(n_potok, Ur/cp, Uthe/cp, Uphi/cp, Vr/cp, Vthe/cp, Vphi/cp, Wr, Wthe, Wphi, mu_, &
					cp, r, II, r_ex(1), r_ex(2), r_ex(3), bb)
				
				Wr = Wr * cp
				Wthe = Wthe * cp
				Wphi = Wphi * cp
				
				do i = 1, II + 1
					if(i == II + 1 .and. bb == .False.) CYCLE  ! Если не запускается основной атом
					call dekard_skorost(r_ex(1), r_ex(2), r_ex(3), Wr(i), Wphi(i), Wthe(i), v1, v2, v3)
					V = (/ v1, v2, v3 /)
					
					call MK_Distination(r_ex, V, to_i, to_j, r_peregel)  ! Находим зону назначения в точке перегелия для рулетки
					
					mu3 = mu_ex * mu_(i)
					call MK_ruletka(n_potok, to_i, to_j, from_i, from_j, area2, r, r_peregel, mu3, bb2)

					if(mu3 < 0.0) then
						print*, "ERROR 209iur92yr8y8p2"
						print*, mu_ex
						print*, i, "  =  ", mu_
						STOP
					end if
					
					if(bb2 == .False.) CYCLE  ! Не запускаем эту частицу, она вырубается
					
					! Запишем новую частицу в стек
					stek(n_potok) = stek(n_potok) + 1
					M_K_particle(1:3, stek(n_potok), n_potok) = r_ex
					M_K_particle(4:6, stek(n_potok), n_potok) = V
					M_K_particle(7, stek(n_potok), n_potok) = mu3
					M_K_particle(8, stek(n_potok), n_potok) = r_peregel

					M_K_particle_2(:,stek(n_potok), n_potok) = (/ cell, area2, to_i, to_j /)
					
					if(MK_Mu_stat == .True.) then
						if(particle_2(2) == area2) then
							M_K_particle_3(:, :, stek(n_potok), n_potok) = particle_3
						else 
							M_K_particle_3(:, :, stek(n_potok), n_potok) = .False.
						end if
					end if
					
					! (par_n_zone + 1, par_m_zone + 1, par_stek, число потоков)
					
				end do
				
				! Проверяем, нужно ли вырубить неперезаряженную часть атома
				call MK_ruletka(n_potok, particle_2(3), particle_2(4), from_i, from_j, particle_2(2), &
						r, particle(8), mu2, bb2)
				if(bb2 == .False.) EXIT loop1  ! Не запускаем эту частицу, она вырубается
				particle(7) = mu2
			
				! Находим следующую ячейку
				
				particle(1:3) = particle(1:3) + time * particle(4:6)
				if(next == 0) then
					if (norm2(particle(1:3)) <= 100.0) then
						print*, "Error 23e2323 ", particle(1:3)
						pause
					end if
					EXIT loop1  ! частица долетела до края области
				end if
				
		11			continue 
				
				
				
				
				call Int2_Get_tetraedron(particle(1), particle(2), particle(3), next)  ! Находим точно следующий тетраэдр
				if (next == -1) then ! ЗАТЫК
					next = 3
					particle(1:3) = particle(1:3) + (time/100.0) * particle(4:6)
					print*, "Zatik"
					pause
					GO TO 11
				end if
				
				particle_2(1) = next
				if(next == 0) then
					if (norm2(particle(1:3)) <= 100.0) then
						print*, "Error wqer2r42  ", particle(1:3)
						pause
					end if
					
					EXIT loop1  ! частица долетела до края области
				end if
				
				if(particle(1) > 0.00001 .and. norm2(particle(1:3)) > par_Rmax + 0.001) then
					EXIT loop1  ! частица долетела до края области
				end if
				
			
				
			
			end do loop1
			
			
		end do
		
	end subroutine M_K_Fly
	
	subroutine M_K_Fly_PUI(n_potok)
		! Функция также учитывает пикапы!
		! Функция запускает все частицы в стеке потока + все дочерние частицы
		
		integer(4), intent(in) :: n_potok  ! Номер потока 
		
		real(8) :: particle(8)
		integer(4):: particle_2(4), i, ijk, jjj
		logical :: particle_3(par_n_zone + 1, par_m_zone + 1)
		
		integer(4) :: num  ! Номер частицы, верхняя в стеке
		integer(4) :: cell ! Номер ячейки-тетраэдра, в которой находится частица
		integer(4) :: cell2 ! Номер ячейки-основной, в которой находится частица
		integer(4) :: next ! Номер ячейки, в которую попадёт частица в следующий раз
		integer(4) :: area2  ! Зона, в которой сейчас находится ячейка
		integer(4) :: new_sort  ! атом какого сорта родился
		integer(4) :: II  ! На сколько атомов расщепляется атом при перезарядке
		integer(4) :: to_i, to_j, from_i, from_j
		logical :: bb, bb2
		
		real(8) :: time ! Оценочное время до вылета частицы из ячейки
		
		real(8) :: cp, vx, vy, vz, ro, PAR(9), MAS_PUI(2), p, ro_pui  ! Параметры плазмы в ячейке
		real(8) :: pui_less
		real(8) :: uz, nu_ex, kappa, ksi, t_ex, t2, mu_ex, mu2, r_ex(3), r, mu, u, V(3), mu3
		real(8) :: nu_ex_pui, kappa_pui, mu_ex_pui
		real(8) :: uz_M, uz_E, k1, k2, k3, u1, u2, u3, skalar
		real(8) :: Ur, Uthe, Uphi, Vr, Vthe, Vphi
		real(8) :: v1, v2, v3, r_peregel, ddt, rho_He
		
		real(8) :: nu_ph, kappa_ph, kappa_all, mu_ph, mu_perez
		
		real(8) :: Wr(par_n_zone + 1), Wthe(par_n_zone + 1), Wphi(par_n_zone + 1), mu_(par_n_zone + 1)
		
		integer(4) :: step
	
		step = 0
		
		
		do while (stek(n_potok) >= 1)
			
			step = step + 1
			
			!if (mod(step, 1) == 0) print*, "step = ", step
			
			!pause
			
			if(stek(n_potok) > par_stek * 0.95) then
				print*, "1234543fj976r  Perepolnen stek", stek(n_potok) 
				!pause
				STOP
			end if
			
				
			num = stek(n_potok)
			stek(n_potok) = stek(n_potok) - 1
			! Берём все параметры частицы
			particle = M_K_particle(:, num, n_potok)
			particle_2 = M_K_particle_2(:, num, n_potok)
			if(MK_Mu_stat) particle_3 = M_K_particle_3(:, :, num, n_potok)
			
			
			loop1: do  ! пока частица не вылетит из области
				cell = particle_2(1)
				mu = particle(7)                                ! Вес частицы
			
				call Int2_Time_fly(particle(1:3), particle(4:6), time, cell, next)  ! Находим время time до вылета из ячейки
				
				time = max(0.00000001_8, time * 1.0001) ! Увеличим время, чтобы частица точно вышла из ячейки
				cell2 = int2_all_tetraendron_point(1, cell) ! INT(cell/6) + 1  ! Номер ячейки, где брать pui

				area2 = int2_Cell_par2(1, int2_all_tetraendron_point(1, cell)) ! Зона рождения

				! if(cell == 334030) then
				! 	print*, "cell = ", cell
				! 	print*, particle(1:3)
				! 	print*, area2
				! 	pause
				! end if

				kappa = 0.0
				kappa_pui = 0.0
				r = norm2(particle(1:3) + time/2.0 * particle(4:6))
				pui_less = 1.0
				MAS_PUI = 0.0

				do ijk = 1, 3  ! Разбиение траектории на 3 части для более точной перезарядки
					pui_less = 1.0
					select case (ijk)
						case(1)
							ddt = 1.0/6.0                                      
						case(2)
							ddt = 5.0/6.0                                    
						case(3)
							ddt = 3.0/6.0                                  
						case default
							print*, "Error uijkhjgfbnbnn hbuhuefw"
							STOP
					end select
					
					if(area2 <= 2) then
						call Int2_Get_par_fast2(particle(1) + time * ddt * particle(4), particle(2)+ time * ddt * particle(5),&
							particle(3) + time * ddt * particle(6), cell, PAR, MAS_PUI = MAS_PUI, rho_He = rho_He)
					else
						call Int2_Get_par_fast2(particle(1) + time * ddt * particle(4), particle(2)+ time * ddt * particle(5),&
							particle(3) + time * ddt * particle(6), cell, PAR, rho_He = rho_He)
					end if
				
					! cp = sqrt(PAR(5)/PAR(1))
					if(area2 <= 2) then
						!! ro = PAR(1) + rho_He
						!! p = PAR(5) * (8.0 * PAR(1)  + 3.0 * rho_He) / (8.0 * PAR(1))
						!! ro_pui = MAS_PUI(1)
						!! p = 4.0 * (ro - ro_pui - rho_He) * (p - MAS_PUI(2) * MAS_PUI(1)) /&
						!! (8 * ro - 5.0 * rho_He - 4.0 * ro_pui)
						!! ro = ro - ro_pui - rho_He

						ro_pui = MAS_PUI(1)

						! if(ro_pui > 0) then
						! 	print*, "ro_pui ", ro_pui
						! 	pause
						! end if

						call Sootnosheniya(PAR(1), PAR(5), rho_He, ro_pui, MAS_PUI(2), 2, rho_Th = ro, p_Th = p)

						if(p < 0.0) p = par_p_0/((r/par_1ae)**(2.0 * ggg))
						
						if(ro < 0.000001) then
							ro = 0.000001
							p = 0.000001
							ro_pui = PAR(1) - ro
							pui_less = ro_pui/MAS_PUI(1)
						end if
						cp = sqrt(2.0 * p/ro)
					else
						!!p = PAR(5)
						!!ro = PAR(1)
						call Sootnosheniya(PAR(1), PAR(5), rho_He, 0.0_8, 0.0_8, 1, rho_Th = ro, p_Th = p)

						ro_pui = 0.0
						cp = sqrt(2.0 * p/ro)
					end if

					vx = PAR(2)
					vy = PAR(3)
					vz = PAR(4)
					
					!! Нужно учесть изменение температуры и плотности протонов из-за отдельного учёта пикапов
				
					if(ro <= 0.0 .or. ro > 1000.0) then
						print*, "ro =", ro
						print*, area2
						print*, "___"
						print*, particle(1:3)
						print*, "___"
						print*, MAS_PUI
						print*, "___"
						print*, PAR
						print*, "___"
						print*, cell, int2_all_tetraendron_point(:, cell)
						print*, "___"
						print*, n_pui(f_pui_num2(104308))
						print*, "___"
						print*, n_pui(f_pui_num2(104368))
						print*, "___"
						print*, n_pui(f_pui_num2(34788))
						print*, "___"
						print*, n_pui(f_pui_num2(104366))
						pause "ERROR ro MK 157 dvrbtrbrv34bt5vestves"
					end if
				
					! Найдём время до перезарядки и веса частиц  ****************************************************************************************
					u = sqrt(kvv(particle(4) - vx, particle(5) - vy, particle(6) - vz))
					u1 =  vx - particle(4)
					u2 =  vy - particle(5)
					u3 =  vz - particle(6)
					skalar = particle(4) * u1 + particle(5) * u2 + particle(6) * u3
				
					if (u / cp > 7.0) then
						uz = MK_Velosity_1(u, cp);
						nu_ex = (ro * uz * MK_sigma(uz)) / par_Kn
					else
						nu_ex = (ro * MK_int_1(u, cp)) / par_Kn        ! Пробуем вычислять интеграллы численно
					end if

					if(area2 <= 2 .and. par_PUI == .True.) kappa_pui = kappa_pui + pui_less * PUI_get_nu_integr(f_pui_num2(cell2), u)/ par_Kn * time/3.0  
			
					kappa = kappa + (nu_ex * time/3.0)  ! по перезарядке
				end do

				
				kappa_ph = 0.0_8
				if(MK_photoionization) then
					nu_ph = par_nu_ph * (par_1ae/r)**2
					kappa_ph = (nu_ph * time)     ! по фотоионизации
				end if
				
				kappa_all = kappa + kappa_pui
				if(MK_photoionization) kappa_all = kappa_all + kappa_ph
				
				call M_K_rand(sensor(1, 2, n_potok), sensor(2, 2, n_potok), sensor(3, 2, n_potok), ksi)
				
				t_ex = -(time / kappa_all) * log(1.0 - ksi * (1.0 - exp(-kappa_all)))  ! Время до перезарядки
				t2 = time - t_ex  ! Время сколько лететь после того, как атом перезарядился
				mu_perez = mu * (1.0 - exp(-kappa_all)) ! вес перезаряженного атома по всем процессам
				mu2 = max(mu - mu_perez, 0.0_8)   !mu * exp(-kappa_all)  ! вес оставшегося неперезаряженного атома 
				mu_ph = 0.0_8
				if(MK_photoionization) mu_ph = (kappa_ph / kappa_all) * mu_perez  ! Вес ионизированного атома
				mu_ex = (kappa / kappa_all) * mu_perez  ! mu_perez - mu_ph      ! Вес перезаряженного на протонах атома
				mu_ex_pui = (kappa_pui / kappa_all) * mu_perez  ! mu_perez - mu_ph      ! Вес перезаряженного на пикапах атома
				
				
				if(mu2 < 0.0) then
					print*, mu2, mu, mu_ex, kappa
					STOP "Eror oiuyyuiojhu987uio9i"
				end if
				

				r_ex = particle(1:3) + t_ex * particle(4:6)   ! Координаты перезарядки
				
				! проверка, если перезарядка произошла за пределами ячейки (нужно немного сдвигать её в этом случае)
				if (Belong_tetraedron(r_ex(1), r_ex(2), r_ex(3), cell) == .False.) then
					t_ex = time * 0.998
					t2 = time - t_ex
					r_ex = particle(1:3) + t_ex * particle(4:6)  
				end if
				
				r = norm2(r_ex)  ! Расстояние от точки перезарядки до Солнца
				
				from_i = MK_geo_zones(r, 1.0_8)     ! Зона по r в точке перезарядки
				from_j = MK_alpha_zones( polar_angle( r_ex(1), sqrt(r_ex(2)**2 + r_ex(3)**2) ) ) ! Зона по углу в точке перезарядки
				
				!print*, from_i, from_j, r_ex, polar_angle( r_ex(1), sqrt(r_ex(2)**2 + r_ex(3)**2) )
				!pause
				
				if(MK_Mu_stat .and. particle_3(from_i, from_j) == .False.) then
					particle_3(from_i, from_j) = .True.
					
					!$omp critical
					MK_Mu_statistic(from_i, from_j, particle_2(2)) = MK_Mu_statistic(from_i, from_j, particle_2(2)) + &
						mu/max(0.3 * MK_SINKR(from_j), sin(polar_angle( r_ex(1), sqrt(r_ex(2)**2 + r_ex(3)**2) )))
					!$omp end critical
				end if
				
				
				
				! Накапливаем моменты и т.д. ______________________________________________________________________________________________________________________________
				if(MK_is_NaN == .True. .and. ieee_is_nan(t_ex * mu + t2 * mu2)) then
					print*, "NaN 789olhgyuimnhyuiolkuytgf"
					pause
				end if
				
				M_K_Moment(1, particle_2(2), cell, n_potok) = M_K_Moment(1, particle_2(2), cell, n_potok) + t_ex * mu + t2 * mu2
				M_K_Moment(2:4, particle_2(2), cell, n_potok) = M_K_Moment(2:4, particle_2(2), cell, n_potok) + (t_ex * mu + t2 * mu2) * particle(4:6)
				M_K_Moment(5, particle_2(2), cell, n_potok) = M_K_Moment(5, particle_2(2), cell, n_potok) + &
					(t_ex * mu + t2 * mu2) * kvv(particle(4), particle(5), particle(6))
				
				if(par_n_moment > 9) then
					M_K_Moment(10, particle_2(2), cell, n_potok) = M_K_Moment(10, particle_2(2), cell, n_potok) + &
						(t_ex * mu + t2 * mu2) * particle(4)**2
					M_K_Moment(11, particle_2(2), cell, n_potok) = M_K_Moment(11, particle_2(2), cell, n_potok) + &
						(t_ex * mu + t2 * mu2) * particle(4)*particle(5)
					M_K_Moment(12, particle_2(2), cell, n_potok) = M_K_Moment(12, particle_2(2), cell, n_potok) + &
						(t_ex * mu + t2 * mu2) * particle(4)*particle(6)
					M_K_Moment(13, particle_2(2), cell, n_potok) = M_K_Moment(13, particle_2(2), cell, n_potok) + &
						(t_ex * mu + t2 * mu2) * particle(5)*particle(5)
					M_K_Moment(14, particle_2(2), cell, n_potok) = M_K_Moment(14, particle_2(2), cell, n_potok) + &
						(t_ex * mu + t2 * mu2) * particle(5)*particle(6)
					M_K_Moment(15, particle_2(2), cell, n_potok) = M_K_Moment(15, particle_2(2), cell, n_potok) + &
						(t_ex * mu + t2 * mu2) * particle(6)*particle(6)
					M_K_Moment(16, particle_2(2), cell, n_potok) = M_K_Moment(16, particle_2(2), cell, n_potok) + &
						(t_ex * mu + t2 * mu2) * particle(4)**3
					M_K_Moment(17, particle_2(2), cell, n_potok) = M_K_Moment(17, particle_2(2), cell, n_potok) + &
						(t_ex * mu + t2 * mu2) * particle(5)**3
					M_K_Moment(18, particle_2(2), cell, n_potok) = M_K_Moment(18, particle_2(2), cell, n_potok) + &
						(t_ex * mu + t2 * mu2) * particle(6)**3
				end if
				
				
				if (u / cp > 7.0) then
					uz_M = MK_Velosity_2(u, cp)/ (uz * (cp**2) * cp * par_pi_8 * par_sqrtpi)
					uz_E = MK_Velosity_3(u, cp)
					
					if(MK_is_NaN == .True. .and. ieee_is_nan(uz_E)) then
						print*, "NaN 34rtyu765rdfgh"
						pause
					end if
					
					M_K_Moment(6, particle_2(2), cell, n_potok) = M_K_Moment(6, particle_2(2), cell, n_potok) - mu_ex * uz_M * u1 / u
					M_K_Moment(7, particle_2(2), cell, n_potok) = M_K_Moment(7, particle_2(2), cell, n_potok) - mu_ex * uz_M * u2 / u
					M_K_Moment(8, particle_2(2), cell, n_potok) = M_K_Moment(8, particle_2(2), cell, n_potok) - mu_ex * uz_M * u3 / u
					M_K_Moment(9, particle_2(2), cell, n_potok) = M_K_Moment(9, particle_2(2), cell, n_potok) + &
						mu_ex * (-0.25 * (3.0 * cp**2 + 2.0 * u**2) * (uz_E / uz) - uz_M * skalar / u)
				else
					k1 = MK_int_1(u, cp)
					k2 = MK_int_2(u, cp)
					k3 = MK_int_3(u, cp)
					
					if( MK_is_NaN == .True. .and. (ieee_is_nan(k1) .or. ieee_is_nan(k2) .or. ieee_is_nan(k3))) then
						print*, "NaN klkjhgfddfujyukjhjhgfghozpem"
						pause
					end if
					
					M_K_Moment(6, particle_2(2), cell, n_potok) = M_K_Moment(6, particle_2(2), cell, n_potok) + mu_ex * (k2/k1) * u1 / u
					M_K_Moment(7, particle_2(2), cell, n_potok) = M_K_Moment(7, particle_2(2), cell, n_potok) + mu_ex * (k2/k1) * u2 / u
					M_K_Moment(8, particle_2(2), cell, n_potok) = M_K_Moment(8, particle_2(2), cell, n_potok) + mu_ex * (k2/k1) * u3 / u
					M_K_Moment(9, particle_2(2), cell, n_potok) = M_K_Moment(9, particle_2(2), cell, n_potok) + &
						mu_ex * (-0.5 * k3/k1 + k2/k1 * skalar / u)
					
				end if

				! Добавим интегралы по перезарядке на пикапах !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
				if(area2 <= 2 .and. par_PUI == .True.) then
					if(f_pui_num2(cell2) < 1) then
						print*, "(__________________________)"
						print*, cell2
						print*, cell
						print*, f_pui_num2(cell2)
						print*, area2
						pause "Error 1344 MK"
					end if

					k1 = PUI_get_nu_integr(f_pui_num2(cell2), u)
					uz_M = PUI_get_Mz_integr(f_pui_num2(cell2), u)/k1
					uz_E = PUI_get_E_integr(f_pui_num2(cell2), u)/k1
					uz_M = 1.0 - uz_M/u
					k1 = (particle(4) - vx) * uz_M
					k2 = (particle(5) - vy) * uz_M
					k3 = (particle(6) - vz) * uz_M
					M_K_Moment(6, particle_2(2), cell, n_potok) = M_K_Moment(6, particle_2(2), cell, n_potok) + k1 * mu_ex_pui
					M_K_Moment(7, particle_2(2), cell, n_potok) = M_K_Moment(7, particle_2(2), cell, n_potok) + k2 * mu_ex_pui
					M_K_Moment(8, particle_2(2), cell, n_potok) = M_K_Moment(8, particle_2(2), cell, n_potok) + k3 * mu_ex_pui
					M_K_Moment(9, particle_2(2), cell, n_potok) = M_K_Moment(9, particle_2(2), cell, n_potok) + mu_ex_pui * (u**2 / 2.0 + &
					(k1 * vx + k2 * vy + k3 * vz) - uz_E)
				end if

				! Добавим интегралы по фотоионизации  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				M_K_Moment(19, particle_2(2), cell, n_potok) = M_K_Moment(19, particle_2(2), cell, n_potok) + mu_ph 
				M_K_Moment(6:8, particle_2(2), cell, n_potok) = M_K_Moment(6:8, particle_2(2), cell, n_potok) + mu_ph * particle(4:6)
				M_K_Moment(9, particle_2(2), cell, n_potok) = M_K_Moment(9, particle_2(2), cell, n_potok) + &
					mu_ph * (0.5 * norm2(particle(4:6)) + par_E_ph)
				
				!! Добавляем для расчёта PUI
				if(area2 <= 2 .and. par_PUI == .True.) then
					!? call PUI_Add(cell, u, kappa/time, nu_ph, mu_ex, mu_ph, t_ex, time)
					!call PUI_Add(cell, u, kappa/time, nu_ph, mu, t_ex)  ! (kappa_all - kappa_ph)/time
					!call PUI_Add(cell, u, kappa/time, nu_ph, mu2, t2)
					
					if(area2 == 1 .and. particle_2(2) /= 1 .and. particle_2(2) /= 5) then  ! От атомов этого-же сорта и тепловой плазмы не рождаются пикапы
						call PUI_Add(cell, u, (kappa_all - kappa_ph)/time, mu, t_ex)
						call PUI_Add(cell, u, (kappa_all - kappa_ph)/time, mu2, t2)
					else
						call PUI_Add(cell, u, 0.0_8, mu, t_ex)
						call PUI_Add(cell, u, 0.0_8, mu2, t2)
					end if
				end if
				!_________________________________________________________________________________________
				
				!! Добавляем расчёт поглощения  ---------------------------------------------------
				if(par_pogloshenie == .True.) then
					uz_M = DOT_PRODUCT(particle(4:6), r_ex)/r   ! Проекция скорости на ось зрения (от звёзды)
					if(uz_M > pogl_v_min .and. uz_M < pogl_v_max) then
						ijk = MK_poglosh_nomer(uz_M)
						call omp_set_lock(pogl_lock(cell))
						pogloshenie(particle_2(2), ijk, cell) = pogloshenie(particle_2(2), ijk, cell) + t_ex * mu + t2 * mu2
						call omp_unset_lock(pogl_lock(cell))
					end if
				end if
				!! --------------------------------------------------------------------------------
				
				
				call spherical_skorost(r_ex(1), r_ex(2), r_ex(3), vx, vy, vz, Ur, Uphi, Uthe)
				call spherical_skorost(r_ex(1), r_ex(2), r_ex(3), particle(4), particle(5), particle(6), Vr, Vphi, Vthe)
			
				! Перезаряжаем
				if (area2 == 1 .or. Ur / cp > 1.8) then   ! Без геометрического расщепления
					II = 0  ! II - это сколько дополнительных атомов запускается (помимо основного)
				else
					II = MK_geo_zones(r, 1.2_8) - 1
				end if

				if(area2 <= 2 .and. par_PUI == .True.) then  !! Нужно перезаряжать с пикапами
					call MK_pui_charge_exchange_velocity(n_potok, vx, vy, vz, particle(4), particle(5), particle(6), f_pui_num2(cell2), k1, k2, k3)
					! Функция возвращает сразу декартовые компоненты новой скорости
					V = (/ k1, k2, k3 /)
					call MK_Distination(r_ex, V, to_i, to_j, r_peregel)
					mu3 = mu_ex_pui
					call MK_ruletka(n_potok, to_i, to_j, from_i, from_j, area2, r, r_peregel, mu3, bb2)

					if(area2 == 1) new_sort = 5
					if(area2 == 2) new_sort = 6

					if(bb2 /= .False.) then  ! Если не надо вырубать атом
						! Запишем новую частицу в стек
						stek(n_potok) = stek(n_potok) + 1
						M_K_particle(1:3, stek(n_potok), n_potok) = r_ex
						M_K_particle(4:6, stek(n_potok), n_potok) = V
						M_K_particle(7, stek(n_potok), n_potok) = mu3
						M_K_particle(8, stek(n_potok), n_potok) = r_peregel

						M_K_particle_2(:,stek(n_potok), n_potok) = (/ cell, new_sort, to_i, to_j /)
						
						if(MK_Mu_stat == .True.) then  !!!  ??????????????????????????????????????? !! !
							if(particle_2(2) == new_sort) then
								M_K_particle_3(:, :, stek(n_potok), n_potok) = particle_3
							else 
								M_K_particle_3(:, :, stek(n_potok), n_potok) = .False.
							end if
						end if
						
					end if
				end if
				
				call M_K_Change_Velosity4(n_potok, Ur/cp, Uthe/cp, Uphi/cp, Vr/cp, Vthe/cp, Vphi/cp, Wr, Wthe, Wphi, mu_, &
					cp, r, II, r_ex(1), r_ex(2), r_ex(3), bb)
				
				Wr = Wr * cp
				Wthe = Wthe * cp
				Wphi = Wphi * cp
				
				do i = 1, II + 1
					if(i == II + 1 .and. bb == .False.) CYCLE  ! Если не запускается основной атом
					call dekard_skorost(r_ex(1), r_ex(2), r_ex(3), Wr(i), Wphi(i), Wthe(i), v1, v2, v3)
					V = (/ v1, v2, v3 /)
					
					call MK_Distination(r_ex, V, to_i, to_j, r_peregel)  ! Находим зону назначения в точке перегелия для рулетки
					
					mu3 = mu_ex * mu_(i)
					call MK_ruletka(n_potok, to_i, to_j, from_i, from_j, area2, r, r_peregel, mu3, bb2)
					
					if(bb2 == .False.) CYCLE  ! Не запускаем эту частицу, она вырубается
					
					! Запишем новую частицу в стек
					stek(n_potok) = stek(n_potok) + 1
					M_K_particle(1:3, stek(n_potok), n_potok) = r_ex
					M_K_particle(4:6, stek(n_potok), n_potok) = V
					M_K_particle(7, stek(n_potok), n_potok) = mu3
					M_K_particle(8, stek(n_potok), n_potok) = r_peregel

					M_K_particle_2(:,stek(n_potok), n_potok) = (/ cell, area2, to_i, to_j /)
					
					if(MK_Mu_stat == .True.) then
						if(particle_2(2) == area2) then
							M_K_particle_3(:, :, stek(n_potok), n_potok) = particle_3
						else 
							M_K_particle_3(:, :, stek(n_potok), n_potok) = .False.
						end if
					end if
					
				end do
				
				! Проверяем, нужно ли вырубить неперезаряженную часть атома
				call MK_ruletka(n_potok, particle_2(3), particle_2(4), from_i, from_j, particle_2(2), &
						r, particle(8), mu2, bb2)
				if(bb2 == .False.) EXIT loop1  ! Не запускаем эту частицу, она вырубается
				particle(7) = mu2
			
				! Находим следующую ячейку
				
				particle(1:3) = particle(1:3) + time * particle(4:6)
				if(next == 0) then
					if (norm2(particle(1:3)) <= 100.0) then
						print*, "Error 23e2323 ", particle(1:3)
						pause
					end if
					EXIT loop1  ! частица долетела до края области
				end if
				
		11			continue 
				
				
				
				
				call Int2_Get_tetraedron(particle(1), particle(2), particle(3), next)  ! Находим точно следующий тетраэдр
				if (next == -1) then ! ЗАТЫК
					next = 3
					particle(1:3) = particle(1:3) + (time/100.0) * particle(4:6)
					print*, "Zatik"
					pause
					GO TO 11
				end if
				
				particle_2(1) = next
				if(next == 0) then
					if (norm2(particle(1:3)) <= 100.0) then
						print*, "Error wqer2r42  ", particle(1:3)
						pause
					end if
					
					EXIT loop1  ! частица долетела до края области
				end if
				
				if(particle(1) > 0.00001 .and. norm2(particle(1:3)) > par_Rmax + 0.001) then
					EXIT loop1  ! частица долетела до края области
				end if
				
			
				
			
			end do loop1
			
			
		end do
	
	end subroutine M_K_Fly_PUI

	subroutine M_K_Calc_Pogloshenie(rx, ry, rz, nomer)
        ! Печатаем поверхности, которые выделяем
        integer, intent(in) :: nomer
        real(8), intent(in) :: rx, ry, rz  !! Направление луча зрения (единичный вектор!)
		real(8) :: r(3), dl
		real(8) :: pogl(par_n_sort, pogl_iter), u
		character(len=5) :: name
		integer(4) :: cell, sort, i

		print*, "Start M_K_Calc_Pogloshenie ", rx, ry, rz, nomer
		if(par_poglosh == .False.) return

		r(1) = rx * par_R0/10.0
		r(2) = ry * par_R0/10.0
		r(3) = rz * par_R0/10.0
        dl = par_R0/2.0
		pogl = 0.0
		write(unit=name,fmt='(i5.5)') nomer
		cell = 3

		do while (.True.)
			r = r + dl * (/ rx, ry, rz /)

			call Int2_Get_tetraedron(r(1), r(2), r(3), cell)
			if(cell < 1) EXIT

			do sort = 1, par_n_sort
				do i = 1, pogl_iter
					pogl(sort, i) = pogl(sort, i) + par_n_p_LISM * dl * pogloshenie(sort, i, cell)/pogl_ddd
				end do
			end do

		end do

		open(1, file = "Pogloshenie_" // name // ".txt")
		if(par_n_sort == 6) write(1,*) "TITLE = 'HP'  VARIABLES = u, f1, f2, f3, f4, f5, f6, ff"
		if(par_n_sort == 4) write(1,*) "TITLE = 'HP'  VARIABLES = u, f1, f2, f3, f4, ff"
		
		pogl = pogl * par_poglosh

		print*, "Size poglosh = ", size(pogloshenie(:, 1, 1)), size(pogloshenie(1, :, 1)), size(pogloshenie(1, 1, :))

		do i = 1, pogl_iter
			u = pogl_v_min + (i + 0.5) * pogl_ddd;
			if(par_n_sort == 6) write(1,*) u, exp(-pogl(1, i)), exp(-pogl(2, i)), exp(-pogl(3, i)), exp(-pogl(4, i)), exp(-pogl(5, i)), exp(-pogl(6, i)), exp(-sum(pogl(:, i)))
			if(par_n_sort == 4) write(1,*) u, exp(-pogl(1, i)), exp(-pogl(2, i)), exp(-pogl(3, i)), exp(-pogl(4, i)), exp(-sum(pogl(:, i)))
		end do

		close(1)
    end subroutine M_K_Calc_Pogloshenie

	integer(4) pure function MK_poglosh_nomer(V)
	    !integer(4) function MK_poglosh_nomer(V)
		!! В какой номер в массиве поглощений нужно записать данную скорость
		real(8), intent(in) :: V

		MK_poglosh_nomer = max(min(INT((V - pogl_v_min)/pogl_ddd), pogl_iter), 1)

		! print*, pogl_v_min + (MK_poglosh_nomer + 0.5) * pogl_ddd, V
		! print*, "VVV", pogl_v_min + (MK_poglosh_nomer - 1 + 0.5) * pogl_ddd, pogl_v_min + (MK_poglosh_nomer + 0.5) * pogl_ddd, pogl_v_min + (MK_poglosh_nomer + 1 + 0.5) * pogl_ddd
		! print*, "dv", V - (pogl_v_min + (MK_poglosh_nomer - 1 + 0.5) * pogl_ddd), V - (pogl_v_min + (MK_poglosh_nomer + 0.5) * pogl_ddd), V - (pogl_v_min + (MK_poglosh_nomer + 1 + 0.5) * pogl_ddd)
		! pause
		return
	end function MK_poglosh_nomer
	
	subroutine M_K_rand(s1, s2, s3, b)
		integer(4), intent(in out) :: s1, s2, s3
		real(8), intent(out) :: b
		integer(4):: ic15 = 32768, ic10 = 1024
		integer(4):: mz = 710, my = 17784, mx = 11973
		real(8):: xi = 9.0949470177292824E-13, c = 1.073741824E9
		integer(4) :: i13, i12, i11, ii
		i13 = mz * s1 + my * s2 + mx * s3
		i12 = my * s1 + mx * s2
		i11 = mx * s1
		ii = i11 / ic15
		i12 = i12 + ii
		s1 = i11 - ic15 * ii
		ii = i12 / ic15
		i13 = i13 + ii
		s2 = i12 - ic15 * ii
		s3 = mod(i13,ic10)
		b = xi * (c * s3 + ic15 * s2 + s1)
	end subroutine M_K_rand
	
	subroutine M_K_Set()
	
		integer(4) :: i, j, n, k, sort, m_zone, n_zone
		real(8) :: Yr
		logical :: exists

		print*, "M_K_Set(), sortov = ", par_n_sort
		
		allocate(M_K_particle(8, par_stek, par_n_potok))
		allocate(M_K_particle_2(4, par_stek, par_n_potok))
		
		if(MK_Mu_stat) then
			allocate(M_K_particle_3(par_n_zone + 1, par_m_zone + 1, par_stek, par_n_potok))
			allocate(MK_Mu_statistic(par_n_zone + 1, par_m_zone + 1, par_n_sort))
			M_K_particle_3 = .False.
			MK_Mu_statistic = 0.0
		end if
		
		allocate(sensor(3, 2, par_n_potok))
		allocate(stek(par_n_potok))
		allocate(MK_Mu(par_n_zone + 1, par_m_zone + 1, par_n_sort))
		allocate(M_K_Moment(par_n_moment, par_n_sort, size(int2_all_tetraendron(1, :)), par_n_potok))
		
	
		MK_Mu = 1.0
		M_K_particle = 0.0
		M_K_particle_2 = 0
		stek = 0
		M_K_Moment = 0.0
		sensor = 1
		
		! Задаём радиусы зон
		!MK_R_zone(1) = 1.0
		MK_R_zone(1) = 2.0
		MK_R_zone(2) = 4.6
		MK_R_zone(3) = 10.0
		MK_R_zone(4) = 25.0
		MK_R_zone(5) = 60.0
		MK_R_zone(6) = 130.0

		
		
		
		! Задаём лучи зон
		do i = 1, par_m_zone
			MK_al_zone(i) = i * par_pi_8/(par_m_zone + 1)
		end do
		
		! Задаём критические веса
		
		do j = 1, par_m_zone + 1
			do i = 1, par_n_zone
				MK_Mu(i, j, :) = 1.0 !(MK_R_zone(i)/par_Rmax)**(1.3_8)
			end do
		end do
		

		if(.False.) then !! Старый файл с весами
			inquire(file="MK_Mu_statistic.txt", exist=exists)
			if (exists == .False.) then
				pause "net faila!!!  cvbdfgertmkopl"
				STOP "net faila!!!"
			end if
			open(2, file = "MK_Mu_statistic.txt", status = 'old')
			
			do k = 1, 4 !!par_n_sort
				do j = 1, par_m_zone + 1
					do i = 1, par_n_zone + 1
						read(2, *) n, n, sort, MK_Mu(i, j, k)
					end do
				end do
			end do
			
			close(2)

			if(par_n_sort == 6) then
				MK_Mu(:, :, 5) = MK_Mu(:, :, 1) !!
				MK_Mu(:, :, 6) = MK_Mu(:, :, 2) !!
			end if
		else
			inquire(file="MK_Mu_statistic2.txt", exist=exists)
			if (exists == .False.) then
				pause "net faila!!!  WCFEWCR2C24R2R235tgwgv3rw3vr4qc2- c"
				STOP "net faila!!!"
			end if
			open(2, file = "MK_Mu_statistic2.txt", status = 'old')
			
			read(2, *) n_zone, m_zone, sort

			if(par_n_sort > sort) then
				print*, "ERROR iji9j9uhcw8ogbor3u4buyv4cwl"
				STOP
			end if

			do k = 1, par_n_sort
				do j = 1, m_zone! par_m_zone + 1
					do i = 1, n_zone!par_n_zone + 1
						read(2, *) n, n, n, MK_Mu(i, j, k)
					end do
				end do
			end do
			
			close(2)

			

			if(par_m_zone + 1 /= m_zone) then
				print*, "ERROR berebt5verby"
				STOP
			end if

			if(par_n_zone + 1 /= n_zone) then
				print*, "ERROR er5buim8ominubyvtc4t5yb6nui6"
				STOP
			end if
		end if
		
		
		MK_Mu(:, :, 2) = MK_Mu(:, :, 2) * 0.8 ! * 0.2
		MK_Mu(:, :, 3) = MK_Mu(:, :, 3) * 2.0! * 0.5
		MK_Mu(:, :, 4) = MK_Mu(:, :, 4) * 2.0! * 0.5

		if(par_n_sort >= 6) then
			MK_Mu(:, :, 6) = MK_Mu(:, :, 6) * 0.5! * 0.5
		end if
		
		
		MK_SINKR(1) = sin(MK_al_zone(1))
		MK_SINKR(par_m_zone + 1) = sin(MK_al_zone(par_m_zone))
		
		do i = 2, par_m_zone
			MK_SINKR(i) = max( dabs(sin(MK_al_zone(i))), dabs(sin(MK_al_zone(i - 1))) )
		end do
		
		do i = 1, par_n_zone
			MK_gam_zone(i) = 1.0 / ((par_Rmax / MK_R_zone(i))**2 - 1.0)
			if (MK_gam_zone(i) < 0.0) pause "ERROR gamma 56789oihgfr6uijt6789"
		end do
		
		Yr = dabs(par_Velosity_inf)
		MK_A0_ = (Yr + 1.0 / (2.0 * Yr)) * erf(Yr) + exp(-(Yr**2)) / par_sqrtpi + Yr
		MK_A1_ = 1.0 + (1.0 + 1.0 / (2.0 * (Yr)**2 )) * erf(Yr) + exp(-(Yr**2) ) / (par_sqrtpi * Yr) 


		!! Создаём массивы для учёта поглощения
		print*, "1 !!!!!!!!!!! - Size poglosh = ", allocated(pogloshenie)

		if(par_pogloshenie == .True. .and. allocated(pogloshenie) == .False.) then
			print*, "Start Set par_pogloshenie ", par_n_sort, pogl_iter, size(int2_all_tetraendron_point(1, :))
			ALLOCATE(pogloshenie(par_n_sort, pogl_iter, size(int2_all_tetraendron_point(1, :))))
			print*, "ALLOCATE"
			pogl_ddd = (pogl_v_max - pogl_v_min)/pogl_iter
			pogloshenie = 0.0
			print*, "End Set par_pogloshenie"
		else
			pogl_ddd = (pogl_v_max - pogl_v_min)/pogl_iter
			pogloshenie = 0.0
		end if

		print*, "!!!!!!!!!!! - Size poglosh = ", size(pogloshenie(:, 1, 1)), size(pogloshenie(1, :, 1)), size(pogloshenie(1, 1, :))
		print*, "2 !!!!!!!!!!! - Size poglosh = ", allocated(pogloshenie)

		if(par_pogloshenie == .True. .and. allocated(pogl_lock) == .False.) then
			ALLOCATE(pogl_lock(size(int2_all_tetraendron_point(1, :))))

			do i = 1, size(int2_all_tetraendron_point(1, :))
				call omp_init_lock(pogl_lock(i))
			end do
		end if
		
	
	end subroutine M_K_Set

	subroutine MK_pui_charge_exchange_velocity(potok, Upx, Upy, Upz, UHx, UHy, UHz, num, VHx, VHy, VHz)
		!? Функция перезарядки - получает новые скорости атома после перезарядки
		USE PUI
		implicit none
		real(8), intent(in) :: Upx, Upy, Upz, UHx, UHy, UHz
		integer, intent(in) :: potok, num ! В какой ячейке сейчас находимся (чтобы брать эту f_pui)
		real(8), intent(out) :: VHx, VHy, VHz
		
		real(8) :: ksi1, ksi2, ksi3
		real(8) :: UH, w, the, u, h0, phi
		real(8) :: ex(3), ey(3), ez(3), vx, vy, vz

		ez(1) = UHx - Upx
		ez(2) = UHy - Upy
		ez(3) = UHz - Upz

		UH = sqrt((Upx - UHx)**2 + (Upy - UHy)**2 + (Upz - UHz)**2)
		ez = ez/UH
		h0 = PUI_get_h0(UH)

		call get_bazis(ez, ex, ey)

		do while(.True.)
			call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi1)
			call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi2)
			call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi3)
			w = PUI_get_F_integer(ksi1, num)
			the = acos(1.0 - 2.0 * ksi2)
			u = sqrt(w**2 * sin(the)**2 + (w * cos(the) - UH)**2)
			if(u * MK_sigma(u)/( (w + pui_h0_wc) * MK_sigma(w + pui_h0_wc) * h0) >= ksi3) EXIT
		end do

		call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi1)
		phi = ksi1 * 2.0 * par_pi_8

		vx = w * sin(the) * cos(phi)
		vy = w * sin(the) * sin(phi)
		vz = w * cos(the)

		VHx = Upx + vx * ex(1) + vy * ey(1) + vz * ez(1)
		VHy = Upy + vx * ex(2) + vy * ey(2) + vz * ez(2)
		VHz = Upz + vx * ex(3) + vy * ey(3) + vz * ez(3)
	end subroutine MK_pui_charge_exchange_velocity
	
	subroutine M_K_Change_Velosity4(potok, Ur, Uthe, Uphi, Vr, Vthe, Vphi, Wr_, Wthe_, Wphi_, mu_, cp, r, I_, x_ex, y_ex, z_ex, bb)
		! Как первая часть, но розыгрышь идёт по-частям
		! Запускаесть один основной и I_ дополнительных атомов
		
		integer(4), intent(in) :: potok, I_
		real(8), intent(in) ::  Ur, Uthe, Uphi, Vr, Vthe, Vphi, cp, r, x_ex, y_ex, z_ex
		real(8), intent(out) ::   Wr_(I_ + 1), Wthe_(I_ + 1), Wphi_(I_ + 1), mu_(I_ + 1)
		logical, intent(out) :: bb
		
		real(8) :: X, uu
		real(8) :: gamma_(I_), Wa_(I_), Mho_(I_)
		real(8) :: ksi, gam1, gam2, Wr1, Wr2, Wr0, ksi1, ksi2, W1, W2, Wa, pp1, pp2, c, p, u
		real(8) :: p4, om1, om2, om3, lo, y1, y2, y3, v1, v2, v3, u1, u2, u3, uuu, yy, h, ksi3, ksi4, ksi5, ksi6, D, ko, gg
		integer(4) :: ii, met, k
		
		X = sqrt( (Vr - Ur)**2 + (Vthe - Uthe)**2 + (Vphi - Uphi)**2 )
		uu = exp(-(X**2)) / par_sqrtpi + (X + 1.0 / (2.0 * X)) * erf(X)
		Wr0 = -1.0
		met = 0

		do ii = 1, I_
			gamma_(ii) = 1.0 / ( (r / MK_R_zone(ii))**2 - 1.0)
		end do

		

		! Разыграем Wr
		
		do ii = 1, I_
			
			if (ii == 1) then
				gam1 = 0.0
				gam2 = gamma_(ii)
			else
				gam1 = gamma_(ii - 1)
				gam2 = gamma_(ii)
			end if
			
			
			
			pp1 = MK_for_Wr_1(0.0_8, gam2, Ur) - MK_for_Wr_1(0.0_8, gam1, Ur)
			pp2 = MK_for_Wr_2(0.0_8, gam2, Ur, Uthe) - MK_for_Wr_2(0.0_8, gam1, Ur, Uthe)
			
			
			call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi)

			if (ksi < pp1 / (pp1 + pp2)) then
				call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi)
				Wr1 = -3.0;
				Wr2 = 0.0;

				do while (MK_H_Wr_1(gam1, gam2, Wr1, Ur, pp1, ksi) >= 0.0)
					Wr1 = Wr1 - 1.0
				end do
					
				k = 0
				do while (dabs(Wr2 - Wr1) > 0.0001)     ! Деление пополам, иначе разваливается
					Wr0 = (Wr1 + Wr2) / 2.0
					if (MK_H_Wr_1(gam1, gam2, Wr0, Ur, pp1, ksi) < 0) then
						Wr1 = Wr0
					else
						Wr2 = Wr0
					end if
					k = k + 1
				end do
				Wr_(ii) = Wr1
			else
				call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi)
				Wr1 = -3.0
				Wr2 = 0.0

				do while (MK_H_Wr_2(gam1, gam2, Wr1, Ur, Uthe, pp2, ksi) >= 0.0)
					Wr1 = Wr1 - 1.0
				end do
				
				k = 0
				do while (dabs(Wr2 - Wr1) > 0.0001)     ! Деление пополам, иначе разваливается
					Wr0 = (Wr1 + Wr2) / 2.0
					if (MK_H_Wr_2(gam1, gam2, Wr0, Ur, Uthe, pp2, ksi) < 0) then
						Wr1 = Wr0
					else
						Wr2 = Wr0
					end if
					k = k + 1
				end do
				Wr_(ii) = Wr1
			end if
		end do

		
		! Разыгрываем  Wa
		do ii = 1, I_
			
			if (ii == 1) then
				gam1 = 0.0
				gam2 = gamma_(ii)
			else
				gam1 = gamma_(ii - 1)
				gam2 = gamma_(ii)
			end if

			do while(.True.)
				call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi1)
				call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi2)
				W1 = sqrt(gam1 * (Wr_(ii)**2))
				W2 = sqrt(gam2 * (Wr_(ii)**2))
				Wa = sqrt(-log(exp(-(W1**2)) - ksi1 * (exp(-(W1**2)) - exp(-(W2**2)))))
				if ((1.0 + (Uthe * Wa)**2) / (1.0 + (Uthe * W2)**2) >= ksi2) EXIT
			end do
			Wa_(ii) = Wa
		end do

		
		! Считаем веса и Mho
		do ii = 1, I_
			c = Uthe * Wa_(ii)
			p = MK_norm_mho(c)
			
			Mho_(ii) = MK_play_mho(potok, c)
			

			Wthe_(ii) = Wa_(ii) * cos(Mho_(ii))
			Wphi_(ii) = Wa_(ii) * sin(Mho_(ii))
			
			if (ii == 1) then
				gam1 = 0.0
				gam2 = gamma_(ii)
			else
				gam1 = gamma_(ii - 1)
				gam2 = gamma_(ii)
			end if

			u = sqrt( (Vr - Wr_(ii))**2 + (Vthe - Wthe_(ii))**2 +(Vphi - Wphi_(ii))**2 )

			if (X > 7.0) then
				mu_(ii) = (u * MK_sigma2(u, cp) / (uu * MK_sigma2(uu, cp))) * (MK_f2(0.0_8, gam1, Ur, Uthe) - MK_f2(0.0_8, gam2, Ur, Uthe)) * exp(-(Uthe**2)) * (p) / (1.0 + (c**2))
			else
				mu_(ii) = (u * MK_sigma2(u, cp) / (MK_int_1(X * cp, cp) / cp)) * (MK_f2(0.0_8, gam1, Ur, Uthe) - MK_f2(0.0_8, gam2, Ur, Uthe)) * exp(-(Uthe**2)) * (p) / (1.0 + (c**2))
			end if
		end do


		! Розыгрыш основного атома
		p4 = 0.5 * par_sqrtpi * X / (1.0 + 0.5 * par_sqrtpi * X)
		

		if (I_ > 0) then
			gg = gamma_(I_)
		else
			gg = 0.0
		end if

		do while(.True.)
			call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi1)
			call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi2)
			call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi3)
			call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi4)
			call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi5)
			call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi6)
			
			if (p4 < ksi1) then
				om1 = 1.0 - 2.0 * ksi4
				om2 = sqrt(1.0 - (om1**2)) * cos(2.0 * par_pi_8 * ksi5)
				om3 = sqrt(1.0 - (om1**2)) * sin(2.0 * par_pi_8 * ksi5)

				lo = sqrt(-log(ksi2 * ksi3))
				y1 = lo * om1
				y2 = lo * om2
				y3 = lo * om3
			else
				y1 = sqrt(-log(ksi2)) * cos(par_pi_8 * ksi3)
				y2 = sqrt(-log(ksi4)) * cos(2.0 * par_pi_8 * ksi5)
				y3 = sqrt(-log(ksi4)) * sin(2.0 * par_pi_8 * ksi5)
			end if
			
			v1 = y1 + Ur
			v2 = y2 + Uthe
			v3 = y3 + Uphi
			u1 = Vr - v1
			u2 = Vthe - v2
			u3 = Vphi - v3
			uuu = sqrt(kvv(u1, u2, u3))
			yy = sqrt(kvv(y1, y2, y3))
			h = ((uuu * MK_sigma2(uuu, cp)) / (MK_sigma2(X, cp) * (X + yy)))
			
			if (h >= ksi6) EXIT
		end do


		Wr_(I_ + 1) = v1
		Wthe_(I_ + 1) = v2
		Wphi_(I_ + 1) = v3


		if (Wr_(I_ + 1) >= 0.0 .or. (Wthe_(I_ + 1)**2) + (Wphi_(I_ + 1)**2) > gg * (Wr_(I_ + 1)**2)) then
			mu_(I_ + 1) = 1.0
			bb = .True.
			return
		else
			mu_(I_ + 1) = 0.0  ! Чтобы не запускать этот атом
			bb = .False.
			return
		end if

		bb = .True.
		return
	end subroutine M_K_Change_Velosity4

	subroutine MK_Velosity_initial2(potok, Vx, Vy, Vz)
		
		integer(4), intent(in) :: potok
		real(8), intent(out) :: Vx, Vy, Vz
		
		real(8) :: ksi1, ksi2, a, ksi3, ksi4, ksi5, ksi6, z, p1
		
		call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi1)
		call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi2)
		
		a = sqrt(-log(1.0 - ksi2))
		Vy = a * cos(2.0 * par_pi_8 * ksi1)
		Vz = a * sin(2.0 * par_pi_8 * ksi1)
		
		z = 0
		p1 = 0.5 * dabs(par_Velosity_inf) * par_sqrtpi / (0.5 + 0.5 * dabs(par_Velosity_inf) * par_sqrtpi);

		do
			call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi3)
			call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi4)
			call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi5)
			call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi6)

			if (p1 > ksi3) then
				z = cos(par_pi_8 * ksi5) * sqrt(-log(ksi4))
			else
				z = sqrt(-log(1.0 - ksi4))
			end if
			
			if((dabs(z + par_Velosity_inf) / (dabs(par_Velosity_inf) + dabs(z)) > ksi6 .and. z >= -par_Velosity_inf)) EXIT
		end do

		Vx = z + par_Velosity_inf
		
		if (Vx <= 0.0) then
			print*, "Error iuygvbnmklo9890pljiouytrtyjhg"
		end if
		
		return
	end subroutine MK_Velosity_initial2
	
	subroutine MK_Velosity_initial(potok, Vx, Vy, Vz)
		
		integer(4), intent(in) :: potok
		real(8), intent(out) :: Vx, Vy, Vz
		
		real(8) :: ksi1, ksi2, a, ksi3, ksi4, ksi5, ksi6, z, p1
		
		call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi1)
		call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi2)
		
		a = sqrt(-log(ksi2))
		Vy = a * cos(2.0 * par_pi_8 * ksi1)
		Vz = a * sin(2.0 * par_pi_8 * ksi1)
		
		z = 0
		p1 = dabs(par_Velosity_inf) * par_sqrtpi / (1.0 + dabs(par_Velosity_inf) * par_sqrtpi)

		do
			call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi3)
			call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi4)
			call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi5)
			call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi6)

			if (p1 > ksi3) then
				z = cos(par_pi_8 * ksi5) * sqrt(-log(ksi4))
			else
				if (ksi4 <= 0.5) then
					z = -sqrt(-log(2.0 * ksi4))
				else
					z = sqrt(-log(2.0 * (1.0 - ksi4)))
				end if
			end if
			
			if(dabs(z + par_Velosity_inf) / (dabs(par_Velosity_inf) + dabs(z)) >= ksi6 .and. z <= -par_Velosity_inf) EXIT
		end do

		Vx = z + par_Velosity_inf
		
		return
	end subroutine MK_Velosity_initial
	
	subroutine MK_Init_Parametrs(potok, mu_, Wt_, Wp_, Wr_, X_, bb)
		! Розыгрышь атомов на начальной сфере (скорости и положения)
		integer(4), intent(in) :: potok
		real(8), intent(out) :: mu_(par_n_zone + 1), Wt_(par_n_zone + 1), Wp_(par_n_zone + 1), Wr_(par_n_zone + 1), X_(par_n_zone + 1)
		logical, intent(out) :: bb
		
		real(8) :: Y, Wr1, Wr2, Wr0, W1, W2, Wa, Yt, c, p
		real(8) :: ksi, ksi1, ksi2
		real(8) :: X1
		real(8) :: X2
		real(8) :: X0                ! для метода хорд
		real(8) :: split, gam1, gam2
		real(8) :: Wa_(par_n_zone + 1)
		real(8) :: Mho_(par_n_zone + 1)
		integer(4) :: i
		real(8) :: p1, ksi3, ksi4, z, h, gg, p4_, ksi5, ksi6
		
		
		X0 = 1.0
		Y = dabs(par_Velosity_inf)
		p1 = erf(Y) / (MK_A1_ * (Y**2))
		
		

		! Разыгрываем  X
		do i = 1, par_n_zone
			
			call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi)
			X1 = 0.0
			X2 = 1.0
			
			if (i == 1) then
				gam1 = 0.0
				gam2 = MK_gam_zone(i)
			else
				gam1 = MK_gam_zone(i - 1)
				gam2 = MK_gam_zone(i)
			end if
			

			do while (dabs(X2 - X1) > 0.0000001)     ! Деление пополам, иначе разваливается
				X0 = (X1 + X2) / 2.0
				if (MK_Hx(gam1, gam2, X0, Y, ksi) < 0.0) then
					X1 = X0
				else
					X2 = X0
				end if
			end do
			X_(i) = X0
		end do

		! Разыгрываем  Wr
		do i = 1, par_n_zone
			
			call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi)
			
			Wr1 = -5.0
			Wr2 = 0.0
			
			if (i == 1) then
				gam1 = 0.0
				gam2 = MK_gam_zone(i)
			else
				gam1 = MK_gam_zone(i - 1)
				gam2 = MK_gam_zone(i)
			end if
			
			do while (MK_Hwr(gam1, gam2, Wr1, X_(i), Y, ksi) >= 0.0)
				Wr1 = Wr1 - 1.0
			end do
				
			do while (dabs(Wr2 - Wr1) > 0.000001)     ! Деление пополам, иначе разваливается
				Wr0 = (Wr1 + Wr2) / 2.0
				if (MK_Hwr(gam1, gam2, Wr0, X_(i), Y, ksi) < 0) then
					Wr1 = Wr0
				else
					Wr2 = Wr0
				end if
			end do

			Wr_(i) = Wr0
		end do

		
		! Разыгрываем  Wa
		
		do i = 1, par_n_zone
			Yt = Y * sqrt(1.0 - (X_(i)**2))
			
			if (i == 1) then
				gam1 = 0.0
				gam2 = MK_gam_zone(i)
			else
				gam1 = MK_gam_zone(i - 1)
				gam2 = MK_gam_zone(i)
			end if

			do while(.True.)
				call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi1)
				call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi2)
				W1 = sqrt(gam1 * (Wr_(i)**2))
				W2 = sqrt(gam2 * (Wr_(i)**2))
				Wa = sqrt(-log(exp(-(W1**2)) - ksi1 * (exp(-(W1**2)) - exp(-(W2**2)))))
				if ((1.0 + (Yt * Wa)**2) / (1.0 + (Yt * W2)**2) >= ksi2) EXIT
			end do


			Wa_(i) = Wa
		end do

		! Считаем веса
		do i = 1, par_n_zone
			
			c = Y * sqrt(1.0 - (X_(i)**2)) * Wa_(i)
			p = MK_norm_mho(c)
			
			
			Mho_(i) = MK_play_mho(potok, c)
			

			Wt_(i) = Wa_(i) * cos(Mho_(i))
			Wp_(i) = Wa_(i) * sin(Mho_(i))
			
			if (i == 1) then
				gam1 = 0.0
				gam2 = MK_gam_zone(i)
			else
				gam1 = MK_gam_zone(i - 1)
				gam2 = MK_gam_zone(i)
			end if
			
			mu_(i) = (MK_F(1.0_8, gam2, Y) - MK_F(1.0_8, gam1, Y)) * (p) / (MK_A0_ * (1.0 + (c**2)))
		end do

		! Разыгрываем основной атом

		! Розыгрыш X целиком из препринта

		
		call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi1)
		
		
		if (p1 < ksi1) then
			
			do while(.True.)
				call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi2)
				call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi3)
				X2 = sqrt(ksi2)
				h = (1.0 + erf(X2 * Y)) / (1.0 + erf(Y))
				if (h > ksi3) EXIT
			end do
		else
			do
				call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi2)
				call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi3)
				X2 = (1.0 / Y) * sqrt(-log(ksi2)) * cos(par_pi_8 * ksi3 / 2.0)
				if(X2 < 1.0) EXIT
			end do
		end if
		

		X_(par_n_zone + 1) = X2

		gg = 0.0
		
		if (par_n_zone > 0) gg = MK_gam_zone(par_n_zone)
		
		p4_ = par_sqrtpi * (Y * X2) / (1.0 + par_sqrtpi * (Y * X2))                 ! Для розыгрыша основного атома на границе по препринту Маламы
		
		do while(.True.)
			call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi1)
			call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi2)
			call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi3)
			if (p4_ > ksi1) then
				z = sqrt(-log(ksi2)) * cos(par_pi_8 * ksi3)
			else
				if (ksi2 <= 0.5) then
					z = -sqrt(-log(2.0 * ksi2))
				else
					z = sqrt(-log(2.0 * (1.0 - ksi2)))
				end if
			end if

			Wr_(par_n_zone + 1) = z - (Y * X2)
			h = dabs(-(Y * X2) + z) / ((Y * X2) + dabs(z))
			call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi4)
			if (h > ksi4 .and. z < (Y * X2)) EXIT
		end do

			
			call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi5)
			call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi6)

			Wt_(par_n_zone + 1) = Y * sqrt(1.0 - (X_(par_n_zone + 1)**2)) + sqrt(-log(ksi5)) * cos(2.0 * par_pi_8 * ksi6)
			Wp_(par_n_zone + 1) = sqrt(-log(ksi5)) * sin(2.0 * par_pi_8 * ksi6)


		if (par_n_zone == 0) then
			mu_(par_n_zone + 1) = 1.0
			bb = .True.
			return
		else
			if (Wr_(par_n_zone + 1) >= 0.0 .or. Wt_(par_n_zone + 1)**2 + Wp_(par_n_zone + 1)**2 > &
				MK_gam_zone(par_n_zone) * Wr_(par_n_zone + 1)**2) then
				mu_(par_n_zone + 1) = 1.0
			else
				mu_(par_n_zone + 1) = 0.0  ! Чтобы не запускать этот атом
				bb = .False.
				return
			end if
		end if

		bb = .True.
		return
	end subroutine MK_Init_Parametrs
	
	subroutine MK_Distination(r, V, to_i, to_j, rr)
		! Variables
		real(8), intent(in) :: r(3), V(3)
		real(8), intent(out) :: rr
		integer(4), intent(out) :: to_i, to_j
		
		real(8) :: time_do_peregel, rk(3)
		
		time_do_peregel = -DOT_PRODUCT(r, V) / kvv(V(1), V(2), V(3))
		rk = r + time_do_peregel * V
		rr = norm2(rk)
		
		to_i = MK_geo_zones(rr, 1.0_8)
		to_j = MK_alpha_zones( polar_angle( rk(1), sqrt(rk(2)**2 + rk(3)**2) ) )
	
	end subroutine MK_Distination
	
	subroutine MK_ruletka(n_potok, to_i, to_j, from_i, from_j, area, r, r_per, mu3, bb2)
		real(8), intent(in out) :: mu3
		real(8), intent(in) :: r, r_per
		integer(4), intent(in) :: n_potok, to_i, to_j, from_i, from_j, area
		logical, intent(out) :: bb2
		
		real(8) :: mu, ksi
		
		mu = min( MK_Mu(to_i, to_j, area) * 0.3 * MK_SINKR(to_j) * (r_per/par_Rmax)**2, &
					MK_Mu(from_i, from_j, area) * 0.3 * MK_SINKR(from_j) * (r/par_Rmax)**2) 
		
		if (mu3 >= mu) then
			bb2 = .True.
			return
		else
			call M_K_rand(sensor(1, 2, n_potok), sensor(2, 2, n_potok), sensor(3, 2, n_potok), ksi)
			if (mu3 >= ksi * mu) then
				mu3 = mu
				bb2 = .True.
				return
			else
				mu3 = 0.0
				bb2 = .False.
				return
			end if
		end if
		
		bb2 = .False.
		return
	
	end subroutine MK_ruletka
	
	integer(4) pure function MK_alpha_zones(al)
		real(8), intent(in) :: al
		integer(4) :: i
		
		do i = 1, par_m_zone
			if(al < MK_al_zone(i)) then
				MK_alpha_zones = i
				return
			end if
		end do
		
		MK_alpha_zones = par_m_zone + 1
		return
		
	end function MK_alpha_zones
	
	integer(4) pure function MK_geo_zones(r, k)
		! В какой геометрической зоне сейчас находится атом (зоны считаются с нуля)
		real(8), intent (in) :: r, k
		integer(4) :: i

		do i = 1, par_n_zone
			if (r < k * MK_R_zone(i)) then
				MK_geo_zones = i
				return
			end if
		end do
		
		MK_geo_zones = par_n_zone + 1
		
	end function MK_geo_zones
	
	real(8) pure function MK_F(X, gam, Y)
		real(8), intent (in) :: X, gam, Y
		real(8) :: bb, gg, Y2, Y4, X2, A, B, C, D
		
		bb = sqrt(1.0 + gam)
		gg = 1.0 + gam
		Y2 = Y**2
		Y4 = Y**4
		X2 = X**2
		A = exp(-Y2) / (2.0 * par_sqrtpi * Y * bb**7)
		B = 2 * X * Y * bb**3 * (2.0 + gam * (3.0 + (X2 - 1.0) * Y2 + gam))
		C = par_sqrtpi * gg * (gg * (gg * (4.0 + gam) + Y2 * (2.0 + 3.0 * gam)) + &
			exp(X2 * Y2 / gg) * (2.0 * X2 * (X2 - 1.0) * Y4 * gam - (gg**2) * (4.0 + gam) + &
				Y2 * gg * (-2.0 - 3.0 * gam + X2 * (2.0 + gam))))
		D = par_sqrtpi * gg * exp(X2 * Y2 / gg) * &
			(-2.0 * X2 * (X2 - 1.0) * Y4 * gam + gg**2 * (4.0 + gam) - Y2 * gg * (-2.0 - 3.0 * gam + X2 * (2.0 + gam))) * &
			erf(X * Y / bb)
			
		MK_F =  A * (B + C - D)
	end function MK_F
	
	real(8) pure function MK_FI(Z, X, gam, Y)
		real(8), intent (in) :: Z, X, gam, Y
		real(8) :: bb, gg, X2, Y2, Y4, Z2, A, B, C, D, ee
		
		bb = sqrt(1.0 + gam)
		gg = 1.0 + gam
		X2 = X**2
		Y2 = Y**2
		Y4 = Y**4
		Z2 = Z**2
		A = exp(-Y2 - 2.0 * X * Y * Z - Z2 * gg) / (par_sqrtpi * bb**7)
		ee = exp((X * Y + Z + Z * gam)**2 / gg)
		B = ee * par_sqrtpi * X * Y * (2.0 * X2 * (X2 - 1.0) * Y4 * gam - 2.0 * gg**2 + (X2 - 1.0) * Y2 * gg * (2.0 + 5.0 * gam))
		C = 2.0 * bb * (X2 * (X2 - 1.0) * Y4 * gam - X * (X2 - 1.0) * Y**3 * Z * gam * gg - gg**2 +  &
			(X2 - 1.0) * Y2 * gg * (1 + gam * (2.0 + Z2 * gg)))
		D = ee * par_sqrtpi * X * Y * (2.0 * X2 * (X2 - 1.0) * Y4 * gam - 2.0 * gg**2 + (X2 - 1.0) * Y2 * gg * (2.0 + 5.0 * gam)) * &
			erf((X * Y + Z + Z * gam) / bb)
			
		MK_FI =  A * (B + C + D)
	end function MK_FI
	
	real(8) pure function MK_Hx(gam1, gam2, X, Y, ksi)
		real(8), intent (in) :: gam1, gam2, X, Y, ksi
		MK_Hx = MK_F(X, gam2, Y) - MK_F(X, gam1, Y) - ksi * (MK_F(1.0_8, gam2, Y) - MK_F(1.0_8, gam1, Y))
	end function MK_Hx
	
	real(8) pure function MK_Hwr(gam1, gam2, Z, X, Y, ksi)
		real(8), intent (in) :: gam1, gam2, Z, X, Y, ksi
		MK_Hwr = MK_FI(Z, X, gam2, Y) - MK_FI(Z, X, gam1, Y) - ksi * (MK_FI(0.0_8, X, gam2, Y) - MK_FI(0.0_8, X, gam1, Y))
	end function MK_Hwr
	
	real(8) pure function MK_norm_mho(c)
		real(8), intent (in) :: c
		real(8) :: s, d, cc, nn
		integer(4) :: i
		
		s = 1.0
		d = 0.0
		cc = 1.0
		nn = 1.0
		
		do i = 1, 99
			cc = cc * c
			nn = nn * i
			d = (cc / nn)**2
			s = s + d
			if (d < 0.00001) then
				EXIT
			end if
		end do
		
		MK_norm_mho = s
	end function MK_norm_mho
	
	real(8) function MK_play_mho(potok, c)
		real(8), intent(in) :: c
		integer(4), intent(in) :: potok
		real(8) :: x, ksi
		
		do while(.True.)
			call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi)
			x = 2.0 * ksi * par_pi_8
			call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi)
			if(exp(2.0 * c * cos(x)) > ksi * exp(2.0 * dabs(c)) ) EXIT
		end do

		MK_play_mho = x
	end function MK_play_mho

	
	real(8) pure function MK_Velosity_1(u, cp)
		real(8), intent (in) :: u, cp
		
		if (u < 0.00001) then
			MK_Velosity_1 = 2.0 * cp / par_sqrtpi + 2.0 * u * u / (3.0 * cp * par_sqrtpi) - &
				u * u * u * u / (15.0 * cp * cp * cp * par_sqrtpi)
		else
			MK_Velosity_1 =  exp(-u * u / cp**2) * cp / par_sqrtpi + (u + (cp**2) / (2.0 * u)) * erf(u / cp)
		end if
		
	end function MK_Velosity_1
	
	real(8) pure function MK_Velosity_2(u, cp)
		real(8), intent (in) :: u, cp
		
		if (u < 0.00001) then
			MK_Velosity_2 = (8.0 / 3.0) * cp**4 * par_pi_8 * u + (8.0 / 15.0) * cp**2 * par_pi_8 * u * u * u - &
				(4.0 / 105.0) * par_pi_8 * u**5
		else
			MK_Velosity_2 =  cp**3 * par_pi_8 * (exp(-u * u / cp**2) * cp * u * 2.0 * (cp**2 + 2.0 * u**2) + &
			par_sqrtpi * (4.0 * u**4 + 4.0 * cp * cp * u**2 - cp**4) * erf(u / cp)) / (4.0 * u * u)
		end if
		
	end function MK_Velosity_2
	
	real(8) pure function MK_Velosity_3(u, cp)
		real(8), intent (in) :: u, cp
		
		if (u < 0.00001) then
			MK_Velosity_3 = 8.0 * cp / (3.0 * par_sqrtpi) + 8.0 * u**2 / (9.0 * cp * par_sqrtpi) - &
				44.0 * u**4 / (135.0 * cp * cp * cp * par_sqrtpi)
		else
			MK_Velosity_3 =  exp(-u**2 / cp**2) * cp * (5.0 * cp**2 + 2.0 * u**2) / (par_sqrtpi * (3.0 * cp**2 + 2.0 * u**2)) + &
			(4.0 * u**4 + 12.0 * cp**2 * u**2 + 3.0 * cp**4) * erf(u / cp) / (2.0 * u * (3.0 * cp**2 + 2.0 * u**2))
		end if
		
	end function MK_Velosity_3
	
	real(8) pure function MK_int_1_f1(x)

		real(8), intent (in) :: x
	
		if (x <= 1.0) then
			MK_int_1_f1 =  6.283185155644284 + 0.000024846677279866114 * x + &
				2.0934329078277405 * x * x + 0.008055998193903208 * x * x * x - &
				0.2355169235647438 * x * x * x * x + 0.03820480582423355 * x * x * x * x * x + &
				0.006992274370591744 * x * x * x * x * x * x
		else if (x <= 3.0) then
			MK_int_1_f1 =  6.437524091973454 - 0.6331520099380095 * x + &
				3.1348881317268997 * x * x - 0.8454201478027856 * x * x * x + &
				0.1004702004260311 * x * x * x * x + 0.0009895488638964746 * x * x * x * x * x - &
				0.000920750276197054 * x * x * x * x * x * x
		else if (x <= 5) then
			MK_int_1_f1 =  4.4920780630505135 + 2.5133093267020654 * x + &
				1.1327223176567935 * x * x - 0.24648691152318875 * x * x * x + &
				0.031326738629523766 * x * x * x * x - 0.0021366031960331384 * x * x * x * x * x + &
				0.00005954097505746697 * x * x * x * x * x * x
		else if (x <= 7) then
			MK_int_1_f1 =  1.9138683588136232 + 5.350374732905213 * x - &
				0.16380205801427633 * x * x + 0.06765898334856263 * x * x * x - &
				0.011071118267864083 * x * x * x * x + 0.0008673476933852199 * x * x * x * x * x - &
				0.00002691859374483661 * x * x * x * x * x * x
		else if (x <= 50.0) then
			MK_int_1_f1 =  1.3138472469154294 + 5.336877156136497 * x + &
				0.020286308991329983 * x * x - 0.9780973533544137 * (x / 10.0)**3 + &
				0.26354051936651874 * (x / 10.0)**4 - 0.03711733070841712 * (x / 10.0)**5 + &
				0.002120935433043921 * (x / 10.0)**6
		else
			MK_int_1_f1 =  1.3138472469154294 + 5.336877156136497 * x + &
				0.020286308991329983 * x * x - 0.9780973533544137 * (x / 10.0)**3 + &
				0.26354051936651874 * (x / 10.0)**4 - 0.03711733070841712 * (x / 10.0)**5 + &
				0.002120935433043921 * (x / 10.0)**6
		end if
				 
		return
	end function MK_int_1_f1
	
	real(8) pure function MK_int_1_f2(x)
		real(8), intent (in) :: x
		
		if (x <= 1.0) then
			MK_int_1_f2 =  1.328216167939543 - 0.000004545681954848391 * x + &
				2.537368073155103 * x * x - 0.0020584991728545624 * x * x * x - &
				0.03742568018912792 * x * x * x * x - 0.010312136385277346 * x * x * x * x * x + &
				0.002767736179209713 * x * x * x * x * x * x
		else if (x <= 3.0) then
			MK_int_1_f2 =  1.2959616295629246 + 0.1533684067037866 * x + &
				2.2354849981206106 * x * x + 0.3113395567715921 * x * x * x - &
				0.21656309882941488 * x * x * x * x + 0.041957500887605075 * x * x * x * x * x - &
				0.0029978773724628604 * x * x * x * x * x * x
		else if (x <= 5.0) then
			MK_int_1_f2 =  1.903643456971281 - 1.4801836911099535 * x + 3.973958664572268 * x * x - &
				0.6482729779428982 * x * x * x + 0.07665007314658864 * x * x * x * x - &
				0.005369758193338703 * x * x * x * x * x + 0.00016605531871992049 * x * x * x * x * x * x
		else if (x <= 7.0) then
			MK_int_1_f2 =  -4.484415105552316 + 5.3747429756690766 * x + &
				0.8892806582308143 * x * x + 0.09767316152573671 * x * x * x - &
				0.025704749778475783 * x * x * x * x + 0.0021937998296249206 * x * x * x * x * x - &
				0.00006928845984076111 * x * x * x * x * x * x
		end if
				 
		return
	end function MK_int_1_f2
	
	real(8) pure function MK_int_1_f3(x)
		real(8), intent (in) :: x
 
		if (x <= 1.0) then
			MK_int_1_f3 = 1.2938345594193854 - 0.000031719847351174835 * x + &
				1.3183710041280094 * x * x - 0.014150512069488197 * x * x * x + &
				0.4226114681928129 * x * x * x * x - 0.06985750969880078 * x * x * x * x * x - &
				0.015347864048406958 * x * x * x * x * x * x
		else if (x <= 3.0) then
			MK_int_1_f3 = 0.9667460440956788 + 1.336271810704016 * x - &
				0.8687355257991665 * x * x + 1.7676868273627229 * x * x * x - &
				0.2731222764016417 * x * x * x * x + 0.004801770033831665 * x * x * x * x * x + &
				0.001780776080720323 * x * x * x * x * x * x
		else if (x <= 5.0) then
			MK_int_1_f3 = 4.760566734123174 - 5.048204299463048 * x + 3.332342585228025 * x * x + &
				0.47584339615235993 * x * x * x - 0.12072786272726124 * x * x * x * x + &
				0.011870955604980658 * x * x * x * x * x - 0.0004580199652402304 * x * x * x * x * x * x
		else if (x <= 7.0) then
			MK_int_1_f3 = 9.370493362469261 - 10.848615619383615 * x + 6.423326878282571 * x * x - &
				0.4148977656870439 * x * x * x + 0.025300923044176957 * x * x * x * x - &
				0.0010108688120876522 * x * x * x * x * x + 0.00001864423130429156 * x * x * x * x * x * x
		end if
			 
		return
	end function MK_int_1_f3
	
	real(8) pure function MK_int_1(x, cp)
		real(8), intent (in) :: x, cp
		real(8) :: b
		
		b = 1.0 - par_a_2 * log(cp)
		MK_int_1 = (cp / (par_sqrtpi**3)) * (b * b * MK_int_1_f1(x / cp) - &
			2.0 * par_a_2 * b * MK_int_1_f2(x / cp) + par_a_2**2 * MK_int_1_f3(x / cp))
	end function MK_int_1
	
	real(8) pure function MK_int_2_f1(x)
		real(8), intent (in) :: x
	
		if (x <= 1.0) then
				MK_int_2_f1 =   8.377571213788123 * x + 0.00047608508679086725 * x * x + &
					1.6710478320575737 * x * x * x + 0.016857530811432916 * x * x * x * x - &
					0.15132474125724812 * x * x * x * x * x + 0.030723378194358945 * x * x * x * x * x * x
				return
		else if (x <= 3.0) then
				MK_int_2_f1 =   -0.11788367995598747 + 8.937936705157014 * x - &
						1.119886471634001 * x * x + 2.8831031948885917 * x * x * x - &
						0.735250146386514 * x * x * x * x + 0.10356311378423572 * x * x * x * x * x - &
						0.006231417172309398 * x * x * x * x * x * x
				return
		else if (x <= 5.0) then
				MK_int_2_f1 =   2.9044497739429858 + 2.757712415967557 * x + 4.239161941189675 * x * x + &
						0.36198838294786784 * x * x * x - 0.05737777787138304 * x * x * x * x + &
						0.004956250079677106 * x * x * x * x * x - 0.0001809238236975877 * x * x * x * x * x * x
				return
		else if (x <= 7.0) then
				MK_int_2_f1 =  41.6323689028892 - 38.118317864344135 * x + 22.211189528076645 * x * x -&
					3.8547348524931246 * x * x * x + 0.5000517174807501 * x * x * x * x -&
					0.03446294709493891 * x * x * x * x * x + 0.0009860204070962582 * x * x * x * x * x * x
				return
		end if
	
		MK_int_2_f1 =  0.0
	end  function MK_int_2_f1
	
	real(8) pure function MK_int_2_f2(x)

		real(8), intent (in) :: x
	
		if (x <= 1.0) then
			MK_int_2_f2 =  3.8653461103376667 * x + 0.0001975300512691014 * x * x + &
				2.4468141895384012 * x * x * x + 0.005987984681429616 * x * x * x * x - &
				0.06453987836713967 * x * x * x * x * x + 0.0066920981111229004 * x * x * x * x * x * x
			return
		else if (x <= 3.0) then
			MK_int_2_f2 =  -0.10983889480661446 + 4.321087890898017 * x - &
				0.7707850845797619 * x * x + 3.1237901158486583 * x * x * x - &
				0.31485222316123385 * x * x * x * x + 0.010270249760261791 * x * x * x * x * x + &
				0.0008259803934338584 * x * x * x * x * x * x
			return
		else if (x <= 5.0) then
			MK_int_2_f2 =  1.8468011847509729 + 1.1986396743254275 * x + &
				1.1421489029509448 * x * x + 2.606316149569781 * x * x * x - &
				0.2788783468089509 * x * x * x * x + 0.019815317035281846 * x * x * x * x * x - &
				0.0006379970557448899 * x * x * x * x * x * x
			return
		else if (x <= 7.0) then
			MK_int_2_f2 =  9.480707804348185 - 8.022988228952784 * x + 5.823555900242488 * x * x + &
				1.3277220473440972 * x * x * x - 0.08074921612981413 * x * x * x * x + &
				0.0033058587723954185 * x * x * x * x * x - 0.00006041810279926061 * x * x * x * x * x * x
			return
		end if
			 
		MK_int_2_f2 = 0.0
	end function MK_int_2_f2
		
	real(8) pure function MK_int_2_f3(x)
        real(8), intent (in) :: x
		if (x <= 1.0) then
			MK_int_2_f3 =  2.6106039258326 * x - 0.0008357997793049243 * x * x + &
				2.0764571907368174 * x * x * x - 0.03182275644841273 * x * x * x * x + &
				0.26310521962808975 * x * x * x * x * x - 0.06034325471643871 * x * x * x * x * x * x
			return
		else if (x <= 3.0) then
			MK_int_2_f3 =   0.20784760901369737 + 1.5920325291316857 * x + &
				2.0985329535259014 * x * x - 0.26286255221171206 * x * x * x + &
				1.4610329096120886 * x * x * x * x - 0.25626862029131897 * x * x * x * x * x + &
				0.01684969647300594 * x * x * x * x * x * x
			return
		else if (x <= 5.0) then
			MK_int_2_f3 =   -6.284115352064703 + 15.665162343948523 * x &
				- 10.766105772158252 * x**2 + 6.0821074614870465 * x**3 - &
				0.3181501403196319 * x**4 + 0.01232319194701587 * x**5 &
				- 0.00017890661550876597 * x**6
			return
		else if (x <= 7.0) then
			MK_int_2_f3 =   -4.355962170454177 + 13.835332665069274 * x &
				- 10.12766071646978 * x**2 + 5.999392227482686 * x**3 - &
				0.32171318647989955 * x**4 + 0.014181987261856027 * x**5 &
				- 0.00030579035551447497 * x**6
			return
		end if
				
		MK_int_2_f3 = 0.0
		
	end function MK_int_2_f3
	
	real(8) pure function MK_int_2(x, cp)
		real(8), intent (in) :: x, cp
		real(8) :: b
		
		b = 1.0 - par_a_2 * log(cp)
		MK_int_2 = -(cp * cp / (par_sqrtpi * par_sqrtpi * par_sqrtpi)) * (b * b * MK_int_2_f1(x / cp) - &
			2.0 * par_a_2 * b * MK_int_2_f2(x / cp) + (par_a_2**2) * MK_int_2_f3(x / cp))
	end function MK_int_2
	
	real(8) pure function MK_int_3(x, cp)
		real(8), intent (in) :: x, cp
		real(8) :: b
		
		b = 1.0 - par_a_2 * log(cp)
		MK_int_3 = (cp**3 / (par_sqrtpi**3)) * (b * b * MK_int_3_f1(x / cp) - 2.0 * par_a_2 * b * MK_int_3_f2(x / cp) + (par_a_2**2) * MK_int_3_f3(x / cp))
	end function MK_int_3

	real(8) pure function MK_int_3_f1(x)
		real(8), intent (in) :: x
		
		if (x <= 1.0) then
			MK_int_3_f1 = 12.566370586001975 - 0.00001944816384202852 * x + 12.567558607381049 * x**2 - &
				0.010507444068349692 * x**3 + 1.2911398125420694 * x**4 - 0.05048482363937502 * x**5 - &
				0.029947937607835744 * x**6
		else if (x <= 3.0) then
			MK_int_3_f1 = 12.451555448799724 + 0.40252674353016715 * x + 12.081033298182223 * x**2 + 0.12939193776331415 * x**3 + &
				1.478876561367302 * x**4 - 0.2237491583356496 * x**5 + 0.014474521138620033 * x**6
		else if (x <= 5.0) then
			MK_int_3_f1 = 8.281962852844913 + 9.783032429527339 * x + 3.165848344887614 * x**2 + 4.711462666119614 * x**3 + &
				0.13739453130827392 * x**4 - 0.012096015315889195 * x**5 + 0.0004514225555943018 * x**6
		else if (x <= 7.0) then
			MK_int_3_f1 = 17.29966669035025 + 2.1457227895928668 * x + 5.572082327920818 * x**2 + 4.4083748449004645 * x**3 + &
				0.13640200155890422 * x**4 - 0.00854302147917508 * x**5 + 0.00022205921430504255 * x**6
		end if
		
	end function MK_int_3_f1

	real(8) pure function MK_int_3_f2(x)
		real(8), intent (in) :: x
		
		if (x <= 1.0) then
			MK_int_3_f2 = 5.798024979296493 - 0.00001772671478406096 * x + 9.98769025405073 * x**2 - 0.0073593516014156535 * x**3 + &
				2.27820901023418 * x**4 - 0.03135086958655956 * x**5 - 0.030403716978821237 * x**6
		else if (x <= 3.0) then
			MK_int_3_f2 = 5.864728705834779 - 0.34799875480550213 * x + 10.760249318358127 * x**2 - 0.9493738978205943 * x**3 + &
				2.948810798239708 * x**4 - 0.29690823082625284 * x**5 + 0.015284639719532296 * x**6
		else if (x <= 5.0) then
			MK_int_3_f2 = -3.21810405152587 + 17.08054504204813 * x - 3.43427364926184 * x**2 + 5.342946243936558 * x**3 + &
				1.3459970589832742 * x**4 - 0.07448103623442687 * x**5 + 0.0021616888853670958 * x**6
		else if (x <= 7.0) then
			MK_int_3_f2 = -59.42860988375287 + 79.24709914283142 * x - 32.304295129821256 * x**2 + 12.558114389999929 * x**3 + &
				0.32138208180756495 * x**4 + 0.003976327853989966 * x**5 - 0.0003703205838879336 * x**6
		end if
		
	end function MK_int_3_f2

	real(8) pure function MK_int_3_f3(x)
		real(8), intent(in) :: x
		
		if (x <= 1.0) then
			MK_int_3_f3 = 3.915885322797866 + 0.000044284982651632276 * x + 7.778946280540595 * x**2 + 0.01990833982643192 * x**3 + &
				2.710531013102172 * x**4 + 0.09480498076917274 * x**5 + 0.026454483470905545 * x**6
		else if (x <= 3.0) then
			MK_int_3_f3 = 4.371710139648201 - 1.787430730965795 * x + 10.530386754306065 * x**2 - 1.9945949141813966 * x**3 + &
				3.310710912254515 * x**4 + 0.1356460090430464 * x**5 - 0.019853464614841342 * x**6
		else if (x <= 5.0) then
			MK_int_3_f3 = 4.6940208896644435 - 6.423081982183987 * x + 18.137713177954524 * x**2 - 7.298291873534797 * x**3 + &
				5.202095367169497 * x**4 - 0.20620743501245353 * x**5 + 0.005094092519200813 * x**6
		else if (x <= 7.0) then
			MK_int_3_f3 = 194.2826492092263 - 195.77327124311523 * x + 95.83253215419927 * x**2 - 23.99281397751645 * x**3 + &
				7.170549820159753 * x**4 - 0.32568903981020303 * x**5 + 0.007955090180048264 * x**6
		end if
		
	end function MK_int_3_f3

	real(8) pure function MK_for_Wr_1(Z, gam, ur)
		!! Для розыгрыша Wr в перезарядке по-частям (первая часть)
		real(8), intent (in) :: Z, gam, ur
		
		MK_for_Wr_1 = -exp(-gam * ur**2 / (gam + 1.0)) * (1.0 + erf((-ur + gam * Z + Z) / sqrt(gam + 1.0))) / (2.0 * sqrt(gam + 1.0))
	end function MK_for_Wr_1

	real(8) pure function MK_for_Wr_2(Z, gam, ur, ut)
		!! Для розыгрыша Wr в перезарядке по-частям 
		real(8), intent (in) :: Z, gam, ur, ut
		real(8) :: gam1, ur2
		
		gam1 = gam + 1.0  ! gam + 1
		ur2 = ur**2
		
		MK_for_Wr_2 =  1.0 / (4.0 * gam1 * gam1 * gam1 * par_sqrtpi) * ut**2 * exp((-1.0 - gam / gam1) * ur2 - gam1 * Z * Z) * &
			(2.0 * exp(ur * (gam * ur / gam1 + 2.0 * Z)) * gam * gam1 * (ur + Z + gam * Z) + &
				exp(ur2 + gam1 * Z * Z) * sqrt(gam1) * par_sqrtpi * (2.0 + gam * (5.0 + 3.0 * gam + 2.0 * ur2)) * &
				(-1.0 - erf((-ur + Z + gam * Z) / sqrt(gam1))))
	end function MK_for_Wr_2

	real(8) pure function MK_H_Wr_1(gam1, gam2,	V, ur, p, ksi)
		real(8), intent (in) :: gam1, gam2, V, ur, p, ksi
		MK_H_Wr_1 = MK_for_Wr_1(V, gam2, ur) - MK_for_Wr_1(V, gam1, ur) - ksi * p
	end function MK_H_Wr_1

	real(8) pure function MK_H_Wr_2(gam1, gam2, V, ur, ut, p, ksi)
		! Для розыгрыша Wr в перезарядке по-частям 
		real(8), intent (in) :: gam1, gam2, V, ur, ut, p, ksi
		MK_H_Wr_2 =  MK_for_Wr_2(V, gam2, ur, ut) - MK_for_Wr_2(V, gam1, ur, ut) - ksi * p
	end function MK_H_Wr_2

	!real(8) pure function MK_sigma(x)
	!	real(8), intent (in) :: x
	!	MK_sigma = (1.0 - par_a_2 * log(x))**2
	!end function MK_sigma
	!
	!real(8) pure function MK_sigma2(x, y)
	!	real(8), intent (in) :: x, y
	!	MK_sigma2 = (1.0 - par_a_2 * log(x * y))**2
	!end function MK_sigma2

	real(8) pure function MK_f2(V, gam, ur, ut)
		real(8), intent (in) :: V, gam, ur, ut
		real(8) :: b, b4, b2, A
		
		b = sqrt(1.0 + gam)
		b4 = b**4
		b2 = b**2

		A = (exp(2.0 * V * ur - (ur**2) - (V**2) * b2) / (2.0 * par_sqrtpi * b4))
		
		if (A <= 0.0000000000001) then
			MK_f2 = 0.0
			return
		end if
			
		MK_f2 = A * (-gam * (ut**2) * (ur + gam * V + V) + &
			(par_sqrtpi / (2.0 * b)) * exp( (V + gam * V - ur)**2 / b2) * (2.0 * b4 + (ut**2) * (b2 * (3.0 * gam + 2.0) + 2.0 * gam * (ur**2))) * &
			(1.0 + erf((-ur + gam * V + V) / b)))
	end function MK_f2

	subroutine Get_sensor(mpi_rank)
		! Считываем датчики случайных чисел из файла
		! Variables
		integer, intent(in) :: mpi_rank
		logical :: exists
		integer(4) :: i, a, b, c, j
		
		
		inquire(file="rnd_my.txt", exist=exists)
		
		if (exists == .False.) then
			pause "net faila!!!  345434wertew21313edftr3e"
			STOP "net faila!!!"
		end if
		
		if (par_n_claster * par_n_potok * 2 > 1021) then
			print*, "NE XVATAET DATCHIKOV 31 miuhi8789pok9"
			pause
		end if
		
		open(1, file = "rnd_my.txt", status = 'old')
		
		do j = 1, mpi_rank + 1
			do i = 1, par_n_potok
				read(1,*) a, b, c
				sensor(:, 1, i) = (/ a, b, c /)
				read(1,*) a, b, c
				sensor(:, 2, i) = (/ a, b, c /)
			end do
		end do
		
		close(1)
	end subroutine Get_sensor
	
	subroutine Get_sensor_sdvig(sdvig)
		! Считываем датчики случайных чисел из файла
		! Variables
		integer, intent(in) :: sdvig
		logical :: exists
		integer(4) :: i, a, b, c, j
		
		
		inquire(file="rnd_my.txt", exist=exists)
		
		if (exists == .False.) then
			pause "net faila!!!  345434wertew21313edftr3e"
			STOP "net faila!!!"
		end if
		
		if (par_n_claster * par_n_potok * 2 > 1021) then
			print*, "NE XVATAET DATCHIKOV 31 miuhi8789pok9"
			pause
		end if
		
		open(1, file = "rnd_my.txt", status = 'old')
		
		do i = 1, sdvig
			read(1,*) a, b, c
			read(1,*) a, b, c
		end do
		
		do i = 1, par_n_potok
			read(1,*) a, b, c
			sensor(:, 1, i) = (/ a, b, c /)
			read(1,*) a, b, c
			sensor(:, 2, i) = (/ a, b, c /)
		end do
		
		close(1)
	end subroutine Get_sensor_sdvig

	subroutine PUI_proverka(x, y, z)
		implicit none
		real(8), intent(in) :: x, y, z
		integer(4) :: n1, cell2
		real(8) :: uH, vH, wH, u, uz, nu_ex, r, u1, u2, u3, skalar
		real(8) :: cp, vx, vy, vz, ro, PAR(9), MAS_PUI(2), p, ro_pui  ! Параметры плазмы в ячейке

		print*, "PUI_proverka = ", x, y, z
		uH = -2.3
		vH = 0.0
		wH = 0.0
		r = sqrt(x**2 + y**2 + z**2)
		n1 = 3
		call Int2_Get_tetraedron_inner(x, y, z, n1)
		! n1 - номер тетраэдра
		call Int2_Get_par_fast2(x, y, z, n1, PAR, MAS_PUI = MAS_PUI)
		cell2 = int2_all_tetraendron_point(1, n1)

		p = PAR(5) - MAS_PUI(2) * MAS_PUI(1)
		if(p < 0.0) p = 4286.72/((r/par_1ae)**(2.0 * ggg))
		ro_pui = MAS_PUI(1)
		ro = PAR(1) - ro_pui
		if(ro < 0.000001) then
			ro = 0.000001
			p = 0.000001
		end if
		cp = sqrt(p/ro)
		vx = PAR(2)
		vy = PAR(3)
		vz = PAR(4)

		! Найдём время до перезарядки и веса частиц  ****************************************************************************************
		u = sqrt(kvv(uH - vx, vH - vy, wH - vz))
		u1 =  vx - uH
		u2 =  vy - vH
		u3 =  vz - wH
		skalar = uH * u1 + vH * u2 + wH * u3

		print*, "ro, p, ro_pui, cp", ro, p, MAS_PUI(1), cp
	
		if (u / cp > 7.0) then
			uz = MK_Velosity_1(u, cp);
			nu_ex = (ro * uz * MK_sigma(uz)) / par_Kn
		else
			nu_ex = (ro * MK_int_1(u, cp)) / par_Kn        ! Пробуем вычислять интеграллы численно
		end if
		print*, "pui_num, u = ", f_pui_num2(cell2), u
		print*, "Chastota ex_p = ", nu_ex
		print*, "Chastota ex_pui = ", PUI_get_nu_integr(f_pui_num2(cell2), u)/ par_Kn

		print*, "END _____________________________________________"

	end subroutine PUI_proverka
	
end module Monte_Karlo 
	
	
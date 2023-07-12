
	
module Monte_Karlo  
	! Модуль Монте-Карло, использует интерполяционную сетку тетраэдров
	USE GEO_PARAM
	USE STORAGE
	USE Interpolate2
	
	implicit none
	
	integer(4), parameter :: par_stek = 1000000  ! Глубина стека (заранее выделяется память под него)
	integer(4), parameter :: par_n_potok = 1  ! Число потоков (у каждого потока свой стек)
	integer(4), parameter :: par_n_zone = 7  ! 
	real(8), parameter :: par_Rmax = 200.0  !  Радиус сферы, с которой запускаем частицы
	
	real(8) :: MK_R_zone(par_n_zone)   ! Радиусы зон
	
	
	real(8) :: MK_gam_zone(par_n_zone)   ! Параметр гамма для зон
	real(8) :: MK_A0_, MK_A1_   ! Параметры для начального запуска
	
	
	
	integer(4), parameter :: par_num_1 = 1   ! Число исходных частиц первого типа (с полусферы)
	
	real(8), allocatable :: M_K_particle(:, :, :)   ! Частицы (7, par_stek, число потоков)
	! (три координаты, три скорости, вес)
	integer(4), allocatable :: M_K_particle_2(:, :, :)  ! Частицы (1, par_stek, число потоков)
	! (в какой ячейке частица)
	
	integer(4), allocatable :: sensor(:, :, :)  !(3, 2, : par_n_potok число потоков)  ! датчики случайных чисел 
	! Каждому потоку по два датчика
	
	integer(4), allocatable :: stek(:)   ! (: число потоков) Переменная чтения и записи в стек
	! Где стоит переменная, там что-то лежит, чтобы записать, нужно увеличить значение на 1
	
	contains
	
	
	subroutine M_K_start()
	! Variables
	integer(4) :: potok, num, i, cell
	real(8) :: mu_(par_n_zone + 1), Wt_(par_n_zone + 1), Wp_(par_n_zone + 1), Wr_(par_n_zone + 1), X_(par_n_zone + 1)
	logical :: bb
	real(8) :: sin_, x, phi, y, z, ksi, Vx, Vy, Vz
	
	call M_K_Set()    ! Создали массивы
	call Get_sensor() ! Считали датчики случайных чисел
	
	
	! Запускаем каждый поток в параллельном цикле
	do potok = 1, par_n_potok    ! Локальные:  num, mu_, Wt_, Wp_, Wr_, X_, bb, i, ksi, Vx, Vy, Vz, cell
		
		! Запускаем частицы первого типа
		do num = 1, par_num_1
			call MK_Init_Parametrs(potok, mu_, Wt_, Wp_, Wr_, X_, bb)
			
			do i = 1, par_n_zone + 1
				sin_ = sqrt(1.0 - (X_(i)**2))
				x = (par_Rmax) * X_(i)
				call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi)
				phi = 2.0 * par_pi_8 * ksi
				y = (par_Rmax) * sin_ * cos(phi)
				z = (par_Rmax) * sin_ * sin(phi)
				
				cell = 3
				call Int2_Get_tetraedron( x, y, z, cell)
				call dekard_skorost(x, y, z, Wr_(i), Wp_(i), Wt_(i), Vx, Vy, Vz)
				
				if(i /= par_n_zone + 1 .or. bb == .True.) then
					! Добавляем частицу в стек
					stek(potok) = stek(potok) + 1
					M_K_particle(:, stek(potok), potok) = (/ x, y, z, Vx, Vy, Vz, mu_(i) /)
					M_K_particle_2(1, stek(potok), potok) = cell
				end if
			end do
			call M_K_Fly(potok)
			
		end do
		
	end do
	
	
	!i = 3
	!
	!
	!
	!M_K_particle(:, 1, 1) = (/ 35.0, 20.0, 3.0, -3.0, 0.01, 0.333, 1.0 /)
	!stek(1) = 1
	!call Int2_Get_tetraedron( M_K_particle(1, 1, 1), M_K_particle(2, 1, 1), M_K_particle(3, 1, 1), i)
	!M_K_particle_2(1, 1, 1) = i
	!
	!call M_K_Fly(1)
	!
	
	end subroutine M_K_start
	

	
	subroutine M_K_Fly(n_potok)
	! Функция запускает все частицы в стеке потока + все дочерние частицы
	
	integer(4), intent(in) :: n_potok  ! Номер потока 
	
	real(8) :: particle(7)
	integer(4):: particle_2(1)
	
	integer(4) :: num  ! Номер частицы, верхняя в стеке
	integer(4) :: cell ! Номер ячейки, в которой находится частица
	integer(4) :: next ! Номер ячейки, в которую попадёт частица в следующий раз
	
	real(8) :: time ! Оценочное время до вылета частицы из ячейки
	real(8) :: time2 ! Время до перезарядки
	
	
	do while (stek(n_potok) >= 1)
		num = stek(n_potok)
		stek(n_potok) = stek(n_potok) - 1
		! Берём все параметры частицы
		particle = M_K_particle(:, num, n_potok)
		particle_2 = M_K_particle_2(:, num, n_potok)
		
		print*, "stek(n_potok) = ", stek(n_potok)
		print*, particle
		print*, particle_2
		pause
		
		
		loop1: do  ! пока частица не вылетит из области
			cell = particle_2(1)
		
			call Int2_Time_fly(particle(1:3), particle(4:6), time, cell, next)
			
			time = max(0.00000001, time * 1.001) ! Увеличим время, чтобы частица точно вышла из ячейки
		
			! Найдём время до перезарядки и веса частиц
		
			! Найдём скорость перезаряженной частицы
		
			! Запишем новую частицу в стек
		
			!stek(n_potok) = stek(n_potok) + 1
			!M_K_particle(:, stek(n_potok), n_potok) = particle ! поменять
		
			! Находим следующую ячайку
			
			particle(1:3) = particle(1:3) + time * particle(4:6)
			if(next == 0) EXIT loop1  ! частица долетела до края области
			
11			continue 
			
			!if(particle(1) < -9.6) then
			!	continue
			!end if
			
			
			
			call Int2_Get_tetraedron(particle(1), particle(2), particle(3), next)  ! Находим точно следующий тетраэдр
			if (next == -1) then ! ЗАТЫК
				next = 3
				particle(1:3) = particle(1:3) + (time/100.0) * particle(4:6)
				print*, "Zatik"
				pause
				GO TO 11
			end if
			
			particle_2(1) = next
			if(next == 0) EXIT loop1
		
		end do loop1
		
		
	end do
	
	
	
	! Variables
	
	! Body of M_K_Start
	
	
	end subroutine M_K_Fly
	
	
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
	
		integer(4) :: i
		real(8) :: Yr
		
		allocate(M_K_particle(7, par_stek, par_n_potok))
		allocate(M_K_particle_2(1, par_stek, par_n_potok))
		allocate(sensor(3, 2, par_n_potok))
		allocate(stek(par_n_potok))
	
		M_K_particle = 0.0
		M_K_particle_2 = 0
		stek = 0
		sensor = 1
		
		MK_R_zone(1) = 1.0
		MK_R_zone(2) = 2.15
		MK_R_zone(3) = 4.6
		MK_R_zone(4) = 10.0
		MK_R_zone(5) = 21.0
		MK_R_zone(6) = 46.0
		MK_R_zone(7) = 99.0
		
		! Задаём начальные параметры
		do i = 1, par_n_zone
			MK_gam_zone(i) = 1.0 / ((par_Rmax / MK_R_zone(i))**2 - 1.0)
			if (MK_gam_zone(i) < 0.0) STOP "ERROR gamma 56789oihgfr6uijt6789"
		end do
		
		Yr = dabs(par_Velosity_inf)
		MK_A0_ = (Yr + 1.0 / (2.0 * Yr)) * erf(Yr) + exp(-(Yr**2)) / par_sqrtpi + Yr
		MK_A1_ = 1.0 + (1.0 + 1.0 / (2.0 * (Yr)**2 )) * erf(Yr) + exp(-(Yr**2) ) / (par_sqrtpi * Yr)
		
	
	end subroutine M_K_Set
	
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
	
	
	p1 = erf(Y) / (MK_A1_ * (Y**2))
	X0 = 1.0
	Y = dabs(par_Velosity_inf)

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
			call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi)
			x = 2.0 * ksi * par_pi_8
			call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi)
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
	
		print*, "Error  int_f1 > 7  =  ", x
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
			 
		print*, "Error  int2_f2 > 7  =  ", x
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
				
		print*, "Error  int2_f3 > 7  =  ", x
		MK_int_2_f3 = 0.0
		
	end function MK_int_2_f3
	
	real(8) pure function MK_int_2(x, cp)
		real(8), intent (in) :: x, cp
		real(8) :: b
		
		b = 1.0 - par_a_2 * log(cp)
		MK_int_2 = -(cp * cp / (par_sqrtpi * par_sqrtpi * par_sqrtpi)) * (b * b * MK_int_2_f1(x / cp) - &
			2.0 * par_a_2 * b * MK_int_2_f2(x / cp) + (par_a_2**2) * MK_int_2_f3(x / cp))
    end function MK_int_2
	
	
	subroutine Get_sensor()
	! Считываем датчики случайных чисел из файла
	! Variables
	logical :: exists
	integer(4) :: i, a, b, c
    
	
    inquire(file="rnd_my.txt", exist=exists)
    
    if (exists == .False.) then
		pause "net faila!!!  345434wertew21313edftr3e"
        STOP "net faila!!!"
    end if
	
	if (par_n_potok * 2 > 1021) then
		print*, "NE XVATAET DATCHIKOV 31 miuhi8789pok9"
		STOP
	end if
	
	open(1, file = "rnd_my.txt", status = 'old')
    
	do i = 1, par_n_potok
		read(1,*) a, b, c
		sensor(:, 1, i) = (/ a, b, c /)
		read(1,*) a, b, c
		sensor(:, 2, i) = (/ a, b, c /)
	end do
	
	close(1)
	
	end subroutine Get_sensor
	
end module Monte_Karlo 
	
	
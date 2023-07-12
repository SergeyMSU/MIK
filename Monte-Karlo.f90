
	
module Monte_Karlo  
	! Модуль Монте-Карло, использует интерполяционную сетку тетраэдров
	USE GEO_PARAM
	USE STORAGE
	USE Interpolate2
	
	implicit none
	
	integer(4), parameter :: par_stek = 1000000  ! Глубина стека (заранее выделяется память под него)
	integer(4), parameter :: par_n_potok = 1  ! Число потоков (у каждого потока свой стек)
	
	real(8), allocatable :: M_K_particle(:, :, :)   ! Частицы (7, par_stek, число потоков)
	! (три координаты, три скорости, вес)
	integer(4), allocatable :: M_K_particle_2(:, :, :)  ! Частицы (1, par_stek, число потоков)
	! (в какой ячейке частица)
	
	integer(4), allocatable :: sensor(:, :)  !(3, : число потоков)  ! датчики случайных чисел
	
	integer(4), allocatable :: stek(:)   ! (число потоков) Переменная чтения и записи в стек
	! Где стоит переменная, там что-то лежит, чтобы записать, нужно увеличить значение на 1
	
	contains
	
	subroutine M_K_start()
	! Variables
	integer(4) :: i
	
	i = 3
	
	call M_K_Set()
	
	M_K_particle(:, 1, 1) = (/ 35.0, 20.0, 3.0, -3.0, 1.01, 0.333, 1.0 /)
	stek(1) = 1
	call Int2_Get_tetraedron( M_K_particle(1, 1, 1), M_K_particle(2, 1, 1), M_K_particle(3, 1, 1), i)
	M_K_particle_2(1, 1, 1) = i
	
	call M_K_Fly(1)
	
	! Body of M_K_start
	
	end subroutine M_K_start
	
	subroutine M_K_Fly(n_potok)
	! Функция запускает все частицы в стеке потока + все новые частицы
	
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
		
		loop1: do  ! пока частица не вылетит из области
			cell = particle_2(1)
			print*, particle(1:3)
		
			call Int2_Time_fly(particle(1:3), particle(4:6), time, cell, next)
			!print*, "A = ", next
			
			time = max(0.00000001, time * 1.001) ! Увеличим время, чтобы частица точно вышла из ячейки
		
			! Найдём время до перезарядки и веса частиц
		
			! Найдём скорость перезаряженной частицы
		
			! Запишем новую частицу в стек
		
			!stek(n_potok) = stek(n_potok) + 1
			!M_K_particle(:, stek(n_potok), n_potok) = particle ! поменять
		
			! Находим следующую ячайку
			
			particle(1:3) = particle(1:3) + time * particle(4:6)
			print*, particle(1:3)
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
			
			!print*, "B = ", next
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
	
		allocate(M_K_particle(7, par_stek, par_n_potok))
		allocate(M_K_particle_2(1, par_stek, par_n_potok))
		allocate(sensor(3, par_n_potok))
		allocate(stek(par_n_potok))
	
		M_K_particle = 0.0
		M_K_particle_2 = 0
		stek = 0
		sensor = 1
	
	end subroutine M_K_Set
	
end module Monte_Karlo 
	
	
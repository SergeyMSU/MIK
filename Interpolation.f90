! Алгоритм интерполяции основан на исходной сетке
! Сначала нужно найти ячейку, в которой находится нижная интерполяционная точка
! Далее определить по каким соседям нужно построить интерполяционную фигуру (тетраэдр)
! Далее нужно интерполировать значения нужных переменных внутри тетраэдра и вывести ответ
	
! ПРОБЛЕМА
! В интерполяции такого рода есть "дыры". А именно точка может лежать возле пересечения 4-х ячеек и в силу
! погрешности и неточности определения грани по четырём точкам алгоритм будет прыгать покругу 
! (по этим четырём ячейкам) и не найдёт, где именно она лежит
	
	module Interpolate2  ! МОДУЛЬ НЕ ДОДЕЛАН
	! Модуль интерполяции по тетраэдрам
    ! Строится двойственная сетка к основной (с дополнительными узлами на разрывах)
	! Строиться сетка тетраэдров на двойственной сетке
	USE GEO_PARAM
	USE STORAGE
	
	implicit none
	
	
	logical :: int2_work             ! Определена ли интерполяционная сетка (выделена ли память и т.д.)
	
	real(8), allocatable :: int2_coord(:, :)    ! (3, :) набор координат двойственной сетки (центры ячеек сетки 1 + дополнительные)
	
	
	
    integer(4), allocatable :: int2_all_Cell(:, :)   ! Весь набор ячеек (8, :) - первая координата массива - это набор узлов ячейки
	integer(4), allocatable :: int2_Cell_A(:, :, :)
	integer(4), allocatable :: int2_Cell_B(:, :, :)
	integer(4), allocatable :: int2_Cell_C(:, :, :)
	
	real(8), allocatable :: int2_Cell_center(:, :)   ! Центр ячеек (3, :) 
	integer(4), allocatable :: int2_all_neighbours(:, :)  ! (6, :)  по 6 соседей на каждую ячейки
	
	integer(4), allocatable :: int2_Point_A(:, :, :)   ! Набор A-точек размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
	! Это точки - центры ячеек основной сетки
	integer(4), allocatable :: int2_Point_B(:, :, :)   ! Набор B-точек размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
	integer(4), allocatable :: int2_Point_C(:, :, :)   ! Набор B-точек размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
	
	contains
	
	subroutine Int2_Initial()
	
	integer(4) :: i, j, k, N1, N2, N3, ijk, N, kk, m, ii, jj
	! Заполняем интерполяционную сетку из основной
	int2_Point_A(1, :, :) = -1  ! Центральная точка
	int2_Point_B(1, :, :) = -1  ! Центральная точка
	
	
	N1 = size(gl_Cell_A(:, 1, 1))
	N2 = size(gl_Cell_A(1, :, 1))
	N3 = size(gl_Cell_A(1, 1, :))
	
	do k = 1, N3
		do j = 1, N2
			do i = 1, par_n_TS - 1
				int2_Point_A(i + 1, j + 1, k) = gl_Cell_A(i, j, k)
				int2_coord(:, gl_Cell_A(i, j, k)) = gl_Cell_center(:, gl_Cell_A(i, j, k))
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2 
			do i = par_n_TS, par_n_HP - 1
				int2_Point_A(i + 3, j + 1, k) = gl_Cell_A(i, j, k)
				int2_coord(:, gl_Cell_A(i, j, k)) = gl_Cell_center(:, gl_Cell_A(i, j, k))
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2
			do i = par_n_HP, par_n_BS - 1
				int2_Point_A(i + 5, j + 1, k) = gl_Cell_A(i, j, k)
				int2_coord(:, gl_Cell_A(i, j, k)) = gl_Cell_center(:, gl_Cell_A(i, j, k))
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2
			do i = par_n_BS, par_n_END - 1
				int2_Point_A(i + 5, j + 1, k) = gl_Cell_A(i, j, k)
				int2_coord(:, gl_Cell_A(i, j, k)) = gl_Cell_center(:, gl_Cell_A(i, j, k))
			end do
		end do
	end do
	
	
	
	! B - ячеейки !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	N1 = size(gl_Cell_B(:, 1, 1))
	N2 = size(gl_Cell_B(1, :, 1))
	N3 = size(gl_Cell_B(1, 1, :))
	
	do k = 1, N3
		do j = 1, N2
			do i = 1, par_n_TS - 1
				int2_Point_B(i + 1, j + 1, k) = gl_Cell_B(i, j, k)
				int2_coord(:, gl_Cell_B(i, j, k)) = gl_Cell_center(:, gl_Cell_B(i, j, k))
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2
			do i = par_n_TS, N1
				int2_Point_B(i + 3, j + 1, k) = gl_Cell_B(i, j, k)
				int2_coord(:, gl_Cell_B(i, j, k)) = gl_Cell_center(:, gl_Cell_B(i, j, k))
			end do
		end do
	end do
	
	
	! C - точки  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	N1 = size(gl_Cell_C(:, 1, 1))
	N2 = size(gl_Cell_C(1, :, 1))
	N3 = size(gl_Cell_C(1, 1, :))
	
	do k = 1, N3
		do j = 1, N2
			do i = 1, par_n_HP - par_n_TS
				int2_Point_C(i, j, k) = gl_Cell_C(i, j, k)
				ijk = gl_Cell_C(i, j, k)
				int2_coord(:, gl_Cell_C(i, j, k)) = gl_Cell_center(:, gl_Cell_C(i, j, k))
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2
			do i = par_n_HP - par_n_TS + 1, par_n_BS - par_n_TS
				int2_Point_C(i + 2, j, k) = gl_Cell_C(i, j, k)
				int2_coord(:, gl_Cell_C(i, j, k)) = gl_Cell_center(:, gl_Cell_C(i, j, k))
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2
			do i = par_n_BS - par_n_TS + 1, par_n_END - par_n_TS
				int2_Point_C(i + 2, j, k) = gl_Cell_C(i, j, k)
				int2_coord(:, gl_Cell_C(i, j, k)) = gl_Cell_center(:, gl_Cell_C(i, j, k))
			end do
		end do
	end do
	
	
	! Заполняем особые точки (поверхности разрывов)
	N = size(gl_Cell_C) + size(gl_Cell_A) + size(gl_Cell_B)
	int2_Point_A(1, :, :) = N + 1
	int2_Point_B(1, :, :) = N + 1
	int2_coord(:, N + 1) = (/ 0.0, 0.0, 0.0 /)
	N = N + 1
	
	N1 = size(gl_Cell_A(:, 1, 1))
	N2 = size(gl_Cell_A(1, :, 1))
	N3 = size(gl_Cell_A(1, 1, :))
	
	! A
	do k = 1, N3
		do j = 1, N2
			int2_Point_A(par_n_TS + 1, j + 1, k) = N + 1
			int2_coord(:, N + 1) = gl_Gran_center(:, gl_Cell_gran(1, gl_Cell_A(par_n_TS - 1, j, k) ))
			N = N + 1
			
			int2_Point_A(par_n_TS + 2, j + 1, k) = N + 1
			int2_coord(:, N + 1) = gl_Gran_center(:, gl_Cell_gran(1, gl_Cell_A(par_n_TS - 1, j, k) ))
			N = N + 1
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2
			int2_Point_A(par_n_HP + 3, j + 1, k) = N + 1
			int2_coord(:, N + 1) = gl_Gran_center(:, gl_Cell_gran(1, gl_Cell_A(par_n_HP - 1, j, k) ))
			N = N + 1
			
			int2_Point_A(par_n_HP + 4, j + 1, k) = N + 1
			int2_coord(:, N + 1) = gl_Gran_center(:, gl_Cell_gran(1, gl_Cell_A(par_n_HP - 1, j, k) ))
			N = N + 1
		end do
	end do
	
	!B
	N1 = size(gl_Cell_B(:, 1, 1))
	N2 = size(gl_Cell_B(1, :, 1))
	N3 = size(gl_Cell_B(1, 1, :))
	
	do k = 1, N3
		do j = 1, N2
			int2_Point_B(par_n_TS + 1, j + 1, k) = N + 1
			int2_coord(:, N + 1) = gl_Gran_center(:, gl_Cell_gran(1, gl_Cell_B(par_n_TS - 1, j, k) ))
			N = N + 1
			
			int2_Point_B(par_n_TS + 2, j + 1, k) = N + 1
			int2_coord(:, N + 1) = gl_Gran_center(:, gl_Cell_gran(1, gl_Cell_B(par_n_TS - 1, j, k) ))
			N = N + 1 
		end do
	end do
	
	!C
	N1 = size(gl_Cell_C(:, 1, 1))
	N2 = size(gl_Cell_C(1, :, 1))
	N3 = size(gl_Cell_C(1, 1, :))
	
	do k = 1, N3
		do j = 1, N2
			int2_Point_C(par_n_HP - par_n_TS + 1, j, k) = N + 1
			int2_coord(:, N + 1) = gl_Gran_center(:, gl_Cell_gran(1, gl_Cell_C(par_n_HP - par_n_TS, j, k) ))
			N = N + 1
			
			int2_Point_C(par_n_HP - par_n_TS + 2, j, k) = N + 1
			int2_coord(:, N + 1) = gl_Gran_center(:, gl_Cell_gran(1, gl_Cell_C(par_n_HP - par_n_TS, j, k) ))
			N = N + 1 
		end do
	end do
	
	! Ячейки на оси
	N1 = size(int2_Point_A(:, 1, 1))
	N2 = size(int2_Point_A(1, :, 1))
	N3 = size(int2_Point_A(1, 1, :))
	
	!A
	do i = 2, N1
		int2_Point_A(i, 1, :) = N + 1
		int2_coord(:, N + 1) = (/ int2_coord(1, int2_Point_A(i, 2, 1)), 0.0_8, 0.0_8 /)
		N = N + 1
	end do
	
	N1 = size(int2_Point_B(:, 1, 1))
	N2 = size(int2_Point_B(1, :, 1))
	N3 = size(int2_Point_B(1, 1, :))
	
	!B
	do i = 2, N1
		int2_Point_B(i, 1, :) = N + 1
		int2_coord(:, N + 1) = (/ int2_coord(1, int2_Point_B(i, 2, 1)), 0.0_8, 0.0_8 /)
		N = N + 1
	end do
	
	! Заполняем ячейки (двойственную сетку)
	
	N = 1
	
	! A - точки 
	N1 = size(int2_Point_A(:, 1, 1))
	N2 = size(int2_Point_A(1, :, 1))
	N3 = size(int2_Point_A(1, 1, :))
	
	do k = 1, N3
		
		kk = k + 1
		if (kk > N3) kk = 1
		
		do j = 1, N2 - 1
			do i = 1, N1 - 1
				
				if (i == par_n_TS + 1) CYCLE
				if (i == par_n_HP + 3) CYCLE
				
				int2_Cell_A(i, j, k) = N
				
				int2_all_Cell(1, N) = int2_Point_A(i, j, k)
				int2_all_Cell(2, N) = int2_Point_A(i + 1, j, k)
				int2_all_Cell(3, N) = int2_Point_A(i + 1, j + 1, k)
				int2_all_Cell(4, N) = int2_Point_A(i, j + 1, k)
				int2_all_Cell(5, N) = int2_Point_A(i, j, kk)
				int2_all_Cell(6, N) = int2_Point_A(i + 1, j, kk)
				int2_all_Cell(7, N) = int2_Point_A(i + 1, j + 1, kk)
				int2_all_Cell(8, N) = int2_Point_A(i, j + 1, kk)
				
				int2_Cell_center(:, N) = 0.0
				do m = 1, 8
					ijk = int2_all_Cell(m, N)
					int2_Cell_center(:, N) = int2_Cell_center(:, N) + int2_coord(:, int2_all_Cell(m, N))
				end do
				int2_Cell_center(:, N) = int2_Cell_center(:, N)/8.0
				
				N = N + 1
			end do
		end do
	end do
	
	
	! B - точки 
	N1 = size(int2_Point_B(:, 1, 1))
	N2 = size(int2_Point_B(1, :, 1))
	N3 = size(int2_Point_B(1, 1, :))
	
	do k = 1, N3
		
		kk = k + 1
		if (kk > N3) kk = 1
		
		do j = 1, N2 - 1
			do i = 1, N1 - 1
				
				if (i == par_n_TS + 1) CYCLE
				
				int2_Cell_B(i, j, k) = N
				int2_all_Cell(1, N) = int2_Point_B(i, j, k)
				int2_all_Cell(4, N) = int2_Point_B(i + 1, j, k)
				int2_all_Cell(3, N) = int2_Point_B(i + 1, j + 1, k)
				int2_all_Cell(2, N) = int2_Point_B(i, j + 1, k)
				int2_all_Cell(5, N) = int2_Point_B(i, j, kk)
				int2_all_Cell(8, N) = int2_Point_B(i + 1, j, kk)
				int2_all_Cell(7, N) = int2_Point_B(i + 1, j + 1, kk)
				int2_all_Cell(6, N) = int2_Point_B(i, j + 1, kk)
				
				int2_Cell_center(:, N) = 0.0
				do m = 1, 8
					ijk = int2_all_Cell(m, N)
					int2_Cell_center(:, N) = int2_Cell_center(:, N) + int2_coord(:, int2_all_Cell(m, N))
				end do
				int2_Cell_center(:, N) = int2_Cell_center(:, N)/8.0
				
				N = N + 1
			end do
		end do
	end do
	
	! C - точки 
	N1 = size(int2_Point_C(:, 1, 1))
	N2 = size(int2_Point_C(1, :, 1))
	N3 = size(int2_Point_C(1, 1, :))
	
	do k = 1, N3
		
		kk = k + 1
		if (kk > N3) kk = 1
		
		do j = 1, N2 - 1
			do i = 1, N1 - 1
				
				if (i == par_n_HP - par_n_TS + 1) CYCLE
				
				int2_Cell_C(i + 1, j, k) = N
				int2_all_Cell(1, N) = int2_Point_C(i, j, k)
				int2_all_Cell(2, N) = int2_Point_C(i + 1, j, k)
				int2_all_Cell(3, N) = int2_Point_C(i + 1, j + 1, k)
				int2_all_Cell(4, N) = int2_Point_C(i, j + 1, k)
				int2_all_Cell(5, N) = int2_Point_C(i, j, kk)
				int2_all_Cell(6, N) = int2_Point_C(i + 1, j, kk)
				int2_all_Cell(7, N) = int2_Point_C(i + 1, j + 1, kk)
				int2_all_Cell(8, N) = int2_Point_C(i, j + 1, kk)
				
				int2_Cell_center(:, N) = 0.0
				do m = 1, 8
					int2_Cell_center(:, N) = int2_Cell_center(:, N) + int2_coord(:, int2_all_Cell(m, N))
				end do
				int2_Cell_center(:, N) = int2_Cell_center(:, N)/8.0
				
				N = N + 1
			end do
		end do
	end do
	
	!Залатываем дырки между группами !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	! AB  - точки 
	N1 = size(int2_Point_A(:, 1, 1))
	N2 = size(int2_Point_A(1, :, 1))
	N3 = size(int2_Point_A(1, 1, :))
	
	do k = 1, N3
		
		kk = k + 1
		if (kk > N3) kk = 1
		
		do j = N2, N2
			do i = 1, par_n_TS
				
				if (i == par_n_TS + 1) CYCLE
				if (i == par_n_HP + 3) CYCLE
				
				int2_Cell_A(i, j, k) = N
				int2_all_Cell(1, N) = int2_Point_A(i, j, k)
				int2_all_Cell(2, N) = int2_Point_A(i + 1, j, k)
				int2_all_Cell(3, N) = int2_Point_B(i + 1, size(int2_Point_B(1, :, 1)), k)
				int2_all_Cell(4, N) = int2_Point_B(i, size(int2_Point_B(1, :, 1)), k)
				int2_all_Cell(5, N) = int2_Point_A(i, j, kk)
				int2_all_Cell(6, N) = int2_Point_A(i + 1, j, kk)
				int2_all_Cell(7, N) = int2_Point_B(i + 1, size(int2_Point_B(1, :, 1)), kk)
				int2_all_Cell(8, N) = int2_Point_B(i, size(int2_Point_B(1, :, 1)), kk)
				
				int2_Cell_center(:, N) = 0.0
				do m = 1, 8
					int2_Cell_center(:, N) = int2_Cell_center(:, N) + int2_coord(:, int2_all_Cell(m, N))
				end do
				int2_Cell_center(:, N) = int2_Cell_center(:, N)/8.0
				
				N = N + 1
			end do
		end do
	end do
	
	! AC - точки 
	N1 = size(int2_Point_A(:, 1, 1))
	N2 = size(int2_Point_A(1, :, 1))
	N3 = size(int2_Point_A(1, 1, :))
	
	do k = 1, N3
		
		kk = k + 1
		if (kk > N3) kk = 1
		
		do j = N2, N2
			do i = par_n_TS + 2, N1 - 1
				
				if (i == par_n_TS + 1) CYCLE
				if (i == par_n_HP + 3) CYCLE
				
				int2_Cell_A(i, j, k) = N
				int2_all_Cell(1, N) = int2_Point_A(i, j, k)
				int2_all_Cell(2, N) = int2_Point_A(i + 1, j, k)
				int2_all_Cell(3, N) = int2_Point_C(i - par_n_TS - 1, 1, k)
				if (i == par_n_TS + 2) then
					int2_all_Cell(4, N) = int2_Point_B(i + 1, size(int2_Point_B(1, :, 1)), k)
				else
					int2_all_Cell(4, N) = int2_Point_C(i - 1 - par_n_TS	- 1, 1, k)
				end if
				int2_all_Cell(5, N) = int2_Point_A(i, j, kk)
				int2_all_Cell(6, N) = int2_Point_A(i + 1, j, kk)
				int2_all_Cell(7, N) = int2_Point_C(i - par_n_TS - 1, 1, kk)
				if (i == par_n_TS + 2) then
					int2_all_Cell(8, N) = int2_Point_B(i + 1, size(int2_Point_B(1, :, 1)), kk)
				else
					int2_all_Cell(8, N) = int2_Point_C(i - 1 - par_n_TS	- 1, 1, kk)
				end if
				
				int2_Cell_center(:, N) = 0.0
				do m = 1, 8
					int2_Cell_center(:, N) = int2_Cell_center(:, N) + int2_coord(:, int2_all_Cell(m, N))
				end do
				int2_Cell_center(:, N) = int2_Cell_center(:, N)/8.0
				
				N = N + 1
			end do
		end do
	end do
	
	N1 = size(int2_Point_A(:, 1, 1))
	N2 = size(int2_Point_A(1, :, 1))
	N3 = size(int2_Point_A(1, 1, :))
	
	! Около тройной точки
	!do k = 1, N3
	!	
	!	kk = k + 1
	!	if (kk > N3) kk = 1
	!	
	!	do j = N2, N2
	!		do i = par_n_TS + 2, par_n_TS + 2
	!			
	!			int2_Cell_A(i, j, k) = N
	!			int2_all_Cell(1, N) = int2_Point_A(i, j, k)
	!			int2_all_Cell(2, N) = int2_Point_A(i + 1, j, k)
	!			int2_all_Cell(3, N) = int2_Point_C(1, 1, k)
	!			int2_all_Cell(4, N) = int2_Point_B(i + 1, size(int2_Point_B(1, :, 1)), k)
	!			int2_all_Cell(5, N) = int2_Point_A(i, j, kk)
	!			int2_all_Cell(6, N) = int2_Point_A(i + 1, j, kk)
	!			int2_all_Cell(7, N) = int2_Point_C(1, 1, kk)
	!			int2_all_Cell(8, N) = int2_Point_B(i + 1, size(int2_Point_B(1, :, 1)), kk)
	!			
	!			int2_Cell_center(:, N) = 0.0
	!			do m = 1, 8
	!				int2_Cell_center(:, N) = int2_Cell_center(:, N) + int2_coord(:, int2_all_Cell(m, N))
	!			end do
	!			int2_Cell_center(:, N) = int2_Cell_center(:, N)/8.0
	!			
	!			N = N + 1
	!		end do
	!	end do
	!end do
	
	do k = 1, N3
		
		kk = k + 1
		if (kk > N3) kk = 1
		
		do j = N2, N2
			do i = par_n_TS + 1, par_n_TS + 1
				
				int2_Cell_A(i, j, k) = N
				int2_all_Cell(1, N) = int2_Point_A(i, j, k)
				int2_all_Cell(2, N) = int2_Point_A(i, j, k)
				int2_all_Cell(3, N) = int2_Point_B(i + 1, size(int2_Point_B(1, :, 1)), k)
				int2_all_Cell(4, N) = int2_Point_B(i, size(int2_Point_B(1, :, 1)), k)
				int2_all_Cell(5, N) = int2_Point_A(i, j, kk)
				int2_all_Cell(6, N) = int2_Point_A(i, j, kk)
				int2_all_Cell(7, N) = int2_Point_B(i + 1, size(int2_Point_B(1, :, 1)), kk)
				int2_all_Cell(8, N) = int2_Point_B(i, size(int2_Point_B(1, :, 1)), kk)
				
				int2_Cell_center(:, N) = 0.0
				do m = 1, 8
					int2_Cell_center(:, N) = int2_Cell_center(:, N) + int2_coord(:, int2_all_Cell(m, N))
				end do
				int2_Cell_center(:, N) = int2_Cell_center(:, N)/8.0
				
				N = N + 1
			end do
		end do
	end do
	
	! Между B и C точками
	
	N1 = size(int2_Point_C(:, 1, 1))
	N2 = size(int2_Point_C(1, :, 1))
	N3 = size(int2_Point_C(1, 1, :))
	
	do k = 1, N3
		
		kk = k + 1
		if (kk > N3) kk = 1
		
		do j = 1, N2 - 1
			do i = 1, 1
				
				int2_Cell_C(i, j, k) = N
				int2_all_Cell(1, N) = int2_Point_B(par_n_TS + 2 + j, size(int2_Point_B(1, :, 1)), k)
				int2_all_Cell(2, N) = int2_Point_C(i, j , k)
				int2_all_Cell(3, N) = int2_Point_C(i, j + 1, k)
				int2_all_Cell(4, N) = int2_Point_B(par_n_TS + 3 + j, size(int2_Point_B(1, :, 1)), k)
				int2_all_Cell(5, N) = int2_Point_B(par_n_TS + 2 + j, size(int2_Point_B(1, :, 1)), kk)
				int2_all_Cell(6, N) = int2_Point_C(i, j , kk)
				int2_all_Cell(7, N) = int2_Point_C(i, j + 1, kk)
				int2_all_Cell(8, N) = int2_Point_B(par_n_TS + 3 + j, size(int2_Point_B(1, :, 1)), kk)
				
				int2_Cell_center(:, N) = 0.0
				do m = 1, 8
					int2_Cell_center(:, N) = int2_Cell_center(:, N) + int2_coord(:, int2_all_Cell(m, N))
				end do
				int2_Cell_center(:, N) = int2_Cell_center(:, N)/8.0
				
				N = N + 1
			end do
		end do
	end do
	
	!print*, " A  ==  B?", N - 1, size(int2_all_Cell(1, :))
	
	! Прописываем связи и соседей всем ячейкам
	
	! связи Для A - ячеек ---------------------------------------------------------------------------------------------------
	N1 = size(int2_Cell_A(:, 1, 1))
	N2 = size(int2_Cell_A(1, :, 1))
	N3 = size(int2_Cell_A(1, 1, :))
	
	do k = 1, N3
		do j = 1, N2
			do i = 1, N1 - 1
				int2_all_neighbours(1, int2_Cell_A(i, j, k)) = int2_Cell_A(i + 1, j, k)
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2
			do i = 2, N1
				int2_all_neighbours(2, int2_Cell_A(i, j, k)) = int2_Cell_A(i - 1, j, k)
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2 - 1
			do i = 1, N1
				int2_all_neighbours(3, int2_Cell_A(i, j, k)) = int2_Cell_A(i, j + 1, k)
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 2, N2
			do i = 1, N1
				int2_all_neighbours(4, int2_Cell_A(i, j, k)) = int2_Cell_A(i, j - 1, k)
			end do
		end do
	end do
	
	do k = 1, N3
		kk = k + 1
		if (kk > N3) kk = 1
		
		do j = 1, N2
			do i = 1, N1
				int2_all_neighbours(5, int2_Cell_A(i, j, k)) = int2_Cell_A(i, j, kk)
			end do
		end do
	end do
	
	do k = 1, N3
		kk = k - 1
		if (kk < 1) kk = N3
		
		do j = 1, N2
			do i = 1, N1
				int2_all_neighbours(6, int2_Cell_A(i, j, k)) = int2_Cell_A(i, j, kk)
			end do
		end do
	end do
	
	! связи Для B - ячеек ---------------------------------------------------------------------------------------------------
	N1 = size(int2_Cell_B(:, 1, 1))
	N2 = size(int2_Cell_B(1, :, 1))
	N3 = size(int2_Cell_B(1, 1, :))
	
	do k = 1, N3
		do j = 1, N2
			do i = 1, N1 - 1
				int2_all_neighbours(1, int2_Cell_B(i, j, k)) = int2_Cell_B(i + 1, j, k)
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2
			do i = 2, N1
				int2_all_neighbours(2, int2_Cell_B(i, j, k)) = int2_Cell_B(i - 1, j, k)
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2 - 1
			do i = 1, N1
				int2_all_neighbours(3, int2_Cell_B(i, j, k)) = int2_Cell_B(i, j + 1, k)
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 2, N2
			do i = 1, N1
				int2_all_neighbours(4, int2_Cell_B(i, j, k)) = int2_Cell_B(i, j - 1, k)
			end do
		end do
	end do
	
	do k = 1, N3
		kk = k + 1
		if (kk > N3) kk = 1
		
		do j = 1, N2
			do i = 1, N1
				int2_all_neighbours(5, int2_Cell_B(i, j, k)) = int2_Cell_B(i, j, kk)
			end do
		end do
	end do
	
	do k = 1, N3
		kk = k - 1
		if (kk < 1) kk = N3
		
		do j = 1, N2
			do i = 1, N1
				int2_all_neighbours(6, int2_Cell_B(i, j, k)) = int2_Cell_B(i, j, kk)
			end do
		end do
	end do
	
	! связи Для C - ячеек ---------------------------------------------------------------------------------------------------
	N1 = size(int2_Cell_C(:, 1, 1))
	N2 = size(int2_Cell_C(1, :, 1))
	N3 = size(int2_Cell_C(1, 1, :))
	
	do k = 1, N3
		do j = 1, N2
			do i = 1, N1 - 1
				int2_all_neighbours(1, int2_Cell_C(i, j, k)) = int2_Cell_C(i + 1, j, k)
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2
			do i = 2, N1
				int2_all_neighbours(2, int2_Cell_C(i, j, k)) = int2_Cell_C(i - 1, j, k)
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2 - 1
			do i = 1, N1
				int2_all_neighbours(3, int2_Cell_C(i, j, k)) = int2_Cell_C(i, j + 1, k)
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 2, N2
			do i = 1, N1
				int2_all_neighbours(4, int2_Cell_C(i, j, k)) = int2_Cell_C(i, j - 1, k)
			end do
		end do
	end do
	
	do k = 1, N3
		kk = k + 1
		if (kk > N3) kk = 1
		
		do j = 1, N2
			do i = 1, N1
				int2_all_neighbours(5, int2_Cell_C(i, j, k)) = int2_Cell_C(i, j, kk)
			end do
		end do
	end do
	
	do k = 1, N3
		kk = k - 1
		if (kk < 1) kk = N3
		
		do j = 1, N2
			do i = 1, N1
				int2_all_neighbours(6, int2_Cell_C(i, j, k)) = int2_Cell_C(i, j, kk)
			end do
		end do
	end do
	
	! Теперь нужно отдельно прописать связи на стыках групп + на поверхностях разрыва + в тройной точке
	
	
	end subroutine Int2_Initial
	
	subroutine Int2_Print_center()
		implicit none
		integer(4) :: i, j, k, N1, N2, N3
		
		open(1, file = 'Int2_center.txt')
		
		N1 = size(int2_Cell_A(:, 1, 1))
		N2 = size(int2_Cell_A(1, :, 1))
		N3 = size(int2_Cell_A(1, 1, :))
		
		do k = 1, 1
			do j = 1, N2
				do i = 1, N1
					if (int2_Cell_A(i, j, k) <= 0) CYCLE
					write(1,*) int2_Cell_center(:, int2_Cell_A(i, j, k))
				end do
			end do
		end do
		
		N1 = size(int2_Cell_B(:, 1, 1))
		N2 = size(int2_Cell_B(1, :, 1))
		N3 = size(int2_Cell_B(1, 1, :))
		
		do k = 1, 1
			do j = 1, N2
				do i = 1, N1
					if (int2_Cell_B(i, j, k) <= 0) CYCLE
					write(1,*) int2_Cell_center(:, int2_Cell_B(i, j, k))
				end do
			end do
		end do
		
		N1 = size(int2_Cell_C(:, 1, 1))
		N2 = size(int2_Cell_C(1, :, 1))
		N3 = size(int2_Cell_C(1, 1, :))
		
		do k = 1, 1
			do j = 1, N2
				do i = 1, N1
					if (int2_Cell_C(i, j, k) <= 0) CYCLE
					write(1,*) int2_Cell_center(:, int2_Cell_C(i, j, k))
				end do
			end do
		end do
		
		
		close(1)
	
	end subroutine Int2_Print_center
	
	
	subroutine Int2_Print_setka_2()
	! печать двойственной сетки
		implicit none
		integer(4) :: i, j, k, N1, N2, N3, N
		
		N = size(int2_Cell_A(1, :, 1)) * size(int2_Cell_A(:, 1, 1)) + size(int2_Cell_B(1, :, 1)) * size(int2_Cell_B(:, 1, 1)) + &
			size(int2_Cell_C(1, :, 1)) * size(int2_Cell_C(:, 1, 1))
		
		open(1, file = 'Int2_setka.txt')
		write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z', 'Q'  ZONE T= 'HP', N= ", 4 * N, ", E =  ", 4 * N , ", F=FEPOINT, ET=LINESEG "
		
		N1 = size(int2_Cell_A(:, 1, 1))
		N2 = size(int2_Cell_A(1, :, 1))
		N3 = size(int2_Cell_A(1, 1, :))
		
		do k = 1, 1
			do j = 1, N2
				do i = 1, N1
					if (int2_Cell_A(i, j, k) <= 0) then
						write(1,*) 0.0, 0.0, 0.0, 1
						write(1,*) 0.0, 0.0, 0.0, 1
						write(1,*) 0.0, 0.0, 0.0, 1
						write(1,*) 0.0, 0.0, 0.0, 1
					else
						write(1,*) int2_coord(:, int2_all_Cell(1, int2_Cell_A(i, j, k))), 1
						write(1,*) int2_coord(:, int2_all_Cell(2, int2_Cell_A(i, j, k))), 1
						write(1,*) int2_coord(:, int2_all_Cell(3, int2_Cell_A(i, j, k))), 1
						write(1,*) int2_coord(:, int2_all_Cell(4, int2_Cell_A(i, j, k))), 1
					end if
					
				end do
			end do
		end do
		
		N1 = size(int2_Cell_B(:, 1, 1))
		N2 = size(int2_Cell_B(1, :, 1))
		N3 = size(int2_Cell_B(1, 1, :))
		
		do k = 1, 1
			do j = 1, N2
				do i = 1, N1
					if (int2_Cell_B(i, j, k) <= 0) then
						write(1,*) 0.0, 0.0, 0.0, 2
						write(1,*) 0.0, 0.0, 0.0, 2
						write(1,*) 0.0, 0.0, 0.0, 2
						write(1,*) 0.0, 0.0, 0.0, 2
					else
						write(1,*) int2_coord(:, int2_all_Cell(1, int2_Cell_B(i, j, k))), 2
						write(1,*) int2_coord(:, int2_all_Cell(2, int2_Cell_B(i, j, k))), 2
						write(1,*) int2_coord(:, int2_all_Cell(3, int2_Cell_B(i, j, k))), 2
						write(1,*) int2_coord(:, int2_all_Cell(4, int2_Cell_B(i, j, k))), 2
					end if
				end do
			end do
		end do
		
		N1 = size(int2_Cell_C(:, 1, 1))
		N2 = size(int2_Cell_C(1, :, 1))
		N3 = size(int2_Cell_C(1, 1, :))
		
		do k = 1, 1
			do j = 1, N2
				do i = 1, N1
					if (int2_Cell_C(i, j, k) <= 0) then
						write(1,*) 0.0, 0.0, 0.0, 3
						write(1,*) 0.0, 0.0, 0.0, 3
						write(1,*) 0.0, 0.0, 0.0, 3
						write(1,*) 0.0, 0.0, 0.0, 3
					else
						write(1,*) int2_coord(:, int2_all_Cell(1, int2_Cell_C(i, j, k))), 3
						write(1,*) int2_coord(:, int2_all_Cell(2, int2_Cell_C(i, j, k))), 3
						write(1,*) int2_coord(:, int2_all_Cell(3, int2_Cell_C(i, j, k))), 3
						write(1,*) int2_coord(:, int2_all_Cell(4, int2_Cell_C(i, j, k))), 3
					end if
				end do
			end do
		end do
		
		! Connectivity list
    do j = 0, N
        write(1,*) 4 * j + 1, 4 * j + 2
        write(1,*) 4 * j + 2, 4 * j + 3
        write(1,*) 4 * j + 3, 4 * j + 4
        write(1,*) 4 * j + 4, 4 * j + 1
    end do
		
		
		close(1)
	
	end subroutine Int2_Print_setka_2
	

	subroutine Int2_Print_point_plane()
	
	integer(4) :: i, j, k, N1, N2, N3
	
	open(1, file = 'Int2_print.txt')
	write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z' , 'T' ZONE T= 'HP'"
	
	N1 = size(int2_Point_A(:, 1, 1))
	N2 = size(int2_Point_A(1, :, 1))
	N3 = size(int2_Point_A(1, 1, :))
	
	do k = 1, 1
		do j = 1, N2
			do i = 1, N1
				if( int2_Point_A(i, j, k) > 0 ) write(1,*) int2_coord(:, int2_Point_A(i, j, k)), 1
			end do
		end do
	end do
	
	!
	N1 = size(int2_Point_B(:, 1, 1))
	N2 = size(int2_Point_B(1, :, 1))
	N3 = size(int2_Point_B(1, 1, :))
	
	do k = 1, 1
		do j = 1, N2
			do i = 1, N1
				if( int2_Point_B(i, j, k) > 0 ) write(1,*) int2_coord(:, int2_Point_B(i, j, k)), 2
			end do
		end do
	end do
	!
	!
	N1 = size(int2_Point_C(:, 1, 1))
	N2 = size(int2_Point_C(1, :, 1))
	N3 = size(int2_Point_C(1, 1, :))
	
	do k = 1, 1
		do j = 1, N2
			do i = 1, N1
				if( int2_Point_C(i, j, k) > 0 ) write(1,*) int2_coord(:, int2_Point_C(i, j, k)), 3
			end do
		end do
	end do
	
	
	close(1)
	
	end subroutine Int2_Print_point_plane
	
	
	subroutine Int2_Set_Interpolate()
	
	implicit none
	integer :: n
	! Выделения памяти и начальная инициализация
	allocate(int2_Point_A( size(gl_Cell_A(:, 1, 1)) + 5, size(gl_Cell_A(1, :, 1)) + 1, size(gl_Cell_A(1, 1, :)) ))
	allocate(int2_Point_B( size(gl_Cell_B(:, 1, 1)) + 3, size(gl_Cell_B(1, :, 1)) + 1, size(gl_Cell_B(1, 1, :)) ))
	allocate( int2_Point_C( size(gl_Cell_C(:, 1, 1)) + 2, size(gl_Cell_C(1, :, 1)), size(gl_Cell_C(1, 1, :)) ) )
	
	n = size(gl_Cell_A(:, 1, 1)) + 2
	n = n * (size(gl_Cell_A(1, :, 1)) + 1) * size(gl_Cell_A(1, 1, :))
	n = n + (size(gl_Cell_B(:, 1, 1)) + 1) * size(gl_Cell_B(1, :, 1)) * size(gl_Cell_B(1, 1, :))
	n = n + (size(gl_Cell_C(:, 1, 1)) + 1) * (size(gl_Cell_C(1, :, 1)) - 1) * size(gl_Cell_C(1, 1, :))
	n = n + size(int2_Point_A(1, 1, :))
	
	allocate(int2_all_Cell(8, n ))  ! Добавили ячейку около тройной точки
	allocate(int2_coord(3, size(int2_Point_A) + size(int2_Point_B) + size(int2_Point_C)))
	allocate(int2_Cell_center(3, n))
	
	allocate(int2_Cell_A(size(gl_Cell_A(:, 1, 1)) + 4, size(gl_Cell_A(1, :, 1)) + 1, size(gl_Cell_A(1, 1, :))))
	allocate(int2_Cell_B(size(gl_Cell_B(:, 1, 1)) + 2, size(gl_Cell_B(1, :, 1)), size(gl_Cell_B(1, 1, :))))
	allocate(int2_Cell_C(size(gl_Cell_C(:, 1, 1)) + 2, size(gl_Cell_C(1, :, 1)) - 1, size(gl_Cell_C(1, 1, :))))
	
	allocate(int2_all_neighbours(6, size(int2_all_Cell(1, :)) ))
	
	int2_coord = 0.0
	int2_all_Cell = -100
	int2_Point_A = -100
	int2_Point_B = -100
	int2_Cell_center(:, :) = 0.0
	int2_Cell_A = 0
	int2_Cell_B = 0
	int2_Cell_C = 0
	int2_all_neighbours = 0 
	
	end subroutine Int2_Set_Interpolate
	
	
	end module Interpolate2
	
	module Interpolate                     ! Модуль геометрических - сеточных параметров - констант программы
	! Модуль НЕ АКТУАЛЕН (интерполяция примитивна, но работает)
	! Интерполяция на основе основной сетки (покрытия пространства не полное, точки могут выпадать)
	! Точка может принадлежать ячейке, но не попадать не в один тетраэдр, из-за этого
	!    интерполированные значения могут быть отрицательны!
    
    implicit none
	
	! Переменные ниже в итоге не понадобились.
	! Но алгоритм их заполнения может быть полезен (только поэтому не удалял их)
	logical, allocatable :: Int1_set_gran(:)   ! Заполнены ли треугольные грани в этой ячейке?
	integer(4), allocatable :: Int1_grans_triangle(:, :, :)  ! (3, 12, :) Треугольные грани в каждой ячейке
	! 12 граней в каждой ячейке, каждая грань состоит из трёх узлов, : - количество ячеек в сетке
	real(8), allocatable :: Int1_grans_normal(:, :, :)  ! (4, 12, :) Внешняя нормаль каждой грани в ячейке + параметр D
	! n1 * x + n2 * y + n3 * z + D = 0  - уравнение грани - нормаль внешняя
    
	contains
	
	subroutine Set_Interpolate_main()
	! Выделяет память и заполняет массивы из основной программы
	
	USE GEO_PARAM
	USE STORAGE
	implicit none
	
	integer :: i, j, k, N1, N2, N3, s, s1, s2, s3
	real(8) :: aa(3), bb(3), normal(3), nn
	
	if (allocated(gl_Cell_gran_inter) == .True.) then
        STOP "Function Set_Interpolate vizvana neskolko raz!!! Dopustimo tolko 1 raz! _9876trtyuiokjhgfde45"    
	end if
	
	if (allocated(gl_Cell_gran) == .False.) then
        STOP "Ne opredeleny osnovnye massivy v programme! _76rfvbnjki87"    
    end if
	
	! Выделяем память
	allocate(gl_Cell_gran_inter, MOLD = gl_Cell_gran)
	allocate(gl_Gran_normal_inter, MOLD = gl_Gran_normal)
	allocate(gl_Gran_center_inter, MOLD = gl_Gran_center)
	allocate(gl_Gran_neighbour_inter, MOLD = gl_Gran_neighbour)
	allocate(gl_Cell_center_inter, MOLD = gl_Cell_center)
	allocate(gl_Cell_neighbour_inter, MOLD = gl_Cell_neighbour)
	allocate(gl_Cell_par_inter, MOLD = gl_Cell_par)
	allocate(gl_Cell_par_MF_inter, MOLD = gl_Cell_par_MF)
	allocate(gl_Cell_dist_inter, MOLD = gl_Cell_dist)
	
	allocate(Int1_grans_triangle(3, 12, size(gl_Cell_par(1, :))))
	allocate(Int1_grans_normal(4, 12, size(gl_Cell_par(1, :))))
	allocate(Int1_set_gran( size(gl_Cell_par(1, :)) ))
	
	gl_Cell_gran_inter = gl_Cell_gran
	gl_Gran_normal_inter = gl_Gran_normal
	gl_Gran_center_inter = gl_Gran_center
	gl_Gran_neighbour_inter = gl_Gran_neighbour
	gl_Cell_center_inter = gl_Cell_center
	gl_Cell_neighbour_inter = gl_Cell_neighbour
	gl_Cell_par_inter = gl_Cell_par
	gl_Cell_par_MF_inter = gl_Cell_par_MF
	gl_Cell_dist_inter = gl_Cell_dist
	Int1_set_gran = .False.
	
	
	! Заполняем треугольные грани
	N1 = size(gl_all_Cell(1, :))
	
	do k = 1, N1
		! 1
		s = gl_Cell_neighbour(1, k)  
		if(s <= 0 .or. Int1_set_gran(max(1, s)) == .False.) then
			Int1_grans_triangle(:, 1, k) = (/ 2, 3, 6  /)
			Int1_grans_triangle(:, 2, k) = (/ 3, 7, 6  /)
		else
			if ( ANY( Int1_grans_triangle(:, 3, s) == 4) .and. ANY( Int1_grans_triangle(:, 3, s) == 5) ) then
				Int1_grans_triangle(:, 1, k) = (/ 2, 3, 6  /)
				Int1_grans_triangle(:, 2, k) = (/ 3, 7, 6  /)
			else
				Int1_grans_triangle(:, 1, k) = (/ 2, 7, 6  /)
				Int1_grans_triangle(:, 2, k) = (/ 2, 3, 7  /)
			end if
		end if
		
		! 2
		s = gl_Cell_neighbour(2, k)  
		if(s <= 0 .or. Int1_set_gran(max(1, s)) == .False.) then
			Int1_grans_triangle(:, 3, k) = (/ 5, 4, 1  /)
			Int1_grans_triangle(:, 4, k) = (/ 5, 8, 4  /)
		else
			if ( ANY( Int1_grans_triangle(:, 1, s) == 3) .and. ANY( Int1_grans_triangle(:, 1, s) == 6) ) then
				Int1_grans_triangle(:, 3, k) = (/ 5, 4, 1  /)
				Int1_grans_triangle(:, 4, k) = (/ 5, 8, 4  /)
			else
				Int1_grans_triangle(:, 3, k) = (/ 8, 4, 1  /)
				Int1_grans_triangle(:, 4, k) = (/ 5, 8, 1  /)
			end if
		end if
		
		
		! 3
		s = gl_Cell_neighbour(3, k)  
		if(s <= 0 .or. Int1_set_gran(max(1, s)) == .False.) then
			Int1_grans_triangle(:, 5, k) = (/ 4, 7, 3  /)
			Int1_grans_triangle(:, 6, k) = (/ 4, 8, 7  /)
		else
			if ( ANY( Int1_grans_triangle(:, 7, s) == 1) .and. ANY( Int1_grans_triangle(:, 7, s) == 6) ) then
				Int1_grans_triangle(:, 5, k) = (/ 4, 7, 3  /)
				Int1_grans_triangle(:, 6, k) = (/ 4, 8, 7  /)
			else
				Int1_grans_triangle(:, 5, k) = (/ 3, 8, 7  /)
				Int1_grans_triangle(:, 6, k) = (/ 4, 8, 3  /)
			end if
		end if
		
		! 4
		s = gl_Cell_neighbour(4, k)  
		if(s <= 0 .or. Int1_set_gran(max(1, s)) == .False.) then
			Int1_grans_triangle(:, 7, k) = (/ 1, 6, 5  /)
			Int1_grans_triangle(:, 8, k) = (/ 1, 2, 6  /)
		else
			if ( ANY( Int1_grans_triangle(:, 5, s) == 4) .and. ANY( Int1_grans_triangle(:, 5, s) == 7) ) then
				Int1_grans_triangle(:, 7, k) = (/ 1, 6, 5  /)
				Int1_grans_triangle(:, 8, k) = (/ 1, 2, 6  /)
			else
				Int1_grans_triangle(:, 7, k) = (/ 2, 6, 5  /)
				Int1_grans_triangle(:, 8, k) = (/ 1, 2, 5  /)
			end if
		end if
		
		
		! 5
		s = gl_Cell_neighbour(5, k)  
		if(s <= 0 .or. Int1_set_gran(max(1, s)) == .False.) then
			Int1_grans_triangle(:, 9, k) = (/ 5, 6, 8  /)
			Int1_grans_triangle(:, 10, k) = (/ 6, 7, 8  /)
		else
			if ( ANY( Int1_grans_triangle(:, 11, s) == 4) .and. ANY( Int1_grans_triangle(:, 11, s) == 2) ) then
				Int1_grans_triangle(:, 9, k) = (/ 5, 6, 8  /)
				Int1_grans_triangle(:, 10, k) = (/ 6, 7, 8  /)
			else
				Int1_grans_triangle(:, 9, k) = (/ 5, 6, 7  /)
				Int1_grans_triangle(:, 10, k) = (/ 5, 7, 8  /)
			end if
		end if
		
		
		! 6
		s = gl_Cell_neighbour(6, k)  
		if(s <= 0 .or. Int1_set_gran(max(1, s)) == .False.) then
			Int1_grans_triangle(:, 11, k) = (/ 1, 4, 2  /)
			Int1_grans_triangle(:, 12, k) = (/ 2, 4, 3  /)
		else
			if ( ANY( Int1_grans_triangle(:, 9, s) == 6) .and. ANY( Int1_grans_triangle(:, 9, s) == 8) ) then
				Int1_grans_triangle(:, 11, k) = (/ 1, 4, 2  /)
				Int1_grans_triangle(:, 12, k) = (/ 2, 4, 3  /)
			else
				Int1_grans_triangle(:, 11, k) = (/ 1, 4, 3  /)
				Int1_grans_triangle(:, 12, k) = (/ 2, 1, 3  /)
			end if
		end if
		
		Int1_set_gran(k) = .True.
		
		!gl_all_Cell(:, k)
		
		do j = 1, 12
			s1 = Int1_grans_triangle(1, j, k)
			s2 = Int1_grans_triangle(2, j, k)
			s3 = Int1_grans_triangle(3, j, k)
			
			aa = (/ gl_x(gl_all_Cell(s2, k)), gl_y(gl_all_Cell(s2, k)), gl_z(gl_all_Cell(s2, k)) /) - (/ gl_x(gl_all_Cell(s1, k)), gl_y(gl_all_Cell(s1, k)), gl_z(gl_all_Cell(s1, k)) /)
			bb = (/ gl_x(gl_all_Cell(s3, k)), gl_y(gl_all_Cell(s3, k)), gl_z(gl_all_Cell(s3, k)) /) - (/ gl_x(gl_all_Cell(s1, k)), gl_y(gl_all_Cell(s1, k)), gl_z(gl_all_Cell(s1, k)) /)
		
			normal(1) = aa(2) * bb(3) - aa(3) * bb(2) 
			normal(2) = aa(3) * bb(1) - aa(1) * bb(3) 
			normal(3) = aa(1) * bb(2) - aa(2) * bb(1) 
			
			nn = norm2(normal)
			normal = normal/nn
			
			Int1_grans_normal(1:3, j, k) = normal
			Int1_grans_normal(4, j, k) = -DOT_PRODUCT(normal, (/ gl_x(gl_all_Cell(s2, k)), gl_y(gl_all_Cell(s2, k)), gl_z(gl_all_Cell(s2, k)) /))
			
			aa = ((/ gl_x(gl_all_Cell(s1, k)), gl_y(gl_all_Cell(s1, k)), gl_z(gl_all_Cell(s1, k)) /) + &
				(/ gl_x(gl_all_Cell(s2, k)), gl_y(gl_all_Cell(s2, k)), gl_z(gl_all_Cell(s2, k)) /) + &
				(/ gl_x(gl_all_Cell(s3, k)), gl_y(gl_all_Cell(s3, k)), gl_z(gl_all_Cell(s3, k)) /))/3.0
			
			bb = gl_Cell_center(:, k)
			
			!if(DOT_PRODUCT(normal, aa - bb) < 0.0) print*, "Error normal 987656789uyghjmnbgyu"
			
		end do
		
		
		! Надо сделать проверку нормалей и граней
		
		
	end do
	
	
	end subroutine Set_Interpolate_main
	
	subroutine Dell_Interpolate()
	! Удалить массивы интерполяции, если больше не нужны
	USE STORAGE
	implicit none
	
	
	DEALLOCATE(gl_Cell_gran_inter)
	DEALLOCATE(gl_Gran_normal_inter)
	DEALLOCATE(gl_Gran_center_inter)
	DEALLOCATE(gl_Gran_neighbour_inter)
	DEALLOCATE(gl_Cell_center_inter)
	DEALLOCATE(gl_Cell_neighbour_inter)
	DEALLOCATE(gl_Cell_par_inter)
	DEALLOCATE(gl_Cell_par_MF_inter)
	DEALLOCATE(Int1_grans_triangle)
	DEALLOCATE(Int1_grans_normal)
	DEALLOCATE(Int1_set_gran)
	
	
	end subroutine Dell_Interpolate
	
	subroutine Save_interpolate_bin(num)  ! Сохранение сетки в бинарном файле
    ! Variables
    use STORAGE
    use GEO_PARAM
    implicit none
    integer, intent(in) :: num
    character(len=5) :: name
    
    write(unit=name, fmt='(i5.5)') num
    
    open(1, file = "interpolate_save_" // name // ".bin", FORM = 'BINARY')
	
    
    write(1) 1
	
	write(1) size(gl_Cell_gran_inter(1, :))
    write(1) gl_Cell_gran_inter
	
	write(1) size(gl_Gran_normal_inter(1, :))
	write(1) gl_Gran_normal_inter
	
	write(1) size(gl_Gran_center_inter(1, :))
	write(1) gl_Gran_center_inter
	
	write(1) size(gl_Gran_neighbour_inter(1, :))
	write(1) gl_Gran_neighbour_inter
	
	write(1) size(gl_Cell_center_inter(1, :))
	write(1) gl_Cell_center_inter
	
	write(1) size(gl_Cell_neighbour_inter(1, :))
	write(1) gl_Cell_neighbour_inter
	
	write(1) size(gl_Cell_par_inter(1, :))
	write(1) gl_Cell_par_inter
	
	write(1) size(gl_Cell_par_MF_inter(1, 1, :))
	write(1) gl_Cell_par_MF_inter
	
    write(1) 1
	write(1) size(gl_Cell_dist_inter(:))
	write(1) gl_Cell_dist_inter
	
	write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    
    
    close(1)
	end subroutine Save_interpolate_bin
	
	subroutine Read_interpolate_bin(num)  
    ! Variables
    use STORAGE
    use GEO_PARAM
    implicit none
    integer, intent(in) :: num
    character(len=5) :: name
    integer :: n
    logical :: exists
    
    write(unit=name,fmt='(i5.5)') num
    
    inquire(file="interpolate_save_" // name // ".bin", exist=exists)
    
    if (exists == .False.) then
        STOP "net faila!!!"
    end if
    
    if (allocated(gl_Cell_gran_inter) == .True.) then
        STOP "Function Read_interpolate_bin vizvana neskolko raz!!! Dopustimo tolko 1 raz! 1qwedsdfr456tfghjnhu876"    
	end if
    
    
    open(1, file = "interpolate_save_" // name // ".bin", FORM = 'BINARY', ACTION = "READ")
    
	
    ! Для  будующих дополнений
    read(1) n
    if (n == 1) then
		
        read(1) n
		allocate(gl_Cell_gran_inter(6, n))
		read(1)	gl_Cell_gran_inter
		
		read(1) n
		allocate(gl_Gran_normal_inter(3, n))
		read(1)	gl_Gran_normal_inter
		
		read(1) n
		allocate(gl_Gran_center_inter(3, n))
		read(1)	gl_Gran_center_inter
		
		read(1) n
		allocate(gl_Gran_neighbour_inter(2, n))
		read(1)	gl_Gran_neighbour_inter
		
		read(1) n
		allocate(gl_Cell_center_inter(3, n))
		read(1)	gl_Cell_center_inter
		
		read(1) n
		allocate(gl_Cell_neighbour_inter(6, n))
		read(1)	gl_Cell_neighbour_inter
		
		read(1) n
		allocate(gl_Cell_par_inter(9, n))
		read(1)	gl_Cell_par_inter
		
		read(1) n
		allocate(gl_Cell_par_MF_inter(5, 4, n))
		read(1)	gl_Cell_par_MF_inter
		
	end if
    
	
    read(1) n
	if (n == 1) then
		read(1) n
		allocate(gl_Cell_dist_inter(n))
		read(1) gl_Cell_dist_inter
	end if
	
	read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    
    
    close(1)
	end subroutine Read_interpolate_bin
	
	subroutine Re_interpolate()
	! Изменяет значения переменных в текущей сетке из интерполированной
	USE STORAGE
	implicit none
	
	integer(4) :: i, zone, number, cell
	real(8) :: F(9), F_mf(5, 4)
	
	number = 1
	
	do i = 1, size(gl_Cell_par(1, :))
		call Interpolate_point(gl_Cell_center(1, i), gl_Cell_center(2, i), gl_Cell_center(3, i), F, F_mf, zone, cell, number)
		number = cell
		gl_Cell_par(:, i) = F
		gl_Cell_par_MF(:, :, i) = F_mf
	end do
	
	end subroutine Re_interpolate
	
	subroutine Interpolate_point(x, y, z, F, F_mf, zone, cell, number)
	! zone - в какой зоне находимся
	! cell - в какой ячейке находимся
	! number - в какой ячейке искать точку (если известно примерно)
	USE STORAGE
	implicit none
	real(8), intent(in) :: x, y, z
	integer(4), intent(in), optional :: number
	integer(4), intent(out) :: zone, cell
	real(8), intent(out) :: F(9), F_mf(5, 4)
	
	integer(4) :: num
	
	if (allocated(gl_Cell_gran_inter) == .False.) then
        STOP "Function Set_Interpolate NE vizvana  789ihj76uhgfrerfdcgt"    
	end if
	
	if(.not. present(number)) then
    num = 1 
    else 
    num = number
	end if
	
	zone = 1
	
	call Get_Cell_Interpolate(x, y, z, num)
	if(num == -1) then
		cell = -1
		return
	end if
	
	cell = num
	call Find_tetraedr_Interpolate(x, y, z, num, F, F_mf)
	
	end subroutine Interpolate_point
		
	subroutine Get_Cell_Interpolate(x, y, z, num) 
	! num - начальная ячейка поиска (если неизвестно, положить равной 1)
	USE GEO_PARAM
	USE STORAGE
	implicit none
	! Variables
	
	real(8), intent(in) :: x, y, z
	integer(4), intent(in out) :: num
	
	integer(4) :: n, i, gr, s, step
	real(8) :: normal(3), center(3), r(3), dd
	
	r(1) = x
	r(2) = y
	r(3) = z
	n = num
	step = 0
	
11	CONTINUE
	step = step + 1
	
	do i = 1, 6
		s = 2
		gr = gl_Cell_gran_inter(i,n)
		if (gr <= 0) CYCLE
		normal = gl_Gran_normal_inter(:, gr)
		center = gl_Gran_center_inter(:, gr)
		
		if(gl_Gran_neighbour_inter(1, gr) /= n) then
			normal = -normal
			s = 1
		end if
		
		if (gl_Gran_neighbour_inter(s, gr) <= 0) CYCLE
		
		dd = DOT_PRODUCT(normal, r - center)
		if( dd > 0.01 * gl_Cell_dist_inter(n)) then
			n = gl_Gran_neighbour_inter(s, gr)
			if (step > 10000) then
				!print*, "gr, n, i = ", gr, n, i
				!print*, "neybor ", gl_Gran_neighbour_inter(1, gr), gl_Gran_neighbour_inter(2, gr)
				!print*, r, dd
				!print*, center
				!print*, normal
				!pause
				print*, "10000 shagov 78987yhjuyhjiu"
				num = n
				EXIT
				!return
			end if
			
			GO TO 11 
		end if
	end do
	
	
	num = n
	
	end subroutine Get_Cell_Interpolate
	
	subroutine Find_tetraedr_Interpolate(x, y, z, num, F, F_mf)
	! num - номер ячейки, в которой находится точка
	USE GEO_PARAM
	USE STORAGE
	USE ieee_arithmetic
	implicit none
	! Variables
	real(8), intent(in) :: x, y, z
	integer(4), intent(in) :: num
	real(8), intent(out) :: F(9), F_mf(5, 4)
	
	integer(4) :: i, j, s1, s2, s3, i_best
	real(8) :: r(3), A(3), B(3), C(3), D(3), M(4, 4), det, aa(3), bb(3), normal(3), dist, dist2
	real(8), dimension(4, 4) :: Minv
	real(8), dimension(4) :: work  ! work array for LAPACK
	real(8), dimension(1, 4) :: vec
	integer, dimension(4) :: ipiv   ! pivot indices
	integer :: n, info
	
	! print*, x, y, z, "START"
	
	!$inter external DGETRF
	!$inter external DGETRI
	
	! print*, "START"
	
	r(1) = x
	r(2) = y
	r(3) = z
	
	A = gl_Cell_center_inter(:, num)
	
	
	do i = 1, 8
		s1 = gl_Cell_neighbour_inter( par_tet(1, i), num)
		s2 = gl_Cell_neighbour_inter( par_tet(2, i), num)
		s3 = gl_Cell_neighbour_inter( par_tet(3, i), num)
		
		if(s1 <= 0 .or. s2 <= 0 .or. s3 <= 0) CYCLE
		
		B = gl_Cell_center_inter(:, s1)
		C = gl_Cell_center_inter(:, s2)
		D = gl_Cell_center_inter(:, s3)
		
		! 1
		aa = B - C
		bb = D - C
		
		normal(1) = aa(2) * bb(3) - aa(3) * bb(2) 
		normal(2) = aa(3) * bb(1) - aa(1) * bb(3) 
		normal(3) = aa(1) * bb(2) - aa(2) * bb(1) 
		
		if ( DOT_PRODUCT(normal, r - C) * DOT_PRODUCT(normal, A - C) < 0) CYCLE
		
		! 2
		aa = C - A
		bb = D - A
		
		normal(1) = aa(2) * bb(3) - aa(3) * bb(2) 
		normal(2) = aa(3) * bb(1) - aa(1) * bb(3) 
		normal(3) = aa(1) * bb(2) - aa(2) * bb(1) 
		
		if ( DOT_PRODUCT(normal, r - A) * DOT_PRODUCT(normal, B - A) < 0) CYCLE
		
		! 3
		aa = B - A
		bb = D - A
		
		normal(1) = aa(2) * bb(3) - aa(3) * bb(2) 
		normal(2) = aa(3) * bb(1) - aa(1) * bb(3) 
		normal(3) = aa(1) * bb(2) - aa(2) * bb(1)  
		
		if ( DOT_PRODUCT(normal, r - A) * DOT_PRODUCT(normal, C - A) < 0) CYCLE
		
		! 4
		aa = B - C
		bb = A - C
		
		normal(1) = aa(2) * bb(3) - aa(3) * bb(2) 
		normal(2) = aa(3) * bb(1) - aa(1) * bb(3) 
		normal(3) = aa(1) * bb(2) - aa(2) * bb(1) 
		
		if ( DOT_PRODUCT(normal, r - C) * DOT_PRODUCT(normal, A - C) < 0) CYCLE
		
		EXIT
		
	end do
	
	if (i == 9) then  ! Если не нашли нужный тетраэдр, значит точка лежит вне тетраэдров из-за кривизны ячейки
		! Будем искать близжайщий тетраэдр
		i_best = 1
		dist2 = 10000000000.0_8
		
		do j = 1, 8
			
		dist = 0.0_8
		s1 = gl_Cell_neighbour_inter( par_tet(1, j), num)
		s2 = gl_Cell_neighbour_inter( par_tet(2, j), num)
		s3 = gl_Cell_neighbour_inter( par_tet(3, j), num)
		
		if(s1 <= 0 .or. s2 <= 0 .or. s3 <= 0) CYCLE
		
		B = gl_Cell_center_inter(:, s1)
		C = gl_Cell_center_inter(:, s2)
		D = gl_Cell_center_inter(:, s3)
		
		! 1
		aa = B - C
		bb = D - C
		
		normal(1) = aa(2) * bb(3) - aa(3) * bb(2) 
		normal(2) = aa(3) * bb(1) - aa(1) * bb(3) 
		normal(3) = aa(1) * bb(2) - aa(2) * bb(1) 
		
		if ( DOT_PRODUCT(normal, r - C) * DOT_PRODUCT(normal, A - C) < 0) then
			if ( dabs(DOT_PRODUCT(normal, r - C)) > dist ) then
				dist = dabs(DOT_PRODUCT(normal, r - C))
			end if
		end if
		
		! 2
		aa = C - A
		bb = D - A
		
		normal(1) = aa(2) * bb(3) - aa(3) * bb(2) 
		normal(2) = aa(3) * bb(1) - aa(1) * bb(3) 
		normal(3) = aa(1) * bb(2) - aa(2) * bb(1) 
		
		if ( DOT_PRODUCT(normal, r - A) * DOT_PRODUCT(normal, B - A) < 0) then
			if ( dabs(DOT_PRODUCT(normal, r - A)) > dist ) then
				dist = dabs(DOT_PRODUCT(normal, r - A))
			end if
		end if
		
		! 3
		aa = B - A
		bb = D - A
		
		normal(1) = aa(2) * bb(3) - aa(3) * bb(2) 
		normal(2) = aa(3) * bb(1) - aa(1) * bb(3) 
		normal(3) = aa(1) * bb(2) - aa(2) * bb(1)  
		
		if ( DOT_PRODUCT(normal, r - A) * DOT_PRODUCT(normal, C - A) < 0) then
			if ( dabs(DOT_PRODUCT(normal, r - A)) > dist ) then
				dist = dabs(DOT_PRODUCT(normal, r - A))
			end if
		end if
		
		! 4
		aa = B - C
		bb = A - C
		
		normal(1) = aa(2) * bb(3) - aa(3) * bb(2) 
		normal(2) = aa(3) * bb(1) - aa(1) * bb(3) 
		normal(3) = aa(1) * bb(2) - aa(2) * bb(1) 
		
		if ( DOT_PRODUCT(normal, r - C) * DOT_PRODUCT(normal, A - C) < 0) then
			if ( dabs(DOT_PRODUCT(normal, r - C)) > dist ) then
				dist = dabs(DOT_PRODUCT(normal, r - C))
			end if
		end if
		
		if(dist < dist2) then
			i_best = j
			dist2 = dist
		end if
	end do
		
		i = i_best
		
	end if
	
	!print*, "i = ", i
	!pause
	
	s1 = gl_Cell_neighbour_inter( par_tet(1, i), num)
	s2 = gl_Cell_neighbour_inter( par_tet(2, i), num)
	s3 = gl_Cell_neighbour_inter( par_tet(3, i), num)
	
	B = gl_Cell_center_inter(:, s1)
	C = gl_Cell_center_inter(:, s2)
	D = gl_Cell_center_inter(:, s3)
	
	!print*, A
	!print*, B
	!print*, C
	!print*, D
	!PAUSE
	
	M(:, 1) = 1.0_8
	M(1, 2:4) = A
	M(2, 2:4) = B
	M(3, 2:4) = C
	M(4, 2:4) = D
	
	Minv = M
	n = size(M,1)
	
	!$inter call DGETRF(n, n, Minv, n, ipiv, info)
	
	if (info /= 0) then
     stop 'Matrix is numerically singular! wfergt56y454erthy5t'
	end if
	
	!$inter call DGETRI(n, Minv, n, ipiv, work, n, info)
	
	if (info /= 0) then
     stop 'Matrix inversion failed! frfg543565tttg'
	end if
	
	
	
	vec(1, 1) = 1.0_8
	vec(1, 2:4) = r
	
	!print*, Minv
	
	vec = MATMUL(vec, Minv)
	
	!print*, "___"
	!print*, Vec
	!print*, "___"
	!print*, gl_Cell_par(1, num), gl_Cell_par(1, s1), gl_Cell_par(1, s2), gl_Cell_par(1, s3)
	
	F = vec(1, 1) * gl_Cell_par_inter(:, num) + vec(1, 2) * gl_Cell_par_inter(:, s1) + vec(1, 3) * gl_Cell_par_inter(:, s2) + vec(1, 4) * gl_Cell_par_inter(:, s3)
	F_mf = vec(1, 1) * gl_Cell_par_MF_inter(:, :, num) + vec(1, 2) * gl_Cell_par_MF_inter(:, :, s1) + vec(1, 3) * gl_Cell_par_MF_inter(:, :, s2) + vec(1, 4) * gl_Cell_par_MF_inter(:, :, s3)
	
	if(ieee_is_nan(F(1))) STOP "NUN 678ojhyuikjhikj"
	if(ieee_is_nan(F(2))) STOP "NUN vbnjko9876t"
	if(ieee_is_nan(F(5))) STOP "NUN 13esasdfg"
	if(ieee_is_nan(F_mf(1, 1))) STOP "NUN 9ol,mnbgtyu"
	if(ieee_is_nan(F_mf(1, 2))) STOP "NUN 876rfghjg"
	if(ieee_is_nan(F_mf(2, 1))) STOP "NUN bnhju78ijh"
	if(ieee_is_nan(F_mf(2, 2))) STOP "NUN zxdse45"
	
	!print*, "___"
	!print*, F(1)
	
	end subroutine Find_tetraedr_Interpolate
	
	subroutine Streem_line(x, y, z, num)
	! Variables
	
	! Body of Streem_line
	
	integer, intent(in) :: num
	real(8), intent(in) :: x, y, z
	
	real(8) :: r(3), V(3)
    character(len=5) :: name
	real(8) :: F(9), F_mf(5, 4)
	integer(4) :: zone, cell, number, N, i
	
	r(1) = x
	r(2) = y
	r(3) = z
	number = 1
	cell = 1
	V = 0.0
	N = 0
    
    write(unit=name,fmt='(i5.5)') num
    
    open(1, file = "streem_line_" // name // ".txt")
	
	
11	continue	
	N = N + 1
	r = r + V
	number = cell
	call Interpolate_point(r(1), r(2), r(3), F, F_mf, zone, cell, number)
	V = F(2:4)
	
	if (r(1) > -400.0 .and. sqrt(r(2)**2 + r(3)**2) < 250.0 .and. r(1) < 150.0 ) then
		GO TO 11
	end if
	
	write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z'  ZONE T= 'HP', N= ", N, ", E =  ", N - 1 , ", F=FEPOINT, ET=LINESEG "
	
	r(1) = x
	r(2) = y
	r(3) = z
	number = 1
	cell = 1
	V = 0.0
	
	12	continue	
	r = r + V
	number = cell
	call Interpolate_point(r(1), r(2), r(3), F, F_mf, zone, cell, number)
	V = F(2:4)
	write(1,*) r(1), r(2), r(3)
	
	if (r(1) > -400.0 .and. sqrt(r(2)**2 + r(3)**2) < 250.0 .and. r(1) < 150.0 ) then
		GO TO 12
	end if
	
	
	do i = 1, N - 1
		write(1,*) i, i + 1
	end do
	
	close(1)
	
	end subroutine Streem_line
	
	end module Interpolate
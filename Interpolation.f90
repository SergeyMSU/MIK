! Алгоритм интерполяции основан на исходной сетке
! Сначала нужно найти ячейку, в которой находится нижная интерполяционная точка
! Далее определить по каким соседям нужно построить интерполяционную фигуру (тетраэдр)
! Далее нужно интерполировать значения нужных переменных внутри тетраэдра и вывести ответ
	
! ПРОБЛЕМА
! В интерполяции такого рода есть "дыры". А именно точка может лежать возле пересечения 4-х ячеек и в силу
! погрешности и неточности определения грани по четырём точкам алгоритм будет прыгать покругу 
! (по этим четырём ячейкам) и не найдёт, где именно она лежит
	
	module Interpolate2  
	! Модуль интерполяции по тетраэдрам
    ! Строится двойственная сетка к основной (с дополнительными узлами на разрывах)
	! Строиться сетка тетраэдров на двойственной сетке
	USE GEO_PARAM
	USE STORAGE
	USE My_func
	USE Solvers
	
	
	logical :: int2_work             ! Определена ли интерполяционная сетка (выделена ли память и т.д.)
	
	real(8), allocatable :: int2_coord(:, :)    ! (3, :) набор координат двойственной сетки (центры ячеек сетки 1 + дополнительные)
	
	real(8), allocatable :: int2_Cell_par(:, :)           ! (9, :) Набор параметров (8 стартовых + Q)
	real(8), allocatable :: int2_Cell_par_2(:, :)           ! (1?, :) (n_He)
	integer(4), allocatable :: int2_Cell_par2(:, :)      ! (1, :) Набор параметров (зона)                        INFO
	real(8), allocatable :: int2_Moment(:, :, :)  ! (par_n_moment, par_n_sort, :)
	real(8), allocatable :: int2_Moment_k(:, :)  ! (5 масса - три импульса - энергия, : число точек)
	! Для массы это просто значение, для остальных это коэффициент отношения
	
    !real(8), allocatable :: int2_Cell_par_MF(:,:,:)           ! Набор параметров (5, 4,:)  Мультифлюид параметры (по 5 для каждой из 4-х жидкостей)
	
	
	
    integer(4), allocatable :: int2_all_Cell(:, :)   ! Весь набор ячеек (8, :) - первая координата массива - это набор узлов ячейки
	integer(4), allocatable :: int2_Cell_A(:, :, :)
	integer(4), allocatable :: int2_Cell_B(:, :, :)
	integer(4), allocatable :: int2_Cell_C(:, :, :)
	
	real(8), allocatable :: int2_Cell_center(:, :)   ! Центр ячеек двойственной сетки (3, :) 
	integer(4), allocatable :: int2_all_neighbours(:, :)  ! (6, :)  по 6 соседей на каждую ячейку
	
	integer(4), allocatable :: int2_Point_A(:, :, :)   ! Набор A-точек размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
	! Это точки - центры ячеек основной сетки
	integer(4), allocatable :: int2_Point_B(:, :, :)   ! Набор B-точек размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
	integer(4), allocatable :: int2_Point_C(:, :, :)   ! Набор C-точек размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
	
	! Тетраэдры
	integer(4), allocatable :: int2_all_tetraendron(:, :)  ! (4 грани, :)
	integer(4), allocatable :: int2_all_tetraendron_point(:, :)  ! (4 точки, :)
	real(8), allocatable :: int2_all_tetraendron_matrix(:, :, :)  ! (4, 4, :)
	integer(4), allocatable :: int2_gran_point(:, :)  ! (3 точки, :)
	integer(4), allocatable :: int2_gran_sosed(:)  ! (:)  номер тетраэдра, соседнего с этой гранью
	!real(8), allocatable :: int2_gran_normal(:, :)  ! (3 точки, :)
	real(8), allocatable :: int2_plane_tetraendron(:, :, :)  ! (4 числа задают плоскость, 4 грани, : тетраэдров)
	! A x + B y + C z + D = 0  при этом ABC - внешняя нормаль к тетраэдру (если > 0, то точка вне тетраэдра)
	real(8), allocatable :: int2_all_Volume(:)  ! Объём тетраэдров
	
	
	contains
	
	
	subroutine Int2_Initial()
	
	integer(4) :: i, j, k, N1, N2, N3, ijk, N, kk, m, ii, jj, s1, s2, ss1, ss2, l
	real(8) :: a1(3), a2(3), a3(3), a4(3), b1(3), b2(3), b3(3), aa(3), bb(3), normal(3), S, di
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
				int2_Cell_par(:, gl_Cell_A(i, j, k)) = gl_Cell_par(:, gl_Cell_A(i, j, k) )
				int2_Cell_par_2(:, gl_Cell_A(i, j, k)) = gl_Cell_par2(:, gl_Cell_A(i, j, k) )
				int2_Cell_par2(1, gl_Cell_A(i, j, k)) = gl_zone_Cell(gl_Cell_A(i, j, k))
				
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2 
			do i = par_n_TS, par_n_HP - 1
				int2_Point_A(i + 3, j + 1, k) = gl_Cell_A(i, j, k)
				int2_coord(:, gl_Cell_A(i, j, k)) = gl_Cell_center(:, gl_Cell_A(i, j, k))
				int2_Cell_par(:, gl_Cell_A(i, j, k)) = gl_Cell_par(:, gl_Cell_A(i, j, k) )
				int2_Cell_par_2(:, gl_Cell_A(i, j, k)) = gl_Cell_par2(:, gl_Cell_A(i, j, k) )
				int2_Cell_par2(1, gl_Cell_A(i, j, k)) = gl_zone_Cell(gl_Cell_A(i, j, k))
					
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2
			do i = par_n_HP, par_n_BS - 1
				int2_Point_A(i + 5, j + 1, k) = gl_Cell_A(i, j, k)
				int2_coord(:, gl_Cell_A(i, j, k)) = gl_Cell_center(:, gl_Cell_A(i, j, k))
				int2_Cell_par(:, gl_Cell_A(i, j, k)) = gl_Cell_par(:, gl_Cell_A(i, j, k) )
				int2_Cell_par_2(:, gl_Cell_A(i, j, k)) = gl_Cell_par2(:, gl_Cell_A(i, j, k) )
				int2_Cell_par2(1, gl_Cell_A(i, j, k)) = gl_zone_Cell(gl_Cell_A(i, j, k))
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2
			do i = par_n_BS, par_n_END - 1
				int2_Point_A(i + 5, j + 1, k) = gl_Cell_A(i, j, k)
				int2_coord(:, gl_Cell_A(i, j, k)) = gl_Cell_center(:, gl_Cell_A(i, j, k))
				int2_Cell_par(:, gl_Cell_A(i, j, k)) = gl_Cell_par(:, gl_Cell_A(i, j, k) )
				int2_Cell_par_2(:, gl_Cell_A(i, j, k)) = gl_Cell_par2(:, gl_Cell_A(i, j, k) )
				int2_Cell_par2(1, gl_Cell_A(i, j, k)) = gl_zone_Cell(gl_Cell_A(i, j, k))
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
				int2_Cell_par(:, gl_Cell_B(i, j, k)) = gl_Cell_par(:, gl_Cell_B(i, j, k) )
				int2_Cell_par_2(:, gl_Cell_B(i, j, k)) = gl_Cell_par2(:, gl_Cell_B(i, j, k) )
				int2_Cell_par2(1, gl_Cell_B(i, j, k)) = gl_zone_Cell(gl_Cell_B(i, j, k))
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2
			do i = par_n_TS, N1
				int2_Point_B(i + 3, j + 1, k) = gl_Cell_B(i, j, k)
				int2_coord(:, gl_Cell_B(i, j, k)) = gl_Cell_center(:, gl_Cell_B(i, j, k))
				int2_Cell_par(:, gl_Cell_B(i, j, k)) = gl_Cell_par(:, gl_Cell_B(i, j, k) )
				int2_Cell_par_2(:, gl_Cell_B(i, j, k)) = gl_Cell_par2(:, gl_Cell_B(i, j, k) )
				int2_Cell_par2(1, gl_Cell_B(i, j, k)) = gl_zone_Cell(gl_Cell_B(i, j, k))
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
				int2_Cell_par(:, gl_Cell_C(i, j, k)) = gl_Cell_par(:, gl_Cell_C(i, j, k) )
				int2_Cell_par_2(:, gl_Cell_C(i, j, k)) = gl_Cell_par2(:, gl_Cell_C(i, j, k) )
				int2_Cell_par2(1, gl_Cell_C(i, j, k)) = gl_zone_Cell(gl_Cell_C(i, j, k))
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2
			do i = par_n_HP - par_n_TS + 1, par_n_BS - par_n_TS
				int2_Point_C(i + 2, j, k) = gl_Cell_C(i, j, k)
				int2_coord(:, gl_Cell_C(i, j, k)) = gl_Cell_center(:, gl_Cell_C(i, j, k))
				int2_Cell_par(:, gl_Cell_C(i, j, k)) = gl_Cell_par(:, gl_Cell_C(i, j, k) )
				int2_Cell_par_2(:, gl_Cell_C(i, j, k)) = gl_Cell_par2(:, gl_Cell_C(i, j, k) )
				int2_Cell_par2(1, gl_Cell_C(i, j, k)) = gl_zone_Cell(gl_Cell_C(i, j, k))
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2
			do i = par_n_BS - par_n_TS + 1, par_n_END - par_n_TS
				int2_Point_C(i + 2, j, k) = gl_Cell_C(i, j, k)
				int2_coord(:, gl_Cell_C(i, j, k)) = gl_Cell_center(:, gl_Cell_C(i, j, k))
				int2_Cell_par(:, gl_Cell_C(i, j, k)) = gl_Cell_par(:, gl_Cell_C(i, j, k) )
				int2_Cell_par_2(:, gl_Cell_C(i, j, k)) = gl_Cell_par2(:, gl_Cell_C(i, j, k) )
				int2_Cell_par2(1, gl_Cell_C(i, j, k)) = gl_zone_Cell(gl_Cell_C(i, j, k))
			end do
		end do
	end do
	
	
	! Заполняем особые точки (поверхности разрывов)
	N = size(gl_Cell_C) + size(gl_Cell_A) + size(gl_Cell_B)
	int2_Point_A(1, :, :) = N + 1
	int2_Point_B(1, :, :) = N + 1
	int2_coord(:, N + 1) = (/ 0.0, 0.0, 0.0 /)
	int2_Cell_par(:, N + 1) = gl_Cell_par(:, gl_Cell_A(1, 1, 1) )
	int2_Cell_par_2(:, N + 1) = gl_Cell_par2(:, gl_Cell_A(1, 1, 1) )
	int2_Cell_par2(:, N + 1) = 1
	N = N + 1
	
	N1 = size(gl_Cell_A(:, 1, 1))
	N2 = size(gl_Cell_A(1, :, 1))
	N3 = size(gl_Cell_A(1, 1, :))
	
	! A
	do k = 1, N3
		do j = 1, N2
			int2_Point_A(par_n_TS + 1, j + 1, k) = N + 1
			int2_coord(:, N + 1) = gl_Gran_center(:, gl_Cell_gran(1, gl_Cell_A(par_n_TS - 1, j, k) ))
			int2_Cell_par(:, N + 1) = gl_Cell_par(:, gl_Cell_A(par_n_TS - 1, j, k) )
			int2_Cell_par_2(:, N + 1) = gl_Cell_par2(:, gl_Cell_A(par_n_TS - 1, j, k) )
			int2_Cell_par2(1, N + 1) = gl_zone_Cell(gl_Cell_A(par_n_TS - 1, j, k) )
			N = N + 1
			
			int2_Point_A(par_n_TS + 2, j + 1, k) = N + 1
			int2_coord(:, N + 1) = gl_Gran_center(:, gl_Cell_gran(1, gl_Cell_A(par_n_TS - 1, j, k) ))
			int2_Cell_par(:, N + 1) = gl_Cell_par(:, gl_Cell_A(par_n_TS, j, k) )
			int2_Cell_par_2(:, N + 1) = gl_Cell_par2(:, gl_Cell_A(par_n_TS, j, k) )
			int2_Cell_par2(1, N + 1) = gl_zone_Cell(gl_Cell_A(par_n_TS, j, k) )
			N = N + 1
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2
			int2_Point_A(par_n_HP + 3, j + 1, k) = N + 1
			int2_coord(:, N + 1) = gl_Gran_center(:, gl_Cell_gran(1, gl_Cell_A(par_n_HP - 1, j, k) ))
			int2_Cell_par(:, N + 1) = gl_Cell_par(:, gl_Cell_A(par_n_HP - 1, j, k) )
			int2_Cell_par_2(:, N + 1) = gl_Cell_par2(:, gl_Cell_A(par_n_HP - 1, j, k) )
			int2_Cell_par2(1, N + 1) = gl_zone_Cell(gl_Cell_A(par_n_HP - 1, j, k) )
			N = N + 1
			
			int2_Point_A(par_n_HP + 4, j + 1, k) = N + 1
			int2_coord(:, N + 1) = gl_Gran_center(:, gl_Cell_gran(1, gl_Cell_A(par_n_HP - 1, j, k) ))
			int2_Cell_par(:, N + 1) = gl_Cell_par(:, gl_Cell_A(par_n_HP, j, k) )
			int2_Cell_par_2(:, N + 1) = gl_Cell_par2(:, gl_Cell_A(par_n_HP, j, k) )
			int2_Cell_par2(1, N + 1) = gl_zone_Cell(gl_Cell_A(par_n_HP, j, k) )
			N = N + 1
		end do
	end do
	
	! Добавим ряд точек на границе полусферы
	
	do k = 1, N3
		do j = 1, N2
			do i = par_n_END, par_n_END
				int2_Point_A(i + 5, j + 1, k) = N + 1
				aa = gl_Cell_center(:, gl_Cell_A(i - 1, j, k))
				
				if(aa(1) >= 0.0) then
					di = norm2(aa)
					int2_coord(:, N + 1) = aa/di * par_R_END
				else
					di = sqrt(aa(2)**2 + aa(3)**2)
					aa(2) = aa(2)/di * par_R_END
					aa(3) = aa(3)/di * par_R_END
					int2_coord(:, N + 1) = aa
				end if
				
				
				int2_Cell_par(:, N + 1) = gl_Cell_par(:, gl_Cell_A(i - 1, j, k) )
				int2_Cell_par_2(:, N + 1) = gl_Cell_par2(:, gl_Cell_A(i - 1, j, k) )
				int2_Cell_par2(1, N + 1) = gl_zone_Cell(gl_Cell_A(i - 1, j, k))
				N = N + 1
			end do
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
			int2_Cell_par(:, N + 1) = gl_Cell_par(:, gl_Cell_B(par_n_TS - 1, j, k) )
			int2_Cell_par_2(:, N + 1) = gl_Cell_par2(:, gl_Cell_B(par_n_TS - 1, j, k) )
			int2_Cell_par2(1, N + 1) = gl_zone_Cell(gl_Cell_B(par_n_TS - 1, j, k) )
			N = N + 1
			
			int2_Point_B(par_n_TS + 2, j + 1, k) = N + 1
			int2_coord(:, N + 1) = gl_Gran_center(:, gl_Cell_gran(1, gl_Cell_B(par_n_TS - 1, j, k) ))
			int2_Cell_par(:, N + 1) = gl_Cell_par(:, gl_Cell_B(par_n_TS, j, k) )
			int2_Cell_par_2(:, N + 1) = gl_Cell_par2(:, gl_Cell_B(par_n_TS, j, k) )
			int2_Cell_par2(1, N + 1) = gl_zone_Cell(gl_Cell_B(par_n_TS, j, k) )
			N = N + 1 
		end do
	end do
	
	! Добавим ряд точек на границе 
	
	do k = 1, N3
		do j = 1, N2
			do i = N1 + 1, N1 + 1
				int2_Point_B(i + 3, j + 1, k) = N + 1
				aa = gl_Cell_center(:, gl_Cell_B(i - 1, j, k))
				int2_coord(:, N + 1) = (/par_R_LEFT, aa(2), aa(3)/)
				int2_Cell_par(:, N + 1) = gl_Cell_par(:, gl_Cell_B(i - 1, j, k) )
				int2_Cell_par_2(:, N + 1) = gl_Cell_par2(:, gl_Cell_B(i - 1, j, k) )
				int2_Cell_par2(1, N + 1) = gl_zone_Cell(gl_Cell_B(i - 1, j, k))
				N = N + 1 
			end do
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
			int2_Cell_par(:, N + 1) = gl_Cell_par(:, gl_Cell_C(par_n_HP - par_n_TS, j, k) )
			int2_Cell_par_2(:, N + 1) = gl_Cell_par2(:, gl_Cell_C(par_n_HP - par_n_TS, j, k) )
			int2_Cell_par2(1, N + 1) = gl_zone_Cell(gl_Cell_C(par_n_HP - par_n_TS, j, k) )
			N = N + 1
			
			int2_Point_C(par_n_HP - par_n_TS + 2, j, k) = N + 1
			int2_coord(:, N + 1) = gl_Gran_center(:, gl_Cell_gran(1, gl_Cell_C(par_n_HP - par_n_TS, j, k) ))
			int2_Cell_par(:, N + 1) = gl_Cell_par(:, gl_Cell_C(par_n_HP - par_n_TS + 1, j, k) )
			int2_Cell_par_2(:, N + 1) = gl_Cell_par2(:, gl_Cell_C(par_n_HP - par_n_TS + 1, j, k) )
			int2_Cell_par2(1, N + 1) = gl_zone_Cell(gl_Cell_C(par_n_HP - par_n_TS + 1, j, k) )
			N = N + 1 
		end do
	end do
	
	! Добавим два ряда на границе
	do k = 1, N3
		do j = 1, N2
			do i = par_n_END - par_n_TS + 1, par_n_END - par_n_TS + 1
				int2_Point_C(i + 2, j, k) = N + 1
				aa = gl_Cell_center(:, gl_Cell_C(i - 1, j, k))
				di = sqrt(aa(2)**2 + aa(3)**2)
				aa(2) = aa(2)/di * par_R_END
				aa(3) = aa(3)/di * par_R_END
				int2_coord(:, N + 1) = aa
				int2_Cell_par(:, N + 1) = gl_Cell_par(:, gl_Cell_C(i - 1, j, k) )
				int2_Cell_par_2(:, N + 1) = gl_Cell_par2(:, gl_Cell_C(i - 1, j, k) )
				int2_Cell_par2(1, N + 1) = gl_zone_Cell(gl_Cell_C(i - 1, j, k))
				N = N + 1 
			end do
		end do
	end do
	
	N1 = size(int2_Point_C(:, 1, 1))
	N2 = size(int2_Point_C(1, :, 1))
	N3 = size(int2_Point_C(1, 1, :))
	
	do k = 1, N3
		do j = N2, N2
			do i = 1, N1
				int2_Point_C(i, j, k) = N + 1
				aa = int2_coord(:, int2_Point_C(i, j - 1, k))
				aa(1) = par_R_LEFT
				int2_coord(:, N + 1) = aa
				int2_Cell_par(:, N + 1) = int2_Cell_par(:, int2_Point_C(i, j - 1, k))
				int2_Cell_par_2(:, N + 1) = int2_Cell_par_2(:, int2_Point_C(i, j - 1, k))
				int2_Cell_par2(1, N + 1) = int2_Cell_par2(1, int2_Point_C(i, j - 1, k))
				N = N + 1 
			end do
		end do
	end do
	
	! Ячейки на оси
	N1 = size(int2_Point_A(:, 1, 1))
	N2 = size(int2_Point_A(1, :, 1))
	N3 = size(int2_Point_A(1, 1, :))
	
	!A
	do i = 2, N1
		int2_Point_A(i, 1, :) = N + 1
		int2_coord(:, N + 1) = 0.0
		int2_Cell_par(:, N + 1) = 0.0
		int2_Cell_par_2(:, N + 1) = 0.0
		int2_Cell_par2(1, N + 1) = int2_Cell_par2(1, int2_Point_A(i, 2, 1))
		
		do k = 1, N3
			int2_coord(:, N + 1) = int2_coord(:, N + 1) + (/ int2_coord(1, int2_Point_A(i, 2, k)), 0.0_8, 0.0_8 /)
			int2_Cell_par(:, N + 1) = int2_Cell_par(:, N + 1) + int2_Cell_par(:, int2_Point_A(i, 2, k))
			int2_Cell_par_2(:, N + 1) = int2_Cell_par_2(:, N + 1) + int2_Cell_par_2(:, int2_Point_A(i, 2, k))
		end do
		int2_coord(:, N + 1) = int2_coord(:, N + 1)/N3
		int2_Cell_par(:, N + 1) = int2_Cell_par(:, N + 1)/N3
		int2_Cell_par_2(:, N + 1) = int2_Cell_par_2(:, N + 1)/N3
		N = N + 1
	end do
	
	N1 = size(int2_Point_B(:, 1, 1))
	N2 = size(int2_Point_B(1, :, 1))
	N3 = size(int2_Point_B(1, 1, :))
	
	!B
	do i = 2, N1
		int2_Point_B(i, 1, :) = N + 1
		int2_coord(:, N + 1) = 0.0
		int2_Cell_par(:, N + 1) = 0.0
		int2_Cell_par_2(:, N + 1) = 0.0
		int2_Cell_par2(1, N + 1) = int2_Cell_par2(1, int2_Point_B(i, 2, 1))
		
		do k = 1, N3
			int2_coord(:, N + 1) = int2_coord(:, N + 1) + (/ int2_coord(1, int2_Point_B(i, 2, k)), 0.0_8, 0.0_8 /)
			int2_Cell_par(:, N + 1) = int2_Cell_par(:, N + 1) + int2_Cell_par(:, int2_Point_B(i, 2, k))
			int2_Cell_par_2(:, N + 1) = int2_Cell_par_2(:, N + 1) + int2_Cell_par_2(:, int2_Point_B(i, 2, k))
		end do
		int2_coord(:, N + 1) = int2_coord(:, N + 1)/N3
		int2_Cell_par(:, N + 1) = int2_Cell_par(:, N + 1)/N3
		int2_Cell_par_2(:, N + 1) = int2_Cell_par_2(:, N + 1)/N3
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
					ijk = int2_all_Cell(m, N)
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
				int2_all_Cell(2, N) = int2_Point_B(i + 2, size(int2_Point_B(1, :, 1)), k)
				int2_all_Cell(3, N) = int2_Point_B(i + 2, size(int2_Point_B(1, :, 1)), k)
				int2_all_Cell(4, N) = int2_Point_B(i + 1, size(int2_Point_B(1, :, 1)), k)
				int2_all_Cell(5, N) = int2_Point_A(i, j, kk)
				int2_all_Cell(6, N) = int2_Point_B(i + 2, size(int2_Point_B(1, :, 1)), kk)
				int2_all_Cell(7, N) = int2_Point_B(i + 2, size(int2_Point_B(1, :, 1)), kk)
				int2_all_Cell(8, N) = int2_Point_B(i + 1, size(int2_Point_B(1, :, 1)), kk)
				
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
				if(int2_Cell_A(i, j, k) == 0) CYCLE
				int2_all_neighbours(1, int2_Cell_A(i, j, k)) = int2_Cell_A(i + 1, j, k)
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2
			do i = 2, N1
				if(int2_Cell_A(i, j, k) == 0) CYCLE
				int2_all_neighbours(2, int2_Cell_A(i, j, k)) = int2_Cell_A(i - 1, j, k)
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2 - 1
			do i = 1, N1
				if(int2_Cell_A(i, j, k) == 0) CYCLE
				int2_all_neighbours(3, int2_Cell_A(i, j, k)) = int2_Cell_A(i, j + 1, k)
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 2, N2
			do i = 1, N1
				if(int2_Cell_A(i, j, k) == 0) CYCLE
				int2_all_neighbours(4, int2_Cell_A(i, j, k)) = int2_Cell_A(i, j - 1, k)
			end do
		end do
	end do
	
	do k = 1, N3
		kk = k + 1
		if (kk > N3) kk = 1
		
		do j = 1, N2
			do i = 1, N1
				if(int2_Cell_A(i, j, k) == 0) CYCLE
				int2_all_neighbours(5, int2_Cell_A(i, j, k)) = int2_Cell_A(i, j, kk)
			end do
		end do
	end do
	
	do k = 1, N3
		kk = k - 1
		if (kk < 1) kk = N3
		
		do j = 1, N2
			do i = 1, N1
				if(int2_Cell_A(i, j, k) == 0) CYCLE
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
				if(int2_Cell_B(i, j, k) == 0) CYCLE
				int2_all_neighbours(1, int2_Cell_B(i, j, k)) = int2_Cell_B(i + 1, j, k)
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2
			do i = 2, N1
				if(int2_Cell_B(i, j, k) == 0) CYCLE
				int2_all_neighbours(2, int2_Cell_B(i, j, k)) = int2_Cell_B(i - 1, j, k)
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2 - 1
			do i = 1, N1
				if(int2_Cell_B(i, j, k) == 0) CYCLE
				int2_all_neighbours(3, int2_Cell_B(i, j, k)) = int2_Cell_B(i, j + 1, k)
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 2, N2
			do i = 1, N1
				if(int2_Cell_B(i, j, k) == 0) CYCLE
				int2_all_neighbours(4, int2_Cell_B(i, j, k)) = int2_Cell_B(i, j - 1, k)
			end do
		end do
	end do
	
	do k = 1, N3
		kk = k + 1
		if (kk > N3) kk = 1
		
		do j = 1, N2
			do i = 1, N1
				if(int2_Cell_B(i, j, k) == 0) CYCLE
				int2_all_neighbours(5, int2_Cell_B(i, j, k)) = int2_Cell_B(i, j, kk)
			end do
		end do
	end do
	
	do k = 1, N3
		kk = k - 1
		if (kk < 1) kk = N3
		
		do j = 1, N2
			do i = 1, N1
				if(int2_Cell_B(i, j, k) == 0) CYCLE
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
				if(int2_Cell_C(i, j, k) == 0) CYCLE
				int2_all_neighbours(1, int2_Cell_C(i, j, k)) = int2_Cell_C(i + 1, j, k)
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2
			do i = 2, N1
				if(int2_Cell_C(i, j, k) == 0) CYCLE
				int2_all_neighbours(2, int2_Cell_C(i, j, k)) = int2_Cell_C(i - 1, j, k)
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2 - 1
			do i = 1, N1
				if(int2_Cell_C(i, j, k) == 0) CYCLE
				int2_all_neighbours(3, int2_Cell_C(i, j, k)) = int2_Cell_C(i, j + 1, k)
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 2, N2
			do i = 1, N1
				if(int2_Cell_C(i, j, k) == 0) CYCLE
				int2_all_neighbours(4, int2_Cell_C(i, j, k)) = int2_Cell_C(i, j - 1, k)
			end do
		end do
	end do
	
	do k = 1, N3
		kk = k + 1
		if (kk > N3) kk = 1
		
		do j = 1, N2
			do i = 1, N1
				if(int2_Cell_C(i, j, k) == 0) CYCLE
				int2_all_neighbours(5, int2_Cell_C(i, j, k)) = int2_Cell_C(i, j, kk)
			end do
		end do
	end do
	
	do k = 1, N3
		kk = k - 1
		if (kk < 1) kk = N3
		
		do j = 1, N2
			do i = 1, N1
				if(int2_Cell_C(i, j, k) == 0) CYCLE
				int2_all_neighbours(6, int2_Cell_C(i, j, k)) = int2_Cell_C(i, j, kk)
			end do
		end do
	end do
	
	! Теперь нужно отдельно прописать связи на стыках групп + на поверхностях разрыва + в тройной точке
	
	! TS в B - группе
	N1 = size(int2_Cell_B(:, 1, 1))
	N2 = size(int2_Cell_B(1, :, 1))
	N3 = size(int2_Cell_B(1, 1, :))
	
	do k = 1, N3
		do j = 1, N2
			do i = par_n_TS, par_n_TS
				int2_all_neighbours(1, int2_Cell_B(i, j, k)) = int2_Cell_B(i + 2, j, k)
				int2_all_neighbours(2, int2_Cell_B(i + 2, j, k)) = int2_Cell_B(i, j, k)
			end do
		end do
	end do
	
	
	! TS в A - группе
	N1 = size(int2_Cell_A(:, 1, 1))
	N2 = size(int2_Cell_A(1, :, 1))
	N3 = size(int2_Cell_A(1, 1, :))
	
	do k = 1, N3
		do j = 1, N2
			do i = par_n_TS, par_n_TS
				int2_all_neighbours(1, int2_Cell_A(i, j, k)) = int2_Cell_A(i + 2, j, k)
				int2_all_neighbours(2, int2_Cell_A(i + 2, j, k)) = int2_Cell_A(i, j, k)
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2
			do i = par_n_HP + 2, par_n_HP + 2
				int2_all_neighbours(1, int2_Cell_A(i, j, k)) = int2_Cell_A(i + 2, j, k)
				int2_all_neighbours(2, int2_Cell_A(i + 2, j, k)) = int2_Cell_A(i, j, k)
			end do
		end do
		
		int2_all_neighbours(1, int2_Cell_A(par_n_TS, N2, k)) = 0
		int2_all_neighbours(1, int2_Cell_A(par_n_TS + 1, N2, k)) = 0
		int2_all_neighbours(2, int2_Cell_A(par_n_TS + 1, N2, k)) = 0
		int2_all_neighbours(2, int2_Cell_A(par_n_TS + 2, N2, k)) = 0
		
		int2_all_neighbours(2, int2_Cell_A(par_n_TS + 2, N2, k)) = int2_Cell_A(par_n_TS + 1, N2, k)
		int2_all_neighbours(1, int2_Cell_A(par_n_TS + 1, N2, k)) = int2_Cell_A(par_n_TS + 2, N2, k)
		
		int2_all_neighbours(2, int2_Cell_A(par_n_TS + 1, N2, k)) = int2_Cell_A(par_n_TS, N2, k)
		int2_all_neighbours(1, int2_Cell_A(par_n_TS, N2, k)) = int2_Cell_A(par_n_TS + 1, N2, k)
			
	end do
	
	
	! HP в C - группе
	N1 = size(int2_Cell_C(:, 1, 1))
	N2 = size(int2_Cell_C(1, :, 1))
	N3 = size(int2_Cell_C(1, 1, :))
	
	do k = 1, N3
		do j = 1, N2
			do i = par_n_HP - par_n_TS + 1, par_n_HP - par_n_TS + 1
				int2_all_neighbours(1, int2_Cell_C(i, j, k)) = int2_Cell_C(i + 2, j, k)
				int2_all_neighbours(2, int2_Cell_C(i + 2, j, k)) = int2_Cell_C(i, j, k)
			end do
		end do
	end do
	
	! AB - связи
	N1 = size(int2_Cell_A(:, 1, 1))
	N2 = size(int2_Cell_A(1, :, 1))
	N3 = size(int2_Cell_A(1, 1, :))
	
	do k = 1, N3
		do j = N2, N2
			do i = 1, par_n_TS
				int2_all_neighbours(3, int2_Cell_A(i, j, k)) = int2_Cell_B(i, size(int2_Cell_B(1, :, 1)), k)
				ijk = int2_Cell_B(i, size(int2_Cell_B(1, :, 1)), k)
				int2_all_neighbours(3, int2_Cell_B(i, size(int2_Cell_B(1, :, 1)), k)) = int2_Cell_A(i, j, k)
			end do
			
			int2_all_neighbours(3, int2_Cell_A(par_n_TS + 1, j, k)) = int2_Cell_B(par_n_TS + 2, size(int2_Cell_B(1, :, 1)), k)
			int2_all_neighbours(3, int2_Cell_B(par_n_TS + 2, size(int2_Cell_B(1, :, 1)), k)) = int2_Cell_A(par_n_TS + 1, j, k)
				
		end do
	end do
	
	! AC - связи
	N1 = size(int2_Cell_A(:, 1, 1))
	N2 = size(int2_Cell_A(1, :, 1))
	N3 = size(int2_Cell_A(1, 1, :))
	
	do k = 1, N3
		do j = N2, N2
			do i = par_n_TS + 2, N1
				if(i == par_n_HP + 3) CYCLE
				int2_all_neighbours(3, int2_Cell_A(i, j, k)) = int2_Cell_C(i - par_n_TS	- 1, 1, k)
				ijk = int2_Cell_C(i - par_n_TS - 1, 1, k)
				int2_all_neighbours(4, int2_Cell_C(i - par_n_TS - 1, 1, k)) = int2_Cell_A(i, j, k)
			end do
		end do
	end do
	
	
	! BC - связи
	N1 = size(int2_Cell_C(:, 1, 1))
	N2 = size(int2_Cell_C(1, :, 1))
	N3 = size(int2_Cell_C(1, 1, :))
	
	do k = 1, N3
		do j = 1, N2
			do i = 1, 1
				int2_all_neighbours(2, int2_Cell_C(i, j, k)) = int2_Cell_B(j + par_n_TS + 2, size(int2_Cell_B(1, :, 1)), k)
				ijk = int2_Cell_B(j + par_n_TS + 2, size(int2_Cell_B(1, :, 1)), k)
				int2_all_neighbours(3, int2_Cell_B(j + par_n_TS + 2, size(int2_Cell_B(1, :, 1)), k)) = int2_Cell_C(i, j, k)
			end do
		end do
	end do
	
	! Делаем проверки связей
	
	print*, "Proverka 1"
	
	N1 = size(int2_all_Cell(1, :))
	do i = 1, N1
		loop1: do j = 1, 6
			s1 = int2_all_neighbours(j, i)
			if (s1 < 1) CYCLE
			do k = 1, 6
				if(int2_all_neighbours(k, s1) == i) CYCLE loop1
			end do
			print*, "No sosed ", i, j
		end do loop1
	end do
	
	print*, "End proverka 1"
	
	! Создаём тетраэдры
	N1 = size(int2_Cell_A(:, 1, 1))
	N2 = size(int2_Cell_A(1, :, 1))
	N3 = size(int2_Cell_A(1, 1, :))
	
	do k = 1, N3
		do j = 1, N2
			do i = 1, N1
				s1 = int2_Cell_A(i, j, k)
				if (s1 == 0) CYCLE
				if(j == N2 .and. i == par_n_TS + 1) CYCLE
				
				int2_gran_point(:, (s1 - 1) * 24 + 1) = (/ int2_all_Cell(1, s1), int2_all_Cell(2, s1), int2_all_Cell(8, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 2) = (/ int2_all_Cell(1, s1), int2_all_Cell(2, s1), int2_all_Cell(5, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 3) = (/ int2_all_Cell(1, s1), int2_all_Cell(8, s1), int2_all_Cell(5, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 4) = (/ int2_all_Cell(2, s1), int2_all_Cell(8, s1), int2_all_Cell(5, s1) /)
				int2_all_tetraendron(:, (s1 - 1) * 6 + 1) = (/ (s1 - 1) * 24 + 1, (s1 - 1) * 24 + 2, (s1 - 1) * 24 + 3, (s1 - 1) * 24 + 4 /)
				
				int2_gran_point(:, (s1 - 1) * 24 + 5) = (/ int2_all_Cell(2, s1), int2_all_Cell(7, s1), int2_all_Cell(8, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 6) = (/ int2_all_Cell(2, s1), int2_all_Cell(7, s1), int2_all_Cell(6, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 7) = (/ int2_all_Cell(2, s1), int2_all_Cell(8, s1), int2_all_Cell(6, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 8) = (/ int2_all_Cell(7, s1), int2_all_Cell(8, s1), int2_all_Cell(6, s1) /)
				int2_all_tetraendron(:, (s1 - 1) * 6 + 2) = (/ (s1 - 1) * 24 + 5, (s1 - 1) * 24 + 6, (s1 - 1) * 24 + 7, (s1 - 1) * 24 + 8 /)
				
				int2_gran_point(:, (s1 - 1) * 24 + 9) = (/ int2_all_Cell(2, s1), int2_all_Cell(3, s1), int2_all_Cell(7, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 10) = (/ int2_all_Cell(2, s1), int2_all_Cell(3, s1), int2_all_Cell(8, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 11) = (/ int2_all_Cell(2, s1), int2_all_Cell(7, s1), int2_all_Cell(8, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 12) = (/ int2_all_Cell(3, s1), int2_all_Cell(7, s1), int2_all_Cell(8, s1) /)
				int2_all_tetraendron(:, (s1 - 1) * 6 + 3) = (/ (s1 - 1) * 24 + 9, (s1 - 1) * 24 + 10, (s1 - 1) * 24 + 11, (s1 - 1) * 24 + 12 /)
			
				int2_gran_point(:, (s1 - 1) * 24 + 13) = (/ int2_all_Cell(1, s1), int2_all_Cell(2, s1), int2_all_Cell(4, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 14) = (/ int2_all_Cell(1, s1), int2_all_Cell(2, s1), int2_all_Cell(8, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 15) = (/ int2_all_Cell(1, s1), int2_all_Cell(4, s1), int2_all_Cell(8, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 16) = (/ int2_all_Cell(2, s1), int2_all_Cell(4, s1), int2_all_Cell(8, s1) /)
				int2_all_tetraendron(:, (s1 - 1) * 6 + 4) = (/ (s1 - 1) * 24 + 13, (s1 - 1) * 24 + 14, (s1 - 1) * 24 + 15, (s1 - 1) * 24 + 16 /)
				
				int2_gran_point(:, (s1 - 1) * 24 + 17) = (/ int2_all_Cell(2, s1), int2_all_Cell(3, s1), int2_all_Cell(4, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 18) = (/ int2_all_Cell(2, s1), int2_all_Cell(3, s1), int2_all_Cell(8, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 19) = (/ int2_all_Cell(2, s1), int2_all_Cell(4, s1), int2_all_Cell(8, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 20) = (/ int2_all_Cell(3, s1), int2_all_Cell(4, s1), int2_all_Cell(8, s1) /)
				int2_all_tetraendron(:, (s1 - 1) * 6 + 5) = (/ (s1 - 1) * 24 + 17, (s1 - 1) * 24 + 18, (s1 - 1) * 24 + 19, (s1 - 1) * 24 + 20 /)
				
				int2_gran_point(:, (s1 - 1) * 24 + 21) = (/ int2_all_Cell(2, s1), int2_all_Cell(5, s1), int2_all_Cell(6, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 22) = (/ int2_all_Cell(2, s1), int2_all_Cell(5, s1), int2_all_Cell(8, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 23) = (/ int2_all_Cell(2, s1), int2_all_Cell(6, s1), int2_all_Cell(8, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 24) = (/ int2_all_Cell(5, s1), int2_all_Cell(6, s1), int2_all_Cell(8, s1) /)
				int2_all_tetraendron(:, (s1 - 1) * 6 + 6) = (/ (s1 - 1) * 24 + 21, (s1 - 1) * 24 + 22, (s1 - 1) * 24 + 23, (s1 - 1) * 24 + 24 /)
			
				int2_gran_sosed((s1 - 1) * 24 + 1) = (s1 - 1) * 6 + 4
				int2_gran_sosed((s1 - 1) * 24 + 4) = (s1 - 1) * 6 + 6
				
				int2_gran_sosed((s1 - 1) * 24 + 5) = (s1 - 1) * 6 + 3
				int2_gran_sosed((s1 - 1) * 24 + 7) = (s1 - 1) * 6 + 6
				
				int2_gran_sosed((s1 - 1) * 24 + 10) = (s1 - 1) * 6 + 5
				int2_gran_sosed((s1 - 1) * 24 + 11) = (s1 - 1) * 6 + 2
				
				int2_gran_sosed((s1 - 1) * 24 + 14) = (s1 - 1) * 6 + 1
				int2_gran_sosed((s1 - 1) * 24 + 16) = (s1 - 1) * 6 + 5
				
				int2_gran_sosed((s1 - 1) * 24 + 18) = (s1 - 1) * 6 + 3
				int2_gran_sosed((s1 - 1) * 24 + 19) = (s1 - 1) * 6 + 4
				
				int2_gran_sosed((s1 - 1) * 24 + 22) = (s1 - 1) * 6 + 1
				int2_gran_sosed((s1 - 1) * 24 + 23) = (s1 - 1) * 6 + 2
				
			end do
		end do
	end do
	
	! Отдельно для тройной точки
	do k = 1, N3
		do j = N2, N2
			do i = par_n_TS + 1, par_n_TS + 1
				s1 = int2_Cell_A(i, j, k)
				if (s1 == 0) CYCLE
				int2_gran_point(:, (s1 - 1) * 24 + 1) = (/ int2_all_Cell(1, s1), int2_all_Cell(5, s1), int2_all_Cell(6, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 2) = (/ int2_all_Cell(1, s1), int2_all_Cell(5, s1), int2_all_Cell(8, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 3) = (/ int2_all_Cell(1, s1), int2_all_Cell(6, s1), int2_all_Cell(8, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 4) = (/ int2_all_Cell(5, s1), int2_all_Cell(6, s1), int2_all_Cell(8, s1) /)
				int2_all_tetraendron(:, (s1 - 1) * 6 + 1) = (/ (s1 - 1) * 24 + 1, (s1 - 1) * 24 + 2, (s1 - 1) * 24 + 3, (s1 - 1) * 24 + 4 /)
				
				int2_gran_point(:, (s1 - 1) * 24 + 5) = (/ int2_all_Cell(1, s1), int2_all_Cell(2, s1), int2_all_Cell(6, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 6) = (/ int2_all_Cell(1, s1), int2_all_Cell(2, s1), int2_all_Cell(8, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 7) = (/ int2_all_Cell(1, s1), int2_all_Cell(6, s1), int2_all_Cell(8, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 8) = (/ int2_all_Cell(2, s1), int2_all_Cell(6, s1), int2_all_Cell(8, s1) /)
				int2_all_tetraendron(:, (s1 - 1) * 6 + 2) = (/ (s1 - 1) * 24 + 5, (s1 - 1) * 24 + 6, (s1 - 1) * 24 + 7, (s1 - 1) * 24 + 8 /)
				
				int2_gran_point(:, (s1 - 1) * 24 + 9) = (/ int2_all_Cell(1, s1), int2_all_Cell(2, s1), int2_all_Cell(4, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 10) = (/ int2_all_Cell(1, s1), int2_all_Cell(2, s1), int2_all_Cell(8, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 11) = (/ int2_all_Cell(1, s1), int2_all_Cell(4, s1), int2_all_Cell(8, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 12) = (/ int2_all_Cell(2, s1), int2_all_Cell(4, s1), int2_all_Cell(8, s1) /)
				int2_all_tetraendron(:, (s1 - 1) * 6 + 3) = (/ (s1 - 1) * 24 + 9, (s1 - 1) * 24 + 10, (s1 - 1) * 24 + 11, (s1 - 1) * 24 + 12 /)
			
				
				int2_gran_sosed((s1 - 1) * 24 + 3) = (s1 - 1) * 6 + 2
				int2_gran_sosed((s1 - 1) * 24 + 7) = (s1 - 1) * 6 + 1
				
				int2_gran_sosed((s1 - 1) * 24 + 6) = (s1 - 1) * 6 + 3
				int2_gran_sosed((s1 - 1) * 24 + 10) = (s1 - 1) * 6 + 2
				
				
			end do
		end do
	end do
	
	N1 = size(int2_Cell_B(:, 1, 1))
	N2 = size(int2_Cell_B(1, :, 1))
	N3 = size(int2_Cell_B(1, 1, :))
	
	do k = 1, N3
		do j = 1, N2
			do i = 1, N1
				s1 = int2_Cell_B(i, j, k)
				if (s1 == 0) CYCLE
				int2_gran_point(:, (s1 - 1) * 24 + 1) = (/ int2_all_Cell(2, s1), int2_all_Cell(3, s1), int2_all_Cell(5, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 2) = (/ int2_all_Cell(2, s1), int2_all_Cell(3, s1), int2_all_Cell(6, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 3) = (/ int2_all_Cell(2, s1), int2_all_Cell(5, s1), int2_all_Cell(6, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 4) = (/ int2_all_Cell(3, s1), int2_all_Cell(5, s1), int2_all_Cell(6, s1) /)
				int2_all_tetraendron(:, (s1 - 1) * 6 + 1) = (/ (s1 - 1) * 24 + 1, (s1 - 1) * 24 + 2, (s1 - 1) * 24 + 3, (s1 - 1) * 24 + 4 /)
				
				int2_gran_point(:, (s1 - 1) * 24 + 5) = (/ int2_all_Cell(3, s1), int2_all_Cell(8, s1), int2_all_Cell(5, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 6) = (/ int2_all_Cell(3, s1), int2_all_Cell(8, s1), int2_all_Cell(7, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 7) = (/ int2_all_Cell(3, s1), int2_all_Cell(5, s1), int2_all_Cell(7, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 8) = (/ int2_all_Cell(8, s1), int2_all_Cell(5, s1), int2_all_Cell(7, s1) /)
				int2_all_tetraendron(:, (s1 - 1) * 6 + 2) = (/ (s1 - 1) * 24 + 5, (s1 - 1) * 24 + 6, (s1 - 1) * 24 + 7, (s1 - 1) * 24 + 8 /)
				
				int2_gran_point(:, (s1 - 1) * 24 + 9) = (/ int2_all_Cell(3, s1), int2_all_Cell(4, s1), int2_all_Cell(8, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 10) = (/ int2_all_Cell(3, s1), int2_all_Cell(4, s1), int2_all_Cell(5, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 11) = (/ int2_all_Cell(3, s1), int2_all_Cell(8, s1), int2_all_Cell(5, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 12) = (/ int2_all_Cell(4, s1), int2_all_Cell(8, s1), int2_all_Cell(5, s1) /)
				int2_all_tetraendron(:, (s1 - 1) * 6 + 3) = (/ (s1 - 1) * 24 + 9, (s1 - 1) * 24 + 10, (s1 - 1) * 24 + 11, (s1 - 1) * 24 + 12 /)
			
				int2_gran_point(:, (s1 - 1) * 24 + 13) = (/ int2_all_Cell(2, s1), int2_all_Cell(3, s1), int2_all_Cell(1, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 14) = (/ int2_all_Cell(2, s1), int2_all_Cell(3, s1), int2_all_Cell(5, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 15) = (/ int2_all_Cell(2, s1), int2_all_Cell(1, s1), int2_all_Cell(5, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 16) = (/ int2_all_Cell(3, s1), int2_all_Cell(1, s1), int2_all_Cell(5, s1) /)
				int2_all_tetraendron(:, (s1 - 1) * 6 + 4) = (/ (s1 - 1) * 24 + 13, (s1 - 1) * 24 + 14, (s1 - 1) * 24 + 15, (s1 - 1) * 24 + 16 /)
				
				int2_gran_point(:, (s1 - 1) * 24 + 17) = (/ int2_all_Cell(3, s1), int2_all_Cell(4, s1), int2_all_Cell(1, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 18) = (/ int2_all_Cell(3, s1), int2_all_Cell(4, s1), int2_all_Cell(5, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 19) = (/ int2_all_Cell(3, s1), int2_all_Cell(1, s1), int2_all_Cell(5, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 20) = (/ int2_all_Cell(4, s1), int2_all_Cell(1, s1), int2_all_Cell(5, s1) /)
				int2_all_tetraendron(:, (s1 - 1) * 6 + 5) = (/ (s1 - 1) * 24 + 17, (s1 - 1) * 24 + 18, (s1 - 1) * 24 + 19, (s1 - 1) * 24 + 20 /)
				
				int2_gran_point(:, (s1 - 1) * 24 + 21) = (/ int2_all_Cell(3, s1), int2_all_Cell(6, s1), int2_all_Cell(7, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 22) = (/ int2_all_Cell(3, s1), int2_all_Cell(6, s1), int2_all_Cell(5, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 23) = (/ int2_all_Cell(3, s1), int2_all_Cell(7, s1), int2_all_Cell(5, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 24) = (/ int2_all_Cell(6, s1), int2_all_Cell(7, s1), int2_all_Cell(5, s1) /)
				int2_all_tetraendron(:, (s1 - 1) * 6 + 6) = (/ (s1 - 1) * 24 + 21, (s1 - 1) * 24 + 22, (s1 - 1) * 24 + 23, (s1 - 1) * 24 + 24 /)
			
				int2_gran_sosed((s1 - 1) * 24 + 1) = (s1 - 1) * 6 + 4
				int2_gran_sosed((s1 - 1) * 24 + 4) = (s1 - 1) * 6 + 6
				
				int2_gran_sosed((s1 - 1) * 24 + 5) = (s1 - 1) * 6 + 3
				int2_gran_sosed((s1 - 1) * 24 + 7) = (s1 - 1) * 6 + 6
				
				int2_gran_sosed((s1 - 1) * 24 + 10) = (s1 - 1) * 6 + 5
				int2_gran_sosed((s1 - 1) * 24 + 11) = (s1 - 1) * 6 + 2
				
				int2_gran_sosed((s1 - 1) * 24 + 14) = (s1 - 1) * 6 + 1
				int2_gran_sosed((s1 - 1) * 24 + 16) = (s1 - 1) * 6 + 5
				
				int2_gran_sosed((s1 - 1) * 24 + 18) = (s1 - 1) * 6 + 3
				int2_gran_sosed((s1 - 1) * 24 + 19) = (s1 - 1) * 6 + 4
				
				int2_gran_sosed((s1 - 1) * 24 + 22) = (s1 - 1) * 6 + 1
				int2_gran_sosed((s1 - 1) * 24 + 23) = (s1 - 1) * 6 + 2
				
			end do
		end do
	end do
	
	
	N1 = size(int2_Cell_C(:, 1, 1))
	N2 = size(int2_Cell_C(1, :, 1))
	N3 = size(int2_Cell_C(1, 1, :))
	
	do k = 1, N3
		do j = 1, N2
			do i = 1, N1
				s1 = int2_Cell_C(i, j, k)
				if (s1 == 0) CYCLE
				int2_gran_point(:, (s1 - 1) * 24 + 1) = (/ int2_all_Cell(2, s1), int2_all_Cell(3, s1), int2_all_Cell(5, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 2) = (/ int2_all_Cell(2, s1), int2_all_Cell(3, s1), int2_all_Cell(6, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 3) = (/ int2_all_Cell(2, s1), int2_all_Cell(5, s1), int2_all_Cell(6, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 4) = (/ int2_all_Cell(3, s1), int2_all_Cell(5, s1), int2_all_Cell(6, s1) /)
				int2_all_tetraendron(:, (s1 - 1) * 6 + 1) = (/ (s1 - 1) * 24 + 1, (s1 - 1) * 24 + 2, (s1 - 1) * 24 + 3, (s1 - 1) * 24 + 4 /)
				
				int2_gran_point(:, (s1 - 1) * 24 + 5) = (/ int2_all_Cell(3, s1), int2_all_Cell(8, s1), int2_all_Cell(5, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 6) = (/ int2_all_Cell(3, s1), int2_all_Cell(8, s1), int2_all_Cell(7, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 7) = (/ int2_all_Cell(3, s1), int2_all_Cell(5, s1), int2_all_Cell(7, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 8) = (/ int2_all_Cell(8, s1), int2_all_Cell(5, s1), int2_all_Cell(7, s1) /)
				int2_all_tetraendron(:, (s1 - 1) * 6 + 2) = (/ (s1 - 1) * 24 + 5, (s1 - 1) * 24 + 6, (s1 - 1) * 24 + 7, (s1 - 1) * 24 + 8 /)
				
				int2_gran_point(:, (s1 - 1) * 24 + 9) = (/ int2_all_Cell(3, s1), int2_all_Cell(4, s1), int2_all_Cell(8, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 10) = (/ int2_all_Cell(3, s1), int2_all_Cell(4, s1), int2_all_Cell(5, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 11) = (/ int2_all_Cell(3, s1), int2_all_Cell(8, s1), int2_all_Cell(5, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 12) = (/ int2_all_Cell(4, s1), int2_all_Cell(8, s1), int2_all_Cell(5, s1) /)
				int2_all_tetraendron(:, (s1 - 1) * 6 + 3) = (/ (s1 - 1) * 24 + 9, (s1 - 1) * 24 + 10, (s1 - 1) * 24 + 11, (s1 - 1) * 24 + 12 /)
			
				int2_gran_point(:, (s1 - 1) * 24 + 13) = (/ int2_all_Cell(2, s1), int2_all_Cell(3, s1), int2_all_Cell(1, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 14) = (/ int2_all_Cell(2, s1), int2_all_Cell(3, s1), int2_all_Cell(5, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 15) = (/ int2_all_Cell(2, s1), int2_all_Cell(1, s1), int2_all_Cell(5, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 16) = (/ int2_all_Cell(3, s1), int2_all_Cell(1, s1), int2_all_Cell(5, s1) /)
				int2_all_tetraendron(:, (s1 - 1) * 6 + 4) = (/ (s1 - 1) * 24 + 13, (s1 - 1) * 24 + 14, (s1 - 1) * 24 + 15, (s1 - 1) * 24 + 16 /)
				
				int2_gran_point(:, (s1 - 1) * 24 + 17) = (/ int2_all_Cell(3, s1), int2_all_Cell(4, s1), int2_all_Cell(1, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 18) = (/ int2_all_Cell(3, s1), int2_all_Cell(4, s1), int2_all_Cell(5, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 19) = (/ int2_all_Cell(3, s1), int2_all_Cell(1, s1), int2_all_Cell(5, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 20) = (/ int2_all_Cell(4, s1), int2_all_Cell(1, s1), int2_all_Cell(5, s1) /)
				int2_all_tetraendron(:, (s1 - 1) * 6 + 5) = (/ (s1 - 1) * 24 + 17, (s1 - 1) * 24 + 18, (s1 - 1) * 24 + 19, (s1 - 1) * 24 + 20 /)
				
				int2_gran_point(:, (s1 - 1) * 24 + 21) = (/ int2_all_Cell(3, s1), int2_all_Cell(6, s1), int2_all_Cell(7, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 22) = (/ int2_all_Cell(3, s1), int2_all_Cell(6, s1), int2_all_Cell(5, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 23) = (/ int2_all_Cell(3, s1), int2_all_Cell(7, s1), int2_all_Cell(5, s1) /)
				int2_gran_point(:, (s1 - 1) * 24 + 24) = (/ int2_all_Cell(6, s1), int2_all_Cell(7, s1), int2_all_Cell(5, s1) /)
				int2_all_tetraendron(:, (s1 - 1) * 6 + 6) = (/ (s1 - 1) * 24 + 21, (s1 - 1) * 24 + 22, (s1 - 1) * 24 + 23, (s1 - 1) * 24 + 24 /)
			
				int2_gran_sosed((s1 - 1) * 24 + 1) = (s1 - 1) * 6 + 4
				int2_gran_sosed((s1 - 1) * 24 + 4) = (s1 - 1) * 6 + 6
				
				int2_gran_sosed((s1 - 1) * 24 + 5) = (s1 - 1) * 6 + 3
				int2_gran_sosed((s1 - 1) * 24 + 7) = (s1 - 1) * 6 + 6
				
				int2_gran_sosed((s1 - 1) * 24 + 10) = (s1 - 1) * 6 + 5
				int2_gran_sosed((s1 - 1) * 24 + 11) = (s1 - 1) * 6 + 2
				
				int2_gran_sosed((s1 - 1) * 24 + 14) = (s1 - 1) * 6 + 1
				int2_gran_sosed((s1 - 1) * 24 + 16) = (s1 - 1) * 6 + 5
				
				int2_gran_sosed((s1 - 1) * 24 + 18) = (s1 - 1) * 6 + 3
				int2_gran_sosed((s1 - 1) * 24 + 19) = (s1 - 1) * 6 + 4
				
				int2_gran_sosed((s1 - 1) * 24 + 22) = (s1 - 1) * 6 + 1
				int2_gran_sosed((s1 - 1) * 24 + 23) = (s1 - 1) * 6 + 2
				
			end do
		end do
	end do
	
	! Ещё не всех соседей тетраэдрам нашли
	! Но сначала нужно обнулить нулевые тетраэдры
	
	N1 = size(int2_all_tetraendron(1, :))
	
	loop2: do i = 1, N1
		do j = 1, 4
			s1 = int2_all_tetraendron(j, i)
			if (s1 == 0) CYCLE
			if(int2_gran_point(1, s1) == int2_gran_point(2, s1) .or. int2_gran_point(1, s1) == int2_gran_point(3, s1) &
				.or. int2_gran_point(2, s1) == int2_gran_point(3, s1)) then
				int2_all_tetraendron(:, i) = 0
				int2_gran_point(:, s1) = 0
				CYCLE loop2
			end if
		end do
	end do loop2
	
	! Находем соседей между ячейками (в А - группе)
	
	N1 = size(int2_Cell_A(:, 1, 1))
	N2 = size(int2_Cell_A(1, :, 1))
	N3 = size(int2_Cell_A(1, 1, :))
	
	do k = 1, N3
		do j = 1, N2 - 2
			do i = 1, N1 - 1
				
				if(i == par_n_TS - 1 .or. i == par_n_TS .or. i == par_n_TS + 1) CYCLE
				if(i == par_n_HP - 1 .or. i == par_n_HP .or. i == par_n_HP + 1) CYCLE
				if(i == par_n_BS - 1 .or. i == par_n_BS .or. i == par_n_BS + 1) CYCLE
				
				s1 = int2_Cell_A(i, j, k)
				s2 = int2_Cell_A(i + 1, j, k)
				
				if(s1 == 0 .or. s2 == 0) CYCLE

				if(j == N2 .and. i == par_n_TS + 1) CYCLE               ! Из-за тройной точки - эти отдельно построим
				if(j == N2 - 1 .and. i == par_n_TS + 1) CYCLE
				if(j == N2 .and. i == par_n_TS) CYCLE
				if(j == N2 - 1 .and. i == par_n_TS) CYCLE
				if(j == N2 .and. i == par_n_TS + 2) CYCLE
				if(j == N2 - 1 .and. i == par_n_TS + 2) CYCLE
				
				if(Int2_Comparison(int2_gran_point(:, (s2 - 1) * 24 + 3), int2_gran_point(:, (s1 - 1) * 24 + 6))) then
					int2_gran_sosed((s2 - 1) * 24 + 3) = (s1 - 1) * 6 + 2
					int2_gran_sosed((s1 - 1) * 24 + 6) = (s2 - 1) * 6 + 1
				end if
				
				
				if(Int2_Comparison(int2_gran_point(:, (s2 - 1) * 24 + 15), int2_gran_point(:, (s1 - 1) * 24 + 9))) then
					int2_gran_sosed((s2 - 1) * 24 + 15) = (s1 - 1) * 6 + 3
					int2_gran_sosed((s1 - 1) * 24 + 9) = (s2 - 1) * 6 + 4
				end if
				
				
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2 - 2
			do i = 1, N1
				
				if(i == par_n_TS - 1 .or. i == par_n_TS .or. i == par_n_TS + 1) CYCLE
				if(i == par_n_HP - 1 .or. i == par_n_HP .or. i == par_n_HP + 1) CYCLE
				if(i == par_n_BS - 1 .or. i == par_n_BS .or. i == par_n_BS + 1) CYCLE
				
				s1 = int2_Cell_A(i, j, k)
				s2 = int2_Cell_A(i, j + 1, k)
				
				if(s1 == 0 .or. s2 == 0) CYCLE
				
				if(Int2_Comparison(int2_gran_point(:, (s2 - 1) * 24 + 2), int2_gran_point(:, (s1 - 1) * 24 + 20))) then
					int2_gran_sosed((s2 - 1) * 24 + 2) = (s1 - 1) * 6 + 5
					int2_gran_sosed((s1 - 1) * 24 + 20) = (s2 - 1) * 6 + 1
				end if
				
				if(Int2_Comparison(int2_gran_point(:, (s2 - 1) * 24 + 21), int2_gran_point(:, (s1 - 1) * 24 + 12) )) then
					int2_gran_sosed((s2 - 1) * 24 + 21) = (s1 - 1) * 6 + 3
					int2_gran_sosed((s1 - 1) * 24 + 12) = (s2 - 1) * 6 + 6
				end if
				
			end do
		end do
	end do
	
	do k = 1, N3 - 1
		do j = 1, N2 - 2
			do i = 1, N1
				
				if(i == par_n_TS - 1 .or. i == par_n_TS .or. i == par_n_TS + 1) CYCLE
				if(i == par_n_HP - 1 .or. i == par_n_HP .or. i == par_n_HP + 1) CYCLE
				if(i == par_n_BS - 1 .or. i == par_n_BS .or. i == par_n_BS + 1) CYCLE
				
				
				s1 = int2_Cell_A(i, j, k)
				s2 = int2_Cell_A(i, j, k + 1)
				
				if(s1 == 0 .or. s2 == 0) CYCLE
				
				
				if(Int2_Comparison(int2_gran_point(:, (s2 - 1) * 24 + 13), int2_gran_point(:, (s1 - 1) * 24 + 24) )) then
					int2_gran_sosed((s2 - 1) * 24 + 13) = (s1 - 1) * 6 + 6
					int2_gran_sosed((s1 - 1) * 24 + 24) = (s2 - 1) * 6 + 4
				end if
				
				
				if(Int2_Comparison(int2_gran_point(:, (s2 - 1) * 24 + 8), int2_gran_point(:, (s1 - 1) * 24 + 17) )) then
					int2_gran_sosed((s2 - 1) * 24 + 8) = (s1 - 1) * 6 + 5
					int2_gran_sosed((s1 - 1) * 24 + 17) = (s2 - 1) * 6 + 2
				end if
				
			end do
		end do
	end do
	
	! Находим соседей между ячейками (в B - группе)
	
	N1 = size(int2_Cell_B(:, 1, 1))
	N2 = size(int2_Cell_B(1, :, 1))
	N3 = size(int2_Cell_B(1, 1, :))
	
	do k = 1, N3
		do j = 1, N2 - 2
			do i = 1, N1 - 1
				
				if(i == par_n_TS - 1 .or. i == par_n_TS .or. i == par_n_TS + 1) CYCLE
				
				s1 = int2_Cell_B(i, j, k)
				s2 = int2_Cell_B(i + 1, j, k)
				
				if(s1 == 0 .or. s2 == 0) CYCLE
				
				
				if(Int2_Comparison( int2_gran_point(1, (s2 - 1) * 24 + 15), int2_gran_point(1, (s1 - 1) * 24 + 9) )) then
					int2_gran_sosed((s2 - 1) * 24 + 15) = (s1 - 1) * 6 + 3
					int2_gran_sosed((s1 - 1) * 24 + 9) = (s2 - 1) * 6 + 4
				end if

				if( Int2_Comparison(int2_gran_point(1, (s2 - 1) * 24 + 3), int2_gran_point(1, (s1 - 1) * 24 + 6) )) then
					int2_gran_sosed((s2 - 1) * 24 + 3) = (s1 - 1) * 6 + 2
					int2_gran_sosed((s1 - 1) * 24 + 6) = (s2 - 1) * 6 + 1
				end if
				
					
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2 - 2
			do i = 1, N1
				
				if(i == par_n_TS - 1 .or. i == par_n_TS .or. i == par_n_TS + 1) CYCLE
				
				s1 = int2_Cell_B(i, j, k)
				s2 = int2_Cell_B(i, j + 1, k)
				
				if(s1 == 0 .or. s2 == 0) CYCLE
				
				if(Int2_Comparison( int2_gran_point(1, (s2 - 1) * 24 + 2), int2_gran_point(1, (s1 - 1) * 24 + 20) )) then
					int2_gran_sosed((s2 - 1) * 24 + 2) = (s1 - 1) * 6 + 5
					int2_gran_sosed((s1 - 1) * 24 + 20) = (s2 - 1) * 6 + 1
				end if
				
				if(Int2_Comparison( int2_gran_point(1, (s2 - 1) * 24 + 21), int2_gran_point(1, (s1 - 1) * 24 + 12) )) then
					int2_gran_sosed((s2 - 1) * 24 + 21) = (s1 - 1) * 6 + 3
					int2_gran_sosed((s1 - 1) * 24 + 12) = (s2 - 1) * 6 + 6
				end if
				
				
			end do
		end do
	end do
	
	do k = 1, N3 - 1
		do j = 1, N2 - 2
			do i = 1, N1
				
				if(i == par_n_TS - 1 .or. i == par_n_TS .or. i == par_n_TS + 1) CYCLE

				s1 = int2_Cell_B(i, j, k)
				s2 = int2_Cell_B(i, j, k + 1)
				
				if(s1 == 0 .or. s2 == 0) CYCLE
				

				if(Int2_Comparison( int2_gran_point(1, (s2 - 1) * 24 + 13), int2_gran_point(1, (s1 - 1) * 24 + 24) )) then
					int2_gran_sosed((s2 - 1) * 24 + 13) = (s1 - 1) * 6 + 6
					int2_gran_sosed((s1 - 1) * 24 + 24) = (s2 - 1) * 6 + 4
				end if
				
				if(Int2_Comparison( int2_gran_point(1, (s2 - 1) * 24 + 17), int2_gran_point(1, (s1 - 1) * 24 + 8) )) then
					int2_gran_sosed((s2 - 1) * 24 + 17) = (s1 - 1) * 6 + 2
					int2_gran_sosed((s1 - 1) * 24 + 8) = (s2 - 1) * 6 + 5
				end if
				
			end do
		end do
	end do
	
	! Находем соседей между ячейками (в C - группе)
	
	N1 = size(int2_Cell_C(:, 1, 1))
	N2 = size(int2_Cell_C(1, :, 1))
	N3 = size(int2_Cell_C(1, 1, :))
	
	do k = 1, N3
		do j = 1, N2
			do i = 1, N1 - 1
				s1 = int2_Cell_C(i, j, k)
				s2 = int2_Cell_C(i + 1, j, k)
				
				if(s1 == 0 .or. s2 == 0) CYCLE
				
				if(Int2_Comparison( int2_gran_point(1, (s2 - 1) * 24 + 15), int2_gran_point(1, (s1 - 1) * 24 + 9) )) then
					int2_gran_sosed((s2 - 1) * 24 + 15) = (s1 - 1) * 6 + 3
					int2_gran_sosed((s1 - 1) * 24 + 9) = (s2 - 1) * 6 + 4
				end if
				
				if(Int2_Comparison( int2_gran_point(1, (s2 - 1) * 24 + 3), int2_gran_point(1, (s1 - 1) * 24 + 6) )) then
					int2_gran_sosed((s2 - 1) * 24 + 3) = (s1 - 1) * 6 + 2
					int2_gran_sosed((s1 - 1) * 24 + 6) = (s2 - 1) * 6 + 1
				end if
				
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2 - 1
			do i = 1, N1
				s1 = int2_Cell_C(i, j, k)
				s2 = int2_Cell_C(i, j + 1, k)
				
				if(s1 == 0 .or. s2 == 0) CYCLE
				
				if(Int2_Comparison( int2_gran_point(1, (s2 - 1) * 24 + 2), int2_gran_point(1, (s1 - 1) * 24 + 20))) then
					int2_gran_sosed((s2 - 1) * 24 + 2) = (s1 - 1) * 6 + 5
					int2_gran_sosed((s1 - 1) * 24 + 20) = (s2 - 1) * 6 + 1
				end if
				
				if(Int2_Comparison( int2_gran_point(1, (s2 - 1) * 24 + 21), int2_gran_point(1, (s1 - 1) * 24 + 12) )) then
					int2_gran_sosed((s2 - 1) * 24 + 21) = (s1 - 1) * 6 + 3
					int2_gran_sosed((s1 - 1) * 24 + 12) = (s2 - 1) * 6 + 6
				end if
				
			end do
		end do
	end do
	
	do k = 1, N3 - 1
		do j = 1, N2
			do i = 1, N1
				s1 = int2_Cell_C(i, j, k)
				s2 = int2_Cell_C(i, j, k + 1)
				
				if(s1 == 0 .or. s2 == 0) CYCLE
				
				if(Int2_Comparison( int2_gran_point(1, (s2 - 1) * 24 + 13), int2_gran_point(1, (s1 - 1) * 24 + 24) )) then
					int2_gran_sosed((s2 - 1) * 24 + 13) = (s1 - 1) * 6 + 6
					int2_gran_sosed((s1 - 1) * 24 + 24) = (s2 - 1) * 6 + 4
				end if
				
				
				if(Int2_Comparison( int2_gran_point(1, (s2 - 1) * 24 + 17), int2_gran_point(1, (s1 - 1) * 24 + 8) )) then
					int2_gran_sosed((s2 - 1) * 24 + 17) = (s1 - 1) * 6 + 2
					int2_gran_sosed((s1 - 1) * 24 + 8) = (s2 - 1) * 6 + 5
				end if
				
				
			end do
		end do
	end do
	
	! Обнулим ссылки не несуществующие тетраэдры:
	N1 = size(int2_all_tetraendron(1, :))
	
	do i = 1, N1
		if(int2_all_tetraendron(1, i) == 0) then
			int2_all_tetraendron(1, i) = 0
			int2_all_tetraendron(2, i) = 0
			int2_all_tetraendron(3, i) = 0
			int2_all_tetraendron(4, i) = 0
			CYCLE !Нет тетраэдра
		end if
		
			
			do k = 1, 4 ! Бежим по граням тетраэдра
			s1 = int2_all_tetraendron(k, i)  ! Номер грани
			ii = int2_gran_sosed(s1)  ! Номер тетраэдра - соседа
			if (ii == 0) CYCLE
			
			if(int2_all_tetraendron(1, ii) == 0) then
				int2_all_tetraendron(1, ii) = 0
				int2_all_tetraendron(2, ii) = 0
				int2_all_tetraendron(3, ii) = 0
				int2_all_tetraendron(4, ii) = 0
				int2_gran_sosed(s1) = 0
				CYCLE
			end if
			end do
	end do 
	
	
! Обнулим грани несуществующих тетраэдров:
	N1 = size(int2_all_tetraendron(1, :))
	
	do i = 1, N1
		if(int2_all_tetraendron(1, i) == 0) then
			int2_gran_point(:, 4 * (i - 1) + 1) = 0
			int2_gran_point(:, 4 * (i - 1) + 2) = 0
			int2_gran_point(:, 4 * (i - 1) + 3) = 0
			int2_gran_point(:, 4 * (i - 1) + 4) = 0
		end if
	end do
	

	
	
	
	! Обнулим ссылки не несуществующие грани:
	N1 = size(int2_gran_point(1, :))
	
	do i = 1, N1
		if(int2_gran_point(1, i) == 0) CYCLE !Нет грани
			
			if(int2_gran_point(1, i) == int2_gran_point(2, i) .or. int2_gran_point(1, i) == int2_gran_point(3, i) .or. &
				int2_gran_point(2, i) == int2_gran_point(3, i)) int2_gran_point(:, i) = 0
	end do 
	
	print*, "Auto connect"
	! Теперь надо соединить между группами и тройную точку (предлагается автоматическое соединение)
	
	N1 = size(int2_all_Cell(1, :))
	if (.True.) then
	do i = 1, N1
		do j = 1, 6
			s1 = (i - 1) * 6 + j  ! Номер тетраэра
			
			
			if(int2_all_tetraendron(1, s1) == 0) CYCLE !Нет тетраэдра
			
			loop3: do k = 1, 4 ! Бежим по граням тетраэдра
				
				!if(s1 == 1195701) then
				!print*, int2_all_tetraendron(:, s1)
				!continue
				!end if
				
				
				s2 = (i - 1) * 24 + (j - 1) * 4 + k  ! Номер грани
				ijk = int2_gran_sosed(s2)   ! для отладки

				if (int2_gran_sosed(s2) /= 0) CYCLE ! Сосед уже найдён
				! Сосед не определён, нужно искать по координатам
				a1 = int2_coord(:, int2_gran_point(1, s2))
				a2 = int2_coord(:, int2_gran_point(2, s2))
				a3 = int2_coord(:, int2_gran_point(3, s2))
				
				do l = 1, 6  ! Бежим по соседям этой ячейки
					ii = int2_all_neighbours(l, i)
					if (ii == 0) CYCLE ! Соседа нет
					
					do jj = 1, 6  ! Тетраэдры соседа
						ss1 = (ii - 1) * 6 + jj  ! Номер тетраэдра
						if(int2_all_tetraendron(1, ss1) == 0) CYCLE !Нет тетраэдра
						
						do kk = 1, 4 ! Бежим по граням тетраэдра
							ss2 = (ii - 1) * 24 + (jj - 1) * 4 + kk  ! Номер грани
							b1 = int2_coord(:, int2_gran_point(1, ss2))
							b2 = int2_coord(:, int2_gran_point(2, ss2))
							b3 = int2_coord(:, int2_gran_point(3, ss2))
							if ( norm2(a1 - b1) < 0.00001 .or. norm2(a1 - b2) < 0.00001 .or. norm2(a1 - b3) < 0.00001) then
								if ( norm2(a2 - b1) < 0.00001 .or. norm2(a2 - b2) < 0.00001 .or. norm2(a2 - b3) < 0.00001) then
									if ( norm2(a3 - b1) < 0.00001 .or. norm2(a3 - b2) < 0.00001 .or. norm2(a3 - b3) < 0.00001) then
										int2_gran_sosed(s2) = ss1
											!if(ss1 == 1075918) then
											!	print*, a1, a2, a3
											!	print*, b1, b2, b3
											!	continue
											!end if
										CYCLE loop3
									else
										CYCLE
									end if
								else
									CYCLE
								end if
							else
								CYCLE
							end if
							
						end do
						
					end do
				end do
				
			end do loop3
		end do
	end do
	end if
	print*, "END auto connect"
	
	! Заполняем плоскости тетраедров
	
	N1 = size(int2_plane_tetraendron(1, 1, :))
	
	do i = 1, N1
		if (int2_all_tetraendron(1, i) == 0) CYCLE  ! Нет тетраэдра
		
		do j = 1, 4
			s1 = int2_all_tetraendron(j, i) ! Номер грани тетраедра
			a1 = int2_coord(:, int2_gran_point(1, s1))
			a2 = int2_coord(:, int2_gran_point(2, s1))
			a3 = int2_coord(:, int2_gran_point(3, s1))
			
			if (j == 1) int2_all_tetraendron_point(1, i) = int2_gran_point(1, s1)
			if (j == 1) int2_all_tetraendron_point(2, i) = int2_gran_point(2, s1)
			if (j == 1) int2_all_tetraendron_point(3, i) = int2_gran_point(3, s1)
			
			aa = a3 - a1
			bb = a2 - a1
		
			normal(1) = aa(2) * bb(3) - aa(3) * bb(2) 
			normal(2) = aa(3) * bb(1) - aa(1) * bb(3) 
			normal(3) = aa(1) * bb(2) - aa(2) * bb(1)  
			
			int2_plane_tetraendron(1:3, j, i) = normal
			int2_plane_tetraendron(4, j, i) = -DOT_PRODUCT(normal, a1)
			
			k = j + 1
			if (k > 4) k = 1
			s2 = int2_all_tetraendron(k, i) ! Номер другой грани тетраедра
			do ii = 1, 3
				if(int2_gran_point(ii, s2) /= int2_gran_point(1, s1) .and. int2_gran_point(ii, s2) /= int2_gran_point(2, s1) .and. &
				  int2_gran_point(ii, s2) /= int2_gran_point(3, s1)) then
					if (j == 1) int2_all_tetraendron_point(4, i) = int2_gran_point(ii, s2)
					if ( DOT_PRODUCT(int2_plane_tetraendron(1:3, j, i), int2_coord(:, int2_gran_point(ii, s2))) + &
						int2_plane_tetraendron(4, j, i) > 0) then
						int2_plane_tetraendron(:, j, i) = -int2_plane_tetraendron(:, j, i)
					end if
					EXIT
				end if
			end do
			
		
		end do
		
	end do
	
	! Вычислим объёмы тетраэдров ******************************************************************
	N1 = size(int2_all_tetraendron(1, :))
	do i = 1, N1
		
		if(int2_all_tetraendron_point(1, i) == 0 .or. int2_all_tetraendron_point(2, i) == 0 .or. &
			int2_all_tetraendron_point(3, i) == 0 .or. int2_all_tetraendron_point(4, i) == 0) then
			int2_all_tetraendron_point(1, i) = 0.0
			int2_all_tetraendron_point(2, i) = 0.0
			int2_all_tetraendron_point(3, i) = 0.0
			int2_all_tetraendron_point(4, i) = 0.0
			int2_all_Volume(i) = 0.0
			CYCLE
		end if
		
		a1 = int2_coord(:, int2_all_tetraendron_point(1, i))
		a2 = int2_coord(:, int2_all_tetraendron_point(2, i))
		a3 = int2_coord(:, int2_all_tetraendron_point(3, i))
		a4 = int2_coord(:, int2_all_tetraendron_point(4, i))
		
		b1 = a2 - a1
		b2 = a3 - a1
		
		b3(1) = b1(2) * b2(3) - b1(3) * b2(2)
        b3(2) = b1(3) * b2(1) - b1(1) * b2(3)
        b3(3) = b1(1) * b2(2) - b1(2) * b2(1)
		S = norm2(b3)
		b3 = b3/S
		S = S/2.0
		
		di = dabs(DOT_PRODUCT(a4, b3) - DOT_PRODUCT(b3, a1))
		
		int2_all_Volume(i) = S * di/3.0
		
		!if(int2_all_Volume(i) > 2.0) then
		!	print*, a1
		!	print*, a2
		!	print*, a3
		!	print*, a4
		!	print*, int2_all_Volume(i)
		!	continue
		!end if
		
		
	end do
	
	
	
	!print*, "proverka sosedey"
	!
	!N1 = size(int2_all_tetraendron(1, :))
	!
	!do i = 1, N1
	!	if(int2_all_tetraendron(1, i) == 0) CYCLE !Нет тетраэдра
	!		
	!	loop4: do k = 1, 4 ! Бежим по граням тетраэдра
	!		s1 = int2_all_tetraendron(k, i)  ! Номер грани
	!		ii = int2_gran_sosed(s1)  ! Номер тетраэдра - соседа
	!		if (ii == 0) CYCLE
	!		a1 = int2_coord(:, int2_gran_point(1, s1))
	!		a2 = int2_coord(:, int2_gran_point(2, s1))
	!		a3 = int2_coord(:, int2_gran_point(3, s1))
	!		do j = 1, 4 ! Бежим по граням тетраэдра - соседа
	!			ss2 = int2_all_tetraendron(j, ii) ! Номер грани
	!			b1 = int2_coord(:, int2_gran_point(1, ss2))
	!			b2 = int2_coord(:, int2_gran_point(2, ss2))
	!			b3 = int2_coord(:, int2_gran_point(3, ss2))
	!			if ( norm2(a1 - b1) < 0.00001 .or. norm2(a1 - b2) < 0.00001 .or. norm2(a1 - b3) < 0.00001) then
	!				if ( norm2(a2 - b1) < 0.00001 .or. norm2(a2 - b2) < 0.00001 .or. norm2(a2 - b3) < 0.00001) then
	!					if ( norm2(a3 - b1) < 0.00001 .or. norm2(a3 - b2) < 0.00001 .or. norm2(a3 - b3) < 0.00001) then
	!					int2_gran_sosed(s2) = ss1
	!					CYCLE loop4
	!					else
	!						CYCLE
	!					end if
	!				else
	!					CYCLE
	!				end if
	!			else
	!				CYCLE
	!			end if
	!			print*, "NET SOSEDA y SOSEDA"
	!		end do
	!	end do loop4
	!end do
	!
	!			
	!print*, "end proverka sosedey"
	
	
	!N1 = size(int2_Cell_A(:, 1, 1))
	!N2 = size(int2_Cell_A(1, :, 1))
	!N3 = size(int2_Cell_A(1, 1, :))
	
	!print*, "Size = ", size(int2_Cell_C) + size(int2_Cell_A) + size(int2_Cell_B), size(int2_all_Cell(1, :))
	
	!do k = 1, 1
	!	do j = N2, N2
	!		do i = par_n_TS + 1, par_n_TS + 1
	!			s1 = int2_Cell_A(i, j, k)
	!			if (s1 == 0) CYCLE
	!			
	!			print*, int2_all_Cell(:, s1)
	!			!if(s1 == 215713) then
	!			!	print*, i, j, k
	!			!	print*, int2_all_neighbours(:, s1)
	!			!else
	!			!	CYCLE
	!			!end if
	!			
	!			print*, "seriya tetraendrov  ", (s1 - 1) * 6
	!			do l = 1, 6
	!				print*, "*******************************************************************"
	!				print*, "l = ", l
	!				print*, int2_all_tetraendron(:, (s1 - 1) * 6 + l)
	!				do ii = 1, 4
	!					print*, "gran = ", ii
	!					ijk = int2_all_tetraendron(ii, (s1 - 1) * 6 + l)
	!					if (ijk /= 0) then
	!						print*, "Sosed = ", int2_gran_sosed(int2_all_tetraendron(ii, (s1 - 1) * 6 + l))
	!						print*, "yzel = ", int2_gran_point(:, ijk)
	!					else
	!						print*, "no gran"
	!					end if
	!				end do
	!				
	!				
	!			end do
	!			
	!		end do
	!	end do
	!end do
	
	end subroutine Int2_Initial
	
	subroutine Int2_Get_tetraedron(x, y, z, num)
	! Найти тетраедр, которому принадлежит точка
	! num по умолчанию должен быть равен 3
	implicit none
	real(8), intent(in) :: x, y, z
	integer(4), intent(in out) :: num  ! Тетраэдр, в котором предположительно находится точка
	integer(4) :: num2
	!integer(4) :: nummm(8), ijk
	integer(4) :: i, j, m, mk, rk
	real(8) :: r(3)
	
	!nummm = 0
	!ijk = 1
	mk = 0  ! Счётчик итераций
	rk = 0
	r = (/ x, y, z /)
	
	if (int2_all_tetraendron(1, num) == 0) then
		
		print*, "Takogo tetraedra net! ERROR  86tryjbyui98765rtyujhvdrtyu8765erthgg"
		return
	end if
	
	
11	CONTINUE	
	
	!nummm(ijk) = num
	!ijk = ijk + 1
	!if(ijk > 8) ijk = 1
	mk = mk + 1
	
	if(mk > 3000) then
		if(rk == 0) then
			r(1) = r(1) + 0.000001
			r(2) = r(2) + 0.000001
			mk = 1
			rk = 1
			!print*, "shevelim tochky"
			GO TO 11
		end if
		num = -1
		print*, "num = -1  oiuygvnmoiuhb"
		return
	end if
	
	m = 0
	do i = 1, 4
		if( DOT_PRODUCT(int2_plane_tetraendron(1:3, i, num), r ) + int2_plane_tetraendron(4, i, num) > 0 ) then
			m = 1
			num2 = int2_gran_sosed(int2_all_tetraendron(i, num))
			if(num2 == 0) CYCLE
			num = num2
			
			!if( ANY(nummm == num) ) then
			!	 num = -1
			!	 return
			!end if
		
			
			GO TO 11
			
		end if
	end do
	
	if(m == 1) then
		num = 0      ! Выход за пределы области
	end if
	
	
	
	end subroutine Int2_Get_tetraedron
	
	subroutine Int2_Get_tetraedron_inner(x, y, z, num)
	! Найти тетраедр, которому принадлежит точка
	! num по умолчанию должен быть равен 3
	! Эта функция вернёт близжайший тетраэдр, даже если точка за пределами области
	implicit none
	real(8), intent(in) :: x, y, z
	integer(4), intent(in out) :: num  ! Тетраэдр, в котором предположительно находится точка
	integer(4) :: num2
	!integer(4) :: nummm(8), ijk
	integer(4) :: i, j, m, mk, rk
	real(8) :: r(3)
	
	!nummm = 0
	!ijk = 1
	mk = 0  ! Счётчик итераций
	rk = 0
	r = (/ x, y, z /)
	
	if (int2_all_tetraendron(1, num) == 0) then
		
		print*, "Takogo tetraedra net! ERROR  86tryjbyui98765rtyujhvdrtyu8765erthgg"
		return
	end if
	
	
11	CONTINUE	
	
	!nummm(ijk) = num
	!ijk = ijk + 1
	!if(ijk > 8) ijk = 1
	mk = mk + 1
	
	if(mk > 3000) then
		if(rk == 0) then
			r(1) = r(1) + 0.000001
			r(2) = r(2) + 0.000001
			mk = 1
			rk = 1
			!print*, "shevelim tochky"
			GO TO 11
		end if
		num = -1
		print*, "num = -1  oiuygvnmoiuhb"
		return
	end if
	
	do i = 1, 4
		
		if( DOT_PRODUCT(int2_plane_tetraendron(1:3, i, num), r ) + int2_plane_tetraendron(4, i, num) > 0 ) then
			num2 = int2_gran_sosed(int2_all_tetraendron(i, num))
			if(num2 == 0) CYCLE
			num = num2
			GO TO 11
			
		end if
	end do
	
	end subroutine Int2_Get_tetraedron_inner
	
	logical pure function Int2_Comparison(a, b)
		! Сравнение наборов
		integer(4), intent(in) :: a(3), b(3)
		
		if(a(1) == b(1) .or. a(1) == b(2) .or. a(1) == b(3)) then
			if(a(2) == b(1) .or. a(2) == b(2) .or. a(2) == b(3)) then
				if(a(3) == b(1) .or. a(3) == b(2) .or. a(3) == b(3)) then
					Int2_Comparison = .True.
					return
				end if
			end if
		end if
		
		Int2_Comparison = .False.
		return
	
	end function Int2_Comparison
	
	
	subroutine Int2_Time_fly(r, V, time, num, next)
	! Variables
	real(8), intent(in) :: r(3), V(3)
	integer(4), intent(in) :: num  ! Номер тетраэдра, в котором находится частица
	integer(4), intent(out) :: next  ! Номер следующего тетраэдра
	real(8), intent(out) :: time

	integer(4) :: i, n2
	real(8) :: t2
	time = 10000000.0
	
	n2 = 1
	
	do i = 1, 4
		t2 = -(DOT_PRODUCT(r, int2_plane_tetraendron(1:3, i, num)) + int2_plane_tetraendron(4, i, num))&
			/(DOT_PRODUCT(V, int2_plane_tetraendron(1:3, i, num)))
		if (t2 < 0) CYCLE
		if(t2 < time) then
			time = t2
			n2 = i
		end if
	end do
	
	next = int2_gran_sosed(int2_all_tetraendron(n2, num))
	
	end subroutine Int2_Time_fly
	
	
	subroutine Int2_Get_par(x, y, z, num, PAR)
	! Найти тетраедр, которому принадлежит точка и получить значения параметров
	! num по умолчанию должен быть равен 3
	implicit none
	real(8), intent(in) :: x, y, z
	real(8), intent(out) :: PAR(9)     ! Выходные параметры
	integer(4), intent(in out) :: num  ! Тетраэдр, в котором предположительно находится точка (num по умолчанию должен быть равен 3)
	integer(4) :: i, j, r(3)
	real(8), dimension(4, 4) :: Minv
	real(8), dimension(4, 4) :: M
	real(8), dimension(4) :: work  ! work array for LAPACK
	real(8), dimension(1, 4) :: vec
	integer, dimension(4) :: ipiv   ! pivot indices
	integer :: n, info
	
	!$inter external DGETRF
	!$inter external DGETRI
	
	info = 0
	
	call Int2_Get_tetraedron(x, y, z, num)
	
	if(num == 0) then
		print*, "Net tetr"
		return
	end if
	
	
	M(:, 1) = 1.0_8
	M(1, 2:4) = int2_coord(:, int2_all_tetraendron_point(1, num) )
	M(2, 2:4) = int2_coord(:, int2_all_tetraendron_point(2, num) )
	M(3, 2:4) = int2_coord(:, int2_all_tetraendron_point(3, num) )
	M(4, 2:4) = int2_coord(:, int2_all_tetraendron_point(4, num) )
	
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
	vec(1, 2:4) = (/ x, y, z /)
	
	vec = MATMUL(vec, Minv)

	
	PAR = vec(1, 1) * int2_Cell_par(:, int2_all_tetraendron_point(1, num) ) + vec(1, 2) * int2_Cell_par(:, int2_all_tetraendron_point(2, num) ) + &
		vec(1, 3) * int2_Cell_par(:, int2_all_tetraendron_point(3, num) ) + vec(1, 4) * int2_Cell_par(:, int2_all_tetraendron_point(4, num) )
	
	
	end subroutine Int2_Get_par
	
	subroutine Int2_Save_bin(num)
	implicit none
    integer, intent(in) :: num
    character(len=5) :: name
	integer(4) :: n1, n2, n3
    
    write(unit=name,fmt='(i5.5)') num
    
    open(1, file = "int2_save_" // name // ".bin", FORM = 'BINARY')
	
	write(1) size(int2_coord(:, 1)), size(int2_coord(1, :))
	write(1) int2_coord
	
	write(1) size(int2_Cell_par(:, 1)), size(int2_Cell_par(1, :))
	write(1) int2_Cell_par
	
	write(1) size(int2_Cell_par2(:, 1)), size(int2_Cell_par2(1, :))
	write(1) int2_Cell_par2
	
	write(1) size(int2_Moment(:, 1, 1)), size(int2_Moment(1, :, 1)), size(int2_Moment(1, 1, :))
	write(1) int2_Moment
	
	write(1) size(int2_all_Cell(:, 1)), size(int2_all_Cell(1, :))
	write(1) int2_all_Cell
	
	write(1) size(int2_Cell_A(:, 1, 1)), size(int2_Cell_A(1, :, 1)), size(int2_Cell_A(1, 1, :))
	write(1) int2_Cell_A
	
	write(1) size(int2_Cell_B(:, 1, 1)), size(int2_Cell_B(1, :, 1)), size(int2_Cell_B(1, 1, :))
	write(1) int2_Cell_B
	
	write(1) size(int2_Cell_C(:, 1, 1)), size(int2_Cell_C(1, :, 1)), size(int2_Cell_C(1, 1, :))
	write(1) int2_Cell_C
	
	write(1) size(int2_Cell_center(:, 1)), size(int2_Cell_center(1, :))
	write(1) int2_Cell_center
	
	write(1) size(int2_all_neighbours(:, 1)), size(int2_all_neighbours(1, :))
	write(1) int2_all_neighbours
	
	write(1) size(int2_Point_A(:, 1, 1)), size(int2_Point_A(1, :, 1)), size(int2_Point_A(1, 1, :))
	write(1) int2_Point_A
	
	write(1) size(int2_Point_B(:, 1, 1)), size(int2_Point_B(1, :, 1)), size(int2_Point_B(1, 1, :))
	write(1) int2_Point_B
	
	write(1) size(int2_Point_C(:, 1, 1)), size(int2_Point_C(1, :, 1)), size(int2_Point_C(1, 1, :))
	write(1) int2_Point_C
	
	write(1) size(int2_all_tetraendron(:, 1)), size(int2_all_tetraendron(1, :))
	write(1) int2_all_tetraendron
	
	write(1) size(int2_all_tetraendron_point(:, 1)), size(int2_all_tetraendron_point(1, :))
	write(1) int2_all_tetraendron_point
	
	write(1) size(int2_all_tetraendron_matrix(:, 1, 1)), size(int2_all_tetraendron_matrix(1, :, 1)), size(int2_all_tetraendron_matrix(1, 1, :))
	write(1) int2_all_tetraendron_matrix
	
	write(1) size(int2_gran_point(:, 1)), size(int2_gran_point(1, :))
	write(1) int2_gran_point
	
	write(1) size(int2_gran_sosed(:))
	write(1) int2_gran_sosed
	
	write(1) size(int2_plane_tetraendron(:, 1, 1)), size(int2_plane_tetraendron(1, :, 1)), size(int2_plane_tetraendron(1, 1, :))
	write(1) int2_plane_tetraendron
	
	write(1) size(int2_all_Volume(:))
	write(1) int2_all_Volume
	
	
	write(1) 1
	write(1) size(int2_Moment_k(:, 1)), size(int2_Moment_k(1, :))
	write(1) int2_Moment_k
	
    write(1) 1
	write(1) size(int2_Cell_par_2(:, 1)), size(int2_Cell_par_2(1, :))
	write(1) int2_Cell_par_2
	
	
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
	
	end subroutine Int2_Save_bin
	
	subroutine Helium_off()
		integer :: zone, i
	
		do i = 1, size(int2_Cell_par(1, :))
			zone = int2_Cell_par2(1, i)
			int2_Cell_par(1, i) = int2_Cell_par(1, i) - int2_Cell_par_2(1, i)
			
			if(zone == 1 .or. zone == 2) then
				int2_Cell_par(5, i) = int2_Cell_par(5, i) / (2.0 * int2_Cell_par(1, i) + 1.5 * int2_Cell_par_2(1, i)) * 2.0 * int2_Cell_par(1, i)
			else
				int2_Cell_par(5, i) = int2_Cell_par(5, i) / (2.0 * int2_Cell_par(1, i) + int2_Cell_par_2(1, i)) * 2.0 * int2_Cell_par(1, i)
			end if
		end do
		
	end subroutine Helium_off
	
	subroutine Helium_on()
		integer :: zone, i
	
		do i = 1, size(int2_Cell_par(1, :))
			zone = int2_Cell_par2(1, i)
			if(zone == 1 .or. zone == 2) then
				int2_Cell_par(5, i) = int2_Cell_par(5, i) * (2.0 * int2_Cell_par(1, i) + 1.5 * int2_Cell_par_2(1, i)) / (2.0 * int2_Cell_par(1, i))
			else
				int2_Cell_par(5, i) = int2_Cell_par(5, i) * (2.0 * int2_Cell_par(1, i) + int2_Cell_par_2(1, i)) / (2.0 * int2_Cell_par(1, i))
			end if
			int2_Cell_par(1, i) = int2_Cell_par(1, i) + int2_Cell_par_2(1, i)
		end do
	end subroutine Helium_on
	
	subroutine Int2_Save_interpol_for_all_MHD(num)
	! Сохраняет интерполяционные файлы для всех желующих
	use STORAGE
    use GEO_PARAM
    implicit none
    integer, intent(in) :: num
    character(len=5) :: name
	integer :: zone, i
    
    write(unit=name, fmt='(i5.5)') num
    
    open(1, file = "int2_save_for_all_MHD" // name // ".bin", FORM = 'BINARY')
	
    
    write(1) size(int2_all_tetraendron(:, 1)), size(int2_all_tetraendron(1, :))
    write(1) int2_all_tetraendron
	
	write(1) size(int2_plane_tetraendron(:, 1, 1)), size(int2_plane_tetraendron(1, :, 1)), size(int2_plane_tetraendron(1, 1, :))
    write(1) int2_plane_tetraendron
	
	write(1) size(int2_gran_sosed(:))
    write(1) int2_gran_sosed
	
	write(1) size(int2_all_tetraendron_matrix(:, 1, 1)), size(int2_all_tetraendron_matrix(1, :, 1)), size(int2_all_tetraendron_matrix(1, 1, :))
    write(1) int2_all_tetraendron_matrix
	
	! He
	do i = 1, size(int2_Cell_par(1, :))
		zone = int2_Cell_par2(1, i)
		int2_Cell_par(1, i) = int2_Cell_par(1, i) - int2_Cell_par_2(1, i)
			
		if(zone == 1 .or. zone == 2) then
			int2_Cell_par(5, i) = int2_Cell_par(5, i) / (2.0 * int2_Cell_par(1, i) + 1.5 * int2_Cell_par_2(1, i)) * 2.0 * int2_Cell_par(1, i)
		else
			int2_Cell_par(5, i) = int2_Cell_par(5, i) / (2.0 * int2_Cell_par(1, i) + int2_Cell_par_2(1, i)) * 2.0 * int2_Cell_par(1, i)
		end if
	end do
	
	
	write(1) size(int2_Cell_par(:, 1)), size(int2_Cell_par(1, :))
    write(1) int2_Cell_par
	
	write(1) size(int2_Cell_par_2(:, 1)), size(int2_Cell_par_2(1, :))
    write(1) int2_Cell_par_2
	
	do i = 1, size(int2_Cell_par(1, :))
		zone = int2_Cell_par2(1, i)
		if(zone == 1 .or. zone == 2) then
			int2_Cell_par(5, i) = int2_Cell_par(5, i) * (2.0 * int2_Cell_par(1, i) + 1.5 * int2_Cell_par_2(1, i)) / (2.0 * int2_Cell_par(1, i))
		else
			int2_Cell_par(5, i) = int2_Cell_par(5, i) * (2.0 * int2_Cell_par(1, i) + int2_Cell_par_2(1, i)) / (2.0 * int2_Cell_par(1, i))
		end if
		int2_Cell_par(1, i) = int2_Cell_par(1, i) + int2_Cell_par_2(1, i)
	end do
	
	write(1) size(int2_all_tetraendron_point(:, 1)), size(int2_all_tetraendron_point(1, :))
    write(1) int2_all_tetraendron_point
	
	write(1) size(int2_Cell_par2(:, 1)), size(int2_Cell_par2(1, :))  ! ЗОНА
    write(1) int2_Cell_par2
	
	close(1)
	
	
	! Body of Save_interpol_for_all
	
	end subroutine Int2_Save_interpol_for_all_MHD
	
	
	subroutine Int2_Save_interpol_for_all_MK(num)
	! Сохраняет интерполяционные файлы для всех желующих
	use STORAGE
    use GEO_PARAM
    implicit none
    integer, intent(in) :: num
    character(len=5) :: name
	integer :: zone, i
    
    write(unit=name, fmt='(i5.5)') num
    
    open(1, file = "int2_save_for_all_MK" // name // ".bin", FORM = 'BINARY')
	
    
    write(1) size(int2_all_tetraendron(:, 1)), size(int2_all_tetraendron(1, :))
    write(1) int2_all_tetraendron
	
	write(1) size(int2_plane_tetraendron(:, 1, 1)), size(int2_plane_tetraendron(1, :, 1)), size(int2_plane_tetraendron(1, 1, :))
    write(1) int2_plane_tetraendron
	
	write(1) size(int2_gran_sosed(:))
    write(1) int2_gran_sosed
	
	write(1) size(int2_all_tetraendron_matrix(:, 1, 1)), size(int2_all_tetraendron_matrix(1, :, 1)), size(int2_all_tetraendron_matrix(1, 1, :))
    write(1) int2_all_tetraendron_matrix
	
	
	write(1) size(int2_Moment(:, 1, 1)), size(int2_Moment(1, :, 1)), size(int2_Moment(1, 1, :))  ! ЗОНА
    write(1) int2_Moment
	
	
	write(1) size(int2_all_tetraendron_point(:, 1)), size(int2_all_tetraendron_point(1, :))
    write(1) int2_all_tetraendron_point
	
	write(1) size(int2_Cell_par2(:, 1)), size(int2_Cell_par2(1, :))  ! ЗОНА
    write(1) int2_Cell_par2
	
	close(1)
	
	
	! Body of Save_interpol_for_all
	
	end subroutine Int2_Save_interpol_for_all_MK
	
	subroutine Int2_Read_bin(num)
	implicit none
	integer, intent(in) :: num
    character(len=5) :: name
    integer(4) :: n
    logical :: exists
	integer(4) :: n1, n2, n3
    
    write(unit=name,fmt='(i5.5)') num
    
    inquire(file="int2_save_" // name // ".bin", exist=exists)
    
    if (exists == .False.) then
		pause "net faila!!!  edhpifs56uhekwkwiuwyvrr"
        STOP "net faila!!!"
    end if
    
    open(1, file = "int2_save_" // name // ".bin", FORM = 'BINARY', ACTION = "READ")
	
	
	read(1) n1, n2
	allocate( int2_coord(n1, n2) )
	read(1) int2_coord
	
	read(1) n1, n2
	allocate( int2_Cell_par(n1, n2) )
	read(1) int2_Cell_par
	
	read(1) n1, n2
	allocate( int2_Cell_par2(n1, n2) )
	read(1) int2_Cell_par2
	
	read(1) n1, n2, n3
	allocate( int2_Moment(n1, n2, n3) )
	read(1) int2_Moment
	par_n_moment = n1
	
	read(1) n1, n2
	allocate( int2_all_Cell(n1, n2) )
	read(1) int2_all_Cell
	
	read(1) n1, n2, n3
	allocate( int2_Cell_A(n1, n2, n3) )
	read(1) int2_Cell_A
	
	read(1) n1, n2, n3
	allocate( int2_Cell_B(n1, n2, n3) )
	read(1) int2_Cell_B
	
	read(1) n1, n2, n3
	allocate( int2_Cell_C(n1, n2, n3) )
	read(1) int2_Cell_C
	
	read(1) n1, n2
	allocate( int2_Cell_center(n1, n2) )
	read(1) int2_Cell_center
	
	read(1) n1, n2
	allocate( int2_all_neighbours(n1, n2) )
	read(1) int2_all_neighbours
	
	read(1) n1, n2, n3
	allocate( int2_Point_A(n1, n2, n3) )
	read(1) int2_Point_A
	
	read(1) n1, n2, n3
	allocate( int2_Point_B(n1, n2, n3) )
	read(1) int2_Point_B
	
	read(1) n1, n2, n3
	allocate( int2_Point_C(n1, n2, n3) )
	read(1) int2_Point_C
	
	read(1) n1, n2
	allocate( int2_all_tetraendron(n1, n2) )
	read(1) int2_all_tetraendron
	
	read(1) n1, n2
	allocate( int2_all_tetraendron_point(n1, n2) )
	read(1) int2_all_tetraendron_point
	
	read(1) n1, n2, n3
	allocate( int2_all_tetraendron_matrix(n1, n2, n3) )
	read(1) int2_all_tetraendron_matrix
	
	read(1) n1, n2
	allocate( int2_gran_point(n1, n2) )
	read(1) int2_gran_point
	
	read(1) n1
	allocate( int2_gran_sosed(n1) )
	read(1) int2_gran_sosed
	
	read(1) n1, n2, n3
	allocate( int2_plane_tetraendron(n1, n2, n3) )
	read(1) int2_plane_tetraendron
	
	read(1) n1
	allocate( int2_all_Volume(n1) )
	read(1) int2_all_Volume
	
	read(1) n
	if(n == 1) then
		read(1) n1, n2
		allocate( int2_Moment_k(n1, n2) )
		read(1) int2_Moment_k
	end if
	
	if( allocated(int2_Moment_k) == .False.) then
		allocate(int2_Moment_k(5, size(int2_Point_A) + size(int2_Point_B) + size(int2_Point_C)))
		int2_Moment_k = 1.0
	end if
	
	read(1) n
	if(n == 1) then
		read(1) n1, n2
		allocate( int2_Cell_par_2(n1, n2) )
		read(1) int2_Cell_par_2
	end if
	
	if( allocated(int2_Cell_par_2) == .False.) then
		allocate(int2_Cell_par_2(1, size(int2_Cell_par(1, :)) ))
		int2_Cell_par_2 = 0.0
	end if
	
	close(1)
	
	end subroutine Int2_Read_bin
	
	subroutine Int2_culc_k()
	! Посчитаем коэффициенты отношения источников мультифлюида к Монте-Карло
	! Также меняет температуру атомов на давление
	integer(4) :: i, j
	real(8) :: sourse(5,5), ss
	
	! Делаем из температуры давление
	do i = 1, size(int2_Moment(1, 1, :))
		int2_Moment(5, 1:4, i) = 0.5 * int2_Moment(5, 1:4, i) * int2_Moment(1, 1:4, i)
	end do
	
	! Бежим по точкам
	do i = 1, size(int2_Moment(1, 1, :))
		call Calc_sourse_MF(int2_Cell_par(:, i), int2_Moment(1:5, 1:4, i), sourse, 1)
		
		!! Считаем коэффициенты (для трёх импульсов и энергии)
		do j = 2, 5
			if( dabs(sourse(j, 1)) > 0.00001) then
				ss = sum(int2_Moment(4 + j, :, i))
				int2_Moment_k(j, i) = ss/sourse(j, 1)
				if(int2_Moment_k(j, i) < 0.3 .or. int2_Moment_k(j, i) > 3.0) int2_Moment_k(j, i) = 1.0
			end if
		end do
		
		int2_Moment_k(1, i) = sum(int2_Moment(19, :, i))
		
	end do
	
	
	! Body of Int2_culc_k
	
	end subroutine Int2_culc_k
	
	subroutine Int2_Set_interpol_matrix()
	! Заполняем интерполяционные матрицы для всех ячееек и записываем в память
	
	integer(4):: num  ! Тетраэдр, в котором предположительно находится точка (num по умолчанию должен быть равен 3)
	integer(4) :: i, j, r(3)
	real(8), dimension(4, 4) :: Minv
	real(8), dimension(4, 4) :: M
	real(8), dimension(4) :: work  ! work array for LAPACK
	real(8), dimension(1, 4) :: vec
	integer, dimension(4) :: ipiv   ! pivot indices
	integer :: n, info
	
	!$inter external DGETRF
	!$inter external DGETRI
	
	info = 0
	
	do num = 1, size(int2_all_tetraendron_point(1, :))
		
		if(int2_all_tetraendron(1, num) == 0) CYCLE
		
		M(:, 1) = 1.0_8
		M(1, 2:4) = int2_coord(:, int2_all_tetraendron_point(1, num) )
		M(2, 2:4) = int2_coord(:, int2_all_tetraendron_point(2, num) )
		M(3, 2:4) = int2_coord(:, int2_all_tetraendron_point(3, num) )
		M(4, 2:4) = int2_coord(:, int2_all_tetraendron_point(4, num) )
	
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
		
		int2_all_tetraendron_matrix(:,:, num) = Minv
	end do
	
	end subroutine Int2_Set_interpol_matrix
	
	subroutine Int2_Get_par_fast(x, y, z, num, PAR, PAR_MOMENT, PAR_k)
	! Найти тетраедр, которому принадлежит точка и получить значения параметров
	! В отличие от медленной версии, эта не вычисляет матрицу интерполяции каждый раз, 
	! предполагается, что матрицы лежат в памяти
	! num по умолчанию должен быть равен 3
	implicit none
	real(8), intent(in) :: x, y, z
	real(8), intent(out) :: PAR(9)     ! Выходные параметры
	integer(4), intent(in out) :: num  ! Тетраэдр, в котором предположительно находится точка (num по умолчанию должен быть равен 3)
	real(8), intent(out), optional :: PAR_MOMENT(par_n_moment, par_n_sort)
	real(8), intent(out), optional :: PAR_k(5)
	
	real(8), dimension(4, 4) :: Minv
	real(8), dimension(1, 4) :: vec
	
	
	call Int2_Get_tetraedron(x, y, z, num)
	
	if(num < 1) then
		print*, "Net tetr"
		return
	end if
	
	Minv = int2_all_tetraendron_matrix(:, :, num)
	
	vec(1, 1) = 1.0_8
	vec(1, 2:4) = (/ x, y, z /)
	
	vec = MATMUL(vec, Minv)

	
	PAR = vec(1, 1) * int2_Cell_par(:, int2_all_tetraendron_point(1, num) ) + vec(1, 2) * int2_Cell_par(:, int2_all_tetraendron_point(2, num) ) + &
		vec(1, 3) * int2_Cell_par(:, int2_all_tetraendron_point(3, num) ) + vec(1, 4) * int2_Cell_par(:, int2_all_tetraendron_point(4, num) )
	
	if(present(PAR_MOMENT)) then
		PAR_MOMENT = vec(1, 1) * int2_Moment(:, :, int2_all_tetraendron_point(1, num) ) + vec(1, 2) * int2_Moment(:, :, int2_all_tetraendron_point(2, num) ) + &
		vec(1, 3) * int2_Moment(:, :, int2_all_tetraendron_point(3, num) ) + vec(1, 4) * int2_Moment(:, :, int2_all_tetraendron_point(4, num) )
	end if
	
	if(present(PAR_k)) then
		PAR_k = vec(1, 1) * int2_Moment_k(:, int2_all_tetraendron_point(1, num) ) + vec(1, 2) * int2_Moment_k(:, int2_all_tetraendron_point(2, num) ) + &
		vec(1, 3) * int2_Moment_k(:, int2_all_tetraendron_point(3, num) ) + vec(1, 4) * int2_Moment_k(:, int2_all_tetraendron_point(4, num) )
	end if
	
	
	end subroutine Int2_Get_par_fast
	
	subroutine Int2_Get_par_fast2(x, y, z, num, PAR)
	! Здесь точно известно в каком тетраэдое точка
	! Найти тетраедр, которому принадлежит точка и получить значения параметров
	! В отличие от медленной версии, эта не вычисляет матрицу интерполяции каждый раз, 
	! предполагается, что матрицы лежат в памяти
	! num по умолчанию должен быть равен 3
	implicit none
	real(8), intent(in) :: x, y, z
	real(8), intent(out) :: PAR(9)     ! Выходные параметры
	integer(4), intent(in) :: num  ! Тетраэдр, в котором предположительно находится точка (num по умолчанию должен быть равен 3)
	
	real(8), dimension(4, 4) :: Minv
	real(8), dimension(1, 4) :: vec
	
	
	
	Minv = int2_all_tetraendron_matrix(:, :, num)
	
	vec(1, 1) = 1.0_8
	vec(1, 2:4) = (/ x, y, z /)
	
	vec = MATMUL(vec, Minv)

	
	PAR = vec(1, 1) * int2_Cell_par(:, int2_all_tetraendron_point(1, num) ) + vec(1, 2) * int2_Cell_par(:, int2_all_tetraendron_point(2, num) ) + &
		vec(1, 3) * int2_Cell_par(:, int2_all_tetraendron_point(3, num) ) + vec(1, 4) * int2_Cell_par(:, int2_all_tetraendron_point(4, num) )
	
	
	end subroutine Int2_Get_par_fast2
	
	logical pure function Belong_tetraedron(x, y, z, num)
	! Принадлежит ли точка тетраедру 
	implicit none
	real(8), intent(in) :: x, y, z
	integer(4), intent(in) :: num  ! Тетраэдр, в котором предположительно находится точка
	integer(4) :: i, j
	do i = 1, 4
		
		if( DOT_PRODUCT(int2_plane_tetraendron(1:3, i, num), (/ x, y, z /) ) + &
						int2_plane_tetraendron(4, i, num) > 0 ) then
			Belong_tetraedron = .False.
			return
		end if
	end do
	
	Belong_tetraedron = .True.
	return
	end function Belong_tetraedron
	
	subroutine Int_2_Cut_tetr(cell, A, B, C, D, bb, CUT)
	! Разрез тетраэдра плоскостью, результат - три точки
		integer(4), intent(in) :: cell    ! Какой тетрадр режим
		real(8), intent(in) :: A, B, C, D  ! Уравнение плоскости
		logical, intent(out) :: bb
		real(8), intent(out) :: CUT(3, 4)  ! (3 координаты, 4 точки)
		
		real(8) :: a1(3), a2(3), a3(3), a4(3), n(3), t, AA(3), p(3, 4)
		integer(4) :: i1, i2, i3, i4, M(2, 6), i, k
		
		if(int2_all_tetraendron_point(1, cell) == 0) then
			bb = .False.   ! Плоскость не пересекает тетраэдр
			return
		end if
		
		if(int2_all_tetraendron_point(2, cell) == 0) then
			bb = .False.   ! Плоскость не пересекает тетраэдр
			return
		end if
		
		if(int2_all_tetraendron_point(3, cell) == 0) then
			bb = .False.   ! Плоскость не пересекает тетраэдр
			return
		end if
		
		if(int2_all_tetraendron_point(4, cell) == 0) then
			bb = .False.   ! Плоскость не пересекает тетраэдр
			return
		end if
		
		
		
		a1 = int2_coord(:, int2_all_tetraendron_point(1, cell))
		a2 = int2_coord(:, int2_all_tetraendron_point(2, cell))
		a3 = int2_coord(:, int2_all_tetraendron_point(3, cell))
		a4 = int2_coord(:, int2_all_tetraendron_point(4, cell))
		
		i1 = signum(a1(1) * A + a1(2) * B + a1(3) * C + D)
		i2 = signum(a2(1) * A + a2(2) * B + a2(3) * C + D)
		i3 = signum(a3(1) * A + a3(2) * B + a3(3) * C + D)
		i4 = signum(a4(1) * A + a4(2) * B + a4(3) * C + D)
		
		if(i1 == 0 .and. i2 == 1 .and. i3 == 1 .and. i4 == 1) GO TO 101
		if(i1 == 0 .and. i2 == -1 .and. i3 == -1 .and. i4 == -1) GO TO 101
		if(i1 == 1 .and. i2 == 0 .and. i3 == 1 .and. i4 == 1) GO TO 101
		if(i1 == -1 .and. i2 == 0 .and. i3 == -1 .and. i4 == -1) GO TO 101
		if(i1 == 1 .and. i2 == 1 .and. i3 == 0 .and. i4 == 1) GO TO 101
		if(i1 == -1 .and. i2 == -1 .and. i3 == 0 .and. i4 == -1) GO TO 101
		if(i1 == 1 .and. i2 == 1 .and. i3 == 1 .and. i4 == 0) GO TO 101
		if(i1 == -1 .and. i2 == -1 .and. i3 == -1 .and. i4 == 0) GO TO 101
		
		
		if(i1 /= 0) then
			if(i1 == i2) then
				i2 = signum(a3(1) * A + a3(2) * B + a3(3) * C + D)
				if(i1 == i2) then
					i2 = signum(a4(1) * A + a4(2) * B + a4(3) * C + D)
					if(i1 == i2) then
						bb = .False.   ! Плоскость не пересекает тетраэдр
						return
					end if
				end if
			end if
		end if
		
		! Плоскость пересекает тетраэдр
		bb = .True.
		AA = (/ A, B, C /)
		p(:, 1) = a1
		p(:, 2) = a2
		p(:, 3) = a3
		p(:, 4) = a4
			
		M(:, 1) = (/ 1, 2 /)
		M(:, 2) = (/ 1, 3 /)
		M(:, 3) = (/ 1, 4 /)
		M(:, 4) = (/ 2, 4 /)
		M(:, 5) = (/ 2, 3 /)
		M(:, 6) = (/ 3, 4 /)
		
		k = 1
		do i = 1, 6
			i1 = M(1, i)
			i2 = M(2, i)
			if(signum(DOT_PRODUCT(AA, p(:, i1)) + D) == signum(DOT_PRODUCT(AA, p(:, i2)) + D)) CYCLE
			n = p(:, i2) - p(:, i1)
			t = (-D - DOT_PRODUCT(AA, p(:, i1)))/DOT_PRODUCT(AA, n)
			CUT(:, k) = p(:, i1) + t * n
			k = k + 1
			if(k == 5) EXIT
		end do
		
		if(k == 4) then
			CUT(:, 4) = CUT(:, 3)
		end if
		
		
		return
		
101		continue
		bb = .False.   ! Плоскость не пересекает тетраэдр
		return
	
	end subroutine Int_2_Cut_tetr
	
	subroutine Int_2_Print_par_2D_set()
		integer :: i, j, N1, N2, p
		real(8) :: PAR(9), PAR_MK(par_n_moment, par_n_sort), PAR_K(5)
		real(8) :: PAR_MK_SUM(par_n_moment)
		character(len=2) :: name, name2

		open(1, file = 'print_par_2D_SET.txt')
		write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z', 'rho', 'u', 'v', 'w', 'p',"
		write(1,*) "'bx', 'by', 'bz', 'Q',"
		
		do i = 1, par_n_moment
			write(unit=name, fmt='(i2.2)') i
			write(1,*) "'SUM_" // name //"',"
		end do
		
		do j = 1, par_n_sort
			write(unit=name2, fmt='(i2.2)') j
			do i = 1, 9
				write(unit=name, fmt='(i2.2)') i
				write(1,*) "'Sort_" // name2 // "_" // name //"',"
			end do
		end do
		
		write(1,*) "'k2', k3', 'k4', 'k5'"
		write(1,*) ", ZONE T= 'HP'"
		
		
		N1 = size(int2_Point_A(:, 1, 1)) 
		N2 = size(int2_Point_A(1, :, 1)) 
		
		do i = 1, N1
			do j = 1, N2
				p = int2_Point_A(i, j, 1)
				
				PAR = int2_Cell_par(:, p)
				PAR_MK = int2_Moment(:, :, p)
				PAR_K = int2_Moment_k(:, p)
				call Int2_sum_moment(PAR_MK, PAR_MK_SUM)
				
				write(1,*) int2_coord(:, p), PAR, PAR_MK_SUM, PAR_MK(1:9, 1), PAR_MK(1:9, 2), PAR_MK(1:9, 3), PAR_MK(1:9, 4), PAR_K(2:5)
				write(1,*) ""
				
			end do
		end do
		
		N1 = size(int2_Point_B(:, 1, 1)) 
		N2 = size(int2_Point_B(1, :, 1)) 
		
		do i = 1, N1
			do j = 1, N2
				p = int2_Point_B(i, j, 1)
				
				PAR = int2_Cell_par(:, p)
				PAR_MK = int2_Moment(:, :, p)
				PAR_K = int2_Moment_k(:, p)
				call Int2_sum_moment(PAR_MK, PAR_MK_SUM)
				
				write(1,*) int2_coord(:, p), PAR, PAR_MK_SUM, PAR_MK(1:9, 1), PAR_MK(1:9, 2), PAR_MK(1:9, 3), PAR_MK(1:9, 4), PAR_K(2:5)
				write(1,*) ""
				
			end do
		end do
		
		N1 = size(int2_Point_C(:, 1, 1)) 
		N2 = size(int2_Point_C(1, :, 1)) 
		
		do i = 1, N1
			do j = 1, N2
				p = int2_Point_C(i, j, 1)
				
				PAR = int2_Cell_par(:, p)
				PAR_MK = int2_Moment(:, :, p)
				PAR_K = int2_Moment_k(:, p)
				call Int2_sum_moment(PAR_MK, PAR_MK_SUM)
				
				write(1,*) int2_coord(:, p), PAR, PAR_MK_SUM, PAR_MK(1:9, 1), PAR_MK(1:9, 2), PAR_MK(1:9, 3), PAR_MK(1:9, 4), PAR_K(2:5)
				write(1,*) ""
				
			end do
		end do

		close(1)
	
	end subroutine Int_2_Print_par_2D_set
	
	subroutine Int_2_Print_par_2D(A, B, C, D, nn)  ! Печатает 2Д сетку с линиями в Техплот
	
		real(8), intent(in) :: A, B, C, D
		integer(4), intent(in) :: nn
		integer :: i, n, j, num
		real(8) :: Mach, PAR(9), PAR_MK(par_n_moment, par_n_sort), a1(3), a2(3), b1(3), b2(3), b3(3), S, PAR_K(5)
		real(8) :: PAR_MK_SUM(par_n_moment)
		real(8), allocatable :: CUT(:, :, :)
		logical :: bb
		character(len=5) :: name

		write(unit=name, fmt='(i5.5)') nn
		
		
		allocate(CUT(3, 4, 500000))
		
		num = 3
		n = 1
		do i = 1, size(int2_all_tetraendron(1, :))
			call Int_2_Cut_tetr(i, A, B, C, D, bb, CUT(:, :, n))
			
			if (bb == .True.) then
				! Если четвёрка неправильно ориентирована
				a1 = CUT(:, 2, n) - CUT(:, 1, n)
				a2 = CUT(:, 4, n) - CUT(:, 3, n)
				
				if(DOT_PRODUCT(a1, a2) > 0) then
					a1 = CUT(:, 4, n)
					CUT(:, 4, n) = CUT(:, 3, n)
					CUT(:, 3, n) = a1
				end if

				b1 = CUT(:, 3, n) - CUT(:, 1, n)
				b2 = CUT(:, 4, n) - CUT(:, 2, n)
				b3(1) = b1(2) * b2(3) - b1(3) * b2(2)
				b3(2) = b1(3) * b2(1) - b1(1) * b2(3)
				b3(3) = b1(1) * b2(2) - b1(2) * b2(1)
				S = norm2(b3)
				if(S < 0.000001) CYCLE
				
				!b1 = (CUT(:, 1, n) + CUT(:, 2, n) + CUT(:, 3, n) + CUT(:, 4, n))/4.0
				!do j = 1, 4
				!	b2 = b1 - CUT(:, j, n)
				!	CUT(:, j, n) = CUT(:, j, n) + b2 * 0.0001 / norm2(b2)
				!end do
				
				
				n = n + 1
				!print*, i
				!pause
			end if
			if(n > 500000) EXIT
		end do
			

		open(1, file = name // '_print_par_2D_interpolate.txt')
		write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z', 'rho', 'u', 'v', 'w', 'p',"
		write(1,*) "'bx', 'by', 'bz', 'Q', 'n_H', 'u_H',"
		write(1,*) "'v_H', 'w_H', 'T_H', 'MK6', 'MK7', 'MK8', 'MK9'"
		write(1,*) ",'n_H1', u_H1', 'v_H1', 'w_H1', 'p_T_H1', 'MK16', 'MK17', 'MK18', 'MK19'"
		write(1,*) ",'n_H2', u_H2', 'v_H2', 'w_H2', 'p_T_H2', 'MK26', 'MK27', 'MK28', 'MK29'"
		write(1,*) ",'n_H3', u_H3', 'v_H3', 'w_H3', 'p_T_H3', 'MK36', 'MK37', 'MK38', 'MK39'"
		write(1,*) ",'n_H4', u_H4', 'v_H4', 'w_H4', 'p_T_H4', 'MK46', 'MK47', 'MK48', 'MK49'"
		write(1,*) ",'k2', k3', 'k4', 'k5'"
		write(1,*) ", ZONE T= 'HP', N= ",  4 * (n - 1) , ", E =  ", (n - 1), ", F=FEPOINT, ET=quadrilateral "

		do i = 1, (n - 1)
			do j = 1, 4
				call Int2_Get_par_fast(CUT(1, j, i), CUT(2, j, i), CUT(3, j, i), num, PAR, PAR_MK, PAR_K)
				call Int2_sum_moment(PAR_MK, PAR_MK_SUM)
				if(num < 1) num = 3
				write(1,*) CUT(:, j, i), PAR, PAR_MK_SUM, PAR_MK(1:9, 1), PAR_MK(1:9, 2), PAR_MK(1:9, 3), PAR_MK(1:9, 4), PAR_K(2:5)
			end do
		end do

		do i = 1, (n - 1)
			write(1,*) 4 * (i - 1) + 1, 4 * (i - 1) + 2, 4 * (i - 1) + 3, 4 * (i - 1) + 4
		end do
		
		deallocate(CUT)
		close(1)


	end subroutine Int_2_Print_par_2D
	
	subroutine Int2_sum_moment(PAR_MK, PAR_MK_SUM)
	! Правильно суммирует моменты, чтобы найти общие параметры
		real(8), intent(in) :: PAR_MK(par_n_moment, par_n_sort)
		real(8), intent(out) :: PAR_MK_SUM(par_n_moment)
		
		integer(4) :: i
		
		PAR_MK_SUM = 0.0
	
		PAR_MK_SUM(1) = SUM(PAR_MK(1, :))
		
		PAR_MK_SUM(6:9) = SUM(PAR_MK(6:9, :), dim = 2)
	
		do i = 1, par_n_sort
			PAR_MK_SUM(2:4) = PAR_MK_SUM(2:4) + PAR_MK(2:4, i) * PAR_MK(1, i)
			PAR_MK_SUM(5) = PAR_MK_SUM(5) + ( PAR_MK(5, i) / (2.0/3.0) + &
				kvv(PAR_MK(2, i), PAR_MK(3, i), PAR_MK(4, i)) ) * PAR_MK(1, i)
			
		end do
		
		if(PAR_MK_SUM(1) > 0.0000001) then
			PAR_MK_SUM(2:4) = PAR_MK_SUM(2:4) / PAR_MK_SUM(1)
			PAR_MK_SUM(5) = (2.0/3.0) * (PAR_MK_SUM(5)/PAR_MK_SUM(1) - kvv(PAR_MK_SUM(2), PAR_MK_SUM(3), PAR_MK_SUM(4)))
		end if
		
	
	end subroutine Int2_sum_moment
	
	subroutine Int2_Print_my()
		implicit none
		integer(4) :: i, j, k, N1, N2, N3, s1, ii
		
		open(1, file = 'Int2_print_my.txt')
		
		N1 = size(int2_gran_point(1, :))
	
		do i = 1, N1
			if(int2_gran_sosed(i) /= 0) CYCLE !Есть сосед
			
			if(int2_gran_point(1, i) == 0) CYCLE  ! Нет грани
			
			!if(norm2(int2_coord(:, int2_gran_point(1, i))) < 25) then
			!	print*, int2_coord(:, int2_gran_point(1, i))
			!	print*, int2_coord(:, int2_gran_point(2, i))
			!	print*, int2_coord(:, int2_gran_point(3, i))
			!	pause
			!end if
			
			
			do k = 1, 3 ! Бежим по граням тетраэдра
				write(1,*) int2_coord(:, int2_gran_point(1, i))
				write(1,*) int2_coord(:, int2_gran_point(2, i))
				write(1,*) int2_coord(:, int2_gran_point(3, i))
			end do
		end do
		
		
	end subroutine Int2_Print_my
	
	
	subroutine Int2_Re_interpol()
	! Меняет значения переменных (основных) основной сетки из интерполированной
	! Variables
	integer(4) :: N1, N2, N3, i, num
	real(8) :: r(3), PAR(9)
	
	if( size(gl_Cell_par(:, 1)) /= size(int2_Cell_par(:, 1)) .or. size(gl_Cell_par(:, 1)) /= 9) then
		print*, "Error ishflspeurbvmdr  ", size(gl_Cell_par(:, 1)), size(int2_Cell_par(:, 1))
		STOP
	end if
	
	N1 = size(gl_Cell_center(1, :))
	num = 3
	
	do i = 1, N1
		r = gl_Cell_center(:, i)
		call Int2_Get_par_fast(r(1), r(2), r(3), num, PAR)
		if(num < 1) then
			num = 3
			call Int2_Get_tetraedron_inner(r(1), r(2), r(3), num)
			if(num < 1) STOP "ERROR  u;dkfg783hju"
			gl_Cell_par(:, i) = int2_Cell_par(:, int2_all_tetraendron_point(1, num))
		else
			gl_Cell_par(:, i) = PAR
		end if
	end do
	
	
	
	! Body of Int2_Re_interpol
	
	end subroutine Int2_Re_interpol
	
	
	subroutine Int2_Print_Cell(n)  ! Печатает ячейку с номером n
    use GEO_PARAM
    use STORAGE

    implicit none
    integer, intent(in) :: n
    integer(4) :: i, k

    open(1, file = 'one_Cell.txt')


    write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z'  ZONE T= 'HP', N= ", 8, ", E =  ", 12 , ", F=FEPOINT, ET=LINESEG "

    do i = 1, 8
            write(1,*) int2_coord(:, int2_all_Cell(i, n))
    end do

    write(1,*) "1 2"
    write(1,*) "2 3"
    write(1,*) "3 4"
    write(1,*) "4 1"
    write(1,*) "1 5"
    write(1,*) "2 6"
    write(1,*) "3 7"
    write(1,*) "4 8"
    write(1,*) "5 6"
    write(1,*) "6 7"
    write(1,*) "7 8"
    write(1,*) "8 5"

    close(1)

    end subroutine Int2_Print_Cell
	
	subroutine Int2_Print_tetraedron(n)
		implicit none
		integer(4), intent(in) :: n
		integer(4) :: i, j, k, N1, N2, N3, s1, ii
		character(len=8) :: name
    
		write(unit=name,fmt='(i8.8)') n
		
		open(1, file = "Int2_print_tetraedron_" // name // ".txt")
		write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z', 'I'  ZONE T= 'HP', N= ", 12, ", E =  ", 12 , ", F=FEPOINT, ET=LINESEG "
	
		i = n
			
		do k = 1, 4 ! Бежим по граням тетраэдра
			s1 = int2_all_tetraendron(k, i)  ! Номер грани
			if (s1 == 0) then
				print*, "NO tetraendron  ", n
				CYCLE
			end if
			
			write(1,*) int2_coord(:, int2_gran_point(1, s1)), n
			write(1,*) int2_coord(:, int2_gran_point(2, s1)), n
			write(1,*) int2_coord(:, int2_gran_point(3, s1)), n
			
		end do
		
			! Connectivity list
		do j = 0, 3
			write(1,*) 3 * j + 1, 3 * j + 2
			write(1,*) 3 * j + 2, 3 * j + 3
			write(1,*) 3 * j + 3, 3 * j + 1
		end do
		
		
	end subroutine Int2_Print_tetraedron
	
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
	
	
	subroutine Int2_Print_sosed()
	! печать двойственной сетки
		implicit none
		integer(4) :: i, j, k, N1, N2, N3, N, l, s1, s2
		
		N = size(int2_Cell_B(1, :, 1)) * size(int2_Cell_B(:, 1, 1)) + size(int2_Cell_A(1, :, 1)) * size(int2_Cell_A(:, 1, 1)) + &
			size(int2_Cell_C(1, :, 1)) * size(int2_Cell_C(:, 1, 1))
		
		open(1, file = 'Int2_sosed.txt')
		write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z', 'Q'  ZONE T= 'HP', N= ", 8 * N, ", E =  ", 4 * N , ", F=FEPOINT, ET=LINESEG "
		
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
						write(1,*) 0.0, 0.0, 0.0, 1
						write(1,*) 0.0, 0.0, 0.0, 1
						write(1,*) 0.0, 0.0, 0.0, 1
						write(1,*) 0.0, 0.0, 0.0, 1
					else
						s1 = int2_Cell_A(i, j, k)
						do l = 1, 4
							s2 = int2_all_neighbours(l, s1)
							if (s2 /= 0) then
								write(1,*) int2_Cell_center(:, s1), 1
								write(1,*) int2_Cell_center(:, s2), 1
							else
								write(1,*) 0.0, 0.0, 0.0, 1
								write(1,*) 0.0, 0.0, 0.0, 1
							end if
							
						end do
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
						write(1,*) 0.0, 0.0, 0.0, 1
						write(1,*) 0.0, 0.0, 0.0, 1
						write(1,*) 0.0, 0.0, 0.0, 1
						write(1,*) 0.0, 0.0, 0.0, 1
						write(1,*) 0.0, 0.0, 0.0, 1
						write(1,*) 0.0, 0.0, 0.0, 1
						write(1,*) 0.0, 0.0, 0.0, 1
						write(1,*) 0.0, 0.0, 0.0, 1
					else
						s1 = int2_Cell_B(i, j, k)
						do l = 1, 4
							s2 = int2_all_neighbours(l, s1)
							if (s2 /= 0) then
								write(1,*) int2_Cell_center(:, s1), 2
								write(1,*) int2_Cell_center(:, s2), 2
							else
								write(1,*) 0.0, 0.0, 0.0, 1
								write(1,*) 0.0, 0.0, 0.0, 1
							end if
						end do
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
						write(1,*) 0.0, 0.0, 0.0, 1
						write(1,*) 0.0, 0.0, 0.0, 1
						write(1,*) 0.0, 0.0, 0.0, 1
						write(1,*) 0.0, 0.0, 0.0, 1
						write(1,*) 0.0, 0.0, 0.0, 1
						write(1,*) 0.0, 0.0, 0.0, 1
						write(1,*) 0.0, 0.0, 0.0, 1
						write(1,*) 0.0, 0.0, 0.0, 1
					else
						s1 = int2_Cell_C(i, j, k)
						do l = 1, 4
							s2 = int2_all_neighbours(l, s1)
							if (s2 /= 0) then
								write(1,*) int2_Cell_center(:, s1), 3
								write(1,*) int2_Cell_center(:, s2), 3
							else
								write(1,*) 0.0, 0.0, 0.0, 3
								write(1,*) 0.0, 0.0, 0.0, 3
							end if
						end do
					end if
				end do
			end do
		end do
		
		! Connectivity list
    do j = 0, N
        write(1,*) 8 * j + 1, 8 * j + 2
        write(1,*) 8 * j + 3, 8 * j + 4
        write(1,*) 8 * j + 5, 8 * j + 6
        write(1,*) 8 * j + 7, 8 * j + 8
    end do
		
		
		close(1)
	
	end subroutine Int2_Print_sosed

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
	USE GEO_PARAM
	implicit none
	integer :: n
	! Выделения памяти и начальная инициализация
	allocate(int2_Point_A( size(gl_Cell_A(:, 1, 1)) + 5 + 1, size(gl_Cell_A(1, :, 1)) + 1, size(gl_Cell_A(1, 1, :)) ))
	allocate(int2_Point_B( size(gl_Cell_B(:, 1, 1)) + 3 + 1, size(gl_Cell_B(1, :, 1)) + 1, size(gl_Cell_B(1, 1, :)) ))
	allocate( int2_Point_C( size(gl_Cell_C(:, 1, 1)) + 2 + 1, size(gl_Cell_C(1, :, 1)) + 1, size(gl_Cell_C(1, 1, :)) ) )
	
	n = size(gl_Cell_A(:, 1, 1)) + 2 + 1
	n = n * (size(gl_Cell_A(1, :, 1)) + 1) * size(gl_Cell_A(1, 1, :))
	n = n + (size(gl_Cell_B(:, 1, 1)) + 1 + 1) * size(gl_Cell_B(1, :, 1)) * size(gl_Cell_B(1, 1, :))
	n = n + (size(gl_Cell_C(:, 1, 1)) + 1 + 1) * (size(gl_Cell_C(1, :, 1)) - 1 + 1) * size(gl_Cell_C(1, 1, :))
	n = n + size(int2_Point_A(1, 1, :))
	
	allocate(int2_all_Cell(8, n ))  ! Добавили ячейку около тройной точки
	allocate(int2_coord(3, size(int2_Point_A) + size(int2_Point_B) + size(int2_Point_C)))
	allocate(int2_Cell_center(3, n))
	
	allocate(int2_Cell_A(size(gl_Cell_A(:, 1, 1)) + 4 + 1, size(gl_Cell_A(1, :, 1)) + 1, size(gl_Cell_A(1, 1, :))))
	allocate(int2_Cell_B(size(gl_Cell_B(:, 1, 1)) + 2 + 1, size(gl_Cell_B(1, :, 1)), size(gl_Cell_B(1, 1, :))))
	allocate(int2_Cell_C(size(gl_Cell_C(:, 1, 1)) + 2 + 1, size(gl_Cell_C(1, :, 1)) - 1 + 1, size(gl_Cell_C(1, 1, :))))
	
	allocate(int2_all_neighbours(6, size(int2_all_Cell(1, :)) ))
	
	allocate(int2_all_tetraendron(4, 6 * n))
	allocate(int2_all_Volume(6 * n))
	allocate(int2_all_tetraendron_point(4, 6 * n))
	allocate(int2_all_tetraendron_matrix(4, 4, 6 * n))
	allocate(int2_gran_point(3, 6 * n * 4))
	allocate(int2_gran_sosed(6 * n * 4))
	!allocate(int2_gran_normal(3, 6 * n * 4))
	allocate(int2_plane_tetraendron(4, 4, 6 * n))
	
	allocate(int2_Cell_par( size(gl_Cell_par(:, 1)), size(int2_Point_A) + size(int2_Point_B) + size(int2_Point_C) ))
	allocate(int2_Cell_par_2( size(gl_Cell_par2(:, 1)), size(int2_Point_A) + size(int2_Point_B) + size(int2_Point_C) ))
	allocate(int2_Cell_par2( 1, size(int2_Point_A) + size(int2_Point_B) + size(int2_Point_C) ))
	allocate(int2_Moment(par_n_moment, par_n_sort, size(int2_Point_A) + size(int2_Point_B) + size(int2_Point_C)))
	allocate(int2_Moment_k(5, size(int2_Point_A) + size(int2_Point_B) + size(int2_Point_C)))
	
	
	int2_Moment = 0.0
	int2_Moment_k = 1.0
	int2_Moment_k(1, :) = 0.0
	int2_all_Volume = 0.0
	int2_plane_tetraendron = 0.0
	int2_coord = 0.0
	int2_all_Cell = -100
	int2_Point_A = -100
	int2_Point_B = -100
	int2_Cell_center(:, :) = 0.0
	int2_Cell_A = 0
	int2_Cell_B = 0
	int2_Cell_C = 0
	int2_all_neighbours = 0 
	int2_all_tetraendron = 0
	int2_all_tetraendron_point = 0
	int2_all_tetraendron_matrix = 0.0
	int2_gran_point = 0
	int2_gran_sosed = 0
	!int2_gran_normal = 0.0
	int2_Cell_par = 0.0
	int2_Cell_par2 = 1
	int2_Cell_par_2 = 0.0
	
	end subroutine Int2_Set_Interpolate
	
	subroutine Int2_Dell_interpolate()
	! Удаляет все массивы интерполяции
	deallocate(int2_Point_A)
	deallocate(int2_Point_B)
	deallocate(int2_Point_C)
	deallocate(int2_all_Cell)
	deallocate(int2_coord)
	deallocate(int2_Cell_center)
	deallocate(int2_Cell_A)
	deallocate(int2_Cell_B)
	deallocate(int2_Cell_C)
	deallocate(int2_all_neighbours)
	deallocate(int2_all_tetraendron)
	deallocate(int2_all_Volume)
	deallocate(int2_all_tetraendron_point)
	deallocate(int2_all_tetraendron_matrix)
	deallocate(int2_gran_point)
	deallocate(int2_gran_sosed)
	deallocate(int2_plane_tetraendron)
	deallocate(int2_Cell_par)
	deallocate(int2_Cell_par2)
	deallocate(int2_Cell_par_2)
	deallocate(int2_Moment)
	deallocate(int2_Moment_k)
	
	! Body of Int2_Dell_interpolate
	
	end subroutine Int2_Dell_interpolate
	
	
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
! Файл для ТВД процедуры
	
	subroutine Find_TVD_sosed()
	use STORAGE
    use GEO_PARAM
	implicit none
	
	integer(4) :: gr, Num, cell, cell2, j, best, cell1, cell3, cell4
	real(8) :: normal(3), c1(3), c2(3), vec(3)
	real(8) :: dotmax, xc
	real(8) :: d1, d2, d3, d4
	real(8) :: c3(3), c4(3), gg(3)
	! Body of Find_TVD_sosed
	
	Num = size(gl_Gran_neighbour(1, :))   ! Число граней в сетке
	
	do gr = 1, Num
		
		normal = gl_Gran_normal(:, gr)
		cell = 0
		cell2 = 0
		dotmax = 0
		best = 0
		
		if(gl_Gran_neighbour(1, gr) > 0) then
			cell = gl_Gran_neighbour(1, gr)
			
			! Не надо делать твд для ячеек около нуля
			if ((gl_Cell_type(cell) == "A" .or. gl_Cell_type(cell) == "B").and.(gl_Cell_number(1, cell) <= 1) ) then
				gl_Gran_neighbour_TVD(1, gr) = 0
				GO TO 11 
			end if
			
			c1 = gl_Cell_center(:, cell)
			do j = 1, 6
				if(gl_Cell_neighbour(j, cell) <= 0) CYCLE
				cell2 = gl_Cell_neighbour(j, cell)
				c2 = gl_Cell_center(:, cell2)
				vec = c2 - c1
				vec = vec/norm2(vec)
				xc = DOT_PRODUCT(vec, -normal)
				if (dotmax < xc) then
					dotmax = xc
					best = j
				end if
			end do
			if (best == 0) then
				gl_Gran_neighbour_TVD(1, gr) = 0
			else
				gl_Gran_neighbour_TVD(1, gr) = gl_Cell_neighbour(best, cell)
			end if
		else
			gl_Gran_neighbour_TVD(1, gr) = 0
		end if
		
11		CONTINUE 
		
		cell = 0
		cell2 = 0
		dotmax = 0
		best = 0
		
		if(gl_Gran_neighbour(2, gr) > 0) then
			cell = gl_Gran_neighbour(2, gr)
			
			if ((gl_Cell_type(cell) == "A" .or. gl_Cell_type(cell) == "B").and.(gl_Cell_number(1, cell) <= 1) ) then
				gl_Gran_neighbour_TVD(2, gr) = 0
				CYCLE
			end if
			
			
			c1 = gl_Cell_center(:, cell)
			do j = 1, 6
				if(gl_Cell_neighbour(j, cell) <= 0) CYCLE
				cell2 = gl_Cell_neighbour(j,cell)
				c2 = gl_Cell_center(:, cell2)
				vec = c2 - c1
				vec = vec/norm2(vec)
				xc = DOT_PRODUCT(vec, normal)
				if (dotmax < xc) then
					dotmax = xc
					best = j
				end if
			end do
			if (best == 0) then
				gl_Gran_neighbour_TVD(2, gr) = 0
			else
				gl_Gran_neighbour_TVD(2, gr) = gl_Cell_neighbour(best, cell)
			end if
		else
			gl_Gran_neighbour_TVD(2, gr) = 0
		end if
		
		CONTINUE 
		
	end do
	
	
	! Можно сделать некоторые проверки правильности расчёта соседей
	do gr = 1, Num
		if(gl_Gran_neighbour(1, gr) == 0) CYCLE
		if(gl_Gran_neighbour(2, gr) == 0) CYCLE
		if(gl_Gran_neighbour_TVD(1, gr) == 0) CYCLE
		if(gl_Gran_neighbour_TVD(2, gr) == 0) CYCLE
		
		cell1 = gl_Gran_neighbour(1, gr)
		cell2 = gl_Gran_neighbour(2, gr)
		cell3 = gl_Gran_neighbour_TVD(1, gr)
		cell4 = gl_Gran_neighbour_TVD(2, gr)
		
		c1 = gl_Cell_center(:, cell1)
		c2 = gl_Cell_center(:, cell2)
		c3 = gl_Cell_center(:, cell3)
		c4 = gl_Cell_center(:, cell4)
		
		if(DOT_PRODUCT(c2-c1, c4-c3) < 0.0) then
			print*, gl_Gran_normal(:, gr)
			print*, "Problem  93 TVD.f90  ", gr
			PAUSE
		end if
		
	end do
	
	print*, "TVD procedure proverena, problem ne obnarugeno"
	
	end subroutine Find_TVD_sosed
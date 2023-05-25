
    
    subroutine Find_Surface()   ! Заполняет массивы поверхностей, которые выделяются
    ! Для BS находятся только поверхности на полусфере! (x > 0)
    use STORAGE
    use GEO_PARAM
    implicit none
    integer :: j, k, node, num, i

    ! HP
    node = 1

    do k = 1, size( gl_Cell_A(par_n_HP - 1, 1, :) )
        do j = 1, size( gl_Cell_A(par_n_HP - 1, :, 1) )
            gl_Contact(node) = gl_Cell_gran(1, gl_Cell_A(par_n_HP - 1, j, k))
			gl_Gran_type(gl_Contact(node)) = 2
            node = node + 1
			
			! Вводим лакса поперёк разрыва 
			gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_A(par_n_HP - 1, j, k))) = 1
			if(gl_Cell_gran(4, gl_Cell_A(par_n_HP - 1, j, k)) /= 0) then
			    gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_A(par_n_HP - 1, j, k)) ) = 1
			end if
			gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_A(par_n_HP - 1, j, k))) = 1
			gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_A(par_n_HP - 1, j, k))) = 1
			
			gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_A(par_n_HP - 2, j, k))) = 1
			if(gl_Cell_gran(4, gl_Cell_A(par_n_HP - 2, j, k)) /= 0) then
			    gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_A(par_n_HP - 2, j, k)) ) = 1
			end if
			gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_A(par_n_HP - 2, j, k))) = 1
			gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_A(par_n_HP - 2, j, k))) = 1
			
			
			gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_A(par_n_HP, j, k))) = 1
			if(gl_Cell_gran(4, gl_Cell_A(par_n_HP, j, k)) /= 0) then
			    gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_A(par_n_HP, j, k))) = 1
			end if
			gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_A(par_n_HP, j, k))) = 1
			gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_A(par_n_HP, j, k))) = 1
			
			gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_A(par_n_HP + 1, j, k))) = 1
			if(gl_Cell_gran(4, gl_Cell_A(par_n_HP + 1, j, k)) /= 0) then
			    gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_A(par_n_HP + 1, j, k))) = 1
			end if
			gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_A(par_n_HP + 1, j, k))) = 1
			gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_A(par_n_HP + 1, j, k))) = 1
			
        end do
    end do

    do k = 1, size( gl_Cell_C(par_n_HP - par_n_TS, 1, :) )
        do j = 1, size( gl_Cell_C(par_n_HP - par_n_TS, :, 1) )
            gl_Contact(node) = gl_Cell_gran(1, gl_Cell_C(par_n_HP - par_n_TS, j, k))
			gl_Gran_type(gl_Contact(node)) = 2
            node = node + 1
			
			! Вводим лакса поперёк разрыва 
			gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_C(par_n_HP - par_n_TS, j, k))) = 1
			if(gl_Cell_gran(4, gl_Cell_C(par_n_HP - par_n_TS, j, k)) /= 0) then
			    gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_C(par_n_HP - par_n_TS, j, k)) ) = 1
			end if
			gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_C(par_n_HP - par_n_TS, j, k))) = 1
			gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_C(par_n_HP - par_n_TS, j, k))) = 1
			
			gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_C(par_n_HP - par_n_TS - 1, j, k))) = 1
			if(gl_Cell_gran(4, gl_Cell_C(par_n_HP - par_n_TS - 1, j, k)) /= 0) then
			    gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_C(par_n_HP - par_n_TS - 1, j, k)) ) = 1
			end if
			gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_C(par_n_HP - par_n_TS - 1, j, k))) = 1
			gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_C(par_n_HP - par_n_TS - 1, j, k))) = 1
			
			
			gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_C(par_n_HP - par_n_TS + 1, j, k))) = 1
			if(gl_Cell_gran(4, gl_Cell_C(par_n_HP - par_n_TS + 1, j, k)) /= 0) then
			    gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_C(par_n_HP - par_n_TS + 1, j, k))) = 1
			end if
			gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_C(par_n_HP - par_n_TS + 1, j, k))) = 1
			gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_C(par_n_HP - par_n_TS + 1, j, k))) = 1
			
			gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_C(par_n_HP - par_n_TS + 2, j, k))) = 1
			if(gl_Cell_gran(4, gl_Cell_C(par_n_HP - par_n_TS + 2, j, k)) /= 0) then
			    gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_C(par_n_HP - par_n_TS + 2, j, k))) = 1
			end if
			gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_C(par_n_HP - par_n_TS + 2, j, k))) = 1
			gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_C(par_n_HP - par_n_TS + 2, j, k))) = 1
        end do
    end do

    node = 1
    ! TS
    do k = 1, size( gl_Cell_A(par_n_TS - 1, 1, :) )
        do j = 1, size( gl_Cell_A(par_n_TS - 1, :, 1) )
            gl_TS(node) = gl_Cell_gran(1, gl_Cell_A(par_n_TS - 1, j, k))
			gl_Gran_type(gl_TS(node)) = 1
            node = node + 1
			
			! Вводим лакса поперёк разрыва 
			gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_A(par_n_TS - 1, j, k))) = 1
			if(gl_Cell_gran(4, gl_Cell_A(par_n_TS - 1, j, k)) /= 0) then
			    gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_A(par_n_TS - 1, j, k)) ) = 1
			end if
			gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_A(par_n_TS - 1, j, k))) = 1
			gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_A(par_n_TS - 1, j, k))) = 1
			!
			!gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_A(par_n_TS, j, k))) = 1
			!if(gl_Cell_gran(4, gl_Cell_A(par_n_TS, j, k)) /= 0) then
			!    gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_A(par_n_TS, j, k))) = 1
			!end if
			!gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_A(par_n_TS, j, k))) = 1
			!gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_A(par_n_TS, j, k))) = 1
			
			gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_A(par_n_TS - 2, j, k))) = 1
			if(gl_Cell_gran(4, gl_Cell_A(par_n_TS - 2, j, k)) /= 0) then
			    gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_A(par_n_TS - 2, j, k)) ) = 1
			end if
			gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_A(par_n_TS - 2, j, k))) = 1
			gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_A(par_n_TS - 2, j, k))) = 1
			
        end do
    end do

    do k = 1, size( gl_Cell_B(par_n_TS - 1, 1, :) )
        do j = 1, size( gl_Cell_B(par_n_TS - 1, :, 1) )
            gl_TS(node) = gl_Cell_gran(1, gl_Cell_B(par_n_TS - 1, j, k))
			gl_Gran_type(gl_TS(node)) = 1
            node = node + 1
			
			! Вводим лакса поперёк разрыва 
			!gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_B(par_n_TS - 1, j, k))) = 1
			!if(gl_Cell_gran(4, gl_Cell_B(par_n_TS - 1, j, k)) /= 0) then
			!    gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_B(par_n_TS - 1, j, k)) ) = 1
			!end if
			!gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_B(par_n_TS - 1, j, k))) = 1
			!gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_B(par_n_TS - 1, j, k))) = 1
			!
			!gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_B(par_n_TS, j, k))) = 1
			!if(gl_Cell_gran(4, gl_Cell_B(par_n_TS, j, k)) /= 0) then
			!    gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_B(par_n_TS, j, k))) = 1
			!end if
			!gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_B(par_n_TS, j, k))) = 1
			!gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_B(par_n_TS, j, k))) = 1
			
        end do
    end do


    node = 1
    ! BS
    num = par_m_A - 1
    do k = 1, size( gl_Cell_A(par_n_BS - 1, 1, :) )
        do j = 1, num
            gl_BS(node) = gl_Cell_gran(1, gl_Cell_A(par_n_BS - 1, j, k))
            node = node + 1
        end do
    end do

    if (par_developer_info) print *, "Find_Surface: Contact", node - 1, size(gl_Contact)
    if (par_developer_info) print *, "Find_Surface: TS", node - 1, size(gl_TS)
    if (par_developer_info) print *, "Find_Surface: BS", node - 1, size(gl_BS)
	
	! Вводим Лакса на оси симметрии
	
	do k = 1, size( gl_Cell_A(1, 1, :) )
        do j = 1, size( gl_Cell_A(:, 1, 1) )
			! Вводим лакса поперёк
			gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_A(j, 1, k))) = 1
			if(gl_Cell_gran(4, gl_Cell_A(j, 1, k)) /= 0) then
			    gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_A(j, 1, k))) = 1
			end if
			gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_A(j, 1, k))) = 1
			gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_A(j, 1, k))) = 1
			
			gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_A(j, 2, k))) = 1
			if(gl_Cell_gran(4, gl_Cell_A(j, 2, k)) /= 0) then
			    gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_A(j, 2, k))) = 1
			end if
			gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_A(j, 2, k))) = 1
			gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_A(j, 2, k))) = 1
			
        end do
	end do
	
	do k = 1, size( gl_Cell_B(1, 1, :) )
        do j = 1, par_n_TS + 5 !size( gl_Cell_B(:, 1, 1) )
			! Вводим лакса поперёк
			gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_B(j, 1, k))) = 1
			if (gl_Cell_gran(4, gl_Cell_B(j, 1, k)) /= 0) then
			    gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_B(j, 1, k))) = 1
			end if
			gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_B(j, 1, k))) = 1
			gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_B(j, 1, k))) = 1
			
			gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_B(j, 2, k))) = 1
			if (gl_Cell_gran(4, gl_Cell_B(j, 2, k)) /= 0) then
			    gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_B(j, 2, k))) = 1
			end if
			gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_B(j, 2, k))) = 1
			gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_B(j, 2, k))) = 1
			
        end do
	end do
	
	
	do k = 1, size( gl_Cell_A(par_n_HP - 1, 1, :) )
        do j = 1, size( gl_Cell_A(par_n_HP - 1, :, 1) )
			do i = 1, par_n_TS - 2
				gl_Gran_scheme(gl_Cell_gran(1, gl_Cell_A(i, j, k))) = 3
			end do
		end do
	end do
			
	

    end subroutine Find_Surface
    
    subroutine Calc_move(now) ! Вычисляем скорости движения узлов (ключевых) с помощью распада на поверхностях разрыва(
    use STORAGE
    use GEO_PARAM
	use Solvers
    implicit none
    
    integer, intent(in) :: now   ! Какиой номер параметров используется в вычислении движения
    
    integer :: i, Num, gr, s1, s2, j, yzel, ndisc
    real(8) :: normal(3), qqq1(8), qqq2(8), dsl, dsp, dsc, POTOK(8), a1, a2, v1, v2, &
        ray(3), norm, b1, b2, c1, c2
    
    !koef1 = 0.3! 0.2
    !koef2 = 1.0 ! 1.0
    !koef3 = 0.7   ! 0.3
	
	ndisc = 0
    
    ! Пробегаемся по всем граням и вычисляем скорости их движения
    
    ! TS
    Num = size(gl_TS)
    
    do i = 1, Num
        
        gr = gl_TS(i)
    	s1 = gl_Gran_neighbour(1, gr)
        s2 = gl_Gran_neighbour(2, gr)
        normal = gl_Gran_normal2(:, gr, now)
        qqq1 = gl_Cell_par(1:8, s1)
        qqq2 = gl_Cell_par(1:8, s2)
        
        call chlld(2, normal(1), normal(2), normal(3), &  ! 2
                0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)
		
        
        dsl = dsl * koef1
		
        !if(i == 50) write(*,*) dsl
		
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
                
                ray = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/) - gl_Gran_center2(:, gr, now)
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
        
    ! Контакт
    Num = size(gl_Contact)
    
    do i = 1, Num
        
        gr = gl_Contact(i)
    	s1 = gl_Gran_neighbour(1, gr)
        s2 = gl_Gran_neighbour(2, gr)
        normal = gl_Gran_normal2(:, gr, now)
        qqq1 = gl_Cell_par(1:8, s1)
        qqq2 = gl_Cell_par(1:8, s2)
		
		! Вычтем нормальную компоненту магнитного поля
		if (gl_Gran_center2(1, gr, now) > -100.0) then
	        qqq1(6:8) = qqq1(6:8) - DOT_PRODUCT(normal, qqq1(6:8)) * normal
	        qqq2(6:8) = qqq2(6:8) - DOT_PRODUCT(normal, qqq2(6:8)) * normal
	        gl_Cell_par(6:8, s1) = qqq1(6:8)
	        gl_Cell_par(6:8, s2) = qqq2(6:8)
		end if
		
		
        
        call chlld(3, normal(1), normal(2), normal(3), &  ! 2
                0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)
        
        !dsc = dsc * koef2
		dsc = (dsc + 1.0 * DOT_PRODUCT(0.5 * (qqq1(2:4) + qqq2(2:4)), normal)) * koef2
		
        !if(i == 2) write(*,*) dsc
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
                
                ray = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/) - gl_Gran_center2(:, gr, now)
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
	if(.True.) then
    Num = size(gl_BS)
    
        do i = 1, Num
        
            gr = gl_BS(i)
    	    s1 = gl_Gran_neighbour(1, gr)
            s2 = gl_Gran_neighbour(2, gr)
            normal = gl_Gran_normal2(:, gr, now)
            qqq1 = gl_Cell_par(1:8, s1)
            qqq2 = gl_Cell_par(1:8, s2)
        
            call chlld(2, normal(1), normal(2), normal(3), &  ! 2
                    0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)
        
            dsp = dsp * koef3
			!if(i == 2) then
		 !       write(*,*) dsp
		 !       write(*,*) normal(1), normal(2), normal(3)
		 !       write(*,*) 2
		 !       write(*,*) qqq1(1), qqq1(2), qqq1(3), qqq1(4), qqq1(5), qqq1(6), qqq1(7), qqq1(8)
		 !       write(*,*) qqq2(1), qqq2(2), qqq2(3), qqq2(4), qqq2(5), qqq2(6), qqq2(7), qqq2(8)
	  !      end if
        
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
                
                    ray = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/) - gl_Gran_center2(:, gr, now)
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
	end if
	
    end subroutine Calc_move   
    
    subroutine Move_all(now, Time)  ! Передвигаем точки в соответствии со скоростью движение узлов на поверхностях разрыва
    use STORAGE
    use GEO_PARAM
    implicit none
    
    integer :: now2             ! Эти параметры мы сейчас меняем на основе now
    real(8), intent(in) :: Time
    integer, intent(in) :: now
    real(8) :: R_TS, proect, vel(3), R_HP, R_BS, Ak(3), Bk(3), Ck(3), Dk(3), Ek(3), ER(3), KORD(3)
    integer :: yzel, N1, N2, N3, i, j, k, yzel2
    real(8) :: the, phi, r, x, y, z, rr, xx, x2, y2, z2, rrr
    
    now2 = mod(now, 2) + 1   
	
	! Поверхностное натяжение ------------------------------------------------------------------------------
	if (.True.) then
	N3 = size(gl_RAY_A(1, 1, :))
    N2 = size(gl_RAY_A(1, :, 1))
    N1 = size(gl_RAY_A(:, 1, 1))

    ! Цикл движения точек на лучах А  
    do k = 1, N3
        do j = 1, N2
			
			if (k /= 1 .and. j == 1) then
                    CYCLE
			end if
			
			! Ударная волна
			
			yzel = gl_RAY_A(par_n_TS, j, k)
			Ak = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/)
			
			if (j < N2) then
			    yzel2 = gl_RAY_A(par_n_TS, j + 1, k)
			else
				yzel2 = gl_RAY_B(par_n_TS, 1, k)
			end if
			Bk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if (k > 1) then
			    yzel2 = gl_RAY_A(par_n_TS, j, k - 1)
			else
			    yzel2 = gl_RAY_A(par_n_TS, j, N3)
			end if
			Ck = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if (j > 1) then
			    yzel2 = gl_RAY_A(par_n_TS, j - 1, k)
			else
				yzel2 = gl_RAY_A(par_n_TS, j, k + ceiling(1.0 * N3/2))
			end if
			Dk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if(k < N3) then
			    yzel2 = gl_RAY_A(par_n_TS, j, k + 1)
			else
				yzel2 = gl_RAY_A(par_n_TS, j, 1)
			end if
			Ek = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
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
			Ak = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/)
			
			if (j < N2) then
			    yzel2 = gl_RAY_A(par_n_HP, j + 1, k)
			else
				yzel2 = gl_RAY_B(par_n_HP, 1, k)
			end if
			Bk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if (k > 1) then
			    yzel2 = gl_RAY_A(par_n_HP, j, k - 1)
			else
			    yzel2 = gl_RAY_A(par_n_HP, j, N3)
			end if
			Ck = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if (j > 1) then
			    yzel2 = gl_RAY_A(par_n_HP, j - 1, k)
			else
				yzel2 = gl_RAY_A(par_n_HP, j, k + ceiling(1.0 * N3/2))
			end if
			Dk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if(k < N3) then
			    yzel2 = gl_RAY_A(par_n_HP, j, k + 1)
			else
				yzel2 = gl_RAY_A(par_n_HP, j, 1)
			end if
			Ek = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
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
			Ak = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/)
			
			if (j < N2) then
			    yzel2 = gl_RAY_A(par_n_BS, j + 1, k)
			else
				CYCLE !yzel2 = yzel                                                    ! Вот тут есть неопределённость что делать с крайним узлом
			end if
			Bk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if (k > 1) then
			    yzel2 = gl_RAY_A(par_n_BS, j, k - 1)
			else
			    yzel2 = gl_RAY_A(par_n_BS, j, N3)
			end if
			Ck = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if (j > 1) then
			    yzel2 = gl_RAY_A(par_n_BS, j - 1, k)
			else
				yzel2 = gl_RAY_A(par_n_BS, j, k + ceiling(1.0 * N3/2))
			end if
			Dk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if(k < N3) then
			    yzel2 = gl_RAY_A(par_n_BS, j, k + 1)
			else
				yzel2 = gl_RAY_A(par_n_BS, j, 1)
			end if
			Ek = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
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
	
	
	N3 = size(gl_RAY_B(1, 1, :))
    N2 = size(gl_RAY_B(1, :, 1))
    N1 = size(gl_RAY_B(:, 1, 1))

    ! Цикл движения точек на лучах B  
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
	
	
	
	N3 = size(gl_RAY_K(1, 1, :))
    N2 = size(gl_RAY_K(1, :, 1))
    N1 = size(gl_RAY_K(:, 1, 1))

    ! Цикл движения точек на лучах K  
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
	
	end if
	
	N3 = size(gl_RAY_O(1, 1, :))
    N2 = size(gl_RAY_O(1, :, 1))
    N1 = size(gl_RAY_O(:, 1, 1))

    ! Цикл движения точек на лучах O  
	if (.True.) then   ! Можно отключить этот цикл
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
	
	!print*, gl_x2(gl_RAY_B(3, 1, 1), now2)
	!print*, gl_y2(gl_RAY_B(6, 4, 4), now2)
	!print*, gl_z2(gl_RAY_B(5, 6, 5), now2)
	!print*, "---------------------------------------------------"
    
		
    N3 = size(gl_RAY_A(1, 1, :))
    N2 = size(gl_RAY_A(1, :, 1))
    N1 = size(gl_RAY_A(:, 1, 1))

    ! Цикл движения точек на лучах А  ************************************************************
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
			R_TS = norm2(KORD + proect * ER)  ! Новое расстояние до TS
			
			!if(j == 1 .and. k==1) write(*,*) "=    ", gl_x2(yzel, now)
            
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
            R_HP = norm2(KORD + proect * ER)  ! Новое расстояние до HP
            
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
            R_BS = norm2(KORD + proect * ER)  ! Новое расстояние до BS
			
            
            ! Далее обычный цикл нахождения координат точек, такой же, как и при построении сетки
            do i = 1, N1
                
                if (i == 1) then
                    CYCLE
                end if

                yzel = gl_RAY_A(i, j, k)
                ! Вычисляем координаты точки на луче

                ! до TS
				if (i <= par_n_IB) then  ! NEW
                    r =  par_R0 + (par_R_inner - par_R0) * (DBLE(i)/(par_n_IB))**par_kk1
                else if (i <= par_n_TS) then  ! До расстояния = R_TS
                    r =  par_R_inner + (R_TS - par_R_inner) * (DBLE(i - par_n_IB)/(par_n_TS - par_n_IB))**par_kk1
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
				
				!if(j == 1 .and. k==1 .and. i==20)  write(*,*) "=    ", gl_x2(yzel, now2)
                
            end do
        end do
	end do
    
		
    N3 = size(gl_RAY_B(1, 1, :))
    N2 = size(gl_RAY_B(1, :, 1))
    N1 = size(gl_RAY_B(:, 1, 1))

    ! Цикл движения точек на лучах B  ************************************************************
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
            R_TS = norm2(KORD + proect * ER)  ! Новое расстояние до TS
            
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
            
            proect = DOT_PRODUCT(vel * Time, ER)
			KORD(1) = gl_x2(yzel, now); KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now) 
            R_HP = norm2(KORD + proect * ER)  ! Новое расстояние до HP
			
                
            do i = 1, N1

                if (i == 1) CYCLE
                
                yzel = gl_RAY_B(i, j, k)
                ! Вычисляем координаты точки на луче

                ! до TS
				if (i <= par_n_IB) then  ! NEW
                    r =  par_R0 + (par_R_inner - par_R0) * (DBLE(i)/(par_n_IB))**par_kk1
                else if (i <= par_n_TS) then  ! До расстояния = R_TS
                    r =  par_R_inner + (R_TS - par_R_inner) * (DBLE(i - par_n_IB)/(par_n_TS - par_n_IB))**par_kk1
                !if (i <= par_n_TS) then  ! До расстояния = R_TS
                !    r =  par_R0 + (R_TS - par_R0) * (REAL(i, KIND = 4)/par_n_TS)**par_kk1
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
    
	
    N3 = size(gl_RAY_C(1, 1, :))
    N2 = size(gl_RAY_C(1, :, 1))
    N1 = size(gl_RAY_C(:, 1, 1))

    ! Цикл движения точек на лучах C  ************************************************************
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

    N3 = size(gl_RAY_O(1, 1, :))
    N2 = size(gl_RAY_O(1, :, 1))
    N1 = size(gl_RAY_O(:, 1, 1))


    ! Цикл движения точек на лучах O   ************************************************************
    do k = 1, N3
        phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
        
        do j = 1, N2
            
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
            R_HP = norm2(KORD + proect * ER )  ! Новое расстояние до HP
			
			! Блокируем схлопывание контакта к оси
			if(R_HP < 15.0_8) then
				R_HP = 15.0_8
			end if
			
            
            xx = gl_x2(gl_RAY_B(par_n_HP, par_m_BC, k), now2)              ! Отталкиваемся от x - координаты крайней точки B на гелиопаузе в этой плоскости (k)
            x = xx - (DBLE(j)/N2)**par_kk3 * (xx - par_R_LEFT)
            
            ! BS     Нужно взять положение BS из её положения на крайнем луче A
            yzel = gl_RAY_A(par_n_BS, size(gl_RAY_A(1, :, 1)), k)
			KORD(1) = 0.0_8; KORD(2) = gl_y2(yzel, now2); KORD(3) = gl_z2(yzel, now2)
            R_BS = norm2(KORD)  ! Новое расстояние до BS
			
			if (R_BS < R_HP + 10.0) then
				R_BS = R_HP + 10.0
				!print*, "R_BS = R_HP + 10.0"
			end if
			
            
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
    
    N3 = size(gl_RAY_K(1, 1, :))
    N2 = size(gl_RAY_K(1, :, 1))
    N1 = size(gl_RAY_K(:, 1, 1))

    ! Цикл движения точек на лучах K  ************************************************************
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
			R_TS = norm2(KORD + proect * ER)  ! Новое расстояние до TS
            
            do i = 1, N1

                if (i == 1) CYCLE
                
                ! Вычисляем координаты точки на луче
                yzel = gl_RAY_K(i, j, k)
				
				if (i <= par_n_IB) then  ! NEW
                    r =  par_R0 + (par_R_inner - par_R0) * (DBLE(i)/(par_n_IB))**par_kk1
                else 
                    r =  par_R_inner + (R_TS - par_R_inner) * (DBLE(i - par_n_IB)/(par_n_TS - par_n_IB))**par_kk1
				end if
					
                !r =  par_R0 + (R_TS - par_R0) * (REAL(i, KIND = 4)/par_n_TS)**par_kk1


                ! Записываем новые координаты
                gl_x2(yzel, now2) = r * cos(the)
                gl_y2(yzel, now2) = r * sin(the) * cos(phi)
                gl_z2(yzel, now2) = r * sin(the) * sin(phi)
                
                
            end do
        end do
    end do
    
    
     N3 = size(gl_RAY_D(1, 1, :))
    N2 = size(gl_RAY_D(1, :, 1))
    N1 = size(gl_RAY_D(:, 1, 1))

    ! Цикл движения точек на лучах D ************************************************************
    do k = 1, N3
        do j = 1, N2
            
            ! Вычисляем координаты текущего луча в пространстве
                phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
				
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
				
				r = sqrt(y**2 + z**2)
				
				
            do i = 1, N1

                if (i == 1) CYCLE

                if (k /= 1 .and. j == 1) CYCLE
				
				
				y = gl_y2(gl_RAY_O(1, i - 1, k), now2)
				z = gl_z2(gl_RAY_O(1, i - 1, k), now2)
				
				rr = sqrt(y**2 + z**2) * 0.3

                yzel = gl_RAY_D(i, j, k)
                
                
                x = xx + (DBLE(i - 1)/(N1 - 1))**par_kk3 * (par_R_LEFT - xx)
				rrr = r * (1.0 - (DBLE(i - 1)**2/(N1 - 1)**2)) + (rr * (DBLE(j - 1)/(N2 - 1))) * ( (DBLE(i - 1)**2) /(N1 - 1)**2)

				
				
                ! Записываем новые координаты
                gl_x2(yzel, now2) = x
                gl_y2(yzel, now2) = rrr * cos(phi)
                gl_z2(yzel, now2) = rrr * sin(phi)
                
            end do
        end do
    end do
    
    N3 = size(gl_RAY_E(1, 1, :))
    N2 = size(gl_RAY_E(1, :, 1))
    N1 = size(gl_RAY_E(:, 1, 1))
    
    ! Цикл движения точек на лучах E  ************************************************************
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
				
				!if(i == 14 .and. j == 12 .and. k == 35) then
				!	print*, "___"
				!	print*, x, x2
				!	print*, y, y2
				!	print*, z, z2
				!	print*, gl_x2(yzel, now2), gl_x2(yzel, now)
				!	print*, gl_y2(yzel, now2), gl_y2(yzel, now)
				!	print*, gl_z2(yzel, now2), gl_z2(yzel, now)
				!	print*, gl_y2(gl_RAY_E(1, j, k), now2), gl_y2(gl_RAY_E(1, j, k), now)
				!	print*, gl_z2(gl_RAY_E(1, j, k), now2), gl_z2(gl_RAY_E(1, j, k), now)
				!	print*, "___"
				!	print*, gl_x2(gl_RAY_O(1, j, k), now2), gl_x2(gl_RAY_O(1, j, k), now)
				!	print*, gl_y2(gl_RAY_O(1, j, k), now2), gl_y2(gl_RAY_O(1, j, k), now)
				!	print*, gl_z2(gl_RAY_O(1, j, k), now2), gl_z2(gl_RAY_O(1, j, k), now)
				!	print*, "___"
				!end if
				
                
            end do
        end do
	end do
 !
	!
	!print*, gl_x2(gl_RAY_A(10, 1, 1), now2)
	!print*, gl_y2(gl_RAY_A(14, 4, 4), now2)
	!print*, gl_z2(gl_RAY_A(18, 8, 5), now2)
	!print*, "---------------------------------------------------"
	!
	!print*, gl_x2(gl_RAY_B(3, 1, 1), now2)
	!print*, gl_y2(gl_RAY_B(6, 4, 4), now2)
	!print*, gl_z2(gl_RAY_B(5, 6, 5), now2)
	!print*, "---------------------------------------------------"
	!
	!print*, gl_x2(gl_RAY_C(3, 1, 1), now2)
	!print*, gl_y2(gl_RAY_C(6, 4, 4), now2)
	!print*, gl_z2(gl_RAY_C(5, 6, 5), now2)
	!print*, "---------------------------------------------------"
	!
	!print*, gl_x2(gl_RAY_O(3, 1, 1), now2)
	!print*, gl_y2(gl_RAY_O(6, 4, 4), now2)
	!print*, gl_z2(gl_RAY_O(5, 6, 5), now2)
	!print*, "---------------------------------------------------"
	!
	!print*, gl_x2(gl_RAY_K(3, 1, 1), now2)
	!print*, gl_y2(gl_RAY_K(6, 4, 4), now2)
	!print*, gl_z2(gl_RAY_K(5, 6, 5), now2)
	!print*, "---------------------------------------------------"
	!
	!print*, gl_x2(gl_RAY_D(3, 1, 1), now2)
	!print*, gl_y2(gl_RAY_D(6, 4, 4), now2)
	!print*, gl_z2(gl_RAY_D(5, 6, 5), now2)
	!print*, "---------------------------------------------------"
	!
	!print*, gl_x2(gl_RAY_E(3, 1, 1), now2)
	!print*, gl_y2(gl_RAY_E(6, 4, 4), now2)
	!print*, gl_z2(gl_RAY_E(5, 6, 5), now2)
	!print*, "---------------------------------------------------"
	!pause
    
    
    end subroutine Move_all
    
    subroutine Get_center_move(n, CENTER, now)  ! Получить массив - центр ячейки
    use STORAGE, only: gl_all_Cell, gl_x2, gl_y2, gl_z2
    implicit none

    integer, intent(in) :: n, now
    real(8), intent(out) :: CENTER(3)
    real(8) :: c(3)
    integer :: i, j

    CENTER = 0.0

    do i = 1, 8
        j = gl_all_Cell(i, n)
        c = (/gl_x2(j, now), gl_y2(j, now), gl_z2(j, now)/)
        CENTER = CENTER + c
    end do

    CENTER = CENTER/8.0

    end subroutine Get_center_move
    
    subroutine calc_all_Gran_move(now)   ! Программа расчёта площадей и нормалей граней и Объёмов ячеек
    ! Версия функции для расчёта на подвижной сетке
    use STORAGE
    use GEO_PARAM
    implicit none
    
    integer, intent(in) :: now           ! Откуда берём координаты узлов и куда записываем вычисленные параметры (площади и объёмы)
    integer :: Ngran, iter
    real(8) :: p(3, 4), Vol, D
    real(8) :: a(3), b(3), c(3), S, node1(3), node2(3)
    real(8) :: dist, di, gr_center(3)
    integer :: i, j, k, ll, grc, ii, now2

    !print*, "calc_all_Gran_move"
    now2 = mod(now, 2) + 1 
    ! Цикл по граням - считаем площадь грани, её нормаль
    Ngran = size(gl_all_Gran(1,:))

    do  iter = 1, Ngran
        grc = 0
        p = 0.0
        a = 0.0
        b = 0.0
        c = 0.0
        gr_center = 0.0
        node1 = 0.0
        node2 = 0.0
        ! Считываем из глобальной памяти координаты точек грани

        !print*, gl_all_Gran(:, iter)
        !pause

        if (gl_all_Gran(1, iter) == gl_all_Gran(3, iter)) then
            print*, "EROROR nierhfue 1", iter
            pause
        end if

        if (gl_all_Gran(2, iter) == gl_all_Gran(4, iter)) then
            print*, "EROROR nierhfue 2", iter
            pause
        end if


        do j = 1, 4
            i = gl_all_Gran(j, iter)
            p(:,j) = (/gl_x2(i, now), gl_y2(i, now), gl_z2(i, now)/)
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

        ! Можно один раз проверить, правильно ли ориентирована нормаль!\

        
        call Get_center_move(gl_Gran_neighbour(1, iter), node1, now)
        node2 = (p(:,1) + p(:,2) + p(:,3) + p(:,4))/4.0
        node1 = node2 - node1
        
        if(DOT_PRODUCT(node1, c) < 0.0) then
            print*, c
            print*, p(1, 1), norm2( (/0.0_8, p(2, 1), p(3, 1) /) )
            print*, "ERROR 13123 move ", DOT_PRODUCT(node1, c)
            print*, iter
            pause
        end if

        ! Нужно записать площадь грани и нормаль в общий массив!
        gl_Gran_normal2(:, iter, now) = c
        gl_Gran_square2(iter, now) = S

        !if (S < 0.000001) then
        !    print *, "ERROR 134443   gran = 0 ", c, a, b
        !    print *, p(:,1), p(:,2), p(:,3), p(:,4)
        !    pause
        !end if


        !print *, gl_Gran_square(iter)
        !pause
    end do

    


    ! Теперь посчитаем объёмы ячеек
    Ngran = size(gl_all_Cell(1,:))

    do  iter = 1, Ngran   ! Пробегаемся по всем ячейкам
        Vol = 0.0
        ll = 0
        dist = 10000.0 * par_R_character
        c = 0.0
        ! Для вычисления правильного центра ячейки, нужно понять какая она (некоторые узлы пропущены)
        if (gl_all_Cell(1, iter) /= gl_all_Cell(5, iter)) then
            do j = 1, 8
                i = gl_all_Cell(j, iter)
                c = c + (/gl_x2(i, now), gl_y2(i, now), gl_z2(i, now)/)
            end do
            c = c/8.0
        else
            if (gl_all_Cell(1, iter) == gl_all_Cell(2, iter)) then
                if (gl_all_Cell(4, iter) == gl_all_Cell(8, iter)) then  ! Класс ячейки = 4
                    do j = 1,8
                        ! if (j == 2 .or. j == 5 .or. j == 6 .or. j == 8) CYCLE
                        i = gl_all_Cell(j, iter)
                        c = c + (/gl_x2(i, now), gl_y2(i, now), gl_z2(i, now)/)
                    end do
                    c = c/8.0  ! 4
                else
                    do j = 1,8
                        ! if (j == 2 .or. j == 5 .or. j == 6) CYCLE
                        i = gl_all_Cell(j, iter)
                        c = c + (/gl_x2(i, now), gl_y2(i, now), gl_z2(i, now)/)
                    end do
                    c = c/8.0
                end if
            else
                if (gl_all_Cell(2, iter) == gl_all_Cell(6, iter)) then
                    if (gl_all_Cell(4, iter) == gl_all_Cell(5, iter)) then ! Класс ячейки = 1
                        do j = 1,8   ! Особая ячейка типа 4
                            ! if (j == 4 .or. j == 5 .or. j == 6 .or. j == 8) CYCLE
                            i = gl_all_Cell(j, iter)
                            c = c + (/gl_x2(i, now), gl_y2(i, now), gl_z2(i, now)/)
                        end do
                        c = c/8.0
                    else
                        do j = 1,8         ! Особая ячейка типа 2
                            ! if (j == 5 .or. j == 6) CYCLE             ! Попробовали убрать такой расчёт центра!
                            i = gl_all_Cell(j, iter)
                            c = c + (/gl_x2(i, now), gl_y2(i, now), gl_z2(i, now)/)
                        end do
                        c = c/8.0    ! 6 ----------------------------------------------------------------------------
                    end if
                else
                    if (gl_all_Cell(4, iter) == gl_all_Cell(5, iter)) then
                        do j = 1,8
                            ! if (j == 4 .or. j == 5 .or. j == 8) CYCLE
                            i = gl_all_Cell(j, iter)
                            c = c + (/gl_x2(i, now), gl_y2(i, now), gl_z2(i, now)/)
                        end do
                        c = c/8.0
                    else            ! Особая ячейка типа 5
                        do j = 1,8
                            ! if (j == 5 .or. j == 8) CYCLE    ! Попробовали убрать такой расчёт центра!
                            i = gl_all_Cell(j, iter)
                            c = c + (/gl_x2(i, now), gl_y2(i, now), gl_z2(i, now)/)
                        end do
                        c = c/8.0   ! 6 ----------------------------------------------------------------------------
                    end if
                end if
            end if
        end if

        gl_Cell_center2(:, iter, now) = c

        !print *, c
        !pause
		
		

        ! Вычислили центр ячейки теперь считаем объём пирамиды на каждую грань
        do j = 1, 6
            i = gl_Cell_gran(j, iter)   ! Берём по очереди все грани ячейки

            !if (iter == gl_Cell_C(1, size(gl_Cell_C(1,:,1)) ,1)) then
            !    print*, gl_Cell_gran(:, iter)
            !    pause
            !end if

            if (i == 0) CYCLE
            k = gl_all_Gran(1, i)  ! Номер первого узла грани
            b = (/gl_x2(k, now), gl_y2(k, now), gl_z2(k, now)/)
            a = gl_Gran_normal2(:,i, now)
            di = dabs(DOT_PRODUCT(c,a) - DOT_PRODUCT(a,b))  ! Расстояние от точки до плоскости
            Vol = Vol + gl_Gran_square2(i, now) * di/3.0
            dist = MIN(dist, di)
            ll = ll + 1
        end do



        !if (iter == gl_Cell_A(5, 1, 1)) then
        !    write(*, *) c
        !    print*, Vol
        !    print*, gl_Cell_gran(:, iter)
        !    pause
        !end if


        gl_Cell_Volume2(iter, now) = Vol
        gl_Cell_dist(iter) = dist
		
		
		

        !        print*, Vol, ll
        !        pause

    end do


    end subroutine calc_all_Gran_move
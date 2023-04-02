!******7**10**********************************************************70
!**                                                                  **
!**    chlld   for                                                   **
!**                                                                  **
!******7**10**********************************************************70
!**    shema HLLD(modification) for 3D MHD                           **
!******7**10**********************************************************70

	
	subroutine Cuda_Calc_move(now)
	
	use MY_CUDA, gl_TS => dev_gl_TS, gl_Gran_neighbour => dev_gl_Gran_neighbour, gl_Gran_normal2 => dev_gl_Gran_normal2, &
		gl_Cell_par => dev_gl_Cell_par, gl_all_Gran => dev_gl_all_Gran, gl_Point_num => dev_gl_Point_num, gl_x2 => dev_gl_x2, &
		gl_y2 => dev_gl_y2, gl_z2 => dev_gl_z2, gl_Gran_center2 => dev_gl_Gran_center2, norm2 => dev_norm2, gl_Vx => dev_gl_Vx, &
		gl_Vy => dev_gl_Vy, gl_Vz => dev_gl_Vz, gl_Contact => dev_gl_Contact, gl_BS => dev_gl_BS
	use GEO_PARAM
	implicit none
	integer, intent(in) :: now
	
	integer, device :: i, Num, gr, s1, s2, j, yzel, metod
    real(8), device :: normal(3), qqq1(8), qqq2(8), dsl, dsp, dsc, POTOK(8), a1, a2, v1, v2, ray(3), norm, b1, b2, c1, c2, koef1, koef2, koef3, www
    
    koef1 = 0.3
    koef2 = 0.7
    koef3 = 0.3
	
	metod = 2
	www = 0.0
    
    ! Пробегаемся по всем граням и вычисляем скорости их движения
    
    ! TS
    Num = size(gl_TS)
    
	!$cuf kernel do <<<*, *>>>
    do i = 1, Num
        
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
        
        call chlld(metod, normal(1), normal(2), normal(3), &
                www, qqq1, qqq2, dsl, dsp, dsc, POTOK)
        
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
        
        call chlld(metod, normal(1), normal(2), normal(3), &
                www, qqq1, qqq2, dsl, dsp, dsc, POTOK)
        
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
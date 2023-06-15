
    
    subroutine Find_Surface()   ! ��������� ������� ������������, ������� ����������
    ! ��� BS ��������� ������ ����������� �� ���������! (x > 0)
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
			
			! ������ ����� ������ ������� 
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
			
			! ������ ����� ������ ������� 
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
			
			! ������ ����� ������ ������� 
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
			
			! ������ ����� ������ ������� 
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
	
	! ������ ����� �� ��� ���������
	
	do k = 1, size( gl_Cell_A(1, 1, :) )
        do j = 1, size( gl_Cell_A(:, 1, 1) )
			! ������ ����� ������
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
			! ������ ����� ������
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
	
	
	! ��������� ��� ������ ������ � ����
	 do k = 1, size( gl_Cell_A(1, 1, :) )
        do j = 1, size( gl_Cell_A(1, :, 1) )
			do i = 1, par_n_TS - 1
				num = gl_Cell_A(i, j, k)
				gl_zone_Cell(num) = 1
			end do
		end do
	 end do
	 
	 do k = 1, size( gl_Cell_A(1, 1, :) )
        do j = 1, size( gl_Cell_A(1, :, 1) )
			do i = par_n_HP, par_n_BS - 1
				num = gl_Cell_A(i, j, k)
				gl_zone_Cell(num) = 3
			end do
		end do
	 end do
	 
	 do k = 1, size( gl_Cell_A(1, 1, :) )
        do j = 1, size( gl_Cell_A(1, :, 1) )
			do i = par_n_BS, size( gl_Cell_A(:, 1, 1) )
				num = gl_Cell_A(i, j, k)
				gl_zone_Cell(num) = 4
			end do
		end do
	 end do
	 
	 do k = 1, size( gl_Cell_B(1, 1, :) )
        do j = 1, size( gl_Cell_B(1, :, 1) )
			do i = 1, par_n_TS - 1
				num = gl_Cell_B(i, j, k)
				gl_zone_Cell(num) = 1
			end do
		end do
	 end do
	 
	
	do k = 1, size( gl_Cell_C(1, 1, :) )
        do j = 1, size( gl_Cell_C(1, :, 1) )
			do i = par_n_HP - par_n_TS + 1, par_n_BS - par_n_TS
				num = gl_Cell_C(i, j, k)
				gl_zone_Cell(num) = 3
			end do
		end do
	end do
	
	do k = 1, size( gl_Cell_C(1, 1, :) )
        do j = 1, size( gl_Cell_C(1, :, 1) )
			do i = par_n_BS - par_n_TS + 1, size( gl_Cell_C(:, 1, 1) )
				num = gl_Cell_C(i, j, k)
				gl_zone_Cell(num) = 4
			end do
		end do
	 end do
			
	

    end subroutine Find_Surface
    
    subroutine Calc_move(now) ! ��������� �������� �������� ����� (��������) � ������� ������� �� ������������ �������(
    use STORAGE
    use GEO_PARAM
	use Solvers
	use My_func
    implicit none
    
    integer, intent(in) :: now   ! ������ ����� ���������� ������������ � ���������� ��������
    
    integer :: i, Num, gr, s1, s2, j, yzel, ndisc
    real(8) :: normal(3), qqq1(8), qqq2(8), dsl, dsp, dsc, POTOK(8), a1, a2, v1, v2, &
        ray(3), norm, b1, b2, c1, c2, vec(3), center(3), center2(3), the1, the2
    
    !koef1 = 0.3! 0.2
    !koef2 = 1.0 ! 1.0
    !koef3 = 0.7   ! 0.3
	
	ndisc = 0
    
    ! ����������� �� ���� ������ � ��������� �������� �� ��������
    
    ! TS
    Num = size(gl_TS)
    
    do i = 1, Num
        
        gr = gl_TS(i)
    	s1 = gl_Gran_neighbour(1, gr)
        s2 = gl_Gran_neighbour(2, gr)
        normal = gl_Gran_normal2(:, gr, now)
        qqq1 = gl_Cell_par(1:8, s1)
        qqq2 = gl_Cell_par(1:8, s2)
        
        call chlld(3, normal(1), normal(2), normal(3), &  ! 2
                0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)
		
		
        
        dsl = dsl * koef1
		
		do j = 1, 4
                yzel = gl_all_Gran(j, gr)
                gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsl
                gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsl
                gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsl
                gl_Point_num(yzel) = gl_Point_num(yzel) + 1
		end do
		CYCLE  ! ����������� � ���� ������, ��������� � ���������
			
		
		if(.False.) then
		    center = gl_Gran_center2(:, gr, now)
		
		    the1 = polar_angle(center(1), sqrt(center(2)**2 + center(3)**2))
		
		    do j = 1, 4
                yzel = gl_all_Gran(j, gr)
			
			    center2(1) = gl_x2(yzel, now)
			    center2(2) = gl_y2(yzel, now)
			    center2(3) = gl_z2(yzel, now)
			    the2 = polar_angle(center2(1), sqrt(center2(2)**2 + center2(3)**2))
			
			    if (the2 < the1 .and. the2 > 0.18) CYCLE
			
			
                gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsl
                gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsl
                gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsl
                gl_Point_num(yzel) = gl_Point_num(yzel) + 1
            end do
            CYCLE  ! ����������� � ���� ������, ��������� � ���������
		end if
		
        !if(i == 50) write(*,*) dsl
		
        a1 = sqrt(ggg * qqq1(5)/qqq1(1))  ! �������� �����
        a2 = sqrt(ggg * qqq2(5)/qqq2(1))
		
        
        if (norm2(qqq1(2:4)) <= a1 .and. norm2(qqq2(2:4)) <= a2) then ! ���������� �������� �� ��� ����
            do j = 1, 4
            	yzel = gl_all_Gran(j, gr)
                gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsl
                gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsl
                gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsl
                gl_Point_num(yzel) = gl_Point_num(yzel) + 1
				
            end do
            CYCLE  ! ����������� � ���� ������, ��������� � ���������
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
					
                    CYCLE   ! ��������� � ���������� ����
                end if
                
                b1 = DOT_PRODUCT(qqq1(6:8), ray)
                b2 = DOT_PRODUCT(qqq2(6:8), ray)
                
                c1 = dabs(b1)/dsqrt(4.0 * par_pi_8 * qqq1(1))
                c2 = dabs(b2)/dsqrt(4.0 * par_pi_8 * qqq2(1))
                
                b1 = norm2(qqq1(6:8))
                b2 = norm2(qqq2(6:8))
                
                ! ���� ������ �������� ������ ��� �������� �� �������� ����� � ������������� �������� ����� �� ����� �����������
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
        
    ! �������
    Num = size(gl_Contact)
    
    do i = 1, Num
        
        gr = gl_Contact(i)
    	s1 = gl_Gran_neighbour(1, gr)
        s2 = gl_Gran_neighbour(2, gr)
        normal = gl_Gran_normal2(:, gr, now)
        qqq1 = gl_Cell_par(1:8, s1)
        qqq2 = gl_Cell_par(1:8, s2)
		
		! ������ ���������� ���������� ���������� ����
		if (gl_Gran_center2(1, gr, now) > -100.0) then
	        qqq1(6:8) = qqq1(6:8) - DOT_PRODUCT(normal, qqq1(6:8)) * normal
	        qqq2(6:8) = qqq2(6:8) - DOT_PRODUCT(normal, qqq2(6:8)) * normal
	        gl_Cell_par(6:8, s1) = qqq1(6:8)
	        gl_Cell_par(6:8, s2) = qqq2(6:8)
		end if
		
		
        
        call chlld(3, normal(1), normal(2), normal(3), &  ! 2
                0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)
		
		vec = 0.5 * (qqq1(2:4) + qqq2(2:4))
		vec = vec - DOT_PRODUCT(vec, normal) * normal
		center = gl_Gran_center2(:, gr, now)
		
        
        !dsc = dsc * koef2
		dsc = (dsc + 0.0 * DOT_PRODUCT(0.5 * (qqq1(2:4) + qqq2(2:4)), normal)) * koef2
		
		!if (gl_Gran_center2(1, gr, now) < -200.0 .and. gl_Gran_center2(1, gr, now) > -220.0 .and. normal(2) > 0 .and. &
		!	dabs(gl_Gran_center2(3, gr, now)) < 5.0) then
		!	print*, gl_Gran_center2(:, gr, now) 
		!	dsc = dsc + 30.0
		!end if
		
		!if(gr == 77) write(*, *) dsc, qqq1(1), qqq2(1), qqq1(5), qqq2(5)
		
		do j = 1, 4
            yzel = gl_all_Gran(j, gr)
			
			
			center2(1) = gl_x2(yzel, now)
			center2(2) = gl_y2(yzel, now)
			center2(3) = gl_z2(yzel, now)
			
			if ( sqrt(center2(2)**2 + center2(3)**2) > 15.0) then
			    if( DOT_PRODUCT(vec, center2 - center) < 0.0) CYCLE
			end if
			
			
            gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsc
            gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsc
            gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsc
            gl_Point_num(yzel) = gl_Point_num(yzel) + 1
        end do
        CYCLE  ! ����������� � ���� ������, ��������� � ���������
		
		
        !if(i == 2) write(*,*) dsc
        a1 = sqrt(ggg * qqq1(5)/qqq1(1))  ! �������� �����
        a2 = sqrt(ggg * qqq2(5)/qqq2(1))
		
        
        if (norm2(qqq1(2:4)) <= a1 .and. norm2(qqq2(2:4)) <= a2) then ! ���������� �������� �� ��� ����
            do j = 1, 4
            	yzel = gl_all_Gran(j, gr)
                gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsc
                gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsc
                gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsc
                gl_Point_num(yzel) = gl_Point_num(yzel) + 1
            end do
            CYCLE  ! ����������� � ���� ������, ��������� � ���������
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
                    CYCLE   ! ��������� � ���������� ����
                end if
                
                b1 = DOT_PRODUCT(qqq1(6:8), ray)
                b2 = DOT_PRODUCT(qqq2(6:8), ray)
                
                c1 = dabs(b1)/dsqrt(4.0 * par_pi_8 * qqq1(1))
                c2 = dabs(b2)/dsqrt(4.0 * par_pi_8 * qqq2(1))
                
                b1 = norm2(qqq1(6:8))
                b2 = norm2(qqq2(6:8))
                
                ! ���� ������ �������� ������ ��� �������� �� �������� ����� � ������������� �������� ����� �� ����� �����������
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
        
            a1 = sqrt(ggg * qqq1(5)/qqq1(1))  ! �������� �����
            a2 = sqrt(ggg * qqq2(5)/qqq2(1))
        
            if (norm2(qqq1(2:4)) <= a1 .and. norm2(qqq2(2:4)) <= a2) then ! ���������� �������� �� ��� ����
                do j = 1, 4
            	    yzel = gl_all_Gran(j, gr)
                    gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsp
                    gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsp
                    gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsp
                    gl_Point_num(yzel) = gl_Point_num(yzel) + 1
                end do
                CYCLE  ! ����������� � ���� ������, ��������� � ���������
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
                        CYCLE   ! ��������� � ���������� ����
                    end if
                
                    b1 = DOT_PRODUCT(qqq1(6:8), ray)
                    b2 = DOT_PRODUCT(qqq2(6:8), ray)
                
                    c1 = dabs(b1)/dsqrt(4.0 * par_pi_8 * qqq1(1))
                    c2 = dabs(b2)/dsqrt(4.0 * par_pi_8 * qqq2(1))
                
                    b1 = norm2(qqq1(6:8))
                    b2 = norm2(qqq2(6:8))
                
                    ! ���� ������ �������� ������ ��� �������� �� �������� ����� � ������������� �������� ����� �� ����� �����������
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
    
    subroutine Move_all(now, Time)  ! ����������� ����� � ������������ �� ��������� �������� ����� �� ������������ �������
    use STORAGE
    use GEO_PARAM
	use My_func
    implicit none
    
    integer :: now2             ! ��� ��������� �� ������ ������ �� ������ now
    real(8), intent(in) :: Time
    integer, intent(in) :: now
    real(8) :: R_TS, proect, vel(3), R_HP, R_BS, Ak(3), Bk(3), Ck(3), Dk(3), Ek(3), ER(3), KORD(3), dist, ddt
    integer :: yzel, N1, N2, N3, i, j, k, yzel2
    real(8) :: the, phi, r, x, y, z, rr, xx, x2, y2, z2, rrr, r1, r2, r3, r4, rd, kk13
    
    now2 = mod(now, 2) + 1   
	
	ddt = Time/0.000127
	
	! ������������� ��������� ------------------------------------------------------------------------------
	if (.True.) then
	N3 = size(gl_RAY_A(1, 1, :))
    N2 = size(gl_RAY_A(1, :, 1))
    N1 = size(gl_RAY_A(:, 1, 1))

    ! ���� �������� ����� �� ����� �  
    do k = 1, N3
        do j = 1, N2
			
			if (j == 1) then ! .and. k /= 1) then  ! k /= 1 .and. 
                    CYCLE
			end if
			
			! ������� �����
			
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
				yzel2 = gl_RAY_A(par_n_TS, j + 1, k + ceiling(1.0 * N3/2))
			end if
			Dk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			if(k < N3) then
			    yzel2 = gl_RAY_A(par_n_TS, j, k + 1)
			else
				yzel2 = gl_RAY_A(par_n_TS, j, 1)
			end if
			Ek = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			
			
			
			if (gl_Point_num(yzel) > 0) then
			    vel = gl_Point_num(yzel) * par_nat_TS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time * ddt
			else
				vel = par_nat_TS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time * ddt
			end if
			
			!if (k == 22) then
			!	print*, "___"
			!	dist = sqrt(gl_Vx(yzel)**2 + gl_Vy(yzel)**2 + gl_Vz(yzel)**2)
			!	print*, norm2(Ak), Ak(1), dist
			!	print*, norm2(vel), dist, 100.0 - (dist - norm2(vel))/dist * 100
			!    end if
			
			
			gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
			gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
			gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
			
			
			! �������
			
			yzel = gl_RAY_A(par_n_HP, j, k)
			Ak = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/)
			r = norm2(Ak)
			
			if (j < N2) then
			    yzel2 = gl_RAY_A(par_n_HP, j + 1, k)
			else
				yzel2 = gl_RAY_B(par_n_HP, 1, k)
			end if
			Bk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			r1 = norm2(Bk)
			
			if (k > 1) then
			    yzel2 = gl_RAY_A(par_n_HP, j, k - 1)
			else
			    yzel2 = gl_RAY_A(par_n_HP, j, N3)
			end if
			Ck = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			r2 = norm2(Ck)
			
			if (j > 1) then
			    yzel2 = gl_RAY_A(par_n_HP, j - 1, k)
			else
				yzel2 = gl_RAY_A(par_n_HP, j, k + ceiling(1.0 * N3/2))
			end if
			Dk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			r3 = norm2(Dk)
			
			if(k < N3) then
			    yzel2 = gl_RAY_A(par_n_HP, j, k + 1)
			else
				yzel2 = gl_RAY_A(par_n_HP, j, 1)
			end if
			Ek = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			r4 = norm2(Ek)
			
			dist = sqrt( (Dk(1) - Ak(1))**2 + (Dk(2) - Ak(2))**2 + (Dk(3) - Ak(3))**2 )
	        dist = max(dist, 1.0_8)
			
			rr = (r1 + r2 + r3 + r4)/4.0
			
			if (gl_Point_num(yzel) > 0) then
			    !vel = gl_Point_num(yzel) * par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time/dist * ddt
				vel = par_nat_HP * 0.0001 * gl_Point_num(yzel) * ((Ak/r) * (rr - r)) * ddt  ! 0.035
			else
				!vel = par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time/dist * ddt
				vel = par_nat_HP * 0.0001 * ((Ak/r) * (rr - r)) * ddt
			end if
			
			!if(k == 2 .and. j == 2) then
			!	print*, gl_Vx(yzel), gl_Vy(yzel), gl_Vz(yzel)
			!	print*, vel
			!	print*, "_____"
			!	pause
			!end if
			
			
			
			gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
			gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
			gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
			
			
			! ������� ������� �����
			
			yzel = gl_RAY_A(par_n_BS, j, k)
			Ak = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/)
			
			if (j < N2) then
			    yzel2 = gl_RAY_A(par_n_BS, j + 1, k)
			else
				CYCLE !yzel2 = yzel                                                    ! ��� ��� ���� ��������������� ��� ������ � ������� �����
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
			
			dist = (sqrt( (Dk(1) - Ak(1))**2 + (Dk(2) - Ak(2))**2 + (Dk(3) - Ak(3))**2 ) + &
				sqrt( (Ek(1) - Ak(1))**2 + (Ek(2) - Ak(2))**2 + (Ek(3) - Ak(3))**2 ))/2.0
	
	        dist = max(dist, 1.0_8)
			
			if (gl_Point_num(yzel) > 0) then
			    vel = gl_Point_num(yzel) * par_nat_BS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time/dist * ddt
			else
				vel = par_nat_BS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time/dist * ddt
			end if
			
			
			
			gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
			gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
			gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
			
		end do
	end do
	
	
	N3 = size(gl_RAY_B(1, 1, :))
    N2 = size(gl_RAY_B(1, :, 1))
    N1 = size(gl_RAY_B(:, 1, 1))

    ! ���� �������� ����� �� ����� B  
    do k = 1, N3
        do j = 1, N2
			
			! ������� �����
			
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
			    vel = gl_Point_num(yzel) * par_nat_TS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time * ddt
			else
				vel = par_nat_TS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time * ddt
			end if
			
			
			
			gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
			gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
			gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
			
			
			! �������
			
			yzel = gl_RAY_B(par_n_HP, j, k)
			Ak(1) = gl_x2(yzel, now); Ak(2) = gl_y2(yzel, now); Ak(3) = gl_z2(yzel, now)
			! Ak = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/)
			r = sqrt(Ak(2)**2 + Ak(3)**2)
			
			if (j < N2) then
			    yzel2 = gl_RAY_B(par_n_HP, j + 1, k)
			else
				!yzel2 = gl_RAY_O(1, 1, k)
				CYCLE
			end if
			Bk(1) = gl_x2(yzel2, now); Bk(2) = gl_y2(yzel2, now); Bk(3) = gl_z2(yzel2, now)
			! Bk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			r1 = sqrt(Bk(2)**2 + Bk(3)**2)
			
			!if (k > 1) then
			!    yzel2 = gl_RAY_B(par_n_HP, j, k - 1)
			!else
			!    yzel2 = gl_RAY_B(par_n_HP, j, N3)
			!end if
			!Ck(1) = gl_x2(yzel2, now); Ck(2) = gl_y2(yzel2, now); Ck(3) = gl_z2(yzel2, now)
			!! Ck = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			!r2 = sqrt(Ck(2)**2 + Ck(3)**2)
			
			if (j > 1) then
			    yzel2 = gl_RAY_B(par_n_HP, j - 1, k)
			else
				yzel2 = gl_RAY_A(par_n_HP, size( gl_RAY_A(par_n_HP, :, k)), k)
			end if
			Dk(1) = gl_x2(yzel2, now); Dk(2) = gl_y2(yzel2, now); Dk(3) = gl_z2(yzel2, now)
			! Dk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			r3 = sqrt(Dk(2)**2 + Dk(3)**2)
			
			!if(k < N3) then
			!    yzel2 = gl_RAY_B(par_n_HP, j, k + 1)
			!else
			!	yzel2 = gl_RAY_B(par_n_HP, j, 1)
			!end if
			!Ek(1) = gl_x2(yzel2, now); Ek(2) = gl_y2(yzel2, now); Ek(3) = gl_z2(yzel2, now)
			!! Ek = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
			!r4 = sqrt(Ek(2)**2 + Ek(3)**2)
			
			!rr = (r1 + r2 + r3 + r4)/4.0
			rr = (r1 + r3)/2.0
			
			!dist = sqrt( (Dk(1) - Ak(1))**2 + (Dk(2) - Ak(2))**2 + (Dk(3) - Ak(3))**2 )
			!dist = max(dist, 1.0_8)
			
			if (gl_Point_num(yzel) > 0) then
			    !vel = gl_Point_num(yzel) * par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time/dist
				vel = par_nat_HP * 0.001 * gl_Point_num(yzel) * ((Ak/r) * (rr - r)) * ddt
				vel(1) = 0.0
			else
				!vel = par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time/dist
				vel = par_nat_HP * 0.001 * ((Ak/r) * (rr - r)) * ddt
				vel(1) = 0.0
			end if
			
			
			
			gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
			gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
			gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
			
		end do
	end do
	
	
	
	N3 = size(gl_RAY_K(1, 1, :))
    N2 = size(gl_RAY_K(1, :, 1))
    N1 = size(gl_RAY_K(:, 1, 1))

    ! ���� �������� ����� �� ����� K  
    do k = 1, N3
        do j = 1, N2
			
			if (j == 1) then
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
			    vel = gl_Point_num(yzel) * par_nat_TS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time * ddt
			else
				vel = par_nat_TS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time * ddt
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

    ! ���� �������� ����� �� ����� O  
	if (.True.) then   ! ����� ��������� ���� ����
        do k = 1, N3
            do j = 1, N2
			
			    ! �������
			    yzel = gl_RAY_O(1, j, k)
				Ak(1) = gl_x2(yzel, now); Ak(2) = gl_y2(yzel, now); Ak(3) = gl_z2(yzel, now)
			    ! Ak = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/)
				r = sqrt(Ak(2)**2 + Ak(3)**2)
			
			    if (j < N2) then
			        yzel2 = gl_RAY_O(1, j + 1, k)
			    else
				    yzel2 = yzel
				end if
				Bk(1) = gl_x2(yzel2, now); Bk(2) = gl_y2(yzel2, now); Bk(3) = gl_z2(yzel2, now)
			    ! Bk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
				r1 = sqrt(Bk(2)**2 + Bk(3)**2)
			
			 !   if (k > 1) then
			 !       yzel2 = gl_RAY_O(1, j, k - 1)
			 !   else
			 !       yzel2 = gl_RAY_O(1, j, N3)
				!end if
				!Ck(1) = gl_x2(yzel2, now); Ck(2) = gl_y2(yzel2, now); Ck(3) = gl_z2(yzel2, now)
			 !   ! Ck = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
				!r2 = sqrt(Ck(2)**2 + Ck(3)**2)
			
			    if (j > 1) then
			        yzel2 = gl_RAY_O(1, j - 1, k)
			    else
				    CYCLE
				    !yzel2 = gl_RAY_C(1, size(gl_RAY_C(1, :, k)), k)
				end if
				Dk(1) = gl_x2(yzel2, now); Dk(2) = gl_y2(yzel2, now); Dk(3) = gl_z2(yzel2, now)
			    ! Dk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
				r3 = sqrt(Dk(2)**2 + Dk(3)**2)
			
			 !   if(k < N3) then
			 !       yzel2 = gl_RAY_O(1, j, k + 1)
			 !   else
				!    yzel2 = gl_RAY_O(1, j, 1)
				!end if
				!Ek(1) = gl_x2(yzel2, now); Ek(2) = gl_y2(yzel2, now); Ek(3) = gl_z2(yzel2, now)
			 !   ! Ek = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
				!r4 = sqrt(Ek(2)**2 + Ek(3)**2)
				
				rr = (r1 + r3)/2.0
			
			    if (gl_Point_num(yzel) > 0) then
			        !vel = gl_Point_num(yzel) * par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time
					vel = par_nat_HP * 0.0001 * gl_Point_num(yzel) * ((Ak/r) * (rr - r)) * ddt ! 0.0003
				    vel(1) = 0.0
				else
				    !vel = par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time  
					vel = par_nat_HP * 0.0001 * ((Ak/r) * (rr - r)) * ddt   ! 0.0003
				    vel(1) = 0.0
				end if
				
				!if (k == 22) then
				!print*, "___"
				!dist = sqrt(gl_Vx(yzel)**2 + gl_Vy(yzel)**2 + gl_Vz(yzel)**2)
				!print*, norm2(Ak), Ak(1), dist
				!print*, norm2(vel), dist, 100.0 - (dist - norm2(vel))/dist * 100
			 !   end if
				
			
			
			    gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
			    gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
			    gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
			
		    end do
	    end do
	end if
	
	
	
	! ����� �������������� ��������� --------------------------------------------------------------------------------------------------------------------
	
	!print*, gl_x2(gl_RAY_B(3, 1, 1), now2)
	!print*, gl_y2(gl_RAY_B(6, 4, 4), now2)
	!print*, gl_z2(gl_RAY_B(5, 6, 5), now2)
	!print*, "---------------------------------------------------"
    
		
    N3 = size(gl_RAY_A(1, 1, :))
    N2 = size(gl_RAY_A(1, :, 1))
    N1 = size(gl_RAY_A(:, 1, 1))

    ! ���� �������� ����� �� ����� �  ************************************************************
    do k = 1, N3
        do j = 1, N2
            
            if (k /= 1 .and. j == 1) then
                    CYCLE
            end if
            
            ! ��������� ���������� �������� ���� � ������������
            the = (j - 1) * par_pi_8/2.0/(N2 - 1)
            phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
            
            ! TS
            yzel = gl_RAY_A(par_n_TS, j, k)
            if(gl_Point_num(yzel) == 0) then
                vel = 0.0
            else
                ! vel = (/gl_Vx(yzel), gl_Vy(yzel), gl_Vz(yzel)/)
				vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
                vel = vel/gl_Point_num(yzel)                       ! ����� �������� �������� ����� ����
            end if
            
            ! ������� ��� ���������� �������������
            gl_Point_num(yzel) = 0
            gl_Vx(yzel) = 0.0
            gl_Vy(yzel) = 0.0
            gl_Vz(yzel) = 0.0
            
			ER(1) = cos(the); ER(2) = sin(the) * cos(phi); ER(3) = sin(the) * sin(phi)
			proect = DOT_PRODUCT(vel * Time, ER)
            ! proect = DOT_PRODUCT(vel * Time, (/cos(the), sin(the) * cos(phi), sin(the) * sin(phi)/))  !  ������� �������� ����������� �� ������ ������ ����
            KORD(1) = gl_x2(yzel, now); KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now) 
			R_TS = norm2(KORD + proect * ER)  ! ����� ���������� �� TS
			
			!if(j == 1 .and. k==1) write(*,*) "=    ", gl_x2(yzel, now)
            
            ! HP
            yzel = gl_RAY_A(par_n_HP, j, k)
            if(gl_Point_num(yzel) == 0) then
                vel = 0.0
            else
                ! vel = (/gl_Vx(yzel), gl_Vy(yzel), gl_Vz(yzel)/)
				vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
                vel = vel/gl_Point_num(yzel)                       ! ����� �������� �������� ����� ����
			end if
            
            ! ������� ��� ���������� �������������
            gl_Point_num(yzel) = 0
            gl_Vx(yzel) = 0.0
            gl_Vy(yzel) = 0.0
            gl_Vz(yzel) = 0.0
            
            proect = DOT_PRODUCT(vel * Time, ER)
			KORD(1) = gl_x2(yzel, now); KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now) 
            R_HP = norm2(KORD + proect * ER)  ! ����� ���������� �� HP
			
		!	if(j == 1 .and. k == 1) then
		!		write(*,*) yzel, now
		!		write(*,*) R_HP
		!		write(*,*) vel(1), vel(2), vel(3)
		!		write(*,*) Time, proect
		!		write(*,*) ER(1), ER(2), ER(3)
		!		write(*, *) KORD(1), KORD(2), KORD(3)
		!		write(*,*) "______________________"
		!end if
			
            
            ! BS
            yzel = gl_RAY_A(par_n_BS, j, k)
            if(gl_Point_num(yzel) == 0) then
                vel = 0.0
            else
                ! vel = (/gl_Vx(yzel), gl_Vy(yzel), gl_Vz(yzel)/)
				vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
                vel = vel/gl_Point_num(yzel)                       ! ����� �������� �������� ����� ����
			end if
			
            
            ! ������� ��� ���������� �������������
            gl_Point_num(yzel) = 0
            gl_Vx(yzel) = 0.0
            gl_Vy(yzel) = 0.0
            gl_Vz(yzel) = 0.0
            
            proect = DOT_PRODUCT(vel * Time, ER)
			KORD(1) = gl_x2(yzel, now); KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now) 
            R_BS = norm2(KORD + proect * ER)  ! ����� ���������� �� BS
			
            
            ! ����� ������� ���� ���������� ��������� �����, ����� ��, ��� � ��� ���������� �����
            do i = 1, N1
                
                if (i == 1) then
                    CYCLE
                end if

                yzel = gl_RAY_A(i, j, k)
                ! ��������� ���������� ����� �� ����
				
				kk13 = par_kk13 * (par_pi_8/2.0 - the)/(par_pi_8/2.0)  +  1.0 * (the/(par_pi_8/2.0))**2

                ! �� TS
				if (i <= par_n_IB) then  ! NEW
						r =  par_R0 + (par_R_inner - par_R0) * (DBLE(i)/(par_n_IB))**par_kk1
				else if (i <= par_n_TS) then  
						r =  par_R_inner + (R_TS - par_R_inner) * (DBLE(i - par_n_IB)/(par_n_TS - par_n_IB))**par_kk12
				!if (i <= par_n_TS) then  ! �� ���������� = R_TS
				!    r =  par_R0 + (R_TS - par_R0) * (DBLE(i)/par_n_TS)**par_kk1
				else if (i <= par_n_HP) then  
					r = R_TS + (i - par_n_TS) * (R_HP - R_TS)/(par_n_HP - par_n_TS)
				else if (i == par_n_HP + 1) then 
					rd = R_TS + (par_n_HP - 1 - par_n_TS) * (R_HP - R_TS)/(par_n_HP - par_n_TS)
					r = R_HP + (R_HP - rd)
				else if (i <= par_n_BS) then 
					rd = R_TS + (par_n_HP - 1 - par_n_TS) * (R_HP - R_TS)/(par_n_HP - par_n_TS)
					rd = R_HP + (R_HP - rd)
					r = rd + (R_BS - rd) * sgushenie_1( 1.0_8 * (i - par_n_HP - 1)/(par_n_BS - par_n_HP - 1), kk13 )
					!r = R_HP + (R_BS - R_HP) * angle_cilindr( 1.0_8 * (i - par_n_HP)/(par_n_BS - par_n_HP), kk13 )
					!r = R_HP + (i - par_n_HP) * (R_BS - R_HP)/(par_n_BS - par_n_HP)
				else
					r = R_BS + (par_R_END - R_BS) * (DBLE(i- par_n_BS)/(par_n_END - par_n_BS))**(par_kk2 * (0.55 + 0.45 * cos(the)) )
				end if
				
				if (i == par_n_TS - 1) then
					r = R_TS - 0.5               
				end if

                ! ���������� ����� ����������
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

    ! ���� �������� ����� �� ����� B  ************************************************************
    do k = 1, N3
        do j = 1, N2
            
            ! ��������� ���������� �������� ���� � ������������
            the = par_pi_8/2.0 + (j) * par_triple_point/(N2)
            phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
            
            ! TS
            yzel = gl_RAY_B(par_n_TS, j, k)
            if(gl_Point_num(yzel) == 0) then
                vel = 0.0
			else
				vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
                vel = vel/gl_Point_num(yzel)                       ! ����� �������� �������� ����� ����
            end if
            
            ! ������� ��� ���������� �������������
            gl_Point_num(yzel) = 0
            gl_Vx(yzel) = 0.0
            gl_Vy(yzel) = 0.0
            gl_Vz(yzel) = 0.0
            
			ER(1) = cos(the); ER(2) = sin(the) * cos(phi); ER(3) = sin(the) * sin(phi)
            proect = DOT_PRODUCT(vel * Time, ER)  !  ������� �������� ����������� �� ������ ������ ����
			KORD(1) = gl_x2(yzel, now); KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now) 
            R_TS = norm2(KORD + proect * ER)  ! ����� ���������� �� TS
            
            ! HP
            yzel = gl_RAY_B(par_n_HP, j, k)
            if(gl_Point_num(yzel) == 0) then
                vel = 0.0
			else
				vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
                vel = vel/gl_Point_num(yzel)                       ! ����� �������� �������� ����� ����
			end if
			
            
            ! ������� ��� ���������� �������������
            gl_Point_num(yzel) = 0
            gl_Vx(yzel) = 0.0
            gl_Vy(yzel) = 0.0
            gl_Vz(yzel) = 0.0
            
            proect = DOT_PRODUCT(vel * Time, ER)
			KORD(1) = gl_x2(yzel, now); KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now) 
            R_HP = norm2(KORD + proect * ER)  ! ����� ���������� �� HP
			
			
			
                
            do i = 1, N1

                if (i == 1) CYCLE
                
                yzel = gl_RAY_B(i, j, k)
                ! ��������� ���������� ����� �� ����

                ! �� TS
				if (i <= par_n_IB) then  ! NEW
                    r =  par_R0 + (par_R_inner - par_R0) * (DBLE(i)/(par_n_IB))**par_kk1
                else if (i <= par_n_TS) then  ! �� ���������� = R_TS
                    r =  par_R_inner + (R_TS - par_R_inner) * (DBLE(i - par_n_IB)/(par_n_TS - par_n_IB))**par_kk12
                !if (i <= par_n_TS) then  ! �� ���������� = R_TS
                !    r =  par_R0 + (R_TS - par_R0) * (REAL(i, KIND = 4)/par_n_TS)**par_kk1
                else if (i <= par_n_HP) then  ! �� ���������� = par_R_character * 1.3
                    r = R_TS + (i - par_n_TS) * (R_HP - R_TS) /(par_n_HP - par_n_TS)
				end if
				
				if (i == par_n_TS - 1) then
					r = R_TS - 0.5               
				end if
				

                ! ���������� ����� ����������
                gl_x2(yzel, now2) = r * cos(the)
                gl_y2(yzel, now2) = r * sin(the) * cos(phi)
                gl_z2(yzel, now2) = r * sin(the) * sin(phi)
                
               
            end do
        end do
	end do
    
	
    N3 = size(gl_RAY_C(1, 1, :))
    N2 = size(gl_RAY_C(1, :, 1))
    N1 = size(gl_RAY_C(:, 1, 1))

    ! ���� �������� ����� �� ����� C  ************************************************************
    do k = 1, N3
        do j = 1, N2
            
            ! ��������� ���������� �������� ���� � ������������
                phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
                
                ! ��������� ���������� ����� �� ����
                x = gl_x2(gl_RAY_B(par_n_HP, j, k), now2)
                y = gl_y2(gl_RAY_B(par_n_HP, j, k), now2)
                z = gl_z2(gl_RAY_B(par_n_HP, j, k), now2)
                rr = (y**2 + z**2)**(0.5)
				
				
				y = gl_y2(gl_RAY_B(par_n_HP - 1, j, k), now2)
                z = gl_z2(gl_RAY_B(par_n_HP - 1, j, k), now2)
                rd = (y**2 + z**2)**(0.5)
				rr = rr + (rr - rd)
                
                ! BS     ����� ����� ��������� BS �� � ��������� �� ������� ���� A
                yzel = gl_RAY_A(par_n_BS, size(gl_RAY_A(1, :, 1)), k)
				ER(1) = 0.0_8; ER(2) = gl_y2(yzel, now2); ER(3) = gl_z2(yzel, now2)
                R_BS = norm2(ER)  ! ����� ���������� �� BS
				
            
            do i = 1, N1
                
                if(i == 1) CYCLE
                
                yzel = gl_RAY_C(i, j, k)

                if(i == 2) then
					r = rr
				else if (i <= par_n_BS - par_n_HP + 1) then
                    r = rr + (i - 2) * (R_BS - rr)/(par_n_BS - par_n_HP - 1)
                else
                    r = R_BS + (DBLE(i - (par_n_BS - par_n_HP + 1))/(N1 - (par_n_BS - par_n_HP + 1) ))**(0.55 * par_kk2) * (par_R_END - R_BS)
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


    ! ���� �������� ����� �� ����� O   ************************************************************
    do k = 1, N3
        phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
        
        do j = 1, N2
            
            yzel = gl_RAY_O(1, j, k)
            if(gl_Point_num(yzel) == 0) then
                vel = 0.0
			else
				vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
                vel = vel/gl_Point_num(yzel)                       ! ����� �������� �������� ����� ����
            end if
            
            ! ������� ��� ���������� �������������
            gl_Point_num(yzel) = 0
            gl_Vx(yzel) = 0.0_8
            gl_Vy(yzel) = 0.0_8
            gl_Vz(yzel) = 0.0_8
            
			ER(1) = 0.0_8; ER(2) = cos(phi); ER(3) = sin(phi)
            proect = DOT_PRODUCT(vel * Time, ER)  !  ������� �������� ����������� �� ������ ������ ����
			KORD(1) = 0.0_8; KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now)
            R_HP = norm2(KORD + proect * ER )  ! ����� ���������� �� HP
			
			
			! ��������� ����������� �������� � ���
			if(R_HP < 20.0_8) then
				R_HP = 20.0_8
			end if
			
            
            xx = gl_x2(gl_RAY_B(par_n_HP, par_m_BC, k), now2)              ! ������������� �� x - ���������� ������� ����� B �� ���������� � ���� ��������� (k)
            x = xx - (DBLE(j)/N2)**par_kk31 * (xx - par_R_LEFT)
            
            ! BS     ����� ����� ��������� BS �� � ��������� �� ������� ���� A
            yzel = gl_RAY_A(par_n_BS, size(gl_RAY_A(1, :, 1)), k)
			KORD(1) = 0.0_8; KORD(2) = gl_y2(yzel, now2); KORD(3) = gl_z2(yzel, now2)
            R_BS = norm2(KORD)  ! ����� ���������� �� BS
			
			if (R_BS < R_HP + 10.0) then
				R_BS = R_HP + 10.0
				!print*, "R_BS = R_HP + 10.0"
			end if
			
            
            do i = 1, N1
                yzel = gl_RAY_O(i, j, k)
                
                if (i <= par_n_BS - par_n_HP + 1) then
                    r = R_HP + (i - 1) * (R_BS - R_HP)/(par_n_BS - par_n_HP)
                else
                    r = R_BS + (DBLE(i - (par_n_BS - par_n_HP + 1))/(N1 - (par_n_BS - par_n_HP + 1) ))**(0.55 * par_kk2) * (par_R_END - R_BS)
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

    ! ���� �������� ����� �� ����� K  ************************************************************
    do k = 1, N3
        do j = 1, N2
            
             ! ��������� ���������� �������� ���� � ������������
            the = par_pi_8/2.0 + par_triple_point + (N2 - j + 1) * (par_pi_8/2.0 - par_triple_point)/(N2)
            phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
            
            if (k /= 1 .and. j == 1) CYCLE
            
            yzel = gl_RAY_K(par_n_TS, j, k)
            if(gl_Point_num(yzel) == 0) then
                vel = 0.0
			else
				vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
                vel = vel/gl_Point_num(yzel)                       ! ����� �������� �������� ����� ����
            end if
            
            ! ������� ��� ���������� �������������
            gl_Point_num(yzel) = 0
            gl_Vx(yzel) = 0.0
            gl_Vy(yzel) = 0.0
            gl_Vz(yzel) = 0.0
            
			ER(1) = cos(the); ER(2) = sin(the) * cos(phi); ER(3) = sin(the) * sin(phi)
            proect = DOT_PRODUCT(vel * Time, ER)  !  ������� �������� ����������� �� ������ ������ ����
            KORD(1) = gl_x2(yzel, now); KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now)
			R_TS = norm2(KORD + proect * ER)  ! ����� ���������� �� TS
            
            do i = 1, N1

                if (i == 1) CYCLE
                
                ! ��������� ���������� ����� �� ����
                yzel = gl_RAY_K(i, j, k)
				
				if (i <= par_n_IB) then  ! NEW
                    r =  par_R0 + (par_R_inner - par_R0) * (DBLE(i)/(par_n_IB))**par_kk1
                else 
                    r =  par_R_inner + (R_TS - par_R_inner) * (DBLE(i - par_n_IB)/(par_n_TS - par_n_IB))**par_kk12
				end if
				
				if (i == par_n_TS - 1) then
					r = R_TS - 0.5              
				end if
					
                !r =  par_R0 + (R_TS - par_R0) * (REAL(i, KIND = 4)/par_n_TS)**par_kk1


                ! ���������� ����� ����������
                gl_x2(yzel, now2) = r * cos(the)
                gl_y2(yzel, now2) = r * sin(the) * cos(phi)
                gl_z2(yzel, now2) = r * sin(the) * sin(phi)
                
                
            end do
        end do
    end do
    
    
     N3 = size(gl_RAY_D(1, 1, :))
    N2 = size(gl_RAY_D(1, :, 1))
    N1 = size(gl_RAY_D(:, 1, 1))

    ! ���� �������� ����� �� ����� D ************************************************************
	
	
	
    do k = 1, N3
        do j = 1, N2
            
            ! ��������� ���������� �������� ���� � ������������
                phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
				
				 ! ��������� ���������� ����� �� ����

                if (j < N2) then
                    xx = gl_x2(gl_RAY_K(par_n_TS, j, k), now2)
                    y = gl_y2(gl_RAY_K(par_n_TS, j, k), now2)
                    z = gl_z2(gl_RAY_K(par_n_TS, j, k), now2)
                else
                    xx = gl_x2(gl_RAY_B(par_n_TS, par_m_BC, k), now2)
                    y = gl_y2(gl_RAY_B(par_n_TS, par_m_BC, k), now2)
                    z = gl_z2(gl_RAY_B(par_n_TS, par_m_BC, k), now2)
				end if
				
				!y = gl_y2(gl_RAY_B(par_n_TS, par_m_BC, k), now2)
				!z = gl_z2(gl_RAY_B(par_n_TS, par_m_BC, k), now2)
				
				r = sqrt(y**2 + z**2)
				
				
            do i = 1, N1

                if (i == 1) CYCLE

                if (k /= 1 .and. j == 1) CYCLE
				
				
				!y = gl_y2(gl_RAY_O(1, i - 1, k), now2)
				!z = gl_z2(gl_RAY_O(1, i - 1, k), now2)
				!
				!rr = sqrt(y**2 + z**2) * 0.4

                yzel = gl_RAY_D(i, j, k)
                
                
                x = xx + (DBLE(i - 1)/(N1 - 1))**par_kk3 * (par_R_LEFT - xx)
				
				!rrr = r * (DBLE(j - 1)/(N2 - 1))
				!rrr = r * (1.0 - (DBLE(i - 1)**2/(N1 - 1)**2)) + (rr * (DBLE(j - 1)/(N2 - 1))) * ( (DBLE(i - 1)**2) /(N1 - 1)**2)

				
				
                ! ���������� ����� ����������
                gl_x2(yzel, now2) = x
                gl_y2(yzel, now2) = r * cos(phi)
                gl_z2(yzel, now2) = r * sin(phi)
                
            end do
        end do
    end do
    
    N3 = size(gl_RAY_E(1, 1, :))
    N2 = size(gl_RAY_E(1, :, 1))
    N1 = size(gl_RAY_E(:, 1, 1))
    
    ! ���� �������� ����� �� ����� E  ************************************************************
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
    
    subroutine Get_center_move(n, CENTER, now)  ! �������� ������ - ����� ������
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
    
    subroutine calc_all_Gran_move(now)   ! ��������� ������� �������� � �������� ������ � ������� �����
    ! ������ ������� ��� ������� �� ��������� �����
    use STORAGE
    use GEO_PARAM
    implicit none
    
    integer, intent(in) :: now           ! ������ ���� ���������� ����� � ���� ���������� ����������� ��������� (������� � ������)
    integer :: Ngran, iter
    real(8) :: p(3, 4), Vol, D
    real(8) :: a(3), b(3), c(3), S, node1(3), node2(3)
    real(8) :: dist, di, gr_center(3)
    integer :: i, j, k, ll, grc, ii, now2

    !print*, "calc_all_Gran_move"
    now2 = mod(now, 2) + 1 
    ! ���� �� ������ - ������� ������� �����, � �������
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
        ! ��������� �� ���������� ������ ���������� ����� �����

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

        ! ����� ��������� ����� �����
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

        ! ����� ���� ��� ���������, ��������� �� ������������� �������!\

        
        call Get_center_move(gl_Gran_neighbour(1, iter), node1, now)
        node2 = (p(:,1) + p(:,2) + p(:,3) + p(:,4))/4.0
        node1 = node2 - node1
        
		if (.False.) then
        !if(DOT_PRODUCT(node1, c) < 0.0) then
            print*, c
            print*, p(1, 1), norm2( (/0.0_8, p(2, 1), p(3, 1) /) )
            print*, "ERROR 13123 move ", DOT_PRODUCT(node1, c)
            print*, iter
            pause
        end if

        ! ����� �������� ������� ����� � ������� � ����� ������!
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

    


    ! ������ ��������� ������ �����
    Ngran = size(gl_all_Cell(1,:))

    do  iter = 1, Ngran   ! ����������� �� ���� �������
        Vol = 0.0
        ll = 0
        dist = 10000.0 * par_R_character
        c = 0.0
        ! ��� ���������� ����������� ������ ������, ����� ������ ����� ��� (��������� ���� ���������)
        if (gl_all_Cell(1, iter) /= gl_all_Cell(5, iter)) then
            do j = 1, 8
                i = gl_all_Cell(j, iter)
                c = c + (/gl_x2(i, now), gl_y2(i, now), gl_z2(i, now)/)
            end do
            c = c/8.0
        else
            if (gl_all_Cell(1, iter) == gl_all_Cell(2, iter)) then
                if (gl_all_Cell(4, iter) == gl_all_Cell(8, iter)) then  ! ����� ������ = 4
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
                    if (gl_all_Cell(4, iter) == gl_all_Cell(5, iter)) then ! ����� ������ = 1
                        do j = 1,8   ! ������ ������ ���� 4
                            ! if (j == 4 .or. j == 5 .or. j == 6 .or. j == 8) CYCLE
                            i = gl_all_Cell(j, iter)
                            c = c + (/gl_x2(i, now), gl_y2(i, now), gl_z2(i, now)/)
                        end do
                        c = c/8.0
                    else
                        do j = 1,8         ! ������ ������ ���� 2
                            ! if (j == 5 .or. j == 6) CYCLE             ! ����������� ������ ����� ������ ������!
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
                    else            ! ������ ������ ���� 5
                        do j = 1,8
                            ! if (j == 5 .or. j == 8) CYCLE    ! ����������� ������ ����� ������ ������!
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
		
		

        ! ��������� ����� ������ ������ ������� ����� �������� �� ������ �����
        do j = 1, 6
            i = gl_Cell_gran(j, iter)   ! ���� �� ������� ��� ����� ������

            !if (iter == gl_Cell_C(1, size(gl_Cell_C(1,:,1)) ,1)) then
            !    print*, gl_Cell_gran(:, iter)
            !    pause
            !end if

            if (i == 0) CYCLE
            k = gl_all_Gran(1, i)  ! ����� ������� ���� �����
            b = (/gl_x2(k, now), gl_y2(k, now), gl_z2(k, now)/)
            a = gl_Gran_normal2(:,i, now)
            di = dabs(DOT_PRODUCT(c,a) - DOT_PRODUCT(a,b))  ! ���������� �� ����� �� ���������
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
	
	
	subroutine Surface_setup()
	! ��� ������� ����������� ���������� �� ������ ���� ��������
	! ����� ��������� ���������� ������ ������������������ ���������� �����
	! Variables
	use GEO_PARAM
	use STORAGE
	use Interpolate
	use My_func
	implicit none
	
	real(8) :: position(10 * par_l_phi, 500)
	integer(4) :: posit_int(10 * par_l_phi, 500)
	real(8) :: Left, Right, x, y, z, r, phi, phi2, dsc
	integer(4) :: i, j, k, num, gr, ip, im
	integer(4) :: zone, cell, number, now, now2, yzel
	real(8) :: F(9), F_mf(5, 4), normal(3)
	
	! ������� ���������� ��� ������� ��� ����������� ��������
    if (allocated(gl_Vx) == .False.) then   ! ���� �������� ����������� ������� ���, �� ����� �������� ������ ��� ������������ ������� � ������ ��������
        allocate(gl_Vx(size(gl_x)))
        allocate(gl_Vy(size(gl_x)))
        allocate(gl_Vz(size(gl_x)))
        allocate(gl_Point_num(size(gl_x)))
        allocate(gl_x2(size(gl_x), 2))
        allocate(gl_y2(size(gl_x), 2))
        allocate(gl_z2(size(gl_x), 2))
        allocate(gl_Cell_Volume2(size(gl_Cell_Volume), 2))
        allocate(gl_Gran_normal2(  size(gl_Gran_normal(:, 1)), size(gl_Gran_normal(1, :)), 2  ))
        allocate(gl_Gran_center2(  size(gl_Gran_center(:, 1)), size(gl_Gran_center(1, :)), 2  ))
        allocate(gl_Cell_center2(size(gl_Cell_center(:, 1)), size(gl_Cell_center(1, :)), 2))
        allocate(gl_Gran_square2(size(gl_Gran_square), 2))
        
        ! ��������� �������������
        gl_Vx = 0.0
        gl_Vy = 0.0
        gl_Vz = 0.0
        gl_Point_num = 0
        gl_x2(:, 1) = gl_x
        gl_x2(:, 2) = gl_x
        gl_y2(:, 1) = gl_y
        gl_y2(:, 2) = gl_y
        gl_z2(:, 1) = gl_z
        gl_z2(:, 2) = gl_z
        gl_Cell_Volume2(:, 1) = gl_Cell_Volume
        gl_Cell_Volume2(:, 2) = gl_Cell_Volume
        gl_Gran_normal2(:, :, 1) = gl_Gran_normal
        gl_Gran_normal2(:, :, 2) = gl_Gran_normal
        gl_Gran_center2(:, :, 1) = gl_Gran_center
        gl_Gran_center2(:, :, 2) = gl_Gran_center
        gl_Cell_center2(:, :, 1) = gl_Cell_center
        gl_Cell_center2(:, :, 2) = gl_Cell_center
        gl_Gran_square2(:, 1) = gl_Gran_square
        gl_Gran_square2(:, 2) = gl_Gran_square
	end if
	
	number = 1
	cell = 1
	
	Left = -450.0_8
	Right = 50.0_8
	position = 0.0_8
	posit_int = 0
	
	do i = 1, 60
		print*, "i = ", i
		do k = 1, 10 * par_l_phi
			number = 1
			x = 25.0
			r = i * 0.2_8
			phi = 2.0 * par_pi_8 * k/(10 * par_l_phi) - par_pi_8/(10 * par_l_phi)
			y = r * cos(phi)
			z = r * sin(phi)
			
			do while (x > -450.0)
				!print*, "x = ", x
				call Interpolate_point(x, y, z, F, F_mf, zone, cell, number)
				if (cell == -1) EXIT
				number = cell
				x = x + 0.01 * F(2)
				y = y + 0.01 * F(3)
				z = z + 0.01 * F(4)
				
				phi2 = polar_angle(y, z)
				
				!print*, phi2, sqrt(y**2 + z**2)
				!pause
				
				position( INT(phi2/(2.0 * par_pi_8/(10 * par_l_phi))) + 1, INT(x + 450.01)) = position( INT(phi2/(2.0 * par_pi_8/(10 * par_l_phi))) + 1, INT(x + 450.01)) + sqrt(y**2 + z**2)
				posit_int(INT(phi2/(2.0 * par_pi_8/(10 * par_l_phi))) + 1, INT(x + 450.01))= posit_int(INT(phi2/(2.0 * par_pi_8/(10 * par_l_phi))) + 1, INT(x + 450.01)) + 1
			end do
			
		end do	
	end do
	
	
	do j = 1, 500
		do i = 1, 10 * par_l_phi
			if(posit_int(i, j) > 0) then
				position(i, j) = position(i, j)/posit_int(i, j)
			end if
			
		end do
	end do
	
	
	do k = 1, 20
		do j = 1, 500
			do i = 1, 10 * par_l_phi
				ip = i + 1
				im = i - 1
				if(i == 1) im = 10 * par_l_phi
				if(i == 10 * par_l_phi) ip = 1
			
				if(position(i, j) <= 0.001) then
					if(position(ip, j) > 0.001 .and. position(im, j) > 0.001) then
						position(i, j) = 0.5 * (position(ip, j) + position(im, j))
					else if (position(ip, j) < 0.001 .and. position(im, j) > 0.001) then
						position(i, j) = position(im, j)
					else if (position(im, j) < 0.001 .and. position(ip, j) > 0.001) then
						position(i, j) = position(ip, j)
					end if
				end if
			
			end do
		end do
	end do
	
	do j = 1, 500
			do i = 1, 10 * par_l_phi
				ip = i + 1
				im = i - 1
				if(i == 1) im = 10 * par_l_phi
				if(i == 10 * par_l_phi) ip = 1
			
				position(i, j) = (position(ip, j) + position(i, j) + position(im, j))/3.0
			
			end do
		end do
	
	
	
	open(2, file = 'save_SERF.txt')
	
	do j = 1, 500
		do i = 1, 10 * par_l_phi
			x = Left + j
			r = position(i, j)
			phi = i * (2.0 * par_pi_8/(10 * par_l_phi))
			y = r * cos(phi)
			z = r * sin(phi)
			write(2,*) x, y, z
		end do
	end do
	
	close(2)
	
	
	
	! ������ ����������� �����������
	
	! �������
    Num = size(gl_Contact)
	now = 2
	
	do k = 1, 600
		if( mod(k, 10) == 0) print*, "k = ", k, "  from 600"
		now2 = now
		now = mod(now, 2) + 1
	
	
		gl_Vx = 0.0
		gl_Vy = 0.0
		gl_Vz = 0.0
		gl_Point_num = 0
    
		do i = 1, Num
			!print*, "i = ", i
			gr = gl_Contact(i)
		
			x = gl_Gran_center2(1, gr, now)
			y = gl_Gran_center2(2, gr, now)
			z = gl_Gran_center2(3, gr, now)
			phi2 = polar_angle(y, z)
		
			if (x > -50) CYCLE
		
			!print*, INT(phi2/(2.0 * par_pi_8/par_l_phi)) + 1, INT(x + 450.01)
			if (sqrt(y**2 + z**2) < position( INT(phi2/(2.0 * par_pi_8/(10 * par_l_phi))) + 1, INT(x + 450.01))) then
				dsc = 0.1
			else
				dsc = -0.1
			end if
		
		
		
			normal = gl_Gran_normal2(:, gr, now)
		
			do j = 1, 4
				yzel = gl_all_Gran(j, gr)
			
				gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsc
				gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsc
				gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsc
				gl_Point_num(yzel) = gl_Point_num(yzel) + 1
			end do
		
		end do
		!print*, "Move_all "
		call Move_all(now, 1.0_8) 
		call calc_all_Gran_move(now2)
	
	
	end do
	
	gl_x = gl_x2(:, now2)
    gl_y = gl_y2(:, now2)
    gl_z = gl_z2(:, now2)
    gl_Cell_Volume = gl_Cell_Volume2(:, now2)
    gl_Gran_normal = gl_Gran_normal2(:, :, now2)
    gl_Gran_center = gl_Gran_center2(:, :, now2)
    gl_Cell_center = gl_Cell_center2(:, :, now2)
    gl_Gran_square = gl_Gran_square2(:, now2)
	
	
	! Body of Surface_setup
	
	end subroutine Surface_setup
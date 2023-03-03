
    
    subroutine Find_Surface()   ! ��������� ������� ������������, ������� ����������
    ! ��� BS ��������� ������ ����������� �� ���������! (x > 0)
    use STORAGE
    use GEO_PARAM
    implicit none
    integer :: j, k, node, num

    ! Body of Find_Surface
    node = 1

    do k = 1, size( gl_Cell_A(par_n_HP - 1, 1, :) )
        do j = 1, size( gl_Cell_A(par_n_HP - 1, :, 1) )
            gl_Contact(node) = gl_Cell_gran(1, gl_Cell_A(par_n_HP - 1, j, k))
            node = node + 1
        end do
    end do

    do k = 1, size( gl_Cell_C(par_n_HP - par_n_TS, 1, :) )
        do j = 1, size( gl_Cell_C(par_n_HP - par_n_TS, :, 1) )
            gl_Contact(node) = gl_Cell_gran(1, gl_Cell_C(par_n_HP - par_n_TS, j, k))
            node = node + 1
        end do
    end do

    node = 1
    ! TS
    do k = 1, size( gl_Cell_A(par_n_TS - 1, 1, :) )
        do j = 1, size( gl_Cell_A(par_n_TS - 1, :, 1) )
            gl_TS(node) = gl_Cell_gran(1, gl_Cell_A(par_n_TS - 1, j, k))
            node = node + 1
        end do
    end do

    do k = 1, size( gl_Cell_B(par_n_TS - 1, 1, :) )
        do j = 1, size( gl_Cell_B(par_n_TS - 1, :, 1) )
            gl_TS(node) = gl_Cell_gran(1, gl_Cell_B(par_n_TS - 1, j, k))
            node = node + 1
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

    end subroutine Find_Surface
    
    subroutine Calc_move(now) ! ��������� �������� �������� ����� (��������) � ������� ������� �� ������������ �������(
    use STORAGE
    use GEO_PARAM
    implicit none
    
    integer, intent(in) :: now   ! ������ ����� ���������� ������������ � ���������� ��������
    
    integer :: i, Num, gr, s1, s2, j, yzel
    real(8) :: normal(3), qqq1(8), qqq2(8), dsl, dsp, dsc, POTOK(8), a1, a2, v1, v2, &
        ray(3), norm, b1, b2, c1, c2, koef1, koef2, koef3
    
    koef1 = 0.3
    koef2 = 0.3
    koef3 = 0.3
    
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
        
        call chlld(2, normal(1), normal(2), normal(3), &
                0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)
        
        dsl = dsl * koef1
        
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
        
        call chlld(2, normal(1), normal(2), normal(3), &
                0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)
        
        dsc = dsc * koef2
        
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
    Num = size(gl_BS)
    
    do i = 1, Num
        
        gr = gl_BS(i)
    	s1 = gl_Gran_neighbour(1, gr)
        s2 = gl_Gran_neighbour(2, gr)
        normal = gl_Gran_normal2(:, gr, now)
        qqq1 = gl_Cell_par(1:8, s1)
        qqq2 = gl_Cell_par(1:8, s2)
        
        call chlld(2, normal(1), normal(2), normal(3), &
                0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)
        
        dsp = dsp * koef3
        
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
    
    end subroutine Calc_move   
    
    subroutine Move_all(now, Time)  ! ����������� ����� � ������������ �� ��������� �������� ����� �� ������������ �������
    use STORAGE
    use GEO_PARAM
    implicit none
    
    integer :: now2             ! ��� ��������� �� ������ ������ �� ������ now
    real(8), intent(in) :: Time
    integer, intent(in) :: now
    real(8) :: R_TS, proect, vel(3), R_HP, R_BS
    integer :: yzel, N1, N2, N3, i, j, k
    real(8) :: the, phi, r, x, y, z, rr, xx, x2, y2, z2
    
    now2 = mod(now, 2) + 1   
    
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
            phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
            
            ! TS
            yzel = gl_RAY_A(par_n_TS, j, k)
            if(gl_Point_num(yzel) == 0) then
                vel = 0.0
            else
                vel = (/gl_Vx(yzel), gl_Vy(yzel), gl_Vz(yzel)/)
                vel = vel/gl_Point_num(yzel)                       ! ����� �������� �������� ����� ����
            end if
            
            ! ������� ��� ���������� �������������
            gl_Point_num(yzel) = 0
            gl_Vx(yzel) = 0.0
            gl_Vy(yzel) = 0.0
            gl_Vz(yzel) = 0.0
            
            proect = DOT_PRODUCT(vel * Time, (/cos(the), sin(the) * cos(phi), sin(the) * sin(phi)/))  !  ������� �������� ����������� �� ������ ������ ����
            R_TS = norm2((/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/) + &
                    proect * (/cos(the), sin(the) * cos(phi), sin(the) * sin(phi)/))  ! ����� ���������� �� TS
            
            ! HP
            yzel = gl_RAY_A(par_n_HP, j, k)
            if(gl_Point_num(yzel) == 0) then
                vel = 0.0
            else
                vel = (/gl_Vx(yzel), gl_Vy(yzel), gl_Vz(yzel)/)
                vel = vel/gl_Point_num(yzel)                       ! ����� �������� �������� ����� ����
            end if
            
            ! ������� ��� ���������� �������������
            gl_Point_num(yzel) = 0
            gl_Vx(yzel) = 0.0
            gl_Vy(yzel) = 0.0
            gl_Vz(yzel) = 0.0
            
            proect = DOT_PRODUCT(vel * Time, (/cos(the), sin(the) * cos(phi), sin(the) * sin(phi)/))
            R_HP = norm2((/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/) + &
                    proect * (/cos(the), sin(the) * cos(phi), sin(the) * sin(phi)/))  ! ����� ���������� �� HP
            
            ! BS
            yzel = gl_RAY_A(par_n_BS, j, k)
            if(gl_Point_num(yzel) == 0) then
                vel = 0.0
            else
                vel = (/gl_Vx(yzel), gl_Vy(yzel), gl_Vz(yzel)/)
                vel = vel/gl_Point_num(yzel)                       ! ����� �������� �������� ����� ����
            end if
            
            ! ������� ��� ���������� �������������
            gl_Point_num(yzel) = 0
            gl_Vx(yzel) = 0.0
            gl_Vy(yzel) = 0.0
            gl_Vz(yzel) = 0.0
            
            proect = DOT_PRODUCT(vel * Time, (/cos(the), sin(the) * cos(phi), sin(the) * sin(phi)/))
            R_BS = norm2((/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/) + &
                    proect * (/cos(the), sin(the) * cos(phi), sin(the) * sin(phi)/))  ! ����� ���������� �� BS
            
            ! ����� ������� ���� ���������� ��������� �����, ����� ��, ��� � ��� ���������� �����
            do i = 1, N1
                
                if (i == 1) then
                    CYCLE
                end if

                yzel = gl_RAY_A(i, j, k)
                ! ��������� ���������� ����� �� ����

                ! �� TS
                if (i <= par_n_TS) then  ! �� ���������� = R_TS
                    r =  par_R0 + (R_TS - par_R0) * (DBLE(i)/par_n_TS)**par_kk1
                    !print *, r
                    !pause
                else if (i <= par_n_HP) then  
                    r = R_TS + (i - par_n_TS) * (R_HP - R_TS)/(par_n_HP - par_n_TS)
                else if (i <= par_n_BS) then 
                    r = R_HP + (i - par_n_HP) * (R_BS - R_HP)/(par_n_BS - par_n_HP)
                else
                    r = R_BS + (par_R_END - R_BS) * (DBLE(i- par_n_BS)/(par_n_END - par_n_BS))**par_kk2
                end if

                ! ���������� ����� ����������
                gl_x2(yzel, now2) = r * cos(the)
                gl_y2(yzel, now2) = r * sin(the) * cos(phi)
                gl_z2(yzel, now2) = r * sin(the) * sin(phi)
                
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
            phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
            
            ! TS
            yzel = gl_RAY_B(par_n_TS, j, k)
            if(gl_Point_num(yzel) == 0) then
                vel = 0.0
            else
                vel = (/gl_Vx(yzel), gl_Vy(yzel), gl_Vz(yzel)/)
                vel = vel/gl_Point_num(yzel)                       ! ����� �������� �������� ����� ����
            end if
            
            ! ������� ��� ���������� �������������
            gl_Point_num(yzel) = 0
            gl_Vx(yzel) = 0.0
            gl_Vy(yzel) = 0.0
            gl_Vz(yzel) = 0.0
            
            proect = DOT_PRODUCT(vel * Time, (/cos(the), sin(the) * cos(phi), sin(the) * sin(phi)/))  !  ������� �������� ����������� �� ������ ������ ����
            R_TS = norm2((/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/) + &
                    proect * (/cos(the), sin(the) * cos(phi), sin(the) * sin(phi)/))  ! ����� ���������� �� TS
            
            ! HP
            yzel = gl_RAY_B(par_n_HP, j, k)
            if(gl_Point_num(yzel) == 0) then
                vel = 0.0
            else
                vel = (/gl_Vx(yzel), gl_Vy(yzel), gl_Vz(yzel)/)
                vel = vel/gl_Point_num(yzel)                       ! ����� �������� �������� ����� ����
            end if
            
            ! ������� ��� ���������� �������������
            gl_Point_num(yzel) = 0
            gl_Vx(yzel) = 0.0
            gl_Vy(yzel) = 0.0
            gl_Vz(yzel) = 0.0
            
            proect = DOT_PRODUCT(vel * Time, (/cos(the), sin(the) * cos(phi), sin(the) * sin(phi)/))
            R_HP = norm2((/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/) + &
                    proect * (/cos(the), sin(the) * cos(phi), sin(the) * sin(phi)/))  ! ����� ���������� �� HP
                
            do i = 1, N1

                if (i == 1) CYCLE
                
                yzel = gl_RAY_B(i, j, k)
                ! ��������� ���������� ����� �� ����

                ! �� TS
                if (i <= par_n_TS) then  ! �� ���������� = R_TS
                    r =  par_R0 + (R_TS - par_R0) * (REAL(i, KIND = 4)/par_n_TS)**par_kk1
                else if (i <= par_n_HP) then  ! �� ���������� = par_R_character * 1.3
                    r = R_TS + (i - par_n_TS) * (R_HP - R_TS) /(par_n_HP - par_n_TS)
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
                phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
                
                ! ��������� ���������� ����� �� ����
                x = gl_x2(gl_RAY_B(par_n_HP, j, k), now2)
                y = gl_y2(gl_RAY_B(par_n_HP, j, k), now2)
                z = gl_z2(gl_RAY_B(par_n_HP, j, k), now2)
                rr = (y**2 + z**2)**(0.5)
                
                ! BS     ����� ����� ��������� BS �� � ��������� �� ������� ���� A
                yzel = gl_RAY_A(par_n_BS, size(gl_RAY_A(1, :, 1)), k)
                R_BS = norm2((/0.0_8, gl_y2(yzel, now2), gl_z2(yzel, now2)/))  ! ����� ���������� �� BS
            
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


    ! ���� �������� ����� �� ����� O   ************************************************************
    do k = 1, N3
        phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
        
        do j = 1, N2
            
            yzel = gl_RAY_O(1, j, k)
            if(gl_Point_num(yzel) == 0) then
                vel = 0.0
            else
                vel = (/gl_Vx(yzel), gl_Vy(yzel), gl_Vz(yzel)/)
                vel = vel/gl_Point_num(yzel)                       ! ����� �������� �������� ����� ����
            end if
            
            ! ������� ��� ���������� �������������
            gl_Point_num(yzel) = 0
            gl_Vx(yzel) = 0.0_8
            gl_Vy(yzel) = 0.0_8
            gl_Vz(yzel) = 0.0_8
            
            proect = DOT_PRODUCT(vel * Time, (/0.0_8, cos(phi), sin(phi)/) )  !  ������� �������� ����������� �� ������ ������ ����
            R_HP = norm2((/0.0_8, gl_y2(yzel, now), gl_z2(yzel, now)/) + &
                    proect * (/0.0_8, cos(phi), sin(phi)/) )  ! ����� ���������� �� HP
            
            xx = gl_x2(gl_RAY_B(par_n_HP, par_m_BC, k), now2)              ! ������������� �� x - ���������� ������� ����� B �� ���������� � ���� ��������� (k)
            x = xx - (DBLE(j)/N2)**par_kk3 * (xx - par_R_LEFT)
            
            ! BS     ����� ����� ��������� BS �� � ��������� �� ������� ���� A
            yzel = gl_RAY_A(par_n_BS, size(gl_RAY_A(1, :, 1)), k)
            R_BS = norm2((/0.0_8, gl_y2(yzel, now2), gl_z2(yzel, now2)/))  ! ����� ���������� �� BS
            
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

    ! ���� �������� ����� �� ����� K  ************************************************************
    do k = 1, N3
        do j = 1, N2
            
             ! ��������� ���������� �������� ���� � ������������
            the = par_pi_8/2.0 + par_triple_point + (N2 - j + 1) * (par_pi_8/2.0 - par_triple_point)/(N2)
            phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
            
            if (k /= 1 .and. j == 1) CYCLE
            
            yzel = gl_RAY_K(par_n_TS, j, k)
            if(gl_Point_num(yzel) == 0) then
                vel = 0.0
            else
                vel = (/gl_Vx(yzel), gl_Vy(yzel), gl_Vz(yzel)/)
                vel = vel/gl_Point_num(yzel)                       ! ����� �������� �������� ����� ����
            end if
            
            ! ������� ��� ���������� �������������
            gl_Point_num(yzel) = 0
            gl_Vx(yzel) = 0.0
            gl_Vy(yzel) = 0.0
            gl_Vz(yzel) = 0.0
            
            proect = DOT_PRODUCT(vel * Time, (/cos(the), sin(the) * cos(phi), sin(the) * sin(phi)/))  !  ������� �������� ����������� �� ������ ������ ����
            R_TS = norm2((/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/) + &
                    proect * (/cos(the), sin(the) * cos(phi), sin(the) * sin(phi)/))  ! ����� ���������� �� TS
            
            do i = 1, N1

                if (i == 1) CYCLE
                
                ! ��������� ���������� ����� �� ����
                yzel = gl_RAY_K(i, j, k)
                r =  par_R0 + (R_TS - par_R0) * (REAL(i, KIND = 4)/par_n_TS)**par_kk1


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
                phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
                
            do i = 1, N1

                if (i == 1) CYCLE

                if (k /= 1 .and. j == 1) CYCLE

                
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

                yzel = gl_RAY_D(i, j, k)
                r = sqrt(y**2 + z**2)
                
                x = xx + (DBLE(i - 1)/(N1 - 1))**par_kk3 * (par_R_LEFT - xx)

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
                
            end do
        end do
    end do

    
    
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
    integer, automatic :: Ngran, iter
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
        
        if(DOT_PRODUCT(node1, c) < 0.0) then
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
                if (gl_all_Cell(4, iter) == gl_all_Cell(8, iter)) then
                    do j = 1,8
                        if (j == 2 .or. j == 5 .or. j == 6 .or. j == 8) CYCLE
                        i = gl_all_Cell(j, iter)
                        c = c + (/gl_x2(i, now), gl_y2(i, now), gl_z2(i, now)/)
                    end do
                    c = c/4.0
                else
                    do j = 1,8
                        if (j == 2 .or. j == 5 .or. j == 6) CYCLE
                        i = gl_all_Cell(j, iter)
                        c = c + (/gl_x2(i, now), gl_y2(i, now), gl_z2(i, now)/)
                    end do
                    c = c/5.0
                end if
            else
                if (gl_all_Cell(2, iter) == gl_all_Cell(6, iter)) then
                    if (gl_all_Cell(4, iter) == gl_all_Cell(5, iter)) then
                        do j = 1,8
                            if (j == 4 .or. j == 5 .or. j == 6 .or. j == 8) CYCLE
                            i = gl_all_Cell(j, iter)
                            c = c + (/gl_x2(i, now), gl_y2(i, now), gl_z2(i, now)/)
                        end do
                        c = c/4.0
                    else
                        do j = 1,8
                            if (j == 5 .or. j == 6) CYCLE
                            i = gl_all_Cell(j, iter)
                            c = c + (/gl_x2(i, now), gl_y2(i, now), gl_z2(i, now)/)
                        end do
                        c = c/6.0
                    end if
                else
                    if (gl_all_Cell(4, iter) == gl_all_Cell(5, iter)) then
                        do j = 1,8
                            if (j == 4 .or. j == 5 .or. j == 8) CYCLE
                            i = gl_all_Cell(j, iter)
                            c = c + (/gl_x2(i, now), gl_y2(i, now), gl_z2(i, now)/)
                        end do
                        c = c/5.0
                    else
                        do j = 1,8
                            if (j == 5 .or. j == 8) CYCLE
                            i = gl_all_Cell(j, iter)
                            c = c + (/gl_x2(i, now), gl_y2(i, now), gl_z2(i, now)/)
                        end do
                        c = c/6.0
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
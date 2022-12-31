!  MIK.f90 
    
!****************************************************************************
!
!  PROGRAM: MIK model - Malama & Izmodenov & Korolkov model
!
!****************************************************************************
    ! �������� �������
    
include "Storage_modul.f90"    
include "Solvers.f90"  
include "Help_func.f90"


    
    module My_func                     ! ������ ����������� ��� ������� �������
    
    interface
        real(8) pure function calc_square(n)
            integer, intent(in) :: n
        end function
        
        real(8) pure function polar_angle(x, y)
            use GEO_PARAM
	        implicit none
            real(8), intent(in) :: x, y
        end function
        
    end interface
    
    end module My_func
    
    ! ��������������� �������
    
    real(8) pure function polar_angle(x, y)
    use GEO_PARAM
	implicit none
    real(8), intent(in) :: x, y
    
	if (abs(x) + abs(y) < 0.00001 / par_R_character) then
        polar_angle = 0.0_8
        return 
    end if
    

	if (x < 0) then
		polar_angle = atan(y / x) + 1.0 * par_pi_8
	elseif (x > 0 .and. y >= 0) then
		polar_angle = atan(y / x)
	elseif (x > 0 .and. y < 0) then
		polar_angle = atan(y / x) + 2.0 * par_pi_8
	elseif (y > 0 .and. x >= 0 .and. x <= 0) then
		polar_angle = par_pi_8 / 2.0
	elseif (y < 0 .and. x >= 0 .and. x <= 0) then
		polar_angle =  3.0 * par_pi_8 / 2.0
	end if
        
	return 
    end function polar_angle

    subroutine spherical_skorost(z, x, y, Vz, Vx, Vy, Vr, Vphi, Vtheta)
    ! Variables
    use My_func
    implicit none 
    real(8), intent(in) :: x, y, z, Vx, Vy, Vz
    real(8), intent(out) :: Vr, Vphi, Vtheta
    real(8) :: r_1, the_1, phi_1
    
    r_1 = sqrt(x * x + y * y + z * z)
	the_1 = acos(z / r_1)
	phi_1 = polar_angle(x, y)

	Vr = Vx * sin(the_1) * cos(phi_1) + Vy * sin(the_1) * sin(phi_1) + Vz * cos(the_1);
	Vtheta = Vx * cos(the_1) * cos(phi_1) + Vy * cos(the_1) * sin(phi_1) - Vz * sin(the_1);
	Vphi = -Vx * sin(phi_1) + Vy * cos(phi_1);
    
    end subroutine spherical_skorost
    
    subroutine dekard_skorost(z, x, y, Vr, Vphi, Vtheta, Vz, Vx, Vy)
    use My_func
    implicit none 
    real(8), intent(in) :: x, y, z,  Vr, Vphi, Vtheta
    real(8), intent(out) :: Vx, Vy, Vz
    real(8) :: r_2, the_2, phi_2
    
    r_2 = sqrt(x * x + y * y + z * z);
	the_2 = acos(z / r_2);
	phi_2 = polar_angle(x, y);

	Vx = Vr * sin(the_2) * cos(phi_2) + Vtheta * cos(the_2) * cos(phi_2) - Vphi * sin(phi_2);
	Vy = Vr * sin(the_2) * sin(phi_2) + Vtheta * cos(the_2) * sin(phi_2) + Vphi * cos(phi_2);
	Vz = Vr * cos(the_2) - Vtheta * sin(the_2);
    
    end subroutine dekard_skorost
    
    ! ****************************************************************************************************************************************************
    ! ���� ������� ��� ���������� ���������� �����
   
    subroutine Set_STORAGE()             ! ������� �������� ������ ��� ���� ��������� ������ STORAGE, ��������� ��������� �� ������ GEO_PARAM
    
    use STORAGE
    use GEO_PARAM
    implicit none
    
    integer :: n1, n2
    
    if (par_developer_info) print *, "START Set_STORAGE"
    ! �������� ������ ��� ����������
    !gl_x = [real(8) ::]
    !gl_y = [real(8) ::]
    !gl_z = [real(8) ::]
    
    
    allocate(gl_RAY_A(par_n_END, par_m_A, par_l_phi))
    allocate(gl_RAY_B(par_n_HP, par_m_BC, par_l_phi))
    allocate(gl_RAY_C(par_n_END - par_n_HP + 1, par_m_BC, par_l_phi))
    allocate(gl_RAY_O(par_n_END - par_n_HP + 1, par_m_O, par_l_phi))
    allocate(gl_RAY_K(par_n_TS, par_m_K, par_l_phi))
    allocate(gl_RAY_D(par_m_O + 1, par_m_K + 1, par_l_phi))
    allocate(gl_RAY_E(par_n_HP - par_n_TS + 1, par_m_O, par_l_phi))
    
    allocate(gl_Cell_A(par_n_END - 1, par_m_A + par_m_BC - 1, par_l_phi))
    allocate(gl_Cell_B( (par_n_TS - 1) + par_m_O, par_m_K, par_l_phi) )
    allocate(gl_Cell_C( par_n_END - par_n_TS, par_m_O, par_l_phi ))
    
    
    allocate(gl_all_Cell(8, size(gl_Cell_A(:,:,:)) + size(gl_Cell_B(:,:,:)) + size(gl_Cell_C(:,:,:)) ) )
    allocate(gl_Cell_neighbour(6, size(gl_Cell_A(:,:,:)) + size(gl_Cell_B(:,:,:)) + size(gl_Cell_C(:,:,:)) ) )
    allocate(gl_Cell_gran(6, size(gl_Cell_A(:,:,:)) + size(gl_Cell_B(:,:,:)) + size(gl_Cell_C(:,:,:))))
    
    allocate(gl_Cell_Volume(size(gl_Cell_A(:,:,:)) + size(gl_Cell_B(:,:,:)) + size(gl_Cell_C(:,:,:))))
    allocate(gl_Cell_dist(size(gl_Cell_A(:,:,:)) + size(gl_Cell_B(:,:,:)) + size(gl_Cell_C(:,:,:))))
    allocate( gl_Cell_center(3, size(gl_Cell_A(:,:,:)) + size(gl_Cell_B(:,:,:)) + size(gl_Cell_C(:,:,:)) ) )
    
    allocate( gl_Cell_par(8, size(gl_Cell_A(:,:,:)) + size(gl_Cell_B(:,:,:)) + size(gl_Cell_C(:,:,:)) ) )
    
    ! ��������� ����� ����� � �����
    par_n_points = par_n_END * par_l_phi * (par_m_A + par_m_BC) + par_m_K * (par_n_TS + par_m_O) * par_l_phi + par_l_phi * (par_n_END - par_n_TS + 1) * par_m_O - &
            (par_m_A + par_m_BC + par_m_K - 1) * par_l_phi - par_n_END * (par_l_phi - 1) - (par_n_TS + par_m_O - 1) * (par_l_phi - 1)  ! ����� ����� � �����
    
    allocate(gl_x(par_n_points))
    allocate(gl_y(par_n_points))
    allocate(gl_z(par_n_points))
    !gl_x = [real(8) ::]
    !gl_y = [real(8) ::]
    !gl_z = [real(8) ::]
    
    n2 =  (par_n_END - 1) * (par_m_A + par_m_BC - 1) + (par_n_TS - 1 + par_m_O) * (par_m_K) + &
        (par_n_END - par_n_TS) * par_m_O ! ����� ����� � 1 ���� �� ����
    
    n1 = 2 * (par_n_END - 1) * (par_m_A + par_m_BC - 1) - (par_n_END - 1) + 2 * (par_n_TS - 1 + par_m_O) * (par_m_K) + &
        2 * (par_n_END - par_n_TS) * par_m_O + (par_n_END - par_n_TS)
    n1 = (n1 + n2) * par_l_phi
    
    if (par_developer_info)print*, "1 sloy CELL ", n2 , " == ", size(gl_Cell_A(:,:,1)) + size(gl_Cell_B(:,:,1)) + size(gl_Cell_C(:,:,1))
    
    allocate(gl_all_Gran(4, n1))
    allocate(gl_Gran_neighbour(2, n1))
    allocate(gl_Gran_normal(3, n1))
    allocate(gl_Gran_square(n1))
    allocate(gl_Gran_POTOK(8, n1))
    allocate(gl_Gran_center(3, n1))
    
    allocate(gl_Contact( (par_m_O + par_m_A + par_m_BC -1) * par_l_phi ))   ! �������� ������ ��� �������
    
    if (par_developer_info) print*, "Set_STORAGE: Ray A contains points: ", size(gl_RAY_A(:,1,1))
    if (par_developer_info) print*, "Set_STORAGE: Ray B contains points: ", size(gl_RAY_B(:,1,1))
    if (par_developer_info) print*, "Set_STORAGE: Ray C contains points: ", size(gl_RAY_C(:,1,1))
    if (par_developer_info) print*, "Set_STORAGE: Ray O contains points: ", size(gl_RAY_O(:,1,1))
    if (par_developer_info) print*, "Set_STORAGE: Ray K contains points: ", size(gl_RAY_K(:,1,1))
    if (par_developer_info) print*, "Set_STORAGE: Ray D contains points: ", size(gl_RAY_D(:,1,1))
    if (par_developer_info) print*, "Set_STORAGE: Ray E contains points: ", size(gl_RAY_E(:,1,1))
    
    if (par_developer_info) print*, "Set_STORAGE: all GRANS: ", n1
    
    if (par_developer_info) print*, "Set_STORAGE: all CELL: ", size(gl_all_Cell(1,:))
    if (par_developer_info) print*, "Set_STORAGE: 1 sloy CELL: ", n2
    
    
    ! ��������� ���������� ����������
    gl_x = 0.0
    gl_y = 0.0
    gl_z = 0.0
    gl_Cell_Volume = 0.0
    gl_RAY_A = -1
    gl_RAY_B = -1
    gl_RAY_C = -1
    gl_RAY_O = -1
    gl_RAY_K = -1
    gl_RAY_D = -1
    gl_RAY_E = -1
    
    gl_Cell_dist = 0.0
    gl_Cell_center = 0.0
    gl_Gran_POTOK = 0.0
    gl_Cell_par = 0.0
    gl_Gran_center = 0.0
    
    gl_Cell_A = -1
    gl_Cell_B = -1
    gl_Cell_C = -1
    
    gl_all_Cell = -1
    gl_Cell_neighbour = 0          ! 0 - ��� ������ (��� ���� ��� �� ������� - � ������� ���� �����)
    gl_Cell_gran = 0
    
    gl_all_Gran = 0
    gl_Gran_neighbour = 0
    
    gl_Contact = 0
    
    gl_Gran_normal = 0.0;
    gl_Gran_square = 0.0;
    
    if (par_developer_info) print *, "END Set_STORAGE"
    
    end subroutine Set_STORAGE

    subroutine Build_Mesh_start()        ! ��������� ���������� ����� 
    
    use STORAGE
    use GEO_PARAM
    implicit none
    
    integer(4), automatic :: i, j, k, N1, N2, N3, i1, kk, node, kk2, ni
    real(8), automatic :: r, phi, the, xx, x, y, z, rr, x2, y2, z2
    
    ! ����� �������� � ���������� ������� ���������, ���� ���� ����� ������������� ������������ 
    ! ����� ����� ��������������� x,y,z <---> x, r, phi  (x - ��� ���������), � ����� ����������� ����������
    
    node = 1
    ! ����� ������ ����
    gl_x(node) = 0.0
    gl_y(node) = 0.0
    gl_z(node) = 0.0
    node = node + 1
    !gl_x = [gl_x, 0.0]
    !gl_y = [gl_y, 0.0]
    !gl_z = [gl_z, 0.0]
    
    N3 = size(gl_RAY_A(1, 1, :))
    N2 = size(gl_RAY_A(1, :, 1))
    N1 = size(gl_RAY_A(:, 1, 1))
    
    ! **************************************************************** ������ � ����� �����
    
    ! ���� ��������� ����� �� ����� � � �� ���������� � ����� ������ ************************************************************
    do k = 1, N3
    	do j = 1, N2
            do i = 1, N1
                
                if (i == 1) then
                    gl_RAY_A(i, j, k) = 1
                    CYCLE
                end if
                
                if (k /= 1 .and. j == 1) then
                    gl_RAY_A(i, j, k) = i
                    CYCLE
                end if
                
    	        ! ��������� ���������� �������� ���� � ������������
                the = (j - 1) * par_pi_8/2.0/(N2 - 1)
                phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
                ! ��������� ���������� ����� �� ����
                
                ! �� TS
                if (i <= par_n_TS) then  ! �� ���������� = par_R_character
                    r =  par_R0 + (par_R_character - par_R0) * (DBLE(i)/par_n_TS)**par_kk1
                    !print *, r
                    !pause
                else if (i <= par_n_HP) then  ! �� ���������� = par_R_character * 1.3
                    r = par_R_character + (i - par_n_TS) * 0.3 * par_R_character/(par_n_HP - par_n_TS)
                else if (i <= par_n_BS) then  ! �� ���������� = par_R_character * 2
                    r = 1.3 * par_R_character + (i - par_n_HP) * 0.7 * par_R_character/(par_n_BS - par_n_HP)
                else  ! �� ���������� = par_R_character * 1.3_8
                    !r = 2.0 * par_R_character + (i - par_n_BS) * (par_R_END - 2.0 * par_R_character)/(par_n_END - par_n_BS)
                    r = 2.0 * par_R_character + (par_R_END - 2.0 * par_R_character) * (DBLE(i- par_n_BS)/(par_n_END - par_n_BS))**par_kk2
                end if
                
                ! ������ �����
                !gl_x = [gl_x, r * cos(the)]
                !gl_y = [gl_y, r * sin(the) * cos(phi)]
                !gl_z = [gl_z, r * sin(the) * sin(phi)]
                gl_x(node) = r * cos(the)
                gl_y(node) = r * sin(the) * cos(phi)
                gl_z(node) = r * sin(the) * sin(phi)
                node = node + 1
                
                gl_RAY_A(i, j, k) = node - 1
            end do
        end do
    end do
    
    N3 = size(gl_RAY_B(1, 1, :))
    N2 = size(gl_RAY_B(1, :, 1))
    N1 = size(gl_RAY_B(:, 1, 1))
    
    ! ���� ��������� ����� �� ����� B � �� ���������� � ����� ������ ************************************************************
    do k = 1, N3
    	do j = 1, N2
            do i = 1, N1
                
                if (i == 1) then
                    gl_RAY_B(i, j, k) = 1
                    CYCLE
                end if
                
    	        ! ��������� ���������� �������� ���� � ������������
                the = par_pi_8/2.0 + (j) * par_triple_point/(N2)
                phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
                ! ��������� ���������� ����� �� ����
                
                ! �� TS
                if (i <= par_n_TS) then  ! �� ���������� = par_R_character
                    r =  par_R0 + (par_R_character - par_R0) * (REAL(i, KIND = 4)/par_n_TS)**par_kk1
                    !print *, r
                    !pause
                else if (i <= par_n_HP) then  ! �� ���������� = par_R_character * 1.3
                    xx = 1.3 * par_R_character * (1 - cos(the - par_pi_8/2.0))/cos(the - par_pi_8/2.0)
                    r = par_R_character + (i - par_n_TS) * (0.3 * par_R_character + xx) /(par_n_HP - par_n_TS)
                end if
                
                ! ������ �����
                !gl_x = [gl_x, r * cos(the)]
                !gl_y = [gl_y, r * sin(the) * cos(phi)]
                !gl_z = [gl_z, r * sin(the) * sin(phi)]
                
                gl_x(node) = r * cos(the)
                gl_y(node) = r * sin(the) * cos(phi)
                gl_z(node) = r * sin(the) * sin(phi)
                node = node + 1
                
                gl_RAY_B(i, j, k) = node - 1
            end do
        end do
    end do
    
    
    N3 = size(gl_RAY_C(1, 1, :))
    N2 = size(gl_RAY_C(1, :, 1))
    N1 = size(gl_RAY_C(:, 1, 1))
    
    ! ���� ��������� ����� �� ����� C � �� ���������� � ����� ������ ************************************************************
    do k = 1, N3
    	do j = 1, N2
            do i = 1, N1
                
                if(i == 1) then
                    gl_RAY_C(i, j, k) = gl_RAY_B(par_n_HP, j, k)
                    CYCLE
                end if
                
    	        ! ��������� ���������� �������� ���� � ������������
                phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
                ! ��������� ���������� ����� �� ����
                
                x = gl_x(gl_RAY_B(par_n_HP, j, k))
                y = gl_y(gl_RAY_B(par_n_HP, j, k))
                z = gl_z(gl_RAY_B(par_n_HP, j, k))
                rr = (y**2 + z**2)**(0.5)
                
                if (i <= par_n_BS - par_n_HP + 1) then  
                    r = rr + (i - 1) * (2.0 * par_R_character - rr)/(par_n_BS - par_n_HP)
                else
                    !r = 2.0 * par_R_character + (i - (par_n_BS - par_n_HP + 1)) * (par_R_END - 2.0 * par_R_character) /(N1 - (par_n_BS - par_n_HP + 1) )
                    r = 2.0 * par_R_character + (DBLE(i - (par_n_BS - par_n_HP + 1))/(N1 - (par_n_BS - par_n_HP + 1) ))**par_kk2 * (par_R_END - 2.0 * par_R_character)
                end if
                
                ! ������ �����
                !gl_x = [gl_x, x]
                !gl_y = [gl_y, r * cos(phi)]
                !gl_z = [gl_z, r * sin(phi)]
                
                gl_x(node) = x
                gl_y(node) = r * cos(phi)
                gl_z(node) = r * sin(phi)
                node = node + 1
                
                gl_RAY_C(i, j, k) = node - 1
            end do
        end do
    end do
    
    N3 = size(gl_RAY_O(1, 1, :))
    N2 = size(gl_RAY_O(1, :, 1))
    N1 = size(gl_RAY_O(:, 1, 1))
    
    xx = gl_x(gl_RAY_B(par_n_HP, par_m_BC, 1))
    
    ! ���� ��������� ����� �� ����� O � �� ���������� � ����� ������ ************************************************************
    do k = 1, N3
    	do j = 1, N2
            !x = xx - j * (xx - par_R_LEFT)/N2
            x = xx - (DBLE(j)/N2)**par_kk3 * (xx - par_R_LEFT)
            do i = 1, N1
                
    	        ! ��������� ���������� �������� ���� � ������������
                phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
                ! ��������� ���������� ����� �� ����
                
                
                if (i <= par_n_BS - par_n_HP + 1) then  
                    r = 1.3 * par_R_character + (i - 1) * par_R_character * (0.7)/(par_n_BS - par_n_HP)
                else
                    !r = 2.0 * par_R_character + (i - (par_n_BS - par_n_HP + 1)) * (par_R_END - 2.0 * par_R_character) /(N1 - (par_n_BS - par_n_HP + 1) )
                    r = 2.0 * par_R_character + (DBLE(i - (par_n_BS - par_n_HP + 1))/(N1 - (par_n_BS - par_n_HP + 1) ))**par_kk2 * (par_R_END - 2.0 * par_R_character)
                end if
                
                ! ������ �����
                !gl_x = [gl_x, x]
                !gl_y = [gl_y, r * cos(phi)]
                !gl_z = [gl_z, r * sin(phi)]
                
                gl_x(node) = x
                gl_y(node) = r * cos(phi)
                gl_z(node) = r * sin(phi)
                node = node + 1
                
                gl_RAY_O(i, j, k) = node - 1
            end do
        end do
    end do
    
    N3 = size(gl_RAY_K(1, 1, :))
    N2 = size(gl_RAY_K(1, :, 1))
    N1 = size(gl_RAY_K(:, 1, 1))
    
    ! ���� ��������� ����� �� ����� K � �� ���������� � ����� ������ ************************************************************
    do k = 1, N3
    	do j = 1, N2
            do i = 1, N1
                
                if (i == 1) then
                    gl_RAY_K(i, j, k) = 1
                    CYCLE
                end if
                
                if (k /= 1 .and. j == 1) then
                    gl_RAY_K(i, j, k) = gl_RAY_K(i, 1, 1)
                    CYCLE
                end if
                
    	        ! ��������� ���������� �������� ���� � ������������
                the = par_pi_8/2.0 + par_triple_point + (N2 - j + 1) * (par_pi_8/2.0 - par_triple_point)/(N2)
                phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
                ! ��������� ���������� ����� �� ����
                
                r =  par_R0 + (par_R_character - par_R0) * (REAL(i, KIND = 4)/par_n_TS)**par_kk1
                
                
                ! ������ �����
                !gl_x = [gl_x, r * cos(the)]
                !gl_y = [gl_y, r * sin(the) * cos(phi)]
                !gl_z = [gl_z, r * sin(the) * sin(phi)]
                
                gl_x(node) = r * cos(the)
                gl_y(node) = r * sin(the) * cos(phi)
                gl_z(node) = r * sin(the) * sin(phi)
                node = node + 1
                
                gl_RAY_K(i, j, k) = node - 1
            end do
        end do
    end do
    
    N3 = size(gl_RAY_D(1, 1, :))
    N2 = size(gl_RAY_D(1, :, 1))
    N1 = size(gl_RAY_D(:, 1, 1))
    
    ! ���� ��������� ����� �� ����� D � �� ���������� � ����� ������ ************************************************************
    do k = 1, N3
    	do j = 1, N2
            do i = 1, N1
                
                if (i == 1) then
                    if (j < N2) then
                        gl_RAY_D(i, j, k) = gl_RAY_K(par_n_TS, j, k)
                    else
                        gl_RAY_D(i, j, k) = gl_RAY_B(par_n_TS, par_m_BC, k)
                    end if
                    CYCLE
                end if
                
                if (k /= 1 .and. j == 1) then
                    gl_RAY_D(i, j, k) = gl_RAY_D(i, 1, 1)
                    CYCLE
                end if
                
    	        ! ��������� ���������� �������� ���� � ������������
                phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
                ! ��������� ���������� ����� �� ����
                
                if (j < N2) then
                        xx = gl_x(gl_RAY_K(par_n_TS, j, k))
                        y = gl_y(gl_RAY_K(par_n_TS, j, k))
                        z = gl_z(gl_RAY_K(par_n_TS, j, k))
                else
                        xx = gl_x(gl_RAY_B(par_n_TS, par_m_BC, k))
                        y = gl_y(gl_RAY_B(par_n_TS, par_m_BC, k))
                        z = gl_z(gl_RAY_B(par_n_TS, par_m_BC, k))
                end if
                
                r = sqrt(y**2 + z**2)
                !x = xx + (i - 1) * (par_R_LEFT - xx)/(N1 - 1)
                x = xx + (DBLE(i - 1)/(N1 - 1))**par_kk3 * (par_R_LEFT - xx)
                
                
                ! ������ �����
                !gl_x = [gl_x, x]
                !gl_y = [gl_y, r * cos(phi)]
                !gl_z = [gl_z, r * sin(phi)]
                
                gl_x(node) = x
                gl_y(node) = r * cos(phi)
                gl_z(node) = r * sin(phi)
                node = node + 1
                
                gl_RAY_D(i, j, k) = node - 1
            end do
        end do
    end do
    
    N3 = size(gl_RAY_E(1, 1, :))
    N2 = size(gl_RAY_E(1, :, 1))
    N1 = size(gl_RAY_E(:, 1, 1))
    
    ! ���� ��������� ����� �� ����� E � �� ���������� � ����� ������ ************************************************************
    do k = 1, N3
    	do j = 1, N2
            do i = 1, N1
                
                if (i == 1) then
                    gl_RAY_E(i, j, k) = gl_RAY_D(j + 1, size(gl_RAY_D(i + 1, :, k)), k)
                    CYCLE
                end if
                
                if (i == N1) then
                    gl_RAY_E(i, j, k) = gl_RAY_O(1, j, k)
                    CYCLE
                end if
                
                x = gl_x(gl_RAY_E(1, j, k))
                y = gl_y(gl_RAY_E(1, j, k))
                z = gl_z(gl_RAY_E(1, j, k))
                
                x2 = gl_x(gl_RAY_O(1, j, k))
                y2 = gl_y(gl_RAY_O(1, j, k))
                z2 = gl_z(gl_RAY_O(1, j, k))
                
                
                ! ������ �����
                !gl_x = [gl_x, x + (x2 - x) * (i - 1)/(N1 - 1)]
                !gl_y = [gl_y, y + (y2 - y) * (i - 1)/(N1 - 1)]
                !gl_z = [gl_z, z + (z2 - z) * (i - 1)/(N1 - 1)]
                
                gl_x(node) = x + (x2 - x) * (i - 1)/(N1 - 1)
                gl_y(node) = y + (y2 - y) * (i - 1)/(N1 - 1)
                gl_z(node) = z + (z2 - z) * (i - 1)/(N1 - 1)
                node = node + 1
                
                gl_RAY_E(i, j, k) = node - 1
            end do
        end do
    end do
    
    ! ������ ���� ������ � ��������� �� � �������
    
    ! � - ������ ����� ************************************************************************************************************************
    N3 = size(gl_Cell_A(1, 1, :))
    N2 = size(gl_Cell_A(1, :, 1))
    N1 = size(gl_Cell_A(:, 1, 1))
    
    i1 = 1
    
    do k = 1, N3
    	do j = 1, N2
            do i = 1, N1
                
                kk = k + 1
                if (kk > par_l_phi) kk = 1
                
                if (j < par_m_A) then
                    gl_all_Cell(1, i1) = gl_RAY_A(i, j, k)
                    gl_all_Cell(2, i1) = gl_RAY_A(i + 1, j, k)
                    gl_all_Cell(3, i1) = gl_RAY_A(i + 1, j + 1, k)
                    gl_all_Cell(4, i1) = gl_RAY_A(i, j + 1, k)
                    gl_all_Cell(5, i1) = gl_RAY_A(i, j, kk)
                    gl_all_Cell(6, i1) = gl_RAY_A(i + 1, j, kk)
                    gl_all_Cell(7, i1) = gl_RAY_A(i + 1, j + 1, kk)
                    gl_all_Cell(8, i1) = gl_RAY_A(i, j + 1, kk)
                else if (j > par_m_A) then
                    if (i < par_n_HP) then
                        gl_all_Cell(1, i1) = gl_RAY_B(i, j - par_m_A, k)
                        gl_all_Cell(2, i1) = gl_RAY_B(i + 1, j - par_m_A, k)
                        gl_all_Cell(3, i1) = gl_RAY_B(i + 1, j - par_m_A + 1, k)
                        gl_all_Cell(4, i1) = gl_RAY_B(i, j - par_m_A + 1, k)
                        gl_all_Cell(5, i1) = gl_RAY_B(i, j - par_m_A, kk)
                        gl_all_Cell(6, i1) = gl_RAY_B(i + 1, j - par_m_A, kk)
                        gl_all_Cell(7, i1) = gl_RAY_B(i + 1, j - par_m_A + 1, kk)
                        gl_all_Cell(8, i1) = gl_RAY_B(i, j - par_m_A + 1, kk)
                    else if (i == par_n_HP) then
                        gl_all_Cell(1, i1) = gl_RAY_B(i, j - par_m_A, k)
                        gl_all_Cell(2, i1) = gl_RAY_C(2, j - par_m_A, k)
                        gl_all_Cell(3, i1) = gl_RAY_C(2, j - par_m_A + 1, k)
                        gl_all_Cell(4, i1) = gl_RAY_B(i, j - par_m_A + 1, k)
                        gl_all_Cell(5, i1) = gl_RAY_B(i, j - par_m_A, kk)
                        gl_all_Cell(6, i1) = gl_RAY_C(2, j - par_m_A, kk)
                        gl_all_Cell(7, i1) = gl_RAY_C(2, j - par_m_A + 1, kk)
                        gl_all_Cell(8, i1) = gl_RAY_B(i, j - par_m_A + 1, kk)
                    else ! i > par_n_HP
                        gl_all_Cell(1, i1) = gl_RAY_C(i - par_n_HP + 1, j - par_m_A, k)
                        gl_all_Cell(2, i1) = gl_RAY_C(i - par_n_HP + 2, j - par_m_A, k)
                        gl_all_Cell(3, i1) = gl_RAY_C(i - par_n_HP + 2, j - par_m_A + 1, k)
                        gl_all_Cell(4, i1) = gl_RAY_C(i - par_n_HP + 1, j - par_m_A + 1, k)
                        gl_all_Cell(5, i1) = gl_RAY_C(i - par_n_HP + 1, j - par_m_A, kk)
                        gl_all_Cell(6, i1) = gl_RAY_C(i - par_n_HP + 2, j - par_m_A, kk)
                        gl_all_Cell(7, i1) = gl_RAY_C(i - par_n_HP + 2, j - par_m_A + 1, kk)
                        gl_all_Cell(8, i1) = gl_RAY_C(i - par_n_HP + 1, j - par_m_A + 1, kk)
                    end if
                else ! j == par_m_A
                    if (i < par_n_HP) then
                        gl_all_Cell(1, i1) = gl_RAY_A(i, j, k)
                        gl_all_Cell(2, i1) = gl_RAY_A(i + 1, j, k)
                        gl_all_Cell(3, i1) = gl_RAY_B(i + 1, 1, k)
                        gl_all_Cell(4, i1) = gl_RAY_B(i, 1, k)
                        gl_all_Cell(5, i1) = gl_RAY_A(i, j, kk)
                        gl_all_Cell(6, i1) = gl_RAY_A(i + 1, j, kk)
                        gl_all_Cell(7, i1) = gl_RAY_B(i + 1, 1, kk)
                        gl_all_Cell(8, i1) = gl_RAY_B(i, 1, kk)
                    else if (i == par_n_HP) then
                        gl_all_Cell(1, i1) = gl_RAY_A(i, j, k)
                        gl_all_Cell(2, i1) = gl_RAY_A(i + 1, j, k)
                        gl_all_Cell(3, i1) = gl_RAY_C(2, 1, k)
                        gl_all_Cell(4, i1) = gl_RAY_B(i, 1, k)
                        gl_all_Cell(5, i1) = gl_RAY_A(i, j, kk)
                        gl_all_Cell(6, i1) = gl_RAY_A(i + 1, j, kk)
                        gl_all_Cell(7, i1) = gl_RAY_C(2, 1, kk)
                        gl_all_Cell(8, i1) = gl_RAY_B(i, 1, kk)
                    else ! i > par_n_HP
                        gl_all_Cell(1, i1) = gl_RAY_A(i, j, k)
                        gl_all_Cell(2, i1) = gl_RAY_A(i + 1, j, k)
                        gl_all_Cell(3, i1) = gl_RAY_C(i - par_n_HP + 2, 1, k)
                        gl_all_Cell(4, i1) = gl_RAY_C(i - par_n_HP + 1, 1, k)
                        gl_all_Cell(5, i1) = gl_RAY_A(i, j, kk)
                        gl_all_Cell(6, i1) = gl_RAY_A(i + 1, j, kk)
                        gl_all_Cell(7, i1) = gl_RAY_C(i - par_n_HP + 2, 1, kk)
                        gl_all_Cell(8, i1) = gl_RAY_C(i - par_n_HP + 1, 1, kk)
                    end if
                end if
                
                gl_Cell_A(i, j, k) = i1
                i1 = i1 + 1
                
            end do
        end do
    end do
    
    ! B - ������ ����� ************************************************************************************************************************
    N3 = size(gl_Cell_B(1, 1, :))
    N2 = size(gl_Cell_B(1, :, 1))
    N1 = size(gl_Cell_B(:, 1, 1))
    
    do k = 1, N3
    	do j = 1, N2
            do i = 1, N1
                
                kk = k + 1
                if (kk > par_l_phi) kk = 1
                
                if (j < par_m_K) then
                    if (i < par_n_TS) then
                        gl_all_Cell(1, i1) = gl_RAY_K(i, j, k)
                        gl_all_Cell(2, i1) = gl_RAY_K(i, j + 1, k)
                        gl_all_Cell(3, i1) = gl_RAY_K(i + 1, j + 1, k)
                        gl_all_Cell(4, i1) = gl_RAY_K(i + 1, j, k)
                        gl_all_Cell(5, i1) = gl_RAY_K(i, j, kk)
                        gl_all_Cell(6, i1) = gl_RAY_K(i, j + 1, kk)
                        gl_all_Cell(7, i1) = gl_RAY_K(i + 1, j + 1, kk)
                        gl_all_Cell(8, i1) = gl_RAY_K(i + 1, j, kk)
                    else if (i == par_n_TS) then
                        gl_all_Cell(1, i1) = gl_RAY_K(i, j, k)
                        gl_all_Cell(2, i1) = gl_RAY_K(i, j + 1, k)
                        gl_all_Cell(3, i1) = gl_RAY_D(2, j + 1, k)
                        gl_all_Cell(4, i1) = gl_RAY_D(2, j, k)
                        gl_all_Cell(5, i1) = gl_RAY_K(i, j, kk)
                        gl_all_Cell(6, i1) = gl_RAY_K(i, j + 1, kk)
                        gl_all_Cell(7, i1) = gl_RAY_D(2, j + 1, kk)
                        gl_all_Cell(8, i1) = gl_RAY_D(2, j, kk)
                    else  ! i > par_n_TS
                        gl_all_Cell(1, i1) = gl_RAY_D(i - par_n_TS + 1, j, k)
                        gl_all_Cell(2, i1) = gl_RAY_D(i - par_n_TS + 1, j + 1, k)
                        gl_all_Cell(3, i1) = gl_RAY_D(i - par_n_TS + 1 + 1, j + 1, k)
                        gl_all_Cell(4, i1) = gl_RAY_D(i - par_n_TS + 1 + 1, j, k)
                        gl_all_Cell(5, i1) = gl_RAY_D(i - par_n_TS + 1, j, kk)
                        gl_all_Cell(6, i1) = gl_RAY_D(i - par_n_TS + 1, j + 1, kk)
                        gl_all_Cell(7, i1) = gl_RAY_D(i - par_n_TS + 1 + 1, j + 1, kk)
                        gl_all_Cell(8, i1) = gl_RAY_D(i - par_n_TS + 1 + 1, j, kk)
                    end if
                else  ! j == par_m_K
                    if (i < par_n_TS) then
                        gl_all_Cell(1, i1) = gl_RAY_K(i, j, k)
                        gl_all_Cell(2, i1) = gl_RAY_B(i, size(gl_RAY_B(i, :, k)), k)
                        gl_all_Cell(3, i1) = gl_RAY_B(i + 1, size(gl_RAY_B(i, :, k)), k)
                        gl_all_Cell(4, i1) = gl_RAY_K(i + 1, j, k)
                        gl_all_Cell(5, i1) = gl_RAY_K(i, j, kk)
                        gl_all_Cell(6, i1) = gl_RAY_B(i, size(gl_RAY_B(i, :, k)), kk)
                        gl_all_Cell(7, i1) = gl_RAY_B(i + 1, size(gl_RAY_B(i, :, k)), kk)
                        gl_all_Cell(8, i1) = gl_RAY_K(i + 1, j, kk)
                    else if (i == par_n_TS) then
                        gl_all_Cell(1, i1) = gl_RAY_K(i, j, k)
                        gl_all_Cell(2, i1) = gl_RAY_B(i, size(gl_RAY_B(i, :, k)), k)
                        gl_all_Cell(3, i1) = gl_RAY_D(2, size(gl_RAY_D(2, :, k)), k)
                        gl_all_Cell(4, i1) = gl_RAY_D(2, size(gl_RAY_D(2, :, k)) - 1, k)
                        gl_all_Cell(5, i1) = gl_RAY_K(i, j, kk)
                        gl_all_Cell(6, i1) = gl_RAY_B(i, size(gl_RAY_B(i, :, k)), kk)
                        gl_all_Cell(7, i1) = gl_RAY_D(2, size(gl_RAY_D(2, :, k)), kk)
                        gl_all_Cell(8, i1) = gl_RAY_D(2, size(gl_RAY_D(2, :, k)) - 1, kk)
                    else   ! i > par_n_TS
                        gl_all_Cell(1, i1) = gl_RAY_D(i - par_n_TS + 1, j, k)
                        gl_all_Cell(2, i1) = gl_RAY_D(i - par_n_TS + 1, j + 1, k)
                        gl_all_Cell(3, i1) = gl_RAY_D(i - par_n_TS + 1 + 1, j + 1, k)
                        gl_all_Cell(4, i1) = gl_RAY_D(i - par_n_TS + 1 + 1, j, k)
                        gl_all_Cell(5, i1) = gl_RAY_D(i - par_n_TS + 1, j, kk)
                        gl_all_Cell(6, i1) = gl_RAY_D(i - par_n_TS + 1, j + 1, kk)
                        gl_all_Cell(7, i1) = gl_RAY_D(i - par_n_TS + 1 + 1, j + 1, kk)
                        gl_all_Cell(8, i1) = gl_RAY_D(i - par_n_TS + 1 + 1, j, kk)
                    end if
                end if
                
                
                gl_Cell_B(i, j, k) = i1
                i1 = i1 + 1
                
            end do
        end do
    end do
    
    ! C - ������ ����� ************************************************************************************************************************
    N3 = size(gl_Cell_C(1, 1, :))
    N2 = size(gl_Cell_C(1, :, 1))
    N1 = size(gl_Cell_C(:, 1, 1))
    
    do k = 1, N3
    	do j = 1, N2
            do i = 1, N1
                
                kk = k + 1
                if (kk > par_l_phi) kk = 1
                
                if (j > 1) then
                    if (i < par_n_HP - par_n_TS + 1) then
                        gl_all_Cell(1, i1) = gl_RAY_E(i, j - 1, k)
                        gl_all_Cell(2, i1) = gl_RAY_E(i + 1, j - 1, k)
                        gl_all_Cell(3, i1) = gl_RAY_E(i + 1, j - 1 + 1, k)
                        gl_all_Cell(4, i1) = gl_RAY_E(i, j - 1 + 1, k)
                        gl_all_Cell(5, i1) = gl_RAY_E(i, j - 1, kk)
                        gl_all_Cell(6, i1) = gl_RAY_E(i + 1, j - 1, kk)
                        gl_all_Cell(7, i1) = gl_RAY_E(i + 1, j - 1 + 1, kk)
                        gl_all_Cell(8, i1) = gl_RAY_E(i, j - 1 + 1, kk)
                    else if (i == par_n_HP - par_n_TS + 1) then
                        gl_all_Cell(1, i1) = gl_RAY_E(i, j - 1, k)
                        gl_all_Cell(2, i1) = gl_RAY_O(2, j - 1, k)
                        gl_all_Cell(3, i1) = gl_RAY_O(2, j - 1 + 1, k)
                        gl_all_Cell(4, i1) = gl_RAY_E(i, j - 1 + 1, k)
                        gl_all_Cell(5, i1) = gl_RAY_E(i, j - 1, kk)
                        gl_all_Cell(6, i1) = gl_RAY_O(2, j - 1, kk)
                        gl_all_Cell(7, i1) = gl_RAY_O(2, j - 1 + 1, kk)
                        gl_all_Cell(8, i1) = gl_RAY_E(i, j - 1 + 1, kk)
                    else  ! i > par_n_HP - par_n_TS + 1
                        gl_all_Cell(1, i1) = gl_RAY_O(i - (par_n_HP - par_n_TS + 1) + 1, j - 1, k)
                        gl_all_Cell(2, i1) = gl_RAY_O(i - (par_n_HP - par_n_TS + 1) + 2, j - 1, k)
                        gl_all_Cell(3, i1) = gl_RAY_O(i - (par_n_HP - par_n_TS + 1) + 2, j - 1 + 1, k)
                        gl_all_Cell(4, i1) = gl_RAY_O(i - (par_n_HP - par_n_TS + 1) + 1, j - 1 + 1, k)
                        gl_all_Cell(5, i1) = gl_RAY_O(i - (par_n_HP - par_n_TS + 1) + 1, j - 1, kk)
                        gl_all_Cell(6, i1) = gl_RAY_O(i - (par_n_HP - par_n_TS + 1) + 2, j - 1, kk)
                        gl_all_Cell(7, i1) = gl_RAY_O(i - (par_n_HP - par_n_TS + 1) + 2, j - 1 + 1, kk)
                        gl_all_Cell(8, i1) = gl_RAY_O(i - (par_n_HP - par_n_TS + 1) + 1, j - 1 + 1, kk)
                    end if
                else  ! j == 1
                    if (i < par_n_HP - par_n_TS + 1) then
                        gl_all_Cell(1, i1) = gl_RAY_B(par_n_TS + i - 1, size(gl_RAY_B(par_n_TS + i - 1, :, k)), k)
                        gl_all_Cell(2, i1) = gl_RAY_B(par_n_TS + i, size(gl_RAY_B(par_n_TS + i - 1, :, k)), k)
                        gl_all_Cell(3, i1) = gl_RAY_E(i + 1, j, k)
                        gl_all_Cell(4, i1) = gl_RAY_E(i, j, k)
                        gl_all_Cell(5, i1) = gl_RAY_B(par_n_TS + i - 1, size(gl_RAY_B(par_n_TS + i - 1, :, k)), kk)
                        gl_all_Cell(6, i1) = gl_RAY_B(par_n_TS + i, size(gl_RAY_B(par_n_TS + i - 1, :, k)), kk)
                        gl_all_Cell(7, i1) = gl_RAY_E(i + 1, j, kk)
                        gl_all_Cell(8, i1) = gl_RAY_E(i, j, kk)
                    else if (i == par_n_HP - par_n_TS + 1) then
                        gl_all_Cell(1, i1) = gl_RAY_B(par_n_TS + i - 1, size(gl_RAY_B(par_n_TS + i - 1, :, k)), k)
                        gl_all_Cell(2, i1) = gl_RAY_C(2, size(gl_RAY_C(2, :, k)), k)
                        gl_all_Cell(3, i1) = gl_RAY_O(2, j, k)
                        gl_all_Cell(4, i1) = gl_RAY_E(i, j, k)
                        gl_all_Cell(5, i1) = gl_RAY_B(par_n_TS + i - 1, size(gl_RAY_B(par_n_TS + i - 1, :, k)), kk)
                        gl_all_Cell(6, i1) = gl_RAY_C(2, size(gl_RAY_C(2, :, k)), kk)
                        gl_all_Cell(7, i1) = gl_RAY_O(2, j, kk)
                        gl_all_Cell(8, i1) = gl_RAY_E(i, j, kk)
                    else  ! i > par_n_HP - par_n_TS + 1
                        gl_all_Cell(1, i1) = gl_RAY_C(i - (par_n_HP - par_n_TS + 1) + 1, size(gl_RAY_C(1, :, k)), k)
                        gl_all_Cell(2, i1) = gl_RAY_C(i - (par_n_HP - par_n_TS + 1) + 2, size(gl_RAY_C(1, :, k)), k)
                        gl_all_Cell(3, i1) = gl_RAY_O(i - (par_n_HP - par_n_TS + 1) + 2, j, k)
                        gl_all_Cell(4, i1) = gl_RAY_O(i - (par_n_HP - par_n_TS + 1) + 1, j, k)
                        gl_all_Cell(5, i1) = gl_RAY_C(i - (par_n_HP - par_n_TS + 1) + 1, size(gl_RAY_C(1, :, k)), kk)
                        gl_all_Cell(6, i1) = gl_RAY_C(i - (par_n_HP - par_n_TS + 1) + 2, size(gl_RAY_C(1, :, k)), kk)
                        gl_all_Cell(7, i1) = gl_RAY_O(i - (par_n_HP - par_n_TS + 1) + 2, j, kk)
                        gl_all_Cell(8, i1) = gl_RAY_O(i - (par_n_HP - par_n_TS + 1) + 1, j, kk)
                    end if
                end if
                
                
                gl_Cell_C(i, j, k) = i1
                i1 = i1 + 1
                
            end do
        end do
    end do
    
    ! ������ ������ � �� �������� (������� ������� �������)
    
    ! ������ � ����� ������ �
    N3 = size(gl_Cell_A(1, 1, :))
    N2 = size(gl_Cell_A(1, :, 1))
    N1 = size(gl_Cell_A(:, 1, 1))
    
    do k = 1, N3
    	do j = 1, N2
            do i = 1, N1
                
                kk = k + 1                      ! ���� ����� ����
                if (kk > par_l_phi) kk = 1
                
                kk2 = k - 1                     ! ���� ����� ����
                if (kk2 < 1) kk2 = par_l_phi
                
                
                ! ��������� ������� ������ (�� ���� �����)
                if (i == N1) then
                    if (j < par_m_A) then
                        gl_Cell_neighbour(1, gl_Cell_A(i, j, k)) = -1   ! ������� (���������� �����)
                    else
                        gl_Cell_neighbour(1, gl_Cell_A(i, j, k)) = -3   ! ������� ������� �������
                    end if  
                else
                    gl_Cell_neighbour(1, gl_Cell_A(i, j, k)) = gl_Cell_A(i + 1, j, k)
                end if
                
                ! ��������� ������� ������ (�� ���� ����)
                if (i == 1) then
                    continue
                else
                    gl_Cell_neighbour(2, gl_Cell_A(i, j, k)) = gl_Cell_A(i - 1, j, k)
                end if
                
                ! ��������� �������� ������ (�� ���� � ���������, � ������� ����������)
                if (j == N2) then
                    if (i < par_n_TS) then
                        gl_Cell_neighbour(3, gl_Cell_A(i, j, k)) = gl_Cell_B(i, par_m_K, k)
                    else
                        gl_Cell_neighbour(3, gl_Cell_A(i, j, k)) = gl_Cell_C(i - par_n_TS + 1, 1, k)
                    end if
                else
                    gl_Cell_neighbour(3, gl_Cell_A(i, j, k)) = gl_Cell_A(i, j + 1, k)
                end if
                
                ! ��������� ��������� ������ (�� ���� � ���������, � ������� ����������)
                if (j == 1) then
                    continue
                else
                    gl_Cell_neighbour(4, gl_Cell_A(i, j, k)) = gl_Cell_A(i, j - 1, k)
                end if
                
                ! ��������� ������ � ������� ������ (����� � ����)
                gl_Cell_neighbour(5, gl_Cell_A(i, j, k)) = gl_Cell_A(i, j, kk)
                gl_Cell_neighbour(6, gl_Cell_A(i, j, k)) = gl_Cell_A(i, j, kk2)
            end do
        end do
    end do
    
    ! ������ ������ B
    N3 = size(gl_Cell_B(1, 1, :))
    N2 = size(gl_Cell_B(1, :, 1))
    N1 = size(gl_Cell_B(:, 1, 1))
    
    do k = 1, N3
    	do j = 1, N2
            do i = 1, N1
                
                kk = k + 1                      ! ���� ����� ����
                if (kk > par_l_phi) kk = 1
                
                kk2 = k - 1                     ! ���� ����� ����
                if (kk2 < 1) kk2 = par_l_phi
                
                
                ! ��������� ������� ������ (�� ���� �����)
                if (i == N1) then
                    gl_Cell_neighbour(1, gl_Cell_B(i, j, k)) = -2   ! �������� �������
                else
                    gl_Cell_neighbour(1, gl_Cell_B(i, j, k)) = gl_Cell_B(i + 1, j, k)
                end if
                
                ! ��������� ������� ������ (�� ���� ����)
                if (i == 1) then
                    continue
                else
                    gl_Cell_neighbour(2, gl_Cell_B(i, j, k)) = gl_Cell_B(i - 1, j, k)
                end if
                
                ! ��������� �������� ������ (�� ���� � ��������� ������ �� �����)
                if (j == N2) then
                    if(i < par_n_TS) then
                        gl_Cell_neighbour(3, gl_Cell_B(i, j, k)) = gl_Cell_A(i, size(gl_Cell_A(i,:, k)), k)
                    else
                        gl_Cell_neighbour(3, gl_Cell_B(i, j, k)) = gl_Cell_C(1, i - par_n_TS + 1, k)
                    end if
                else
                    gl_Cell_neighbour(3, gl_Cell_B(i, j, k)) = gl_Cell_B(i, j + 1, k)
                end if
                
                ! ��������� ��������� ������ (�� ���� � ���������)
                if (j == 1) then
                    continue
                else
                    gl_Cell_neighbour(4, gl_Cell_B(i, j, k)) = gl_Cell_B(i, j - 1, k)
                end if
                
                ! ��������� ������ � ������� ������ (����� � ����)
                gl_Cell_neighbour(5, gl_Cell_B(i, j, k)) = gl_Cell_B(i, j, kk)
                gl_Cell_neighbour(6, gl_Cell_B(i, j, k)) = gl_Cell_B(i, j, kk2)
            end do
        end do
    end do
    
    ! ������ ������ C
    N3 = size(gl_Cell_C(1, 1, :))
    N2 = size(gl_Cell_C(1, :, 1))
    N1 = size(gl_Cell_C(:, 1, 1))
    
    do k = 1, N3
    	do j = 1, N2
            do i = 1, N1
                
                kk = k + 1                      ! ���� ����� ����
                if (kk > par_l_phi) kk = 1
                
                kk2 = k - 1                     ! ���� ����� ����
                if (kk2 < 1) kk2 = par_l_phi
                
                
                ! ��������� ������� ������ (�� ���� �����)
                if (i == N1) then
                    gl_Cell_neighbour(1, gl_Cell_C(i, j, k)) = -3   ! ������� �������
                else
                    gl_Cell_neighbour(1, gl_Cell_C(i, j, k)) = gl_Cell_C(i + 1, j, k)
                end if
                
                ! ��������� ������� ������ (�� ���� ����)
                if (i == 1) then
                    gl_Cell_neighbour(2, gl_Cell_C(i, j, k)) = gl_Cell_B(j + par_n_TS - 1, size(gl_Cell_B(j + par_n_TS - 1, :, k)), k)
                else
                    gl_Cell_neighbour(2, gl_Cell_C(i, j, k)) = gl_Cell_C(i - 1, j, k)
                end if
                
                ! ��������� �������� ������ (�� ���� � ���������, � ������� ���������� ����!)
                if (j == N2) then
                    gl_Cell_neighbour(3, gl_Cell_C(i, j, k)) = -2  ! �������� �������
                else
                    gl_Cell_neighbour(3, gl_Cell_C(i, j, k)) = gl_Cell_C(i, j + 1, k)
                end if
                
                ! ��������� ��������� ������ (�� ���� � ���������, � ������� ����������)
                if (j == 1) then
                    gl_Cell_neighbour(4, gl_Cell_C(i, j, k)) = gl_Cell_A(i + par_n_TS - 1, size(gl_Cell_A(i + par_n_TS - 1, : , k)), k)
                else
                    gl_Cell_neighbour(4, gl_Cell_C(i, j, k)) = gl_Cell_C(i, j - 1, k)
                end if
                
                ! ��������� ������ � ������� ������ (����� � ����)
                gl_Cell_neighbour(5, gl_Cell_C(i, j, k)) = gl_Cell_C(i, j, kk)
                gl_Cell_neighbour(6, gl_Cell_C(i, j, k)) = gl_Cell_C(i, j, kk2)
            end do
        end do
    end do
    
    if (par_developer_info) print*, "Build_Mesh_start: Number all points: ", size(gl_x), " == ", node - 1
    
    ! �������� �����
    
    node = 1
    
    ! ������ � ����� ������ �
    N3 = size(gl_Cell_A(1, 1, :))
    N2 = size(gl_Cell_A(1, :, 1))
    N1 = size(gl_Cell_A(:, 1, 1))
    
    do k = 1, N3
    	do j = 1, N2
            do i = 1, N1
                
                kk = k + 1                      ! ���� ����� ����
                if (kk > par_l_phi) kk = 1
                
                ni = gl_Cell_A(i, j, k)
                
                ! ��������� 1-�� ����� ��� ���� ������
                if (i < N1) then
                    gl_all_Gran(1, node) = gl_all_Cell(2, ni)
                    gl_all_Gran(2, node) = gl_all_Cell(3, ni)
                    gl_all_Gran(3, node) = gl_all_Cell(7, ni)
                    gl_all_Gran(4, node) = gl_all_Cell(6, ni)
                    
                    gl_Gran_neighbour(1, node) = ni
                    gl_Gran_neighbour(2, node) = gl_Cell_A(i + 1, j, k)
                    
                    gl_Cell_gran(1, ni) = node
                    gl_Cell_gran(2, gl_Cell_A(i + 1, j, k)) = node
                    
                    node = node + 1
                else
                    gl_all_Gran(1, node) = gl_all_Cell(2, ni)
                    gl_all_Gran(2, node) = gl_all_Cell(3, ni)
                    gl_all_Gran(3, node) = gl_all_Cell(7, ni)
                    gl_all_Gran(4, node) = gl_all_Cell(6, ni)
                    
                    gl_Gran_neighbour(1, node) = ni
                    
                    if (j < par_m_A) then
                        gl_Gran_neighbour(2, node) = -1
                    else
                        gl_Gran_neighbour(2, node) = -3
                    end if
                    
                    gl_Cell_gran(1, ni) = node
                    node = node + 1
                end if
                
                ! ��������� 4-�� ����� ��� ���� ������
                if (j /= 1) then
                    gl_all_Gran(1, node) = gl_all_Cell(1, ni)
                    gl_all_Gran(2, node) = gl_all_Cell(2, ni)
                    gl_all_Gran(3, node) = gl_all_Cell(6, ni)
                    gl_all_Gran(4, node) = gl_all_Cell(5, ni)
                    
                    gl_Gran_neighbour(1, node) = ni
                    gl_Gran_neighbour(2, node) = gl_Cell_A(i, j - 1, k)
                    
                    gl_Cell_gran(4, ni) = node
                    gl_Cell_gran(3, gl_Cell_A(i, j - 1, k)) = node
                    
                    node = node + 1
                end if
                
                ! ��������� 5-�� ����� ��� ���� ������
                gl_all_Gran(1, node) = gl_all_Cell(5, ni)
                gl_all_Gran(2, node) = gl_all_Cell(6, ni)
                gl_all_Gran(3, node) = gl_all_Cell(7, ni)
                gl_all_Gran(4, node) = gl_all_Cell(8, ni)
                    
                gl_Gran_neighbour(1, node) = ni
                gl_Gran_neighbour(2, node) = gl_Cell_A(i, j, kk)
                    
                gl_Cell_gran(5, ni) = node
                gl_Cell_gran(6, gl_Cell_A(i, j, kk)) = node
                    
                node = node + 1
            end do
        end do
    end do
    
     ! ������ ������ B
    N3 = size(gl_Cell_B(1, 1, :))
    N2 = size(gl_Cell_B(1, :, 1))
    N1 = size(gl_Cell_B(:, 1, 1))
    
    do k = 1, N3
    	do j = 1, N2
            do i = 1, N1
                
                kk = k + 1                      ! ���� ����� ����
                if (kk > par_l_phi) kk = 1
                
                ni = gl_Cell_B(i, j, k)
                
                ! ��������� 1-�� ����� ��� ���� ������
                if (i < N1) then
                    gl_all_Gran(1, node) = gl_all_Cell(3, ni)
                    gl_all_Gran(2, node) = gl_all_Cell(4, ni)
                    gl_all_Gran(3, node) = gl_all_Cell(8, ni)
                    gl_all_Gran(4, node) = gl_all_Cell(7, ni)
                    
                    gl_Gran_neighbour(1, node) = ni
                    gl_Gran_neighbour(2, node) = gl_Cell_B(i + 1, j, k)
                    
                    gl_Cell_gran(1, ni) = node
                    gl_Cell_gran(2, gl_Cell_B(i + 1, j, k)) = node
                    
                    node = node + 1
                else
                    gl_all_Gran(1, node) = gl_all_Cell(3, ni)
                    gl_all_Gran(2, node) = gl_all_Cell(4, ni)
                    gl_all_Gran(3, node) = gl_all_Cell(8, ni)
                    gl_all_Gran(4, node) = gl_all_Cell(7, ni)
                    
                    gl_Gran_neighbour(1, node) = ni
                    gl_Gran_neighbour(2, node) = -2
                    
                    gl_Cell_gran(1, ni) = node
                    
                    node = node + 1
                end if
                
                ! ��������� 4-�� ����� ��� ���� ������
                if (j /= 1) then
                    gl_all_Gran(1, node) = gl_all_Cell(4, ni)
                    gl_all_Gran(2, node) = gl_all_Cell(1, ni)
                    gl_all_Gran(3, node) = gl_all_Cell(5, ni)
                    gl_all_Gran(4, node) = gl_all_Cell(8, ni)
                    
                    gl_Gran_neighbour(1, node) = ni
                    gl_Gran_neighbour(2, node) = gl_Cell_B(i, j - 1, k)
                    
                    gl_Cell_gran(4, ni) = node
                    gl_Cell_gran(3, gl_Cell_B(i, j - 1, k)) = node
                    
                    node = node + 1
                end if
                
                ! ��������� 5-�� ����� ��� ���� ������
                gl_all_Gran(1, node) = gl_all_Cell(5, ni)
                gl_all_Gran(2, node) = gl_all_Cell(6, ni)
                gl_all_Gran(3, node) = gl_all_Cell(7, ni)
                gl_all_Gran(4, node) = gl_all_Cell(8, ni)
                    
                gl_Gran_neighbour(1, node) = ni
                gl_Gran_neighbour(2, node) = gl_Cell_B(i, j, kk)
                    
                gl_Cell_gran(5, ni) = node
                gl_Cell_gran(6, gl_Cell_B(i, j, kk)) = node
                    
                node = node + 1
                
                if (j == N2) then  ! 3-�� ����� ��� �������
                    gl_all_Gran(1, node) = gl_all_Cell(2, ni)
                    gl_all_Gran(2, node) = gl_all_Cell(3, ni)
                    gl_all_Gran(3, node) = gl_all_Cell(7, ni)
                    gl_all_Gran(4, node) = gl_all_Cell(6, ni)
                    
                    gl_Gran_neighbour(1, node) = ni
                    
                    if(i < par_n_TS) then
                        gl_Gran_neighbour(2, node) = gl_Cell_A(i, size(gl_Cell_A(i,:, k)), k)
                        gl_Cell_gran(3, gl_Cell_A(i, size(gl_Cell_A(i,:, k)), k)) = node
                    else
                        gl_Gran_neighbour(2, node) = gl_Cell_C(1, i - par_n_TS + 1, k)
                        gl_Cell_gran(2, gl_Cell_C(1, i - par_n_TS + 1, k)) = node
                    end if
                    
                    gl_Cell_gran(3, ni) = node
                    
                    node = node + 1
                end if
                
                
            end do
        end do
    end do
    
    ! ������ ������ C
    N3 = size(gl_Cell_C(1, 1, :))
    N2 = size(gl_Cell_C(1, :, 1))
    N1 = size(gl_Cell_C(:, 1, 1))
    
    do k = 1, N3
    	do j = 1, N2
            do i = 1, N1
                
                kk = k + 1                      ! ���� ����� ����
                if (kk > par_l_phi) kk = 1
                
                ni = gl_Cell_C(i, j, k)
                
                ! ��������� 1-�� ����� ��� ���� ������
                if (i < N1) then
                    gl_all_Gran(1, node) = gl_all_Cell(2, ni)
                    gl_all_Gran(2, node) = gl_all_Cell(3, ni)
                    gl_all_Gran(3, node) = gl_all_Cell(7, ni)
                    gl_all_Gran(4, node) = gl_all_Cell(6, ni)
                    
                    gl_Gran_neighbour(1, node) = ni
                    gl_Gran_neighbour(2, node) = gl_Cell_C(i + 1, j, k)
                    
                    gl_Cell_gran(1, ni) = node
                    gl_Cell_gran(2, gl_Cell_C(i + 1, j, k)) = node
                    
                    node = node + 1
                else
                    gl_all_Gran(1, node) = gl_all_Cell(2, ni)
                    gl_all_Gran(2, node) = gl_all_Cell(3, ni)
                    gl_all_Gran(3, node) = gl_all_Cell(7, ni)
                    gl_all_Gran(4, node) = gl_all_Cell(6, ni)
                    
                    gl_Gran_neighbour(1, node) = ni
                    gl_Gran_neighbour(2, node) = -3    ! ������� �������
                    
                    gl_Cell_gran(1, ni) = node
                    
                    node = node + 1
                end if
                
                ! ��������� 4-�� ����� ��� ���� ������
                if (j /= 1) then
                    gl_all_Gran(1, node) = gl_all_Cell(1, ni)
                    gl_all_Gran(2, node) = gl_all_Cell(2, ni)
                    gl_all_Gran(3, node) = gl_all_Cell(6, ni)
                    gl_all_Gran(4, node) = gl_all_Cell(5, ni)
                    
                    gl_Gran_neighbour(1, node) = ni
                    gl_Gran_neighbour(2, node) = gl_Cell_C(i, j - 1, k)
                    
                    gl_Cell_gran(4, ni) = node
                    gl_Cell_gran(3, gl_Cell_C(i, j - 1, k)) = node
                    
                    node = node + 1
                end if
                
                if (j == N2) then  ! ������ ����� �� �������
                    gl_all_Gran(1, node) = gl_all_Cell(3, ni)
                    gl_all_Gran(2, node) = gl_all_Cell(4, ni)
                    gl_all_Gran(3, node) = gl_all_Cell(8, ni)
                    gl_all_Gran(4, node) = gl_all_Cell(7, ni)
                    
                    gl_Gran_neighbour(1, node) = ni
                    gl_Gran_neighbour(2, node) = -2
                    
                    gl_Cell_gran(3, ni) = node
                    node = node + 1
                end if
                
                
                ! ��������� 5-�� ����� ��� ���� ������
                gl_all_Gran(1, node) = gl_all_Cell(5, ni)
                gl_all_Gran(2, node) = gl_all_Cell(6, ni)
                gl_all_Gran(3, node) = gl_all_Cell(7, ni)
                gl_all_Gran(4, node) = gl_all_Cell(8, ni)
                    
                gl_Gran_neighbour(1, node) = ni
                gl_Gran_neighbour(2, node) = gl_Cell_C(i, j, kk)
                    
                gl_Cell_gran(5, ni) = node
                gl_Cell_gran(6, gl_Cell_C(i, j, kk)) = node
                    
                node = node + 1
                
                if (j == 1) then  ! 4-�� ����� ��� �������
                    gl_all_Gran(1, node) = gl_all_Cell(1, ni)
                    gl_all_Gran(2, node) = gl_all_Cell(2, ni)
                    gl_all_Gran(3, node) = gl_all_Cell(6, ni)
                    gl_all_Gran(4, node) = gl_all_Cell(5, ni)
                    
                    gl_Gran_neighbour(1, node) = ni
                    
                    gl_Gran_neighbour(2, node) = gl_Cell_A(i + par_n_TS - 1, size(gl_Cell_A(i + par_n_TS - 1, : , k)), k)
                    gl_Cell_gran(3, gl_Cell_A(i + par_n_TS - 1, size(gl_Cell_A(i + par_n_TS - 1, : , k)), k)) = node
                    
                    gl_Cell_gran(4, ni) = node
                    
                    node = node + 1
                end if
                
            end do
        end do
    end do
    
    
    if (par_developer_info) print *, "Build_Mesh_start: All Cell = ", i1 - 1
    
    if (par_developer_info) print*, "Build_Mesh_start: Number all gran: ", size(gl_all_Gran(1,:)), " == ", node - 1
    
    !print*, par_n_END * par_l_phi * (par_m_A + par_m_BC) + par_m_K * (par_n_TS + par_m_O) * par_l_phi + par_l_phi * (par_n_END - par_n_TS + 1) * par_m_O - &
    !        (par_m_A + par_m_BC + par_m_K - 1) * par_l_phi - par_n_END * (par_l_phi - 1) - (par_n_TS + par_m_O - 1) * (par_l_phi - 1)
    
    end subroutine Build_Mesh_start
    
    subroutine Find_Surface()   ! ��������� ������� ������������, ������� ����������
    use STORAGE
    use GEO_PARAM
    implicit none
    integer :: j, k, node
    
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
    
    if (par_developer_info) print *, "Find_Surface: ", node - 1, size(gl_Contact)
    
    end subroutine Find_Surface
    
    ! ************************************************************************************************************************************************
    ! ���� ������
    
    subroutine Initial_conditions()  ! ����� ��������� �������
    use STORAGE
    use GEO_PARAM
    implicit none
    
    integer(4) :: ncell, gr
    real(8) :: r
    real(8) :: ro, P_E
    real(8) :: c(3)
    
    ncell = size(gl_all_Cell(1, :))
    
    !P_E = par_kk/(par_chi_real**2 * par_R0**2) * par_chi_real**2 / (ggg * 5.0**2)
    
    do gr = 1, ncell
        c = gl_Cell_center(:, gr)
        r = norm2(c)
        if (r <= par_R0 * par_R_int * 1.01) then   !(.FALSE.) then!
            ro = par_kk/(par_chi_real**2 * r**2)
            P_E = ro * par_chi_real**2 / (ggg * 7.0**2)
            c = DBLE(c) * par_chi_real/DBLE(r)
            gl_Cell_par(:, gr) = (/ro, DBLE(c(1)), DBLE(c(2)), DBLE(c(3)), P_E, 0.0_8, 0.0_8, 0.0_8/)
        else
            !gl_Cell_par(:, gr) = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 1.0_8, 0.0_8, 0.0_8, 0.0_8/)
        end if
    end do

    
    end subroutine Initial_conditions
    
    subroutine calc_all_Gran()   ! ��������� ������� �������� � �������� ������ � ������� �����
    ! ��� ������� ����� ������ � ����������� ��������� (�.�. � ��������� ������ ������� ����� ��� �� ������� "������"
    ! ���� � ����� ������������ � �������� ��������� �����
    use STORAGE
    use GEO_PARAM
    implicit none
    
    integer, automatic :: Ngran, iter
    real(8) :: p(3, 4), Vol, D
    real(8) :: a(3), b(3), c(3), S, node1(3), node2(3)
    real(8) :: dist, di, gr_center(3)
    integer :: i, j, k, ll, grc
    
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
            p(:,j) = (/gl_x(i), gl_y(i), gl_z(i)/)
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
        
        gl_Gran_center(:, iter) = gr_center
        
        c(1) = a(2) * b(3) - a(3) * b(2) 
        c(2) = a(3) * b(1) - a(1) * b(3) 
        c(3) = a(1) * b(2) - a(2) * b(1) 
        
        S = norm2(C)  ! S = S/2
        c = c/S
        S = S/2.0
        
        ! ����� ���� ��� ���������, ��������� �� ������������� �������!\
        
        call Get_center(gl_Gran_neighbour(1, iter), node1)
        node2 = (p(:,1) + p(:,2) + p(:,3) + p(:,4))/4.0
        node1 = node2 - node1
        if(DOT_PRODUCT(node1, c) < 0.0) then
            print*, "ERROR 1333 ", DOT_PRODUCT(node1, c)
            pause
        end if
        
        ! ����� �������� ������� ����� � ������� � ����� ������!
        gl_Gran_normal(:, iter) = c
        gl_Gran_square(iter) = S
        
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
                c = c + (/gl_x(i), gl_y(i), gl_z(i)/)
            end do
            c = c/8.0
        else
            if (gl_all_Cell(1, iter) == gl_all_Cell(2, iter)) then
                if (gl_all_Cell(4, iter) == gl_all_Cell(8, iter)) then
                    do j = 1,8
                        if (j == 2 .or. j == 5 .or. j == 6 .or. j == 8) CYCLE
        	            i = gl_all_Cell(j, iter)
                        c = c + (/gl_x(i), gl_y(i), gl_z(i)/)
                    end do
                    c = c/4.0
                else
                    do j = 1,8
                        if (j == 2 .or. j == 5 .or. j == 6) CYCLE
        	            i = gl_all_Cell(j, iter)
                        c = c + (/gl_x(i), gl_y(i), gl_z(i)/)
                    end do
                    c = c/5.0
                end if
            else
                if (gl_all_Cell(2, iter) == gl_all_Cell(6, iter)) then
                    if (gl_all_Cell(4, iter) == gl_all_Cell(5, iter)) then
                        do j = 1,8
                            if (j == 4 .or. j == 5 .or. j == 6 .or. j == 8) CYCLE
        	                i = gl_all_Cell(j, iter)
                            c = c + (/gl_x(i), gl_y(i), gl_z(i)/)
                        end do
                        c = c/4.0
                    else
                        do j = 1,8
                            if (j == 5 .or. j == 6) CYCLE
        	                i = gl_all_Cell(j, iter)
                            c = c + (/gl_x(i), gl_y(i), gl_z(i)/)
                        end do
                        c = c/6.0
                    end if
                else
                    if (gl_all_Cell(4, iter) == gl_all_Cell(5, iter)) then
                        do j = 1,8
                            if (j == 4 .or. j == 5 .or. j == 8) CYCLE
        	                i = gl_all_Cell(j, iter)
                            c = c + (/gl_x(i), gl_y(i), gl_z(i)/)
                        end do
                        c = c/5.0
                    else
                        do j = 1,8
                            if (j == 5 .or. j == 8) CYCLE
        	                i = gl_all_Cell(j, iter)
                            c = c + (/gl_x(i), gl_y(i), gl_z(i)/)
                        end do
                        c = c/6.0
                    end if
                end if
            end if
        end if
        
        gl_Cell_center(:, iter) = c
        
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
            b = (/gl_x(k), gl_y(k), gl_z(k)/)
            a = gl_Gran_normal(:,i)
            di = abs(DOT_PRODUCT(c,a) - DOT_PRODUCT(a,b))  ! ���������� �� ����� �� ���������
            Vol = Vol + gl_Gran_square(i) * di/3.0
            dist = MIN(dist, di)
            ll = ll + 1
        end do
        
        
        
        !if (iter == gl_Cell_A(5, 1, 1)) then
        !    write(*, *) c
        !    print*, Vol
        !    print*, gl_Cell_gran(:, iter)
        !    pause
        !end if
        
        
        gl_Cell_Volume(iter) = Vol
        gl_Cell_dist(iter) = dist
        
!        print*, Vol, ll
!        pause
        
    end do
    
    
    end subroutine calc_all_Gran
    
    subroutine Start_GD_1(steps)
    ! ��������� ������ ������� �������� (��� ������������, ��� ��������� �����)
    ! ��� ������ � �.�. ��� ���
    ! ������� �������������� �� OPEN_MP
    use STORAGE
    use GEO_PARAM
    USE OMP_LIB
    implicit none
    
    integer(4), intent(in) :: steps
    integer(4) :: st, gr, ngran, ncell, s1, s2, i, j, k
    real(8) :: qqq1(8), qqq2(8), qqq(8)  ! ���������� � ������
    real(8) :: dist, dsl, dsc, dsp
    real(8) :: POTOK(8)
    real(8) :: time, Volume, TT, U8, rad1, rad2, aa, bb, cc
    
    real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3
    
    ngran = size(gl_all_Gran(1, :))
    ncell = size(gl_all_Cell(1, :))
    
    call calc_all_Gran()   ! ��������� ��� ������, ������, ������� � �.�.
    call Initial_conditions()  ! ����� ��������� �������
    
    time = 0.00000001
    
    !k = gl_Cell_C(1, size(gl_Cell_C(1,:,1)) ,1)
    !print*, gl_Cell_Volume(k)
    !
    !k = gl_Cell_C(1, size(gl_Cell_C(1,:,1)) - 1 ,1)
    !print*, gl_Cell_Volume(k)
    !
    !pause
    
    
    do st = 1, steps
        ! ������ ���� �� ������ � ������� ������ ����� ���
        TT = time
        time = 100000.0
        if (mod(st, 50) == 0)  print*, st, TT
        
        !$omp parallel
        
        !$omp do private(POTOK, s1, s2, qqq1, qqq2, dist, dsl, dsp, dsc, rad1, rad2, aa, bb, cc) &
        !$omp & reduction(min:time)
    	    do gr = 1, ngran
            POTOK = 0.0
    	    s1 = gl_Gran_neighbour(1, gr)
            if(s1 < 1) print*, "ERROR 1465gdfdyhe"
            s2 = gl_Gran_neighbour(2, gr)
            qqq1 = gl_Cell_par(:, s1)
            
            ! ��������� ������ ��������� ��������������� ��������
            if(norm2(qqq1(2:4))/sqrt(ggg*qqq1(5)/qqq1(1)) > 5.0) then
                rad1 = norm2(gl_Cell_center(:, s1))
                rad2 = norm2(gl_Gran_center(:, gr))
                qqq1(1) = qqq1(1) * rad1**2 / rad2**2
                qqq1(5) = qqq1(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
                ! �������� ������ � ����������� �.�.
                call spherical_skorost(gl_Cell_center(1, s1), gl_Cell_center(2, s1), gl_Cell_center(3, s1), &
                    qqq1(2), qqq1(3), qqq1(4), aa, bb, cc)
                call dekard_skorost(gl_Gran_center(1, gr), gl_Gran_center(2, gr), gl_Gran_center(3, gr), &
                    aa, bb, cc, qqq1(2), qqq1(3), qqq1(4))
                
            end if
            
            if (s2 >= 1) then
                if ( norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int .and. norm2(gl_Cell_center(:, s2)) <= par_R0 * par_R_int) CYCLE
                qqq2 = gl_Cell_par(:, s2 )
                dist = min(gl_Cell_dist(s1), gl_Cell_dist(s2))
                
                ! ��������� ������ ��������� ��������������� ��������
                if(norm2(qqq2(2:4))/sqrt(ggg*qqq2(5)/qqq2(1)) > 5.0) then
                    rad1 = norm2(gl_Cell_center(:, s2))
                    rad2 = norm2(gl_Gran_center(:, gr))
                    qqq2(1) = qqq2(1) * rad1**2 / rad2**2
                    qqq2(5) = qqq2(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
                    call spherical_skorost(gl_Cell_center(1, s2), gl_Cell_center(2, s2), gl_Cell_center(3, s2), &
                        qqq2(2), qqq2(3), qqq2(4), aa, bb, cc)
                    call dekard_skorost(gl_Gran_center(1, gr), gl_Gran_center(2, gr), gl_Gran_center(3, gr), &
                        aa, bb, cc, qqq2(2), qqq2(3), qqq2(4))
                end if
                
            else  ! � ������ ��������� ����� - ��������� �������
                if (norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int) CYCLE
                if(s2 == -1) then  ! ���������� �����
                    dist = gl_Cell_dist(s1)
                    qqq2 = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 1.0_8, 0.0_8, 0.0_8, 0.0_8/)
                else  ! ����� ����� ������ �������
                    dist = gl_Cell_dist(s1)
                    qqq2 = qqq1
                    qqq2(5) = 1.0_8
                    if(qqq2(2) > par_Velosity_inf) then
                        qqq2(2) = par_Velosity_inf ! ����� �������� 
                    end if
                end if
            end if
            
            !print*, POTOK
            call chlld(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                           0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)
            !print*, POTOK
            !pause
            
            time = min(time, 0.9 * dist/max(abs(dsl), abs(dsp)) )   ! REDUCTION
            gl_Gran_POTOK(:, gr) = POTOK * gl_Gran_square(gr)
        end do
    !$omp end do
    
    !$omp barrier ! ������������� �����
        ! ������ ���� �� �������
    !$omp do private(POTOK, Volume, qqq, i, j, ro3, u3, v3, w3, p3)
        do gr = 1, ncell
            if (norm2(gl_Cell_center(:, gr)) <= par_R0 * par_R_int) CYCLE    ! �� ������� ������ �����
            POTOK = 0.0
            Volume = gl_Cell_Volume(gr)
            qqq = gl_Cell_par(:, gr)            
            ! ������������ ������ ����� �����
            do i = 1, 6
            	j = gl_Cell_gran(i, gr)
                if (j == 0) CYCLE
                if (gl_Gran_neighbour(1, j) == gr) then
                    POTOK = POTOK + gl_Gran_POTOK(:, j)
                else
                    POTOK = POTOK - gl_Gran_POTOK(:, j)
                end if
            end do
            
            ro3 = qqq(1) - TT * POTOK(1) / Volume
            if (ro3 <= 0.0_8) then
                print*, "Ro < 0  1490 ", ro3
                pause
            end if
            u3 = (qqq(1) * qqq(2) - TT * POTOK(2) / Volume) / ro3
            v3 = (qqq(1) * qqq(3) - TT * POTOK(3) / Volume) / ro3
            w3 = (qqq(1) * qqq(4) - TT * POTOK(4) / Volume) / ro3
            !bx3 = (bx * Volume_do / Volume - T[now1] * (K->Potok[4] + qqq(2) * K->Potok[8]) / Volume)
            !by3 = (by * Volume_do / Volume - T[now1] * (K->Potok[5] + qqq(3) * K->Potok[8]) / Volume)
            !bz3 = (bz * Volume_do / Volume - T[now1] * (K->Potok[6] + qqq(4) * K->Potok[8]) / Volume)
            !U8 = ( qqq(5) / (ggg - 1.0) + 0.5 * qqq(1) * norm2(qqq(2:4))**2 )
            p3 = ((  ( qqq(5) / (ggg - 1.0) + 0.5 * qqq(1) * norm2(qqq(2:4))**2 ) &
                - TT * POTOK(5)/ Volume) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) ) * (ggg - 1.0)
            
            if (p3 <= 0.0_8) then
                print*, "p < 0  1490 ", p3 , gl_Cell_center(:, gr)
                p3 = 0.000001
                !pause
            end if
            
            gl_Cell_par(:, gr) = (/ro3, u3, v3, w3, p3, 0.0_8, 0.0_8, 0.0_8/)
            
        end do
    !$omp end do
        
    !$omp end parallel
    
    end do
    
    end subroutine Start_GD_1
    
    
    ! ************************************************************************************************************************************************
    ! ���� ���������, ����������� ��� ������ - ������� �������, ������
    
    real(8) pure function calc_square(n) ! ��������� ������� ����� ��� ������� n � ����� ������ ������
    use STORAGE, only: gl_all_Gran, gl_x, gl_y, gl_z
    implicit none
    
    integer, intent (in) :: n
    real(8) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4
    real(8) :: v1(3), v2(3), d1, d2
    
    x1 = gl_x(gl_all_Gran(1, n))
    y1 = gl_y(gl_all_Gran(1, n))
    z1 = gl_z(gl_all_Gran(1, n))
    
    x2 = gl_x(gl_all_Gran(2, n))
    y2 = gl_y(gl_all_Gran(2, n))
    z2 = gl_z(gl_all_Gran(2, n))
    
    x3 = gl_x(gl_all_Gran(3, n))
    y3 = gl_y(gl_all_Gran(3, n))
    z3 = gl_z(gl_all_Gran(3, n))
    
    x4 = gl_x(gl_all_Gran(4, n))
    y4 = gl_y(gl_all_Gran(4, n))
    z4 = gl_z(gl_all_Gran(4, n))
    
    v1 = (/x3-x1, y3-y1, z3-z1/)
    v2 = (/x4-x2, y4-y2, z4-z2/)
    
    d1 = norm2(v1)
    d2 = norm2(v2)
    
    calc_square = 0.5 * d1 * d2 * sqrt(1.0 - (DOT_PRODUCT(v1, v2)/d1/d2)**2 )
    end function calc_square
    
    pure subroutine Get_center(n, CENTER)  ! �������� ������ - ����� ������
    use STORAGE, only: gl_all_Cell, gl_x, gl_y, gl_z
    implicit none
    
    integer, intent(in) :: n
    real(8), intent(out) :: CENTER(3)
    real(8) :: c(3)
    integer :: i, j
    
    CENTER = 0.0
    
    do i = 1, 8
        j = gl_all_Cell(i, n)
    	c = (/gl_x(j), gl_y(j), gl_z(j)/)
        CENTER = CENTER + c
    end do
    
    CENTER = CENTER/8.0
    
    end subroutine Get_center
    
    ! ****************************************************************************************************************************************************
    ! ���� ������� ��� ������ � �������� ��������� �����
    
    subroutine Print_all_surface(str)  ! ������ ������������, ������� ���������� (TS, HP, BS) 
    use GEO_PARAM
    use STORAGE
    use My_func
    implicit none
    
    !interface
    !    real(8) pure function calc_square(n)
    !    integer, intent(in) :: n
    !    end function
    !end interface
    
    character, intent(in) :: str
    integer :: m, n, j, i
    
    if (str == "C") then ! Contact
         m = size(gl_Contact)
         open(1, file = 'Contact_geometry.txt')  
         write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z', 'M'  ZONE T= 'HP', N= ", 4 * (m), ", E =  ", 4 * (m) , ", F=FEPOINT, ET=LINESEG "
         do j = 1, m
             do i = 1, 4
                 n = gl_Contact(j)
                 write(1,*) gl_x(gl_all_Gran(i, n)), gl_y(gl_all_Gran(i, n)), gl_z(gl_all_Gran(i, n)), calc_square(n)
             end do
         end do
         
         do j = 0, m
             write(1,*) 4 * j + 1, 4 * j + 2
             write(1,*) 4 * j + 2, 4 * j + 3
             write(1,*) 4 * j + 3, 4 * j + 4
             write(1,*) 4 * j + 4, 4 * j + 1
         end do
         close(1)
    else
        print*, "Print_all_surface: ERROR 1304"
    end if
    
    end subroutine Print_all_surface
    
    subroutine Print_Setka_2D()  ! �������� 2� ����� � ������� � �������
    use GEO_PARAM
    use STORAGE
    implicit none
    
    integer, automatic :: N1, N2, kk, k, i, j, N
    
    
    N = size(gl_Cell_A(1, :, 1)) * size(gl_Cell_A(:, 1, 1)) + & 
        size(gl_Cell_B(1, :, 1)) * size(gl_Cell_B(:, 1, 1)) + &
        size(gl_Cell_C(1, :, 1)) * size(gl_Cell_C(:, 1, 1))
    N = N * 2
    
    open(1, file = 'print_setka_2D.txt')  
    write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y'  ZONE T= 'HP', N= ", 4 * N, ", E =  ", 4 * N , ", F=FEPOINT, ET=LINESEG "
    
    
    kk = 1
    N2 = size(gl_Cell_A(1, :, 1))
    N1 = size(gl_Cell_A(:, 1, 1))
    do j = 1, N2
        do i = 1, N1
            do k = 1, 4
                write(1,*) gl_x(gl_all_Cell(k, gl_Cell_A(i, j, kk) )), gl_y(gl_all_Cell(k, gl_Cell_A(i, j, kk) ))
            end do
        end do
    end do
    N2 = size(gl_Cell_B(1, :, 1))
    N1 = size(gl_Cell_B(:, 1, 1))
    do j = 1, N2
        do i = 1, N1
            do k = 1, 4
                write(1,*) gl_x(gl_all_Cell(k, gl_Cell_B(i, j, kk) )), gl_y(gl_all_Cell(k, gl_Cell_B(i, j, kk) ))
            end do
        end do
    end do
    N2 = size(gl_Cell_C(1, :, 1))
    N1 = size(gl_Cell_C(:, 1, 1))
    do j = 1, N2
        do i = 1, N1
            do k = 1, 4
                write(1,*) gl_x(gl_all_Cell(k, gl_Cell_C(i, j, kk) )), gl_y(gl_all_Cell(k, gl_Cell_C(i, j, kk) ))
            end do
        end do
    end do
    
    ! ������ ����� �����
    
    kk = par_l_phi/2 + 1
    N2 = size(gl_Cell_A(1, :, 1))
    N1 = size(gl_Cell_A(:, 1, 1))
    do j = 1, N2
        do i = 1, N1
            do k = 1, 4
                write(1,*) gl_x(gl_all_Cell(k, gl_Cell_A(i, j, kk) )), gl_y(gl_all_Cell(k, gl_Cell_A(i, j, kk) ))
            end do
        end do
    end do
    N2 = size(gl_Cell_B(1, :, 1))
    N1 = size(gl_Cell_B(:, 1, 1))
    do j = 1, N2
        do i = 1, N1
            do k = 1, 4
                write(1,*) gl_x(gl_all_Cell(k, gl_Cell_B(i, j, kk) )), gl_y(gl_all_Cell(k, gl_Cell_B(i, j, kk) ))
            end do
        end do
    end do
    N2 = size(gl_Cell_C(1, :, 1))
    N1 = size(gl_Cell_C(:, 1, 1))
    do j = 1, N2
        do i = 1, N1
            do k = 1, 4
                write(1,*) gl_x(gl_all_Cell(k, gl_Cell_C(i, j, kk) )), gl_y(gl_all_Cell(k, gl_Cell_C(i, j, kk) ))
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
        
    
    end subroutine Print_Setka_2D
    
    subroutine Print_par_2D()  ! �������� 2� ����� � ������� � �������
    use GEO_PARAM
    use STORAGE
    implicit none
    
    integer, automatic :: N1, N2, kk, k, i, j, N, m
    real(8) :: c(3), Mach
    
    
    N = size(gl_Cell_A(1, :, 1)) * size(gl_Cell_A(:, 1, 1)) + & 
        size(gl_Cell_B(1, :, 1)) * size(gl_Cell_B(:, 1, 1)) + &
        size(gl_Cell_C(1, :, 1)) * size(gl_Cell_C(:, 1, 1))
    N = N * 2
    
    open(1, file = 'print_par_2D.txt')  
    write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z', 'rho', 'u', 'v', 'w', 'p', 'bx', 'by', 'bz', 'Volume', 'Mach'  ZONE T= 'HP'"
    
    
    kk = 1
    N2 = size(gl_Cell_A(1, :, 1))
    N1 = size(gl_Cell_A(:, 1, 1))
    do j = 1, N2
        do i = 1, N1
            c = gl_Cell_center(:, gl_Cell_A(i, j, kk))
            m = gl_Cell_A(i, j, kk)
            Mach = norm2(gl_Cell_par(2:4, m ))/sqrt(ggg*gl_Cell_par(5, m )/gl_Cell_par(1, m ))
            write(1,*) c, gl_Cell_par(:, m ), gl_Cell_Volume(m), Mach
        end do
    end do
    N2 = size(gl_Cell_B(1, :, 1))
    N1 = size(gl_Cell_B(:, 1, 1))
    do j = 1, N2
        do i = 1, N1
            m = gl_Cell_B(i, j, kk)
            c = gl_Cell_center(:, m)
            Mach = norm2(gl_Cell_par(2:4, m ))/sqrt(ggg*gl_Cell_par(5, m )/gl_Cell_par(1, m ))
            write(1,*) c, gl_Cell_par(:, m ), gl_Cell_Volume(m), Mach
        end do
    end do
    N2 = size(gl_Cell_C(1, :, 1))
    N1 = size(gl_Cell_C(:, 1, 1))
    do j = 1, N2
        do i = 1, N1
            m = gl_Cell_C(i, j, kk)
            c = gl_Cell_center(:, m)
            Mach = norm2(gl_Cell_par(2:4, m ))/sqrt(ggg*gl_Cell_par(5, m )/gl_Cell_par(1, m ))
            write(1,*) c, gl_Cell_par(:, m ), gl_Cell_Volume(m), Mach
        end do
    end do
    
    ! ������ ����� �����
    
    kk = par_l_phi/2 + 1
    N2 = size(gl_Cell_A(1, :, 1))
    N1 = size(gl_Cell_A(:, 1, 1))
    do j = 1, N2
        do i = 1, N1
            m = gl_Cell_A(i, j, kk)
            c = gl_Cell_center(:, m)
            Mach = norm2(gl_Cell_par(2:4, m ))/sqrt(ggg*gl_Cell_par(5, m )/gl_Cell_par(1, m ))
            write(1,*) c, gl_Cell_par(:, m ), gl_Cell_Volume(m), Mach
        end do
    end do
    N2 = size(gl_Cell_B(1, :, 1))
    N1 = size(gl_Cell_B(:, 1, 1))
    do j = 1, N2
        do i = 1, N1
            m = gl_Cell_B(i, j, kk)
            c = gl_Cell_center(:, m)
            Mach = norm2(gl_Cell_par(2:4, m ))/sqrt(ggg*gl_Cell_par(5, m )/gl_Cell_par(1, m ))
            write(1,*) c, gl_Cell_par(:, m ), gl_Cell_Volume(m), Mach
        end do
    end do
    N2 = size(gl_Cell_C(1, :, 1))
    N1 = size(gl_Cell_C(:, 1, 1))
    do j = 1, N2
        do i = 1, N1
            m = gl_Cell_C(i, j, kk)
            c = gl_Cell_center(:, m)
            Mach = norm2(gl_Cell_par(2:4, m ))/sqrt(ggg*gl_Cell_par(5, m )/gl_Cell_par(1, m ))
            write(1,*) c, gl_Cell_par(:, m ), gl_Cell_Volume(m), Mach
        end do
    end do
    
    
    close(1)
        
    
    end subroutine Print_par_2D
    
    subroutine Print_Point_Plane()
    ! ����� � ��������� XOY - ����� ��� ����������
    ! ������� ��� �������� ���������
    use GEO_PARAM
    use STORAGE
    implicit none
    
    integer :: N1, N2, k, l, i, j
    
    open(1, file = 'all_point_plane.txt')
    
    N2 = size(gl_RAY_A(1, :, 1))
    N1 = size(gl_RAY_A(:, 1, 1))
    
    k = 1
    do j = 1, N2
        do i = 1, N1
            l = gl_RAY_A(i, j, k)
            write(1,*) gl_x(l), gl_y(l)
        end do
    end do
    
    k = par_l_phi/2 + 1
    do j = 1, N2
        do i = 1, N1
            l = gl_RAY_A(i, j, k)
            write(1,*) gl_x(l), gl_y(l)
        end do
    end do
    
    ! ----------------------------------------------------- B
    N2 = size(gl_RAY_B(1, :, 1))
    N1 = size(gl_RAY_B(:, 1, 1))
    
    k = 1
    do j = 1, N2
        do i = 1, N1
            l = gl_RAY_B(i, j, k)
            write(1,*) gl_x(l), gl_y(l)
        end do
    end do
    
    k = par_l_phi/2 + 1
    do j = 1, N2
        do i = 1, N1
            l = gl_RAY_B(i, j, k)
            write(1,*) gl_x(l), gl_y(l)
        end do
    end do
    
    ! ----------------------------------------------------- C
    N2 = size(gl_RAY_C(1, :, 1))
    N1 = size(gl_RAY_C(:, 1, 1))
    
    k = 1
    do j = 1, N2
        do i = 1, N1
            l = gl_RAY_C(i, j, k)
            write(1,*) gl_x(l), gl_y(l)
        end do
    end do
    
    k = par_l_phi/2 + 1
    do j = 1, N2
        do i = 1, N1
            l = gl_RAY_C(i, j, k)
            write(1,*) gl_x(l), gl_y(l)
        end do
    end do
    
    ! ----------------------------------------------------- O
    N2 = size(gl_RAY_O(1, :, 1))
    N1 = size(gl_RAY_O(:, 1, 1))
    
    k = 1
    do j = 1, N2
        do i = 1, N1
            l = gl_RAY_O(i, j, k)
            write(1,*) gl_x(l), gl_y(l)
        end do
    end do
    
    k = par_l_phi/2 + 1
    do j = 1, N2
        do i = 1, N1
            l = gl_RAY_O(i, j, k)
            write(1,*) gl_x(l), gl_y(l)
        end do
    end do
    
    ! ----------------------------------------------------- K
    N2 = size(gl_RAY_K(1, :, 1))
    N1 = size(gl_RAY_K(:, 1, 1))
    
    k = 1
    do j = 1, N2
        do i = 1, N1
            l = gl_RAY_K(i, j, k)
            write(1,*) gl_x(l), gl_y(l)
        end do
    end do
    
    k = par_l_phi/2 + 1
    do j = 1, N2
        do i = 1, N1
            l = gl_RAY_K(i, j, k)
            write(1,*) gl_x(l), gl_y(l)
        end do
    end do
    
    ! ----------------------------------------------------- D
    N2 = size(gl_RAY_D(1, :, 1))
    N1 = size(gl_RAY_D(:, 1, 1))
    
    k = 1
    do j = 1, N2
        do i = 1, N1
            l = gl_RAY_D(i, j, k)
            write(1,*) gl_x(l), gl_y(l)
        end do
    end do
    
    k = par_l_phi/2 + 1
    do j = 1, N2
        do i = 1, N1
            l = gl_RAY_D(i, j, k)
            write(1,*) gl_x(l), gl_y(l)
        end do
    end do
    
    ! ----------------------------------------------------- E
    N2 = size(gl_RAY_E(1, :, 1))
    N1 = size(gl_RAY_E(:, 1, 1))
    
    k = 1
    do j = 1, N2
        do i = 1, N1
            l = gl_RAY_E(i, j, k)
            write(1,*) gl_x(l), gl_y(l)
        end do
    end do
    
    k = par_l_phi/2 + 1
    do j = 1, N2
        do i = 1, N1
            l = gl_RAY_E(i, j, k)
            write(1,*) gl_x(l), gl_y(l)
        end do
    end do
    
    close(1)
    
    end subroutine Print_Point_Plane
    
    subroutine Print_All_Points()   ! ������� ������ �������� ��� ����� ����� (����� ����� ����������� �� � ��������)
    ! Variables
    use STORAGE
    
    open(1, file = 'all_point.txt')  
    
    do i=1,size(gl_x)  
      write(1,*) gl_x(i), gl_y(i), gl_z(i)   
    end do  
    
    close(1)
    
    end subroutine Print_All_Points
    
    subroutine Print_A_Ray(n1, n2)   ! ������� ������ �������� ���� ��� (� ���������, � ������������) + ������ �������� �����
    
    use GEO_PARAM
    use STORAGE
    
    implicit none
    integer, intent(in) :: n1, n2
    integer(4), automatic :: i, k
    
    open(1, file = 'Ray.txt')  
    
    do i = 1, size(gl_RAY_A(:, n1, n2))
        k = 0
        if (i == par_n_TS .or. i == par_n_HP .or. i == par_n_BS .or. i == par_n_END) k = 1
        write(1,*) gl_x(gl_RAY_A(i, n1, n2)), gl_y(gl_RAY_A(i, n1, n2)), gl_z(gl_RAY_A(i, n1, n2)),  k 
    end do  
    
    close(1)
    
    end subroutine Print_A_Ray
    
    subroutine Print_Cell(n1, n2, n3, str)  ! �������� ������ � �������� n1, n2, n3 �� ��������� ����� str
    use GEO_PARAM
    use STORAGE
    
    implicit none
    integer, intent(in) :: n1, n2, n3
    character, intent(in) :: str
    integer(4), automatic :: i, k
    
     open(1, file = 'one_Cell.txt')  
    
    
    write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z'  ZONE T= 'HP', N= ", 8, ", E =  ", 12 , ", F=FEPOINT, ET=LINESEG "
    
    do i = 1, 8
        if (str == "A") then
            if (n1 > size(gl_Cell_A(:, 1, 1))) PAUSE "n1 too mach! ERROR 822"
            if (n2 > size(gl_Cell_A(1, :, 1))) PAUSE "n2 too mach! ERROR 822"
            if (n3 > size(gl_Cell_A(1, 1, :))) PAUSE "n3 too mach! ERROR 822"
            write(1,*) gl_x(gl_all_Cell(i, gl_Cell_A(n1, n2, n3))), gl_y(gl_all_Cell(i, gl_Cell_A(n1, n2, n3))), gl_z(gl_all_Cell(i, gl_Cell_A(n1, n2, n3)))
        else if (str == "B") then
            if (n1 > size(gl_Cell_B(:, 1, 1))) PAUSE "n1 too mach! ERROR 827"
            if (n2 > size(gl_Cell_B(1, :, 1))) PAUSE "n2 too mach! ERROR 827"
            if (n3 > size(gl_Cell_B(1, 1, :))) PAUSE "n3 too mach! ERROR 827"
            write(1,*) gl_x(gl_all_Cell(i, gl_Cell_B(n1, n2, n3))), gl_y(gl_all_Cell(i, gl_Cell_B(n1, n2, n3))), gl_z(gl_all_Cell(i, gl_Cell_B(n1, n2, n3)))
        else
            if (n1 > size(gl_Cell_C(:, 1, 1))) PAUSE "n1 too mach! ERROR 827"
            if (n2 > size(gl_Cell_C(1, :, 1))) PAUSE "n2 too mach! ERROR 827"
            if (n3 > size(gl_Cell_C(1, 1, :))) PAUSE "n3 too mach! ERROR 827"
            write(1,*) gl_x(gl_all_Cell(i, gl_Cell_C(n1, n2, n3))), gl_y(gl_all_Cell(i, gl_Cell_C(n1, n2, n3))), gl_z(gl_all_Cell(i, gl_Cell_C(n1, n2, n3)))
        end if
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
    
    end subroutine Print_Cell
       
    subroutine Print_all_Cell()  ! �������� ��� ������ ����� (� ������) - ���� ������ � �����-�� ��������
    use GEO_PARAM
    use STORAGE
    implicit none
    
    integer(4) :: N1, i, j
    
    N1 = size(gl_Cell_A(:, :, :)) + size(gl_Cell_B(:, :, :)) + size(gl_Cell_C(:, :, :))
    
    open(1, file = 'all_Cell.txt')  
    write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z'  ZONE T= 'HP', N= ", 8 * N1, ", E =  ", 12 * N1 , ", F=FEPOINT, ET=LINESEG "
    
    
    
    do j = 1, N1
            do i = 1, 8
                write(1,*) gl_x(gl_all_Cell(i,j)), gl_y(gl_all_Cell(i,j)), gl_z(gl_all_Cell(i,j))
            end do
    end do
    
    do j = 0, N1-1
        write(1,*) 8 * j + 1, 8 * j + 2
        write(1,*) 8 * j + 2, 8 * j + 3
        write(1,*) 8 * j + 3, 8 * j + 4
        write(1,*) 8 * j + 4, 8 * j + 1
        write(1,*) 8 * j + 1, 8 * j + 5
        write(1,*) 8 * j + 2, 8 * j + 6
        write(1,*) 8 * j + 3, 8 * j + 7
        write(1,*) 8 * j + 4, 8 * j + 8
        write(1,*) 8 * j + 5, 8 * j + 6
        write(1,*) 8 * j + 6, 8 * j + 7
        write(1,*) 8 * j + 7, 8 * j + 8
        write(1,*) 8 * j + 8, 8 * j + 5
    end do
    
    
    close(1)
                
    
    end subroutine Print_all_Cell
    
    subroutine Print_cell_and_neighbour(n1, n2, n3)  ! �������� ���� ������� ������ � ���� ������
    use GEO_PARAM
    use STORAGE
    implicit none
    
    integer, intent(in) :: n1, n2, n3
    integer, automatic :: n, m, i, j, k
    
    n = gl_Cell_C(n1, n2, n3)  ! �������� ����� ������ � ����� ������
    
    m = COUNT(gl_Cell_neighbour(:, n) > 0)
    
    if (par_developer_info) print *, "Kolichestvo sosedey = ", m
    
    open(1, file = 'Cell+neighbour.txt')  
    
    
    write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z', 'M'  ZONE T= 'HP', N= ", 8 * (m + 1), ", E =  ", 12 * (m + 1) , ", F=FEPOINT, ET=LINESEG "
    
    ! �������� ���� ������
    do i = 1, 8
        write(1,*) gl_x(gl_all_Cell(i, n)), gl_y(gl_all_Cell(i, n)), gl_z(gl_all_Cell(i, n)), 0
    end do
    
    ! �������� �������
    do j = 1, 6
        if (gl_Cell_neighbour(j, n) <= 0) CYCLE
        k = gl_Cell_neighbour(j, n)
        do i = 1, 8
        write(1,*) gl_x(gl_all_Cell(i, k)), gl_y(gl_all_Cell(i, k)), gl_z(gl_all_Cell(i, k)), j
        end do
    end do
    
    
    do j = 0, m
        write(1,*) 8 * j + 1, 8 * j + 2
        write(1,*) 8 * j + 2, 8 * j + 3
        write(1,*) 8 * j + 3, 8 * j + 4
        write(1,*) 8 * j + 4, 8 * j + 1
        write(1,*) 8 * j + 1, 8 * j + 5
        write(1,*) 8 * j + 2, 8 * j + 6
        write(1,*) 8 * j + 3, 8 * j + 7
        write(1,*) 8 * j + 4, 8 * j + 8
        write(1,*) 8 * j + 5, 8 * j + 6
        write(1,*) 8 * j + 6, 8 * j + 7
        write(1,*) 8 * j + 7, 8 * j + 8
        write(1,*) 8 * j + 8, 8 * j + 5
    end do
    
    close(1)
    
    
    end subroutine Print_cell_and_neighbour
    
    subroutine Print_gran()  ! �������� ����� ����� (����� � �����-�� ��������)
    use GEO_PARAM
    use STORAGE
    implicit none

    integer, automatic :: n, m, i, j, k
    
    m = 0
    
    !do i = 1, size(gl_all_Gran(1,:))
    !    if ( ANY( gl_Gran_neighbour(:,i) == -2 ) ) m = m + 1   ! ��� ������� ����� �������� � ���� (� ����� ����� ������)
    !end do
    m = 1
    
    if (par_developer_info) print *, "Kolichestvo graney = ", m
    
    open(1, file = 'Grans.txt')  
    
    
    write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z', 'M'  ZONE T= 'HP', N= ", 4 * (m), ", E =  ", 4 * (m) , ", F=FEPOINT, ET=LINESEG "
    
    ! ��������
    do n = 3953, 3953 !1, size(gl_all_Gran(1,:))
        !if ( ANY( gl_Gran_neighbour(:,n) == -2 ) == .FALSE. ) CYCLE
        do i = 1, 4
            write(1,*) gl_x(gl_all_Gran(i, n)), gl_y(gl_all_Gran(i, n)), gl_z(gl_all_Gran(i, n)), n
        end do
    end do
    
    
    do j = 0, m-1
        write(1,*) 4 * j + 1, 4 * j + 2
        write(1,*) 4 * j + 2, 4 * j + 3
        write(1,*) 4 * j + 3, 4 * j + 4
        write(1,*) 4 * j + 4, 4 * j + 1
    end do
    
    close(1)
    
    end subroutine Print_gran
    
    subroutine Geometry_check()
    use GEO_PARAM
    use STORAGE
    
    implicit none
    integer(4) :: n1, n2, n3, i, j, k, j1
    integer(4) :: a1, a2, a3
    
    if (par_developer_info) print *, "Geometry_check: Proverka 1"
    ! ��������� ������� �� ���� ������ ��� ����� (�� �������� �� ��������������������)
    ! ----------------------------------------------------- A
    N2 = size(gl_RAY_A(1, :, 1))
    N1 = size(gl_RAY_A(:, 1, 1))
    N3 = size(gl_RAY_A(1, 1, :))
    do k = 1, N3
        do j = 1, N2
            do i = 1, N1
                if (gl_RAY_A(i, j, k) < 1) STOP "Error gl_RAY_A 618"
            end do
        end do
    end do
    
    ! ----------------------------------------------------- B
    N2 = size(gl_RAY_B(1, :, 1))
    N1 = size(gl_RAY_B(:, 1, 1))
    N3 = size(gl_RAY_B(1, 1, :))
    do k = 1, N3
        do j = 1, N2
            do i = 1, N1
                if (gl_RAY_B(i, j, k) < 1) STOP "Error gl_RAY_B 618"
            end do
        end do
    end do
    
    ! ----------------------------------------------------- C
    N2 = size(gl_RAY_C(1, :, 1))
    N1 = size(gl_RAY_C(:, 1, 1))
    N3 = size(gl_RAY_C(1, 1, :))
    do k = 1, N3
        do j = 1, N2
            do i = 1, N1
                if (gl_RAY_C(i, j, k) < 1) STOP "Error gl_RAY_C 618"
            end do
        end do
    end do
    
    ! ----------------------------------------------------- O
    N2 = size(gl_RAY_O(1, :, 1))
    N1 = size(gl_RAY_O(:, 1, 1))
    N3 = size(gl_RAY_O(1, 1, :))
    do k = 1, N3
        do j = 1, N2
            do i = 1, N1
                if (gl_RAY_O(i, j, k) < 1) STOP "Error gl_RAY_O 618"
            end do
        end do
    end do
    
    ! ----------------------------------------------------- K
    N2 = size(gl_RAY_K(1, :, 1))
    N1 = size(gl_RAY_K(:, 1, 1))
    N3 = size(gl_RAY_K(1, 1, :))
    do k = 1, N3
        do j = 1, N2
            do i = 1, N1
                if (gl_RAY_K(i, j, k) < 1) STOP "Error gl_RAY_K 618"
            end do
        end do
    end do
    
    ! ----------------------------------------------------- D
    N2 = size(gl_RAY_D(1, :, 1))
    N1 = size(gl_RAY_D(:, 1, 1))
    N3 = size(gl_RAY_D(1, 1, :))
    do k = 1, N3
        do j = 1, N2
            do i = 1, N1
                if (gl_RAY_D(i, j, k) < 1) STOP "Error gl_RAY_D 618"
            end do
        end do
    end do
    
    ! ----------------------------------------------------- E
    N2 = size(gl_RAY_E(1, :, 1))
    N1 = size(gl_RAY_E(:, 1, 1))
    N3 = size(gl_RAY_E(1, 1, :))
    do k = 1, N3
        do j = 1, N2
            do i = 1, N1
                if (gl_RAY_E(i, j, k) < 1) STOP "Error gl_RAY_E 618"
            end do
        end do
    end do
    
    if (par_developer_info) print *, "Geometry_check: Proverka 2"
    ! ��������� ������������ ������� (����� ���� ����������)
    do i = 1, size(gl_all_Cell(1, :))
        do j = 1, 6
            if (gl_Cell_neighbour(j, i) > 0) then
                j1 = gl_Cell_neighbour(j, i)
                if( ANY(gl_Cell_neighbour(:, j1) == i) == .FALSE.) then
                    pause "ERROR 1763"
                end if
            end if
        end do
    end do
    
    if (par_developer_info) print *, "Geometry_check: Proverka 3"
    ! ��������� ������������ ������� ������ (����� ��� ������� ��������� �������, � �� ��� ����)
    do i = 1, size(gl_Gran_neighbour(1, :))
        if (gl_Gran_neighbour(1, i) < 1 .or. gl_Gran_neighbour(2, i) < 1) CYCLE
        j1 = gl_Gran_neighbour(1, i)
        if( ANY(gl_Cell_neighbour(:, j1) == gl_Gran_neighbour(2, i)) == .FALSE.) then
            pause "ERROR 1775"
        end if
    end do
    
    
    if (par_developer_info) print *, "Geometry_check: OK"
    
    ! �����-�� ��������� �������� ����� ���� ���������
    !print *, "B111:"
    !a1 = gl_Cell_B(1,1,1)
    !do a2 = 1, 8 
    !	print*, gl_all_Cell(a2, a1)
    !end do
    
    
    end subroutine Geometry_check
    
    ! ****************************************************************************************************************************************************
    ! �������� ���������
    
    program MIK

    implicit none
    
    real(16) :: fg(5), dd
    
    fg = 10
    fg(1) = 2
    dd = 124532
    
    !call Set_STORAGE()
    !call Build_Mesh_start()
    
    call Read_setka_bin(1)
    
    call Find_Surface()
    call Print_All_Points()
    call Print_A_Ray(1,3)
    call Print_Point_Plane()
    call Geometry_check()
    !call Print_Cell(1, 2, 1, "C")
    !call Print_all_Cell()
    call Print_Setka_2D()
    !call Print_cell_and_neighbour(1,2,1)
    !call Print_gran()
    !call Print_all_surface("C")
    call calc_all_Gran()
    
    print *, "Start_GD_1"
    call Start_GD_1(5000)
    call Print_par_2D()
    call Save_setka_bin(2)
    ! Variables
    
    pause
    end program MIK


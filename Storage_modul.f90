
    
    module GEO_PARAM                     ! ������ �������������� - �������� ���������� - �������� ���������
    
    implicit none
    
    logical, parameter :: par_developer_info = .True.   ! �������� ������������ ��� ������ �������������� ���������
    real(8), parameter :: par_R_character = 35.6505         ! ����������� ������ � ������ (���������� �� TS �� ��������� ����� ���������� �����)
    real(8), parameter :: par_koeff_HP = 1.3            ! �� ������� �������� ����������� ������ ��� ��������� HP
    real(8), parameter :: par_koeff_BS = 2.0            ! �� ������� �������� ����������� ������ ��� ��������� BS
    real(8), parameter :: par_R0 = 0.256505         ! ����������� ������ 1 �.�. (���������� �����)
    ! �������� ����� ��� �������������� �������� ����� eps = par_R_character/10000
    real(8), parameter :: par_R_END = 300.0         ! 
    real(8), parameter :: par_R_LEFT = -390.0         ! 
    real(8), parameter :: par_pi_8 = acos(-1.0_8)         ! ����������� ������ � ������ (��������� ���������� �� TS)
    real(8), parameter :: par_pi_4 = acos(-1.0_4)         ! ����������� ������ � ������ (��������� ���������� �� TS)
    real(8), parameter :: cpi4 = 12.56637061435917295384
    real(8), parameter :: ggg = (5.0/3.0)
    real(8), parameter :: par_kk = 10000.0                ! ������� ������������ ������� ������������
    
    real(8), parameter :: par_a_2 = 0.1307345665         ! �������� � ������� �����������
    real(8), parameter :: par_n_p_LISM = 1.0_8         ! � �����������
    real(8), parameter :: par_n_H_LISM_ = 3.0_8
    real(8), parameter :: par_Kn = 43.3   !0.4326569808         ! � �����������
    
    real(8), parameter :: par_chi_real = 1.0_8      ! � ����� �� ������� �������
    real(8), parameter :: par_chi = 36.1275_8      ! � ����� �� ���� ���� �� �������
    real(8), parameter :: par_Velosity_inf = -2.54338_8
    
    integer, parameter :: par_R_int = 70  ! ������� �.�. �� ������� ������
    
    real(8), parameter :: par_kk1 = 1.7_8     ! ������� �������� ����� � ���� � ������� �� TS: 1 - ��������, 2 - ������������ � �.�.
    real(8), parameter :: par_kk2 = 2.0_8     ! ������� �������� � �������� ������� �� �������������
    real(8), parameter :: par_kk3 = 1.8_8     ! ������� �������� � ������

    integer :: par_l_phi = 48 !4   ! ���������� �������� ��������� ����� ��� ��������� ��������� (���������� �� ���� phi)
                                ! ������ �������� �� 4 ��� �������� ������ ����������� � ����������
    integer :: par_m_A = 30 !5     ! ���������� ����� A � ���������
    integer :: par_m_BC = 18 !3     ! ���������� ����� B/C � ���������
    integer :: par_m_O = 17 !7 !3     ! ���������� ����� O � ���������
    integer :: par_m_K = 7 !2     ! ���������� ����� K � ���������
    real(8) :: par_triple_point = 13.0 * par_pi_4/40.0     ! �� ������ ���� ������� �� pi/2 (� �������������� x) ������� �����
    
    ! ���������� ����� �� ����� A
    integer :: par_n_TS =  26! 3                 ! ���������� ����� �� TS (TS ����������)
    integer :: par_n_HP =  40! 4                 ! ���������� ����� HP (HP ����������)  �� �� 0 ���������
    integer :: par_n_BS =  60! 5                 ! ���������� ����� BS (BS ����������)
    integer :: par_n_END = 72! 6                ! ���������� ����� �� ����� ����� (����� ����������)
    
    integer :: par_n_points  ! ����� ����� � �����
    
    end module GEO_PARAM
    
    
    module STORAGE                       ! ������ ���������� ������ � ����� (��� ���������� ���������� �� gl - global)
    implicit none
    
    real(8), allocatable :: gl_x(:)   ! ����� x-��������� ����� ����� 
    real(8), allocatable :: gl_y(:)   ! ����� y-��������� ����� ����� 
    real(8), allocatable :: gl_z(:)   ! ����� z-��������� ����� ����� 
    ! ���� ����� �������� ��������, �� ��� ����� ��� ������� � ������� ��������� ������ � �������������, ����� �� ������
    
    ! ����, �� ������� ������������� ����� ����� - ��� 
    integer(4), allocatable :: gl_RAY_A(:,:,:)   ! ����� �-����� ����������� 3 (�� ���� ����, � ���� ���������, �� ���� � ������������)
    integer(4), allocatable :: gl_RAY_B(:,:,:)   ! ����� B-����� ����������� 3 (�� ���� ����, � ���� ���������, �� ���� � ������������)
    integer(4), allocatable :: gl_RAY_C(:,:,:)   ! ����� C-����� ����������� 3 (�� ���� ����, � ���� ���������, �� ���� � ������������)
    integer(4), allocatable :: gl_RAY_O(:,:,:)   ! ����� O-����� ����������� 3 (�� ���� ����, � ���� ���������, �� ���� � ������������)
    integer(4), allocatable :: gl_RAY_K(:,:,:)   ! ����� K-����� ����������� 3 (�� ���� ����, � ���� ���������, �� ���� � ������������)
    integer(4), allocatable :: gl_RAY_D(:,:,:)   ! ����� D-����� ����������� 3 (�� ���� ����, � ���� ���������, �� ���� � ������������)
    integer(4), allocatable :: gl_RAY_E(:,:,:)   ! ����� E-����� ����������� 3 (�� ���� ����, � ���� ���������, �� ���� � ������������)
    
    ! ������
    integer(4), allocatable :: gl_Cell_A(:,:,:)   ! ����� A-����� ����������� 4 (�� ���� ����, � ���� ���������, �� ���� � ������������)
    integer(4), allocatable :: gl_Cell_B(:,:,:)   ! ����� B-����� ����������� 4 (�� ���� ����, � ���� ���������, �� ���� � ������������)
    integer(4), allocatable :: gl_Cell_C(:,:,:)   ! ����� C-����� ����������� 4 (�� ���� ����, � ���� ���������, �� ���� � ������������)
    ! ���������� ������ ����� ��� ������� �������������� �����, �������� ������� � �.�.
    ! ������ ��������� �� ������
    
    integer(4), allocatable :: gl_all_Cell(:,:)   ! ���� ����� ����� (8, :) - ������ ���������� ������� - ��� ����� ����� ������
    
    integer(4), allocatable :: gl_all_Cell_inner(:)   ! ������ �����, ������� ��������� ������ (��������� ��������)
    
    integer(4), allocatable :: gl_Cell_neighbour(:,:)   ! ����� �� 6 ������� ��� ������ ������ (���� ����� ������ = 0, �� ��� ��� � ���� �����������)
    integer(4), allocatable :: gl_Cell_gran(:,:)        ! ����� �� 6 ������ ��� ������ ������ (���� ����� = 0, �� ����� ��� � ���� �����������)
    integer(4), allocatable :: gl_Cell_info(:)        ! ����� �� 6 ������ ��� ������ ������ (���� ����� = 0, �� ����� ��� � ���� �����������)
    
    real(8), allocatable :: gl_Cell_Volume(:)           ! ����� ������� �����
    real(8), allocatable :: gl_Cell_dist(:)             ! ����������� ���������� �� ����� � ������ ������  4444444444444444
    real(8), allocatable :: gl_Cell_center(:, :)             ! (3, :) ����� ������ ������  4444444444444444
    real(8), allocatable :: gl_Cell_par(:, :)           ! ����� ���������� (8 ���������)
    real(8), allocatable :: gl_Cell_par_MF(:,:,:)           ! ����� ���������� (5,4,:)  ����������� ��������� (�� 5 ��� ������ �� 4-� ���������)
    
    
    ! ���� ������ ----------------- ����� ��� ������� ��������, �������, �������, �������� � �.�.
    integer(4), allocatable :: gl_all_Gran(:,:)       ! ��� ����� (4,:) ����� �� 4 ����
    integer(4), allocatable :: gl_Gran_neighbour(:,:) ! ������ ������ ����� (2,:) ����� �� 2 ������, ������� ���� �� ������� �� �������
    
    real(8), allocatable :: gl_Gran_normal(:,:)       ! (3, :) ������� �����   4444444444444444
    real(8), allocatable :: gl_Gran_square(:)         ! (:) ������� �����
    real(8), allocatable :: gl_Gran_POTOK(:,:)       ! (8, :) ����� �����
    real(8), allocatable :: gl_Gran_POTOK_MF(:,:,:)       ! (5,4, :) ����� ����� �������������� ���������
    real(8), allocatable :: gl_Gran_center(:,:)       ! (3, :) ����� �����
    integer(4), allocatable :: gl_Gran_info(:)         ! (:) ������� �����
    integer(4), allocatable :: gl_all_Gran_inner(:)   ! ������ �����, ������� ��������� ������ (��������� ��������)
    
    ! ����������� ���������
    integer(4), allocatable :: gl_Contact(:)       ! ������� (�� ������ � ��������� ���������� �����, ������ � ��������� ������� Find_Surface
    
    end module STORAGE
    
    
    
    
    
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
    allocate(gl_Cell_info(size(gl_Cell_A(:,:,:)) + size(gl_Cell_B(:,:,:)) + size(gl_Cell_C(:,:,:))))
    
    allocate( gl_Cell_par(9, size(gl_Cell_A(:,:,:)) + size(gl_Cell_B(:,:,:)) + size(gl_Cell_C(:,:,:)) ) )
    allocate(gl_Cell_par_MF(5, 4, size(gl_Cell_A(:,:,:)) + size(gl_Cell_B(:,:,:)) + size(gl_Cell_C(:,:,:))))
    
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
    allocate(gl_Gran_info(n1))
    allocate(gl_Gran_POTOK(9, n1))
    allocate(gl_Gran_POTOK_MF(5, 4, n1))
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
    gl_Gran_info = 0
    gl_Cell_info = 0
    
    
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

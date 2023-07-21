
    
    module GEO_PARAM                     ! ������ �������������� - �������� ���������� - �������� ���������
    
    implicit none
    
    logical, parameter :: par_developer_info = .True.   ! �������� ������������ ��� ������ �������������� ���������
	logical, parameter :: par_TVD = .True.				! ������ �� ���
	logical, parameter :: par_null_bn = .True.          ! �������� �� bn �� ��������
	
	 real(8), parameter :: par_null_bn_x = -2500.0_8   ! �� ������ ���������� �������� ��������� bn
	
	
    real(8), parameter :: par_R_character = 35.6505         ! ����������� ������ � ������ (���������� �� TS �� ��������� ����� ���������� �����)
    real(8), parameter :: par_koeff_HP = 1.3            ! �� ������� �������� ����������� ������ ��� ��������� HP
    real(8), parameter :: par_koeff_BS = 2.0            ! �� ������� �������� ����������� ������ ��� ��������� BS
    real(8), parameter :: par_R0 = 0.197035         ! ����������� ������ 1 �.�. (���������� �����) ��� ��������� ������ ����� �� ����� �� ����� (������ ��������� � ����)
    ! �������� ����� ��� �������������� �������� ����� eps = par_R_character/10000
    real(8) :: par_R_END = 300.0         !  
    real(8) :: par_R_LEFT = -240.0 ! -390.0         !  ����� �������
    real(8), parameter :: par_pi_8 = acos(-1.0_8)         
    real(8), parameter :: par_pi_4 = acos(-1.0_4)      
	real(8), parameter :: par_sqrtpi = sqrt(par_pi_8)
    real(8), parameter :: cpi4 = 12.56637061435917295384
    real(8), parameter :: ggg = (5.0/3.0)
    real(8), parameter :: par_kk = 10000.0                ! ������� ������������ ������� ������������
	real(8) :: par_R_inner = 9.0! 5.0_8     ! �� ������ ���������� ���������� �����
	
	
	real(8), parameter :: lock_move = 1.0_8 !1.0_8 
	real(8), parameter :: par_nat_TS = lock_move * 0.3 * 0.004_8 ! 0.002_8 !0.0000001_8 !0.003_8                ! ����������� ��������� ������� �����  0.002
	real(8), parameter :: par_nat_HP = lock_move * 0.30_8 ! 0.1  0.8                 ! ����������� ��������� ��������  0.0001
	real(8), parameter :: par_nat_BS = lock_move * 0.00004_8                ! ����������� ��������� ������� ������� ����� 0.0002
	
	real(8), parameter :: koef1 = lock_move * 0.1 * 0.4_8! 0.3      ! ����������� ������������ �������� ������� �����
    real(8), parameter :: koef2 = lock_move * 1.0_8 ! 1.0  0.5  0.01
    real(8), parameter :: koef3 = lock_move * 0.7_8   ! 0.3
	
    
    real(8), parameter :: par_a_2 = 0.130738_8        ! �������� � ������� �����������
    real(8), parameter :: par_n_p_LISM = 3.5_8         ! � �����������
    real(8), parameter :: par_n_H_LISM_ = 1.0_8
    real(8), parameter :: par_Kn = 49.9018   !0.4326569808         ! � �����������
    
    real(8), parameter :: par_chi_real = 1.0_8      ! � ����� �� ������� �������
    real(8), parameter :: par_chi = 41.6479_8      ! � ����� �� ���� ���� �� �������
    real(8), parameter :: par_Velosity_inf = -2.54279_8
	real(8), parameter :: par_Mach_alf = 12.8816_8
	real(8), parameter :: par_Mach_0 = 6.44_8
	real(8), parameter :: par_B_inf = 13.9666_8
	real(8), parameter :: par_alphaB_inf = 1.04719755_8   ! 60 ��������
	real(8), parameter :: par_k_Br = 0.00197035_8
	
	! ��������� ��� �����-�����
	integer(4), parameter :: par_n_potok = 32  ! ����� ������� (� ������� ������ ���� ����)
	integer(4), parameter :: par_n_parallel = 20  ! ��� ����������������� ����� (�.�. ������ ����� ����� � ������� ������������
	! ����� ����� ��������
	integer(4), parameter :: par_n_zone = 6!7  !  ���������� �������� (�� ���� ��� ������� ����)
	integer(4), parameter :: par_m_zone = 7! 6  !  ���������� ����� �� ���� (�� 0 �� 180)
	integer(4), parameter :: par_n_sort = 4  !  ���������� ������ ������
	integer(4), parameter :: par_n_moment = 9  !  ������� ��������� �������� ������� (������ �������)
	real(8), parameter :: par_Rmax = 220.0  !  ������ �����, � ������� ��������� �������
	real(8), parameter :: par_Rleft = -220.0 !-400.0 + 0.01  !  ������ ������
	real(8), parameter :: par_Rup = 250.0 - 0.01  !  ������� ������
	
	! ����� ������ � ������� ������!
	! ����� ������ ���� ������ par_n_parallel
	integer(4), parameter :: MK_N1 = 6 * 600/par_n_parallel   ! ����� �������� ������ ������� ���� (� ���������)
	integer(4), parameter :: MK_N2 = 6 * 100/par_n_parallel  
	integer(4), parameter :: MK_N3 = 0/par_n_parallel     ! (����� �����)
	integer(4), parameter :: MK_N4 = 6 * 100/par_n_parallel   ! (����� ������� � ��������)
	
    
    integer, parameter :: par_R_int = 70  ! ������� �.�. �� ������� ������
    
    real(8) :: par_kk1 = 1.7_8     ! ������� �������� ����� � ���� � ������� �� TS: 1 - ��������, 2 - ������������ � �.�.
	real(8) :: par_kk12 = 1.4_8 ! 1.7_8     ! ������� �������� ����� � ���� � ������� �� TS �� ����� ���������� �����: 1 - ��������, 2 - ������������ � �.�.
    real(8) :: par_kk2 = 2.0_8     ! ������� �������� � �������� ������� �� �������������
    real(8) :: par_kk3 = 1.8_8     ! ������� �������� � ������
	real(8) :: par_kk31 = 1.0_8     ! ������� �������� � ������ ��� ����� �� �������� (������ ����� � � - ����)
	real(8) :: par_kk13 = 0.5_8     ! ������� �������� ����� � �������� ������� �� ������� ������� ����  �� 0 �� 1
	

    integer :: par_l_phi = 80! 48 !4   ! ���������� �������� ��������� ����� ��� ��������� ��������� (���������� �� ���� phi)
                                ! ������ �������� �� 4 ��� �������� ������ ����������� � ����������
	real(8) :: par_al1 = 0.35! 1.0   ! ��� �������� ����� �� ���� (���������������)
    
    integer(4) :: par_m_A = 40! 30      ! ���������� ����� A � ���������
    integer(4) :: par_m_BC = 35! 18      ! ���������� ����� B/C � ���������
    integer(4) :: par_m_O = 20! 17      ! ���������� ����� O � ���������
    integer(4) :: par_m_K = 14! 7      ! ���������� ����� K � ���������
    !integer :: par_m_A = 5      ! ���������� ����� A � ���������
    !integer :: par_m_BC = 5      ! ���������� ����� B/C � ���������
    !integer :: par_m_O = 5      ! ���������� ����� O � ���������
    !integer :: par_m_K = 5      ! ���������� ����� K � ���������
    real(8) :: par_triple_point = 13.0 * par_pi_4/40.0     ! �� ������ ���� ������� �� pi/2 (� �������������� x) ������� �����
    
    ! ���������� ����� �� ����� A
    integer(4) :: par_n_TS =  35! 26                    ! ���������� ����� �� TS (TS ����������)
    integer(4) :: par_n_HP =  65! 40                 ! ���������� ����� HP (HP ����������)  �� �� 0 ���������
    integer(4) :: par_n_BS =  135! 60! 5                 ! ���������� ����� BS (BS ����������)
    integer(4) :: par_n_END = 150! 72! 6                ! ���������� ����� �� ����� ����� (����� ����������)
    integer(4) :: par_n_IA =  17! 12                   ! ���������� �����, ������� ������ �� ���������� �������
	integer(4) :: par_n_IB =  19! 14                   ! ���������� �����, ������� ������ �� ���������� ������� (� �������)
    
    integer :: par_n_points  ! ����� ����� � �����
	
	NAMELIST /SETKA_PARAM/ par_n_TS, par_n_HP, par_n_BS, par_n_END, par_n_IA, par_n_IB, par_triple_point, par_l_phi, &
		par_m_A, par_m_BC, par_m_O, par_m_K, &
		par_R_inner, par_kk1, par_kk12, par_kk2, par_kk3, par_kk31, par_kk13, par_al1, par_n_points
	
	integer, parameter :: par_tet(3, 8) = reshape( (/ 1, 3, 5, 2, 3, 5, 1, 3, 6, 2, 3, 6, 1, 4, 5, 2, 4, 5, 1, 4, 6, 2, 4, 6/), (/3, 8/) )
	
	real(8), parameter :: par_velocity_in(19) = (/ 630.2394085106936, 640.4105750846768, 630.8115381444824, &
		599.5063133133626, 578.7699879954978, 548.7649580013795, 526.4677621421822, &
		491.8843643204692, 452.6566774983212, 427.7441748498429, 432.2370567789446, 475.3568579689629, &
		510.0483064012807, 537.8387913811598, 565.7919325692278, 596.3368522862609, &
		617.5387454952307, 626.1195404367945, 621.216192571027/)  ! �������� �� �� 1 �.�. � �������������� �� ����
	
	real(8), parameter :: par_density_in(19) = (/ 3.273539443111209, 3.396055899930126, 3.822089066410917, &
		4.345206094170982, 4.810128100113988, 5.298016202758408, 5.663744567157964, 5.946408496653754, &
		6.255793360396144, 6.553349257233052, 6.561404691492232, 6.192779851956457, 5.891359370071873, &
		5.660302584234213, 5.303056768788817, 4.642336212635667, 3.865584371689210, 3.394746713363845, &
		3.603241022818391/)  ! / 0.04
    
	end module GEO_PARAM
    
	
    
    module STORAGE                       ! ������ ���������� ������ � ����� (��� ���������� ���������� �� gl - global)
    implicit none

    ! ���� ����� �������� ��������, �� ��� ����� ��� ������� � ������� ��������� ������ � �������������, �� ����� �� ������
    
    real(8), allocatable :: gl_x(:)   ! ����� x-��������� ����� ����� 
    real(8), allocatable :: gl_y(:)   ! ����� y-��������� ����� ����� 
    real(8), allocatable :: gl_z(:)   ! ����� z-��������� ����� ����� 
	
    
    ! �������� MOVE - ��������, ��� ��� ������� ������������ (� ����������������) ������ � ������ �������� �����
    ! �������� �������� �����
    real(8), allocatable :: gl_Vx(:)   !     !MOVE
    real(8), allocatable :: gl_Vy(:)   !     !MOVE
    real(8), allocatable :: gl_Vz(:)   !     !MOVE
    integer(4), allocatable :: gl_Point_num(:)   ! ������� ������ �������� ���� �������� �������� � ������ ����      !MOVE    
    
    ! ���������� ����� (��������� ��� ������ �������� �����)
    real(8), allocatable :: gl_x2(:, :)   ! (:, 2) ����� x-��������� ����� �����     !MOVE
    real(8), allocatable :: gl_y2(:, :)   ! (:, 2) ����� y-��������� ����� �����     !MOVE
    real(8), allocatable :: gl_z2(:, :)   ! (:, 2) ����� z-��������� ����� �����     !MOVE
    
    ! ����, �� ������� ������������� ����� ����� - ��� 
    integer(4), allocatable :: gl_RAY_A(:,:,:)   ! ����� �-����� ����������� 3 (�� ���� ����, � ���� ���������, �� ���� � ������������)
    integer(4), allocatable :: gl_RAY_B(:,:,:)   ! ����� B-����� ����������� 3 (�� ���� ����, � ���� ���������, �� ���� � ������������)
    integer(4), allocatable :: gl_RAY_C(:,:,:)   ! ����� C-����� ����������� 3 (�� ���� ����, � ���� ���������, �� ���� � ������������)
    integer(4), allocatable :: gl_RAY_O(:,:,:)   ! ����� O-����� ����������� 3 (�� ���� ����, � ���� ���������, �� ���� � ������������)
    integer(4), allocatable :: gl_RAY_K(:,:,:)   ! ����� K-����� ����������� 3 (�� ���� ����, � ���� ���������, �� ���� � ������������)
    integer(4), allocatable :: gl_RAY_D(:,:,:)   ! ����� D-����� ����������� 3 (�� ���� ����, � ���� ���������, �� ���� � ������������)
    integer(4), allocatable :: gl_RAY_E(:,:,:)   ! ����� E-����� ����������� 3 (�� ���� ����, � ���� ���������, �� ���� � ������������)
    
    ! ������
    integer(4), allocatable :: gl_Cell_A(:,:,:)   ! ����� A-����� ����������� 3 (�� ���� ����, � ���� ���������, �� ���� � ������������)
    integer(4), allocatable :: gl_Cell_B(:,:,:)   ! ����� B-����� ����������� 3 (�� ���� ����, � ���� ���������, �� ���� � ������������)
    integer(4), allocatable :: gl_Cell_C(:,:,:)   ! ����� C-����� ����������� 3 (�� ���� ����, � ���� ���������, �� ���� � ������������)
    ! ���������� ������ ����� ��� ������� �������������� �����, �������� ������� � �.�.
    ! ������ ��������� �� ������
    
    integer(4), allocatable :: gl_all_Cell(:,:)   ! ���� ����� ����� (8, :) - ������ ���������� ������� - ��� ����� ����� ������
    integer(4), allocatable :: gl_zone_Cell(:)   ! ����� ���� ����������� ������ (1 - 4)
	! 1 - 4
	
    integer(4), allocatable :: gl_all_Cell_inner(:)   ! ������ �����, ������� ��������� ������ (��������� ��������)
    
    integer(4), allocatable :: gl_Cell_neighbour(:,:)   ! (6, :) ����� �� 6 ������� ��� ������ ������ 
	! (���� ����� ������ = 0, �� ��� ��� � ���� �����������)
    ! -1   ! ������� (���������� �����)
	! -3   ! ������� ������� �������
	!  -2  ! �������� �������
	
	integer(4), allocatable :: gl_Cell_gran(:,:)        ! (6, :) ����� �� 6 ������ ��� ������ ������ (���� ����� = 0, �� ����� ��� � ���� �����������)
    integer(4), allocatable :: gl_Cell_info(:)        ! 
    
    real(8), allocatable :: gl_Cell_Volume(:)           ! ����� ������� �����
    real(8), allocatable :: gl_Cell_dist(:)             ! ����������� ���������� �� ����� � ������ ������  4444444444444444
    real(8), allocatable :: gl_Cell_center(:, :)             ! (3, :) ����� ������ ������  4444444444444444
    real(8), allocatable :: gl_Cell_par(:, :)           ! (9, :) ����� ���������� (8 ��������� + Q)
    real(8), allocatable :: gl_Cell_par_MF(:,:,:)           ! ����� ���������� (5, 4,:)  ����������� ��������� (�� 5 ��� ������ �� 4-� ���������)
	real(8), allocatable :: gl_Cell_par_MK(:,:,:)           ! ����� ���������� (9, ������,:)  ����������� ��������� (�� 5 ��� ������ �� 4-� ���������)
    character, allocatable :: gl_Cell_type(:)           ! ��� ������ ������ �, �, �
    integer(4), allocatable :: gl_Cell_number(:, :)     ! (3, :) ����� ������ ������ ������ ������ ����
    
    ! ����� ���������� ��� ��������� ����� (��� ����� ������ ������, ���������� ��� �� ������, ��������� ��� ����������)
    real(8), allocatable :: gl_Cell_Volume2(:, :)           ! (����� �����, 2) ����� ������� �����   !MOVE
    real(8), allocatable :: gl_Cell_center2(:, :, :)             ! (3, :) ����� ������ ������  4444444444444444  !MOVE
    
    
    ! ���� ������ ----------------- ����� ��� ������� ��������, �������, �������, �������� � �.�.
    integer(4), allocatable :: gl_all_Gran(:,:)       ! ��� ����� (4,:) ����� �� 4 ����
    integer(4), allocatable :: gl_Gran_neighbour(:,:) ! ������ ������ ����� (2,:) ����� �� 2 ������, ������� ���� �� ������� �� �������
	! -1  -2  -3  ������
	! -1 - ������� �����
	! -2 - �������� ����� ������
	! -3  -  ������� �������
	
	integer(4), allocatable :: gl_Gran_neighbour_TVD(:,:) ! TVD-������ ������ ����� (2,:) ����� �� 2 ������
	! 0 - ������ ������ ���
	! ��� ���� ������ TVD-����� - ��� ����� ������� �������� ������. �.�. ������� ����� ���� ���� �� ������� �� ������� TVD-������
    
    real(8), allocatable :: gl_Gran_normal(:,:)       ! (3, :) ������� �����   4444444444444444
    real(8), allocatable :: gl_Gran_square(:)         ! (:) ������� �����
    real(8), allocatable :: gl_Gran_POTOK(:,:)       ! (10, :) ����� �����    ��������� - ����������� ���������� ���� ��� �������
    real(8), allocatable :: gl_Gran_POTOK_MF(:,:,:)       ! (5, 4, :) ����� ����� �������������� ���������
    real(8), allocatable :: gl_Gran_center(:,:)       ! (3, :) 
	integer(4), allocatable :: gl_Gran_type(:)      ! ���������� ��� ����� (0 - �������, 1 - TS, 2 - HP)
	integer(4), allocatable :: gl_Gran_scheme(:)      ! ���������� ��� ����� (0 - �������, 1 - TS, 2 - HP)
    integer(4), allocatable :: gl_Gran_info(:)         ! (:) ������������� ����� (��. �����) 
	! ��� ������� ����� ����� ���������� ��� ���������� ����� (������� �������, ���� ������������)
	
    integer(4), allocatable :: gl_all_Gran_inner(:)   ! ������ ������, ������� ��������� ������ (��������� ��������)
    
    ! ����� ���������� ��� ��������� ����� (��� ����� ������ ������, ���������� ��� �� ������, ��������� ��� ����������)
    real(8), allocatable :: gl_Gran_normal2(:, :, :)       ! (3, :, 2) ������� �����                       !MOVE
    real(8), allocatable :: gl_Gran_center2(:, :, :)       ! (3, :, 2)                          !Move
    real(8), allocatable :: gl_Gran_square2(:, :)         ! (:, 2) ������� �����          !MOVE
    
    ! ����������� ���������
    integer(4), allocatable :: gl_Contact(:)       ! ������� (�� ������ � ��������� ���������� �����, ������ � ��������� ������� Find_Surface
    ! ������� ������� �� ������� ������
    integer(4), allocatable :: gl_TS(:)
    integer(4), allocatable :: gl_BS(:)
	
	
	! ����� �������� ��� ��������� ������������ (����������� �������� �����)
	integer(4), allocatable :: gl_Cell_gran_inter(:,:)  ! (6, :)
	real(8), allocatable :: gl_Gran_normal_inter(:,:)   ! (3, :)
	real(8), allocatable :: gl_Gran_center_inter(:,:)  ! (3, :)
	integer(4), allocatable :: gl_Gran_neighbour_inter(:,:)   ! (2, :)
	real(8), allocatable :: gl_Cell_center_inter(:, :)        ! (3, :)
	integer(4), allocatable :: gl_Cell_neighbour_inter(:,:)   ! (6, :)
	real(8), allocatable :: gl_Cell_par_inter(:, :)			  ! (9, :)
    real(8), allocatable :: gl_Cell_par_MF_inter(:,:,:)		  ! (5, 4,:) 
	real(8), allocatable :: gl_Cell_dist_inter(:)
    
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
    
    
    if (allocated(gl_RAY_A) == .True.) then
        STOP "Function Set_STORAGE vizvana neskolko raz!!! Dopustimo tolko 1 raz!"    
    end if
    
    
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
    allocate(gl_Cell_type(size(gl_Cell_A(:,:,:)) + size(gl_Cell_B(:,:,:)) + size(gl_Cell_C(:,:,:))))
    allocate(gl_Cell_number(3, size(gl_Cell_A(:,:,:)) + size(gl_Cell_B(:,:,:)) + size(gl_Cell_C(:,:,:))))
	allocate(gl_zone_Cell(size(gl_Cell_A(:,:,:)) + size(gl_Cell_B(:,:,:)) + size(gl_Cell_C(:,:,:))))
    
    allocate( gl_Cell_par(9, size(gl_Cell_A(:,:,:)) + size(gl_Cell_B(:,:,:)) + size(gl_Cell_C(:,:,:)) ) )
    allocate(gl_Cell_par_MF(5, 4, size(gl_Cell_A(:,:,:)) + size(gl_Cell_B(:,:,:)) + size(gl_Cell_C(:,:,:))))
	allocate(gl_Cell_par_MK(9, par_n_sort, size(gl_Cell_A(:,:,:)) + size(gl_Cell_B(:,:,:)) + size(gl_Cell_C(:,:,:))))
    
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
	allocate(gl_Gran_neighbour_TVD(2, n1))
    allocate(gl_Gran_normal(3, n1))
    allocate(gl_Gran_square(n1))
    allocate(gl_Gran_info(n1))  
	allocate(gl_Gran_type(n1))
	allocate(gl_Gran_scheme(n1))
    allocate(gl_Gran_POTOK(10, n1))
    allocate(gl_Gran_POTOK_MF(5, 4, n1))
    allocate(gl_Gran_center(3, n1))
    
    allocate(gl_Contact( (par_m_O + par_m_A + par_m_BC -1) * par_l_phi ))   ! �������� ������ ��� �������
    allocate(gl_TS( (par_m_A + par_m_BC + par_m_K -1) * par_l_phi ))   ! �������� ������ ��� TS
    allocate(gl_BS( (par_m_A - 1) * par_l_phi ))   ! �������� ������ ��� BS
    
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
    
	gl_Cell_par_MK = 0.0
    gl_Cell_dist = 0.0
    gl_Cell_center = 0.0
    gl_Gran_POTOK = 0.0
    gl_Cell_par = 0.0
    gl_Gran_center = 0.0
    gl_Gran_info = 0
    gl_Cell_info = 2
	gl_Gran_type = 0
	gl_Gran_scheme = 2
	gl_zone_Cell = 2
    
    
    gl_Cell_A = -1
    gl_Cell_B = -1
    gl_Cell_C = -1
    
    gl_all_Cell = -1
    gl_Cell_neighbour = 0          ! 0 - ��� ������ (��� ���� ��� �� ������� - � ������� ���� �����)
    gl_Cell_gran = 0
    
    gl_all_Gran = 0
    gl_Gran_neighbour = 0
	gl_Gran_neighbour_TVD = 0
    
    gl_Contact = 0
    
    gl_Gran_normal = 0.0;
    gl_Gran_square = 0.0;
    
    if (par_developer_info) print *, "END Set_STORAGE"
    
    end subroutine Set_STORAGE

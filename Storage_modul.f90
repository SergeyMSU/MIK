
    
module GEO_PARAM                     ! ������ �������������� - �������� ���������� - �������� ���������
    
    implicit none
    
    logical, parameter :: par_developer_info = .True.   ! �������� ������������ ��� ������ �������������� ���������
	logical, parameter :: par_TVD = .True.				! ������ �� ���
	logical, parameter :: par_null_bn = .True.          ! �������� �� bn �� ��������
	logical, parameter :: par_helium = .True.          ! �������� �� �����? (��� ���� ��������� ������� � �.�.)
	logical, parameter :: par_PUI = .True.          ! �������� �� �����? (��� ���� ��������� ������� � �.�.)
	
	 real(8), parameter :: par_null_bn_x = -10000.0_8 !-2500.0_8   ! �� ������ ���������� �������� ��������� bn
	
	
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
    real(8), parameter :: cpi4 = 12.56637061435917295384_8
    real(8), parameter :: ggg = (5.0/3.0)
    real(8), parameter :: par_kk = 10000.0                ! ������� ������������ ������� ������������
	real(8) :: par_R_inner = 9.0! 5.0_8     ! �� ������ ���������� ���������� �����
	
	
	real(8), parameter :: lock_move = 1.0_8 !1.0_8 
	real(8), parameter :: par_nat_TS = 0.3 * lock_move * 0.3 * 0.01_8 ! 0.7 * 0.002_8 !0.0000001_8 !0.003_8                ! ����������� ��������� ������� �����  0.002
	real(8), parameter :: par_nat_HP = 0.8 * 0.1 * lock_move! 0.1  0.8                 ! ����������� ��������� ��������  0.0001
	real(8), parameter :: par_nat_BS = lock_move * 0.00004_8                ! ����������� ��������� ������� ������� ����� 0.0002
	
	real(8), parameter :: koef1 = lock_move * 0.1! 0.2  � 5 ��� ��������     ! ����������� ������������ �������� ������� �����
    real(8), parameter :: koef2 = 0.1 * lock_move ! 1.0  0.5  0.01
    real(8), parameter :: koef3 = lock_move * 0.7_8   ! 0.3
	
    
    real(8), parameter :: par_a_2 = 0.130738_8        ! �������� � ������� �����������
    real(8), parameter :: par_n_p_LISM = 3.5_8         ! � �����������
    real(8), parameter :: par_n_H_LISM_ = 1.0_8
    real(8), parameter :: par_Kn = 49.9018   !0.4326569808         ! � �����������
    real(8), parameter :: par_nu_ph = 12.2125 
    real(8), parameter :: par_E_ph = 0.10878
	
    
    real(8), parameter :: par_chi_real = 41.6479_8! 1.0_8      ! � ����� �� ������� �������
    real(8), parameter :: par_chi = 41.6479_8      ! � ����� �� ���� ���� �� �������
    real(8), parameter :: par_Velosity_inf = -2.54279_8
	real(8), parameter :: par_Mach_alf = 12.8816_8
	real(8), parameter :: par_Mach_0 = 6.44_8
	real(8), parameter :: par_p_0 = 4624.57_8 ! 4790.19! 4624.57_8
	real(8), parameter :: par_B_inf = 13.9666_8
	real(8), parameter :: par_alphaB_inf = 1.04719755_8   ! 60 ��������
	real(8), parameter :: par_k_Br = 0.00197035_8
	real(8), parameter :: par_1ae = 0.197035_8
	
	! ��������� ��� �����-�����
	integer(4), parameter :: par_n_potok = 32! 32  ! ����� ������� (� ������� ������ ���� ����)
	integer(4), parameter :: par_n_parallel = 20  ! ��� ����������������� ����� (�.�. ������ ����� ����� � ������� ������������
	integer(4), parameter :: par_n_claster = 1  ! ����� ����������� (��� MPI)
	! ����� ����� ��������
	integer(4), parameter :: par_n_zone = 6!7  !  ���������� �������� (�� ���� ��� ������� ����)
	integer(4), parameter :: par_m_zone = 7! 6  !  ���������� ����� �� ���� (�� 0 �� 180)
	integer(4), parameter :: par_n_sort = 4  !  ���������� ������ ������
	integer(4) :: par_n_moment = 19 !9  !  ������� ��������� �������� ������� (������ �������)
	real(8), parameter :: par_Rmax = 220.0  !  ������ �����, � ������� ��������� �������
	real(8) :: par_Rleft! = par_R_LEFT + 0.0001 !-400.0 + 0.01  !  ������ ������
	real(8) :: par_Rup! = par_R_END - 0.001  !  ������� ������
	
	! ����� ������ � ������� ������!
	! ����� ������ ���� ������ par_n_parallel
	integer(4), parameter :: MK_k_multiply = 3 * 18!6 * 11! 17   ! 1 = 10 ����� ����� (� �������� 18 �����)
	integer(4), parameter :: MK_k_mul1 = 6 * MK_k_multiply! 6
	integer(4), parameter :: MK_k_mul2 = 1 * MK_k_multiply! 
	integer(4), parameter :: MK_N1 = MK_k_mul1 * 60/par_n_parallel   ! 600 ����� �������� ������ ������� ���� (� ���������)
	integer(4), parameter :: MK_N2 = MK_k_mul1 * 20/par_n_parallel  ! 200
	integer(4), parameter :: MK_N3 = MK_k_mul2 * 20/par_n_parallel     ! (����� �����)
	integer(4), parameter :: MK_N4 = MK_k_mul1 * 20/par_n_parallel   ! (����� ������� � ��������)
	
    
    integer, parameter :: par_R_int = 70  ! ������� �.�. �� ������� ������
    
    real(8) :: par_kk1 = 1.7_8     ! ������� �������� ����� � ���� � ������� �� TS: 1 - ��������, 2 - ������������ � �.�.
	!real(8) :: par_kk12 = 1.4_8 ! 1.7_8     ! ������� �������� ����� � ���� � ������� �� TS �� ����� ���������� �����: 1 - ��������, 2 - ������������ � �.�.
    real(8) :: par_kk2 = 2.0_8     ! ������� �������� � �������� ������� �� �������������
    real(8) :: par_kk3 = 1.8_8     ! ������� �������� � ������
	real(8) :: par_kk31 = 1.0_8     ! ������� �������� � ������ ��� ����� �� �������� (������ ����� � � - ����)
	real(8) :: par_kk13 = 1.8_8     ! ������� �������� ����� � �������� ������� �� ������� ������� ����  �� 0 �� 1
	                                ! �������� ����? 
	real(8) :: par_kk131 = 0.1_8
	real(8) :: par_kk132 = 1.5_8
	
	real(8) :: par_kk14 = 1.0_8     ! ������� �������� ����� � �������� ������� �� ���������� ������� ����  �� 0 �� 1
	!                                 (�������� ����� � TS � HP)  
	real(8) :: par_kk12 = 1.4_8     ! ������� �������� ����� �� TS � ������� �����  >= 1

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
    real(8) :: par_triple_point = 13.0 * par_pi_8/40.0     ! �� ������ ���� ������� �� pi/2 (� �������������� x) ������� �����
    real(8) :: par_triple_point_2 = 7.0 * par_pi_8/40.0     ! ��� ����� ����� ������� ��� ����� ������� ����� ������� �� pi/2 (� �������������� x) 
    
    ! ���������� ����� �� ����� A
    integer(4) :: par_n_TS =  35! 26                    ! ���������� ����� �� TS (TS ����������)
    integer(4) :: par_n_HP =  65! 40                 ! ���������� ����� HP (HP ����������)  �� �� 0 ���������
    integer(4) :: par_n_BS =  135! 60! 5                 ! ���������� ����� BS (BS ����������)
    integer(4) :: par_n_END = 150! 72! 6                ! ���������� ����� �� ����� ����� (����� ����������)
    integer(4) :: par_n_IA =  17! 12                   ! ���������� �����, ������� ������ �� ���������� �������
	integer(4) :: par_n_IB =  19! 14                   ! ���������� �����, ������� ������ �� ���������� ������� (� �������)
    
    integer :: par_n_points  ! ����� ����� � �����
	
	NAMELIST /SETKA_PARAM/ par_n_TS, par_n_HP, par_n_BS, par_n_END, par_n_IA, par_n_IB, par_triple_point, par_triple_point_2, par_l_phi, &
		par_m_A, par_m_BC, par_m_O, par_m_K, &
		par_R_inner, par_kk1, par_kk12, par_kk2, par_kk3, par_kk31, par_kk13, par_al1, par_n_points
	
	integer, parameter :: par_tet(3, 8) = reshape( (/ 1, 3, 5, 2, 3, 5, 1, 3, 6, 2, 3, 6, 1, 4, 5, 2, 4, 5, 1, 4, 6, 2, 4, 6/), (/3, 8/) )
	
	real(8), parameter :: par_velocity_in(19) = (/ 630.728222, 640.919966, 631.236137, &
		599.740090, 578.917817, 548.910412, 526.687826, &
		492.094696, 452.743990, 427.788554, 432.315053, 475.512527, &
		510.250092, 538.045884, 565.994056, 596.540539, &
		617.781334, 626.338464, 621.444686/)  ! �������� �� �� 1 �.�. � ��c����������� �� ����
	
	real(8), parameter :: par_density_in(19) = (/ 3.264516, 3.387225, 3.812428, &
		4.334442, 4.798314, 5.285154, 5.649792, 5.931632, &
		6.240762, 6.537483, 6.545379, 6.177649, 5.877042, &
		5.646527, 5.289991, 4.630639, 3.855446, 3.386001, &
		3.594525/)  ! / 0.04
    
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
    integer(4), allocatable :: gl_zone_Cell(:)   ! ����� ���� ����������� ������ (1 - 4)   INFO
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
    real(8), allocatable :: gl_Cell_par2(:, :)           ! (1?, :) (n_He, div(V))
    real(8), allocatable :: gl_Cell_par_MF(:,:,:)           ! ����� ���������� (5, 4,:)  ����������� ��������� (�� 5 ��� ������ �� 4-� ���������)
	real(8), allocatable :: gl_Cell_par_MK(:,:,:)           ! ����� ���������� (10, ������,:)  ����������� ��������� (�� 5 ��� ������ �� 4-� ���������)
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
    real(8), allocatable :: gl_Gran_POTOK2(:,:)       ! (10, :) ����� �����    ��������� - ����������� ���������� ���� ��� �������
    real(8), allocatable :: gl_Gran_POTOK_MF(:,:,:)       ! (5, 4, :) ����� ����� �������������� ���������
    real(8), allocatable :: gl_Gran_center(:,:)       ! (3, :) 
	integer(4), allocatable :: gl_Gran_type(:)      ! ���������� ��� ����� (0 - �������, 1 - TS, 2 - HP)      INFO
	integer(4), allocatable :: gl_Gran_scheme(:)      ! ���������� ��� ����� (0 - �������, 1 - TS, 2 - HP)    INFO
    integer(4), allocatable :: gl_Gran_info(:)         ! (:) ������������� ����� (��. ����� ����)                  INFO
	! 0, 1, 2   
	! 2 - ��� ����� ����� �������� ������ ����-����� (��� ��������� ����, 50 �.�. ��������)
	! 1 - ��� ����������� �����
	! 0 - ��� ����� ����� �������� �������� �� �������� ��������� �������
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


    !! ����� �������� ��� ��������� PUI 
    ! �� ����������� � ���������� �������� �������� � ������ PUI, ����� ��� �����
    ! ������ ��� ����, ����� �� ����� ���� ������������ � ������ �������

    integer :: pui_nW = 40      ! 50
	real(8) :: pui_wR = 100.0    ! 150.0
    real(8), allocatable :: f_pui(:, :)           ! (pui_nW, :)
	integer, allocatable :: f_pui_num(:)           ! �� ������ � ������� ���, ���������� ����� ���� � ���������������� �����
	integer, allocatable :: f_pui_num2(:)		   ! �� ������ ���� � ���������������� �����, ���������� ����� � ������� PUI
    real(8), allocatable :: n_pui(:)           ! (:)	   ������������ �������
	real(8), allocatable :: T_pui(:)           ! (:)	   ����������� �������
    
    contains
    
    
    subroutine Set_STORAGE()             ! ������� �������� ������ ��� ���� ��������� ������ STORAGE, ��������� ��������� �� ������ GEO_PARAM
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
	
    if(par_helium) allocate( gl_Cell_par2(1, size(gl_Cell_A(:,:,:)) + size(gl_Cell_B(:,:,:)) + size(gl_Cell_C(:,:,:)) ) )
	
    allocate(gl_Cell_par_MF(5, 4, size(gl_Cell_A(:,:,:)) + size(gl_Cell_B(:,:,:)) + size(gl_Cell_C(:,:,:))))
	allocate(gl_Cell_par_MK(10, par_n_sort, size(gl_Cell_A(:,:,:)) + size(gl_Cell_B(:,:,:)) + size(gl_Cell_C(:,:,:))))
    
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
    allocate(gl_Gran_POTOK2(size(gl_Cell_par2(:, 1)), n1))
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
    if (par_developer_info) print*, "Kol-vo ploskostey po ygly =  ", size(gl_RAY_A(1, 1, :))
    
    
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
	gl_Cell_par_MK(1, :, :) = 0.000001
	
    gl_Cell_dist = 0.0
    gl_Cell_center = 0.0
    gl_Gran_POTOK = 0.0
	gl_Gran_POTOK2 = 0.0
    gl_Cell_par = 0.0
    gl_Gran_center = 0.0
    gl_Gran_info = 0
    gl_Cell_info = 2
	gl_Gran_type = 0
	!gl_Gran_scheme = 2
	gl_Gran_scheme = 3
	gl_zone_Cell = 2
	if(par_helium) gl_Cell_par2 = 0.0
    
    
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

	subroutine Dell_STORAGE()
	use GEO_PARAM
	deallocate(gl_RAY_A)
    deallocate(gl_RAY_B)
    deallocate(gl_RAY_C)
    deallocate(gl_RAY_O)
    deallocate(gl_RAY_K)
    deallocate(gl_RAY_D)
    deallocate(gl_RAY_E)
    
    deallocate(gl_Cell_A)
    deallocate(gl_Cell_B)
    deallocate(gl_Cell_C)
    
    
    deallocate(gl_all_Cell)
    deallocate(gl_Cell_neighbour)
    deallocate(gl_Cell_gran)
    
    deallocate(gl_Cell_Volume)
    deallocate(gl_Cell_dist)
    deallocate( gl_Cell_center)
    deallocate(gl_Cell_info)
    deallocate(gl_Cell_type)
    deallocate(gl_Cell_number)
	deallocate(gl_zone_Cell)
    
    deallocate(gl_Cell_par)
	if(par_helium) deallocate(gl_Cell_par2)
    deallocate(gl_Cell_par_MF)
	deallocate(gl_Cell_par_MK)
    
    deallocate(gl_x)
    deallocate(gl_y)
    deallocate(gl_z)
    
    deallocate(gl_all_Gran)
    deallocate(gl_Gran_neighbour)
	deallocate(gl_Gran_neighbour_TVD)
    deallocate(gl_Gran_normal)
    deallocate(gl_Gran_square)
    deallocate(gl_Gran_info)
	deallocate(gl_Gran_type)
	deallocate(gl_Gran_scheme)
    deallocate(gl_Gran_POTOK)
	deallocate(gl_Gran_POTOK2)
    deallocate(gl_Gran_POTOK_MF)
    deallocate(gl_Gran_center)
	
    
    deallocate(gl_Contact)
    deallocate(gl_TS)
    deallocate(gl_BS)
	
	end subroutine Dell_STORAGE
	
end module STORAGE
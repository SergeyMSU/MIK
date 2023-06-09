! Для запуска на Линукс нужно изменить расширение файла на cuf

	
module MY_CUDA
	 use cudafor
	 use Solvers
	 use My_func
	 
	 ! Константы
	 real(8), constant :: dev_ggg
	 real(8), constant :: dev_par_Velosity_inf
	 
	 integer(4), constant :: dev_Ngran
	 integer(4), constant :: dev_Ncell
	 real(8), constant :: dev_par_R0
	 integer(4), constant :: dev_Ngran_inner
	 integer(4), constant :: dev_Ncell_inner
	
	 real(8), constant :: dev_par_pi_8  
	 real(8), constant :: dev_par_a_2
	 real(8), constant :: dev_par_n_p_LISM
    real(8), constant :: dev_par_n_H_LISM_
    real(8), constant :: dev_par_Kn
	real(8), constant :: dev_par_R_inner
	
	integer(4), constant :: dev_par_n_TS
	integer(4), constant :: dev_par_n_HP
	integer(4), constant :: dev_par_n_BS
	integer(4), constant ::	dev_par_n_END
	real(8), constant ::	dev_par_triple_point
	integer(4), constant ::	dev_par_m_BC
	integer(4), constant :: dev_par_n_IA
	integer(4), constant :: dev_par_n_IB
	real(8), constant ::  dev_par_kk1     ! Степень сгущения сетки к нулю в области до TS: 1 - линейное, 2 - квадратичное и т.д.
	real(8), constant ::  dev_par_kk12     ! Степень сгущения сетки к нулю в области до TS но после внутренней сферы: 1 - линейное, 2 - квадратичное и т.д.
    real(8), constant ::  dev_par_kk2     ! Степень сгущения в головной области на бесконечности
    real(8), constant ::  dev_par_kk3     ! Степень сгущения в хвосте
	real(8), constant ::  dev_par_kk31
	real(8), constant ::  dev_par_R_END
	real(8), constant ::  dev_par_R_LEFT
	real(8), constant ::  dev_par_al1
	
	
	integer(4), device :: dev_mutex_1
	integer(4), device :: dev_mutex_2
	integer(4), device :: dev_mutex_3
	integer(4), device :: dev_mutex_4
	 
	! Создаём набор необходимых массивов для неподвижной сетки, аналогичных массивам на хосте
	 integer(4), device, allocatable :: dev_gl_Gran_neighbour(:,:)
	 integer(4), device, allocatable :: dev_gl_Gran_info(:)         ! (:) классификатор грани (см. схему)
	 integer(4), device, allocatable :: dev_gl_Gran_scheme(:)
	 integer(4), device, allocatable :: dev_gl_Gran_type(:)         ! (:) классификатор грани (см. схему)
	 integer(4), device, allocatable :: dev_gl_Cell_info(:)
	 real(8), device, allocatable :: dev_gl_Cell_par(:, :)
	 real(8), device, allocatable :: dev_gl_Cell_par_MF(:,:,:)
	 real(8), device, allocatable :: dev_gl_Cell_center(:, :)
	 real(8), device, allocatable :: dev_gl_Gran_normal(:,:)       ! (3, :) Нормаль грани   4444444444444444
     real(8), device, allocatable :: dev_gl_Gran_square(:)         ! (:) Площадь грани
	 real(8), device, allocatable :: dev_gl_Gran_center(:,:)
	 real(8), device, allocatable :: dev_gl_Cell_dist(:)
	 real(8), device, allocatable :: dev_gl_Cell_Volume(:)
	 real(8), device, allocatable :: dev_gl_Gran_POTOK(:,:)       ! (8, :) поток грани
	 real(8), device, allocatable :: dev_gl_Gran_POTOK_MF(:,:,:)       ! (5, 4, :) поток грани мультифлюидных жидкостей
	 integer(4), device, allocatable :: dev_gl_Cell_gran(:,:)
	 integer(4), device, allocatable :: dev_gl_all_Gran_inner(:)
	 integer(4), device, allocatable :: dev_gl_all_Cell_inner(:)
	 character,  device, allocatable :: dev_gl_Cell_type(:)
	 integer(4),  device, allocatable :: dev_gl_Cell_number(:, :)
	 integer(4), device, allocatable :: dev_gl_Gran_neighbour_TVD(:,:)
	 real(8), device, allocatable :: dev_gl_x(:)   ! набор z-координат узлов сетки    !MOVE
     real(8), device, allocatable :: dev_gl_y(:)   ! набор z-координат узлов сетки    !MOVE
     real(8), device, allocatable :: dev_gl_z(:)   ! набор z-координат узлов сетки    !MOVE
	 
	 ! Создаём набор массивов для подвижной сетки, аналогичных массивам на хосте
	 
	 real(8), device, allocatable :: dev_gl_Vx(:)   ! набор z-координат узлов сетки    !MOVE
     real(8), device, allocatable :: dev_gl_Vy(:)   ! набор z-координат узлов сетки    !MOVE
     real(8), device, allocatable :: dev_gl_Vz(:)   ! набор z-координат узлов сетки    !MOVE
     integer(4), device, allocatable :: dev_gl_Point_num(:)   ! Сколько граней записали свою скорость движения в данный узел      !MOVE    
     real(8), device, allocatable :: dev_gl_x2(:, :)   ! (:, 2) набор x-координат узлов сетки     !MOVE
     real(8), device, allocatable :: dev_gl_y2(:, :)   ! (:, 2) набор y-координат узлов сетки     !MOVE
     real(8), device, allocatable :: dev_gl_z2(:, :)   ! (:, 2) набор z-координат узлов сетки     !MOVE
	 integer(4), device, allocatable :: dev_gl_all_Cell(:,:)   ! Весь набор ячеек (8, :) - первая координата массива - это набор узлов ячейки
	 real(8), device, allocatable :: dev_gl_Cell_Volume2(:, :)           ! (всего ячеек, 2) Набор объёмов ячеек   !MOVE
     real(8), device, allocatable :: dev_gl_Cell_center2(:, :, :)             ! (3, :) Центр каждой ячейки  4444444444444444  !MOVE
	 real(8), device, allocatable :: dev_gl_Gran_normal2(:, :, :)       ! (3, :, 2) Нормаль грани                       !MOVE
     real(8), device, allocatable :: dev_gl_Gran_center2(:, :, :)       ! (3, :, 2)                          !Move
     real(8), device, allocatable :: dev_gl_Gran_square2(:, :)         ! (:, 2) Площадь грани          !MOVE
	 integer(4), device, allocatable :: dev_gl_Contact(:)       ! Контакт (не входит в начальное построение сетки, ищется в отдельной функции Find_Surface
     integer(4), device, allocatable :: dev_gl_TS(:)
     integer(4), device, allocatable :: dev_gl_BS(:)
	 real(8), device, allocatable :: dev_gl_Vel_gran(:)        ! Скорость граней (для движения граней - которые выделяются) 
	 ! нужно было для избежания атомарных операции и безпроблемного доступа к памяти
	 integer(4), device, allocatable :: dev_gl_zone_Cell(:)   ! Какой зоне принадлежит ячейка
	 
	 integer(4), device, allocatable :: dev_gl_all_Gran(:,:)       ! Все грани (4,:) имеют по 4 узла
	 
	 ! Лучи, на которых распологаются точки сетки - так 
    integer(4), device, allocatable :: dev_gl_RAY_A(:,:,:)   ! Набор А-лучей размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
    integer(4), device, allocatable :: dev_gl_RAY_B(:,:,:)   ! Набор B-лучей размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
    integer(4), device, allocatable :: dev_gl_RAY_C(:,:,:)   ! Набор C-лучей размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
    integer(4), device, allocatable :: dev_gl_RAY_O(:,:,:)   ! Набор O-лучей размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
    integer(4), device, allocatable :: dev_gl_RAY_K(:,:,:)   ! Набор K-лучей размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
    integer(4), device, allocatable :: dev_gl_RAY_D(:,:,:)   ! Набор D-лучей размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
    integer(4), device, allocatable :: dev_gl_RAY_E(:,:,:)   ! Набор E-лучей размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
	 
	 
	 real(8), device :: time_all
	 real(8), device :: time_step
	 real(8), device :: time_step2
	 real(8), device :: dev_time_step_inner
	 
	 interface
	 
	 
	 attributes(device) subroutine spherical_skorost(z, x, y, Vz, Vx, Vy, Vr, Vphi, Vtheta)
    use My_func
    implicit none
    real(8), intent(in) :: x, y, z, Vx, Vy, Vz
    real(8), intent(out) :: Vr, Vphi, Vtheta
	
	 end subroutine spherical_skorost
	 
	 attributes(device) subroutine dekard_skorost(z, x, y, Vr, Vphi, Vtheta, Vz, Vx, Vy)
    use My_func
    implicit none
    real(8), intent(in) :: x, y, z,  Vr, Vphi, Vtheta
    real(8), intent(out) :: Vx, Vy, Vz
	end subroutine dekard_skorost
	
	!attributes(device) subroutine chlld_Q(n_state, al, be, ge, &
 !                                w, qqq1, qqq2, &
 !                                dsl, dsp, dsc, &
 !                                qqq, null_bn1, n_disc)
 !     implicit real*8 (a-h,o-z)
 !     real(8), intent(out) :: dsl, dsp, dsc
 !     real(8), intent(in) :: al, be, ge, w
 !     integer(4), intent(in) :: n_state
	!  logical, intent(in), optional :: null_bn1
	!  integer, intent(in), optional :: n_disc
 !     dimension qqq(9),qqq1(9),qqq2(9)
	! end subroutine chlld_Q
	 
	 attributes(device) subroutine chlld_gd(n_state, al, be, ge, &  
                                 w, qqq1, qqq2, &
                                 dsl, dsp, dsc, &
                                 qqq)
      implicit real*8 (a-h,o-z)
      
      real(8), intent(out) :: dsl, dsp, dsc
      real(8), intent(in) :: al, be, ge, w
      integer(4), intent(in) :: n_state
      
      dimension qqq(5),qqq1(5),qqq2(5)
	end subroutine chlld_gd
	
	!attributes(device) subroutine chlld(n_state, al, be, ge, &
 !                                w, qqq1, qqq2, &
 !                                dsl, dsp, dsc, &
 !                                qqq, n_disc)
	!
 !     implicit real*8 (a-h,o-z)
 !     
 !     real(8), intent(out) :: dsl, dsp, dsc
 !     real(8), intent(in) :: al, be, ge, w
 !     integer(4), intent(in) :: n_state
 !     dimension qqq(8),qqq1(8),qqq2(8)
	!  integer, intent(in), optional :: n_disc
	!  
	!  end subroutine chlld
	
    end interface
	
	contains 
	
	subroutine Set_CUDA()
		use GEO_PARAM
		use STORAGE
		integer(4) :: Ncell, Ngran
	
		 if (allocated(dev_gl_Gran_neighbour) == .True.) then
		 	STOP "Function Set_CUDA vizvana neskolko raz!!! Dopustimo tolko 1 raz!"    
		 end if
		
		 if (allocated(gl_Cell_par) == .False.) then
		 	STOP "Function Set_storage ne vizvana, Set_CUDA precrashaet rabotu"    
		 end if
	
		 Ncell = size(gl_Cell_par(1, :))
		 Ngran = size(gl_Gran_neighbour(1, :))
		
		 allocate(dev_gl_x(size(gl_x(:))))
		 allocate(dev_gl_y(size(gl_x(:))))
		 allocate(dev_gl_z(size(gl_x(:))))
		 
		 allocate(dev_gl_Gran_neighbour(2, Ngran))
		 allocate(dev_gl_Gran_info(Ngran))
		 allocate(dev_gl_Gran_scheme(Ngran))
		 allocate(dev_gl_Gran_type(Ngran))
		 allocate(dev_gl_Cell_par(9, Ncell))
		 allocate(dev_gl_Cell_par_MF(5, 4, Ncell))
		 allocate(dev_gl_Cell_center(3, Ncell))
		 allocate(dev_gl_Gran_normal(3, Ngran))
		 allocate(dev_gl_Gran_square(Ngran))
		 allocate(dev_gl_Gran_center(3, Ngran))
		 allocate(dev_gl_Cell_dist(Ncell))
		 allocate(dev_gl_Cell_Volume(Ncell))
		 allocate(dev_gl_Gran_POTOK( size(gl_Gran_POTOK(:, 1))  , Ngran))
		 allocate(dev_gl_Gran_POTOK_MF(5, 4, Ngran))
		 allocate(dev_gl_Cell_info(Ncell))
		 allocate(dev_gl_Cell_gran(6, Ncell))
		 allocate(dev_gl_all_Gran_inner( size(gl_all_Gran_inner(:) )))
		 allocate(dev_gl_all_Cell_inner( size(gl_all_Cell_inner(:) )))
		 allocate(dev_gl_Cell_type( size(gl_Cell_type(:))  ))
		 allocate(dev_gl_Cell_number(3,  size(gl_Cell_number(1, :))  ))
		 
		 allocate(dev_gl_Gran_neighbour_TVD, mold = gl_Gran_neighbour_TVD)
		 
		 dev_mutex_1 = 0
		 dev_mutex_2 = 0
		 dev_mutex_3 = 0
		 dev_mutex_4 = 0
		 time_all = 0.0
		 time_step = 1.0 ! 0.00000000001
		 
		 dev_gl_Gran_POTOK = 0.0
		 dev_gl_Gran_POTOK_MF = 0.0
	
	end subroutine Set_CUDA
	
	subroutine Alloc_CUDA_move()
		use GEO_PARAM
		use STORAGE
		implicit none
		integer(4) :: Ncell, Ngran, Npoint
	
		 if (allocated(dev_gl_Vx) == .True.) then
		 	STOP "Function Set_CUDA_move vizvana neskolko raz!!! Dopustimo tolko 1 raz!"    
		 end if
		
		 if (allocated(gl_Cell_par) == .False.) then
		 	STOP "Function Set_storage ne vizvana, Set_CUDA_move precrashaet rabotu"    
		 end if
	
		 Ncell = size(gl_Cell_par(1, :))
		 Ngran = size(gl_Gran_neighbour(1, :))
		 Npoint = size(gl_x(:))
		
		 
		 allocate(dev_gl_Vx(Npoint))
		 allocate(dev_gl_Vy(Npoint))
		 allocate(dev_gl_Vz(Npoint))
		 
		 allocate(dev_gl_Point_num(Npoint))
		 allocate(dev_gl_x2(Npoint, 2))
		 allocate(dev_gl_y2(Npoint, 2))
		 allocate(dev_gl_z2(Npoint, 2))
		 allocate(dev_gl_all_Cell(8, size(gl_all_Cell(1, :))  ))
		 
		 allocate(dev_gl_Cell_Volume2(Ncell, 2))
		 allocate(dev_gl_Cell_center2(3, Ncell, 2))
		 allocate(dev_gl_Gran_normal2(3, Ngran, 2))
		 allocate(dev_gl_Gran_center2(3, Ngran, 2))
		 allocate(dev_gl_Gran_square2(Ngran, 2))
		 allocate(dev_gl_Vel_gran(Ngran))
		 allocate(dev_gl_zone_Cell(Ncell))
		 
		 allocate(dev_gl_Contact(size(gl_Contact(:))))
		 allocate(dev_gl_TS(size(gl_TS(:))))
		 allocate(dev_gl_BS(size(gl_BS(:))))
		 allocate(dev_gl_all_Gran(4, Ngran))
		 
		 
		 allocate(dev_gl_RAY_A, MOLD = gl_RAY_A)
		 allocate(dev_gl_RAY_B, MOLD = gl_RAY_B)
		 allocate(dev_gl_RAY_C, MOLD = gl_RAY_C)
		 allocate(dev_gl_RAY_O, MOLD = gl_RAY_O)
		 allocate(dev_gl_RAY_K, MOLD = gl_RAY_K)
		 allocate(dev_gl_RAY_D, MOLD = gl_RAY_D)
		 allocate(dev_gl_RAY_E, MOLD = gl_RAY_E)
		 
		 dev_gl_Gran_square2 = 0.0
		 dev_gl_Gran_center2 = 0.0
		 dev_gl_Gran_normal2 = 0.0
		 dev_gl_Cell_center2 = 0.0
		 dev_gl_Cell_Volume2 = 0.0
		 
		 
	end subroutine Alloc_CUDA_move
	
	subroutine Set_CUDA_move()
		use STORAGE
		implicit none
	
		dev_gl_Vx = 0.0
        dev_gl_Vy = 0.0
        dev_gl_Vz = 0.0
        dev_gl_Point_num = 0
        dev_gl_x2(:, 1) = dev_gl_x
        dev_gl_x2(:, 2) = dev_gl_x
        dev_gl_y2(:, 1) = dev_gl_y
        dev_gl_y2(:, 2) = dev_gl_y
        dev_gl_z2(:, 1) = dev_gl_z
        dev_gl_z2(:, 2) = dev_gl_z
		dev_gl_all_Cell = gl_all_Cell
        dev_gl_Cell_Volume2(:, 1) = dev_gl_Cell_Volume
        dev_gl_Cell_Volume2(:, 2) = dev_gl_Cell_Volume
        dev_gl_Gran_normal2(:, :, 1) = dev_gl_Gran_normal
        dev_gl_Gran_normal2(:, :, 2) = dev_gl_Gran_normal
        dev_gl_Gran_center2(:, :, 1) = dev_gl_Gran_center
        dev_gl_Gran_center2(:, :, 2) = dev_gl_Gran_center
        dev_gl_Cell_center2(:, :, 1) = dev_gl_Cell_center
        dev_gl_Cell_center2(:, :, 2) = dev_gl_Cell_center
        dev_gl_Gran_square2(:, 1) = dev_gl_Gran_square
        dev_gl_Gran_square2(:, 2) = dev_gl_Gran_square
		dev_gl_Contact = gl_Contact
		dev_gl_TS = gl_TS
		dev_gl_BS = gl_BS
		dev_gl_all_Gran = gl_all_Gran
		dev_gl_RAY_A  = gl_RAY_A
		dev_gl_RAY_B  = gl_RAY_B
		dev_gl_RAY_C = gl_RAY_C
		dev_gl_RAY_O = gl_RAY_O
		dev_gl_RAY_K = gl_RAY_K
		dev_gl_RAY_D = gl_RAY_D
		dev_gl_RAY_E = gl_RAY_E
		dev_gl_Vel_gran = 0.0
		dev_gl_zone_Cell = gl_zone_Cell
		 
	end subroutine Set_CUDA_move
	
	subroutine Set_CUDA_move_reverse(now2)
	! Для чего эта функция????????????????????????????????????????????????????????????????????????
		implicit none
		integer(4), intent(in) :: now2
		
		dev_gl_x = dev_gl_x2(:, now2)
        dev_gl_y = dev_gl_y2(:, now2)
        dev_gl_z = dev_gl_z2(:, now2)
        dev_gl_Cell_Volume = dev_gl_Cell_Volume2(:, now2)
        dev_gl_Gran_normal = dev_gl_Gran_normal2(:, :, now2)
        dev_gl_Gran_center = dev_gl_Gran_center2(:, :, now2)
        dev_gl_Cell_center = dev_gl_Cell_center2(:, :, now2)
        dev_gl_Gran_square = dev_gl_Gran_square2(:, now2)
		
	end subroutine Set_CUDA_move_reverse
		
	
	subroutine Send_data_to_Cuda()
		use GEO_PARAM
		use STORAGE
		 dev_gl_Gran_neighbour = gl_Gran_neighbour
		 dev_gl_Cell_par = gl_Cell_par
		 dev_gl_Cell_par_MF = gl_Cell_par_MF
		 dev_gl_Cell_center = gl_Cell_center
		 dev_gl_Gran_normal = gl_Gran_normal
		 dev_gl_Gran_square = gl_Gran_square
		 dev_gl_Gran_center = gl_Gran_center
		 dev_gl_Cell_dist = gl_Cell_dist
		 dev_gl_Cell_Volume = gl_Cell_Volume
		 dev_gl_Gran_info = gl_Gran_info
		 dev_gl_Gran_scheme = gl_Gran_scheme
		 dev_gl_Gran_type = gl_Gran_type
		 dev_gl_Cell_info = gl_Cell_info
		 dev_gl_Cell_gran = gl_Cell_gran
		 dev_gl_all_Gran_inner = gl_all_Gran_inner
		 dev_gl_all_Cell_inner = gl_all_Cell_inner
		 dev_gl_Cell_type = gl_Cell_type
		 dev_gl_Cell_number = gl_Cell_number
		 dev_gl_Gran_neighbour_TVD = gl_Gran_neighbour_TVD
		 
		 dev_gl_x = gl_x
		 dev_gl_y = gl_y
		 dev_gl_z = gl_z
		 
		 dev_ggg = ggg
		 dev_par_Velosity_inf = par_Velosity_inf
		 dev_Ngran = size(gl_all_Gran(1, :))
		 dev_Ncell = size(gl_Cell_par(1, :))
		 dev_Ngran_inner = size(gl_all_Gran_inner(:))
		 dev_Ncell_inner = size(gl_all_Cell_inner(:))
		 dev_par_R0 = par_R0
		 dev_par_a_2 = par_a_2
		 dev_par_n_p_LISM = par_n_p_LISM
		 dev_par_n_H_LISM_ = par_n_H_LISM_
		 dev_par_Kn = par_Kn
		 dev_par_R_inner = par_R_inner
		 dev_par_pi_8 = par_pi_8
		 dev_par_n_TS = par_n_TS
		 dev_par_n_BS = par_n_BS
		 dev_par_n_HP = par_n_HP
		 dev_par_n_END = par_n_END
		 dev_par_triple_point = par_triple_point
		 dev_par_m_BC = par_m_BC
		 dev_par_n_IA = par_n_IA
		 dev_par_n_IB = par_n_IB
		 dev_par_kk1 = par_kk1    
	     dev_par_kk12 = par_kk12    
         dev_par_kk2 = par_kk2     
         dev_par_kk3 = par_kk3    
	     dev_par_kk31 = par_kk31
		 dev_par_R_END = par_R_END
		 dev_par_R_LEFT = par_R_LEFT
		 dev_par_al1 = par_al1
		 
	end subroutine Send_data_to_Cuda
	
	subroutine Send_data_to_Host()
		use GEO_PARAM
		use STORAGE
		gl_Cell_par = dev_gl_Cell_par
		gl_Cell_par_MF = dev_gl_Cell_par_MF
	end subroutine Send_data_to_Host
	
	subroutine Send_data_to_Host_move(now)
		use GEO_PARAM
		use STORAGE
		implicit none
		integer, intent(in) :: now
		gl_x = dev_gl_x2(:, now)
        gl_y = dev_gl_y2(:, now)
        gl_z = dev_gl_z2(:, now)
        gl_Cell_Volume = dev_gl_Cell_Volume2(:, now)
        gl_Gran_normal = dev_gl_Gran_normal2(:, :, now)
        gl_Gran_center = dev_gl_Gran_center2(:, :, now)
        gl_Cell_center = dev_gl_Cell_center2(:, :, now)
        gl_Gran_square = dev_gl_Gran_square2(:, now)
		par_al1 = dev_par_al1
		gl_Gran_neighbour_TVD = dev_gl_Gran_neighbour_TVD
		
	end subroutine Send_data_to_Host_move
	
	
	attributes(device) real(8) function dev_norm2(x)
	! РАБОТАЕТ ТОЛЬКО ДЛЯ ТРЁХМЕРНЫХ ВЕКТОРОВ
    implicit none
	real(8), device :: x(*)

    dev_norm2 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)
	return
	end function dev_norm2
	
	attributes(device) subroutine dev_Calc_sourse_MF(plasma, fluid, sourse, zone)  ! Считаются мультифлюидные источники
    use GEO_PARAM
	implicit none
    real(8), intent(in) :: plasma(9)
    real(8), intent(in) :: fluid(5,4)
    real(8), intent(out) :: sourse(5,5)
    integer(4), intent(in) :: zone
    
    integer(4) :: i, kk(4)
    real(8) :: U_M_H(4), U_H(4), sigma(4), nu(4), S1, S2
	
    sourse = 0.0
    kk = 0
    kk(zone) = 1
    S1 = 0.0
    S2 = 0.0
    
    ! Body of Calc_sourse_MF
    do i = 1, 4
    	U_M_H(i) = sqrt( (plasma(2) - fluid(2, i))**2 + (plasma(3) - fluid(3, i))**2 + (plasma(4) - fluid(4, i))**2 + &
           (64.0 / (9.0 * dev_par_pi_8)) * (plasma(5) / plasma(1) + 2.0 * fluid(5, i) / fluid(1, i)) )
        U_H(i) = sqrt( (plasma(2) - fluid(2, i))**2 + (plasma(3) - fluid(3, i))**2 + (plasma(4) - fluid(4, i))**2 + &
           (4.0 / dev_par_pi_8) * (plasma(5) / plasma(1) + 2.0 * fluid(5, i) / fluid(1, i)) )
        sigma(i) = (1.0 - par_a_2 * log(U_M_H(i)))**2
        nu(i) = plasma(1) * fluid(1, i) * U_M_H(i) * sigma(i)
    end do
    
    do i = 1, 4
        sourse(2, 1) =  sourse(2, 1) + nu(i) * (fluid(2, i) - plasma(2))
        sourse(3, 1) =  sourse(3, 1) + nu(i) * (fluid(3, i) - plasma(3))
        sourse(4, 1) =  sourse(4, 1) + nu(i) * (fluid(4, i) - plasma(4))
        sourse(5, 1) = sourse(5, 1) + nu(i) * ( (fluid(2, i)**2 + fluid(3, i)**2 + fluid(4, i)**2 - &
            plasma(2)**2 - plasma(3)**2 - plasma(4)**2)/2.0 + (U_H(i)/U_M_H(i)) * ( 2.0 * fluid(5, i)/fluid(1, i) - plasma(5)/plasma(1) ) )
    end do
    
    sourse(2:5, 1) =  sourse(2:5, 1) * (par_n_p_LISM/par_Kn)
    
    do i = 1, 4
        S1 = S1 + nu(i)
        S2 = S2 + nu(i) * ( (plasma(2)**2 + plasma(3)**2 + plasma(4)**2)/2.0 + (U_H(i)/U_M_H(i)) * (plasma(5)/plasma(1)) )
    end do
    
    do i = 1, 4
        sourse(1, i + 1) = (par_n_H_LISM_/par_Kn) * (kk(i) * S1 - nu(i))
        sourse(2, i + 1) = (par_n_H_LISM_/par_Kn) * (kk(i) * S1 * plasma(2) - nu(i) * fluid(2, i))
        sourse(3, i + 1) = (par_n_H_LISM_/par_Kn) * (kk(i) * S1 * plasma(3) - nu(i) * fluid(3, i))
        sourse(4, i + 1) = (par_n_H_LISM_/par_Kn) * (kk(i) * S1 * plasma(4) - nu(i) * fluid(4, i))
        sourse(5, i + 1) = (par_n_H_LISM_/par_Kn) * (kk(i) * S2 - nu(i) * ( (fluid(2, i)**2 + fluid(3, i)**2 + fluid(4, i)**2)/2.0 + &
              (U_H(i)/U_M_H(i)) * 2.0 * (fluid(5, i) / fluid(1, i)) ) )
	end do
    
	end subroutine dev_Calc_sourse_MF
	
	attributes(global) subroutine dev_retime()
	    time_step = 10000.0
	end subroutine dev_retime
	
	end module MY_CUDA
	
	include "cuf_Solvers.cuf"
	 
	attributes(global) subroutine CUF_hellow()
	
	print*, "Hellow from CUDA"
	
	end subroutine CUF_hellow
	
	
	
	subroutine CUDA_START_MGD_move()
	! На подвижной сетке! с мультифлюидом
	use STORAGE
    use GEO_PARAM
	use MY_CUDA
	implicit none
    integer :: step, now, now2, step2, i, alla2(100), Num
	integer(4):: ierrSync, ierrAsync, nx, ny, ijk, istat
	integer(4), device :: dev_now, dev_now2
	real :: time_work
	real(8) :: local1
	type(dim3) :: grid, tBlock
	type(cudaEvent) :: startEvent, stopEvent
	
	call Set_CUDA()
	call Send_data_to_Cuda()
	call Alloc_CUDA_move()
	call Set_CUDA_move()
	
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc start 0: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync start 0: ', cudaGetErrorString(cudaGetLastError())
	
	tBlock = dim3(32, 8, 1)
	
	now = 2
	time_all = 0.0_8                 ! Глобальное время (хранится на девайсе)
	time_step2 = 0.00002_8              ! Шаг по времени (хранится на девайсе)
	
	istat = cudaEventCreate(startEvent)
	istat = cudaEventCreate(stopEvent)
	istat = cudaEventRecord(startEvent, 0)
	
	step = dev_Ngran
	print*, step, "   ", size(gl_all_Gran(1, :))
	step = dev_Ncell
	print*, step, "   ", size(gl_all_Cell(1, :))
	
	
	dev_gl_Point_num = 0.0
	
	! Главный цикл
	do step = 1,  50000  ! ---------------------------------------------------------------------------------------------------
		ierrAsync = cudaDeviceSynchronize()
		if (mod(step, 1000) == 0) then
			local1 = time_step2
			print*, "Step = ", step , "  step_time = ", local1
		end if
		
		
		
	!if(par_al1 > 0.3) par_al1 = par_al1 - 0.0000002
	!dev_par_al1 = par_al1
		
	
	time_step = time_step2
	!$cuf kernel do <<<1,1>>>
	do ijk = 1, 1
		time_all =  time_all + time_step
	end do
	
	time_step2 = 1000000_8       ! Шаг по времени (хранится на девайсе)
	
	now2 = now
	now = mod(now, 2) + 1
	ierrAsync = cudaDeviceSynchronize()
	dev_now = now
	dev_now2 = now2
	ierrAsync = cudaDeviceSynchronize()

	
	if(.False.) then  ! Есть ли вообще движение сетки
	! Сначала вычисляем скорости движения поверхностей
	
	
	! Вычисляем движения узлов на каждой поверхности (из задачи о распаде разрыва) 
	! *** аналог функции Calc_move на хосте ***
	Num = size(gl_TS)
	call Cuda_Calc_move_TS<<<ceiling(real(Num)/256), 256>>>(dev_now)
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) write (*,*) 'Error Sinc start 1: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) write(*,*) 'Error ASync start 1: ', cudaGetErrorString(cudaGetLastError())
	
	Num = size(gl_Contact)
	call Cuda_Calc_move_HP<<<ceiling(real(Num)/256), 256>>>(dev_now)
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) write (*,*) 'Error Sinc start 2: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) write(*,*) 'Error ASync start 2: ', cudaGetErrorString(cudaGetLastError())
	
	
	Num = size(gl_BS)
	call Cuda_Calc_move_BS<<<ceiling(real(Num)/256), 256>>>(dev_now)
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) write (*,*) 'Error Sinc start 3: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) write(*,*) 'Error ASync start 3: ', cudaGetErrorString(cudaGetLastError())
	
	!dev_gl_Vx = 0.0
	!dev_gl_Vy = 0.0
	!dev_gl_Vz = 0.0
	
	
	! Делаем три цикла поверхностного натяжения (по разным лучам, как на хосте) -------------------------------------------------------------------------
	! *** аналог функции Move_all на хосте ***
	if (.True.) then
		nx = size(gl_RAY_A(1, 1, :))
		ny = size(gl_RAY_A(1, :, 1))
		grid = dim3( ceiling(real(nx)/tBlock%x), &
			ceiling(real(ny)/tBlock%y), 1)
	
		call Cuda_Move_all_1 <<< grid, tBlock>>> (dev_now)
		ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
			write (*,*) 'Error Sinc start 4: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
			write(*,*) 'Error ASync start 4: ', cudaGetErrorString(cudaGetLastError())
	
	
		nx = size(gl_RAY_B(1, 1, :))
		ny = size(gl_RAY_B(1, :, 1))
		grid = dim3( ceiling(real(nx)/tBlock%x), &
			ceiling(real(ny)/tBlock%y), 1)
	
		call Cuda_Move_all_2 <<< grid, tBlock>>> (dev_now)
		ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
			write (*,*) 'Error Sinc start 5: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
			write(*,*) 'Error ASync start 5: ', cudaGetErrorString(cudaGetLastError())
	
		nx = size(gl_RAY_K(1, 1, :))
		ny = size(gl_RAY_K(1, :, 1))
		grid = dim3( ceiling(real(nx)/tBlock%x), &
			ceiling(real(ny)/tBlock%y), 1)
	
		call Cuda_Move_all_3 <<< grid, tBlock>>> (dev_now)
		ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
			write (*,*) 'Error Sinc start 6: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
			write(*,*) 'Error ASync start 6: ', cudaGetErrorString(cudaGetLastError())
	
	
		nx = size(gl_RAY_O(1, 1, :))
		ny = size(gl_RAY_O(1, :, 1))
		grid = dim3( ceiling(real(nx)/tBlock%x), &
			ceiling(real(ny)/tBlock%y), 1)
	 
		call Cuda_Move_all_4 <<< grid, tBlock>>> (dev_now)
		ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
			write (*,*) 'Error Sinc start 6: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
			write(*,*) 'Error ASync start 6: ', cudaGetErrorString(cudaGetLastError())
	
	end if
	
	! Здесь нужно подвинуть узел на TS на оси симметрии (найти его скорость скорее)
	
	!call Cuda_Move_all_5 <<< 1, 1>>> (dev_now)
	!	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
	!		write (*,*) 'Error Sinc start 6: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
	!		write(*,*) 'Error ASync start 6: ', cudaGetErrorString(cudaGetLastError())
	
	! Теперь собственно циклы движения точек ---------------------------------------------------------------------------------------------------------
	
	nx = size(gl_RAY_A(1, 1, :))
    ny = size(gl_RAY_A(1, :, 1))
	grid = dim3( ceiling(real(nx)/tBlock%x), &
		ceiling(real(ny)/tBlock%y), 1)
	
	call Cuda_Move_all_A <<< grid, tBlock>>> (dev_now)
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc start 7: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync start 7: ', cudaGetErrorString(cudaGetLastError())
	
	
	nx = size(gl_RAY_B(1, 1, :))
    ny = size(gl_RAY_B(1, :, 1))
	grid = dim3( ceiling(real(nx)/tBlock%x), &
		ceiling(real(ny)/tBlock%y), 1)
	
	call Cuda_Move_all_B <<< grid, tBlock>>> (dev_now)
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc start 8: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync start 8: ', cudaGetErrorString(cudaGetLastError())
	
	nx = size(gl_RAY_C(1, 1, :))
    ny = size(gl_RAY_C(1, :, 1))
	grid = dim3( ceiling(real(nx)/tBlock%x), &
		ceiling(real(ny)/tBlock%y), 1)
	
	call Cuda_Move_all_C <<< grid, tBlock>>> (dev_now)
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc start 9: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync start 9: ', cudaGetErrorString(cudaGetLastError())
	
	nx = size(gl_RAY_O(1, 1, :))
    ny = size(gl_RAY_O(1, :, 1))
	grid = dim3( ceiling(real(nx)/tBlock%x), &
		ceiling(real(ny)/tBlock%y), 1)
	
	call Cuda_Move_all_O <<< grid, tBlock>>> (dev_now)
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc start 10: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync start 10: ', cudaGetErrorString(cudaGetLastError())
	
	
	nx = size(gl_RAY_K(1, 1, :))
    ny = size(gl_RAY_K(1, :, 1))
	grid = dim3( ceiling(real(nx)/tBlock%x), &
		ceiling(real(ny)/tBlock%y), 1)
	
	call Cuda_Move_all_K <<< grid, tBlock>>> (dev_now)
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc start 11: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync start 11: ', cudaGetErrorString(cudaGetLastError())
	
	
	nx = size(gl_RAY_D(1, 1, :))
    ny = size(gl_RAY_D(1, :, 1))
	grid = dim3( ceiling(real(nx)/tBlock%x), &
		ceiling(real(ny)/tBlock%y), 1)
	
	call Cuda_Move_all_D <<< grid, tBlock>>> (dev_now)
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc start 12: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync start 12: ', cudaGetErrorString(cudaGetLastError())
	
	
	nx = size(gl_RAY_E(1, 1, :))
    ny = size(gl_RAY_E(1, :, 1))
	grid = dim3( ceiling(real(nx)/tBlock%x), &
		ceiling(real(ny)/tBlock%y), 1)
	
	call Cuda_Move_all_E <<< grid, tBlock>>> (dev_now)
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc start 13: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync start 13: ', cudaGetErrorString(cudaGetLastError())
	
	
	! Теперь посчитаем новые объёмы и площади граней  ---------------------------------------------------------------------------------------------------------
	! *** аналог функции calc_all_Gran_move на хосте ***
	
	Num = size(gl_all_Gran(1,:))
	call Cuda_calc_all_Gran_move_1 <<< ceiling(real(Num)/256), 256>>> (dev_now2)  ! цикл по граням
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc start 14: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync start 14: ', cudaGetErrorString(cudaGetLastError())
	
	Num = size(gl_all_Cell(1,:))
	call Cuda_calc_all_Gran_move_2 <<< ceiling(real(Num)/256), 256>>> (dev_now2)  ! цикл по ячейкам
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc start 15: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync start 15: ', cudaGetErrorString(cudaGetLastError())
	
	
	! Для сравнения результатов с хостом
	!gl_Gran_square = dev_gl_Gran_square2(:, now2) 
	!gl_Gran_normal = dev_gl_Gran_normal2(:, :, now2) 
	!gl_Cell_Volume = dev_gl_Cell_Volume2(:, now2) 
	!gl_Gran_center = dev_gl_Gran_center2(:, :, now2) 
	!
	!print*, gl_Gran_center(:, 30)
	!print*, gl_Gran_center(:, 90)
	!print*, gl_Gran_center(:, 830)
	!print*, "_______"
	!print*, gl_Gran_square(10)
	!print*, gl_Gran_square(50)
	!print*, gl_Gran_square(300)
	!print*, "_______"
	!print*, gl_Cell_Volume(10)
	!print*, gl_Cell_Volume(50)
	!print*, gl_Cell_Volume(300)
	!print*, "_______"
	!print*, gl_Gran_normal(:, 10)
	!print*, gl_Gran_normal(:, 50)
	!print*, gl_Gran_normal(:, 300)
	!print*, "_______"
	
	end if
	
	! Теперь нужен цикл по граням для решения задачи о распаде произвольного разрыва на них
	
	
	Num = size(gl_all_Gran(1, :))
	call CUF_MGD_grans_MF <<< ceiling(real(Num)/256), 256>>> (dev_now)  ! цикл по граням
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc start 16: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync start 16: ', cudaGetErrorString(cudaGetLastError())

	
	Num = size(gl_all_Cell(1, :))
	call CUF_MGD_cells_MF <<< ceiling(real(Num)/256), 256>>> (dev_now)  ! цикл по ячейкам
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc start 17: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync start 17: ', cudaGetErrorString(cudaGetLastError())
	
	!gl_Cell_par = dev_gl_Cell_par
	
	!print*, "____"
	!print*, gl_Cell_par(:, 10)
	!print*, "____"
	!print*, gl_Cell_par(:, 1000)
	!print*, "____"
	!print*, "____"
	
	
	! Здесь нужно два цикла по внутренней области (сфере) ячеек и граней
	
	dev_gl_x = dev_gl_x2(:, now2)
    dev_gl_y = dev_gl_y2(:, now2)
    dev_gl_z = dev_gl_z2(:, now2)
    dev_gl_Cell_Volume = dev_gl_Cell_Volume2(:, now2)
    dev_gl_Gran_normal = dev_gl_Gran_normal2(:, :, now2)
    dev_gl_Gran_center = dev_gl_Gran_center2(:, :, now2)
    dev_gl_Cell_center = dev_gl_Cell_center2(:, :, now2)
    dev_gl_Gran_square = dev_gl_Gran_square2(:, now2)
	
	! Нужно обновить граничные условия на внутренней сфере (так как она немного подвигалась)
	
	!if (mod(step, 100) == 0) then
	!		gl_Cell_center = dev_gl_Cell_center
	!		gl_Cell_par = dev_gl_Cell_par
	!		gl_Cell_par_MF = dev_gl_Cell_par_MF
	!		call Initial_conditions()
	!		dev_gl_Cell_par = gl_Cell_par
	!		dev_gl_Cell_par_MF = gl_Cell_par_MF
	!end if
	
	if (.True.) then
		do ijk = 1, 5  ! Несколько раз просчитываем внутреннюю область
			dev_time_step_inner = 1000000.0
			ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
				write (*,*) 'Error Sinc start 18: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
				write(*,*) 'Error ASync start 18: ', cudaGetErrorString(cudaGetLastError())
		
			Num = size(gl_all_Gran_inner(:))
			call CUF_MGD_grans_MF_inner <<< ceiling(real(Num)/256), 256>>> ()  ! цикл по граням
			ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
				write (*,*) 'Error Sinc start 19: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
				write(*,*) 'Error ASync start 19: ', cudaGetErrorString(cudaGetLastError())
 
			Num = size(gl_all_Cell_inner(:))
			call CUF_MGD_cells_MF_inner <<< ceiling(real(Num)/256), 256>>> ()  ! цикл по ячейкам
			ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
				write (*,*) 'Error Sinc start 20: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
				write(*,*) 'Error ASync start 20: ', cudaGetErrorString(cudaGetLastError())
		end do
	end if
	
	
	dev_gl_Vx = 0.0
    dev_gl_Vy = 0.0
    dev_gl_Vz = 0.0
    dev_gl_Point_num = 0
	
	
	if (mod(step, 20000) == 0 .or. step == 2000 .or. step == 10000 .or. step == 5000 .or. step == 15000) then
		print*, "PECHAT ", step
		par_al1 = dev_par_al1
		print*, "par_al1 = ", par_al1
		call Send_data_to_Host_move(now2)
		call Send_data_to_Host()
		call Print_surface_2D()
		call Print_Setka_2D()
		call Print_par_2D()
		call Print_par_y_2D()
		call Print_Contact_3D()
		call Print_TS_3D()
		call Print_Setka_y_2D()
		call Print_surface_y_2D()
		call Print_Cell(gl_Cell_number(1, 207606), gl_Cell_number(2, 207606), gl_Cell_number(3, 207606), gl_Cell_type(207606))
	end if
	
	if (mod(step, 100000) == 0) then
		print*, "PECHAT 2 ", step
		call Send_data_to_Host_move(now2)
		call Send_data_to_Host()
		call Save_setka_bin(79)
	end if
	
	if (mod(step, 5000) == 0) then
		print*, "Renew TVD ", step
		call Send_data_to_Host_move(now2)
		call Send_data_to_Host()
		call Find_TVD_sosed()
		dev_gl_Gran_neighbour_TVD = gl_Gran_neighbour_TVD
	end if
	
	
	end do
	
	istat = cudaEventRecord(stopEvent, 0)
	istat = cudaEventSynchronize(stopEvent)
	istat = cudaEventElapsedTime(time_work, startEvent, stopEvent)
	print *, "CUDA Time work: ", (time_work)/(60*1000.0), "   in minutes"
	
	par_al1 = dev_par_al1
	print*, "par_al1 = ", par_al1
	call Send_data_to_Host_move(now2)
	call Send_data_to_Host()
	call Print_surface_2D()
    call Print_Setka_2D()
    call Print_par_2D()
	call Print_par_y_2D()
	call Print_Setka_y_2D()
	call Print_TS_3D()
	call Print_Contact_3D
	call Print_surface_y_2D()
	
	
	end subroutine CUDA_START_MGD_move
	
	
	subroutine CUDA_START_GD_3()
	! На неподвижной сетке! Без магнитных полей, с мультифлюидом
	use STORAGE
    use GEO_PARAM
	use MY_CUDA
	
	implicit none
	
	integer(4):: NGRAN, NCELL, NGRAN_INNER, NCELL_INNER, ierrSync, ierrAsync, step, st, step_all
	
	NGRAN = size(gl_all_Gran(1, :))
	NCELL = size(gl_Cell_Volume(:))
	NGRAN_INNER = size(gl_all_Gran_inner(:))
    NCELL_INNER = size(gl_all_Cell_inner(:))
	
	print*, NGRAN, NCELL, NGRAN_INNER, NCELL_INNER
	print*, ceiling(real(NGRAN)/256), ceiling(real(NCELL)/256), ceiling(real(NGRAN_INNER)/256), ceiling(real(NCELL_INNER)/256)
	
	print*, "DO = ", gl_Cell_par(1, 100)
	
	call Set_CUDA()
	call Send_data_to_Cuda()
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc start: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync start: ', cudaGetErrorString(cudaGetLastError())
	
	step_all = 20000 * 2 * 6
	
	do step = 1, step_all
		
		if(mod(step, 1000) == 0) then
			print*, "step = " , step , "   of   ", step_all
		end if
	
	    call CUF_GD_3_grans<<<ceiling(real(NGRAN)/256), 256>>>()
		ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc 1: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync 1: ', cudaGetErrorString(cudaGetLastError())
		
		call CUF_GD_3_cells<<<ceiling(real(NCELL)/256), 256>>>()
		ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc 2: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync 2: ', cudaGetErrorString(cudaGetLastError())
		
		call dev_retime<<<1, 1>>>()
		ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		    write (*,*) 'Error Sinc 3: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		    write(*,*) 'Error ASync 3: ', cudaGetErrorString(cudaGetLastError())
		
		do st = 1, 4
			
			call CUF_GD_3_grans_inner<<<ceiling(real(NGRAN_INNER)/256), 256>>>()
		    ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		    write (*,*) 'Error Sinc 4: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		    write(*,*) 'Error ASync 4: ', cudaGetErrorString(cudaGetLastError())
		
		    call CUF_GD_3_cells_inner<<<ceiling(real(NCELL_INNER)/256), 256>>>()
		    ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		    write (*,*) 'Error Sinc 5: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		    write(*,*) 'Error ASync 5: ', cudaGetErrorString(cudaGetLastError())
		
		    call dev_retime<<<1, 1>>>()
			ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		    write (*,*) 'Error Sinc 6: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		    write(*,*) 'Error ASync 6: ', cudaGetErrorString(cudaGetLastError())
			
		end do
		
	
	end do
	
	call Send_data_to_Host()
	ierrSync = cudaGetLastError(); ierrAsync = cudaDeviceSynchronize(); if (ierrSync /= cudaSuccess) &
		write (*,*) 'Error Sinc end: ', cudaGetErrorString(ierrSync); if(ierrAsync /= cudaSuccess) & 
		write(*,*) 'Error ASync end: ', cudaGetErrorString(cudaGetLastError())
	
	print*, "POSLE = ", gl_Cell_par(1, 100)
	
	end subroutine CUDA_START_GD_3
	
	attributes(global) subroutine CUF_GD_3_grans_inner()
	! Просчёт потоков через грани и запись этих потоков в массивы для использования далее
	! Неподвижная сетка, газовая динамика + мультифлюид
	! Программа, вызываемая с хоста, содержащая cuf-ядро
	use MY_CUDA, gl_Gran_info => dev_gl_Gran_info, gl_Gran_neighbour => dev_gl_Gran_neighbour, &
		gl_Cell_par => dev_gl_Cell_par, gl_Cell_par_MF => dev_gl_Cell_par_MF, gl_Cell_center => dev_gl_Cell_center, &
		gl_Gran_center => dev_gl_Gran_center, norm2 => dev_norm2, ggg => dev_ggg, par_Velosity_inf => dev_par_Velosity_inf, &
		gl_Cell_dist => dev_gl_Cell_dist, gl_Gran_POTOK => dev_gl_Gran_POTOK, gl_Gran_POTOK_MF => dev_gl_Gran_POTOK_MF, &
		gl_Gran_square => dev_gl_Gran_square, gl_Gran_normal => dev_gl_Gran_normal, gl_all_Gran_inner => dev_gl_all_Gran_inner
	!use cudadevice
	!use libm
	implicit none
	
	integer(4) :: gr
	
    integer(4) :: st, s1, s2, i, j, k, zone, iter
    real(8) :: qqq1(9), qqq2(9), qqq(9)  ! Переменные в ячейке
    real(8) :: fluid1(5, 4), fluid2(5, 4)
    real(8) :: dist, dsl, dsc, dsp
    real(8) :: POTOK(9), ttest(3)
    real(8) :: POTOK_MF(5)
    real(8) :: POTOK_MF_all(5, 4)
    real(8) :: time, Volume, TT, U8, rad1, rad2, aa, bb, cc
    real(8) :: SOURSE(5,5)  ! Источники массы, импульса и энергии для плазмы и каждого сорта мультифлюида
    real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3
	
	time = 100000.0
	iter = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	
	if (iter > dev_Ngran_inner) return
	
	!if(gr == 1) then
	!	print*, "Hellow from potok (1, 1) ", dev_Ngran
	!end if
	
	
	! ----------------------------------------------------------- можно просто скопировать код из хоста ----------------------
	gr = gl_all_Gran_inner(iter)
	
            POTOK = 0.0
            s1 = gl_Gran_neighbour(1, gr)
            s2 = gl_Gran_neighbour(2, gr)
            qqq1 = gl_Cell_par(:, s1)
            fluid1 = gl_Cell_par_MF(:, :, s1)   ! Загрузили параметры жидкостей для мультифлюида

            ! Попробуем снести плотность пропорционально квадрату
            if(norm2(qqq1(2:4))/sqrt(ggg*qqq1(5)/qqq1(1)) > 5.0) then
                rad1 = norm2(gl_Cell_center(:, s1))
                rad2 = norm2(gl_Gran_center(:, gr))
                qqq1(1) = qqq1(1) * rad1**2 / rad2**2
                qqq1(9) = qqq1(9) * rad1**2 / rad2**2
                qqq1(5) = qqq1(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
                ! Скорости сносим в сферической С.К.
                call spherical_skorost(gl_Cell_center(1, s1), gl_Cell_center(2, s1), gl_Cell_center(3, s1), &
                    qqq1(2), qqq1(3), qqq1(4), aa, bb, cc)
                call dekard_skorost(gl_Gran_center(1, gr), gl_Gran_center(2, gr), gl_Gran_center(3, gr), &
                    aa, bb, cc, qqq1(2), qqq1(3), qqq1(4))

                call spherical_skorost(gl_Cell_center(1, s1), gl_Cell_center(2, s1), gl_Cell_center(3, s1), &
                    fluid1(2, 1), fluid1(3, 1), fluid1(4, 1), aa, bb, cc)
                call dekard_skorost(gl_Gran_center(1, gr), gl_Gran_center(2, gr), gl_Gran_center(3, gr), &
                    aa, bb, cc, fluid1(2, 1), fluid1(3, 1), fluid1(4, 1))
            end if

            if (s2 >= 1) then
                !if ( norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int .and. norm2(gl_Cell_center(:, s2)) <= par_R0 * par_R_int) CYCLE
                qqq2 = gl_Cell_par(:, s2)
                fluid2 = gl_Cell_par_MF(:, :, s2)   ! Загрузили параметры жидкостей для мультифлюида
                dist = min(gl_Cell_dist(s1), gl_Cell_dist(s2))

                ! Попробуем снести плотность пропорционально квадрату
                if(norm2(qqq2(2:4))/sqrt(ggg*qqq2(5)/qqq2(1)) > 5.0) then
                    rad1 = norm2(gl_Cell_center(:, s2))
                    rad2 = norm2(gl_Gran_center(:, gr))
                    qqq2(1) = qqq2(1) * rad1**2 / rad2**2
                    qqq2(9) = qqq2(9) * rad1**2 / rad2**2
                    qqq2(5) = qqq2(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
                    call spherical_skorost(gl_Cell_center(1, s2), gl_Cell_center(2, s2), gl_Cell_center(3, s2), &
                        qqq2(2), qqq2(3), qqq2(4), aa, bb, cc)
                    call dekard_skorost(gl_Gran_center(1, gr), gl_Gran_center(2, gr), gl_Gran_center(3, gr), &
                        aa, bb, cc, qqq2(2), qqq2(3), qqq2(4))

                    call spherical_skorost(gl_Cell_center(1, s2), gl_Cell_center(2, s2), gl_Cell_center(3, s2), &
                        fluid2(2, 1), fluid2(3, 1), fluid2(4, 1), aa, bb, cc)
                    call dekard_skorost(gl_Gran_center(1, gr), gl_Gran_center(2, gr), gl_Gran_center(3, gr), &
                        aa, bb, cc, fluid2(2, 1), fluid2(3, 1), fluid2(4, 1))
                end if

            else  ! В случае граничных ячеек - граничные условия
                !if (norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int) CYCLE
                if(s2 == -1) then  ! Набегающий поток
                    dist = gl_Cell_dist(s1)
                    qqq2 = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 1.0_8, 0.0_8, 0.0_8, 0.0_8, 100.0_8/)
                    fluid2(:, 1) = (/0.000001_8, 0.0_8, 0.0_8, 0.0_8, 0.000001_8/)
                    fluid2(:, 2) = (/0.000001_8, 0.0_8, 0.0_8, 0.0_8, 0.000001_8/)
                    fluid2(:, 3) = (/0.000001_8, 0.0_8, 0.0_8, 0.0_8, 0.000001_8/)
                    fluid2(:, 4) = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 0.5_8/)
                else  ! Здесь нужны мягкие условия (это задняя стенка)
                    dist = gl_Cell_dist(s1)
                    qqq2 = qqq1
                    fluid2 = fluid1
                    qqq2(5) = 1.0_8
                    if(qqq2(2) > par_Velosity_inf) then
                        qqq2(2) = par_Velosity_inf ! Отсос жидкости
                    end if

                end if
            end if


            call chlld_Q(2, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)
            time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )   ! REDUCTION
            gl_Gran_POTOK(:, gr) = POTOK * gl_Gran_square(gr)

            call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                0.0_8, fluid1(:, 1), fluid2(:, 1), dsl, dsp, dsc, POTOK_MF)
            time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )
            gl_Gran_POTOK_MF(:, 1, gr) = POTOK_MF * gl_Gran_square(gr)

            call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                0.0_8, fluid1(:, 2), fluid2(:, 2), dsl, dsp, dsc, POTOK_MF)
            time = min(time, 0.8 * dist/max(dabs(dsl), dabs(dsp)) )
            gl_Gran_POTOK_MF(:, 2, gr) = POTOK_MF * gl_Gran_square(gr)

            call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                0.0_8, fluid1(:, 3), fluid2(:, 3), dsl, dsp, dsc, POTOK_MF)
            time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )
            gl_Gran_POTOK_MF(:, 3, gr) = POTOK_MF * gl_Gran_square(gr)

            call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                0.0_8, fluid1(:, 4), fluid2(:, 4), dsl, dsp, dsc, POTOK_MF)
            time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )
            gl_Gran_POTOK_MF(:, 4, gr) = POTOK_MF * gl_Gran_square(gr)
			
			
			if (time_step > time) then 
			time =  atomicmin(time_step, time)   ! Атомарная операция взятия минимального значения
			end if

	
	
	end subroutine CUF_GD_3_grans_inner
	
	attributes(global) subroutine CUF_GD_3_grans()
	! Просчёт потоков через грани и запись этих потоков в массивы для использования далее
	! Неподвижная сетка, газовая динамика + мультифлюид
	! Программа, вызываемая с хоста, содержащая cuf-ядро
	use MY_CUDA, gl_Gran_info => dev_gl_Gran_info, gl_Gran_neighbour => dev_gl_Gran_neighbour, &
		gl_Cell_par => dev_gl_Cell_par, gl_Cell_par_MF => dev_gl_Cell_par_MF, gl_Cell_center => dev_gl_Cell_center, &
		gl_Gran_center => dev_gl_Gran_center, norm2 => dev_norm2, ggg => dev_ggg, par_Velosity_inf => dev_par_Velosity_inf, &
		gl_Cell_dist => dev_gl_Cell_dist, gl_Gran_POTOK => dev_gl_Gran_POTOK, gl_Gran_POTOK_MF => dev_gl_Gran_POTOK_MF, &
		gl_Gran_square => dev_gl_Gran_square, gl_Gran_normal => dev_gl_Gran_normal
	!use cudadevice
	!use libm
	implicit none
	
	integer(4) :: gr
	
    integer(4) :: st, s1, s2, i, j, k, zone
    real(8) :: qqq1(9), qqq2(9), qqq(9)  ! Переменные в ячейке
    real(8) :: fluid1(5, 4), fluid2(5, 4)
    real(8) :: dist, dsl, dsc, dsp
    real(8) :: POTOK(9), ttest(3)
    real(8) :: POTOK_MF(5)
    real(8) :: POTOK_MF_all(5, 4)
    real(8) :: time, Volume, TT, U8, rad1, rad2, aa, bb, cc
    real(8) :: SOURSE(5,5)  ! Источники массы, импульса и энергии для плазмы и каждого сорта мультифлюида
    real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3
	
	time = 100000.0
	gr = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	
	!if(gr == 1) then
	!	print*, "Hellow from potok (1, 1) ", dev_Ngran
	!end if
	
	
	! ----------------------------------------------------------- можно просто скопировать код из хоста ----------------------
	if (gr > dev_Ngran) return
	
	if(gl_Gran_info(gr) == 2) return
	
	s1 = gl_Gran_neighbour(1, gr)
	s2 = gl_Gran_neighbour(2, gr)
	qqq1 = gl_Cell_par(:, s1)
	fluid1 = gl_Cell_par_MF(:, :, s1)   ! Загрузили параметры жидкостей для мультифлюида
	
	! Попробуем снести плотность пропорционально квадрату
    if(norm2(qqq1(2:4))/sqrt(ggg*qqq1(5)/qqq1(1)) > 2.2) then
        rad1 = norm2(gl_Cell_center(:, s1))
        rad2 = norm2(gl_Gran_center(:, gr))
        qqq1(1) = qqq1(1) * rad1**2 / rad2**2
        qqq1(9) = qqq1(9) * rad1**2 / rad2**2
        qqq1(5) = qqq1(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
        ! Скорости сносим в сферической С.К.
        call spherical_skorost(gl_Cell_center(1, s1), gl_Cell_center(2, s1), gl_Cell_center(3, s1), &
            qqq1(2), qqq1(3), qqq1(4), aa, bb, cc)
        call dekard_skorost(gl_Gran_center(1, gr), gl_Gran_center(2, gr), gl_Gran_center(3, gr), &
            aa, bb, cc, qqq1(2), qqq1(3), qqq1(4))

        call spherical_skorost(gl_Cell_center(1, s1), gl_Cell_center(2, s1), gl_Cell_center(3, s1), &
            fluid1(2, 1), fluid1(3, 1), fluid1(4, 1), aa, bb, cc)
        call dekard_skorost(gl_Gran_center(1, gr), gl_Gran_center(2, gr), gl_Gran_center(3, gr), &
            aa, bb, cc, fluid1(2, 1), fluid1(3, 1), fluid1(4, 1))
	end if
	
	if (s2 >= 1) then
                !if ( norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int .and. norm2(gl_Cell_center(:, s2)) <= par_R0 * par_R_int) CYCLE
                qqq2 = gl_Cell_par(:, s2)
                fluid2 = gl_Cell_par_MF(:, :, s2)   ! Загрузили параметры жидкостей для мультифлюида
                dist = min(gl_Cell_dist(s1), gl_Cell_dist(s2))

                ! Попробуем снести плотность пропорционально квадрату
                if(norm2(qqq2(2:4))/sqrt(ggg*qqq2(5)/qqq2(1)) > 2.2) then
                    rad1 = norm2(gl_Cell_center(:, s2))
                    rad2 = norm2(gl_Gran_center(:, gr))
                    qqq2(1) = qqq2(1) * rad1**2 / rad2**2
                    qqq2(9) = qqq2(9) * rad1**2 / rad2**2
                    qqq2(5) = qqq2(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
                    call spherical_skorost(gl_Cell_center(1, s2), gl_Cell_center(2, s2), gl_Cell_center(3, s2), &
                        qqq2(2), qqq2(3), qqq2(4), aa, bb, cc)
                    call dekard_skorost(gl_Gran_center(1, gr), gl_Gran_center(2, gr), gl_Gran_center(3, gr), &
                        aa, bb, cc, qqq2(2), qqq2(3), qqq2(4))

                    call spherical_skorost(gl_Cell_center(1, s2), gl_Cell_center(2, s2), gl_Cell_center(3, s2), &
                        fluid2(2, 1), fluid2(3, 1), fluid2(4, 1), aa, bb, cc)
                    call dekard_skorost(gl_Gran_center(1, gr), gl_Gran_center(2, gr), gl_Gran_center(3, gr), &
                        aa, bb, cc, fluid2(2, 1), fluid2(3, 1), fluid2(4, 1))
                end if

            else  ! В случае граничных ячеек - граничные условия
                !if (norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int) CYCLE
                if(s2 == -1) then  ! Набегающий поток
                    dist = gl_Cell_dist(s1)
                    qqq2 = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 1.0_8, 0.0_8, 0.0_8, 0.0_8, 100.0_8/)
                    fluid2(:, 1) = fluid1(:, 1)
                    fluid2(:, 2) = fluid1(:, 2)
                    fluid2(:, 3) = fluid1(:, 3)
                    fluid2(:, 4) = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 0.5_8/)
                else  ! Здесь нужны мягкие условия
                    dist = gl_Cell_dist(s1)
                    qqq2 = qqq1
                    fluid2 = fluid1
                    qqq2(5) = 1.0_8
                    if(qqq2(2) > par_Velosity_inf) then
                        qqq2(2) = par_Velosity_inf ! Отсос жидкости
                    end if

                end if
	end if
	
	call chlld_Q(2, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)
    time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )   ! REDUCTION
    gl_Gran_POTOK(:, gr) = POTOK * gl_Gran_square(gr)

    call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
        0.0_8, fluid1(:, 1), fluid2(:, 1), dsl, dsp, dsc, POTOK_MF)
    time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )
    gl_Gran_POTOK_MF(:, 1, gr) = POTOK_MF * gl_Gran_square(gr)

    call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
        0.0_8, fluid1(:, 2), fluid2(:, 2), dsl, dsp, dsc, POTOK_MF)
    time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )
    gl_Gran_POTOK_MF(:, 2, gr) = POTOK_MF * gl_Gran_square(gr)

    call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
        0.0_8, fluid1(:, 3), fluid2(:, 3), dsl, dsp, dsc, POTOK_MF)
    time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )
    gl_Gran_POTOK_MF(:, 3, gr) = POTOK_MF * gl_Gran_square(gr)

    call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
        0.0_8, fluid1(:, 4), fluid2(:, 4), dsl, dsp, dsc, POTOK_MF)
    time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )
    gl_Gran_POTOK_MF(:, 4, gr) = POTOK_MF * gl_Gran_square(gr)
	
	! ----------------------------------------------------------- копируем до этого момента ----------------------
	
	if (time_step > time) then 
		time =  atomicmin(time_step, time)   ! Атомарная операция взятия минимального значения
	end if
	
	!if(gr <= 200) then
	!	print*, gr, time, time_step
	!end if
	
	
	end subroutine CUF_GD_3_grans
	
	
	attributes(global) subroutine CUF_GD_3_cells_inner()
	! Неподвижная сетка, газовая динамика + мультифлюид
	use MY_CUDA, gl_Gran_info => dev_gl_Gran_info, gl_Gran_neighbour => dev_gl_Gran_neighbour, &
		gl_Cell_par => dev_gl_Cell_par, gl_Cell_par_MF => dev_gl_Cell_par_MF, gl_Cell_center => dev_gl_Cell_center, &
		gl_Gran_center => dev_gl_Gran_center, norm2 => dev_norm2, ggg => dev_ggg, par_Velosity_inf => dev_par_Velosity_inf, &
		gl_Cell_dist => dev_gl_Cell_dist, gl_Gran_POTOK => dev_gl_Gran_POTOK, gl_Gran_POTOK_MF => dev_gl_Gran_POTOK_MF, &
		gl_Gran_square => dev_gl_Gran_square, gl_Gran_normal => dev_gl_Gran_normal, &
		gl_Cell_info => dev_gl_Cell_info, gl_Cell_Volume => dev_gl_Cell_Volume, gl_Cell_gran => dev_gl_Cell_gran, &
		par_R0 => dev_par_R0, gl_all_Cell_inner => dev_gl_all_Cell_inner, gl_Cell_number=>dev_gl_Cell_number, &
		gl_Cell_type=>dev_gl_Cell_type
	!use cudadevice
	!use libm
	implicit none
	
	integer(4) :: gr, iter
	
    integer(4) :: st, s1, s2, i, j, k, zone
    real(8) :: qqq1(9), qqq2(9), qqq(9)  ! Переменные в ячейке
    real(8) :: fluid1(5, 4), fluid2(5, 4)
    real(8) :: dist, dsl, dsc, dsp
    real(8) :: POTOK(9), ttest(3)
    real(8) :: POTOK_MF(5)
    real(8) :: POTOK_MF_all(5, 4)
    real(8) :: time, Volume, TT, U8, rad1, rad2, aa, bb, cc
    real(8) :: SOURSE(5,5)  ! Источники массы, импульса и энергии для плазмы и каждого сорта мультифлюида
    real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3
	logical :: l_1
	
	iter = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	time = time_step
	
	if (iter > dev_Ncell_inner) return
	
	gr = gl_all_Cell_inner(iter)
            !if(gl_Cell_info(gr) == 2) CYCLE
            l_1 = .TRUE.
            if ((gl_Cell_type(gr) == "A" .or. gl_Cell_type(gr) == "B").and.(gl_Cell_number(1, gr) <= 2) ) l_1 = .FALSE.    ! Не считаем в первых двух ячейках
            POTOK = 0.0
            SOURSE = 0.0
            POTOK_MF_all = 0.0
            Volume = gl_Cell_Volume(gr)
            qqq = gl_Cell_par(:, gr)
            fluid1 = gl_Cell_par_MF(:, :, gr)
            ! Просуммируем потоки через грани
            do i = 1, 6
                j = gl_Cell_gran(i, gr)
                if (j == 0) CYCLE
                if (gl_Gran_neighbour(1, j) == gr) then
                    POTOK = POTOK + gl_Gran_POTOK(:, j)
                    POTOK_MF_all(:, 1) = POTOK_MF_all(:, 1) +  gl_Gran_POTOK_MF(:, 1, j)
                    POTOK_MF_all(:, 2) = POTOK_MF_all(:, 2) +  gl_Gran_POTOK_MF(:, 2, j)
                    POTOK_MF_all(:, 3) = POTOK_MF_all(:, 3) +  gl_Gran_POTOK_MF(:, 3, j)
                    POTOK_MF_all(:, 4) = POTOK_MF_all(:, 4) +  gl_Gran_POTOK_MF(:, 4, j)
                else
                    POTOK = POTOK - gl_Gran_POTOK(:, j)
                    POTOK_MF_all(:, 1) = POTOK_MF_all(:, 1) -  gl_Gran_POTOK_MF(:, 1, j)
                    POTOK_MF_all(:, 2) = POTOK_MF_all(:, 2) -  gl_Gran_POTOK_MF(:, 2, j)
                    POTOK_MF_all(:, 3) = POTOK_MF_all(:, 3) -  gl_Gran_POTOK_MF(:, 3, j)
                    POTOK_MF_all(:, 4) = POTOK_MF_all(:, 4) -  gl_Gran_POTOK_MF(:, 4, j)
                end if
            end do


            ! Определяем зону в которой находимся

            zone = 1


            call dev_Calc_sourse_MF(qqq, fluid1, SOURSE, zone)  ! Вычисляем источники

            if (l_1 == .TRUE.) then
                ro3 = qqq(1) - time * POTOK(1) / Volume
                Q3 = qqq(9) - time * POTOK(9) / Volume
                if (ro3 <= 0.0_8) then
					write(*, *) "Ro < 0  1490    ", qqq(1), Volume, POTOK(1), time
                end if
                u3 = (qqq(1) * qqq(2) - time * POTOK(2) / Volume + time * SOURSE(2, 1)) / ro3
                v3 = (qqq(1) * qqq(3) - time * POTOK(3) / Volume + time * SOURSE(3, 1)) / ro3
                w3 = (qqq(1) * qqq(4) - time * POTOK(4) / Volume + time * SOURSE(4, 1)) / ro3
                p3 = ((  ( qqq(5) / (ggg - 1.0) + 0.5 * qqq(1) * norm2(qqq(2:4))**2 ) &
                    - time * POTOK(5)/ Volume + time * SOURSE(5, 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) ) * (ggg - 1.0)

                if (p3 <= 0.0_8) then
                    p3 = 0.000001
                end if

                gl_Cell_par(:, gr) = (/ro3, u3, v3, w3, p3, 0.0_8, 0.0_8, 0.0_8, Q3/)

            end if

            ! Теперь посчитаем законы сохранения для остальных жидкостей

            do i = 1, 4
                if (i == 1 .and. l_1 == .FALSE.) CYCLE
                
                if (l_1 == .FALSE.) SOURSE(:, i + 1) = 0.0       ! Не перезаряжаем жидкости внутри сферы
                ro3 = fluid1(1, i) - time * POTOK_MF_all(1, i) / Volume + time * SOURSE(1, i + 1)
                if (ro3 <= 0.0_8) then
                    write(*, *) "Ro < 0  in ", ro3,  i, time, fluid1(1, i), POTOK_MF_all(1, i), Volume, SOURSE(1, i + 1)
                end if
                u3 = (fluid1(1, i) * fluid1(2, i) - time * POTOK_MF_all(2, i) / Volume + time * SOURSE(2, i + 1)) / ro3
                v3 = (fluid1(1, i) * fluid1(3, i) - time * POTOK_MF_all(3, i) / Volume + time * SOURSE(3, i + 1)) / ro3
                w3 = (fluid1(1, i) * fluid1(4, i) - time * POTOK_MF_all(4, i) / Volume + time * SOURSE(4, i + 1)) / ro3
                p3 = ((  ( fluid1(5, i) / (ggg - 1.0) + 0.5 * fluid1(1, i) * norm2(fluid1(2:4, i))**2 ) &
                    - time * POTOK_MF_all(5, i)/ Volume + time * SOURSE(5, i + 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) ) * (ggg - 1.0)
                if (p3 <= 0.0_8) then
                    p3 = 0.000001
                end if

                gl_Cell_par_MF(:, i, gr) = (/ro3, u3, v3, w3, p3/)
            end do
	
	! ----------------------------------------------------------- копируем до этого момента ----------------------
	
	end subroutine CUF_GD_3_cells_inner
	
	
	attributes(global) subroutine CUF_GD_3_cells()
	! Неподвижная сетка, газовая динамика + мультифлюид
	use MY_CUDA, gl_Gran_info => dev_gl_Gran_info, gl_Gran_neighbour => dev_gl_Gran_neighbour, &
		gl_Cell_par => dev_gl_Cell_par, gl_Cell_par_MF => dev_gl_Cell_par_MF, gl_Cell_center => dev_gl_Cell_center, &
		gl_Gran_center => dev_gl_Gran_center, norm2 => dev_norm2, ggg => dev_ggg, par_Velosity_inf => dev_par_Velosity_inf, &
		gl_Cell_dist => dev_gl_Cell_dist, gl_Gran_POTOK => dev_gl_Gran_POTOK, gl_Gran_POTOK_MF => dev_gl_Gran_POTOK_MF, &
		gl_Gran_square => dev_gl_Gran_square, gl_Gran_normal => dev_gl_Gran_normal, &
		gl_Cell_info => dev_gl_Cell_info, gl_Cell_Volume => dev_gl_Cell_Volume, gl_Cell_gran => dev_gl_Cell_gran, &
		par_R0 => dev_par_R0
	!use cudadevice
	!use libm
	implicit none
	
	integer(4) :: gr
	
    integer(4) :: st, s1, s2, i, j, k, zone
    real(8) :: qqq1(9), qqq2(9), qqq(9)  ! Переменные в ячейке
    real(8) :: fluid1(5, 4), fluid2(5, 4)
    real(8) :: dist, dsl, dsc, dsp
    real(8) :: POTOK(9), ttest(3)
    real(8) :: POTOK_MF(5)
    real(8) :: POTOK_MF_all(5, 4)
    real(8) :: time, Volume, TT, U8, rad1, rad2, aa, bb, cc
    real(8) :: SOURSE(5,5)  ! Источники массы, импульса и энергии для плазмы и каждого сорта мультифлюида
    real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3
	logical :: l_1
	
	gr = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! Номер потока
	time = time_step
	
	if (gr > dev_Ncell) return
	
	! ----------------------------------------------------------- можно просто скопировать код из хоста ----------------------
	
	if(gl_Cell_info(gr) == 0) return
	
    l_1 = .TRUE.
    !if (norm2(gl_Cell_center(:, gr)) <= 1.0 * par_R0)   write(*,*) "NONE ERROR  27465678uhgfdr"  !l_1 = .FALSE.    ! Не считаем внутри сферы   здесь такого не должно быть
    POTOK = 0.0
    SOURSE = 0.0
    POTOK_MF_all = 0.0
    Volume = gl_Cell_Volume(gr)
    qqq = gl_Cell_par(:, gr)
    fluid1 = gl_Cell_par_MF(:, :, gr)
    ! Просуммируем потоки через грани
    do i = 1, 6
        j = gl_Cell_gran(i, gr)
        if (j == 0) CYCLE
        if (gl_Gran_neighbour(1, j) == gr) then
            POTOK = POTOK + gl_Gran_POTOK(:, j)
            POTOK_MF_all(:, 1) = POTOK_MF_all(:, 1) +  gl_Gran_POTOK_MF(:, 1, j)
            POTOK_MF_all(:, 2) = POTOK_MF_all(:, 2) +  gl_Gran_POTOK_MF(:, 2, j)
            POTOK_MF_all(:, 3) = POTOK_MF_all(:, 3) +  gl_Gran_POTOK_MF(:, 3, j)
            POTOK_MF_all(:, 4) = POTOK_MF_all(:, 4) +  gl_Gran_POTOK_MF(:, 4, j)
        else
            POTOK = POTOK - gl_Gran_POTOK(:, j)
            POTOK_MF_all(:, 1) = POTOK_MF_all(:, 1) -  gl_Gran_POTOK_MF(:, 1, j)
            POTOK_MF_all(:, 2) = POTOK_MF_all(:, 2) -  gl_Gran_POTOK_MF(:, 2, j)
            POTOK_MF_all(:, 3) = POTOK_MF_all(:, 3) -  gl_Gran_POTOK_MF(:, 3, j)
            POTOK_MF_all(:, 4) = POTOK_MF_all(:, 4) -  gl_Gran_POTOK_MF(:, 4, j)
        end if
    end do
 
 
    ! Определяем зону в которой находимся
    if(qqq(9)/qqq(1) < 50.0) then
 
        if(norm2(qqq(2:4))/sqrt(ggg*qqq(5)/qqq(1)) > 3.0) then
            zone = 1
        else
            zone = 2
        end if
    else
        if(norm2(qqq(2:4))/sqrt(ggg*qqq(5)/qqq(1)) > 1.1) then
            zone = 4
        else
            zone = 3
        end if
    end if
 
 
    call dev_Calc_sourse_MF(qqq, fluid1, SOURSE, zone)  ! Вычисляем источники
 
 
 
    if (l_1 == .TRUE.) then
        ro3 = qqq(1) - time * POTOK(1) / Volume
        Q3 = qqq(9) - time * POTOK(9) / Volume
        if (ro3 <= 0.0_8) then
            write(*,*) "Ro < 0  1490  ", qqq(1), Volume
        end if
        u3 = (qqq(1) * qqq(2) - time * POTOK(2) / Volume + time * SOURSE(2, 1)) / ro3
        v3 = (qqq(1) * qqq(3) - time * POTOK(3) / Volume + time * SOURSE(3, 1)) / ro3
        w3 = (qqq(1) * qqq(4) - time * POTOK(4) / Volume + time * SOURSE(4, 1)) / ro3
        p3 = ((  ( qqq(5) / (ggg - 1.0) + 0.5 * qqq(1) * norm2(qqq(2:4))**2 ) &
            - time * POTOK(5)/ Volume + time * SOURSE(5, 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) ) * (ggg - 1.0)
        
        if (p3 <= 0.0_8) then
            p3 = 0.000001
        end if
    
        gl_Cell_par(:, gr) = (/ro3, u3, v3, w3, p3, 0.0_8, 0.0_8, 0.0_8, Q3/)
    
    end if
 
    ! Теперь посчитаем законы сохранения для остальных жидкостей
 
    do i = 1, 4
        if (i == 1 .and. l_1 == .FALSE.) CYCLE       ! Пропускаем внутреннюю сферу для сорта 1
        if (l_1 == .FALSE.) SOURSE(:, i + 1) = 0.0       ! Пропускаем внутреннюю сферу для сорта 1
        ro3 = fluid1(1, i) - time * POTOK_MF_all(1, i) / Volume + time * SOURSE(1, i + 1)
        if (ro3 <= 0.0_8) then
            print*, "Ro < 0  in "
        end if
        u3 = (fluid1(1, i) * fluid1(2, i) - time * POTOK_MF_all(2, i) / Volume + time * SOURSE(2, i + 1)) / ro3
        v3 = (fluid1(1, i) * fluid1(3, i) - time * POTOK_MF_all(3, i) / Volume + time * SOURSE(3, i + 1)) / ro3
        w3 = (fluid1(1, i) * fluid1(4, i) - time * POTOK_MF_all(4, i) / Volume + time * SOURSE(4, i + 1)) / ro3
        p3 = ((  ( fluid1(5, i) / (ggg - 1.0) + 0.5 * fluid1(1, i) * norm2(fluid1(2:4, i))**2 ) &
            - time * POTOK_MF_all(5, i)/ Volume + time * SOURSE(5, i + 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) ) * (ggg - 1.0)
        if (p3 <= 0.0_8) then
            p3 = 0.000001
        end if
 
        gl_Cell_par_MF(:, i, gr) = (/ro3, u3, v3, w3, p3/)
	end do
	
	! ----------------------------------------------------------- копируем до этого момента ----------------------
	
	end subroutine CUF_GD_3_cells
	
	
	subroutine CUDA_info()
	! Функция, печатающая информацию о видеокарте на экран (только о первой видеокарте)
	    use OMP_lib
         use cudafor
        
        implicit none
        
         type (cudaDeviceProp) :: prop
        integer :: nDevices = 0, i, ierr
        
         ierr = cudaGetDeviceCount(nDevices)
		
		 print*, "  "
		 print*, "  "
		print*, "CUDA_info()"
		print*, "  "
        
        write(*, "('Chislo videokard: ', i0)")  nDevices
		
		if (nDevices == 0) then
			STOP "NET CUDA!!!"
		end if
		
		 ierr = cudaGetDeviceProperties(prop, 0)
		
		 write(*, "('Version: ', i0,'.',i0)") prop%major, prop%minor
		 write(*, "('Chislo mulitiprocessorov: ', i0)") prop%multiProcessorCount
		 write(*, "('Max chislo potokov/multiprocessor: ', i0)") prop%maxThreadsPerMultiprocessor
		 write(*, "('Global memory (GB): ', f9.3)") prop%totalGlobalMem / 1024.0**3
		
		 print*, "Konfiguration of working"
		 write(*, "('max grid size: ', 2(i0, ' x '), i0 )") prop%maxGridSize
		 write(*, "('max block size: ', 2(i0, ' x '), i0 )") prop%maxThreadsDim
		 write(*, "('Max chislo potokov v blocke: ', i0)") prop%maxThreadsPerBlock
		print*, "                "
		print*, "                "
	end subroutine CUDA_info

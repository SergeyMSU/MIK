
    
    module GEO_PARAM                     ! Модуль геометрических - сеточных параметров - констант программы
    
    implicit none
    
    logical, parameter :: par_developer_info = .True.   ! Параметр разработчика для вывода информационных сообщений
    real(8), parameter :: par_R_character = 35.6505         ! Характерный размер в задаче (расстояние до TS на начальном этапе построения сетки)
    real(8), parameter :: par_koeff_HP = 1.3            ! На сколько умножить характерный размер для получения HP
    real(8), parameter :: par_koeff_BS = 2.0            ! На сколько умножить характерный размер для получения BS
    real(8), parameter :: par_R0 = 0.197035         ! Характерный размер 1 а.е. (внутренней сферы) Там находится вторая точка на лучах от цетра (первая находится в нуле)
    ! параметр нужен для геометрической точности сетки eps = par_R_character/10000
    real(8) :: par_R_END = 300.0         !  
    real(8) :: par_R_LEFT = -390.0         !  Левая граница
    real(8), parameter :: par_pi_8 = acos(-1.0_8)         
    real(8), parameter :: par_pi_4 = acos(-1.0_4)        
    real(8), parameter :: cpi4 = 12.56637061435917295384
    real(8), parameter :: ggg = (5.0/3.0)
    real(8), parameter :: par_kk = 10000.0                ! Масштаб характерного размера регулируется
	real(8) :: par_R_inner = 5.0_8     ! До какого расстояния внутренняя сфера
	
	
	
	real(8), parameter :: par_nat_TS = 0.002_8 ! 0.003_8 !0.0000001_8 !0.003_8                ! Коэффициент натяжения ударной волны  0.002
	real(8), parameter :: par_nat_HP = 1.5_8 ! 1.5                 ! Коэффициент натяжения контакта  0.0001
	real(8), parameter :: par_nat_BS = 0.00003_8                ! Коэффициент натяжения внешней ударной волны 0.0002
	
	real(8), parameter :: koef1 = 0.4_8! 0.3      ! Коэффицинт запаздывания скорости ударной волны
    real(8), parameter :: koef2 = 0.6_8 ! 0.5
    real(8), parameter :: koef3 = 0.7_8   ! 0.3
	
    
    real(8), parameter :: par_a_2 = 0.130738_8        ! Параметр в сечении перезарядки
    real(8), parameter :: par_n_p_LISM = 3.5_8         ! в перезарядке
    real(8), parameter :: par_n_H_LISM_ = 1.0_8
    real(8), parameter :: par_Kn = 49.9018   !0.4326569808         ! в перезарядке
    
    real(8), parameter :: par_chi_real = 1.0_8      ! С каким хи считаем реально
    real(8), parameter :: par_chi = 41.6479_8      ! С каким хи надо было бы считать
    real(8), parameter :: par_Velosity_inf = -2.54279_8
	real(8), parameter :: par_Mach_alf = 12.8816_8
	real(8), parameter :: par_Mach_0 = 6.44_8
	real(8), parameter :: par_B_inf = 13.9666_8
	real(8), parameter :: par_alphaB_inf = 1.04719755_8   ! 60 градусов
	real(8), parameter :: par_k_Br = 0.00197035_8
    
    integer, parameter :: par_R_int = 70  ! Сколько а.е. не считаем внутри
    
    real(8) :: par_kk1 = 1.7_8     ! Степень сгущения сетки к нулю в области до TS: 1 - линейное, 2 - квадратичное и т.д.
	real(8) :: par_kk12 = 1.7_8     ! Степень сгущения сетки к нулю в области до TS но после внутренней сферы: 1 - линейное, 2 - квадратичное и т.д.
    real(8) :: par_kk2 = 2.0_8     ! Степень сгущения в головной области на бесконечности
    real(8) :: par_kk3 = 1.8_8     ! Степень сгущения в хвосте
	real(8) :: par_kk31 = 1.1_8     ! Степень сгущения в хвосте для точек на контакте (первая точка в О - луче)
	

    integer :: par_l_phi = 48 !4   ! Количество вращений двумерной сетки для получения трёхмерной (разрешение по углу phi)
                                ! Должно делиться на 4 для удобного вывода результатов в плоскостях
    
    integer(4) :: par_m_A = 30      ! Количество лучей A в плоскости
    integer(4) :: par_m_BC = 18      ! Количество лучей B/C в плоскости
    integer(4) :: par_m_O = 17      ! Количество лучей O в плоскости
    integer(4) :: par_m_K = 7      ! Количество лучей K в плоскости
    !integer :: par_m_A = 5      ! Количество лучей A в плоскости
    !integer :: par_m_BC = 5      ! Количество лучей B/C в плоскости
    !integer :: par_m_O = 5      ! Количество лучей O в плоскости
    !integer :: par_m_K = 5      ! Количество лучей K в плоскости
    real(8) :: par_triple_point = 13.0 * par_pi_4/40.0     ! До какого угла начиная от pi/2 (с положительного x) тройная точка
    
    ! Количество точек по лучам A
    integer(4) :: par_n_TS =  26! 3                 ! Количество точек до TS (TS включается)
    integer(4) :: par_n_HP =  40! 4                 ! Количество точек HP (HP включается)  всё от 0 считается
    integer(4) :: par_n_BS =  60! 5                 ! Количество точек BS (BS включается)
    integer(4) :: par_n_END = 72! 6                ! Количество точек до конца сетки (конец включается)
    integer(4) :: par_n_IA =  12                   ! Количество точек, которые входят во внутреннюю область
	integer(4) :: par_n_IB =  14                   ! Количество точек, которые входят во внутреннюю область (с зазором)
    
    integer :: par_n_points  ! Всего точек в сетке
	
	integer, parameter :: par_tet(3, 8) = reshape( (/ 1, 3, 5, 2, 3, 5, 1, 3, 6, 2, 3, 6, 1, 4, 5, 2, 4, 5, 1, 4, 6, 2, 4, 6/), (/3, 8/) )
    
    end module GEO_PARAM
    
    
    module STORAGE                       ! Модуль глобальных данных и типов (все переменные начинаются на gl - global)
    implicit none

    ! Если нужно поменять точность, то это нужно ещё сделать в функции выделения памяти и инициализации, но лучше не менять
    
    real(8), allocatable :: gl_x(:)   ! набор x-координат узлов сетки 
    real(8), allocatable :: gl_y(:)   ! набор y-координат узлов сетки 
    real(8), allocatable :: gl_z(:)   ! набор z-координат узлов сетки 
	
    
    ! Параметр MOVE - означает, что эти массивы используются (и инициализируются) только в случае движения сетки
    ! Скорости движения узлов
    real(8), allocatable :: gl_Vx(:)   ! набор z-координат узлов сетки    !MOVE
    real(8), allocatable :: gl_Vy(:)   ! набор z-координат узлов сетки    !MOVE
    real(8), allocatable :: gl_Vz(:)   ! набор z-координат узлов сетки    !MOVE
    integer(4), allocatable :: gl_Point_num(:)   ! Сколько граней записали свою скорость движения в данный узел      !MOVE    
    
    ! Координаты узлов (сдвоенные для случая движения сетки)
    real(8), allocatable :: gl_x2(:, :)   ! (:, 2) набор x-координат узлов сетки     !MOVE
    real(8), allocatable :: gl_y2(:, :)   ! (:, 2) набор y-координат узлов сетки     !MOVE
    real(8), allocatable :: gl_z2(:, :)   ! (:, 2) набор z-координат узлов сетки     !MOVE
    
    ! Лучи, на которых распологаются точки сетки - так 
    integer(4), allocatable :: gl_RAY_A(:,:,:)   ! Набор А-лучей размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
    integer(4), allocatable :: gl_RAY_B(:,:,:)   ! Набор B-лучей размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
    integer(4), allocatable :: gl_RAY_C(:,:,:)   ! Набор C-лучей размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
    integer(4), allocatable :: gl_RAY_O(:,:,:)   ! Набор O-лучей размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
    integer(4), allocatable :: gl_RAY_K(:,:,:)   ! Набор K-лучей размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
    integer(4), allocatable :: gl_RAY_D(:,:,:)   ! Набор D-лучей размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
    integer(4), allocatable :: gl_RAY_E(:,:,:)   ! Набор E-лучей размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
    
    ! Ячейки
    integer(4), allocatable :: gl_Cell_A(:,:,:)   ! Набор A-ечеек размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
    integer(4), allocatable :: gl_Cell_B(:,:,:)   ! Набор B-ечеек размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
    integer(4), allocatable :: gl_Cell_C(:,:,:)   ! Набор C-ечеек размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
    ! Предыдущие наборы нужны для большей структуризации сетки, удобство доступа и т.д.
    ! Ячейки разделены на группы
    
    integer(4), allocatable :: gl_all_Cell(:,:)   ! Весь набор ячеек (8, :) - первая координата массива - это набор узлов ячейки
    
	
    integer(4), allocatable :: gl_all_Cell_inner(:)   ! номера ячеек, которые находятся внутри (считаются отдельно)
    
    integer(4), allocatable :: gl_Cell_neighbour(:,:)   ! (6, :) Набор из 6 соседей для каждой ячейки 
	! (если номер соседа = 0, то его нет в этом направлении)
    ! -1   ! Граница (набегающий поток)
	! -3   ! Граница верхний цилиндр
	!  -2  ! Выходная граница
	
	integer(4), allocatable :: gl_Cell_gran(:,:)        ! (6, :) Набор из 6 граней для каждой ячейки (если номер = 0, то грани нет в этом направлении)
    integer(4), allocatable :: gl_Cell_info(:)        ! 
    
    real(8), allocatable :: gl_Cell_Volume(:)           ! Набор объёмов ячеек
    real(8), allocatable :: gl_Cell_dist(:)             ! Минимальное расстояние до грани в каждой ячейки  4444444444444444
    real(8), allocatable :: gl_Cell_center(:, :)             ! (3, :) Центр каждой ячейки  4444444444444444
    real(8), allocatable :: gl_Cell_par(:, :)           ! (9, :) Набор параметров (8 стартовых + Q)
    real(8), allocatable :: gl_Cell_par_MF(:,:,:)           ! Набор параметров (5, 4,:)  Мультифлюид параметры (по 5 для каждой из 4-х жидкостей)
    character, allocatable :: gl_Cell_type(:)           ! Тип каждой ячейки А, Б, С
    integer(4), allocatable :: gl_Cell_number(:, :)     ! (3, :) номер каждой ячейки внутри своего типа
    
    ! Набор переменных для подвижной сетки (они имеют второй индекс, обознающий шаг по времни, следующий или предыдущий)
    real(8), allocatable :: gl_Cell_Volume2(:, :)           ! (всего ячеек, 2) Набор объёмов ячеек   !MOVE
    real(8), allocatable :: gl_Cell_center2(:, :, :)             ! (3, :) Центр каждой ячейки  4444444444444444  !MOVE
    
    
    ! Блок граней ----------------- нужны для расчёта площадей, потоков, объёмов, нормалей и т.д.
    integer(4), allocatable :: gl_all_Gran(:,:)       ! Все грани (4,:) имеют по 4 узла
    integer(4), allocatable :: gl_Gran_neighbour(:,:) ! Соседи каждой грани (2,:) имеют по 2 соседа, нормаль ведёт от первого ко второму
	! -1  -2  -3  бывает
	! -1 - входная сфера
	! -2 - выходная заняя стенка
	! -3  -  верхний цилиндр
	
	integer(4), allocatable :: gl_Gran_neighbour_TVD(:,:) ! TVD-Соседи каждой грани (2,:) имеют по 2 соседа
	! 0 - значит соседа нет
	! При этом первый TVD-сосед - это сосед первого обычного соседа. Т.е. нормаль грани тоже ведёт от первого ко второму
    
    real(8), allocatable :: gl_Gran_normal(:,:)       ! (3, :) Нормаль грани   4444444444444444
    real(8), allocatable :: gl_Gran_square(:)         ! (:) Площадь грани
    real(8), allocatable :: gl_Gran_POTOK(:,:)       ! (10, :) поток грани    последний - дивергенция магнитного поля для очистки
    real(8), allocatable :: gl_Gran_POTOK_MF(:,:,:)       ! (5, 4, :) поток грани мультифлюидных жидкостей
    real(8), allocatable :: gl_Gran_center(:,:)       ! (3, :) 
	integer(4), allocatable :: gl_Gran_type(:)      ! Показывает тип грани (0 - обычная, 1 - TS, 2 - HP)
	integer(4), allocatable :: gl_Gran_scheme(:)      ! Показывает тип грани (0 - обычная, 1 - TS, 2 - HP)
    integer(4), allocatable :: gl_Gran_info(:)         ! (:) классификатор грани (см. схему) 
	! для расчёта какие грани пропускать для внутренней сферы (бредово сделано, надо переделывать)
	
    integer(4), allocatable :: gl_all_Gran_inner(:)   ! номера граней, которые находятся внутри (считаются отдельно)
    
    ! Набор переменных для подвижной сетки (они имеют второй индекс, обознающий шаг по времни, следующий или предыдущий)
    real(8), allocatable :: gl_Gran_normal2(:, :, :)       ! (3, :, 2) Нормаль грани                       !MOVE
    real(8), allocatable :: gl_Gran_center2(:, :, :)       ! (3, :, 2)                          !Move
    real(8), allocatable :: gl_Gran_square2(:, :)         ! (:, 2) Площадь грани          !MOVE
    
    ! Поверхности выделения
    integer(4), allocatable :: gl_Contact(:)       ! Контакт (не входит в начальное построение сетки, ищется в отдельной функции Find_Surface
    ! Контакт состоит из номеров граней
    integer(4), allocatable :: gl_TS(:)
    integer(4), allocatable :: gl_BS(:)
	
	
	! Набор массивов для программы интерполяции (аналогичных основной сетке)
	integer(4), allocatable :: gl_Cell_gran_inter(:,:)  ! (6, :)
	real(8), allocatable :: gl_Gran_normal_inter(:,:)   ! (3, :)
	real(8), allocatable :: gl_Gran_center_inter(:,:)  ! (3, :)
	integer(4), allocatable :: gl_Gran_neighbour_inter(:,:)   ! (2, :)
	real(8), allocatable :: gl_Cell_center_inter(:, :)        ! (3, :)
	integer(4), allocatable :: gl_Cell_neighbour_inter(:,:)   ! (6, :)
	real(8), allocatable :: gl_Cell_par_inter(:, :)			  ! (9, :)
    real(8), allocatable :: gl_Cell_par_MF_inter(:,:,:)		  ! (5, 4,:) 
    
    end module STORAGE
    
    
    subroutine Set_STORAGE()             ! Функция выделяет память для всех элементов модуля STORAGE, используя параметры из модуля GEO_PARAM
    
    use STORAGE
    use GEO_PARAM
    implicit none
    
    integer :: n1, n2
    
    if (par_developer_info) print *, "START Set_STORAGE"
    ! Выделяем память под переменные
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
    
    allocate( gl_Cell_par(9, size(gl_Cell_A(:,:,:)) + size(gl_Cell_B(:,:,:)) + size(gl_Cell_C(:,:,:)) ) )
    allocate(gl_Cell_par_MF(5, 4, size(gl_Cell_A(:,:,:)) + size(gl_Cell_B(:,:,:)) + size(gl_Cell_C(:,:,:))))
    
    ! Посчитаем число узлов в сетке
    par_n_points = par_n_END * par_l_phi * (par_m_A + par_m_BC) + par_m_K * (par_n_TS + par_m_O) * par_l_phi + par_l_phi * (par_n_END - par_n_TS + 1) * par_m_O - &
            (par_m_A + par_m_BC + par_m_K - 1) * par_l_phi - par_n_END * (par_l_phi - 1) - (par_n_TS + par_m_O - 1) * (par_l_phi - 1)  ! Всего точек в сетке
    
    allocate(gl_x(par_n_points))
    allocate(gl_y(par_n_points))
    allocate(gl_z(par_n_points))
    !gl_x = [real(8) ::]
    !gl_y = [real(8) ::]
    !gl_z = [real(8) ::]
    
    n2 =  (par_n_END - 1) * (par_m_A + par_m_BC - 1) + (par_n_TS - 1 + par_m_O) * (par_m_K) + &
        (par_n_END - par_n_TS) * par_m_O ! Число ячеек в 1 слое по углу
    
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
    
    allocate(gl_Contact( (par_m_O + par_m_A + par_m_BC -1) * par_l_phi ))   ! Выделяем память под контакт
    allocate(gl_TS( (par_m_A + par_m_BC + par_m_K -1) * par_l_phi ))   ! Выделяем память под TS
    allocate(gl_BS( (par_m_A - 1) * par_l_phi ))   ! Выделяем память под BS
    
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
    
    
    ! Заполняем начальными значениями
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
    gl_Cell_info = 2
	gl_Gran_type = 0
	gl_Gran_scheme = 2
    
    
    gl_Cell_A = -1
    gl_Cell_B = -1
    gl_Cell_C = -1
    
    gl_all_Cell = -1
    gl_Cell_neighbour = 0          ! 0 - нет соседа (при этом это НЕ граница - у границы свой номер)
    gl_Cell_gran = 0
    
    gl_all_Gran = 0
    gl_Gran_neighbour = 0
	gl_Gran_neighbour_TVD = 0
    
    gl_Contact = 0
    
    gl_Gran_normal = 0.0;
    gl_Gran_square = 0.0;
    
    if (par_developer_info) print *, "END Set_STORAGE"
    
    end subroutine Set_STORAGE

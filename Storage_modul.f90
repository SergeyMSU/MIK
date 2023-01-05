
    
    module GEO_PARAM                     ! Модуль геометрических - сеточных параметров - констант программы
    
    implicit none
    
    logical, parameter :: par_developer_info = .True.   ! Параметр разработчика для вывода информационных сообщений
    real(8), parameter :: par_R_character = 35.6505         ! Характерный размер в задаче (расстояние до TS на начальном этапе построения сетки)
    real(8), parameter :: par_koeff_HP = 1.3            ! На сколько умножить характерный размер для получения HP
    real(8), parameter :: par_koeff_BS = 2.0            ! На сколько умножить характерный размер для получения BS
    real(8), parameter :: par_R0 = 0.256505         ! Характерный размер 1 а.е. (внутренней сферы)
    ! параметр нужен для геометрической точности сетки eps = par_R_character/10000
    real(8), parameter :: par_R_END = 300.0         ! 
    real(8), parameter :: par_R_LEFT = -390.0         ! 
    real(8), parameter :: par_pi_8 = acos(-1.0_8)         ! Характерный размер в задаче (примерное расстояние до TS)
    real(8), parameter :: par_pi_4 = acos(-1.0_4)         ! Характерный размер в задаче (примерное расстояние до TS)
    real(8), parameter :: cpi4 = 12.56637061435917295384
    real(8), parameter :: ggg = (5.0/3.0)
    real(8), parameter :: par_kk = 10000.0                ! Масштаб характерного размера регулируется
    
    real(8), parameter :: par_a_2 = 0.1307345665         ! Параметр в сечении перезарядки
    real(8), parameter :: par_n_p_LISM = 1.0_8         ! в перезарядке
    real(8), parameter :: par_n_H_LISM_ = 3.0_8
    real(8), parameter :: par_Kn = 43.3   !0.4326569808         ! в перезарядке
    
    real(8), parameter :: par_chi_real = 1.0_8      ! С каким хи считаем реально
    real(8), parameter :: par_chi = 36.1275_8      ! С каким хи надо было бы считать
    real(8), parameter :: par_Velosity_inf = -2.54338_8
    
    integer, parameter :: par_R_int = 70  ! Сколько а.е. не считаем внутри
    
    real(8), parameter :: par_kk1 = 1.7_8     ! Степень сгущения сетки к нулю в области до TS: 1 - линейное, 2 - квадратичное и т.д.
    real(8), parameter :: par_kk2 = 2.0_8     ! Степень сгущения в головной области на бесконечности
    real(8), parameter :: par_kk3 = 1.8_8     ! Степень сгущения в хвосте

    integer :: par_l_phi = 48 !4   ! Количество вращений двумерной сетки для получения трёхмерной (разрешение по углу phi)
                                ! Должно делиться на 4 для удобного вывода результатов в плоскостях
    integer :: par_m_A = 30 !5     ! Количество лучей A в плоскости
    integer :: par_m_BC = 18 !3     ! Количество лучей B/C в плоскости
    integer :: par_m_O = 17 !7 !3     ! Количество лучей O в плоскости
    integer :: par_m_K = 7 !2     ! Количество лучей K в плоскости
    real(8) :: par_triple_point = 13.0 * par_pi_4/40.0     ! До какого угла начиная от pi/2 (с положительного x) тройная точка
    
    ! Количество точек по лучам A
    integer :: par_n_TS =  26! 3                 ! Количество точек до TS (TS включается)
    integer :: par_n_HP =  40! 4                 ! Количество точек HP (HP включается)  всё от 0 считается
    integer :: par_n_BS =  60! 5                 ! Количество точек BS (BS включается)
    integer :: par_n_END = 72! 6                ! Количество точек до конца сетки (конец включается)
    
    integer :: par_n_points  ! Всего точек в сетке
    
    end module GEO_PARAM
    
    
    module STORAGE                       ! Модуль глобальных данных и типов (все переменные начинаются на gl - global)
    implicit none
    
    real(8), allocatable :: gl_x(:)   ! набор x-координат узлов сетки 
    real(8), allocatable :: gl_y(:)   ! набор y-координат узлов сетки 
    real(8), allocatable :: gl_z(:)   ! набор z-координат узлов сетки 
    ! Если нужно поменять точность, то это нужно ещё сделать в функции выделения памяти и инициализации, лучше не менять
    
    ! Лучи, на которых распологаются точки сетки - так 
    integer(4), allocatable :: gl_RAY_A(:,:,:)   ! Набор А-лучей размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
    integer(4), allocatable :: gl_RAY_B(:,:,:)   ! Набор B-лучей размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
    integer(4), allocatable :: gl_RAY_C(:,:,:)   ! Набор C-лучей размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
    integer(4), allocatable :: gl_RAY_O(:,:,:)   ! Набор O-лучей размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
    integer(4), allocatable :: gl_RAY_K(:,:,:)   ! Набор K-лучей размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
    integer(4), allocatable :: gl_RAY_D(:,:,:)   ! Набор D-лучей размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
    integer(4), allocatable :: gl_RAY_E(:,:,:)   ! Набор E-лучей размерности 3 (на этом луче, в этой плоскости, по углу в пространстве)
    
    ! Ячейки
    integer(4), allocatable :: gl_Cell_A(:,:,:)   ! Набор A-ечеек размерности 4 (на этом луче, в этой плоскости, по углу в пространстве)
    integer(4), allocatable :: gl_Cell_B(:,:,:)   ! Набор B-ечеек размерности 4 (на этом луче, в этой плоскости, по углу в пространстве)
    integer(4), allocatable :: gl_Cell_C(:,:,:)   ! Набор C-ечеек размерности 4 (на этом луче, в этой плоскости, по углу в пространстве)
    ! Предыдущие наборы нужны для большей структуризации сетки, удобство доступа и т.д.
    ! Ячейки разделены на группы
    
    integer(4), allocatable :: gl_all_Cell(:,:)   ! Весь набор ячеек (8, :) - первая координата массива - это набор узлов ячейки
    
    integer(4), allocatable :: gl_all_Cell_inner(:)   ! номера ячеек, которые находятся внутри (считаются отдельно)
    
    integer(4), allocatable :: gl_Cell_neighbour(:,:)   ! Набор из 6 соседей для каждой ячейки (если номер соседа = 0, то его нет в этом направлении)
    integer(4), allocatable :: gl_Cell_gran(:,:)        ! Набор из 6 граней для каждой ячейки (если номер = 0, то грани нет в этом направлении)
    integer(4), allocatable :: gl_Cell_info(:)        ! Набор из 6 граней для каждой ячейки (если номер = 0, то грани нет в этом направлении)
    
    real(8), allocatable :: gl_Cell_Volume(:)           ! Набор объёмов ячеек
    real(8), allocatable :: gl_Cell_dist(:)             ! Минимальное расстояние до грани в каждой ячейки  4444444444444444
    real(8), allocatable :: gl_Cell_center(:, :)             ! (3, :) Центр каждой ячейки  4444444444444444
    real(8), allocatable :: gl_Cell_par(:, :)           ! Набор параметров (8 стартовых)
    real(8), allocatable :: gl_Cell_par_MF(:,:,:)           ! Набор параметров (5,4,:)  Мультифлюид параметры (по 5 для каждой из 4-х жидкостей)
    
    
    ! Блок граней ----------------- нужны для расчёта площадей, потоков, объёмов, нормалей и т.д.
    integer(4), allocatable :: gl_all_Gran(:,:)       ! Все грани (4,:) имеют по 4 узла
    integer(4), allocatable :: gl_Gran_neighbour(:,:) ! Соседи каждой грани (2,:) имеют по 2 соседа, нормаль ведёт от первого ко второму
    
    real(8), allocatable :: gl_Gran_normal(:,:)       ! (3, :) Нормаль грани   4444444444444444
    real(8), allocatable :: gl_Gran_square(:)         ! (:) Площадь грани
    real(8), allocatable :: gl_Gran_POTOK(:,:)       ! (8, :) поток грани
    real(8), allocatable :: gl_Gran_POTOK_MF(:,:,:)       ! (5,4, :) поток грани мультифлюидных жидкостей
    real(8), allocatable :: gl_Gran_center(:,:)       ! (3, :) поток грани
    integer(4), allocatable :: gl_Gran_info(:)         ! (:) Площадь грани
    integer(4), allocatable :: gl_all_Gran_inner(:)   ! номера ячеек, которые находятся внутри (считаются отдельно)
    
    ! Поверхности выделения
    integer(4), allocatable :: gl_Contact(:)       ! Контакт (не входит в начальное построение сетки, ищется в отдельной функции Find_Surface
    
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
    allocate(gl_Gran_normal(3, n1))
    allocate(gl_Gran_square(n1))
    allocate(gl_Gran_info(n1))
    allocate(gl_Gran_POTOK(9, n1))
    allocate(gl_Gran_POTOK_MF(5, 4, n1))
    allocate(gl_Gran_center(3, n1))
    
    allocate(gl_Contact( (par_m_O + par_m_A + par_m_BC -1) * par_l_phi ))   ! Выделяем память под контакт
    
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
    gl_Cell_info = 0
    
    
    gl_Cell_A = -1
    gl_Cell_B = -1
    gl_Cell_C = -1
    
    gl_all_Cell = -1
    gl_Cell_neighbour = 0          ! 0 - нет соседа (при этом это НЕ граница - у границы свой номер)
    gl_Cell_gran = 0
    
    gl_all_Gran = 0
    gl_Gran_neighbour = 0
    
    gl_Contact = 0
    
    gl_Gran_normal = 0.0;
    gl_Gran_square = 0.0;
    
    if (par_developer_info) print *, "END Set_STORAGE"
    
    end subroutine Set_STORAGE

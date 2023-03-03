    !  MIK.f90

    !****************************************************************************
    !
    !  PROGRAM: MIK model - Malama & Izmodenov & Korolkov model
    !
    !****************************************************************************
    ! Описание модулей

    include "Storage_modul.f90"
    include "Solvers.f90"
    include "Help_func.f90"
    include "Move_func.f90"

    ! ceiling(a) возвращает наименьшее целое число, большее или равное a. Тип - integer по умолчанию

    module My_func                     ! Модуль интерфейсов для внешних функций

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

    ! Вспомогательные функции

    real(8) pure function polar_angle(x, y)
    use GEO_PARAM
    implicit none
    real(8), intent(in) :: x, y

    if (dabs(x) + dabs(y) < 0.00001 / par_R_character) then
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
    ! Блок функций для начального построения сетки

    subroutine Build_Mesh_start()        ! Начальное построение сетки

    use STORAGE
    use GEO_PARAM
    implicit none

    integer(4), automatic :: i, j, k, N1, N2, N3, i1, kk, node, kk2, ni
    real(8), automatic :: r, phi, the, xx, x, y, z, rr, x2, y2, z2

    ! Сетка строится в декартовой системе координат, хотя сама сетка цилиндрически симметричная
    ! Нужно уметь преобразовывать x,y,z <---> x, r, phi  (x - ось симметрии), а также сферические координаты

    node = 1
    ! Задаём первый узел
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

    ! **************************************************************** Начнём с узлов сетки

    ! Цикл генерации точек на лучах А и их связывание с этими лучами ************************************************************
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

                ! Вычисляем координаты текущего луча в пространстве
                the = (j - 1) * par_pi_8/2.0/(N2 - 1)
                phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
                ! Вычисляем координаты точки на луче

                ! до TS
                if (i <= par_n_TS) then  ! До расстояния = par_R_character
                    r =  par_R0 + (par_R_character - par_R0) * (DBLE(i)/par_n_TS)**par_kk1
                    !print *, r
                    !pause
                else if (i <= par_n_HP) then  ! До расстояния = par_R_character * 1.3
                    r = par_R_character + (i - par_n_TS) * 0.3 * par_R_character/(par_n_HP - par_n_TS)
                else if (i <= par_n_BS) then  ! До расстояния = par_R_character * 2
                    r = 1.3 * par_R_character + (i - par_n_HP) * 0.7 * par_R_character/(par_n_BS - par_n_HP)
                else  ! До расстояния = par_R_character * 1.3_8
                    !r = 2.0 * par_R_character + (i - par_n_BS) * (par_R_END - 2.0 * par_R_character)/(par_n_END - par_n_BS)
                    r = 2.0 * par_R_character + (par_R_END - 2.0 * par_R_character) * (DBLE(i- par_n_BS)/(par_n_END - par_n_BS))**par_kk2
                end if

                ! Создаём точку
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

    ! Цикл генерации точек на лучах B и их связывание с этими лучами ************************************************************
    do k = 1, N3
        do j = 1, N2
            do i = 1, N1

                if (i == 1) then
                    gl_RAY_B(i, j, k) = 1
                    CYCLE
                end if

                ! Вычисляем координаты текущего луча в пространстве
                the = par_pi_8/2.0 + (j) * par_triple_point/(N2)
                phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
                ! Вычисляем координаты точки на луче

                ! до TS
                if (i <= par_n_TS) then  ! До расстояния = par_R_character
                    r =  par_R0 + (par_R_character - par_R0) * (REAL(i, KIND = 4)/par_n_TS)**par_kk1
                    !print *, r
                    !pause
                else if (i <= par_n_HP) then  ! До расстояния = par_R_character * 1.3
                    xx = 1.3 * par_R_character * (1.0_8 - cos(the - par_pi_8/2.0_8))/cos(the - par_pi_8/2.0_8)
                    r = par_R_character + (i - par_n_TS) * (0.3 * par_R_character + xx) /(par_n_HP - par_n_TS)
                end if

                ! Создаём точку
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

    ! Цикл генерации точек на лучах C и их связывание с этими лучами ************************************************************
    do k = 1, N3
        do j = 1, N2
            do i = 1, N1

                if(i == 1) then
                    gl_RAY_C(i, j, k) = gl_RAY_B(par_n_HP, j, k)
                    CYCLE
                end if

                ! Вычисляем координаты текущего луча в пространстве
                phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
                ! Вычисляем координаты точки на луче

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

                ! Создаём точку
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

    ! Цикл генерации точек на лучах O и их связывание с этими лучами ************************************************************
    do k = 1, N3
        do j = 1, N2
            !x = xx - j * (xx - par_R_LEFT)/N2
            x = xx - (DBLE(j)/N2)**par_kk3 * (xx - par_R_LEFT)
            do i = 1, N1

                ! Вычисляем координаты текущего луча в пространстве
                phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
                ! Вычисляем координаты точки на луче


                if (i <= par_n_BS - par_n_HP + 1) then
                    r = 1.3 * par_R_character + (i - 1) * par_R_character * (0.7)/(par_n_BS - par_n_HP)
                else
                    !r = 2.0 * par_R_character + (i - (par_n_BS - par_n_HP + 1)) * (par_R_END - 2.0 * par_R_character) /(N1 - (par_n_BS - par_n_HP + 1) )
                    r = 2.0 * par_R_character + (DBLE(i - (par_n_BS - par_n_HP + 1))/(N1 - (par_n_BS - par_n_HP + 1) ))**par_kk2 * (par_R_END - 2.0 * par_R_character)
                end if

                ! Создаём точку
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

    ! Цикл генерации точек на лучах K и их связывание с этими лучами ************************************************************
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

                ! Вычисляем координаты текущего луча в пространстве
                the = par_pi_8/2.0 + par_triple_point + (N2 - j + 1) * (par_pi_8/2.0 - par_triple_point)/(N2)
                phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
                ! Вычисляем координаты точки на луче

                r =  par_R0 + (par_R_character - par_R0) * (REAL(i, KIND = 4)/par_n_TS)**par_kk1


                ! Создаём точку
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

    ! Цикл генерации точек на лучах D и их связывание с этими лучами ************************************************************
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

                ! Вычисляем координаты текущего луча в пространстве
                phi = (k - 1) * 2.0_8 * par_pi_8/(N3)
                ! Вычисляем координаты точки на луче

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


                ! Создаём точку
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

    ! Цикл генерации точек на лучах E и их связывание с этими лучами ************************************************************
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


                ! Создаём точку
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

    ! Строим сами ячейки и связываем их с точками

    ! А - группа ячеек ************************************************************************************************************************
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

    ! B - группа ячеек ************************************************************************************************************************
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

    ! C - группа ячеек ************************************************************************************************************************
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

    ! Свяжем ячейки с их соседями (имеется строгий порядок)

    ! Начнём с ячеек группы А
    N3 = size(gl_Cell_A(1, 1, :))
    N2 = size(gl_Cell_A(1, :, 1))
    N1 = size(gl_Cell_A(:, 1, 1))

    do k = 1, N3
        do j = 1, N2
            do i = 1, N1

                kk = k + 1                      ! Слой ячеек выше
                if (kk > par_l_phi) kk = 1

                kk2 = k - 1                     ! Слой ячеек ниже
                if (kk2 < 1) kk2 = par_l_phi


                ! Заполняем первого соседа (по лучу вверх)
                if (i == N1) then
                    if (j < par_m_A) then
                        gl_Cell_neighbour(1, gl_Cell_A(i, j, k)) = -1   ! Граница (набегающий поток)
                    else
                        gl_Cell_neighbour(1, gl_Cell_A(i, j, k)) = -3   ! Граница верхний цилиндр
                    end if
                else
                    gl_Cell_neighbour(1, gl_Cell_A(i, j, k)) = gl_Cell_A(i + 1, j, k)
                end if

                ! Заполняем второго соседа (по лучу вниз)
                if (i == 1) then
                    continue
                else
                    gl_Cell_neighbour(2, gl_Cell_A(i, j, k)) = gl_Cell_A(i - 1, j, k)
                end if

                ! Заполняем третьего соседа (по углу в плоскости, в сторону увеличения)
                if (j == N2) then
                    if (i < par_n_TS) then
                        gl_Cell_neighbour(3, gl_Cell_A(i, j, k)) = gl_Cell_B(i, par_m_K, k)
                    else
                        gl_Cell_neighbour(3, gl_Cell_A(i, j, k)) = gl_Cell_C(i - par_n_TS + 1, 1, k)
                    end if
                else
                    gl_Cell_neighbour(3, gl_Cell_A(i, j, k)) = gl_Cell_A(i, j + 1, k)
                end if

                ! Заполняем четвёртого соседа (по углу в плоскости, в сторону уменьшения)
                if (j == 1) then
                    continue
                else
                    gl_Cell_neighbour(4, gl_Cell_A(i, j, k)) = gl_Cell_A(i, j - 1, k)
                end if

                ! Заполняем пятого и шестого соседа (вверх и вниз)
                gl_Cell_neighbour(5, gl_Cell_A(i, j, k)) = gl_Cell_A(i, j, kk)
                gl_Cell_neighbour(6, gl_Cell_A(i, j, k)) = gl_Cell_A(i, j, kk2)
            end do
        end do
    end do

    ! Ячееки группы B
    N3 = size(gl_Cell_B(1, 1, :))
    N2 = size(gl_Cell_B(1, :, 1))
    N1 = size(gl_Cell_B(:, 1, 1))

    do k = 1, N3
        do j = 1, N2
            do i = 1, N1

                kk = k + 1                      ! Слой ячеек выше
                if (kk > par_l_phi) kk = 1

                kk2 = k - 1                     ! Слой ячеек ниже
                if (kk2 < 1) kk2 = par_l_phi


                ! Заполняем первого соседа (по лучу вверх)
                if (i == N1) then
                    gl_Cell_neighbour(1, gl_Cell_B(i, j, k)) = -2   ! Выходная граница
                else
                    gl_Cell_neighbour(1, gl_Cell_B(i, j, k)) = gl_Cell_B(i + 1, j, k)
                end if

                ! Заполняем второго соседа (по лучу вниз)
                if (i == 1) then
                    continue
                else
                    gl_Cell_neighbour(2, gl_Cell_B(i, j, k)) = gl_Cell_B(i - 1, j, k)
                end if

                ! Заполняем третьего соседа (по углу в плоскости наверх на схеме)
                if (j == N2) then
                    if(i < par_n_TS) then
                        gl_Cell_neighbour(3, gl_Cell_B(i, j, k)) = gl_Cell_A(i, size(gl_Cell_A(i,:, k)), k)
                    else
                        gl_Cell_neighbour(3, gl_Cell_B(i, j, k)) = gl_Cell_C(1, i - par_n_TS + 1, k)
                    end if
                else
                    gl_Cell_neighbour(3, gl_Cell_B(i, j, k)) = gl_Cell_B(i, j + 1, k)
                end if

                ! Заполняем четвёртого соседа (по углу в плоскости)
                if (j == 1) then
                    continue
                else
                    gl_Cell_neighbour(4, gl_Cell_B(i, j, k)) = gl_Cell_B(i, j - 1, k)
                end if

                ! Заполняем пятого и шестого соседа (вверх и вниз)
                gl_Cell_neighbour(5, gl_Cell_B(i, j, k)) = gl_Cell_B(i, j, kk)
                gl_Cell_neighbour(6, gl_Cell_B(i, j, k)) = gl_Cell_B(i, j, kk2)
            end do
        end do
    end do

    ! Ячееки группы C
    N3 = size(gl_Cell_C(1, 1, :))
    N2 = size(gl_Cell_C(1, :, 1))
    N1 = size(gl_Cell_C(:, 1, 1))

    do k = 1, N3
        do j = 1, N2
            do i = 1, N1

                kk = k + 1                      ! Слой ячеек выше
                if (kk > par_l_phi) kk = 1

                kk2 = k - 1                     ! Слой ячеек ниже
                if (kk2 < 1) kk2 = par_l_phi


                ! Заполняем первого соседа (по лучу вверх)
                if (i == N1) then
                    gl_Cell_neighbour(1, gl_Cell_C(i, j, k)) = -3   ! Верхний цилиндр
                else
                    gl_Cell_neighbour(1, gl_Cell_C(i, j, k)) = gl_Cell_C(i + 1, j, k)
                end if

                ! Заполняем второго соседа (по лучу вниз)
                if (i == 1) then
                    gl_Cell_neighbour(2, gl_Cell_C(i, j, k)) = gl_Cell_B(j + par_n_TS - 1, size(gl_Cell_B(j + par_n_TS - 1, :, k)), k)
                else
                    gl_Cell_neighbour(2, gl_Cell_C(i, j, k)) = gl_Cell_C(i - 1, j, k)
                end if

                ! Заполняем третьего соседа (по углу в плоскости, в сторону увеличения угла!)
                if (j == N2) then
                    gl_Cell_neighbour(3, gl_Cell_C(i, j, k)) = -2  ! Выходная граница
                else
                    gl_Cell_neighbour(3, gl_Cell_C(i, j, k)) = gl_Cell_C(i, j + 1, k)
                end if

                ! Заполняем четвёртого соседа (по углу в плоскости, в сторону уменьшения)
                if (j == 1) then
                    gl_Cell_neighbour(4, gl_Cell_C(i, j, k)) = gl_Cell_A(i + par_n_TS - 1, size(gl_Cell_A(i + par_n_TS - 1, : , k)), k)
                else
                    gl_Cell_neighbour(4, gl_Cell_C(i, j, k)) = gl_Cell_C(i, j - 1, k)
                end if

                ! Заполняем пятого и шестого соседа (вверх и вниз)
                gl_Cell_neighbour(5, gl_Cell_C(i, j, k)) = gl_Cell_C(i, j, kk)
                gl_Cell_neighbour(6, gl_Cell_C(i, j, k)) = gl_Cell_C(i, j, kk2)
            end do
        end do
    end do

    if (par_developer_info) print*, "Build_Mesh_start: Number all points: ", size(gl_x), " == ", node - 1

    ! Создадим грани

    node = 1

    ! Начнём с ячеек группы А
    N3 = size(gl_Cell_A(1, 1, :))
    N2 = size(gl_Cell_A(1, :, 1))
    N1 = size(gl_Cell_A(:, 1, 1))

    do k = 1, N3
        do j = 1, N2
            do i = 1, N1

                kk = k + 1                      ! Слой ячеек выше
                if (kk > par_l_phi) kk = 1

                ni = gl_Cell_A(i, j, k)

                ! Заполняем 1-ую грань для этой ячейки
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

                ! Заполняем 4-ую грань для этой ячейки
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

                ! Заполняем 5-ую грань для этой ячейки
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

    ! Ячееки группы B
    N3 = size(gl_Cell_B(1, 1, :))
    N2 = size(gl_Cell_B(1, :, 1))
    N1 = size(gl_Cell_B(:, 1, 1))

    do k = 1, N3
        do j = 1, N2
            do i = 1, N1

                kk = k + 1                      ! Слой ячеек выше
                if (kk > par_l_phi) kk = 1

                ni = gl_Cell_B(i, j, k)

                ! Заполняем 1-ую грань для этой ячейки
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

                ! Заполняем 4-ую грань для этой ячейки
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

                ! Заполняем 5-ую грань для этой ячейки
                gl_all_Gran(1, node) = gl_all_Cell(5, ni)
                gl_all_Gran(2, node) = gl_all_Cell(6, ni)
                gl_all_Gran(3, node) = gl_all_Cell(7, ni)
                gl_all_Gran(4, node) = gl_all_Cell(8, ni)

                gl_Gran_neighbour(1, node) = ni
                gl_Gran_neighbour(2, node) = gl_Cell_B(i, j, kk)

                gl_Cell_gran(5, ni) = node
                gl_Cell_gran(6, gl_Cell_B(i, j, kk)) = node

                node = node + 1

                if (j == N2) then  ! 3-яя грань для границы
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

    ! Ячейки группы C
    N3 = size(gl_Cell_C(1, 1, :))
    N2 = size(gl_Cell_C(1, :, 1))
    N1 = size(gl_Cell_C(:, 1, 1))

    do k = 1, N3
        do j = 1, N2
            do i = 1, N1

                kk = k + 1                      ! Слой ячеек выше
                if (kk > par_l_phi) kk = 1

                ni = gl_Cell_C(i, j, k)

                ! Заполняем 1-ую грань для этой ячейки
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
                    gl_Gran_neighbour(2, node) = -3    ! Верхний цилиндр

                    gl_Cell_gran(1, ni) = node

                    node = node + 1
                end if

                ! Заполняем 4-ую грань для этой ячейки
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

                if (j == N2) then  ! Третья грань на границе
                    gl_all_Gran(1, node) = gl_all_Cell(3, ni)
                    gl_all_Gran(2, node) = gl_all_Cell(4, ni)
                    gl_all_Gran(3, node) = gl_all_Cell(8, ni)
                    gl_all_Gran(4, node) = gl_all_Cell(7, ni)

                    gl_Gran_neighbour(1, node) = ni
                    gl_Gran_neighbour(2, node) = -2

                    gl_Cell_gran(3, ni) = node
                    node = node + 1
                end if


                ! Заполняем 5-ую грань для этой ячейки
                gl_all_Gran(1, node) = gl_all_Cell(5, ni)
                gl_all_Gran(2, node) = gl_all_Cell(6, ni)
                gl_all_Gran(3, node) = gl_all_Cell(7, ni)
                gl_all_Gran(4, node) = gl_all_Cell(8, ni)

                gl_Gran_neighbour(1, node) = ni
                gl_Gran_neighbour(2, node) = gl_Cell_C(i, j, kk)

                gl_Cell_gran(5, ni) = node
                gl_Cell_gran(6, gl_Cell_C(i, j, kk)) = node

                node = node + 1

                if (j == 1) then  ! 4-яя грань для границы
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


    subroutine Find_inner()  ! Находит внутренние ячейки и грани (действует после calc_all_Gran)
    use STORAGE
    use GEO_PARAM
    implicit none

    integer :: n, Ngran, iter

    n = 0
    Ngran = size(gl_all_Gran(1,:))

    do  iter = 1, Ngran
        if (gl_Gran_info(iter) /= 0) n = n + 1
    end do

    allocate(gl_all_Gran_inner(n))
    n = 1

    do  iter = 1, Ngran
        if (gl_Gran_info(iter) /= 0) then
            gl_all_Gran_inner(n) = iter
            n = n + 1
        end if
    end do


    n = 0
    Ngran = size(gl_all_Cell(1,:))

    do  iter = 1, Ngran
        if (gl_Cell_info(iter) /= 2) n = n + 1
    end do

    allocate(gl_all_Cell_inner(n))
    n = 1

    do  iter = 1, Ngran
        if (gl_Cell_info(iter) /= 2) then
            gl_all_Cell_inner(n) = iter
            n = n + 1
        end if
    end do



    end subroutine Find_inner
    ! ************************************************************************************************************************************************
    ! Блок физики

    subroutine Initial_conditions()  ! Задаём начальные условия
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
            !c = DBLE(c) * par_chi_real/DBLE(r)
            !gl_Cell_par(:, gr) = (/ro, DBLE(c(1)), DBLE(c(2)), DBLE(c(3)), P_E, 0.0_8, 0.0_8, 0.0_8, ro * 1.0_8/)
            c = DBLE(c) * par_chi/DBLE(r)
            ! Задаём мультифлюид
            gl_Cell_par_MF(:, 1, gr) = (/0.000001_8, c(1), c(2), c(3), 0.000001_8/)

            !gl_Cell_par_MF(:, 2, gr) = (/0.001_8, 0.0_8, 0.0_8, 0.0_8, 0.00001_8/)
            !gl_Cell_par_MF(:, 3, gr) = (/0.001_8, 0.0_8, 0.0_8, 0.0_8, 0.00001_8/)
            !gl_Cell_par_MF(:, 4, gr) = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 0.5_8/)
        else
            !gl_Cell_par(:, gr) = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 1.0_8, 0.0_8, 0.0_8, 0.0_8, 100.0_8/)

            ! Задаём мультифлюид
            !gl_Cell_par_MF(:, 1, gr) = (/0.000001_8, 0.0_8, 0.0_8, 0.0_8, 0.000001_8/)
            !gl_Cell_par_MF(:, 2, gr) = (/0.001_8, 0.0_8, 0.0_8, 0.0_8, 0.00001_8/)
            !gl_Cell_par_MF(:, 3, gr) = (/0.001_8, 0.0_8, 0.0_8, 0.0_8, 0.00001_8/)
            !gl_Cell_par_MF(:, 4, gr) = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 0.5_8/)
        end if
    end do


    end subroutine Initial_conditions

    subroutine calc_all_Gran()   ! Программа расчёта площадей и нормалей граней и Объёмов ячеек
    ! Эта функция нужна только в неподвижной программе (т.к. в подвижном случае быстрее будет это всё считать "налету"
    ! Либо её можно использовать в проверке геометрии сетки
    use STORAGE
    use GEO_PARAM
    implicit none

    integer, automatic :: Ngran, iter
    real(8) :: p(3, 4), Vol, D
    real(8) :: a(3), b(3), c(3), S, node1(3), node2(3)
    real(8) :: dist, di, gr_center(3)
    integer :: i, j, k, ll, grc

    ! Цикл по граням - считаем площадь грани, её нормаль
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
        ! Считываем из глобальной памяти координаты точек грани

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

        ! Нужно сохранить центр грани
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

        ! Можно один раз проверить, правильно ли ориентирована нормаль!\
        
        

        call Get_center(gl_Gran_neighbour(1, iter), node1)
        node2 = (p(:,1) + p(:,2) + p(:,3) + p(:,4))/4.0
        node1 = node2 - node1
        
        
        if(DOT_PRODUCT(node1, c) < 0.0) then
            print*, "ERROR 1333 ", DOT_PRODUCT(node1, c)
            pause
        end if

        ! Нужно записать площадь грани и нормаль в общий массив!
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




    ! Теперь посчитаем объёмы ячеек
    Ngran = size(gl_all_Cell(1,:))

    do  iter = 1, Ngran   ! Пробегаемся по всем ячейкам
        Vol = 0.0
        ll = 0
        dist = 10000.0 * par_R_character
        c = 0.0
        ! Для вычисления правильного центра ячейки, нужно понять какая она (некоторые узлы пропущены)
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

        ! Вычислили центр ячейки теперь считаем объём пирамиды на каждую грань
        do j = 1, 6
            i = gl_Cell_gran(j, iter)   ! Берём по очереди все грани ячейки

            !if (iter == gl_Cell_C(1, size(gl_Cell_C(1,:,1)) ,1)) then
            !    print*, gl_Cell_gran(:, iter)
            !    pause
            !end if

            if (i == 0) CYCLE
            k = gl_all_Gran(1, i)  ! Номер первого узла грани
            b = (/gl_x(k), gl_y(k), gl_z(k)/)
            a = gl_Gran_normal(:,i)
            di = dabs(DOT_PRODUCT(c,a) - DOT_PRODUCT(a,b))  ! Расстояние от точки до плоскости
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


    ! Вычислим, что за грань и ячейки       для разделения сетки на внутреннюю и внешнюю области, которые считаются отдельно
    Ngran = size(gl_all_Gran(1,:))
    gl_Gran_info = 0

    do  iter = 1, Ngran
        i = gl_Gran_neighbour(1, iter)
        if (norm2(gl_Cell_center(:, i)) <= par_R0 * 45) then
            gl_Cell_info(i) = 0                                        ! Внутренняя сфера
            gl_Gran_info(iter) = gl_Gran_info(iter) + 1
            if (norm2(gl_Cell_center(:, i)) >= par_R0 * 40) then
                gl_Cell_info(i) = 1                                        ! Пограничный слой
            end if
        else
            gl_Cell_info(i) = 2
        end if

        i = gl_Gran_neighbour(2, iter)
        if (i < 1) CYCLE
        if (norm2(gl_Cell_center(:, i)) <= par_R0 * 45) then
            gl_Cell_info(i) = 0                                        ! Внутренняя сфера
            gl_Gran_info(iter) = gl_Gran_info(iter) + 1
            if (norm2(gl_Cell_center(:, i)) >= par_R0 * 40) then
                gl_Cell_info(i) = 1                                        ! Пограничный слой
            end if
        else
            gl_Cell_info(i) = 2
        end if

    end do


    end subroutine calc_all_Gran


    subroutine Start_GD_1(steps)
    ! Считается просто газовая динамика (без мультифлюида, без магнитных полей)
    ! Без атомов и т.д. Без ТВД
    ! Функция паспараллелена на OPEN_MP
    use STORAGE
    use GEO_PARAM
    USE OMP_LIB
    implicit none

    integer(4), intent(in) :: steps
    integer(4) :: st, gr, ngran, ncell, s1, s2, i, j, k
    real(8) :: qqq1(8), qqq2(8), qqq(8)  ! Переменные в ячейке
    real(8) :: dist, dsl, dsc, dsp
    real(8) :: POTOK(8)
    real(8) :: time, Volume, TT, U8, rad1, rad2, aa, bb, cc

    real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3

    ngran = size(gl_all_Gran(1, :))
    ncell = size(gl_all_Cell(1, :))

    call calc_all_Gran()   ! Посчитаем все объёмы, центры, площади и т.д.
    call Initial_conditions()  ! Задаём начальные условия

    time = 0.00000001

    !k = gl_Cell_C(1, size(gl_Cell_C(1,:,1)) ,1)
    !print*, gl_Cell_Volume(k)
    !
    !k = gl_Cell_C(1, size(gl_Cell_C(1,:,1)) - 1 ,1)
    !print*, gl_Cell_Volume(k)
    !
    !pause


    do st = 1, steps
        ! Делаем цикл по граням и считаем потоки через них
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

            ! Попробуем снести плотность пропорционально квадрату
            if(norm2(qqq1(2:4))/sqrt(ggg*qqq1(5)/qqq1(1)) > 5.0) then
                rad1 = norm2(gl_Cell_center(:, s1))
                rad2 = norm2(gl_Gran_center(:, gr))
                qqq1(1) = qqq1(1) * rad1**2 / rad2**2
                qqq1(5) = qqq1(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
                ! Скорости сносим в сферической С.К.
                call spherical_skorost(gl_Cell_center(1, s1), gl_Cell_center(2, s1), gl_Cell_center(3, s1), &
                    qqq1(2), qqq1(3), qqq1(4), aa, bb, cc)
                call dekard_skorost(gl_Gran_center(1, gr), gl_Gran_center(2, gr), gl_Gran_center(3, gr), &
                    aa, bb, cc, qqq1(2), qqq1(3), qqq1(4))

            end if

            if (s2 >= 1) then
                if ( norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int .and. norm2(gl_Cell_center(:, s2)) <= par_R0 * par_R_int) CYCLE
                qqq2 = gl_Cell_par(:, s2 )
                dist = min(gl_Cell_dist(s1), gl_Cell_dist(s2))

                ! Попробуем снести плотность пропорционально квадрату
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

            else  ! В случае граничных ячеек - граничные условия
                if (norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int) CYCLE
                if(s2 == -1) then  ! Набегающий поток
                    dist = gl_Cell_dist(s1)
                    qqq2 = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 1.0_8, 0.0_8, 0.0_8, 0.0_8/)
                else  ! Здесь нужны мягкие условия
                    dist = gl_Cell_dist(s1)
                    qqq2 = qqq1
                    qqq2(5) = 1.0_8
                    if(qqq2(2) > par_Velosity_inf) then
                        qqq2(2) = par_Velosity_inf ! Отсос жидкости
                    end if
                end if
            end if

            !print*, POTOK
            call chlld(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)
            !print*, POTOK
            !pause

            time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )   ! REDUCTION
            gl_Gran_POTOK(:, gr) = POTOK * gl_Gran_square(gr)
        end do
        !$omp end do

        !$omp barrier ! Синхронизация нитей
        ! Теперь цикл по ячейкам
        !$omp do private(POTOK, Volume, qqq, i, j, ro3, u3, v3, w3, p3)
        do gr = 1, ncell
            if (norm2(gl_Cell_center(:, gr)) <= par_R0 * par_R_int) CYCLE    ! Не считаем внутри сферы
            POTOK = 0.0
            Volume = gl_Cell_Volume(gr)
            qqq = gl_Cell_par(:, gr)
            ! Просуммируем потоки через грани
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

    subroutine Start_GD_2(steps)
    ! Считается просто газовая динамика + параметр Q (без мультифлюида, без магнитных полей)
    ! Без атомов и т.д. Без ТВД
    ! Функция паспараллелена на OPEN_MP
    use STORAGE
    use GEO_PARAM
    USE OMP_LIB
    implicit none

    integer(4), intent(in) :: steps
    integer(4) :: st, gr, ngran, ncell, s1, s2, i, j, k
    real(8) :: qqq1(9), qqq2(9), qqq(9)  ! Переменные в ячейке
    real(8) :: dist, dsl, dsc, dsp
    real(8) :: POTOK(9)
    real(8) :: time, Volume, TT, U8, rad1, rad2, aa, bb, cc

    real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3

    ngran = size(gl_all_Gran(1, :))
    ncell = size(gl_all_Cell(1, :))

    call calc_all_Gran()   ! Посчитаем все объёмы, центры, площади и т.д.
    call Initial_conditions()  ! Задаём начальные условия

    time = 0.00000001

    !k = gl_Cell_C(1, size(gl_Cell_C(1,:,1)) ,1)
    !print*, gl_Cell_Volume(k)
    !
    !k = gl_Cell_C(1, size(gl_Cell_C(1,:,1)) - 1 ,1)
    !print*, gl_Cell_Volume(k)
    !
    !pause


    do st = 1, steps
        ! Делаем цикл по граням и считаем потоки через них
        TT = time
        time = 100000.0
        if (mod(st, 100) == 0)  print*, st, TT

        !$omp parallel

        !$omp do private(POTOK, s1, s2, qqq1, qqq2, dist, dsl, dsp, dsc, rad1, rad2, aa, bb, cc) &
        !$omp & reduction(min:time)
        do gr = 1, ngran
            POTOK = 0.0
            s1 = gl_Gran_neighbour(1, gr)
            if(s1 < 1) print*, "ERROR 1465gdfdyhe"
            s2 = gl_Gran_neighbour(2, gr)
            qqq1 = gl_Cell_par(:, s1)

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

            end if

            if (s2 >= 1) then
                if ( norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int .and. norm2(gl_Cell_center(:, s2)) <= par_R0 * par_R_int) CYCLE
                qqq2 = gl_Cell_par(:, s2 )
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
                end if

            else  ! В случае граничных ячеек - граничные условия
                if (norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int) CYCLE
                if(s2 == -1) then  ! Набегающий поток
                    dist = gl_Cell_dist(s1)
                    qqq2 = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 1.0_8, 0.0_8, 0.0_8, 0.0_8, 100.0_8/)
                else  ! Здесь нужны мягкие условия
                    dist = gl_Cell_dist(s1)
                    qqq2 = qqq1
                    qqq2(5) = 1.0_8
                    if(qqq2(2) > par_Velosity_inf) then
                        qqq2(2) = par_Velosity_inf ! Отсос жидкости
                    end if
                end if
            end if

            !print*, POTOK
            !print*, qqq1, qqq2
            call chlld_Q(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)
            !print*, POTOK
            !pause
            !
            !call chlld(2, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
            !               0.0_8, qqq1(1:8), qqq2(1:8), dsl, dsp, dsc, POTOK(1:8))
            !print*, POTOK
            !pause

            time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )   ! REDUCTION
            gl_Gran_POTOK(:, gr) = POTOK * gl_Gran_square(gr)
        end do
        !$omp end do

        !$omp barrier ! Синхронизация нитей
        ! Теперь цикл по ячейкам
        !$omp do private(POTOK, Volume, qqq, i, j, ro3, u3, v3, w3, p3, Q3)
        do gr = 1, ncell
            if (norm2(gl_Cell_center(:, gr)) <= par_R0 * par_R_int) CYCLE    ! Не считаем внутри сферы
            POTOK = 0.0
            Volume = gl_Cell_Volume(gr)
            qqq = gl_Cell_par(:, gr)
            ! Просуммируем потоки через грани
            do i = 1, 6
                j = gl_Cell_gran(i, gr)
                if (j == 0) CYCLE
                if (gl_Gran_neighbour(1, j) == gr) then
                    POTOK = POTOK + gl_Gran_POTOK(:, j)
                else
                    POTOK = POTOK - gl_Gran_POTOK(:, j)
                end if
            end do

            ro3 = qqq(1) - time * POTOK(1) / Volume
            Q3 = qqq(9) - time * POTOK(9) / Volume
            if (ro3 <= 0.0_8) then
                print*, "Ro < 0  1490 ", ro3
                pause
            end if
            u3 = (qqq(1) * qqq(2) - time * POTOK(2) / Volume) / ro3
            v3 = (qqq(1) * qqq(3) - time * POTOK(3) / Volume) / ro3
            w3 = (qqq(1) * qqq(4) - time * POTOK(4) / Volume) / ro3
            !bx3 = (bx * Volume_do / Volume - T[now1] * (K->Potok[4] + qqq(2) * K->Potok[8]) / Volume)
            !by3 = (by * Volume_do / Volume - T[now1] * (K->Potok[5] + qqq(3) * K->Potok[8]) / Volume)
            !bz3 = (bz * Volume_do / Volume - T[now1] * (K->Potok[6] + qqq(4) * K->Potok[8]) / Volume)
            !U8 = ( qqq(5) / (ggg - 1.0) + 0.5 * qqq(1) * norm2(qqq(2:4))**2 )
            p3 = ((  ( qqq(5) / (ggg - 1.0) + 0.5 * qqq(1) * norm2(qqq(2:4))**2 ) &
                - time * POTOK(5)/ Volume) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) ) * (ggg - 1.0)

            if (p3 <= 0.0_8) then
                print*, "p < 0  1490 ", p3 , gl_Cell_center(:, gr)
                p3 = 0.000001
                !pause
            end if

            gl_Cell_par(:, gr) = (/ro3, u3, v3, w3, p3, 0.0_8, 0.0_8, 0.0_8, Q3/)

        end do
        !$omp end do

        !$omp end parallel

    end do

    end subroutine Start_GD_2

    subroutine Start_GD_3(steps)
    ! Считается газовая динамика + параметр Q + мультифлюид (без магнитных полей)
    ! Без атомов и т.д. Без ТВД
    ! Функция паспараллелена на OPEN_MP
    use STORAGE
    use GEO_PARAM
    USE OMP_LIB
    implicit none

    integer(4), intent(in) :: steps
    integer(4) :: st, gr, ngran, ncell, s1, s2, i, j, k, zone
    real(8) :: qqq1(9), qqq2(9), qqq(9)  ! Переменные в ячейке
    real(8) :: fluid1(5, 4), fluid2(5, 4)
    real(8) :: dist, dsl, dsc, dsp
    real(8) :: POTOK(9)
    real(8) :: POTOK_MF(5)
    real(8) :: POTOK_MF_all(5, 4)
    real(8) :: time, Volume, TT, U8, rad1, rad2, aa, bb, cc
    real(8) :: SOURSE(5,5)  ! Источники массы, импульса и энергии для плазмы и каждого сорта мультифлюида

    real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3

    logical :: l_1

    ngran = size(gl_all_Gran(1, :))
    ncell = size(gl_all_Cell(1, :))

    !call calc_all_Gran()   ! Посчитаем все объёмы, центры, площади и т.д.
    !call Initial_conditions()  ! Задаём начальные условия

    time = 0.000000001

    do st = 1, steps
        ! Делаем цикл по граням и считаем потоки через них
        TT = time
        time = 100000.0

        !$omp parallel

        !$omp do private(POTOK, s1, s2, qqq1, qqq2, dist, dsl, dsp, dsc, rad1, rad2, aa, bb, cc, fluid1, fluid2, POTOK_MF) &
        !$omp & reduction(min:time)
        do gr = 1, ngran
            if(gl_Gran_info(gr) == 2) CYCLE
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
                    fluid2(:, 1) = (/0.0001_8, 0.0_8, 0.0_8, 0.0_8, 0.0001_8/)
                    fluid2(:, 2) = (/0.0001_8, 0.0_8, 0.0_8, 0.0_8, 0.0001_8/)
                    fluid2(:, 3) = (/0.0001_8, 0.0_8, 0.0_8, 0.0_8, 0.0001_8/)
                    fluid2(:, 4) = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 0.5_8/)
                else  ! Здесь нужны мягкие условия
                    dist = gl_Cell_dist(s1)
                    qqq2 = qqq1
                    fluid2 = fluid1
                    qqq2(5) = 1.0_8
                    if(qqq2(2) > par_Velosity_inf) then
                        qqq2(2) = par_Velosity_inf ! Отсос жидкости
                    end if

                    !if(fluid2(2,1) > par_Velosity_inf) then
                    !    fluid2(2,1) = par_Velosity_inf ! Отсос жидкости
                    !end if
                    !if(fluid2(2,2) > par_Velosity_inf) then
                    !    fluid2(2,2) = par_Velosity_inf ! Отсос жидкости
                    !end if
                    !if(fluid2(2,3) > par_Velosity_inf) then
                    !    fluid2(2,3) = par_Velosity_inf ! Отсос жидкости
                    !end if
                    !if(fluid2(2,4) > par_Velosity_inf) then
                    !    fluid2(2,4) = par_Velosity_inf ! Отсос жидкости
                    !end if
                end if
            end if


            call chlld_Q(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)
            time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )   ! REDUCTION
            gl_Gran_POTOK(:, gr) = POTOK * gl_Gran_square(gr)

            !print*, time
            !pause

            call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                0.0_8, fluid1(:, 1), fluid2(:, 1), dsl, dsp, dsc, POTOK_MF)
            time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )
            gl_Gran_POTOK_MF(:, 1, gr) = POTOK_MF * gl_Gran_square(gr)

            !print*, time
            !print*, fluid1(:, 1)
            !print*, fluid2(:, 1)
            !print*, dsl, dsp, dsc
            !pause

            call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                0.0_8, fluid1(:, 2), fluid2(:, 2), dsl, dsp, dsc, POTOK_MF)
            time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )
            gl_Gran_POTOK_MF(:, 2, gr) = POTOK_MF * gl_Gran_square(gr)

            !print*, time
            !pause

            call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                0.0_8, fluid1(:, 3), fluid2(:, 3), dsl, dsp, dsc, POTOK_MF)
            time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )
            gl_Gran_POTOK_MF(:, 3, gr) = POTOK_MF * gl_Gran_square(gr)

            !print*, time
            !pause

            call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                0.0_8, fluid1(:, 4), fluid2(:, 4), dsl, dsp, dsc, POTOK_MF)
            time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )
            gl_Gran_POTOK_MF(:, 4, gr) = POTOK_MF * gl_Gran_square(gr)

            !if (dabs(POTOK_MF(1)) > 10.0) then
            !    print*, "POTOK ", POTOK_MF
            !    print*, "fluid1 ",fluid1(:, 4)
            !    print*, "fluid2 ",fluid2(:, 4)
            !    pause
            !end if

        end do
        !$omp end do



        !$omp barrier ! Синхронизация нитей

        !!$omp single
        !    if (mod(st, 5) == 0)  print*, st, time
        !!$omp end single

        ! Теперь цикл по ячейкам
        !$omp do private(POTOK, Volume, qqq, i, j, ro3, u3, v3, w3, p3, Q3, POTOK_MF_all, zone, l_1, fluid1, SOURSE)
        do gr = 1, ncell
            if(gl_Cell_info(gr) == 0) CYCLE
            l_1 = .TRUE.
            if (norm2(gl_Cell_center(:, gr)) <= par_R0 + (par_R_character - par_R0) * (3.0_8/par_n_TS)**par_kk1) l_1 = .FALSE.    ! Не считаем внутри сферы
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

                ! Перенормируем параметры
                qqq(2:4) = qqq(2:4) * (par_chi/par_chi_real)
                qqq(1) = qqq(1) / (par_chi/par_chi_real)**2

                if(norm2(qqq(2:4))/sqrt(ggg*qqq(5)/qqq(1)) > 2.0) then
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

            ! Перенормируем первую жидкость
            fluid1(2:4, 1) = fluid1(2:4, 1) * (par_chi/par_chi_real)
            fluid1(1, 1) = fluid1(1, 1) / (par_chi/par_chi_real)**2

            fluid1(2:4, 2) = fluid1(2:4, 2) * (3.0_8)
            fluid1(1, 2) = fluid1(1, 2) / (3.0_8)**2

            call Calc_sourse_MF(qqq, fluid1, SOURSE, zone)  ! Вычисляем источники

            ! Перенормируем первую жидкость обратно
            fluid1(2:4, 1) = fluid1(2:4, 1) / (par_chi/par_chi_real)
            fluid1(1, 1) = fluid1(1, 1) * (par_chi/par_chi_real)**2
            SOURSE(5, 2) = SOURSE(5, 2)/ (par_chi/par_chi_real)
            SOURSE(1, 2) = SOURSE(1, 2)* (par_chi/par_chi_real)

            fluid1(2:4, 2) = fluid1(2:4, 2) / (3.0_8)
            fluid1(1, 2) = fluid1(1, 2) * (3.0_8)**2
            SOURSE(5, 3) = SOURSE(5, 3)/ (3.0_8)
            SOURSE(1, 3) = SOURSE(1, 3)* (3.0_8)

            ! Перенормируем обратно
            if(zone == 1 .or. zone == 2) then
                qqq(2:4) = qqq(2:4) / (par_chi/par_chi_real)
                qqq(1) = qqq(1) * (par_chi/par_chi_real)**2
                SOURSE(5, 1) = SOURSE(5, 1)/ (par_chi/par_chi_real)
            end if


            if (l_1 == .TRUE.) then
                ro3 = qqq(1) - time * POTOK(1) / Volume
                Q3 = qqq(9) - time * POTOK(9) / Volume
                if (ro3 <= 0.0_8) then
                    print*, "Ro < 0  1490 ", ro3, gl_Cell_center(:, gr)
                    pause
                end if
                u3 = (qqq(1) * qqq(2) - time * POTOK(2) / Volume + time * SOURSE(2, 1)) / ro3
                v3 = (qqq(1) * qqq(3) - time * POTOK(3) / Volume + time * SOURSE(3, 1)) / ro3
                w3 = (qqq(1) * qqq(4) - time * POTOK(4) / Volume + time * SOURSE(4, 1)) / ro3
                p3 = ((  ( qqq(5) / (ggg - 1.0) + 0.5 * qqq(1) * norm2(qqq(2:4))**2 ) &
                    - time * POTOK(5)/ Volume + time * SOURSE(5, 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) ) * (ggg - 1.0)

                if (p3 <= 0.0_8) then
                    !print*, "p < 0  plasma 2028 ", p3 , gl_Cell_center(:, gr)
                    p3 = 0.000001
                    !pause
                end if

                gl_Cell_par(:, gr) = (/ro3, u3, v3, w3, p3, 0.0_8, 0.0_8, 0.0_8, Q3/)

            end if

            ! Теперь посчитаем законы сохранения для остальных жидкостей

            do i = 1, 4
                if (i == 1 .and. l_1 == .FALSE.) CYCLE       ! Пропускаем внутреннюю сферу для сорта 1
                if (l_1 == .FALSE.) SOURSE(:, i + 1) = 0.0       ! Пропускаем внутреннюю сферу для сорта 1
                ro3 = fluid1(1, i) - time * POTOK_MF_all(1, i) / Volume + time * SOURSE(1, i + 1)
                if (ro3 <= 0.0_8) then
                    print*, "Ro < 0  in ", i,  ro3
                    print*, "gl_Cell_center(:, gr)", gl_Cell_center(:, gr)
                    print*, "gl_Cell_info(gr)", gl_Cell_info(gr)
                    print*, fluid1(:, i)
                    print*, l_1
                    print*, SOURSE(:, i + 1)
                    print*, time
                    print*, "POTOK ", POTOK_MF_all(:, i)
                    pause
                end if
                u3 = (fluid1(1, i) * fluid1(2, i) - time * POTOK_MF_all(2, i) / Volume + time * SOURSE(2, i + 1)) / ro3
                v3 = (fluid1(1, i) * fluid1(3, i) - time * POTOK_MF_all(3, i) / Volume + time * SOURSE(3, i + 1)) / ro3
                w3 = (fluid1(1, i) * fluid1(4, i) - time * POTOK_MF_all(4, i) / Volume + time * SOURSE(4, i + 1)) / ro3
                p3 = ((  ( fluid1(5, i) / (ggg - 1.0) + 0.5 * fluid1(1, i) * norm2(fluid1(2:4, i))**2 ) &
                    - time * POTOK_MF_all(5, i)/ Volume + time * SOURSE(5, i + 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) ) * (ggg - 1.0)
                if (p3 <= 0.0_8) then
                    !print*, "p < 0  in ", i, p3, gl_Cell_center(:, gr)
                    !print*, fluid1(:, i)
                    !print*, l_1
                    !print*, "Sourse", SOURSE(:, i + 1)
                    !print*, time
                    !print*, "zone", zone, norm2(qqq(2:4))/sqrt(ggg*qqq(5)/qqq(1))
                    !print*, "POTOK ", POTOK_MF_all(:, i)
                    !print*, "qqq", qqq
                    !pause
                    p3 = 0.000001
                    !pause
                end if

                gl_Cell_par_MF(:, i, gr) = (/ro3, u3, v3, w3, p3/)
            end do

        end do
        !$omp end do

        !$omp end parallel

    end do

    end subroutine Start_GD_3

    subroutine Start_GD_3_inner(steps)
    ! Считается газовая динамика + параметр Q + мультифлюид (без магнитных полей)
    ! Без атомов и т.д. Без ТВД
    ! Функция паспараллелена на OPEN_MP
    use STORAGE
    use GEO_PARAM
    USE OMP_LIB
    implicit none

    integer(4), intent(in) :: steps
    integer(4) :: st, gr, ngran, ncell, s1, s2, i, j, k, zone, iter
    real(8) :: qqq1(9), qqq2(9), qqq(9)  ! Переменные в ячейке
    real(8) :: fluid1(5, 4), fluid2(5, 4)
    real(8) :: dist, dsl, dsc, dsp
    real(8) :: POTOK(9)
    real(8) :: POTOK_MF(5)
    real(8) :: POTOK_MF_all(5, 4)
    real(8) :: time, Volume, TT, U8, rad1, rad2, aa, bb, cc
    real(8) :: SOURSE(5,5)  ! Источники массы, импульса и энергии для плазмы и каждого сорта мультифлюида

    real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3

    logical :: l_1

    ngran = size(gl_all_Gran_inner(:))
    ncell = size(gl_all_Cell_inner(:))

    !call calc_all_Gran()   ! Посчитаем все объёмы, центры, площади и т.д.
    !call Initial_conditions()  ! Задаём начальные условия

    time = 0.000000001

    do st = 1, steps
        ! Делаем цикл по граням и считаем потоки через них
        TT = time
        time = 100000.0
        !if (mod(st, 10) == 0)  print*, "inner",st, TT

        !$omp parallel

        !$omp do private(POTOK, s1, s2, qqq1, qqq2, dist, dsl, dsp, dsc, rad1, rad2, aa, bb, cc, fluid1, fluid2, POTOK_MF, gr) &
        !$omp & reduction(min:time)
        do iter = 1, ngran
            gr = gl_all_Gran_inner(iter)
            !if(gl_Gran_info(gr) == 0) CYCLE
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
                else  ! Здесь нужны мягкие условия
                    dist = gl_Cell_dist(s1)
                    qqq2 = qqq1
                    fluid2 = fluid1
                    qqq2(5) = 1.0_8
                    if(qqq2(2) > par_Velosity_inf) then
                        qqq2(2) = par_Velosity_inf ! Отсос жидкости
                    end if

                    !if(fluid2(2,1) > par_Velosity_inf) then
                    !    fluid2(2,1) = par_Velosity_inf ! Отсос жидкости
                    !end if
                    !if(fluid2(2,2) > par_Velosity_inf) then
                    !    fluid2(2,2) = par_Velosity_inf ! Отсос жидкости
                    !end if
                    !if(fluid2(2,3) > par_Velosity_inf) then
                    !    fluid2(2,3) = par_Velosity_inf ! Отсос жидкости
                    !end if
                    !if(fluid2(2,4) > par_Velosity_inf) then
                    !    fluid2(2,4) = par_Velosity_inf ! Отсос жидкости
                    !end if
                end if
            end if


            call chlld_Q(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)
            time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )   ! REDUCTION
            gl_Gran_POTOK(:, gr) = POTOK * gl_Gran_square(gr)

            !print*, time
            !pause

            call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                0.0_8, fluid1(:, 1), fluid2(:, 1), dsl, dsp, dsc, POTOK_MF)
            time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )
            gl_Gran_POTOK_MF(:, 1, gr) = POTOK_MF * gl_Gran_square(gr)

            !print*, time
            !print*, fluid1(:, 1)
            !print*, fluid2(:, 1)
            !print*, dsl, dsp, dsc
            !pause

            call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                0.0_8, fluid1(:, 2), fluid2(:, 2), dsl, dsp, dsc, POTOK_MF)
            time = min(time, 0.8 * dist/max(dabs(dsl), dabs(dsp)) )
            gl_Gran_POTOK_MF(:, 2, gr) = POTOK_MF * gl_Gran_square(gr)

            !print*, time
            !pause

            call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                0.0_8, fluid1(:, 3), fluid2(:, 3), dsl, dsp, dsc, POTOK_MF)
            time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )
            gl_Gran_POTOK_MF(:, 3, gr) = POTOK_MF * gl_Gran_square(gr)

            !print*, time
            !pause

            call chlld_gd(1, gl_Gran_normal(1, gr), gl_Gran_normal(2, gr), gl_Gran_normal(3, gr), &
                0.0_8, fluid1(:, 4), fluid2(:, 4), dsl, dsp, dsc, POTOK_MF)
            time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )
            gl_Gran_POTOK_MF(:, 4, gr) = POTOK_MF * gl_Gran_square(gr)

            !if (dabs(POTOK_MF(1)) > 10.0) then
            !    print*, "POTOK ", POTOK_MF
            !    print*, "fluid1 ",fluid1(:, 4)
            !    print*, "fluid2 ",fluid2(:, 4)
            !    pause
            !end if

        end do
        !$omp end do

        !$omp barrier ! Синхронизация нитей

        !!$omp single
        !    if (mod(st, 50) == 0)  print*, st, time
        !!$omp end single

        ! Теперь цикл по ячейкам
        !$omp do private(POTOK, Volume, qqq, i, j, ro3, u3, v3, w3, p3, Q3, POTOK_MF_all, zone, l_1, fluid1, SOURSE, gr)
        do iter = 1, ncell
            gr = gl_all_Cell_inner(iter)
            !if(gl_Cell_info(gr) == 2) CYCLE
            l_1 = .TRUE.
            if (norm2(gl_Cell_center(:, gr)) <= par_R0 + (par_R_character - par_R0) * (2.0_8/par_n_TS)**par_kk1) l_1 = .FALSE.    ! Не считаем внутри сферы
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

                ! Перенормируем параметры
                qqq(2:4) = qqq(2:4) * (par_chi/par_chi_real)
                qqq(1) = qqq(1) / (par_chi/par_chi_real)**2

                if(norm2(qqq(2:4))/sqrt(ggg*qqq(5)/qqq(1)) > 2.0) then
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


            ! Перенормируем первую жидкость

            fluid1(2:4, 1) = fluid1(2:4, 1) * (par_chi/par_chi_real)
            fluid1(1, 1) = fluid1(1, 1) / (par_chi/par_chi_real)**2

            fluid1(2:4, 2) = fluid1(2:4, 2) * (3.0)
            fluid1(1, 2) = fluid1(1, 2) / (3.0)**2

            call Calc_sourse_MF(qqq, fluid1, SOURSE, zone)  ! Вычисляем источники

            ! Перенормируем первую жидкость обратно

            fluid1(2:4, 1) = fluid1(2:4, 1) / (par_chi/par_chi_real)
            fluid1(1, 1) = fluid1(1, 1) * (par_chi/par_chi_real)**2
            SOURSE(5, 2) = SOURSE(5, 2)/ (par_chi/par_chi_real)
            SOURSE(1, 2) = SOURSE(1, 2)* (par_chi/par_chi_real)

            fluid1(2:4, 2) = fluid1(2:4, 2) / (3.0)
            fluid1(1, 2) = fluid1(1, 2) * (3.0)**2
            SOURSE(5, 3) = SOURSE(5, 3)/ (3.0)
            SOURSE(1, 3) = SOURSE(1, 3)* (3.0)


            ! Перенормируем обратно
            if(zone == 1 .or. zone == 2) then
                qqq(2:4) = qqq(2:4) / (par_chi/par_chi_real)
                qqq(1) = qqq(1) * (par_chi/par_chi_real)**2
                SOURSE(5, 1) = SOURSE(5, 1)/ (par_chi/par_chi_real)
            end if


            if (l_1 == .TRUE.) then
                ro3 = qqq(1) - time * POTOK(1) / Volume
                Q3 = qqq(9) - time * POTOK(9) / Volume
                if (ro3 <= 0.0_8) then
                    print*, "Ro < 0  1490 ", ro3, gl_Cell_center(:, gr)
                    pause
                end if
                u3 = (qqq(1) * qqq(2) - time * POTOK(2) / Volume + time * SOURSE(2, 1)) / ro3
                v3 = (qqq(1) * qqq(3) - time * POTOK(3) / Volume + time * SOURSE(3, 1)) / ro3
                w3 = (qqq(1) * qqq(4) - time * POTOK(4) / Volume + time * SOURSE(4, 1)) / ro3
                p3 = ((  ( qqq(5) / (ggg - 1.0) + 0.5 * qqq(1) * norm2(qqq(2:4))**2 ) &
                    - time * POTOK(5)/ Volume + time * SOURSE(5, 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) ) * (ggg - 1.0)

                if (p3 <= 0.0_8) then
                    !print*, "p < 0  plasma 2028 ", p3 , gl_Cell_center(:, gr)
                    p3 = 0.000001
                    !pause
                end if

                gl_Cell_par(:, gr) = (/ro3, u3, v3, w3, p3, 0.0_8, 0.0_8, 0.0_8, Q3/)

            end if

            ! Теперь посчитаем законы сохранения для остальных жидкостей

            do i = 1, 4
                if (i == 1 .and. l_1 == .FALSE.) then
                    gl_Cell_par_MF(1, i, gr) = 0.0000001
                    gl_Cell_par_MF(5, i, gr) = 0.0000001
                    CYCLE       ! Пропускаем внутреннюю сферу для сорта 1
                end if
                if (l_1 == .FALSE.) SOURSE(:, i + 1) = 0.0       ! Пропускаем внутреннюю сферу для сорта не 1
                ro3 = fluid1(1, i) - time * POTOK_MF_all(1, i) / Volume + time * SOURSE(1, i + 1)
                if (ro3 <= 0.0_8) then
                    print*, "Ro < 0  in ", i,  ro3, gl_Cell_center(:, gr)
                    print*, fluid1(:, i)
                    print*, l_1
                    print*, SOURSE(:, i + 1)
                    print*, time
                    print*, "POTOK ", POTOK_MF_all(:, i)
                    pause
                end if
                u3 = (fluid1(1, i) * fluid1(2, i) - time * POTOK_MF_all(2, i) / Volume + time * SOURSE(2, i + 1)) / ro3
                v3 = (fluid1(1, i) * fluid1(3, i) - time * POTOK_MF_all(3, i) / Volume + time * SOURSE(3, i + 1)) / ro3
                w3 = (fluid1(1, i) * fluid1(4, i) - time * POTOK_MF_all(4, i) / Volume + time * SOURSE(4, i + 1)) / ro3
                p3 = ((  ( fluid1(5, i) / (ggg - 1.0) + 0.5 * fluid1(1, i) * norm2(fluid1(2:4, i))**2 ) &
                    - time * POTOK_MF_all(5, i)/ Volume + time * SOURSE(5, i + 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) ) * (ggg - 1.0)
                if (p3 <= 0.0_8) then
                    !print*, "p < 0  in ", i, p3, gl_Cell_center(:, gr)
                    !print*, fluid1(:, i)
                    !print*, l_1
                    !print*, "Sourse", SOURSE(:, i + 1)
                    !print*, time
                    !print*, "zone", zone, norm2(qqq(2:4))/sqrt(ggg*qqq(5)/qqq(1))
                    !print*, "POTOK ", POTOK_MF_all(:, i)
                    !print*, "qqq", qqq
                    !pause
                    p3 = 0.000001
                    !pause
                end if

                gl_Cell_par_MF(:, i, gr) = (/ro3, u3, v3, w3, p3/)
            end do

        end do
        !$omp end do

        !$omp end parallel

    end do

    end subroutine Start_GD_3_inner

    subroutine Start_GD_move()
    ! Функция газовой динамики + мультифлюид + подвижная сетка
    ! Это управляющая функция, в отличие от предыдущих версий, здесь показано как правильно вызывать процедуры для движения сетки
    use GEO_PARAM
    use STORAGE
    use My_func
    implicit none
    integer :: step, now, now2, step2
    integer(4) :: st, gr, ngran, ncell, s1, s2, i, j, k, zone, iter
    real(8) :: qqq1(9), qqq2(9), qqq(9)  ! Переменные в ячейке
    real(8) :: fluid1(5, 4), fluid2(5, 4)
    real(8) :: dist, dsl, dsc, dsp
    real(8) :: POTOK(9)
    real(8) :: POTOK_MF(5)
    real(8) :: POTOK_MF_all(5, 4)
    real(8) :: time, Volume, Volume2, TT, U8, rad1, rad2, aa, bb, cc, wc
    real(8) :: SOURSE(5,5)  ! Источники массы, импульса и энергии для плазмы и каждого сорта мультифлюида

    real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3

    logical :: l_1


    ! Сначала подготовим все массивы для правильного движения
    if (allocated(gl_Vx) == .False.) then   ! Если движение запускается впервый раз, то нужно выделить память под используемые массивы и задать значения
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
        
        ! Начальная инициализвация
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
    
    
    
    ! Нужно посчитать новые нормали граней, их площади, центры, объёмы ячеек, центры ячеек (зачем?)

    ! Запускаем глобальный цикл
    now = 2                           ! Какие параметры сейчас будут считаться (1 или 2). Они меняются по очереди
    time = 0.00000000001               ! Начальная инициализация шага по времени (в данной программе это не нужно, так как шаг вычисляется налету)
    do step = 1, 5000 * 6 * 5   ! Нужно чтобы это число было чётным!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        if (mod(step, 1000) == 0) print*, "Step = ", step , "  step_time = ", time
        now2 = now
        now = mod(now, 2) + 1
        
        TT = time
        time = 100000.0
        
        ! Вычисляем новые скорости управляющих узлов
        call Calc_move(now)   ! Записали скорости каждого узла в этот узел для последующего движения
        ! Пока только Контакт сделал
        
        ! Двигаем все узлы сетки в соответствии с расчитанными скоростями в предыдущей функции
        call Move_all(now, TT)    ! Пока только А-лучи сделал
        
        call calc_all_Gran_move(now2)   ! Расчитываются новые объёмы\площади\нормали и т.д.
        
        ngran = size(gl_all_Gran_inner(:))
        ncell = size(gl_all_Cell_inner(:))
        
        !$omp parallel

        ! Теперь по основным ячейкам
        ngran = size(gl_all_Gran(1, :))
        ncell = size(gl_all_Cell(1, :))
        
        !$omp do private(POTOK, s1, s2, qqq1, qqq2, dist, dsl, dsp, dsc, rad1, rad2, aa, bb, cc, fluid1, fluid2, POTOK_MF, wc) &
        !$omp & reduction(min:time)
        do gr = 1, ngran
            if(gl_Gran_info(gr) == 2) CYCLE
            POTOK = 0.0
            s1 = gl_Gran_neighbour(1, gr)
            s2 = gl_Gran_neighbour(2, gr)
            qqq1 = gl_Cell_par(:, s1)
            fluid1 = gl_Cell_par_MF(:, :, s1)   ! Загрузили параметры жидкостей для мультифлюида

            ! Попробуем снести плотность пропорционально квадрату
            if(norm2(qqq1(2:4))/sqrt(ggg*qqq1(5)/qqq1(1)) > 5.0) then
                rad1 = norm2(gl_Cell_center2(:, s1, now))
                rad2 = norm2(gl_Gran_center2(:, gr, now))
                qqq1(1) = qqq1(1) * rad1**2 / rad2**2
                qqq1(9) = qqq1(9) * rad1**2 / rad2**2
                qqq1(5) = qqq1(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
                ! Скорости сносим в сферической С.К.
                call spherical_skorost(gl_Cell_center2(1, s1, now), gl_Cell_center2(2, s1, now), gl_Cell_center2(3, s1, now), &  
                    qqq1(2), qqq1(3), qqq1(4), aa, bb, cc)
                call dekard_skorost(gl_Gran_center2(1, gr, now), gl_Gran_center2(2, gr, now), gl_Gran_center2(3, gr, now), &  
                    aa, bb, cc, qqq1(2), qqq1(3), qqq1(4))

                call spherical_skorost(gl_Cell_center2(1, s1, now), gl_Cell_center2(2, s1, now), gl_Cell_center2(3, s1, now), & 
                    fluid1(2, 1), fluid1(3, 1), fluid1(4, 1), aa, bb, cc)
                call dekard_skorost(gl_Gran_center2(1, gr, now), gl_Gran_center2(2, gr, now), gl_Gran_center2(3, gr, now), & 
                    aa, bb, cc, fluid1(2, 1), fluid1(3, 1), fluid1(4, 1))
            end if

            if (s2 >= 1) then
                !if ( norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int .and. norm2(gl_Cell_center(:, s2)) <= par_R0 * par_R_int) CYCLE
                qqq2 = gl_Cell_par(:, s2)
                fluid2 = gl_Cell_par_MF(:, :, s2)   ! Загрузили параметры жидкостей для мультифлюида
                dist = min(gl_Cell_dist(s1), gl_Cell_dist(s2))                                              

                ! Попробуем снести плотность пропорционально квадрату
                if(norm2(qqq2(2:4))/sqrt(ggg*qqq2(5)/qqq2(1)) > 5.0) then 
                    rad1 = norm2(gl_Cell_center2(:, s2, now))                              
                    rad2 = norm2(gl_Gran_center2(:, gr, now))
                    qqq2(1) = qqq2(1) * rad1**2 / rad2**2
                    qqq2(9) = qqq2(9) * rad1**2 / rad2**2
                    qqq2(5) = qqq2(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
                    call spherical_skorost(gl_Cell_center2(1, s2, now), gl_Cell_center2(2, s2, now), gl_Cell_center2(3, s2, now), &
                        qqq2(2), qqq2(3), qqq2(4), aa, bb, cc)
                    call dekard_skorost(gl_Gran_center2(1, gr, now), gl_Gran_center2(2, gr, now), gl_Gran_center2(3, gr, now), &
                        aa, bb, cc, qqq2(2), qqq2(3), qqq2(4))

                    call spherical_skorost(gl_Cell_center2(1, s2, now), gl_Cell_center2(2, s2, now), gl_Cell_center2(3, s2, now), &
                        fluid2(2, 1), fluid2(3, 1), fluid2(4, 1), aa, bb, cc)
                    call dekard_skorost(gl_Gran_center2(1, gr, now), gl_Gran_center2(2, gr, now), gl_Gran_center2(3, gr, now), &
                        aa, bb, cc, fluid2(2, 1), fluid2(3, 1), fluid2(4, 1))
                end if

            else  ! В случае граничных ячеек - граничные условия
                !if (norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int) CYCLE
                if(s2 == -1) then  ! Набегающий поток
                    dist = gl_Cell_dist(s1)
                    qqq2 = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 1.0_8, 0.0_8, 0.0_8, 0.0_8, 100.0_8/)
                    fluid2(:, 1) = (/0.0001_8, 0.0_8, 0.0_8, 0.0_8, 0.0001_8/)
                    fluid2(:, 2) = (/0.0001_8, 0.0_8, 0.0_8, 0.0_8, 0.0001_8/)
                    fluid2(:, 3) = (/0.0001_8, 0.0_8, 0.0_8, 0.0_8, 0.0001_8/)
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
            
            ! Нужно вычислить скорость движения грани
            wc = DOT_PRODUCT((gl_Gran_center2(:, gr, now2) -  gl_Gran_center2(:, gr, now))/TT, gl_Gran_normal2(:, gr, now))
            
            call chlld_Q(1, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
                wc, qqq1, qqq2, dsl, dsp, dsc, POTOK)
            time = min(time, 0.9 * dist/max(dabs(dsl), dabs(dsp)) )   ! REDUCTION
            gl_Gran_POTOK(:, gr) = POTOK * gl_Gran_square2(gr, now)

            call chlld_gd(1, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
                wc, fluid1(:, 1), fluid2(:, 1), dsl, dsp, dsc, POTOK_MF)
            time = min(time, 0.9 * dist/(max(dabs(dsl), dabs(dsp)) + dabs(wc)) )
            gl_Gran_POTOK_MF(:, 1, gr) = POTOK_MF * gl_Gran_square2(gr, now)

            call chlld_gd(1, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
                wc, fluid1(:, 2), fluid2(:, 2), dsl, dsp, dsc, POTOK_MF)
            time = min(time, 0.9 * dist/(max(dabs(dsl), dabs(dsp)) + dabs(wc)) )
            gl_Gran_POTOK_MF(:, 2, gr) = POTOK_MF * gl_Gran_square2(gr, now)
            
            call chlld_gd(1, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
                wc, fluid1(:, 3), fluid2(:, 3), dsl, dsp, dsc, POTOK_MF)
            time = min(time, 0.9 * dist/(max(dabs(dsl), dabs(dsp)) + dabs(wc)) )
            gl_Gran_POTOK_MF(:, 3, gr) = POTOK_MF * gl_Gran_square2(gr, now)

            call chlld_gd(1, gl_Gran_normal2(1, gr, now), gl_Gran_normal2(2, gr, now), gl_Gran_normal2(3, gr, now), &
                wc, fluid1(:, 4), fluid2(:, 4), dsl, dsp, dsc, POTOK_MF)
            time = min(time, 0.9 * dist/(max(dabs(dsl), dabs(dsp)) + dabs(wc)) )
            gl_Gran_POTOK_MF(:, 4, gr) = POTOK_MF * gl_Gran_square2(gr, now)


        end do
        !$omp end do



        !$omp barrier ! Синхронизация нитей

        !!$omp single
        !    if (mod(st, 5) == 0)  print*, st, time
        !!$omp end single

        ! Теперь цикл по ячейкам
        !$omp do private(POTOK, Volume, Volume2, qqq, i, j, ro3, u3, v3, w3, p3, Q3, POTOK_MF_all, zone, l_1, fluid1, SOURSE)
        do gr = 1, ncell
            if(gl_Cell_info(gr) == 0) CYCLE
            l_1 = .TRUE.
            if (norm2(gl_Cell_center2(:, gr, now)) <= par_R0 + (par_R_character - par_R0) * (3.0_8/par_n_TS)**par_kk1) l_1 = .FALSE.    ! Не считаем внутри сферы
            POTOK = 0.0
            SOURSE = 0.0
            POTOK_MF_all = 0.0
            Volume = gl_Cell_Volume2(gr, now)
            Volume2 = gl_Cell_Volume2(gr, now2)
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

                ! Перенормируем параметры
                qqq(2:4) = qqq(2:4) * (par_chi/par_chi_real)
                qqq(1) = qqq(1) / (par_chi/par_chi_real)**2

                if(norm2(qqq(2:4))/sqrt(ggg*qqq(5)/qqq(1)) > 2.0) then
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

            ! Перенормируем первую жидкость
            fluid1(2:4, 1) = fluid1(2:4, 1) * (par_chi/par_chi_real)
            fluid1(1, 1) = fluid1(1, 1) / (par_chi/par_chi_real)**2

            fluid1(2:4, 2) = fluid1(2:4, 2) * (3.0_8)
            fluid1(1, 2) = fluid1(1, 2) / (3.0_8)**2

            call Calc_sourse_MF(qqq, fluid1, SOURSE, zone)  ! Вычисляем источники

            ! Перенормируем первую жидкость обратно
            fluid1(2:4, 1) = fluid1(2:4, 1) / (par_chi/par_chi_real)
            fluid1(1, 1) = fluid1(1, 1) * (par_chi/par_chi_real)**2
            SOURSE(5, 2) = SOURSE(5, 2)/ (par_chi/par_chi_real)
            SOURSE(1, 2) = SOURSE(1, 2)* (par_chi/par_chi_real)

            fluid1(2:4, 2) = fluid1(2:4, 2) / (3.0_8)
            fluid1(1, 2) = fluid1(1, 2) * (3.0_8)**2
            SOURSE(5, 3) = SOURSE(5, 3)/ (3.0_8)
            SOURSE(1, 3) = SOURSE(1, 3)* (3.0_8)

            ! Перенормируем обратно
            if(zone == 1 .or. zone == 2) then
                qqq(2:4) = qqq(2:4) / (par_chi/par_chi_real)
                qqq(1) = qqq(1) * (par_chi/par_chi_real)**2
                SOURSE(5, 1) = SOURSE(5, 1)/ (par_chi/par_chi_real)
            end if


            if (l_1 == .TRUE.) then
                ro3 = qqq(1)* Volume / Volume2 - time * POTOK(1) / Volume2
                Q3 = qqq(9)* Volume / Volume2 - time * POTOK(9) / Volume2
                if (ro3 <= 0.0_8) then
                    print*, "Ro < 0  1490 ", ro3, gl_Cell_center(:, gr)
                    pause
                end if
                u3 = (qqq(1) * qqq(2)* Volume / Volume2 - time * POTOK(2) / Volume2 + time * SOURSE(2, 1)) / ro3
                v3 = (qqq(1) * qqq(3)* Volume / Volume2 - time * POTOK(3) / Volume2 + time * SOURSE(3, 1)) / ro3
                w3 = (qqq(1) * qqq(4)* Volume / Volume2 - time * POTOK(4) / Volume2 + time * SOURSE(4, 1)) / ro3
                p3 = ((  ( qqq(5) / (ggg - 1.0) + 0.5 * qqq(1) * norm2(qqq(2:4))**2 )* Volume / Volume2 &
                    - time * POTOK(5)/ Volume2 + time * SOURSE(5, 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) ) * (ggg - 1.0)

                if (p3 <= 0.0_8) then
                    !print*, "p < 0  plasma 2028 ", p3 , gl_Cell_center(:, gr)
                    p3 = 0.000001
                    !pause
                end if

                gl_Cell_par(:, gr) = (/ro3, u3, v3, w3, p3, 0.0_8, 0.0_8, 0.0_8, Q3/)

            end if

            ! Теперь посчитаем законы сохранения для остальных жидкостей

            do i = 1, 4
                if (i == 1 .and. l_1 == .FALSE.) CYCLE       ! Пропускаем внутреннюю сферу для сорта 1
                if (l_1 == .FALSE.) SOURSE(:, i + 1) = 0.0       ! Пропускаем внутреннюю сферу для сорта 1
                ro3 = fluid1(1, i)* Volume / Volume2 - time * POTOK_MF_all(1, i) / Volume2 + time * SOURSE(1, i + 1)
                if (ro3 <= 0.0_8) then
                    print*, "Ro < 0  in ", i,  ro3
                    print*, "gl_Cell_center(:, gr)", gl_Cell_center(:, gr)
                    print*, "gl_Cell_info(gr)", gl_Cell_info(gr)
                    print*, fluid1(:, i)
                    print*, l_1
                    print*, SOURSE(:, i + 1)
                    print*, time
                    print*, "POTOK ", POTOK_MF_all(:, i)
                    pause
                end if
                u3 = (fluid1(1, i) * fluid1(2, i)* Volume / Volume2 - time * POTOK_MF_all(2, i) / Volume2 + time * SOURSE(2, i + 1)) / ro3
                v3 = (fluid1(1, i) * fluid1(3, i)* Volume / Volume2 - time * POTOK_MF_all(3, i) / Volume2 + time * SOURSE(3, i + 1)) / ro3
                w3 = (fluid1(1, i) * fluid1(4, i)* Volume / Volume2 - time * POTOK_MF_all(4, i) / Volume2 + time * SOURSE(4, i + 1)) / ro3
                p3 = ((  ( fluid1(5, i) / (ggg - 1.0) + 0.5 * fluid1(1, i) * norm2(fluid1(2:4, i))**2 )* Volume / Volume2 &
                    - time * POTOK_MF_all(5, i)/ Volume2 + time * SOURSE(5, i + 1)) - 0.5 * ro3 * (u3**2 + v3**2 + w3**2) ) * (ggg - 1.0)
                if (p3 <= 0.0_8) then
                    !print*, "p < 0  in ", i, p3, gl_Cell_center(:, gr)
                    !print*, fluid1(:, i)
                    !print*, l_1
                    !print*, "Sourse", SOURSE(:, i + 1)
                    !print*, time
                    !print*, "zone", zone, norm2(qqq(2:4))/sqrt(ggg*qqq(5)/qqq(1))
                    !print*, "POTOK ", POTOK_MF_all(:, i)
                    !print*, "qqq", qqq
                    !pause
                    p3 = 0.000001
                    !pause
                end if

                gl_Cell_par_MF(:, i, gr) = (/ro3, u3, v3, w3, p3/)
            end do

        end do
        !$omp end do
        
        !$omp end parallel
        
        gl_x = gl_x2(:, now2)
        gl_y = gl_y2(:, now2)
        gl_z = gl_z2(:, now2)
        gl_Cell_Volume = gl_Cell_Volume2(:, now2)
        gl_Gran_normal = gl_Gran_normal2(:, :, now2)
        gl_Gran_center = gl_Gran_center2(:, :, now2)
        gl_Cell_center = gl_Cell_center2(:, :, now2)
        gl_Gran_square = gl_Gran_square2(:, now2)
        
        call Start_GD_3_inner(3)
        
    end do


        gl_x = gl_x2(:, 2)
        gl_y = gl_y2(:, 2)
        gl_z = gl_z2(:, 2)
        gl_Cell_Volume = gl_Cell_Volume2(:, 2)
        gl_Gran_normal = gl_Gran_normal2(:, :, 2)
        gl_Gran_center = gl_Gran_center2(:, :, 2)
        gl_Cell_center = gl_Cell_center2(:, :, 2)
        gl_Gran_square = gl_Gran_square2(:, 2)



    end subroutine Start_GD_move
                             
      
    ! ************************************************************************************************************************************************
    ! Блок геометрии, необходимой для физики - Считаем площади, объёмы

    real(8) pure function calc_square(n) ! Вычисляет площадь грани под номером n в общем списке граней
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

    pure subroutine Get_center(n, CENTER)  ! Получить массив - центр ячейки
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
    ! Блок функций для печати и проверки геометрии сетки

    subroutine Print_all_surface(str)  ! Печать поверхностей, которые выделяются (TS, HP, BS)
    
    use GEO_PARAM
    use STORAGE
    use My_func
    implicit none

    !interface
    !    real(8) pure function calc_square(n)
    !    integer, intent(in) :: n
    !    end function
    !end interface

    character, intent(in) :: str   ! Задаёт поверхность, которую надо напечатать буквой (см. ниже)
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
    elseif(str == "B") then
        m = size(gl_BS)
        open(1, file = 'BS_geometry.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z', 'M'  ZONE T= 'HP', N= ", 4 * (m), ", E =  ", 4 * (m) , ", F=FEPOINT, ET=LINESEG "
        do j = 1, m
            do i = 1, 4
                n = gl_BS(j)
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
    elseif(str == "T") then
        m = size(gl_TS)
        open(1, file = 'TS_geometry.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z', 'M'  ZONE T= 'HP', N= ", 4 * (m), ", E =  ", 4 * (m) , ", F=FEPOINT, ET=LINESEG "
        do j = 1, m
            do i = 1, 4
                n = gl_TS(j)
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

    subroutine Print_Setka_2D()  ! Печатает 2Д сетку с линиями в Техплот
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

    ! Нижняя часть сетки

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

    subroutine Print_surface_2D()
    use GEO_PARAM
    use STORAGE
    implicit none
    integer :: kk, j, yzel, N, M, L
    ! Body of Print_surface
    kk = 1
    N = size(gl_RAY_A(1, :, kk)) + size(gl_RAY_B(1, :, kk)) + size(gl_RAY_K(1, :, kk))
    M = size(gl_RAY_A(1, :, kk)) + size(gl_RAY_C(1, :, kk)) + size(gl_RAY_O(1, :, kk))
    L = M

    open(1, file = 'print_surface_2D.txt')
    write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y'  ZONE T= 'HP', N= ", 2 * N + 2 * M + 2 * L, ", E =  ", 2 * (N - 1) + 2 * (M - 1) + &
        + 2 * (L - 1), ", F=FEPOINT, ET=LINESEG "

    kk = 1
    do j = 1, size(gl_RAY_A(1, :, kk))
        yzel = gl_RAY_A(par_n_TS, j, kk)
        write(1,*) gl_x(yzel), gl_y(yzel)
    end do

    do j = 1, size(gl_RAY_B(1, :, kk))
        yzel = gl_RAY_B(par_n_TS, j, kk)
        write(1,*) gl_x(yzel), gl_y(yzel)
    end do

    do j = size(gl_RAY_K(1, :, kk)), 1, -1
        yzel = gl_RAY_K(par_n_TS, j, kk)
        write(1,*) gl_x(yzel), gl_y(yzel)
    end do

    kk = par_l_phi/2 + 1
    do j = 1, size(gl_RAY_A(1, :, kk))
        yzel = gl_RAY_A(par_n_TS, j, kk)
        write(1,*) gl_x(yzel), gl_y(yzel)
    end do

    do j = 1, size(gl_RAY_B(1, :, kk))
        yzel = gl_RAY_B(par_n_TS, j, kk)
        write(1,*) gl_x(yzel), gl_y(yzel)
    end do

    do j = size(gl_RAY_K(1, :, kk)), 1, -1
        yzel = gl_RAY_K(par_n_TS, j, kk)
        write(1,*) gl_x(yzel), gl_y(yzel)
    end do

    ! Контакт
    kk = 1
    do j = 1, size(gl_RAY_A(1, :, kk))
        yzel = gl_RAY_A(par_n_HP, j, kk)
        write(1,*) gl_x(yzel), gl_y(yzel)
    end do

    do j = 1, size(gl_RAY_C(1, :, kk))
        yzel = gl_RAY_C(1, j, kk)
        write(1,*) gl_x(yzel), gl_y(yzel)
    end do

    do j =  1, size(gl_RAY_O(1, :, kk))
        yzel = gl_RAY_O(1, j, kk)
        write(1,*) gl_x(yzel), gl_y(yzel)
    end do

    kk = par_l_phi/2 + 1
    do j = 1, size(gl_RAY_A(1, :, kk))
        yzel = gl_RAY_A(par_n_HP, j, kk)
        write(1,*) gl_x(yzel), gl_y(yzel)
    end do

    do j = 1, size(gl_RAY_C(1, :, kk))
        yzel = gl_RAY_C(1, j, kk)
        write(1,*) gl_x(yzel), gl_y(yzel)
    end do

    do j =  1, size(gl_RAY_O(1, :, kk))
        yzel = gl_RAY_O(1, j, kk)
        write(1,*) gl_x(yzel), gl_y(yzel)
    end do

    ! BS
    kk = 1
    do j = 1, size(gl_RAY_A(1, :, kk))
        yzel = gl_RAY_A(par_n_BS, j, kk)
        write(1,*) gl_x(yzel), gl_y(yzel)
    end do

    do j = 1, size(gl_RAY_C(1, :, kk))
        yzel = gl_RAY_C(par_n_BS - par_n_HP + 1, j, kk)
        write(1,*) gl_x(yzel), gl_y(yzel)
    end do

    do j =  1, size(gl_RAY_O(1, :, kk))
        yzel = gl_RAY_O(par_n_BS - par_n_HP + 1, j, kk)
        write(1,*) gl_x(yzel), gl_y(yzel)
    end do

    kk = par_l_phi/2 + 1
    do j = 1, size(gl_RAY_A(1, :, kk))
        yzel = gl_RAY_A(par_n_BS, j, kk)
        write(1,*) gl_x(yzel), gl_y(yzel)
    end do

    do j = 1, size(gl_RAY_C(1, :, kk))
        yzel = gl_RAY_C(par_n_BS - par_n_HP + 1, j, kk)
        write(1,*) gl_x(yzel), gl_y(yzel)
    end do

    do j =  1, size(gl_RAY_O(1, :, kk))
        yzel = gl_RAY_O(par_n_BS - par_n_HP + 1, j, kk)
        write(1,*) gl_x(yzel), gl_y(yzel)
    end do


    ! Connectivity list
    do j = 1, N - 1
        write(1,*) j, j + 1
    end do

    do j = N + 1, 2 * N - 1
        write(1,*) j, j + 1
    end do


    do j = 2 * N + 1, 2 * N + 1 + M - 2
        write(1,*) j, j + 1
    end do

    do j = 2 * N + 1 + M, 2 * N + 1 + 2 * M - 2
        write(1,*) j, j + 1
    end do

    do j = 2 * N + 1 + 2 * M, 2 * N + 1 + 2 * M + L - 2
        write(1,*) j, j + 1
    end do

    do j = 2 * N + 1 + 2 * M + L, 2 * N + 1 + 2 * M + L + L - 2
        write(1,*) j, j + 1
    end do

    close(1)


    end subroutine Print_surface_2D

    subroutine Print_par_2D()  ! Печатает 2Д сетку с линиями в Техплот
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
    write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z', 'rho', 'u', 'v', 'w', 'p', 'bx', 'by', 'bz', 'Volume', 'Mach', 'Q',"
    write(1,*) "'T','rho1', 'u1', 'v1', 'w1', 'p1', 'rho2',"
    write(1,*)" 'u2', 'v2', 'w2', 'p2', 'rho3', 'u3', 'v3', 'w3', 'p3', "
    write(1,*) "'rho4', 'u4', 'v4', 'w4', 'p4', ZONE T= 'HP'"


    kk = 1
    N2 = size(gl_Cell_A(1, :, 1))
    N1 = size(gl_Cell_A(:, 1, 1))
    do j = 1, N2
        do i = 1, N1
            c = gl_Cell_center(:, gl_Cell_A(i, j, kk))
            m = gl_Cell_A(i, j, kk)
            Mach = norm2(gl_Cell_par(2:4, m ))/sqrt(ggg*gl_Cell_par(5, m )/gl_Cell_par(1, m ))
            write(1,*) c, gl_Cell_par(1:8, m ), gl_Cell_Volume(m), Mach, gl_Cell_par(9, m )/gl_Cell_par(1, m ), gl_Cell_par(5, m )/gl_Cell_par(1, m ), gl_Cell_par_MF(:, :, m)
        end do
    end do
    N2 = size(gl_Cell_B(1, :, 1))
    N1 = size(gl_Cell_B(:, 1, 1))
    do j = 1, N2
        do i = 1, N1
            m = gl_Cell_B(i, j, kk)
            c = gl_Cell_center(:, m)
            Mach = norm2(gl_Cell_par(2:4, m ))/sqrt(ggg*gl_Cell_par(5, m )/gl_Cell_par(1, m ))
            write(1,*)  c, gl_Cell_par(1:8, m ), gl_Cell_Volume(m), Mach, gl_Cell_par(9, m )/gl_Cell_par(1, m ), gl_Cell_par(5, m )/gl_Cell_par(1, m ), gl_Cell_par_MF(:, :, m)
        end do
    end do
    N2 = size(gl_Cell_C(1, :, 1))
    N1 = size(gl_Cell_C(:, 1, 1))
    do j = 1, N2
        do i = 1, N1
            m = gl_Cell_C(i, j, kk)
            c = gl_Cell_center(:, m)
            Mach = norm2(gl_Cell_par(2:4, m ))/sqrt(ggg*gl_Cell_par(5, m )/gl_Cell_par(1, m ))
            write(1,*)  c, gl_Cell_par(1:8, m ), gl_Cell_Volume(m), Mach, gl_Cell_par(9, m )/gl_Cell_par(1, m ), gl_Cell_par(5, m )/gl_Cell_par(1, m ), gl_Cell_par_MF(:, :, m)
        end do
    end do

    ! Нижняя часть сетки

    kk = par_l_phi/2 + 1
    N2 = size(gl_Cell_A(1, :, 1))
    N1 = size(gl_Cell_A(:, 1, 1))
    do j = 1, N2
        do i = 1, N1
            m = gl_Cell_A(i, j, kk)
            c = gl_Cell_center(:, m)
            Mach = norm2(gl_Cell_par(2:4, m ))/sqrt(ggg*gl_Cell_par(5, m )/gl_Cell_par(1, m ))
            write(1,*)  c, gl_Cell_par(1:8, m ), gl_Cell_Volume(m), Mach, gl_Cell_par(9, m )/gl_Cell_par(1, m ), gl_Cell_par(5, m )/gl_Cell_par(1, m ), gl_Cell_par_MF(:, :, m)
        end do
    end do
    N2 = size(gl_Cell_B(1, :, 1))
    N1 = size(gl_Cell_B(:, 1, 1))
    do j = 1, N2
        do i = 1, N1
            m = gl_Cell_B(i, j, kk)
            c = gl_Cell_center(:, m)
            Mach = norm2(gl_Cell_par(2:4, m ))/sqrt(ggg*gl_Cell_par(5, m )/gl_Cell_par(1, m ))
            write(1,*)  c, gl_Cell_par(1:8, m ), gl_Cell_Volume(m), Mach, gl_Cell_par(9, m )/gl_Cell_par(1, m ), gl_Cell_par(5, m )/gl_Cell_par(1, m ), gl_Cell_par_MF(:, :, m)
        end do
    end do
    N2 = size(gl_Cell_C(1, :, 1))
    N1 = size(gl_Cell_C(:, 1, 1))
    do j = 1, N2
        do i = 1, N1
            m = gl_Cell_C(i, j, kk)
            c = gl_Cell_center(:, m)
            Mach = norm2(gl_Cell_par(2:4, m ))/sqrt(ggg*gl_Cell_par(5, m )/gl_Cell_par(1, m ))
            write(1,*)  c, gl_Cell_par(1:8, m ), gl_Cell_Volume(m), Mach, gl_Cell_par(9, m )/gl_Cell_par(1, m ), gl_Cell_par(5, m )/gl_Cell_par(1, m ), gl_Cell_par_MF(:, :, m)
        end do
    end do


    close(1)


    end subroutine Print_par_2D

    subroutine Print_Point_Plane()
    ! Точки в плоскости XOY - имеют две координаты
    ! Функция для проверки геометрии
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

    subroutine Print_All_Points()   ! Функция просто печатает все точки сетки (далее можно просмотреть их в текплоте)
    ! Variables
    use STORAGE

    open(1, file = 'all_point.txt')

    do i=1,size(gl_x)
        write(1,*) gl_x(i), gl_y(i), gl_z(i)
    end do

    close(1)

    end subroutine Print_All_Points

    subroutine Print_A_Ray(n1, n2)   ! Функция просто печатает один луч (в плоскости, в пространстве) + красит ключевые точки

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

    subroutine Print_Cell(n1, n2, n3, str)  ! Печатает ячейку с номерами n1, n2, n3 из множества ячеек str
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

    subroutine Print_all_Cell()  ! Печатает все ячейки сетки (с линями) - Либо ячейки с каким-то условием
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

    subroutine Print_cell_and_neighbour(n1, n2, n3)  ! Печатает всех соседей ячейки и саму ячейку
    use GEO_PARAM
    use STORAGE
    implicit none

    integer, intent(in) :: n1, n2, n3
    integer, automatic :: n, m, i, j, k

    n = gl_Cell_C(n1, n2, n3)  ! Получили номер ячейки в общем списке

    m = COUNT(gl_Cell_neighbour(:, n) > 0)

    if (par_developer_info) print *, "Kolichestvo sosedey = ", m

    open(1, file = 'Cell+neighbour.txt')


    write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z', 'M'  ZONE T= 'HP', N= ", 8 * (m + 1), ", E =  ", 12 * (m + 1) , ", F=FEPOINT, ET=LINESEG "

    ! Печатаем саму ячейку
    do i = 1, 8
        write(1,*) gl_x(gl_all_Cell(i, n)), gl_y(gl_all_Cell(i, n)), gl_z(gl_all_Cell(i, n)), 0
    end do

    ! Печатаем соседей
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

    subroutine Print_gran()  ! Печатает грани сетки (можно с каким-то условием)
    use GEO_PARAM
    use STORAGE
    implicit none

    integer, automatic :: n, m, i, j, k

    m = 0

    !do i = 1, size(gl_all_Gran(1,:))
    !    if ( ANY( gl_Gran_neighbour(:,i) == -2 ) ) m = m + 1   ! Это условие нужно задавать и ниже (в самом цикле печати)
    !end do
    m = 1

    if (par_developer_info) print *, "Kolichestvo graney = ", m

    open(1, file = 'Grans.txt')


    write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z', 'M'  ZONE T= 'HP', N= ", 4 * (m), ", E =  ", 4 * (m) , ", F=FEPOINT, ET=LINESEG "

    ! Печатаем
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
    ! Проверяем каждому ли лучу задали все точки (не осталось ли неинициализированных)
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
    ! Проверяем правильность соседей (чтобы была обоюдность)
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
    ! Проверяем правильность соседей граней (чтобы они реально связывали соседей, а не абы кого)
    do i = 1, size(gl_Gran_neighbour(1, :))
        if (gl_Gran_neighbour(1, i) < 1 .or. gl_Gran_neighbour(2, i) < 1) CYCLE
        j1 = gl_Gran_neighbour(1, i)
        if( ANY(gl_Cell_neighbour(:, j1) == gl_Gran_neighbour(2, i)) == .FALSE.) then
            pause "ERROR 1775"
        end if
    end do


    if (par_developer_info) print *, "Geometry_check: OK"

    ! Какие-то временные проверки можно сюда вставлять
    !print *, "B111:"
    !a1 = gl_Cell_B(1,1,1)
    !do a2 = 1, 8
    !	print*, gl_all_Cell(a2, a1)
    !end do


    end subroutine Geometry_check

    ! ****************************************************************************************************************************************************
    ! Основная программа

    program MIK
    use STORAGE
    use GEO_PARAM
    implicit none

    integer :: i

    ! Процесс построения сетки (не менять, все шаги необходимы для корректной работы)
    !call Set_STORAGE()                 ! Выделяем память под все массимы рограммы
    !call Build_Mesh_start()            ! Запускаем начальное построение сетки (все ячейки связываются, но поверхности не выделены)
    
    call Read_setka_bin(12)            ! Либо считываем сетку с файла (при этом всё равно вызывается предыдущие функции под капотом)
    
    call Find_Surface()                ! Ищем поверхности, которые будем выделять (вручную)
    call calc_all_Gran()               ! Программа расчёта объёмов ячеек, площадей и нормалей граней (обязательна здесь)
    call Find_inner()                  ! Находит ячейки внутри небольшой сферы, в которых счёт будет происходить отдельно (обязательно после 
                                       ! предыдущей функции)
    
    
    !call Print_All_Points()
    !call Print_A_Ray(1,3)
    !call Print_Point_Plane()
    call Geometry_check()              ! Проверка геометрии сетки, чтобы не было ошибок в построении
    !call Print_Cell(1, 2, 1, "C")
    !call Print_all_Cell()
    call Print_Setka_2D()
    !call Print_cell_and_neighbour(1,2,1)
    !call Print_gran()
    !call Print_all_surface("C")
    !print *, "Start_GD_2"
    !call Start_GD_2(10000)
    print *, "Start_GD_3"

    !call Start_GD_2(100000)
    !call Initial_conditions()

    print*, "Size = " , size(gl_all_Gran_inner), size(gl_all_Cell_inner)
    !call Initial_conditions()


    ! Перенормировка сортов для более быстрого счёта. Перенормируем первый и второй сорта. Делаем это ВЕЗДЕ
    gl_Cell_par_MF(2:4, 1, :) = gl_Cell_par_MF(2:4, 1, :) / (par_chi/par_chi_real)
    gl_Cell_par_MF(1, 1, :) =  gl_Cell_par_MF(1, 1, :) * (par_chi/par_chi_real)**2

    gl_Cell_par_MF(2:4, 2, :) = gl_Cell_par_MF(2:4, 2, :) / (3.0)
    gl_Cell_par_MF(1, 2, :) =  gl_Cell_par_MF(1, 2, :) * (3.0)**2

    !call Print_par_2D()
    call Start_GD_move()
    
    !do i = 1, 112000 * 3
    !    if(mod(i, 1000) == 0) print*, i
    !    call Start_GD_3(1)
    !    call Start_GD_3_inner(5)
    !end do

    ! Перенормировка сортов обратно
    gl_Cell_par_MF(2:4, 1, :) = gl_Cell_par_MF(2:4, 1, :) * (par_chi/par_chi_real)
    gl_Cell_par_MF(1, 1, :) =  gl_Cell_par_MF(1, 1, :) / (par_chi/par_chi_real)**2

    gl_Cell_par_MF(2:4, 2, :) = gl_Cell_par_MF(2:4, 2, :) * (3.0)
    gl_Cell_par_MF(1, 2, :) =  gl_Cell_par_MF(1, 2, :) / (3.0)**2


    call Print_surface_2D()
    call Print_Setka_2D()
    !call Print_all_surface("C")
    !call Print_all_surface("B")
    !call Print_all_surface("T")
    call Print_par_2D()
    call Save_setka_bin(13)
    ! Variables

    !pause
    end program MIK


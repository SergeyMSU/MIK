! �������� ������������ ������� �� �������� �����
! ������� ����� ����� ������, � ������� ��������� ������ ���������������� �����
! ����� ���������� �� ����� ������� ����� ��������� ���������������� ������ (��������)
! ����� ����� ��������������� �������� ������ ���������� ������ ��������� � ������� �����
	
! ��������
! � ������������ ������ ���� ���� "����". � ������ ����� ����� ������ ����� ����������� 4-� ����� � � ����
! ����������� � ���������� ����������� ����� �� ������ ������ �������� ����� ������� ������� 
! (�� ���� ������ �������) � �� �����, ��� ������ ��� �����
	
	module Interpolate2  ! ������ �� �������
	! ������ ������������ �� ����������
    ! �������� ������������ ����� � �������� (� ��������������� ������ �� ��������)
	! ��������� ����� ���������� �� ������������ �����
	USE GEO_PARAM
	USE STORAGE
	
	implicit none
	
	
	logical :: int2_work             ! ���������� �� ���������������� ����� (�������� �� ������ � �.�.)
	
	real(8), allocatable :: int2_coord(:, :)    ! (3, :) ����� ��������� ������������ ����� (������ ����� ����� 1)
	
	
	
    integer(4), allocatable :: int2_all_Cell(:, :)   ! ���� ����� ����� (8, :) - ������ ���������� ������� - ��� ����� ����� ������
	
	integer(4), allocatable :: int2_Cell_A(:, :, :)   ! ����� A-����� ����������� 3 (�� ���� ����, � ���� ���������, �� ���� � ������������)
	
	
	contains
	
	subroutine Int2_Initial()
	
	integer(4) :: i, j, k, N1, N2, N3
	! ��������� ���������������� ����� �� ��������
	int2_Cell_A(1, :, :) = -1  ! ����������� �����
	
	N1 = size(gl_Cell_A(:, 1, 1))
	N2 = size(gl_Cell_A(1, :, 1))
	N3 = size(gl_Cell_A(1, 1, :))
	
	do k = 1, N3
		do j = 1, N2
			do i = 1, par_n_TS - 1
				int2_Cell_A(i + 1, j, k) = gl_Cell_A(i, j, k)
				int2_coord(:, gl_Cell_A(i, j, k)) = gl_Cell_center(:, gl_Cell_A(i, j, k))
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2
			do i = par_n_TS, par_n_HP - 1
				int2_Cell_A(i + 3, j, k) = gl_Cell_A(i, j, k)
				int2_coord(:, gl_Cell_A(i, j, k)) = gl_Cell_center(:, gl_Cell_A(i, j, k))
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2
			do i = par_n_HP, par_n_BS - 1
				int2_Cell_A(i + 5, j, k) = gl_Cell_A(i, j, k)
				int2_coord(:, gl_Cell_A(i, j, k)) = gl_Cell_center(:, gl_Cell_A(i, j, k))
			end do
		end do
	end do
	
	do k = 1, N3
		do j = 1, N2
			do i = par_n_BS, par_n_END - 1
				int2_Cell_A(i + 6, j, k) = gl_Cell_A(i, j, k)
				int2_coord(:, gl_Cell_A(i, j, k)) = gl_Cell_center(:, gl_Cell_A(i, j, k))
			end do
		end do
	end do
	
	
	
	end subroutine Int2_Initial
	
	
	subroutine Int2_Print_point_plane()
	
	integer(4) :: i, j, k, N1, N2, N3
	
	open(1, file = 'Int2_print.txt')
	
	N1 = size(int2_Cell_A(:, 1, 1))
	N2 = size(int2_Cell_A(1, :, 1))
	N3 = size(int2_Cell_A(1, 1, :))
	
	do k = 1, 1
		do j = 1, N2
			do i = 1, N1
				if( int2_Cell_A(i, j, k) > 0 ) write(1,*) int2_coord(:, int2_Cell_A(i, j, k))
			end do
		end do
	end do
	
	write(1,*)
	
	close(1)
	
	end subroutine Int2_Print_point_plane
	
	
	subroutine Int2_Set_Interpolate()
	
	allocate(int2_Cell_A( size(gl_Cell_A(:, 1, 1)) + 7, size(gl_Cell_A(1, :, 1)), size(gl_Cell_A(1, 1, :)) ))
	allocate(int2_all_Cell(8, size(int2_Cell_A)))
	allocate(int2_coord(3, size(gl_x)))
	
	int2_coord = 0.0
	int2_all_Cell = -100
	int2_Cell_A = -100
	
	end subroutine Int2_Set_Interpolate
	
	
	end module Interpolate2
	
	module Interpolate                     ! ������ �������������� - �������� ���������� - �������� ���������
	! ������������ �� ������ �������� ����� (�������� ������������ �� ������, ����� ����� ��������)
    
    implicit none
	
	logical, allocatable :: Int1_set_gran(:)   ! ��������� �� ����������� ����� � ���� ������?
	integer(4), allocatable :: Int1_grans_triangle(:, :, :)  ! (3, 12, :) ����������� ����� � ������ ������
	! 12 ������ � ������ ������, ������ ����� ������� �� ��� �����, : - ���������� ����� � �����
	real(8), allocatable :: Int1_grans_normal(:, :, :)  ! (4, 12, :) ������� ������� ������ ����� � ������ + �������� D
	! n1 * x + n2 * y + n3 * z + D = 0  - ��������� ����� - ������� �������
    
	contains
	
	subroutine Set_Interpolate_main()
	! �������� ������ � ��������� ������� �� �������� ���������
	
	USE GEO_PARAM
	USE STORAGE
	implicit none
	
	integer :: i, j, k, N1, N2, N3, s, s1, s2, s3
	real(8) :: aa(3), bb(3), normal(3), nn
	
	if (allocated(gl_Cell_gran_inter) == .True.) then
        STOP "Function Set_Interpolate vizvana neskolko raz!!! Dopustimo tolko 1 raz! _9876trtyuiokjhgfde45"    
	end if
	
	if (allocated(gl_Cell_gran) == .False.) then
        STOP "Ne opredeleny osnovnye massivy v programme! _76rfvbnjki87"    
    end if
	
	! �������� ������
	allocate(gl_Cell_gran_inter, MOLD = gl_Cell_gran)
	allocate(gl_Gran_normal_inter, MOLD = gl_Gran_normal)
	allocate(gl_Gran_center_inter, MOLD = gl_Gran_center)
	allocate(gl_Gran_neighbour_inter, MOLD = gl_Gran_neighbour)
	allocate(gl_Cell_center_inter, MOLD = gl_Cell_center)
	allocate(gl_Cell_neighbour_inter, MOLD = gl_Cell_neighbour)
	allocate(gl_Cell_par_inter, MOLD = gl_Cell_par)
	allocate(gl_Cell_par_MF_inter, MOLD = gl_Cell_par_MF)
	allocate(gl_Cell_dist_inter, MOLD = gl_Cell_dist)
	
	allocate(Int1_grans_triangle(3, 12, size(gl_Cell_par(1, :))))
	allocate(Int1_grans_normal(4, 12, size(gl_Cell_par(1, :))))
	allocate(Int1_set_gran( size(gl_Cell_par(1, :)) ))
	
	gl_Cell_gran_inter = gl_Cell_gran
	gl_Gran_normal_inter = gl_Gran_normal
	gl_Gran_center_inter = gl_Gran_center
	gl_Gran_neighbour_inter = gl_Gran_neighbour
	gl_Cell_center_inter = gl_Cell_center
	gl_Cell_neighbour_inter = gl_Cell_neighbour
	gl_Cell_par_inter = gl_Cell_par
	gl_Cell_par_MF_inter = gl_Cell_par_MF
	gl_Cell_dist_inter = gl_Cell_dist
	Int1_set_gran = .False.
	
	
	! ��������� ����������� �����
	N1 = size(gl_all_Cell(1, :))
	
	do k = 1, N1
		! 1
		s = gl_Cell_neighbour(1, k)  
		if(s <= 0 .or. Int1_set_gran(max(1, s)) == .False.) then
			Int1_grans_triangle(:, 1, k) = (/ 2, 3, 6  /)
			Int1_grans_triangle(:, 2, k) = (/ 3, 7, 6  /)
		else
			if ( ANY( Int1_grans_triangle(:, 3, s) == 4) .and. ANY( Int1_grans_triangle(:, 3, s) == 5) ) then
				Int1_grans_triangle(:, 1, k) = (/ 2, 3, 6  /)
				Int1_grans_triangle(:, 2, k) = (/ 3, 7, 6  /)
			else
				Int1_grans_triangle(:, 1, k) = (/ 2, 7, 6  /)
				Int1_grans_triangle(:, 2, k) = (/ 2, 3, 7  /)
			end if
		end if
		
		! 2
		s = gl_Cell_neighbour(2, k)  
		if(s <= 0 .or. Int1_set_gran(max(1, s)) == .False.) then
			Int1_grans_triangle(:, 3, k) = (/ 5, 4, 1  /)
			Int1_grans_triangle(:, 4, k) = (/ 5, 8, 4  /)
		else
			if ( ANY( Int1_grans_triangle(:, 1, s) == 3) .and. ANY( Int1_grans_triangle(:, 1, s) == 6) ) then
				Int1_grans_triangle(:, 3, k) = (/ 5, 4, 1  /)
				Int1_grans_triangle(:, 4, k) = (/ 5, 8, 4  /)
			else
				Int1_grans_triangle(:, 3, k) = (/ 8, 4, 1  /)
				Int1_grans_triangle(:, 4, k) = (/ 5, 8, 1  /)
			end if
		end if
		
		
		! 3
		s = gl_Cell_neighbour(3, k)  
		if(s <= 0 .or. Int1_set_gran(max(1, s)) == .False.) then
			Int1_grans_triangle(:, 5, k) = (/ 4, 7, 3  /)
			Int1_grans_triangle(:, 6, k) = (/ 4, 8, 7  /)
		else
			if ( ANY( Int1_grans_triangle(:, 7, s) == 1) .and. ANY( Int1_grans_triangle(:, 7, s) == 6) ) then
				Int1_grans_triangle(:, 5, k) = (/ 4, 7, 3  /)
				Int1_grans_triangle(:, 6, k) = (/ 4, 8, 7  /)
			else
				Int1_grans_triangle(:, 5, k) = (/ 3, 8, 7  /)
				Int1_grans_triangle(:, 6, k) = (/ 4, 8, 3  /)
			end if
		end if
		
		! 4
		s = gl_Cell_neighbour(4, k)  
		if(s <= 0 .or. Int1_set_gran(max(1, s)) == .False.) then
			Int1_grans_triangle(:, 7, k) = (/ 1, 6, 5  /)
			Int1_grans_triangle(:, 8, k) = (/ 1, 2, 6  /)
		else
			if ( ANY( Int1_grans_triangle(:, 5, s) == 4) .and. ANY( Int1_grans_triangle(:, 5, s) == 7) ) then
				Int1_grans_triangle(:, 7, k) = (/ 1, 6, 5  /)
				Int1_grans_triangle(:, 8, k) = (/ 1, 2, 6  /)
			else
				Int1_grans_triangle(:, 7, k) = (/ 2, 6, 5  /)
				Int1_grans_triangle(:, 8, k) = (/ 1, 2, 5  /)
			end if
		end if
		
		
		! 5
		s = gl_Cell_neighbour(5, k)  
		if(s <= 0 .or. Int1_set_gran(max(1, s)) == .False.) then
			Int1_grans_triangle(:, 9, k) = (/ 5, 6, 8  /)
			Int1_grans_triangle(:, 10, k) = (/ 6, 7, 8  /)
		else
			if ( ANY( Int1_grans_triangle(:, 11, s) == 4) .and. ANY( Int1_grans_triangle(:, 11, s) == 2) ) then
				Int1_grans_triangle(:, 9, k) = (/ 5, 6, 8  /)
				Int1_grans_triangle(:, 10, k) = (/ 6, 7, 8  /)
			else
				Int1_grans_triangle(:, 9, k) = (/ 5, 6, 7  /)
				Int1_grans_triangle(:, 10, k) = (/ 5, 7, 8  /)
			end if
		end if
		
		
		! 6
		s = gl_Cell_neighbour(6, k)  
		if(s <= 0 .or. Int1_set_gran(max(1, s)) == .False.) then
			Int1_grans_triangle(:, 11, k) = (/ 1, 4, 2  /)
			Int1_grans_triangle(:, 12, k) = (/ 2, 4, 3  /)
		else
			if ( ANY( Int1_grans_triangle(:, 9, s) == 6) .and. ANY( Int1_grans_triangle(:, 9, s) == 8) ) then
				Int1_grans_triangle(:, 11, k) = (/ 1, 4, 2  /)
				Int1_grans_triangle(:, 12, k) = (/ 2, 4, 3  /)
			else
				Int1_grans_triangle(:, 11, k) = (/ 1, 4, 3  /)
				Int1_grans_triangle(:, 12, k) = (/ 2, 1, 3  /)
			end if
		end if
		
		Int1_set_gran(k) = .True.
		
		!gl_all_Cell(:, k)
		
		do j = 1, 12
			s1 = Int1_grans_triangle(1, j, k)
			s2 = Int1_grans_triangle(2, j, k)
			s3 = Int1_grans_triangle(3, j, k)
			
			aa = (/ gl_x(gl_all_Cell(s2, k)), gl_y(gl_all_Cell(s2, k)), gl_z(gl_all_Cell(s2, k)) /) - (/ gl_x(gl_all_Cell(s1, k)), gl_y(gl_all_Cell(s1, k)), gl_z(gl_all_Cell(s1, k)) /)
			bb = (/ gl_x(gl_all_Cell(s3, k)), gl_y(gl_all_Cell(s3, k)), gl_z(gl_all_Cell(s3, k)) /) - (/ gl_x(gl_all_Cell(s1, k)), gl_y(gl_all_Cell(s1, k)), gl_z(gl_all_Cell(s1, k)) /)
		
			normal(1) = aa(2) * bb(3) - aa(3) * bb(2) 
			normal(2) = aa(3) * bb(1) - aa(1) * bb(3) 
			normal(3) = aa(1) * bb(2) - aa(2) * bb(1) 
			
			nn = norm2(normal)
			normal = normal/nn
			
			Int1_grans_normal(1:3, j, k) = normal
			Int1_grans_normal(4, j, k) = -DOT_PRODUCT(normal, (/ gl_x(gl_all_Cell(s2, k)), gl_y(gl_all_Cell(s2, k)), gl_z(gl_all_Cell(s2, k)) /))
			
			aa = ((/ gl_x(gl_all_Cell(s1, k)), gl_y(gl_all_Cell(s1, k)), gl_z(gl_all_Cell(s1, k)) /) + &
				(/ gl_x(gl_all_Cell(s2, k)), gl_y(gl_all_Cell(s2, k)), gl_z(gl_all_Cell(s2, k)) /) + &
				(/ gl_x(gl_all_Cell(s3, k)), gl_y(gl_all_Cell(s3, k)), gl_z(gl_all_Cell(s3, k)) /))/3.0
			
			bb = gl_Cell_center(:, k)
			
			!if(DOT_PRODUCT(normal, aa - bb) < 0.0) print*, "Error normal 987656789uyghjmnbgyu"
			
		end do
		
		
		! ���� ������� �������� �������� � ������
		
		
	end do
	
	
	end subroutine Set_Interpolate_main
	
	subroutine Dell_Interpolate()
	! ������� ������� ������������, ���� ������ �� �����
	USE STORAGE
	implicit none
	
	
	DEALLOCATE(gl_Cell_gran_inter)
	DEALLOCATE(gl_Gran_normal_inter)
	DEALLOCATE(gl_Gran_center_inter)
	DEALLOCATE(gl_Gran_neighbour_inter)
	DEALLOCATE(gl_Cell_center_inter)
	DEALLOCATE(gl_Cell_neighbour_inter)
	DEALLOCATE(gl_Cell_par_inter)
	DEALLOCATE(gl_Cell_par_MF_inter)
	DEALLOCATE(Int1_grans_triangle)
	DEALLOCATE(Int1_grans_normal)
	DEALLOCATE(Int1_set_gran)
	
	
	end subroutine Dell_Interpolate
	
	subroutine Save_interpolate_bin(num)  ! ���������� ����� � �������� �����
    ! Variables
    use STORAGE
    use GEO_PARAM
    implicit none
    integer, intent(in) :: num
    character(len=5) :: name
    
    write(unit=name, fmt='(i5.5)') num
    
    open(1, file = "interpolate_save_" // name // ".bin", FORM = 'BINARY')
	
    
    write(1) 1
	
	write(1) size(gl_Cell_gran_inter(1, :))
    write(1) gl_Cell_gran_inter
	
	write(1) size(gl_Gran_normal_inter(1, :))
	write(1) gl_Gran_normal_inter
	
	write(1) size(gl_Gran_center_inter(1, :))
	write(1) gl_Gran_center_inter
	
	write(1) size(gl_Gran_neighbour_inter(1, :))
	write(1) gl_Gran_neighbour_inter
	
	write(1) size(gl_Cell_center_inter(1, :))
	write(1) gl_Cell_center_inter
	
	write(1) size(gl_Cell_neighbour_inter(1, :))
	write(1) gl_Cell_neighbour_inter
	
	write(1) size(gl_Cell_par_inter(1, :))
	write(1) gl_Cell_par_inter
	
	write(1) size(gl_Cell_par_MF_inter(1, 1, :))
	write(1) gl_Cell_par_MF_inter
	
    write(1) 1
	write(1) size(gl_Cell_dist_inter(:))
	write(1) gl_Cell_dist_inter
	
	write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    write(1) 0
    
    
    close(1)
	end subroutine Save_interpolate_bin
	
	subroutine Read_interpolate_bin(num)  
    ! Variables
    use STORAGE
    use GEO_PARAM
    implicit none
    integer, intent(in) :: num
    character(len=5) :: name
    integer :: n
    logical :: exists
    
    write(unit=name,fmt='(i5.5)') num
    
    inquire(file="interpolate_save_" // name // ".bin", exist=exists)
    
    if (exists == .False.) then
        STOP "net faila!!!"
    end if
    
    if (allocated(gl_Cell_gran_inter) == .True.) then
        STOP "Function Read_interpolate_bin vizvana neskolko raz!!! Dopustimo tolko 1 raz! 1qwedsdfr456tfghjnhu876"    
	end if
    
    
    open(1, file = "interpolate_save_" // name // ".bin", FORM = 'BINARY', ACTION = "READ")
    
	
    ! ���  �������� ����������
    read(1) n
    if (n == 1) then
		
        read(1) n
		allocate(gl_Cell_gran_inter(6, n))
		read(1)	gl_Cell_gran_inter
		
		read(1) n
		allocate(gl_Gran_normal_inter(3, n))
		read(1)	gl_Gran_normal_inter
		
		read(1) n
		allocate(gl_Gran_center_inter(3, n))
		read(1)	gl_Gran_center_inter
		
		read(1) n
		allocate(gl_Gran_neighbour_inter(2, n))
		read(1)	gl_Gran_neighbour_inter
		
		read(1) n
		allocate(gl_Cell_center_inter(3, n))
		read(1)	gl_Cell_center_inter
		
		read(1) n
		allocate(gl_Cell_neighbour_inter(6, n))
		read(1)	gl_Cell_neighbour_inter
		
		read(1) n
		allocate(gl_Cell_par_inter(9, n))
		read(1)	gl_Cell_par_inter
		
		read(1) n
		allocate(gl_Cell_par_MF_inter(5, 4, n))
		read(1)	gl_Cell_par_MF_inter
		
	end if
    
	
    read(1) n
	if (n == 1) then
		read(1) n
		allocate(gl_Cell_dist_inter(n))
		read(1) gl_Cell_dist_inter
	end if
	
	read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    read(1) n
    
    
    close(1)
	end subroutine Read_interpolate_bin
	
	subroutine Re_interpolate()
	! �������� �������� ���������� � ������� ����� �� �����������������
	USE STORAGE
	implicit none
	
	integer(4) :: i, zone, number, cell
	real(8) :: F(9), F_mf(5, 4)
	
	number = 1
	
	do i = 1, size(gl_Cell_par(1, :))
		if (gl_Cell_center(1, i) < -300.0 .or. norm2(gl_Cell_center(:, i)) > 200.0) CYCLE
		
		call Interpolate_point(gl_Cell_center(1, i), gl_Cell_center(2, i), gl_Cell_center(3, i), F, F_mf, zone, cell, number)
		number = cell
		gl_Cell_par(:, i) = F
		gl_Cell_par_MF(:, :, i) = F_mf
	end do
	
	end subroutine Re_interpolate
	
	subroutine Interpolate_point(x, y, z, F, F_mf, zone, cell, number)
	! zone - � ����� ���� ���������
	! cell - � ����� ������ ���������
	! number - � ����� ������ ������ ����� (���� �������� ��������)
	USE STORAGE
	implicit none
	real(8), intent(in) :: x, y, z
	integer(4), intent(in), optional :: number
	integer(4), intent(out) :: zone, cell
	real(8), intent(out) :: F(9), F_mf(5, 4)
	
	integer(4) :: num
	
	if (allocated(gl_Cell_gran_inter) == .False.) then
        STOP "Function Set_Interpolate NE vizvana  789ihj76uhgfrerfdcgt"    
	end if
	
	if(.not. present(number)) then
    num = 1 
    else 
    num = number
	end if
	
	zone = 1
	
	call Get_Cell_Interpolate(x, y, z, num)
	if(num == -1) then
		cell = -1
		return
	end if
	
	cell = num
	call Find_tetraedr_Interpolate(x, y, z, num, F, F_mf)
	
	end subroutine Interpolate_point
		
	subroutine Get_Cell_Interpolate(x, y, z, num) 
	! num - ��������� ������ ������ (���� ����������, �������� ������ 1)
	USE GEO_PARAM
	USE STORAGE
	implicit none
	! Variables
	
	real(8), intent(in) :: x, y, z
	integer(4), intent(in out) :: num
	
	integer(4) :: n, i, gr, s, step
	real(8) :: normal(3), center(3), r(3), dd
	
	r(1) = x
	r(2) = y
	r(3) = z
	n = num
	step = 0
	
11	CONTINUE
	step = step + 1
	
	do i = 1, 6
		s = 2
		gr = gl_Cell_gran_inter(i,n)
		if (gr <= 0) CYCLE
		normal = gl_Gran_normal_inter(:, gr)
		center = gl_Gran_center_inter(:, gr)
		
		if(gl_Gran_neighbour_inter(1, gr) /= n) then
			normal = -normal
			s = 1
		end if
		
		if (gl_Gran_neighbour_inter(s, gr) <= 0) CYCLE
		
		dd = DOT_PRODUCT(normal, r - center)
		if( dd > 0.01 * gl_Cell_dist_inter(n)) then
			n = gl_Gran_neighbour_inter(s, gr)
			if (step > 10000) then
				!print*, "gr, n, i = ", gr, n, i
				!print*, "neybor ", gl_Gran_neighbour_inter(1, gr), gl_Gran_neighbour_inter(2, gr)
				!print*, r, dd
				!print*, center
				!print*, normal
				!pause
				print*, "10000 shagov 78987yhjuyhjiu"
				num = n
				EXIT
				!return
			end if
			
			GO TO 11 
		end if
	end do
	
	
	num = n
	
	end subroutine Get_Cell_Interpolate
	
	subroutine Find_tetraedr_Interpolate(x, y, z, num, F, F_mf)
	! num - ����� ������, � ������� ��������� �����
	USE GEO_PARAM
	USE STORAGE
	USE ieee_arithmetic
	implicit none
	! Variables
	real(8), intent(in) :: x, y, z
	integer(4), intent(in) :: num
	real(8), intent(out) :: F(9), F_mf(5, 4)
	
	integer(4) :: i, j, s1, s2, s3, i_best
	real(8) :: r(3), A(3), B(3), C(3), D(3), M(4, 4), det, aa(3), bb(3), normal(3), dist, dist2
	real(8), dimension(4, 4) :: Minv
	real(8), dimension(4) :: work  ! work array for LAPACK
	real(8), dimension(1, 4) :: vec
	integer, dimension(4) :: ipiv   ! pivot indices
	integer :: n, info
	
	! print*, x, y, z, "START"
	
	!$inter external DGETRF
	!$inter external DGETRI
	
	! print*, "START"
	
	r(1) = x
	r(2) = y
	r(3) = z
	
	A = gl_Cell_center_inter(:, num)
	
	
	do i = 1, 8
		s1 = gl_Cell_neighbour_inter( par_tet(1, i), num)
		s2 = gl_Cell_neighbour_inter( par_tet(2, i), num)
		s3 = gl_Cell_neighbour_inter( par_tet(3, i), num)
		
		if(s1 <= 0 .or. s2 <= 0 .or. s3 <= 0) CYCLE
		
		B = gl_Cell_center_inter(:, s1)
		C = gl_Cell_center_inter(:, s2)
		D = gl_Cell_center_inter(:, s3)
		
		! 1
		aa = B - C
		bb = D - C
		
		normal(1) = aa(2) * bb(3) - aa(3) * bb(2) 
		normal(2) = aa(3) * bb(1) - aa(1) * bb(3) 
		normal(3) = aa(1) * bb(2) - aa(2) * bb(1) 
		
		if ( DOT_PRODUCT(normal, r - C) * DOT_PRODUCT(normal, A - C) < 0) CYCLE
		
		! 2
		aa = C - A
		bb = D - A
		
		normal(1) = aa(2) * bb(3) - aa(3) * bb(2) 
		normal(2) = aa(3) * bb(1) - aa(1) * bb(3) 
		normal(3) = aa(1) * bb(2) - aa(2) * bb(1) 
		
		if ( DOT_PRODUCT(normal, r - A) * DOT_PRODUCT(normal, B - A) < 0) CYCLE
		
		! 3
		aa = B - A
		bb = D - A
		
		normal(1) = aa(2) * bb(3) - aa(3) * bb(2) 
		normal(2) = aa(3) * bb(1) - aa(1) * bb(3) 
		normal(3) = aa(1) * bb(2) - aa(2) * bb(1)  
		
		if ( DOT_PRODUCT(normal, r - A) * DOT_PRODUCT(normal, C - A) < 0) CYCLE
		
		! 4
		aa = B - C
		bb = A - C
		
		normal(1) = aa(2) * bb(3) - aa(3) * bb(2) 
		normal(2) = aa(3) * bb(1) - aa(1) * bb(3) 
		normal(3) = aa(1) * bb(2) - aa(2) * bb(1) 
		
		if ( DOT_PRODUCT(normal, r - C) * DOT_PRODUCT(normal, A - C) < 0) CYCLE
		
		EXIT
		
	end do
	
	if (i == 9) then  ! ���� �� ����� ������ ��������, ������ ����� ����� ��� ���������� ��-�� �������� ������
		! ����� ������ ���������� ��������
		i_best = 1
		dist2 = 10000000000.0_8
		
		do j = 1, 8
			
		dist = 0.0_8
		s1 = gl_Cell_neighbour_inter( par_tet(1, j), num)
		s2 = gl_Cell_neighbour_inter( par_tet(2, j), num)
		s3 = gl_Cell_neighbour_inter( par_tet(3, j), num)
		
		if(s1 <= 0 .or. s2 <= 0 .or. s3 <= 0) CYCLE
		
		B = gl_Cell_center_inter(:, s1)
		C = gl_Cell_center_inter(:, s2)
		D = gl_Cell_center_inter(:, s3)
		
		! 1
		aa = B - C
		bb = D - C
		
		normal(1) = aa(2) * bb(3) - aa(3) * bb(2) 
		normal(2) = aa(3) * bb(1) - aa(1) * bb(3) 
		normal(3) = aa(1) * bb(2) - aa(2) * bb(1) 
		
		if ( DOT_PRODUCT(normal, r - C) * DOT_PRODUCT(normal, A - C) < 0) then
			if ( dabs(DOT_PRODUCT(normal, r - C)) > dist ) then
				dist = dabs(DOT_PRODUCT(normal, r - C))
			end if
		end if
		
		! 2
		aa = C - A
		bb = D - A
		
		normal(1) = aa(2) * bb(3) - aa(3) * bb(2) 
		normal(2) = aa(3) * bb(1) - aa(1) * bb(3) 
		normal(3) = aa(1) * bb(2) - aa(2) * bb(1) 
		
		if ( DOT_PRODUCT(normal, r - A) * DOT_PRODUCT(normal, B - A) < 0) then
			if ( dabs(DOT_PRODUCT(normal, r - A)) > dist ) then
				dist = dabs(DOT_PRODUCT(normal, r - A))
			end if
		end if
		
		! 3
		aa = B - A
		bb = D - A
		
		normal(1) = aa(2) * bb(3) - aa(3) * bb(2) 
		normal(2) = aa(3) * bb(1) - aa(1) * bb(3) 
		normal(3) = aa(1) * bb(2) - aa(2) * bb(1)  
		
		if ( DOT_PRODUCT(normal, r - A) * DOT_PRODUCT(normal, C - A) < 0) then
			if ( dabs(DOT_PRODUCT(normal, r - A)) > dist ) then
				dist = dabs(DOT_PRODUCT(normal, r - A))
			end if
		end if
		
		! 4
		aa = B - C
		bb = A - C
		
		normal(1) = aa(2) * bb(3) - aa(3) * bb(2) 
		normal(2) = aa(3) * bb(1) - aa(1) * bb(3) 
		normal(3) = aa(1) * bb(2) - aa(2) * bb(1) 
		
		if ( DOT_PRODUCT(normal, r - C) * DOT_PRODUCT(normal, A - C) < 0) then
			if ( dabs(DOT_PRODUCT(normal, r - C)) > dist ) then
				dist = dabs(DOT_PRODUCT(normal, r - C))
			end if
		end if
		
		if(dist < dist2) then
			i_best = j
			dist2 = dist
		end if
	end do
		
		i = i_best
		
	end if
	
	!print*, "i = ", i
	!pause
	
	s1 = gl_Cell_neighbour_inter( par_tet(1, i), num)
	s2 = gl_Cell_neighbour_inter( par_tet(2, i), num)
	s3 = gl_Cell_neighbour_inter( par_tet(3, i), num)
	
	B = gl_Cell_center_inter(:, s1)
	C = gl_Cell_center_inter(:, s2)
	D = gl_Cell_center_inter(:, s3)
	
	!print*, A
	!print*, B
	!print*, C
	!print*, D
	!PAUSE
	
	M(:, 1) = 1.0_8
	M(1, 2:4) = A
	M(2, 2:4) = B
	M(3, 2:4) = C
	M(4, 2:4) = D
	
	Minv = M
	n = size(M,1)
	
	!$inter call DGETRF(n, n, Minv, n, ipiv, info)
	
	if (info /= 0) then
     stop 'Matrix is numerically singular! wfergt56y454erthy5t'
	end if
	
	!$inter call DGETRI(n, Minv, n, ipiv, work, n, info)
	
	if (info /= 0) then
     stop 'Matrix inversion failed! frfg543565tttg'
	end if
	
	
	
	vec(1, 1) = 1.0_8
	vec(1, 2:4) = r
	
	!print*, Minv
	
	vec = MATMUL(vec, Minv)
	
	!print*, "___"
	!print*, Vec
	!print*, "___"
	!print*, gl_Cell_par(1, num), gl_Cell_par(1, s1), gl_Cell_par(1, s2), gl_Cell_par(1, s3)
	
	F = vec(1, 1) * gl_Cell_par_inter(:, num) + vec(1, 2) * gl_Cell_par_inter(:, s1) + vec(1, 3) * gl_Cell_par_inter(:, s2) + vec(1, 4) * gl_Cell_par_inter(:, s3)
	F_mf = vec(1, 1) * gl_Cell_par_MF_inter(:, :, num) + vec(1, 2) * gl_Cell_par_MF_inter(:, :, s1) + vec(1, 3) * gl_Cell_par_MF_inter(:, :, s2) + vec(1, 4) * gl_Cell_par_MF_inter(:, :, s3)
	
	if(ieee_is_nan(F(1))) STOP "NUN 678ojhyuikjhikj"
	if(ieee_is_nan(F(2))) STOP "NUN vbnjko9876t"
	if(ieee_is_nan(F(5))) STOP "NUN 13esasdfg"
	if(ieee_is_nan(F_mf(1, 1))) STOP "NUN 9ol,mnbgtyu"
	if(ieee_is_nan(F_mf(1, 2))) STOP "NUN 876rfghjg"
	if(ieee_is_nan(F_mf(2, 1))) STOP "NUN bnhju78ijh"
	if(ieee_is_nan(F_mf(2, 2))) STOP "NUN zxdse45"
	
	!print*, "___"
	!print*, F(1)
	
	end subroutine Find_tetraedr_Interpolate
	
	subroutine Streem_line(x, y, z, num)
	! Variables
	
	! Body of Streem_line
	
	integer, intent(in) :: num
	real(8), intent(in) :: x, y, z
	
	real(8) :: r(3), V(3)
    character(len=5) :: name
	real(8) :: F(9), F_mf(5, 4)
	integer(4) :: zone, cell, number, N, i
	
	r(1) = x
	r(2) = y
	r(3) = z
	number = 1
	cell = 1
	V = 0.0
	N = 0
    
    write(unit=name,fmt='(i5.5)') num
    
    open(1, file = "streem_line_" // name // ".txt")
	
	
11	continue	
	N = N + 1
	r = r + V
	number = cell
	call Interpolate_point(r(1), r(2), r(3), F, F_mf, zone, cell, number)
	V = F(2:4)
	
	if (r(1) > -400.0 .and. sqrt(r(2)**2 + r(3)**2) < 250.0 .and. r(1) < 150.0 ) then
		GO TO 11
	end if
	
	write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Z'  ZONE T= 'HP', N= ", N, ", E =  ", N - 1 , ", F=FEPOINT, ET=LINESEG "
	
	r(1) = x
	r(2) = y
	r(3) = z
	number = 1
	cell = 1
	V = 0.0
	
	12	continue	
	r = r + V
	number = cell
	call Interpolate_point(r(1), r(2), r(3), F, F_mf, zone, cell, number)
	V = F(2:4)
	write(1,*) r(1), r(2), r(3)
	
	if (r(1) > -400.0 .and. sqrt(r(2)**2 + r(3)**2) < 250.0 .and. r(1) < 150.0 ) then
		GO TO 12
	end if
	
	
	do i = 1, N - 1
		write(1,*) i, i + 1
	end do
	
	close(1)
	
	end subroutine Streem_line
	
	end module Interpolate
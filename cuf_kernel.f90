! ��� ������� �� ������ ����� �������� ���������� ����� �� cuf

	 include "cuf_Solvers.cuf"
	
module MY_CUDA
	 use cudafor
	 
	 ! ���������
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
	 
	! ������ ����� ����������� ��������, ����������� �������� �� �����
	 integer(4), device, allocatable :: dev_gl_Gran_neighbour(:,:)
	 integer(4), device, allocatable :: dev_gl_Gran_info(:)         ! (:) ������������� ����� (��. �����)
	 integer(4), device, allocatable :: dev_gl_Cell_info(:)
	 real(8), device, allocatable :: dev_gl_Cell_par(:, :)
	 real(8), device, allocatable :: dev_gl_Cell_par_MF(:,:,:)
	 real(8), device, allocatable :: dev_gl_Cell_center(:, :)
	 real(8), device, allocatable :: dev_gl_Gran_normal(:,:)       ! (3, :) ������� �����   4444444444444444
     real(8), device, allocatable :: dev_gl_Gran_square(:)         ! (:) ������� �����
	 real(8), device, allocatable :: dev_gl_Gran_center(:,:)
	 real(8), device, allocatable :: dev_gl_Cell_dist(:)
	 real(8), device, allocatable :: dev_gl_Cell_Volume(:)
	 real(8), device, allocatable :: dev_gl_Gran_POTOK(:,:)       ! (8, :) ����� �����
	 real(8), device, allocatable :: dev_gl_Gran_POTOK_MF(:,:,:)       ! (5, 4, :) ����� ����� �������������� ���������
	 integer(4), device, allocatable :: dev_gl_Cell_gran(:,:)
	 integer(4), device, allocatable :: dev_gl_all_Gran_inner(:)
	 integer(4), device, allocatable :: dev_gl_all_Cell_inner(:)
	 character,  device, allocatable :: dev_gl_Cell_type(:)
	 integer(4),  device, allocatable :: dev_gl_Cell_number(:, :)
	 
	 real(8), device :: time_all
	 real(8), device :: time_step
	 
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
	
	attributes(device) subroutine chlld_Q(n_state, al, be, ge, &
                                 w, qqq1, qqq2, &
                                 dsl, dsp, dsc, &
                                 qqq)
      implicit real*8 (a-h,o-z)
      real(8), intent(out) :: dsl, dsp, dsc
      real(8), intent(in) :: al, be, ge, w
      integer(4), intent(in) :: n_state
      dimension qqq(9),qqq1(9),qqq2(9)
	 end subroutine chlld_Q
	 
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
		
		 allocate(dev_gl_Gran_neighbour(2, Ngran))
		 allocate(dev_gl_Gran_info(Ngran))
		 allocate(dev_gl_Cell_par(9, Ncell))
		 allocate(dev_gl_Cell_par_MF(5, 4, Ncell))
		 allocate(dev_gl_Cell_center(3, Ncell))
		 allocate(dev_gl_Gran_normal(3, Ngran))
		 allocate(dev_gl_Gran_square(Ngran))
		 allocate(dev_gl_Gran_center(3, Ngran))
		 allocate(dev_gl_Cell_dist(Ncell))
		 allocate(dev_gl_Cell_Volume(Ncell))
		 allocate(dev_gl_Gran_POTOK(9, Ngran))
		 allocate(dev_gl_Gran_POTOK_MF(5, 4, Ngran))
		 allocate(dev_gl_Cell_info(Ncell))
		 allocate(dev_gl_Cell_gran(6, Ncell))
		 allocate(dev_gl_all_Gran_inner( size(gl_all_Gran_inner(:) )))
		 allocate(dev_gl_all_Cell_inner( size(gl_all_Cell_inner(:) )))
		 allocate(dev_gl_Cell_type( size(gl_Cell_type(:))  ))
		 allocate(dev_gl_Cell_number(3,  size(gl_Cell_number(1, :))  ))
		 
		 time_all = 0.0
		 time_step = 10000.0
	
	end subroutine Set_CUDA
	
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
		 dev_gl_Cell_info = gl_Cell_info
		 dev_gl_Cell_gran = gl_Cell_gran
		 dev_gl_all_Gran_inner = gl_all_Gran_inner
		 dev_gl_all_Cell_inner = gl_all_Cell_inner
		 dev_gl_Cell_type = gl_Cell_type
		 dev_gl_Cell_number = gl_Cell_number
		 
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
		 dev_par_pi_8 = par_pi_8
		 
		 
	
	end subroutine Send_data_to_Cuda
	
	subroutine Send_data_to_Host()
		use GEO_PARAM
		use STORAGE
		gl_Cell_par = dev_gl_Cell_par
		gl_Cell_par_MF = dev_gl_Cell_par_MF
	end subroutine Send_data_to_Host
	
	
	attributes(device) real(8) function dev_norm2(x)
    implicit none
	real(8), device :: x(*)

    dev_norm2 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)
	return
	end function dev_norm2
	
	attributes(device) subroutine dev_Calc_sourse_MF(plasma, fluid, sourse, zone)  ! ��������� �������������� ���������
    implicit none
    real(8), intent(in) :: plasma(9)
    real(8), intent(in) :: fluid(5,4)
    real(8), intent(out) :: sourse(5,5)
    integer(4), intent(in) :: zone
    
    integer(4) :: i, kk(4)
    real(8) :: U_M_H(4), U_H(4), sigma(4), nu(4), S1, S2
	
	associate(par_a_2=>dev_par_a_2)
	associate(par_n_p_LISM=>dev_par_n_p_LISM)
	associate(par_Kn=>dev_par_Kn)
	associate(par_n_H_LISM_=>dev_par_n_H_LISM_)
    
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
    
	
	end associate
	end associate
	end associate
	end associate
    
	end subroutine dev_Calc_sourse_MF
	
	attributes(global) subroutine dev_retime()
	    time_step = 10000.0
	end subroutine dev_retime
	
	end module MY_CUDA
	
	attributes(global) subroutine CUF_hellow()
	
	print*, "Hellow from CUDA"
	
	end subroutine CUF_hellow
	
	
	subroutine CUDA_START_GD_move()
	! �� ��������� �����! ��� ��������� �����, � �������������
	use STORAGE
    use GEO_PARAM
	use MY_CUDA
	
	
	end subroutine CUDA_START_GD_move
	
	
	subroutine CUDA_START_GD_3()
	! �� ����������� �����! ��� ��������� �����, � �������������
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
	
	step_all = 20000 * 5 * 6 * 2
	
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
	! ������� ������� ����� ����� � ������ ���� ������� � ������� ��� ������������� �����
	! ����������� �����, ������� �������� + �����������
	! ���������, ���������� � �����, ���������� cuf-����
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
    real(8) :: qqq1(9), qqq2(9), qqq(9)  ! ���������� � ������
    real(8) :: fluid1(5, 4), fluid2(5, 4)
    real(8) :: dist, dsl, dsc, dsp
    real(8) :: POTOK(9), ttest(3)
    real(8) :: POTOK_MF(5)
    real(8) :: POTOK_MF_all(5, 4)
    real(8) :: time, Volume, TT, U8, rad1, rad2, aa, bb, cc
    real(8) :: SOURSE(5,5)  ! ��������� �����, �������� � ������� ��� ������ � ������� ����� ������������
    real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3
	
	time = 100000.0
	iter = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! ����� ������
	
	if (iter > dev_Ngran_inner) return
	
	!if(gr == 1) then
	!	print*, "Hellow from potok (1, 1) ", dev_Ngran
	!end if
	
	
	! ----------------------------------------------------------- ����� ������ ����������� ��� �� ����� ----------------------
	gr = gl_all_Gran_inner(iter)
	
            POTOK = 0.0
            s1 = gl_Gran_neighbour(1, gr)
            s2 = gl_Gran_neighbour(2, gr)
            qqq1 = gl_Cell_par(:, s1)
            fluid1 = gl_Cell_par_MF(:, :, s1)   ! ��������� ��������� ��������� ��� ������������

            ! ��������� ������ ��������� ��������������� ��������
            if(norm2(qqq1(2:4))/sqrt(ggg*qqq1(5)/qqq1(1)) > 5.0) then
                rad1 = norm2(gl_Cell_center(:, s1))
                rad2 = norm2(gl_Gran_center(:, gr))
                qqq1(1) = qqq1(1) * rad1**2 / rad2**2
                qqq1(9) = qqq1(9) * rad1**2 / rad2**2
                qqq1(5) = qqq1(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
                ! �������� ������ � ����������� �.�.
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
                fluid2 = gl_Cell_par_MF(:, :, s2)   ! ��������� ��������� ��������� ��� ������������
                dist = min(gl_Cell_dist(s1), gl_Cell_dist(s2))

                ! ��������� ������ ��������� ��������������� ��������
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

            else  ! � ������ ��������� ����� - ��������� �������
                !if (norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int) CYCLE
                if(s2 == -1) then  ! ���������� �����
                    dist = gl_Cell_dist(s1)
                    qqq2 = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 1.0_8, 0.0_8, 0.0_8, 0.0_8, 100.0_8/)
                    fluid2(:, 1) = (/0.000001_8, 0.0_8, 0.0_8, 0.0_8, 0.000001_8/)
                    fluid2(:, 2) = (/0.000001_8, 0.0_8, 0.0_8, 0.0_8, 0.000001_8/)
                    fluid2(:, 3) = (/0.000001_8, 0.0_8, 0.0_8, 0.0_8, 0.000001_8/)
                    fluid2(:, 4) = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 0.5_8/)
                else  ! ����� ����� ������ ������� (��� ������ ������)
                    dist = gl_Cell_dist(s1)
                    qqq2 = qqq1
                    fluid2 = fluid1
                    qqq2(5) = 1.0_8
                    if(qqq2(2) > par_Velosity_inf) then
                        qqq2(2) = par_Velosity_inf ! ����� ��������
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
			time =  atomicmin(time_step, time)   ! ��������� �������� ������ ������������ ��������
			end if

	
	
	end subroutine CUF_GD_3_grans_inner
	
	attributes(global) subroutine CUF_GD_3_grans()
	! ������� ������� ����� ����� � ������ ���� ������� � ������� ��� ������������� �����
	! ����������� �����, ������� �������� + �����������
	! ���������, ���������� � �����, ���������� cuf-����
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
    real(8) :: qqq1(9), qqq2(9), qqq(9)  ! ���������� � ������
    real(8) :: fluid1(5, 4), fluid2(5, 4)
    real(8) :: dist, dsl, dsc, dsp
    real(8) :: POTOK(9), ttest(3)
    real(8) :: POTOK_MF(5)
    real(8) :: POTOK_MF_all(5, 4)
    real(8) :: time, Volume, TT, U8, rad1, rad2, aa, bb, cc
    real(8) :: SOURSE(5,5)  ! ��������� �����, �������� � ������� ��� ������ � ������� ����� ������������
    real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3
	
	time = 100000.0
	gr = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! ����� ������
	
	!if(gr == 1) then
	!	print*, "Hellow from potok (1, 1) ", dev_Ngran
	!end if
	
	
	! ----------------------------------------------------------- ����� ������ ����������� ��� �� ����� ----------------------
	if (gr > dev_Ngran) return
	
	if(gl_Gran_info(gr) == 2) return
	
	s1 = gl_Gran_neighbour(1, gr)
	s2 = gl_Gran_neighbour(2, gr)
	qqq1 = gl_Cell_par(:, s1)
	fluid1 = gl_Cell_par_MF(:, :, s1)   ! ��������� ��������� ��������� ��� ������������
	
	! ��������� ������ ��������� ��������������� ��������
    if(norm2(qqq1(2:4))/sqrt(ggg*qqq1(5)/qqq1(1)) > 2.2) then
        rad1 = norm2(gl_Cell_center(:, s1))
        rad2 = norm2(gl_Gran_center(:, gr))
        qqq1(1) = qqq1(1) * rad1**2 / rad2**2
        qqq1(9) = qqq1(9) * rad1**2 / rad2**2
        qqq1(5) = qqq1(5) * rad1**(2 * ggg) / rad2**(2 * ggg)
        ! �������� ������ � ����������� �.�.
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
                fluid2 = gl_Cell_par_MF(:, :, s2)   ! ��������� ��������� ��������� ��� ������������
                dist = min(gl_Cell_dist(s1), gl_Cell_dist(s2))

                ! ��������� ������ ��������� ��������������� ��������
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

            else  ! � ������ ��������� ����� - ��������� �������
                !if (norm2(gl_Cell_center(:, s1)) <= par_R0 * par_R_int) CYCLE
                if(s2 == -1) then  ! ���������� �����
                    dist = gl_Cell_dist(s1)
                    qqq2 = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 1.0_8, 0.0_8, 0.0_8, 0.0_8, 100.0_8/)
                    fluid2(:, 1) = (/0.0001_8, 0.0_8, 0.0_8, 0.0_8, 0.0001_8/)
                    fluid2(:, 2) = (/0.0001_8, 0.0_8, 0.0_8, 0.0_8, 0.0001_8/)
                    fluid2(:, 3) = (/0.0001_8, 0.0_8, 0.0_8, 0.0_8, 0.0001_8/)
                    fluid2(:, 4) = (/1.0_8, par_Velosity_inf, 0.0_8, 0.0_8, 0.5_8/)
                else  ! ����� ����� ������ �������
                    dist = gl_Cell_dist(s1)
                    qqq2 = qqq1
                    fluid2 = fluid1
                    qqq2(5) = 1.0_8
                    if(qqq2(2) > par_Velosity_inf) then
                        qqq2(2) = par_Velosity_inf ! ����� ��������
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
	
	! ----------------------------------------------------------- �������� �� ����� ������� ----------------------
	
	if (time_step > time) then 
		time =  atomicmin(time_step, time)   ! ��������� �������� ������ ������������ ��������
	end if
	
	!if(gr <= 200) then
	!	print*, gr, time, time_step
	!end if
	
	
	end subroutine CUF_GD_3_grans
	
	
	attributes(global) subroutine CUF_GD_3_cells_inner()
	! ����������� �����, ������� �������� + �����������
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
    real(8) :: qqq1(9), qqq2(9), qqq(9)  ! ���������� � ������
    real(8) :: fluid1(5, 4), fluid2(5, 4)
    real(8) :: dist, dsl, dsc, dsp
    real(8) :: POTOK(9), ttest(3)
    real(8) :: POTOK_MF(5)
    real(8) :: POTOK_MF_all(5, 4)
    real(8) :: time, Volume, TT, U8, rad1, rad2, aa, bb, cc
    real(8) :: SOURSE(5,5)  ! ��������� �����, �������� � ������� ��� ������ � ������� ����� ������������
    real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3
	logical :: l_1
	
	iter = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! ����� ������
	time = time_step
	
	if (iter > dev_Ncell_inner) return
	
	gr = gl_all_Cell_inner(iter)
            !if(gl_Cell_info(gr) == 2) CYCLE
            l_1 = .TRUE.
            if ((gl_Cell_type(gr) == "A" .or. gl_Cell_type(gr) == "B").and.(gl_Cell_number(1, gr) <= 2) ) l_1 = .FALSE.    ! �� ������� � ������ ���� �������
            POTOK = 0.0
            SOURSE = 0.0
            POTOK_MF_all = 0.0
            Volume = gl_Cell_Volume(gr)
            qqq = gl_Cell_par(:, gr)
            fluid1 = gl_Cell_par_MF(:, :, gr)
            ! ������������ ������ ����� �����
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


            ! ���������� ���� � ������� ���������

            zone = 1


            call dev_Calc_sourse_MF(qqq, fluid1, SOURSE, zone)  ! ��������� ���������

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

            ! ������ ��������� ������ ���������� ��� ��������� ���������

            do i = 1, 4
                if (i == 1 .and. l_1 == .FALSE.) CYCLE
                
                if (l_1 == .FALSE.) SOURSE(:, i + 1) = 0.0       ! �� ������������ �������� ������ �����
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
	
	! ----------------------------------------------------------- �������� �� ����� ������� ----------------------
	
	end subroutine CUF_GD_3_cells_inner
	
	
	attributes(global) subroutine CUF_GD_3_cells()
	! ����������� �����, ������� �������� + �����������
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
    real(8) :: qqq1(9), qqq2(9), qqq(9)  ! ���������� � ������
    real(8) :: fluid1(5, 4), fluid2(5, 4)
    real(8) :: dist, dsl, dsc, dsp
    real(8) :: POTOK(9), ttest(3)
    real(8) :: POTOK_MF(5)
    real(8) :: POTOK_MF_all(5, 4)
    real(8) :: time, Volume, TT, U8, rad1, rad2, aa, bb, cc
    real(8) :: SOURSE(5,5)  ! ��������� �����, �������� � ������� ��� ������ � ������� ����� ������������
    real(8) :: ro3, u3, v3, w3, p3, bx3, by3, bz3, Q3
	logical :: l_1
	
	gr = blockDim%x * (blockIdx%x - 1) + threadIdx%x   ! ����� ������
	time = time_step
	
	if (gr > dev_Ncell) return
	
	! ----------------------------------------------------------- ����� ������ ����������� ��� �� ����� ----------------------
	
	if(gl_Cell_info(gr) == 0) return
	
    l_1 = .TRUE.
    !if (norm2(gl_Cell_center(:, gr)) <= 1.0 * par_R0)   write(*,*) "NONE ERROR  27465678uhgfdr"  !l_1 = .FALSE.    ! �� ������� ������ �����   ����� ������ �� ������ ����
    POTOK = 0.0
    SOURSE = 0.0
    POTOK_MF_all = 0.0
    Volume = gl_Cell_Volume(gr)
    qqq = gl_Cell_par(:, gr)
    fluid1 = gl_Cell_par_MF(:, :, gr)
    ! ������������ ������ ����� �����
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
 
 
    ! ���������� ���� � ������� ���������
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
 
 
    call dev_Calc_sourse_MF(qqq, fluid1, SOURSE, zone)  ! ��������� ���������
 
 
 
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
 
    ! ������ ��������� ������ ���������� ��� ��������� ���������
 
    do i = 1, 4
        if (i == 1 .and. l_1 == .FALSE.) CYCLE       ! ���������� ���������� ����� ��� ����� 1
        if (l_1 == .FALSE.) SOURSE(:, i + 1) = 0.0       ! ���������� ���������� ����� ��� ����� 1
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
	
	! ----------------------------------------------------------- �������� �� ����� ������� ----------------------
	
	end subroutine CUF_GD_3_cells
	
	
	subroutine CUDA_info()
	! �������, ���������� ���������� � ���������� �� ����� (������ � ������ ����������)
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

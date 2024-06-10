
	module MY_CUDA_smooth
		USE My_func
		contains 
		!@cuf attributes(host, device) & 
		subroutine Smooth_kvadr(x1, y1, z1, x2, y2, z2, x3, y3, z3, x, y, z, xx, yy, zz)
		! Функция определяющая новвое положение точки x, y, z
		! В соответствии с тремя другими точками, строя квадратичный сплайн
		real(8), intent(in) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x, y, z
		real(8), intent(out) :: xx, yy, zz
		real(8) :: ex(3), ey(3), c(3), xx2, yy2, xx3, xx4, a, b, yy4, nn
		
		ex(1) = x3 - x1
		ex(2) = y3 - y1
		ex(3) = z3 - z1
		
		ey(1) = x2 - x1
		ey(2) = y2 - y1
		ey(3) = z2 - z1
		
		nn = sqrt(ex(1)**2 + ex(2)**2 + ex(3)**2)
		ex = ex/nn
		
		nn = sqrt(ey(1)**2 + ey(2)**2 + ey(3)**2)
		ey = ey/nn
		
		!print*, "ex = ", ex
		
		! Если они не сонаправлены
		if( dabs(DOT_PRODUCT(ex, ey)) < 0.99 ) then
			c(1) = ex(2) * ey(3) - ex(3) * ey(2)
			c(2) = ex(3) * ey(1) - ex(1) * ey(3)
			c(3) = ex(1) * ey(2) - ex(2) * ey(1)
			
			ey(1) = c(2) * ex(3) - c(3) * ex(2)
			ey(2) = c(3) * ex(1) - c(1) * ex(3)
			ey(3) = c(1) * ex(2) - c(2) * ex(1)
			
			
			nn = sqrt(ey(1)**2 + ey(2)**2 + ey(3)**2)
			ey = ey/nn
			
			!print*, "ey = ", ey
			
			c(1) = x2 - x1; c(2) = y2 - y1; c(3) = z2 - z1
			xx2 = DOT_PRODUCT(ex, c)
			yy2 = DOT_PRODUCT(ey, c)
			c(1) = x3 - x1; c(2) = y3 - y1; c(3) = z3 - z1
			xx3 = DOT_PRODUCT(ex, c)
			c(1) = x - x1; c(2) = y - y1; c(3) = z - z1
			xx4 = DOT_PRODUCT(ex, c)
			
			a = -yy2/((xx3 - xx2) * xx2)
			b = -xx3 * yy2/((xx2 - xx3) * xx2)
			yy4 = a * xx4**2 + b * xx4
			
			xx = x1 + xx4 * ex(1) + yy4 * ey(1)
			yy = y1 + xx4 * ex(2) + yy4 * ey(2)
			zz = z1 + xx4 * ex(3) + yy4 * ey(3)
			
			return
		else
			xx = x
			yy = y
			zz = z
		end if
		
		
		! Body of Smooth_kvadr
		
		end subroutine Smooth_kvadr
		
		
		!@cuf attributes(host) & 
		subroutine Smooth_kvadr2(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5, xx, yy, zz)
		!@cuf use cublas
		! Функция определяющая новвое положение точки x, y, z
		! В соответствии с тремя другими точками, строя квадратичный сплайн по методу наименьших квадратов
			real(8), intent(in) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5
			real(8), intent(out) :: xx, yy, zz
			real(8) :: ex(3), ey(3), c(3), xx2, yy2, xx3, yy3, xx4, yy4, xx5, nn
		
		Integer, Parameter :: nb = 64
		!     .. Local Scalars ..
		Real(8) :: rcond
		Integer :: i, info, lda, lwork,rank
		Integer, parameter :: n = 3
		Integer, parameter :: m = 5
		!     .. Local Arrays ..
		Real(8) :: a(5, 3), b(5), work(265)
		Integer :: jpvt(3)
		
		
		!$inter external dgelsy
		
		ex(1) = x5 - x1
		ex(2) = y5 - y1
		ex(3) = z5 - z1
		
		ey(1) = x3 - x1
		ey(2) = y3 - y1
		ey(3) = z3 - z1
		
		nn = sqrt(ex(1)**2 + ex(2)**2 + ex(3)**2)
		ex = ex/nn
		
		nn = sqrt(ey(1)**2 + ey(2)**2 + ey(3)**2)
		ey = ey/nn
		
		!print*, "ex = ", ex
		
		! Если они не сонаправлены
		if( dabs(DOT_PRODUCT(ex, ey)) < 0.99 ) then
			c(1) = ex(2) * ey(3) - ex(3) * ey(2)
			c(2) = ex(3) * ey(1) - ex(1) * ey(3)
			c(3) = ex(1) * ey(2) - ex(2) * ey(1)
			
			ey(1) = c(2) * ex(3) - c(3) * ex(2)
			ey(2) = c(3) * ex(1) - c(1) * ex(3)
			ey(3) = c(1) * ex(2) - c(2) * ex(1)
			
			
			nn = sqrt(ey(1)**2 + ey(2)**2 + ey(3)**2)
			ey = ey/nn
			
			!print*, "ey = ", ey
			
			c(1) = x2 - x1; c(2) = y2 - y1; c(3) = z2 - z1
			xx2 = DOT_PRODUCT(ex, c)
			yy2 = DOT_PRODUCT(ey, c)
			
			c(1) = x3 - x1; c(2) = y3 - y1; c(3) = z3 - z1
			xx3 = DOT_PRODUCT(ex, c)
			yy3 = DOT_PRODUCT(ey, c)
			
			c(1) = x4 - x1; c(2) = y4 - y1; c(3) = z4 - z1
			xx4 = DOT_PRODUCT(ex, c)
			yy4 = DOT_PRODUCT(ey, c)
			
			c(1) = x5 - x1; c(2) = y5 - y1; c(3) = z5 - z1
			xx5 = DOT_PRODUCT(ex, c)
			
			!     Skip heading in data file
				
		
				lda = m
				lwork = 3*n + nb*(n+1)

			!     Read A and B from data file

				a(1, 1) = 0
				a(1, 2) = 0
				a(1, 3) = 1
		
				a(2, 1) = xx2**2
				a(2, 2) = xx2
				a(2, 3) = 1
		
				a(3, 1) = xx3**2
				a(3, 2) = xx3
				a(3, 3) = 1
		
				a(4, 1) = xx4**2
				a(4, 2) = xx4
				a(4, 3) = 1
		
				a(5, 1) = xx5**2
				a(5, 2) = xx5
				a(5, 3) = 1.0
		
				b(1) = 0.0
				b(2) = yy2
				b(3) = yy3
				b(4) = yy4
				b(5) = 0.0
		

			!     Initialize JPVT to be zero so that all columns are free

				jpvt(1:n) = 0

			!     Choose RCOND to reflect the relative accuracy of the input data

				rcond = 0.000001_8

			!     Solve the least squares problem min( norm2(b - Ax) ) for the x
			!     of minimum norm.

				!$inter  Call dgelsy(m, n, 1, a, lda, b, m, jpvt, rcond, rank, work, lwork, info)

			!     Print solution
		
				!  print*, "solution 1 = ", b(1:3)
			
				
			yy3 = b(1) * xx3**2 + b(2) * xx3 + b(3)
			
			xx = x1 + xx3 * ex(1) + yy3 * ey(1)
			yy = y1 + xx3 * ex(2) + yy3 * ey(2)
			zz = z1 + xx3 * ex(3) + yy3 * ey(3)
			
			return
		else
			xx = x3
			yy = y3
			zz = z3
			
			return
		end if

		
		end subroutine Smooth_kvadr2
		
		!@cuf attributes(host, device) & 
		subroutine Smooth_kvadr3(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5, xx, yy, zz)
		!@cuf use cublas
		! Функция определяющая новвое положение точки x, y, z
		! В соответствии с тремя другими точками, строя квадратичный сплайн по методу наименьших квадратов
			real(8), intent(in) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5
			real(8), intent(out) :: xx, yy, zz
			real(8) :: ex(3), ey(3), c(3), xx2, yy2, xx3, yy3, xx4, yy4, xx5, nn
			real(8) :: a, b, cc, g2, g3, g4, h2, h4
		
		
		ex(1) = x5 - x1
		ex(2) = y5 - y1
		ex(3) = z5 - z1
		
		ey(1) = x3 - x1
		ey(2) = y3 - y1
		ey(3) = z3 - z1
		
		nn = sqrt(ex(1)**2 + ex(2)**2 + ex(3)**2)
		ex = ex/nn
		
		nn = sqrt(ey(1)**2 + ey(2)**2 + ey(3)**2)
		ey = ey/nn
		
		
		! Если они не сонаправлены
		if( dabs(DOT_PRODUCT(ex, ey)) < 0.99 ) then
			c(1) = ex(2) * ey(3) - ex(3) * ey(2)
			c(2) = ex(3) * ey(1) - ex(1) * ey(3)
			c(3) = ex(1) * ey(2) - ex(2) * ey(1)
			
			ey(1) = c(2) * ex(3) - c(3) * ex(2)
			ey(2) = c(3) * ex(1) - c(1) * ex(3)
			ey(3) = c(1) * ex(2) - c(2) * ex(1)
			
			
			nn = sqrt(ey(1)**2 + ey(2)**2 + ey(3)**2)
			ey = ey/nn
			
			
			c(1) = x2 - x1; c(2) = y2 - y1; c(3) = z2 - z1
			xx2 = DOT_PRODUCT(ex, c)
			yy2 = DOT_PRODUCT(ey, c)
			
			c(1) = x3 - x1; c(2) = y3 - y1; c(3) = z3 - z1
			xx3 = DOT_PRODUCT(ex, c)
			yy3 = DOT_PRODUCT(ey, c)
			
			c(1) = x4 - x1; c(2) = y4 - y1; c(3) = z4 - z1
			xx4 = DOT_PRODUCT(ex, c)
			yy4 = DOT_PRODUCT(ey, c)
			
			c(1) = x5 - x1; c(2) = y5 - y1; c(3) = z5 - z1
			xx5 = DOT_PRODUCT(ex, c)
			
			g2 = xx2/xx5
			g3 = xx3/xx5
			g4 = xx4/xx5
			
			if(dabs(y3) < 0.0000001) GO TO 11
			
			h2 = yy2/yy3
			h4 = yy4/yy3
			
			
			a = ((-1 + g4)*g4* ((-1 + g4)*(1 + h2) + 3*h4) +  g3**3*(h2 + g4*h2 + h4 - 3*g4*h4) +  &
			g2**3*(1 + g4 + g3*(-3 + h4) +  h4 - 3*g4*h4) + g3*(-3 + h2 + h4 + g4*(1 + g4* &
					(1 + g4*(-3 + h2) - 2*h4) + h4)) +   g2*(1 - 3*h2 +  &
				g3*(1 + g4**2)*(1 + h2) + h4 - 2*g3*g4**2*h4 + g3**3*(-3*h2 + h4) + &
				g4*(h2 +  g4* (g4 + h2 - 3*g4*h2 - 2*h4) + h4) +  g3**2* &
				(-2 + h2 + g4*(-2 + h2 + h4)) ) + g3**2* (3 - 2*h2 - 2*h4 +  &
				g4*(-2 + h4 + g4*(3 - 2*h2 + 3*h4))) + g2**2*(-2 + 3*h2 + &
				g3**2*(3 + 3*h2 - 2*h4) - 2*h4 + g3*(1 + g4 - 2*h2 - 2*g4*h2 + g4*h4) + g4*(-2*h2 + h4 + &
					g4*(-2 + 3*h2 + 3*h4))))/ (3*(-1 + g4)**2*g4**2 - &
			2*g3*(-1 + g4)**2*g4*(1 + g4) + 2*g3**3* (-3 + g4 + g4**2 - 3*g4**3) +  &
			g3**4*(3 + g4*(-2 + 3*g4)) + g2**4*(3 + 3*g3**2 -  &
				2*g3*(1 + g4) + g4*(-2 + 3*g4)  ) + 2*g2**3* (-3 - 3*g3**3 + g4 + g4**2 -  &
				3*g4**3 + g3**2*(1 + g4) +  g3*(1 + g4**2)) +  2*g2*(-(g3**4*(1 + g4)) - &
				(-1 + g4)**2*g4*(1 + g4) +g3**3*(1 + g4**2) + g3**2*(1 + g4**3) - &
				g3*(1 + g4**4)) +  g3**2*(3 + g4*(2 + &
					g4*(-6 + g4*(2 + 3*g4)))) + g2**2* (3 + 3*g3**4 +   2*g3**3*(1 + g4) - 6*g3**2*(1 + g4**2) +  &
				2*g3*(1 + g4**3) +  g4*(2 +  g4*(-6 + g4*(2 + 3*g4)))))
			
			b = ((-1 + g4)*g4*   (1 + h2 - g4**2*(1 + h2) - 3*h4)  - g3*(-3 + g4**4*(-3 + h2) + &
				h2 - g4**2*(-2 + h4) + h4) -  g3**4*(h2 + g4*h2 + h4 -	  3*g4*h4) -  g2**4*(1 + g4 + g3*(-3 + h4) + &
				h4 - 3*g4*h4) +  g3**3*(h2 + g4**2*h2 + h4 -  3*g4**2*h4) + g2**3*(1 + g4**2 +  &
				g3**2*(-3 + h4) + h4 - 3*g4**2*h4) +  g3**2*(-3 + h2 + h4 + g4*(1 - 2*h4 +  g4*(1 + g4*(-3 + h2) + h4)) ) + g2* &
			(-1 + 3*h2 +  g3**4*(3*h2 - h4) - h4 +  g4**2* (-2*h2 + g4**2*(-1 + 3*h2) +  h4) +  g3**2*(1 - 2*h2 + &
					g4**2*(1 - 2*h2 + h4))) +  g2**2*(1 +  g3*(1 + g4**2)*(-2 + h2) -  3*h2 + h4 + g3*g4**2*h4 + g3**3*(-3*h2 + h4) + &
				g3**2* (1 + g4 + h2 + g4*h2 -   2*g4*h4) +  g4*(h2 - 2*h4 +  g4*(g4 + h2 - 3*g4*h2 + h4)  )))/ &
			(3*(-1 + g4)**2*g4**2 -  2*g3*(-1 + g4)**2*g4*(1 + g4) + 2*g3**3*(-3 + g4 + g4**2 - 3*g4**3) + &
			g3**4*(3 + g4*(-2 + 3*g4)) + g2**4*(3 + 3*g3**2 - 2*g3*(1 + g4) + g4*(-2 + 3*g4)) + 2*g2**3* &
			(-3 - 3*g3**3 + g4 + g4**2 - 3*g4**3 + g3**2*(1 + g4) +  g3*(1 + g4**2)) +  2*g2*(-(g3**4*(1 + g4)) - &
				(-1 + g4)**2*g4*(1 + g4) +  g3**3*(1 + g4**2) +  g3**2*(1 + g4**3) - g3*(1 + g4**4)) + &
			g3**2*(3 + g4*(2 +  g4*(-6 + g4*(2 + 3*g4)))) + g2**2*(3 + 3*g3**4 + 2*g3**3*(1 + g4) - &
				6*g3**2*(1 + g4**2) + 2*g3*(1 + g4**3) +  g4*(2 +  g4*(-6 + g4*(2 + 3*g4)))))
			
			cc = ((-1 + g4)**2*g4**2*(1 + h2) -  g3*(-1 + g4)*g4* (-1 + g4**2 - h4) +  g3**4*(h2 + g4**2*h2 + h4 -  g4*h4) + &
			g3**3*(-2*(1 + g4**3)*h2 +  (-2 + g4 + g4**2)*h4) +  g3**2*(h2 + h4 +  g4*(1 +  g4* (-2 + g4 + g4**2*h2 - &
					2*h4) + h4)) +  g2**4*(1 + g4**2 + h4 +  g3**2*h4 - g4*h4 -   g3*(1 + g4 + g4*h4)) +  g2*(g3**2*(1 + g4**3)*(1 + h2) - &
				g3*(1 + g4**4)*(1 + h2) + (-1 + g4)*g4*(h2 - g4**2*h2 + h4) - g3**4*(h2 + g4*h2 + g4*h4) + g3**3* &
				(h2 + g4**2*h2 + g4**2*h4)) + g2**3* (-2*(1 + g4**3) - 2*g3**3*h4 +  (-2 + g4 + g4**2)*h4 +   g3**2*(1 + g4 + g4*h4) + &
				g3*(1 + g4**2*(1 + h4))) +  g2**2*(1 +   g3*(1 + g4**3)*(1 + h2) +  h4 + g3**4*h4 +   g3**3*(h2 + g4*h2 + g4*h4) - &
				2*g3**2*  (1 + h2 +  g4**2*(1 + h2 + h4)) +  g4*(h2 + h4 + g4*  (g4*(g4 + h2) -   2*(h2 + h4)))))/ &
			(3*(-1 + g4)**2*g4**2 -  2*g3*(-1 + g4)**2*g4*(1 + g4) + 2*g3**3* (-3 + g4 + g4**2 - 3*g4**3) + &
			g3**4*(3 + g4*(-2 + 3*g4)) +  g2**4*(3 + 3*g3**2 -    2*g3*(1 + g4) + g4*(-2 + 3*g4)  ) + 2*g2**3*  (-3 - 3*g3**3 + g4 + g4**2 - &
				3*g4**3 + g3**2*(1 + g4) +  g3*(1 + g4**2)) +  2*g2*(-(g3**4*(1 + g4)) -   (-1 + g4)**2*g4*(1 + g4) + &
				g3**3*(1 + g4**2) +  g3**2*(1 + g4**3) -  g3*(1 + g4**4)) +  g3**2*(3 +   g4*(2 + &
					g4*(-6 + g4*(2 + 3*g4))))  + g2**2*(3 + 3*g3**4 +   2*g3**3*(1 + g4) -   6*g3**2*(1 + g4**2) + &
				2*g3*(1 + g4**3) +   g4*(2 +    g4*(-6 + g4*(2 + 3*g4)))))
			
				
			yy3 = (a * (xx3/xx5)**2 + b * (xx3/xx5) + cc) * yy3
			
			
			xx = x1 + xx3 * ex(1) + yy3 * ey(1)
			yy = y1 + xx3 * ex(2) + yy3 * ey(2)
			zz = z1 + xx3 * ex(3) + yy3 * ey(3)
			
			return
		else		
			xx = x3
			yy = y3
			zz = z3
			
			return
		end if

		11      continue		
			xx = x3
			yy = y3
			zz = z3
			
			return
		
		end subroutine Smooth_kvadr3
		
		!@cuf attributes(host, device) & 
		subroutine Smooth_kvadr4(x1, y1, z1,  x2, y2, z2,  x3, y3, z3,  x4, y4, z4,  x5, y5, z5,  x6, y6, z6,  x7, y7, z7,  xx, yy, zz)
		! Функция определяющая новвое положение точки x, y, z
		! В соответствии с тремя другими точками, строя квадратичный сплайн по методу наименьших квадратов
			real(8), intent(in) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5, x6, y6, z6, x7, y7, z7
			real(8), intent(out) :: xx, yy, zz
			real(8) :: ex(3), ey(3), c(3), xx2, yy2, xx3, yy3, xx4, yy4, xx5, yy5, xx6, yy6, xx7, nn
			real(8) :: AA(7, 3), AAT(3, 7), BB(7, 1)
			real(8) :: CC(3, 1), DDR(3, 3), EE(3, 1)
			real(8) :: DD(3, 3) 
			integer :: i
		
		
			ex(1) = x7 - x1
			ex(2) = y7 - y1
			ex(3) = z7 - z1
		
			ey(1) = x4 - x1
			ey(2) = y4 - y1
			ey(3) = z4 - z1
		
			nn = sqrt(ex(1)**2 + ex(2)**2 + ex(3)**2)
			ex = ex/nn
		
			nn = sqrt(ey(1)**2 + ey(2)**2 + ey(3)**2)
			ey = ey/nn
		
		!print*, dabs(DOT_PRODUCT(ex, ey))
		! Если они не сонаправлены
		if( dabs(DOT_PRODUCT(ex, ey)) < 0.9999 ) then
			c(1) = ex(2) * ey(3) - ex(3) * ey(2)
			c(2) = ex(3) * ey(1) - ex(1) * ey(3)
			c(3) = ex(1) * ey(2) - ex(2) * ey(1)
			
			ey(1) = c(2) * ex(3) - c(3) * ex(2)
			ey(2) = c(3) * ex(1) - c(1) * ex(3)
			ey(3) = c(1) * ex(2) - c(2) * ex(1)
			
			
			nn = sqrt(ey(1)**2 + ey(2)**2 + ey(3)**2)
			ey = ey/nn
			
			
			c(1) = x2 - x1; c(2) = y2 - y1; c(3) = z2 - z1
			xx2 = DOT_PRODUCT(ex, c)
			yy2 = DOT_PRODUCT(ey, c)
			
			c(1) = x3 - x1; c(2) = y3 - y1; c(3) = z3 - z1
			xx3 = DOT_PRODUCT(ex, c)
			yy3 = DOT_PRODUCT(ey, c)
			
			c(1) = x4 - x1; c(2) = y4 - y1; c(3) = z4 - z1
			xx4 = DOT_PRODUCT(ex, c)
			yy4 = DOT_PRODUCT(ey, c)
			
			c(1) = x5 - x1; c(2) = y5 - y1; c(3) = z5 - z1
			xx5 = DOT_PRODUCT(ex, c)
			yy5 = DOT_PRODUCT(ey, c)
			
			c(1) = x6 - x1; c(2) = y6 - y1; c(3) = z6 - z1
			xx6 = DOT_PRODUCT(ex, c)
			yy6 = DOT_PRODUCT(ey, c)
			
			c(1) = x7 - x1; c(2) = y7 - y1; c(3) = z7 - z1
			xx7 = DOT_PRODUCT(ex, c)
			
			AA(1, 1) = 0
			AA(1, 2) = 0
			AA(1, 3) = 1
		
			AA(2, 1) = xx2**2
			AA(2, 2) = xx2
			AA(2, 3) = 1
		
			AA(3, 1) = xx3**2
			AA(3, 2) = xx3
			AA(3, 3) = 1
		
			AA(4, 1) = xx4**2
			AA(4, 2) = xx4
			AA(4, 3) = 1
		
			AA(5, 1) = xx5**2
			AA(5, 2) = xx5
			AA(5, 3) = 1.0
			
			AA(6, 1) = xx6**2
			AA(6, 2) = xx6
			AA(6, 3) = 1.0
			
			AA(7, 1) = xx7**2
			AA(7, 2) = xx7
			AA(7, 3) = 1.0
		
			BB(1, 1) = 0.0
			BB(2, 1) = yy2
			BB(3, 1) = yy3
			BB(4, 1) = yy4
			BB(5, 1) = yy5
			BB(6, 1) = yy6
			BB(7, 1) = 0.0
			
			do i = 1, 7
				AAT(:, i) = AA(i, :)
			end do
			
			CC = MATMUL(AAT, BB)
			
			DD = MATMUL(AAT, AA)
			
			call Smooth_inverse(DD, DDR)
			
			EE = MATMUL(DDR, CC)
			
			!print*, "EE = ", EE
				
			yy4 = (EE(1, 1) * xx4**2 + EE(2, 1) * xx4 + EE(3, 1))
			
			
			xx = x1 + xx4 * ex(1) + yy4 * ey(1)
			yy = y1 + xx4 * ex(2) + yy4 * ey(2)
			zz = z1 + xx4 * ex(3) + yy4 * ey(3)
			
			return
		else		
			xx = x4
			yy = y4
			zz = z4
			
			return
		end if

		11      continue		
			xx = x4
			yy = y4
			zz = z4
			
			return
		
		end subroutine Smooth_kvadr4
		
		
		!@cuf attributes(host, device) & 
		subroutine Smooth_kvadr5(x1, y1, z1,  x2, y2, z2,  x3, y3, z3,  x4, y4, z4,  x5, y5, z5,  x6, y6, z6,  x7, y7, z7, &
			x8, y8, z8, x9, y9, z9, xx, yy, zz)
		! Функция определяющая новвое положение точки x, y, z
		! В соответствии с тремя другими точками, строя квадратичный сплайн по методу наименьших квадратов
			real(8), intent(in) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5, x6, y6, z6, x7, y7, z7, x8, y8, z8, x9, y9, z9
			real(8), intent(out) :: xx, yy, zz
			real(8) :: ex(3), ey(3), c(3), xx2, yy2, xx3, yy3, xx4, yy4, xx5, yy5, xx6, yy6, xx7, yy7, xx8, yy8, xx9, nn
			real(8) :: AA(9, 3), AAT(3, 9), BB(9, 1)
			real(8) :: CC(3, 1), DDR(3, 3), EE(3, 1)
			real(8) :: DD(3, 3) 
			integer :: i
		
		
			ex(1) = x9 - x1
			ex(2) = y9 - y1
			ex(3) = z9 - z1
		
			ey(1) = x5 - x1
			ey(2) = y5 - y1
			ey(3) = z5 - z1
		
			nn = sqrt(ex(1)**2 + ex(2)**2 + ex(3)**2)
			ex = ex/nn
		
			nn = sqrt(ey(1)**2 + ey(2)**2 + ey(3)**2)
			ey = ey/nn
		
		!print*, dabs(DOT_PRODUCT(ex, ey))
		! Если они не сонаправлены
		if( dabs(DOT_PRODUCT(ex, ey)) < 0.9999 ) then
			c(1) = ex(2) * ey(3) - ex(3) * ey(2)
			c(2) = ex(3) * ey(1) - ex(1) * ey(3)
			c(3) = ex(1) * ey(2) - ex(2) * ey(1)
			
			ey(1) = c(2) * ex(3) - c(3) * ex(2)
			ey(2) = c(3) * ex(1) - c(1) * ex(3)
			ey(3) = c(1) * ex(2) - c(2) * ex(1)
			
			
			nn = sqrt(ey(1)**2 + ey(2)**2 + ey(3)**2)
			ey = ey/nn
			
			
			c(1) = x2 - x1; c(2) = y2 - y1; c(3) = z2 - z1
			xx2 = DOT_PRODUCT(ex, c)
			yy2 = DOT_PRODUCT(ey, c)
			
			c(1) = x3 - x1; c(2) = y3 - y1; c(3) = z3 - z1
			xx3 = DOT_PRODUCT(ex, c)
			yy3 = DOT_PRODUCT(ey, c)
			
			c(1) = x4 - x1; c(2) = y4 - y1; c(3) = z4 - z1
			xx4 = DOT_PRODUCT(ex, c)
			yy4 = DOT_PRODUCT(ey, c)
			
			c(1) = x5 - x1; c(2) = y5 - y1; c(3) = z5 - z1
			xx5 = DOT_PRODUCT(ex, c)
			yy5 = DOT_PRODUCT(ey, c)
			
			c(1) = x6 - x1; c(2) = y6 - y1; c(3) = z6 - z1
			xx6 = DOT_PRODUCT(ex, c)
			yy6 = DOT_PRODUCT(ey, c)
			
			c(1) = x7 - x1; c(2) = y7 - y1; c(3) = z7 - z1
			xx7 = DOT_PRODUCT(ex, c)
			yy7 = DOT_PRODUCT(ey, c)
			
			c(1) = x8 - x1; c(2) = y8 - y1; c(3) = z8 - z1
			xx8 = DOT_PRODUCT(ex, c)
			yy8 = DOT_PRODUCT(ey, c)
			
			c(1) = x9 - x1; c(2) = y9 - y1; c(3) = z9 - z1
			xx9 = DOT_PRODUCT(ex, c)
			
			AA(1, 1) = 0
			AA(1, 2) = 0
			AA(1, 3) = 1
		
			AA(2, 1) = xx2**2
			AA(2, 2) = xx2
			AA(2, 3) = 1
		
			AA(3, 1) = xx3**2
			AA(3, 2) = xx3
			AA(3, 3) = 1
		
			AA(4, 1) = xx4**2
			AA(4, 2) = xx4
			AA(4, 3) = 1
		
			AA(5, 1) = xx5**2
			AA(5, 2) = xx5
			AA(5, 3) = 1.0
			
			AA(6, 1) = xx6**2
			AA(6, 2) = xx6
			AA(6, 3) = 1.0
			
			AA(7, 1) = xx7**2
			AA(7, 2) = xx7
			AA(7, 3) = 1.0
			
			AA(8, 1) = xx8**2
			AA(8, 2) = xx8
			AA(8, 3) = 1.0
			
			AA(9, 1) = xx9**2
			AA(9, 2) = xx9
			AA(9, 3) = 1.0
		
			BB(1, 1) = 0.0
			BB(2, 1) = yy2
			BB(3, 1) = yy3
			BB(4, 1) = yy4
			BB(5, 1) = yy5
			BB(6, 1) = yy6
			BB(7, 1) = yy7
			BB(8, 1) = yy8
			BB(9, 1) = 0.0
			
			do i = 1, 9
				AAT(:, i) = AA(i, :)
			end do
			
			CC = MATMUL(AAT, BB)
			
			DD = MATMUL(AAT, AA)
			
			call Smooth_inverse(DD, DDR)
			
			EE = MATMUL(DDR, CC)
			
			!print*, "EE = ", EE
				
			yy5 = (EE(1, 1) * xx5**2 + EE(2, 1) * xx5 + EE(3, 1))
			
			
			xx = x1 + xx5 * ex(1) + yy5 * ey(1)
			yy = y1 + xx5 * ex(2) + yy5 * ey(2)
			zz = z1 + xx5 * ex(3) + yy5 * ey(3)
			
			return
		else		
			xx = x5
			yy = y5
			zz = z5
			
			return
		end if

		11      continue		
			xx = x5
			yy = y5
			zz = z5
			
			return
		
		end subroutine Smooth_kvadr5
		
		!@cuf attributes(host, device) & 
		subroutine Smooth_inverse(a, c)
		!============================================================
		! Inverse matrix
		! Method: Based on Doolittle LU factorization for Ax=b
		! Alex G. December 2009
		!-----------------------------------------------------------
		! input ...
		! a(n,n) - array of coefficients for matrix A
		! n      - dimension
		! output ...
		! c(n,n) - inverse matrix of A
		! comments ...
		! the original matrix a(n,n) will be destroyed 
		! during the calculation
		!===========================================================
		implicit none 
		real(8), intent(in out) :: a(3,3)
		real(8), intent(out) :: c(3,3)
		integer, parameter :: n = 3
		real(8) :: L(n,n), U(n,n), b(n), d(n), x(n)
		real(8) :: coeff
		integer :: i, j, k

		! step 0: initialization for matrices L and U and b
		! Fortran 90/95 aloows such operations on matrices
		L=0.0
		U=0.0
		b=0.0

		! step 1: forward elimination
		do k=1, n-1
		do i=k+1,n
			coeff = a(i,k) / a(k,k)
			L(i,k) = coeff
			do j=k+1,n
				a(i,j) = a(i,j) - coeff * a(k,j)
			end do
		end do
		end do

		! Step 2: prepare L and U matrices 
		! L matrix is a matrix of the elimination coefficient
		! + the diagonal elements are 1.0
		do i=1,n
		L(i,i) = 1.0
		end do
		! U matrix is the upper triangular part of A
		do j=1,n
		do i=1,j
			U(i,j) = a(i,j)
		end do
		end do

		! Step 3: compute columns of the inverse matrix C
		do k=1,n
		b(k)=1.0
		d(1) = b(1)
		! Step 3a: Solve Ld=b using the forward substitution
		do i=2,n
			d(i)=b(i)
			do j=1,i-1
			d(i) = d(i) - L(i,j)*d(j)
			end do
		end do
		! Step 3b: Solve Ux=d using the back substitution
		x(n)=d(n)/U(n,n)
		do i = n-1,1,-1
			x(i) = d(i)
			do j=n,i+1,-1
			x(i)=x(i)-U(i,j)*x(j)
			end do
			x(i) = x(i)/u(i,i)
		end do
		! Step 3c: fill the solutions x(n) into column k of C
		do i=1,n
			c(i,k) = x(i)
		end do
		b(k)=0.0
		end do
		
		end subroutine Smooth_inverse
		
		
		
		! Добавил в этот модуль движение сетки
		
		!@cuf attributes(host, device) & 
		real(8) pure function Setka_A(i, R_TS, R_HP, R_BS, the, dk13, par_R0, par_R_inner, par_n_IB, par_kk1, par_kk12, &
			par_n_TS, par_n_HP, par_n_BS, par_kk14, par_R_END, par_kk2, par_n_END)
			! Variables
			implicit none
			real(8), intent(in) :: R_TS, R_HP, R_BS, dk13, the, par_kk1, par_kk12, par_R0, par_R_inner, par_kk14, par_R_END, par_kk2
			integer, intent(in) :: i, par_n_IB, par_n_TS, par_n_HP, par_n_BS, par_n_END
			real(8) :: r, rr, rrr, dr, ddr
			
			if (i <= par_n_IB) then  ! NEW
					if(i == 2) then
						r =  par_R0 - (par_R_inner - par_R0) * (DBLE(3 - 2)/(par_n_IB - 2))**par_kk1
					else
						r =  par_R0 + (par_R_inner - par_R0) * (DBLE(i - 2)/(par_n_IB - 2))**par_kk1
					end if
			else if (i <= par_n_TS) then  
				r =  par_R_inner + (R_TS - par_R_inner) * sgushenie_3( (DBLE(i - par_n_IB)/(par_n_TS - par_n_IB)) , par_kk12)
			else if (i <= par_n_HP) then  
				r = R_TS + (R_HP - R_TS) * sgushenie_2(DBLE(i - par_n_TS)/(par_n_HP - par_n_TS), par_kk14)
			else if (i <= par_n_BS) then 
				ddr = R_HP - (R_TS + (R_HP - R_TS) * sgushenie_2(DBLE(par_n_HP - 1 - par_n_TS)/(par_n_HP - par_n_TS), par_kk14))
				r = R_HP + max((R_BS - R_HP) * (DBLE(i - par_n_HP)/(par_n_BS - par_n_HP))**1.6, &
						ddr * (i - par_n_HP))
			else
				ddr = R_HP - (R_TS + (R_HP - R_TS) * sgushenie_2(DBLE(par_n_HP - 1 - par_n_TS)/(par_n_HP - par_n_TS), par_kk14))
				dr = R_BS - (max(R_HP + (R_BS - R_HP) * (DBLE(par_n_BS - 1 - par_n_HP)/(par_n_BS - par_n_HP))**1.6, &
						ddr * (i - par_n_HP)))
				r = R_BS + max((par_R_END - R_BS) * (DBLE(i- par_n_BS)/(par_n_END - par_n_BS))**3.0, &
								dr * (i - par_n_BS))
			end if
			
			Setka_A = r
			return
		
		end function Setka_A
			
			!@cuf attributes(host, device) & 
		real(8) pure function Setka_C(i, R_BS_, dk13, par_n_HP, par_n_BS, par_kk2, par_R_END, N1, rr, dr, ddr, rr0)
			! Variables
			! ddr - расстояние крайней ячеки до BS изнутри
			! rr0 - высота HP на крайнем A-луче
			implicit none
			real(8), intent(in) :: R_BS_, dk13, par_kk2, par_R_END, rr, dr, ddr, rr0
			integer, intent(in) :: i, par_n_HP, par_n_BS, N1
			real(8) :: r, rrr, r1, R_BS

			!R_BS = R_BS_
			R_BS = min(rr + (R_BS_ - rr0), 0.92 * par_R_END)
			!R_BS = max(R_BS, dr * (par_n_BS - par_n_HP + 1) )
			
			!if (i <= 11) then
				!r = rr + dk13 * (i - 1)
			!else if (i <= par_n_BS - par_n_HP + 1) then
			! if (i == 2) then
			! 	r = rr + dr
			! else 
			if (i <= par_n_BS - par_n_HP + 1) then
				!rrr = rr + dk13 * (10)
				!r1 = log((1.5 * dk13)/(R_BS - rr))/log(DBLE(1.0)/(par_n_BS - par_n_HP - 10))
				!r = rrr + (R_BS - rrr) * (DBLE(i - 11)/(par_n_BS - par_n_HP - 10))**r1
				r = rr + max((R_BS - rr) * (DBLE(i - 1)/(par_n_BS - par_n_HP))**1.6, &
							dr * (i - 1.0)/( (1.0 * (i - 1))**(0.2) ))
			else
				!ddr = R_BS - (rr + dr +  (R_BS - rr - dr) * (DBLE(par_n_BS - par_n_HP - 2)/(par_n_BS - par_n_HP - 1))**1.5)
				
				r = max(R_BS + (DBLE(i - (par_n_BS - par_n_HP + 1))/(N1 - (par_n_BS - par_n_HP + 1) ))**3.0 * (par_R_END - R_BS), &
						R_BS + ddr * (i - (par_n_BS - par_n_HP + 1)) )

				if(r > par_R_END) then
					r = par_R_END - ddr + i * ddr/300
				end if
			end if
			
			Setka_C = r
			return
		
		end function Setka_C
		
		!@cuf attributes(host, device) &
		real(8) pure function Setka_O(i, r1, R_HP, R_BS_, dk13, par_n_HP, par_n_BS, par_kk2, par_R_END, N1, ddr, rr0, dddr)
			! Variables
			! dddr - высота крайней ячейки перед HP
			! rr0 - высота HP на крайнем A-луче
			implicit none
			real(8), intent(in) :: R_BS_, dk13, par_kk2, par_R_END, r1, R_HP, ddr, rr0, dddr
			integer, intent(in) :: i, par_n_HP, par_n_BS, N1
			real(8) :: r, rr, R_BS

			R_BS = min(R_HP + (R_BS_ - rr0), 0.92 * par_R_END)
			
			!if (i <= 10) then
			!    r = R_HP + (i - 1) * dk13 * (R_BS - R_HP)
			!else if (i <= par_n_BS - par_n_HP + 1) then
			if (i <= par_n_BS - par_n_HP + 1) then
				!rr = R_HP + 9 * dk13 * (R_BS - R_HP)
				!r = rr + (R_BS - rr) * (DBLE(i - 10)/(par_n_BS - par_n_HP - 9))**r1

				!r = R_HP + (R_BS - R_HP) * (DBLE(i - 1)/(par_n_BS - par_n_HP))**1.6
				r = max(R_HP + (R_BS - R_HP) * (DBLE(i - 1)/(par_n_BS - par_n_HP))**1.6,&
					R_HP + dddr * (i - 1.0)/( (1.0 * (i))**(0.2)))
				!print*, r, R_HP, dddr, i
				!pause
			else
				!r = R_BS + (DBLE(i - (par_n_BS - par_n_HP + 1))/(N1 - (par_n_BS - par_n_HP + 1) ))**(0.55 * par_kk2) * (par_R_END - R_BS)
				r = max(R_BS + (DBLE(i - (par_n_BS - par_n_HP + 1))/(N1 - (par_n_BS - par_n_HP + 1) ))**3.0 * (par_R_END - R_BS), &
						R_BS + ddr * (i - (par_n_BS - par_n_HP + 1)))
			end if
			
			Setka_O = r
			return
		
		end function Setka_O
	
	end module MY_CUDA_smooth
    
    subroutine Find_Surface()   ! Заполняет массивы поверхностей, которые выделяются
		! Для BS находятся только поверхности на полусфере! (x > 0)
		use STORAGE
		use GEO_PARAM
		implicit none
		integer :: j, k, node, num, i

		! HP
		node = 1

		do k = 1, size( gl_Cell_A(par_n_HP - 1, 1, :) )
			do j = 1, size( gl_Cell_A(par_n_HP - 1, :, 1) )
				gl_Contact(node) = gl_Cell_gran(1, gl_Cell_A(par_n_HP - 1, j, k))
				gl_Gran_type(gl_Contact(node)) = 2
				node = node + 1
				

			end do
		end do

		do k = 1, size( gl_Cell_C(par_n_HP - par_n_TS, 1, :) )
			do j = 1, size( gl_Cell_C(par_n_HP - par_n_TS, :, 1) )
				gl_Contact(node) = gl_Cell_gran(1, gl_Cell_C(par_n_HP - par_n_TS, j, k))
				gl_Gran_type(gl_Contact(node)) = 2
				node = node + 1
				
			end do
		end do

		node = 1
		! TS
		do k = 1, size( gl_Cell_A(par_n_TS - 1, 1, :) )
			do j = 1, size( gl_Cell_A(par_n_TS - 1, :, 1) )
				gl_TS(node) = gl_Cell_gran(1, gl_Cell_A(par_n_TS - 1, j, k))
				gl_Gran_type(gl_TS(node)) = 1
				node = node + 1
				
			end do
		end do

		do k = 1, size( gl_Cell_B(par_n_TS - 1, 1, :) )
			do j = 1, size( gl_Cell_B(par_n_TS - 1, :, 1) )
				gl_TS(node) = gl_Cell_gran(1, gl_Cell_B(par_n_TS - 1, j, k))
				gl_Gran_type(gl_TS(node)) = 1
				node = node + 1
				
				
			end do
		end do


		node = 1
		! BS
		num = par_m_A - 1
		do k = 1, size( gl_Cell_A(par_n_BS - 1, 1, :) )
			do j = 1, num
				gl_BS(node) = gl_Cell_gran(1, gl_Cell_A(par_n_BS - 1, j, k))
				gl_Gran_type(gl_BS(node)) = 3
				node = node + 1
			end do
		end do

		if (par_developer_info) print *, "Find_Surface: Contact", node - 1, size(gl_Contact)
		if (par_developer_info) print *, "Find_Surface: TS", node - 1, size(gl_TS)
		if (par_developer_info) print *, "Find_Surface: BS", node - 1, size(gl_BS)
		
		
		
		! Определим для каждой ячейки её зону
		do k = 1, size( gl_Cell_A(1, 1, :) )
			do j = 1, size( gl_Cell_A(1, :, 1) )
				do i = 1, par_n_TS - 1
					num = gl_Cell_A(i, j, k)
					gl_zone_Cell(num) = 1
				end do
			end do
		end do
		
		do k = 1, size( gl_Cell_A(1, 1, :) )
			do j = 1, size( gl_Cell_A(1, :, 1) )
				do i = par_n_HP, par_n_BS - 1
					num = gl_Cell_A(i, j, k)
					gl_zone_Cell(num) = 3
				end do
			end do
		end do
		
		do k = 1, size( gl_Cell_A(1, 1, :) )
			do j = 1, size( gl_Cell_A(1, :, 1) )
				do i = par_n_BS, size( gl_Cell_A(:, 1, 1) )
					num = gl_Cell_A(i, j, k)
					gl_zone_Cell(num) = 4
				end do
			end do
		end do
		
		do k = 1, size( gl_Cell_B(1, 1, :) )
			do j = 1, size( gl_Cell_B(1, :, 1) )
				do i = 1, par_n_TS - 1
					num = gl_Cell_B(i, j, k)
					gl_zone_Cell(num) = 1
				end do
			end do
		end do
		
		
		do k = 1, size( gl_Cell_C(1, 1, :) )
			do j = 1, size( gl_Cell_C(1, :, 1) )
				do i = par_n_HP - par_n_TS + 1, par_n_BS - par_n_TS
					num = gl_Cell_C(i, j, k)
					gl_zone_Cell(num) = 3
				end do
			end do
		end do
		
		do k = 1, size( gl_Cell_C(1, 1, :) )
			do j = 1, size( gl_Cell_C(1, :, 1) )
				do i = par_n_BS - par_n_TS + 1, size( gl_Cell_C(:, 1, 1) )
					num = gl_Cell_C(i, j, k)
					gl_zone_Cell(num) = 4
				end do
			end do
		end do
		
		gl_Gran_scheme = 0

		! Считаем HLLD везде кроме 2 области
		
		do k = 1, size(gl_Gran_scheme(:))
			
			if(gl_Gran_neighbour(1, k) > 0 .and. gl_Gran_neighbour(2, k) > 0) then
			
				if(gl_zone_Cell(gl_Gran_neighbour(1, k)) == 1 .and. gl_zone_Cell(gl_Gran_neighbour(2, k)) == 1) then
					gl_Gran_scheme(k) = 3
				end if
			
				if(gl_zone_Cell(gl_Gran_neighbour(1, k)) == 4 .and. gl_zone_Cell(gl_Gran_neighbour(2, k)) == 4) then
					gl_Gran_scheme(k) = 1!3
				end if
			
				if(gl_zone_Cell(gl_Gran_neighbour(1, k)) == 3 .and. gl_zone_Cell(gl_Gran_neighbour(2, k)) == 3) then
					gl_Gran_scheme(k) = 0!3
				end if
				!
				if(gl_zone_Cell(gl_Gran_neighbour(1, k)) == 2 .and. gl_zone_Cell(gl_Gran_neighbour(2, k)) == 2) then
					gl_Gran_scheme(k) = 0!3
				end if
			end if
			
			
			! if(gl_Gran_center(1, k) < -40) then
			! 	gl_Gran_scheme(k) = 3
			! end if
			
		end do
		
		
		! HP
		node = 1

		
		do k = 1, size( gl_Cell_A(par_n_HP - 1, 1, :) )
			do j = 1, size( gl_Cell_A(par_n_HP - 1, :, 1) )
				gl_Contact(node) = gl_Cell_gran(1, gl_Cell_A(par_n_HP - 1, j, k))
				gl_Gran_type(gl_Contact(node)) = 2
				node = node + 1
				
				! Вводим лакса поперёк разрыва 
				! gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_A(par_n_HP - 1, j, k))) = 1!1
				! if(gl_Cell_gran(4, gl_Cell_A(par_n_HP - 1, j, k)) /= 0) then
				! gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_A(par_n_HP - 1, j, k)) ) = 1!1
				! end if
				! gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_A(par_n_HP - 1, j, k))) = 1!1
				! gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_A(par_n_HP - 1, j, k))) = 1!1
				! gl_Gran_scheme(gl_Cell_gran(2, gl_Cell_A(par_n_HP - 1, j, k))) = 1!1
				
				
				! gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_A(par_n_HP - 2, j, k))) = 1!1
				! if(gl_Cell_gran(4, gl_Cell_A(par_n_HP - 2, j, k)) /= 0) then
				! gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_A(par_n_HP - 2, j, k)) ) = 1!1
				! end if
				! gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_A(par_n_HP - 2, j, k))) = 1!1
				! gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_A(par_n_HP - 2, j, k))) = 1!1
				! gl_Gran_scheme(gl_Cell_gran(2, gl_Cell_A(par_n_HP - 2, j, k))) = 1!1
				
				! gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_A(par_n_HP - 3, j, k))) = 1
				! if(gl_Cell_gran(4, gl_Cell_A(par_n_HP - 3, j, k)) /= 0) then
				! gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_A(par_n_HP - 3, j, k)) ) = 1
				! end if
				! gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_A(par_n_HP - 3, j, k))) = 1
				! gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_A(par_n_HP - 3, j, k))) = 1
				! gl_Gran_scheme(gl_Cell_gran(2, gl_Cell_A(par_n_HP - 3, j, k))) = 1
				
				
				! gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_A(par_n_HP, j, k))) = 1!1
				! if(gl_Cell_gran(4, gl_Cell_A(par_n_HP, j, k)) /= 0) then
				! gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_A(par_n_HP, j, k))) = 1!1
				! end if
				! gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_A(par_n_HP, j, k))) = 1!1
				! gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_A(par_n_HP, j, k))) = 1!1
				! gl_Gran_scheme(gl_Cell_gran(1, gl_Cell_A(par_n_HP, j, k))) = 1!1
				
				! gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_A(par_n_HP + 1, j, k))) = 1!1
				! if(gl_Cell_gran(4, gl_Cell_A(par_n_HP + 1, j, k)) /= 0) then
				! gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_A(par_n_HP + 1, j, k))) = 1!1
				! end if
				! gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_A(par_n_HP + 1, j, k))) = 1!1
				! gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_A(par_n_HP + 1, j, k))) = 1!1
				! gl_Gran_scheme(gl_Cell_gran(1, gl_Cell_A(par_n_HP + 1, j, k))) = 1!1
				
				! gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_A(par_n_HP + 2, j, k))) = 1!1
				! if(gl_Cell_gran(4, gl_Cell_A(par_n_HP + 2, j, k)) /= 0) then
				! gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_A(par_n_HP + 2, j, k))) = 1!1
				! end if
				! gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_A(par_n_HP + 2, j, k))) = 1!1
				! gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_A(par_n_HP + 2, j, k))) = 1!1
				! gl_Gran_scheme(gl_Cell_gran(1, gl_Cell_A(par_n_HP + 2, j, k))) = 1!1
				
			end do
		end do

		
		do k = 1, size( gl_Cell_C(par_n_HP - par_n_TS, 1, :) )
			do j = 1, size( gl_Cell_C(par_n_HP - par_n_TS, :, 1) )
				gl_Contact(node) = gl_Cell_gran(1, gl_Cell_C(par_n_HP - par_n_TS, j, k))
				gl_Gran_type(gl_Contact(node)) = 2
				node = node + 1
				
				! Вводим лакса поперёк разрыва 
				! gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_C(par_n_HP - par_n_TS, j, k))) = 1!1
				! if(gl_Cell_gran(4, gl_Cell_C(par_n_HP - par_n_TS, j, k)) /= 0) then
				! gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_C(par_n_HP - par_n_TS, j, k)) ) = 1!1
				! end if
				! gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_C(par_n_HP - par_n_TS, j, k))) = 1!1
				! gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_C(par_n_HP - par_n_TS, j, k))) = 1!1
				! gl_Gran_scheme(gl_Cell_gran(1, gl_Cell_C(par_n_HP - par_n_TS, j, k))) = 1!1
				
				! gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_C(par_n_HP - par_n_TS - 1, j, k))) = 1
				! if(gl_Cell_gran(4, gl_Cell_C(par_n_HP - par_n_TS - 1, j, k)) /= 0) then
				! gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_C(par_n_HP - par_n_TS - 1, j, k)) ) = 1
				! end if
				! gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_C(par_n_HP - par_n_TS - 1, j, k))) = 1
				! gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_C(par_n_HP - par_n_TS - 1, j, k))) = 1
				! gl_Gran_scheme(gl_Cell_gran(1, gl_Cell_C(par_n_HP - par_n_TS - 1, j, k))) = 1
				
				
				! gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_C(par_n_HP - par_n_TS + 1, j, k))) = 1!1
				! if(gl_Cell_gran(4, gl_Cell_C(par_n_HP - par_n_TS + 1, j, k)) /= 0) then
				! gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_C(par_n_HP - par_n_TS + 1, j, k))) = 1!1
				! end if
				! gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_C(par_n_HP - par_n_TS + 1, j, k))) = 1!1
				! gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_C(par_n_HP - par_n_TS + 1, j, k))) = 1!1
				! gl_Gran_scheme(gl_Cell_gran(1, gl_Cell_C(par_n_HP - par_n_TS + 1, j, k))) = 1!1
				
				! gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_C(par_n_HP - par_n_TS + 2, j, k))) = 1
				! if(gl_Cell_gran(4, gl_Cell_C(par_n_HP - par_n_TS + 2, j, k)) /= 0) then
				! gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_C(par_n_HP - par_n_TS + 2, j, k))) = 1
				! end if
				! gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_C(par_n_HP - par_n_TS + 2, j, k))) = 1
				! gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_C(par_n_HP - par_n_TS + 2, j, k))) = 1
				! gl_Gran_scheme(gl_Cell_gran(1, gl_Cell_C(par_n_HP - par_n_TS + 2, j, k))) = 1
			end do
		end do
		
		node = 1
		! TS
		do k = 1, size( gl_Cell_A(par_n_TS - 1, 1, :) )
			do j = 1, size( gl_Cell_A(par_n_TS - 1, :, 1) )
				gl_TS(node) = gl_Cell_gran(1, gl_Cell_A(par_n_TS - 1, j, k))
				gl_Gran_type(gl_TS(node)) = 1
				node = node + 1
				
				! Вводим лакса поперёк разрыва 
				!gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_A(par_n_TS - 1, j, k))) = 0!1
				!if(gl_Cell_gran(4, gl_Cell_A(par_n_TS - 1, j, k)) /= 0) then
				!    gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_A(par_n_TS - 1, j, k)) ) = 0!1
				!end if
				!gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_A(par_n_TS - 1, j, k))) = 0!1
				!gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_A(par_n_TS - 1, j, k))) = 0!1
				!
				!!
				!gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_A(par_n_TS, j, k))) = 0!1
				!if(gl_Cell_gran(4, gl_Cell_A(par_n_TS, j, k)) /= 0) then
				!    gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_A(par_n_TS, j, k))) = 0!1
				!end if
				!gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_A(par_n_TS, j, k))) = 0!1
				!gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_A(par_n_TS, j, k))) = 0!1
				
				
				!gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_A(par_n_TS - 2, j, k))) = 1
				!if(gl_Cell_gran(4, gl_Cell_A(par_n_TS - 2, j, k)) /= 0) then
				!    gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_A(par_n_TS - 2, j, k)) ) = 1
				!end if
				!gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_A(par_n_TS - 2, j, k))) = 1
				!gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_A(par_n_TS - 2, j, k))) = 1
				
				!gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_A(par_n_TS + 1, j, k))) = 0
				!if(gl_Cell_gran(4, gl_Cell_A(par_n_TS + 1, j, k)) /= 0) then
				!    gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_A(par_n_TS + 1, j, k))) = 0
				!end if
				!gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_A(par_n_TS + 1, j, k))) = 0
				!gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_A(par_n_TS + 1, j, k))) = 0
				
				!gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_A(par_n_TS + 2, j, k))) = 0
				!if(gl_Cell_gran(4, gl_Cell_A(par_n_TS + 2, j, k)) /= 0) then
				!    gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_A(par_n_TS + 2, j, k))) = 0
				!end if
				!gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_A(par_n_TS + 2, j, k))) = 0
				!gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_A(par_n_TS + 2, j, k))) = 0
				!
				!gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_A(par_n_TS + 3, j, k))) = 0
				!if(gl_Cell_gran(4, gl_Cell_A(par_n_TS + 3, j, k)) /= 0) then
				!    gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_A(par_n_TS + 3, j, k))) = 0
				!end if
				!gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_A(par_n_TS + 3, j, k))) = 0
				!gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_A(par_n_TS + 3, j, k))) = 0
				
			end do
		end do

		do k = 1, size( gl_Cell_B(par_n_TS - 1, 1, :) )
			do j = 1, size( gl_Cell_B(par_n_TS - 1, :, 1) )
				gl_TS(node) = gl_Cell_gran(1, gl_Cell_B(par_n_TS - 1, j, k))
				gl_Gran_type(gl_TS(node)) = 1
				node = node + 1
				
				! Вводим лакса поперёк разрыва 
				!gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_B(par_n_TS - 1, j, k))) = 1
				!if(gl_Cell_gran(4, gl_Cell_B(par_n_TS - 1, j, k)) /= 0) then
				!    gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_B(par_n_TS - 1, j, k)) ) = 1
				!end if
				!gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_B(par_n_TS - 1, j, k))) = 1
				!gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_B(par_n_TS - 1, j, k))) = 1
				!
				!gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_B(par_n_TS, j, k))) = 1
				!if(gl_Cell_gran(4, gl_Cell_B(par_n_TS, j, k)) /= 0) then
				!    gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_B(par_n_TS, j, k))) = 1
				!end if
				!gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_B(par_n_TS, j, k))) = 1
				!gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_B(par_n_TS, j, k))) = 1
				
			end do
		end do


		node = 1
		! BS
		num = par_m_A - 1
		do k = 1, size( gl_Cell_A(par_n_BS - 1, 1, :) )
			do j = 1, num
				gl_BS(node) = gl_Cell_gran(1, gl_Cell_A(par_n_BS - 1, j, k))
				gl_Gran_type(gl_BS(node)) = 3
				node = node + 1
			end do
		end do

		if (par_developer_info) print *, "Find_Surface: Contact", node - 1, size(gl_Contact)
		if (par_developer_info) print *, "Find_Surface: TS", node - 1, size(gl_TS)
		if (par_developer_info) print *, "Find_Surface: BS", node - 1, size(gl_BS)
		
		! Вводим Лакса на оси симметрии
		if(.False.) then
			do k = 1, size( gl_Cell_A(1, 1, :) )
				do j = par_n_HP + 10, size( gl_Cell_A(:, 1, 1) )
					!Вводим лакса поперёк
					gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_A(j, 1, k))) = 1!1
					if(gl_Cell_gran(4, gl_Cell_A(j, 1, k)) /= 0) then
						gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_A(j, 1, k))) = 1!1
					end if
					gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_A(j, 1, k))) = 1!1
					gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_A(j, 1, k))) = 1!1
				
					gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_A(j, 2, k))) = 1!1
					if(gl_Cell_gran(4, gl_Cell_A(j, 2, k)) /= 0) then
						gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_A(j, 2, k))) = 1!1
					end if
					gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_A(j, 2, k))) = 1!1
					gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_A(j, 2, k))) = 1!1
				
				end do
			end do
		
			do k = 1, size( gl_Cell_B(1, 1, :) )
				do j = 1, par_n_TS + 1 !size( gl_Cell_B(:, 1, 1) )
					! Вводим лакса поперёк
					gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_B(j, 1, k))) = 1!1
					if (gl_Cell_gran(4, gl_Cell_B(j, 1, k)) /= 0) then
						gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_B(j, 1, k))) = 1!1
					end if
					gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_B(j, 1, k))) = 1!1
					gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_B(j, 1, k))) = 1!1
				!
					gl_Gran_scheme(gl_Cell_gran(3, gl_Cell_B(j, 2, k))) = 1!1
					if (gl_Cell_gran(4, gl_Cell_B(j, 2, k)) /= 0) then
						gl_Gran_scheme(gl_Cell_gran(4, gl_Cell_B(j, 2, k))) = 1!1
					end if
					gl_Gran_scheme(gl_Cell_gran(5, gl_Cell_B(j, 2, k))) = 1!1
					gl_Gran_scheme(gl_Cell_gran(6, gl_Cell_B(j, 2, k))) = 1!1
				
				end do
			end do
		end if
		
		!do k = 1, size( gl_Cell_A(par_n_HP - 1, 1, :) )
		!       do j = 1, size( gl_Cell_A(par_n_HP - 1, :, 1) )
		!		do i = 1, par_n_TS - 2
		!			gl_Gran_scheme(gl_Cell_gran(1, gl_Cell_A(i, j, k))) = 3
		!		end do
		!	end do
		!end do
		
    end subroutine Find_Surface
    
    subroutine Calc_move(now) ! Вычисляем скорости движения узлов (ключевых) с помощью распада на поверхностях разрыва(
		use STORAGE
		use GEO_PARAM
		use Solvers
		use My_func
		use cgod
		implicit none
		
		integer, intent(in) :: now   ! Какиой номер параметров используется в вычислении движения
		
		integer :: i, Num, gr, s1, s2, j, yzel, ndisc, metod, ss1, ss2
		real(8) :: normal(3), qqq1(8), qqq2(8), dsl, dsp, dsc, POTOK(8), a1, a2, v1, v2, &
			ray(3), norm, b1, b2, c1, c2, vec(3), center(3), center2(3), the1, the2, www, rast(3), df1, df2, dff2, dff1
		real(8) :: rad1, rad5, qqq1_TVD(8), aa, bb, cc, k1
		integer(4) :: kdir, KOBL, idgod
		
		KOBL = 0
		kdir = 0
		idgod = 0
		
		!koef1 = 0.3! 0.2
		!koef2 = 1.0 ! 1.0
		!koef3 = 0.7   ! 0.3
		
		ndisc = 0
		
		! Пробегаемся по всем граням и вычисляем скорости их движения
		
		! TS
		Num = size(gl_TS)
		
		do i = 1, Num
			
			metod = 3 
		www = 0.0
		KOBL = 0
		kdir = 0
		idgod = 0
		
			
		gr = gl_TS(i)
		s1 = gl_Gran_neighbour(1, gr)
		s2 = gl_Gran_neighbour(2, gr)
		ss1 = gl_Gran_neighbour_TVD(1, gr)
		ss2 = gl_Gran_neighbour_TVD(2, gr)
		normal = gl_Gran_normal2(:, gr, now)
		qqq1 = gl_Cell_par(1:8, s1)
		qqq2 = gl_Cell_par(1:8, s2)
		
		if(gl_Gran_center2(1, gr, now) > -5.0) then
		
			rast = gl_Gran_center2(:, gr, now) - gl_Cell_center2(:, s1, now)
			df1 = norm2(rast)
			rast = gl_Gran_center2(:, gr, now) - gl_Cell_center2(:, s2, now)
			df2 = norm2(rast)
			rast = gl_Gran_center2(:, gr, now) - gl_Cell_center2(:, ss2, now)
			dff2 = norm2(rast)
		
		
			rad1 = norm2(gl_Cell_center2(:, s1, now))                                
			rad5 = norm2(gl_Gran_center2(:, gr, now))
		
			qqq1_TVD = qqq1
						
			qqq1_TVD(1) = qqq1_TVD(1) * rad1**2 / rad5**2
			qqq1_TVD(5) = qqq1_TVD(5) * rad1**(2 * ggg) / rad5**(2 * ggg)
						
			call spherical_skorost(gl_Cell_center2(1, s1, now), gl_Cell_center2(2, s1, now), gl_Cell_center2(3, s1, now), &
				qqq1_TVD(2), qqq1_TVD(3), qqq1_TVD(4), aa, bb, cc)
						
			call dekard_skorost(gl_Gran_center2(1, gr, now), gl_Gran_center2(2, gr, now), gl_Gran_center2(3, gr, now), &
				aa, bb, cc, qqq1_TVD(2), qqq1_TVD(3), qqq1_TVD(4))
					
		
			qqq1 = qqq1_TVD
		
		end if
			
		call chlld(metod, normal(1), normal(2), normal(3), &
				www, qqq1, qqq2, dsl, dsp, dsc, POTOK)
			
			
			
			dsl = dsl * koef1
			
			center = gl_Gran_center2(:, gr, now)
			
			the1 = polar_angle(center(1), sqrt(center(2)**2 + center(3)**2))
			
			do j = 1, 4
					yzel = gl_all_Gran(j, gr)
					
					center2(1) = gl_x2(yzel, now)
					center2(2) = gl_y2(yzel, now)
					center2(3) = gl_z2(yzel, now)
					the2 = polar_angle(center2(1), sqrt(center2(2)**2 + center2(3)**2))
				
					if (the2 > 0.1 .and. the2 < the1 .and. the1 < par_pi_8/3.0) CYCLE
			
					gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsl
					gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsl
					gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsl
					gl_Point_num(yzel) = gl_Point_num(yzel) + 1
			end do
			CYCLE  ! Заканчиваем с этой гранью, переходим к следующей
				
			
			if(.False.) then
				center = gl_Gran_center2(:, gr, now)
			
				the1 = polar_angle(center(1), sqrt(center(2)**2 + center(3)**2))
			
				do j = 1, 4
					yzel = gl_all_Gran(j, gr)
				
					center2(1) = gl_x2(yzel, now)
					center2(2) = gl_y2(yzel, now)
					center2(3) = gl_z2(yzel, now)
					the2 = polar_angle(center2(1), sqrt(center2(2)**2 + center2(3)**2))
				
					if (the2 < the1 .and. the2 > 0.18) CYCLE
				
				
					gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsl
					gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsl
					gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsl
					gl_Point_num(yzel) = gl_Point_num(yzel) + 1
				end do
				CYCLE  ! Заканчиваем с этой гранью, переходим к следующей
			end if
			
			!if(i == 50) write(*,*) dsl
			
			a1 = sqrt(ggg * qqq1(5)/qqq1(1))  ! Скорости звука
			a2 = sqrt(ggg * qqq2(5)/qqq2(1))
			
			
			if (norm2(qqq1(2:4)) <= a1 .and. norm2(qqq2(2:4)) <= a2) then ! Записываем скорость во все узлы
				do j = 1, 4
					yzel = gl_all_Gran(j, gr)
					gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsl
					gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsl
					gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsl
					gl_Point_num(yzel) = gl_Point_num(yzel) + 1
					
				end do
				CYCLE  ! Заканчиваем с этой гранью, переходим к следующей
			end if
			
			do j = 1, 4
					yzel = gl_all_Gran(j, gr)
					
					ray = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/) - gl_Gran_center2(:, gr, now)
					norm = norm2(ray)
					ray = ray/norm
					v1 = DOT_PRODUCT(qqq1(2:4), ray)
					v2 = DOT_PRODUCT(qqq2(2:4), ray)
					
					
					if (v1 >= 0 .or. v2 >= 0) then
						
						gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsl
						gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsl
						gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsl
						gl_Point_num(yzel) = gl_Point_num(yzel) + 1
						
						CYCLE   ! Переходим к следующему узлу
					end if
					
					b1 = DOT_PRODUCT(qqq1(6:8), ray)
					b2 = DOT_PRODUCT(qqq2(6:8), ray)
					
					c1 = dabs(b1)/dsqrt(4.0 * par_pi_8 * qqq1(1))
					c2 = dabs(b2)/dsqrt(4.0 * par_pi_8 * qqq2(1))
					
					b1 = norm2(qqq1(6:8))
					b2 = norm2(qqq2(6:8))
					
					! Если модуль скорости меньше чем максимум из скорости звука и альфвеновской скорости звука по этому направлению
					if( dabs(v1) <= max(a1, 0.5 * (sqrt(b1 * b1/(4.0 * par_pi_8 * qqq1(1)) + a1 * a1+ 2.0 * a1 * c1) &
						+ sqrt(b1 * b1/(4.0 * par_pi_8 * qqq1(1)) + a1 * a1 - 2.0 * a1 * c1)  ) ) &
						.and. dabs(v2) <= max(a2, 0.5 * (sqrt(b2 * b2/(4.0 * par_pi_8 * qqq2(1)) + a2 * a2+ 2.0 * a2 * c2) &
						+ sqrt(b2 * b2/(4.0 * par_pi_8 * qqq2(1)) + a2 * a2 - 2.0 * a2 * c2)  ) ) ) then
						gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsl
						gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsl
						gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsl
						gl_Point_num(yzel) = gl_Point_num(yzel) + 1
						
						CYCLE
					end if
					
			end do
			
		end do
			
		! Контакт
		Num = size(gl_Contact)
		
		do i = 1, Num
			
			KOBL = 0
		kdir = 0
		idgod = 0
			metod = 3
		www = 0.0
		

		gr = gl_Contact(i)
		s1 = gl_Gran_neighbour(1, gr)
		s2 = gl_Gran_neighbour(2, gr)
		normal = gl_Gran_normal2(:, gr, now)
		qqq1 = gl_Cell_par(1:8, s1)
		qqq2 = gl_Cell_par(1:8, s2)
		
		k1 = 1.0
		!if (gl_Gran_center2(1, gr, now) >= 2.0) k1 = 0.7 !0.14
		
		! Вычтем нормальную компоненту магнитного поля для значений в самой ячейке
		if (gl_Gran_center2(1, gr, now) >= par_null_bn_x .and. par_null_bn == .True.) then
			qqq1(6:8) = qqq1(6:8) - DOT_PRODUCT(normal, qqq1(6:8)) * normal
			qqq2(6:8) = qqq2(6:8) - DOT_PRODUCT(normal, qqq2(6:8)) * normal
			gl_Cell_par(6:8, s1) = qqq1(6:8)
			gl_Cell_par(6:8, s2) = qqq2(6:8)
			
			
		end if
			
			vec = 0.5 * (qqq1(2:4) + qqq2(2:4))
		vec = vec - DOT_PRODUCT(vec, normal) * normal
		center = gl_Gran_center2(:, gr, now)
			
		dsc = 0.0_8
		
		if (par_null_bn == .True.) then
			! Теперь сделаем для газовой динамики
			qqq1(5) = qqq1(5) + (norm2(qqq1(6:8))**2)/(8.0 * par_pi_8)
			qqq2(5) = qqq2(5) + (norm2(qqq2(6:8))**2)/(8.0 * par_pi_8)
			call cgod3d(KOBL, 0, 0, 0, kdir, idgod, &
									normal(1), normal(2), normal(3), 1.0_8, &
									www, qqq1(1:8), qqq2(1:8), &
									dsl, dsp, dsc, 1.0_8, 1.66666666666666_8, &
									POTOK)
			if (idgod == 2) then
			call chlld_gd(3, normal(1), normal(2), normal(3), &
					www, qqq1(1:5), qqq2(1:5), dsl, dsp, dsc, POTOK)
			end if
			
		else
			call chlld(metod, normal(1), normal(2), normal(3), &
					www, qqq1, qqq2, dsl, dsp, dsc, POTOK)
		end if
			
		!if(gr == 31449) then
		!	print*, "DSC = ", dsc
		!	PAUSE
		!end if
		
		!print*, "1 = ", dsc
		dsc = k1 * (dsc) * koef2
		
		
		!print*, "2 = ", dsc
		!PAUSE
			
			!if (gl_Gran_center2(1, gr, now) < -200.0 .and. gl_Gran_center2(1, gr, now) > -220.0 .and. normal(2) > 0 .and. &
			!	dabs(gl_Gran_center2(3, gr, now)) < 5.0) then
			!	print*, gl_Gran_center2(:, gr, now) 
			!	dsc = dsc + 30.0
			!end if
			
			!if(gr == 77) write(*, *) dsc, qqq1(1), qqq2(1), qqq1(5), qqq2(5)
			
			do j = 1, 4
				yzel = gl_all_Gran(j, gr)
				gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsc
				gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsc
				gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsc
				gl_Point_num(yzel) = gl_Point_num(yzel) + 1
			end do
			CYCLE  ! Заканчиваем с этой гранью, переходим к следующей
			
			
			!if(i == 2) write(*,*) dsc
			a1 = sqrt(ggg * qqq1(5)/qqq1(1))  ! Скорости звука
			a2 = sqrt(ggg * qqq2(5)/qqq2(1))
			
			
			if (norm2(qqq1(2:4)) <= a1 .and. norm2(qqq2(2:4)) <= a2) then ! Записываем скорость во все узлы
				do j = 1, 4
					yzel = gl_all_Gran(j, gr)
					gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsc
					gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsc
					gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsc
					gl_Point_num(yzel) = gl_Point_num(yzel) + 1
				end do
				CYCLE  ! Заканчиваем с этой гранью, переходим к следующей
			end if
			
			do j = 1, 4
					yzel = gl_all_Gran(j, gr)
					
					ray = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/) - gl_Gran_center2(:, gr, now)
					norm = norm2(ray)
					ray = ray/norm
					v1 = DOT_PRODUCT(qqq1(2:4), ray)
					v2 = DOT_PRODUCT(qqq2(2:4), ray)
					
					if (v1 >= 0 .or. v2 >= 0) then
						gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsc
						gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsc
						gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsc
						gl_Point_num(yzel) = gl_Point_num(yzel) + 1
						CYCLE   ! Переходим к следующему узлу
					end if
					
					b1 = DOT_PRODUCT(qqq1(6:8), ray)
					b2 = DOT_PRODUCT(qqq2(6:8), ray)
					
					c1 = dabs(b1)/dsqrt(4.0 * par_pi_8 * qqq1(1))
					c2 = dabs(b2)/dsqrt(4.0 * par_pi_8 * qqq2(1))
					
					b1 = norm2(qqq1(6:8))
					b2 = norm2(qqq2(6:8))
					
					! Если модуль скорости меньше чем максимум из скорости звука и альфвеновской скорости звука по этому направлению
					if( dabs(v1) <= max(a1, 0.5 * (sqrt(b1 * b1/(4.0 * par_pi_8 * qqq1(1)) + a1 * a1+ 2.0 * a1 * c1) &
						+ sqrt(b1 * b1/(4.0 * par_pi_8 * qqq1(1)) + a1 * a1 - 2.0 * a1 * c1)  ) ) &
						.and. dabs(v2) <= max(a2, 0.5 * (sqrt(b2 * b2/(4.0 * par_pi_8 * qqq2(1)) + a2 * a2+ 2.0 * a2 * c2) &
						+ sqrt(b2 * b2/(4.0 * par_pi_8 * qqq2(1)) + a2 * a2 - 2.0 * a2 * c2)  ) ) ) then
						gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsc
						gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsc
						gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsc
						gl_Point_num(yzel) = gl_Point_num(yzel) + 1
						CYCLE
					end if
					
			end do
			
		end do
		
		
		! BS
		if(.True.) then
		Num = size(gl_BS)
		
			do i = 1, Num
			
				gr = gl_BS(i)
				s1 = gl_Gran_neighbour(1, gr)
				s2 = gl_Gran_neighbour(2, gr)
				normal = gl_Gran_normal2(:, gr, now)
				qqq1 = gl_Cell_par(1:8, s1)
				qqq2 = gl_Cell_par(1:8, s2)
				KOBL = 0
				kdir = 0
				idgod = 0
			
				call chlld(2, normal(1), normal(2), normal(3), &  ! 2
						0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)
			
				dsp = dsp * koef3
				!if(i == 2) then
			!       write(*,*) dsp
			!       write(*,*) normal(1), normal(2), normal(3)
			!       write(*,*) 2
			!       write(*,*) qqq1(1), qqq1(2), qqq1(3), qqq1(4), qqq1(5), qqq1(6), qqq1(7), qqq1(8)
			!       write(*,*) qqq2(1), qqq2(2), qqq2(3), qqq2(4), qqq2(5), qqq2(6), qqq2(7), qqq2(8)
		!      end if
			
				a1 = sqrt(ggg * qqq1(5)/qqq1(1))  ! Скорости звука
				a2 = sqrt(ggg * qqq2(5)/qqq2(1))
			
				if (norm2(qqq1(2:4)) <= a1 .and. norm2(qqq2(2:4)) <= a2) then ! Записываем скорость во все узлы
					do j = 1, 4
						yzel = gl_all_Gran(j, gr)
						gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsp
						gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsp
						gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsp
						gl_Point_num(yzel) = gl_Point_num(yzel) + 1
					end do
					CYCLE  ! Заканчиваем с этой гранью, переходим к следующей
				end if
			
				do j = 1, 4
						yzel = gl_all_Gran(j, gr)
					
						ray = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/) - gl_Gran_center2(:, gr, now)
						norm = norm2(ray)
						ray = ray/norm
						v1 = DOT_PRODUCT(qqq1(2:4), ray)
						v2 = DOT_PRODUCT(qqq2(2:4), ray)
					
						if (v1 >= 0 .or. v2 >= 0) then
							gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsp
							gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsp
							gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsp
							gl_Point_num(yzel) = gl_Point_num(yzel) + 1
							CYCLE   ! Переходим к следующему узлу
						end if
					
						b1 = DOT_PRODUCT(qqq1(6:8), ray)
						b2 = DOT_PRODUCT(qqq2(6:8), ray)
					
						c1 = dabs(b1)/dsqrt(4.0 * par_pi_8 * qqq1(1))
						c2 = dabs(b2)/dsqrt(4.0 * par_pi_8 * qqq2(1))
					
						b1 = norm2(qqq1(6:8))
						b2 = norm2(qqq2(6:8))
					
						! Если модуль скорости меньше чем максимум из скорости звука и альфвеновской скорости звука по этому направлению
						if( dabs(v1) <= max(a1, 0.5 * (sqrt(b1 * b1/(4.0 * par_pi_8 * qqq1(1)) + a1 * a1+ 2.0 * a1 * c1) &
							+ sqrt(b1 * b1/(4.0 * par_pi_8 * qqq1(1)) + a1 * a1 - 2.0 * a1 * c1)  ) ) &
							.and. dabs(v2) <= max(a2, 0.5 * (sqrt(b2 * b2/(4.0 * par_pi_8 * qqq2(1)) + a2 * a2+ 2.0 * a2 * c2) &
							+ sqrt(b2 * b2/(4.0 * par_pi_8 * qqq2(1)) + a2 * a2 - 2.0 * a2 * c2)  ) ) ) then
							gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsp
							gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsp
							gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsp
							gl_Point_num(yzel) = gl_Point_num(yzel) + 1
							CYCLE
						end if
					
				end do
			
			end do
		end if
	
    end subroutine Calc_move   
    
	subroutine Move_all(now, Time)  ! Передвигаем точки в соответствии со скоростью движение узлов на поверхностях разрыва
		use STORAGE
		use GEO_PARAM
		use My_func
		use MY_CUDA_smooth
		implicit none
		
		integer :: now2             ! Эти параметры мы сейчас меняем на основе now
		real(8), intent(in) :: Time
		integer, intent(in) :: now
		real(8) :: R_TS, proect, vel(3), R_HP, R_BS, Ak(3), Bk(3), Ck(3), Dk(3), Ek(3), ER(3), ER2(3), KORD(3), dist, ddt
		integer :: yzel, N1, N2, N3, i, j, k, yzel2, yzel3, yzel22, yzel33, yzel333, yzel222, yzel2222, yzel3333
		real(8) :: the, the2, phi, r, x, y, z, rr, xx, x2, y2, z2, rrr, r1, r2, r3, r4, rd, kk13, kkk, k1
		
		now2 = mod(now, 2) + 1   
		
		ddt = Time/0.000127
		
		
		! Поверхностное натяжение ------------------------------------------------------------------------------
		if (.True.) then
		N3 = size(gl_RAY_A(1, 1, :))
		N2 = size(gl_RAY_A(1, :, 1))
		N1 = size(gl_RAY_A(:, 1, 1))

		! Цикл движения точек на лучах А  
		do k = 1, N3
			do j = 1, N2
				
				k1 = 0.002! 0.001
		if(j <= 18) k1 = 0.03   ! 0.03
		
		
		if (j == 1 .and. k/= 1) then
						CYCLE
		end if
		
		ddt = Time/0.000127
				
		! Ударная волна
		yzel = gl_RAY_A(par_n_TS, j, k)
			Ak(1) = gl_x2(yzel, now); Ak(2) = gl_y2(yzel, now); Ak(3) = gl_z2(yzel, now)
				
		if (j < N2) then
			yzel2 = gl_RAY_A(par_n_TS, j + 1, k)
		else
			yzel2 = gl_RAY_B(par_n_TS, 1, k)
		end if
		!Bk(1) = gl_x2(yzel2, now); Bk(2) = gl_y2(yzel2, now); Bk(3) = gl_z2(yzel2, now)
		
		if (j < N2 - 1) then
			yzel22 = gl_RAY_A(par_n_TS, j + 2, k)
		else if (j == N2 - 1) then
			yzel22 = gl_RAY_B(par_n_TS, 1, k)
		else
			yzel22 = gl_RAY_B(par_n_TS, 2, k)
		end if
		
		if (j < N2 - 2) then
			yzel222 = gl_RAY_A(par_n_TS, j + 3, k)
		else if (j == N2 - 2) then
			yzel222 = gl_RAY_B(par_n_TS, 1, k)
		else if (j == N2 - 1) then
			yzel222 = gl_RAY_B(par_n_TS, 2, k)
		else
			yzel222 = gl_RAY_B(par_n_TS, 3, k)
		end if
				
		if (j > 1) then
			yzel3 = gl_RAY_A(par_n_TS, j - 1, k)
		else
			yzel3 = gl_RAY_A(par_n_TS, j + 1, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
		end if
		!Dk(1) = gl_x2(yzel3, now); Dk(2) = gl_y2(yzel3, now); Dk(3) = gl_z2(yzel3, now)
		
		if (j > 2) then
			yzel33 = gl_RAY_A(par_n_TS, j - 2, k)
		else if(j == 2) then
			yzel33 = gl_RAY_A(par_n_TS, 2, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
		else
			yzel33 = gl_RAY_A(par_n_TS, 3, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
		end if
		
		if (j > 3) then
			yzel333 = gl_RAY_A(par_n_TS, j - 3, k)
		else if(j == 3) then
			yzel333 = gl_RAY_A(par_n_TS, 2, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
		else if(j == 2) then
			yzel333 = gl_RAY_A(par_n_TS, 3, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
		else
			yzel333 = gl_RAY_A(par_n_TS, 4, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
		end if
		
		call Smooth_kvadr4(gl_x2(yzel333, now), gl_y2(yzel333, now), gl_z2(yzel333, now), &
				gl_x2(yzel33, now), gl_y2(yzel33, now), gl_z2(yzel33, now), &
				gl_x2(yzel3, now), gl_y2(yzel3, now), gl_z2(yzel3, now), &
				Ak(1), Ak(2), Ak(3), &
				gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now), &
				gl_x2(yzel22, now), gl_y2(yzel22, now), gl_z2(yzel22, now), &
				gl_x2(yzel222, now), gl_y2(yzel222, now), gl_z2(yzel222, now), &
				Ek(1), Ek(2), Ek(3))
		
		
		if (k > 1) then
			yzel2 = gl_RAY_A(par_n_TS, j, k - 1)
		else
			yzel2 = gl_RAY_A(par_n_TS, j, N3)
		end if
		!Ck(1) = gl_x2(yzel2, now); Ck(2) = gl_y2(yzel2, now); Ck(3) = gl_z2(yzel2, now)
		
		if (k > 2) then
			yzel33 = gl_RAY_A(par_n_TS, j, k - 2)
		else if(k == 2) then
			yzel33 = gl_RAY_A(par_n_TS, j, N3)
		else
			yzel33 = gl_RAY_A(par_n_TS, j, N3 - 1)
		end if
		
		if (k > 3) then
			yzel333 = gl_RAY_A(par_n_TS, j, k - 3)
		else if(k == 3) then
			yzel333 = gl_RAY_A(par_n_TS, j, N3)
		else if(k == 2) then
			yzel333 = gl_RAY_A(par_n_TS, j, N3 - 1)
		else
			yzel333 = gl_RAY_A(par_n_TS, j, N3 - 2)
		end if
		
		if(k < N3) then
			yzel2 = gl_RAY_A(par_n_TS, j, k + 1)
		else
			yzel2 = gl_RAY_A(par_n_TS, j, 1)
		end if
		!Ek(1) = gl_x2(yzel2, now); Ek(2) = gl_y2(yzel2, now); Ek(3) = gl_z2(yzel2, now)
		
		if(k < N3 - 1) then
			yzel22 = gl_RAY_A(par_n_TS, j, k + 2)
		else if(k == N3 - 1) then
			yzel22 = gl_RAY_A(par_n_TS, j, 1)
		else
			yzel22 = gl_RAY_A(par_n_TS, j, 2)
		end if
		
		if(k < N3 - 2) then
				yzel222 = gl_RAY_A(par_n_TS, j, k + 2)
			else if(k == N3 - 2) then
				yzel222 = gl_RAY_A(par_n_TS, j, 1)
			else if(k == N3 - 1) then
				yzel222 = gl_RAY_A(par_n_TS, j, 2)
			else
				yzel222 = gl_RAY_A(par_n_TS, j, 3)
			end if
		
		
		call Smooth_kvadr4(gl_x2(yzel333, now), gl_y2(yzel333, now), gl_z2(yzel333, now), &
				gl_x2(yzel33, now), gl_y2(yzel33, now), gl_z2(yzel33, now), &
				gl_x2(yzel3, now), gl_y2(yzel3, now), gl_z2(yzel3, now), &
				Ak(1), Ak(2), Ak(3), &
				gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now), &
				gl_x2(yzel22, now), gl_y2(yzel22, now), gl_z2(yzel22, now), &
				gl_x2(yzel222, now), gl_y2(yzel222, now), gl_z2(yzel222, now), &
				Ck(1), Ck(2), Ck(3))
				
		
		k1 = 30.0 !10
		if(j <= 5) k1 = 150.0   ! 0.03
		
		if (gl_Point_num(yzel) > 0) then
			vel = k1 * par_nat_TS * 0.006 * (Ck/2.0 + Ek/2.0 - Ak) * gl_Point_num(yzel)/Time * ddt  ! 0.003
		else
			vel = k1 * par_nat_TS * 0.006 * (Ck/2.0 + Ek/2.0 - Ak)/Time * ddt   ! 0.003
		end if
			
		if(j >= N2 - 1) then
			if (gl_Point_num(yzel) > 0) then
				vel = k1 * par_nat_TS * 0.006 * ( (Ck + 6.0 * Ek)/7.0 - Ak) * gl_Point_num(yzel)/Time * ddt  ! 0.003
			else
				vel = k1 *  par_nat_TS * 0.006 * ( (Ck + 6.0 * Ek)/7.0 - Ak)/Time * ddt   ! 0.003
			end if
		end if
		
				
		gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
		gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
		gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
		
		
		! Контакт
			yzel = gl_RAY_A(par_n_HP, j, k)
			Ak(1) = gl_x2(yzel, now); Ak(2) = gl_y2(yzel, now); Ak(3) = gl_z2(yzel, now)
			
			if (j < N2) then
				yzel2 = gl_RAY_A(par_n_HP, j + 1, k)
			else
				yzel2 = gl_RAY_B(par_n_HP, 1, k)
			end if
		
			if (j < N2 - 1) then
				yzel22 = gl_RAY_A(par_n_HP, j + 2, k)
			else if (j == N2 - 1) then
				yzel22 = gl_RAY_B(par_n_HP, 1, k)
			else
				yzel22 = gl_RAY_B(par_n_HP, 2, k)
			end if
		
			if (j < N2 - 2) then
				yzel222 = gl_RAY_A(par_n_HP, j + 3, k)
			else if (j == N2 - 2) then
				yzel222 = gl_RAY_B(par_n_HP, 1, k)
			else if (j == N2 - 1) then
				yzel222 = gl_RAY_B(par_n_HP, 2, k)
			else
				yzel222 = gl_RAY_B(par_n_HP, 3, k)
			end if
			
			if (j > 1) then
				yzel3 = gl_RAY_A(par_n_HP, j - 1, k)
			else
				yzel3 = gl_RAY_A(par_n_HP, 1, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
			end if
			
			if (j > 2) then
				yzel33 = gl_RAY_A(par_n_HP, j - 2, k)
			else if(j == 2) then
				yzel33 = gl_RAY_A(par_n_HP, 1, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
			else
				yzel33 = gl_RAY_A(par_n_HP, 2, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
			end if
		
			if (j > 3) then
				yzel333 = gl_RAY_A(par_n_HP, j - 3, k)
			else if(j == 3) then
				yzel333 = gl_RAY_A(par_n_HP, 1, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
			else if(j == 2) then
				yzel333 = gl_RAY_A(par_n_HP, 2, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
			else
				yzel333 = gl_RAY_A(par_n_HP, 3, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
			end if
		
			call Smooth_kvadr4(gl_x2(yzel333, now), gl_y2(yzel333, now), gl_z2(yzel333, now), &
				gl_x2(yzel33, now), gl_y2(yzel33, now), gl_z2(yzel33, now), &
				gl_x2(yzel3, now), gl_y2(yzel3, now), gl_z2(yzel3, now), &
				Ak(1), Ak(2), Ak(3), &
				gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now), &
				gl_x2(yzel22, now), gl_y2(yzel22, now), gl_z2(yzel22, now), &
				gl_x2(yzel222, now), gl_y2(yzel222, now), gl_z2(yzel222, now), &
				Ek(1), Ek(2), Ek(3))
		
			
			if(k < N3) then
				yzel2 = gl_RAY_A(par_n_HP, j, k + 1)
			else
				yzel2 = gl_RAY_A(par_n_HP, j, 1)
			end if
			
			if(k < N3 - 1) then
				yzel22 = gl_RAY_A(par_n_HP, j, k + 2)
			else if(k == N3 - 1) then
				yzel22 = gl_RAY_A(par_n_HP, j, 1)
			else
				yzel22 = gl_RAY_A(par_n_HP, j, 2)
			end if
		
			if(k < N3 - 2) then
				yzel222 = gl_RAY_A(par_n_HP, j, k + 2)
			else if(k == N3 - 2) then
				yzel222 = gl_RAY_A(par_n_HP, j, 1)
			else if(k == N3 - 1) then
				yzel222 = gl_RAY_A(par_n_HP, j, 2)
			else
				yzel222 = gl_RAY_A(par_n_HP, j, 3)
			end if
		
			if (k > 1) then
				yzel3 = gl_RAY_A(par_n_HP, j, k - 1)
			else
				yzel3 = gl_RAY_A(par_n_HP, j, N3)
			end if
		
			if (k > 2) then
				yzel33 = gl_RAY_A(par_n_HP, j, k - 2)
			else if(k == 2) then
				yzel33 = gl_RAY_A(par_n_HP, j, N3)
			else
				yzel33 = gl_RAY_A(par_n_HP, j, N3 - 1)
			end if
		
			if (k > 3) then
				yzel333 = gl_RAY_A(par_n_HP, j, k - 3)
			else if(k == 3) then
				yzel333 = gl_RAY_A(par_n_HP, j, N3)
			else if(k == 2) then
				yzel333 = gl_RAY_A(par_n_HP, j, N3 - 1)
			else
				yzel333 = gl_RAY_A(par_n_HP, j, N3 - 2)
			end if
		
			call Smooth_kvadr4(gl_x2(yzel333, now), gl_y2(yzel333, now), gl_z2(yzel333, now), &
				gl_x2(yzel33, now), gl_y2(yzel33, now), gl_z2(yzel33, now), &
				gl_x2(yzel3, now), gl_y2(yzel3, now), gl_z2(yzel3, now), &
				Ak(1), Ak(2), Ak(3), &
				gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now), &
				gl_x2(yzel22, now), gl_y2(yzel22, now), gl_z2(yzel22, now), &
				gl_x2(yzel222, now), gl_y2(yzel222, now), gl_z2(yzel222, now), &
				Ck(1), Ck(2), Ck(3))
		
			
			if (gl_Point_num(yzel) > 0) then
				!vel = par_nat_HP * 0.006 * (Bk/4.0 + Ck/4.0 + Dk/4.0 + Ek/4.0 - Ak) * gl_Point_num(yzel)/Time  ! 0.00006
				vel = par_nat_HP * 0.003 * (Ck/2.0 + Ek/2.0 - Ak) * gl_Point_num(yzel)/Time  ! 0.00006
			else
				!vel = par_nat_HP * 0.006 * (Bk/4.0 + Ck/4.0 + Dk/4.0 + Ek/4.0 - Ak)/Time   !  0.00006
				vel = par_nat_HP * 0.003 * (Ck/2.0 + Ek/2.0 - Ak)/Time   ! 0.001     0.00006
			end if
			
				
		gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
		gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
		gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
				
				
		! Внешняя ударная волна
				
		yzel = gl_RAY_A(par_n_BS, j, k)
		Ak(1) = gl_x2(yzel, now); Ak(2) = gl_y2(yzel, now); Ak(3) = gl_z2(yzel, now)
				
		if (j < N2) then
			yzel2 = gl_RAY_A(par_n_BS, j + 1, k)
		else
			CYCLE !yzel2 = yzel                                                    ! Вот тут есть неопределённость что делать с крайним узлом
		end if
		Bk(1) = gl_x2(yzel2, now); Bk(2) = gl_y2(yzel2, now); Bk(3) = gl_z2(yzel2, now)
				
		if (k > 1) then
			yzel2 = gl_RAY_A(par_n_BS, j, k - 1)
		else
			yzel2 = gl_RAY_A(par_n_BS, j, N3)
		end if
		Ck(1) = gl_x2(yzel2, now); Ck(2) = gl_y2(yzel2, now); Ck(3) = gl_z2(yzel2, now)
				
		if (j > 1) then
			yzel2 = gl_RAY_A(par_n_BS, j - 1, k)
		else
			yzel2 = gl_RAY_A(par_n_BS, j, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
		end if
		Dk(1) = gl_x2(yzel2, now); Dk(2) = gl_y2(yzel2, now); Dk(3) = gl_z2(yzel2, now)
				
		if(k < N3) then
			yzel2 = gl_RAY_A(par_n_BS, j, k + 1)
		else
			yzel2 = gl_RAY_A(par_n_BS, j, 1)
		end if
		Ek(1) = gl_x2(yzel2, now); Ek(2) = gl_y2(yzel2, now); Ek(3) = gl_z2(yzel2, now)
		
		dist = (sqrt( (Dk(1) - Ak(1))**2 + (Dk(2) - Ak(2))**2 + (Dk(3) - Ak(3))**2 ) + &
					sqrt( (Ek(1) - Ak(1))**2 + (Ek(2) - Ak(2))**2 + (Ek(3) - Ak(3))**2 ))/2.0
		
		dist = max(dist, 1.0_8)
				
		if (gl_Point_num(yzel) > 0) then
			vel = gl_Point_num(yzel) * par_nat_BS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time/dist * ddt
		else
			vel = par_nat_BS * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time/dist * ddt
		end if
				
				
		gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
		gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
		gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
				
			end do
		end do
		
		
		N3 = size(gl_RAY_B(1, 1, :))
		N2 = size(gl_RAY_B(1, :, 1))
		N1 = size(gl_RAY_B(:, 1, 1))

		! Цикл движения точек на лучах B  
		do k = 1, N3
			do j = 1, N2
				
				! Ударная волна
				
				yzel = gl_RAY_B(par_n_TS, j, k)
				Ak(1) = gl_x2(yzel, now); Ak(2) = gl_y2(yzel, now); Ak(3) = gl_z2(yzel, now)
				! Ak = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/)
				
				if (j < N2) then
					yzel2 = gl_RAY_B(par_n_TS, j + 1, k)
				else
					yzel2 = gl_RAY_K(par_n_TS, size(gl_RAY_K(par_n_TS, :, k)), k)
				end if
				!Bk(1) = gl_x2(yzel2, now); Bk(2) = gl_y2(yzel2, now); Bk(3) = gl_z2(yzel2, now)
				! Bk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
				
				if (j < N2 - 1) then
					yzel22 = gl_RAY_B(par_n_TS, j + 2, k)
				else if (j == N2 - 1) then
					yzel22 = gl_RAY_K(par_n_TS, size(gl_RAY_K(par_n_TS, :, k)), k)
				else
					yzel22 = gl_RAY_K(par_n_TS, size(gl_RAY_K(par_n_TS, :, k)) - 1, k)
				end if
		
				if (j < N2 - 2) then
					yzel222 = gl_RAY_B(par_n_TS, j + 3, k)
				else if (j == N2 - 2) then
					yzel222 = gl_RAY_K(par_n_TS, size(gl_RAY_K(par_n_TS, :, k)), k)
				else if (j == N2 - 1) then
					yzel222 = gl_RAY_K(par_n_TS, size(gl_RAY_K(par_n_TS, :, k)) - 1, k)
				else
					yzel222 = gl_RAY_K(par_n_TS, size(gl_RAY_K(par_n_TS, :, k))	- 2, k)
				end if
				
				if (j < N2 - 3) then
					yzel2222 = gl_RAY_B(par_n_TS, j + 4, k)
				else if (j == N2 - 3) then
					yzel2222 = gl_RAY_K(par_n_TS, size(gl_RAY_K(par_n_TS, :, k)) - 1, k)
				else if (j == N2 - 2) then
					yzel2222 = gl_RAY_K(par_n_TS, size(gl_RAY_K(par_n_TS, :, k)) - 2, k)
				else if (j == N2 - 1) then
					yzel2222 = gl_RAY_K(par_n_TS, size(gl_RAY_K(par_n_TS, :, k)) - 3, k)
				else
					yzel2222 = gl_RAY_K(par_n_TS, size(gl_RAY_K(par_n_TS, :, k)) - 4, k)
				end if
		
				if (j > 1) then
					yzel3 = gl_RAY_B(par_n_TS, j - 1, k)
				else
					yzel3 = gl_RAY_A(par_n_TS, size( gl_RAY_A(par_n_TS, :, k)), k)
				end if
				!Dk(1) = gl_x2(yzel2, now); Dk(2) = gl_y2(yzel2, now); Dk(3) = gl_z2(yzel2, now)
				! Dk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
		
				if (j > 2) then
					yzel33 = gl_RAY_B(par_n_TS, j - 2, k)
				else if(j == 2) then
					yzel33 = gl_RAY_A(par_n_TS, size( gl_RAY_A(par_n_TS, :, k)), k)
				else
					yzel33 = gl_RAY_A(par_n_TS, size( gl_RAY_A(par_n_TS, :, k)) - 1, k)
				end if
		
				if (j > 3) then
					yzel333 = gl_RAY_B(par_n_TS, j - 3, k)
				else if(j == 3) then
					yzel333 = gl_RAY_A(par_n_TS, size( gl_RAY_A(par_n_TS, :, k)), k)
				else if(j == 2) then
					yzel333 = gl_RAY_A(par_n_TS, size( gl_RAY_A(par_n_TS, :, k)) - 1, k)
				else
					yzel333 = gl_RAY_A(par_n_TS, size( gl_RAY_A(par_n_TS, :, k)) - 2, k)
				end if
				
				if (j > 4) then
					yzel3333 = gl_RAY_B(par_n_TS, j - 4, k)
				else if(j == 4) then
					yzel3333 = gl_RAY_A(par_n_TS, size( gl_RAY_A(par_n_TS, :, k)), k)
				else if(j == 3) then
					yzel3333 = gl_RAY_A(par_n_TS, size( gl_RAY_A(par_n_TS, :, k)) - 1, k)
				else if(j == 2) then
					yzel3333 = gl_RAY_A(par_n_TS, size( gl_RAY_A(par_n_TS, :, k)) - 2, k)
				else
					yzel3333 = gl_RAY_A(par_n_TS, size( gl_RAY_A(par_n_TS, :, k)) - 3, k)
				end if
				
				call Smooth_kvadr5(gl_x2(yzel3333, now), gl_y2(yzel3333, now), gl_z2(yzel3333, now), &
					gl_x2(yzel333, now), gl_y2(yzel333, now), gl_z2(yzel333, now), &
				gl_x2(yzel33, now), gl_y2(yzel33, now), gl_z2(yzel33, now), &
				gl_x2(yzel3, now), gl_y2(yzel3, now), gl_z2(yzel3, now), &
				Ak(1), Ak(2), Ak(3), &
				gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now), &
				gl_x2(yzel22, now), gl_y2(yzel22, now), gl_z2(yzel22, now), &
				gl_x2(yzel222, now), gl_y2(yzel222, now), gl_z2(yzel222, now), &
					gl_x2(yzel2222, now), gl_y2(yzel2222, now), gl_z2(yzel2222, now), &
				Ek(1), Ek(2), Ek(3))
				
				if (k > 1) then
					yzel3 = gl_RAY_B(par_n_TS, j, k - 1)
				else
					yzel3 = gl_RAY_B(par_n_TS, j, N3)
				end if
				!Ck(1) = gl_x2(yzel2, now); Ck(2) = gl_y2(yzel2, now); Ck(3) = gl_z2(yzel2, now)
				! Ck = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
				
				if (k > 2) then
					yzel33 = gl_RAY_B(par_n_TS, j, k - 2)
				else if(k == 2) then
					yzel33 = gl_RAY_B(par_n_TS, j, N3)
				else
					yzel33 = gl_RAY_B(par_n_TS, j, N3 - 1)
				end if
		
				if (k > 3) then
					yzel333 = gl_RAY_B(par_n_TS, j, k - 3)
				else if(k == 3) then
					yzel333 = gl_RAY_B(par_n_TS, j, N3)
				else if(k == 2) then
					yzel333 = gl_RAY_B(par_n_TS, j, N3 - 1)
				else
					yzel333 = gl_RAY_B(par_n_TS, j, N3 - 2)
				end if
				
				if(k < N3) then
					yzel2 = gl_RAY_B(par_n_TS, j, k + 1)
				else
					yzel2 = gl_RAY_B(par_n_TS, j, 1)
				end if
				!Ek(1) = gl_x2(yzel2, now); Ek(2) = gl_y2(yzel2, now); Ek(3) = gl_z2(yzel2, now)
				! Ek = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
				
				if(k < N3 - 1) then
					yzel22 = gl_RAY_B(par_n_TS, j, k + 2)
				else if(k == N3 - 1) then
					yzel22 = gl_RAY_B(par_n_TS, j, 1)
				else
					yzel22 = gl_RAY_B(par_n_TS, j, 2)
				end if
		
				if(k < N3 - 2) then
					yzel222 = gl_RAY_B(par_n_TS, j, k + 2)
				else if(k == N3 - 2) then
					yzel222 = gl_RAY_B(par_n_TS, j, 1)
				else if(k == N3 - 1) then
					yzel222 = gl_RAY_B(par_n_TS, j, 2)
				else
					yzel222 = gl_RAY_B(par_n_TS, j, 3)
				end if
		
				call Smooth_kvadr4(gl_x2(yzel333, now), gl_y2(yzel333, now), gl_z2(yzel333, now), &
				gl_x2(yzel33, now), gl_y2(yzel33, now), gl_z2(yzel33, now), &
				gl_x2(yzel3, now), gl_y2(yzel3, now), gl_z2(yzel3, now), &
				Ak(1), Ak(2), Ak(3), &
				gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now), &
				gl_x2(yzel22, now), gl_y2(yzel22, now), gl_z2(yzel22, now), &
				gl_x2(yzel222, now), gl_y2(yzel222, now), gl_z2(yzel222, now), &
				Ck(1), Ck(2), Ck(3))
				
				if (gl_Point_num(yzel) > 0) then
					vel = 10.0 * gl_Point_num(yzel) * par_nat_TS * (Ck/2.0 + Ek/2.0 - Ak)/Time * ddt
				else
					vel = 10.0 * par_nat_TS * (Ck/2.0 + Ek/2.0 - Ak)/Time * ddt
				end if
				
				
				gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
				gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
				gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
				
				!return
				kkk = 6.0 !0.2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				!if(j < 15) kkk = 6.0
				! Контакт
				if (.True.) then !(j < N2 - 3) then
				yzel = gl_RAY_B(par_n_HP, j, k)
				Ak(1) = gl_x2(yzel, now); Ak(2) = gl_y2(yzel, now); Ak(3) = gl_z2(yzel, now)
				
				!if(Ak(2) < 0 .and. Ak(3) > 0) kkk = 6.0
			
				if (j < N2) then
					yzel2 = gl_RAY_B(par_n_HP, j + 1, k)
				else
					yzel2 = gl_RAY_O(1, 1, k)
				end if
		
				if (j < N2 - 1) then
					yzel22 = gl_RAY_B(par_n_HP, j + 2, k)
				else if (j == N2 - 1) then
					yzel22 = gl_RAY_O(1, 1, k)
				else
					yzel22 = gl_RAY_O(1, 2, k)
				end if
		
				if (j < N2 - 2) then
					yzel222 = gl_RAY_B(par_n_HP, j + 3, k)
				else if (j == N2 - 2) then
					yzel222 = gl_RAY_O(1, 1, k)
				else if (j == N2 - 1) then
					yzel222 = gl_RAY_O(1, 2, k)
				else
					yzel222 = gl_RAY_O(1, 3, k)
				end if
				
				if (j < N2 - 3) then
					yzel2222 = gl_RAY_B(par_n_HP, j + 4, k)
				else if (j == N2 - 3) then
					yzel2222 = gl_RAY_O(1, 1, k)
				else if (j == N2 - 2) then
					yzel2222 = gl_RAY_O(1, 2, k)
				else if (j == N2 - 1) then
					yzel2222 = gl_RAY_O(1, 3, k)
				else
					yzel2222 = gl_RAY_O(1, 4, k)
				end if
			
				if (j > 1) then
					yzel3 = gl_RAY_B(par_n_HP, j - 1, k)
				else
					yzel3 = gl_RAY_A(par_n_HP, size( gl_RAY_A(par_n_HP, :, k)), k)
				end if
			
				if (j > 2) then
					yzel33 = gl_RAY_B(par_n_HP, j - 2, k)
				else if(j == 2) then
					yzel33 = gl_RAY_A(par_n_HP, size( gl_RAY_A(par_n_HP, :, k)), k)
				else
					yzel33 = gl_RAY_A(par_n_HP, size( gl_RAY_A(par_n_HP, :, k)) - 1, k)
				end if
		
				if (j > 3) then
					yzel333 = gl_RAY_B(par_n_HP, j - 3, k)
				else if(j == 3) then
					yzel333 = gl_RAY_A(par_n_HP, size( gl_RAY_A(par_n_HP, :, k)), k)
				else if(j == 2) then
					yzel333 = gl_RAY_A(par_n_HP, size( gl_RAY_A(par_n_HP, :, k)) - 1, k)
				else
					yzel333 = gl_RAY_A(par_n_HP, size( gl_RAY_A(par_n_HP, :, k)) - 2, k)
				end if
				
				if (j > 4) then
					yzel3333 = gl_RAY_B(par_n_HP, j - 4, k)
				else if(j == 4) then
					yzel3333 = gl_RAY_A(par_n_HP, size( gl_RAY_A(par_n_HP, :, k)), k)
				else if(j == 3) then
					yzel3333 = gl_RAY_A(par_n_HP, size( gl_RAY_A(par_n_HP, :, k)) - 1, k)
				else if(j == 2) then
					yzel3333 = gl_RAY_A(par_n_HP, size( gl_RAY_A(par_n_HP, :, k)) - 2, k)
				else
					yzel3333 = gl_RAY_A(par_n_HP, size( gl_RAY_A(par_n_HP, :, k)) - 3, k)
				end if
				
				call Smooth_kvadr5(gl_x2(yzel3333, now), gl_y2(yzel3333, now), gl_z2(yzel3333, now), &
					gl_x2(yzel333, now), gl_y2(yzel333, now), gl_z2(yzel333, now), &
				gl_x2(yzel33, now), gl_y2(yzel33, now), gl_z2(yzel33, now), &
				gl_x2(yzel3, now), gl_y2(yzel3, now), gl_z2(yzel3, now), &
				Ak(1), Ak(2), Ak(3), &
				gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now), &
				gl_x2(yzel22, now), gl_y2(yzel22, now), gl_z2(yzel22, now), &
				gl_x2(yzel222, now), gl_y2(yzel222, now), gl_z2(yzel222, now), &
					gl_x2(yzel2222, now), gl_y2(yzel2222, now), gl_z2(yzel2222, now), &
				Ek(1), Ek(2), Ek(3))
		
			
				if(k < N3) then
					yzel2 = gl_RAY_B(par_n_HP, j, k + 1)
				else
					yzel2 = gl_RAY_B(par_n_HP, j, 1)
				end if
			
				if(k < N3 - 1) then
					yzel22 = gl_RAY_B(par_n_HP, j, k + 2)
				else if(k == N3 - 1) then
					yzel22 = gl_RAY_B(par_n_HP, j, 1)
				else
					yzel22 = gl_RAY_B(par_n_HP, j, 2)
				end if
		
				if(k < N3 - 2) then
					yzel222 = gl_RAY_B(par_n_HP, j, k + 2)
				else if(k == N3 - 2) then
					yzel222 = gl_RAY_B(par_n_HP, j, 1)
				else if(k == N3 - 1) then
					yzel222 = gl_RAY_B(par_n_HP, j, 2)
				else
					yzel222 = gl_RAY_B(par_n_HP, j, 3)
				end if
		
				if (k > 1) then
					yzel3 = gl_RAY_B(par_n_HP, j, k - 1)
				else
					yzel3 = gl_RAY_B(par_n_HP, j, N3)
				end if
		
				if (k > 2) then
					yzel33 = gl_RAY_B(par_n_HP, j, k - 2)
				else if(k == 2) then
					yzel33 = gl_RAY_B(par_n_HP, j, N3)
				else
					yzel33 = gl_RAY_B(par_n_HP, j, N3 - 1)
			end if
		
			if (k > 3) then
				yzel333 = gl_RAY_B(par_n_HP, j, k - 3)
			else if(k == 3) then
				yzel333 = gl_RAY_B(par_n_HP, j, N3)
			else if(k == 2) then
				yzel333 = gl_RAY_B(par_n_HP, j, N3 - 1)
			else
				yzel333 = gl_RAY_B(par_n_HP, j, N3 - 2)
			end if
		
			call Smooth_kvadr4(gl_x2(yzel333, now), gl_y2(yzel333, now), gl_z2(yzel333, now), &
				gl_x2(yzel33, now), gl_y2(yzel33, now), gl_z2(yzel33, now), &
				gl_x2(yzel3, now), gl_y2(yzel3, now), gl_z2(yzel3, now), &
				Ak(1), Ak(2), Ak(3), &
				gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now), &
				gl_x2(yzel22, now), gl_y2(yzel22, now), gl_z2(yzel22, now), &
				gl_x2(yzel222, now), gl_y2(yzel222, now), gl_z2(yzel222, now), &
				Ck(1), Ck(2), Ck(3))
		
			
				if (gl_Point_num(yzel) > 0) then
					!vel = par_nat_HP * 0.006 * (Bk/4.0 + Ck/4.0 + Dk/4.0 + Ek/4.0 - Ak) * gl_Point_num(yzel)/Time  ! 0.00006
					vel = kkk * par_nat_HP * 0.0002 * (Ck/2.0 + Ek/2.0 - Ak) * gl_Point_num(yzel)/Time  ! 0.00006
				else
					!vel = par_nat_HP * 0.006 * (Bk/4.0 + Ck/4.0 + Dk/4.0 + Ek/4.0 - Ak)/Time   !  0.00006
					vel = kkk * par_nat_HP * 0.0002 * (Ck/2.0 + Ek/2.0 - Ak)/Time   !  0.00006
				end if
				
				else  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					yzel = gl_RAY_B(par_n_HP, j, k)
					Ak(1) = gl_x2(yzel, now); Ak(2) = gl_y2(yzel, now); Ak(3) = gl_z2(yzel, now)
					! Ak = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/)
					r = sqrt(Ak(2)**2 + Ak(3)**2)
				
					if (j < N2) then
						yzel2 = gl_RAY_B(par_n_HP, j + 1, k)
					else
						yzel2 = gl_RAY_O(1, 1, k)
						!return
					end if
					Bk(1) = gl_x2(yzel2, now); Bk(2) = gl_y2(yzel2, now); Bk(3) = gl_z2(yzel2, now)
					! Bk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
					r1 = sqrt(Bk(2)**2 + Bk(3)**2)
				
					if (k > 1) then
						yzel2 = gl_RAY_B(par_n_HP, j, k - 1)
					else
						yzel2 = gl_RAY_B(par_n_HP, j, N3)
					end if
					Ck(1) = gl_x2(yzel2, now); Ck(2) = gl_y2(yzel2, now); Ck(3) = gl_z2(yzel2, now)
					! Ck = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
					r2 = sqrt(Ck(2)**2 + Ck(3)**2)
				
					if (j > 1) then
						yzel2 = gl_RAY_B(par_n_HP, j - 1, k)
					else
						yzel2 = gl_RAY_A(par_n_HP, size( gl_RAY_A(par_n_HP, :, k)), k)
					end if
					Dk(1) = gl_x2(yzel2, now); Dk(2) = gl_y2(yzel2, now); Dk(3) = gl_z2(yzel2, now)
					! Dk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
					r3 = sqrt(Dk(2)**2 + Dk(3)**2)
				
					if(k < N3) then
						yzel2 = gl_RAY_B(par_n_HP, j, k + 1)
					else
						yzel2 = gl_RAY_B(par_n_HP, j, 1)
					end if
					Ek(1) = gl_x2(yzel2, now); Ek(2) = gl_y2(yzel2, now); Ek(3) = gl_z2(yzel2, now)
					! Ek = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
					r4 = sqrt(Ek(2)**2 + Ek(3)**2)
				
					rr = (r1 + r2 + r3 + r4)/4.0
					!rr = (r1 + r3)/2.0
				
					!dist = sqrt( (Dk(1) - Ak(1))**2 + (Dk(2) - Ak(2))**2 + (Dk(3) - Ak(3))**2 )
					!dist = max(dist, 1.0_8)
				
					if (gl_Point_num(yzel) > 0) then
						!vel = gl_Point_num(yzel) * par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time/dist
					
						!vel = par_nat_HP * 0.0001 * gl_Point_num(yzel) * ((Ak/r) * (rr - r)) * ddt  !0.001
						!vel(1) = 0.0
					
						vel = par_nat_HP * 0.00001 * (Bk/4.0 + Ck/4.0 + Dk/4.0 + Ek/4.0 - Ak) * gl_Point_num(yzel)/Time  !0.000003 
					else
						!vel = par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time/dist
					
						!vel = par_nat_HP * 0.0001 * ((Ak/r) * (rr - r)) * ddt
						!vel(1) = 0.0
			
						vel = par_nat_HP * 0.00001 * (Bk/4.0 + Ck/4.0 + Dk/4.0 + Ek/4.0 - Ak)/Time  !0.000003
					end if
				end if
		
				
				
				gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
				gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
				gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
				
			end do
		end do
		
		
		
		N3 = size(gl_RAY_K(1, 1, :))
		N2 = size(gl_RAY_K(1, :, 1))
		N1 = size(gl_RAY_K(:, 1, 1))

		! Цикл движения точек на лучах K  
		do k = 1, N3
			do j = 1, N2
				
				if (j == 1) then
						CYCLE
				end if
				
				if (j == 1 .and. k /= 1) then
				CYCLE
				end if
		
				
		yzel = gl_RAY_K(par_n_TS, j, k)
		Ak(1) = gl_x2(yzel, now); Ak(2) = gl_y2(yzel, now); Ak(3) = gl_z2(yzel, now)
		
		
				
		if (j < N2) then
			yzel2 = gl_RAY_K(par_n_TS, j + 1, k)
		else
			yzel2 = gl_RAY_B(par_n_TS, size(gl_RAY_B(par_n_TS, :, k)), k)
		end if
		!Bk(1) = gl_x2(yzel2, now); Bk(2) = gl_y2(yzel2, now); Bk(3) = gl_z2(yzel2, now)
		
		if (j < N2 - 1) then
			yzel22 = gl_RAY_K(par_n_TS, j + 2, k)
		else if (j == N2 - 1) then
			yzel22 = gl_RAY_B(par_n_TS, size(gl_RAY_B(par_n_TS, :, k)), k)
		else
			yzel22 = gl_RAY_B(par_n_TS, size(gl_RAY_B(par_n_TS, :, k)) - 1, k)
		end if
		
		if(j < N2 - 2) then
				yzel222 = gl_RAY_K(par_n_TS, j + 3, k)
			else if(j == N2 - 2) then
				yzel222 = gl_RAY_B(par_n_TS, size(gl_RAY_B(par_n_TS, :, k)), k)
			else if(j == N2 - 1) then
				yzel222 = gl_RAY_B(par_n_TS, size(gl_RAY_B(par_n_TS, :, k)) - 1, k)
			else
				yzel222 = gl_RAY_B(par_n_TS, size(gl_RAY_B(par_n_TS, :, k)) - 2, k)
		end if
		
		if (j > 1) then
			yzel3 = gl_RAY_K(par_n_TS, j - 1, k)
		else
			yzel3 = gl_RAY_K(par_n_TS, 1, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
		end if
		!Dk(1) = gl_x2(yzel2, now); Dk(2) = gl_y2(yzel2, now); Dk(3) = gl_z2(yzel2, now)
		
		if (j > 2) then
			yzel33 = gl_RAY_K(par_n_TS, j - 2, k)
		else if(j == 2) then
			yzel33 = gl_RAY_K(par_n_TS, 1, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
		else
			yzel33 = gl_RAY_K(par_n_TS, 2, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
		end if
		
		if (j > 3) then
			yzel333 = gl_RAY_K(par_n_TS, j - 3, k)
		else if(j == 3) then
			yzel333 = gl_RAY_K(par_n_TS, 1, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
		else if(j == 2) then
			yzel333 = gl_RAY_K(par_n_TS, 2, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
		else
			yzel333 = gl_RAY_K(par_n_TS, 3, mod(k + ceiling(1.0 * N3/2) - 1, N3) + 1)
		end if
		
		call Smooth_kvadr4(gl_x2(yzel333, now), gl_y2(yzel333, now), gl_z2(yzel333, now), &
				gl_x2(yzel33, now), gl_y2(yzel33, now), gl_z2(yzel33, now), &
				gl_x2(yzel3, now), gl_y2(yzel3, now), gl_z2(yzel3, now), &
				Ak(1), Ak(2), Ak(3), &
				gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now), &
				gl_x2(yzel22, now), gl_y2(yzel22, now), gl_z2(yzel22, now), &
				gl_x2(yzel222, now), gl_y2(yzel222, now), gl_z2(yzel222, now), &
				Ek(1), Ek(2), Ek(3))
		
				
		if (k > 1) then
			yzel3 = gl_RAY_K(par_n_TS, j, k - 1)
		else
			yzel3 = gl_RAY_K(par_n_TS, j, N3)
		end if
		!Ck(1) = gl_x2(yzel2, now); Ck(2) = gl_y2(yzel2, now); Ck(3) = gl_z2(yzel2, now)
				
		if (k > 2) then
			yzel33 = gl_RAY_K(par_n_TS, j, k - 2)
		else if(k == 2) then
			yzel33 = gl_RAY_K(par_n_TS, j, N3)
		else
			yzel33 = gl_RAY_K(par_n_TS, j, N3 - 1)
		end if
		
		if (k > 3) then
			yzel333 = gl_RAY_K(par_n_TS, j, k - 3)
		else if(k == 3) then
			yzel333 = gl_RAY_K(par_n_TS, j, N3)
		else if(k == 2) then
			yzel333 = gl_RAY_K(par_n_TS, j, N3 - 1)
		else
			yzel333 = gl_RAY_K(par_n_TS, j, N3 - 2)
		end if
				
		if(k < N3) then
			yzel2 = gl_RAY_K(par_n_TS, j, k + 1)
		else
			yzel2 = gl_RAY_K(par_n_TS, j, 1)
		end if
		!Ek(1) = gl_x2(yzel2, now); Ek(2) = gl_y2(yzel2, now); Ek(3) = gl_z2(yzel2, now)
		
		if(k < N3 - 1) then
			yzel22 = gl_RAY_K(par_n_TS, j, k + 2)
		else if(k == N3 - 1) then
			yzel22 = gl_RAY_K(par_n_TS, j, 1)
		else
			yzel22 = gl_RAY_K(par_n_TS, j, 2)
		end if
		
		if(k < N3 - 2) then
			yzel222 = gl_RAY_K(par_n_TS, j, k + 2)
		else if(k == N3 - 2) then
			yzel222 = gl_RAY_K(par_n_TS, j, 1)
		else if(k == N3 - 1) then
			yzel222 = gl_RAY_K(par_n_TS, j, 2)
		else
			yzel222 = gl_RAY_K(par_n_TS, j, 3)
		end if
		
		call Smooth_kvadr4(gl_x2(yzel333, now), gl_y2(yzel333, now), gl_z2(yzel333, now), &
		gl_x2(yzel33, now), gl_y2(yzel33, now), gl_z2(yzel33, now), &
		gl_x2(yzel3, now), gl_y2(yzel3, now), gl_z2(yzel3, now), &
		Ak(1), Ak(2), Ak(3), &
		gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now), &
		gl_x2(yzel22, now), gl_y2(yzel22, now), gl_z2(yzel22, now), &
		gl_x2(yzel222, now), gl_y2(yzel222, now), gl_z2(yzel222, now), &
		Ck(1), Ck(2), Ck(3))
				
		if (gl_Point_num(yzel) > 0) then
			vel = 10.0 * gl_Point_num(yzel) * par_nat_TS * (Ck/2.0 + Ek/2.0 - Ak)/Time * ddt
		else
			vel = 10.0 * par_nat_TS * (Ck/2.0 + Ek/2.0 - Ak)/Time * ddt
		end if
			
				
		gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
		gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
		gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
				
			end do
		end do
		
		end if
		
		N3 = size(gl_RAY_O(1, 1, :))
		N2 = size(gl_RAY_O(1, :, 1))
		N1 = size(gl_RAY_O(:, 1, 1))

		! Цикл движения точек на лучах O  
		if (.True.) then   ! Можно отключить этот цикл
			do k = 1, N3
				do j = 1, N2
				
				! Контакт
		yzel = gl_RAY_O(1, j, k)
		Ak(1) = gl_x2(yzel, now); Ak(2) = gl_y2(yzel, now); Ak(3) = gl_z2(yzel, now)
		! Ak = (/gl_x2(yzel, now), gl_y2(yzel, now), gl_z2(yzel, now)/)
		r = sqrt(Ak(2)**2 + Ak(3)**2)
				
		if (j < N2) then
			yzel2 = gl_RAY_O(1, j + 1, k)
		else
			yzel2 = yzel
		end if
		Bk(1) = gl_x2(yzel2, now); Bk(2) = gl_y2(yzel2, now); Bk(3) = gl_z2(yzel2, now)
		! Bk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
		r1 = sqrt(Bk(2)**2 + Bk(3)**2)
				
		if (k > 1) then
			yzel2 = gl_RAY_O(1, j, k - 1)
		else
			yzel2 = gl_RAY_O(1, j, N3)
		end if
		Ck(1) = gl_x2(yzel2, now); Ck(2) = gl_y2(yzel2, now); Ck(3) = gl_z2(yzel2, now)
		! Ck = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
		r2 = sqrt(Ck(2)**2 + Ck(3)**2)
				
		if (j > 1) then
			yzel2 = gl_RAY_O(1, j - 1, k)
		else
			yzel2 = gl_RAY_C(1, size(gl_RAY_C(1, :, k)), k)
		end if
		Dk(1) = gl_x2(yzel2, now); Dk(2) = gl_y2(yzel2, now); Dk(3) = gl_z2(yzel2, now)
		! Dk = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
		r3 = sqrt(Dk(2)**2 + Dk(3)**2)
				
		if(k < N3) then
			yzel2 = gl_RAY_O(1, j, k + 1)
		else
			yzel2 = gl_RAY_O(1, j, 1)
		end if
		Ek(1) = gl_x2(yzel2, now); Ek(2) = gl_y2(yzel2, now); Ek(3) = gl_z2(yzel2, now)
		! Ek = (/gl_x2(yzel2, now), gl_y2(yzel2, now), gl_z2(yzel2, now)/)
		r4 = sqrt(Ek(2)**2 + Ek(3)**2)
		
		!dist = sqrt( (Dk(1) - Ak(1))**2 + (Dk(2) - Ak(2))**2 + (Dk(3) - Ak(3))**2 )
		!dist = max(dist, 1.0_8)
		
		rr = (r1 + r2 + r3 + r4)/4.0
		
		!rr = (r1 + r3)/2.0
		
		
		
				
		if (gl_Point_num(yzel) > 0) then
			!vel = gl_Point_num(yzel) * par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time
			
			!vel = par_nat_HP * 0.0001 * gl_Point_num(yzel) * ((Ak/r) * (rr - r)) * ddt   ! 0.0003
			!vel(1) = 0.0
			
			vel = par_nat_HP * 0.0006 * (Bk/4.0 + Ck/4.0 + Dk/4.0 + Ek/4.0 - Ak) * gl_Point_num(yzel)/Time  ! надо ещё уменьшать
			
		else
			!vel = par_nat_HP * (Bk/8.0 + Ck/8.0 + Dk/8.0 + Ek/8.0 - Ak/2.0)/Time
			
			!vel = par_nat_HP * 0.0001 * ((Ak/r) * (rr - r)) * ddt
			!vel(1) = 0.0
			
			vel = par_nat_HP * 0.0006 * (Bk/4.0 + Ck/4.0 + Dk/4.0 + Ek/4.0 - Ak)/Time  !0.00006
		end if
		
		!Ak = Bk/4.0 + Ck/4.0 + Dk/4.0 + Ek/4.0 - Ak
		!nd = norm2(Ak)
		!vel = 10.0 * gl_Vn * (Ak)/nd
				
		gl_Vx(yzel) = gl_Vx(yzel) + vel(1)
		gl_Vy(yzel) = gl_Vy(yzel) + vel(2)
		gl_Vz(yzel) = gl_Vz(yzel) + vel(3)
				
				end do
			end do
		end if
		
		
		! Конец поверхностного натяжения --------------------------------------------------------------------------------------------------------------------
		
		!print*, gl_x2(gl_RAY_B(3, 1, 1), now2)
		!print*, gl_y2(gl_RAY_B(6, 4, 4), now2)
		!print*, gl_z2(gl_RAY_B(5, 6, 5), now2)
		!print*, "---------------------------------------------------"
		
			
		N3 = size(gl_RAY_A(1, 1, :))
		N2 = size(gl_RAY_A(1, :, 1))
		N1 = size(gl_RAY_A(:, 1, 1))

		! Цикл движения точек на лучах А  ************************************************************
		do k = 1, N3
			do j = 1, N2
				
				if (k /= 1 .and. j == 1) then
						CYCLE
				end if
				
				! Вычисляем координаты текущего луча в пространстве
				the = (j - 1) * par_pi_8/2.0/(N2 - 1)
				phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
				
				! TS
				yzel = gl_RAY_A(par_n_TS, j, k)
				if(gl_Point_num(yzel) == 0) then
					vel = 0.0
				else
					! vel = (/gl_Vx(yzel), gl_Vy(yzel), gl_Vz(yzel)/)
					vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
					vel = vel/gl_Point_num(yzel)                       ! Нашли скорость движения этого узла
				end if
				
				! Обнулим для следующего использования
				gl_Point_num(yzel) = 0
				gl_Vx(yzel) = 0.0
				gl_Vy(yzel) = 0.0
				gl_Vz(yzel) = 0.0
				
				ER(1) = cos(the); ER(2) = sin(the) * cos(phi); ER(3) = sin(the) * sin(phi)
				proect = DOT_PRODUCT(vel * Time, ER)
				! proect = DOT_PRODUCT(vel * Time, (/cos(the), sin(the) * cos(phi), sin(the) * sin(phi)/))  !  Находим проекцию перемещения на радиус вектор луча
				KORD(1) = gl_x2(yzel, now); KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now) 
				R_TS = norm2(KORD + proect * ER)  ! Новое расстояние до TS
				
				!if(j == 1 .and. k==1) write(*,*) "=    ", gl_x2(yzel, now)
				
				! HP
				yzel = gl_RAY_A(par_n_HP, j, k)
				if(gl_Point_num(yzel) == 0) then
					vel = 0.0
				else
					! vel = (/gl_Vx(yzel), gl_Vy(yzel), gl_Vz(yzel)/)
					vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
					vel = vel/gl_Point_num(yzel)                       ! Нашли скорость движения этого узла
				end if
				
				! Обнулим для следующего успользования
				gl_Point_num(yzel) = 0
				gl_Vx(yzel) = 0.0
				gl_Vy(yzel) = 0.0
				gl_Vz(yzel) = 0.0
				
				proect = DOT_PRODUCT(vel * Time, ER)
				KORD(1) = gl_x2(yzel, now); KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now) 
				R_HP = norm2(KORD + proect * ER)  ! Новое расстояние до HP
				
				
				
				! BS
				yzel = gl_RAY_A(par_n_BS, j, k)
				if(gl_Point_num(yzel) == 0) then
					vel = 0.0
				else
					! vel = (/gl_Vx(yzel), gl_Vy(yzel), gl_Vz(yzel)/)
					vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
					vel = vel/gl_Point_num(yzel)                       ! Нашли скорость движения этого узла
				end if
				
				
				! Обнулим для следующего использования
				gl_Point_num(yzel) = 0
				gl_Vx(yzel) = 0.0
				gl_Vy(yzel) = 0.0
				gl_Vz(yzel) = 0.0
				
				proect = DOT_PRODUCT(vel * Time, ER)
				KORD(1) = gl_x2(yzel, now); KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now) 
				R_BS = norm2(KORD + proect * ER)  ! Новое расстояние до BS
				
				
				! Далее обычный цикл нахождения координат точек, такой же, как и при построении сетки
				do i = 1, N1
					
					if (i == 1) CYCLE

					yzel = gl_RAY_A(i, j, k)
					! Вычисляем координаты точки на луче
					
					kk13 = par_kk13 * (par_pi_8/2.0 - the)/(par_pi_8/2.0)  +  1.0 * (the/(par_pi_8/2.0))**2

					! до TS
					if (i <= par_n_IB) then  ! NEW
							if(i == 2) then
								r =  par_R0 - (par_R_inner - par_R0) * (DBLE(3 - 2)/(par_n_IB - 2))**par_kk1
								if(r < 0.0) then
									print*, "Error iouihjgfdcydygy  ", r
									STOP
								end if
							else
								r =  par_R0 + (par_R_inner - par_R0) * (DBLE(i - 2)/(par_n_IB - 2))**par_kk1
							end if
							!r =  par_R0 + (par_R_inner - par_R0) * (DBLE(i - 2)/(par_n_IB - 2))**par_kk1
					else if (i <= par_n_TS) then  
							r =  par_R_inner + (R_TS - par_R_inner) * (DBLE(i - par_n_IB)/(par_n_TS - par_n_IB))**par_kk12
					!if (i <= par_n_TS) then  ! До расстояния = R_TS
					!    r =  par_R0 + (R_TS - par_R0) * (DBLE(i)/par_n_TS)**par_kk1
					else if (i <= par_n_HP) then  
						r = R_TS + (i - par_n_TS) * (R_HP - R_TS)/(par_n_HP - par_n_TS)
					else if (i == par_n_HP + 1) then 
						rd = R_TS + (par_n_HP - 1 - par_n_TS) * (R_HP - R_TS)/(par_n_HP - par_n_TS)
						r = R_HP + (R_HP - rd)
					else if (i <= par_n_BS) then 
						rd = R_TS + (par_n_HP - 1 - par_n_TS) * (R_HP - R_TS)/(par_n_HP - par_n_TS)
						rd = R_HP + (R_HP - rd)
						r = rd + (R_BS - rd) * sgushenie_1( 1.0_8 * (i - par_n_HP - 1)/(par_n_BS - par_n_HP - 1), kk13 )
						!r = R_HP + (R_BS - R_HP) * angle_cilindr( 1.0_8 * (i - par_n_HP)/(par_n_BS - par_n_HP), kk13 )
						!r = R_HP + (i - par_n_HP) * (R_BS - R_HP)/(par_n_BS - par_n_HP)
					else
						r = R_BS + (par_R_END - R_BS) * (DBLE(i- par_n_BS)/(par_n_END - par_n_BS))**(par_kk2 * (0.55 + 0.45 * cos(the)) )
					end if
					
					if (i == par_n_TS - 1) then
						r = R_TS - 0.5               
					end if

					! Записываем новые координаты
					gl_x2(yzel, now2) = r * cos(the)
					gl_y2(yzel, now2) = r * sin(the) * cos(phi)
					gl_z2(yzel, now2) = r * sin(the) * sin(phi)
					
					!if(j == 1 .and. k==1 .and. i==20)  write(*,*) "=    ", gl_x2(yzel, now2)
					
				end do
			end do
		end do
		
			
		N3 = size(gl_RAY_B(1, 1, :))
		N2 = size(gl_RAY_B(1, :, 1))
		N1 = size(gl_RAY_B(:, 1, 1))

		! Цикл движения точек на лучах B  ************************************************************
		do k = 1, N3
			do j = 1, N2
				
				! Вычисляем координаты текущего луча в пространстве
				the = par_pi_8/2.0 + (j) * par_triple_point/(N2)
				the2 = par_pi_8/2.0 + (j) * par_triple_point_2/(N2)
				phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
				
				! TS
				yzel = gl_RAY_B(par_n_TS, j, k)
				if(gl_Point_num(yzel) == 0) then
					vel = 0.0
				else
					vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
					vel = vel/gl_Point_num(yzel)                       ! Нашли скорость движения этого узла
				end if
				
				! Обнулим для следующего использования
				gl_Point_num(yzel) = 0
				gl_Vx(yzel) = 0.0
				gl_Vy(yzel) = 0.0
				gl_Vz(yzel) = 0.0
				
				ER(1) = cos(the); ER(2) = sin(the) * cos(phi); ER(3) = sin(the) * sin(phi)
				proect = DOT_PRODUCT(vel * Time, ER)  !  Находим проекцию перемещения на радиус вектор луча
				KORD(1) = gl_x2(yzel, now); KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now) 
				R_TS = norm2(KORD + proect * ER)  ! Новое расстояние до TS
				
				! HP
				yzel = gl_RAY_B(par_n_HP, j, k)
				if(gl_Point_num(yzel) == 0) then
					vel = 0.0
				else
					vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
					vel = vel/gl_Point_num(yzel)                       ! Нашли скорость движения этого узла
				end if
				
				
				! Обнулим для следующего использования
				gl_Point_num(yzel) = 0
				gl_Vx(yzel) = 0.0
				gl_Vy(yzel) = 0.0
				gl_Vz(yzel) = 0.0
				
				ER2(1) = cos(the2); ER2(2) = sin(the2) * cos(phi); ER2(3) = sin(the2) * sin(phi)
				proect = DOT_PRODUCT(vel * Time, ER2)
				KORD(1) = gl_x2(yzel, now); KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now) 
				R_HP = norm2(KORD + proect * ER2 - (R_TS * ER))  ! Новое расстояние до HP от края TS
				
				
					
				do i = 1, N1

					if (i == 1) CYCLE
					
					yzel = gl_RAY_B(i, j, k)
					! Вычисляем координаты точки на луче

					
					! до TS
					if (i <= par_n_IB) then  ! NEW
						if(i == 2) then
							r =  par_R0 - (par_R_inner - par_R0) * (DBLE(3 - 2)/(par_n_IB - 2))**par_kk1
							if(r < 0.0) then
								print*, "Error iouihjgfdcydygy  ", r
								STOP
							end if
						else
							r =  par_R0 + (par_R_inner - par_R0) * (DBLE(i - 2)/(par_n_IB - 2))**par_kk1
						end if
						!r =  par_R0 + (par_R_inner - par_R0) * (DBLE(i - 2)/(par_n_IB - 2))**par_kk1
					else if (i <= par_n_TS) then  ! До расстояния = R_TS
						r =  par_R_inner + (R_TS - par_R_inner) * (DBLE(i - par_n_IB)/(par_n_TS - par_n_IB))**par_kk12
					!if (i <= par_n_TS) then  ! До расстояния = R_TS
					!    r =  par_R0 + (R_TS - par_R0) * (REAL(i, KIND = 4)/par_n_TS)**par_kk1
					else if (i <= par_n_HP) then  ! До расстояния = par_R_character * 1.3
						!r = R_TS + (i - par_n_TS) * (R_HP - R_TS) /(par_n_HP - par_n_TS)
						
						r1 = par_R_inner + (R_TS - par_R_inner)
						r =  (i - par_n_TS) * (R_HP) /(par_n_HP - par_n_TS)
						gl_x2(yzel, now2) = r1 * cos(the) + r * cos(the2)
						gl_y2(yzel, now2) = r1 * sin(the) * cos(phi) + r * sin(the2) * cos(phi)
						gl_z2(yzel, now2) = r1 * sin(the) * sin(phi) + r * sin(the2) * sin(phi)
						CYCLE
					end if
					
					if (i == par_n_TS - 1) then
						r = R_TS - 0.5               
					end if
					

					! Записываем новые координаты
					gl_x2(yzel, now2) = r * cos(the)
					gl_y2(yzel, now2) = r * sin(the) * cos(phi)
					gl_z2(yzel, now2) = r * sin(the) * sin(phi)
					
				
				end do
			end do
		end do
		
		
		N3 = size(gl_RAY_C(1, 1, :))
		N2 = size(gl_RAY_C(1, :, 1))
		N1 = size(gl_RAY_C(:, 1, 1))

		! Цикл движения точек на лучах C  ************************************************************
		do k = 1, N3
			do j = 1, N2
				
				! Вычисляем координаты текущего луча в пространстве
					phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
					
					! Вычисляем координаты точки на луче
					x = gl_x2(gl_RAY_B(par_n_HP, j, k), now2)
					y = gl_y2(gl_RAY_B(par_n_HP, j, k), now2)
					z = gl_z2(gl_RAY_B(par_n_HP, j, k), now2)
					rr = (y**2 + z**2)**(0.5)
					
					
					y = gl_y2(gl_RAY_B(par_n_HP - 1, j, k), now2)
					z = gl_z2(gl_RAY_B(par_n_HP - 1, j, k), now2)
					rd = (y**2 + z**2)**(0.5)
					rr = rr + (rr - rd)
					
					! BS     Нужно взять положение BS из её положения на крайнем луче A
					yzel = gl_RAY_A(par_n_BS, size(gl_RAY_A(1, :, 1)), k)
					ER(1) = 0.0_8; ER(2) = gl_y2(yzel, now2); ER(3) = gl_z2(yzel, now2)
					R_BS = norm2(ER)  ! Новое расстояние до BS
					
				
				do i = 1, N1
					
					if(i == 1) CYCLE
					
					yzel = gl_RAY_C(i, j, k)

					if(i == 2) then
						r = rr
					else if (i <= par_n_BS - par_n_HP + 1) then
						r = rr + (i - 2) * (R_BS - rr)/(par_n_BS - par_n_HP - 1)
					else
						r = R_BS + (DBLE(i - (par_n_BS - par_n_HP + 1))/(N1 - (par_n_BS - par_n_HP + 1) ))**(0.55 * par_kk2) * (par_R_END - R_BS)
					end if
					
					gl_x2(yzel, now2) = x
					gl_y2(yzel, now2) = r * cos(phi)
					gl_z2(yzel, now2) = r * sin(phi)
					

				end do
			end do
		end do

		N3 = size(gl_RAY_O(1, 1, :))
		N2 = size(gl_RAY_O(1, :, 1))
		N1 = size(gl_RAY_O(:, 1, 1))


		! Цикл движения точек на лучах O   ************************************************************
		do k = 1, N3
			phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
			
			do j = 1, N2
				
				yzel = gl_RAY_O(1, j, k)
				if(gl_Point_num(yzel) == 0) then
					vel = 0.0
				else
					vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
					vel = vel/gl_Point_num(yzel)                       ! Нашли скорость движения этого узла
				end if
				
				! Обнулим для следующего использования
				gl_Point_num(yzel) = 0
				gl_Vx(yzel) = 0.0_8
				gl_Vy(yzel) = 0.0_8
				gl_Vz(yzel) = 0.0_8
				
				ER(1) = 0.0_8; ER(2) = cos(phi); ER(3) = sin(phi)
				proect = DOT_PRODUCT(vel * Time, ER)  !  Находим проекцию перемещения на радиус вектор луча
				KORD(1) = 0.0_8; KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now)
				R_HP = norm2(KORD + proect * ER )  ! Новое расстояние до HP
				
				
				! Блокируем схлопывание контакта к оси
				if(R_HP < 20.0_8) then
					R_HP = 20.0_8
				end if
				
				
				xx = gl_x2(gl_RAY_B(par_n_HP, par_m_BC, k), now2)              ! Отталкиваемся от x - координаты крайней точки B на гелиопаузе в этой плоскости (k)
				x = xx - (DBLE(j)/N2)**par_kk31 * (xx - par_R_LEFT)
				
				! BS     Нужно взять положение BS из её положения на крайнем луче A
				yzel = gl_RAY_A(par_n_BS, size(gl_RAY_A(1, :, 1)), k)
				KORD(1) = 0.0_8; KORD(2) = gl_y2(yzel, now2); KORD(3) = gl_z2(yzel, now2)
				R_BS = norm2(KORD)  ! Новое расстояние до BS
				
				if (R_BS < R_HP + 10.0) then
					R_BS = R_HP + 10.0
					!print*, "R_BS = R_HP + 10.0"
				end if
				
				
				do i = 1, N1
					yzel = gl_RAY_O(i, j, k)
					
					if (i <= par_n_BS - par_n_HP + 1) then
						r = R_HP + (i - 1) * (R_BS - R_HP)/(par_n_BS - par_n_HP)
					else
						r = R_BS + (DBLE(i - (par_n_BS - par_n_HP + 1))/(N1 - (par_n_BS - par_n_HP + 1) ))**(0.55 * par_kk2) * (par_R_END - R_BS)
					end if


					gl_x2(yzel, now2) = x
					gl_y2(yzel, now2) = r * cos(phi)
					gl_z2(yzel, now2) = r * sin(phi)
					
				end do
			end do
		end do
		
		N3 = size(gl_RAY_K(1, 1, :))
		N2 = size(gl_RAY_K(1, :, 1))
		N1 = size(gl_RAY_K(:, 1, 1))

		! Цикл движения точек на лучах K  ************************************************************
		do k = 1, N3
			do j = 1, N2
				
				! Вычисляем координаты текущего луча в пространстве
				the = par_pi_8/2.0 + par_triple_point + (N2 - j + 1) * (par_pi_8/2.0 - par_triple_point)/(N2)
				phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
				
				if (k /= 1 .and. j == 1) CYCLE
				
				yzel = gl_RAY_K(par_n_TS, j, k)
				if(gl_Point_num(yzel) == 0) then
					vel = 0.0
				else
					vel(1) = gl_Vx(yzel); vel(2) = gl_Vy(yzel); vel(3) = gl_Vz(yzel)
					vel = vel/gl_Point_num(yzel)                       ! Нашли скорость движения этого узла
				end if
				
				! Обнулим для следующего использования
				gl_Point_num(yzel) = 0
				gl_Vx(yzel) = 0.0
				gl_Vy(yzel) = 0.0
				gl_Vz(yzel) = 0.0
				
				ER(1) = cos(the); ER(2) = sin(the) * cos(phi); ER(3) = sin(the) * sin(phi)
				proect = DOT_PRODUCT(vel * Time, ER)  !  Находим проекцию перемещения на радиус вектор луча
				KORD(1) = gl_x2(yzel, now); KORD(2) = gl_y2(yzel, now); KORD(3) = gl_z2(yzel, now)
				R_TS = norm2(KORD + proect * ER)  ! Новое расстояние до TS
				
				do i = 1, N1

					if (i == 1) CYCLE
					
					! Вычисляем координаты точки на луче
					yzel = gl_RAY_K(i, j, k)
					
					if (i <= par_n_IB) then  ! NEW
						if(i == 2) then
							r =  par_R0 - (par_R_inner - par_R0) * (DBLE(3 - 2)/(par_n_IB - 2))**par_kk1
							if(r < 0.0) then
								print*, "Error iouihjgfdcydygy  ", r
								STOP
							end if
						else
							r =  par_R0 + (par_R_inner - par_R0) * (DBLE(i - 2)/(par_n_IB - 2))**par_kk1
						end if
						!r =  par_R0 + (par_R_inner - par_R0) * (DBLE(i - 2)/(par_n_IB - 2))**par_kk1
					else 
						r =  par_R_inner + (R_TS - par_R_inner) * (DBLE(i - par_n_IB)/(par_n_TS - par_n_IB))**par_kk12
					end if
					
					if (i == par_n_TS - 1) then
						r = R_TS - 0.5              
					end if
						
					!r =  par_R0 + (R_TS - par_R0) * (REAL(i, KIND = 4)/par_n_TS)**par_kk1


					! Записываем новые координаты
					gl_x2(yzel, now2) = r * cos(the)
					gl_y2(yzel, now2) = r * sin(the) * cos(phi)
					gl_z2(yzel, now2) = r * sin(the) * sin(phi)
					
					
				end do
			end do
		end do
		
		
		N3 = size(gl_RAY_D(1, 1, :))
		N2 = size(gl_RAY_D(1, :, 1))
		N1 = size(gl_RAY_D(:, 1, 1))

		! Цикл движения точек на лучах D ************************************************************
		
		
		
		do k = 1, N3
			do j = 1, N2
				
				! Вычисляем координаты текущего луча в пространстве
					phi = 2.0_8 * par_pi_8 * angle_cilindr((k - 1.0_8)/(N3), par_al1)  !(k - 1) * 2.0_8 * par_pi_8/(N3)
					
					! Вычисляем координаты точки на луче

					if (j < N2) then
						xx = gl_x2(gl_RAY_K(par_n_TS, j, k), now2)
						y = gl_y2(gl_RAY_K(par_n_TS, j, k), now2)
						z = gl_z2(gl_RAY_K(par_n_TS, j, k), now2)
					else
						xx = gl_x2(gl_RAY_B(par_n_TS, par_m_BC, k), now2)
						y = gl_y2(gl_RAY_B(par_n_TS, par_m_BC, k), now2)
						z = gl_z2(gl_RAY_B(par_n_TS, par_m_BC, k), now2)
					end if
					
					!y = gl_y2(gl_RAY_B(par_n_TS, par_m_BC, k), now2)
					!z = gl_z2(gl_RAY_B(par_n_TS, par_m_BC, k), now2)
					
					r = sqrt(y**2 + z**2)
					
					
				do i = 1, N1

					if (i == 1) CYCLE

					if (k /= 1 .and. j == 1) CYCLE
					
					
					!y = gl_y2(gl_RAY_O(1, i - 1, k), now2)
					!z = gl_z2(gl_RAY_O(1, i - 1, k), now2)
					!
					!rr = sqrt(y**2 + z**2) * 0.4

					yzel = gl_RAY_D(i, j, k)
					
					
					x = xx + (DBLE(i - 1)/(N1 - 1))**par_kk3 * (par_R_LEFT - xx)
					
					!rrr = r * (DBLE(j - 1)/(N2 - 1))
					!rrr = r * (1.0 - (DBLE(i - 1)**2/(N1 - 1)**2)) + (rr * (DBLE(j - 1)/(N2 - 1))) * ( (DBLE(i - 1)**2) /(N1 - 1)**2)

					
					
					! Записываем новые координаты
					gl_x2(yzel, now2) = x
					gl_y2(yzel, now2) = r * cos(phi)
					gl_z2(yzel, now2) = r * sin(phi)
					
				end do
			end do
		end do
		
		N3 = size(gl_RAY_E(1, 1, :))
		N2 = size(gl_RAY_E(1, :, 1))
		N1 = size(gl_RAY_E(:, 1, 1))
		
		! Цикл движения точек на лучах E  ************************************************************
		do k = 1, N3
			do j = 1, N2
				do i = 1, N1

					if (i == 1) CYCLE

					if (i == N1) CYCLE

					x = gl_x2(gl_RAY_E(1, j, k), now2)
					y = gl_y2(gl_RAY_E(1, j, k), now2)
					z = gl_z2(gl_RAY_E(1, j, k), now2)

					x2 = gl_x2(gl_RAY_O(1, j, k), now2)
					y2 = gl_y2(gl_RAY_O(1, j, k), now2)
					z2 = gl_z2(gl_RAY_O(1, j, k), now2)

					yzel = gl_RAY_E(i, j, k)

					gl_x2(yzel, now2) = x + (x2 - x) * (i - 1)/(N1 - 1)
					gl_y2(yzel, now2) = y + (y2 - y) * (i - 1)/(N1 - 1)
					gl_z2(yzel, now2) = z + (z2 - z) * (i - 1)/(N1 - 1)
					
					!if(i == 14 .and. j == 12 .and. k == 35) then
					!	print*, "___"
					!	print*, x, x2
					!	print*, y, y2
					!	print*, z, z2
					!	print*, gl_x2(yzel, now2), gl_x2(yzel, now)
					!	print*, gl_y2(yzel, now2), gl_y2(yzel, now)
					!	print*, gl_z2(yzel, now2), gl_z2(yzel, now)
					!	print*, gl_y2(gl_RAY_E(1, j, k), now2), gl_y2(gl_RAY_E(1, j, k), now)
					!	print*, gl_z2(gl_RAY_E(1, j, k), now2), gl_z2(gl_RAY_E(1, j, k), now)
					!	print*, "___"
					!	print*, gl_x2(gl_RAY_O(1, j, k), now2), gl_x2(gl_RAY_O(1, j, k), now)
					!	print*, gl_y2(gl_RAY_O(1, j, k), now2), gl_y2(gl_RAY_O(1, j, k), now)
					!	print*, gl_z2(gl_RAY_O(1, j, k), now2), gl_z2(gl_RAY_O(1, j, k), now)
					!	print*, "___"
					!end if
					
					
				end do
			end do
		end do
		!
		!
		!print*, gl_x2(gl_RAY_A(10, 1, 1), now2)
		!print*, gl_y2(gl_RAY_A(14, 4, 4), now2)
		!print*, gl_z2(gl_RAY_A(18, 8, 5), now2)
		!print*, "---------------------------------------------------"
		!
		!print*, gl_x2(gl_RAY_B(3, 1, 1), now2)
		!print*, gl_y2(gl_RAY_B(6, 4, 4), now2)
		!print*, gl_z2(gl_RAY_B(5, 6, 5), now2)
		!print*, "---------------------------------------------------"
		!
		!print*, gl_x2(gl_RAY_C(3, 1, 1), now2)
		!print*, gl_y2(gl_RAY_C(6, 4, 4), now2)
		!print*, gl_z2(gl_RAY_C(5, 6, 5), now2)
		!print*, "---------------------------------------------------"
		!
		!print*, gl_x2(gl_RAY_O(3, 1, 1), now2)
		!print*, gl_y2(gl_RAY_O(6, 4, 4), now2)
		!print*, gl_z2(gl_RAY_O(5, 6, 5), now2)
		!print*, "---------------------------------------------------"
		!
		!print*, gl_x2(gl_RAY_K(3, 1, 1), now2)
		!print*, gl_y2(gl_RAY_K(6, 4, 4), now2)
		!print*, gl_z2(gl_RAY_K(5, 6, 5), now2)
		!print*, "---------------------------------------------------"
		!
		!print*, gl_x2(gl_RAY_D(3, 1, 1), now2)
		!print*, gl_y2(gl_RAY_D(6, 4, 4), now2)
		!print*, gl_z2(gl_RAY_D(5, 6, 5), now2)
		!print*, "---------------------------------------------------"
		!
		!print*, gl_x2(gl_RAY_E(3, 1, 1), now2)
		!print*, gl_y2(gl_RAY_E(6, 4, 4), now2)
		!print*, gl_z2(gl_RAY_E(5, 6, 5), now2)
		!print*, "---------------------------------------------------"
		!pause
		
		
    end subroutine Move_all
    
    subroutine Get_center_move(n, CENTER, now)  ! Получить массив - центр ячейки
		use STORAGE, only: gl_all_Cell, gl_x2, gl_y2, gl_z2
		implicit none

		integer, intent(in) :: n, now
		real(8), intent(out) :: CENTER(3)
		real(8) :: c(3)
		integer :: i, j

		CENTER = 0.0

		do i = 1, 8
			j = gl_all_Cell(i, n)
			c = (/gl_x2(j, now), gl_y2(j, now), gl_z2(j, now)/)
			CENTER = CENTER + c
		end do

		CENTER = CENTER/8.0

    end subroutine Get_center_move
    
    subroutine calc_all_Gran_move(now)   ! Программа расчёта площадей и нормалей граней и Объёмов ячеек
    ! Версия функции для расчёта на подвижной сетке
    use STORAGE
    use GEO_PARAM
    implicit none
    
    integer, intent(in) :: now           ! Откуда берём координаты узлов и куда записываем вычисленные параметры (площади и объёмы)
    integer :: Ngran, iter
    real(8) :: p(3, 4), Vol, D
    real(8) :: a(3), b(3), c(3), S, node1(3), node2(3)
    real(8) :: dist, di, gr_center(3)
    integer :: i, j, k, ll, grc, ii, now2

    !print*, "calc_all_Gran_move"
    now2 = mod(now, 2) + 1 
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
            p(:,j) = (/gl_x2(i, now), gl_y2(i, now), gl_z2(i, now)/)
        end do
        a = p(:,3) - p(:,1)
        b = p(:,4) - p(:,2)

        ! Нужно сохранить центр грани
		gr_center = p(:, 1) + p(:, 2) + p(:, 3) + p(:, 4)
		grc = 4
        !gr_center = p(:, 1)
        !grc = 1
        !if (gl_all_Gran(2, iter) /= gl_all_Gran(1, iter)) then
        !    grc = grc + 1
        !    gr_center = gr_center + p(:,2)
        !end if
        !
        !if (gl_all_Gran(3, iter) /= gl_all_Gran(2, iter) .and. gl_all_Gran(3, iter) /= gl_all_Gran(1, iter)) then
        !    grc = grc + 1
        !    gr_center = gr_center + p(:,3)
        !end if
        !
        !if (gl_all_Gran(4, iter) /= gl_all_Gran(1, iter) .and. gl_all_Gran(4, iter) /= gl_all_Gran(2, iter) &
        !    .and. gl_all_Gran(4, iter) /= gl_all_Gran(3, iter)) then
        !    grc = grc + 1
        !    gr_center = gr_center + p(:,4)
        !end if

        gr_center = gr_center/grc

        gl_Gran_center2(:, iter, now) = gr_center

        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)

        S = norm2(c)  ! S = S/2
        c = c/S
        S = S/2.0

        ! Можно один раз проверить, правильно ли ориентирована нормаль!\

        
        call Get_center_move(gl_Gran_neighbour(1, iter), node1, now)
        node2 = (p(:,1) + p(:,2) + p(:,3) + p(:,4))/4.0
        node1 = node2 - node1
        
		if (.False.) then
        !if(DOT_PRODUCT(node1, c) < 0.0) then
            print*, c
            print*, p(1, 1), norm2( (/0.0_8, p(2, 1), p(3, 1) /) )
            print*, "ERROR 13123 move ", DOT_PRODUCT(node1, c)
            print*, iter
            pause
        end if

        ! Нужно записать площадь грани и нормаль в общий массив!
        gl_Gran_normal2(:, iter, now) = c
        gl_Gran_square2(iter, now) = S

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
                c = c + (/gl_x2(i, now), gl_y2(i, now), gl_z2(i, now)/)
            end do
            c = c/8.0
        else
            if (gl_all_Cell(1, iter) == gl_all_Cell(2, iter)) then
                if (gl_all_Cell(4, iter) == gl_all_Cell(8, iter)) then  ! Класс ячейки = 4
                    do j = 1,8
                        ! if (j == 2 .or. j == 5 .or. j == 6 .or. j == 8) CYCLE
                        i = gl_all_Cell(j, iter)
                        c = c + (/gl_x2(i, now), gl_y2(i, now), gl_z2(i, now)/)
                    end do
                    c = c/8.0  ! 4
                else
                    do j = 1,8
                        ! if (j == 2 .or. j == 5 .or. j == 6) CYCLE
                        i = gl_all_Cell(j, iter)
                        c = c + (/gl_x2(i, now), gl_y2(i, now), gl_z2(i, now)/)
                    end do
                    c = c/8.0
                end if
            else
                if (gl_all_Cell(2, iter) == gl_all_Cell(6, iter)) then
                    if (gl_all_Cell(4, iter) == gl_all_Cell(5, iter)) then ! Класс ячейки = 1
                        do j = 1,8   ! Особая ячейка типа 4
                            ! if (j == 4 .or. j == 5 .or. j == 6 .or. j == 8) CYCLE
                            i = gl_all_Cell(j, iter)
                            c = c + (/gl_x2(i, now), gl_y2(i, now), gl_z2(i, now)/)
                        end do
                        c = c/8.0
                    else
                        do j = 1,8         ! Особая ячейка типа 2
                            ! if (j == 5 .or. j == 6) CYCLE             ! Попробовали убрать такой расчёт центра!
                            i = gl_all_Cell(j, iter)
                            c = c + (/gl_x2(i, now), gl_y2(i, now), gl_z2(i, now)/)
                        end do
                        c = c/8.0    ! 6 ----------------------------------------------------------------------------
                    end if
                else
                    if (gl_all_Cell(4, iter) == gl_all_Cell(5, iter)) then
                        do j = 1,8
                            ! if (j == 4 .or. j == 5 .or. j == 8) CYCLE
                            i = gl_all_Cell(j, iter)
                            c = c + (/gl_x2(i, now), gl_y2(i, now), gl_z2(i, now)/)
                        end do
                        c = c/8.0
                    else            ! Особая ячейка типа 5
                        do j = 1,8
                            ! if (j == 5 .or. j == 8) CYCLE    ! Попробовали убрать такой расчёт центра!
                            i = gl_all_Cell(j, iter)
                            c = c + (/gl_x2(i, now), gl_y2(i, now), gl_z2(i, now)/)
                        end do
                        c = c/8.0   ! 6 ----------------------------------------------------------------------------
                    end if
                end if
            end if
        end if

        gl_Cell_center2(:, iter, now) = c

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
            b = (/gl_x2(k, now), gl_y2(k, now), gl_z2(k, now)/)
            a = gl_Gran_normal2(:,i, now)
            di = dabs(DOT_PRODUCT(c,a) - DOT_PRODUCT(a,b))  ! Расстояние от точки до плоскости
            Vol = Vol + gl_Gran_square2(i, now) * di/3.0
            dist = MIN(dist, di)
            ll = ll + 1
        end do



        !if (iter == gl_Cell_A(5, 1, 1)) then
        !    write(*, *) c
        !    print*, Vol
        !    print*, gl_Cell_gran(:, iter)
        !    pause
        !end if


        gl_Cell_Volume2(iter, now) = Vol
        gl_Cell_dist(iter) = dist
		
		
		

        !        print*, Vol, ll
        !        pause

    end do


	end subroutine calc_all_Gran_move
	
	
	subroutine Surface_setup()
	! Эта функция передвигает гелиопаузу по линиям тока жидкости
	! Также заполняет подвинутые ячейки интерполированными значениями полей
	! Variables
	use GEO_PARAM
	use STORAGE
	use Interpolate
	use My_func
	implicit none
	
	real(8) :: position(10 * par_l_phi, 500)
	integer(4) :: posit_int(10 * par_l_phi, 500)
	real(8) :: Left, Right, x, y, z, r, phi, phi2, dsc
	integer(4) :: i, j, k, num, gr, ip, im
	integer(4) :: zone, cell, number, now, now2, yzel
	real(8) :: F(9), F_mf(5, 4), normal(3)
	
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
        
        ! Начальная инициализация
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
	
	number = 1
	cell = 1
	
	Left = -450.0_8
	Right = 50.0_8
	position = 0.0_8
	posit_int = 0
	
	do i = 1, 60
		print*, "i = ", i
		do k = 1, 10 * par_l_phi
			number = 1
			x = 25.0
			r = i * 0.2_8
			phi = 2.0 * par_pi_8 * k/(10 * par_l_phi) - par_pi_8/(10 * par_l_phi)
			y = r * cos(phi)
			z = r * sin(phi)
			
			do while (x > -450.0)
				!print*, "x = ", x
				call Interpolate_point(x, y, z, F, F_mf, zone, cell, number)
				if (cell == -1) EXIT
				number = cell
				x = x + 0.01 * F(2)
				y = y + 0.01 * F(3)
				z = z + 0.01 * F(4)
				
				phi2 = polar_angle(y, z)
				
				!print*, phi2, sqrt(y**2 + z**2)
				!pause
				
				position( INT(phi2/(2.0 * par_pi_8/(10 * par_l_phi))) + 1, INT(x + 450.01)) = position( INT(phi2/(2.0 * par_pi_8/(10 * par_l_phi))) + 1, INT(x + 450.01)) + sqrt(y**2 + z**2)
				posit_int(INT(phi2/(2.0 * par_pi_8/(10 * par_l_phi))) + 1, INT(x + 450.01))= posit_int(INT(phi2/(2.0 * par_pi_8/(10 * par_l_phi))) + 1, INT(x + 450.01)) + 1
			end do
			
		end do	
	end do
	
	
	do j = 1, 500
		do i = 1, 10 * par_l_phi
			if(posit_int(i, j) > 0) then
				position(i, j) = position(i, j)/posit_int(i, j)
			end if
			
		end do
	end do
	
	
	do k = 1, 20
		do j = 1, 500
			do i = 1, 10 * par_l_phi
				ip = i + 1
				im = i - 1
				if(i == 1) im = 10 * par_l_phi
				if(i == 10 * par_l_phi) ip = 1
			
				if(position(i, j) <= 0.001) then
					if(position(ip, j) > 0.001 .and. position(im, j) > 0.001) then
						position(i, j) = 0.5 * (position(ip, j) + position(im, j))
					else if (position(ip, j) < 0.001 .and. position(im, j) > 0.001) then
						position(i, j) = position(im, j)
					else if (position(im, j) < 0.001 .and. position(ip, j) > 0.001) then
						position(i, j) = position(ip, j)
					end if
				end if
			
			end do
		end do
	end do
	
	do j = 1, 500
			do i = 1, 10 * par_l_phi
				ip = i + 1
				im = i - 1
				if(i == 1) im = 10 * par_l_phi
				if(i == 10 * par_l_phi) ip = 1
			
				position(i, j) = (position(ip, j) + position(i, j) + position(im, j))/3.0
			
			end do
		end do
	
	
	
	open(2, file = 'save_SERF.txt')
	
	do j = 1, 500
		do i = 1, 10 * par_l_phi
			x = Left + j
			r = position(i, j)
			phi = i * (2.0 * par_pi_8/(10 * par_l_phi))
			y = r * cos(phi)
			z = r * sin(phi)
			write(2,*) x, y, z
		end do
	end do
	
	close(2)
	
	
	
	! Теперь передвигаем поверхность
	
	! Контакт
    Num = size(gl_Contact)
	now = 2
	
	do k = 1, 600
		if( mod(k, 10) == 0) print*, "k = ", k, "  from 600"
		now2 = now
		now = mod(now, 2) + 1
	
	
		gl_Vx = 0.0
		gl_Vy = 0.0
		gl_Vz = 0.0
		gl_Point_num = 0
    
		do i = 1, Num
			!print*, "i = ", i
			gr = gl_Contact(i)
		
			x = gl_Gran_center2(1, gr, now)
			y = gl_Gran_center2(2, gr, now)
			z = gl_Gran_center2(3, gr, now)
			phi2 = polar_angle(y, z)
		
			if (x > -50) CYCLE
		
			!print*, INT(phi2/(2.0 * par_pi_8/par_l_phi)) + 1, INT(x + 450.01)
			if (sqrt(y**2 + z**2) < position( INT(phi2/(2.0 * par_pi_8/(10 * par_l_phi))) + 1, INT(x + 450.01))) then
				dsc = 0.1
			else
				dsc = -0.1
			end if
		
		
		
			normal = gl_Gran_normal2(:, gr, now)
		
			do j = 1, 4
				yzel = gl_all_Gran(j, gr)
			
				gl_Vx(yzel) = gl_Vx(yzel) + normal(1) * dsc
				gl_Vy(yzel) = gl_Vy(yzel) + normal(2) * dsc
				gl_Vz(yzel) = gl_Vz(yzel) + normal(3) * dsc
				gl_Point_num(yzel) = gl_Point_num(yzel) + 1
			end do
		
		end do
		!print*, "Move_all "
		call Move_all(now, 1.0_8) 
		call calc_all_Gran_move(now2)
	
	
	end do
	
	gl_x = gl_x2(:, now2)
    gl_y = gl_y2(:, now2)
    gl_z = gl_z2(:, now2)
    gl_Cell_Volume = gl_Cell_Volume2(:, now2)
    gl_Gran_normal = gl_Gran_normal2(:, :, now2)
    gl_Gran_center = gl_Gran_center2(:, :, now2)
    gl_Cell_center = gl_Cell_center2(:, :, now2)
    gl_Gran_square = gl_Gran_square2(:, now2)
	
	
	! Body of Surface_setup
	
	end subroutine Surface_setup
	

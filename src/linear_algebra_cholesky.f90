program linear_algebra_cholesky
	use omp_lib
	implicit none
	! logical, parameter :: isCheck = .false., isShow = .false.
	integer, parameter :: n = 500, m = 10, cutoff = 0
	real(4), dimension(1:n, 1:n) :: a, l, lT, id
	real(4), dimension(1:n, 1:m) :: b, lb, x, ax
	real(4), dimension(1:n) :: y
	integer :: i, j, ip
	real(8) :: time, s
	! Initialization of the matrix A
	call random_number(a)
	ForAll(i = 1:n, j = 1:n) id(i,j) = (i/j)*(j/i)
	a = (a + transpose(a)) / 2 + n * id
	! Matrix b initialization
	b(:,:) = 1.0
	! Initialization of a time counter
	time = omp_get_wtime()
	! Cholesky decomposition
	l = a
	do  i = 1, n
		s = l(i,i)
		do  ip = 1, i - 1
			s = s - dprod(l(i,ip), l(i,ip))
		end do
		l(i,i) = sqrt(s)
		!$omp parallel do private(j,ip,s) shared(l,i) schedule(static) if (j < n - cutoff)
		do  j = i + 1, n
			s = l(j,i)
			do  ip = 1, i-1
				s = s - dprod(l(i,ip), l(j,ip))
			end do
			l(j,i) = s / l(i,i)
		end do
		!$omp end parallel do
	end do
	! back substitution
	lT = transpose(l)
	!$omp parallel do private(i,j,y,s) shared(x,l,lT) schedule(static)
	do i = 1, m
		y(1) = b(1,i) / l(1,1)
		do j = 2, n
			s = b(j,i)
			do ip = 1, j - 1
				s = s - dprod(l(j,ip), y(ip))
			end do
			y(j) = s / l(j,j)
		end do
		x(n,i) = y(n) / lT(n,n)
		do j = n - 1, 1, -1
			s = y(j)
			do ip = j + 1, n
				s = s - dprod(lT(j,ip), x(ip,i)) 
			end do
			x(j,i) = s / lT(j,j)
		end do
	end do
	!$omp end parallel do
	! Time
	time = omp_get_wtime() - time
	print *, 'Execution time: ', time
	! ! Show matrices
	! if (isShow) then
	! 	print *, 'A:'
	! 	do i=1, n
	! 		print *, a(i,:)
	! 	end do
	! 	print *, 'X:'
	! 	do i=1, n
	! 		print *, x(i,:)
	! 	end do
	! 	print *, 'AX:'
	! 	ax = matmul(a,x)
	! 	do i=1, n
	! 		print *, ax(i,:)
	! 	end do
	! 	print *, 'B:'
	! 	do i=1, n
	! 		print *, b(i,:)
	! 	end do
	! end if
	! ! Check
	! if (isCheck) then
	! 	print *, 'Checking...'
	! 	ax = matmul(a,x)
	! 	outer: do i = 1, n
	! 		do j = 1, m
	! 			if (abs(ax(i,j) - b(i,j)) .gt. 1e-6) then
	! 				print *, 'There is a mistake', i, j, ax(i,j), b(i,j)
	! 				exit outer
	! 			end if
	! 		end do 
	! 	end do outer
	! 	print *, 'Сheck has been passed successfully'
	! end if
end program linear_algebra_cholesky
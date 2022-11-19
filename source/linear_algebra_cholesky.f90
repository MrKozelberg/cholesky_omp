program linear_algebra_cholesky
	use omp_lib
	implicit none
	! logical, parameter :: isCheck = .false., isShow = .false.
	integer, parameter :: m = 10
	integer :: n, cutoff
	real(4), allocatable :: a(:,:), l(:,:), lT(:,:), id(:,:)
	real(4), allocatable :: b(:,:), x(:,:), y(:)
	integer :: i, j, ip
	real(8) :: time, s
	! Reading n and cutoff
	read(*, *) n, cutoff
	! Id matrix
	allocate(id(n,n))
	ForAll(i = 1:n, j = 1:n) id(i,j) = (i/j)*(j/i)
	! Initialization of the matrix A
	allocate(a(n,n))
	call random_number(a)	
	a = (a + transpose(a)) / 2 + n * id
	! Matrix b initialization
	allocate(b(n,m))
	b(:,:) = 1.0
	! Initialization of a time counter
	time = omp_get_wtime()
	! Cholesky decomposition
	allocate(l(n,n))
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
	allocate(lT(n,n))
	lT = transpose(l)
	allocate(y(n))
	allocate(x(n,m))
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
	write(*,*) n, time
	deallocate(id)
	deallocate(a)
	deallocate(b)
	deallocate(l)
	deallocate(lT)
	deallocate(y)
	deallocate(x)
end program linear_algebra_cholesky
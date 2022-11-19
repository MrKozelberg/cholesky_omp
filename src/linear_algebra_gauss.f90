program linear_algebra_gauss
	integer, parameter :: n = 3
	real, dimension(1:n, 1:n) :: a, left
	real, dimension(1:n) :: x, b
	data left/9 * 0./
	data a/ .471, 4.27, .012, 3.21, -.513, 1.273, -1.307, 1.102, -4.175 /
	data b/ 2.425, -.176, 1.423 /
	! решение: 0.07535443 0.6915624 -0.1297573
	! прямой ход метода гаусса
	do k = 1, n - 1
		do i = k + 1, n
			left(i, k) = a(i, k) / a(k, k)
			b(i) = b(i) - left(i, k) * b(k)
		end do
		do j = k + 1, n
			do i = k + 1, n
				a(j, i) = a(j, i) - left(j, k) * a(k, i)
			end do
		end do
	end do
	! обратная подстановка
	x(n) = b(n) / a(n, n)
	do i = n - 1, 1, -1
		x(i) = b(i)
		do j = i + 1, n
			x(i) = x(i) - a(i, j) * x(j)
		end do
		x(i) = x(i) / a(i, i)
	end do
	do i=1, n 
		print *, a(i,1), a(i,2), a(i, 3)
	end do
	print *, ''
	print *, (x(i), i = 1, n)
end program linear_algebra_gauss
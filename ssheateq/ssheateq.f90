program steadystateheateq
implicit none
real*8, dimension(:,:), allocatable :: T
real*8, dimension(:), allocatable :: b, x
integer, dimension(:), allocatable :: p, q
integer :: allocatestatus, k, l, m, i, j, meth, iterations, allow, alhigh, n = 1, start, finish
real*8 :: Tal = 100, Tsp = 2.725, Tsun = 340.31

n = 10*n+1

allow = int(0.3*n+0.7)
alhigh = int(0.7*n+0.3)

allocate(T(n**3,n**3), STAT = allocatestatus)
if (AllocateStatus /= 0) stop "***Insufficient memory***"
allocate(b(n**3), STAT = allocatestatus)
if (allocatestatus /= 0) stop "***Insufficient memory***"
allocate(x(n**3), STAT = allocatestatus)
if (allocatestatus /= 0) stop "***Insufficient memory***"
allocate(p(n**3), STAT = allocatestatus)
if (allocatestatus /= 0) stop "***Insufficient memory***"
allocate(q(n**3), STAT = allocatestatus)
if (allocatestatus /= 0) stop "***Insufficient memory***"

T = 0
b = 0

! INITIALIZING FACES

do k = 1,n
	do l = 1,n
	T((k-1)*n + l, (k-1)*n + l) = 1
	b((k-1)*n + l) = Tsun ! i = 1 FACE TO SUN
	T((n-1)*n*n + (k-1)*n + l, (n-1)*n*n + (k-1)*n + l) = 1
	b((n-1)*n*n + (k-1)*n + l) = Tsp ! i = n FACE FARTHEST FROM SUN
	T((k-1)*n*n + l, (k-1)*n*n + l)= 1
	b((k-1)*n*n + l) = (Tsp - Tsun)/(n-1)*(k-1) + Tsun ! GRADIENT AWAY FROM SUN, j = 1
	T((k-1)*n*n + (n-1)*n + l, (k-1)*n*n + (n-1)*n + l) = 1
	b((k-1)*n*n + (n-1)*n + l) = (Tsp - Tsun)/(n-1)*(k-1) + Tsun ! GRADIENT, j = n 
	T((k-1)*n*n + (l-1)*n + 1, (k-1)*n*n + (l-1)*n + 1) = 1 
	b((k-1)*n*n + (l-1)*n + 1) = (Tsp - Tsun)/(n-1)*(k-1) + Tsun ! GRADIENT, k = 1
	T((k-1)*n*n + (l-1)*n + n, (k-1)*n*n + (l-1)*n + n) = 1
	b((k-1)*n*n + (l-1)*n + n) = (Tsp - Tsun)/(n-1)*(k-1) + Tsun ! GRADIENT, k = n
	end do
end do

! FINITE DIFFERENCES WEIGHTING FOR ALL INTERIOR POINTS

do k = 2, n-1
	do l = 2, n-1
		do m = 2, n-1
		T((k-1)*n*n + (l-1)*n + m, (k-1)*n*n + (l-1)*n + m) = -6
		T((k-1)*n*n + (l-1)*n + m, (k-1 + 1)*n*n + (l-1)*n + m) = 1
		T((k-1)*n*n + (l-1)*n + m, (k-1 - 1)*n*n + (l-1)*n + m) = 1
		T((k-1)*n*n + (l-1)*n + m, (k-1)*n*n + (l-1 + 1)*n + m) = 1
		T((k-1)*n*n + (l-1)*n + m, (k-1)*n*n + (l-1 - 1)*n + m) = 1
		T((k-1)*n*n + (l-1)*n + m, (k-1)*n*n + (l-1)*n + m + 1) = 1
		T((k-1)*n*n + (l-1)*n + m, (k-1)*n*n + (l-1)*n + m - 1) = 1
		end do
	end do
end do

! ALUMINUM BOX CONDITIONS
! Removes finite differences placed where unnecessary, reassigns as definite solutions

do k = allow,alhigh
	do l = allow,alhigh
		do m = allow,alhigh
		T((k-1)*n*n + (l-1)*n + m, (k-1 + 1)*n*n + (l-1)*n + m) = 0
		T((k-1)*n*n + (l-1)*n + m, (k-1 - 1)*n*n + (l-1)*n + m) = 0
		T((k-1)*n*n + (l-1)*n + m, (k-1)*n*n + (l-1 + 1)*n + m) = 0
		T((k-1)*n*n + (l-1)*n + m, (k-1)*n*n + (l-1 - 1)*n + m) = 0
		T((k-1)*n*n + (l-1)*n + m, (k-1)*n*n + (l-1)*n + m + 1) = 0
		T((k-1)*n*n + (l-1)*n + m, (k-1)*n*n + (l-1)*n + m - 1) = 0
		T((k-1)*n*n + (l-1)*n + m, (k-1)*n*n + (l-1)*n + m) = 1
		b((k-1)*n*n + (l-1)*n + m) = Tal
		end do
	end do
end do

print*,'Array initialized...'

!open(unit=100, file = 'array.txt')
!open(unit=200, file = 'b.txt')
!do k = 1,n**3
!	write(100,*) (T(k,l), l = 1,n**3)
!	write(200,*) b(k)
!end do
!close(100)
!close(200)

print*,'Which method?'
print*,'Self-written (1), LAPACK (2)'
read(*,*) meth

if (meth == 1) then
call system_clock(count=start)
call computepaq(size(T,1),T,p,q)
print*,'PAQ factorization calculated...'
call paqsolve(size(T,1),T,b,p,q,x)
print*,'Solution found!'
call system_clock(count=finish)
print*,'Calculation time:', float(finish-start)/1000
open(unit=300, file = 'ans.txt')
print*,'Writing data...'

do l = 1,n**3
	k = mod(l,n)
	if(k==0) then
		k = n
	end if
	j = mod((l-k)/n + 1, n)
	if(j==0) then
		j = n
	end if
	i = (l-k)/(n**2) + (1-j)/n + 1
	if(i == 0) then
		i = n
	end if
	write(300,*) i, j ,k , x(l)
end do

else if (meth==2) then
call system_clock(count=start)
call dgesv(size(T,1),1,T,size(T,1),p,b,size(b,1),allocatestatus)
call system_clock(count=finish)
print*,'Calculation time:', float(finish-start)/1000
print*,'Solution found!'
open(unit=300, file = 'ans.txt')
print*,'Writing data...'
do l = 1,n**3
	k = mod(l,n)
	if(k==0) then
		k = n
	end if
	j = mod((l-k)/n + 1, n)
	if(j==0) then
		j = n
	end if
	i = (l-k)/(n**2) + (1-j)/n + 1
	if(i == 0) then
		i = n
	end if
	write(300,*) i, j ,k , b(l)
end do

end if

close(300)
deallocate(T,p,q,b,x)

end program steadystateheateq





subroutine computepaq(n,a,p,q)
double precision :: a(n,n), mult, maxvalue, temp(n)
integer :: n, i, j, k, l, maxrow, maxcol, p(n), q(n)

do i = 1,n
	p(i) = i
	q(i) = i	!initialize the row and column trackers
end do
do i = 1,n-1
	maxvalue = abs(a(i,i))
	maxrow = i
	maxcol = i
	do j = i,n
		do k = i,n
			if (abs(a(j,k)) > maxvalue) then
				maxrow = j
				maxcol = k
				maxvalue = abs(a(j,k))
			end if
		end do
	end do

	if (maxrow > i) then		! Reorganize a and record shifts in p and q
		p(i) = maxrow
		temp(i:n) = a(maxrow,i:n)
		a(maxrow,i:n) = a(i,i:n)
		a(i,i:n) = temp(i:n)
	end if
	if (maxcol > i) then
		q(i) = maxcol
		temp = a(1:n,maxcol)
		a(1:n,maxcol) = a(1:n,i)
		a(1:n,i) = temp(1:n)
	end if

	a(i+1:n,i) = a(i+1:n,i)/a(i,i)	!multipliers in the usual way
	
	do j = i+1,n
		a(j,i+1:n) = a(j,i+1:n) - a(j,i)*a(i,i+1:n)
	end do
end do

return
end

subroutine paqsolve(n,a,b,p,q,x)
implicit none
double precision :: a(n,n), b(n), w(n), x(n), y(n), s
integer :: n, i, j, k, p(n), q(n)

! First solve Ly = Pb

y = b

do k = 1,n-1
	if (p(k) > k) then
		s = y(k)
		y(k) = y(p(k))
		y(p(k)) = s
	end if
	do i = k+1,n
		y(i) = y(i) - a(i,k)*y(k)
	end do
end do

! Now solve Uy = Qx

x(n) = y(n)/a(n,n)

do i = n-1,1,-1
	s = y(i)
	do j = i+1,n
		s = s - a(i,j)*x(j)
	end do

	x(i) = s/a(i,i)
end do

! Next properly permute the rows again to solve x = Qy

do k = n-1,1,-1
	if(q(k) > k) then
		s = x(k)
		x(k) = x(q(k))
		x(q(k)) = s
	end if
end do

return
end

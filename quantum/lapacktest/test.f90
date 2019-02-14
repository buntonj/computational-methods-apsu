program lapacktest
real*8, dimension(:,:), allocatable :: A
integer :: numel = 5, info
real*8, allocatable, dimension(:) :: eigenvalues, work
integer :: lwork, i
character*1 :: jobz = 'V', uplo = 'U'

allocate(A(numel,numel), eigenvalues(numel),work(lwork))
lwork = 3*numel

do i = 1, numel
	A(i,i) = 2*i
end do

call dsyev(jobz,uplo,numel,A,numel,eigenvalues,work,lwork,info)

do i = 1,numel
	print*,(A(i,j), j = 1,numel)
end do 

end program lapacktest

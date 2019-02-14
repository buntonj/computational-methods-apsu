module constants

integer, parameter :: nmax = 5 ! number of basis function terms to include
integer, parameter :: elevels = 10
real*8, parameter :: pi = dacos(-1.0d0)
real*8, parameter :: sd = 4!0.1!1000!1D10	! well width
real*8, parameter :: depth = -2.179e-17!-3! -2.179D-5 ! well depth
real*8, parameter :: l = 1D-10	! period of basis functions and width of well
real*8, parameter :: hbar = 1.0545718D-34
real*8, parameter :: m = 9.1094D-31

!V(x) = depth*exp(-sd*(x/l)**2)

end module constants

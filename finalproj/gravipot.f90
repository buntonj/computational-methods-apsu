module gravipot
use constants
implicit none
contains

subroutine gforce(r, v, f, m, stage)
real*8, dimension(:) :: r, v, F
real*8 :: m, length, fmag
integer :: stage

length = norm2(r(:))
fmag = -GG*Me*m/(length**2)
F(:) = fmag*r(:)/length
if (stage == 0) then
	!F(:) = F(:) + v(:) * 2.0 * stage0fuelmass / tstage0
	!F(:) = F(:) + v(:)/norm2(v) * stage0thrust
	F(:) = F(:) + v(:)/norm2(v) * isp0mult* 2.0 * isp0 * g* stage0fuelmass / tstage0
else if (stage == 1) then
	F(:) = F(:) + v(:)/norm2(v) * isp1mult*isp1 * g * stage1fuelmass / tstage1
end if

end subroutine gforce


subroutine energy(r, v, E, m)
real*8, dimension(:) :: r, v
real*8 :: E, m

E = -GG*Me*m/norm2(r(:)) + 0.5*norm2(v(:))**2

end subroutine energy

end module gravipot

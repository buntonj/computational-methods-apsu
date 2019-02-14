module constants

real*8, parameter :: Re = 6.371D6 ! meters
real*8, parameter :: Me = 5.972D24 ! kg
real*8, parameter :: GG = 6.67408D-11 ! m^3 kg^-1 s^-2
real*8, parameter :: payloadmass = 1360.0D0 ! kg
real*8, parameter :: pi = dacos(-1.0D0)
real*8, parameter :: dt = 1D-2 ! seconds
real*8, parameter :: g = 9.8D0 ! m/s^2

! ATLAS LV-3B ROCKET PARAMETERS
real*8, parameter :: isp0 = 282 ! seconds
real*8, parameter :: isp1 = 309 ! seconds
real*8, parameter :: totalmass = 120000.0D0 ! kg
real*8, parameter :: stage1fuelmass = 110703.0D0 ! kg
real*8, parameter :: stage0fuelmass = 3350.0D0 ! kg
real*8, parameter :: stage0thrust = 1517.4D3 ! Newtons
real*8, parameter :: stage1mass = 116100.0D0 ! kg
real*8, parameter :: tstage0 = 134.0D0 ! seconds
real*8, parameter :: tstage1 = 303.0D0 ! seconds
real*8, parameter :: emptymass = totalmass - stage0fuelmass - stage1fuelmass
real*8, parameter :: isp0mult = 30
real*8, parameter :: isp1mult = 3
end module constants

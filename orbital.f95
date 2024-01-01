program orbital

implicit none

real::pi, t, tstep, ke, pe, te, mass
real, dimension(6)::state, k1, k2, k3, k4, temp
real, dimension(3)::r, v, a

pi = 4.0*atan(1.0)
t = 0.0
tstep = 0.00001
mass = 1.0

open(unit = 100, file = 'state')
open(unit = 110, file = 'r')
open(unit = 120, file = 'v')
open(unit = 130, file = 'ke')
open(unit = 140, file = 'pe')
open(unit = 150, file = 'te')

state = 0.0
state(1) = 1.0
state(5) = 2*pi
r = state(1:3)
v = state(4:6)
a = -4*pi**2/((norm2(r))**3) * r

write(100,*) state
write(110,*) r
write(120,*) v

do while (t<1.0)

	state(1:3) = r(1:3)
	state(4:6) = v(1:3)
	
	a = -4.0*pi**2.0/((norm2(r))**3) * r
	k1(1:3) = v(1:3)
	k1(4:6) = a(1:3)
	temp = state + 0.5 * tstep * k1
	r = temp(1:3)
	v = temp(4:6)
	
	a = -4.0*pi**2.0/((norm2(r))**3) * r
	k2(1:3) = v(1:3)
	k2(4:6) = a(1:3)
	temp = state + 0.5 * tstep * k2
	r = temp(1:3)
	v = temp(4:6)
	
	a = -4.0*pi**2.0/((norm2(r))**3) * r
	k3(1:3) = v(1:3)
	k3(4:6) = a(1:3)
	temp = state + tstep * k3
	r = temp(1:3)
	v = temp(4:6)
	
	a = -4.0*pi**2.0/((norm2(r))**3) * r
	k4(1:3) = v(1:3)
	k4(4:6) = a(1:3)
	
	state = state + tstep/6.0 * (k1 + 2.0*k2 + 2.0*k3 + k4)
	r = state(1:3)
	v = state(4:6) 
	t = t + tstep
	
	ke = 0.5 * mass * norm2(v)**2
	pe = -4*pi*mass/norm2(r)
	te = ke + pe
	
	write(100,*) state
	write(110,*) r
	write(120,*) v
	write(130,*) t, ke
	write(140,*) t, pe
	write(150,*) t, te
end do

end program orbital
! Code works, not sure what issues you had, but I ran it and got the correct graphs. Nice job.

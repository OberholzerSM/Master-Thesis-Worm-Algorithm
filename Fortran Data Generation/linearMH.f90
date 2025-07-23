program linearMH

	implicit none
	
	integer :: i, datasize, samplesize, u
	real (kind = 16) :: om, ob, ln_f
	real (kind = 16), dimension(201) :: x, y
	real (kind = 16), dimension(10**5,2) :: samples
	
	call random_seed()
	
	datasize = 201
	samplesize = 10**5
	
	call DatenGenerator(datasize, x, y)
	
	call MH_Sampler(datasize, samplesize, x, y)
	

end program linearMH


subroutine GaussSample(mu, sigma, x) !Generiert ein Sample aus einer Gauss-Verteilung. Inputs: mu, sigma. Output: x

	implicit none

	real (kind = 16) :: u1,u2,z, mu, sigma, x
	real (kind = 16), parameter :: pi = 4.0*atan(1.0)
	
	call random_number(u1)
	call random_number(u2)

	z = sqrt( -2.0*log(u1) )*cos(2.0*pi*u2)
	x = sigma*z + mu
	
end subroutine GaussSample

subroutine DatenGenerator(n, x, y) !Generiert n Daten y, die um x verteilt sind.

	implicit none
	
	integer :: n, i
	real (kind = 16) :: sigma
	real (kind = 16), dimension(n) :: x, y
	
	sigma = 0.5
	
	do i = 1, n
	
		x(i) = i - ((n-1)/2) - 1
		call GaussSample(2.0*x(i) - 1.0, sigma, y(i))
	
	end do

end subroutine DatenGenerator

real (kind = 16) function ln_f(m, b, n, x, y) !Berechnet den ln der Wahrscheinlichkeit, dass die Parameter m und b bei gegeben Daten x und y der Länge n stimmen.

	integer :: i,n
	real (kind = 16) :: ln_L, m, b, two, nreal
	real (kind = 16), parameter :: pi = 4*atan(1.0)
	real (kind = 16), dimension(201) :: x,y
	
	two = real(2.0,16)
	nreal = real(n,16)
	
	ln_L = 0
	do i = 1,n
	
		ln_L = ln_L + ( ( m*x(i) + b - y(i) )**2 )
		
	end do
	
	ln_L = ( -two * ln_L) - ( (nreal/two) * log(two / pi) )
	
	ln_f = ln_L - real(log(4.0),16)

end function ln_f

subroutine MH_Sampler(datasize, samplesize, x, y) !Erstellt zwei .txt Files mit den m und b samples

	implicit none
	
	integer :: i, datasize, samplesize, rate1, rate2, test
	real :: realrate1, realrate2 
	real (kind = 16) :: ln_f, om, ob, r1, r2
	real (kind = 16), dimension(201) :: x, y
	real (kind = 16), dimension(2) :: theta1, theta2  
	
	om = real(0.001,16)
	ob = real(0.006,16)
	theta1(1) = real(0,16)
	theta1(2) = real(0,16)
	theta2(1) = real(0,16)
	theta2(2) = real(0,16)
	rate1 = 0
	rate2 = 0
	test = 0
	
	open(1, file = 'm_samplings.txt', status = 'old')
	open(2, file = 'b_samplings.txt', status = 'old')
	
	do i = 1,samplesize
		
		call random_number(r1)
		call random_number(r2)
			
			
		call GaussSample(theta1(1), om, theta2(1))
			
			
		if ( ln_f(theta2(1), theta2(2), datasize, x, y) - ln_f(theta1(1), theta1(2), datasize, x, y)  > log(r1) ) then
			
			theta1(1) = theta2(1)
			rate1 = rate1 + 1
			
		end if
		
		write(1,*) theta1(1)
		
		call GaussSample(theta1(2), ob, theta2(2))
			
		if ( ln_f(theta2(1), theta2(2), datasize, x, y) - ln_f(theta1(1), theta1(2), datasize, x, y)  > log(r2) ) then
			
			theta1(2) = theta2(2)
			rate2 = rate2 + 1
			
		end if
			
		write(2,*) theta1(2)
		
		if (test == 0 .and. mod(i,1000) == 0 ) then
			
			if (rate1 < real(600,16)) then
				om = om - real(0.00001,16)
			end if
		
			if (rate1 > real(800,16)) then
				om = om + real(0.00001,16)
			end if
		
			if (rate2 < real(600,16)) then
				ob = ob - real(0.00001,16)
			end if
		
			if (rate2 > real(800,16)) then
				ob = ob + real(0.00001,16)
			end if

			if( min(rate1,rate2) > 600 .and. max(rate1,rate2) < 800 ) then
				
				test = 1
				print *, 'Daten koennen ab dem ', i, '. Schritt genutzt werden.'
				print *, 'om :', om
				print *, 'ob :', ob
			
			end if
			
			rate1 = 0
			rate2 = 0
			
		end if
		
	end do
	
	close(1)
	close(2)
	
	realrate1 = real(rate1)
	realrate1 = realrate1 / real(samplesize)
	print *, 'Rate m: ', realrate1
	
	realrate2 = real(rate2)
	realrate2 = realrate2 / real(samplesize)
	print *, 'Rate b: ', realrate2

end subroutine MH_Sampler
program linearMH

	implicit none
	
	integer :: datasize,samplesize,t,i,j,k,l_Block,n_Block
	real :: chi0,chi,tau
	real(kind=8) :: summe1,summe2,summe3,sd,sd_naiv
	real, dimension(201) :: x,y
	real, dimension(:),allocatable :: m_list,b_list,Block_list
	
	call random_seed()
	
	datasize = 201
	samplesize = 10**5
	
	allocate(m_list(samplesize))
	allocate(b_list(samplesize))
	
	call DatenGenerator(datasize,x,y)
	call MCMC_sampler(samplesize,datasize,x,y,m_list,b_list)
	
	!Autokorrelation für b, da für m Fortran nicht funktioniert
	chi0 = chi(0,samplesize,b_list)
	tau = ( ( chi0 + chi(1,samplesize,b_list) ) / (2.0*chi0) )
	t = 1
	do while( 6.0*tau > real(t)  )
		tau = tau + ( ( chi(t,samplesize,b_list) + chi(t+1,samplesize,b_list) ) / (2.0*chi0) )
		t = t + 1
	end do
	
	!Naiver Fehler
	print*,"Berechne naiven Fehler..."
	call bootstrap(samplesize,b_list,sd)
	sd_naiv = tau*sd
	print*,"Naiver Fehler",sd_naiv
	
	!Blockfehler
	print*,"Berechne Blockfehler..."
	l_Block = ceiling(6.0*tau) + 10 - modulo(ceiling(6.0*tau),10)!Länge eines Blocks: 6*tau aufgerundet auf die nächste durch 10 teilbare Zahl.
	
	do k = l_Block-20,l_Block
	
		n_Block = samplesize / k !Anzahl Blöcke
		summe1 = 0.0_8
		allocate(Block_list(k))
		
		do i = 1,n_Block
		
			!Erstelle die Liste innerhalb des Blocks
			do j = 1,k
				
				Block_list(j) = b_list(j+(i-1)*k)
				
			end do
		
			!Bestimme den Fehler und summiere ihn auf
			call bootstrap(k,Block_list,sd)
			summe1 = summe1 + sd**2
		end do
		
		print*,"Blockfehler",k,"/",l_Block,",",sqrt(summe1)/real(n_Block),"/",sd_naiv
		
		deallocate(Block_list)
	
	end do
	
	deallocate(m_list)
	deallocate(b_list)

end program linearMH

subroutine GaussSample(mu,sigma,x) !Generiert ein Sample aus einer Gauss-Verteilung. Inputs: mu, sigma. Output: x

	implicit none

	real :: u1,u2,z,mu,sigma,x
	real, parameter :: pi = 4.0*atan(1.0)
	
	call random_number(u1)
	call random_number(u2)

	z = sqrt( -2.0*log(u1) )*cos(2.0*pi*u2)
	x = sigma*z + mu
	
end subroutine GaussSample

subroutine DatenGenerator(n,x,y) !Generiert n Daten y, die um x verteilt sind.

	implicit none
	
	integer :: n,i
	real:: sigma
	real, dimension(n) :: x,y
	
	sigma = 0.5
	do i = 1, n
	
		x(i) = i - ((n-1)/2) - 1
		call GaussSample(2.0*x(i) - 1.0, sigma, y(i))
	
	end do

end subroutine DatenGenerator

function ln_f(m,b,datasize,x,y) !Berechnet den ln der Wahrscheinlichkeit, dass die Parameter m und b bei gegeben Daten x und y der Länge n stimmen.

	integer :: i,datasize
	real :: ln_L,ln_f,m,b,nreal
	real, parameter :: pi = 4*atan(1.0)
	real, dimension(201) :: x,y
	
	ln_L = 0.0
	do i = 1,datasize
		ln_L = ln_L + ( ( m*x(i) + b - y(i) )**2 )
	end do
	
	ln_L = -2.0*ln_L - (real(datasize)/2.0) * log(2.0/pi)
	
	ln_f = ln_L - log(4.0)

end function ln_f

subroutine MCMC_sampler(samplesize,datasize,x,y,m_list,b_list) !Erzeugt die Liste der samples

	implicit none
	integer :: samplesize,datasize,i
	real :: ln_f,u,proposal_m,proposal_b,sigma_m,sigma_b,p_A
	real, dimension(datasize) :: x,y
	real, dimension(samplesize) :: m_list,b_list

	!Startwerte
	m_list(1) = 2.0
	b_list(1) = -1.0
	
	sigma_m = 0.001
	sigma_b = 0.006
	
	open(1,file="m_samplings.dat")
	open(2,file="b_samplings.dat")
	
	write(1,*) m_list(1)
	write(2,*) b_list(1)
	
	do i = 2,samplesize
	
		!Erzeuge einen neuen Vorschlag
		call GaussSample(m_list(i-1),sigma_m,proposal_m)
		call GaussSample(b_list(i-1),sigma_b,proposal_b)
		
		!Bestimme die Akzepanzwahrscheinlichkeit
		p_A = exp(ln_f(proposal_m,proposal_b,datasize,x,y) - ln_f(m_list(i-1),b_list(i-1),datasize,x,y))
		call random_number(u)
		if (p_A > u) then
			m_list(i) = proposal_m
			b_list(i) = proposal_b
		else
			m_list(i) = m_list(i-1)
			b_list(i) = b_list(i-1)
		end if
		
		write(1,*) m_list(i)
		write(2,*) b_list(i)
		
	end do

	close(1)
	close(2)

end subroutine MCMC_sampler

function chi(t,samplesize,list) !Autokorrealtionsfunktion chi(t)

	implicit none
	integer :: t,samplesize,i
	real, dimension(samplesize) :: list
	real :: chi
	real:: summe1,summe2,summe3
	
	summe1 = 0.0
	summe2 = 0.0
	summe3 = 0.0
	
	do i = 1,samplesize-t
		summe1 = summe1 + list(i)*list(i+t)
		summe2 = summe2 + list(i)
		summe3 = summe3 + list(i+t)
	end do
	
	summe1 = summe1/real(samplesize-t)
	summe2 = summe2/real(samplesize-t)
	summe3 = summe3/real(samplesize-t)
	chi = summe1 - summe2*summe3

end function chi

function randomint(a,b) !Zufälliges Integer zwischen a und b (inklusive a und b)

	implicit none
	
	integer :: randomint,a,b
	real :: u
	
	call random_number(u)
	randomint = a + floor((b-a+1)*u)
	
end function randomint

subroutine bootstrap(samplesize,list,sd)!Fehlerberechnung. Inputs: samplesize,list. Output: sd

	implicit none
	integer :: samplesize,n_tries,i,j,t,randomint
	real(kind=8) :: sd,summe,summe2,average
	real, dimension(samplesize) :: list

	n_tries = 10**3
	summe = 0.0_8
	summe2 = 0.0_8
	
	do i = 1,n_tries
		
		!Bestimme den Durchschnittswert einer Liste mit zufälligen Elementen
		average = 0.0_8
		do j = 1,samplesize
			t = randomint(1,samplesize)
			average = average + real(list(t),8)
		end do
		average = average / real(samplesize,8)
		
		summe = summe + average
		summe2 = summe2 + average**2
		
	end do

	sd = sqrt( summe2/real(n_tries,8) - (summe/real(n_tries,8))**2 )

end subroutine bootstrap
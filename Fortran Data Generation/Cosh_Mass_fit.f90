program Cosh_Mass_fit
	
	implicit none
	integer :: N0,N1,n_tries,n_masse,n_Gitter,i,j,k,model,Gitter_modus 
	real :: M_bar,m,m_sd,C,C_sd
	
	n_tries = 10**5 !Anzahl Durchgänge pro Fit
	n_masse = 5 !Anzahl gemessene bare Massen
	n_Gitter = 6 !Anzahl gemessene Gitter
	Gitter_modus = 1 !Ob N0 oder N1 variiert wird
	model = 7 !7 Vertex oder 163 Vertex Model
	
	open(10,file="Mass_Fit.dat")
	write(10,*) "#m, m_sd, C, C_sd, M_bar, N0, N1"
	
	if (Gitter_modus < 3) then
		do i = 1,n_Gitter
			select case(Gitter_modus)
			case(1)
				N0 = 2**i
				N1 = 2
			case(2)
				N0 = 64
				N1 = 2**i
			end select
			do j = 1,n_masse
				M_bar = real(j)*2.0/real(n_masse)
				call Mass_fit(N0,N1,M_bar,n_tries,m,m_sd,C,C_sd,model)
				write(10,*) m,m_sd,C,C_sd,M_bar,N0,N1
				print*,100*(j + (i-1)*n_masse)/(n_masse*n_Gitter),"%"
			end do
		end do
	else !Variiere zuerst N0 und dann N1
		do i = 1,n_Gitter
			N0 = 2**i
			N1 = 2
			do j = 1,n_masse
				M_bar = real(j)*2.0/real(n_masse)
				call Mass_fit(N0,N1,M_bar,n_tries,m,m_sd,C,C_sd,model)
				write(10,*) m,m_sd,C,C_sd,M_bar,N0,N1
				print*,50*(j + (i-1)*n_masse)/(n_masse*n_Gitter),"%"
			end do
		end do
		do i = 1,n_Gitter
			N0 = 64
			N1 = 2**i
			do j = 1,n_masse
				M_bar = real(j)*2.0/real(n_masse)
				call Mass_fit(N0,N1,M_bar,n_tries,m,m_sd,C,C_sd,model)
				write(10,*) m,m_sd,C,C_sd,M_bar,N0,N1
				print*,50+50*(j + (i-1)*n_masse)/(n_masse*n_Gitter),"%"
			end do
		end do
	end if
	
	close(10)
	
end program Cosh_Mass_fit

subroutine Mass_fit(N0,N1,M_bar,n_tries,m,sigma_m,C,sigma_C,model) !Fittet die Masse für eine fixe bare Masse und Gittergrösse

	implicit none
	integer, parameter :: qp = selected_real_kind(33, 4931)
	integer :: i,j,k,t,N0,N1,n_tries,model
	real :: m,m_old,C,C_old,M_bar,sigma_m,sigma_C,u,pA
	real(kind=qp) :: c0,mu
	real :: random_Gauss,random_Cauchy,cosh_m
	real(kind=qp) :: log_weight,log_weight_old,log_Gauss
	integer(kind=qp), dimension(2*N0-1) :: dt_list
	real, dimension(2*N0-1) :: sd_list
	real(kind=qp), dimension(int(N0/2)+1) :: y,sigma_y
	real, dimension(n_tries) :: m_list,C_list
	character(len=10) :: file_M,file_N0,file_N1
	character(len=100) :: file_name
	
	!Lese die Daten ein
	write(file_M,"(f10.2)") M_bar
	write(file_N0,"(i3)") N0
	write(file_N1,"(i3)") N1
	
	if (model==7) then
		file_name = "7Vertex_dt_number_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	else if (model==163) then
		file_name = "163Vertex_dt_number_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	end if
	open(1,file=file_name)
	
	if (model==7) then
		file_name = "7Vertex_dt_sd_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	else if (model==163) then
		file_name = "163Vertex_dt_sd_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	end if
	open(2,file=file_name)
	
	do i = 1,2*N0-1
		read(1,*) dt_list(i)
		read(2,*) sd_list(i)
	end do
	
	close(1)
	close(2)
	
	!Fasse die äquivalenten Daten zusammen	
	c0 = real(dt_list(N0),qp)
	y(1) = 1.0_qp
	sigma_y(1) = real(sd_list(N0),qp)/c0
	
	do i = 1,int(N0/2)
		y(i+1) = real( dt_list(i) + dt_list(N0-i) + dt_list(N0+i) + dt_list(2*N0-i),qp ) / (2.0_qp*c0)
		sigma_y(i+1) = sqrt( real(sd_list(i)**2 + sd_list(N0-i)**2 + sd_list(N0+i)**2 + sd_list(2*N0-i)**2,16) ) / (2.0_qp*c0)
	end do
	
	!Startwere
	m_old = 0.0
	m = m_old
	m_list(1) = m_old
	C_old = 1.0
	C = C_old
	C_list(1) = C_old
	!Bestimme das Gewicht
	log_weight_old = 0.0_qp
	do t = 1,int(N0/2)+1
		if (y(t) > epsilon(0.0)) then
			mu = real(cosh_m(t-1,m,C,N0),qp)
			log_weight_old = log_weight_old + log_Gauss(y(t),mu,sigma_y(t))
		end if
	end do

	!MCMC Algorithmus
	do k = 2,n_tries
		
		!Schlage eine zufällige neue Masse vor
		sigma_m = random_Cauchy(0.0,1.0)
		m = abs(random_Gauss(m_old,sigma_m))
		sigma_C = random_Cauchy(0.0,1.0)
		C = random_Gauss(C_old,sigma_m)
		
		!Bestimme das Gewicht
		log_weight = 0.0_qp
		do t = 1,int(N0/2)+1
			if (y(t) > epsilon(0.0)) then !Nur falls eine Messung vorhanden ist
				mu = real(cosh_m(t-1,m,C,N0),qp)
				log_weight = log_weight + log_Gauss(y(t),mu,sigma_y(t))
			end if
		end do
		
		!Teste, ob das neue Gewicht akzeptiert wird
		call random_number(u)
		if ( log_weight - log_weight_old > real(log(u),qp) ) then
			m_old = m
			C_old = C
			log_weight_old = log_weight
		end if
		
		m_list(k) = m_old
		C_list(k) = C_old
		
	end do
	
	!Berechne m und C
	m = sum(m_list)/real(n_tries)
	C = sum(C_list)/real(n_tries)
	
	!Berechne den Fehler
	call Blockfehler(m_list,n_tries,sigma_m)
	call Blockfehler(C_list,n_tries,sigma_C)

end subroutine Mass_fit

subroutine Blockfehler(liste,n,sd) !Nimmt eine Liste von reals der Länge n und berechnet den Fehler des Durchschnitts der Liste.

	implicit none
	integer :: n,n_Block,l_Block,i,j,k,t,randomint,n_tries
	real(kind=16) :: summe0,summe1,summe2
	real,dimension(n) :: liste
	real,dimension(100) :: Block_list !Länge ist "n_Block".
	real :: sd
	
	n_Block=100
	l_Block=n/n_Block
	
	!Block_list: m für jeden Block
	do i=1,n_Block
		Block_list(i)=0.0
		do j=1,l_Block
			Block_list(i)=Block_list(i)+liste(j+(i-1)*l_Block)
		end do
		Block_list(i) = Block_list(i)/l_Block
	end do
	
	!Bootstrappe die Blockliste, um den Fehler abzuschätzen
	n_tries=10**3
	summe1=0_16
	summe2=0_16
	do k=1,n_tries
		summe0 = 0_16
		do i=1,n_Block
			t=randomint(1,n_Block)
			summe0 = summe0 + Block_list(t)
		end do
		summe1 = summe1 + summe0
		summe2 = summe2 + summe0**2
	end do
	
	sd = sqrt( real(summe2,16)/real(n_tries,16) - (real(summe1,16)/real(n_tries,16))**2 )

end subroutine Blockfehler

function cosh_m(t,m,C,N0) !Pseudo cosh(t)

	implicit none
	integer :: t,N0
	real :: cosh_m,x,m,A,B,C
	
	x = real(t)
	A = C / ( 1.0 + exp(-m*real(N0)) )
	B = C / ( 1.0 + exp(m*real(N0)) )
	
	cosh_m = A*exp(-m*x) + B*exp(m*x)

end function cosh_m

function inverse_cosh(t,m,C,N0) !Inverse Funktion von A exp[-mt] + B exp[+mt]

	implicit none
	integer :: t,N0
	real :: inverse_cosh,x,m,A,B,C,d,e
	
	x = real(t)
	A = C / ( 1.0 + exp(-m*real(N0)) )
	B = C / ( 1.0 + exp(m*real(N0)) )
	e = x / (2.0*B)
	d = x**2/(4.0*B**2) - A/B
	
	if (d > 0.0) then
		inverse_cosh = log(e+sqrt(d))
	else
		inverse_cosh = -1.0
	end if

end function inverse_cosh

function log_Gauss(x,mu,sigma) !Gibt den ln der Normalverteilungsdichte zurück

	implicit none
	integer, parameter :: qp = selected_real_kind(33, 4931)
	real(kind=qp) :: x,mu,sigma,pi
	real(kind=qp) :: log_Gauss
	
	pi=real(2.0*acos(0.0),qp)
	log_Gauss = real(log(1.0_qp/(sqrt(2.0_qp*pi)*sigma)) - (x-mu)**2/(2.0_qp*sigma**2),qp)

end function log_Gauss

function random_Gauss(mu,sigma) !Erzeugt eine normalverteilte Zufallsvariable

	implicit none
	integer :: test
	real :: u1,u2,p,q,mu,sigma,random_Gauss
	
	test = 0
	do while(test==0)
	
		call random_number(u1)
		u1 = 2.0*u1 - 1.0
		call random_number(u2)
		u2 = 2.0*u2 - 1.0
		
		q = u1**2 + u2**2
		if (q > 0.0 .and. q < 1.0) then
			test = 1
		end if
	end do
	
	p = sqrt( -2.0*log(q)/q )
	random_Gauss = sigma*p*u1 + mu

end function random_Gauss

function random_Cauchy(t,s) !Gibt eine zufällige Cauchy Zufallsvariable zurück

	implicit none
	real :: random_Cauchy,t,s,pi,cauchy,u
	
	call random_number(u)
	pi=2.0*acos(0.0)
	cauchy = 1.0 / tan(pi*u)
	random_Cauchy = s*cauchy + t
	
end function random_Cauchy

function randomint(a,b) !Zufälliges Integer zwischen a und b (inklusive a und b)

	implicit none
	integer :: randomint,a,b
	real :: u
	
	call random_number(u)
	randomint = a + floor(real(b-a+1)*u)
	
end function randomint
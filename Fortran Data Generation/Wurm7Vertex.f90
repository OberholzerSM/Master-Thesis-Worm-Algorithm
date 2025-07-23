program Wurm7Vertex

	implicit none
	integer :: N0,N1,n_data,i,j,k0,k1,n_masse,option,n_Block
	real :: M,M_max,M_min,start_time,end_time
	character(len=10) :: file_k,file_M,file_N0,file_N1
	character(len=100) :: file_name

	call cpu_time(start_time)
	
	option=1
	if (option<3 .and. option/=0) then !Lese die Gittergrösse ein
		call getarg(1,file_k)
		read(file_k,*) k0
		call getarg(2,file_k)
		read(file_k,*) k1
	end if
	if (option==1) then !Lese die Masse ein
		call getarg(3,file_M)
		read(file_M,*) M
	end if
	N0 = 2**k0	!Gittergrösse
	N1 = 2**k1
	n_data = 10**5 !Anzahl Daten pro Massenpunkt
	n_Block = 100 !In wie viele Blöcke die Daten zusammengefasst werden
	n_masse=20 !Anzahl Massenpunkte zwischen M=0.0 und M=2.0
	
	select case(option)
	case(0) !Erstelle .txt File mit Parametern
		call make_txt_file(n_masse)
	case(1) !Fixe Gittergrösse und Masse
		call Auswertung(N0,N1,M,n_data,n_Block)
	case(2) !Fixe Gittergrösse, variiere Masse
		do i = 1,n_masse
			M = i*2.0/real(n_masse)
			call Auswertung(N0,N1,M,n_data,n_Block)
		end do
	case(3) !Variiere N0 und Masse
		do j = 1,6
			N0 = 2**j
			N1 = 2 
			do i = 1,n_masse
				M = i*2.0/real(n_masse)
				call Auswertung(N0,N1,M,n_data,n_Block)
			end do
		end do
	case(4) !Variiere N1 und Masse
		do j = 1,6
			N0 = 64
			N1 = 2 **j
			do i = 1,n_masse
				M = i*2.0/real(n_masse)
				call Auswertung(N0,N1,M,n_data,n_Block)
			end do
		end do
	end select

	call cpu_time(end_time)
	print*,"Time: ",end_time-start_time

end program Wurm7Vertex

subroutine make_txt_file(n_masse)

	implicit none
	integer :: n_masse,i,j,k0,k1
	real :: M,M_max,M_min
	
	M_min = 0.0
	M_max = 2.0
	open(1,file="params_51.txt")
	write(1,*) 0
	!do k0 = 1,6
		k0=5
		!do k1 = 1,5
			k1=6
			do i = 1,n_masse
				M = real(i)*(M_max-M_min)/real(n_masse) + M_min
				write(1,*) k0,k1,M
			end do
		!end do
	!end do
	close(1)

end subroutine make_txt_file

subroutine Auswertung(N0,N1,M,n_data,n_Block) !Wurm-Algorithmus für eine fixe Masse und Gittergrösse. Beinhaltet sowohl den Wurm als auch dessen Auswertung.

	implicit none
	integer :: N0,N1,n_data,topology,n_Block,l_Block
	integer :: i,j,k,t,Block_index
	real :: M
	real(kind=16) :: c0
	integer(kind=8), dimension(4) :: topology_number
	integer(kind=8), dimension(n_Block,4) :: topology_list
	real, dimension(4,4) :: topology_sd
	integer(kind=8), dimension(2*N0-1) :: dt_number,dt_list_Wurm
	integer(kind=8), dimension(n_Block,2*N0-1) :: dt_list
	real, dimension(2*N0-1,2*N0-1) :: dt_sd
	integer,dimension(N0,N1,4) :: lattice
	character(len=10) :: file_M,file_N0,file_N1
	character(len=100) :: file_name
	
	write(file_M,"(f10.2)") M
	write(file_N0,"(i3)") N0
	write(file_N1,"(i3)") N1
	
	l_Block=n_data/n_Block
	
	file_name = "7Vertex_dt_Blocklist_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	open(1,file=file_name)
	file_name = "7Vertex_topology_Blocklist_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	open(2,file=file_name)
	file_name = "7Vertex_last_lattice_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	open(3,file=file_name)
	
	!Startkonfiguration
	do i = 1,N0
		do j = 1,N1
			do k = 1,4
				lattice(i,j,k) = 0
			end do
			write(3,'(999999999I3)') lattice(i,j,:)
		end do
	end do
	
	!Topologie der Startkonfiguration
	do i = 1,4
		topology_number(i) = 0_8 !Wie häufig jede Topologieklasse vorkommt
		topology_list(1,i)=0_8 !Binäre Liste für die Fehlerbestimmung. Gibt an, welche Topologieklasse an der k-ten Stelle vorkommt.
	end do
	t = topology(lattice,N0,N1)
	topology_number(t) = 1_8
	topology_list(1,t)=1_8
	
	!Propagator Zeitlänge
	do i = 1,2*N0-1
		dt_number(i)=0_8
		dt_list(1,i)=0_8
	end do
	
	!Hauptloop: Erzeuge n_data viele "wurmlose" Konfigurationen
	mainloop: do k=2,n_data
		
		!Drucke den Fortschritt aus
		do i = 1,100
			if ( i*n_data/100==k ) then
				print*,i,"%, M=",file_M,", N=",file_N0,"x",file_N1
			end if
		end do
		
		!Lasse den Wurm laufen
		call Wurm(N0,N1,M,lattice,dt_list_Wurm)
		
		!Block-Index
		Block_index = 1+(k-1)/l_Block
		
		!Erfasse die Topologie-Statistik der wurmlosen Konfigurationen
		t = topology(lattice,N0,N1)
		topology_number(t) = topology_number(t)+1_8
		topology_list(Block_index,t)=topology_list(Block_index,t)+1_8
		
		!Erfasse die dt Statistik des vorherigen Wurm
		dt_list(Block_index,:) = dt_list(Block_index,:) + dt_list_Wurm
		dt_number = dt_number + dt_list_Wurm
		
		!Drucke den Blockwert aus, sobald er fertig ist
		if ( Block_index*l_Block==k ) then
		
			write(1,'(999999999I16)') dt_list(Block_index,:)
			write(2,'(999999999I9)') topology_list(Block_index,:)
			
			!Setze die Zähler für den nächsten Block auf 0
			if (Block_index /= n_block) then
				do i=1,4
					topology_list(Block_index+1,i)=0_8
				end do
				do i = 1,2*N0-1
					dt_list(Block_index+1,i)=0_8
				end do
			end if
			
			!Drucke die letzte generierte Lattice aus
			do i = 1,N0
				do j = 1,N1
					write(3,'(999999999I3)') lattice(i,j,:)
				end do
			end do
			rewind(3)
			
		end if
		
	end do mainloop
	close(1)
	close(2)
	close(3)
	
	!Drucke aus, wie häufig die Topologischen Klassen vorgekommen sind
	file_name = "7Vertex_topology_number_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	open(1,file=file_name)
	do i = 1,4
		write(1,'(I9)') topology_number(i)
	end do
	close(1)

	!Bestimme den Fehler der Topologie-Klassen
	print*,"Berechne Fehler Topologie..."
	file_name = "7Vertex_topology_sd_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	open(1,file=file_name)
	file_name = "7Vertex_topology_cov_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	open(2,file=file_name)
	call Bootstrap(real(topology_list,16)/real(n_data,16),n_Block,4,10**3,topology_sd)
	do i = 1,4
		write(1,'(999999999F0.9)') sqrt(topology_sd(i,i))
		write(2,'(999999999F30.9)') topology_sd(i,:)
	end do
	close(1)
	close(2)
	
	!Drucke aus, wie häufig die dt Abstände vorgekommen sind
	file_name = "7Vertex_dt_number_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	open(1,file=file_name)
	do i = 1,2*N0-1
		write(1,'(I16)') dt_number(i)
	end do
	close(1)
	
	!Normierung dt
	c0 = real(sum(dt_list(:,N0)),16)
	!Bestimme den Fehler der dt Abstände
	print*,"Berechne Fehler dt..."
	file_name = "7Vertex_dt_sd_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	open(1,file=file_name)
	file_name = "7Vertex_dt_cov_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	open(2,file=file_name)
	call Bootstrap(real(dt_list,16)/c0,n_Block,2*N0-1,10**3,dt_sd)
	do i = 1,2*N0-1
		write(1,'(999999999F0.9)') sqrt(dt_sd(i,i))
		write(2,'(999999999F30.9)') dt_sd(i,:)
	end do
	close(1)
	close(2)
	
end subroutine Auswertung

subroutine Wurm(N0,N1,M,lattice,dt_list_Wurm)

	implicit none
	integer :: N0,N1,n_data,randomint,check_end,check_end_prop,check_delete
	integer :: i,j,k,t,check_valid,check_valid_new,dt
	integer :: x0,x1,y0,y1,start0,start1,direction
	integer :: n_bonds_x_new,n_bonds_y_new,n_bonds_start,n_bonds_start_new
	integer, dimension(4) :: oldvertex_x,oldvertex_y,newvertex_x,newvertex_y
	integer,dimension(N0,N1,4) :: lattice,lattice_prop
	integer(kind=8), dimension(2*N0-1) :: dt_list_Wurm
	real :: M,p0,proposal,q1,q2,pA,u
	real :: weight,weightx,weighty,weightx_new,weighty_new
	
	!Setze die dt_liste auf 0
	do i = 1,2*N0-1
		dt_list_Wurm(i) = 0_8
	end do
	
	!Schlage einen Startort für den neuen Wurm vor
	start0 = randomint(1,N0)
	start1= randomint(1,N1)
	
	!Wahrscheinlichkeit, eine Quelle/Senke wegzunehmen
	p0 = 1.0/real(N0*N1)
	
	!Bestimme, ob man die Startkonfiguration akzeptiert
	check_end = 1
	q1 = 1.0/real(N0*N1)
	q2 = p0
	weightx = weight(M,lattice(start0,start1,:),0)
	weightx_new = weight(M,lattice(start0,start1,:),1)
	pA = (weightx_new/weightx)*(q2/q1)
	
	call random_number(u)
	if (pA>u) then !Falls die Startkonfiguration akzeptiert wird
		check_end=0
		x0=start0
		x1=start1
		
		!Bestimme die Anzahl Bonds beim Start
		n_bonds_start=0
		do i=1,4
			if (lattice(start0,start1,i)==1) then
				n_bonds_start=n_bonds_start+1
			end if
		end do
		
		!Falls n_bonds_start=0 hat man eine physikalisch valide Startconfig.
		check_valid=0
		if (n_bonds_start==0) then
			check_valid=1
			dt_list_Wurm(N0) = 1_8
		end if
		
	end if
	
	!Wurmloop: Führe den Wurm weiter, bis man zu einer validen "wurmlosen" Kofniguration zurückwechselt
	Wurmloop: do while(check_end==0)
		
		!Falls man einen Loop gemacht hat, schlage mit p0 Wahrscheinlichkeit vor zurück in die wurmlose config. zu wechseln
		oldvertex_x=lattice(x0,x1,:)
		check_end_prop=0
		if (x0==start0 .and. x1==start1) then
			call random_number(u)
			if (p0>u) then
				check_end_prop=1
				y0=0
				y1=0
				newvertex_x = oldvertex_x
				oldvertex_y = oldvertex_x
				newvertex_y = oldvertex_x
				lattice_prop=lattice
			end if
		end if
		
		!Falls man nicht vorschlägt aufzuhören, schlage vor zum Punkt y weiter zu gehen
		check_delete=0 !Testet, ob man einen Bond löscht
		if (check_end_prop==0) then
		
			!Schlage eine neue zufällige Richtung vor
			direction = randomint(1,4)
			
			!Definiere den neuen potentiellen Wurmkopf y sowie die neuen Vertex bei x und y
			select case(direction)
				
			case(1)
			!Definiere den neuen Punkt (y0,y1)
			if (x0 < N0) then
				y0 = x0 + 1
			else !Randbedingung
				y0 = 1
			end if
			y1 = x1
			oldvertex_y = lattice(y0,y1,:)
			!Falls kein Hop in der Richtung existiert, erstelle einen neuen Hop
			if (oldvertex_x(direction)==0) then
				newvertex_x = oldvertex_x + (/ 1,0,0,0 /)
				newvertex_y = oldvertex_y + (/ 0,0,1,0 /)
			else !Falls bereits ein Hop in der Richtung existiert, lösche den Hop
				check_delete = 1
				newvertex_x = oldvertex_x + (/ -1,0,0,0 /)
				newvertex_y = oldvertex_y + (/ 0,0,-1,0 /)
			end if
					
			case(2)
			!Definiere den neuen Punkt (y0,y1)
			y0 = x0
			if (x1 < N1) then
				y1 = x1 + 1
			else !Randbedingung
				y1 = 1
			end if
			oldvertex_y = lattice(y0,y1,:)
			!Falls kein Hop in der Richtung existiert, erstelle einen neuen Hop
			if (oldvertex_x(direction)==0) then
				newvertex_x = oldvertex_x + (/ 0,1,0,0 /)
				newvertex_y = oldvertex_y + (/ 0,0,0,1 /)
			else !Falls bereits ein Hop in der Richtung existiert, lösche den Hop
				check_delete = 1
				newvertex_x = oldvertex_x + (/ 0,-1,0,0 /)
				newvertex_y = oldvertex_y + (/ 0,0,0,-1 /)
			end if
				
			case(3)
			!Definiere den neuen Punkt (y0,y1)
			if (x0 > 1) then
				y0 = x0 - 1
			else !Randbedingung
				y0 = N0
			end if
			y1 = x1
			oldvertex_y = lattice(y0,y1,:)
			!Falls kein Hop in der Richtung existiert, erstelle einen neuen Hop
			if (oldvertex_x(direction)==0) then
				newvertex_x = oldvertex_x + (/ 0,0,1,0 /)
				newvertex_y = oldvertex_y + (/ 1,0,0,0 /)
			else !Falls bereits ein Hop in der Richtung existiert, lösche den Hop
				check_delete = 1
				newvertex_x = oldvertex_x + (/ 0,0,-1,0 /)
				newvertex_y = oldvertex_y + (/ -1,0,0,0 /)
			end if

			case(4)
			!Definiere den neuen Punkt (y0,y1)
			y0 = x0
			if (x1 > 1) then
				y1 = x1 - 1
			else !Randbedingung
				y1 = N1
			end if
			oldvertex_y = lattice(y0,y1,:)
			!Falls kein Hop in der Richtung existiert, erstelle einen neuen Hop
			if (oldvertex_x(direction)==0) then
				newvertex_x = oldvertex_x + (/ 0,0,0,1 /)
				newvertex_y = oldvertex_y + (/ 0,1,0,0 /)
			else !Falls bereits ein Hop in der Richtung existiert, lösche den Hop
				check_delete = 1
				newvertex_x = oldvertex_x + (/ 0,0,0,-1 /)
				newvertex_y = oldvertex_y + (/ 0,-1,0,0 /)
			end if
			
			end select
			
			lattice_prop=lattice
			lattice_prop(x0,x1,:)=newvertex_x
			lattice_prop(y0,y1,:)=newvertex_y
		end if
		
		!Bestimme die Anzahl Bonds bei x,y und start
		n_bonds_x_new = 0
		n_bonds_y_new = 0
		n_bonds_start_new = 0
		do i = 1,4
			if (newvertex_x(i)==1) then
				n_bonds_x_new = n_bonds_x_new+1
			end if
			if (newvertex_y(i)==1) then
				n_bonds_y_new = n_bonds_y_new+1
			end if
			if (lattice_prop(start0,start1,i)==1) then
				n_bonds_start_new = n_bonds_start_new+1
			end if
		end do
		
		!Bestimme die Gewichte bei x
		if (y0==0 .and. y1==0) then !Falls man aufhört, nimm den Wurmkopf weg
			weightx = weight(M,oldvertex_x,1)
			weightx_new = weight(M,newvertex_x,0)
		else if (x0==start0 .and. x1==start1) then !Falls man am Start ist, bleibt der Wurmkopf
			weightx = weight(M,oldvertex_x,1)
			weightx_new = weight(M,newvertex_x,1)
			if (n_bonds_x_new==3) then !Hinterlasse keinen 3-er Vertex beim Start
				weightx_new = 0.0
			end if
		else !Sonst bewegt sich der Wurmkopf weiter zu y
			weightx = weight(M,oldvertex_x,1)
			weightx_new = weight(M,newvertex_x,0)
		end if
		
		!Bestimme die Gewichte bei y
		if (y0==0 .and. y1==0) then !Falls man aufhört, ist y=0.
			weighty = 1.0
			weighty_new = 1.0
		else if (y0==start0 .and. y1==start1) then !Falls man zum Start geht, bleibt der Wurmkopf
			weighty = weight(M,oldvertex_y,1)
			weighty_new = weight(M,newvertex_y,1)
		else !Sonst bewegt sich der Wurmkopf weiter zu y
			weighty = weight(M,oldvertex_y,0)
			weighty_new = weight(M,newvertex_y,1)
		end if
		
		!Bestimme, ob die neue config valide ist
		check_valid_new=0
		if (y0==start0 .and. y1==start1) then
			if (n_bonds_start_new==0) then
				if (n_bonds_x_new==0 .or. n_bonds_x_new==2) then
					check_valid_new=1
				end if
			end if
		else if (x0==start0 .and. x1==start1) then
			if (y0==0 .and. y1==0 .and. n_bonds_start_new==0) then
				check_valid_new=1
			else
				if (n_bonds_start_new==1 .and. n_bonds_y_new==1) then
					check_valid_new=1
				end if
			end if
		else
			if (n_bonds_start_new==1 .and. n_bonds_y_new==1) then
				if (n_bonds_x_new==2 .or. n_bonds_x_new==0) then
					check_valid_new=1
				end if
			end if
		end if
		
		!Bestimme die Proposal-Wahrscheinlichkeiten, von x->y zu gehen, sowie die Akzeptanzwahrscheinlichkeit
		q1 = proposal(N0,N1,x0,x1,y0,y1,p0,oldvertex_x,newvertex_x)
		q2 = proposal(N0,N1,y0,y1,x0,x1,p0,newvertex_y,oldvertex_y)
		pA = (weightx_new/weightx) * (weighty_new/weighty) * (q2/q1)
		
		!Bestimme, ob der Schritt akzeptiert wird
		call random_number(u)
		if ( pA > u ) then
			
			if (y0==0 .and. y1==0) then !Falls man aufhört
				check_end=1
				exit Wurmloop
			else !Falls man mit dem Wurm weitermacht
		
				lattice(x0,x1,:) = newvertex_x
				lattice(y0,y1,:) = newvertex_y
				x0 = y0
				x1 = y1
				
				n_bonds_start=n_bonds_start_new
				check_valid = check_valid_new
				
			end if
		end if
		
		!Führe die dt_liste weiter, falls man eine physikalisch valide config hat
		if (check_valid==1) then
			dt=x0-start0
			dt_list_Wurm(dt+N0) = dt_list_Wurm(dt+N0)+1_8
		end if
		
	end do Wurmloop

end subroutine Wurm

subroutine Bootstrap(liste,n_Block,n_parameters,n_tries,corr) !Nimmt eine Liste von n_Block (real) Messungen von n_parameters Messgrössen und schätze die Kovarianzmatrix ab

	implicit none
	integer :: n_Block,n_parameters,i,j,k,t,randomint,n_tries
	real(kind=16),dimension(n_parameters) :: summe1,bootstrap_liste
	real(kind=16),dimension(n_parameters,n_parameters) :: summe2
	real(kind=16),dimension(n_Block,n_parameters) :: liste
	real, dimension(n_parameters,n_parameters) :: corr
	
	!Setzte die Summen auf 0
	do i = 1,n_parameters
		summe1(i) = 0.0_16
		do j = 1,n_parameters
			summe2(i,j) = 0.0_16
		end do
	end do
	
	!Erstelle n_tries viele Bootstrap-Listen
	do k = 1,n_tries
	
		!Bootstrap-Liste: Wähle für jeden Eintrag einen zufälligen Block und summiere sie auf
		do i = 1,n_parameters
			bootstrap_liste(i) = 0.0_16
		end do
		do i = 1,n_Block
			t = randomint(1,n_Block)
			bootstrap_liste(:) = bootstrap_liste(:) + liste(t,:)
		end do
		
		!summe1: Summiere jeden Parameter einzeln auf
		do i = 1,n_parameters
			summe1(i) = summe1(i) + bootstrap_liste(i)
			!summe2: Summiere die Parameter miteinander multipliziert auf
			do j = 1,n_parameters
				summe2(i,j) = summe2(i,j) + bootstrap_liste(i)*bootstrap_liste(j)
			end do
		end do
		
	end do
	
	!Bestimme die Kovarianzmatrix
	do i = 1,n_parameters
		do j = 1,n_parameters
			corr(i,j) = summe2(i,j)/real(n_tries,16) - (summe1(i)/real(n_tries,16))*(summe1(j)/real(n_tries,16))
		end do
	end do
	
end subroutine Bootstrap

function weight(M,mu,check_x)!Gewicht von einem Vertex. check_x=1: Quelle/Senke am Vertex.

	implicit none
	integer :: summe,check_x
	real :: M,weight
	integer, dimension(4) :: mu
	
	summe = sum(mu)
	
	if (check_x==1) then !Falls man eine Quelle/Senke am Vertex hat
	
		select case (summe)
		
			case(0) !Quelle und Senke am gleichen Ort
			weight = 4.0 !w=4, da man verbundene und unverbundene zusammenzählt
			
			case(1)
			weight = 1.0
			
			case(2)
			if (mu(1)==1 .and. mu(3)==1) then !Gerade horizontale Linie
				weight = 1.0
			else if(mu(2)==1 .and. mu(4)==1) then !Gerade vertikale Linie
				weight = 1.0
			else !Ecke
				weight = 0.5
			end if
			
			case(3)
			weight = 0.25**(1.0/3.0)
			
			case default
			weight = 0.0
		
		end select
	
	else !Falls man keine Quellen/Senken hat
	
		select case (summe)
		
			case(0)
			weight=M**2
			
			case(1)
			weight=0.0
			
			case(2)
			if (mu(1)==1 .and. mu(3)==1) then !Gerade horizontale Linie
				weight = 1.0
			else if(mu(2)==1 .and. mu(4)==1) then !Gerade vertikale Linie
				weight = 1.0
			else !Ecke
				weight = 0.5
			end if
			
			case default
			weight = 0.0
		
		end select
	
	end if

end function weight

function proposal(N0,N1,x0,x1,y0,y1,p0,oldvertex_x,newvertex_x) !Gibt die proposal Wahrscheinlichkeit, von oldvertex_x -> newvertex_x zu gehen.

	implicit none
	integer :: N0,N1,x0,x1,y0,y1,n_bonds,i
	integer, dimension(4) :: oldvertex_x,newvertex_x
	real :: proposal,p0
	
	if ( x0 == 0 .and. x1 == 0 ) then
		proposal = 1.0 / real(N0*N1)
	else
		!Bestimme die Anzahl Bonds beim Wurmkopf
		n_bonds = 0
		do i = 1,4
			if ( oldvertex_x(i) /= 0 ) then
				n_bonds = n_bonds + 1
			end if
		end do
	
		select case(n_bonds)
		
		case(0)!Valider Vertex
		if ( y0 == 0 .and. y1 == 0 ) then
			proposal = p0
		else
			proposal = (1.0-p0) / 4.0
		end if
		
		case(1) !1-point Vertex, alle Richtungen möglich
		proposal = 1.0 / 4.0
		
		
		case(2)!Valider Vertex
		if ( y0 == 0 .and. y1 == 0 ) then
			proposal = p0
		else
			proposal = (1.0-p0) / 4.0
		end if
		
		case(3) !3-point Vertex, nur drei Richtungen möglich (aber alle werden vorgeschlagen)
		proposal = 1.0 / 4.0
		
		case default
		proposal = 0.0
		
		end select
	end if

end function proposal

function topology(lattice,N0,N1) !Gibt die Topologie-Klasse einer (wurmlosen) Konfiguration zurück

	implicit none
	integer :: topology,N0,N1,i,n_0,n_1
	integer, dimension(N0,N1,4) :: lattice
	
	n_0 = 0
	n_1 = 0
	
	!Anzahl Hops über den unteren Rand
	do i = 1,N0
	
		if ( lattice(i,1,4)==1 ) then
			n_1 = n_1 + 1
		end if
	
	end do

	!Anzahl Hops über den linken Rand
	do i = 1,N1
	
		if ( lattice(1,i,3)==1 ) then
			n_0 = n_0 + 1
		end if
		
	end do

	n_0 = modulo(n_0,2)
	n_1 = modulo(n_1,2)
	!Topologische Klasse = (n_0,n_1)
	
	if ( n_0 == 0 .and. n_1 == 0) then
		topology = 1
	else if ( n_0 == 1 .and. n_1 == 0) then
		topology = 2
	else if ( n_0 == 0 .and. n_1 == 1) then
		topology = 3
	else if ( n_0 == 1 .and. n_1 == 1) then
		topology = 4
	end if

end function topology

function randomint(a,b) !Zufälliges Integer zwischen a und b (inklusive a und b)

	implicit none
	integer :: randomint,a,b
	real :: u
	
	call random_number(u)
	randomint = a + floor(real(b-a+1)*u)
	
end function randomint
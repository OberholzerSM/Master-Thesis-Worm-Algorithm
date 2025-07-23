program Wurm163Vertex
	
	implicit none
	integer :: N0,N1,n_data,option,n_masse,n_Gitter,n_Block
	integer :: i,j,k,k0,k1
	real :: M,M_min,M_max,start_time,end_time
	real :: Z,Z_abs
	real,dimension(2) :: sign_weights
	real,dimension(5) :: isospin_weights
	real,dimension(5,3) :: dt_weights
	real(kind=16), dimension(36,36) :: T,T_abs
	real(kind=16), dimension(5,36,36) :: S_source,S_sink,S0,S_source_abs,S_sink_abs,S0_abs
	character(len=10) :: file_M,file_k,file_N0,file_N1
	character(len=100) :: file_name
	
	call cpu_time(start_time)
	option=1
	if (option<3 .and. option>0) then !Lese die Gittergrösse ein
		call getarg(1,file_k)
		read(file_k,*) k0
		call getarg(2,file_k)
		read(file_k,*) k1
	end if
	if (option==1) then !Lese die Masse ein
		call getarg(3,file_M)
		read(file_M,*) M
	end if
	N0=2**k0
	N1=2**k1
	n_masse = 20
	M_min = 0.0
	M_max = 4.0
	n_Gitter = 5
	n_Block = 100
	n_data = 10**5

	select case(option)
	case(-2) !Parameter File
		call make_txt_file()
	
	case(-1) !Transfermatrix 2x1
		
		do i = 1,n_masse
			print*,real(100*i)/n_masse,"%"
			M = real(i)*(M_max-M_min)/real(n_masse) + M_min
			write(file_M,"(f10.2)") M
			file_name = "163Vertex_T_M="// trim(adjustl(file_M)) //".dat"
			open(100,file=file_name)
			file_name = "163Vertex_S_source_M="// trim(adjustl(file_M)) //".dat"
			open(101,file=file_name)
			file_name = "163Vertex_S_sink_M="// trim(adjustl(file_M)) //".dat"
			open(102,file=file_name)
			file_name = "163Vertex_S0_M="// trim(adjustl(file_M)) //".dat"
			open(103,file=file_name)
			file_name = "163Vertex_T_abs_M="// trim(adjustl(file_M)) //".dat"
			open(104,file=file_name)
			file_name = "163Vertex_S_source_abs_M="// trim(adjustl(file_M)) //".dat"
			open(105,file=file_name)
			file_name = "163Vertex_S_sink_abs_M="// trim(adjustl(file_M)) //".dat"
			open(106,file=file_name)
			file_name = "163Vertex_S0_abs_M="// trim(adjustl(file_M)) //".dat"
			open(107,file=file_name)
			
			call Transfermatrix(M,T,T_abs,S_source,S_sink,S0,S_source_abs,S_sink_abs,S0_abs)
			do j = 1,36
				write(100,'(36(F50.9,:,1X),/)') (T(j,:))
				write(104,'(36(F50.9,:,1X),/)') (T_abs(j,:))
			end do
			write(100,*) " "
			write(104,*) " "
			
			do k=1,5
				do j = 1,36
					write(101,'(36(F50.9,:,1X),/)') (S_source(k,j,:))
					write(102,'(36(F50.9,:,1X),/)') (S_sink(k,j,:))
					write(103,'(36(F50.9,:,1X),/)') (S0(k,j,:))
					write(105,'(36(F50.9,:,1X),/)') (S_source_abs(k,j,:))
					write(106,'(36(F50.9,:,1X),/)') (S_sink_abs(k,j,:))
					write(107,'(36(F50.9,:,1X),/)') (S0_abs(k,j,:))
				end do
				write(101,*) " "
				write(102,*) " "
				write(103,*) " "
				write(105,*) " "
				write(106,*) " "
				write(107,*) " "
			end do
			
			close(100)
			close(101)
			close(102)
			close(103)
			close(104)
			close(105)
			close(106)
			close(107)
			
		end do
		
	case(0) !Erstelle die Normierung für das 2x2 Gittergrösse
		
		file_name = "163Vertex_Z_N=2x2.dat"
		open(100,file=file_name)
		file_name = "163Vertex_weights_isospin_N=2x2.dat"
		open(101,file=file_name)
		file_name = "163Vertex_weights_sign_N=2x2.dat"
		open(102,file=file_name)
		file_name = "163Vertex_weights_dt_N=2x2.dat"
		open(103,file=file_name)
		
		do i = 1,n_masse
			M = real(i)*2.0/real(n_masse)
			call Normierung_2x2(M,Z,Z_abs,dt_weights,isospin_weights,sign_weights)	
			write(100,*) Z,Z_abs
			write(101,*) isospin_weights
			write(102,*) sign_weights
			do j = 1,5
				write(103,*) dt_weights(j,:)
			end do
			write(103,*)
		end do
		
		close(100)
		close(101)
		close(102)
		close(103)
		
	case(1) !Wurm Fixe Gittergrösse und Masse
		call Auswertung(N0,N1,M,n_data,n_Block)
	case(2) !Wurm Fixe Gittergrösse und variable Masse
		do i = 1,n_masse
			M = real(i)*2.0/real(n_masse)
			call Auswertung(N0,N1,M,n_data,n_Block)
		end do
	case(3) !Wurm Variable Gittergrösse und Masse
	
		do j = 1,n_Gitter
			N0 = 2**j
			N1 = 2
			do i = 1,n_masse
				M = real(i)*2.0/real(n_masse)
				call Auswertung(N0,N1,M,n_data,n_Block)
			end do
		end do

	end select
	
	call cpu_time(end_time)
	print*,"Time: ",end_time-start_time
	
end program Wurm163Vertex

subroutine make_txt_file() !Gibt ein txt File mit den Parametern für SLURM aus

	implicit none
	integer :: n_masse,n_Gitter,i,j,k
	real :: M,M_min,M_max
	
	n_Gitter = 5
	n_masse = 20
	
	M_min=0.0
	M_Max=2.0
	open(1,file="params_N0x2_Mk.txt")
	write(1,*) 0
	do j = 1,n_Gitter
		do i = 1,n_masse
			M = real(i)*(M_max-M_min)/real(n_masse) + M_min
			write(1,*) j,1,M
		end do
	end do
	close(1)
	
	M_min=2.0
	M_Max=4.0
	open(1,file="params_N0x2_MG.txt")
	write(1,*) 0
	do j = 1,n_Gitter
		do i = 1,n_masse
			M = real(i)*(M_max-M_min)/real(n_masse) + M_min
			write(1,*) j,1,M
		end do
	end do
	close(1)
	
	M_min=0.0
	M_Max=2.0
	open(1,file="params_32xN1_Mk.txt")
	write(1,*) 0
	do j = 2,n_Gitter
		do i = 1,n_masse
			M = real(i)*(M_max-M_min)/real(n_masse) + M_min
			write(1,*) n_Gitter,j,M
		end do
	end do
	close(1)
	
	M_min=2.0
	M_Max=4.0
	open(1,file="params_32xN1_MG.txt")
	write(1,*) 0
	do j = 2,n_Gitter
		do i = 1,n_masse
			M = real(i)*(M_max-M_min)/real(n_masse) + M_min
			write(1,*) n_Gitter,j,M
		end do
	end do
	close(1)
	
	M_min=0.0
	M_Max=2.0
	open(1,file="params_32x64_Mk.txt")
	write(1,*) 0
	do i = 1,n_masse
		M = real(i)*(M_max-M_min)/real(n_masse) + M_min
		write(1,*) 5,6,M
	end do
	close(1)
	
	M_min=2.0
	M_Max=4.0
	open(1,file="params_32x64_MG.txt")
	write(1,*) 0
	do i = 1,n_masse
		M = real(i)*(M_max-M_min)/real(n_masse) + M_min
		write(1,*) 5,6,M
	end do
	close(1)

end subroutine make_txt_file

subroutine Normierung_2x2(M,Z,Z_abs,dt_weights,isospin_weights,sign_weights) !Gibt die Normierung Z fuer das 2x2 Gitter fuer eine gegeben Masse an

	implicit none
	integer :: i,j,isospin,t,check_x,check_start,dt,counter,counter_closed
	integer :: i1,i2,i3,i4,i5,i6,i7,i8,i_source,j_source,i_sink,j_sink,hope_mode
	real :: M,weight_163Vertex,weight_vertex
	real :: Z,Z_abs,weight_lattice
	real,dimension(2) :: sign_weights
	real,dimension(5) :: isospin_weights
	real,dimension(5,3) :: dt_weights
	integer, dimension(2,2,4) :: lattice
	integer, dimension(4251,2,2,4) :: closed_configs_list
	integer, dimension(5) :: counter_propagator
	integer, dimension(15792,2,2,4) :: hopetype1_list
	integer, dimension(15792,2,2,4) :: hopetype2_list
	integer, dimension(21524,2,2,4) :: hopetype3_list
	integer, dimension(21524,2,2,4) :: hopetype4_list
	integer, dimension(7224,2,2,4) :: hopetype5_list
	integer, dimension(15792,2,2) :: positions_1
	integer, dimension(15792,2,2) :: positions_2
	integer, dimension(21524,2,2) :: positions_3
	integer, dimension(21524,2,2) :: positions_4
	integer, dimension(7224,2,2) :: positions_5
	
	do i = 1,3
		do j = 1,5
			dt_weights(j,i) = 0.0
		end do
	end do
	do i = 1,5
		isospin_weights(i) = 0.0
	end do
	do i = 1,2
		sign_weights(i) = 0.0
	end do

	counter = 0
	counter_closed = 0
	do i=1,5
		counter_propagator(i)=0
	end do

	Z = 0.0_8
	Z_abs = 0.0_8
	do i1 = 0,5
		do i2 = 0,5
			do i3 = 0,5
				do i4 = 0,5
					do i5 = 0,5
						do i6 = 0,5
							do i7 = 0,5
								do i8 = 0,5
									
									!Drucke den Fortschritt aus
									counter = counter + 1
									do i = 1,100
										if ( i*(6**8)/100==counter ) then
											print*,i,"%, M=",M
										end if
									end do
									
									!Definiere die Konfiguration
									lattice(1,1,:) = (/i1,i2,i3,i4/)
									lattice(2,2,:) = (/i5,i6,i7,i8/)
									lattice(1,2,:) = (/i7,i4,i5,i2/)
									lattice(2,1,:) = (/i3,i8,i1,i6/)
									do i = 1,4
										if ( lattice(1,2,i) == 3 ) then
											lattice(1,2,i) = 4
										else if ( lattice(1,2,i) == 4 ) then
											lattice(1,2,i) = 3
										end if
										if ( lattice(2,1,i) == 3 ) then
											lattice(2,1,i) = 4
										else if ( lattice(2,1,i) == 4 ) then
											lattice(2,1,i) = 3
										end if
									end do
									
									!Gewicht der Konfiguration ohne Quellen/Senken
									weight_lattice = 1.0
									do i = 1,2
										do j = 1,2
											weight_vertex = weight_163Vertex(M,lattice(i,j,:),1,0,0)
											weight_lattice = weight_lattice * weight_vertex
										end do
									end do
									
									if ( abs(weight_lattice) > 0.0 ) then
										counter_closed = counter_closed + 1
										closed_configs_list(counter_closed,:,:,:) = lattice
										Z = Z + weight_lattice
										Z_abs = Z_abs + abs(weight_lattice)
										
										!Isospin
										t = isospin(lattice,2,2) + 3 !3=N1+1
										isospin_weights(t) = isospin_weights(t) + weight_lattice
										
										!Sign
										if (weight_lattice > 0.0) then
											sign_weights(1) = sign_weights(1) + weight_lattice
										else
											sign_weights(2) = sign_weights(2) + weight_lattice
										end if
										
									end if
									
									!Gewicht der Konfiguration mit Quellen und Senken
									do hope_mode = 1,5
										do i_sink = 1,2
											do j_sink = 1,2
												do i_source = 1,2
													do j_source = 1,2
														
														weight_lattice = 1.0
														do i = 1,2
															do j = 1,2
																check_start=0
																check_x = 0
																if (i==i_source .and. j==j_source) then
																	check_start=1
																end if
																if (i==i_sink .and. j==j_sink) then
																	check_x=1
																end if
																weight_vertex = weight_163Vertex(M,lattice(i,j,:),hope_mode,check_x,check_start)
																weight_lattice = weight_lattice * weight_vertex
															end do
														end do
														
														dt = i_sink - i_source
														dt_weights(hope_mode,dt+2) = dt_weights(hope_mode,dt+2) + weight_lattice
														
														if ( abs(weight_lattice) > epsilon(0.0) ) then
															counter_propagator(hope_mode)=counter_propagator(hope_mode)+1
															
															select case(hope_mode)
															case(1)
																hopetype1_list(counter_propagator(hope_mode),:,:,:) = lattice
																positions_1(counter_propagator(hope_mode),1,:) = (/i_source,j_source/)
																positions_1(counter_propagator(hope_mode),2,:) = (/i_sink,j_sink/)
															case(2)
																hopetype2_list(counter_propagator(hope_mode),:,:,:) = lattice
																positions_2(counter_propagator(hope_mode),1,:) = (/i_source,j_source/)
																positions_2(counter_propagator(hope_mode),2,:) = (/i_sink,j_sink/)
															case(3)
																hopetype3_list(counter_propagator(hope_mode),:,:,:) = lattice
																positions_3(counter_propagator(hope_mode),1,:) = (/i_source,j_source/)
																positions_3(counter_propagator(hope_mode),2,:) = (/i_sink,j_sink/)
															case(4)
																hopetype4_list(counter_propagator(hope_mode),:,:,:) = lattice
																positions_4(counter_propagator(hope_mode),1,:) = (/i_source,j_source/)
																positions_4(counter_propagator(hope_mode),2,:) = (/i_sink,j_sink/)
															case(5)
																hopetype5_list(counter_propagator(hope_mode),:,:,:) = lattice
																positions_5(counter_propagator(hope_mode),1,:) = (/i_source,j_source/)
																positions_5(counter_propagator(hope_mode),2,:) = (/i_sink,j_sink/)
															end select
										
														end if
														
													end do
												end do
											end do
										end do
									end do
									
								end do
							end do
						end do
					end do
				end do
			end do
		end do
	end do
	print*,""
	!print*,counter_closed
	!print*,counter_propagator
	
end subroutine Normierung_2x2

subroutine Transfermatrix(M,T,T_abs,S_source,S_sink,S0,S_source_abs,S_sink_abs,S0_abs) !Berechnet die Transfermatrizen für ein gegebenes M für das 2x1 Gitter

	implicit none
	integer :: i,j,k,i_sink,i_source,hope_mode
	integer :: i1,i2,i3,i4,i5,i6,i7,i8,i3_l,i7_l
	integer,dimension(2,4) :: lattice
	real :: M,weight_163Vertex,weight_lattice
	real(kind=16), dimension(36,36) :: T,T_abs
	real(kind=16), dimension(5,36,36) :: S_source,S_sink,S0,S_source_abs,S_sink_abs,S0_abs
	
	do i=1,36
		do j = 1,36
			T(i,j)=0.0_16
			T_abs(i,j)=0.0_16
			do k = 1,5
				S_source(k,i,j)=0.0_16
				S_sink(k,i,j)=0.0_16
				S0_abs(k,i,j)=0.0_16
				S_source_abs(k,i,j)=0.0_16
				S_sink_abs(k,i,j)=0.0_16
				S0_abs(k,i,j)=0.0_16
			end do
		end do
	end do
	
	do i1 = 0,5
		do i2 = 0,5
			do i3 = 0,5
				do i4 = 0,5
					do i5 = 0,5
						do i7 = 0,5
							
							!Definiere die 1x2 Lattice
							lattice(2,:) = (/i1,i2,i3,i4/)
							if (i2==3) then
								i8=4
							else if (i2==4) then
								i8=3
							else
								i8=i2
							end if
							if (i4==3) then
								i6=4
							else if (i4==4) then
								i6=3
							else
								i6=i4
							end if
							lattice(1,:) = (/i5,i6,i7,i8/)
							
							!Index linke und rechte Seite (vertausche für eine Seite bonds 3 und 4)
							if (i7==3) then
								i7_l=4
							else if (i7==4) then
								i7_l=3
							else
								i7_l=i7
							end if
							if (i3==3) then
								i3_l=4
							else if (i3==4) then
								i3_l=3
							else
								i3_l=i3
							end if
							i = (i7_l+1) + 6*i3_l
							j = (i5+1) + 6*i1
							
							!Bestimme das Gewicht ohne Quellen
							weight_lattice = 1.0
							do k = 1,2
								weight_lattice = weight_lattice * weight_163Vertex(M,lattice(k,:),0,0,0)
							end do
							T(i,j)=T(i,j)+real(weight_lattice,16)
							T_abs(i,j)=T_abs(i,j)+real(abs(weight_lattice),16)
							
							do hope_mode = 1,5
								!Bestimme das Gewicht mit einer Quelle
								weight_lattice = 1.0
								weight_lattice = weight_lattice * weight_163Vertex(M,lattice(2,:),hope_mode,0,1)
								weight_lattice = weight_lattice * weight_163Vertex(M,lattice(1,:),hope_mode,0,0)
								S_source(hope_mode,i,j) = S_source(hope_mode,i,j) + real(weight_lattice,16)
								S_source_abs(hope_mode,i,j) = S_source_abs(hope_mode,i,j) + real(abs(weight_lattice),16)
								
								weight_lattice = 1.0
								weight_lattice = weight_lattice * weight_163Vertex(M,lattice(2,:),hope_mode,0,0)
								weight_lattice = weight_lattice * weight_163Vertex(M,lattice(1,:),hope_mode,0,1)
								S_source(hope_mode,i,j) = S_source(hope_mode,i,j) + real(weight_lattice,16)
								S_source_abs(hope_mode,i,j) = S_source_abs(hope_mode,i,j) + real(abs(weight_lattice),16)
								
								!Bestimme das Gewicht mit einer Senke
								weight_lattice = 1.0
								weight_lattice = weight_lattice * weight_163Vertex(M,lattice(2,:),hope_mode,1,0)
								weight_lattice = weight_lattice * weight_163Vertex(M,lattice(1,:),hope_mode,0,0)
								S_sink(hope_mode,i,j) = S_sink(hope_mode,i,j) + real(weight_lattice,16)
								S_sink_abs(hope_mode,i,j) = S_sink_abs(hope_mode,i,j) + real(abs(weight_lattice),16)
								
								weight_lattice = 1.0
								weight_lattice = weight_lattice * weight_163Vertex(M,lattice(2,:),hope_mode,0,0)
								weight_lattice = weight_lattice * weight_163Vertex(M,lattice(1,:),hope_mode,1,0)
								S_sink(hope_mode,i,j) = S_sink(hope_mode,i,j) + real(weight_lattice,16)
								S_sink_abs(hope_mode,i,j) = S_sink_abs(hope_mode,i,j) + real(abs(weight_lattice),16)
								
								!Bestimme das Gewicht mit einer Quelle und einer Senke
								weight_lattice = 1.0
								weight_lattice = weight_lattice * weight_163Vertex(M,lattice(2,:),hope_mode,1,0)
								weight_lattice = weight_lattice * weight_163Vertex(M,lattice(1,:),hope_mode,0,1)
								S0(hope_mode,i,j) = S0(hope_mode,i,j) + real(weight_lattice,16)
								S0_abs(hope_mode,i,j) = S0_abs(hope_mode,i,j) + real(abs(weight_lattice),16)
								
								weight_lattice = 1.0
								weight_lattice = weight_lattice * weight_163Vertex(M,lattice(2,:),hope_mode,0,1)
								weight_lattice = weight_lattice * weight_163Vertex(M,lattice(1,:),hope_mode,1,0)
								S0(hope_mode,i,j) = S0(hope_mode,i,j) + real(weight_lattice,16)
								S0_abs(hope_mode,i,j) = S0_abs(hope_mode,i,j) + real(abs(weight_lattice),16)
								
								weight_lattice = 1.0
								weight_lattice = weight_lattice * weight_163Vertex(M,lattice(2,:),hope_mode,1,1)
								weight_lattice = weight_lattice * weight_163Vertex(M,lattice(1,:),hope_mode,0,0)
								S0(hope_mode,i,j) = S0(hope_mode,i,j) + real(weight_lattice,16)
								S0_abs(hope_mode,i,j) = S0_abs(hope_mode,i,j) + real(abs(weight_lattice),16)
								
								weight_lattice = 1.0
								weight_lattice = weight_lattice * weight_163Vertex(M,lattice(2,:),hope_mode,0,0)
								weight_lattice = weight_lattice * weight_163Vertex(M,lattice(1,:),hope_mode,1,1)
								S0(hope_mode,i,j) = S0(hope_mode,i,j) + real(weight_lattice,16)
								S0_abs(hope_mode,i,j) = S0_abs(hope_mode,i,j) + real(abs(weight_lattice),16)
			
							end do
							
						end do
					end do
				end do
			end do
		end do
	end do

end subroutine Transfermatrix

subroutine Auswertung(N0,N1,M,n_data,n_Block)
	
	implicit none
	integer :: N0,N1,n_data,n_Block,l_Block,Block_index,sign_lattice,sign_weight
	integer :: randomint,isospin,t_isospin,i,j,k,hope_mode
	integer, dimension(N0,N1,4) :: lattice
	integer(kind=8), dimension(5,2*N0-1) :: dt_number,dt_number_abs
	integer(kind=8), dimension(2*N0-1) :: dt_list_Wurm,dt_list_Wurm_abs
	integer(kind=8), dimension(n_Block,5*(2*N0-1)) :: dt_list,dt_list_abs
	real, dimension(5*(2*N0-1),5*(2*N0-1)) :: dt_sd,dt_sd_abs
	integer(kind=8), dimension(2*N1+1) :: isospin_number,isospin_number_abs
	integer(kind=8), dimension(n_Block,2*N1+1) :: isospin_list,isospin_list_abs
	real, dimension(2*N1+1,2*N1+1) :: isospin_sd,isospin_sd_abs
	integer(kind=8), dimension(2) :: sign_number
	integer(kind=8), dimension(n_Block,2) :: sign_list
	real, dimension(2,2) :: sign_sd
	real :: M,weight,weight_lattice
	real(kind=16) :: c0
	character(len=10) :: file_M,file_N0,file_N1,file_i
	character(len=100) :: file_name
	
	write(file_M,"(f10.2)") M
	write(file_N0,"(i3)") N0
	write(file_N1,"(i3)") N1

	do i=1,5
		write(file_i,"(i3)") i
		file_name = "163Vertex_dt_Blocklist_Meson"// trim(adjustl(file_i)) // "_M=" // trim(adjustl(file_M)) //&
		",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
		open(i,file=file_name)
	end do
	file_name = "163Vertex_isospin_Blocklist_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	open(6,file=file_name)
	file_name = "163Vertex_sign_Blocklist_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	open(7,file=file_name)
	file_name = "163Vertex_last_lattice_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	open(8,file=file_name)
	
	do i=1,5
		write(file_i,"(i3)") i
		file_name = "163Vertex_dt_Blocklist_abs_Meson"// trim(adjustl(file_i)) // "_M=" // trim(adjustl(file_M)) //&
		",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
		open(i+8,file=file_name)
	end do
	file_name = "163Vertex_isospin_Blocklist_abs_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	open(14,file=file_name)
	
	l_Block=n_data/n_Block
	
	!Startkonfiguration
	do i = 1,N0
		do j = 1,N1
			do k = 1,4
				lattice(i,j,k) = 0
			end do
			write(4,'(2000000000I9)') lattice(i,j,:)
		end do
	end do
	sign_lattice = 1
	
	!Zeitabstände Propagator
	do i = 1,5
		do j = 1,2*N0-1
			dt_number(i,j) = 0_8
			dt_list(1,j+(i-1)*(2*N0-1)) = 0_8
			dt_number_abs(i,j) = 0_8
			dt_list_abs(i,j+(i-1)*(2*N0-1)) = 0_8
		end do
	end do
	
	!Isospin
	t_isospin = isospin(lattice,N0,N1)+1+N1
	do i = 1,2*N1+1
		isospin_number(i) = 0_8
		isospin_list(1,i) = 0_8
		isospin_number_abs(i) = 0_8
		isospin_list_abs(1,i) = 0_8
	end do
	isospin_number(t_isospin) = 1_8
	isospin_list(1,t_isospin) = 1_8
	isospin_number_abs(t_isospin) = 1_8
	isospin_list_abs(1,t_isospin) = 1_8
	
	!Vorzeichen
	if (sign_lattice==1) then
		sign_number(1) = 1_8
		sign_number(2) = 0_8
		sign_list(1,1) = 1_8
		sign_list(1,2) = 0_8
	else if (sign_lattice==-1) then
		sign_number(1) = 0_8
		sign_number(2) = 1_8
		sign_list(1,1) = 0_8
		sign_list(1,2) = 1_8
	end if
	
	!Erzeuge n_data viele Konfigurationen
	do k = 2,n_data
		
		!Drucke den Fortschritt aus
		do i = 1,100
			if ( i*n_data/100==k ) then
				print*,i,"%, M=",file_M,", N=",file_N0,"x",file_N1
			end if
		end do
		
		!Generiere eine neue geschlossene Konfiguration
		call Wurm(N0,N1,M,lattice,sign_lattice,hope_mode,dt_list_Wurm,dt_list_Wurm_abs)
		
		!Block-Index
		Block_index = 1+(k-1)/l_Block
		
		!dt-Statistik
		dt_number(hope_mode,:) = dt_number(hope_mode,:) + dt_list_Wurm
		dt_list(Block_index,((hope_mode-1)*(2*N0-1)):(hope_mode*(2*N0-1))) = &
			dt_list(Block_index,((hope_mode-1)*(2*N0-1)):(hope_mode*(2*N0-1))) + dt_list_Wurm
		dt_number_abs(hope_mode,:) = dt_number_abs(hope_mode,:) + dt_list_Wurm_abs
		dt_list_abs(Block_index,((hope_mode-1)*(2*N0-1)):(hope_mode*(2*N0-1))) = &
			dt_list_abs(Block_index,((hope_mode-1)*(2*N0-1)):(hope_mode*(2*N0-1))) + dt_list_Wurm_abs
		
		!Isospin Statistik
		t_isospin = isospin(lattice,N0,N1)+1+N1
		isospin_number(t_isospin) = isospin_number(t_isospin) + sign_lattice
		isospin_list(Block_index,t_isospin) = isospin_list(Block_index,t_isospin) + sign_lattice
		isospin_number_abs(t_isospin) = isospin_number_abs(t_isospin) + 1_8
		isospin_list_abs(Block_index,t_isospin) = isospin_list_abs(Block_index,t_isospin) + 1_8
		
		!Vorzeichen Statistik
		if (sign_lattice==1) then
			sign_number(1) = sign_number(1)+1_8
			sign_list(Block_index,1) = sign_list(Block_index,1)+1_8
		else if (sign_lattice==-1) then
			sign_number(2) = sign_number(2)+1_8
			sign_list(Block_index,2) = sign_list(Block_index,2)+1_8
		end if
		
		!Drucke den Blockwert aus, sobald er fertig ist
		if ( Block_index*l_Block==k ) then
			
			do i = 1,5
				write(i,'(2000000000I16)') dt_list(Block_index,((i-1)*(2*N0-1)):(i*(2*N0-1)))
			end do
			write(6,'(2000000000I9)') isospin_list(Block_index,:)
			write(7,'(2000000000I9)') sign_list(Block_index,:)
			do i = 1,5
				write(i+8,'(2000000000I16)') dt_list_abs(Block_index,((i-1)*(2*N0-1)):(i*(2*N0-1)))
			end do
			write(14,'(2000000000I9)') isospin_list_abs(Block_index,:)
			
			!Setze die Zähler für den nächsten Block auf 0
			if (Block_index /= n_block) then
				do i=1,2*N1+1
					isospin_list(Block_index+1,i)=0_8
					isospin_list_abs(Block_index+1,i)=0_8
				end do
				do i = 1,2
					sign_list(Block_index+1,i)=0_8
				end do
				do i = 1,5
					do j = 1,2*N0-1
						dt_list(Block_index+1,j+(i-1)*(2*N0-1))=0_8
						dt_list_abs(Block_index+1,j+(i-1)*(2*N0-1))=0_8
					end do
				end do
			end if
			
			!Drucke die letzte generierte Lattice aus
			do i = 1,N0
				do j = 1,N1
					write(8,'(2000000000I9)') lattice(i,j,:)
				end do
			end do
			rewind(8)
			
		end if
	
	end do
	
	do i = 1,14
		close(i)
	end do
	
	!Drucke aus, wie häufig die Vorzeichen vorgekommen sind
	file_name = "163Vertex_sign_number_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	open(1,file=file_name)
	do i = 1,2
		write(1,'(2000000000I9)') sign_number(i)
	end do
	close(1)
	
	!Bestimme den Fehler der Vorzeichen
	print*,"Berechne Fehler Vorzeichen..."
	file_name = "163Vertex_sign_sd_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	open(1,file=file_name)
	file_name = "163Vertex_sign_cov_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	open(2,file=file_name)
	call Bootstrap(real(sign_list,16)/real(n_data,16),n_Block,2,10**3,sign_sd)
	do i = 1,2
		write(1,'(2000000000F0.16)') sqrt(sign_sd(i,i))
		write(2,'(2000000000F30.16)') sign_sd(i,:)
	end do
	close(1)
	close(2)
	
	!Drucke aus, wie häufig die Isospins vorgekommen sind
	file_name = "163Vertex_isospin_number_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	open(1,file=file_name)
	do i = 1,2*N1+1
		write(1,'(2000000000I9)') isospin_number(i)
	end do
	close(1)
	
	!Bestimme den Fehler von isospin_number
	print*,"Berechne Fehler Isospin..."
	file_name = "163Vertex_isospin_sd_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	open(1,file=file_name)
	file_name = "163Vertex_isospin_cov_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	open(2,file=file_name)
	call Bootstrap(real(isospin_list,16)/real(n_data,16),n_Block,2*N1+1,10**3,isospin_sd)
	do i = 1,2*N1+1
		write(1,'(2000000000F0.16)') sqrt(isospin_sd(i,i))
		write(2,'(2000000000F30.16)') isospin_sd(i,:)
	end do
	close(1)
	close(2)
	
	!Drucke aus, wie häufig die Isospins im abs vorgekommen sind
	file_name = "163Vertex_isospin_number_abs_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	open(1,file=file_name)
	do i = 1,2*N1+1
		write(1,'(2000000000I9)') isospin_number_abs(i)
	end do
	close(1)
	
	!Bestimme den Fehler von isospin_number_abs
	print*,"Berechne Fehler Isospin abs..."
	file_name = "163Vertex_isospin_sd_abs_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	open(1,file=file_name)
	file_name = "163Vertex_isospin_cov_abs_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	open(2,file=file_name)
	call Bootstrap(real(isospin_list_abs,16)/real(n_data,16),n_Block,2*N1+1,10**3,isospin_sd_abs)
	do i = 1,2*N1+1
		write(1,'(2000000000F0.16)') sqrt(isospin_sd_abs(i,i))
		write(2,'(2000000000F30.16)') isospin_sd_abs(i,:)
	end do
	close(1)
	close(2)
	
	!Drucke aus, wie häufig die dt Abstände vorgekommen sind
	do i = 1,5
		write(file_i,"(i3)") i
		file_name = "163Vertex_dt_number_Meson"// trim(adjustl(file_i)) // "_M=" // trim(adjustl(file_M)) //&
		",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
		open(1,file=file_name)
		do j = 1,2*N0-1
			write(1,'(2000000000I16)') dt_number(i,j)
		end do
		close(1)
	end do
	
	!Bestimme den Fehler der dt Abstände für alle fünf Mesonen
	print*,"Berechne Fehler dt..."
	file_name = "163Vertex_dt_cov_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	open(1,file=file_name)
	call Bootstrap(real(dt_list(:,:),16)/real(n_data,16),n_Block,5*(2*N0-1),10**3,dt_sd(:,:))
	do j = 1,5*(2*N0-1)
		write(1,'(2000000000F30.16)') dt_sd(j,:)
	end do
	close(1)
	
	!Drucke aus, wie häufig die dt Abstände im abs vorgekommen sind
	do i = 1,5
		write(file_i,"(i3)") i
		file_name = "163Vertex_dt_number_abs_Meson"// trim(adjustl(file_i)) // "_M=" // trim(adjustl(file_M)) //&
		",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
		open(1,file=file_name)
		do j = 1,2*N0-1
			write(1,'(2000000000I16)') dt_number_abs(i,j)
		end do
		close(1)
	end do
	
	!Bestimme den Fehler der dt abs Abstände für alle fünf Mesonen
	print*,"Berechne Fehler dt abs..."
	file_name = "163Vertex_dt_cov_abs_M=" // trim(adjustl(file_M)) //&
	",N="// trim(adjustl(file_N0)) // "x"// trim(adjustl(file_N1)) //".dat"
	open(1,file=file_name)
	call Bootstrap(real(dt_list_abs(:,:),16)/real(n_data,16),n_Block,5*(2*N0-1),10**3,dt_sd_abs(:,:))
	do j = 1,5*(2*N0-1)
		write(1,'(2000000000F30.16)') dt_sd_abs(j,:)
	end do
	close(1)
	
end subroutine Auswertung

subroutine Wurm(N0,N1,M,lattice,sign_lattice,hope_mode,dt_list_Wurm,dt_list_Wurm_abs) !Bond-Typ wird zufällig vorgeschlagen
	
	implicit none
	integer :: N0,N1,hope_mode,randomint,check_end,check_end_prop,sign_lattice,check_valid,check_valid_old
	integer :: i,j,k,dt
	integer :: x0,x1,y0,y1,start0,start1,direction_x,direction_y,old_bond,new_bond_x,new_bond_y
	integer, dimension(4) :: oldvertex_x,oldvertex_y,newvertex_x,newvertex_y
	integer, dimension(N0,N1,4) :: lattice
	integer(kind=8), dimension(2*N0-1) :: dt_list_Wurm,dt_list_Wurm_abs
	real :: M,weight,weightx,weighty,weightx_new,weighty_new
	real :: pA,p0,q1,q2,proposal,u
	
	!Setze die dt_liste auf 0
	do i = 1,2*N0-1
		dt_list_Wurm(i) = 0_8
		dt_list_Wurm_abs(i) = 0_8
	end do
	
	!Schlage einen neuen Startpunkt und Hoptype vor
	start0 = randomint(1,N0)
	start1= randomint(1,N1)
	hope_mode = randomint(1,5)
	
	!Wahrscheinlichkeit, eine Quelle/Senke wegzunehmen
	p0 = 1.0/real(N0*N1*5)
	
	!Bestimme, ob man die Startkonfiguration akzeptiert
	check_end = 1
	q1 = 1.0/real(N0*N1*5)
	q2 = p0
	weightx = weight(M,lattice(start0,start1,:),hope_mode,0,0,1)
	weightx_new = weight(M,lattice(start0,start1,:),hope_mode,1,1,0)
	pA = abs(weightx_new/weightx)*(q2/q1)
	
	call random_number(u)
	if (pA > u) then !Falls man die Startkonfiguration akzeptiert
		check_end = 0
		x0 = start0
		x1 = start1
		
		!Vorzeichen der neuen Konfiguration
		if ((weightx_new/weightx) < 0.0) then
			sign_lattice = -1*sign_lattice
		end if
		
		!Teste, ob die neue Konfiguration valide ist
		if ( abs(weight(M,lattice(start0,start1,:),hope_mode,1,1,1)) > epsilon(0.0) ) then
			check_valid_old = 1
			dt_list_Wurm(N0) = sign_lattice
			dt_list_Wurm_abs(N0) = 1_8
		else
			check_valid_old = 0
		end if
	end if
	
	Wurmloop : do while(check_end==0) !Führe den Wurm weiter, bis man wieder zu einer geschlossenen wurmlosen Konfiguration zurückkehrt.
		
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
			end if
		end if
		
		!Falls man nicht vorschlägt aufzuhören, schlage vor zum Punkt y weiter zu gehen
		if (check_end_prop==0) then
		
			!Schlage eine neue zufällige Richtung vor und bestimme den Punkt y
			direction_x = randomint(1,4)
			select case(direction_x)
			case(1)
				direction_y = 3
				!Definiere den neuen Punkt (y0,y1)
				if (x0 < N0) then
					y0 = x0 + 1
				else !Randbedingung
					y0 = 1
				end if
				y1 = x1
				
			case(2)
				direction_y = 4
				!Definiere den neuen Punkt (y0,y1)
				y0 = x0
				if (x1 < N1) then
					y1 = x1 + 1
				else !Randbedingung
					y1 = 1
				end if
				
			case(3)
				direction_y = 1
				!Definiere den neuen Punkt (y0,y1)
				if (x0 > 1) then
					y0 = x0 - 1
				else !Randbedingung
					y0 = N0
				end if
				y1 = x1
				
			case(4)
				direction_y = 2
				!Definiere den neuen Punkt (y0,y1)
				y0 = x0
				if (x1 > 1) then
					y1 = x1 - 1
				else !Randbedingung
					y1 = N1
				end if
				
			end select
			oldvertex_y = lattice(y0,y1,:)
			
			!Wähle einen zufälligen neuen Bond
			new_bond_x = randomint(0,5)
			if (new_bond_x == 3) then
				new_bond_y = 4
			else if(new_bond_x == 4) then
				new_bond_y = 3
			else
				new_bond_y = new_bond_x
			end if
			
			!Bestimme die neuen Vertex
			newvertex_x = oldvertex_x
			newvertex_x(direction_x) = new_bond_x
			newvertex_y = oldvertex_y
			newvertex_y(direction_y) = new_bond_y
			
		end if
		
		!Bestimme die Gewichte bei x
		if (y0==0 .and. y1==0) then !Falls man vorschlägt aufzuhören, nimm den Wurmkopf weg
			weightx = weight(M,oldvertex_x,hope_mode,1,1,0)
			weightx_new = weight(M,newvertex_x,hope_mode,0,0,1)
		else if (x0==start0 .and. x1==start1) then !Falls man am Start ist
			weightx = weight(M,oldvertex_x,hope_mode,1,1,0)
			weightx_new = weight(M,newvertex_x,hope_mode,0,1,1)
		else !Sonst bewege den Wurmkopf zu y
			weightx = weight(M,oldvertex_x,hope_mode,1,0,0)
			weightx_new = weight(M,newvertex_x,hope_mode,0,0,1)
		end if
		
		!Bestimme die Gewichte bei y
		if (y0==0 .and. y1==0) then !Falls man aufhört, ist y=0.
			weighty = 1.0
			weighty_new = 1.0
			if ( abs(weightx_new)>0.0 ) then
				check_valid = 1
			else
				check_valid = 0
			end if
		else if (y0==start0 .and. y1==start1) then !Falls man zum Start geht, schliesst man die config.
			weighty = weight(M,oldvertex_y,hope_mode,0,1,1)
			weighty_new = weight(M,newvertex_y,hope_mode,1,1,0)
			!Bestimme, ob die neue Konfiguration valide ist
			if ( abs(weightx_new)>0.0 .and. abs(weight(M,newvertex_y,hope_mode,1,1,1))>0.0 ) then
				check_valid = 1
			else
				check_valid = 0
			end if
		else !Sonst bewegt sich der Wurmkopf weiter zu y
			weighty = weight(M,oldvertex_y,hope_mode,0,0,1)
			weighty_new = weight(M,newvertex_y,hope_mode,1,0,0)
			!Bestimme, ob die neue Konfiguration valide ist
			if ( abs(weightx_new)>0.0 .and. abs(weight(M,newvertex_y,hope_mode,1,0,1))>0.0 ) then
				check_valid = 1
			else
				check_valid = 0
			end if
		end if
		
		!Bestimme die Proposal-Wahrscheinlichkeiten, von x->y zu gehen, sowie die Akzeptanzwahrscheinlichkeit
		q1 = proposal(N0,N1,x0,x1,y0,y1,start0,start1,p0)
		q2 = proposal(N0,N1,y0,y1,x0,x1,start0,start1,p0)
		pA = abs(weightx_new/weightx) * abs(weighty_new/weighty) * (q2/q1)
		
		call random_number(u)
		if (pA > u) then !Falls man die neue Konfiguration akzeptiert
			
			check_valid_old = check_valid
			!Vorzeichen der neuen Konfiguration
			if ((weightx_new/weightx) * (weighty_new/weighty) < 0.0) then
				sign_lattice = -1*sign_lattice
			end if
			
			if (y0==0 .and. y1==0) then !Falls man aufhört
				check_end = 1
				exit Wurmloop
			else
				lattice(x0,x1,:) = newvertex_x
				lattice(y0,y1,:) = newvertex_y
				x0 = y0
				x1 = y1
			end if
		end if
		
		!Führe die dt-Liste weiter, falls man eine valide Konfiguration hat
		if (check_valid_old==1) then
			dt=x0-start0
			dt_list_Wurm(dt+N0) = dt_list_Wurm(dt+N0)+sign_lattice
			dt_list_Wurm_abs(dt+N0) = dt_list_Wurm_abs(dt+N0)+1_8			
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

function proposal(N0,N1,x0,x1,y0,y1,start0,start1,p0) !Gibt die Wahrscheinlichkeit, eine Konfiguration vorgeschlagen zu haben.

	implicit none
	integer :: N0,N1,x0,x1,y0,y1,start0,start1
	real :: proposal,p0
	
	if (x0==0 .and. x1==0) then !Falls man auf eine wurmlose config. einen neuen Wurm setzt
		proposal = 1.0 / real(N0*N1*5)
	else if (y0==0 .and. y1==0) then !Falls man in die wurmlose config. zurückwechselt
		proposal = p0
	else if (x0==start0 .and. x1==start1) then !Falls man vom Start aus den Wurm fortsetzt
		proposal = (1.0-p0)/(4.0*6.0)
	else !Default: Wähle eine zufällige Richtung und Bondtyp
		proposal = 1.0/(4.0*6.0)
	end if
	
end function proposal

function weight(M,mu,hope_mode,check_x,check_start,onlyreal) !Gibt das Gewicht des Vertex mu. hope_mode: Typ der Quelle/Senke. check_x=0: Keine Quelle/Senke. check_x=1: Senke. check_start=1: Quelle. onlyreal: Ob "verbotene Vertex" zulässig sind.

	implicit none
	integer :: hope_mode,check_x,check_start,onlyreal,direction,bond,n_valid
	integer, dimension(4) :: mu,mu_test
	real :: M,weight,test_weight,weight_vertex,weight_163Vertex
	
	test_weight = weight_163Vertex(M,mu,hope_mode,check_x,check_start)

	if (onlyreal==0 .and. abs(test_weight)<epsilon(0.0)) then !Falls das Gewicht 0 ist und "verbotene Vertex" zulässig sind
		!Teste, ob der Vertex nach einem Updateschritt valide wird
		test_weight = 1.0
		n_valid = 0
		
		do direction = 1,4
			do bond = 0,5
				
				mu_test = mu
				mu_test(direction) = bond
				if (check_x==1 .and. check_start==0) then !Wurmkopf wird wegbewegt
					weight_vertex = abs(weight_163Vertex(M,mu_test,hope_mode,0,0))
				else if (check_x==1 .and. check_start==1) then !Wurmkopf wird vom Start aus wegbewegt
					weight_vertex = abs(weight_163Vertex(M,mu_test,hope_mode,0,1))
				else
					weight_vertex = 0.0
				end if
				
				if ( weight_vertex > epsilon(0.0) ) then !Falls der Vertex, sobald man die Quelle bewegt, valide wird
					test_weight = test_weight * weight_vertex
					n_valid = n_valid + 1
				end if
				
			end do
		end do
		
		if (n_valid == 0) then
			test_weight = 0.0
		else
			test_weight = test_weight**(1.0 / real(n_valid))
		end if
	end if
	
	weight = test_weight
	
end function weight

function weight_163Vertex(M,mu,hope_mode,check_x,check_start) !Gibt das 163-Vertex Gewicht.

	implicit none
	real:: M,weight_163Vertex,weight_blue,weight_orange,weight_FreeDirac
	integer :: i,hope_mode,check_x,check_start,n_sources_blue,n_sinks_blue,n_sources_orange,n_sinks_orange
	integer, dimension(4) :: mu,blue_out,blue_in,orange_out,orange_in
	
	!Bestimme die Position der blauen und orangen ausgehenden und eingehenden Hops
	do i = 1,4
		select case(mu(i))
		
			case(1)
			blue_out(i) = 1
			blue_in(i) = 1
			orange_out(i) = 0
			orange_in(i) = 0
			
			case(2)
			blue_out(i) = 0
			blue_in(i) = 0
			orange_out(i) = 1
			orange_in(i) = 1
			
			case(3)
			blue_out(i) = 1
			blue_in(i) = 0
			orange_out(i) = 0
			orange_in(i) = 1
			
			case(4)
			blue_out(i) = 0
			blue_in(i) = 1
			orange_out(i) = 1
			orange_in(i) = 0
			
			case(5)
			blue_out(i) = 1
			blue_in(i) = 1
			orange_out(i) = 1
			orange_in(i) = 1
			
			case default
			blue_out(i) = 0
			blue_in(i) = 0
			orange_out(i) = 0
			orange_in(i) = 0
		
		end select
	end do
	
	!Bestimme, ob und wie viele Quellen und Senken vorkommen
	if (check_x==0 .and. check_start==0) then !Keine Quellen oder Senken
		n_sources_blue = 0
		n_sinks_blue = 0
		n_sources_orange = 0
		n_sinks_orange = 0

	else if (check_x==1 .and. check_start==0) then !Wurmkopf
	
		select case(hope_mode)
		case(1)
			n_sources_blue = 1
			n_sinks_blue = 1
			n_sources_orange = 0
			n_sinks_orange = 0
		case(2)
			n_sources_blue = 0
			n_sinks_blue = 0
			n_sources_orange = 1
			n_sinks_orange = 1
		case(3)
			n_sources_blue = 0
			n_sinks_blue = 1
			n_sources_orange = 1
			n_sinks_orange = 0
		case(4)
			n_sources_blue = 1
			n_sinks_blue = 0
			n_sources_orange = 0
			n_sinks_orange = 1
		case(5)
			n_sources_blue = 1
			n_sinks_blue = 1
			n_sources_orange = 1
			n_sinks_orange = 1
		end select
	
	else if (check_x==0 .and. check_start==1) then !Start
		
		select case(hope_mode)
		case(1)
			n_sources_blue = 1
			n_sinks_blue = 1
			n_sources_orange = 0
			n_sinks_orange = 0
		case(2)
			n_sources_blue = 0
			n_sinks_blue = 0
			n_sources_orange = 1
			n_sinks_orange = 1
		case(3)
			n_sources_blue = 1
			n_sinks_blue = 0
			n_sources_orange = 0
			n_sinks_orange = 1
		case(4)
			n_sources_blue = 0
			n_sinks_blue = 1
			n_sources_orange = 1
			n_sinks_orange = 0
		case(5)
			n_sources_blue = 1
			n_sinks_blue = 1
			n_sources_orange = 1
			n_sinks_orange = 1
		end select
		
	else if (check_x==1 .and. check_start==1) then !Wurmkopf ist am Start
	
		select case(hope_mode)
		case(1)
			n_sources_blue = 2
			n_sinks_blue = 2
			n_sources_orange = 0
			n_sinks_orange = 0
		case(2)
			n_sources_blue = 0
			n_sinks_blue = 0
			n_sources_orange = 2
			n_sinks_orange = 2
		case(3)
			n_sources_blue = 1
			n_sinks_blue = 1
			n_sources_orange = 1
			n_sinks_orange = 1
		case(4)
			n_sources_blue = 1
			n_sinks_blue = 1
			n_sources_orange = 1
			n_sinks_orange = 1
		case(5)
			n_sources_blue = 2
			n_sinks_blue = 2
			n_sources_orange = 2
			n_sinks_orange = 2
		end select
	
	end if
	
	weight_blue = weight_FreeDirac(M,blue_out,blue_in,n_sources_blue,n_sinks_blue)
	weight_orange = weight_FreeDirac(M,orange_out,orange_in,n_sources_orange,n_sinks_orange)
	
	weight_163Vertex = weight_blue * weight_orange
	
end function weight_163Vertex

function weight_FreeDirac(M,mu_out,mu_in,n_sources,n_sinks) !Gibt das 49-Vertex Gewicht.

	implicit none
	
	real:: weight_FreeDirac,M,weight_out,weight_in
	integer :: i, summe_out,summe_in,n_sources,n_sinks
	integer, dimension(4) :: mu_out,mu_in
	
	summe_out = 0
	summe_in = 0
	
	do i = 1,4
		summe_out = summe_out + mu_out(i)
		summe_in = summe_in + mu_in(i)
	end do
	
	if (n_sinks==0 .and. n_sources==0) then !Keine Quellen/Senken
		!Leere Konfiguration
		if ( summe_out == 0 .and. summe_in == 0 ) then
			weight_FreeDirac = M**2
		
		else if( summe_out == 1 .and. summe_in == 1 ) then
			if ( all( mu_out == mu_in ) ) then
				weight_FreeDirac = 0.0
			!Gerade einfache Hops
			else if ( (mu_out(1)==1 .and. mu_in(3)==1) .or. &
			(mu_out(3)==1 .and. mu_in(1)==1) .or. &
			(mu_out(2)==1 .and. mu_in(4)==1) .or. &
			(mu_out(4)==1 .and. mu_in(2)==1) ) then
				weight_FreeDirac = M
			else
				weight_FreeDirac = M / sqrt(2.0)
			end if
		
		else if( summe_out == 2 .and. summe_in == 2 ) then
			!Doppelpfeil-Terme
			if ( all( mu_out == mu_in ) ) then
				if ( mu_out(1) == mu_out(3) .or. mu_out(2) == mu_out(4) ) then
					weight_FreeDirac = 1.0
				else
					weight_FreeDirac = 0.5
				end if
			!Kreuzterme
			else if ( all( mu_out /= mu_in ) ) then
				if ( mu_out(1) == mu_out(3) .or. mu_out(2) == mu_out(4) ) then
					weight_FreeDirac = 1.0
				else
					weight_FreeDirac = -0.5
				end if

			!3-point terme
			else if ( ( all( mu_out == (/0,1,1,0/) ) .and. all( mu_in == (/0,0,1,1/) ) ) .or. &
			( all( mu_out == (/1,0,0,1/) ) .and. all( mu_in == (/1,1,0,0/) ) ) .or. &
			( all( mu_out == (/0,0,1,1/) ) .and. all( mu_in == (/1,0,0,1/) ) ) .or. &
			( all( mu_out == (/1,1,0,0/) ) .and. all( mu_in == (/0,1,1,0/) ) )  ) then
				weight_FreeDirac = 0.5
				
			else if ( ( all( mu_out == (/1,1,0,0/) ) .and. all( mu_in == (/1,0,0,1/) ) ) .or. &
			( all( mu_out == (/0,0,1,1/) ) .and. all( mu_in == (/0,1,1,0/) ) ) .or. &
			( all( mu_out == (/0,1,1,0/) ) .and. all( mu_in == (/1,1,0,0/) ) ) .or. &
			( all( mu_out == (/1,0,0,1/) ) .and. all( mu_in == (/0,0,1,1/) ) )  ) then
				weight_FreeDirac = -0.5
			
			else if ( ( all( mu_out == (/0,1,1,0/) ) .and. all( mu_in == (/1,0,1,0/) ) ) .or. &
			( all( mu_out == (/1,0,1,0/) ) .and. all( mu_in == (/1,1,0,0/) ) ) .or. &
			( all( mu_out == (/0,1,0,1/) ) .and. all( mu_in == (/1,0,0,1/) ) ) .or. &
			( all( mu_out == (/1,1,0,0/) ) .and. all( mu_in == (/0,1,0,1/) ) ) .or. &
			( all( mu_out == (/1,0,1,0/) ) .and. all( mu_in == (/0,0,1,1/) ) ) .or. &
			( all( mu_out == (/1,0,0,1/) ) .and. all( mu_in == (/1,0,1,0/) ) ) .or. &
			( all( mu_out == (/0,0,1,1/) ) .and. all( mu_in == (/0,1,0,1/) ) ) .or. &
			( all( mu_out == (/0,1,0,1/) ) .and. all( mu_in == (/0,1,1,0/) ) ) ) then
				weight_FreeDirac = 1.0/sqrt(2.0)
			
			else if ( ( all( mu_out == (/1,1,0,0/) ) .and. all( mu_in == (/1,0,1,0/) ) ) .or. &
			( all( mu_out == (/1,0,1,0/) ) .and. all( mu_in == (/0,1,1,0/) ) ) .or. &
			( all( mu_out == (/0,1,0,1/) ) .and. all( mu_in == (/1,1,0,0/) ) ) .or. &
			( all( mu_out == (/1,0,0,1/) ) .and. all( mu_in == (/0,1,0,1/) ) ) .or. &
			( all( mu_out == (/1,0,1,0/) ) .and. all( mu_in == (/1,0,0,1/) ) ) .or. &
			( all( mu_out == (/0,0,1,1/) ) .and. all( mu_in == (/1,0,1,0/) ) ) .or. &
			( all( mu_out == (/0,1,1,0/) ) .and. all( mu_in == (/0,1,0,1/) ) ) .or. &
			( all( mu_out == (/0,1,0,1/) ) .and. all( mu_in == (/0,0,1,1/) ) ) ) then
				weight_FreeDirac = -1.0/sqrt(2.0)
			
			end if
		else
			weight_FreeDirac = 0.0
		end if
		
	else if (n_sources==1 .and. n_sinks==0) then !psi(x)
		
		if (summe_out==1 .and.summe_in==0) then !psi(x) mit Massenterm
			weight_FreeDirac = -M
		else if (summe_out==2 .and. summe_in==1) then !psi(x) mit zwei Hoptermen
			if (mu_out(1)==1 .and. mu_out(3)==1) then
				weight_FreeDirac = 1.0
			else if(mu_out(2)==1 .and. mu_out(4)==1) then
				weight_FreeDirac = 1.0
			else
				weight_FreeDirac = 1.0/sqrt(2.0)
			end if
		else
			weight_FreeDirac = 0.0
		end if
	
	else if (n_sources==0 .and. n_sinks==1) then !bar{psi}(x)
		
		if (summe_out==0 .and. summe_in==1) then !bar{psi}(x) mit Massenterm
			weight_FreeDirac = -M
		else if (summe_out==1 .and. summe_in==2) then !bar{psi}(x) mit zwei Hoptermen
			if (mu_in(1)==1 .and. mu_in(3)==1) then
				weight_FreeDirac = 1.0
			else if(mu_in(2)==1 .and. mu_in(4)==1) then
				weight_FreeDirac = 1.0
			else !Ecke
				weight_FreeDirac = 1.0/sqrt(2.0)
			end if	
		else
			weight_FreeDirac = 0.0
		end if
		
	else if (n_sinks==1 .and. n_sources==1) then !bar{psi}(x) psi(x)
		
		if (summe_out == 0 .and. summe_in == 0) then !Quelle/Senke am Start
			weight_FreeDirac = -2.0*M
		else if (summe_out == 1 .and. summe_in == 1) then !Doppel Quelle/Senke
			if ( all( mu_out == mu_in ) ) then !Quelle/Senke für Wurmtyp 1/2
				weight_FreeDirac = -1.0
			else if ( (mu_out(1)==1 .and. mu_in(3)==1) .or. &
			(mu_out(3)==1 .and. mu_in(1)==1) .or. &
			(mu_out(2)==1 .and. mu_in(4)==1) .or. &
			(mu_out(4)==1 .and. mu_in(2)==1) ) then
				weight_FreeDirac = -1.0
			else
				weight_FreeDirac = -1.0 / sqrt(2.0)
			end if
		else
			weight_FreeDirac = 0.0
		end if
		
	else if (n_sinks==2 .and. n_sources==2) then !bar{psi}(x) psi(x) bar{psi}(x) psi(x)
		
		if (summe_out == 0 .and. summe_in == 0) then !Quelle/Senke am Start für Wurmtyp 1/2
			weight_FreeDirac = 4.0
		else
			weight_FreeDirac = 0.0
		end if
		
	else
		weight_FreeDirac = 0.0
	
	end if

end function weight_FreeDirac

function isospin(lattice,N0,N1) !Gibt den Isospin einer Konfiguration zurück

	implicit none
	integer :: isospin,N0,N1,i,test
	integer, dimension(N0,N1,4) :: lattice
	
	test = 0

	!Isospin über den rechten Rand (Nur die Zeit zählt)
	do i = 1,N1
	
		if ( lattice(N0,i,1)==3 ) then
			test = test + 1
		else if ( lattice(N0,i,1)==4 ) then
			test = test - 1
		end if
		
	end do

	isospin = test

end function isospin

function randomint(a,b) !Zufälliges Integer zwischen a und b (inklusive a und b)

	implicit none
	integer :: randomint,a,b
	real :: u
	
	call random_number(u)
	randomint = a + floor(real(b-a+1)*u)
	
end function randomint
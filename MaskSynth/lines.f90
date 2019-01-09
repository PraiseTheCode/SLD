subroutine lines ! Calculate individual lines profiles
use dfport
use synth_data
use lines_data
use chemical_elements
use radiation_field, only : mu_fix,Inu,Hnu
use mol_dat, only : Mol_id, Nmol
integer integer_wavelength,integer_E ! Data in binary file are placed as a integer(4) values
integer(2) integer_gf,integer_code,integer_dumping1,integer_dumping2,integer_dumping3 ! Data in binary file are placed as a integer(2) values
real(8) E,limb_darkening,a,b
character(5) new_VALD
character(100)string
rewind(3)


  if(atomic_lines == 1)then ! Take atomic lines data from binary file
  nrec=first_record_atomic-1; mode=1
  do
   read(3,rec=nrec)integer_wavelength,integer_gf,integer_E,integer_code,integer_dumping1,integer_dumping2,integer_dumping3
   central_wavelength=integer_wavelength/1000.d0 ! (A)
   if(central_wavelength > wavelength_last+2.5d0)exit ! Complete lines selection
   nrec=nrec+1
   gf=10**(integer_gf/1000.d0); if(gf > 5)gf=-gf ! For prevent possible bugs in source VALD list
   code=integer_code/100.d0; E=integer_E/1000.d0; Elow=E*2.997925d10 !Convert integers to float (true) values
   gamma_Rad=10**(integer_dumping1/100.d0) ;gamma_Stark=10**(integer_dumping2/100.d0) ;gamma_WdW=10**(integer_dumping3/100.d0)
   atomic_number=code+0.5; ionization=integer_code-atomic_number*100+1
   write(line_identifier,'(a2,i2)')Elem_id(atomic_number),ionization
   if(ionization > 6)cycle ! Maximum 5th ions
   call calculate_profile(mode) 
  enddo
 endif

  if(atomic_lines == 2)then ! Take atomic lines data from VALD file
   read(3,*,iostat=ios)wstart,wmax,Nlines
   if(ios /= 0)stop 'Wrong data in VALD file'
   read(3,*) ! Skip head's records
   read(3,*)
   do nl=1,Nlines ! Read lines data
    mode=1; new_VALD=' '; VALD_string(1:200)=' '
	read(3,'((a))',iostat=ios)VALD_string
	if(ios /= 0)stop 'Number of lines in VALD file is incorrect'
	read(VALD_string,*,iostat=ios)new_VALD,central_wavelength,E,dummy,gf,gamma_Rad,gamma_Stark,gamma_WdW
	if(ios /= 0)then
	 write(*,*)'Wrong record number ',nl,' in VALD file'; stop
    endif
	if(new_VALD == 'H  1')cycle  ! Skip Hydrogen lines
	if(central_wavelength < wavelength_first-2.5d0)cycle ! skip lines outside wavelengths region
    if(central_wavelength > wavelength_last+2.5d0)exit
    Elow=E*2.997925d10*8066.d0
    gf=10**gf
	line_identifier=new_VALD(1:4)
	if(molecular_lines == 2)then
     do l=1,Nmol
	  k=index(new_VALD,trim(Mol_id(l)))
	  if(k /= 0)then
       mode = 2; exit
      endif
     enddo
	endif
	if(mode == 2)cycle ! Skip molecular lines
!	 if(new_VALD(3:3) == ' '.and.new_VALD(4:4) == '1')line_identifier(4:4)=' '  ! CO, C2 etc.!
!	 if(new_VALD(3:3) == ' '.and.new_VALD(4:4) == '2')line_identifier(3:4)='+ ' ! CO+,CH+ etc.
!	 if(new_VALD(3:3) /= ' '.and.new_VALD(5:5) == '2')line_identifier(4:4)='+'  ! MgH+ etc.
!	 if(gamma_Rad == 0.d0)then
!	  gamma_Rad=2.223d13/(central_wavelength**2)
!     else
!	  gamma_Rad=10**gamma_Rad
!     endif
!     gamma_Stark=1.d-5   ! Stark and Wan der Waals constant fixed for molecular lines
!     gamma_WdW=1.d-7
!	else  
     read(line_identifier(3:4),*)ionization
	 if(ionization > 6)cycle ! Maximum 5th ions
     do i=1,99  ! Search the atomic number
      if(line_identifier(1:2) == Elem_id(i))exit
     enddo
     atomic_number=i
	 if(atomic_number == 1)cycle  ! Skip Hydrogen lines  
     code=atomic_number+(ionization-1)/100.d0
  	 call calculate_line_constants(central_wavelength,E,code,gamma_Rad,gamma_Stark,gamma_WdW) ! Calculate dumping constants if in VALD file it's equial zero
     gamma_Rad=10**gamma_Rad; gamma_Stark=10**gamma_Stark; if(gamma_WdW < 10)gamma_WdW=10**gamma_WdW
!	endif
    call calculate_profile(mode) 
   enddo

  endif

  close(223)

end


subroutine get_alphas
	use synth_data
	real(8), allocatable :: vprof(:), rprofs(:,:), alphas(:,:), wl_alphas(:,:), wlines(:), simple_sum(:), waves(:), vels(:), interp_alph(:)
	!real(8), allocatable :: vprof(:), rprofs(:,:), alphas(:,:), wl_alphas(:), wlines(:), simple_sum(:), waves(:), vels(:), interp_alph(:)
	integer n_lines, n_points, n_middle, n_blend, cont_count
	character(4), allocatable :: temp_idents(:)
	real(8) profint, vp, alphaint, contr
	real(8), allocatable :: contribs(:)
	data vc/2.997925d5/

	synth_R = 1.d0 - synth_R

	open(666, file='profiles.bnr',status='old',form='binary',access='direct',recl=12,iostat=ios)
	if(ios /= 0)stop 'cannot open profiles.bnr file'
	read(666, rec=1) n_lines,n_points,n_middle
	allocate(temp_idents(n_lines),vprof(n_points),rprofs(n_points, n_lines), wlines(n_lines), simple_sum(number_of_wavelengths), waves(number_of_wavelengths))
	allocate(alphas(n_points,n_lines),wl_alphas(number_of_wavelengths,n_lines),vels(number_of_wavelengths),interp_alph(n_points))
	close (666)
	open(666,file='profiles.bnr',form='binary',access='direct',recl=n_points*8+12)
	read(666,rec=2) vprof

	do i = 1,n_lines
		read(666,rec=i+2) wlines(i), temp_idents(i), rprofs(1:n_points,i)
		rprofs(1:n_points,i) = 1.d0 - rprofs(1:n_points,i)
	enddo
	close(666)
	
	wl_alphas = 1.d0
	simple_sum = 0.d0
	do i = 1, number_of_wavelengths
		waves(i) = synth_wl(i)
		n_blend = 0
		do j = 1, n_lines
			vp=((wlines(j)-waves(i))/waves(i))*vc 
			if(vp < vprof(1))cycle 
			if(vp > vprof(n_points))exit
			n_blend = n_blend + 1
			id=map1(vprof,rprofs(1:n_points,j),n_points,vp,profint,1)
			simple_sum(i)=simple_sum(i)+profint
		enddo
		cont_count = 0
		if (n_blend > 0) then
			!allocate(contribs(n_blend))
			do j = 1, n_lines
				vp=((wlines(j)-waves(i))/waves(i))*vc 
				if(vp < vprof(1))cycle 
				if(vp > vprof(n_points))exit
				cont_count = cont_count + 1
				id=map1(vprof,rprofs(1:n_points,j),n_points,vp,profint,1)
				!contribs(cont_count) = profint / simple_sum(i) 
				contr = profint / simple_sum(i) 

				wl_alphas(i,j) = (synth_R(i)-simple_sum(i))/simple_sum(i)
				wl_alphas(i,j) = wl_alphas(i,j) + 1.d0
			enddo
			!deallocate(contribs)
		endif
	enddo

	do i = 1, n_lines
		do j = 1, number_of_wavelengths
			vels(j) = ((waves(j)-wlines(i))/wlines(i))*vc 
		enddo
		id = map1(vels,wl_alphas(1:number_of_wavelengths,i),number_of_wavelengths,vprof,interp_alph,n_points)
		alphas(1:n_points,i) = interp_alph
	enddo

	open(666,file='testtest.dat')
	simple_sum = 0.d0
	do i = 1, number_of_wavelengths
		waves(i) = synth_wl(i)
		do j = 1, n_lines
			vp=((waves(i)-wlines(j))/wlines(j))*vc 
			if(vp < vprof(1))cycle 
			if(vp > vprof(n_points))cycle
			id=map1(vprof,rprofs(1:n_points,j),n_points,vp,profint,1)
			id=map1(vprof,alphas(1:n_points,j),n_points,vp,alphaint,1)
			simple_sum(i)=simple_sum(i)+profint*alphaint !wl_alphas(i,j) 
			!print *, alphaint, profint
		enddo
		!write(*,*) (synth_R(i)-simple_sum(i))/synth_R(i)
		write(666,*) synth_wl(i), 1.d0-synth_R(i), 1.d0-simple_sum(i)
		!pause
	enddo
	close(666)
	
	open(666,file='profiles.bnr',form='binary',access='direct',recl=n_points*8*2+12)
	write(666,rec=2)(vprof(i),i=1,n_points)
	write(666,rec=1)n_lines,n_points,n_middle
	do 	nrecord=3,n_lines+2
		write(666,rec=nrecord)wlines(nrecord-2),temp_idents(nrecord-2),(1.d0-rprofs(i,nrecord-2),i=1,n_points),(alphas(i,nrecord-2),i=1,n_points)
		!print *, temp_idents(nrecord-2)
		!write(*,*) wlines(nrecord-2),temp_idents(nrecord-2),(1.d0-rprofs(i,nrecord-2),i=1,n_points),(alphas(i,nrecord-2),i=1,n_points)
		!pause
	enddo
	close(666)

end subroutine get_alphas

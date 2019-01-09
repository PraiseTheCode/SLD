subroutine fun(x,f,jacobian,nx,mode,nspec) !get_model_spectrum
use common_data
use dfport
implicit real(8) (a-h,o-z)
character(4) line_ident,mode
real(8) x(nx),f(nobs),jacobian(nobs,nx)
data c/2.997925d5/

 !SLD_prof=x(2:nx)
 model_spectrum=0.d0; f=0.d0
 if(mode == 'grad')jacobian=0.d0
 Vrad=x(1); dVrad=5.d0

if (nSLD == 1) then
	 do n=1,nlines
	  !read(2,rec=n+2)central_wavelength,line_ident,rprof
	  central_wavelength = wlines(n)
	  rprof = all_rprofs(1:npoints,n)
	  !rprof=1.d0-rprof
	  loc_alphas = alphas(1:npoints,n)
	  SLD_prof = x(2:(npoints+1))
	  ifspecial = 0
	  do j = 1, nspecial
		if (line_idents(n) == special_idents(j)) then
			ifspecial = 1
			exit
		endif
	  enddo
	  !print *, ifspecial
	  !if (line_idents(n) == 'Fe 1' .or. line_idents(n) == 'Fe 2') then
	  if (ifspecial == 1) then
		SLD_prof = x((npoints+2):nx) 
		call convolve
		  do i=1,nobs
		   vp=((wl_obs(i)-central_wavelength)/central_wavelength)*c !  Corresponded for given wavelength velocity
		   if(vp < vprof(1)+Vrad)cycle 
		   if(vp > vprof(npoints)+Vrad)exit
		   id=map1(vprof+Vrad,ySLD,npoints,vp,profint,1)
		   model_spectrum(i)=model_spectrum(i)+profint

		! Jacobian
		  if(mode == 'grad')then
			if(vp < vprof(2)+Vrad)cycle 
			if(vp > vprof(npoints-1)+Vrad)exit 
			do k=1,npoints-1
			 if(vp >= vprof(k)+Vrad.and.vp <= vprof(k+1)+Vrad)then
			  do l=2,npoints+1
			  jacobian(i,l+npoints)=jacobian(i,l+npoints)+loc_grad(k,l-1)+(loc_grad(k+1,l-1)-loc_grad(k,l-1))*(vp-vprof(k)-Vrad)/(vprof(k+1)-vprof(k)-Vrad)
			  enddo
			  exit
			 endif
			enddo
		   endif
		  enddo
		  
	   else
		  call convolve
		  do i=1,nobs
		   vp=((wl_obs(i)-central_wavelength)/central_wavelength)*c !  Corresponded for given wavelength velocity
		   if(vp < vprof(1)+Vrad)cycle 
		   if(vp > vprof(npoints)+Vrad)exit
		   id=map1(vprof+Vrad,ySLD,npoints,vp,profint,1)
		   model_spectrum(i)=model_spectrum(i)+profint

		! Jacobian
		  if(mode == 'grad')then
			if(vp < vprof(2)+Vrad)cycle 
			if(vp > vprof(npoints-1)+Vrad)exit 
			do k=1,npoints-1
			 if(vp >= vprof(k)+Vrad.and.vp <= vprof(k+1)+Vrad)then
			  do l=2,npoints+1
			  jacobian(i,l)=jacobian(i,l)+loc_grad(k,l-1)+(loc_grad(k+1,l-1)-loc_grad(k,l-1))*(vp-vprof(k)-Vrad)/(vprof(k+1)-vprof(k)-Vrad)
			  enddo
			  exit
			 endif
			enddo
		   endif
		  enddo
		endif

	 enddo
else
	!print *, 1
	do n=1,nlines
	  !read(2,rec=n+2)central_wavelength,line_ident,rprof
	  central_wavelength = wlines(n)
	  rprof = all_rprofs(1:npoints,n)
	  !rprof=1.d0-rprof
	  loc_alphas = alphas(1:npoints,n)
	  !SLD_prof = x(2:(npoints+1))
	  ifspecial = 0
	  do j = 1, nspecial
		if (line_idents(n) == special_idents(j)) then
			ifspecial = 1
			exit
		endif
	  enddo
	!print *, 2
	  call interp_SLD(ifspecial+1, x, nx, rprof(nmiddle+1), kk)
	!print *, 3
	  SLD_prof = flux_sld
	!print *, 4
	  call convolve
	!print *, 5
	  do i=1,nobs
		   vp=((wl_obs(i)-central_wavelength)/central_wavelength)*c !  Corresponded for given wavelength velocity
		   if(vp < vprof(1)+Vrad)cycle 
		   if(vp > vprof(npoints)+Vrad)exit
		   id=map1(vprof+Vrad,ySLD,npoints,vp,profint,1)
		   model_spectrum(i)=model_spectrum(i)+profint

		! Jacobian
		  if(mode == 'grad')then
			if(vp < vprof(2)+Vrad)cycle 
			if(vp > vprof(npoints-1)+Vrad)exit 
			do k=1,npoints-1
			 if(vp >= vprof(k)+Vrad.and.vp <= vprof(k+1)+Vrad)then
			  do l=2,npoints+1
			  jacobian(i,l+ifspecial*npoints*nSLD + npoints*(kk-1))=jacobian(i,l+ifspecial*npoints*nSLD + npoints*(kk-1))+loc_grad(k,l-1)+(loc_grad(k+1,l-1)-loc_grad(k,l-1))*(vp-vprof(k)-Vrad)/(vprof(k+1)-vprof(k)-Vrad)
			  enddo
			  exit
			 endif
			enddo
		   endif
	   enddo
	 enddo
endif

if(mode == 'grad')then ! df/dx(1) calculated numerically by finite differencies
 do n=1,nlines
  central_wavelength = wlines(n)
  rprof = all_rprofs(1:npoints,n)
  loc_alphas = alphas(1:npoints,n)
  !read(2,rec=n+2)central_wavelength,line_ident,rprof
  !rprof=1.d0-rprof
  ySLD = x(2:npoints)
  call convolve
  do i=1,nobs
   vp=((wl_obs(i)-central_wavelength)/central_wavelength)*c
     if(vp < vprof(1)+Vrad+dVrad)cycle 
   if(vp > vprof(npoints)+Vrad+dVrad)exit 
   id=map1(vprof+Vrad+dVrad,ySLD,npoints,vp,profint,1)
   f(i)=f(i)+profint  
  enddo
 enddo
 jacobian(1:nobs,1)=(model_spectrum-f)/dVrad
endif

 f=model_spectrum

! print *, nobs, size(jacobian(1,:)), npoints
! print *, jacobian

end


subroutine interp_SLD(icomp,x,nx,flux_line,kk)	! interpolate SLD profile for a given line strength
	use common_data
	implicit real(8) (a-h,o-z)
	real(8) x(nx)

	lsd_temp = 0.d0; flux_sld = 0.d0
	do iv = 1, npoints
		do ilsd = 1, nSLD
			lsd_temp(iv,ilsd) = x(1+(npoints*nSLD*(icomp-1)+npoints*(ilsd-1))+iv)
		enddo
	enddo

	do iv = 1, npoints
		ios = map1(LS_mean,lsd_temp(iv,1:nSLD),nSLD,flux_line,flux_sld(iv),1)
	enddo

	if(flux_line <= LS_limits(1)) then
		kk = 1
	elseif(flux_line > LS_limits(nSLD-1)) then
		kk = nSLD
	else
		do ilimit = 1, nSLD-2
			if(flux_line > LS_limits(ilimit) .and. flux_line <= LS_limits(ilimit+1)) then
				kk = ilimit + 1
				exit
			endif
		enddo
	endif

end subroutine interp_SLD

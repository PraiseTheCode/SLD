subroutine SLD_profile ! Search convolution profile by adjust model spectrum to obsrved
use common_data
implicit real(8) (a-h,o-z)
real(8), allocatable :: x(:),weight(:)
real(8) gauss_pars(3)
external fun, gauss_fun_finite

nns = 1
if(nspecial > 0) nns = 2

if(.not.allocated(x))&
  allocate (x(nns*nSLD*npoints+1),weight(nobs),stat=ios)
   if(ios /= 0)stop 'Memory allocation failed'
weight=1.d0

Vrad = 0.d0
if (cross_corr) call get_Vrad
print *, Vrad

gauss_pars(1) = 0.01d0
gauss_pars(2) = Vrad
gauss_pars(3) = 80.d0

call mini(gauss_pars,r_obs,weight,3,nobs,1,gauss_fun_finite)
print *, gauss_pars
 
 weight=1.d0; x=0.d0; 
 x(1)=gauss_pars(2)
 gauss_pars(2) = 0.d0
 !x(1) = Vrad
! x(nmiddle+1)=1.d0
 do i = 1, nns
	do j = 1, nSLD
		do k = 1, npoints
			x(1+(i-1)*nSLD*npoints + (j-1)*npoints + k) = G(vprof(k)+Vrad,gauss_pars)
		enddo
	enddo
 enddo
!model_spectrum = 0.d0
!do j = 1,nlines
!  do i=1,nobs
!   vp=((wl_obs(i)-wlines(j))/wlines(j))*2.997925d5 !  Corresponded for given wavelength velocity
!   if(vp < vprof(1))cycle 
!   if(vp > vprof(npoints))exit
!   id=map1(vprof,all_rprofs(1:npoints,j),npoints,vp,profint,1)
!   id=map1(vprof,alphas(1:npoints,j),npoints,vp,al,1)
!   model_spectrum(i)=model_spectrum(i)+profint*al
!   enddo
!enddo

!call mini(x,r_obs,weight,npoints+1,nobs,isp,fun)

!regpar(1) = 0.1d0
!regpar(2) = 0.1d0
if (regpar(1) == 0.d0) then
	call mini(x,r_obs,weight,nns*nSLD*npoints+1,nobs,isp,fun)
else 
	allocate(absc(nns*nSLD*npoints+1))
	absc(1) = 0.d0
	!absc(2:npoints+1) = vprof
	!absc(npoints+2:2*npoints+1) = vprof 
	do i = 1, nns
		do j = 1, nSLD
			do k = 1, npoints
				absc(1+(i-1)*nSLD*npoints + (j-1)*npoints + k) = vprof(k) + Vrad
			enddo
		enddo
	 enddo
	call mini_reg(x,absc,r_obs,weight,regpar,nns*nSLD*npoints+1,nobs,isp,fun, npoints)
endif

write(*,*)x(1) ! Vrad
Vrad = x(1)

open(100,file = 'model.dat')
do i=1,nobs
write(100,*)wl_obs(i),1.d0-model_spectrum(i),1.d0-r_obs(i)
enddo
close(100)

open(101,file = 'prof.dat')
do k=1,npoints
	write(101, '(1x,<nSLD*nns+1>f10.5)') vprof(k)+Vrad, ((x(1+(j-1)*nSLD*npoints + (i-1)*npoints + k) , i = 1, nSLD), j = 1, nns)
enddo
close(101)

end

subroutine get_Vrad
	use common_data
	real(8) CCF1, CCF2, Rv1, Rv2, i
	external CCF_max
	!tolerant=1.d-5; ! Code should correct radial velocity by search maximum CCF in region Rv_old +- dv

	!Rv1=fmin(-200.d0,0.d0,CCF_max,tolerant,CCF1, 1)
	!Rv2=fmin(0.d0,200.d0,CCF_max,tolerant,CCF2, 1)
	!Vrad=Rv1
	!if (CCF1 > CCF2) Vrad=Rv2

	CCF1 = CCF_max(0.d0, 1)
	Vrad = 0.d0
	do i = -200.d0, 200.d0, 1.d0
		CCF2 = CCF_max(i, 1)
		if (CCF1 < CCF2) then
			Vrad=i
			CCF1 = CCF2
		endif
		!print *, CCF2, CCF1, i
	enddo

end subroutine get_Vrad


real(8) function CCF_max(v, n)
	use common_data
	implicit real(8) (a-h,o-z)

	allocate(xint(nobs), yint(nobs), stat=ios)
	if(ios /= 0) stop 'Memory allocation failed (xint/yint in subroutine correction)'
	xint = wl_obs*(1-v/vc)
	id=map1(wl_synth_full,r_synth_full,nsynth_full,xint,yint,nobs)

	sx=sum(r_obs**2); sy=sum(yint**2); sxy=sum(r_obs*yint)

	if(sx <= 0.or.sy <= 0.or.sxy <= 0)then
		CCF_max=100
	else  
		CCF_max=sqrt(sx)*sqrt(sy)/sxy ! really 1/CCF - fmin search minimum of function, so CCF_max=1/CCF
	endif
	deallocate (xint,yint)
end


subroutine gauss_fun_finite(x,f,jacobian,nx,mode,nspec) !get_model_spectrum
use common_data
use dfport
implicit real(8) (a-h,o-z)
character(4) line_ident,mode
real(8) x(nx),xf(nx),f(nobs),jacobian(nobs,nx)
data c/2.997925d5/
	
	SLD_prof=0.d0
	model_spectrum=0.d0
	f=0.d0
	if(mode == 'grad')jacobian=0.d0
	!Vrad=0.d0;
	Vrad = x(2) 
	dVrad=5.d0
	do i = 1,size(SLD_prof)
		SLD_prof(i) = G(vprof(i),x)
	enddo

 do n=1,nlines
  !read(2,rec=n+2)central_wavelength,line_ident,rprof
  central_wavelength = wlines(n)
  rprof = all_rprofs(1:npoints,n)
  !rprof=1.d0-rprof
  loc_alphas = alphas(1:npoints,n)
  call convolve
  do i=1,nobs
   vp=((wl_obs(i)-central_wavelength)/central_wavelength)*c !  Corresponded for given wavelength velocity
   if(vp < vprof(1))cycle 
   if(vp > vprof(npoints))exit
   id=map1(vprof,ySLD,npoints,vp,profint,1)
   model_spectrum(i)=model_spectrum(i)+profint
  enddo
 enddo


 f=model_spectrum
! Jacobian
	
	if (mode == 'grad') then
		
		! x1
		model_spectrum = 0.d0
		xf = x
		xf(1) = xf(1)+ 0.001d0
		do i = 1,size(SLD_prof)
			SLD_prof(i) = G(vprof(i),xf)
		enddo
		 do n=1,nlines
		  !read(2,rec=n+2)central_wavelength,line_ident,rprof
		  central_wavelength = wlines(n)
		  rprof = all_rprofs(1:npoints,n)
		  !rprof=1.d0-rprof
		  loc_alphas = alphas(1:npoints,n)
		  call convolve
		  do i=1,nobs
		   vp=((wl_obs(i)-central_wavelength)/central_wavelength)*c !  Corresponded for given wavelength velocity
		   if(vp < vprof(1))cycle 
		   if(vp > vprof(npoints))exit
		   id=map1(vprof,ySLD,npoints,vp,profint,1)
		   model_spectrum(i)=model_spectrum(i)+profint
		  enddo
		 enddo
		 do i = 1,nobs
		  jacobian(i,1) = (model_spectrum(i) - f(i))/0.001d0
		 enddo

		! x2
		model_spectrum = 0.d0
		xf = x
		xf(2) = xf(2)+ dVrad
		do i = 1,size(SLD_prof)
			SLD_prof(i) = G(vprof(i),xf)
		enddo
		 do n=1,nlines
		  !read(2,rec=n+2)central_wavelength,line_ident,rprof
		  central_wavelength = wlines(n)
		  rprof = all_rprofs(1:npoints,n)
		  !rprof=1.d0-rprof
		  loc_alphas = alphas(1:npoints,n)
		  call convolve
		  do i=1,nobs
		   vp=((wl_obs(i)-central_wavelength)/central_wavelength)*c !  Corresponded for given wavelength velocity
		   if(vp < vprof(1))cycle 
		   if(vp > vprof(npoints))exit
		   id=map1(vprof,ySLD,npoints,vp,profint,1)
		   model_spectrum(i)=model_spectrum(i)+profint
		  enddo
		 enddo
		 do i = 1,nobs
		  jacobian(i,2) = (model_spectrum(i) - f(i))/dVrad
		 enddo

		 
		! x3
		model_spectrum = 0.d0
		xf = x
		xf(3) = xf(3)+ 0.001d0
		do i = 1,size(SLD_prof)
			SLD_prof(i) = G(vprof(i),xf)
		enddo
		 do n=1,nlines
		  !read(2,rec=n+2)central_wavelength,line_ident,rprof
		  central_wavelength = wlines(n)
		  rprof = all_rprofs(1:npoints,n)
		  !rprof=1.d0-rprof
		  loc_alphas = alphas(1:npoints,n)
		  call convolve
		  do i=1,nobs
		   vp=((wl_obs(i)-central_wavelength)/central_wavelength)*c !  Corresponded for given wavelength velocity
		   if(vp < vprof(1))cycle 
		   if(vp > vprof(npoints))exit
		   id=map1(vprof,ySLD,npoints,vp,profint,1)
		   model_spectrum(i)=model_spectrum(i)+profint
		  enddo
		 enddo
		 do i = 1,nobs
		  jacobian(i,3) = (model_spectrum(i) - f(i))/0.001d0
		 enddo

	endif

end

function G(x, gx)
	real(8) G, gx(3), x

	G = gx(1) * exp(-(x-gx(2))**2/gx(3)**2)
end function

subroutine init
use common_data
implicit real(8) (a-h,o-z)
character(100) fn, synth_name, tempstr
character(4) line_ident

vc = 299792.d0

 open(1,file='SLD.config',status='old',iostat=ios)
 if(ios /= 0)stop 'SLD.config missed' 
 read(1,'((a))') fn
 read(1,'((a))') synth_name
 cross_corr = (synth_name(1:4) /= 'skip')
 read(1,*) br
 read(1,*) regpar


 read(1, *) nSLD
 if(nSLD > 1) then
 	allocate(LS_limits(nSLD-1), LS_mean(nSLD))
 	read(1,*) (LS_limits(i), i = 1, nSLD-1)
 else
	read(1,*)
 endif

 if(nSLD > 1) then
	do ilsd = 1, nSLD
		if(ilsd == 1) then
			LS_mean(ilsd) = (0.d0 + LS_limits(ilsd))*0.5d0
		elseif(ilsd == nSLD) then
			LS_mean(ilsd) = (1.d0 + LS_limits(ilsd-1))*0.5d0
		else
			LS_mean(ilsd) = (LS_limits(ilsd-1) + LS_limits(ilsd))*0.5d0
		endif
	enddo
 endif

 read(1,*) nspecial
 if (nspecial > 0) then
  allocate(special_idents(nspecial))
  read(1,'((a))') tempstr
  do i = 1, nspecial
   special_idents(i) = tempstr(1+(i-1)*4:4+(i-1)*4)
  enddo
 endif

 open(2,file='profiles.bnr',status='old',form='binary',access='direct',recl=12,iostat=ios) !Initialize lines profile database
 if(ios /= 0)stop 'cannot open profiles.bnr file'
 read(2,rec=1)nlines,npoints,nmiddle
 print *, nlines, npoints, nmiddle
 allocate(vprof_ext(npoints))
 close (2)

 open(2,file='profiles.bnr',form='binary',access='direct',recl=npoints*8*2+12)
 read(2,rec=2) vprof_ext

if (br == 0.d0) then
 allocate(vprof(npoints))
 vprof = vprof_ext

 allocate(alphas(npoints,nlines))
 allocate(rprof(npoints))
 allocate(wlines(nlines), line_idents(nlines), all_rprofs(npoints,nlines),loc_alphas(npoints)) 
 do i = 1,nlines
	read(2,rec=i+2) wlines(i), line_idents(i), all_rprofs(1:npoints,i), alphas(1:npoints,i)
	all_rprofs(1:npoints,i) = 1.d0 - all_rprofs(1:npoints,i)
 enddo

else
 hdv = abs(vprof_ext(2)-vprof_ext(1))
 nvs = 2*(br+5.d0)/hdv - 1
 allocate(vprof(nvs))
 vprof((nvs-1)/2+1) = 0.d0
 do i = 1,(nvs-1)/2	!	extended velocity grid
  vprof((nvs-1)/2+1 - i) = - i * hdv
  vprof((nvs-1)/2+1 + i) = i * hdv
 enddo

 nmiddle = (nvs-1)/2+1
 nbuff = npoints
 npoints = nvs
 nvs = nbuff

 allocate(alphas(npoints,nlines))
 allocate(rprof(npoints))
 allocate(wlines(nlines), line_idents(nlines), all_rprofs(npoints,nlines),loc_alphas(npoints))

 all_rprofs = 0.d0
 alphas = 1.d0
 do i = 1,nlines
	read(2,rec=i+2) wlines(i), line_idents(i), all_rprofs(nmiddle-nvs/2+1 : nmiddle+nvs/2,i), alphas(nmiddle-nvs/2+1 : nmiddle+nvs/2,i)
	all_rprofs(nmiddle-nvs/2+1 : nmiddle+nvs/2,i) = 1.d0 - all_rprofs(nmiddle-nvs/2+1 : nmiddle+nvs/2,i)
 enddo

 close(2)
 nmiddle = nmiddle - 1
endif

 k= index(fn,'*')
 if(k == 0)then
  call load_obs(fn)
 else             ! Cycle from all files with observed spectra
 endif

 allocate (model_spectrum(nobs),SLD_prof(npoints),ySLD(npoints))

 allocate (loc_grad(npoints,npoints))

if (cross_corr) then
 open(10,file=synth_name)
 nsynth=0
 do
  read(10,*,iostat=ios)
  if(ios /= 0)exit
  nsynth=nsynth+1
 enddo
 rewind (10)
 nsynth_full = nsynth
 if(nsynth > 1000) nsynth = 1000
 allocate (wl_synth(nsynth),r_synth(nsynth))
 do i=1,nsynth
  read(10,*)wl_synth(i),r_synth(i)
 enddo
 r_synth=1.d0-r_synth
 rewind(10)
 allocate(wl_synth_full(nsynth_full), r_synth_full(nsynth_full))
 do i=1,nsynth_full
  read(10,*)wl_synth_full(i),r_synth_full(i)
 enddo
 close (10)
endif

 close(1)

 allocate(lsd_temp(npoints,nSLD), flux_sld(npoints))

end

subroutine load_obs(fn) ! Load given observed spectrum
use common_data
character(80) fn

 open(10,file=fn)
 nobs=0
 do
  read(10,*,iostat=ios)
  if(ios /= 0)exit
  nobs=nobs+1
 enddo
 rewind (10)
 allocate (wl_obs(nobs),r_obs(nobs))
 do i=1,nobs
  read(10,*)wl_obs(i),r_obs(i)
 enddo
 r_obs=1.d0-r_obs
 close (10)

end

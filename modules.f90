module common_data
real(8), allocatable :: xint(:), yint(:)

real(8), allocatable :: wl_obs(:),r_obs(:),wl_synth(:),r_synth(:), wl_synth_full(:), r_synth_full(:) ! Observed spectrum
integer nobs, nsynth, nsynth_full                      ! number of points of observed spectrum

real(8), allocatable :: SLD_prof(:),ySLD(:) ! SLD profile and same convolved
real(8) Vrad
real(8), allocatable :: model_spectrum(:)
integer nlines,npoints,nmiddle ! fixed length SLD database data: number of lines, number of line profile points
real(8), allocatable :: vprof(:),rprof(:),loc_alphas(:), vprof_ext(:), absc(:) ! single line profile
real(8), allocatable :: alphas(:,:), wlines(:), all_rprofs(:,:)
character(4), allocatable :: line_idents(:)
real(8) regpar(2)
! Jacobian data
real(8), allocatable :: loc_grad(:,:) ! gradients of single line profile
!real(8), allocatable :: temp(:) ! temporary array which need in calculations of jacobian

integer nmasks, nspecial, nSLD
character(4), allocatable :: special_idents(:)
real(8), allocatable :: LS_limits(:), LS_mean(:), flux_sld(:), lsd_temp(:,:)

real(8) vc 

logical cross_corr

end module

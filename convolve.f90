subroutine convolve ! Convolution of SLD and line profile
use common_data
implicit real(8) (a-h,o-z)

!print *, rprof
!print *, '!!'
!SLD_prof = 0.d0
!SLD_prof(nmiddle) = 1.d0
 !snorm=sum(rprof)

!print *, rprof

 loc_grad=0.d0 
 do i=1,npoints
  ss=0.d0
   do k=1,npoints
    ik=i-k
	if(ik <= 0.and.abs(ik) < nmiddle)m=nmiddle+ik
	if(abs(ik) > nmiddle)cycle
	if(ik > 0.and.ik <= nmiddle)m=nmiddle+ik
    ss=ss+SLD_prof(k)*rprof(m)*loc_alphas(m)
    !loc_grad(i,k)=rprof(m)/snorm*loc_alphas(m)
	loc_grad(i,k)=rprof(m)*loc_alphas(m)
   enddo
  !ySLD(i)=ss/snorm
  ySLD(i)=ss
 enddo

end
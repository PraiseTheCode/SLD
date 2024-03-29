!
! modified Levenberg-Marquardt minimization routine MINI
!
! REAL*8 XOPT(N)     : initial guess for the vector of unknowns
! REAL*8 ABSC(N)     : abscissas (e.g. log tau) for the vector of unknowns
! REAL*8 OBS(NWL)    : observations
! REAL*8 WEIGHT(NWL) : normalized weights, computed from S/N as follows
!                      WEIGHT(IWL)=SNR(IWL)*SNR(IWL); MAXWEIGHT=MAXVAL(WEIGHT)
!                      WEIGHT(IWL)=WEIGHT(IWL)/MAXWEIGHT
! REAL*8 REGPAR(2)   : Tikhonov regularization parameter
!                      REGPAR(1) gives absolute scale of the regularization
!                      REGPAR(2) gives fractional reduction of regularization
!                      according to total sensitivity function computed
!                      using array of the first derivatives
! INTEGER N          : number of unknowns
! INTEGER NWL        : number of spectral points
!
! on exit MINI returns the vector of uknowns in XOPT
!
! MINI calls a subroutine FUNC which should return theoretical
! profiles and their first derivaties
!
! CALL FUNC(X,PRF,DERIVA,TYPE)
!
! where
! 
! REAL*8 X(N)          : vector of unknowns
! REAL*8 PRF(NWL)      : theoretical lines profiles computed for the same set
!                        of points as observations
! REAL*8 DERIVA(NWL,N) : array of first derivatives; DERIVA(I,J) contains the
!                        value of derivative of spectral point I with repect to
!                        parameter J
! CHARACTER*4 TYPE     : if TYPE='GRAD', FUNC should compute PRF and DERIVA
!                        if TYPE='FUNC', FUNC should only compute PRF
!

subroutine mini_reg(xopt,absc,obs,weight,regpar,n,nwl,nnspec,fun,npoints) ! Tichonov's regularization included

real(8) xopt(n), obs(nwl), weight(nwl),absc(n),regpar(2),admin,admax,regul 
real(8) ftoler, freduc, fvalue, fopt, dev, oldfun, d, cormax, corlim, corscl, di, qlambda, qscale
logical failure, done
integer info, wrksiz, iterf, maxitr, iterc, nspec, nnspec
parameter (ftoler = 0.0005d0, freduc = 0.9995d0, corlim = 3.d0)
real(8), allocatable :: x(:), prf(:), work(:), sigc(:,:), diag(:), rhsold(:), deriva(:,:), derivc(:,:), rhs(:),aderiv(:)
integer, allocatable :: ipiv(:)


allocate (x(n), prf(nwl), ipiv(n), work(45*n), sigc(n,n), diag(n), rhsold(n), deriva(nwl,n), derivc(n,n), rhs(n),aderiv(n),stat=ios)
if(ios /= 0)stop 'Cannot allocate arrays in mini_reg subroutine'
work=0.d0

!  Set initial values to counters and flags. ITERC and NCALLS
!  count the number of iterations and calls of GRFLUX.

iterc   = 0
deriva  = 0.d0
wrksiz  = 45 * n
maxitr  = 50
ncalls  = 0
failure = .false.
qlambda = 10.0
qscale  = 10.0
fopt    = 1.d30
iterf   = iterc

!  The modified Marquardt loop stars here. The big loop computes the
!  function and derivatives using FUN(... ,'GRAD') call.
!  The inner loop adjusts the Marquardt parameter QLAMBDA and compares
!  the function values using the call FUN(... ,'FUNC').
!  So we start with the main loop:

!  Copy the starting guess

do 
 iterc = iterc + 1

 call dcopy(n,xopt,1,x,1)

 oldfun = fopt
  
!  Compute function and straight derivative. Rescale derivative to account
!  for weights and normalization of misfit function.
 
 call fun(x,prf,deriva,n,'grad',nnspec)

 do i = 1, n
  do iwl = 1, nwl
   deriva(iwl,i) = deriva(iwl,i) * sqrt(weight(iwl) / nwl)
  enddo
 enddo  
      
!  Compute misfit function and average deviation
 
 call deviat(obs,prf,weight,nwl,d,dev)

!  After evaluation of gradient compute total relative sensitivity function
!  which will be used for rescaling regularization

  do i=1,n
   aderiv(i)=0.d0
   do iwl=1,nwl    
    aderiv(i)=aderiv(i)+deriva(iwl,i)**2
   enddo
  enddo  	 
  aderiv=sqrt(aderiv)
  admin=minval(aderiv)  
  admax=maxval(aderiv)  
  do i=1,n
   aderiv(i)=(aderiv(i)-admin)/(admax-admin) 
  enddo

!  Compute regularization and total function value
  call regfunq(x,absc,aderiv,regpar,n,regul,npoints)
  fvalue=d+regul 

!  Copy the best guess and function value; ignore the first iteration

 if(fvalue < fopt) then
  call dcopy(n,x,1,xopt,1)
  fopt = fvalue
 endif
 
!  Compute the symmetric curvature matrix DERIVC = DERIVA' # DERIVA 
 
 call dsyrk('u','t',n,nwl,1.d0,deriva,nwl,0.d0,derivc,n)
 do i1 = 2, n
  do i2 = 1, i1-1
   derivc(i1,i2) = derivc(i2,i1)
  enddo
 enddo

!  Compute the Right Hand Side RHS and add the second derivatives to the
!  main diagonal of trace(DERIVC)=trace(DERIVC)+sum_wl d^2C/dx^2 * [O-C]

 do i = 1, n
  di = 0.d0
  do iwl = 1, nwl
   di = di + sqrt(weight(iwl)) * (obs(iwl) - prf(iwl)) * deriva(iwl,i)
  enddo
  rhs(i) = di / sqrt(dble(nwl))
 enddo

!  Derivatives of the regularization functions
 call regder(x,absc,aderiv,regpar,n,derivc,rhs)

!  Make backup copies of the main diagonal of the curvatire matrix and
!  of the right hand side in case the current QLAMBDA fails.

 call dcopy(n,derivc,n+1,diag,1)
 call dcopy(n,rhs,1,rhsold,1)

!  Internal loop: apply Lambda correction

 done = .false.
 ifail = 0
 itry = 0
 oldfun = fvalue
  
 do
  call dscal(n,1.d0+qlambda,derivc,n+1)
  ifail = ifail + 1
  itry = itry + 1 

!  Solve the system of linear equations
  
  call dsysv('u',n,1,derivc,n,ipiv,rhs,n,work,wrksiz,info)
  if(info /= 0) failure = .true.

!  Check if correction was properly computed

  if(failure) then
   write(*,'('' Degenerate matrix DERIVC.'','' Cannot solve SLE.'',T95,''*** (MINI  ) ***'')')
   write(*,*) 'INFO=', info  
   deallocate (x, prf, ipiv, work, sigc, diag, rhsold, deriva, derivc, rhs)
   stop
  endif

!  Find maximum correction
  
  cormax = 0.d0
  do i = 1, n
   if(cormax < abs(rhs(i))) cormax = abs(rhs(i))
  enddo

!  Apply the correction to X(I)

  corscl = 1.d0
  if(cormax > corlim) corscl = corlim / cormax
  do i = 1, n
   x(i) = x(i) + rhs(i) * corscl
  enddo

!  Evaluate the function only

  call fun(x,prf,deriva,n,'func',nnspec)
  call deviat(obs,prf,weight,nwl,d,dev)
  call regfunq(x,absc,aderiv,regpar,n,regul,npoints)
  fvalue = d+regul; dev = sqrt(dev)

!  Try changing QLAMBDA. If the last attempt was not useful, change
!  direction. The first QSCALE in case of change of direction is to
!  undo the previous iteration.

  if(itry == 2 .and. fvalue > oldfun*freduc) then
   if(qscale > 1.) then
    qscale = 0.1
   else
    qscale = 10.
   endif
   qlambda = qlambda * qscale * qscale
   call dcopy(n,xopt,1,x,1)
   fvalue = fopt
   done = .false.
  elseif(itry > 2 .and. fvalue > oldfun*freduc) then
   qlambda = qlambda / qscale
   done = .true.   
  else
   qlambda = qlambda * qscale
   if(qlambda > 1.e10) then
    qlambda = 1.e9
	done = .true.
   elseif(qlambda < 1.e-10) then
    qlambda = 1.e-9
	done = .true.
   else
    done = .false.
   endif
  endif
     
!  Copy the best guess and function value

  if(fvalue < fopt*(1.d0-ftoler)) then
   call deviat(obs,prf,weight,nwl,di,dev)
   dev = sqrt(dev)
   call dcopy(n,x,1,xopt,1)
   fopt = fvalue
   ifail = 0
  endif

!  Restore curvature matrix and the right hand side

  do i1 = 2, n
   do i2 = 1, i1-1
    derivc(i2,i1) = derivc(i1,i2)
   enddo
  enddo
  call dcopy(n,diag,1,derivc,n+1)
  call dcopy(n,rhsold,1,rhs,1)

!  Continue internal loop
  
  if(.not. done .or. itry <= 2 .or. fvalue < oldfun*freduc) then
   oldfun = fvalue
   cycle
  else
   exit ! From internal loop
  endif
 
 enddo

!  Continue external loop

 write(*,"(' ====  ITER:',I3,' MEAN DEV.=',F7.3,'% FUNCT.=',G12.5,' FOPT=',G12.5)") iterc, dev*100, fvalue, fopt
 write(500,"(' ====  ITER:',I3,' MEAN DEV.=',F7.3,'% FUNCT.=',G12.5,' FOPT=',G12.5)") iterc, dev*100, fvalue, fopt
 if(ifail < itry .or. (iterc - iterf) < 2) then
  call dcopy(n,xopt,1,x,1)
  fvalue = fopt
  oldfun = fopt 

!  Check if the limit of iteration is reached. If it is, save and exit
  
  if(iterc >= maxitr) then
   write(*,'(''The limit of iteration number has been reached'',t95,''*** (mini  ) ***'')')
   call deviat(obs,prf,weight,nwl,di,dev)     
   dev = sqrt(dev)
   exit ! From main loop
  endif          
 else
  exit ! From main loop
 endif 

enddo ! End of main loop     

!  Normal exit door is here

call dcopy(n,xopt,1,x,1)
  
deallocate (x, prf, ipiv, work, sigc, diag, rhsold, deriva, derivc, rhs,aderiv)

end



subroutine regfunq(x,absc,aderiv,regpar,n,regul, npoints)
real(8) x(n),absc(n),aderiv(n),regpar(2),regul,wt,da,df

!  Tikhonov regularization, quadratic form
 regul=0.d0
  do i=3,n
   i1=i
   i2=i-1
   wt=1.d0-regpar(2)*0.5d0*(aderiv(i1)+aderiv(i2))
   da=absc(i1)-absc(i2)
   df=x(i1)-x(i2)
   if(mod((i-1),npoints) /= 1) regul=regul+wt*(df/da)**2
  enddo 

! Scale by regularization parameter 
  regul=regpar(1)*regul/(n-1) 

end

subroutine regder(x,absc,aderiv,regpar,n,derivc,rhs)
real(8) x(n),absc(n),aderiv(n),regpar(2),rgn,wt,da,df,derivc(n,n),rhs(n)
!  Derivatives of Tikhonov regularization, quadratic form

  do i=3,n
   i1=i
   i2=i-1
   wt=1.d0-regpar(2)*0.5d0*(aderiv(i1)+aderiv(i2))
   da=absc(i1)-absc(i2)
   df=x(i1)-x(i2)
   rgn=regpar(1)*wt/(n-1)/(da*da)
   rhs(i1)=rhs(i1)-rgn*df
   rhs(i2)=rhs(i2)+rgn*df
   derivc(i1,i1)=derivc(i1,i1)+rgn
   derivc(i2,i2)=derivc(i2,i2)+rgn
   derivc(i1,i2)=derivc(i1,i2)-rgn
   derivc(i2,i1)=derivc(i2,i1)-rgn
  enddo
end
 
      

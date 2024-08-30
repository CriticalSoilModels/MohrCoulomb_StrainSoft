! Functions that are used by incrmental driver


module incrementalDriver_funcs
   implicit none
   private
   public :: splitaLine, ReadStepCommons, PARSER, get_increment, USOLVER, EXITNOW
contains


! ==========================================================================
!     basing on an input command with parameters converts  deltaLoad or deltaLoadCirc
!     to the canonical three lists:  dstress(), dstrain(), ifstress()
!     get\_increment is called in each increment (and not once per step )
   subroutine get_increment(keywords, time, deltaTime,ifstress,ninc,      &
      deltaLoadCirc,phase0,deltaLoad,deltaTemp, &  ! AN 2023 temperat
      dtime, ddstress,  dstran , dTemp, Qb33,   &! AN 2023 temperat
      dfgrd0, dfgrd1,drot )
      implicit none
      character(40):: keywords(10)
      integer, intent(in)  :: ifstress(6),ninc
      real(8), intent(in) :: time(2), deltaTime, deltaLoadCirc(6),phase0(6), deltaLoad(9), deltaTemp
      real(8), intent(out) ::  dtime, ddstress(6), dstran(6), Qb33(3,3), dTemp
      real(8), intent(in out) ::  dfgrd0(3,3), dfgrd1(3,3), drot(3,3)


      real(8), parameter :: Pi = 3.1415926535897932385d0
      real(8),parameter,dimension(3,3):: delta = reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))
      real(8),dimension(3,3):: Fb,Fbb, dFb,aux33,dLb,depsb,dOmegab
      real(8):: wd(6),  & ! angular velocity (in future individual for each component)
         w0(6),  & ! initial phase shift for a component
         t         ! step time
      real(8) :: arandom
      integer(4) :: i
      logical :: ok

      dtime =  deltaTime/ ninc
      dTemp = deltaTemp / ninc    ! AN 2023 temperat   (perturbations and random walk may need some programming)
      dstran= 0
      ddstress=0
      Qb33 = delta
      drot = delta
      dfgrd0=delta
      dfgrd1=delta

!------------------------------------------------------
      if(keywords(2) == '*LinearLoad') then                               !  proportional loading
         do i=1,6
            if (ifstress(i)==1)   ddstress(i) = deltaLoad(i)/ ninc
            if (ifstress(i)==0)    dstran(i) = deltaLoad(i)/ ninc              ! log strain -> corresp. displac. inc. not constant
         enddo
! here dfgrd0 and dfgrd1   can be defined from stran assuming polar decomposition F=V.R with R=1  and V = exp(stran)
! for dfgrd0 use stran
! for dfgrd1 use stran-dstran
      endif
!--------------------------------------------------
      if(keywords(2) == '*DeformationGradient') then                     ! full deformation gradient.
! finite rotations calculated after Hughes+Winget 1980
         Fb = reshape((/deltaLoad(1), deltaLoad(5), deltaLoad(7),  &
                      deltaLoad(4), deltaLoad(2), deltaLoad(9),    &
                      deltaLoad(6), deltaLoad(8), deltaLoad(3)/),  &
                      (/3,3/))
                      
         Fbb = delta + (Fb-delta)*(time(1)/deltaTime)
         dfgrd0  = Fbb
         dFb = (Fb-delta)/ninc
         aux33 =  Fbb + dFb/2.0d0
         dfgrd1   = Fbb  + dFb

!  call matrix('inverse', aux33, 3, ok )
         aux33 = inv33(aux33)
         dLb =  matmul(dFb,aux33)
         depsb = 0.5d0*(dLb + transpose(dLb))
         dstran=(/depsb(1,1), depsb(2,2),depsb(3,3), 2.0d0*depsb(1,2),2.0d0*depsb(1,3),2.0d0*depsb(2,3)/)
         dOmegab =    0.5d0*(dLb - transpose(dLb))
         aux33 =  delta - 0.5d0*dOmegab
!     call matrix('inverse', aux33, 3, ok )
         aux33 = inv33(aux33)
         Qb33 = matmul(aux33, (delta+0.5d0*dOmegab))
         drot=Qb33
      endif
!------------------------------------------------------
      if(keywords(2) == '*CirculatingLoad' )then                         !  harmonic oscillation
         wd(:) = 2*Pi/deltaTime
         w0 = phase0
         t= time(1)  + dtime/2   ! step time in the middle of the increment
         do i=1,6
            if(ifstress(i)==1) ddstress(i)= dtime * deltaLoadCirc(i) * wd(i) * Cos(wd(i) * t + w0(i)) + deltaLoad(i)/ ninc
            if(ifstress(i)==0) dstran(i)  = dtime * deltaLoadCirc(i) * wd(i) * Cos(wd(i) * t + w0(i)) + deltaLoad(i)/ ninc
         enddo
! here dfgrd0 and dfgrd1   can be defined from stran assuming polar decomposition F=V.R with R=1  and V = exp(stran)
! for dfgrd0 use stran
! for dfgrd1 use stran-dstran
      endif

!--------------------------------------------------------
      if(keywords(2) == '*PerturbationsS' )then
         ddstress(1)= deltaLoad(1)*cos( time(1)*2*Pi/deltaTime )
         ddstress(2)= deltaLoad(1)*sin( time(1)*2*Pi/deltaTime  )
! here dfgrd0 and dfgrd1   can be defined from stran assuming polar decomposition F=V.R with R=1  and V = exp(stran)
! for dfgrd0 use stran
! for dfgrd1 use stran-dstran
      endif

!--------------------------------------------------------
      if(keywords(2) == '*PerturbationsE' )then
         dstran(1)= deltaLoad(1)*cos( time(1)*2*Pi/deltaTime )
         dstran(2)= deltaLoad(1)*sin( time(1)*2*Pi/deltaTime  )
! here dfgrd0 and dfgrd1   can be defined from stran assuming polar decomposition F=V.R with R=1  and V = exp(stran)
! for dfgrd0 use stran
! for dfgrd1 use stran-dstran
      endif

      if(keywords(2) == '*RandomWalk' )then
         call random_seed
         do i =1,6
            call random_number(arandom)
            if( ifstress(i)== 1) ddstress(i)= 2*(arandom-0.5d0)*deltaLoad(i)
            if( ifstress(i)== 0) dstran(i)= 2*(arandom-0.5d0)*deltaLoad(i)
         enddo
      endif

      return

   contains !=======================================================

!  contained in get\_increment inverts a 3x3 matrix
      function inv33( a )  !==================contained in get\_increment
         implicit none
         real(8), dimension(3,3), intent(in) :: a
         real(8), dimension(3,3) :: b
         real(8), dimension(3,3) :: inv33
         real(8) :: det
         det = - a(1,3)*a(2,2)*a(3,1) + a(1,2)*a(2,3)*a(3,1) &
            + a(1,3)*a(2,1)*a(3,2) - a(1,1)*a(2,3)*a(3,2) &
            - a(1,2)*a(2,1)*a(3,3) + a(1,1)*a(2,2)*a(3,3)

         b= reshape( [-a(2,3)*a(3,2) + a(2,2)*a(3,3), a(1,3)*a(3,2) - a(1,2)*a(3,3), &
            -a(1,3)*a(2,2) + a(1,2)*a(2,3), a(2,3)*a(3,1) - a(2,1)*a(3,3), &
            -a(1,3)*a(3,1) + a(1,1)*a(3,3), a(1,3)*a(2,1) - a(1,1)*a(2,3), &
            -a(2,2)*a(3,1) + a(2,1)*a(3,2), a(1,2)*a(3,1) - a(1,1)*a(3,2), &
            -a(1,2)*a(2,1) + a(1,1)*a(2,2)],                               &
            [3,3])
         inv33 = transpose(b)/det
      end function inv33

   end subroutine get_increment



!    Imitation of utility routine provided by abaqus for people writing  umats
!     rotates a tensor input as vector : if  LSTR == 1 $\to$  stress  or    LSTR == 0 $\to$  strain
   SUBROUTINE ROTSIG(S,R,SPRIME,LSTR,NDI,NSHR)
      implicit none
      integer, intent(in) ::  LSTR,NDI,NSHR
      integer :: ntens
      real(8), dimension(3,3),intent(in) ::  R
      real(8), dimension(1:NDI+NSHR), intent(in) :: S
      real(8), dimension(1:NDI+NSHR) , intent(out):: SPRIME
      real(8)::  a(6), b(3,3)
      ntens = ndi+nshr
      a(:) = 0
      a(1:ntens) = S(:)

      if(LSTR==1) b = reshape( [a(1),a(4),a(5),a(4),a(2),a(6), a(5),a(6),a(3)], [3,3] )
      if(LSTR==0) b = reshape([a(1),a(4)/2,a(5)/2,a(4)/2,a(2),a(6)/2, a(5)/2, a(6)/2, a(3) ],[3,3] )

      b = matmul( matmul(R,b),transpose(R))
      if(LSTR==1) a = [b(1,1),b(2,2),b(3,3),b(1,2),b(1,3),b(2,3)]
      if(LSTR==0) a = [b(1,1),b(2,2),b(3,3),2*b(1,2),2*b(1,3),2*b(2,3)]
      SPRIME = a(1:ntens)
      return
   END  SUBROUTINE ROTSIG


!    Imitation of utility routine provided by abaqus for people writing  umats
!    returns  two  stress invariants
   subroutine SINV(STRESS,SINV1,SINV2,NDI,NSHR)
      implicit none
      real(8),intent(in) :: STRESS(NDI+NSHR)
      real(8),intent(out) ::  SINV1,SINV2
      integer, intent(in) ::  NDI,NSHR
      real(8) :: devia(NDI+NSHR)
      real(8), parameter :: sq2 = 1.4142135623730950488d0
      if(NDI /= 3) stop 'stopped because ndi/=3 in sinv'
      sinv1 = (stress(1) + stress(2) + stress(3) )/3.0d0
      devia(1:3) = stress(1:3) - sinv1
      devia(3+1:3+nshr) = stress(3+1:3+nshr) * sq2
      sinv2 = sqrt(1.5d0 *  dot_product(devia, devia)  )
   end subroutine SINV

!    Imitation of utility routine provided by abaqus for people writing  umats
!    returns  principal values if  LSTR == 1 ->  for stress  or    LSTR == 2 ->   for strain
   subroutine SPRINC(S,PS,LSTR,NDI,NSHR)
      integer, intent(in) :: LSTR,NDI,NSHR
      real(8),intent(in) :: S(NDI+NSHR)
      real(8),intent(out) :: PS(NDI+NSHR)
      real(8):: A(3,3),AN(3,3)
      real(8) :: r(6)
      if(NDI /= 3) stop 'stopped because ndi/=3 in sprinc'
      r(1:3) = s(1:3)
      if(LSTR == 1 .and. nshr > 0) r(4:3+nshr) = s(4:3+nshr)
      if(LSTR == 2 .and. nshr > 0) r(4:3+nshr) = s(4:3+nshr)/2
      A= reshape([r(1),r(4),r(5),r(4),r(2),r(6),r(5),r(6),r(3)],[3,3])
      call spectral_decomposition_of_symmetric(A, PS, AN, 3)
      return
   end subroutine SPRINC

!    Imitation of utility routine provided by abaqus for people writing  umats
!     returns principal directions LSTR == 1 ->   stress  or    LSTR == 2 ->   strain
   subroutine SPRIND(S,PS,AN,LSTR,NDI,NSHR)
      implicit none
      real(8),intent(in) :: S(NDI+NSHR)
      real(8),intent(out) :: PS(3),AN(3,3)
      integer, intent(in) :: LSTR,NDI,NSHR
      real(8):: A(3,3)
      real(8) :: r(6)
      if(NDI /= 3) stop 'stopped because ndi/=3 in sprind'
      r(1:3) = s(1:3)
      if(LSTR == 1 .and. nshr > 0) r(4:3+nshr) = s(4:3+nshr)
      if(LSTR == 2 .and. nshr > 0) r(4:3+nshr) = s(4:3+nshr)/2
      A= reshape([r(1),r(4),r(5),r(4),r(2),r(6),r(5),r(6),r(3)],[3,3])
      call spectral_decomposition_of_symmetric(A, PS, AN, 3)
      return
   end subroutine SPRIND

!    Imitation of quit utility routine provided by abaqus for people writing  umats
   subroutine XIT
      stop 'stopped because umat called XIT'
   end subroutine XIT

!    used by  utility routine SPRINC  or SPRIND
   SUBROUTINE  spectral_decomposition_of_symmetric(A, Lam, G, n)
      implicit none
      integer, intent(in) :: n                                         ! size of the matrix
      real(8), INTENT(in)  :: A(n,n)                                   ! symmetric input matrix  n x n   (not destroyed in this routine)
      real(8), INTENT(out)  :: Lam(n)                                  ! eigenvalues
      real(8), INTENT(out)  :: G(n,n)                                  ! corresponding eigenvectors in columns of G
      integer ::  iter,i, p,q
      real(8) ::   cosine, sine
      real(8), dimension(:), allocatable :: pcol ,qcol
      real(8), dimension(:,:), allocatable :: x

      allocate(pcol(n) ,qcol(n), x(n,n) )
      x = A
      G=0.0d0
      do i=1,n
         G(i,i) = 1.0d0
      enddo

      do  iter = 1,30
         call  get_jacobian_rot(x, p ,q, cosine, sine, n)                !  find how to apply  optimal similarity  mapping
         call  app_jacobian_similarity(x, p,q, cosine, sine, n)          !   perform mapping

         pcol = G(:,p)                                                   !  collect rotations  to global similarity matrix
         qcol = G(:,q)
         G(:,p) =   pcol*cosine - qcol*sine
         G(:,q) =   pcol* sine + qcol *cosine

! here write a problem-oriented accuracy test max\_off\_diagonal < something
! but 30 iterations are usually ok for 3x3 stress or 6x6 stiffness matrix
      enddo

      do i=1,n
         Lam(i) = x(i,i)                                                  !  eigenvalues
      enddo
      deallocate( pcol ,qcol, x )
      return
   end

!    used by  utility routine SPRINC  or SPRIND
   SUBROUTINE  app_jacobian_similarity(A, p,q, c, s, n)              !  jacobian similarity tranformation of a square symmetric matrix A
      implicit none                                                     !  ( $ A : =  G^T .A . G $ with    Givens   rotation  G\_pq = $\{\{c,s\},\{-s,c\}\}$  )
      INTEGER, INTENT(IN)        :: p,q                                 !   G is an identity n x n matrix overridden with  values  \{\{c,s\},\{-s,c\}\}  )
      real(8), INTENT(IN)        :: c ,s                                !   in cells $\{\{pp, pq\},\{qp,qq\}\}$  algorithm according to Kielbasinski  p.385
      integer , INTENT(IN)       :: n
      real(8), dimension(n,n),intent(inout) :: A
      real(8), dimension(n)  :: prow ,qrow
      real(8) :: App, Apq, Aqq

      if(p == q)  stop 'error: jacobian_similarity  p == q'
      if(p<1 .or. p>n) stop 'error: jacobian_similarity p out of range'
      if(q<1 .or. q>n) stop 'error: jacobian_similarity q out of range'

      prow(1:n) = c*A(1:n,p) - s*A(1:n,q)
      qrow(1:n) = s*A(1:n,p) + c*A(1:n,q)
      App = c*c*A(p,p) -2*c*s*A(p,q) + s*s*A(q,q)
      Aqq =  s*s*A(p,p) +2*c*s*A(p,q) + c*c*A(q,q)
      Apq = c*s*(A(p,p) - A(q,q)) + (c*c - s*s)* A(p,q)
      A(p,1:n) =   prow(1:n)
      A(1:n,p) =   prow(1:n)
      A(q,1:n) =   qrow(1:n)
      A(1:n,q) =   qrow(1:n)
      A(p,p) =   App
      A(q,q) =   Aqq
      A(p,q) =   Apq
      A(q,p) =   Apq

   END SUBROUTINE  app_jacobian_similarity


!    used by  utility routine SPRINC  or SPRIND for iterative diagonalization
   SUBROUTINE  get_jacobian_rot(A, p,q, c, s, n)          !   \com returns jacobian similarity  tranformation param.
      implicit none                                          !  \com  for iterative diagonalization of  a square symm.  A
      integer , INTENT(IN)                :: n               !  \com  algorithm according to Kielbasinski 385-386
      real(8), dimension(n,n),intent(in) :: A
      INTEGER, INTENT(OUT)                :: p,q
      real(8), INTENT(OUT)                :: c ,s
      real(8) :: App, Apq, Aqq, d, t, maxoff
      integer ::   i,j

      p = 0
      q = 0
      maxoff  = tiny(maxoff)
      do i=1,n-1
         do j=i+1,n
            if( abs(A(i,j)) > maxoff ) then
               maxoff = abs(A(i,j))
               p=i
               q=j
            endif
         enddo
      enddo
      if (p > 0) then
         App = A(p,p)
         Apq = A(p,q)
         Aqq = A(q,q)
         d = (Aqq - App)/ (2.0d0*Apq)
         t = 1.0d0/ sign(abs(d) + sqrt(1.0d0 + d*d) , d )
         c = 1.0d0/sqrt(1.0d0 + t*t)
         s = t*c
      else                                                               ! \com  no rotation
         p=1
         q=2
         c=1
         s=0
      endif
   end subroutine get_jacobian_rot

   subroutine ReadStepCommons(from, ninc, maxiter,deltaTime, deltaTemp, every)   ! AN 2023 temperat
      integer, intent(in) :: from
      real(8), intent(out) :: deltaTime, deltaTemp ! increase of time and temperature within the whole step
      integer, intent(out) :: ninc,maxiter,every
      logical ::  okSplit
      character(Len=40)   aShortLine, leftLine, rightLine
      read(from,'(a)') aShortLine
      call splitaLine(aShortLine,':',leftLine,rightLine,okSplit )     !  if ':' is absent,   okSplit=False and the whole aShortline is copied to the leftLine
      read(leftLine,*)  ninc, maxiter, deltaTime
      every = 1
      if(okSplit)  read(rightLine,*)  every
      if (every > ninc) every=ninc
      if (every < 1 ) every= 1

      call splitaLine(aShortLine,'T',leftLine,rightLine,okSplit)  ! AN2023  read the step increase of temerature if any
      if(okSplit)  read(rightLine,*)  deltaTemp


   end subroutine ReadStepCommons

!------------------------------------------------------------------------------------------
!  SplitaLine gets aLinie and returns two portions left of the separator sep and right of the separator if sep is found
!  then ok is set to .true.
   subroutine splitaLine(aLine,  sep, left,  right , ok )
      implicit none
      character(len=40), intent(in) ::  aLine
      character(len=40), intent(out):: left,  right
      character(len=40) :: tmp
      character(len=1), intent(in):: sep

      integer:: iSep
      logical:: ok
      ok=.False.
      isep = index(aLine,sep);
      if(isep==0) then
         ok=.False.
         left = trim(adjustl( aLine))
         right = '  '
      endif
      if(isep > 0) then
         ok=.True.
         tmp  =  aLine(:isep-1)
         right = aLine(isep+1:)
         left = tmp
      endif
   end subroutine splitaLine

!------------------------------------------------------------------------------------------
!  reads a condition (= string cond) and returns true if stress stran and statev satisfy this condition
!  it is used after each increment of a step. If cond == true then the remaining increments of a step are skipped
   function EXITNOW(cond, stress,stran,statev,nstatv) !--AN 2016-------------------->
      implicit none
      integer, parameter:: ntens=6, mSummands=5
      integer, intent(in):: nstatv
      real(8), intent(in) :: stress(ntens), stran(ntens),statev(nstatv)
      character(len=40), intent(in) :: cond
      logical:: EXITNOW
      integer:: i,igt, ilt,iis,imin,iplus,iminus,Nsummands,itimes
      character(len=40) :: inp, rhs, summand(mSummands), aux
      real(8):: factor(mSummands),fac,x,y
      real(8), parameter :: sq3  = 1.7320508075689d0, &
         sq23 = 0.81649658092773d0


      exitnow = .False.
      igt = index(cond,'>'); ilt = index(cond,'<');iis = max(igt,ilt)     !  look for a < > sign
      if (iis == 0) goto 555                                              ! correct condition  must contain < or >
      inp = adjustl(cond(:iis)); rhs =  trim(adjustl( cond (iis+1:)))


      factor(1) = 1;
      if(inp(1:1) == '-') then       ! do not treat the first minus as a separator
         factor(1) = -1;  inp=inp(2:)  ! remove the first character = '-' from inp
      endif
      do i=1,mSummands ! loop over all possible summands
         iplus = index(inp,'+');  if(iplus==0) iplus=200      ! position of an operator in the string  set to 200 if this operator is absent
         iminus = index(inp,'-');  if(iminus==0) iminus=200
         igt= index(inp,'>');  if(igt==0) igt=200          ! actually inp cannot contain > or <
         ilt= index(inp,'<');  if(ilt==0) ilt=200

         imin = min(iplus,iminus,igt,ilt) ! choose the first separator
         if(imin==200)  exit               ! no more summands  encountered
         if(imin==iplus) then               ! separator= '+'  everything left from + save as summand  and positive sign for the next summand
            summand(i) = inp(:imin-1) ; factor(i+1) = 1
            inp = inp(imin+1:)
         endif
         if(imin==iminus) then ! separator= '+'  everything left from + save as summand
            summand(i) = inp(:imin-1);  factor(i+1) = -1
            inp = inp(imin+1:)
         endif
         if(imin==ilt .or. imin == igt) then
            summand(i) = inp(:imin-1);  exit
         endif
      enddo
      Nsummands = i  ! last factor(i)*summand(i) was encountered before exit
      x = 0;
      do i=1,Nsummands   ! for each summand  on the LHS
         aux =   adjustl( summand(i) )
         itimes = index(aux,'*')
         if(itimes /= 0) then !  '*' exists: split the summand into factor and component
            read(aux(:itimes-1),*) fac   ! numeric factor  of the summand
            factor(i) = factor(i)*fac   ! the signed numeric factor  of the summand
            aux=trim(adjustl(aux(itimes+1:)))
         endif
         if(itimes == 0)  aux=trim(aux)  ! no '*' aux == component
         select case(aux)
          case ('s1');x = x + factor(i)*stress(1)
          case ('s2');x = x + factor(i)*stress(2)
          case ('s3');x = x + factor(i)*stress(3)
          case ('s12');x = x + factor(i)*stress(4)   ! AN 2020
          case ('s13');x = x + factor(i)*stress(5)   ! AN 2020
          case ('s23');x = x + factor(i)*stress(6)   ! AN 2020


          case ('v1');x = x + factor(i)*statev(1)
          case ('v2');x = x + factor(i)*statev(2)
          case ('v3');x = x + factor(i)*statev(3)
          case ('v4');x = x + factor(i)*statev(4)
          case ('v5');x = x + factor(i)*statev(5)
          case ('v6');x = x + factor(i)*statev(6)
          case ('v7');x = x + factor(i)*statev(7)
          case ('v8');x = x + factor(i)*statev(8)
          case ('v9');x = x + factor(i)*statev(9)

          case ('p');x=x-factor(i)*(stress(1)+stress(2)+stress(3))/3
          case ('q');x = x - factor(i)*(stress(1) - stress(3) )
          case ('P');x=x-factor(i)*(stress(1)+stress(2)+stress(3))/sq3
          case ('Q');x = x - factor(i)*(stress(1) - stress(3) )

          case ('e1');x = x + factor(i)*stran(1)  ! AN 2017
          case ('e2');x = x + factor(i)*stran(2)
          case ('e3');x = x + factor(i)*stran(3)
          case ('g12');x = x + factor(i)*stran(4)   ! AN 2020
          case ('g13');x = x + factor(i)*stran(5)   ! AN 2020
          case ('g23');x = x + factor(i)*stran(6)   ! AN 2020

          case ('ev');x = x - factor(i)* (stran(1)+stran(2)+stran(3))
          case ('eq'); x = x - 2* (stran(1)- stran(3))/3                  ! AN 2017
          case ('eP');x =x- factor(i)*(stran(1)+stran(2)+stran(3))/sq3
          case ('eQ'); x = x - factor(i)* sq23* (stran(1)- stran(3))
          case DEFAULT; goto 555
         end select
      enddo

      read(rhs,*) y
      igt = index(cond,'>'); ilt = index(cond,'<')
      if(igt /= 0)   exitnow  = ( x > y )
      if(ilt /= 0)   exitnow  = ( x < y )
      return
555   write(*,*) 'inp syntax error: ',cond,' exit condition ignored'
      EXITNOW = .False.
   end function  EXITNOW     !<--AN 2016--------------------

!------------------------------------------------------------------------------------------
!  used to read test.inp when the option *ObeyRestrictions is used
   subroutine  PARSER(inputline, Mt,Me,mb)
      implicit none
      character(260), intent(in) ::  inputline(6)
      real(8), dimension(6,6), intent(out) :: Mt , Me
      real(8), dimension(6),intent(out) :: mb

      character(len=260) ::  inp, aux,aux3
      character(40) ::  summand(13)
      integer :: iis,i,iplus,iminus,iequal,imin,iex,itimes,Irestr,ihash, Nsummands
      real(8) :: factor(13),fac


      Mt = 0; Me= 0; mb= 0

      Do Irestr = 1,6 ! Irestr loop over restriction lines

         inp = trim(adjustl(inputline(Irestr)))
         ihash = index(inp,'#')
         if(ihash /= 0) inp = inp(:ihash-1)
         iis = index(inp,'=')
         if(iis==0) stop 'parser error: no = in restriction'

         factor(1) =1;
         if(inp(1:1) == '-') then ! do not treat the first minus as a separator
            factor(1) = -1
            inp=inp(2:)            ! remove the first character = '-' from inp
         endif
         do i=1,13 ! loop over possible summands
            iplus = index(inp,'+');  if(iplus==0) iplus=200
            iminus = index(inp,'-');  if(iminus==0) iminus=200
            iequal = index(inp,'=');  if(iequal==0) iequal=200
            imin = min(iplus,iminus,iequal) ! choose the first separator
            if(imin==200)  stop 'parser err: no +,-,= in restric'
            if(imin==iplus) then ! separator= '+'  everything left from + save as summand
               summand(i) = inp(:imin-1) ; factor(i+1) = 1
               inp = inp(imin+1:)
            endif
            if(imin==iminus) then ! separator= '+'  everything left from + save as summand
               summand(i) = inp(:imin-1);  factor(i+1) = -1
               inp = inp(imin+1:)
            endif
            if(imin==iequal) then ! separator= '='  everything left from + save as summand
               summand(i) = inp(:imin-1);
               inp = inp(imin+1:)  !  rhs possibly with sign
               iminus = index(inp,'-');  if(iminus==0) iminus=200
               iex = index(inp,'!');   if(iex==0) iex=len(inp)+1  ! right limit = comment or EOL
               if(iminus == 200) then     ! '=' is not followed by  '-'
                  factor(i+1) = 1
                  summand(i+1) = inp(:iex-1)
               else                       ! double separator: '=' followed by '-'
                  factor(i+1) = -1
                  summand(i+1) = inp(iminus+1:iex-1)
               endif
               exit   ! reading a single summand after '=' ends reading of the line
            endif
         enddo ! i-loop
         Nsummands=i+1  ! summand()=LHS, summand(Nsummands)=RHS, signs in factor()

         Do i=1,Nsummands-1  ! for summands on the LHS
            aux =   adjustl( summand(i) )
            itimes = index(aux,'*')
            if(itimes /= 0) then              ! if exists '*' then split the summand into factor and component
               read(aux(:itimes-1),*) fac    ! numeric factor  of the summand  ------------------------  TODO it need not be a number it can be a stress s1,s2,s3,s4,s5,s6
               factor(i) = factor(i)*fac     ! the signed numeric factor  of the summand
               aux=adjustl(aux(itimes+1:))
            endif
            aux3 = aux(1:3)
            select case(aux3)
             case ('sd1') ;   Mt(Irestr,1) = factor(i)
             case ('sd2') ;   Mt(Irestr,2) = factor(i)
             case ('sd3') ;   Mt(Irestr,3) = factor(i)
             case ('sd4') ;   Mt(Irestr,4) = factor(i)
             case ('sd5') ;   Mt(Irestr,5) = factor(i)
             case ('sd6') ;   Mt(Irestr,6) = factor(i)
             case ('ed1') ;   Me(Irestr,1) = factor(i)
             case ('ed2') ;   Me(Irestr,2) = factor(i)
             case ('ed3') ;   Me(Irestr,3) = factor(i)
             case ('ed4') ;   Me(Irestr,4) = factor(i)
             case ('ed5') ;   Me(Irestr,5) = factor(i)
             case ('ed6') ;   Me(Irestr,6) = factor(i)
            end select
         enddo
         read(summand(Nsummands) ,*) mb(Irestr)     ! \com RHS numeric without sign
         mb(Irestr) =  mb(Irestr)*factor(Nsummands)  ! \com RHS numeric with sign
      enddo ! Irestr
   end subroutine PARSER

!  solver for unsymmetric matrix and  unknowns on both sides of equation
   subroutine USOLVER(KK,u,rhs,is,ntens) ! 23.7.2008  new usolver with improvement after numerical recipes
!  \com    KK - stiffness  is not spoiled  within the subroutine
!  \com    u - strain   rhs - stress
!  \com    is(i)= 1 means rhs(i) is prescribed,
!  \com    is(i)= 0  means u(i) is prescribed

      implicit none
      integer, intent(in):: ntens
      integer, dimension(1:ntens), intent(in):: is
      real(8), dimension(1:ntens,1:ntens), intent(in):: KK
      real(8), dimension(1:ntens), intent(inout)::  u,rhs
      real(8), dimension(1:ntens):: rhs1
      real(8), allocatable :: rhsPrim(:), KKprim(:,:), uprim(:)
      integer ::  i,j,ii,nis
      integer,allocatable :: is1(:)

      nis = sum(is)                                                     ! \com number of prescribed stress components

      if (all( is(1:ntens)== 0) ) then
         rhs =  matmul(KK,u)
         return
      endif

      if (all(is(1:ntens) == 1)) then                                   ! \com a special case with full stress control
         u =xLittleUnsymmetricSolver(KK,rhs)
         return
      endif

      rhs1 = rhs  ! \com modify the rhs  to rhs1
      do i=1,ntens
         if (is(i) == 0) rhs1 = rhs1 - u(i)*KK(:,i)                        ! \com  modify rhs wherever strain control
      enddo

      allocate(KKprim(nis,nis), rhsprim(nis), uprim(nis), is1(nis))     ! \com re-dimension  stiffness and rhs

      ii=0
      do i=1,ntens
         if(is(i)==1) then
            ii = ii+1
            is1(ii) = i                                                    ! list with positions  of is(i) == 1
         endif
      enddo


      do i=1,nis
         rhsPrim(i) = rhs1( is1(i) )
         do j=1,nis
            KKprim(i,j) =  KK(is1(i),is1(j))
         enddo
      enddo

      if (nis ==1) uprim = rhsprim / KKprim(1,1)
      if (nis > 1) uprim =xLittleUnsymmetricSolver(KKprim,rhsprim)
      do i=1,nis
         u(is1(i)) = uprim(i)
      enddo
      do i=1,ntens
         if ( is(i) == 0 ) rhs(i) = dot_product( KK(i,:), u)             ! \com   calculate rhs where u prescribed
      enddo
      deallocate(KKprim,rhsprim,uprim,is1)


   CONTAINS  !===================================================

!  contained in USOLVER LU-decomposition from NR
      SUBROUTINE ludcmp(a,indx,d)
         IMPLICIT NONE
         REAL(8), DIMENSION(:,:), INTENT(INOUT) :: a
         INTEGER, DIMENSION(:), INTENT(OUT) :: indx
         REAL(8), INTENT(OUT) :: d
         REAL(8), DIMENSION(size(a,1)) :: vv ,aux
         integer, dimension(1) :: imaxlocs
         REAL(8), PARAMETER :: TINY=1.0d-20
         INTEGER  :: j,n,imax
         n = size(a,1)
         d=1.0
         vv=maxval(abs(a),dim=2)
         if (any(vv == 0.0)) stop 'singular matrix in ludcmp'
         vv=1.0d0/vv
         do j=1,n
            imaxlocs=maxloc(  vv(j:n)*abs( a(j:n,j) ) )
            imax=(j-1)+imaxlocs(1)
            if (j /= imax) then
               aux = a(j,:)      ! call swap(a(imax,:),a(j,:))
               a(j,:) = a(imax,:)
               a(imax,:) = aux
               d=-d
               vv(imax)=vv(j)
            end if
            indx(j)=imax
            if (a(j,j) == 0.0) a(j,j)=TINY
            a(j+1:n,j)=a(j+1:n,j)/a(j,j)
            a(j+1:n,j+1:n)=a(j+1:n,j+1:n)- spread(a(j+1:n,j),2,n-j )* spread(a(j,j+1:n),1, n-j)   ! outerprod
         end do
      END SUBROUTINE ludcmp

!  contained in USOLVER LU-back substitution from NR
      SUBROUTINE lubksb(a,indx,b)
         IMPLICIT NONE
         REAL(8), DIMENSION(:,:), INTENT(IN) :: a
         INTEGER, DIMENSION(:), INTENT(IN) :: indx
         REAL(8), DIMENSION(:), INTENT(INOUT) :: b
         INTEGER :: i,n,ii,ll
         REAL(8) :: summ
         n=size(a,1)
         ii=0
         do i=1,n
            ll=indx(i)
            summ=b(ll)
            b(ll)=b(i)
            if (ii /= 0) then
               summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
            else if (summ /= 0.0) then
               ii=i
            end if
            b(i)=summ
         end do
         do i=n,1,-1
            b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
         end do
      END SUBROUTINE lubksb

!  contained in USOLVER  improvement of the accuracy
      SUBROUTINE mprove(a,alud,indx,b,x)
         IMPLICIT NONE
         REAL(8), DIMENSION(:,:), INTENT(IN) :: a,alud
         INTEGER, DIMENSION(:), INTENT(IN) :: indx
         REAL(8), DIMENSION(:), INTENT(IN) :: b
         REAL(8), DIMENSION(:), INTENT(INOUT) :: x
         REAL(8), DIMENSION(size(a,1)) :: r
         r=matmul(a,x)-b
         call lubksb(alud,indx,r)
         x=x-r
      END SUBROUTINE mprove

!  solver contained in USOLVER  for problems with unknowns on the left-hand side
      function xLittleUnsymmetricSolver(a,b)
         IMPLICIT NONE                !==== solves $a . x = b$ \& doesn't spoil a or b
         REAL(8), DIMENSION(:), intent(inout) :: b
         REAL(8), DIMENSION(:,:), intent(in) ::  a
         REAL(8), DIMENSION(size(b,1)) :: x
         REAL(8), DIMENSION(size(b,1),size(b,1)) ::  aa
         INTEGER, DIMENSION(1:size(b,1)) :: indx
         real(8), DIMENSION(1:size(b,1)):: xLittleUnsymmetricSolver
         REAL(8) :: d
         x(:)=b(:)
         aa(:,:)=a(:,:)
         call ludcmp(aa,indx,d)
         call lubksb(aa,indx,x)
         call mprove(a,aa,indx,b,x)
         xLittleUnsymmetricSolver= x(:)
      end function xLittleUnsymmetricSolver


   end subroutine USOLVER



   subroutine stopp(i, whyStopText)                      ! AN 2016
      USE ISO_FORTRAN_ENV  ! , ONLY : ERROR\_UNIT           ! AN 2016
      implicit none                                         ! AN 2016
      integer, intent(in) :: i                             ! AN 2016
      character(*)      :: whyStopText                    ! AN 2016
      stop  'whyStopText'                                   ! AN 2016
      WRITE(ERROR_UNIT,*)   whyStopText                  ! AN 2016
      CALL EXIT(5)                                       ! AN 2016
   end subroutine stopp                             ! AN 2016


end module incrementalDriver_funcs

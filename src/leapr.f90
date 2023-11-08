module leapm
   ! Provides leapr for NJOY2016.
   use locale
   implicit none
   private
   public leapr

   ! global variables

   ! user's input
   integer::nout
   integer::iprint
   integer::nphon
   integer::mat
   real(kr)::smin
   real(kr)::za
   real(kr)::awr
   real(kr)::spr
   integer::npr
   integer::iel
   integer::nss
   integer::ncold,nsk
   real(kr)::b7
   real(kr)::aws
   real(kr)::sps
   integer::mss
   integer::nalpha,nbeta,lat
   real(kr),dimension(:),allocatable::alpha,beta
   real(kr)::twt,c,tbeta
   real(kr),dimension(:),allocatable::tempr
   integer::np1
   real(kr)::delta1
   real(kr),dimension(:),allocatable::p1
   integer::nd
   real(kr),dimension(:),allocatable::bdel,adel
   integer::nka
   real(kr)::dka
   real(kr),dimension(:),allocatable::ska
   real(kr)::cfrac

   ! other global variables for module
   integer::naint,nbint
   real(kr)::f0,tbar
   real(kr)::arat
   real(kr)::tev,deltab
   real(kr),dimension(:,:,:),allocatable::ssm,ssp
   real(kr),dimension(:),allocatable::dwpix,dwp1
   real(kr),dimension(:),allocatable::tempf,tempf1

   ! min phonon expansion for time warning message
   integer,parameter,public::maxnphon=250

contains

   subroutine leapr
   !--------------------------------------------------------------------
   !
   ! Calculate S(alpha,beta)
   !
   ! Calculates the thermal scattering law, S(alpha,beta), in the
   ! incoherent and gaussian approximations.  The scattering law
   ! for solid-type frequency distributions is calculated using
   ! the phonon expansion method without recourse to the usual
   ! edgewood and sct approximations.  If desired, an analytic
   ! representation of diffusion or free-gas scattering can be
   ! convolved with the solid-type scattering law.  In addition,
   ! up to 50 discrete oscillators can be convolved with the
   ! continuous scattering law.  The results of the calculation
   ! are written out in ENDF-6 File 7 format, ready to be
   ! processed by the thermr module of NJOY.
   !
   ! It is possible to generate S(alpha,beta) for composite
   ! moderators like BeO, where Be in BeO is combined with O in
   ! BeO and normalized to be used with the Be cross section.
   !
   ! Incoherent elastic or coherent elastic scattering functions
   ! can also be included using the ENDF-6 format.  The incoherent
   ! result depends on the Debye-Waller factor computed during the
   ! S(alpha,beta) calculation.  The coherent result is computed
   ! using the methods developed for the thermr module of NJOY
   ! (which were based on the HEXSCAT code).  This scattering
   ! law depends on the Debye-Waller factor from the S(alpha,beta)
   ! calculation, on lattice parameters that are built in to data
   ! statements in the code, and on the coherent scattering cross
   ! section (which is also built in).
   !
   ! A special option exists for liquid hydrogen and deuterium.
   ! A solid-type spectrum and a diffusive spectrum can be given
   ! in the normal way.  The resulting S(alpha,beta) is then
   ! convolved with rotational modes calculated using the method
   ! of Young and Koppel.  Because of the inclusion of spin
   ! correlations, the resulting S(alpha,beta) is not symmetric in
   ! beta, and the lasym option is used in MF7.
   !
   ! This module is loosly based on the British code 'LEAP+ADDELT',
   ! originally written by R.C.F.McLatchie at Harwell (1962,
   ! unpublished), then implemented by A.T.D.Butland at Winfrith
   ! (AEEW 1200, 1973), and finally modified to work better for
   ! cold moderators as part of the thesis of D.J.Picton, now
   ! at the University of Birmingham.  The first ENDF and NJOY
   ! compatible version was prepared by R.E.MacFarlane at
   ! Los Alamos in 1987.  The main changes to the original code
   ! were: 1) the change to NJOY style, 2) the addition of ENDF-6
   ! output, 3) the addition of incoherent elastic output, 4) the
   ! addition of a coherent elastic calculation, 5) a major
   ! speed up of the diffusion calculation by using interpolation
   ! instead of direct recalculation of S-solid(alpha,beta),
   ! and 6) the liquid hydrogen and deuterium treatments.
   ! A second version was prepared by R.E.MacFarlane in 1989 by
   ! removing the Edgewood and SCT approximations in favor of
   ! direct use of the phonon expansion for all phonon orders.
   ! In addition, free gas scattering was added, the code was
   ! simplified and scratch tapes were eliminated.  Thus, the
   ! code takes advantage of the capabilities of large, fast
   ! computers that weren't available to the designers of the
   ! original lEAP code.  This 1992 version changed to using
   ! the asymmetric S(alpha,beta) for better numerics on
   ! short-word machines, added the mixed moderator capability,
   ! rebuilt the discrete-oscillator calculation for better
   ! accuracy, and made many other smaller improvements.
   !
   !----- user input (free format) ---------------------------------
   !
   ! card 1 - units
   !    nout     endf output unit for thermal file
   !
   ! card 2 - title
   !
   ! card 3 - run control
   !    ntempr  number of temperatures (def=1)
   !    iprint  print control (0=min, 1=more, 2=most, def=1)
   !    nphon   phonon-expansion order (def=100)
   !
   ! card 4 - endf output control
   !    mat     endf mat number
   !    za      1000*z+a for principal scatterer
   !    isabt   sab type (0=symmetric, 1=asymmetric, def=0)
   !    ilog    log flag (0=s, 1=log10(s), def=0)
   !    smin    minimum S(alpha, beta) stored in file (def=1e-75)
   !
   ! card 5 - principal scatterer control
   !    awr     weight ratio to neutron for principal scatterer
   !    spr     free atom cross section for principal scatterer
   !    npr     number of principal scattering atoms in compound
   !    iel     coherent elastic option
   !                   0  none (default)
   !                   1  graphite
   !                   2  beryllium
   !                   3  beryllium oxide
   !                   4  aluminum
   !                   5  lead
   !                   6  iron
   !    ncold   cold hydrogen option
   !                   0   none (default)
   !                   1   ortho hydrogen
   !                   2   para hydrogen
   !                   3   otho deuterium
   !                   4   para deuterium
   !    nsk            0   none (default)
   !                   1   vineyard
   !                   2   skold
   !
   ! card 6 - secondary scatterer control
   !    nss     number of secondary scatterers (0 or 1)
   !    b7      secondary scatterer type
   !             (0=sct only, 1=free, 2=diffusion)
   !    aws     weight ratio to neutron for secondary scatterer
   !    sps     free atoms cross section for secondary scatterer
   !    mss     number of atoms of this type in the compound
   !
   ! card 7 - alpha, beta control
   !    nalpha   number of alpha values
   !    nbeta    number of beta values
   !    lat      if lat.eq.1, alpha and beta values are scaled
   !               by .0253/tev, where tev is temp in ev.  (def=0)
   !
   ! card 8 - alpha values (increasing order)
   ! card 9 - beta values (increasing order)
   !
   ! scatterer loop, do temperature loop for principal scatterer.
   !         repeat for secondary scatterer (if any) if b7=0.
   !
   ! temperature loop, repeat cards 10 to 18 for each temperature
   !
   !    card 10 - temperature (k)
   !       a negative value means skip cards 11 to 18,
   !          thereby using previous parameters for this temp.
   !
   !    card 11 -- continuous distribution control
   !       delta    interval in ev
   !       ni       number of points
   !
   !    card 12 -- rho(energy) (order of increasing ev)
   !
   !    card 13 - continuous distribution parameters
   !       twt       translational weight
   !       c         diffusion constant (zero for free gas)
   !       tbeta     normalization for continuous part
   !
   !    card 14 - discrete oscillator control
   !       nd     number of discrete oscillators
   !
   !    card 15 - oscillator energies (ev)
   !    card 16 - oscillator weights (sum to 1.-tbeta-twt)
   !
   !    card 17 - pair correlation control (nsk.ne.0 .or. ncold.ne.0)
   !       nka     number of kappa values
   !       dka     kappa increment (inv. angstroms)
   !
   !    card 18  skappa values in increasing order (inv. ang.)
   !
   !    card 19 - coherent scattering fraction (nsk.ne.0)
   !       cfrac   coherent fraction
   !
   ! card 20 - file 1 comments, repeat until blank line is read.
   !
   !--------------------------------------------------------------------
   use mainio  ! provides nysi,nsyso,nsyse
   use physics ! provides bk (boltzmann constant)
   use util    ! provides timer,openz,error,mess
   ! internals
   integer::itemp,idone,i,ntempr,isabt,ilog,ni,nedge
   integer::isym,mscr,maxb,isecs
   real(kr)::time
   character(4)::title(20)
   character(60)::strng
   real(kr)::temp,emax
   character::text*80
   real(kr),dimension(:),allocatable::bragg
   real(kr),dimension(:),allocatable::scr
   real(kr),parameter::zero=0

   !--initialize
   call timer(time)
   write(nsyso,'(/'' leapr...'',&
     &''compute thermal scattering law'',30x,f8.1,''s'')') time
   write(nsyse,'(/'' leapr...'',60x,f8.1,''s'')') time

   !--read user's global input
   read(nsysi,*) nout
   read(nsysi,*) text
   read(text,'(20a4)') (title(i),i=1,20)
   write(nsyso,'(/5x,20a4)') (title(i),i=1,20)
   ntempr=1
   iprint=1
   nphon=100
   read(nsysi,*) ntempr,iprint,nphon
   isabt=0
   ilog=0
   smin=1.0e-75_kr
   read(nsysi,*) mat,za,isabt,ilog,smin
   write(nsyso,'(/&
     &  '' no. of temperatures .................. '',i10/&
     &  '' print flag ........................... '',i10/&
     &  '' phonon-expansion order ............... '',i10/&
     &  '' endf mat number ...................... '',i10/&
     &  '' za ................................... '',i10/&
     &  '' isabt ................................ '',i10/&
     &  '' ilog ................................. '',i10/&
     &  '' smin.................. ............... '',es10.3)')&
     &  ntempr,iprint,nphon,mat,nint(za),isabt,ilog,smin
   if (isabt.ne.0) write(nsyso,'(/&
     &''*** Warning.  isabt=1 pendf tapes CANNOT be processed '',&
     &''by the NJOY THERMR module ***'')')
   if (isabt.ne.0) call mess('leapr','isabt=1 pendf tapes CANNOT be',&
     'processed by thermr.')
   iel=0
   ncold=0
   nsk=0
   read(nsysi,*) awr,spr,npr,iel,ncold,nsk
   write(nsyso,'(/&
     &  '' awr for principal scatterer .......... '',f10.3/&
     &  '' free xsec for principal scatterer .... '',f10.3/&
     &  '' number of principal atoms ............ '',i10/&
     &  '' elastic option ....................... '',i10/&
     &  '' cold moderator option ................ '',i10/&
     &  '' s(kappa) option ...................... '',i10)')&
     &  awr,spr,npr,iel,ncold,nsk
   b7=0
   aws=0
   sps=0
   mss=0
   read(nsysi,*) nss,b7,aws,sps,mss
   if (nss.gt.0) write(nsyso,'(/&
     &  '' secondary scatterer type ............. '',f10.3/&
     &  '' awr for secondary scatterer .......... '',f10.3/&
     &  '' free xsec for secondary scatterer .... '',f10.3/&
     &  '' number of secondary atoms ............ '',i10)')&
     &  b7,aws,sps,mss
   lat=0
   read(nsysi,*) nalpha,nbeta,lat
   naint=0
   nbint=0
   if (naint.eq.0) naint=1
   if (nbint.eq.0) nbint=1
   allocate(alpha(nalpha))
   allocate(beta(nbeta))
   read(nsysi,*) (alpha(i),i=1,nalpha)
   read(nsysi,*) (beta(i),i=1,nbeta)

   ! warn for excessive computation time
   if (nphon.gt.maxnphon) then
      write(strng,'('' phonon expansion order is larger than '',i3)') maxnphon
      call mess('leapr',strng,'calculation time may be excessive')
   endif

   !--open the output unit
   call openz(nout,1)

   !--allocate storage for ssm (and ssp if needed)
   allocate(ssm(nbeta,nalpha,ntempr))
   if (ncold.ne.0) allocate(ssp(nbeta,nalpha,ntempr))

   !--allocate storage for tempr, dwpix, dwp1, tempf and tempf1
   allocate(tempr(ntempr))
   allocate(dwpix(ntempr))
   allocate(dwp1(ntempr))
   allocate(tempf(ntempr))
   allocate(tempf1(ntempr))

   !--loop over desired scatterers and temperatures
   isecs=0
   arat=1
   idone=0
   do while (idone.eq.0)
      if (isecs.eq.0) write(nsyso,'(/'' principal scatterer...'')')
      if (isecs.gt.0) then
         arat=aws/awr
         write(nsyso,'(/'' secondary scatterer...''/&
           &'' input alpha values divided by'',f7.3)') arat
      endif
      do itemp=1,ntempr
         read(nsysi,*) temp
         tempr(itemp)=abs(temp)
         write(nsyso,'(/'' doing temp ='',f10.2)') temp
         tev=bk*abs(temp)
         if (itemp.eq.1.or.temp.ge.zero) then

            !--read in parameters for a continuous distribution
            read(nsysi,*) delta1,ni
            if (allocated(p1)) deallocate(p1)
            allocate(p1(ni))
            read(nsysi,*) (p1(i),i=1,ni)
            np1=ni
            read(nsysi,*) twt,c,tbeta

            !--read in oscillator data
            read(nsysi,*) nd
            if (nd.gt.0) then
               if (allocated(bdel)) deallocate(bdel)
               if (allocated(adel)) deallocate(adel)
               allocate(bdel(nd))
               allocate(adel(nd))
               read(nsysi,*) (bdel(i),i=1,nd)
               read(nsysi,*) (adel(i),i=1,nd)
            endif

            !--read in pair correlation function, s_kappa
            if ((nsk.gt.0).or.(ncold.gt.0)) then
               read(nsysi,*) nka,dka
               if (allocated(ska)) deallocate(ska)
               allocate(ska(nka))
               read(nsysi,*) (ska(i),i=1,nka)
               if (nsk.eq.1) write(nsyso,'(/'' s(kappa) for vineyard method'')')
               if (nsk.eq.2) write(nsyso,'(/'' s(kappa) for skold method'')')
               do i=1,nka
                  write(nsyso,'(1p,2e12.4)') dka*i,ska(i)
               enddo
            endif

            !--read in coherent fraction for skold method
            if (nsk.gt.0) read(nsysi,*) cfrac

         endif

         !--continuous part of distribution
         call contin(temp,itemp,np1,nphon)

         !--translational part, if any
         if (twt.gt.zero) call trans(itemp)

         !--discrete oscillators, if any
         if (nd.gt.0) call discre(itemp)

         !--check for special hydrogen and deuterium options
         if (ncold.gt.0) call coldh(itemp,temp)

         !--check for skold option for correlations
         if ((nsk.eq.2) .and. (ncold.eq.0))&
            call skold(itemp,temp,ssm,nalpha,nbeta,ntempr)

      !--continue temperature loop
      enddo
      if (allocated(adel)) deallocate(adel)
      if (allocated(bdel)) deallocate(bdel)

      !--save ssm for principal scatterer on scratch file
      if (nss.eq.0.or.b7.gt.zero.or.isecs.gt.0) then
         idone=1
      else
         isecs=isecs+1
         call copys(ssm,nbeta,nalpha,ntempr)
         do itemp=1,ntempr
            tempf1(itemp)=tempf(itemp)
            dwp1(itemp)=dwpix(itemp)
         enddo
      endif
   enddo

   !--construct bragg edges if coherent was requested.
   nedge=0
   if (iel.gt.0) then
      maxb=60000
      allocate(bragg(maxb))
      emax=5
      call coher(iel,npr,bragg,nedge,maxb,emax)
   else
      maxb=1
      allocate(bragg(maxb))
   endif

   !--write output in endf format
   isym=0
   if (ncold.ne.0) isym=1
   if (isabt.eq.1) isym=isym+2
   
   ! Based on endout, to write the actual TSL data, the max number of entries
   ! needed in scr is either 8+2*nalpha, or 8+2*nedge. However, we have no way
   ! of knowing how many comment lines were added to the leaper input. The
   ! previous hard coded limit of 4000 is also used as a possible max as this
   ! has apparently been sufficient to hold all comments in the past.
   mscr = max(8 + 2*nalpha, 8 + 2*nedge, 4000)
   allocate(scr(mscr))
   call endout(ntempr,bragg,nedge,maxb,scr,mscr,isym,ilog)

   deallocate(tempr)
   deallocate(p1)
   deallocate(ssm)
   deallocate(scr)
   if (allocated(ssp)) deallocate(ssp)
   deallocate(alpha)
   deallocate(beta)
   deallocate(dwpix)
   deallocate(dwp1)
   deallocate(tempf)
   deallocate(tempf1)

   !--finished
   call closz(nout)
   call timer(time)
   write(nsyso,'(69x,f8.1,''s''/1x,77(''*''))') time
   return
   end subroutine leapr

   subroutine contin(temp,itemp,np,maxn)
   !--------------------------------------------------------------------
   ! Main routine for calculating S(alpha,beta) at temp
   ! for continuous phonon frequency distributions.
   ! Called by leapr.
   ! Uses start, terpt, convol.
   !--------------------------------------------------------------------
   use physics ! provides pi
   use mainio  ! provides nsyso
   use util    ! provides timer
   ! externals
   real(kr)::temp
   integer::itemp,np,maxn
   ! internals
   integer::i,j,k,n,npn,npl,iprt,jprt
   integer,allocatable,dimension(:)::maxt
   character(3)::tag
   real(kr)::al,be,bel,ex,exx,st,add,sc,alp,alw,ssct,ckk,time
   real(kr)::ff0,ff1,ff2,ff1l,ff2l,sum0,sum1
   real(kr),dimension(:),allocatable::p,tlast,tnow,xa
   real(kr),parameter::therm=0.0253e0_kr
   real(kr),parameter::tiny=1.e-30_kr
   real(kr),parameter::explim=-250.e0_kr
   real(kr),parameter::zero=0

   !--write heading
   write(nsyso,'(/'' solid-type contributions to scattering law'')')

   !--allocate temporary arrays
   allocate(p(np1))
   allocate(tlast(nphon*np1))
   allocate(tnow(nphon*np1))
   allocate(xa(nalpha))
   allocate(maxt(nbeta))

   !--calculate various parameters for this temperature
   call start(itemp,p,np,deltab,tev)
   sc=1
   if (lat.eq.1) sc=therm/tev

   !--start the phonon expansion sum with t1
   do i=1,np
      tlast(i)=p(i)
   enddo
   do j=1,nalpha
      al=alpha(j)*sc/arat
      xa(j)=log(al*f0)
      ex=-f0*al+xa(j)
      exx=0
      if (ex.gt.explim) exx=exp(ex)
      do k=1,nbeta
         be=beta(k)*sc
         st=terpt(p,np,deltab,be)
         add=st*exx
         if (add.lt.tiny) add=0
         ssm(k,j,itemp)=add
      enddo
   enddo
   npl=np

   !--do the phonon expansion sum
   do j=1,nbeta
      maxt(j)=nalpha+1
   enddo
   if (iprint.eq.2) then
      write(nsyso,'(/'' normalization check for phonon expansion'')')
   endif
   if (maxn.gt.maxnphon) then
      call timer(time)
      write(nsyse,'(/'' performing phonon expansion sum'',&
            &37x,f8.1,''s'')'),time
   endif
   do n=2,maxn
      npn=np+npl-1
      call convol(p,tlast,tnow,np,npl,npn,deltab,ckk)
      if (iprint.eq.2) write(nsyso,'(5x,i5,f12.5)') n,ckk
      do j=1,nalpha
         al=alpha(j)*sc/arat
         xa(j)=xa(j)+log(al*f0/n)
         ex=-f0*al+xa(j)
         exx=0
         if (ex.gt.explim) exx=exp(ex)
         do k=1,nbeta
            be=beta(k)*sc
            st=terpt(tnow,npn,deltab,be)
            add=st*exx
            if (add.lt.tiny) add=0
            ssm(k,j,itemp)=ssm(k,j,itemp)+add
            if (ssm(k,j,itemp).ne.zero.and.n.ge.maxn) then
            if (add.gt.ssm(k,j,itemp)/1000.and.j.lt.maxt(k)) maxt(k)=j
            endif
         enddo
      enddo
      do i=1,npn
         tlast(i)=tnow(i)
      enddo
      npl=npn
      if (mod(n,maxnphon).eq.0) then
         call timer(time)
         write(nsyse,'(2x,i5,'' of '',i5,&
               &'' loops done for phonon expansion sum'',17x,f8.1,''s'')')&
               &n,maxn,time
      endif
   enddo
   if (maxn.gt.maxnphon) then
      call timer(time)
      write(nsyse,'(/'' done with phonon expansion sum'',&
            &38x,f8.1,''s'')'),time
   endif

   !--print out start of sct range for each beta
   if (iprint.ne.0) then
      write(nsyso,'(/''         beta   alpha sct'')')
      do i=1,nbeta
         if (i.gt.1) then
            if (maxt(i).gt.maxt(i-1)) maxt(i)=maxt(i-1)
         endif
         if (maxt(i).gt.nalpha) then
            write(nsyso,'(1x,f12.4,''      none'')') beta(i)*sc
         else
            write(nsyso,'(1x,2f12.4)')&
              beta(i)*sc,alpha(maxt(i))*sc/arat
         endif
      enddo
   endif

   !---check the moments of s(alpha,beta)
   if (iprint.eq.1) write(nsyso,&
     '(/'' sab checks''/''      alpha      norm   sum rule'')')
   do j=1,nalpha
      iprt=mod(j-1,naint)+1
      if (j.eq.nalpha) iprt=1
      if (iprt.eq.1) then
         al=alpha(j)*sc/arat
         if (iprint.eq.2) then
            write(nsyso,'(/''  alpha='',f10.4)') al
            write(nsyso,'(5x,''beta'',7x,''s(alpha,beta)'',&
              &6x,''ss(alpha,beta)'',5x,''ss(alpha,-beta)'')')
         endif
         bel=0
         ff1l=0
         ff2l=0
         sum0=0
         sum1=0
         do k=1,nbeta
            jprt=mod(k-1,nbint)+1
            if (k.eq.nbeta) jprt=1
            be=beta(k)*sc
            alw=al*tbeta
            alp=alw*tbar
            ex=-(alw-be)**2/(4*alp)
            ssct=0
            if (ex.gt.explim) ssct=exp(ex)/sqrt(4*pi*alp)
            tag='   '
            if (j.ge.maxt(k)) then
               tag='sct'
               ssm(k,j,itemp)=ssct
            endif
            ff2=ssm(k,j,itemp)
            ff1=ssm(k,j,itemp)*exp(-be)
            ff0=ssm(k,j,itemp)*exp(-be/2)
            if (jprt.eq.1.and.iprint.eq.2.and.ff2.gt.zero)&
              write(nsyso,'(f10.4,1p,e18.5,2e20.5,4x,a3)')&
              be,ff0,ff1,ff2,tag
            if (k.gt.1) then
               sum0=sum0+(be-bel)*(ff1l+ff2l+ff1+ff2)/2
               sum1=sum1+(be-bel)*(ff2l*bel+ff2*be-ff1l*bel-ff1*be)/2
               ff1l=ff1
               ff2l=ff2
               bel=be
            else
               bel=be
               ff1l=ff1
               ff2l=ff2
               sum0=0
               sum1=0
            endif
         enddo
         sum0=sum0/(1-exp(-al*f0))
         sum1=sum1/al/tbeta
         if (iprint.eq.2) then
            write(nsyso,'(''   normalization check ='',f8.4)') sum0
            write(nsyso,'(''        sum rule check ='',f8.4)') sum1
         else if (iprint.eq.1) then
            write(nsyso,'(1x,3f10.4)') al,sum0,sum1
         endif
      endif
   enddo

   !--finished with continuous distribution
   end subroutine contin

   subroutine start(itemp,p,np,deltab,tev)
   !--------------------------------------------------------------------
   ! Computes several integral functions of the
   ! phonon frequency distribution.
   ! Called by contin.
   ! Uses fsum.
   !--------------------------------------------------------------------
   use mainio ! provides nsyso
   ! externals
   integer::itemp,np
   real(kr)::deltab,tev
   real(kr)::p(np)
   ! internals
   integer::i,j,npt
   real(kr)::u,v,vv,tau,an,be,z,ee,bigp,rho,rhoe

   !--copy input spectrum into p array
   deltab=delta1/tev
   do i=1,np1
      p(i)=p1(i)
   enddo
   npt=np1
   u=deltab
   v=exp(deltab/2)
   p(1)=p(2)/deltab**2
   vv=v
   do j=2,np1
      p(j)=p(j)/(u*(vv-1/vv))
      vv=v*vv
      u=u+deltab
   enddo

   !--calculate normalizing constant, an
   tau=1
   tau=tau/2
   an=fsum(1,p,npt,tau,deltab)
   an=an/tbeta
   do i=1,npt
      p(i)=p(i)/an
      be=deltab*(i-1)
      z=exp(be/2)
      rho=p(i)*be*(z-1/z)
   enddo

   !--calculate debye-waller lambda and effective temperature
   f0=fsum(0,p,npt,tau,deltab)
   tbar=fsum(2,p,npt,tau,deltab)/(2*tbeta)

   !--convert p(beta) into t1(beta)
   do i=1,npt
      be=deltab*(i-1)
      p(i)=p(i)*exp(be/2)/f0
   enddo

   !--print out the results
   write(nsyso,'(/'' frequency distribution'')')
   write(nsyso,'(8x,''e'',8x,''rho(e)'',6x,''beta'',&
     & 5x,''rho(beta)'',10x,''bigp'',6x,''t1(beta)'')')
   do i=1,npt
      be=deltab*(i-1)
      ee=be*tev
      z=exp(be/2)
      bigp=p(i)*f0/exp(be/2)
      rho=bigp*be*(z-1/z)
      rhoe=rho*an
      write(nsyso,'(f9.5,1p,e14.4,0p,f10.4,1p,3e14.4)')&
        ee,rhoe,be,rho,bigp,p(i)
   enddo
   dwpix(itemp)=f0
   tempf(itemp)=tbar*tempr(itemp)
   write(nsyso,'(/'' p(beta) is scaled to'',f9.5//&
     &  '' for p-bound only''/&
     &  ''             effective temp = '',f10.3/&
     &  ''        debye-waller lambda='',f10.6)')&
     tbeta,tempf(itemp),f0
   write(nsyso,'(''          factor_ph = '',f10.6)') tbar
   return
   end subroutine start

   real(kr) function fsum(n,p,np,tau,deltab)
   !--------------------------------------------------------------------
   ! Computes integrals over the phonon frequency
   ! of the form
   !    integral 0 to infinity of
   !       2*p*beta**n*hyperbolic
   !    dbeta
   ! where
   !    p is p/beta**2, or rho/(2*beta*sinh(beta/2)), and
   !    'hyperbolic' is cosh(tau*beta) for n even
   !       and sinh(tau*beta) for n odd.
   ! Called by start.
   !--------------------------------------------------------------------
   ! externals
   integer::n,np
   real(kr)::tau,deltab
   real(kr)::p(np)
   ! internals
   integer::ij
   real(kr)::arg,edsq,v,an,be,fs,w,ff

   arg=deltab*tau/2
   edsq=exp(arg)
   v=1
   an=1-2*mod(n,2)
   be=0
   fs=0
   w=1
   do ij=1,np
      if (n.gt.0) w=be**n
      ff=((p(ij)*v)*v+(p(ij)*an/v)/v)*w
      if (ij.eq.1.or.ij.eq.np) ff=ff/2
      fs=fs+ff
      be=be+deltab
      v=v*edsq
   enddo
   fsum=fs*deltab
   return
   end function fsum

   real(kr) function terpt(tn,ntn,delta,be)
   !--------------------------------------------------------------------
   ! Interpolate in a table of t_n(beta) for a required beta.
   ! Called by contin.
   !--------------------------------------------------------------------
   ! externals
   integer::ntn
   real(kr)::delta,be
   real(kr)::tn(ntn)
   ! internals
   integer::i
   real(kr)::bt,btp
   terpt=0
   if (be.gt.ntn*delta) return
   i=int(be/delta)
   if (i.lt.ntn-1) then
      bt=i*delta
      btp=bt+delta
      i=i+1
      terpt=tn(i)+(be-bt)*(tn(i+1)-tn(i))/(btp-bt)
   else
      terpt=0
   endif
   return
   end function terpt

   subroutine convol(t1,tlast,tnext,n1,nl,nn,delta,ckk)
   !--------------------------------------------------------------------
   ! Calculate the next term in the phonon expansion by
   ! convolving t1 with tlast and writing the result in tnext.
   ! The integral of tnext is also checked.
   ! Trapazoidal integration.
   ! Called by start.
   !--------------------------------------------------------------------
   ! externals
   integer::n1,nl,nn
   real(kr)::delta,ckk
   real(kr)::t1(n1),tlast(nl),tnext(nn)
   ! internals
   integer::j,k,i1,i2
   real(kr)::f1,f2,cc,be
   real(kr),parameter::tiny=1.e-30_kr
   real(kr),parameter::zero=0

   ckk=0
   do k=1,nn
      tnext(k)=0
      do j=1,n1
         i1=k+j-2
         i2=k-j
         f1=0
         be=(j-1)*delta
         if (t1(j).gt.zero) then
            if (i1+1.le.nl) f1=tlast(i1+1)*exp(-be)
            f2=0
            if (i2.ge.0.and.i2+1.le.nl) then
               f2=tlast(i2+1)
            else if (i2.lt.0.and.1-i2.le.nl) then
               be=-i2*delta
               f2=tlast(1-i2)*exp(-be)
            endif
            cc=t1(j)*(f1+f2)
            if (j.eq.1.or.j.eq.n1) cc=cc/2
            tnext(k)=tnext(k)+cc
         endif
      enddo
      tnext(k)=tnext(k)*delta
      if (tnext(k).lt.tiny) tnext(k)=0
      cc=tnext(k)
      be=(k-1)*delta
      cc=cc+tnext(k)*exp(-be)
      if (k.eq.1.or.k.eq.nn) cc=cc/2
      ckk=ckk+cc
   enddo
   ckk=ckk*delta
   return
   end subroutine convol

   subroutine trans(itemp)
   !--------------------------------------------------------------------
   ! Controls the addition of a translational contribution
   ! to a continuous S(alpha,beta).  The translational component
   ! can be either diffusion, or a free gas.  The values of the input
   ! s(alpha,beta) for the convolution are obtained by interpolation.
   ! Called by leapr.
   ! Uses stable, sbfill, terps.
   !--------------------------------------------------------------------
   use mainio ! provides nsyso
   use util   ! provides timer
   ! externals
   integer::itemp
   ! internals
   integer::ialpha,ibeta,i,nbt,iprt,jprt,nsd,nu,ndmax
   real(kr)::time,s1,s2,sum0,sum1,ff1,ff2,ff1l,ff2l
   real(kr)::sc,al,deb,ded,delta,f,s,bb,st,be,bel
   real(kr),dimension(:),allocatable::betan,ap,sd,sb
   real(kr),parameter::therm=.0253e0_kr
   real(kr),parameter::c0=.4e0_kr
   real(kr),parameter::c1=1.e0_kr
   real(kr),parameter::c2=1.42e0_kr
   real(kr),parameter::c3=.2e0_kr
   real(kr),parameter::c4=10.e0_kr
   real(kr),parameter::tiny=1.e-30_kr
   real(kr),parameter::zero=0

   !--write heading for translational calculation
   call timer(time)
   write(nsyso,'(/'' translational part of scattering law'',&
     &32x,f8.1,''s'')') time
   sc=1
   if (lat.eq.1) sc=therm/tev

   !---allocate scratch storage
   ndmax=max(nbeta,1000000)
   allocate(betan(nbeta))
   allocate(ap(ndmax))
   allocate(sd(ndmax))
   allocate(sb(ndmax))

   !--alpha loop
   if (iprint.eq.1) write(nsyso,&
     '(/'' sab checks''/''      alpha      norm   sum rule'')')
   do ialpha=1,nalpha
      iprt=mod(ialpha-1,naint)+1
      if (ialpha.eq.nalpha) iprt=1
      al=alpha(ialpha)*sc/arat
      if (iprt.eq.1.and.iprint.eq.2)&
        write(nsyso,'(/'' alpha='',f10.4)') al

      !--choose beta interval for convolution
      ded=c0*(twt*c*al)/sqrt(c1+c2*(twt*c*al)*c)
      if (ded.eq.zero) ded=c3*sqrt(twt*al)
      deb=c4*al*deltab
      delta=deb
      if (ded.lt.delta) delta=ded
      nu=1
      if (iprt.eq.1.and.iprint.eq.2) write(nsyso,&
        '(/'' delta d='',e18.5,5x,''delta b='',e18.5,&
        &10x,''delta='',e18.5)') ded,deb,delta

      !--make table of s-diffusion or s-free on this interval
      call stable(ap,sd,nsd,al,delta,iprt,nu,ndmax)
      if (nsd.gt.1) then

         !--copy original ss(-beta) to a temporary array
         do i=1,nbeta
            betan(i)=beta(i)*sc
            ap(i)=ssm(i,ialpha,itemp)
         enddo

         !--loop over beta values
         if (iprt.eq.1.and.iprint.eq.2) write(nsyso,&
           '(/'' results after convolution ''/&
           & 4x,'' beta'',7x,''s(alpha,beta)'',6x,''ss(alpha,beta)'',&
           &5x,''ss(alpha,-beta)'')')
         do ibeta=1,nbeta
            jprt=mod(ibeta-1,nbint)+1
            if (ibeta.eq.nbeta) jprt=1
            s=0
            be=betan(ibeta)

            !--prepare table of continuous ss on new interval
            nbt=nsd
            call sbfill(sb,nbt,delta,be,ap,betan,nbeta,ibeta,ndmax)

            !--convolve s-transport with s-continuous
            do i=1,nbt
               f=2*(mod(i-1,2)+1)
               if (i.eq.1.or.i.eq.nbt) f=1
               s=s+f*sd(i)*sb(nbt+i-1)
               bb=(i-1)*delta
               s=s+f*sd(i)*sb(nbt-i+1)*exp(-bb)
            enddo
            s=s*delta/3
            if (s.lt.tiny) s=0
            st=terps(sd,nbt,delta,be)
            if (st.gt.zero) s=s+exp(-al*f0)*st

            !--store results
            ssm(ibeta,ialpha,itemp)=s
            if (s.ne.zero) then
               s1=s*exp(-be/2)
               s2=s*exp(-be)
               if (iprt.eq.1.and.jprt.eq.1.and.iprint.eq.2)&
                 write(nsyso,'(f10.4,1p,e18.5,2e20.5)') be,s1,s2,s
            endif

         !--continue beta loop
         enddo

         !--check moments of calculated s(alpha,beta).
         if (iprt.eq.1) then
            sum0=0
            sum1=0
            ff1l=0
            ff2l=0
            bel=0
            do ibeta=1,nbeta
               be=betan(ibeta)
               ff2=ssm(ibeta,ialpha,itemp)
               ff1=ssm(ibeta,ialpha,itemp)*exp(-be)
               if (ibeta.gt.1) then
                  sum0=sum0+(be-bel)*(ff1l+ff2l+ff1+ff2)/2
                  sum1=sum1+(be-bel)&
                    *(ff2l*bel+ff2*be-ff1l*bel-ff1*be)/2
                  ff1l=ff1
                  ff2l=ff2
                  bel=be
               else
                  bel=be
                  ff1l=ff1
                  ff2l=ff2
                  sum0=0
                  sum1=0
               endif
            enddo
            sum1=sum1/al/(tbeta+twt)
            if (iprint.eq.2) then
               write(nsyso,'(&
                 &''     normalization check ='',f8.4)') sum0
               write(nsyso,'(&
                 &''          sum rule check ='',f8.4)') sum1
            else if (iprint.eq.1) then
               write(nsyso,'(1x,f10.4,2f10.4)') al,sum0,sum1
            endif
         endif
      endif

   !--continue alpha loop
   enddo

   !--update effective temperature
   tempf(itemp)=(tbeta*tempf(itemp)+twt*tempr(itemp))/(tbeta+twt)
   write(nsyso,'(/''     new effective temp = '',f10.3)') tempf(itemp)

   !--deallocate scratch storage
   deallocate(sb)
   deallocate(sd)
   deallocate(ap)
   deallocate(betan)
   return
   end subroutine trans

   subroutine stable(ap,sd,nsd,al,delta,iprt,nu,ndmax)
   !--------------------------------------------------------------------
   ! Sets up table of S-diffusion or S-free in the array sd,
   ! evaluated at intervals delta determined by trans.
   ! Tabulation is continued until
   !       sd(j) is less than 1e-7*sd(1)
   !       or nsd is 1999
   ! Here, nsd is always odd for use with Simpson's rule, and
   ! ap is used for temporary storage of beta values.
   ! This routine returns the -beta side of
   ! the asymmetric S(alpha,beta).
   ! Called by trans.
   ! Uses besk1.
   !--------------------------------------------------------------------
   use mainio ! provide nsyso
   use physics ! provides pi
   ! externals
   integer::nsd,iprt,nu,ndmax
   real(kr)::al,delta,check0,check1,f,bb,sfree
   real(kr)::ap(ndmax),sd(ndmax)
   ! internals
   integer::i,j,idone,icheck
   real(kr)::d,c2,c3,c4,c5,c6,c7,c8,be,ex,wal
   real(kr),parameter::quart=0.25e0_kr
   real(kr),parameter::eps=1.e-7_kr
   real(kr),parameter::zero=0
   real(kr),parameter::one=1
   icheck=0

   !--diffusion branch
   if (c.ne.zero) then
      if (iprt.eq.1.and.iprint.eq.2) write(nsyso,&
        '(/4(4x,''beta'',4x,''ssdiff(-beta)''))')
      d=twt*c
      c2=sqrt(c*c+quart)
      c3=2*d*al
      c4=c3*c3
      c8=c2*c3/pi
      c3=2*d*c*al
      be=0
      j=1
      idone=0
      do while (idone.eq.0)
         c6=sqrt(be*be+c4)
         c7=c6*c2
         if (c7.le.one) c5=c8*exp(c3+be/2)
         if (c7.gt.one) then
            c5=0
            ex=c3-c7+be/2
            c5=c8*exp(ex)
         endif
         sd(j)=c5*besk1(c7)/c6
         ap(j)=be
         be=be+delta
         j=j+1
         if (mod(j,2).eq.0) then
            if (j.ge.ndmax) idone=1
            if (eps*sd(1).ge.sd(j-1)) idone=1
         endif
      enddo
      j=j-1
      nsd=j
      if (iprt.eq.1.and.iprint.eq.2)&
        write(nsyso,'(0pf10.5,1pe15.7,0pf10.5,1pe15.7,0pf10.5,&
        &1pe15.7,0pf10.5,1pe15.7)') (ap(i),sd(i),i=1,j,nu)

   !--free-gas branch
   else
      if (iprt.eq.1.and.iprint.eq.2) write(nsyso,&
        '(/4(4x,''beta'',4x,''ssfree(-beta)''))')
      be=0
      j=1
      wal=twt*al
      idone=0
      do while (idone.eq.0)
         ex=-(wal-be)**2/(4*wal)
         sfree=exp(ex)/sqrt(4*pi*wal)
         sd(j)=sfree
         ap(j)=be
         be=be+delta
         j=j+1
         if (mod(j,2).eq.0) then
            if (j.ge.ndmax) idone=1
            if (eps*sd(1).ge.sd(j-1)) idone=1
         endif
      enddo
      j=j-1
      nsd=j
      if (iprt.eq.1.and.iprint.eq.2)&
        write(nsyso,'(0pf10.5,1pe15.7,0pf10.5,1pe15.7,0pf10.5,&
        &1pe15.7,0pf10.5,1pe15.7)') (ap(i),sd(i),i=1,j,nu)
   endif

   !--check the moments of the distribution
   if (iprt.ne.1.or.iprint.ne.2) return
   if (icheck.eq.0) return
   check0=0
   check1=0
   do i=1,nsd
      f=2*(mod(i-1,2)+1)
      if (i.eq.1.or.i.eq.nsd) f=1
      bb=(i-1)*delta
      check0=check0+f*sd(i)
      check0=check0+sd(i)*exp(-bb)
      check1=check1+f*bb*sd(i)
      check1=check1+f*bb*sd(i)*exp(-bb)
   enddo
   check0=check0*delta/3
   check1=check1*delta/3
   check1=check1/(al*twt)
   write(nsyso,'(/'' check0='',f11.7,''   check1='',f11.7)')&
     check0,check1
   return
   end subroutine stable

   real(kr) function terps(sd,nsd,delta,be)
   !--------------------------------------------------------------------
   ! Interpolate in a table of S(alpha,beta) for a required beta.
   ! Used in trans.
   !--------------------------------------------------------------------
   ! externals
   integer::nsd
   real(kr)::delta,be
   real(kr)::sd(nsd)
   ! internals
   integer::i
   real(kr)::bt,btp,st,stp,stt
   real(kr),parameter::slim=-225.e0_kr
   real(kr),parameter::zero=0

   terps=0
   if (be.gt.delta*nsd) return
   i=int(be/delta)
   if (i.lt.nsd-1) then
      bt=i*delta
      btp=bt+delta
      i=i+1
      if (sd(i).le.zero) then
         st=slim
      else
         st=log(sd(i))
      endif
      if (sd(i+1).le.zero) then
         stp=slim
      else
         stp=log(sd(i+1))
      endif
      stt=st+(be-bt)*(stp-st)/(btp-bt)
      terps=0
      if (stt.gt.slim) terps=exp(stt)
      return
   endif
   return
   end function terps

   subroutine sbfill(sb,nbt,delta,be,s,betan,nbeta,nbe,ndmax)
   !--------------------------------------------------------------------
   ! For translational cases only.
   ! Generates s(beta) on a new energy grid for convolution
   ! with a diffusion or free-gas shape.  Interpolation is used.
   ! Called by trans.
   !--------------------------------------------------------------------
   use util    ! provides error
   ! externals
   integer::nbt,nbeta,nbe,ndmax
   real(kr)::delta,be
   real(kr)::sb(ndmax),s(ndmax),betan(nbeta)
   ! internals
   character(60)::strng
   integer::i,j,idone
   real(kr)::bmin,bmax,b,st,stm,arg,bet
   real(kr),parameter::shade=1.00001e0_kr
   real(kr),parameter::slim=-225.e0_kr
   real(kr),parameter::zero=0

   bmin=-be-(nbt-1)*delta
   bmax=-be+(nbt-1)*delta+delta/100
   if ((1+int((bmax-bmin)/delta)).gt.ndmax) then
      write(strng,'(''ndmax needs to be at least '',i6)')&
                    1+int((bmax-bmin)/delta)
      call error('sbfill',strng,' ')
   endif
   j=nbeta
   i=0
   bet=bmin
   do while (bet.le.bmax)
      i=i+1
      b=abs(bet)
      ! search for correct beta range
      idone=0
      do while (idone.eq.0)
         if (b.gt.betan(j)) then
            if (j.eq.nbeta.and.b.lt.shade*betan(j)) then
               idone=1
            else
               if (j.eq.nbeta) then
                  idone=2
               else
                  ! move up
                  j=j+1
               endif
            endif
         else
            if (b.gt.betan(j-1)) then
               idone=1
            else
               if (j.eq.2) then
                  idone=1
               else
                  ! move down
                  j=j-1
               endif
            endif
         endif
      enddo
      ! interpolate in this range
      if (idone.eq.1) then
         if (s(j).le.zero) then
            st=slim
         else
            st=log(s(j))
         endif
         if (s(j-1).le.zero) then
            stm=slim
         else
            stm=log(s(j-1))
         endif
         sb(i)=st+(b-betan(j))*(stm-st)/(betan(j-1)-betan(j))
         if (bet.gt.zero) sb(i)=sb(i)-bet
         arg=sb(i)
         sb(i)=0
         if (arg.gt.slim) sb(i)=exp(arg)
      else
         sb(i)=0
      endif
      ! if delta is to small for the current value of beta, increase it
      do while (bet.eq.(bet+delta))
         delta=delta*10
      end do
      bet=bet+delta
   enddo
   return
   end subroutine sbfill

   real(kr) function besk1(x)
   !--------------------------------------------------------------------
   ! Computes modified Bessel function, K1.
   ! The exponential part for x>1 is omitted (see stable).
   ! Called by stable.
   !--------------------------------------------------------------------
   ! externals
   real(kr)::x
   ! internals
   real(kr)::test,v,u,bi1,bi3
   real(kr),parameter::c0=.125e0_kr
   real(kr),parameter::c1=.442850424e0_kr
   real(kr),parameter::c2=.584115288e0_kr
   real(kr),parameter::c3=6.070134559e0_kr
   real(kr),parameter::c4=17.864913364e0_kr
   real(kr),parameter::c5=48.858995315e0_kr
   real(kr),parameter::c6=90.924600045e0_kr
   real(kr),parameter::c7=113.795967431e0_kr
   real(kr),parameter::c8=85.331474517e0_kr
   real(kr),parameter::c9=32.00008698e0_kr
   real(kr),parameter::c10=3.999998802e0_kr
   real(kr),parameter::c11=1.304923514e0_kr
   real(kr),parameter::c12=1.47785657e0_kr
   real(kr),parameter::c13=16.402802501e0_kr
   real(kr),parameter::c14=44.732901977e0_kr
   real(kr),parameter::c15=115.837493464e0_kr
   real(kr),parameter::c16=198.437197312e0_kr
   real(kr),parameter::c17=222.869709703e0_kr
   real(kr),parameter::c18=142.216613971e0_kr
   real(kr),parameter::c19=40.000262262e0_kr
   real(kr),parameter::c20=1.999996391e0_kr
   real(kr),parameter::c21=1.e0_kr
   real(kr),parameter::c22=.5e0_kr
   real(kr),parameter::c23=.5772156649e0_kr
   real(kr),parameter::c24=1.e0_kr
   real(kr),parameter::c25=.0108241775e0_kr
   real(kr),parameter::c26=.0788000118e0_kr
   real(kr),parameter::c27=.2581303765e0_kr
   real(kr),parameter::c28=.5050238576e0_kr
   real(kr),parameter::c29=.663229543e0_kr
   real(kr),parameter::c30=.6283380681e0_kr
   real(kr),parameter::c31=.4594342117e0_kr
   real(kr),parameter::c32=.2847618149e0_kr
   real(kr),parameter::c33=.1736431637e0_kr
   real(kr),parameter::c34=.1280426636e0_kr
   real(kr),parameter::c35=.1468582957e0_kr
   real(kr),parameter::c36=.4699927013e0_kr
   real(kr),parameter::c37=1.2533141373e0_kr

   test=1
   if (x.le.test) then
      v=c0*x
      u=v*v
      bi1=(((((((((c1*u+c2)*u+c3)*u+c4)*u+c5)*u+c6)*u+c7)*u&
        +c8)*u+c9)*u+c10)*v
      bi3=(((((((((c11*u+c12)*u+c13)*u+c14)*u+c15)*u+c16)*u&
        +c17)*u+c18)*u+c19)*u+c20)
      besk1=c21/x+bi1*(log(c22*x)+c23)-v*bi3
   else
      u=c24/x
      bi3=((((((((((((-c25*u+c26)*u-c27)*u+c28)*u-c29)*u+c30)*u&
        -c31)*u+c32)*u-c33)*u+c34)*u-c35)*u+c36)*u+c37)
      besk1=sqrt(u)*bi3
   endif
   return
   end function besk1

   subroutine discre(itemp)
   !--------------------------------------------------------------------
   ! Controls the convolution of discrete oscillators with
   ! the continuous S(alpha,beta) computed in contin.
   ! Called by leapr.
   ! Uses bfact, bfill, exts, sint.
   !--------------------------------------------------------------------
   use mainio  ! provides nsyso
   use util    ! provides timer
   use physics ! provides bk (boltzmann constant)
   ! externals
   integer::itemp
   ! internals
   real(kr)::time,ss,s1,s2,sc,sumn,sumr,st,add,besn,wtsn
   real(kr)::sum0,sum1,bel,ff1,ff2,ff1l,ff2l
   real(kr)::x,db,save,dwc,dwf,al,tbart,wt,be,cn,sn
   real(kr)::tsave,dw0,dwt
   integer::nal,iprt,i,j,k,m,n,ibeta,idone
   integer::jj,nn,jprt,nbx,maxbb,maxdd
   real(kr)::bdeln(50),eb(50),dbw(50),ar(50),dist(50)
   real(kr)::bzero,bplus(50),bminus(50)
   real(kr),dimension(:),allocatable::betan,exb,sexpb
   real(kr),dimension(:),allocatable::bex,rdbex,sex
   real(kr),dimension(:),allocatable::bes,wts,ben,wtn
   real(kr),parameter::therm=.0253e0_kr
   real(kr),parameter::small=1.e-8_kr
   real(kr),parameter::vsmall=1.e-10_kr
   real(kr),parameter::tiny=1.e-20_kr
   real(kr),parameter::zero=0

   !--write heading
   call timer(time)
   write(nsyso,&
     '(/'' discrete-oscillator part of scattering law'',&
     &26x,f8.1,''s'')') time
   sc=1
   if (lat.eq.1) sc=therm/tev

   !--allocate scratch storage
   allocate(betan(nbeta))
   allocate(exb(nbeta))
   allocate(sexpb(nbeta))
   maxbb=2*nbeta+1
   allocate(bex(maxbb))
   allocate(rdbex(maxbb))
   allocate(sex(maxbb))
   maxdd=500
   allocate(bes(maxdd))
   allocate(wts(maxdd))
   allocate(ben(maxdd))
   allocate(wtn(maxdd))

   !--set up oscillator parameters
   dwt=0
   do i=1,nd
      bdeln(i)=bdel(i)/tev
      dwt=dwt+adel(i)
   enddo
   tsave=0
   dw0=dwpix(itemp)
   do i=1,nd
      eb(i)=exp(bdeln(i)/2)
      sn=(eb(i)-1/eb(i))/2
      cn=(eb(i)+1/eb(i))/2
      ar(i)=adel(i)/(sn*bdeln(i))
      dist(i)=adel(i)*bdel(i)*cn/(2*sn)
      tsave=tsave+dist(i)/bk
      dbw(i)=ar(i)*cn
      if (dwpix(itemp).gt.zero) dwpix(itemp)=dwpix(itemp)+dbw(i)
   enddo
   write(nsyso,'(/'' add delta functions''//(5x,i3,1p,2e14.4))')&
     (i,bdeln(i),adel(i),i=1,nd)

   !--prepare functions of beta
   do i=1,nbeta
      be=beta(i)*sc
      exb(i)=exp(-be/2)
      betan(i)=be
   enddo
   call bfill(bex,rdbex,nbx,betan,nbeta,maxbb)
   wt=tbeta
   tbart=tempf(itemp)/tempr(itemp)

   !--main alpha loop
   if (iprint.eq.1) write(nsyso,&
     '(/'' sab checks''/''      alpha      norm   sum rule'')')
   do nal=1,nalpha
      iprt=mod(nal-1,naint)+1
      if (nal.eq.nalpha) iprt=1
      al=alpha(nal)*sc/arat
      if (iprt.eq.1.and.iprint.eq.2) write(nsyso,&
        '(/3x,''alpha='',f10.5)') al
      dwf=exp(-al*dw0)
      if (iprt.eq.1.and.iprint.eq.2) write(nsyso,&
        '(/''      debye-waller factor='',1p,e12.4)') dwf
      call exts(ssm(1,nal,itemp),sex,exb,betan,nbeta,maxbb)
      do j=1,nbeta
         sexpb(j)=0
      enddo

      !---initialize for delta function calculation
      ben(1)=0
      wtn(1)=1
      nn=1
      n=0

      !--loop over all oscillators
      do i=1,nd
         dwc=al*dbw(i)
         x=al*ar(i)
         call bfact(x,bzero,bplus,bminus,dwc,bdeln(i))

         !--do convolution for the delta functions
         !--n=0 term
         do m=1,nn
            besn=ben(m)
            wtsn=wtn(m)*bzero
            if (besn.le.zero.or.wtsn.ge.small) then
               if (n.lt.maxdd) then
                  n=n+1
                  bes(n)=besn
                  wts(n)=wtsn
               endif
            endif
         enddo

         !--negative n terms
         k=0
         idone=0
         do while (k.lt.50.and.idone.eq.0)
            k=k+1
            if (bminus(k).le.zero) then
               idone=1
            else
               do m=1,nn
                  besn=ben(m)-k*bdeln(i)
                  wtsn=wtn(m)*bminus(k)
                  if (wtsn.ge.small.and.n.lt.maxdd) then
                     n=n+1
                     bes(n)=besn
                     wts(n)=wtsn
                  endif
               enddo
            endif
         enddo

         !--positive n terms
         k=0
         idone=0
         do while (k.lt.50.and.idone.eq.0)
            k=k+1
            if (bplus(k).le.zero) then
               idone=1
            else
               do m=1,nn
                  besn=ben(m)+k*bdeln(i)
                  wtsn=wtn(m)*bplus(k)
                  if (wtsn.ge.small.and.n.lt.maxdd) then
                     n=n+1
                     bes(n)=besn
                     wts(n)=wtsn
                  endif
               enddo
            endif
         enddo

         !--continue oscillator loop
         nn=n
         do m=1,nn
            ben(m)=bes(m)
            wtn(m)=wts(m)
         enddo
         n=0
         wt=wt+adel(i)
         tbart=tbart+dist(i)/bk/tempr(itemp)
      enddo
      n=nn

      !--sort the discrete lines
      !--and throw out the smallest ones
      nn=n-1
      do i=2,n-1
         do j=i+1,n
            if (wts(j).ge.wts(i)) then
               save=wts(j)
               wts(j)=wts(i)
               wts(i)=save
               save=bes(j)
               bes(j)=bes(i)
               bes(i)=save
            endif
         enddo
      enddo
      i=0
      idone=0
      do while (i.lt.nn.and.idone.eq.0)
         i=i+1
         n=i
         if (wts(i).lt.100*small.and.i.gt.5) idone=1
      enddo

      !--report the discrete lines
      if (iprint.ge.2) then
         write(nsyso,'(/''      discrete lines''/&
           &14x,''beta'',6x,''weight'')')
         sumn=0
         sumr=0
         do m=1,n
            write(nsyso,'(6x,f12.4,1p,e12.4)') bes(m),wts(m)
            sumn=sumn+wts(m)
            sumr=sumr-bes(m)*wts(m)
         enddo
         sumr=sumr/al/dwt
         write(nsyso,'(8x,''norm check'',f10.4/&
           &8x,''rule check'',f10.4)')&
           sumn,sumr
      endif

      !--add the continuum part to the scattering law
      do m=1,n
         do j=1,nbeta
            be=-betan(j)-bes(m)
            st=sint(be,bex,rdbex,sex,nbx,al,tbeta+twt,tbart,betan,&
              nbeta,maxbb)
            add=wts(m)*st
            if (add.ge.tiny) sexpb(j)=sexpb(j)+add
         enddo
      enddo

      !--add the delta functions to the scattering law
      !--delta(0.) is saved for the incoherent elastic
      if (twt.le.zero) then
         m=0
         idone=0
         do while (m.lt.n.and.idone.eq.0)
            m=m+1
            if (dwf.lt.vsmall) then
               idone=1
            else
               if (bes(m).lt.zero) then
                  be=-bes(m)
                  if (be.le.betan(nbeta-1)) then
                     db=1000
                     idone=0
                     j=0
                     do while (j.lt.nbeta.and.idone.eq.0)
                        j=j+1
                        jj=j
                        if (abs(be-betan(j)).gt.db) then
                           idone=1
                        else
                           db=abs(be-betan(j))
                        endif
                     enddo
                     if (jj.le.2) then
                        add=wts(m)/betan(jj)
                     else
                        add=2*wts(m)/(betan(jj)-betan(jj-2))
                     endif
                     add=add*dwf
                     if (add.ge.tiny) sexpb(jj-1)=sexpb(jj-1)+add
                  endif
               endif
            endif
         enddo
      endif

      !--record the results
      do j=1,nbeta
         ssm(j,nal,itemp)=sexpb(j)
      enddo
      if (iprt.eq.1.and.iprint.eq.2) write(nsyso,&
        '(/4x,'' beta'',7x,''s(alpha,beta)'',7x,''ss(alpha,beta)'',&
        &5x,''ss(alpha,-beta)'')')
      do i=1,nbeta
         be=beta(i)*sc
         ss=ssm(i,nal,itemp)
         s1=ss*exp(-be/2)
         s2=ss*exp(-be)
         jprt=mod(i-1,nbint)+1
         if (i.eq.nbeta) jprt=1
         if (iprt.eq.1.and.jprt.eq.1.and.iprint.eq.2)&
           write(nsyso,'(f10.4,1pe18.5,1p,2e20.5)') betan(i),s1,s2,ss
         enddo

      !--check moments of calculated s(alpha,beta).
      if (iprt.eq.1) then
         sum0=0
         sum1=0
         ff1l=0
         ff2l=0
         bel=0
         do ibeta=1,nbeta
            be=betan(ibeta)
            ff2=ssm(ibeta,nal,itemp)
            ff1=ssm(ibeta,nal,itemp)*exp(-be)
            if (ibeta.gt.1) then
               sum0=sum0+(be-bel)*(ff1l+ff2l+ff1+ff2)/2
               sum1=sum1+(be-bel)*(ff2l*bel+ff2*be-ff1l*bel-ff1*be)/2
               ff1l=ff1
               ff2l=ff2
               bel=be
            else
               bel=be
               ff1l=ff1
               ff2l=ff2
               sum0=0
               sum1=0
            endif
         enddo
         if (twt.eq.zero) sum0=sum0/(1-exp(-al*dwpix(itemp)))
         sum1=sum1/al
         if (iprint.eq.2) then
            write(nsyso,'(''     normalization check ='',f8.4)') sum0
            write(nsyso,'(''          sum rule check ='',f8.4)') sum1
         else if (iprint.eq.1) then
            write(nsyso,'(1x,f10.4,2f10.4)') al,sum0,sum1
         endif
      endif

   !--continue the alpha loop
   enddo

   !--finished
   tempf(itemp)=(tbeta+twt)*tempf(itemp)+tsave
   write(nsyso,'(/&
     &  ''       new effective temp = '',f10.3/&
     &  ''  new debye-waller lambda='',f10.6)')&
     tempf(itemp),dwpix(itemp)
   write(nsyso,'('' discr.-oscill. part of eff. temp = '',f10.3)') tsave
   deallocate(wtn)
   deallocate(ben)
   deallocate(wts)
   deallocate(bes)
   deallocate(sex)
   deallocate(rdbex)
   deallocate(bex)
   deallocate(sexpb)
   deallocate(exb)
   deallocate(betan)
   return
   end subroutine discre

   subroutine bfact(x,bzero,bplus,bminus,dwc,betai)
   !--------------------------------------------------------------------
   ! Calculates the Bessel function terms for discrete oscillators.
   ! Called by discre.
   !--------------------------------------------------------------------
   ! externals
   real(kr)::x,dwc,betai
   real(kr)::bzero,bplus(50),bminus(50)
   ! internals
   integer::i,imax,j
   real(kr)::bn(50)
   real(kr)::y,u,v,bessi0,bessi1,rat,arg
   real(kr),parameter::c0=3.75e0_kr
   real(kr),parameter::c1=1.e0_kr
   real(kr),parameter::c2=3.5156229e0_kr
   real(kr),parameter::c3=3.0899424e0_kr
   real(kr),parameter::c4=1.2067492e0_kr
   real(kr),parameter::c5=0.2659732e0_kr
   real(kr),parameter::c6=0.0360768e0_kr
   real(kr),parameter::c7=0.0045813e0_kr
   real(kr),parameter::c8=0.39894228e0_kr
   real(kr),parameter::c9=0.01328592e0_kr
   real(kr),parameter::c10=0.00225319e0_kr
   real(kr),parameter::c11=0.00157565e0_kr
   real(kr),parameter::c12=0.00916281e0_kr
   real(kr),parameter::c13=0.02057706e0_kr
   real(kr),parameter::c14=0.02635537e0_kr
   real(kr),parameter::c15=0.01647633e0_kr
   real(kr),parameter::c16=0.00392377e0_kr
   real(kr),parameter::c17=0.5e0_kr
   real(kr),parameter::c18=0.87890594e0_kr
   real(kr),parameter::c19=0.51498869e0_kr
   real(kr),parameter::c20=0.15084934e0_kr
   real(kr),parameter::c21=0.02658733e0_kr
   real(kr),parameter::c22=0.00301532e0_kr
   real(kr),parameter::c23=0.00032411e0_kr
   real(kr),parameter::c24=0.02282967e0_kr
   real(kr),parameter::c25=0.02895312e0_kr
   real(kr),parameter::c26=0.01787654e0_kr
   real(kr),parameter::c27=0.00420059e0_kr
   real(kr),parameter::c28=0.39894228e0_kr
   real(kr),parameter::c29=0.03988024e0_kr
   real(kr),parameter::c30=0.00362018e0_kr
   real(kr),parameter::c31=0.00163801e0_kr
   real(kr),parameter::c32=0.01031555e0_kr
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::tiny=1.e-30_kr
   real(kr),parameter::zero=0
   real(kr),parameter::one=1

   !--compute bessi0
   y=x/c0
   if (y.le.one) then
      u=y*y
      bessi0=c1+u*(c2+u*(c3+u*(c4+u*(c5+u*(c6+u*c7)))))
   else
      v=1/y
      bessi0=(c8+v*(c9+v*(c10+v*(-c11+v*(c12+v*(-c13&
        +v*(c14+v*(-c15+v*c16))))))))/sqrt(x)
   endif

   !--compute bessi1
   if (y.le.one) then
      u=y*y
      bessi1=(c17+u*(c18+u*(c19+u*(c20+u*(c21+u*(c22+u*c23))))))*x
   else
      v=1/y
      bessi1=c24+v*(-c25+v*(c26-v*c27))
      bessi1=c28+v*(-c29+v*(-c30+v*(c31+v*(-c32+v*bessi1))))
      bessi1=bessi1/sqrt(x)
   endif

   !--generate higher orders by reverse recursion
   imax=50
   bn(imax)=0
   bn(imax-1)=1
   i=imax-1
   do while (i.gt.1)
      bn(i-1)=bn(i+1)+i*(2/x)*bn(i)
      i=i-1
      if (bn(i).ge.big) then
         do j=i,imax
            bn(j)=bn(j)/big
         enddo
      endif
   enddo
   rat=bessi1/bn(1)
   do i=1,imax
      bn(i)=bn(i)*rat
      if (bn(i).lt.tiny) bn(i)=0
   enddo

   !--apply exponential terms to bessel functions
   if (y.le.one) then
      bzero=bessi0*exp(-dwc)
      do i=1,imax
         bminus(i)=0
         bplus(i)=0
         if (bn(i).ne.zero) then
            arg=-dwc-i*betai/2
            bplus(i)=0
            bplus(i)=exp(arg)*bn(i)
            if (bplus(i).lt.tiny) bplus(i)=0
            bminus(i)=0
            arg=-dwc+i*betai/2
            bminus(i)=exp(arg)*bn(i)
            if (bminus(i).lt.tiny) bminus(i)=0
         else
            bplus(i)=0
            bminus(i)=0
         endif
      enddo
   else
      bzero=bessi0*exp(-dwc+x)
      do i=1,imax
         bminus(i)=0
         bplus(i)=0
         if (bn(i).ne.zero) then
            bplus(i)=0
            arg=-dwc-i*betai/2+x
            bplus(i)=exp(arg)*bn(i)
            if (bplus(i).lt.tiny) bplus(i)=0
            bminus(i)=0
            arg=-dwc+i*betai/2+x
            bminus(i)=exp(arg)*bn(i)
            if (bminus(i).lt.tiny) bminus(i)=0
         else
            bplus(i)=0
            bminus(i)=0
         endif
      enddo
   endif
   return
   end subroutine bfact

   subroutine bfill(bex,rdbex,nbx,betan,nbetan,maxbb)
   !--------------------------------------------------------------------
   ! Sets up the arrays bex and rdbex used by sint.
   ! Called by discre.
   !--------------------------------------------------------------------
   ! externals
   integer::nbx,nbetan,maxbb
   real(kr)::bex(maxbb),rdbex(maxbb),betan(nbetan)
   ! internals
   integer::i,k,nbxm
   real(kr),parameter::small=1.e-9_kr

   k=nbeta
   do i=1,nbeta
      bex(i)=-betan(k)
      k=k-1
   enddo
   if (betan(1).le.small) then
      bex(nbeta)=0
      k=nbeta+1
   else
      k=nbeta+2
      bex(nbeta+1)=betan(1)
   endif
   do i=2,nbeta
      bex(k)=betan(i)
      k=k+1
   enddo
   nbxm=k-2
   nbx=k-1
   do i=1,nbxm
      rdbex(i)=1/(bex(i+1)-bex(i))
   enddo
   return
   end subroutine bfill

   subroutine exts(sexpb,sex,exb,betan,nbetan,maxbb)
   !--------------------------------------------------------------------
   ! Sets up the array sex for sint.
   ! Here sexpb contains the asymmetric sab for negative beta, and
   ! sex contains the asymmetric sab extended to plus and minus beta.
   ! Called by discre.
   !--------------------------------------------------------------------
   ! externals
   integer::nbetan,maxbb
   real(kr)::sexpb(nbetan),sex(maxbb),exb(maxbb),betan(nbetan)
   ! internals
   integer::i,k
   real(kr),parameter::small=1.e-9_kr

   k=nbeta
   do i=1,nbeta
      sex(i)=sexpb(k)
      k=k-1
   enddo
   if (betan(1).le.small) then
      sex(nbeta)=sexpb(1)
      k=nbeta+1
   else
      k=nbeta+2
      sex(nbeta+1)=sexpb(1)
   endif
   do i=2,nbeta
      sex(k)=sexpb(i)*exb(i)*exb(i)
      k=k+1
   enddo
   return
   end subroutine exts

   real(kr) function sint(x,bex,rdbex,sex,nbx,alph,wt,tbart,&
     betan,nbetan,maxbb)
   !--------------------------------------------------------------------
   ! Interpolates in scattering function, or uses SCT approximation
   ! to extrapolate outside the range in memory.
   ! Called by discre.
   !--------------------------------------------------------------------
   use physics ! provides pi
   ! externals
   integer::nbx,nbetan,maxbb
   real(kr)::x,alph,wt,tbart
   real(kr)::bex(maxbb),rdbex(maxbb),sex(maxbb),betan(nbetan)
   ! internals
   integer::k1,k2,k3,idone
   real(kr)::sv,ex,ss1,ss3
   real(kr),parameter::slim=-225.e0_kr
   real(kr),parameter::zero=0

   !--sct approx
   if (abs(x).gt.betan(nbeta)) then
      if (alph.le.zero) then
         sv=0
      else
         ex=-(wt*alph-abs(x))**2/(4*wt*alph*tbart)
         if (x.gt.zero) ex=ex-x
         sv=exp(ex)/(4*pi*wt*alph*tbart)
      endif
      sint=sv
      return
   endif

   !--interpolation
   k1=1
   k2=nbeta
   k3=nbx
   !  bisect for x
   idone=0
   do while (idone.eq.0)
      if (x.eq.bex(k2)) then
         sv=sex(k2)
         sint=sv
         return
      else if (x.gt.bex(k2)) then
         k1=k2
         k2=(k3-k2)/2+k2
         if (k3-k1.le.1) idone=1
      else
         k3=k2
         k2=(k2-k1)/2+k1
         if (k3-k1.le.1) idone=1
      endif
   enddo
   if (sex(k1).le.zero) then
      ss1=slim
   else
      ss1=log(sex(k1))
   endif
   if (sex(k3).le.zero) then
      ss3=slim
   else
      ss3=log(sex(k3))
   endif
   ex=((bex(k3)-x)*ss1+(x-bex(k1))*ss3)*rdbex(k1)
   sv=0
   if (ex.gt.slim) sv=exp(ex)
   sint=sv
   return
   end function sint

   subroutine coldh(itemp,temp)
   !--------------------------------------------------------------------
   ! Convolve the current solid-type and/or diffusive S(alpha,beta)
   ! with discrete rotational modes for ortho or para hydrogen or
   ! deuterium.   The discrete modes are computed using the formulas
   ! of Young and Koppel for the vibrational ground state with
   ! coding based on contributions from Robert (Grenoble) and
   ! Neef (Julich).  The approach of using solid/diffusive modes
   ! with discrete rotations is based on the work of Keinert and
   ! Sax.  Note that the final S(alpha,beta) is not symmetric in beta.
   !--------------------------------------------------------------------
   use physics ! provides pi,bk,hbar,ev
   use mainio  ! provides nsyso
   use util    ! provides timer
   ! externals
   integer::itemp
   real(kr)::temp
   ! internals
   real(kr)::time,tev,sc,de,x,amassm,bp,sampc,sampi
   real(kr)::snlg,betap,bn,ex,add,snlk,up,down,sn,snorm
   real(kr)::sum0,bel,ff1,ff2,ff1l,ff2l,tmp,total,be
   real(kr)::al,alp,waven,y,sk,swe,swo,wt,tbart,pj
   integer::i,j,k,l,jj,jjmax,jprt,nbx,maxbb
   integer::law,nal,iprt,ipo,jt1,lp,jp,nbe
   real(kr),dimension(:),allocatable::betan,exb
   real(kr),dimension(:),allocatable::bex,rdbex,sex
   real(kr),parameter::pmass=1.6726231e-24_kr
   real(kr),parameter::dmass=3.343586e-24_kr
   real(kr),parameter::deh=0.0147e0_kr
   real(kr),parameter::ded=0.0074e0_kr
   real(kr),parameter::sampch=0.356e0_kr
   real(kr),parameter::sampcd=0.668e0_kr
   real(kr),parameter::sampih=2.526e0_kr
   real(kr),parameter::sampid=0.403e0_kr
   real(kr),parameter::small=1.e-6_kr
   real(kr),parameter::therm=.0253e0_kr
   real(kr),parameter::angst=1.e-8_kr
   integer::ifree=0
   integer::nokap=0
   integer::jterm=3
   real(kr),parameter::zero=0

   !--write header
   call timer(time)
   write(nsyso,'(/'' cold hydrogen or deuterium scattering'',&
     &31x,f8.1,''s'')') time

   !--allocate scratch storage
   allocate(betan(nbeta))
   allocate(exb(nbeta))
   maxbb=2*nbeta+1
   allocate(bex(maxbb))
   allocate(rdbex(maxbb))
   allocate(sex(maxbb))

   !--set up constants
   tev=bk*abs(temp)
   sc=1
   if (lat.eq.1) sc=therm/tev
   law=ncold+1
   de=deh
   if (law.gt.3) de=ded
   x=de/tev
   if (law.gt.3) then
     amassm=6.69E-24_kr
     ! amassm=2*(amassd+amasse)*amu*ev/(clight*clight)
     sampc=sampcd
     bp=hbar/2*sqrt(2/ded/ev/dmass)/angst
     sampi=sampid
   else
     amassm=3.3464e-24_kr
     ! amassm=2*(amassp+amasse)*amu*ev/(clight*clight)
     sampc=sampch
     bp=hbar/2*sqrt(2/deh/ev/pmass)/angst
     sampi=sampih
   endif
   wt=twt+tbeta
   tbart=tempf(itemp)/tempr(itemp)

   !--main alpha loop
   do nal=1,nalpha
      iprt=mod(nal-1,naint)+1
      if (nal.eq.nalpha) iprt=1
      al=alpha(nal)*sc/arat
      alp=wt*al
      waven=angst*sqrt(amassm*tev*ev*al)/hbar
      y=bp*waven
      if (iprt.eq.1) write(nsyso,'(//i4,3x,''alpha='',f10.5)') nal,al
      sk=terpk(ska,nka,dka,waven)
      if (nokap.eq.1) sk=1
      if (iprt.eq.1.and.iprint.eq.2) then
         write(nsyso,'(/'' wave number             ='',f10.4)') waven
         write(nsyso,'('' static structure factor ='',f10.4)') sk
         write(nsyso,'('' oscillators:'')')
      endif

      !--spin-correlation factors
      if (law.eq.2) swe=sampi**2/3
      if (law.eq.2) swo=sk*sampc**2+2*sampi**2/3
      if (law.eq.3) swe=sk*sampc**2
      if (law.eq.3) swo=sampi**2
      if (law.eq.4) swe=sk*sampc**2+5*sampi**2/8
      if (law.eq.4) swo=3*sampi**2/8
      if (law.eq.5) swe=3*sampi**2/4
      if (law.eq.5) swo=sk*sampc**2+sampi**2/4
      snorm=sampi**2+sampc**2
      swe=swe/snorm
      swo=swo/snorm

      !--prepare arrays for sint
      if (nal.eq.1) then
         do i=1,nbeta
            be=beta(i)
            if (lat.eq.1) be=be*therm/tev
            exb(i)=exp(-be/2)
            betan(i)=be
         enddo
         call bfill(bex,rdbex,nbx,betan,nbeta,maxbb)
      endif
      call exts(ssm(1,nal,itemp),sex,exb,betan,nbeta,maxbb)

      !--loop over all beta values
      !    results for positive beta go into ssp
      !    results for negative beta go into ssm
      jjmax=2*nbeta-1
      do jj=1,jjmax
         if (jj.lt.nbeta) k=nbeta-jj+1
         if (jj.ge.nbeta) k=jj-nbeta+1
         be=betan(k)
         if (jj.lt.nbeta) be=-be
         sn=0
         total=0

         !--loop over all oscillators
         ! para-h2: j=0,2,....; ortho-h2: j=1,3,....
         ! ortho-d2: j=0,2,....; para-d2: j=1,3,....
         ipo=1
         if (law.eq.2.or.law.eq.5) ipo=2
         jt1=2*jterm
         if (ipo.eq.2) jt1=jt1+1
         do l=ipo,jt1,2
            j=l-1
            call bt(j,pj,x)

            !--sum over even values of j-prime
            snlg=0
            do lp=1,10,2
               jp=lp-1
               betap=(-j*(j+1)+jp*(jp+1))*x/2
               tmp=(2*jp+1)*pj*swe*4*sumh(j,jp,y)
               if (jj.eq.1.and.tmp.ge.small) then
                  write(nsyso,'(5x,f10.4,2i4,f10.6)') betap,j,jp,tmp
                  total=total+tmp
               endif
               bn=be+betap
               if (ifree.eq.1) then
                  ex=-(alp-abs(bn))**2/(4*alp)
                  if (bn.gt.zero) ex=ex-bn
                  add=exp(ex)/sqrt(4*pi*alp)
               else
                  add=sint(bn,bex,rdbex,sex,nbx,al,wt,tbart,&
                    betan,nbeta,maxbb)
               endif
               snlg=snlg+tmp*add
            enddo

            !--sum over the odd values of j-prime
            snlk=0
            do lp=2,10,2
               jp=lp-1
               betap=(-j*(j+1)+jp*(jp+1))*x/2
               tmp=(2*jp+1)*pj*swo*4*sumh(j,jp,y)
               if (jj.eq.1.and.tmp.ge.small) then
                  write(nsyso,'(5x,f10.4,2i4,f10.6)') betap,j,jp,tmp
                  total=total+tmp
               endif
               bn=be+betap
               if (ifree.eq.1) then
                  ex=-(alp-abs(bn))**2/(4*alp)
                  if (bn.gt.zero) ex=ex-bn
                  add=exp(ex)/sqrt(4*pi*alp)
               else
                  add=sint(bn,bex,rdbex,sex,nbx,al,wt,tbart,&
                    betan,nbeta,maxbb)
               endif
               snlk=snlk+tmp*add
            enddo

            !--continue the j loop
            sn=sn+snlg+snlk
         enddo
         if (jj.eq.1.and.iprt.eq.1) then
            write(nsyso,'(5x,''total'',5x,f10.6)') total
         endif

         !--continue the beta loop
         if (jj.le.nbeta) ssm(k,nal,itemp)=sn
         if (jj.ge.nbeta) ssp(k,nal,itemp)=sn
      enddo

      !--record the results
      if (iprt.eq.1) write(nsyso,&
        '(/4x,'' beta'',7x,''s(alpha,beta)'',7x,&
        &''s(alpha,-beta)'',7x,''ss(alpha,beta)'',&
        &5x,''ss(alpha,-beta)'')')
      do i=1,nbeta
         jprt=mod(i-1,nbint)+1
         if (i.eq.nbeta) jprt=1
         down=ssm(i,nal,itemp)*exb(i)
         up=0
         if (exb(i).ne.zero) up=ssp(i,nal,itemp)/exb(i)
         if (iprt.eq.1.and.jprt.eq.1) write(nsyso,&
           '(f10.4,1p,e18.5,3e20.5)') betan(i),&
           up,down,ssp(i,nal,itemp),ssm(i,nal,itemp)
      enddo

      !--check moments of calculated s(alpha,beta).
      sum0=0
      bel=0
      ff1l=0
      ff2l=0
      do nbe=1,nbeta
         be=betan(nbe)
         ff2=ssm(nbe,nal,itemp)
         ff1=ssp(nbe,nal,itemp)
         if (nbe.ne.1) then
            sum0=sum0+(be-bel)*(ff1l+ff2l+ff1+ff2)/2
            ff1l=ff1
            ff2l=ff2
            bel=be
         else
            bel=be
            ff1l=ff1
            ff2l=ff2
            sum0=0
         endif
      enddo
      write(nsyso,'(''     normalization check ='',f8.4)') sum0

   !--continue the alpha loop
   enddo
   deallocate(sex)
   deallocate(rdbex)
   deallocate(bex)
   deallocate(exb)
   deallocate(betan)
   return
   end subroutine coldh

   subroutine bt(j,pj,x)
   !--------------------------------------------------------------------
   ! Statistical weight factor
   ! for cold hydrogen or deuterium calculation
   !--------------------------------------------------------------------
   ! externals
   integer::j
   real(kr)::pj,x
   ! internals
   integer::i,k
   real(kr)::yy,a,b
   real(kr),parameter::half=0.5e0_kr

   yy=half*j*(j+1)
   a=(2*j+1)*exp(-yy*x)
   b=0
   do i=1,10
      k=2*i-2
      if (mod(j,2).eq.1) k=k+1
      yy=half*k*(k+1)
      b=b+(2*k+1)*exp(-yy*x)
   enddo
   pj=a/(2*b)
   return
   end subroutine bt

   real(kr) function sumh(j,jp,y)
   !--------------------------------------------------------------------
   ! Does sum over Bessel functions and Clebsch-Gordon coefficients
   ! for cold hydrogen or deuterium calculation.
   !--------------------------------------------------------------------
   ! externals
   integer::j,jp
   real(kr)::y
   ! internals
   integer::imk,ipk1,mpk,ipk,n,n1
   real(kr)::sum1,sum2

   if (j.eq.0) then
      sum2=(sjbes(jp,y)*cn(j,jp,jp))**2
   else if (jp.eq.0) then
      sum2=(sjbes(j,y)*cn(j,0,j))**2
   else
      sum1=0
      imk=iabs(j-jp)+1
      ipk1=j+jp+1
      mpk=ipk1-imk
      if (mpk.le.9) then
         ipk=ipk1
      else
         ipk=imk+9
      endif
      do n=imk,ipk
         n1=n-1
         sum1=sum1+(sjbes(n1,y)*cn(j,jp,n1))**2
      enddo
      sum2=sum1
   endif
   sumh=sum2
   return
   end function sumh

   real(kr) function cn(jj,ll,nn)
   !--------------------------------------------------------------------
   ! Clebsch-Gordon coefficients
   ! for cold hydrogen or deuterium calculation
   !--------------------------------------------------------------------
   ! externals
   integer::jj,ll,nn
   ! internals
   integer::i,kdet,kdel,ka1,ka2,ka3,ka4,kb1,kb2,kb3,kb4,iwign
   real(kr)::s,fact,zi,a1,a2,a3,a4,b1,b2,b3,b4,rat,wign
   real(kr),parameter::zero=0

   kdet=(jj+ll+nn)/2
   kdel=jj+ll+nn-2*kdet
   if (kdel.eq.0) then
      ka1=jj+ll+nn
      ka2=jj+ll-nn
      ka3=jj-ll+nn
      ka4=ll-jj+nn
      kb1=ka1/2
      kb2=ka2/2
      kb3=ka3/2
      kb4=ka4/2
      s=0
      fact=1
      do i=1,ka1
         zi=i
         s=s+log(zi)
      enddo
      if (s.gt.zero) fact=exp(s)
      a1=sqrt(fact)
      s=0
      fact=1
      do i=1,ka2
         zi=i
         s=s+log(zi)
      enddo
      if (s.gt.zero) fact=exp(s)
      a2=sqrt(fact)
      s=0
      fact=1
      do i=1,ka3
         zi=i
         s=s+log(zi)
      enddo
      if (s.gt.zero) fact=exp(s)
      a3=sqrt(fact)
      s=0
      fact=1
      do i=1,ka4
         zi=i
         s=s+log(zi)
      enddo
      if (s.gt.zero) fact=exp(s)
      a4=sqrt(fact)
      s=0
      b1=1
      do i=1,kb1
         zi=i
         s=s+log(zi)
      enddo
      if (s.gt.zero) b1=exp(s)
      s=0
      b2=1
      do i=1,kb2
         zi=i
         s=s+log(zi)
      enddo
      if (s.gt.zero) b2=exp(s)
      s=0
      b3=1
      do i=1,kb3
         zi=i
         s=s+log(zi)
      enddo
      if (s.gt.zero) b3=exp(s)
      s=0
      b4=1
      do i=1,kb4
         zi=i
         s=s+log(zi)
      enddo
      if (s.gt.zero) b4=exp(s)
      rat=2*nn+1
      rat=rat/(jj+ll+nn+1)
      iwign=(jj+ll-nn)/2
      wign=(-1)**iwign
      wign=wign*sqrt(rat)*b1/a1*a2/b2*a3/b3*a4/b4
   else
      wign=0
   endif
   cn=wign
   return
   end function cn

   real(kr) function sjbes(n,x)
   !--------------------------------------------------------------------
   ! Bessel functions for cold hydrogen or deuterium calculation
   !--------------------------------------------------------------------
   use util ! provides error
   ! externals
   integer::n
   real(kr)::x
   ! internals
   integer::i,k,l,iii,kmax,nm
   real(kr)::w,bessel,y,z,sj,t1,t2,t3
   character(60)::strng
   real(kr),parameter::huge=1.e25_kr
   real(kr),parameter::small=2.e-38_kr
   real(kr),parameter::break1=3.e4_kr
   real(kr),parameter::break2=7.e-4_kr
   real(kr),parameter::break3=0.2e0_kr
   real(kr),parameter::one=1
   real(kr),parameter::ten=10
   real(kr),parameter::hund=100
   real(kr),parameter::zero=0

   !--check for large arguments
   if (n.ge.30000.or.x.gt.break1) then
      write(strng,'(&
        &''value is not accurate  n = '',i7,10x,''x = '',e14.7)') n,x
      call mess('sjbes',strng,' ')
      sjbes=0
      return
   endif

   !--check for bad arguments
   if (x.lt.zero.or.n.lt.0) then
      write(strng,'(&
        &''argument is invalid  n = '',i7,10x,''x = '',e14.7)') n,x
      call error('sjbes',strng,' ')
      sjbes=0
      return
   endif

   !--compute normal values
   if (x.le.break2) then
      w=1
      if (n.eq.0) then
         bessel=w
      else if (n.gt.10) then
         bessel=0
      else
         t1=3
         t2=1
         do i=1,n
            t3=t2*x/t1
            t1=t1+2
            t2=t3
         enddo
         bessel=t3
      endif
   else
      if (x.lt.break3) then
         y=x**2
         w=1-y*(1-y/20)/6
      else
         w=sin(x)/x
      endif
      if (n.eq.0) then
         bessel=w
      else
         if (x.ge.hund) then
            l=int(x/50+18)
         else if (x.ge.ten) then
            l=int(x/10+10)
         else if (x.gt.one) then
            l=int(x/2+5)
         else
            l=5
         endif
         iii=int(x)
         kmax=n
         if (iii.gt.n) kmax=iii
         nm=kmax+l
         z=1/x
         t3=0
         t2=small
         do i=1,nm
            k=nm-i
            t1=(2*k+3)*z*t2-t3
            if (n.eq.k) sj=t1
            if (abs(t1).ge.huge) then
               t1=t1/huge
               t2=t2/huge
               sj=sj/huge
            endif
            t3=t2
            t2=t1
         enddo
         bessel=w*sj/t1
      endif
   endif
   sjbes=bessel
   return
   end function sjbes

   real(kr) function terpk(ska,nka,delta,be)
   !--------------------------------------------------------------------
   ! Interpolate in a table of ska(kappa) for a required kappa.
   ! Called by coldh.
   !--------------------------------------------------------------------
   ! externals
   integer::nka
   real(kr)::ska(nka),delta,be
   ! internals
   integer::i
   real(kr)::bt,btp

   terpk=1
   if (be.gt.nka*delta) return
   i=int(be/delta)
   if (i.lt.nka-1) then
      bt=i*delta
      btp=bt+delta
      i=i+1
      terpk=ska(i)+(be-bt)*(ska(i+1)-ska(i))/(btp-bt)
   endif
   return
   end function terpk

   subroutine copys(sab,nbeta,nalpha,ntempr)
   !--------------------------------------------------------------------
   ! Copy sab for principal scatterer to scratch tape on unit 10.
   !--------------------------------------------------------------------
   use util ! provides openz
   ! externals
   integer::nbeta,nalpha,ntempr
   real(kr)::sab(nbeta,nalpha,ntempr)
   ! internals
   integer::i,j,k,nscr

   nscr=-10
   call openz(nscr,1)
   do k=1,ntempr
      do j=1,nalpha
         write(-nscr) (sab(i,j,k),i=1,nbeta)
      enddo
   enddo
   return
   end subroutine copys

   subroutine coher(lat,natom,b,nbe,maxb,emax)
   !--------------------------------------------------------------------
   ! Compute Bragg energies and associated structure factors
   ! for coherent elastic scattering from graphite, Be, or BeO.
   !--------------------------------------------------------------------
   use mainio  ! provides nsyso
   use physics ! provides pi,hbar,ev,amu,amassn
   use util    ! provides timer,error
   ! externals
   integer::lat,natom,nbe,maxb
   real(kr)::b(maxb),emax
   ! internals
   integer::i,j,k,imax,jmin,idone,ifl,i1m,nw
   integer::i1,i2,i3,l1,l2,l3,i2m,i3m
   real(kr)::time,twopis,amne,econ,tsqx
   real(kr)::a,c,amsc,scoh,c1,c2
   real(kr)::recon,scon,wint,t2,ulim,phi
   real(kr)::w1,w2,w3,tsq,tau,w,f
   real(kr)::x,st,sf,bel,be,bs
   real(kr),parameter::gr1=2.4573e-8_kr
   real(kr),parameter::gr2=6.700e-8_kr
   real(kr),parameter::gr3=12.011e0_kr
   real(kr),parameter::gr4=5.50e0_kr
   real(kr),parameter::be1=2.2856e-8_kr
   real(kr),parameter::be2=3.5832e-8_kr
   real(kr),parameter::be3=9.01e0_kr
   real(kr),parameter::be4=7.53e0_kr
   real(kr),parameter::beo1=2.695e-8_kr
   real(kr),parameter::beo2=4.39e-8_kr
   real(kr),parameter::beo3=12.5e0_kr
   real(kr),parameter::beo4=1.0e0_kr
   real(kr),parameter::al1=4.04e-8_kr
   real(kr),parameter::al3=26.7495e0_kr
   real(kr),parameter::al4=1.495e0_kr
   real(kr),parameter::pb1=4.94e-8_kr
   real(kr),parameter::pb3=207.e0_kr
   real(kr),parameter::pb4=1.e0_kr
   real(kr),parameter::fe1=2.86e-8_kr
   real(kr),parameter::fe3=55.454e0_kr
   real(kr),parameter::fe4=12.9e0_kr
   real(kr),parameter::twothd=0.666666666667e0_kr
   real(kr),parameter::sqrt3=1.732050808e0_kr
   real(kr),parameter::toler=1.e-6_kr
   real(kr),parameter::eps=.05e0_kr
   real(kr),parameter::zero=0

   !--write header
   call timer(time)
   write(nsyso,'(/'' bragg edges for coherent elastic scattering'',&
     &25x,f8.1,''s'')') time

   !--initialize.
   twopis=(2*pi)**2
   amne=amassn*amu
   econ=ev*8*(amne/hbar)/hbar
   recon=1/econ
   tsqx=econ/20
   if (lat.eq.1) then
      ! graphite constants.
      a=gr1
      c=gr2
      amsc=gr3
      scoh=gr4/natom
   else if (lat.eq.2) then
      !  beryllium constants
      a=be1
      c=be2
      amsc=be3
      scoh=be4/natom
   else if (lat.eq.3) then
      !  beryllium oxide constants
      a=beo1
      c=beo2
      amsc=beo3
      scoh=beo4/natom
   else if (lat.eq.4) then
      ! aluminum constants
      a=al1
      amsc=al3
      scoh=al4/natom
   else if (lat.eq.5) then
      ! lead constants
      a=pb1
      amsc=pb3
      scoh=pb4/natom
   else if (lat.eq.6) then
      ! iron constants
      a=fe1
      amsc=fe3
      scoh=fe4/natom
   else
      call error('coh','illegal lat.',' ')
   endif
   if (lat.lt.4) then
      c1=4/(3*a*a)
      c2=1/(c*c)
      scon=scoh*(4*pi)**2/(2*a*a*c*sqrt3*econ)
   else if (lat.ge.4.and.lat.le.5) then
      c1=3/(a*a)
      scon=scoh*(4*pi)**2/(16*a*a*a*econ)
   else if (lat.eq.6) then
      c1=2/(a*a)
      scon=scoh*(4*pi)**2/(8*a*a*a*econ)
   endif
   wint=0
   t2=hbar/(2*amu*amsc)
   ulim=econ*emax
   ifl=1
   nw=maxb

   !--compute lattice factors for hexagonal lattices
   if (lat.gt.3) go to 210
   phi=ulim/twopis
   i1m=int(a*sqrt(phi))
   i1m=i1m+1
   k=0
   do i1=1,i1m
      l1=i1-1
      i2m=int((l1+sqrt(3*(a*a*phi-l1*l1)))/2)
      i2m=i2m+1
      do i2=i1,i2m
         l2=i2-1
         x=phi-c1*(l1*l1+l2*l2-l1*l2)
         i3m=0
         if (x.gt.zero) i3m=int(c*sqrt(x))
         i3m=i3m+1
         do i3=1,i3m
            l3=i3-1
            w1=2
            if (l1.eq.l2) w1=1
            w2=2
            if (l1.eq.0.or.l2.eq.0) w2=1
            if (l1.eq.0.and.l2.eq.0) w2=1
            if (l1.eq.0.and.l2.eq.0) w2=w2/2
            w3=2
            if (l3.eq.0) w3=1
            tsq=tausq(l1,l2,l3,c1,c2,twopis)
            if (tsq.gt.zero.and.tsq.le.ulim) then
               tau=sqrt(tsq)
               w=exp(-tsq*t2*wint)*w1*w2*w3/tau
               f=w*formf(lat,l1,l2,l3)
               if (k.le.0.or.tsq.le.tsqx) then
                  k=k+1
                  if ((2*k).gt.nw) call error('coh',&
                    'storage exceeded',' ')
                  b(ifl+2*k-2)=tsq
                  b(ifl+2*k-1)=f
               else
                  i=0
                  idone=0
                  do while (i.lt.k.and.idone.eq.0)
                     i=i+1
                     if (tsq.ge.b(ifl+2*i-2).and.&
                       tsq.lt.(1+eps)*b(ifl+2*i-2)) then
                           b(ifl+2*i-1)=b(ifl+2*i-1)+f
                        idone=1
                     endif
                  enddo
                  if (idone.eq.0) then
                     k=k+1
                     if ((2*k).gt.nw) call error('coh',&
                       'storage exceeded',' ')
                     b(ifl+2*k-2)=tsq
                     b(ifl+2*k-1)=f
                  endif
               endif
            endif
            tsq=tausq(l1,-l2,l3,c1,c2,twopis)
            if (tsq.gt.zero.and.tsq.le.ulim) then
               tau=sqrt(tsq)
               w=exp(-tsq*t2*wint)*w1*w2*w3/tau
               f=w*formf(lat,l1,-l2,l3)
               if (k.le.0.or.tsq.le.tsqx) then
                  k=k+1
                  if ((2*k).gt.nw) call error('coh',&
                    'storage exceeded',' ')
                  b(ifl+2*k-2)=tsq
                  b(ifl+2*k-1)=f
               else
                  i=0
                  idone=0
                  do while (i.lt.k.and.idone.eq.0)
                     i=i+1
                     if (tsq.ge.b(ifl+2*i-2).and.&
                       tsq.lt.(1+eps)*b(ifl+2*i-2)) then
                        b(ifl+2*i-1)=b(ifl+2*i-1)+f
                        idone=1
                     endif
                  enddo
                  if (idone.eq.0) then
                     k=k+1
                     if ((2*k).gt.nw) call error('coh',&
                       'storage exceeded',' ')
                     b(ifl+2*k-2)=tsq
                     b(ifl+2*k-1)=f
                  endif
               endif
            endif
         enddo
      enddo
   enddo
   imax=k-1
   go to 220

   !--compute lattice factors for fcc lattices
  210 continue
   if (lat.gt.5) go to 215
   phi=ulim/twopis
   i1m=int(a*sqrt(phi))
   i1m=15
   k=0
   do i1=-i1m,i1m
      i2m=i1m
      do i2=-i2m,i2m
         i3m=i1m
         do i3=-i3m,i3m
            tsq=taufcc(i1,i2,i3,c1,twothd,twopis)
            if (tsq.gt.zero.and.tsq.le.ulim) then
               tau=sqrt(tsq)
               w=exp(-tsq*t2*wint)/tau
               f=w*formf(lat,i1,i2,i3)
               k=k+1
               if ((2*k).gt.nw) call error('coh','storage exceeded',' ')
               b(ifl+2*k-2)=tsq
               b(ifl+2*k-1)=f
            endif
         enddo
      enddo
   enddo
   imax=k-1
   go to 220

   !--compute lattice factors for bcc lattices
  215 continue
   phi=ulim/twopis
   i1m=int(a*sqrt(phi))
   i1m=15
   k=0
   do i1=-i1m,i1m
      i2m=i1m
      do i2=-i2m,i2m
         i3m=i1m
         do i3=-i3m,i3m
            tsq=taubcc(i1,i2,i3,twopis)
            if (tsq.gt.zero.and.tsq.le.ulim) then
               tau=sqrt(tsq)
               w=exp(-tsq*t2*wint)/tau
               f=w*formf(lat,i1,i2,i3)
               k=k+1
               if ((2*k).gt.nw) call error('coh','storage exceeded',' ')
               b(ifl+2*k-2)=tsq
               b(ifl+2*k-1)=f
            endif
         enddo
      enddo
   enddo
   imax=k-1

   !--sort lattice factors

  220 continue
   do i=1,imax
      jmin=i+1
      do j=jmin,k
         if (b(ifl+2*j-2).lt.b(ifl+2*i-2)) then
            st=b(ifl+2*i-2)
            sf=b(ifl+2*i-1)
            b(ifl+2*i-2)=b(ifl+2*j-2)
            b(ifl+2*i-1)=b(ifl+2*j-1)
            b(ifl+2*j-2)=st
            b(ifl+2*j-1)=sf
         endif
      enddo
   enddo
   k=k+1
   b(ifl+2*k-2)=ulim
   b(ifl+2*k-1)=b(ifl+2*k-3)
   nw=2*k

   !--convert to practical units
   !--and combine duplicate bragg edges.
   bel=-1
   j=0
   do i=1,k
      be=b(ifl+2*i-2)*recon
      bs=b(ifl+2*i-1)*scon
      if (be-bel.lt.toler) then
         b(ifl+2*j-1)=b(ifl+2*j-1)+bs
      else
         j=j+1
         b(ifl+2*j-2)=be
         b(ifl+2*j-1)=bs
         bel=be
      endif
   enddo
   nbe=j
   maxb=2*nbe
   write(nsyso,'(/''   found'',i5,'' edges below'',&
     &f6.2,'' ev'')') nbe,emax
   return

   contains

      real(kr) function tausq(m1,m2,m3,c1,c2,twopis)
      integer::m1,m2,m3
      real(kr)::c1,c2,twopis
      tausq=(c1*(m1*m1+m2*m2+m1*m2)+(m3*m3*c2))*twopis
      return
      end function tausq

      real(kr) function taufcc(m1,m2,m3,c1,twothd,twopis)
      integer::m1,m2,m3
      real(kr)::c1,twothd,twopis
      taufcc=c1*(m1*m1+m2*m2+m3*m3+twothd*m1*m2&
        +twothd*m1*m3-twothd*m2*m3)*twopis
      return
      end function taufcc

      real(kr) function taubcc(m1,m2,m3,twopis)
      integer::m1,m2,m3
      real(kr)::twopis
      taubcc=c1*(m1*m1+m2*m2+m3*m3+m1*m2+m2*m3+m1*m3)*twopis
      return
      end function taubcc

   end subroutine coher

   subroutine skold(itemp,temp,ssm,nalpha,nbeta,ntempr)
   !--------------------------------------------------------------------
   ! use skold approximation to add in the effects
   ! of intermolecular coherence.
   !--------------------------------------------------------------------
   use mainio  ! provides nsyso
   use physics ! provides bk,ev,hbar,amassn,amu
   use endf    ! provides terp1
   ! externals
   integer::itemp,nalpha,nbeta,ntempr
   real(kr)::temp
   real(kr)::ssm(nbeta,nalpha,ntempr)
   ! internals
   integer::i,j,k,kk,nal,ibeta,iprt,jprt
   real(kr)::tev,sc,amass,al,sk,ap,be,ss,s1,s2
   real(kr)::sum0,sum1,ff1l,ff2l,bel,ff1,ff2,waven
   real(kr)::scoh(1000)
   real(kr),parameter::angst=1.e-8_kr
   real(kr),parameter::therm=.0253e0_kr
   real(kr),parameter::zero=0.0

   !--apply the skold approximation
   tev=bk*abs(temp)
   sc=1
   if (lat.eq.1) sc=therm/tev
   amass=awr*amassn*amu
   do i=1,nbeta
      do j=1,nalpha
         al=alpha(j)*sc/arat
         waven=angst*sqrt(2*amass*tev*ev*al)/hbar
         sk=terpk(ska,nka,dka,waven)
         ap=alpha(j)/sk
         do k=1,nalpha
            kk=k
            if (ap.lt.alpha(k)) exit
         enddo
         if (kk.eq.1) kk=2
         if (ssm(i,kk-1,itemp).eq.zero.or.ssm(i,kk,itemp).eq.zero) then
            scoh(j)=zero
         else
         call terp1(alpha(kk-1),ssm(i,kk-1,itemp),&
           alpha(kk),ssm(i,kk,itemp),ap,scoh(j),5)
         endif
         scoh(j)=scoh(j)*sk
      enddo
      do j=1,nalpha
         ssm(i,j,itemp)=(1-cfrac)*ssm(i,j,itemp)+cfrac*scoh(j)
      enddo
   enddo

   !--report the results
   if (iprint.eq.2) write(nsyso,&
     '(/'' results after applying skold approximation'')')
   do nal=1,nalpha
      iprt=mod(nal-1,naint)+1
      if (nal.eq.nalpha) iprt=1
      al=alpha(nal)*sc/arat
      if (iprt.eq.1.and.iprint.eq.2) write(nsyso,&
        '(/3x,''alpha='',f10.5)') al
      if (iprt.eq.1.and.iprint.eq.2) write(nsyso,&
        '(/4x,'' beta'',7x,''s(alpha,beta)'',7x,''ss(alpha,beta)'',&
        &5x,''ss(alpha,-beta)'')')
      do i=1,nbeta
         be=beta(i)*sc
         ss=ssm(i,nal,itemp)
         s1=ss*exp(-be/2)
         s2=ss*exp(-be)
         jprt=mod(i-1,nbint)+1
         if (i.eq.nbeta) jprt=1
         if (iprt.eq.1.and.jprt.eq.1.and.iprint.eq.2)&
           write(nsyso,'(f10.4,1pe18.5,1p,2e20.5)') beta(i),s1,s2,ss
      enddo
      if (iprt.eq.1) then
         sum0=0
         sum1=0
         ff1l=0
         ff2l=0
         bel=0
         do ibeta=1,nbeta
            be=beta(ibeta)
            ff2=ssm(ibeta,nal,itemp)
            ff1=ssm(ibeta,nal,itemp)*exp(-be)
            if (ibeta.gt.1) then
               sum0=sum0+(be-bel)*(ff1l+ff2l+ff1+ff2)/2
               sum1=sum1+(be-bel)*(ff2l*bel+ff2*be-ff1l*bel-ff1*be)/2
               ff1l=ff1
               ff2l=ff2
               bel=be
            else
               bel=be
               ff1l=ff1
               ff2l=ff2
               sum0=0
               sum1=0
            endif
         enddo
         sum1=sum1/al
         if (iprint.eq.2) then
            write(nsyso,'(''     normalization check ='',f8.4)') sum0
            write(nsyso,'(''          sum rule check ='',f8.4)') sum1
         else if (iprint.eq.1) then
            write(nsyso,'(1x,f10.4,2f10.4)') al,sum0,sum1
         endif
      endif
   enddo
   return
   end subroutine skold

   real(kr) function formf(lat,l1,l2,l3)
   !--------------------------------------------------------------------
   ! Compute form factors for the specified lattice.
   !       lat=1    graphite
   !       lat=2    Be
   !       lat=3    BeO
   !       lat=4,5  fcc lattice (aluminum, lead)
   !       lat=6    bcc lattice (iron)
   !--------------------------------------------------------------------
   use physics ! provides pi
   ! externals
   integer::lat,l1,l2,l3
   ! internals
   integer::i
   real(kr)::e1,e2,e3
   real(kr),parameter::c1=7.54e0_kr
   real(kr),parameter::c2=4.24e0_kr
   real(kr),parameter::c3=11.31e0_kr

   if (lat.eq.1) then
      ! graphite.
      i=l3/2
      if ((2*i).ne.l3) then
      formf=sin(pi*(l1-l2)/3)**2
      else
         formf=(6+10*cos(2*pi*(l1-l2)/3))/4
      endif
   else if (lat.eq.2) then
      ! beryllium.
      formf=1+cos(2*pi*(2*l1+4*l2+3*l3)/6)
   else if (lat.eq.3) then
      ! beryllium oxide.
      formf=(1+cos(2*pi*(2*l1+4*l2+3*l3)/6))&
        *(c1+c2+c3*cos(3*pi*l3/4))
   else if (lat.eq.4.or.lat.eq.5) then
      ! fcc lattices.
      e1=2*pi*l1
      e2=2*pi*(l1+l2)
      e3=2*pi*(l1+l3)
      formf=(1+cos(e1)+cos(e2)+cos(e3))**2+(sin(e1)+sin(e2)+sin(e3))**2
   else if (lat.eq.6) then
      ! bcc lattices.
      e1=2*pi*(l1+l2+l3)
      formf=(1+cos(e1))**2+(sin(e1))**2
   endif
   return
   end function formf

   subroutine endout(ntempr,bragg,nedge,maxb,scr,mscr,isym,ilog)
   !--------------------------------------------------------------------
   ! ENDF output routine
   !--------------------------------------------------------------------
   use mainio  ! provide nsyso
   use endf    ! provides endf routines, npage, iverf
   use physics ! provide bk (boltzmann constant)
   use util    ! provides timer,repoz,error,sigfig
   ! externals
   integer::ntempr,maxb,mscr,nedge,isym,ilog
   real(kr)::bragg(maxb)
   real(kr)::scr(mscr)
   ! internals
   real(kr)::time,sum,suml,w,e,sb,sbs,srat,sc,be
   integer::i,j,k,l,nt,jmax,ii,jj,j1,ndw,nbt,ntf,ll
   integer::ncards,idone,nc,n,nb,nw,nprnt,nscr
   character(66)::text
   character(4)::t(17)
   real(kr)::z(17)
   equivalence(t(1),z(1))
   real(kr),parameter::small=1.e-9_kr
   real(kr),parameter::tiny=-999.e0_kr
   real(kr),parameter::tol=0.9e-7_kr
   real(kr),parameter::up=1.01e0_kr
   real(kr),parameter::therm=.0253e0_kr
   real(kr),parameter::zero=0

   !--write header
   call timer(time)
   write(nsyso,'(//'' endf output'',57x,f8.1,''s'')') time
   if (isym.eq.0) then
      write(nsyso,'(/'' endf output is s'')')
   else if (isym.eq.1) then
      write(nsyso,'(/'' endf output is s for +&- beta'')')
   else if (isym.eq.2) then
      write(nsyso,'(/'' endf output is ss for - beta'')')
   else if (isym.eq.3) then
      write(nsyso,'(/'' endf output is ss for +&- beta'')')
   endif
   if (ilog.eq.1) write(nsyso,'('' log10 s or ss is written'')')

   !---compute bound scattering cross sections
   sb=spr*((1+awr)/awr)**2
   if (aws.ne.zero) sbs=sps*((1+aws)/aws)**2

   !--for mixed moderators, merge ssm results
   if (nss.ne.0.and.b7.le.zero) then
      srat=sbs/sb
      nscr=-10
      call repoz(nscr)
      do k=1,ntempr
         do j=1,nalpha
            read(-nscr) (scr(i),i=1,nbeta)
            do i=1,nbeta
               ssm(i,j,k)=srat*ssm(i,j,k)+scr(i)
            enddo
         enddo
      enddo
   endif

   !---display endf t-effective and debye-waller integral
   do i=1,ntempr
      if (nss.eq.0.or.b7.gt.zero) then
         dwpix(i)=dwpix(i)/(awr*tempr(i)*bk)
      else
         dwpix(i)=dwpix(i)/(aws*tempr(i)*bk)
         dwp1(i)=dwp1(i)/(awr*tempr(i)*bk)
      endif
   enddo
   if (nss.eq.0.or.b7.gt.zero) then
      write(nsyso,'(/''    i     tempr       tempf       awr*w''/&
        &(i5,1p,3e12.3))') (i,tempr(i),tempf(i),awr*dwpix(i),i=1,ntempr)
   else
      write(nsyso,'(/26x,''principal'',15x,''secondary''/&
        &''    i     tempr       tempf       awr*w'',&
        &''       tempf       aws*w''/(i5,1p,5e12.3))')&
        (i,tempr(i),tempf1(i),awr*dwp1(i),&
        tempf(i),aws*dwpix(i),i=1,ntempr)
   endif

   !--write endf-6 file 1.
   if (iel.eq.0.and.twt.eq.zero) iel=-1
   write(nsyso,'(//)')
   nprnt=0
   nsp=0
   nsc=0
   math=1
   mfh=0
   mth=0
   nsh=0
   text=' '
   read(text,'(16a4,a2)') (t(i),i=1,17)
   call tpidio(0,nout,nprnt,z,nb,nw)
   math=mat
   mfh=1
   mth=451
   scr(1)=za
   scr(2)=awr
   scr(3)=-1
   scr(4)=0
   scr(5)=0
   scr(6)=0
   call contio(0,nout,nprnt,scr(1),nb,nw)
   scr(1)=0
   scr(2)=0
   scr(3)=0
   scr(4)=0
   scr(5)=0
   scr(6)=6
   call contio(0,nout,nprnt,scr(1),nb,nw)
   scr(1)=1
   scr(2)=0
   scr(3)=0
   scr(4)=0
   scr(5)=12
   scr(6)=6
   call contio(0,nout,nprnt,scr(1),nb,nw)
   scr(1)=0
   scr(2)=0
   scr(3)=0
   scr(4)=0
   scr(5)=0
   scr(6)=2
   if (iel.ne.0) scr(6)=3
   n=6
   nc=0
   idone=0
   do while (idone.eq.0)
      text='$'
      read(nsysi,*) text
      if (text.eq.'$') then
         idone=1
      else
         nc=nc+1
         read(text,'(16a4,a2)') (t(j),j=1,17)
         do j=1,17
            scr(n+j)=z(j)
         enddo
         n=n+17
         if (n.gt.mscr) call error('endout',&
           &'scratch storage exceeded for hollerith data',' ')
      endif
   enddo
   scr(5)=17*nc
   l=1
   call hdatio(0,nout,nprnt,scr(l),nb,nw)
   do while (nb.ne.0)
      l=l+nw
      call moreio(0,nout,nprnt,scr(l),nb,nw)
   enddo
   scr(1)=0
   scr(2)=0
   scr(3)=1
   scr(4)=451
   scr(5)=5+nc+1
   if (iel.ne.0) scr(5)=scr(5)+1
   scr(6)=0
   ii=6
   if (iel.ne.0) then
      scr(ii+1)=0
      scr(ii+2)=0
      scr(ii+3)=7
      scr(ii+4)=2
      if (iel.lt.0) ncards=3+(2*ntempr+4)/6
      if (iel.gt.0) ncards=3+(2*nedge+4)/6
      if (iel.gt.0.and.ntempr.gt.1)&
        ncards=ncards+(ntempr-1)*(1+(nedge+5)/6)
      scr(ii+5)=ncards
      scr(ii+6)=0
      ii=ii+6
   endif
   scr(ii+1)=0
   scr(ii+2)=0
   scr(ii+3)=7
   scr(ii+4)=4
   ncards=2+(2*nalpha+4)/6
   if (ntempr.gt.1) ncards=ncards+(ntempr-1)*(1+(nalpha+5)/6)
   ncards=5+nbeta*ncards
   scr(ii+5)=ncards
   scr(ii+6)=0
   nw=2
   if (iel.ne.0) nw=3
   call dictio(0,nout,nprnt,scr(1),nb,nw)
   call asend(nout,nprnt)
   call afend(nout,nprnt)

   !--write incoherent elastic part
   if (iel.lt.0) then
      mfh=7
      mth=2
      scr(1)=za
      scr(2)=awr
      scr(3)=2
      scr(4)=0
      scr(5)=0
      scr(6)=0
      call contio(0,nout,nprnt,scr(1),nb,nw)
      scr(1)=sb*npr
      scr(2)=0
      scr(3)=0
      scr(4)=0
      scr(5)=1
      ndw=ntempr
      if (ndw.eq.1) ndw=2
      scr(6)=ndw
      scr(7)=ndw
      scr(8)=2
      do i=1,ndw
         if (i.le.ntempr) then
            scr(2*i+7)=tempr(i)
            scr(2*i+8)=sigfig(dwpix(i),7,0)
         else
            scr(2*i+7)=scr(2*i+5)
            scr(2*i+8)=scr(2*i+6)
         endif
      enddo
      call tab1io(0,nout,nprnt,scr(1),nb,nw)
      call asend(nout,nprnt)
   endif

   !--write coherent elastic part
   if (iel.ge.1) then
      mfh=7
      mth=2
      scr(1)=za
      scr(2)=awr
      scr(3)=1
      scr(4)=0
      scr(5)=0
      scr(6)=0
      call contio(0,nout,nprnt,scr(1),nb,nw)

      !--thin out the 1/e part at high energies.
      sum=0
      suml=0
      w=dwpix(1)
      if (nss.gt.0.and.b7.eq.zero) w=(dwpix(1)+dwp1(1))/2
      do j=1,nedge
         e=bragg(1+2*j-2)
         sum=sum+exp(-4*w*e)*bragg(1+2*j-1)
         if (sum-suml.gt.tol*sum) then
            jmax=j
            suml=sum
         endif
      enddo

      !--now output the records for each temperature
      do i=1,ntempr
         if (i.eq.1) then
            scr(1)=tempr(i)
            scr(2)=0
            scr(3)=ntempr-1
            scr(4)=0
            scr(5)=1
            scr(6)=jmax
            scr(7)=jmax
            scr(8)=1
            ii=8
            jj=0
            j1=0
            sum=0
            w=dwpix(i)
            if (nss.gt.0.and.b7.eq.zero) w=(dwpix(i)+dwp1(i))/2
            j=0
            do while (j.lt.nedge)
               j=j+1
               e=bragg(1+2*j-2)
               if (j.le.jmax) jj=jj+2
               scr(ii+jj-1-j1)=sigfig(e,7,0)
               scr(ii+jj-j1)=sum+exp(-4*w*e)*bragg(1+2*j-1)
               sum=scr(ii+jj-j1)
               scr(ii+jj-j1)=sigfig(scr(ii+jj-j1),7,0)
               if (j.ge.nedge.or.jj-j1.ge.npage) then
                  if (ii.ge.8) then
                     call tab1io(0,nout,nprnt,scr(1),nb,nw)
                     ii=0
                     j1=npage
                     idone=1
                  else
                     call moreio(0,nout,nprnt,scr(1),nb,nw)
                     j1=j1+npage
                  endif
               endif
            enddo
         else
            scr(1)=tempr(i)
            scr(2)=0
            scr(3)=2
            scr(4)=0
            scr(5)=jmax
            scr(6)=0
            ii=6
            jj=0
            j1=0
            sum=0
            w=dwpix(i)
            if (nss.gt.0.and.b7.eq.zero) w=(dwpix(i)+dwp1(i))/2
            do j=1,nedge
               if (j.le.jmax) jj=jj+1
               e=sigfig(bragg(1+2*jj-2),7,0)
               scr(ii+jj-j1)=sum+exp(-4*w*e)*bragg(1+2*jj-1)
               sum=scr(ii+jj-j1)
               scr(ii+jj-j1)=sigfig(scr(ii+jj-j1),7,0)
               if (j.ge.nedge.or.jj-j1.ge.npage) then
                  if (ii.ge.6) then
                     call listio(0,nout,nprnt,scr(1),nb,nw)
                     ii=0
                     j1=npage
                  else
                     call moreio(0,nout,nprnt,scr(1),nb,nw)
                     j1=j1+npage
                  endif
               endif
            enddo
         endif
      enddo
      call asend(nout,nprnt)
   endif

   !--write inelastic part
   mfh=7
   mth=4
   scr(1)=za
   scr(2)=awr
   scr(3)=0
   scr(4)=lat
   scr(5)=isym
   scr(6)=0
   call contio(0,nout,nprnt,scr(1),nb,nw)
   scr(1)=0
   scr(2)=0
   scr(3)=0
   if (ilog.ne.0) scr(3)=1
   scr(4)=0
   scr(5)=6
   if (nss.gt.0) scr(5)=6*(nss+1)
   scr(6)=nss
   scr(7)=npr*spr
   scr(8)=beta(nbeta)
   scr(9)=awr
   scr(10)=sigfig(therm*beta(nbeta),7,0)
   scr(11)=0
   scr(12)=npr
   if (nss.ne.0) then
      scr(13)=b7
      scr(14)=mss*sps
      scr(15)=aws
      scr(16)=0
      scr(17)=0
      scr(18)=mss
   endif
   call listio(0,nout,nprnt,scr(1),nb,nw)
   scr(1)=0
   scr(2)=0
   scr(3)=0
   scr(4)=0
   scr(5)=1
   nbt=nbeta
   if (isym.eq.1.or.isym.eq.3) nbt=2*nbeta-1
   scr(6)=nbt
   scr(7)=nbt
   scr(8)=4
   call tab2io(0,nout,nprnt,scr(1),nb,nw)
   ii=0
   do i=1,nbt
      do nt=1,ntempr
         sc=1
         if (lat.eq.1) sc=therm/(bk*tempr(nt))
         if (nt.eq.1) then
            scr(1)=tempr(nt)
            if (mod(isym,2).eq.0) scr(2)=beta(i)
            if (mod(isym,2).eq.1.and.i.lt.nbeta)&
              scr(2)=-beta(nbeta-i+1)
            if (mod(isym,2).eq.1.and.i.ge.nbeta)&
              scr(2)=beta(i-nbeta+1)
            be=scr(2)*sc
            scr(3)=ntempr-1
            scr(4)=0
            scr(5)=1
            scr(6)=nalpha
            scr(7)=nalpha
            scr(8)=4
            do j=1,nalpha
               scr(7+2*j)=alpha(j)
               if (isym.eq.0) then
                  if (ilog.eq.0) then
                     scr(8+2*j)=ssm(i,j,nt)*exp(-be/2)
                     if (scr(8+2*j).ge.small) then
                        scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                     else
                        scr(8+2*j)=sigfig(scr(8+2*j),6,0)
                     endif
                  else
                     scr(8+2*j)=tiny
                     if (ssm(i,j,nt).gt.zero) then
                        scr(8+2*j)=log(ssm(i,j,nt))-be/2
                        scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                     endif
                  endif
               else if (isym.eq.1) then
                  if (i.lt.nbeta) then
                     if (ilog.eq.0) then
                        scr(8+2*j)=ssm(nbeta-i+1,j,nt)*exp(be/2)
                        if (scr(8+2*j).ge.small) then
                           scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                        else
                           scr(8+2*j)=sigfig(scr(8+2*j),6,0)
                        endif
                     else
                        scr(8+2*j)=tiny
                        if (ssm(nbeta-i+1,j,nt).gt.zero) then
                           scr(8+2*j)=log(ssm(nbeta-i+1,j,nt))+be/2
                           scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                        endif
                     endif
                  else
                     if (ilog.eq.0) then
                        scr(8+2*j)=ssp(i-nbeta+1,j,nt)*exp(be/2)
                        if (scr(8+2*j).ge.small) then
                           scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                        else
                           scr(8+2*j)=sigfig(scr(8+2*j),6,0)
                        endif
                     else
                        scr(8+2*j)=tiny
                        if (ssp(i-nbeta+1,j,nt).gt.zero) then
                           scr(8+2*j)=log(ssp(i-nbeta+1,j,nt))+be/2
                           scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                        endif
                     endif
                  endif
               else if (isym.eq.2) then
                  if (ilog.eq.0) then
                     scr(8+2*j)=ssm(i,j,nt)
                    if (scr(8+2*j).ge.small) then
                        scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                     else
                        scr(8+2*j)=sigfig(scr(8+2*j),6,0)
                     endif
                  else
                     scr(8+2*j)=tiny
                     if (ssm(i,j,nt).gt.zero) then
                        scr(8+2*j)=log(ssm(i,j,nt))
                        scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                     endif
                  endif
               else if (isym.eq.3) then
                  if (i.lt.nbeta) then
                     if (ilog.eq.0) then
                        scr(8+2*j)=ssm(nbeta-i+1,j,nt)
                        if (scr(8+2*j).ge.small) then
                           scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                        else
                           scr(8+2*j)=sigfig(scr(8+2*j),6,0)
                        endif
                     else
                        scr(8+2*j)=tiny
                        if (ssm(nbeta-i+1,j,nt).gt.zero) then
                           scr(8+2*j)=log(ssm(nbeta-i+1,j,nt))
                           scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                        endif
                     endif
                  else
                     if (ilog.eq.0) then
                        scr(8+2*j)=ssp(i-nbeta+1,j,nt)
                        if (scr(8+2*j).ge.small) then
                           scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                        else
                           scr(8+2*j)=sigfig(scr(8+2*j),6,0)
                        endif
                     else
                        scr(8+2*j)=tiny
                        if (ssp(-nbeta+1,j,nt).gt.zero) then
                           scr(8+2*j)=log(ssp(i-nbeta+1,j,nt))
                           scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                        endif
                     endif
                  endif
               endif
               if (ilog.eq.0.and.scr(8+2*j).lt.smin) scr(8+2*j)=0
            enddo
            call tab1io(0,nout,nprnt,scr(1),nb,nw)
            ll=1+nw
            do while (nb.ne.0)
               call moreio(0,nout,nprnt,scr(ll),nb,nw)
               ll=ll+nw
            enddo
         else
            scr(1)=tempr(nt)
            if (mod(isym,2).eq.0) scr(2)=beta(i)
            if (mod(isym,2).eq.1.and.i.lt.nbeta)&
              scr(2)=-beta(nbeta-i+1)
            if (mod(isym,2).eq.1.and.i.ge.nbeta)&
              scr(2)=beta(i-nbeta+1)
            be=scr(2)*sc
            scr(3)=4
            scr(4)=0
            scr(5)=nalpha
            scr(6)=0
            do j=1,nalpha
               if (isym.eq.0) then
                  if (ilog.eq.0) then
                     scr(6+j)=ssm(i,j,nt)*exp(-be/2)
                     if (scr(6+j).ge.small) then
                        scr(6+j)=sigfig(scr(6+j),7,0)
                     else
                        scr(6+j)=sigfig(scr(6+j),6,0)
                     endif
                  else
                     scr(6+j)=0
                     if (ssm(i,j,nt).gt.zero) then
                        scr(6+j)=log(ssm(i,j,nt))-be/2
                        scr(6+j)=sigfig(scr(6+j),7,0)
                     endif
                  endif
               else if (isym.eq.1) then
                  if (i.lt.nbeta) then
                     if (ilog.eq.0) then
                        scr(6+j)=ssm(nbeta-i+1,j,nt)*exp(be/2)
                        if (scr(6+j).ge.small) then
                           scr(6+j)=sigfig(scr(6+j),7,0)
                        else
                           scr(6+j)=sigfig(scr(6+j),6,0)
                        endif
                     else
                        scr(6+j)=tiny
                        if (ssm(nbeta-i+1,j,nt).gt.zero) then
                           scr(6+j)=log(ssm(nbeta-i+1,j,nt))+be/2
                           scr(6+j)=sigfig(scr(6+j),7,0)
                        endif
                     endif
                  else
                     if (ilog.eq.0) then
                        scr(6+j)=ssp(i-nbeta+1,j,nt)*exp(be/2)
                        if (scr(6+j).ge.small) then
                           scr(6+j)=sigfig(scr(6+j),7,0)
                        else
                           scr(6+j)=sigfig(scr(6+j),6,0)
                        endif
                     else
                        scr(6+j)=tiny
                        if (ssp(i-nbeta+1,j,nt).gt.zero) then
                           scr(6+j)=log(ssp(i-nbeta+1,j,nt))+be/2
                           scr(6+j)=sigfig(scr(6+j),7,0)
                        endif
                     endif
                  endif
               else if (isym.eq.2) then
                  if (ilog.eq.0) then
                     scr(6+j)=ssm(i,j,nt)
                     if (scr(6+j).ge.small) then
                        scr(6+j)=sigfig(scr(6+j),7,0)
                     else
                        scr(6+j)=sigfig(scr(6+j),6,0)
                     endif
                  else
                     scr(6+j)=tiny
                     if (ssm(i,j,nt).gt.zero) then
                        scr(6+j)=log(ssm(i,j,nt))
                        scr(6+j)=sigfig(scr(6+j),7,0)
                     endif
                  endif
               else if (isym.eq.3) then
                  if (i.lt.nbeta) then
                     if (ilog.eq.0) then
                        scr(6+j)=ssm(nbeta-i+1,j,nt)
                        if (scr(6+j).ge.small) then
                           scr(6+j)=sigfig(scr(6+j),7,0)
                        else
                           scr(6+j)=sigfig(scr(6+j),6,0)
                        endif
                     else
                        scr(6+j)=tiny
                        if (ssm(nbeta-i+1,j,nt).gt.zero) then
                           scr(6+j)=log(ssm(nbeta-i+1,j,nt))
                           scr(6+j)=sigfig(scr(6+j),7,0)
                        endif
                     endif
                  else
                     if (ilog.eq.0) then
                        scr(6+j)=ssp(i-nbeta+1,j,nt)
                        if (scr(6+j).ge.small) then
                           scr(6+j)=sigfig(scr(6+j),7,0)
                        else
                           scr(6+j)=sigfig(scr(6+j),6,0)
                        endif
                     else
                        scr(6+j)=tiny
                        if (ssp(-nbeta+1,j,nt).gt.zero) then
                           scr(6+j)=log(ssp(i-nbeta+1,j,nt))
                           scr(6+j)=sigfig(scr(6+j),7,0)
                        endif
                     endif
                  endif
               endif
               if (ilog.eq.0.and.scr(6+j).lt.smin) scr(6+j)=0
            enddo
            call listio(0,nout,nprnt,scr(1),nb,nw)
            ll=1
            do while (nb.ne.0)
               ll=ll+nw
               call moreio(0,nout,nprnt,scr(ll),nb,nw)
            enddo
         endif
      enddo
   enddo
   if (nss.ne.0.and.b7.le.zero) then
      scr(1)=0
      scr(2)=0
      scr(3)=0
      scr(4)=0
      scr(5)=1
      ntf=ntempr
      scr(6)=ntf
      scr(7)=ntf
      scr(8)=2
      do i=1,ntf
         if (i.le.ntempr) then
            scr(2*i+7)=sigfig(tempr(i),7,0)
            scr(2*i+8)=sigfig(tempf1(i),7,0)
         else
            scr(2*i+7)=sigfig(up*scr(2*i+5),7,0)
            scr(2*i+8)=scr(2*i+6)
         endif
      enddo
      call tab1io(0,nout,nprnt,scr(1),nb,nw)
   endif
   scr(1)=0
   scr(2)=0
   scr(3)=0
   scr(4)=0
   scr(5)=1
   ntf=ntempr
   scr(6)=ntf
   scr(7)=ntf
   scr(8)=2
   do i=1,ntf
      if (i.le.ntempr) then
         scr(2*i+7)=sigfig(tempr(i),7,0)
         scr(2*i+8)=sigfig(tempf(i),7,0)
      else
         scr(2*i+7)=sigfig(up*scr(2*i+5),7,0)
         scr(2*i+8)=scr(2*i+6)
      endif
   enddo
   call tab1io(0,nout,nprnt,scr(1),nb,nw)
   call asend(nout,nprnt)
   call afend(nout,nprnt)
   call amend(nout,nprnt)
   call atend(nout,nprnt)
   return
   end subroutine endout

end module leapm

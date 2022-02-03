module groupm
   ! provides subroutine groupr for NJOY2016
   use locale
   implicit none
   private
   public groupr

   ! global variables

   ! particle group structure
   integer::ign,ngn
   real(kr),dimension(:),allocatable::egn

   ! gamma group structure
   integer::igg,ngg
   real(kr),dimension(:),allocatable::egg

   ! weighting function
   integer::iwt,jsigz,jtemp
   real(kr),dimension(:),allocatable::wght
   real(kr),dimension(:),allocatable::wtbuf
   integer::nfp,nfv,nfscr,nflmax
   real(kr)::felo,fehi,sigpot,alpha2,alpha3,beta,sam,gamma
   integer,parameter::nbuf=1000

   ! temperatures and background cross sections
   real(kr),dimension(:),allocatable::temp
   real(kr),dimension(:),allocatable::sigz
   integer::ntemp,nsigz

   ! input mat, legendre order, print control, smoothing option and run title
   integer::matb,lord,iprint,ismooth
   real(kr)::rtitle(17)
   character(4)::title(17)
   equivalence(title(1),rtitle(1))

   ! unit numbers
   integer::nendf,npend,nend2,nend3,nflx,npend2,ninwt
   integer::ngout1,ngout2,ntw

   ! auto reaction processing lists
   ! - lfs8(i) points to the "level number" from mf8.
   ! - mlfs8(i) is calculated and corresponds to NJOY's assumption
   !   of the ground state or isomer number.
   integer,parameter::maxr1=500
   integer,parameter::maxr2=500
   integer::mf4(maxr1),mf6(maxr1),mf12(maxr1),mf13(maxr1),mf18(maxr1),&
     mf4r(6,maxr1),mf6p(6,maxr1),mf10f(maxr2),mf10s(maxr2),mf10i(maxr2),&
     lfs8(maxr2),mlfs8(maxr2)
   integer::lastza,izatest

   ! reaction and kinematics parameters
   integer::matd,mfd,mtd
   real(kr)::econst
   integer::jconst
   integer::ndelg
   real(kr)::dntc(8)
   real(kr)::awr,q,thresh,alpha
   integer::lrflag
   real(kr)::awrp
   integer::izap
   real(kr)::xtot
   real(kr)::zap,awp,aprime
   integer::lip,law,lct
   integer::izat
   real(kr)::spi
   real(kr)::emaxx,ebeg
   integer::lfs,isom
   integer::izar

   ! unresolved resonance parameters
   integer::nunr,intunr,lrp,lssf
   real(kr),dimension(:),allocatable::unr

   ! temporary arrays for reaction data
   real(kr),dimension(:),allocatable::yield
   real(kr),dimension(:),allocatable::sigma
   real(kr),dimension(:),allocatable::fls
   real(kr),dimension(:),allocatable::flgn
   real(kr),dimension(:),allocatable::gyln
   real(kr),dimension(:),allocatable::gyle
   real(kr),dimension(:),allocatable::eyle
   real(kr),dimension(:,:),allocatable::gfle
   real(kr),dimension(:),allocatable::sedist
   real(kr),dimension(:,:),allocatable::sede
   real(kr),dimension(:),allocatable::ddmf6
   real(kr),dimension(:),allocatable::aes
   real(kr),dimension(:,:),allocatable::stmp
   real(kr),dimension(:,:),allocatable::slst
   real(kr),dimension(:,:),allocatable::ftmp
   real(kr),dimension(:,:),allocatable::flst
   real(kr),dimension(:),allocatable::falo,fahi
   integer::ipan

contains

   subroutine groupr
   !-------------------------------------------------------------------
   !
   ! compute self-shielded group-averaged cross sections
   !
   ! Produces self-shielded cross sections, neutron scattering
   ! matrices, and photon production matrices.  Scattering and
   ! photon matrices may be self-shielded if desired (see init).
   ! Bondarenko weighting is normally used.   Optionally, the flux
   ! can be computed for an infinite mixture of heavy absorber
   ! and light moderator.  Delayed neutron data and thermal
   ! scattering matrices are handled specially.
   !
   ! The integration over initial energy is handled in the same
   ! way for all reaction types by using the integrand
   !                   feed*xsec*flux                   .
   ! Feed is the source into final energy group gprime and
   ! Legendre order l from initial energy e (see getff).  For
   ! vectors, the feed is 1. or a yield (nubar, mubar).  For two
   ! body scattering, a center-of-mass Gaussian integration is used
   ! to obtain accurate results even for small Legendre components
   ! of the group-to-group scattering.  Additional initial energy
   ! quadrature points are added to integrate the known polynomial
   ! order of this feed function.  Feed for tabulated continuum
   ! reactions is computed exactly on the endf grid points and
   ! then interpolated at e.  A special projection interpolation
   ! scheme is used for thermal matrices (see getaed).  The feed
   ! for analytic continuum reactions is exact.
   !
   !---input specifications (free format)---------------------------
   !
   ! card1
   !    nendf   unit for endf tape
   !    npend   unit for pendf tape
   !    ngout1  unit for input gout tape (default=0)
   !    ngout2  unit for output gout tape (default=0)
   ! card2
   !    matb    material to be processed
   !             if ngout=0, matb<0 is an option to automatically
   !              process all the mats on the endf input tape.
   !             otherwise, matb<0 is a flag to add mts to and/or
   !              replace individual mts on gout1.
   !    ign     neutron group structure option
   !    igg     gamma group structure option
   !    iwt     weight function option
   !    lord    legendre order
   !    ntemp   number of temperatures (default=1)
   !    nsigz   number of sigma zeroes (default=1)
   !    iprint  long print option (0/1=minimum/maximum)
   !            (default=1)
   !    ismooth switch on/off smoothing operation (1/0, default=1=on)
   !            set ismooth to 1 to enable sqrt(e) smoothing for
   !            mf6 cm emission spectra at low energies and for
   !            histogram delayed neutron spectra at low energies.
   ! card3
   !    title   run label (up to 80 characters delimited by quotes,
   !            ended with /)  (default=blank)
   ! card4
   !    temp    temperatures in kelvin
   ! card5
   !    sigz    sigma zero values (including infinity)
   !
   !          if ign=1, read neutron group structure (6a and 6b)
   ! card6a
   !    ngn     number of groups
   ! card6b
   !    egn     ngn+1 group breaks (ev)
   !
   !          if igg=1, read gamma group structure (7a and 7b)
   ! card7a
   !    ngg     number of groups
   ! card7b
   !    egg     ngg+1 group breaks (ev)
   !
   !          weight function options (8a,8b,8c,8d)
   ! card8a     flux calculator parameters (iwt.lt.0 only)
   !    fehi    break between computed flux and bondarenko flux
   !            (must be in the resolved resonance range)
   !    sigpot  estimate of potential scattering cross section
   !    nflmax  maximum number of computed flux points
   !    ninwt   tape unit for new flux parameters (default=0)
   !            note: weighting flux file is always written binary
   !    jsigz   index of reference sigma zero in sigz array
   !            (default=0)
   !    alpha2   alpha for admixed moderator (def=o=none)
   !    sam      admixed moderator xsec in barns per absorber
   !             atom (def=0=none)
   !    beta     heterogeneity parameter (def=0=none)
   !    alpha3   alpha for external moderator (def=0=none)
   !    gamma    fraction of admixed moderator cross section in
   !              external moderator cross section (def=0)
   ! card8b     tabulated (iwt=1 or -1 only)
   !    wght    read weight function as tab1 record,
   !            this may span multiple lines and ends with a /.
   ! card8c     analytic flux parameters (iwt=4 or -4 only)
   !    eb      thermal break (ev)
   !    tb      thermal temperature (ev)
   !    ec      fission break (ev)
   !    tc      fission temperature (ev)
   ! card8d     input resonance flux (iwt=0 only)
   !    ninwt   tape unit for flux parameters (binary)
   !
   ! card9
   !    mfd     file to be processed
   !    mtd     section to be processed
   !    mtname  description of section to be processed
   !          repeat for all reactions desired
   !          mfd=0/ terminates this temperature/material.
   ! card10
   !    matd    next mat number to be processed
   !            matd=0/ terminates groupr run.
   !
   !---options for input variables----------------------------------
   !
   !     ign          meaning
   !     ---          -------
   !      1           arbitrary structure (read in)
   !      2           csewg 239-group structure
   !      3           lanl 30-group structure
   !      4           anl 27-group structure
   !      5           rrd 50-group structure
   !      6           gam-i 68-group structure
   !      7           gam-ii 100-group structure
   !      8           laser-thermos 35-group structure
   !      9           epri-cpm 69-group structure
   !     10           lanl 187-group structure
   !     11           lanl 70-group structure
   !     12           sand-ii 620-group structure
   !     13           lanl 80-group structure
   !     14           eurlib 100-group structure
   !     15           sand-iia 640-group structure
   !     16           vitamin-e 174-group structure
   !     17           vitamin-j 175-group structure
   !     18           xmas nea-lanl
   !     all new additional group structure with 7 significant
   !     decimal digits compatible with calendf
   !     19           ecco  33-group structure
   !     20           ecco 1968-group structure
   !     21           tripoli 315-group structure
   !     22           xmas lwpc 172-group structure
   !     23           vit-j lwpc 175-group structure
   !     24           shem cea 281-group structure
   !     25           shem epm 295-group structure
   !     26           shem cea/epm 361-group structure
   !     27           shem epm 315-group structure
   !     28           rahab aecl 89-group structure
   !     29           ccfe   660-group structure  (30 MeV)
   !     30           ukaea 1025-group structure  (30 MeV)
   !     31           ukaea 1067-group structure (200 MeV)
   !     32           ukaea 1102-group structure   (1 GeV)
   !     33           ukaea  142-group structure (200 MeV)
   !     34           lanl 618-group structure
   !
   !     igg          meaning
   !     ---          -------
   !      0           none
   !      1           arbitrary structure (read in)
   !      2           csewg 94-group structure
   !      3           lanl 12-group structure
   !      4           steiner 21-group gamma-ray structure
   !      5           straker 22-group structure
   !      6           lanl 48-group structure
   !      7           lanl 24-group structure
   !      8           vitamin-c 36-group structure
   !      9           vitamin-e 38-group structure
   !     10           vitamin-j 42-group structure
   !
   !     iwt          meaning
   !     ---          -------
   !      1           read in smooth weight function
   !      2           constant
   !      3           1/e
   !      4           1/e + fission spectrum + thermal maxwellian
   !      5           epri-cell lwr
   !      6           (thermal) -- (1/e) -- (fission + fusion)
   !      7           same with t-dep thermal part
   !      8           thermal--1/e--fast reactor--fission + fusion
   !      9           claw weight function
   !     10           claw with t-dependent thermal part
   !     11           vitamin-e weight function (ornl-5505)
   !     12           vit-e with t-dep thermal part
   !     -n           compute flux with weight n
   !      0           read in resonance flux from ninwt
   !
   !     mfd          meaning
   !     ---          -------
   !      3           cross section or yield vector
   !      5           fission chi by short-cut method
   !      6           neutron-neutron matrix (mf4/5)
   !      8           neutron-neutron matrix (mf6)
   !     12           photon prod. xsec (photon yields given, mf12)
   !     13           photon prod. xsec (photon xsecs given, mf13)
   !     16           neutron-gamma matrix (photon yields given)
   !     17           neutron-gamma matrix (photon xsecs given)
   !     18           neutron-gamma matrix (mf6)
   !         note: if necessary, mfd=13 will automatically change
   !         to 12 and mfd=16 will automatically change to 17 or 18.
   !     21           proton production matrix (mf6)
   !     22           deuteron production (mf6)
   !     23           triton production (mf6)
   !     24           he-3 production (mf6)
   !     25           alpha production (mf6)
   !     26           residual nucleus (a>4) production (mf6)
   !     31           proton production matrix (mf4)
   !     32           deuteron production (mf4)
   !     33           triton production (mf4)
   !     34           he-3 production (mf4)
   !     35           alpha production (mf4)
   !     36           residual nucleus (a>4) production (mf4)
   !          note: if necessary, mfd=21-26 will
   !          automatically change to 31-36.
   !    1zzzaaam       nuclide production for zzzaaam
   !                     subsection from file 3
   !    2zzzaaam       nuclide production for zzzaaam
   !                     subsection from file 6
   !    3zzzaaam       nuclide production for zzzaaam
   !                     subsection from file 9
   !    4zzzaaam       nuclide production for zzzaaam
   !                     subsection from file 10
   !    40000000       fission product production (mtd=18 only)
   !                     subsection from file 10
   !
   !     mtd          meaning
   !     ---          -------
   !     -n           process all mt numbers from the previous
   !                          entry to n inclusive
   !     221-250      reserved for thermal scattering
   !     257          average energy
   !     258          average lethargy
   !     259          average inverse velocity (m/sec)
   !
   !     automatic reaction processing options
   !     -------------------------------------
   !        3/        do all reactions in file3 on input pendf
   !        6/        do all matrix reactions in endf dictionary
   !       10/        do all isotope productions using mf8
   !       13/        do all photon production cross sections
   !       16/        do all photon production matrices
   !       21/        do all proton production matrices
   !       22/        do all deuteron production matrices
   !       23/        do all triton production matrices
   !       24/        do all he-3 production matrices
   !       25/        do all alpha production matrices
   !       26/        do all a>4 production matrices
   !
   !-------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso,nsyse
   use endf   ! provides endf routines and variables
   use util   ! provides openz,timer,repoz,skiprz,error,mess,closz
   ! internals
   integer::nwscr,nb,nw,itend,nwds,i,itemp,nz
   integer::ngnp1,nggp1,np,loc,ngi,mtdp,iauto,izam,j,ibase
   integer::iaddmt,iglo,nq,ig1,ig,nll,ig2lo,it,iz,il,igzero
   integer::naed,ng2,nl,idis,lim,nlg,ng2g,mfdn,jzam
   integer::itmp
   real(kr)::time,za,tempin,eps,diff,ee,first,en
   real(kr)::enext,elo,ehi,yld,test
   real(kr),dimension(:,:),allocatable::flux,sig
   character(4)::mtname(15)
   character(60)::strng
   character(66)::text
   character(4)::tz(17)
   real(kr)::z(17)
   equivalence(tz(1),z(1))
   real(kr),dimension(:),allocatable::scr
   real(kr),dimension(:),allocatable::flxc
   real(kr),dimension(:,:,:),allocatable::ans
   real(kr),dimension(:,:),allocatable::ff
   real(kr),dimension(:,:,:),allocatable::prod
   real(kr),parameter::small=1.e-6_kr
   real(kr),parameter::zero=0

   !--initialize
   nwscr=10000
   nfscr=10
   nflx=11
   nend2=12
   nend3=13
   npend2=0
   ninwt=0
   ipan=0
   call openz(-nflx,1)
   allocate(wtbuf(nbuf))
   ebeg=10*small

   !--write heading and read user input
   call timer(time)
   write(nsyso,'(/'' groupr...'',&
     &''compute self-shielded group-averaged cross-sections'',&
     &8x,f8.1,''s'')') time
   write(nsyse,'(/'' groupr...'',59x,f8.1,''s'')') time
   call ruinb(iaddmt)
   nwscr=max(nwscr,nsigz+ngn+ngg+9)
   allocate(scr(nwscr))
   if (nendf.lt.0) nend2=-nend2
   if (nendf.lt.0) nend3=-nend3
   call openz(nend2,1)
   call openz(nend3,1)
   if (npend.lt.0) nfscr=-nfscr
   call openz(nfscr,1)
   naed=14
   if (npend.lt.0) naed=-naed
   call openz(naed,1)

   !--determine endf format being used
   call tpidio(nendf,0,0,scr,nb,nw)
   call contio(nendf,0,0,scr,nb,nw)
   call contio(nendf,0,0,scr,nb,nw)
   if (n1h.ne.0) then
      iverf=4
   else if (n2h.eq.0) then
      iverf=5
   else
      iverf=6
   endif
   write(nsyso,'(/'' using endf-'',i1,'' format'')') iverf
   if (iverf.le.4) then
      emaxx=1.5e7_kr
   else if (iverf.eq.5) then
      emaxx=2.0e7_kr
   else
      call contio(nendf,0,0,scr,nb,nw)
      emaxx=2.0e7_kr
      if (c2h.gt.zero) emaxx=c2h
   endif

   !--locate position for new material on old gout tape.
   call repoz(ngout1)
   call repoz(ngout2)
   math=0
   mfh=0
   mth=0
   nsh=0
   itend=0
   nwds=17
   if (ngout1.eq.0) call tpidio(0,ngout2,0,rtitle,nb,nwds)
   ntw=1
   rtitle(1)=0
   call repoz(nendf)
   call repoz(npend)
   call tpidio(nendf,0,0,scr,nb,nw)
   call tpidio(npend,0,0,scr,nb,nw)
   if (ngout1.eq.0.or.ngout2.eq.0) go to 170
   call tpidio(ngout1,ngout2,0,scr,nb,nwds)
   if (iaddmt.gt.0) go to 170
  120 continue
   call contio(ngout1,0,0,scr,nb,nw)
   if (math.eq.-1) go to 160
   if (math.eq.1) go to 120
   if (math.eq.matb) go to 140
   if (math.gt.matb) go to 160
   call contio(0,ngout2,0,scr,nb,nw)
   call tomend(ngout1,ngout2,0,scr)
   go to 120
  140 continue
   call tomend(ngout1,0,0,scr)
   go to 170
  160 continue
   call skiprz(ngout1,-1)
  170 continue

   !--locate desired material on endf tape
   !--make two copies for use by retrieval routines.
   call repoz(nendf)
   if (iaddmt.lt.0) then
      call contio(nendf,0,0,scr,nb,nw)
      if (math.lt.0) go to 650
      matb=math
      call skiprz(nendf,-1)
   else
      call findf(matb,1,451,nendf)
   endif
   call repoz(nend2)
   call repoz(nend3)
   nsh=1
   nsc=1
   math=1
   call afend(nend2,nend3)
   call conver(nendf,nend2,nend3)
   call atend(nend2,nend3)
   call repoz(nend2)
   call repoz(nend3)

   !--search for desired mat and temperatures on pendf tape
   !--store unresolved resonance cross sections if present.
   do 630 itemp=1,ntemp
   call findf(matb,1,451,npend)
   call contio(npend,0,0,scr,nb,nw)
   lrp=l1h
   jtemp=itemp
   tempin=temp(itemp)
   eps=tempin/10000
   diff=100000000
   do while (diff.gt.1+eps)
      za=c1h
      awr=c2h
      if (iverf.ge.5) call contio(npend,0,0,scr,nb,nw)
      if (iverf.ge.6) call contio(npend,0,0,scr,nb,nw)
      call hdatio(npend,0,0,scr,nb,nw)
      diff=abs(c1h-tempin)
      if (diff.gt.1+eps) then
         if (c1h.gt.tempin) then
            write(strng,'(''unable to find mat='',i4,'' t='',&
              &1p,e12.4, 1x, e12.4)') matb,tempin, c1h
            call error('groupr',strng,' ')
         endif
         call tomend(npend,0,0,scr)
         call contio(npend,0,0,scr,nb,nw)
         if (math.ne.matb) then
            write(strng,'(&
              &''unable to find mat='',i4,'' t='',1p,e12.4, 1x, i4)')&
              matb,tempin, math
            call error('groupr',strng,' ')
         endif
      else
         if (diff.gt.small) then
            write(nsyso,'(&
              &''   difference between temperatures desired'',&
              &''  and found is '',1p,e10.2)') diff
            write(nsyse,'(&
              &''   difference between temperatures desired'',&
              &''  and found is '',1p,e10.2)') diff
         endif
         do i=1,17
            z(i)=scr(6+i)
         enddo
         write(nsyso,'(/'' processing mat '',i6/&
           &1x,9(''----''),''--''/2x,17a4)') matb,(tz(i),i=1,17)
         write(nsyso,'('' '')')
         call stounr(matb,tempin,npend)
      endif
   enddo

   !--generate weighting flux on nflx.
   if (nsigz.gt.1) then
      nw=nflmax*(3+nsigz)
      if (nw.eq.0) nw=1
      allocate(flxc(nw))
      call genflx(flxc)
      deallocate(flxc)
   endif

   !--write head record for this material on gout tape.
   if (ngout2.ne.0) then
      if (iaddmt.eq.0) then
         nz=nsigz
         nsh=1
         math=matb
         mfh=1
         mth=451
         scr(1)=za
         scr(2)=awr
         scr(3)=0
         scr(4)=nz
         scr(5)=-1
         scr(6)=ntw
         call contio(0,ngout2,0,scr,nb,nw)
         scr(1)=tempin
         scr(2)=0
         scr(3)=ngn
         scr(4)=ngg
         scr(5)=0
         scr(6)=0
         nw=6
         do i=1,ntw
            nw=nw+1
            scr(nw)=rtitle(i)
         enddo
         do i=1,nz
            nw=nw+1
            scr(nw)=sigz(i)
         enddo
         ngnp1=ngn+1
         do i=1,ngnp1
            nw=nw+1
            scr(nw)=egn(i)
         enddo
         nggp1=ngg+1
         do i=1,nggp1
            nw=nw+1
            scr(nw)=egg(i)
         enddo
         if (ngg.eq.0) scr(nw)=0
         np=nw-6
         scr(5)=np
         loc=1
         call listio(0,ngout2,0,scr(loc),nb,nw)
         do while (nb.ne.0)
            loc=loc+nw
            call moreio(0,ngout2,0,scr(loc),nb,nw)
         enddo
         call afend(ngout2,0)
      else
         call contio(ngout1,ngout2,0,scr,nb,nw)
         call tofend(ngout1,ngout2,0,scr)
      endif
   endif

   !--main entry point for loop over reactions
   !--initialize and write appropriate heading.
   !--write head record on output tape.
   ngi=ngn
   if (izap.eq.0) ngi=ngg
   matd=matb
   mtd=1
   mtdp=1
   iauto=0
  360 continue
   lfs=0
   isom=0
   izam=0
   if (mtd+mtdp.lt.0) go to 400
   if (iauto.eq.0) go to 365
   call nextr(iauto,matd,mfd,mtdp,scr)
   if (mtdp.ne.0) go to 390
  365 continue
   mtdp=-1000
   strng=' '
   read(nsysi,*) mfd,mtdp,strng
   if (mfd.lt.0.or.mfd.eq.1.or.mfd.eq.2.or.mfd.eq.4) go to 381
   if (mfd.eq.7.or.mfd.eq.9.or.mfd.eq.11) go to 381
   if (mfd.eq.14.or.mfd.eq.15) go to 381
   if (mfd.gt.18.and.mfd.lt.21) go to 381
   if (mfd.gt.26.and.mfd.lt.31) go to 381
   if (mfd.gt.36.and.mfd.lt.10000000) go to 381
   if (mfd.ge.12.and.mfd.le.18.and.igg.eq.0)&
     call error('groupr','photons not allowed with igg=0.',' ')
   if (mfd.eq.0) go to 590
   if (mtdp.eq.-1000) go to 382
   read(strng,'(15a4)') (mtname(i),i=1,15)
   go to 390
  381 continue
   write(strng,'(''illegal mfd ='',i3)') mfd
   call error('groupr',strng,' ')
  382 continue
   iauto=1
   mtdp=0
   call nextr(iauto,matd,mfd,mtdp,scr)
   if (mtdp.ne.0) go to 390
   write(strng,'(''auto finds no reactions for mf='',i3)') mfd
   call mess('groupr',strng,' ')
   go to 365
  390 continue
   if (mtdp.ge.0) mtd=mtdp
   ! -- if auto processing, already have level number (lfs) and isomer
   !    number (isom) from nextr; if user input may need to decode the
   !    user mfd.  then calculate izam
   if (iauto.eq.0.and.mfd.ge.10000000) then
      ! -- decode from user mfd input
      if (mfd.eq.40000000) then ! fission special case for mf10
         lfs=0
         isom=0
         izam=-1
      else
         itmp=mfd/10000000
         itmp=(mfd-10000000*itmp)/10
         izar = itmp
         lfs=mfd-(10000000*(mfd/10000000)+10*itmp)
         isom=lfs
         if (lfs.lt.10) then
            izam=mod(mfd,10000000)
         else
            izam=10*mod(mfd,10000000)
         endif
      endif
   elseif (iauto.eq.1.and.mfd.ge.10000000) then
      if (mfd.eq.40000000) then ! fission special case for mf10
         izam=-1
      else
         izar=(mfd-((mfd/1000000)*1000000))/10
         if (lfs.lt.10) then
            izam=mod(mfd,10000000)+lfs
         else
            izam=10*mod(mfd,10000000)+lfs
         endif
      endif
   endif
  400 continue
   if (mtdp.lt.0) mtd=mtd+1
   if (iauto.eq.0.and.(mfd.eq.13.or.mfd.eq.16)) call mfchk(mfd,mtd)
   if (iauto.eq.0.and.mfd.ge.21.and.mfd.le.26) call mfchk2(mfd,mtd)
   if (iauto.eq.1) call namer(izap,izam,mfd,mtd,mtname)
   if (iaddmt.le.0) go to 490
  410 continue
   call contio(ngout1,0,0,scr,nb,nw)
   if (math.eq.0.and.mfh.eq.0.and.mth.eq.0) go to 430
   if (mfh.eq.0.or.mth.eq.0) go to 410
   if (mfh.ne.mfd.or.mth.ne.mtd) go to 420
   call tosend(ngout1,0,0,scr)
   go to 490
  420 continue
   call contio(0,ngout2,0,scr,nb,nw)
   call tosend(ngout1,ngout2,0,scr)
   go to 410
  430 continue
   call skiprz(ngout1,-1)
  490 continue
   call timer(time)
   call skiprz(npend,-1)
   call init(ng2g,nlg,nz)
   allocate(ans(nlg,nz,ng2g))
   ng2=ng2g
   nl=nlg
   allocate(ff(nlg,ng2g))
   allocate(prod(nlg,nz,ngi))
   allocate(flux(nz,nl))
   allocate(sig(nz,nl))
   ee=0
   if (mfd.ne.5) then
      call getsig(ee,first,idis,sig,nl,nz)
      call getflx(ee,en,idis,flux,nl,nz)
   endif
   call getff(ee,en,idis,yld,ff,nl,ng2,iglo,nq,nlg,ng2g)
   test=-1
   if (en.eq.test) go to 582
   ig1=0
   ig=0

   !--spectrum calculation (mfd=5 or constant spectra)
   nll=1
   if (mfd.eq.5) go to 502
   if (econst.eq.zero) go to 510
   go to 505
  502 continue
   if (mtd.eq.455) then
      do i=1,ndelg
         do iz=1,nz
            ans(i,iz,1)=dntc(i)
         enddo
      enddo
      nll=ndelg
   endif
  505 continue
   call getff(en,enext,idis,yld,ff,nll,ng2,ig2lo,nq,nlg,ng2g)
   do il=1,nll
      do i=1,ng2
         j=i
         if (mtd.eq.455) j=i+1
         do iz=1,nz
            ans(il,iz,j)=ff(il,i)
         enddo
      enddo
   enddo
   if (mfd.eq.5) ig=ngi
   if (mtd.eq.455) ng2=ng2+1
   go to 540

   !--loop over initial energy groups
  510 continue
   ig=ig+1
   if (ngi.eq.ngn) then
      elo=egn(ig)
      ehi=egn(ig+1)
   else
      elo=egg(ig)
      ehi=egg(ig+1)
   endif
   ig2lo=0
   ng2=2
   if (ehi.le.first) go to 580
   if (elo.ge.emaxx.and.ig.ne.ngi) go to 580
   enext=ehi
   do it=1,ng2g
      do iz=1,nz
         do il=1,nl
            ans(il,iz,it)=0
         enddo
      enddo
   enddo
  530 continue
   call panel(elo,enext,ans,ff,nl,nz,ng2,ig2lo,nlg,ng2g)
   if (enext.eq.ehi) go to 540
   elo=enext
   enext=ehi
   go to 530
  540 continue

   !--accumulate production below econst, and print it
   if (econst.gt.0.and.mfd.ne.5.and.ig.ne.0.and.ig.le.jconst) then
      do iz=1,nz
         do il=1,nl
            ans(il,iz,2)=ans(il,iz,2)/ans(il,iz,1)
            prod(il,iz,ig)=ans(il,iz,2)
         enddo
      enddo
      ig2lo=0
      if (ig.ge.jconst) then
         call displa(ig,prod,nl,nz,jconst,ig2lo,igzero,nlg,ngi)
      endif

   !--print results for this initial energy group
   else
      if (ig1.le.0) then
         ig1=1
         if (tempin.ne.0) write(nsyso,'(&
           &'' group constants at t='',1p,e9.3,'' deg k'',&
           &32x,0p,f8.1,''s'')') tempin,time
         if (tempin.eq.0) write(nsyso,'(/&
           &'' group constants at t=zero deg k'',37x,f8.1,''s'')')&
           time
         if (mfd.lt.1000000) then
            write(nsyso,'('' for mf'',i3,'' and mt'',i3,1x,15a4)')&
              mfd,mtd,(mtname(i),i=1,15)
         else
            mfdn=mfd/10000000
            jzam=izam
            if (mfdn.eq.1) then
               write(nsyso,'(&
                 &'' for mf3 mt'',i3,'' zam'',i8,1x,15a4)')&
                 mtd,jzam,(mtname(i),i=1,15)
            else if (mfdn.eq.2) then
               write(nsyso,'(&
                 &'' for mf3*mf6 mt'',i3,'' zam'',i8,1x,15a4)')&
                 mtd,jzam,(mtname(i),i=1,15)
            else if (mfdn.eq.3) then
               write(nsyso,'(&
                 &'' for mf3*mf9 mt'',i3,'' zam'',i8,1x,15a4)')&
                 mtd,jzam,(mtname(i),i=1,15)
            else if (mfdn.eq.4) then
               write(nsyso,'(&
                 &'' for mf10 mt'',i3,'' zam'',i8,1x,15a4)')&
                 mtd,jzam,(mtname(i),i=1,15)
            endif
         endif
         if (lrflag.gt.0) write(nsyso,'(14x,''lr'',i3,&
           &'' particle emission'')') lrflag
         if (econst.gt.zero) write(nsyso,'(/&
           &'' spectrum constant below '',1p,e10.4,'' ev'',&
           &''  ('',i4,'' groups)'')') econst,jconst
         math=matb
         mfh=mfd
         if (mfd.eq.8) mfh=6
         if (mfd.eq.17) mfh=16
         if (mfd.eq.18) mfh=16
         if (mfd.ge.31.and.mfd.le.36) mfh=mfd-10
         if (mfd.ge.10000000) mfh=3
         mth=mtd
         scr(1)=za
         scr(2)=izam
         scr(3)=nl
         scr(4)=nz
         scr(5)=lrflag
         scr(6)=ngi
         call contio(0,ngout2,0,scr,nb,nw)
         call displa(-1,ans,nl,nz,ng2,ig2lo,igzero,nlg,ng2g)
      endif
      call displa(ig,ans,nl,nz,ng2,ig2lo,igzero,nlg,ng2g)
   endif

   !--write group constants on output tape
   lim=nl*nz*ng2
   if (ngout2.ne.0) then
      if (igzero.ne.0.or.ig.eq.ngi) then
         math=matb
         mfh=mfd
         if (mfd.eq.8) mfh=6
         if (mfd.eq.17) mfh=16
         if (mfd.eq.18) mfh=16
         if (mfd.ge.31.and.mfd.le.36) mfh=mfd-10
         if (mfd.ge.10000000) mfh=3
         mth=mtd
         scr(1)=tempin
         scr(2)=0
         scr(3)=ng2
         scr(4)=ig2lo
         scr(5)=lim
         scr(6)=ig
         i=0
         j=6
         ibase=6
         do it=1,ng2
            do iz=1,nz
               do il=1,nl
                  i=i+1
                  j=j+1
                  scr(j)=ans(il,iz,it)
                  if (j.ge.(npage+ibase).or.i.eq.lim) then
                     if (ibase.ne.0) then
                        call listio(0,ngout2,0,scr,nb,j)
                        ibase=0
                        j=0
                     else
                        call moreio(0,ngout2,0,scr,nb,nw)
                        j=0
                     endif
                  endif
               enddo
            enddo
         enddo
      endif
   endif
  580 continue
   if (ig.lt.ngi.and.mfd.ne.5) go to 510
   if (ig1.eq.0) write(nsyso,'(/&
     &'' threshold is above highest energy bound for mt'',i3)') mtd
   write(nsyso,'('' '')')
   go to 586
  582 continue
   write(nsyso,'(/'' no photons in mf6, mt'',i3)') mth

   !--loop over other desired reactions
  586 continue
   deallocate(ans)
   deallocate(ff)
   deallocate(prod)
   deallocate(flux)
   deallocate(sig)
   if (allocated(sigma)) deallocate(sigma)
   if (allocated(yield)) deallocate(yield)
   if (allocated(fls)) deallocate(fls)
   if (allocated(flgn)) deallocate(flgn)
   if (allocated(gyln)) deallocate(gyln)
   if (allocated(gyle)) deallocate(gyle)
   if (allocated(eyle)) deallocate(eyle)
   if (allocated(gfle)) deallocate(gfle)
   if (allocated(sedist)) deallocate(sedist)
   if (allocated(sede)) deallocate(sede)
   if (allocated(ddmf6)) deallocate(ddmf6)
   if (allocated(aes)) deallocate(aes)
   if (mth.gt.0) call tosend(npend,0,0,scr)
   if (ig1.ne.0) call asend(ngout2,0)
   go to 360

   !--loop over desired temperatures
  590 continue
   if (allocated(unr)) deallocate(unr)
   if (iaddmt.gt.0) go to 600
   call amend(ngout2,0)
   if (itemp.eq.ntemp) go to 630
   call tomend(npend,0,0,scr)
   call contio(npend,0,0,scr,nb,nw)
   if (math.ne.matb)&
     call error('groupr','unable to find next temp.',' ')
   go to 630
   ! copy rest of material if necessary
  600 continue
   call contio(ngout1,0,0,scr,nb,nw)
   if (math.le.0) go to 610
   if (math.ne.matd) go to 620
  610 continue
   call contio(0,ngout2,0,scr,nb,nw)
   if (math.lt.0) itend=1
   if (math.le.0) go to 630
   call tomend(ngout1,ngout2,0,scr)
   go to 630
  620 continue
   call skiprz(ngout1,-2)
  630 continue
   call tomend(nendf,0,0,scr)

   !--loop over desired materials.
   if (iaddmt.lt.0) go to 160
   read(nsysi,*) matb
   iaddmt=0
   if (matb.lt.0.and.ngout1.ne.0) iaddmt=1
   if (matb.lt.0.and.ngout1.eq.0) iaddmt=-1
   matb=iabs(matb)
   if (matb.eq.0) go to 650
   if (ngout1.eq.0) go to 160
   go to 120

   !--groupr is complete
  650 continue
   if (ngout1.eq.0) call atend(ngout2,0)
   if (ngout1.ne.0.and.itend.eq.0)&
     call totend(ngout1,ngout2,0,scr)
   deallocate(egn)
   deallocate(egg)
   deallocate(wtbuf)
   if (allocated(wght)) deallocate(wght)
   if (allocated(temp)) deallocate(temp)
   if (allocated(sigz)) deallocate(sigz)
   if (allocated(stmp)) then
      deallocate(stmp)
      deallocate(slst)
      deallocate(ftmp)
      deallocate(flst)
   endif
   if (allocated(falo)) then
       deallocate(falo)
       deallocate(fahi)
   endif
   call repoz(nendf)
   call repoz(npend)
   call repoz(ngout1)
   call repoz(ngout2)
   call repoz(ninwt)
   call closz(nendf)
   call closz(npend)
   call closz(ngout1)
   call closz(ngout2)
   call closz(nend2)
   call closz(nend3)
   call closz(ninwt)
   call closz(naed)
   call closz(nfscr)
   call closz(-nflx)
   call timer(time)
   write(nsyso,'(69x,f8.1,''s''/&
     &1x,7(''**********''),''*******'')') time
   return
   end subroutine groupr

   subroutine ruinb(iaddmt)
   !-------------------------------------------------------------------
   ! Read user input for group averaging.
   !-------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use util   ! provides openz,error
   ! externals
   integer::iaddmt
   ! internals
   integer::i
   character(66)::text
   real(kr),parameter::sigzmx=1.e10_kr
   real(kr),parameter::zero=0

   !--read and display user input
   ngout1=0
   ngout2=0
   read(nsysi,*) nendf,npend,ngout1,ngout2
   call openz(nendf,0)
   call openz(npend,0)
   call openz(ngout1,0)
   call openz(ngout2,1)
   iprint=1
   ismooth=1
   ntemp=1
   nsigz=1
   read(nsysi,*) matb,ign,igg,iwt,lord,ntemp,nsigz,iprint,ismooth
   iaddmt=0
   if (matb.lt.0.and.ngout1.ne.0) iaddmt=1
   if (matb.lt.0.and.ngout1.eq.0) iaddmt=-1
   ntw=17
   text=' '
   read(nsysi,*) text
   read(text,'(16a4,a2)') (title(i),i=1,ntw)
   write(nsyso,'(/&
     &'' unit for endf tape ................... '',i10/&
     &'' unit for pendf tape .................. '',i10/&
     &'' unit for input gout tape ............. '',i10/&
     &'' unit for output gout tape ............ '',i10)')&
     nendf,npend,ngout1,ngout2
   write(nsyso,'(&
     &'' mat to be processed .................. '',i10/&
     &'' neutron group structure option ....... '',i10/&
     &'' gamma group option ................... '',i10/&
     &'' weight function option ............... '',i10/&
     &'' legendre order ....................... '',i10/&
     &'' print option (0 min, 1 max) .......... '',i10/&
     &'' smoothing option (0 off, 1 on) ....... '',i10)')&
     matb,ign,igg,iwt,lord,iprint,ismooth
   if (ismooth.ne.0.and.ismooth.ne.1) then
      call error('ruinb','illegal ismooth.',' ')
   endif
   write(nsyso,'(/'' run title''/&
     &1x,3(''----------''),''--------''/&
     &6x,16a4,a2)') (title(i),i=1,ntw)
   write(nsyso,'('' '')')
   if (allocated(temp)) then
      deallocate(temp)
      deallocate(sigz)
   endif
   allocate(temp(ntemp))
   allocate(sigz(nsigz))
   read(nsysi,*) (temp(i),i=1,ntemp)
   if (temp(1).eq.zero) write(nsyso,'(&
     &'' temperatures (kelvin) ................ '',6x,''zero'')')
   if (temp(1).ne.zero) write(nsyso,'(&
     &'' temperatures (kelvin) ................ '',1p,e10.2)')&
     temp(1)
   if (ntemp.gt.1) write(nsyso,'(40x,1p,e10.2)')&
     (temp(i),i=2,ntemp)
   read(nsysi,*) (sigz(i),i=1,nsigz)
   sigz(1)=sigzmx
   if (nsigz.gt.0) write(nsyso,'(&
     &'' sigma zeroes ......................... '',2x,''infinity'')')
   if (nsigz.gt.1) write(nsyso,'(40x,1p,e10.2)')&
     (sigz(i),i=2,nsigz)
   if (iaddmt.gt.0) write(nsyso,'(/&
     &'' replacing and/or adding mts'')')
   call gengpn
   call gengpg
   call genwtf
   return
   end subroutine ruinb

   subroutine nextr(iauto,matd,mfd,mtd,scr)
   !-------------------------------------------------------------------
   ! Find next reaction in this file of an endf or pendf tape.
   ! If mtd=0, first reaction is found.  Returns iauto=0 if
   ! file end is found.
   !-------------------------------------------------------------------
   use util ! provides skiprz
   use endf ! provides endf routines and variables
   ! externals
   integer::iauto,matd,mfd,mtd
   real(kr)::scr(*)
   ! internals
   integer::nin,mft,nb,nw,ir,idone
   save ir

   !--for mf3, process the reactions on pendf for this temp.
   !--exclude thermal reactions
   if (mfd.le.3) then
      nin=npend
      if (mtd.le.0) then
         mft=mfd
         call findf(matd,mft,mtd,nin)
      endif
      idone=0
      do while (idone.eq.0)
         call contio(nin,0,0,scr,nb,nw)
         if (mfh.le.0.or.mth.ne.0) then
            mtd=mth
            if (mtd.le.200) idone=1
            if (mtd.ge.203.and.mtd.le.207) idone=1
            if (mtd.gt.300) idone=1
            if (idone.eq.0) call tosend(nin,0,0,scr)
         endif
      enddo
      if (mtd.le.0) iauto=0
      call skiprz(nin,-1)
      return
   endif

   !--for other reactions, use lists from conver.
   if (mtd.eq.0) ir=0
   ir=ir+1
   ! neutrons
   if (mfd.ne.6) go to 220
   mtd=mf4(ir)
   if (mtd.ne.0.and.iabs(mtd).le.91) go to 290
   ir=1
   mfd=8
  220 continue
   if (mfd.ne.8) go to 225
   mtd=mf6(ir)
   if (mtd.ne.0) go to 290
   go to 280
   ! nuclide production (first time)
  225 continue
   if (mfd.ne.10) go to 230
   ir=1
   if (mf10f(ir).eq.0) go to 280
   mfd=mf10i(ir)
   lfs=lfs8(ir)
   isom=mlfs8(ir)
   if (mf10f(ir).eq.3) mfd=mfd+10000000
   if (mf10f(ir).eq.6) mfd=mfd+20000000
   if (mf10f(ir).eq.9) mfd=mfd+30000000
   if (mf10f(ir).eq.10) mfd=mfd+40000000
   mtd=mf10s(ir)
   go to 290
   ! photons
  230 continue
   if (mfd.eq.16) then
      mtd=mf12(ir)
      if (mtd.ne.0) return
      ir=1
      mfd=17
   endif
   if (mfd.eq.17) then
      mtd=mf13(ir)
      if (mtd.ne.0) return
      ir=1
      mfd=18
   endif
   if (mfd.eq.18) then
      mtd=mf18(ir)
      if (mtd.ne.0) return
      go to 280
   endif
   ! protons
   if (mfd.ne.21) go to 320
   mtd=mf6p(1,ir)
   if (mtd.ne.0) return
   ir=1
   mfd=31
  320 continue
   if (mfd.ne.31) go to 330
   mtd=mf4r(1,ir)
   if (mtd.ne.0) return
   go to 280
   ! deuterons
  330 continue
   if (mfd.ne.22) go to 340
   mtd=mf6p(2,ir)
   if (mtd.ne.0) return
   ir=1
   mfd=32
  340 continue
   if (mfd.ne.32) go to 350
   mtd=mf4r(2,ir)
   if (mtd.ne.0) return
   go to 280
   ! tritons
  350 continue
   if (mfd.ne.23) go to 360
   mtd=mf6p(3,ir)
   if (mtd.ne.0) return
   ir=1
   mfd=33
  360 continue
   if (mfd.ne.33) go to 370
   mtd=mf4r(3,ir)
   if (mtd.ne.0) return
   go to 280
   ! he-3 nuclei
  370 continue
   if (mfd.ne.24) go to 380
   mtd=mf6p(4,ir)
   if (mtd.ne.0) return
   ir=1
   mfd=34
  380 continue
   if (mfd.ne.34) go to 390
   mtd=mf4r(4,ir)
   if (mtd.ne.0) return
   go to 280
   ! alphas
  390 continue
   if (mfd.ne.25) go to 400
   mtd=mf6p(5,ir)
   if (mtd.ne.0) return
   ir=1
   mfd=35
  400 continue
   if (mfd.ne.35) go to 410
   mtd=mf4r(5,ir)
   if (mtd.ne.0) return
   go to 280
   ! heavy particles (a>4)
  410 continue
   if (mfd.ne.26) go to 420
   mtd=mf6p(6,ir)
   if (mtd.ne.0) return
   ir=1
   mfd=36
  420 continue
   if (mfd.ne.36) go to 430
   mtd=mf4r(6,ir)
   if (mtd.ne.0) return
   go to 280
   ! production (after first 10/ entry)
  430 continue
   if (mfd.lt.10000000) go to 280
   if (mf10f(ir).eq.0) go to 280
   mfd=mf10i(ir)
   lfs=lfs8(ir)
   isom=mlfs8(ir)
   if (mf10f(ir).eq.3) mfd=mfd+10000000
   if (mf10f(ir).eq.6) mfd=mfd+20000000
   if (mf10f(ir).eq.9) mfd=mfd+30000000
   if (mf10f(ir).eq.10) mfd=mfd+40000000
   mtd=mf10s(ir)
   go to 290

   !--return next mt, or zero when finished.
  280 continue
   mtd=0
   iauto=0
   if (mfd.eq.8) mfd=6
   if (mfd.eq.17.or.mfd.eq.18) mfd=16
   if (mfd.eq.31) mfd=21
   if (mfd.eq.32) mfd=22
   if (mfd.eq.33) mfd=23
   if (mfd.eq.34) mfd=24
   if (mfd.eq.35) mfd=25
   if (mfd.eq.36) mfd=26
 290 continue
   return
   end subroutine nextr

   subroutine namer(izad,izam,mfd,mtd,mtname)
   !-------------------------------------------------------------------
   ! Generates standard names for ENDF reactions in 4 char. words.
   ! The letter 'h' is used for He-3.
   !-------------------------------------------------------------------
   use endf ! provides iverf
   ! externals
   integer::izad,izam,mfd,mtd
   character(4)::mtname(15)
   ! internals
   integer::i,j,izaa,imm
   character(60)::dummy,result
   character(1)::proj
   character(7)::reac
   character(8)::part
   character(8)::azam
   integer,dimension(7),parameter::ip=(/&
     0,1,1001,1002,1003,2003,2004/)
   character(1),dimension(7),parameter::np=(/&
     'g','n','p','d','t','h','a'/)
   integer,dimension(17),parameter::i2=(/&
     16,17,18,6,8,21,22,23,24,25,26,31,32,33,34,35,36/)
   character(8),dimension(17),parameter::n2=(/&
     'gamma   ','gamma   ','gamma   ','neutron ','neutron ',&
     'proton  ','deuteron','triton  ','he3     ','alpha   ',&
     'recoil  ','proton  ','deuteron','triton  ','he3     ',&
     'alpha   ','recoil  '/)
   integer,dimension(113),parameter::ir=(/&
     11,16,17,22,23,24,25,28,29,30,32,33,34,35,36,37,41,42,44,45,&
     108,109,111,112,113,114,115,116,117,&
         152,153,154,155,156,157,158,159,160,&
     161,162,163,164,165,166,167,168,169,170,&
     171,172,173,174,175,176,177,178,179,180,&
     181,182,183,184,185,186,187,188,189,190,&
     191,192,193,194,195,196,197,198,199,200,&
     18,19,20,21,38,&
     102,103,104,105,106,107,201,202,203,204,205,206,207,&
     501,502,504,515,516,517,1,2,3,4,5,10,101,301,443,444,599/)
   character(7),dimension(113),parameter::nr=(/&
     '2nd    ','2n     ','3n     ','na     ','n3a    ','2na    ',&
     '3na    ','np     ','n2a    ','2n2a   ','nd     ','nt     ',&
     'nh     ','nd2a   ','nt2a   ','4n     ','2np    ','3np    ',&
     'n2p    ','npa    ','2a     ','3a     ','2p     ','pa     ',&
     't2a    ','d2a    ','pd     ','pt     ','da     ',&
     '5n     ','6n     ','2nt    ','ta     ','4np    ','3nd    ',&
     'nda    ','2npa   ','7n     ','8n     ','5np    ','6np    ',&
     '7np    ','4na    ','5na    ','6na    ','7na    ','4nd    ',&
     '5nd    ','6nd    ','3nt    ','4nt    ','5nt    ','6nt    ',&
     '2nh    ','3nh    ','4nh    ','3n2p   ','3n2a   ','3npa   ',&
     'dt     ','npd    ','npt    ','ndt    ','nph    ','ndh    ',&
     'nth    ','nta    ','2n2p   ','ph     ','dh     ','ha     ',&
     '4n2p   ','4n2a   ','4npa   ','3p     ','n3p    ','3n2pa  ',&
     '5n2p   ','fission',&
     'f      ','nf     ','2nf    ','3nf    ','g      ','p      ',&
     'd      ','t      ','h      ','a      ','xn     ','xg     ',&
     'xp     ','xd     ','xt     ','xh     ','xa     ','total  ',&
     'coh    ','incoh  ','paire  ','pair   ','pairn  ','total  ',&
     'elastic','nonel  ','inel   ','x      ','cont   ','dis    ',&
     'heat   ','kerma  ','damage ','stop   '/)
   integer,parameter::nreac=113
   integer,parameter::npart=17
   integer,parameter::nproj=7

   !--get projectile name
   proj='?'
   do i=1,nproj
      if (izad.eq.ip(i)) proj=np(i)
   enddo

   !--*get reaction name
   reac='?'
   if (iverf.lt.6) then
      if (mtd.ge.51.and.mtd.lt.60) then
         write(reac,'(''n0'',i1)') mtd-50
      else if (mtd.ge.60.and.mtd.le.90) then
         write(reac,'(''n'',i2)') mtd-50
      else if (mtd.eq.91) then
         reac='nc'
      else if (mtd.ge.700.and.mtd.le.709) then
         write(reac,'(''p0'',i1)') mtd-700
      else if (mtd.ge.710.and.mtd.le.718) then
         write(reac,'(''p'',i2)') mtd-700
      else if (mtd.eq.719) then
         reac='pc'
      else if (mtd.ge.720.and.mtd.le.729) then
         write(reac,'(''d0'',i1)') mtd-720
      else if (mtd.ge.730.and.mtd.le.738) then
         write(reac,'(''d'',i2)') mtd-720
      else if (mtd.eq.739) then
         reac='dc'
      else if (mtd.ge.740.and.mtd.le.749) then
         write(reac,'(''t0'',i1)') mtd-740
      else if (mtd.ge.750.and.mtd.le.758) then
         write(reac,'(''t'',i2)') mtd-740
      else if (mtd.eq.759) then
         reac='tc'
      else if (mtd.ge.760.and.mtd.le.769) then
         write(reac,'(''h0'',i1)') mtd-760
      else if (mtd.ge.770.and.mtd.le.778) then
         write(reac,'(''h'',i2)') mtd-760
      else if (mtd.eq.779) then
         reac='hc'
      else if (mtd.ge.780.and.mtd.le.789) then
         write(reac,'(''a0'',i1)') mtd-780
      else if (mtd.ge.790.and.mtd.le.798) then
         write(reac,'(''a'',i2)') mtd-780
      else if (mtd.eq.799) then
         reac='ac'
      else if (mtd.eq.602) then
         reac='photo'
      endif
   else
      if (mtd.ge.50.and.mtd.le.59) then
         write(reac,'(''n0'',i1)') mtd-50
      else if (mtd.ge.60.and.mtd.le.90) then
         write(reac,'(''n'',i2)') mtd-50
      else if (mtd.eq.91) then
         reac='nc'
      else if (mtd.ge.600.and.mtd.le.609) then
         write(reac,'(''p0'',i1)') mtd-600
      else if (mtd.ge.610.and.mtd.le.648) then
         write(reac,'(''p'',i2)') mtd-600
      else if (mtd.eq.649) then
         reac='pc'
      else if (mtd.ge.650.and.mtd.le.659) then
         write(reac,'(''d0'',i1)') mtd-650
      else if (mtd.ge.660.and.mtd.le.698) then
         write(reac,'(''d'',i2)') mtd-650
      else if (mtd.eq.699) then
         reac='dc'
      else if (mtd.ge.700.and.mtd.le.709) then
         write(reac,'(''t0'',i1)') mtd-700
      else if (mtd.ge.710.and.mtd.le.748) then
         write(reac,'(''t'',i2)') mtd-700
      else if (mtd.eq.749) then
         reac='tc'
      else if (mtd.ge.750.and.mtd.le.759) then
         write(reac,'(''h0'',i1)') mtd-750
      else if (mtd.ge.760.and.mtd.le.798) then
         write(reac,'(''h'',i2)') mtd-750
      else if (mtd.eq.799) then
         reac='hc'
      else if (mtd.ge.800.and.mtd.le.809) then
         write(reac,'(''a0'',i1)') mtd-800
      else if (mtd.ge.810.and.mtd.le.848) then
         write(reac,'(''a'',i2)') mtd-800
      else if (mtd.eq.849) then
         reac='ac'
      else if (mtd.ge.875.and.mtd.le.884) then
         write(reac,'(''2n0'',i1)') mtd-875
      else if (mtd.ge.885.and.mtd.le.890) then
         write(reac,'(''2n'',i2)') mtd-875
      else if (mtd.eq.891) then
         reac='2nc'
      else if (mtd.eq.522) then
         reac='photo'
      else if (mtd.eq.500) then
         reac='stop'
      endif
   endif
   if (reac.eq.'?') then
      do i=1,nreac
         if (mtd.eq.ir(i)) reac=nr(i)
      enddo
   endif

   !--get outgoing particle name
   if (mfd.eq.3) then
      dummy='('//proj//','//reac//')-cross-section.'
   else if (mfd.ge.10000000) then
      izaa=izam/10
      imm=isom
      if (imm.eq.0) then
         write(azam,'(i5)') izaa
         dummy='('//proj//','//reac//')-'//azam(1:5)//'-production.'
      else if (imm.lt.10) then
         write(azam,'(i5,''m'',i1)') izaa,imm
         dummy='('//proj//','//reac//')-'//azam//'-production.'
      else
         izaa=izaa/10
         write(azam,'(i5,''m'',i2)') izaa,imm
         dummy='('//proj//','//reac//')-'//azam//'-production.'
      endif
   else
      part='?'
      do i=1,npart
         if (mfd.eq.i2(i)) part=n2(i)
      enddo
      dummy='('//proj//','//reac//')-'//part//'-matrix.'
   endif

   !--strip unneeded spaces
   i=1
   j=0
   result=' '
   do while (dummy(i:i).ne.'.')
      do while (dummy(i:i).eq.' ')
         i=i+1
      enddo
      j=j+1
      if (dummy(i:i).eq.'-') dummy(i:i)=' '
      result(j:j)=dummy(i:i)
      i=i+1
   enddo

   !--return final result in hollerith form
   read(result,'(15a4)') mtname
   return
   end subroutine namer

   subroutine mfchk(mf,mt)
   !-------------------------------------------------------------------
   ! Check for 12 required when 13 is given, or
   ! for 17 required when 16 is given.
   !-------------------------------------------------------------------
   ! externals
   integer::mf,mt
   ! internals
   integer::max,i

   max=maxr1
   if (mf.eq.13) then
      i=1
      do while (i.le.max)
         if (mf12(i).eq.0) return
         if (mf12(i).eq.mt) then
            mf=12
            return
         endif
         if (i.gt.1) then
            if (mf12(i).le.0.and.&
              mf12(i-1).le.mt.and.(-mf12(i)).ge.mt) then
               mf=12
               return
            endif
         endif
         i=i+1
      enddo
   else if (mf.eq.16) then
      i=1
      do while (i.le.max)
         if (mf13(i).eq.0) return
         if (mf13(i).eq.mt) then
            mf=17
            return
         endif
         if (i.gt.1) then
            if (mf13(i).le.0.and.&
              mf13(i-1).le.mt.and.(-mf13(i)).ge.mt) then
               mf=17
               return
            endif
         endif
         i=i+1
      enddo
   endif
   return
   end subroutine mfchk

   subroutine mfchk2(mf,mt)
   !-------------------------------------------------------------------
   ! Check for charged particle production from mf4
   ! recoils when mf=21-26 is requested.
   !-------------------------------------------------------------------
   ! externals
   integer::mf,mt
   ! internals
   integer::max,i,mtd

   max=maxr1
   i=1
   do while (i.le.max)
      mtd=mf4r(mf-20,i)
      if (mtd.eq.0) return
      if (mt.eq.mtd) then
         mf=mf+10
         return
      endif
      i=i+1
   enddo
   return
   end subroutine mfchk2

   subroutine gengpn
   !-------------------------------------------------------------------
   ! Generate requested neutron group structure or read in from
   ! the system input file in the form of an ENDF list record
   !
   !    ign     meaning
   !    ---     ---------------------------------------
   !     1      arbitrary structure (read in)
   !     2      CSEWG 239 group structure
   !     3      LANL 30 group structure
   !     4      ANL 27 group structure
   !     5      RRD 50 group structure
   !     6      GAM-I 68 group structure
   !     7      GAM-II 100 group structure
   !     8      LASER-THERMOS 35 group
   !     9      EPRI-CPM/WIMS 69 group structure
   !    10      LANL 187-group structure
   !    11      LANL 70-group structure
   !    12      SAND-II 620-group structure
   !    13      LANL 80-group structure
   !    14      EURLIB 100-group structure
   !    15      SAND-IIA 640-group structure
   !    16      VITAMIN-E 174-group structure
   !    17      VITAMIN-J 175-group structure
   !    18      XMAS 172-group structure
   !    19      ECCO  33-group structure
   !    20      ECCO 1968-group structure
   !    21      TRIPOLI 315-group structure
   !    22      XMASLWC 172-group structure
   !    23      VIT-J lwpc 175-group structure
   !    24      SHEM CEA 281-group structure
   !    25      SHEM EPM 295-group structure
   !    26      SHEM CEA/EPM 361-group structure
   !    27      SHEM EPM 315-group structure
   !    28      RAHAB AECL 89-group structure
   !    29      CCFE   660-group structure
   !    30      UKAEA 1025-group structure
   !    31      UKAEA 1067-group structure
   !    32      UKAEA 1102-group structure
   !    33      UKAEA  142-group structure
   !    34      LANL 618 group structure
   !
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use util   ! provides error
   ! internals
   integer::lflag,ig,ngp,n1,n2,n,i,ic
   integer::ngp_hold
   real(kr)::u,du,delta
   real(kr),dimension(241),parameter::gl2=(/&
     27.631e0_kr,17.0e0_kr,16.75e0_kr,16.588e0_kr,16.5e0_kr,16.3e0_kr,&
     16.25e0_kr,16.0e0_kr,15.75e0_kr,15.5e0_kr,15.25e0_kr,15.e0_kr,&
     14.75e0_kr,14.5e0_kr,14.25e0_kr,14.e0_kr,13.75e0_kr,13.5e0_kr,&
     13.25e0_kr,13.e0_kr,12.75e0_kr,12.5e0_kr,12.25e0_kr,12.e0_kr,&
     11.75e0_kr,11.5e0_kr,11.25e0_kr,11.e0_kr,10.75e0_kr,10.5e0_kr,&
     10.25e0_kr,10.e0_kr,9.75e0_kr,9.5e0_kr,9.25e0_kr,9.e0_kr,&
     8.9e0_kr,8.8e0_kr,8.75e0_kr,8.7e0_kr,8.6e0_kr,8.5e0_kr,8.4e0_kr,&
     8.3e0_kr,8.25e0_kr,8.2e0_kr,8.1583e0_kr,8.1e0_kr,8.e0_kr,&
     7.9e0_kr,7.8e0_kr,7.75e0_kr,7.7e0_kr,7.6e0_kr,7.5e0_kr,&
     7.375e0_kr,7.25e0_kr,7.125e0_kr,7.e0_kr,6.875e0_kr,6.75e0_kr,&
     6.625e0_kr,6.5e0_kr,6.375e0_kr,6.25e0_kr,6.15e0_kr,6.125e0_kr,&
     6.05e0_kr,6.025e0_kr,6.e0_kr,5.95e0_kr,5.875e0_kr,5.75e0_kr,&
     5.675e0_kr,5.65e0_kr,5.625e0_kr,5.5e0_kr,5.375e0_kr,5.25e0_kr,&
     5.175e0_kr,5.125e0_kr,5.075e0_kr,5.e0_kr,4.875e0_kr,4.75e0_kr,&
     4.625e0_kr,4.5e0_kr,4.45e0_kr,4.4e0_kr,4.35e0_kr,4.3e0_kr,&
     4.25e0_kr,4.2e0_kr,4.15e0_kr,4.125e0_kr,4.1e0_kr,4.075e0_kr,&
     4.05e0_kr,4.e0_kr,3.95e0_kr,3.9e0_kr,3.85e0_kr,3.8e0_kr,&
     3.75e0_kr,3.7e0_kr,3.65e0_kr,3.6e0_kr,3.575e0_kr,3.55e0_kr,&
     3.525e0_kr,3.5e0_kr,3.475e0_kr,3.45e0_kr,3.4e0_kr,3.35e0_kr,&
     3.3e0_kr,3.25e0_kr,3.2e0_kr,3.15e0_kr,3.1e0_kr,3.05e0_kr,&
     3.e0_kr,2.975e0_kr,2.95e0_kr,2.925e0_kr,2.9e0_kr,2.85e0_kr,&
     2.8e0_kr,2.75e0_kr,2.7e0_kr,2.65e0_kr,2.6e0_kr,2.55e0_kr,&
     2.5e0_kr,2.45e0_kr,2.4e0_kr,2.35e0_kr,2.3417e0_kr,2.325e0_kr,&
     2.3e0_kr,2.25e0_kr,2.2e0_kr,2.15e0_kr,2.125e0_kr,2.1e0_kr,&
     2.05e0_kr,2.e0_kr,1.95e0_kr,1.9e0_kr,1.875e0_kr,1.85e0_kr,&
     1.825e0_kr,1.8e0_kr,1.75e0_kr,1.7e0_kr,1.675e0_kr,1.65e0_kr,&
     1.625e0_kr,1.6e0_kr,1.55e0_kr,1.5e0_kr,1.4833e0_kr,1.4667e0_kr,&
     1.45e0_kr,1.4417e0_kr,1.4333e0_kr,1.4167e0_kr,1.4e0_kr,1.35e0_kr,&
     1.3e0_kr,1.25e0_kr,1.2e0_kr,1.175e0_kr,1.15e0_kr,1.125e0_kr,&
     1.1e0_kr,1.05e0_kr,1.e0_kr,.95e0_kr,.9e0_kr,.85e0_kr,.8e0_kr,&
     .775e0_kr,.75e0_kr,.725e0_kr,.7e0_kr,.65e0_kr,.6e0_kr,.55e0_kr,&
     .525e0_kr,.5e0_kr,.475e0_kr,.45e0_kr,.425e0_kr,.41667e0_kr,&
     .40833e0_kr,.4e0_kr,.375e0_kr,.35e0_kr,.325e0_kr,.3e0_kr,&
     .275e0_kr,.25e0_kr,.225e0_kr,.2e0_kr,.175e0_kr,.15e0_kr,&
     .125e0_kr,.1e0_kr,.075e0_kr,.05e0_kr,.025e0_kr,0.e0_kr,&
     -.025e0_kr,-.05e0_kr,-.075e0_kr,-.1e0_kr,-.125e0_kr,-.15e0_kr,&
     -.175e0_kr,-.2e0_kr,-.225e0_kr,-.25e0_kr,-.275e0_kr,-.3e0_kr,&
     -.325e0_kr,-.35e0_kr,-.375e0_kr,-.4e0_kr,-.425e0_kr,-.45e0_kr,&
     -.475e0_kr,-.5e0_kr,-.525e0_kr,-.55e0_kr,-.575e0_kr,-.6e0_kr,&
     -.625e0_kr,-.65e0_kr,-.675e0_kr,-.69167e0_kr/)
   real(kr),dimension(31),parameter::eg3=(/&
     1.39e-4_kr,1.52e-1_kr,4.14e-1_kr,1.13e0_kr,3.06e0_kr,8.32e0_kr,&
     2.26e1_kr,6.14e1_kr,1.67e2_kr,4.54e2_kr,1.235e3_kr,3.35e3_kr,&
     9.12e3_kr,2.48e4_kr,6.76e4_kr,1.84e5_kr,3.03e5_kr,5.00e5_kr,&
     8.23e5_kr,1.353e6_kr,1.738e6_kr,2.232e6_kr,2.865e6_kr,3.68e6_kr,&
     6.07e6_kr,7.79e6_kr,1.00e7_kr,1.20e7_kr,1.35e7_kr,1.50e7_kr,&
     1.70e7_kr/)
   real(kr),dimension(28),parameter::gl4=(/&
      14.5e0_kr,13.0e0_kr,12.5e0_kr,12.0e0_kr,11.5e0_kr,11.0e0_kr,&
      10.5e0_kr,10.0e0_kr,9.5e0_kr,9.0e0_kr,8.5e0_kr,8.0e0_kr,&
      7.5e0_kr,7.0e0_kr,6.5e0_kr,6.0e0_kr,5.5e0_kr,5.0e0_kr,4.5e0_kr,&
      4.0e0_kr,3.5e0_kr,3.0e0_kr,2.5e0_kr,2.0e0_kr,1.5e0_kr,1.0e0_kr,&
      0.5e0_kr,0.0e0_kr/)
   real(kr),dimension(51),parameter::gl5=(/&
      27.631e0_kr,16.5e0_kr,16.e0_kr,15.5e0_kr,15.e0_kr,14.5e0_kr,&
      14.e0_kr,13.5e0_kr,13.e0_kr,12.5e0_kr,12.e0_kr,11.5e0_kr,&
      11.e0_kr,10.5e0_kr,10.25e0_kr,10.e0_kr,9.75e0_kr,9.5e0_kr,&
      9.25e0_kr,9.e0_kr,8.75e0_kr,8.5e0_kr,8.25e0_kr,8.e0_kr,&
      7.75e0_kr,7.5e0_kr,7.25e0_kr,7.e0_kr,6.75e0_kr,6.5e0_kr,&
      6.25e0_kr,6.e0_kr,5.75e0_kr,5.5e0_kr,5.25e0_kr,5.e0_kr,&
      4.75e0_kr,4.5e0_kr,4.25e0_kr,4.e0_kr,3.75e0_kr,3.5e0_kr,&
      3.25e0_kr,3.e0_kr,2.5e0_kr,2.e0_kr,1.5e0_kr,1.e0_kr,.5e0_kr,&
      0.e0_kr,-.6917e0_kr/)
   real(kr),dimension(36),parameter::eg6=(/&
     .253e-3_kr,.2277e-2_kr,.6325e-2_kr,.12397e-1_kr,.20493e-1_kr,&
     .30613e-1_kr,.42757e-1_kr,.56925e-1_kr,.81972e-1_kr,.11159e0_kr,&
     .14573e0_kr,.18444e0_kr,.2277e0_kr,.25104e0_kr,.27053e0_kr,&
     .29075e0_kr,.30113e0_kr,.32064e0_kr,.35768e0_kr,.41704e0_kr,&
     .50326e0_kr,.62493e0_kr,.78211e0_kr,.95070e0_kr,.10137e+1_kr,&
     .10428e+1_kr,.10525e+1_kr,.10624e+1_kr,.10722e+1_kr,&
     .10987e+1_kr,.11664e+1_kr,.13079e+1_kr,.14575e+1_kr,.1595e+1_kr,&
     .17262e+1_kr,.1855e+1_kr/)
   real(kr),dimension(70),parameter::eg9=(/&
     1.e-5_kr,.005e0_kr,.01e0_kr,.015e0_kr,.02e0_kr,.025e0_kr,&
     .03e0_kr,.035e0_kr,.042e0_kr,.05e0_kr,.058e0_kr,.067e0_kr,&
     .08e0_kr,.1e0_kr,.14e0_kr,.18e0_kr,.22e0_kr,.25e0_kr,&
     .28e0_kr,.3e0_kr,.32e0_kr,.35e0_kr,.4e0_kr,.5e0_kr,.625e0_kr,&
     .78e0_kr,.85e0_kr,.91e0_kr,.95e0_kr,.972e0_kr,.996e0_kr,&
     1.02e0_kr,1.045e0_kr,1.071e0_kr,1.097e0_kr,1.123e0_kr,1.15e0_kr,&
     1.3e0_kr,1.5e0_kr,2.1e0_kr,2.6e0_kr,3.3e0_kr,4.e0_kr,9.877e0_kr,&
     15.968e0_kr,27.7e0_kr,48.052e0_kr,75.501e0_kr,148.728e0_kr,&
     367.262e0_kr,906.898e0_kr,1425.1e0_kr,2239.45e0_kr,3519.1e0_kr,&
     5530.e0_kr,9118.e0_kr,1.503e4_kr,2.478e4_kr,4.085e4_kr,&
     6.734e4_kr,1.11e5_kr,1.83e5_kr,3.025e5_kr,5.e5_kr,8.21e5_kr,&
     1.353e6_kr,2.231e6_kr,3.679e6_kr,6.0655e6_kr,1.e7_kr/)
   real(kr),dimension(48),parameter::eg10a=(/&
     1.e-5_kr,2.5399e-4_kr,7.6022e-4_kr,2.2769e-3_kr,6.3247e-3_kr,&
     .012396e0_kr,.020492e0_kr,.0255e0_kr,.030612e0_kr,.0355e0_kr,&
     .042755e0_kr,.05e0_kr,.056922e0_kr,.067e0_kr,.081968e0_kr,&
     .11157e0_kr,.14572e0_kr,.1523e0_kr,.18443e0_kr,.22769e0_kr,&
     .25103e0_kr,.27052e0_kr,.29074e0_kr,.30112e0_kr,.32063e0_kr,&
     .35767e0_kr,.41499e0_kr,.50323e0_kr,.62506e0_kr,.78208e0_kr,&
     .83368e0_kr,.87642e0_kr,.91e0_kr,.95065e0_kr,.971e0_kr,&
     .992e0_kr,1.0137e0_kr,1.0427e0_kr,1.0525e0_kr,1.0623e0_kr,&
     1.0722e0_kr,1.0987e0_kr,1.1254e0_kr,1.1664e0_kr,1.3079e0_kr,&
     1.4574e0_kr,1.5949e0_kr,1.7261e0_kr/)
   real(kr),dimension(13),parameter::eg10b=(/&
     1.1e7_kr,1.2e7_kr,1.3e7_kr,1.35e7_kr,1.375e7_kr,1.394e7_kr,&
     1.42e7_kr,1.442e7_kr,1.464e7_kr,1.5e7_kr,1.6e7_kr,1.7e7_kr,&
     2.e7_kr/)
   real(kr),dimension(71),parameter::eg11=(/10.677e0_kr,61.4421e0_kr,&
     101.301e0_kr,130.073e0_kr,167.017e0_kr,214.454e0_kr,275.365e0_kr,&
     353.575e0_kr,453.999e0_kr,582.947e0_kr,748.518e0_kr,961.117e0_kr,&
     1089.09e0_kr,1234.1e0_kr,1398.42e0_kr,1584.61e0_kr,1795.6e0_kr,&
     2034.68e0_kr,2305.6e0_kr,2612.59e0_kr,2960.45e0_kr,3354.63e0_kr,&
     3801.29e0_kr,4307.43e0_kr,4880.95e0_kr,5530.84e0_kr,6267.27e0_kr,&
     7101.74e0_kr,8047.33e0_kr,9118.82e0_kr,10333.3e0_kr,11708.8e0_kr,&
     13267.8e0_kr,15034.4e0_kr,17036.2e0_kr,19304.5e0_kr,21874.9e0_kr,&
     24787.5e0_kr,28087.9e0_kr,31827.8e0_kr,40867.7e0_kr,52475.2e0_kr,&
     67379.5e0_kr,86517.e0_kr,111090.e0_kr,142642.e0_kr,183156.e0_kr,&
     235178.e0_kr,301974.e0_kr,387742.e0_kr,439369.e0_kr,497871.e0_kr,&
     564161.e0_kr,639279.e0_kr,724398.e0_kr,820850.e0_kr,930145.e0_kr,&
     1053990.e0_kr,1194330.e0_kr,1353350.e0_kr,1737740.e0_kr,&
     2231300.e0_kr,2865050.e0_kr,3678790.e0_kr,4723670.e0_kr,&
     6065310.e0_kr,7788010.e0_kr,1.e7_kr,1.28403e7_kr,1.64872e7_kr,&
     2.e7_kr/)
   real(kr),dimension(8),parameter::deltl=(/&
     5.e0_kr,7.5e0_kr,10.e0_kr,15.e0_kr,20.e0_kr,25.e0_kr,30.e0_kr,&
     40.e0_kr/)
   integer,dimension(9),parameter::ndelta=(/2,6,10,19,23,28,36,40,46/)
   real(kr),dimension(80),parameter::u80=(/&
     .1681472e0_kr,.125e0_kr,.1e0_kr,.125e0_kr,.175e0_kr,.25e0_kr,&
     .25e0_kr,.25e0_kr,.25e0_kr,.25e0_kr,.25e0_kr,.25e0_kr,.25e0_kr,&
     .125e0_kr,.125e0_kr,.125e0_kr,.125e0_kr,.125e0_kr,.125e0_kr,&
     .125e0_kr,.125e0_kr,.125e0_kr,.125e0_kr,.25e0_kr,.25e0_kr,&
     .25e0_kr,.25e0_kr,.25e0_kr,.25e0_kr,.25e0_kr,.25e0_kr,.25e0_kr,&
     .25e0_kr,.125e0_kr,.075e0_kr,.05e0_kr,.125e0_kr,.125e0_kr,&
     .125e0_kr,.125e0_kr,.125e0_kr,.125e0_kr,.125e0_kr,.125e0_kr,&
     .125e0_kr,.125e0_kr,.125e0_kr,.125e0_kr,.125e0_kr,.125e0_kr,&
     .125e0_kr,.125e0_kr,.125e0_kr,.125e0_kr,.125e0_kr,.125e0_kr,&
     .125e0_kr,.125e0_kr,.125e0_kr,.125e0_kr,.125e0_kr,.125e0_kr,&
     .25e0_kr,.25e0_kr,.25e0_kr,.25e0_kr,.25e0_kr,.5e0_kr,.5e0_kr,&
     .5e0_kr,.5e0_kr,.5e0_kr,.5e0_kr,.5e0_kr,.5e0_kr,.5e0_kr,1.e0_kr,&
     1.e0_kr,1.e0_kr,7.e0_kr/)
   integer,dimension(19),parameter::ig14=(/&
     2,9,13,15,17,23,25,55,60,61,63,64,65,93,94,95,99,100,101/)
   real(kr),dimension(19),parameter::gl14=(/&
     .1e0_kr,.05e0_kr,.1e0_kr,.05e0_kr,.1e0_kr,.05e0_kr,.1e0_kr,&
     .25e0_kr,.2e0_kr,.05e0_kr,.075e0_kr,.125e0_kr,.25e0_kr,&
     .5e0_kr,.25e0_kr,.5e0_kr,.588e0_kr,.412e0_kr,10.631e0_kr/)
   real(kr),dimension(84),parameter::eg15a=(/&
     1.0e-5_kr,1.0e-1_kr,4.1399e-1_kr,5.3158e-1_kr,6.8256e-1_kr,&
     8.7642e-1_kr,1.1254e0_kr,1.4450e0_kr,1.8554e0_kr,2.3824e0_kr,&
     3.0590e0_kr,3.9279e0_kr,5.0435e0_kr,6.4760e0_kr,8.3153e0_kr,&
     1.0677e1_kr,1.3710e1_kr,1.7603e1_kr,2.2603e1_kr,2.9023e1_kr,&
     3.7267e1_kr,4.7851e1_kr,6.1442e1_kr,7.8893e1_kr,1.0130e2_kr,&
     1.3007e2_kr,1.6702e2_kr,2.1445e2_kr,2.7536e2_kr,3.5358e2_kr,&
     4.5400e2_kr,5.8295e2_kr,7.4852e2_kr,9.6112e2_kr,1.2341e3_kr,&
     1.5846e3_kr,2.0347e3_kr,2.2487e3_kr,2.4852e3_kr,2.6126e3_kr,&
     2.7465e3_kr,3.0354e3_kr,3.3546e3_kr,3.7074e3_kr,4.3074e3_kr,&
     5.5308e3_kr,7.1017e3_kr,9.1188e3_kr,1.0595e4_kr,1.1709e4_kr,&
     1.5034e4_kr,1.9305e4_kr,2.1875e4_kr,2.3579e4_kr,2.4176e4_kr,&
     2.4788e4_kr,2.6058e4_kr,2.7000e4_kr,2.8500e4_kr,3.1828e4_kr,&
     3.4307e4_kr,4.0868e4_kr,4.6309e4_kr,5.2475e4_kr,5.6562e4_kr,&
     6.7379e4_kr,7.2000e4_kr,7.9500e4_kr,8.2500e4_kr,8.6517e4_kr,&
     9.8037e4_kr,1.1109e5_kr,1.1679e5_kr,1.2277e5_kr,1.2907e5_kr,&
     1.3569e5_kr,1.4264e5_kr,1.4996e5_kr,1.5764e5_kr,1.6573e5_kr,&
     1.7422e5_kr,1.8316e5_kr,1.9255e5_kr,2.0242e5_kr/)
   real(kr),dimension(91),parameter::eg15b=(/&
     2.1280e5_kr,2.2371e5_kr,2.3518e5_kr,2.4724e5_kr,2.7324e5_kr,&
     2.8725e5_kr,2.9452e5_kr,2.9720e5_kr,2.9850e5_kr,3.0197e5_kr,&
     3.3373e5_kr,3.6883e5_kr,3.8774e5_kr,4.0762e5_kr,4.5049e5_kr,&
     4.9787e5_kr,5.2340e5_kr,5.5023e5_kr,5.7844e5_kr,6.0810e5_kr,&
     6.3928e5_kr,6.7206e5_kr,7.0651e5_kr,7.4274e5_kr,7.8082e5_kr,&
     8.2085e5_kr,8.6294e5_kr,9.0718e5_kr,9.6164e5_kr,1.0026e6_kr,&
     1.1080e6_kr,1.1648e6_kr,1.2246e6_kr,1.2873e6_kr,1.3534e6_kr,&
     1.4227e6_kr,1.4957e6_kr,1.5724e6_kr,1.6530e6_kr,1.7377e6_kr,&
     1.8268e6_kr,1.9205e6_kr,2.0190e6_kr,2.1225e6_kr,2.2313e6_kr,&
     2.3069e6_kr,2.3457e6_kr,2.3653e6_kr,2.3852e6_kr,2.4660e6_kr,&
     2.5924e6_kr,2.7253e6_kr,2.8650e6_kr,3.0119e6_kr,3.1664e6_kr,&
     3.3287e6_kr,3.6788e6_kr,4.0657e6_kr,4.4933e6_kr,4.7237e6_kr,&
     4.9659e6_kr,5.2205e6_kr,5.4881e6_kr,5.7695e6_kr,6.0653e6_kr,&
     6.3763e6_kr,6.5924e6_kr,6.7032e6_kr,7.0469e6_kr,7.4082e6_kr,&
     7.7880e6_kr,8.1873e6_kr,8.6071e6_kr,9.0484e6_kr,9.5123e6_kr,&
     1.0000e7_kr,1.0513e7_kr,1.1052e7_kr,1.1618e7_kr,1.2214e7_kr,&
     1.2523e7_kr,1.3499e7_kr,1.3840e7_kr,1.4191e7_kr,1.4550e7_kr,&
     1.4918e7_kr,1.5683e7_kr,1.6487e7_kr,1.6905e7_kr,1.7333e7_kr,&
     1.9640e7_kr/)
   real(kr),dimension(173),parameter::eg18=(/&
     1.96403e7_kr,1.73325e7_kr,1.49182e7_kr,1.38403e7_kr,1.16183e7_kr,&
     1.00000e7_kr,8.18731e6_kr,6.70320e6_kr,6.06531e6_kr,5.48812e6_kr,&
     4.49329e6_kr,3.67879e6_kr,3.01194e6_kr,2.46597e6_kr,2.23130e6_kr,&
     2.01897e6_kr,1.65299e6_kr,1.35335e6_kr,1.22456e6_kr,1.10803e6_kr,&
     1.00259e6_kr,9.07180e5_kr,8.20850e5_kr,6.08101e5_kr,5.50232e5_kr,&
     4.97871e5_kr,4.50492e5_kr,4.07622e5_kr,3.01974e5_kr,2.73237e5_kr,&
     2.47235e5_kr,1.83156e5_kr,1.22773e5_kr,1.11090e5_kr,8.22975e4_kr,&
     6.73795e4_kr,5.51656e4_kr,4.08677e4_kr,3.69786e4_kr,2.92830e4_kr,&
     2.73944e4_kr,2.47875e4_kr,1.66156e4_kr,1.50344e4_kr,1.11378e4_kr,&
     9.11882e3_kr,7.46586e3_kr,5.53084e3_kr,5.00451e3_kr,3.52662e3_kr,&
     3.35463e3_kr,2.24867e3_kr,2.03468e3_kr,1.50733e3_kr,1.43382e3_kr,&
     1.23410e3_kr,1.01039e3_kr,9.14242e2_kr,7.48518e2_kr,6.77287e2_kr,&
     4.53999e2_kr,3.71703e2_kr,3.04325e2_kr,2.03995e2_kr,1.48625e2_kr,&
     1.36742e2_kr,9.16609e1_kr,7.56736e1_kr,6.79041e1_kr,5.55951e1_kr,&
     5.15780e1_kr,4.82516e1_kr,4.55174e1_kr,4.01690e1_kr,3.72665e1_kr,&
     3.37201e1_kr,3.05113e1_kr,2.76077e1_kr,2.49805e1_kr,2.26033e1_kr,&
     1.94548e1_kr,1.59283e1_kr,1.37096e1_kr,1.12245e1_kr,9.90555e0_kr,&
     9.18981e0_kr,8.31529e0_kr,7.52398e0_kr,6.16012e0_kr,5.34643e0_kr,&
     5.04348e0_kr,4.12925e0_kr,4.00000e0_kr,3.38075e0_kr,3.30000e0_kr,&
     2.76792e0_kr,2.72000e0_kr,2.60000e0_kr,2.55000e0_kr,2.36000e0_kr,&
     2.13000e0_kr,2.10000e0_kr,2.02000e0_kr,1.93000e0_kr,1.84000e0_kr,&
     1.75500e0_kr,1.67000e0_kr,1.59000e0_kr,1.50000e0_kr,1.47500e0_kr,&
     1.44498e0_kr,1.37000e0_kr,1.33750e0_kr,1.30000e0_kr,1.23500e0_kr,&
     1.17000e0_kr,1.15000e0_kr,1.12535e0_kr,1.11000e0_kr,1.09700e0_kr,&
     1.07100e0_kr,1.04500e0_kr,1.03500e0_kr,1.02000e0_kr,&
     9.96000e-1_kr,9.86000e-1_kr,9.72000e-1_kr,9.50000e-1_kr,&
     9.30000e-1_kr,9.10000e-1_kr,8.60000e-1_kr,8.50000e-1_kr,&
     7.90000e-1_kr,7.80000e-1_kr,7.05000e-1_kr,6.25000e-1_kr,&
     5.40000e-1_kr,5.00000e-1_kr,4.85000e-1_kr,4.33000e-1_kr,&
     4.00000e-1_kr,3.91000e-1_kr,3.50000e-1_kr,3.20000e-1_kr,&
     3.14500e-1_kr,3.00000e-1_kr,2.80000e-1_kr,2.48000e-1_kr,&
     2.20000e-1_kr,1.89000e-1_kr,1.80000e-1_kr,1.60000e-1_kr,&
     1.40000e-1_kr,1.34000e-1_kr,1.15000e-1_kr,1.00001e-1_kr,&
     9.50000e-2_kr,8.00000e-2_kr,7.70000e-2_kr,6.70000e-2_kr,&
     5.80000e-2_kr,5.00000e-2_kr,4.20000e-2_kr,3.50000e-2_kr,&
     3.00000e-2_kr,2.50000e-2_kr,2.00000e-2_kr,1.50000e-2_kr,&
     1.00000e-2_kr,6.90000e-3_kr,5.00000e-3_kr,3.00000e-3_kr,&
     1.00001e-5_kr/)
   real(kr),dimension(34),parameter::eg19=(/&
     1.000010e-05_kr,1.000000e-01_kr,5.400000e-01_kr,4.000000e+00_kr,&
     8.315287e+00_kr,1.370959e+01_kr,2.260329e+01_kr,4.016900e+01_kr,&
     6.790405e+01_kr,9.166088e+01_kr,1.486254e+02_kr,3.043248e+02_kr,&
     4.539993e+02_kr,7.485183e+02_kr,1.234098e+03_kr,2.034684e+03_kr,&
     3.354626e+03_kr,5.530844e+03_kr,9.118820e+03_kr,1.503439e+04_kr,&
     2.478752e+04_kr,4.086771e+04_kr,6.737947e+04_kr,1.110900e+05_kr,&
     1.831564e+05_kr,3.019738e+05_kr,4.978707e+05_kr,8.208500e+05_kr,&
     1.353353e+06_kr,2.231302e+06_kr,3.678794e+06_kr,6.065307e+06_kr,&
     1.000000e+07_kr,1.964033e+07_kr/)
   real(kr),dimension(392),parameter::eg20a=(/&
     1.000010e-05_kr,3.000000e-03_kr,5.000000e-03_kr,6.900000e-03_kr,&
     1.000000e-02_kr,1.500000e-02_kr,2.000000e-02_kr,2.500000e-02_kr,&
     3.000000e-02_kr,3.500000e-02_kr,4.200000e-02_kr,5.000000e-02_kr,&
     5.800000e-02_kr,6.700000e-02_kr,7.700000e-02_kr,8.000000e-02_kr,&
     9.500000e-02_kr,1.000000e-01_kr,1.150000e-01_kr,1.340000e-01_kr,&
     1.400000e-01_kr,1.463700e-01_kr,1.530300e-01_kr,1.600000e-01_kr,&
     1.697100e-01_kr,1.800000e-01_kr,1.890000e-01_kr,1.988100e-01_kr,&
     2.091400e-01_kr,2.200000e-01_kr,2.335800e-01_kr,2.480000e-01_kr,&
     2.635100e-01_kr,2.800000e-01_kr,3.000000e-01_kr,3.145000e-01_kr,&
     3.200000e-01_kr,3.346600e-01_kr,3.500000e-01_kr,3.699300e-01_kr,&
     3.910000e-01_kr,4.000000e-01_kr,4.139900e-01_kr,4.330000e-01_kr,&
     4.496800e-01_kr,4.670100e-01_kr,4.850000e-01_kr,5.000000e-01_kr,&
     5.196200e-01_kr,5.315800e-01_kr,5.400000e-01_kr,5.669600e-01_kr,&
     5.952800e-01_kr,6.250000e-01_kr,6.531500e-01_kr,6.825600e-01_kr,&
     7.050000e-01_kr,7.415500e-01_kr,7.800000e-01_kr,7.900000e-01_kr,&
     8.194500e-01_kr,8.500000e-01_kr,8.600000e-01_kr,8.764250e-01_kr,&
     9.100000e-01_kr,9.300000e-01_kr,9.500000e-01_kr,9.720000e-01_kr,&
     9.860000e-01_kr,9.960000e-01_kr,1.020000e+00_kr,1.035000e+00_kr,&
     1.045000e+00_kr,1.071000e+00_kr,1.080000e+00_kr,1.097000e+00_kr,&
     1.110000e+00_kr,1.123000e+00_kr,1.150000e+00_kr,1.170000e+00_kr,&
     1.202060e+00_kr,1.235000e+00_kr,1.267080e+00_kr,1.300000e+00_kr,&
     1.337500e+00_kr,1.370000e+00_kr,1.404560e+00_kr,1.440000e+00_kr,&
     1.475000e+00_kr,1.500000e+00_kr,1.544340e+00_kr,1.590000e+00_kr,&
     1.629510e+00_kr,1.670000e+00_kr,1.711970e+00_kr,1.755000e+00_kr,&
     1.797000e+00_kr,1.840000e+00_kr,1.855390e+00_kr,1.884460e+00_kr,&
     1.930000e+00_kr,1.974490e+00_kr,2.020000e+00_kr,2.059610e+00_kr,&
     2.100000e+00_kr,2.130000e+00_kr,2.185310e+00_kr,2.242050e+00_kr,&
     2.300270e+00_kr,2.360000e+00_kr,2.382370e+00_kr,2.421710e+00_kr,&
     2.485030e+00_kr,2.550000e+00_kr,2.600000e+00_kr,2.659320e+00_kr,&
     2.720000e+00_kr,2.767920e+00_kr,2.837990e+00_kr,2.909830e+00_kr,&
     2.983490e+00_kr,3.059020e+00_kr,3.137330e+00_kr,3.217630e+00_kr,&
     3.300000e+00_kr,3.380750e+00_kr,3.466330e+00_kr,3.554080e+00_kr,&
     3.644050e+00_kr,3.736300e+00_kr,3.830880e+00_kr,3.927860e+00_kr,&
     4.000000e+00_kr,4.129250e+00_kr,4.233782e+00_kr,4.340961e+00_kr,&
     4.450853e+00_kr,4.563526e+00_kr,4.679053e+00_kr,4.797503e+00_kr,&
     4.918953e+00_kr,5.043477e+00_kr,5.085681e+00_kr,5.128239e+00_kr,&
     5.171153e+00_kr,5.214426e+00_kr,5.258061e+00_kr,5.302061e+00_kr,&
     5.346430e+00_kr,5.391169e+00_kr,5.436284e+00_kr,5.481775e+00_kr,&
     5.527647e+00_kr,5.573904e+00_kr,5.620547e+00_kr,5.667581e+00_kr,&
     5.715008e+00_kr,5.762832e+00_kr,5.811056e+00_kr,5.859684e+00_kr,&
     5.908719e+00_kr,5.958164e+00_kr,6.008022e+00_kr,6.058298e+00_kr,&
     6.108995e+00_kr,6.160116e+00_kr,6.211665e+00_kr,6.263645e+00_kr,&
     6.316060e+00_kr,6.368914e+00_kr,6.422210e+00_kr,6.475952e+00_kr,&
     6.530144e+00_kr,6.584789e+00_kr,6.639892e+00_kr,6.695455e+00_kr,&
     6.751484e+00_kr,6.807981e+00_kr,6.864952e+00_kr,6.922399e+00_kr,&
     6.980326e+00_kr,7.038739e+00_kr,7.097640e+00_kr,7.157034e+00_kr,&
     7.216925e+00_kr,7.277317e+00_kr,7.338215e+00_kr,7.399622e+00_kr,&
     7.461544e+00_kr,7.523983e+00_kr,7.586945e+00_kr,7.650434e+00_kr,&
     7.714454e+00_kr,7.779009e+00_kr,7.844105e+00_kr,7.909746e+00_kr,&
     7.975936e+00_kr,8.042680e+00_kr,8.109982e+00_kr,8.177848e+00_kr,&
     8.246281e+00_kr,8.315287e+00_kr,8.384871e+00_kr,8.455037e+00_kr,&
     8.525790e+00_kr,8.597135e+00_kr,8.669077e+00_kr,8.741621e+00_kr,&
     8.814772e+00_kr,8.888536e+00_kr,8.962916e+00_kr,9.037919e+00_kr,&
     9.113550e+00_kr,9.189814e+00_kr,9.266715e+00_kr,9.344261e+00_kr,&
     9.422455e+00_kr,9.501303e+00_kr,9.580812e+00_kr,9.660985e+00_kr,&
     9.741830e+00_kr,9.823351e+00_kr,9.905554e+00_kr,9.988446e+00_kr,&
     1.007203e+01_kr,1.015631e+01_kr,1.024130e+01_kr,1.032701e+01_kr,&
     1.041342e+01_kr,1.050056e+01_kr,1.058843e+01_kr,1.067704e+01_kr,&
     1.076639e+01_kr,1.085648e+01_kr,1.094733e+01_kr,1.103894e+01_kr,&
     1.113132e+01_kr,1.122446e+01_kr,1.131839e+01_kr,1.141311e+01_kr,&
     1.150861e+01_kr,1.160492e+01_kr,1.170203e+01_kr,1.179995e+01_kr,&
     1.189870e+01_kr,1.199827e+01_kr,1.209867e+01_kr,1.219991e+01_kr,&
     1.230201e+01_kr,1.240495e+01_kr,1.250876e+01_kr,1.261343e+01_kr,&
     1.271898e+01_kr,1.282542e+01_kr,1.293274e+01_kr,1.304097e+01_kr,&
     1.315010e+01_kr,1.326014e+01_kr,1.337110e+01_kr,1.348299e+01_kr,&
     1.359582e+01_kr,1.370959e+01_kr,1.382431e+01_kr,1.394000e+01_kr,&
     1.405665e+01_kr,1.417428e+01_kr,1.429289e+01_kr,1.441250e+01_kr,&
     1.453310e+01_kr,1.465472e+01_kr,1.477735e+01_kr,1.490101e+01_kr,&
     1.502570e+01_kr,1.515144e+01_kr,1.527823e+01_kr,1.540608e+01_kr,&
     1.553500e+01_kr,1.566500e+01_kr,1.579609e+01_kr,1.592827e+01_kr,&
     1.606156e+01_kr,1.619597e+01_kr,1.633150e+01_kr,1.646816e+01_kr,&
     1.660597e+01_kr,1.674493e+01_kr,1.688506e+01_kr,1.702635e+01_kr,&
     1.716883e+01_kr,1.731250e+01_kr,1.745738e+01_kr,1.760346e+01_kr,&
     1.775077e+01_kr,1.789931e+01_kr,1.804910e+01_kr,1.820013e+01_kr,&
     1.835244e+01_kr,1.850601e+01_kr,1.866087e+01_kr,1.881703e+01_kr,&
     1.897449e+01_kr,1.913328e+01_kr,1.929339e+01_kr,1.945484e+01_kr,&
     1.961764e+01_kr,1.978180e+01_kr,1.994734e+01_kr,2.011426e+01_kr,&
     2.028258e+01_kr,2.045231e+01_kr,2.062345e+01_kr,2.079603e+01_kr,&
     2.097006e+01_kr,2.114554e+01_kr,2.132249e+01_kr,2.150092e+01_kr,&
     2.168084e+01_kr,2.186227e+01_kr,2.204522e+01_kr,2.222969e+01_kr,&
     2.241572e+01_kr,2.260329e+01_kr,2.279244e+01_kr,2.298317e+01_kr,&
     2.317550e+01_kr,2.336944e+01_kr,2.356499e+01_kr,2.376219e+01_kr,&
     2.396104e+01_kr,2.416154e+01_kr,2.436373e+01_kr,2.456761e+01_kr,&
     2.477320e+01_kr,2.498050e+01_kr,2.518954e+01_kr,2.540033e+01_kr,&
     2.561289e+01_kr,2.582722e+01_kr,2.604335e+01_kr,2.626128e+01_kr,&
     2.648104e+01_kr,2.670264e+01_kr,2.692609e+01_kr,2.715141e+01_kr,&
     2.737862e+01_kr,2.760773e+01_kr,2.783875e+01_kr,2.807171e+01_kr,&
     2.830662e+01_kr,2.854349e+01_kr,2.878235e+01_kr,2.902320e+01_kr,&
     2.926607e+01_kr,2.951098e+01_kr,2.975793e+01_kr,3.000695e+01_kr,&
     3.025805e+01_kr,3.051126e+01_kr,3.076658e+01_kr,3.102404e+01_kr,&
     3.128365e+01_kr,3.154544e+01_kr,3.180942e+01_kr,3.207560e+01_kr,&
     3.234401e+01_kr,3.261467e+01_kr,3.288760e+01_kr,3.316281e+01_kr,&
     3.344032e+01_kr,3.372015e+01_kr,3.400233e+01_kr,3.428686e+01_kr,&
     3.457378e+01_kr,3.486310e+01_kr,3.515484e+01_kr,3.544902e+01_kr,&
     3.574566e+01_kr,3.604479e+01_kr,3.634642e+01_kr,3.665057e+01_kr,&
     3.695727e+01_kr,3.726653e+01_kr,3.757838e+01_kr,3.789285e+01_kr,&
     3.820994e+01_kr,3.852969e+01_kr,3.885211e+01_kr,3.917723e+01_kr,&
     3.950507e+01_kr,3.983565e+01_kr,4.016900e+01_kr,4.050514e+01_kr/)
   real(kr),dimension(392),parameter::eg20b=(/&
     4.084410e+01_kr,4.118589e+01_kr,4.153054e+01_kr,4.187807e+01_kr,&
     4.222851e+01_kr,4.258189e+01_kr,4.293822e+01_kr,4.329753e+01_kr,&
     4.365985e+01_kr,4.402521e+01_kr,4.439361e+01_kr,4.476511e+01_kr,&
     4.513971e+01_kr,4.551744e+01_kr,4.589834e+01_kr,4.628243e+01_kr,&
     4.666972e+01_kr,4.706026e+01_kr,4.745407e+01_kr,4.785117e+01_kr,&
     4.825160e+01_kr,4.865538e+01_kr,4.906253e+01_kr,4.947309e+01_kr,&
     4.988709e+01_kr,5.030456e+01_kr,5.072551e+01_kr,5.114999e+01_kr,&
     5.157802e+01_kr,5.200963e+01_kr,5.244486e+01_kr,5.288373e+01_kr,&
     5.332626e+01_kr,5.377251e+01_kr,5.422248e+01_kr,5.467623e+01_kr,&
     5.513376e+01_kr,5.559513e+01_kr,5.606036e+01_kr,5.652948e+01_kr,&
     5.700253e+01_kr,5.747954e+01_kr,5.796053e+01_kr,5.844556e+01_kr,&
     5.893464e+01_kr,5.942781e+01_kr,5.992511e+01_kr,6.042657e+01_kr,&
     6.093223e+01_kr,6.144212e+01_kr,6.195628e+01_kr,6.247474e+01_kr,&
     6.299754e+01_kr,6.352471e+01_kr,6.405630e+01_kr,6.459233e+01_kr,&
     6.513285e+01_kr,6.567789e+01_kr,6.622749e+01_kr,6.678169e+01_kr,&
     6.734053e+01_kr,6.790405e+01_kr,6.847228e+01_kr,6.904527e+01_kr,&
     6.962305e+01_kr,7.020566e+01_kr,7.079316e+01_kr,7.138556e+01_kr,&
     7.198293e+01_kr,7.258529e+01_kr,7.319270e+01_kr,7.380518e+01_kr,&
     7.442280e+01_kr,7.504558e+01_kr,7.567357e+01_kr,7.630682e+01_kr,&
     7.694537e+01_kr,7.758926e+01_kr,7.823854e+01_kr,7.889325e+01_kr,&
     7.955344e+01_kr,8.021915e+01_kr,8.089044e+01_kr,8.156734e+01_kr,&
     8.224991e+01_kr,8.293819e+01_kr,8.363223e+01_kr,8.433208e+01_kr,&
     8.503778e+01_kr,8.574939e+01_kr,8.646695e+01_kr,8.719052e+01_kr,&
     8.792015e+01_kr,8.865588e+01_kr,8.939776e+01_kr,9.014586e+01_kr,&
     9.090021e+01_kr,9.166088e+01_kr,9.242791e+01_kr,9.320136e+01_kr,&
     9.398128e+01_kr,9.476773e+01_kr,9.556076e+01_kr,9.636043e+01_kr,&
     9.716679e+01_kr,9.797990e+01_kr,9.879981e+01_kr,9.962658e+01_kr,&
     1.004603e+02_kr,1.013009e+02_kr,1.021486e+02_kr,1.030034e+02_kr,&
     1.038654e+02_kr,1.047345e+02_kr,1.056110e+02_kr,1.064947e+02_kr,&
     1.073859e+02_kr,1.082845e+02_kr,1.091907e+02_kr,1.101044e+02_kr,&
     1.110258e+02_kr,1.119548e+02_kr,1.128917e+02_kr,1.138364e+02_kr,&
     1.147890e+02_kr,1.157496e+02_kr,1.167182e+02_kr,1.176949e+02_kr,&
     1.186798e+02_kr,1.196729e+02_kr,1.206744e+02_kr,1.216842e+02_kr,&
     1.227024e+02_kr,1.237292e+02_kr,1.247646e+02_kr,1.258087e+02_kr,&
     1.268615e+02_kr,1.279231e+02_kr,1.289935e+02_kr,1.300730e+02_kr,&
     1.311615e+02_kr,1.322590e+02_kr,1.333658e+02_kr,1.344818e+02_kr,&
     1.356072e+02_kr,1.367420e+02_kr,1.378862e+02_kr,1.390401e+02_kr,&
     1.402036e+02_kr,1.413768e+02_kr,1.425599e+02_kr,1.437529e+02_kr,&
     1.449558e+02_kr,1.461688e+02_kr,1.473920e+02_kr,1.486254e+02_kr,&
     1.498691e+02_kr,1.511232e+02_kr,1.523879e+02_kr,1.536631e+02_kr,&
     1.549489e+02_kr,1.562456e+02_kr,1.575531e+02_kr,1.588715e+02_kr,&
     1.602010e+02_kr,1.615415e+02_kr,1.628933e+02_kr,1.642565e+02_kr,&
     1.656310e+02_kr,1.670170e+02_kr,1.684146e+02_kr,1.698239e+02_kr,&
     1.712451e+02_kr,1.726781e+02_kr,1.741231e+02_kr,1.755802e+02_kr,&
     1.770494e+02_kr,1.785310e+02_kr,1.800250e+02_kr,1.815315e+02_kr,&
     1.830505e+02_kr,1.845823e+02_kr,1.861269e+02_kr,1.876845e+02_kr,&
     1.892551e+02_kr,1.908388e+02_kr,1.924358e+02_kr,1.940461e+02_kr,&
     1.956699e+02_kr,1.973073e+02_kr,1.989584e+02_kr,2.006233e+02_kr,&
     2.023021e+02_kr,2.039950e+02_kr,2.057021e+02_kr,2.074234e+02_kr,&
     2.091592e+02_kr,2.109095e+02_kr,2.126744e+02_kr,2.144541e+02_kr,&
     2.162487e+02_kr,2.180583e+02_kr,2.198830e+02_kr,2.217230e+02_kr,&
     2.235784e+02_kr,2.254494e+02_kr,2.273360e+02_kr,2.292384e+02_kr,&
     2.311567e+02_kr,2.330910e+02_kr,2.350416e+02_kr,2.370084e+02_kr,&
     2.389917e+02_kr,2.409917e+02_kr,2.430083e+02_kr,2.450418e+02_kr,&
     2.470924e+02_kr,2.491601e+02_kr,2.512451e+02_kr,2.533476e+02_kr,&
     2.554676e+02_kr,2.576054e+02_kr,2.597611e+02_kr,2.619348e+02_kr,&
     2.641267e+02_kr,2.663370e+02_kr,2.685657e+02_kr,2.708131e+02_kr,&
     2.730793e+02_kr,2.753645e+02_kr,2.776688e+02_kr,2.799924e+02_kr,&
     2.823354e+02_kr,2.846980e+02_kr,2.870804e+02_kr,2.894827e+02_kr,&
     2.919052e+02_kr,2.943479e+02_kr,2.968110e+02_kr,2.992948e+02_kr,&
     3.017993e+02_kr,3.043248e+02_kr,3.068715e+02_kr,3.094394e+02_kr,&
     3.120288e+02_kr,3.146399e+02_kr,3.172729e+02_kr,3.199279e+02_kr,&
     3.226051e+02_kr,3.253047e+02_kr,3.280269e+02_kr,3.307719e+02_kr,&
     3.335398e+02_kr,3.363309e+02_kr,3.391454e+02_kr,3.419834e+02_kr,&
     3.448452e+02_kr,3.477309e+02_kr,3.506408e+02_kr,3.535750e+02_kr,&
     3.565338e+02_kr,3.595173e+02_kr,3.625258e+02_kr,3.655595e+02_kr,&
     3.686185e+02_kr,3.717032e+02_kr,3.748137e+02_kr,3.779502e+02_kr,&
     3.811129e+02_kr,3.843021e+02_kr,3.875180e+02_kr,3.907608e+02_kr,&
     3.940308e+02_kr,3.973281e+02_kr,4.006530e+02_kr,4.040057e+02_kr,&
     4.073865e+02_kr,4.107955e+02_kr,4.142332e+02_kr,4.176995e+02_kr,&
     4.211949e+02_kr,4.247195e+02_kr,4.282736e+02_kr,4.318575e+02_kr,&
     4.354713e+02_kr,4.391154e+02_kr,4.427900e+02_kr,4.464953e+02_kr,&
     4.502317e+02_kr,4.539993e+02_kr,4.577984e+02_kr,4.616294e+02_kr,&
     4.654923e+02_kr,4.693877e+02_kr,4.733156e+02_kr,4.772763e+02_kr,&
     4.812703e+02_kr,4.852976e+02_kr,4.893587e+02_kr,4.934537e+02_kr,&
     4.975830e+02_kr,5.017468e+02_kr,5.059455e+02_kr,5.101793e+02_kr,&
     5.144486e+02_kr,5.187536e+02_kr,5.230946e+02_kr,5.274719e+02_kr,&
     5.318859e+02_kr,5.363368e+02_kr,5.408249e+02_kr,5.453506e+02_kr,&
     5.499142e+02_kr,5.545160e+02_kr,5.591563e+02_kr,5.638354e+02_kr,&
     5.685536e+02_kr,5.733114e+02_kr,5.781089e+02_kr,5.829466e+02_kr,&
     5.878248e+02_kr,5.927438e+02_kr,5.977040e+02_kr,6.027057e+02_kr,&
     6.077492e+02_kr,6.128350e+02_kr,6.179633e+02_kr,6.231345e+02_kr,&
     6.283489e+02_kr,6.336071e+02_kr,6.389092e+02_kr,6.442557e+02_kr,&
     6.496469e+02_kr,6.550832e+02_kr,6.605651e+02_kr,6.660928e+02_kr,&
     6.716668e+02_kr,6.772874e+02_kr,6.829550e+02_kr,6.886701e+02_kr,&
     6.944330e+02_kr,7.002441e+02_kr,7.061038e+02_kr,7.120126e+02_kr,&
     7.179709e+02_kr,7.239790e+02_kr,7.300373e+02_kr,7.361464e+02_kr,&
     7.423066e+02_kr,7.485183e+02_kr,7.547820e+02_kr,7.610981e+02_kr,&
     7.674671e+02_kr,7.738894e+02_kr,7.803654e+02_kr,7.868957e+02_kr,&
     7.934805e+02_kr,8.001205e+02_kr,8.068160e+02_kr,8.135676e+02_kr,&
     8.203756e+02_kr,8.272407e+02_kr,8.341631e+02_kr,8.411435e+02_kr,&
     8.481824e+02_kr,8.552801e+02_kr,8.624372e+02_kr,8.696542e+02_kr,&
     8.769316e+02_kr,8.842699e+02_kr,8.916696e+02_kr,8.991312e+02_kr,&
     9.066553e+02_kr,9.142423e+02_kr,9.218928e+02_kr,9.296074e+02_kr,&
     9.373865e+02_kr,9.452307e+02_kr,9.531405e+02_kr,9.611165e+02_kr,&
     9.691593e+02_kr,9.772694e+02_kr,9.854473e+02_kr,9.936937e+02_kr,&
     1.002009e+03_kr,1.010394e+03_kr,1.018849e+03_kr,1.027375e+03_kr,&
     1.035972e+03_kr,1.044641e+03_kr,1.053383e+03_kr,1.062198e+03_kr/)
   real(kr),dimension(392),parameter::eg20c=(/&
     1.071087e+03_kr,1.080050e+03_kr,1.089088e+03_kr,1.098201e+03_kr,&
     1.107391e+03_kr,1.116658e+03_kr,1.126002e+03_kr,1.135425e+03_kr,&
     1.144926e+03_kr,1.154507e+03_kr,1.164168e+03_kr,1.173910e+03_kr,&
     1.183734e+03_kr,1.193639e+03_kr,1.203628e+03_kr,1.213700e+03_kr,&
     1.223857e+03_kr,1.234098e+03_kr,1.244425e+03_kr,1.254839e+03_kr,&
     1.265339e+03_kr,1.275928e+03_kr,1.286605e+03_kr,1.297372e+03_kr,&
     1.308228e+03_kr,1.319176e+03_kr,1.330215e+03_kr,1.341346e+03_kr,&
     1.352571e+03_kr,1.363889e+03_kr,1.375303e+03_kr,1.386811e+03_kr,&
     1.398416e+03_kr,1.410118e+03_kr,1.421919e+03_kr,1.433817e+03_kr,&
     1.445816e+03_kr,1.457915e+03_kr,1.470115e+03_kr,1.482417e+03_kr,&
     1.494822e+03_kr,1.507331e+03_kr,1.519944e+03_kr,1.532663e+03_kr,&
     1.545489e+03_kr,1.558422e+03_kr,1.571463e+03_kr,1.584613e+03_kr,&
     1.597874e+03_kr,1.611245e+03_kr,1.624728e+03_kr,1.638324e+03_kr,&
     1.652034e+03_kr,1.665858e+03_kr,1.679798e+03_kr,1.693855e+03_kr,&
     1.708030e+03_kr,1.722323e+03_kr,1.736735e+03_kr,1.751268e+03_kr,&
     1.765923e+03_kr,1.780701e+03_kr,1.795602e+03_kr,1.810628e+03_kr,&
     1.825780e+03_kr,1.841058e+03_kr,1.856464e+03_kr,1.871999e+03_kr,&
     1.887665e+03_kr,1.903461e+03_kr,1.919389e+03_kr,1.935451e+03_kr,&
     1.951647e+03_kr,1.967979e+03_kr,1.984447e+03_kr,2.001053e+03_kr,&
     2.017798e+03_kr,2.034684e+03_kr,2.051710e+03_kr,2.068879e+03_kr,&
     2.086192e+03_kr,2.103650e+03_kr,2.121253e+03_kr,2.139004e+03_kr,&
     2.156904e+03_kr,2.174953e+03_kr,2.193153e+03_kr,2.211506e+03_kr,&
     2.230012e+03_kr,2.248673e+03_kr,2.267490e+03_kr,2.286465e+03_kr,&
     2.305599e+03_kr,2.324892e+03_kr,2.344347e+03_kr,2.363965e+03_kr,&
     2.383747e+03_kr,2.403695e+03_kr,2.423809e+03_kr,2.444092e+03_kr,&
     2.464545e+03_kr,2.485168e+03_kr,2.505965e+03_kr,2.526935e+03_kr,&
     2.548081e+03_kr,2.569403e+03_kr,2.590904e+03_kr,2.612586e+03_kr,&
     2.634448e+03_kr,2.656494e+03_kr,2.678723e+03_kr,2.701139e+03_kr,&
     2.723743e+03_kr,2.746536e+03_kr,2.769519e+03_kr,2.792695e+03_kr,&
     2.816065e+03_kr,2.839630e+03_kr,2.863392e+03_kr,2.887354e+03_kr,&
     2.911515e+03_kr,2.935879e+03_kr,2.960447e+03_kr,2.985221e+03_kr,&
     3.010202e+03_kr,3.035391e+03_kr,3.060792e+03_kr,3.086405e+03_kr,&
     3.112233e+03_kr,3.138276e+03_kr,3.164538e+03_kr,3.191019e+03_kr,&
     3.217722e+03_kr,3.244649e+03_kr,3.271800e+03_kr,3.299179e+03_kr,&
     3.326787e+03_kr,3.354626e+03_kr,3.382698e+03_kr,3.411005e+03_kr,&
     3.439549e+03_kr,3.468332e+03_kr,3.497355e+03_kr,3.526622e+03_kr,&
     3.556133e+03_kr,3.585891e+03_kr,3.615898e+03_kr,3.646157e+03_kr,&
     3.676668e+03_kr,3.707435e+03_kr,3.738460e+03_kr,3.769744e+03_kr,&
     3.801290e+03_kr,3.833099e+03_kr,3.865175e+03_kr,3.897520e+03_kr,&
     3.930135e+03_kr,3.963023e+03_kr,3.996186e+03_kr,4.029627e+03_kr,&
     4.063347e+03_kr,4.097350e+03_kr,4.131637e+03_kr,4.166211e+03_kr,&
     4.201075e+03_kr,4.236230e+03_kr,4.271679e+03_kr,4.307425e+03_kr,&
     4.343471e+03_kr,4.379817e+03_kr,4.416468e+03_kr,4.453426e+03_kr,&
     4.490693e+03_kr,4.528272e+03_kr,4.566165e+03_kr,4.604375e+03_kr,&
     4.642906e+03_kr,4.681758e+03_kr,4.720936e+03_kr,4.760441e+03_kr,&
     4.800277e+03_kr,4.840447e+03_kr,4.880952e+03_kr,4.921797e+03_kr,&
     4.962983e+03_kr,5.004514e+03_kr,5.046393e+03_kr,5.088622e+03_kr,&
     5.131204e+03_kr,5.174143e+03_kr,5.217441e+03_kr,5.261101e+03_kr,&
     5.305127e+03_kr,5.349521e+03_kr,5.394287e+03_kr,5.439427e+03_kr,&
     5.484945e+03_kr,5.530844e+03_kr,5.577127e+03_kr,5.623797e+03_kr,&
     5.670858e+03_kr,5.718312e+03_kr,5.766164e+03_kr,5.814416e+03_kr,&
     5.863072e+03_kr,5.912135e+03_kr,5.961609e+03_kr,6.011496e+03_kr,&
     6.061802e+03_kr,6.112528e+03_kr,6.163678e+03_kr,6.215257e+03_kr,&
     6.267267e+03_kr,6.319712e+03_kr,6.372597e+03_kr,6.425924e+03_kr,&
     6.479697e+03_kr,6.533920e+03_kr,6.588597e+03_kr,6.643731e+03_kr,&
     6.699327e+03_kr,6.755388e+03_kr,6.811918e+03_kr,6.868921e+03_kr,&
     6.926401e+03_kr,6.984362e+03_kr,7.042809e+03_kr,7.101744e+03_kr,&
     7.161172e+03_kr,7.221098e+03_kr,7.281525e+03_kr,7.342458e+03_kr,&
     7.403901e+03_kr,7.465858e+03_kr,7.528334e+03_kr,7.591332e+03_kr,&
     7.654857e+03_kr,7.718914e+03_kr,7.783507e+03_kr,7.848641e+03_kr,&
     7.914319e+03_kr,7.980548e+03_kr,8.047330e+03_kr,8.114671e+03_kr,&
     8.182576e+03_kr,8.251049e+03_kr,8.320095e+03_kr,8.389719e+03_kr,&
     8.459926e+03_kr,8.530719e+03_kr,8.602106e+03_kr,8.674090e+03_kr,&
     8.746676e+03_kr,8.819869e+03_kr,8.893675e+03_kr,8.968099e+03_kr,&
     9.043145e+03_kr,9.118820e+03_kr,9.195127e+03_kr,9.272074e+03_kr,&
     9.349664e+03_kr,9.427903e+03_kr,9.506797e+03_kr,9.586352e+03_kr,&
     9.666572e+03_kr,9.747463e+03_kr,9.829031e+03_kr,9.911282e+03_kr,&
     9.994221e+03_kr,1.007785e+04_kr,1.016219e+04_kr,1.024723e+04_kr,&
     1.033298e+04_kr,1.041944e+04_kr,1.050664e+04_kr,1.059456e+04_kr,&
     1.068321e+04_kr,1.077261e+04_kr,1.086276e+04_kr,1.095366e+04_kr,&
     1.104532e+04_kr,1.113775e+04_kr,1.123095e+04_kr,1.132494e+04_kr,&
     1.141970e+04_kr,1.151527e+04_kr,1.161163e+04_kr,1.170880e+04_kr,&
     1.180678e+04_kr,1.190558e+04_kr,1.200521e+04_kr,1.210567e+04_kr,&
     1.220697e+04_kr,1.230912e+04_kr,1.241212e+04_kr,1.251599e+04_kr,&
     1.262073e+04_kr,1.272634e+04_kr,1.283283e+04_kr,1.294022e+04_kr,&
     1.304851e+04_kr,1.315770e+04_kr,1.326780e+04_kr,1.337883e+04_kr,&
     1.349079e+04_kr,1.360368e+04_kr,1.371752e+04_kr,1.383231e+04_kr,&
     1.394806e+04_kr,1.406478e+04_kr,1.418247e+04_kr,1.430116e+04_kr,&
     1.442083e+04_kr,1.454151e+04_kr,1.466319e+04_kr,1.478590e+04_kr,&
     1.490963e+04_kr,1.503439e+04_kr,1.516020e+04_kr,1.528706e+04_kr,&
     1.541499e+04_kr,1.554398e+04_kr,1.567406e+04_kr,1.580522e+04_kr,&
     1.593748e+04_kr,1.607085e+04_kr,1.620533e+04_kr,1.634094e+04_kr,&
     1.647768e+04_kr,1.661557e+04_kr,1.675461e+04_kr,1.689482e+04_kr,&
     1.703620e+04_kr,1.717876e+04_kr,1.732251e+04_kr,1.746747e+04_kr,&
     1.761364e+04_kr,1.776104e+04_kr,1.790966e+04_kr,1.805953e+04_kr,&
     1.821066e+04_kr,1.836305e+04_kr,1.851671e+04_kr,1.867166e+04_kr,&
     1.882791e+04_kr,1.898547e+04_kr,1.914434e+04_kr,1.930454e+04_kr,&
     1.946608e+04_kr,1.962898e+04_kr,1.979324e+04_kr,1.995887e+04_kr,&
     2.012589e+04_kr,2.029431e+04_kr,2.046413e+04_kr,2.063538e+04_kr,&
     2.080806e+04_kr,2.098218e+04_kr,2.115777e+04_kr,2.133482e+04_kr,&
     2.151335e+04_kr,2.169338e+04_kr,2.187491e+04_kr,2.205796e+04_kr,&
     2.224255e+04_kr,2.242868e+04_kr,2.261636e+04_kr,2.280562e+04_kr,&
     2.299646e+04_kr,2.318890e+04_kr,2.338295e+04_kr,2.357862e+04_kr,&
     2.377593e+04_kr,2.397489e+04_kr,2.417552e+04_kr,2.437782e+04_kr,&
     2.458182e+04_kr,2.478752e+04_kr,2.499495e+04_kr,2.520411e+04_kr,&
     2.541502e+04_kr,2.562770e+04_kr,2.584215e+04_kr,2.605841e+04_kr,&
     2.627647e+04_kr,2.649635e+04_kr,2.671808e+04_kr,2.694166e+04_kr,&
     2.700000e+04_kr,2.716711e+04_kr,2.739445e+04_kr,2.762369e+04_kr/)
   real(kr),dimension(392),parameter::eg20d=(/&
     2.785485e+04_kr,2.808794e+04_kr,2.832299e+04_kr,2.850000e+04_kr,&
     2.856000e+04_kr,2.879899e+04_kr,2.903999e+04_kr,2.928300e+04_kr,&
     2.952804e+04_kr,2.977514e+04_kr,3.002430e+04_kr,3.027555e+04_kr,&
     3.052890e+04_kr,3.078437e+04_kr,3.104198e+04_kr,3.130174e+04_kr,&
     3.156368e+04_kr,3.182781e+04_kr,3.209415e+04_kr,3.236272e+04_kr,&
     3.263353e+04_kr,3.290662e+04_kr,3.318198e+04_kr,3.345965e+04_kr,&
     3.373965e+04_kr,3.402199e+04_kr,3.430669e+04_kr,3.459377e+04_kr,&
     3.488326e+04_kr,3.517517e+04_kr,3.546952e+04_kr,3.576633e+04_kr,&
     3.606563e+04_kr,3.636743e+04_kr,3.667176e+04_kr,3.697864e+04_kr,&
     3.728808e+04_kr,3.760011e+04_kr,3.791476e+04_kr,3.823203e+04_kr,&
     3.855196e+04_kr,3.887457e+04_kr,3.919988e+04_kr,3.952791e+04_kr,&
     3.985869e+04_kr,4.019223e+04_kr,4.052857e+04_kr,4.086771e+04_kr,&
     4.120970e+04_kr,4.155455e+04_kr,4.190229e+04_kr,4.225293e+04_kr,&
     4.260651e+04_kr,4.296305e+04_kr,4.332257e+04_kr,4.368510e+04_kr,&
     4.405066e+04_kr,4.441928e+04_kr,4.479099e+04_kr,4.516581e+04_kr,&
     4.554376e+04_kr,4.592488e+04_kr,4.630919e+04_kr,4.669671e+04_kr,&
     4.708747e+04_kr,4.748151e+04_kr,4.787884e+04_kr,4.827950e+04_kr,&
     4.868351e+04_kr,4.909090e+04_kr,4.950170e+04_kr,4.991594e+04_kr,&
     5.033364e+04_kr,5.075484e+04_kr,5.117957e+04_kr,5.160785e+04_kr,&
     5.203971e+04_kr,5.247518e+04_kr,5.291430e+04_kr,5.335710e+04_kr,&
     5.380360e+04_kr,5.425384e+04_kr,5.470784e+04_kr,5.516564e+04_kr,&
     5.562728e+04_kr,5.609278e+04_kr,5.656217e+04_kr,5.703549e+04_kr,&
     5.751277e+04_kr,5.799405e+04_kr,5.847935e+04_kr,5.896871e+04_kr,&
     5.946217e+04_kr,5.995976e+04_kr,6.046151e+04_kr,6.096747e+04_kr,&
     6.147765e+04_kr,6.199211e+04_kr,6.251086e+04_kr,6.303396e+04_kr,&
     6.356144e+04_kr,6.409333e+04_kr,6.462968e+04_kr,6.517051e+04_kr,&
     6.571586e+04_kr,6.626579e+04_kr,6.682031e+04_kr,6.737947e+04_kr,&
     6.794331e+04_kr,6.851187e+04_kr,6.908519e+04_kr,6.966330e+04_kr,&
     7.024626e+04_kr,7.083409e+04_kr,7.142684e+04_kr,7.202455e+04_kr,&
     7.262726e+04_kr,7.323502e+04_kr,7.384786e+04_kr,7.446583e+04_kr,&
     7.508897e+04_kr,7.571733e+04_kr,7.635094e+04_kr,7.698986e+04_kr,&
     7.763412e+04_kr,7.828378e+04_kr,7.893887e+04_kr,7.950000e+04_kr,&
     7.959944e+04_kr,8.026554e+04_kr,8.093721e+04_kr,8.161451e+04_kr,&
     8.229747e+04_kr,8.250000e+04_kr,8.298615e+04_kr,8.368059e+04_kr,&
     8.438084e+04_kr,8.508695e+04_kr,8.579897e+04_kr,8.651695e+04_kr,&
     8.724094e+04_kr,8.797098e+04_kr,8.870714e+04_kr,8.944945e+04_kr,&
     9.019798e+04_kr,9.095277e+04_kr,9.171388e+04_kr,9.248135e+04_kr,&
     9.325525e+04_kr,9.403563e+04_kr,9.482253e+04_kr,9.561602e+04_kr,&
     9.641615e+04_kr,9.722297e+04_kr,9.803655e+04_kr,9.885694e+04_kr,&
     9.968419e+04_kr,1.005184e+05_kr,1.013595e+05_kr,1.022077e+05_kr,&
     1.030630e+05_kr,1.039254e+05_kr,1.047951e+05_kr,1.056720e+05_kr,&
     1.065563e+05_kr,1.074480e+05_kr,1.083471e+05_kr,1.092538e+05_kr,&
     1.101681e+05_kr,1.110900e+05_kr,1.120196e+05_kr,1.129570e+05_kr,&
     1.139022e+05_kr,1.148554e+05_kr,1.158165e+05_kr,1.167857e+05_kr,&
     1.177629e+05_kr,1.187484e+05_kr,1.197421e+05_kr,1.207441e+05_kr,&
     1.217545e+05_kr,1.227734e+05_kr,1.238008e+05_kr,1.248368e+05_kr,&
     1.258814e+05_kr,1.269348e+05_kr,1.279970e+05_kr,1.290681e+05_kr,&
     1.301482e+05_kr,1.312373e+05_kr,1.323355e+05_kr,1.334429e+05_kr,&
     1.345596e+05_kr,1.356856e+05_kr,1.368210e+05_kr,1.379660e+05_kr,&
     1.391205e+05_kr,1.402847e+05_kr,1.414586e+05_kr,1.426423e+05_kr,&
     1.438360e+05_kr,1.450396e+05_kr,1.462533e+05_kr,1.474772e+05_kr,&
     1.487113e+05_kr,1.499558e+05_kr,1.512106e+05_kr,1.524760e+05_kr,&
     1.537519e+05_kr,1.550385e+05_kr,1.563359e+05_kr,1.576442e+05_kr,&
     1.589634e+05_kr,1.602936e+05_kr,1.616349e+05_kr,1.629875e+05_kr,&
     1.643514e+05_kr,1.657268e+05_kr,1.671136e+05_kr,1.685120e+05_kr,&
     1.699221e+05_kr,1.713441e+05_kr,1.727779e+05_kr,1.742237e+05_kr,&
     1.756817e+05_kr,1.771518e+05_kr,1.786342e+05_kr,1.801291e+05_kr,&
     1.816364e+05_kr,1.831564e+05_kr,1.846891e+05_kr,1.862346e+05_kr,&
     1.877930e+05_kr,1.893645e+05_kr,1.909491e+05_kr,1.925470e+05_kr,&
     1.941583e+05_kr,1.957830e+05_kr,1.974214e+05_kr,1.990734e+05_kr,&
     2.007393e+05_kr,2.024191e+05_kr,2.041130e+05_kr,2.058210e+05_kr,&
     2.075434e+05_kr,2.092801e+05_kr,2.110314e+05_kr,2.127974e+05_kr,&
     2.145781e+05_kr,2.163737e+05_kr,2.181844e+05_kr,2.200102e+05_kr,&
     2.218512e+05_kr,2.237077e+05_kr,2.255797e+05_kr,2.274674e+05_kr,&
     2.293709e+05_kr,2.312903e+05_kr,2.332258e+05_kr,2.351775e+05_kr,&
     2.371455e+05_kr,2.391299e+05_kr,2.411310e+05_kr,2.431488e+05_kr,&
     2.451835e+05_kr,2.472353e+05_kr,2.493042e+05_kr,2.513904e+05_kr,&
     2.534941e+05_kr,2.556153e+05_kr,2.577544e+05_kr,2.599113e+05_kr,&
     2.620863e+05_kr,2.642794e+05_kr,2.664910e+05_kr,2.687210e+05_kr,&
     2.709697e+05_kr,2.732372e+05_kr,2.755237e+05_kr,2.778293e+05_kr,&
     2.801543e+05_kr,2.824986e+05_kr,2.848626e+05_kr,2.872464e+05_kr,&
     2.896501e+05_kr,2.920740e+05_kr,2.945181e+05_kr,2.969826e+05_kr,&
     2.972000e+05_kr,2.985000e+05_kr,2.994678e+05_kr,3.019738e+05_kr,&
     3.045008e+05_kr,3.070489e+05_kr,3.096183e+05_kr,3.122093e+05_kr,&
     3.148219e+05_kr,3.174564e+05_kr,3.201129e+05_kr,3.227916e+05_kr,&
     3.254928e+05_kr,3.282166e+05_kr,3.309631e+05_kr,3.337327e+05_kr,&
     3.365254e+05_kr,3.393415e+05_kr,3.421812e+05_kr,3.450446e+05_kr,&
     3.479320e+05_kr,3.508435e+05_kr,3.537795e+05_kr,3.567399e+05_kr,&
     3.597252e+05_kr,3.627354e+05_kr,3.657708e+05_kr,3.688317e+05_kr,&
     3.719181e+05_kr,3.750304e+05_kr,3.781687e+05_kr,3.813333e+05_kr,&
     3.845243e+05_kr,3.877421e+05_kr,3.909868e+05_kr,3.942586e+05_kr,&
     3.975578e+05_kr,4.008846e+05_kr,4.042393e+05_kr,4.076220e+05_kr,&
     4.110331e+05_kr,4.144727e+05_kr,4.179410e+05_kr,4.214384e+05_kr,&
     4.249651e+05_kr,4.285213e+05_kr,4.321072e+05_kr,4.357231e+05_kr,&
     4.393693e+05_kr,4.430460e+05_kr,4.467535e+05_kr,4.504920e+05_kr,&
     4.542618e+05_kr,4.580631e+05_kr,4.618963e+05_kr,4.657615e+05_kr,&
     4.696591e+05_kr,4.735892e+05_kr,4.775523e+05_kr,4.815485e+05_kr,&
     4.855782e+05_kr,4.896416e+05_kr,4.937390e+05_kr,4.978707e+05_kr,&
     5.020369e+05_kr,5.062381e+05_kr,5.104743e+05_kr,5.147461e+05_kr,&
     5.190535e+05_kr,5.233971e+05_kr,5.277769e+05_kr,5.321934e+05_kr,&
     5.366469e+05_kr,5.411377e+05_kr,5.456660e+05_kr,5.502322e+05_kr,&
     5.548366e+05_kr,5.594796e+05_kr,5.641614e+05_kr,5.688824e+05_kr,&
     5.736429e+05_kr,5.784432e+05_kr,5.832837e+05_kr,5.881647e+05_kr,&
     5.930866e+05_kr,5.980496e+05_kr,6.030542e+05_kr,6.081006e+05_kr,&
     6.131893e+05_kr,6.183206e+05_kr,6.234948e+05_kr,6.287123e+05_kr,&
     6.339734e+05_kr,6.392786e+05_kr,6.446282e+05_kr,6.500225e+05_kr,&
     6.554620e+05_kr,6.609470e+05_kr,6.664779e+05_kr,6.720551e+05_kr,&
     6.776790e+05_kr,6.833499e+05_kr,6.890683e+05_kr,6.948345e+05_kr/)
   real(kr),dimension(392),parameter::eg20e=(/&
     7.006490e+05_kr,7.065121e+05_kr,7.124243e+05_kr,7.183860e+05_kr,&
     7.243976e+05_kr,7.304594e+05_kr,7.365720e+05_kr,7.427358e+05_kr,&
     7.489511e+05_kr,7.552184e+05_kr,7.615382e+05_kr,7.679109e+05_kr,&
     7.743369e+05_kr,7.808167e+05_kr,7.873507e+05_kr,7.939393e+05_kr,&
     8.005831e+05_kr,8.072825e+05_kr,8.140380e+05_kr,8.208500e+05_kr,&
     8.277190e+05_kr,8.346455e+05_kr,8.416299e+05_kr,8.486728e+05_kr,&
     8.557746e+05_kr,8.629359e+05_kr,8.701570e+05_kr,8.774387e+05_kr,&
     8.847812e+05_kr,8.921852e+05_kr,8.996511e+05_kr,9.071795e+05_kr,&
     9.147709e+05_kr,9.224259e+05_kr,9.301449e+05_kr,9.379285e+05_kr,&
     9.457772e+05_kr,9.536916e+05_kr,9.616723e+05_kr,9.697197e+05_kr,&
     9.778344e+05_kr,9.860171e+05_kr,9.942682e+05_kr,1.002588e+06_kr,&
     1.010978e+06_kr,1.019438e+06_kr,1.027969e+06_kr,1.036571e+06_kr,&
     1.045245e+06_kr,1.053992e+06_kr,1.062812e+06_kr,1.071706e+06_kr,&
     1.080674e+06_kr,1.089717e+06_kr,1.098836e+06_kr,1.108032e+06_kr,&
     1.117304e+06_kr,1.126654e+06_kr,1.136082e+06_kr,1.145588e+06_kr,&
     1.155175e+06_kr,1.164842e+06_kr,1.174589e+06_kr,1.184418e+06_kr,&
     1.194330e+06_kr,1.204324e+06_kr,1.214402e+06_kr,1.224564e+06_kr,&
     1.234812e+06_kr,1.245145e+06_kr,1.255564e+06_kr,1.266071e+06_kr,&
     1.276666e+06_kr,1.287349e+06_kr,1.298122e+06_kr,1.308985e+06_kr,&
     1.319938e+06_kr,1.330984e+06_kr,1.342122e+06_kr,1.353353e+06_kr,&
     1.364678e+06_kr,1.376098e+06_kr,1.387613e+06_kr,1.399225e+06_kr,&
     1.410934e+06_kr,1.422741e+06_kr,1.434646e+06_kr,1.446652e+06_kr,&
     1.458758e+06_kr,1.470965e+06_kr,1.483274e+06_kr,1.495686e+06_kr,&
     1.508202e+06_kr,1.520823e+06_kr,1.533550e+06_kr,1.546383e+06_kr,&
     1.559323e+06_kr,1.572372e+06_kr,1.585530e+06_kr,1.598797e+06_kr,&
     1.612176e+06_kr,1.625667e+06_kr,1.639271e+06_kr,1.652989e+06_kr,&
     1.666821e+06_kr,1.680770e+06_kr,1.694834e+06_kr,1.709017e+06_kr,&
     1.723318e+06_kr,1.737739e+06_kr,1.752281e+06_kr,1.766944e+06_kr,&
     1.781731e+06_kr,1.796640e+06_kr,1.811675e+06_kr,1.826835e+06_kr,&
     1.842122e+06_kr,1.857538e+06_kr,1.873082e+06_kr,1.888756e+06_kr,&
     1.904561e+06_kr,1.920499e+06_kr,1.936570e+06_kr,1.952776e+06_kr,&
     1.969117e+06_kr,1.985595e+06_kr,2.002210e+06_kr,2.018965e+06_kr,&
     2.035860e+06_kr,2.052897e+06_kr,2.070076e+06_kr,2.087398e+06_kr,&
     2.104866e+06_kr,2.122480e+06_kr,2.140241e+06_kr,2.158151e+06_kr,&
     2.176211e+06_kr,2.194421e+06_kr,2.212785e+06_kr,2.231302e+06_kr,&
     2.249973e+06_kr,2.268802e+06_kr,2.287787e+06_kr,2.306932e+06_kr,&
     2.326237e+06_kr,2.345703e+06_kr,2.365332e+06_kr,2.385126e+06_kr,&
     2.405085e+06_kr,2.425211e+06_kr,2.445505e+06_kr,2.465970e+06_kr,&
     2.486605e+06_kr,2.507414e+06_kr,2.528396e+06_kr,2.549554e+06_kr,&
     2.570889e+06_kr,2.592403e+06_kr,2.614096e+06_kr,2.635971e+06_kr,&
     2.658030e+06_kr,2.680272e+06_kr,2.702701e+06_kr,2.725318e+06_kr,&
     2.748124e+06_kr,2.771121e+06_kr,2.794310e+06_kr,2.817693e+06_kr,&
     2.841272e+06_kr,2.865048e+06_kr,2.889023e+06_kr,2.913199e+06_kr,&
     2.937577e+06_kr,2.962159e+06_kr,2.986947e+06_kr,3.011942e+06_kr,&
     3.037147e+06_kr,3.062562e+06_kr,3.088190e+06_kr,3.114032e+06_kr,&
     3.140091e+06_kr,3.166368e+06_kr,3.192864e+06_kr,3.219583e+06_kr,&
     3.246525e+06_kr,3.273692e+06_kr,3.301087e+06_kr,3.328711e+06_kr,&
     3.356566e+06_kr,3.384654e+06_kr,3.412978e+06_kr,3.441538e+06_kr,&
     3.470337e+06_kr,3.499377e+06_kr,3.528661e+06_kr,3.558189e+06_kr,&
     3.587965e+06_kr,3.617989e+06_kr,3.648265e+06_kr,3.678794e+06_kr,&
     3.709579e+06_kr,3.740621e+06_kr,3.771924e+06_kr,3.803488e+06_kr,&
     3.835316e+06_kr,3.867410e+06_kr,3.899773e+06_kr,3.932407e+06_kr,&
     3.965314e+06_kr,3.998497e+06_kr,4.031957e+06_kr,4.065697e+06_kr,&
     4.099719e+06_kr,4.134026e+06_kr,4.168620e+06_kr,4.203504e+06_kr,&
     4.238679e+06_kr,4.274149e+06_kr,4.309916e+06_kr,4.345982e+06_kr,&
     4.382350e+06_kr,4.419022e+06_kr,4.456001e+06_kr,4.493290e+06_kr,&
     4.530890e+06_kr,4.568805e+06_kr,4.607038e+06_kr,4.645590e+06_kr,&
     4.684465e+06_kr,4.723666e+06_kr,4.763194e+06_kr,4.803053e+06_kr,&
     4.843246e+06_kr,4.883775e+06_kr,4.924643e+06_kr,4.965853e+06_kr,&
     5.007408e+06_kr,5.049311e+06_kr,5.091564e+06_kr,5.134171e+06_kr,&
     5.177135e+06_kr,5.220458e+06_kr,5.264143e+06_kr,5.308195e+06_kr,&
     5.352614e+06_kr,5.397406e+06_kr,5.442572e+06_kr,5.488116e+06_kr,&
     5.534042e+06_kr,5.580351e+06_kr,5.627049e+06_kr,5.674137e+06_kr,&
     5.721619e+06_kr,5.769498e+06_kr,5.817778e+06_kr,5.866462e+06_kr,&
     5.915554e+06_kr,5.965056e+06_kr,6.014972e+06_kr,6.065307e+06_kr,&
     6.116062e+06_kr,6.167242e+06_kr,6.218851e+06_kr,6.270891e+06_kr,&
     6.323367e+06_kr,6.376282e+06_kr,6.429639e+06_kr,6.483443e+06_kr,&
     6.537698e+06_kr,6.592406e+06_kr,6.647573e+06_kr,6.703200e+06_kr,&
     6.759294e+06_kr,6.815857e+06_kr,6.872893e+06_kr,6.930406e+06_kr,&
     6.988401e+06_kr,7.046881e+06_kr,7.105850e+06_kr,7.165313e+06_kr,&
     7.225274e+06_kr,7.285736e+06_kr,7.346704e+06_kr,7.408182e+06_kr,&
     7.470175e+06_kr,7.532687e+06_kr,7.595721e+06_kr,7.659283e+06_kr,&
     7.723377e+06_kr,7.788008e+06_kr,7.853179e+06_kr,7.918896e+06_kr,&
     7.985162e+06_kr,8.051983e+06_kr,8.119363e+06_kr,8.187308e+06_kr,&
     8.255820e+06_kr,8.324906e+06_kr,8.394570e+06_kr,8.464817e+06_kr,&
     8.535652e+06_kr,8.607080e+06_kr,8.679105e+06_kr,8.751733e+06_kr,&
     8.824969e+06_kr,8.898818e+06_kr,8.973284e+06_kr,9.048374e+06_kr,&
     9.124092e+06_kr,9.200444e+06_kr,9.277435e+06_kr,9.355070e+06_kr,&
     9.433354e+06_kr,9.512294e+06_kr,9.591895e+06_kr,9.672161e+06_kr,&
     9.753099e+06_kr,9.834715e+06_kr,9.917013e+06_kr,1.000000e+07_kr,&
     1.008368e+07_kr,1.016806e+07_kr,1.025315e+07_kr,1.033895e+07_kr,&
     1.042547e+07_kr,1.051271e+07_kr,1.060068e+07_kr,1.068939e+07_kr,&
     1.077884e+07_kr,1.086904e+07_kr,1.095999e+07_kr,1.105171e+07_kr,&
     1.114419e+07_kr,1.123745e+07_kr,1.133148e+07_kr,1.142631e+07_kr,&
     1.152193e+07_kr,1.161834e+07_kr,1.171557e+07_kr,1.181360e+07_kr,&
     1.191246e+07_kr,1.201215e+07_kr,1.211267e+07_kr,1.221403e+07_kr,&
     1.231624e+07_kr,1.241930e+07_kr,1.252323e+07_kr,1.262802e+07_kr,&
     1.273370e+07_kr,1.284025e+07_kr,1.294770e+07_kr,1.305605e+07_kr,&
     1.316531e+07_kr,1.327548e+07_kr,1.338657e+07_kr,1.349859e+07_kr,&
     1.361155e+07_kr,1.372545e+07_kr,1.384031e+07_kr,1.395612e+07_kr,&
     1.407291e+07_kr,1.419068e+07_kr,1.430943e+07_kr,1.442917e+07_kr,&
     1.454991e+07_kr,1.467167e+07_kr,1.479444e+07_kr,1.491825e+07_kr,&
     1.504309e+07_kr,1.516897e+07_kr,1.529590e+07_kr,1.542390e+07_kr,&
     1.555297e+07_kr,1.568312e+07_kr,1.581436e+07_kr,1.594670e+07_kr,&
     1.608014e+07_kr,1.621470e+07_kr,1.635039e+07_kr,1.648721e+07_kr,&
     1.662518e+07_kr,1.676430e+07_kr,1.690459e+07_kr,1.704605e+07_kr,&
     1.718869e+07_kr,1.733253e+07_kr,1.747757e+07_kr,1.762383e+07_kr,&
     1.777131e+07_kr,1.792002e+07_kr,1.806998e+07_kr,1.822119e+07_kr/)
   real(kr),dimension(9),parameter::eg20f=(/&
     1.837367e+07_kr,1.852742e+07_kr,1.868246e+07_kr,1.883880e+07_kr,&
     1.899644e+07_kr,1.915541e+07_kr,1.931570e+07_kr,1.947734e+07_kr,&
     1.964033e+07_kr/)
   real(kr),parameter::eg21(316)=(/&
     1.000010e-05_kr,1.100000e-04_kr,3.000000e-03_kr,5.500100e-03_kr,&
     1.000000e-02_kr,1.500000e-02_kr,2.000000e-02_kr,3.000000e-02_kr,&
     3.200000e-02_kr,3.238000e-02_kr,4.300000e-02_kr,5.900100e-02_kr,&
     7.700100e-02_kr,9.500000e-02_kr,1.000000e-01_kr,1.150000e-01_kr,&
     1.340000e-01_kr,1.600000e-01_kr,1.890000e-01_kr,2.200000e-01_kr,&
     2.480000e-01_kr,2.825000e-01_kr,3.145000e-01_kr,3.520000e-01_kr,&
     3.910100e-01_kr,4.139900e-01_kr,4.330000e-01_kr,4.850100e-01_kr,&
     5.315800e-01_kr,5.400100e-01_kr,6.250100e-01_kr,6.825600e-01_kr,&
     7.050000e-01_kr,7.900100e-01_kr,8.600100e-01_kr,8.764200e-01_kr,&
     9.300100e-01_kr,9.860100e-01_kr,1.010000e+00_kr,1.035000e+00_kr,&
     1.070000e+00_kr,1.080000e+00_kr,1.090000e+00_kr,1.110000e+00_kr,&
     1.125400e+00_kr,1.170000e+00_kr,1.235000e+00_kr,1.305000e+00_kr,&
     1.370000e+00_kr,1.440000e+00_kr,1.445000e+00_kr,1.510000e+00_kr,&
     1.590000e+00_kr,1.670000e+00_kr,1.755000e+00_kr,1.840000e+00_kr,&
     1.855400e+00_kr,1.930000e+00_kr,2.020000e+00_kr,2.130000e+00_kr,&
     2.360000e+00_kr,2.372400e+00_kr,2.767900e+00_kr,3.059000e+00_kr,&
     3.380700e+00_kr,3.927900e+00_kr,4.129200e+00_kr,4.470000e+00_kr,&
     4.670000e+00_kr,5.043500e+00_kr,5.623000e+00_kr,6.160100e+00_kr,&
     6.476000e+00_kr,7.079000e+00_kr,7.524000e+00_kr,7.943000e+00_kr,&
     8.315300e+00_kr,8.913000e+00_kr,9.189800e+00_kr,1.000000e+01_kr,&
     1.067700e+01_kr,1.122400e+01_kr,1.259000e+01_kr,1.371000e+01_kr,&
     1.522700e+01_kr,1.674500e+01_kr,1.760300e+01_kr,1.902800e+01_kr,&
     2.045200e+01_kr,2.260300e+01_kr,2.498000e+01_kr,2.791800e+01_kr,&
     2.920300e+01_kr,3.051100e+01_kr,3.388900e+01_kr,3.726700e+01_kr,&
     3.981000e+01_kr,4.551700e+01_kr,4.785100e+01_kr,5.012000e+01_kr,&
     5.559500e+01_kr,6.144200e+01_kr,6.310000e+01_kr,6.790400e+01_kr,&
     7.079000e+01_kr,7.889300e+01_kr,8.527700e+01_kr,9.166100e+01_kr,&
     1.013000e+02_kr,1.122000e+02_kr,1.300700e+02_kr,1.367400e+02_kr,&
     1.585000e+02_kr,1.670200e+02_kr,1.778000e+02_kr,2.039900e+02_kr,&
     2.144500e+02_kr,2.430100e+02_kr,2.753600e+02_kr,3.043200e+02_kr,&
     3.535800e+02_kr,3.981000e+02_kr,4.540000e+02_kr,5.144600e+02_kr,&
     5.829500e+02_kr,6.310000e+02_kr,6.772900e+02_kr,7.079000e+02_kr,&
     7.485200e+02_kr,8.482000e+02_kr,9.611200e+02_kr,1.010400e+03_kr,&
     1.116700e+03_kr,1.234100e+03_kr,1.363900e+03_kr,1.507300e+03_kr,&
     1.584600e+03_kr,1.795600e+03_kr,2.034700e+03_kr,2.113000e+03_kr,&
     2.248700e+03_kr,2.371000e+03_kr,2.485200e+03_kr,2.612600e+03_kr,&
     2.661000e+03_kr,2.746500e+03_kr,2.818000e+03_kr,3.035400e+03_kr,&
     3.162000e+03_kr,3.354600e+03_kr,3.548000e+03_kr,3.707400e+03_kr,&
     3.981000e+03_kr,4.307400e+03_kr,4.642900e+03_kr,5.004500e+03_kr,&
     5.530800e+03_kr,6.267300e+03_kr,7.101700e+03_kr,7.465900e+03_kr,&
     8.251000e+03_kr,9.118800e+03_kr,1.007800e+04_kr,1.113800e+04_kr,&
     1.170900e+04_kr,1.272600e+04_kr,1.383200e+04_kr,1.503400e+04_kr,&
     1.585000e+04_kr,1.661600e+04_kr,1.778000e+04_kr,1.930500e+04_kr,&
     1.995000e+04_kr,2.054000e+04_kr,2.113000e+04_kr,2.187500e+04_kr,&
     2.239000e+04_kr,2.304000e+04_kr,2.357900e+04_kr,2.417600e+04_kr,&
     2.441000e+04_kr,2.478800e+04_kr,2.512000e+04_kr,2.585000e+04_kr,&
     2.605800e+04_kr,2.661000e+04_kr,2.700000e+04_kr,2.738000e+04_kr,&
     2.818000e+04_kr,2.850000e+04_kr,2.901000e+04_kr,2.985000e+04_kr,&
     3.073000e+04_kr,3.162000e+04_kr,3.182800e+04_kr,3.430700e+04_kr,&
     3.697900e+04_kr,4.086800e+04_kr,4.358900e+04_kr,4.630900e+04_kr,&
     4.939200e+04_kr,5.247500e+04_kr,5.516600e+04_kr,5.656200e+04_kr,&
     6.172500e+04_kr,6.737900e+04_kr,7.200000e+04_kr,7.499000e+04_kr,&
     7.950000e+04_kr,8.229700e+04_kr,8.250000e+04_kr,8.651700e+04_kr,&
     9.803700e+04_kr,1.110900e+05_kr,1.167900e+05_kr,1.227700e+05_kr,&
     1.290700e+05_kr,1.356900e+05_kr,1.426400e+05_kr,1.499600e+05_kr,&
     1.576400e+05_kr,1.657300e+05_kr,1.742200e+05_kr,1.831600e+05_kr,&
     1.925500e+05_kr,2.024200e+05_kr,2.128000e+05_kr,2.237100e+05_kr,&
     2.351800e+05_kr,2.472400e+05_kr,2.732400e+05_kr,2.872500e+05_kr,&
     2.945200e+05_kr,2.972000e+05_kr,2.985000e+05_kr,3.019700e+05_kr,&
     3.337300e+05_kr,3.688300e+05_kr,3.877400e+05_kr,4.076200e+05_kr,&
     4.504900e+05_kr,5.234000e+05_kr,5.502300e+05_kr,5.784400e+05_kr,&
     6.081000e+05_kr,6.392800e+05_kr,6.720600e+05_kr,7.065100e+05_kr,&
     7.427400e+05_kr,7.808200e+05_kr,8.208500e+05_kr,8.629400e+05_kr,&
     9.071800e+05_kr,9.616400e+05_kr,1.002600e+06_kr,1.108000e+06_kr,&
     1.164800e+06_kr,1.224600e+06_kr,1.287300e+06_kr,1.353400e+06_kr,&
     1.422700e+06_kr,1.495700e+06_kr,1.572400e+06_kr,1.653000e+06_kr,&
     1.737700e+06_kr,1.826800e+06_kr,1.920500e+06_kr,2.019000e+06_kr,&
     2.122500e+06_kr,2.231300e+06_kr,2.306900e+06_kr,2.345700e+06_kr,&
     2.365300e+06_kr,2.385200e+06_kr,2.466000e+06_kr,2.592400e+06_kr,&
     2.725300e+06_kr,2.865000e+06_kr,3.011900e+06_kr,3.166400e+06_kr,&
     3.328700e+06_kr,3.678800e+06_kr,4.065700e+06_kr,4.493300e+06_kr,&
     4.723700e+06_kr,4.965900e+06_kr,5.220500e+06_kr,5.488100e+06_kr,&
     5.769500e+06_kr,6.065300e+06_kr,6.376300e+06_kr,6.592400e+06_kr,&
     6.703200e+06_kr,7.046900e+06_kr,7.408200e+06_kr,7.788000e+06_kr,&
     8.187300e+06_kr,8.607100e+06_kr,9.048400e+06_kr,9.512300e+06_kr,&
     1.000000e+07_kr,1.051300e+07_kr,1.105200e+07_kr,1.161800e+07_kr,&
     1.221400e+07_kr,1.284000e+07_kr,1.349900e+07_kr,1.384000e+07_kr,&
     1.419100e+07_kr,1.455000e+07_kr,1.491800e+07_kr,1.568300e+07_kr,&
     1.648700e+07_kr,1.690500e+07_kr,1.733300e+07_kr,1.964000e+07_kr/)
   real(kr),parameter::eg22(173)=(/&
     1.000010e-05_kr,3.000000e-03_kr,5.000000e-03_kr,6.900000e-03_kr,&
     1.000000e-02_kr,1.500000e-02_kr,2.000000e-02_kr,2.500000e-02_kr,&
     3.000000e-02_kr,3.500000e-02_kr,4.200000e-02_kr,5.000000e-02_kr,&
     5.800000e-02_kr,6.700000e-02_kr,7.700000e-02_kr,8.000000e-02_kr,&
     9.500000e-02_kr,1.000000e-01_kr,1.150000e-01_kr,1.340000e-01_kr,&
     1.400000e-01_kr,1.600000e-01_kr,1.800000e-01_kr,1.890000e-01_kr,&
     2.200000e-01_kr,2.480000e-01_kr,2.800000e-01_kr,3.000000e-01_kr,&
     3.145000e-01_kr,3.200000e-01_kr,3.500000e-01_kr,3.910000e-01_kr,&
     4.000000e-01_kr,4.330000e-01_kr,4.850000e-01_kr,5.000000e-01_kr,&
     5.400000e-01_kr,6.250000e-01_kr,7.050000e-01_kr,7.800000e-01_kr,&
     7.900000e-01_kr,8.500000e-01_kr,8.600000e-01_kr,9.100000e-01_kr,&
     9.300000e-01_kr,9.500000e-01_kr,9.720000e-01_kr,9.860000e-01_kr,&
     9.960000e-01_kr,1.020000e+00_kr,1.035000e+00_kr,1.045000e+00_kr,&
     1.071000e+00_kr,1.097000e+00_kr,1.110000e+00_kr,1.123000e+00_kr,&
     1.150000e+00_kr,1.170000e+00_kr,1.235000e+00_kr,1.300000e+00_kr,&
     1.337500e+00_kr,1.370000e+00_kr,1.440000e+00_kr,1.475000e+00_kr,&
     1.500000e+00_kr,1.590000e+00_kr,1.670000e+00_kr,1.755000e+00_kr,&
     1.840000e+00_kr,1.930000e+00_kr,2.020000e+00_kr,2.100000e+00_kr,&
     2.130000e+00_kr,2.360000e+00_kr,2.550000e+00_kr,2.600000e+00_kr,&
     2.720000e+00_kr,2.767920e+00_kr,3.300000e+00_kr,3.380750e+00_kr,&
     4.000000e+00_kr,4.129250e+00_kr,5.043477e+00_kr,5.346430e+00_kr,&
     6.160116e+00_kr,7.523983e+00_kr,8.315287e+00_kr,9.189814e+00_kr,&
     9.905554e+00_kr,1.122446e+01_kr,1.370959e+01_kr,1.592827e+01_kr,&
     1.945484e+01_kr,2.260329e+01_kr,2.498050e+01_kr,2.760773e+01_kr,&
     3.051126e+01_kr,3.372015e+01_kr,3.726653e+01_kr,4.016900e+01_kr,&
     4.551744e+01_kr,4.825160e+01_kr,5.157802e+01_kr,5.559513e+01_kr,&
     6.790405e+01_kr,7.567357e+01_kr,9.166088e+01_kr,1.367420e+02_kr,&
     1.486254e+02_kr,2.039950e+02_kr,3.043248e+02_kr,3.717032e+02_kr,&
     4.539993e+02_kr,6.772874e+02_kr,7.485183e+02_kr,9.142423e+02_kr,&
     1.010394e+03_kr,1.234098e+03_kr,1.433817e+03_kr,1.507331e+03_kr,&
     2.034684e+03_kr,2.248673e+03_kr,3.354626e+03_kr,3.526622e+03_kr,&
     5.004514e+03_kr,5.530844e+03_kr,7.465858e+03_kr,9.118820e+03_kr,&
     1.113775e+04_kr,1.503439e+04_kr,1.661557e+04_kr,2.478752e+04_kr,&
     2.739445e+04_kr,2.928300e+04_kr,3.697864e+04_kr,4.086771e+04_kr,&
     5.516564e+04_kr,6.737947e+04_kr,8.229747e+04_kr,1.110900e+05_kr,&
     1.227734e+05_kr,1.831564e+05_kr,2.472353e+05_kr,2.732372e+05_kr,&
     3.019738e+05_kr,4.076220e+05_kr,4.504920e+05_kr,4.978707e+05_kr,&
     5.502322e+05_kr,6.081006e+05_kr,8.208500e+05_kr,9.071795e+05_kr,&
     1.002588e+06_kr,1.108032e+06_kr,1.224564e+06_kr,1.353353e+06_kr,&
     1.652989e+06_kr,2.018965e+06_kr,2.231302e+06_kr,2.465970e+06_kr,&
     3.011942e+06_kr,3.678794e+06_kr,4.493290e+06_kr,5.488116e+06_kr,&
     6.065307e+06_kr,6.703200e+06_kr,8.187308e+06_kr,1.000000e+07_kr,&
     1.1618343e+07_kr,1.3840307e+07_kr,1.4918247e+07_kr,1.733253e+07_kr,&
     1.964033e+07_kr/)
   real(kr),parameter::eg23(176)=(/&
     1.000010e-05_kr,1.000010e-01_kr,4.139940e-01_kr,5.315790e-01_kr,&
     6.825600e-01_kr,8.764250e-01_kr,1.123000e+00_kr,1.440000e+00_kr,&
     1.855390e+00_kr,2.382370e+00_kr,3.059020e+00_kr,3.927860e+00_kr,&
     5.043480e+00_kr,6.475950e+00_kr,8.315290e+00_kr,1.067700e+01_kr,&
     1.370960e+01_kr,1.760350e+01_kr,2.260330e+01_kr,2.902320e+01_kr,&
     3.726650e+01_kr,4.785120e+01_kr,6.144210e+01_kr,7.889320e+01_kr,&
     1.013010e+02_kr,1.300730e+02_kr,1.670170e+02_kr,2.144540e+02_kr,&
     2.753640e+02_kr,3.535750e+02_kr,4.539990e+02_kr,5.829470e+02_kr,&
     7.485180e+02_kr,9.611170e+02_kr,1.234100e+03_kr,1.584610e+03_kr,&
     2.034680e+03_kr,2.248670e+03_kr,2.485170e+03_kr,2.612590e+03_kr,&
     2.746540e+03_kr,3.035390e+03_kr,3.354630e+03_kr,3.707440e+03_kr,&
     4.307420e+03_kr,5.530840e+03_kr,7.101740e+03_kr,9.118820e+03_kr,&
     1.059460e+04_kr,1.170880e+04_kr,1.503440e+04_kr,1.930450e+04_kr,&
     2.187490e+04_kr,2.357860e+04_kr,2.417550e+04_kr,2.478750e+04_kr,&
     2.605840e+04_kr,2.700010e+04_kr,2.850110e+04_kr,3.182780e+04_kr,&
     3.430670e+04_kr,4.086770e+04_kr,4.630920e+04_kr,5.247520e+04_kr,&
     5.656220e+04_kr,6.737950e+04_kr,7.202450e+04_kr,7.949870e+04_kr,&
     8.250340e+04_kr,8.651700e+04_kr,9.803650e+04_kr,1.110900e+05_kr,&
     1.167860e+05_kr,1.227730e+05_kr,1.290680e+05_kr,1.356860e+05_kr,&
     1.426420e+05_kr,1.499560e+05_kr,1.576440e+05_kr,1.657270e+05_kr,&
     1.742240e+05_kr,1.831560e+05_kr,1.925470e+05_kr,2.024190e+05_kr,&
     2.127970e+05_kr,2.237080e+05_kr,2.351770e+05_kr,2.472350e+05_kr,&
     2.732370e+05_kr,2.872460e+05_kr,2.945180e+05_kr,2.972110e+05_kr,&
     2.984910e+05_kr,3.019740e+05_kr,3.337330e+05_kr,3.688320e+05_kr,&
     3.877420e+05_kr,4.076220e+05_kr,4.504920e+05_kr,4.978710e+05_kr,&
     5.233970e+05_kr,5.502320e+05_kr,5.784430e+05_kr,6.081010e+05_kr,&
     6.392790e+05_kr,6.720550e+05_kr,7.065120e+05_kr,7.427360e+05_kr,&
     7.808170e+05_kr,8.208500e+05_kr,8.629360e+05_kr,9.071800e+05_kr,&
     9.616720e+05_kr,1.002590e+06_kr,1.108030e+06_kr,1.164840e+06_kr,&
     1.224560e+06_kr,1.287350e+06_kr,1.353350e+06_kr,1.422740e+06_kr,&
     1.495690e+06_kr,1.572370e+06_kr,1.652990e+06_kr,1.737740e+06_kr,&
     1.826840e+06_kr,1.920500e+06_kr,2.018970e+06_kr,2.122480e+06_kr,&
     2.231300e+06_kr,2.306930e+06_kr,2.345700e+06_kr,2.365330e+06_kr,&
     2.385130e+06_kr,2.465970e+06_kr,2.592400e+06_kr,2.725320e+06_kr,&
     2.865050e+06_kr,3.011940e+06_kr,3.166370e+06_kr,3.328710e+06_kr,&
     3.678790e+06_kr,4.065700e+06_kr,4.493290e+06_kr,4.723670e+06_kr,&
     4.965850e+06_kr,5.220460e+06_kr,5.488120e+06_kr,5.769500e+06_kr,&
     6.065310e+06_kr,6.376280e+06_kr,6.592410e+06_kr,6.703200e+06_kr,&
     7.046880e+06_kr,7.408180e+06_kr,7.788010e+06_kr,8.187310e+06_kr,&
     8.607080e+06_kr,9.048370e+06_kr,9.512290e+06_kr,1.000000e+07_kr,&
     1.051270e+07_kr,1.105170e+07_kr,1.161830e+07_kr,1.221400e+07_kr,&
     1.252320e+07_kr,1.284030e+07_kr,1.349860e+07_kr,1.384030e+07_kr,&
     1.419070e+07_kr,1.454990e+07_kr,1.491820e+07_kr,1.568310e+07_kr,&
     1.648720e+07_kr,1.690460e+07_kr,1.733250e+07_kr,1.964030e+07_kr/)
   real(kr),parameter::eg24(282)=(/&
     1.100027E-04_kr,2.499897E-03_kr,4.556021E-03_kr,7.145263E-03_kr,&
     1.045050E-02_kr,1.482996E-02_kr,2.001035E-02_kr,2.493942E-02_kr,&
     2.929889E-02_kr,3.439976E-02_kr,4.029993E-02_kr,4.730186E-02_kr,&
     5.549815E-02_kr,6.519936E-02_kr,7.649686E-02_kr,8.979683E-02_kr,&
     1.042977E-01_kr,1.199949E-01_kr,1.379994E-01_kr,1.618953E-01_kr,&
     1.900049E-01_kr,2.096102E-01_kr,2.311923E-01_kr,2.549965E-01_kr,&
     2.799888E-01_kr,3.050115E-01_kr,3.250079E-01_kr,3.529935E-01_kr,&
     3.900011E-01_kr,4.315786E-01_kr,4.750165E-01_kr,5.200108E-01_kr,&
     5.549897E-01_kr,5.949930E-01_kr,6.249987E-01_kr,7.199989E-01_kr,&
     8.200371E-01_kr,8.800244E-01_kr,9.199779E-01_kr,9.440222E-01_kr,&
     9.639598E-01_kr,9.819591E-01_kr,9.965005E-01_kr,1.009035E+00_kr,&
     1.021012E+00_kr,1.034993E+00_kr,1.077986E+00_kr,1.091982E+00_kr,&
     1.103950E+00_kr,1.116049E+00_kr,1.129974E+00_kr,1.147969E+00_kr,&
     1.169989E+00_kr,1.213968E+00_kr,1.250939E+00_kr,1.293038E+00_kr,&
     1.330952E+00_kr,1.380981E+00_kr,1.410007E+00_kr,1.443967E+00_kr,&
     1.519976E+00_kr,1.588030E+00_kr,1.668949E+00_kr,1.779966E+00_kr,&
     1.900077E+00_kr,1.989920E+00_kr,2.070095E+00_kr,2.156948E+00_kr,&
     2.217087E+00_kr,2.272986E+00_kr,2.330061E+00_kr,2.469941E+00_kr,&
     2.550003E+00_kr,2.590094E+00_kr,2.620053E+00_kr,2.640041E+00_kr,&
     2.700115E+00_kr,2.719898E+00_kr,2.740922E+00_kr,2.775121E+00_kr,&
     2.884047E+00_kr,3.142109E+00_kr,3.543073E+00_kr,3.712087E+00_kr,&
     3.882170E+00_kr,4.000000E+00_kr,4.219828E+00_kr,4.309812E+00_kr,&
     4.419800E+00_kr,4.767845E+00_kr,4.933232E+00_kr,5.109974E+00_kr,&
     5.210076E+00_kr,5.320112E+00_kr,5.380032E+00_kr,5.410245E+00_kr,&
     5.488167E+00_kr,5.530036E+00_kr,5.619790E+00_kr,5.720146E+00_kr,&
     5.800211E+00_kr,5.960142E+00_kr,6.059906E+00_kr,6.160108E+00_kr,&
     6.280153E+00_kr,6.359784E+00_kr,6.432057E+00_kr,6.481775E+00_kr,&
     6.514916E+00_kr,6.539066E+00_kr,6.556090E+00_kr,6.571843E+00_kr,&
     6.588293E+00_kr,6.606106E+00_kr,6.631257E+00_kr,6.716683E+00_kr,&
     6.742254E+00_kr,6.759807E+00_kr,6.776050E+00_kr,6.791653E+00_kr,&
     6.810696E+00_kr,6.835259E+00_kr,6.870208E+00_kr,6.917776E+00_kr,&
     6.994292E+00_kr,7.139869E+00_kr,7.380153E+00_kr,7.600350E+00_kr,&
     7.739943E+00_kr,7.839651E+00_kr,7.970079E+00_kr,8.130272E+00_kr,&
     8.300322E+00_kr,8.524074E+00_kr,8.673690E+00_kr,8.800375E+00_kr,&
     8.979950E+00_kr,9.140311E+00_kr,9.500024E+00_kr,1.057925E+01_kr,&
     1.080376E+01_kr,1.105292E+01_kr,1.126944E+01_kr,1.158944E+01_kr,&
     1.170943E+01_kr,1.181529E+01_kr,1.197947E+01_kr,1.213015E+01_kr,&
     1.230855E+01_kr,1.247210E+01_kr,1.259997E+01_kr,1.332970E+01_kr,&
     1.354604E+01_kr,1.404961E+01_kr,1.425053E+01_kr,1.447024E+01_kr,&
     1.459522E+01_kr,1.473012E+01_kr,1.486626E+01_kr,1.577923E+01_kr,&
     1.604977E+01_kr,1.655014E+01_kr,1.683053E+01_kr,1.744572E+01_kr,&
     1.756476E+01_kr,1.775903E+01_kr,1.795905E+01_kr,1.908484E+01_kr,&
     1.919969E+01_kr,1.939265E+01_kr,1.959735E+01_kr,2.007338E+01_kr,&
     2.027512E+01_kr,2.041754E+01_kr,2.051988E+01_kr,2.060213E+01_kr,&
     2.068470E+01_kr,2.076761E+01_kr,2.097632E+01_kr,2.106040E+01_kr,&
     2.114481E+01_kr,2.122956E+01_kr,2.133597E+01_kr,2.148585E+01_kr,&
     2.170178E+01_kr,2.200114E+01_kr,2.215569E+01_kr,2.237836E+01_kr,&
     2.253556E+01_kr,2.460856E+01_kr,2.760769E+01_kr,3.372011E+01_kr,&
     4.016895E+01_kr,4.399581E+01_kr,4.579131E+01_kr,5.267255E+01_kr,&
     6.144204E+01_kr,7.504548E+01_kr,8.895177E+01_kr,1.086459E+02_kr,&
     1.327005E+02_kr,1.620807E+02_kr,1.979658E+02_kr,2.417960E+02_kr,&
     2.837502E+02_kr,3.199275E+02_kr,3.535746E+02_kr,4.107950E+02_kr,&
     5.017462E+02_kr,6.128342E+02_kr,7.485173E+02_kr,9.075007E+02_kr,&
     1.064962E+03_kr,1.135007E+03_kr,1.345061E+03_kr,1.614038E+03_kr,&
     1.910451E+03_kr,2.219627E+03_kr,2.578838E+03_kr,2.996183E+03_kr,&
     3.481068E+03_kr,4.097345E+03_kr,5.004508E+03_kr,6.112520E+03_kr,&
     7.465848E+03_kr,9.118808E+03_kr,1.113774E+04_kr,1.360366E+04_kr,&
     1.489967E+04_kr,1.620045E+04_kr,1.858471E+04_kr,2.269941E+04_kr,&
     2.499908E+04_kr,2.610010E+04_kr,2.739441E+04_kr,2.928101E+04_kr,&
     3.345961E+04_kr,3.697859E+04_kr,4.086766E+04_kr,4.991587E+04_kr,&
     5.516557E+04_kr,6.737938E+04_kr,8.229736E+04_kr,9.466450E+04_kr,&
     1.156235E+05_kr,1.227732E+05_kr,1.399995E+05_kr,1.649989E+05_kr,&
     1.950077E+05_kr,2.300137E+05_kr,2.678264E+05_kr,3.206464E+05_kr,&
     3.838835E+05_kr,4.125012E+05_kr,4.560211E+05_kr,4.940018E+05_kr,&
     5.784425E+05_kr,7.065112E+05_kr,8.600058E+05_kr,9.511189E+05_kr,&
     1.051149E+06_kr,1.162048E+06_kr,1.286961E+06_kr,1.336941E+06_kr,&
     1.405768E+06_kr,1.636539E+06_kr,1.901387E+06_kr,2.231299E+06_kr,&
     2.725314E+06_kr,3.328707E+06_kr,4.065691E+06_kr,4.965847E+06_kr,&
     6.065299E+06_kr,6.703192E+06_kr,7.408173E+06_kr,8.187297E+06_kr,&
     9.048363E+06_kr,9.999987E+06_kr,1.161833E+07_kr,1.384029E+07_kr,&
     1.491823E+07_kr,1.964030E+07_kr/)
   real(kr),parameter::eg25(296)=(/&
     1.100027E-04_kr,2.499897E-03_kr,4.556021E-03_kr,7.145263E-03_kr,&
     1.045050E-02_kr,1.482996E-02_kr,2.001035E-02_kr,2.493942E-02_kr,&
     2.929889E-02_kr,3.439976E-02_kr,4.029993E-02_kr,4.730186E-02_kr,&
     5.549815E-02_kr,6.519936E-02_kr,7.649686E-02_kr,8.979683E-02_kr,&
     1.042977E-01_kr,1.199949E-01_kr,1.379994E-01_kr,1.618953E-01_kr,&
     1.900049E-01_kr,2.096102E-01_kr,2.311923E-01_kr,2.549965E-01_kr,&
     2.799888E-01_kr,3.050115E-01_kr,3.250079E-01_kr,3.529935E-01_kr,&
     3.900011E-01_kr,4.315786E-01_kr,4.750165E-01_kr,5.200108E-01_kr,&
     5.549897E-01_kr,5.949930E-01_kr,6.249987E-01_kr,7.199989E-01_kr,&
     8.200371E-01_kr,8.800244E-01_kr,9.199779E-01_kr,9.440222E-01_kr,&
     9.639598E-01_kr,9.819591E-01_kr,9.965005E-01_kr,1.009035E+00_kr,&
     1.021012E+00_kr,1.034993E+00_kr,1.077986E+00_kr,1.091982E+00_kr,&
     1.103950E+00_kr,1.116049E+00_kr,1.129974E+00_kr,1.147969E+00_kr,&
     1.169989E+00_kr,1.213968E+00_kr,1.250939E+00_kr,1.293038E+00_kr,&
     1.330952E+00_kr,1.380981E+00_kr,1.410007E+00_kr,1.443967E+00_kr,&
     1.519976E+00_kr,1.588030E+00_kr,1.668949E+00_kr,1.779966E+00_kr,&
     1.900077E+00_kr,1.989920E+00_kr,2.070095E+00_kr,2.156948E+00_kr,&
     2.217087E+00_kr,2.272986E+00_kr,2.330061E+00_kr,2.469941E+00_kr,&
     2.550003E+00_kr,2.590094E+00_kr,2.620053E+00_kr,2.640041E+00_kr,&
     2.700115E+00_kr,2.719898E+00_kr,2.740922E+00_kr,2.775121E+00_kr,&
     2.884047E+00_kr,3.142109E+00_kr,3.543073E+00_kr,3.712087E+00_kr,&
     3.882170E+00_kr,4.000000E+00_kr,4.219828E+00_kr,4.309812E+00_kr,&
     4.402156E+00_kr,4.632489E+00_kr,4.958456E+00_kr,5.317984E+00_kr,&
     5.618667E+00_kr,5.918567E+00_kr,6.228244E+00_kr,6.508405E+00_kr,&
     6.828427E+00_kr,7.214511E+00_kr,7.622423E+00_kr,7.965298E+00_kr,&
     8.415661E+00_kr,8.855993E+00_kr,9.310049E+00_kr,9.787385E+00_kr,&
     1.023788E+01_kr,1.070911E+01_kr,1.120202E+01_kr,1.174109E+01_kr,&
     1.223248E+01_kr,1.288539E+01_kr,1.357316E+01_kr,1.421211E+01_kr,&
     1.497069E+01_kr,1.573825E+01_kr,1.656173E+01_kr,1.737608E+01_kr,&
     1.830354E+01_kr,1.903148E+01_kr,2.010753E+01_kr,2.053425E+01_kr,&
     2.122319E+01_kr,2.200114E+01_kr,2.267118E+01_kr,2.465783E+01_kr,&
     2.788515E+01_kr,3.169295E+01_kr,3.308547E+01_kr,3.453918E+01_kr,&
     3.569799E+01_kr,3.605676E+01_kr,3.641914E+01_kr,3.685880E+01_kr,&
     3.730377E+01_kr,3.779188E+01_kr,3.878736E+01_kr,3.972951E+01_kr,&
     4.122704E+01_kr,4.214409E+01_kr,4.312463E+01_kr,4.417214E+01_kr,&
     4.529037E+01_kr,4.620529E+01_kr,4.751732E+01_kr,4.925911E+01_kr,&
     5.178468E+01_kr,5.298953E+01_kr,5.405999E+01_kr,5.705949E+01_kr,&
     5.992503E+01_kr,6.230828E+01_kr,6.363059E+01_kr,6.459225E+01_kr,&
     6.504598E+01_kr,6.550290E+01_kr,6.583123E+01_kr,6.616121E+01_kr,&
     6.649285E+01_kr,6.682614E+01_kr,6.906820E+01_kr,7.188692E+01_kr,&
     7.355948E+01_kr,7.633216E+01_kr,7.936793E+01_kr,8.393934E+01_kr,&
     8.877405E+01_kr,9.332559E+01_kr,9.732874E+01_kr,1.005942E+02_kr,&
     1.010984E+02_kr,1.016052E+02_kr,1.021145E+02_kr,1.030376E+02_kr,&
     1.056461E+02_kr,1.102879E+02_kr,1.128539E+02_kr,1.154797E+02_kr,&
     1.165237E+02_kr,1.175771E+02_kr,1.205536E+02_kr,1.262286E+02_kr,&
     1.327005E+02_kr,1.395042E+02_kr,1.466567E+02_kr,1.541759E+02_kr,&
     1.630561E+02_kr,1.675186E+02_kr,1.752291E+02_kr,1.832945E+02_kr,&
     1.849516E+02_kr,1.862508E+02_kr,1.875592E+02_kr,1.888767E+02_kr,&
     1.902035E+02_kr,1.930780E+02_kr,1.959960E+02_kr,2.009577E+02_kr,&
     2.121077E+02_kr,2.243247E+02_kr,2.355903E+02_kr,2.417960E+02_kr,&
     2.567478E+02_kr,2.682969E+02_kr,2.764678E+02_kr,2.848875E+02_kr,&
     2.883267E+02_kr,2.959215E+02_kr,3.199275E+02_kr,3.353230E+02_kr,&
     3.535746E+02_kr,3.717027E+02_kr,3.907603E+02_kr,4.190936E+02_kr,&
     4.539987E+02_kr,5.017462E+02_kr,5.392042E+02_kr,5.771455E+02_kr,&
     5.929407E+02_kr,6.000988E+02_kr,6.128342E+02_kr,6.468370E+02_kr,&
     6.772865E+02_kr,7.485173E+02_kr,8.322179E+02_kr,9.096813E+02_kr,&
     9.824941E+02_kr,1.064323E+03_kr,1.134667E+03_kr,1.343582E+03_kr,&
     1.586197E+03_kr,1.811833E+03_kr,2.084104E+03_kr,2.397290E+03_kr,&
     2.700236E+03_kr,2.996183E+03_kr,3.481068E+03_kr,4.097345E+03_kr,&
     5.004508E+03_kr,6.112520E+03_kr,7.465848E+03_kr,9.118808E+03_kr,&
     1.113774E+04_kr,1.360366E+04_kr,1.489967E+04_kr,1.620045E+04_kr,&
     1.858471E+04_kr,2.269941E+04_kr,2.499908E+04_kr,2.610010E+04_kr,&
     2.739441E+04_kr,2.928101E+04_kr,3.345961E+04_kr,3.697859E+04_kr,&
     4.086766E+04_kr,4.991587E+04_kr,5.516557E+04_kr,6.737938E+04_kr,&
     8.229736E+04_kr,9.466450E+04_kr,1.156235E+05_kr,1.227732E+05_kr,&
     1.400976E+05_kr,1.650650E+05_kr,1.950662E+05_kr,2.300598E+05_kr,&
     2.678264E+05_kr,3.206464E+05_kr,3.838835E+05_kr,4.125012E+05_kr,&
     4.560211E+05_kr,4.940018E+05_kr,5.784425E+05_kr,7.065112E+05_kr,&
     8.600058E+05_kr,9.511189E+05_kr,1.051149E+06_kr,1.162048E+06_kr,&
     1.286961E+06_kr,1.336941E+06_kr,1.405768E+06_kr,1.636539E+06_kr,&
     1.901387E+06_kr,2.231299E+06_kr,2.725314E+06_kr,3.328707E+06_kr,&
     4.065691E+06_kr,4.965847E+06_kr,6.065299E+06_kr,6.703192E+06_kr,&
     7.408173E+06_kr,8.187297E+06_kr,9.048363E+06_kr,9.999987E+06_kr,&
     1.161833E+07_kr,1.384029E+07_kr,1.491823E+07_kr,1.964030E+07_kr/)
   real(kr),parameter::eg26(362)=(/&
     1.100027E-04_kr,2.499897E-03_kr,4.556021E-03_kr,7.145263E-03_kr,&
     1.045050E-02_kr,1.482996E-02_kr,2.001035E-02_kr,2.493942E-02_kr,&
     2.929889E-02_kr,3.439976E-02_kr,4.029993E-02_kr,4.730186E-02_kr,&
     5.549815E-02_kr,6.519936E-02_kr,7.649686E-02_kr,8.979683E-02_kr,&
     1.042977E-01_kr,1.199949E-01_kr,1.379994E-01_kr,1.618953E-01_kr,&
     1.900049E-01_kr,2.096102E-01_kr,2.311923E-01_kr,2.549965E-01_kr,&
     2.799888E-01_kr,3.050115E-01_kr,3.250079E-01_kr,3.529935E-01_kr,&
     3.900011E-01_kr,4.315786E-01_kr,4.750165E-01_kr,5.200108E-01_kr,&
     5.549897E-01_kr,5.949930E-01_kr,6.249987E-01_kr,7.199989E-01_kr,&
     8.200371E-01_kr,8.800244E-01_kr,9.199779E-01_kr,9.440222E-01_kr,&
     9.639598E-01_kr,9.819591E-01_kr,9.965005E-01_kr,1.009035E+00_kr,&
     1.021012E+00_kr,1.034993E+00_kr,1.077986E+00_kr,1.091982E+00_kr,&
     1.103950E+00_kr,1.116049E+00_kr,1.129974E+00_kr,1.147969E+00_kr,&
     1.169989E+00_kr,1.213968E+00_kr,1.250939E+00_kr,1.293038E+00_kr,&
     1.330952E+00_kr,1.380981E+00_kr,1.410007E+00_kr,1.443967E+00_kr,&
     1.519976E+00_kr,1.588030E+00_kr,1.668949E+00_kr,1.779966E+00_kr,&
     1.900077E+00_kr,1.989920E+00_kr,2.070095E+00_kr,2.156948E+00_kr,&
     2.217087E+00_kr,2.272986E+00_kr,2.330061E+00_kr,2.469941E+00_kr,&
     2.550003E+00_kr,2.590094E+00_kr,2.620053E+00_kr,2.640041E+00_kr,&
     2.700115E+00_kr,2.719898E+00_kr,2.740922E+00_kr,2.775121E+00_kr,&
     2.884047E+00_kr,3.142109E+00_kr,3.543073E+00_kr,3.712087E+00_kr,&
     3.882170E+00_kr,4.000000E+00_kr,4.219828E+00_kr,4.309812E+00_kr,&
     4.419800E+00_kr,4.767845E+00_kr,4.933232E+00_kr,5.109974E+00_kr,&
     5.210076E+00_kr,5.320112E+00_kr,5.380032E+00_kr,5.410245E+00_kr,&
     5.488167E+00_kr,5.530036E+00_kr,5.619790E+00_kr,5.720146E+00_kr,&
     5.800211E+00_kr,5.960142E+00_kr,6.059906E+00_kr,6.160108E+00_kr,&
     6.280153E+00_kr,6.359784E+00_kr,6.432057E+00_kr,6.481775E+00_kr,&
     6.514916E+00_kr,6.539066E+00_kr,6.556090E+00_kr,6.571843E+00_kr,&
     6.588293E+00_kr,6.606106E+00_kr,6.631257E+00_kr,6.716683E+00_kr,&
     6.742254E+00_kr,6.759807E+00_kr,6.776050E+00_kr,6.791653E+00_kr,&
     6.810696E+00_kr,6.835259E+00_kr,6.870208E+00_kr,6.917776E+00_kr,&
     6.994292E+00_kr,7.139869E+00_kr,7.380153E+00_kr,7.600350E+00_kr,&
     7.739943E+00_kr,7.839651E+00_kr,7.970079E+00_kr,8.130272E+00_kr,&
     8.300322E+00_kr,8.524074E+00_kr,8.673690E+00_kr,8.800375E+00_kr,&
     8.979950E+00_kr,9.140311E+00_kr,9.500024E+00_kr,1.057925E+01_kr,&
     1.080376E+01_kr,1.105292E+01_kr,1.126944E+01_kr,1.158944E+01_kr,&
     1.170943E+01_kr,1.181529E+01_kr,1.197947E+01_kr,1.213015E+01_kr,&
     1.230855E+01_kr,1.247210E+01_kr,1.259997E+01_kr,1.332970E+01_kr,&
     1.354604E+01_kr,1.404961E+01_kr,1.425053E+01_kr,1.447024E+01_kr,&
     1.459522E+01_kr,1.473012E+01_kr,1.486626E+01_kr,1.577923E+01_kr,&
     1.604977E+01_kr,1.655014E+01_kr,1.683053E+01_kr,1.744572E+01_kr,&
     1.756476E+01_kr,1.775903E+01_kr,1.795905E+01_kr,1.908484E+01_kr,&
     1.919969E+01_kr,1.939265E+01_kr,1.959735E+01_kr,2.007338E+01_kr,&
     2.027512E+01_kr,2.041754E+01_kr,2.051988E+01_kr,2.060213E+01_kr,&
     2.068470E+01_kr,2.076761E+01_kr,2.097632E+01_kr,2.106040E+01_kr,&
     2.114481E+01_kr,2.122956E+01_kr,2.133597E+01_kr,2.148585E+01_kr,&
     2.170178E+01_kr,2.200114E+01_kr,2.215569E+01_kr,2.237836E+01_kr,&
     2.253556E+01_kr,2.465783E+01_kr,2.788515E+01_kr,3.169295E+01_kr,&
     3.308547E+01_kr,3.453918E+01_kr,3.569799E+01_kr,3.605676E+01_kr,&
     3.641914E+01_kr,3.685880E+01_kr,3.730377E+01_kr,3.779188E+01_kr,&
     3.878736E+01_kr,3.972951E+01_kr,4.122704E+01_kr,4.214409E+01_kr,&
     4.312463E+01_kr,4.417214E+01_kr,4.529037E+01_kr,4.620529E+01_kr,&
     4.751732E+01_kr,4.925911E+01_kr,5.178468E+01_kr,5.298953E+01_kr,&
     5.405999E+01_kr,5.705949E+01_kr,5.992503E+01_kr,6.230828E+01_kr,&
     6.363059E+01_kr,6.459225E+01_kr,6.504598E+01_kr,6.550290E+01_kr,&
     6.583123E+01_kr,6.616121E+01_kr,6.649285E+01_kr,6.682614E+01_kr,&
     6.906820E+01_kr,7.188692E+01_kr,7.355948E+01_kr,7.633216E+01_kr,&
     7.936793E+01_kr,8.393934E+01_kr,8.877405E+01_kr,9.332559E+01_kr,&
     9.732874E+01_kr,1.005942E+02_kr,1.010984E+02_kr,1.016052E+02_kr,&
     1.021145E+02_kr,1.030376E+02_kr,1.056461E+02_kr,1.102879E+02_kr,&
     1.128539E+02_kr,1.154797E+02_kr,1.165237E+02_kr,1.175771E+02_kr,&
     1.205536E+02_kr,1.262286E+02_kr,1.327005E+02_kr,1.395042E+02_kr,&
     1.466567E+02_kr,1.541759E+02_kr,1.630561E+02_kr,1.675186E+02_kr,&
     1.752291E+02_kr,1.832945E+02_kr,1.849516E+02_kr,1.862508E+02_kr,&
     1.875592E+02_kr,1.888767E+02_kr,1.902035E+02_kr,1.930780E+02_kr,&
     1.959960E+02_kr,2.009577E+02_kr,2.121077E+02_kr,2.243247E+02_kr,&
     2.355903E+02_kr,2.417960E+02_kr,2.567478E+02_kr,2.682969E+02_kr,&
     2.764678E+02_kr,2.848875E+02_kr,2.883267E+02_kr,2.959215E+02_kr,&
     3.199275E+02_kr,3.353230E+02_kr,3.535746E+02_kr,3.717027E+02_kr,&
     3.907603E+02_kr,4.190936E+02_kr,4.539987E+02_kr,5.017462E+02_kr,&
     5.392042E+02_kr,5.771455E+02_kr,5.929407E+02_kr,6.000988E+02_kr,&
     6.128342E+02_kr,6.468370E+02_kr,6.772865E+02_kr,7.485173E+02_kr,&
     8.322179E+02_kr,9.096813E+02_kr,9.824941E+02_kr,1.064323E+03_kr,&
     1.134667E+03_kr,1.343582E+03_kr,1.586197E+03_kr,1.811833E+03_kr,&
     2.084104E+03_kr,2.397290E+03_kr,2.700236E+03_kr,2.996183E+03_kr,&
     3.481068E+03_kr,4.097345E+03_kr,5.004508E+03_kr,6.112520E+03_kr,&
     7.465848E+03_kr,9.118808E+03_kr,1.113774E+04_kr,1.360366E+04_kr,&
     1.489967E+04_kr,1.620045E+04_kr,1.858471E+04_kr,2.269941E+04_kr,&
     2.499908E+04_kr,2.610010E+04_kr,2.739441E+04_kr,2.928101E+04_kr,&
     3.345961E+04_kr,3.697859E+04_kr,4.086766E+04_kr,4.991587E+04_kr,&
     5.516557E+04_kr,6.737938E+04_kr,8.229736E+04_kr,9.466450E+04_kr,&
     1.156235E+05_kr,1.227732E+05_kr,1.400976E+05_kr,1.650650E+05_kr,&
     1.950662E+05_kr,2.300598E+05_kr,2.678264E+05_kr,3.206464E+05_kr,&
     3.838835E+05_kr,4.125012E+05_kr,4.560211E+05_kr,4.940018E+05_kr,&
     5.784425E+05_kr,7.065112E+05_kr,8.600058E+05_kr,9.511189E+05_kr,&
     1.051149E+06_kr,1.162048E+06_kr,1.286961E+06_kr,1.336941E+06_kr,&
     1.405768E+06_kr,1.636539E+06_kr,1.901387E+06_kr,2.231299E+06_kr,&
     2.725314E+06_kr,3.328707E+06_kr,4.065691E+06_kr,4.965847E+06_kr,&
     6.065299E+06_kr,6.703192E+06_kr,7.408173E+06_kr,8.187297E+06_kr,&
     9.048363E+06_kr,9.999987E+06_kr,1.161833E+07_kr,1.384029E+07_kr,&
     1.491823E+07_kr,1.964030E+07_kr/)
   real(kr),parameter::eg27(316)=(/&
     1.100028E-04_kr,2.499900E-03_kr,4.556027E-03_kr,7.145272E-03_kr,&
     1.045051E-02_kr,1.482998E-02_kr,2.001038E-02_kr,2.493945E-02_kr,&
     2.929893E-02_kr,3.439981E-02_kr,4.029998E-02_kr,4.730192E-02_kr,&
     5.549822E-02_kr,6.519944E-02_kr,7.649695E-02_kr,8.979694E-02_kr,&
     1.000001E-01_kr,1.199951E-01_kr,1.379996E-01_kr,1.618955E-01_kr,&
     1.900051E-01_kr,2.096105E-01_kr,2.311926E-01_kr,2.549968E-01_kr,&
     2.799892E-01_kr,3.050119E-01_kr,3.250084E-01_kr,3.529940E-01_kr,&
     3.900016E-01_kr,4.315792E-01_kr,4.750172E-01_kr,5.200115E-01_kr,&
     5.400010E-01_kr,5.949938E-01_kr,6.249996E-01_kr,7.199998E-01_kr,&
     8.200382E-01_kr,8.800255E-01_kr,9.199791E-01_kr,9.440234E-01_kr,&
     9.639611E-01_kr,9.819603E-01_kr,9.965018E-01_kr,1.009036E+00_kr,&
     1.021013E+00_kr,1.034994E+00_kr,1.077987E+00_kr,1.091983E+00_kr,&
     1.103951E+00_kr,1.116050E+00_kr,1.129975E+00_kr,1.147971E+00_kr,&
     1.169991E+00_kr,1.213970E+00_kr,1.250940E+00_kr,1.293040E+00_kr,&
     1.330954E+00_kr,1.380983E+00_kr,1.410008E+00_kr,1.443969E+00_kr,&
     1.519978E+00_kr,1.588032E+00_kr,1.668952E+00_kr,1.779968E+00_kr,&
     1.900079E+00_kr,1.989922E+00_kr,2.070097E+00_kr,2.156951E+00_kr,&
     2.217090E+00_kr,2.272989E+00_kr,2.330064E+00_kr,2.469944E+00_kr,&
     2.550006E+00_kr,2.590098E+00_kr,2.620056E+00_kr,2.640044E+00_kr,&
     2.700118E+00_kr,2.719901E+00_kr,2.740926E+00_kr,2.775125E+00_kr,&
     2.884050E+00_kr,3.142113E+00_kr,3.543077E+00_kr,3.712092E+00_kr,&
     3.882175E+00_kr,4.000000E+00_kr,4.219833E+00_kr,4.309818E+00_kr,&
     4.402161E+00_kr,4.632495E+00_kr,4.958462E+00_kr,5.317991E+00_kr,&
     5.618674E+00_kr,5.918575E+00_kr,6.228252E+00_kr,6.508413E+00_kr,&
     6.828436E+00_kr,7.214520E+00_kr,7.622433E+00_kr,7.965308E+00_kr,&
     8.315287E+00_kr,8.856004E+00_kr,9.310061E+00_kr,9.787398E+00_kr,&
     1.023789E+01_kr,1.070912E+01_kr,1.120204E+01_kr,1.174110E+01_kr,&
     1.223249E+01_kr,1.288541E+01_kr,1.370959E+01_kr,1.421213E+01_kr,&
     1.497071E+01_kr,1.573828E+01_kr,1.656175E+01_kr,1.737610E+01_kr,&
     1.830356E+01_kr,1.903150E+01_kr,2.010756E+01_kr,2.053428E+01_kr,&
     2.122322E+01_kr,2.200117E+01_kr,2.260329E+01_kr,2.465786E+01_kr,&
     2.788519E+01_kr,3.169299E+01_kr,3.308552E+01_kr,3.453923E+01_kr,&
     3.569804E+01_kr,3.605681E+01_kr,3.641918E+01_kr,3.685885E+01_kr,&
     3.730382E+01_kr,3.779193E+01_kr,3.878741E+01_kr,4.016900E+01_kr,&
     4.122709E+01_kr,4.214414E+01_kr,4.312469E+01_kr,4.417220E+01_kr,&
     4.529043E+01_kr,4.620535E+01_kr,4.751739E+01_kr,4.925918E+01_kr,&
     5.178475E+01_kr,5.298960E+01_kr,5.406006E+01_kr,5.705956E+01_kr,&
     5.992511E+01_kr,6.230836E+01_kr,6.363067E+01_kr,6.459233E+01_kr,&
     6.504606E+01_kr,6.550298E+01_kr,6.583132E+01_kr,6.616130E+01_kr,&
     6.649293E+01_kr,6.682623E+01_kr,6.790405E+01_kr,6.906828E+01_kr,&
     7.188702E+01_kr,7.355958E+01_kr,7.633226E+01_kr,7.936803E+01_kr,&
     8.393945E+01_kr,8.877416E+01_kr,9.166088E+01_kr,9.332571E+01_kr,&
     9.732887E+01_kr,1.005943E+02_kr,1.010985E+02_kr,1.016053E+02_kr,&
     1.021146E+02_kr,1.030378E+02_kr,1.056462E+02_kr,1.102881E+02_kr,&
     1.128541E+02_kr,1.154798E+02_kr,1.165238E+02_kr,1.175773E+02_kr,&
     1.205537E+02_kr,1.262287E+02_kr,1.327006E+02_kr,1.395043E+02_kr,&
     1.485759E+02_kr,1.541761E+02_kr,1.630563E+02_kr,1.675188E+02_kr,&
     1.752293E+02_kr,1.832948E+02_kr,1.849519E+02_kr,1.862511E+02_kr,&
     1.875594E+02_kr,1.888769E+02_kr,1.902037E+02_kr,1.930783E+02_kr,&
     1.959963E+02_kr,2.009579E+02_kr,2.121080E+02_kr,2.243249E+02_kr,&
     2.355906E+02_kr,2.417963E+02_kr,2.567482E+02_kr,2.682973E+02_kr,&
     2.723521E+02_kr,2.764682E+02_kr,2.848879E+02_kr,2.883271E+02_kr,&
     2.959219E+02_kr,3.043248E+02_kr,3.199279E+02_kr,3.353235E+02_kr,&
     3.535750E+02_kr,3.717032E+02_kr,3.907608E+02_kr,4.190942E+02_kr,&
     4.539993E+02_kr,5.017468E+02_kr,5.392049E+02_kr,5.771462E+02_kr,&
     5.929414E+02_kr,6.000996E+02_kr,6.128350E+02_kr,6.468379E+02_kr,&
     6.772874E+02_kr,7.485183E+02_kr,8.322190E+02_kr,9.096825E+02_kr,&
     9.824954E+02_kr,1.064325E+03_kr,1.134668E+03_kr,1.234098E+03_kr,&
     1.343584E+03_kr,1.586199E+03_kr,1.811835E+03_kr,2.034684E+03_kr,&
     2.284941E+03_kr,2.423809E+03_kr,2.571117E+03_kr,2.768596E+03_kr,&
     2.951579E+03_kr,3.146656E+03_kr,3.354626E+03_kr,3.707435E+03_kr,&
     4.097350E+03_kr,4.528272E+03_kr,5.004514E+03_kr,5.530844E+03_kr,&
     6.112528E+03_kr,6.755388E+03_kr,7.465858E+03_kr,8.251049E+03_kr,&
     9.118820E+03_kr,1.007785E+04_kr,1.113775E+04_kr,1.360368E+04_kr,&
     1.503439E+04_kr,1.620047E+04_kr,1.858473E+04_kr,2.269944E+04_kr,&
     2.478752E+04_kr,2.610013E+04_kr,2.739445E+04_kr,2.928104E+04_kr,&
     3.345965E+04_kr,3.697864E+04_kr,4.086771E+04_kr,4.991594E+04_kr,&
     5.516564E+04_kr,6.737947E+04_kr,8.229747E+04_kr,9.466462E+04_kr,&
     1.110900E+05_kr,1.227734E+05_kr,1.400977E+05_kr,1.650652E+05_kr,&
     1.831564E+05_kr,1.950665E+05_kr,2.300601E+05_kr,2.678268E+05_kr,&
     3.019738E+05_kr,3.206469E+05_kr,3.838840E+05_kr,4.125017E+05_kr,&
     4.560217E+05_kr,4.978707E+05_kr,5.784432E+05_kr,7.065121E+05_kr,&
     7.808167E+05_kr,8.208500E+05_kr,9.508348E+05_kr,1.051150E+06_kr,&
     1.162049E+06_kr,1.286963E+06_kr,1.353353E+06_kr,1.408584E+06_kr,&
     1.636541E+06_kr,1.901390E+06_kr,2.231302E+06_kr,2.528396E+06_kr,&
     2.865048E+06_kr,3.246525E+06_kr,3.678794E+06_kr,4.168620E+06_kr,&
     4.723666E+06_kr,5.352614E+06_kr,6.065307E+06_kr,6.703200E+06_kr,&
     7.408182E+06_kr,8.187308E+06_kr,9.048374E+06_kr,1.000000E+07_kr,&
     1.161834E+07_kr,1.384031E+07_kr,1.491825E+07_kr,1.964033E+07_kr/)
   real(kr),parameter::eg28(90)=(/&
     2.000000E-04_kr,5.000000E-03_kr,1.000000E-02_kr,1.500000E-02_kr,&
     2.000000E-02_kr,2.500000E-02_kr,3.000000E-02_kr,3.500000E-02_kr,&
     4.200000E-02_kr,5.000000E-02_kr,5.800000E-02_kr,6.700000E-02_kr,&
     8.000000E-02_kr,1.000000E-01_kr,1.400000E-01_kr,1.800000E-01_kr,&
     2.200000E-01_kr,2.500000E-01_kr,2.800000E-01_kr,3.000000E-01_kr,&
     3.200000E-01_kr,3.500000E-01_kr,4.000000E-01_kr,5.000000E-01_kr,&
     6.250000E-01_kr,7.800000E-01_kr,8.500000E-01_kr,9.100000E-01_kr,&
     9.500000E-01_kr,9.720000E-01_kr,9.960000E-01_kr,1.020000E+00_kr,&
     1.045000E+00_kr,1.071000E+00_kr,1.097000E+00_kr,1.123000E+00_kr,&
     1.150000E+00_kr,1.300000E+00_kr,1.500000E+00_kr,2.100000E+00_kr,&
     2.600000E+00_kr,3.300000E+00_kr,4.000000E+00_kr,5.043500E+00_kr,&
     6.476000E+00_kr,8.315300E+00_kr,1.067700E+01_kr,1.371000E+01_kr,&
     1.760300E+01_kr,2.260300E+01_kr,2.902300E+01_kr,3.726700E+01_kr,&
     4.785100E+01_kr,6.144200E+01_kr,7.889300E+01_kr,1.013000E+02_kr,&
     1.300700E+02_kr,1.670200E+02_kr,2.753600E+02_kr,4.540000E+02_kr,&
     7.485200E+02_kr,1.234100E+03_kr,2.034700E+03_kr,3.354600E+03_kr,&
     5.530800E+03_kr,9.118800E+03_kr,1.503400E+04_kr,2.478800E+04_kr,&
     4.086800E+04_kr,6.737900E+04_kr,8.651700E+04_kr,1.110900E+05_kr,&
     1.426400E+05_kr,1.831600E+05_kr,2.351800E+05_kr,3.019700E+05_kr,&
     3.877300E+05_kr,4.978700E+05_kr,6.392800E+05_kr,8.208500E+05_kr,&
     1.054000E+06_kr,1.353400E+06_kr,1.737700E+06_kr,2.231300E+06_kr,&
     2.865000E+06_kr,3.678800E+06_kr,4.723700E+06_kr,6.065300E+06_kr,&
     7.788000E+06_kr,1.000000E+07_kr/)
   real(kr),dimension(661),parameter::eg29=(/&
     1.000000E-05_kr,1.047129E-05_kr,1.096478E-05_kr,1.148154E-05_kr,&
     1.202264E-05_kr,1.258925E-05_kr,1.318257E-05_kr,1.380384E-05_kr,&
     1.445440E-05_kr,1.513561E-05_kr,1.584893E-05_kr,1.659587E-05_kr,&
     1.737801E-05_kr,1.819701E-05_kr,1.905461E-05_kr,1.995262E-05_kr,&
     2.089296E-05_kr,2.187762E-05_kr,2.290868E-05_kr,2.398833E-05_kr,&
     2.511886E-05_kr,2.630268E-05_kr,2.754229E-05_kr,2.884032E-05_kr,&
     3.019952E-05_kr,3.162278E-05_kr,3.311311E-05_kr,3.467369E-05_kr,&
     3.630781E-05_kr,3.801894E-05_kr,3.981072E-05_kr,4.168694E-05_kr,&
     4.365158E-05_kr,4.570882E-05_kr,4.786301E-05_kr,5.011872E-05_kr,&
     5.248075E-05_kr,5.495409E-05_kr,5.754399E-05_kr,6.025596E-05_kr,&
     6.309573E-05_kr,6.606934E-05_kr,6.918310E-05_kr,7.244360E-05_kr,&
     7.585776E-05_kr,7.943282E-05_kr,8.317638E-05_kr,8.709636E-05_kr,&
     9.120108E-05_kr,9.549926E-05_kr,1.000000E-04_kr,1.047129E-04_kr,&
     1.096478E-04_kr,1.148154E-04_kr,1.202264E-04_kr,1.258925E-04_kr,&
     1.318257E-04_kr,1.380384E-04_kr,1.445440E-04_kr,1.513561E-04_kr,&
     1.584893E-04_kr,1.659587E-04_kr,1.737801E-04_kr,1.819701E-04_kr,&
     1.905461E-04_kr,1.995262E-04_kr,2.089296E-04_kr,2.187762E-04_kr,&
     2.290868E-04_kr,2.398833E-04_kr,2.511886E-04_kr,2.630268E-04_kr,&
     2.754229E-04_kr,2.884032E-04_kr,3.019952E-04_kr,3.162278E-04_kr,&
     3.311311E-04_kr,3.467369E-04_kr,3.630781E-04_kr,3.801894E-04_kr,&
     3.981072E-04_kr,4.168694E-04_kr,4.365158E-04_kr,4.570882E-04_kr,&
     4.786301E-04_kr,5.011872E-04_kr,5.248075E-04_kr,5.495409E-04_kr,&
     5.754399E-04_kr,6.025596E-04_kr,6.309573E-04_kr,6.606934E-04_kr,&
     6.918310E-04_kr,7.244360E-04_kr,7.585776E-04_kr,7.943282E-04_kr,&
     8.317638E-04_kr,8.709636E-04_kr,9.120108E-04_kr,9.549926E-04_kr,&
     1.000000E-03_kr,1.047129E-03_kr,1.096478E-03_kr,1.148154E-03_kr,&
     1.202264E-03_kr,1.258925E-03_kr,1.318257E-03_kr,1.380384E-03_kr,&
     1.445440E-03_kr,1.513561E-03_kr,1.584893E-03_kr,1.659587E-03_kr,&
     1.737801E-03_kr,1.819701E-03_kr,1.905461E-03_kr,1.995262E-03_kr,&
     2.089296E-03_kr,2.187762E-03_kr,2.290868E-03_kr,2.398833E-03_kr,&
     2.511886E-03_kr,2.630268E-03_kr,2.754229E-03_kr,2.884032E-03_kr,&
     3.019952E-03_kr,3.162278E-03_kr,3.311311E-03_kr,3.467369E-03_kr,&
     3.630781E-03_kr,3.801894E-03_kr,3.981072E-03_kr,4.168694E-03_kr,&
     4.365158E-03_kr,4.570882E-03_kr,4.786301E-03_kr,5.011872E-03_kr,&
     5.248075E-03_kr,5.495409E-03_kr,5.754399E-03_kr,6.025596E-03_kr,&
     6.309573E-03_kr,6.606934E-03_kr,6.918310E-03_kr,7.244360E-03_kr,&
     7.585776E-03_kr,7.943282E-03_kr,8.317638E-03_kr,8.709636E-03_kr,&
     9.120108E-03_kr,9.549926E-03_kr,1.000000E-02_kr,1.047129E-02_kr,&
     1.096478E-02_kr,1.148154E-02_kr,1.202264E-02_kr,1.258925E-02_kr,&
     1.318257E-02_kr,1.380384E-02_kr,1.445440E-02_kr,1.513561E-02_kr,&
     1.584893E-02_kr,1.659587E-02_kr,1.737801E-02_kr,1.819701E-02_kr,&
     1.905461E-02_kr,1.995262E-02_kr,2.089296E-02_kr,2.187762E-02_kr,&
     2.290868E-02_kr,2.398833E-02_kr,2.511886E-02_kr,2.630268E-02_kr,&
     2.754229E-02_kr,2.884032E-02_kr,3.019952E-02_kr,3.162278E-02_kr,&
     3.311311E-02_kr,3.467369E-02_kr,3.630781E-02_kr,3.801894E-02_kr,&
     3.981072E-02_kr,4.168694E-02_kr,4.365158E-02_kr,4.570882E-02_kr,&
     4.786301E-02_kr,5.011872E-02_kr,5.248075E-02_kr,5.495409E-02_kr,&
     5.754399E-02_kr,6.025596E-02_kr,6.309573E-02_kr,6.606934E-02_kr,&
     6.918310E-02_kr,7.244360E-02_kr,7.585776E-02_kr,7.943282E-02_kr,&
     8.317638E-02_kr,8.709636E-02_kr,9.120108E-02_kr,9.549926E-02_kr,&
     1.000000E-01_kr,1.047129E-01_kr,1.096478E-01_kr,1.148154E-01_kr,&
     1.202264E-01_kr,1.258925E-01_kr,1.318257E-01_kr,1.380384E-01_kr,&
     1.445440E-01_kr,1.513561E-01_kr,1.584893E-01_kr,1.659587E-01_kr,&
     1.737801E-01_kr,1.819701E-01_kr,1.905461E-01_kr,1.995262E-01_kr,&
     2.089296E-01_kr,2.187762E-01_kr,2.290868E-01_kr,2.398833E-01_kr,&
     2.511886E-01_kr,2.630268E-01_kr,2.754229E-01_kr,2.884032E-01_kr,&
     3.019952E-01_kr,3.162278E-01_kr,3.311311E-01_kr,3.467369E-01_kr,&
     3.630781E-01_kr,3.801894E-01_kr,3.981072E-01_kr,4.168694E-01_kr,&
     4.365158E-01_kr,4.570882E-01_kr,4.786301E-01_kr,5.011872E-01_kr,&
     5.248075E-01_kr,5.495409E-01_kr,5.754399E-01_kr,6.025596E-01_kr,&
     6.309573E-01_kr,6.606934E-01_kr,6.918310E-01_kr,7.244360E-01_kr,&
     7.585776E-01_kr,7.943282E-01_kr,8.317638E-01_kr,8.709636E-01_kr,&
     9.120108E-01_kr,9.549926E-01_kr,1.000000E+00_kr,1.047129E+00_kr,&
     1.096478E+00_kr,1.148154E+00_kr,1.202264E+00_kr,1.258925E+00_kr,&
     1.318257E+00_kr,1.380384E+00_kr,1.445440E+00_kr,1.513561E+00_kr,&
     1.584893E+00_kr,1.659587E+00_kr,1.737801E+00_kr,1.819701E+00_kr,&
     1.905461E+00_kr,1.995262E+00_kr,2.089296E+00_kr,2.187762E+00_kr,&
     2.290868E+00_kr,2.398833E+00_kr,2.511886E+00_kr,2.630268E+00_kr,&
     2.754229E+00_kr,2.884032E+00_kr,3.019952E+00_kr,3.162278E+00_kr,&
     3.311311E+00_kr,3.467369E+00_kr,3.630781E+00_kr,3.801894E+00_kr,&
     3.981072E+00_kr,4.168694E+00_kr,4.365158E+00_kr,4.570882E+00_kr,&
     4.786301E+00_kr,5.011872E+00_kr,5.248075E+00_kr,5.495409E+00_kr,&
     5.754399E+00_kr,6.025596E+00_kr,6.309573E+00_kr,6.606934E+00_kr,&
     6.918310E+00_kr,7.244360E+00_kr,7.585776E+00_kr,7.943282E+00_kr,&
     8.317638E+00_kr,8.709636E+00_kr,9.120108E+00_kr,9.549926E+00_kr,&
     1.000000E+01_kr,1.047129E+01_kr,1.096478E+01_kr,1.148154E+01_kr,&
     1.202264E+01_kr,1.258925E+01_kr,1.318257E+01_kr,1.380384E+01_kr,&
     1.445440E+01_kr,1.513561E+01_kr,1.584893E+01_kr,1.659587E+01_kr,&
     1.737801E+01_kr,1.819701E+01_kr,1.905461E+01_kr,1.995262E+01_kr,&
     2.089296E+01_kr,2.187762E+01_kr,2.290868E+01_kr,2.398833E+01_kr,&
     2.511886E+01_kr,2.630268E+01_kr,2.754229E+01_kr,2.884032E+01_kr,&
     3.019952E+01_kr,3.162278E+01_kr,3.311311E+01_kr,3.467369E+01_kr,&
     3.630781E+01_kr,3.801894E+01_kr,3.981072E+01_kr,4.168694E+01_kr,&
     4.365158E+01_kr,4.570882E+01_kr,4.786301E+01_kr,5.011872E+01_kr,&
     5.248075E+01_kr,5.495409E+01_kr,5.754399E+01_kr,6.025596E+01_kr,&
     6.309573E+01_kr,6.606934E+01_kr,6.918310E+01_kr,7.244360E+01_kr,&
     7.585776E+01_kr,7.943282E+01_kr,8.317638E+01_kr,8.709636E+01_kr,&
     9.120108E+01_kr,9.549926E+01_kr,1.000000E+02_kr,1.047129E+02_kr,&
     1.096478E+02_kr,1.148154E+02_kr,1.202264E+02_kr,1.258925E+02_kr,&
     1.318257E+02_kr,1.380384E+02_kr,1.445440E+02_kr,1.513561E+02_kr,&
     1.584893E+02_kr,1.659587E+02_kr,1.737801E+02_kr,1.819701E+02_kr,&
     1.905461E+02_kr,1.995262E+02_kr,2.089296E+02_kr,2.187762E+02_kr,&
     2.290868E+02_kr,2.398833E+02_kr,2.511886E+02_kr,2.630268E+02_kr,&
     2.754229E+02_kr,2.884032E+02_kr,3.019952E+02_kr,3.162278E+02_kr,&
     3.311311E+02_kr,3.467369E+02_kr,3.630781E+02_kr,3.801894E+02_kr,&
     3.981072E+02_kr,4.168694E+02_kr,4.365158E+02_kr,4.570882E+02_kr,&
     4.786301E+02_kr,5.011872E+02_kr,5.248075E+02_kr,5.495409E+02_kr,&
     5.754399E+02_kr,6.025596E+02_kr,6.309573E+02_kr,6.606934E+02_kr,&
     6.918310E+02_kr,7.244360E+02_kr,7.585776E+02_kr,7.943282E+02_kr,&
     8.317638E+02_kr,8.709636E+02_kr,9.120108E+02_kr,9.549926E+02_kr,&
     1.000000E+03_kr,1.047129E+03_kr,1.096478E+03_kr,1.148154E+03_kr,&
     1.202264E+03_kr,1.258925E+03_kr,1.318257E+03_kr,1.380384E+03_kr,&
     1.445440E+03_kr,1.513561E+03_kr,1.584893E+03_kr,1.659587E+03_kr,&
     1.737801E+03_kr,1.819701E+03_kr,1.905461E+03_kr,1.995262E+03_kr,&
     2.089296E+03_kr,2.187762E+03_kr,2.290868E+03_kr,2.398833E+03_kr,&
     2.511886E+03_kr,2.630268E+03_kr,2.754229E+03_kr,2.884032E+03_kr,&
     3.019952E+03_kr,3.162278E+03_kr,3.311311E+03_kr,3.467369E+03_kr,&
     3.630781E+03_kr,3.801894E+03_kr,3.981072E+03_kr,4.168694E+03_kr,&
     4.365158E+03_kr,4.570882E+03_kr,4.786301E+03_kr,5.011872E+03_kr,&
     5.248075E+03_kr,5.495409E+03_kr,5.754399E+03_kr,6.025596E+03_kr,&
     6.309573E+03_kr,6.606934E+03_kr,6.918310E+03_kr,7.244360E+03_kr,&
     7.585776E+03_kr,7.943282E+03_kr,8.317638E+03_kr,8.709636E+03_kr,&
     9.120108E+03_kr,9.549926E+03_kr,1.000000E+04_kr,1.047129E+04_kr,&
     1.096478E+04_kr,1.148154E+04_kr,1.202264E+04_kr,1.258925E+04_kr,&
     1.318257E+04_kr,1.380384E+04_kr,1.445440E+04_kr,1.513561E+04_kr,&
     1.584893E+04_kr,1.659587E+04_kr,1.737801E+04_kr,1.819701E+04_kr,&
     1.905461E+04_kr,1.995262E+04_kr,2.089296E+04_kr,2.187762E+04_kr,&
     2.290868E+04_kr,2.398833E+04_kr,2.511886E+04_kr,2.630268E+04_kr,&
     2.754229E+04_kr,2.884032E+04_kr,3.019952E+04_kr,3.162278E+04_kr,&
     3.311311E+04_kr,3.467369E+04_kr,3.630781E+04_kr,3.801894E+04_kr,&
     3.981072E+04_kr,4.168694E+04_kr,4.365158E+04_kr,4.570882E+04_kr,&
     4.786301E+04_kr,5.011872E+04_kr,5.248075E+04_kr,5.495409E+04_kr,&
     5.754399E+04_kr,6.025596E+04_kr,6.309573E+04_kr,6.606934E+04_kr,&
     6.918310E+04_kr,7.244360E+04_kr,7.585776E+04_kr,7.943282E+04_kr,&
     8.317638E+04_kr,8.709636E+04_kr,9.120108E+04_kr,9.549926E+04_kr,&
     1.000000E+05_kr,1.047129E+05_kr,1.096478E+05_kr,1.148154E+05_kr,&
     1.202264E+05_kr,1.258925E+05_kr,1.318257E+05_kr,1.380384E+05_kr,&
     1.445440E+05_kr,1.513561E+05_kr,1.584893E+05_kr,1.659587E+05_kr,&
     1.737801E+05_kr,1.819701E+05_kr,1.905461E+05_kr,1.995262E+05_kr,&
     2.089296E+05_kr,2.187762E+05_kr,2.290868E+05_kr,2.398833E+05_kr,&
     2.511886E+05_kr,2.630268E+05_kr,2.754229E+05_kr,2.884032E+05_kr,&
     3.019952E+05_kr,3.162278E+05_kr,3.311311E+05_kr,3.467369E+05_kr,&
     3.630781E+05_kr,3.801894E+05_kr,3.981072E+05_kr,4.168694E+05_kr,&
     4.365158E+05_kr,4.570882E+05_kr,4.786301E+05_kr,5.011872E+05_kr,&
     5.248075E+05_kr,5.495409E+05_kr,5.754399E+05_kr,6.025596E+05_kr,&
     6.309573E+05_kr,6.606934E+05_kr,6.918310E+05_kr,7.244360E+05_kr,&
     7.585776E+05_kr,7.943282E+05_kr,8.317638E+05_kr,8.709636E+05_kr,&
     9.120108E+05_kr,9.549926E+05_kr,1.000000E+06_kr,1.047129E+06_kr,&
     1.096478E+06_kr,1.148154E+06_kr,1.202264E+06_kr,1.258925E+06_kr,&
     1.318257E+06_kr,1.380384E+06_kr,1.445440E+06_kr,1.513561E+06_kr,&
     1.584893E+06_kr,1.659587E+06_kr,1.737801E+06_kr,1.819701E+06_kr,&
     1.905461E+06_kr,1.995262E+06_kr,2.089296E+06_kr,2.187762E+06_kr,&
     2.290868E+06_kr,2.398833E+06_kr,2.511886E+06_kr,2.630268E+06_kr,&
     2.754229E+06_kr,2.884032E+06_kr,3.019952E+06_kr,3.162278E+06_kr,&
     3.311311E+06_kr,3.467369E+06_kr,3.630781E+06_kr,3.801894E+06_kr,&
     3.981072E+06_kr,4.168694E+06_kr,4.365158E+06_kr,4.570882E+06_kr,&
     4.786301E+06_kr,5.011872E+06_kr,5.248075E+06_kr,5.495409E+06_kr,&
     5.754399E+06_kr,6.025596E+06_kr,6.309573E+06_kr,6.606934E+06_kr,&
     6.918310E+06_kr,7.244360E+06_kr,7.585776E+06_kr,7.943282E+06_kr,&
     8.317638E+06_kr,8.709636E+06_kr,9.120108E+06_kr,9.549926E+06_kr,&
     1.000000E+07_kr,1.020000E+07_kr,1.040000E+07_kr,1.060000E+07_kr,&
     1.080000E+07_kr,1.100000E+07_kr,1.120000E+07_kr,1.140000E+07_kr,&
     1.160000E+07_kr,1.180000E+07_kr,1.200000E+07_kr,1.220000E+07_kr,&
     1.240000E+07_kr,1.260000E+07_kr,1.280000E+07_kr,1.300000E+07_kr,&
     1.320000E+07_kr,1.340000E+07_kr,1.360000E+07_kr,1.380000E+07_kr,&
     1.400000E+07_kr,1.420000E+07_kr,1.440000E+07_kr,1.460000E+07_kr,&
     1.480000E+07_kr,1.500000E+07_kr,1.520000E+07_kr,1.540000E+07_kr,&
     1.560000E+07_kr,1.580000E+07_kr,1.600000E+07_kr,1.620000E+07_kr,&
     1.640000E+07_kr,1.660000E+07_kr,1.680000E+07_kr,1.700000E+07_kr,&
     1.720000E+07_kr,1.740000E+07_kr,1.760000E+07_kr,1.780000E+07_kr,&
     1.800000E+07_kr,1.820000E+07_kr,1.840000E+07_kr,1.860000E+07_kr,&
     1.880000E+07_kr,1.900000E+07_kr,1.920000E+07_kr,1.940000E+07_kr,&
     1.960000E+07_kr,1.980000E+07_kr,2.000000E+07_kr,2.100000E+07_kr,&
     2.200000E+07_kr,2.300000E+07_kr,2.400000E+07_kr,2.500000E+07_kr,&
     2.600000E+07_kr,2.700000E+07_kr,2.800000E+07_kr,2.900000E+07_kr,&
     3.000000E+07_kr/)
   real(kr),dimension(1026),parameter::eg30=(/&
     1.000000e-05_kr,1.047129e-05_kr,1.096478e-05_kr,1.148154e-05_kr,&
     1.202264e-05_kr,1.258925e-05_kr,1.318257e-05_kr,1.380384e-05_kr,&
     1.445440e-05_kr,1.513561e-05_kr,1.584893e-05_kr,1.659587e-05_kr,&
     1.737801e-05_kr,1.819701e-05_kr,1.905461e-05_kr,1.995262e-05_kr,&
     2.089296e-05_kr,2.187762e-05_kr,2.290868e-05_kr,2.398833e-05_kr,&
     2.511886e-05_kr,2.630268e-05_kr,2.754229e-05_kr,2.884032e-05_kr,&
     3.019952e-05_kr,3.162278e-05_kr,3.311311e-05_kr,3.467369e-05_kr,&
     3.630781e-05_kr,3.801894e-05_kr,3.981072e-05_kr,4.168694e-05_kr,&
     4.365158e-05_kr,4.570882e-05_kr,4.786301e-05_kr,5.011872e-05_kr,&
     5.248075e-05_kr,5.495409e-05_kr,5.754399e-05_kr,6.025596e-05_kr,&
     6.309573e-05_kr,6.606934e-05_kr,6.918310e-05_kr,7.244360e-05_kr,&
     7.585776e-05_kr,7.943282e-05_kr,8.317638e-05_kr,8.709636e-05_kr,&
     9.120108e-05_kr,9.549926e-05_kr,1.000000e-04_kr,1.047129e-04_kr,&
     1.096478e-04_kr,1.148154e-04_kr,1.202264e-04_kr,1.258925e-04_kr,&
     1.318257e-04_kr,1.380384e-04_kr,1.445440e-04_kr,1.513561e-04_kr,&
     1.584893e-04_kr,1.659587e-04_kr,1.737801e-04_kr,1.819701e-04_kr,&
     1.905461e-04_kr,1.995262e-04_kr,2.089296e-04_kr,2.187762e-04_kr,&
     2.290868e-04_kr,2.398833e-04_kr,2.511886e-04_kr,2.630268e-04_kr,&
     2.754229e-04_kr,2.884032e-04_kr,3.019952e-04_kr,3.162278e-04_kr,&
     3.311311e-04_kr,3.467369e-04_kr,3.630781e-04_kr,3.801894e-04_kr,&
     3.981072e-04_kr,4.168694e-04_kr,4.365158e-04_kr,4.570882e-04_kr,&
     4.786301e-04_kr,5.011872e-04_kr,5.248075e-04_kr,5.495409e-04_kr,&
     5.754399e-04_kr,6.025596e-04_kr,6.309573e-04_kr,6.606934e-04_kr,&
     6.918310e-04_kr,7.244360e-04_kr,7.585776e-04_kr,7.943282e-04_kr,&
     8.317638e-04_kr,8.709636e-04_kr,9.120108e-04_kr,9.549926e-04_kr,&
     1.000000e-03_kr,1.047129e-03_kr,1.096478e-03_kr,1.148154e-03_kr,&
     1.202264e-03_kr,1.258925e-03_kr,1.318257e-03_kr,1.380384e-03_kr,&
     1.445440e-03_kr,1.513561e-03_kr,1.584893e-03_kr,1.659587e-03_kr,&
     1.737801e-03_kr,1.819701e-03_kr,1.905461e-03_kr,1.995262e-03_kr,&
     2.089296e-03_kr,2.187762e-03_kr,2.290868e-03_kr,2.398833e-03_kr,&
     2.511886e-03_kr,2.630268e-03_kr,2.754229e-03_kr,2.884032e-03_kr,&
     3.019952e-03_kr,3.162278e-03_kr,3.311311e-03_kr,3.467369e-03_kr,&
     3.630781e-03_kr,3.801894e-03_kr,3.981072e-03_kr,4.168694e-03_kr,&
     4.365158e-03_kr,4.570882e-03_kr,4.786301e-03_kr,5.011872e-03_kr,&
     5.248075e-03_kr,5.495409e-03_kr,5.754399e-03_kr,6.025596e-03_kr,&
     6.309573e-03_kr,6.606934e-03_kr,6.918310e-03_kr,7.244360e-03_kr,&
     7.585776e-03_kr,7.943282e-03_kr,8.317638e-03_kr,8.709636e-03_kr,&
     9.120108e-03_kr,9.549926e-03_kr,1.000000e-02_kr,1.047129e-02_kr,&
     1.096478e-02_kr,1.148154e-02_kr,1.202264e-02_kr,1.258925e-02_kr,&
     1.318257e-02_kr,1.380384e-02_kr,1.445440e-02_kr,1.513561e-02_kr,&
     1.584893e-02_kr,1.659587e-02_kr,1.737801e-02_kr,1.819701e-02_kr,&
     1.905461e-02_kr,1.995262e-02_kr,2.089296e-02_kr,2.187762e-02_kr,&
     2.290868e-02_kr,2.398833e-02_kr,2.511886e-02_kr,2.630268e-02_kr,&
     2.754229e-02_kr,2.884032e-02_kr,3.019952e-02_kr,3.162278e-02_kr,&
     3.311311e-02_kr,3.467369e-02_kr,3.630781e-02_kr,3.801894e-02_kr,&
     3.981072e-02_kr,4.168694e-02_kr,4.365158e-02_kr,4.570882e-02_kr,&
     4.786301e-02_kr,5.011872e-02_kr,5.248075e-02_kr,5.495409e-02_kr,&
     5.754399e-02_kr,6.025596e-02_kr,6.309573e-02_kr,6.606934e-02_kr,&
     6.918310e-02_kr,7.244360e-02_kr,7.585776e-02_kr,7.943282e-02_kr,&
     8.317638e-02_kr,8.709636e-02_kr,9.120108e-02_kr,9.549926e-02_kr,&
     1.000000e-01_kr,1.047129e-01_kr,1.096478e-01_kr,1.148154e-01_kr,&
     1.202264e-01_kr,1.258925e-01_kr,1.318257e-01_kr,1.380384e-01_kr,&
     1.445440e-01_kr,1.513561e-01_kr,1.584893e-01_kr,1.659587e-01_kr,&
     1.737801e-01_kr,1.819701e-01_kr,1.905461e-01_kr,1.995262e-01_kr,&
     2.089296e-01_kr,2.187762e-01_kr,2.290868e-01_kr,2.398833e-01_kr,&
     2.511886e-01_kr,2.630268e-01_kr,2.754229e-01_kr,2.884032e-01_kr,&
     3.019952e-01_kr,3.162278e-01_kr,3.311311e-01_kr,3.467369e-01_kr,&
     3.630781e-01_kr,3.801894e-01_kr,3.981072e-01_kr,4.168694e-01_kr,&
     4.365158e-01_kr,4.570882e-01_kr,4.786301e-01_kr,5.011872e-01_kr,&
     5.248075e-01_kr,5.500000e-01_kr,5.750000e-01_kr,6.000000e-01_kr,&
     6.250000e-01_kr,6.500000e-01_kr,6.750000e-01_kr,7.000000e-01_kr,&
     7.250000e-01_kr,7.500000e-01_kr,7.750000e-01_kr,8.000000e-01_kr,&
     8.250000e-01_kr,8.500000e-01_kr,8.750000e-01_kr,9.000000e-01_kr,&
     9.250000e-01_kr,9.500000e-01_kr,9.750000e-01_kr,1.000000e+00_kr,&
     1.025000e+00_kr,1.050000e+00_kr,1.075000e+00_kr,1.100000e+00_kr,&
     1.125000e+00_kr,1.150000e+00_kr,1.175000e+00_kr,1.200000e+00_kr,&
     1.225000e+00_kr,1.250000e+00_kr,1.275000e+00_kr,1.300000e+00_kr,&
     1.325000e+00_kr,1.350000e+00_kr,1.375000e+00_kr,1.400000e+00_kr,&
     1.425000e+00_kr,1.450000e+00_kr,1.475000e+00_kr,1.500000e+00_kr,&
     1.525000e+00_kr,1.550000e+00_kr,1.575000e+00_kr,1.600000e+00_kr,&
     1.625000e+00_kr,1.650000e+00_kr,1.675000e+00_kr,1.700000e+00_kr,&
     1.725000e+00_kr,1.750000e+00_kr,1.775000e+00_kr,1.800000e+00_kr,&
     1.825000e+00_kr,1.850000e+00_kr,1.875000e+00_kr,1.900000e+00_kr,&
     1.925000e+00_kr,1.950000e+00_kr,1.975000e+00_kr,2.000000e+00_kr,&
     2.025000e+00_kr,2.050000e+00_kr,2.075000e+00_kr,2.100000e+00_kr,&
     2.125000e+00_kr,2.150000e+00_kr,2.175000e+00_kr,2.200000e+00_kr,&
     2.225000e+00_kr,2.250000e+00_kr,2.275000e+00_kr,2.300000e+00_kr,&
     2.325000e+00_kr,2.350000e+00_kr,2.375000e+00_kr,2.400000e+00_kr,&
     2.425000e+00_kr,2.450000e+00_kr,2.475000e+00_kr,2.500000e+00_kr,&
     2.525000e+00_kr,2.550000e+00_kr,2.575000e+00_kr,2.600000e+00_kr,&
     2.625000e+00_kr,2.650000e+00_kr,2.675000e+00_kr,2.700000e+00_kr,&
     2.725000e+00_kr,2.750000e+00_kr,2.775000e+00_kr,2.800000e+00_kr,&
     2.825000e+00_kr,2.850000e+00_kr,2.875000e+00_kr,2.900000e+00_kr,&
     2.925000e+00_kr,2.950000e+00_kr,2.975000e+00_kr,3.000000e+00_kr,&
     3.025000e+00_kr,3.050000e+00_kr,3.075000e+00_kr,3.100000e+00_kr,&
     3.125000e+00_kr,3.150000e+00_kr,3.175000e+00_kr,3.200000e+00_kr,&
     3.225000e+00_kr,3.250000e+00_kr,3.275000e+00_kr,3.300000e+00_kr,&
     3.325000e+00_kr,3.350000e+00_kr,3.375000e+00_kr,3.400000e+00_kr,&
     3.425000e+00_kr,3.450000e+00_kr,3.475000e+00_kr,3.500000e+00_kr,&
     3.525000e+00_kr,3.550000e+00_kr,3.575000e+00_kr,3.600000e+00_kr,&
     3.625000e+00_kr,3.650000e+00_kr,3.675000e+00_kr,3.700000e+00_kr,&
     3.725000e+00_kr,3.750000e+00_kr,3.775000e+00_kr,3.800000e+00_kr,&
     3.825000e+00_kr,3.850000e+00_kr,3.875000e+00_kr,3.900000e+00_kr,&
     3.925000e+00_kr,3.950000e+00_kr,3.975000e+00_kr,4.000000e+00_kr,&
     4.025000e+00_kr,4.050000e+00_kr,4.075000e+00_kr,4.100000e+00_kr,&
     4.125000e+00_kr,4.150000e+00_kr,4.175000e+00_kr,4.200000e+00_kr,&
     4.225000e+00_kr,4.250000e+00_kr,4.275000e+00_kr,4.300000e+00_kr,&
     4.325000e+00_kr,4.350000e+00_kr,4.375000e+00_kr,4.400000e+00_kr,&
     4.425000e+00_kr,4.450000e+00_kr,4.475000e+00_kr,4.500000e+00_kr,&
     4.525000e+00_kr,4.550000e+00_kr,4.575000e+00_kr,4.600000e+00_kr,&
     4.625000e+00_kr,4.650000e+00_kr,4.675000e+00_kr,4.700000e+00_kr,&
     4.725000e+00_kr,4.750000e+00_kr,4.775000e+00_kr,4.800000e+00_kr,&
     4.825000e+00_kr,4.850000e+00_kr,4.875000e+00_kr,4.900000e+00_kr,&
     4.925000e+00_kr,4.950000e+00_kr,4.975000e+00_kr,5.000000e+00_kr,&
     5.025000e+00_kr,5.050000e+00_kr,5.075000e+00_kr,5.100000e+00_kr,&
     5.125000e+00_kr,5.150000e+00_kr,5.175000e+00_kr,5.200000e+00_kr,&
     5.225000e+00_kr,5.250000e+00_kr,5.275000e+00_kr,5.300000e+00_kr,&
     5.325000e+00_kr,5.350000e+00_kr,5.375000e+00_kr,5.400000e+00_kr,&
     5.425000e+00_kr,5.450000e+00_kr,5.475000e+00_kr,5.500000e+00_kr,&
     5.525000e+00_kr,5.550000e+00_kr,5.575000e+00_kr,5.600000e+00_kr,&
     5.625000e+00_kr,5.650000e+00_kr,5.675000e+00_kr,5.700000e+00_kr,&
     5.725000e+00_kr,5.750000e+00_kr,5.775000e+00_kr,5.800000e+00_kr,&
     5.825000e+00_kr,5.850000e+00_kr,5.875000e+00_kr,5.900000e+00_kr,&
     5.925000e+00_kr,5.950000e+00_kr,5.975000e+00_kr,6.000000e+00_kr,&
     6.025000e+00_kr,6.050000e+00_kr,6.075000e+00_kr,6.100000e+00_kr,&
     6.125000e+00_kr,6.150000e+00_kr,6.175000e+00_kr,6.200000e+00_kr,&
     6.225000e+00_kr,6.250000e+00_kr,6.275000e+00_kr,6.300000e+00_kr,&
     6.325000e+00_kr,6.350000e+00_kr,6.375000e+00_kr,6.400000e+00_kr,&
     6.425000e+00_kr,6.450000e+00_kr,6.475000e+00_kr,6.500000e+00_kr,&
     6.525000e+00_kr,6.550000e+00_kr,6.575000e+00_kr,6.600000e+00_kr,&
     6.625000e+00_kr,6.650000e+00_kr,6.675000e+00_kr,6.700000e+00_kr,&
     6.725000e+00_kr,6.750000e+00_kr,6.775000e+00_kr,6.800000e+00_kr,&
     6.825000e+00_kr,6.850000e+00_kr,6.875000e+00_kr,6.900000e+00_kr,&
     6.925000e+00_kr,6.950000e+00_kr,6.975000e+00_kr,7.000000e+00_kr,&
     7.025000e+00_kr,7.050000e+00_kr,7.075000e+00_kr,7.100000e+00_kr,&
     7.125000e+00_kr,7.150000e+00_kr,7.175000e+00_kr,7.200000e+00_kr,&
     7.225000e+00_kr,7.250000e+00_kr,7.275000e+00_kr,7.300000e+00_kr,&
     7.325000e+00_kr,7.350000e+00_kr,7.375000e+00_kr,7.400000e+00_kr,&
     7.425000e+00_kr,7.450000e+00_kr,7.475000e+00_kr,7.500000e+00_kr,&
     7.525000e+00_kr,7.550000e+00_kr,7.575000e+00_kr,7.600000e+00_kr,&
     7.625000e+00_kr,7.650000e+00_kr,7.675000e+00_kr,7.700000e+00_kr,&
     7.725000e+00_kr,7.750000e+00_kr,7.775000e+00_kr,7.800000e+00_kr,&
     7.825000e+00_kr,7.850000e+00_kr,7.875000e+00_kr,7.900000e+00_kr,&
     7.925000e+00_kr,7.950000e+00_kr,7.975000e+00_kr,8.000000e+00_kr,&
     8.025000e+00_kr,8.050000e+00_kr,8.075000e+00_kr,8.100000e+00_kr,&
     8.125000e+00_kr,8.150000e+00_kr,8.175000e+00_kr,8.200000e+00_kr,&
     8.225000e+00_kr,8.250000e+00_kr,8.275000e+00_kr,8.300000e+00_kr,&
     8.325000e+00_kr,8.350000e+00_kr,8.375000e+00_kr,8.400000e+00_kr,&
     8.425000e+00_kr,8.450000e+00_kr,8.475000e+00_kr,8.500000e+00_kr,&
     8.525000e+00_kr,8.550000e+00_kr,8.575000e+00_kr,8.600000e+00_kr,&
     8.625000e+00_kr,8.650000e+00_kr,8.675000e+00_kr,8.700000e+00_kr,&
     8.725000e+00_kr,8.750000e+00_kr,8.775000e+00_kr,8.800000e+00_kr,&
     8.825000e+00_kr,8.850000e+00_kr,8.875000e+00_kr,8.900000e+00_kr,&
     8.925000e+00_kr,8.950000e+00_kr,8.975000e+00_kr,9.000000e+00_kr,&
     9.025000e+00_kr,9.050000e+00_kr,9.075000e+00_kr,9.100000e+00_kr,&
     9.125000e+00_kr,9.150000e+00_kr,9.175000e+00_kr,9.200000e+00_kr,&
     9.225000e+00_kr,9.250000e+00_kr,9.275000e+00_kr,9.300000e+00_kr,&
     9.325000e+00_kr,9.350000e+00_kr,9.375000e+00_kr,9.400000e+00_kr,&
     9.425000e+00_kr,9.450000e+00_kr,9.475000e+00_kr,9.500000e+00_kr,&
     9.525000e+00_kr,9.550000e+00_kr,9.575000e+00_kr,9.600000e+00_kr,&
     9.625000e+00_kr,9.650000e+00_kr,9.675000e+00_kr,9.700000e+00_kr,&
     9.725000e+00_kr,9.750000e+00_kr,9.775000e+00_kr,9.800000e+00_kr,&
     9.825000e+00_kr,9.850000e+00_kr,9.875000e+00_kr,9.900000e+00_kr,&
     9.925000e+00_kr,9.950000e+00_kr,9.975000e+00_kr,1.000000e+01_kr,&
     1.047129e+01_kr,1.096478e+01_kr,1.148154e+01_kr,1.202264e+01_kr,&
     1.258925e+01_kr,1.318257e+01_kr,1.380384e+01_kr,1.445440e+01_kr,&
     1.513561e+01_kr,1.584893e+01_kr,1.659587e+01_kr,1.737801e+01_kr,&
     1.819701e+01_kr,1.905461e+01_kr,1.995262e+01_kr,2.089296e+01_kr,&
     2.187762e+01_kr,2.290868e+01_kr,2.398833e+01_kr,2.511886e+01_kr,&
     2.630268e+01_kr,2.754229e+01_kr,2.884032e+01_kr,3.019952e+01_kr,&
     3.162278e+01_kr,3.311311e+01_kr,3.467369e+01_kr,3.630781e+01_kr,&
     3.801894e+01_kr,3.981072e+01_kr,4.168694e+01_kr,4.365158e+01_kr,&
     4.570882e+01_kr,4.786301e+01_kr,5.011872e+01_kr,5.248075e+01_kr,&
     5.495409e+01_kr,5.754399e+01_kr,6.025596e+01_kr,6.309573e+01_kr,&
     6.606934e+01_kr,6.918310e+01_kr,7.244360e+01_kr,7.585776e+01_kr,&
     7.943282e+01_kr,8.317638e+01_kr,8.709636e+01_kr,9.120108e+01_kr,&
     9.549926e+01_kr,1.000000e+02_kr,1.047129e+02_kr,1.096478e+02_kr,&
     1.148154e+02_kr,1.202264e+02_kr,1.258925e+02_kr,1.318257e+02_kr,&
     1.380384e+02_kr,1.445440e+02_kr,1.513561e+02_kr,1.584893e+02_kr,&
     1.659587e+02_kr,1.737801e+02_kr,1.819701e+02_kr,1.905461e+02_kr,&
     1.995262e+02_kr,2.089296e+02_kr,2.187762e+02_kr,2.290868e+02_kr,&
     2.398833e+02_kr,2.511886e+02_kr,2.630268e+02_kr,2.754229e+02_kr,&
     2.884032e+02_kr,3.019952e+02_kr,3.162278e+02_kr,3.311311e+02_kr,&
     3.467369e+02_kr,3.630781e+02_kr,3.801894e+02_kr,3.981072e+02_kr,&
     4.168694e+02_kr,4.365158e+02_kr,4.570882e+02_kr,4.786301e+02_kr,&
     5.011872e+02_kr,5.248075e+02_kr,5.495409e+02_kr,5.754399e+02_kr,&
     6.025596e+02_kr,6.309573e+02_kr,6.606934e+02_kr,6.918310e+02_kr,&
     7.244360e+02_kr,7.585776e+02_kr,7.943282e+02_kr,8.317638e+02_kr,&
     8.709636e+02_kr,9.120108e+02_kr,9.549926e+02_kr,1.000000e+03_kr,&
     1.047129e+03_kr,1.096478e+03_kr,1.148154e+03_kr,1.202264e+03_kr,&
     1.258925e+03_kr,1.318257e+03_kr,1.380384e+03_kr,1.445440e+03_kr,&
     1.513561e+03_kr,1.584893e+03_kr,1.659587e+03_kr,1.737801e+03_kr,&
     1.819701e+03_kr,1.905461e+03_kr,1.995262e+03_kr,2.089296e+03_kr,&
     2.187762e+03_kr,2.290868e+03_kr,2.398833e+03_kr,2.511886e+03_kr,&
     2.630268e+03_kr,2.754229e+03_kr,2.884032e+03_kr,3.019952e+03_kr,&
     3.162278e+03_kr,3.311311e+03_kr,3.467369e+03_kr,3.630781e+03_kr,&
     3.801894e+03_kr,3.981072e+03_kr,4.168694e+03_kr,4.365158e+03_kr,&
     4.570882e+03_kr,4.786301e+03_kr,5.011872e+03_kr,5.248075e+03_kr,&
     5.495409e+03_kr,5.754399e+03_kr,6.025596e+03_kr,6.309573e+03_kr,&
     6.606934e+03_kr,6.918310e+03_kr,7.244360e+03_kr,7.585776e+03_kr,&
     7.943282e+03_kr,8.317638e+03_kr,8.709636e+03_kr,9.120108e+03_kr,&
     9.549926e+03_kr,1.000000e+04_kr,1.047129e+04_kr,1.096478e+04_kr,&
     1.148154e+04_kr,1.202264e+04_kr,1.258925e+04_kr,1.318257e+04_kr,&
     1.380384e+04_kr,1.445440e+04_kr,1.513561e+04_kr,1.584893e+04_kr,&
     1.659587e+04_kr,1.737801e+04_kr,1.819701e+04_kr,1.905461e+04_kr,&
     1.995262e+04_kr,2.089296e+04_kr,2.187762e+04_kr,2.290868e+04_kr,&
     2.398833e+04_kr,2.511886e+04_kr,2.630268e+04_kr,2.754229e+04_kr,&
     2.884032e+04_kr,3.019952e+04_kr,3.162278e+04_kr,3.311311e+04_kr,&
     3.467369e+04_kr,3.630781e+04_kr,3.801894e+04_kr,3.981072e+04_kr,&
     4.168694e+04_kr,4.365158e+04_kr,4.570882e+04_kr,4.786301e+04_kr,&
     5.011872e+04_kr,5.248075e+04_kr,5.495409e+04_kr,5.754399e+04_kr,&
     6.025596e+04_kr,6.309573e+04_kr,6.606934e+04_kr,6.918310e+04_kr,&
     7.244360e+04_kr,7.585776e+04_kr,7.943282e+04_kr,8.317638e+04_kr,&
     8.709636e+04_kr,9.120108e+04_kr,9.549926e+04_kr,1.000000e+05_kr,&
     1.047129e+05_kr,1.096478e+05_kr,1.148154e+05_kr,1.202264e+05_kr,&
     1.258925e+05_kr,1.318257e+05_kr,1.380384e+05_kr,1.445440e+05_kr,&
     1.513561e+05_kr,1.584893e+05_kr,1.659587e+05_kr,1.737801e+05_kr,&
     1.819701e+05_kr,1.905461e+05_kr,1.995262e+05_kr,2.089296e+05_kr,&
     2.187762e+05_kr,2.290868e+05_kr,2.398833e+05_kr,2.511886e+05_kr,&
     2.630268e+05_kr,2.754229e+05_kr,2.884032e+05_kr,3.019952e+05_kr,&
     3.162278e+05_kr,3.311311e+05_kr,3.467369e+05_kr,3.630781e+05_kr,&
     3.801894e+05_kr,3.981072e+05_kr,4.168694e+05_kr,4.365158e+05_kr,&
     4.570882e+05_kr,4.786301e+05_kr,5.011872e+05_kr,5.248075e+05_kr,&
     5.495409e+05_kr,5.754399e+05_kr,6.025596e+05_kr,6.309573e+05_kr,&
     6.606934e+05_kr,6.918310e+05_kr,7.244360e+05_kr,7.585776e+05_kr,&
     7.943282e+05_kr,8.317638e+05_kr,8.709636e+05_kr,9.120108e+05_kr,&
     9.549926e+05_kr,1.000000e+06_kr,1.047129e+06_kr,1.096478e+06_kr,&
     1.148154e+06_kr,1.202264e+06_kr,1.258925e+06_kr,1.318257e+06_kr,&
     1.380384e+06_kr,1.445440e+06_kr,1.513561e+06_kr,1.584893e+06_kr,&
     1.659587e+06_kr,1.737801e+06_kr,1.819701e+06_kr,1.905461e+06_kr,&
     1.995262e+06_kr,2.089296e+06_kr,2.187762e+06_kr,2.290868e+06_kr,&
     2.398833e+06_kr,2.511886e+06_kr,2.630268e+06_kr,2.754229e+06_kr,&
     2.884032e+06_kr,3.019952e+06_kr,3.162278e+06_kr,3.311311e+06_kr,&
     3.467369e+06_kr,3.630781e+06_kr,3.801894e+06_kr,3.981072e+06_kr,&
     4.168694e+06_kr,4.365158e+06_kr,4.570882e+06_kr,4.786301e+06_kr,&
     5.000000e+06_kr,5.200000e+06_kr,5.400000e+06_kr,5.600000e+06_kr,&
     5.800000e+06_kr,6.000000e+06_kr,6.200000e+06_kr,6.400000e+06_kr,&
     6.600000e+06_kr,6.800000e+06_kr,7.000000e+06_kr,7.200000e+06_kr,&
     7.400000e+06_kr,7.600000e+06_kr,7.800000e+06_kr,8.000000e+06_kr,&
     8.200000e+06_kr,8.400000e+06_kr,8.600000e+06_kr,8.800000e+06_kr,&
     9.000000e+06_kr,9.200000e+06_kr,9.400000e+06_kr,9.600000e+06_kr,&
     9.800000e+06_kr,1.000000e+07_kr,1.020000e+07_kr,1.040000e+07_kr,&
     1.060000e+07_kr,1.080000e+07_kr,1.100000e+07_kr,1.120000e+07_kr,&
     1.140000e+07_kr,1.160000e+07_kr,1.180000e+07_kr,1.200000e+07_kr,&
     1.220000e+07_kr,1.240000e+07_kr,1.260000e+07_kr,1.280000e+07_kr,&
     1.300000e+07_kr,1.320000e+07_kr,1.340000e+07_kr,1.360000e+07_kr,&
     1.380000e+07_kr,1.400000e+07_kr,1.420000e+07_kr,1.440000e+07_kr,&
     1.460000e+07_kr,1.480000e+07_kr,1.500000e+07_kr,1.520000e+07_kr,&
     1.540000e+07_kr,1.560000e+07_kr,1.580000e+07_kr,1.600000e+07_kr,&
     1.620000e+07_kr,1.640000e+07_kr,1.660000e+07_kr,1.680000e+07_kr,&
     1.700000e+07_kr,1.720000e+07_kr,1.740000e+07_kr,1.760000e+07_kr,&
     1.780000e+07_kr,1.800000e+07_kr,1.820000e+07_kr,1.840000e+07_kr,&
     1.860000e+07_kr,1.880000e+07_kr,1.900000e+07_kr,1.920000e+07_kr,&
     1.940000e+07_kr,1.960000e+07_kr,1.980000e+07_kr,2.000000e+07_kr,&
     2.020000e+07_kr,2.040000e+07_kr,2.060000e+07_kr,2.080000e+07_kr,&
     2.100000e+07_kr,2.120000e+07_kr,2.140000e+07_kr,2.160000e+07_kr,&
     2.180000e+07_kr,2.200000e+07_kr,2.220000e+07_kr,2.240000e+07_kr,&
     2.260000e+07_kr,2.280000e+07_kr,2.300000e+07_kr,2.320000e+07_kr,&
     2.340000e+07_kr,2.360000e+07_kr,2.380000e+07_kr,2.400000e+07_kr,&
     2.420000e+07_kr,2.440000e+07_kr,2.460000e+07_kr,2.480000e+07_kr,&
     2.500000e+07_kr,2.520000e+07_kr,2.540000e+07_kr,2.560000e+07_kr,&
     2.580000e+07_kr,2.600000e+07_kr,2.620000e+07_kr,2.640000e+07_kr,&
     2.660000e+07_kr,2.680000e+07_kr,2.700000e+07_kr,2.720000e+07_kr,&
     2.740000e+07_kr,2.760000e+07_kr,2.780000e+07_kr,2.800000e+07_kr,&
     2.820000e+07_kr,2.840000e+07_kr,2.860000e+07_kr,2.880000e+07_kr,&
     2.900000e+07_kr,2.920000e+07_kr,2.940000e+07_kr,2.960000e+07_kr,&
     2.980000e+07_kr,3.000000e+07_kr/)
   real(kr),dimension(1068),parameter::eg31=(/&
     1.000000e-05_kr,1.047129e-05_kr,1.096478e-05_kr,1.148154e-05_kr,&
     1.202264e-05_kr,1.258925e-05_kr,1.318257e-05_kr,1.380384e-05_kr,&
     1.445440e-05_kr,1.513561e-05_kr,1.584893e-05_kr,1.659587e-05_kr,&
     1.737801e-05_kr,1.819701e-05_kr,1.905461e-05_kr,1.995262e-05_kr,&
     2.089296e-05_kr,2.187762e-05_kr,2.290868e-05_kr,2.398833e-05_kr,&
     2.511886e-05_kr,2.630268e-05_kr,2.754229e-05_kr,2.884032e-05_kr,&
     3.019952e-05_kr,3.162278e-05_kr,3.311311e-05_kr,3.467369e-05_kr,&
     3.630781e-05_kr,3.801894e-05_kr,3.981072e-05_kr,4.168694e-05_kr,&
     4.365158e-05_kr,4.570882e-05_kr,4.786301e-05_kr,5.011872e-05_kr,&
     5.248075e-05_kr,5.495409e-05_kr,5.754399e-05_kr,6.025596e-05_kr,&
     6.309573e-05_kr,6.606934e-05_kr,6.918310e-05_kr,7.244360e-05_kr,&
     7.585776e-05_kr,7.943282e-05_kr,8.317638e-05_kr,8.709636e-05_kr,&
     9.120108e-05_kr,9.549926e-05_kr,1.000000e-04_kr,1.047129e-04_kr,&
     1.096478e-04_kr,1.148154e-04_kr,1.202264e-04_kr,1.258925e-04_kr,&
     1.318257e-04_kr,1.380384e-04_kr,1.445440e-04_kr,1.513561e-04_kr,&
     1.584893e-04_kr,1.659587e-04_kr,1.737801e-04_kr,1.819701e-04_kr,&
     1.905461e-04_kr,1.995262e-04_kr,2.089296e-04_kr,2.187762e-04_kr,&
     2.290868e-04_kr,2.398833e-04_kr,2.511886e-04_kr,2.630268e-04_kr,&
     2.754229e-04_kr,2.884032e-04_kr,3.019952e-04_kr,3.162278e-04_kr,&
     3.311311e-04_kr,3.467369e-04_kr,3.630781e-04_kr,3.801894e-04_kr,&
     3.981072e-04_kr,4.168694e-04_kr,4.365158e-04_kr,4.570882e-04_kr,&
     4.786301e-04_kr,5.011872e-04_kr,5.248075e-04_kr,5.495409e-04_kr,&
     5.754399e-04_kr,6.025596e-04_kr,6.309573e-04_kr,6.606934e-04_kr,&
     6.918310e-04_kr,7.244360e-04_kr,7.585776e-04_kr,7.943282e-04_kr,&
     8.317638e-04_kr,8.709636e-04_kr,9.120108e-04_kr,9.549926e-04_kr,&
     1.000000e-03_kr,1.047129e-03_kr,1.096478e-03_kr,1.148154e-03_kr,&
     1.202264e-03_kr,1.258925e-03_kr,1.318257e-03_kr,1.380384e-03_kr,&
     1.445440e-03_kr,1.513561e-03_kr,1.584893e-03_kr,1.659587e-03_kr,&
     1.737801e-03_kr,1.819701e-03_kr,1.905461e-03_kr,1.995262e-03_kr,&
     2.089296e-03_kr,2.187762e-03_kr,2.290868e-03_kr,2.398833e-03_kr,&
     2.511886e-03_kr,2.630268e-03_kr,2.754229e-03_kr,2.884032e-03_kr,&
     3.019952e-03_kr,3.162278e-03_kr,3.311311e-03_kr,3.467369e-03_kr,&
     3.630781e-03_kr,3.801894e-03_kr,3.981072e-03_kr,4.168694e-03_kr,&
     4.365158e-03_kr,4.570882e-03_kr,4.786301e-03_kr,5.011872e-03_kr,&
     5.248075e-03_kr,5.495409e-03_kr,5.754399e-03_kr,6.025596e-03_kr,&
     6.309573e-03_kr,6.606934e-03_kr,6.918310e-03_kr,7.244360e-03_kr,&
     7.585776e-03_kr,7.943282e-03_kr,8.317638e-03_kr,8.709636e-03_kr,&
     9.120108e-03_kr,9.549926e-03_kr,1.000000e-02_kr,1.047129e-02_kr,&
     1.096478e-02_kr,1.148154e-02_kr,1.202264e-02_kr,1.258925e-02_kr,&
     1.318257e-02_kr,1.380384e-02_kr,1.445440e-02_kr,1.513561e-02_kr,&
     1.584893e-02_kr,1.659587e-02_kr,1.737801e-02_kr,1.819701e-02_kr,&
     1.905461e-02_kr,1.995262e-02_kr,2.089296e-02_kr,2.187762e-02_kr,&
     2.290868e-02_kr,2.398833e-02_kr,2.511886e-02_kr,2.630268e-02_kr,&
     2.754229e-02_kr,2.884032e-02_kr,3.019952e-02_kr,3.162278e-02_kr,&
     3.311311e-02_kr,3.467369e-02_kr,3.630781e-02_kr,3.801894e-02_kr,&
     3.981072e-02_kr,4.168694e-02_kr,4.365158e-02_kr,4.570882e-02_kr,&
     4.786301e-02_kr,5.011872e-02_kr,5.248075e-02_kr,5.495409e-02_kr,&
     5.754399e-02_kr,6.025596e-02_kr,6.309573e-02_kr,6.606934e-02_kr,&
     6.918310e-02_kr,7.244360e-02_kr,7.585776e-02_kr,7.943282e-02_kr,&
     8.317638e-02_kr,8.709636e-02_kr,9.120108e-02_kr,9.549926e-02_kr,&
     1.000000e-01_kr,1.047129e-01_kr,1.096478e-01_kr,1.148154e-01_kr,&
     1.202264e-01_kr,1.258925e-01_kr,1.318257e-01_kr,1.380384e-01_kr,&
     1.445440e-01_kr,1.513561e-01_kr,1.584893e-01_kr,1.659587e-01_kr,&
     1.737801e-01_kr,1.819701e-01_kr,1.905461e-01_kr,1.995262e-01_kr,&
     2.089296e-01_kr,2.187762e-01_kr,2.290868e-01_kr,2.398833e-01_kr,&
     2.511886e-01_kr,2.630268e-01_kr,2.754229e-01_kr,2.884032e-01_kr,&
     3.019952e-01_kr,3.162278e-01_kr,3.311311e-01_kr,3.467369e-01_kr,&
     3.630781e-01_kr,3.801894e-01_kr,3.981072e-01_kr,4.168694e-01_kr,&
     4.365158e-01_kr,4.570882e-01_kr,4.786301e-01_kr,5.011872e-01_kr,&
     5.248075e-01_kr,5.500000e-01_kr,5.750000e-01_kr,6.000000e-01_kr,&
     6.250000e-01_kr,6.500000e-01_kr,6.750000e-01_kr,7.000000e-01_kr,&
     7.250000e-01_kr,7.500000e-01_kr,7.750000e-01_kr,8.000000e-01_kr,&
     8.250000e-01_kr,8.500000e-01_kr,8.750000e-01_kr,9.000000e-01_kr,&
     9.250000e-01_kr,9.500000e-01_kr,9.750000e-01_kr,1.000000e+00_kr,&
     1.025000e+00_kr,1.050000e+00_kr,1.075000e+00_kr,1.100000e+00_kr,&
     1.125000e+00_kr,1.150000e+00_kr,1.175000e+00_kr,1.200000e+00_kr,&
     1.225000e+00_kr,1.250000e+00_kr,1.275000e+00_kr,1.300000e+00_kr,&
     1.325000e+00_kr,1.350000e+00_kr,1.375000e+00_kr,1.400000e+00_kr,&
     1.425000e+00_kr,1.450000e+00_kr,1.475000e+00_kr,1.500000e+00_kr,&
     1.525000e+00_kr,1.550000e+00_kr,1.575000e+00_kr,1.600000e+00_kr,&
     1.625000e+00_kr,1.650000e+00_kr,1.675000e+00_kr,1.700000e+00_kr,&
     1.725000e+00_kr,1.750000e+00_kr,1.775000e+00_kr,1.800000e+00_kr,&
     1.825000e+00_kr,1.850000e+00_kr,1.875000e+00_kr,1.900000e+00_kr,&
     1.925000e+00_kr,1.950000e+00_kr,1.975000e+00_kr,2.000000e+00_kr,&
     2.025000e+00_kr,2.050000e+00_kr,2.075000e+00_kr,2.100000e+00_kr,&
     2.125000e+00_kr,2.150000e+00_kr,2.175000e+00_kr,2.200000e+00_kr,&
     2.225000e+00_kr,2.250000e+00_kr,2.275000e+00_kr,2.300000e+00_kr,&
     2.325000e+00_kr,2.350000e+00_kr,2.375000e+00_kr,2.400000e+00_kr,&
     2.425000e+00_kr,2.450000e+00_kr,2.475000e+00_kr,2.500000e+00_kr,&
     2.525000e+00_kr,2.550000e+00_kr,2.575000e+00_kr,2.600000e+00_kr,&
     2.625000e+00_kr,2.650000e+00_kr,2.675000e+00_kr,2.700000e+00_kr,&
     2.725000e+00_kr,2.750000e+00_kr,2.775000e+00_kr,2.800000e+00_kr,&
     2.825000e+00_kr,2.850000e+00_kr,2.875000e+00_kr,2.900000e+00_kr,&
     2.925000e+00_kr,2.950000e+00_kr,2.975000e+00_kr,3.000000e+00_kr,&
     3.025000e+00_kr,3.050000e+00_kr,3.075000e+00_kr,3.100000e+00_kr,&
     3.125000e+00_kr,3.150000e+00_kr,3.175000e+00_kr,3.200000e+00_kr,&
     3.225000e+00_kr,3.250000e+00_kr,3.275000e+00_kr,3.300000e+00_kr,&
     3.325000e+00_kr,3.350000e+00_kr,3.375000e+00_kr,3.400000e+00_kr,&
     3.425000e+00_kr,3.450000e+00_kr,3.475000e+00_kr,3.500000e+00_kr,&
     3.525000e+00_kr,3.550000e+00_kr,3.575000e+00_kr,3.600000e+00_kr,&
     3.625000e+00_kr,3.650000e+00_kr,3.675000e+00_kr,3.700000e+00_kr,&
     3.725000e+00_kr,3.750000e+00_kr,3.775000e+00_kr,3.800000e+00_kr,&
     3.825000e+00_kr,3.850000e+00_kr,3.875000e+00_kr,3.900000e+00_kr,&
     3.925000e+00_kr,3.950000e+00_kr,3.975000e+00_kr,4.000000e+00_kr,&
     4.025000e+00_kr,4.050000e+00_kr,4.075000e+00_kr,4.100000e+00_kr,&
     4.125000e+00_kr,4.150000e+00_kr,4.175000e+00_kr,4.200000e+00_kr,&
     4.225000e+00_kr,4.250000e+00_kr,4.275000e+00_kr,4.300000e+00_kr,&
     4.325000e+00_kr,4.350000e+00_kr,4.375000e+00_kr,4.400000e+00_kr,&
     4.425000e+00_kr,4.450000e+00_kr,4.475000e+00_kr,4.500000e+00_kr,&
     4.525000e+00_kr,4.550000e+00_kr,4.575000e+00_kr,4.600000e+00_kr,&
     4.625000e+00_kr,4.650000e+00_kr,4.675000e+00_kr,4.700000e+00_kr,&
     4.725000e+00_kr,4.750000e+00_kr,4.775000e+00_kr,4.800000e+00_kr,&
     4.825000e+00_kr,4.850000e+00_kr,4.875000e+00_kr,4.900000e+00_kr,&
     4.925000e+00_kr,4.950000e+00_kr,4.975000e+00_kr,5.000000e+00_kr,&
     5.025000e+00_kr,5.050000e+00_kr,5.075000e+00_kr,5.100000e+00_kr,&
     5.125000e+00_kr,5.150000e+00_kr,5.175000e+00_kr,5.200000e+00_kr,&
     5.225000e+00_kr,5.250000e+00_kr,5.275000e+00_kr,5.300000e+00_kr,&
     5.325000e+00_kr,5.350000e+00_kr,5.375000e+00_kr,5.400000e+00_kr,&
     5.425000e+00_kr,5.450000e+00_kr,5.475000e+00_kr,5.500000e+00_kr,&
     5.525000e+00_kr,5.550000e+00_kr,5.575000e+00_kr,5.600000e+00_kr,&
     5.625000e+00_kr,5.650000e+00_kr,5.675000e+00_kr,5.700000e+00_kr,&
     5.725000e+00_kr,5.750000e+00_kr,5.775000e+00_kr,5.800000e+00_kr,&
     5.825000e+00_kr,5.850000e+00_kr,5.875000e+00_kr,5.900000e+00_kr,&
     5.925000e+00_kr,5.950000e+00_kr,5.975000e+00_kr,6.000000e+00_kr,&
     6.025000e+00_kr,6.050000e+00_kr,6.075000e+00_kr,6.100000e+00_kr,&
     6.125000e+00_kr,6.150000e+00_kr,6.175000e+00_kr,6.200000e+00_kr,&
     6.225000e+00_kr,6.250000e+00_kr,6.275000e+00_kr,6.300000e+00_kr,&
     6.325000e+00_kr,6.350000e+00_kr,6.375000e+00_kr,6.400000e+00_kr,&
     6.425000e+00_kr,6.450000e+00_kr,6.475000e+00_kr,6.500000e+00_kr,&
     6.525000e+00_kr,6.550000e+00_kr,6.575000e+00_kr,6.600000e+00_kr,&
     6.625000e+00_kr,6.650000e+00_kr,6.675000e+00_kr,6.700000e+00_kr,&
     6.725000e+00_kr,6.750000e+00_kr,6.775000e+00_kr,6.800000e+00_kr,&
     6.825000e+00_kr,6.850000e+00_kr,6.875000e+00_kr,6.900000e+00_kr,&
     6.925000e+00_kr,6.950000e+00_kr,6.975000e+00_kr,7.000000e+00_kr,&
     7.025000e+00_kr,7.050000e+00_kr,7.075000e+00_kr,7.100000e+00_kr,&
     7.125000e+00_kr,7.150000e+00_kr,7.175000e+00_kr,7.200000e+00_kr,&
     7.225000e+00_kr,7.250000e+00_kr,7.275000e+00_kr,7.300000e+00_kr,&
     7.325000e+00_kr,7.350000e+00_kr,7.375000e+00_kr,7.400000e+00_kr,&
     7.425000e+00_kr,7.450000e+00_kr,7.475000e+00_kr,7.500000e+00_kr,&
     7.525000e+00_kr,7.550000e+00_kr,7.575000e+00_kr,7.600000e+00_kr,&
     7.625000e+00_kr,7.650000e+00_kr,7.675000e+00_kr,7.700000e+00_kr,&
     7.725000e+00_kr,7.750000e+00_kr,7.775000e+00_kr,7.800000e+00_kr,&
     7.825000e+00_kr,7.850000e+00_kr,7.875000e+00_kr,7.900000e+00_kr,&
     7.925000e+00_kr,7.950000e+00_kr,7.975000e+00_kr,8.000000e+00_kr,&
     8.025000e+00_kr,8.050000e+00_kr,8.075000e+00_kr,8.100000e+00_kr,&
     8.125000e+00_kr,8.150000e+00_kr,8.175000e+00_kr,8.200000e+00_kr,&
     8.225000e+00_kr,8.250000e+00_kr,8.275000e+00_kr,8.300000e+00_kr,&
     8.325000e+00_kr,8.350000e+00_kr,8.375000e+00_kr,8.400000e+00_kr,&
     8.425000e+00_kr,8.450000e+00_kr,8.475000e+00_kr,8.500000e+00_kr,&
     8.525000e+00_kr,8.550000e+00_kr,8.575000e+00_kr,8.600000e+00_kr,&
     8.625000e+00_kr,8.650000e+00_kr,8.675000e+00_kr,8.700000e+00_kr,&
     8.725000e+00_kr,8.750000e+00_kr,8.775000e+00_kr,8.800000e+00_kr,&
     8.825000e+00_kr,8.850000e+00_kr,8.875000e+00_kr,8.900000e+00_kr,&
     8.925000e+00_kr,8.950000e+00_kr,8.975000e+00_kr,9.000000e+00_kr,&
     9.025000e+00_kr,9.050000e+00_kr,9.075000e+00_kr,9.100000e+00_kr,&
     9.125000e+00_kr,9.150000e+00_kr,9.175000e+00_kr,9.200000e+00_kr,&
     9.225000e+00_kr,9.250000e+00_kr,9.275000e+00_kr,9.300000e+00_kr,&
     9.325000e+00_kr,9.350000e+00_kr,9.375000e+00_kr,9.400000e+00_kr,&
     9.425000e+00_kr,9.450000e+00_kr,9.475000e+00_kr,9.500000e+00_kr,&
     9.525000e+00_kr,9.550000e+00_kr,9.575000e+00_kr,9.600000e+00_kr,&
     9.625000e+00_kr,9.650000e+00_kr,9.675000e+00_kr,9.700000e+00_kr,&
     9.725000e+00_kr,9.750000e+00_kr,9.775000e+00_kr,9.800000e+00_kr,&
     9.825000e+00_kr,9.850000e+00_kr,9.875000e+00_kr,9.900000e+00_kr,&
     9.925000e+00_kr,9.950000e+00_kr,9.975000e+00_kr,1.000000e+01_kr,&
     1.047129e+01_kr,1.096478e+01_kr,1.148154e+01_kr,1.202264e+01_kr,&
     1.258925e+01_kr,1.318257e+01_kr,1.380384e+01_kr,1.445440e+01_kr,&
     1.513561e+01_kr,1.584893e+01_kr,1.659587e+01_kr,1.737801e+01_kr,&
     1.819701e+01_kr,1.905461e+01_kr,1.995262e+01_kr,2.089296e+01_kr,&
     2.187762e+01_kr,2.290868e+01_kr,2.398833e+01_kr,2.511886e+01_kr,&
     2.630268e+01_kr,2.754229e+01_kr,2.884032e+01_kr,3.019952e+01_kr,&
     3.162278e+01_kr,3.311311e+01_kr,3.467369e+01_kr,3.630781e+01_kr,&
     3.801894e+01_kr,3.981072e+01_kr,4.168694e+01_kr,4.365158e+01_kr,&
     4.570882e+01_kr,4.786301e+01_kr,5.011872e+01_kr,5.248075e+01_kr,&
     5.495409e+01_kr,5.754399e+01_kr,6.025596e+01_kr,6.309573e+01_kr,&
     6.606934e+01_kr,6.918310e+01_kr,7.244360e+01_kr,7.585776e+01_kr,&
     7.943282e+01_kr,8.317638e+01_kr,8.709636e+01_kr,9.120108e+01_kr,&
     9.549926e+01_kr,1.000000e+02_kr,1.047129e+02_kr,1.096478e+02_kr,&
     1.148154e+02_kr,1.202264e+02_kr,1.258925e+02_kr,1.318257e+02_kr,&
     1.380384e+02_kr,1.445440e+02_kr,1.513561e+02_kr,1.584893e+02_kr,&
     1.659587e+02_kr,1.737801e+02_kr,1.819701e+02_kr,1.905461e+02_kr,&
     1.995262e+02_kr,2.089296e+02_kr,2.187762e+02_kr,2.290868e+02_kr,&
     2.398833e+02_kr,2.511886e+02_kr,2.630268e+02_kr,2.754229e+02_kr,&
     2.884032e+02_kr,3.019952e+02_kr,3.162278e+02_kr,3.311311e+02_kr,&
     3.467369e+02_kr,3.630781e+02_kr,3.801894e+02_kr,3.981072e+02_kr,&
     4.168694e+02_kr,4.365158e+02_kr,4.570882e+02_kr,4.786301e+02_kr,&
     5.011872e+02_kr,5.248075e+02_kr,5.495409e+02_kr,5.754399e+02_kr,&
     6.025596e+02_kr,6.309573e+02_kr,6.606934e+02_kr,6.918310e+02_kr,&
     7.244360e+02_kr,7.585776e+02_kr,7.943282e+02_kr,8.317638e+02_kr,&
     8.709636e+02_kr,9.120108e+02_kr,9.549926e+02_kr,1.000000e+03_kr,&
     1.047129e+03_kr,1.096478e+03_kr,1.148154e+03_kr,1.202264e+03_kr,&
     1.258925e+03_kr,1.318257e+03_kr,1.380384e+03_kr,1.445440e+03_kr,&
     1.513561e+03_kr,1.584893e+03_kr,1.659587e+03_kr,1.737801e+03_kr,&
     1.819701e+03_kr,1.905461e+03_kr,1.995262e+03_kr,2.089296e+03_kr,&
     2.187762e+03_kr,2.290868e+03_kr,2.398833e+03_kr,2.511886e+03_kr,&
     2.630268e+03_kr,2.754229e+03_kr,2.884032e+03_kr,3.019952e+03_kr,&
     3.162278e+03_kr,3.311311e+03_kr,3.467369e+03_kr,3.630781e+03_kr,&
     3.801894e+03_kr,3.981072e+03_kr,4.168694e+03_kr,4.365158e+03_kr,&
     4.570882e+03_kr,4.786301e+03_kr,5.011872e+03_kr,5.248075e+03_kr,&
     5.495409e+03_kr,5.754399e+03_kr,6.025596e+03_kr,6.309573e+03_kr,&
     6.606934e+03_kr,6.918310e+03_kr,7.244360e+03_kr,7.585776e+03_kr,&
     7.943282e+03_kr,8.317638e+03_kr,8.709636e+03_kr,9.120108e+03_kr,&
     9.549926e+03_kr,1.000000e+04_kr,1.047129e+04_kr,1.096478e+04_kr,&
     1.148154e+04_kr,1.202264e+04_kr,1.258925e+04_kr,1.318257e+04_kr,&
     1.380384e+04_kr,1.445440e+04_kr,1.513561e+04_kr,1.584893e+04_kr,&
     1.659587e+04_kr,1.737801e+04_kr,1.819701e+04_kr,1.905461e+04_kr,&
     1.995262e+04_kr,2.089296e+04_kr,2.187762e+04_kr,2.290868e+04_kr,&
     2.398833e+04_kr,2.511886e+04_kr,2.630268e+04_kr,2.754229e+04_kr,&
     2.884032e+04_kr,3.019952e+04_kr,3.162278e+04_kr,3.311311e+04_kr,&
     3.467369e+04_kr,3.630781e+04_kr,3.801894e+04_kr,3.981072e+04_kr,&
     4.168694e+04_kr,4.365158e+04_kr,4.570882e+04_kr,4.786301e+04_kr,&
     5.011872e+04_kr,5.248075e+04_kr,5.495409e+04_kr,5.754399e+04_kr,&
     6.025596e+04_kr,6.309573e+04_kr,6.606934e+04_kr,6.918310e+04_kr,&
     7.244360e+04_kr,7.585776e+04_kr,7.943282e+04_kr,8.317638e+04_kr,&
     8.709636e+04_kr,9.120108e+04_kr,9.549926e+04_kr,1.000000e+05_kr,&
     1.047129e+05_kr,1.096478e+05_kr,1.148154e+05_kr,1.202264e+05_kr,&
     1.258925e+05_kr,1.318257e+05_kr,1.380384e+05_kr,1.445440e+05_kr,&
     1.513561e+05_kr,1.584893e+05_kr,1.659587e+05_kr,1.737801e+05_kr,&
     1.819701e+05_kr,1.905461e+05_kr,1.995262e+05_kr,2.089296e+05_kr,&
     2.187762e+05_kr,2.290868e+05_kr,2.398833e+05_kr,2.511886e+05_kr,&
     2.630268e+05_kr,2.754229e+05_kr,2.884032e+05_kr,3.019952e+05_kr,&
     3.162278e+05_kr,3.311311e+05_kr,3.467369e+05_kr,3.630781e+05_kr,&
     3.801894e+05_kr,3.981072e+05_kr,4.168694e+05_kr,4.365158e+05_kr,&
     4.570882e+05_kr,4.786301e+05_kr,5.011872e+05_kr,5.248075e+05_kr,&
     5.495409e+05_kr,5.754399e+05_kr,6.025596e+05_kr,6.309573e+05_kr,&
     6.606934e+05_kr,6.918310e+05_kr,7.244360e+05_kr,7.585776e+05_kr,&
     7.943282e+05_kr,8.317638e+05_kr,8.709636e+05_kr,9.120108e+05_kr,&
     9.549926e+05_kr,1.000000e+06_kr,1.047129e+06_kr,1.096478e+06_kr,&
     1.148154e+06_kr,1.202264e+06_kr,1.258925e+06_kr,1.318257e+06_kr,&
     1.380384e+06_kr,1.445440e+06_kr,1.513561e+06_kr,1.584893e+06_kr,&
     1.659587e+06_kr,1.737801e+06_kr,1.819701e+06_kr,1.905461e+06_kr,&
     1.995262e+06_kr,2.089296e+06_kr,2.187762e+06_kr,2.290868e+06_kr,&
     2.398833e+06_kr,2.511886e+06_kr,2.630268e+06_kr,2.754229e+06_kr,&
     2.884032e+06_kr,3.019952e+06_kr,3.162278e+06_kr,3.311311e+06_kr,&
     3.467369e+06_kr,3.630781e+06_kr,3.801894e+06_kr,3.981072e+06_kr,&
     4.168694e+06_kr,4.365158e+06_kr,4.570882e+06_kr,4.786301e+06_kr,&
     5.000000e+06_kr,5.200000e+06_kr,5.400000e+06_kr,5.600000e+06_kr,&
     5.800000e+06_kr,6.000000e+06_kr,6.200000e+06_kr,6.400000e+06_kr,&
     6.600000e+06_kr,6.800000e+06_kr,7.000000e+06_kr,7.200000e+06_kr,&
     7.400000e+06_kr,7.600000e+06_kr,7.800000e+06_kr,8.000000e+06_kr,&
     8.200000e+06_kr,8.400000e+06_kr,8.600000e+06_kr,8.800000e+06_kr,&
     9.000000e+06_kr,9.200000e+06_kr,9.400000e+06_kr,9.600000e+06_kr,&
     9.800000e+06_kr,1.000000e+07_kr,1.020000e+07_kr,1.040000e+07_kr,&
     1.060000e+07_kr,1.080000e+07_kr,1.100000e+07_kr,1.120000e+07_kr,&
     1.140000e+07_kr,1.160000e+07_kr,1.180000e+07_kr,1.200000e+07_kr,&
     1.220000e+07_kr,1.240000e+07_kr,1.260000e+07_kr,1.280000e+07_kr,&
     1.300000e+07_kr,1.320000e+07_kr,1.340000e+07_kr,1.360000e+07_kr,&
     1.380000e+07_kr,1.400000e+07_kr,1.420000e+07_kr,1.440000e+07_kr,&
     1.460000e+07_kr,1.480000e+07_kr,1.500000e+07_kr,1.520000e+07_kr,&
     1.540000e+07_kr,1.560000e+07_kr,1.580000e+07_kr,1.600000e+07_kr,&
     1.620000e+07_kr,1.640000e+07_kr,1.660000e+07_kr,1.680000e+07_kr,&
     1.700000e+07_kr,1.720000e+07_kr,1.740000e+07_kr,1.760000e+07_kr,&
     1.780000e+07_kr,1.800000e+07_kr,1.820000e+07_kr,1.840000e+07_kr,&
     1.860000e+07_kr,1.880000e+07_kr,1.900000e+07_kr,1.920000e+07_kr,&
     1.940000e+07_kr,1.960000e+07_kr,1.980000e+07_kr,2.000000e+07_kr,&
     2.020000e+07_kr,2.040000e+07_kr,2.060000e+07_kr,2.080000e+07_kr,&
     2.100000e+07_kr,2.120000e+07_kr,2.140000e+07_kr,2.160000e+07_kr,&
     2.180000e+07_kr,2.200000e+07_kr,2.220000e+07_kr,2.240000e+07_kr,&
     2.260000e+07_kr,2.280000e+07_kr,2.300000e+07_kr,2.320000e+07_kr,&
     2.340000e+07_kr,2.360000e+07_kr,2.380000e+07_kr,2.400000e+07_kr,&
     2.420000e+07_kr,2.440000e+07_kr,2.460000e+07_kr,2.480000e+07_kr,&
     2.500000e+07_kr,2.520000e+07_kr,2.540000e+07_kr,2.560000e+07_kr,&
     2.580000e+07_kr,2.600000e+07_kr,2.620000e+07_kr,2.640000e+07_kr,&
     2.660000e+07_kr,2.680000e+07_kr,2.700000e+07_kr,2.720000e+07_kr,&
     2.740000e+07_kr,2.760000e+07_kr,2.780000e+07_kr,2.800000e+07_kr,&
     2.820000e+07_kr,2.840000e+07_kr,2.860000e+07_kr,2.880000e+07_kr,&
     2.900000e+07_kr,2.920000e+07_kr,2.940000e+07_kr,2.960000e+07_kr,&
     2.980000e+07_kr,3.000000e+07_kr,3.019952e+07_kr,3.162278e+07_kr,&
     3.311311e+07_kr,3.467369e+07_kr,3.630781e+07_kr,3.801894e+07_kr,&
     3.981072e+07_kr,4.168694e+07_kr,4.365158e+07_kr,4.570882e+07_kr,&
     4.786301e+07_kr,5.011872e+07_kr,5.248075e+07_kr,5.495409e+07_kr,&
     5.754399e+07_kr,6.025596e+07_kr,6.309573e+07_kr,6.606934e+07_kr,&
     6.918310e+07_kr,7.244360e+07_kr,7.585776e+07_kr,7.943282e+07_kr,&
     8.317638e+07_kr,8.709636e+07_kr,9.120108e+07_kr,9.549926e+07_kr,&
     1.000000e+08_kr,1.047129e+08_kr,1.096478e+08_kr,1.148154e+08_kr,&
     1.202264e+08_kr,1.258925e+08_kr,1.318257e+08_kr,1.380384e+08_kr,&
     1.445440e+08_kr,1.513561e+08_kr,1.584893e+08_kr,1.659587e+08_kr,&
     1.737801e+08_kr,1.819701e+08_kr,1.905461e+08_kr,1.995262e+08_kr/)
   real(kr),dimension(1103),parameter::eg32=(/&
     1.000000e-05_kr,1.047129e-05_kr,1.096478e-05_kr,1.148154e-05_kr,&
     1.202264e-05_kr,1.258925e-05_kr,1.318257e-05_kr,1.380384e-05_kr,&
     1.445440e-05_kr,1.513561e-05_kr,1.584893e-05_kr,1.659587e-05_kr,&
     1.737801e-05_kr,1.819701e-05_kr,1.905461e-05_kr,1.995262e-05_kr,&
     2.089296e-05_kr,2.187762e-05_kr,2.290868e-05_kr,2.398833e-05_kr,&
     2.511886e-05_kr,2.630268e-05_kr,2.754229e-05_kr,2.884032e-05_kr,&
     3.019952e-05_kr,3.162278e-05_kr,3.311311e-05_kr,3.467369e-05_kr,&
     3.630781e-05_kr,3.801894e-05_kr,3.981072e-05_kr,4.168694e-05_kr,&
     4.365158e-05_kr,4.570882e-05_kr,4.786301e-05_kr,5.011872e-05_kr,&
     5.248075e-05_kr,5.495409e-05_kr,5.754399e-05_kr,6.025596e-05_kr,&
     6.309573e-05_kr,6.606934e-05_kr,6.918310e-05_kr,7.244360e-05_kr,&
     7.585776e-05_kr,7.943282e-05_kr,8.317638e-05_kr,8.709636e-05_kr,&
     9.120108e-05_kr,9.549926e-05_kr,1.000000e-04_kr,1.047129e-04_kr,&
     1.096478e-04_kr,1.148154e-04_kr,1.202264e-04_kr,1.258925e-04_kr,&
     1.318257e-04_kr,1.380384e-04_kr,1.445440e-04_kr,1.513561e-04_kr,&
     1.584893e-04_kr,1.659587e-04_kr,1.737801e-04_kr,1.819701e-04_kr,&
     1.905461e-04_kr,1.995262e-04_kr,2.089296e-04_kr,2.187762e-04_kr,&
     2.290868e-04_kr,2.398833e-04_kr,2.511886e-04_kr,2.630268e-04_kr,&
     2.754229e-04_kr,2.884032e-04_kr,3.019952e-04_kr,3.162278e-04_kr,&
     3.311311e-04_kr,3.467369e-04_kr,3.630781e-04_kr,3.801894e-04_kr,&
     3.981072e-04_kr,4.168694e-04_kr,4.365158e-04_kr,4.570882e-04_kr,&
     4.786301e-04_kr,5.011872e-04_kr,5.248075e-04_kr,5.495409e-04_kr,&
     5.754399e-04_kr,6.025596e-04_kr,6.309573e-04_kr,6.606934e-04_kr,&
     6.918310e-04_kr,7.244360e-04_kr,7.585776e-04_kr,7.943282e-04_kr,&
     8.317638e-04_kr,8.709636e-04_kr,9.120108e-04_kr,9.549926e-04_kr,&
     1.000000e-03_kr,1.047129e-03_kr,1.096478e-03_kr,1.148154e-03_kr,&
     1.202264e-03_kr,1.258925e-03_kr,1.318257e-03_kr,1.380384e-03_kr,&
     1.445440e-03_kr,1.513561e-03_kr,1.584893e-03_kr,1.659587e-03_kr,&
     1.737801e-03_kr,1.819701e-03_kr,1.905461e-03_kr,1.995262e-03_kr,&
     2.089296e-03_kr,2.187762e-03_kr,2.290868e-03_kr,2.398833e-03_kr,&
     2.511886e-03_kr,2.630268e-03_kr,2.754229e-03_kr,2.884032e-03_kr,&
     3.019952e-03_kr,3.162278e-03_kr,3.311311e-03_kr,3.467369e-03_kr,&
     3.630781e-03_kr,3.801894e-03_kr,3.981072e-03_kr,4.168694e-03_kr,&
     4.365158e-03_kr,4.570882e-03_kr,4.786301e-03_kr,5.011872e-03_kr,&
     5.248075e-03_kr,5.495409e-03_kr,5.754399e-03_kr,6.025596e-03_kr,&
     6.309573e-03_kr,6.606934e-03_kr,6.918310e-03_kr,7.244360e-03_kr,&
     7.585776e-03_kr,7.943282e-03_kr,8.317638e-03_kr,8.709636e-03_kr,&
     9.120108e-03_kr,9.549926e-03_kr,1.000000e-02_kr,1.047129e-02_kr,&
     1.096478e-02_kr,1.148154e-02_kr,1.202264e-02_kr,1.258925e-02_kr,&
     1.318257e-02_kr,1.380384e-02_kr,1.445440e-02_kr,1.513561e-02_kr,&
     1.584893e-02_kr,1.659587e-02_kr,1.737801e-02_kr,1.819701e-02_kr,&
     1.905461e-02_kr,1.995262e-02_kr,2.089296e-02_kr,2.187762e-02_kr,&
     2.290868e-02_kr,2.398833e-02_kr,2.511886e-02_kr,2.630268e-02_kr,&
     2.754229e-02_kr,2.884032e-02_kr,3.019952e-02_kr,3.162278e-02_kr,&
     3.311311e-02_kr,3.467369e-02_kr,3.630781e-02_kr,3.801894e-02_kr,&
     3.981072e-02_kr,4.168694e-02_kr,4.365158e-02_kr,4.570882e-02_kr,&
     4.786301e-02_kr,5.011872e-02_kr,5.248075e-02_kr,5.495409e-02_kr,&
     5.754399e-02_kr,6.025596e-02_kr,6.309573e-02_kr,6.606934e-02_kr,&
     6.918310e-02_kr,7.244360e-02_kr,7.585776e-02_kr,7.943282e-02_kr,&
     8.317638e-02_kr,8.709636e-02_kr,9.120108e-02_kr,9.549926e-02_kr,&
     1.000000e-01_kr,1.047129e-01_kr,1.096478e-01_kr,1.148154e-01_kr,&
     1.202264e-01_kr,1.258925e-01_kr,1.318257e-01_kr,1.380384e-01_kr,&
     1.445440e-01_kr,1.513561e-01_kr,1.584893e-01_kr,1.659587e-01_kr,&
     1.737801e-01_kr,1.819701e-01_kr,1.905461e-01_kr,1.995262e-01_kr,&
     2.089296e-01_kr,2.187762e-01_kr,2.290868e-01_kr,2.398833e-01_kr,&
     2.511886e-01_kr,2.630268e-01_kr,2.754229e-01_kr,2.884032e-01_kr,&
     3.019952e-01_kr,3.162278e-01_kr,3.311311e-01_kr,3.467369e-01_kr,&
     3.630781e-01_kr,3.801894e-01_kr,3.981072e-01_kr,4.168694e-01_kr,&
     4.365158e-01_kr,4.570882e-01_kr,4.786301e-01_kr,5.011872e-01_kr,&
     5.248075e-01_kr,5.500000e-01_kr,5.750000e-01_kr,6.000000e-01_kr,&
     6.250000e-01_kr,6.500000e-01_kr,6.750000e-01_kr,7.000000e-01_kr,&
     7.250000e-01_kr,7.500000e-01_kr,7.750000e-01_kr,8.000000e-01_kr,&
     8.250000e-01_kr,8.500000e-01_kr,8.750000e-01_kr,9.000000e-01_kr,&
     9.250000e-01_kr,9.500000e-01_kr,9.750000e-01_kr,1.000000e+00_kr,&
     1.025000e+00_kr,1.050000e+00_kr,1.075000e+00_kr,1.100000e+00_kr,&
     1.125000e+00_kr,1.150000e+00_kr,1.175000e+00_kr,1.200000e+00_kr,&
     1.225000e+00_kr,1.250000e+00_kr,1.275000e+00_kr,1.300000e+00_kr,&
     1.325000e+00_kr,1.350000e+00_kr,1.375000e+00_kr,1.400000e+00_kr,&
     1.425000e+00_kr,1.450000e+00_kr,1.475000e+00_kr,1.500000e+00_kr,&
     1.525000e+00_kr,1.550000e+00_kr,1.575000e+00_kr,1.600000e+00_kr,&
     1.625000e+00_kr,1.650000e+00_kr,1.675000e+00_kr,1.700000e+00_kr,&
     1.725000e+00_kr,1.750000e+00_kr,1.775000e+00_kr,1.800000e+00_kr,&
     1.825000e+00_kr,1.850000e+00_kr,1.875000e+00_kr,1.900000e+00_kr,&
     1.925000e+00_kr,1.950000e+00_kr,1.975000e+00_kr,2.000000e+00_kr,&
     2.025000e+00_kr,2.050000e+00_kr,2.075000e+00_kr,2.100000e+00_kr,&
     2.125000e+00_kr,2.150000e+00_kr,2.175000e+00_kr,2.200000e+00_kr,&
     2.225000e+00_kr,2.250000e+00_kr,2.275000e+00_kr,2.300000e+00_kr,&
     2.325000e+00_kr,2.350000e+00_kr,2.375000e+00_kr,2.400000e+00_kr,&
     2.425000e+00_kr,2.450000e+00_kr,2.475000e+00_kr,2.500000e+00_kr,&
     2.525000e+00_kr,2.550000e+00_kr,2.575000e+00_kr,2.600000e+00_kr,&
     2.625000e+00_kr,2.650000e+00_kr,2.675000e+00_kr,2.700000e+00_kr,&
     2.725000e+00_kr,2.750000e+00_kr,2.775000e+00_kr,2.800000e+00_kr,&
     2.825000e+00_kr,2.850000e+00_kr,2.875000e+00_kr,2.900000e+00_kr,&
     2.925000e+00_kr,2.950000e+00_kr,2.975000e+00_kr,3.000000e+00_kr,&
     3.025000e+00_kr,3.050000e+00_kr,3.075000e+00_kr,3.100000e+00_kr,&
     3.125000e+00_kr,3.150000e+00_kr,3.175000e+00_kr,3.200000e+00_kr,&
     3.225000e+00_kr,3.250000e+00_kr,3.275000e+00_kr,3.300000e+00_kr,&
     3.325000e+00_kr,3.350000e+00_kr,3.375000e+00_kr,3.400000e+00_kr,&
     3.425000e+00_kr,3.450000e+00_kr,3.475000e+00_kr,3.500000e+00_kr,&
     3.525000e+00_kr,3.550000e+00_kr,3.575000e+00_kr,3.600000e+00_kr,&
     3.625000e+00_kr,3.650000e+00_kr,3.675000e+00_kr,3.700000e+00_kr,&
     3.725000e+00_kr,3.750000e+00_kr,3.775000e+00_kr,3.800000e+00_kr,&
     3.825000e+00_kr,3.850000e+00_kr,3.875000e+00_kr,3.900000e+00_kr,&
     3.925000e+00_kr,3.950000e+00_kr,3.975000e+00_kr,4.000000e+00_kr,&
     4.025000e+00_kr,4.050000e+00_kr,4.075000e+00_kr,4.100000e+00_kr,&
     4.125000e+00_kr,4.150000e+00_kr,4.175000e+00_kr,4.200000e+00_kr,&
     4.225000e+00_kr,4.250000e+00_kr,4.275000e+00_kr,4.300000e+00_kr,&
     4.325000e+00_kr,4.350000e+00_kr,4.375000e+00_kr,4.400000e+00_kr,&
     4.425000e+00_kr,4.450000e+00_kr,4.475000e+00_kr,4.500000e+00_kr,&
     4.525000e+00_kr,4.550000e+00_kr,4.575000e+00_kr,4.600000e+00_kr,&
     4.625000e+00_kr,4.650000e+00_kr,4.675000e+00_kr,4.700000e+00_kr,&
     4.725000e+00_kr,4.750000e+00_kr,4.775000e+00_kr,4.800000e+00_kr,&
     4.825000e+00_kr,4.850000e+00_kr,4.875000e+00_kr,4.900000e+00_kr,&
     4.925000e+00_kr,4.950000e+00_kr,4.975000e+00_kr,5.000000e+00_kr,&
     5.025000e+00_kr,5.050000e+00_kr,5.075000e+00_kr,5.100000e+00_kr,&
     5.125000e+00_kr,5.150000e+00_kr,5.175000e+00_kr,5.200000e+00_kr,&
     5.225000e+00_kr,5.250000e+00_kr,5.275000e+00_kr,5.300000e+00_kr,&
     5.325000e+00_kr,5.350000e+00_kr,5.375000e+00_kr,5.400000e+00_kr,&
     5.425000e+00_kr,5.450000e+00_kr,5.475000e+00_kr,5.500000e+00_kr,&
     5.525000e+00_kr,5.550000e+00_kr,5.575000e+00_kr,5.600000e+00_kr,&
     5.625000e+00_kr,5.650000e+00_kr,5.675000e+00_kr,5.700000e+00_kr,&
     5.725000e+00_kr,5.750000e+00_kr,5.775000e+00_kr,5.800000e+00_kr,&
     5.825000e+00_kr,5.850000e+00_kr,5.875000e+00_kr,5.900000e+00_kr,&
     5.925000e+00_kr,5.950000e+00_kr,5.975000e+00_kr,6.000000e+00_kr,&
     6.025000e+00_kr,6.050000e+00_kr,6.075000e+00_kr,6.100000e+00_kr,&
     6.125000e+00_kr,6.150000e+00_kr,6.175000e+00_kr,6.200000e+00_kr,&
     6.225000e+00_kr,6.250000e+00_kr,6.275000e+00_kr,6.300000e+00_kr,&
     6.325000e+00_kr,6.350000e+00_kr,6.375000e+00_kr,6.400000e+00_kr,&
     6.425000e+00_kr,6.450000e+00_kr,6.475000e+00_kr,6.500000e+00_kr,&
     6.525000e+00_kr,6.550000e+00_kr,6.575000e+00_kr,6.600000e+00_kr,&
     6.625000e+00_kr,6.650000e+00_kr,6.675000e+00_kr,6.700000e+00_kr,&
     6.725000e+00_kr,6.750000e+00_kr,6.775000e+00_kr,6.800000e+00_kr,&
     6.825000e+00_kr,6.850000e+00_kr,6.875000e+00_kr,6.900000e+00_kr,&
     6.925000e+00_kr,6.950000e+00_kr,6.975000e+00_kr,7.000000e+00_kr,&
     7.025000e+00_kr,7.050000e+00_kr,7.075000e+00_kr,7.100000e+00_kr,&
     7.125000e+00_kr,7.150000e+00_kr,7.175000e+00_kr,7.200000e+00_kr,&
     7.225000e+00_kr,7.250000e+00_kr,7.275000e+00_kr,7.300000e+00_kr,&
     7.325000e+00_kr,7.350000e+00_kr,7.375000e+00_kr,7.400000e+00_kr,&
     7.425000e+00_kr,7.450000e+00_kr,7.475000e+00_kr,7.500000e+00_kr,&
     7.525000e+00_kr,7.550000e+00_kr,7.575000e+00_kr,7.600000e+00_kr,&
     7.625000e+00_kr,7.650000e+00_kr,7.675000e+00_kr,7.700000e+00_kr,&
     7.725000e+00_kr,7.750000e+00_kr,7.775000e+00_kr,7.800000e+00_kr,&
     7.825000e+00_kr,7.850000e+00_kr,7.875000e+00_kr,7.900000e+00_kr,&
     7.925000e+00_kr,7.950000e+00_kr,7.975000e+00_kr,8.000000e+00_kr,&
     8.025000e+00_kr,8.050000e+00_kr,8.075000e+00_kr,8.100000e+00_kr,&
     8.125000e+00_kr,8.150000e+00_kr,8.175000e+00_kr,8.200000e+00_kr,&
     8.225000e+00_kr,8.250000e+00_kr,8.275000e+00_kr,8.300000e+00_kr,&
     8.325000e+00_kr,8.350000e+00_kr,8.375000e+00_kr,8.400000e+00_kr,&
     8.425000e+00_kr,8.450000e+00_kr,8.475000e+00_kr,8.500000e+00_kr,&
     8.525000e+00_kr,8.550000e+00_kr,8.575000e+00_kr,8.600000e+00_kr,&
     8.625000e+00_kr,8.650000e+00_kr,8.675000e+00_kr,8.700000e+00_kr,&
     8.725000e+00_kr,8.750000e+00_kr,8.775000e+00_kr,8.800000e+00_kr,&
     8.825000e+00_kr,8.850000e+00_kr,8.875000e+00_kr,8.900000e+00_kr,&
     8.925000e+00_kr,8.950000e+00_kr,8.975000e+00_kr,9.000000e+00_kr,&
     9.025000e+00_kr,9.050000e+00_kr,9.075000e+00_kr,9.100000e+00_kr,&
     9.125000e+00_kr,9.150000e+00_kr,9.175000e+00_kr,9.200000e+00_kr,&
     9.225000e+00_kr,9.250000e+00_kr,9.275000e+00_kr,9.300000e+00_kr,&
     9.325000e+00_kr,9.350000e+00_kr,9.375000e+00_kr,9.400000e+00_kr,&
     9.425000e+00_kr,9.450000e+00_kr,9.475000e+00_kr,9.500000e+00_kr,&
     9.525000e+00_kr,9.550000e+00_kr,9.575000e+00_kr,9.600000e+00_kr,&
     9.625000e+00_kr,9.650000e+00_kr,9.675000e+00_kr,9.700000e+00_kr,&
     9.725000e+00_kr,9.750000e+00_kr,9.775000e+00_kr,9.800000e+00_kr,&
     9.825000e+00_kr,9.850000e+00_kr,9.875000e+00_kr,9.900000e+00_kr,&
     9.925000e+00_kr,9.950000e+00_kr,9.975000e+00_kr,1.000000e+01_kr,&
     1.047129e+01_kr,1.096478e+01_kr,1.148154e+01_kr,1.202264e+01_kr,&
     1.258925e+01_kr,1.318257e+01_kr,1.380384e+01_kr,1.445440e+01_kr,&
     1.513561e+01_kr,1.584893e+01_kr,1.659587e+01_kr,1.737801e+01_kr,&
     1.819701e+01_kr,1.905461e+01_kr,1.995262e+01_kr,2.089296e+01_kr,&
     2.187762e+01_kr,2.290868e+01_kr,2.398833e+01_kr,2.511886e+01_kr,&
     2.630268e+01_kr,2.754229e+01_kr,2.884032e+01_kr,3.019952e+01_kr,&
     3.162278e+01_kr,3.311311e+01_kr,3.467369e+01_kr,3.630781e+01_kr,&
     3.801894e+01_kr,3.981072e+01_kr,4.168694e+01_kr,4.365158e+01_kr,&
     4.570882e+01_kr,4.786301e+01_kr,5.011872e+01_kr,5.248075e+01_kr,&
     5.495409e+01_kr,5.754399e+01_kr,6.025596e+01_kr,6.309573e+01_kr,&
     6.606934e+01_kr,6.918310e+01_kr,7.244360e+01_kr,7.585776e+01_kr,&
     7.943282e+01_kr,8.317638e+01_kr,8.709636e+01_kr,9.120108e+01_kr,&
     9.549926e+01_kr,1.000000e+02_kr,1.047129e+02_kr,1.096478e+02_kr,&
     1.148154e+02_kr,1.202264e+02_kr,1.258925e+02_kr,1.318257e+02_kr,&
     1.380384e+02_kr,1.445440e+02_kr,1.513561e+02_kr,1.584893e+02_kr,&
     1.659587e+02_kr,1.737801e+02_kr,1.819701e+02_kr,1.905461e+02_kr,&
     1.995262e+02_kr,2.089296e+02_kr,2.187762e+02_kr,2.290868e+02_kr,&
     2.398833e+02_kr,2.511886e+02_kr,2.630268e+02_kr,2.754229e+02_kr,&
     2.884032e+02_kr,3.019952e+02_kr,3.162278e+02_kr,3.311311e+02_kr,&
     3.467369e+02_kr,3.630781e+02_kr,3.801894e+02_kr,3.981072e+02_kr,&
     4.168694e+02_kr,4.365158e+02_kr,4.570882e+02_kr,4.786301e+02_kr,&
     5.011872e+02_kr,5.248075e+02_kr,5.495409e+02_kr,5.754399e+02_kr,&
     6.025596e+02_kr,6.309573e+02_kr,6.606934e+02_kr,6.918310e+02_kr,&
     7.244360e+02_kr,7.585776e+02_kr,7.943282e+02_kr,8.317638e+02_kr,&
     8.709636e+02_kr,9.120108e+02_kr,9.549926e+02_kr,1.000000e+03_kr,&
     1.047129e+03_kr,1.096478e+03_kr,1.148154e+03_kr,1.202264e+03_kr,&
     1.258925e+03_kr,1.318257e+03_kr,1.380384e+03_kr,1.445440e+03_kr,&
     1.513561e+03_kr,1.584893e+03_kr,1.659587e+03_kr,1.737801e+03_kr,&
     1.819701e+03_kr,1.905461e+03_kr,1.995262e+03_kr,2.089296e+03_kr,&
     2.187762e+03_kr,2.290868e+03_kr,2.398833e+03_kr,2.511886e+03_kr,&
     2.630268e+03_kr,2.754229e+03_kr,2.884032e+03_kr,3.019952e+03_kr,&
     3.162278e+03_kr,3.311311e+03_kr,3.467369e+03_kr,3.630781e+03_kr,&
     3.801894e+03_kr,3.981072e+03_kr,4.168694e+03_kr,4.365158e+03_kr,&
     4.570882e+03_kr,4.786301e+03_kr,5.011872e+03_kr,5.248075e+03_kr,&
     5.495409e+03_kr,5.754399e+03_kr,6.025596e+03_kr,6.309573e+03_kr,&
     6.606934e+03_kr,6.918310e+03_kr,7.244360e+03_kr,7.585776e+03_kr,&
     7.943282e+03_kr,8.317638e+03_kr,8.709636e+03_kr,9.120108e+03_kr,&
     9.549926e+03_kr,1.000000e+04_kr,1.047129e+04_kr,1.096478e+04_kr,&
     1.148154e+04_kr,1.202264e+04_kr,1.258925e+04_kr,1.318257e+04_kr,&
     1.380384e+04_kr,1.445440e+04_kr,1.513561e+04_kr,1.584893e+04_kr,&
     1.659587e+04_kr,1.737801e+04_kr,1.819701e+04_kr,1.905461e+04_kr,&
     1.995262e+04_kr,2.089296e+04_kr,2.187762e+04_kr,2.290868e+04_kr,&
     2.398833e+04_kr,2.511886e+04_kr,2.630268e+04_kr,2.754229e+04_kr,&
     2.884032e+04_kr,3.019952e+04_kr,3.162278e+04_kr,3.311311e+04_kr,&
     3.467369e+04_kr,3.630781e+04_kr,3.801894e+04_kr,3.981072e+04_kr,&
     4.168694e+04_kr,4.365158e+04_kr,4.570882e+04_kr,4.786301e+04_kr,&
     5.011872e+04_kr,5.248075e+04_kr,5.495409e+04_kr,5.754399e+04_kr,&
     6.025596e+04_kr,6.309573e+04_kr,6.606934e+04_kr,6.918310e+04_kr,&
     7.244360e+04_kr,7.585776e+04_kr,7.943282e+04_kr,8.317638e+04_kr,&
     8.709636e+04_kr,9.120108e+04_kr,9.549926e+04_kr,1.000000e+05_kr,&
     1.047129e+05_kr,1.096478e+05_kr,1.148154e+05_kr,1.202264e+05_kr,&
     1.258925e+05_kr,1.318257e+05_kr,1.380384e+05_kr,1.445440e+05_kr,&
     1.513561e+05_kr,1.584893e+05_kr,1.659587e+05_kr,1.737801e+05_kr,&
     1.819701e+05_kr,1.905461e+05_kr,1.995262e+05_kr,2.089296e+05_kr,&
     2.187762e+05_kr,2.290868e+05_kr,2.398833e+05_kr,2.511886e+05_kr,&
     2.630268e+05_kr,2.754229e+05_kr,2.884032e+05_kr,3.019952e+05_kr,&
     3.162278e+05_kr,3.311311e+05_kr,3.467369e+05_kr,3.630781e+05_kr,&
     3.801894e+05_kr,3.981072e+05_kr,4.168694e+05_kr,4.365158e+05_kr,&
     4.570882e+05_kr,4.786301e+05_kr,5.011872e+05_kr,5.248075e+05_kr,&
     5.495409e+05_kr,5.754399e+05_kr,6.025596e+05_kr,6.309573e+05_kr,&
     6.606934e+05_kr,6.918310e+05_kr,7.244360e+05_kr,7.585776e+05_kr,&
     7.943282e+05_kr,8.317638e+05_kr,8.709636e+05_kr,9.120108e+05_kr,&
     9.549926e+05_kr,1.000000e+06_kr,1.047129e+06_kr,1.096478e+06_kr,&
     1.148154e+06_kr,1.202264e+06_kr,1.258925e+06_kr,1.318257e+06_kr,&
     1.380384e+06_kr,1.445440e+06_kr,1.513561e+06_kr,1.584893e+06_kr,&
     1.659587e+06_kr,1.737801e+06_kr,1.819701e+06_kr,1.905461e+06_kr,&
     1.995262e+06_kr,2.089296e+06_kr,2.187762e+06_kr,2.290868e+06_kr,&
     2.398833e+06_kr,2.511886e+06_kr,2.630268e+06_kr,2.754229e+06_kr,&
     2.884032e+06_kr,3.019952e+06_kr,3.162278e+06_kr,3.311311e+06_kr,&
     3.467369e+06_kr,3.630781e+06_kr,3.801894e+06_kr,3.981072e+06_kr,&
     4.168694e+06_kr,4.365158e+06_kr,4.570882e+06_kr,4.786301e+06_kr,&
     5.000000e+06_kr,5.200000e+06_kr,5.400000e+06_kr,5.600000e+06_kr,&
     5.800000e+06_kr,6.000000e+06_kr,6.200000e+06_kr,6.400000e+06_kr,&
     6.600000e+06_kr,6.800000e+06_kr,7.000000e+06_kr,7.200000e+06_kr,&
     7.400000e+06_kr,7.600000e+06_kr,7.800000e+06_kr,8.000000e+06_kr,&
     8.200000e+06_kr,8.400000e+06_kr,8.600000e+06_kr,8.800000e+06_kr,&
     9.000000e+06_kr,9.200000e+06_kr,9.400000e+06_kr,9.600000e+06_kr,&
     9.800000e+06_kr,1.000000e+07_kr,1.020000e+07_kr,1.040000e+07_kr,&
     1.060000e+07_kr,1.080000e+07_kr,1.100000e+07_kr,1.120000e+07_kr,&
     1.140000e+07_kr,1.160000e+07_kr,1.180000e+07_kr,1.200000e+07_kr,&
     1.220000e+07_kr,1.240000e+07_kr,1.260000e+07_kr,1.280000e+07_kr,&
     1.300000e+07_kr,1.320000e+07_kr,1.340000e+07_kr,1.360000e+07_kr,&
     1.380000e+07_kr,1.400000e+07_kr,1.420000e+07_kr,1.440000e+07_kr,&
     1.460000e+07_kr,1.480000e+07_kr,1.500000e+07_kr,1.520000e+07_kr,&
     1.540000e+07_kr,1.560000e+07_kr,1.580000e+07_kr,1.600000e+07_kr,&
     1.620000e+07_kr,1.640000e+07_kr,1.660000e+07_kr,1.680000e+07_kr,&
     1.700000e+07_kr,1.720000e+07_kr,1.740000e+07_kr,1.760000e+07_kr,&
     1.780000e+07_kr,1.800000e+07_kr,1.820000e+07_kr,1.840000e+07_kr,&
     1.860000e+07_kr,1.880000e+07_kr,1.900000e+07_kr,1.920000e+07_kr,&
     1.940000e+07_kr,1.960000e+07_kr,1.980000e+07_kr,2.000000e+07_kr,&
     2.020000e+07_kr,2.040000e+07_kr,2.060000e+07_kr,2.080000e+07_kr,&
     2.100000e+07_kr,2.120000e+07_kr,2.140000e+07_kr,2.160000e+07_kr,&
     2.180000e+07_kr,2.200000e+07_kr,2.220000e+07_kr,2.240000e+07_kr,&
     2.260000e+07_kr,2.280000e+07_kr,2.300000e+07_kr,2.320000e+07_kr,&
     2.340000e+07_kr,2.360000e+07_kr,2.380000e+07_kr,2.400000e+07_kr,&
     2.420000e+07_kr,2.440000e+07_kr,2.460000e+07_kr,2.480000e+07_kr,&
     2.500000e+07_kr,2.520000e+07_kr,2.540000e+07_kr,2.560000e+07_kr,&
     2.580000e+07_kr,2.600000e+07_kr,2.620000e+07_kr,2.640000e+07_kr,&
     2.660000e+07_kr,2.680000e+07_kr,2.700000e+07_kr,2.720000e+07_kr,&
     2.740000e+07_kr,2.760000e+07_kr,2.780000e+07_kr,2.800000e+07_kr,&
     2.820000e+07_kr,2.840000e+07_kr,2.860000e+07_kr,2.880000e+07_kr,&
     2.900000e+07_kr,2.920000e+07_kr,2.940000e+07_kr,2.960000e+07_kr,&
     2.980000e+07_kr,3.000000e+07_kr,3.019952e+07_kr,3.162278e+07_kr,&
     3.311311e+07_kr,3.467369e+07_kr,3.630781e+07_kr,3.801894e+07_kr,&
     3.981072e+07_kr,4.168694e+07_kr,4.365158e+07_kr,4.570882e+07_kr,&
     4.786301e+07_kr,5.011872e+07_kr,5.248075e+07_kr,5.495409e+07_kr,&
     5.754399e+07_kr,6.025596e+07_kr,6.309573e+07_kr,6.606934e+07_kr,&
     6.918310e+07_kr,7.244360e+07_kr,7.585776e+07_kr,7.943282e+07_kr,&
     8.317638e+07_kr,8.709636e+07_kr,9.120108e+07_kr,9.549926e+07_kr,&
     1.000000e+08_kr,1.047129e+08_kr,1.096478e+08_kr,1.148154e+08_kr,&
     1.202264e+08_kr,1.258925e+08_kr,1.318257e+08_kr,1.380384e+08_kr,&
     1.445440e+08_kr,1.513561e+08_kr,1.584893e+08_kr,1.659587e+08_kr,&
     1.737801e+08_kr,1.819701e+08_kr,1.905461e+08_kr,1.995262e+08_kr,&
     2.089296e+08_kr,2.187762e+08_kr,2.290868e+08_kr,2.398833e+08_kr,&
     2.511886e+08_kr,2.630268e+08_kr,2.754229e+08_kr,2.884032e+08_kr,&
     3.019952e+08_kr,3.162278e+08_kr,3.311311e+08_kr,3.467369e+08_kr,&
     3.630781e+08_kr,3.801894e+08_kr,3.981072e+08_kr,4.168694e+08_kr,&
     4.365158e+08_kr,4.570882e+08_kr,4.786301e+08_kr,5.011872e+08_kr,&
     5.248075e+08_kr,5.495409e+08_kr,5.754399e+08_kr,6.025596e+08_kr,&
     6.309573e+08_kr,6.606934e+08_kr,6.918310e+08_kr,7.244360e+08_kr,&
     7.585776e+08_kr,7.943282e+08_kr,8.317638e+08_kr,8.709636e+08_kr,&
     9.120108e+08_kr,9.549926e+08_kr,1.000000e+09_kr/)
   real(kr),dimension(143),parameter::eg33=(/&
     5.000000e+03_kr,1.000000e+04_kr,1.500000e+04_kr,2.000000e+04_kr,&
     2.500000e+04_kr,3.000000e+04_kr,3.500000e+04_kr,4.000000e+04_kr,&
     4.500000e+04_kr,5.000000e+04_kr,5.500000e+04_kr,6.000000e+04_kr,&
     6.500000e+04_kr,7.000000e+04_kr,7.500000e+04_kr,8.000000e+04_kr,&
     8.500000e+04_kr,9.000000e+04_kr,9.500000e+04_kr,1.000000e+05_kr,&
     1.200000e+05_kr,1.400000e+05_kr,1.600000e+05_kr,1.800000e+05_kr,&
     2.000000e+05_kr,2.200000e+05_kr,2.400000e+05_kr,2.600000e+05_kr,&
     2.800000e+05_kr,3.000000e+05_kr,3.250000e+05_kr,3.500000e+05_kr,&
     3.750000e+05_kr,4.000000e+05_kr,4.250000e+05_kr,4.500000e+05_kr,&
     4.750000e+05_kr,5.000000e+05_kr,5.250000e+05_kr,5.500000e+05_kr,&
     5.750000e+05_kr,6.000000e+05_kr,6.250000e+05_kr,6.500000e+05_kr,&
     6.750000e+05_kr,7.000000e+05_kr,7.250000e+05_kr,7.500000e+05_kr,&
     7.750000e+05_kr,8.000000e+05_kr,8.250000e+05_kr,8.500000e+05_kr,&
     8.750000e+05_kr,9.000000e+05_kr,9.250000e+05_kr,9.500000e+05_kr,&
     9.750000e+05_kr,1.000000e+06_kr,1.125000e+06_kr,1.250000e+06_kr,&
     1.375000e+06_kr,1.500000e+06_kr,1.625000e+06_kr,1.750000e+06_kr,&
     1.875000e+06_kr,2.000000e+06_kr,2.125000e+06_kr,2.250000e+06_kr,&
     2.375000e+06_kr,2.500000e+06_kr,2.625000e+06_kr,2.750000e+06_kr,&
     2.875000e+06_kr,3.000000e+06_kr,3.125000e+06_kr,3.250000e+06_kr,&
     3.375000e+06_kr,3.500000e+06_kr,3.625000e+06_kr,3.750000e+06_kr,&
     3.875000e+06_kr,4.000000e+06_kr,4.200000e+06_kr,4.400000e+06_kr,&
     4.600000e+06_kr,4.800000e+06_kr,5.000000e+06_kr,5.200000e+06_kr,&
     5.400000e+06_kr,5.600000e+06_kr,5.800000e+06_kr,6.000000e+06_kr,&
     6.200000e+06_kr,6.400000e+06_kr,6.600000e+06_kr,6.800000e+06_kr,&
     7.000000e+06_kr,7.200000e+06_kr,7.400000e+06_kr,7.600000e+06_kr,&
     7.800000e+06_kr,8.000000e+06_kr,8.200000e+06_kr,8.400000e+06_kr,&
     8.600000e+06_kr,8.800000e+06_kr,9.000000e+06_kr,9.200000e+06_kr,&
     9.400000e+06_kr,9.600000e+06_kr,9.800000e+06_kr,1.000000e+07_kr,&
     1.100000e+07_kr,1.200000e+07_kr,1.300000e+07_kr,1.400000e+07_kr,&
     1.500000e+07_kr,1.600000e+07_kr,1.700000e+07_kr,1.800000e+07_kr,&
     1.900000e+07_kr,2.000000e+07_kr,2.200000e+07_kr,2.400000e+07_kr,&
     2.600000e+07_kr,2.800000e+07_kr,3.000000e+07_kr,3.500000e+07_kr,&
     4.000000e+07_kr,4.500000e+07_kr,5.000000e+07_kr,5.400000e+07_kr,&
     5.500000e+07_kr,6.000000e+07_kr,7.000000e+07_kr,8.000000e+07_kr,&
     9.000000e+07_kr,1.000000e+08_kr,1.200000e+08_kr,1.400000e+08_kr,&
     1.600000e+08_kr,1.800000e+08_kr,2.000000e+08_kr/)
   real(kr),dimension(619),parameter::eg618=(/&
     1.00000000000000e-05_kr,2.56901129797510e-05_kr,4.23558357164050e-05_kr,&
     6.98329672839171e-05_kr,1.15135098557100e-04_kr,1.39000000000000e-04_kr,&
     1.89825685995247e-04_kr,2.43741005558083e-04_kr,3.12969646225607e-04_kr,&
     4.01860980405450e-04_kr,5.15999712815652e-04_kr,6.62556746258873e-04_kr,&
     8.50739702194323e-04_kr,1.09237140060287e-03_kr,1.23781896276759e-03_kr,&
     1.40263264283687e-03_kr,1.58939100945164e-03_kr,1.80101596367844e-03_kr,&
     2.04081845319089e-03_kr,2.31255027322349e-03_kr,2.62046276474246e-03_kr,&
     2.96937332818714e-03_kr,3.36474079341315e-03_kr,3.81275082502696e-03_kr,&
     4.32041269930856e-03_kr,4.89566896683177e-03_kr,5.54751971649269e-03_kr,&
     6.28616338510141e-03_kr,7.12315631555298e-03_kr,8.07159355992206e-03_kr,&
     9.14631375620984e-03_kr,1.03641312841130e-02_kr,1.10325603234355e-02_kr,&
     1.17440993319742e-02_kr,1.25015286638674e-02_kr,1.33078079906897e-02_kr,&
     1.41660878664320e-02_kr,1.50797220383603e-02_kr,1.60522805518561e-02_kr,&
     1.70875637004458e-02_kr,1.81896168755305e-02_kr,1.93627463738410e-02_kr,&
     2.06115362243856e-02_kr,2.19408661006432e-02_kr,2.33559303879934e-02_kr,&
     2.48622584808902e-02_kr,2.64657363890912e-02_kr,2.81726297373683e-02_kr,&
     2.99896082485731e-02_kr,3.19237718057234e-02_kr,3.39826781949507e-02_kr,&
     3.61743726377138e-02_kr,3.85074192276762e-02_kr,4.09909343950883e-02_kr,&
     4.36346225294370e-02_kr,4.64488138995581e-02_kr,4.94445050193864e-02_kr,&
     5.26334016170732e-02_kr,5.60279643753727e-02_kr,5.96414576220314e-02_kr,&
     6.34880011604368e-02_kr,6.75826254430556e-02_kr,7.19413303032538e-02_kr,&
     7.65811474749932e-02_kr,8.15202071447017e-02_kr,8.67778087953710e-02_kr,&
     9.23744966197059e-02_kr,9.83321397970035e-02_kr,1.04674017947447e-01_kr,&
     1.11424912097725e-01_kr,1.18611201513438e-01_kr,1.26260966776645e-01_kr,&
     1.34404099511350e-01_kr,1.43072419185677e-01_kr,1.52000000000000e-01_kr,&
     1.62122290476778e-01_kr,1.72578279879602e-01_kr,1.83708622661412e-01_kr,&
     1.95556810878505e-01_kr,2.08169141583816e-01_kr,2.21594897733660e-01_kr,&
     2.35886540761951e-01_kr,2.51099915574398e-01_kr,2.67294468763689e-01_kr,&
     2.84533480898340e-01_kr,3.02884313792894e-01_kr,3.22418673725673e-01_kr,&
     3.43212891632624e-01_kr,3.65348221372105e-01_kr,3.88911157226101e-01_kr,&
     4.14000000000000e-01_kr,4.40694076191185e-01_kr,4.69116402183442e-01_kr,&
     4.99371810711756e-01_kr,5.31578525442442e-01_kr,5.65862394813206e-01_kr,&
     6.02357383788648e-01_kr,6.41206097331274e-01_kr,6.82560337633487e-01_kr,&
     7.26581697287950e-01_kr,7.73442190714156e-01_kr,8.23324926308510e-01_kr,&
     8.76424821944364e-01_kr,9.32949366617847e-01_kr,9.93119431215624e-01_kr,&
     1.05717013157269e+00_kr,1.13000000000000e+00_kr,1.19793069922005e+00_kr,&
     1.27519059148733e+00_kr,1.35743331870246e+00_kr,1.44498024610924e+00_kr,&
     1.53817346522906e+00_kr,1.63737713059081e+00_kr,1.74297888267274e+00_kr,&
     1.85539136261598e+00_kr,1.97505382462877e+00_kr,2.10243385238185e+00_kr,&
     2.23802918610180e+00_kr,2.38236966750182e+00_kr,2.53601931014967e+00_kr,&
     2.69957850336301e+00_kr,2.87368635824370e+00_kr,3.06000000000000e+00_kr,&
     3.25631325144309e+00_kr,3.46632741266196e+00_kr,3.68988632357374e+00_kr,&
     3.92786354548104e+00_kr,4.18118897955002e+00_kr,4.45085250041942e+00_kr,&
     4.73790782415717e+00_kr,5.04347662567888e+00_kr,5.36875292171691e+00_kr,&
     5.71500773646672e+00_kr,6.08359406814152e+00_kr,6.47595217584221e+00_kr,&
     6.89361520740109e+00_kr,7.33821519019035e+00_kr,7.81148940830449e+00_kr,&
     8.32000000000000e+00_kr,8.85157713916813e+00_kr,9.42245481732848e+00_kr,&
     1.00301509424501e+01_kr,1.06770401003478e+01_kr,1.13656500244640e+01_kr,&
     1.20986714730416e+01_kr,1.28789687433204e+01_kr,1.37095908638408e+01_kr,&
     1.45937835085895e+01_kr,1.55350016795403e+01_kr,1.65369232071503e+01_kr,&
     1.76034631215617e+01_kr,1.87387889506673e+01_kr,1.99473370048166e+01_kr,&
     2.12338297117944e+01_kr,2.26000000000000e+01_kr,2.40610812906042e+01_kr,&
     2.56128877094204e+01_kr,2.72647770435633e+01_kr,2.90232040865040e+01_kr,&
     3.08950399301257e+01_kr,3.28875988136648e+01_kr,3.50086667042598e+01_kr,&
     3.72665317207867e+01_kr,3.96700165198641e+01_kr,4.22285127705753e+01_kr,&
     4.49520178526194e+01_kr,4.78511739212901e+01_kr,5.09373094919281e+01_kr,&
     5.42224837063415e+01_kr,5.77195334541645e+01_kr,6.14000000000000e+01_kr,&
     6.54048000453254e+01_kr,6.96230472348795e+01_kr,7.41133479945056e+01_kr,&
     7.88932482720022e+01_kr,8.39814256315774e+01_kr,8.93977622368364e+01_kr,&
     9.51634225407686e+01_kr,1.01300935986307e+02_kr,1.07834285040617e+02_kr,&
     1.14788998907105e+02_kr,1.22192253281342e+02_kr,1.30072976540676e+02_kr,&
     1.38461962782503e+02_kr,1.47391992152865e+02_kr,1.56897958935589e+02_kr,&
     1.67000000000000e+02_kr,1.77788679457205e+02_kr,1.89255064140519e+02_kr,&
     2.01460967099726e+02_kr,2.14454083165892e+02_kr,2.28285183222401e+02_kr,&
     2.43008312593295e+02_kr,2.58681002226541e+02_kr,2.75364493497472e+02_kr,&
     2.93123977510781e+02_kr,3.12028849836190e+02_kr,3.32152981673137e+02_kr,&
     3.53575008504100e+02_kr,3.76378637364449e+02_kr,4.00652973929511e+02_kr,&
     4.26492870696926e+02_kr,4.54000000000000e+02_kr,4.83279736674251e+02_kr,&
     5.14448601797023e+02_kr,5.47627686010971e+02_kr,5.82946637308688e+02_kr,&
     6.20543465259898e+02_kr,6.60565080286848e+02_kr,7.03167867719981e+02_kr,&
     7.48518298877006e+02_kr,7.96793581553195e+02_kr,8.48182352464692e+02_kr,&
     9.02885414350579e+02_kr,9.61116520613947e+02_kr,1.02310321056796e+03_kr,&
     1.05557999269466e+03_kr,1.08908769855066e+03_kr,1.12365905316802e+03_kr,&
     1.15932782038279e+03_kr,1.19612883581024e+03_kr,1.23500000000000e+03_kr,&
     1.27327251787187e+03_kr,1.31369052626409e+03_kr,1.35539153996700e+03_kr,&
     1.39841628594101e+03_kr,1.44280678395902e+03_kr,1.48860638764470e+03_kr,&
     1.53585982681347e+03_kr,1.58461325115751e+03_kr,1.63491427531748e+03_kr,&
     1.68681202538499e+03_kr,1.74035718688117e+03_kr,1.79560205425833e+03_kr,&
     1.85260058197288e+03_kr,1.91140843717952e+03_kr,1.97208305409813e+03_kr,&
     2.03468369010644e+03_kr,2.09927148361327e+03_kr,2.16590951376885e+03_kr,&
     2.23466286207060e+03_kr,2.30559867592442e+03_kr,2.37878623422368e+03_kr,&
     2.45429701500989e+03_kr,2.53220476528119e+03_kr,2.61258557301668e+03_kr,&
     2.69551794148722e+03_kr,2.78108286592499e+03_kr,2.86936391262682e+03_kr,&
     2.96044730056855e+03_kr,3.05442198561012e+03_kr,3.15137974737356e+03_kr,&
     3.25141527887886e+03_kr,3.35000000000000e+03_kr,3.46111354800741e+03_kr,&
     3.57098108576248e+03_kr,3.68433619353942e+03_kr,3.80128957869464e+03_kr,&
     3.92195546281323e+03_kr,4.04645169326264e+03_kr,4.17489985828732e+03_kr,&
     4.30742540575688e+03_kr,4.44415776568380e+03_kr,4.58523047663021e+03_kr,&
     4.73078131612718e+03_kr,4.88095243523415e+03_kr,5.03589049736953e+03_kr,&
     5.19574682154838e+03_kr,5.36067753016696e+03_kr,5.53084370147834e+03_kr,&
     5.70641152690821e+03_kr,5.88755247336443e+03_kr,6.07444345069879e+03_kr,&
     6.26726698448458e+03_kr,6.46621139427874e+03_kr,6.67147097754267e+03_kr,&
     6.88324619940125e+03_kr,7.10174388842549e+03_kr,7.32717743863004e+03_kr,&
     7.55976701788271e+03_kr,7.79973978292963e+03_kr,8.04733010124613e+03_kr,&
     8.30277977992978e+03_kr,8.56633830185941e+03_kr,8.83826306935050e+03_kr,&
     9.12000000000000e+03_kr,9.40828206378196e+03_kr,9.70693299519909e+03_kr,&
     1.00150641248322e+04_kr,1.03329763864764e+04_kr,1.06609802665909e+04_kr,&
     1.09993961075332e+04_kr,1.13485544204187e+04_kr,1.17087962079117e+04_kr,&
     1.20804732972634e+04_kr,1.24639486839205e+04_kr,1.28595968860421e+04_kr,&
     1.32678043102699e+04_kr,1.36889696291092e+04_kr,1.41235041702888e+04_kr,&
     1.45718323184816e+04_kr,1.50343919297757e+04_kr,1.55116347593038e+04_kr,&
     1.60040269024456e+04_kr,1.65120492500366e+04_kr,1.70361979580257e+04_kr,&
     1.75769849320427e+04_kr,1.81349383273462e+04_kr,1.87106030646422e+04_kr,&
     1.93045413622771e+04_kr,1.99173332853231e+04_kr,2.05495773120946e+04_kr,&
     2.12018909186467e+04_kr,2.18749111818289e+04_kr,2.25692954014803e+04_kr,&
     2.32857217423771e+04_kr,2.40248898965561e+04_kr,2.48000000000000e+04_kr,&
     2.55743621709957e+04_kr,2.63861795709192e+04_kr,2.72237668213834e+04_kr,&
     2.80879419452551e+04_kr,2.89795489322345e+04_kr,2.98994585631306e+04_kr,&
     3.08485692603026e+04_kr,3.18278079650967e+04_kr,3.28381310431359e+04_kr,&
     3.38805252183471e+04_kr,3.49560085366367e+04_kr,3.60656313601573e+04_kr,&
     3.72104773931352e+04_kr,3.83916647402616e+04_kr,3.96103469986807e+04_kr,&
     4.08677143846407e+04_kr,4.21649948959093e+04_kr,4.35034555110877e+04_kr,&
     4.48844034269952e+04_kr,4.63091873353325e+04_kr,4.77791987398702e+04_kr,&
     4.92958733154505e+04_kr,5.08606923101270e+04_kr,5.24751839918138e+04_kr,&
     5.41409251408564e+04_kr,5.58595425899810e+04_kr,5.76327148131282e+04_kr,&
     5.94621735647209e+04_kr,6.13497055709682e+04_kr,6.32971542748575e+04_kr,&
     6.53064216365378e+04_kr,6.76000000000000e+04_kr,6.95183239638479e+04_kr,&
     7.17250724500870e+04_kr,7.40018706527728e+04_kr,7.63509421885996e+04_kr,&
     7.87745812594328e+04_kr,8.12751548929221e+04_kr,8.38551052542408e+04_kr,&
     8.65169520312063e+04_kr,8.92632948951132e+04_kr,9.20968160396814e+04_kr,&
     9.50202828005989e+04_kr,9.80365503582183e+04_kr,1.01148564526046e+05_kr,&
     1.04359364627745e+05_kr,1.07672086465471e+05_kr,1.11089965382423e+05_kr,&
     1.14616339422619e+05_kr,1.18254652590966e+05_kr,1.22008458216826e+05_kr,&
     1.25881422424340e+05_kr,1.29877327712922e+05_kr,1.34000076651408e+05_kr,&
     1.38253695689465e+05_kr,1.42642339089993e+05_kr,1.47170292986351e+05_kr,&
     1.51841979568379e+05_kr,1.56661961401289e+05_kr,1.61634945881659e+05_kr,&
     1.66765789834876e+05_kr,1.72059504258514e+05_kr,1.77521259216285e+05_kr,&
     1.84000000000000e+05_kr,1.88970396775857e+05_kr,1.94968961085980e+05_kr,&
     2.01157940267409e+05_kr,2.07543378736997e+05_kr,2.14131512781987e+05_kr,&
     2.20928776650624e+05_kr,2.27941808836123e+05_kr,2.35177458560091e+05_kr,&
     2.42642792461767e+05_kr,2.50345101499601e+05_kr,2.58291908071908e+05_kr,&
     2.66490973363555e+05_kr,2.74950304925867e+05_kr,2.83678164497131e+05_kr,&
     2.92683076071361e+05_kr,3.03000000000000e+05_kr,3.11559512696998e+05_kr,&
     3.21449473268761e+05_kr,3.31653374889103e+05_kr,3.42181183116660e+05_kr,&
     3.53043179850858e+05_kr,3.64249973373642e+05_kr,3.75812508709979e+05_kr,&
     3.87742078317220e+05_kr,4.00050333113793e+05_kr,4.12749293857975e+05_kr,&
     4.25851362887876e+05_kr,4.39369336234074e+05_kr,4.53316416116763e+05_kr,&
     4.67706223839590e+05_kr,4.82552813092796e+05_kr,5.00000000000000e+05_kr,&
     5.13674795672507e+05_kr,5.29980584033558e+05_kr,5.46803973679148e+05_kr,&
     5.64161395037774e+05_kr,5.82069800095720e+05_kr,6.00546678953079e+05_kr,&
     6.19610076905320e+05_kr,6.39278612067076e+05_kr,6.59571493555382e+05_kr,&
     6.80508540250102e+05_kr,7.02110200149880e+05_kr,7.24397570342515e+05_kr,&
     7.47392417609257e+05_kr,7.71117199683167e+05_kr,7.95595087182277e+05_kr,&
     8.23000000000000e+05_kr,8.46906561847805e+05_kr,8.73790261954204e+05_kr,&
     9.01527342308164e+05_kr,9.30144892106635e+05_kr,9.59670860449985e+05_kr,&
     9.90134083638263e+05_kr,1.02156431333394e+06_kr,1.05399224561864e+06_kr,&
     1.08744955097221e+06_kr,1.12196890520344e+06_kr,1.15758402136263e+06_kr,&
     1.19432968266720e+06_kr,1.23224177647237e+06_kr,1.27135732932036e+06_kr,&
     1.31171454310194e+06_kr,1.35300000000000e+06_kr,1.39631286281399e+06_kr,&
     1.44063659101453e+06_kr,1.48636730538123e+06_kr,1.53354966844928e+06_kr,&
     1.58222976049498e+06_kr,1.63245512453958e+06_kr,1.68427481278184e+06_kr,&
     1.73800000000000e+06_kr,1.79290120550119e+06_kr,1.84981399907304e+06_kr,&
     1.90853339864316e+06_kr,1.96911675204194e+06_kr,2.03162322751532e+06_kr,&
     2.09611387151098e+06_kr,2.16265166829887e+06_kr,2.23200000000000e+06_kr,&
     2.30213071747361e+06_kr,2.37520819095458e+06_kr,2.45060539245526e+06_kr,&
     2.52839595804746e+06_kr,2.60865586126285e+06_kr,2.69146348729184e+06_kr,&
     2.77689970953790e+06_kr,2.86500000000000e+06_kr,2.95599435377371e+06_kr,&
     3.04982768711059e+06_kr,3.14663961018459e+06_kr,3.24652467358350e+06_kr,&
     3.34958042925295e+06_kr,3.45590752576975e+06_kr,3.56560980663947e+06_kr,&
     3.68000000000000e+06_kr,3.79557188183090e+06_kr,3.91605626676799e+06_kr,&
     4.04036523663342e+06_kr,4.10399173096370e+06_kr,4.16862019678508e+06_kr,&
     4.23426641285263e+06_kr,4.30094640640062e+06_kr,4.36867645705557e+06_kr,&
     4.43747310081080e+06_kr,4.50735313406362e+06_kr,4.57833361771614e+06_kr,&
     4.65043188134056e+06_kr,4.72366552741015e+06_kr,4.79805243559677e+06_kr,&
     4.87361076713619e+06_kr,4.95035896926199e+06_kr,5.02831577970941e+06_kr,&
     5.10750023129011e+06_kr,5.18793165653889e+06_kr,5.26962969243371e+06_kr,&
     5.35261428518990e+06_kr,5.43690569513000e+06_kr,5.52252450163020e+06_kr,&
     5.60949160814471e+06_kr,5.69782824730923e+06_kr,5.78755598612484e+06_kr,&
     5.87869673122347e+06_kr,5.97127273421627e+06_kr,6.07000000000000e+06_kr,&
     6.16082127790678e+06_kr,6.25784009604591e+06_kr,6.35638673826052e+06_kr,&
     6.45648526427892e+06_kr,6.55816011271502e+06_kr,6.66143610703488e+06_kr,&
     6.76633846161729e+06_kr,6.87289278790972e+06_kr,6.98112510068126e+06_kr,&
     7.09106182437398e+06_kr,7.20272979955440e+06_kr,7.31615628946642e+06_kr,&
     7.43136898668758e+06_kr,7.54839601989007e+06_kr,7.66726596070820e+06_kr,&
     7.79000000000000e+06_kr,7.91065110850296e+06_kr,8.03522573689061e+06_kr,&
     8.16176213022340e+06_kr,8.29029118180400e+06_kr,8.42084427143382e+06_kr,&
     8.55345327307423e+06_kr,8.68815056262843e+06_kr,8.82496902584595e+06_kr,&
     8.96394206635151e+06_kr,9.10510361380034e+06_kr,9.24848813216205e+06_kr,&
     9.39413062813476e+06_kr,9.54206665969188e+06_kr,9.69233234476344e+06_kr,&
     9.84496437005408e+06_kr,1.00000000000000e+07_kr,1.01250000000000e+07_kr,&
     1.02500000000000e+07_kr,1.03750000000000e+07_kr,1.05000000000000e+07_kr,&
     1.06250000000000e+07_kr,1.07500000000000e+07_kr,1.08750000000000e+07_kr,&
     1.10000000000000e+07_kr,1.11250000000000e+07_kr,1.12500000000000e+07_kr,&
     1.13750000000000e+07_kr,1.15000000000000e+07_kr,1.16250000000000e+07_kr,&
     1.17500000000000e+07_kr,1.18750000000000e+07_kr,1.20000000000000e+07_kr,&
     1.21250000000000e+07_kr,1.22500000000000e+07_kr,1.23750000000000e+07_kr,&
     1.25000000000000e+07_kr,1.26250000000000e+07_kr,1.27500000000000e+07_kr,&
     1.28750000000000e+07_kr,1.30000000000000e+07_kr,1.31250000000000e+07_kr,&
     1.32500000000000e+07_kr,1.33750000000000e+07_kr,1.35000000000000e+07_kr,&
     1.36250000000000e+07_kr,1.37500000000000e+07_kr,1.38750000000000e+07_kr,&
     1.40000000000000e+07_kr,1.41250000000000e+07_kr,1.42500000000000e+07_kr,&
     1.43750000000000e+07_kr,1.45000000000000e+07_kr,1.46250000000000e+07_kr,&
     1.47500000000000e+07_kr,1.48750000000000e+07_kr,1.50000000000000e+07_kr,&
     1.51250000000000e+07_kr,1.52500000000000e+07_kr,1.53750000000000e+07_kr,&
     1.55000000000000e+07_kr,1.56250000000000e+07_kr,1.57500000000000e+07_kr,&
     1.58750000000000e+07_kr,1.60000000000000e+07_kr,1.61250000000000e+07_kr,&
     1.62500000000000e+07_kr,1.63750000000000e+07_kr,1.65000000000000e+07_kr,&
     1.66250000000000e+07_kr,1.67500000000000e+07_kr,1.68750000000000e+07_kr,&
     1.70000000000000e+07_kr,1.71250000000000e+07_kr,1.72500000000000e+07_kr,&
     1.73750000000000e+07_kr,1.75000000000000e+07_kr,1.76250000000000e+07_kr,&
     1.77500000000000e+07_kr,1.78750000000000e+07_kr,1.80000000000000e+07_kr,&
     1.81250000000000e+07_kr,1.82500000000000e+07_kr,1.83750000000000e+07_kr,&
     1.85000000000000e+07_kr,1.86250000000000e+07_kr,1.87500000000000e+07_kr,&
     1.88750000000000e+07_kr,1.90000000000000e+07_kr,1.91250000000000e+07_kr,&
     1.92500000000000e+07_kr,1.93750000000000e+07_kr,1.95000000000000e+07_kr,&
     1.96250000000000e+07_kr,1.97500000000000e+07_kr,1.98750000000000e+07_kr,&
     2.00000000000000e+07_kr/)
   real(kr),parameter::ezero=1.e7_kr
   real(kr),parameter::tenth=0.10e0_kr
   real(kr),parameter::eighth=0.125e0_kr
   real(kr),parameter::quart=0.25e0_kr
   real(kr),parameter::bgam2=27.631021e0_kr
   real(kr),parameter::tgam2=-0.53062825e0_kr
   real(kr),parameter::u187a=-15.5e0_kr
   real(kr),parameter::u187b=-14.125e0_kr
   real(kr),parameter::u187c=-5.875e0_kr
   real(kr),parameter::e187d=2.6058e4_kr
   real(kr),parameter::e187e=6.868e0_kr
   real(kr),parameter::sanda=1.e-4_kr
   real(kr),parameter::sandb=1.e-6_kr
   real(kr),parameter::sandc=2.8e-4_kr
   real(kr),parameter::sandd=1.e6_kr
   real(kr),parameter::sande=1.e5_kr
   real(kr),parameter::sandf=1.e6_kr
   real(kr),parameter::uu80=.6931472e0_kr
   real(kr),parameter::e175=1.284e7_kr

   !--choose option according to ign
   lflag=0

   !--group structure is read in (free format)
   if (ign.eq.1) then
      read(nsysi,*) ngn
      ngp=ngn+1
      allocate(egn(ngp))
      read(nsysi,*) (egn(ig),ig=1,ngp)
      do ig=1,ngn
         if (egn(ig).ge.egn(ig+1)) call error('gengpn',&
           'read-in group structure is out of order.',' ')
      enddo

   !--csewg 239 group structure
   else if (ign.eq.2) then
      ngn=240
      ngp=ngn+1
      allocate(egn(ngp))
      do ig=1,ngp
         egn(ig)=gl2(ig)
      enddo
      lflag=1

   !--lanl 30 group structure
   else if (ign.eq.3) then
      ngn=30
      ngp=ngn+1
      allocate(egn(ngp))
      do ig=1,ngp
         egn(ig)=eg3(ig)
      enddo

   !--anl 27 group structure
   else if (ign.eq.4) then
      ngn=27
      ngp=ngn+1
      allocate(egn(ngp))
      do ig=1,ngp
         egn(ig)=gl4(ig)
      enddo
      lflag=1

   !--rrd 50 group structure
   else if (ign.eq.5) then
      ngn=50
      ngp=ngn+1
      allocate(egn(ngp))
      do ig=1,ngp
         egn(ig)=gl5(ig)
      enddo
      lflag=1

   !--gam-i 68 group structure
   else if (ign.eq.6) then
      ngn=68
      u=-quart
      du=quart
      ngp=ngn+1
      allocate(egn(ngp))
      do ig=1,ngp
         u=u+du
         egn(70-ig)=u
      enddo
      lflag=1

   !--gam-ii 100 group structure
   else if (ign.eq.7) then
      ngn=100
      ngp=ngn+1
      allocate(egn(ngp))
      u=-4*tenth
      du=tenth
      do ig=1,99
         u=u+du
         egn(101-ig)=u
         if (ig.eq.49) du=quart
      enddo
      egn(1)=bgam2
      ! upper limit changed to 17 mev
      egn(101)=tgam2
      lflag=1

   !--laser-thermos 35 group structure
   else if (ign.eq.8) then
      ngn=35
      ngp=ngn+1
      allocate(egn(ngp))
      do ig=1,ngp
         egn(ig)=eg6(ig)
      enddo

   !--epri-cpm 69 group structure
   else if (ign.eq.9) then
      ngn=69
      ngp=ngn+1
      allocate(egn(ngp))
      do ig=1,ngp
         egn(ig)=eg9(ig)
      enddo

   !--lanl 187-group structure
   else if (ign.eq.10) then
      ngn=187
      ngp=ngn+1
      allocate(egn(ngp))
      do ig=1,48
         egn(ig)=eg10a(ig)
      enddo
      u=u187a
      do ig=49,59
         egn(ig)=ezero*exp(u)
         u=u+eighth
      enddo
      egn(60)=e187e
      u=u187b
      do ig=61,126
         egn(ig)=ezero*exp(u)
         u=u+eighth
      enddo
      egn(127)=e187d
      u=u187c
      do ig=128,175
         egn(ig)=ezero*exp(u)
         u=u+eighth
      enddo
      do ig=176,188
         egn(ig)=eg10b(ig-175)
      enddo

   !--lanl 70 group structure
   else if (ign.eq.11) then
      ngn=70
      ngp=ngn+1
      allocate(egn(ngp))
      do ig=1,ngp
         egn(ig)=eg11(ig)
      enddo

   !--sand-ii 620- and 640-group structures
   else if (ign.eq.12.or.ign.eq.15 .or. ign.eq.24) then
      ngn=620
      if (ign.eq.24) ngn=770
      if (ign.eq.15) ngn=640
      ngp=ngn+1
      allocate(egn(ngp))
      egn(1)=sanda
      ! generate the first 45 boundaries
      do i=1,8
         delta=deltl(i)*sandb
         n1=ndelta(i)
         n2=ndelta(i+1)-1
         do n=n1,n2
            egn(n)=egn(n-1)+delta
         enddo
      enddo
      ! correct group 21
      egn(21)=sandc
      ! groups 46 to 450 are multiples of previous groups
      do i=46,450
         egn(i)=egn(i-45)*10
      enddo
      ! groups 451 through 620 have constant spacing of 1.e5
      egn(451)=sandd
      ngp_hold = min(ngp,641)
      do i=452,ngp_hold 
         egn(i)=egn(i-1)+sande
      enddo
      if ( ngp.eq. 771) then
          do i=642,ngp                                                                     
             egn(i)=egn(i-1)+sandf                                                         
          enddo   
      endif

   !--lanl 80-group structure
   else if (ign.eq.13) then
      ngn=80
      u=uu80
      ngp=ngn+1
      allocate(egn(ngp))
      do ig=1,ngp
         egn(82-ig)=ezero*exp(u)
         if (ig.le.ngn) u=u-u80(ig)
      enddo
      egn(81)=2*ezero

   !--eurlib 100-group structure
   else if (ign.eq.14) then
      ngn=100
      ngp=ngn+1
      allocate(egn(ngp))
      egn(101)=-4
      egn(101)=egn(101)/10
      ic=0
      do ig=2,101
         if (ig.eq.ig14(ic+1)) ic=ic+1
         egn(102-ig)=egn(103-ig)+gl14(ic)
      enddo
      lflag=1

   !--vitamin-e 174- and vitamin-j 175-group structures (ornl-5510)
   else if (ign.eq.16.or.ign.eq.17) then
      if (ign.eq.16) ngn=174
      if (ign.eq.17) ngn=175
      ngp=ngn+1
      allocate(egn(ngp))
      do ig=1,84
         egn(ig)=eg15a(ig)
      enddo
      do ig=85,175
         egn(ig)=eg15b(ig-84)
      enddo
      if (ign.ne.16) then
         egn(166)=e175
         do ig=167,176
            egn(ig)=eg15b(ig-85)
         enddo
      endif

   !--xmas 172-group structure
   else if (ign.eq.18) then
      ngn=172
      ngp=ngn+1
      allocate(egn(ngp))
      do ig=1,ngp
         egn(ig)=eg18(174-ig)
      enddo

   !--ecco  33-group structure
   else if (ign.eq.19) then
      ngn=33
      ngp=ngn+1
      allocate(egn(ngp))
      do ig=1,34
         egn(ig)=eg19(ig)
      enddo

   !--ecco 1968-group structure
   else if (ign.eq.20) then
      ngn=1968
      ngp=ngn+1
      allocate(egn(ngp))
      do ig=1,392
         egn(ig)=eg20a(ig)
      enddo
      do ig=1,392
         egn(ig+392)=eg20b(ig)
      enddo
      do ig=1,392
         egn(ig+784)=eg20c(ig)
      enddo
      do ig=1,392
         egn(ig+1176)=eg20d(ig)
      enddo
      do ig=1,392
         egn(ig+1568)=eg20e(ig)
      enddo
      do ig=1,9
         egn(ig+1960)=eg20f(ig)
      enddo

   !--tripoli 315-group structure
   else if (ign.eq.21) then
      ngn=315
      ngp=ngn+1
      allocate(egn(ngp))
      do ig=1,316
         egn(ig)=eg21(ig)
      enddo

   !--xmas lwpc 172-group structure
   else if (ign.eq.22) then
      ngn=172
      ngp=ngn+1
      allocate(egn(ngp))
      do ig=1,173
         egn(ig)=eg22(ig)
      enddo

   !--vit-j lwpc 175-group structure
   else if (ign.eq.23) then
      ngn=175
      ngp=ngn+1
      allocate(egn(ngp))
      do ig=1,176
         egn(ig)=eg23(ig)
      enddo

   !--illegal ign
   else
      call error('gengpn','illegal group structure.',' ')
   endif

   !--convert lethargy grid to energies
   if (lflag.eq.1) then
      ngp=ngn+1
      do ig=1,ngp
         egn(ig)=sigfig(ezero*exp(-egn(ig)),7,0)
      enddo
   endif

   !--display group structure
   if (ign.eq.1) write(nsyso,'(/&
     &'' neutron group structure......read in'')')
   if (ign.eq.2) write(nsyso,'(/&
     &'' neutron group structure......csewg 240 group'')')
   if (ign.eq.3) write(nsyso,'(/&
     &'' neutron group structure......lanl 30 group'')')
   if (ign.eq.4) write(nsyso,'(/&
     &'' neutron group structure......anl 27 group'')')
   if (ign.eq.5) write(nsyso,'(/&
     &'' neutron group structure......rrd 50 group'')')
   if (ign.eq.6) write(nsyso,'(/&
     &'' neutron group structure......gam-i 68 group'')')
   if (ign.eq.7) write(nsyso,'(/&
     &'' neutron group structure......gam-ii 100 group'')')
   if (ign.eq.8) write(nsyso,'(/&
     &'' neutron group structure......laser-thermos 35 group'')')
   if (ign.eq.9) write(nsyso,'(/&
     &'' neutron group structure......epri-cpm 69 group'')')
   if (ign.eq.10) write(nsyso,'(/&
     &'' neutron group structure......lanl 187-group'')')
   if (ign.eq.11) write(nsyso,'(/&
     &'' neutron group structure......lanl 70-group'')')
   if (ign.eq.12) write(nsyso,'(/&
     &'' neutron group structure......sand-ii 620 group'')')
   if (ign.eq.13) write(nsyso,'(/&
     &'' neutron group structure......lanl 80-group'')')
   if (ign.eq.14) write(nsyso,'(/&
     &'' neutron group structure......eurlib 100-group'')')
   if (ign.eq.15) write(nsyso,'(/&
     &'' neutron group structure......sand-iia 640-group'')')
   if (ign.eq.16) write(nsyso,'(/&
     &'' neutron group structure......vitamin-e 174-group'')')
   if (ign.eq.17) write(nsyso,'(/&
     &'' neutron group structure......vitamin-j 175-group'')')
   if (ign.eq.18) write(nsyso,'(/&
     &'' neutron group structure......xmas 172-group'')')
   if (ign.eq.19) write(nsyso,'(/&
     &  '' neutron group structure......ecco  33-group'')')
   if (ign.eq.20) write(nsyso,'(/&
     &  '' neutron group structure......ecco  1968-group'')')
   if (ign.eq.21) write(nsyso,'(/&
     &  '' neutron group structure......tripoli 315-group'')')
   if (ign.eq.22) write(nsyso,'(/&
     &  '' neutron group structure......xmas lwpc 172-group'')')
   if (ign.eq.23) write(nsyso,'(/&
     &  '' neutron group structure......vit-j lwpc 175-group'')')
   if (ign.eq.24) write(nsyso,'(/&                                                     
   &  '' neutron group structure......SAND-IV 770-group'')')   
   do ig=1,ngn
      write(nsyso,'(1x,i5,2x,1p,e12.5,''  - '',e12.5)')&
        ig,egn(ig),egn(ig+1)
   enddo
   return
   end subroutine gengpn

   subroutine gengpg
   !-------------------------------------------------------------------
   ! Generate requested gamma group structure or read in from
   ! the system input file in the form of an ENDF list record
   !
   !    igg     meaning
   !    ---     --------------------------------------
   !     0      none
   !     1      arbitrary structure (read in)
   !     2      CSEWG 94-group structure
   !     3      LANL 12-group structure
   !     4      Steiner 21-group gamma-ray structure (ORNL-TM-2564)
   !     5      Straker 22-group structure
   !     6      LANL 48-group structure
   !     7      LANL 24-group structure
   !     8      VITAMIN-C 36-group structure
   !     9      VITAMIN-E 38-group structure (R. Roussin, Feb 86)
   !    10      VITAMIN-J 42-group structure
   !
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use util   ! provides error
   ! internals
   integer::ngp,ig
   real(kr),dimension(95),parameter::eg2=(/&
     .005e0_kr,.01e0_kr,.015e0_kr,.02e0_kr,.03e0_kr,.035e0_kr,&
     .04e0_kr,.045e0_kr,.055e0_kr,.06e0_kr,.065e0_kr,.075e0_kr,&
     .08e0_kr,.09e0_kr,.1e0_kr,.12e0_kr,.14e0_kr,.15e0_kr,&
     .16e0_kr,.19e0_kr,.22e0_kr,.26e0_kr,.3e0_kr,.325e0_kr,.35e0_kr,&
     .375e0_kr,.4e0_kr,.425e0_kr,.45e0_kr,.5e0_kr,.525e0_kr,&
     .55e0_kr,.575e0_kr,.6e0_kr,.625e0_kr,.65e0_kr,.675e0_kr,&
     .7e0_kr,.75e0_kr,.8e0_kr,.825e0_kr,.865e0_kr,.9e0_kr,1.e0_kr,&
     1.125e0_kr,1.2e0_kr,1.25e0_kr,1.33e0_kr,1.42e0_kr,1.5e0_kr,&
     1.6e0_kr,1.66e0_kr,1.75e0_kr,1.875e0_kr,2.e0_kr,2.166e0_kr,&
     2.333e0_kr,2.5e0_kr,2.666e0_kr,2.833e0_kr,3.e0_kr,3.166e0_kr,&
     3.333e0_kr,3.5e0_kr,3.65e0_kr,3.8e0_kr,3.9e0_kr,4.e0_kr,&
     4.2e0_kr,4.4e0_kr,4.5e0_kr,4.7e0_kr,5.e0_kr,5.2e0_kr,5.4e0_kr,&
     5.5e0_kr,5.75e0_kr,6.e0_kr,6.25e0_kr,6.5e0_kr,6.75e0_kr,7.e0_kr,&
     7.25e0_kr,7.5e0_kr,7.75e0_kr,8.e0_kr,8.5e0_kr,9.e0_kr,9.5e0_kr,&
     10.e0_kr,10.6e0_kr,11.e0_kr,12.e0_kr,14.e0_kr,20.e0_kr/)
   real(kr),dimension(13),parameter::eg3=(/&
     .01e0_kr,.10e0_kr,.50e0_kr,1.0e0_kr,2.0e0_kr,3.0e0_kr,4.0e0_kr,&
     5.0e0_kr,6.0e0_kr,7.0e0_kr,8.0e0_kr,9.0e0_kr,20.0e0_kr/)
   real(kr),dimension(22),parameter::eg4=(/&
     .01e0_kr,.1e0_kr,.2e0_kr,.4e0_kr,1.e0_kr,1.5e0_kr,2.e0_kr,&
     2.5e0_kr,3.e0_kr,3.5e0_kr,4.e0_kr,4.5e0_kr,5.e0_kr,5.5e0_kr,&
     6.e0_kr,6.5e0_kr,7.e0_kr,7.5e0_kr,8.e0_kr,10.e0_kr,12.e0_kr,&
     14.e0_kr/)
   real(kr),dimension(23),parameter::eg5=(/&
     .01e0_kr,.03e0_kr,.06e0_kr,.10e0_kr,.15e0_kr,.30e0_kr,&
     .45e0_kr,.60e0_kr,.80e0_kr,1.0e0_kr,1.33e0_kr,1.66e0_kr,&
     2.0e0_kr,2.5e0_kr,3.0e0_kr,3.5e0_kr,4.0e0_kr,5.0e0_kr,6.0e0_kr,&
     7.0e0_kr,8.0e0_kr,10.0e0_kr,14.0e0_kr/)
   real(kr),dimension(49),parameter::eg6=(/&
     .001e0_kr,.01e0_kr,.02e0_kr,.03e0_kr,.045e0_kr,.06e0_kr,&
     .08e0_kr,.1e0_kr,.15e0_kr,.2e0_kr,.3e0_kr,.4e0_kr,.45e0_kr,&
     .5e0_kr,.525e0_kr,.6e0_kr,.7e0_kr,.8e0_kr,.9e0_kr,1.e0_kr,&
     1.125e0_kr,1.2e0_kr,1.33e0_kr,1.5e0_kr,1.66e0_kr,1.875e0_kr,&
     2.e0_kr,2.333e0_kr,2.5e0_kr,2.666e0_kr,3.e0_kr,3.5e0_kr,&
     4.e0_kr,4.5e0_kr,5.e0_kr,5.5e0_kr,6.e0_kr,6.5e0_kr,7.e0_kr,&
     7.5e0_kr,8.e0_kr,9.e0_kr,10.e0_kr,12.e0_kr,14.e0_kr,17.e0_kr,&
     20.e0_kr,30.e0_kr,50.e0_kr/)
   real(kr),dimension(25),parameter::eg7=(/&
     1.e4_kr,3.e4_kr,6.e4_kr,1.e5_kr,2.e5_kr,3.e5_kr,5.e5_kr,&
     5.25e5_kr,7.5e5_kr,1.e6_kr,1.33e6_kr,1.66e6_kr,2.e6_kr,&
     2.5e6_kr,3.e6_kr,4.e6_kr,5.e6_kr,6.e6_kr,7.e6_kr,8.e6_kr,&
     9.e6_kr,1.e7_kr,1.2e7_kr,1.7e7_kr,3.e7_kr/)
   real(kr),dimension(39),parameter::eg8=(/&
     .01e0_kr,.02e0_kr,.03e0_kr,.045e0_kr,.06e0_kr,.07e0_kr,&
     .075e0_kr,.10e0_kr,.15e0_kr,.20e0_kr,.30e0_kr,.40e0_kr,&
     .45e0_kr,.51e0_kr,.512e0_kr,.60e0_kr,.70e0_kr,.80e0_kr,1.0e0_kr,&
     1.33e0_kr,1.5e0_kr,1.66e0_kr,2.0e0_kr,2.5e0_kr,3.0e0_kr,&
     3.5e0_kr,4.0e0_kr,4.5e0_kr,5.0e0_kr,5.5e0_kr,6.0e0_kr,&
     6.5e0_kr,7.0e0_kr,7.5e0_kr,8.0e0_kr,10.e0_kr,12.e0_kr,&
     14.e0_kr,20.e0_kr/)
   real(kr),dimension(43),parameter::eg10=(/&
     1.0e3_kr,1.0e4_kr,2.0e4_kr,3.0e4_kr,4.5e4_kr,6.0e4_kr,7.0e4_kr,&
     7.5e4_kr,1.0e5_kr,1.50e5_kr,2.00e5_kr,3.00e5_kr,4.00e5_kr,&
     4.50e5_kr,5.10e5_kr,5.12e5_kr,6.00e5_kr,7.00e5_kr,8.00e5_kr,&
     1.00e6_kr,1.33e6_kr,1.34e6_kr,1.50e6_kr,1.66e6_kr,2.00e6_kr,&
     2.50e6_kr,3.00e6_kr,3.50e6_kr,4.00e6_kr,4.50e6_kr,5.00e6_kr,&
     5.50e6_kr,6.00e6_kr,6.50e6_kr,7.00e6_kr,7.50e6_kr,8.00e6_kr,&
     1.00e7_kr,1.20e7_kr,1.40e7_kr,2.00e7_kr,3.00e7_kr,5.00e7_kr/)
   real(kr),parameter::emev=1.e6_kr

   !--select structure
   if (igg.eq.0) then
      ngg=0
      allocate(egg(1))
      egg(1)=0
      return

   !--group structure is read in.
   else if (igg.eq.1) then
      read(nsysi,*) ngg
      ngp=ngg+1
      allocate(egg(ngp))
      read(nsysi,*) (egg(ig),ig=1,ngp)

   !--csewg 94 group structure
   else if (igg.eq.2) then
      ngg=94
      ngp=ngg+1
      allocate(egg(ngp))
      do ig=1,ngp
         egg(ig)=eg2(ig)*emev
      enddo

   !--lanl 12 group structure
   else if (igg.eq.3) then
      ngg=12
      ngp=ngg+1
      allocate(egg(ngp))
      do ig=1,ngp
         egg(ig)=eg3(ig)*emev
      enddo

   !--steiner 21-group gamma structure (ornl-tm-2564)
   else if (igg.eq.4) then
      ngg=21
      ngp=ngg+1
      allocate(egg(ngp))
      do ig=1,ngp
         egg(ig)=eg4(ig)*emev
      enddo

   !--straker 22 group structure
   else if (igg.eq.5) then
      ngg=22
      ngp=ngg+1
      allocate(egg(ngp))
      do ig=1,ngp
         egg(ig)=eg5(ig)*emev
      enddo

   !--lanl 48-group structure
   else if (igg.eq.6) then
      ngg=48
      ngp=ngg+1
      allocate(egg(ngp))
      do ig=1,ngp
         egg(ig)=eg6(ig)*emev
      enddo

   !--lanl 24-group structure
   else if (igg.eq.7) then
      ngg=24
      ngp=ngg+1
      allocate(egg(ngp))
      do ig=1,ngp
         egg(ig)=eg7(ig)
      enddo

   !--vitamin-series 36- and 38-group structures
   else if (igg.eq.8.or.igg.eq.9) then
      ngg=38
      if (igg.eq.8) ngg=36
      ngp=ngg+1
      allocate(egg(ngp))
      do ig=1,ngp
         egg(ig)=eg8(ig)*emev
      enddo
      if (igg.ne.9) then
         ! remove group bounds eg8(7) and eg8(39) if igg=8
         do ig=7,ngp
            egg(ig)=eg8(ig+1)*emev
         enddo
      endif

   !--vitamin-j 42-group structure
   else if (igg.eq.10) then
      ngg=42
      ngp=ngg+1
      allocate(egg(ngp))
      do ig=1,ngp
         egg(ig)=eg10(ig)
      enddo

   !--illegal igg
   else
      call error('gengpg','illegal group structure.',' ')
   endif

   !--display group structure
   if (igg.eq.1) write(nsyso,'(/&
     &'' gamma group structure......read in'')')
   if (igg.eq.2) write(nsyso,'(/&
     &'' gamma group structure......csewg 94 group'')')
   if (igg.eq.3) write(nsyso,'(/&
     &'' gamma group structure......lanl 12 group'')')
   if (igg.eq.4) write(nsyso,'(/&
     &'' gamma group structure......steiner 21-group'')')
   if (igg.eq.5) write(nsyso,'(/&
     &'' gamma group structure......straker 22 group'')')
   if (igg.eq.6) write(nsyso,'(/&
     &'' gamma group structure......lanl 48-group'')')
   if (igg.eq.7) write(nsyso,'(/&
     &'' gamma group structure......lanl 24-group'')')
   if (igg.eq.8) write(nsyso,'(/&
     &'' gamma group structure......vitamin-c 36-group'')')
   if (igg.eq.9) write(nsyso,'(/&
     &'' gamma group structure......vitamin-e 38-group'')')
   if (igg.eq.10) write(nsyso,'(/&
     &'' gamma group structure......vitamin-j 42-group'')')
   do ig=1,ngg
      write(nsyso,'(1x,i5,2x,1p,e12.5,''  - '',e12.5)')&
        ig,egg(ig),egg(ig+1)
   enddo
   return
   end subroutine gengpg

   subroutine genwtf
   !-------------------------------------------------------------------
   ! Set up calculation of weight functions or read in arbitary
   ! function in the form of an ENDF tab1 record or
   ! read in parameters for an analytic weight function.
   ! Negative values trigger flux calculation.
   !
   !    iwt     meaning
   !    ---     -------
   !     1      read in
   !     2      constant
   !     3      1/e
   !     4      1/e + fission spectrum + thermal Maxwellian
   !     5      EPRI-CELL LWR
   !     6      (thermal) -- (1/e) -- (fission + fusion)
   !     7      same with t-dep thermal part
   !     8      thermal--1/e--fast reactor--fission + fusion
   !     9      claw weight function
   !    10      claw with t-dependent thermal part
   !    11      VITAMIN-E weight function (ORNL-5505)
   !    12      VITAMIN-E with t-dep thermal part
   !     0      read resonance flux for flux calculator
   !    -n      read flux calc. input, and set up iwt=n.
   !
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use util   ! provides error
   ! internals
   integer::iwtt,i,nr,np,ntmp,iw
   real(kr)::ehi,eb,tb,ec,tc,ab,ac
   real(kr),dimension(:),allocatable::tmp
   real(kr),dimension(92),parameter::w1=(/&
     0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,1.e0_kr,88.e0_kr,88.e0_kr,&
     5.e0_kr,1.e-5_kr,5.25e-4_kr,.009e0_kr,.355e0_kr,.016e0_kr,&
     .552e0_kr,.024e0_kr,.712e0_kr,.029e0_kr,.785e0_kr,.033e0_kr,&
     .829e0_kr,.043e0_kr,.898e0_kr,.05e0_kr,.918e0_kr,.054e0_kr,&
     .921e0_kr,.059e0_kr,.918e0_kr,.07e0_kr,.892e0_kr,.09e0_kr,&
     .799e0_kr,.112e0_kr,.686e0_kr,.14e0_kr,.52e0_kr,.17e0_kr,&
     .383e0_kr,.21e0_kr,.252e0_kr,.3e0_kr,.108e0_kr,.4e0_kr,&
     .0687e0_kr,.49e0_kr,.051e0_kr,.57e0_kr,.0437e0_kr,.6e0_kr,&
     .0413e0_kr,1.e0_kr,.024914e0_kr,1.01e3_kr,3.7829e-5_kr,2.e4_kr,&
     2.2257e-6_kr,3.07e4_kr,1.5571e-6_kr,6.07e4_kr,9.1595e-7_kr,&
     1.2e5_kr,5.7934e-7_kr,2.01e5_kr,4.3645e-7_kr,2.83e5_kr,&
     3.8309e-7_kr,3.56e5_kr,3.6926e-7_kr,3.77e5_kr,3.4027e-7_kr,&
     3.99e5_kr,2.7387e-7_kr,4.42e5_kr,1.0075e-7_kr,4.74e5_kr,&
     2.1754e-7_kr,5.02e5_kr,2.6333e-7_kr,5.4e5_kr,3.0501e-7_kr,&
     6.5e5_kr,2.9493e-7_kr,7.7e5_kr,2.5005e-7_kr,9.e5_kr,&
     2.1479e-7_kr,9.41e5_kr,1.7861e-7_kr,1.e6_kr,9.1595e-8_kr,&
     1.05e6_kr,1.1518e-7_kr/)
   real(kr),dimension(92),parameter::w2=(/&
     1.12e6_kr,1.3648e-7_kr,1.19e6_kr,1.5479e-7_kr,1.21e6_kr,&
     1.5022e-7_kr,1.31e6_kr,6.8696e-8_kr,1.4e6_kr,1.2182e-7_kr,&
     2.22e6_kr,5.9033e-8_kr,2.35e6_kr,9.1595e-8_kr,2.63e6_kr,&
     3.9981e-8_kr,3.e6_kr,3.1142e-8_kr,4.e6_kr,1.7073e-8_kr,&
     5.e6_kr,9.0679e-9_kr,6.e6_kr,4.7153e-9_kr,8.e6_kr,1.2276e-9_kr,&
     1.e7_kr,3.0953e-10_kr,1.257e7_kr,2.4619e-10_kr,1.26e7_kr,&
     3.4731e-10_kr,1.27e7_kr,1.0357e-9_kr,1.28e7_kr,2.8436e-9_kr,&
     1.29e7_kr,7.191e-9_kr,1.3e7_kr,1.6776e-8_kr,1.31e7_kr,&
     3.6122e-8_kr,1.32e7_kr,7.1864e-8_kr,1.33e7_kr,1.3222e-7_kr,&
     1.34e7_kr,2.2511e-7_kr,1.35e7_kr,3.5512e-7_kr,1.36e7_kr,&
     5.1946e-7_kr,1.37e7_kr,7.0478e-7_kr,1.38e7_kr,8.8825e-7_kr,&
     1.39e7_kr,1.0408e-6_kr,1.407e7_kr,1.154e-6_kr,1.42e7_kr,&
     1.087e-6_kr,1.43e7_kr,9.5757e-7_kr,1.44e7_kr,7.7804e-7_kr,&
     1.45e7_kr,6.0403e-7_kr,1.46e7_kr,4.3317e-7_kr,1.47e7_kr,&
     2.9041e-7_kr,1.48e7_kr,1.8213e-7_kr,1.49e7_kr,1.0699e-7_kr,&
     1.5e7_kr,5.8832e-8_kr,1.51e7_kr,3.0354e-8_kr,1.52e7_kr,&
     1.4687e-8_kr,1.53e7_kr,6.6688e-9_kr,1.54e7_kr,2.845e-9_kr,&
     1.55e7_kr,1.1406e-9_kr,1.5676e7_kr,1.978e-10_kr,2.e7_kr,&
     1.5477e-10_kr/)
   real(kr),dimension(66),parameter::w8=(/&
     0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,1.e0_kr,29.e0_kr,29.e0_kr,&
     5.e0_kr,.139000e-03_kr,.751516e-03_kr,.100000e-01_kr,&
     .497360e-01_kr,.200000e-01_kr,.754488e-01_kr,.400000e-01_kr,&
     .107756e+00_kr,.600000e-01_kr,.110520e+00_kr,.800000e-01_kr,&
     .101542e+00_kr,.100000e+00_kr,.884511e-01_kr,.614000e+02_kr,&
     .144057e-03_kr,.788930e+02_kr,.217504e-03_kr,.312030e+03_kr,&
     .127278e-02_kr,.179560e+04_kr,.236546e-02_kr,.804730e+04_kr,&
     .114311e-02_kr,.463090e+05_kr,.387734e-03_kr,.161630e+06_kr,&
     .125319e-03_kr,.639280e+06_kr,.207541e-04_kr,.286500e+07_kr,&
     .216111e-05_kr,.472370e+07_kr,.748998e-06_kr,.100000e+08_kr,&
     .573163e-07_kr,.127900e+08_kr,.940528e-08_kr,.129000e+08_kr,&
     .973648e-08_kr,.135500e+08_kr,.985038e-07_kr,.137500e+08_kr,&
     .176388e-06_kr,.139500e+08_kr,.239801e-06_kr,.140700e+08_kr,&
     .251963e-06_kr,.141900e+08_kr,.239298e-06_kr,.143900e+08_kr,&
     .176226e-06_kr,.145900e+08_kr,.992422e-07_kr,.155500e+08_kr,&
     .150737e-08_kr,.200000e+08_kr,.725000e-10_kr/)
   real(kr),dimension(102),parameter::w9=(/&
     0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,1e0_kr,47e0_kr,47e0_kr,5e0_kr,&
      1.39e-4_kr,3.019e6_kr,5.e-4_kr,1.07e7_kr,1.e-3_kr,2.098e7_kr,&
      5.e-3_kr,8.939e7_kr,1.e-2_kr,1.4638e8_kr,2.5e-2_kr,2.008e8_kr,&
      4.e-2_kr,1.7635e8_kr,5.e-2_kr,1.478e8_kr,1.e-1_kr,4.e7_kr,&
      1.4e-1_kr,1.13e7_kr,1.5e-1_kr,7.6e6_kr,4.14e-1_kr,2.79e6_kr,&
      1.13e0_kr,1.02e6_kr,3.06e0_kr,3.77e5_kr,8.32e0_kr,1.39e5_kr,&
      2.26e1_kr,5.11e4_kr,6.14e1_kr,1.88e4_kr,1.67e2_kr,6.91e3_kr,&
      4.54e2_kr,2.54e3_kr,1.235e3_kr,9.35e2_kr,3.35e3_kr,3.45e2_kr,&
      9.12e3_kr,1.266e2_kr,2.48e4_kr,4.65e1_kr,6.76e4_kr,1.71e1_kr,&
      1.84e5_kr,6.27e0_kr,3.03e5_kr,3.88e0_kr,5.e5_kr,3.6e0_kr,&
      8.23e5_kr,2.87e0_kr,1.353e6_kr,1.75e0_kr,1.738e6_kr,1.13e0_kr,&
      2.232e6_kr,0.73e0_kr,2.865e6_kr,0.4e0_kr,3.68e6_kr,2.05e-1_kr,&
      6.07e6_kr,3.9e-2_kr,7.79e6_kr,1.63e-2_kr,1.e7_kr,6.5e-3_kr,&
      1.2e7_kr,7.6e-3_kr,1.3e7_kr,1.23e-2_kr,1.35e7_kr,2.64e-2_kr,&
      1.4e7_kr,1.14e-1_kr,1.41e7_kr,1.14e-1_kr,1.42e7_kr,1.01e-1_kr,&
      1.43e7_kr,6.5e-2_kr,1.46e7_kr,1.49e-2_kr,1.5e7_kr,4.e-3_kr,&
      1.6e7_kr,1.54e-3_kr,1.7e7_kr,0.85e-3_kr/)
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::zero=0.e0_kr
   real(kr),parameter::onep5=1.5e0_kr

   !--read flux calculator input, if any.
   iwtt=iabs(iwt)
   nflmax=0
   sigpot=0
   if (iwt.lt.0) then
      jsigz=1
      alpha2=0
      sam=0
      beta=0
      alpha3=0
      gamma=0
      read(nsysi,*) fehi,sigpot,nflmax,ninwt,jsigz,&
        alpha2,sam,beta,alpha3,gamma
      ninwt=iabs(ninwt)
      call openz(-ninwt,1)
      write(nsyso,'(/&
        &'' compute flux...fehi, sigpot, nflmax ='',f9.1,f9.2,i8/&
        &4x,''ninwt, jsigz ='',i3,i4)')&
        fehi,sigpot,nflmax,-ninwt,jsigz
      if (alpha3.ne.zero.and.alpha2.eq.zero) alpha2=small
      if (alpha2.gt.zero) write(nsyso,'(&
        &7x,''alpha2, sam ='',1p,e12.4,0p,f8.3/&
        &7x,''beta, alpha3, gamma ='',f6.3,1p,e12.4,0p,f8.4)')&
        alpha2,sam,beta,alpha3,gamma
   endif

   !--compute flux from data on ninwt
   if (iwtt.eq.0) then
      write(nsyso,'(/'' compute flux......from data on ninwt'')')
      read(nsysi,*) ninwt
      ninwt=iabs(ninwt)
      call openz(-ninwt,0)
      write(nsyso,'(/'' ninwt......'',i4)') -ninwt
      return
   endif

   !--arbitary
   if (iwtt.eq.1) then
      write(nsyso,'(/'' weight function......read in'')')
      ntmp=10000
      allocate(tmp(ntmp))
      read(nsysi,*) (tmp(i),i=1,ntmp)
      nr=nint(tmp(5))
      np=nint(tmp(6))
      iw=6+2*nr+2*np
      if (iw.gt.ntmp) call error('genwtf',&
        'exceeded storage reading user weight function',' ')
      allocate(wght(iw))
      do i=1,iw
         wght(i)=tmp(i)
      enddo
      deallocate(tmp)

   !--constant
   else if (iwtt.eq.2) then
      write(nsyso,'(/'' weight function......constant for all l'')')

   !--1/e
   else if (iwtt.eq.3) then
      write(nsyso,'(/'' weight function......1/e for all l'')')

   !--1/e+fission+thermal
   else if (iwtt.eq.4) then
      read(nsysi,*) eb,tb,ec,tc
      if (eb.gt.50*tb) then
         ab=1
         ac=0
      else
         ab=1/(exp(-eb/tb)*eb**2)
         ac=1/(exp(-ec/tc)*ec**onep5)
      endif
      iw=6
      allocate(wght(iw))
      wght(1)=eb
      wght(2)=tb
      wght(3)=ab
      wght(4)=ec
      wght(5)=tc
      wght(6)=ac
      write(nsyso,'(/&
        &'' weight function......thermal + 1/e + fission''/&
        &'' thermal breakpoint and temperature  '',1p,2e12.4/&
        &'' fission breakpoint and temperature  '',2e12.4)')&
        eb,tb,ec,tc

   !--epri-cell light water reactor weight.
   else if (iwtt.eq.5) then
      write(nsyso,'(/'' weight function......epri-cell lwr'')')
      iw=184
      allocate(wght(iw))
      do i=1,92
         wght(i)=w1(i)
      enddo
      do i=1,92
         wght(i+92)=w2(i)
      enddo

   !--(thermal) -- (1/e) -- (fission + fusion)
   else if (iwtt.eq.6.or.iwtt.eq.7) then
      write(nsyso,'(/&
        &'' weight function......(thermal) -- (1/e) -- '',&
        &''(fission + fusion)'')')
      if (iwtt.gt.6) write(nsyso,'(22x,''temperature dependent'')')

   !--thermal--1/e--fast reactor--fission + fusion
   else if (iwtt.eq.8) then
      write(nsyso,'(/&
        &'' weight function...thermal--1/e--fast reactor--'',&
        &''fission + fusion'')')
      iw=66
      allocate(wght(iw))
      do i=1,66
         wght(i)=w8(i)
      enddo

   !--claw weight function
   else if (iwtt.eq.9.or.iwtt.eq.10) then
      write(nsyso,&
        &'(/'' weight function......claw weight function'')')
      if (iwtt.gt.9) then
         write(nsyso,'(22x,''temperature dependent'')')
      endif
      iw=102
      allocate(wght(iw))
      do i=1,102
         wght(i)=w9(i)
      enddo

   !--vitamin-e weight function
   else if (iwtt.eq.11.or.iwtt.eq.12) then
      write(nsyso,'(/'' weight function......vitamin-e'')')
      if (iwtt.gt.11) write(nsyso,'(22x,''temperature dependent'')')

   !--illegal iwt
   else
      call error('genwtf','illegal weight function.',' ')
   endif
   return
   end subroutine genwtf

   subroutine getwtf(e,enext,idis,lord,wtf)
   !-------------------------------------------------------------------
   ! Retrieve or compute required Legendre component of the
   ! weight function constructed or read in by genwtf.
   !-------------------------------------------------------------------
   use physics ! provides bk
   use endf    ! provides terpa
   use util    ! provides sigfig
   ! externals
   integer::idis,lord
   real(kr)::e,enext,wtf
   ! internals
   integer::iwtt,ir,ip,ipl
   real(kr)::step,tt,bb,cc,pow,test,ea,eb,enxt
   real(kr),parameter::con1=7.45824e+07_kr
   real(kr),parameter::con2=1.e0_kr
   real(kr),parameter::con3=1.44934e-09_kr
   real(kr),parameter::con4=3.90797e-02_kr
   real(kr),parameter::con5=2.64052e-05_kr
   real(kr),parameter::con6=6.76517e-02_kr
   real(kr),parameter::en1=.414e0_kr
   real(kr),parameter::en2=2.12e6_kr
   real(kr),parameter::en3=1.e7_kr
   real(kr),parameter::en4=1.252e7_kr
   real(kr),parameter::en5=1.568e7_kr
   real(kr),parameter::therm=.0253e0_kr
   real(kr),parameter::theta=1.415e6_kr
   real(kr),parameter::fusion=2.5e4_kr
   real(kr),parameter::ep=1.407e7_kr
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::s110=1.10e0_kr
   real(kr),parameter::s101=1.01e0_kr
   real(kr),parameter::s1002=1.002e0_kr
   real(kr),parameter::s1005=1.005e0_kr
   real(kr),parameter::s10001=1.0001e0_kr
   real(kr),parameter::tenth=0.1e0_kr
   real(kr),parameter::half=0.5e0_kr
   real(kr),parameter::veb=5.e5_kr
   real(kr),parameter::wt6a=.054e0_kr
   real(kr),parameter::wt6b=1.578551e-3_kr
   real(kr),parameter::wt6c=2.1e6_kr
   real(kr),parameter::wt6d=2.32472e-12_kr
   real(kr),parameter::wt6e=1.4e6_kr
   real(kr),parameter::wt6f=2.5e4_kr
   real(kr),parameter::wt6g=1.407e7_kr
   real(kr),parameter::wt6h=2.51697e-11_kr
   real(kr),parameter::wt6i=1.6e6_kr
   real(kr),parameter::wt6j=3.3e5_kr
   real(kr),parameter::wt10a=.15e0_kr
   real(kr),parameter::wt10b=300.e0_kr
   real(kr),parameter::wt10c=1.15e6_kr
   real(kr),parameter::exmin=-89.e0_kr
   real(kr),parameter::two=2
   real(kr),parameter::zero=0
   save ip,ir,ipl,step

   !--initialize
   iwtt=iabs(iwt)
   idis=0
   if (e.eq.0) then
      ip=2
      ir=1
      ipl=0
      enext=emax
      step=s110
      return
   endif

   !--branch to desired method
   iwtt=iabs(iwt)

   !--tabulated
   if (iwtt.eq.1.or.iwtt.eq.5.or.iwtt.eq.8.or.iwtt.eq.9) then
      call terpa(wtf,e,enext,idis,wght,ip,ir)
      if (wtf.ne.zero) then
         if (ip.ne.ipl) then
            step=s10001*(enext/e)**tenth
            ipl=ip
         endif
         enxt=step*e
         if (enxt.gt.s101*e) enxt=s101*e
         if (enext.gt.enxt) idis=0
         if (enext.gt.enxt) enext=enxt
      endif

   !--constant for all orders
   else if (iwtt.eq.2) then
      wtf=1
      enext=s101*e

   !--1/e for all orders
   else if (iwtt.eq.3) then
      wtf=1/e
      enext=s101*e

   !--thermal + 1/e + fission
   !--wght(1) to wght(6) are eb, tb, ab, ec, tc, ac
   else if (iwtt.eq.4) then
      if (e.le.wght(1)) then
         wtf=wght(3)*e*exp(-e/wght(2))
         enext=s101*e
         if (e.lt.wght(1).and.enext.gt.wght(1)) enext=wght(1)
      else if (e.le.wght(4)) then
         wtf=1/e
         enext=s101*e
         if (e.lt.wght(4).and.enext.gt.wght(4)) enext=wght(4)
      else
         wtf=wght(6)*sqrt(e)*exp(-e/wght(5))
         enext=s101*e
      endif

   !--(thermal) -- (1/e) -- (fission + fusion)
   !--with optional t dependence
   else if (iwtt.eq.6.or.iwtt.eq.7) then
      tt=wt6a
      if (iwtt.gt.6) tt=temp(jtemp)*bk
      bb=2*tt
      cc=1
      if (iwtt.gt.6) cc=wt6b*exp(two)/bb**2
      if (e.le.bb) then
         wtf=cc*e*exp(-e/tt)
         enext=s101*e
      else if (e.le.wt6c) then
         wtf=wt6b/e
         enext=s101*e
      else
         wtf=wt6d*sqrt(e)*exp(-e/wt6e)
         pow=-(sqrt(e/wt6f)-sqrt(wt6g/wt6f))**2/2
         if (pow.gt.exmin) wtf=wtf+wt6h*exp(pow)
         enext=s101*e
         test=wt6i
         if (abs(e-wt6g).le.test) enext=s1005*e
         test=wt6j
         if (abs(e-wt6g).le.test) enext=s1002*e
      endif

   !--temperature-dependent thermal part
   else if (iwtt.eq.10) then
      ea=bk*temp(jtemp)
      eb=wt10a*temp(jtemp)/wt10b
      if (e.lt.eb) then
         wtf=wt10c*(e/eb**2)*exp(-(e-eb)/ea)
         enext=s101*e
         if (enext.gt.eb) enext=eb
      else
         call terpa(wtf,e,enext,idis,wght,ip,ir)
         if (wtf.eq.zero) then
            enext=emax
         else
            if (ip.ne.ipl) then
               step=s10001*(enext/e)**tenth
               ipl=ip
            endif
            enxt=step*e
            if (enxt.gt.s101*e) enxt=s101*e
            if (enext.gt.enxt) idis=0
            if (enext.gt.enxt) enext=enxt
         endif
      endif

   !--vitamin-e weight function (ornl-5510)
   !--with optional t dependence
   else if (iwtt.eq.11.or.iwtt.eq.12) then
      enext=s101*e
      if (e.lt.en1) then
         tt=therm
         if (iwtt.gt.11) tt=temp(jtemp)*bk
         cc=con1
         if (iwtt.gt.11) cc=con2*exp(en1/tt)/en1**2
         wtf=cc*e*exp(-e/tt)
         if (enext.gt.en1) enext=en1
      else if (e.lt.en2) then
         wtf=con2/e
         if (enext.gt.en2) enext=en2
      else if (e.lt.en3) then
         wtf=con3*e**half*exp(-e/theta)
         if (enext.gt.en3) enext=en3
      else if (e.lt.en4) then
         wtf=con4/e
         if (enext.gt.en4) enext=en4
      else if (e.lt.en5) then
         wtf=con5*exp(-5*(e**half-ep**half)**2/fusion)
         if (abs(e-ep).le.veb) enext=s1002*e
         if (enext.gt.en5) enext=en5
      else
         wtf=con6/e
      endif
   endif

   !--return enext on an even grid
   enext=sigfig(enext,6,1)
   return
   end subroutine getwtf

   subroutine genflx(b)
   !-------------------------------------------------------------------
   ! Flux calculator.
   ! If iwt.gt.0, the Bondarenko model is used.
   ! If iwt.lt.0, the flux is calculated for an infinite mixture
   ! of a heavy absorber and a light moderator.
   ! The method used is an iterative solution of the integral
   ! slowing down equation assuming
   !   1) isotropic scattering in cm system
   !   2) scattering from background atoms (sigma zero) gives
   !      a smooth asymptotic spectrum equal to the specified
   !      weight function (see iwt in genwtf)
   ! The flux is computed between felo (lowest group bound) and
   ! fehi (must be in resolved range) or until nflmax points
   ! have been generated.  The flux at higher energies is extended
   ! using the bondarenko narrow resonance model.
   ! Legendre moments are computed using the B0 large system
   ! approximation--phi(l)=phi(0)/sigma**(l+1).
   !
   ! Special features are available (see ir) for heterogeneity
   ! adjustment or for two moderators (external and/or mixed with
   ! fuel).
   !   alpha(2) scattering alpha for second moderator
   !            (may be admixed or external)
   !   alpha(3) scattering alpha for third moderator
   !            (must be externa)
   !   beta     heterogeneity parameter
   !   sam      second (or admixed) moderator barns/atom of absorber
   !   gamma    fraction of second moderator in exerna monerator
   ! An option for writing and reading a resonance flux to include
   ! interference effects is also available.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use endf   ! provides endf routines and variables
   use util  ! provides timer,repoz,loada,error
   ! externals
   real(kr)::b(*)
   ! internals
   integer::idis,ie,li,iz,lim,k,je,lj,l,il,indx,mt,nscr
   integer::nl,nw,nb,net2,nin,nalph,nemax,nx,net,ne,net1,nett
   real(kr)::e,enext,en,ebot,f1,ej,ejp,f2,f3,f4
   real(kr)::elim,g3,g4,ei,enxt,fac,flag,t,wtf,xc,time
   real(kr)::tt,ss,s1,s2,sigzj,ttt
   real(kr)::alpha(3)
   real(kr),dimension(:),allocatable::fout
   real(kr),dimension(:),allocatable::denom,add,factor,tot
   real(kr),dimension(:),allocatable::scr1,scr2
   real(kr),parameter::up=1.001e0_kr
   real(kr),parameter::dn=0.999e0_kr
   real(kr),parameter::tenth=0.1e0_kr
   real(kr),parameter::sigzmx=1.e10_kr
   real(kr),parameter::etop=1.e10_kr
   real(kr),parameter::zero=0

   !--allocate temporary arrays
   if (allocated(denom)) then
      deallocate(denom)
      deallocate(add)
      deallocate(factor)
      deallocate(tot)
   endif
   allocate(denom(nsigz))
   allocate(add(nsigz))
   allocate(factor(nsigz))
   allocate(tot(nsigz))

   nl=lord+1
   if (nl.eq.1.and.nsigz.gt.1) nl=2
   nfv=1+nl*nsigz+1
   allocate(fout(nfv))
   !--bondarenko flux model.
   if (nflmax.eq.0) then
      nw=npage+50
      allocate(scr1(nw))
      scr1=0
      nscr=npend
      call findf(matb,3,1,nscr)
      call contio(nscr,0,0,scr1,nb,nw)
      e=0
      call gety1(e,enext,idis,t,nscr,scr1)
      if (enext.gt.up*egn(1)) call error('genflx',&
        'total not defined over energy range.',' ')
      if (iwt.ne.0) call getwtf(e,en,idis,0,wtf)
      if (iwt.eq.0) call getfwt(e,en,wtf,xc)
      ie=0
      net2=0

   !--compute flux from slowing down equation.
   else
      call timer(time)
      write(nsyso,'(68x,f9.1,''s'')') time
      nin=npend
      nscr=nfscr
      alpha(2)=alpha2
      alpha(3)=alpha3
      nalph=1
      if (alpha2.ne.zero) nalph=2
      if (alpha3.ne.zero) nalph=3
      nemax=nflmax
      nx=3+nsigz
      felo=egn(1)
      ebot=felo
      if (felo.lt.tenth) felo=tenth
      nw=npage+50
      allocate(scr1(nw))
      allocate(scr2(nw))
      scr1=0
      scr2=0
      e=0
      call getwtf(e,en,idis,0,wtf)

      !--read required cross sections.
      call repoz(nscr)
      call findf(matb,3,1,nin)
      nsh=1
      call tosend(nin,nscr,0,scr1)
      call repoz(nscr)
      call contio(nscr,0,0,scr1,nb,nw)
      awr=scr1(2)
      alpha(1)=((awr-1)/(awr+1))**2
      e=0
      call gety1(e,enext,idis,tt,nscr,scr1)
      call findf(matb,3,2,nin)
      call contio(nin,0,0,scr2,nb,nw)
      call gety2(e,en,idis,ss,nin,scr2)
      enext=felo
      net=0
      ne=0
      do while (e.lt.fehi.and.ne.lt.nemax)
         ne=ne+1
         e=enext
         if (e.gt.fehi) e=fehi
         call gety1(e,enext,idis,tt,nscr,scr1)
         call gety2(e,en,idis,ss,nin,scr2)
         call getwtf(e,en,idis,0,wtf)
         li=nx*(ne-1)
         b(li+1)=e
         b(li+2)=tt
         b(li+3)=ss
         do iz=1,nsigz
            b(iz+3+li)=(sigz(iz)-sam)*wtf*(1-beta)
         enddo
      enddo
      fehi=e

      !--assume nr flux at fehi and above
      !--calculate contributions below fehi for nr flux
      do iz=1,nsigz
         b(iz+3+li)=(sigz(iz)+sigpot)*wtf/(sigz(iz)+tt)
      enddo
      lim=ne-1
      do k=1,nalph
         do je=1,lim
            lj=nx*(je-1)
            e=b(lj+1)
            if (e.ge.alpha(k)*fehi) then
               f1=(1-alpha(k)*fehi/e)*wtf/(1-alpha(k))
               do iz=1,nsigz
                  if (k.eq.1) b(lj+3+iz)=b(lj+3+iz)+f1*sigpot
                  if (k.eq.2) b(lj+3+iz)=b(lj+3+iz)&
                    +f1*(sam+beta*gamma*(sigz(iz)-sam))
                  if (k.eq.3) b(lj+3+iz)=b(lj+3+iz)&
                    +f1*beta*(1-gamma)*(sigz(iz)-sam)
               enddo
            endif
         enddo
      enddo

      !--solve slowing down equation from high to low energies
      je=ne-1
      do while (je.gt.0)
         lj=nx*(je-1)
         ej=b(lj+1)
         ejp=b(lj+1+nx)
         f1=ejp/(ejp-ej)
         f2=ej/(ejp-ej)
         f3=f1*log(ejp/ej)-1
         f4=1-f2*log(ejp/ej)
         do iz=1,nsigz
            denom(iz)=sigz(iz)+b(lj+2)
         enddo
         do k=1,nalph
            elim=ej/alpha(k)
            if (elim.gt.ejp) elim=ejp
            do iz=1,nsigz
               if (k.eq.1) s1=b(lj+3)
               if (k.eq.1) s2=b(lj+3+nx)
               if (k.eq.2) s1=sam+beta*gamma*(sigz(iz)-sam)
               if (k.eq.2) s2=s1
               if (k.eq.3) s1=beta*(1-gamma)*(sigz(iz)-sam)
               if (k.eq.3) s2=s1
               g3=f1*log(elim/ej)-(elim-ej)/(ejp-ej)
               g4=(elim-ej)/(ejp-ej)-f2*log(elim/ej)
               b(lj+3+iz)=b(lj+3+iz)+g4*s2*b(lj+3+iz+nx)/(1-alpha(k))
               denom(iz)=denom(iz)-g3*s1/(1-alpha(k))
            enddo
         enddo
         do iz=1,nsigz
            b(iz+3+lj)=b(iz+3+lj)/denom(iz)
         enddo
         do k=1,nalph
            do iz=1,nsigz
               if (k.eq.1) s1=b(lj+3)
               if (k.eq.1) s2=b(lj+3+nx)
               if (k.eq.2) s1=sam+beta*gamma*(sigz(iz)-sam)
               if (k.eq.2) s2=s1
               if (k.eq.3) s1=beta*(1-gamma)*(sigz(iz)-sam)
               if (k.eq.3) s2=s1
               add(iz)=(f3*s1*b(lj+3+iz)&
                 +f4*s2*b(lj+3+iz+nx))/(1-alpha(k))
            enddo
            elim=2*ej
            ie=je
            do while (ie.gt.1.and.elim.gt.ej)
               ie=ie-1
               li=nx*(ie-1)
               ei=b(li+1)
               elim=ei/alpha(k)
               if (elim.gt.ej) then
                  if (elim.ge.ejp) then
                     do iz=1,nsigz
                        b(iz+3+li)=b(iz+3+li)+add(iz)
                     enddo
                  else
                     do iz=1,nsigz
                        if (k.eq.1) s1=b(lj+3)
                        if (k.eq.1) s2=b(lj+3+nx)
                        if (k.eq.2) s1=sam+beta*gamma*(sigz(iz)-sam)
                        if (k.eq.2) s2=s1
                        if (k.eq.3) s1=beta*(1-gamma)*(sigz(iz)-sam)
                        if (k.eq.3) s2=s1
                        g3=f1*log(elim/ej)-(elim-ej)/(ejp-ej)
                        g4=(elim-ej)/(ejp-ej)-f2*log(elim/ej)
                        b(li+3+iz)=b(li+3+iz)+(g3*s1*b(lj+3+iz)&
                          +g4*s2*b(lj+3+iz+nx))/(1-alpha(k))
                     enddo
                  endif
               endif
            enddo
         enddo
         je=je-1
      enddo

      !--flux is converged
      !--use specified weight function below elow
      e=0
      call getwtf(e,enxt,idis,0,wtf)
      e=b(1)
      call getwtf(e,enxt,idis,0,wtf)
      do iz=1,nsigz
         factor(iz)=b(iz+3)/wtf
      enddo
      e=0
      call getwtf(e,enxt,idis,0,wtf)
      net=0
      enxt=ebot
      do while (enxt.lt.b(1))
         net=net+1
         e=enxt
         call getwtf(e,enxt,idis,0,wtf)
         fout(1)=e
         do iz=1,nsigz
            l=(iz-1)*nl+1
            do il=1,nl
               fout(il+l)=factor(iz)*wtf
            enddo
         enddo
         call loada(net,fout,nfv,nflx,wtbuf,nbuf)
      enddo
      write(nsyso,'(&
        &'' flux calculator used weight function from 1e-5 to'',&
        &1p,e12.4,'' ev''/1x,i6,'' points'')') fout(1),net
      net1=net

      !--output converged flux
      sigzj=sigzmx
      indx=1
      if (ninwt.ne.0.and.iwt.ne.0) then
         if (jtemp.eq.1) call repoz(ninwt)
         write(ninwt) temp(jtemp)
         sigzj=sigz(jsigz)
         indx=2+nl*(jsigz-1)
         call repoz(nscr)
         call contio(nscr,0,0,scr1,nb,nw)
         e=0
         call gety1(e,enext,idis,tt,nscr,scr1)
      endif
      do ie=1,ne
         net=net+1
         li=nx*(ie-1)
         fout(1)=b(li+1)
         do iz=1,nsigz
            fac=(sigpot+sigz(iz))/(b(li+2)+sigz(iz))
            do il=1,nl
               l=1+il+nl*(iz-1)
               if (il.eq.1) fout(l)=b(li+3+iz)
               if (il.gt.1) fout(l)=fout(l-1)*fac
            enddo
         enddo
         call loada(net,fout,nfv,nflx,wtbuf,nbuf)
         if (ninwt.ne.0.and.iwt.ne.0) then
            e=fout(1)
            call gety1(e,enext,idis,tt,nscr,scr1)
            xc=1+tt/sigzj
            write(ninwt) e,fout(indx),xc
         endif
      enddo
      ie=net
      write(nsyso,'(&
        &'' computed flux from'',1p,e12.4,'' to '',1p,e12.4,'' ev''/&
        &1x,i6,'' points'')') felo,fout(1),net-net1
      net2=net
   endif

   !--complete flux using nr approximation
   nett=1
   do while (nett.gt.0)
      ie=ie+1
      e=enext
      if (iwt.ne.0) call getwtf(e,enext,idis,0,wtf)
      if (iwt.eq.0) call getfwt(e,enext,wtf,xc)
      call gety1(e,en,idis,tt,nscr,scr1)
      if (en.lt.enext) enext=en
      if (en.gt.dn*etop) enext=etop
      fout(1)=e
      ttt=tt
      do il=1,nl
         if (il.eq.1) mt=1
         if (il.gt.1) mt=261
         do iz=1,nsigz
            tot(iz)=ttt
         enddo
         xtot=0
         call getunr(mt,e,en,tot)
         if (en.ge.zero) then
            if (en.lt.enext) idis=0
            if (en.lt.enext) enext=en
         endif
         do iz=1,nsigz
            tt=tot(iz)
            fac=0
            if (iwt.ne.0) fac=(sigpot+sigz(iz))/(tt+sigz(iz))
            if (iwt.eq.0) fac=sigz(iz)*xc/(sigz(iz)*xc+tt)
            l=1+il+nl*(iz-1)
            if (il.eq.1) fout(l)=wtf*fac
            if (il.gt.1) fout(l)=fout(l-1)*fac
         enddo
      enddo
      fout(nfv)=xtot
      net=ie
      nett=net
      if (enext.gt.dn*etop) nett=-nett
      call loada(nett,fout,nfv,nflx,wtbuf,nbuf)
      if (ninwt.ne.0.and.iwt.ne.0) then
         write(ninwt) e,wtf,xc
      endif
   enddo
   nfp=net
   write(nsyso,'('' finished with narrow-resonance flux to'',&
     & 1p,e12.4,'' ev''/1x,i6,'' points'')') fout(1),net-net2

   deallocate(fout)
   if (allocated(denom)) then
      deallocate(denom)
      deallocate(add)
      deallocate(factor)
      deallocate(tot)
   endif

   if (allocated(scr1)) deallocate(scr1)
   if (allocated(scr2)) deallocate(scr2)
   if (ninwt.eq.0.or.iwt.eq.0) return
   flag=-1
   write(ninwt) flag,flag,flag
   return
   end subroutine genflx

   subroutine getfwt(e,enext,w,xc)
   !-------------------------------------------------------------------
   ! Read ninwt for the weight function and xc value
   ! corresponding to e.
   !-------------------------------------------------------------------
   use endf ! provides terp1
   use util ! provides repoz,error
   ! externals
   real(kr)::e,enext,w,xc
   ! internals
   real(kr)::el,en,wl,wn,xl,xn,flag,tempn
   character(60)::strng
   real(kr),parameter::eps=1.e-6_kr
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::zero=0
   real(kr)::templ=1.e10_kr
   save el,en,xl,xn,wl,wn

   !--initialize when e=0
   if (e.gt.zero) go to 200
   if (jtemp.eq.1.or.temp(jtemp).lt.templ) call repoz(ninwt)
   el=0
   wl=0
   xl=0
   en=0
   wn=0
   xn=0
   ! search tape for desired temperature
   if (jtemp.eq.1) go to 110
   read(ninwt) flag
   if (flag.lt.zero) go to 110
   if (abs(temp(jtemp)-flag).gt.eps) go to 130
   tempn=flag
   go to 150
  110 continue
   read(ninwt,end=145) tempn
   if (abs(tempn-temp(jtemp)).le.eps) go to 150
  130 continue
   read(ninwt,end=145) flag
   if (flag.ge.zero) go to 130
   go to 110
  145 continue
   write(strng,'(''temperature '',1p,e11.4)') temp(jtemp)
   call error('getfwt',strng,' ')
  150 continue
   templ=tempn
   enext=emax
   return

   !--search for e on ninwt
  200 continue
   if (abs(e-el).lt.el*small.and.en.lt.0) then
      w=wl
      xc=xl
      enext=en
      go to 250
   endif
   if (abs(e-en).lt.en*small) go to 240
   if (e.gt.el*(1-small).and.e.lt.en*(1+small)) go to 230
   if (e.lt.el*(1-small))&
     call error('getfwt','requested e is out of order',' ')
  210 continue
   el=en
   wl=wn
   xl=xn
   read(ninwt,end=225) en,wn,xn
   if (abs(e-en).lt.en*small) go to 240
   if (e.gt.el*(1-small).and.e.lt.en*(1+small)) go to 230
   go to 210
  225 continue
   call error('getfwt','e outside range of data',' ')
   ! interpolate for requested e
  230 continue
   call terp1(el,wl,en,wn,e,w,2)
   call terp1(el,xl,en,xn,e,xc,2)
   enext=en
   go to 250
  240 continue
   w=wn
   xc=xn
   el=en
   wl=wn
   xl=xn
   read(ninwt) en,wn,xn
   enext=en
   enext=emax
  250 continue
   return
   end subroutine getfwt

   subroutine init(ng,nl,nz)
   !-------------------------------------------------------------------
   ! Select parameters appropriate for this reaction type.
   ! Set nz=nsigz for mfd=6 if shielded neutron scattering
   ! matrices are desired.
   ! Set nz=nsigz for mfd=16 and 17 if shielded photon production
   ! matrices are desired.
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use util ! provides skiprz
   ! externals
   integer::ng,nl,nz
   ! internals
   real(kr),dimension(:),allocatable::tmp
   integer,parameter::naed=14
   real(kr),parameter::ebig=1.e8_kr

   !--select control parameters
   ng=2
   if (mfd.eq.6) ng=ngn+1
   if (mfd.eq.8) ng=ngn+1
   if (mtd.eq.251.or.mtd.eq.252.or.mtd.eq.253) ng=3
   if (mtd.eq.452.or.mtd.eq.455.or.mtd.eq.456) ng=3
   if (mfd.eq.5) ng=ngn+1
   if (mfd.eq.16.or.mfd.eq.17) ng=ngg+1
   if (mfd.eq.18) ng=ngg+1
   if (mfd.ge.21.and.mfd.le.26) ng=ngn+1
   if (mfd.ge.31.and.mfd.le.36) ng=ngn+1
   nl=1
   if (mtd.eq.1.and.nsigz.gt.1) nl=2
   if (mfd.eq.5.and.mtd.eq.455) nl=ndelg
   if (mfd.eq.6) nl=lord+1
   if (mfd.eq.8) nl=lord+1
   if (mfd.eq.16.or.mfd.eq.17) nl=lord+1
   if (mfd.eq.18) nl=lord+1
   if (mfd.ge.21.and.mfd.le.26) nl=lord+1
   if (mfd.ge.31.and.mfd.le.36) nl=lord+1
   if (mtd.eq.18) nl=1
   if (mtd.eq.19.or.mtd.eq.20.or.mtd.eq.21.or.mtd.eq.38) nl=1
   nz=1
   if (mfd.eq.3.and.mtd.eq.1) nz=nsigz
   if (mfd.eq.3.and.mtd.eq.2) nz=nsigz
   if (mfd.eq.6.and.mtd.eq.2) nz=nsigz
   if (mfd.eq.3.and.mtd.eq.18) nz=nsigz
   if (mfd.eq.6.and.mtd.eq.18) nz=nsigz
   if (mfd.eq.3.and.mtd.eq.19) nz=nsigz
   if (mfd.eq.3.and.mtd.eq.51) nz=nsigz
   if (mfd.eq.3.and.mtd.eq.102) nz=nsigz
   if (mfd.eq.3.and.mtd.eq.301) nz=nsigz
   if (mfd.eq.3.and.mtd.eq.302) nz=nsigz
   if (mfd.eq.3.and.mtd.eq.318) nz=nsigz
   if (mfd.eq.3.and.mtd.eq.402) nz=nsigz
   if (mfd.eq.3.and.mtd.eq.443) nz=nsigz
   if (mfd.eq.3.and.mtd.eq.444) nz=nsigz
   econst=0
   if (mfd.ne.3.and.mfd.ne.8.and.mfd.ne.18.and.mfd.lt.10000000) then
      if (mtd.eq.18.or.mtd.eq.19) econst=ebig
      if (mtd.eq.452.or.mtd.eq.455.or.mtd.eq.456) econst=ebig
      if (mtd.eq.102) econst=ebig
   endif

   !--copy file 6 to npend2 for use by retrieval routines.
   if (mfd.eq.6.and.mtd.ge.221.and.mtd.le.250) then
      allocate(tmp(npage+50))
      call findf(matb,6,mtd,npend)
      npend2=naed
      if (npend.lt.0) npend2=-npend2
      call repoz(npend2)
      nsc=1
      math=1
      call afend(0,npend2)
      call tosend(npend,0,npend2,tmp)
      call amend(0,npend2)
      call atend(0,npend2)
      call repoz(npend2)
      call skiprz(npend,-4)
      deallocate(tmp)
   endif
   return
   end subroutine init

   subroutine panel(elo,ehi,ans,ff,nl,nz,ng,iglo,nlg,ng2g)
   !-------------------------------------------------------------------
   ! Perform generalized group constant integrals for one panel.
   ! The upper boundry of the panel is chosen to be the smallest
   ! of ehi, the next cross section point, the next flux point,
   ! and the next feed function point.  Use Lobatto quadrature
   ! with order two larger than that used for the feed function.
   ! Note that three integrals are accumulated for ratio quantities
   ! (e.g. mt251, mt452)--flux, production rate, and reaction rate.
   !-------------------------------------------------------------------
   use util ! provides error,mess
   ! externals
   integer::nl,nz,ng,iglo,nlg,ng2g
   real(kr)::elo,ehi,ans(nlg,nz,ng2g),ff(nlg,ng2g)
   ! internals
   integer::ig1w,iglow,idiscf,ng1,ig1,nq,iz,il
   integer::nqp,iq,ig,igt
   real(kr)::elow,ehigh,enext,en,yld,aq,bq,eq,a,b,rr
   real(kr)::eqw,wq,t1
   character(60)::strng1,strng2
   real(kr),dimension(2),parameter::qp2=(/-1.e0_kr,1.e0_kr/)
   real(kr),dimension(2),parameter::qw2=(/1.e0_kr,1.e0_kr/)
   real(kr),dimension(6),parameter::qp6=(/&
     -1.e0_kr,-.76505532e0_kr,-.28523152e0_kr,.28523152e0_kr,&
     .76505532e0_kr,1.e0_kr/)
   real(kr),dimension(6),parameter::qw6=(/&
     .06666667e0_kr,.37847496e0_kr,.55485838e0_kr,.55485838e0_kr,&
     .37847496e0_kr,.06666667e0_kr/)
   real(kr),dimension(10),parameter::qp10=(/&
     -1.e0_kr,-.9195339082e0_kr,-.7387738651e0_kr,&
     -.4779249498e0_kr,-.1652789577e0_kr,.1652789577e0_kr,&
     .4779249498e0_kr,.7387738651e0_kr,.9195339082e0_kr,1.e0_kr/)
   real(kr),dimension(10),parameter::qw10=(/&
     .0222222222e0_kr,.1333059908e0_kr,.2248893420e0_kr,&
     .2920426836e0_kr,.3275397612e0_kr,.3275397612e0_kr,&
     .2920426836e0_kr,.2248893420e0_kr,.1333059908e0_kr,&
     .0222222222e0_kr/)
   real(kr),parameter::rndoff=1.000002e0_kr
   real(kr),parameter::delta=0.999995e0_kr
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::zero=0.e0_kr
   real(kr),parameter::smin=1.e-9_kr
   real(kr),parameter::ttest=0.1e0_kr
   real(kr),parameter::trange=1.e-6_kr
   real(kr)::elast=0.e0_kr
   integer::idisc=0
   integer::matl=0
   integer::mfl=0
   integer::mtl=0
   integer::allerr
   save ig1w,iglow,eqw,nq,nqp,ig1,ng1
   save enext

   !--retrieve factors in integrands at lower boundary
   if (matd.ne.matl.or.mfd.ne.mfl.or.mtd.ne.mtl.or.ipan.eq.0) then
      ipan=1
      matl=matd
      mfl=mfd
      mtl=mtd
      ig1w=0
      iglow=0
      eqw=0
      nq=0
      elast=0
      idisc=0
      !--reset these arrays whenever matd, mfd or mtd changes
      if (allocated(stmp)) then
         deallocate(stmp)
         deallocate(slst)
         deallocate(ftmp)
         deallocate(flst)
      endif
      allocate(stmp(nz,nl))
      allocate(slst(nz,nl))
      allocate(ftmp(nz,nl),stat=allerr)
      allocate(flst(nz,nl))
   endif
   elow=elo
   if (elo.gt.ehi*(1+small)) call error('panel','elo.gt.ehi.',' ')
   if (abs(elo-elast).ge.elast*small) then
      if (elo*rndoff.lt.ehi) elo=elo*rndoff
      elast=elo
      call getflx(elo,en,idiscf,flst,nl,nz)
      call getsig(elo,enext,idisc,slst,nl,nz)
      if (enext.ge.emax*(1-small)) go to 115
      if (abs(en-enext).lt.enext*small.and.idiscf.gt.idisc)&
        idisc=idiscf
      if (en.lt.enext*(1-small)) idisc=idiscf
      if (en.lt.enext*(1-small)) enext=en
      call getff(elo,en,idiscf,yld,ff,nl,ng1,ig1,nq,nlg,ng2g)
      if (abs(en-enext).lt.enext*small.and.idiscf.gt.idisc)&
        idisc=idiscf
      if (en.lt.enext*(1-small)) idisc=idiscf
      if (en.lt.enext*(1-small)) enext=en
      nq=nq+2
      if (nq.gt.10) nq=10
   endif

   !--retrieve cross section and flux at upper boundary
   if (enext.lt.delta*ehi) then
      ehi=enext
      ehigh=ehi
      if (idisc.gt.0.and.ehi*delta.gt.elo) ehigh=ehi*delta
   else
      ehigh=delta*ehi
   endif
   call getflx(ehigh,en,idiscf,ftmp,nl,nz)
   call getsig(ehigh,enext,idisc,stmp,nl,nz)
   if (enext.gt.emax*(1-small)) go to 110
   if (abs(en-enext).lt.enext*small.and.idiscf.gt.idisc) idisc=idiscf
   if (en.lt.enext*(1-small)) idisc=idiscf
   if (en.lt.enext*(1-small)) enext=en
   go to 125

   !--check for zero cross section over entire panel
  110 continue
   if (stmp(1,1).ne.zero) go to 125
  115 continue
   do iz=1,nz
      do il=1,nl
         ans(il,iz,1)=1
         ans(il,iz,2)=0
      enddo
   enddo
   ng=2
   iglo=1
   en=emax
   nqp=0
   go to 220
  125 continue

   !--compute group fluxes assuming flux is linear over panel.
   aq=(ehi+elow)/2
   bq=(ehi-elow)/2
   do iz=1,nz
      do il=1,nl
         ans(il,iz,1)=ans(il,iz,1)+(ftmp(iz,il)+flst(iz,il))*bq
      enddo
   enddo

   !--integrate over panel using lobatto quadrature.
   do iq=1,nq
      if (nq.eq.2) then
         eq=aq+bq*qp2(iq)
         wq=bq*qw2(iq)
      else if (nq.eq.6) then
         eq=aq+bq*qp6(iq)
         wq=bq*qw6(iq)
      else if (nq.eq.10) then
         eq=aq+bq*qp10(iq)
         wq=bq*qw10(iq)
      else
         call error('panel','bad nq in panel',' ')
      endif
      eq=sigfig(eq,9,0)
      t1=(eq-elow)/(ehi-elow)

      !--retrieve the feed function at the quadrature point.
      !--first point was last point of previous panel.
      if (iq.gt.1) then
         if (eq.gt.ehigh) eq=ehigh
         call getff(eq,en,idiscf,yld,ff,nl,ng1,ig1,nqp,nlg,ng2g)
      endif

      !--accumulate the ng*nl*nz integrals simultaneously.
      !--assuming that the reaction rate is linear across the panel.
      if (iglo.eq.0) iglo=ig1
      do iz=1,nz
         do il=1,nl
            a=stmp(iz,il)*ftmp(iz,il)
            b=slst(iz,il)*flst(iz,il)
            rr=b+(a-b)*t1
            rr=rr*wq
            do ig=1,ng1
               if (ff(il,ig).ne.zero) then
                  igt=ig1+ig-iglo+1
                  if (igt.le.1) then
                     if (ig1.ne.ig1w.or.iglo.ne.iglow&
                       .or.abs(eq-eqw).ge.trange) then
                        write(strng1,'(''thermal range problem at'',&
                          &1p,e12.4)') eq
                        write(strng2,&
                          '(''scattering to group '',i3,&
                          &''.  expected '',i3)') ig1,iglo
                        call mess('panel',strng1,strng2)
                        ig1w=ig1
                        iglow=iglo
                        eqw=eq
                     endif
                  endif
                  if (igt.gt.1) then
                     ans(il,iz,igt)=ans(il,iz,igt)+rr*ff(il,ig)
                     if (igt.gt.ng) ng=igt
                  endif
               endif
            enddo
         enddo
      enddo
   enddo

   !--remove partial thermal group.
   if (mtd.ge.221.and.mtd.le.250) then
      if (elo.ge.ttest) then
         if (abs(stmp(1,1)).le.smin) then
            do iz=1,nz
               do il=1,nl
                  ans(il,iz,2)=0
               enddo
            enddo
            ng=2
         endif
      endif
   endif

   !--save cross section and flux for next panel.
   !--determine next point and quadrature order from last ff.
  220 continue
   elast=ehigh
   do iz=1,nz
      do il=1,nl
         slst(iz,il)=stmp(iz,il)
         flst(iz,il)=ftmp(iz,il)
      enddo
   enddo
   if (abs(en-enext).lt.enext*small.and.idiscf.gt.idisc)&
     idisc=idiscf
   if (en.lt.enext*(1-small)) idisc=idiscf
   if (en.lt.enext*(1-small)) enext=en
   if (enext.le.ehi) enext=rndoff*ehi
   nq=nqp+2
   if (nq.gt.10) nq=10
   return
   end subroutine panel

   subroutine displa(ig,ans,nl,nz,ng2,ig2lo,igzero,nlg,ng2g)
   !-------------------------------------------------------------------
   ! Display generalized group constants generated by panel.
   ! Remove small values.  sum over final groups.
   ! Ratio records are transformed to flux, yield, cross section.
   ! Write an appropriate heading if called with ig.lt.0.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use util   ! provides sigfig
   ! external
   integer::ig,nl,nz,ng2,ig2lo,igzero,nlg,ng2g
   real(kr)::ans(nlg,nz,ng2g)
   ! internals
   integer::i,ig2,igj,irat,itr,lt,l,nd,ilo,ihi,igt,izero,il
   integer::j,iz,locd,lim,i2
   integer::nlnz10
   real(kr),dimension(:),allocatable::result
   character(10),dimension(:),allocatable::field
   real(kr),parameter::smin=1.e-9_kr
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::zero=0
   save irat,itr
   igzero=0

   nlnz10=max(nl,nz,10)
   if (allocated(result)) deallocate(result)
   allocate(result(nlnz10))
   if (allocated(field)) deallocate(field)
   allocate(field(nlnz10))

   !--write appropriate heading.
   if (ig.ge.0) go to 210
   write(nsyso,'('' '')')
   itr=0
   if (mfd.eq.6.or.mfd.eq.8&
     .or.mfd.eq.16.or.mfd.eq.17.or.mfd.eq.18&
     .or.(mfd.ge.21.and.mfd.le.26)&
     .or.(mfd.ge.31.and.mfd.le.36)) go to 130
   if (iprint.ne.1) go to 150
   if (mfd.eq.5) go to 140
   if (nz.gt.1) go to 120
   if (nl.eq.1) write(nsyso,'(&
     &'' enrgy  group constants at''/&
     &'' group  infinite dilution'')')
   if (nl.gt.1) write(nsyso,'(&
     &'' enrgy  lgend  group constants at''/&
     &'' group  order  infinite dilution'')')
   go to 150
  120 continue
   if (nl.eq.1) then
      do i=2,nz
         call a10(sigz(i),field(i))
      enddo
      write(nsyso,'(&
        &'' enrgy  group constants vs sigma zero''/&
        &'' group  infinity '',1p,9a11,(/,17x,9a11))') (field(i),i=2,nz)
   endif
   if (nl.gt.1) then
      do i=2,nz
         call a10(sigz(i),field(i))
      enddo
      write(nsyso,'(&
        &'' enrgy  lgend  group constants vs sigma zero''/&
        &'' group  order  infinity '',1p,9a11,(/,24x,9a11))')&
        & (field(i),i=2,nz)
   endif
   go to 150
  130 continue
   itr=1
   if (iprint.ne.1) go to 150
   if (nl.eq.1.and.nz.eq.1) go to 135
   lt=nl-1
   if (lt.gt.4) lt=4
   if (nz.eq.1) write(nsyso,'(&
     &'' initl  final  group constants vs legendre order''/&
     &'' group  group  0'',4(10x,i1))') (l,l=1,lt)
   if (nz.gt.1) then
      do i=2,nz
         call a10(sigz(i),field(i))
      enddo
      write(nsyso,'(&
        &'' initl  final  lgend  group constants vs sigma zero''/&
        &'' group  group  order  infinity '',1p,9a11)')&
        (field(i),i=2,nz)
   endif
   go to 150
  135 continue
   write(nsyso,'(&
     &'' initl  final  isotropic matrix vs final group''/&
     &'' group  group  +0         +1         +2         +3'')')
   go to 150
  140 continue
   nd=nl
   if (nd.eq.1) write(nsyso,'('' normalized fission spectrum'')')
   if (nd.gt.1) then
      do i=1,nd
         call a10(dntc(i),field(i))
      enddo
      write(nsyso,'(&
        &'' normalized delayed neutron spectra vs time constant''/&
        &'' group'',1p,9a11)') (field(i),i=1,nd)
   endif
  150 continue
   write(nsyso,'('' '')')
   irat=0
   if (mtd.eq.251.or.mtd.eq.252.or.mtd.eq.253) irat=1
   if (mtd.eq.452.or.mtd.eq.455.or.mtd.eq.456) irat=1
   return

   !--write out results in appropriate format.
  210 continue
   if (mfd.eq.5) go to 400
   if (ig.eq.0) go to 460
   if (ig2lo.eq.0) go to 470
   if (nz.gt.1) go to 310
   if (itr.eq.0) go to 310
   if (nl.eq.1) go to 250

   !--infinite dilution transfer matrices (iz=1).
   ilo=0
   ihi=ng2
   do ig2=2,ng2
      igt=ig2lo+ig2-2
      izero=0
      do il=1,nl
         result(il)=0
         if (ans(il,1,ig2).ne.zero) then
            if (ans(il,1,1).eq.zero) ans(il,1,1)=big
            ans(il,1,1)=sigfig(ans(il,1,1),7,0)
            result(il)=ans(il,1,ig2)/ans(il,1,1)
            result(il)=sigfig(result(il),7,0)
            if (abs(result(il)).ge.smin) izero=1
         endif
      enddo
      if (ilo.eq.0.and.(izero.gt.0.or.ig2.eq.ng2)) ilo=ig2
      if (ilo.ne.0) then
         do il=1,nl
            ans(il,1,ig2-ilo+2)=result(il)
         enddo
         if (izero.gt.0) ihi=ig2
         if (izero.eq.1) igzero=1
         if (iprint.eq.1.and.izero.gt.0) then
            do i=1,nl
               call a10(result(i),field(i))
            enddo
            write(nsyso,'(1x,i5,2x,i5,1x,1p,10a11)')&
              &ig,igt,(field(i),i=1,nl)
         endif
      endif
   enddo
   ig2lo=ig2lo+ilo-2
   ng2=ihi-ilo+2
   return

   !--isotropic infinite dilution transfer matrix (il=1,iz=1).
  250 continue
   ilo=0
   ihi=ng2
   j=2
   do ig2=2,ng2
      if (ilo.gt.0) j=j+1
      if (ans(1,1,1).eq.zero) ans(1,1,1)=big
      ans(1,1,1)=sigfig(ans(1,1,1),7,0)
      ans(1,1,j)=ans(1,1,ig2)/ans(1,1,1)
      ans(1,1,j)=sigfig(ans(1,1,j),7,0)
      if (abs(ans(1,1,j)).lt.smin/1000) ans(1,1,j)=0
      if (ans(1,1,j).ne.zero) then
         if (ilo.eq.0) ilo=ig2
         ihi=ig2
         igzero=1
      endif
   enddo
   if (igzero.eq.0) ilo=ng2
   ig2lo=ig2lo+ilo-2
   ng2=ihi-ilo+2
   if (igzero.eq.0) return
   if (iprint.ne.1) return
   do ig2=2,ng2,6
      igt=ig2lo+ig2-2
      igj=ig2+5
      if (igj.gt.ng2) igj=ng2
      do j=ig2,igj
         call a10(ans(1,1,j),field(j-ig2+1))
      enddo
      write(nsyso,'(1x,i5,2x,i5,1x,1p,10a11)')&
        ig,igt,(field(j),j=1,igj-ig2+1)
   enddo
   return

   !--self-shielded cross-sections and transfer matrices.
  310 continue
   ilo=0
   ihi=ng2
   do ig2=2,ng2
      igt=ig2lo+ig2-2
      izero=0
      do il=1,nl
         do iz=1,nz
            result(iz)=0
            if (ig2.eq.2) then
               ans(il,iz,1)=sigfig(ans(il,iz,1),7,0)
               if (irat.eq.1) ans(il,iz,3)=sigfig(ans(il,iz,3),7,0)
            endif
            if (ans(il,iz,ig2).ne.zero) then
               locd=1
               if (irat.eq.1.and.ig2.eq.2) locd=3
               if (ans(il,iz,locd).eq.zero) ans(il,iz,locd)=big
               result(iz)=ans(il,iz,ig2)/ans(il,iz,locd)
               result(iz)=sigfig(result(iz),7,0)
               izero=1
            endif
         enddo
         if (ilo.eq.0.and.(izero.gt.0.or.ig2.eq.ng2)) ilo=ig2
         if (ilo.eq.0.and.irat.eq.1) ilo=ig2
         if (ilo.ne.0) then
            do iz=1,nz
               ans(il,iz,ig2-ilo+2)=result(iz)
            enddo
            if (izero.gt.0) ihi=ig2
            if (izero.eq.1) igzero=1
            if (izero.ne.0.and.iprint.eq.1) then
               if (irat.ne.1.or.ig2.ne.3) then
                  l=il-1
                  if (itr.eq.1) then
                     do i=1,nz
                        call a10(result(i),field(i))
                     enddo
                     write(nsyso,'(1x,i5,2x,i5,6x,i1,2x,1p,10a11,&
                       &(/,22x,10a11))')ig,igt,l,(field(i),i=1,nz)
                  else if (nl.gt.1) then
                     do i=1,nz
                        call a10(result(i),field(i))
                     enddo
                     write(nsyso,'(1x,i5,6x,i1,2x,1p,10a11,&
                        &(/,15x,10a11))')ig,l,(field(i),i=1,nz)
                     if (mfd.eq.3.and.mtd.eq.1) then
                        do i=1,nz
                           call a10(ans(il,i,1),field(i))
                        enddo
                        write(nsyso,'(5x,''flux'',1x,i1,2x,1p,10a11,&
                          &(/,13x,10a11))')l,(field(iz),iz=1,nz)
                     endif
                  else
                     do i=1,nz
                        call a10(ans(il,i,2),field(i))
                     enddo
                     write(nsyso,'(1x,i5,2x,1p,10a11,&
                       &(/,8x,10a11))')ig,(field(i),i=1,nz)
                     if (mfd.eq.3.and.mtd.eq.1) then
                        do i=1,nz
                           call a10(ans(1,i,1),field(i))
                        enddo
                        write(nsyso,'(2x,''flx'',1x,1p,10a11,&
                          &(/,6x,10a11))')(field(iz),iz=1,nz)
                     endif
                  endif
               endif
            endif
         endif
      enddo
   enddo
   ig2lo=ig2lo+ilo-2
   ng2=ihi-ilo+2
   return

   !--fission spectrum (mfd=5/mtd=18)
  400 continue
   if (mtd.ne.455) then
      ! prompt
      do j=1,ng2
         ans(1,1,j)=sigfig(ans(1,1,j),7,0)
      enddo
      igzero=1
      if (iprint.ne.1) return
      do j=1,ng2,10
         lim=j+9
         if (lim.gt.ng2) lim=ng2
         do i=j,lim
            call a10(ans(1,1,i),field(i-j+1))
         enddo
         write(nsyso,'(1x,i5,2x,1p,10a11)')&
           j,(field(i),i=1,lim-j+1)
      enddo
   else
      !--delayed
      igzero=1
      do j=1,ng2
         do il=1,nl
            ans(il,1,j)=sigfig(ans(il,1,j),7,0)
         enddo
      enddo
      if (iprint.ne.1) return
      do j=2,ng2
         i2=j-1
         do il=1,nl
            call a10(ans(il,1,j),field(il))
         enddo
         write(nsyso,'(1x,i5,2x,1p,10a11)')&
           i2,(field(il),il=1,nl)
      enddo
   endif
   return

   !--constant spectra
  460 continue
   do i=1,ng2
      ans(1,1,i)=sigfig(ans(1,1,i),7,0)
   enddo
   igzero=1
   if (iprint.ne.1) return
   do i=1,ng2,6
      lim=i+5
      if (lim.gt.ng2) lim=ng2
      i2=ig2lo-1+i
      do j=i,lim
         call a10(ans(1,1,j),field(j-i+1))
      enddo
      write(nsyso,'('' spec'',3x,i5,2x,1p,6a11)')&
        i2,(field(j),j=1,lim-i+1)
   enddo
   write(nsyso,'('' '')')
   return

   !--constant productions
  470 continue
   igzero=1
   do i=1,ng2
      ans(1,1,i)=sigfig(ans(1,1,i),7,0)
   enddo
   if (iprint.ne.1) return
   do i=1,ng2,6
      lim=i+5
      if (lim.gt.ng2) lim=ng2
      do j=i,lim
         call a10(ans(1,1,j),field(j-i+1))
      enddo
      write(nsyso,'(1x,i5,''   prod'',2x,1p,6a11)')&
        i,(field(j),j=1,lim-i+1)
   enddo
   write(nsyso,'('' '')')

   deallocate(result)

   return
   end subroutine displa

   subroutine getflx(e,enext,idis,flux,nl,nz)
   !-------------------------------------------------------------------
   ! Retrieve weighting flux.  For nz=1 (infinite dilution) use
   ! getwtf.  For self-shielded cases, use flux components from
   ! tape prepared by genflx.
   !-------------------------------------------------------------------
   use util ! provides finda,scana
   use endf ! provides terp1
   ! externals
   integer::idis,nl,nz
   real(kr)::e,enext,flux(nz,nl)
   ! internals
   integer::ip,i,llord,iz,il,l
   real(kr)::wtf
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::zero=0
   save ip

   !--initialize.
   if (e.eq.zero) then
      if (nsigz.gt.1) then
         if (allocated(falo)) then
             deallocate(falo)
             deallocate(fahi)
         endif
         allocate(falo(nfv))
         allocate(fahi(nfv))
         call finda(1,falo,nfv,nflx,wtbuf,nbuf)
         call finda(2,fahi,nfv,nflx,wtbuf,nbuf)
         enext=falo(1)
         ip=0
      else
         call getwtf(e,enext,idis,0,wtf)
      endif
      return
   endif

   !--normal entry
   if (nsigz.gt.1) then

      !--branch for nz>1 (self shielded)
      if (ip.eq.0) then
         call scana(e,ip,nfp,nfv,nflx,wtbuf,nbuf)
         call finda(ip,falo,nfv,nflx,wtbuf,nbuf)
         call finda(ip+1,fahi,nfv,nflx,wtbuf,nbuf)
      endif

      !--check energy range
      !--read next block of data if necessary
      do while (e.ge.fahi(1).and.ip.lt.nfp)
         do i=1,nfv
            falo(i)=fahi(i)
         enddo
         ip=ip+1
         call finda(ip,fahi,nfv,nflx,wtbuf,nbuf)
      enddo

      !--in the right energy range
      !--interpolate for desired flux components.
      llord=lord
      if (lord.eq.0.and.nsigz.gt.1) llord=1
      do iz=1,nz
         do il=1,nl
            l=1+il+(llord+1)*(iz-1)
            call terp1(falo(1),falo(l),fahi(1),fahi(l),e,flux(iz,il),2)
         enddo
      enddo
      call terp1(falo(1),falo(nfv),fahi(1),fahi(nfv),e,xtot,2)
      enext=fahi(1)
      idis=0
      if (ip.eq.nfp) enext=emax
   else

      !--branch for nz=1 (infinite dilution only)
      call getwtf(e,enext,idis,0,wtf)
      do il=1,nl
         flux(1,il)=wtf
      enddo
   endif
   return
   end subroutine getflx

   subroutine getyld(e,enext,idis,yld,mat,mf,mt,lfs,itape)
   !-------------------------------------------------------------------
   ! Retrieve or compute yield for requested mat,mf,mt on itape.
   ! If mt=455, retrieve delayed neutron time constants.
   ! Initialize if e=0.
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use util ! provides error
   ! externals
   integer::idis,mat,mf,mt,lfs,itape,i,nr,nc
   real(kr)::e,enext,yld,term
   ! internals
   integer::mft,nb,nw,nk,ik,lnu,iza,lnd,izn,lfn,ip,ir,na,loc
   integer::ntmp
   real(kr),dimension(:),allocatable::tmp
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::zero=0
   save lnu,ip,ir
   character(60)::strng

   !--initialize
   if (e.gt.zero) go to 200
   ntmp=10000
   allocate(tmp(ntmp))
   mft=mf
   if (mft.ge.40000000) mft=10
   if (mft.ge.30000000) mft=9
   if (mft.ge.20000000) mft=6
   call findf(mat,mft,mt,itape)
   call contio(itape,0,0,tmp,nb,nw)
   nk=0
   if (mf.gt.99) nk=nint(tmp(5))
   ik=0
   lnu=0
   if (mt.ge.452.and.mt.le.456) lnu=nint(tmp(4))
   iza=0
   if (mf.ge.10000000) iza=mod(mf/10,1000000)
   loc=1
   if (mt.eq.455) go to 110
   lnd=0
   if (lnu.ne.1) go to 130
   call listio(itape,0,0,tmp(loc),nb,nw)
   na=nw
   enext=emax
   go to 190
  110 continue
   call listio(itape,0,0,tmp(loc),nb,nw)
   lnd=nint(tmp(loc+4))
   if (lnd.gt.8) call error('getyld','illegal lnd.',' ')
   do i=1,lnd
      dntc(i)=tmp(i+5+loc)
   enddo
  130 continue
   ik=ik+1
   loc=1
   call tab1io(itape,0,0,tmp(loc),nb,nw)
   loc=loc+nw
   do while (nb.ne.0)
      call moreio(itape,0,0,tmp(loc),nb,nw)
      loc=loc+nw
   enddo
   if (nk.eq.0) go to 180
   izn=0
   if (mft.eq.9.or.mft.eq.10) izn=nint(tmp(3))
   if (mft.eq.6) izn=nint(tmp(1))
   if (mft.eq.9.or.mft.eq.10) lfn=nint(tmp(4))
   if (mft.eq.6) lfn=nint(tmp(3))
   if (iza.eq.izn.and.lfs.eq.lfn) go to 180
   if (mft.gt.6.and.izn.eq.0.and.lfs.eq.lfn) go to 180
   if (mft.eq.9.or.mft.eq.10) go to 170
   law=nint(tmp(4))
   call skip6(itape,0,0,tmp(loc),law)
  170 continue
   if (ik.lt.nk) go to 130
   write(strng,'(''unable to find nuclide for iza='',i7,'' lfs='',i3)') iza,lfs
   call error('getyld',strng,' ')
  180 continue
   na=loc-1
   ip=2
   ir=1
   nr=nint(tmp(5))
   enext=tmp(7+2*nr)
  190 continue
   idis=0
   allocate(yield(na))
   do i=1,na
      yield(i)=tmp(i)
   enddo
   deallocate(tmp)
   return

   !--tabulated data
  200 continue
   loc=1
   if (lnu.eq.1) go to 300
   call terpa(yld,e,enext,idis,yield(loc),ip,ir)
   return

   !--polynomial data
  300 continue
   yld=yield(loc+6)
   term=1
   nc=nint(yield(loc+4))
   do i=2,nc
      term=term*e
      yld=yld+yield(loc+5+i)*term
   enddo
   enext=emax
   idis=0
   return
   end subroutine getyld

   subroutine getsig(e,enext,idis,sig,nl,nz)
   !-------------------------------------------------------------------
   ! Retrieve the reaction cross-section defined by mfd and mtd.
   ! Remove discontinuities by moving second point up by eps.
   ! Use unresolved cross-section if needed.  Initialize if e=0.
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use util ! provides error
   ! externals
   integer::idis,nl,nz
   real(kr)::e,enext,sig(nz,nl)
   ! internals
   integer::nsig,mf,mt,nb,nw,i,il,iz,mtt,nskip,nfs,jfs
   real(kr)::en,s
   character(60)::strng
   real(kr),parameter::emin=1.e-5_kr
   real(kr),parameter::step=1.01e0_kr
   real(kr),parameter::vc=1.919e8_kr
   real(kr),parameter::ezero=1.e7_kr
   real(kr),parameter::step2=1.10e0_kr
   real(kr),parameter::zero=0
   save mt,nsig

   !--initialize
   if (e.eq.zero) then
      nsig=npend
      if (mfd.eq.5) then
         ! spectra
         enext=emin
         thresh=0
         lrflag=0
         allocate(sigma(1))
         return
      endif
      mf=3
      if (mfd.eq.13.or.mfd.eq.17) mf=13
      if (mfd.ge.10000000) mf=3
      if (mfd.ge.20000000) mf=3
      if (mfd.ge.30000000) mf=3
      if (mfd.ge.40000000) mf=10
      mt=0
      if (mtd.eq.259.or.mtd.eq.258.or.mtd.eq.257) then
         ! non cross section quantities
         enext=emin
         thresh=0
         lrflag=0
         allocate(sigma(1))
         return
      endif
      if (mtd.le.150) mt=mtd
      if (mtd.ge.152.and.mtd.le.200) mt=mtd
      if (mtd.ge.201.and.mtd.le.250) mt=mtd
      if (mtd.gt.299.and.mtd.lt.451) mt=mtd
      if (iverf.ge.6.and.mtd.eq.500) mt=mtd
      if (iverf.lt.6.and.mtd.ge.700.and.mtd.le.799) mt=mtd
      if (iverf.ge.6.and.mtd.ge.600.and.mtd.le.849) mt=mtd
      if (iverf.ge.6.and.mtd.ge.850.and.mtd.le.874) mt=mtd
      if (iverf.ge.6.and.mtd.ge.875.and.mtd.le.891) mt=mtd
      if (iverf.ge.6.and.mtd.eq.599) mt=mtd
      if (mtd.eq.251.or.mtd.eq.252.or.mtd.eq.253) mt=2
      if (mtd.eq.261) mt=261
      if (mtd.eq.452.or.mtd.eq.455.or.mtd.eq.456) mt=18
      if (mt.eq.0) then
         write(strng,'(i4,'' invalid in endf'')') mtd
         call error('getsig','illegal mt.',strng)
      endif
      call findf(matd,mf,mt,nsig)
      nw=npage+50
      if (allocated(sigma)) deallocate(sigma)
      allocate(sigma(nw))
      call contio(nsig,0,0,sigma,nb,nw)
      awr=sigma(2)
      if (awrp.ne.zero) awr=awr/awrp
      if (mf.eq.10) then
         nfs=n1h                                    !# of tab1's to check
         jfs=-1
         do i=1,nfs
            call tab1io(nsig,0,0,sigma,nb,nw)
            if (l1h.eq.izar .and. l2h.eq.lfs) jfs=i  !id the one we want ...
            do while (nb.ne.0)
               call moreio(nsig,0,0,sigma,nb,nw)
            enddo
         enddo
         if (jfs.lt.0) then
            write(strng,'("can''t find mf,mt,izar,lfs = ",3i9,i5)')&

                                       mf,mt,izar,lfs
            call error('getsig',strng,' ')
         endif
         nskip=jfs-1
         call skiprz(nsig,-1)
         call findf(matd,mf,mt,nsig)
         call contio(nsig,0,0,sigma,nb,nw)
         if (nskip.gt.0) then
            do i=1,nskip
               call tab1io(nsig,0,0,sigma,nb,nw)
               do while (nb.ne.0)
                  call moreio(nsig,0,0,sigma,nb,nw)
               enddo
            enddo
         endif
      endif
      call gety1(e,enext,idis,s,nsig,sigma)
      q=c2h
      lrflag=l2h
      if (mf.ne.3) lrflag=0
      thresh=(awr+1)*(-q)/awr
      alpha=(awr-1)**2/(awr+1)**2
      return
   endif

   !--retrieve point cross sections.
   !--check for unresolved cross sections.
   if (mtd.eq.259) then
      ! average reciprocal velocity
      idis=0
      enext=step*e
      sig(1,1)=1/sqrt(vc*e)
   else if (mtd.eq.258) then
      ! average lethargy
      idis=0
      enext=step*e
      sig(1,1)=log(ezero/e)
   else if (mtd.eq.257) then
      ! average energy
      idis=0
      enext=step*e
      sig(1,1)=e
   else if (mtd.eq.2.and.izap.gt.1) then
      ! cp elastic
      idis=0
      enext=step2*e
      do il=1,nl
         do iz=1,nz
            sig(iz,il)=1
         enddo
      enddo
   else
      ! other reactions
      call gety1(e,enext,idis,s,nsig,sigma)
      do il=1,nl
         do iz=1,nz
            sig(iz,il)=s
         enddo
         mtt=mt
         if (il.gt.1.and.mt.eq.1) mtt=261
         en=enext
         if (nz.gt.1) then
            call getunr(mtt,e,en,sig(1,il))
         endif
      enddo
      if (en.lt.zero) return
      if (en.lt.enext) idis=0
      if (en.lt.enext) enext=en
   endif
   end subroutine getsig

   subroutine stounr(mat,temp,ntape)
   !-------------------------------------------------------------------
   ! Read self-shielded unresolved resonance cross-sections
   ! from ntape and store them for retrieval with getunr.
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use util ! provides error,mess
   ! externals
   integer::mat,ntape
   real(kr)::temp
   ! internals
   integer::mfd,mtd,nb,nw,i,j,nsig0,l,ntmp
   real(kr)::temz
   character(60)::strng
   real(kr),dimension(:),allocatable::tmp

   !--initialize
   nunr=0
   if (ntape.eq.0) return

   !--locate mat and temp
   if (lrp.ne.-1) go to 100
   allocate(unr(1))
   return
  100 continue
   ntmp=10000
   allocate(tmp(ntmp))
   mfd=2
   mtd=152
   call findf(mat,2,0,ntape)
   call contio(ntape,0,0,tmp,nb,nw)
   if (mfh.gt.mfd) go to 105
   if (mth.eq.mtd) go to 110
   call tosend(ntape,0,0,tmp)
   call contio(ntape,0,0,tmp,nb,nw)
   if (mfh.gt.mfd) go to 105
   if (mth.ne.mtd) go to 105
   go to 110
  105 continue
   deallocate(tmp)
   allocate(unr(1))
   return

   !--read in desired data
  110 continue
   l=1
   intunr=n2h
   lssf=l1h
   call listio(ntape,0,0,tmp(l),nb,nw)
   temz=tmp(l)
   if (abs(temz-temp).le.1.) go to 120
   if (temz.gt.temp) then
      write(strng,'(''cannot find temp='',1p,2e12.3)') temp,temz
      call error('stounr',strng,' ')
   endif
   call tomend(ntape,0,0,tmp(l))
   go to 100
  120 continue
   nsig0=nint(tmp(4))
   nunr=nint(tmp(6))
   l=l+nw
   do while (nb.ne.0)
      if (l.gt.ntmp)&
        call error('stounr','storage exceeded.',' ')
      call moreio(ntape,0,0,tmp(l),nb,nw)
      l=l+nw
   enddo
   if (nsig0.gt.ntmp)&
     call error('stounr','storage exceeded.',' ')
   allocate(unr(l))
   do i=1,l
      unr(i)=tmp(i)
   enddo
   deallocate(tmp)

   !--check for sigmazero conflicts
   if (nsigz.eq.1) return
   if (nsig0.gt.1) go to 160
   call mess('stounr','no unresolved sigma zero data',&
     'will use infinite dilution only')
   go to 190
  160 continue
   do 170 i=1,nsigz
   do 180 j=1,nsig0
   if (abs(sigz(i)-unr(6+j)).lt.sigz(i)/100) go to 170
  180 continue
   call mess('stounr','sigma zero grids do not match',&
     &  'will interpolate for desired values')
   go to 190
  170 continue
  190 continue
   return
   end subroutine stounr

   subroutine getunr(mt,e,enext,sig)
   !-------------------------------------------------------------------
   ! Retrieve unresolved resonance cross-sections stored by getunr
   ! and correct sig for self-shielding effects.  If unresolved
   ! has resolved resonance overlap, add it to the background.
   ! If e is outside the unresolved energy range or if this mat has
   ! no unresolved parameters (nunr=0), return enext=-1.
   !-------------------------------------------------------------------
   use endf ! provides terp1
   use util ! provides mess
   ! externals
   integer::mt
   real(kr)::e,enext,sig(*)
   ! internals
   integer::ibase,ix,lim,iflag,i,is,ib,nx,nsig0,ncyc
   integer::locl,locs,loc,loc1,iovl
   real(kr)::enunr,s,sl,sn
   real(kr)::e1,elo,ehi,el,en,sigb,sinf
   real(kr),parameter::eps=1.e-4_kr
   real(kr),parameter::del=1.e-10_kr
   real(kr),parameter::zero=0

   !--check nunr and energy range
   enext=-1
   iovl=0
   if (nunr.eq.0) return
   nx=nint(unr(3))
   nsig0=nint(unr(4))
   if (nsig0.le.1) return
   ncyc=1+nx*nsig0
   ibase=7+nsig0
   e1=abs(unr(ibase))
   if (e.lt.e1-del*e1) return
   enunr=abs(unr(ibase+(nunr-1)*ncyc))
   if (e.gt.enunr+del*enunr) return

   !--determine reaction type
   ix=0
   if (mt.eq.1) ix=1
   if (mt.eq.2) ix=2
   if (mt.eq.18) ix=3
   if (mt.eq.19) ix=3
   if (mt.eq.102) ix=4
   if (mt.eq.261) ix=5
   if (ix.eq.0) return
   locl=ibase+(nunr-1)*ncyc+nsig0*(ix-1)
   locs=7

   !--find interpolation range
   lim=nunr-1
   iflag=0
   i=0
   do while (iflag.eq.0.and.i.lt.lim)
      i=i+1
      loc=ibase+(i-1)*ncyc
      loc1=loc+ncyc
      elo=abs(unr(loc))
      elo=elo-del*elo
      ehi=abs(unr(loc1))
      if (e.gt.elo.and.e.lt.ehi) iflag=1
   enddo

   !--special exit for last point
   if (iflag.eq.0) then
      sinf=unr(locl+1)
      if (unr(locl).lt.zero) iovl=1
      if (ix.eq.1.and.iovl.eq.1) xtot=sig(1)-sinf
      do is=1,nsigz
         sigb=sigz(is)
         if (iovl.eq.1) sigb=sigb+xtot
         call terpu(s,sigb,unr(locl+1),unr(locs),nsig0)
         if (lssf.eq.0) sig(is)=sig(is)+s-sinf
         if (lssf.gt.0.and.sinf.ne.zero) sig(is)=sig(is)*s/sinf
      enddo
      enext=enunr*(1+eps)

   !--interpolate for cross-sections requested
   else
      el=abs(unr(loc))
      en=abs(unr(loc1))
      if (unr(loc).lt.zero.and.unr(loc1).lt.zero) iovl=1
      ib=nsig0*(ix-1)+1
      call terp1(el,unr(loc+ib),en,unr(loc1+ib),e,sinf,intunr)
      if (ix.eq.1.and.iovl.eq.1) xtot=sig(1)-sinf
      do is=1,nsigz
         sigb=sigz(is)
         if (iovl.eq.1) sigb=sigb+xtot
         call terpu(sl,sigb,unr(loc+ib),unr(locs),nsig0)
         call terpu(sn,sigb,unr(loc1+ib),unr(locs),nsig0)
         if (sl.lt.zero.or.sn.lt.zero) call mess('getunr',&
           ' Warning, negative URR cross sections found, check unresr',' ')
         call terp1(el,sl,en,sn,e,s,intunr)
         if (lssf.eq.0) sig(is)=sig(is)+s-sinf
         if (lssf.gt.0.and.sinf.ne.zero) sig(is)=sig(is)*s/sinf
      enddo
      enext=en
   endif
   return
   end subroutine getunr

   subroutine terpu(sig,sigz,xs,sigs,ns)
   !-------------------------------------------------------------------
   ! Interpolate for the cross section at a desired sigma zero value
   ! for unresolved data.
   !-------------------------------------------------------------------
   ! externals
   integer::ns
   real(kr)::sig,sigz,xs(ns),sigs(ns)
   ! internals
   integer::is,js
   real(kr)::terp,sm,st
   real(kr),parameter::rlo=0.99e0_kr
   real(kr),parameter::rhi=10.1e0_kr

   sig=0
   do 110 is=1,ns
   terp=1
   sm=sigs(is)
   if (sigz.lt.sigs(ns).and.sm.gt.rhi*sigs(ns)) go to 110
   if (sigz.lt.sigs(2).and.sm.gt.rhi*sigs(2)) go to 110
   if (sigz.ge.sigs(2).and.sm.lt.rlo*sigs(2)) go to 110
   do 120 js=1,ns
   st=sigs(js)
   if (abs(sm-st).lt.sm/1000) go to 120
   if (sigz.ge.sigs(ns)) go to 130
   if (st.gt.rhi*sigs(ns)) go to 120
   terp=terp*(sigz**2-st**2)/(sm**2-st**2)
   go to 120
  130 continue
   if (sigz.ge.sigs(2)) go to 140
   if (st.gt.rhi*sigs(2)) go to 120
   terp=terp*log(sigz/st)/log(sm/st)
   go to 120
  140 continue
   if (st.lt.rlo*sigs(2)) go to 120
   terp=terp*((st/sigz)-1)/((st/sm)-1)
  120 continue
   sig=sig+terp*xs(is)
  110 continue
   return
   end subroutine terpu

   subroutine getff(e,enext,idisc,yld,ff,nl,ng,iglo,nq,nlg,ng2g)
   !-------------------------------------------------------------------
   ! Compute feed function or yield for desired reaction type.
   !-------------------------------------------------------------------
   use util ! provides error
   use endf ! provides iverf
   ! externals
   integer::idisc,nl,ng,iglo,nq,nlg,ng2g
   real(kr)::e,enext,yld,ff(nlg,ng2g)
   ! internals
   integer::il,ifirst,nyl,igmin,idis,nle,lcd,nk
   integer::ngn1,ig,ighi,ik,jg,mft,iglo1
   real(kr)::en,econ,epl,eph,tempo,temp
   real(kr)::y(1),eccn(1)
   real(kr)::s(1,1)
   character(40)::com
   real(kr)::fle(21)
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::smin=1.e-9_kr
   real(kr),parameter::zero=0
   save ifirst,nyl,igmin

   !--select option according to reaction type requested
   if (mtd.eq.251.or.mtd.eq.252.or.mtd.eq.253) go to 400
   if (izap.gt.1.and.mfd.eq.3.and.mtd.eq.2) go to 800
   if (mfd.eq.3) go to 100
   if (mfd.eq.5) go to 200
   if (mfd.eq.12.or.mfd.eq.13) go to 100
   if (mfd.ge.10000000) go to 100
   if (mfd.eq.6.and.mtd.eq.2) go to 400
   if (mfd.eq.6.and.mtd.ge.221.and.mtd.le.250) go to 700
   if (mfd.eq.6.and.(mtd.ge.51.and.mtd.le.90)) go to 400
   if (mfd.eq.6.and.(mtd.ge.6.and.mtd.le.9)) go to 400
   if (mfd.eq.6) go to 200
   if (mfd.eq.16.or.mfd.eq.17) go to 300
   if (mfd.eq.8) go to 800
   if (mfd.eq.18) go to 810
   if (mfd.ge.21.and.mfd.le.26) go to 800
   if (mfd.ge.31.and.mfd.le.36.and.mtd.eq.2) go to 400
   if (mfd.ge.31.and.mfd.le.36.and.&
     (mtd.ge.51.and.mtd.lt.91)) go to 400
   if (mfd.ge.31.and.mfd.le.36.and.iverf.lt.6.and.&
     (mtd.ge.700.and.mtd.lt.719)) go to 400
   if (mfd.ge.31.and.mfd.le.36.and.iverf.lt.6.and.&
     (mtd.ge.720.and.mtd.lt.739)) go to 400
   if (mfd.ge.31.and.mfd.le.36.and.iverf.lt.6.and.&
     (mtd.ge.740.and.mtd.lt.759)) go to 400
   if (mfd.ge.31.and.mfd.le.36.and.iverf.lt.6.and.&
     (mtd.ge.760.and.mtd.lt.779)) go to 400
   if (mfd.ge.31.and.mfd.le.36.and.iverf.lt.6.and.&
     (mtd.ge.780.and.mtd.lt.799)) go to 400
   if (mfd.ge.31.and.mfd.le.36.and.iverf.ge.6.and.&
     (mtd.ge.600.and.mtd.lt.649)) go to 400
   if (mfd.ge.31.and.mfd.le.36.and.iverf.ge.6.and.&
     (mtd.ge.650.and.mtd.lt.699)) go to 400
   if (mfd.ge.31.and.mfd.le.36.and.iverf.ge.6.and.&
     (mtd.ge.700.and.mtd.lt.749)) go to 400
   if (mfd.ge.31.and.mfd.le.36.and.iverf.ge.6.and.&
     (mtd.ge.750.and.mtd.lt.799)) go to 400
   if (mfd.ge.31.and.mfd.le.36.and.iverf.ge.6.and.&
     (mtd.ge.800.and.mtd.lt.849)) go to 400
   if (mfd.ge.31.and.mfd.le.36.and.iverf.ge.6.and.&
     (mtd.ge.875.and.mtd.lt.891)) go to 400
   write(com,'(''do not know how to handle mf,mt: '',&
     &i2,'','',i3)') mfd,mtd
   call error('getff',com,' ')

   !--self-shielded group cross sections and ratios.
  100 continue
   enext=emax
   idisc=0
   en=enext
   nq=0
   ng=1
   if (mtd.eq.452.or.mtd.eq.455.or.mtd.eq.456) ng=2
   iglo=1
   yld=1
   if (mfd.eq.12.or.(mfd.ge.30000000.and.mfd.lt.40000000))&
     call getyld(e,en,idis,yld,matd,mfd,mtd,lfs,nend3)
   if (mfd.ge.20000000.and.mfd.lt.30000000)&
     call getyld(e,en,idis,yld,matd,mfd,mtd,isom,nend3)
   if (mtd.eq.452.or.mtd.eq.455.or.mtd.eq.456)&
     call getyld(e,en,idis,yld,matd,1,mtd,0,nend3)
   if (en.lt.enext*(1-small)) idisc=idis
   if (en.lt.enext*(1-small)) enext=en
   if (e.eq.zero) return
   do il=1,nl
      ff(il,1)=yld
      if (ng.eq.2) ff(il,2)=1
   enddo
   return

   !--neutron continuum transfer matrices
  200 continue
   enext=emax
   idisc=0
   en=enext
   nq=0
   yld=1
   if (e.eq.zero) ifirst=0
   if (e.gt.zero) ifirst=ifirst+1
   if (mfd.eq.5) go to 202
   if (e.lt.econst*(1-small).and.ifirst.eq.1) go to 202
   if (mtd.eq.16.or.mtd.eq.24.or.mtd.eq.26.or.mtd.eq.30) yld=2
   if (mtd.eq.11.or.mtd.eq.41) yld=2
   if (mtd.eq.17.or.mtd.eq.25) yld=3
   if (mtd.eq.42) yld=3
   if (mtd.eq.37) yld=4
   if ((mtd.ge.18.and.mtd.le.21).or.mtd.eq.38)&
     call getyld(e,en,idis,yld,matd,1,456,0,nend3)
   if (en.lt.enext*(1-small)) idisc=idis
   if (en.lt.enext*(1-small)) enext=en
   nle=nl
   lcd=1
   call getfle(e,en,idis,fle,nle,lcd,matd,4,mtd,nend2)
   if (e.eq.0.and.nle.eq.1) nl=1
   if (abs(en-enext).lt.enext*small.and.idis.gt.idisc) idisc=idis
   if (en.lt.enext*(1-small)) idisc=idis
   if (en.lt.enext*(1-small)) enext=en
  202 continue
   if (e.lt.econst*(1-small).and.ifirst.gt.1) go to 230
   nk=1
   if (mtd.eq.455) nk=ndelg
   if (e.eq.zero) allocate(sede(nk,ngn))
   call getsed(e,en,idis,sede,egn,ngn,nk,matd,5,mtd,nendf)
   if (abs(en-enext).lt.enext*small.and.idis.gt.idisc) idisc=idis
   if (en.lt.enext*(1-small)) idisc=idis
   if (en.lt.enext*(1-small)) enext=en
   if (e.gt.zero) go to 205
   if (econst.eq.zero) return
   ! find group where e-dependence begins
   ngn1=ngn+1
   do ig=1,ngn1
      if (econst.gt.egn(ig)*(1+small)) jconst=ig-1
   enddo
   econst=egn(jconst+1)
   if (jconst.le.2) econst=0
   return
   ! for e.gt.econst, return matrix
  205 continue
   if (mfd.eq.5) go to 260
   if (e.lt.econst*(1-small)) go to 230
   iglo=1
   ighi=1
   do ig=1,ngn
      temp=yld*sede(1,ig)
      if (temp.gt.0) ighi=ig
      do il=1,nl
         ff(il,ig)=temp*fle(il)
      enddo
   enddo
   ng=ighi-iglo+1
   return
   ! for e.lt.econst, return yield only
  230 continue
   if (ifirst.eq.1) go to 260
   iglo=1
   ng=1
   ff(1,1)=yld
   enext=econst
   idisc=1
   return
   ! spectrum calculation (e.lt.econst or mfd=5)
  260 continue
   iglo=0
   if (mfd.eq.5) iglo=1
   do 270 ig=1,ngn
   if (iglo.gt.0) go to 280
   do 275 ik=1,nk
   if (sede(ik,ig).gt.smin/1000.or.mtd.eq.18) go to 280
  275 continue
   go to 270
  280 continue
   if (iglo.eq.0) iglo=ig
   if (mfd.ne.5) ighi=ig
   do ik=1,nk
      jg=ig-iglo+1
      ff(ik,jg)=sede(ik,ig)
      if (mfd.eq.5) then
         if (ff(ik,jg).ne.zero) ighi=ig
      endif
   enddo
  270 continue
   if (iglo.eq.0) ighi=1
   if (iglo.eq.0) iglo=1
   ng=ighi-iglo+1
   return

   !--gamma production matrices
  300 continue
   enext=emax
   idisc=0
   nq=0
   if (mfd.eq.16) mft=12
   if (mfd.eq.17) mft=13
   ! initialize retrieval routines and assign storage.
   if (e.eq.zero) then
      call getgyl(e,en,idis,y,eccn,nyl,matd,mft,mtd,nend2)
      econ=eccn(1)
      if (en.lt.enext*(1-small)) enext=en
      allocate(gyle(nyl))
      allocate(eyle(nyl))
      allocate(gfle(nlg,nyl))
      call getgfl(e,en,idis,gfle,nl,nlg,nyl,matd,14,mtd,nend2)
      if (econ.eq.zero) then
         nk=1
         call getsed(e,en,idis,s,egg,ngg,nk,matd,15,mtd,nendf)
         if (en.lt.enext*(1-small)) enext=en
         allocate(sede(nk,ngg))
      endif
      ! find group containing lowest energy photon
      igmin=1
      do while (econ.ge.egg(igmin+1))
         igmin=igmin+1
      enddo
      if (econst.ne.zero) then
         ! find group where e dependence starts
         ngn1=ngn+1
         do ig=1,ngn1
            if (econst.gt.egn(ig)) jconst=ig-1
         enddo
         econst=egn(jconst+1)
         if (jconst.le.2) econst=0
      endif
      ifirst=0
      return
   endif
   ! retrieve factors of feed function.
   ifirst=ifirst+1
   call getgyl(e,en,idis,gyle,eyle,nyl,matd,mft,mtd,nend2)
   if (abs(en-enext).lt.enext*small.and.idis.gt.idisc) idisc=idis
   if (en.lt.enext*(1-small)) idisc=idis
   if (en.lt.enext*(1-small)) enext=en
   if (e.lt.econst*(1-small).and.ifirst.gt.1) go to 395
   call getgfl(e,en,idis,gfle,nl,nlg,nyl,matd,14,mtd,nend2)
   if (abs(en-enext).lt.enext*small.and.idis.gt.idisc) idisc=idis
   if (en.lt.enext*(1-small)) idisc=idis
   if (en.lt.enext*(1-small)) enext=en
   if (eyle(nyl).ne.zero) go to 330
   nk=1
   call getsed(e,en,idis,sede,egg,ngg,nk,matd,15,mtd,nendf)
   if (abs(en-enext).lt.enext*small.and.idis.gt.idisc) idisc=idis
   if (en.lt.enext*(1-small)) idisc=idis
   if (en.lt.enext*(1-small)) enext=en
   if (e.gt.zero) go to 330
   return
  330 continue
   iglo=igmin
   ighi=1
   do 360 ig=1,ngg
   epl=egg(ig)
   eph=egg(ig+1)
   do il=1,nl
      ff(il,ig)=0
   enddo
   do 350 ik=1,nyl
   if (gyle(ik).eq.zero) go to 350
   if (eyle(ik).eq.zero) go to 340
   if (eyle(ik).lt.epl*(1-small)) go to 350
   if (eyle(ik).ge.eph*(1-small)) go to 350
  340 continue
   tempo=gyle(ik)
   if (eyle(ik).eq.zero) tempo=tempo*sede(1,ig)
   do il=1,nl
      ff(il,ig)=ff(il,ig)+tempo*gfle(il,ik)
      if (ig.gt.ighi) ighi=ig
   enddo
  350 continue
  360 continue
   yld=0
   do ik=1,nyl
      yld=gyle(ik)+yld
   enddo
   ng=ighi-iglo+1
   if (e.ge.econst.or.ifirst.ne.1) then
      if (iglo.eq.1) return
      iglo1=iglo-1
      do ig=1,ng
         do il=1,nl
            ff(il,ig)=ff(il,ig+iglo1)
         enddo
      enddo
      return
   endif
   iglo=0
   ighi=0
   do ig=1,ngg
      ff(1,ig)=ff(1,ig)/yld
      if (ff(1,ig).ge.small) then
         if (iglo.eq.0) iglo=ig
         ighi=ig
      endif
   enddo
   ng=ighi-iglo+1
   do ig=1,ng
      ff(1,ig)=ff(1,iglo-1+ig)
   enddo
   return
  395 continue
   iglo=-1
   ng=1
   yld=0
   do ik=1,nyl
      yld=yld+gyle(ik)
   enddo
   ff(1,1)=yld
   enext=econst
   idisc=1
   return

   !--discrete channel scattering
   !--by the center-of-mass method
  400 continue
   call parts (mfd,mtd)
   call getdis(e,enext,idisc,yld,ff,nl,ng,iglo,nq,nend2,nlg)
   return

   !--angle-energy distribution matrices (mf6).
   !--special for thermal data (mt221 to 250).
  700 continue
   idisc=0
   call getaed(ff,e,enext,idisc,nl,nq,matd,mfd,mtd,npend2,nlg)
   yld=1
   ng=ngn
   if (e.eq.zero) return
   iglo=ngn
   ighi=1
   do ig=1,ngn
      do il=1,nl
         if (abs(ff(il,ig)).gt.small) then
            if (ig.lt.iglo) iglo=ig
            if (ig.gt.ighi) ighi=ig
         endif
      enddo
   enddo
   if (iglo.gt.ighi) iglo=ighi
   ! iglo=iglo-2
   iglo=1
   if (iglo.lt.1) iglo=1
   ng=ighi-iglo+1
   if (iglo.eq.1) return
   do ig=1,ng
      do il=1,nl
         ff(il,ig)=ff(il,iglo+ig-1)
      enddo
   enddo
   return

   !--energy-angle distributions for particle
   !--from endf-6 file 6
  800 continue
   call getmf6(ff,e,enext,idisc,yld,egn,ngn,nl,iglo,ng,nq,&
     matd,mfd,mtd,nend3,nlg)
   if ((mtd.ge.18.and.mtd.le.21).or.mtd.eq.38) then
      call getyld(e,en,idis,yld,matd,1,456,0,nend3)
      do ig=1,ng
         do il=1,nl
            ff(il,ig)=ff(il,ig)*yld
         enddo
      enddo
      if (en.lt.enext*(1-small)) then
         enext=en
         idisc=idis
      endif
   endif
   return

   !--energy-angle distributions for gamma
   !--from endf-6 file 6
  810 continue
   call getmf6(ff,e,enext,idisc,yld,egg,ngg,nl,iglo,ng,nq,&
     matd,mfd,mtd,nend3,nlg)
   return
   end subroutine getff

   subroutine parts(mfd,mtd)
   !-------------------------------------------------------------------
   ! Set up particles for reactions that use file 4.
   !-------------------------------------------------------------------
   use physics !get global physics and light particle mass constants
   use endf ! provides iverf
   ! externals
   integer::mfd,mtd

   law=0
   zap=1
   aprime=1
   if (iverf.lt.6) then
      if (mfd.eq.3) return
      if (mtd.eq.2.or.(mtd.ge.51.and.mtd.le.91)) then
         zap=1
         aprime=1
         if (mfd.ne.6) then
            aprime=awr
            law=4
         endif
      else if (mtd.ge.700.and.mtd.lt.720) then
         zap=1001
         aprime=pnratio
         if (mfd.ne.31) then
            aprime=awr+1-aprime
            law=4
         endif
      else if (mtd.ge.720.and.mtd.lt.740) then
         zap=1002
         aprime=dnratio
         if (mfd.ne.32) then
            aprime=awr+1-aprime
            law=4
         endif
      else if (mtd.ge.740.and.mtd.lt.760) then
         zap=1003
         aprime=tnratio
         if (mfd.ne.33) then
            aprime=awr+1-aprime
            law=4
         endif
      else if (mtd.ge.760.and.mtd.lt.780) then
         zap=2003
         aprime=hnratio
         if (mfd.ne.34) then
            aprime=awr+1-aprime
            law=4
         endif
      else if (mtd.ge.780.and.mtd.lt.800) then
         zap=2004
         aprime=anratio
         if (mfd.ne.35) then
            aprime=awr+1-aprime
            law=4
         endif
      endif
   else
      if (mfd.eq.3.and.mtd.ne.2) return
      if (mtd.eq.2) then
         zap=izap
         aprime=1
         if (mfd.ne.6) then
            aprime=awr
            law=4
         endif
      else if (mtd.ge.50.and.mtd.le.91) then
         zap=1
         aprime=1
         if (mfd.ne.6) then
            aprime=awr
            law=4
         endif
      else if (mtd.ge.600.and.mtd.lt.650) then
         zap=1001
         aprime=pnratio
         if (mfd.ne.31) then
            aprime=awr+1-aprime
            law=4
         endif
      else if (mtd.ge.650.and.mtd.lt.700) then
         zap=1002
         aprime=dnratio
         if (mfd.ne.32) then
            aprime=awr+1-aprime
            law=4
         endif
      else if (mtd.ge.700.and.mtd.lt.750) then
         zap=1003
         aprime=tnratio
         if (mfd.ne.33) then
            aprime=awr+1-aprime
            law=4
         endif
      else if (mtd.ge.750.and.mtd.lt.800) then
         zap=2003
         aprime=hnratio
         if (mfd.ne.34) then
            aprime=awr+1-aprime
            law=4
         endif
      else if(mtd.ge.800.and.mtd.lt.850) then
         zap=2004
         aprime=anratio
         if (mfd.ne.35) then
            aprime=awr+1-aprime
            law=4
         endif
      endif
   endif
   return
   end subroutine parts

   subroutine getmf6(ans,ed,enext,idisc,yld,eg,ng,nl,iglo,ng2,nq,&
     matd,mfd,mtd,nin,nlg)
   !-------------------------------------------------------------------
   ! Compute feed function from ENDF-6 File 6 data.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use endf   ! provides endf routines and variables
   use util   ! provides error,mess,skiprz
   ! externals
   integer::idisc,ng,nl,iglo,ng2,nq,matd,mfd,mtd,nin,nlg
   real(kr)::ed,enext,yld
   real(kr)::ans(nlg,*),eg(*)
   ! internals
   integer::mfn,nb,nw,lct,lct3,ik,nne,ne,int,nss
   integer::ie,ilo,jlo,jhi,ii,nn,nnn,langn,lepn,idis,jzap
   integer::nk,jzad,lang,lep,i,npsx,irr,npp,nmu,l1
   integer::j,iss,ip,ir,jgmax,jj,jg,ndlo,nplo,nclo,nphi,nchi
   integer::llo,lhi,iz,l,iy,max,nc,lf
   real(kr)::zad,elo,ehi,apsx,enow,eihi,ep,epnext,en
   real(kr)::pspmax,yldd,el,eh,e0,g0,e1,e2,test,pe,disc102
   real(kr)::val,fx,ex,cx,cxx,sum,rn,dx
   integer(kr)::nx,ncyc,n,ix
   integer,parameter::mxlg=65
   real(kr)::term(mxlg),terml(mxlg)
   character(60)::strng
   integer,parameter::maxss=500
   integer,parameter::nssm=9
   integer,dimension(nssm)::iyss,izss,jjss
   integer,dimension(maxss)::jloss
   real(kr),dimension(:),allocatable::tmp
   integer,parameter::ncmax=350
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::shade=1.1999e0_kr
   real(kr),parameter::step=1.2e0_kr
   real(kr),parameter::up=1.00001e0_kr
   real(kr),parameter::dn=0.99999e0_kr
   real(kr),parameter::eps=0.02e0_kr
   real(kr),parameter::zero=0
   real(kr),parameter::alight=5
   integer,parameter::ntmp=990000
   save nne,ne,int
   save jlo,elo,jhi,ehi,terml
   save pspmax,langn,lepn,disc102
   save idis,iyss,izss,jjss,jloss,nss,jzap,lct3,lct

   !--initialize
   if (ed.gt.zero) go to 200
   idis=0
   allocate(tmp(ntmp))
   max=ntmp
   mfn=6
   call findf(matd,mfn,mtd,nin)
   call contio(nin,0,0,tmp,nb,nw)
   nk=nint(tmp(5))
   lct=nint(tmp(4))
   lct3=lct
   izat=nint(tmp(1))
   zad=1
   if (mfd.eq.3) zad=izap
   if (mfd.eq.18) zad=0
   if (mfd.eq.18) lct=1
   if (mfd.eq.21) zad=1001
   if (mfd.eq.22) zad=1002
   if (mfd.eq.23) zad=1003
   if (mfd.eq.24) zad=2003
   if (mfd.eq.25) zad=2004
   if (mfd.eq.26) zad=9999
   jzad=nint(zad)
   iy=7
   nss=1
   jjss(1)=1
   iyss(1)=iy
   elo=emax
   ehi=0
   ik=0
  100 continue
   ik=ik+1
   if (ik.gt.nk.and.mfd.eq.18) go to 199
   if (ik.gt.nk)&
     call error('getmf6','desired particle not found.',' ')
   l=iy
   call tab1io(nin,0,0,tmp(l),nb,nw)
   l=l+nw
   do while (nb.ne.0)
      call moreio(nin,0,0,tmp(l),nb,nw)
      l=l+nw
   enddo
   disc102=0
   jzap=nint(tmp(iy))
   zap=jzap
   awp=tmp(iy+1)
   lip=nint(tmp(iy+2))
   law=nint(tmp(iy+3))
   if (jzap.eq.0.and.law.eq.2) then
      disc102=awp
      awp=0
   endif
   if (lct3.eq.3) then
      lct=1
      if (jzap.ge.1.and.jzap.le.2004) lct=2
   endif
   aprime=awp
   if (awrp.ne.0) aprime=awp/awrp
   if (jzap.eq.jzad) go to 130
   if (mfd.eq.26.and.jzap.gt.2004) go to 130
   call skip6(nin,0,0,tmp(l),law)
   go to 100
  130 continue

   !--for discrete recoil, back up to particle distribution,
   !--which is assumed to be the first subsection.
   if (law.eq.4) then
      call skiprz(nin,-2)
      call findf(matd,mfn,mtd,nin)
      l=1
      call contio(nin,0,0,tmp(l),nb,nw)
      l=l+nw
      call tab1io(nin,0,0,tmp(l),nb,nw)
      tmp(l)=zap
      tmp(l+1)=awp
      lf=nint(tmp(l+3))
      tmp(l+3)=law
      if (lf.eq.3) law=-4
      l=l+nw
      do while (nb.ne.0)
         call moreio(nin,0,0,tmp(l),nb,nw)
         l=l+nw
      enddo
   endif

   !--check law
  140 continue
   if (disc102.gt.zero) go to 195
   if (law.ge.2.and.law.le.5) go to 194
   if (law.eq.-4) go to 194
   if (law.gt.7) call error('getmf6','illegal law.',' ')

   !--read in data for law 1
   if (law.eq.1) then
      izss(nss)=l
      iz=l
      call tab2io(nin,0,0,tmp(l),nb,nw)
      lang=nint(tmp(l+2))
      lep=nint(tmp(l+3))
      ne=nint(tmp(l+5))
      int=nint(tmp(l+7))
      ! force unit base interpolation for smoother scattering source
      if (int.eq.2) int=22
      nn=0
      l=l+nw
      nnn=0
      do ie=1,ne
         ilo=l
         jlo=ilo
         ii=jjss(nss)+ie-1
         if (ii.gt.maxss) call error('getmf6',&
           'too many subsection energy points',' ')
         jloss(ii)=jlo
         call listio(nin,0,0,tmp(l),nb,nw)
         if (tmp(l+1).lt.elo) elo=tmp(l+1)
         if (tmp(l+1).gt.ehi) ehi=tmp(l+1)
         nn=nint(tmp(l+5))
         ncyc=nint(tmp(l+4)/nn)
         if (ie.gt.1.and.nn.ne.nnn.and.int.ge.11.and.int.le.20) then
            int=int+10
            tmp(iz+7)=int
            call mess('getmf6',&
              'bad grids for corresponding-point interpolation-',&
              'changing to unit-base interpolation')
         endif
         nnn=nn
         l=l+nw
         do while (nb.ne.0)
            if (l.gt.ntmp) call error('getmf6',&
              'storage exceeded',' ')
            call moreio(nin,0,0,tmp(l),nb,nw)
            l=l+nw
         enddo
         if (nn.eq.2.and.lep.eq.2) then
            j=0
            if (tmp(ilo+6+ncyc).lt.tmp(ilo+1)/(awr+1)**2) then
               j=1
               tmp(ilo+7)=0
               tmp(ilo+8)=tmp(ilo+1)/(awr+1)**2
               tmp(ilo+9)=2/tmp(ilo+8)
            elseif (tmp(ilo+7).ne.0.and.tmp(ilo+9).eq.0) then
               j=1
               tmp(ilo+9)=tmp(ilo+7)
               tmp(ilo+7)=0
            endif
            if (j.ne.0) then
               write(nsyso,'('' patching low-energy distribution at'',&
                 &1p,e10.3)')tmp(ilo+1)
            endif
         endif
         if (ismooth.gt.0.and.jzap.eq.1.and.lep.eq.1) then
            fx=.8409
            ex=40
            ncyc=nint(tmp(ilo+3))+2
            cx=tmp(ilo+6+ncyc)*tmp(ilo+7)
            nx=nint(tmp(ilo+4))
            n=nint(tmp(ilo+5))
            do while (n.gt.2)
               cxx=cx+tmp(ilo+7+ncyc)*(tmp(ilo+6+2*ncyc)-tmp(ilo+6+ncyc))
               if (abs(cxx/tmp(ilo+6+2*ncyc)**1.5&
                 -cx/tmp(ilo+6+ncyc)**1.5)&
                 .gt.cx/tmp(ilo+6+ncyc)**1.5/50) exit
               tmp(ilo+7)=(tmp(ilo+7)*tmp(ilo+6+ncyc)&
                 +tmp(ilo+7+ncyc)*(tmp(ilo+6+2*ncyc)&
                 -tmp(ilo+6+ncyc)))/tmp(ilo+6+2*ncyc)
               do ix=1,nx-2*ncyc
                  tmp(ilo+5+ix+ncyc)=tmp(ilo+5+ix+2*ncyc)
               enddo
               cx=cxx
               nx=nx-ncyc
               n=n-1
            enddo
            write(nsyso,'('' extending histogram as sqrt(E) below'',&
              &1p,e10.2,'' eV for E='',e10.2,'' eV'')')&
              tmp(ilo+6+ncyc),tmp(ilo+1)
            do while (tmp(ilo+6+ncyc).gt.ex)
               do ix=nx,1,-1
                  tmp(ilo+5+ncyc+ix)=tmp(ilo+5+ix)
                  tmp(ilo+5+ncyc+ix)=sigfig(tmp(ilo+5+ncyc+ix),7,0)
               enddo
               tmp(ilo+6+ncyc)=fx*tmp(ilo+6+2*ncyc)
               tmp(ilo+6+ncyc)=sigfig(tmp(ilo+6+ncyc),7,0)
               val=tmp(ilo+7)
               tmp(ilo+7)=sqrt(fx)*val
               tmp(ilo+7)=sigfig(tmp(ilo+7),7,0)
               tmp(ilo+7+ncyc)=(1-fx*sqrt(fx))*val/(1-fx)
               tmp(ilo+7+ncyc)=sigfig(tmp(ilo+7+ncyc),7,0)
               nx=nx+ncyc
               n=n+1
               tmp(ilo+4)=nx
               tmp(ilo+5)=n
            enddo
            l=ilo+6+nx
         else if (ismooth.gt.0.and.jzap.eq.1.and.lep.eq.2) then
             ncyc=nint(tmp(ilo+3))+2
             nx=nint(tmp(ilo+4))
             n=nint(tmp(ilo+5))
             write(nsyso,'('' extending lin-lin as sqrt(E) below'',&
               &1p,e10.2,'' eV for E='',e10.2,'' eV'')')&
               tmp(ilo+6+ncyc),tmp(ilo+1)
             ex=40
             fx=0.50
             nn=0
             cx=(tmp(ilo+6+ncyc)-tmp(ilo+6))*(tmp(ilo+7+ncyc)+tmp(ilo+7))/2
             cxx=0
             do i=1,ncyc-1
                tmp(ilo+6+i)=0
             enddo
             do while (tmp(ilo+6+ncyc).gt.ex)
                nn=nn+1
                do ix=nx,ncyc,-1
                   tmp(ilo+5+ncyc+ix)=tmp(ilo+5+ix)
                enddo
                tmp(ilo+6+ncyc)=fx*tmp(ilo+6+2*ncyc)
                tmp(ilo+6+ncyc)=sigfig(tmp(ilo+6+ncyc),6,0)
                do i=1,ncyc-1
                   tmp(ilo+6+ncyc+i)=tmp(ilo+6+2*ncyc+i)*&
                     sqrt(tmp(ilo+6+ncyc)/tmp(ilo+6+2*ncyc))
                enddo
                if (nn.gt.1) then
                   cxx=cxx+(tmp(ilo+6+2*ncyc)-tmp(ilo+6+ncyc))*&
                     (tmp(ilo+7+2*ncyc)+tmp(ilo+7+ncyc))/2
                endif
                nx=nx+ncyc
                n=n+1
                tmp(ilo+4)=nx
                tmp(ilo+5)=n
             enddo
             cxx=cxx+tmp(ilo+6+ncyc)*tmp(ilo+7+ncyc)/2
             dx=tmp(ilo+6+(nn+1)*ncyc)-tmp(ilo+6+nn*ncyc)
             rn=1
             if (cxx+dx*tmp(ilo+7+nn*ncyc)/2.ne.zero) then
                rn=(cx-dx*tmp(ilo+7+(nn+1)*ncyc)/2)/&
                 (cxx+dx*tmp(ilo+7+nn*ncyc)/2)
             endif
             do i=1,nn
                do j=1,ncyc-1
                   tmp(ilo+6+i*ncyc+j)=rn*tmp(ilo+6+i*ncyc+j)
                enddo
             enddo
             l=ilo+6+nx
         endif
         if (lct.eq.2.or.(lct.eq.3.and.awp.le.alight)) then
            call cm2lab(ilo,jlo,l,tmp,nl,lang,lep,max)
            do i=jlo,l
               tmp(ilo+i-jlo)=tmp(i)
            enddo
            l=ilo+l-jlo
         endif
      enddo

   !--read in data for law 6.  convert them to law 1.
   else if (law.eq.6) then
      call contio(nin,0,0,tmp(l),nb,nw)
      apsx=tmp(l)
      npsx=nint(tmp(l+5))
      tmp(l)=0
      tmp(l+1)=0
      tmp(l+2)=0
      tmp(l+3)=0
      lang=0
      lep=2
      lct=2
      irr=nint(tmp(iy+4))
      npp=nint(tmp(iy+5))
      enow=tmp(iy+6+2*irr)/shade
      pspmax=tmp(iy+4+2*irr+2*npp)
      izss(nss)=l
      l=l+8
      ne=0
      do while (enow.lt.pspmax)
         ne=ne+1
         ii=jjss(nss)+ne-1
         if (ii.gt.maxss) call error('getmf6',&
           'too many subsection energy points',' ')
         ilo=l
         jloss(ii)=ilo
         tmp(ilo)=apsx
         enow=enow*step
         if (enow.gt.pspmax) enow=pspmax
         tmp(ilo+1)=enow
         if (enow.lt.elo) elo=enow
         tmp(ilo+2)=0
         tmp(ilo+3)=0
         tmp(ilo+4)=0
         tmp(ilo+5)=npsx
         l=l+6
         call cm2lab(ilo,jlo,l,tmp,nl,lang,lep,max)
         do i=jlo,l
            tmp(ilo+i-jlo)=tmp(i)
         enddo
         l=ilo+l-jlo
      enddo
      i=izss(nss)
      lang=1
      tmp(i+2)=lang
      tmp(i+3)=lep
      tmp(i+4)=1
      tmp(i+5)=ne
      tmp(i+6)=ne
      tmp(i+7)=22
      law=1
      tmp(iy+3)=law

   !--read in data for law 7
   else if (law.eq.7) then
      izss(nss)=l
      call tab2io(nin,0,0,tmp(l),nb,nw)
      lang=nint(tmp(l+2))
      lep=nint(tmp(l+3))
      ne=nint(tmp(l+5))
      l=l+nw
      do ie=1,ne
         ilo=l
         jlo=ilo
         ii=jjss(nss)+ie-1
         if (ii.gt.maxss) call error('getmf6',&
           'too many subsection energy points',' ')
         jloss(ii)=ilo
         call tab2io(nin,0,0,tmp(l),nb,nw)
         if (tmp(l+1).lt.elo) elo=tmp(l+1)
         if (tmp(l+1).gt.ehi) ehi=tmp(l+1)
         nmu=nint(tmp(l+5))
         l=l+nw
         l1=l
         do i=1,nmu
            call tab1io(nin,0,0,tmp(l),nb,nw)
            l=l+nw
            do while (nb.ne.0)
               if (l.gt.ntmp) call error('getmf6',&
                 'storage exceeded',' ')
               call moreio(nin,0,0,tmp(l),nb,nw)
               l=l+nw
            enddo
         enddo
         call ll2lab(ilo,jlo,l,tmp,nl,max)
         do i=jlo,l
            tmp(ilo+i-jlo)=tmp(i)
         enddo
         l=ilo+l-jlo
      enddo
   endif

   !--check for multiple particle subsections
   if (ik.eq.nk) go to 195
   ik=ik+1
   call tab1io(nin,0,0,tmp(l),nb,nw)
   l1=l
   l=l+nw
   do while (nb.ne.0)
      if (l.gt.ntmp) call error('getmf6','storage exceeded',' ')
      call moreio(nin,0,0,tmp(l),nb,nw)
      l=l+nw
   enddo
   if (nint(tmp(l1)).ne.jzap) go to 195
   law=nint(tmp(l1+3))
   nss=nss+1
   if (nss.gt.nssm) call error('getmf6',&
     'too many subsections for one particle',' ')
   iyss(nss)=l1
   jjss(nss)=jjss(nss-1)+ne
   go to 140

   !--discrete two-body scattering
  194 continue
   izss(nss)=l
   call getdis(ed,elo,idisc,yldd,ans,nl,ng2,iglo,nq,nin,nlg)
   l=l+ncmax

   !--initialization complete
  195 continue
   if (nss.gt.1) call mess('getmf6',&
     'there are multiple subsections in mf6',&
     'for this emitted particle')
   nc=l-1
   allocate(ddmf6(nc))
   do i=1,nc
      ddmf6(i)=tmp(i)
   enddo
   deallocate(tmp)
   enext=elo
   idisc=0
   return

   !--special exit. no mf6 gammas.
  199 continue
   enext=-1
   idisc=0
   return

   !--normal entry
  200 continue
   do j=1,ng
      do l=1,nl
         ans(l,j)=0
      enddo
   enddo
   enext=emax
   yld=0
   iss=1
  207 continue
   iy=iyss(iss)
   iz=izss(iss)
   lip=nint(ddmf6(iy+2))
   law=nint(ddmf6(iy+3))
   if (disc102.gt.zero) go to 600
   if (law.ge.2.and.law.le.5) go to 500

   !--interpolate for fractional probability
   ip=2
   ir=1
   call terpa(pe,ed,eihi,idisc,ddmf6(iy),ip,ir)
   if (abs(eihi-enext).lt.enext*small.and.idisc.gt.idis) idis=idisc
   if (eihi.lt.enext*(1-small)) idis=idisc
   if (eihi.lt.enext*(1-small)) enext=eihi
   if (pe.eq.zero) go to 450
   yld=yld+pe
   lang=nint(ddmf6(iz+2))
   lep=nint(ddmf6(iz+3))
   lepn=lep
   if (law.eq.1.and.lct.eq.2) lepn=2
   if (law.eq.7) lepn=2
   langn=lang
   if (law.eq.1.and.lct.eq.2) langn=1
   if (law.eq.1.and.lct.eq.3.and.awp.lt.alight) langn=1
   if (law.eq.7) langn=1

   !--find jgmax
   jgmax=ng+1
   if (q.le.zero.and.jzap.ne.0) then
      jgmax=1
      do while (eg(jgmax+1).lt.ed*(1-small))
         jgmax=jgmax+1
      enddo
   endif

   !--find desired energy panel
   nne=nint(ddmf6(iz+5))
   int=nint(ddmf6(iz+7))
   ! force unit base interpolation for smoother scattering source
   if (int.eq.2) int=22
   jj=jjss(iss)
   do 230 ie=1,nne
      jlo=jloss(jj+ie-1)
      jhi=jloss(jj+ie)
      if (ed.ge.ddmf6(jlo+1)*(1-small).and.&
        ed.le.ddmf6(jhi+1)*(1+small)) then
         elo=ddmf6(jlo+1)
         ehi=ddmf6(jhi+1)
         go to 300
      endif
  230 continue
   go to 450

   !--set up loop over secondary energy
  300 continue
   ep=0
   call f6lab(ep,epnext,terml,nl,law,int,langn,lepn,&
     ddmf6(jlo),ddmf6(jhi),ed)
   jg=1

   !--find next energy grid point
   !--and retrieve the legendre components.
  320 continue
   en=epnext
   if (en.ge.emax*(1-small)) go to 420
   if (jg.lt.ng.and.ep.ge.eg(jg+1)*(1-small)) jg=jg+1
   if (jg.lt.ng.and.eg(jg+1).lt.en*(1-small)) en=eg(jg+1)
   call f6lab(en,epnext,term,nl,law,int,langn,lepn,ddmf6(jlo),&
     ddmf6(jhi),ed)

   !--trapazoidal integration over panel.
   do l=1,nl
      ans(l,jg)=ans(l,jg)+pe*(en-ep)*(term(l)+terml(l))/2
      terml(l)=term(l)
   enddo
   ep=en
   go to 320

   !--add in discrete energies for lab, if any.
  420 continue
   if (lct.ne.2) then
      if (lct.ne.3.or.awp.ge.alight) then
         ndlo=nint(ddmf6(jlo+2))
         nplo=nint(ddmf6(jlo+5))
         nclo=nint(ddmf6(jlo+4))/nplo
         nphi=nint(ddmf6(jhi+5))
         nchi=nint(ddmf6(jhi+4))/nphi
         if (ndlo.gt.0) then
            do i=1,ndlo
               llo=jlo+6+nclo*(i-1)
               lhi=jhi+6+nchi*(i-1)
               el=ddmf6(llo)
               eh=ddmf6(lhi)
               if (el.lt.zero) el=ed*awr/(awr+1)-el
               if (eh.lt.zero) eh=ed*awr/(awr+1)-eh
               call terp1(elo,el,ehi,eh,ed,e0,2)
               call terp1(elo,ddmf6(llo+1),ehi,ddmf6(lhi+1),ed,g0,2)
               do j=1,ng
                  e1=eg(j)
                  if (j.eq.1) e1=0
                  e2=eg(j+1)
                  if (j.eq.ng) e2=emax
                  if (e0.ge.e1*(1-small).and.e0.lt.e2*(1-small))&
                    ans(1,j)=ans(1,j)+pe*g0
               enddo
            enddo
         endif
      endif
   endif

   !--finished with energy loop.
  450 continue
   iglo=1
   ng2=ng
   nq=0
   if (ed.lt.ehi*(1-small).and.ehi.lt.enext*(1-small)) enext=ehi

   !--loop back to process additional particles, if any.
   iss=iss+1
   if (iss.le.nss) go to 207
   go to 700

   !--laws 2 thru 5 -- discrete two-body scattering
  500 continue
   ilo=iz
   call getdis(ed,enext,idisc,yldd,ans,nl,ng2,iglo,nq,nin,nlg)
   return

   !--discrete relativistic capture gamma
  600 continue
   call gam102(ans,ed,enext,disc102,law,nl,iglo,ng2,nq)
   return

   !--check normalization
   !--move illegal upscatters to ingroup
  700 continue
   if (ed.le.zero) return
   if (yld.eq.zero) return
   test=0
   jj=jgmax-iglo+1
   do i=1,ng2
      test=test+ans(1,i)
      if (i.gt.jj.and.nint(zap).eq.izap.and.q.le.zero) then
         do l=1,nl
            ans(l,jj)=ans(l,jj)+ans(l,i)
            ans(l,i)=0
         enddo
      endif
   enddo
   if (test.eq.zero) return
   if (abs(test-yld).gt.eps.and.(ed.lt.up*elo.or.ed.gt.dn*ehi))&
     write(nsyso,'(''     normalization '',1p,e12.4,1p,e15.6)')&
     ed,test
   test=yld/test
   do l=1,nl
      do i=1,ng2
         ans(l,i)=test*ans(l,i)
      enddo
   enddo
   return
   end subroutine getmf6

   subroutine cm2lab(inow,jnow,lnow,c,nl,lang,lep,max)
   !-------------------------------------------------------------------
   ! Convert the cm distribution starting at inow into the laboratory
   ! frame.  The lab distribution starts at jnow.
   !-------------------------------------------------------------------
   use endf ! provides mth
   use util ! provides error,mess
   ! externals
   integer::inow,jnow,lnow,nl,lang,lep,max
   real(kr)::c(*)
   ! internals
   integer::l,i,iflag,nw,np
   real(kr)::sum,dy,da,xm,ym,test,dely,epnext,epn,rat
   integer,parameter::imax=15
   integer,parameter::nlmax=65
   real(kr)::term(nlmax),x(imax),y(imax,nlmax)
   character(60)::strng
   real(kr),parameter::tol=0.01e0_kr
   real(kr),parameter::eps=1.e-8_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::zero=0

   !--jnow points to the new lab distribution
   np=nint(c(inow+4))
   jnow=lnow
   c(jnow)=0
   c(jnow+1)=c(inow+1)
   lnow=jnow+6

   !--prime the stack with first point
   sum=0
   x(2)=0
   call f6cm(x(2),epnext,term,nl,lang,lep,c(inow))
   do l=1,nl
      y(2,l)=term(l)
   enddo

   !--prime the stack with the top of the next panel
   do while (epnext.lt.emax*(1-small))
      x(1)=epnext
      call f6cm(x(1),epnext,term,nl,lang,lep,c(inow))
      do l=1,nl
         y(1,l)=term(l)
      enddo

      !--adaptive integration over panel for cm.
      !--this gives a good energy grid in the lab.
      i=2
      dy=(y(1,1)+y(2,1))/2
      do while (i.gt.1)

         !--test for convergence.
         iflag=0
         if (i.gt.1.and.i.lt.imax) then
            da=(x(i-1)-x(i))*(y(i-1,1)+y(i,1))/2
            rat=2
            if (x(i).ne.zero) rat=x(i-1)/x(i)
            if (rat.gt.1+tol) then
               if (da.ge.eps) then
                  xm=(x(i-1)+x(i))/2
                  if (xm.ne.x(i-1).and.xm.ne.x(i)) then
                     epn=x(i)
                     call f6cm(xm,epn,term,nl,lang,lep,c(inow))
                        do l=1,nl
                        ym=(y(i-1,l)+y(i,l))/2
                           test=l*tol*(abs(ym))
                        dely=abs(term(l)-ym)
                        if (dely.gt.test) iflag=1
                     enddo
                  endif
               endif
            endif
         endif

         !--fails tests.
         !--add midpoint to stack and continue.
         if (iflag.eq.1) then
            i=i+1
            x(i)=x(i-1)
            do l=1,nl
               y(i,l)=y(i-1,l)
            enddo
            x(i-1)=xm
            do l=1,nl
               y(i-1,l)=term(l)
            enddo

         !--passes tests.
         !--take top point off of the stack.
         else
            c(lnow)=x(i)
            do l=1,nl
               c(lnow+l)=y(i,l)
            enddo
            lnow=lnow+1+nl
            if (lnow+nl+1.gt.max)&
              &call error('cm2lab','storage exceeded.',' ')
            sum=sum+(x(i-1)-x(i))*(y(i-1,1)+y(i,1))/2
            i=i-1
         endif
      enddo

      !--continue looping over panels.
      x(2)=x(1)
      do l=1,nl
         y(2,l)=y(1,l)
      enddo
   enddo
   if (abs(sum-1).gt.tol) then
      write(strng,'(''mt='',i3,'' e='',1p,e10.3,'' lab sum='',0p,f6.3)')&
        mth,c(inow+1),sum
      call mess('cm2lab','lab normalization problem ',strng)
   endif

   !--linearization and conversion complete.
   nw=lnow-jnow-6
   np=nw/(nl+1)
   c(jnow+4)=nw
   c(jnow+5)=np
   c(jnow+2)=0
   c(jnow+3)=nl-1
   return
   end subroutine cm2lab

   subroutine f6cm(ep,epnext,term,nl,lang,lep,cnow)
   !-------------------------------------------------------------------
   ! Compute the Legendre coefficients for the double-differential
   ! cross section for secondary energy ep (in the lab system)
   ! using the ENDF-6 File 6 data in the cm system located
   ! in the array cnow.  Call with ep=0 to initialize.
   ! Thereafter, ep can be requested in any order.
   !-------------------------------------------------------------------
   use util ! provides error
   use mathm ! provides legndr
   use endf
   ! externals
   integer::nl,lang,lep
   real(kr)::ep,epnext,term(*),cnow(*)
   ! internals
   integer::ndnow,npnow,ncnow,i,l,na,j,iflag,ll
   real(kr)::xc,elmax,epnn,epmax,w,xx,cc,c,umin,da
   real(kr)::dm,umax,un,u,yy,epp,ym,dy,wmin,wmax
   real(kr)::eps,e,t,epn,epm,epx,s,dx,epnxt,f,test
   integer,parameter::imax=10
   integer,parameter::mxlg=65
   real(kr)::p(mxlg)
   real(kr)::x(imax),y(imax,mxlg),yt(mxlg)
   real(kr),parameter::tol=0.005e0_kr
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::check=.99999e0_kr
   real(kr),parameter::rdn=.999995e0_kr
   real(kr),parameter::rup=1.000005e0_kr
   real(kr),parameter::amil=1.e-6_kr
   real(kr),parameter::athou=1.e-3_kr
   real(kr),parameter::dn=.999999e0_kr
   real(kr),parameter::up=1.000001e0_kr
   real(kr),parameter::step=.05e0_kr
   real(kr),parameter::one=1.e0_kr
   real(kr),parameter::tiny=1.e-4_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::zero=0
   save na,eps,xc,ndnow,npnow,ncnow,elmax,e,epmax

   !--initialize.
   if (nl.gt.mxlg) call error('f6cm','nl>mxlg',' ')
   if (ep.eq.zero) then
      e=cnow(2)
      ndnow=nint(cnow(3))
      npnow=nint(cnow(6))
      ncnow=nint(cnow(5))/npnow
      xc=aprime/(awr+1)**2
      epnext=emax
      elmax=0
      epnn=emax
      epmax=0
      if (npnow.ne.ndnow) then
         if (lang.gt.0) t=f6ddx(ep,epn,epm,w,cnow,lang,lep)
         if (lang.eq.0) t=f6psp(ep,epn,epm,w,e,cnow)
         if (epn.lt.epnn*(1-small)) epnn=epn
         if (epm.gt.epmax*(1+small)) epmax=epm
         elmax=e*(sqrt(epmax/e)+sqrt(xc))**2
         elmax=up*elmax
         epn=step*elmax
         if(epmax*up.lt.xc*e) epn=e*(sqrt(epmax/e)-sqrt(xc))**2
         if (epn.lt.epnext) epnext=epn
         eps=tiny/elmax
         epnext=dn*epnext
      endif
      if (ndnow.ne.0) then
         t=f6dis(0,epn,epm,w,cnow,lang)
         do i=1,ndnow
            t=f6dis(i,epn,epm,w,cnow,lang)
            if (epn.lt.epnn*(1-small)) epnn=epn
            if (epn.gt.epmax*(1+small)) epmax=epn
            epx=dn*e*(sqrt(epn/e)-sqrt(xc))**2
            if (epx.lt.epnext*(1-small)) epnext=epx
            epx=up*e*(sqrt(epn/e)+sqrt(xc))**2
            if (epx.gt.elmax*(1+small)) elmax=epx
         enddo
      endif
      na=nl-1
      do l=1,nl
         term(l)=0
      enddo
      return
   endif

   !--normal entry.
   do l=1,nl
      term(l)=0
   enddo
   if (ndnow.eq.npnow) go to 330
   xx=ep/e
   cc=xc/xx
   c=sqrt(cc)
   umin=(1+cc-epmax/ep)/(2*c)
   if (umin.lt.check) go to 215
   if (c.le.one) go to 420
   if (c.gt.one) go to 320
  215 continue
   if (umin.lt.-one) umin=-one
   dm=(1-umin)*athou+amil
   epnn=emax

   !--for the continuum part of the distribution,
   !--integrate over the trajectory ep=constant
   !--using adaptive integration on lab cosine.
   umax=1
   un=umax
   x(2)=un
   do while (x(2).gt.umin)

      !--load the boundaries for this cosine panel.
      j=2
      do while (j.eq.2)
         u=un
         yy=1+cc-2*c*u
         if (yy.lt.amil) then
            yy=amil
            c=u-athou
            cc=c*c
         endif
         epp=yy*ep
         w=(u-c)/sqrt(yy)
         if (lang.gt.0) s=f6ddx(epp,epn,epm,w,cnow,lang,lep)
         if (lang.eq.0) s=f6psp(epp,epn,epm,w,e,cnow)
         if (u.eq.umax.and.epn.lt.epnn) epnn=epn
         call legndr(u,p,na)
         j=1
         if (u.eq.umax) j=2
         x(j)=u
         do l=1,nl
            y(j,l)=p(l)*s/sqrt(yy)
         enddo
         un=u-(epn-epp)/(2*c*ep)
         if (un.gt.u-tiny) un=u-tiny
         if (un.lt.umin+tiny/10) un=umin
      enddo

      !--do an adaptive reconstruction for this panel.
      i=2
      do while (i.gt.1)

         !--test for convergence.
         iflag=0
         if (i.gt.1.and.i.lt.imax) then
            dx=x(i)-x(i-1)
            if (dx.ge.dm) then
               da=dx*(y(i,1)+y(i-1,1))/2
               if (da.ge.eps) then
                  u=(x(i-1)+x(i))/2
                  yy=1+cc-2*c*u
                  epp=yy*ep
                  w=(u-c)/sqrt(yy)
                  if (lang.gt.0) s=f6ddx(epp,epn,epm,w,cnow,lang,lep)
                  if (lang.eq.0) s=f6psp(epp,epn,epm,w,e,cnow)
                  call legndr(u,p,na)
                  do l=1,nl
                     yt(l)=p(l)*s/sqrt(yy)
                  enddo
                  do l=1,nl
                     ym=(y(i-1,l)+y(i,l))/2
                     dy=abs(yt(l)-ym)
                     if (dy.gt.l*tol*abs(ym)+eps) iflag=1
                  enddo
               endif
            endif
         endif

         !--fails test.
         !--add midpoint to stack.
         if (iflag.eq.1) then
            i=i+1
            x(i)=x(i-1)
            do l=1,nl
               y(i,l)=y(i-1,l)
            enddo
            x(i-1)=u
            do l=1,nl
               y(i-1,l)=yt(l)
            enddo

         !--passes test.
         !--take top point off of stack.
         else
            do l=1,nl
               term(l)=term(l)+(x(i)-x(i-1))*(y(i,l)+y(i-1,l))/2
            enddo
            i=i-1
         endif
      enddo

      !--continue loop over cosine panels
      x(2)=x(1)
      do l=1,nl
         y(2,l)=y(1,l)
      enddo
   enddo

   !--select next energy point
   epnxt=e*(sqrt(epnn/e)-sqrt(xc))**2
   if (epnxt.le.ep*(1+small)) epnxt=e*(sqrt(epnn/e)+sqrt(xc))**2
   if (epnxt.le.ep*(1+small)) epnxt=3*ep/2
   if (epnxt.gt.elmax*(1+small).and.elmax.gt.ep*(1+small))&
     epnxt=elmax
   go to 380
  320 continue
   epnxt=(epnext+elmax)/2
   go to 380

   !--check for contributions from delta functions
  330 continue
   if (ndnow.ne.0) then
      epnxt=emax
      do i=1,ndnow
         ll=7+ncnow*(i-1)
         epp=cnow(ll)
         w=(ep-epp-xc*e)/(2*sqrt(xc*e)*sqrt(epp))
         u=(xc*e+ep-epp)/(2*sqrt(xc*e)*sqrt(ep))
         wmin=-1
         wmax=1
         if (w.lt.wmin) then
            epn=e*(sqrt(epp/e)-sqrt(xc))**2
            if (ep.lt.check*epn) then
               epn=rdn*epn
            else
               epn=rup*epn
            endif
            if (epn.lt.epnxt*(1-small)) epnxt=epn
         else if (w.le.wmax) then
            call legndr(u,p,na)
            f=f6dis(i,epn,epm,w,cnow,lang)
            f=f/(2*sqrt(xc*e*epp))
            epn=e*(sqrt(epp/e)+sqrt(xc))**2
            if (ep.lt.check*epn) then
               epn=rdn*epn
            else
               epn=rup*epn
            endif
            if (epn.lt.epnxt*(1-small)) epnxt=epn
            do l=1,nl
               term(l)=term(l)+f*p(l)
            enddo
         endif
      enddo
   endif

   !--remove small, inaccurate moments
  380 continue
   if (nl.gt.1) then
      test=tol*abs(term(1))
      do l=2,nl
         if (abs(term(l)).lt.test) term(l)=0
      enddo
   endif

   !--select epnext
   epnext=epnxt
   return
  420 continue
   epnext=emax
   return
   end subroutine f6cm

   real(kr) function f6ddx(ep,epnext,epmax,w,cnow,lang,lep)
   !-------------------------------------------------------------------
   ! Retrieve the double differential cross section for
   ! secondary energy ep and cosine w (in the cm system) from a
   ! subsection in ENDF-6 File 6 format stored in cnow.  Call with
   ! ep=0 to initialize for each new e. Thereafter, ep can be
   ! requested in any order.  Only the continuum part is
   ! returned.  Delta functions, if given, must be handled
   ! separately.
   !-------------------------------------------------------------------
   use mathm ! provides legndr
   use endf  ! provides terp1
   use util  ! error,mess
   ! externals
   integer::lang,lep
   real(kr)::ep,epnext,epmax,w,cnow(*)
   ! internals
   integer::nl,inow,lnow,mnow,ncnow,na,ndnow,npnow,idone,illdef
   integer::l,iza2,int,ii,jj,lll,ia
   real(kr)::efirst,enow,t,eplast,s,r,aa,ss,bb,sa,tii,tjj,tt
   real(kr)::x1,x2,y1,y2
   integer,parameter::mxlg=65
   real(kr)::p(mxlg)
   real(kr),parameter::up=1.00001e0_kr
   real(kr),parameter::dn=.99999e0_kr
   real(kr),parameter::off=.999995e0_kr
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::step=0.05e0_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::tomev=1.e-6_kr
   real(kr),parameter::half=0.5e0_kr
   real(kr),parameter::zero=0
   save enow,efirst,nl,inow,lnow,mnow,ncnow,na,illdef

   !--initialize
   if (ep.eq.zero) then
      enow=cnow(2)
      ndnow=nint(cnow(3))
      npnow=nint(cnow(6))
      ncnow=nint(cnow(5))/npnow
      mnow=7+(npnow-1)*ncnow
      inow=7+ncnow*ndnow
      lnow=inow
      if (cnow(lnow).le.zero) lnow=lnow+ncnow
      epnext=cnow(lnow)
      epmax=cnow(mnow)
      if (epnext.gt.step*epmax) epnext=step*epmax
      nl=ncnow-1
      if (nl.gt.mxlg) call error('f6ddx','nl>mxlg',' ')
      na=nl-1
      efirst=0
      f6ddx=0
      illdef=0
      return
   endif

   !--use analytic form below efirst.
   if (ep.lt.efirst*(1-small)) then
      t=0

   !--compute the result for higher energies
   else

      !--check for desired interpolation range
      idone=0
      do while (idone.eq.0)
         epnext=cnow(lnow)
         if (ep.lt.off*epnext) then
            eplast=cnow(lnow-ncnow)
            if (ep.ge.off*eplast) then
               idone=1
            else
               if (lnow.le.inow+ncnow) then
                  idone=1
               else
                  lnow=lnow-ncnow
               endif
            endif
         else if (lnow.ge.mnow) then
            epnext=emax
            idone=1
         else
            lnow=lnow+ncnow
         endif
      enddo

      !--legendre coefficients
      if (lang.eq.1) then
         call legndr(w,p,na)
         t=0
         do l=1,nl
            if (l.le.ncnow-1.and.ep.ge.cnow(inow)*(1-small).and.&
              ep.le.cnow(mnow)*(1+small)) then
               lll=lnow+l
               x1=cnow(lnow-ncnow)
               x2=cnow(lnow)
               y1=cnow(lll-ncnow)
               y2=cnow(lll)
               if (x1.eq.x2) then
                  tt=y1
               else
                  call terp1(x1,y1,x2,y2,ep,tt,lep)
               endif
               t=t+(2*l-1)*tt*p(l)/2
            endif
         enddo
         if (t.lt.zero) t=0

      !--kalbach systematics
      else if (lang.eq.2) then
         s=0
         r=0
         if (ep.ge.cnow(inow)*(1-small).and.&
           ep.le.cnow(mnow)*(1+small)) then
            x1=cnow(lnow-ncnow)
            x2=cnow(lnow)
            if (x1.eq.x2.and.lep.gt.1) then
               x2=sigfig(x2,6,1)
               if (illdef.eq.0) then
                  call mess('f6ddx',&
                    'vertical segment(s) in distribution',&
                    'y(x) is ill defined')
                  illdef=1
               endif
            endif
            y1=cnow(lnow-ncnow+1)
            y2=cnow(lnow+1)
            call terp1(x1,y1,x2,y2,ep,s,lep)
            y1=cnow(lnow-ncnow+2)
            y2=cnow(lnow+2)
            call terp1(x1,y1,x2,y2,ep,r,lep)
         endif
         ! kalbach-86 (obsolete kalbach-mann coding has been deleted)
         iza2=nint(zap)
         if (na.eq.2) then
            y1=cnow(lnow-ncnow+3)
            y2=cnow(lnow+3)
            call terp1(x1,y1,x2,y2,ep,aa,lep)
         else
            aa=bach(izap,iza2,izat,enow,ep)
         endif

         t=aa*(cosh(aa*w)+r*sinh(aa*w))/(2*sinh(aa))
         t=t*s
         if (t.lt.zero) t=0

      !--tabulated angular distribution
      else
         if (lang.lt.11.or.lang.gt.15)&
           call error('f6ddx','illegal lang.',' ')
         int=lang-10
         ii=lnow-ncnow
         jj=lnow
         x1=cnow(ii)
         x2=cnow(jj)
         if (x1.eq.x2.and.lep.gt.1) then
            x2=sigfig(x2,6,1)
            if (illdef.eq.0) then
               call mess('f6ddx',&
                 'vertical segment(s) in distribution',&
                 'y(x) is ill defined')
               illdef=1
            endif
         endif
         y1=cnow(ii+1)
         y2=cnow(jj+1)
         call terp1(x1,y1,x2,y2,ep,s,lep)
         tii=0
         do ia=1,na,2
            if (w.ge.cnow(ii+1+ia).and.w.le.cnow(ii+3+ia)) then
               call terp1(cnow(ii+1+ia),cnow(ii+2+ia),cnow(ii+3+ia),&
                 cnow(ii+4+ia),w,tii,int)
            endif
         enddo
         tjj=0
         do ia=1,na,2
            if (w.ge.cnow(jj+1+ia).and.w.le.cnow(jj+3+ia)) then
               call terp1(cnow(jj+1+ia),cnow(jj+2+ia),cnow(jj+3+ia),&
                 cnow(jj+4+ia),w,tjj,int)
            endif
         enddo
         call terp1(x1,tii,x2,tjj,ep,t,lep)
         t=t*s
      endif
   endif

   !--return final value
   if (lep.gt.1) then
      f6ddx=t
   else if (abs(epnext-emax).lt.emax*small) then
      f6ddx=t
   else if (ep.ge.dn*dn*epnext) then
      epnext=up*epnext
      f6ddx=t
   else
      epnext=dn*epnext
      f6ddx=t
   endif
   return
   end function f6ddx

   real(kr) function f6dis(i,epnext,epmax,w,cnow,lang)
   !-------------------------------------------------------------------
   ! Retrieve the double differential cross section for
   ! discrete-energy i and cosine w (in the cm system) from a
   ! subsection in ENDF-6 File 6 format stored in cnow.  Call with
   ! i=0 to initialize for each new e. Thereafter, i can be
   ! requested in any order.
   !-------------------------------------------------------------------
   use mathm ! provides legndr
   use util  ! provides error
   use endf  ! provides terp1
   ! externals
   integer::i,lang
   real(kr)::epnext,epmax,w,cnow(*)
   ! internals
   integer::nl,inow,lnow,mnow,ncnow,na,ndnow,npnow,l
   integer::iza2,int,ia
   real(kr)::enow,t,s,r,ep,aa,ss,bb,sa,tt
   integer,parameter::mxlg=65
   real(kr)::p(mxlg)
   real(kr),parameter::tomev=1.e-6_kr
   real(kr),parameter::half=0.5e0_kr
   real(kr),parameter::zero=0.e0_kr
   save nl,inow,lnow,mnow,ncnow,na,enow

   !--initialize
   if (i.eq.0) then
      enow=cnow(2)
      ndnow=nint(cnow(3))
      npnow=nint(cnow(6))
      ncnow=nint(cnow(5))/npnow
      inow=7
      lnow=inow
      mnow=lnow+(ndnow-1)*ncnow
      epnext=cnow(lnow)
      epmax=cnow(mnow)
      nl=ncnow-1
      na=nl-1
      f6dis=0
      return
   endif

   !--select desired discrete line
   inow=7+(i-1)*ncnow
   epnext=cnow(inow)

   !--legendre coefficients
   if (lang.eq.1) then
      call legndr(w,p,na)
      t=0
      do l=1,nl
         if (l.le.ncnow-1) then
            t=t+(2*l-1)*p(l)*cnow(inow+1+l)/2
         endif
      enddo
      t=t*cnow(inow+1)
      if (t.lt.zero) t=0

   !--kalbach systematics
   else if (lang.eq.2) then
      s=cnow(inow+1)
      r=cnow(inow+2)
      ! kalbach-86 (obsolete kalbach-mann coding has been deleted)
      iza2=nint(zap)
      if(na.eq.2) then
         aa=cnow(inow+3)
      else
         aa=bach(izap,iza2,izat,enow,epnext)
      endif
      t=aa*(cosh(aa*w)+r*sinh(aa*w))/(2*sinh(aa))
      t=t*s
      if (t.lt.zero) t=0

   !--tabulated angular distribution
   else if (lang.ge.11.and.lang.le.15) then
      int=lang-10
      t=0
      do ia=1,na-2,2
         if (w.ge.cnow(inow+1+ia).and.w.le.cnow(inow+3+ia)) then
            call terp1(cnow(inow+1+ia),cnow(inow+2+ia),&
              cnow(inow+3+ia),cnow(inow+4+ia),w,t,int)
         endif
      enddo
      t=t*cnow(inow+1)

   !--illegal lang
   else
      call error('f6dis','illegal lang.',' ')
   endif

   !--return final value
   f6dis=t
   return
   end function f6dis

   real(kr) function bach(iza1i,iza2,izat,e,ep)
   !-------------------------------------------------------------------
   ! Compute the Kalbach a parameter.
   !-------------------------------------------------------------------
   use util    ! provides error
   use physics ! provides amassn,amu,ev,clight
   ! externals
   integer::iza1i,iza2,izat
   real(kr)::e,ep
   ! internals
   integer::iza1,iza,na
   real(kr)::emc2,aa,za,ac,ab,zb,zc,sa,sb,nc,nb
   real(kr)::ecm,ea,eb,x1,x3,fa,fb,bb,fact,test
   character(60)::strng
   real(kr),parameter::third=.333333333e0_kr
   real(kr),parameter::twoth=.666666667e0_kr
   real(kr),parameter::fourth=1.33333333e0_kr
   real(kr),parameter::c1=15.68e0_kr
   real(kr),parameter::c2=-28.07e0_kr
   real(kr),parameter::c3=-18.56e0_kr
   real(kr),parameter::c4=33.22e0_kr
   real(kr),parameter::c5=-0.717e0_kr
   real(kr),parameter::c6=1.211e0_kr
   real(kr),parameter::s2=2.22e0_kr
   real(kr),parameter::s3=8.48e0_kr
   real(kr),parameter::s4=7.72e0_kr
   real(kr),parameter::s5=28.3e0_kr
   real(kr),parameter::brk1=130.e0_kr
   real(kr),parameter::brk2=41.e0_kr
   real(kr),parameter::half=0.5e0_kr
   real(kr),parameter::b1=0.04e0_kr
   real(kr),parameter::b2=1.8e-6_kr
   real(kr),parameter::b3=6.7e-7_kr
   real(kr),parameter::d1=9.3e0_kr
   real(kr),parameter::tomev=1.e-6_kr
   real(kr),parameter::zero=0
   emc2=tomev*amassn*amu*clight*clight/ev

   iza1=iza1i
   if (iza1i.eq.0) iza1=1
   iza=izat
   if (iza.eq.6000) iza=6012
   if (iza.eq.12000) iza=12024
   if (iza.eq.14000) iza=14028
   if (iza.eq.16000) iza=16032
   if (iza.eq.17000) iza=17035
   if (iza.eq.19000) iza=19039
   if (iza.eq.20000) iza=20040
   if (iza.eq.22000) iza=22048
   if (iza.eq.23000) iza=23051
   if (iza.eq.24000) iza=24052
   if (iza.eq.26000) iza=26056
   if (iza.eq.28000) iza=28058
   if (iza.eq.29000) iza=29063
   if (iza.eq.31000) iza=31069
   if (iza.eq.40000) iza=40090
   if (iza.eq.42000) iza=42096
   if (iza.eq.48000) iza=48112
   if (iza.eq.49000) iza=49115
   if (iza.eq.50000) iza=50120
   if (iza.eq.63000) iza=63151
   if (iza.eq.72000) iza=72178
   if (iza.eq.74000) iza=74184
   if (iza.eq.82000) iza=82208
   aa=mod(iza,1000)
   if (aa.eq.zero) then
      write(strng,'(''dominant isotope not known for '',i8)') iza
      call error('bach',strng,' ')
   endif
   za=int(iza/1000)
   ac=aa+mod(iza1,1000)
   zc=za+int(iza1/1000)
   ab=ac-mod(iza2,1000)
   zb=zc-int(iza2/1000)
   na=nint(aa-za)
   nb=nint(ab-zb)
   nc=nint(ac-zc)
   sa=c1*(ac-aa)&
     +c2*((nc-zc)**2/ac-(na-za)**2/aa)&
     +c3*(ac**twoth-aa**twoth)&
     +c4*((nc-zc)**2/ac**fourth-(na-za)**2/aa**fourth)&
     +c5*(zc**2/ac**third-za**2/aa**third)&
     +c6*(zc**2/ac-za**2/aa)
   if (iza1.eq.1002) sa=sa-s2
   if (iza1.eq.1003) sa=sa-s3
   if (iza1.eq.2003) sa=sa-s4
   if (iza1.eq.2004) sa=sa-s5
   sb=c1*(ac-ab)&
     +c2*((nc-zc)**2/ac-(nb-zb)**2/ab)&
     +c3*(ac**twoth-ab**twoth)&
     +c4*((nc-zc)**2/ac**fourth-(nb-zb)**2/ab**fourth)&
     +c5*(zc**2/ac**third-zb**2/ab**third)&
     +c6*(zc**2/ac-zb**2/ab)
   if (iza2.eq.1002) sb=sb-s2
   if (iza2.eq.1003) sb=sb-s3
   if (iza2.eq.2003) sb=sb-s4
   if (iza2.eq.2004) sb=sb-s5
   ecm=aa*e/ac
   ea=ecm*tomev+sa
   eb=tomev*ep*ac/ab+sb
   x1=eb
   if (ea.gt.brk1) x1=brk1*eb/ea
   x3=eb
   if (ea.gt.brk2) x3=brk2*eb/ea
   fa=1
   if (iza1.eq.2004) fa=0
   fb=1
   if (iza2.eq.1) fb=half
   if (iza2.eq.2004) fb=2
   bb=b1*x1+b2*x1**3+b3*fa*fb*x3**4
   if (iza1i.eq.0) then
      fact=d1/sqrt(ep*tomev)
      test=1
      if (fact.lt.test) fact=test
      test=4
      if (fact.gt.test) fact=test
      bb=bb*sqrt((tomev*e)/(2*emc2))*fact
   endif
   bach=bb
   return
   end function bach

   subroutine ll2lab(inow,jnow,lnow,c,nl,max)
   !-------------------------------------------------------------------
   ! Change the angle-energy lab distribution in law 7 format into
   ! a laboratory Legendre coefficient representation in law 1 format.
   !-------------------------------------------------------------------
   use util  ! provides error
   use endf  ! provides terpa
   use mathm ! provides legndr
   ! externals
   integer::inow,jnow,lnow,nl,max
   real(kr)::c(*)
   ! internals
   integer::nll,ie,ip,ir,nmu,nr,next,i,ll,np,idis,l,j,nmm,iflag,nw
   real(kr)::epnext,ep,epn,f,ff,sum,fact,f1,f2,tmid,test
   integer,parameter::mxlg=65
   real(kr)::term(mxlg),p(mxlg),amu(50),fmu(50)
   integer,parameter::nqp=8
   real(kr),dimension(8),parameter::qp=(/&
     -.9602898565e0_kr,-.7966664774e0_kr,-.5255324099e0_kr,&
     -.1834346425e0_kr,.1834346425e0_kr,.5255324099e0_kr,&
     .7966664774e0_kr,.9602898565e0_kr/)
   real(kr),dimension(8),parameter::qw=(/&
     .1012285362e0_kr,.2223810345e0_kr,.3137066459e0_kr,&
     .3626837834e0_kr,.3626837834e0_kr,.3137066459e0_kr,&
     .2223810345e0_kr,.1012285363e0_kr/)
   real(kr),parameter::tol=.005e0_kr
   real(kr),parameter::shade=.99999e0_kr
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::zero=0

   !--initialize output record
   jnow=lnow
   c(jnow)=0
   c(jnow+1)=c(inow+1)
   c(jnow+2)=0
   c(jnow+3)=nl-1
   nll=nl-1
   ie=0
   lnow=jnow+6

   !--loop over union grid of secondary energies
   epnext=0
   ip=2
   ir=1
   do while (epnext.lt.emax*(1-small))
      ep=epnext
      epnext=emax
      if (lnow+nl+1.gt.max)&
        call error('ll2lab','storage exceeded.',' ')

      !--retrieve the angular distribution
      nmu=nint(c(inow+5))
      nr=nint(c(inow+4))
      next=inow+6+2*nr
      do i=1,nmu
         ll=next
         amu(i)=c(ll+1)
         np=nint(c(ll+5))
         nr=nint(c(ll+4))
         next=ll+6+2*nr+2*np
         call terpa(fmu(i),ep,epn,idis,c(ll),ip,ir)
         test=shade*epn
         if (idis.gt.0.and.ep.lt.test) epn=test
         if (epn.lt.epnext*(1-small)) epnext=epn
      enddo

      !--compute legendre coefficients for this distribution.
      do l=1,nl
         term(l)=0
      enddo
      do i=1,nqp
         call legndr(qp(i),p,nll)
         j=1
         do while (qp(i).ge.amu(j))
            j=j+1
         enddo
         do l=1,nl
            f=(qp(i)-amu(j-1))/(amu(j)-amu(j-1))
            ff=(1-f)*fmu(j-1)+f*fmu(j)
            term(l)=term(l)+ff*p(l)*qw(i)
         enddo
      enddo

      !--check normalization
      nmm=nmu-1
      sum=0
      do i=1,nmm
         sum=sum+(amu(i+1)-amu(i))*(fmu(i+1)+fmu(i))/2
      enddo
      fact=1
      if (term(1).ne.zero) fact=sum/term(1)
      c(lnow)=ep
      do i=1,nl
         c(lnow+i)=term(i)*fact
      enddo

      !--is the previous point still needed?
      if (ie.lt.3) then
         ie=ie+1
         lnow=lnow+nl+1
      else
         f1=(c(lnow-nl-1)-c(lnow-2*nl-2))/(c(lnow)-c(lnow-2*nl-2))
         f2=1-f1
         iflag=0
         do i=1,nl
            tmid=f1*c(lnow-2*nl-2+i)+f2*c(lnow+i)
            test=tol*abs(tmid)+small/100
            if (abs(tmid-c(lnow-nl-1+i)).gt.test) iflag=1
         enddo
         if (iflag.eq.0) then
            c(lnow-nl-1)=c(lnow)
            do i=1,nl
               c(lnow-nl-1+i)=c(lnow+i)
            enddo
         else
            ie=ie+1
            lnow=lnow+nl+1
         endif
      endif
   enddo

   !--finished.
   nw=lnow-jnow-6
   c(jnow+4)=nw
   c(jnow+5)=ie
   return
   end subroutine ll2lab

   subroutine f6lab(ep,epnext,term,nl,law,int,lang,lep,clo,chi,e)
   !-------------------------------------------------------------------
   ! Retrieve the Legendre coefficients of the double-differential
   ! cross section at ep in the laboratory system from a part
   ! of File 6 in law 1 format.  Only the continuum part is
   ! returned by this subroutine.  Delta functions, if given,
   ! must be handled separately.  Call with ep=0 to initialize.
   ! Subsequent values of ep must be requested in increasing order.
   !-------------------------------------------------------------------
   use util  ! provides error
   use endf  ! provides terp1
   use mathm ! provides legndr
   ! externals
   integer::nl,law,int,lang,lep
   real(kr)::ep,epnext,term(nl),clo(*),chi(*),e
   ! internals
   integer::ilo,llo,mlo,nclo,ihi,lhi,mhi,nchi,intt,l,idone
   integer::ndlo,ndhi,nplo,nphi,int2,na,i,j,jj,nll
   real(kr)::xend,xlo,xhi,elo,ehi,epn,eplast,f1,f2,epp
   real(kr)::t0,t1,t2
   integer,parameter::mxlg=65
   real(kr)::term1(mxlg),term2(mxlg),p(mxlg)
   integer,parameter::nqp=8
   real(kr),dimension(8),parameter::qp=(/&
     -.9602898565e0_kr,-.7966664774e0_kr,-.5255324099e0_kr,&
     -.1834346425e0_kr,.1834346425e0_kr,.5255324099e0_kr,&
     .7966664774e0_kr,.9602898565e0_kr/)
   real(kr),dimension(8),parameter::qw=(/&
     .1012285362e0_kr,.2223810345e0_kr,.3137066459e0_kr,&
     .3626837834e0_kr,.3626837834e0_kr,.3137066459e0_kr,&
     .2223810345e0_kr,.1012285363e0_kr/)
   real(kr),parameter::up=1.00001e0_kr
   real(kr),parameter::dn=.99999e0_kr
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::emin=1.e-5_kr
   real(kr),parameter::zero=0
   save xend
   save ilo,llo,mlo,nclo,xlo
   save ihi,lhi,mhi,nchi,xhi
   save elo,ehi,intt

   !--initialization
   if (ep.eq.zero) then
      xlo=1
      xhi=1
      xend=1
      elo=clo(2)
      ehi=chi(2)
      ndlo=nint(clo(3))
      ndhi=nint(chi(3))
      intt=int
      if (intt.gt.20) intt=int-20
      if (intt.gt.10) intt=int-10
      nplo=nint(clo(6))
      nclo=nint(clo(5))/nplo
      mlo=7+(nplo-1)*nclo
      ilo=7+nclo*ndlo
      llo=ilo
      if (clo(llo).le.zero) llo=llo+nclo
      if (int.gt.20) xlo=clo(mlo)
      nphi=nint(chi(6))
      nchi=nint(chi(5))/nphi
      mhi=7+(nphi-1)*nchi
      ihi=7+nchi*ndhi
      lhi=ihi
      if (chi(lhi).le.zero) lhi=lhi+nchi
      if (int.gt.20) xhi=chi(mhi)
      if (int.lt.11.or.int.gt.15) then
         xend=xlo+(e-elo)*(xhi-xlo)/(ehi-elo)
         if (clo(ilo+1).gt.0.0.or.chi(ihi+1).gt.zero) then
            epnext=emin
         else
            epnext=clo(llo)/xlo
            epn=chi(lhi)/xhi
            if (epn.lt.epnext*(1-small)) epnext=epn
            epnext=epnext*xend
         endif
      else
         if (clo(ilo+1).gt.0.0.or.chi(ihi+1).gt.zero) then
            epnext=emin
         else
            call terp1(elo,clo(llo),ehi,chi(lhi),e,epnext,2)
         endif
      endif
      do l=1,nl
         term(l)=0
      enddo
      return
   endif

   !--normal entry
   !--law 1.  find energy range.
   if (int.lt.11.or.int.gt.15) then
      idone=0
      do while (idone.eq.0)
         if (llo.gt.mlo) then
            epnext=emax
            idone=1
         else
            epnext=clo(llo)*xend/xlo
            if (ep.lt.epnext*(1-small)) then
               idone=1
            else
               llo=llo+nclo
            endif
         endif
      enddo
      idone=0
      do while (idone.eq.0)
         if (lhi.gt.mhi) then
            epn=emax
            idone=1
         else
            epn=chi(lhi)*xend/xhi
            if (ep.lt.epn*(1-small)) then
               idone=1
            else
               lhi=lhi+nchi
            endif
         endif
      enddo
      if (epn.lt.epnext*(1-small)) epnext=epn
   else
      idone=0
      do while (idone.eq.0)
         call terp1(elo,clo(llo),ehi,chi(lhi),e,epnext,2)
         if (ep.lt.epnext*(1-small)) then
            idone=1
         else
            if (llo.ge.mlo) then
               epnext=emax
               idone=1
            else
               llo=llo+nclo
               lhi=lhi+nchi
            endif
         endif
      enddo
      if (llo.ne.ilo.and.llo.le.mlo) then
         call terp1(elo,clo(llo-nclo),ehi,chi(lhi-nchi),e,eplast,2)
         f1=(epnext-ep)/(epnext-eplast)
         f2=1-f1
      endif
   endif

   !--interpolate for coefficients.
   if (lang.eq.1) then
      do l=1,nl
         if (int.lt.11.or.int.gt.15) then
            term1(l)=0
            if (llo.ne.ilo.and.llo.le.mlo) then
               if (l.le.nclo-1) then
                  call terp1(clo(llo-nclo),clo(llo-nclo+l),&
                    clo(llo),clo(llo+l),ep*xlo/xend,term1(l),lep)
                  term1(l)=term1(l)*xlo
               endif
            endif
            term2(l)=0
            if (lhi.ne.ihi.and.lhi.le.mhi) then
               if (l.le.nchi-1) then
                  call terp1(chi(lhi-nchi),chi(lhi-nchi+l),&
                    chi(lhi),chi(lhi+l),ep*xhi/xend,term2(l),lep)
                  term2(l)=term2(l)*xhi
               endif
            endif
            call terp1(elo,term1(l),ehi,term2(l),&
              e,term(l),intt)
            term(l)=term(l)/xend
         else
            term1(l)=0
            if (llo.ne.ihi.and.llo.le.mlo) then
               if (l.le.nclo-1) then
                  epp=f1*clo(llo)+f2*clo(llo-nclo)
                  call terp1(clo(llo-nclo),clo(llo-nclo+l),&
                    clo(llo),clo(llo+l),epp,term1(l),lep)
               endif
            endif
            term2(l)=0
            if (lhi.ne.ihi.and.lhi.le.mhi) then
               if (l.le.nchi-1) then
                  epp=f1*chi(lhi)+f2*chi(lhi-nchi)
                  call terp1(chi(lhi-nchi),chi(lhi-nchi+l),&
                    chi(lhi),chi(lhi+l),epp,term2(l),lep)
               endif
            endif
            call terp1(elo,term1(l),ehi,term2(l),e,term(l),intt)
         endif
      enddo

   !--compute coefficients from tabulated distribution.
   else
      if (lang.lt.11.or.lang.gt.15)&
        call error('f6lab','illegal lang.',' ')
      int2=lang-10
      nll=nl-1
      na=(nclo-2)/2
      do l=1,nl
         term1(l)=0
      enddo
      do i=1,nqp
         call legndr(qp(i),p,nll)
         t0=0
         if (llo.ne.ilo.and.llo.le.mlo) then
            j=1
            jj=llo-nclo+2
            do while (clo(jj).le.qp(i).and.j.lt.na)
               j=j+1
               jj=jj+2
            enddo
            call terp1(clo(jj-2),clo(jj-1),clo(jj),clo(jj+1),&
              qp(i),t1,int2)
            j=1
            jj=llo+2
            do while (clo(jj).le.qp(i).and.j.lt.na)
               j=j+1
               jj=jj+2
            enddo
            call terp1(clo(jj-2),clo(jj-1),clo(jj),clo(jj+1),&
              qp(i),t2,int2)
            call terp1(clo(llo-nclo),t1,clo(llo),t2,&
              ep*xlo/xhi,t0,lep)
         endif
         t0=t0*qw(i)
         do l=1,nl
            term1(l)=term1(l)+t0*p(l)
         enddo
      enddo
      do l=1,nl
         term2(l)=0
      enddo
      do i=1,nqp
         call legndr(qp(i),p,nll)
         t0=0
         if (lhi.ne.ihi.and.lhi.le.mhi) then
            j=1
            jj=lhi-nchi+2
            do while (chi(jj).le.qp(i).and.j.lt.na)
               j=j+1
               jj=jj+2
            enddo
            call terp1(chi(jj-2),chi(jj-1),chi(jj),chi(jj+1),&
              qp(i),t1,int2)
            j=1
            jj=lhi+2
            do while (chi(jj).le.qp(i).and.j.lt.na)
               j=j+1
               jj=jj+2
            enddo
            call terp1(chi(jj-2),chi(jj-1),chi(jj),chi(jj+1),&
              qp(i),t2,int2)
            call terp1(chi(lhi-nchi),t1,chi(lhi),t2,&
              ep*xhi/xhi,t0,lep)
         endif
         t0=t0*qw(i)
         do l=1,nl
            term2(l)=term2(l)+t0*p(l)
         enddo
      enddo
      do l=1,nl
         call terp1(elo,term1(l),ehi,term2(l),e,term(l),intt)
      enddo
   endif

   !--finished with law 1.
   if (lep.gt.1) return
   if (abs(epnext-emax).lt.small*emax) return
   if (ep.ge.dn*dn*epnext) then
      epnext=up*epnext
   else
      epnext=dn*epnext
   endif
   return
   end subroutine f6lab

   subroutine getdis(e,enext,idisc,yld,ff,nl,ng,iglo,nq,nin,nlg)
   !-------------------------------------------------------------------
   ! Computes feed function for elastic or discrete inelastic
   ! scattering using data from either File 4 or File 6.
   !-------------------------------------------------------------------
   use util    ! provides error
   use mathm   ! provides legndr
   use physics ! provides amassn,amu,hbar,ev,clight
   ! externals
   integer::idisc,nl,ng,iglo,nq,nin,nlg
   real(kr)::e,enext,yld,ff(nlg,*)
   ! internals
   integer::iecl,iech,nqp0,nld,lidp,lcd,mft,mtt,idis
   integer::npo,ig,il,nqp,iqp,i2s,l,ndig,iii,ii,ld
   integer::ignow,ngnow,ngn1
   real(kr)::ecl,ech,ecn,en,awr2,ast,whi,ep,aa,b,prob
   real(kr)::ai,at,cc1,ee,cc2,eta,sigc,wlab,fact,test
   real(kr)::ef,disc,af,wmin,wmax,wlo,wqp,wqw,zt,zi,wn
   real(kr)::fle(21),p(21),flt(20)
   real(kr),dimension(4),parameter::qp4=(/&
     -.8611363116e0_kr,-.3399810436e0_kr,.3399810436e0_kr,&
    .8611363116e0_kr/)
   real(kr),dimension(4),parameter::qw4=(/&
     .3478548451e0_kr,.6521451549e0_kr,.6521451549e0_kr,&
     .3478548451e0_kr/)
   real(kr),dimension(8),parameter::qp8=(/&
     -.9602898565e0_kr,-.7966664774e0_kr,-.5255324099e0_kr,&
     -.1834346425e0_kr,.1834346425e0_kr,.5255324099e0_kr,&
     .7966664774e0_kr,.9602898565e0_kr/)
   real(kr),dimension(8),parameter::qw8=(/&
     .1012285362e0_kr,.2223810345e0_kr,.3137066459e0_kr,&
     .3626837834e0_kr,.3626837834e0_kr,.3137066459e0_kr,&
     .2223810345e0_kr,.1012285363e0_kr/)
   real(kr),dimension(12),parameter::qp12=(/&
     -.9815606342e0_kr,-.9041172564e0_kr,-.7699026742e0_kr,&
     -.5873179543e0_kr,-.3678314990e0_kr,-.1252334085e0_kr,&
     .1252334085e0_kr,.3678314990e0_kr,.5873179543e0_kr,&
     .7699026742e0_kr,.9041172564e0_kr,.9815606342e0_kr/)
   real(kr),dimension(12),parameter::qw12=(/&
     .0471753364e0_kr,.1069393260e0_kr,.1600783285e0_kr,&
     .2031674267e0_kr,.2334925365e0_kr,.2491470458e0_kr,&
     .2491470458e0_kr,.2334925365e0_kr,.2031674267e0_kr,&
     .1600783285e0_kr,.1069393260e0_kr,.0471753364e0_kr/)
   real(kr),dimension(20),parameter::qp20=(/&
     -.9931286e0_kr,-.9639719e0_kr,-.9122344e0_kr,&
     -.839117e0_kr,-.7463319e0_kr,-.6360537e0_kr,&
     -.510867e0_kr,-.3737061e0_kr,-.2277859e0_kr,&
     -.0765265e0_kr,.0765265e0_kr,.2277859e0_kr,&
     .3737061e0_kr,.5108670e0_kr,.6360537e0_kr,&
     .7463319e0_kr,.8391170e0_kr,.9122344e0_kr,&
     .9639719e0_kr,.9931286e0_kr/)
   real(kr),dimension(20),parameter::qw20=(/&
     .0176140e0_kr,.0406014e0_kr,.0626720e0_kr,&
     .0832767e0_kr,.1019301e0_kr,.1181945e0_kr,&
     .1316886e0_kr,.1420961e0_kr,.1491730e0_kr,&
     .1527534e0_kr,.1527534e0_kr,.1491730e0_kr,&
     .1420961e0_kr,.1316886e0_kr,.1181945e0_kr,&
     .1019301e0_kr,.0832767e0_kr,.0626720e0_kr,&
     .0406014e0_kr,.0176140e0_kr/)
   real(kr),parameter::fm=1.e-12_kr
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::wcut=.94e0_kr
   real(kr),parameter::shade=.99999995e0_kr
   real(kr),parameter::zero=0
   real(kr),parameter::ten=10
   save iecl,iech,ecl,ech,ecn,nqp0

   nld=21
   ngn1=ngn+1
   enext=emax
   idisc=0
   lidp=0
   if (izap.eq.izat) lidp=1
   yld=1
   if (lrflag.eq.16.or.lrflag.eq.21.or.lrflag.eq.24&
     .or.lrflag.eq.26.or.lrflag.eq.30) yld=2
   if (lrflag.eq.17.or.lrflag.eq.25.or.lrflag.eq.38) yld=3
   if (lrflag.eq.37) yld=4
   lcd=2
   mft=4
   if (mfd.eq.8) mft=8
   if (mfd.ge.21.and.mfd.le.26) mft=8
   if (mfd.eq.3.and.mtd.eq.2.and.izap.gt.1) mft=8
   mtt=mtd
   if (mtd.eq.251.or.mtd.eq.252.or.mtd.eq.253) mtt=2
   call getfle(e,en,idis,fle,nld,lcd,matd,mft,mtt,nin)
   if (en.lt.enext*(1-small)) idisc=idis
   if (en.lt.enext*(1-small)) enext=en
   awr2=awr*(awr+1-aprime)/aprime
   if (e.eq.0) go to 600
   ld=nld-1
   npo=ld+nl+int(log(300/awr))
   if (mfd.eq.3) npo=npo+lord
   if (izap.gt.1) npo=npo+6
   nqp0=4
   if (npo.gt.8) nqp0=8
   if (npo.gt.12) nqp0=12
   if (npo.gt.16) nqp0=20
   ast=0
   if (e.gt.thresh) ast=sqrt(awr2)*sqrt(1-thresh/e)
   ng=0
   ig=0
   ngnow=1
   ignow=0
   iglo=iech-1
   if (iglo.lt.1) iglo=1
   wmin=-1
   wmax=1
   if (izap.gt.1.and.mtd.eq.2) then
      wmax=wcut
      if (lidp.eq.1) wmin=-wcut
   endif
   whi=wmin
  410 continue
   ng=ng+1
   if (ng.gt.ngn1) then
      ng=ngn1
      ig=ignow
      do il=1,nl
         ff(il,ng)=0
      enddo
      go to 465
   endif
   ngnow=ng
   do il=1,nl
      ff(il,ng)=0
   enddo
   if (ast.eq.0) go to 465
  425 continue
   wlo=whi
   ig=ig+1
   if (iglo+ig.gt.ngn1) then
      ng=ngnow
      ig=ig-1
      go to 465
   endif
   ignow=ig
   ep=egn(iglo+ig)
   if ((mtd.eq.252.or.mtd.eq.253).and.ep.ge.e) go to 465
   whi=((ep/e)*(1+awr)**2/aprime-(1+ast*ast))/(2*ast)
   if (whi.lt.shade*wmin) whi=wmin
   if (whi.gt.shade*wmax) whi=wmax
   if ((mtd.eq.251.or.mtd.eq.252.or.mtd.eq.253)&
     .and.wlo.eq.whi) go to 425
   if (mfd.eq.3.and.wlo.eq.whi) go to 425
   if (wlo.eq.whi) go to 410
   if ((mtd.eq.252.or.mtd.eq.253).and.&
     whi.ge.wmax) go to 464

   !--integrate between wlo and whi using gaussian quadrature
   aa=(whi+wlo)/2
   b=(whi-wlo)/2
   nqp=int(npo*2*b)
   if (nqp.gt.npo) nqp=npo
   if (nqp.le.8) nqp=4
   if (nqp.gt.8.and.nqp.le.12) nqp=8
   if (nqp.gt.12.and.nqp.le.16) nqp=12
   if (nqp.gt.16) nqp=20
   do il=1,nl
      flt(il)=0
   enddo
   do iqp=1,nqp
      if (nqp.eq.4) then
         wqp=aa+b*qp4(iqp)
         wqw=b*qw4(iqp)
      else if (nqp.eq.8) then
         wqp=aa+b*qp8(iqp)
         wqw=b*qw8(iqp)
      else if (nqp.eq.12) then
         wqp=aa+b*qp12(iqp)
         wqw=b*qw12(iqp)
      else if (nqp.eq.20) then
         wqp=aa+b*qp20(iqp)
         wqw=b*qw20(iqp)
      else
         call error('getdis','illegal nqp',' ')
      endif
      if (law.eq.4) wqp=-wqp

      !--compute scattering probability at this cm cosine
      call legndr(wqp,p,ld)
      prob=0
      do il=1,nld
         prob=prob+fle(il)*p(il)*(2*il-1)/2
      enddo

      !--add in coulomb part for charged particles
      if (mtd.eq.2.and.izap.gt.1) then
         i2s=nint(2*spi)
         ai=awrp*amassn
         at=awr*amassn
         zt=int(izat/1000)
         zi=int(izap/1000)
         cc1=2*amu*ev*fm**2/hbar**2
         ee=(ev/10000000)*(clight/10)
         cc2=ee**4*amu/(2*hbar**2*ev)
         wn=at*sqrt(cc1*e*ai)/(ai+at)
         eta=zt*zi*sqrt(cc2*ai/e)
         sigc=0
         if (lidp.eq.0) sigc=(eta**2/wn**2)/(1-wqp)**2
         if (lidp.eq.1) sigc=((2*eta**2/wn**2)&
           /(1-wqp**2))*((1+wqp**2)/(1-wqp**2)&
           +(-1**i2s)*cos(eta*log((1+wqp)/(1-wqp)))/(2*spi+1))
         if (lidp.eq.0) prob=prob/(1-wqp)
         if (lidp.eq.1) prob=prob/(1-wqp*wqp)
         prob=prob+sigc
      endif

      !--multiply scattering probability by legendre polynomials at
      !--corresponding lab cosine and accumulate integrals
      wlab=(1+ast*wqp)/sqrt(1+ast*ast+2*ast*wqp)
      l=nl-1
      call legndr(wlab,p,l)
      prob=yld*wqw*prob
      if (mtd.eq.251) prob=wlab*prob
      do il=1,nl
         flt(il)=flt(il)+prob*p(il)
      enddo

   !--end of quadrature loop
   enddo
   do il=1,nl
      ff(il,ng)=ff(il,ng)+flt(il)
   enddo
   if (whi.ge.wmax) go to 465
   if (iglo+ig.gt.ngn) go to 465
   if (mtd.eq.251.or.mtd.eq.252.or.mtd.eq.253) go to 425
   if (mfd.eq.3) go to 425
   go to 410
  464 continue
   ff(1,1)=1
   ig=ig+1
  465 continue

   !--reduce significant figures for small values
   ndig=7
   fact=ten**ndig
   do il=1,nl
      do ii=1,ng
         iii=nint(fact*ff(il,ii)+ten**(ndig-11))
         ff(il,ii)=iii/fact
      enddo
   enddo
   if (mfd.eq.3.and.mtd.eq.2) iglo=1
   if (mtd.eq.251.or.mtd.eq.252.or.mtd.eq.253) then
      do il=1,nl
         if (mtd.eq.251.or.mtd.eq.252) ff(il,ng+1)=1
         if (iglo+ig-1.ge.1) then
            if (mtd.eq.252) ff(il,ng)=ff(il,ng)&
              *log(egn(iglo+ig)/egn(iglo+ig-1))
            if (mtd.eq.253) then
               ff(il,ng+1)=ff(il,ng)
               ff(il,ng)=ff(il,ng)*log(e/egn(iglo+ig-1))
            endif
         endif
      enddo
      ng=ng+1
      iglo=1
   endif
   go to 610

   !--check next critical point for
   !--discrete channel scattering
  600 continue
   if (e.gt.zero) go to 610
   iecl=0
   iech=0
   ecl=0
   ech=0
   ecn=0
  610 continue
   if (e.lt.shade*ecn) go to 680
  620 continue
   if (e.lt.shade*ecl) go to 640
   ecl=emax
   if (iecl.eq.ngn+1) go to 640
   iecl=iecl+1
   if (thresh.gt.zero) go to 630
   test=1
   if (aprime.gt.test) go to 625
   ecl=egn(iecl)
   go to 620
  625 continue
   ecl=egn(iecl)*(awr+1)**2/(4*awr)
   go to 620
  630 continue
   ef=(awr+1)*egn(iecl)/(awr+1-aprime)
   aa=1+ef/(-q)
   disc=(awr2*aa-1)*ef/(-q)
   if (disc.lt.zero) go to 620
   disc=sqrt(disc)
   af=(1-disc)/aa
   if (af**2.eq.awr2) then
      ecl=emax
   else
      ecl=thresh/(1-af**2/awr2)
   endif
   go to 620
  640 continue
   if (e.lt.shade*ech) go to 670
   ech=emax
   if (iech.eq.ngn+1) go to 670
   iech=iech+1
   if (thresh.gt.zero) go to 650
   test=1
   if (aprime.gt.test) go to 645
   test=test/10
   if (alpha.lt.test) go to 640
   ech=egn(iech)/alpha
   go to 640
  645 continue
   ech=emax
   go to 640
  650 continue
   ef=(awr+1)*egn(iech)/(awr+1-aprime)
   aa=1+ef/(-q)
   disc=(awr2*aa-1)*ef/(-q)
   if (disc.lt.zero) go to 640
   disc=sqrt(disc)
   af=(1+disc)/aa
   if (af**2.eq.awr2) then
      ech=emax
   else
      ech=thresh/(1-af**2/awr2)
   endif
   go to 640
  670 continue
   ecn=ecl
   if (ech.lt.ecn*(1-small)) ecn=ech
  680 continue
   if (ecn.lt.enext*(1-small)) idisc=1
   if (ecn.lt.enext*(1-small)) enext=ecn
   nq=nqp0
   if (iecl.eq.iech) nq=0
   return
   end subroutine getdis

   subroutine getfle(e,enext,idis,fle,nle,lcd,matd,mfd,mtd,nin)
   !-------------------------------------------------------------------
   ! Retrieve or compute Legendre coefficients at e.
   ! Return isotropic distribution below first point.
   !-------------------------------------------------------------------
   use util ! provides error,mess
   use endf ! provides endf routines and variables
   ! externals
   integer::idis,nle,lcd,matd,mfd,mtd,nin
   real(kr)::e,enext,fle(*)
   ! internals
   integer::iso,lidp,ltt,iraw,ir,nne,ne,int,nlo,nhi,ltt3,lttn
   integer::nb,nwc,lvt,i,nlmax,il
   real(kr)::elo,ehi
   character(60)::strng
   integer,parameter::mxlg=65
   real(kr)::flo(mxlg),fhi(mxlg)
   integer,parameter::ncmax=350
   real(kr)::fls(ncmax)
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::over=1.01e0_kr
   real(kr),parameter::zero=0
   save iso,lidp,ltt,iraw,ir,nne,ne,int
   save elo,ehi,nlo,nhi,flo,fhi
   save ltt3,lttn,fls

   !--initialize
   idis=0
   if (e.gt.zero) go to 200
   iso=0
   if (mfd.ne.8) then
      call findf(matd,mfd,mtd,nin)
      call contio(nin,0,0,fls,nb,nwc)
      awr=fls(2)
      lvt=nint(fls(3))
      ltt=nint(fls(4))
      ltt3=ltt
      if (ltt.eq.3) then
         ltt=1
         lttn=1
      endif
      if (lvt.eq.0) then
         call contio(nin,0,0,fls,nb,nwc)
      else
         call listio(nin,0,0,fls,nb,nwc)
         do while (nb.ne.0)
            call moreio(nin,0,0,fls(7),nb,nwc)
         enddo
      endif
      iso=nint(fls(3))
      lct=nint(fls(4))
   else
      if (iabs(law).eq.3) iso=1
      if (law.eq.-4) iso=1
      ltt=1
   endif
   ! check for lab distribution for two-body scattering.
   if (lct.eq.1.and.awr.ge.10) then
      if (mtd.ge.51.and.mtd.le.90) then
         write(strng,&
           '(''lab distribution changed to cm for mt='',i3)') mtd
         call mess('getfle',strng,' ')
         lct=2
      endif
   endif
   elo=emax
   if (iso.ne.1) then
      call tab2io(nin,0,0,fls,nb,nwc)
      lidp=-1
      if (iabs(law).eq.5) lidp=nint(fls(3))
      ne=nint(fls(6))
      ! read in raw data at first two incident energies.
      ! retrieve or compute legendre coefficients.
      iraw=1+nwc
      if (ltt.eq.1) call listio(nin,0,0,fls(iraw),nb,nwc)
      if (ltt.eq.2) call tab1io(nin,0,0,fls(iraw),nb,nwc)
      elo=fls(iraw+1)
      nlo=nle
      if (lidp.ge.0) fls(iraw+3)=lidp
      call getco(flo,nlo,lcd,fls(iraw),lct,ltt,idis)
      if (ltt.eq.1) call listio(nin,0,0,fls(iraw),nb,nwc)
      if (ltt.eq.2) call tab1io(nin,0,0,fls(iraw),nb,nwc)
      ehi=fls(iraw+1)
      nhi=nle
      if (lidp.ge.0) fls(iraw+3)=lidp
      call getco(fhi,nhi,lcd,fls(iraw),lct,ltt,idis)
      ! check for total isotropy.
      if (ne.eq.2.and.nlo.eq.1.and.nhi.eq.1) iso=1
      int=nint(fls(8))
   endif
   if (iso.eq.1) nle=1
   enext=elo
   nne=2
   ir=1
   return

   !--normal entry
   !--is desired energy in current panel
  200 continue
   if (nle.eq.1) go to 440
   if (iso.eq.1) go to 440
   if (e.lt.ehi*(1-small)) go to 300
   if (nne.eq.ne.and.e.lt.over*ehi) then
      if (ltt3.eq.3.and.lttn.eq.1) go to 210
      go to 300
   endif

   !--no.  slide raw high energy data into low energy positions.
  210 continue
   do i=1,nhi
      flo(i)=fhi(i)
   enddo
   nlo=nhi
   elo=ehi

   !--read in new raw data at high energy.
   if (nne.eq.ne.and.ltt3.eq.3.and.lttn.eq.1) then
      call tab2io(nin,0,0,fls(1),nb,nwc)
      ne=nint(fls(6))
      nne=0
      ir=1
      ltt=2
      lttn=2
   else if (nne.eq.ne) then
      call error('getfle','desired energy above highest given.',' ')
   endif
   if (ltt.eq.1) call listio(nin,0,0,fls(iraw),nb,nwc)
   if (ltt.eq.2) call tab1io(nin,0,0,fls(iraw),nb,nwc)
   ehi=fls(iraw+1)
   nhi=nle
   if (lidp.ge.0) fls(iraw+3)=lidp
   call getco(fhi,nhi,lcd,fls(iraw),lct,ltt,idis)
   nne=nne+1
   if (nne.gt.fls(5+2*ir)) ir=ir+1
   int=nint(fls(6+2*ir))
   if (ehi.le.e*(1-small).and.nne.lt.ne) go to 210

   !--yes.  interpolate for coefficients at desired energy.
  300 continue
   if (e.ge.elo*(1-small)) then
      ! normal energies
      nlmax=nlo
      if (nlmax.lt.nhi) nlmax=nhi
      do i=1,nle
         if (i.le.nlmax) then
            call terp1(elo,flo(i),ehi,fhi(i),e,fle(i),int)
         else
            fle(i)=0
         endif
      enddo
      nle=nlmax
      enext=ehi
      if (int.eq.1) idis=1
   else
      !--return isotropic distribution below first point.
      fle(1)=1
      do il=2,nle
         fle(il)=0
      enddo
      nle=1
      enext=ehi
      idis=1
   endif
   return

   !--isotropic distribution or request.
  440 continue
   fle(1)=1
   if (nle.gt.1) then
      do il=2,nle
         fle(il)=0
      enddo
      nle=1
   endif
   enext=emax
   idis=0
   return
   end subroutine getfle

   subroutine getaed(aed,e,enext,idis,nl,nq,mat,mf,mt,nin,nlg)
   !-------------------------------------------------------------------
   ! Retrieve or compute Legendre coefficients at e for thermal
   ! angle and energy distributions in File 6.  A special incident
   ! energy interpolation scheme is used for incoherent inelastic
   ! data.  Angular data for coherent elastic scattering is
   ! reconstructed using the Bragg edges in the cross section.
   !-------------------------------------------------------------------
   use util  ! provides error
   use endf  ! provides endf routines and variables
   use mathm ! provides legndr
   ! externals
   integer::idis,nl,nq,mat,mf,mt,nin,nlg
   real(kr)::aed(nlg,*),e,enext
   ! internals
   integer::itt,law,ir,ne,nwt,ncyc,nu,nlo,ie,nw1,nw2,l1,idone
   integer::l2,l3,llo,lhi,nhi,jbrag,nbrag,nb,nw
   integer::i,il,ig,k1,k2,iu,ib,loc
   real(kr)::elo,ehi,clast,ebrag,xs,cnow,eg,eb,egp,egp1,egp2
   real(kr)::eg1,eg2,ek1,ei1,ek2,ei2,ei,ee,f1,f2,x1,x2
   real(kr)::u1,u2,u
   real(kr)::p(20)
   real(kr)::fl1(20),fl2(20),fi(20),fl(20)
   real(kr),parameter::emin=1.e-5_kr
   real(kr),parameter::emax=10.e0_kr
   real(kr),parameter::etop=1.e10_kr
   real(kr),parameter::eps=1.e-5_kr
   real(kr),parameter::shade=.999999e0_kr
   real(kr),parameter::step=1.1e0_kr
   real(kr),parameter::up=1.00001e0_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::zero=0
   integer,parameter::maxaes=200000
   save nwt,ir,ncyc,nu,ne,ie,nw1,nw2
   save l2,l3,llo,nlo,elo,lhi,nhi,ehi
   save itt,law,jbrag,nbrag

   !--initialize.
   if (e.eq.zero) then
      call findf(mat,mf,mt,nin)
      call contio(nin,0,0,p,nb,nw)
      call tab1io(nin,0,0,p,nb,nw)
      itt=-l1h
      law=l2h

      !--incoherent elastic or inelastic with E-E'-mu ordering
      if (law.eq.1) then
         allocate(aes(maxaes))
         l1=1
         call tab2io(nin,0,0,aes(l1),nb,nw)
         ir=1
         ne=nint(aes(l1+5))
         nwt=nw
         l2=l1+nwt
         call listio(nin,0,0,aes(l2),nb,nw)
         elo=aes(l2+1)
         ncyc=nint(aes(l2+5))
         nu=ncyc-2
         nlo=nint(aes(l2+4))/ncyc
         loc=l2+nw
         do while (nb.ne.0)
            call moreio(nin,0,0,aes(loc),nb,nw)
            loc=loc+nw
            if (loc.gt.maxaes)&
              call error('getaed','storage exceeded.',' ')
         enddo
         llo=loc
         if (loc.gt.maxaes)&
           call error('getaed','storage exceeded.',' ')
         lhi=llo
         nw1=loc-l2
         ehi=0
         enext=elo
         ie=1

      !--incoherent inelastic with E-mu-E' ordering
      else if (law.eq.7) then
         call error('getaed','thermal mf6/law7 not coded',' ')

      !--coherent elastic
      else
         nbrag=itt
         allocate(aes(2*nbrag))
         enext=emin
         clast=0
         jbrag=0

         !--store bragg edges
         idone=0
         do while (idone.eq.0)
            ebrag=enext
            call gety1(ebrag,enext,idis,xs,npend,sigma)
            cnow=ebrag*xs
            if (enext.gt.emax*(1+small)) then
               idone=1
            else
               if (abs(enext-emax).le.eps) then
                  idone=1
               else
                  if ((cnow-clast).gt.(eps*clast+eps)) then
                     jbrag=jbrag+1
                     aes(jbrag)=ebrag
                     aes(nbrag+jbrag)=cnow-clast
                     clast=cnow
                  endif
               endif
            endif
         enddo

         !--reposition npend for getsig
         call findf(mat,3,mt,npend)
         call contio(npend,0,0,sigma,nb,nw)
         call gety1(e,enext,idis,xs,npend,sigma)
         enext=etop
      endif
      return
   endif

   !--normal entry for incoherent inelastic
   if (itt.ne.1) go to 490
   if (e.lt.ehi*(1-small)) go to 300
   if (ie.eq.ne.and.e.gt.ehi) go to 410
   if (ie.eq.ne) go to 300
   l1=1
   l2=l1+nwt
   l3=l2+nw1
   if (ie.eq.1) go to 220

   !--slide high data into low positions.
  205 continue
   do i=1,nw2
      aes(i-1+l2)=aes(i-1+l3)
   enddo
   llo=l2+lhi-l3
   nw1=nw2
   l3=l2+nw1
   elo=ehi
   nlo=nhi

   !--read new high data.
  220 continue
   call listio(nin,0,0,aes(l3),nb,nw)
   ehi=aes(l3+1)
   nhi=nint(aes(l3+4))/ncyc
   loc=l3+nw
   do while (nb.ne.0)
      call moreio(nin,0,0,aes(loc),nb,nw)
      loc=loc+nw
      if (loc.gt.maxaes)&
        call error('getaed','storage exceeded.',' ')
   enddo
   if (loc.gt.maxaes)&
     call error('getaed','storage exceeded.',' ')
   nw2=loc-l3
   ie=ie+1
   if (e.ge.ehi*(1-small).and.ie.lt.ne) go to 205
   if (ie.le.aes(l1+4+2*ir)) go to 300
   ir=ir+1

   !--begin loop over energy groups.
  300 continue
   do il=1,nl
      do ig=1,ngn
         aed(il,ig)=0
      enddo
   enddo
   k1=1
   k2=1
   do il=1,nl
      fl(il)=0
   enddo
   eg=0
   eb=elo/10
   i=0
   egp=0
   egp1=0
   egp2=0
  330 continue
   i=i+1
   egp=egn(i+1)
   if (i.eq.ngn) egp=egp+1
   eg1=egp1
   eg2=egp2
   egp1=egp-e+elo
   egp2=egp+ehi-e
   if (egp1.lt.eb*(1-small)) then
      egp1=egp*eb/(e-elo+eb)
      egp2=egp*(ehi-elo+eb)/(e-elo+eb)
   endif

   !--get next point projected from low side
  345 continue
   if (k1.gt.nlo) go to 360
   ek1=aes(l2+6+ncyc*(k1-1))
   if (ek1.gt.up*eg1) go to 355
   k1=k1+1
   go to 345
  360 continue
   ek1=aes(l3+6+ncyc*(nhi-1))
  355 continue
   ei1=ek1+e-elo
   if (ek1.lt.elo/10) ei1=ek1*(e-elo+eb)/eb

   !--get next point projected from high side
  375 continue
   if (k2.gt.nhi) go to 450
   ek2=aes(l3+6+ncyc*(k2-1))
   if (ek2.gt.up*eg2) go to 385
   k2=k2+1
   go to 375
  385 continue
   ei2=ek2-ehi+e
   if (ek2.lt.ehi-elo+eb)&
     ei2=ek2*(e-elo+eb)/(ehi-elo+eb)

   !--do integrals to next point
   ei=egp
   if (ei1.lt.ei) ei=ei1
   if (ei2.lt.ei) ei=ei2
   if (abs(ei-egp).lt.egp*small) then
      call aedi(egp1,fl1,nl,aes(l2))
      call aedi(egp2,fl2,nl,aes(l3))
   else if (abs(ei-ei1).lt.ei1*small) then
      call aedi(ek1,fl1,nl,aes(l2))
      ee=ek1+ehi-elo
      if (ek1.lt.eb*(1-small)) ee=ek1*(ehi-elo+eb)/eb
      call aedi(ee,fl2,nl,aes(l3))
      eg1=ek1
      eg2=ee
   else if (abs(ei-ei2).lt.ei2*small) then
      ee=ek2-ehi+elo
      if (ee.lt.eb*(1-small)) ee=ek2*eb/(ehi-elo+eb)
      call aedi(ee,fl1,nl,aes(l2))
      call aedi(ek2,fl2,nl,aes(l3))
      eg2=ek2
      eg1=ee
   endif
   f1=(ehi-e)/(ehi-elo)
   f2=(e-elo)/(ehi-elo)
   do il=1,nl
      fi(il)=f1*fl1(il)+f2*fl2(il)
   enddo
   ! write(6,'(1p,3e12.4)') ei,fi(1),fi(2)
   do il=1,nl
      aed(il,i)=aed(il,i)+(fi(il)+fl(il))*(ei-eg)/2
      fl(il)=fi(il)
   enddo
   eg=ei
   if (ei.lt.egp*(1-small)) go to 345

   !--close loop over energy groups
  450 continue
   if (i.lt.ngn.and.(k1.lt.nlo.or.k2.lt.nhi)) go to 330
   enext=ehi
   if (enext.le.e*(1+small)) enext=etop
   nq=0
   return

   !--upper end of thermal table
  410 continue
   do il=1,nl
      do ig=1,ngn
         aed(il,ig)=0
      enddo
   enddo
   enext=etop
   nq=0
   return

   !--normal entry for elastic
   !--find group index of in-group
  490 continue
   ig=1
   do while (e.ge.egn(ig+1)*(1-small).and.ig.ne.ngn)
      ig=ig+1
   enddo
   do il=1,nl
      do i=1,ngn
         aed(il,i)=0
      enddo
   enddo
   if (itt.gt.2) go to 600

   !--incoherent elastic
   if (e.lt.ehi*(1-small)) go to 560
   if (ie.eq.ne.and.e.gt.step*ehi) go to 410
   if (ie.eq.ne) go to 560
   if (ie.eq.1) go to 530

   !--slide high data into low positions
   do i=1,nw2
      aes(i-1+l2)=aes(i-1+l3)
   enddo
   nw1=nw2
   l3=l2+nw1
   elo=ehi
   llo=lhi

   !--read new high data
  530 continue
   ehi=0
   do while (e.ge.ehi*(1-small).and.ie.lt.ne)
      call listio(nin,0,0,aes(l3),nb,nw)
      ehi=aes(l3+1)
      nhi=1
      loc=l3+nw
      do while (nb.ne.0)
         call moreio(nin,0,0,aes(loc),nb,nw)
         loc=loc+nw
      enddo
      nw2=loc-l3
      ie=ie+1
   enddo

   !--compute legendre components
  560 continue
   x1=aes(l2+7)
   x2=aes(l3+7)
   xs=x1+(x2-x1)*(e-elo)/(ehi-elo)
   xs=xs/nu
   do iu=1,nu
      u1=aes(l2+7+iu)
      u2=aes(l3+7+iu)
      u=u1+(u2-u1)*(e-elo)/(ehi-elo)
      call legndr(u,p,nl)
      do il=1,nl
         aed(il,ig)=aed(il,ig)+xs*p(il)
      enddo
   enddo
   enext=ehi
   if (enext.gt.egn(ig+1)*(1+small)) enext=egn(ig+1)
   idis=1
   if (abs(enext-ehi).lt.ehi*small) idis=0
   if (enext.le.e*(1+small)) enext=etop
   nq=0
   return

   !--coherent elastic
  600 continue
   if (jbrag.ne.0) then
      ib=1
      do while (e.ge.shade*aes(ib).and.ib.le.jbrag)
         u=1-2*aes(ib)/e
         call legndr(u,p,nl)
         do il=1,nl
            aed(il,ig)=aed(il,ig)+aes(nbrag+ib)*p(il)/e
         enddo
         ib=ib+1
      enddo
      if (aed(1,ig).ne.zero) then
         if (nl.ge.2) then
            do il=2,nl
               aed(il,ig)=aed(il,ig)/aed(1,ig)
            enddo
         endif
      endif
   endif
   aed(1,ig)=1
   enext=egn(ig+1)
   idis=1
   nq=0
   return
   end subroutine getaed

   subroutine aedi(ee,fl,nl,aa)
   !-------------------------------------------------------------------
   ! Interpolate for Legendre components of thermal scattering at ee.
   !-------------------------------------------------------------------
   use mathm ! provides legndr
   ! externals
   integer::nl
   real(kr)::ee,fl(nl),aa(*)
   ! internals
   integer::ncyc,nw,np,nu,il,ip,i,iu
   real(kr)::f1,f2,u
   real(kr)::p(20)

   ncyc=nint(aa(6))
   nw=nint(aa(5))
   np=nw/ncyc
   nu=ncyc-2
   do il=1,nl
      fl(il)=0
   enddo
   ip=0
   i=1
   do while (ip.eq.0.and.i.lt.np)
      if (ee.le.aa(7+ncyc*i)) ip=i
      i=i+1
   enddo
   if (ip.gt.0) then
      f1=(aa(7+ncyc*ip)-ee)/(aa(7+ncyc*ip)-aa(7+ncyc*(ip-1)))
      f2=(ee-aa(7+ncyc*(ip-1)))/(aa(7+ncyc*ip)-aa(7+ncyc*(ip-1)))
      do iu=1,nu
         u=aa(8+iu+ncyc*(ip-1))
         call legndr(u,p,nl)
         do il=1,nl
            fl(il)=fl(il)+f1*aa(8+ncyc*(ip-1))*p(il)/nu
         enddo
         u=aa(8+iu+ncyc*ip)
         call legndr(u,p,nl)
         do il=1,nl
            fl(il)=fl(il)+f2*aa(8+ncyc*ip)*p(il)/nu
         enddo
      enddo
   endif
   return
   end subroutine aedi

   subroutine getgfl(ed,enext,idis,gfl,nl,nlg,ng,mat,mf,mt,nin)
   !-------------------------------------------------------------------
   ! Retrieve gamma angular distributions as Legendre coefficients.
   ! Routine returns coefficients for all gamma rays (discrete and
   ! continuum) simultaneously.  Coded for coefficient data only.
   ! Initialize if ed=0.
   !-------------------------------------------------------------------
   use endf ! proves endf routines and variables
   use util ! provides error
   ! externals
   integer::idis,nl,nlg,ng,mat,mf,mt,nin
   real(kr)::ed,enext,gfl(nlg,ng)
   ! internals
   integer::li,ltt,lcd,ig,il,j,i,im,m
   integer::ip1,k,l,lint,lnow,in,int,nb,nw,ni,np,nnp
   integer::nh,nipj,lnext,nbt,nld,loc,ntmp,na
   real(kr)::elo,ehi
   integer,parameter::mxlg=65
   real(kr)::b(6),alo(mxlg),ahi(mxlg)
   integer,parameter::maxgfl=500
   integer::loca(maxgfl)
   real(kr),dimension(:),allocatable::tmp
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::up=1.1e0_kr
   real(kr),parameter::zero=0
   save li,loca,ltt

   if (ed.gt.zero) go to 200

   !--initialize
   lct=1
   lcd=1
   call findf(mat,mf,mt,nin)
   call contio(nin,0,0,b,nb,nw)
   li=nint(b(3))
   ltt=nint(b(4))
   ni=nint(b(6))
   ntmp=10000
   allocate(tmp(ntmp))
   enext=emax

   !--reduce legendre order if possible
   if (li.eq.1) nl=1
   na=1
   if (li.eq.1) then
      allocate(flgn(1))
      deallocate(tmp)
      return
   endif

   !--find starting location of each sub-section
   ig=0
   il=1
   if (ni.gt.0) then
      do while (ig.lt.ni)
         call contio(nin,0,0,tmp(il),nb,nw)
         ig=ig+1
         if (ig.gt.maxgfl) call error('getgfl','too many gammas',' ')
         loca(ig)=il-1
         il=il+nw
      enddo
   endif
   do while (ig.lt.ng)
      call tab2io(nin,0,0,tmp(il),nb,nw)
      ig=ig+1
      if (ig.gt.maxgfl) call error('getgfl','too many gammas',' ')
      loca(ig)=il-1
      np=nint(tmp(il+5))
      nnp=0
      il=il+nw
      loc=il
      do while (nnp.lt.np)
         if (loc+npage.gt.ntmp)&
           call error('getgfl','storage exceeded',' ')
         if (ltt.eq.1) call listio(nin,0,0,tmp(loc),nb,nw)
         if (ltt.eq.2) call tab1io(nin,0,0,tmp(loc),nb,nw)
         loc=loc+nw
         do while (nb.ne.0)
            if (loc+npage.gt.ntmp)&
              call error('getgfl','storage exceeded',' ')
            call moreio(nin,0,0,tmp(loc),nb,nw)
            loc=loc+nw
         enddo
         nnp=nnp+1
         if (nnp.eq.1.and.tmp(il+1).lt.enext) enext=tmp(il+1)
         il=loc
      enddo
   enddo
   na=il-1
   allocate(flgn(na))
   do i=1,na
      flgn(i)=tmp(i)
   enddo
   deallocate(tmp)

   !--sort locations into descending order
   if (ng.gt.1) then
      nh=ng-ni
      do j=1,nh
         nipj=ni+j
         if (ni.eq.0) nipj=ng
         do i=1,nipj
            il=1+loca(i)
            im=1+loca(nipj)
            if (flgn(il).le.flgn(im)) then
               m=-1
               ip1=i+1
               do k=ip1,nipj
                  m=m+1
                  loca(nipj-m)=loca(nipj-m-1)
               enddo
               loca(i)=im-1
            endif
         enddo
      enddo
   endif
   return

   !--normal entry
  200 continue
   enext=emax
   idis=0
   do 260 j=1,ng
   if (li.eq.1) go to 205
   l=1+loca(j)
   lint=l
   lnow=nint(flgn(l+4))
   if (lnow.gt.0) go to 220
   ! for isotropic photons and and energies less than threshold
  205 continue
   gfl(1,j)=1
   if (nl.eq.1) go to 260
   do i=2,nl
      gfl(i,j)=0
   enddo
   go to 260
   ! for anisotropic photons
  220 continue
   np=nint(flgn(l+5))
   im=1
   in=1
   lnext=lnow*2
  240 continue
   l=l+lnext+6
   elo=flgn(l+1)
   if (elo.gt.ed*(1+small).and.im.eq.1) go to 205
   if (enext.lt.elo*(1-small)) enext=elo
   lnext=nint(flgn(l+4))
   if (ltt.eq.2) lnext=2*nint(flgn(l+4))+2*nint(flgn(l+5))
   ehi=flgn(l+lnext+7)
   im=im+1
   if (im.gt.np) call error('getgfl',&
     'desired energy at highest given energy.',' ')
   nbt=nint(flgn(lint+4+2*in))
   if (im.gt.nbt) in=in+1
   if (ehi.le.ed*(1+small).and.im.lt.np) go to 240
   if (im.eq.np.and.ed.gt.up*ehi) call error('getgfl',&
     'desired energy at highest given energy.',' ')
   int=nint(flgn(lint+5+2*in))
   nld=nl
   call getco(alo,nld,lcd,flgn(l),lct,ltt,idis)
   l=l+lnext+6
   nld=nl
   call getco(ahi,nld,lcd,flgn(l),lct,ltt,idis)
   do i=1,nl
      call terp1(elo,alo(i),ehi,ahi(i),ed,gfl(i,j),int)
   enddo
   if (ehi.gt.enext*(1+small)) go to 260
   idis=0
   if (int.eq.1) idis=1
   enext=ehi
  260 continue
   return
   end subroutine getgfl

   subroutine getco(fl,nl,lcd,c,lct,ltt,idis)
   !-------------------------------------------------------------------
   ! Retrieve or compute Legendre coefficients using raw data in
   ! ENDF format.  Conversion between reference systems and
   ! conversion of tabulated data to coefficients is performed by
   ! direct Gaussian quadrature of order 64. Significance of the
   ! results is reduced to reflect precision of calculation.
   !      nl   on entry, maximum number of coefficients desired
   !      nl   on return, number of coefficients required
   !      lcd  system desired (1 for lab, 2 for cm)
   !      lct  system of raw data
   !      ltt  format of raw data (1 for coefficients, 2 for tabulated)
   !      awr  atomic weight ratio
   ! For charged-particle elastic scattering, getco returns the
   ! "nuclear+interference" cross section only.  The Coulomb
   ! part is added in getdis to get a better energy-dependence.
   !-------------------------------------------------------------------
   use util    ! provides error
   use endf    ! provides terpa,terp1,endf variables
   use mathm   ! provideds legndr
   use physics ! provides amassn,amu,ev,clight,hbar
   ! externals
   integer::nl,lcd,lct,ltt,idis
   real(kr)::fl(*),c(*)
   ! internals
   real(kr)::e,ai,x,ee,cc2,eta,sigr,sigi,sgn,xn,y
   real(kr)::xnext,zt,zi
   integer::il,ipc,irc,ipm,intn,iq,ip,l,lang,np,nlz,nmu
   integer::lidp,nt,it,j
   complex(kr)::carg1,carg2,cs1,cs2
   real(kr)::p(65)
   real(kr),dimension(64),parameter::qp=(/&
     -9.99305042e-01_kr,-9.96340117e-01_kr,-9.91013371e-01_kr,&
     -9.83336254e-01_kr,-9.73326828e-01_kr,-9.61008800e-01_kr,&
     -9.46411375e-01_kr,-9.29569172e-01_kr,-9.10522137e-01_kr,&
     -8.89315446e-01_kr,-8.65999398e-01_kr,-8.40629296e-01_kr,&
     -8.13265315e-01_kr,-7.83972359e-01_kr,-7.52819907e-01_kr,&
     -7.19881850e-01_kr,-6.85236313e-01_kr,-6.48965471e-01_kr,&
     -6.11155355e-01_kr,-5.71895646e-01_kr,-5.31279464e-01_kr,&
     -4.89403146e-01_kr,-4.46366017e-01_kr,-4.02270158e-01_kr,&
     -3.57220158e-01_kr,-3.11322872e-01_kr,-2.64687162e-01_kr,&
     -2.17423644e-01_kr,-1.69644420e-01_kr,-1.21462819e-01_kr,&
     -7.29931218e-02_kr,-2.43502927e-02_kr, 2.43502927e-02_kr,&
      7.29931218e-02_kr, 1.21462819e-01_kr, 1.69644420e-01_kr,&
      2.17423644e-01_kr, 2.64687162e-01_kr, 3.11322872e-01_kr,&
      3.57220158e-01_kr, 4.02270158e-01_kr, 4.46366017e-01_kr,&
      4.89403146e-01_kr, 5.31279464e-01_kr, 5.71895646e-01_kr,&
      6.11155355e-01_kr, 6.48965471e-01_kr, 6.85236313e-01_kr,&
      7.19881850e-01_kr, 7.52819907e-01_kr, 7.83972359e-01_kr,&
      8.13265315e-01_kr, 8.40629296e-01_kr, 8.65999398e-01_kr,&
      8.89315446e-01_kr, 9.10522137e-01_kr, 9.29569172e-01_kr,&
      9.46411375e-01_kr, 9.61008800e-01_kr, 9.73326828e-01_kr,&
      9.83336254e-01_kr, 9.91013371e-01_kr, 9.96340117e-01_kr,&
      9.99305042e-01_kr/)
   real(kr),dimension(64),parameter::qw=(/&
      1.78328072e-03_kr, 4.14703326e-03_kr, 6.50445797e-03_kr,&
      8.84675983e-03_kr, 1.11681395e-02_kr, 1.34630479e-02_kr,&
      1.57260305e-02_kr, 1.79517158e-02_kr, 2.01348232e-02_kr,&
      2.22701738e-02_kr, 2.43527026e-02_kr, 2.63774697e-02_kr,&
      2.83396726e-02_kr, 3.02346571e-02_kr, 3.20579284e-02_kr,&
      3.38051618e-02_kr, 3.54722133e-02_kr, 3.70551285e-02_kr,&
      3.85501532e-02_kr, 3.99537411e-02_kr, 4.12625632e-02_kr,&
      4.24735151e-02_kr, 4.35837245e-02_kr, 4.45905582e-02_kr,&
      4.54916279e-02_kr, 4.62847966e-02_kr, 4.69681828e-02_kr,&
      4.75401657e-02_kr, 4.79993886e-02_kr, 4.83447622e-02_kr,&
      4.85754674e-02_kr, 4.86909570e-02_kr, 4.86909570e-02_kr,&
      4.85754674e-02_kr, 4.83447622e-02_kr, 4.79993886e-02_kr,&
      4.75401657e-02_kr, 4.69681828e-02_kr, 4.62847966e-02_kr,&
      4.54916279e-02_kr, 4.45905582e-02_kr, 4.35837245e-02_kr,&
      4.24735151e-02_kr, 4.12625632e-02_kr, 3.99537411e-02_kr,&
      3.85501532e-02_kr, 3.70551285e-02_kr, 3.54722133e-02_kr,&
      3.38051618e-02_kr, 3.20579284e-02_kr, 3.02346571e-02_kr,&
      2.83396726e-02_kr, 2.63774697e-02_kr, 2.43527026e-02_kr,&
      2.22701738e-02_kr, 2.01348232e-02_kr, 1.79517158e-02_kr,&
      1.57260305e-02_kr, 1.34630479e-02_kr, 1.11681395e-02_kr,&
      8.84675983e-03_kr, 6.50445797e-03_kr, 4.14703326e-03_kr,&
      1.78328072e-03_kr/)
   real(kr),parameter::toler=1.e-8_kr
   integer,parameter::nqp=64
   integer,parameter::nlmax=65
   real(kr),parameter::uno=1
   real(kr),parameter::zero=0

   l=nl-1
   if (nl.gt.nlmax)&
     call error('getco','limited to 64 legendre coefficients.',' ')
   lang=-1
   if (mfh.eq.6) lang=l1h
   if (lcd.ne.lct) go to 300
   if (lang.eq.1.or.lang.gt.2) go to 300
   if (ltt.eq.2) go to 300

   !--retrieve coefficients from raw data
   if (lang.ne.2) then
      np=nint(c(5))+1
      fl(1)=1
      nlz=1
      do il=2,nl
         if (il.le.np) then
            fl(il)=c(il+5)
            if (fl(il).ne.zero) nlz=il
         else
            fl(il)=0
         endif
      enddo
      nl=nlz
   else
      np=nint(c(5))
      nlz=1
      do il=1,nl
         if (il.le.np) then
            fl(il)=c(il+6)
            if (fl(il).ne.zero) nlz=il
         else
            fl(il)=0
         endif
      enddo
      nl=nlz
   endif
   return

   !--integrate for fl in desired system
  300 continue
   do il=1,nl
      fl(il)=0
   enddo
   ipc=2
   irc=1
   ipm=0
   if (lang.gt.2) then
      nmu=n2h
      ipm=6+2*nmu
      intn=l1h-10
   endif
   do 390 iq=1,nqp
   x=qp(iq)
   if (ltt.eq.2) go to 350
   if (lang.gt.2) go to 350
   if (lang.eq.1) go to 360

   !--compute scattering probability from raw coefficients
   np=nint(c(5))
   call legndr(x,p,np)
   y=1
   y=y/2
   do ip=1,np
      y=y+(2*ip+1)*p(ip+1)*c(ip+6)/2
   enddo
   go to 370

   !--retrieve scattering probability from raw tabular data
  350 continue
   if (mfh.eq.6) go to 352
   call terpa(y,x,xnext,idis,c,ipc,irc)
   go to 370
  352 continue
   if (x.lt.c(7)) go to 356
   if (ipc.gt.ipm) go to 356
   if (x.le.c(ipc+7)) go to 354
   ipc=ipc+2
   go to 352
  354 continue
   call terp1(c(ipc+5),c(ipc+6),c(ipc+7),c(ipc+8),x,y,intn)
   go to 370
  356 continue
   y=0
   go to 370

   !--coulomb nuclear+interference by amplitude expansion
  360 continue
   if (izap.le.1.or.mth.ne.2) go to 370
   if (lang.ne.1) go to 370
   e=c(2)
   lidp=nint(c(4))
   ai=awrp*amassn
   zt=int(izat/1000)
   zi=int(izap/1000)
   ee=(ev/10000000)*(clight/10)
   cc2=ee**4*amu/(2*hbar**2*ev)
   eta=zt*zi*sqrt(cc2*ai/e)
   nt=nint(c(6))
   np=2*nt
   call legndr(x,p,np)
   if (lidp.ne.1) then
      sigr=c(7)/2
      do ip=1,np
         sigr=sigr+(2*ip+1)*p(ip+1)*c(ip+7)/2
      enddo
      cs1=cmplx(c(8+np),c(9+np))/2
      do it=1,nt
         cs1=cs1+(2*it+1)*p(it+1)&
           *cmplx(c(8+np+2*it),c(9+np+2*it))/2
      enddo
      carg1=cmplx(zero,uno)*eta*log((1-x)/2)
      sigi=(-2*eta/(1-x))*real(cs1*exp(carg1))
      y=(1-x)*(sigr+sigi)
   else
      sigr=c(7)/2
      do it=1,nt
         sigr=sigr+(4*it+1)*p(2*it+1)*c(it+7)/2
      enddo
      cs1=cmplx(c(8+nt),c(9+nt))/2
      cs2=cs1
      sgn=-1
      do it=1,nt
         cs1=cs1+(2*it+1)*p(it+1)&
           *cmplx(c(8+nt+2*it),c(9+nt+2*it))/2
         cs2=cs2+sgn*(2*it+1)*p(it+1)&
           *cmplx(c(8+nt+2*it),c(9+nt+2*it))/2
         sgn=-sgn
      enddo
      carg1=cmplx(zero,uno)*eta*log((1-x)/2)
      carg2=cmplx(zero,uno)*eta*log((1+x)/2)
      sigi=(-2*eta/(1-x*x))*real(cs1*(1+x)*exp(carg1)&
         +cs2*(1-x)*exp(carg2))
      y=(1-x*x)*(sigr+sigi)
   endif

   !--multiply scattering probability by legendre polynomials in
   !--desired system and accumulate
  370 continue
   xn=x
   if (lcd.eq.2.and.lct.eq.1)&
     call error('getco','lab to cm conversion not coded.',' ')
   if (lcd.eq.1.and.lct.eq.2) xn=(1+awr*x)/sqrt(1+awr*awr+2*awr*x)
   call legndr(xn,p,l)
   y=y*qw(iq)
   do il=1,nl
      fl(il)=fl(il)+y*p(il)
   enddo
  390 continue

   !--reduce number of significant figures
   !--and count number of coefficients required
   nlz=1
   do il=2,nl
      j=nint(fl(il)/toler)
      fl(il)=j*toler
      if (j.ne.0) nlz=il
   enddo
   nl=nlz
   return
   end subroutine getco

   subroutine getgyl(ed,enext,idis,gyl,eyl,nyl,matd,mfd,mtd,nin)
   !-------------------------------------------------------------------
   ! Retrieve gamma yields from mf12 or calculate them from mf13.
   ! Routine returns vectors of yields and associated gamma energy.
   ! For fission and radiative capture, the initialization call
   ! returns the energy where significant energy-dependence starts.
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use util ! provides error
   ! externals
   integer::idis,nyl,matd,mfd,mtd,nin
   real(kr)::ed,enext,gyl(*),eyl(*)
   ! internals
   integer::il,ik,ipd,ird,idone,i,ip,ir,idisc,nb,nw,lo,ntmp,na
   integer::nr,l,lp,np,ll,j,ngi,loc
   real(kr)::awr,en,ey,e,sig,rsig
   integer,parameter::nylmax=550
   integer::loca(nylmax)
   real(kr),dimension(:),allocatable::tmp
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::zero=0
   save awr,loca,ird,ipd

   !--initialize
   if (ed.eq.zero) then
      ntmp=99000
      allocate(tmp(ntmp))
      call findf(matd,mfd,mtd,nin)
      call contio(nin,0,0,tmp,nb,nw)
      awr=tmp(2)
      lo=nint(tmp(3))
      if (lo.eq.2) call error('getgyl','lo=2 not coded.',' ')
      nyl=nint(tmp(5))
      if (nyl.gt.nylmax) call error('getgyl','too many gammas.',' ')
      ! find starting location for each sub section
      loc=1
      call tab1io(nin,0,0,tmp(loc),nb,nw)
      if (nyl.ne.1.or.mfd.ne.13) then
         loc=loc+nw
         do while (nb.ne.0)
            call moreio(nin,0,0,tmp(loc),nb,nw)
            loc=loc+nw
         enddo
      endif
      il=loc
      ik=0
      eyl(1)=tmp(1)
      nr=nint(tmp(5))
      enext=tmp(7+2*nr)
      loca(1)=0
      if (nyl.eq.1) then
         na=loc-1
         if (mfd.eq.13) na=1
      else
         do while (ik.lt.nyl)
            if (loc+npage.gt.ntmp)&
              call error('getgyl','storage exceeded.',' ')
            call tab1io(nin,0,0,tmp(loc),nb,nw)
            loc=loc+nw
            do while (nb.ne.0)
               if (loc+npage.gt.ntmp)&
                 call error('getgyl','storage exceeded.',' ')
               call moreio(nin,0,0,tmp(loc),nb,nw)
               loc=loc+nw
            enddo
            ik=ik+1
            loca(ik)=il-1
            nr=nint(tmp(il+4))
            en=tmp(il+6+2*nr)
            if (en.lt.enext) enext=en
            il=loc
         enddo
         il=1+loca(ik)
         eyl(1)=tmp(il)
         ipd=2
         ird=1
         na=loc-1
      endif
      ! check for constant range at low energies
      if (mfd.eq.12) then
         idone=0
         i=0
         do while (i.lt.nyl.and.idone.eq.0)
            i=i+1
            l=1+loca(i)
            ey=tmp(l)
            lp=nint(tmp(l+2))
            nr=nint(tmp(l+4))
            np=nint(tmp(l+5))
            ll=l+5+2*nr
            e=0
            !--if histogram interpolation, start with e equal
            !  to the second energy point
            if (nint(tmp(l+7)).eq.1) e=tmp(ll+3)
            do j=2,np
               if (tmp(ll+2*j).eq.tmp(ll+2)) e=tmp(ll+2*j-1)
            enddo
            if (lp.eq.2.and.ey/100.lt.e) e=ey/100
            if (e.lt.econst) econst=e
            if (econst.eq.zero) idone=1
         enddo
      else
         econst=zero
      endif
      allocate(gyln(na))
      do i=1,na
         gyln(i)=tmp(i)
      enddo
      deallocate(tmp)
      return
   endif

   !--normal entry
   enext=emax
   idis=0
   if (nyl.ne.1.or.mfd.ne.13) then
      do i=1,nyl
         l=1+loca(i)
         eyl(i)=gyln(l)
         lp=nint(gyln(l+2))
         if (lp.eq.2) eyl(i)=eyl(i)+ed*awr/(awr+1)
         ip=2
         ir=1
         call terpa(gyl(i),ed,en,idisc,gyln(l),ip,ir)
         if (en.eq.enext.and.idisc.gt.idis) idis=idisc
         if (en.lt.enext) idis=idisc
         if (en.lt.enext) enext=en
         if (lp.eq.2) then
            ngi=ngg+1
            if (gyln(l).ne.zero) then
               do j=1,ngi
                  e=(egg(j)-gyln(l))*(awr+1)/awr
                  if (e.ge.ed.and.e.lt.enext) then
                     idis=1
                     enext=e
                  endif
               enddo
            endif
         endif
      enddo
      if (mfd.ne.12) then
         ! divide by cross section to get yield
         call terpa(sig,ed,en,idis,gyln,ipd,ird)
         if (en.eq.enext.and.idisc.gt.idis) idis=idisc
         if (en.lt.enext) idis=idisc
         if (en.lt.enext) enext=en
         if (sig.ne.zero) then
            rsig=1/sig
            do i=1,nyl
               gyl(i)=gyl(i)*rsig
            enddo
         endif
      endif
   else
      gyl(1)=1
      eyl(1)=gyln(1)
   endif
   return
   end subroutine getgyl

   subroutine gam102(ans,ed,enext,disc102,law,nl,iglo,ng2,nq)
   !-------------------------------------------------------------------
   ! Process the relativistic discrete gamma or its recoil as
   ! given in mf6/mt102 for ENDF/B-VII neutron + H-1.
   !-------------------------------------------------------------------
   ! externals
   integer::law,nl,iglo,ng2,nq
   real(kr)::ed,enext,disc102
   real(kr)::ans(nl,*)
   ! internals
   integer::i
   real(kr)::edis
   real(kr),parameter::big=1.e10_kr

   !--approximate using discrete gamma for now
   if (law.eq.2) then
      edis=disc102+ed*awr/(awr+1)
      do i=1,ngg
         if (edis.ge.egg(i).and.edis.lt.egg(i+1)) iglo=i
      enddo
      ng2=1
      do i=1,nl
         ans(i,1)=0
      enddo
      ans(1,1)=1
      nq=0
      if (ed.lt.egg(iglo+1)) then
         enext=egg(iglo+1)
      else if (ed.lt.egg(iglo)) then
         enext=egg(iglo)
      else
         enext=big
      endif
   else if (law.eq.4) then
      edis=ed/(awr+1)
      do i=1,ngg
         if (edis.ge.egn(i).and.edis.lt.egn(i+1)) iglo=i
      enddo
      ng2=1
      do i=1,nl
         ans(i,1)=1
      enddo
      nq=0
      if (ed.lt.egn(iglo+1)) then
         enext=egn(iglo+1)
      else if (ed.lt.egn(iglo)) then
         enext=egn(iglo)
      else
         enext=big
      endif
   endif
   return
   end subroutine gam102

   subroutine conver(nin,nout,nscr)
   !-------------------------------------------------------------------
   ! Convert photon transition probability arrays (lo=2), if any,
   ! to photon yields (lo=1).  Add mt456 if necessary.
   ! copy all other sections in this material.
   !-------------------------------------------------------------------
   use util ! provides error,mess
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nout,nscr
   ! internals
   integer::l2flg,imax,lmax,i,istor,ii,idone
   integer::ik,idis,ie,il,iis,nnth,jdis
   integer::imf4,imf6,imf10,imf12,imf13,imf18,imf21,imf22,imf23
   integer::imf24,imf25,imf26,irr21,irr22,irr23,irr24,irr25,irr26
   integer::jza2,itest,imax2,j,jm1,ja,jb,ip1,im,izan,imf
   integer::mt0,mt0old,nb,nw,jzar,lg,n,kk,k,kp1,l,lm1,mtl,no455
   integer::nnu,lnu,mf,mt,l1,l2,n1,n2,nk,jzap,ne,ltp,nm
   integer::m1,m2,mtnow,mttst,nn
   integer::ltt,lcd,nl
   real(kr)::g,ei,eja,ejb,e1,enext,p,za2,za,ysum,yy,tsave,w,elow,etop,enxt,x
   character(60)::strng
   integer::ngam(440)
   real(kr)::sig(1,1)   !matches parameters passed via call getsig below
   integer,parameter::mxlg=65
   integer,parameter::mxnnth=450
   real(kr)::fl(mxlg)
   integer::mtth(mxnnth)
   real(kr)::eeth(mxnnth)
   real(kr),dimension(:),allocatable::scr
   real(kr),dimension(:),allocatable::e,eg,es,y,aa,r
   real(kr),dimension(:),allocatable::nu
   real(kr),parameter::zero=0

   !--initialize
   l2flg=0
   ndelg=0
!   imax=49
   imax=50
   lmax=500
   do i=1,imax
      ngam(i)=0
   enddo
   istor=0
   allocate(scr(npage+50))
   nsh=0
   nsc=0
   imf4=1
   imf6=1
   imf10=1
   imf12=1
   imf13=1
   imf18=1
   imf21=1
   imf22=1
   imf23=1
   imf24=1
   imf25=1
   imf26=1
   irr21=1
   irr22=1
   irr23=1
   irr24=1
   irr25=1
   irr26=1
   mt0=0
   mt0old=0
   mtnow=-1
   nnth=0
   izatest=0
   lastza=0

   !--loop over all sections of this material
  110 continue
   call contio(nin,0,0,scr,nb,nw)
   if (mfh.eq.0.and.math.ne.0) go to 130
   if (mfh.eq.1) go to 550

   !--get thresholds vs mt number
   if (mfh.eq.3.and.mth.ne.0) then
      ei=0
      call gety1(ei,enxt,jdis,x,nin,scr)
      nnth=nnth+1
      if (nnth.gt.mxnnth) call error('conver','nnth too large',' ')
      eeth(nnth)=enxt
      mtth(nnth)=mth
      call contio(0,nout,nscr,scr,nb,nw)
      call tosend(nin,nout,nscr,scr)
      go to 110
   endif

   if (mfh.eq.12.and.mth.eq.460) then
      call tosend(nin,0,0,scr)
      go to 110
   endif

   if (mfh.eq.6.and.mth.eq.18.and.l1h.ne.0) then
      call mess('conver','skipping new mf6/mt18 multiplicity section',&
                '')
      call tosend(nin,0,0,scr)
      go to 110
   endif

   if (mfh.ne.4) go to 111
   if (imf4.eq.1) go to 211
   if (mth.eq.18) go to 211
   if (mth.eq.19.and.mf4(imf4-1).eq.18) go to 311
   if (mth.gt.iabs(mf4(imf4-1))+1) go to 211
   if (mf4(imf4-1).lt.0) imf4=imf4-1
   mf4(imf4)=-mth
   imf4=imf4+1
   go to 312
  311 continue
   imf4=imf4-1
  211 continue
   mf4(imf4)=mth
   imf4=imf4+1
  312 continue
   za2=-1
   if (mth.eq.2.or.(mth.ge.50.and.mth.lt.91)) za2=1
   if (iverf.lt.6) then
      if (mth.ge.700.and.mth.lt.719) za2=1001
      if (mth.ge.720.and.mth.lt.739) za2=1002
      if (mth.ge.740.and.mth.lt.759) za2=1003
      if (mth.ge.760.and.mth.lt.779) za2=2003
      if (mth.ge.780.and.mth.lt.799) za2=2004
   else
      if (mth.ge.600.and.mth.lt.649) za2=1001
      if (mth.ge.650.and.mth.lt.699) za2=1002
      if (mth.ge.700.and.mth.lt.749) za2=1003
      if (mth.ge.750.and.mth.lt.799) za2=2003
      if (mth.ge.800.and.mth.lt.849) za2=2004
      if (mth.ge.875.and.mth.lt.891) za2=1
   endif
   if (za2.lt.0..and.&
     mth.ne.18.and.mth.ne.19.and.mth.ne.20.and.&
     mth.ne.21.and.mth.ne.38) then
      write(strng,&
        '(''cannot do complete particle production for mt='',&
        &i3)') mth
      call mess('conver',strng,'only mf4/mf5 provided')
      go to 119
   endif
   izat=nint(c1h)
   jza2=nint(za2)
   jzar=izat+izap-jza2
   if (jza2.eq.1001) then
      itest=0
      if (irr21.ne.1) then
         if (mth.le.iabs(mf4r(1,irr21-1))+1) itest=1
      endif
      if (itest.eq.1) then
         if (mf4r(1,irr21-1).lt.0) irr21=irr21-1
         mf4r(1,irr21)=-mth
         irr21=irr21+1
      else
         mf4r(1,irr21)=mth
         irr21=irr21+1
      endif
   else if (jza2.eq.1002) then
      itest=0
      if (irr22.ne.1) then
         if (mth.le.iabs(mf4r(2,irr22-1))+1) itest=1
      endif
      if (itest.eq.1) then
         if (mf4r(2,irr22-1).lt.0) irr22=irr22-1
         mf4r(2,irr22)=-mth
         irr22=irr22+1
      else
         mf4r(2,irr22)=mth
         irr22=irr22+1
      endif
   else if (jza2.eq.1003) then
      itest=0
      if (irr23.ne.1) then
         if (mth.le.iabs(mf4r(3,irr23-1))+1) itest=1
      endif
      if (itest.eq.1) then
         if (mf4r(3,irr23-1).lt.0) irr23=irr23-1
         mf4r(3,irr23)=-mth
         irr23=irr23+1
      else
         mf4r(3,irr23)=mth
         irr23=irr23+1
      endif
   else if (jza2.eq.2003) then
      itest=0
      if (irr24.ne.1) then
         if (mth.le.iabs(mf4r(4,irr24-1))+1) itest=1
      endif
      if (itest.eq.1) then
         if (mf4r(4,irr24-1).lt.0) irr24=irr24-1
         mf4r(4,irr24)=-mth
         irr24=irr24+1
      else
         mf4r(4,irr24)=mth
         irr24=irr24+1
      endif
   else if (jza2.eq.2004) then
      itest=0
      if (irr25.ne.1) then
         if (mth.le.iabs(mf4r(5,irr25-1))+1) itest=1
      endif
      if (itest.eq.1) then
         if (mf4r(5,irr25-1).lt.0) irr25=irr25-1
         mf4r(5,irr25)=-mth
         irr25=irr25+1
      else
         mf4r(5,irr25)=mth
         irr25=irr25+1
      endif
   endif
   if (jzar.gt.2004) then
      itest=0
      if (irr26.ne.1) then
         if (mth.le.iabs(mf4r(6,irr26-1))+1) itest=1
      endif
      if (itest.eq.1) then
         if (mf4r(6,irr26-1).lt.0) irr26=irr26-1
         mf4r(6,irr26)=-mth
         irr26=irr26+1
      else
         mf4r(6,irr26)=mth
         irr26=irr26+1
      endif
   endif
   go to 119
  111 continue
   if (mfh.eq.6) go to 620
   if (mfh.eq.8) go to 820
   if (mfh.eq.12) then
      itest=0
      if (imf12.ne.1) then
         if (mth.le.iabs(mf12(imf12-1))+1) itest=1
      endif
      if (itest.eq.1) then
         if (mf12(imf12-1).lt.0) imf12=imf12-1
         mf12(imf12)=-mth
         imf12=imf12+1
      else
         mf12(imf12)=mth
         imf12=imf12+1
      endif
   else if (mfh.eq.13) then
      itest=0
      if (imf13.ne.1) then
         if (mth.le.iabs(mf13(imf13-1))+1) itest=1
      endif
      if (itest.eq.1) then
         if (mf13(imf13-1).lt.0) imf13=imf13-1
         mf13(imf13)=-mth
         imf13=imf13+1
      else
         mf13(imf13)=mth
         imf13=imf13+1
      endif
   endif
  119 continue
   if (mfh.eq.12) go to 120
   if (math.eq.1149.and.mfh.eq.13.and.mth.ge.51.and.iverf.lt.5)&
     go to 140
   if (math.eq.1150.and.mfh.eq.13.and.mth.le.54.and.iverf.lt.5)&
     go to 140
   if (mfh.eq.14.and.l2flg.eq.1) go to 500
   call contio(0,nout,nscr,scr,nb,nw)
   if (math.eq.0) go to 600
   call tosend(nin,nout,nscr,scr)
   go to 110
  120 continue
   if (l1h.eq.2) go to 150
   call contio(0,nout,nscr,scr,nb,nw)
   call tosend(nin,nout,nscr,scr)
   go to 110
  130 continue
   call contio(0,nout,nscr,scr,nb,nw)
   go to 110
  140 continue
   call tosend(nin,0,0,scr)
   write(strng,'(''gamma production patch made for '',i3)') mth
   call mess('conver',strng,' ')
   go to 110

   !--allocate storage
  150 continue
   if (istor.eq.0) then
      istor=1
      imax2=imax*imax
      allocate(e(imax))
      allocate(eg(lmax))
      allocate(es(lmax))
      allocate(y(lmax))
      allocate(aa(imax2))
      allocate(r(imax2))
      e(1)=0
      do i=1,imax2
         r(i)=0
         aa(i)=0
      enddo
      do i=1,imax
         r(imax*(i-1)+i)=1
      enddo
   endif

   !--check for missing mf12 mt's ... advise user if found
   if (mth.eq.51.or.                                                &
       (iverf.ge.6.and.(mth.eq.601.or.mth.eq.651.or.mth.eq.701.or.  &
                        mth.eq.751.or.mth.eq.801.or.mth.eq.876)).or.&
       (iverf.le.5.and.(mth.eq.701.or.mth.eq.721.or.mth.eq.741.or.  &
                        mth.eq.761.or.mth.eq.781))) then
      mtnow=mth
   elseif ((mth.gt.51.and.mth.lt.91).or.                 &
           (iverf.ge.6.and.mth.gt.601.and.mth.lt.648).or.&
           (iverf.ge.6.and.mth.gt.651.and.mth.lt.698).or.&
           (iverf.ge.6.and.mth.gt.701.and.mth.lt.748).or.&
           (iverf.ge.6.and.mth.gt.751.and.mth.lt.798).or.&
           (iverf.ge.6.and.mth.gt.801.and.mth.lt.848).or.&
           (iverf.ge.6.and.mth.gt.876.and.mth.lt.890).or.&
           (iverf.le.5.and.mth.gt.701.and.mth.lt.718).or.&
           (iverf.le.5.and.mth.gt.721.and.mth.lt.738).or.&
           (iverf.le.5.and.mth.gt.741.and.mth.lt.758).or.&
           (iverf.le.5.and.mth.gt.761.and.mth.lt.778).or.&
           (iverf.le.5.and.mth.gt.781.and.mth.lt.798)) then
      if (mtnow+1.ne.mth) then
         write(strng,'(''mf12, mt'',i2,'' may be missing'')')mth-1
         call mess('conver',strng,&
                            'discrete photon data may be incomplete')
      endif
      mtnow=mth
   endif

   !--convert transition probability array to yields
   za=c1h
   awr=c2h
   lg=l2h
   g=1
   l2flg=1
   call listio(nin,0,0,scr,nb,nw)

   !--make sure mt0 is correct for this range of mt's.
   if (mth.ge.51.and.mth.le.90.and.mt0.ne.49) mt0=49
   if (iverf.ge.6) then
      if (mth.ge.601.and.mth.le.649.and.mt0.ne.599) mt0=599
      if (mth.ge.651.and.mth.le.699.and.mt0.ne.649) mt0=649
      if (mth.ge.701.and.mth.le.749.and.mt0.ne.699) mt0=699
      if (mth.ge.751.and.mth.le.799.and.mt0.ne.749) mt0=749
      if (mth.ge.801.and.mth.le.849.and.mt0.ne.799) mt0=799
      if (mth.ge.876.and.mth.le.891.and.mt0.ne.874) mt0=874
   elseif (iverf.le.5) then
      if (mth.ge.701.and.mth.le.719.and.mt0.ne.699) mt0=699
      if (mth.ge.721.and.mth.le.739.and.mt0.ne.719) mt0=719
      if (mth.ge.741.and.mth.le.759.and.mt0.ne.739) mt0=739
      if (mth.ge.761.and.mth.le.779.and.mt0.ne.759) mt0=759
      if (mth.ge.781.and.mth.le.799.and.mt0.ne.779) mt0=779
   endif

   !--load the e(_) array with mf3 -q.  Will overwrite with mf12
   !  mt data when possible (which should be always but if there
   !  is a missing mf12 mt, we're covered).
   if (mt0.ne.mt0old) then
      e=0
      mt0old=mt0
      m1=mt0+2
      if (mt0.eq.49) then
         m2=91
      elseif (iverf.ge.6.and.mt0.ne.874) then
         m2=m1+48
      elseif (iverf.ge.6.and.mt0.eq.874) then
         m2=m1+15
      elseif (iverf.le.5) then
         m2=m1+17
      endif
      nn=1
      mttst=-1
      do while (mttst.lt.m2.and.nn.le.nnth)
         mttst=mtth(nn)
         if (mttst.ge.m1.and.mttst.le.m2)&
                                  e(mttst-mt0)=eeth(nn)*awr/(awr+1)
         nn=nn+1
      enddo
   endif

   j=mth-mt0
   e(j)=scr(1)
   n=n2h
   jm1=j-1
   do kk=1,jm1
      k=jm1-kk+1
      ii=0
      idone=0
      do while (ii.lt.n.and.idone.eq.0)
         ii=ii+1
         i=ii
         ei=scr(6+(lg+1)*i-lg)
         if (ei.eq.zero.and.e(k).eq.zero) idone=1
         if (ei.ne.zero) then
            if (abs(ei-e(k))/ei.lt.0.0001) idone=1
         endif
      enddo
      if (idone.eq.0) then
         aa((k-1)*imax+j)=0
         r((k-1)*imax+j)=0
      else
         p=scr(7+(lg+1)*i-lg)
         g=1
         if (lg.eq.2) g=scr(6+(lg+1)*i)
         aa((k-1)*imax+j)=p*g
         r((k-1)*imax+j)=p
      endif
      if (k.ne.jm1) then
         kp1=k+1
         do ii=kp1,jm1
            r((k-1)*imax+j)=r((k-1)*imax+j)&
              +r((ii-1)*imax+j)*aa((k-1)*imax+ii)
         enddo
      endif
   enddo
   l=0
   ysum=0
   do i=2,j
      ja=j+2-i
      eja=e(ja)
      do ii=1,jm1
         jb=j-ii
         ejb=e(jb)
         yy=aa((jb-1)*imax+ja)*r((ja-1)*imax+j)
         if (yy.ne.zero) then
            l=l+1
            if (l.gt.lmax) call error('conver',&
              'too many lo=2 gammas',' ')
            eg(l)=eja-ejb
            es(l)=eja
            y(l)=yy
            ysum=ysum+y(l)
         endif
      enddo
   enddo
   e(1)=0

   !--arrange gamma energies in descending order
   if (l.gt.1) then
      lm1=l-1
      do i=1,lm1
         ip1=i+1
         do ii=ip1,l
            if (eg(i).lt.eg(ii)) then
               tsave=eg(i)
               eg(i)=eg(ii)
               eg(ii)=tsave
               tsave=y(i)
               y(i)=y(ii)
               y(ii)=tsave
               tsave=es(i)
               es(i)=es(ii)
               es(ii)=tsave
            endif
         enddo
      enddo
   endif

   !--output computed yields in endf lo=1 format
   scr(1)=za
   scr(2)=awr
   scr(3)=1
   scr(4)=0
   scr(5)=l
   scr(6)=0
   mtl=49
   if (iverf.le.5.and.mth.ge.700) mtl=649
   if (iverf.ge.6.and.mth.ge.600) mtl=549
   ngam(mth-mtl)=l
   call contio(0,nout,nscr,scr,nb,nw)
   etop=emaxx
   elow=ebeg
   do i=1,nnth
      if (mtth(i).eq.mth) elow=eeth(i)
   enddo
   ! output tab1 sum record, if needed
   if (l.gt.1) then
      scr(1)=0
      scr(2)=0
      scr(3)=0
      scr(5)=1
      scr(6)=2
      scr(7)=2
      scr(8)=2
      scr(9)=elow
      scr(10)=ysum
      scr(11)=etop
      scr(12)=ysum
      nw=12
      call tab1io(0,nout,nscr,scr,nb,nw)
   endif
   ! output rest of tab1 records
   do i=1,l
      scr(1)=eg(i)
      scr(2)=es(i)
      scr(3)=0
      scr(4)=2
      scr(5)=1
      scr(6)=2
      scr(7)=2
      scr(8)=2
      scr(9)=elow
      scr(10)=y(i)
      scr(11)=etop
      scr(12)=y(i)
      nw=12
      call tab1io(0,nout,nscr,scr,nb,nw)
   enddo

   !--copy send record and loop over remaining sections.
   call contio(nin,nout,nscr,scr,nb,nw)
   go to 110

   !--angular distributions
  500 continue
   mt0=0
   if (mth.ge.50.and.mth.le.90) mt0=50
   if (iverf.le.5.and.mth.ge.700.and.mth.le.717) mt0=700
   if (iverf.le.5.and.mth.ge.720.and.mth.le.737) mt0=720
   if (iverf.le.5.and.mth.ge.740.and.mth.le.757) mt0=740
   if (iverf.le.5.and.mth.ge.760.and.mth.le.777) mt0=760
   if (iverf.le.5.and.mth.ge.780.and.mth.le.797) mt0=780
   if (iverf.ge.6.and.mth.ge.600.and.mth.le.648) mt0=600
   if (iverf.ge.6.and.mth.ge.650.and.mth.le.698) mt0=650
   if (iverf.ge.6.and.mth.ge.700.and.mth.le.748) mt0=700
   if (iverf.ge.6.and.mth.ge.750.and.mth.le.798) mt0=750
   if (iverf.ge.6.and.mth.ge.800.and.mth.le.848) mt0=800
   if (iverf.ge.6.and.mth.ge.875.and.mth.le.890) mt0=875
   if (mt0.ne.0) then
      mtl=49
      if (iverf.le.5.and.mth.ge.700) mtl=649
      if (iverf.ge.6.and.mth.ge.600) mtl=549
      scr(3)=1
      scr(4)=0
      scr(6)=0
      scr(5)=ngam(mth-mtl)
      call contio(0,nout,nscr,scr,nb,nw)
      call tosend(nin,0,0,scr)
      call asend(nout,nscr)
   else
      call contio(0,nout,nscr,scr,nb,nw)
      call tosend(nin,nout,nscr,scr)
   endif
   go to 110

   !--file1.  add 456=452 if necessary.
  550 continue
   call contio(0,nout,nscr,scr,nb,nw)
   izap=1
   awrp=1
   if (iverf.ge.6) then
      call contio(nin,nout,nscr,scr,nb,nw)
      call contio(nin,nout,nscr,scr,nb,nw)
      izap=n1h/10
      awrp=c1h
   endif
   call tosend(nin,nout,nscr,scr)
   no455=0
   call contio(nin,nout,nscr,scr,nb,nw)
   if (mfh.eq.0) go to 110
   if (mth.ne.452) go to 595
   nnu=8000
   allocate(nu(nnu))
   l=1
   lnu=l2h
   if (lnu.eq.1) call listio(nin,nout,nscr,nu(l),nb,nw)
   if (lnu.eq.2) call tab1io(nin,nout,nscr,nu(l),nb,nw)
   do while (nb.ne.0)
      if (l+nw.gt.nnu) call error('conver',&
        'storage for fission nu exceeded',' ')
      l=l+nw
      call moreio(nin,nout,nscr,nu(l),nb,nw)
   enddo
   call tosend(nin,nout,nscr,scr)
  565 continue
   call contio(nin,0,0,scr,nb,nw)
   if (mth.eq.456) go to 590
   if (mth.eq.0.and.no455.eq.0) go to 570
   if (mfh.eq.0) go to 590
   if (mth.eq.455) no455=1
   if (mth.gt.456.and.no455.eq.0) go to 570
   call contio(0,nout,nscr,scr,nb,nw)
   if (mth.eq.455) then
      call listio(nin,nout,nscr,scr,nb,nw)
      ndelg=n1h
   endif
   call tosend(nin,nout,nscr,scr)
   go to 565
  570 continue
   if (no455.eq.1) go to 590
   mf=mfh
   mt=mth
   l1=l1h
   l2=l2h
   n1=n1h
   n2=n2h
   mfh=1
   mth=456
   scr(3)=0
   scr(4)=lnu
   scr(5)=0
   scr(6)=0
   l=1
   call contio(0,nout,nscr,scr,nb,nw)
   if (lnu.eq.1) call listio(0,nout,nscr,nu(l),nb,nw)
   if (lnu.eq.2) call tab1io(0,nout,nscr,nu(l),nb,nw)
   do while (nb.ne.0)
      l=l+nw
      call moreio(0,nout,nscr,nu(l),nb,nw)
   enddo
   call asend(nout,nscr)
   mfh=mf
   mth=mt
   scr(3)=l1
   scr(4)=l2
   scr(5)=n1
   scr(6)=n2
  590 continue
   call contio(0,nout,nscr,scr,nb,nw)
   deallocate(nu)
   if (mfh.eq.0) go to 110
  595 continue
   call tofend(nin,nout,nscr,scr)
   go to 110

   !--examine contents of file 6.
  620 continue
   call contio(0,nout,nscr,scr,nb,nw)
   nk=n1h
   ik=0
  625 continue
   ik=ik+1
   call tab1io(nin,nout,nscr,scr,nb,nw)
   do while (nb.ne.0)
      call moreio(nin,nout,nscr,scr,nb,nw)
   enddo
   jzap=nint(c1h)
   law=l2h

   !--convert cp elastic ni format
   if (law.eq.5) then
      izap=-izap
      matd=math
      mfd=mfh
      mtd=mth
      e1=0
      call getsig(e1,enext,idis,sig,1,1)
      call tab2io(nin,nout,nscr,scr,nb,nw)
      spi=c1h
      ne=n2h
      do ie=1,ne
         call listio(nin,0,0,scr,nb,nw)
         ltp=l1h
         if (ltp.le.2) then
            call listio(0,nout,nscr,scr,nb,nw)
         else
            e1=c2h
            call getsig(e1,enext,idis,sig,1,1)
            nm=n2h
            do im=1,nm
               w=scr(5+2*im)
               scr(6+2*im)=sig(1,1)*(1-w)*scr(6+2*im)
            enddo
            ltt=2
            lcd=2
            lct=2
            nl=21
            call getco(fl,nl,lcd,scr,lct,ltt,idis)
            scr(3)=2
            scr(5)=nl
            scr(6)=nl
            do il=1,nl
               scr(6+il)=fl(il)
            enddo
            call listio(0,nout,nscr,scr,nb,nw)
         endif
      enddo
      izap=-izap

   !--just copy this subsection
   else
      call skip6(nin,nout,nscr,scr,law)
   endif

   !--record reaction information
   if (jzap.ne.0) go to 725
   if (imf18.eq.1) go to 721
   if (mth.eq.iabs(mf18(imf18-1))) go to 790
   if (mth.gt.iabs(mf18(imf18-1))+1) go to 721
   if (mf18(imf18-1).lt.0) imf18=imf18-1
   mf18(imf18)=-mth
   imf18=imf18+1
   go to 790
  721 continue
   mf18(imf18)=mth
   imf18=imf18+1
   go to 790
  725 continue
   if (jzap.ne.1) go to 730
   if (imf6.eq.1) go to 726
   if (mth.eq.iabs(mf6(imf6-1))) go to 790
   if (mth.gt.iabs(mf6(imf6-1))+1) go to 726
   if (mf6(imf6-1).lt.0) imf6=imf6-1
   mf6(imf6)=-mth
   imf6=imf6+1
   go to 790
  726 continue
   mf6(imf6)=mth
   imf6=imf6+1
   go to 790
  730 continue
   if (jzap.ne.1001) go to 735
   if (imf21.eq.1) go to 731
   if (mth.eq.iabs(mf6p(1,imf21-1))) go to 790
   if (mth.gt.iabs(mf6p(1,imf21-1))+1) go to 731
   if (mf6p(1,imf21-1).lt.0) imf21=imf21-1
   mf6p(1,imf21)=-mth
   imf21=imf21+1
   go to 790
  731 continue
   mf6p(1,imf21)=mth
   imf21=imf21+1
   go to 790
  735 continue
   if (jzap.ne.1002) go to 740
   if (imf22.eq.1) go to 736
   if (mth.eq.iabs(mf6p(2,imf22-1))) go to 790
   if (mth.gt.iabs(mf6p(2,imf22-1))+1) go to 736
   if (mf6p(2,imf22-1).lt.0) imf22=imf22-1
   mf6p(2,imf22)=-mth
   imf22=imf22+1
   go to 790
  736 continue
   mf6p(2,imf22)=mth
   imf22=imf22+1
   go to 790
  740 continue
   if (jzap.ne.1003) go to 745
   if (imf23.eq.1) go to 741
   if (mth.eq.iabs(mf6p(3,imf23-1))) go to 790
   if (mth.gt.iabs(mf6p(3,imf23-1))+1) go to 741
   if (mf6p(3,imf23-1).lt.0) imf23=imf23-1
   mf6p(3,imf23)=-mth
   imf23=imf23+1
   go to 790
  741 continue
   mf6p(3,imf23)=mth
   imf23=imf23+1
   go to 790
  745 continue
   if (jzap.ne.2003) go to 750
   if (imf24.eq.1) go to 746
   if (mth.eq.iabs(mf6p(4,imf24-1))) go to 790
   if (mth.gt.iabs(mf6p(4,imf24-1))+1) go to 746
   if (mf6p(4,imf24-1).lt.0) imf24=imf24-1
   mf6p(4,imf24)=-mth
   imf24=imf24+1
   go to 790
  746 continue
   mf6p(4,imf24)=mth
   imf24=imf24+1
   go to 790
  750 continue
   if (jzap.ne.2004) go to 755
   if (imf25.eq.1) go to 751
   if (mth.eq.iabs(mf6p(5,imf25-1))) go to 790
   if (mth.gt.iabs(mf6p(5,imf25-1))+1) go to 751
   if (mf6p(5,imf25-1).lt.0) imf25=imf25-1
   mf6p(5,imf25)=-mth
   imf25=imf25+1
   go to 790
  751 continue
   mf6p(5,imf25)=mth
   imf25=imf25+1
   go to 790
  755 continue
   if (imf26.eq.1) go to 756
   if (mth.eq.iabs(mf6p(6,imf26-1))) go to 790
   if (mth.gt.iabs(mf6p(6,imf26-1))+1) go to 756
   if (mf6p(6,imf26-1).lt.0) imf26=imf26-1
   mf6p(6,imf26)=-mth
   imf26=imf26+1
   go to 790
  756 continue
   mf6p(6,imf26)=mth
   imf26=imf26+1
   go to 790

   !--examine contents of file 8
  820 continue
   nk=n1h
   if (mth.le.200.or.mth.ge.600) then
      ik=0
      do while (ik.lt.nk)
         ik=ik+1
         call listio(nin,0,0,scr,nb,nw)
         izan=nint(c1h)
         imf=l1h
         iis=l2h
         mf10f(imf10)=imf
         mf10s(imf10)=mth
         mf10i(imf10)=10*izan
         if (izan.eq.-1.and.mth.eq.18) mf10i(imf10)=0
         lfs8(imf10)=iis
         if (izan.ne.lastza) then
            if (iis.eq.0) then
               mlfs8(imf10)=0
            else
               mlfs8(imf10)=1
            endif
         else
            mlfs8(imf10)=mlfs8(imf10-1)+1
         endif
         lastza=izan
         imf10=imf10+1
      enddo
   endif
   call tosend(nin,0,0,scr)
   go to 110

   !--continue loop over subsections
  790 continue
   if (ik.lt.nk) go to 625
   call tosend(nin,nout,nscr,scr)
   go to 110

   !--conver is finished.
  600 continue
   if (istor.gt.0) then
      deallocate(e)
      deallocate(eg)
      deallocate(es)
      deallocate(y)
      deallocate(aa)
      deallocate(r)
   endif
   if (imf4.gt.maxr1) call error('conver','imf4 too big',' ')
   if (irr21.gt.maxr1) call error('conver','irr21 too big',' ')
   if (irr22.gt.maxr1) call error('conver','irr22 too big',' ')
   if (irr23.gt.maxr1) call error('conver','irr23 too big',' ')
   if (irr24.gt.maxr1) call error('conver','irr24 too big',' ')
   if (irr24.gt.maxr1) call error('conver','irr24 too big',' ')
   if (irr25.gt.maxr1) call error('conver','irr25 too big',' ')
   if (irr26.gt.maxr1) call error('conver','irr26 too big',' ')
   if (imf12.gt.maxr1) call error('conver','imf12 too big',' ')
   if (imf13.gt.maxr1) call error('conver','imf13 too big',' ')
   if (imf18.gt.maxr1) call error('conver','imf18 too big',' ')
   if (imf6.gt.maxr1) call error('conver','imf6 too big',' ')
   if (imf21.gt.maxr1) call error('conver','imf21 too big',' ')
   if (imf22.gt.maxr1) call error('conver','imf22 too big',' ')
   if (imf23.gt.maxr1) call error('conver','imf23 too big',' ')
   if (imf24.gt.maxr1) call error('conver','imf24 too big',' ')
   if (imf25.gt.maxr1) call error('conver','imf25 too big',' ')
   if (imf26.gt.maxr1) call error('conver','imf26 too big',' ')
   if (imf26.gt.maxr1) call error('conver','imf26 too big',' ')
   if (imf10.gt.maxr2) call error('conver','imf10 too big',' ')
   mf4(imf4)=0
   mf6(imf6)=0
   mf10f(imf10)=0
   mf10s(imf10)=0
   mf10i(imf10)=0
   lfs8(imf10)=0
   mf12(imf12)=0
   mf13(imf13)=0
   mf18(imf18)=0
   mf6p(1,imf21)=0
   mf6p(2,imf22)=0
   mf6p(3,imf23)=0
   mf6p(4,imf24)=0
   mf6p(5,imf25)=0
   mf6p(6,imf26)=0
   mf4r(1,irr21)=0
   mf4r(2,irr22)=0
   mf4r(3,irr23)=0
   mf4r(4,irr24)=0
   mf4r(5,irr25)=0
   mf4r(6,irr26)=0
   return
   end subroutine conver

   subroutine getsed(ed,enext,idis,sed,eg,ng,nk,matd,mfd,mtd,nin)
   !-------------------------------------------------------------------
   ! Compute secondary energy distribution for all sink groups
   ! simultaneously.  Laws 1, 3, 5, 7, 9, and 11 are coded.
   ! Initialize if ed=0 and return number of subsections in nk.
   ! On a normal entry if nk=1, return the sum of all the subsections,
   ! but if nk.gt.1, return the contributions separately.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use util   ! provides mess
   use endf   ! provides endf routines and variables
   ! externals
   integer::idis,ng,nk,matd,mfd,mtd,nin,nc
   real(kr)::ed,enext,sed(nk,*),eg(*)
   ! internals
   integer::nktot,nupm,nb,nw,l,il,ik,lf,ln,ip,ne,nbt,int,nnow,idisc
   integer::nne,iraw,ir,nr,np,nnt,ig,lnow,i,llo,lhi,iout,ikt,mnow
   integer::ntmp,m,m1,klo,nrlo,nplo,khi,nrhi,nphi,jnt
   integer::ier
   real(kr)::elo,delta,ee,e1,e2,ehi,pe,eihi,test,sigup,ethi,s
   real(kr)::xlo,xhi,xend,e1lo,e1hi,e2lo,e2hi,flo,fhi,fe
   real(kr)::val,fx,ex
   integer::mm,ix,l1
   character(60)::strng
   integer,parameter::nkmax=20
   integer::loc(nkmax)
   real(kr),dimension(:),allocatable::tmp
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::eps=.001e0_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::ebig=1.e8_kr
   real(kr),parameter::tenth=0.1e0_kr
   real(kr),parameter::zero=0
   save nktot,nupm,loc

   !--initialize
   if (ed.eq.zero) then
      ier=1
      ntmp=250000
      do while (ier.ne.0)
         if (allocated(tmp)) deallocate(tmp)
         allocate(tmp(ntmp),stat=ier)
         if (ier.ne.0) then
            ntmp=ntmp/2
         endif
      enddo
      call findf(matd,mfd,mtd,nin)
      call contio(nin,0,0,tmp,nb,nw)
      nktot=nint(tmp(5))
      if (nktot.gt.nkmax)&
        call error('getsed','too many subsections.',' ')
      nk=nktot
      l=1
      il=0
      ik=0
      enext=emax
      do while (ik.lt.nk)
         ik=ik+1
         il=l-1
         loc(ik)=il
         call tab1io(nin,0,0,tmp(l),nb,nw)
         l=l+nw
         do while (nb.ne.0)
            if (l+nw.gt.ntmp) call error('getsed',&
            'storage for tmp exceeded1',' ')
            call moreio(nin,0,0,tmp(l),nb,nw)
            l=l+nw
         enddo
         lf=nint(tmp(il+4))
         nr=nint(tmp(il+5))
         elo=tmp(il+1+6+2*nr)
         if (elo.lt.enext) enext=elo
         il=l-1

         !--analytic subsection.
         if (lf.gt.1) then
            if (lf.ne.3) then
               call tab1io(nin,0,0,tmp(l),nb,nw)
               nr=nint(tmp(l+4))
               np=nint(tmp(l+5))
               ln=l+2*nr+5
               l=l+nw
               do while (nb.ne.0)
                  call moreio(nin,0,0,tmp(l),nb,nw)
                  l=l+nw
               enddo
               do ip=2,np
                  delta=abs(tmp(ln+2*ip)-tmp(ln+2*ip-2))
                  if (delta.ge.eps*tmp(ln+2*ip-2)) then
                     ee=tmp(ln+2*ip-3)+eps*tmp(ln+2*ip-2)*&
                       (tmp(ln+2*ip-1)-tmp(ln+2*ip-3))/delta
                     if (ee.lt.econst) econst=ee
                  endif
               enddo
               if (lf.eq.5.or.lf.eq.11) then
                  l1=l
                  call tab1io(nin,0,0,tmp(l),nb,nw)
                  nr=nint(tmp(l+4))
                  np=nint(tmp(l+5))
                  ln=l+2*nr+5
                  l=l+nw
                  do while (nb.ne.0)
                     call moreio(nin,0,0,tmp(l),nb,nw)
                     l=l+nw
                  enddo
                  !extend lowest delayed bin using sqrt(e) shape
                  if (ismooth.gt.0.and.mtd.eq.455.and.&
                    nint(tmp(l1+7)).eq.1) then
                     ex=40
                     fx=.8409
                     write(nsyso,'('' extending lowest delayed bin'',&
                      &'' using sqrt(E)'')')
                     mm=nint(tmp(l1+5))
                     do while (tmp(l1+10).gt.ex)
                        if (l1+9+2*mm.gt.ntmp) call error('getsed',&
                        'storage for tmp exceeded3',' ')
                        do ix=2*mm,1,-1
                           tmp(l1+9+ix)=tmp(l1+7+ix)
                        enddo
                        tmp(l1+10)=fx*tmp(l1+12)
                        val=tmp(l1+9)
                        tmp(l1+9)=sqrt(fx)*val
                        tmp(l1+9)=sigfig(tmp(l1+9),7,0)
                        tmp(l1+11)=(1-fx*sqrt(fx))*val/(1-fx)
                        tmp(l1+11)=sigfig(tmp(l1+11),7,0)
                        mm=mm+1
                     enddo
                     tmp(l1+5)=mm
                     tmp(l1+6)=mm
                     l=l1+8+2*mm
                  endif
                  if (lf.ne.5) then
                     do ip=2,np
                        delta=abs(tmp(ln+2*ip)-tmp(ln+2*ip-2))
                        if (delta.ge.eps*tmp(ln+2*ip-2)) then
                           ee=tmp(ln+2*ip-3)+eps*tmp(ln+2*ip-2)&
                             *(tmp(ln+2*ip-1)-tmp(ln+2*ip-3))/delta
                           if (ee.lt.econst) econst=ee
                        endif
                     enddo
                  endif
               endif
            endif
            il=l-1
            m=l

         !--tabulated subsection
         else
            call tab2io(nin,0,0,tmp(l),nb,nw)
            ne=nint(tmp(l+5))
            nnt=l+5
            nbt=nint(tmp(nnt+1))
            int=nint(tmp(nnt+2))
            if (int.eq.2) int=22 !force unit base interpolation if lin-lin
            l=l+nw
            iraw=l+ne*(ng+1)
            m=iraw
            nne=0

            !--read and average spectrum for each incident energy
            do while (nne.lt.ne)
               m1=m
               call tab1io(nin,0,0,tmp(m),nb,nw)
               m=m+nw
               do while (nb.ne.0)
                  if (m.gt.ntmp) call error('getsed',&
                  'storage for tmp exceeded2',' ')
                  call moreio(nin,0,0,tmp(m),nb,nw)
                  m=m+nw
               enddo
               ! extend lowest delayed bin using sqrt(e) shape
               if (ismooth.gt.0.and.mtd.eq.455.and.&
                 nint(tmp(m1+7)).eq.1) then
                  ex=40
                  fx=.8409
                  write(nsyso,'('' extending lowest delayed bin'',&
                   &'' using sqrt(E)'')')
                  mm=nint(tmp(m1+5))
                  do while (tmp(m1+10).gt.ex)
                     if (m1+9+2*mm.gt.ntmp) call error('getsed',&
                     'storage for tmp exceeded4',' ')
                     do ix=2*mm,1,-1
                        tmp(m1+9+ix)=tmp(m1+7+ix)
                     enddo
                     tmp(m1+10)=fx*tmp(m1+12)
                     val=tmp(m1+9)
                     tmp(m1+9)=sqrt(fx)*val
                     tmp(m1+9)=sigfig(tmp(m1+9),7,0)
                     tmp(m1+11)=(1-fx*sqrt(fx))*val/(1-fx)
                     tmp(m1+11)=sigfig(tmp(m1+11),7,0)
                     mm=mm+1
                  enddo
                  tmp(m1+5)=mm
                  tmp(m1+6)=mm
                  m=m1+8+2*mm
               endif
               nne=nne+1
               tmp(l)=tmp(m1+1)
               l=l+1
               ip=2
               ir=1
               do ig=1,ng
                  e1=eg(ig)
                  if (ig.eq.1) e1=0
                  e2=eg(ig+1)
                  if (ig.eq.ng) e2=ebig
                  call intega(tmp(l),e1,e2,tmp(m1),ip,ir)
                  l=l+1
               enddo
            enddo

            !--determine econst
            l=loc(ik)+1
            nr=nint(tmp(l+4))
            np=nint(tmp(l+5))
            l=l+6+2*nr+2*np
            nr=nint(tmp(l+4))
            ne=nint(tmp(l+5))
            ir=1
            nbt=nint(tmp(l+4+2*ir))
            int=nint(tmp(l+5+2*ir))
            if (int.eq.2) int=22 !force unit base interpolation if lin-lin
            if (int.ge.11.and.int.le.15) call mess('getsed',&
              'corresponding point interpolation not available',' ')
            lnow=l+6+2*nr
            do i=2,ne
               if (i.gt.nbt) then
                  ir=ir+1
                  nbt=nint(tmp(l+4+2*ir))
                  int=nint(tmp(l+5+2*ir))
                  if (int.eq.2) int=22 !force unit base interpolation if lin-lin
               endif
               llo=lnow+(i-2)*(ng+1)
               elo=tmp(llo)
               lhi=lnow+(i-1)*(ng+1)
               ehi=tmp(lhi)
               do ig=1,ng
                  delta=abs(tmp(lhi+ig)-tmp(llo+ig))
                  if (delta.gt.eps*tmp(llo+ig)) then
                     if (int.eq.1.and.ehi.lt.econst) econst=ehi
                     if (int.ne.1) then
                        ee=elo+eps*tmp(llo+ig)*(ehi-elo)/delta
                        if (ee.lt.econst) econst=ee
                     endif
                  endif
               enddo
            enddo
            l=m
         endif
      enddo

      !--initialization complete
      nc=m-1
      allocate(sedist(nc))
      do i=1,nc
         sedist(i)=tmp(i)
      enddo
      deallocate(tmp)
      nupm=0
      return
   endif

   !--normal entry.
   enext=emax
   idis=0
   do ik=1,nk
      do ig=1,ng
         sed(ik,ig)=0
      enddo
   enddo
   iout=0
   ik=0
   do while (ik.lt.nktot.and.iout.eq.0)
      ik=ik+1
      ikt=ik
      if (nk.eq.1) ikt=1
      lnow=loc(ik)+1
      lf=nint(sedist(lnow+3))

      !--interpolate for fractional probability.
      ip=2
      ir=1
      call terpa(pe,ed,eihi,idisc,sedist(lnow),ip,ir)
      if (abs(eihi-enext).lt.enext*small.and.idisc.gt.idis)&
        idis=idisc
      if (eihi.lt.enext*(1-small)) idis=idisc
      if (eihi.lt.enext*(1-small)) enext=eihi
      if (pe.gt.zero) then
         nr=nint(sedist(lnow+4))
         np=nint(sedist(lnow+5))
         mnow=lnow+6+2*nr+2*np

         !--analytic subsection.  compute sed at e.
         if (lf.gt.1) then
            ip=2
            ir=1
            do ig=1,ng
               e1=eg(ig)
               if (ig.eq.1) e1=0
               e2=eg(ig+1)
               if (ig.eq.ng) e2=ebig
               call anased(s,ed,ethi,idisc,e1,e2,sedist(lnow),ip,ir)
               s=s*pe
               sed(ikt,ig)=sed(ikt,ig)+s
            enddo
            if (abs(ethi-enext).lt.enext*small.and.idisc.gt.idis)&
               idis=idisc
            if (ethi.lt.enext*(1-small)) idis=idisc
            if (ethi.lt.enext*(1-small)) enext=ethi

         !--tabulated subsection.  interpolate for sed at e.
         else
            nne=0
            nr=nint(sedist(mnow+4))
            ne=nint(sedist(mnow+5))
            nnow=mnow+6+2*nr
            nnt=6
            nbt=nint(sedist(mnow+nnt))
            int=nint(sedist(mnow+nnt+1))
            if (int.eq.2) int=22 !force unit base interpolation if lin-lin
            ehi=0
            do while (nne.lt.ne.and.ed.gt.ehi*(1+small).and.iout.eq.0)
               nne=nne+1
               if (nne.gt.ne) then
                  iout=1
               else
                  if (nne.gt.nbt) then
                     nnt=nnt+2
                     nbt=nint(sedist(mnow+nnt))
                     int=nint(sedist(mnow+nnt+1))
                     if (int.eq.2) int=22 !force unit base interpolation
                  endif
                  llo=nnow+(ng+1)*(nne-1)
                  elo=sedist(llo)
                  lhi=nnow+(ng+1)*nne
                  ehi=sedist(lhi)
               endif
            enddo
            if (int.gt.5) then
               klo=nnow+(ng+1)*ne
               elo=sedist(klo+1)
               nrlo=nint(sedist(klo+4))
               nplo=nint(sedist(klo+5))
               xlo=sedist(klo+4+2*nrlo+2*nplo)
               khi=klo+6+2*nrlo+2*nplo
               ehi=sedist(khi+1)
               nrhi=nint(sedist(khi+4))
               nphi=nint(sedist(khi+5))
               xhi=sedist(khi+4+2*nrhi+2*nphi)
               do while (nne.lt.ne.and.ed.gt.ehi*(1+small).and.iout.eq.0)
                  klo=khi
                  elo=ehi
                  nrlo=nrhi
                  nplo=nphi
                  xlo=xhi
                  khi=klo+6+2*nrlo+2*nplo
                  ehi=sedist(khi+1)
                  nrhi=nint(sedist(khi+4))
                  nphi=nint(sedist(khi+5))
                  xhi=sedist(khi+4+2*nrhi+2*nphi)
               enddo
            endif
            if (iout.eq.0) then
               !--unit base.
               if (int.ge.21) then
                  jnt=int-20
                  call terp1(elo,xlo,ehi,xhi,ed,xend,jnt)
                  do ig=1,ng
                     e1=eg(ig)
                     if (ig.eq.1) e1=0
                     e1lo=e1*xlo/xend
                     e1hi=e1*xhi/xend
                     e2=eg(ig+1)
                     if (ig.eq.ng) e2=ebig
                     e2lo=e2*xlo/xend
                     e2hi=e2*xhi/xend
                     ip=2
                     ir=1
                     call intega(flo,e1lo,e2lo,sedist(klo),ip,ir)
                     ip=2
                     ir=1
                     call intega(fhi,e1hi,e2hi,sedist(khi),ip,ir)
                     call terp1(elo,flo,ehi,fhi,ed,fe,jnt)
                     sed(ikt,ig)=sed(ikt,ig)+pe*fe
                  enddo
               !--corresponding points.  not implemented.
               else if (int.ge.11) then
                  do ig=1,ng
                     sed(ikt,ig)=0
                  enddo
               !--cartesion.
               else
                  do ig=1,ng
                     call terp1(elo,sedist(llo+ig),ehi,sedist(lhi+ig),ed,s,int)
                     sed(ikt,ig)=sed(ikt,ig)+s*pe
                  enddo
               endif
               !--upscatter is not allowed in secondary energy
               !--spectra, except for fission and gamma-ray
               !--production. put the upscatters into the in-group.
               sigup=0
               if (mfd.ne.15.and.mtd.ne.455) then
                  if (mtd.lt.18.or.(mtd.gt.21.and.mtd.ne.38)) then
                     ig=ng
                     do while (ed.le.eg(ig)*(1+eps).and.ig.gt.1)
                        sed(ikt,ig-1)=sed(ikt,ig-1)+sed(ikt,ig)
                        sigup=sed(ikt,ig)
                        sed(ikt,ig)=0
                        ig=ig-1
                     enddo
                  endif
               endif
               test=tenth
               if (sigup.gt.test.and.nupm.lt.10) then
                  nupm=nupm+1
                  write(strng,'(''upscatter correction '',&
                    &1p,e12.4)') sigup
                  call mess('getsed',strng,' ')
                  if (nupm.eq.10) call mess('getsed',&
                    'additional messages supressed',' ')
               endif
               if (int.eq.1.and.ehi.lt.enext) idisc=1
               if (ehi.lt.enext*(1-small)) enext=ehi
            endif
         endif
      endif
   enddo

   !--return the normal results
   if (iout.eq.0) then
      return

   !--return zeroes outside range of table
   else
      if (int.ne.1) then
         do ig=1,ng
            sed(1,ig)=0
         enddo
         enext=ehi
         if (nne.eq.ne) enext=emax
         idis=1
      else
         do ig=1,ng
            sed(1,ig)=ehi
         enddo
         enext=emax
      endif
   endif
   return
   end subroutine getsed

   subroutine anased(g,e,enext,idis,ep1,ep2,a,ip,ir)
   !-------------------------------------------------------------------
   ! Analytic secondary energy distributions.  Compute the integral
   ! between sink energies epl and eph for source energy e.  The
   ! parameters are in a in packed endf format.  Laws 5, 7, 9,
   ! 11, and 12 only.
   ! Here, gami is the incomplete gamma function (law 12 only).
   ! Here, e1 is the first order exponential integral function
   ! (law 12 only).
   !-------------------------------------------------------------------
   use physics ! provides pi
   use util    ! provides error
   use endf    ! provides terpa,intega
   use mathm   ! provides e1,gami,erfc
   ! externals
   integer::idis,ip,ir
   real(kr)::g,e,enext,ep1,ep2,a(*)
   ! internals
   integer::new,nr,np,ip2,ir2,i,idisc,lf,loct,locg,locb
   real(kr)::rp4,epl,eph,u,de,xone,xtwo,theta,xc,r1,r2,temp1
   real(kr)::xlo,xhi,rc,r4,expa,bot,rl,rh,r3,expb,expc,b,en
   real(kr)::top,ca,ef,alpha,sa,sb,aa,bb,ab,fact,ap,bp,ans
   real(kr)::h(2)
   real(kr),parameter::xmin=1.e-3_kr
   real(kr),parameter::thrhaf=1.5e0_kr
   real(kr),parameter::fivhaf=2.5e0_kr
   real(kr),parameter::brk=.01e0_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::etest=1.1e-5_kr
   real(kr),parameter::zero=0
   save new,theta,xc,rc,bot,ca,loct

   rp4=sqrt(pi/4)
   epl=ep1
   eph=ep2
   new=0
   if (ep1.lt.etest) new=1

   ! check limit on integration.
   enext=emax
   idis=0
   g=0
   lf=nint(a(4))
   if (lf.ne.12) then
      u=a(1)
      de=e-u
      if (epl.gt.de*(1+small).or.de.le.zero) return
   endif

   !--retrieve theta by interpolation.
   if (new.ne.0) then
      nr=nint(a(5))
      np=nint(a(6))
      ip2=2
      ir2=1
      loct=6+2*nr+2*np
      call terpa(temp1,e,en,idisc,a(loct+1),ip2,ir2)
      theta=temp1
      if (lf.ne.12) xc=de/theta
      if (en.lt.enext*(1-small)) idis=idisc
      if (en.lt.enext*(1-small)) enext=en
   endif

   !--compute common quantities and integration limits.
   if (lf.ne.12) then
      xlo=epl/theta
      if (eph.le.de*(1+small)) xhi=eph/theta
      if (eph.gt.de*(1+small)) xhi=xc
   endif

   !--complete processing according to law.

   !--law 5.  general evaporation spectrum.
   if (lf.eq.5) then
      nr=nint(a(loct+5))
      np=nint(a(loct+6))
      locg=loct+6+2*nr+2*np
      xone=xlo
      xtwo=xhi
      call intega(g,xone,xtwo,a(locg+1),ip,ir)
      g=g*theta

   !--law 7.  simple fission spectrum.
   else if (lf.eq.7) then
      if (new.ne.0) then
         rc=sqrt(xc)
         r4=-xc
         bot=rp4*(1-erfc(rc))-rc*exp(r4)
      endif
      rl=sqrt(xlo)
      rh=sqrt(xhi)
      r3=-xlo
      r4=-xhi
      top=rl*exp(r3)-rh*exp(r4)-rp4*(erfc(rh)-erfc(rl))
      g=top/bot

   !--law 9.  evaporation spectrum.
   else if (lf.eq.9) then
      rc=xc
      if (new.ne.0) then
         if (xc.lt.xmin) bot=rc*rc/2
         if (xc.ge.xmin) bot=1-exp(-rc)*(1+rc)
      endif
      r1=xlo
      r2=xhi
      if (xhi.lt.xmin) top=(r2*r2-r1*r1)/2
      if (xhi.ge.xmin) top=exp(-r1)*(r1+1)-exp(-r2)*(r2+1)
      g=top/bot

   !--law 11. energy dependent watt spectrum.
   else if (lf.eq.11) then
      if (new.ne.0) then
         nr=nint(a(loct+5))
         np=nint(a(loct+6))
         ip=2
         ir=1
         locb=loct+6+2*nr+2*np
         call terpa(b,e,en,idisc,a(locb+1),ip,ir)
         rc=theta*b/4
         rc=sqrt(rc)
         ca=4*rc*exp(-rc*rc)/(3*rp4)
      endif
      rh=sqrt(xc)
      r1=-rc
      r2=rh-rc
      r3=rc
      r4=rh+rc
      call hnab(h,r1,r2)
      bot=h(2)+rc*h(1)
      call hnab(h,r3,r4)
      bot=bot-h(2)+rc*h(1)
      rl=sqrt(xlo)
      rh=sqrt(xhi)
      if (rh.ge.brk) then
         r1=rl-rc
         r2=rh-rc
         r3=rl+rc
         r4=rh+rc
         call hnab(h,r1,r2)
         top=h(2)+rc*h(1)
         call hnab(h,r3,r4)
         top=top-h(2)+rc*h(1)
      else
         top=ca*(rh**3-rl**3)
      endif
      g=top/bot

   !--law 12.  madland-nix fission spectrum
   else if (lf.eq.12) then
      ef=a(loct+1)
      i=0
      alpha=sqrt(theta)
      sa=sqrt(epl)
      sb=sqrt(eph)
      do while (i.lt.2)
         beta=sqrt(ef)
         aa=(sa+beta)*(sa+beta)/theta
         bb=(sb+beta)*(sb+beta)/theta
         ab=alpha*beta
         fact=1/(3*ab)
         ap=(sa-beta)*(sa-beta)/theta
         bp=(sb-beta)*(sb-beta)/theta
         ! region i
         if (epl.ge.ef*(1+small).and.eph.gt.ef*(1+small)) then
            ans=fact*(((4*theta*bb**fivhaf/10-ab*bb*bb/2)*e1(bb)-&
              (4*theta*aa**fivhaf/10-ab*aa*aa/2)*e1(aa))-&
              ((4*theta*bp**fivhaf/10+ab*bp*bp/2)*e1(bp)-&
                (4*theta*ap**fivhaf/10+ab*ap*ap/2)*e1(ap))+&
              ((theta*bb-2*ab*sqrt(bb))*gami(thrhaf,bb)-&
                (theta*aa-2*ab*sqrt(aa))*gami(thrhaf,aa))-&
              ((theta*bp+2*ab*sqrt(bp))*gami(thrhaf,bp)-&
                (theta*ap+2*ab*sqrt(ap))*gami(thrhaf,ap))-&
              (6*theta*(gami(fivhaf,bb)-gami(fivhaf,aa)&
                -gami(fivhaf,bp)+gami(fivhaf,ap))/10)&
              -(thrhaf*ab*(exp(-bb)*(1+bb)-exp(-aa)*(1+aa)+&
                exp(-bp)*(1+bp)-exp(-ap)*(1+ap))))
         ! region ii
         else if (epl.lt.ef*(1-small).and.eph.le.ef*(1+small)) then
            ans=fact*(((4*theta*bb**fivhaf/10-ab*bb*bb/2)*e1(bb)-&
              (4*theta*aa**fivhaf/10-ab*aa*aa/2)*e1(aa))-&
              ((4*theta*bp**fivhaf/10-ab*bp*bp/2)*e1(bp)-&
              (4*theta*ap**fivhaf/10-ab*ap*ap/2)*e1(ap))+&
              ((theta*bb-2*ab*sqrt(bb))*gami(thrhaf,bb)-&
               (theta*aa-2*ab*sqrt(aa))*gami(thrhaf,aa))-&
              ((theta*bp-2*ab*sqrt(bp))*gami(thrhaf,bp)-&
               (theta*ap-2*ab*sqrt(ap))*gami(thrhaf,ap))-&
               (6*theta*(gami(fivhaf,bb)-gami(fivhaf,aa)&
                -gami(fivhaf,bp)+gami(fivhaf,ap))/10)&
               -(thrhaf*ab*(exp(-bb)*(1+bb)-exp(-aa)*(1+aa)-&
               exp(-bp)*(1+bp)+exp(-ap)*(1+ap))))
         ! region iii
         else if (epl.lt.ef*(1-small).and.eph.gt.ef*(1+small)) then
            ans=fact*(((4*theta*bb**fivhaf/10-ab*bb*bb/2)*e1(bb)-&
               (4*theta*aa**fivhaf/10-ab*aa*aa/2)*e1(aa))-&
             ((4*theta*bp**fivhaf/10+ab*bp*bp/2)*e1(bp)-&
               (4*theta*ap**fivhaf/10-ab*ap*ap/2)*e1(ap))+&
             ((theta*bb-2*ab*sqrt(bb))*gami(thrhaf,bb)-&
               (theta*aa-2*ab*sqrt(aa))*gami(thrhaf,aa))-&
             ((theta*bp+2*ab*sqrt(bp))*gami(thrhaf,bp)-&
               (theta*ap-2*ab*sqrt(ap))*gami(thrhaf,ap))-&
             (6*theta*(gami(fivhaf,bb)-gami(fivhaf,aa)-&
               gami(fivhaf,bp)+gami(fivhaf,ap))/10)&
             -(thrhaf*ab*(exp(-bb)*(1+bb)-exp(-aa)*(1+aa)&
              +exp(-bp)*(1+bp)+exp(-ap)*(1+ap)-2)))
         endif
         g=g+ans
         i=i+1
         if (i.lt.2) then
            ef=a(loct+2)
         endif
      enddo
      g=g/2

   !--illegal law
   else
      call error('anased','illegal lf.',' ')
   endif
   return

   end subroutine anased

   subroutine hnab(hh,aa,bb)
   !-------------------------------------------------------------------
   ! Compute integral from a to b of
   ! (1/sqrt(pi))*(u**n)*exp(-u**2)
   ! for n=0 and n=1.
   ! For b-a large, a difference is used.
   ! For b-a small, a direct Taylor expansion of the integral is used.
   ! Change toler to control the choice of method.
   !-------------------------------------------------------------------
   use physics ! provides pi
   use mathm   ! provides erfc
   ! externals
   real(kr)::hh(2),aa,bb
   ! internals
   integer::k,kd,idone,m,n,j,mflag,kstar,kdstar,kk,jalpha
   real(kr)::a,b,asq,con,expa,expb,h,fact,resqpi,sgn,x,xk,s,qmn
   real(kr)::term,test,xn1,xn2,xx
   real(kr)::cm(50),cmstar(50),fa(2),fb(2)
   real(kr),parameter::toler=1.e-5_kr
   real(kr),dimension(5),parameter::pow2=(/1.4142135623731e0_kr,&
     2.0e0_kr,2.8284271247462e0_kr,4.0e0_kr,5.6568542494924e0_kr/)
   real(kr),parameter::one=1
   real(kr),parameter::epsrel=1.e-9_kr
   real(kr),parameter::epsabs=1.e-30_kr
   real(kr),parameter::zero=0
   resqpi=1/sqrt(pi)

   expa=exp(-aa*aa)
   expb=exp(-bb*bb)
   fa(1)=1-erfc(abs(aa))
   if (aa.lt.zero) fa(1)=-fa(1)
   fb(1)=1-erfc(abs(bb))
   if (bb.lt.zero) fb(1)=-fb(1)
   fa(2)=(1-expa)*resqpi
   fb(2)=(1-expb)*resqpi
   hh(1)=fb(1)-fa(1)
   if (abs(hh(1)).gt.toler*abs(fb(1)).and.&
     abs(hh(1)).gt.toler*abs(fa(1))) go to 110
   if (abs(aa-bb).ge.one) go to 110
   n=0
   go to 120
  110 continue
   hh(2)=fb(2)-fa(2)
   if (abs(hh(2)).gt.toler*abs(fb(2)).and.&
     abs(hh(2)).gt.toler*abs(fa(2))) go to 130
   if (abs(aa-bb).ge.one) go to 130
   n=1

   !--compute as a taylor expansion
  120 continue
   sgn=1
   if (bb.lt.aa) then
      a=abs(bb)
      b=abs(aa)
      sgn=-sgn
   else
      a=abs(aa)
      b=abs(bb)
   endif
   if (bb.lt.0.and.mod(n,2).ne.0) sgn=-sgn
   h=(b-a)*pow2(1)
   x=pow2(1)*a
   xx=x*x
   asq=a*a
   con=exp(-asq)*resqpi/pow2(n+1)
   mflag=0
   k=n
   kd=0
   cm(1)=1
   xk=1
   if (k.ne.0) xk=x**k
   s=h*xk
   fact=h
   idone=0
   m=1
   do while (idone.eq.0.and.m.lt.50)
   m=m+1
      fact=fact*h/m
      kstar=k
      kdstar=kd+1
      do j=1,kdstar
         cmstar(j)=cm(j)
      enddo
      k=n-m+1
      if (k.lt.0) then
         kk=mod(k,2)
         k=0
         if (kk.ne.0) k=1
      endif
      kd=(n+m-1-k)/2
      cm(kd+1)=-cmstar(kdstar)
      qmn=cm(kd+1)
      if (kd.ne.0) then
         do j=1,kd
            jalpha=(2*j+k-1-kstar)/2
            beta=0
            if (jalpha.ne.0) beta=cmstar(jalpha)
            cm(j)=(2*j+k-1)*cmstar(jalpha+1)-beta
         enddo
         do j=1,kd
            qmn=qmn*xx+cm(kd+1-j)
         enddo
      endif
      xk=1
      if (k.ne.0) xk=x**k
      term=fact*xk*qmn
      s=s+term
      xn1=n+1
      xn2=abs(h*x)
      if (m.ge.max(xn1,xn2)) then
         test=epsabs+epsrel*abs(s)
         if (abs(term).le.test) then
            if (mflag.eq.1) idone=1
            mflag=1
         else
            mflag=0
         endif
      endif
   enddo
   hh(n+1)=2*con*s*sgn
   if (n.eq.0) go to 110
  130 continue
   return
   end subroutine hnab

   real(kr) function f6psp(ep,epnext,epmax,w,e,c)
   !-------------------------------------------------------------------
   ! Compute the double differential cross section for this incident
   ! energy (in the lab system), and secondary energy ep and cosine
   ! w (in the cm system) using the n-body phase-space formula.
   ! Call with ep=0 to initialize the retrieval routine for each
   ! new e.  Thereater, ep can be requested in any order.
   !-------------------------------------------------------------------
   use util ! provides error
   use endf ! provides contio
   ! externals
   real(kr)::ep,epnext,epmax,w,e,c(*)
   ! internals
   integer::npsx
   real(kr)::apsx,f1,f2,ex,eimax,cn,s
   real(kr),parameter::c3=1.2732e0_kr
   real(kr),parameter::c4=3.2813e0_kr
   real(kr),parameter::c5=5.8205e0_kr
   real(kr),parameter::thrhaf=1.5e0_kr
   real(kr),parameter::sevhaf=3.5e0_kr
   real(kr),parameter::step=.05e0_kr
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::up=1.000001e0_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::zero=0
   save cn,eimax,ex

   !--initialize
   if (ep.eq.zero) then
      apsx=c(1)
      npsx=nint(c(6))
      f1=(apsx-awp)/apsx
      f2=awr/(awr+1)
      ex=thrhaf*npsx-4
      eimax=f1*(f2*e+q)
      if (eimax.le.zero) eimax=1
      if (npsx.eq.3) then
         cn=c3/eimax**2
      else if (npsx.eq.4) then
         cn=c4/eimax**sevhaf
      else if (npsx.eq.5) then
         cn=c5/eimax**5
      else
         call error('f6psp','3, 4, or 5 particles only.',' ')
      endif
      f6psp=0
      epnext=step*eimax
      epmax=eimax

   !--compute cross section
   else
      s=0
      if (ep.lt.eimax*(1-small)) s=cn*sqrt(ep)*(eimax-ep)**ex
      f6psp=s
      epnext=ep+step*eimax
      if (epnext.gt.eimax*(1+small)) epnext=up*eimax
      if (ep.ge.eimax*(1-small)) epnext=emax
   endif
   return
   end function f6psp

end module groupm

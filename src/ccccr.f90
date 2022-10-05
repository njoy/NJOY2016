module ccccm
   ! provides subroutine ccccr for NJOY2016
   use locale
   implicit none
   private
   public ccccr

   ! equivalenced arrays for CCCC input and output
   integer,parameter:: isiza4=100000, isiza8=isiza4/2
   real(k4)::a(isiza4)      ! reals are 4-byte
   integer(k4)::ia(isiza4)  ! integers are 4-byte
   real(k8)::ha(isiza8)     ! Hollerith data are 8-byte
   character(8)::ta(isiza8) ! text equivalent to Hollerith
   integer::mult=2         ! used for counting 8-byte entries
   equivalence(a(1),ia(1),ha(1),ta(1))

   ! i/o units
   integer::nin,nisot,nbrks,ndlay,nscrt1,nscrt2,nscrt3

   ! common variables
   real(8)::huse(2)
   character(8)::tuse(2)
   equivalence(huse(1),tuse(1))
   integer::ivers

   ! cccc user and version values

   ! array for reading gendf data
   integer,parameter::maxe=8000
   real(kr),dimension(:),allocatable::e

   ! hollerith set id record
   character(6)::hsetid(12)

   ! common parameters
   integer::ngps,nggrup,niso,maxord,ifopt,nsblk,nscmax,ichix

   ! isotxs material identity and classification arrays
   integer::imat(200),kbr(200)
   character(6)::hisonm(200),habsid(200),hident(200),hmat(200)

   ! isotxs isotope data
   real(kr)::amass(200),efiss(200),ecapt(200),temp(200),&
     sigpot(200),adens(200)

   integer::iflags(15)
   integer::ltot,ltrn

   ! isotxs matrix parameters
   integer::irsize,ng2z,jlz,ieof,lord1,nsiza

   ! user input for brkoxs file
   integer::nti,nzi,nreact,nzt,ntj,nzj

   ! temperature and sigmazero areas for brkoxs data
   real(kr)::atem(20),ctem(20),ntat(20),item(20)
   real(kr)::asig(20),csig(20),ntap(20),isig(20)
   real(kr)::xspo,xspot(200)
   real(kr)::tzro

   ! brkoxs parameters
   integer::next,isopec,l1,l2,l13
   integer::jbl,jbh
   integer::mt2tem(10),ntfl
   integer::matd,nsb,lrsize

   ! tolerance for self-shielding factors in brkoxs data
   real(kr),parameter::sstol=.00001e0_kr

   ! delayed parameters
   integer::nisod,nfam

contains

   subroutine ccccr
   !--------------------------------------------------------------------
   !
   !  Produce CCCC-IV files from njoy intermediate cross section
   !  library.
   !
   !  Working from a groupr output tape, this module produces
   !  the following three standard interface files,
   !
   !           ISOTXS       BRKOXS       DLAYXS,
   !
   !  as specified by the committee for computer code coordination
   !  (CCCC), to facilitate the exchange of nuclear data for reactor
   !  calculations (Reference 1).
   !      In a given run, all three files can be produced using the
   !  same user-specified list of isotopes.  The code will ignore
   !  isotopes which are not present on the groupr tape (and in the
   !  case of DLAYXS, isotopes without delayed neutron data).
   !      The ISOTXS coding allows for NSBLK equal to one or ngroup.
   !  In addition, files with higher order matrices can be produced
   !  with a separate block for each l-order (ifopt=2) or with all
   !  orders in one block (ifopt=1).  This flexibility accommodates
   !  large group structures.  Fission vectors or fission
   !  matrices can be produced.
   !      In BRKOXS, the potential scattering cross section for all
   !  energy groups is equal to the user-input value (xspo).
   !  The elastic removal f-factor is supplied as the sixth reaction.
   !
   !  1. R.D.Odell. Standard Interface Files and Procedures
   !     for Reactor Physics Codes, Version IV,
   !     LANL report LA-6941-MS (Sept.77)
   !
   !
   !---input specifications (free format)---------------------------
   !
   !-ccccr-
   ! card 1 units
   !    nin      input unit for data from groupr
   !    nisot    output unit for isotxs (0 if isotxs not wanted)
   !    nbrks    output unit for brkoxs (0 if brkoxs not wanted)
   !    ndlay    output unit for dlayxs (0 if dlayxs not wanted)
   ! card 2 identification
   !    lprint    print flag (0/1=not print/printed)
   !    ivers     file version number (default=0)
   !    huse      user identification (12 characters)
   !              delimited by *, ended by /.
   !              (default=blank)
   ! card 3
   !    hsetid    hollerith identification of set (12 characters)
   !              delimited by *, ended by /.
   !              (default=blank)
   ! card 4 file control
   !    ngroup    number of neutron energy groups
   !    nggrup    number of gamma energy groups
   !    niso      number of isotopes desire
   !    maxord    maximum legendre order
   !    ifopt     matrix blocking option (1/2=blocking by
   !                               reaction/legendre order)
   ! card 5 isotope parameters (one card per isotope)
   !       (first four words are hollerith, up to six characters
   !       each, delimited by *)
   !    hisnm     hollerith isotope label
   !    habsid    hollerith absolute isotope label
   !    hident    identifier of data source library (endf)
   !    hmat      isotope identification
   !    imat      numerical isotope identifier (endf mat number)
   !    xspo      average potential scattering cross sect. (brkoxs)
   !
   !-cisotx- (only if nisot.gt.0)
   ! card 1 file control
   !    nsblok    subblocking option for scattering matrix
   !              (1 or ngrup sub-blocks allowed)
   !    maxup     maximum number of upscatter groups (always zero)
   !    maxdn     maximum number of downscatter groups
   !    ichix     fission chi representation
   !                   -1   vector (using groupr flux)
   !                    0   none
   !                   +1   vector (using input flux)
   !                .gt.1   matrix
   ! card 2 chi vector control (ichix=1 only)
   !    spec      ngroup flux values used to collapse the groupr
   !              fission matrix into a chi vector
   ! card 3 chi matrix control (ichix.gt.1 only)
   !    spec      ngroup values of spec(i)=k define the range of
   !              groups i to be averaged to obtain spectrum k.
   !              index k ranges from 1 to ichi.
   !              the model flux is used to weight each group i.
   ! card 4 isotope control (one card per isotope)
   !    kbr       isotope classification
   !    amass     gram atomic weight
   !    efiss     total thermal energy/fission
   !    ecapt     total thermal energy/capture
   !    temp      isotope temperature
   !    sigpot     average effective potential scattering
   !    adens     density of isotope in mixture
   !
   !-cbrkxs- (only if nbrks.gt.0)
   ! card 1 (2i6) file data
   !    nti       number of temperatures desired
   !              (-n means accept first n temperatures)
   !    nzi       number of sigpo values desire
   !              (-n means accept first n dilution factors)
   ! card 2 (not needed if nti.lt.0)
   !    atem(nti) values of desired temperatures
   ! card 3 (not needed if nzi.lt.0)
   !    asig(nzi) values of desired sigpo
   !
   !-cdlayx-- no input required
   !
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use util ! provides timer,openz,repoz,closz
   ! internals
   integer::lprint,i
   real(kr)::time
   character(72)::text
   character(12)::word

   !--initialize.
   call timer(time)
   write(nsyso,'(/'' ccccr...'',&
     &''produce cccc format output files'',28x,f8.1,''s'')')&
     time
   write(nsyse,'(/'' ccccr...'',60x,f8.1,''s'')') time
   allocate(e(maxe))
   next=1
   nscrt1=10
   nscrt2=11
   nscrt3=12
   call openz(-nscrt2,1)
   call openz(-nscrt3,1)

   !--user input for units to be used.
   nisot=0
   nbrks=0
   ndlay=0
   read(nsysi,*) nin,nisot,nbrks,ndlay
   write(nsyso,'(/&
     &'' input gendf unit ..................... '',i10/&
     &'' output isotxs unit ................... '',i10/&
     &'' output brkoxs unit ................... '',i10/&
     &'' output dlayxs unit ................... '',i10)')&
     nin,nisot,nbrks,ndlay
   if (nin.lt.0) nscrt1=-nscrt1
   call openz(nin,0)
   call openz(-nisot,1)
   call openz(-nbrks,1)
   call openz(-ndlay,1)
   call openz(nscrt1,1)

   !--user input for file identification
   ivers=0
   word=' '
   read(nsysi,*) lprint,ivers,word
   read(word,'(2a6)') tuse(1),tuse(2)

   !--read other user input
   text=' '
   read(nsysi,*) text
   read(text,'(12a6)') (hsetid(i),i=1,12)
   read(nsysi,*) ngps,nggrup,niso,maxord,ifopt
   do i=1,niso
      read(nsysi,*) hisonm(i),habsid(i),hident(i),&
        hmat(i),imat(i),xspot(i)
   enddo

   !--prepare isotxs file if desired
   if (nisot.ne.0) then
      call cisotx
      if (lprint.eq.1) call pisotx(nisot)
   endif

   !--prepare brkoxs file if desired
   if (nbrks.ne.0) then
      call cbrkxs
      if (lprint.eq.1) call pbrkxs(nbrks)
   endif

   !--prepare dlayxs file if desired
   if (ndlay.ne.0) then
      call cdlyxs(ndlay,nin)
      if (lprint.eq.1.and.nisod.ne.0) call pdlyxs(ndlay)
   endif

   !--finished
   deallocate(e)
   call repoz(-nisot)
   call repoz(-nbrks)
   call repoz(-ndlay)
   call repoz(nin)
   call closz(nin)
   call closz(-nisot)
   call closz(-nbrks)
   call closz(-ndlay)
   call closz(nscrt1)
   call closz(-nscrt2)
   call closz(-nscrt3)
   call timer(time)
   write(nsyso,'(69x,f8.1,''s''/1x,&
     &7(''**********''),''*******'')') time
   return
   end subroutine ccccr

   subroutine cisotx
   !--------------------------------------------------------------------
   ! Write CCCC-IV ISOTXS interface file.  The n2n matrix is production
   ! based for version IV.  To get reaction based matrix as in
   ! version III, set n2niv to 0.
   !--------------------------------------------------------------------
   use util ! provides repoz
   ! internals
   integer::n4,irec,nwds,l1h,ngroup,ichist,nsblok,irec2,irec3
   integer::jlord,jband,jtemp,iso,lnext,ichi,ifis,ialf,inp,in2n,ind,int
   integer::ltot,ltrn,istrpd,inc,n,lordn,lord,kmax,k,ji,jl,ju,jj
   integer::l,ndex,mdex,j,n1,ibw,jm,isnk,isrce,jk
   integer::l1r,l1i,l13h,nwh,l13r,nwr,l13i,nwi,i
   character(6)::hname='isotxs'
   integer::n2niv=1

   n4=nisot

   !--read user input
   call ruinis

   !--create and write isotxs data sets on scratch tape
   call isxdat

   !--file identification
   irec=1
   nwds=3*mult+1
   l1h=(l1-1)/mult+1
   read(hname,'(a6)') ta(l1h)
   ha(l1h+1)=huse(1)
   ha(l1h+2)=huse(2)
   ia(3*mult+l1)=ivers
   call repoz(-n4)
   write(n4)(ha(l1h+i-1),i=1,3),ia(l1+3*mult)

   !--file control
   irec=2
   l2=l1+nwds
   nwds=8
   write(n4)(ia(l2+i-1),i=1,nwds)
   ngroup=ia(l2)
   niso=ia(l2+1)
   ichist=ia(l2+5)
   nscmax=ia(l2+6)
   nsblok=ia(l2+7)

   !--file data
   irec=3
   l13h=(l13-1)/mult
   nwh=12+niso
   l13r=l13-1+nwh*mult
   nwr=ngroup*(2+ichist*(2/(ichist+1)))+1
   l13i=l13r+nwr
   nwi=niso
   nwds=nwh*mult+nwr+nwi
   write(n4)(ha(l13h+i),i=1,nwh),(a(l13r+i),i=1,nwr),&
     (ia(l13i+i),i=1,nwi)

   !--set chi data
   ! no set chi data provided

   !--loop over isotopes
   irec2=0
   irec3=0
   call repoz(-nscrt3)
   call repoz(-nscrt2)
   jlord=l1
   jband=jlord+2*nscmax
   jtemp=jband+nscmax*ngroup
   do 110 iso=1,niso

   !--isotope control and group independent data
   irec=irec+1
   l1h=(l1-1)/mult+1
   nwh=3
   l1r=l1+mult*nwh
   nwr=6
   l1i=l1r+nwr
   nwi=11+2*nscmax+2*ngroup*nscmax
   nwds=mult*nwh+nwr+nwi
   irec3=irec3+1
   read(nscrt3)(ha(l1h+i-1),i=1,nwh),(a(l1r-1+i),i=1,nwr),&
     (ia(l1i+i-1),i=1,nwi)
   write(n4)(ha(l1h+i-1),i=1,nwh),(a(l1r-1+i),i=1,nwr),&
     (ia(l1i+i-1),i=1,nwi)
   lnext=l1+3*mult
   ichi=ia(lnext+7)
   ifis=ia(lnext+8)
   ialf=ia(lnext+9)
   inp=ia(lnext+10)
   in2n=ia(lnext+11)
   ind=ia(lnext+12)
   int=ia(lnext+13)
   ltot=ia(lnext+14)
   ltrn=ia(lnext+15)
   istrpd=ia(lnext+16)
   call stow(a(lnext+17),a(jlord),2*nscmax)
   call stow(a(lnext+17+2*nscmax),a(jband),nscmax*ngroup)

   !--principal cross sections
   irec=irec+1
   nwr=(1+ltrn+ltot+ialf+inp+in2n+ind+int+istrpd+2*ifis+&
     ichi*(2/(ichi+1)))*ngroup
   nwds=nwr
   irec2=irec2+1
   read(nscrt2)(a(jtemp+i-1),i=1,nwr)
   write(n4)(a(jtemp+i-1),i=1,nwr)

   !--isotope chi data
   if (ichi.gt.1) then
      nwr=ngroup*ichi
      nwi=ngroup
      nwds=nwr+nwi
      irec2=irec2+1
      read(nscrt2)(a(jtemp+i-1),i=1,nwr),(ia(jtemp+nwr+i-1),i=1,nwi)
      irec=irec+1
      write(n4)(a(jtemp+i-1),i=1,nwr),(ia(jtemp+nwr+i-1),i=1,nwi)
   endif

   !--scattering data
   inc=(ngroup-1)/nsblok+1
   do n=1,nscmax
      lordn=ia(jlord+nscmax+n-1)
      if (lordn.ne.0) then
         lord=maxord+1
         if (ifopt.eq.2) lord=1
         kmax=0
         k=jlord-1+2*nscmax+(n-1)*ngroup
         do ji=1,nsblok
            jl=(ji-1)*inc+1
            ju=ji*inc
            nwds=0
            do jj=jl,ju
               l=k+jj
               nwds=nwds+ia(l)
                  kmax=kmax+ia(l)
            enddo
            nwds=nwds*lord
            if (nwds.gt.0) then
               irec2=irec2+1
               read(nscrt2)(a(jtemp+i-1),i=1,nwds)
               if (lordn.ne.0) then
                  ndex=jtemp-1
                  mdex=jtemp-1
                  do j=jl,ju
                     n1=jlord-1+2*nscmax+(n-1)*ngroup+j
                     ibw=ia(n1)
                     if (ibw.ne.0) then
                        do jm=1,lordn
                           isnk=ndex+(jm-1)*ibw
                           isrce=mdex+(jm-1)*ibw
                           do jk=1,ibw
                              a(jk+isnk)=a(jk+isrce)
                           enddo
                           isnk=isnk+ibw
                        enddo
                        ndex=isnk
                     endif
                     mdex=mdex+ibw*lord
                  enddo
                  nwds=ndex-jtemp+1
                  if (nwds.gt.0) then
                     if (ia(jlord+n-1).ge.300.and.n2niv.gt.1) then
                        do j=1,nwds
                           jj=jtemp-1+j
                           a(jj)=a(jj)/2
                        enddo
                     endif
                     irec=irec+1
                     write(n4)(a(jtemp+i-1),i=1,nwds)
                  endif
               endif
            endif
         enddo
      endif
   enddo

   !--close isotope loop
  110 continue
   next=isopec
   return
   end subroutine cisotx

   subroutine ruinis
   !--------------------------------------------------------------------
   ! Reads user input pertinent to ISOTXS file
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   ! internals
   integer::maxup,maxdn,ichist,nwds,i,ndex,j

   !--isotope independent data
   read(nsysi,*) nsblk,maxup,maxdn,ichix
   ichist=0
   nscmax=4+4*maxord*(ifopt-1)
   nwds=ngps
   if (2*(nwds/2).ne.nwds) nwds=nwds+1
   isopec=next
   next=isopec+nwds
   if (ichix.gt.0) then
      read(nsysi,*) (a(isopec+i-1),i=1,ngps)
   endif
   nwds=3*mult+8+mult
   l1=next
   next=l1+nwds
   nwds=12*mult
   l13=next
   next=l13+nwds
   ndex=l1+3*mult
   ia(ndex+1)=ngps
   ia(ndex+2)=niso
   ia(ndex+3)=maxup
   ia(ndex+4)=maxdn
   ia(ndex+5)=maxord
   ia(ndex+6)=ichist
   ia(ndex+7)=nscmax
   ia(ndex+8)=nsblk
   j=(l13-1)/mult
   do i=1,12
      read(hsetid(i),'(a6)') ta(i+j)
   enddo

   !--isotope dependent data
   do i=1,niso
      read(nsysi,*) kbr(i),amass(i),efiss(i),ecapt(i),&
        temp(i),sigpot(i),adens(i)
   enddo
   return
   end subroutine ruinis

   subroutine isxdat
   !--------------------------------------------------------------------
   ! Processes all ISOTXS data and writes on scratch tapes.
   ! Data stored in a(i) as follows:
   !
   !    location                 variable                    length
   !    --------                 --------                    ------
   !      l1                    *isotxs*                       1
   !      l2                     huse(1)                       1
   !      l3                     huse(2)                       1
   !      l4                      ivers                        1
   !      l5                     ngroup                        1
   !      l6                     niso                          1
   !      l7                     maxup                         1
   !      l8                     maxdn                         1
   !      l9                     maxord                        1
   !      l10                    ichist                        1
   !      l11                    nscmax                        1
   !      l12                    nsblok                        1
   !      l13                    hsetid                       12
   !      l14                    hisonm                     niso
   !
   !      l16                     vel                     ngroup
   !      l17                    emax                     ngroup
   !      l18                    emin                          1
   !      l19                    loca                       niso
   !      l20            isotope dependent data               -
   !      l21            idsct,lord,jband,ijj                 -
   !      l22      principal x-sections or scatt. matrix      -
   !
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use util ! provides repoz,error
   use endf ! provides endf routines and variables
   ! internals
   integer::nwds,l14,l16,l17,nwords,l18,l19,l20,l21,l21a,l22
   integer::nb,nw,nz,nt,ntw,loc,iqd,nenpts,i,jiso,loca,irec2,irec3
   integer::j,iin,jj,nsz,ichid,l20m,loci,k,nl,lk,l,ijz,nlgt,isol
   integer::nwh,nwr,nwi,llh,llr,lli
   integer::nrct(4)
   integer::imusd(200)

   !--assign storage
   nwds=mult*niso
   l14=next
   next=l14+nwds
   l16=next
   next=l16+ngps
   l17=next
   next=l17+ngps
   nwords=1
   l18=next
   next=l18+nwords
   nwds=niso+1
   if (((l18+nwds)/mult)*mult.eq.(l18+niso)) nwds=niso
   l19=next
   next=l19+nwds
   nwds=3*mult+17
   l20=next
   next=l20+nwds
   nwds=2*nscmax*(1+ngps)
   l21=next
   next=l21+nwds
   nwds=ngps
   l21a=next
   next=l21a+nwds
   l22=next

   !--read header for first material
   call repoz(nin)
   call tpidio(nin,0,0,e(1),nb,nw)
  100 continue
   call contio(nin,0,0,e(1),nb,nw)
   if (mfh.eq.0.or.mth.eq.0) go to 100
   nz=l2h
   nt=n2h
   ntw=nt
   call listio(nin,0,0,e(1),nb,nw)
   ngps=l1h
   loc=l1+1+3*mult
   if (ngps.ne.ia(loc))&
     call error('isxdat','incompatible group structures.',' ')
   loc=1+nw
   do while (nb.ne.0)
      call moreio(nin,0,0,e(loc),nb,nw)
      loc=loc+nw
   enddo
   iqd=1+6+ntw+nz+ngps
   ! rearrange group bounds, store at l17
   nenpts=ngps+1
   do i=1,nenpts
      a(i-1+l17)=e(-i+1+iqd)
   enddo

   !--begin isotope loop
   jiso=0
   loca=0
   ia(l19)=0
   do i=1,200
      imusd(i)=0
   enddo
   call repoz(nin)
   irec2=0
   irec3=0
   call repoz(-nscrt2)
   call repoz(-nscrt3)
   do 300 i=1,niso
   call tpidio(nin,0,0,e(1),nb,nw)
  120 continue
   call contio(nin,0,0,e(1),nb,nwds)
   if (math.eq.-1) go to 350
   if (math.eq.0) go to 120
   if (mfh.eq.0.or.mth.eq.0) go to 120
   if (mfh.gt.1) go to 130
   call tosend(nin,0,0,e(1))
   go to 120
  130 continue
   do 140 j=1,niso
   if (math.ne.imat(j)) go to 140
   if (imusd(j).eq.1) go to 140
   jiso=jiso+1
   iin=j
   imusd(j)=1
   jj=(l14-1)/mult+jiso
   read(hisonm(j),'(a6)') ta(jj)
   go to 160
  140 continue
   call tomend(nin,0,0,e(1))
   go to 120

   !--copy this material to a scratch tape
  160 continue
   nsz=l2h
   call repoz(nscrt1)
   nsh=1
   call contio(0,nscrt1,0,e(1),nb,nw)
   call tomend(nin,nscrt1,0,e(1))
   call atend(nscrt1,0)
   call repoz(nscrt1)

   !--process principal cross section file
   call prinxs(a(isopec),a(l16),a(l22),ia(l22),a(l21a),i,nsz,ichid,&
     irec2)

   !--write isotope independent data
   jj=(l20-1)/mult+1
   read(habsid(iin),'(a6)') ta(jj)
   read(hident(iin),'(a6)') ta(jj+1)
   read(hmat(iin),'(a6)') ta(jj+2)
   l20m=l20+3*mult-1
   a(l20m+1)=amass(iin)
   a(l20m+2)=efiss(iin)
   a(l20m+3)=ecapt(iin)
   a(l20m+4)=temp(iin)
   a(l20m+5)=sigpot(iin)
   a(l20m+6)=adens(iin)
   ia(l20m+7)=kbr(iin)
   ia(l20m+8)=ichid
   ia(l20m+9)=iflags(4)
   ia(l20m+10)=iflags(7)
   ia(l20m+11)=iflags(8)
   ia(l20m+12)=iflags(9)
   ia(l20m+13)=iflags(10)
   ia(l20m+14)=iflags(11)
   ia(l20m+15)=ltrn
   ia(l20m+16)=ltot
   ia(l20m+17)=0

   !--process scattering matrix file
   call isomtx(ia(l21),a(l22),l21,l22,loci,nsz,irec2)
   ia(l1+3*mult+8)=nsblk
   nrct(1)=200
   nrct(2)=100
   nrct(3)=300
   nrct(4)=000
   if (ifopt.ne.2) then
      do k=1,4
         ia(l21-1+k)=nrct(k)
      enddo
   else
      nl=maxord+1
      do k=1,4
         lk=l21-1+(k-1)*nl
         do l=1,nl
            ia(l+lk)=l-1+nrct(k)
         enddo
      enddo
   endif
   ijz=l21-1+nscmax*(2+ngps)
   do k=1,ngps
      lk=ijz+(k-1)*nscmax
      do l=1,nscmax
         ia(l+lk)=1
      enddo
   enddo
   ! calculate loca(i)
   if (i.ne.niso) then
      ia(l19+i)=loca+loci+2
      loca=ia(l19+i)
   endif

   !--write isotope-dependent data to nscrt3
   nlgt=l22-l20
   nwh=3
   nwr=6
   nwi=11+2*nscmax+2*ngps*nscmax
   nlgt=mult*nwh+nwr+nwi
   llh=(l20-1)/mult
   llr=l20+mult*nwh-1
   lli=llr+nwr
   irec3=irec3+1
   write(nscrt3)(ha(llh+j),j=1,nwh),(a(llr+j),j=1,nwr),&
     (ia(lli+j),j=1,nwi)

   !--all data read
  300 continue
   go to 360
  350 write(nsyso,'(/&
     &'' input tape does not contain all specified isotopes.'')')
  360 ia(l1+3*mult+2)=jiso

   !--thin out the file data (jiso not equal to niso)
   isol=niso-jiso
   do i=l16,l20
      a(i-isol)=a(i)
   enddo
   return
   end subroutine isxdat

   subroutine prinxs(spec,b,a,ia,cspc,iiso,nsz,ichid,irec2)
   !--------------------------------------------------------------------
   ! Creates principal cross section file.
   !--------------------------------------------------------------------
   use util ! provides repoz,error
   use endf ! provides endf routines and variables
   ! externals
   integer::iiso,nsz,ichid,irec2
   real(4)::spec(*),b(*),a(*),cspc(*)
   integer(4)::ia(*)
   ! internals
   integer::i,mt19,j,nptr,mchi,nwd,ig,nb,nwds,nl,nz,ng2,ig2lo
   integer::lz,loc,indx,mptr,ispec,jg2,locsg,ig2,jg,locch,loca
   integer::id,locn,locd,irct,ndex,mdex,jndx,jmdx
   real(kr)::flux,sig
   real(kr)::cnorm(20),dnorm(20)
   integer::iptr(12)
   integer,dimension(15),parameter::nsm=(/0,1,102,0,0,0,107,103,&
     16,104,105,6,7,8,9/)
   integer,dimension(11),parameter::nstx=(/0,1,102,0,0,0,107,103,&
     0,104,105/)
   integer,dimension(5),parameter::ntn=(/6,7,8,9,16/)
   real(kr),parameter::epsn=1.e-10_kr
   real(kr),parameter::zero=0

   ltot=1
   ltrn=1
   do i=1,20
      cnorm(i)=0
      dnorm(i)=0
   enddo
   do i=1,ngps
      cspc(i)=0
   enddo
   mt19=0

   !--write principal cross sections
   ! calculate position in a(i) where each mt file begins
   do i=1,12
      j=i-1
      iptr(i)=j*ngps+1
   enddo
   nptr=iptr(12)+ngps-1+ngps*nsz*(maxord+1)
   if (ichix.gt.1) mchi=nptr+1
   if (ichix.gt.1) nptr=nptr+ichix*ngps
   ! clear iflags(i),a(i)
   do i=1,nptr
      a(i)=0
   enddo
   do i=1,15
      iflags(i)=0
   enddo
   ichid=0
   ! read input data
   nwd=nsz*(maxord+1)*ngps
   call repoz(nscrt1)
   ig=ngps
  115 continue
   if (ig.lt.ngps) go to 120
   call contio(nscrt1,0,0,e(1),nb,nwds)
   nl=l1h
   nz=l2h
   if (math.eq.0) go to 240
   if (mfh.eq.0.or.mth.eq.0) go to 115
   if (mfh.gt.1) go to 120
   call tofend(nscrt1,0,0,e(1))
   go to 115
  120 continue
   call listio(nscrt1,0,0,e(1),nb,nwds)
   ng2=l1h
   ig2lo=l2h
   ig=n2h
   lz=6
   loc=1+nwds
   do while (nb.ne.0)
      if (loc+302.gt.maxe) call error('prinxs',&
        'endf input size exceeded',' ')
      call moreio(nscrt1,0,0,e(loc),nb,nwds)
      loc=loc+nwds
   enddo
   if (mfh.ne.3) go to 205

   !--process mf=3 data
   ! read average velocities.
   if (mth.ne.259) go to 125
   if (iiso.ne.1) go to 115
   indx=lz+nl*nz+1
   mptr=ngps-ig+1
   ! assume reciprocal average velocities.
   b(mptr)=sigfig(100/e(indx),7,0)
   go to 115
  125 continue
   ! check mt numbers and set flags
   do 130 i=1,15
   if (mth.ne.nsm(i)) go to 130
   iflags(i)=1
   go to 140
  130 continue
   go to 190
   ! check for (n,2n) and sum
  140 continue
   do 150 i=1,5
   if (mth.eq.ntn(i)) go to 160
  150 continue
   go to 170
  160 continue
   mptr=iptr(9)-1+ngps-ig+1+nwd
   indx=lz+nl*nz+1
   a(mptr)=a(mptr)+e(indx)
   iflags(9)=1
   go to 115
   ! check for total, capture, (n,p), (n,d), (n,t), (n,alpha)
  170 continue
   do 180 i=1,11
   if (mth.ne.nstx(i)) go to 180
   mptr=iptr(i)-1+ngps-ig+1+nwd
   indx=lz+nl*nz+1
   a(mptr)=e(indx)
   go to 190
  180 continue
   ! store parameters for transport cross section
  190 continue
   if (mth.eq.1) then
      mptr=iptr(1)-1+ngps-ig+1+nwd
      indx=lz+nl*(nz+1)
      a(mptr)=e(indx)
      mptr=iptr(12)-1+ngps-ig+1+nwd
      indx=lz+nl
      a(mptr)=e(indx)
   ! fission cross section
   else if (mth.eq.18.or.mth.eq.19.or.mth.eq.20.or.mth.eq.21&
     .or.mth.eq.38) then
      if (mth.eq.19) mt19=1
      if (mth.eq.18) then
         iflags(4)=1
         iflags(5)=1
         if (ichix.ne.0) ichid=iabs(ichix)
         iflags(6)=ichid
         mptr=iptr(4)-1+ngps-ig+1+nwd
         indx=lz+nl*nz+1
         a(mptr)=e(indx)
      endif
   ! delay nu
   else if (mth.eq.455) then
      mptr=iptr(5)-1+ngps-ig+1+nwd
      a(mptr)=a(mptr)+e(lz+2)*e(lz+3)
      if (ichid.eq.1) then
         ispec=1
         if (ichix.lt.0) flux=e(lz+1)
         if (ichix.gt.0) flux=spec(ngps-ig+1)
      else
         ispec=nint(spec(ngps-ig+1))
         flux=e(lz+1)
      endif
      dnorm(ispec)=dnorm(ispec)+flux*e(lz+2)*e(lz+3)
   endif

   !--process mf=6 data
  205 continue
   if (mfh.ne.6) go to 230
   if (mth.eq.18.or.mth.eq.19.or.mth.eq.20.or.mth.eq.21&
     .or.mth.eq.38) go to 207
   ! transport cross section
   do j=2,ng2
      jg2=ngps+1-ig2lo-j+2
      if (nl.ge.2) then
         locsg=lz+nl*nz*(j-1)+2
         mptr=iptr(1)-1+ngps+1-ig+nwd
         a(mptr)=a(mptr)-e(locsg)
      endif
   enddo
   go to 230
   ! snutot and chi
  207 continue
   if (mt19.eq.0.and.mth.ne.18) go to 230
   if (mt19.eq.1.and.mth.eq.18) go to 230
   if (ig.ne.0) then
      do j=2,ng2
         jg2=ngps+1-ig2lo-j+2
         locsg=lz+nl*nz*(j-1)+1
         if (ig2lo.gt.0) then
            !-- matrix part
            mptr=iptr(5)-1+ngps-ig+1+nwd
            a(mptr)=a(mptr)+e(locsg)
            if (ichid.eq.1) then
               ispec=1
               if (ichix.lt.0) flux=e(lz+1)
               if (ichix.gt.0) flux=spec(ngps-ig+1)
               mptr=iptr(6)-1+jg2+nwd
            else
               ispec=nint(spec(ngps-ig+1))
               flux=a(lz+1)
               mptr=mchi+(jg2-1)*ichid+ispec-1
            endif
            a(mptr)=a(mptr)+e(locsg)*flux
            cnorm(ispec)=cnorm(ispec)+e(locsg)*flux
         else
            !--production part
            do i=1,ngps
               sig=e(locsg)*cspc(i)
               mptr=iptr(5)-1+ngps-ig+1+nwd
               a(mptr)=a(mptr)+sig
               if (ichid.eq.1) then
                  ispec=1
                  if (ichix.lt.0) flux=e(lz+1)
                  if (ichix.gt.0) flux=spec(ngps-ig+1)
                  mptr=iptr(6)-1+i+nwd
               else
                  ispec=nint(spec(i-1)*ichid)+ispec-1
                  flux=e(lz+1)
                  mptr=mchi+(i-1)*ichid+ispec-1
               endif
               a(mptr)=a(mptr)+sig*flux
               cnorm(ispec)=cnorm(ispec)+sig*flux
            enddo
         endif
      enddo
   else
      !--spectrum part
      do j=1,ng2
         jg2=ngps-ig2lo-j+2
         locsg=lz+nl*nz*(j-1)+1
         cspc(jg2)=e(locsg)
      enddo
   endif

   !--process mf=5 data
   ! complete accumulation of data for snutot and chiso
  230 continue
   if (mfh.eq.5.and.mth.eq.455) then
      do ig2=2,ng2
         jg=ngps+1-ig2lo-ig2+2
         do ispec=1,ichid
            if (ichid.eq.1) locch=iptr(6)-1+jg+nwd
            if (ichid.gt.1) locch=mchi+(jg-1)*ichid+ispec-1
            loca=lz+nl*(ig2-1)
            do id=1,nl
               a(locch)=e(id+loca)*dnorm(ispec)+a(locch)
               cnorm(ispec)=e(id+loca)*dnorm(ispec)+cnorm(ispec)
            enddo
         enddo
      enddo
   endif
   go to 115

   !--complete calculations of fission nu and chi
  240 continue
   if (iflags(4).ne.0) then
      do j=1,ngps
         !--fission nu
         locn=iptr(5)-1+j+nwd
         locd=iptr(4)-1+j+nwd
         if (a(locd).le.epsn) then
            a(locn)=0
         else
            a(locn)=a(locn)/a(locd)
         endif
         !--fission chi
         do ispec=1,ichid
            if (ichid.eq.1) locn=iptr(6)-1+j+nwd
            if (ichid.gt.1) locn=mchi+(j-1)*ichid+ispec-1
            if (cnorm(ispec).ne.zero) a(locn)=a(locn)/cnorm(ispec)
         enddo
      enddo
   endif

   !--collapse principal cross section file to eliminate reactions
   !--which have zero cross section
   irct=ltot+ltrn+1
   ndex=irct*ngps+1+nwd
   mdex=ndex
   do i=4,11
      if (iflags(i).eq.1) then
         jndx=ndex-1
         jmdx=mdex-1
         do j=1,ngps
            a(j+jndx)=a(j+jmdx)
         enddo
         irct=irct+1
         ndex=ndex+ngps
      endif
      mdex=mdex+ngps
   enddo
   nwds=irct*ngps

   !--write principal cross sections on scratch file
   ndex=nwd+1
   irec2=irec2+1
   write(nscrt2)(a(ndex+i-1),i=1,nwds)
   if (ichid.le.1) return

   !--write isotope chi data
   nwds=ngps*(ichid+1)
   irec2=irec2+1
   do i=1,ngps
      ia(nptr+i)=nint(spec(i))
   enddo
   write(nscrt2)(a(mchi+i-1),i=1,nwds)
   return
   end subroutine prinxs

   subroutine isomtx(ia,b,l21,l22,loci,nsz,irec2)
   !--------------------------------------------------------------------
   ! Supervises the changing of matrix format to standard 4C.
   !--------------------------------------------------------------------
   use util ! provides error,repoz
   use endf ! provides endf routines and variables
   ! externals
   integer::l21,l22,loci,nsz,irec2
   integer(4)::ia(*)
   real(4)::b(*)
   ! internals
   integer::nosiza,lrsize,nrmx,nrec,npass,iasiz,j,i,nj,nb,nwds
   integer::nl,m,mn,l

   !--assign storage

   ! nsiza--space available in b(i)
   ! irsize--length of one input record
   ! nosiza--space available in b(i) for output matrix
   ! lrsize--maximum length of output record
   ! nrmx--total number of records in matrix
   ! nrec--number of output records which fit in nosiza
   ! npass--number of tape passes necessary to fill entire matrix
   ! for ifopt=2 and nsblk=ngps, nrec and npass are not used

   !--set up sizes
   nsiza=isiza4-l22+1
   irsize=6+(ngps+1)*(maxord+1)*nsz
   if (irsize.lt.6*ngps) irsize=6*ngps
   if (irsize.gt.nsiza)&
     call error('isomtx','input record too large.',' ')
   nosiza=nsiza-irsize
   if (ifopt.ne.2) then
      lrsize=((maxord+1)*ngps*(ngps+1))/2
      if (nosiza.lt.lrsize) nsblk=ngps
      if (nsblk.eq.1) then
         nrmx=1
      else
         lrsize=ngps*(maxord+1)
         nrmx=ngps
      endif
   else
      lrsize=(ngps*(ngps+1))/2
      if (nosiza.lt.lrsize) nsblk=ngps
      if (nsblk.eq.1) then
         nrmx=(maxord+1)
      else
         nrmx=ngps*(maxord+1)
         lrsize=ngps
      endif
   endif
   nrec=nosiza/lrsize
   if (nrec.ge.nrmx) nrec=nrmx
   if (nrec.lt.1)&
     call error('isomtx','output record too large.',' ')
   npass=(nrmx+nrec-1)/nrec

   !--clear ia(i)
   iasiz=l22-l21
   do j=1,iasiz
      ia(j)=0
   enddo

   !--loop over data for inelastic, elastic, (n,2n) and total
   do 500 i=1,4
   lord1=0
   ng2z=0
   jlz=0
   ! loop through input tape
   do 400 nj=1,npass
   ! clear b(i)
   do j=1,nsiza
      b(j)=0
   enddo
   call repoz(nscrt1)
  220 continue
   call contio(nscrt1,0,0,e(1),nb,nwds)
   nl=l1h
   if (math.eq.0) go to 300
   if (mfh.eq.0.or.mth.eq.0) go to 220
   if (mfh.eq.6) go to 225
   call tosend(nscrt1,0,0,e(1))
   go to 220
  225 continue
   if (mfh.ne.6) go to 220
   ! check for end of block
   if (ng2z.ge.ngps) go to 405
   if (jlz.gt.maxord) go to 405
   ! check mt number
   go to (230,240,250,260),i
   ! inelastic
  230 if (mth.eq.2.or.mth.eq.16) go to 220
   if (mth.ge.6.and.mth.le.9) go to 220
   if (mth.ge.46.and.mth.le.49) go to 220
   if (mth.eq.18.or.mth.eq.19.or.mth.eq.20.or.mth.eq.21.or.mth.eq.38)&
    go to 220
   if (lord1.lt.1) lord1=1
   go to 270
   ! elastic
  240 continue
   if (mth.gt.2) go to 300
   lord1=1
   go to 270
   ! (n,2n)
  250 continue
   if (mth.gt.49) go to 300
   if (mth.lt.6) go to 220
   if (mth.gt.9.and.mth.lt.46.and.mth.ne.16) go to 220
   if (lord1.lt.1) lord1=1
   go to 270
   ! total
  260 continue
   if (mth.eq.18.or.mth.eq.19.or.mth.eq.20.or.mth.eq.21.or.mth.eq.38)&
    go to 220
   if (lord1.lt.1) lord1=1
   go to 270

   !--rearrange matrix elements
  270 continue
   ieof=0
   call shuffl(b,nl)
   if (ieof.eq.1) go to 300
   go to 220

   !--calculate bandwidths and write matrix elements to nscrt2
  300 continue
   call wrtmtx(ia(1),b(1),i,nrec,irec2)
  400 continue

   !--calculate lord(n)
  405 continue
   if (ifopt.ne.2) then
      ia(4+i)=lord1
   else if (lord1.ne.0) then
      do m=1,lord1
         mn=nscmax+(maxord+1)*(i-1)+m
         ia(mn)=1
      enddo
   endif
  500 continue

   !--calculate loci
   loci=0
   if (nsblk.eq.1) then
      do l=1,nscmax
         mn=nscmax+l
         if (ia(mn).ne.0) then
            loci=loci+1
         endif
      enddo
   else
      do l=1,nscmax
         do m=1,ngps
            mn=2*nscmax+(l-1)*ngps+m
            if (ia(mn).ne.0) then
               loci=loci+1
            endif
         enddo
      enddo
   endif
   return
   end subroutine isomtx

   subroutine shuffl(a,nl)
   !--------------------------------------------------------------------
   ! Reads one entire MT file and transfers data.
   !--------------------------------------------------------------------
   use util ! provides error
   use endf ! provides endf routines and variables
   ! externals
   integer::nl
   real(4)::a(*)
   ! internals
   integer::nz,ig,nb,nw,ng2,ig2lo,loc,jg1,il,ik,niloc,jg2,noloc,inl

   nl=l1h
   nz=l2h
   ig=0
   do while (ig.lt.ngps)
      call listio(nscrt1,0,0,e,nb,nw)
      ng2=l1h
      ig2lo=l2h
      ig=n2h
      loc=1+nw
      do while (nb.ne.0)
         if (loc+302.gt.maxe) call error('shuffl',&
           'size of endf input array exceeded',' ')
         call moreio(nscrt1,0,0,e(loc),nb,nw)
         loc=loc+nw
      enddo
      !--calculate lord1
      if (nl.gt.lord1) then
         lord1=nl
         if (lord1.gt.maxord+1) lord1=maxord+1
      endif
      jg1=ngps+1-ig
      !--step through the input record
      do il=1,nl
         if (il.le.maxord+1) then
            do ik=2,ng2
               niloc=6+il+nz*nl*(ik-1)
               jg2=ngps+1-(ig2lo+ik-2)
               !--calculate noloc, the corresponding
               !--word in the output matrix
               !--ifopt=1
               if (ifopt.eq.1) then
                  noloc=irsize+(maxord+1)*(jg2*(jg2-1)-ng2z&
                    *(ng2z+1))/2+(il-1)*jg2+jg2-jg1+1
               !--ifopt=2     nsblk=1
               else if (nsblk.eq.1) then
                  noloc=irsize+(jg2-jg1)+(jg2*(jg2-1))/2+&
                    ((il-jlz-1)*ngps*(ngps+1))/2+1
               !--ifopt=2     nsblk=ngps
               else
                  inl=jlz+1
                  if (inl.ge.il) then
                     noloc=irsize+1+(jg2-jg1)+(jg2*(jg2-1)&
                       -ng2z*(ng2z+1))/2
                  else
                     noloc=irsize+1+(jg2-jg1)+(jg2*(jg2-1))/2+&
                       ((il-inl)*ngps*(ngps+1)-ng2z*(ng2z+1))/2
                  endif
               endif
               if (noloc.gt.irsize.and.noloc.le.nsiza) then
                  a(noloc)=a(noloc)+e(niloc)
               endif
            enddo
         endif
      enddo
   enddo
   return
   end subroutine shuffl

   subroutine wrtmtx(ia,a,i,nrec,irec2)
   !--------------------------------------------------------------------
   ! Calculates bandwidths and writes matrix elements to nscrt2.
   !--------------------------------------------------------------------
   ! externals
   integer::i,nrec,irec2
   integer(4)::ia(*)
   real(4)::a(*)
   ! internals
   integer::nzero,nupp,nrc,nlow,nrnge,inc,jg,jg2,ibw,nrd,noloca
   integer::jbloc,lord,ndex,mdex,jj,loc,jm,jk,lout,lin,nwds,inl
   integer::narc,ltst,ii
   real(kr),parameter::eps=1.e-10_kr

   nzero=irsize
   if (ifopt.eq.2) go to 225

   !--ifopt=2, nsblk=ngps (calculate only pzero bandwidths)
   nupp=nzero
   do 210 nrc=1,nrec
   if (nsblk.eq.1) go to 561
   nlow=nupp+1
   if ((ng2z+nrc).gt.nsblk) go to 210
   nrnge=(maxord+1)*(ng2z+nrc)
   go to 562
  561 continue

   !--ifopt=1, nsblk=1
   nlow=nzero+1
   nrnge=(ngps*(ngps+1)*(maxord+1))/2
  562 continue
   nupp=nlow+nrnge-1
   inc=(ngps-1)/nsblk+1
   do jg=1,inc
      jg2=ng2z+(nrc-1)+jg
      ibw=0
      do nrd=1,jg2
         noloca=nlow+(nrd-1)+((maxord+1)*jg*(jg-1))/2
         if (abs(a(noloca)).ge.eps) then
            ibw=nrd
         endif
      enddo
      jbloc=2*nscmax+jg2+(i-1)*ngps
      ia(jbloc)=ibw
   enddo

   !--thin and copy matrix records
   lord=maxord+1
   ndex=nlow-1
   mdex=ndex
   do jg=1,inc
      jj=ng2z+(nrc-1)+jg
      loc=2*nscmax+(i-1)*ngps+jj
      ibw=ia(loc)
      if (ibw.gt.0) then
         do jm=1,lord
            do jk=1,ibw
               lout=ndex+(jm-1)*ibw+jk
               lin=mdex+(jm-1)*jj+jk
               a(lout)=a(lin)
            enddo
         enddo
         ndex=lout
      endif
      mdex=mdex+jj*lord
   enddo
   nwds=ndex-nlow+1
   if (nwds.gt.0) then
      irec2=irec2+1
      write(nscrt2)(a(nlow+ii-1),ii=1,nwds)
   endif
  210 continue
   ng2z=ng2z+nrec
   jlz=0
   go to 300
  225 continue
   if (nsblk.gt.1) go to 250

   !--ifopt=2, nsblk=1
   do nrc=1,nrec
      nlow=nzero+1+((nrc-1)*ngps*(ngps+1))/2
      nrnge=(ngps*(ngps+1))/2
      nupp=nlow+nrnge-1
      do jg2=1,ngps
         ibw=0
         do nrd=1,jg2
            noloca=nlow-1+nrd+(jg2*(jg2-1))/2
            if (abs(a(noloca)).ge.eps) then
               ibw=nrd
            endif
         enddo
         jbloc=2*nscmax+(i-1)*ngps*(maxord+1)+(nrc+jlz-1)*ngps+jg2
         ia(jbloc)=ibw
      enddo

      !--thin and copy matrix records
      lord=1
      ndex=nlow-1
      mdex=ndex
      do jj=1,ngps
         loc=2*nscmax+(i-1)*ngps*(maxord+1)+(nrc+jlz-1)*ngps+jj
         ibw=ia(loc)
         if (ibw.gt.0) then
            jm=1
            do jk=1,ibw
               lout=ndex+(jm-1)*ibw+jk
               lin=mdex+(jm-1)*jj+jk
               a(lout)=a(lin)
            enddo
            ndex=lout
         endif
         mdex=mdex+jj*lord
      enddo
      nwds=ndex-nlow+1
      if (nwds.gt.0) then
         irec2=irec2+1
         write(nscrt2)(a(nlow+ii-1),ii=1,nwds)
      endif
   enddo
   ng2z=0
   jlz=jlz+nrec
   go to 300
  250 continue

   !--ifopt=2, nsblk=ngps
   inl=1
   nupp=nzero
   nrnge=ng2z+1
   narc=nrec*(maxord+1)
   do 290 nrc=1,narc
   jg2=nrc+ng2z-(inl-1)*ngps
   if (jg2.le.ngps) go to 260
   jg2=jg2-ngps
   inl=inl+1
   ! do not write beyond the end of the scattering matrix
   ltst=jlz+inl-maxord-1
   if (ltst.gt.0) jlz=maxord+1
   if (jlz.gt.maxord) go to 300
  260 continue
   nlow=nupp+1
   nrnge=jg2
   nupp=nlow+nrnge-1
   if (nupp.le.nsiza) go to 270
   jlz=jlz+inl-1
   ng2z=jg2-1
   if (ng2z.le.0) ng2z=ngps
   go to 300
  270 continue
   ibw=0
   do nrd=1,jg2
      noloca=nlow+nrd-1
      if (abs(a(noloca)).ge.eps) then
         ibw=nrd
      endif
   enddo
   jbloc=2*nscmax+(i-1)*ngps*(maxord+1)+(jlz+inl-1)*ngps+jg2
   ia(jbloc)=ibw

   !--thin and copy matrix records
   lord=1
   ndex=nlow-1
   mdex=ndex
   jj=jg2
   loc=2*nscmax+(i-1)*ngps*(maxord+1)+(jlz+inl-1)*ngps+jj
   ibw=ia(loc)
   if (ibw.gt.0) then
      jm=1
      do jk=1,ibw
         lout=ndex+(jm-1)*ibw+jk
         lin=mdex+(jm-1)*jj+jk
         a(lout)=a(lin)
      enddo
      ndex=lout
   endif
   mdex=mdex+jj*lord
   nwds=ndex-nlow+1
   if (nwds.gt.0) then
      irec2=irec2+1
      write(nscrt2)(a(nlow+ii-1),ii=1,nwds)
   endif
  290 continue
  300 continue
   return
   end subroutine wrtmtx

   subroutine pisotx(n4)
   !--------------------------------------------------------------------
   ! Controls printing of the ISOTXS file.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use util ! provides error,a10,repoz
   ! externals
   integer::n4
   ! internals
   integer::i,irec1,nwds,iver,ngroup,niso,ichist,nscmax,nsblok
   integer::j,next,npt,jread,l,ichi,ifis,ialf,inp,in2n,ind,int
   integer::ltot,ltrn,istrpd,l1,l2,j1,iblk,jblk,mblk,k,nblk,n
   integer::ntype,idsct,lonow,lordn,kmax,inc,ii,lordn1,jl,ju,jj
   integer::ndn,jk,np,lll,ip,ip1,ip2,jk1ip,nwh,llr,nwr,lli,nwi,nn
   real(kr)::ad
   character(8)::hn,hu,hs,hblock
   integer::ird(10)
   character(6)::htype(20)
   character(10)::field(10)
   character(6)::htotal=' total'
   character(6)::helas='elastc'
   character(6)::hinel='inelas'
   character(6)::hn2n='   n2n'
   character(6)::hstrpl=' strpl'
   character(6)::hsttpl='stotpl'
   character(6)::hsngam=' sngam'
   character(6)::hsfis='  sfis'
   character(6)::hsnutt='snutot'
   character(6)::hchiso=' chiso'
   character(6)::hsnalf=' snalf'
   character(6)::hsnp='   snp'
   character(6)::hsn2n='  sn2n'
   character(6)::hsnd='   snd'
   character(6)::hsnt='   snt'
   character(6)::hingp='in gp '
   character(6)::houtgp='out gp'
   character(6)::hblank='      '
   character(6)::horder=' order'
   character(6)::hundl='------'
   integer::m1=1
   integer::izero=0
   do i=1,10
      ird(i)=1
   enddo
   irec1=0
   call repoz(-n4)

   !--file identification
   nwds=1+3*mult
   irec1=irec1+1
   read(n4)(ha(i),i=1,3),ia(nwds)
   hn=ta(1)
   hu=ta(2)
   hs=ta(3)
   iver=ia(nwds)
   write(nsyso,'(/'' ***file '',a6,'' -- version'',i3,&
     &'' -- unit'',i3,''***'')') hn,iver,n4
   write(nsyso,'('' **user identification**'',a6,3x,a6)') hu,hs

   !--file control
   irec1=irec1+1
   read(n4)(ia(i),i=1,8)
   ngroup=ia(1)
   niso=ia(2)
   ichist=ia(6)
   nscmax=ia(7)
   nsblok=ia(8)
   if (ird(1).eq.1) write(nsyso,'(/&
     &'' file control parameters''//&
     &''  ngroup    number of energy groups in set         '',i12/&
     &''  niso      number of isotopes in set              '',i12/&
     &''  maxup     maximum number of upscatter groups     '',i12/&
     &''  maxdn     maximum number of downscatter groups   '',i12/&
     &''  maxord    maximum scattering order               '',i12/&
     &''  ichist    set fission spectrum flag              '',i12/&
     &14x,''ichist=1'',7x,''set vector''/&
     &20x,''=ngroup, set matrix''/&
     &''  nscmax    maximum number of blocks of            '',i12/&
     &14x,''scattering data''/&
     &''  nsblok    blocking control for scattering data   '',i12)')&
     (ia(i),i=1,8)

   !--file data
   nwh=12+niso
   llr=nwh*mult
   nwr=1+ngroup*(2+ichist*(2/(ichist+1)))
   lli=llr+nwr
   nwi=niso
   nwds=mult*nwh+nwr+nwi
   irec1=irec1+1
   if (nwds.ge.isiza4)&
     call error('pisotx','input record too large.',' ')
   read(n4)(ha(i),i=1,nwh),(a(llr+i),i=1,nwr),(ia(lli+i),i=1,nwi)
   if (ird(2).eq.1) then
      write(nsyso,'(/1x,12a6)') (ta(i),i=1,12)
      write(nsyso,'(/4x,''isotope'',4x,''name'')')
      do i=1,niso
         j=12+i
         write(nsyso,'(6x,i3,6x,a6)') i,ta(j)
      enddo
      next=(12+niso)*mult+1
      if (ichist.eq.1) then
         write(nsyso,'(/'' fission spectrum'')')
         call wot(a(next),ichist,ngroup,m1,hingp,houtgp,hblank,nsyso)
         next=next+ngroup*ichist
      endif
      write(nsyso,'(/11x,''neutron velocity'',9x,''upper energy''/&
        &'' group'',13x,''(cm/sec)'',17x,''(ev)''//)')
      write(nsyso,'(1x,i5,1p,2e21.6)')&
        (j,a(next-1+j),a(next-1+ngroup+j),j=1,ngroup)
      next=next+2*ngroup
      write(nsyso,'(27x,1p,e21.6)') a(next)
      next=next+1
      write(nsyso,'(/'' number of records to be skipped''//&
        &5x,''isotope'',2x,''number'')')
      do i=1,niso
         write(nsyso,'(i10,i8)') i,ia(next+i-1)
      enddo
   endif

   !--set chi data
   if (ichist.gt.1) then
      nwr=ngroup*ichist
      nwi=ngroup
      nwds=nwr+nwi
      npt=1
      llr=npt-1
      lli=llr+nwr
      irec1=irec1+1
      if (nwds.ge.isiza4)&
        call error('pisotx','input record too large.',' ')
      read(n4)(a(llr+i),i=1,nwr),(ia(lli+i),i=1,nwi)
      if (ird(3).eq.1) then
         write(nsyso,'(/'' set fission spectrum'')')
         call wot(a(1),ichist,ngroup,m1,hingp,houtgp,hblank,nsyso)
      endif
   endif

   !--isotope control and group independent data
   do 220 i=1,niso
   nwh=3
   llr=mult*nwh
   nwr=6
   lli=llr+nwr
   nwi=11+2*nscmax+2*ngroup*nscmax
   nwds=mult*nwh+nwr+nwi
   irec1=irec1+1
   if (nwds.ge.isiza4)&
     call error('pisotx','input record too large.',' ')
   read(n4)(ha(j),j=1,nwh),(a(llr+j),j=1,nwr),(ia(lli+j),j=1,nwi)
   jread=1+nwds
   l=3*mult+7
   ichi=ia(l+1)
   ifis=ia(l+2)
   ialf=ia(l+3)
   inp=ia(l+4)
   in2n=ia(l+5)
   ind=ia(l+6)
   int=ia(l+7)
   ltot=ia(l+8)
   ltrn=ia(l+9)
   istrpd=ia(l+10)
   if (ird(4).eq.1) then
      write(nsyso,'(/'' isotope'',i3)') i
      l1=3*mult
      l2=l1+6
      write(nsyso,'(/&
        &'' isotope control parameters''//&
        &''  habsid  absolute isotope label'',26x,a6/&
        &''  hident  library identifier'',30x,a6/&
        &''  hmat    isotope identification'',26x,a6/&
        &''  amass   gram atomic weight'',24x,1p,e12.5/&
        &''  efiss   thermal energy/fission (w*sec/fiss)'',7x,e12.5/&
        &''  ecapt   thermal energy/capture (w*sec/capt)'',7x,e12.5/&
        &''  temp    isotope temperature (deg k)'',15x,e12.5/&
        &''  sigpot  ave. potential scattering (barns/atom)'',&
        &4x,e12.5/&
        &''  adens   reference atom density (a/b*cm)'',11x,e12.5/&
        &''  kbr     isotope classification'',20x,i12/&
        &''  ichi    fission spectrum flag'',&
        &'' (0/1/n=set chi/vector/matrix)'',i3)')&
        (ta(j),j=1,3),(a(l1+j),j=1,6),(ia(l2+j),j=1,2)
      l1=l2+2
      write(nsyso,'(&
        &''  ifis    (n,f) x-sec flag (0/1=no/yes)'',13x,i12/&
        &''  ialf    (n,a) x-sec flag (0/1=no/yes)'',13x,i12/&
        &''  inp     (n,p) x-sec flag (0/1=no/yes)'',13x,i12/&
        &''  in2n    (n,2n) x-sec flag (0/1=no/yes)'',12x,i12/&
        &''  ind     (n,d) x-sec flag (0/1=no/yes)'',13x,i12/&
        &''  int     (n,t) x-sec flag (0/1=no/yes)'',13x,i12/&
        &''  ltot    number of total x-sec moments'',13x,i12/&
        &''  ltrn    number of transport x-sec moments'',9x,i12/&
        &''  istrpd  number of transport x-sec directions'',i18/)')&
        (ia(l1+j),j=1,9)
      l1=l1+9
      write(nsyso,&
        '(1x,''block'',5x,''type'',4x,''ident'',3x,''orders''/)')
      do j=1,nscmax
         read(htotal,'(a6)') htype(j)
         if (ia(l1+j).ge.100) read(helas,'(a6)') htype(j)
         if (ia(l1+j).ge.200) read(hinel,'(a6)') htype(j)
         if (ia(l1+j).ge.300) read(hn2n,'(a6)') htype(j)
      enddo
      write(nsyso,'(1x,i5,3x,a6,2i9)')&
        (j,htype(j),ia(l1+j),ia(l1+nscmax+j),j=1,nscmax)
      j1=l1+2*nscmax
      write(nsyso,'(/&
        &'' scattering bandwidth and in-group'',&
        &'' scattering position'')')
      iblk=(nscmax+1)/10+1
      jblk=nscmax/iblk+1
      mblk=1
      do k=1,iblk
         nblk=mblk+jblk-1
         if (nblk.gt.nscmax) nblk=nscmax
         write(nsyso,'(/'' group/block'',3x,20i5)')&
           (j,j=mblk,nblk),(j,j=mblk,nblk)
         do n=1,ngroup
            write(nsyso,'(1x,i5,9x,20i5)')&
              n,(ia(j1+j*ngroup-ngroup+n),j=mblk,nblk),&
              (ia(j1+ngroup*nscmax+j*ngroup-ngroup+n),j=mblk,nblk)
         enddo
         mblk=nblk+1
      enddo
   endif

   !--principal cross sections
   nwds=(1+ltrn+ltot+ialf+inp+in2n+ind+int+istrpd+2*ifis+&
     ichi*(2/(ichi+1)))*ngroup
   irec1=irec1+1
   read(hstrpl,'(a6)') htype(1)
   read(hsttpl,'(a6)') htype(2)
   read(hsngam,'(a6)') htype(3)
   j=4
   if (ifis.gt.0) read(hsfis,'(a6)') htype(j)
   if (ifis.gt.0) read(hsnutt,'(a6)') htype(j+1)
   if (ifis.gt.0) j=j+2
   if (ichi.eq.1) read(hchiso,'(a6)') htype(j)
   if (ichi.eq.1) j=j+1
   if (ialf.eq.1) read(hsnalf,'(a6)') htype(j)
   if (ialf.eq.1) j=j+1
   if (inp.eq.1) read(hsnp,'(a6)') htype(j)
   if (inp.eq.1) j=j+1
   if (in2n.eq.1) read(hsn2n,'(a6)') htype(j)
   if (in2n.eq.1) j=j+1
   if (ind.eq.1) read(hsnd,'(a6)') htype(j)
   if (ind.eq.1) j=j+1
   if (int.eq.1) read(hsnt,'(a6)') htype(j)
   if (nwds.ge.isiza4)&
     call error('pisotx','input record too large.',' ')
   read(n4)(a(jread+k-1),k=1,nwds)
   if (ird(5).eq.1) then
      write(nsyso,'(/'' principal cross-sections'')')
      ntype=nwds/ngroup
      iblk=(ntype-1)/8+1
      jblk=ntype/iblk+1
      mblk=1
      do k=1,iblk
         nblk=mblk+jblk-1
         if (nblk.gt.ntype) nblk=ntype
         write(nsyso,'(/'' group'',9(3x,a8,''  ''))')&
           (htype(j),j=mblk,nblk)
         do n=1,ngroup
            write(nsyso,'(1x,i3,2x,1p,9e13.6)')&
              n,(a(jread-1+n+j*ngroup-ngroup),j=mblk,nblk)
         enddo
         mblk=nblk+1
      enddo
   endif

   !--isotope chi data
   if (ichi.gt.1) then
      nwr=ngroup*ichi
      nwi=ngroup
      nwds=nwr+nwi
      npt=jread
      llr=npt-1
      lli=llr+nwr
      irec1=irec1+1
      if (nwds.ge.isiza4)&
         call error('pisotx','input record too large.',' ')
      read(n4)(a(llr+i),nn=1,nwr),(ia(lli+nn),nn=1,nwi)
      if (ird(6).eq.1) then
         write(nsyso,'('' isotope fission spectrum'')')
         call wot(a(jread),ichi,ngroup,m1,hingp,houtgp,hblank,nsyso)
      endif
   endif

   !--scattering cross sections
   write(nsyso,'(/'' scattering matrices'')')
   do 221 n=1,nscmax
   k=n+3*mult+17
   idsct=ia(k)
   idsct=idsct/100
   lonow=ia(k)-100*idsct
   k=k+nscmax
   lordn=ia(k)
   idsct=idsct+1
   if (idsct.eq.1) then
      read(htotal,'(a6)') hblock
   else if (idsct.eq.2) then
      read(helas,'(a6)') hblock
   else if (idsct.eq.3) then
      read(hinel,'(a6)') hblock
   else
      read(hn2n,'(a6)') hblock
   endif
   if (lordn.ne.0) go to 160
   write(nsyso,'(/'' no data for '',a6,''  p'',i1)') hblock,lonow
   go to 221
  160 continue
   kmax=0
   k=3*mult+17+2*nscmax+(n-1)*ngroup
   inc=(ngroup-1)/nsblok+1
   npt=jread
   kmax=0
   if (lonow.eq.0) write(nsyso,&
     '(/'' matrix ** '',a6,'' **'')') hblock
   if (lonow.gt.0) write(nsyso,&
     '(/'' matrix ** '',a6,'' ** p'',i1)') hblock,lonow
   if (lordn.eq.1) write(nsyso,'(/&
     &''  final  initl    xsec vs initl group''/&
     &''  group  group'',5x,&
     &''-0'',10x,''-1'',10x,''-2'',10x,''-3'')')
   if (lordn.eq.1) write(nsyso,'(''  -----  -----  '',9(2a6))')&
     (hundl,hundl,ii=1,5)
   lordn1=lordn-1
   if (lordn.gt.1) write(nsyso,'(/&
     &''  final  initl    xsec vs legendre order''/&
     &''  group  group'',4x,9(a6,1x,i1,4x))')&
     (horder,ii,ii=izero,lordn1)
   if (lordn.gt.1) write(nsyso,'(''  -----  -----  '',9(2a6))')&
     (hundl,hundl,ii=1,lordn)
   do 180 j=1,nsblok
   jl=(j-1)*inc+1
   ju=j*inc
   nwds=0
   do jj=jl,ju
      l=k+jj
      nwds=nwds+ia(l)
      kmax=ia(l)+kmax
   enddo
   nwds=nwds*lordn
   if (nwds.eq.0) go to 180
   irec1=irec1+1
   if (nwds.ge.isiza4)&
     call error('pisotx','input record too large.',' ')
   read(n4)(a(npt+ii-1),ii=1,nwds)
   if (ird(7).ne.1) go to 180
   ndn=npt
   do 175 jj=jl,ju
   l=k+jj
   nwds=lordn*ia(l)
   if (nwds.eq.0) go to 175
   jk=j+jj/jl-1
   np=ia(l)
   if (lordn.gt.1) go to 178
   lll=ndn-1
   ip=1
   ip1=ip
  171 continue
   if (a(lll+ip1).ne.0.) go to 172
   ip1=ip1+1
   if (ip1.gt.np) go to 176
   go to 171
  172 continue
   ip2=ip1+5
   if (ip2.gt.np) ip2=np
  173 continue
   if (a(lll+ip2).ne.0.) go to 174
   ip2=ip2-1
   if (ip2.eq.ip1) go to 174
   go to 173
  174 continue
   jk1ip=jk+1-ip1
   do ii=ip1,ip2
      ad=a(lll+ii)
      call a10(ad,field(ii-ip1+1))
   enddo
   if (ip.eq.1) write(nsyso,'(3x,i3,4x,i3,2x,1p,9a12)')&
     jk,jk1ip,(field(ii),ii=1,ip2-ip1+1)
   if (ip.gt.1) write(nsyso,'(10x,i3,2x,1p,9a12)')&
     jk1ip,(field(ii),ii=1,ip2-ip1+1)
   ip=ip+1
   ip1=ip2+1
   if (ip1.le.np) go to 171
  176 continue
   ndn=ndn+nwds
   go to 175
  178 continue
   ip1=1
   do 179 ip=1,np
   lll=ndn-1+ip-np
   do 181 ii=1,lordn
   if (a(np*ii+lll).ne.0.) go to 182
  181 continue
   go to 179
  182 continue
   jk1ip=jk+1-ip
   do ii=1,lordn
      ad=a(np*ii+lll)
      call a10(ad,field(ii))
   enddo
   if (ip1.eq.1) write(nsyso,'(3x,i3,4x,i3,2x,1p,9a12)')&
     jk,jk1ip,(field(ii),ii=1,lordn)
   if (ip1.gt.1) write(nsyso,'(10x,i3,2x,1p,9a12)')&
     jk1ip,(field(ii),ii=1,lordn)
   ip1=ip1+1
  179 continue
   ndn=ndn+nwds
  175 continue
  180 continue
  221 continue
  220 continue
   return
   end subroutine pisotx

   subroutine cbrkxs
   !--------------------------------------------------------------------
   ! Write a CCCC-IV BRKOXS interface file.
   !--------------------------------------------------------------------
   use util ! provides repoz
   ! internals
   integer::nwds,irec,l1h,ngps,nisosh,nsigpt,ntempt,nreact
   integer::iblk,l3,l4,irec3,i,loc,jbl,jbh,ntap,ntat,nblok
   integer::ml,mu,j,l3h,nwh,l3r,nwr,l3i,nwi,k
   character(6)::hname='brkoxs'

   !--set first pointer
   nwds=3*mult+1+6
   next=1
   l1=next
   next=l1+nwds
   next=next+1
   irec=0

   !--read user input
   call ruinbr

   !--create and write brkoxs data sets on scratch tapes
   call brkdat

   !--file identification
   nwds=3*mult+1
   l1h=(l1-1)/mult+1
   read(hname,'(a6)') ta(l1h)
   ha(l1h+1)=huse(1)
   ha(l1h+2)=huse(2)
   ia(3*mult+l1)=ivers
   irec=irec+1
   call repoz(-nbrks)
   write(nbrks)(ha(i),i=1,3),ia(nwds)

   !--file control
   l2=l1+3*mult+1
   ngps=ia(l2)
   nisosh=ia(l2+1)
   nsigpt=ia(l2+2)
   ntempt=ia(l2+3)
   nreact=ia(l2+4)
   iblk=ia(l2+5)
   irec=irec+1
   write(nbrks)(ia(l2-1+i),i=1,6)

   !--file data
   l3=l2+6
   l3=l3+1
   l3h=(l3-1)/mult
   nwh=mult*nisosh
   l3r=l3+mult*nwh-1
   nwr=nsigpt+ntempt+ngps+1
   l3i=l3r+nwr
   nwi=4*nisosh
   nwds=mult*nwh+nwr+nwi
   irec=irec+1
   write(nbrks)(ha(l3h+i),i=1,nwh),(a(l3r+i),i=1,nwr),&
     (ia(l3i+i),i=1,nwi)
   l4=l3+nwds

   !--loop over isotopes
   irec3=0
   call repoz(-nscrt3)
   do i=1,nisosh

      !--self-shielding factors
      loc=l3+mult*nisosh+nsigpt+ntempt+ngps+1+i-1
      jbl=ia(loc)
      jbh=ia(loc+nisosh)
      ntap=ia(loc+2*nisosh)
      ntat=ia(loc+3*nisosh)
      if (iblk.eq.0) then
         nblok=1
         ml=1
         mu=nreact
      else
         nblok=nreact
         ml=1
         mu=1
      endif
      nwds=ntap*ntat*(jbh-jbl+1)*(mu-ml+1)
      do j=1,nblok
         irec3=irec3+1
         read(nscrt3)(a(l4+k-1),k=1,nwds)
         irec=irec+1
         write(nbrks)(a(l4+k-1),k=1,nwds)
      enddo

      !--cross sections
      nwds=6*ngps
      irec3=irec3+1
      read(nscrt3)(a(l4+j-1),j=1,nwds)
      irec=irec+1
      write(nbrks)(a(l4+j-1),j=1,nwds)
   enddo
   next=l1
   return
   end subroutine cbrkxs

   subroutine ruinbr
   !--------------------------------------------------------------------
   ! Reads user input pertinent to the BRKOXS file.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   ! internals
   integer::i

   !--read sig0 and temperature data
   read(nsysi,*) nti,nzi
   nreact=6
   if (nti.gt.0) read(nsysi,*) (atem(i),i=1,nti)
   if (nzi.gt.0) read(nsysi,*) (asig(i),i=1,nzi)
   end subroutine ruinbr

   subroutine brkdat
   !--------------------------------------------------------------------
   ! Processes all BRKOXS data and writes on scratch tape.
   ! Data stored in a(i) as follows
   !
   ! location               variable                 length
   ! --------               --------                 ------
   !    l1                  *brkoxs*                  mult
   !                         huse(1)                  mult
   !                         huse(2)                  mult
   !                         ivers                       1
   !                         ngroup                      1
   !                         niso                        1
   !                         nsigpt                      1
   !                         ntempt                      1
   !                         nreact                      1
   !                         iblk                        1
   !                         dummy                       1
   !    l9                   hisonm              niso*mult
   !    l10                  emax                   ngroup
   !    l11                  emin                        1
   !    l12                   jbfl                    niso
   !    l13                   jbfh                    niso
   !    l14                  ntabp                    niso
   !    l15                  ntabt                    niso
   !    l16                    x                  niso*nzj
   !    l17                    tb                 niso*ntj
   !    l18               input record                  -
   !    l19      infinite dilution cross section 13*ngroup
   !    l20          self shielding factors             -
   !
   !--------------------------------------------------------------------
   use util ! provides repoz,error,mess
   use endf ! provides endf routines and variables
   ! internals
   integer::ntl,nzl,i,nwds,l9,l10,l11,nwords,l12,l14,l15,l16,l17,l18
   integer::nb,nw,nz,lz,nt,ntw,ngpn,loc,iqd,nenpts,l19,l20,nosize
   integer::iffile,inv,jiso,irec2,j,nrec,nlines,iin,nwl,k,lcsg0
   integer::l,itzro,jj,nwd,locsg0,locsg,mtc,nwon,nrnge,iv,noloc
   integer::ndex,nzf,jz,mdex,ntf,jt,irec3
   real(kr)::c,abe,test,amd
   integer::imusd(200),ntabt(200),ntabp(200)
   character(60)::strng
   real(kr),parameter::eps=1.e-04_kr
   real(kr),parameter::eps2=1.e-10_kr

   !--file data
   ia(l1+3*mult+1)=ngps
   ia(l1+3*mult+2)=niso
   ia(l1+3*mult+5)=nreact
   ntl=iabs(nti)
   nzl=iabs(nzi)
   do i=1,20
      item(i)=0
      isig(i)=0
   enddo

   !--assign pointers
   nwds=mult*niso
   l9=next
   l10=l9+nwds
   l11=l10+ngps
   nwords=1
   l12=l11+nwords
   l13=l12+niso
   l14=l13+niso
   l15=l14+niso
   l16=l15+niso
   nwds=niso*nzl
   l17=l16+nwds
   nwds=niso*ntl
   next=l17+nwds
   l18=next

   !--process header for first material
   call repoz(nin)
   call tpidio(nin,0,0,e(1),nb,nw)
  100 continue
   call contio(nin,0,0,e(1),nb,nw)
   if (mfh.eq.0.or.mth.eq.0) go to 100
   nz=l2h
   lz=6
   nt=n2h
   ntw=nt
   call listio(nin,0,0,e(1),nb,nw)
   ngpn=l1h
   loc=l1+3*mult+1
   if (ngpn.ne.ia(loc))&
     call error('brkdat','incompatible group structures.',' ')
   loc=1+nw
   do while (nb.ne.0)
      if (loc+302.gt.maxe) call error('brkdat',&
        'max size of endf record exceeded',' ')
      call moreio(nin,0,0,e(loc),nb,nw)
      loc=loc+nw
   enddo
   iqd=1+lz+ntw+nz+ngpn
   ! reorder group bounds, store in a(l10)
   nenpts=ngpn+1
   do i=1,nenpts
      a(i-1+l10)=e(-i+1+iqd)
   enddo

   !--more pointers
   irsize=200
   irsize=(maxord+1)*(ngpn+1)*7+100
   if (irsize.gt.2000) irsize=2000
   next=l18+irsize
   nwds=(13+nzl)*ngps
   l19=next
   next=l19+nwds
   l20=next

   !--determine sub blocking for self shielding file
   !--nosize--space available in a(i) for output file
   !--iffile--length of entire self shielding file
   !--nsblk--number of s-s records to be written on scratch tape
   !--lrsize--length of one s-s record
   nosize=isiza4-l20+1
   iffile=nreact*ngpn*ntl*nzl
   nsblk=1
   if (iffile.gt.nosize) nsblk=nreact
   lrsize=iffile/nsblk
   inv=nreact/nsblk
   if (nsblk.gt.1) ia(l1+3*mult+6)=1
   if (nsblk.eq.1) ia(l1+3*mult+6)=0

   !--begin isotope loop
   do i=1,200
      imusd(i)=0
   enddo
   ia(l1+3*mult+3)=0
   ia(l1+3*mult+4)=0
   jiso=0
   irec3=0
   call repoz(-nscrt3)
   do 370 i=1,niso
   call repoz(nin)
   irec2=0
   call repoz(-nscrt2)
   jbl=ngpn
   jbh=0
   do j=l18,isiza4
      a(j)=0
   enddo
   do j=1,10
      mt2tem(j)=0
   enddo
   ntfl=0
   ! find a desired mat on nin file tape
   call tpidio(nin,0,0,e(1),nb,nw)
  150 continue
   call contio(nin,0,0,e(1),nb,nw)
   if (mfh.eq.0.or.mth.eq.0) go to 150
   if (math.eq.-1) go to 180
   nrec=1
   nlines=1
   do 170 j=1,niso
   if (math.ne.imat(j)) go to 170
   if (imusd(j).eq.1) go to 170
   nt=n2h
   ntw=nt
   iin=j
   matd=imat(j)
   xspo=xspot(j)
   imusd(j)=1
   nzj=0
   nzt=l2h
   call listio(nin,0,0,e(1),nb,nw)
   nwl=nint(e(5))
   nlines=nlines+(nwl+5)/6
   nrec=nrec+1
   loc=1+nw
   do while (nb.ne.0)
      if (loc+302.gt.maxe) call error('brkdat',&
        'max size of endf record exceeded',' ')
      call moreio(nin,0,0,e(loc),nb,nw)
      nrec=nrec+1
      loc=loc+nw
   enddo
   ! check which values of sig0 are present
   do 165 k=1,nzt
   lcsg0=6+ntw+k
   if (nzi.lt.0) go to 160
   do 156 l=1,nzl
   c=abs(asig(l)-e(lcsg0))
   if (c.le.eps*asig(l)) go to 160
  156 continue
   go to 165
  160 continue
   nzj=nzj+1
   isig(nzj)=k
   csig(nzj)=e(lcsg0)
  165 continue
   go to 190
  170 continue
   call tomend(nin,0,0,e(1))
   go to 150
  180 continue
   call mess('brkdat','all available mats have been processed',' ')
   go to 450
  190 continue
   call repoz(nscrt1)
   nsh=1
   ! copy desired temperatures to nscrt1
   if (nin.gt.0) nrec=nlines+1
   call skiprz(nin,-nrec)
   ntj=0
   itzro=0
   do 250 k=1,ntl
  210 continue
   call contio(nin,0,0,e(1),nb,nw)
   if (math.ne.matd) go to 255
   call listio(nin,0,0,e(1),nb,nw)
   loc=1+nw
   do while (nb.ne.0)
      if (loc+302.gt.maxe) call error('brkdat',&
        'max size of endf record exceeded',' ')
      call moreio(nin,0,0,e(loc),nb,nw)
      loc=loc+nw
   enddo
   if (itzro.gt.0) go to 235
   tzro=c1h
   itzro=1
   call tofend(nin,0,0,e(1))
   call tomend(nin,nscrt1,0,e(1))
   if (nti.lt.0) go to 230
   do 225 l=1,ntl
   if (abs(atem(l)-tzro).le.eps*tzro) go to 230
  225 continue
   go to 210
  230 continue
   ntj=ntj+1
   ctem(ntj)=tzro
   item(ntj)=ntj
   go to 250
  235 continue
   if (nti.lt.0) go to 245
   do 240 l=1,ntl
   if (abs(atem(l)-c1h).le.eps*c1h) go to 245
  240 continue
   call tomend(nin,0,0,e(1))
   go to 210
  245 continue
   call tofend(nin,0,0,e(1))
   ntj=ntj+1
   ctem(ntj)=c1h
   item(ntj)=ntj
   call tomend(nin,nscrt1,0,e(1))
  250 continue
  255 continue
   call atend(nscrt1,0)
   if (ntj.gt.0) go to 260
   write(strng,'(''no temperatures for mat='',i4)') matd
   call mess('brkdat',strng,' ')
   go to 370
  260 continue
   jiso=jiso+1
   jj=(l9-1)/mult+jiso
   read(hisonm(iin),'(a6)') ta(jj)
   ia(l1+3*mult+2)=jiso
   next=l18+irsize
   nwd=(13+nzj)*ngpn
   l19=next
   next=l19+nwd
   l20=next
   do 340 nsb=1,nsblk
   call repoz(nscrt1)
  265 continue
   call contio(nscrt1,0,0,e(1),nb,nw)
   if (math.eq.-1) go to 270
   if (mfh.eq.0.or.mth.eq.0) go to 265

   !--process cross sections (mf=3)
   if (mfh.eq.3) call xsproc(a(l18),l18,l19,l20)

   !--process scattering matrix (mf=6)
   if (mfh.eq.6) call mxproc(a(l18),l18,l20)
   call tosend(nscrt1,0,0,e(1))
   go to 265
  270 continue

   !--calculate transport self-shielding factors
   if (nsblk.eq.1.or.nsb.eq.4) then
      do j=1,ngpn
         locsg0=l19-1+3*ngpn+j
         do k=1,ntj
            do l=1,nzj
               locsg=l20-1+(3/nsb)*ngpn*nzj*ntj+(j-1)*nzj*ntj&
                 +nzj*(k-1)+l
               if (a(locsg0).ge.eps2) then
                  a(locsg)=a(locsg)/a(locsg0)
               else
                  a(locsg)=1
               endif
            enddo
         enddo
      enddo
   endif

   !--check whether transport cross section properly calculated
   if (ntfl.ne.1) then
      mtc=0
      do k=1,10
         mtc=mtc+mt2tem(k)
      enddo
      if (mtc.ne.ntj) then
         call mess('brkdat',&
           'need elastic matrices at higher temperatures',' ')
      endif
   endif
   ! write data to scratch tape
   if (nsb.eq.1) then
      nwon=l19+7*ngpn
      nrnge=6*ngpn
      irec2=irec2+1
      write(nscrt2)(a(nwon+j-1),j=1,nrnge)
   endif
   do iv=1,inv
      do j=1,ngpn
         do k=1,ntj
            do l=1,nzj
               noloc=l20-1+(iv-1)*ngpn*nzj*ntj+(j-1)*ntj*nzj&
                 +(k-1)*nzj+l
               abe=abs(a(noloc))
               if (abe.lt.eps2) a(noloc)=1
               test=1
               test=test/100
               if (a(noloc).lt.test) a(noloc)=test
               if (abs(a(noloc)-1.0).ge.sstol) then
                  if (jbh.lt.j) jbh=j
                  if (jbl.gt.j) jbl=j
                  ntat(k)=1
                  ntap(l)=1
               endif
            enddo
         enddo
      enddo
   enddo
   ia(l12-1+i)=jbl
   ia(l13-1+i)=jbh
   ! write f-factors to nscrt2
   irec2=irec2+1
   write(nscrt2)(a(l20+j-1),j=1,lrsize)
  340 continue

   !--calculate isotope flags
   ntabt(i)=0
   do j=1,ntj
      if (ntat(j).ne.0) then
         ntabt(i)=ntabt(i)+1
         ndex=l17+(i-1)*ntl+ntabt(i)-1
         a(ndex)=ctem(j)
      endif
   enddo
   ntabp(i)=0
   do j=1,nzj
      if (ntap(j).ne.0) then
         ntabp(i)=ntabp(i)+1
         ndex=l16+(i-1)*nzl+ntabp(i)-1
         a(ndex)=csig(j)
      endif
   enddo
   ia(l14+i-1)=ntabp(i)
   ia(l15+i-1)=ntabt(i)

   !--thin data and write on scratch tape
   call thnwrt(a(l18),irec3)

   !--close isotope loop
  370 continue

   !--rearrange and thin file data

   !--build file data record
   nrnge=ngpn+1
   do i=1,nrnge
      a(l20+i-1)=a(l10+i-1)
      a(l10+i-1)=0
   enddo
   nrnge=4*niso
   do i=1,nrnge
      ndex=ngpn+i
      ia(l20+ndex)=ia(l10+ndex)
      ia(l10+ndex)=0
   enddo
   ! thin x and tb, and store in l10
   ndex=l10-1
   ia(l1+3*mult+3)=0
   do i=1,jiso
      nzf=ia(l20-1+ngpn+1+2*niso+i)
      do jz=1,nzf
         ndex=ndex+1
         mdex=l16+(i-1)*nzl+jz-1
         amd=a(mdex)
         a(ndex)=log10(amd)
         test=1
         test=test/100000
         if (abs(a(ndex)).lt.test) a(ndex)=0
         a(mdex)=0
         ia(l1+3*mult+3)=ia(l1+3*mult+3)+1
      enddo
   enddo
   ia(l1+3*mult+4)=0
   do i=1,jiso
      ntf=ia(l20-1+ngpn+1+3*niso+i)
      do jt=1,ntf
         ndex=ndex+1
         mdex=l17+(i-1)*ntl+jt-1
         a(ndex)=a(mdex)-273.16
         ia(l1+3*mult+4)=ia(l1+3*mult+4)+1
         a(mdex)=0
      enddo
   enddo
   ! move up remaining file data
   nrnge=ngpn+1
   do i=1,nrnge
      a(ndex+i)=a(l20+i-1)
      a(l20+i-1)=0
   enddo
   ndex=ndex+ngpn+1
   do i=1,4
      k=(i-1)*jiso
      l=ngpn+(i-1)*niso
      do j=1,jiso
         ia(j+k+ndex)=ia(j+l+l20)
         ia(j+l+l20)=0
      enddo
   enddo
  450 return
   end subroutine brkdat

   subroutine xsproc(a,l18,l19,l20)
   !--------------------------------------------------------------------
   ! Derives information needed from cross section file (MF=3)
   !--------------------------------------------------------------------
   use util ! provides error
   use endf ! provides endf routines and variables
   ! externals
   integer::l18,l19,l20
   real(4)::a(*)
   ! internals
   integer::inv,nl,nz,lz,nwon,nb,nw,jmt,ig,loc,noloc,i,isig2
   integer::none,j,itm,jg,imt,idum,locsg0,jz,locsg,locss
   character(60)::strng
   real(kr)::temp,abtm,sg0
   real(kr),parameter::eps=1.e-06_kr
   real(kr),parameter::zero=0

   inv=nreact/nsblk
   nl=l1h
   nz=l2h
   lz=6
   nwon=l19-l18+1

   !--read all records for this reaction
  120 continue
   call listio(nscrt1,0,0,e(1),nb,nw)
   temp=c1h
   abtm=abs(temp-tzro)-0.01e0_kr
   jmt=0
   ig=n2h
   loc=1+nw
   do while (nb.ne.0)
      if (loc+302.gt.maxe) call error('xsproc',&
         'max endf input record exceeded',' ')
      call moreio(nscrt1,0,0,e(loc),nb,nw)
      loc=loc+nw
   enddo

   !--store flux for this temperature
   if (mth.eq.1) then
      noloc=13*ngps+nzj*(ngps-ig)+irsize
      do i=1,nzj
         isig2=isig(i)*isig(i)
         a(i+noloc)=e(lz+isig2)
      enddo
   endif

   !--store infinite dilution cross sections in a(l19)
   ! total
   if (mth.ne.1) go to 140
   jmt=1
   if (nsb.ne.1) go to 200
   if (abtm.gt.eps*temp) go to 200
   a(nwon+ngps-ig)=e(lz+nl*nz+1)
   a(nwon+4*ngps-ig)=e(lz+nl*nz+nl)
   a(nwon+6*ngps-ig)=e(lz+nl*nz+nl)
   a(nwon+7*ngps-ig)=e(lz+nl)
   go to 200
   ! capture
  140 continue
   if (mth.ne.102) go to 150
   jmt=2
   if (nsb.ne.1) go to 200
   if (abtm.gt.eps*temp) go to 200
   a(nwon+2*ngps-ig)=e(lz+nl*nz+1)
   go to 200
   ! fission
  150 continue
    if (mth.ne.18.and.mth.ne.19.and.mth.ne.20.and.mth.ne.21&
      .and.mth.ne.38) go to 160
   if (mth.ne.18) go to 160
   jmt=3
   if (nsb.ne.1) go to 200
   if (abtm.gt.eps*temp) go to 200
   a(nwon+3*ngps-ig)=e(lz+nl*nz+1)
   go to 200
   ! elastic
  160 continue
   if (mth.ne.2) go to 180
   jmt=5
   if (nsb.ne.1) go to 200
   if (abtm.gt.eps*temp) go to 200
   a(nwon+5*ngps-ig)=e(lz+nl*nz+1)
   ! load geometric hard core scattering cross section into xspo
   a(nwon+8*ngps-ig)=xspo
   go to 200
   ! mubar
  180 continue
   if (mth.ne.251) go to 190
   if (nsb.ne.1) go to 195
   if (abtm.gt.eps*temp) go to 195
   a(nwon+11*ngps-ig)=e(lz+nl*nz+1)
   go to 195
   ! xi
  190 continue
   if (mth.ne.252) go to 192
   if (nsb.ne.1) go to 195
   if (abtm.gt.eps*temp) go to 195
   a(nwon+13*ngps-ig)=e(lz+nl*nz+1)
  192 continue
  195 continue
   if (ig.lt.ngps) go to 120

   !--calculate and store self-shielding factors
   !--(for transport cross section store sigma total)
  200 continue
   none=l20-l18+1
   do 330 j=1,ntj
   itm=j
   if (abs(temp-ctem(j)).le.eps*temp) go to 335
  330 continue
   go to 395
  335 jg=(ngps+1)-ig
   do 390 j=1,inv
   imt=j+nsb-1
   if (imt.eq.jmt) go to 340
   if (imt.eq.4.and.mth.eq.1) go to 345
   go to 390
  340 continue
   idum=1
      if (mth.eq.1.and.nl.eq.2) idum=2
      locsg0=(l19-l18)+(imt-1)*ngps+jg
      sg0=a(locsg0)
      go to 350
  345 idum=nl
      sg0=1
  350 do jz=1,nzj
      locsg=lz+nl*nz+nl*(isig(jz)-1)+idum
      if (e(locsg).ne.zero) then
         locss=(none-1)+(imt-nsb)*ntj*nzj*ngps+(jg-1)*nzj*ntj&
           +(itm-1)*nzj+jz
         if (sg0.eq.0.) then
            write(strng,'(''infinite f-factor '',3i6,1p,e12.3)')&
              mth,jg,jz,temp
            call mess('xsproc',strng,' ')
         endif
         if (sg0.ne.zero) then
            a(locss)=a(locss)+e(locsg)/sg0
         endif
      endif
   enddo
  390 continue
  395 continue
   if (ig.lt.ngps) go to 120
   return
   end subroutine xsproc

   subroutine mxproc(a,l18,l20)
   !--------------------------------------------------------------------
   ! Derives information needed from matrix file (MF=6)
   !--------------------------------------------------------------------
   use util ! provides error
   use endf ! provides endf routines and variables
   ! externals
   integer::l18,l20
   real(4)::a(*)
   ! internals
   integer::nl,nz,lz,nb,nw,ng2,ig2lo,ig,jg1,loc,j,niloc,jg2,noloc
   integer::inc,k,locsg,l,itemp,lrem
   real(kr)::temp,abtm
   real(kr),parameter::eps=1.e-6_kr

   nl=l1h
   nz=l2h
   lz=6

   !--process all records for this reaction
  120 continue
   call listio(nscrt1,0,0,e(1),nb,nw)
   temp=c1h
   abtm=abs(temp-tzro)-0.01e0_kr
   ng2=l1h
   ig2lo=l2h
   ig=n2h
   jg1=ngps+1-ig
   loc=1+nw
   do while (nb.ne.0)
      if (loc+302.gt.maxe)&
        call error('mxproc','max endf size exceeded',' ')
      call moreio(nscrt1,0,0,e(loc),nb,nw)
      loc=loc+nw
   enddo
   if (mth.gt.91) go to 140
   if ((mth.ge.18.and.mth.le.21).or.mth.eq.38) go to 140
   if (abtm.gt.eps*temp) go to 200
   if (nsb.gt.1) go to 200

   !--complete the cross section block
   ! elastic xsec and elastic removal
   if (mth.ne.2) go to 160
   do j=2,ng2
      niloc=lz+1+nl*nz*(j-1)
      jg2=(ngps+1)-(ig2lo+j-2)
      noloc=irsize+9*ngps+jg1
      a(noloc)=a(noloc)+e(niloc)
      if (jg2.ne.jg1) then
         noloc=irsize+11*ngps+jg1
         a(noloc)=a(noloc)+e(niloc)
      endif
   enddo
   ! inelastic scattering cross section
   ! i. e., everything but elastic and n2n including multiplicity
  160 continue
   if (mth.eq.2) go to 180
   if (mth.eq.16) go to 180
   if (mth.ge.6.and.mth.le.9) go to 180
   if (mth.ge.46.and.mth.le.49) go to 180
   do j=2,ng2
      niloc=lz+1+nl*nz*(j-1)
      jg2=(ngps+1)-(ig2lo+j-2)
      noloc=irsize+8*ngps+jg1
      a(noloc)=a(noloc)+e(niloc)
   enddo
   ! transport xsec (total-p1 scattering)
  180 continue
   if (nl.lt.2) go to 200
   do j=2,ng2
      niloc=lz+1+nl*nz*(j-1)+1
      noloc=irsize+3*ngps+jg1
      a(noloc)=a(noloc)-e(niloc)
   enddo

   !--transport shielding factor
  200 continue
  if (nl.lt.2) go to 300
   inc=4
   if (nsblk.gt.1) inc=1
   if (inc+nsb.ne.5) go to 300
   if (abtm.gt.eps*temp) go to 250
   ! finite dilution for non-self-shielded reactions
   if (nz.gt.1.and.mth.eq.2) go to 250
   if (mth.eq.2) ntfl=1
   do j=1,nzj
      noloc=l20-l18+(inc-1)*ngps*nzj*ntj+(jg1-1)*nzj*ntj-nzj+j
      do k=2,ng2
         jg2=(ngps+1)-(ig2lo+k-2)
         locsg=lz+nl*nz*(k-1)+2
         do l=1,ntj
            a(nzj*l+noloc)=a(nzj*l+noloc)-e(locsg)
         enddo
      enddo
   enddo
   go to 300
   ! finite dilution for elastic scattering
  250 continue
   if (mth.ne.2.or.ntfl.eq.1.or.nz.eq.1) go to 300
   do 260 j=1,ntj
   itemp=j
   if (abs(temp-ctem(j)).le.eps*temp) go to 270
  260 continue
   go to 300
  270 continue
   if (mth.eq.2) mt2tem(itemp)=1
   do j=1,nzj
      noloc=l20-l18+(inc-1)*ngps*nzj*ntj+(jg1-1)*nzj*ntj&
        +(itemp-1)*nzj+j
      locsg=lz-nl*nz+nl*(isig(j)-1)+2
      do k=2,ng2
         a(noloc)=a(noloc)-e(nl*nz*k+locsg)
      enddo
   enddo

   !--elastic removal shielding factor
  300 continue
   if (mth.ne.2) go to 140
   inc=6
   if (nsblk.gt.1) inc=1
   if (inc+nsb.ne.7) go to 140
   do 310 j=1,ntj
   itemp=j
   if (abs(temp-ctem(j)).le.eps*temp) go to 320
  310 continue
   go to 140
  320 continue
   do j=1,nzj
      do k=2,ng2
         jg2=(ngps+1)-(ig2lo+k-2)
         if (jg2.ne.jg1) then
            locsg=lz+nl*nz*(k-1)+(isig(j)-1)*nl+1
            noloc=l20-l18+(inc-1)*ngps*nzj*ntj+(jg1-1)*nzj*ntj&
              +(itemp-1)*nzj+j
            lrem=irsize+11*ngps+jg1
            a(noloc)=a(noloc)+e(locsg)/a(lrem)
         endif
      enddo
   enddo

  140 continue
   if (ig.lt.ngps) go to 120
   return
   end subroutine mxproc

   subroutine thnwrt(a,irec3)
   !--------------------------------------------------------------------
   ! Thins data and writes on nscrt1 in proper order.
   !--------------------------------------------------------------------
   use util ! provides repoz
   ! externals
   integer::irec3
   real(4)::a(*)
   ! internals
   integer::irec2,nrnge,i,nlow,ndex,inv,ir,jg,jt,jz,mdex,nup,j

   !--read self-shielding factors
   irec2=0
   call repoz(-nscrt2)
   nrnge=6*ngps
   irec2=irec2+1
   read(nscrt2)(a(i),i=1,nrnge)
   do i=1,nsblk
      nlow=nrnge+ngps*ntj*nzj*(i-1)+1
      irec2=irec2+1
      read(nscrt2)(a(nlow+j-1),j=1,lrsize)

      !--thin data (to (jbh-jbl+1)*ntap*ntat)
      ndex=nrnge
      inv=nreact/nsblk
      do ir=1,inv
         do jg=1,ngps
            if (jg.ge.jbl.and.jg.le.jbh) then
               do jt=1,ntj
                  if (ntat(jt).ne.0) then
                     do jz=1,nzj
                        if (ntap(jz).ne.0) then
                           mdex=(ir-1)*ngps*ntj*nzj&
                             +(jg-1)*ntj*nzj+(jt-1)*nzj+jz+nrnge
                           ndex=ndex+1
                           !--set s-s factors equal to 1.0
                           !--if reaction is not present
                           if (a(mdex).eq.0.0) a(mdex)=1
                           a(ndex)=a(mdex)
                        endif
                     enddo
                  endif
               enddo
            endif
         enddo
      enddo
      !--write on nscrt3
      nlow=nrnge+1
      nup=ndex-nrnge
      irec3=irec3+1
      write(nscrt3)(a(nlow+j-1),j=1,nup)
   enddo

   !--write cross sections to nscrt3
   irec3=irec3+1
   write(nscrt3)(a(i),i=1,nrnge)
   return
   end subroutine thnwrt

   subroutine pbrkxs(n4)
   !--------------------------------------------------------------------
   ! Controls printing of the BRKOXS file.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use util ! provides error,repoz
   ! externals
   integer::n4
   ! internals
   integer::irec1,i,nwds,iver,ngroup,nisosh,nsigpt,ntempt,nreact
   integer::iblk,next,nex1,ntp,j,ntt,nex2,nex3,nex4,jssf,k,nbint
   integer::nbtem,jbfli,jbfhi,ndiff,npro,j2,j3,j4,j5,j6,n
   integer::nwh,llr,nwr,lli,nwi
   character(8)::hn,hu,hs
   integer::ird(10)
   character(6),parameter::hgroup='index'
   character(6),parameter::htemp='temp'
   character(6),parameter::hsigo='sigo'

   irec1=0
   call repoz(-n4)
   do i=1,10
      ird(i)=1
   enddo

   !--file identification
   nwds=1+3*mult
   irec1=irec1+1
   read(n4)(ha(i),i=1,3),ia(nwds)
   hn=ta(1)
   hu=ta(2)
   hs=ta(3)
   iver=ia(nwds)
   write(nsyso,'(''1''/'' ***file '',a6,'' -- version'',i3,&
     &'' -- unit'',i3,''***'')') hn,iver,n4
   write(nsyso,'('' **user identification**'',a6,3x,a6)') hu,hs

   !--file control
   irec1=irec1+1
   read(n4)(ia(i),i=1,6)
   ngroup=ia(1)
   nisosh=ia(2)
   nsigpt=ia(3)
   ntempt=ia(4)
   nreact=ia(5)
   iblk=ia(6)
   if (ird(1).eq.1) write(nsyso,'(/&
     &'' file control parameters''//&
     &''  ngroup    number of energy groups in set         '',i12/&
     &''  nisosh    number of isotopes with self-''/&
     &''              shielding factors                    '',i12/&
     &''  nsigpt    total number of values of variable x''/&
     &14x,''which are given.  nsigpt is equal to''/&
     &''              the sum from 1 to nisosh of ntabp(i) '',i12/&
     &''  ntempt    total number of values of variable tb''/&
     &14x,''which are given.  ntempt is equal to''/&
     &''            the sum from 1 to nisosh of ntabt(i)   '',i12/&
     &''  nreact    number of reactions                    '',i12/&
     &''  iblk      blocking option                        '',i12)')&
     (ia(i),i=1,6)

   !--file data
   nwh=nisosh
   llr=mult*nwh
   nwr=nsigpt+ntempt+ngroup+1
   lli=llr+nwr
   nwi=4*nisosh
   nwds=mult*nwh+nwr+nwi
   irec1=irec1+1
   if (nwds.ge.isiza4)&
     call error('pbrkxs','input record too large.',' ')
   read(n4)(ha(i),i=1,nwh),(a(llr+i),i=1,nwr),(ia(lli+i),i=1,nwi)
   if (ird(2).eq.1) then
      write(nsyso,'(/4x,''isotope'',4x,''name'')')
      write(nsyso,'(5x,i3,7x,a6)') (i,ta(i),i=1,nisosh)
      next=nisosh*mult+1
      write(nsyso,'(/&
        &'' ln(sigp0)/ln(10) values for all isotopes''/5x,&
        &''isotope'',6x,''1st value'',6x,''2nd value'',&
        &10x,''. . .'')')
      nex1=(mult+2)*nisosh+nsigpt+ntempt+ngroup+1
      do i=1,nisosh
         nex1=nex1+1
         ntp=next+ia(nex1)-1
         write(nsyso,'(7x,i3,1p,8e15.5)') i,(a(j),j=next,ntp)
         next=ntp+1
      enddo
      write(nsyso,'(/&
        &'' temperatures (deg c) for all isotopes''//'' isotope'',&
        &6x,''1st value'',6x,''2nd value'',10x,''. . .'')')
      do i=1,nisosh
         nex1=nex1+1
         ntt=next+ia(nex1)-1
         write(nsyso,'(7x,i3,1p,8e15.5)') i,(a(j),j=next,ntt)
         next=ntt+1
      enddo
      write(nsyso,'(/'' maximum energy bound''/3x,''group j'')')
      do i=1,ngroup
         write(nsyso,'(i10,1p,e13.5)') i,a(next+i-1)
      enddo
      next=next+ngroup
      write(nsyso,'(/&
        &'' minimum energy bound of set''//11x,1p,e12.5)') a(next)
      write(nsyso,'(/&
        &'' f-factor start and stop groups''/&
        &'' and number of sig0 and temperature values''//&
        &1x,''isotope'',11x,''jbfh'',11x,''jbfl'',10x,''ntabf'',&
        &10x,''ntabt'')')
      write(nsyso,'(1x,i7,4i15)')&
        (i,ia(next+i),ia(next+i+nisosh),ia(i+next+2*nisosh),&
        ia(next+i+3*nisosh),i=1,nisosh)
   endif
   next=nisosh*mult+nsigpt+ntempt+ngroup+1
   nex1=next+1
   nex2=nex1+nisosh
   nex3=nex2+nisosh
   nex4=nex3+nisosh

   !--self-shielding factors
   jssf=nex4+nisosh
   do 110 i=1,nisosh
   k=nex3-1+i
   nbint=ia(k)
   k=nex4-1+i
   nbtem=ia(k)
   k=nex1-1+i
   jbfli=ia(k)
   k=nex2-1+i
   jbfhi=ia(k)
   ndiff=jbfhi-jbfli+1
   npro=nbint*nbtem*ndiff
   if (iblk.gt.1) go to 82
   nwds=nreact*npro
   if (nwds.le.0) go to 100
   irec1=irec1+1
   if (nwds.ge.isiza4)&
     call error('pbrkxs','input record too large.',' ')
   read(n4)(a(jssf+j-1),j=1,nwds)
   if (ird(3).ne.1) go to 100
   j2=jssf+npro
   j3=j2+npro
   j4=j3+npro
   j5=j4+npro
   j6=j5+npro
   go to 83
   82 nwds=npro
   if (nwds.le.0) go to 100
   j2=jssf
   j3=j2
   j4=j3
   j5=j4
   j6=j5
   irec1=irec1+1
   read(n4)(a(jssf+j-1),j=1,nwds)
   83 if (ird(3).eq.1) then
   write(nsyso,'(/&
        &'' total self-shielding factors'',10x,''isotope'',i3)') i
      call wot(a(jssf),nbint,nbtem,ndiff,hsigo,htemp,hgroup,nsyso)
   endif
   if (iblk.ne.0) then
      irec1=irec1+1
      read(n4)(a(j2+j-1),j=1,nwds)
   endif
   if (ird(3).eq.1) then
      write(nsyso,'(/&
        &'' capture self-shielding factors'',10x,''isotope'',i3)') i
      call wot(a(j2),nbint,nbtem,ndiff,hsigo,htemp,hgroup,nsyso)
   endif
   if (iblk.ne.0) then
      irec1=irec1+1
      read(n4)(a(j3+j-1),j=1,nwds)
   endif
   if (ird(3).eq.1) then
      write(nsyso,'(/&
        &'' fission self-shielding factors'',10x,''isotope'',i3)') i
      call wot(a(j3),nbint,nbtem,ndiff,hsigo,htemp,hgroup,nsyso)
   endif
   if (iblk.ne.0) then
      irec1=irec1+1
      read(n4)(a(j4+j-1),j=1,nwds)
   endif
   if (ird(3).eq.1) then
      write(nsyso,'(/&
        &'' transport self-shielding factors'',10x,''isotope'',i3)') i
      call wot(a(j4),nbint,nbtem,ndiff,hsigo,htemp,hgroup,nsyso)
   endif
   if (iblk.ne.0) then
      irec1=irec1+1
      read(n4)(a(j5+j-1),j=1,nwds)
   endif
   if (ird(3).eq.1) then
      write(nsyso,'(/&
        &'' elastic self-shielding factors'',10x,''isotope'',i3)') i
      call wot(a(j5),nbint,nbtem,ndiff,hsigo,htemp,hgroup,nsyso)
   endif
   if (nreact.lt.6) go to 100
   if (iblk.ne.0) then
      irec1=irec1+1
      read(n4)(a(j6+j-1),j=1,nwds)
   endif
   if (ird(3).eq.1) then
      write(nsyso,'(/&
        &'' removal self-shielding factors'',10x,''isotope'',i3)') i
      call wot(a(j6),nbint,nbtem,ndiff,hsigo,htemp,hgroup,nsyso)
   endif
100   continue

   !--cross sections
   nwds=6*ngroup
   irec1=irec1+1
   if (nwds.ge.isiza4)&
     call error('pbrkxs','input record too large.',' ')
   read(n4)(a(jssf+n-1),n=1,nwds)
   if (ird(4).eq.1) then
      j2=jssf+ngroup
      j3=j2+ngroup
      j4=j3+ngroup
      j5=j4+ngroup
      j6=j5+ngroup
      write(nsyso,'(/&
        &'' group'',9x,''xspo'',9x,''xsin'',10x,''xse'',9x,''xsmu'',&
        &9x,''xsed'',9x,''xsxi'')')
      do n=1,ngroup
         write(nsyso,'(1x,i5,1p,6e13.6)')&
           n,(a(jssf-1+n+j*ngroup-ngroup),j=1,6)
      enddo
   endif

   !--continue loop over shielding factors
  110 continue
   return
   end subroutine pbrkxs

   subroutine cdlyxs(idlay,idlayt)
   !--------------------------------------------------------------------
   ! Write a CCCC-III DLAYXS interface file.
   !--------------------------------------------------------------------
   use util ! provides mess,repoz
   ! externals
   integer::idlay,idlayt
   ! internals
   integer::lzero,irec,nwds,lzeroh,l3,n,l4,nkfam
   integer::nwh,nwr,nwi,lli,l1h,llr,l1i,l1r,i
   character(6)::hname='dlayxs'

   !--read delayed neutron data
   lzero=next
   call dldata(idlayt)
   if (nisod.ne.0) then

      !--file identification
      irec=1
      nwds=1+3*mult
      lzeroh=(lzero-1)/mult+1
      read(hname,'(a6)') ta(lzeroh)
      ha(lzeroh+1)=huse(1)
      ha(lzeroh+2)=huse(2)
      ia(3*mult+lzero)=ivers
      call repoz(idlay)
      write(idlay)(ha(lzeroh+i-1),i=1,3),ia(lzeroh+nwds-1)

      !--file control
      irec=2
      nwds=4
      ia(lzero)=ngps
      ia(lzero+1)=nisod
      ia(lzero+2)=nfam
      ia(lzero+3)=nfam
      write(idlay)(ia(lzero+i-1),i=1,4)

      !--file data, decay constants, and emission spectra
      irec=3
      l1h=(l1-1)/mult
      nwh=nisod
      l1r=l1+mult*nwh-1
      nwr=(ngps+1)*(nfam+1)
      l1i=l1r+nwr
      nwi=2*nisod
      nwds=mult*nwh+nwr+nwi
      write(idlay)(ha(l1h+i),i=1,nwh),(a(l1r+i),i=1,nwr),&
        (ia(l1i+i),i=1,nwi)

      !--delayed neutron precursor yield data by isotope
      l3=l1+nwds
      do n=1,nisod
         irec=4
         l4=l1+nisod*mult+nfam+ngps*nfam+ngps
         nkfam=ia(l4+n)
         llr=l3-1
         nwr=ngps*nkfam
         lli=llr+nwr
         nwi=nkfam
         nwds=nwr+nwi
         write(idlay)(a(llr+i),i=1,nwr),(ia(lli+i),i=1,nwi)
         l3=l3+nwds
      enddo
      next=l1

   !--no delayed neutron data found
   else
      call mess('cdlyxs','no delayed neutron data found.',' ')
      next=l1
   endif
   return
   end subroutine cdlyxs

   subroutine dldata(nin)
   !--------------------------------------------------------------------
   ! Read delayed neutron data from input tape.
   ! Store all data in a according to the following assignments:
   !
   ! location  variable   length
   ! --------  --------   ------
   !   l1      isonm      nisod
   !   l2      flam       nfam
   !   l3      chid       ngroup*nfam
   !   l4      emax       ngroup
   !   l5      emin       1
   !   l6      nfami      nisod
   !   l7      loca       nisod
   !   l8      snudel     ngroup*nfam
   !            +numfam     +nfam
   !   l9      next available location
   !
   ! In this routine, "family" means an isotope and time group
   ! (there are typically nisod*ndg families, where ndg is the
   ! number of delayed neutron groups for this isotope).
   !--------------------------------------------------------------------
   use util ! provides repox,error
   use endf ! provides endf routines and variables
   ! externals
   integer::nin
   ! internals
   integer::nisoo,nwds,lzero,l3,l4,l5,l6,l7,l8,l9,la,last
   integer::i,ipr,nb,nw,ig,ig1,nl,nz,ng,ng2,loc,ngn,ngp1,loca
   integer::ifam,nd,jj,locb,j,ndg
   integer::imusd(200)
   integer,parameter::ndmax=8
   real(kr)::fract(ndmax)
   character(8)::hisnm
   nisoo=niso

   !--get the number of delayed neutron groups for this nuclide
   !--from groupr's mf5, mt455 head record.
   call repoz(nin)
   call tpidio(nin,0,0,e(1),nb,nw)
   do while (mfh.lt.5)
      call contio(nin,0,0,e(1),nb,nw)
   enddo
   if (mth.eq.455) then
      ndg=nint(e(3))
   else if (mth.lt.455) then
      do while (mfh.eq.5.and.mth.lt.455)
         call contio(nin,0,0,e(1),nb,nw)
      enddo
      if (mfh.eq.5) then
         ! number of delayed neutron groups
         ndg=nint(e(3))
      else
         nisod=0
         return
      endif
   else
      nisod=0
      return
   endif
   if (ndg.eq.0) then
       nisod=0
       return
   elseif (ndg.gt.ndmax) then
       call mess ('dldata','too many delayed neutron groups',&
         'dlayxs request ignored')
       nisod=0
       return
   endif

   !--assign storage
   nfam=ndg*nisoo
   nwds=3*mult+mult
   lzero=next
   next=lzero+nwds+4
   nwds=nisoo*mult
   l1=next
   next=l1+nwds
   l2=next
   next=l2+nfam
   nwds=ngps*nfam
   l3=next
   next=l3+nwds
   l4=next
   next=l4+ngps
   l5=next
   next=l5+1
   l6=next
   next=l6+nisoo
   l7=next
   next=l7+nisoo
   nwds=ngps*nfam+nfam
   l8=next
   next=l8+nwds
   l9=next
   la=l9+5
   last=la+10*ngps

   if (last.gt.isiza4) call error('dldata',&
        'max4a for container a exceeded',' ')

   !--read through input tape for desired data
   do i=1,200
      imusd(i)=0
   enddo
   ipr=0
   nisod=0
   do i=l1,last
      a(i)=0
   enddo
   call repoz(nin)
   call tpidio(nin,0,0,e(1),nb,nw)
  110 call contio(nin,0,0,e(1),nb,nw)
   if (math.eq.-1) return
   if (mfh.eq.0) go to 110
   do 115 i=1,nisoo
   read(hisonm(i),'(a6)') hisnm
   if (math.ne.imat(i)) go to 115
   if (imusd(i).gt.0) go to 115
   imusd(i)=1
   go to 150
  115 continue
   call tomend(nin,0,0,e(1))
   go to 110
  130 if (ig.lt.ngps) go to 170
   call contio(nin,0,0,e(1),nb,nw)
   if (math.eq.0) go to 110
   if (mfh.eq.0.or.mth.eq.0) go to 130
  150 nl=l1h
   nz=l2h
   ng=n2h
  170 call listio(nin,0,0,e(1),nb,nw)
   ng2=l1h
   ig=n2h
   loc=1+nw
   do while (nb.ne.0)
      if (loc+302.gt.maxe) call error('dldata',&
        'max for endf input data exceeded',' ')
      call moreio(nin,0,0,e(loc),nb,nw)
      loc=loc+nw
   enddo
   if (mth.eq.451) go to 210
   if (mfh.eq.3.and.mth.eq.455) go to 310
   if (mfh.eq.5.and.mth.eq.455) go to 410
   go to 130

   !--process material header
  210 continue
   nisod=nisod+1
   if (ipr.gt.0) go to 130
   ngn=ng2
   if (ngn.ne.ngps)&
     call error('dldata','incompatible group structures.',' ')
   ngp1=ngn+1
   do i=1,ngp1
      loca=l4+ngn-i+1
      a(loca)=e(6+ng+nz+i)
   enddo
   ipr=1
   ig=ngps
   go to 130

   !--process delayed neutron yield records
  310 continue
   do i=1,ndg
      ifam=ndg*(nisod-1)+i
      loca=l8+ngn-ig+ngn*(ifam-1)+ndg*(nisod-1)
      a(loca)=e(8)
   enddo
   go to 130

   !--process delayed neutron spectra record
  410 nd=nl
   ! lowest energy group that holds delayed data
   ! subtraction of -1 because first position in gendf mf5mt455
   ! is decay constant of that precursor family
   ig1=l1h-1
   ! hisnm--isotope name
   jj=(l1-1)/mult+nisod
   ta(jj)=hisnm
   ! loca--number of records before desired isotope.
   ia(l7+nisod-1)=nisod-1
   ! nkfam--number of families for this isotope.
   loca=l6+nisod-1
   ia(loca)=nd
   ! flam--decay constants for these families.
   ifam=ndg*(nisod-1)-1+l2
   do i=1,nd
      a(i+ifam)=e(i+6)
   enddo
   ! chid--delayed neutron spectra for these families.
   do i=1,nd
      ifam=ndg*(nisod-1)+i
      loca=l3+ngn+ngn*(ifam-1)
      locb=6+i
      fract(i)=0
      do j=1,ig1
         a(loca-j)=e(nl*j+locb)
         fract(i)=e(nl*j+locb)+fract(i)
      enddo
   enddo
   ! snudel--yield of delayed neutrons for each family.
   do i=1,nd
      ifam=i+ndg*(nisod-1)
      loca=l8+ngn+ngn*(ifam-1)+ndg*(nisod-1)
   !--numfam--family number for each yield vector
      locb=l8-1+ndg*ngn*nisod+ndg*(nisod-1)+i
      ia(locb)=ifam
      do j=1,ngn
         a(loca-j)=fract(i)*a(loca-j)
      enddo
   enddo
   go to 130
   end subroutine dldata

   subroutine pdlyxs(idlay)
   !--------------------------------------------------------------------
   ! Print the CCCC DLAYXS file
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use util ! provides error,repoz
   ! externals
   integer::idlay
   ! internals
   integer::n4,irec1,nwds,i4,ngroup,nisod,nfam,i,nsp,nint
   integer::k,nex1,nex2,jwd1,kk,jwd2,nwh,llr,lli,j
   character(6),parameter::h1=' group'
   character(6),parameter::h2='family'
   character(6),parameter::h3='isotop'
   character(6),parameter::h4='number'
   character(6),parameter::h5='vector'
   character(6),parameter::hblank='      '
   integer::m1=1
   integer::m4=4
   n4=idlay

   !--file identification
   irec1=1
   nwds=1+3*mult
   call repoz(n4)
   read(n4)(ha(i),i=1,3),ia(nwds)
   i4=1+3*mult
   write(nsyso,'(&
     &''1delayed neutron precursor data '',a6/&
     &'' user identification'',1x,2a6/&
     &'' file version number'',i6)')&
     ta(1),ta(2),ta(3),ia(i4)

   !--file control
   irec1=irec1+1
   read(n4)(ia(i),i=1,m4)
   ngroup=ia(1)
   nisod=ia(2)
   nfam=ia(3)
   write(nsyso,'(/&
     &'' file control parameters''//&
     &''  ngroup    number of neutron energy groups in set '',i12/&
     &''  nisod     number of isotopes in delayed          '',i12/&
     &''              neutron set                          ''/&
     &''  nfam      number of delayed neutron families     '',i12/&
     &''              in set                               ''/&
     &''  idum      dummy to make up 4-word record         '',i12)')&
     (ia(i),i=1,4)

   !--file data
   nsp=(nfam+1)*(ngroup+1)
   nint=2*nisod
   nwh=nisod
   llr=mult*nwh
   lli=llr+nsp
   nwds=nwh*mult+nsp+nint
   irec1=irec1+1
   if (nwds.ge.isiza4)&
     call error('pdlyxs','input record too large.',' ')
   read(n4)(ha(i),i=1,nwh),(a(llr+i),i=1,nsp),(ia(lli+i),i=1,nint)
   write(nsyso,'(/3x,''isotope'',5x,''name'')')
   next=nisod
   do i=1,next
      write(nsyso,'(5x,i3,7x,a6)') i,ta(i)
   enddo
   next=mult*nisod
   write(nsyso,'(/&
     &'' delayed neutron precursor decay constant for family n'',&
     &/6x,''family'',9x,''n'')')
   do i=1,nfam
      k=next+i
      write(nsyso,'(6x,i3,6x,1p,e12.5)') i,a(k)
   enddo
   next=next+nfam+1
   write(nsyso,'(/&
     &'' fraction of delayed neutrons emitted into neutron '',&
     &''energy group from precursor family'')')
   call wot(a(next),ngroup,nfam,m1,h1,h2,hblank,nsyso)
   next=next+ngroup*nfam
   write(nsyso,'(/&
     &'' maximum energy bound''/7x,''group'',5x,''j'')')
   do i=1,ngroup
      write(nsyso,'(i10,1p,e13.5)') i,a(next+i-1)
   enddo
   next=next+ngroup
   write(nsyso,'(/&
     &'' minimum energy bound of set''//11x,1p,e12.5)')&
     a(next)
   nex1=next+1
   write(nsyso,'(/&
     &'' number of families to which fission in isotope '',&
     &''contributes delayed neutron precursors'')')
   call woti(ia(nex1),nisod,m1,m1,h3,h4,hblank,nsyso)
   nex2= nex1+nisod
   write(nsyso,'(/&
     &'' number of records to be skipped to read data for '',&
     &''isotope'')')
   call woti(ia(nex2),nisod,m1,m1,h3,h4,hblank,nsyso)
   next=mult*nisod+(ngroup+1)*(nfam+1)

   !--delayed neutron procursor yield data
   jwd1=nwds+1
   do i=1,nisod
      k=next+i
      kk=ia(k)
      nsp=ngroup*kk
      nint=kk
      nwds=nsp+nint
      llr=jwd1-1
      lli=llr+nsp
      irec1=irec1+1
      if (nwds.ge.isiza4)&
        call error('pdlyxs','input record too large.',' ')
      read(n4)(a(llr+j),j=1,nsp),(ia(lli+j),j=1,nint)
      jwd2=jwd1+ngroup*kk
      write(nsyso,'(/&
        &''number of delayed neutron precursors produced in '',&
        &''family per fission in group'')')
      call wot(a(jwd1),ngroup,kk,m1,h1,h2,hblank,nsyso)
      write(nsyso,'(/&
        &'' family number of k-th yield vector'')')
      call woti(ia(jwd2),m1,kk,m1,h2,h5,hblank,nsyso)
   enddo
   return
   end subroutine pdlyxs

   subroutine wot(x,ii,jj,kk,top1,top2,top3,nou)
   !--------------------------------------------------------------------
   ! Prints a single one-, two-, or three-dimensional real array.
   !--------------------------------------------------------------------
   ! externals
   integer::ii,jj,kk,nou
   real(4)::x(ii,jj,kk)
   character(6)::top1,top2,top3
   ! internals
   integer::k,jo1,jo2,jo3,jk,j,nflag,i,n,im1
   real(kr)::ssum
   real(kr),parameter::zero=0

   do 110 k=1,kk
   jo2=0
   jo3=(jj+8)/9
   if (kk.gt.1) then
      write(nou,'('' '',6x,a6,i5)') top3,k
   endif
   do 100 jk=1,jo3
   jo1=jo2+1
   jo2=min0(jo1+8,jj)
   write(nou,'('' '',6x,a6,9(4x,a6,i3))')&
     top1,(top2,j,j=jo1,jo2)
   nflag=0
   do 90 i=1,ii
   ssum=0
   do j=jo1,jo2
      ssum=ssum+abs(x(i,j,k))
   enddo
   if (ssum.ne.zero) go to 30
   nflag=nflag+1
   if (i.lt.ii) go to 90
  30 if (nflag.eq.1) go to 60
   if (nflag.lt.1) go to 80
   n=i-nflag
   im1=i-1
   if (i.ge.ii.and.ssum.eq.zero) then
   n=n+1
   im1=i
   endif
   write(nou,'('' '',i7,'' to'',i4,1p,9e13.5)')&
     n,im1,(x(im1,j,k),j=jo1,jo2)
   go to 70
   60 if (ssum.le.zero) go to 80
   im1=i-1
   write(nou,'('' '',i10,4x,1p,9e13.5)')&
     im1,(x(im1,j,k),j=jo1,jo2)
   70 continue
   nflag=0
   if (ssum.le.zero) go to 90
   80 write(nou,'('' '',i10,4x,1p,9e13.5)')&
        i,(x(i,j,k),j=jo1,jo2)
   90 continue
  100 continue
  110 continue
   return
   end subroutine wot

   subroutine woti(ix,ii,jj,kk,top1,top2,top3,nou)
   !--------------------------------------------------------------------
   ! Prints a single one-, two-, or three-dimensional integer array.
   !--------------------------------------------------------------------
   ! externals
   integer::ii,jj,kk,nou
   integer(4)::ix(ii,jj,kk)
   character(6)::top1,top2,top3
   ! internals
   integer::i,j,k,jk,jo1,jo2,jo3

   do k=1,kk
      jo2=0
      jo3=(jj+24)/25
      if (kk.gt.1) then
         write(nou,'(/''0'',6x,a6,i5)') top3,k
      endif
      do jk=1,jo3
         jo1=jo2+1
         jo2=min0(jo1+24,jj)
         write(nou,'('' '',a6,''/'',a6,25i4)')&
           top1,top2,(j,j=jo1,jo2)
         do i=1,ii
            write(nou,'('' '',i6,7x,25i4)')&
              i,(ix(i,j,k),j=jo1,jo2)
         enddo
      enddo
   enddo
   return
   end subroutine woti

   subroutine stow(a,b,nwds)
   !--------------------------------------------------------------------
   ! Move data from one region of core to another.
   !--------------------------------------------------------------------
   ! externals
   integer::nwds
   real(4)::a(*),b(*)
   ! internals
   integer::n

   if (nwds.le.0) return
   do n=1,nwds
      b(n)=a(n)
   enddo
   return
   end subroutine stow

end module ccccm

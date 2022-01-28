module dtfm
   ! Module to provide dtfr for NJOY2016
   use locale
   implicit none
   private
   public dtfr

   ! global variables for dtfr

   ! units
   integer::nout,npend,nin,nplot

   ! neutron group structure
   real(kr)::egn(2501)
   integer::ngn

   ! gamma group structure
   real(kr)::egg(400)
   integer::ngp

   ! dtfr parameters
   integer::iprint,ifilm,iedit,nlmax,ng,iptotl,ipingp
   integer::itabl,ned,nptabl,jgthrm
   integer,parameter::nedmax=50
   character(6)::hednam(3+nedmax)
   integer::jped(3+nedmax),mted(3+nedmax)
   integer::multed(3+nedmax),ids(3+nedmax)
   character(6)::hisnam
   integer::matd,jz
   real(kr)::dtemp
   integer::kpos,ked,lpos,led,iphph
   integer::mti,mtc,nlc
   real(kr)::bot,axl,axd,top
   integer::ixo,iyo
   integer::nplt
   real(kr)::xmin,xmax,ymin,ymax

   ! sigma-zero storage
   real(kr)::sigz(10)
   integer::nsigz

   ! storage arrays
   integer,parameter::ndim=7000
   real(kr)::x(ndim),y(ndim),z(ndim)
   integer,parameter::nwamax=40000
   real(kr)::a(nwamax)
   integer,parameter::nwsmax=500000
   real(kr)::sig(nwsmax)

contains

   subroutine dtfr
   !-------------------------------------------------------------------
   !
   ! convert output of groupr to dtf format.
   !
   ! Processes neutron and gamma production cross sections and
   ! matrices.  The neutron tables can have reduced table length.
   ! Up-scatter is allowed.  The absorption reaction is computed
   ! from the total cross section and total scattering.  Any edits
   ! can be produced which are either given in the ENDF file
   ! or are linear combinations of ENDF cross sections.  The
   ! fission nu*sigf and chi are computed from the fission matrices
   ! for all partial fission reactions.  Chi includes source
   ! weighting.  The Pl tables for l.gt.0 contain the Pl weighted
   ! total in the total position and the pl transport cross section
   ! in the absorption position.  The gamma tables have gamma group
   ! 1 in position 1, 2 in position 2, etc, with a table length
   ! equal to the number of gamma groups.
   !
   ! Plots can be prepared in viewr format.
   !
   !---input specifications (free format)---------------------------
   !
   ! card 1       units
   !    nin       input unit with data from groupr (binary).
   !    nout      output unit containing dtf tables (coded).
   !              (default=0=none)
   !    npend     input unit with pendf tape for point plots.
   !              (default=0=none)
   !    nplot     output plot info for plotr module
   !              (default=0=none)
   ! card 2       options
   !    iprint    print control (0 minimum, 1 maximum)
   !    ifilm     film control (0/1/2=no/yes with 1 plot per frame/
   !              yes with 4 plots per frame (default=0)
   !    iedit     edit control (0/1=in table/separate) (default=0)
   !
   !       cards 3 through 5 only for iedit=0
   !
   ! card 3       neutron tables
   !    nlmax     number of neutron tables desired.
   !    ng        number of neutron groups
   !    iptotl    position of total cross section
   !    ipingp    position of in-group scattering cross section.
   !    itabl     neutron table length desired.
   !    ned       number of entries in edit table (default=0).
   !    ntherm    number of thermal groups (default=0).
   !  card 3a only for ntherm ne 0
   ! card 3a      thermal incoherent and coherent mts
   !    mti       mt for thermal incoherent data
   !    mtc       mt for thermal coherent data (default=0)
   !    nlc       no. coherent legendre orders (default=0)
   ! card 4       edit names
   !       six character hollerith names for edits for as many
   !       cards as needed.  there will be iptotl-3 names read.
   !       each name is delimited with *.
   ! card 5       edit specifications
   !       ned triplets of numbers on as many cards as needed.
   !       positions can appear more than once.
   !       reaction types can appear more than once.
   !    jpos      position of edit quantity.
   !    mt        endf reaction number.
   !    mult      multiplicity to be used when adding this mt.
   !
   !       card 6 for iedit=1
   !
   ! card 6       claw-format tables
   !    nlmax     number of neutron tables (def=5)
   !    ng        number of neutron groups (def=30)
   !              (number of thermal groups is zero)
   !
   ! card 7       gamma ray tables
   !    nptabl    number of gamma tables desired (default=0)
   !    ngp       number of gamma groups (default=0)
   ! card 8       material description
   !       one card for each table set desired.
   !       empty card (/) terminates execution of dtfr.
   !    hisnam    6-character isotope name
   !    mat       material number as in endf (default=0)
   !    jsigz     index number of sigma-zero desired (default=1)
   !    dtemp     temperature desired (default=300)
   !
   !-------------------------------------------------------------------
   use util ! provides timer,openz,repoz,error,closz
   use endf ! provides endf routines and variables
   use mainio ! provides nsysi,nsyso
   ! internals
   integer::nscr,nb,nw,i,nl,nz,ng2,ig2lo,lz,nwa
   integer::ntw,ngg,jbase,j,ngnp1,nggp1,il,ip,i618,ig,nws,n3
   integer::jg,loca,locs,ied,jpos,k,jg2,locb,locf,kg,id,jgp
   integer::jzd,ilmax,ipmax,loc,idone
   real(kr)::time,csz,temp,test,t,ff,sss,cnorm,cnm,dnorm
   real(kr)::spect(2501),fcap(2501),ffis(2501)
   real(kr),parameter::size=.30e0_kr
   real(kr),parameter::zero=0

   !--read  user input and write heading
   call timer(time)
   write(nsyso,'(/&
     &'' dtfr...produce dtf format from groupr output'',&
     &24x,f8.1,''s'')') time
   write(nsyse,'(/'' dtfr...'',61x,f8.1,''s'')') time
   call ruin
   nscr=11
   if (nin.lt.0) nscr=-nscr
   call openz(nscr,1)

   !--initialize plotting
   csz=size
   if (ifilm.gt.1) csz=csz/2
   write(nplot,'(''1 2'',f6.3,'' /'')') csz

   !--select next material, dilution, and temperature.
  105 continue
   hisnam=' '
   matd=0
   jzd=1
   dtemp=300
   read(nsysi,*) hisnam,matd,jzd,dtemp
   cnorm=0
   cnm=0
   dnorm=0
   ilmax=0
   ipmax=0
   ixo=0
   iyo=0
   if (matd.eq.0) go to 900
   write(nsyso,'(//&
     &'' mat='',i4,''  iso='',a6,''  sigzero no='',i2,&
     &''  temp='',1p,e10.3)') matd,hisnam,jzd,dtemp

   !--find desired material on input tape, process header section,
   !--and copy rest of material to scratch tape
   call repoz(nin)
   call tpidio(nin,0,0,a,nb,nw)
   call repoz(nscr)
   nsc=0
  130 continue
   call contio(nin,0,0,a,nb,nw)
   if (math.eq.matd) go to 135
   if (math.eq.-1) go to 132
   if (mfh.eq.0) go to 130
   call tomend(nin,0,0,a)
   go to 130
  132 continue
   write(nsyso,'(/'' ***material not found'')')
   go to 105
  135 continue
   loc=7
   call listio(nin,0,0,a(loc),nb,nw)
   temp=a(7)
   test=1
   test=test/1000
   if (abs(temp-dtemp).le.test) go to 136
   call tomend(nin,0,0,a)
   go to 130
  136 continue
   do while (nb.ne.0)
      loc=loc+nw
      call moreio(nin,0,0,a(loc),nb,nw)
   enddo
   nsigz=nint(a(4))
   ntw=nint(a(6))
   ngn=nint(a(9))
   if (ngn.ne.ng) call error('dtfr',&
     'number of neutron groups disagrees with number requested.',&
     ' ')
   ngg=nint(a(10))
   if (ngp.gt.0.and.ngg.ne.ngp) call error('dtfr',&
     'number of gamma groups disagrees with number requested.',&
     ' ')
   jbase=12
   jbase=jbase+ntw
   do j=1,nsigz
      sigz(j)=a(j+jbase)
   enddo
   jbase=jbase+nsigz
   ngnp1=ngn+1
   do j=1,ngnp1
      egn(j)=a(j+jbase)
   enddo
   jbase=jbase+ngnp1
   nggp1=ngg+1
   do j=1,nggp1
      egg(j)=a(j+jbase)
   enddo
   call tofend(nin,0,0,a)
   call tomend(nin,0,nscr,a)

   !--find desired material and temperature on pendf tape
   if (npend.eq.0) go to 129
   call repoz(npend)
   call tpidio(npend,0,0,a,nb,nw)
   call findf(matd,1,451,npend)
   test=dtemp/1000+1
   t=0
   idone=0
   do while (idone.eq.0)
      call contio(npend,0,0,a,nb,nw)
      if (math.ne.matd) then
         call error('dtfr','desired temperature not on pendf',' ')
      endif
      call contio(npend,0,0,a,nb,nw)
      if (n1h.ne.0) then
         iverf=4
      else if (n2h.eq.0) then
         iverf=5
      else
         iverf=6
      endif
      call skiprz(npend,-1)
      if (iverf.ge.5) call contio(npend,0,0,a,nb,nw)
      if (iverf.ge.6) call contio(npend,0,0,a,nb,nw)
      call hdatio(npend,0,0,a,nb,nw)
      t=a(1)
      if (abs(t-dtemp).ge.test) then
         call tomend(npend,0,0,a)
      else
         idone=1
      endif
   enddo

   !--generate requested neutron and photon tables.
  129 continue
   il=1
   if (nlmax.eq.0) il=0
   ip=0
   if (nlmax.eq.0) ip=1
   i618=0
  110 continue
   call repoz(nscr)
   ig=ngn
   if (il.gt.0) nws=itabl*ng
   if (ip.gt.0) nws=ngp*ng
   if (nws.gt.nwsmax)&
     call error('dtfr','not enough storage for table.',' ')
   n3=ned+3
   do i=1,nwsmax
      if (i.le.n3) ids(i)=0
      sig(i)=0
   enddo

   !--read next record.  process according to type.
  150 continue
   if (ig.lt.ngn) go to 160
   call contio(nscr,0,0,a,nb,nw)
   if (math.eq.0) go to 800
   if (mfh.eq.0.or.mth.eq.0) go to 150
   if (mth.eq.501) mted(ned+3)=501
   if (mfh.eq.23) iphph=1
   nl=nint(a(3))
   if (il.gt.nl.or.ip.gt.nl) go to 155
   nz=nint(a(4))
   jz=jzd
   if (jz.gt.nz) jz=1
   go to 160
  155 continue
   call tosend(nscr,0,0,a)
   go to 150
  160 continue
   call listio(nscr,0,0,a,nb,nw)
   ng2=nint(a(3))
   ig2lo=nint(a(4))
   ig=nint(a(6))
   lz=6
   nwa=1+nw
   do while (nb.ne.0)
      call moreio(nscr,0,0,a(nwa),nb,nw)
      nwa=nwa+nw
      if (nwa.gt.nwamax)&
        call error('dtfr','not enough storage for record.',' ')
   enddo
   if (il.eq.1.and.mth.eq.455) go to 400
   if (il.gt.0.and.mfh.eq.3) go to 200
   if (il.gt.0.and.mfh.eq.23) go to 200
   if (il.gt.0.and.mfh.eq.6) go to 300
   if (il.gt.0.and.mfh.eq.26) go to 300
   if (ip.gt.0.and.(mfh.eq.16.or.mfh.eq.17)) go to 600
   go to 150

   !--process neutron cross sections.
  200 continue
   if (mth.eq.1) go to 210
   if (mth.eq.501) go to 210
   if (mth.eq.102.or.mth.eq.18) go to 220
   if (ned.gt.0) go to 240
   go to 260
   ! save total.
  210 continue
   if (il.gt.nl) go to 150
   jg=ng-ig+1
   loca=lz+il+nl*((jz-1)+nz)
   locs=iptotl+itabl*(jg-1)
   ids(ned+3)=1
   sig(locs)=a(loca)
   if (il.gt.1) go to 150
   locs=locs-2
   sig(locs)=sig(locs)+a(loca)
   ids(ned+1)=1
   if (ned.eq.0) go to 150
   go to 240
   ! save capture and fission self-shielding factors
  220 continue
   jg=ng-ig+1
   loca=lz+il+nl*nz
   if (mth.eq.18) ffis(jg)=1
   if (mth.eq.18.and.a(loca).ne.zero)&
     ffis(jg)=a(loca+nl*(jz-1))/a(loca)
   if (mth.eq.102) fcap(jg)=1
   if (mth.eq.102.and.a(loca).ne.zero)&
     fcap(jg)=a(loca+nl*(jz-1))/a(loca)
   if (ned.gt.0) go to 240
   go to 260
   ! save edit cross sections
  240 continue
   if (il.gt.1) go to 150
   do ied=1,ned
      if (mth.eq.1.and.mted(ied).eq.300) then
         loca=lz+il+nl*(jz-1)
         jpos=jped(ied)
         locs=jpos+itabl*(jg-1)
         sig(locs)=a(loca)*multed(ied)+sig(locs)
         ids(ied)=ied
      else if (mth.eq.mted(ied)) then
         jg=ng-ig+1
         loca=lz+il+nl*((jz-1)+nz)
         locs=jped(ied)+itabl*(jg-1)
         sig(locs)=sig(locs)+multed(ied)*a(loca)
         ids(ied)=ied
      endif
   enddo
   ! thermal correction for total and absorption
  260 continue
   if (il.gt.1) go to 150
   if (mth.ne.2.and.mth.ne.mti.and.mth.ne.mtc) go to 150
   jg=ng-ig+1
   if (jg.lt.jgthrm) go to 150
   loca=lz+il+nl*((jz-1)+nz)
   locs=iptotl+itabl*(jg-1)
   if (mth.eq.2) sig(locs)=sig(locs)-a(loca)
   if (mth.eq.mti.or.mth.eq.mtc) sig(locs)=sig(locs)+a(loca)
   if (mth.eq.2) sig(locs-2)=sig(locs-2)-a(loca)
   if (mth.eq.mti.or.mth.eq.mtc) sig(locs-2)=sig(locs-2)+a(loca)
   go to 150

   !--accumulate neutron transfer matrices.
  300 continue
   if (mth.eq.18) go to 340
   if (mth.eq.19.or.mth.eq.20.or.mth.eq.21.or.mth.eq.38) go to 340
   ilmax=il
   jg=ng-ig+1
   if (mth.eq.2.and.jg.ge.jgthrm) go to 150
   if (mth.lt.200.or.mth.gt.250) go to 305
   if (jg.lt.jgthrm) go to 150
   if (mth.eq.mti) go to 305
   if (mth.eq.mtc.and.mth-mtc+1.eq.il) go to 305
   go to 150
  305 continue
   do k=2,ng2
      jg2=ng-ig2lo-k+3
      jpos=ipingp+(jg2-jg)
      if (jpos.gt.itabl) jg2=jg2-jpos+itabl
      if (jpos.gt.itabl) jpos=itabl
      if (jpos.le.iptotl) jg2=jg2+iptotl+1-jpos
      if (jpos.le.iptotl) jpos=iptotl+1
      if (jpos.le.itabl) then
         loca=lz+il+nl*((jz-1)+nz*(k-1))
         locs=jpos+itabl*(jg2-1)
         sig(locs)=sig(locs)+a(loca)
         if (il.le.1) then
            locb=iptotl-2+itabl*(jg-1)
            sig(locb)=sig(locb)-a(loca)
         endif
      endif
   enddo
   go to 150

   !--special branch for fission reactions.
   !--accumulate nu*sigf and chi.
  340 continue
   if (il.gt.1) go to 150
   if (mth.eq.18) i618=1
   if (mth.eq.19.and.i618.gt.0)&
     call error('dtfr','mt18 already processed, mt19 not allowed.',&
     ' ')
   ids(ned+2)=1
   if (ig.eq.0) go to 360
   if (ig2lo.eq.0) go to 380
   jg=ng-ig+1
   do k=2,ng2
      loca=lz+il+nl*((jz-1)+nz*(k-1))
      locs=iptotl-1+itabl*(jg-1)
      sss=a(loca)
      if (mth.eq.18.or.mth.eq.19) sss=sss*ffis(jg)
      sig(locs)=sig(locs)+sss
      if (kpos.gt.0) then
         locf=lz+il+nl*(jz-1)
         jg2=ng-ig2lo-k+3
         locs=kpos+itabl*(jg2-1)
         sig(locs)=sig(locs)+a(locf)*sss
         cnorm=cnorm+a(locf)*sss
         ids(ked)=ked
      endif
   enddo
   go to 390
  360 continue
   do k=1,ng
      jg=ng-k+1
      spect(jg)=0
      if (k.ge.ig2lo.and.k.le.(ig2lo+ng2-1)) then
         loca=lz+k-ig2lo+1
         spect(jg)=a(loca)
      endif
   enddo
   go to 150
  380 continue
   if (kpos.eq.0) go to 150
   jg=ng-ig+1
   loca=lz+il+nl*(jz-1)+jz
   locf=lz+il
   locs=iptotl-1+itabl*(jg-1)
   sss=a(loca)
   if (mth.eq.18.or.mth.eq.19) sss=sss*ffis(jg)
   sig(locs)=sig(locs)+sss
   cnm=cnm+sss*a(locf)
  390 continue
   if (ig.ge.ng.and.cnm.ne.zero) then
      do k=1,ng
         locs=kpos+itabl*(k-1)
         sig(locs)=sig(locs)+cnm*spect(k)
         cnorm=cnorm+cnm*spect(k)
         ids(ked)=ked
      enddo
   endif
   go to 150

   !--special branch for delayed neutrons.
  400 continue
   if (mfh.eq.5) go to 440
   jg=ng-ig+1
   locs=iptotl-1+itabl*(jg-1)
   loca=lz+1
   sig(locs)=sig(locs)+a(loca+1)*a(loca+2)*ffis(jg)
   ids(ned+2)=1
   dnorm=dnorm+a(loca)*a(loca+1)*a(loca+2)*ffis(jg)
   go to 200
  440 continue
   do kg=2,ng2
      jg=ng-ig2lo-kg+3
      do id=1,nl
         loca=lz+id+nl*(kg-1)
         if (lpos.gt.0) then
            locs=lpos+itabl*(jg-1)
            sig(locs)=sig(locs)+a(loca)
            ids(led)=led
         endif
         if (kpos.gt.0) then
            if (dnorm.eq.zero) call error('dtfr',&
              'delayed nubar required',&
              'to add delayed chi to total.')
            locs=kpos+itabl*(jg-1)
            sig(locs)=sig(locs)+dnorm*a(loca)
            cnorm=cnorm+dnorm*a(loca)
            ids(ked)=ked
         endif
      enddo
   enddo
   go to 150

   !--accumulate photon production matrices.
  600 continue
   jg=ng-ig+1
   ff=1
   if (mth.eq.102) ff=fcap(jg)
   if (mth.eq.18) ff=ffis(jg)
   if (ig.eq.0) then
      do k=1,ngp
         jgp=ngp-k+1
         spect(jgp)=0
         if (k.ge.ig2lo.and.k.le.(ig2lo+ng2-1)) then
            loca=lz+k-ig2lo+1
            spect(jgp)=a(loca)
         endif
      enddo
   else if (ig2lo.eq.0) then
      loca=lz+ip+nl*((jz-1)+nz)
      do k=1,ngp
         locs=k+ngp*(jg-1)
         sig(locs)=sig(locs)+a(loca)*ff*spect(k)
      enddo
   else
      ipmax=ip
      do k=2,ng2
         jpos=ngp-ig2lo-k+3
         if (jpos.ge.1.and.jpos.le.ngp) then
            loca=lz+ip+nl*((jz-1)+nz*(k-1))
            locs=jpos+ngp*(jg-1)
            sig(locs)=sig(locs)+a(loca)*ff
         endif
      enddo
   endif
   go to 150

   !--accumulation complete
   !--make transport corrections if desired.
   !--normalize fission chi vector if desired.
   !--print out dtf tables if desired
  800 continue
   if (il.eq.1.and.kpos.ne.0.and.cnorm.ne.zero) then
      do k=1,ng
         locs=kpos+itabl*(k-1)
         sig(locs)=sig(locs)/cnorm
      enddo
   endif
   call dtfout(sig,il,ip,ilmax,ipmax)
   if (ip.eq.0) il=il+1
   if (il.gt.nlmax) il=0
   if (il.eq.0) ip=ip+1
   if (ip.gt.nptabl) ip=0
   if (il.gt.0) go to 110
   if (ip.gt.0) go to 110
   go to 105

   !--dtfr is finished.
  900 continue
   ! terminate plotting system
   write(nplot,'('' 99/'')')
   call closz(nplot)
   call closz(nin)
   call closz(npend)
   call closz(nout)
   call closz(nscr)
   call timer(time)
   write(nsyso,'(69x,f8.1,''s''/&
     &1x,7(''**********''),''*******'')') time
   return
   end subroutine dtfr

   subroutine ruin
   !-------------------------------------------------------------------
   ! Read user input.
   !-------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use util   ! provides error,openz,repoz
   ! internals
   integer::i,j,jlast,jedit,nthrm,nedpos
   character(6)::hword
   character(6),parameter::hnusf='nusigf'
   character(6),parameter::hnabs='absorp'
   character(6),parameter::hntotl='total '
   character(6),parameter::hblank='      '
   integer,dimension(50),parameter::kmted=(/2,4,16,6,7,8,9,17,&
     102,107,103,19,20,21,18,104,105,106,111,112,113,114,108,109,&
     32,35,34,22,24,23,25,28,29,30,33,36,37,38,470,471,455,300,&
     301,443,443,444,452,1,0,0/)
   integer,dimension(50),parameter::kjped=(/1,2,3,3,3,3,3,4,5,&
     6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,&
     27,28,29,30,31,32,33,34,35,36,37,38,39,39,40,41,42,43,44,45/)
   integer,dimension(50),parameter::kmultd=(/1,1,1,1,1,1,1,1,1,1,&
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,&
     1,1,1,0,1,1,1,1,1,1/)
   character(5),dimension(50),parameter::hmtid=(/&
     '  els','  ins','  n2n','  n3n','  ngm','  nal','   np',&
     ' fdir','  nnf',' n2nf',' ftot','   nd','   nt',' nhe3',&
     '  n2p','  npa',' nt2a',' nd2a','  n2a','  n3a','  nnd',&
     'nnd2a','nnhe3','  nna',' n2na',' nn3a',' n3na','  nnp',&
     ' nn2a','n2n2a','  nnt','nnt2a','  n4n',' n3nf','  chi',&
     ' chid','  nud','  phi','theat','kerma','tdame',' nusf',&
     ' totl','     ','     ','     ','     ','     ','     ',&
     '     '/)

   !--read initial parameters
   mti=0
   mtc=0
   nlc=0
   nout=0
   npend=0
   nplot=0
   read(nsysi,*) nin,nout,npend,nplot
   call openz(nin,0)
   call openz(nout,1)
   call openz(npend,0)
   call openz(nplot,1)
   ifilm=0
   iedit=0
   read(nsysi,*) iprint,ifilm,iedit
   if (nplot.eq.0) ifilm=0
   write(nsyso,'(/&
     &'' input gendf unit ..................... '',i10/&
     &'' output unit .......................... '',i10/&
     &'' input pendf unit ..................... '',i10/&
     &'' output plot data ..................... '',i10)')&
     nin,nout,npend,nplot
   write(nsyso,'(&
     &'' print option (0 min, 1 max) .......... '',i10/&
     &'' film option (none, 1/plot, 4/plot ) .. '',i10/&
     &'' edit option (0 table, 1 sep.)......... '',i10)')&
     iprint, ifilm, iedit

   !--read parameter values for iedit=0
   if (iedit.eq.0) then
      ned=0
      nthrm=0
      read(nsysi,*) nlmax,ng,iptotl,ipingp,itabl,ned,nthrm
      if (ipingp.le.iptotl) call error('dtfr','iping.le.iptotl.',' ')
      jgthrm=ng+1-nthrm
      if (nthrm.ne.0) then
         mtc=0
         nlc=0
         read(nsysi,*) mti,mtc,nlc
      endif

   !--set standard values for iedit=1 (claw format)
   else
      nlmax=5
      ng=30
      read(nsysi,*) nlmax,ng
      iptotl=43+3
      ned=48
      jgthrm=ng+1
      ipingp=iptotl+1
      itabl=iptotl+ng
      nedpos=iptotl-3
      nptabl=0
      do i=1,nedpos
      hednam(i)=hmtid(i)
      enddo
      do i=1,ned
         jped(i)=kjped(i)
         mted(i)=kmted(i)
         multed(i)=kmultd(i)
      enddo
   endif

   !--write out neutron table specifications
   write(nsyso,'(&
     &'' number of neutron tables ............. '',i10/&
     &'' number of neutron groups ............. '',i10)')&
     nlmax, ng
   if (nlmax.gt.0) then
      write(nsyso,'(&
        &'' position of total .................... '',i10/&
        &'' position of in-group ................. '',i10/&
        &'' table length ......................... '',i10)')&
        iptotl, ipingp, itabl
      if (jgthrm.le.ng) write(nsyso,'(&
        &'' thermal break group .................. '',i10)')&
        jgthrm
      if (jgthrm.le.ng) write(nsyso,'(&
        &'' thermal incoherent mt ................ '',i10/&
        &'' thermal coherent mt .................. '',i10/&
        &'' no. coherent legendre orders ......... '',i10)')&
        mti,mtc,nlc

      !--read edit descriptions for iedit=0
      if (iedit.eq.0) then
         if (ned+3.gt.nedmax)&
           call error('dtfr','not enough storage for edits.',' ')
         if (ned.ne.0) then
            nedpos=iptotl-3
            if (nedpos.ne.0) then
              read(nsysi,*) (hednam(i),i=1,nedpos)
            endif
            read(nsysi,*) (jped(j),mted(j),multed(j),j=1,ned)
         endif
      endif

      !--write out edit specifications
      if (ned.gt.0) then
         write(nsyso,'(&
           &'' edit cross sections .................. '',&
           &'' name   position reaction multiplicity'')')
         jlast=0
         jedit=0
         kpos=0
         lpos=0
         do i=1,ned
            hword=hblank
            if (nedpos.ne.0.and.jped(i).ne.jlast) then
               jedit=jedit+1
               jlast=jped(i)
               hword=hednam(jedit)
            endif
            if (mted(i).eq.470) kpos=jped(i)
            if (mted(i).eq.470) ked=i
            if (mted(i).eq.471) lpos=jped(i)
            if (mted(i).eq.471) led=i
            write(nsyso,'(40x,a6,3x,i2,7x,i3,8x,i2)')&
              hword,jped(i),mted(i),multed(i)
         enddo
      endif
   endif

   !--add three standard edits
   hednam(iptotl-1)=hnusf
   jped(ned+2)=iptotl-1
   mted(ned+2)=1000
   hednam(iptotl-2)=hnabs
   jped(ned+1)=iptotl-2
   mted(ned+1)=1000
   hednam(iptotl)=hntotl
   jped(ned+3)=iptotl
   mted(ned+3)=1

   !--read and display photon table specifications
   nptabl=0
   ngp=0
   read(nsysi,*) nptabl,ngp
   write(nsyso,'(&
     &'' number of photon tables .............. '',i10)') nptabl
   if (nptabl.gt.0) write(nsyso,'(&
     &'' number of photon groups .............. '',i10)') ngp
   write(nsyso,'(/)')
   call repoz(nout)
   return
   end subroutine ruin

   subroutine dtfout(sig,il,ip,ilmax,ipmax)
   !-------------------------------------------------------------------
   ! Write and/or plot group constants in DTF format.
   !-------------------------------------------------------------------
   use util ! provides dater
   use mainio ! provides nsyso
   ! externals
   integer::il,ip,ilmax,ipmax
   real(kr)::sig(*)
   ! internals
   integer::j,nj,ifis,i,jpos,iseq,k,ig,locs,i1,ii,l
   integer::med,nws,ltabn,igp,nw
   real(kr)::dat(6)
   character(6)::hmti
   character(8)::hdat
   character(8),parameter::id='njoy    '

   !--special format 0 -- output dtf tables with internal edits.
   if (iedit.ne.0) go to 130
   if (il.ne.0) then
      nw=itabl*ng
      if (nout.gt.0) write(nout,'(/'' il= '',i1,'' table'',i3,&
        &'' gp'',i3,'' pos, mat='',i5,'' iz='',i2,&
        &'' temp='',1p,e12.5)') il,ng,itabl,matd,jz,dtemp
      if (iprint.eq.1) write(nsyso,'(/'' il= '',i1,'' table'',i3,&
        &'' gp'',i3,'' pos, mat='',i5,'' iz='',i2,&
        &'' temp='',1p,e12.5)') il,ng,itabl,matd,jz,dtemp
   else
      nw=ngp*ng
      if (nout.gt.0) write(nout,'(/'' ip= '',i1,'' table'',i3,&
        &'' gp'',i3,'' pos, mat='',i5,'' iz='',i2,&
        &'' temp='',1p,e12.5)') ip,ng,ngp,matd,jz,dtemp
      if (iprint.eq.1) write(nsyso,'(/'' ip= '',i1,'' table'',i3,&
        &'' gp'',i3,'' pos, mat='',i5,'' iz='',i2,&
        &'' temp='',1p,e12.5)') ip,ng,ngp,matd,jz,dtemp
   endif
   if (nout.gt.0) write(nout,'(1p,6e12.4)') (sig(j),j=1,nw)
   if (iprint.eq.1) write(nsyso,'(1x,1p,6e12.4)') (sig(j),j=1,nw)
   go to 300

   !--special format 1 -- td6 format.  write separate edits
  130 continue
   if (il.eq.0) go to 250
   if (il.gt.1) go to 200
   call dater(hdat)
   if (nout.gt.0) write(nout,'(a6,''     edit xsec ('',&
     &i3,''x'',i3,'') proc by '',a6,'' on '',2a8)')&
     hisnam,ng,ned,id,hdat
   if (iprint.eq.1) write(nsyso,'(1x,a6,''     edit xsec ('',&
     &i3,''x'',i3,'') proc by '',a6,'' on '',2a8)')&
     hisnam,ng,ned,id,hdat
   ifis=ids(ned+2)
   nj=iptotl-3
   do 140 j=1,nj
   if (j.le.7) go to 150
   if (ifis.gt.0.and.j.le.11) go to 150
   if (j.ge.nj-1) go to 150
   do 145 i=1,ned
   if (jped(i).eq.j.and.ids(i).gt.0) go to 150
  145 continue
   go to 140
  150 continue
   jpos=j
   if (j.eq.nj-1.and.ifis.eq.0) go to 140
   if (j.eq.nj-1) jpos=iptotl-1
   if (j.eq.nj) jpos=iptotl
   hmti=hednam(j)
   iseq=0
   k=0
   do ig=1,ng
      locs=jpos+itabl*(ig-1)
      k=k+1
      dat(k)=sig(locs)
      if (k.ge.6.or.ig.ge.ng) then
         iseq=iseq+1
         if (k.ne.6) then
            i1=k+1
            do ii=i1,6
               dat(ii)=0
            enddo
         endif
         if (nout.gt.0) write(nout,'(1p,6e12.5,a2,a5,i1)')&
           (dat(i),i=1,6),hisnam,hmti,iseq
         if (iprint.eq.1) write(nsyso,'(1x,1p,6e12.5,a2,a5,i1)')&
           (dat(i),i=1,6),hisnam,hmti,iseq
         k=0
      endif
   enddo
  140 continue

   !--write reduced neutron table with edits removed.
  200 continue
   l=il-1
   if (il.le.ilmax) then
      med=iptotl-3
      nws=ng*itabl
      iseq=0
      ltabn=itabl-iptotl+3
      if (nout.gt.0) write(nout,'(a6,'' l='',i1,&
        &'' n-n table ('',i3,''x'',i3,'')'')')&
        hisnam,l,ltabn,ng
      if (iprint.eq.1) write(nsyso,'(1x,a6,'' l='',i1,&
        &'' n-n table ('',i3,''x'',i3,'')'')')&
        hisnam,l,ltabn,ng
      k=0
      do ig=1,ng
         do jpos=1,itabl
            if (jpos.gt.med) then
               k=k+1
               locs=jpos+itabl*(ig-1)
               dat(k)=sig(locs)
               if (k.ge.6.or.locs.ge.nws) then
                  iseq=iseq+1
                  if (k.ne.6) then
                     i1=k+1
                     do i=i1,6
                        dat(i)=0
                     enddo
                  endif
                  if (nout.gt.0) write(nout,&
                    '(1p,6e12.5,a2,i2,i4)')&
                    (dat(i),i=1,6),hisnam,l,iseq
                  if (iprint.eq.1) write(nsyso,&
                    '(1p,6e12.5,a2,i2,i4)')&
                    (dat(i),i=1,6),hisnam,l,iseq
                  k=0
               endif
            endif
         enddo
      enddo
   endif
   go to 300

   !--write photon table
  250 continue
   l=ip-1
   if (ip.le.ipmax) then
      iseq=0
      if (nout.gt.0) write(nout,'(1x,a6,'' l='',i1,&
        &''h n-p table ('',i3,''x'',i3,'')'')')&
        hisnam,l,ngp,ng
      if (iprint.eq.1) write(nsyso,'(a6,'' l='',i1,&
        &'' n-p table ('',i3,''x'',i3,'')'')')&
        hisnam,l,ngp,ng
      k=0
      nws=ng*ngp
      do ig=1,ng
         do igp=1,ngp
            k=k+1
            locs=igp+ngp*(ig-1)
            dat(k)=sig(locs)
            if (k.ge.6.or.locs.ge.nws) then
               iseq=iseq+1
               if (k.ne.6) then
                  i1=k+1
                  do i=i1,6
                     dat(i)=0
                  enddo
               endif
               if (nout.gt.0) write(nout,&
                 '(1p,6e12.5,a2,i2,i4)')&
                 (dat(i),i=1,6),hisnam,l,iseq
               if (iprint.eq.1) write(nsyso,&
                 '(1p,6e12.5,a2,i2,i4)')&
                 (dat(i),i=1,6),hisnam,l,iseq
               k=0
            endif
         enddo
      enddo
   endif

   !--plot dtf edit quantities and scattering tables
  300 continue
   if (ifilm.gt.0.and.il.eq.1) call ploted(sig)
   if (ifilm.gt.0.and.il.eq.1.and.il.le.ilmax) call plotnn(sig,il)
   if (ifilm.gt.0.and.ip.eq.1.and.ip.le.ipmax) call plotnp(sig,ip)
   return
   end subroutine dtfout

   subroutine ploted(sig)
   !-------------------------------------------------------------------
   ! Generate plots of all edit and cross section positions.
   ! Overlay corresponding pointwise cross sections from pendf.
   !-------------------------------------------------------------------
   use endf ! provides iverf
   ! externals
   real(kr)::sig(*)
   ! internals
   integer::i,nk,k,ksave,jpos,mt,locs
   integer::nh,ih,j,iy1,lin,nsym,ndash,ipt,npts
   real(kr)::deltay,xll,yll,xur,yur,xtag,ytag
   character(6)::hedn
   character(16)::labelz
   character(8)::l1
   character(1)::qu=''''
   character(1),dimension(5),parameter::nchar=(/' ','x','+','*','0'/)
   real(kr),parameter::one=1
   real(kr),parameter::ten=10
   real(kr),parameter::ww=5.e0_kr
   real(kr),parameter::hh=3.75e0_kr
   real(kr),parameter::zero=0

   !--plot cross sections (non-zero reactions only)
   nplt=iptotl
   deltay=16
   do 100 i=1,nplt
   nk=ned+3
   do 110 k=1,nk
   ksave=k
   if (jped(k).eq.i.and.ids(k).gt.0) go to 112
  110 continue
   go to 100
  112 continue
   jpos=i
   mt=mted(ksave)
   if (mt.eq.300) go to 100
   if (mt.eq.1.and.i.ne.iptotl) go to 100
   hedn=hednam(jpos)

   !--draw histogram.
   write(labelz,'(a6,4x,a6)') hisnam,hedn
   do k=1,ngn
      locs=jpos+itabl*(ngn-k)
      y(k)=sig(locs)
   enddo
   call histod(x,y,z,nh,ngn,egn)
   if (ifilm.eq.1) then
      write(nplot,'('' 1/'')')
   else
      xll=ww*ixo
      yll=hh-hh*iyo
      xur=ww
      yur=hh
      if (ixo.eq.0.and.iyo.eq.0) then
         write(nplot,'('' 1 0 1. 1.'',5f6.2,''/'')')&
           xll,yll,xur,yur,zero
      else
         write(nplot,'('' -1 0 1. 1.'',5f6.2,''/'')')&
           xll,yll,xur,yur,zero
      endif
   endif
   if (mt.ge.500.and.mt.le.522) then
      if (ifilm.eq.1) xtag=ten**xmax/10
      if (ifilm.eq.2) xtag=3*ten**xmax/100
   else
      xtag=xmin+(xmax-xmin)/50
      xtag=ten**xtag
   endif
   ytag=ymax-(ymax-ymin)/50
   ytag=ten**ytag
   write(nplot,'(1x,a,''<'',a20,a,''/'')') qu,labelz,qu
   write(nplot,'('' /'')')
   write(nplot,'('' 4 0 3 1'',2e10.2,''/'')') xtag,ytag
   write(nplot,'(1p,3e10.2,''/'')') ten**xmin,ten**xmax,one
   write(nplot,'(1x,a,''<e>nergy (e<v>)'',a,''/'')') qu,qu
   write(nplot,'(1p,3e10.2,''/'')') ten**ymin,ten**ymax,one
   write(nplot,'(1x,a,''<c>ross <s>ection (barns)'',a,''/'')') qu,qu
   write(nplot,'(''/'')')
   write(nplot,'('' 0 0 0 0 1/'')')
   write(nplot,'(1x,a,''<m.g.>'',a,''/'')') qu,qu
   write(nplot,'('' 0'')')
   do ih=1,nh
      write(nplot,'(1p,2e13.4,''/'')') x(ih),z(ih)
   enddo
   write(nplot,'('' /'')')
   if (ifilm.ne.1) then
      ixo=ixo+1
      if (ixo.gt.1) ixo=0
      if (ixo.eq.0) iyo=iyo+1
      if (iyo.gt.1) iyo=0
   endif
   if (iphph.eq.1) then
      if (iverf.ge.6) then
         if (mt.lt.501.or.mt.gt.522) go to 100
      else
         if (mt.lt.501.or.mt.gt.602) go to 100
      endif
   else
      if (iverf.ge.6) then
         if (mt.gt.150.and.mt.lt.301) go to 100
         if (mt.gt.450.and.mt.lt.600) go to 100
         if (mt.gt.849) go to 100
      else
         if (mt.gt.150.and.mt.lt.301) go to 100
         if (mt.gt.450.and.mt.lt.700) go to 100
         if (mt.gt.799) go to 100
      endif
   endif

   !--draw smooth plot(s)
   j=0
   k=ksave
   iy1=10
  140 continue
   j=j+1
   write(l1,'(''mt'',i3)') mt
   nsym=ichar(nchar(j))
   call dpend(a,dtemp,matd,mt,x,y,npts)
   if (j.lt.5) then
      lin=0
      nsym=0
      ndash=j
   else if (j.ge.5.and.j.lt.9) then
      lin=1+npts/5
      nsym=4
      ndash=j-4
   else
      lin=1+npts/5
      nsym=3
      ndash=j-8
   endif
   write(nplot,'(i3,''/ pendf plot for mt='',i4)') j+1,mt
   write(nplot,'(''/'')')
   write(nplot,'(i5,i3,i3,'' 0 1/'')') lin,nsym,ndash
   write(nplot,'(1x,a,''<'',a8,a,''/'')') qu,l1,qu
   write(nplot,'('' 0'')')
   do ipt=1,npts
      write(nplot,'(1p,2e13.4,''/'')') x(ipt),y(ipt)
   enddo
   write(nplot,'('' /'')')
   iy1=iy1+nint(deltay)
  150 continue
   k=k+1
   if (k.gt.ned) go to 100
   if (ids(k).eq.0.or.jped(k).ne.i) go to 150
   mt=mted(k)
   go to 140
  100 continue
   return
   end subroutine ploted

   subroutine histod(x,y,z,n,ng,eg)
   !-------------------------------------------------------------------
   ! Prepare histogram data for plotting and choose scales for plots.
   !-------------------------------------------------------------------
   ! externals
   integer::n,ng
   real(kr)::x(*),y(*),z(*),eg(*)
   ! internals
   integer::i,iph,ipl,ipx
   real(kr)::yh,yl,ayh,ayl,ynow,axh,xnow
   real(kr),parameter::zero=0
   real(kr),parameter::two=2
   real(kr),parameter::five=5
   real(kr),parameter::ten=10
   real(kr),parameter::ptwo=.2e0_kr
   real(kr),parameter::pfive=.5e0_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::big=1.e+10_kr
   real(kr),parameter::fact=1.7e0_kr

   !--construct histogram removing leading and trailing zeros
   yh=small
   yl=big
   n=1
   z(n)=0
   do i=1,ng
      if (y(i).ne.zero.or.n.ne.1) then
         n=n+2
         z(n-1)=y(i)
         z(n)=y(i)
         x(n-2)=eg(i)
         x(n-1)=eg(i)
         if (y(i).gt.yh) yh=y(i)
         if (y(i).gt.zero.and.y(i).lt.yl) yl=y(i)
      endif
   enddo
   x(n)=eg(ng+1)
   n=n+1
   x(n)=eg(ng+1)
   z(n)=0

   !--set scales for plotting
   z(n)=0
   yh=fact*yh
   ayh=log10(yh)
   iph=int(ayh)
   if (iph.ge.0) iph=iph+1
   ayl=log10(yl)
   ipl=int(ayl-1)
   if ((iph-ipl).gt.6) ipl=iph-6
   if (iph.eq.ipl) iph=iph+1
   ymin=ipl
   ymax=iph
   if (iph-ipl.le.3) then
      ynow=ymin
      if (yl.gt.two*ten**ynow) ymin=ynow+log10(two)
      if (yl.gt.five*ten**ynow) ymin=ynow+log10(five)
      ynow=ymax
      if (yh.lt.pfive*ten**ynow) ymax=ynow+log10(pfive)
      if (yh.lt.ptwo*ten**ynow) ymax=ynow+log10(ptwo)
   endif
   bot=ten**ymin
   top=ten**ymax
   axh=log10(x(n))
   iph=int(axh)
   if (axh-iph.ne.zero) iph=iph+1
   xmax=iph
   ipl=int(log10(x(1))+small)
   ipx=-4
   if (ipl.gt.ipx) ipx=-2
   if (ipl.gt.ipx) ipx=4
   if (ipl.gt.ipx) ipx=5
   if (ipl.gt.ipx) ipx=6
   xmin=ipx
   if (iph-ipx.le.2) then
      xnow=xmax
      if (x(n).lt.pfive*ten**xnow) xmax=xnow+log10(pfive)
      if (x(n).lt.ptwo*ten**xnow) xmax=xnow+log10(ptwo)
      xnow=xmin
      if (x(1).gt.two*ten**xnow) xmin=xnow+log10(two)
      if (x(1).gt.five*ten**xnow) xmin=xnow+log10(five)
   endif
   axl=xmin
   axd=xmax-xmin
   do i=1,n
      if (z(i).le.bot) z(i)=bot
   enddo
   return
   end subroutine histod

   subroutine dpend(a,tempr,mat,mt,x,y,npts)
   !-------------------------------------------------------------------
   ! Retrieve point cross sections from pendf file and thin for
   ! plotting if necessary.  Add points for plotting lin-lin
   ! functions on a log-log graph if necessary.
   !-------------------------------------------------------------------
   use util ! provides error
   use endf ! provides findf,contio,gety1
   ! externals
   integer::mat,mt,npts
   real(kr)::tempr,a(*),x(*),y(*)
   ! internals
   integer::mft,nb,nw,idis,ixlast,ix,ns
   real(kr)::e,enext,fact,step,elow,elast,s,s1,s2
   real(kr),parameter::ehigh=5.e7_kr
   real(kr),parameter::one=1
   real(kr),parameter::ten=10
   integer,parameter::nstep=50
   real(kr),parameter::zero=0

   !--find desired reaction.
   npts=0
   mft=3
   if (iphph.eq.1) mft=23
   call findf(mat,mft,mt,npend)
   call contio(npend,0,0,a,nb,nw)
   e=0
   call gety1(e,enext,idis,s,npend,a)
   fact=2000/axd
   step=ten**(one/nstep)
   elow=ten**axl

   !--read this reaction, thinning or thickening as necessary.
   ixlast=-100
   elast=enext
  210 continue
   e=enext
   if (e.gt.step*elast) e=step*elast
   elast=e
   call gety1(e,enext,idis,s,npend,a)
   if (enext.lt.elow) go to 210
   if (s.lt.zero) s=-s
   if (s.gt.top) s=top
   if (s.lt.bot) s=bot
   if (e.lt.elow) e=elow
   ix=int((log10(e)-axl)*fact)
   if (ix.eq.ixlast) go to 220
   ixlast=ix
   ns=0
   go to 230
  220 continue
   if (s.eq.y(npts)) go to 260
   go to (230,240,240,250,270),ns
  240 continue
   if ((s-y(npts))*(y(npts)-y(npts-1)).gt.zero) go to 260
  230 continue
   npts=npts+1
   if (npts.gt.ndim) call error('dpend','npts exceeds ndim.',' ')
   ns=ns+1
  260 continue
   x(npts)=e
   y(npts)=s
   if (enext.le.ehigh) go to 210
   return
  250 continue
   ns=5
   s2=y(npts-2)
   s1=y(npts-1)
   if (s2.gt.s1) go to 270
   y(npts-2)=s1
   y(npts-1)=s2
  270 continue
   if (s.le.y(npts-2)) go to 280
   y(npts-2)=s
   go to 260
  280 continue
   if (s.ge.y(npts-1)) go to 260
   y(npts-1)=s
   go to 260
   end subroutine dpend

   subroutine plotnn(sig,il)
   !-------------------------------------------------------------------
   ! Generate an isometric plot of a neutron table.
   ! This routine destroys the table.
   !-------------------------------------------------------------------
   ! externals
   real(kr)::sig(*)
   integer::il
   ! internals
   integer::ngn,ltabn,j,k,ig,ig2,i,lim,ip,l,ig1,ig2l
   real(kr)::top,bot,dl,sigj,xll,yll,xur,yur,zlast,ss
   character(26)::ititle*26
   character(1),parameter::qu=''''
   real(kr),parameter::emin=1.e3_kr
   real(kr),parameter::emax=1.e8_kr
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::ten=10
   real(kr),parameter::ww=5.0e0_kr
   real(kr),parameter::hh=3.75e0_kr
   real(kr),parameter::zero=0

   !--remove edits and cross sections.
   ngn=ng
   ltabn=itabl-iptotl
   j=0
   k=0
   top=-big
   do ig=1,ngn
      j=itabl*(ig-1)+iptotl
      ig2=ngn-ig+1
      dl=log(egn(ig2+1))-log(egn(ig2))
      do ip=1,ltabn
         k=k+1
         j=j+1
         sigj=sig(j)
         if (egn(ig2).lt.emin) sigj=0
         if (egn(ig2+1).gt.emax) sigj=0
         sig(k)=sigj/dl
         if (sig(k).gt.top) top=sig(k)
      enddo
   enddo
   i=int(log10(top))
   if (i.ge.0) i=i+1
   bot=ten**(i-5)
   top=ten**i
   lim=ngn*ltabn
   do k=1,lim
      if (sig(k).lt.bot) sig(k)=bot
   enddo

   !--plot the modified table
   l=il-1
   if (iphph.eq.0) write(ititle,&
     '(a6,'' l='',i1,'' neut-neut table'')') hisnam,l
   if (iphph.eq.1) write(ititle,&
     '(a6,'' l='',i1,'' phot-phot table'')') hisnam,l
   xmin=egn(1)
   if (xmin.lt.emin) xmin=emin
   i=int(log10(xmin))
   if (i.lt.0) i=i-1
   xmin=ten**i
   xmax=egn(ngn+1)
   if (xmax.gt.emax) xmax=emax
   i=int(log10(xmax))
   if (i.ge.0) i=i+1
   xmax=ten**i
   if (ifilm.eq.1) then
      write(nplot,'('' 1/ 3d data'')')
   else
      xll=ww*ixo
      yll=hh-hh*iyo
      xur=ww
      yur=hh
      if (ixo.eq.0.and.iyo.eq.0) then
         write(nplot,'('' 1 0 1. 1.'',5f6.2,''/ 3d data'')')&
           xll,yll,xur,yur,zero
      else
         write(nplot,'('' -1 0 1. 1.'',5f6.2,''/ 3d data'')')&
           xll,yll,xur,yur,zero
      endif
   endif
   write(nplot,'(1x,a,''<'',a26,a,''/'')') qu,ititle,qu
   write(nplot,'('' /'')')
   write(nplot,'('' -4 2/'')')
   write(nplot,'(1p,2e12.4,'' 1./'')') xmin,xmax
   write(nplot,'(1x,a,''<s>ec. <e>nergy'',a,''/'')') qu,qu
   write(nplot,'(1p,2e12.4,'' 1./'')') xmin,xmax
   write(nplot,'(1x,a,''<e>nergy (e<v>)'',a,''/'')') qu,qu
   write(nplot,'(1p,2e12.4,'' 1./'')') bot,top
   write(nplot,'(1x,a,''<x>sec/leth'',a,''/'')') qu,qu
   write(nplot,'(''/'')')
   if (ifilm.eq.1) then
      write(nplot,'('' 15. -15. 15. -2.5 6.5 2.5/'')')
   else
      write(nplot,'('' 12. -12. 12. -2. 5.25 2./'')')
   endif
   write(nplot,'('' 1/ 3d data'')')
   do ig1=1,ngn
      if (egn(ig1).ge.xmin.and.egn(ig1+1).le.xmax) then
         i=ngn-ig1+1
         write(nplot,'(1x,1p,e12.4,''/'')') egn(ig1)
         zlast=bot
         do ig2=1,ngn
            if (egn(ig2).ge.xmin.and.egn(ig2+1).le.xmax) then
               j=ngn-ig2+1
               k=j-i+ipingp-iptotl
               if (k.ge.1) then
                  ss=sig(k+ltabn*(j-1))
                  if (k.gt.ltabn) ss=bot
                  if (ss.lt.bot) ss=bot
                  write(nplot,'(1x,1p,2e12.4,''/'')') egn(ig2),zlast
                  write(nplot,'(1x,1p,2e12.4,''/'')') egn(ig2),ss
                  zlast=ss
                  ig2l=ig2
               endif
            endif
         enddo
         write(nplot,'(1x,1p,2e12.4,''/'')') egn(ig2l+1),zlast
         write(nplot,'(1x,1p,2e12.4,''/'')') egn(ig2l+1),bot
         write(nplot,'('' /'')')
      endif
   enddo
   write(nplot,'('' /'')')
   if (ifilm.eq.1) return
   ixo=ixo+1
   if (ixo.gt.1) ixo=0
   if (ixo.eq.0) iyo=iyo+1
   if (iyo.gt.1) iyo=0
   return
   end subroutine plotnn

   subroutine plotnp(sig,ip)
   !-------------------------------------------------------------------
   ! Generate an isometric plot of a photon table.
   ! This routine destroys the table.
   !-------------------------------------------------------------------
   ! externals
   real(kr)::sig(*)
   integer::ip
   ! internals
   integer::ign,i,igp,j,igpl
   real(kr)::xll,yll,xur,yur,top,bot,ss,xmin,xmax,zlast
   character(26)::ititle
   character(1),parameter::qu=''''
   real(kr),parameter::emin=1.e3_kr
   real(kr),parameter::emax=1.e8_kr
   real(kr),parameter::ten=10
   real(kr),parameter::ww=5
   real(kr),parameter::hh=3.75e0_kr
   real(kr),parameter::zero=0

   write(ititle,'(a6,'' l='',i1,'' neut-phot table'')') hisnam,ip-1
   if (ifilm.eq.1) then
      write(nplot,'('' 1/ 3d data'')')
   else
      xll=ww*ixo
      yll=hh-hh*iyo
      xur=ww
      yur=hh
      if (ixo.eq.0.and.iyo.eq.0) then
         write(nplot,'('' 1 0 1. 1.'',5f6.2,''/ 3d data'')')&
           xll,yll,xur,yur,zero
      else
         write(nplot,'('' -1 0 1. 1.'',5f6.2,''/ 3d data'')')&
           xll,yll,xur,yur,zero
      endif
   endif
   top=0
   do ign=1,ngn
      i=ngn-ign+1
      do igp=1,ngp
         j=ngp-igp+1
         ss=sig(j+ngp*(i-1))/(egg(igp+1)-egg(igp))
         if (egn(ign).ge.emin.and.egn(ign+1).le.emax) then
            if (egg(igp).ge.emin.and.egg(igp+1).le.emax) then
               if (ss.gt.top) top=ss
            endif
         endif
      enddo
   enddo
   i=int(log10(top))
   if (i.ge.0) i=i+1
   bot=ten**(i-5)
   top=ten**i
   xmin=egn(1)
   if (xmin.lt.emin) xmin=emin
   i=int(log10(xmin))
   if (i.lt.0) i=i-1
   xmin=ten**i
   xmax=egn(ngn+1)
   if (xmax.gt.emax) xmax=emax
   i=int(log10(xmax))
   if (i.ge.0) i=i+1
   xmax=ten**i
   write(nplot,'(1x,a,''<'',a26,a,''/'')') qu,ititle,qu
   write(nplot,'('' /'')')
   write(nplot,'('' -4 2/'')')
   write(nplot,'(1p,2e12.4,'' 1./'')') xmin,xmax
   write(nplot,'(1x,a,''<g>amma <e>nergy'',a,''/'')') qu,qu
   write(nplot,'(1p,2e12.4,'' 1./'')') xmin,xmax
   write(nplot,'(1x,a,''<e>nergy (e<v>)'',a,''/'')') qu,qu
   write(nplot,'(1p,2e12.4,'' 1./'')') bot,top
   write(nplot,'(1x,a,''<x>sec/e<v>'',a,''/'')') qu,qu
   write(nplot,'(''/'')')
   if (ifilm.eq.1) then
      write(nplot,'('' 15. -15. 15. -2.5 6.5 2.5/'')')
   else
      write(nplot,'('' 12. -12. 12. -2. 5.25 2./'')')
   endif
   write(nplot,'('' 1/ 3d data'')')
   do ign=1,ngn
      if (egn(ign).ge.xmin.and.egn(ign+1).le.xmax) then
         i=ngn-ign+1
         write(nplot,'(1x,1p,e12.4,''/'')') egn(ign)
         zlast=bot
         do igp=1,ngp
            j=ngp-igp+1
            ss=sig(j+ngp*(i-1))/(egg(igp+1)-egg(igp))
            if (ss.lt.bot) ss=bot
            write(nplot,'(1x,1p,2e12.4,''/'')') egg(igp),zlast
            write(nplot,'(1x,1p,2e12.4,''/'')') egg(igp),ss
            zlast=ss
            igpl=igp
         enddo
         write(nplot,'(1x,1p,2e12.4,''/'')') egg(igpl+1),zlast
         write(nplot,'(1x,1p,2e12.4,''/'')') egg(igpl+1),ss
         write(nplot,'('' /'')')
      endif
   enddo
   write(nplot,'('' /'')')
   if (ifilm.eq.1) return
   ixo=ixo+1
   if (ixo.gt.1) ixo=0
   if (ixo.eq.0) iyo=iyo+1
   if (iyo.gt.1) iyo=0
   return
   end subroutine plotnp

end module dtfm


module purm
   ! provides subroutine purr for NJOY2016
   use locale
   implicit none
   private
   public purr

   ! units
   integer::nendf,nin,nout

   ! control parameters
   integer::init
   integer::nunr,lssf,iinel,iabso
   integer::iprint,nermax,nladr,nmode

   ! arrays for the w function
   real(kr)::tr(41,27),ti(41,27),trs(41,27),tis(41,27)

   ! arrays for unresolved resonance parameter sections
   integer::nsect
   integer::isot(20),modet(20),ibaset(20)
   real(kr)::abnt(20),elt(20),eht(20)

   ! resonance data
   integer::nseq0,nro,naps,intunr
   integer,parameter::mxns0=100
   real(kr)::e,cth(mxns0),csz(mxns0),cc2p(mxns0),cs2p(mxns0),&
     cgn(mxns0),cgg(mxns0),cgf(mxns0),cgx(mxns0),cgt(mxns0),&
     dbar(mxns0),spot,dbarin,sigi(4)
   integer::ndfn(mxns0),ndff(mxns0),ndfx(mxns0)

   ! array for unresolved energy grid
   integer,parameter::meunr=150
   real(kr)::eunr(meunr)

   ! array for background cross sections
   real(kr),dimension(:),allocatable::sb

   ! probability table globals
   real(kr),dimension(:,:,:),allocatable::bval
   real(kr),dimension(:),allocatable::tmin,tmax,tsum

   ! storage array for unresolved resonance parameters
   integer,parameter::jx=10000
   real(kr),dimension(:),allocatable::arry

   ! optional plots
   integer ipl
   real(kr)::epl(200)
   real(kr),dimension(:,:,:),allocatable::sigpl

   ! random number generator
   integer::kk
   real(kr)::rr

   ! message flags
   integer::mflg1,mflg2

contains

   subroutine purr
   !-------------------------------------------------------------------
   !
   !  Probabalistic unresolved calculation of
   !  Bondarenko moments and probability tables
   !
   !-------------------------------------------------------------------
   !
   !  Purr constructs a series of resonance ladders that obey the
   !  distributions given in MT151 of the ENDF tape.  Each ladder
   !  is sampled randomly to produce contributions to a probability
   !  table and a set of bondarenko moments.  When the table is
   !  complete, Bondarenko moments are computed from the table to
   !  provide a convergence check.  All temperatures are computed
   !  simultaneously to preserve temperature correlations.  The
   !  Bondarenko tables are written on the pendf tape using MT152,
   !  and the probability tables are written using MT153.
   !  A conditional probability for heating is added to the table.
   !  If partial heating cross sections for elastic (302), fission
   !  (318), and capture (402) are available from heatr, full
   !  fluctuations will be provided for the total heating.
   !  Otherwise, the same value will be provided for each bin.
   !
   !---input data cards---------------------------------------------
   !
   ! card 1
   !   nendf   unit for endf tape
   !   nin     unit for input pendf tape
   !   nout    unit for output pendf tape
   ! card 2
   !   matd    material to be processed
   !           matd=0 terminates purr
   !   ntemp   no. of temperatures (default=1)
   !   nsigz   no of sigma zeros (default=1)
   !   nbin    no. of probability bins (15 or more)
   !   nladr   no. of resonance ladders
   !   iprint  print option (0=min, 1=max, def=1)
   !   nunx    no. of energy points desired (def=0=all)
   ! card 3
   !   temp    temperatures in kelvin (including zero)
   ! card 4
   !   sigz    sigma zero values (including infinity)
   !
   !-------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use util ! provides timer,error,openz,repoz,closz,sigfig
   use endf ! provides endf routines and variables
   ! internals
   integer::lrp,nstep,ie,i,j,k,new152,new153,it,nx,newmat
   integer::matd,mfd,mtd,ix,is,n1,nsamp,nscr,maxscr,iscr,nb,nw
   integer::nunx,ntemp,nsigz,nbin,nwds,n,lstart,l
   integer::nwd,nc153,ncds,ihave,nc,nd,maxtst
   real(kr)::time,za,awr,ez,sigx,temz,h
   real(kr)::bkgz(4)
   character(60)::strng1,strng2
   real(kr),dimension(:),allocatable::a
   real(kr),dimension(:,:,:),allocatable::tabl
   real(kr),dimension(:,:),allocatable::tval
   real(kr),dimension(:),allocatable::er,gnr,gfr,ggr,gxr,gt
   real(kr),dimension(:),allocatable::es,xs
   real(kr),dimension(:,:),allocatable::fis,cap,els
   real(kr),dimension(:,:,:),allocatable::heat
   real(kr),dimension(:),allocatable::c,d
   real(kr),dimension(:),allocatable::temp,sigz
   real(kr),dimension(:,:,:),allocatable::sigu
   real(kr),parameter::tref=300.e0_kr
   real(kr),parameter::sigmin=1.e-5_kr
   real(kr),parameter::zero=0

   !--initialize
   ipl=0
   nmode=0
   nermax=1000
   nsamp=10000
   nsh=0
   call timer(time)
   write(nsyso,'(/'' purr...'',&
     &''probabalistic unresolved calculation'',&
     &25x,f8.1,''s'')') time
   write(nsyse,'(/'' purr...'',61x,f8.1,''s'')') time
   read(nsysi,*) nendf,nin,nout
   if ((nin.lt.0.and.nout.gt.0).or.(nin.gt.0.and.nout.lt.0))&
     call error('purr',&
     'mode conversion between nin and nout not allowed',' ')
   call openz(nendf,0)
   call openz(nin,0)
   call openz(nout,1)
   nscr=10
   call openz(-nscr,1)
   maxscr=20000
   allocate(a(maxscr))
   iscr=1
   write(nsyso,'(&
     &'' unit for input endf tape ............. '',i10/&
     &'' unit for input pendf tape ............ '',i10/&
     &'' unit for output pendf tape ........... '',i10)')&
     nendf,nin,nout
   call uwtab2
   call repoz(nin)
   call repoz(nendf)
   call repoz(nout)
   call tpidio(nin,nout,0,a(iscr),nb,nw)
   call tpidio(nendf,0,0,a(iscr),nb,nw)

   ! set random number sequence
   kk=-101
   ez=rann(kk)

   !--loop over requested materials
  110 continue
   iprint=1
   nunx=0
   ntemp=1
   nsigz=1
   init=0
   read(nsysi,*) matd,ntemp,nsigz,nbin,nladr,iprint,nunx
   if (matd.eq.0) go to 400
   call repoz(-nscr)
   if (allocated(temp)) then
      deallocate(temp)
      deallocate(sigz)
      deallocate(sigu)
      deallocate(sigpl)
   endif
   allocate(temp(ntemp))
   allocate(sigz(nsigz))
   allocate(sigu(5,nsigz,ntemp))
   allocate(sigpl(200,5,nsigz))
   read(nsysi,*) (temp(i),i=1,ntemp)
   read(nsysi,*) (sigz(i),i=1,nsigz)
   write(nsyso,'(&
     &'' temperatures ......................... '',1p,e10.3)')&
     temp(1)
   if (ntemp.gt.1) write(nsyso,'(40x,1p,e10.3)') (temp(i),i=2,ntemp)
   write(nsyso,'(&
     &'' sigma zero values .................... '',&
     &2x,''infinity'')')
   if (nsigz.gt.1) write(nsyso,'(40x,1p,e10.3)') (sigz(i),i=2,nsigz)
   write(nsyso,'(&
     &'' number of probability bins ........... '',i10/&
     &'' number of resonance ladders .......... '',i10/&
     &'' print option (0=min, 1=max) .......... '',i10/&
     &'' no. of energy points (0=all) ......... '',i10)')&
     nbin,nladr,iprint,nunx
   if (nbin.lt.15) call error('purr','nbin should be 15 or more',' ')

   !--allocate storage for ladders and tables
   allocate(arry(jx))
   if (allocated(tabl)) then
      deallocate(tabl)
      deallocate(tval)
      deallocate(fis)
      deallocate(cap)
      deallocate(els)
   endif
   allocate(tabl(nbin,5,ntemp))
   allocate(tval(nbin,ntemp))
   allocate(er(nermax))
   allocate(gnr(nermax))
   allocate(gfr(nermax))
   allocate(ggr(nermax))
   allocate(gxr(nermax))
   allocate(gt(nermax))
   allocate(es(nsamp))
   allocate(xs(nsamp))
   allocate(fis(ntemp,nsamp))
   allocate(cap(ntemp,nsamp))
   allocate(els(ntemp,nsamp))
   arry=0.
   xs=0

   !--process this material
   call timer(time)
   write(nsyso,'(/'' mat = '',i4,58x,f8.1,''s'')') matd,time
   write(nsyse,'(/'' mat = '',i4,58x,f8.1,''s'')') matd,time

   !--read resonance parameters from file 2 on endf tape.
   !--read background cross sections from file 3
   nunr=0
   call findf(matd,1,451,nendf)
   call contio(nendf,0,0,a(iscr),nb,nw)
   za=c1h
   awr=c2h
   lrp=l1h
   if (lrp.le.0) go to 370
   call contio(nendf,0,0,a(iscr),nb,nw)
   if (n1h.ne.0) then
      iverf=4
   else if (n2h.eq.0) then
      iverf=5
   else
      iverf=6
   endif
   call findf(matd,2,151,nendf)
   call contio(nendf,0,0,a(iscr),nb,nw)
   call rdf2un(a)
   if (nunr.eq.0) go to 380
   maxtst=max(12+nsigz+nunr*(1+5*nsigz),12+(1+6*nbin)*nunr)
   if (maxtst.gt.maxscr) then
      write(strng1,'(''maxscr is too small, increase to at least '',&
                     &i7,'' words'')')maxtst
      call error('purr',strng1,' ')
   endif
   call findf(matd,3,1,nendf)
   call rdf3un(a)

   !--read in the total and partial heating cross sections
   if (allocated(heat)) deallocate(heat)
   h=0
   allocate(heat(4,nunr,ntemp))
   heat(:,:,:)=zero
   ihave=0
   call rdheat(a,heat,eunr,temp,ntemp,nunr,ihave,matd)
   if (ihave.eq.0) call mess('purr',&
     'no heating found on pendf',&
     'ur heating set to zero')
   if (ihave.eq.1) call mess('purr',&
     'no partial heating xsecs found on pendf',&
     'ur heating will not selfshield')

   !--compute unresolved resonance cross-sections
   !--for all grid points, temperatures, and sigma zero values.
   if ((nunx.eq.0).or.(nunx.gt.nunr)) nunx=nunr
   if (nunx.eq.1) nunx=2
   nstep=nunr/(nunx-1)
   if (nstep.lt.1) nstep=1
   ie=1-nstep
   mflg1=0
   mflg2=0
   do while (ie.lt.nunr)
      ie=ie+nstep
      if (ie.gt.nunr) ie=nunr
      ez=eunr(ie)
      bkgz(1)=sb(ie)
      bkgz(2)=sb(nunr+ie)
      bkgz(3)=sb(nunr*2+ie)
      bkgz(4)=sb(nunr*3+ie)
      sigx=bkgz(1)-bkgz(2)-bkgz(3)-bkgz(4)
      if (sigx.lt.sigmin) sigx=0
      call unresx(ez,tref)
      call unrest(bkgz,sigz,sigu,tabl,tval,er,gnr,gfr,ggr,gxr,gt,&
        es,xs,cap,fis,els,temp,nsigz,ntemp,nbin,nsamp)

      !--purr can run for a long time ... keep the user informed
      call timer(time)
      write(nsyse,'(2x,i5," of ",i5,&
            &" loops done for all temps & sig0s.",19x,f8.1,"s")')&
            &ie,nunx,time

      !--write bondarenko data and probability table to scratch file.
      write(nscr) ez
      write(nscr) (((tabl(i,j,k),i=1,nbin),j=1,5),k=1,ntemp)
      write(nscr) (((sigu(i,j,k),i=1,5),j=1,nsigz),k=1,ntemp)
   enddo

   !--check for previous mt152 and mt153
   new152=1
   new153=1
   call findf(matd,2,151,nin)
   call tosend(nin,0,0,a(iscr))
   call contio(nin,0,0,a(iscr),nb,nw)
   if (mth.eq.152) then
      new152=0
      call tosend(nin,0,0,a(iscr))
      call contio(nin,0,0,a(iscr),nb,nw)
   endif
   if (mth.eq.153) new153=0
   call findf(matd,1,451,nin)

   !--write data to new pendf tape for all temperatures
   do 350 it=1,ntemp
   temz=temp(it)
   nx=0
   call contio(nin,0,0,a(iscr),nb,nw)
   if (iverf.lt.5) then
      nx=n2h
      if (nx.gt.0.and.new152.gt.0) a(iscr+5)=a(iscr+5)+1
      if (nx.gt.0.and.new153.gt.0) a(iscr+5)=a(iscr+5)+1
   endif
   call contio(0,nout,0,a(iscr),nb,nw)
   if (iverf.ge.5) call contio(nin,nout,0,a(iscr),nb,nw)
   if (iverf.ge.6) call contio(nin,nout,0,a(iscr),nb,nw)
   call hdatio(nin,0,0,a(iscr),nb,nw)
   a(iscr)=temz
   if (iverf.ge.5) then
      nx=n2h
      if (nx.gt.0.and.new152.gt.0) a(iscr+5)=a(iscr+5)+1
      if (nx.gt.0.and.new153.gt.0) a(iscr+5)=a(iscr+5)+1
   endif
   call hdatio(0,nout,0,a(iscr),nb,nw)
   do while (nb.ne.0)
      call moreio(nin,nout,0,a(iscr),nb,nw)
   enddo
   nw=nx
   if (nw.eq.0) go to 250
   nc=6*(nx+new152+new153)
   if (allocated(c)) deallocate(c)
   allocate(c(nc))
   nd=6*nx
   if (allocated(d)) deallocate(d)
   allocate(d(nd))
   call dictio(nin,0,0,d,nb,nw)
   j=0
   newmat=0
   ncds=nsigz+nunx*(1+5*nsigz)
   ncds=2+(ncds-1)/6
   nc153=(1+6*nbin)*nunx
   nc153=2+(nc153-1)/6
   do 240 i=1,nw,6
   if (newmat.gt.0) go to 220
   mfd=nint(d(i+2))
   mtd=nint(d(i+3))
   if (mfd.lt.2) go to 220
   if (mfd.eq.2.and.mtd.eq.151) go to 220
   if (mfd.eq.2.and.mtd.eq.152) go to 240
   if (mfd.eq.2.and.mtd.eq.153) go to 240
   newmat=1
   c(1+j)=0
   c(2+j)=0
   c(3+j)=2
   c(4+j)=152
   c(5+j)=ncds
   c(6+j)=0
   c(7+j)=0
   c(8+j)=0
   c(9+j)=2
   c(10+j)=153
   c(11+j)=nc153
   c(12+j)=0
   j=j+12
  220 continue
   do k=1,6
      c(j+k)=d(k+i-1)
   enddo
   j=j+6
  240 continue
   nwd=j/6
   call dictio(0,nout,0,c,nb,nwd)
   deallocate(c)
   deallocate(d)
   ! copy through file 2, mt 151
  250 continue
   call tofend(nin,nout,0,a(iscr))
   call tosend(nin,nout,0,a(iscr))
   call contio(nin,0,0,a(iscr),nb,nw)
   if (mth.eq.152) then
      call tosend(nin,0,0,a(iscr))
      call contio(nin,0,0,a(iscr),nb,nw)
   endif
   if (mth.eq.153) then
      call tosend(nin,0,0,a(iscr))
      call contio(nin,0,0,a(iscr),nb,nw)
   endif

   !--write new mt152 on output tape
   !--for bondarenko self-shielded cross sections
   math=matd
   mfh=2
   mth=152
   l=iscr
   a(l)=za
   a(l+1)=awr
   a(l+2)=lssf
   a(l+3)=0
   a(l+4)=0
   a(l+5)=intunr
   a(l+6)=temz
   a(l+7)=0
   a(l+8)=5
   a(l+9)=nsigz
   a(l+10)=nsigz+nunx*(1+5*nsigz)
   a(l+11)=nunx
   l=l+11
   do i=1,nsigz
      l=l+1
      a(l)=sigz(i)
   enddo
   call repoz(-nscr)
   ie=1-nstep
   do while (ie.lt.nunr)
      ie=ie+nstep
      if (ie.gt.nunr) ie=nunr
      read(nscr) ez
      read(nscr) (((tabl(i,j,k),i=1,nbin),j=1,5),k=1,ntemp)
      read(nscr) (((sigu(i,j,k),i=1,5),j=1,nsigz),k=1,ntemp)
      l=l+1
      a(l)=eunr(ie)
      do ix=1,5
         do is=1,nsigz
            l=l+1
            a(l)=sigfig(sigu(ix,is,it),7,0)
         enddo
      enddo
   enddo
   call contio(0,nout,0,a(iscr),nb,nw)
   lstart=iscr+6
   nwds=l-iscr-5
   nw=nwds
   if (nw.gt.npage+6) nw=npage+6
   call listio(0,nout,0,a(lstart),nb,nw)
  290 continue
   lstart=lstart+nw
   nwds=nwds-nw
   nw=nwds
   if (nw.gt.npage) nw=npage
   if (nb.eq.0) go to 300
   call moreio(0,nout,0,a(lstart),nb,nw)
   go to 290
  300 continue
   call asend(nout,0)

   !--write new mt153 on output tape
   !--for probability tables
   math=matd
   mfh=2
   mth=153
   n=iscr
   a(n)=za
   a(n+1)=awr
   a(n+2)=iinel
   a(n+3)=iabso
   a(n+4)=intunr
   a(n+5)=nbin
   a(n+6)=temz
   a(n+7)=0
   a(n+8)=lssf
   a(n+9)=0
   a(n+10)=(1+6*nbin)*nunx
   a(n+11)=nunx
   n=n+11
   call repoz(-nscr)
   ie=1-nstep
   do while (ie.lt.nunr)
      ie=ie+nstep
      if (ie.gt.nunr) ie=nunr
      read(nscr) ez
      read(nscr) (((tabl(i,j,k),i=1,nbin),j=1,5),k=1,ntemp)
      read(nscr) (((sigu(i,j,k),i=1,5),j=1,nsigz),k=1,ntemp)
      n=n+1
      a(n)=sigfig(eunr(ie),7,0)
      n1=n
      do i=1,5
         do j=1,nbin
            n=n+1
            a(n)=sigfig(tabl(j,i,it),7,0)
         enddo
      enddo
      do j=1,nbin
         n=n+1
         a(n)=0
      enddo
      do i=2,5
         do j=1,nbin
            l=n1+j+nbin*(i-1)
            if (lssf.eq.1) then
               if (sigu(i-1,1,1).ne.0) then
                  a(l)=a(l)/sigu(i-1,1,1)
               else
                  a(l)=1
               endif
            endif
            a(l)=sigfig(a(l),7,0)
         enddo
      enddo
      do j=1,nbin
         l=n1+j+5*nbin
         if (ihave.eq.1) then
            if (lssf.eq.1) then
               a(l)=1
            else
               a(l)=heat(1,ie,it)/sigu(1,1,1)
               endif
         else if (ihave.eq.2) then
            a(l)=a(l)+(heat(1,ie,it)-heat(2,ie,it)-heat(3,ie,it)-heat(4,ie,it))
            h=heat(2,ie,it)
            if (lssf.eq.1) then
               h=h*a(n1+j+2*nbin)
            else if (sigu(2,1,1).ne.zero) then
               h=h*a(n1+j+2*nbin)/sigu(2,1,1)
            endif
            a(l)=a(l)+h
            h=heat(3,ie,it)
            if (lssf.eq.1) then
               h=h*a(n1+j+3*nbin)
            else if (sigu(3,1,1).ne.zero) then
               h=h*a(n1+j+3*nbin)/sigu(3,1,1)
            endif
            a(l)=a(l)+h
            h=heat(4,ie,it)
            if (lssf.eq.1) then
               h=h*a(n1+j+4*nbin)
            elseif (sigu(4,1,1).ne.zero) then
               h=h*a(n1+j+4*nbin)/sigu(4,1,1)
            endif
            a(l)=a(l)+h
            if (a(n1+j+nbin).ne.zero) a(l)=a(l)/a(n1+j+nbin)
            if (lssf.eq.1) a(l)=a(l)/heat(1,ie,it)
         endif
      enddo
   enddo
   call contio(0,nout,0,a(iscr),nb,nw)
   lstart=iscr+6
   nwds=6*nbin*nunx
   nw=nwds
   if (nw.gt.npage+6) nw=npage+6
   call listio(0,nout,0,a(lstart),nb,nw)
  330 continue
   lstart=lstart+nw
   nwds=nwds-nw
   nw=nwds
   if (nw.gt.npage) nw=npage
   if (nb.eq.0) go to 340
   call moreio(0,nout,0,a(lstart),nb,nw)
   go to 330
  340 continue
   call asend(nout,0)
   ! write file end card for file 2
   call afend(nout,0)
   ! copy rest of material to output tape
   call contio(nin,nout,0,a(iscr),nb,nw)
   call tomend(nin,nout,0,a(iscr))
  350 continue

   !--write report of calculation
   call timer(time)
   write(nsyso,'(/'' generated cross sections at '',i3,&
     &'' points'',30x,f8.1,''s'')') nunx,time
   deallocate(sb)
   go to 110

   !--mat has no resonance parameters. copy as is to nout.
  370 continue
   write(strng1,&
     '(''mat'',i5,'' has no resonance parameters'')') matd
   write(strng2,'(''copy as is to nout'')')
   call mess('purr',strng1,strng2)
   go to 385

   !--mat has no unresolved.  copy as is to nout
  380 continue
   write(strng1,&
     '(''mat'',i5,'' has no unresolved parameters'')') matd
   write(strng2,'(''copy as is to nout'')')
   call mess('purr',strng1,strng2)
  385 continue
   call findf(matd,1,451,nin)
   it=0
  390 continue
   call contio(nin,0,0,a(iscr),nb,nw)
   if (math.ne.matd) go to 110
   call contio(0,nout,0,a(iscr),nb,nw)
   call tomend(nin,nout,0,a(iscr))
   it=it+1
   if (it.lt.ntemp) go to 390
   go to 110

   !--finished with purr
  400 continue
   !--print bondarenko cross sections for plotting
   !--uncomment the following lines to activate
    !do k=1,5
    !   do j=1,nsigz
    !      if (k.eq.1) write(nsyso,&
    !        '(/'' p0 total, sigz='',1p,e10.2)') sigz(j)
    !      if (k.eq.2) write(nsyso,&
    !        '(/'' elastic, sigz='',1p,e10.2)') sigz(j)
    !      if (k.eq.3) write(nsyso,&
    !        '(/'' fission, sigz='',1p,e10.2)') sigz(j)
    !      if (k.eq.4) write(nsyso,&
    !        '(/'' capture, sigz='',1p,e10.2)') sigz(j)
    !      if (k.eq.5) write(nsyso,&
    !        '(/'' p1 total, sigz='',1p,e10.2)') sigz(j)
    !      do i=1,ipl
    !         write(nsyso,'(1p,2e11.4)') epl(i),sigpl(i,k,j)
    !      enddo
    !   enddo
    !enddo

   deallocate(arry)
   deallocate(temp)
   deallocate(sigz)
   deallocate(sigu)
   deallocate(sigpl)
   if (allocated(bval)) then
      deallocate(bval)
      deallocate(tmin)
      deallocate(tmax)
      deallocate(tsum)
   endif
   if (allocated(tabl)) then
      deallocate(tabl)
      deallocate(tval)
      deallocate(er)
      deallocate(gnr)
      deallocate(gfr)
      deallocate(ggr)
      deallocate(gxr)
      deallocate(gt)
      deallocate(es)
      deallocate(xs)
      deallocate(fis)
      deallocate(cap)
      deallocate(els)
   endif
   if (allocated(heat)) deallocate(heat)
   if (allocated(c)) deallocate(c)
   if (allocated(d)) deallocate(d)

   call atend(nout,0)
   call repoz(nout)
   call closz(nendf)
   call closz(nin)
   call closz(nout)
   call closz(-nscr)
   call timer(time)
   write(nsyso,'(/69x,f8.1,''s''/&
     &7(''**********''),''*******'')') time
   return
   end subroutine purr

   subroutine rdf2un(a)
   !-------------------------------------------------------------------
   ! Read and/or copy unresolved resonance parameters from File 2.
   !-------------------------------------------------------------------
   use util ! provides error,mess,sigfig
   use endf ! provides endf routines and variables
   ! externals
   real(kr)::a(*)
   ! internals
   integer::indep,nis,inow,is,nb,nw,lfw,ner,ier,lru,lrf
   integer::nls,l,njs,j,mode,ll,jnow,n,i,ne,inow1,ir,k,jj
   integer::ist,loc,nwds,ione,itwo,ithree,ii,iovlp,iscr
   real(kr)::elr,ehr,abn,el,eh,et,spi,ay,awri,enow,elast,enext
   real(kr)::test,enut
   integer,parameter::ngridu=78
   real(kr),dimension(78)::egridu=(/&
     1.0e1_kr,1.25e1_kr,1.5e1_kr,1.7e1_kr,2.0e1_kr,2.5e1_kr,3.0e1_kr,&
     3.5e1_kr,4.0e1_kr,5.0e1_kr,6.0e1_kr,7.2e1_kr,8.5e1_kr,1.0e2_kr,&
     1.25e2_kr,1.5e2_kr,1.7e2_kr,2.0e2_kr,2.5e2_kr,3.0e2_kr,3.5e2_kr,&
     4.0e2_kr,5.0e2_kr,6.0e2_kr,7.2e2_kr,8.5e2_kr,1.0e3_kr,1.25e3_kr,&
     1.5e3_kr,1.7e3_kr,2.0e3_kr,2.5e3_kr,3.0e3_kr,3.5e3_kr,4.0e3_kr,&
     5.0e3_kr,6.0e3_kr,7.2e3_kr,8.5e3_kr,1.0e4_kr,1.25e4_kr,1.5e4_kr,&
     1.7e4_kr,2.0e4_kr,2.5e4_kr,3.0e4_kr,3.5e4_kr,4.0e4_kr,5.0e4_kr,&
     6.0e4_kr,7.2e4_kr,8.5e4_kr,1.0e5_kr,1.25e5_kr,1.5e5_kr,1.7e5_kr,&
     2.0e5_kr,2.5e5_kr,3.0e5_kr,3.5e5_kr,4.0e5_kr,5.0e5_kr,6.0e5_kr,&
     7.2e5_kr,8.5e5_kr,1.e6_kr,1.25e6_kr,1.5e6_kr,1.7e6_kr,2.0e6_kr,&
     2.5e6_kr,3.0e6_kr,3.5e6_kr,4.0e6_kr,5.0e6_kr,6.0e6_kr,7.2e6_kr,&
     8.5e6_kr/)
   real(kr),parameter::etop=5.e6_kr
   real(kr),parameter::small=1.e-12_kr
   real(kr),parameter::wide=1.26e0_kr

   iscr=1
   indep=0
   elr=0
   ehr=etop
   nunr=1
   eunr(1)=etop

   !--head card read in calling program
   nis=n1h
   nsect=0
   inow=1

   !--loop over isotopes
   do 110 is=1,nis
   call contio(nendf,0,0,a(iscr),nb,nw)
   abn=c2h
   lfw=l2h
   ner=n1h
   ier=0

   !--loop over all energy ranges
  115 continue
   call contio(nendf,0,0,a(iscr),nb,nw)
   if (mth.eq.0) go to 300
   ier=ier+1
   el=sigfig(c1h,7,0)
   eh=sigfig(c2h,7,0)
   lru=l1h
   lrf=l2h
   nro=n1h
   naps=n2h
   if (lru.eq.2) go to 140
   if (eh.gt.elr) elr=eh

   !--read through this sub-section
   call contio(nendf,0,0,a(iscr),nb,nw)
   nls=n1h
   if (lrf.eq.4.or.lrf.eq.7) then
      call listio(nendf,0,0,a(iscr),nb,nw)
      do while (nb.ne.0)
         call moreio(nendf,0,0,a(iscr),nb,nw)
      enddo
   endif
   do l=1,nls
      if (lrf.ne.4) then
         call listio(nendf,0,0,a(iscr),nb,nw)
         do while (nb.ne.0)
            call moreio(nendf,0,0,a(iscr),nb,nw)
         enddo
         if (lrf.eq.7) then
            call listio(nendf,0,0,a(iscr),nb,nw)
            do while (nb.ne.0)
               call moreio(nendf,0,0,a(iscr),nb,nw)
            enddo
         endif
      else
         call contio(nendf,0,0,a(iscr),nb,nw)
         njs=n1h
         do j=1,njs
            call listio(nendf,0,0,a(iscr),nb,nw)
            do while (nb.ne.0)
               call moreio(nendf,0,0,a(iscr),nb,nw)
            enddo
         enddo
      endif
   enddo
   if (ier.eq.ner) go to 110
   go to 115

   !--read and copy unresolved resonance parameters
  140 continue
   nsect=nsect+1
   isot(nsect)=is
   abnt(nsect)=abn
   elt(nsect)=el
   eht(nsect)=eh
   mode=lrf+4*(lru-1)
   modet(nsect)=mode
   ibaset(nsect)=inow
   if (eh.lt.ehr) ehr=eh
   ! add el and eh to list of nodes
   ! shade nodes around discontinuities
   et=sigfig(el,7,-1)
   call ilist2(et,eunr,nunr)
   et=sigfig(el,7,+1)
   call ilist2(et,eunr,nunr)
   et=sigfig(eh,7,-1)
   call ilist2(et,eunr,nunr)
   et=sigfig(eh,7,+1)
   call ilist2(et,eunr,nunr)

   ! if present, read and store the energy-dependent
   ! scattering radius tab1 data.
   if (nro.eq.1) then
      call tab1io(nendf,0,0,arry(inow),nb,nw)
      jj=inow+nw
      do while (nb.ne.0)
         call moreio(nendf,0,0,arry(jj),nb,nw)
         jj=jj+nw
         if (jj.gt.jx) then
            call error('rdf2un','storage in a exceeded',' ')
         endif
      enddo
      inow=jj
   endif

   !--branch to specified representation
   if (lrf.ne.2) then

      !--all parameters independent of energy
      if (lfw.ne.1) then
         call contio(nendf,0,0,a(iscr),nb,nw)
         spi=c1h
         ay=c2h
         lssf=l1h
         nls=n1h
         ! loop over l states
         do l=1,nls
            call listio(nendf,0,0,a(iscr),nb,nw)
            awri=c1h
            ll=l1h
            njs=n2h
            if (l.eq.1) then
               arry(inow)=awri
               arry(inow+1)=ay
               arry(inow+2)=spi
               arry(inow+3)=lfw
               arry(inow+4)=nls
               arry(inow+5)=0
               inow=inow+6
            endif
            arry(inow)=njs
            arry(inow+1)=ll
            inow=inow+2
            jnow=iscr+5
            do n=1,njs
               do i=1,5
                  arry(inow)=a(jnow+i)
                  inow=inow+1
               enddo
               jnow=jnow+6
            enddo
         enddo
         if (inow.gt.jx) call error('rdf2un','storage exceeded',' ')
         indep=1
         intunr=5

      !--fission widths energy dependent
      else
         call listio(nendf,0,0,a(iscr),nb,nw)
         spi=c1h
         ay=c2h
         lssf=l1h
         ne=n1h
         nls=n2h
         arry(inow+1)=ay
         arry(inow+2)=spi
         arry(inow+3)=lfw
         arry(inow+4)=nls
         arry(inow+5)=ne
         ! store fission width energy points
         inow1=inow+5
         do i=1,ne
            arry(inow1+i)=a(iscr+i+5)
            enow=a(iscr+i+5)
            arry(inow1+i)=enow
            ! use as energy nodes
            if (i.gt.1.and.i.lt.ne) then
               if (nunr.ge.meunr) call error('rdunf2',&
                 'too many ur energy points',' ')
               if (enow.ge.el.and.enow.le.eh) call ilist2(enow,eunr,nunr)
            endif
         enddo
         ! loop over l states
         do l=1,nls
            call contio(nendf,0,0,a(iscr),nb,nw)
            if (l.eq.1) then
               arry(inow)=c1h
               inow=inow+6+ne
            endif
            njs=n1h
            arry(inow)=njs
            arry(inow+1)=l1h
            inow=inow+2
            ! loop over j states
            do j=1,njs
               call listio(nendf,0,0,a(iscr),nb,nw)
               ne=n1h-6
               arry(inow)=l2h
               ! store parameters
               do ir=1,5
                  arry(inow+ir)=a(iscr+5+ir)
               enddo
               inow1=inow+5
               ! store fission widths
               do k=1,ne
                  arry(inow1+k)=a(iscr+11+k)
               enddo
               inow=inow+6+ne
            enddo
         enddo
         if (inow.gt.jx) call error('rdf2un','storage exceeded',' ')
         intunr=5
      endif

   !--all parameters energy dependent
   else
      ist=inow
      ! first cont record
      call contio(nendf,0,0,a(iscr),nb,nw)
      lssf=l1h
      nls=n1h
      arry(inow+1)=c2h
      arry(inow+2)=c1h
      arry(inow+3)=n1h
      inow=inow+4
      ! loop over l states
      do l=1,nls
         call contio(nendf,0,0,a(iscr),nb,nw)
         if (l.eq.1) arry(ist)=c1h
         njs=n1h
         arry(inow)=njs
         arry(inow+1)=l1h
         inow=inow+2
         ! loop over j states
         do j=1,njs
            call listio(nendf,0,0,a(iscr),nb,nw)
            loc=iscr+nw
            do while (nb.ne.0)
               call moreio(nendf,0,0,a(loc),nb,nw)
               loc=loc+nw
            enddo
            ne=n2h
            intunr=l1h
            arry(inow)=c1h
            arry(inow+1)=l1h
            arry(inow+2)=ne
            ! store numbers of degrees of freedom
            do k=3,6
               arry(inow+k)=a(iscr+5+k)
            enddo
            inow=inow+7
            jnow=iscr+11
            ! store parameters
            do n=1,ne
               do k=1,6
                  if (a(jnow+k).le.0.) a(jnow+k)=small
                  arry(inow+k-1)=a(jnow+k)
               enddo
               inow=inow+6
               if (n.ne.1.and.n.ne.ne) then
                  if (l.eq.1) then
                     if (j.eq.1) then
                        ! store energy as a node
                        enow=a(jnow+1)
                        if (nunr.ge.meunr) call error('rdunf2',&
                          'too many ur energy points',' ')
                        if (enow.ge.el.and.enow.le.eh)&
                          call ilist2(enow,eunr,nunr)
                     endif
                  endif
               endif
               jnow=jnow+6
            enddo
         enddo
      enddo
      if (inow.gt.jx) call error('rdf2un','storage exceeded',' ')
   endif
   if (ier.lt.ner) go to 115
  110 continue
  300 continue
   nwds=inow
   nunr=nunr-1
   if (nunr.eq.0) return

   !--remove first and last energy nodes
   !--flag resolved-unresolved overlap energies
   nunr=nunr+1
   i=1
   elast=eunr(2)
   ione=0
   do while (ione.eq.0)
      i=i+1
      enext=eunr(1+i)
      test=etop
      if (enext.ge.test) then
         ione=1
      else
         if (enext.ge.wide*elast.or.indep.ne.0) then
            et=elast
            itwo=0
            do while (itwo.eq.0)
               ithree=0
               ii=0
               do while (ii.lt.ngridu.and.ithree.eq.0)
                  ii=ii+1
                  enut=egridu(ii)
                  if (enut.gt.et+et/100) ithree=1
               enddo
               if (ithree.eq.0) enut=enext
               et=enut
               if (et.ge.enext) then
                  itwo=1
               else
                  if (nunr.ge.meunr) call error('rdunf2',&
                    'too many ur energy points',' ')
                  call ilist2(et,eunr,nunr)
                  i=i+1
               endif
            enddo
         endif
         elast=eunr(1+i)
      endif
   enddo
   iovlp=0
   nunr=nunr-2
   do i=2,nunr
      et=eunr(i)
      if (et.lt.elr) et=-et
      if (et.gt.ehr) et=-et
      if (et.lt.0.) iovlp=1
      eunr(i-1)=et
   enddo
   nunr=nunr-1
   if (iovlp.eq.1) call mess('rdf2un',&
     'resolved-unresolved overlap energies',&
     'marked with minus signs')
   return
   end subroutine rdf2un

   subroutine ilist2(e,elist,nlist)
   !-------------------------------------------------------------------
   ! Insert a new energy into a list of enegies in ascending order.
   ! Omit duplicate values.  Initial list must be primed with one
   ! energy larger than any others to be added.
   !-------------------------------------------------------------------
   ! externals
   integer::nlist
   real(kr)::e,elist(*)
   ! internals
   integer::jl,i,idone,j

   jl=0
   i=0
   idone=0
   do while (i.lt.nlist.and.idone.eq.0)
      i=i+1
      if (e.le.elist(i)) then
         if (e.ne.elist(i)) then
            jl=nlist-i+1
            do j=1,jl
               elist(nlist-j+2)=elist(nlist-j+1)
            enddo
            elist(i)=e
         endif
         idone=1
      endif
   enddo
   if (idone.eq.0) then
      jl=1
      elist(nlist+1)=e
   endif
   if (jl.gt.0) nlist=nlist+1
   return
   end subroutine ilist2

   subroutine rdf3un(a)
   !-------------------------------------------------------------------
   ! Read unresolved region background cross sections from File 3.
   ! Compute the background cross sections on the unresolved
   ! energy grid and store for use in unresx.  Also, analyze the
   ! competitive reactions, if any.
   !-------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use util ! provides mess
   use endf ! provides endf routines and variables
   ! externals
   real(kr)::a(*)
   ! internals
   integer::iscr,nthr,nb,nw,idis,ix,jthr,ibase,ie,icx,kthr,mtc,mtmax
   real(kr)::e,enext,bkg,sigx,ecomp,total
   real(kr),dimension(:),allocatable::thr
   character(60)::strng1,strng2
   integer,dimension(4),parameter::mtx=(/1,2,18,102/)
   real(kr),parameter::up=1.00001e0_kr
   real(kr),parameter::dn=.99999e0_kr
   real(kr),parameter::small=1.e-5_kr
   real(kr),parameter::tol=1.e-6_kr

   !--allocate storage
   iscr=1
   nthr=1000
   nw=2*nthr
   allocate(thr(nw))
   nw=nunr*4
   allocate(sb(nunr*4))

   !--initialise background array for total, elastic, fission and capture
   sb(1:nunr*4)=0

   !--loop over reactions
   !--saving resonance cross sections
   !--and reaction thresholds
   jthr=0
   ix=1
   ibase=1+nunr*(ix-1)-1
   mtmax=891 ! continuum n,2n
   call contio(nendf,0,0,a(iscr),nb,nw)
   do while (mfh.eq.3.and.mth.le.mtmax)
      if (mth.lt.6.or.(mth.gt.10.and.mth.lt.12).or. &
         & (mth.gt.15.and.mth.lt.43).or. &
         & (mth.gt.43.and.mth.lt.46).or. &
         & (mth.gt.49.and.mth.lt.92).or. &
         & (mth.gt.100.and.mth.lt.110).or. &
         & (mth.gt.110.and.mth.lt.118).or. &
         & (mth.gt.150.and.mth.le.200).or.mth.ge.600) then
         !--get reaction threshold and store data
         e=0
         call gety1(e,enext,idis,bkg,nendf,a(iscr))
         if (mth.ge.4) then
            jthr=jthr+1
            thr(1+2*(jthr-1))=mth
            thr(2+2*(jthr-1))=enext
         endif
         !--skip unresolved resonance reaction if required
         do while (mtx(ix).ne.102.and.mtx(ix).lt.mth)
            ix=ix+1
         enddo
         !--store data is this is one we need
         if (mth.eq.mtx(ix)) then
            ibase=1+nunr*(ix-1)-1
            do ie=1,nunr
               e=abs(eunr(ie))
               if (ie.eq.1) e=up*e
               if (ie.eq.nunr) e=dn*e
               call gety1(e,enext,idis,sb(ibase+ie),nendf,a(iscr))
            enddo
         endif
      endif
      call tosend(nendf,0,0,a(iscr))
      call contio(nendf,0,0,a(iscr),nb,nw)
   enddo

   !--output lssf flag
   write(nsyso,'('' evaluation lssf equal to'',i2)') lssf

   !--check for competitive reactions
   icx=0
   do ie=1,nunr
      sigx=sb(ie)-sb(nunr+ie)-sb(2*nunr+ie)-sb(3*nunr+ie)
      if (icx.eq.0.and.sigx.gt.small) icx=ie
   enddo
   iinel=-1 ! default is no competition
   iabso=-1 ! default is no competition
   if (icx.gt.0) then
      ecomp=eunr(icx)
      write(nsyso,'(/'' competition starts at'',1p,e12.4)') ecomp
      do kthr=1,jthr
         mtc=nint(thr(1+2*(kthr-1)))
         if (mtc.gt.4.and.mtc.ne.18.and.mtc.ne.19.and.mtc.ne.102) then
            if (up*thr(2+2*(kthr-1)).lt.eunr(nunr)) then
               write(nsyso,'(''   ur competes with mt'',i3)') mtc
               if (mtc.ge.51.and.mtc.le.91) then
                  if (iinel.lt.0) then
                     iinel=mtc ! use mtc as flag if it is the only one
                  else
                     iinel=4   ! more than one in competition, use mt4 instead
                  endif
               else
                  if (iabso.lt.0) then
                     iabso=mtc ! use mtc as flag if it is the only one
                  else
                     iabso=0   ! more than one in competition, use 0
                  endif
               endif
            endif
         endif
      enddo
      write(nsyso,'('' inelastic competition flag set to '',i3/&
        &'' absorption competition flag set to '',i3)') iinel,iabso
   endif

   !--sanity check for lssf=1, derive corresponding sb(i) array
   if (lssf.gt.0) then
      do ie=1,nunr
         total=sb(ie)
         sb(ie)=sb(ie)-sb(ie+nunr)-sb(ie+2*nunr)-sb(ie+3*nunr)
         if (sb(ie).gt.tol*total) then
            !--this can only happen if there is competition and above ecomp
            if (iinel.lt.0.and.iabso.lt.0) then
               write(strng1,&
                 '(''total xs greater than its components at e=''&
                 &,1p,e12.4)') eunr(ie)
               write(strng2,'(''check evaluation file'')')
               call mess('purr',strng1,strng2)
               sb(ie)=0
            elseif (eunr(ie).lt.ecomp) then
               write(strng1,&
                 '(''total xs greater than its components at e=''&
                 &,1p,e12.4)') eunr(ie)
               write(strng2,'(''check evaluation file below competition'')')
               call mess('purr',strng1,strng2)
               sb(ie)=0
            endif
         else
            !--either roundoff or total is smaller than its components
            if (sb(ie).lt.-tol*total) then
               write(strng1,&
                 '(''total xs less than its components at e=''&
                 &,1p,e12.4)') eunr(ie)
               write(strng2,'(''check evaluation file'')')
               call mess('purr',strng1,strng2)
            endif
            sb(ie)=0
         endif
         sb(ie+nunr)=0
         sb(ie+2*nunr)=0
         sb(ie+3*nunr)=0
      enddo
   endif

   return
   end subroutine rdf3un

   subroutine rdheat(a,heat,eunr,temp,ntemp,nunr,ihave,matd)
   !-------------------------------------------------------------------
   ! Read the total heating (MT=301) and partial heating cross
   ! sections for elastic (302), fission (318), and capture (402)
   ! from the pendf tape on the unresolved energy grid.
   !-------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use util ! provides repoz
   use endf ! provides endf routines and variables
   ! externals
   integer::ntemp,nunr,ihave,matd
   real(kr)::a(*),heat(4,nunr,ntemp),eunr(nunr),temp(ntemp)
   ! internals
   integer::iscr,nb,nw,itemp,idis,ix,ie
   real(kr)::enext,h
   iscr=1

   !--loop over the temperatures
   call repoz(nin)
   call tpidio(nin,0,0,a(iscr),nb,nw)
   itemp=0
   do while (itemp.lt.ntemp)
      itemp=itemp+1
      call findf(matd,3,1,nin)
      ihave=0
      mth=1
      do while (mth.ne.0)
         call contio(nin,0,0,a(iscr),nb,nw)
         if (mth.ne.0) then
            ix=0
            if (mth.eq.301) ix=1
            if (mth.eq.302) ix=2
            if (mth.eq.318) ix=3
            if (mth.eq.402) ix=4
            if (ix.gt.0) then
               e=0
               call gety1(e,enext,idis,h,nin,a(iscr))
               if (ix.eq.1) ihave=1
               if (ix.gt.1) ihave=2
               do ie=1,nunr
                  e=abs(eunr(ie))
                  call gety1(e,enext,idis,h,nin,a(iscr))
                  heat(ix,ie,itemp)=h
               enddo
            endif
            call tosend(nin,0,0,a(iscr))
         endif
      enddo
      call tomend(nin,0,0,a(iscr))
   enddo
   call repoz(nin)
   return
   end subroutine rdheat

   subroutine unresx(ee,tt)
   !-------------------------------------------------------------------
   ! Construct ladder parameters, infinite-dilution cross sections,
   ! and potential scattering for this energy.
   !-------------------------------------------------------------------
   use util    ! provides error
   use endf    ! provides terpa
   use physics ! provides pi,amassn,amu,hbar,ev
   ! externals
   real(kr)::ee,tt
   ! internals
   integer::i,isect,mode,inow,lfw,nls,ne,iest,l,njs,ll,j
   integer::ifst,int,lu,mu,nu,ip,ir,iro,idx
   real(kr)::cwaven,k,t,abn,ay,spi,awri,aaa,enx
   real(kr)::rat,aw,aa,e2,ab,rho,rhoc,amux,gxx,amuf,dx,aj
   real(kr)::amun,gnox,ggx,gfx,vl,ps,gj,gnx,tp,gs,gc,gf,temp
   real(kr),parameter::rc1=.123e0_kr
   real(kr),parameter::rc2=.08e0_kr
   real(kr),parameter::third=.333333333e0_kr
   real(kr),parameter::con1=2901.34e0_kr
   cwaven=sqrt(2*amassn*amu*ev)*1.e-12_kr/hbar

   !--initialize
   e=abs(ee)
   t=tt
   spot=0
   dbarin=0
   nseq0=0
   sigi(1)=0
   sigi(2)=0
   sigi(3)=0
   sigi(4)=0
   lfw=0
   vl=0.
   ps=0.
   ne=0
   iest=0
   aa=0.

   !--find sections of resonance parameters which contribute
   do i=1,nsect
      isect=i
      abn=abnt(i)
      mode=modet(i)
      if ((e.ge.elt(i).and.e.le.eht(i)).and.&
        (mode.eq.5.or.mode.eq.6)) then

         !--retrieve starting location in storage for this section
         inow=ibaset(isect)

         !--retrieve nuclide information
         if (nro.eq.1) then
            iro=inow
            ip=2
            ir=1
            call terpa(ay,e,enx,idx,arry(iro),ip,ir)
            inow=inow+6+2*nint(arry(iro+4))+2*nint(arry(iro+5))
         else
            ay=arry(inow+1)
         endif
         awri=arry(inow)
         aaa=0
         if (naps.eq.2.and.nro.eq.1) aaa=arry(inow+1)
         spi=arry(inow+2)
         if (mode.ne.6) then
            lfw=nint(arry(inow+3))
            nls=nint(arry(inow+4))
            ne=nint(arry(inow+5))
            inow=inow+6
            ! save starting location for fission width energies
            if (lfw.ne.0) then
               iest=inow
               inow=inow+ne
            endif
         else
            nls=nint(arry(inow+3))
            inow=inow+4
         endif

         !--calculate nuclide dependent parameters
         rat=awri/(awri+1)
         aw=awri*amassn
         if (naps.eq.0) then
            aa=rc1*aw**third+rc2
         elseif (naps.eq.1) then
            aa=ay
         else if (naps.eq.2 .and. nro.eq.1) then
            aa=aaa
         else
            call error('unresx','illegal naps',' ')
         endif
         e2=sqrt(e)
         k=cwaven*rat*e2
         ab=4*pi/k**2
         rho=k*aa
         rhoc=k*ay

         !--loop over all sequences (l,j) in this section
         do l=1,nls
            njs=nint(arry(inow))
            ll=nint(arry(inow+1))
            inow=inow+2
            do j=1,njs

               !--retrieve parameters for this sequence
               if (mode.ne.6) then
                  amux=0
                  gxx=0
                  if (lfw.eq.0) then
                     amuf=0
                     dx=arry(inow)
                     aj=arry(inow+1)
                     amun=arry(inow+2)
                     gnox=arry(inow+3)
                     ggx=arry(inow+4)
                     gfx=0
                     inow=inow+5
                  else
                     amuf=arry(inow)
                     dx=arry(inow+1)
                     aj=arry(inow+2)
                     amun=arry(inow+3)
                     gnox=arry(inow+4)
                     ggx=arry(inow+5)
                     inow=inow+6
                     ! interpolate for fission width at this energy
                     ifst=inow
                     call intrf2(e,gfx,iest,ifst,ne,2,arry)
                     inow=inow+ne
                  endif
               else
                  aj=arry(inow)
                  int=nint(arry(inow+1))
                  ne=nint(arry(inow+2))
                  amux=arry(inow+3)
                  amun=arry(inow+4)
                  amuf=arry(inow+6)
                  inow=inow+7
                  ! interpolate for parameters at this energy
                  call intr2(e,dx,gxx,gnox,ggx,gfx,inow,ne,int,arry)
                  inow=inow+6*ne
               endif

               !--calculate penetrability and phase shifts
               call unfac2(ll,rho,rhoc,amun,vl,ps)
               gj=(2*aj+1)/(4*spi+2)
               gnx=gnox*vl*e2*amun

               !--compute potential scattering
               if (j.eq.1) then
                  spot=spot+abn*ab*(2*ll+1)*sin(ps)**2
               endif

               !--store ladder parameters
               dbarin=dbarin+1/dx
               lu=nint(amux)
               mu=nint(amuf)
               nu=nint(amun)
               gnx=gnx/nu
               if (mu.le.0) gfx=0
               if (lu.le.0) gxx=0
               nseq0=nseq0+1
               if (nseq0.gt.mxns0) call error('unresx',&
                 'too many sequences, increase mxns0',' ')
               csz(nseq0)=abn*ab*gj
               tp=t
               if (tp.eq.0) tp=1
               cth(nseq0)=sqrt(con1*awri/(e*tp))
               cc2p(nseq0)=cos(2*ps)
               cs2p(nseq0)=sin(2*ps)
               cgn(nseq0)=gnx
               cgg(nseq0)=ggx
               cgf(nseq0)=gfx
               cgx(nseq0)=gxx
               cgt(nseq0)=gnx+ggx+gfx+gxx
               ndfn(nseq0)=nu
               ndff(nseq0)=mu
               ndfx(nseq0)=lu
               dbar(nseq0)=dx

               !--compute infinite dilution cross sections
               call gnrx(gnx,gfx,ggx,nu,mu,lu,gs,gxx,1)
               call gnrx(gnx,gfx,ggx,nu,mu,lu,gc,gxx,2)
               call gnrx(gnx,gfx,ggx,nu,mu,lu,gf,gxx,3)
               temp=abn*pi*ab*gj*gnx/(2*dx)
               sigi(2)=sigi(2)+temp*(gs*gnx-2*sin(ps)**2)
               sigi(3)=sigi(3)+temp*gf*gfx
               sigi(4)=sigi(4)+temp*gc*ggx
            enddo
         enddo
      endif
   enddo

   return
   end subroutine unresx

   subroutine unfac2(l,rho,rhoc,amun,vl,ps)
   !-------------------------------------------------------------------
   ! Calculates the penetrability factor (vl) and phase shift (ps)
   !-------------------------------------------------------------------
   ! externals
   integer::l
   real(kr)::rho,rhoc,amun,vl,ps
   ! internals
   real(kr)::r2,r4

   if (l.eq.0) then
      vl=amun
      ps=rhoc
   else if (l.eq.1) then
      r2=rho*rho
      vl=amun*r2/(1+r2)
      ps=rhoc-atan(rhoc)
   else if (l.eq.2) then
      r2=rho*rho
      r4=r2*r2
      vl=amun*r4/(9+3*r2+r4)
      ps=rhoc-atan(3*rhoc/(3-rhoc*rhoc))
   endif
   return
   end subroutine unfac2

   subroutine intrf2(e,gfx,iest,ifst,ne,int,a)
   !-------------------------------------------------------------------
   ! Interpolates fission widths for unresolved representation 1.
   !-------------------------------------------------------------------
   use endf ! provides terp1
   ! externals
   integer::iest,ifst,ne,int
   real(kr)::e,gfx,a(*)
   ! internals
   integer::i,idone,i1,i2

   i=1
   idone=0
   i1=0
   i2=0
   do while (i.lt.ne.and.idone.eq.0)
      i=i+1
      i1=iest+i-1
      i2=ifst+i-1
      if (e.lt.a(i1)) idone=1
   enddo
   call terp1(a(i1-1),a(i2-1),a(i1),a(i2),e,gfx,int)
   return
   end subroutine intrf2

   subroutine intr2(e,dx,gxx,gnox,ggx,gfx,inow,ne,int,a)
   !-------------------------------------------------------------------
   ! Interpolates energy dependent unresolved resonance parameters.
   !-------------------------------------------------------------------
   use endf ! provides terp
   ! externals
   integer::inow,ne,int
   real(kr)::e,dx,gxx,gnox,ggx,gfx,a(*)
   ! internals
   integer::i1,i,idone,i2
   real(kr),parameter::small=1.e-8_kr

   i1=inow
   i2=0
   i=1
   idone=0
   do while (i.lt.ne.and.idone.eq.0)
      i=i+1
      i2=i1+6
      if (e.lt.a(i2)) then
         idone=1
      else
         if (i.lt.ne) i1=i2
      endif
   enddo
   call terp1(a(i1),a(i1+1),a(i2),a(i2+1),e,dx,int)
   call terp1(a(i1),a(i1+2),a(i2),a(i2+2),e,gxx,int)
   call terp1(a(i1),a(i1+3),a(i2),a(i2+3),e,gnox,int)
   call terp1(a(i1),a(i1+4),a(i2),a(i2+4),e,ggx,int)
   call terp1(a(i1),a(i1+5),a(i2),a(i2+5),e,gfx,int)
   if (gxx.lt.small) gxx=0
   if (gfx.lt.small) gfx=0
   return
   end subroutine intr2

   subroutine gnrx(galpha,gbeta,gamma,mu,nu,lamda,s,df,id)
   !-------------------------------------------------------------------
   ! Calculates fluctuation integrals for unresolved resonances
   ! using MC2-2 quadrature scheme.
   !-------------------------------------------------------------------
   ! externals
   integer::mu,nu,lamda,id
   real(kr)::galpha,gbeta,gamma,s,df
   ! internals
   integer::j,k,l
   real(kr)::xj,wj,xk,wk,xl,wl
   real(kr),dimension(10,4),parameter::qw=reshape((/&
     1.1120413e-1_kr,2.3546798e-1_kr,2.8440987e-1_kr,2.2419127e-1_kr,&
     .10967668e0_kr,.030493789e0_kr,.0042930874e0_kr,2.5827047e-4_kr,&
     4.9031965e-6_kr,1.4079206e-8_kr,.033773418e0_kr,.079932171e0_kr,&
     .12835937e0_kr,.17652616e0_kr,.21347043e0_kr,.21154965e0_kr,&
     .13365186e0_kr,.022630659e0_kr,1.6313638e-5_kr,2.745383e-31_kr,&
     3.3376214e-4_kr,.018506108e0_kr,.12309946e0_kr,.29918923e0_kr,&
     .33431475e0_kr,.17766657e0_kr,.042695894e0_kr,4.0760575e-3_kr,&
     1.1766115e-4_kr,5.0989546e-7_kr,1.7623788e-3_kr,.021517749e0_kr,&
     .080979849e0_kr,.18797998e0_kr,.30156335e0_kr,.29616091e0_kr,&
     .10775649e0_kr,2.5171914e-3_kr,8.9630388e-10_kr,0.e0_kr/),&
     (/10,4/))
   real(kr),dimension(10,4),parameter::qp=reshape((/&
     3.0013465e-3_kr,7.8592886e-2_kr,4.3282415e-1_kr,1.3345267e0_kr,&
     3.0481846e0_kr,5.8263198e0_kr,9.9452656e0_kr,1.5782128e+1_kr,&
     23.996824e0_kr,36.216208e0_kr,1.3219203e-2_kr,7.2349624e-2_kr,&
     .19089473e0_kr,.39528842e0_kr,.74083443e0_kr,1.3498293e0_kr,&
     2.5297983e0_kr,5.2384894e0_kr,13.821772e0_kr,75.647525e0_kr,&
     1.0004488e-3_kr,.026197629e0_kr,.14427472e0_kr,.44484223e0_kr,&
     1.0160615e0_kr,1.9421066e0_kr,3.3150885e0_kr,5.2607092e0_kr,&
     7.9989414e0_kr,12.072069e0_kr,.013219203e0_kr,.072349624e0_kr,&
     .19089473e0_kr,.39528842e0_kr,.74083443e0_kr,1.3498293e0_kr,&
     2.5297983e0_kr,5.2384894e0_kr,13.821772e0_kr,75.647525e0_kr/),&
     (/10,4/))
   real(kr),parameter::zero=0

   s=0
   if (galpha.le.zero) return
   if (gamma.le.zero) return
   if (gbeta.lt.zero) return
   if (gbeta.le.zero.and.df.lt.zero) return
   if (gbeta.le.zero) then
      if (df.le.zero) then
         do j=1,10
            xj=qp(j,mu)
            wj=qw(j,mu)
            if (id.eq.1) then
               s=s+wj*xj*xj/(galpha*xj+gamma)
            else if (id.eq.2) then
               s=s+wj*xj/(galpha*xj+gamma)
            endif
         enddo
      else
         do j=1,10
            xj=qp(j,mu)
            wj=qw(j,mu)
            do k=1,10
               xk=qp(k,lamda)
               wk=qw(k,lamda)
               if (id.eq.1) then
                  s=s+wj*wk*xj*xj/(galpha*xj+gamma+df*xk)
               else if (id.eq.2) then
                  s=s+wj*wk*xj/(galpha*xj+gamma+df*xk)
               endif
            enddo
         enddo
      endif
   else
      if (df.le.zero) then
         if (df.ge.zero) then
            do j=1,10
               xj=qp(j,mu)
               wj=qw(j,mu)
               do k=1,10
                  xk=qp(k,nu)
                  wk=qw(k,nu)
                  if (id.eq.1) then
                     s=s+wj*wk*xj*xj/(galpha*xj+gbeta*xk+gamma)
                  else if (id.eq.2) then
                     s=s+wj*wk*xj/(galpha*xj+gbeta*xk+gamma)
                  else if (id.eq.3) then
                     s=s+wj*wk*xj*xk/(galpha*xj+gbeta*xk+gamma)
                  endif
               enddo
            enddo
         endif
      else
         do j=1,10
            xj=qp(j,mu)
            wj=qw(j,mu)
            do k=1,10
               xk=qp(k,nu)
               wk=qw(k,nu)
               do l=1,10
                  xl=qp(l,lamda)
                  wl=qw(l,lamda)
                  if (id.eq.1) then
                     s=s+wj*wk*wl*xj*xj&
                       /(galpha*xj+gbeta*xk+gamma+df*xl)
                  else if (id.eq.2) then
                     s=s+wj*wk*wl*xj&
                       /(galpha*xj+gbeta*xk+gamma+df*xl)
                  else if (id.eq.3) then
                     s=s+wj*wk*wl*xj*xk&
                       /(galpha*xj+gbeta*xk+gamma+df*xl)
                  endif
               enddo
            enddo
         enddo
      endif
   endif
   end subroutine gnrx

   subroutine ladr2(nseqz,elow,ehigh,nr,er,gnr,gfr,ggr,gxr,gt)
   !-------------------------------------------------------------------
   ! Generate ladders of statistically chosen resonances from
   ! average resonance data from ENDF data.
   !-------------------------------------------------------------------
   use physics ! provides pi
   use util ! provides error
   ! externals
   integer::nseqz,nr
   real(kr)::elow,ehigh,gg,gf,gx,gn,xr,dn
   real(kr)::er(*),gnr(*),gfr(*),ggr(*),gxr(*),gt(*)
   ! internals
   integer::nr0,ir,idone,ndf,i,n
   real(kr)::dcon
   real(kr),dimension(20,4),parameter::chisq=reshape((/&
     1.31003e-3_kr,9.19501e-3_kr,.0250905e0_kr,.049254e0_kr,&
     .0820892e0_kr,.124169e0_kr,.176268e0_kr,.239417e0_kr,&
     .314977e0_kr,.404749e0_kr,.511145e0_kr,.637461e0_kr,&
     .788315e0_kr,.970419e0_kr,1.194e0_kr,1.47573e0_kr,&
     1.84547e0_kr,2.36522e0_kr,3.20371e0_kr,5.58201e0_kr,&
     .0508548e0_kr,.156167e0_kr,.267335e0_kr,.38505e0_kr,&
     .510131e0_kr,.643564e0_kr,.786543e0_kr,.940541e0_kr,&
     1.1074e0_kr,1.28947e0_kr,1.48981e0_kr,1.71249e0_kr,&
     1.96314e0_kr,2.24984e0_kr,2.58473e0_kr,2.98744e0_kr,&
     3.49278e0_kr,4.17238e0_kr,5.21888e0_kr,7.99146e0_kr,&
     .206832e0_kr,.470719e0_kr,.691933e0_kr,.901674e0_kr,&
     1.10868e0_kr,1.31765e0_kr,1.53193e0_kr,1.75444e0_kr,&
     1.98812e0_kr,2.23621e0_kr,2.50257e0_kr,2.79213e0_kr,&
     3.11143e0_kr,3.46967e0_kr,3.88053e0_kr,4.36586e0_kr,&
     4.96417e0_kr,5.75423e0_kr,6.94646e0_kr,10.0048e0_kr,&
     .459462e0_kr,.893735e0_kr,1.21753e0_kr,1.50872e0_kr,&
     1.78605e0_kr,2.05854e0_kr,2.33194e0_kr,2.61069e0_kr,&
     2.89878e0_kr,3.20032e0_kr,3.51995e0_kr,3.86331e0_kr,&
     4.23776e0_kr,4.65345e0_kr,5.12533e0_kr,5.67712e0_kr,&
     6.35044e0_kr,7.22996e0_kr,8.541e0_kr,11.8359e0_kr/),(/20,4/))
   real(kr),parameter::start=19.9999e0_kr
   real(kr),parameter::zero=0.e0_kr
     real(kr)::rnsum
     rnsum=0

   !--select resonance parameters
   !--until resonance energy is above ehi
   dcon=dbar(nseqz)*sqrt(4/pi)
   nr0=0
   ir=nr0
   idone=0
   do while (idone.eq.0)
      ir=ir+1
      if (ir.gt.nermax) call error('ladr2',&
        'too many resonances in ladder',' ')
      if (ir.ne.nr0+1) then
         ! select next resonance location from wigner distribution
         rr=rann(kk)
         er(ir)=er(ir-1)+dcon*sqrt(-log(rr))
      else
         ! choose starting point from uniform distribution
         rr=rann(kk)
         er(ir)=elow+dcon*rr
      endif

      !--select level widths
      i=nseqz
      gg=cgg(i)
      gf=cgf(i)
      gx=cgx(i)
      gn=cgn(i)
      ! choose neutron width for ndn degrees of freedom
      ndf=ndfn(i)
      if (ndf.gt.0) gn=gn/ndf
      n=int(1+start*rann(kk))
         rnsum=rnsum+n
      xr=chisq(n,ndf)
      gn=gn*xr
      ! choose fission width for ndf degrees of freedom
      if (cgf(i).ne.zero.and.ndff(i).ne.0) then
         ndf=ndff(i)
         if (ndf.gt.0) gf=gf/ndf
         n=int(1+start*rann(kk))
         xr=chisq(n,ndf)
         gf=gf*xr
      endif
      ! choose competitive width for ndx degrees of freedom
      if (cgx(i).ne.zero.and.ndfx(i).ne.0) then
         ndf=ndfx(i)
         if (ndf.gt.0) gx=gx/ndf
         dn=20
         dn=dn-dn/10000
         n=int(1+dn*rann(kk))
         xr=chisq(n,ndf)
         gx=gx*xr
      endif
      gt(ir)=gn+gf+gg+gx
      gnr(ir)=gn/gt(ir)
      gfr(ir)=gf/gt(ir)
      ggr(ir)=gg/gt(ir)
      gxr(ir)=gx/gt(ir)
      if (er(ir).gt.ehigh) idone=1
   enddo
   nr=ir-1
   return
   end subroutine ladr2

   subroutine unrest(bkg,sig0,sigf,tabl,tval,er,gnr,gfr,ggr,gxr,gt,&
     es,xs,cap,fis,els,temp,nsig0,ntemp,nbin,nsamp)
   !-------------------------------------------------------------------
   ! Probability tables and Bondarenko moments.
   !-------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use physics ! provides pi
   use util ! provides timer
   ! externals
   integer::nsig0,ntemp,nbin,nsamp
   real(kr)::bkg(4),sig0(nsig0),tval(nbin,ntemp)
   real(kr)::sigf(5,nsig0,ntemp),tabl(nbin,5,ntemp)
   real(kr)::er(*),gnr(*),gfr(*),ggr(*),gxr(*),gt(*)
   real(kr)::es(nsamp),xs(nsamp)
   real(kr)::cap(ntemp,nsamp),fis(ntemp,nsamp),els(ntemp,nsamp)
   real(kr)::temp(ntemp)
   ! internals
   integer::navoid,i,nres,nrest,iladr,ne,ie,itemp,k,nr,jr
   integer::i7,i0,it,i6,i1,i5,i2,is,i4,i3,jj,j,n,ii
   integer::ixx,l,mfl,nebin,ibin,izeroprob,ibadxs
   real(kr)::rpi,binmin,elow,dmin,erange,ehigh,emin,emax,espan
   real(kr)::dbart,sigx,ctx,chek1,chek2,chekn,delr,elo,ehi
   real(kr)::y,yy,szy,cc2,cs2,ccg,ccf,test,x,a1,rew,aimw
   real(kr)::a2,a3,temp1,temp2,a4,a5,f1,f2,a6,e1,e2,e3
   real(kr)::tempor,q,q2,hq,hq2,ax,aki,p,p2,hp,hp2,pq,tav
   real(kr)::tvar,eav,evar,fav,fvar,cav,cvar,totf,elsf,capf,fisf
   real(kr)::efact,cfact,ffact
   real(kr)::tot,tem,argt,arge,argf,argc,tnorm,enorm
   real(kr)::fnorm,cnorm,denom,den,ttt
   character(3),dimension(4),parameter::nmr=(/&
     'tot','els','fis','cap'/)
   character(60)::strng
   real(kr),parameter::c1=0.5641895835e0_kr
   real(kr),parameter::c2=0.2752551e0_kr
   real(kr),parameter::c3=2.724745e0_kr
   real(kr),parameter::c4=0.5124242e0_kr
   real(kr),parameter::c5=0.05176536e0_kr
   real(kr),parameter::half=0.5e0_kr
   real(kr),parameter::ten=10.0e0_kr
   real(kr),parameter::d1=0.4613135e0_kr
   real(kr),parameter::d2=0.1901635e0_kr
   real(kr),parameter::d3=0.09999216e0_kr
   real(kr),parameter::d4=1.7844927e0_kr
   real(kr),parameter::d5=0.002883894e0_kr
   real(kr),parameter::d6=5.5253437e0_kr
   real(kr),parameter::fifty=50
   real(kr),parameter::con1=31.83e0_kr
   real(kr),parameter::con2=20.e0_kr
   real(kr),parameter::break1=.5e0_kr
   real(kr),parameter::break2=3.9e0_kr
   real(kr),parameter::tref=300.e0_kr
   real(kr),parameter::big=1.e6_kr
   real(kr),parameter::zero=0
   save navoid,binmin,elow
   rpi=sqrt(pi)
   izeroprob=0
   ibadxs=0
   ! data for optional bondarenko plots
   if (ipl.lt.200) ipl=ipl+1
   epl(ipl)=e

   !--initialize
   if (init.eq.0) then
      init=1
      navoid=300
      binmin=spot/10
      elow=10
      if (allocated(bval)) then
         deallocate(bval)
         deallocate(tmin)
         deallocate(tmax)
         deallocate(tsum)
      endif
      allocate(bval(8,nsig0,ntemp))
      allocate(tmin(ntemp))
      allocate(tmax(ntemp))
      allocate(tsum(ntemp))
   endif
   dmin=100000
   do i=1,nseq0
      if (dbar(i).lt.dmin) dmin=dbar(i)
   enddo
   erange=9*nermax*dmin/10
   nres=int(erange*dbarin)
   nrest=nres-2*navoid
   ehigh=elow+erange
   emin=elow+navoid/dbarin
   emax=ehigh-navoid/dbarin
   espan=emax-emin
   if ((nres.le.zero).or.(emax.lt.emin)) then
      call error('unrest','bad value for nres or emin>emax, increase dmin','')
   endif
   dbart=1/dbarin
   sigx=bkg(1)-bkg(2)-bkg(3)-bkg(4)
   do itemp=1,ntemp
      do i=nsig0,1,-1
         do j=1,8
            bval(j,i,itemp)=0
         enddo
      enddo
   enddo

   if (iprint.gt.0) then
      write(nsyso,'(/,1p,''e='',e11.4,&
        &2x,''spot='',e11.4,2x,''dbar='',e11.4,2x,''sigx='',e11.4/&
        &13x,''total'',5x,''elastic'',5x,''fission'',5x,&
        &''capture'')') e,spot,dbart,sigx
   endif

   !--loop over ladders
   !--use first pass to set prob table limits
   do 150 iladr=1,nladr
   ne=nsamp

   !--choose random energy grid
   do ie=1,ne
      rr=rann(kk)
      es(ie)=emin+espan*rr
      do itemp=1,ntemp
         cap(itemp,ie)=0
         fis(itemp,ie)=0
         els(itemp,ie)=spot
      enddo
   enddo
   call fsort(es,xs,ne,1)

   !--loop over sequences
   do 170 k=1,nseq0
   call ladr2(k,elow,ehigh,nr,er,gnr,gfr,ggr,gxr,gt)

   !--loop over temperatures
   do 180 itemp=1,ntemp
   ctx=cth(k)*sqrt(tref/temp(itemp))

   !--accumulate contributions to cross sections
   do 210 jr=1,nr
   chek1=con1*gt(jr)
   chek2=con2/ctx
   chekn=chek1
   if (chek2.gt.chek1) chekn=chek2
   delr=chek1+chekn
   elo=er(jr)-delr
   ehi=er(jr)+delr
   call fsrch(ehi,es,ne,i7,mfl)
   call fsrch(elo,es,i7,i0,mfl)
   if (i0.eq.i7) go to 210
   it=i7-i0+1
   y=ctx*gt(jr)/2
   yy=y*y
   szy=csz(k)*gnr(jr)*rpi*y
   cc2=szy*(cc2p(k)-1+gnr(jr))
   cs2=szy*cs2p(k)
   ccg=szy*ggr(jr)
   ccf=szy*gfr(jr)
   do ie=i0,i7
      xs(ie)=ctx*(es(ie)-er(jr))
   enddo

   !--accumulate cross sections using an
   !--appropriate approximation in each range

    !  asymptotic term for x.gt.100 or y.gt.100

   i6=i0-1
   test=100
   if (y.gt.test) go to 340
   call fsrch(-test,xs(i0),it,i1,mfl)
   if (i1.eq.1.and.mfl.eq.2) i1=0
   i1=i1+i0-1
   if (i1.lt.i0) i1=i0-1
   if (i1.le.i0) go to 240
   do ie=i0,i1
      x=xs(ie)
      test=x*x+yy
      a1=c1/test
      rew=y*a1
      aimw=x*a1
      cap(itemp,ie)=cap(itemp,ie)+ccg*rew
      fis(itemp,ie)=fis(itemp,ie)+ccf*rew
      els(itemp,ie)=els(itemp,ie)+cc2*rew+cs2*aimw
   enddo
   if (i1.ge.i7) go to 370

     ! rational approximation for x.gt.6 or y.gt.6 but x and y .le.100

  240 continue
   i5=i1
   test=6
   if (y.gt.test) go to 330
   test=-6
   call fsrch(test,xs(i0),it,i2,mfl)
   if (i2.eq.1.and.mfl.eq.2) i2=0
   i2=i2+i0-1
   if (i2.lt.i1) i2=i1
   if (i2.le.i1) go to 250
   is=i1+1
   do ie=is,i2
      x=xs(ie)
      a1=x*x-yy
      a2=2*x*y
      a3=a2*a2
      temp1=a2*x
      temp2=a2*y
      a4=a1-c2
      a5=a1-c3
      f1=c4/(a4*a4+a3)
      f2=c5/(a5*a5+a3)
      rew=f1*(temp1-a4*y)+f2*(temp1-a5*y)
      aimw=f1*(a4*x+temp2)+f2*(a5*x+temp2)
      cap(itemp,ie)=cap(itemp,ie)+ccg*rew
      fis(itemp,ie)=fis(itemp,ie)+ccf*rew
      els(itemp,ie)=els(itemp,ie)+cc2*rew+cs2*aimw
   enddo
   if (i2.ge.i7) go to 370

     ! rational approximation for x.gt.3.9 or y.gt.3.0 but x and y .le.6

  250 continue
   i4=i2
   test=3
   if (y.gt.test) go to 320
   test=-break2
   call fsrch(test,xs(i0),it,i3,mfl)
   if (i3.eq.1.and.mfl.eq.2) i3=0
   i3=i3+i0-1
   if (i3.lt.i2) i3=i2
   if (i3.le.i2) go to 260
   is=i2+1
   do ie=is,i3
      x=xs(ie)
      a1=x*x-yy
      a2=2*x*y
      a3=a2*a2
      temp1=a2*x
      temp2=a2*y
      a4=a1-d2
      a5=a1-d4
      a6=a1-d6
      e1=d1/(a4*a4+a3)
      e2=d3/(a5*a5+a3)
      e3=d5/(a6*a6+a3)
      rew=e1*(temp1-a4*y)+e2*(temp1-a5*y)+e3*(temp1-a6*y)
      aimw=e1*(a4*x+temp2)+e2*(a5*x+temp2)+e3*(a6*x+temp2)
      cap(itemp,ie)=cap(itemp,ie)+ccg*rew
      fis(itemp,ie)=fis(itemp,ie)+ccf*rew
      els(itemp,ie)=els(itemp,ie)+cc2*rew+cs2*aimw
   enddo
   if (i3.ge.i7) go to 370

     ! six point table interpolation for 0.le.abs(x).le.3.9

  260 continue
   test=break2
   call fsrch(test,xs(i0),it,i4,mfl)
   if (i4.eq.1.and.mfl.eq.2) i4=0
   i4=i4+i0-1
   if (i4.lt.i3) i4=i3
   if (i4.le.i3) go to 320
   is=i3+1

     ! y.ge.0.5

   if (y.lt.break1) go to 310
   tempor=ten*y
   jj=int(tempor)
   j=jj-3
   n=j-1
   q=tempor-jj
   q2=q*q
   hq=half*q
   hq2=half*q2
   a1=hq2-hq
   do ie=is,i4
      x=xs(ie)
      ax=abs(x)
      aki=1
      if (x.lt.zero) aki=-1
      tempor=ten*ax
      ii=int(tempor)
      i=ii+2
      p=tempor-ii
      p2=p*p
      hp=half*p
      hp2=half*p2
      a2=hp2-hp
      pq=p*q
      a3=1+pq-p2-q2
      a4=hp2-pq+hp
      a5=hq2-pq+hq
      rew=a1*tr(i,n)+a2*tr(i-1,j)+a3*tr(i,j)+a4*tr(i+1,j)&
        +a5*tr(i,j+1)+pq*tr(i+1,j+1)
      aimw=a1*ti(i,n)+a2*ti(i-1,j)+a3*ti(i,j)+a4*ti(i+1,j)&
        +a5*ti(i,j+1)+pq*ti(i+1,j+1)
      aimw=aimw*aki
      cap(itemp,ie)=cap(itemp,ie)+ccg*rew
      fis(itemp,ie)=fis(itemp,ie)+ccf*rew
      els(itemp,ie)=els(itemp,ie)+cc2*rew+cs2*aimw
   enddo
   if (i4.ge.i7) go to 370
   go to 320

     ! y.lt.0.5

  310 continue
   tempor=fifty*y
   jj=int(tempor)
   j=jj+2
   n=j-1
   q=tempor-jj
   q2=q*q
   hq=half*q
   hq2=half*q2
   a1=hq2-hq
   do ie=is,i4
      x=xs(ie)
      ax=abs(x)
      aki=1
      if (x.lt.zero) aki=-1
      tempor=ten*ax
      ii=int(tempor)
      i=ii+2
      p=tempor-ii
      p2=p*p
      hp=half*p
      hp2=half*p2
      a2=hp2-hp
      pq=p*q
      a3=1+pq-p2-q2
      a4=hp2-pq+hp
      a5=hq2-pq+hq
      rew=a1*trs(i,n)+a2*trs(i-1,j)+a3*trs(i,j)+a4*trs(i+1,j)+&
        a5*trs(i,j+1)+pq*trs(i+1,j+1)
      aimw=a1*tis(i,n)+a2*tis(i-1,j)+a3*tis(i,j)+a4*tis(i+1,j)+&
        a5*tis(i,j+1)+pq*tis(i+1,j+1)
      aimw=aimw*aki
      cap(itemp,ie)=cap(itemp,ie)+ccg*rew
      fis(itemp,ie)=fis(itemp,ie)+ccf*rew
      els(itemp,ie)=els(itemp,ie)+cc2*rew+cs2*aimw
   enddo
   if (i4.ge.i7) go to 370

   ! rational approximation for x.gt.3.9 or y.gt.3.0 but x and y .le.6

  320 continue
   test=6
   call fsrch(test,xs(i0),it,i5,mfl)
   if (i5.eq.1.and.mfl.eq.2) i5=0
   i5=i5+i0-1
   if (i5.lt.i4) i5=i4
   if (i5.le.i4) go to 330
   is=i4+1
   do ie=is,i5
      x=xs(ie)
      a1=x*x-yy
      a2=2*x*y
      a3=a2*a2
      temp1=a2*x
      temp2=a2*y
      a4=a1-d2
      a5=a1-d4
      a6=a1-d6
      e1=d1/(a4*a4+a3)
      e2=d3/(a5*a5+a3)
      e3=d5/(a6*a6+a3)
      rew=e1*(temp1-a4*y)+e2*(temp1-a5*y)+e3*(temp1-a6*y)
      aimw=e1*(a4*x+temp2)+e2*(a5*x+temp2)+e3*(a6*x+temp2)
      cap(itemp,ie)=cap(itemp,ie)+ccg*rew
      fis(itemp,ie)=fis(itemp,ie)+ccf*rew
      els(itemp,ie)=els(itemp,ie)+cc2*rew+cs2*aimw
   enddo
   if (i5.ge.i7) go to 370

   ! rational approximation for x.gt.6 or y.gt.6 but x and y .le.100

  330 continue
   test=100
   call fsrch(test,xs(i0),it,i6,mfl)
   if (i6.eq.1.and.mfl.eq.2) i6=0
   i6=i6+i0-1
   if (i6.lt.i5) i6=i5
   if (i6.le.i5) go to 340
   is=i5+1
   do ie=is,i6
      x=xs(ie)
      a1=x*x-yy
      a2=2*x*y
      a3=a2*a2
      temp1=a2*x
      temp2=a2*y
      a4=a1-c2
      a5=a1-c3
      f1=c4/(a4*a4+a3)
      f2=c5/(a5*a5+a3)
      rew=f1*(temp1-a4*y)+f2*(temp1-a5*y)
      aimw=f1*(a4*x+temp2)+f2*(a5*x+temp2)
      cap(itemp,ie)=cap(itemp,ie)+ccg*rew
      fis(itemp,ie)=fis(itemp,ie)+ccf*rew
      els(itemp,ie)=els(itemp,ie)+cc2*rew+cs2*aimw
   enddo
   if (i6.ge.i7) go to 370

    ! asymptotic term for x.gt.100 or y.gt.100

  340 continue
   is=i6+1
   do ie=is,i7
      x=xs(ie)
      test=x*x+yy
      a1=c1/test
      rew=y*a1
      aimw=x*a1
      cap(itemp,ie)=cap(itemp,ie)+ccg*rew
      fis(itemp,ie)=fis(itemp,ie)+ccf*rew
      els(itemp,ie)=els(itemp,ie)+cc2*rew+cs2*aimw
   enddo
  370 continue

   !--close loops over resonances and sequences
  210 continue
  180 continue
  170 continue

   !--eliminate negative elastic cross sections; reset to 1 microbarn
   do itemp=1,ntemp
      do ie=1,ne
         if (els(itemp,ie).lt.-bkg(2)) els(itemp,ie)=-bkg(2)+1/big
      enddo
   enddo

   !--compute infinitely-dilute cross sections
   if (iladr.eq.1) then
      tav=0
      tvar=0
      eav=0
      evar=0
      fav=0
      fvar=0
      cav=0
      cvar=0
   endif
   totf=0
   elsf=0
   capf=0
   fisf=0
   do ie=1,ne
      totf=totf+els(1,ie)+fis(1,ie)+cap(1,ie)+bkg(1)
      elsf=elsf+els(1,ie)+bkg(2)
      capf=capf+cap(1,ie)+bkg(4)
      fisf=fisf+fis(1,ie)+bkg(3)
   enddo
   totf=totf/ne
   elsf=elsf/ne
   capf=capf/ne
   fisf=fisf/ne
   tav=tav+totf
   tvar=tvar+totf*totf
   eav=eav+elsf
   evar=evar+elsf*elsf
   fav=fav+fisf
   fvar=fvar+fisf*fisf
   cav=cav+capf
   cvar=cvar+capf*capf
   if (iprint.gt.0) write(nsyso,'(i6,1p,4e12.4)')&
     iladr,totf,elsf,fisf,capf

   !--renormalize reaction cross sections
   !--compute total cross section
   if (nmode.eq.1) then
      efact=(sigi(2)+bkg(2))/(elsf-spot)
      cfact=(sigi(4)+bkg(4))/capf
      ffact=1
      if (fisf.ne.zero) ffact=(sigi(3)+bkg(3))/fisf
      do ie=1,ne
         do itemp=1,ntemp
            els(itemp,ie)=spot+(els(itemp,ie)-spot)*efact
            cap(itemp,ie)=cap(itemp,ie)*cfact
            fis(itemp,ie)=fis(itemp,ie)*ffact
         enddo
      enddo
   endif

   !--loop over temperatures
   !--using a different total cross section bin structure for each
   do 140 itemp=1,ntemp

   !--choose probability table bounds
   !--and zero accumulators
   if (iladr.eq.1) then
      do ie=1,ne
         es(ie)=els(itemp,ie)+fis(itemp,ie)+cap(itemp,ie)+bkg(1)
      enddo
      call fsort(es,xs,ne,1)
      tmin(itemp)=es(1)
      tmax(itemp)=es(ne)
      nebin=int(nsamp/(nbin-10+1.76))
      ibin=nebin/200
      if (ibin.le.0) then
         if (mflg1.eq.0) then
            mflg1=1
            write(strng,'(''reset ibin=1, consider larger nsamp'',&
                         &'' or smaller nbin'')')
            call mess('purr',strng,'')
         endif
         ibin=1
      endif
      do i=1,nbin-1
         if (ibin.gt.nsamp) then
            if (mflg2.eq.0) then
               mflg2=1
               write(strng,'(''reset ibin=nsamp,'',&
                            &'' consider smaller nbin'')')
               call mess('purr',strng,'')
            endif
            ibin=nsamp
         endif
         tval(i,itemp)=es(ibin)
         if (i.gt.1) then
            if (tval(i,itemp).le.tval(i-1,itemp))&
              tval(i,itemp)=tval(i-1,itemp)+tval(i-1,itemp)/20
         endif
         if (i.eq.1) ibin=ibin+nebin/40
         if (i.eq.2) ibin=ibin+nebin/10
         if (i.eq.3) ibin=ibin+nebin/4
         if (i.eq.4) ibin=ibin+nebin/2
         if (i.gt.4.and.i.lt.nbin-5) ibin=ibin+nebin
         if (i.eq.nbin-5) ibin=ibin+nebin/2
         if (i.eq.nbin-4) ibin=ibin+nebin/4
         if (i.eq.nbin-3) ibin=ibin+nebin/10
         if (i.eq.nbin-2) ibin=ibin+nebin/40
         if (i.eq.nbin-1) ibin=ibin+nebin/200
      enddo
      tval(nbin,itemp)=big
      do i=1,nbin
         do j=1,5
            tabl(i,j,itemp)=0
         enddo
      enddo
      tsum(itemp)=0
   endif

   !--accumulate cross sections into tables
   do ie=1,ne
      tot=els(itemp,ie)+fis(itemp,ie)+cap(itemp,ie)+bkg(1)
      if (tot.lt.tmin(itemp)) tmin(itemp)=tot
      if (tot.gt.tmax(itemp)) tmax(itemp)=tot
      call fsrch(tot,tval(1,itemp),nbin,ii,mfl)
      if (mfl.ne.2.and.ii.lt.nbin) ii=ii+1
      if (ii.lt.1) ii=1
      tsum(itemp)=tsum(itemp)+1
      tabl(ii,1,itemp)=tabl(ii,1,itemp)+1
      tabl(ii,2,itemp)=tabl(ii,2,itemp)+tot
      tabl(ii,3,itemp)=tabl(ii,3,itemp)+els(itemp,ie)+bkg(2)
      tabl(ii,4,itemp)=tabl(ii,4,itemp)+fis(itemp,ie)+bkg(3)
      tabl(ii,5,itemp)=tabl(ii,5,itemp)+cap(itemp,ie)+bkg(4)
      do i=1,nsig0
         tem=sig0(i)/(sig0(i)+tot)
         bval(1,i,itemp)=bval(1,i,itemp)+tot*tem
         bval(2,i,itemp)=bval(2,i,itemp)+(els(itemp,ie)+bkg(2))*tem
         bval(3,i,itemp)=bval(3,i,itemp)+(fis(itemp,ie)+bkg(3))*tem
         bval(4,i,itemp)=bval(4,i,itemp)+(cap(itemp,ie)+bkg(4))*tem
         bval(5,i,itemp)=bval(5,i,itemp)+tot*tem*tem
         bval(6,i,itemp)=bval(6,i,itemp)+tem
         bval(7,i,itemp)=bval(7,i,itemp)+tem*tem
      enddo
   enddo

   !--close loop over temperatures
  140 continue

   !--close loop over ladders
  150 continue

   !--write overall average cross sections
   tav=tav/nladr
   argt=tvar/nladr-tav*tav
   if (argt.lt.zero) argt=0
   eav=eav/nladr
   arge=evar/nladr-eav*eav
   if (arge.lt.zero) arge=0
   fav=fav/nladr
   argf=fvar/nladr-fav*fav
   if (argf.lt.zero) argf=0
   cav=cav/nladr
   argc=cvar/nladr-cav*cav
   if (argc.lt.zero) argc=0
   tvar=100*sqrt(argt)/tav
   evar=100*sqrt(arge)/eav
   if (fav.ne.zero) fvar=100*sqrt(argf)/fav
   cvar=100*sqrt(argc)/cav
   sigi(1)=sigi(2)+sigi(3)+sigi(4)+spot+bkg(1)
   sigi(2)=sigi(2)+spot+bkg(2)
   sigi(3)=sigi(3)+bkg(3)
   sigi(4)=sigi(4)+bkg(4)
   if (iprint.gt.0) write(nsyso,'(&
     &''  bkgd'',1p,4e12.4/''  infd'',1p,4e12.4/&
     &''  aver'',1p,4e12.4/''  pcsd'',0p,4f12.2/&
     &''  nres'',i12)')&
     (bkg(i),i=1,4),(sigi(i),i=1,4),&
     tav,eav,fav,cav,tvar,evar,fvar,cvar,nrest

   !--compute and write bondarenko table
   if (iprint.gt.0) write(nsyso,'(/&
     &'' bondarenko cross sections by direct sampling''/&
     &9x,''temp'',6x,''sig0'',4x,''p0 total'',5x,&
     &''elastic'',5x,''fission'',5x,''capture'',4x,&
     &''p1 total'')')
   do i=1,nsig0
      do itemp=1,ntemp
         sigf(1,i,itemp)=bval(1,i,itemp)/bval(6,i,itemp)
         sigf(2,i,itemp)=bval(2,i,itemp)/bval(6,i,itemp)
         sigf(3,i,itemp)=bval(3,i,itemp)/bval(6,i,itemp)
         sigf(4,i,itemp)=bval(4,i,itemp)/bval(6,i,itemp)
         sigf(5,i,itemp)=bval(5,i,itemp)/bval(7,i,itemp)
      enddo
   enddo
   if (iprint.gt.0) then
      do itemp=1,ntemp
         do i=1,nsig0
            write(nsyso,'(3x,1p,2e10.3,5e12.4)')&
              temp(itemp),sig0(i),(sigf(j,i,itemp),j=1,5)
         enddo
      enddo
   endif

   !--normalize and write probability table
   if (iprint.gt.0) write(nsyso,'(/'' probability table'')')
   tnorm=sigi(1)-sigf(1,1,1)
   enorm=sigi(2)-sigf(2,1,1)
   fnorm=sigi(3)-sigf(3,1,1)
   cnorm=sigi(4)-sigf(4,1,1)
   tnorm=0
   enorm=0
   fnorm=0
   cnorm=0
   do i=1,nbin
      do itemp=1,ntemp
         tmin(itemp)=tmin(itemp)+tnorm
         tval(nbin,itemp)=tmax(itemp)
         denom=tabl(i,1,itemp)
         if (denom.eq.zero) denom=1
         tabl(i,1,itemp)=tabl(i,1,itemp)/tsum(itemp)
         tabl(i,2,itemp)=tabl(i,2,itemp)/denom
         tabl(i,3,itemp)=tabl(i,3,itemp)/denom
         tabl(i,4,itemp)=tabl(i,4,itemp)/denom
         tabl(i,5,itemp)=tabl(i,5,itemp)/denom
         tabl(i,2,itemp)=tabl(i,2,itemp)-tnorm
         tabl(i,3,itemp)=tabl(i,3,itemp)-enorm
         tabl(i,4,itemp)=tabl(i,4,itemp)-fnorm
         tabl(i,5,itemp)=tabl(i,5,itemp)-cnorm
         if (tabl(i,1,itemp).le.zero) izeroprob=izeroprob+1
         if (tabl(i,2,itemp).lt.zero) ibadxs=ibadxs+1
         if (tabl(i,3,itemp).lt.zero) ibadxs=ibadxs+1
         if (tabl(i,4,itemp).lt.zero) ibadxs=ibadxs+1
         if (tabl(i,5,itemp).lt.zero) ibadxs=ibadxs+1
      enddo
   enddo
   do itemp=1,ntemp
      do ixx=1,4
         if (ixx.eq.1) then
            write(nsyso,'('' tmin'',1p,e11.3,1p,10e11.3)')&
              temp(itemp),tmin(itemp)
            write(nsyso,'('' tmax'',1p,e11.3,1p,10e11.3/(16x,10e11.3))')&
              temp(itemp),(tval(i,itemp),i=1,nbin)
            write(nsyso,'('' prob'',1p,e11.3,1p,10e11.3/(16x,10e11.3))')&
              temp(itemp),(tabl(i,1,itemp),i=1,nbin)
         endif
         write(nsyso,'(1x,a,1x,1p,e11.3,10e11.3/(16x,10e11.3))')&
           nmr(ixx),temp(itemp),(tabl(i,ixx+1,itemp),i=1,nbin)
      enddo
   enddo
   if (izeroprob.gt.0) then
      write(strng,'(''ptable has '',i3,'' zero probability bins'')') izeroprob
      call mess('purr',strng,' ')
   endif
   if (ibadxs.gt.0) then
      write(strng,'(''ptable has '',i3,'' negative xs values'')') ibadxs
      call mess('purr',strng,' ')
   endif

   !--compute bondarenko cross sections from prob. table
   if (iprint.gt.0) write(nsyso,'(/&
     &'' bondarenko cross sections from probability table''/&
     &9x,''temp'',6x,''sig0'',4x,''p0 total'',5x,''elastic'',&
     &5x,''fission'',5x,''capture'',4x,''p1 total'')')
   do i=1,nsig0
      do j=1,7
         do itemp=1,ntemp
            bval(j,i,itemp)=0
         enddo
      enddo
      do itemp=1,ntemp
         do j=1,nbin
            if (tabl(j,1,itemp).ne.zero) then
               den=sig0(i)/(sig0(i)+tabl(j,2,itemp))
               ttt=tabl(j,1,itemp)
               bval(1,i,itemp)=bval(1,i,itemp)+ttt*tabl(j,2,itemp)*den
               bval(2,i,itemp)=bval(2,i,itemp)+ttt*tabl(j,3,itemp)*den
               bval(3,i,itemp)=bval(3,i,itemp)+ttt*tabl(j,4,itemp)*den
               bval(4,i,itemp)=bval(4,i,itemp)+ttt*tabl(j,5,itemp)*den
               bval(5,i,itemp)=bval(5,i,itemp)+ttt*tabl(j,2,itemp)*den*den
               bval(6,i,itemp)=bval(6,i,itemp)+ttt*den
               bval(7,i,itemp)=bval(7,i,itemp)+ttt*den*den
            endif
         enddo
         sigf(1,i,itemp)=bval(1,i,itemp)/bval(6,i,itemp)
         sigf(2,i,itemp)=bval(2,i,itemp)/bval(6,i,itemp)
         sigf(3,i,itemp)=bval(3,i,itemp)/bval(6,i,itemp)
         sigf(4,i,itemp)=bval(4,i,itemp)/bval(6,i,itemp)
         sigf(5,i,itemp)=bval(5,i,itemp)/bval(7,i,itemp)
      enddo
   enddo
   if (iprint.gt.0) then
      do itemp=1,ntemp
         do i=1,nsig0
            write(nsyso,'(3x,1p,2e10.3,5e12.4)')&
              temp(itemp),sig0(i),(sigf(j,i,itemp),j=1,5)
         enddo
      enddo
   endif

   !--renormalize probability table and bondarenko
   !--cross sections to the computed infinitely-dilute values
   do itemp=1,ntemp
      do i=2,5
         do j=1,nbin
            l=4*(itemp-1)
            if (sigf(i-1,1,itemp).ne.0)&
              tabl(j,i,itemp)=tabl(j,i,itemp)*sigi(i-1)/sigf(i-1,1,itemp)
         enddo
      enddo
   enddo
   do itemp=1,ntemp
      do i=nsig0,1,-1
         do j=1,5
            k=j
            if (j.eq.5) k=1
            if (sigf(j,1,itemp).ne.0)&
              sigf(j,i,itemp)=sigf(j,i,itemp)*sigi(k)/sigf(j,1,itemp)
            sigpl(ipl,j,i)=sigf(j,i,1)
         enddo
      enddo
   enddo
     !optional printout for plotting probability per barn with viewr
     !uncomment the following lines to activate
     !write(nsyso,'(/'' probability per barn versus total'')')
     !write(nsyso,'('' e='',1p,e12.4)') e
     !tnorm=0
     !do i=1,nbin
     !   denom=tval(i,1)-tnorm
     !   denom=tabl(i,1,1)/denom
     !   write(nsyso,'(1p,2e11.4,''/'')') tabl(i,2,1),denom
     !   tnorm=tval(i,1)
     !enddo
   return
   end subroutine unrest

   subroutine uwtab2
   !-------------------------------------------------------------------
   ! Subroutine wtable controls the calculation and writing of the
   ! coarse and fine, real and imaginary parts of the w table.
   !-------------------------------------------------------------------
   ! internals
   integer::i,j
   real(kr)::delx,dely,x,y,rew,aimw
   real(kr)::ax(41),ay(27)
   real(kr),parameter::tenth=.1e0_kr
   real(kr),parameter::oh2=.02e0_kr

   ax(1)=-tenth
   ax(2)=0
   delx=tenth
   do i=3,41
      ax(i)=ax(i-1)+delx
   enddo
   ay(1)=4*tenth
   ay(2)=5*tenth
   dely=tenth
   do j=3,27
      ay(j)=ay(j-1)+dely
   enddo
   do i=2,41
      x=ax(i)
      do j=1,27
         y=ay(j)
         call uw2(x,y,rew,aimw)
         tr(i,j)=rew
         ti(i,j)=aimw
      enddo
   enddo
   do j=1,27
      tr(1,j)=tr(3,j)
      ti(1,j)=-ti(3,j)
      ti(2,j)=0
   enddo
   ay(1)=-oh2
   ay(2)=0
   dely=oh2
   do j=3,27
      ay(j)=ay(j-1)+dely
   enddo
   do i=2,41
      x=ax(i)
      do j=1,27
         y=ay(j)
         call uw2(x,y,rew,aimw)
         trs(i,j)=rew
         tis(i,j)=aimw
      enddo
   enddo
   do j=1,27
      trs(1,j)=trs(3,j)
      tis(1,j)=-tis(3,j)
      tis(2,j)=0
   enddo
   return
   end subroutine uwtab2

   subroutine uw2(rez,aim1,rew,aimw)
   !-------------------------------------------------------------------
   ! Compute complex probability integral w.
   !-------------------------------------------------------------------
   use physics ! provides pi
   ! externals
   real(kr)::rez,aim1,rew,aimw
   ! internals
   integer::kw,k,jsig
   real(kr)::rpi,aimz,abrez,ai2,r2,rv,ak,el,h,b,a,tempm,temel
   real(kr)::g,c,d,am,aak,ajtemp,temp4,ajp,pr,pim,amagn,tempc,tempd
   real(kr)::temp1,temp2,aj,ajsig,sigp,expon,expc,exps,sig2p
   real(kr)::aj4sig,aj4sm1,temp3,tt4,temp7,ref,aimf
   real(kr),parameter::c1=1.1283792e0_kr
   real(kr),parameter::c2=1.5e0_kr
   real(kr),parameter::up=1.e15_kr
   real(kr),parameter::dn=1.e-15_kr
   real(kr),parameter::brk1=1.25e0_kr
   real(kr),parameter::brk2=5.0e0_kr
   real(kr),parameter::brk3=1.863636e0_kr
   real(kr),parameter::brk4=4.1e0_kr
   real(kr),parameter::brk5=1.71e0_kr
   real(kr),parameter::brk6=2.89e0_kr
   real(kr),parameter::brk7=1.18e0_kr
   real(kr),parameter::brk8=5.76e0_kr
   real(kr),parameter::brk9=1.5e0_kr
   real(kr),parameter::eps=1.e-7_kr
   real(kr),parameter::zero=0

   rpi=sqrt(pi)
   rew=0
   aimw=0
   aimz=abs(aim1)
   abrez=abs(rez)
   if (abrez+aimz.ne.zero) go to 20
   rew=1
   aimw=0
   return

  20 continue
   r2=rez*rez
   ai2=aimz*aimz
   if (abrez+brk1*aimz-brk2.gt.zero) go to 350
   if (abrez+brk3*aimz-brk4.gt.zero) go to 350
   if (r2+brk5*ai2-brk6.lt.zero) go to 340
   if (r2+brk7*ai2-brk8.ge.zero) go to 340
   if (aimz-brk9.ge.zero) go to 350
  340 continue
   kw=2
   aimz=aim1
   go to 420
  350 continue
   kw=1
   if (aim1.ge.zero) go to 370
   kw=2
   aimz=aim1
   go to 420
   ! wa is obtained from asymtotic series
  370 continue
   rv=2*(r2-ai2)
   ak=4*rez*aimz
   el=ak
   h=0
   b=0
   a=0
   tempm=0
   temel=0
   g=1
   c=-c1*aimz
   d=c1*rez
   am=rv-1
   aak=1
   k=0
  380 continue
   ajtemp=2*aak
   temp4=(1-ajtemp)*ajtemp
   ajp=rv-(4*aak+1)
   go to 480
  390 continue
   aak=aak+1
   k=k+1
   pr=rew
   pim=aimw
   amagn=tempm**2+temel**2
   rew=(tempc*tempm+tempd*temel)/amagn
   aimw=(tempm*tempd-temel*tempc)/amagn
   if (abs(rew-pr).ge.eps) go to 380
   if (rez.eq.zero) go to 550
   if (abs(aimw-pim).ge.eps) go to 380
   return
   ! wt is obtained from taylor series
  420 continue
   temp1=r2+ai2
   temp2=2*temp1*temp1
   aj=-(r2-ai2)/temp2
   ak=2*rez*aimz/temp2
   c=0
   b=0
   ajsig=0
   d=0
   jsig=0
   g=0
   h=0
   el=0
   a=1
   am=1
   sigp=c2
   expon=exp(temp2*aj)
   expc=expon*cos(temp2*ak)
   exps=-expon*sin(temp2*ak)
   sig2p=2*sigp
  430 continue
   aj4sig=4*ajsig
   aj4sm1=aj4sig-1
   temp3=1/(aj4sm1*(aj4sig+3))
   tt4=sig2p*(2*ajsig-1)
   temp4=tt4/(aj4sm1*(aj4sig+1)*(aj4sig-3)*aj4sm1)
   ajp=aj+temp3
   go to 480
  440 continue
   ajsig=ajsig+1
   jsig=jsig+1
   temp7=rpi*(am*am+el*el)
   ref=(aimz*(c*am+d*el)-rez*(am*d-c*el))/temp7/temp1
   aimf=(aimz*(am*d-c*el)+rez*(c*am+d*el))/temp7/temp1
   pr=rew
   pim=aimw
   rew=expc-ref
   aimw=exps-aimf
   if (abs(rew-pr).ge.eps) go to 470
   if (rez.eq.zero) go to 550
   if (abs(aimw-pim).ge.eps) go to 470
   return
  470 continue
   sig2p=2*ajsig
   go to 430
  480 continue
   tempc=ajp*c+temp4*a-ak*d
   tempd=ajp*d+temp4*b+ak*c
   temel=ajp*el+temp4*h+ak*am
   tempm=ajp*am+temp4*g-ak*el
   a=c
   b=d
   g=am
   h=el
   c=tempc
   d=tempd
   am=tempm
   el=temel
   if (abs(tempm)+abs(temel)-up.lt.zero) go to 500
   c=dn*c
   d=dn*d
   am=dn*am
   el=dn*el
   tempc=dn*tempc
   tempd=dn*tempd
   tempm=dn*tempm
   temel=dn*temel
   go to 520
  500 continue
   if (abs(tempm)+abs(temel)-dn.gt.zero) go to 520
   c=up*c
   d=up*d
   am=up*am
   el=up*el
   tempc=up*tempc
   tempd=up*tempd
   tempm=up*tempm
   temel=up*temel
  520 continue
   go to (390,440,530),kw
  530 return
  550 continue
   aimw=0
   return
   end subroutine uw2

   subroutine fsort(x,y,n,i)
   !-------------------------------------------------------------------
   ! Floating-point sort routine.
   ! Sort x and y into increasing x order.
   !-------------------------------------------------------------------
   ! externals
   integer::n,i
   real(kr)::x(n),y(n)
   ! internals
   integer::k,j
   real(kr)::xt,yt

   do  k=1,n-1
      do j=k+1,n
         if (x(k).gt.x(j)) then
            xt=x(k)
            yt=y(k)
            x(k)=x(j)
            y(k)=y(j)
            x(j)=xt
            y(j)=yt
         endif
      enddo
   enddo
   return
   end subroutine fsort

   subroutine fsrch(x,xarray,n,i,k)
   !-------------------------------------------------------------------
   ! Floating-point search routine.
   ! Search the xarray for x.ge.xarray(i) and x.lt.xarray(i+1).
   ! Return i=1 and k=2 if x is below the lower limit.
   ! Return i=n and k=3 if x is above the upper limit.
   !-------------------------------------------------------------------
   ! externals
   integer::n,i,k
   real(kr)::x,xarray(n)
   ! internals
   integer::i1,i2,idone

   if (x.lt.xarray(1)) then
      i=1
      k=2
      return
   endif
   if (x.gt.xarray(n)) then
      i=n
      k=3
      return
   endif
   i1=1
   i2=n
   idone=0
   do while (idone.eq.0)
      if (i1+1.eq.i2) then
         idone=1
      else
         if (x.ge.xarray((i1+i2)/2)) then
            i1=(i1+i2)/2
         else
            i2=(i1+i2)/2
         endif
      endif
   enddo
   i=i1
   k=1
   return
   end subroutine fsrch

   real(kr) function rann(idum)
   !-------------------------------------------------------------------
   ! Random number generator.
   ! Set idum negative to reset the seed.
   !-------------------------------------------------------------------
   use util ! provides error
   ! externals
   integer::idum
   ! internals
   integer,parameter::m=714025
   integer,parameter::ia=1366
   integer,parameter::ic=150889
   real(kr),parameter::zero=0
   real(kr),parameter::one=1
   real(kr),parameter::rm=one/m
   integer::j,iy
   integer::ir(97)
   integer::iff=0
   save iff,iy,ir

   ! Do not allow zero
100 continue
   if (idum.lt.0.or.iff.eq.0) then
      iff=1
      idum=mod(ic-idum,m)
      do j=1,97
         idum=mod(ia*idum+ic,m)
         ir(j)=idum
      enddo
      idum=mod(ia*idum+ic,m)
      iy=idum
   endif
   j=1+(97*iy)/m
   if (j.gt.97.or.j.lt.1) call error('rann','failed',' ')
   iy=ir(j)
   rann=iy*rm
   idum=mod(ia*idum+ic,m)
   ir(j)=idum
   if (rann.eq.zero) go to 100
   return
   end function rann

end module purm

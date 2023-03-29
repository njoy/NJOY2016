module covm
   ! provides subroutine covr for NJOY2016
   use locale
   implicit none
   private
   public covr

   ! global variables

   ! user's input
   integer::nin,nout
   integer::nplot
   integer::icolor
   real(kr)::epmin
   real(kr)::einc
   integer::irelco
   integer::ncase
   integer::noleg
   integer::nstart
   integer::ndiv
   integer::matype
   character(12)::hlibid
   character(21)::hdescr
   integer::mat,mt,mat1,mt1
   integer::mfflg,mf3,mf5,mf35

   ! other common variables

   integer,parameter::ncasemx=100
   integer::nscr,nscr1,nscr2
   integer::nrow,ncol,itype,nvf,ncf
   integer::ixmin,ixmax
   integer::izero,ismall
   integer::ishade
   integer::nlev
   integer::nexp
   integer::iza,iza1,izap
   integer::nwscr,max,nwcf,nwig
   real(kr),dimension(:),allocatable::scr
   real(kr),dimension(:),allocatable::xlev
   real(kr),dimension(:),allocatable::cf
   real(kr),dimension(:),allocatable::x,y
   real(kr),dimension(:),allocatable::xx,xy
   real(kr),dimension(:),allocatable::rsdx,rsdy
   integer,dimension(:),allocatable::imtx,imtx1,imatx1

contains

   subroutine covr
   !--------------------------------------------------------------------
   !
   ! Plot covariance data from errorr or make a condensed library.
   !
   ! In the plot option, covr plots a matrix of correlation
   ! coefficients and an associated pair of standard deviation
   ! vectors, i.e.,a covariance matrix.  The correlation
   ! matrix is plotted as a shaded contour plot and the vectors
   ! are plotted as semi-log plots, one rotated by 90 degrees.
   ! The log energy grids for the vector plots are identical
   ! to the grids for the matrix plot.  This version plots
   ! through viewr.
   !
   ! In the library option, covr produces a condensed bcd
   ! covariance library in the boxer format.  This format is
   ! efficient for matrices of simple blocks.
   !
   !---input specifications (free format)---------------------------
   !
   !  card 1
   !     nin            input tape unit
   !     nout           output tape unit
   !                    (default=0=none)
   !     nplot          viewr output unit
   !                    (default=0=none)
   !
   !   ---cards 2, 2', 2a, and 3a for nout.le.0 only (plot option)
   !
   !  card 2
   !     icolor         select color or monochrome style
   !                      0=monochrome (uses cross hatching)
   !                      1=color background and contours
   !                      2=color background and contours plus
   !                        card 2' follows.
   !                      (default=0)
   !  card 2' (only when icolor=2)
   !     nlev,(tlev(i),i=1,nlev)
   !                    defines the number of correlation matrix
   !                    intervals and their boundaries.  Zero is
   !                    assumed as the lower limit of the first
   !                    boundary, but the User must specify unity
   !                    as the upper limit of the last boundary.
   !                    nlev is a positive integer .le. 9.
   !                    default values (when icolor=1) are:
   !                      6,0.001,0.1,0.2,0.3,0.6,1.0
   !  card 2a
   !     epmin          lowest energy of interest (default=0.)
   !  card 3a
   !     irelco         type of covariances present on nin
   !                    0/1=absolute/relative covariances
   !                    (default=1)
   !     ncase          no. cases to be run (maximum=60)
   !                    (default=1)
   !     noleg          plot legend option
   !                    -1/0/1=legend for first subcase only/
   !                    legend for all plots/no legends
   !                    (default=0)
   !     nstart         sequential figure number
   !                    0/n=not needed/first figure is figure n.
   !                    (default=1)
   !     ndiv           no. of subdivisions of each of the
   !                    gray shades (default=1)
   !
   !   ---cards 2b, 3b, and 3c for nout gt 0 (library option) only--
   !
   !  card 2b
   !     matype         output library matrix option
   !                    3/4=covariances/correlations
   !                    (default=3)
   !     ncase          no. cases to be run (maximum=ncasemx=100)
   !                    (default=1)
   !  card 3b
   !     hlibid         up to 6 characters for identification
   !  card 3c
   !     hdescr         up to 21 characters of descriptive
   !                    information
   !
   !   ---cards 4 for both options---
   !
   !  card 4
   !     mat            desired mat number
   !     mt             desired mt number
   !     mat1           desired mat1 number
   !     mt1            desired mt1 number
   !                    (default for mt, mat1 and mt1 are 0,0,0
   !                    meaning process all mts for this mat
   !                    with mat1=mat)
   !                    (neg. values for mt, mat1, and mt1 mean
   !                    process all mts for this mat, except for
   !                    the mt-numbers -mt, -mat1, and -mt1.  in
   !                    general, -n will strip both mt=1 and mt=n.
   !                    -4 will strip mt=1, mt=3, and mt=4, and
   !                    -62, for example, will strip mt=1, mt=62,
   !                    mt=63, ... up to and incl. mt=90.)
   !          repeat card 4 ncase times
   !
   ! note---if more than one material appears on the input tape,
   ! the mat numbers must be in ascending order.
   !
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso,nsyse
   use util ! provides timer,openz,error,repoz,mess
   use endf ! provides endf routines and variables
   ! internals
   integer,parameter::nlevmx=9
   integer::iomit,nfigmx,ncamx,ne,i,nfig,nwlev,ilev,j,n,nb,nw
   integer::nwy,ixpn
   real(kr)::sec,da,tlow,tx
   character(60)::strng
   character(6)::hinpid
   real(kr),dimension(:),allocatable::xig,yig
   real(kr),dimension(nlevmx)::tlev=(/.001e0_kr,.1e0_kr,.2e0_kr,&
             .3e0_kr,.6e0_kr,1.0e0_kr,0.0e0_kr,0.0e0_kr,0.0e0_kr/)
   integer::imat(ncasemx),imt(ncasemx),imat1(ncasemx),imt1(ncasemx)
   real(kr),parameter::rdn=0.999998e0_kr
   character(5),dimension(3)::hreas=(/'null ','small','empty'/)
   integer,parameter::ntics3=600

   !--initialize
   call timer(sec)
   write(nsyso,'(/&
     &'' covr...process covariance data'',38x,f8.1,''s'')') sec
   write(nsyse,'(/'' covr...'',61x,f8.1,''s'')') sec
   nfigmx=ncasemx
   max=1+nfigmx*(nfigmx+1)/2
   mf3=3
   mf5=5
   ncamx=ncasemx
   iomit=0

   !--read and write out user input
   nout=0
   nplot=0
   read(nsysi,*) nin,nout,nplot
   call openz(nin,0)
   call openz(nout,1)
   if (nout.le.0) then
      nscr1=11
      nscr2=0
      call openz(nscr1,1)
      call repoz(-nscr1)
      write(nscr1,'(/&
        &'' the following plots were omitted or were plotted but '',&
        &''empty'',21x/4x,''mat'',5x,''mt'',3x,''mat1'',4x,''mt1'',3x,&
        &''reason'',43x/3x,''----'',4x,''---'',3x,''----'',4x,&
        &''---'',3x,''------'',43x)')
      ! input for plot option
      icolor=0
      read(nsysi,*) icolor
      nlev=6
      if (icolor.eq.2) then
         read(nsysi,*)nlev,(tlev(i),i=1,nlev)
         if (nlev.gt.nlevmx) then
           write(strng,'('' nlev>'',i2,'' not allowed'')')nlevmx
            call error('covr',strng,' ')
         endif
         do i=2,nlev
            if (tlev(i).le.tlev(i-1)) then
               call error('covr',&
                          'tlev array must sequentially increase',&
                          ' ')
            endif
         enddo
         if (tlev(nlev).ne.1) then
            write(strng,'(''reset tlev('',i2,'') from '',1pe11.4,&
                          &'' to 1.0'')')nlev,tlev(nlev)
            call mess('covr',strng,' ')
            tlev(nlev)=1
         endif
      endif
      epmin=0
      read(nsysi,*) epmin
      epmin=epmin*rdn
      irelco=1
      ncase=1
      noleg=0
      nstart=1
      ndiv=1
      read(nsysi,*) irelco,ncase,noleg,nstart,ndiv
      if (ndiv.eq.0) ndiv=1
      if (ncase.gt.ncamx)&
        call error('covr','requested too many cases.',' ')
   else
      ! input for library option
      nscr1=11
      nscr2=12
      call openz(-nscr1,1)
      call openz(-nscr2,1)
      call repoz(-nscr1)
      call repoz(-nscr2)
      ! do not translate abs. covariances to rel. in library option,
      ! thus allowing stand. devs. and covariances to track the
      ! errorr output.
      irelco=1
      matype=3
      ncase=1
      read(nsysi,*) matype,ncase
      if (matype.ne.4) matype=3
      if (ncase.le.0) ncase=1
      if (ncase.gt.ncamx)&
        call error('covr','requested too many cases.',' ')
      read(nsysi,*) hinpid
      hdescr=' '
      read(nsysi,*) hdescr
   endif
   do i=1,ncase
      imt(i)=0
      imat1(i)=0
      imt1(i)=0
      read(nsysi,*) imat(i),imt(i),imat1(i),imt1(i)
   enddo
   if (nout.eq.0) write(nsyso,'(/&
     &'' unit for input covariance tape ....... '',i10/&
     &'' unit for output covariance tape ...... '',i10/&
     &'' unit for plot output ................. '',i10/&
     &'' icolor ............................... '',i10/&
     &'' rel. cov. option (0=abs/1=rel) ....... '',i10/&
     &'' ncase ................................ '',i10/&
     &'' legend option (-1=first/0=all/1=none)  '',i10/&
     &'' starting fig. no. (0=omit fig. nos.) . '',i10/&
     &'' no. subdivisions of each shade ....... '',i10/&
     &'' minimum energy to be plotted ......... '',1p,e10.3)')&
     nin,nout,nplot,icolor,irelco,ncase,noleg,nstart,ndiv,epmin
   if (nout.gt.0) write(nsyso,'(/&
     &'' unit for input covariance tape ....... '',i10/&
     &'' unit for output covariance tape ...... '',i10/&
     &'' output matrix option (3=cov./4=corr.)  '',i10/&
     &'' no. cases to run ..................... '',i10/&
     &'' library identification ............... ''/6x,a6/&
     &'' library description .................. ''/6x,a21)')&
     nin,nout,matype,ncase,hinpid,hdescr
   write(nsyso,'(/12x,''mat'',6x,''mt'',6x,''mat1'',6x,''mt1''/&
     &11x,''----'',5x,''---'',6x,''----'',6x,''---''/&
     &(11x,i4,5x,i3,6x,i4,6x,i3))')&
     (imat(i),imt(i),imat1(i),imt1(i),i=1,ncase)

   !--initialize parameters for plots
   if (nout.le.0) then
      nfig=nstart-1
      ! set up shade levels
      nwlev=nlev*ndiv
      allocate(xlev(nwlev))
      ilev=0
      tlow=0
      do i=1,nlev
         tx=(tlev(i)-tlow)/ndiv
         do j=1,ndiv
            ilev=ilev+1
            xlev(ilev)=tlow+tx*j
         enddo
         tlow=tlev(i)
      enddo
      da=1
      da=da/1000
      xlev(ilev)=xlev(ilev)+da
      nlev=ilev
   endif
   call repoz(nin)
   nscr=10
   nsc=0
   if (nin.lt.0) nscr=-nscr
   call openz(nscr,1)
   nwscr=npage*2+50

   !--read the first two records to verify that this is a legal tape
   !--set the mf35 flag and determine the number of groups for scr
   allocate(scr(17))
   call tpidio(nin,0,0,scr,nb,nw)
   call contio(nin,0,0,scr,nb,nw)
   if (mfh.ne.1.or.mth.ne.451) call error('covr','illegal input tape',' ')
   mfflg=nint(scr(5))
   if (mfflg.eq.-11.or.mfflg.eq.-14) then
      mf35=mf3
   else if (mfflg.eq.-12) then
      mf35=mf5
   else
      call error('covr','illegal errorr output tape for covr',' ')
   endif
   call contio(nin,0,0,scr,nb,nw)
   nwscr=l1h+7
   if (nwscr.lt.17) nwscr=17
   deallocate(scr)
   allocate(scr(nwscr))
   call repoz(nin)
   call openz(nplot,1)
   if (nplot.ne.0) write(nplot,'(''1 2 .22'',i3,''/'')') icolor

   !--loop over cases
   do 130 n=1,ncase
   mat=imat(n)

   ! copy this mat from nin to nscr
   call repoz(nin)
   call repoz(nscr)
   call tpidio(nin,0,nscr,scr,nb,nw)
   call finds(mat,0,0,nin)
   call contio(nin,0,nscr,scr,nb,nw)
   iverf=l1h
   iza=nint(scr(1))
   iza1=iza
   call tomend(nin,0,nscr,scr)
   allocate(imtx(max))
   allocate(imtx1(max))
   allocate(imatx1(max))
   nexp=1
   imtx(1)=imt(n)
   imatx1(1)=imat1(n)
   imtx1(1)=imt1(n)
   if (imat1(n).ne.imat(n).and.imat1(n).gt.0) then
      ! if needed, also copy mat1 from nin to nscr
      call finds(imat1(n),0,0,nin)
      call contio(nin,0,nscr,scr,nb,nw)
      iza1=nint(scr(1))
      call tomend(nin,0,nscr,scr)
   endif
   call atend(0,nscr)
   call repoz(nscr)

   !--expand the mt-mt1 list
   if (imt(n).le.0) call expndo(nscr)

   !--loop over mt-mt1 pairs
   do 191 ne=1,nexp
   mt=imtx(ne)
   mat1=imatx1(ne)
   if (mat1.le.0) mat1=mat
   mt1=imtx1(ne)
   if (mt1.eq.0) mt1=mt
   write(nsyso,'(/&
     &'' reaction-pairs processed '',i4,2x,i3,3x,i4,2x,i3)')&
     mat,mt,mat1,mt1

   !--read through input errorr file, calculating relative standard
   !--deviations and correlation coefficients
   call corr(nscr,mat,mt,mat1,mt1,izap)
   if (nout.gt.0) go to 210
   if (ismall.gt.0) go to 155
  150 continue
   iomit=iomit+1
   ! document null and small matrices on nscr1, then skip plot
   write(nscr1,'(4i7,3x,a5,44x)') mat,mt,mat1,mt1,hreas(izero+1)
   go to 190
  155 continue
   nwcf=ixmax*ixmax
   nwig=2*(2*(ixmax+1)+ntics3)
   allocate(xig(nwig))
   nwy=nwig
   if (nwig.lt.nwcf) nwy=nwcf
   allocate(yig(nwy))
   ! cf storage is in the middle.  move it to the top
   do i=1,nwcf
      yig(i)=cf(i)
   enddo
   deallocate(cf)
   call truncg(epmin,x,y,xx,xy,ixmin,ixmax)
   ! cf storage goes into highest part
   allocate(cf(nwcf))
   do i=1,nwcf
      cf(i)=yig(i)
      yig(i)=0
   enddo
   if (izero.eq.0) go to 190
   ixpn=ixmax-ixmin+1
   ishade=0
   call plotit(x(ixmin:),y(ixmin:),xig,yig,ixpn,ixmax,&
     rsdx(ixmin:),rsdy(ixmin:),xlev,noleg,ne,izap)
   if (ishade.eq.0) izero=2
   if (ishade.eq.0) go to 150
   if (nfig-nstart+1.lt.nfigmx) go to 190
   write(strng,'(''mat='',i4,'' mt='',i3,'' mat1='',i4,'' mt1='',i3,&
     &'' nfig='',i3)') mat,mt,mat1,mt1,nfig
   call mess('covr','have plotted all that fit.',strng)
   go to 310

   !--process library option
  210 continue
   deallocate(y)
   deallocate(xy)
   if (izero.ne.0) go to 220
   write(nsyso,'(/&
     &3x,''null covariance matrix.  output suppressed.'')')
   go to 190
  220 continue
   if (ne.le.1) then
      if (matype.eq.3) write(hlibid,'(a6,''-a-'',i3)') hinpid,ixmax
      if (matype.eq.4) write(hlibid,'(a6,''-b-'',i3)') hinpid,ixmax
   endif
   nrow=ixmax+1
   ncol=1
   nvf=10
   ncf=3
   itype=0
   ! for first sub-case of first case, write the group structure
   if (n.eq.1.and.ne.eq.1) call press(0,nout,x,nrow,ncol)
   nrow=ixmax
   if (mat1.eq.mat.and.mt1.eq.mt) then
      itype=1
      call press(0,nout,xx,ixmax,1)
      itype=2
      call press(0,nout,rsdx,ixmax,1)
   endif
   ncol=nrow
   if (mat1.eq.mat.and.mt1.eq.mt) ncol=0
   nvf=7
   if (matype.eq.3) nvf=10
   ncf=6
   if (ixmax.le.100) ncf=5
   if (ixmax.le.30) ncf=4
   itype=matype
   call press(0,nout,cf,2,ixmax)
  190 continue
   if (allocated(x))deallocate(x)
   if (allocated(y))deallocate(y)
   if (allocated(rsdx)) deallocate(rsdx)
   if (allocated(rsdy)) deallocate(rsdy)
   if (allocated(xig)) deallocate(xig)
   if (allocated(yig)) deallocate(yig)
   if (allocated(xx)) deallocate(xx)
   if (allocated(xy)) deallocate(xy)
   if (allocated(cf)) deallocate(cf)
  191 continue
   deallocate(imtx)
   deallocate(imtx1)
   deallocate(imatx1)
  130 continue
   if (nout.eq.0) go to 310
   nscr1=-nscr1
   nscr2=-nscr2
   go to 330

   !--plot the table of contents
  310 continue
   write(nplot,'(''99/'')')

   !--copy the statistics file to output
   if (iomit.ne.0) then
      write(nsyso,'(/)')
      call copyst(nscr1,nsyso)
   endif

   !--finished
  330 continue
   if (allocated(scr)) deallocate(scr)
   if (allocated(xlev)) deallocate(xlev)
   call repoz(nin)
   call repoz(nout)
   call closz(nin)
   call closz(nout)
   call closz(nplot)
   call closz(nscr)
   call closz(nscr1)
   call closz(nscr2)
   call timer(sec)
   write(nsyso,'(69x,f8.1,''s''/&
     &1x,7(''**********''),''****'')') sec
   return
   end subroutine covr

   subroutine expndo(nin)
   !--------------------------------------------------------------------
   ! Expand the MT-MT1 list to include all possible combinations
   ! in ENDF order.  Obtain the set of MT-numbers present by
   ! examining the cross section file on unit nin.
   !--------------------------------------------------------------------
   use util ! provides error,mess,repoz
   use endf ! provides endf routines and variables
   ! externals
   integer::nin
   ! internals
   integer::nlstm,nb,nw,nmt,l,mtn,istr,im,jm
   character(60)::strng
   integer::mstrip(3)
   integer,dimension(:),allocatable::lstm

   nlstm=ncasemx
   allocate(lstm(nlstm))
   call repoz(nin)
   call tpidio(nin,0,0,scr,nb,nw)
   nmt=0
   mstrip(1)=-imtx(1)
   mstrip(2)=-imatx1(1)
   mstrip(3)=-imtx1(1)
   call finds(mat,mf35,0,nin)
  130 continue
   call listio(nin,0,0,scr,nb,nw)
   if (mfh.ne.mf35) go to 160
   l=1
   mtn=mth
   do while (nb.ne.0)
      l=l+nw
      if (l.gt.nwscr)&
        call error('expndo','storage exceeded.',' ')
      call moreio(nin,0,0,scr(l),nb,nw)
   enddo
   ! read send card
   call contio(nin,0,0,scr,nb,nw)
   do 155 istr=1,3
   if (mstrip(istr).eq.0) go to 155
   if (mtn.eq.1) go to 180
   if (mtn.eq.mstrip(istr)) go to 180
   if (mtn.eq.3.and.mstrip(istr).eq.4) go to 180
   if (mstrip(istr).lt.51.or.mstrip(istr).gt.90) go to 155
   if (mtn.gt.mstrip(istr).and.mtn.le.90) go to 180
  155 continue
   nmt=nmt+1
   lstm(nmt)=mtn
   go to 130
  160 continue
   nexp=0
   do im=1,nmt
      do jm=im,nmt
         nexp=nexp+1
         if (nexp.gt.max) then
            call error('expndo','storage exceeded.',' ')
         endif
         imtx(nexp)=lstm(im)
         imatx1(nexp)=0
         imtx1(nexp)=lstm(jm)
      enddo
   enddo
   deallocate(lstm)
   return
  180 continue
   write(strng,'(''mt='',i3,'' stripped from list.'')') mtn
   call mess('expndo',strng,' ')
   go to 130
   end subroutine expndo

   subroutine corr(nin,mat,mt,mat1,mt1,izap)
   !--------------------------------------------------------------------
   ! Convert relative covariances to standard deviations
   ! and correlations.
   !--------------------------------------------------------------------
   use util ! provides error,repoz
   ! externals
   integer::nin,mat,mt,mat1,mt1,izap
   ! internals
   integer::icall,ixtest,ixp,i,ig,ii,j,nswap
   real(kr)::epmn,xcycle
   character(60)::strng
   real(kr),parameter::zp4=0.4e0_kr
   real(kr),parameter::xsize=4.25e0_kr
   real(kr),parameter::epmn0=0.499999e-4_kr
   real(kr),parameter::zero=0

   ismall=0

   !--for cross-reaction matrices, check to see if the matrix
   !--is null before proceeding with the correlation calculation
   icall=0
   if (mat.ne.mat1.or.mt.ne.mt1) then
      call covard(nin,mat,mt,mat1,mt1,izap,icall)
      icall=1
      if (izero.eq.0) return
   endif

  110 continue
   call covard(nin,mat1,mt1,mat1,mt1,izap,icall)
   ixtest=ixmax
   if (icall.eq.0) icall=1
   allocate(rsdx(ixmax))
   allocate(rsdy(ixmax))
   icall=iabs(icall)
   go to (130,165,190),icall

  130 continue
   ixp=ixmax+1
   do i=1,ixp
      x(i)=y(i)
   enddo
   if (nout.le.0) then
      epmn=epmn0
      xcycle=0
      do while (xcycle.lt.zp4)
         epmn=epmn*2
         xcycle=xsize/log10(x(1+ixmax)/epmn)
      enddo
      ! have found epmin value that satifies criterion
      if (epmn.gt.epmin) then
         epmin=epmn
         write(strng,&
           '(''epmin reset to'',1p,e12.4,'', xcycle='',1p,e12.4)')&
           epmin,xcycle
         call mess('corr',strng,' ')
      endif
      do i=1,ixmax
         if (cf(ixmax*(i-1)+i).gt.0) then
            rsdy(i)=sqrt(cf(ixmax*(i-1)+i))
         else
            rsdy(i)=0
         endif
         rsdx(i)=rsdy(i)
      enddo
   else
      call repoz(-nscr2)
      do i=1,ixmax
         read(nscr2) ig,(cf(ii),ii=1,ixmax)
         rsdy(i)=sqrt(cf(i))
         rsdx(i)=rsdy(i)
      enddo
   endif
   if (mat.eq.mat1.and.mt.eq.mt1) go to 190
   icall=icall+1
   call covard(nin,mat,mt,mat,mt,izap,icall)
   if (icall.lt.0) go to 110

  165 continue
   if (nout.le.0) then
      do i=1,ixmax
         rsdx(i)=sqrt(cf(ixmax*(i-1)+i))
      enddo
   else
      call repoz(-nscr2)
      do i=1,ixmax
         read(nscr2) ig,(cf(ii),ii=1,ixmax)
         rsdx(i)=sqrt(cf(i))
      enddo
   endif
   icall=icall+1
   call covard(nin,mat,mt,mat1,mt1,izap,icall)
   if (icall.lt.0) go to 110

  190 continue
   if (izero.ne.0) then
      if (nout.le.0) then
         do i=1,ixmax
            do j=1,ixmax
               if (cf(ixmax*(i-1)+j).ne.zero.and.&
                 rsdx(i)*rsdy(j).ne.zero) then
                  cf(ixmax*(i-1)+j)=cf(ixmax*(i-1)+j)/&
                    (rsdx(i)*rsdy(j))
                  ! test for the presence of any significant (i.e.,
                  ! plotable) values of the correlation coefficient
                  if (abs(cf(ixmax*(i-1)+j)).ge.xlev(1)) ismall=1
               else
                  cf(ixmax*(i-1)+j)=0
               endif
            enddo
         enddo
      else
         call repoz(-nscr2)
         call repoz(-nscr1)
         if (matype.ne.3) then
            do i=1,ixmax
               read(nscr2) ig,(cf(ii),ii=1,ixmax)
               do j=1,ixmax
                  if (cf(j).ne.zero.and.&
                    rsdx(i)*rsdy(j).ne.zero) then
                     cf(j)=cf(j)/(rsdx(i)*rsdy(j))
                  else
                     cf(j)=0
                  endif
               enddo
               write(nscr1) ig,(cf(ii),ii=1,ixmax)
            enddo
            call repoz(-nscr2)
         else
            nswap=nscr1
            nscr1=nscr2
            nscr2=nswap
         endif
      endif
   endif

   ! finished
   if (ixmax.ne.ixtest)&
     call error('corr','group structures do not agree.',' ')
   return
   end subroutine corr

   subroutine covard(nin,mat,mt,mat1,mt1,izap,icall)
   !--------------------------------------------------------------------
   ! Read covariance data in the ENDF-like format output by
   ! njoy/errorr and transform it into engineering format,
   ! that is, with full matrices (zeroes included).
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso,nsyse
   use util ! provides repoz,finds,error,mess
   use endf ! provides ENDF routines and variables
   ! externals
   integer::nin,mat,mt,mat1,mt1,izap,icall
   ! internals
   integer::nb,nw,l,i,k,nmt,imt,k1,mat1x,mtx,kgp,kl,lgp1
   integer::j1,j2,kj,j,ir,ipflag,ind,kk,n,mf3x,ixp,ixpo
   character(60)::strng,strng2
   integer::igprnt=0
   real(kr),parameter::zero=0
   save ixpo,igprnt

   if (icall.le.1) ixpo=0
   call repoz(nin)
   call tpidio(nin,0,0,scr,nb,nw)
   call finds(mat,1,451,nin)
   if (nout.gt.0) call repoz(-nscr1)
   if (nout.gt.0) call repoz(-nscr2)

   !--read group structure
   call contio(nin,0,0,scr,nb,nw)
   call listio(nin,0,0,scr,nb,nw)
   ixmax=l1h
   ixp=n1h
   if (icall.eq.0) go to 120
   if (icall.eq.1) go to 130
   if (ixpo.eq.0.and.ixp.eq.igprnt) go to 120
   if (ixp.le.ixpo) go to 130
   deallocate(x)
   deallocate(y)
   deallocate(cf)
   deallocate(xx)
   deallocate(xy)
   icall=-icall
  120 continue
   allocate(x(ixp))
   allocate(y(ixp))
   allocate(xx(ixmax))
   allocate(xy(ixmax))
   nwcf=ixmax*(ixmax+3)
   if (nout.gt.0) nwcf=ixmax*2
   allocate(cf(nwcf))
  130 continue
   l=1
   do while (nb.ne.0)
      l=l+nw
      if (l.gt.nwscr)&
        call error('covard','storage exceeded.',' ')
      call moreio(nin,0,0,scr(l),nb,nw)
   enddo
   do i=1,ixp
      y(i)=scr(i+6)
   enddo
   if (igprnt.ne.ixp) write(nsyso,'(/&
     &'' group structure''/'' ---------------''/&
     &(1x,1p,6e12.4))') (y(i),i=1,ixp)
   igprnt=ixp
   ixpo=ixp

   !--read cross sections for mt and mt1
   k=1
   call finds(mat,mf35,mt,nin)
  200 continue
   call listio(nin,0,0,scr,nb,nw)
   if (mf35.eq.5) einc=scr(2)
   l=1
   do while (nb.ne.0)
      l=l+nw
      if (l.gt.nwscr)&
        call error('covard','storage exceeded.',' ')
      call moreio(nin,0,0,scr(l),nb,nw)
   enddo
   if (k.gt.1) go to 240
   do i=1,ixmax
      xx(i)=scr(i+6)
   enddo
   if (mth.ne.mt1) go to 260
  240 continue
   do i=1,ixmax
      xy(i)=scr(i+6)
   enddo
  260 continue
   if (k.gt.1) go to 270
   if (mat.eq.mat1.and.mt.eq.mt1) go to 270
   call finds(mat1,mf35,mt1,nin)
   k=k+1
   go to 200
  270 continue
   do i=1,nwcf
      cf(i)=0
   enddo

   !--read covariances.
   mf3x=33
   if (mfflg.eq.-14) mf3x=40
   if (mf35.eq.5) mf3x=35
   if (mt.eq.251) mf3x=34
   call finds(mat,mf3x,mt,nin)
   call contio(nin,0,0,scr,nb,nw)
   nmt=n2h
   if (nmt.eq.0) go to 450
   imt=0
   k1=0
   ! search for desired subsection of this mt
   do
      imt=imt+1
      if (imt.gt.nmt) then
         write(strng,'(''did not find file '',i2,'' subsection'')')mf3x
         write(strng2,'(''for mt='',i3,'' mat='',i4)') mt1,mat1
         call error('covard',strng,strng2)
      endif
      izap=0
      if (mfflg.eq.-14) then
         if (mat.eq.mat1.and.mt.eq.mt1) then
            izap=nint(c2h)
         endif
      endif
      call contio(nin,0,0,scr,nb,nw)
      mat1x=l1h
      if (mf3x.eq.34) mat1x=math
      if (mat1x.eq.0) mat1x=math
      mtx=l2h
      if (mf3x.eq.34) mtx=l1h
      kgp=n2h
      if (mat1x.eq.mat1.and.mtx.eq.mt1) exit
      do
         call listio(nin,0,0,scr,nb,nw)
         kl=n2h
         do while (nb.ne.0)
            call moreio(nin,0,0,scr,nb,nw)
         enddo
         if (kl.ge.kgp) exit
      enddo
   enddo
   do
      call listio(nin,0,0,scr,nb,nw)
      lgp1=l2h
      kl=n2h
      l=1
      do while (nb.ne.0)
         l=l+nw
         if (l.gt.nwscr)&
           call error('covard','storage exceeded.',' ')
         call moreio(nin,0,0,scr(l),nb,nw)
      enddo
      j1=7
      j2=l+nw-1
      if (nout.le.0) then
         kj=0
         do j=j1,j2
            kj=kj+1
            cf(ixmax*(kl-1)+lgp1-1+kj)=scr(j)
         enddo
      else
         nw=j2-j1+1
         write(nscr1) kl,lgp1,nw,(scr(j),j=j1,j2)
         if (k1.eq.0) k1=kl
      endif
      if (kl.ge.kgp) exit
   enddo
   if (nout.eq.0) go to 410
   call repoz(-nscr1)
   kl=k1
   ir=0

   !--eliminate spurious covariances in regions of zero xsec,
   !--convert absolute covariances to relative,
   !--and test for a null covariance matrix
  410 continue
   ipflag=0
   izero=0
   do 420 k=1,ixmax
   if (nout.eq.0) ind=ixmax*(k-1)
   if (nout.eq.0) go to 425
   ind=0
   do kk=1,ixmax
      cf(kk)=0
   enddo
   if (k.lt.kl) go to 425
   if (ir.eq.0) read(nscr1) kl,lgp1,nw,(scr(i),i=1,nw)
   ir=ir+1
   if (k.lt.kl) go to 425
   ir=0
   do kk=1,nw
      cf(lgp1+kk-1)=scr(kk)
   enddo
  425 continue
   do 428 n=1,ixmax
   if (xx(k)*xy(n).ne.zero) go to 430
   if (cf(ind+n).eq.zero) go to 428
   cf(ind+n)=0
   ipflag=ipflag+1
   if (ipflag.eq.1) then
      write(strng,&
        '(''mt='',i3,'' ig='',i3,'' mt1='',i3,'' ig1='',i3)')&
        mt,k,mt1,n
      call mess('covard',&
        'due to zero xsec, covariance zeroed for',strng)
   endif
   if (ipflag.eq.2) call mess('covard','additional warnings suppressed',' ')
   go to 428
  430 continue
   if (irelco.ne.1) cf(ind+n)=cf(ind+n)/(xx(k)*xy(n))
   if (cf(ind+n).ne.zero) izero=1
  428 continue
   if (nout.ne.0) then
      write(nscr2) k,(cf(kk),kk=1,ixmax)
   endif
  420 continue
  450 return
   end subroutine covard

   subroutine truncg(epmin,x,y,xx,xy,ixmin,ixmax)
   !--------------------------------------------------------------------
   ! Truncate the lower end of the MT and MT1 energy group
   ! structures to eliminate zero cross section regions from plots.
   ! Define ethresh to be the lower of either the actual threshold
   ! or 1 MeV (emin).  Then the plots will be truncated at ethresh
   ! or epmin, whichever is higher.
   !--------------------------------------------------------------------
   use util ! provides error
   ! externals
   integer::ixmin,ixmax
   real(kr)::epmin,x(*),y(*),xx(*),xy(*)
   ! internals
   integer::i,iymin,ielo
   real(kr)::rlimx,rlimy,ethrx,ethry,elo
   real(kr),parameter::emin=0.9999e6_kr
   real(kr),parameter::xslim=1.e-4_kr
   real(kr),parameter::zero=0
   real(kr),parameter::ten=10.e0_kr

   !--define alternative limit as average/xslim
   rlimx=0
   rlimy=0
   ethrx=x(1)
   ethry=y(1)
   do i=2,ixmax
      rlimx=rlimx+xx(i-1)*(x(i)-x(i-1))
      rlimy=rlimy+xy(i-1)*(y(i)-y(i-1))
      if (xx(i-1).le.zero) ethrx=x(i)
      if (xy(i-1).le.zero) ethry=y(i)
   enddo
   if (rlimx.ne.zero.and.rlimy.ne.zero) then
      rlimx=xslim*rlimx/(y(ixmax)-ethry)
      rlimy=xslim*rlimy/(x(ixmax)-ethrx)
   endif
   ixmin=1
   do i=1,ixmax
      if (x(1+i).gt.epmin) then
         if (xx(i).ge.xslim.or.xx(i).ge.rlimx) go to 120
         if (x(1+i).gt.emin) go to 120
      endif
      ixmin=i+1
   enddo
   call error('truncg','bad data.',' ')
  120 continue
   iymin=1
   do i=1,ixmax
      if (y(1+i).gt.epmin) then
         if (xy(i).ge.xslim.or.xy(i).ge.rlimy) go to 140
         if (y(1+i).gt.emin) go to 140
      endif
      iymin=i+1
   enddo
   call error('truncg','bad data.',' ')
  140 continue
   if (iymin.lt.ixmin) ixmin=iymin

   ! limit bin width of thermal groupr to one decade
   if (ixmin.eq.1) then
      if (10*x(1).lt.x(2).and.10*y(1).lt.y(2)) then
         elo=log10(x(2)/10)
         ielo=nint(elo)
         if (x(1).lt.ten**ielo) x(1)=ten**ielo
         elo=log10(y(2)/10)
         ielo=nint(elo)
         if (y(1).lt.ten**ielo) y(1)=ten**ielo
      endif
   endif
   if (x(ixmin).lt.epmin) x(ixmin)=epmin
   if (y(ixmin).lt.epmin) y(ixmin)=epmin

   ! finished
   return
   end subroutine truncg

   subroutine plotit(x,y,xig,yig,ixn,ixmax,rsdx,rsdy,xlev,noleg,ne,izap)
   !--------------------------------------------------------------------
   ! Make desired plots.
   !--------------------------------------------------------------------
   ! externals
   integer::ixn,ixmax,noleg,ne,izap
   real(kr)::x(*),y(*),xig(*),yig(*),rsdx(*),rsdy(*),xlev(*)
   ! internals
   integer::i,ii,j,jj,npts,mtflg,iwarn
   real(kr)::xmin,xmax,ymin,ymax,t,t1,xpos,ypos,wa,yyy
   real(kr)::ydec,yymin,yrange,de,xmin1
   character(80)::strng
   character(1)::qu=''''
   real(kr),parameter::two=2
   real(kr),parameter::five=5
   real(kr),parameter::ten=10
   real(kr),parameter::dx1=1.e0_kr
   real(kr),parameter::dy1=.75e0_kr
   real(kr),parameter::dx2=-.25e0_kr
   real(kr),parameter::dy2=.5e0_kr
   real(kr),parameter::dx3=-.75e0_kr
   real(kr),parameter::dy3=1.9e0_kr
   real(kr),parameter::dy4=.75e0_kr
   real(kr),parameter::yyym=59.999e0_kr
   real(kr),parameter::xsize=5.00e0_kr
   real(kr),parameter::ysize=3.38e0_kr
   real(kr),parameter::yrtest=10
   real(kr),parameter::ylogmn=0.1001e0_kr
   real(kr),parameter::ylogmx=99.9e0_kr
   real(kr),parameter::zero=0

   !--check if mat=mat1 and mt=mt1
   !--if true, move nu-bar, few-group cross section or spectrum/eV into rsdx
   mtflg=0
   if (mat.eq.mat1.and.mt.eq.mt1) then
      mtflg=1
      if (mf35.ne.5) then
         do i=1,ixn
            rsdx(i)=xx(i+(ixmax-ixn))
         enddo
      else
         do i=1,ixn
            de=y(i+1)-y(i)
            rsdx(i)=xx(i+(ixmax-ixn))/de
         enddo
      endif
   endif

   !--determine first rsdx array index with non-zero data
   ii=1
   do while (ii.le.ixn.and.rsdx(ii).le.zero)
      ii=ii+1
   enddo
   !if (ii.gt.1) ii=ii-1

   !--determine first rsdy array index with non-zero data
   jj=1
   do while (jj.le.ixn.and.rsdy(jj).le.zero)
      jj=jj+1
   enddo
   !if (jj.gt.1) jj=jj-1

   !--only plot from higher energy group
   !  (makes correlation matrix appear square)
   if (jj.gt.ii)ii=jj
   j=0
   xig(1)=x(ii)
   yig(1)=0
   ymin=1.0e20_kr
   ymax=-1.0e20_kr
   do i=ii,ixn
      j=j+1
      xig(2*j)=x(i)
      yig(2*j)=rsdx(i)
      if (rsdx(i).gt.zero.and.rsdx(i).lt.ymin) ymin=rsdx(i)
      if (rsdx(i).gt.ymax) ymax=rsdx(i)
      xig(2*j+1)=x(i+1)
      yig(2*j+1)=rsdx(i)
   enddo
   if (ymin.ne.zero) then
      yrange=ymax/ymin
   else
      yrange=0
   endif
   ydec=log10(ymin)
   if (ydec.lt.zero) ydec=ydec-1
   yymin=ten**int(ydec)
   npts=2*j+2
   xig(npts)=x(ixn+1)
   yig(npts)=0
   xmin=xig(1)
   xmin1=xmin
   t=log10(xmin)
   if (t.lt.zero) then
      t1=t-0.99999e0_kr
   else
      t1=t+0.00001e0_kr
   endif
   i=int(t1)
   xmin=ten**i
   if (t-i.gt.log10(two)) xmin=two*ten**i
   if (t-i.gt.log10(five)) xmin=five*ten**i
   xmax=xig(npts)
   t=log10(xmax)
   t=t-t/200
   i=int(t)
   xmax=ten**(i+1)
   if (t-i.lt.log10(five)) xmax=five*ten**i
   if (t-i.lt.log10(two)) xmax=two*ten**i
   xpos=ysize-dy1
   ypos=xsize-dx1
   wa=0
   write(nplot,'(''1 0 1. 1.'',5f8.3,''/'')')&
     xpos,ypos,xsize,ysize,wa
   call smilab(iza,mat,mt,izap,mtflg,strng)
   write(nplot,'(a,a,a,''/'')') qu,strng,qu
   write(nplot,'(''/'')')
   if (yrange.gt.yrtest) then
      write(nplot,'(''4/'')')
   else
      write(nplot,'(''3/'')')
   endif
   write(nplot,'(1p,2e12.4,''/'')') xmin,xmax
   write(nplot,'(a,''.'',a,''/'')') qu,qu
   write(nplot,'(''/'')')
   write(nplot,'(a,''.'',a,''/'')') qu,qu
   write(nplot,'(''/'')')
   write(nplot,'(''0 0 0 3 2/'')')
   write(nplot,'(''0'')')

   !--Write (e,std.dev.), or if the mtflg flag is set, write
   !--(e,few-group-nu-bar) or (e,few-group-xsec), to the plot file.
   !--If plotting the relative standard deviation,
   !--restrict the axis limits to 0.1% to 100% for a log scale
   !--or less than 60% for a linear scale.
   !--For few-group data, there is no axis limitation.
   !--Also set a flag to warn the user when data are modified
   !--to fit on the plot (but don't get confused by zero data).
   iwarn=0
   do i=1,npts
      if (mtflg.eq.0) then
         yyy=100*yig(i)
         if (yrange.gt.yrtest.and.yyy.lt.ylogmn) then
            if (yyy.ne.zero) iwarn=1
            yyy=ylogmn
         endif
         if (yrange.gt.yrtest.and.yyy.gt.ylogmx) then
            iwarn=1
            yyy=ylogmx
         else if (yrange.le.yrtest.and.yyy.gt.yyym) then
            iwarn=1
            yyy=yyym
         endif
      else
         yyy=yig(i)
         if (yrange.gt.yrtest) then
            if (yyy.le.yymin) then
               if (yyy.ne.0) iwarn=1
               yyy=1.0001e0_kr*yymin
            endif
         endif
      endif
      write(nplot,'(1p,2e13.4,''/'')') xig(i),yyy
   enddo
   write(nplot,'(''/'')')

   j=0
   xig(1)=y(ii)
   yig(1)=0
   ymin=1.e20_kr
   ymax=-1.e20_kr
   do i=ii,ixn
      j=j+1
      xig(2*j)=y(i)
      yig(2*j)=rsdy(i)
      if (rsdy(i).gt.zero .and. rsdy(i).lt.ymin) ymin=rsdy(i)
      if (rsdy(i).gt.ymax) ymax=rsdy(i)
      xig(2*j+1)=y(i+1)
      yig(2*j+1)=rsdy(i)
   enddo
   if (ymin.ne.zero) then
      yrange=ymax/ymin
   else
      yrange=0
   endif
   npts=2*j+2
   xig(npts)=y(ixn+1)
   yig(npts)=0
   ymin=xig(1)
   ! set uncertainty plot lower limit to match the cross
   ! section plot so get a square correlation matrix.
   if (mtflg.eq.1) ymin=xmin1
   t=log10(ymin)
   if (t.lt.zero) then
      t1=t-0.99999e0_kr
   else
      t1=t+0.00001e0_kr
   endif
   i=int(t1)
   if (t-i.gt.zero) ymin=ten**i
   if (t-i.gt.log10(two)) ymin=two*ten**i
   if (t-i.gt.log10(five)) ymin=five*ten**i
   ymax=xig(npts)
   t=log10(ymax)
   t=t-t/200
   i=int(t)
   ymax=ten**(i+1)
   if (t-i.lt.log10(five)) ymax=five*ten**i
   if (t-i.lt.log10(two)) ymax=two*ten**i
   xpos=ysize+dy2
   ypos=dx2
   wa=90
   write(nplot,'(''-1 0 1. 1.'',5f8.3,''/'')')&
     xpos,ypos,xsize,ysize,wa
   call smilab(iza1,mat1,mt1,izap,0,strng)
   write(nplot,'(a,a,a,''/'')') qu,strng,qu
   write(nplot,'(''/'')')
   if (yrange.gt.yrtest) then
      write(nplot,'(''4/'')')
   else
      write(nplot,'(''3/'')')
   endif
   write(nplot,'(1p,2e12.4,''/'')') ymin,ymax
   write(nplot,'(a,''.'',a,''/'')') qu,qu
   write(nplot,'(''/'')')
   write(nplot,'(a,''.'',a,''/'')') qu,qu
   write(nplot,'(''/'')')
   write(nplot,'(''0 0 0 3 2/'')')
   write(nplot,'(''0'')')

   !--Write (e, std.dev.) to the plot file.  Restrict the axis limits
   !--to 0.1% to 100% for a log scale or less than 60% for a linear scale.
   do i=1,npts
      yyy=100*yig(i)
      if (yrange.gt.yrtest.and.yyy.lt.ylogmn) then
         if (yyy.ne.zero) iwarn=1
         yyy=ylogmn
      endif
      if (yrange.gt.yrtest.and.yyy.gt.ylogmx) then
         iwarn=1
         yyy=ylogmx
      else if (yrange.le.yrtest.and.yyy.gt.yyym) then
         iwarn=1
         yyy=yyym
      endif
      write(nplot,'(1p,2e13.4,''/'')') xig(i),yyy
   enddo
   write(nplot,'(''/'')')

   !--write the text describing the axes
   xpos=dy3
   ypos=xsize+dx3
   wa=90
   write(nplot,'(''-1 0 1. 1.'',5f8.3,''/'')')&
     xpos,ypos,ysize,1.,wa
   if (mtflg.eq.0) then
      write(nplot,'(a,''#H.75<o>rdinate scale is %'',a,''/'')') qu,qu
      write(nplot,'(a,''#H.75<>relative standard deviation.'',a,''/'')') qu,qu
   else if (mt.eq.251) then
      write(nplot,'(a,''#H.75<o>rdinate scales are % relative'',a,''/'')') qu,qu
      write(nplot,'(a,''#H.75<>standard deviation and mu-bar.'',a,''/'')') qu,qu
   else if (mt.ge.452.and.mt.le.456) then
      write(nplot,'(a,''#H.75<o>rdinate scales are % relative'',a,''/'')') qu,qu
      write(nplot,'(a,''#H.75<>standard deviation and nu-bar.'',a,''/'')') qu,qu
   else if (mf35.eq.5) then
      write(nplot,'(a,''#H.75<o>rdinate scales are % standard'',a,''/'')') qu,qu
      write(nplot,'(a,''#H.75<>deviation and spectrum/eV.'',a,''/'')') qu,qu
   else
      write(nplot,'(a,''#H.75<o>rdinate scales are % relative'',a,''/'')') qu,qu
      write(nplot,'(a,''#H.75<>standard deviation and barns.'',a,''/'')') qu,qu
   endif
   write(nplot,'(''0/'')')
   xpos=xpos+dy4
   write(nplot,'(''-1 0 1. 1.'',5f8.3,''/'')')&
     xpos,ypos,ysize,1.,wa
   write(nplot,'(a,''#H.75<a>bscissa scales are energy (e<v>).'',&
     &a,''/'')') qu,qu
   write(nplot,'(''/'')')
   if (iwarn.ne.0) then
       xpos=xpos+0.45e0_kr
       write(nplot,'(''0/'')')
       write(nplot,'(''-1 0 1. 1.'',5f8.3,''/'')')xpos,ypos,ysize,1.,wa
       write(nplot,'(a,''#H.75<w>arning:  some uncertainty'',a,''/'')') qu,qu
       write(nplot,'(a,''#H.75<>data were suppressed.'',a,''/'')') qu,qu
   endif
   write(nplot,'(''0/'')')

   !--plot the correlation matrix as a shaded contour plot
   call matshd(xig,yig,x,y,ixn,xmin,xmax,ymin,ymax)

   if (allocated(xx)) deallocate(xx)
   if (allocated(xy)) deallocate(xy)

   return
   end subroutine plotit

   subroutine matshd(xig,yig,x,y,ixnow,xmin,xmax,ymin,ymax)
   !--------------------------------------------------------------------
   ! Define regions with similar correlations
   ! and shade them accordingly
   !--------------------------------------------------------------------
   use util ! provides error
   ! externals
   integer::ixnow
   real(kr)::xig(*),yig(*),x(*),y(*)
   real(kr)::xmin,xmax,ymin,ymax,cof
   ! internals
   integer::ipat,inext,jnext,nwlcf,ii,i,jj,j
   integer::il,jl,ilevel,jfirst,jlast,i1,i2,istart,ilim
   integer::iloc,ilast,ifirst,iloop,jpat,none,ntwo,ixmx,jxmx
   real(kr)::xpos,ypos,wa,tlow,thi,cofm,cofa
   character(60)::strn1,strng
   character(1)::qu=''''
   integer,dimension(:),allocatable::ixmip,ixmap,ilcf
   real(kr),parameter::xsize=5.0e0_kr
   real(kr),parameter::ysize=3.38e0_kr
   real(kr),parameter::dx1=-.25e0_kr
   real(kr),parameter::dy1=-.75e0_kr
   real(kr),parameter::dx2=.75e0_kr
   real(kr),parameter::dy2=.60e0_kr
   real(kr),parameter::dx3=.625e0_kr
   real(kr),parameter::dy3=1.75e0_kr
   real(kr),parameter::dy4=1.375e0_kr
   real(kr),parameter::eps=1.e-6_kr
   real(kr),parameter::zero=0
   real(kr),parameter::one=1
   real(kr),parameter::two=2

   !--initialize pattern search.
   ipat=0
   inext=0
   jnext=1
   allocate(ixmip(ixnow))
   allocate(ixmap(ixnow))
   nwlcf=ixnow*ixnow
   allocate(ilcf(nwlcf))
   ii=0
   none=0
   ntwo=0
   cofm=0
   do i=ixmin,ixmax
      ii=ii+1
      jj=0
      do j=ixmin,ixmax
         jj=jj+1
         cof=cf(ixmax*(i-1)+j)
         cofa=abs(cof)-eps
         if (cofa.gt.abs(cofm)) then
            cofm=cof
            ixmx=i
            jxmx=j
         endif
         if (cofa.gt.two) then
            ntwo=ntwo+1
         else if (cofa.gt.one) then
            none=none+1
         endif
         if (cof.gt. one) cof= one
         if (cof.lt.-one) cof=-one
         ilcf(ixnow*(ii-1)+jj)=level(cof,xlev,nlev)
      enddo
   enddo
   if (none.gt.0) then
      write(strn1,'(''processing of mat/mt'',i5,''/'',i3,&
        &'' vs. mat1/mt1'',i5,''/'',i3)') mat,mt,mat1,mt1
      write(strng,'(''largest coefficient='',1p,e13.5&
        &,'' at index'',2i4)') cofm,ixmx,jxmx
      call mess('matshd',strn1,strng)
      write(strng,'(i4,'' coefficients > 1'')') none
      call mess('matshd',strng,'reset and continue.')
   endif
   if (ntwo.gt.0) then
      write(strng,'(i4,'' coefficients > 2'')') ntwo
      call  mess('matshd',strng,'reset and continue')
   endif
   deallocate(cf)

   !--loop over patterns
  200 continue
   ipat=ipat+1
   if (ipat.gt.99999) call error('matshd','ipat gt 99999.',' ')

   !--find next un-allocated square
  210 continue
   inext=inext+1
   if (inext.gt.ixnow) jnext=jnext+1
   if (jnext.gt.ixnow) go to 400
   if (inext.gt.ixnow) inext=1
   if (ilcf(ixnow*(inext-1)+jnext).eq.0) go to 210
   il=inext
   jl=jnext
   ilevel=ilcf(ixnow*(il-1)+jl)
   ilcf(ixnow*(il-1)+jl)=0

   !--find the extent of this pattern in its first row
   ixmip(jl)=il
   ixmap(jl)=il
   jfirst=jl
   jlast=jl
  220 continue
   il=il+1
   if (il.gt.ixnow) go to 230
   if (ilcf(ixnow*(il-1)+jl).ne.ilevel) go to 230
   ilcf(ixnow*(il-1)+jl)=0
   ixmap(jl)=il
   inext=il
   go to 220
  230 continue
   jl=jl+1
   if (jl.gt.ixnow) go to 290
   i1=ixmip(jl-1)
   i2=ixmap(jl-1)
   do 240 i=i1,i2
   if (ilcf(ixnow*(i-1)+jl).ne.ilevel) go to 240

   !--starter square found for this row
   jlast=jl
   ixmip(jl)=i
   ixmap(jl)=i
   istart=i
   ilcf(ixnow*(i-1)+jl)=0
   go to 250
  240 continue

   !--no starter is found.  this is the end of current pattern.
   go to 290

   !--find maximum left extent of pattern in this row
  250 continue
   if (istart.ne.1) then
      ilim=istart-1
      do i=1,ilim
         if (ilcf(ixnow*(istart-i-1)+jl).ne.ilevel) exit
         ixmip(jl)=istart-i
         ilcf(ixnow*(istart-i-1)+jl)=0
      enddo
   endif

   !--now find maximum right extent of pattern in this row
   if (istart.ne.ixnow) then
      ilim=ixnow-istart
      do i=1,ilim
         if (ilcf(ixnow*(istart+i-1)+jl).ne.ilevel) exit
         ixmap(jl)=istart+i
         ilcf(ixnow*(istart+i-1)+jl)=0
      enddo
   endif
   go to 230

   !--determine the boundary curve (xig,yig) for this pattern
  290 continue
   iloc=0
   jl=jfirst-1
   do
      jl=jl+1
      iloc=iloc+1
      if (iloc.gt.nwig) call error('matshd','storage exceeded.',' ')
      ilast=ixmap(jl)
      xig(iloc)=x(ilast+1)
      yig(iloc)=y(jl)
      iloc=iloc+1
      if (iloc.gt.nwig) call error('matshd','storage exceeded.',' ')
      xig(iloc)=x(ilast+1)
      yig(iloc)=y(jl+1)
      if (jl.ge.jlast) exit
   enddo
   jl=jlast+1
   do
      jl=jl-1
      iloc=iloc+1
      if (iloc.gt.nwig) call error('matshd','storage exceeded.',' ')
      ifirst=ixmip(jl)
      xig(iloc)=x(ifirst)
      yig(iloc)=y(jl+1)
      iloc=iloc+1
      if (iloc.gt.nwig) call error('matshd','storage exceeded.',' ')
      xig(iloc)=x(ifirst)
      yig(iloc)=y(jl)
      if (jl.le.jfirst) exit
   enddo

   !--finished with pattern contour logic
   !--ready to draw the pattern
   if (ipat.eq.1) then
      xpos=ysize+dy1
      ypos=dx1
      wa=0
      write(nplot,'(''-1 0 1. 1.'',5f8.3,''/'')')&
        xpos,ypos,xsize,xsize,wa
      write(nplot,'(''/'')')
      write(nplot,'(''/'')')
      write(nplot,'(''4 0 0/'')')
      write(nplot,'(1p,2e12.4,''/'')') xmin,xmax
      write(nplot,'(''/'')')
      write(nplot,'(1p,2e12.4,''/'')') ymin,ymax
      write(nplot,'(''/'')')
   else
      write(nplot,'(''2/'')')
   endif
   write(nplot,'(''/'')')
   call patlev(jpat,ilevel)
   if (jpat.ne.0) ishade=ishade+1
   write(nplot,'(''0 0 0 0 0'',i3,''/'')') jpat
   write(nplot,'(''0'')')
   do i=1,iloc
      write(nplot,'(1p,2e13.4,''/'')') xig(i),yig(i)
   enddo
   write(nplot,'(''/'')')
   go to 200

   !--draw the legend for the contour plot
  400 continue
   xpos=ysize-dx2+xsize+dy2
   ypos=dx2
   wa=90
   write(nplot,'(''-1 0 1. 1.'',5f8.3,''/'')')&
     xpos,ypos,xpos,1.,wa
   write(nplot,'(a,''<c>orrelation <m>atrix'',a,''/'')') qu,qu
   write(nplot,'(''/'')')
   write(nplot,'(''0/'')')
   xpos=xpos+dy3
   ypos=dx3
   write(nplot,'(''-1 0 1. 1.'',5f8.3,''/'')')&
     xpos,ypos,1.75,2.50,wa
   write(nplot,'(''/'')')
   write(nplot,'(''/'')')
   write(nplot,'(''1 0 0/'')')
   write(nplot,'(1p,3e12.4,''/'')') 0.,1.,1.
   write(nplot,'(''/'')')
   write(nplot,'(1p,3e12.4,''/'')') 0.,1.,.2
   write(nplot,'(a,''.'',a,''/'')') qu,qu
   write(nplot,'(''/'')')
   ! positive part
   tlow=0
   do ilevel=1,nlev
      call patlev(jpat,ilevel)
      ! switch orientation for positive part of legend on b-w plots
      if (icolor.eq.0.and.jpat.ne.0) jpat=jpat+10
      write(nplot,'(''0 0 0 0 0'',i3,''/'')') jpat
      write(nplot,'(''0'')')
      thi=xlev(ilevel)
      write(nplot,'(''0. '',f6.3,''/'')') tlow
      write(nplot,'(''1. '',f6.3,''/'')') tlow
      write(nplot,'(''1. '',f6.3,''/'')') thi
      write(nplot,'(''0. '',f6.3,''/'')') thi
      write(nplot,'(''0. '',f6.3,''/'')') tlow
      write(nplot,'(''/'')')
      if (ilevel.lt.nlev) write(nplot,'(''2/'')')
      if (ilevel.lt.nlev) write(nplot,'(''/'')')
      tlow=thi
   enddo
   ypos=ypos+dy4
   write(nplot,'(''-1 0 1. 1.'',5f8.3,''/'')') xpos,ypos,1.75,2.50,wa
   write(nplot,'(''/'')')
   write(nplot,'(''/'')')
   write(nplot,'(''1 0 0/'')')
   write(nplot,'(1p,3e12.4,''/'')') 0.,1.,1.
   write(nplot,'(''/'')')
   write(nplot,'(1p,3e12.4,''/'')') 0.,-1.,-.2
   write(nplot,'(a,''.'',a,''/'')') qu,qu
   write(nplot,'(''/'')')
   ! negative part
   tlow=0
   do iloop=1,nlev
      ilevel=-iloop
      call patlev(jpat,ilevel)
      write(nplot,'(''0 0 0 0 0'',i3,''/'')') jpat
      write(nplot,'(''0'')')
      thi=-xlev(iloop)
      write(nplot,'(''0. '',f6.3,''/'')') tlow
      write(nplot,'(''1. '',f6.3,''/'')') tlow
      write(nplot,'(''1. '',f6.3,''/'')') thi
      write(nplot,'(''0. '',f6.3,''/'')') thi
      write(nplot,'(''0. '',f6.3,''/'')') tlow
      write(nplot,'(''/'')')
      if (iloop.lt.nlev) write(nplot,'(''2/'')')
      if (iloop.lt.nlev) write(nplot,'(''/'')')
      tlow=thi
   enddo

   !--finished
   deallocate(ixmip)
   deallocate(ixmap)
   deallocate(ilcf)
   return
   end subroutine matshd

   integer function level(c,xlev,nlev)
   !--------------------------------------------------------------------
   ! Index the intensity c on scale xlev.
   !--------------------------------------------------------------------
   ! externals
   integer::nlev
   real(kr)::c,xlev(*)
   ! internals
   integer::i,ilev
   real(kr),parameter::zero=0

   do i=1,nlev
      ilev=i
      if (abs(c).lt.xlev(i)) exit
   enddo
   if (c.lt.zero) ilev=-ilev
   level=ilev
   return
   end function level

   subroutine patlev(jpat,ilevel)
   !--------------------------------------------------------------------
   ! Convert the correlation level to a color or b-w shading pattern
   ! following the conventions of viewr.
   !--------------------------------------------------------------------
   ! externals
   integer::jpat,ilevel
   ! internals
   real(kr)::scale,t

   jpat=0
   if (icolor.eq.0) then
      if (ilevel.gt.1) jpat=20-(nlev-ilevel)*2/ndiv
      if (ilevel.lt.-1) jpat=40-(nlev+ilevel)*2/ndiv
   else
      scale=10./(nlev-1)
      if (ilevel.gt.1) then
         t=float(nlev-ilevel)
         jpat=50-nint(scale*t)
      endif
      if (ilevel.lt.-1) then
         t=float(nlev+ilevel)
         jpat=60-nint(scale*t)
      endif
   endif
   return
   end subroutine patlev

   subroutine smilab(iza,mat,mt,izap,mtflg,strng)
   !--------------------------------------------------------------------
   ! Prepare a label for the semi-log plots.
   !--------------------------------------------------------------------
   ! externals
   integer::iza,mat,mt,izap,mtflg
   character strng*(*)
   ! internals
   integer::l1,l2,lstr3
   character str1*14,str2*36,str3*39,str4*16

   l1=6
   if (mt.ge.452.and.mt.le.456.and.mtflg.eq.0) then
       write(str1,'(''[d]n/n'')')
   else if (mt.ge.452.and.mt.le.456.and.mtflg.ne.0) then
       write(str1,'(''    ]n'')')
   else if (mt.eq.251.and.mtflg.ne.0) then
       write(str1,'(''    ]m'')')
   else if (mtflg.ne.0.and.mf35.eq.5) then
       write(str1,'(''Grp-average ]f'')')
       l1=14
   else if (mtflg.eq.0.and.mf35.eq.5) then
       write(str1,'(''[d]f/f'')')
   else if (mtflg.eq.0.and.mt.eq.251) then
          write(str1,'(''[d]m/m'')')
   else if (mtflg.ne.0) then
       write(str1,'(''    ]s'')')
   else
       write(str1,'(''[d]s/s'')')
   endif

   if (mtflg.ne.0.and.mf35.eq.5.and.einc.gt.9.999e4_kr) then
      write(str2,'(''>(<e>#LH>in#HXLX>='',f5.2,&
        &'' <m>e<v>), '')') einc/1.e6
      l2=34
   else if (mtflg.ne.0.and.mf35.eq.5.and.einc.le.9.999e4_kr) then
      write(str2,'(''>(<e>#LH>in#HXLX>= '',1pe9.2,&
        &'' e<v>), '')') einc
      l2=35
   else
      write(str2,'(''> vs. <e> for '')')
      l2=14
   endif
   call matmes(iza,mat,mt,str3,lstr3)
   if (mfflg.ne.-14) then
      strng=str1(1:l1)//str2(1:l2)//str3(1:lstr3)
   elseif (mfflg.eq.-14.and.izap.eq.0) then
      strng=str1(1:l1)//str2(1:l2)//str3(1:lstr3)//', MF40'
   else
      write(str4,'('', izap = '',i7)')izap
      strng=str1(1:l1)//str2(1:l2)//str3(1:lstr3)//str4
   endif

   return
   end subroutine smilab

   subroutine matmes(iza,mat,mt,strng,lstrng)
   !--------------------------------------------------------------------
   ! Generate a reaction name in nuclear physics notation.
   !--------------------------------------------------------------------
   ! externals
   integer::iza,mat,mt,lstrng
   character strng*(*)
   ! internals
   integer::iz,ia,iname,niso,nmat,nnam,inamel,ivl,inamer
   character lnamel*16,lnamer*8,lname*4
   character liso*12,lmat*16,lnam*20

   iz=iza/1000
   ia=iza-iz*1000
   call elem(iz,lname,iname)
   if (ia.gt.0.and.ia.lt.10) then
      write(liso,'(''#EH<'',i1,''#HXEX'')') ia
      niso=10
   else if (ia.ge.10.and.ia.lt.100) then
      write(liso,'(''#EH<'',i2,''#HXEX'')') ia
      niso=11
   else if (ia.ge.100) then
      write(liso,'(''#EH<'',i3,''#HXEX'')') ia
      niso=12
   endif
   if (ia.eq.0) then
      lmat=lname(1:iname)
      nmat=iname
   else
      lmat=liso(1:niso)//lname(1:iname)
      nmat=niso+iname
   endif
   call mtno(mt,lnamel,inamel,ivl,lnamer,inamer)
   if (ivl.le.0) then
      lnam=lnamer
      nnam=inamer
   else if (ivl.lt.10) then
      write(lnam,'(''#L.25H.75<'',i1,''#HXLX<)'')') ivl
      nnam=18
   else if (ivl.lt.100) then
      write(lnam,'(''#L.25H.75<'',i2,''#HXLX<)'')') ivl
      nnam=19
   endif
   strng=lmat(1:nmat)//lnamel(1:inamel)//lnam(1:nnam)
   lstrng=nmat+inamel+nnam
   return
   end subroutine matmes

   subroutine elem(iz,lname,iname)
   !--------------------------------------------------------------------
   ! Convert a 'Z' value into a character element name.
   !--------------------------------------------------------------------
   ! externals
   integer::iz,iname
   character(4)::lname
   ! internals
   character(4),dimension(103)::ename=(/&
     '<h> ','<h>e','<l>i','<b>e','<b> ','<c> ','<n> ','<o> ','<f> ',&
     '<n>e','<n>a','<m>g','<a>l','<s>i','<p> ','<s> ','<c>l','<a>r',&
     '<k> ','<c>a','<s>c','<t>i','<v> ','<c>r','<m>n','<f>e','<c>o',&
     '<n>i','<c>u','<z>n','<g>a','<g>e','<a>s','<s>e','<b>r','<k>r',&
     '<r>b','<s>r','<y> ','<z>r','<n>b','<m>o','<t>c','<r>u','<r>h',&
     '<p>d','<a>g','<c>d','<i>n','<s>n','<s>b','<t>e','<i> ','<x>e',&
     '<c>s','<b>a','<l>a','<c>e','<p>r','<n>d','<p>m','<s>m','<e>u',&
     '<g>d','<t>b','<d>y','<h>o','<e>r','<t>m','<y>b','<l>u','<h>f',&
     '<t>a','<w> ','<r>e','<o>s','<i>r','<p>t','<a>u','<h>g','<t>l',&
     '<p>b','<b>i','<p>o','<a>t','<r>n','<f>r','<r>a','<a>c','<t>h',&
     '<p>a','<u> ','<n>p','<p>u','<a>m','<c>m','<b>k','<c>f','<e>s',&
     '<f>m','<m>d','<n>o','<l>w'/)

   lname=ename(iz)
   iname=4
   if (lname(4:4).eq.' ') iname=3
   return
   end subroutine elem

   subroutine mtno(mt,lnamel,inamel,ivl,lnamer,inamer)
   !--------------------------------------------------------------------
   ! Convert an MT-number to a Hollerith name.
   ! At present, only those MT-numbers associated with File-31
   ! or File-33 data in ENDF/B-5 are recognized.  For other MTs,
   ! the MT-number is used in place of a reaction name.
   !--------------------------------------------------------------------
   use endf ! provides iverf
   ! externals
   integer::mt,inamel,ivl,inamer
   character(16)::lnamel
   character(8)::lnamer
   ! internals
   integer::j,jloc,idone
   integer,dimension(33)::ira1=(/&
     1,2,3,4,16,17,18,22,28,37,51,91,102,103,104,105,106,107,111,&
     207,780,781,452,455,456,25,24,32,33,41,112,115,116/)
   integer,dimension(33)::ira2=(/&
     5,4,7,6,3,3,2,5,3,3,1,6,4,2,2,2,6,4,3,5,5,5,5,5,5,6,6,3,3,4,5,3,3/)
   character(8),dimension(33)::hira=(/&
     'tot.)   ','el.)    ','nonel.) ','inel.)  ','2n)     ',&
     '3n)     ','f)      ','n]a<)   ','np)     ','4n)     ',&
     ')       ','cont.)  ',']g<)    ','p)      ','d)      ',&
     't)      ','<h>e3)  ',']a<)    ','2p      ','<h>e)   ',&
     ']a<0)   ',']a<1)   ',' ]n<)   ',' ]n<)   ',' ]n<)   ',&
     '3n]a<)  ','2n]a<)  ','nd)     ','nt)     ','2np)    ',&
     'p]a<)   ','pd)     ','pt)     '/)
   character(8)::blank='        '
   character(8)::nmea1='(n,     '
   character(8)::nmeh1='(mt     '
   character(8)::nmed1='(total  '
   character(8)::nmeb1='(n,n    '
   character(8)::nmep1='(n,p    '
   character(8)::nmez1='(n,d    '
   character(8)::nmet1='(n,t    '
   character(8)::nmex1='(n,h    '
   character(8)::nmey1='(n,a    '
   character(8)::nmee1='(delayed'
   character(8)::nmef1='(prompt '
   character(8)::nmeg1='(spectr.'
   lnamel=blank

   !--set the left name of this reaction
   if (mt.gt.50.and.mt.lt.92) then
      lnamel=nmeb1
      inamel=4
   else if (mt.ge.600.and.mt.lt.650.and.iverf.gt.5) then
      lnamel=nmep1
      inamel=4
   else if (mt.ge.650.and.mt.lt.700.and.iverf.gt.5) then
      lnamel=nmez1
      inamel=4
   else if (mt.ge.700.and.mt.lt.750.and.iverf.gt.5) then
      lnamel=nmet1
      inamel=4
   else if (mt.ge.750.and.mt.lt.800.and.iverf.gt.5) then
      lnamel=nmex1
      inamel=4
   else if (mt.ge.800.and.mt.lt.850.and.iverf.gt.5) then
      lnamel=nmey1
      inamel=4
   else if (mt.eq.452) then
      lnamel=nmed1
      inamel=6
   else if (mt.eq.455) then
      lnamel=nmee1
      inamel=8
   else if (mt.eq.456) then
      lnamel=nmef1
      inamel=7
   else if (mt.eq.251) then
      lnamel=nmeh1
      inamel=3
   else if (mt.eq.261) then
      lnamel=nmeg1
      inamel=8
   else
      lnamel=nmea1
      inamel=3
   endif
   jloc=11

   !--set the discrete level number
   ivl=-1
   if (mt.gt.50.and.mt.lt.91) ivl=mt-50
   if (mt.ge.600.and.mt.lt.650.and.iverf.gt.5) ivl=mt-600
   if (mt.ge.650.and.mt.lt.700.and.iverf.gt.5) ivl=mt-650
   if (mt.ge.700.and.mt.lt.750.and.iverf.gt.5) ivl=mt-700
   if (mt.ge.750.and.mt.lt.800.and.iverf.gt.5) ivl=mt-750
   if (mt.ge.800.and.mt.lt.850.and.iverf.gt.5) ivl=mt-800

   !--set the right name of this reaction
   j=0
   idone=0
   do while (j.lt.33.and.idone.eq.0)
      j=j+1
      jloc=j
      if (ira1(j).eq.mt) idone=1
      if (ivl.ge.0.and.ivl.le.99.and.ira1(j).eq.51) idone=1
   enddo
   if (idone.eq.0) then
      write(nmeh1(4:6),'(i3)') mt
      lnamel=nmeh1
      inamel=6
      jloc=11
   endif
   lnamer=hira(jloc)
   inamer=ira2(jloc)
   return
   end subroutine mtno

   subroutine copyst(nin,nout)
   !--------------------------------------------------------------------
   ! Copy nin to nout.
   !--------------------------------------------------------------------
   use util ! provides repoz
   ! externals
   integer::nin,nout
   ! internals
   integer::i
   character(4)::a(20)

   call repoz(nin)
  110 continue
   read(nin,'(20a4)',end=130) (a(i),i=1,20)
   write(nout,'(20a4)') (a(i),i=1,20)
   go to 110
  130 continue
   return
   end subroutine copyst

   subroutine finds(mat,mf,mt,ntape)
   !--------------------------------------------------------------------
   ! Search through the errorr tape ntape for the desired section.
   ! If MF=0, find the first occurrence of that MAT.
   ! If MT=0, find the first section of that MAT with the file MF.
   !--------------------------------------------------------------------
   use util ! provides error,repoz,skiprz
   ! externals
   integer::mat,mf,mt,ntape
   ! internals
   integer::nt,mat1,mf1,mt1,irew,i,math,mfh,mth
   character(60)::strng

   !--initialize
   nt=iabs(ntape)
   mat1=-2
   mf1=-2
   mt1=-2
   irew=0
   i=0

   !--read through ntape
  110 continue
   if (ntape.ge.0) then
      read(ntape,'(66x,i4,i2,i3)') math,mfh,mth
   else
      read(nt) math,mfh,mth
   endif
   if (math.eq.-1) go to 160
   if (math.eq.0.or.mfh.eq.0.or.mth.eq.0) go to 110
   if (mat1.ne.-2) go to 140
   mat1=math
   mf1=mfh
   mt1=mth
  140 continue
   i=i+1
   if (math.ne.mat) go to 150
   if (mf.eq.0) go to 210
   if (mfh.ne.mf) go to 150
   if (mth.eq.mt.or.mt.eq.0) go to 210
  150 continue
   if (irew.eq.0) go to 110
   if (mat1.ne.math.or.mf1.ne.mfh.or.mt1.ne.mth) go to 110
  160 continue
   if (irew.gt.0) then
      write(strng,&
        '(''mat '',i4,'' mf '',i2,'' mt '',i3,'' not on tape '',i3)')&
        mat,mf,mt,ntape
      call error('finds',strng,' ')
   endif
   irew=irew+1
   call repoz(ntape)
   go to 110

   !--finished
  210 continue
   call skiprz(ntape,-1)
   ! check to see we are at beginning of section
   if (i.gt.1) return
  220 continue
   if (ntape.ge.0) then
      read(ntape,'(66x,i4,i2,i3)') math,mfh,mth
   else
      read(nt) math,mfh,mth
   endif
   if (math.eq.1) return
   if (mat.ne.math) go to 250
   if (mf.eq.0.or.mfh.eq.mf) go to 260
  250 continue
   call skiprz(ntape,-2)
   go to 220
  260 continue
   if (mt.eq.0.or.mth.eq.mt) go to 270
   go to 220
  270 continue
   i=i+1
   go to 210
   end subroutine finds

   subroutine press(nscr,nout,xa,nr,nc)
   !--------------------------------------------------------------------
   ! Convert the data matrix xa(i,j) to the compressed "boxer" format
   ! and write the result to unit nout.
   !
   ! Definitions:
   !
   !   nrow      number of rows in matrix xa,
   !   ncol      number of columns in xa,
   !                =1 if xa is a vector.
   !                =0 if xa is a symmetric matrix, represented in the
   !                   *press* format by just the upper right triangle
   !                   (i.e., data begins at the diagonal).
   !   i         row index of matrix xa(i,j), normally the energy group
   !                of the reaction (mat,mt).
   !   j         column index of xa(i,j), normally the energy group of
   !                the reaction (mat1,mt1).
   !   xval(iv)  data value array in the *press* format (iv=1,nval).
   !   nvmax     maximum allowable value of nval.
   !   icon(ic)  control parameter array in the *press* format
   !                (ic=1,ncon).
   !   ncmax     maximum allowable value of ncon.
   !   nvf       format-type for xval(iv).
   !   ncf       format-type for icon(ic).
   !   itype     data type
   !                =0 for energy-group boundaries
   !                =1 for cross sections
   !                =2 for relative standard deviations
   !                =3 for relative covariance matrix
   !                =4 for correlation matrix
   !
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso,nsyse
   use util ! provides error,repoz,mess,sigfig
   ! externals
   integer::nscr,nout,nr,nc
   real(kr)::xa(nr,nc)
   ! internals
   integer::ndimx,ndig,lvmax,lcmax,i,istart
   integer::lv,lc,ivopt,icopt,iv,ic,nsym,ig,j
   integer::ii,jj,jlow,im1,im1ind,iind,icd
   integer::leng,nval,ncon,nrowm
   real(kr)::xvtest,crit,test,xctest
   real(kr)::xval(880)
   integer::icon(900)
   character(60)::strng
   character(12)::ivft,icft
   integer,parameter::nvmax=880
   integer,parameter::ncmax=900
   real(kr),parameter::eps=1.e-10_kr
   real(kr),parameter::zero=0

   ndimx=max0(nr,nc)+1
   if (nrow.gt.ndimx.or.ncol.gt.ndimx)&
     call error('press','storage exceeded.',' ')

   !--output
   if (nout.eq.0) go to 510
   ! set formats
   call setfor(ivft,icft,nvf,ncf)
   ! initialize
   ndig=nvf-6
   if (nvf.lt.9) ndig=nvf-3
   lvmax=10**(ncf-1)-1
   lcmax=10**ncf-1
   i=0
   istart=1
   call repoz(-nscr1)
  480 continue
   lv=0
   lc=0
   ivopt=1
   icopt=1
   iv=0
   ic=0
   nsym=0
   if (ncol.eq.0) nsym=1
   if (ncol.eq.0) ncol=nrow
   if (nrow.eq.ncol) then
      read(nscr1) ig,(scr(j),j=1,ncol)
      do j=1,ncol
         xa(2,j)=sigfig(scr(j),ndig,0)
      enddo
      xvtest=xa(2,1)
      if (nsym.eq.1) xvtest=xa(2,i+1)
   else
      do ii=1,nrow
         do jj=1,ncol
            xa(ii,jj)=sigfig(xa(ii,jj),ndig,0)
         enddo
      enddo
      xvtest=xa(istart,1)
   endif

   !--build xval and icon from a matrix xa(i,j)
   ! begin loop over rows
  470 continue
   i=i+1
   jlow=1
   im1=i-1
   im1ind=im1
   iind=i
   if (ncol.ne.nrow) go to 250
   im1ind=1
   iind=2
   if (nsym.eq.1) jlow=i
   if (i.eq.istart) go to 250
   do j=1,ncol
      xa(1,j)=xa(2,j)
   enddo
   do while (ig.lt.i)
      read(nscr1) ig,(scr(j),j=1,ncol)
   enddo
   do j=1,ncol
      xa(iind,j)=sigfig(scr(j),ndig,0)
   enddo
   ! check for symmetry in the matrix
   if (nsym.ne.1) go to 250
   ! begin loop over columns
   if (xa(im1ind,i).eq.zero) go to 250
   crit=(xa(im1ind,i)-xa(iind,im1))/xa(im1ind,i)
   test=1
   test=test/1000
   if (abs(crit).gt.test) then
      write(strng,&
        '(''i'',i3,'' j'',i3,'' xa(i,j)='',1p,e12.4,'' xa(j,i)='',&
       &1p,e12.4)') ii,im1,xa(im1ind,ii),xa(iind,im1)
      call error('press','matrix not symmetric',strng)
   endif
   test=1
   test=test/1000000
   if (abs(crit).lt.test) go to 250
   write(strng,&
     '(''i'',i3,'' j'',i3,'' xa(i,j)='',1p,e12.4,'' xa(j,i)='',&
     &1p,e12.4)') ii,im1,xa(im1ind,ii),xa(iind,im1)
   call mess('press','matrix not symmetric',strng)
  250 continue
   do 210 j=jlow,ncol

   !--begin pattern-search logic
  310 continue
   if (ivopt.ne.0) then
      test=eps*abs(xvtest)
      if (abs(xa(iind,j)-xvtest).le.test.and.lv.ne.lvmax) then
         ! passes repeated-value test
         lv=lv+1
         icd=0
      else
         ! fails test
         ivopt=0
      endif
   endif
   if (icopt.ne.0) then
      ! test carry-down option
      xctest=0
      if (i.gt.istart) xctest=xa(im1ind,j)
      test=eps*abs(xctest)
      if (abs(xa(iind,j)-xctest).le.test.and.lc.ne.lcmax) then
         ! passes carry-down test
         lc=lc+1
         icd=1
      else
         ! fails test
         icopt=0
      endif
   endif
   if (icopt.ne.0.or.ivopt.ne.0) go to 210

   !--have found end of data pattern.  store data and re-initialize
   ic=ic+1
   if (icd.ne.0) then
      ! carry-down option
      icon(ic)=lc
   else
      ! repeated-value option
      icon(ic)=-lv
      iv=iv+1
      xval(iv)=xvtest
   endif
   lv=0
   lc=0
   ivopt=1
   icopt=1
   xvtest=xa(iind,j)
   go to 310
  210 continue
   ! if container is nearly full, discontinue pattern searches,
   ! write out the data processed up to this point, and then
   ! begin building a new xval and icon array (ie., a new page).
   leng=ncol+1
   if (nsym.eq.1) leng=ncol+1-i
   if (iv+leng.gt.nvmax.or.ic+leng.gt.ncmax) go to 440
   if (i.lt.nrow) go to 470
   ! store last pattern
  440 ic=ic+1
   if (icd.ne.0) then
      icon(ic)=lc
   else
      icon(ic)=-lv
      iv=iv+1
      xval(iv)=xvtest
   endif
   nval=iv
   ncon=ic
   nrowm=nrow-i
   if (nsym.eq.1) ncol=0

   !--finished building xval and icon.  write results
   if (istart.eq.1) write(nout,&
     '(i1,1x,a12,1x,a21,2(i5,i4),2(i4,i3),3i4)')&
     itype,hlibid,hdescr,mat,mt,mat1,&
     mt1,nval,nvf,ncon,ncf,nrowm,nrow,ncol
   if (istart.gt.1) write(nout,&
     '(i1,1x,34(''-''),2(i5,i4),2(i4,i3),3i4)')&
     itype,mat,mt,mat1,mt1,nval,nvf,ncon,ncf,nrowm,nrow,ncol
   if (nval.gt.0) write(nout,ivft) (xval(iv),iv=1,nval)
   if (ncon.gt.0) write(nout,icft) (icon(ic),ic=1,ncon)
   if (nrowm.eq.0) go to 500
   ! write new page
   istart=i+1
   go to 480
  500 continue
   if (itype.gt.2) then
      write(nsyso,'(/3x,''matrix data written.'')')
   endif
  510 continue
   return
   end subroutine press

   subroutine setfor(ivft,icft,nvf,ncf)
   !--------------------------------------------------------------------
   ! Set formats for input and output.
   !--------------------------------------------------------------------
   use util ! provides error
   ! externals
   integer::nvf,ncf
   character(12)::ivft,icft
   ! internals
   character(60)::strng
   character(12),dimension(14)::ift=(/&
     '(80i1)      ','(40i2)      ','(26i3)      ','(20i4)      ',&
     '(16i5)      ','(13i6)      ','(11f7.4)    ','(10f8.5)    ',&
     '(1p8e9.2)   ','(1p8e10.3)  ','(1p7e11.4)  ','(1p6e12.5)  ',&
     '(1p6e13.6)  ','(1p5e14.7)  '/)

   if (nvf.lt.7.or.nvf.gt.14.or.ncf.lt.1.or.ncf.gt.6) then
      write(strng,&
        '(''nvf (='',i3,'') or ncf (='',i3,'' is illegal.'')')&
        nvf,ncf
      call error('setfor',strng,' ')
   endif

   !--set formats
   ivft=ift(nvf)
   icft=ift(ncf)
   return
   end subroutine setfor

end module covm


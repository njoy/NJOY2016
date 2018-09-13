module gaminm
   ! provides subroutine gaminr for NJOY2016
   use locale
   implicit none
   private
   public gaminr

   ! global variables
   integer::igg,ngg
   integer,parameter::ngmax=400
   real(kr)::egg(ngmax)
   integer::iwt
   integer,parameter::iwmax=200
   real(kr)::wght(iwmax)
   integer::matb,lord
   integer::nendf,npend
   integer::ngam1,ngam2,ntw
   character(17)::title(17)
   integer::matd,mfd,mtd
   integer::ig2pp
   integer::iprint

contains

   subroutine gaminr
   !------------------------------------------------------------------
   !
   ! compute multigroup photon cross sections
   !
   ! Produce multigroup photon interaction cross sections
   ! and heating kerma factors using ENDF cross sections
   ! and coherent and incoherent form factors.  Initial energy
   ! quadrature techiques are identical to those used in groupr.
   ! Secondary energy-angle quadrature is performed using Gaussian
   ! integration.
   !
   !---input specifications (free format)---------------------------
   !
   ! card1
   !    nendf   unit for endf tape
   !    npend   unit for pendf tape
   !    ngam1   unit for input ngam tape (default=0)
   !    ngam2   unit for output ngam tape (default=0)
   ! card2
   !    matb    material to be processed
   !            input materials in ascending order
   !    igg     gamma group structure option
   !    iwt     weight function option
   !    lord    legendre order
   !    iprint  print option (0/1=minimum/maximum) (default=1)
   ! card3
   !    title   run label up to 80 characters (delimited by ',
   !            ended with /)
   ! card4      (igg=1 only)
   !    ngg     number of groups
   !    egg     ngg+1 group bounds (ev)
   ! card5      (iwt=1 only)
   !    wght    weight function as tab1 record
   ! card6
   !    mfd     file to be processed
   !    mtd     section to be processed
   !    mtname  description of section to be processed
   !            repeat for all reactions desired
   !            mfd=0/ terminates this material
   !            mfd=-1/ is a flag to process all sections present
   !            for this material  (termination is automatic)
   ! card7
   !    matd    next mat number to be processed
   !            terminate gaminr run with matd=0.
   !
   !---options for input variables----------------------------------
   !
   !        igg     meaning
   !        ---     -------
   !         0      none
   !         1      arbitrary structure (read in)
   !         2      csewg 94-group structure
   !         3      lanl 12-group structure
   !         4      steiner 21-group gamma-ray structure
   !         5      straker 22-group structure
   !         6      lanl 48-group structure
   !         7      lanl 24-group structure
   !         8      vitamin-c 36-group structure
   !         9      vitamin-e 38-group structure
   !         10      vitamin-j 42-group structure
   !
   !        iwt     meaning
   !        ---     -------
   !         1      read in
   !         2      constant
   !         3      1/e + rolloffs
   !
   !------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso,nsyse
   use endf   ! provides endf routines and variables
   use util   ! provides timer,repoz,skiprz,closz
   use physics ! provides electron mass
   integer::nw,nb,nwds,nggp1,np,loc,nleft,icnt,mflg,nl
   integer::idis,ng2,iglo,nq,ig1,ig,ig2lo,it,iz,il
   integer::idone,igzero,j,ibase,i,lim,l
   integer::ntoth,ii,jj
   real(kr)::time,zref,znow,za,awr,e,en,elo,ehi,enext
   real(kr)::thresh
   character(4)::mtname(15)
   integer::ng2s(ngmax+1),ig2s(ngmax+1)
   real(kr)::z(10,10)
   character(66)::text
   character(4)::tt(17)
   real(kr)::rt(17)
   equivalence(tt(1),rt(1))
   integer,dimension(9),parameter::mflst=(/&
     23,23,26,23,26,23,26,23,23/)
   integer,dimension(9),parameter::mtlst=(/&
     501,502,502,504,504,516,516,602,621/)
   integer,dimension(9),parameter::mtlst6=(/&
     501,502,502,504,504,516,516,522,525/)
   character(4),dimension(9)::nmlst=(/&
     'totl','coht','coht','inch','inch','pair','pair','abst','heat'/)
   real(kr),dimension(:),allocatable::scr
   real(kr),dimension(:),allocatable::sdat
   real(kr),dimension(:),allocatable::pff
   real(kr),dimension(:),allocatable::toth
   real(kr),dimension(:,:),allocatable::ff
   real(kr),dimension(:,:,:),allocatable::akn
   real(kr),dimension(:,:,:),allocatable::ans
   real(kr),parameter::ekn=12.4e3_kr
   real(kr),parameter::emax=1.e12_kr
   real(kr),parameter::zero=0
   integer,parameter::ncnt=9
   integer,parameter::nwpff=5000
   integer::nz=1

   !--write heading and read user input
   call timer(time)
   write(nsyso,'(/&
     &'' gaminr...produce photon interaction cross sections'',&
     &18x,f8.1,''s'')') time
   write(nsyse,'(/'' gaminr...'',59x,f8.1,''s'')') time
   call ruing

   !--allocate storage arrays
   nw=npage+50
   allocate(scr(nw))
   allocate(sdat(nw))
   allocate(pff(nwpff))
   ntoth=ngg*2
   allocate(toth(ntoth))
   zref=101
   znow=0
   allocate(akn(lord+1,ngg+3,ngg))
   call repoz(nendf)

   !--determine endf format being used
   call repoz(nendf)
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
   call repoz(nendf)
   call repoz(npend)
   call tpidio(nendf,0,0,scr,nb,nw)
   call tpidio(npend,0,0,scr,nb,nwds)

   !--locate position for new material on old gamout tape.
   call repoz(ngam1)
   call repoz(ngam2)
   nsh=0
   if (ngam1.eq.0) call tpidio(0,ngam2,0,scr,nb,nwds)
   if (ngam1.eq.0.or.ngam2.eq.0) go to 150
   call tpidio(ngam1,ngam2,0,scr,nb,nwds)
  110 continue
   call contio(ngam1,0,0,scr,nb,nw)
   if (math.eq.-1) go to 140
   if (math.eq.1) go to 110
   if (math.eq.matb) go to 130
   if (math.gt.matb) go to 140
   call contio(0,ngam2,0,scr,nb,nw)
   call tomend(ngam1,ngam2,0,scr)
   go to 110
  130 continue
   call tomend(ngam1,0,0,scr)
   go to 150
  140 continue
   call skiprz(ngam1,-1)
  150 continue

   !--search for desired material on pendf tape
  160 continue
   call findf(matb,1,451,npend)
   call contio(npend,0,0,scr,nb,nw)
   za=c1h
   awr=c2h
   if (iverf.ge.5) call contio(npend,0,0,scr,nb,nw)
   if (iverf.ge.6) call contio(npend,0,0,scr,nb,nw)
   call hdatio(npend,0,0,scr,nb,nw)
   do i=1,17
      rt(i)=scr(6+i)
   enddo
   write(nsyso,'(/&
     &'' processing mat '',i4/1x,66(''-'')/1x,17a4)')&
     matb,(tt(i),i=1,17)

   !--write head record for this material on gamout tape.
   if (ngam2.ne.0) then
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
      call contio(0,ngam2,0,scr,nb,nw)
      scr(1)=0
      scr(2)=0
      scr(3)=ngg
      scr(4)=0
      scr(5)=0
      scr(6)=0
      nw=6
      do i=1,ntw
         nw=nw+1
         tt(i)=title(i)
         scr(nw)=rt(i)
      enddo
      nw=nw+1
      scr(nw)=emax
      nggp1=ngg+1
      do i=1,nggp1
         nw=nw+1
         if (egg(i).lt.epair) ig2pp=i
         scr(nw)=egg(i)
      enddo
      nw=nw+1
      scr(nw)=0
      np=nw-6
      scr(5)=np
      loc=1
      nw=np
      if (nw.gt.npage) nw=npage
      nw=nw+6
      call listio(0,ngam2,0,scr(loc),nb,nw)
      nleft=np-nw+6
      do while (nb.ne.0)
         loc=loc+nw
         nw=nleft
         if (nw.gt.npage) nw=npage
         call moreio(0,ngam2,0,scr(loc),nb,nw)
      enddo
      call afend(ngam2,0)
   endif

   !--main entry point for loop over reactions
   !--initialize and write appropriate heading.
   !--write head record on output tape.
   matd=matb
   do i=1,ntoth
      toth(i)=0
   enddo
   icnt=0
   mflg=0
  220 continue
   if (mflg.eq.-1) go to 250
   text=' '
   read(nsysi,*) mfd,mtd,text
   mflg=mfd
   if (mflg.eq.-1) go to 250
   if (mfd.eq.0) go to 420
   read(text,'(15a4)') (mtname(i),i=1,15)
   go to 270
  250 continue
   icnt=icnt+1
   if (icnt.gt.ncnt) go to 420
   text=' '
   read(text,'(15a4)') (mtname(i),i=1,15)
   mfd=mflst(icnt)
   mtd=mtlst(icnt)
   if (iverf.ge.6) mtd=mtlst6(icnt)
   read(nmlst(icnt),'(a4)') mtname(1)
  270 continue
   call timer(time)

   ! special branch for total heating assigned to mt=621.
   ! for version 6 format, heating is mt=525
   if (mtd.eq.621) go to 400
   if (mtd.eq.525) go to 400

   ! branch for normal reactions
   nl=1
   if (mfd.eq.26) nl=lord+1
   if (mtd.eq.516) nl=1
   nz=1
   if (allocated(ans)) deallocate(ans)
   allocate(ans(nl,nz,ngg+3))
   allocate(ff(nl,ngg+3))
   e=0
   call gtsig(e,thresh,idis,z,sdat)
   call gtflx(e,en,idis,z,nl,nz)
   call gtff(e,en,idis,ff,nl,ng2,iglo,nq,pff,nwpff)
   ig1=0
   ig=0
   if (mfd.eq.26.and.mtd.eq.504) then
      znow=nint(za/1000)
      if (znow.le.zref) then
         zref=znow
         znow=-znow
      endif
   endif

   !--loop over initial energy groups
  300 continue
   ig=ig+1
   elo=egg(ig)
   ehi=egg(ig+1)
   if (znow.gt.zero.and.elo.ge.znow*ekn) go to 325
   ig2lo=0
   ng2=2
   if (ehi.le.thresh) go to 380
   enext=ehi
   do it=1,ngg+3
      do iz=1,nz
         do il=1,nl
            ans(il,iz,it)=0
         enddo
      enddo
   enddo
   idone=0
   do while (idone.eq.0)
      call gpanel(elo,enext,ans,ff,nl,nz,ng2,ig2lo,sdat,pff,nwpff)
      if (enext.eq.ehi) then
         idone=1
      else
         elo=enext
         enext=ehi
      endif
   enddo
   go to 330
  325 continue
   do i=1,3+ngg
      do j=1,nl
         ans(j,1,i)=akn(j,i,ig)
         if (i.gt.1) ans(j,1,i)=ans(j,1,i)*znow
      enddo
   enddo
   ng2=ng2s(ig)
   ig2lo=ig2s(ig)
  330 continue

   !--print results for this initial energy group
   if (ig1.eq.0) then
      ig1=1
      write(nsyso,'(/'' group constants'',53x,f8.1,''s'')') time
      write(nsyso,'('' for mf'',i2,'' and mt'',i3,1x,15a4)')&
        mfd,mtd,(mtname(i),i=1,15)
      math=matb
      mfh=mfd
      mth=mtd
      scr(1)=za
      scr(2)=awr
      scr(3)=nl
      scr(4)=nz
      scr(5)=0
      scr(6)=ngg
      call contio(0,ngam2,0,scr,nb,nwds)
      call dspla(-1,ans,nl,nz,ng2,ig2lo,igzero)
   endif
   if (znow.lt.zero) then
      do i=1,3+ngg
         do j=1,nl
            akn(j,i,ig)=ans(j,1,i)
            if (i.gt.1) akn(j,i,ig)=akn(j,i,ig)/abs(znow)
         enddo
      enddo
      ng2s(ig)=ng2
      ig2s(ig)=ig2lo
   endif
   call dspla(ig,ans,nl,nz,ng2,ig2lo,igzero)

   !--accumulate total heating.
   if (ng2.ne.2.and.mtd.ne.502) then
      l=1+2*(ig-1)
      toth(l)=ans(1,1,1)
      toth(l+1)=toth(l+1)+ans(1,1,ng2)*ans(1,1,1)
      ng2=ng2-1
      if (mtd.eq.504) ng2=ng2-1
   endif

   !--write results on output tape
   nw=nl*nz*ng2
   lim=nw
   if (ngam2.ne.0.and.(igzero.ne.0.or.ig.eq.ngg)) then
      math=matb
      mfh=mfd
      mth=mtd
      scr(1)=0
      scr(2)=0
      scr(3)=ng2
      scr(4)=ig2lo
      scr(5)=nw
      scr(6)=ig
      j=6
      ibase=6
      do i=1,lim
         j=j+1
         jj=(i-1)/nl
         ii=i-nl*jj
         jj=jj+1
         scr(j)=ans(ii,1,jj)
         if (j.ge.(npage+ibase).or.i.eq.lim) then
            if (ibase.ne.0) then
               call listio(0,ngam2,0,scr,nb,j)
               ibase=0
               j=0
            else
               call moreio(0,ngam2,0,scr,nb,j)
               j=0
            endif
         endif
      enddo
      if (ig.eq.ngg) call asend(ngam2,0)
   endif

   !--continue the loop over groups
  380 continue
   if (ig.lt.ngg) go to 300
   if (ig1.eq.0) write(nsyso,'(/&
     &'' threshold above highest energy bound for mt='',i3)') mtd
   znow=0

   !--loop over other desired reactions
   deallocate(ff)
   deallocate(ans)
   go to 220

   !--edit out total heating.
  400 continue
   ig=0
   ig1=0
   nl=1
   do while (ig.lt.ngg)
      ig=ig+1
      ng2=2
      ig2lo=1
      allocate(ans(1,1,2))
      l=1+2*(ig-1)
      ans(1,1,1)=toth(l)
      ans(1,1,2)=toth(l+1)
      if (ig1.eq.0) then
         ig1=1
         write(nsyso,'(/'' group constants'',53x,f8.1,''s'')') time
         write(nsyso,'('' for mf'',i2,'' and mt'',i3,1x,15a4)')&
           mfd,mtd,(mtname(i),i=1,15)
         math=matb
         mfh=mfd
         mth=mtd
         scr(1)=za
         scr(2)=awr
         scr(3)=nl
         scr(4)=nz
         scr(5)=0
         scr(6)=ngg
         call contio(0,ngam2,0,scr,nb,nwds)
         call dspla(-1,ans,nl,nz,ng2,ig2lo,igzero)
      endif
      call dspla(ig,ans,nl,nz,ng2,ig2lo,igzero)
      nw=nl*nz*ng2
      lim=nw
      if (ngam2.ne.0.and.(igzero.ne.0.or.ig.eq.ngg)) then
         math=matb
         mfh=mfd
         mth=mtd
         scr(1)=0
         scr(2)=0
         scr(3)=ng2
         scr(4)=ig2lo
         scr(5)=nw
         scr(6)=ig
         j=6
         ibase=6
         do i=1,lim
            j=j+1
            jj=(i-1)/nl
            ii=i-nl*jj
            jj=jj+1
            scr(j)=ans(ii,1,jj)
            if (j.ge.(npage+ibase).or.i.eq.lim) then
               if (ibase.ne.0) then
                  call listio(0,ngam2,0,scr,nb,j)
                  ibase=0
                  j=0
               else
                  call moreio(0,ngam2,0,scr,nb,j)
                  j=0
               endif
            endif
         enddo
         if (ig.eq.ngg) call asend(ngam2,0)
      endif
      deallocate(ans)
   enddo
   go to 220

   !--loop over desired materials.
  420 continue
   call amend(ngam2,0)
   matb=0
   read(nsysi,*) matb
   if (matb.eq.0) go to 500
   if (ngam1.eq.0) go to 160
   go to 110

   !--gaminr is complete.
  500 continue
   if (ngam1.eq.0) call atend(ngam2,0)
   if (ngam1.ne.0) call totend(ngam1,ngam2,0,scr)
   call repoz(ngam1)
   call repoz(ngam2)
   call repoz(nendf)
   call repoz(npend)
   call closz(ngam1)
   call closz(ngam2)
   call closz(nendf)
   call closz(npend)
   call timer(time)
   write(nsyso,'(69x,f8.1,''s''/1x,77(''*''))') time
   return
   end subroutine gaminr

   subroutine ruing
   !------------------------------------------------------------------
   ! Read user input for group averaging.
   !------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use util   ! provides openz
   ! internals
   integer::i
   character(66)::text

   !--read and display user input
   ngam1=0
   ngam2=0
   read(nsysi,*) nendf,npend,ngam1,ngam2
   call openz(nendf,0)
   call openz(npend,0)
   call openz(ngam1,0)
   call openz(ngam2,1)
   iprint=1
   read(nsysi,*) matb,igg,iwt,lord,iprint
   read(nsysi,*) text
   read(text,'(16a4,a2)') (title(i),i=1,17)
   write(nsyso,'(/&
     &'' unit for endf tape ................... '',i10/&
     &'' unit for pendf tape .................. '',i10/&
     &'' unit for input gamout tape ........... '',i10/&
     &'' unit for output gamout tape .......... '',i10)')&
     nendf,npend,ngam1,ngam2
   write(nsyso,'(&
     &'' mat to be processed .................. '',i10/&
     &'' gamma group option ................... '',i10/&
     &'' weight function option ............... '',i10/&
     &'' legendre order ....................... '',i10/&
     &'' print option (0 min, 1 max) .......... '',i10)')&
     matb,igg,iwt,lord,iprint
   write(nsyso,'(/'' run title''/1x,38(''-'')/&
     &6x,16a4,a2)') (title(i),i=1,17)
   write(nsyso,'('' '')')
   if (ngam2.gt.0) then
      ntw=1
      title(1)=''
   endif
   call genggp
   call gnwtf
   return
   end subroutine ruing

   subroutine genggp
   !------------------------------------------------------------------
   ! Generate requested gamma group structure or read in from
   ! the system input file in the form of an ENDF/B list record.
   !
   !    igg     meaning
   !    ---     --------------------------------------
   !     0      none
   !     1      arbitrary structure (read in)
   !     2      csewg 94-group structure
   !     3      lanl 12-group structure
   !     4      steiner 21-group gamma-ray structure (ornl-tm-2564)
   !     5      straker 22-group structure
   !     6      lanl 48-group structure
   !     7      lanl 24-group structure
   !     8      vitamin-c 36-group structure
   !     9      vitamin-e 38-group structure (r. roussin, feb 86)
   !    10      vitamin-j 42-group structure
   !
   !------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use util   ! provides error
   ! internals
   integer::ig,ngm,i,ngp
   real(kr),dimension(95),parameter::eg2=(/&
     .005e0_kr,.01e0_kr,.015e0_kr,.02e0_kr,.03e0_kr,.035e0_kr,&
     .04e0_kr,.045e0_kr,.055e0_kr,.06e0_kr,.065e0_kr,.075e0_kr,&
     .08e0_kr,.09e0_kr,.1e0_kr,.12e0_kr,.14e0_kr,.15e0_kr,.16e0_kr,&
     .19e0_kr,.22e0_kr,.26e0_kr,.3e0_kr,.325e0_kr,.35e0_kr,&
     .375e0_kr,.4e0_kr,.425e0_kr,.45e0_kr,.5e0_kr,.525e0_kr,&
     .55e0_kr,.575e0_kr,.6e0_kr,.625e0_kr,.65e0_kr,.675e0_kr,&
     .7e0_kr,.75e0_kr,.8e0_kr,.825e0_kr,.865e0_kr,.9e0_kr,1.e0_kr,&
     1.125e0_kr,1.2e0_kr,1.25e0_kr,1.33e0_kr,1.42e0_kr,1.5e0_kr,&
     1.6e0_kr,1.66e0_kr,1.75e0_kr,1.875e0_kr,2.e0_kr,2.166e0_kr,&
     2.333e0_kr,2.5e0_kr,2.666e0_kr,2.833e0_kr,3.e0_kr,3.166e0_kr,&
     3.333e0_kr,3.5e0_kr,3.65e0_kr,3.8e0_kr,3.9e0_kr,4.e0_kr,&
     4.2e0_kr,4.4e0_kr,4.5e0_kr,4.7e0_kr,5.e0_kr,5.2e0_kr,5.4e0_kr,&
     5.5e0_kr,5.75e0_kr,6.e0_kr,6.25e0_kr,6.5e0_kr,6.75e0_kr,&
     7.e0_kr,7.25e0_kr,7.5e0_kr,7.75e0_kr,8.e0_kr,8.5e0_kr,9.e0_kr,&
     9.5e0_kr,10.e0_kr,10.6e0_kr,11.e0_kr,12.e0_kr,14.e0_kr,20.e0_kr/)
   real(kr),dimension(13),parameter::eg3=(/&
     .01e0_kr,.10e0_kr,.50e0_kr,1.0e0_kr,2.0e0_kr,3.0e0_kr,4.0e0_kr,&
     5.0e0_kr,6.0e0_kr,7.0e0_kr,8.0e0_kr,9.0e0_kr,20.0e0_kr/)
   real(kr),dimension(22),parameter::eg4=(/&
     .01e0_kr,.1e0_kr,.2e0_kr,.4e0_kr,1.e0_kr,1.5e0_kr,2.e0_kr,&
     2.5e0_kr,3.e0_kr,3.5e0_kr,4.e0_kr,4.5e0_kr,5.e0_kr,5.5e0_kr,&
     6.e0_kr,6.5e0_kr,7.e0_kr,7.5e0_kr,8.e0_kr,10.e0_kr,12.e0_kr,&
     14.e0_kr/)
   real(kr),dimension(23),parameter::eg5=(/&
     .01e0_kr,.03e0_kr,.06e0_kr,.10e0_kr,.15e0_kr,.30e0_kr,.45e0_kr,&
     .60e0_kr,.80e0_kr,1.0e0_kr,1.33e0_kr,1.66e0_kr,2.0e0_kr,&
     2.5e0_kr,3.0e0_kr,3.5e0_kr,4.0e0_kr,5.0e0_kr,6.0e0_kr,7.0e0_kr,&
     8.0e0_kr,10.0e0_kr,14.0e0_kr/)
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
     5.25e5_kr,7.5e5_kr,1.e6_kr,1.33e6_kr,1.66e6_kr,2.e6_kr,2.5e6_kr,&
     3.e6_kr,4.e6_kr,5.e6_kr,6.e6_kr,7.e6_kr,8.e6_kr,9.e6_kr,1.e7_kr,&
     1.2e7_kr,1.7e7_kr,3.e7_kr/)
   real(kr),dimension(39),parameter::eg8=(/&
     .01e0_kr,.02e0_kr,.03e0_kr,.045e0_kr,.06e0_kr,.07e0_kr,&
     .075e0_kr,.10e0_kr,.15e0_kr,.20e0_kr,.30e0_kr,.40e0_kr,.45e0_kr,&
     .510e0_kr,.512e0_kr,.60e0_kr,.70e0_kr,.80e0_kr,1.0e0_kr,&
     1.33e0_kr,1.50e0_kr,1.66e0_kr,2.0e0_kr,2.5e0_kr,3.0e0_kr,&
     3.5e0_kr,4.0e0_kr,4.5e0_kr,5.0e0_kr,5.5e0_kr,6.0e0_kr,6.5e0_kr,&
     7.0e0_kr,7.5e0_kr,8.0e0_kr,10.e0_kr,12.e0_kr,14.e0_kr,20.e0_kr/)
   real(kr),dimension(43),parameter::eg10=(/&
     1.0e3_kr,1.0e4_kr,2.0e4_kr,3.0e4_kr,4.5e4_kr,6.0e4_kr,7.0e4_kr,&
     7.5e4_kr,1.0e5_kr,1.50e5_kr,2.00e5_kr,3.00e5_kr,4.00e5_kr,&
     4.50e5_kr,5.10e5_kr,5.12e5_kr,6.00e5_kr,7.00e5_kr,8.00e5_kr,&
     1.00e6_kr,1.33e6_kr,1.34e6_kr,1.50e6_kr,1.66e6_kr,2.00e6_kr,&
     2.50e6_kr,3.00e6_kr,3.50e6_kr,4.00e6_kr,4.50e6_kr,5.00e6_kr,&
     5.50e6_kr,6.00e6_kr,6.50e6_kr,7.00e6_kr,7.50e6_kr,8.00e6_kr,&
     1.00e7_kr,1.20e7_kr,1.40e7_kr,2.00e7_kr,3.00e7_kr,5.00e7_kr/)
   real(kr),parameter::mev=1.e6_kr

   !--select structure
   ngg=0
   egg(1)=0
   if (igg.eq.0) return

   !--group structure is read in.
   if (igg.eq.1) then
      read(nsysi,*) ngg
      ngp=ngg+1
      if (ngp.gt.ngmax) call error('genggp','too many groups.',' ')
      read(nsysi,*) (egg(i),i=1,ngp)

   !--csewg 94 group structure
   else if (igg.eq.2) then
      ngg=94
      do ig=1,95
         egg(ig)=eg2(ig)*mev
      enddo

   !--lanl 12 group structure
   else if (igg.eq.3) then
      ngg=12
      do ig=1,13
         egg(ig)=eg3(ig)*mev
      enddo

   !--steiner 21-group gamma structure (ornl-tm-2564)
   else if (igg.eq.4) then
      ngg=21
      do ig=1,22
         egg(ig)=eg4(ig)*mev
      enddo

   !--straker 22 group structure
   else if (igg.eq.5) then
      ngg=22
      do ig=1,23
         egg(ig)=eg5(ig)*mev
      enddo

   !--lanl 48-group structure
   else if (igg.eq.6) then
      ngg=48
      do ig=1,49
         egg(ig)=eg6(ig)*mev
      enddo

   !--lanl 24-group structure
   else if (igg.eq.7) then
      ngg=24
      do ig=1,25
         egg(ig)=eg7(ig)
      enddo

   !--vitamin-series 36- and 38-group structures
   else if (igg.eq.8.or.igg.eq.9) then
      ngg=38
      if (igg.eq.8) ngg=36
      ngm=ngg+1
      do ig=1,ngm
         egg(ig)=eg8(ig)*mev
      enddo
      if (igg.ne.9) then
         ! remove group bounds eg8(7) and eg8(39) if igg=8
         do ig=7,ngm
            egg(ig)=eg8(ig+1)*mev
         enddo
      endif

   !--vitamin-j 42-group structure
   else if (igg.eq.10) then
      ngg=42
      do ig=1,43
         egg(ig)=eg10(ig)
      enddo

   !--illegal value for igg
   else
      call error('genggp','illegal group structure.',' ')
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
   end subroutine genggp

   subroutine gnwtf
   !------------------------------------------------------------------
   ! Set up calculation of weight functions or read in arbitary
   ! function in the form of an ENDF TAB1 record or
   ! read in parameters for an analytic weight function.
   !
   !    iwt     meaning
   !    ---     -------
   !     1      read in
   !     2      constant
   !     3      1/e + rolloffs
   !
   !------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use util   ! provides error
   ! internals
   integer::i
   real(kr),dimension(16),parameter::wt1=(/&
     0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,1.e0_kr,4.e0_kr,4.e0_kr,&
     5.e0_kr,1.e3_kr,1.e-4_kr,1.e5_kr,1.e0_kr,1.e7_kr,1.e-2_kr,&
     3.e7_kr,1.e-4_kr/)

   !--branch on weight function option

   !--arbitary
   if (iwt.eq.1) then
      write(nsyso,'(/'' weight function......read in'')')
      read(nsysi,*) (wght(i),i=1,iwmax)

   !--constant
   else if (iwt.eq.2) then
      write(nsyso,'(/'' weight function......constant for all l'')')

   !--1/e with high and low energy rolloffs
   else if (iwt.eq.3) then
      write(nsyso,'(/'' weight function......1/e with rolloffs'')')
      do i=1,16
         wght(i)=wt1(i)
      enddo

   !--illegal value for iwt
   else
      call error('gnwtf','illegal iwt',' ')
   endif
   return
   end subroutine gnwtf

   subroutine gtflx(e,enext,idis,flux,nl,nz)
   !------------------------------------------------------------------
   ! Retrieve or compute required legendre component of the
   ! weight function constructed or read in by genwtf.
   !------------------------------------------------------------------
   use endf ! provides terpa
   ! externals
   integer::idis,nl,nz
   real(kr)::e,enext,flux(10,10)
   ! internals
   integer::ip,ir,il
   real(kr)::wtf,enxt
   real(kr),parameter::emax=1.e12_kr
   real(kr),parameter::step=1.05e0_kr
   real(kr),parameter::zero=0
   save ip,ir

   !--initialize
   idis=0
   if (e.eq.zero) then
      ip=2
      ir=1
      enext=emax

   !--branch to desired method
   else

      !--tabulated
      if (iwt.ne.2) then
         call terpa(wtf,e,enext,idis,wght,ip,ir)
         do il=1,nl
            flux(1,il)=wtf
         enddo
         enxt=step*e
         if (enext.gt.enxt) idis=0
         if (enext.gt.enxt) enext=enxt

     !--constant for all orders
      else
         do il=1,nl
            flux(1,il)=1
         enddo
         enext=emax
         idis=0
      endif
   endif
   return
   end subroutine gtflx

   subroutine gpanel(elo,ehi,ans,ff,nl,nz,ng,iglo,sdat,pff,nwpff)
   !------------------------------------------------------------------
   ! Perform generalized group constant integrals for one panel.
   ! The upper boundry of the panel is chosen to be the smallest
   ! of ehi, the next cross section point, the next flux point,
   ! and the next feed function point.  Use Lobatto quadrature
   ! with order two larger than that used for the feed function.
   !------------------------------------------------------------------
   use util ! provides error
   ! externals
   integer::nl,nz,ng,iglo,nwpff
   real(kr)::elo,ehi,ans(nl,nz,*),ff(nl,*),sdat(*),pff(nwpff)
   ! internals
   integer::idiscf,ng1,ig1,nq,iz,il,iq,nqp,ig,igt
   real(kr)::en,ehigh,aq,bq,eq,wq,t1,a,b,rr,enext
   real(kr)::sig(10,10),slst(10,10),flux(10,10),flst(10,10)
   real(kr),dimension(2),parameter::qp2=(/-1.e0_kr,1.e0_kr/)
   real(kr),dimension(2),parameter::qw2=(/1.e0_kr,1.e0_kr/)
   real(kr),dimension(6),parameter::qp6=(/&
     -1.e0_kr,-.76505532e0_kr,-.28523152e0_kr,.28523152e0_kr,&
     .76505532e0_kr,1.e0_kr/)
   real(kr),dimension(6),parameter::qw6=(/&
     .06666667e0_kr,.37847496e0_kr,.55485838e0_kr,.55485838e0_kr,&
     .37847496e0_kr,.06666667e0_kr/)
   real(kr),dimension(10),parameter::qp10=(/&
     -1.e0_kr,-.9195339082e0_kr,-.7387738651e0_kr,-.4779249498e0_kr,&
     -.1652789577e0_kr,.1652789577e0_kr,.4779249498e0_kr,&
     .7387738651e0_kr,.9195339082e0_kr,1.e0_kr/)
   real(kr),dimension(10),parameter::qw10=(/&
     .0222222222e0_kr,.1333059908e0_kr,.2248893420e0_kr,&
     .2920426836e0_kr,.3275397612e0_kr,.3275397612e0_kr,&
     .2920426836e0_kr,.2248893420e0_kr,.1333059908e0_kr,&
     .0222222222e0_kr/)
   real(kr),parameter::rndoff=1.000002e0_kr
   real(kr),parameter::delta=0.999995e0_kr
   integer::idisc=0
   real(kr)::elast=0
   save nq,enext,elast,slst,flst,idisc,ng1,ig1

   !--retrieve factors in integrands at lower boundry.
   if (elo.gt.ehi) call error('gpanel','elo gt ehi.',' ')
   if (elo.ne.elast) then
      if (elo*rndoff.lt.ehi) elo=elo*rndoff
      elast=elo
      call gtsig(elo,enext,idisc,slst,sdat)
      call gtflx(elo,en,idiscf,flst,nl,nz)
      if (en.eq.enext.and.idiscf.gt.idisc) idisc=idiscf
      if (en.lt.enext) idisc=idiscf
      if (en.lt.enext) enext=en
      call gtff(elo,en,idiscf,ff,nl,ng1,ig1,nq,pff,nwpff)
      if (en.eq.enext.and.idiscf.gt.idisc) idisc=idiscf
      if (en.lt.enext) idisc=idiscf
      if (en.lt.enext) enext=en
      nq=nq+2
      if (nq.gt.10) nq=10
   endif

   !--retrieve cross section and flux at upper boundry.
   if (enext.lt.delta*ehi) then
      ehi=enext
      ehigh=ehi
      if (idisc.gt.0.and.ehi*delta.gt.elo) ehigh=ehi*delta
   else
      ehigh=delta*ehi
   endif
   call gtsig(ehigh,enext,idisc,sig,sdat)
   call gtflx(ehigh,en,idiscf,flux,nl,nz)
   if (en.eq.enext.and.idiscf.gt.idisc) idisc=idiscf
   if (en.lt.enext) idisc=idiscf
   if (en.lt.enext) enext=en

   !--compute group fluxes assuming flux is linear over panel.
   aq=(ehigh+elo)/2
   bq=(ehigh-elo)/2
   do iz=1,nz
      do il=1,nl
         ans(il,iz,1)=ans(il,iz,1)+(flux(iz,il)+flst(iz,il))*bq
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
      endif
      t1=(eq-elo)/(ehi-elo)

      !--retrieve the feed function at the quadrature point.
      !--first point was last point of previous panel.
      if (iq.gt.1) then
         if (eq.gt.ehigh) eq=ehigh
         call gtff(eq,en,idiscf,ff,nl,ng1,ig1,nqp,pff,nwpff)
      endif

      !--accumulate the ng*nl*nz integrals simultaneously.
      !--assuming that the reaction rate is linear across the panel.
      if (iglo.eq.0) iglo=ig1
      do iz=1,nz
         do il=1,nl
            a=sig(iz,1)*flux(iz,il)
            b=slst(iz,1)*flst(iz,il)
            rr=b+(a-b)*t1
            rr=rr*wq
            do ig=1,ng1
               igt=ig1+ig-iglo+1
               if (igt.gt.1) then
                  ans(il,iz,igt)=ans(il,iz,igt)+rr*ff(il,ig)
                  if (igt.gt.ng) ng=igt
               endif
            enddo
         enddo
      enddo
   enddo

   !--save cross section and flux for next panel.
   !--determine next point and quadrature order from last ff.
   elast=ehigh
   do iz=1,nz
      slst(iz,1)=sig(iz,1)
      do il=1,nl
         flst(iz,il)=flux(iz,il)
      enddo
   enddo
   if (en.eq.enext.and.idiscf.gt.idisc) idisc=idiscf
   if (en.lt.enext) idisc=idiscf
   if (en.lt.enext) enext=en
   if (enext.le.ehi) enext=rndoff*ehi
   nq=nqp+2
   if (nq.gt.10) nq=10
   return
   end subroutine gpanel

   subroutine dspla(ig,ans,nl,nz,ng2,ig2lo,igzero)
   !------------------------------------------------------------------
   ! Display generalized group constants generated by panel.
   !------------------------------------------------------------------
   use mainio ! provides nsyso
   use util ! provides a10
   ! externals
   integer::ig,nl,nz,ng2,ig2lo,igzero
   real(kr)::ans(nl,nz,*)
   ! internals
   integer::i,ig2,igt,max,il
   real(kr)::result(20)
   character(10)::field(20)
   real(kr),parameter::small=1.e-9_kr
   real(kr),parameter::zero=0

   !--write appropriate heading.
   igzero=0
   if (ig.lt.0) then
      if (iprint.eq.1) then
         write(nsyso,'('' '')')
         if (mfd.ne.26) then
            if (mtd.eq.501.or.mtd.eq.502.or.mtd.eq.504)&
              write(nsyso,'(&
              &'' gamma  sigma  ''/'' group  (barns)'')')
            if (mtd.eq.516) write(nsyso,'(&
              &'' gamma  sigma  ''/'' group  (barns)'')')
            if (mtd.eq.602) write(nsyso,'(&
              &'' gamma  sigma      heating''/&
              &'' group  (barns)    (ev-barns)'')')
            if (mtd.eq.621) write(nsyso,'(&
              &'' gamma  heating''/&
              &'' group  (ev-barns)'')')
            if (mtd.eq.522) write(nsyso,'(&
              &'' gamma  sigma      heating''/&
              &'' group  (barns)    (ev-barns)'')')
            if (mtd.eq.525) write(nsyso,'(&
              &'' gamma  heating''/&
              &'' group  (ev-barns)'')')
         else
            if (mtd.eq.502.or.mtd.eq.504) write(nsyso,'(&
              &'' initl  final  cross sections vs''/&
              &'' group  group  legendre order   '')')
            if (mtd.eq.516) write(nsyso,'(&
              &'' initl  final  cross sections vs''/&
              &'' group  group  legendre order   '')')
         endif
      endif
      return
   endif

   !--write out results in appropriate format.
   !--cross sections and heating values.
   if (mfd.eq.23) then
      do i=2,ng2
         result(i-1)=0
         if (ans(1,1,1).ne.zero.and.abs(ans(1,1,i)).ge.small) then
            result(i-1)=ans(1,1,i)/ans(1,1,1)
            ans(1,1,i)=result(i-1)
            if (ans(1,1,i).ge.small) igzero=1
         endif
      enddo
      if (result(1).ne.zero.and.iprint.eq.1) then
         do i=2,ng2
            call a10(result(i-1),field(i-1))
         enddo
         write(nsyso,'(1x,i3,2x,1p,2a11)') ig,(field(i-1),i=2,ng2)
      endif

   !--scattering matrices.
   else if (mfd.eq.26) then
      do ig2=2,ng2
         igt=ig2lo+ig2-2
         if (mtd.eq.516.and.ig2.gt.2) then
            result(1)=ans(1,1,ig2)/ans(1,1,1)
            ans(1,1,ig2)=result(1)
            if (ans( 1,1,ig2).ge.small) igzero=1
            if (iprint.eq.1) then
               write(nsyso,'(1x,i3,3x,''heat'',2x,1p,e11.3)')&
                 ig,result(1)
            endif
         else if (igt.le.ig) then
            max=nl
            if (max.gt.6) max=6
            do il=1,max
               result(il)=ans(il,1,ig2)/ans(il,1,1)
               ans(il,1,ig2)=result(il)
               if (ans(il,1,ig2).ge.small) igzero=1
               call a10(result(il),field(il))
            enddo
            if (iprint.eq.1)&
              write(nsyso,'(1x,i3,4x,i3,2x,1p,6a11)')&
              ig,igt,(field(i),i=1,max)
            if (nl.gt.6) then
               do il=7,nl
                  result(il)=ans(il,1,ig2)/ans(il,1,1)
                  ans(il,1,ig2)=result(il)
                  if (ans(il,1,ig2).ge.small) igzero=1
                  call a10(result(il),field(il))
               enddo
               if (iprint.eq.1)&
                 write(nsyso,'(13x,1p,6a11)') (field(i),i=7,nl)
            endif
         else
            result(1)=ans(1,1,ig2)/ans(1,1,1)
            ans(1,1,ig2)=result(1)
            if (ans( 1,1,ig2).ge.small) igzero=1
            call a10(result(1),field(1))
            if (iprint.eq.1) then
               if (igt.eq.ig+1) write(nsyso,'(&
                 &1x,i3,3x,''xsec'',2x,1p,a11)') ig,field(1)
               if (igt.eq.ig+2) write(nsyso,'(&
                 &1x,i3,3x,''heat'',2x,1p,a11)') ig,field(1)
            endif
         endif
      enddo
   endif
   return
   end subroutine dspla

   subroutine gtsig(e,enext,idis,sig,sdat)
   !------------------------------------------------------------------
   ! Retrieve the reaction cross-section defined by mfd and mtd.
   !------------------------------------------------------------------
   use endf ! provides endf routines and variables
   ! externals
   integer::idis
   real(kr)::e,enext,sig(10,10),sdat(*)
   ! internals
   integer::nb,nw,mf,mt
   real(kr)::s
   real(kr),parameter::zero=0

   !--initialize
   if (e.eq.zero) then
      mf=23
      mt=mtd
      call findf(matd,mf,mt,npend)
      call contio(npend,0,0,sdat,nb,nw)
      call gety1(e,enext,idis,s,npend,sdat)

   !--normal entry
   else
      call gety1(e,enext,idis,s,npend,sdat)
      sig(1,1)=s
   endif
   return
   end subroutine gtsig

   subroutine gtff(e,enext,idisc,ff,nl,ng,iglo,nq,pff,nwpff)
   !------------------------------------------------------------------
   ! Compute feed function or yield for desired reaction type.
   !------------------------------------------------------------------
   use util ! provides error
   use endf ! provides endf routines and variables
   use mathm ! provides legndr
   use physics ! provides electron mass
   ! externals
   integer::idisc,nl,ng,iglo,nq,nwpff
   real(kr)::e,enext,ff(nl,*),pff(nwpff)
   ! internals
   integer::igp,il,ip,idis,idone,ig,ifini,i,ir,iq,l
   integer::lim,j,nb,nw
   real(kr)::c1e,aq,bq,fact,enow,enowi,enow2,ebar,dk,zz
   real(kr)::pnow,pnext,pnowi,xnow,unow,snow,xlim,stest,wq,test
   real(kr)::xnext,xq,xzz
   real(kr)::sigcoh,q2m,siginc,q2,unext,px,uq,rm2,rm,rt
   real(kr)::arg(12),pl(12)
   save zz

   real(kr),dimension(6),parameter::qp6=(/&
      -1.e0_kr,-.76505532e0_kr,-.28523152e0_kr,.28523152e0_kr,&
     .76505532e0_kr,1.e0_kr/)
   real(kr),dimension(6),parameter::qw6=(/&
     .06666667e0_kr,.37847496e0_kr,.55485838e0_kr,.55485838e0_kr,&
     .37847496e0_kr,.06666667e0_kr/)
   real(kr),dimension(10),parameter::qp10=(/&
     -1.e0_kr,-.9195339082e0_kr,-.7387738651e0_kr,-.4779249498e0_kr,&
     -.1652789577e0_kr,.1652789577e0_kr,.4779249498e0_kr,&
     .7387738651e0_kr,.9195339082e0_kr,1.e0_kr/)
   real(kr),dimension(10),parameter::qw10=(/&
     .0222222222e0_kr,.1333059908e0_kr,.2248893420e0_kr,&
     .2920426836e0_kr,.3275397612e0_kr,.3275397612e0_kr,&
     .2920426836e0_kr,.2248893420e0_kr,.1333059908e0_kr,&
     .0222222222e0_kr/)
   real(kr),parameter::c1=57.03156e-6_kr
   real(kr),parameter::c2=0.249467e0_kr
   real(kr),parameter::c3=1.95693e-6_kr
   real(kr),parameter::c4=0.0485262e0_kr
   real(kr),parameter::c5=20.60744e0_kr
   real(kr),parameter::toler=1.e-6_kr
   real(kr),parameter::rndoff=1.0000001e0_kr
   real(kr),parameter::emax=1.e12_kr
   real(kr),parameter::two=2.e0_kr
   real(kr),parameter::half=0.5e0_kr
   real(kr),parameter::close=.99999e0_kr
   real(kr),parameter::zero=0

   !--photon interaction cross sections and heating.
   if (mfd.eq.23) then
      nq=0
      iglo=1
      ng=0
      if (mtd.eq.501.or.mtd.eq.502.or.mtd.eq.504) ng=1
      if (mtd.eq.516) ng=1
      if (mtd.eq.602) ng=2
      if (mtd.eq.522) ng=2
      if (ng.eq.0)&
        call error('gtff','illegal reaction for cross section.',' ')
      enext=emax
      idisc=0
      if (e.eq.zero) return
      ff(1,1)=1
      if (mtd.eq.602) ff(1,2)=e
      if (mtd.eq.522) ff(1,2)=e
      return

   !--photon scattering cross sections.
   !--initialize.
   else if (mfd.eq.26) then
      if (e.eq.zero.or.mtd.ne.502) then
         if (e.eq.zero.or.mtd.ne.504) then
            if (mtd.ne.516) then
               call findf(matd,27,mtd,nendf)
               l=1
               call contio(nendf,0,0,pff(l),nb,nw)
               call tab1io(nendf,0,0,pff(l),nb,nw)
               zz=pff(l+1)
               l=l+nw
               do while (nb.ne.0)
                  if (l.gt.nwpff) call error('gtff',&
                    'insufficient storage for form factor.',' ')
                  call moreio(nendf,0,0,pff(l),nb,nw)
                  l=l+nw
               enddo
               ng=1
               if (mtd.eq.504) ng=ngg+2
               enext=emax
               idisc=0
               return
            endif
         endif
      endif

      !--photon coherent scattering.
      !--compute all legendre components at this energy.
      if (mtd.eq.502) then
         igp=1
         do while (e.ge.egg(igp+1))
            igp=igp+1
         enddo
         do il=1,nl
            ff(il,1)=0
         enddo
         ng=1
         iglo=igp
         xnow=0
         unow=1
         ip=2
         ir=1
         call terpa(snow,xnow,xnext,idis,pff,ip,ir)
         fact=2*c2/(c1*c1*e*e)
         xlim=sqrt(two)*c1*e
         c1e=1/(c1*e)**2
         do il=1,nl
            arg(il)=2*snow**2
         enddo
         stest=toler*snow

         !--loop over panels of coherent form factor.
         idone=0
         do while (idone.eq.0)
            aq=(xnext+xnow)/2
            bq=(xnext-xnow)/2
            nq=6
            if (nl.gt.4) nq=10
            do iq=1,nq
               if (nq.eq.6) then
                  xq=aq+bq*qp6(iq)
                  wq=bq*qw6(iq)
               else
                  xq=aq+bq*qp10(iq)
                  wq=bq*qw10(iq)
               endif
               xnow=xq*rndoff
               if (iq.gt.1) then
                  unow=1-c1e*xnow*xnow
                  test=-1
                  if (unow.lt.test) unow=-1
                  call terpa(snow,xnow,xnext,idis,pff,ip,ir)
                  call legndr(unow,pl,nl)
                  do il=1,nl
                     arg(il)=(1+unow*unow)*snow*snow*pl(il)
                  enddo
               endif
               do il=1,nl
                  ff(il,1)=ff(il,1)+wq*fact*xnow*arg(il)
               enddo
            enddo
            test=-1
            if (unow.le.test) then
               idone=1
            else
               if (xnow.ge.xlim) then
                  idone=1
               else
                  if (snow.lt.stest) then
                     idone=1
                  else
                     if (xnext.gt.xlim) xnext=xlim
                  endif
               endif
            endif
         enddo
         sigcoh=ff(1,1)
         do il=1,nl
            ff(il,1)=ff(il,1)/sigcoh
         enddo
         enext=egg(igp+1)
         idisc=1
         nq=0
         return

      !--photon incoherent scattering and heating.
      !--compute for all legendre orders and all sink groups.
      else if (mtd.eq.504) then
         igp=1
         do while (e.ge.egg(igp+1))
            igp=igp+1
         enddo
         enow=c3*e
         enowi=1/enow
         enow2=enow*enow
         pnow=enow
         xnow=0
         xzz=c5*sqrt(enow/500)
         q2m=(2*enow*(1+enow)/(1+2*enow))**2
         if (xzz.le.zz) then
            ip=2
            ir=1
            call terpa(snow,xnow,xnext,idis,pff,ip,ir)
         else
            snow=zz
            xnext=emax
         endif
         do il=1,nl
            arg(il)=0
            lim=igp+2
            do ig=1,lim
               ff(il,ig)=0
            enddo
         enddo
         ig=igp
         siginc=0
         ebar=0
         if (e.eq.egg(igp)) igp=igp-1

         !--loop over panels defined by group and form factor breaks.
         ifini=0
         do while (ifini.eq.0)
            idone=0
            do while (idone.eq.0)
               q2=(c4*xnext)**2
               unext=-1
               if (q2.gt.q2m) then
                  idone=1
               else
                  unext=1-((1-q2*enowi)-sqrt(1+q2))/(q2-enow2-2*enow)
                  test=-1
                  if (unext.lt.test) unext=-1
                  if (unext.lt.close) then
                     idone=1
                  else
                     xnow=xnext*rndoff
                     call terpa(snow,xnow,xnext,idis,pff,ip,ir)
                  endif
               endif
            enddo
            pnext=enow/(1+2*enow)
            if (igp.gt.0) pnext=c3*egg(igp)
            px=enow/(1+enow*(1-unext))
            if (px.gt.pnext) pnext=px
            if (pnext.gt.pnow/rndoff) pnext=pnow/rndoff
            aq=(pnext+pnow)/2
            bq=(pnext-pnow)/2
            nq=6
            if (nl.gt.6) nq=10
            do iq=1,nq
               if (nq.eq.6) then
                  uq=aq+bq*qp6(iq)
                  wq=-c2*bq*qw6(iq)
               else
                  uq=aq+bq*qp10(iq)
                  wq=-c2*bq*qw10(iq)
               endif
               pnow=uq
               if (iq.ne.1) then
                  pnowi=1/pnow
                  unow=1+enowi-pnowi
                  test=1
                  if (unow.gt.test) unow=1
                  if (xzz.le.zz) then
                     rm2=(1-unow)/2
                     rm=sqrt(rm2)
                     rt=1+2*enow*rm2
                     xnow=c5*2*enow*rm*sqrt(rt+enow2*rm2)/rt
                     xnow=xnow*rndoff
                     call terpa(snow,xnow,xnext,idis,pff,ip,ir)
                  endif
                  call legndr(unow,pl,nl)
                  dk=unow-1
                  fact=snow*(enow*pnowi+pnow*enowi+dk*(2+dk))/enow2
                  do il=1,nl
                     arg(il)=fact*pl(il)
                  enddo
               endif
               if (igp.ne.0) then
                  do il=1,nl
                     ff(il,igp)=ff(il,igp)+wq*arg(il)
                  enddo
               endif
               siginc=siginc+wq*arg(1)
               ebar=ebar+wq*arg(1)*pnow/c3
            enddo
            if (unow.lt.-close) then
               ifini=1
            else
               if (igp.gt.0) then
                  if (pnext.le.c3*egg(igp)) igp=igp-1
               endif
            endif
         enddo

         !--strip sink groups with zero cross section.
         if (igp.le.0) igp=1
         iglo=igp
         ng=ig-igp+1
         do i=1,ng
            j=iglo+i-1
            do il=1,nl
               ff(il,i)=ff(il,j)
            enddo
         enddo
         ebar=ebar/siginc
         ff(1,ng+1)=siginc
         ff(1,ng+2)=(e-ebar)*siginc
         ng=ng+2
         do il=1,nl
            do i=1,ng
               ff(il,i)=ff(il,i)/siginc
            enddo
         enddo

         !--determine next critical point for incoherent scattering.
         if (c3*e.ge.half) then
            enext=emax
         else
            idone=0
            i=0
            do while (i.lt.ngg.and.idone.eq.0)
               i=i+1
               if (c3*egg(i).ge.test) then
                  idone=1
               else
                  enext=egg(i)/(1-2*c3*egg(i))
                  if (enext.gt.e) idone=2
               endif
            enddo
            if (idone.lt.2) enext=emax
         endif
         if (egg(ig+1).lt.enext) enext=egg(ig+1)
         idisc=1
         nq=nq+2
         return

      !--pair production matrix
      else if (mtd.eq.516) then
         enext=emax
         idisc=0
         if (e.gt.zero) then
            igp=1
            do while (e.ge.egg(igp+1))
               igp=igp+1
            enddo
            do il=1,nl
               ff(il,1)=0
            enddo
            ng=2
            iglo=ig2pp
            nq=0
            ff(1,1)=2
            ff(1,2)=e-2.0*epair
         endif
         return
      endif

   !--bad file type
   else
      call error('gtff','illegal file type.',' ')
   endif
   return
   end subroutine gtff

end module gaminm


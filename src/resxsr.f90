module resxsm
   ! provides subroutine resxsr for NJOY2016
   use locale
   implicit none
   private
   public resxsr

contains

   subroutine resxsr
   !-------------------------------------------------------------------
   !
   !  Construct a RESXS resonance cross section
   !  file from NJOY pendf cross sections.
   !
   !  User input:
   !
   !   card 1   units
   !      nout     output unit
   !
   !   card 2   control
   !      nmat     number of materials
   !      maxt     max. number of temperatures
   !      nholl    number of lines of descriptive comments
   !      efirst   lower energy limit (ev)
   !      elast    upper energy limit
   !      eps      thinning tolerance
   !
   !   card 3   user id
   !      huse     hollerith user identification (up to 16 chars)
   !      ivers    file version number
   !
   !   card 4   descriptive data (repeat nholl times)
   !      holl     line of hollerith data (72 chars max)
   !
   !   card 5   material specifications (repeat nmat times)
   !      hmat     hollerith name for material (up to 8 chars)
   !      mat      endf mat number for material
   !      unit     njoy unit number for pendf data
   !
   !  The RESXS format specification follows:
   !
!***********************************************************************
!               proposed 09/24/90                                      -
!                                                                      -
!f           resxs                                                     -
!e           resonance cross section file                              -
!                                                                      -
!n                       this file contains pointwise cross            -
!n                       sections for the epithermal resonance         -
!n                       range to be used for hyper-fine flux          -
!n                       calculations.  elastic, fission, and          -
!n                       capture cross sections are given vs           -
!n                       temperature.  linear interpolation is         -
!n                       assumed.                                      -
!                                                                      -
!n           formats given are for file exchange only                  -
!                                                                      -
!***********************************************************************
!
!
!-----------------------------------------------------------------------
!s          file structure                                             -
!s                                                                     -
!s              record type                       present if           -
!s              ==============================    ===============      -
!s              file identification                 always             -
!s              file control                        always             -
!s              set hollerith identification        always             -
!s              file data                           always             -
!s                                                                     -
!s   *************(repeat for all materials)                           -
!s   *          material control                    always             -
!s   *                                                                 -
!s   * ***********(repeat for all cross section blocks)                -
!s   * *        cross section block                 always             -
!s   * ***********                                                     -
!s   *************                                                     -
!                                                                      -
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!r           file identification                                       -
!                                                                      -
!l    hname,(huse(i),i=1,2),ivers                                      -
!                                                                      -
!w    1+3*mult                                                         -
!                                                                      -
!b    format(4h ov ,a8,1h*,2a8,1h*,i6)                                 -
!                                                                      -
!d    hname         hollerith file name  - resxs -  (a8)               -
!d    huse          hollerith user identifiation    (a8)               -
!d    ivers         file version number                                -
!d    mult          double precision parameter                         -
!d                       1- a8 word is single word                     -
!d                       2- a8 word is double precision word           -
!                                                                      -
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!r           file control                                              -
!                                                                      -
!l    efirst,elast,nholl,nmat,nblok
!                                                                      -
!w    5                                                                -
!                                                                      -
!b    format(4h 1d ,2i6)                                               -
!                                                                      -
!d    efirst      lowest energy on file (ev)                           -
!d    elast       highest enery on file (ev)                           -
!d    nholl       number of words in set hollerith                     -
!d                    identification record                            -
!d    nmat        number of materials on file                          -
!d    nblok       energy blocking factor                               -
!                                                                      -
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!r           set hollerith identification                              -
!                                                                      -
!l    (hsetid(i),i=1,nholl)                                            -
!                                                                      -
!w    nholl*mult                                                       -
!                                                                      -
!b    format(4h 2d ,8a8/(9a8))                                         -
!                                                                      -
!d    hsetid      hollerith identification of set (a8)                 -
!d                 (to be edited out 72 characters per line)           -
!                                                                      -
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!r          file data                                                  -
!                                                                      -
!l    (hmatn(i),i=1,nmat),(ntemp(i),i=1,nmat),(locm(i),i=1,nmat)       -
!                                                                      -
!w    (mult+2)*nmat                                                    -
!                                                                      -
!b    format(4h 3d ,8a8/(9a8))      hmatn                              -
!b    format(12i6)                  ntemp,locm                         -
!                                                                      -
!d    hmatn(i)    hollerith identification for material i              -
!d    ntemp(i)    number of temperatures for material i                -
!d    locm(i)     location of material i                               -
!                                                                      -
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!r          material control                                           -
!                                                                      -
!l    hmat,amass,(temp(i),i=1,ntemp),nreac,nener                       -
!                                                                      -
!w    mult+3+ntemp                                                     -
!                                                                      -
!b    format(4h 6d ,a8,1h*,1p1e12.5/(6e12.5))     hmat,temp            -
!b    format(2i6)                                 nener,blok           -
!                                                                      -
!d    hmat        hollerith material identifier                        -
!d    amass       atomic weight ratio                                  -
!d    temp        temperature values for this material                 -
!d    nreac       number of reactions for this material                -
!d                  (3 for fissionable, 2 for nonfissionable)          -
!d    nener       number of energies for this material                 -
!                                                                      -
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!r          cross section block                                        -
!                                                                      -
!l    (xsb(i),i=1,imax)                                                -
!                                                                      -
!c    imax=3*ntemp*(number of energies in the block)                   -
!                                                                      -
!w    imax                                                             -
!                                                                      -
!b    format(4h 8d ,1p5e12.5/(6e12.5))                                 -
!                                                                      -
!d    xsb(i)      data for a block of nblok or fewer point energy      -
!d                values.  the data values given for each energy       -
!d                are nelas, nfis, and ng at temp(1), followed by      -
!d                nelas, nfis, and ng at temp(2), and so on.           -
!                                                                      -
!-----------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use endf ! provides endf routines and variables
   use util ! provides openz,loada,finda
   ! internals
   integer::iscr,iscrh,nscr,inew,iold,irec
   integer::nout,nmat,maxt,nholl,ivers,i,j,k,im,matd
   integer::nin,nb,nw,ix,itemp,jx,idis,isave,ie,je
   integer::ne,n,lim,nwds,nn,jrec,krec,nwr,nreac,nener
   real(kr)::e,efirst,elast,eps,awr,enext,s,ei,test
   real(kr)::temp(10),sigs(31),sigi(31),sigl(31)
   character(8)::hmat(20)
   integer::mat(20),ninm(20),ntemp(20),locm(20)
   real(kr)::awrm(20)
   real(kr)::huse(2)
   character(8)::holl(9,10)
   real(kr)::stack(50,31)
   character(12)::uid
   character(72)::text
   real(kr)::edat(5000)
   real(k4)::a(10000)
   integer(k4)::ia(10000)
   real(k8)::ha(5000)
   equivalence (a(1),ia(1),ha(1))
   real(kr),dimension(:),allocatable::bufo,bufn
   character(6)::hfile='resxsr'
   character(6)::hblank='      '
   integer::mult=2
   integer::nbuf=5000
   integer::nx=31
   integer::nblok=5000

   !--assign storage pointers
   iscr=1
   iscrh=1
   allocate(bufo(nbuf))
   allocate(bufn(nbuf))

   !--assign i/o units
   nscr=10
   inew=11
   iold=12
   call openz(-nscr,0)
   call openz(-inew,0)
   call openz(-iold,0)
   irec=0

   !--read user input
   read(nsysi,*) nout
   read(nsysi,*) nmat,maxt,nholl,efirst,elast,eps
   read(nsysi,*) uid,ivers
   read(uid,'(2a6)') huse(1),huse(2)
   do i=1,nholl
      text=' '
      read(nsysi,*) text
      read(text,'(9a8)') (holl(j,i),j=1,9)
   enddo
   do i=1,nmat
      read(nsysi,*) hmat(i),mat(i),ninm(i)
   enddo

   !--loop over desired materials
   do 400 im=1,nmat
   matd=mat(im)
   locm(im)=irec
   nin=ninm(im)
   call openz(nin,0)
   call tpidio(nin,0,0,edat,nb,nw)

   !--initialize loada/finda file
   do ix=1,nx
      sigs(ix)=0
   enddo
   sigs(1)=efirst
   call loada(1,sigs,nx,inew,bufn,nbuf)
   sigs(1)=elast
   call loada(-2,sigs,nx,inew,bufn,nbuf)

   !--loop over temperatures
   itemp=0
   jx=0
  110 continue
   call contio(nin,0,0,edat,nb,nw)
   if (math.ne.matd) go to 250
   if (itemp.eq.maxt) go to 250
   awr=c2h
   awrm(im)=awr
   call contio(nin,0,0,edat,nb,nw)
   if (n1h.ne.0) then
      iverf=4
   else if (n2h.eq.0) then
      iverf=5
   else
      iverf=6
   endif
   call skiprz(nin,-1)
   if (iverf.eq.4) nx=n2h
   if (iverf.ge.5) call contio(nin,0,0,edat,nb,nw)
   if (iverf.ge.6) call contio(nin,0,0,edat,nb,nw)
   call hdatio(nin,0,0,edat,nb,nw)
   itemp=itemp+1
   ntemp(im)=itemp
   temp(itemp)=edat(1)
   write(nsyso,'('' found: '',a8,i5,1p,e12.4)')&
     hmat(im),math,temp(itemp)

   !--loop over resonance reactions
   call findf(matd,3,0,nin)
  120 continue
   call contio(nin,0,0,edat,nb,nw)
   if (mfh.eq.0) go to 210
   if (mth.eq.2.or.mth.eq.18.or.mth.eq.102) go to 130
   call tosend(nin,0,0,edat)
   go to 120
  130 continue
   jx=jx+1

   !--loop through the energy grid for this reaction
   !--and add the cross sections to the union set
   e=0
   call gety1(e,enext,idis,s,nin,edat)
   isave=iold
   iold=inew
   inew=isave
   ie=1
   call finda(ie,sigs,nx,iold,bufo,nbuf)
   e=sigs(1)
   call gety1(e,enext,idis,s,nin,edat)
   sigs(1+jx)=s
   je=1
   call loada(je,sigs,nx,inew,bufn,nbuf)
  150 continue
   ie=ie+1
   do ix=1,nx
      sigl(ix)=sigs(ix)
   enddo
   call finda(ie,sigi,nx,iold,bufo,nbuf)
   ei=sigi(1)
  160 continue
   e=ei
   if (e.gt.enext) e=enext
   call gety1(e,enext,idis,s,nin,edat)
   if (jx.eq.1) go to 180
   do ix=2,jx
      call terp1(sigl(1),sigl(ix),sigi(1),sigi(ix),e,sigs(ix),2)
   enddo
  180 continue
   sigs(1)=e
   sigs(1+jx)=s
   je=je+1
   if (e.gt.elast-elast/200000) go to 190
   call loada(je,sigs,nx,inew,bufn,nbuf)
   if (e.eq.ei) go to 150
   go to 160
  190 continue
   call loada(-je,sigs,nx,inew,bufn,nbuf)
   write(nsyso,'(''  mt='',i4,4x,''ne='',i6)') mth,je
   call tosend(nin,0,0,edat)
   go to 120

   !--continue temperature loop for this material
  210 continue
   call tomend(nin,0,0,edat)
   go to 110

   !--thin the cross sections for this material
  250 continue
   ne=je
   isave=inew
   inew=iold
   iold=isave
   ie=0
   je=0
   n=0
  260 continue
   ie=ie+1
   call finda(ie,sigs,nx,iold,bufo,nbuf)
   if (ie.eq.ne) go to 350
   n=n+1
   do j=1,nx
      stack(n,j)=sigs(j)
   enddo
   if (n.lt.3) go to 260
   lim=n-1
   do 290 i=2,lim
   do 280 j=2,nx
   call terp1(stack(1,1),stack(1,j),stack(n,1),stack(n,j),&
     stack(i,1),test,2)
   if (abs(test-stack(i,j)).gt.eps*stack(i,j)) go to 300
  280 continue
  290 continue
   go to 260
  300 continue
   do j=1,nx
      sigs(j)=stack(1,j)
   enddo
   je=je+1
   call loada(je,sigs,nx,inew,bufn,nbuf)
   do i=1,2
      do j=1,nx
         stack(i,j)=stack(lim+i-1,j)
      enddo
   enddo
   n=2
   go to 260
  350 continue
   je=je+1
   call loada(-je,sigs,nx,inew,bufn,nbuf)
   write(nsyso,'('' after thinning, ne='',i6)') je

   !--write the resxs material control record
   n=ntemp(im)
   nwds=mult+n+3
   irec=irec+1
   read(hmat(im),'(a6)') ha(iscrh)
   a(iscr+mult)=awrm(im)
   do i=1,n
      a(iscr+mult+i)=temp(i)
   enddo
   ia(iscr+mult+n+1)=jx/n
   ia(iscr+mult+n+2)=je
   write(nscr) ha(1),(a(mult+i),i=iscr,iscr+n),ia(nwds-1),ia(nwds)

   !--write the thinned cross-section blocks
   nn=1+jx
   nwds=0
   do i=1,je
      call finda(i,sigs,nx,inew,bufn,nbuf)
      if (nwds+nn.gt.nblok) then
         irec=irec+1
         write(nscr)(a(iscr-1+j),j=1,nwds)
         nwds=0
      endif
      do j=1,nn
         a(iscr-1+nwds+j)=sigs(j)
      enddo
      nwds=nwds+nn
   enddo
   if (nwds.ne.0) then
      irec=irec+1
      write(nscr)(a(iscr-1+i),i=1,nwds)
   endif

   !--continue the material loop
  400 continue

   !--write the output in resxs format
   !--start with the file identification record
   call openz(-nout,1)
   nwds=3*mult+1
   read(hfile,'(a6)') ha(iscrh)
   ha(iscrh+1)=huse(1)
   ha(iscrh+2)=huse(2)
   ia(iscr+3*mult)=ivers
   jrec=1
   call repoz(-nout)
   write(nout)(ha(i),i=1,3),ia(nwds)

   !--file control record
   nwds=5
   jrec=2
   a(iscr)=efirst
   a(iscr+1)=elast
   ia(iscr+2)=nholl
   ia(iscr+3)=nmat
   ia(iscr+4)=nblok
   write(nout)(a(i),i=1,2),(ia(2+i),i=1,3)

   !--set hollerith identification
   nwds=mult*nholl
   jrec=3
   do i=1,nwds
      read(hblank,'(a8)') ha(i)
   enddo
   write(nout)(ha(i),i=1,nholl)

   !--file data
   nwds=(mult+2)*nmat
   jrec=4
   do i=1,nmat
      read(hmat(i),'(a6)') ha(iscrh-1+i)
      ia(iscr-1+nmat*mult+i)=ntemp(i)
      ia(iscr-1+nmat*(mult+1)+i)=locm(i)
   enddo
   write(nout)(ha(i),i=1,nmat),(ia(mult*nmat+i),i=1,2*nmat)

   !--copy the material control and cross-section blocks
   !--from the scratch file to nout
   krec=0
   call repoz(-nscr)
   do i=1,nmat
      krec=krec+1
      nwr=1+ntemp(i)
      nwds=mult+ntemp(i)+3
      read(nscr) ha(1),(a(mult+j),j=1,nwr),ia(nwds-1),ia(nwds)
      jrec=jrec+1
      write(nout) ha(1),(a(mult+j),j=1,nwr),ia(nwds-1),ia(nwds)
      nreac=ia(iscr+mult+ntemp(i)+1)
      nener=ia(iscr+mult+ntemp(i)+2)
      nn=1+nreac*ntemp(i)
      nw=nener*nn
      nwds=(nblok/nn)*nn
      nb=nw/nwds
      if (nw.gt.nb*nwds) nb=nb+1
      do j=1,nb
         if (nw.lt.nwds) nwds=nw
         krec=krec+1
         read(nscr)(a(k),k=1,nwds)
         jrec=jrec+1
         write(nout)(a(k),k=1,nwds)
         nw=nw-nwds
      enddo
   enddo
   end subroutine resxsr

end module resxsm


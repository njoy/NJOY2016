module matxsm
   ! provides subroutine matxsr for NJOY2016
   use locale
   implicit none
   private
   public matxsr
   ! global variables

   ! equivalenced arrays for CCCC input and output
   integer::isiza=200000    ! container array size
   real(k4)::a(200000)      ! reals are 4-byte
   integer(k4)::ia(200000)  ! integers are 4-byte
   real(k8)::ha(100000)     ! Hollerith data are 8-byte
   character(8)::ta(100000) ! text equivalent to Hollerith
   integer::mult=2          ! used for counting 8-byte entries
   equivalence(a(1),ia(1),ha(1),ta(1))

   ! unit numbers
   integer::ngen1,ngen2,nmatx,nscrt1,nscrt2,nscrt3,nscrt4,&
     nscrt5,nscrt6,nscrt7,nscrt8,ngen3,ngen4,ngen5,ngen6,ngen7,ngen8

   ! pointers
   integer::next,icont,iholl,ifild,igrup

   ! file identification
   real(kr)::huse(2)
   character(8)::tuse(2)
   equivalence(huse(1),tuse(1))
   integer::ivers

   ! file control
   integer::ntype,npart,nholl,nmat,maxw,length

   ! file data
   character(8)::hprt(10),htype(10),hmatn(250)
   integer::ngrp(10),jinp(10),joutp(10),matno(250),matgg(250)

   ! group structure sizes
   integer::ning,noutg

   ! reaction names
   character(8)::hvps(3000),hmtx(3000)

   ! self-shielding parameters
   integer::iz,iglt

   ! matrix parameters
   integer::imcon,icdat,ijgll,imdat,ng2z,lord1,jconst
   integer::maxord,nsubmx
   integer::nritev,nritem
   integer::ntw

contains

   subroutine matxsr
   !-------------------------------------------------------------------
   !
   !  Produce MATXS interface file from NJOY intermediate cross
   !  section data from group or gaminr.
   !
   !  The matxs file uses a generalized, flexible format based on
   !  the CCCC interface conventions.  Working from groupr and/or
   !  gaminr output tapes, this module can process neutron,
   !  thermal, photon, and charged-particle data into a
   !  single output file.  This file can then be used by the
   !  TRANSX code to prepare data libraries for transport codes.
   !
   !  A MATXS file specification may be found following the
   !  input instructions.
   !
   !---input specifications (free format)---------------------------
   !
   ! card 1 units
   !   ngen1     input unit for data from groupr
   !   ngen2     input unit for data from gaminr
   !   nmatx     output unit for matxs
   !   ngen3     incident proton data from groupr (default=0)
   !   ngen4     incident deuteron data from groupr (default=0)
   !   ngen5     incident triton data from groupr (default=0)
   !   ngen6     incident he3 data from groupr (default=0)
   !   ngen7     incident alpha data from groupr (default=0)
   !   ngen8     photonuclear data from groupr (default=0)
   ! card 2 user identification
   !   ivers     file version number (default=0)
   !   huse      user id (up to 16 characters, delimited by *,
   !             ended by /) (default=blank)
   ! card 3 file control
   !   npart     number of particles for which group
   !                structures are given
   !   ntype     number of data types in set
   !   nholl     number of cards to be read for hollerith
   !             identification record.
   !   nmat      number of materials desired
   ! card 4 set hollerith identification
   !   hsetid    hollerith identification of set
   !             (each line can be up to 72 characters,
   !             delimited with *, ended by /)
   ! card 5 particle identifiers
   !   hpart     hollerith identifiers for particles
   !             (up to 8 characters each)
   ! card 6 energy groups
   !   ngrp      number of groups for each particle
   ! card 7 data type identifiers
   !   htype     hollerith identifiers for data types
   !             (up to 8 characters each)
   ! card 8 input particle ids
   !   jinp     input particle id for each data type
   ! card 9 output particle ids
   !   joutp    output particle id for each data type
   ! card 10 material data (one card per material)
   !   hmat      hollerith material identifier
   !             (up to 8 characters each)
   !   matno     integer material identifier
   !             (endf mat number)
   !   matgg     mat number for photoatomic data
   !             (default=100*(matno/100) as in endf-6)
   !
   !-------------------------------------------------------------------

!          Standardized CCCC format listing for MATXS file
!c
!c**********************************************************************
!c               proposed 09/09/77
!c                       (modified 09/80)
!c                       (nomenclature changed 06/88)
!c                       (modified for const sub-blocks 06/90)
!c                       (ordering changed 10/90)
!c     c                       (bcd format changed 12/21/91)
!c
!cf           matxs
!ce           material cross section file
!c
!cn                       this file contains cross section
!cn                       vectors and matrices for all
!cn                       particles, materials, and reactions;
!cn                       delayed neutron spectra by time group;
!cn                       and decay heat and photon spectra.
!c
!cn           formats given are for file exchange only
!c
!c**********************************************************************
!c
!c
!c----------------------------------------------------------------------
!cs          file structure
!cs
!cs              record type                       present if
!cs              ==============================    ===============
!cs              file identification                 always
!cs              file control                        always
!cs              set hollerith identification        always
!cs              file data                           always
!cs
!cs   *************(repeat for all particles)
!cs   *          group structures                    always
!cs   *************
!cs
!cs   *************(repeat for all materials)
!cs   *          material control                    always
!cs   *
!cs   * ***********(repeat for all submaterials)
!cs   * *        vector control                      n1db.gt.0
!cs   * *
!cs   * * *********(repeat for all vector blocks)
!cs   * * *      vector block                        n1db.gt.0
!cs   * * *********
!cs   * *
!cs   * * *********(repeat for all matrix blocks)
!cs   * * *      matrix control                      n2d.gt.0
!cs   * * *
!cs   * * * *******(repeat for all sub-blocks)
!cs   * * * *    matrix sub-block                    n2d.gt.0
!cs   * * * *******
!cs   * * *
!cs   * * *      constant sub-block                  jconst.gt.0
!cs   * * *
!cs   *************
!c
!c----------------------------------------------------------------------
!c
!c
!c----------------------------------------------------------------------
!cr           file identification
!c
!cl    hname,(huse(i),i=1,2),ivers
!c
!cw    1+3*mult
!c
!cb    format(4h 0v ,a8,1h*,2a8,1h*,i6)
!c
!cd    hname         hollerith file name  - matxs -  (a8)
!cd    huse          hollerith user identifiation    (a8)
!cd    ivers         file version number
!cd    mult          double precision parameter
!cd                       1- a8 word is single word
!cd                       2- a8 word is double precision word
!c
!c----------------------------------------------------------------------
!c
!c
!c----------------------------------------------------------------------
!cr           file control
!c
!cl    npart,ntype,nholl,nmat,maxw,length
!c
!cw    6
!c
!cb    format(6h 1d   ,6i6)
!c
!cd    npart       number of particles for which group
!cd                   structures are given
!cd    ntype       number of data types present in set
!cd    nholl       number of words in set hollerith
!cd                    identification record
!cd    nmat        number of materials on file
!cd    maxw        maximum record size for sub-blocking
!cd    length      length of file
!c
!c----------------------------------------------------------------------
!c
!c
!c----------------------------------------------------------------------
!cr           set hollerith identification
!c
!cl    (hsetid(i),i=1,nholl)
!c
!cw    nholl*mult
!c
!cb    format(4h 2d /(9a8))
!c
!cd    hsetid      hollerith identification of set (a8)
!cd                 (to be edited out 72 characters per line)
!c
!c----------------------------------------------------------------------
!c
!c
!c----------------------------------------------------------------------
!cr          file data
!c
!cl    (hprt(j),j=1,npart),(htype(k),k=1,ntype),(hmatn(i),i=1,nmat),
!cl   1(ngrp(j),j=1,npart),(jinp(k),k=1,ntype,(joutp(k),k=1,ntype),
!cl   2(nsubm(i)i=1,nmat),(locm(i),i=1,nmat)
!c
!cw    (npart+ntype+nmat)*mult+2*ntype+npart+2*nmat
!c
!cb    format(4h 3d ,4x,8a8/(9a8))    hprt,htype,hmatn
!cb    format(12i6)                  ngrp,jinp,joutp,nsubm,locm
!c
!cd    hprt(j)     hollerith identification for particle j
!cd                     n         neutron
!cd                     g         gamma
!cd                     p         proton
!cd                     d         deuteron
!cd                     t         triton
!cd                     h         he-3 nucleus
!cd                     a         alpha (he-4 nucleus)
!cd                     b         beta
!cd                     r         residual or recoil
!cd                               (heavier than alpha)
!cd    htype(k)     hollerith identification for data type k
!cd                     nscat     neutron scattering
!cd                     ng        neutron induced gamma production
!cd                     gscat     gamma scattering (atomic)
!cd                     gg        gamma scattering (photonuclear)
!cd                     pn        proton induced neutron production
!cd                       .          .
!cd                       .          .
!cd                       .          .
!cd                     dkn       delayed neutron data
!cd                     dkhg      decay heat and gamma data
!cd                     dkb       decay beta data
!cd    hmatn(i)    hollerith identification for material i
!cd    ngrp(j)      number of energy groups for particle j
!cd    jinp(k)     type of incident particle associated with
!cd                   data type k.  for dk data types, jinp is 0.
!cd    joutp(k)    type of outgoing particle associated with
!cd                   data type k
!cd    nsubm(i)    number of submaterials for material i
!cd    locm(i)     location of material i
!c
!c----------------------------------------------------------------------
!c
!c
!c----------------------------------------------------------------------
!cr          group structure
!c
!cl    (gpb(i),i=1,ngr),emin
!c
!cc    ngr=ngrp(j)
!c
!cw    ngrp(j)+1
!c
!cb    format(4h 4d ,8x,1p,5e12.5/(6e12.5))
!c
!cd    gpb(i)      maximum energy bound for group i for particle j
!cd    emin        minimum energy bound for particle j
!c
!c----------------------------------------------------------------------
!c
!c
!c----------------------------------------------------------------------
!cr          material control
!c
!cl    hmat,amass,(temp(i),sigz(i),itype(i),n1d(i),n2d(i),
!cl   1locs(i),i=1,nsubm)
!c
!cw    mult+1+6*nsubm
!c
!cb    format(4h 5d ,a8,1p,2e12.5/(2e12.5,5i6))
!c
!cd    hmat        hollerith material identifier
!cd    amass       atomic weight ratio
!cd    temp        ambient temperature or other parameters for
!cd                    submaterial i
!cd    sigz        dilution factor or other parameters for
!cd                    submaterial i
!cd    itype       data type for submaterial i
!cd    n1d         number of vectors for submaterial i
!cd    n2d         number of matrix blocks for submaterial i
!cd    locs        location of submaterial i
!c
!c----------------------------------------------------------------------
!c
!c
!c----------------------------------------------------------------------
!cr          vector control
!c
!cl    (hvps(i),i=1,n1d),(nfg(i),i=1,n1d),(nlg(i),i=1,n1d)
!c
!cw    (mult+2)*n1d
!c
!cb    format(4h 6d ,4x,8a8/(9a8))      hvps
!cb    format(12i6)                    iblk,nfg,nlg
!c
!cd    hvps(i)     hollerith identifier of vector
!cd                      nelas     neutron elastic scattering
!cd                      n2n       (n,2n)
!cd                      nnf       second chance fission
!cd                      gabs      gamma absorption
!cd                      p2n       protons in, 2 neutrons out
!cd                         .          .
!cd                         .          .
!cd                         .          .
!cd    nfg(i)      number of first group in band for vector i
!cd    nlg(i)      number of last group in band for vector i
!c
!c----------------------------------------------------------------------
!c
!c
!c----------------------------------------------------------------------
!cr          vector block
!c
!cl    (vps(i),i=1,kmax)
!c
!cc    kmax=sum over group band for each vector in block j
!c
!cw    kmax
!c
!cb    format(4h 7d ,8x,1p,5e12.5/(6e12.5))
!c
!cd    vps(i)      data for group bands for vectors in block j.
!cd                block size is determined by taking all the group
!cd                bands that have a total length less than or equal
!cd                to maxw.
!c
!c----------------------------------------------------------------------
!c
!c
!c----------------------------------------------------------------------
!cr        scattering matrix control
!c
!cl    hmtx,lord,jconst,
!cl   1(jband(l),l=1,noutg(k)),(ijj(l),l=1,noutg(k))
!c
!cw    mult+2+2*noutg(k)
!c
!cb    format(4h 8d ,4x,a8/(12i6))      hmtx,lord,jconst,
!cb                                     jband,ijj
!c
!cd    hmtx        hollerith identification of block
!cd    lord        number of orders present
!cd    jconst      number of groups with constant spectrum
!cd    jband(l)    bandwidth for group l
!cd    ijj(l)      lowest group in band for group l
!c
!c----------------------------------------------------------------------
!c
!c
!c----------------------------------------------------------------------
!cr          scattering sub-block
!c
!cl    (scat(k),k=1,kmax)
!c
!cc    kmax=lord times the sum over all jband in the group range of
!cc            this sub-block
!c
!cb    format(4h 9d ,8x,1p,5e12.5/(6e12.5))
!c
!cw    kmax
!c
!cd    scat(k)     matrix data given as bands of elements for initial
!cd                groups that lead to each final group.  the order
!cd                of the elements is as follows:  band for p0 of
!cd                group i, band for p1 of group i, ... , band for p0
!cd                of group i+1, band for p1 of group i+1, etc.  the
!cd                groups in each band are given in descending order.
!cd                the size of each sub-block is determined by the
!cd                total length of a group of bands that is less than
!cd                or equal to maxw.
!cd
!cd                if jconst.gt.0, the contributions from the jconst
!cd                low-energy groups are given separately.
!c
!c----------------------------------------------------------------------
!c
!c
!c----------------------------------------------------------------------
!cr          constant sub-block
!c
!cl    (spec(l),l=1,noutg(k)),(prod(l),l=l1,ning(k))
!c
!cc    l1=ning(k)-jconst+1
!c
!cw    noutg(k)+jconst
!c
!cb    format(4h10d ,8x,1p,5e12.5/(6e12.5))
!c
!cd    spec        normalized spectrum of final particles for initial
!cd                particles in groups l1 to ning(k)
!cd    prod        production cross section (e.g., nu*sigf) for
!cd                initial groups l1 through ning(k)
!cd
!cd         this option is normally used for the energy-independent
!cd         neutron and photon spectra from fission and radiative
!cd         capture usually seen at low energies.
!c
!c----------------------------------------------------------------------
!c

   !-------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso,nsyse
   use util   ! provides timer,openz,closz
   ! internals
   real(kr)::time
   character(16)::word

   !--initialize.
   call timer(time)
   write(nsyso,'(/&
     &'' matxsr...produce a matxs format output file'',&
     &25x,f8.1,''s'')') time
   write(nsyse,'(/'' matxsr...'',59x,f8.1,''s'')') time

   !--read users input/output unit numbers.
   ngen1=0
   ngen2=0
   nmatx=0
   ngen3=0
   ngen4=0
   ngen5=0
   ngen6=0
   ngen7=0
   ngen8=0
   read(nsysi,*) ngen1,ngen2,nmatx,ngen3,ngen4,ngen5,ngen6,ngen7,ngen8
   write(nsyso,'(/&
     &'' input gendf unit ..................... '',i10/&
     &'' input gamout unit .................... '',i10/&
     &'' output matxs unit .................... '',i10)')&
     ngen1,ngen2,nmatx
   if (ngen3.ne.0) write(nsyso,'(&
     &'' incident proton unit ................. '',i10/&
     &'' incident deuteron unit ............... '',i10/&
     &'' incident triton unit ................. '',i10/&
     &'' incident he3 unit .................... '',i10/&
     &'' incident alpha unit .................. '',i10/&
     &'' photonuclear unit .................... '',i10)')&
     ngen3,ngen4,ngen5,ngen6,ngen7,ngen8
   nscrt1=10
   nscrt2=11
   nscrt3=12
   nscrt4=13
   nscrt5=14
   nscrt6=15
   nscrt7=16
   nscrt8=17
   call openz(ngen1,0)
   call openz(ngen2,0)
   call openz(ngen3,0)
   call openz(ngen4,0)
   call openz(ngen5,0)
   call openz(ngen6,0)
   call openz(ngen7,0)
   call openz(ngen8,0)
   call openz(nmatx,1)
   call openz(-nscrt2,1)
   call openz(-nscrt3,1)
   call openz(-nscrt6,1)

   !--user input for file identification
   ivers=0
   word=' '
   read(nsysi,*) ivers,word
   read(word,'(a8,a8)') tuse(1),tuse(2)
   write(nsyso,'(&
     &'' file version number .................. '',i10/&
     &'' user id .............................. '',2a8)')&
     ivers,tuse

   !--generate matxs file
   call cmatxs

   !--finished
   call closz(nmatx)
   call closz(ngen1)
   call closz(ngen2)
   call closz(ngen3)
   call closz(ngen4)
   call closz(ngen5)
   call closz(ngen6)
   call closz(ngen7)
   call closz(ngen8)
   call closz(-nscrt2)
   call closz(-nscrt3)
   call closz(-nscrt6)
   call timer(time)
   write(nsyso,'(69x,f8.1,''s''/1x,7(''**********''),&
     &''*******'')') time
   return
   end subroutine matxsr

   subroutine cmatxs
   !-------------------------------------------------------------------
   ! Writes matxs interface file to nmatx.
   !------------------------------------------------------------------
   use util ! provides repoz
   ! internals
   integer::irec2,irec3,irec6,iscrt,iscrh,nwds,irec,i
   integer::ihholl,ndr,nir,ihfild,lout,loc,j,imatc,nsubj,ihmatc
   integer::ll,ij,itype,n1d,n2d,ivcon,nw,ihvcon,ivdat,l
   integer::nfg,nlg,nwmcr,nwdc,n,ihmcon,lord,m,jj,k
   character(6),parameter::hfile='matxs '
   nsubmx=100
   maxord=5
   maxw=5000

   !--read user input
   call ruinm

   !--initiate reformatting of data
   call mtxdat
   irec2=0
   irec3=0
   irec6=0
   call repoz(-nscrt2)
   call repoz(-nscrt3)
   call repoz(-nscrt6)

   !--point to scratch area
   iscrt=next
   iscrt=(iscrt/2)*2+1
   iscrh=(iscrt-1)/mult+1

   !--write file identification to nmatx
   nwds=3*mult+1
 ! read(hfile,'(a)') ha(iscrh)
   ha(iscrh+1)=huse(1)
   ha(iscrh+2)=huse(2)
   ia(iscrt+3*mult)=ivers
   irec=1
   if (nmatx.lt.0)&
     write(-nmatx)(ha(iscrh+i-1),i=1,3),ia(iscrt+3*mult)
   if (nmatx.gt.0) write(nmatx,'('' 0v '',a8,''*'',2a8,''*'',i6)')&
     hfile,tuse(1),tuse(2),ivers

   !--write file control record to nmatx
   nwds=6
   irec=2
   ia(icont+4)=maxw
   ia(icont+5)=length+5+npart
   if (nmatx.lt.0) write(-nmatx)(ia(icont+i-1),i=1,nwds)
   if (nmatx.gt.0) write(nmatx,'('' 1d   '',6i6)')&
     (ia(icont+i-1),i=1,nwds)

   !--write set hollerith id to nmatx
   nwds=mult*nholl
   irec=3
   ihholl=(iholl-1)/mult
   if (nmatx.lt.0) write(-nmatx) (ha(ihholl+i),i=1,nholl)
   if (nmatx.gt.0) write(nmatx,'('' 2d ''/(9a8))')&
     (ta(ihholl+i),i=1,nholl)

   !--write file data
   ndr=npart+ntype+nmat
   nir=npart+2*ntype+2*nmat
   nwds=ndr*mult+nir
   irec=4
   ihfild=(ifild-1)/mult
   if (nmatx.lt.0) write(-nmatx) (ha(ihfild+i),i=1,ndr),&
     (ia(ifild+ndr*mult+i-1),i=1,nir)
   if (nmatx.gt.0) then
      write(nmatx,'('' 3d '',4x,8a8/(9a8))') (ta(ihfild+i),i=1,ndr)
      write(nmatx,'(12i6)') (ia(ifild+ndr*mult+i-1),i=1,nir)
   endif

   !--write group structures
   lout=igrup
   do i=1,npart
      loc=ifild-1+(npart+ntype+nmat)*mult+i
      nwds=ia(loc)+1
      irec=irec+1
      if (nmatx.lt.0) write(-nmatx) (a(lout+j-1),j=1,nwds)
      if (nmatx.gt.0) write(nmatx,&
        '('' 4d '',1p,8x,5e12.5/(6e12.5))') (a(lout+j-1),j=1,nwds)
      lout=lout+nwds
   enddo

   !--loop over materials
   do j=1,nmat

      !--copy material control
      imatc=iscrt
      loc=ifild-1+(npart+ntype+nmat)*mult+npart+2*ntype+j
      nsubj=ia(loc)
      nwds=mult+1+6*nsubj
      irec6=irec6+1
      ihmatc=(imatc-1)/mult+1
      ll=imatc+mult
      read(nscrt6) ha(ihmatc),a(imatc+mult),&
        ((a(ll+6*(i-1)+k),k=1,2),(ia(ll+6*(i-1)+k),k=3,6),i=1,nsubj)
      irec=irec+1
      if (nmatx.lt.0) write(-nmatx) ha(ihmatc),a(imatc+mult),&
        ((a(ll+6*(i-1)+k),k=1,2),(ia(ll+6*(i-1)+k),k=3,6),i=1,nsubj)
      if (nmatx.gt.0) then
         write(nmatx,'('' 5d '',a8,1p,e12.5)')&
           ta(ihmatc),a(imatc+mult)
         do i=1,nsubj
            ll=imatc+mult+1+6*(i-1)
            write(nmatx,'(1p,2e12.5,4i6)') a(ll),a(ll+1),&
              ia(ll+2),ia(ll+3),ia(ll+4),ia(ll+5)
         enddo
      endif
      next=imatc+nwds+mult-1

      !--loop over all submaterials
      do ij=1,nsubj
         itype=ia(imatc+mult+6*(ij-1)+3)
         n1d=ia(imatc+mult+6*(ij-1)+4)
         n2d=ia(imatc+mult+6*(ij-1)+5)
         ning=ngrp(jinp(itype))
         noutg=ngrp(joutp(itype))
         if (n1d.ne.0) then

            !--copy vector control
            ivcon=next
            nwds=(mult+2)*n1d
            nw=2*n1d
            irec2=irec2+1
            ihvcon=(ivcon-1)/mult
            read(nscrt2)&
              (ha(ihvcon+i),i=1,n1d),(ia(ivcon+n1d*mult+i-1),i=1,nw)
            irec=irec+1
            if (nmatx.lt.0) then
               write(-nmatx)&
                 (ha(ihvcon+i),i=1,n1d),(ia(ivcon+n1d*mult+i-1),i=1,nw)
            else if (nmatx.gt.0) then
               write(nmatx,'('' 6d '',4x,8a8/(9a8))')&
                 (ta(ihvcon+i),i=1,n1d)
               write(nmatx,'(12i6)') (ia(ivcon+n1d*mult+i-1),i=1,nw)
            endif
            ivdat=ivcon+nwds

            !--copy vector partial blocks
            nwds=0
            do l=1,n1d
               loc=ivcon-1+mult*n1d+l
               nfg=ia(loc)
               nlg=ia(loc+n1d)
               nw=nlg-nfg+1
               if (nwds+nw.ge.maxw) then
                  irec3=irec3+1
                  read(nscrt3) (a(ivdat+i-1),i=1,nwds)
                  irec=irec+1
                  if (nmatx.lt.0) then
                     write(-nmatx) (a(ivdat+i-1),i=1,nwds)
                  else if (nmatx.gt.0) then
                     write(nmatx,&
                       '('' 7d '',8x,1p,5e12.5/(6e12.5))')&
                       (a(ivdat+i-1),i=1,nwds)
                  endif
                  nwds=0
               endif
               nwds=nwds+nw
            enddo
            if (nwds.gt.0) then
               irec3=irec3+1
               read(nscrt3) (a(ivdat+i-1),i=1,nwds)
               irec=irec+1
               if (nmatx.lt.0) then
                  write(-nmatx) (a(ivdat+i-1),i=1,nwds)
               else if (nmatx.gt.0) then
                  write(nmatx,&
                    '('' 7d '',8x,1p,5e12.5/(6e12.5))')&
                    (a(ivdat+i-1),i=1,nwds)
               endif
            endif
         endif

         !--loop over matrix blocks
         if (n2d.ne.0) then
            imcon=next
            nwmcr=mult+2+2*noutg
            icdat=imcon+nwmcr
            nwdc=noutg+ning
            imdat=icdat+nwdc
            do n=1,n2d

               !--copy matrix control record
               irec2=irec2+1
               ihmcon=(imcon-1)/mult+1
               nw=nwmcr-mult
               read(nscrt2) ha(ihmcon),(ia(imcon+mult+i-1),i=1,nw)
               irec=irec+1
               if (nmatx.lt.0) then
                  write(-nmatx) ha(ihmcon),(ia(imcon+mult+i-1),i=1,nw)
               else if (nmatx.gt.0) then
                  write(nmatx,'('' 8d '',4x,a8)') ta(ihmcon)
                  write(nmatx,'(12i6)') (ia(imcon+mult+i-1),i=1,nw)
               endif

               !--copy matrix sub-blocks
               loc=imcon+mult
               lord=ia(loc)
               loc=loc+1
               jconst=ia(loc)
               nwdc=noutg+jconst
               if (lord.ne.0) then
                  m=loc
                  nwds=0
                  do jj=1,noutg
                     l=m+jj
                     nw=ia(l)
                     if (nwds+nw*lord.ge.maxw.and.nwds.gt.0) then
                        irec3=irec3+1
                        read(nscrt3) (a(imdat+i-1),i=1,nwds)
                        irec=irec+1
                        if (nmatx.lt.0) then
                           write(-nmatx) (a(imdat+i-1),i=1,nwds)
                        else if (nmatx.gt.0) then
                           write(nmatx,&
                             '('' 9d '',8x,1p,5e12.5/(6e12.5))')&
                             (a(imdat+i-1),i=1,nwds)
                        endif
                        nwds=0
                     endif
                     nwds=nwds+nw*lord
                  enddo
                  if (nwds.gt.0) then
                     irec3=irec3+1
                     read(nscrt3) (a(imdat+i-1),i=1,nwds)
                     irec=irec+1
                     if (nmatx.lt.0) then
                        write(-nmatx) (a(imdat+i-1),i=1,nwds)
                     else if (nmatx.gt.0) then
                        write(nmatx,&
                          '('' 9d '',8x,1p,5e12.5/(6e12.5))')&
                          (a(imdat+i-1),i=1,nwds)
                     endif
                  endif
               endif
               if (jconst.gt.0) then
                  irec3=irec3+1
                  read(nscrt3) (a(icdat+i-1),i=1,nwdc)
                  irec=irec+1
                  if (nmatx.lt.0) then
                     write(-nmatx) (a(icdat+i-1),i=1,nwdc)
                  else if (nmatx.gt.0) then
                     write(nmatx,&
                       '(''10d '',8x,1p,5e12.5/(6e12.5))')&
                       (a(icdat+i-1),i=1,nwdc)
                  endif
               endif

            !--end of loop over matrix blocks
            enddo
         endif

      !--end of submaterial loop
      enddo

   !--end of material loop
   enddo
   return
   end subroutine cmatxs

   subroutine ruinm
   !-------------------------------------------------------------------
   ! Read user input.
   !-------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   ! internals
   integer::lim,i,j,ifih,nwds,loc
   character(72)::text

   !--read file control parameters
   read(nsysi,*) npart,ntype,nholl,nmat
   write(nsyso,'(&
     &'' npart ................................ '',i10/&
     &'' ntype ................................ '',i10/&
     &'' nholl ................................ '',i10/&
     &'' nmat ................................. '',i10/&
     &'' maxw ................................. '',i10)')&
     npart,ntype,nholl,nmat,maxw

   !--load file identification record
   lim=9
   icont=1
   next=icont+6
   ia(icont)=npart
   ia(icont+1)=ntype
   ia(icont+2)=nholl*lim
   ia(icont+3)=nmat
   ia(icont+4)=0
   ia(icont+5)=0

   !--read and load hollerith set identification
   write(nsyso,&
     '(/'' set identification''/'' ------------------'')')
   lim=9
   iholl=next
   nwds=nholl*mult*lim
   next=iholl+nwds
   do i=1,nholl
      text=' '
      read(nsysi,*) text
      loc=(iholl-1)/mult+lim*(i-1)
      read(text,'(9a8)') (ta(loc+j),j=1,lim)
      write(nsyso,'(6x,9a8)') (ta(j+loc),j=1,lim)
   enddo
   nholl=nholl*lim

   !--read particle identifiers
   read(nsysi,*) (hprt(i),i=1,npart)
   write(nsyso,'(/&
     &'' names of particles ................... '',&
     &6x,5(a8,2x))')&
     (hprt(i),i=1,npart)

   !--read number of groups for each particle
   read(nsysi,*) (ngrp(i),i=1,npart)
   write(nsyso,'(&
     &'' no. groups for each particle ......... '',&
     &2x,12i10)')&
     (ngrp(i),i=1,npart)

   !--read hollerith names for data types
   read(nsysi,*) (htype(i),i=1,ntype)
   write(nsyso,'(&
     &'' names of data types .................. '',&
     &6x,12(a8,2x))')&
     (htype(i),i=1,ntype)

   !--read input and output particles for data types
   read(nsysi,*) (jinp(i),i=1,ntype)
   read(nsysi,*) (joutp(i),i=1,ntype)
   write(nsyso,'(&
     &'' input particles ...................... '',2x,12i10)')&
     &(jinp(i),i=1,ntype)
   write(nsyso,'(&
     &'' output particles ..................... '',2x,12i10)')&
     (joutp(i),i=1,ntype)

   !--read in material names and id numbers
   do i=1,nmat
      matgg(i)=0
      read(nsysi,*) hmatn(i),matno(i),matgg(i)
      if (matgg(i).eq.0) matgg(i)=100*(matno(i)/100)
   enddo

   !--set up file data record
   ifild=next
   nwds=(npart+ntype+nmat)*mult+2*ntype+2*nmat+npart
   nwds=((nwds+1)/2)*2
   next=ifild+nwds
   ifih=(ifild-1)/mult+1
   loc=ifih-1
   do i=1,npart
      read(hprt(i),'(a8)') ta(loc+i)
   enddo
   loc=loc+npart
   do i=1,ntype
      read(htype(i),'(a8)') ta(loc+i)
   enddo
   loc=loc+ntype
   do i=1,nmat
      read(hmatn(i),'(a8)') ta(loc+i)
   enddo
   loc=ifild-1+(npart+ntype+nmat)*mult
   do i=1,npart
      ia(loc+i)=ngrp(i)
   enddo
   do i=1,ntype
      ia(loc+npart+i)=jinp(i)
      ia(loc+npart+ntype+i)=joutp(i)
   enddo
   return
   end subroutine ruinm

   subroutine mtxdat
   !-------------------------------------------------------------------
   ! Processes matxs data and writes to scratch tapes.
   !-------------------------------------------------------------------
   use endf   ! provides endf routines and variables
   use util   ! provides error,repoz,openz,timer,closz
   use mainio ! provides nsyso
   ! internals
   integer::ngpb,i,nin,nb,nw,ng1,ng2,iqd,iqg,ndex
   integer::ip,imatc,loca,irec2,irec3,irec6,im,k,imloc
   integer::imatch,isubm,locs,imat,n2d,itemp,isigz
   integer::nsigz,ll,lz,nwds,nz,nqd,nqg,nwmc,j
   integer::ip1,ip2,nscr,nref,mfv,mfm,nsz,jmat
   integer::nnmat,iref,mz,nl,ng,max,ig,loc,n1d,ihmatc
   real(kr)::temp,awr,sec
   character(60)::strng
   integer,parameter::nbmax=2000
   real(kr)::b(nbmax)
   real(kr)::sigz(15)
   integer::noned(15),ntwod(15)
   character(8)::htyp
   character(1)::hp
   character(8),parameter::hnt='n     '
   character(8),parameter::hgm='g     '
   character(8),parameter::hprot='p     '
   character(8),parameter::hdeut='d     '
   character(8),parameter::htrit='t     '
   character(8),parameter::hhe3='h     '
   character(8),parameter::halph='a     '
   character(8),parameter::hgsct='gscat '
   character(8),parameter::hgg='gg '
   character(8),parameter::hnthr='ntherm'
   integer,parameter::mtt1=221
   integer,parameter::mtt2=250
   real(kr),parameter::eps=1.0e-06_kr
   lz=6

   !--assign storage for group structure records
   ngpb=0
   do i=1,npart
      ngpb=ngrp(i)+1+ngpb
   enddo
   nwds=((ngpb-1)/2+1)*2
   igrup=next
   next=igrup+nwds

   !--choose input file to read group structures
   !--for now, all particles are assumed to use the neutron structur
   !--this is a limitation of the current implementation in groupr
   nin=ngen1
   if (nin.eq.0) nin=ngen3
   if (nin.eq.0) nin=ngen4
   if (nin.eq.0) nin=ngen5
   if (nin.eq.0) nin=ngen6
   if (nin.eq.0) nin=ngen7
   if (nin.eq.0) nin=ngen8
   if (nin.eq.0) nin=ngen2
   if (nin.eq.0) call error('mtxdat','input error (nin=0)',' ')

   !--read head record and store groups bounds
   call repoz(nin)
   call tpidio(nin,0,0,b(1),nb,nw)
   mfh=0
   do while (mfh.eq.0.or.mth.eq.0)
      call contio(nin,0,0,b(1),nb,nw)
   enddo
   nz=l2h
   ntw=n2h
   call lst1io(nin,0,0,b(1),nb,nw,nbmax)
   ng1=l1h
   ng2=l2h
   iqd=1+6+ntw+nz
   iqg=iqd+ng1+1
   ndex=igrup
   nqd=ng1+1
   nqg=ng2+1
   do ip=1,npart
      if (ip.eq.1.or.hprt(ip).ne.hgm) then
         do i=1,nqd
            a(ndex+nqd-i)=b(iqd-1+i)
         enddo
         ndex=ndex+nqd
      else
         do i=1,nqg
            a(ndex+nqg-i)=b(iqg-1+i)
         enddo
         ndex=ndex+nqg
      endif
   enddo

   !--material loop
   nwmc=2*mult+6*nsubmx
   nwmc=(((nwmc-1)/2)+1)*2
   imatc=next
   next=imatc+nwmc
   loca=0
   length=0
   write(nsyso,'(/)')
   irec2=0
   irec3=0
   irec6=0
   do 700 im=1,nmat
   do k=1,nwmc
      a(k-1+imatc)=0.
   enddo
   imatch=(imatc-1)/mult+1
   read(hmatn(im),'(a8)') ta(imatch)
   isubm=0
   locs=0

   !--data type loop
   do 600 i=1,ntype
   imat=matno(im)
   if (htype(i).eq.hgsct) then
      imat=matgg(im)
      if (hmatn(im)(1:2).eq.'fm') imat=9920
   endif
   ning=ngrp(jinp(i))
   noutg=ngrp(joutp(i))

   !--choose input unit and data files for this data type
   htyp=htype(i)
   ip1=jinp(i)
   ip2=joutp(i)
   if (hprt(ip1).eq.hnt) nin=ngen1
   if (hprt(ip1).eq.hgm) nin=ngen8
   if (hprt(ip1).eq.hgm.and.htyp.eq.hgsct) nin=ngen2
   if (hprt(ip1).eq.hprot) nin=ngen3
   if (hprt(ip1).eq.hdeut) nin=ngen4
   if (hprt(ip1).eq.htrit) nin=ngen5
   if (hprt(ip1).eq.hhe3) nin=ngen6
   if (hprt(ip1).eq.halph) nin=ngen7
   call repoz(nin)
   nscr=nscrt1
   if (nin.lt.0) nscr=-nscrt1
   nref=nscrt5
   if (nin.lt.0) nref=-nscrt5
   call openz(nscr,1)
   call openz(nref,1)
   mfv=3
   if (ip1.ne.ip2) mfv=0
   if (hprt(ip2).eq.hgm.and.htyp.ne.hgg) mfv=13
   if (hprt(ip1).eq.hgm.and.htyp.eq.hgsct) mfv=23
   if (mfv.ne.0) mfm=mfv+3
   if (hprt(ip2).eq.hnt) mfm=6
   if (hprt(ip2).eq.hprot) mfm=21
   if (hprt(ip2).eq.hdeut) mfm=22
   if (hprt(ip2).eq.htrit) mfm=23
   if (hprt(ip2).eq.hhe3) mfm=24
   if (hprt(ip2).eq.halph) mfm=25
   if (hprt(ip2).eq.hgm.and.htyp.eq.hgg) mfm=16
   hp=hprt(ip1)(1:1)

   !--find the material
   call tpidio(nin,0,0,b(1),nb,nw)
   nsz=0
   jmat=0
   nnmat=0
   isigz=0
   nsigz=0
   call findg(imat,1,451,nin)

   !--read head record for material
   call contio(nin,0,0,b(1),nb,nw)
   ll=1+6
   if (i.eq.1) awr=c2h
   nz=l2h
   ntw=n2h
   nsh=1
   call lst1io(nin,0,0,b(ll),nb,nw,nbmax)
   temp=c1h
   nnmat=nnmat+1
   nsigz=nz
   do k=1,nsigz
      sigz(k)=b(k+lz+ntw-1+ll)
   enddo
   isigz=1
   if (imat.eq.jmat) go to 300
   itemp=1
   jmat=imat
   go to 320
  300 continue
   call contio(nin,0,0,b(1),nb,nw)
   if (math.ne.imat) go to 475
   itemp=itemp+1
   isigz=1
   if (mfh.ne.1.or.mth.ne.451) then
      call listio(nin,0,0,b(ll),nb,nw)
      temp=c1h
      do while (nb.ne.0)
         call moreio(nin,0,0,b(ll),nb,nw)
      enddo
   else
      ntw=n2h
      call lst1io(nin,0,0,b(ll),nb,nw,nbmax)
      temp=c1h
   endif

   !--copy this material to scratch files
  320 continue
   call repoz(nscr)
   iref=0
   if (itemp.eq.1.and.isigz.eq.1) iref=nref
   call repoz(iref)
   nsh=0
   nsc=0
   call amend(nscr,iref)
   iglt=1000
   do mz=1,nsigz
      noned(mz)=0
      ntwod(mz)=0
   enddo
  340 continue
   call contio(nin,0,0,b(1),nb,nw)
   if (math.le.0) go to 430
   if (mfh.eq.0.or.mth.eq.0) go to 340
   nl=l1h
   nz=l2h
   ng=n2h
   max=6+nl*nz*(ng+1)
   if (max.gt.nsz) nsz=max
   call listio(nin,0,0,b(ll),nb,nw)
   ng2=l1h
   ig=n2h
   if (ng2.ne.2.or.ig.ne.ng) go to 360
   loc=ll+6+nl*nz
   if (abs(b(loc)).ge.eps) go to 360
   go to 420
  360 continue
   if (ig.gt.1.and.ig.lt.iglt) iglt=ig
   if (htyp.eq.hnthr.and.(mth.lt.mtt1.or.mth.gt.mtt2)) go to 420
   if (htyp.ne.hnthr.and.(mth.ge.mtt1.and.mth.le.mtt2)) go to 420
   if (htyp.eq.hnthr.and.mfh.ge.90) go to 420
   if (mfv.eq.3.and.mfh.gt.90) go to 365
   if (mfh.ne.mfv.and.mfh.ne.(mfv+2).and.&
     mfh.ne.mfm.and.mfh.ne.(mfm+1)) go to 420
  365 continue
   do 410 mz=1,nsigz
   if (nz.lt.mz) go to 410
   if (mfh.gt.90) go to 380
   if (mfh.ne.mfv) go to 390
   if (mth.ne.1) go to 370
   ! flux and ntot for all l-orders
   noned(mz)=noned(mz)+2*nl
   go to 410
   ! flux and gtot
  370 continue
   if (mth.ne.501) go to 375
   noned(mz)=noned(mz)+2
   go to 410
  375 continue
   if (mfh.gt.90) go to 380
   if (mfh.ne.2.or.mfm.lt.21) go to 380
   noned(mz)=noned(mz)+3
   go to 410
   ! other vector partials
  380 continue
   noned(mz)=noned(mz)+1
   go to 410
   ! search for delayed fission neutrons
  390 continue
   if (mfh.ne.(mfv+2)) go to 400
   if (mth.eq.455) noned(mz)=noned(mz)+1
   go to 410
   ! matrices (2-d arrays)
  400 continue
   if (mfh.eq.mfm) ntwod(mz)=ntwod(mz)+1
   if (mfh.eq.17.and.mfm.eq.16) ntwod(mz)=ntwod(mz)+1
  410 continue
   call contio(0,nscr,iref,b(1),nb,nw)
   call listio(0,nscr,iref,b(ll),nb,nw)
   call tosend(nin,nscr,iref,b(1))
   go to 340
  420 continue
   call tosend(nin,0,0,b(1))
   go to 340
  430 continue
   call amend(nscr,iref)
   call atend(nscr,iref)
   if (noned(1).gt.0.or.ntwod(1).gt.0) go to 435

   !--no data found
   write(strng,&
     '(''no data for itype='',i2,'' mat='',i4,'' itemp='',i1)')&
     i,imat,itemp
   call mess('mtxdat',strng,' ')

   !--label what we are now doing
  435 continue
   if (itemp.eq.1) then
      call timer(sec)
      write(nsyso,'(10x,''doing '',a8,'' for '',a8,'' ... '',&
        &f7.1,''s'')') hmatn(im),htype(i),sec
   endif

   !--loop over sigzeroes for this mat and temperature
   do while (isigz.le.nsigz)
      call repoz(nscr)
      iref=nref
      if (isigz.eq.1.and.itemp.eq.1) iref=0
      call repoz(iref)
      iz=isigz
      n1d=noned(isigz)
      nritev=0
      n2d=ntwod(isigz)
      nritem=0
      if (n1d.ne.0.or.n2d.ne.0) then

         !--produce vector cross sections
         if (n1d.ne.0)&
           call vector(mfv,hp,n1d,nscr,iref,irec2,irec3)

         !--produce matrix cross sections
         if (n2d.ne.0)&
           call matrix(irec3,irec2,mfm,hp,n2d,nscr,nscrt2,iref)

         !--enter submaterial data into material control
         isubm=isubm+1
         if (isubm.gt.nsubmx)&
            call error('mtxdat','too many submaterials',' ')
         imloc=imatc+mult+6*(isubm-1)
         a(imloc+1)=temp
         a(imloc+2)=sigz(isigz)
         ia(imloc+3)=i
         ia(imloc+4)=n1d
         ia(imloc+5)=n2d
         ia(imloc+6)=locs
         locs=locs+nritev+nritem
      endif
      isigz=isigz+1
   enddo
   go to 300

   !--close the scratch files
  475 continue
   call closz(nref)
   call closz(nscr)

   !--continue data type loop
  600 continue

   !--write material control and continue material loop
   a(imatc+mult)=awr
   nw=mult+1+6*isubm
   irec6=irec6+1
   ihmatc=(imatc-1)/mult+1
   ll=imatc+mult
   write(nscrt6) ha(ihmatc),a(imatc+mult),&
     ((a(ll+6*(i-1)+j),j=1,2),(ia(ll+6*(i-1)+j),j=3,6),i=1,isubm)
   loc=ifild-1+(npart+ntype+nmat)*mult+npart+2*ntype+im
   ia(loc)=isubm
   loc=loc+nmat
   ia(loc)=loca
   loca=loca+1+locs
   length=length+1+locs

   !--continue material loop
  700 continue

   !--finished
   next=imatc
   return
   end subroutine mtxdat

   subroutine findg(mat,mf,mt,ntape)
   !-------------------------------------------------------------------
   ! Search a gendf tape for the desired section.
   !-------------------------------------------------------------------
   use util ! provides error,repoz,skiprz
   ! externals
   integer::mat,mf,mt,ntape
   ! internals
   integer::itend,nt,ifirst,mat1,mf1,mt1,math,mfh,mth
   character(60)::strng

   if (mat.le.0.or.mf.le.0.or.mt.le.0) call error('findg',&
     'mat or mf or mt le 0 not allowed',' ')
   itend=0
   nt=iabs(ntape)
   ifirst=0
   mat1=-10
   mf1=-10
   mt1=-10
  110 continue
   if (ntape.ge.0) then
      read(ntape,'(66x,i4,i2,i3)',end=160) math,mfh,mth
   else
      read(nt,end=160) math,mfh,mth
   endif
   if (ifirst.gt.0) go to 150
   if (math.eq.0.and.mfh.eq.0.and.mth.eq.0) go to 110
   mat1=math
   mf1=mfh
   mt1=mth
   ifirst=1
  150 continue
   if (math.eq.mat) go to 170
   if (math.eq.-1) go to 160
   if (itend.eq.0) go to 110
   if (math.eq.mat1.and.mfh.eq.mf1.and.mth.eq.mt1) go to 165
   go to 110
   ! at end of tape.  rewind and try beginning
  160 continue
   itend=itend+1
   if (itend.gt.1) go to 165
   call repoz(ntape)
   go to 110
   ! desired section not on tape
  165 continue
   write(strng,'(''mat='',i4,'' mf='',i2,'' mt='',i3,&
     &'' not on tape '',i3)') mat,mf,mt,ntape
   call error('findg',strng,' ')
   ! desired mat found.  check for mf and mt
  170 continue
   if (mfh.ne.mf.or.mth.ne.mt) go to 110
   ! desired section found.  backspace one record
   call skiprz(ntape,-1)
   return
   end subroutine findg

   subroutine hname(hreact,hp,mt,lr,izam)
   !-------------------------------------------------------------------
   ! Assigns hollerith names to reactions.  If lr flag is set,
   ! an additional identifier is placed in the last three
   ! (or six) characters.
   !-------------------------------------------------------------------
   use endf ! provides iverf
   ! externals
   integer::mt,lr,izam
   character::hreact*(*),hp*(*)
   ! internals
   integer::i,jmt,jlr,ii,mm,jj
   character(6)::strng,temp
   integer,parameter::num1=33
   integer,dimension(33),parameter::i1=(/&
     2,3,4,5,6,7,8,9,16,17,18,19,20,21,22,23,24,25,&
     26,28,29,30,32,33,34,35,36,37,38,41,42,44,45/)
   character(5),dimension(33),parameter::h1=(/&
     'elas ','nonel','inel ','x    ','2n1  ','2n2  ','2n3  ','2n4  ',&
     '2n   ','3n   ','ftot ','f    ','nf   ','2nf  ','na   ','n3a  ',&
     '2na  ','3na  ','2ni  ','np   ','n2a  ','2n2a ','nd   ','nt   ',&
     'nh   ','nd2a ','nt2a ','4n   ','3nf  ','2np  ','3np  ','n2p  ',&
     'npa  '/)
   integer,parameter::num2=16
   integer,dimension(16),parameter::i2=(/&
     101,102,103,104,105,106,107,108,109,111,112,113,114,115,116,117/)
   character(5),dimension(16),parameter::h2=(/&
     'abs  ','g    ','p    ','d    ','t    ','h    ','a    ',&
     '2a   ','h3a  ','2p   ','pa   ','t2a  ','d2a  ','pd   ','pt   ',&
     'da   '/)
   character(5),dimension(7),parameter::h3=(/&
     '.neut','.gam ','.h1  ','.h2  ','.h3  ','.he3 ','.he4 '/)
   integer,parameter::num4=246
   character(6),dimension(26),parameter::h4=(/&
     'free  ','hh2o  ','poly  ','poly$ ','hzrh  ','hzrh$ ',&
     'benz  ','dd2o  ','graph ','graph$','be    ','be$   ',&
     'bebeo ','bebeo$','zrzrh ','zrzrh$','obeo  ','obeo$ ',&
     'ouo2  ','ouo2$ ','uuo2  ','uuo2$ ','al    ','al$   ',&
     'fe    ','fe$   '/)
   integer,parameter::num5=10
   integer,dimension(10),parameter::i5=(/&
     251,252,253,259,301,443,444,452,455,456/)
   character(6),dimension(10),parameter::h5=(/&
     'mubar ','xi    ','gamma ','invel ','heat  ','kerma ','dame  ',&
     'nutot ','nudel ','nupmt '/)
   integer,parameter::num6=10
   integer,dimension(10),parameter::i6=(/&
     500,502,504,515,516,517,522,525,602,621/)
   character(6),dimension(10),parameter::h6=(/&
     'stop  ','gcoh  ','ginch ','gpaire','gpair ','gpairn',&
     'gabs  ','gheat ','gabs  ','gheat '/)
   integer,parameter::num7=19
   integer,dimension(19),parameter::ilr=(/&
     31,16,17,22,23,24,25,26,28,29,30,32,33,34,35,36,37,39,40/)
   character(3),dimension(19),parameter::hlr=(/&
     'g  ','n  ','nn ','a  ','aaa','na ','nna','ni ','p  ','aa ',&
     'naa','d  ','t  ','h  ','daa','taa','nnn','ic ','ee '/)

!!!!!!!!!!!!!!!!!!!!!!
   iverf=6
!!!!!!!!!!!!!!!!!!!!!!
   !--determine name from mt number or radionuclide
   if (izam.gt.0) go to 340

   !--neutron-emitting reactions
   if (mt.gt.49) go to 130
   do 110 i=1,num1
   jmt=i
   if (mt.eq.i1(i)) go to 120
  110 continue
   go to 300
  120 continue
   write(strng,'(a1,a5)') hp,h1(jmt)
   go to 350

   !--neutron-emitting levels
  130 continue
   if (mt.gt.91) go to 150
   if (mt.eq.91) go to 140
   if (mt.lt.60) write(strng,'(''n0'',i1)') mt-50
   if (mt.ge.60) write(strng,'(''n'',i2)') mt-50
   go to 310
  140 continue
   write(strng,'(''ncn'')')
   go to 310

   !--particle-emitting reactions
  150 continue
   if (mt.gt.150) go to 180
   do 160 i=1,num2
   jmt=i
   if (mt.eq.i2(i)) go to 170
  160 continue
   go to 300
  170 continue
   write(strng,'(a1,a5)') hp,h2(jmt)
   go to 350

   !--gas production
  180 continue
   if (mt.gt.207) go to 190
   if (mt.lt.201) go to 300
   write(strng,'(a6)') h3(mt-200)
   go to 350

   !--special njoy thermal mt numbers
  190 continue
   if (mt.gt.num4) go to 200
   if (mt.lt.221) go to 300
   write(strng,'(a6)') h4(mt-220)
   go to 350

   !--miscellaneous quantities
  200 continue
   if (mt.gt.499) go to 250
   do 210 i=1,num5
   jmt=i
   if (mt.eq.i5(i)) go to 220
  210 continue
   go to 300
  220 continue
   write(strng,'(a6)') h5(jmt)
   go to 350

   !--photon interactions
  250 continue
   if (iverf.le.5.and.mt.gt.699) go to 280
   if (iverf.ge.6.and.mt.gt.599) go to 290
   do 260 i=1,num6
   jmt=i
   if (mt.eq.i6(i)) go to 270
  260 continue
   go to 300
  270 continue
   write(strng,'(a6)') h6(jmt)
   go to 350

   !--particle-emitting levels (endf-5)
  280 continue
   if (mt.gt.799) go to 300
   if (mt.ge.700.and.mt.lt.710) write(strng,'(''p0'',i1)') mt-700
   if (mt.ge.710.and.mt.lt.719) write(strng,'(''p'',i2)') mt-700
   if (mt.eq.719) write(strng,'(''pcn'')')
   if (mt.ge.720.and.mt.lt.730) write(strng,'(''d0'',i1)') mt-720
   if (mt.ge.730.and.mt.lt.739) write(strng,'(''d'',i2)') mt-720
   if (mt.eq.739) write(strng,'(''dcn'')')
   if (mt.ge.740.and.mt.lt.750) write(strng,'(''t0'',i1)') mt-740
   if (mt.ge.750.and.mt.lt.759) write(strng,'(''t'',i2)') mt-740
   if (mt.eq.759) write(strng,'(''tcd'')')
   if (mt.ge.760.and.mt.lt.770) write(strng,'(''h0'',i1)') mt-760
   if (mt.ge.770.and.mt.lt.779) write(strng,'(''h'',i2)') mt-760
   if (mt.eq.779) write(strng,'(''hcn'')')
   if (mt.ge.780.and.mt.lt.790) write(strng,'(''a0'',i1)') mt-780
   if (mt.ge.790.and.mt.lt.799) write(strng,'(''a'',i2)') mt-780
   if (mt.eq.799) write(strng,'(''acn'')')
   go to 310

   !--particle-emitting levels (endf-6)
  290 continue
   if (mt.gt.891) go to 300
   if (mt.ge.600.and.mt.lt.610) write(strng,'(''p0'',i1)') mt-600
   if (mt.ge.610.and.mt.lt.649) write(strng,'(''p'',i2)') mt-600
   if (mt.eq.649) write(strng,'(''pcn'')')
   if (mt.ge.650.and.mt.lt.660) write(strng,'(''d0'',i1)') mt-650
   if (mt.ge.660.and.mt.lt.699) write(strng,'(''d'',i2)') mt-650
   if (mt.eq.699) write(strng,'(''dcn'')')
   if (mt.ge.700.and.mt.lt.710) write(strng,'(''t0'',i1)') mt-700
   if (mt.ge.710.and.mt.lt.749) write(strng,'(''t'',i2)') mt-700
   if (mt.eq.749) write(strng,'(''tcd'')')
   if (mt.ge.750.and.mt.lt.760) write(strng,'(''h0'',i1)') mt-750
   if (mt.ge.760.and.mt.lt.799) write(strng,'(''h'',i2)') mt-750
   if (mt.eq.799) write(strng,'(''hcn'')')
   if (mt.ge.800.and.mt.lt.810) write(strng,'(''a0'',i1)') mt-800
   if (mt.ge.810.and.mt.lt.849) write(strng,'(''a'',i2)') mt-800
   if (mt.eq.849) write(strng,'(''acn'')')
   if (mt.ge.875.and.mt.lt.885) write(strng,'(''2n0'',i1)') mt-875
   if (mt.ge.885.and.mt.lt.891) write(strng,'(''2n'',i2)') mt-875
   if (mt.eq.891) write(strng,'(''2ncn'')')
   go to 310

   !--encode mt number into name
  300 write(strng,'(''mt'',i3)') mt
   go to 350

   !--encode lr flag into name
  310 if (lr.eq.0) go to 350
   temp=strng(1:3)
   do 320 i=1,num7
   jlr=i
   if (lr.eq.ilr(i)) go to 330
  320 continue
   go to 350
  330 write(strng,'(a3,a3)') temp,hlr(jlr)
   go to 350

   !--radionuclide production
  340 ii=izam/10
   mm=izam-10*ii
   jj=ii+100*mm
   if (mt.eq.102.and.jj.ge.10000) write(strng,'(''c'',i5)') jj
   if (mt.eq.102.and.jj.lt.10000) write(strng,'(''c0'',i4)') jj
   if (mt.ne.102.and.jj.ge.10000) write(strng,'(''r'',i5)') jj
   if (mt.ne.102.and.jj.lt.10000) write(strng,'(''r0'',i4)') jj

   !--name is finished
  350 continue
   read(strng,'(a6)') hreact
   return
   end subroutine hname

   subroutine vector(mfd,hp,n1d,nscr,iref,irec2,irec3)
   !-------------------------------------------------------------------
   ! Loads vector control and vector cross sections onto
   ! scratch files.
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use util ! provides repoz,error
   ! externals
   integer::mfd,n1d,nscr,iref,irec2,irec3
   character::hp*(*)
   ! internals
   integer::k004,k016,k103,k104,k105,k106,k107,lz,n1i,nwds
   integer::ivcon,ivdat,i,ig,igr,igw,ig2lo,idref,il
   integer::izam,nl,nz,matr,mfr,mtr,loc,jg2,lout,lin
   integer::j,k,mdex,ivnow,nb,nw,ltot,ng2,n1j,nmove
   integer::ndex,nfg,nlg,nxsd,n,ihvcon
   integer,parameter::maxb=30000
   real(kr)::b(maxb)
   character(8)::hchid='chid'
   character(8)::hgtot='gtot0'
   character(8)::hgflx='gwt0'
   character(8)::hnt='n'
   character(8)::hga='g'
   real(kr),parameter::stol=1.e-3_kr
   real(kr),parameter::small=1.e-30_kr

   !--initialize.
   k004=0
   k016=0
   k103=0
   k104=0
   k105=0
   k106=0
   k107=0
   lz=6
   n1i=0
   nwds=(3+mult)*n1d
   nwds=((nwds-1)/2+1)*2
   ivcon=next
   next=ivcon+nwds
   nwds=ning*n1d
   ivdat=next
   next=ivdat+nwds
   do i=ivcon,isiza
      a(i)=0.0
   enddo

   !--read input data
   call repoz(nscr)
   call repoz(iref)
   call contio(nscr,0,0,b(1),nb,nw)
   if (iref.ne.0) call contio(iref,0,0,b(1),nb,nw)
   ig=ning

   !--set up next reaction
  115 continue
   if (ig.lt.ning) go to 155
   call contio(nscr,0,0,b(1),nb,nw)
   if (math.le.0) go to 400
   if (mfh.eq.0.or.mth.eq.0) go to 115
   if (mfh.eq.12) mfh=13
   if (mfh.le.1) go to 145
   izam=0
   if (mfh.gt.90) izam=nint(c2h)
   nl=l1h
   nz=l2h
   if (nz.lt.iz) go to 145
   ! delayed fission chi
   if (mfh.ne.(mfd+2)) go to 120
   n1i=n1i+1
   hvps(n1i)=hchid
   go to 140
  120 continue
   if (mfh.ne.mfd.and.mfh.le.90) go to 145
   if (mth.ne.1.and.mth.ne.501) go to 134
   ! flux and total cross section
   if (mth.eq.501) go to 130
   ltot=l1h
   do i=1,ltot
      write(hvps(i),'(a1,''wt'',i1)') hp,i-1
      write(hvps(i+ltot),'(a1,''tot'',i1)') hp,i-1
   enddo
   n1i=n1i+2*ltot
   go to 140
  130 continue
   hvps(1)=hgflx
   hvps(2)=hgtot
   n1i=n1i+2
   go to 140
  134 continue
   if (mth.eq.2.and.hp.ne.hnt) then
      write(hvps(1),'(a1,''wt0'')') hp
      write(hvps(2),'(a1,''tot0'')') hp
      n1i=n1i+2
   endif
   ! other vector partials (including delayed nu)
   n1i=n1i+1
   ! find the hollerith identifier
   call hname(hvps(n1i),hp,mth,n1h,izam)
  140 continue
   if (iref.eq.0) go to 150
   matr=math
   mfr=mfh
   mtr=mth
   call findg(matr,mfr,mtr,iref)
   call contio(iref,0,0,b(1),nb,nw)
   go to 150
  145 continue
   call tosend(nscr,0,0,b(1))
   go to 115
  150 continue
   ig=0
   igr=0
   igw=0

   !--read data for this group
  155 continue
   if (iref.eq.0.or.ig.le.igr) then
      call listio(nscr,0,0,b(1),nb,nw)
      if (mfh.eq.12) mfh=13
      ng2=l1h
      ig2lo=l2h
      ig=n2h
      loc=1+nw
      do while (nb.ne.0)
         if (loc+302.gt.maxb) call error('vector',&
           'exceeded input data array size',' ')
         call moreio(nscr,0,0,b(loc),nb,nw)
         if (mfh.eq.12) mfh=13
         loc=loc+nw
      enddo
   endif
   if (iref.ne.0.and.igr.le.ig) then
      idref=loc
      call listio(iref,0,0,b(idref),nb,nw)
      igr=n2h
      loc=loc+nw
      do while (nb.ne.0)
         if (loc+302.gt.maxb) call error('vector',&
           'exceeded input data array size',' ')
         call moreio(iref,0,0,b(loc),nb,nw)
         loc=loc+nw
      enddo
   endif
   if (mfh.eq.(mfd+2)) go to 200
   go to 300

   !--delayed fission spectrum(chi)
  200 continue
   ng2=l1h
   ig2lo=l2h
   do i=2,ng2
      jg2=noutg+1-ig2lo-i+2
      lout=ivdat-1+(n1i-1)*ning+jg2
      lin=lz+nl*(i-1)
      do j=1,nl
         a(lout)=b(j+lin)+a(lout)
      enddo
   enddo
   go to 115

   !--all other reactions
  300 continue
   if (mth.eq.1.or.mth.eq.501) go to 320
   if (mth.eq.2.and.hp.ne.hnt) go to 350
   if (iref.ne.0.and.ig.gt.igr) go to 310
   lout=ivdat+n1i*ning-ig
   lin=lz+nl*(nz+iz-1)+1
   a(lout)=b(lin)
   if (hp.eq.hnt) go to 310
   if (hp.eq.hga) go to 310
   if (mth.eq.500) go to 310
   if (mth.ge.201.and.mth.le.207) go to 310
   if (mth.eq.4) k004=1
   if (k004.eq.1.and.mth.ge.50.and.mth.le.91) go to 310
   if (mth.eq.103) k103=1
   if (k103.eq.1.and.mth.ge.600.and.mth.le.649) go to 310
   if (mth.eq.104) k104=1
   if (k104.eq.1.and.mth.ge.650.and.mth.le.699) go to 310
   if (mth.eq.105) k105=1
   if (k105.eq.1.and.mth.ge.700.and.mth.le.749) go to 310
   if (mth.eq.106) k106=1
   if (k106.eq.1.and.mth.ge.750.and.mth.le.799) go to 310
   if (mth.eq.107) k107=1
   if (k107.eq.1.and.mth.ge.800.and.mth.le.849) go to 310
   if (mth.eq.16) k016=1
   if (k016.eq.1.and.mth.ge.875.and.mth.le.891) go to 310
   lout=ivdat+2*ning-ig
   a(lout)=a(lout)+b(lin)
  310 continue
   if (iref.eq.0) go to 115
   if (igr.gt.ig) go to 115
   lout=ivdat+n1i*ning-igr
   lin=idref-1+lz+nl*nz+1
   a(lout)=a(lout)-b(lin)
   if (abs(a(lout)).lt.stol*b(lin)) a(lout)=0
   go to 115
   ! flux and total cross section
  320 continue
   do il=1,nl
      ! flux
      lout=ivdat+il*ning-ig
      lin=lz+nl*(iz-1)+il
      a(lout)=b(lin)
      ! total cross sections
      if (iref.eq.0.or.ig.le.igr) then
         lout=ivdat+(nl+il)*ning-ig
         lin=lz+nz*nl+nl*(iz-1)+il
         a(lout)=b(lin)
      endif
      if (iref.ne.0.and.igr.le.ig) then
         lout=ivdat+(nl+il)*ning-igr
         lin=idref-1+lz+nz*nl+il
         a(lout)=a(lout)-b(lin)
         if (abs(a(lout)).ge.stol*b(lin)) igw=1
         if (abs(a(lout)).lt.stol*b(lin).and.igw.gt.0) a(lout)=0
      endif
   enddo
   go to 115
   ! charged-particle elastic and flux
  350 continue
   lout=ivdat+ning-ig
   lin=lz+nl*(iz-1)+1
   a(lout)=b(lin)
   lout=lout+ning
   lin=lin+nl*nz
   a(lout)=b(lin)
   lout=lout+ning
   a(lout)=b(lin)
   go to 115
  400 continue

   !--add up separate contributions to production reactions
  410 continue
   do 430 i=1,n1d
   n1i=n1d-i+1
   do 420 j=1,n1i-1
   n1j=j
   if (hvps(n1j).eq.hvps(n1i)) go to 440
  420 continue
  430 continue
   go to 480
  440 continue
   do k=1,ning
      lout=ivdat+n1j*ning-k
      lin=ivdat+n1i*ning-k
      a(lout)=a(lout)+a(lin)
   enddo
   if (n1i.ne.n1d) then
      nmove=(n1d-n1i)*ning
      do k=1,nmove
         lout=ivdat+(n1i-1)*ning+k-1
         lin=ivdat+n1i*ning+k-1
         a(lout)=a(lin)
      enddo
      nmove=n1d-n1i
      do k=1,nmove
         hvps(n1i+k-1)=hvps(n1i+k)
      enddo
   endif
   n1d=n1d-1
   go to 410
  480 continue

   !--output the vector data to the scratch files
   ndex=ivdat
   mdex=ivdat
   do i=1,n1d

      !--calculate nfg and nlg
      nfg=ning+1
      nlg=0
      do k=1,ning
         nxsd=ndex-1+k
         if (abs(a(nxsd)).ge.small) then
            if (nfg.gt.k) nfg=k
            if (nlg.lt.k) nlg=k
         endif
      enddo

      !--if entire partial is zero, keep group 1 only
      if (nlg.eq.0) then
         nfg=1
         nlg=1
      endif

      !--crunch data using nfg and nlg
      do k=1,ning
         if (nfg.le.k.and.nlg.ge.k) then
            a(mdex)=a(ndex)
            mdex=mdex+1
         endif
         ndex=ndex+1
      enddo

      !--write vector control data
      n=(ivcon-1)/mult+i
      read(hvps(i),'(a8)') ta(n)
      ia(ivcon+n1d*mult-1+i)=nfg
      ia(ivcon+n1d*mult+n1d-1+i)=nlg
   enddo

   !--write vector blocks to nscrt3
   nwds=0
   ivnow=ivdat
   do i=1,n1d
      loc=ivcon-1+mult*n1d+i
      nw=ia(loc+n1d)-ia(loc)+1
      if (nwds+nw.ge.maxw) then
         irec3=irec3+1
         write(nscrt3) (a(ivnow+j-1),j=1,nwds)
         ivnow=ivnow+nwds
         nritev=nritev+1
         nwds=0
      endif
      nwds=nwds+nw
   enddo
   if (nwds.gt.0) then
      irec3=irec3+1
      write(nscrt3) (a(ivnow+j-1),j=1,nwds)
      ivnow=ivnow+nwds
      nritev=nritev+1
   endif

   !--write vector control record
   irec2=irec2+1
   nw=2*n1d
   ihvcon=(ivcon-1)/mult
   write(nscrt2)&
     (ha(ihvcon+i),i=1,n1d),(ia(ivcon+n1d*mult+i-1),i=1,nw)
   nritev=nritev+1

   !--vector production complete
   next=ivcon
   return
   end subroutine vector

   subroutine matrix(irec3,irec2,mfd,hp,n2d,nscr,nscr2,iref)
   !-------------------------------------------------------------------
   ! Produces matrix data and writes it to scratch tape.
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use util ! provides openz,repoz,closz
   ! externals
   integer::irec3,irec2,mfd,n2d,nscr,nscr2,iref
   character::hp*(*)
   ! internals
   integer::nscrr,irefr,nb,nw,i,j,irite,iskip,lr
   integer::nlord,njcon,lordn,ioff,ik,loca
   integer::nwmcr,nwm,nwdc,nwds,nhmtx,nz,ihmcon
   real(kr)::b(350)

   !--open scratch files for reaction data
   nscrr=nscrt4
   if (nscr.lt.0) nscrr=-nscrt4
   irefr=nscrt7
   if (iref.lt.0) irefr=-nscrt7
   if (iref.eq.0) irefr=0
   call openz(nscrr,1)
   call openz(irefr,1)

   !--set pointers for matrix control and matrix data blocks
   nwmcr=mult+2+2*noutg
   nwm=nwmcr
   imcon=next
   next=imcon+nwm
   nwdc=noutg+ning
   icdat=next
   next=icdat+nwdc
   nwds=2*noutg
   ijgll=next
   next=ijgll+nwds
   imdat=next
   nhmtx=(imcon-1)/mult+1

   !--loop over matrix partials
   call repoz(nscr)
   call repoz(iref)
   call contio(nscr,0,0,b(1),nb,nw)
   if (iref.ne.0) call contio(iref,0,0,b(1),nb,nw)
   nritem=0
   do i=1,n2d

      !--process this reaction
      do j=1,nwm
         a(j-1+imcon)=0
      enddo
      do j=1,nwdc
         a(icdat-1+j)=0
      enddo
      lord1=0
      ng2z=0
      irite=0
      jconst=0
      iskip=0
      do while (iskip.eq.0)
         if (iref.ne.0) call contio(iref,0,0,b(1),nb,nw)
         call contio(nscr,0,0,b(1),nb,nw)
         nz=l2h
         lr=n1h
         if (nz.ge.iz) then
            if (mfh.eq.mfd) iskip=1
            if (mfh.eq.17.and.mfd.eq.16) iskip=1
         endif
         if (iskip.eq.0) then
           if (iref.ne.0) call tosend(iref,0,0,b(1))
           call tosend(nscr,0,0,b(1))
         endif
      enddo
      call findf(math,mfh,mth,nscr)
      call repoz(nscrr)
      call repoz(irefr)
      call contio(nscr,nscrr,0,b(1),nb,nw)
      if (iref.ne.0) then
         call findf(math,mfh,mth,iref)
         call contio(iref,irefr,0,b(1),nb,nw)
      endif
      do j=imdat,isiza
         a(j)=0
      enddo

      !--find hollerith names
      call hname(hmtx(i),hp,mth,lr,0)

      !--read and copy this reaction
      !--while computing banding information
      call band(nscr,iref,ia(ijgll),ia(ijgll+noutg),nscrr,irefr)

      !--read reaction again from scratch file
      !--and shuffle the numbers into the bands
      call shufl(nscrr,irefr,ia(ijgll),ia(ijgll+noutg),&
        nscrt3,irec3,irite)

      !--calculate lord(n),  hmtx(n)
      read(hmtx(i),'(a8)') ta(nhmtx)
      nlord=imcon+mult
      njcon=nlord+1
      if (irite.ne.0) then
         ia(nlord)=lord1
         ia(njcon)=jconst
         call tosend(nscr,0,0,b(1))
         irec2=irec2+1
         ihmcon=(imcon-1)/mult+1
         nw=nwmcr-mult
         write(nscr2) ha(ihmcon),(ia(imcon+mult+j-1),j=1,nw)
         irite=irite+1
      else
         call tosend(nscr,0,0,b(1))
         if (lord1.eq.0) lord1=maxord+1
         ia(nlord)=lord1
         ia(njcon)=jconst
         ia(njcon+1)=1
         ia(njcon+1+noutg)=1
         ihmcon=(imcon-1)/mult+1
         nw=nwmcr-mult
         write(nscr2) ha(ihmcon),(ia(imcon+mult+j-1),j=1,nw)
         irite=irite+1
         lordn=maxord+1
         do j=1,lordn
            a(j-1+imdat)=0.
         enddo
         irec3=irec3+1
         write(nscrt3) (a(imdat+j-1),j=1,lordn)
         irite=irite+1
      endif

      !--write constant sub-block if present
      if (jconst.ne.0) then
         ioff=ning-jconst
         do ik=1,jconst
            loca=icdat-1+noutg+ik
            a(loca)=a(loca+ioff)
         enddo
         nwds=noutg+jconst
         irec3=irec3+1
         write(nscrt3) (a(icdat+j-1),j=1,nwds)
         irite=irite+1
      endif

      !--close loop over matrix partials
      nritem=nritem+irite
   enddo

   !--finished
   next=imcon
   call closz(nscrr)
   call closz(irefr)
   return
   end subroutine matrix

   subroutine band(nscr,iref,jg1lo,jg1hi,nscrr,irefr)
   !-------------------------------------------------------------------
   ! Reads one entire mt file, copies it to scratch files,
   ! and computes the banding parameters.
   !-------------------------------------------------------------------
   use endf ! provide endf routines and variables
   use util ! provides error
   ! externals
   integer::nscr,iref,nscrr,irefr
   integer(k4)::jg1lo(*),jg1hi(*)
   ! internals
   integer::i,irinp,ig,igr,ig2lo,idref,ig2lor,jg1,ik
   integer::nb,nw,ng2,loc,ng2r,jg2,n,nl
   integer,parameter::maxb=30000
   real(kr)::b(maxb)

   !--initialize banding data
   do i=1,noutg
      jg1lo(i)=ning+1
      jg1hi(i)=0
   enddo

   !--loop over records on input tape
   irinp=1
   igr=0
   nl=l1h
   maxord=nl-1

   !--main group loop
   ig=0
   do while (ig.lt.ning)
      call listio(nscr,nscrr,0,b(irinp),nb,nw)
      ng2=l1h
      ig2lo=l2h
      ig=n2h
      loc=irinp+nw
      do while (nb.ne.0)
         if (loc+302.gt.maxb) call error('band',&
           'input too large',' ')
         call moreio(nscr,nscrr,0,b(loc),nb,nw)
         loc=loc+nw
      enddo

      !--reference group loop
     120 continue
         if (iref.ne.0.and.ig.gt.igr) then
            idref=loc
            if (idref+302.gt.maxb) call error('band',&
              'input too large',' ')
            call listio(iref,irefr,0,b(idref),nb,nw)
            ng2r=l1h
            ig2lor=l2h
            igr=n2h
            loc=idref+nw
            do while (nb.ne.0)
               if (loc+302.gt.maxb) call error('band',&
                 'input too large',' ')
               call moreio(iref,irefr,0,b(loc),nb,nw)
               loc=loc+nw
            enddo
         endif

         !--analyze banding for matrix sub-blocks
         if (ig.ne.0.and.ig2lo.ne.0) then
            if (iref.eq.0.or.ig.le.igr) then
              jg1=ning+1-ig
              do ik=2,ng2
                 jg2=noutg+1-(ig2lo+ik-2)
                 if (jg1.lt.jg1lo(jg2)) jg1lo(jg2)=jg1
                 if (jg1.gt.jg1hi(jg2)) jg1hi(jg2)=jg1
               enddo
            endif
            if (iref.ne.0.and.igr.le.ig) then
               jg1=ning+1-igr
               do ik=2,ng2r
                  jg2=noutg+1-(ig2lor+ik-2)
                  if (jg1.lt.jg1lo(jg2)) jg1lo(jg2)=jg1
                  if (jg1.gt.jg1hi(jg2)) jg1hi(jg2)=jg1
               enddo
            endif
         endif
      if (iref.ne.0.and.igr.lt.ig) go to 120
   enddo

   !--load jband and ijj for the matrix control record
   loc=imcon+mult+1
   do i=1,noutg
      n=jg1hi(i)-jg1lo(i)+1
      if (n.lt.0) n=0
      ia(loc+i)=n
      ia(loc+noutg+i)=jg1hi(i)
   enddo
   return
   end subroutine band

   subroutine shufl(nscr,iref,jg1lo,jg1hi,nout3,irec3,irite)
   !-------------------------------------------------------------------
   ! Reads one entire mt file and rearranges data
   ! using the banding information computed by subroutine band.
   ! Make multiple passes through the files if necessary.
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use util ! provides repoz,error
   ! externals
   integer::nscr,iref,nout3,irec3,irite
   integer(k4)::jg1lo(*),jg1hi(*)
   ! internals
   integer::i0,num,nwds,nsb,i,idone,n,imax,irinp,nz
   integer::ig,igr,ig2lo,idref,ig2lor,ik,jg2,jg1,il
   integer::i1,imax1,nb,nw,lz,nl,ng2,loc,ng2r,niloc
   integer::kp,jp,np,noloc,j
   integer,parameter::maxb=30000
   real(kr)::b(maxb)

   !--loop over records on input tape
   ng2z=0
  100 continue
   i0=ng2z+1

   !--how many groups and sub blocks will fit?
   num=0
   nwds=0
   nsb=0
   i=i0-1
   idone=0
   do while (i.lt.noutg.and.idone.eq.0)
      i=i+1
      n=jg1hi(i)-jg1lo(i)+1
      if (n.lt.0) n=0
      n=n*(maxord+1)
      nwds=nwds+n
      if (imdat+nwds.gt.isiza) then
         idone=1
      else
         if (num+n.ge.maxw) then
            imax=i-1
            nsb=nsb+1
            num=0
         endif
         num=num+n
      endif
   enddo
   if (idone.eq.0) imax=noutg

   !--read the head record
   irinp=1
   call repoz(nscr)
   call repoz(iref)
   call contio(nscr,0,0,b(irinp),nb,nw)
   lz=6
   igr=0
   nl=l1h
   nz=l2h
   if (iref.ne.0) call contio(iref,0,0,b(irinp),nb,nw)

   !--read the group records
  110 call listio(nscr,0,0,b(irinp),nb,nw)
   ng2=l1h
   ig2lo=l2h
   ig=n2h
   loc=irinp+nw
   do while (nb.ne.0)
      if (loc+302.ge.maxb) call error('shufl','input too large',' ')
      call moreio(nscr,0,0,b(loc),nb,nw)
      loc=loc+nw
   enddo
  120 if (iref.ne.0.and.ig.gt.igr) then
      idref=loc
      if (idref+302.ge.maxb) call error('shufl',&
        'input too large',' ')
      call listio(iref,0,0,b(idref),nb,nw)
      ng2r=l1h
      ig2lor=l2h
      igr=n2h
      loc=idref+nw
      do while (nb.ne.0)
         if (loc+302.ge.maxb) call error('shufl',&
           'input too large',' ')
         call moreio(iref,0,0,b(loc),nb,nw)
         loc=loc+nw
      enddo
   endif

   !--store constant sub-block if present
   if (ig.eq.0) go to 135
   if (ig2lo.eq.0) go to 140
   go to 145
  135 do ik=1,ng2
      jg2=noutg-ig2lo-ik+2
      a(icdat-1+jg2)=b(irinp+lz+ik-1)
   enddo
   go to 240
  140 jconst=ig
   jg1=ning+1-ig
   a(icdat-1+noutg+jg1)=b(irinp+lz+1)
   go to 240

   !--store matrix sub-blocks
  145 if (iref.eq.0.or.ig.le.igr) then
      jg1=ning+1-ig
      do ik=2,ng2
         niloc=irinp+lz+nl*nz*(ik-1)+nl*(iz-1)
         jg2=noutg+1-(ig2lo+ik-2)
         if (jg2.ge.i0.and.jg2.le.imax) then
            kp=jg1hi(jg2)
            jp=0
            i=i0-1
            idone=0
            do while (i.lt.noutg.and.idone.eq.0)
               i=i+1
               np=jg1hi(i)-jg1lo(i)+1
               if (np.lt.0) np=0
               if (i.eq.jg2) then
                  idone=1
               else
                  jp=jp+np
               endif
            enddo
            noloc=imdat+nl*jp+kp-jg1
            do il=1,nl
               a(noloc+np*(il-1))=b(niloc+il-1)
            enddo
         endif
      enddo
   endif
   if (iref.ne.0.and.igr.le.ig) then
      jg1=ning+1-igr
      do ik=2,ng2r
         niloc=idref+lz+nl*nz*(ik-1)
         jg2=noutg+1-(ig2lor+ik-2)
         if (jg2.ge.i0.and.jg2.le.imax) then
            kp=jg1hi(jg2)
            jp=0
            i=i0-1
            idone=0
            do while (i.lt.noutg.and.idone.eq.0)
               i=i+1
               np=jg1hi(i)-jg1lo(i)+1
               if (np.lt.0) np=0
               if (i.eq.jg2) then
                  idone=1
               else
                  jp=jp+np
               endif
            enddo
            noloc=imdat+nl*jp+kp-jg1
            do il=1,nl
               a(noloc+np*(il-1))=a(noloc+np*(il-1))-b(niloc+il-1)
            enddo
         endif
      enddo
   endif
  240 continue
   if (iref.ne.0.and.igr.lt.ig) go to 120
   if (ig.lt.ning) go to 110
   if (nl+1.gt.lord1) lord1=nl+1
   if (lord1.gt.maxord+1) lord1=maxord+1

   !--write out this range of sub-block data
   nwds=0
   i1=imdat
   imax1=imax+1
   do i=i0,imax1
      if (i.le.imax) then
         np=jg1hi(i)-jg1lo(i)+1
         if (np.lt.0) np=0
         np=np*nl
         if (nwds+np.lt.maxw) then
            nwds=nwds+np
         else
            if (nwds.eq.0) then
               nwds=nwds+np
            else
               irec3=irec3+1
               write(nout3) (a(i1+j-1),j=1,nwds)
               irite=irite+1
               i1=i1+nwds
               nwds=0
               nwds=nwds+np
            endif
         endif
      else
         if (nwds.eq.0) then
            nwds=nwds+np
         else
            irec3=irec3+1
            write(nout3) (a(i1+j-1),j=1,nwds)
            irite=irite+1
            i1=i1+nwds
            nwds=0
            nwds=nwds+np
         endif
      endif
   enddo
   ng2z=imax
   if (ng2z.lt.noutg) go to 100
   return
   end subroutine shufl

   subroutine lst1io(nin,nout,nscr,a,nb,nw,namax)
   !-------------------------------------------------------------------
   ! Convert a gendf file 1 mt 451 list record to opposite format.
   ! Routine lst1io processes the entire record.
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use util ! provides error
   ! externals
   integer::nin,nout,nscr,nb,nw,namax
   real(kr)::a(namax)
   ! internals
   integer::nb1,nw1,ng,nwh,nx,nwd,i,lim,l,nt,ns,iwd,mtw
   integer::j
   real(kr)::z(100)
   character(4)::tz(100)
   equivalence(z(1),tz(1))
   character(11)::field(6)
   character(40)::hfmt
   character(4)::blank='    '
   nwd=0

   !--input
   if (nin.eq.0) go to 200
   if (nin.lt.0) go to 140
   call contio(nin,0,0,a,nb1,nw1)
   nw=6+n1h
   ng=l1h
   if (ng.eq.0) ng=l2h
   if (ntw.eq.0) go to 130
   if (ntw.ge.14) go to 120
   nwh=(ntw*4+10)/11
   nx=nwh*11-ntw*4
   nwd=6-nwh
   if (nwd.le.0) go to 120
   if (nx.ne.0) write(hfmt,&
     &'(''('',i2,''a4,'',i2,''x,1p'',i1,''e11.0)'')')&
     ntw,nx,nwd
   if (nx.eq.0) write(hfmt,&
     '(''('',i2,''a4,'',''1p'',i1,''e11.0)'')')&
     ntw,nwd
   read(nin,hfmt) (tz(i),i=1,ntw),(a(i+nw1+ntw),i=1,nwd)
   do i=1,ntw
      a(i+nw1)=z(i)
   enddo
   nw1=nw1+ntw+nwd
   go to 130
  120 continue
   read(nin,'(16a4,a2)') (tz(i),i=1,ntw)
   do i=1,ntw
      a(i+nw1)=z(i)
   enddo
   nw1=nw1+ntw
  130 continue
   lim=nw-nw1
   read(nin,'(6e11.0)') (a(i+nw1),i=1,lim)
   go to 200
  140 continue
   call listio(nin,0,0,a,nb,nw)
   ng=l1h
   if (ng.eq.0) ng=l2h
   l=1
  150 continue
   if (nb.eq.0) go to 200
   l=l+nw
   call moreio(nin,0,0,a(l),nb,nw)
   if ((l-1+nw).gt.namax)&
     call error('lst1io','storage exceeded',' ')
   go to 150

   !--output
  200 continue
   c1h=a(1)
   c2h=a(2)
   l1h=nint(a(3))
   l2h=nint(a(4))
   n1h=nint(a(5))
   n2h=nint(a(6))
   nw=6+n1h
   nb=0
   if (nout.eq.0.and.nscr.eq.0) go to 400
   if (nout.le.0.and.nscr.le.0) go to 300
   nt=nout
   if (nout.le.0) nt=0
   ns=nscr
   if (nscr.le.0) ns=0
   call contio(0,nt,ns,a,nb1,nw1)
   nwh=(ntw*4+10)/11
   nx=nwh*11-ntw*4
   nwd=6-nwh
   if (ntw.ge.14) go to 240
   if (nwd.le.0) go to 240
   do iwd=1,nwd
      if (iwd+nw1+ntw.gt.nw) a(nwd+nw1+ntw)=0
      call a11(a(iwd+nw1+ntw),field(iwd))
   enddo
   if (nx.ne.0) write(hfmt,&
     '(''('',i1,''a4,'',i2,''x,'',i1,&
     &''(a11),i4,i2,i3,i5)'')')&
     ntw,nx,nwd
   if (nx.eq.0) write(hfmt,&
     &'(''('',i1,''a4,'',i1,''(a11),i4,i2,i3,i5)'')')&
     nwd,ntw
   if (nout.gt.0) then
      nsh=nsh+1
      do i=1,ntw
         z(i)=a(i+nw1)
      enddo
      write(nt,hfmt) (tz(i),i=1,ntw),(field(i),i=1,nwd),&
        math,mfh,mth,nsh
   endif
   if (nscr.le.0) go to 260
   nsc=nsc+1
   do i=1,ntw
      z(i)=a(i+nw1)
   enddo
   write(ns,hfmt) (tz(i),i=1,ntw),(field(i),i=1,nwd),&
     math,mfh,mth,nsc
  240 continue
   mtw=17-ntw
   if (nout.gt.0) then
      nsh=nsh+1
      do i=1,ntw
         z(i)=a(i+nw1)
      enddo
      if (ntw.eq.17) write(nt,'(16a4,a2,i4,i2,i3,i5)')&
        (tz(i),i=1,ntw),math,mfh,mth,nsh
      if (ntw.lt.17) write(nt,'(16a4,a2,i4,i2,i3,i5)')&
        (tz(i),i=1,ntw),('    ',i=1,mtw),math,mfh,mth,nsh
   endif
   if (nscr.gt.0) then
      nsc=nsc+1
      do i=1,ntw
         z(i)=a(i+nw1)
      enddo
      if (ntw.eq.17) write(ns,'(16a4,a2,i4,i2,i3,i5)')&
        (tz(i),i=1,ntw),math,mfh,mth,nsc
      if (ntw.lt.17) write(ns,'(16a4,a2,i4,i2,i3,i5)')&
        (tz(i),i=1,ntw),(blank,i=1,mtw),math,mfh,mth,nsc
   endif
  260 continue
   if (nwd.lt.0) nwd=0
   nw1=nw1+nwd+ntw
   lim=nw-nw1
   do i=1,lim,6
      do j=1,6
         if ((j+nw1).gt.nw) a(j+nw1)=0
         call a11(a(j+nw1),field(j))
      enddo
      if (nout.gt.0) then
         nsh=nsh+1
         write(nt,'(6(a11),i4,i2,i3,i5)')&
           (field(j),j=1,6),math,mfh,mth,nsh
      endif
      if (nscr.gt.0) then
         nsc=nsc+1
         write(ns,'(6(a11),i4,i2,i3,i5)')&
           (field(j),j=1,6),math,mfh,mth,nsc
      endif
      nw1=nw1+6
   enddo
   if (nout.ge.0.and.nscr.ge.0) go to 400
   nt=nout
   if (nout.ge.0) nt=0
   ns=nscr
   if (nscr.ge.0) ns=0
   go to 310
  300 continue
   nt=nout
   ns=nscr
  310 continue
   call listio(0,nt,ns,a,nb,nw)
   l=1
   do while (nb.ne.0)
      l=l+nw
      call moreio(0,nt,ns,a(l),nb,nw)
      if ((l-1+nw).gt.namax) call error('lst1io',&
        'storage exceeded',' ')
   enddo

   !--finished
  400 continue
   return
   end subroutine lst1io

end module matxsm


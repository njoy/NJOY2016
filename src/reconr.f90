module reconm
   ! Module to provide reconr for NJOY2016
   use locale
   implicit none
   private
   public reconr

   ! global variables for reconr

   integer::nin,nout
   real(kr)::zain,awin
   real(kr)::zai,el,eh,err,errmax,errint,za,awr,tempr,q18
   real(kr)::tempi
   real(kr)::elis,sta,efmax
   integer::lis,lis0,nfor,lrel,nver
   integer::lfw,mata,itype,lrp,lfi,lssf,lrx
   integer,parameter::nmtmax=10
   integer::mtr4,mtr18,mtr(nmtmax),mtrt(nmtmax),nmtr
   integer::mt103,mt104,mt105,mt106,mt107
   integer::mpmin,mpmax,mdmin,mdmax,mtmin,mtmax,m3min,m3max,m4min,m4max
   integer::nxc,ngo,mtr522,ncards
   integer::isot(20),modet(20),ibaset(20),isect,nsect
   real(kr)::abnt(20),elt(20),eht(20)
   real(kr)::eresl,eresr,eresh,eresu,eresm
   integer::nodes,nunr,maxres
   real(kr)::elst(20),slst(20,4),enxt(20),snxt(20,4)
   integer::ilast(20)
   integer::mmtres(10),mcards(10)
   integer::nmtres,ncoef,nsig,nresp
   real(kr)::spin,ascat,thr6
   character(4),dimension(:),allocatable::card
   real(kr),dimension(:),allocatable::enode
   real(kr),dimension(:),allocatable::eunr
   real(kr),dimension(:),allocatable::res
   real(kr),dimension(:,:),allocatable::tr,ti
   real(kr),dimension(:),allocatable::dict
   integer,dimension(:),allocatable::mfs,mts,ncs
   real(kr),dimension(:),allocatable::sunr

   integer,parameter::maxunr=500
   integer,parameter::nodmax=800000
   integer,parameter::nbufg=2000
   integer,parameter::nbufr=2000
   integer,parameter::nbuf=2000
   integer,parameter::nbufl=2000

   !--option for sammy method
   !integer::isammy=0  ! no sammy calculation
   integer::isammy=1  ! sammy calculation

contains

   subroutine reconr
   !-------------------------------------------------------------------
   !
   ! Reconstruct pointwise cross sections
   !
   ! Resonance cross sections are calculated using methods from the
   ! earlier versions of NJOY (which are descended from those of
   ! the early RESEND code) or the methods developed for the SAMMY
   ! R-matrix fitting code, depending on the format used.
   !
   ! This program generates an energy grid which is the union of an
   ! input grid (if any), the resonance energies (if any), and the
   ! energies of cross sections in mf3 and mf13 (or mf23). The
   ! pointwise cross sections are then computed on this grid and
   ! points are added so that the resonance cross sections and any
   ! cross sections represented by non-linear interpolation are
   ! reproduced within a specified tolerance by linear interpolation.
   ! Psi-Chi reconstruction can be used if desired.  Sections which
   ! are not cross sections (mubar,nubar) and photon multiplicities
   ! (MF12) are not processed.  Redundant reactions are reconstructed
   ! to be the sum of their parts.  The pendf tape contains point
   ! cross sections in MF3 and MF13 (or MF23) and a description of
   ! the processing in MF1.  The MF1 dictionary is updated.  The
   ! C1 and C2 fields of the second card in MF1 contain the temper-
   ! ature and reconstruction tolerance, respectively.  An MF2
   ! appropriate to no resonance parameters is constructed with
   ! the potential scattering length added.
   !
   ! If unresolved parameters are present, the infinitely dilute
   ! cross sections are computed on a special energy grid chosen to
   ! preserve the required interpolation properties.  This table is
   ! added to the pendf tape using a special format in MF2/MT152,
   ! and the table is also used to compute the unresolved contribu-
   ! tions in MF3.  This allows resolved resonance cross sections
   ! which overlap the resolved range to be recovered by subtraction.
   !
   !---input specifications (free format)---------------------------
   !
   ! card 1
   !    nendf    unit for endf tape
   !    npend    unit for pendf tape
   ! card 2
   !    tlabel   66 character label for new pendf tape
   !             delimited with quotes, ended with /.
   ! card 3
   !    mat      material to be reconstructed
   !    ncards   number of cards of descriptive data for new mf1
   !             (default=0)
   !    ngrid    number of user energy grid points to be added.
   !             (default=0)
   ! card 4
   !    err      fractional reconstruction tolerance used when
   !             resonance-integral error criterion (see errint)
   !             is not satisfied.
   !    tempr    reconstruction temperature (deg kelvin)
   !             (default=0)
   !    errmax   fractional reconstruction tolerance used when
   !             resonance-integral error criterion is satisfied
   !             (errmax.ge.err, default=10*err)
   !    errint   maximum resonance-integral error (in barns)
   !             per grid point (default=err/20000)
   !             (note: the max cross section difference for
   !             linearization, errlim, and for reconstruction,
   !             errmin, are also tied to errint.  to get maximum
   !             accuracy, set errint to a very small number.
   !             for economical production, use the defaults.)
   ! card 5
   !    cards    ncards of descriptive comments for mt451
   !             each card delimited with quotes, ended with /.
   ! card 6
   !    enode    users energy grid points
   !
   !     cards 3, 4, 5, 6 must be input for each material desired
   !     mat=0/ terminates execution of reconr.
   !
   !-------------------------------------------------------------------
   use util   ! provides timer,repoz,closz,error
   use endf   ! provides endf routines and variables
   use samm   ! provides s2sammy,desammy
   use mainio ! provides nsysi,nsyso,nsyse
   ! internals
   integer::i,nb,nw,nrtot,nwscr,n6,nx,nsub,iold,inew,ngrid
   integer::nendf,npend,nscr1,nscr2,nscr3,nscr4,intunr
   real(kr)::time
   real(kr)::rlabel(17),z(17)
   character(4)::tlabel(17),tz(17)
   equivalence(tlabel(1),rlabel(1))
   equivalence(z(1),tz(1))
   character(66)::text
   real(kr),dimension(:),allocatable::scr
   real(kr),parameter::elow=1.e-5_kr
   real(kr),parameter::ehigh=20.e6_kr
   real(kr),parameter::elarge=9.e9_kr
   logical there

   !--set samrml options
   Want_Partial_Derivs=.false.
   Want_Angular_Dist=.false.
   Want_SAMRML_BW=.false.
   Want_SAMRML_RM=.false.

   !--define local i/o units
   nscr1=-10
   nscr2=11
   nscr3=12

   !--write heading and read user input
   call timer(time)
   write(nsyso,'(/'' reconr...'',&
     &''reconstruct pointwise cross sections in pendf format'',&
     &7x,f8.1,''s'')') time
   write(nsyse,'(/'' reconr...'',59x,f8.1,''s'')') time
   read(nsysi,*) nendf,npend
   call openz(nendf,0)
   call openz(npend,1)
   read(nsysi,*) text
   read(text,'(16a4,a2)') (tlabel(i),i=1,17)
   write(nsyso,'(/&
     &'' unit for endf tape ................... '',i10/&
     &'' unit for pendf tape .................. '',i10)')&
     nendf,npend
   write(nsyso,'(/&
     &'' label for pendf tape''/1x,66(''-'')/1x,17a4)')&
     (tlabel(i),i=1,17)

   !--read tape identification record from input tape.
   !--write new tape identification label on output tape.
   nin=nendf
   nout=npend
   call repoz(nin)
   call repoz(nout)
   call tpidio(nin,0,0,z,nb,nw)
   write(nsyso,'(/'' tape label''/1x,66(''-'')/1x,17a4)')&
     (tz(i),i=1,nw)
   math=1
   mfh=0
   mth=0
   nsh=0
   nw=17
   call tpidio(0,nout,0,rlabel,nb,nw)
   if (npend.gt.0) nscr1=10
   mata=1

   !--loop through users materials
   do while (mata.ne.0)
      allocate(enode(nodmax))

      !--read users input for next material.
      call ruina
      if (mata.ne.0) then

         !--open files for this material
         call openz(-nscr2,1)
         call openz(-nscr3,1)
         call openz(nscr1,1)
         call repoz(nscr1)
         call repoz(nscr2)
         call repoz(nscr3)
         nscr4=13
         nwscr=npage+50
         allocate(scr(nwscr))

         !--position to required material
         call findf(mata,1,451,nin)

         !--process file 1
         tempi=0
         call contio(nin,0,0,scr,nb,nw)
         za=c1h
         awr=c2h
         lrp=l1h
         lfi=l2h
         n6=n2h
         call contio(nin,0,0,scr,nb,nw)
         if (n1h.ne.0) then
            iverf=4
            nx=n6
         else if (n2h.eq.0) then
            iverf=5
         else
            iverf=6
         endif
         call skiprz(nin,-1)
         if (iverf.ge.5) then
            call contio(nin,0,0,scr,nb,nw)
            elis=c1h
            sta=c2h
            lis=l1h
            lis0=l2h
            nfor=n2h
         endif
         nsub=10
         zain=1
         awin=1
         if (iverf.ge.6) then
            call contio(nin,0,0,scr,nb,nw)
            nsub=n1h
            zain=int(nsub/10)
            awin=c1h
            efmax=c2h
            lrel=l1h
            nver=n2h
         endif
         if (10*nint(zain).ne.nsub.and.nsub.ne.3) then
            call error('reconr','illegal nsub for reconr',' ')
         endif
         call hdatio(nin,0,0,scr,nb,nw)
         if (iverf.ge.6) tempi=scr(1)
         do i=1,17
            z(i)=scr(6+i)
         enddo
         write(nsyso,'(/&
           &'' processing mat '',i4,&
           &'' in endf-'',i1,'' format''/&
           &1x,66(''-'')/1x,17a4/)')&
           mata,iverf,(tz(i),i=1,17)
         do while (nb.ne.0)
            call moreio(nin,0,0,scr,nb,nw)
         enddo
         if (iverf.ne.4) nx=n2h
         call anlyzd(nin,nx)
         call tofend(nin,0,0,scr)

         !--scan file 2 to check whether sammy method is needed.
         !--get sizes for sammy arrays and allocate some arrays.
         lrx=0
         nmtres=0
         if (lrp.eq.1) then
            if (isammy.gt.0) then
              allocate(res(maxres))
              call s2sammy(nin,res,maxres,mmtres,nmtres)
              if (nmtres.gt.0) then
                 do i=1,nmtres
                    if (mmtres(i).eq.103) mmtres(i)=600
                    if (mmtres(i).eq.104) mmtres(i)=650
                    if (mmtres(i).eq.105) mmtres(i)=700
                    if (mmtres(i).eq.106) mmtres(i)=750
                    if (mmtres(i).eq.107) mmtres(i)=800
                 enddo
              endif
              call repoz(nin)
              call findf(mata,2,0,nin)
              deallocate(res)
            endif
         endif

         !--process file 2
         allocate(eunr(maxunr))
         nunr=0
         if (lrp.ne.1) then
            write(nsyso,'(/'' mat has no resonance parameters'')')
            nrtot=0
            eresl=elow
            eresr=ehigh
            eresh=ehigh
            eresu=elarge
            eresm=elarge
            if (lrp.eq.0 .or. lrp.gt.1) then
               call contio(nin,0,0,scr,nb,nw)
               call contio(nin,0,0,scr,nb,nw)
               call contio(nin,0,0,scr,nb,nw)
               call contio(nin,0,0,scr,nb,nw)
               spin=c1h
               ascat=c2h
            elseif (lrp.lt.0) then
               spin=0
               ascat=0
            endif
         else
            nrtot=maxres
            allocate(res(nrtot))
            ncoef=1
            do i=1,10
               mcards(i)=0
            enddo
            call rdfil2(nrtot,intunr)
            if (nmtres.eq.0) then
               nsig=5
            else
               nsig=4
               write(nsyso,'(/''   samm resonance reactions:'',10i5)')&
                 (mmtres(i),i=1,nmtres)
               if (nmtres.gt.2) then
                  do i=3,nmtres
                     if (mmtres(i).ne.18) nsig=nsig+1
                  enddo
               endif
            endif
            write(nsyso,'(''   samm max legendre order:'',i3)') ncoef-1
            if (ncoef.gt.1)&
              write(nsyso,&
               '(''   generating File 4 for resonance angular distributions'')')
            if (lrp.eq.1) write(nsyso,'(/&
              &'' mat has no unresolved resonance parameters'')')
            ! temporary meaning of lrp------
            !   1 means resolved only
            !   3 means resolved and unresolved
            if (lrp.eq.3) call genunr(nin,intunr)
            if (lssf.gt.0) eresh=eresr
         endif
         deallocate(eunr)

         !--process file 3 and file 13
         call lunion(nin,nscr1,nscr2,ngrid)
         deallocate(enode)

         !--compute resonance cross sections
         if (lrp.eq.1.or.lrp.eq.3) then
            call resxs(nscr2,nscr3,nscr4,nrtot)
            deallocate(res)
         endif

         !--clean up after sammy calculation
         if (nmtres.gt.0) call desammy

         !--do energy merge and add resonance cross sections
         if (npend.lt.0) nscr4=-nscr4
         call emerge(nscr1,nscr2,ngrid,nscr3,nrtot,iold,inew,nscr4)

         !--write output tape.
         call recout(iold,nscr4,nrtot)

         !--deallocate arrays for this material
         deallocate(card)
         deallocate(scr)
         deallocate(mfs)
         deallocate(mts)
         deallocate(ncs)
         if (allocated(tr)) deallocate(tr)
         if (allocated(ti)) deallocate(ti)
         if (allocated(sunr)) deallocate(sunr)
         !-- careful, unit number must match nscrl from resxs.
         inquire(unit=16,exist=there)
         if (there) call closz(16)
      else
         deallocate(enode)
      endif
   enddo

   !--reconstruction complete.
   call atend(nout,0)
   call repoz(nout)
   call closz(nout)
   call closz(nin)
   call timer(time)
   write(nsyso,'(69x,f8.1,''s''/&
     &1x,7(''**********''),''*******'')') time
   return
   end subroutine reconr

   subroutine ruina
   !-------------------------------------------------------------------
   ! Read user input for reconr
   !-------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso,nsyse
   use util   ! provides sigfig
   ! internals
   integer::nw,ic,i,j,nwords,ngrid,n2
   real(kr)::z(17)
   character(4)::tz(17)
   equivalence(tz(1),z(1))
   character(66)::text
   real(kr),parameter::zero=0

   !--read user input
   ncards=0
   ngrid=0
   read(nsysi,*) mata,ncards,ngrid
   if (mata.eq.0) return
   tempr=0
   errmax=-1
   errint=-1
   read(nsysi,*) err,tempr,errmax,errint
   if (errmax.le.zero) errmax=10*err
   if (errmax.lt.err) errmax=err
   if (errint.le.zero) errint=err/20000
   if (ncards.gt.0) then
      nw=17*ncards
      allocate(card(nw))
      ic=-17
      do i=1,ncards
         text=' '
         read(nsysi,*) text
         read(text,'(16a4,a2)') (tz(j),j=1,17)
         ic=ic+17
         do j=1,17
            card(j+ic)=tz(j)
         enddo
      enddo
   else
      ncards=1
      nwords=17
      allocate(card(nwords))
      text=' '
      read(text,'(16a4,a2)') (tz(j),j=1,17)
      do i=1,17
         card(i)=tz(i)
      enddo
   endif
   if (ngrid.gt.0) then
      read(nsysi,*) (enode(i),i=1,ngrid)
      do i=1,ngrid
         enode(i)=sigfig(enode(i),7,0)
      enddo
   endif

   !--write out user input
   write(nsyso,'(/&
     &'' material to be processed ............. '',i10/&
     &'' reconstruction tolerance ............. '',f10.3/&
     &'' reconstruction temperature ........... '',f10.2,''k''/&
     &'' resonance-integral-check tolerance ... '',f10.3/&
     &'' max resonance-integral error ......... '',1p,e10.3)')&
     mata,err,tempr,errmax,errint
   if (ncards.gt.0) then
      write(nsyso,'(/&
        &'' descriptive cards for pendf tape''/1x,66(''-''))')
      n2=17*ncards
      write(nsyso,'(1x,17a4)') (card(i),i=1,n2)
   endif
   write(nsyso,'('' '')')
   if (ngrid.gt.0) then
      write(nsyso,'(&
        &'' no. users energy grid points ......... '',i10)') ngrid
      write(nsyso,'(3x,1p,6e12.5)') (enode(i),i=1,ngrid)
   endif
   nodes=ngrid

   !--construct w table if needed for psi-chi.
   if (tempr.gt.zero) then
      allocate(tr(62,62))
      allocate(ti(62,62))
      call wtab(tr,ti)
   endif
   return
   end subroutine ruina

   subroutine anlyzd(nin,nx)
   !-------------------------------------------------------------------
   ! Analyze the dictionary, initialize quantities for new
   ! dictionary, and select desired redundant reactions.
   !-------------------------------------------------------------------
   use endf ! provides dictio
   use util ! provides error
   ! externals
   integer::nin,nx
   ! internals
   integer::nb,nw,i,j,igam,nxn,mfi,mti,i2,j1,mtsave

   !--read in old dictionary.
   nw=6*nx
   allocate(dict(nw))
   nw=nx
   call dictio(nin,0,0,dict,nb,nw)

   !--select desired redundant reactions.
   igam=0
   nxn=0
   nmtr=0
   mtr4=0
   mtr18=0
   mt103=0
   mt104=0
   mt105=0
   mt106=0
   mt107=0
   if (iverf.ge.6) then
       mpmin=600
       mpmax=649
       mdmin=650
       mdmax=699
       mtmin=700
       mtmax=749
       m3min=750
       m3max=799
       m4min=800
       m4max=849
   else
       mpmin=700
       mpmax=718
       mdmin=720
       mdmax=738
       mtmin=740
       mtmax=758
       m3min=760
       m3max=768
       m4min=780
       m4max=798
   endif
   mtr522=0
   do i=1,nw,6
      if (nmtr.ge.nmtmax)&
        call error('anlyzd','too many redundant reactions.',' ')
      mfi=nint(dict(i+2))
      if (mfi.eq.2) then
         mti=nint(dict(i+3))
         if (mti.eq.151) maxres=12*nint(dict(i+4))
      else if (mfi.eq.3) then
         nxn=nxn+1
         if (nmtr.eq.0) mtr(1)=1
         if (nmtr.eq.0) nmtr=1
         mti=nint(dict(i+3))
         if (mti.eq.19) then
            nmtr=nmtr+1
            mtr(nmtr)=18
            mtr18=1
         endif
         if (mti.ge.51.and.mti.le.91.and.mtr4.eq.0) then
            nmtr=nmtr+1
            mtr(nmtr)=4
            mtr4=1
         endif
         if (mti.ge.mpmin.and.mti.le.mpmax.and.mt103.eq.0) then
            nmtr=nmtr+1
            mtr(nmtr)=103
            mt103=1
         endif
         if (mti.ge.mdmin.and.mti.le.mdmax.and.mt104.eq.0) then
            nmtr=nmtr+1
            mtr(nmtr)=104
            mt104=1
         endif
         if (mti.ge.mtmin.and.mti.le.mtmax.and.mt105.eq.0) then
            nmtr=nmtr+1
            mtr(nmtr)=105
            mt105=1
         endif
         if (mti.ge.m3min.and.mti.le.m3max.and.mt106.eq.0) then
            nmtr=nmtr+1
            mtr(nmtr)=106
            mt106=1
         endif
         if (mti.ge.m4min.and.mti.le.m4max.and.mt107.eq.0) then
            nmtr=nmtr+1
            mtr(nmtr)=107
            mt107=1
         endif
      else if (mfi.eq.10) then
         nxn=nxn+1
      else if (mfi.eq.12) then
         nxn=nxn+1
         mti=nint(dict(i+3))
         if (mti.eq.3) nmtr=nmtr+1
         if (mti.eq.3) mtr(nmtr)=3
      else if (mfi.eq.13) then
         nxn=nxn+1
      else if (mfi.eq.23) then
         nxn=nxn+1
         igam=1
         if (iverf.lt.6) then
            zain=0
            awin=0
         endif
         if (nmtr.eq.0) mtr(1)=501
         if (nmtr.eq.0) nmtr=1
         if (mtr522.eq.0) then
            mti=nint(dict(i+3))
            if (mti.ge.534.and.mti.le.572) then
               nmtr=nmtr+1
               mtr(nmtr)=522
               mtr522=1
            endif
         endif
      endif
   enddo
   do i=1,nmtr
      mtrt(i)=1000000
   enddo
   do i=1,nmtr
      if (mtr(i).eq.1) mtrt(i)=1
      if (mtr(i).eq.501) mtrt(i)=1
   enddo

   !--sort reactions into increasing order.
   if (nmtr.gt.2) then
      i2=nmtr-1
      do i=2,i2
         j1=i+1
         do j=j1,nmtr
            if (mtr(j).le.mtr(i)) then
               mtsave=mtr(i)
               mtr(i)=mtr(j)
               mtr(j)=mtsave
            endif
         enddo
      enddo
   endif

   !--initialize new dictionary.
   deallocate(dict)
   nxn=nxn+4
   allocate(mfs(nxn))
   allocate(mts(nxn))
   allocate(ncs(nxn))
   mfs(1)=1
   mts(1)=451
   ncs(1)=ncards+2
   nxc=1
   if (igam.eq.1) return
   mfs(2)=2
   mts(2)=151
   ncs(2)=4
   nxc=2
   return
   end subroutine anlyzd

   subroutine rdfil2(nrtot,intunr)
   !-------------------------------------------------------------------
   ! Read ENDF File 2 data and store it into res, or for the sammy
   ! options, read the data into special arrays.  Nodes are added
   ! just below and above the boundaries of energy ranges, at the
   ! resonance energies, and at values one half-width above and
   ! below the resonance energies.  In the unresolved range, the
   ! scheme depends on the representation.  For energy independent
   ! parameters, nodes are added at 10 standard values per decade.
   ! for sections with only fission widths dependent on energy, the
   ! given values are used as nodes unless an energy step exceeds a
   ! ratio of 1.26 (see wide).  This is then assumed to be a carry-over
   ! evaluation which violates the current rules, and nodes are
   ! added at the 13-per-decade standard values.  A similar procedure
   ! is followed for energy-dependent parameters.  At the same time,
   ! a special grid of unresolved energies is constructed using the
   ! same rules.  This grid is chosen to preserve the correct inter-
   ! polation of unresolved cross sections for MF2, MT153.
   !-------------------------------------------------------------------
   use endf   ! provides endf routines and variable
   use samm   ! provides rdsammy,ppsammy
   use util   ! provides error
   use mainio ! provides nsyso
   ! externals
   integer::nrtot,intunr
   ! internals
   integer::i,j,nb,nw,nis,mode,ner,lru,lrf,nro,naps,jnow
   real(kr)::e,abn
   real(kr),parameter::therm=.0253e0_kr
   real(kr),parameter::eps=1.e-10_kr
   real(kr),parameter::elow=1.e-5_kr
   real(kr),parameter::elarge=9.e9_kr
   real(kr),parameter::small=1.e-6_kr

   !--initialize.
   ! eresl is lowest resonance lower bound
   ! eresh is highest resonance upper bound
   ! eresr is highest resolved upper bound
   ! eresu is lowest unresolved lower bound
   ! eresm is lowest unresolved upper bound
   eresl=elarge
   eresh=0
   eresr=0
   eresu=elarge
   eresm=elarge

   !--read head record.
   jnow=1
   call contio(nin,0,0,res(jnow),nb,nw)
   if (mfh.ne.2) then
      write(nsyso,'(/'' mat has no file 2.'')')
      nrtot=0
      lrp=0
      call tosend(nin,0,0,res)
      return
   endif
   nis=n1h
   ! preset counters
   nsect=0

   !--do loop over all isotopes.
   do i=1,nis
      call contio(nin,0,0,res(jnow),nb,nw)
      jnow=jnow+nw
      if (jnow.gt.maxres) call error('rdfil2',&
        'res storage exceeded',' ')
      zai=c1h
      abn=c2h
      lfw=l2h
      ner=n1h

      !--do loop over all energy ranges.
      do j=1,ner
         call contio(nin,0,0,res(jnow),nb,nw)
         el=c1h
         eh=c2h
         lru=l1h
         lrf=l2h

         !--add el and eh to list of nodes.
         !--shade nodes to prevent discontinuities.
         if (el.lt.eresl) eresl=el
         if (eh.gt.eresh) eresh=eh
         if (lru.le.1.and.eh.gt.eresr) eresr=eh
         if (lru.ne.0) then
            if (abs(el-elow).le.small) then
               nodes=nodes+1
               enode(nodes)=el
            else
               nodes=nodes+1
               enode(nodes)=sigfig(el,7,-1)
               nodes=nodes+1
               enode(nodes)=sigfig(el,7,+1)
               if (lru.eq.2) then
                  if (el.lt.eresu) eresu=el
                  if (nunr.ge.maxunr) call error('rdfil2',&
                    'storage in eunr exceeded.',' ')
                  nunr=nunr+1
                  eunr(nunr)=sigfig(el,7,-1)
                  nunr=nunr+1
                  eunr(nunr)=sigfig(el,7,+1)
               endif
            endif
            nodes=nodes+1
            enode(nodes)=sigfig(eh,7,-1)
            nodes=nodes+1
            enode(nodes)=sigfig(eh,7,+1)
            if (lru.eq.2) then
               if (eh.lt.eresm) eresm=eh
               if (nunr.ge.maxunr) call error('rdfil2',&
                 &'storage in eunr exceeded.',' ')
               nunr=nunr+1
               eunr(nunr)=sigfig(eh,7,-1)
               nunr=nunr+1
               eunr(nunr)=sigfig(eh,7,+1)
            endif
         endif

         !--check for energy-dependent scattering radius
         !--or energy-dependent fission widths
         nro=n1h
         if (lru.eq.1) then
            res(jnow+4)=nro
         else
            res(jnow+4)=lfw
         endif
         ! get scattering radius flag
         naps=n2h
         res(jnow+5)=naps
         ! test formalism
         mode=lrf+10*(lru-1)
         if (lru.eq.0) mode=0
         if (lru.eq.2) lrp=3
         ! store isotope-section table data
         nsect=nsect+1
         isot(nsect)=i
         abnt(nsect)=abn
         elt(nsect)=el
         eht(nsect)=eh
         modet(nsect)=mode
         ! ibaset stores location of isotope card
         ibaset(nsect)=jnow
         jnow=jnow+nw
         if (jnow.gt.maxres) call error('rdfil2',&
           'res storage exceeded',' ')

         !--process the rest of this part
         !--according to its formalism
         intunr=5
         if (mode.le.3.and.nmtres.eq.0) then
            call rdf2bw(nin,jnow,lrf,nro,naps,mode)
         else if (mode.le.3.and.nmtres.gt.0) then
            call rdsammy(nin,j,jnow,nro,naps,mode,el,eh,&
              enode,nodes,nodmax,spin,ascat,res,maxres)
            call ppsammy(j,ncoef,nresp)
         else if (mode.eq.4) then
            call rdf2aa(nin,jnow,mode)
         else if (mode.eq.6) then
            call rdf2hy(nin,jnow)
         else if (mode.eq.7.and.nmtres.gt.0) then
            call rdsammy(nin,j,jnow,nro,naps,mode,el,eh,&
              enode,nodes,nodmax,spin,ascat,res,maxres)
            call ppsammy(j,ncoef,nresp)
         else if (mode.eq.11.and.lfw.eq.0) then
            call rdf2u0(nin,jnow,nro)
         else if (mode.eq.11.and.lfw.eq.1) then
            call rdf2u1(nin,jnow,nro)
         else if (mode.eq.12) then
            call rdf2u2(nin,jnow,nro,intunr)
         else
            call error('rdfil2','illegal resonance mode',' ')
         endif

         !--add .0253 ev to the list of nodes
         if (el.lt.therm) then
            nodes=nodes+1
            enode(nodes)=therm
         endif

      !--continue loop over energy ranges
      enddo

   !--continue loop over isotopes
   enddo
   if (eresr.lt.eresl) eresr=eresl
   if (eresr.gt.eresh) eresr=eresh

   !--finish off file and sort nodes into order.
   call tofend(nin,0,0,res(jnow))
   nrtot=jnow
   call order(enode,nodes)
   if (nunr.eq.0) return

   !--sort unresolved nodes into order.
   !--mark resolved-unresolved overlaps.
   call order(eunr,nunr)
   j=0
   el=0
   do i=1,nunr
      e=eunr(i)
      if (e.gt.(1-eps)*eresu.and.e.lt.(1+eps)*eresh&
        .and.e.ge.el) then
         j=j+1
         if (e.lt.eresr) e=-e
         eunr(j)=e
         el=sigfig(abs(e),7,+1)
      endif
   enddo
   nunr=j
   if (eresu.lt.eresr) then
      write(nsyso,'(/&
        &'' mat has resolved-unresolved overlap in range''/&
        &5x,1p,2e12.3)') eresu,eresr
   endif
   if (eresm.lt.eresh) then
      write(nsyso,'(/&
        &'' mat has unresolved-smooth overlap in range''/&
        &5x,1p,2e12.3)') eresm,eresh
   endif
   return
   end subroutine rdfil2

   subroutine rdf2bw(nin,jnow,lrf,nro,naps,mode)
   !-------------------------------------------------------------------
   ! Read in resonance data for single-level Breit-Wigner,
   ! multi-level Breit-Wigner, and Reich-Moore formats.
   !-------------------------------------------------------------------
   use util    ! provides error,mess
   use endf    ! provides endf routines and variables
   use physics ! provides amassn,amu,hbar,ev
   ! externals
   integer::nin,jnow,lrf,nro,naps,mode
   ! internals
   integer::nb,nw,nls,lad,nrs,ncyc,ll,lp,n,jj,l,ndig
   real(kr)::cwaven,awri,ascatl,arat,aw,ra,ral,qx
   real(kr)::hw,ehalf,rho,rhoc,ser,per,sec,pec
   real(kr),parameter::rc1=.123e0_kr
   real(kr),parameter::rc2=.08e0_kr
   real(kr),parameter::third=.333333333e0_kr
   real(kr),parameter::gxmin=1.0e-5_kr
   real(kr),parameter::zero=0
   cwaven=sqrt(2*amassn*amu*ev)*1.e-12_kr/hbar

   !--check for energy-dependent scattering radius
   if (nro.ne.0) then
      if (mode.ne.2) call error('rdf2bw',&
        &'energy-dep scattering length',&
        &'only works for mlbw')
      call tab1io(nin,0,0,res(jnow),nb,nw)
      jnow=jnow+nw
      if (jnow.gt.maxres) then
         call error('rdfil2','res storage exceeded',' ')
      endif
   endif

   !--read some parameters
   call contio(nin,0,0,res(jnow),nb,nw)
   jnow=jnow+nw
   if (jnow.gt.maxres) call error('rdf2bw','res storage exceeded',' ')
   nls=n1h
   lad=0
   if (mode.eq.3) lad=l1h
   if (lad.ne.0) call mess('rdf2bw',&
     &'calculation of angular distribution not installed.',' ')
   lssf=0
   spin=c1h
   ascat=c2h

   !--do loop over all l states
   do l=1,nls
      call listio(nin,0,0,res(jnow),nb,nw)
      awri=res(jnow)
      ascatl=ascat
      if (lrf.eq.3) ascatl=res(jnow+1)
      if (ascatl.eq.zero) ascatl=ascat
      arat=awri/(awri+1)
      aw=amassn*awri
      ra=rc1*aw**third+rc2
      ral=ra
      if (naps.eq.1) then
         ral=ascatl
         ra=ascat
      endif
      jj=jnow+nw
      do while (nb.ne.0)
         call moreio(nin,0,0,res(jj),nb,nw)
         jj=jj+nw
         if (jj.gt.maxres) call error('rdf2bw',&
           'res storage exceeded',' ')
      enddo
      nrs=nint(res(jnow+5))
      ncyc=nint(res(jnow+4))/nrs
      ll=nint(res(jnow+2))
      qx=res(jnow+1)
      lrx=nint(res(jnow+3))
      if (lrx.ne.0) then
         nsig=5
      endif
      jnow=jnow+6

      !--loop over resonances in this sequence
      do n=1,nrs
         if (res(jnow).gt.el.and.res(jnow).le.eh) then
            if (nodes.ge.nodmax) call error('rdf2bw',&
              'storage in enode exceeded.',' ')
            nodes=nodes+1
            hw=res(jnow+2)/2
            if (mode.eq.3)&
              hw=hw+(res(jnow+3)+abs(res(jnow+4))+abs(res(jnow+5)))/2
            ndig=5
            if (res(jnow).gt.zero) ndig=2+nint(log10(res(jnow)/(hw/10)))
            if (ndig.lt.5) ndig=5
            if (ndig.gt.9) ndig=9
            enode(nodes)=sigfig(res(jnow),ndig,0)
            if ((res(jnow)+hw).lt.eh) then
               if (nodes.ge.nodmax) call error('rdf2bw',&
                 &'storage in enode exceeded.',' ')
               ehalf=res(jnow)+hw
               nodes=nodes+1
               enode(nodes)=sigfig(ehalf,ndig,0)
            endif
            if ((res(jnow)-hw).gt.el) then
               if (nodes.ge.nodmax) call error('rdf2bw',&
                 &'storage in enode exceeded.',' ')
               ehalf=res(jnow)-hw
               nodes=nodes+1
               enode(nodes)=sigfig(ehalf,ndig,0)
            endif
         endif
         rho=cwaven*arat*sqrt(abs(res(jnow)))*ral
         call facts(ll,rho,ser,per)
         res(jj)=ser
         res(jj+1)=per
         res(jj+2)=0
         res(jj+3)=0
         res(jj+4)=0
         res(jj+5)=0
         if (lrx.ne.0) then
            rhoc=cwaven*arat*sqrt(abs(res(jnow)+qx/arat))*ra
            lp=ll
            ! set competing l' value for a 2+ inelastic residual (u238)
              if (ll.eq.0) lp=2
              if (ll.eq.2) lp=0
            ! use ll to compute gx0
            call facts(lp,rhoc,sec,pec)
            res(jj+2)=res(jnow+2)-res(jnow+3)-res(jnow+4)-res(jnow+5)
            if (res(jj+2).lt.gxmin*res(jnow+2)) res(jj+2)=0
            if (res(jnow).lt.-qx/arat) res(jj+2)=0
            res(jj+3)=pec
         endif
         jj=jj+6
         if (jj.gt.maxres) call error('rdf2bw',&
           'res storage exceeded',' ')
         jnow=jnow+ncyc
      enddo
      jnow=jj
   enddo
   return
   end subroutine rdf2bw

   subroutine rdf2aa(nin,jnow,mode)
   !-------------------------------------------------------------------
   ! Read in resonance data for the Adler-Adler format.
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variable
   use util ! provides error
   ! externals
   integer::nin,jnow,mode
   ! internals
   integer::nb,nw,nls,lssf,l,njs,n,nen,jen,nloc,ien,ndig
   real(kr)::ener,hw,ehalf
   real(kr),parameter::zero=0

   !--read some parameters
   call contio(nin,0,0,res(jnow),nb,nw)
   jnow=jnow+nw
   if (jnow.gt.maxres) call error('rdf2aa','res storage exceeded',' ')
   nls=n1h
   lssf=0
   spin=c1h
   ascat=c2h

   !--read adler-adler background corrections.
   call listio(nin,0,0,res(jnow),nb,nw)
   jnow=jnow+nw
   do while (nb.ne.0)
      if (jnow.gt.maxres) call error('rdf2aa',&
        'res storage exceeded',' ')
      call moreio(nin,0,0,res(jnow),nb,nw)
      jnow=jnow+nw
   enddo

   !--do loop over all l states
   do l=1,nls
      call contio(nin,0,0,res(jnow),nb,nw)
      jnow=jnow+nw
      if (jnow.gt.maxres) call error('rdf2aa',&
        &'res storage exceeded',' ')
      njs=n1h
      do n=1,njs
         call listio(nin,0,0,res(jnow),nb,nw)
         nen=n2h
         jen=12
         nloc=jnow-6
         jnow=jnow+nw
         do while (nb.ne.0)
            if (jnow.gt.maxres) call error('rdf2aa',&
              &'res storage exceeded',' ')
            call moreio(nin,0,0,res(jnow),nb,nw)
            jnow=jnow+nw
         enddo
         if (n.eq.1.and.l.eq.1) then
            do ien=1,nen
               ener=res(nloc+ien*jen)
               if (ener.gt.el.and.ener.le.eh) then
                  if (nodes.ge.nodmax) call error('rdf2aa',&
                    &'storage in enode exceeded.',' ')
                  nodes=nodes+1
                  hw=res(nloc+ien*jen+1)
                  ndig=5
                  if (ener.gt.zero) ndig=2+nint(log10(ener/(hw/10)))
                  if (ndig.lt.5) ndig=5
                  if (ndig.gt.9) ndig=9
                  enode(nodes)=sigfig(ener,ndig,0)
                  if ((ener+hw).le.eh) then
                     if (nodes.ge.nodmax) call error('rdf2aa',&
                       &'storage in enode exceeded.',' ')
                     ehalf=ener+hw
                     nodes=nodes+1
                     enode(nodes)=sigfig(ehalf,ndig,0)
                  endif
                  if (ener-hw.gt.el) then
                     if (nodes.ge.nodmax) call error('rdf2aa',&
                       &'storage in enode exceeded.',' ')
                     ehalf=ener-hw
                     nodes=nodes+1
                     enode(nodes)=sigfig(ehalf,ndig,0)
                  endif
               endif
            enddo
         endif
      enddo
   enddo
   return
   end subroutine rdf2aa

   subroutine rdf2hy(nin,jnow)
   !-------------------------------------------------------------------
   ! Read in resonance data for the hybrid R-matrix format.
   !-------------------------------------------------------------------
   use util ! provides error,mess
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,jnow
   ! internals
   integer::nb,nw,nls,lad,l,ns,nj,ir,ll,ig,ndig
   integer::nire,ncre,ncr,lil,nss,njs,lbk,lps,nlsj
   real(kr)::hw,ehalf,enow
   real(kr),parameter::zero=0

   !--read some parameters
   call contio(nin,0,0,res(jnow),nb,nw)
   jnow=jnow+nw
   if (jnow.gt.maxres) call error('rdf2hy','res storage exceeded',' ')
   nls=n1h
   lad=l1h
   if (lad.ne.0) call mess('rdf2hy',&
     &'calculation of angular distribution not installed.',' ')
   lssf=0
   spin=c1h
   ascat=c2h

   !--read additional hybrid-r parameters
   call contio(nin,0,0,res(jnow),nb,nw)
   nire=n1h
   ncre=n2h
   if (nire+ncre.gt.0) call mess('rdf2hy',&
     &'hybrid competing reactions not yet added to file 3.',' ')
   jnow=jnow+nw
   call contio(nin,0,0,res(jnow),nb,nw)
   jnow=jnow+nw
   call listio(nin,0,0,res(jnow),nb,nw)
   jnow=jnow+nw
   if (ncre.ne.0) then
      do ncr=1,ncre
         do lil=1,4
            call tab1io(nin,0,0,res(jnow),nb,nw)
            jnow=jnow+nw
         enddo
      enddo
   endif
   if (nls.eq.0) return

   !--do loop over all l states
   do l=1,nls
      call contio(nin,0,0,res(jnow),nb,nw)
      nss=n1h
      jnow=jnow+nw
      do ns=1,nss
         call contio(nin,0,0,res(jnow),nb,nw)
         njs=n1h
         jnow=jnow+nw
         do nj=1,njs
            call listio(nin,0,0,res(jnow),nb,nw)
            lbk=l1h
            lps=l2h
            nlsj=n2h
            ll=jnow+6
            do ir=1,nlsj
               hw=0
               do ig=1,7
                  hw=hw+res(ll+ig)
               enddo
               hw=hw/2
               enow=res(ll)
               ndig=5
               if (enow.gt.zero) ndig=2+nint(log10(enow/(hw/10)))
               if (ndig.lt.5) ndig=5
               if (ndig.gt.9) ndig=9
               if (nodes+3.ge.nodmax) call error('rdf2hy',&
                 &'storage in enode exceeded.',' ')
               if (enow.gt.el.and.enow.le.eh) then
                  nodes=nodes+1
                  enode(nodes)=sigfig(enow,ndig,0)
                  if (enow-hw.gt.el) then
                     ehalf=enow-hw
                     nodes=nodes+1
                     enode(nodes)=sigfig(ehalf,ndig,0)
                  endif
                  if (enow+hw.lt.eh) then
                     ehalf=enow+hw
                     nodes=nodes+1
                     enode(nodes)=sigfig(ehalf,ndig,0)
                  endif
               endif
               ll=ll+12
            enddo
            jnow=jnow+nw
            if (lbk.ne.0) then
               call tab1io(nin,0,0,res(jnow),nb,nw)
               jnow=jnow+nw
               call tab1io(nin,0,0,res(jnow),nb,nw)
               jnow=jnow+nw
            endif
            if (lps.ne.0) then
               call tab1io(nin,0,0,res(jnow),nb,nw)
               jnow=jnow+nw
               call tab1io(nin,0,0,res(jnow),nb,nw)
               jnow=jnow+nw
            endif
         enddo
      enddo
   enddo
   return
   end subroutine rdf2hy

   subroutine rdf2u0(nin,jnow,nro)
   !-------------------------------------------------------------------
   ! Read in resonance data for the energy-independent
   ! unresolved format (mode=11, lfw=0).
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use util ! provides error
   ! externals
   integer::nin,jnow,nro
   ! internals
   integer::nb,nw,nls,l,jj,ii
   real(kr)::ener,enut
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

   !--read some parameters
   if (nro.eq.1) then
      call tab1io(nin,0,0,res(jnow),nb,nw)
      res(jnow+1)=-float(nro)
      jj=jnow+nw
      do while (nb.ne.0)
         call moreio(nin,0,0,res(jj),nb,nw)
         jj=jj+nw
         if (jj.gt.maxres) then
            call error('rdf2u0','storage in res exceeded',' ')
         endif
      enddo
      jnow=jj
   endif
   call contio(nin,0,0,res(jnow),nb,nw)
   jnow=jnow+nw
   if (jnow.gt.maxres) call error('rdf2u0','res storage exceeded',' ')
   nls=n1h
   lssf=l1h
   spin=c1h
   ascat=c2h

   !--do loop over all l states
   do l=1,nls
      call listio(nin,0,0,res(jnow),nb,nw)
      jj=jnow+nw
      do while (nb.ne.0)
         call moreio(nin,0,0,res(jj),nb,nw)
         jj=jj+nw
         if (jj.gt.maxres) call error('rdf2u0',&
           &'res storage exceeded',' ')
      enddo
      jnow=jj

      !--add extra unresolved nodes
      ener=el
      do while (ener.lt.eh)
         enut=0
         ii=0
         do while (ii.le.ngridu.and.enut.lt.ener+ener/100)
            ii=ii+1
            enut=egridu(ii)
         enddo
         ener=enut
         if (ener.lt.eh) then
            if (nodes.ge.nodmax)&
              call error('rdf2u0','storage in enode exceeded.',' ')
            nodes=nodes+1
            enode(nodes)=sigfig(ener,7,0)
            if (nunr.ge.maxunr)&
              call error('rdf2u0','storage in eunr exceeded.',' ')
            nunr=nunr+1
            eunr(nunr)=sigfig(ener,7,0)
         endif
      enddo
   enddo
   return
   end subroutine rdf2u0

   subroutine rdf2u1(nin,jnow,nro)
   !-------------------------------------------------------------------
   ! Read in resonance data for the energy-independent
   ! unresolved format with energy-dependent fission widths
   ! (mode=11, lfw=1).
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use util ! provides error
   ! externals
   integer::nin,jnow,nro
   ! internals
   integer::nb,nw,nloc,ifill,ii,l,njs,ne,nls,n,jj
   real(kr)::enut,ener,enex
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
   real(kr),parameter::wide=1.26e0_kr

   !--read some parameters
   if (nro.eq.1) then
      call tab1io(nin,0,0,res(jnow),nb,nw)
      res(jnow+1)=-float(nro)
      jj=jnow+nw
      do while (nb.ne.0)
         call moreio(nin,0,0,res(jj),nb,nw)
         jj=jj+nw
         if (jj.gt.maxres) then
            call error('rdf2u1','storage in res exceeded',' ')
         endif
      enddo
      jnow=jj
   endif
   call listio(nin,0,0,res(jnow),nb,nw)
   nloc=jnow+5
   jnow=jnow+nw
   if (jnow.gt.maxres) call error('rdf2u1','res storage exceeded',' ')
   do while (nb.ne.0)
      call moreio(nin,0,0,res(jnow),nb,nw)
   enddo
   lssf=l1h
   ne=n1h
   nls=n2h
   spin=c1h
   ascat=c2h

   !--read in the energy-dependent fission widths
   do n=1,ne
      ener=res(nloc+n)
      if (ener.ge.el.and.ener.lt.eh) then
         enex=res(nloc+n+1)
         ifill=0
         if (enex.gt.wide*ener) ifill=1
         if (nodes.ge.nodmax)&
           &call error('rdf2u1','storage in enode exceeded.',' ')
         nodes=nodes+1
         enode(nodes)=sigfig(ener,7,0)
         if (nunr.ge.maxunr)&
           &call error('rdf2u1','storage in eunr exceeded.',' ')
         nunr=nunr+1
         eunr(nunr)=sigfig(ener,7,0)
         do while (ifill.ne.0)
            enut=0
            ii=0
            do while (ii.lt.ngridu.and.enut.lt.ener+ener/100)
               ii=ii+1
               enut=egridu(ii)
            enddo
            ener=enut
            if (ener.ge.enex) ifill=0
            if (ifill.eq.1) then
               if (nodes.ge.nodmax) call error('rdf2u1',&
                 &'storage in enode exceeded.',' ')
               nodes=nodes+1
               enode(nodes)=sigfig(ener,7,0)
               if (nunr.ge.maxunr) call error('rdf2u1',&
                 &'storage in eunr exceeded.',' ')
               nunr=nunr+1
                  eunr(nunr)=sigfig(ener,7,0)
            endif
         enddo
      endif
   enddo

   !--do loop over all l states
   do l=1,nls
      call contio(nin,0,0,res(jnow),nb,nw)
      jnow=jnow+nw
      if (jnow.gt.maxres) call error('rdf2u1',&
        &'res storage exceeded',' ')
      njs=n1h
      do n=1,njs
         call listio(nin,0,0,res(jnow),nb,nw)
         nloc=jnow+6
         jnow=jnow+nw
         do while (nb.ne.0)
            if (jnow.gt.maxres) call error('rdf2u1',&
              'res storage exceeded',' ')
            call moreio(nin,0,0,res(jnow),nb,nw)
            jnow=jnow+nw
         enddo
      enddo
   enddo
   return
   end subroutine rdf2u1

   subroutine rdf2u2(nin,jnow,nro,intunr)
   !-------------------------------------------------------------------
   ! Read in resonance data for the energy-dependent
   ! unresolved format (mode=12).
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use util ! provides error
   ! externals
   integer::nin,jnow,nro,intunr
   ! internals
   integer::nb,nw,nls,l,njs,n,nen,jen,nloc,ien,ifill,ii,jj
   real(kr)::ener,enex,enut
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
   real(kr),parameter::wide=1.26e0_kr

   !--read some parameters
   if (nro.eq.1) then
      call tab1io(nin,0,0,res(jnow),nb,nw)
      res(jnow+1)=-float(nro)
      jj=jnow+nw
      do while (nb.ne.0)
         call moreio(nin,0,0,res(jj),nb,nw)
         jj=jj+nw
         if (jj.gt.maxres) then
            call error('rdf2u2','storage in res exceeded',' ')
         endif
      enddo
      jnow=jj
   endif
   call contio(nin,0,0,res(jnow),nb,nw)
   jnow=jnow+nw
   if (jnow.gt.maxres) call error('rdf2u2','res storage exceeded',' ')
   nls=n1h
   lssf=l1h
   spin=c1h
   ascat=c2h

   !--do loop over all l states
   do l=1,nls
      call contio(nin,0,0,res(jnow),nb,nw)
      jnow=jnow+nw
      if (jnow.gt.maxres) call error('rdf2u2',&
        &'res storage exceeded',' ')
      njs=n1h
      do n=1,njs
         call listio(nin,0,0,res(jnow),nb,nw)
         intunr=l1h
         nen=n2h
         jen=6
         nloc=jnow+6
         jnow=jnow+nw
         do while (nb.ne.0)
            if (jnow.gt.maxres) call error('rdf2u2',&
              &'res storage exceeded',' ')
            call moreio(nin,0,0,res(jnow),nb,nw)
            jnow=jnow+nw
         enddo
         if (n.eq.1.and.l.eq.1) then
            do ien=1,nen
               ener=res(nloc+(ien-1)*jen)
               if (ener.ge.el.and.ener.lt.eh) then
                  enex=res(nloc+(ien-1)*jen+jen)
                  ifill=0
                  if (enex.gt.wide*ener) ifill=1
                  if (nodes.ge.nodmax) call error('rdf2u2',&
                    &'storage in enode exceeded.',' ')
                  nodes=nodes+1
                  enode(nodes)=sigfig(ener,7,0)
                  if (nunr.ge.maxunr) call error('rdf2u2',&
                    &'storage in eunr exceeded.',' ')
                  nunr=nunr+1
                  eunr(nunr)=sigfig(ener,7,0)
                  do while (ifill.ne.0)
                     enut=0
                     ii=0
                     do while (ii.lt.ngridu.and.&
                       &enut.lt.ener+ener/1000)
                        ii=ii+1
                        enut=egridu(ii)
                     enddo
                     ener=enut
                     if (ener.ge.enex) ifill=0
                     if (ifill.eq.1) then
                        if (nodes.ge.nodmax) call error('rdf2u2',&
                          &'storage in enode exceeded.',' ')
                        nodes=nodes+1
                        enode(nodes)=sigfig(ener,7,0)
                        if (nunr.ge.maxunr) call error('rdf2u2',&
                          &'storage in eunr exceeded.',' ')
                        nunr=nunr+1
                        eunr(nunr)=sigfig(ener,7,0)
                     endif
                  enddo
               endif
            enddo
         endif
      enddo
   enddo
   return
   end subroutine rdf2u2

   subroutine order(x,n)
   !-------------------------------------------------------------------
   ! Sort the n elements of x into ascending order
   ! removing any duplicate elements
   !-------------------------------------------------------------------
   ! externals
   integer::n
   real(kr)::x(*)
   ! internals
   integer::m,i,j,k
   real(kr)::tsave
   real(kr),parameter::small=1.e-10_kr

   if (n.le.2) return
   m=n
   i=0
   do while (i.lt.m-1)
      i=i+1
      j=i
      do while (j.lt.m)
         j=j+1
         if (x(j).lt.x(i)) then
            tsave=x(j)
            x(j)=x(i)
            x(i)=tsave
         endif
      enddo
      if (i.gt.1) then
         if (abs(x(i)-x(i-1)).le.small*x(i)) then
            m=m-1
            if (i.ge.m) then
               n=m
               return
            endif
            do k=i,m
               x(k)=x(k+1)
            enddo
            i=i-1
         endif
      endif
   enddo
   if (abs(x(m)-x(m-1)).le.small*x(m)) m=m-1
   n=m
   return
   end subroutine order

   subroutine facts (l,rho,se,pe)
   !-------------------------------------------------------------------
   ! Calculates penetration and shift factors
   !-------------------------------------------------------------------
   ! externals
   integer::l
   real(kr)::rho,se,pe
   ! internals
   real(kr)::den,r2,r4,r6,r8

   r2=rho*rho
   if (l.eq.0) then
      se=0
      pe=rho
   else if (l.eq.1) then
      den=1+r2
      pe=r2*rho/den
      se=-1/den
   else if (l.eq.2) then
      r4=r2*r2
      den=3*r2+r4+9
      pe=r4*rho/den
      se=-(18+3*r2)/den
   else if (l.eq.3) then
      r4=r2*r2
      r6=r4*r2
      den=225+45*r2+6*r4+r6
      pe=r6*rho/den
      se=-(675+90*r2+6*r4)/den
   else if (l.eq.4) then
      r4=r2*r2
      r6=r4*r2
      r8=r4*r4
      den=11025+1575*r2+135*r4+10*r6+r8
      pe=r8*rho/den
      se=-(44100+4725*r2+270*r4+10*r6)/den
   endif
   return
   end subroutine facts

   subroutine genunr(nin,intunr)
   !-------------------------------------------------------------------
   ! Compute infinitely dilute cross sections on the unresolved energy
   ! grid for later retrieval with sigunr.  The table is prepared in
   ! the format of MF2/MT152.
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use util ! provides sigfig
   ! externals
   integer::nin,intunr
   !internals
   real(kr)::sigu(4)
   integer::nb,nw,l,i,j,k,ix,idis,mode,mt1
   real(kr)::e,enext,abn,bkg
   real(kr),dimension(:),allocatable::scr
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::zero=0

   !--reserve space for unresolved table
   nw=6*nunr+13
   allocate(sunr(nw))

   !--store cont and list constants
   sunr(1)=za
   sunr(2)=awr
   sunr(3)=lssf
   sunr(4)=0
   sunr(5)=0
   sunr(6)=intunr
   sunr(7)=tempr
   sunr(8)=0
   sunr(9)=5
   sunr(10)=1
   sunr(11)=1+6*nunr
   sunr(12)=nunr
   sunr(13)=big

   !--compute and store unresolved part of cross sections
   e=0
   call csunr1(e,sigu)
   call csunr2(e,sigu)
   l=14
   do i=1,nunr
      e=eunr(i)
      sunr(l)=e
      e=abs(e)
      do k=1,4
         sunr(l+k)=0
      enddo
      do j=1,nsect
         mode=modet(j)
         if (mode.ge.11.and.e.ge.elt(j).and.e.lt.eht(j)) then
            abn=abnt(j)
            isect=j
            if (mode.eq.11) call csunr1(e,sigu)
            if (mode.eq.12) call csunr2(e,sigu)
            do k=1,4
               sunr(l+k)=sunr(l+k)+abn*sigu(k)
               sunr(l+k)=sigfig(sunr(l+k),7,0)
            enddo
         endif
      enddo
      sunr(l+5)=sunr(l+1)
      l=l+6
   enddo
   if (lssf.ne.0) return

   !--add on unresolved background from mf3 on endf tape
   !--any mf3 background in a range of resolved-unresolved
   !--overlap is omitted.
   allocate(scr(npage+50))
   mt1=0
   call findf(mata,3,0,nin)
   mfh=3
   do while (mfh.ne.0)
      call contio(nin,0,0,scr,nb,nw)
      if (mfh.ne.0) then
         if (mt1.eq.0) mt1=mth
         ix=0
         if (mth.eq.1) ix=1
         if (mth.eq.2) ix=2
         if (mth.eq.18) ix=3
         if (mth.eq.102) ix=4
         if (ix.ne.0) then
            e=0
            call gety1(e,enext,idis,bkg,nin,scr)
            l=14
            do i=1,nunr
               e=sunr(l)
               if (e.ge.zero) then
                  call gety1(e,enext,idis,bkg,nin,scr)
                  sunr(l+ix)=sunr(l+ix)+bkg
                  sunr(l+ix)=sigfig(sunr(l+ix),7,0)
                  if (ix.eq.1) then
                     sunr(l+5)=sunr(l+5)+bkg
                     sunr(l+5)=sigfig(sunr(l+5),7,0)
                  endif
               endif
               l=l+6
            enddo
         endif
         call tosend(nin,0,0,scr)
      endif
   enddo
   call findf(mata,3,mt1,nin)
   deallocate(scr)
   return
   end subroutine genunr

   subroutine sigunr(e,sigu)
   !-------------------------------------------------------------------
   ! Retrieve infinitely dilute unresolved cross sections
   ! from table produced by genunr.
   !-------------------------------------------------------------------
   use endf ! provides terp1
   ! externals
   real(kr)::e,sigu(4)
   ! internals
   integer::intunr,l,l1,l2,i
   real(kr)::en,e1,e2

   intunr=nint(sunr(6))
   l=14
   en=0
   i=1
   do while (i.lt.nunr.and.e.ge.en)
      i=i+1
      l2=l+6*(i-1)
      en=abs(sunr(l2))
   enddo
   l1=l2-6
   e1=abs(sunr(l1))
   e2=abs(sunr(l2))
   do i=1,4
      if (e.ge.e1.and.e.le.e2) then
         call terp1(e1,sunr(l1+i),e2,sunr(l2+i),e,sigu(i),intunr)
      else
         sigu(i)=0
      endif
   enddo
   return
   end subroutine sigunr

   subroutine lunion(nin,nout,ngrid,ngr)
   !-------------------------------------------------------------------
   ! Copy desired reactions from MF3 and MF13 on nin to nout.
   ! Check each reaction against the energy grid and add any points
   ! needed to unionize and linearize all desired reactions to
   ! the given tolerance.  Remove singularities by moving first
   ! point down and second point up in last sig fig.
   !-------------------------------------------------------------------
   use util   ! provides openz,repoz,error,timer
   use endf   ! provides endf routines and variables
   use mainio ! provides nsyso
   ! externals
   integer::nin,nout,ngrid,ngr
   ! internals
   integer::ndim,iold,inew,i,j,mfl,nb,nw,nss,iss
   integer::l,lsave,in,ig,ir,jr,lr,ibase,nbta,inta,intn
   integer::ipwr,inx,ngn,isave,ngneg,jg,igt,ngpos,npr
   character(40)::text
   real(kr)::stpmax,awrx,eg,qx,thrx,thrxx,et
   real(kr)::er,sr,ernext,srnext,enl,test,snl,erl,srl
   real(kr)::en,sn,xm,ym,yl,errn,ent,egl,time
   real(kr)::aa(1)
   real(kr),dimension(:),allocatable::bufn,bufo
   real(kr),dimension(:),allocatable::x,y
   real(kr),dimension(:),allocatable::scr
   real(kr),parameter::elow=1.e-5_kr
   real(kr)::elim=.99e6_kr
   real(kr),parameter::up=1.001e0_kr
   real(kr),parameter::dn=0.999e0_kr
   real(kr),parameter::therm=.0253e0_kr
   real(kr),parameter::small=1.e-9_kr
   real(kr),parameter::ten=10.e0_kr
   real(kr),parameter::ssmall=1.e-30_kr
   real(kr),parameter::stpmin=1.001e0_kr
   real(kr),parameter::ovfact=5.3e0_kr
   real(kr),parameter::trange=.4999e0_kr
   real(kr),parameter::emax=19.e6_kr
   real(kr),parameter::zero=0

   !--initialize.
   if (eresr.lt.elim) elim=eresr
   ndim=50
   allocate(bufo(nbuf))
   allocate(bufn(nbuf))
   allocate(x(ndim))
   allocate(y(ndim))
   allocate(scr(npage+50))

   ! this value fits 1/v to within err
   stpmax=1+sqrt(ovfact*err)

   !--copy nodes from the global area to *old* tape
   iold=14
   inew=15
   call openz(-iold,1)
   call openz(-inew,1)
   call repoz(-iold)
   call repoz(-inew)
   if (nodes.gt.0) then
      do i=1,nodes
         j=i
         if (i.eq.nodes) j=-i
         aa(1)=enode(i)
         call loada(j,aa,1,iold,bufo,nbuf)
      enddo
   endif
   ngo=nodes

   !--find next reaction to be processed
   !--last card read from nin was fend card
   call repoz(nout)
   nsh=0
   mfl=0
  110 continue
   call contio(nin,0,0,scr,nb,nw)
   nss=1
   if (mfh.eq.10) nss=n1h
   iss=nss
   if (math.gt.0) go to 120
   call amend(nout,0)
   call atend(nout,0)
   go to 410
  120 continue
   if (mfh.gt.0) go to 130
   if (mfl.eq.3.or.mfl.eq.10.or.mfl.eq.12.or.mfl.eq.13.or.mfl.eq.23)&
     call afend(nout,0)
   mfl=0
   go to 110
   ! process mf3, mf10, mf13, and mf23 only.
  130 continue
   if (mfh.eq.3.or.mfh.eq.10.or.mfh.eq.13.or.mfh.eq.23) go to 140
   if (mfh.eq.12) go to 140
   call tofend(nin,0,0,scr)
   go to 110
   ! skip redundant reactions and non-cross-sections in n files
   ! remove discrete photon records from file 13.
  140 continue
   mfl=mfh
   if (mth.eq.501) go to 150
   if (mth.eq.522.and.mtr522.gt.0) go to 150
   if (mth.eq.460) go to 150
   if (mfh.eq.12.and.nint(scr(3)).ne.1) go to 150
   if (mfh.eq.12) scr(5)=1
   if (mfh.eq.13) scr(5)=1
   if (mfh.ne.3) go to 180
   if (mth.eq.1.or.mth.eq.3) go to 150
   if (mth.eq.4.and.mtr4.gt.0) go to 150
   if (mth.eq.103.and.mt103.gt.0) go to 150
   if (mth.eq.104.and.mt104.gt.0) go to 150
   if (mth.eq.105.and.mt105.gt.0) go to 150
   if (mth.eq.106.and.mt106.gt.0) go to 150
   if (mth.eq.107.and.mt107.gt.0) go to 150
   if (mth.eq.101) go to 150
   if (mth.eq.120) go to 150
   if (mth.eq.151) go to 150
   if ((mth.ge.251.and.mth.le.300).and.mth.ne.261) go to 150
   if (mth.eq.18.and.mtr18.gt.0) go to 150
   if (iverf.lt.6) then
      if (mth.ge.451.and.mth.lt.700) go to 150
      if (mth.ge.800) go to 150
   else
      if (mth.ge.451.and.mth.lt.600) go to 150
      if (mth.gt.850.and.mth.le.870) go to 150
      if (mth.gt.891) go to 150
   endif
   go to 180
  150 continue
   call tosend(nin,0,0,scr)
   go to 110

   !--process this reaction
  180 continue
   call contio(0,nout,0,scr,nb,nw)
   if (mfh.ne.23.and.awin.ne.zero) awrx=c2h/awin
  181 continue
   call tab1io(nin,0,0,scr,nb,nw)
   if (mfh.eq.23) go to 190
   qx=c2h
   if (mth.eq.19.and.mtr18.gt.0) q18=qx
   if (awin.ne.0) then
       thrx=-qx*(awrx+1)/awrx
   else
       thrx=-qx
   endif
   thr6=thrx
   if (thr6.lt.zero) thr6=0
   if (thrx.le.zero) go to 190
   thrxx=sigfig(thrx,7,+1)
   l=7+2*nint(scr(5))
   if (scr(l+1).ne.zero) then
      write(text,'(''xsec nonzero at threshold for mt='',i3)') mth
      call mess('lunion',text,'adjusted using jump in xsec')
   endif
   if (scr(l).ge.thrxx) go to 190
   thrx=thrxx
   write(nsyso,'(/&
     &'' changed threshold from'',1p,e13.6,'' to'',e13.6,&
     &'' for mt'',i3,''.'')') scr(l),thrx,mth
   scr(l)=thrx
   lsave=l
   do while (scr(l+2).le.scr(l))
      if (l.gt.lsave+20) call error('lunion',&
         'ill-behaved threshold.',' ')
      scr(l+2)=sigfig(scr(l),7,+1)
      l=l+2
   enddo
  190 continue
   call tab1io(0,nout,0,scr,nb,nw)

   !--set up points for first panel
   in=0
   ig=1
   ir=0
   jr=1
   lr=0
   npr=nint(scr(6))
   ibase=6+2*nint(scr(5))
   nbta=1
  205 continue
   ir=ir+1
   if (lr+ir.gt.nbta) jr=jr+1
   nbta=nint(scr(5+2*jr))
   inta=nint(scr(6+2*jr))
   er=scr(ibase+ir*2-1)
   sr=scr(ibase+ir*2)
   eg=er
   if (mth.eq.2) go to 210
   if (mth.eq.18.or.mth.eq.19) go to 210
   if (mth.eq.102) go to 210
   ! check for pseudo threshold
   if ((ibase+2*ir).lt.nw.or.nb.eq.0) go to 207
   call moreio(nin,nout,0,scr(ibase+1),nb,nw)
   lr=lr+ir
   ir=0
   nw=ibase+nw
  207 continue
   ernext=scr(ibase+ir*2+1)
   srnext=scr(ibase+ir*2+2)
   if (sr.lt.ssmall.and.srnext.lt.ssmall.and.ir.lt.npr-1) go to 205
   ! check for initial discontinuity
   if (abs(er-ernext).gt.small*er) go to 210
   er=sigfig(er,7,0)
   ernext=sigfig(ernext,7,+1)
   scr(ibase+ir*2+1)=ernext
   210 if (ig.gt.ngo) go to 220
   call finda(ig,aa,1,iold,bufo,nbuf)
   eg=aa(1)
   if (abs(eg).ge.er*(1-small)) go to 220
   enl=eg
   in=in+1
   aa(1)=enl
   call loada(in,aa,1,inew,bufn,nbuf)
   ig=ig+1
   go to 210
  220 continue
   test=elow
   test=sigfig(test,7,+1)
   if (er.lt.test) go to 225
   if ((er-sigfig(abs(eg),7,-1)).gt.1.e-8_kr*er) go to 222
   if (sr.eq.zero) go to 225
   enl=sigfig(er,7,0)
   er=sigfig(er,7,+1)
   in=in+1
   aa(1)=enl
   call loada(in,aa,1,inew,bufn,nbuf)
   go to 225
  222 continue
   er=abs(eg)
  225 continue
   enl=er
   snl=sr
   erl=er
   srl=sr
   if (ig.lt.ngo.and.abs(eg).le.er*(1+small)) then
      ig=ig+1
      call finda(ig,aa,1,iold,bufo,nbuf)
      eg=aa(1)
   endif
   ir=ir+1
   er=scr(ibase+ir*2-1)
   sr=scr(ibase+ir*2)
   if (lr+ir.gt.nbta) then
      jr=jr+1
      nbta=nint(scr(5+2*jr))
      inta=nint(scr(6+2*jr))
   endif
   ! check ahead for discontinuity
   if (ir.lt.nbta) then
      ernext=scr(ibase+2*ir+1)
      srnext=scr(ibase+2*ir+2)
      if (abs(er-ernext).le.small*er) then
         er=sigfig(er,7,-1)
         ernext=sigfig(ernext,7,+1)
         scr(ibase+2*ir+1)=ernext
      endif
   endif
   et=0

   !--determine upper limit of this panel
  240 continue
   if (lr+ir.gt.n2h) go to 350
   if (ig.ge.ngo) go to 250
   if (abs(eg).ge.er*(1-small)) go to 250
   en=abs(eg)
   call terp1(abs(erl),srl,er,sr,en,sn,inta)
   if (eg.lt.zero) en=-en
   intn=inta
   ig=ig+1
   call finda(ig,aa,1,iold,bufo,nbuf)
   eg=aa(1)
   go to 300
   ! histogram interpolation
  250 continue
   if (inta.ne.1) go to 255
   if (lr+ir.ge.n2h) go to 255
   if (er.lt.(1+small)*enl) go to 255
   if (abs(et-enl).lt.small*et) go to 252
   en=er
   et=sigfig(en,7,-1)
   en=et
   sn=snl
   intn=inta
   erl=en
   srl=sr
   go to 300
  252 continue
   er=sigfig(er,7,+1)
   ! normal interpolation
  255 continue
   en=er
   sn=sr
   intn=inta
   erl=er
   srl=sr
   260 ir=ir+1
   if (lr+ir.gt.n2h) go to 280
   er=scr(ibase+ir*2-1)
   sr=scr(ibase+ir*2)
   if (lr+ir.gt.nbta) then
      jr=jr+1
      nbta=nint(scr(5+2*jr))
      inta=nint(scr(6+2*jr))
   endif
   if ((ibase+2*ir).ge.nw.and.nb.ne.0) then
      call moreio(nin,nout,0,scr(ibase+1),nb,nw)
      lr=lr+ir
      ir=0
      nw=ibase+nw
   endif
   ! check ahead for a discontinuity.
   if (lr+ir.ge.n2h) go to 280
   ernext=scr(ibase+2*ir+1)
   srnext=scr(ibase+2*ir+2)
   if (abs(ernext-er).gt.small*er) go to 280
   if (abs(srnext-sr).lt.small*sr) go to 260
   er=sigfig(er,7,-1)
   ernext=sigfig(ernext,7,+1)
   scr(ibase+2*ir+1)=ernext
  280 continue
   if (abs(eg).gt.erl*(1+small)) go to 300
   ig=ig+1
   if (ig.gt.ngo) go to 300
   call finda(ig,aa,1,iold,bufo,nbuf)
   eg=aa(1)
  300 continue
   if (abs(en).le.abs(enl)*(1+small)) go to 240

   !--check panel for linearity and add points if needed.
   i=2
   x(2)=abs(enl)
   y(2)=snl
   x(1)=abs(en)
   y(1)=sn
  310 continue
   if (x(i-1)/x(i).lt.stpmin) go to 325
   ! force the 1,2,5,10,... points and the .0253 ev
   ! thermal point to appear in the grid.
   if (x(i).gt.elim) go to 320
   do 312 ipwr=-4,5
      xm=ten**ipwr
      if (x(i-1).gt.up*xm.and.x(i).lt.dn*xm) go to 315
      xm=xm/2
      xm=sigfig(xm,1,0)
      if (x(i-1).gt.up*xm.and.x(i).lt.dn*xm) go to 315
      xm=4*xm/10
      xm=sigfig(xm,1,0)
      if (x(i-1).gt.up*xm.and.x(i).lt.dn*xm) go to 315
  312 continue
   xm=therm
   if (x(i-1).gt.up*xm.and.x(i).lt.dn*xm) go to 315
   go to 320
   315 call terp1(x(i),y(i),x(i-1),y(i-1),xm,ym,intn)
   go to 330
   320 if (intn.gt.2) go to 322
   ! linear interpolation.
   ! for energies less than elim, add points by dividing
   ! panels in half to make sure that a 1/v shape is represented
   ! to the desired accuracy.
   if (x(i).gt.elim) go to 325
   if (x(i-1)/x(i).le.stpmax) go to 325
   xm=(x(i-1)+x(i))/2
   xm=sigfig(xm,7,0)
   call terp1(x(i),y(i),x(i-1),y(i-1),xm,ym,2)
   go to 330
   ! nonlinear interpolation.
   ! add points by dividing panels in half to make sure that
   ! the cross section is represented to the desired accuracy.
   322 xm=(x(i-1)+x(i))/2
   xm=sigfig(xm,7,0)
   call terp1(x(i),y(i),x(i-1),y(i-1),xm,ym,intn)
   call terp1(x(i),y(i),x(i-1),y(i-1),xm,yl,2)
   errn=err
   if (x(i-1).lt.trange) errn=err/5
   test=errn*ym+ssmall
   if (abs(ym-yl).lt.test) go to 325
   go to 330
  325 continue
   in=in+1
   ent=-x(i)
   if (abs(abs(enl)-x(i)).lt.small*x(i)) ent=enl
   aa(1)=ent
   call loada(in,aa,1,inew,bufn,nbuf)
   i=i-1
   if (i.gt.1) go to 310
   enl=en
   snl=sn
   go to 240
  330 continue
   i=i+1
   if (i.le.ndim) go to 340
   call error('lunion','exceeded stack.',' ')
  340 continue
   x(i)=x(i-1)
   y(i)=y(i-1)
   x(i-1)=xm
   y(i-1)=ym
   go to 310
  350 continue
   in=in+1
   inx=-in
   aa(1)=en
   call loada(inx,aa,1,inew,bufn,nbuf)

   !--finished with this reaction
   ngn=in
   ngo=ngn
   isave=iold
   iold=inew
   inew=isave
   iss=iss-1
   if (iss.gt.0) go to 181

   !--write send record and find next reaction.
   call tosend(nin,0,0,scr)
   call asend(nout,0)
   go to 110

   !--linearization and unionization complete
  410 continue
   ngr=ngo
   ngneg=0
   jg=0
   egl=0
   do 420 ig=1,ngo
   call finda(ig,aa,1,iold,bufo,nbuf)
   eg=aa(1)
   if (abs(abs(eg)-egl).le.small*egl) go to 430
   if (eg.lt.zero) ngneg=ngneg+1
   eg=abs(eg)
   if (ig.eq.1.or.ig.eq.ngo) go to 415
   if (eg.gt.emax) go to 415
   if (abs(eg-eresl).le.small*eresl) go to 430
   if (abs(eg-eresr).le.small*eresr) go to 430
   if (abs(eg-eresu).le.small*eresu) go to 430
   if (abs(eg-eresm).le.small*eresm) go to 430
   if (abs(eg-eresh).le.small*eresh) go to 430
  415 continue
   egl=eg
   jg=jg+1
   igt=jg
   if (ig.eq.ngo) igt=-igt
   aa(1)=eg
   call loada(igt,aa,1,ngrid,bufn,nbuf)
   go to 420
  430 continue
   ngr=ngr-1
  420 continue
   ngo=ngr
   ngpos=ngo-ngneg
   call timer(time)
   write(nsyso,'(/&
     &'' number of user and resonance nodes           ='',i8,/&
     &'' points in initial unionized grid             ='',i8,/&
     &'' points added by linearization                ='',i8,14x,f8.1,''s'')')&
     nodes,ngpos,ngneg,time
   deallocate(y)
   deallocate(x)
   deallocate(bufn)
   deallocate(bufo)
   return
   end subroutine lunion

   subroutine resxs(ngrid,nout,nscr,nrtot)
   !-------------------------------------------------------------------
   ! Compute resonance cross sections on existing energy grid.
   ! Reconstruct the resonance cross sections between the grid
   ! points to the specified tolerance.  The final unionized
   ! energy grid is on nout with the resonance cross sections.
   !-------------------------------------------------------------------
   use util   ! provides error,timer,releas,closz,openz
   use mainio ! provides nsyso
   ! externals
   integer::ngrid,nout,nscr,nrtot
   ! internals
   integer::ngneg,nmax,iii,in,ig,i,j,k,l,inf,ip,ndig,no
   real(kr)::cint,cmax,eint,emax,fint,fmax
   real(kr)::elo,ehi,xxx,e,egl,eg
   real(kr)::xm,dx,c1,c2
   real(kr)::fr1,fr2,sl,dm1,dm2,dm3,errm,errn
   real(kr)::tsti,est,time
   real(kr)::dm(10)
   real(kr)::sig(nsig),sigl(nsig),res(1+nsig)
   real(kr)::sigd(10,10),resl(101)
   real(kr)::aa(1)
   real(kr),dimension(:),allocatable::bufr,bufg,bufl
   real(kr),dimension(:),allocatable::x,y
   real(kr),dimension(:,:),allocatable::sigs
   integer,parameter::ndim=30
   real(kr),parameter::half=0.5e0_kr
   real(kr),parameter::estp=4.1e0_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::one=1,ten=10
   real(kr),parameter::trange=.4999e0_kr
   real(kr),parameter::zero=0
   real(kr),parameter::tenth=0.1e0_kr
   integer::nscrl=16

   !--if urr region only and lssf=1, set nrtot=0 & return
   if (lssf.eq.1.and.eresl.eq.eresu) then
      nrtot=0
      return
   endif

   !--initialize.
   call openz(-nscr,1)
   call openz(-nscrl,1)
   allocate(bufr(nbufr))
   allocate(bufg(nbufg))
   allocate(bufl(nbufl))
   allocate(x(ndim))
   allocate(y(ndim))
   allocate(sigs(nsig-1,ndim))
   ngneg=0
   nmax=0
   cmax=0
   cint=0
   eint=0
   emax=0
   fint=0
   fmax=0
   do i=1,nsig+1
      res(i)=0
   enddo
   elo=eresl
   iii=int(log10(elo*(1+small))+small)
   if (log10(elo*(1+small)).lt.zero) iii=iii-1
   xxx=ten**iii
   if (2*xxx.gt.elo*(1+small).and.elo.ge.one*(1-small)) then
      ehi=2*xxx
   else if (5*xxx.gt.elo*(1+small).and.elo.ge.one*(1-small)) then
      ehi=5*xxx
   else
      ehi=10*xxx
   endif
   if (ehi.gt.eresr) ehi=eresr
   e=0
   call csunr1(e,sig)
   call csunr2(e,sig)
   write(nsyso,'(/&
     &15x,''estimated maximum error due to''/&
     &15x,''resonance integral check (errmax,errint)''//&
     &4x,''upper'',6x,''elastic'',3x,''percent'',3x,&
     &''capture'',3x,''percent'',3x,''fission'',3x,''percent''/&
     &4x,''energy'',5x,''integral'',3x,''error'',4x,&
     &''integral'',3x,''error'',4x,''integral'',3x,''error''/&
     &1p,e10.2)') elo

   !--set up the first grid point in the resonance range.
   in=0
   ig=0
  100 continue
   ig=ig+1
   call finda(ig,aa,1,ngrid,bufg,nbufg)
   egl=aa(1)
   if (egl.lt.eresl) go to 100
   call sigma(egl,sigl,sigd)

   !--set up the next grid point
  105 continue
   ig=ig+1
   call finda(ig,aa,1,ngrid,bufg,nbufg)
   eg=aa(1)
   call sigma(eg,sig,sigd)
   if (eg.ge.eresh) go to 190

   !--add points between nodes if needed.
   i=2
   x(2)=egl
   y(2)=sigl(1)
   x(1)=eg
   y(1)=sig(1)
   do j=1,nsig-1
      sigs(j,2)=sigl(j+1)
      sigs(j,1)=sig(j+1)
   enddo
   if (x(2).lt.eresr.and.x(1).gt.eresr) go to 150
   if (x(2).lt.eresu.and.x(1).gt.eresu) go to 150
   if (x(2).lt.eresm.and.x(1).gt.eresm) go to 150
   ! compare the true function to linear interpolation at the mid point
  115 continue
   xm=half*(x(i)+x(i-1))
   dx=x(i-1)-x(i)
   ndig=9
   if (xm.gt.tenth.and.xm.lt.one) ndig=8
   if (xm.gt.sigfig(x(i),ndig,+1).and.&
     xm.lt.sigfig(x(i-1),ndig,-1)) go to 135
   ! convergence forced by significant figures check
   go to 150
  135 continue
   if (xm.gt.sigfig(x(i),7,+1).and.&
     xm.lt.sigfig(x(i-1),7,-1)) then
      xm=sigfig(xm,7,0)
   else
      xm=sigfig(xm,ndig,0)
   endif
   call sigma(xm,sig,sigd)
   fr2=(xm-x(i))/dx
   fr1=1-fr2
   do j=1,nsig-1
      ! temporary patch to cope with the shift of l>0 resonances
      ! when RML (LRF=7) is used.  the initial grid points are
      ! wrong, resulting in sometimes missing part of the peak.
      ! if (sigs(j,i-1).ne.zero) then
      !    if (sigs(j,i)/sigs(j,i-1).gt.1.25e0_kr) go to 175
      ! endif
      ! if (sigs(j,i).ne.zero) then
      !    if (sigs(j,i-1)/sigs(j,i).gt.1.25e0_kr) go to 175
      ! endif
      sl=fr1*sigs(j,i)+fr2*sigs(j,i-1)
      dm(j)=abs(sig(j+1)-sl)
   enddo
   errn=err
   if (x(i-1).lt.trange) errn=err/5
   errm=errmax
   if (x(i-1).lt.trange) errm=errmax/5
   no=0
   do j=1,nsig-1
      if (dm(j).gt.errn*sig(j+1)) no=1
   enddo
   if (no.eq.0) go to 145
   do j=1,nsig-1
      if (dm(j).gt.errm*sig(j+1)) go to 175
   enddo
   tsti=2*errint*xm/dx
   do j=1,nsig-1
      if (dm(j).ge.tsti) go to 175
   enddo
   dm1=dm(1)
   dm2=dm(2)
   dm3=dm(3)

   ! convergence forced by resonance integral check
   nmax=nmax+1
   cmax=cmax+dm3*dx/(2*xm)
   emax=emax+dm1*dx/(2*xm)
   fmax=fmax+dm2*dx/(2*xm)
   go to 150
   ! don't allow big increases in the energy step
   ! they may be misconvergences
  145 continue
   est=estp*(x(i)-res(1))
   if (in.gt.3.and.dx.gt.est.and.xm.lt.eresr) go to 175
   ! converged.  write out top point in stack.
  150 continue
   in=in+1
   res(1)=x(i)
   if (abs(egl-x(i)).le.small*egl) then
      res(1)=egl
   else
      ngneg=ngneg+1
   endif
   res(2)=y(i)
   do j=1,nsig-1
      res(j+2)=sigs(j,i)
   enddo
   call loada(in,res,nsig+1,nscr,bufr,nbufr)
   if (ncoef.gt.1) then
      call sigma(res(1),sig,sigd)
      do j=1,nmtres
         if (mcards(j).eq.0) then
            if (sigd(1,j).gt.zero) then
               mcards(j)=in
            endif
         endif
      enddo
      k=1
      resl(k)=res(1)
      do l=1,ncoef
         do j=1,nmtres
            k=k+1
            resl(k)=sigd(l,j)
         enddo
      enddo
      call loada(in,resl,1+ncoef*nmtres,nscrl,bufl,nbufl)
   endif
   c1=sigs(3,i)
   c2=sigs(3,i-1)
   cint=cint+(c1+c2)*dx/(2*xm)
   c1=sigs(1,i)
   c2=sigs(1,i-1)
   eint=eint+(c1+c2)*dx/(2*xm)
   c1=sigs(2,i)
   c2=sigs(2,i-1)
   fint=fint+(c1+c2)*dx/(2*xm)
   if (x(i-1).lt.ehi*(1-small)) go to 165
   if (ehi.gt.eresh) go to 165
   if (elo.ge.ehi*(1-small)) go to 165
   ! print resonance integral diagnostics
   if (cint.ne.zero) cmax=100*cmax/cint
   if (eint.ne.zero) emax=100*emax/eint
   if (fint.ne.zero) fmax=100*fmax/fint
   write(nsyso,'(1p,e10.2,1x,3(1p,e12.2,f8.3))')&
     x(i-1),eint,emax,cint,cmax,fint,fmax
   elo=x(i-1)
   iii=int(log10(elo*(1+small))+small)
   if (log10(elo*(1+small)).lt.0) iii=iii-1
   xxx=ten**iii
   if (2*xxx.gt.elo*(1+small).and.elo.ge.one*(1-small)) then
      ehi=2*xxx
   else if (5*xxx.gt.elo*(1+small).and.elo.ge.one*(1-small)) then
      ehi=5*xxx
   else
      ehi=10*xxx
   endif
   if (ehi.gt.eresr) ehi=eresr
   if (elo.ge.eresr) ehi=eresh
   cint=0
   eint=0
   fint=0
   cmax=0
   emax=0
   fmax=0
  165 continue
   i=i-1
   if (i.gt.1) go to 115
   ! finished with stack.  go to next panel.
   egl=eg
   sigl(1)=y(1)
   do j=1,nsig-1
      sigl(j+1)=sigs(j,1)
   enddo
   go to 105
   ! not converged.  add mid point to stack and continue.
  175 continue
   i=i+1
   if (i.gt.ndim) call error('resxs','stack exceeded.',' ')
   x(i)=x(i-1)
   y(i)=y(i-1)
   x(i-1)=xm
   y(i-1)=sig(1)
   do j=1,nsig-1
      sigs(j,i)=sigs(j,i-1)
      sigs(j,i-1)=sig(j+1)
   enddo
   go to 115
   ! write out last point in grid.
  190 continue
   in=in+1
   inf=-in
   res(1)=egl
   res(2)=sigl(1)
   do j=1,nsig-1
      res(j+2)=sigs(j,1)
   enddo
   nrtot=in
   call loada(inf,res,nsig+1,nscr,bufr,nbufr)
   if (ncoef.gt.1) then
      call sigma(res(1),sig,sigd)
      k=1
      resl(k)=res(1)
      do l=1,ncoef
         do j=1,nmtres
            k=k+1
            resl(k)=sigd(l,j)
         enddo
      enddo
      call loada(inf,resl,1+ncoef*nmtres,nscrl,bufl,nbufl)
      do i=1,nmtres
         if (mcards(i).gt.0) mcards(i)=3+2*(nrtot-mcards(i)+1)
      enddo
   endif

   !--copy resonance data to output
   ip=1
   call finda(ip,x,nsig+1,nscr,bufr,nbufr)
   in=1
   call loada(in,x,nsig+1,nout,bufg,nbufg)
   do while (ip.lt.nrtot)
      ip=ip+1
      call finda(ip,x,nsig+1,nscr,bufr,nbufr)
      in=in+1
      inf=in
      if (ip.eq.nrtot) inf=-in
      call loada(inf,x,nsig+1,nout,bufg,nbufg)
   enddo

   !--finished writing resonance tape.
   nrtot=in
   call timer(time)
   write(nsyso,'(/&
     &'' points added by resonance reconstruction     ='',i8/&
     &'' points affected by resonance integral check  ='',i8/&
     &'' final number of resonance points             ='',i8)')&
     ngneg,nmax,nrtot
   deallocate(sigs)
   deallocate(y)
   deallocate(x)
   deallocate(bufg)
   deallocate(bufr)
   call closz(-nscr)
   return
   end subroutine resxs

   subroutine sigma(e,sig,sigd)
   !-------------------------------------------------------------------
   ! Calculate cross section at energy e from resonance parameters.
   !-------------------------------------------------------------------
   use samm ! provides cssammy
   use util ! provides error
   ! externals
   real(kr)::e,sig(nsig),sigd(10,10)
   ! internals
   integer::i,j,mode
   real(kr)::abn
   real(kr)::sigp(nsig)
   real(kr)::siga(ncoef,nmtres),sigx(1,1)
   real(kr),parameter::zero=0

   !--initialize reaction cross sections to zero
   do i=1,nsig
      sig(i)=0
   enddo
   if (ncoef.gt.1) then
      do i=1,ncoef
         do j=1,nmtres
            siga(i,j)=0
         enddo
      enddo
   endif

   !--calculate cross sections from data in a for each section
   do i=1,nsect
      isect=i
      abn=abnt(i)
      ! test if present section contributes
      if (e.ge.elt(i).and.e.lt.eht(i)) then
         ! determine formalism
         mode=modet(i)
         if (mode.ne.11.and.mode.ne.12) then
            if (mode.eq.0) then
               ! no res parameters
               call csnorp(e,sigp)
            else if (mode.eq.1) then
               ! single level breit-wigner
               call csslbw(e,sigp)
            else if (mode.eq.2.and.tempr.eq.zero) then
               ! multilevel breit-wigner
               call csmlbw(e,sigp)
            else if (mode.eq.2.and.tempr.gt.zero) then
               ! multilevel breit-wigner using g/h method
               call csmlbw2(e,sigp)
            else if (mode.eq.3.and.(isammy.eq.0.or.nmtres.eq.0)) then
               ! multi-level reich-moore (r-matrix)
               call csrmat(e,sigp)
            else if (mode.eq.3.and.isammy.eq.1.and.nmtres.gt.0) then
               ! sammy calculation
               call cssammy(e,sigp,siga,sigx,mmtres,nmtres,ncoef,nresp,i)
            else if (mode.eq.4) then
               ! multilevel adler-adler
                  call csaa(e,sigp)
            else if (mode.eq.5) then
               ! general r-matrix
               call error('sigma',&
                 'general r-matrix not installed.',' ')
            else if (mode.eq.6) then
               ! hybrid r-function
               call cshybr(e,sigp)
            else if (mode.eq.7) then
               ! sammy calculation
               call cssammy(e,sigp,siga,sigx,mmtres,nmtres,ncoef,nresp,i)
            else
               ! illegal mode
               call error('sigma','illegal option.',' ')
            endif
            do j=1,nsig
               if (sigp(j).lt.zero) sigp(j)=0
               sig(j)=sig(j)+sigp(j)*abn
            enddo
         endif
      endif
   enddo

   !--output angular coefficients
   if (ncoef.gt.1) then
      do i=1,ncoef
         do j=1,nmtres
            sigd(i,j)=siga(i,j)
         enddo
      enddo
    endif

   !--add unresolved contribution, if any.
   if (e.lt.eresu) return
   if (e.ge.eresh) return
   call sigunr(e,sigp)
   do j=1,4
      sig(j)=sig(j)+sigp(j)
   enddo
   return
   end subroutine sigma

   subroutine csnorp(e,sigp)
   !-------------------------------------------------------------------
   ! Calculates potential scattering cross section for an isotope
   ! with a section without res paramters (lru=0, lrf=0)
   !-------------------------------------------------------------------
   use physics ! provides pi
   ! externals
   real(kr)::e,sigp(nsig)
   ! internals
   integer::i,inow
   real(kr)::ap

   ! zero out cross section contribution from this section (sigp)
   do i=1,nsig
      sigp(i)=0
   enddo
   ! retrieve starting location for data in a
   inow=ibaset(isect)
   ap=res(inow+7)
   sigp(2)=4*pi*ap*ap
   return
   end subroutine csnorp

   subroutine csslbw(e,sigp)
   !-------------------------------------------------------------------
   ! Calculates single level Breit-Wigner cross sections at energy e
   ! for one section (one isotope-one energy range)
   !-------------------------------------------------------------------
   use physics ! provides pi,bk,amassn,amu,hbar,ev
   use util    ! provides sigfig
   ! externals
   real(kr)::e,sigp(nsig)
   ! internals
   integer::i,ki,inow,naps,nls,l,nrs,ll,lp,in
   real(kr)::rpi,cwaven
   real(kr)::tbk,awri,ap,aw,ra,spi,spifac,k,pifac,rho,rhoc
   real(kr)::qx,se,pe,sec,pec,pex,phi,cos2p,sin2p,sinsq,spot
   real(kr)::er,aj,gn,gg,gf,ser,per,rper,gc,gx,gj
   real(kr)::erp,edelt,gne,gtt,ex,delta,theta,ax,y
   real(kr)::rew,aimw,psi,chi,smax,comfac,arat,rhop,add
   real(kr),parameter::rc1=.123e0_kr
   real(kr),parameter::rc2=.08e0_kr
   real(kr),parameter::third=.333333333e0_kr
   real(kr),parameter::zero=0
   rpi=sqrt(pi)
   cwaven=sqrt(2*amassn*amu*ev)*1.e-12_kr/hbar

   !--compute cross sections for this energy
   do i=1,nsig
      sigp(i)=0
   enddo
   tbk=tempr*bk
   ki=1
   ! retrieve starting location for data in res
   inow=ibaset(isect)
   ! retrieve nuclide information
   naps=nint(res(inow+5))
   awri=res(inow+12)
   ap=res(inow+7)
   ! calculate channel radius (ra)
   aw=amassn*awri
   ra=rc1*aw**third+rc2
   if (naps.eq.1) ra=ap
   spi=res(inow+6)
   spifac=1/(2*spi+1)
   nls=nint(res(inow+10))
   inow=inow+12
   ! calculate wave number(k),rho and rhocap at energy (e)
   arat=awri/(awri+1)
   k=cwaven*arat*sqrt(abs(e))
   pifac=pi/(k*k)
   rho=k*ra
   rhoc=k*ap

   !--loop over l states
   do l=1,nls
      nrs=nint(res(inow+5))
      ll=nint(res(inow+2))
      qx=res(inow+1)
      lrx=nint(res(inow+3))
      in=inow+nrs*6+6

      !--calculate shift and penetration factors
      !--at cross section energy
      call facts(ll,rho,se,pe)
      pec=0
      if (lrx.ne.0) then
         rhop=cwaven*arat*sqrt(abs(e+qx/arat))*ra
         lp=ll
         ! set competing l' value for a 2+ inelastic residual (u238)
         if (ll.eq.0) lp=2
         if (ll.eq.2) lp=0
         call facts(lp,rhop,sec,pec)
         if (e+qx/arat.lt.zero) pec=0
      endif
      call facphi(ll,rhoc,phi)

      !--constants independent of res. energy
      cos2p=cos(2*phi)
      sin2p=sin(2*phi)
      sinsq=(sin(phi))**2

      ! calculate potential scattering
      spot=4*(2*ll+1)*pifac*sinsq

      !--use doppler-broadened line shapes
      if (tempr.gt.zero) then
         do i=1,nrs
            inow=inow+6
            er=res(inow)
            aj=res(inow+1)
            gn=res(inow+3)
            gg=res(inow+4)
            gf=res(inow+5)
            ser=res(in)
            per=res(in+1)
            rper=1/per
            gc=res(in+2)
            pex=res(in+3)
            in=in+6
            gx=gg+gf
            gj=(2*aj+1)*spifac/2
            erp=er+gn*(ser-se)*rper/2
            edelt=e-erp
            gne=gn*pe*rper
            gtt=gne+gx
            if (gc.ne.zero) gtt=gtt+gc*pec/pex
            ex=2*(e-erp)/gtt
            delta=sqrt(4*tbk*e/awri)
            theta=gtt/delta
            ax=theta*ex/2
            y=theta/2
            call quickw(ax,y,rew,aimw,ki,tr,ti)
            psi=rpi*theta*rew/2
            chi=rpi*theta*aimw/2
            smax=4*pifac*gj*gne/gtt**2
            sigp(2)=sigp(2)+smax*((cos2p*gtt-gx)*psi+sin2p*chi*gtt)
            sigp(3)=sigp(3)+smax*gf*psi
            sigp(4)=sigp(4)+smax*gg*psi
            if (lrx.ne.0.and.gc.ne.zero) sigp(5)=sigp(5)+smax*gc*pec*psi/pex
         enddo
         inow=in

      !--use zero-temperature line shapes
      else
         do i=1,nrs
            inow=inow+6
            er=res(inow)
            aj=res(inow+1)
            gn=res(inow+3)
            gg=res(inow+4)
            gf=res(inow+5)
            ser=res(in)
            per=res(in+1)
            rper=1/per
            gc=res(in+2)
            pex=res(in+3)
            in=in+6
            gx=gg+gf
            gj=(2*aj+1)*spifac/2
            erp=er+gn*(ser-se)*rper/2
            edelt=e-erp
            gne=gn*pe*rper
            gtt=gne+gx
            if (gc.ne.zero) gtt=gtt+gc*pec/pex
            comfac=pifac*gj*gne/(edelt**2+gtt*gtt/4)
            add=comfac*(gne*cos2p-2*gx*sinsq+2*edelt*sin2p)
            sigp(2)=sigp(2)+add
            sigp(3)=sigp(3)+comfac*gf
            sigp(4)=sigp(4)+comfac*gg
            if (lrx.ne.0.and.gc.ne.zero) sigp(5)=sigp(5)+comfac*gc*pec/pex
         enddo
         inow=in
      endif

      !--add potential scattering
      sigp(2)=sigfig(sigp(2),8,0)
      spot=sigfig(spot,8,0)
      sigp(2)=sigp(2)+spot
   enddo

   !--construct the total cross section
   sigp(1)=sigp(2)+sigp(3)+sigp(4)
   if (lrx.ne.0) sigp(1)=sigp(1)+sigp(5)
   return
   end subroutine csslbw

   subroutine csmlbw(e,sigp)
   !-------------------------------------------------------------------
   ! Calculates multilevel Breit-Wigner cross sections at energy e
   ! for one section (one isotope-one energy range)
   !-------------------------------------------------------------------
   use physics ! provides pi,amassn,amu,hbar,ev
   use util    ! provides error
   use endf    ! provides terpa
   ! externals
   real(kr)::e,sigp(nsig)
   ! internals
   integer::i,inow,nro,naps,iro,ip,ir,idx,nls,l,ll,lp,nrs
   integer::nj,ii,in,j
   real(kr)::cwaven,ape,enx,awri,ap,aw,ra,spi,den
   real(kr)::arat,k,pifac,rho,rhoc,qx,se,pe,rhop,sec,pec,pex
   real(kr)::phi,cos2p,sin2p,sum,fl,ajmin,ajmax,aj,diff
   real(kr)::er,gn,gg,gf,ser,per,gc,erp,edelt,gne,gx,gtt
   real(kr)::x,comfac,rper,add
   real(kr)::sigj(10,2),gj(10)
   real(kr),parameter::rc1=.123e0_kr
   real(kr),parameter::rc2=.08e0_kr
   real(kr),parameter::third=.333333333e0_kr
   real(kr),parameter::half=.5e0_kr
   real(kr),parameter::zero=0
   cwaven=sqrt(2*amassn*amu*ev)*1.e-12_kr/hbar

   !--doppler broadening not allowed
   if (tempr.gt.zero) call error('csmlbw',&
     'not coded for temperature gt 0 deg k.',' ')

   !--compute cross sections for this energy
   do i=1,nsig
      sigp(i)=0
   enddo
   ! retrieve starting location for data in a
   inow=ibaset(isect)
   ! retrieve nuclide information
   nro=nint(res(inow+4))
   naps=nint(res(inow+5))
   if (nro.ne.0) then
      iro=inow+6
      inow=inow+6+2*nint(res(iro+4))+2*nint(res(iro+5))
      ip=2
      ir=1
      call terpa(ape,e,enx,idx,res(iro),ip,ir)
   endif
   awri=res(inow+12)
   ap=res(inow+7)
   aw=amassn*awri
   ra=rc1*aw**third+rc2
   if (nro.eq.0) then
      if (naps.eq.1) ra=ap
   else
      if (naps.eq.0) then
         ap=ape
      else if (naps.eq.1) then
         ap=ape
         ra=ape
      else
         ra=ap
         ap=ape
      endif
   endif
   spi=res(inow+6)
   den=4*spi+2
   nls=nint(res(inow+10))
   ! calculate wave number(k),rho and rhocap at energy (e)
   arat=awri/(awri+1)
   k=cwaven*arat*sqrt(abs(e))
   pifac=pi/(k*k)
   rho=k*ra
   rhoc=k*ap
   inow=inow+12

   !--loop over l states
   do l=1,nls
      nrs=nint(res(inow+5))
      ll=nint(res(inow+2))
      qx=res(inow+1)
      lrx=nint(res(inow+3))
      call facts(ll,rho,se,pe)
      pec=0
      if (lrx.ne.0) then
         rhop=cwaven*arat*sqrt(abs(e+qx/arat))*ra
         lp=ll
         ! set competing l' value for a 2+ inelastic residual (u238)
         if (ll.eq.0) lp=2
         if (ll.eq.2) lp=0
         call facts(lp,rhop,sec,pec)
         if (e+qx/arat.lt.zero) pec=0
      endif
      call facphi(ll,rhoc,phi)
      cos2p=1-cos(2*phi)
      sin2p=sin(2*phi)
      sum=0
      fl=ll
      ajmin=abs(abs(spi-fl)-half)
      ajmax=spi+fl+half
      nj=nint(ajmax-ajmin+1)
      aj=ajmin
      do i=1,nj
         gj(i)=(2*aj+1)/den
         aj=aj+1
         sum=sum+gj(i)
      enddo
      diff=2*fl+1-sum
      do ii=1,2
         do i=1,nj
            sigj(i,ii)=0
         enddo
      enddo
      inow=inow+6
      in=inow+nrs*6

      !--loop over all resonances
      do i=1,nrs
         er=res(inow)
         j=nint(res(inow+1)-ajmin+1)
         gn=res(inow+3)
         gg=res(inow+4)
         gf=res(inow+5)
         ser=res(in)
         per=res(in+1)
         rper=1/per
         gc=res(in+2)
         pex=res(in+3)
         in=in+6
         inow=inow+6
         erp=er+gn*(ser-se)*rper/2
         edelt=e-erp
         gne=gn*pe*rper
         gx=gg+gf
         gtt=gne+gx
         if (gc.ne.zero) gtt=gtt+gc*pec/pex
         x=2*edelt/gtt
         comfac=2*gne/gtt/(1+x*x)
         sigj(j,1)=sigj(j,1)+comfac
         sigj(j,2)=sigj(j,2)+comfac*x
         comfac=comfac*gj(j)/gtt
         sigp(3)=sigp(3)+comfac*gf
         sigp(4)=sigp(4)+comfac*gg
         if (lrx.ne.0.and.gc.ne.zero) sigp(5)=sigp(5)+comfac*gc*pec/pex
      enddo
      do j=1,nj
         add=gj(j)*((cos2p-sigj(j,1))**2+(sin2p+sigj(j,2))**2)
         sigp(2)=sigp(2)+add
      enddo
      sigp(2)=sigp(2)+2*diff*cos2p
      inow=in
   enddo

   !--construct the final cross sections
   sigp(2)=sigp(2)*pifac
   sigp(3)=sigp(3)*2*pifac
   sigp(4)=sigp(4)*2*pifac
   sigp(1)=sigp(2)+sigp(3)+sigp(4)
   if (lrx.ne.0) sigp(1)=sigp(1)+sigp(5)
   return
   end subroutine csmlbw

   subroutine csmlbw2(e,sigp)
   !-------------------------------------------------------------------
   ! Calculates multi level Breit-Wigner cross sections at energy e
   ! using the g/h method
   ! for one section (one isotope-one energy range)
   !-------------------------------------------------------------------
   use physics ! provides pi,bk,amassn,amu,hbar,ev
   use util    ! provides sigfig
   ! externals
   real(kr)::e,sigp(nsig)
   ! internals
   integer::i,j,ki,inow,naps,nls,l,nrs,ll,in,igh
   real(kr)::rpi,cwaven
   real(kr)::tbk,awri,ap,aw,ra,spi,spifac,k,pifac,rho,rhoc
   real(kr)::qx,se,pe,sec,pec,phi,cos2p,sin2p,sinsq,spot
   real(kr)::er,aj,gt,gn,gg,gf,ser,per,rper,gc,gx,gj,gne
   real(kr)::g,h,ajj,erj,gnj,ggj,gfj,gtj,sej,pej,gcj,sml,ssl
   real(kr)::erp,gtt,ex,delta,theta,ax,y,pex,rpex,pcj,enow,test
   real(kr)::rew,aimw,psi,chi,smax,comfac,arat,rhop,add
   real(kr),parameter::rc1=.123e0_kr
   real(kr),parameter::rc2=.08e0_kr
   real(kr),parameter::third=.333333333e0_kr
   real(kr),parameter::tmin=10e0_kr
   real(kr),parameter::zero=0
   save enow
   rpi=sqrt(pi)
   cwaven=sqrt(2*amassn*amu*ev)*1.e-12_kr/hbar

   !--compute cross sections for this energy
   do i=1,nsig
      sigp(i)=0
   enddo
   tbk=tempr*bk
   ki=1
   ! retrieve starting location for data in res
   inow=ibaset(isect)
   ! retrieve nuclide information
   naps=nint(res(inow+5))
   awri=res(inow+12)
   ap=res(inow+7)
   ! calculate channel radius (ra)
   aw=amassn*awri
   ra=rc1*aw**third+rc2
   if (naps.eq.1) ra=ap
   spi=res(inow+6)
   spifac=1/(2*spi+1)
   nls=nint(res(inow+10))
   inow=inow+12
   ! calculate wave number(k),rho and rhocap at energy (e)
   arat=awri/(awri+1)
   k=cwaven*arat*sqrt(abs(e))
   pifac=pi/(k*k)
   rho=k*ra
   rhoc=k*ap

   !--accelerating the gh calculation
   igh=0
   if (abs(e-1.e-5_kr).lt.1.e-7_kr) then
      enow=e
      igh=1
   endif
   test=enow+enow/20
   if (enow+10.lt.test) test=enow+10
   if (e.ge.test) then
      enow=e
      igh=1
   endif

   !--loop over l states
   do l=1,nls
      nrs=nint(res(inow+5))
      ll=nint(res(inow+2))
      qx=res(inow+1)
      lrx=nint(res(inow+3))
      in=inow+nrs*6+6

      !--calculate shift and penetration factors
      !--at cross section energy
      call facts(ll,rho,se,pe)
      pec=0
      if (lrx.ne.0) then
         rhop=cwaven*arat*sqrt(abs(e+qx))*ra
         call facts(ll,rhop,sec,pec)
         if (e+qx/arat.lt.zero) pec=0
      endif
      call facphi(ll,rhoc,phi)

      !--constants independent of res. energy
      cos2p=cos(2*phi)
      sin2p=sin(2*phi)
      sinsq=(sin(phi))**2

      ! calculate potential scattering
      spot=4*(2*ll+1)*pifac*sinsq
      sml=0

      !--loop over resonances for this sequence
      inow=inow+6
      do i=1,nrs

         !--get parameters for this resonance
         er=res(inow+6*(i-1))
         aj=res(inow+6*(i-1)+1)
         gt=res(inow+6*(i-1)+2)
         gn=res(inow+6*(i-1)+3)
         gg=res(inow+6*(i-1)+4)
         gf=res(inow+6*(i-1)+5)
         ser=res(in+6*(i-1))
         per=res(in+6*(i-1)+1)
         rper=1/per
         gc=res(in+6*(i-1)+2)
         pex=res(in+6*(i-1)+3)
         rpex=1/pex
         gx=gg+gf
         gj=(2*aj+1)*spifac/2
         erp=er+gn*(ser-se)*rper/2
         gne=gn*pe*rper
         gtt=gne+gx
         if (gc.ne.zero) gtt=gtt+gc*pec*rpex

         !--compute g and h for multi level correction
         !--at the resonance centers
         if (igh.eq.1) then
            g=0
            h=0
            do j=1,nrs
               ajj=res(inow+6*(j-1)+1)
               if (j.ne.i.and.ajj.eq.aj.and.abs(i-j).lt.200) then
                  erj=res(inow+6*(j-1))
                  sej=res(in+6*(j-1))
                  pej=res(in+6*(j-1)+1)
                  gcj=res(in+6*(j-1)+2)
                  pcj=res(in+6*(j-1)+3)
                  gnj=res(inow+6*(j-1)+3)
                  ggj=res(inow+6*(j-1)+4)
                  gfj=res(inow+6*(j-1)+5)
                  erj=erj+(sej-se)*gnj/(2*pej)
                  gnj=gnj*pe/pej
                  gtj=gnj+ggj+gfj
                  if (gcj.ne.zero) gtj=gtj+gcj*pec/pcj
                  g=g+gne*gnj*(gtt+gtj)/((erp-erj)**2+(gtt+gtj)**2/4)
                  h=h+gne*gnj*(erp-erj)/((erp-erj)**2+(gtt+gtj)**2/4)
               endif
            enddo
            g=g/2
            res(in+6*(i-1)+4)=g/e
            res(in+6*(i-1)+5)=h/e
         endif

         !--finish cross section calculation
         ex=2*(e-erp)/gtt
         if (tempr.gt.tmin) then
            delta=sqrt(4*tbk*e/awri)
            theta=gtt/delta
            ax=theta*ex/2
            y=theta/2
            call quickw(ax,y,rew,aimw,ki,tr,ti)
            psi=rpi*theta*rew/2
            chi=rpi*theta*aimw/2
         else
            psi=1/(1+ex**2)
            chi=ex/(1+ex**2)
         endif
         g=res(in+6*(i-1)+4)
         h=res(in+6*(i-1)+5)
         smax=4*pifac*gj*gne/gtt**2
         sigp(2)=sigp(2)+smax*((cos2p*gtt-gx)*psi+sin2p*gtt*chi)
         sml=sml+smax*e*(g*psi+h*chi)*gtt/gne
         sigp(3)=sigp(3)+smax*gf*psi
         sigp(4)=sigp(4)+smax*gg*psi
      enddo
      inow=inow+12*nrs

      !--add potential scattering
      ssl=sigp(2)+spot
      sigp(2)=ssl+sml
   enddo

   !--construct the total cross section
   sigp(1)=sigp(2)+sigp(3)+sigp(4)
   return
   end subroutine csmlbw2

   subroutine csrmat(e,sigp)
   !-------------------------------------------------------------------
   ! Calculates R-matrix (Reich-Moore) cross sections at energy e
   ! for one section  (one isotope-one energy range)
   !-------------------------------------------------------------------
   use physics ! provides pi,amassn,amu,hbar,ev
   use util    ! provides error
   ! externals
   real(kr)::e,sigp(nsig)
   ! internals
   integer::i,inow,naps,nls,l,inowb,nrs,ncyc,ll,numj
   integer::jj,jjl,j,in,kchanl,idone,kpstv,kngtv,iskip,kkkkkk
   real(kr)::cwaven,awri,ap,aw,ra,spi,arat,k,pifac,gjd
   real(kr)::rho,rhoc,gfa,gfb,gf,apl,se,pe,dum1
   real(kr)::phi,phid,p1,p2,fl,ajmin,ajmax,ajc,gj,aj
   real(kr)::er,gn,gg,per,a1,a2,a3,diff,den,de2,gg4
   real(kr)::t1,t2,t3,t4,u11r,u11i,termt,termn,termf,termg
   real(kr)::dd,rr,ss,amag,rri,ssi,uur,uui,xx
   real(kr)::r(3,3),s(3,3),ri(3,3),si(3,3)
   real(kr),parameter::rc1=.123e0_kr
   real(kr),parameter::rc2=.08e0_kr
   real(kr),parameter::third=.333333333e0_kr
   real(kr),parameter::quar=.25e0_kr
   real(kr),parameter::haf=.5e0_kr
   real(kr),parameter::uno=1.e0_kr
   real(kr),parameter::two=2.e0_kr
   real(kr),parameter::four=4.0e0_kr
   real(kr),parameter::small=3.e-4_kr
   real(kr),parameter::zero=0
   cwaven=sqrt(2*amassn*amu*ev)*1.e-12_kr/hbar

   !--doppler broadening not provided.
   if (tempr.gt.zero) call error('csrmat',&
     'not coded for temperatures gt 0 deg k.',' ')

   !--compute cross sections at this energy
   do i=1,nsig
      sigp(i)=0
   enddo
   ! retrieve starting location for data in a
   inow=ibaset(isect)
   ! retrieve nuclide information
   naps=nint(res(inow+5))
   awri=res(inow+12)
   ap=res(inow+7)
   aw=amassn*awri
   ra=rc1*aw**third+rc2
   if (naps.eq.1) ra=ap
   spi=res(inow+6)
   gjd=2*(2*spi+1)
   nls=nint(res(inow+10))
   ! calculate wave number(k),rho and rhocap at energy (e)
   arat=awri/(awri+1)
   k=cwaven*arat*sqrt(abs(e))
   pifac=pi/(k*k)
   rho=k*ra
   rhoc=k*ap
   gfa=0
   gfb=0
   gf=0
   inow=inow+12

   !--loop over l states
   do l=1,nls
      inowb=inow
      nrs=nint(res(inow+5))
      ncyc=nint(res(inow+4))/nrs
      ll=nint(res(inow+2))
      apl=res(inow+1)
      rhoc=k*ap
      rho=k*ra
      if (apl.ne.zero) rhoc=k*apl
      if (apl.ne.zero.and.naps.eq.1) rho=k*apl
      ! calculate shift and penetration factors
      call facts(ll,rho,se,dum1)
      pe=dum1
      call facphi(ll,rhoc,phi)
      ! constants independent of res. energy
      phid=phi
      p1=cos(2*phid)
      p2=sin(2*phid)

      !--loop over possible j values
      fl=ll
      ajmin=abs(abs(spi-fl)-haf)
      ajmax=spi+fl+haf
      numj=nint(ajmax-ajmin+1)
      ajc=ajmin-1
      if (ll.ne.0.and.(fl.gt.spi-haf.and.fl.le.spi)) then
         jjl=0
      else
         jjl=1
      endif
      do jj=1,numj
         inow=inowb
         ajc=ajc+1
         gj=(2*ajc+1)/gjd

         !--loop over possible channel spins
         kchanl=0
         idone=0
         do while (kchanl.lt.2.and.idone.eq.0)
            kchanl=kchanl+1
            inow=inowb
            kpstv=0
            kngtv=0
            !initialize matrix
               do j=1,3
                  do i=1,3
                     s(j,i)=0
                     r(j,i)=0
                  enddo
               enddo
            !--loop over resonances
            inow=inow+6
            in=inow+nrs*6
            do i=1,nrs
               aj=abs(res(inow+1))
               !select only resonances with current j value
               if (abs(aj-ajc).le.quar) then
                  if (res(inow+1).lt.zero) kngtv=kngtv+1
                  if (res(inow+1).gt.zero) kpstv=kpstv+1
                  iskip=0
                  if (kchanl.eq.1.and.res(inow+1).lt.zero) iskip=1
                  if (kchanl.eq.2.and.res(inow+1).gt.zero) iskip=1
                  if (iskip.eq.0) then
                     !retrieve parameters
                     er=res(inow)
                     gn=res(inow+2)
                     gg=res(inow+3)
                     gfa=res(inow+4)
                     gfb=res(inow+5)
                     per=res(in+1)
                     !gc=res(in+2)
                     a1=sqrt(gn*pe/per)
                     a2=0
                     if (gfa.ne.zero) a2=sqrt(abs(gfa))
                     if (gfa.lt.zero) a2=-a2
                     a3=0
                     if (gfb.ne.zero) a3=sqrt(abs(gfb))
                     if (gfb.lt.zero) a3=-a3
                     !compute energy factors
                     diff=er-e
                     den=diff*diff+quar*gg*gg
       if(den.eq.0)write(*,'("er,e,diff,gg=",4(1pe11.4))')&
                              er,e,diff,gg
                     de2=haf*diff/den
                     gg4=quar*gg/den
                     !calculate r-function, or
                     !calculate upper triangular matrix terms
                     r(1,1)=r(1,1)+gg4*a1*a1
                     s(1,1)=s(1,1)-de2*a1*a1
                     if (gfa.ne.zero.or.gfb.ne.zero) then
                        r(1,2)=r(1,2)+gg4*a1*a2
                        s(1,2)=s(1,2)-de2*a1*a2
                        r(1,3)=r(1,3)+gg4*a1*a3
                        s(1,3)=s(1,3)-de2*a1*a3
                        r(2,2)=r(2,2)+gg4*a2*a2
                        s(2,2)=s(2,2)-de2*a2*a2
                        r(3,3)=r(3,3)+gg4*a3*a3
                        s(3,3)=s(3,3)-de2*a3*a3
                        r(2,3)=r(2,3)+gg4*a2*a3
                        s(2,3)=s(2,3)-de2*a2*a3
                        gf=1
                     endif
                  endif
               endif
               inow=inow+ncyc
               in=in+6
            enddo

            !--take care of extra channel spin as defined
            !--by the sign of aj:
            !-- kkkkkk = 0 => do not add anything in here
            !-- kkkkkk = 1 => add resonance contribution but
            !--  not extra hard-sphere
            !-- kkkkkk = 2 => add resonance plus hard-sphere
            !--  phase shift contribution
            kkkkkk = 0
            if (kchanl.eq.1) then
               if (kpstv.gt.0) then
                  if (kngtv.eq.0) then
                     if (jj.gt.jjl.and.jj.lt.numj) then
                        kkkkkk=2
                     else
                        kkkkkk=1
                     endif
                  else if (kngtv.gt.0) then
                     kkkkkk=1
                  endif
               else if (kpstv.eq.0) then
                  if (kngtv.eq.0) then
                     if (jj.gt.jjl.and.jj.lt.numj) then
                        kkkkkk=2
                     else
                        kkkkkk=1
                     endif
                  else if (kngtv.gt.0) then
                     kkkkkk=0
                  endif
               endif
            else if (kchanl.eq.2) then
               if (kpstv.gt.0) then
                  if (kngtv.eq.0) then
                  else if (kngtv.gt.0) then
                     kkkkkk=1
                  endif
               else if (kpstv.eq.0) then
                  if (kngtv.eq.0) then
                  else if (kngtv.gt.0) then
                     if (jj.gt.jjl.and.jj.lt.numj) then
                        kkkkkk=2
                     else
                        kkkkkk=1
                     endif
                  endif
               endif
            endif
            if (kkkkkk.ne.0) then

               !--r-matrix path -- make symmetric matrix
               if (gf.ne.zero) then
                  r(1,1)=uno+r(1,1)
                  r(2,2)=uno+r(2,2)
                  r(3,3)=uno+r(3,3)
                  r(2,1)=r(1,2)
                  s(2,1)=s(1,2)
                  r(3,1)=r(1,3)
                  s(3,1)=s(1,3)
                  r(3,2)=r(2,3)
                  s(3,2)=s(2,3)
                  !invert the complex matrix
                  call frobns(r,s,ri,si)
                  !fission term for r-matrix path
                  t1=ri(1,2)
                  t2=si(1,2)
                  t3=ri(1,3)
                  t4=si(1,3)
                  termf=four*gj*(t1*t1+t2*t2+t3*t3+t4*t4)
                  u11r=p1*(two*ri(1,1)-uno)+two*p2*si(1,1)
                  u11i=p2*(uno-two*ri(1,1))+two*p1*si(1,1)
                  termt=two*gj*(uno-u11r)
                  termn=gj*((uno-u11r)**2+u11i**2)

               !--r-function path
               else
                  dd=r(1,1)
                  rr=uno+dd
                  ss=s(1,1)
                  amag=rr**2+ss**2
                  rri=rr/amag
                  ssi=-ss/amag
                  uur=p1*(two*rri-uno)+two*p2*ssi
                  uui=p2*(uno-two*rri)+two*p1*ssi
                  if (abs(dd).lt.small.and.&
                    abs(phid).lt.small) then
                     xx=2*dd
                     xx=xx+2*(dd*dd+ss*ss+phid*phid+p2*ss)
                     xx=xx-2*phid*phid*(dd*dd+ss*ss)
                     xx=xx/amag
                     termt=two*gj*xx
                     termn=gj*(xx**2+uui**2)
                  else
                     termt=two*gj*(uno-uur)
                     termn=gj*((uno-uur)**2+uui**2)
                  endif
                  termf=0
               endif

               !--cross sections contributions
               if (kkkkkk.eq.2) then
                  termn=termn+two*gj*(1-p1)
                  termt=termt+two*gj*(1-p1)
               endif
               termg=termt-termf-termn
               sigp(2)=sigp(2)+termn
               sigp(4)=sigp(4)+termg
               sigp(3)=sigp(3)+termf
               sigp(1)=sigp(1)+termt
            endif
         enddo
      enddo
      inow=in

   !--continue the loop over l values
   enddo

   !--calculate final cross sections and store for return
   sigp(1)=pifac*sigp(1)
   sigp(2)=pifac*sigp(2)
   sigp(3)=pifac*sigp(3)
   sigp(4)=pifac*sigp(4)
   return
   end subroutine csrmat

   subroutine frobns(a,b,c,d)
   !-------------------------------------------------------------------
   ! This subroutine inverts a complex matrix with real and imaginary
   ! parts a and b and gives c and d the real and imaginary parts of
   ! the inverse. Frobenius-Schur method of inversion.
   !-------------------------------------------------------------------
   ! externals
   real(kr)::a(3,3),b(3,3),c(3,3),d(3,3),q(3,3)
   ! internals
   integer::i,j,ind

   do i=1,3
      do j=1,3
         c(i,j)=a(i,j)
      enddo
   enddo
   call thrinv(a,3,ind)
   if (ind.ne.1) then
      call abcmat(a,b,q)
      call abcmat(b,q,d)
      do i=1,3
         do j=1,3
            c(i,j)=c(i,j)+d(i,j)
         enddo
      enddo
      call thrinv(c,3,ind)
      call abcmat(q,c,d)
      do i=1,3
         do j=1,3
            d(i,j)=-d(i,j)
         enddo
      enddo
   endif
   return
   end subroutine frobns

   subroutine thrinv(d,n,kimerr)
   !-------------------------------------------------------------------
   ! Inverts symmetric matrix (d(i,j),j=1,n,i=1,j)
   !-------------------------------------------------------------------
   ! externals
   integer::n,kimerr
   real(kr)::d(3,3)
   ! internals
   integer::i,j,lr
   real(kr)::fooey
   real(kr)::s(3)
   real(kr),parameter::zero=0
   real(kr),parameter::uno=1

   kimerr=0
   do j=1,n
      do i=1,j
         d(i,j)=-d(i,j)
         d(j,i)=d(i,j)
      enddo
      d(j,j)=uno+d(j,j)
   enddo
   do lr=1,n
      fooey=uno-d(lr,lr)
      if (fooey.eq.zero) then
         kimerr=1
         return
      endif
      d(lr,lr)=uno/fooey
      do j=1,n
         s(j)=d(lr,j)
         if (j.ne.lr) then
            d(j,lr)=d(j,lr)*d(lr,lr)
            d(lr,j)=d(j,lr)
         endif
      enddo
      do j=1,n
         if (j.ne.lr) then
            do i=1,j
               if (i.ne.lr) then
                  d(i,j)=d(i,j)+d(i,lr)*s(j)
                     d(j,i)=d(i,j)
               endif
            enddo
         endif
      enddo
   enddo
   return
   end subroutine thrinv

   subroutine abcmat(a,b,c)
   !-------------------------------------------------------------------
   ! Routine to do a matrix multiplication
   !-------------------------------------------------------------------
   ! externals
   real(kr)::a(3,3),b(3,3),c(3,3)
   ! internals
   integer::i,j,k

   do i=1,3
      do j=1,3
         c(i,j)=0
         do k=1,3
            c(i,j)=c(i,j)+a(i,k)*b(k,j)
         enddo
      enddo
   enddo
   return
   end subroutine abcmat

   subroutine cshybr(e,sigp)
   !-------------------------------------------------------------------
   ! Calculates hybrid R-function cross sections at energy e
   ! for one section (one isotope-one energy range)
   !-------------------------------------------------------------------
   use physics ! provides pi,amassn,amu,hbar,ev
   use util    ! provides error
   use endf    ! provides terpa
   ! externals
   real(kr)::e,sigp(nsig)
   ! internals
   integer::i,inow,nls,ngre,nfre,nire,ncre,nrchan,np,nr,itpt
   integer::n,l,nw,nss,ns,njs,nj,ll,lbk,lps,nrs,ii,ip,ir,idis,iii
   real(kr)::cwaven,spi,awri,arat,k,kfac,pifac
   real(kr)::aj,gj,ac,rho,phir,se,pe
   real(kr)::er,gn,gg,gf,gelim,eeff,ereff,phres,pd,phen
   real(kr)::sn,pn,eenn,rhres,ser,per,gne,gnr,gtt
   real(kr)::edelt,den,comfac,sd,rr0,ri0,phro,phri,rpart,sabs
   complex(kr)::one,eye,rfcnn,slsjnn,phi0,phi,r0,exphi
   integer::lc(4),icrpt(4,4)
   real(kr)::qre(4),gc(4),sigc(4)
   real(kr)::zero=0
   real(kr)::uno=1

   eye=dcmplx(zero,uno)
   one=dcmplx(uno,zero)
   cwaven=sqrt(2*amassn*amu*ev)*1.e-12_kr/hbar

   !--no doppler broadening
   if (tempr.gt.zero) call error('cshybr',&
     'doppler broadening not provided for hybrid r-mtrix.',' ')

   !--zero out cross section contribution from this section (sigp)
   do i=1,4
      sigp(i)=0
   enddo

   !--retrieve nuclide information
   inow=ibaset(isect)
   spi=res(inow+6)
   nls=nint(res(inow+10))

   !--retrieve competing channel descriptions
   ngre=nint(res(inow+14))
   nfre=nint(res(inow+15))
   nire=nint(res(inow+16))
   ncre=nint(res(inow+17))
   nrchan=nire+ncre
   qre(1)=res(inow+29)
   qre(2)=res(inow+30)
   qre(3)=res(inow+31)
   qre(4)=res(inow+32)
   if (nrchan.gt.0) then
      do n=1,nrchan
         sigc(n)=0
      enddo
   endif

   !--make pointer index for charged particle penetrabilities
   inow=inow+34
   if (ncre.ne.0) then
      do n=1,ncre
         do l=1,4
            icrpt(l,n)=inow
            nw=6+2*nint(res(inow+4))+2*nint(res(inow+5))
            inow=inow+nw
         enddo
      enddo
   endif

   !--calculate wave number(k) at energy(e)
   awri=res(inow)
   arat=awri/(awri+1)
   kfac=cwaven*arat
   k=kfac*sqrt(abs(e))
   pifac=pi/(k*k)

   !--loop over l states
   do 140 l=1,nls
   ll=nint(res(inow+2))
   nss=nint(res(inow+4))
   inow=inow+6

   !--loop over channel spin states
   do 130 ns=1,nss
   njs=nint(res(inow+4))
   inow=inow+6

   !--loop over j states
   do 120 nj=1,njs
   aj=res(inow)
   gj=(2*aj+1)/(4*spi+2)
   ac=res(inow+1)
   lbk=nint(res(inow+2))
   lps=nint(res(inow+3))
   nrs=nint(res(inow+5))
   inow=inow+6
   rho=k*ac
   ! calculate hard sphere phase shift
   call facphi(ll,rho,phir)
   phi0=dcmplx(phir,zero)
   ! calculate shift and penetration factors at cross section energy
   call facts(ll,rho,se,pe)

   !--loop over all resonances
   rfcnn=dcmplx(zero,zero)
   if (nrs.eq.0) go to 140
   do 110 i=1,nrs
   er=res(inow)
   gn=res(inow+1)
   gg=res(inow+2)
   gf=res(inow+3)
   inow=inow+3
   gelim=gg+gf
   if (nrchan.ne.0) then
      ! determine value of competitive widths
      do ii=1,nrchan
         gc(ii)=0
         lc(ii)=nint(res(inow+ii+4))
         eeff=e+(qre(ii)/arat)
         if (eeff.gt.zero) then
            ereff=abs(er+(qre(ii)/arat))
            ! inelastic penetrability
            if (ii.le.nire) then
               phres=kfac*sqrt(ereff)*ac
               call facts(lc(ii),phres,sd,pd)
               phen=kfac*sqrt(eeff)*ac
               call facts(lc(ii),phen,sn,pn)
               gc(ii)=pn*res(inow+ii)/pd
            ! charged particle penetrability
            else if (lc(ii)+1.le.4) then
               itpt=icrpt(lc(ii)+1,ii)
               ip=2
               ir=1
               call terpa(pd,ereff,eenn,idis,res(itpt),ip,ir)
               call terpa(pn,eeff,eenn,idis,res(itpt),ip,ir)
               gc(ii)=pn*res(inow+ii)/pd
            endif
         endif
         gelim=gelim+gc(ii)
      enddo
   endif
   inow=inow+9

   !--cross section calculations
   rhres=kfac*sqrt(abs(er))*ac
   call facts(ll,rhres,ser,per)
   gne=gn*pe/per
   gnr=gn/(2*per)
   gtt=gne+gelim
   edelt=er-e
   den=edelt**2+gtt*gtt/4
   comfac=pifac*gj*gne/den
   rfcnn=rfcnn+gnr/(edelt-eye*gelim/2)
   if (ngre.ne.0) sigp(4)=sigp(4)+comfac*gg
   if (nfre.ne.0) sigp(3)=sigp(3)+comfac*gf
   if (nrchan.gt.0) then
      do iii=1,nrchan
         sigc(iii)=sigc(iii)+comfac*gc(iii)
      enddo
   endif

   !--end of loop over resonances
  110 continue

   !--background r-function
   r0=dcmplx(zero,zero)
   if (lbk.ne.0)  then
      ip=2
      ir=1
      call terpa(rr0,eeff,eenn,idis,res(inow),ip,ir)
      np=nint(res(inow+5))
      nr=nint(res(inow+4))
      inow=inow+6+2*nr+2*np
      call terpa(ri0,eeff,eenn,idis,res(inow),ip,ir)
      np=nint(res(inow+5))
      nr=nint(res(inow+4))
      inow=inow+6+2*nr+2*np
      r0=dcmplx(rr0,ri0)
   endif
   rfcnn=rfcnn+r0

   !--optical model phase shift
   phi=phi0
   if (lps.ne.0)  then
      ip=2
      ir=1
      call terpa(phro,eeff,eenn,idis,res(inow),ip,ir)
      np=nint(res(inow+5))
      nr=nint(res(inow+4))
      inow=inow+6+2*nr+2*np
      call terpa(phri,eeff,eenn,idis,res(inow),ip,ir)
      np=nint(res(inow+5))
      nr=nint(res(inow+4))
      inow=inow+6+2*nr+2*np
      phi=dcmplx(phro,phri)
   endif

   !--r-function calculations
   exphi=exp(-2*eye*phi)
   slsjnn=(one-exphi)-eye*pe*rfcnn*(one+exphi)
   slsjnn=slsjnn/(one-eye*pe*rfcnn)
   rpart=dble(slsjnn)
   sabs=abs(slsjnn)
   sigp(1)=sigp(1)+2*pifac*gj*rpart
   sigp(2)=sigp(2)+pifac*gj*sabs*sabs
   sigp(1)=sigp(2)+sigp(3)+sigp(4)

   !--continue l, s, and j loops
  120 continue
  130 continue
  140 continue

   !--finished.
   return
   end subroutine cshybr

   subroutine csunr1(e,sigp)
   !-------------------------------------------------------------------
   ! Unresolved resonance region (format 1).
   ! Single-level Breit-Wigner formalism.
   ! Energy independent parameters:
   !    lfw=0  no fission widths given
   !    lfw=1  fission widths tabulated
   ! For lfw=0, parameter interpolation is always used.
   ! For lfw=1, linear cross section interpolation is used
   ! unless the energy step is too large (see wide), then
   ! parameter interpolation is used.  Call csunr1 with
   ! e=0 to intialize flags for each section.
   !-------------------------------------------------------------------
   use physics ! provides pi,amassn,amu,hbar,ev
   use endf    ! provides terp1
   use util    ! provides error
   ! externals
   real(kr)::e,sigp(4)
   ! internals
   integer::ne,ie,i,iest,inow,lfw,nls,l,ll,j,muf,nep6,ifst
   integer::mu,nu,lamda,i1,i2,njs,nro,ip,ir,iro,idx,naps
   real(kr)::ee,enext,elast,spi,ay,aaa,awri,rat,aw,aa,const
   real(kr)::gfx,gxx,dx,aj,amun,gnox,ggx,gj
   real(kr)::e2,k,rho,rhoc,vl,ps,gnx,diff,den,enx
   real(kr)::cwaven,spot,y1,y2
   real(kr)::temp,terg,ters,terf,gs,gc,gff,add
   integer::ifirst=1
   real(kr),parameter::wide=1.26e0_kr
   real(kr),parameter::rc1=.123e0_kr
   real(kr),parameter::rc2=.08e0_kr
   real(kr),parameter::third=.333333333e0_kr
   real(kr),parameter::zero=0
   cwaven=sqrt(2*amassn*amu*ev)*1.e-12_kr/hbar

   !--initialize (e=0.)
   ne=0
   if (e.gt.zero) go to 125
   do i=1,20
      ilast(i)=0
   enddo
   go to 220
  125 continue
   ee=e
   ifirst=ilast(isect)
   inow=ibaset(isect)
   lfw=nint(res(inow+4))

   ! check sign of res(inow+7).  If < 0, it equals -nro and
   ! res(inow+6) marks the beginning of an energy dependent tab1
   ! record for the scattering radius; if > 0 then res(inow+6)
   ! is the cont record beginning with spin data.
   nro=0
   if (res(inow+7).lt.zero) then
      nro=-nint(res(inow+7))
      inow=inow+6+2*nint(res(inow+10))+2*nint(res(inow+11))
   endif

   if (ifirst.gt.0) go to 100
   enext=e
   elast=0
   enxt(isect)=e
   elst(isect)=0
   ie=1
   go to 120

   !--check energy range
  100 continue
   enext=enxt(isect)
   elast=elst(isect)
   ie=ilast(isect)
   if (lfw.eq.0) go to 130
   iest=inow+12
   if (e.lt.enext.and.enext.ge.wide*elast) go to 130
   if (e.le.enext) go to 210
   if (ie.eq.ne) go to 120
  110 continue
   elast=enext
   elst(isect)=elast
   do i=1,4
      slst(isect,i)=snxt(isect,i)
   enddo
   ie=ie+1
   enext=res(iest+ie-1)
   enxt(isect)=enext
   ilast(isect)=ie
  120 continue
   e=enext

   !--compute unresolved cross sections
  130 continue
   do i=2,4
      sigp(i)=0
   enddo
   ! retrieve starting for location for data in a for this section
   inow=ibaset(isect)
   ! retrieve nuclide information
   lfw=nint(res(inow+4))
   naps=nint(res(inow+5))
   if (res(inow+7).gt.0.) then
      ay=res(inow+7)
   else
      iro=inow+6
      ip=2
      ir=1
      call terpa(ay,e,enx,idx,res(iro),ip,ir)
      inow=inow+6+2*nint(res(iro+4))+2*nint(res(iro+5))
   endif
   spi=res(inow+6)
   aaa=res(inow+7)
   if (lfw.ne.1) then
      ne=0
      nls=nint(res(inow+10))
      inow=inow+12
      awri=res(inow)
   else
      ne=nint(res(inow+10))
      nls=nint(res(inow+11))
      ! save starting location for fission width energies
      iest=inow+12
      inow=inow+12+ne
      awri=res(inow)
   endif
   ! calculate channel radius
   rat=awri/(awri+1)
   aw=awri*amassn
   if (naps.eq.0) then
      aa=rc1*aw**third+rc2
   else if (naps.eq.1) then
      aa=ay
   else if (naps.eq.2.and.nro.eq.1) then
      aa=aaa
   else
      call error('csunr1','illegal naps',' ')
   endif
   const=(2*pi**2)/(cwaven*rat)**2

   !--do loop over all l states
   do l=1,nls
      if (lfw.eq.0) njs=nint(res(inow+5))
      if (lfw.eq.1) njs=nint(res(inow+4))
      ll=nint(res(inow+2))
      inow=inow+6

      !do loop  over all j states
      do j=1,njs
         gfx=0
         gxx=0
         if (lfw.le.0) then
            muf=0
            dx=res(inow)
            aj=res(inow+1)
            amun=res(inow+2)
            gnox=res(inow+3)
            ggx=res(inow+4)
            inow=inow+6
         else
            muf=nint(res(inow+3))
            nep6=nint(res(inow+4))
            dx=res(inow+6)
            aj=res(inow+7)
            amun=res(inow+8)
            gnox=res(inow+9)
            ggx=res(inow+10)
            ! save starting location for fission widths
            ifst=inow+12
            inow=inow+6+nep6
         endif
         gj=(2*aj+1)/(4*spi+2)
         mu=int(amun)
         nu= muf
         lamda=0
         ! interpolate parameters at energy of interest
         if (lfw.ne.0) then
            i1=iest
            if (ie.gt.1) i1=iest+ie-2
            i2=ifst
            if (ie.gt.1) i2=ifst+ie-2
            call terp1(res(i1),res(i2),res(i1+1),res(i2+1),e,gfx,2)
         endif
         e2=sqrt(e)
         k=rat*e2*cwaven
         rho=k*aa
         rhoc=k*ay
         ! calculate penetrability (vl) and phase shift(ps)
         call unfac(ll,rho,rhoc,amun,vl,ps)
         vl=vl*e2
         ! calculate potential scattering
         if (j.eq.1) spot=4*pi*(2*ll+1)*(sin(ps)/k)**2
         ! compute cross section contributions
         gnx=gnox*vl
         diff=gxx
         den=e*dx
         temp=const*gj*gnx/den
         terg=temp*ggx
         ters=temp*gnx
         terf=temp*gfx
         ! calculate fluctuation integrals
         call gnrl(gnx,gfx,ggx,mu,nu,lamda,gs ,diff,1)
         call gnrl(gnx,gfx,ggx,mu,nu,lamda,gc ,diff,2)
         call gnrl(gnx,gfx,ggx,mu,nu,lamda,gff,diff,3)
         gc=gc*terg
         gff=gff*terf
         gs=gs*ters
         ! add interference correction
         add=const*gj*2*gnx*sin(ps)**2
         add=add/(e*dx)
         gs=gs-add
         ! cross sections
         sigp(2)=sigp(2)+gs
         sigp(3)=sigp(3)+gff
         sigp(4)=sigp(4)+gc

      !continue j loop
      enddo
      sigp(2)=sigp(2)+spot

   !--continue l loop
   enddo
   sigp(1)=sigp(2)+sigp(3)+sigp(4)
   if (ifirst.eq.0.and.lfw.eq.0) ifirst=1
   if (lfw.eq.0) go to 220
   if (ifirst.gt.0.and.e.lt.enext.and.&
     enext.gt.wide*elast) go to 220

   !--replace snext with sigp.
   do i=1,4
      snxt(isect,i)=sigp(i)
   enddo
   if (ifirst.gt.0.and.e.eq.enext.and.&
     enext.gt.wide*elast) go to 200
   if (ifirst.gt.0) go to 210
   if (ifirst.lt.0) go to 205
   ifirst=-1
   go to 110
  200 continue
   e=ee
   go to 130
  205 continue
   ifirst=1

   !--interpolate.
  210 continue
   e=ee
   do i=2,4
      y1=slst(isect,i)
      y2=snxt(isect,i)
      call terp1(elast,y1,enext,y2,e,sigp(i),2)
   enddo
   sigp(1)=sigp(2)+sigp(3)+sigp(4)
  220 continue
   return
   end subroutine csunr1

   subroutine csunr2(e,sigp)
   !-------------------------------------------------------------------
   ! Unresolved resonance region (format 2).
   ! Single-level Breit-Wigner formalism.
   ! Energy-dependent parameters.
   ! On initial entry (e=0.), initialization flags are set for all
   ! sections.  On subsequent entries (e.gt.0.), cross sections
   ! are either computed at energy grid points or interpolated
   ! from previously computed values.  The energy grid is the
   ! grid given by the evaluator.  If the energy step is too large
   ! (see wide), parameter interpolation is used. in that panel.
   !-------------------------------------------------------------------
   use physics ! provides pi,amassn,amu,hbar,ev
   use endf    ! provides terp1
   use util    ! provides error
   ! externals
   real(kr)::e,sigp(4)
   ! internals
   integer::i,inow,inoww,int,ne,iloc,ie,nls,l,ll,njs,j,jnow
   integer::mu,nu,lamda,i1,i2,nro,ip,ir,iro,idx,naps
   real(kr)::ee,enext,elast,awri,ay,aaa,spi,rat,aw,aa,const,e2,k
   real(kr)::rho,rhoc,aj,amux,amun,amuf,gj,vl,ps,spot
   real(kr)::gnx,gamma,galpha,gbeta,diff,den,enx
   real(kr)::cwaven,dx,gxx,gnox,ggx,gfx
   real(kr)::temp,terg,ters,terf,gs,gc,gff,add,y1,y2
   integer::ifirst=1
   real(kr),parameter::wide=1.26e0_kr
   real(kr),parameter::rc1=.123e0_kr
   real(kr),parameter::rc2=.08e0_kr
   real(kr),parameter::third=.333333333e0_kr
   real(kr),parameter::small=1.e-8_kr
   real(kr),parameter::zero=0
   cwaven=sqrt(2*amassn*amu*ev)*1.e-12_kr/hbar

   !--initialize (e=0.)
   if (e.gt.zero) go to 125
   do i=1,20
      ilast(i)=0
   enddo
   go to 210
  125 continue
   ee=e
   ifirst=ilast(isect)
   inow=ibaset(isect)

   ! check sign of res(inow+7).  If < 0, it equals -nro and
   ! res(inow+6) marks the beginning of an energy dependent tab1
   ! record for the scattering radius; if > 0 then res(inow+6)
   ! is the cont record beginning with spin data.
   nro=0
   if (res(inow+7).lt.0.) then
      nro=-nint(res(inow+7))
      inoww=inow+6+2*nint(res(inow+10))+2*nint(res(inow+11))
   else
      inoww=inow
   endif
   int=nint(res(inoww+20))
   ne=nint(res(inoww+23))
   iloc=inoww+30
   if (ifirst.gt.0) go to 100
   enext=e
   elast=0
   enxt(isect)=e
   elst(isect)=0
   ie=1
   go to 120

   !--check energy range
  100 continue
   enext=enxt(isect)
   elast=elst(isect)
   ie=ilast(isect)
   if (e.lt.enext.and.enext.ge.wide*elast) go to 130
   if (e.le.enext) go to 200
   if (ie.eq.ne) go to 120
  110 continue
   elast=enext
   elst(isect)=elast
   do i=1,4
      slst(isect,i)=snxt(isect,i)
   enddo
   ie=ie+1
   enext=res(iloc+(ie-1)*6)
   enxt(isect)=enext
   ilast(isect)=ie
  120 continue
   e=enext

   !--compute unresolved cross sections
  130 continue
   do i=2,4
      sigp(i)=0
   enddo
   ! retrieve starting location for data in a (for this section)
   inow=ibaset(isect)
   ! retrieve nuclide information
   naps=nint(res(inow+5))
   if (res(inow+7).gt.0.) then
      ay=res(inow+7)
   else
      iro=inow+6
      ip=2
      ir=1
      call terpa(ay,e,enx,idx,res(iro),ip,ir)
      inow=inow+6+2*nint(res(iro+4))+2*nint(res(iro+5))
   endif
   spi=res(inow+6)
   aaa=res(inow+7)
   nls=nint(res(inow+10))
   awri=res(inow+12)
   spi=res(inow+6)
   nls=nint(res(inow+10))
   inow=inow+12
   ! calculate channel radius
   rat=awri/(awri+1)
   aw=awri*amassn
   if (naps.eq.0) then
      aa=rc1*aw**third+rc2
   else if (naps.eq.1) then
      aa=ay
   else if (naps.eq.2.and.nro.eq.1) then
      aa=aaa
   else
      call error('csunr2','illegal naps',' ')
   endif
   const=(2*pi**2)/(cwaven*rat)**2
   e2=sqrt(e)
   k=cwaven*rat*e2
   rho=k*aa
   rhoc=k*ay

   !--do loop over all l states
   do l=1,nls
      njs=nint(res(inow+4))
      ll=nint(res(inow+2))
      inow=inow+6

      !--do loop  over all j states
      do j=1,njs
         aj=res(inow)
         int=nint(res(inow+2))
         ne=nint(res(inow+5))
         amux=res(inow+8)
         amun=res(inow+9)
         amuf=res(inow+11)
         jnow=nint(res(inow+4))+inow+6
         inow=inow+12
         gj=(2*aj+1)/(4*spi+2)
         mu=nint(amun)
         nu=nint(amuf)
         lamda=nint(amux)
         ! interpolate parameters at energy of interest
         i1=inow
         if ((i1+6).lt.jnow) then
            do while (res(i1+6).le.e-small)
               i1=i1+6
            enddo
         endif
         i2=i1+6
         call terp1(res(i1),res(i1+1),res(i2),res(i2+1),e,dx,2)
         call terp1(res(i1),res(i1+2),res(i2),res(i2+2),e,gxx,2)
         call terp1(res(i1),res(i1+3),res(i2),res(i2+3),e,gnox,2)
         call terp1(res(i1),res(i1+4),res(i2),res(i2+4),e,ggx,2)
         call terp1(res(i1),res(i1+5),res(i2),res(i2+5),e,gfx,2)
         if (gxx.lt.small) gxx=0
         if (gfx.lt.small) gfx=0
         ! calculate penetrability (vl) and phase shift(ps)
         call unfac(ll,rho,rhoc,amun,vl,ps)
         vl=vl*e2
         ! calculate potential scattering
         if (j.le.1) then
            spot=4*pi*(2*ll+1)*(sin(ps)/k)**2
         endif
         ! compute cross section contributions
         gnx=gnox*vl
         gamma=ggx
         galpha=gnx
         gbeta=gfx
         diff=gxx
         den=e*dx
         temp=const*gj*gnx/den
         terg=temp*gamma
         ters=temp*gnx
         terf=temp*gbeta
         ! calculate fluctuation integrals
         call gnrl(galpha,gbeta,gamma,mu,nu,lamda,gs ,diff,1)
         call gnrl(galpha,gbeta,gamma,mu,nu,lamda,gc ,diff,2)
         call gnrl(galpha,gbeta,gamma,mu,nu,lamda,gff,diff,3)
         gc=gc*terg
         gff=gff*terf
         gs=gs*ters
         ! add interference correction term
         add=const*gj*2*gnx*(sin(ps))**2
         add=add/(e*dx)
         gs=gs-add
         ! cross sections
         sigp(2)=sigp(2)+gs
         sigp(3)=sigp(3)+gff
         sigp(4)=sigp(4)+gc
         inow=jnow

      !--continue j loop
      enddo
      sigp(2)=sigp(2)+spot

   !--continue l loop
   enddo
   sigp(1)=sigp(2)+sigp(3)+sigp(4)
   if (ifirst.gt.0.and.e.lt.enext.and.&
     enext.ge.wide*elast) go to 210

   !--replace snext with sigp
   do i=1,4
      snxt(isect,i)=sigp(i)
   enddo
   if (ifirst.gt.0.and.e.eq.enext.and.&
     enext.ge.wide*elast) go to 185
   if (ifirst.gt.0) go to 200
   if (ifirst.lt.0) go to 190
   ifirst=-1
   go to 110
  185 continue
   e=ee
   go to 130
  190 continue
   ifirst=1

   !--interpolate.
  200 continue
   e=ee
   do i=2,4
      y1=slst(isect,i)
      y2=snxt(isect,i)
      call terp1(elast,y1,enext,y2,e,sigp(i),int)
   enddo
   sigp(1)=sigp(2)+sigp(3)+sigp(4)
  210 continue
   return
   end subroutine csunr2

   subroutine csaa(e,sigp)
   !-------------------------------------------------------------------
   ! Calculates multilevel Adler-Adler cross sections at energy e
   ! for one section (one isotope-one energy range)
   !-------------------------------------------------------------------
   use physics ! provides pi,amassn,amu,hbar,ev,bk
   use util    ! provides error
   ! externals
   real(kr)::e,sigp(4)
   ! internals
   integer::i,inow,nls,l,njs,j,nlj,n,ist,li
   real(kr)::rpi,cwaven,awri,ap,delta,arat,eroot,esq,ecub,k
   real(kr)::c,omg,snf,csf,bakt,bakf,bakc
   real(kr)::de,dw,gr,gi,x,psi,chi,theta,ax,y,rew,aimw
   integer::ki=1
   real(kr)::cadler=6.5099897e5_kr
   real(kr),parameter::zero=0
   rpi=sqrt(pi)
   cwaven=sqrt(2*amassn*amu*ev)*1.e-12_kr/hbar

   !--compute cross sections for this energy
   do i=1,4
      sigp(i)=0
   enddo
   ! retrieve starting location for data in a
   inow=ibaset(isect)
   ! retrieve nuclide information
   awri=res(inow+12)
   ap=res(inow+7)
   nls=nint(res(inow+10))
   li=nint(res(inow+14))
   inow=inow+18
   ! calculate constants
   if (tempr.gt.zero) delta=1/sqrt(4*bk*tempr*e/awri)
   arat=awri/(awri+1)
   eroot=sqrt(abs(e))
   esq=e*e
   ecub=esq*e
   k=cwaven*arat*eroot
   c=cadler/(arat*arat)
   omg=2*k*ap
   snf=sin(omg)
   csf=cos(omg)
   ! calculate background
   bakt=0
   bakf=0
   bakc=0
   if (li.lt.5) call error('csaa','bad li value',' ')
   if (li.ne.6) then
      bakt=res(inow)+res(inow+1)/e+res(inow+2)/esq
      bakt=bakt+res(inow+3)/ecub+res(inow+4)*e+res(inow+5)*esq
      bakt=bakt*c/eroot
   endif
   inow=inow+6
   if (li.ne.5) then
      bakf=res(inow)+res(inow+1)/e+res(inow+2)/esq
      bakf=bakf+res(inow+3)/ecub+res(inow+4)*e+res(inow+5)*esq
      bakf=bakf*c/eroot
   endif
   inow=inow+6
   bakc=res(inow)+res(inow+1)/e+res(inow+2)/esq
   bakc=bakc+res(inow+3)/ecub+res(inow+4)*e+res(inow+5)*esq
   bakc=bakc*c/eroot
   inow=inow+6

   !--calculate resonance contribution
   do l=1,nls
      njs=nint(res(inow+4))
      inow=inow+6
      do j=1,njs
         nlj=nint(res(inow+5))
         inow=inow+6
         ! loop over individual resonances for this l-j state
         do n=1,nlj
            ! loop over reaction types
            ist=inow
            do i=1,4
               if (i.ne.2) then
                  de=res(ist)
                  dw=res(ist+1)
                  gr=res(ist+2)
                  gi=res(ist+3)
                  ist=ist+4
                  x=(de-e)/dw
                  if (tempr.le.zero) then
                     psi=1/(1+x*x)
                     chi=x*psi
                  else
                     theta=2*dw*delta
                     ax=theta*x/2
                     y=theta/2
                     call quickw(ax,y,rew,aimw,ki,tr,ti)
                     psi=rpi*theta*rew/2
                     chi=rpi*theta*aimw/2
                  endif
                  ! compute cross section contributions
                  if (i.eq.1) sigp(i)=sigp(i)+((gr*csf+gi*snf)*psi+&
                    (gi*csf-gr*snf)*chi)/dw
                  if (i.gt.1) sigp(i)=sigp(i)+(gr*psi+gi*chi)/dw
               endif
            enddo
            inow=ist
         enddo
      enddo
   enddo

   !--add background and factors
   sigp(1)=(sigp(1)+bakt)*c/eroot
   ! total
   sigp(1)=sigp(1)+2*c*(1-csf)/e
   ! fission
   sigp(3)=(sigp(3)+bakf)*c/eroot
   ! capture
   sigp(4)=(sigp(4)+bakc)*c/eroot
   if (li.ne.6) then
      sigp(2)=sigp(1)-(sigp(3)+sigp(4))
   else
      sigp(1)=sigp(3)+sigp(4)
   endif
   return
   end subroutine csaa

   subroutine facphi(l,rho,phi)
   !-------------------------------------------------------------------
   ! Calculates phase shift.
   !-------------------------------------------------------------------
   ! externals
   integer::l
   real(kr)::rho,phi
   ! internals
   real(kr)::r2,r4,top,bot
   real(kr),parameter::test=1.e-6_kr

   r2=rho*rho
   if (l.eq.0) then
      phi=rho
   else if (l.eq.1) then
      phi=rho-atan(rho)
   else if (l.eq.2)then
      phi=rho-atan(3*rho/(3-r2))
      if ((phi/rho).lt.test) phi=0
   else if (l.eq.3) then
      phi=rho-atan((15*rho-rho*r2)/(15-6*r2))
      if ((phi/rho).lt.test) phi=0
   else
      r4=r2*r2
      top=105*rho-10*r2*rho
      bot=105-45*r2+r4
      phi=rho-atan(top/bot)
      if ((phi/rho).lt.test) phi=0
   endif
   return
   end subroutine facphi

   subroutine unfac(l,rho,rhoc,amun,vl,ps)
   !-------------------------------------------------------------------
   ! Calculates the penetrability factor (vl) and phase shift (ps).
   !-------------------------------------------------------------------
   ! externals
   integer::l
   real(kr)::rho,rhoc,amun,vl,ps
   ! internals
   real(kr)::r2,r4

   r2=rho*rho
   if (l.eq.0) then
      vl=amun
      ps=rhoc
   else if (l.eq.1) then
      vl=amun*r2/(1+r2)
      ps=rhoc-atan(rhoc)
   else if (l.eq.2) then
      r4=r2*r2
      vl=amun*r4/(9+3*r2+r4)
      ps=rhoc-atan(3*rhoc/(3-rhoc*rhoc))
   endif
   return
   end subroutine unfac

   subroutine gnrl(galpha,gbeta,gamma,mu,nu,lamda,s,df,id)
   !-------------------------------------------------------------------
   ! Calculates fluctuation integrals for unresolved resonances.
   !-------------------------------------------------------------------
   ! externals
   integer::mu,nu,lamda,id
   real(kr)::galpha,gbeta,gamma,s,df
   ! internals
   integer::j,k,l
   real(kr)::xj
   real(kr),dimension(10,4)::qw=reshape((/&
     1.1120413e-1_kr,2.3546798e-1_kr,2.8440987e-1_kr,&
     2.2419127e-1_kr,0.10967668e0_kr,.030493789e0_kr,&
     0.0042930874e0_kr,2.5827047e-4_kr,4.9031965e-6_kr,&
     1.4079206e-8_kr,0.033773418e0_kr,0.079932171e0_kr,&
     0.12835937e0_kr,0.17652616e0_kr,0.21347043e0_kr,&
     0.21154965e0_kr,0.13365186e0_kr,0.022630659e0_kr,&
     1.6313638e-5_kr,2.745383e-31_kr,3.3376214e-4_kr,&
     0.018506108e0_kr,0.12309946e0_kr,0.29918923e0_kr,&
     0.33431475e0_kr,0.17766657e0_kr,0.042695894e0_kr,&
     4.0760575e-3_kr,1.1766115e-4_kr,5.0989546e-7_kr,&
     1.7623788e-3_kr,0.021517749e0_kr,0.080979849e0_kr,&
     0.18797998e0_kr,0.30156335e0_kr,0.29616091e0_kr,&
     0.10775649e0_kr,2.5171914e-3_kr,8.9630388e-10_kr,&
     0.e0_kr/),(/10,4/))
   real(kr),dimension(10,4)::qp=reshape((/&
     3.0013465e-3_kr,7.8592886e-2_kr,4.3282415e-1_kr,&
     1.3345267e0_kr,3.0481846e0_kr,5.8263198e0_kr,&
     9.9452656e0_kr,1.5782128e1_kr,23.996824e0_kr,&
     36.216208e0_kr,1.3219203e-2_kr,7.2349624e-2_kr,&
     0.19089473e0_kr,0.39528842e0_kr,0.74083443e0_kr,&
     1.3498293e0_kr,2.5297983e0_kr,5.2384894e0_kr,&
     13.821772e0_kr,75.647525e0_kr,1.0004488e-3_kr,&
     0.026197629e0_kr,0.14427472e0_kr,0.44484223e0_kr,&
     1.0160615e0_kr,1.9421066e0_kr,3.3150885e0_kr,&
     5.2607092e0_kr,7.9989414e0_kr,12.072069e0_kr,&
     0.013219203e0_kr,0.072349624e0_kr,0.19089473e0_kr,&
     0.39528842e0_kr,0.74083443e0_kr,1.3498293e0_kr,&
     2.5297983e0_kr,5.2384894e0_kr,13.821772e0_kr,&
     75.647525e0_kr/),(/10,4/))
   real(kr),parameter::zero=0

   s=0
   if (galpha.le.zero) return
   if (gamma.le.zero) return
   if (gbeta.lt.zero) return
   if (gbeta.gt.zero.and.df.lt.zero) return

   if (gbeta.eq.zero.and.df.eq.zero) then
      if (id.eq.1) then
         do j=1,10
            s=s+qw(j,mu)*qp(j,mu)*qp(j,mu)/(galpha*qp(j,mu)+gamma)
         enddo
      else if (id.eq.2) then
         do j=1,10
            s=s+qw(j,mu)*qp(j,mu)/(galpha*qp(j,mu)+gamma)
         enddo
      endif
      return
   endif

   if (gbeta.eq.zero.and.df.gt.zero) then
      if (id.eq.1) then
         do j=1,10
            xj=qp(j,mu)
            do k=1,10
               s=s+qw(j,mu)*qw(k,lamda)*xj*xj/(galpha*xj+gamma&
                 +df*qp(k,lamda))
            enddo
         enddo
      else if (id.eq.2) then
         do j=1,10
            xj=qp(j,mu)
            do k=1,10
               s=s+qw(j,mu)*qw(k,lamda)*xj/(galpha*xj+gamma&
                 +df*qp(k,lamda))
            enddo
         enddo
      endif
      return
   endif

   if (gbeta.gt.zero.and.df.eq.0) then
      if (id.eq.1) then
         do j=1,10
            xj=qp(j,mu)
            do k=1,10
               s=s+qw(j,mu)*qw(k,nu)*xj*xj/(galpha*xj&
                 +gbeta*qp(k,nu)+gamma)
            enddo
         enddo
      else if (id.eq.2) then
         do j=1,10
            xj=qp(j,mu)
            do k=1,10
               s=s+qw(j,mu)*qw(k,nu)*xj/(galpha*xj&
                 +gbeta*qp(k,nu)+gamma)
            enddo
         enddo
      else if (id.eq.3) then
         do j=1,10
            do k=1,10
               s=s+qw(j,mu)*qw(k,nu)*qp(j,mu)*qp(k,nu)&
                 /(galpha*qp(j,mu)+gbeta*qp(k,nu)+gamma)
            enddo
         enddo
      endif
      return
   endif

   if (gbeta.gt.zero.and.df.gt.zero) then
      if (id.eq.1) then
         do j=1,10
            xj=qp(j,mu)
            do k=1,10
               do l=1,10
                  s=s+qw(j,mu)*qw(k,nu)*qw(l,lamda)*xj*xj&
                    /(galpha*xj+gbeta*qp(k,nu)+gamma&
                    +df*qp(l,lamda))
               enddo
            enddo
         enddo
      else if (id.eq.2) then
         do j=1,10
            do k=1,10
               do l=1,10
                  s=s+qw(j,mu)*qw(k,nu)*qw(l,lamda)*qp(j,mu)/&
                    (galpha*qp(j,mu)+gbeta*qp(k,nu)+gamma&
                    +df*qp(l,lamda))
               enddo
            enddo
         enddo
      else if (id.eq.3) then
         do j=1,10
            do k=1,10
               do l=1,10
                  s=s+qw(j,mu)*qw(k,nu)*qw(l,lamda)*qp(j,mu)&
                    *qp(k,nu)/(galpha*qp(j,mu)+gbeta*qp(k,nu)&
                    +gamma+df*qp(l,lamda))
               enddo
            enddo
         enddo
      endif
      return
   endif
   return
   end subroutine gnrl

   subroutine emerge(nin,ngrid,ngp,nres,nrtot,iold,inew,nscr)
   !-------------------------------------------------------------------
   ! Merge energies from input grid with resonance grid (if any).
   ! Compute all cross sections on full grid and add resonance cross
   ! sections where needed.  Accumulate the total cross sections.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use endf   ! provides endf routines and variables
   use util   ! provides loada,finda,openz,repoz,closz
   ! externals
   integer::nin,ngrid,ngp,nres,nrtot,iold,inew,nscr
   ! internals
   integer::nneg,ntot,i,in,ig,inn,nss,iss,nb,nw,idis,it
   integer::imtr,itt,k,istart,iend,j,ib,isave,ir,ith
   real(kr)::er,eg,en,e,thresh,sn,sg,enext
   real(kr)::res(nsig+1),tot(10)
   real(kr)::aa(1)
   real(kr),dimension(:),allocatable::bufo,bufn,bufg,bufr
   real(kr),dimension(:),allocatable::scr
   real(kr),parameter::test=1.e-10_kr
   real(kr),parameter::finity=.99e12_kr
   real(kr),parameter::small=1.e-8_kr
   real(kr),parameter::zero=0
   real(kr),parameter::one=1

   !--set up storage allocations.
   allocate(bufo(nbuf))
   allocate(bufn(nbuf))
   allocate(bufg(nbufg))
   allocate(bufr(nbufr))
   allocate(scr(npage+50))
   nneg=0
   ntot=nmtr+1
   do i=1,nmtr
      tot(i+1)=0
   enddo

   !--assign scratch units.
   iold=14
   inew=15

   !--merge background grid with resonance grid
   in=0
   ig=1
   call finda(ig,aa,1,ngrid,bufg,nbufg)
   eg=aa(1)
   er=finity
   if (nrtot.gt.0) then
      ir=1
      call finda(ir,res,nsig+1,nres,bufr,nbufr)
      er=res(1)
   endif
  130 continue
   in=in+1
   en=eg
   if (er.lt.eg) en=er
   tot(1)=en
   inn=in
   if (ig.eq.ngp.and.en.eq.eg) inn=-in
   call loada(inn,tot,ntot,iold,bufo,nbuf)
  135 continue
   if (en.ne.er) go to 150
   if (ir.lt.nrtot) go to 140
   er=finity
   go to 150
  140 continue
   ir=ir+1
   call finda(ir,res,nsig+1,nres,bufr,nbufr)
   er=res(1)
   go to 135
  150 continue
   if (en.ne.eg) go to 130
   if (ig.eq.ngp) go to 160
   ig=ig+1
   call finda(ig,aa,1,ngrid,bufg,nbufg)
   eg=aa(1)
   go to 150
  160 continue
   ngo=in

   !--find next reaction on input file
   call openz(nscr,1)
   call repoz(nin)
   call repoz(nscr)
   nsc=0
  210 continue
   call contio(nin,0,nscr,scr,nb,nw)
   nss=1
   if (mfh.eq.10) nss=n1h
   iss=nss
   if (math.gt.0) go to 220
   call atend(0,nscr)
   go to 510
  220 continue
   if (mth.eq.0) go to 210
   e=0
   call gety1(e,thresh,idis,sn,nin,scr)
   thresh=sigfig(thresh,7,0)

   !--identify resonance reactions
   itype=0
   if (nrtot.ne.0.and.mfh.eq.3) then
      if (mth.eq.2) itype=2
      if (mth.eq.18) itype=3
      if (mth.eq.19) itype=3
      if (mth.eq.102) itype=4
      if (nmtres.gt.2) then
         do i=3,nmtres
            if (mth.ne.18.and.mth.ne.19.and.mth.eq.mmtres(i)) then
               itype=2+i
            endif
         enddo
      endif
      if (lrx.ne.0.and.mth.eq.51) itype=5
   endif

   !--set up first resonance point.
   er=finity
   if (itype.ne.0) then
      ir=1
      call finda(ir,res,nsig+1,nres,bufr,nbufr)
      er=res(1)
      if (er.lt.thresh) thresh=sigfig(er,7,0)
   endif
   in=0
   ith=0
   it=0
   ig=0

   !--compute cross section at next grid point
  340 continue
   ig=ig+1
   call finda(ig,tot,ntot,iold,bufo,nbuf)
   eg=tot(1)
   sg=tot(2)
   sn=0
   if (thresh-eg.gt.test*thresh.and.itype.eq.0) go to 370
   call gety1(eg,enext,idis,sn,nin,scr)
   if (thresh.gt.one.and.abs(thresh-eg).lt.test*thresh) sn=0
   ! backgrounds in a range of unresolved-smooth overlap
   ! are arbitrarily assigned to the unresolved component
   if (eg.ge.eresr.and.eg.lt.eresh.and.itype.gt.0) sn=0
  345 continue
   if (er-eg.gt.test*eg) go to 360
   if (abs(eg-er).lt.test*eg) go to 355
   ir=ir+1
   call finda(ir,res,nsig+1,nres,bufr,nbufr)
   er=res(1)
   go to 345
  355 continue
   sn=sn+res(1+itype)
   if (ir.ge.nrtot) then
      er=finity
   else
      ir=ir+1
      call finda(ir,res,nsig+1,nres,bufr,nbufr)
      er=res(1)
   endif
  360 continue
   in=in+1
   if (itype.eq.2.and.sn.le.small) then
      nneg=nneg+1
      if (nneg.eq.1) then
         if (sn .le. 0.0) then
            call mess('emerge',&
            'nonpositive elastic cross sections found.',' ')
         else
            call mess('emerge',&
            'tiny elastic cross sections replaced by the lower limit.',' ')
         end if
      endif
      sn=small
   endif
   sn=sigfig(sn,7,0)
   tot(2)=sn
   if (ith.eq.0.and.sn.gt.zero) ith=in
   inn=in
   if (ig.eq.ngo) inn=-in
   call loada(inn,tot,2,ngrid,bufg,nbufg)
   tot(2)=sg

   !--accumulate redundant reactions.
  370 continue
   it=it+1
   if (thresh-eg.gt.test*thresh.and.itype.eq.0) go to 440
   if (mth.eq.10) go to 440
   if (mfh.ne.3) go to 430
   if (mth.ge.46.and.mth.le.49) go to 440
   if (mth.gt.200.and.mth.lt.mpmin) go to 440
   do 410 imtr=1,nmtr
   if (mtr(imtr).ne.1) go to 375
   go to 400
  375 continue
   if (mtr(imtr).ne.3) go to 380
   if (mth.eq.2) go to 410
   go to 400
  380 continue
   if (mtr(imtr).ne.4) go to 381
   if (mth.lt.51.or.mth.gt.91) go to 410
   if (sn.eq.zero) go to 410
   go to 400
  381 continue
   if (mtr(imtr).ne.103) go to 382
   if (mth.lt.mpmin.or.mth.gt.mpmax) go to 410
   if (sn.eq.zero) go to 410
   go to 400
  382 continue
   if (mtr(imtr).ne.104) go to 383
   if (mth.lt.mdmin.or.mth.gt.mdmax) go to 410
   if (sn.eq.zero) go to 410
   go to 400
  383 continue
   if (mtr(imtr).ne.105) go to 384
   if (mth.lt.mtmin.or.mth.gt.mtmax) go to 410
   if (sn.eq.zero) go to 410
   go to 400
  384 continue
   if (mtr(imtr).ne.106) go to 385
   if (mth.lt.m3min.or.mth.gt.m3max) go to 410
   if (sn.eq.zero) go to 410
   go to 400
  385 continue
   if (mtr(imtr).ne.107) go to 386
   if (mth.lt.m4min.or.mth.gt.m4max) go to 410
   if (sn.eq.zero) go to 410
   go to 400
  386 continue
   if (mtr(imtr).ne.18) go to 410
   if (mth.lt.19.or.mth.gt.21.and.mth.ne.38) go to 410
  400 continue
   if (it.lt.mtrt(imtr)) then
        mtrt(imtr)=it
       if (mtrt(imtr).gt.1) mtrt(imtr)=mtrt(imtr)-1
   endif
   tot(1+imtr)=tot(1+imtr)+sn
  410 continue
  430 continue
   if (mfh.ne.23) go to 440
   do 435 imtr=1,nmtr
   if (mtr(imtr).ne.501) go to 433
   if (mth.eq.515.or.mth.eq.517) go to 435
   go to 434
  433 continue
   if (mtr(imtr).ne.522) go to 435
   if (mth.lt.534.or.mth.gt.572) go to 435
  434 continue
   if (it.lt.mtrt(imtr)) then
       mtrt(imtr)=it
       if (mtrt(imtr).gt.1) mtrt(imtr)=mtrt(imtr)-1
   endif
   tot(1+imtr)=tot(1+imtr)+sn
  435 continue
  440 continue
   itt=it
   if (ig.eq.ngo) itt=-it
   call loada(itt,tot,ntot,inew,bufn,nbuf)
   if (ig.lt.ngo) go to 340

   !--finished with this reaction
   if (ith.gt.1) ith=ith-1
   if (iss.eq.nss) nxc=nxc+1
   mfs(nxc)=mfh
   mts(nxc)=mth
   if (iss.eq.nss) ncs(nxc)=1
   ncs(nxc)=ncs(nxc)+2+int(((in-ith+1)+2)/3)
   scr(5)=1
   scr(6)=in-ith+1
   scr(7)=in-ith+1
   scr(8)=2
   k=8
   if (ith.gt.0) then
      do ib=1,ith
         call finda(ib,tot,2,ngrid,bufg,nbufg)
      enddo
   endif
   istart=ith
   nb=1
   do while (nb.ne.0)
      iend=in
      if ((iend-istart).ge.npage/2) iend=istart+npage/2-1
      j=k-1
      ib=istart-1
      do while (ib.lt.iend)
         j=j+2
         ib=ib+1
         call finda(ib,tot,2,ngrid,bufg,nbufg)
         scr(j)=tot(1)
         scr(j+1)=tot(2)
      enddo
      nw=j+1
      if (k.ne.0) then
         k=0
         call tab1io(0,0,nscr,scr,nb,nw)
      else
         call moreio(0,0,nscr,scr,nb,nw)
      endif
      if (nb.ne.0) istart=iend+1
   enddo
   iss=iss-1
   if (iss.gt.0) go to 220
   call tosend(nin,0,0,scr)
   call asend(0,nscr)
   isave=iold
   iold=inew
   inew=isave
   go to 210

   !--merge is finished.
  510 continue

   if (nneg.gt.0) write(nsyso,'(/&
     &'' number of nonpositive cross sections removed ='',i8)') nneg
   write(nsyso,'(&
     &'' number of points in final unionized grid     ='',i8)') ngo
   deallocate(bufr)
   deallocate(bufg)
   deallocate(bufn)
   deallocate(bufo)
   call closz(-inew)
   call closz(nin)
   call closz(-ngrid)
   call closz(-nres)
   return
   end subroutine emerge

   subroutine recout(iold,nscr,nrtot)
   !-------------------------------------------------------------------
   ! Add a new material to the output pendf tape.
   !-------------------------------------------------------------------
   use endf   ! provides endf routines and variable
   use util   ! provides error,closz,repoz
   ! externals
   integer::iold,nscr,nrtot
   ! internals
   integer::i152,nb,nw,nwd,i,j,nc,no2,no3,imtr,np,nxcc
   integer::l,ntot,lis3,lfs,mtl,istart,k,last,iend,ib,mtd
   integer::mfl,n1l,n2l
   integer::iang,imt,idone
   real(kr)::resl(1+ncoef*nmtres)
   real(kr),dimension(:),allocatable::bufo,bufn,bufl
   real(kr),dimension(:),allocatable::scr,dict
   real(kr)::tot(10)
   character::strng*60
   character(4)::tz(17)
   real(kr)::z(17)
   equivalence(tz(1),z(1))
   real(kr),parameter::emin=1.e-5_kr
   real(kr),parameter::zero=0

   !--initialize.
   allocate(scr(npage+50))
   allocate(bufo(nbuf))
   allocate(bufn(nbuf))
   allocate(bufl(nbufl))
   i152=0
   if (nunr.gt.0) i152=1
   iang=16

!  --set flag for presence of file 3 with incident neutrons.  if
!    none found increment nxc since we'll insert a dummy mf3, mt1
!    section into the output tape later.
   no3=0
   if (nint(zain).eq.1) then
      do i=1,nxc
         if (mfs(i).eq.3)no3=1
      enddo
      if (no3.eq.0) nxc=nxc+1
   endif

   !--write control records.
   scr(1)=za
   scr(2)=awr
   if (iverf.eq.6) lrp=2
   scr(3)=lrp
   scr(4)=lfi
   scr(5)=0
   scr(6)=0
   if (iverf.eq.4) scr(6)=nxc+nmtr+i152
   math=mata
   mfh=1
   mth=451
   nsh=1
   nw=6
   call contio(0,nout,0,scr,nb,nw)
   if (iverf.ne.4) then
      scr(1)=elis
      scr(2)=sta
      scr(3)=lis
      scr(4)=lis0
      scr(5)=0
      scr(6)=nfor
      if (iverf.eq.6) scr(6)=6
      call contio(0,nout,0,scr,nb,nw)
   endif
   if (iverf.eq.6) then
      scr(1)=awin
      scr(2)=efmax
      scr(3)=lrel
      scr(4)=0
      scr(5)=int(10*zain)
      scr(6)=nver
      call contio(0,nout,0,scr,nb,nw)
   endif
   scr(1)=tempr
   if (tempr.eq.0) scr(1)=tempi
   scr(2)=err
   scr(3)=0
   if (iverf.eq.6) scr(3)=1
   scr(4)=0
   nwd=ncards*17
   scr(5)=nwd
   scr(6)=0
   if (iverf.ge.5) scr(6)=nxc+nmtr+i152
   ncs(1)=ncs(1)+nxc+nmtr+i152
   if (ncoef.gt.1) then
      scr(6)=scr(6)+nmtres-1
      ncs(1)=ncs(1)+nmtres-1
   endif

   !--write descriptive information and dictionary.
   nw=6
   nc=0
   do i=1,ncards
      do j=1,17
         tz(j)=card(j+nc)
      enddo
      do j=1,17
         scr(j+nw)=z(j)
      enddo
      nw=nw+17
      nc=nc+17
   enddo
   nxcc=nxc
   if (nint(zain).eq.1.and.no3.eq.0) nxcc=nxcc-1
   scr(6)=nxcc+nmtr+i152
   call hdatio(0,nout,0,scr,nb,nw)
   no2=1
   if (nxc.ne.0) then
      nw=(nxc+nmtr+i152)*6
      if (ncoef.gt.1) nw=nw+(nmtres-1)*6
      allocate(dict(nw))
      j=0
      imtr=1
      idone=0
      do i=1,nxcc
         if (mfs(i).eq.2) no2=0
         if (mfs(i).eq.3.or.mfs(i).eq.23) then
            do while (imtr.le.nmtr.and.mts(i).ge.mtr(imtr))
               dict(j+1)=0
               dict(j+2)=0
               dict(j+3)=mfs(i)
               if (mfs(i).eq.3.and.mtr(imtr).eq.1.and.awin.eq.0) then
                  dict(j+4)=3
               else
                  dict(j+4)=mtr(imtr)
               endif
               np=ngo-mtrt(imtr)+1
               dict(j+5)=3+int((np+2)/3)
               dict(j+6)=0
               j=j+6
               imtr=imtr+1
            enddo
         endif
         if (no3.eq.0.and.nint(zain).eq.1.and.mfs(i).gt.3)then
             no3=-1
             dict(j+1)=0
             dict(j+2)=0
             dict(j+3)=3
             dict(j+4)=1
             dict(j+5)=4
             dict(j+6)=0
             j=j+6
        endif
        if (ncoef.gt.1.and.(i.eq.nxc.or.mfs(i).gt.4).and.idone.eq.0) then
            idone=1
            do k=1,nmtres
               if (k.eq.1.or.k.gt.2) then
                  dict(j+1)=0
                  dict(j+2)=0
                  dict(j+3)=4
                  if (k.eq.1) then
                     dict(j+4)=mmtres(1)
                     dict(j+5)=mcards(1)
                  endif
                  if (k.gt.2) then
                     dict(j+4)=mmtres(k)
                     dict(j+5)=mcards(k)
                  endif
                  dict(j+6)=0
                  j=j+6
               endif
            enddo
         endif
         dict(j+1)=0
         dict(j+2)=0
         dict(j+3)=mfs(i)
         dict(j+4)=mts(i)
         dict(j+5)=ncs(i)
         dict(j+6)=0
         j=j+6
         if (mfs(i).eq.2.and.i152.ne.0) then
            dict(j+1)=0
            dict(j+2)=0
            dict(j+3)=2
            dict(j+4)=152
            dict(j+5)=3+nunr
            dict(j+6)=0
            j=j+6
         endif
      enddo
      nw=nxcc+nmtr+i152
      if (ncoef.gt.1) then
         nw=nw+nmtres-1
      endif
      call dictio(0,nout,0,dict,nb,nw)
      deallocate(dict)
   endif
   call asend(nout,0)
   call afend(nout,0)

   !--write default mf2/mt151.
   !--if unresolved present, write mf2/mt152.
   if (no2.le.0) then
      scr(1)=za
      scr(2)=awr
      scr(3)=0
      scr(4)=0
      scr(5)=1
      scr(6)=0
      mfh=2
      mth=151
      call contio(0,nout,0,scr,nb,nw)
      scr(2)=1
      call contio(0,nout,0,scr,nb,nw)
      scr(1)=emin
      ! use top of resolved range for upper limit
      scr(2)=eresh
      if (eresr.lt.eresh) scr(2)=eresr
      scr(5)=0
      call contio(0,nout,0,scr,nb,nw)
      scr(1)=spin
      scr(2)=ascat
      call contio(0,nout,0,scr,nb,nw)
      call asend(nout,0)
      if (i152.gt.0) then
         mfh=2
         mth=152
         call contio(0,nout,0,sunr,nb,nw)
         l=7
         call listio(0,nout,0,sunr(l),nb,nw)
         do while (nb.ne.0)
            l=l+nw
            call moreio(0,nout,0,sunr(l),nb,nw)
         enddo
         call asend(nout,0)
      endif
      call afend(nout,0)
   endif

   !--write file 3 using reactions from scratch tape
   !--and adding desired redundant reactions.
   ntot=nmtr+1
   call repoz(nscr)
   call contio(nscr,0,0,scr,nb,nw)
   za=c1h
   awr=c2h
   lis3=l1h
   lfs=l2h
   n1l=n1h
   n2l=n2h
   mfl=mfh
   mtl=mth
   !--test mfl=0 (occurs when all mf3 data are redundant and removed)
   if (nint(zain).eq.1.and.(mfl.eq.0.or.mfl.gt.3)) then
      scr(3)=0
      scr(4)=99
      scr(5)=0
      scr(6)=0
      mfh=3
      mth=1
      call contio(0,nout,0,scr,nb,nw)
      scr(1)=tempr
      scr(2)=0
      scr(3)=0
      scr(4)=0
      scr(5)=1
      scr(6)=2
      scr(7)=2
      scr(8)=2
      scr(9)=emin
      scr(10)=0
      scr(11)=eresh
      scr(12)=0
      nw=12
      call tab1io(0,nout,0,scr,nb,nw)
      call asend(nout,0)
      call afend(nout,0)
   !--if no other mf3 sections on scratch, copy the rest
      if(mfl.eq.0) go to 272
      goto 270
   endif
   mth=1
   if (awin.eq.0) then
      mth=3
   endif
   if (mfh.eq.23) mth=501
   scr(3)=0
   scr(4)=99
   imtr=1
   call contio(0,nout,0,scr,nb,nw)

   !--create a new section for a redundant reaction.
  210 continue
   istart=mtrt(imtr)
   np=ngo-istart+1
   if (np.le.0) then
      write(strng,'(''for mf'',i2,'' mt'',i3,&
        &'' np:'',i6,''-'',i6,''+1='',i6)') mfh,mth,ngo,istart,np
      call error('recout',strng,' ')
   endif
   scr(1)=tempr
   scr(2)=0
   if (mth.eq.18) scr(2)=sigfig(q18,7,0)
   scr(3)=0
   scr(4)=0
   scr(5)=1
   scr(6)=np
   scr(7)=np
   scr(8)=2
   k=8
   if (istart.gt.1) then
      last=istart-1
      i=0
      do while (i.lt.last)
         i=i+1
         call finda(i,tot,ntot,iold,bufo,nbuf)
      enddo
   endif
   nb=1
   do while (nb.ne.0)
      iend=ngo
      if ((iend-istart).ge.npage/2) iend=istart+npage/2-1
      j=k-1
      ib=istart-1
      do while (ib.lt.iend)
         j=j+2
         ib=ib+1
         call finda(ib,tot,ntot,iold,bufo,nbuf)
         scr(j)=tot(1)
         scr(j+1)=sigfig(tot(imtr+1),7,0)
      enddo
      nw=j+1
      if (k.ne.0) then
         k=0
         call tab1io(0,nout,0,scr,nb,nw)
      else
         call moreio(0,nout,0,scr,nb,nw)
      endif
      if (nb.ne.0) istart=iend+1
   enddo
   call asend(nout,0)
   imtr=imtr+1
   if (imtr.gt.nmtr) go to 270
   mtd=mtr(imtr)
  250 continue
   if (mtl.gt.mtd) go to 260

   !--this is a non-redundant reaction.
   mth=mtl
   scr(1)=za
   scr(2)=awr
   scr(3)=lis3
   scr(4)=lfs
   scr(5)=0
   scr(6)=0
   call contio(0,nout,0,scr,nb,nw)
   call tosend(nscr,nout,0,scr)
   call contio(nscr,0,0,scr,nb,nw)
   lis3=l1h
   lfs=l2h
   mtl=mth
   go to 250

   !--this is a redundant reaction.
  260 continue
   mth=mtd
   scr(1)=za
   scr(2)=awr
   scr(3)=0
   scr(4)=99
   scr(5)=0
   scr(6)=0
   call contio(0,nout,0,scr,nb,nw)
   go to 210

   !--all remaining reactions are non-redundant.
  270 continue
   scr(1)=za
   scr(2)=awr
   scr(3)=lis3
   scr(4)=lfs
   scr(5)=n1l
   scr(6)=n2l
   mfh=mfl
   mth=mtl
   call contio(0,nout,0,scr,nb,nw)
   call tofend(nscr,nout,0,scr)

   !--process angular coefficients
  272 continue
   if (ncoef.gt.1) then
      do imt=2,nmtres
         if (imt.eq.2) l=1
         if (imt.gt.2) l=imt
         k=0
         do i=1,nrtot
            call finda(i,resl,1+ncoef*nmtres,iang,bufl,nbufl)
            if (resl(1+l+nmtres).ne.zero) exit
            k=i
         enddo
         mfh=4
         mth=mmtres(imt)
         if (imt.eq.2) mth=mmtres(1)
         scr(1)=za
         scr(2)=awr
         scr(3)=0
         scr(4)=1
         scr(5)=0
         scr(6)=0
         call contio(0,nout,0,scr,nb,nw)
         scr(1)=0
         scr(2)=awr
         scr(3)=0
         scr(4)=2
         scr(5)=0
         scr(6)=0
         call contio(0,nout,0,scr,nb,nw)
         scr(1)=0
         scr(2)=0
         scr(3)=0
         scr(4)=0
         scr(5)=1
         scr(6)=nrtot-k
         scr(7)=nrtot-k
         scr(8)=2
         nw=8
         call tab2io(0,nout,0,scr,nb,nw)
         if (imt.eq.2) l=1
         if (imt.gt.2) l=imt
         do i=k+1,nrtot
            call finda(i,resl,1+ncoef*nmtres,iang,bufl,nbufl)
            scr(1)=0
            scr(2)=resl(1)
            scr(3)=0
            scr(4)=0
            scr(5)=0
            scr(6)=0
            if (resl(1+l).ne.zero) then
               do j=2,ncoef
                  if (resl(1+l+(j-1)*nmtres).ne.zero) then
                     scr(5+j)=resl(1+l+(j-1)*nmtres)
                     scr(5)=scr(5)+1
                  endif
               enddo
            else
               scr(5)=1
               scr(7)=0
            endif
            call listio(0,nout,0,scr,nb,nw)
         enddo
         call asend(nout,0)
      enddo
      call afend(nout,0)
   endif
   call tomend(nscr,nout,0,scr)

   call closz(-iold)
   call closz(nscr)
   deallocate(scr)
   deallocate(bufn)
   deallocate(bufo)
   return
   end subroutine recout

   subroutine quickw(ax,y,rew,aimw,ki,tr,ti)
   !-------------------------------------------------------------------
   ! Used to calculate Psi and Chi line shape functions.
   !-------------------------------------------------------------------
   use physics ! provides pi
   ! externals
   integer::ki
   real(kr)::ax,y,rew,aimw
   real(kr)::tr(62,62),ti(62,62)
   ! internals
   integer::ii,jj,i,j,n
   real(kr)::rpi,aki,x,test,p,q,p2,q2,pq,hp,hq,hq2,hp2
   real(kr)::a1,a2,a3,a4,a5,d1,d2
   real(kr),parameter::break1=36.e0_kr
   real(kr),parameter::break2=144.e0_kr
   real(kr),parameter::break3=10000.e0_kr
   real(kr),parameter::c1=.2752551e0_kr
   real(kr),parameter::c2=2.724745e0_kr
   real(kr),parameter::c3=.5124242e0_kr
   real(kr),parameter::c4=.05176536e0_kr
   real(kr),parameter::c5=1.1283792e0_kr
   real(kr),parameter::zero=0
   rpi=sqrt(pi)

   aki=1
   if (ax.lt.zero) aki=-1
   x=abs(ax)
   test=x*x+y*y
   if (test.lt.break1) then
      ii=int(x*10)
      jj=int(y*10)
      i=ii+2
      j=jj+2
      n=j-1
      p=10*x-ii
      q=10*y-jj
      p2=p*p
      q2=q*q
      pq=p*q
      hp=p/2
      hq=q/2
      hq2=q2/2
      hp2=p2/2
      a1=hq2-hq
      a2=hp2-hp
      a3=1+pq-p2-q2
      a4=hp2-pq+hp
      a5=hq2-pq+hq
      rew=a1*tr(i,n)+a2*tr(i-1,j)+a3*tr(i,j)+a4*tr(i+1,j)&
        +a5*tr(i,j+1)+pq*tr(i+1,j+1)
      if (ki.gt.0) then
         aimw=a1*ti(i,n)+a2*ti(i-1,j)+a3*ti(i,j)+a4*ti(i+1,j)&
           +a5*ti(i,j+1)+pq*ti(i+1,j+1)
         aimw=aimw*aki
      endif
   else if (test.lt.break2) then
      a1=x**2-y**2
      a2=2*x*y
      a3=a2**2
      a4=a1-c1
      a5=a1-c2
      d1=c3/(a4**2+a3)
      d2=c4/(a5**2+a3)
      rew=d1*(a2*x-a4*y)+d2*(a2*x-a5*y)
      if (ki.gt.0) then
         aimw=d1*(a4*x+a2*y)+d2*(a5*x+a2*y)
         aimw=aimw*aki
      endif
   else if (test.lt.break3) then
      a1=(x**2-y**2)*2
      a2=4*x*y
      a4=a1-1
      d1=c5/(a4**2+a2**2)
      rew=d1*(a2*x-a4*y)
      if (ki.gt.0) then
         aimw=d1*(a4*x+a2*y)
         aimw=aimw*aki
      endif
   else
      a1=1/(rpi*test)
      rew=y*a1
      if (ki.gt.0) aimw=x*a1*aki
   endif
   return
   end subroutine quickw

   subroutine wtab(rw,aimw)
   !-------------------------------------------------------------------
   ! Generate table of complex probability integral w.
   !-------------------------------------------------------------------
   ! externals
   real(kr)::rw(62,62),aimw(62,62)
   ! internals
   integer::nx,ny,nx1,ny1,i,j
   real(kr)::x0,y0,dx,dy,xi,yj,rwt,aimt
   real(kr)::x(62),y(62)
   real(kr),parameter::tenth=.1e0_kr

   x0=-tenth
   y0=-tenth
   dx=tenth
   dy=tenth
   nx=62
   ny=62
   y(1)=y0
   nx1=nx-1
   ny1=ny-1
   do i=1,ny1
      y(i+1)=i*dy+y0
   enddo
   x(1)=x0
   do i=1,nx1
      x(i+1)=i*dx+x0
   enddo
   do i=1,nx
      do j=1,ny
         xi=x(i)
         yj=y(j)
         call w (xi,yj,rwt,aimt)
         rw(i,j)=rwt
         aimw(i,j)=aimt
      enddo
   enddo
   return
   end subroutine wtab

   subroutine w(rez,aim1,rew,aimw)
   !-------------------------------------------------------------------
   ! Compute complex probability integral w.
   !-------------------------------------------------------------------
   use physics ! provides pi
   ! externals
   real(kr)::rez,aim1,rew,aimw
   ! internals
   integer::kw,k,jsig
   real(kr)::tempc,tempd,r2,ai2,aimz
   real(kr)::rv,ak,el,a,b,h,g,am,aak
   real(kr)::prr,pii,amagn,temp1,temp2,aj,c,d
   real(kr)::ajtemp,temp4,ajp,ajsig,sigp,expon,expc,exps
   real(kr)::sig2p,aj4sig,temp3,tt4,temp7
   real(kr)::ref,aimf,rpi,abrez,tempm,temel,aj4sm1
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
   tempc=0
   tempd=0
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

   !--wa is obtained from asymtotic series
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
   prr=rew
   pii=aimw
   amagn=tempm**2+temel**2
   rew=(tempc*tempm+tempd*temel)/amagn
   aimw=(tempm*tempd-temel*tempc)/amagn
   if (abs(rew-prr).ge.eps) go to 380
   if (abs(aimw-pii).ge.eps) go to 380
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
   prr=rew
   pii=aimw
   rew=expc-ref
   aimw=exps-aimf
   if (abs(rew-prr).ge.eps) go to 470
   if (abs(aimw-pii).ge.eps) go to 470
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
   if ((abs(tempm)+abs(temel)-up).lt.zero) go to 500
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
   if ((abs(tempm)+abs(temel)-dn).gt.zero) go to 520
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
  530 continue
   return
   end subroutine w

end module reconm


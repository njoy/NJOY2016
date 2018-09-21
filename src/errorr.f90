module errorm
   ! provides error for NJOY2016
   use locale
   implicit none
   private
   public errorr

   ! global variables

   ! user input quantities
   integer::nendf,npend,nout,nin,ngout,nstan
   integer::ngoutu,nscr5
   integer::matd
   integer::ign
   integer::iprint
   integer::irelco
   integer::iread
   integer::mfcov
   integer::irespr
   integer::legord
   integer::ifissp

   real(kr)::tempin

   integer::nitape,notape

   integer::mfd,mtd

   ! scratch files
   integer::nscr1,nscr2,nscr3

   ! group structure
   integer::ngn,nwgp
   real(kr),dimension(:),allocatable::egn

   ! weight function
   integer::iwt
   integer::mwtf
   real(kr) eb,tb,ab,ec,tc,ac
   real(kr),dimension(:),allocatable::wght

   real(kr)::za,awr
   integer::imode
   integer::mf32,mf33
   integer::nunion
   integer::nresg
   integer::nga
   integer::lct4,lct34,ltt4,ltt34
   real(kr)::ety1,ety2
   real(kr)::emaxx,eltt4

   ! fission spectrum processing
   real(kr)::eclo,echi,efmean
   integer::mffis,igflag,ncove,ncovl
   integer::nendf2

   ! scattering radius uncertainty
   real(kr)::dap,dap3(5)
   integer::isr,isru

   integer,parameter::nenimx=5000
   integer::neni
   real(kr),dimension(:),allocatable::eni

   integer::nlump
   integer::nlmt=50
   integer,dimension(:),allocatable::lmt
   integer,dimension(:),allocatable::lump
   integer,dimension(:),allocatable::iga

   integer,parameter::nmtmax=150
   integer::nmt,nmt1,mats(nmtmax),mts(nmtmax)

   integer::mzap(nmtmax)
   integer::lfs

   integer,parameter::nkmax=50
   integer::nek
   real(kr)::ek(nkmax)

   integer,parameter::nasmax=5
   integer::nas
   integer,dimension(nasmax)::matb,mtb,matc,mtc

   real(kr)::sigz(10)
   integer::nsigz

   real(kr),dimension(:),allocatable::flx,sig,cov
   real(kr),dimension(:),allocatable::un
   real(kr),dimension(:,:,:),allocatable::akxy
   real(kr),dimension(:),allocatable::cff,cfg,cgg
   real(kr),dimension(:),allocatable::cee,cef,ceg,ctt
   real(kr),dimension(:),allocatable::ufg,uef,ueg,uff,ugg,uee,utt
   real(kr),dimension(:),allocatable::cflx
   real(kr),dimension(:,:),allocatable::csig

   real(kr)::u1lele(10),plele(2501,10)
   integer::nr,nbt(20),jnt(20)

   ! globals for resonance processing
   integer::ifresr,ifunrs
   integer::nj,nis,lfw,nls,lrf,nro,naps,nsrs,nvs1
   integer::mpar,npar,il2,il3,nrb,lru,nlrs,lcomp,nls1,lptr
   integer::isrr,ier
   integer::nscr6
   integer::ll
   integer::nlspepi(20)
   real(kr)::ap,arat,ra,spifac,cwaven
   real(kr)::abn,ehr,elr,ehg,elg,spi,ral,apl
   real(kr)::ajmin,gj(10),diff

   integer,parameter::ndig=6

   integer::nscr4

   integer::nmtres
   integer::mmtres(10)

   real(kr),dimension(:,:,:,:),allocatable::crr

   !--speed up options
   real(kr)::eskip1,eskip2,eskip3,eskip4

   !--control for sammy resonance processing
   !integer::isammy=0 ! sammy processing off
   integer::isammy=1 ! sammy processing on

contains

   subroutine errorr
   !--------------------------------------------------------------------
   !
   ! Produce cross section covariances from error files in ENDF format
   !
   ! First, the union energy grid of the users group structure
   ! and the ENDF covariance energies is determined.  The array
   ! of coefficients for derived cross sections is also constructed.
   ! Then multigroup cross sections are computed on the union
   ! grid (see grpav), or they are read from a multigroup cross
   ! section library and then collapsed to the union grid.  The
   ! methods of groupr are used for cross section averaging.  ENDF
   ! covariances and the group cross sections are then combined
   ! to get the basic covariance matrices (see covcal).  Finally,
   ! the basic matrices are combined to get covariances for
   ! derived reactions, the matrices are collapsed to the user's
   ! group structure, and the results are printed and/or written
   ! onto an output gendf tape for later use (see covout).
   !
   ! If resonance parameter covariances are present (MF32), they
   ! are processed and combined with the MF33 covariances.  This
   ! coding was taken from ERRORJ developed in Japan.
   !
   ! Covariances for angular distributions (MF34) and secondary
   ! energy distributions (MF35) can also be processed.
   !
   ! 9/17/2012 ... important changes
   ! - ngout = input gendf input tape
   !           -  this tape may contain multiple temperatures,
   !              multiple Legendre moments and multiple sigma-0
   !              data, but ...
   !              - only the first (infinitely dilute) sigma-0
   !                data will be used.
   !
   ! - card 3 is now REQUIRED and the specified temperature must
   !   be one of those on the gendf tape.
   !   - default values are mprint=1, tempin=300
   !
   !---input specifications (free format)---------------------------
   !
   !  card 1
   !    nendf   unit for endf tape
   !    npend   unit for pendf tape
   !    ngout   unit for input group xsec (gendf) tape
   !            (if zero, group xsecs will be calculated)
   !            (if iread eq 2 or if mfcov eq 31, 35 or 40 (see
   !             card 7), then ngout cannot be zero)
   !            (if mfcov eq 35 (see card 7),
   !              ngout cannot be zero)
   !            (default=0)
   !    nout    unit for output covariance tape (default=0)
   !    nin     unit for input covariance tape (default=0)
   !            (nin and nout must be both coded or both binary)
   !    nstan   unit for ratio-to-standard tape (default=0)
   !
   !  card 2
   !    matd    material to be processed
   !    ign     neutron group option
   !            (ign definition same as groupr, except ign=-1,
   !            which means read in an energy grid, as in ign=1,
   !            and supplement this with the endf covariance grid
   !            within the range of the user-specified energies)
   !            (default=1)
   !    iwt     weight function option (default=6)
   !    iprint  print option (0/1=minimum/maximum) (default=1)
   !    irelco  covariance form (0/1=absolute/relative) (default=1)
   !
   !  card 3    *** REQIURED for njoy2012_0917 & later ***
   !    mprint  print option for group averaging (0=min., 1=max.)
   !    tempin  temperature (default=300)
   !
   !---for endf/b version 4 (iverf=4) only--------------------------
   !
   !  card 4
   !    nek     number of derived xsec energy ranges
   !            (if zero, all xsecs are independent)
   !  card 5    (omit if nek=0)
   !    ek      nek+1 derived xsec energy bounds
   !  card 6    (omit if nek=0)
   !    akxy    derived cross section coefficients, one row/line
   !
   !---for endf/b version 5 or 6 (iverf=5 or 6) only------------------------
   !
   !  card 7
   !    iread   0/1/2=program calculated mts/input mts and eks/
   !            calculated mts plus extra mat1-mt1 pairs from input
   !            (default=0)
   !    mfcov   endf covariance file (31, 33, 34, 35 or 40) to be
   !            processed (default=33).
   !            note--contribution to group cross section
   !            covariances from resonance-parameter uncertainties
   !            (mf=32) is included when mfcov=33 is specified.
   !            (mf=-33) high speed calc for test case
   !            (mf=333) hight speed calc for test case (faster)
   !    irespr  processing for resonance parameter covariances
   !            (mf=32) (default=1)
   !            0 = area sensitivity method
   !            1 = 1% sensitivity method
   !    legord  legendre order for calculating covariances (default=1)
   !            (if mfcov is not 34, legord is ignored)
   !    ifissp  subsection of the fission spectrum covariance
   !            matrix to process (default=-1 which means process
   !            the subsection that includes efmean).  The value
   !            for ifissp that appears in njoy's standard output
   !            will equal the subsection containing efmean.
   !            (if mfcov is not 35, ifissp is ignored)
   !    efmean  incident neutron energy (eV).  Process the covar-
   !            iance matrix subsection whose energy interval in-
   !            cludes efmean.  if ifissp=-1 and efmean is not
   !            specified, set efmean=2.e6_kr.  But if there is only
   !            one subsection, process it and if the input efmean
   !            was not within this subsection's energy range then
   !            redefine efmean to equal the average energy for
   !            this subsection.
   !            (if mfcov is not 35, efmean is ignored)
   !    dap     user specified scattering radius uncertainty, given
   !            as a fraction (i.e., dap=0.1 means 10% uncertainty
   !            in the scattering radius).  The default value is
   !            zero.  This variable is only defined for mfcov=33
   !            and if non-zero will be used in lieu of any data
   !            that might have been read from the nendf tape.
   !
   !  following cards only if iread eq 1
   !  card 8
   !    nmt     no. mts to be processed
   !    nek     no. derived cross section energy ranges
   !            (if zero, all xsecs are independent)
   !  card 8a
   !    mts     nmt mts
   !  card 8b   (omit if nek=0)
   !    ek      nek+1 derived cross section energy bounds
   !  card 9    (omit if nek=0)
   !    akxy    derived cross section coefficients, one row/line
   !
   !  following card only if iread eq 2
   !  card 10
   !    mat1    cross-material reaction to be added to
   !    mt1         covariance reaction list.
   !            repeat for all mat1-mt1 pairs desired
   !            terminate with mat1=0.
   !
   !  following card only if nstan ne 0
   !  card 11
   !    matb    standards reaction referenced
   !    mtb         in matd.
   !    matc    standards reaction to be
   !    mtc         used instead.
   !            repeat for all standard reactions to be redefined.
   !            terminate with matb=0.
   !  note.  if matb(1) and mtb(1) are negative, then matc(1) and
   !    mtc(1) identify a third reaction, correlated with matd thru
   !    the use of the same standard.  covariances of all reactions
   !    in matd (which reference the standard) with the reaction
   !    matc(1)-mtc(1) will be produced.  the standard reaction
   !    must be identified on card 10 and repeated as the negative
   !    entries on card 11.  the group xsec tape ngout must include
   !    all covariance reactions in matd, plus matc(1)-mtc(1).
   !
   !  card 12a (for ign eq 1 or ign eq -1)
   !    ngn     number of groups
   !  card 12b
   !    egn     ngn+1 group bounds (ev)
   !  card 13 (for iwt eq 1 only)
   !    wght    weight function as a tab1 record
   !  card 13b  analytic flux parameters (iwt=4 only)
   !    eb      thermal break (ev)
   !    tb      thermal temperature (ev)
   !    ec      fission break (ev)
   !    tc      fission temperature (ev)
   !
   ! ---options for input variables------------------------------------
   !
   !      ign          meaning
   !      ---          -------
   !       1           arbitrary structure (read in)
   !      -1           arbitrary structure (read in), supplemented with endf
   !                   covariance grid
   !       2           csewg 239-group structure
   !       3           lanl 30-group structure
   !       4           anl 27-group structure
   !       5           rrd 50-group structure
   !       6           gam-i 68-group structure
   !       7           gam-ii 100-group structure
   !       8           laser-thermos 35-group structure
   !       9           epri-cpm 69-group structure
   !      10           lanl 187-group structure
   !      11           lanl 70-group structure
   !      12           sand-ii 620-group structure
   !      13           lanl 80-group structure
   !      14           eurlib 100-group structure
   !      15           sand-iia 640-group structure
   !      16           vitamin-e 174-group structure
   !      17           vitamin-j 175-group structure
   !      18           xmas nea-lanl
   !      19           ecco  33-group structure
   !      20           ecco 1968-group structure
   !      21           tripoli 315-group structure
   !      22           xmas lwpc 172-group structure
   !      23           vit-j lwpc 175-group structure
   !      24           shem cea 281-group structure
   !      25           shem epm 295-group structure
   !      26           shem cea/epm 361-group structure
   !      27           shem epm 315-group structure
   !      28           rahab aecl 89-group structure
   !      29           ccfe   660-group structure  (30 MeV)
   !      30           ukaea 1025-group structure  (30 MeV)
   !      31           ukaea 1067-group structure (200 MeV)
   !      32           ukaea 1102-group structure   (1 GeV)
   !      33           ukaea  142-group structure (200 MeV)   
   !
   !      iwt          meaning
   !      ---          -------
   !       1           read in smooth weight function
   !       2           constant
   !       3           1/e
   !       4           1/e + fission spectrum + thermal maxwellian
   !       5           epri-cell lwr
   !       6           (thermal) -- (1/e) -- (fission + fusion)
   !       7           same with t-dep thermal part
   !       8           thermal--1/e--fast reactor--fission + fusion
   !       9           claw weight function
   !      10           claw with t-dependent thermal part
   !      11           vitamin-e weight function (ornl-5505)
   !      12           vit-e with t-dep thermal part
   !
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso,nsyse
   use util ! provides timer,openz,error,repoz,closz,sigfig
   use endf ! provides endf routines and variables
   use samm ! provides s2sammy
   use physics ! provides amassn, amu, ev, hbar
   ! internals
   integer::i,icov,nb,nw,mprint,nwi,nx,ndictm,nwl,ix,l
   integer::mf,mfi,mt,neki,j,k,ii1,ii2,ii3,ii4
   integer::lim,ii,ntape,nek1,idone,iadd,nwscr
   integer::ng,ngp,iwtt,iw,np,nr
   integer::mat,nfissp,is
   integer::iaddmt(5)
   real(kr)::dmy,time
   character(60)::strng
   real(kr)::b(17),z(10)
   character(4)::id(17)
   real(kr)::rd(17)
   equivalence(id(1),rd(1))
   real(kr),dimension(:),allocatable::dict
   real(kr),dimension(:),allocatable::temp
   real(kr),dimension(:),allocatable::scr
   integer,parameter::maxa=50000
   real(kr),dimension(:),allocatable::a
   real(kr),parameter::eps=1.e-9_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::elo=1.e-5_kr
   real(kr),parameter::zero=0
   logical::need1,need2

   !--set samrml options
   Want_Partial_Derivs=.true.
   Want_Angular_Dist=.false.
   Want_SAMRML_BW=.false.
   Want_SAMRML_RM=.false.

   cwaven=sqrt(2*amassn*amu*ev)*1.e-12_kr/hbar

   !--read user input and write header.
   imode=-1
   iread=0
   mfcov=33
   iwt=0
   nas=0
   eclo=1.0e-5_kr
   echi=2.0e7_kr
   efmean=2.0e6_kr
   mffis=0
   igflag=0
   ncove=0
   ncovl=0
   call timer(time)
   write(nsyso,'(/&
     &'' errorr...produce cross section covariances'',&
     &26x,f8.1,''s'')') time
   write(nsyse,'(/'' errorr...'',59x,f8.1,''s'')') time
   ngout=0
   nout=0
   nin=0
   nstan=0
   nscr5=0
   read(nsysi,*) nendf,npend,ngout,nout,nin,nstan
   ngoutu=ngout

   if (nendf.eq.999) then
      read(nsysi,*) nitape,notape
      iadd=0
 1000 continue
      read(nsysi,*) ii
      if (ii.ne.0) then
         iadd=iadd+1
         if (iadd.gt.5) then
            call error('errorr','errorr in 999 option',' ')
         endif
         iaddmt(iadd)=ii
         goto 1000
      else
         call covadd(iadd,iaddmt,5,nitape,notape)
         return
      endif
      call error('errorr','errorr in 999 option',' ')
   endif

   if (nendf.eq.nstan)&
      call error('errorr','nstan should be different from nendf','')

   call openz(nendf,0)
   call openz(npend,0)
   call openz(ngout,0)
   call openz(nout,1)
   call openz(nin,0)
   call openz(nstan,0)
   call repoz(nendf)
   call tpidio(nendf,0,0,b,nb,nw)
   call contio(nendf,0,0,b,nb,nw)
   call contio(nendf,0,0,b,nb,nw)
   emaxx=2.0e7_kr
   if (n1h.ne.0.and.n2h.eq.0) then
      iverf=4
      emaxx=1.5e7_kr
   else if (n2h.eq.0) then
      iverf=5
   else
      iverf=6
      call contio(nendf,0,0,b,nb,nw)
      if (c2h.gt.zero) emaxx=c2h
   endif
   nmt=0
   nmt1=0
   ign=1
   iwt=6
   iprint=1
   irelco=1
   mwtf=0
   read(nsysi,*) matd,ign,iwt,iprint,irelco
   if (iwt.le.0) then
      call mess ('errorj','input weighting function not supported',&
                 'switching to default, iwt=6')
      iwt=6
   endif
   if (ign.eq.19) then
      call mess ('errorr','neutron group structure option 19 has been ',&
                 'changed, you may want to use -1 instead')
      iwt=6
   endif
   mprint=0
   tempin=300
   read(nsysi,*) mprint,tempin   ! ***Card 3 ... now required***
   allocate(eni(nenimx))
   write(nsyso,'(/&
     &'' unit for endf tape ................... '',i10/&
     &'' unit for pendf tape .................. '',i10/&
     &'' unit for input gendf tape ............ '',i10/&
     &'' unit for output covariance tape ...... '',i10/&
     &'' unit for input covariance tape ....... '',i10/&
     &'' unit for ratio-to-standard tape ...... '',i10)')&
     nendf,npend,ngoutu,nout,nin,nstan
   write(nsyso,'(&
     &'' material to be processed ............. '',i10/&
     &'' neutron group option ................. '',i10/&
     &'' weight function option ............... '',i10/&
     &'' print option (0 min, 1 max) .......... '',i10/&
     &'' rel. cov. option (0 abs, 1 rel) ...... '',i10)')&
     matd,ign,iwt,iprint,irelco
   write(nsyso,'(&
     &'' group av. print option (0 min, 1 max)  '',i10)') mprint
   if (tempin.eq.0.) write(nsyso,'(&
     &'' temperature .......................... '',6x,''zero'')')
   if (tempin.ne.0.) write(nsyso,'(&
     &'' temperature .......................... '',f10.1)') tempin
   if (iverf.ne.4) then
      iread=0
      mfcov=33
      irespr=1
      legord=1
      ifissp=-1
      isru=0
      dap=0
      read(nsysi,*) iread,mfcov,irespr,legord,ifissp,efmean,dap

      !--only allow legord=1 at this time
      if (legord.ne.1) then
         write(strng,'(''reset legord from '',i2,'' to 1'')')legord
         call mess('errorr',strng,'')
         legord=1
      endif

      !--input check for legal and/or consistent options
      if (iread.lt.0.or.iread.gt.2) then
         write(strng,'(''illegal iread='',i3)') iread
         call error('errorr',strng,' ')
      endif

      if (mfcov.eq.33.and.dap.ne.0) isru=1

      if ((mfcov.eq.31.or.mfcov.eq.35.or.mfcov.eq.40).and.ngout.eq.0)&
      then
         write(strng,'('' ngout must be nonzero when mfcov ='',i2)')&
               mfcov
         call error('errorr',strng,' ')
      endif

      !--for mf35 processing, check for spectrum data on the
      !--gendf tape, check for covariance matrix data on
      !--the endf tape and make a copy of this mat on nendf2.
      if (mfcov.eq.35) then
         call ngchk
         ngout=nscr5  !make ngout point to condensed gendf tape
         nwscr=500000
         allocate(scr(nwscr))
         ! check for spectrum in the gendf input file
         mat=0
         mf=0
         do while (mffis.le.6.and.mf.lt.7.and.mat.ne.-1)
            if (ngout.gt.0) then
               read(ngout,'(66x,i4,i2,i3)') mat,mf,mt
            else
               read(-ngout) mat,mf,mt
            endif
            if ((mf.eq.5.or.mf.eq.6).and.mt.eq.18) then
               mffis=mf
            endif
         enddo
         call repoz(ngout)
         if (mffis.eq.0) then
            strng='mfcov=35 requested but no spectra present'
            call error('errorr',strng,' ')
         else if(mffis.le.5) then
           strng='mf35 processing with the mf5 spectrum vector'
           call mess('errorr',strng,' ')
         else
           strng='mf35 processing with the mf6 spectrum matrix'
           call mess('errorr',strng,' ')
         endif
         ! check for the presence of spectrum covariances and for
         ! legal ifissp and/or efmean input.  if ifissp>0, set
         ! efmean=average for this interval; if ifissp=-1,
         ! redefine ifissp to whatever interval includes efmean.
         call findf(matd,35,18,nendf)
         call contio(nendf,0,0,scr,nb,nw)
         nfissp=n1h
         if (ifissp.gt.nfissp) then
            write(strng,'(''User ifissp ('',i2,'') not found.  Max on '',&
               & ''file is '',i2,''.'')') ifissp,nfissp
            call error('errorr',strng,' ')
         endif
         i=0
         idone=0
         do while (i.lt.nfissp)
            i=i+1
            call listio(nendf,0,0,scr,nb,nw)
            if (idone.eq.0) then
               eclo=c1h
               echi=c2h
            endif
            if (n1h.gt.ncovl) ncovl=n1h+6
            is=1
            do while (nb.ne.0)
               is=is+nw
               if (is.gt.nwscr) call error('errorj','storage exceeded.',' ')
               call moreio(nendf,0,0,scr(is),nb,nw)
            enddo
            if (i.eq.ifissp) then
               efmean=(eclo+echi)/2
               idone=1
            else if (ifissp.le.0) then
               if (eclo.le.efmean .and. echi.ge.efmean) then
                  ifissp=i
                  igflag=1
                  idone=1
               endif
               if (nfissp.eq.1) then
                  if (efmean.lt.eclo.or.efmean.gt.echi) then
                     efmean=(eclo+echi)/2
                     ifissp=1
                     write(strng,'('' reset efmean to '',1pe11.4,&
                        & '' eV.'')')efmean
                     call mess('errorr','only one subsection',strng)
                  else
                     ifissp=1
                     igflag=1
                  endif
               endif
            endif
         enddo
         if (ifissp.le.0) call error('errorj',&
           'no covariance data found for user ifissp/efmean',' ')
         call repoz(nendf)
         ! copy matd from nendf to nendf2
         nendf2=19
         if (nendf.lt.0)nendf2=-nendf2
         call openz(nendf2,0)
         call tpidio(nendf,0,0,b,nb,nw)
         math=0
         mfh=0
         mth=0
         call tpidio(0,0,nendf2,b,nb,nw)
         math=0
         do while (math.le.matd)
            call contio(nendf,0,0,b,nb,nw)
            if (math.eq.matd) then
               call contio(0,0,nendf2,b,nb,nw)
               call tomend(nendf,0,nendf2,scr)
               call atend(0,nendf2)
               math=10000
            else
               call tomend(nendf,0,0,b)
            endif
         enddo
         if (math.ne.10000) call error('errorj',&
           'nendf-to-nendf2 copy for matd failed',' ')
         deallocate(scr)
         call repoz(nendf)
      endif

      !--high speed test options
      eskip1=1.00002e0_kr
      eskip2=1.0003e0_kr
      eskip3=1.005e0_kr
      eskip4=1.05e0_kr
      if (mfcov.eq.-33) then
         eskip1=1.0002e0_kr
         eskip2=1.0003e0_kr
         eskip3=1.005e0_kr
         mfcov=33
      endif
      if (mfcov.eq.333) then
         eskip1=1.002e0_kr
         eskip2=1.003e0_kr
         eskip3=1.005e0_kr
         mfcov=33
      endif

      write(nsyso,'(&
        &'' read option (0 calc, 1 read, 2 combo)  '',i10/&
        &'' endf covariance file to be processed . '',i10)')&
        iread,mfcov
      if (mfcov.eq.33) write(nsyso,&
        '('' irespr processing option for mf=33 ... '',i10)') irespr
      if (isru.ne.0) write(nsyso,&
        '('' user dap for scat radius unc override '',f11.4)') dap
      if (mfcov.eq.34) write(nsyso,&
        '('' legendre order for mf=34 ............. '',i10)') legord
      if (mfcov.eq.35) then
         write(nsyso,&
           '('' igflag ............................... '',i10/&
           &'' ifissp ............................... '',i10/&
           &'' efmean ...............................'',1p,e11.4)')&
           igflag,ifissp,efmean
         write(nsyso,'('' covariance matrix energy range .......'',&
           &1p,e11.4,'' to'',e11.4,'' eV'')') eclo,echi
      endif

      if (mfcov.ne.31.and.mfcov.ne.33.and.mfcov.ne.34.and.&
          mfcov.ne.35.and.mfcov.ne.40) then
         write(strng,'(''not coded for mfcov='',i3)') mfcov
         call error('errorr',strng,' ')
      endif
      write(nsyso,'(/'' using endf-'',i1,'' format'')') iverf
   endif

   !--read covariance reaction types from end/b dictionary
   !--and set file 32 flag
   nscr2=0
   nwi=4000
   allocate(dict(nwi))
   call repoz(nendf)
   call tpidio(nendf,0,0,dict,nb,nw)
   call findf(matd,1,451,nendf)
   call contio(nendf,0,0,dict,nb,nw)
   nx=nint(dict(6))
   if (iverf.gt.4) call contio(nendf,0,0,dict,nb,nw)
   if (iverf.gt.5) call contio(nendf,0,0,dict,nb,nw)
   call hdatio(nendf,0,0,dict,nb,nw)
   if (iverf.gt.4) nx=nint(dict(6))
   do i=1,17
      rd(i)=dict(i+6)
   enddo
   do while (nb.ne.0)
      call moreio(nendf,0,0,dict,nb,nw)
   enddo
   ndictm=6*nx
   if (ndictm.gt.nwi) call error('errorr',&
     'not enough space for endf dictionary',' ')
   nw=nx
   call dictio(nendf,0,0,dict,nb,nw)
   nmt=0
   mf32=0
   mf33=0
   nga=0
   nwi=200
   if (ngout.eq.0.or.mfcov.eq.34) allocate(iga(nwi))
   nlump=0
   nwl=nlmt*2
   allocate(lump(nwl))
   icov=0
   do ix=1,nx
      l=6*ix-3
      mf=nint(dict(l))
      if (mf.eq.32) mf32=1
      if (mf.eq.33) mf33=1
      mt=nint(dict(l+1))

      !--check dictionary for required files
      if (mfcov.eq.30.and.(mf.ge.30 .and. mf.le.33)) icov=icov+1
      if (mfcov.eq.31.and.mf.eq.31) icov=icov+1
      if (mfcov.eq.33.and.(mf.eq.32.or.mf.eq.33)) icov=icov+1
      if (mfcov.eq.34.and.mf.eq.34) icov=icov+1
      if (mfcov.eq.35.and.mf.eq.35) icov=icov+1
      if (mfcov.eq.40.and.mf.eq.40) icov=icov+1
      if (mf.eq.mfcov) then
         if (mt.eq.850.or.(mt.gt.870.and.mt.lt.875).or.mt.gt.891) then
            write(strng,'(''ignoring unknown mt = '',i5)')mt
            call mess('errorr',strng,'')
            cycle
         endif
         if (mt.gt.850.and.mt.le.870) go to 121
         if (ngout.ne.0.and.mfcov.ne.34) go to 125
         nga=nga+1
         if (nga.gt.nwi) call error('errorr','too many reactions for mf34.',' ')
         iga(nga)=mt
         go to 125
     121 continue
         nlump=nlump+1
         if (nlump.gt.nlmt)&
           call error('errorr','too many lumped reaction types',' ')
         lump(1+2*(nlump-1))=mt
         lump(1+2*(nlump-1)+1)=0
     125 continue
         if (iverf.le.4) then
            nmt=nmt+1
            if (nmt.gt.nmtmax)&
              call error('errorr','too many reaction types.',' ')
            nmt1=nmt
            mats(nmt)=0
            mts(nmt)=mt
         endif
      endif
   enddo
   call tofend(nendf,0,0,dict)
   deallocate(dict)
   nwl=nlump*2

   !--check MF2 for the use of the sammy method
   nmtres=0
   if (isammy.gt.0) then
      allocate(a(maxa))
      call s2sammy(nendf,a,maxa,mmtres,nmtres)
      if (nmtres.gt.0) then
         do i=1,nmtres
            if (mmtres(i).eq.103) mmtres(i)=600
            if (mmtres(i).eq.104) mmtres(i)=650
            if (mmtres(i).eq.105) mmtres(i)=700
            if (mmtres(i).eq.106) mmtres(i)=750
            if (mmtres(i).eq.107) mmtres(i)=800
         enddo
      endif
      call findf(matd,2,0,nendf)
      deallocate(a)
   endif

   if (iverf.gt.4) go to 200

   !--set up coefficients for derived cross sections
   !--for covariance data in the endf/b-iv format
   read(nsysi,*) neki
   if (neki.gt.nkmax) then
      write(strng,'(''only'',i3,'' ek energies allowed'')') nkmax
      call error('errorr',strng,' ')
   endif
   write(nsyso,'(&
     &'' no. of derived xsec energy ranges .... '',i10)') neki
   nek=neki
   if (neki.eq.0) nek=1
   allocate(akxy(nmt,nmt,nek))
   if (neki.eq.0) then
      ek(1)=small
      ek(2)=big
      do i=1,nmt
         do j=1,nmt
            akxy(j,i,1)=0
            if (i.eq.j) akxy(j,i,1)=1
         enddo
      enddo
   else
     nek1=nek+1
     read(nsysi,*) (ek(i),i=1,neki)
     do i=1,nek1
        ek(i)=sigfig(ek(i),ndig,0)
     enddo
     do i=1,nek
        do j=1,nmt
           nw=nmt
           do k=1,nw
              z(k)=0
           enddo
             read(nsysi,*) (z(k),k=1,nw)
           do k=1,nmt
              akxy(k,j,i)=z(k)
           enddo
        enddo
     enddo
   endif
   go to 280

   !--set up coefficients for derived cross sections
   !--for covariance data in the endf/b-v, vi, and vii formats
  200 continue
   if (iread.eq.1) then
      nek=0
      read(nsysi,*) nmt,nek
      if (nmt.gt.nmtmax)&
        call error('errorr','too many reaction types.',' ')
      nmt1=nmt
      neki=nek
      if (nek.eq.0) nek=1
      nek1=nek+1
      write(nsyso,'(&
        &'' no. of mts to be processed ........... '',i10)') nmt
      write(nsyso,'(&
        &'' no. of derived xsec energy ranges .... '',i10)') nek
      allocate(akxy(nmt,nmt,nek))
      nw=nmtmax
      allocate(temp(nw))
      read(nsysi,*) (temp(i),i=1,nmt)
      do i=1,nmt
         mats(i)=0
         mts(i)=nint(temp(i))
      enddo
      ek(1)=small
      ek(2)=big
      if (neki.eq.0) then
         do j=1,nmt
            do k=1,nmt
               akxy(k,j,1)=0
               if (j.eq.k) akxy(k,j,1)=1
            enddo
         enddo
      else
         read(nsysi,*) (ek(i),i=1,nek1)
         do i=1,nek1
            ek(i)=sigfig(ek(i),ndig,0)
         enddo
         do i=1,nek
            do j=1,nmt
               read(nsysi,*) (akxy(k,j,i),k=1,nmt)
            enddo
         enddo
      endif
      deallocate(temp)
   endif

   !--read additional user-supplied mat1-mt1 pairs
   if (iread.eq.2) then
      idone=0
      do while (idone.eq.0)
         ii1=0
         read(nsysi,*) ii1,ii2
         if (ii1.eq.0) idone=1
         if (idone.eq.0) then
            nmt1=nmt1+1
            mats(nmt1)=ii1
            mts(nmt1)=ii2
         endif
      enddo
   endif

   !--read input for redefining the standard
   if (nstan.ne.0) then
      idone=0
      do while (idone.eq.0)
         ii1=0
         read(nsysi,*) ii1,ii2,ii3,ii4
         if (ii1.eq.0) idone=1
         if (idone.eq.0) then
            nas=nas+1
            if (nas.gt.nasmax)&
              call error('errorr','too many standards redefined.',' ')
            matb(nas)=ii1
            mtb(nas)=ii2
            matc(nas)=ii3
            mtc(nas)=ii4
         endif
      enddo
   endif

   !--check if there are data on file to process
   if (icov.eq.0) then
      write(strng,'(''no data on file for mfcov='',i3)') mfcov
      call mess('errorr',strng,'processing terminated')
      !--skip remaining errorr input (if any)
      if (abs(ign).eq.1) then
         read(nsysi,*) ng
         ngp=ng+1
         read(nsysi,*) (dmy,i=1,ngp)
      endif
      iwtt=iabs(iwt)
      if (iwtt.eq.1) then
         iw=1000
         allocate(wght(iw))
         read(nsysi,*)(wght(i),i=1,iw)
         nr=nint(wght(5))
         np=nint(wght(6))
         iw=6+2*nr+2*np
         deallocate(wght)
      else if (iwtt.eq.4) then
         read(nsysi,*) eb,tb,ec,tc
      endif
      go to 330
   endif

   !--if present, condense the input GENDF tape:
   !  - will only contain matd plus (if iread=2) the nmt1 list of
   !    materials from card 10 at tempin, infinitely dilute and p1.
   !  - if mfcov=35, we already did this.
   if (ngout.ne.0 .and. mfcov.ne.35) then
      call ngchk
      ngout=nscr5  !make ngout point to condensed gendf tape
   endif

   !--construct the union energy grid and akxy array
   call gridd(neki)
   if (nlump.gt.0) call lumpmt

   !--print the akxy array
  280 continue
   if (neki.ne.0) then
      write(nsyso,'(/'' coefficients for derived cross sections'')')
      lim=nmt
      if (nmt.gt.10) lim=10
      do i=1,nek
         write(nsyso,'(/'' for'',1p,e12.4,'' to '',e12.4,'' ev''/)')&
           ek(i),ek(i+1)
         write(nsyso,'(3x,''mt -'',10(4x,i3))') (mts(ii),ii=1,lim)
         if (nmt.gt.lim) write(nsyso,'(7x,10i7)') (mts(ii),ii=11,nmt)
         write(nsyso,'(''  --- -'',10(4x,a3))') ('---',ii=1,lim)
         do j=1,nmt
            write(nsyso,'(2x,i3,'' -'',10f7.1/(5x,'' -'',10f7.1))')&
              mts(j),(akxy(k,j,i),k=1,nmt)
         enddo
      enddo
   endif

   !--print out resonance method used
   if (mf33.eq.1.and.nmtres.gt.0) then
      write(nsyso,'(/'' using sammy resonance method'')')
   else if (mf33.eq.1) then
      write(nsyso,'(/'' using errorj resonance method'')')
   endif

   !--prepare a scratch tape for use by grpav or colaps & define
   !  the user multigroup array (egn(1:ngn+1))

   ntape=-10 !careful ... must match ngout (grpav) and ntp (colaps)!
   call openz(ntape,1)
   call egngpn

   !--if an mfcov=34 job, check various mf4 and mf34 flags to see
   !  whether a gendf tape is required.
   if (mfcov.eq.34) then
      need1=.false.           !Assumed default, no input gendf tape
      need2=.false.           !Assumed default, no gendf tape required
      if (ngout.ne.0) need1=.true. !User has supplied a gendf tape
      if (ltt4.eq.0.or.ltt4.eq.2) need2=.true.     !gendf tape required
      if (ltt4.eq.3.and.&
          (eltt4+elo).lt.egn(ngn+1)) need2=.true.  !gendf tape required
      if (need2 .and. .not.need1) call error('errorr',&
                                 'ngout input tape is required','')
   endif

   !--compute group constants on union grid, either from pointwise
   !  input (npend) or multigroup input (ngout).
   !  - net result is an ngout=-10 "gendf" tape with multigroup data
   !    collapsed to the union energy grid.
   if (ngout.eq.0) then
      call grpav(mprint,tempin)
   else
      call colaps
   endif

   !--if processing mf=34 data (mt2/mt251, or mu-bar only at the moment)
   !  then call grpav4 and create a gendf-like tape with a mf4, mt2
   !  header record followed by union grid multigroup mu-bar identified
   !  as mf3,mt2.
   if (mfcov.eq.34) then
      call grpav4(mprint)
   endif

   write(nsyso,'(/&
     &'' processing mat  '',i4/&
     &'' ---------------------''/&
     &1x,17a4)') matd,(id(i),i=1,17)
   write(nsyse,'(/&
     &'' processing mat  '',i4/&
     &'' ---------------------''/&
     &1x,17a4)') matd,(id(i),i=1,17)

   !--compute MF33 covariance matrices
   ek(1)=sigfig(ek(1),ndig,0)
   if (abs(egn(1)-elo).le.eps) egn(1)=elo
   call covcal

   !--add MF32 covariance matrices and write output tape.
   call covout

   !--errorr is finished.
   call atend(nout,0)
  330 continue
   if (allocated(flx))  deallocate(flx)
   if (allocated(sig))  deallocate(sig)
   if (allocated(cov))  deallocate(cov)
   if (allocated(egn))  deallocate(egn)
   if (allocated(lump)) deallocate(lump)
   if (allocated(iga))  deallocate(iga)
   if (allocated(un))   deallocate(un)
   if (allocated(akxy)) deallocate(akxy)
   if (allocated(eni))  deallocate(eni)
   if (allocated(wght)) deallocate(wght)
   call repoz(nout)
   call repoz(nin)
   call repoz(nendf)
   call repoz(ngout)
   call repoz(npend)
   call repoz(nstan)
   call closz(nstan)
   call closz(nendf)
   if (mfcov.eq.35) call closz(nendf2)
   call closz(npend)
   call closz(ngout)
   call closz(nin)
   call closz(nout)
   call timer(time)
   write(nsyso,'(69x,f8.1,''s''/1x,77(''*''))') time
   return
   end subroutine errorr

   subroutine gridd(neki)
   !--------------------------------------------------------------------
   ! Read through mfcov in version 5 or 6 format and extract the union
   ! energy grid for the derivation relations (NC-type sub-subsections
   ! with LTY=0), construct the matrix of derivation coefficients
   ! and extract the union energy grid from NI-type sub-subsections
   ! for fine-group covariance calculations.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso,nsyse
   use util ! provides repoz,error,sigfig
   use endf ! provides endf routines and variables
   ! externals
   integer::neki
   ! internals
   integer::nb,nw,nwscr,ir,nmtt,nsub,i,isub,izap
   integer::mat1,mt1,iok,nmtp,nc,ni,ic,lty,l,nt,matstd,mtstd
   integer::jtop,j,ii,nx,lb,nec,iloc,nl
   integer::ik,nr,irs,ilo,ihi,ntr,nim,ixm,ider,nder
   integer::jakj,mat2,ijk,ld1,ld,mt2,nl1
   integer::ilfs,nfs
   real(kr)::elh,ehh
   character(60)::strng1,strng2
   integer,parameter::irmax=80
   real(kr)::el(irmax),eh(irmax)
   integer::nmtr(irmax),imtr(irmax),jak(irmax)
   integer,parameter::nxmax=5000
   real(kr)::x(nxmax)
   real(kr),dimension(:),allocatable::scr
   real(kr),dimension(:,:),allocatable::ak
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::small5=1.e-5_kr
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::zero=0

   !--allocate storage and initialize.
   if (iread.ne.1) then
      allocate(ak(nmtmax*2,irmax))
   endif
   nwscr=1400000
   allocate(scr(nwscr))
   neki=0
   neni=0
   ir=0
   nmtt=nmt1
   call repoz(nendf)
   call tpidio(nendf,0,0,scr,nb,nw)

   if (mfcov.eq.34) then
   !--get ltt, lct and (maybe) transition energy, eltt4, between
   !  Legendre and probability distribution representations.
      call findf(matd,4,2,nendf)
      call contio(nendf,0,0,scr,nb,nw)
      ltt4=l2h
      call contio(nendf,0,0,scr,nb,nw)
      lct4=l2h
      eltt4=emaxx
      if (ltt4.eq.3) then
         call tab2io(nendf,0,0,scr,nb,nw)
         ii=n2h
         do i=1,ii
            call listio(nendf,0,0,scr,nb,nw)
         enddo
         eltt4=c2h
      endif
   endif

   call findf(matd,mfcov,0,nendf)
   ld=0
   ld1=0
   mat2=0
   mt2=0

   !--loop over sections.
  110 continue
   call contio(nendf,0,0,scr,nb,nw)
   if (mfh.eq.0.and.iread.eq.1) then
      deallocate(scr)
      go to 610
   endif
   if (mfcov.eq.34 .and. mth.ne.2) then
     call tofend(nendf,0,0,scr)
     go to 450
   elseif (mfcov.eq.34) then
     ltt34=l2h
   endif
   if (mfh.eq.0) go to 450
   ! ignore components of a lumped reaction
   izap=0
   nfs=1
   if (mfcov.eq.35) then
      nsub=n1h
   else
      nsub=n2h  !NMT1 in the manual
      if (mfcov.eq.40) then
         za=c1h
         nfs=n1h
         nsub=1
         call contio(nendf,0,0,scr,nb,nw)
         lfs=l2h
         izap=10*l1h+l2h
      endif
   endif
   if (nsub.eq.0) go to 410
   if (iread.ne.1) go to 130
   do 120 i=1,nmt
      if (mth.eq.mts(i)) go to 160
  120 continue
   go to 420
  130 continue
   nmt=nmt+1
   nmt1=nmt1+1
   if (nmt1.gt.nmtmax)&
     call error('gridd','too many reaction types.',' ')
   if (iread.eq.2) then
      do i=1,nmtt
         mats(nmt+nmtt+1-i)=mats(nmt+nmtt-i)
         mts(nmt+nmtt+1-i)=mts(nmt+nmtt-i)
      enddo
   endif
   mats(nmt)=0
   mts(nmt)=mth
   mzap(nmt)=izap
  160 continue

   !--angular distribution covariances
   if (mfcov.eq.34) then
      if (nsub.gt.1) call error('gridd','not coded for nmt1>1 of mf=34',' ')
      call contio(nendf,0,0,scr,nb,nw)
      mat2=l1h  !MAT1 in the manual
      mt2=l2h   !MT1  in the manual
      nl=n1h
      nl1=n2h
      nsub=nl*nl1
      if (mth.eq.mt2) nsub=nl*(nl+1)/2
   endif

   !--loop over subsections
   do 310 ilfs=1,nfs
   do 300 isub=1,nsub
   if (mfcov.ne.35) call contio(nendf,0,0,scr,nb,nw)
   if (mfcov.eq.34.and.mth.eq.0) go to 110
   iok=1
   if (mfcov.eq.34) then
      mat1=mat2
      mt1=mt2
      ld=l1h    !L1  in the manual
      ld1=l2h   !L11 in the manual
      lct34=n1h
      ni=n2h
      nc=0
      go to 290
   else if (mfcov.eq.35) then
      mat1=0
      mt1=mth
   else
      mat1=l1h
      mt1=l2h
   endif
   if (mt1.eq.0) call error('gridd','illegal mt1=0.',' ')
   if (iread.gt.0) go to 161
   if (mat1.gt.0) go to 175
   go to 180
  161 continue
   if (iread.gt.1) go to 165
   if (mat1.gt.0) go to 175
   do 162 i=1,nmt
   if (mt1.eq.mts(i)) go to 180
  162 continue
   go to 175
  165 continue
   if (mat1.eq.0) go to 180
   nmtp=nmt+1
   do 170 i=nmtp,nmt1
   if (mat1.eq.mats(i).and.mt1.eq.mts(i)) go to 180
  170 continue
   ! covariance matrix for mat1-mt1 is present in mfcov , but is
   ! not wanted by user.  flag this case by setting iok=0, in order
   ! to avoid adding unnecessary points to the union energy grid.
  175 continue
   iok=0
  180 continue
   if (mfcov.eq.35) then
      nc=0
      ni=1
   else
      nc=n1h
      ni=n2h
   endif
   if (nc.eq.0) go to 290

   !--loop over nc sub-subsections
   do 210 ic=1,nc
      call contio(nendf,0,0,scr,nb,nw)
      lty=l2h
      if (lty.gt.3) then
         write(nsyso,'(/'' ***error in gridd***covariances of reaction mt='',&
           &i3,'' with (mat1='',i4,'' mt1='',i3,'')''/&
           &4x,''cannot be calculated.  not coded for lty='',i3)')&
           mth,mat1,mt1,lty
         write(nsyse,'(/'' ***error in gridd***covariances of reaction mt='',&
           &i3,'' with (mat1='',i4,'' mt1='',i3,'')''/&
           &4x,''cannot be calculated.  not coded for lty='',i3)')&
           mth,mat1,mt1,lty
         call error('gridd',' ',' ')
      endif
      call listio(nendf,0,0,scr,nb,nw)
      l=1
      do while (nb.ne.0)
         l=l+nw
         if (l.gt.nwscr) call error('gridd','nc subsection too big',&
           'see nwscr')
         call moreio(nendf,0,0,scr(l),nb,nw)
      enddo
      continue
      if (iok.eq.0) go to 210
      if (mfcov.eq.34) then
         if (ld.gt.legord.or.ld1.gt.legord) go to 210
      endif
      nt=2
      call merge(scr(1),nt,nxmax,eni,neni,nenimx,ndig,zero,zero)
      elh=sigfig(c1h,ndig,0)
      ehh=sigfig(c2h,ndig,0)
      if (lty.eq.0) go to 260
      matstd=l1h
      mtstd=l2h
      if (nstan.eq.0) then
         write(nsyso,'(/&
           &'' ***error in gridd***cannot calculate covariances of '',&
           &''reaction mt='',i3,'' with'',/4x,''mt1='',i3,&
           &'' because nstan=0.  to proceed, mount an endf tape '',&
           &''containing'',/4x,''the standard'',&
           &'' reaction (matstd=,i4,8h, mtstd=,i3,16h) on unit nstan.''/&
           &4x,''if necessary matstd and mtstd can be redefined on input '',&
           &''card 11.'')') mth,mt1,matstd,mtstd
         write(nsyse,'(/&
           &'' ***error in gridd***cannot calculate covariances of '',&
           &''reaction mt='',i3,'' with'',/4x,''mt1='',i3,&
           &'' because nstan=0.  to proceed, mount an endf tape '',&
           &''containing'',/4x,''the standard'',&
           &'' reaction (matstd=,i4,8h, mtstd=,i3,16h) on unit nstan.''/&
           &4x,''if necessary matstd and mtstd can be redefined on input '',&
           &''card 11.'')') mth,mt1,matstd,mtstd
         call error('gridd',' ',' ')
      endif
      if (lty.eq.1) call grist(matstd,mtstd,nxmax,elh,ehh,scr,x)
      if (lty.eq.2) call grist(matstd,mtstd,nxmax,zero,zero,scr,x)
      go to 210
  260 continue
      if (iread.eq.1) go to 210
      ir=ir+1
      if (ir.gt.irmax) call error('gridd',&
        'too many formulas in nc-type sub-subsections with lty=0.',' ')
      imtr(ir)=nmt
      el(ir)=elh
      eh(ir)=ehh
      nmtr(ir)=n2h
      if (n2h.gt.nmtmax) call error('gridd',&
        'too many mt-numbers in nc-type subsections with lty=0.',' ')

      !--save the derivation formula for later processing.
      jtop=n1h
      do j=1,jtop
         ak(j,ir)=scr(j+6)
      enddo
      nt=2
      call merge(scr(1),nt,nxmax,ek,neki,nkmax,ndig,zero,zero)
  210 continue
  290 continue
   if (ni.eq.0) go to 300

   !--loop over ni sub-subsections
   do 350 ii=1,ni
      call listio(nendf,0,0,scr,nb,nw)
      if (mfh.eq.35) then
         if (n2h.gt.ncove) ncove=n2h
      endif
      l=1
      do while (nb.ne.0)
         l=l+nw
         if (l.gt.nwscr) call error('gridd','ni subsection too big',&
           'see nwscr')
         call moreio(nendf,0,0,scr(l),nb,nw)
      enddo
      if (iok.eq.0) go to 350
      if (mfcov.eq.34) then
         if (ld.gt.legord.or.ld1.gt.legord) go to 350
      endif
      nx=n2h
      if (nx.gt.nxmax) call error('gridd','nx is too large',&
            'increase nxmax')
      lb=l2h
      if (lb.lt.5.or.lb.eq.8) then
         do i=1,nx
            x(i)=scr(2*i+5)
         enddo
         nl=l1h
         nx=nx-nl
         call merge(x,nx,nxmax,eni,neni,nenimx,ndig,zero,zero)
         call merge(x(1+nx),nl,nxmax,eni,neni,nenimx,ndig,zero,zero)
      else
         call merge(scr(7),nx,nxmax,eni,neni,nenimx,ndig,zero,zero)
         if (lb.eq.5) go to 350
         if (mfcov.eq.35.and.lb.eq.7) go to 350
         nec=(n1h-1)/nx
         iloc=7+nx
         call merge(scr(iloc),nec,nxmax,eni,neni,nenimx,ndig,zero,zero)
      endif
  350 continue
  300 continue
  310 continue
  410 continue
   call contio(nendf,0,0,scr,nb,nw)
   go to 110
  420 continue
   call tosend(nendf,0,0,scr)
   go to 110

   !--set up coefficients for derived cross sections.
  450 continue
   call repoz(nendf)
   call tpidio(nendf,0,0,scr,nb,nw)
   deallocate(scr)
   if (neki.gt.0.and.small5+small5/1000.lt.ek(1)) then
      ! if necessary, increment neki so that a default coefficient
      ! table is produced for the energy range 10**-5 eV to ek(1),
      ! plus redefine the ek(neki) array to include this interval.
      neki=neki+1
      do ijk=neki,2,-1
         ek(ijk)=ek(ijk-1)
      enddo
      ek(1)=small5
   endif
   nek=neki-1
   if (neki.eq.0) nek=1
   allocate(akxy(nmt1,nmt1,nek))
   do ik=1,nek
      do i=1,nmt1
         do j=1,nmt1
            akxy(j,i,ik)=0
            if (i.eq.j) akxy(j,i,ik)=1
         enddo
      enddo
   enddo
   if (neki.gt.0) go to 500
   ek(1)=small
   ek(2)=big
   go to 600

   !--reconstruct full akxy table.
  500 continue
   nr=ir
   do 590 ir=1,nr
      irs=ir
      ilo=0
      do
         ilo=ilo+1
         if (el(ir).eq.ek(ilo)) exit
      enddo
      ihi=ilo
      do
         ihi=ihi+1
         if (eh(ir).eq.ek(ihi)) exit
      enddo
      ihi=ihi-1
      ntr=nmtr(ir)
      do 560 nim=1,ntr
         ixm=nint(ak(2*nim,ir))
         do 550 i=1,nmt
            jak(nim)=i
            if (ixm.eq.mts(i)) go to 560
  550    continue
         write(strng1,&
           '(''mt'',i3,'' referenced in derivation formula'')') ixm
         write(strng2,&
           '(''for range '',i2,'' does not appear in mfcov'')') irs
         call error('gridd',strng1,strng2)
  560 continue
      ider=imtr(ir)
      nder=nmtr(ir)
      do i=ilo,ihi
         do j=1,nder
            jakj=jak(j)
            akxy(jakj,ider,i)=ak(2*j-1,ir)
         enddo
         akxy(ider,ider,i)=0
      enddo
  590 continue
  600 continue
   deallocate(ak)
  610 continue
   if (allocated(scr)) deallocate(scr)
   return
   end subroutine gridd

   subroutine merge(x,nx,nxmax,y,ny,nymax,ndig,e1,e2)
   !--------------------------------------------------------------------
   ! Merge an energy grid in x with a previously existing one in y.
   ! Both grids are assumed to be in increasing order.
   !--------------------------------------------------------------------
   use util ! provides error,sigfig
   ! externals
   integer::nx,nxmax,ny,nymax,ndig
   real(kr)::x(*),y(*),e1,e2
   ! internals
   integer::j,i,jsave,isave,k,loc,ny1
   character(60)::strng
   real(kr),parameter::eps=1.e-5_kr
   real(kr),parameter::zero=0

   if (nx.eq.0) return
   j=0
   do 100 i=1,nx
   x(i)=sigfig(x(i),ndig,0)
   if (e1.eq.zero.and.e2.eq.zero) go to 110
   if (x(i).le.e1) go to 100
   if (x(i).ge.e2) go to 120
  110 continue
   j=j+1
   x(j)=x(i)
  100 continue
  120 continue
   nx=j
   if (ny.gt.0) go to 140
   if (nx.gt.nxmax) call error('merge','x storage exceeded.',' ')
   do i=1,nx
      y(i)=x(i)
   enddo
   ny=nx
   go to 200

  140 continue
   j=0
   do 180 i=1,nx
  150 continue
   j=j+1
   if (j.gt.ny) go to 170
   if (y(j).lt.x(i)) go to 150
   if (y(j).eq.x(i)) go to 180
   if (abs(y(j)-x(i)).le.eps*y(j)) go to 180
   ! insert x(i) in the y array.
   if (ny.eq.nymax) call error('merge','y storage exceeded.',' ')
   do k=j,ny
      loc=ny+j-k
      y(loc+1)=y(loc)
   enddo
  170 continue
   ny=ny+1
   y(j)=x(i)
  180 continue

   !--check grid for data inconsistencies
  200 continue
   ny1=ny-1
   do i=1,ny1
      isave=i
      do j=i,ny
         jsave=j
         if (y(i).gt.y(j)) then
            write(strng,&
              '(''y('',i4,'')='',1p,e12.4,'' lt y('',i4,'')='',&
              &1p,e12.4)') jsave,y(jsave),isave,y(isave)
            call error('merge',strng,' ')
         endif
      enddo
   enddo
   return
   end subroutine merge

   subroutine grist(matstd,mtstd,nxmax,el,eh,scr,x)
   !--------------------------------------------------------------------
   ! Merge the energy grid from the standard tape nstan into the union
   ! energy grid from nendf.
   !--------------------------------------------------------------------
   use util ! provides repoz,mess,error,sigfig
   use endf ! provides endf routines and variables
   ! externals
   integer::matstd,mtstd,nxmax
   real(kr)::el,eh,scr(*),x(*)
   ! internals
   integer::nb,nw,i,nsub,il,mat1,mt1,nc,ni,ic,lty,nt,ii
   integer::l,lb,nx,nec,iloc,nl,idone
   real(kr)::zero
   character(60)::strng
   real(kr),parameter::small=1.e-10_kr

   call repoz(nstan)
   call tpidio(nstan,0,0,scr,nb,nw)

   !--redefine standard if necessary
   if (nas.gt.0) then
      i=0
      idone=0
      do while (i.lt.nas.and.idone.eq.0)
         i=i+1
         if (matstd.eq.matb(i).and.mtstd.eq.mtb(i)) then
            write(strng,'(''standards reaction (,'',i4,'','',i3,&
              &'') replaced by ('',i4,'','',i3,'')'')')&
              matstd,mtstd,matc(i),mtc(i)
            call mess('grist',strng,' ')
            matstd=matc(i)
            mtstd=mtc(i)
            idone=1
         endif
      enddo
   endif
   call findf(matstd,mfcov,mtstd,nstan)
   call contio(nstan,0,0,scr,nb,nw)
   nsub=n2h
   if (nsub.le.0) call error('grist','standards tape bad.',' ')

   !--loop over subsections
   il=0
   idone=0
   do while (il.lt.nsub.and.idone.eq.0)
      il=il+1
      call contio(nstan,0,0,scr,nb,nw)
      mat1=l1h
      mt1=l2h
      nc=n1h
      ni=n2h

      !--read and merge energies from the nc-type sub-subsections
      if (nc.gt.0) then
         do ic=1,nc
            call contio(nstan,0,0,scr,nb,nw)
            lty=l2h
            call listio(nstan,0,0,scr,nb,nw)
            do while (nb.ne.0)
               call moreio(nstan,0,0,scr(3),nb,nw)
            enddo
            if (il.ne.1) then
               if (mat1.eq.matc(1).and.mt1.eq.mtc(1).and.lty.eq.3) then
                  ety1=sigfig(c1h,ndig,0)
                  ety2=sigfig(c2h,ndig,0)
                  nt=2
                  zero=0
                  call merge(scr,nt,nxmax,eni,neni,nenimx,ndig,&
                    zero,zero)
               endif
            endif
         enddo
      endif

      !--read and merge the energies from the ni sub-subsections
      if (ni.eq.0.and.il.eq.1) then
         write(strng,'(''matstd='',i4,'', mtstd='',i3)') matstd,mtstd
         call error('grist','illegal ni=0 in the standard',strng)
      endif
      if (ni.ne.0) then
         do ii=1,ni
            call listio(nstan,0,0,scr,nb,nw)
            l=1
            do while (nb.ne.0)
               l=l+nw
               call moreio(nstan,0,0,scr(l),nb,nw)
            enddo
            continue
            if (il.le.1) then
               lb=l2h
               nx=n2h
               if (lb.eq.0) call error('grist','illegal lb=0.',' ')
               if (lb.ge.5.and.lb.ne.8) then
                  call merge(scr(7),nx,nxmax,eni,neni,nenimx,ndig,el,eh)
                  if (lb.ne.5) then
                     nec=(n1h-1)/nx
                     iloc=7+nx
                     call merge(scr(iloc),nec,nxmax,eni,neni,nenimx,&
                       ndig,el,eh)
                  endif
               else
                  do i=1,nx
                     x(i)=scr(2*i+5)
                  enddo
                  nl=l1h
                  nx=nx-nl
                  call merge(x(1),nx,nxmax,eni,neni,nenimx,ndig,el,eh)
                  call merge(x(1+nx),nl,nxmax,eni,neni,nenimx,ndig,&
                    el,eh)
               endif
            endif
         enddo
         if (el.gt.small) idone=1
         if (nas.eq.0) idone=1
         if (matb(1).ge.0) idone=1
         if (matstd.ne.-matb(1).or.mtstd.ne.-mtb(1)) idone=1
      endif
   enddo
   return
   end subroutine grist

   subroutine lumpmt
   !--------------------------------------------------------------------
   ! Read through the File 33 MTs and store the list of component MTs
   ! making up the lumped MTs.
   !--------------------------------------------------------------------
   use util ! provides error,openz,repoz
   use endf ! provides endf routines and variables
   ! internals
   integer::nb,nw,i,nwl,max,mt1,l,mtl,k,j,loc1,loc
   character(4)::cid(17),bl='    '
   real(kr)::rid(17)
   real(kr),dimension(:),allocatable::scr
   equivalence (cid(1),rid(1))

   nwl=nlump*nlmt
   allocate(lmt(nwl))
   nw=2*npage+50
   allocate(scr(nw))
   do i=1,nwl
      lmt(i)=0
   enddo
   call findf(matd,mfcov,0,nendf)
   max=0

   !--loop over mts
  110 continue
   call contio(nendf,0,0,scr,nb,nw)
   if (mfh.eq.0) go to 200
   mt1=l2h
   if (mt1.lt.851.or.mt1.gt.870) go to 140
   do 120 l=1,nlump
   mtl=lump(1+2*(l-1))
   if (mt1.ne.mtl) go to 120
   lump(1+2*(l-1)+1)=lump(1+2*(l-1)+1)+1
   k=lump(1+2*(l-1)+1)
   if (k.gt.nlmt) call error('lumpmt','storage exceeded.',' ')
   if (k.gt.max) max=k
   lmt(nlmt*(l-1)+k)=mth
   ! set this mth in mts negative
   do 130 j=1,nmt
   if (mts(j).ne.mth) go to 130
   mts(j)=-mts(j)
   if (mats(j).eq.0) mats(j)=-1
   if (mats(j).gt.0) mats(j)=-mats(j)
   go to 140
  130 continue
  120 continue
  140 continue
   call tosend(nendf,0,0,scr)
   go to 110

   !--determine the maximum no. of words needed
  200 continue
   if (max.ne.nlmt) then
      ! squeeze storage
      loc1=max
      do l=2,nlump
         loc=nlmt*(l-1)
         do j=1,max
            lmt(j+loc1)=lmt(j+loc)
         enddo
         loc1=loc1+max
      enddo
      nwl=max*nlump
      nlmt=max
   endif
   ! copy mfcov to nscr3 for use in lumpxs
   nscr3=15
   if (nendf.lt.0) nscr3=-nscr3
   call openz(nscr3,1)
   call repoz(nscr3)
   do i=1,17
      read(bl,'(a4)') cid(i)
   enddo
   math=0
   mfh=0
   mth=0
   nsc=0
   call tpidio(0,0,nscr3,rid,nb,nw)
   call findf(matd,mfcov,0,nendf)
   call contio(nendf,0,nscr3,scr,nb,nw)
   call tofend(nendf,0,nscr3,scr)
   call amend(0,nscr3)
   call atend(0,nscr3)
   call repoz(nscr3)
   deallocate(scr)
   return
   end subroutine lumpmt

   subroutine covcal
   !--------------------------------------------------------------------
   ! Calculate absolute covariances in the union-group structure.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso,nsyse
   use util ! provides openz,repoz,error,mess,timer
   use endf ! provides endf routines and variables
   ! internals
   integer::nmts,nb,nw,nwi,mfd,mtd,is,nl,i,mtl,il,mat1,mt1
   integer::nc,ni,iok,kmt1,li,l,lty,ic,ltyi,np,locli,jg
   integer::jh,loci,lt,lb,locip4,locip6,nk1,k,k2,locnec,nl1
   integer::ifloc,nlt,nk,nlt1,locl,lend,loclp4,loclp6
   integer::l2,m,m2,jgend,ih,ibase,ij,kmtb,isrrr,namx
   integer::mat2,mt2,nlg1,nlg2,ld,ld1,izap,izero
   real(kr)::eg,xcv,time,de,flux,dne
   character(70)::strng,strng2
   integer,parameter::locm=30
   integer::loc(locm)
   real(kr),dimension(:),allocatable::sig1,scr2,scr
   real(kr),dimension(:),allocatable::alp1,alp2
   real(kr),dimension(:),allocatable::egt
   real(kr),dimension(:),allocatable::b
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::zero=0

   !--initialize
   izero=0
   nscr2=13
   if (ngout.lt.0) nscr2=-nscr2
   nscr1=11*imode
   call openz(nscr1,1)
   call openz(nscr2,1)
   call repoz(nscr1)
   nmts=0
   za=0
   nwi=nunion+1
   nw=nwi
   if ((npage+50).gt.nw) nw=npage+50
   allocate(b(nw))
   allocate(egt(nwi))
   mfd=1
   mtd=451
   kmt1=0
   mat2=0
   mt2=0
   nl=0
   ! skip over energy group bounds on ngout
   call rdgout(ngout,matd,mfd,mtd,izero,b,egt)
   if (mtd.gt.nwi)&
     call error('covcal','storage exceeded in egt.',' ')
   deallocate(egt)
   ! assign storage.
   allocate(flx(nunion))
   allocate(sig(nunion+1))
   allocate(cov(nunion))
   allocate(sig1(nunion+1))
   allocate(scr2(nunion+1))
   if (mfcov.eq.34) then
      allocate(alp1(nunion))
      allocate(alp2(nunion))
   endif
   if (mfcov.eq.35) then
      namx=max(npage+6,ncovl)
   else
      namx=2000000
   endif
   allocate(scr(namx))
   mfd=3
   mtd=-1
   ! store flux array for later use
   call rdgout(ngout,matd,mfd,mtd,izero,b,flx)
   nsc=0
   call rdsig(matd,izero,izero,b,scr)

   !--temporary patch for missing mf33
   if (mf33.eq.0.and.mfcov.eq.33) then
      go to 1000
   endif

   !--if the total cross section is absent, estimate the
   !--fine-group fluxes in the sub-threshold energy region
   !--by assuming dn/de is constant
   if (mtd.eq.-2) then
      dne=small
      do is=1,nunion
         isrrr=nunion+1-is
         flux=flx(isrrr)
         de=un(1+isrrr)-un(isrrr)
         if (flux.ne.zero.and.de.ne.zero) then
            dne=flux/de
         else
            flx(isrrr)=dne*de
         endif
      enddo
      mtd=-1
   endif
   kmtb=0

   !--loop over reactions in mfcov
   call findf(matd,mfcov,0,nendf)
  140 continue
   call contio(nendf,0,0,scr,nb,nw)
   if (math.lt.1) go to 700
   if (math.ne.matd) go to 700
   if (mth.eq.0) go to 140
   if (mfh.lt.mfcov) go to 140
   if (mfh.gt.mfcov) go to 700
   if (mfcov.eq.34 .and. mth.ne.2) go to 700
   if (iverf.eq.4) then
      nl=l2h
   else if (iverf.gt.4) then
      nl=n2h
   endif
   izap=0
   if (mfcov.eq.34) then
      nl=1
   else if (mfcov.eq.35) then
      nl=n1h
   else if (mfcov.eq.40) then
      za=c1h
      awr=c2h
      nl=n1h
      call contio(nendf,0,0,scr,nb,nw)
      if (l1h.eq.0) then
         write(strng,'(''WARNING!  izap=0 for mf40, mt'',i3,&
                      &'' is nonstandard'')')mth
         if (l2h.eq.0) then
            write(strng2,'(''if produced, covr plot title will '',&
                          &''be ambiguous'')')
            call mess('covcal',strng,strng2)
         else
            write(strng2,'(''covcal may skip this mt, or covr '',&
                          &''plot title will be ambiguous'')')
            call mess('covcal',strng,strng2)
         endif
      endif
      izap=10*l1h+l2h
   endif

   ! ignore components of a lumped reaction
   if (nl.eq.0) go to 140

   if (iread.ne.1) go to 170
   do 150 i=1,nmt
   if (mth.eq.mts(i)) go to 170
  150 continue
   call tosend(nendf,0,0,scr)
   go to 140
  170 continue
   if (za.eq.zero) then
      za=c1h
      awr=c2h
   endif
   nmts=nmts+1
   if (mts(nmts).ne.mth) call error('covcal',&
     'mfcov mt found not equal to input mt.',' ')
   if (mfcov.eq.34) then
      call contio(nendf,0,0,scr,nb,nw)
      mat2=l1h
      mt2=l2h
      nlg1=min(n1h,legord)
      nlg2=min(n2h,legord)
      nl=nlg1*(nlg2+1)/2
      b(1)=za
      b(2)=awr
      b(3)=nlg1
      b(4)=nlg2
      b(5)=nl
      b(6)=nunion
      nl=n1h*(n2h+1)/2
   else
      b(1)=za
      b(2)=awr
      b(3)=0
      b(4)=nl
      b(5)=0
      b(6)=nunion
   endif
   call contio(0,0,nscr1,b,nb,nw)
   mtl=mth
   if (mth.le.850.or.mth.gt.870) then
      if (mfcov.eq.35) then
         call rdchi(math,b,sig)
      else
         call rdsig(math,mth,izap,b,sig)
      endif
   else
      call lumpxs(mtl,mtl,sig,sig1,b,scr2)
   endif

   !--loop over different covariance matrices for this reaction
   do 650 il=1,nl
   if (mfcov.eq.35) then
      mat1=0
      mt1=mth
      nc=0
      ni=1
   else
      call contio(nendf,0,0,scr,nb,nw)
      if (mfcov.eq.34.and.mth.eq.0) go to 660
      if (mfcov.eq.34) then
         mat1=mat2
         mt1=mt2
         ld=l1h
         ld1=l2h
      else
         mat1=l1h
         mt1=l2h
      endif
      if (mt1.eq.0) call error('covcal','illegal mt1=0.',' ')
      nc=n1h
      ni=n2h
   endif
   if (ni.gt.locm) call error('covcal','storage exceeded in loc.',' ')
   iok=1
   do 210 i=1,nmt1
   kmt1=i
   if (mt1.eq.mts(i).and.mat1.eq.mats(i)) go to 220
  210 continue
   ! covariance matrix for mat1-mt1 is present in mfcov, but is
   ! not wanted by user.  flag this by setting iok=0, and later
   ! write a null matrix on the output file.
   iok=0

   !--if necessary, redefine mat1 and mt1
  220 continue
   if (nas.eq.0) go to 230
   if (mat1.ne.-matb(1).or.mt1.ne.-mtb(1)) go to 230
   if (iok.eq.0) then
      write(strng,&
        '(''must request mat1='',i3,'' and mt1='',i3)') mat1,mt1
      call error('covcal',strng,'on card 10.')
   endif
   mat1=matc(1)
   mt1=mtc(1)
   kmtb=kmt1

   !--read in all sub-subsections for this matrix.
  230 continue
   li=0
   l=1
   if (nc.gt.0) then
      lty=0
      do ic=1,nc
         if (iverf.gt.4) call contio(nendf,0,0,scr(l),nb,nw)
         if (iverf.gt.4) lty=l2h
         call listio(nendf,0,0,scr(l),nb,nw)
         do while (nb.ne.0)
            call moreio(nendf,0,0,scr(l),nb,nw)
         enddo
         if (iok.ne.0) then
            if (mfcov.eq.34) then
               if (ld.gt.legord.or.ld1.gt.legord) cycle
            endif
            if (lty.gt.0.and.lty.lt.4) call stand(li,l,loc,lty,scr)
         endif
      enddo
   endif
   if (ni.eq.0.and.li.eq.0) go to 600
   if (ni.gt.0) go to 285
   ni=li
   go to 320
  285 continue
   ltyi=0
   ni=ni+li
  290 continue
   li=li+1
   loc(li)=l
   call listio(nendf,0,0,scr(l),nb,nw)
   np=n1h
   if (l2h.eq.6) scr(l+2)=(n1h-1)/n2h
   scr(l+4)=ltyi
   l=l+nw
   do while (nb.ne.0)
      if (l.gt.namx) call error('covcal',&
        'storage exceeded in scr.','see namx.')
      call moreio(nendf,0,0,scr(l),nb,nw)
      l=l+nw
   enddo
   locli=loc(li)+5
   if (mfcov.eq.35) then
      call sumchk(scr(loc(li)))
      call covbin(scr(loc(li)))
   endif
   do i=1,np
      scr(i+locli)=sigfig(scr(i+locli),ndig,0)
   enddo
   if (li.lt.ni) go to 290
  320 continue
   if (iok.eq.0) go to 600
   if (mfcov.eq.34) then
      if (ld.gt.legord.or.ld1.gt.legord) go to 650
   endif

   !--retrieve sigma for mt1, either from ngout or sig.
   if (kmt1.ne.nmts) then
      if (mt1.lt.851.or.mt1.gt.870) then
         call rdsig(mat1,mt1,izero,b,sig1)
      else
         call lumpxs(mt1,mtl,sig,sig1,b,scr2)
      endif
   else
      do jg=1,nunion
         sig1(jg)=sig(jg)
      enddo
   endif

   if (mfcov.eq.34) then
      call rdlgnd(nscr4,matd,mth,ld,b,alp1)
      call rdlgnd(nscr4,matd,mt1,ld1,b,alp2)
   endif

   !--generate covariance matrix using specified laws.
   do 570 jg=1,nunion
   eg=un(jg)
   do 520 jh=1,nunion
   ehr=un(jh)
   cov(jh)=0
   do 510 i=1,ni
   loci=loc(i)
   lt=nint(scr(loci+2))
   lb=nint(scr(loci+3))
   ltyi=nint(scr(loci+4))
   np=nint(scr(loci+5))
   if (mfcov.eq.34.and.&
     (lb.lt.0.or.lb.eq.3.or.lb.eq.4.or.lb.eq.7.or.lb.gt.8)) then
      write(strng,'(''unpermitted for lb='',i2)') lb
      call error('covcal',strng,'in mf=34.')
   else if (mfcov.eq.35.and.lb.ne.7) then
      write(strng,'("unpermitted for lb=",i2)') lb
      call error('covcal',strng,'in mf=35.')
   endif
   if (ltyi.eq.0) go to 345
   if (ltyi.lt.1.or.ltyi.gt.3.or.&
     scr(loci).le.0..or.scr(loci+1).le.scr(loci))&
     call error('covcal','data in scr(loci) is illegal.',' ')
   if (ltyi.eq.3) go to 340
   ! for lty = 1 and 2, apply energy window to mt groups
   if (eg.lt.scr(loci).or.eg.ge.scr(loci+1)) go to 510
   if (ltyi.eq.2) go to 345
   ! for lty = 1 and 3, apply energy window to mt1 groups
  340 continue
   if (ehr.lt.scr(loci).or.ehr.ge.scr(loci+1)) go to 510
   go to 346
   ! if necessary, apply ety energy window to mt1 groups
  345 continue
   if (nas.eq.0) go to 346
   if (matb(1).ge.0) go to 346
   if (mat1.ne.matc(1).or.mt1.ne.mtc(1)) go to 346
   if (ehr.lt.ety1.or.ehr.ge.ety2) go to 510
  346 continue
   if (lb.eq.7.or.lb.gt.8) then
      if (mfcov.eq.35.and.lb.eq.7) go to 347
      write(strng,'(''not coded for lb='',i2)') lb
      call error('covcal',strng,' ')
   endif
   if (lb.lt.3.and.lt.gt.0) then
      write(strng,'(''lb='',i2,'' when lt='',i2)') lb,lt
      call error('covcal',strng,' ')
   endif
  347 continue
   locip4=loci+4
   locip6=loci+6
   nk1=np-1
   if (lb.ne.8) go to 880

   !--separate treatment for lb=8
   if (jh.ne.jg) go to 510
   k=0
  860 continue
   k=k+1
   k2=k*2
   if (eg.ge.scr(locip4+k2).and.eg.lt.scr(locip6+k2)) go to 870
   if (k.lt.nk1) go to 860
   go to 510
   ! assume flux is constant within a union group
  870 continue
   xcv=(scr(locip6+k2)-scr(locip4+k2))/(un(1+jg)-un(jg))
   cov(jh)=cov(jh)+scr(loci+5+k2)*xcv
   go to 510
  880 if (lb.ne.7) go to 800

   !--separate treatment for lb=7 (mf=35)
   k=0
  890 continue
   k=k+1
   if (eg.ge.scr(locip6+k-1).and.eg.lt.scr(locip6+k)) go to 900
   if (k.lt.nk1) go to 890
   go to 510
  900 continue
   l=0
  910 continue
   l=l+1
   if (ehr.ge.scr(locip6+l-1).and.ehr.lt.scr(locip6+l)) go to 920
   if (l.lt.nk1) go to 910
   go to 510
  920 continue
   if (l.ge.k) then
      ifloc=locip6+nk1*np/2-(np-k+1)*(np-k)/2+l-k
   else
      ifloc=locip6+nk1*np/2-(np-l+1)*(np-l)/2+k-l
   endif
   cov(jh)=cov(jh)+scr(ifloc+np)
   go to 510
  800 continue
   if (lb.ne.6) go to 850

   !--separate treatment for lb=6
   k=0
  810 continue
   k=k+1
   if (eg.ge.scr(locip6+k-1).and.eg.lt.scr(locip6+k)) go to 820
   if (k.lt.nk1) go to 810
   go to 510
  820 continue
   locnec=locip6+np
   nl1=lt-1
   l=0
  830 continue
   l=l+1
   if (ehr.ge.scr(locnec+l-1).and.ehr.lt.scr(locnec+l)) go to 840
   if (l.lt.nl1) go to 830
   go to 510
  840 continue
   ifloc=locnec+lt+(k-1)*nl1+l-1
   if (mfcov.eq.34) then
      cov(jh)=cov(jh)+scr(ifloc)*alp1(jg)*alp2(jh)
   else
      cov(jh)=cov(jh)+scr(ifloc)*sig(jg)*sig1(jh)
   endif
   go to 510
  850 continue
   if (lb.ne.5) go to 410

   !--separate treatment for lb=5.
   k=0
  350 continue
   k=k+1
   if (eg.ge.scr(locip6+k-1).and.eg.lt.scr(locip6+k)) go to 360
   if (k.lt.nk1) go to 350
   go to 510
  360 continue
   l=0
  370 continue
   l=l+1
   if (ehr.ge.scr(locip6+l-1).and.ehr.lt.scr(locip6+l)) go to 380
   if (l.lt.nk1) go to 370
   go to 510
  380 continue
   if (lt.eq.1) go to 390
   ifloc=locip6+(k-1)*nk1+l-1
   go to 400
  390 continue
   ifloc=locip6+nk1*np/2-(np-l+1)*(np-l)/2+k-l
   if (l.ge.k) ifloc=locip6+nk1*np/2-(np-k+1)*(np-k)/2+l-k
  400 continue
   if (mfcov.eq.34) then
      cov(jh)=cov(jh)+scr(ifloc+np)*alp1(jg)*alp2(jh)
   else
      cov(jh)=cov(jh)+scr(ifloc+np)*sig(jg)*sig1(jh)
   endif
   go to 510

   !--integrated treatment for lb=0 thru lb=4.
  410 continue
   nlt=lt
   nk=np-nlt
   nk1=nk-1
   nlt1=nlt-1
   k=0
  420 continue
   k=k+1
   k2=k*2
   if (eg.lt.scr(locip4+k2).or.eg.ge.scr(locip6+k2)) go to 430
   if (lb.eq.2.or.lb.eq.3) go to 440
   if (ehr.ge.scr(locip4+k2).and.ehr.lt.scr(locip6+k2)) go to 490
   go to 510
  430 continue
   if (k.lt.nk1) go to 420
   go to 510
  440 continue
   if (lb.gt.2) go to 450
   locl=loci
   lend=nk1
   go to 460
  450 continue
   locl=loci+nk*2
   lend=nlt1
  460 continue
   loclp4=locl+4
   loclp6=locl+6
   l=0
  470 continue
   l=l+1
   l2=l*2
   if (ehr.ge.scr(loclp4+l2).and.ehr.lt.scr(loclp6+l2)) go to 480
   if (l.lt.lend) go to 470
   go to 510
  480 continue
   if (lb.ne.4) go to 486
   m=0
  482 continue
   m=m+1
   m2=m*2
   if (eg.ge.scr(loclp4+m2).and.ehr.lt.scr(loclp6+m2)) go to 484
   if (m.lt.lend) go to 482
   go to 510
  484 continue
   if (mfcov.eq.34) then
      cov(jh)=cov(jh)+scr(loci+5+k2)*scr(locl+5+m2)*scr(locl+5+l2)&
        *alp1(jg)*alp2(jh)
   else
      cov(jh)=cov(jh)+scr(loci+5+k2)*scr(locl+5+m2)*scr(locl+5+l2)&
        *sig(jg)*sig1(jh)
   endif
   go to 510
  486 continue
   if (mfcov.eq.34) then
      cov(jh)=cov(jh)+scr(loci+5+k2)*scr(locl+5+l2)*alp1(jg)*alp2(jh)
   else
      cov(jh)=cov(jh)+scr(loci+5+k2)*scr(locl+5+l2)*sig(jg)*sig1(jh)
   endif
   go to 510
  490 continue
   if (lb.eq.4) go to 450
   if (lb.eq.0) go to 500
   if (mfcov.eq.34) then
      cov(jh)=cov(jh)+scr(loci+5+k2)*alp1(jg)*alp2(jh)
   else
      cov(jh)=cov(jh)+scr(loci+5+k2)*sig(jg)*sig1(jh)
   endif
   go to 510
  500 continue
   cov(jh)=cov(jh)+scr(loci+5+k2)
  510 continue
  520 continue

   !--write one row of the covariance matrix on scratch tape.
   jgend=0
   do ih=1,nunion
      if (cov(ih).ne.zero) jgend=ih
   enddo
   if (jgend.gt.0) go to 540
   if (jg.lt.nunion) go to 570
   jgend=1
  540 continue
   mfh=mfcov
   math=matd
   mth=mts(nmts)
   if (mfcov.eq.34) then
      b(1)=ld
      b(2)=ld1
   else
      b(1)=0
      b(2)=0
   endif
   b(3)=mat1
   b(4)=mt1
   b(5)=jgend
   b(6)=jg
   ibase=6
   ic=ibase
   do ij=1,jgend
      ic=ic+1
      if (mfcov.eq.34) then
         b(ic)=cov(ij)*(sig(jg)*flx(jg))*(sig1(ij)*flx(ij))
      else
         b(ic)=cov(ij)*flx(jg)*flx(ij)
      endif
      if ((ic-ibase).ge.npage.or.ij.eq.jgend) then
         if (ibase.ne.0) then
            call listio(0,0,nscr1,b,nb,ic)
            ibase=0
            ic=0
         else
            call moreio(0,0,nscr1,b,nb,ic)
            ic=0
         endif
      endif
   enddo
  570 continue
   go to 650

   !--write a null covariance matrix on scratch tape.
  600 continue
   mth=mts(nmts)
   math=matd
   mfh=mfcov
   b(1)=0
   b(2)=0
   b(3)=mat1
   b(4)=mt1
   b(5)=1
   b(6)=nunion
   b(7)=0
   nw=7
   call listio(0,0,nscr1,b,nb,nw)

   !--close loop over subsections of mfcov
  650 continue
  660 continue
   call asend(0,nscr1)

   !--close loop over sections of mfcov
   go to 140
  700 continue
   ! if necessary, redefine one mats(i)-mts(i) pair
   if (kmtb.gt.0) then
      mats(kmtb)=matc(1)
      mts(kmtb)=mtc(1)
   endif

   !--covcal is finished
  1000 continue
   call afend(0,nscr1)
   call amend(0,nscr1)
   call atend(0,nscr1)
   call timer(time)
   if (mf33.gt.0 .or. mfcov.ne.33) then
      write(nsyso,'(/&
        &'' covariances calculated for '',i2,'' reactions and '',&
        &i4,'' groups'',13x,f8.1,''s'')') nmts,nunion,time
      write(nsyse,'(/&
        &'' covariances calculated for '',i2,'' reactions and '',&
        &i4,'' groups'',13x,f8.1,''s'')') nmts,nunion,time
   endif
   deallocate(sig1)
   deallocate(scr2)
   deallocate(scr)
   if (mfcov.eq.34) then
      deallocate(alp1)
      deallocate(alp2)
   endif
   deallocate(b)
   ! call closz(nscr2)  !leave open for covout.

   return
   end subroutine covcal

   subroutine sumchk(covl)
   !--------------------------------------------------------------------
   ! compute row (or column) sums for the input covariance matrix and
   ! the multigroup spectrum integral using the covariance matrix
   ! energy grid.  Then perform the endf specified "zero-sum" rule
   ! test.  If necessary, correct the individual matrix elements per
   ! the endf manual (section 35.3 in the april, 2001 edition).
   !  covl(1:ncovl)  = mf35 list record (header, matrix energy grid
   !                   and the triangular covariance matrix elements).
   !--------------------------------------------------------------------
   use util ! provides error,mess
   ! externals
   real(kr)::covl(*)
   ! internals
   integer::isym,lb,ibase,ne,ne1,i,j,is
   real(kr)::sumt
   character(66)::c
   real(kr)::b(17)
   real(kr),dimension(:),allocatable::rcs,spc
   real(kr),parameter::stst=1.0e-5_kr
   real(kr),parameter::sml=1.0e-30_kr

   !--check for required symmetry and lb=7 flags
   isym=nint(covl(3))
   lb=nint(covl(4))
   if (isym.ne.1) call error('sumchk','endf file error',&
                             'mf35 matrix must be symmetric')
   if (lb.ne.7) call error('sumchk','endf file error',&
                           'mf35 matrix lb flag must be 7')
   ibase=6
   ne=nint(covl(ibase))
   if (ne.gt.ncove) call error('sumchk','ne, ncove mismatch',' ')
   ne1=ne-1

   !--reserve space for ne1 row (column) sums
   !--and ne1 spectrum integrals.
   allocate(rcs(ne1))
   allocate(spc(ne1))

   !--set pointer to matrix elements (a triangular matrix given in
   !--a vector array starting at cov(ibase+1).  also initialize the
   !--summation variables and then accumulate the sum over rows (or
   !--columns) and a cummulative sum of all matrix elements.
   ibase=ibase+ne
   do i=1,ne1
      rcs(i)=0
   enddo
   sumt=0
   do i=1,ne1
      do j=i,ne1
         ibase=ibase+1
         rcs(i)=rcs(i)+covl(ibase)
         sumt=sumt+covl(ibase)
         if (j.gt.i) then
            rcs(j)=rcs(j)+covl(ibase)
            sumt=sumt+covl(ibase)
         endif
      enddo
   enddo

   !--get multigroup spectrum integrals on the covariance matrix
   !--energy grid (needed for zero-sum rule test).
   call spcint(covl,spc,ne1)

   !--perform zero-sum test on the covariance matrix and apply
   !--correction if needed.
   is=0
   i=1
   do while (is.eq.0.and.i.lt.ne1)
      if ((rcs(i)/spc(i)).gt.stst) is=is+1
      i=i+1
   enddo
   if (is.eq.0) then
      call mess('sumchk','zero-sum test passed',' ')
   else
      call mess('sumchk','zero-sum test failed',&
                'applying normalization correction')
      ibase=ne+6
      do i=1,ne1
         do j=i,ne1
            ibase=ibase+1
            covl(ibase)=covl(ibase)-spc(i)*rcs(j)&
                                   -spc(j)*rcs(i)&
                                   +spc(i)*spc(j)*sumt
            if (abs(covl(ibase)).lt.sml) covl(ibase)=0
         enddo
      enddo
   endif

   !--release spc and rcs (they will be re-assigned, if
   !--necessary, upon the next entry into sumchk).
   deallocate(spc)
   deallocate(rcs)

   return

   end subroutine sumchk

   subroutine spcint(covm,spc,ne1)
   !-------------------------------------------------------------------
   ! routine to obtain the spectrum from the input tape (nendf2) and
   ! compute multigroup integrals on the covariance matrix energy grid
   !  covm(1:ncovl)  = mf35 list record (header, matrix energy grid
   !                   and the triangular covariance matrix elements).
   !
   ! the spectrum is read from either file5 or file6.  if from file5
   ! can use intega directly since the data are in a tab1 record.  if
   ! from file6 must convert the list record into an equivalent tab1
   ! before using intega for the integration.
   !
   ! in this initial version, only support file5, lf=1 and file6,
   ! law=lang=1.
   !-------------------------------------------------------------------
   use util ! provides error,mess
   use endf ! provides intega,terpa
   ! externals
   integer::ne1
   real(kr)::covm(*),spc(ne1)
   ! internals
   character(60)::strng
   character(66)::c
   integer,parameter::nnw=10000
   integer::mtt,mf56,mat,mf,mt,nk
   integer::i,nb,nw,k,lf,nr1,nr2,np1,np2,ib,ib2,iloop,ibx,nr12,np12
   integer::ir,ip,idis,izap,law,lang,lep,ib2x,na,nep,indx1,nnk,nd
   real(kr)::c1,c2,esp,pe,enext,elow,ehigh,aaa
   real(kr),dimension(:),allocatable::scr3
   real(kr)::b(nnw)

   !--allocate scratch storage and initialize the spectrum
   !--integral array
   allocate(scr3(nnw))
   do i=1,ne1
      spc(i)=0
   enddo
   mtt=mth

   !--read the endf2 input tape until reach the required mf/mtt
   call repoz(nendf2)
   call tpidio(nendf2,0,0,b,nb,nw)
   call contio(nendf2,0,0,b,nb,nw)
   call findf(math,3,0,nendf2)
   mf56=0
   do
      call contio(nendf2,0,0,b,nb,nw)
      if (mfh.eq.0.and.mth.eq.0)call contio(nendf2,0,0,b,nb,nw)
      if ((mfh.eq.5.or.mfh.eq.6).and.mth.eq.mtt) mf56=mfh
      if (mfh.gt.6) mf56=999
      if (mf56.gt.0) exit
      call tosend(nendf2,0,0,b)
   enddo
   if (mf56.eq.999) then
      write(strng,'(''no mf5 or mf6, mt='',i3,'' spectrum on nendf2'')') mtt
      call error('spcint',strng,' ')
   endif
   nk=n1h
   if (mf56.eq.5) then

      !--loop over the nk subsections.  First subsection is always a
      !--tab1 regardless of subsequent lf flag.
      do k=1,nk
         call tab1io(nendf2,0,0,scr3,nb,nw)
         lf=nint(scr3(4))
         if (lf.ne.1) then
            write(strng,'(''not ready for lf = '',i2)') lf
            call error('spcint',strng,' ')
         endif
         nr1=nint(scr3(5))
         np1=nint(scr3(6))
         ib=6+2*(nr1+np1)
         if (lf.eq.1) then
            call tab2io(nendf2,0,0,scr3(ib),nb,nw)
            nr2=nint(scr3(ib+4))
            np2=nint(scr3(ib+5))
            ib2=ib+6+2*nr2
            iloop=0
            ibx=0

            !--loop over tab1 spectrum data.  Use the first
            !--spectrum whose incident energy equals or exceeds
            !--the covariance matrix lower energy bound.
            do while (iloop.le.np2)
               iloop=iloop+1
               call tab1io(nendf2,0,0,scr3(ib2),nb,nw)
               if (nb.ne.0) ibx=ib2+nw
               do while (nb.ne.0)
                  call moreio(nendf2,0,0,scr3(ibx),nb,nw)
                  ibx=ibx+nw
                  if (ibx.ge.nnw) call error('spcint','array overflow',' ')
               enddo
               esp=scr3(ib2+1)
               nr12=nint(scr3(ib2+4))
               np12=nint(scr3(ib2+5))
               if (esp.ge.covm(1)) then
                  ir=1
                  ip=2
                  call terpa(pe,esp,enext,idis,scr3,ip,ir)
                  ir=1
                  ip=2
                  do i=1,ne1
                     elow=covm(i+6)
                     ehigh=covm(i+7)
                     call intega(aaa,elow,ehigh,scr3(ib2),ip,ir)
                     spc(i)=spc(i)+pe*aaa
                  enddo
                  iloop=np2+1
               endif
            enddo
         endif
      enddo

   else if (mf56.eq.6) then

      !--general format is a loop over nk, but the neutron (zap=1)
      !--should be the first particle defined, therefore restrict
      !--the do loop to a single iteration.
      nnk=1
      do k=1,nnk
         call tab1io(nendf2,0,0,scr3,nb,nw)
         izap=nint(scr3(1))
         if (izap.ne.1) then
            write(strng,'(''looking for mf=6,mt='',i3,'',izap=1 but '',&
                         &''found izap ='',i5)') mtt,nint(scr3(1))
            call error('spcint',strng,' ')
         endif
         law=nint(scr3(4))
         if (law.ne.1) then
            write(strng,'(''not ready for mf=6, mt='',i3,'', law = '',&
                          &i2)') mtt,law
            call error('spcint',strng,' ')
         endif
         nr1=nint(scr3(5))
         np1=nint(scr3(6))
         ib=6+2*(nr1+np1)
         if (law.eq.1) then
            call tab2io(nendf2,0,0,scr3(ib),nb,nw)
            lang=nint(scr3(ib+2))
            lep=nint(scr3(ib+3))
            nr2=nint(scr3(ib+4))
            np2=nint(scr3(ib+5))
            ib2=ib+6+2*nr2
            iloop=0

            !--loop over list spectrum data.  Use the first
            !--spectrum whose incident energy equals or exceeds
            !--the covariance matrix lower energy bound.  When
            !--found, convert into an equivalent tab1 record
            !--(overwrite the tab2 record starting at a(ib) since
            !--we don't need these data any longer) and use
            !--intega to get the spectrum integral.
            do while (iloop.le.np2)
               iloop=iloop+1
               call listio(nendf2,0,0,scr3(ib2),nb,nw)
               ib2x=ib2+nw
               do while (nb.ne.0)
                  call moreio(nendf2,0,0,scr3(ib2x),nb,nw)
                  ib2x=ib2x+nw
               enddo
               esp=scr3(ib2+1)
               nd=nint(scr3(ib2+2))
               na=nint(scr3(ib2+3))
               nw=nint(scr3(ib2+4))
               nep=nint(scr3(ib2+5))
               if (esp.ge.covm(1)) then
                  ir=1
                  ip=2
                  call terpa(pe,esp,enext,idis,scr3,ip,ir)
                  indx1=ib+6
                  scr3(ib)=0
                  scr3(ib+1)=0
                  scr3(ib+2)=0
                  scr3(ib+3)=0
                  scr3(ib+4)=1
                  scr3(ib+5)=nep
                  scr3(ib+6)=nep
                  scr3(ib+7)=lep
                  do i=1,nep,nw/nep
                     indx1=indx1+2
                     scr3(indx1)=scr3(ib2+i+6)
                     scr3(indx1+1)=scr3(ib2+i+7)
                  enddo
                  ir=1
                  ip=2
                  do i=1,ne1
                     elow=covm(i+6)
                     ehigh=covm(i+7)
                     call intega(aaa,elow,ehigh,scr3(ib),ip,ir)
                     spc(i)=spc(i)+pe*aaa
                  enddo
               endif
               iloop=np2+1
            enddo
         endif
      enddo
   endif

   !--release scratch storage
   deallocate(scr3)

   return
   end subroutine spcint

   subroutine covbin(covl)
   !---------------------------------------------------------------
   ! errorj works with absolute covariances of the probability
   ! distribution function while the endf-6 format specifies
   ! covariances of the bin probabilities. Divide covariances
   ! by the bin width.
   !  - this routine modeled after covbin in chmf35 by A.Trkov.
   !---------------------------------------------------------------
   use util ! provides mess
   ! externals
   real(kr)::covl(*)
   ! internals
   character(60)::strng
   integer::ibase,ng,midx,i,j
   real(kr)::dei,dej

   write(strng,&
     '(''converting mf35 convariance data to errorj format'')')
   call mess('covbin',strng,' ')
   ibase=6
   ng=nint(covl(6))-1
   midx=ibase+ng+1
   do i=1,ng
      dei=covl(ibase+i+1)-covl(ibase+i)
      do j=i,ng
         dej=covl(ibase+j+1)-covl(ibase+j)
         midx=midx+1
         covl(midx)=covl(midx)/(dei*dej)
     enddo
   enddo

   return
   end subroutine covbin

   subroutine rdsig(mat,mt,izap,b,sig)
   !--------------------------------------------------------------------
   ! Read cross sections from ngout tape using subroutine rdgout,
   ! re-initializing that subroutine for each new value of mat.
   !--------------------------------------------------------------------
   ! externals
   integer::mat,mt,izap
   real(kr)::b(*),sig(*)
   ! internals
   integer::matd,mfrd,matlst,matrd,mtrd,mfri,mtri,izero
   save matd,mfrd,matlst

   izero=0
   if (mt.eq.0) then
      matd=mat
      mfrd=3
      matlst=10000
   else
      matrd=mat
      if (mat.eq.0) matrd=matd
      mtrd=mt
      if (mat.ne.matlst) then
         mfri=1
         mtri=451
         call rdgout(ngout,matrd,mfri,mtri,izero,b,sig)
         matlst=mat
      endif
      call rdgout(ngout,matrd,mfrd,mtrd,izero,b,sig)
   endif

   return
   end subroutine rdsig

   subroutine rdgout(ngout,matd,mfd,mti,izap,b,sig)
   !--------------------------------------------------------------------
   ! Find the desired information from a groupr-type output tape.
   ! The flag mti=-1 is used to retrieve flux from mt=1 records.
   !--------------------------------------------------------------------
   use util ! provides repoz,error
   use endf ! provides endf routines and variables
   ! externals
   integer::ngout,matd,mfd,mti,izap
   real(kr)::b(*),sig(*)
   ! internals
   integer::mtd,nb,nwds,nz,ntw,ngtp1,ibase,isave
   integer::i,is,nl,jg,if,ib,nwds2,ibn,mtlast,ngt,iz,jzap
   integer::mtsig0(100)
   integer::mtsig=0
   character(60)::strng
   save mtlast,ngt,iz

   mtd=iabs(mti)
   if (mfd.gt.1) go to 200

   !--copy rest of this mat to a scratch file.
   call repoz(nscr2)
   call repoz(ngout)
   call tpidio(ngout,0,0,b,nb,nwds)
  100 continue
   call contio(ngout,0,0,b,nb,nwds)
   if (math.eq.matd) go to 120
   if (math.ne.-1) go to 110
   write(strng,'(''mat'',i4,'' not found.'')') matd
   call error('rdgout',strng,' ')
  110 continue
   call tomend(ngout,0,0,b)
   go to 100
  120 continue
   if (mfh.ne.mfd.or.mth.ne.mtd) then
      write(strng,'(''mf'',i2,'' mt'',i3,'' not found.'')') mfd,mtd
      call error('rdgout',strng,' ')
   endif
   nz=l2h
   ntw=n2h
   call listio(ngout,0,0,b,nb,nwds)
   ngt=l1h
   ngtp1=ngt+1
   ibase=ntw+nz+6
   isave=0
   do i=1,ngtp1
      isave=isave+1
      sig(i)=b(ibase+isave)
      if (nb.gt.0.and.ibase+isave.ge.nwds) then
         call moreio(ngout,0,0,b,nb,nwds)
         ibase=0
         isave=0
      endif
   enddo
   mti=ngt
   iz=1
   nsc=0
   call tofend(ngout,0,0,b)
   ! copy rest of material to scratch tape
   call contio(ngout,0,nscr2,b,nb,nwds)
   call tomend(ngout,0,nscr2,b)
   call atend(0,nscr2)
   mtlast=1000
   return

   !--retrieve desired cross section or flux.
   !--construct mt=3 from total minus elastic
  200 continue
   if (mtd.le.mtlast) call repoz(nscr2)
   if (mtd.eq.3) call repoz(nscr2)
   mtlast=mtd
   if (mtd.eq.3) mtd=1
   do is=1,ngt
      sig(is)=0
   enddo
  210 continue
   call contio(nscr2,0,0,b,nb,nwds)
   if (math.lt.1) then
      if (mti.ne.-2) then
         if (mtsig.gt.0) then
            do i=1,mtsig
               if (mtd.eq.mtsig0(i)) go to 212
            enddo
            mtsig=mtsig+1
            mtsig0(mtsig)=mtd
         else
            mtsig0(1)=mtd
            mtsig=1
         endif
         write(strng,'(''mf'',i2,'' mt'',i3,'' not found.'')') mfd,mtd
         call mess('rdgout',strng,'calculation continued with sigma=0')
  212    continue
      endif
      call repoz(nscr2)
      return
   endif
   jzap=0
   if (izap.ne.0) jzap=nint(c2h)
   if (mfh.eq.0.or.mth.eq.0) go to 210
   if (mfh.eq.mfd.and.mth.eq.mtd.and.jzap.eq.izap) go to 230
   if (mfh.eq.mfd.and.mth.gt.mtd.and.mti.lt.0) go to 220
   call tosend(nscr2,0,0,b)
   go to 210

   !--if the total cross-section is absent, construct the
   !--flux vector as the union of the fluxes
   !--from all partials present
  220 continue
   mti=-2
   mtd=mth
  230 continue
   nl=l1h
   nz=l2h
  240 continue
   call listio(nscr2,0,0,b,nb,nwds)
   jg=n2h
   if=1+nl*(iz-1)
   is=nl*nz+1+nl*(iz-1)
   ib=is
   if (mti.lt.0) ib=if
   if ((ib+6).gt.nwds) go to 250
   if (mtd.eq.1.or.mtlast.ne.3) sig(jg)=b(ib+6)
   if (mtd.eq.2.and.mtlast.eq.3) sig(jg)=sig(jg)-b(ib+6)
   go to 270
  250 continue
   if (nb.eq.0) call error('rdgout',&
     'bad index for b equivalent to sig(ig).',' ')
   call moreio(nscr2,0,0,b,nb,nwds2)
   ibn=ib+6-nwds
   if (ibn.le.nwds2) go to 260
   nwds=nwds+nwds2
   go to 250
  260 continue
   if (mtd.eq.1.or.mtlast.ne.3) sig(jg)=b(ibn)
   if (mtd.eq.2.and.mtlast.eq.3) sig(jg)=sig(jg)-b(ibn)
  270 continue
   do while (nb.ne.0)
      call moreio(nscr2,0,0,b,nb,nwds)
   enddo
   if (jg.lt.ngt) go to 240
   if (mti.eq.-2) go to 210
   if (mtlast.eq.3.and.mtd.eq.2) mtd=3
   if (mtlast.ne.3.or.mtd.eq.3) return
   mtd=2
   go to 210
   end subroutine rdgout

   subroutine stand(li,l,loc,lty,scr)
   !--------------------------------------------------------------------
   ! Read and store the appropriate data from nstan.
   !--------------------------------------------------------------------
   use util ! provides sigfig
   use endf ! provides endf routines and variables
   ! externals
   integer::li,l,loc(*),lty
   real(kr)::scr(*)
   ! internals
   integer::matstd,mtstd,i,nb,nw,nc,ni,ic,ii,np,locli,idone

   elr=sigfig(c1h,ndig,0)
   ehr=sigfig(c2h,ndig,0)
   matstd=l1h
   mtstd=l2h

   !--redefine standard if necessary
   if (nas.gt.0) then
      idone=0
      i=0
      do while (i.lt.nas.and.idone.eq.0)
         i=i+1
         if (matstd.eq.matb(i).and.mtstd.eq.mtb(i)) then
            matstd=matc(i)
            mtstd=mtc(i)
            idone=1
         endif
      enddo
   endif
   call findf(matstd,mfcov,mtstd,nstan)
   call contio(nstan,0,0,scr(l),nb,nw)

   !--first subsection is the one we want
   call contio(nstan,0,0,scr(l),nb,nw)
   nc=n1h
   ni=n2h

   !--skip over nc sub-subsections
   if (nc.ne.0) then
      do ic=1,nc
         call contio(nstan,0,0,scr(l),nb,nw)
         call listio(nstan,0,0,scr(l),nb,nw)
         do while (nb.ne.0)
            call moreio(nstan,0,0,scr(l),nb,nw)
         enddo
      enddo
   endif

   !--loop over ni sub-subsections
   do ii=1,ni
      call listio(nstan,0,0,scr(l),nb,nw)
      li=li+1
      loc(li)=l
      scr(l)=elr
      scr(1+l)=ehr
      if (l2h.eq.6) scr(l+2)=int((n1h-1)/n2h)
      scr(l+4)=lty
      np=n1h
      l=l+nw
      do while (nb.ne.0)
         call moreio(nstan,0,0,scr(l),nb,nw)
         l=l+nw
      enddo
      locli=loc(li)+5
      do i=1,np
         scr(i+locli)=sigfig(scr(i+locli),ndig,0)
      enddo
   enddo
   return
   end subroutine stand

   subroutine resprx(nwscr,a)
   !--------------------------------------------------------------------
   ! prepare tables containing the resonance-parameter contributions
   ! to coarse-group covariances.
   !--------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use util ! provides error,openz,repoz,closz
   use samm ! provides desammy
   use mainio ! provides nsyso
   ! externals
   integer::nwscr
   real(kr)::a(nwscr)
   ! internals
   integer::indx
   integer::nb,nw,nngn,ig,ner
   integer::ip1,ip2,i,iest,ieed,iscr,jscr,nro,mls
   integer,parameter::mxlru2=100
   real(kr)::za,awr,ee
   real(kr)::amur(3,mxlru2)
   character(60)::strng1,strng2

   !--initialize
   nresg=0
   ifresr=0
   ifunrs=0
   nscr6=0
   indx=0
   iscr=1
   if (mfcov.eq.33.and.mf32.ne.0) then
      call repoz(nendf)
      call tpidio(nendf,0,0,a(iscr),nb,nw)
      nngn=ngn*(ngn+1)/2
      do ig=1,nngn
         cff(ig)=0
         cee(ig)=0
         cgg(ig)=0
         ctt(ig)=0
         uff(ig)=0
         ugg(ig)=0
         uee(ig)=0
         utt(ig)=0
      enddo
      nngn=ngn*ngn
      do ig=1,nngn
         cfg(ig)=0
         cef(ig)=0
         ceg(ig)=0
         ufg(ig)=0
         uef(ig)=0
         ueg(ig)=0
      enddo
      nscr6=16
      if (nendf.lt.0) nscr6=-nscr6
      call openz(nscr6,1)
      call repoz(nscr6)
      call rdumrd2(matd,nendf,nscr6,a(iscr),amur,mxlru2,nwscr)
      call repoz(nscr6)

      !--start reading the covariance files
      call findf(matd,32,151,nendf)
      call contio(nendf,0,0,a(iscr),nb,nw)
      za=c1h
      awr=c2h
      nis=n1h
   endif

   !--loop over isotopes
   do isrr=1,nis
      call contio(nendf,0,0,a(iscr),nb,nw)
      abn=c2h
      lfw=l2h
      ner=n1h

      !--loop over energy ranges
      do ier=1,ner
         indx=indx+1
         write (*,'(/'' Energy range : '',i5,''/'',i5)') ier, ner
         call contio(nendf,0,0,a(iscr),nb,nw)
         elr=c1h
         ehr=c2h
         ehg=ehr
         elg=elr
         ip1=0
         ip2=0
         do i=2,ngn+1
            ee=egn(i)
            if (ip1.eq.0.and.ee.gt.elr) then
               ip1=1
               elg=egn(i-1)
               iest=i-1
            endif
            if (ip2.eq.0.and.ee.ge.ehr) then
               ip2=1
               ehg=egn(i)
               ieed=i-1
            endif
         enddo
         if (elg.ge.ehr) ieed=0
         lru=l1h
         lrf=l2h
         nro=n1h
         naps=n2h
         if (lru.eq.1.and.lrf.ge.1.and.lrf.le.3) go to 130
         if (lru.eq.1.and.lrf.ge.1.and.lrf.le.7) go to 130
         if (lru.eq.2.and.lrf.ge.1.and.lrf.le.2) go to 130
            write(strng2,'(''lrf='',i4,''  lru='',i4)') lrf,lru
            call error('resprx',&
              'illegal or no coding data structure in mf32',strng2)
  130    continue
         if (nro.ne.0) then
            write(strng2,'(''nro='',i4)') nro
            call error('resprx',&
              'illegal or unrecognized data structure in mf32',strng2)
         endif
         if (lrf.eq.7.and.isammy.eq.0) call error('resprx',&
            'cannot handle lrf=7 RML resonance representation',&
            'with isammy=0')
         call contio(nendf,0,0,a(iscr),nb,nw)
         spi=c1h
         spifac=1/(2*spi+1)
         ap=c2h
         lcomp=l2h
         if (lrf.ne.7) then
            nls=nlspepi(indx)
            if (n1h.gt.nls) then
               call error('resprx',&
                          'mf2/mf32 l-state mis-match',&
                          'probable evaluation file error')
            else if (n1h.lt.nls) then
               write(strng1,'(''mf2 nls='',i1,'', but mf32 nls='',i1)')&
                                           nls,                   n1h
               call mess('resprx',strng1,&
                         'continue with partial urr covariance data')
               nls=n1h
            endif
         else
            nls=n1h
         endif

         !--Scattering radius uncertainty processing.  May be no data,
         !  may be data on file, may be user data.
         isr=n2h
         if (lrf.ne.7) then
            if (isr.eq.1) then
               jscr=iscr+6
               if (lrf.eq.1.or.lrf.eq.2) then
                  call contio(nendf,0,0,a(jscr),nb,nw)
                  if (isru.eq.0) then
                     dap=a(jscr+1)
                  else
                     dap=ap*dap
                  endif
                  do i=1,nls
                     dap3(i)=dap
                  enddo
               else if (lrf.eq.3) then
                  call listio(nendf,0,0,a(jscr),nb,nw)
                  if (isru.eq.0) then
                     mls=nint(a(jscr+4))
                     dap=a(jscr+6)
                  else
                     mls=1
                     dap=ap*dap
                  endif
                  if (mls.eq.1) then
                     do i=1,nls
                        dap3(i)=dap
                     enddo
                  else if (mls.gt.1.and.mls.eq.nls+1) then
                     do i=1,nls
                        dap3(i)=a(jscr+6+i)
                     enddo
                  else if (mls.gt.1.and.mls.lt.nls+1) then
                     do i=1,mls-1
                        dap3(i)=a(jscr+6+i)
                     enddo
                     do i=mls,nls
                        dap3(i)=dap
                     enddo
                  else
                     write(strng1,'(''mls='',i1,'', nls='',i1,&
                                   &''are inconsistent'')')mls,nls
                     write(strng2,'(''will ignore scattering '',&
                                   &''radius uncertainty'')')
                     call mess('resprx',strng1,strng2)
                     isr=0
                     dap=0
                  endif
               else
                  write(strng1,'(''not ready for isr=1, lrf='',i1)') lrf
                  call error('resprx',strng1,' ')
               endif
            else if (isr.eq.0.and.isru.ne.0) then
               isr=1
               dap=ap*dap
               if (nls.gt.0) then
                  do i=1,nls
                     dap3(i)=dap
                  enddo
               endif
            else if (isr.ne.0) then
               write(strng1,*) 'illegal isr',isr
               call error('resprx',strng1,' ')
            endif
            if (isr.ne.0.and.isru.ne.0) then
               write(strng1,'(''user override for scattering radius '',&
                             &''uncertainty'')')
               write(strng2,'(''use DAP ='',f6.3,''*AP for all L'')')dap/ap
               call mess('resprx',strng1,strng2)
            endif
         else if (lrf.eq.7.and.(isr.ne.0.or.isru.ne.0)) then
            call mess('resprx','scat. radius unc not ready for lrf=7',' ')
         endif

         lptr=iscr+6

         !--Unresolved
         if (lru.eq.2) then
            call rpxunr(a,amur,mxlru2,iest,ieed,nwscr)

         !--Resolved with sammy method
         else if (nmtres.gt.0) then
            call rpxsamm(nwscr,a,ier)

         !--Resolved with errorj method
         else
            if (lcomp.eq.0) then
               call rpxlc0(nwscr,a)
            else if (lcomp.eq.1.or.lcomp.eq.2) then
               call rpxlc12(nwscr,a,iest,ieed)
            endif
         endif

         !--close loops on isotope and energy range
      enddo
   enddo
   call closz(nscr6)
   if (nmtres.gt.0) call desammy

   return
   end subroutine resprx

   subroutine rpxsamm(nwscr,a,ier)
   !--------------------------------------------------------------------
   ! Resolved resonances with sammy method: lrf=3 or 7, lcomp=1 or 2.
   !--------------------------------------------------------------------
   use endf ! provides listio, moreio, intgio
   use util ! provides error
   use samm ! provides rdsammy
   use physics ! provides amassn
   use mainio
   ! externals
   integer::nwscr,ier
   real(kr)::a(nwscr)
   ! internals
   integer::jnow,nodes,nnn,nm,nx
   integer::i1,i2,n1,n2,nn1,nn2,nn2p,n3,j,ig,ig2
   integer::l1,l2,l3,k2,k,nind,lbg,nn
   integer::nch,ip,ncoef,nresp
   integer::nsrs,nlrs,njsx,nparb,idis,lord
   integer::i,nb,nw
   integer::isrrr
   real(kr)::el,eh,spin,ap,awri,aw,fact,en,rat,frac
   real(kr)::e,enext,wt,wtl,elo,ehi,enxt,ee,eel,tmp
   integer,parameter::nodmax=30000
   real(kr)::enode(nodmax)
   integer,parameter::maxb=90000
   real(kr)::b(maxb)

   real(kr),dimension(:),allocatable::sdev
   real(kr),dimension(:,:),allocatable::cov
   real(kr),dimension(:,:,:),allocatable::sens
   real(kr),dimension(:),allocatable::sigp,sigpn,sigpl
   real(kr),dimension(:,:),allocatable::siga,sigan,sigal
   real(kr),dimension(:,:),allocatable::sigx,sigxn,sigxl
   real(kr),dimension(:),allocatable::sflx
   real(kr),dimension(:,:),allocatable::sigs
   character(60)::strng1,strng2
   real(kr),parameter::eps=0.01e0_kr
   real(kr),parameter::epm=1.0e-7_kr
   real(kr),parameter::rc1=.123e0_kr
   real(kr),parameter::rc2=.08e0_kr
   real(kr),parameter::third=0.333333333e0_kr
   real(kr),parameter::zero=0
   real(kr),parameter::ten=10

   !--read in the resonance parameters for computing cross sections
   isrrr=1
   ier=1
   call rskiprp(nscr6,b,isrrr,ier,maxb)
   jnow=1
   call contio(nscr6,0,0,b(jnow),nb,nw)
   el=c1h
   enode(1)=sigfig(el,7,+1)
   eh=c2h
   enode(2)=sigfig(eh,7,-1)
   nodes=2
   lru=l1h
   lrf=l2h
   jnow=jnow+6
   call rdsammy(nscr6,ier,jnow,nro,naps,lrf,el,eh,&
     enode,nodes,nodmax,spin,ap,b,maxb)
   nresp=1
   call ppsammy(1,ncoef,nresp)

   !--read the parameter covariances
   if (lrf.ne.7.or.lcomp.ne.2) then
      allocate(sdev(nresp))
      allocate(cov(nresp,nresp))
   else
      allocate(sdev(2*nresp))
      allocate(cov(2*nresp,2*nresp))
   endif
   cov=0

   !--lrf=3 and lcomp=1
   if (lrf.eq.3.and.lcomp.eq.1) then
      l1=lptr
      call contio(nendf,0,0,a(l1),nb,nw)
      awri=c1h
      nsrs=n1h
      nlrs=n2h
      do i=1,nsrs
         l1=lptr
         call listio(nendf,0,0,a(l1),nb,nw)
         l1=l1+nw
         do while (nb.gt.0)
            if (l1.gt.nwscr) call error('rpxsamm','storage exceeded',' ')
            call moreio(nendf,0,0,a(l1),nb,nw)
            l1=l1+nw
         enddo
         mpar=nint(a(lptr+2))
         nrb=nint(a(lptr+5))
         nvs1=6*nrb
         npar=mpar*nrb
      enddo

   !--lrf=3 and lcomp=2
   else if (lrf.eq.3.and.lcomp.eq.2) then
      nrb=0
      nind=1
      lbg=lptr
      do nn=1,nls
         l1=lptr
         call listio(nendf,0,0,a(l1),nb,nw)
         l1=l1+nw
         do while (nb.gt.0)
            call moreio(nendf,0,0,a(l1),nb,nw)
            l1=l1+nw
         enddo
         if (l1.gt.nwscr) then
            write(strng2,'(''require='',i8,''  supply='',i8,&
              &'' for nwds given in s.covout'')') l1,nwscr
            call error('rpxsamm','storage exceeded',strng2)
         endif
         if (nn.eq.1) then
            awri=a(lptr)
            arat=awri/(awri+1)
            aw=amassn*awri
            ra=rc1*aw**third+rc2
            ral=ra
            apl=ap
         endif
         mpar=nint(a(lptr+2))
         nrb=nrb+nint(a(lptr+5))
         l3=lptr
         do n2=1,nint(a(lptr+5))
            l3=l3+12
            cov(nind,nind)=a(l3)
            cov(nind+1,nind+1)=a(l3+2)
            cov(nind+2,nind+2)=a(l3+3)
            nind=nind+3
         enddo
         l3=lbg+(nrb-nint(a(lptr+5))+1)*6
         l2=lptr+6
         do n2=1,nint(a(lptr+5))
            do n3=1,6
               a(l3+n3-1)=a(l2+n3-1)
            enddo
            l3=l3+6
            l2=l2+12
         enddo
         lptr=l1
      enddo
      call contio(nendf,0,0,a(l3),nb,nw)
      nx=l1h
      if (nx.eq.0) nx=2
      nnn=l2h
      nm=n1h
      do n2=1,nm
         nw=nx
         call intgio(nendf,0,0,a(l3),nb,nw)
         nn1=nint(a(l3))
         nn2=nint(a(l3+1))
         nn2p=nn2-1
         fact=2*10**nx
         do n3=3,nw
            nn2p=nn2p+1
            if (nn2p.ge.nn1) exit
            l2=l3-1+n3
            if (nint(a(l2)).gt.0) then
               cov(nn2p,nn1)=((2*a(l2)+1)/fact)*cov(nn2p,nn2p)*cov(nn1,nn1)
            else if (nint(a(l2)).lt.0) then
               cov(nn2p,nn1)=(-(-2*a(l2)+1)/fact)*cov(nn2p,nn2p)*cov(nn1,nn1)
            endif
            cov(nn1,nn2p)=cov(nn2p,nn1)
         enddo
      enddo

   !--lrf=7 and lcomp=1
   else if (lrf.eq.7.and.lcomp.eq.1) then
      call contio(nendf,0,0,a(lptr),nb,nw)
      nsrs=n1h
      nlrs=n2h
      call contio(nendf,0,0,a(lptr),nb,nw)
      njsx=l1h
      do i=1,nls
         l1=lptr
         call listio(nendf,0,0,a(l1),nb,nw)
         nch=l1h
         nrb=l2h
         nx=n2h
         do while (nb.ne.0)
            l1=l1+nw
            if (l1.gt.nwscr) call error('rpxsam','storage exceeded','in a')
            call moreio(nendf,0,0,a(l1),nb,nw)
         enddo
      enddo
      l1=lptr
      call listio(nendf,0,0,a(l1),nb,nw)
      nparb=n2h
      do while (nb.ne.0)
         l1=l1+nw
         if (l1.gt.nwscr) call error('rpxsam','storage exceeded','in a')
         call moreio(nendf,0,0,a(l1),nb,nw)
      enddo
      l2=lptr+5
      do i1=1,nparb
         do i2=i1,nparb
            l2=l2+1
            cov(i1,i2)=a(l2)
            if (i1.ne.i2) cov(i2,i1)=cov(i1,i2)
         enddo
      enddo

   !--lrf=7 and lcomp=2
   else if (lrf.eq.7.and.lcomp.eq.2) then
      call listio(nendf,0,0,a,nb,nw)
      ip=0
      do i=1,nls
         l1=lptr
         call listio(nendf,0,0,a(l1),nb,nw)
         nch=n2h
         call listio(nendf,0,0,a(l1),nb,nw)
         nx=n2h
         do while (nb.ne.0)
            l1=l1+nw
            if (l1.gt.nwscr) call error('rpxsam','storage exceeded','in a')
            call moreio(nendf,0,0,a(l1),nb,nw)
         enddo

         !--extract the standard deviations (uncertainties)
         do i1=1,nx
            do i2=1,1+nch
               ip=ip+1
               sdev(ip)=a(lptr+6+(i1-1)*12+6+i2-1)
            enddo
         enddo
      enddo
      npar=ip

      !--read the compact correlation matrix
      l3=lptr
      call contio(nendf,0,0,a(l3),nb,nw)
      nx=l1h
      if (nx.eq.0) nx=2
      nnn=l2h
      nm=n1h
      do i=1,npar
         cov(i,i)=1
      enddo
      do n2=1,nm
         nw=nx
         call intgio(nendf,0,0,a(l3),nb,nw)
         nn1=nint(a(l3))
         nn2=nint(a(l3+1))
         nn2p=nn2-1
         fact=2*ten**nx
         do n3=3,nw
            nn2p=nn2p+1
            if (nn2p.ge.nn1) exit
            l2=l3-1+n3
            if (nint(a(l2)).gt.0) then
               cov(nn2p,nn1)=((2*a(l2)+1)/fact)
            else if (nint(a(l2)).lt.0) then
               cov(nn2p,nn1)=(-(-2*a(l2)+1)/fact)
            endif
            cov(nn1,nn2p)=cov(nn2p,nn1)
         enddo
      enddo

      !--convert the correlation matrix to a covariance matrix
      do i=1,npar
         do j=1,npar
            cov(i,j)=cov(i,j)*sdev(i)*sdev(j)
         enddo
      enddo
   endif

   !--calculate the group-wise sensitivities using sammy derivatives
   nnn=nmt
   if (nmtres+2.gt.nnn) nnn=nmtres+2
   allocate(sens(nmt,ngn,nresp))
   allocate(sigp(nnn))
   allocate(sigpn(nnn))
   allocate(sigpl(nnn))
   allocate(siga(ncoef,nmtres))
   allocate(sigan(ncoef,nmtres))
   allocate(sigal(ncoef,nmtres))
   allocate(sigx(nresp,nmtres))
   allocate(sigxn(nresp,nmtres))
   allocate(sigxl(nresp,nmtres))
   allocate(sflx(ngn))
   allocate(sigs(ngn,nnn))
   e=0
   call egtwtf(e,enext,idis,lord,wtl)
   i2=1
   do i1=1,ngn
      sflx(i1)=0
      do n1=1,nmtres+2
         sigs(i1,n1)=0
      enddo
      do n1=1,nmt
         do n3=1,nresp
            sens(n1,i1,n3)=0
         enddo
      enddo
   enddo
   do i1=1,ngn
      frac=0
      elo=egn(i1)
      ehi=egn(i1+1)
      e=elo
      call cssammy(e,sigpl,sigal,sigxl,mmtres,nmtres,ncoef,nresp,ier)
      if (e.gt.enode(nodes)) then
         do n1=1,nmtres+2
            sigpl(n1)=0
         enddo
         do n1=1,nmtres
            do n2=1,nresp
               sigxl(n2,n1)=0
            enddo
         enddo
      endif
      call egtwtf(e,enext,idis,lord,wtl)
      ee=elo
      eel=ee
      enext=ehi
   100 continue
      if (ee.gt.enode(nodes)) then
         enext=ehi
      else if (ee.eq.enode(nodes)) then
         enext=sigfig(enode(nodes),7,+1)
      else
         if ((1+eps)*ee.lt.enext) enext=(1+eps)*ee
         i2=1
         do while (enode(i2).le.ee.and.i2.lt.nodes)
            i2=i2+1
         enddo
         if (i2.le.nodes.and.enode(i2).lt.enext) enext=enode(i2)
      endif
      ee=enext
   110 continue
      call cssammy(ee,sigp,siga,sigx,mmtres,nmtres,ncoef,nresp,ier)
      if (ee.gt.enode(nodes)) then
         do n1=1,nmtres+2
            sigp(n1)=0
         enddo
         do n1=1,nmtres
            do n2=1,nresp
               sigx(n2,n1)=0
            enddo
         enddo
      endif
      call egtwtf(ee,enext,idis,lord,wt)
      e=(ee+eel)/2
      call cssammy(e,sigpn,sigan,sigxn,mmtres,nmtres,ncoef,nresp,ier)
      if (e.gt.enode(nodes)) then
         do n1=1,nmtres+2
            sigpn(n1)=0
         enddo
         do n1=1,nmtres
            do n2=1,nresp
               sigxn(n2,n1)=0
            enddo
         enddo
      endif
      if (eel.lt.e) then
         do n1=1,nmtres+2
            if (ee.ge.enode(nodes)) exit
            if (abs(sigpn(n1)-(sigp(n1)+sigpl(n1))/2).gt.eps*sigpn(n1)+epm) then
               ee=e
               go to 110
            endif
         enddo
      else
         write(strng1,'(''convergence issue for e='',1p,e10.3)') eel
         call mess('rpxsamm',strng1,'check reconstructed xs for steep increase')
      endif
      k=0
      do i=1,nek
         if (ee.ge.ek(i).and.ee.lt.ek(i+1)) k=i
      enddo
      do n1=1,nmt
         do n2=1,nmt
            k2=0
            if (mts(n2).eq.2) k2=2
            if (mts(n2).eq.18) k2=3
            if (mts(n2).eq.102) k2=4
            do i=3,nmtres
               if (mts(n2).eq.mmtres(i)) k2=i+2
            enddo
            if (k2.ne.0) then
               sigs(i1,n1)=sigs(i1,n1)+akxy(n2,n1,k)*&
                 (ee-eel)*(sigp(k2)*wt+sigpl(k2)*wtl)/2
            endif
         enddo
      enddo
      do n1=1,nmt
         do n2=1,nmt
            k2=0
            do i=1,nmtres
               if (mts(n2).eq.mmtres(i)) k2=i
            enddo
            if (k2.ne.0) then
               do n3=1,nresp
                  sens(n1,i1,n3)=sens(n1,i1,n3)+akxy(n2,n1,k)*&
                    (ee-eel)*(sigx(n3,k2)*wt+sigxl(n3,k2)*wtl)/2
               enddo
            endif
         enddo
      enddo
      sflx(i1)=sflx(i1)+(ee-eel)*(wt+wtl)/2
      if (sigp(1).ne.zero) frac=frac+(ee-eel)*(wt+wtl)/2
      do n1=1,nmtres+2
         sigpl(n1)=sigp(n1)
      enddo
      eel=ee
      wtl=wt
      do n1=1,nmtres
         do n3=1,ncoef
            sigal(n3,n1)=siga(n3,n1)
         enddo
         do n3=1,nresp
            sigxl(n3,n1)=sigx(n3,n1)
         enddo
      enddo
      if (ee.lt.ehi) go to 100

      !--adjust to groupr cross sections
      frac=frac/sflx(i1)
      do n1=1,nmt
         rat=1
         if (sigs(i1,n1).ne.zero) rat=cflx(i1)*csig(i1,n1)/sigs(i1,n1)
         do n3=1,nresp
            sens(n1,i1,n3)=sens(n1,i1,n3)*rat*frac
         enddo
      enddo

      !--continue the group loop
   enddo

   deallocate(sigp)
   deallocate(sigpl)
   deallocate(siga)
   deallocate(sigal)
   deallocate(sigx)
   deallocate(sigxl)
   deallocate(sigs)

   !--fold sensitivities with covariances to get cross reaction values
   allocate(crr(ngn,ngn,nmt,nmt))
   do ig=1,ngn
      do ig2=1,ngn
         do n1=1,nmt
            do n2=1,nmt
               crr(ig,ig2,n1,n2)=0
            enddo
         enddo
      enddo
   enddo
   do ig=1,ngn
      do ig2=1,ngn
         do i=1,nresp
            do j=1,nresp
               tmp=cov(i,j)
               if (tmp.ne.zero) then
                  do n1=1,nmt
                     do n2=1,nmt
                        crr(ig,ig2,n1,n2)=crr(ig,ig2,n1,n2)+&
                          tmp*sens(n1,ig,i)*sens(n2,ig2,j)
                     enddo
                  enddo
               endif
            enddo
         enddo
      enddo
   enddo

   !--clean up allocated arrays
   deallocate(sens)
   deallocate(cov)

   return
   end subroutine rpxsamm

   subroutine rpxlc0(nwscr,a)
   !--------------------------------------------------------------------
   ! Resolved resonances: lru=1, lcomp=0
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use physics ! provides amassn
   use util ! provides error,mess
   ! externals
   integer::nwscr
   real(kr)::a(nwscr)
   ! internals
   integer::nl,iloc,nrs,i,nr,ig,j,loop,igind,inow,id,ig1,ig2,ij,ii
   integer::nb,nw,jd
   integer,parameter::maxnls=10
   integer,parameter::maxe=600000
   integer::istloc(maxnls)
   real(kr)::awri,aw,sum,den,fl,ajmax,aj
   real(kr)::er,rho,ser,per,corr,ebc,e1,ekp
   real(kr)::er1,er2,er3,er4,er5,er6,s
   real(kr)::sig(maxe,5),gsig(4,2501,6),sig1(4)
   real(kr)::sens(4,6,2501),cov(5,5)
   real(kr)::ag(6),aa(3),aa2(3)
   character(60)::strng1,strng2
   logical::lneger
   real(kr),parameter::rc1=.123e0_kr
   real(kr),parameter::rc2=.08e0_kr
   real(kr),parameter::third=0.333333333e0_kr
   real(kr),parameter::half=0.5e0_kr
   real(kr),parameter::zero=0

   lneger=.false.

   if (isr.ne.0) then
      write(strng1,'(''lcomp=0, scattering radius uncertainty is '',&
                    &''not included'')')
      call mess('rpxlc0',strng1,' ')
   endif

   !--resolved resonance parameters (lru=1)
   !--compatible resolved resonance subsection format (lcomp=0)
   if (nls.gt.maxnls) then
      write(strng2,'(''nls='',i8,''  maxnls='',i8)') nls,maxnls
      call error('rpxlc0','storage exceeded.',strng2)
   endif
   do nl=1,nls
      istloc(nl)=lptr
      call listio(nendf,0,0,a(lptr),nb,nw)
      lptr=lptr+nw
      if (nb.gt.0) then
         if ((lptr+nb-1).gt.nwscr) then
            write(strng2,'(''need '',i7,'' words but only '',&
              &i7,'' allocated'')') lptr+nb-1,nwscr
            call error('rpxlc0','scr storage exceeded',strng2)
         endif
         do while (nb.ne.0)
            call moreio(nendf,0,0,a(lptr),nb,nw)
            lptr=lptr+nw
         enddo
      endif
   enddo

   !--loop over l states
   do 160 nl=1,nls
      iloc=istloc(nl)
      awri=a(iloc)
      ll=nint(a(iloc+2))
      nrs=nint(a(iloc+5))
      arat=awri/(awri+1)
      aw=amassn*awri
      ra=rc1*aw**third+rc2
      ral=ra
      apl=ap
      if (naps.eq.1) then
         ral=apl
         ra=ap
      endif
      if (lrf.eq.2) then
         sum=0
         den=4*spi+2
         fl=ll
         ajmin=abs(abs(spi-fl)-half)
         ajmax=spi+fl+half
         nj=nint(ajmax-ajmin+1)
         aj=ajmin
         do i=1,min(10,nj)
            gj(i)=(2*aj+1)/den
            aj=aj+1
            sum=sum+gj(i)
         enddo
         diff=2*fl+1-sum
      endif
      inow=iloc+6

      !--loop over all resonances
      do 170 nr=1,nrs
         er=abs(a(inow))
         rho=cwaven*arat*sqrt(er)*ral
         call efacts(ll,rho,ser,per)
         aa(1)=ser
         aa(2)=per
         aa(3)=0

         do ig=1,ngn
            do j=1,6
               do i=1,4
                  gsig(i,ig,j)=0
                  sens(i,j,ig)=0
              enddo
            enddo
         enddo

         do 180 loop=1,6
            do j=1,6
               ag(j)=a(inow+j-1)
            enddo
            do j=1,3
               aa2(j)=aa(j)
            enddo
            if (loop.eq.2) then
               if (ag(1).lt.zero) then
                  ag(1)=ag(1)*1.0001e0_kr
                  lneger=.false.
               else
                  do ig1=1,ngn
                     if (ag(1).ge.egn(ig1).and.ag(1).lt.egn(ig1+1)) exit
                  enddo
                  e1=ag(1)*1.0001e0_kr
                  do ig2=ig1,ngn
                     if (e1.ge.egn(ig2).and.e1.lt.egn(ig2+1)) exit
                  enddo
                  if (ig1.eq.ig2) then
                     ag(1)=e1
                     lneger=.false.
                  else
                     ag(1)=ag(1)*0.9999e0_kr
                     lneger=.true.
                  endif
               endif
               rho=cwaven*arat*sqrt(abs(ag(1)))*ral
               call efacts(ll,rho,ser,per)
               aa2(1)=ser
               aa2(2)=per
            else if (loop.eq.3) then
               go to 180
            else if (loop.ge.4.and.loop.le.6) then
               if (ag(loop).eq.0.) go to 180
               ag(loop)=ag(loop)*1.01e0_kr
               ag(3)=ag(4)+ag(5)+ag(6)
            endif
            e1=elg
            ii=0
            if (nr.eq.1) then
               er1=er/10
               er2=er*10
            else if (er.le.1.e+2_kr) then
               er1=er/4
               er2=er*4
            else
               er1=er/2.5e0_kr
               er2=er*2.5e0_kr
            endif
            if (loop.eq.2) then
               er3=abs(ag(1))*0.995e0_kr
               er4=abs(ag(1))*1.005e0_kr
               er5=abs(ag(1))*0.9992e0_kr
               er6=abs(ag(1))*1.0008e0_kr
            else
               er3=er*0.995e0_kr
               er4=er*1.005e0_kr
               er5=er*0.9992e0_kr
               er6=er*1.0008e0_kr
            endif
            go to 220

  210       continue
            if (e1.ge.er5.and.e1.le.er6) then
               ekp=1.000001e0_kr
            else if (e1.ge.er3.and.e1.le.er4) then
               ekp=1.00001e0_kr
            else if (e1.ge.er1.and.e1.le.er2) then
               ekp=1.0018e0_kr
            else
               ekp=1.02e0_kr
            endif
            e1=e1*ekp
            ebc=e1/ekp
            if (ebc.lt.elr.and.e1.gt.elr) e1=elr
            if (ebc.lt.ehr.and.e1.gt.ehr) e1=ehr

  220       continue
            if (e1.gt.ehg) e1=ehg
            if (e1.ge.elr.and.e1.le.ehr) then
               if (lrf.eq.1) then
                  call ssslbw(e1,sig1,ag,aa2)
               else if (lrf.eq.2) then
                  call ssmlbw(e1,sig1,ag,aa2)
               else
                  write(strng2,'(''lrf='',i4,'' for lcomp=0'')')lrf
                 call error('rpxlc0','not allowed lrf.',strng2)
               endif
            else
               do i=1,4
                  sig1(i)=0
               enddo
            endif
            ii=ii+1
            if (ii.gt.maxe) call error('rpxlc0',&
              'number of pointwise xsec of resonance exceeded.',&
     &         'please increase the maxe parameter.')
            do i=1,4
               sig(ii,i)=sig1(i)
            enddo
            sig(ii,5)=e1
            if (e1.ge.ehg) go to 230
            go to 210

  230       continue
            call rpxgrp(ngn,egn,sig,ii,gsig(1,1,loop),a,nwscr)
  180    continue

         do ig=1,ngn
            if (gsig(1,ig,1).gt.zero) then
               do j=1,4
                  do i=1,5
                     if (gsig(j,ig,i+1).le.zero) cycle
                     s=gsig(j,ig,i+1)-gsig(j,ig,1)
                     if (i.eq.1) then
                        ! Parameter is resonance energy
                        if (lneger) s=-s
                        s=10000*s/er
                        ij=1
                     else if (i.eq.2) then
                        ! Parameter is total width (irrelevant)
                        ij=5
                     else
                        ! Parameter is width
                        if (a(inow+i).eq.zero) cycle
                        s=100*s/a(inow+i)
                        ij=i-1
                     endif
                     if (abs(s).ge.1.e-10_kr) sens(j,ij,ig)=s*cflx(ig)*abn
                  enddo
               enddo
            endif
         enddo
         do j=1,5
            do i=1,5
               cov(i,j)=0
            enddo
         enddo

         inow=inow+6
         cov(1,1)=a(inow)
         do i=2,5
            do j=2,i
               inow=inow+1
               cov(i,j)=a(inow)
               cov(j,i)=a(inow)
            enddo
         enddo

         !--check for bad covariance data
         inow=inow+2
         if (iverf.eq.6) then
            ! For ENDFB/6 only four parameters per resonance
            jd=5
            do i=1,5
               id=i
               if (cov(i,5).ne.zero) then
                  write(strng2,'('' resonance parameters '',i1,'' and '',i1,&
                    &'' at er='',1pe12.4)') id,jd,er
                  call error('rpxlc0','bad covariance data for',strng2)
               endif
            enddo
         endif
         do i=1,5
            id=i
            jd=i
            if (cov(i,i).lt.zero) then
               write(strng2,'('' resonance parameters '',i1,'' and '',i1,&
                 &'' at er='',1pe12.4)') id,jd,er
               call error('rpxlc0','bad covariance data for',strng2)
            endif
         enddo
         do i=1,5
            id=i
            do j=1,5
               jd=j
               if (cov(i,i).gt.zero.and.cov(j,j).gt.zero) then
                  corr=cov(i,j)/sqrt(cov(i,i)*cov(j,j))
                  if (abs(corr).ge.1.0001e0_kr) then
                     if (abs(corr).gt.2.e0_kr) then
                        write(strng2,'('' resonance parameters '',i1,&
                          &'' and '',i1,'' at er='',1pe12.4)') id,jd,er
                     endif
                     call error('rpxlc0','bad covariance data for',strng2)
                     write(strng1,'(''correlation coeff.='',f8.4)') corr
                     write(strng2,'(''for resonance parameters '',i1,&
                       &'' and '',i1,'' at er='',1pe12.4)') i,j,er
                     call mess('rpxlc0',strng1,strng2)
                  endif
               else
                  if (cov(i,j).ne.zero) then
                     write(strng2,'('' resonance parameters '',i1,&
                       &'' and '',i1,'' at er='',1pe12.4)') id,jd,er
                     call error('rpxlc0','bad covariance data for',strng2)
                  endif
               endif
            enddo
         enddo

         igind=0
         do ig=1,ngn
            if (gsig(1,ig,1).le.zero) then
               igind=igind+(ngn-ig+1)
               cycle
            endif
            do ig2=ig,ngn
               igind=igind+1
               do i=1,5
                  if (sens(1,i,ig).eq.zero.and.sens(2,i,ig).eq.zero&
                    .and.sens(3,i,ig).eq.zero.and.sens(4,i,ig).eq.zero) cycle
                  do j=1,5
                     if (abs(cov(i,j)).le.zero) cycle
                     cff(igind)=cff(igind)+cov(i,j)*sens(3,i,ig)*sens(3,j,ig2)
                     cgg(igind)=cgg(igind)+cov(i,j)*sens(4,i,ig)*sens(4,j,ig2)
                     cee(igind)=cee(igind)+cov(i,j)*sens(2,i,ig)*sens(2,j,ig2)
                     ctt(igind)=ctt(igind)+cov(i,j)*sens(1,i,ig)*sens(1,j,ig2)
                  enddo
               enddo
            enddo
            if (ig.gt.nresg) nresg=ig
         enddo

         igind=0
         do ig=1,ngn
            do ig2=1,ngn
               igind=igind+1
               do i=1,5
                  if (sens(1,i,ig).eq.zero.and.sens(2,i,ig).eq.zero&
                   .and.sens(3,i,ig).eq.zero.and.sens(4,i,ig).eq.zero) cycle
                  do j=1,5
                     if (abs(cov(i,j)).le.zero) cycle
                     cef(igind)=cef(igind)+cov(i,j)*sens(2,i,ig)*sens(3,j,ig2)
                     ceg(igind)=ceg(igind)+cov(i,j)*sens(2,i,ig)*sens(4,j,ig2)
                     cfg(igind)=cfg(igind)+cov(i,j)*sens(3,i,ig)*sens(4,j,ig2)
                  enddo
               enddo
            enddo
         enddo

  170    continue
  160 continue

   !--end of do loops
   ifresr=1

   return
   end subroutine rpxlc0

   subroutine rpxlc12(nwscr,a,iest,ieed)
   !--------------------------------------------------------------------
   ! Resolved resonances: lru=1, lcomp=1 or 2
   !--------------------------------------------------------------------
   use mainio ! provides nsyso
   use endf ! provides endf routines and variables
   use physics ! provides amassn
   use util ! provides error,mess
   ! externals
   integer::nwscr,iest,ieed
   real(kr)::a(nwscr)
   ! internals
   integer::nb,nw,lru1,lrf1,lb,i,ig,ig2,il,ii2,ii1,ipp,j
   integer::loopm,loopn,loop,ns,nrs1,nsmax,ipos,l1,l2,l3
   integer::lb2,lldum,igind,ipara,itmp,ilnum,ind,ii
   integer::imess,inow,jj,lll,nrs
   integer,parameter::maxe=600000
   integer,parameter::mxnpar=7000
   integer,parameter::maxb=30000
   real(kr)::awri,aw,eres,ajres,backdt,backdt2,backdt3,ajres2
   real(kr)::eres2,gwidth,rho,rr,rr2,ser,per,tmp
   real(kr)::dap2
   real(kr)::b(maxb)
   real(kr)::sigr(maxe,5),sigp(maxe,5),gsig(4,2501)
   real(kr)::sens(4,mxnpar,2501)
   real(kr)::cov(mxnpar,mxnpar)
   real(kr)::pneorg(10000)
   real(kr)::time
   integer::llmat(5)
   character(60)::strng1,strng2
   real(kr),parameter::rc1=.123e0_kr
   real(kr),parameter::rc2=.08e0_kr
   real(kr),parameter::third=0.333333333e0_kr
   real(kr),parameter::zero=0

   !--general resolved resonance subsection formats (lcomp=1)
   !--compact resolved resonance subsection formats (lcomp=2)
   imess=0
   if (lcomp.eq.1) then
      call contio(nendf,0,0,a(lptr),nb,nw)
      awri=c1h
      nsrs=n1h
      nlrs=n2h
      if (nsrs.gt.0) then
         arat=awri/(awri+1)
         aw=amassn*awri
         ra=rc1*aw**third+rc2
         ral=ra
         apl=ap
      endif
      if (nsrs.le.0) go to 600
   else
      call rpxlc2(nwscr,cov,mxnpar,a)
   endif

   !--write MF=2 data on b array
   call rskiprp(nscr6,b,isrr,ier,maxb)
   call contio(nscr6,0,0,b,nb,nw)
   lru1=l1h
   lrf1=l2h
   if (lrf1.eq.7) call error('rpxlc12','cannot handle RML with isammy=0',' ')
   lb=7
   call contio(nscr6,0,0,b(lb),nb,nw)
   nls1=n1h
   if (lru.ne.lru1.or.lrf.ne.lrf1) then
      write(strng2,'(''lru/lrf(mf=32)='',i3,''/'',i3,&
        &''  vs.  lru/lrf(mf=2)='',i3,''/'',i3)')lru,lrf,lru1,lrf1
      call error('rpxlc12','different type of resonance for lcomp=1',strng2)
   endif
   lb=lb+6
   l2=lb

   do il=1,nls1
      itmp=l2
      call listio(nscr6,0,0,b(l2),nb,nw)
      l2=l2+nw
      if (nb.gt.0) then
         if ((l2+nb-1).gt.maxb) then
            write(strng2,'(''need '',i7,'' words but only '',&
              &i7,'' allocated'')') l2+nb-1,maxb
            call error('rpxlc12','b array storage exceeded',strng2)
         endif
         do while (nb.ne.0)
            call moreio(nscr6,0,0,b(l2),nb,nw)
            l2=l2+nw
         enddo
      endif

      ind=itmp
      if (lrf1.eq.3) then
         apl=b(ind+1)
         if (apl.eq.zero) apl=ap
      endif
      if (naps.eq.1) then
         ral=apl
         ra=ap
      endif
      ll=nint(b(ind+2))
      llmat(il)=ll
      nrs1=nint(b(ind+5))
      do nr=1,nrs1
         rho=cwaven*arat*sqrt(abs(b(ind+6*nr)))*ral
         call efacts(ll,rho,ser,per)
         b(l2)=ser
         b(l2+1)=per
         b(l2+2)=0
         l2=l2+3
      enddo
   enddo

   lb2=l2
   l3=lb2

   !--end of writing MF=2 data to b array

   if (lcomp.eq.1) then
      lptr=lptr+6
      nsmax=nsrs
   else
      nsmax=1
   endif

   !--loop over the number of "sections" of covariance matrix
   !--store that information in array "a"
   do ns=1,nsmax
      if (lcomp.eq.1) then
         l1=lptr
         call listio(nendf,0,0,a(l1),nb,nw)
         l1=l1+nw
         if (nb.gt.0) then
            if ((l1+nb-1).gt.nwscr) then
               write(strng2,'(''need '',i7,'' words but only '',&
                 &i7,'' allocated'')') l1+nb-1,nwscr
               call error('rpxlc12','a array storage exceeded',strng2)
            endif
            do while (nb.ne.0)
               call moreio(nendf,0,0,a(l1),nb,nw)
               l1=l1+nw
            enddo
         endif
         mpar=nint(a(lptr+2))
         nrb=nint(a(lptr+5))
         nvs1=6*nrb
         npar=mpar*nrb
         if (npar+1.gt.mxnpar) then
            write(strng2,'(''npar='',i8, '' +1 >  mxnpar='',i8)')&
              npar,mxnpar
            call error('rpxlc12','storage exceeded.', strng2)
         endif
      endif

      if (mpar.gt.4.and.lrf.le.2) then
         call error('rpxlc12','mpar.gt.4.and.lrf.le.2 not coded',' ')
      endif
      if (lrf.le.0.or.lrf.gt.3) then
         write(strng2,'(''lrf='',i3,'' is not coded.'')') lrf
         call error('rpxlc12','lcomp=1 general form.',strng2)
      endif

      !--check for presence or absence of scattering radius uncertainty
      !  data.  Advise user via mess and proceed.
      if (isr.eq.0.and.imess.eq.0) then
         imess=1
         write(strng1,'(''no scattering radius uncertainty'')')
         call mess('rpxlc12',strng1,' ')
      else if (isr.eq.1.and.imess.eq.0) then
         imess=1
         write(strng1,'(''include scattering radius uncertainty'')')
         call mess('rpxlc12',strng1,' ')
      endif

      if (isr.eq.1) then

         !--start with reference (no perturbation) calculation
         call rpendf(ii,99,0._kr,a,nwscr,sigr,-1._kr,b,maxb)

         !--now perturb the system, including
         !  - perturbation of scattering radius from mf=2
         !  - perturbation of penetration factor
         inow=1
         ap=b(inow+7)
         ap=ap+dap
         b(inow+7)=ap
         nls=nint(b(inow+10))
         inow=inow+12
         itmp=1
         do lll=1,nls
            apl=b(inow+1)
            if (apl.eq.0) then
               apl=ap
            else
               apl=apl+dap3(lll)
            endif
            b(inow+1)=apl
            ll=nint(b(inow+2))
            nrs=nint(b(inow+5))
            inow=inow+6
            do jj=1,nrs
               rho=cwaven*arat*sqrt(abs(b(inow+6*(jj-1))))*apl
               call efacts(ll,rho,ser,per)
               pneorg(itmp)=b(inow+6*nrs+3*(jj-1))
               pneorg(itmp+1)=b(inow+6*nrs+3*(jj-1)+1)
               itmp=itmp+2
               b(inow+6*nrs+3*(jj-1))=ser
               b(inow+6*nrs+3*(jj-1)+1)=per
            enddo
            inow=inow+6*nrs+3*nrs
         enddo
         call rpendf(ii,99,0._kr,a,nwscr,sigp,-1._kr,b,maxb)

         !--sensitivity calculation
         do ii1=1,4
            do ii2=1,ii
               tmp=((sigp(ii2,ii1)-sigr(ii2,ii1)))/dap
               sigp(ii2,ii1)=tmp
            enddo
         enddo
         call rpxgrp(ngn,egn,sigp,ii,gsig(1,1),a,nwscr)
         do ii1=1,ngn
            tmp=cflx(ii1)*abn
            do j=1,4
               gsig(j,ii1)=gsig(j,ii1)*tmp
           enddo
         enddo

         !--to get absolute standard deviation for scattering radius
         dap2=dap*dap

         !--error propagation
         igind=0
         do ig=1,ngn
            do ig2=ig,ngn
               igind=igind+1
               cff(igind)=cff(igind)+dap2*gsig(3,ig)*gsig(3,ig2)
               cgg(igind)=cgg(igind)+dap2*gsig(4,ig)*gsig(4,ig2)
               cee(igind)=cee(igind)+dap2*gsig(2,ig)*gsig(2,ig2)
               ctt(igind)=ctt(igind)+dap2*gsig(1,ig)*gsig(1,ig2)
            enddo
         enddo

         igind=0
         do ig=1,ngn
            do ig2=1,ngn
               igind=igind+1
               cef(igind)=cef(igind)+dap2*gsig(2,ig)*gsig(3,ig2)
               ceg(igind)=ceg(igind)+dap2*gsig(2,ig)*gsig(4,ig2)
               cfg(igind)=cfg(igind)+dap2*gsig(3,ig)*gsig(4,ig2)
            enddo
         enddo

         !--restore reference data
         inow=1
         ap=b(inow+7)
         ap=ap-dap
         b(inow+7)=ap
         nls=nint(b(inow+10))
         itmp=1
         inow=inow+12
         do lll=1,nls
            apl=b(inow+1)
            apl=apl/(1+dap3(lll))
            b(inow+1)=apl
            nrs=nint(b(inow+5))
            inow=inow+6
            do jj=1,nrs
               b(inow+6*nrs+3*(jj-1))=pneorg(itmp)
               b(inow+6*nrs+3*(jj-1)+1)=pneorg(itmp+1)
               itmp=itmp+2
            enddo
            inow=inow+6*nrs+3*nrs
         enddo
      !--end of scattering radius uncertainty treatment
      endif

      ipos=0
      loop=0
      do loopm=1,nrb
         do loopn=1,mpar
            loop=loop+1
            if (loopn.eq.1) then
               write (*,'(''Resonance number'',i5,''(/'',i5,&
                 &'')   Resonance energy'',1p,e15.7)') loopm,nrb,a(lptr+6*loopm)
            else
               write (*,'(''Resonance number'',i5,''(/'',i5,&
                 &'')   Width number'', i2)') loopm,nrb,loopn-1
            endif

            !--search aimed mf32 resonance in mf=2
            if (loopn.eq.1) then
               eres =a(lptr+6*loopm)
               ajres=a(lptr+6*loopm+1)
               il2=lb
               do il=1,nls1
                  itmp=il2
                  ipara=nint(b(il2+5))
                  if (ipara.ne.0) then
                     do ipp=1,ipara
                        il2=il2+6
                        eres2 =b(il2)
                        ajres2=b(il2+1)
                        if (eres*eres2.gt.zero) then
                           rr=abs(eres/eres2-1)
                           rr2=abs(ajres-ajres2)
                           if (rr.lt.1e-6_kr.and.rr2.lt.1e-4_kr) then
                              ipos=itmp+6+ipara*6+(ipp-1)*3
                              goto 461
                           endif
                        endif
                     enddo
                     il2=il2+6
                     il2=il2+ipara*3
                  endif
               enddo
               write(strng1,'(''E:'',1p,e12.4,'' ajres:'',e12.4)')
               call error('rpxlc12','problem',strng1)
 461           continue
               ajres=abs(ajres)
               ilnum=il
            endif

            !--perturbed(-)
            if (loopn.eq.1) then
               il3=il2
               backdt=b(il2)
               b(il2)=backdt*0.9999e0_kr
               gwidth=backdt*0.0001e0_kr
               backdt2=b(ipos)
               backdt3=b(ipos+1)
               rho=cwaven*arat*sqrt(abs(b(il2)))*ral
               lldum=llmat(il)
               call efacts(lldum,rho,ser,per)
               b(ipos)=ser
               b(ipos+1)=per
            else
               if (lrf.eq.1.or.lrf.eq.2) then
                  il3=il2+loopn+1
               else if (lrf.eq.3) then
                  il3=il2+loopn
               endif
               backdt=b(il3)
               gwidth=backdt*0.01e0_kr
               b(il3)=backdt*0.99e0_kr
            endif
            if (gwidth.ne.zero)&
              call rpendf(ii,ilnum,ajres,a,nwscr,sigr,eres,b,maxb)

            b(il3)=backdt
            if (loopn.eq.1) then
               b(ipos)=backdt2
               b(ipos+1)=backdt3
            endif

            !--perturbed(+)
            if (loopn.eq.1) then
               il3=il2
               b(il2)=backdt*1.0001e0_kr
               rho=cwaven*arat*sqrt(abs(b(il2)))*ral
               lldum=llmat(il)
               call efacts(lldum,rho,ser,per)
               b(ipos)=ser
               b(ipos+1)=per
            else
               b(il3)=backdt*1.01e0_kr
            endif
            if (gwidth.ne.zero)&
              call rpendf(ii,ilnum,ajres,a,nwscr,sigp,eres,b,maxb)

            b(il3)=backdt
            if (loopn.eq.1) then
               b(ipos)=backdt2
               b(ipos+1)=backdt3
            endif
            if (gwidth.ne.zero) then

               !--differencing
               do ii1=1,4
                  do ii2=1,ii
                    tmp=(sigp(ii2,ii1)-sigr(ii2,ii1))/(gwidth*2)
                    sigp(ii2,ii1)=tmp
                  enddo
               enddo

               !--integration
               call rpxgrp(ngn,egn,sigp,ii,gsig(1,1),a,nwscr)
               do ig=iest,ieed
                  tmp=cflx(ig)*abn
                  do j=1,4
                    sens(j,loop,ig)=gsig(j,ig)*tmp
                  enddo
               enddo
            else
               do ig=iest,ieed
                  do j=1,4
                    sens(j,loop,ig)=0
                  enddo
               enddo
            endif

         enddo !--end of do-loop over number of parameters per resonance
      enddo !--end of do-loop over number of resonances

      if (lcomp.eq.1) then
         l3=lptr+5+nvs1
         do i=1,npar
            do j=i,npar
               l3=l3+1
               tmp=a(l3)
               cov(i,j)=tmp
               cov(j,i)=tmp
            enddo
         enddo
      endif
      call timer(time)
      write(strng1,'("resonance parameter loop done",13x,f8.1,"s")')time
      call mess('rpxlc12',strng1,'')

      igind=0
      do ig=1,ieed
         do ig2=ig,ngn
            igind=igind+1
            if (ig.ge.iest.and.ig.le.ieed.and.ig2.ge.iest.and.ig2.le.ieed) then
               do i=1,npar
                  do j=i,npar
                    tmp=cov(i,j)
                    if (tmp.ne.zero) then
                       cff(igind)=cff(igind)+tmp*sens(3,i,ig)*sens(3,j,ig2)
                       cgg(igind)=cgg(igind)+tmp*sens(4,i,ig)*sens(4,j,ig2)
                       cee(igind)=cee(igind)+tmp*sens(2,i,ig)*sens(2,j,ig2)
                       ctt(igind)=ctt(igind)+tmp*sens(1,i,ig)*sens(1,j,ig2)
                       if (i.ne.j) then
                          cff(igind)=cff(igind)+tmp*sens(3,j,ig)*sens(3,i,ig2)
                          cgg(igind)=cgg(igind)+tmp*sens(4,j,ig)*sens(4,i,ig2)
                          cee(igind)=cee(igind)+tmp*sens(2,j,ig)*sens(2,i,ig2)
                          ctt(igind)=ctt(igind)+tmp*sens(1,j,ig)*sens(1,i,ig2)
                       endif
                    endif
                  enddo
               enddo
            endif
         enddo
         if (ig.gt.nresg) nresg=ig
      enddo
      call timer(time)
      write(strng1,'("sensitivity calculation continues",9x,f8.1,"s")')time
      call mess('rpxlc12',strng1,'')

      igind=0
      do ig=1,ieed
         do ig2=1,ngn
            igind=igind+1
            if (ig.ge.iest.and.ig.le.ieed.and.ig2.ge.iest.and.ig2.le.ieed) then
               do i=1,npar
                  do j=i,npar
                     tmp=cov(i,j)
                     if (tmp.ne.zero) then
                        cef(igind)=cef(igind)+tmp*sens(2,i,ig)*sens(3,j,ig2)
                        ceg(igind)=ceg(igind)+tmp*sens(2,i,ig)*sens(4,j,ig2)
                        cfg(igind)=cfg(igind)+tmp*sens(3,i,ig)*sens(4,j,ig2)
                        if (i.ne.j) then
                           cef(igind)=cef(igind)+tmp*sens(2,j,ig)*sens(3,i,ig2)
                           ceg(igind)=ceg(igind)+tmp*sens(2,j,ig)*sens(4,i,ig2)
                           cfg(igind)=cfg(igind)+tmp*sens(3,j,ig)*sens(4,i,ig2)
                        endif
                     endif
                  enddo
               enddo
            endif
         enddo
      enddo
      call timer(time)
      write(strng1,'("sensitivity calculation completed",9x,f8.1,"s")')time
      call mess('rpxlc12',strng1,'')

   enddo
   !--end of "sections of covariance matrix" from ENDF File32

  600 continue
   if (nlrs.gt.0) then
      do ns=1,nlrs
         call listio(nendf,0,0,a(lptr),nb,nw)
         lb=nint(a(lptr+3))
         l1=lptr+nw
         if (nb.gt.0) then
            if ((l1+nb-1).gt.nwscr) then
               write(strng2,'(''need '',i7,'' words but only '',&
                 &i7,'' allocated'')') l1+nb-1,nwscr
               call error('rpxlc12','a array for nlrs storage exceeded',strng2)
            endif
            do while (nb.ne.0)
               call moreio(nendf,0,0,a(l1),nb,nw)
               l1=l1+nw
            enddo
         endif
         call error('rpxlc12','nlrs>0 not coded.',' ')
      enddo
   endif
   ifresr=1
   return
   end subroutine rpxlc12

   subroutine rpxlc2(nwscr,cov,mxnpar,a)
   !--------------------------------------------------------------------
   ! Resolved resonances: lru=1, lcomp=2
   !--------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use physics ! provides ammasn
   use util ! provides error
   ! externals
   integer::nwscr,mxnpar
   real(kr)::cov(mxnpar,mxnpar),a(nwscr)
   ! internals
   integer::nb,nw,i1,i2,nind,lbg,l1,l2,l3,n1,n2,n3
   integer::nn1,nn2,nnn,nx,nn2p,nm,i,m,ii,mm,mmm,ndigit,mbase
   integer,dimension(6)::mpid
   integer,dimension(6),parameter::mpidbw=(/1,4,5,6,0,0/)
   integer,dimension(6),parameter::mpidrm=(/1,3,4,5,6,0/)
   real(kr)::aw,awri,fd,std1,std2
   character(60)::strng2
   real(kr),parameter::rc1=.123e0_kr
   real(kr),parameter::rc2=.08e0_kr
   real(kr),parameter::third=0.333333333e0_kr
   real(kr),parameter::half=0.5e0_kr
   real(kr),parameter::zero=0

   if (lrf.eq.1.or.lrf.eq.2) then
      mpid=mpidbw
   elseif (lrf.eq.3) then
      mpid=mpidrm
   else
      write(strng2,'(''not ready for lrf='',i1,'', lcomp=2'')')lrf
      call error('rpxlc2',strng2,'')
   endif
   do i1=1,mxnpar
      do i2=1,mxnpar
         cov(i1,i2)=0
      enddo
   enddo
   nind=1
   lbg=lptr
   l1=lptr
   call listio(nendf,0,0,a(l1),nb,nw)
   l1=l1+nw
   if (nb.gt.0) then
      if ((l1+nb-1).gt.nwscr) then
         write(strng2,'(''need '',i7,'' words but only '',&
           &i7,'' allocated'')') l1+nb-1,nwscr
         call error('rpxlc2','a array storage exceeded',strng2)
      endif
      do while (nb.ne.0)
         call moreio(nendf,0,0,a(l1),nb,nw)
         l1=l1+nw
      enddo
   endif
   awri=a(lptr)
   arat=awri/(awri+1)
   aw=amassn*awri
   ra=rc1*aw**third+rc2
   ral=ra
   apl=ap

   !--number of resonances (nrsa in endf manual)
   nrb=nint(a(lbg+5))
   l3=lbg

   !--pack uncertainties of selected parameters into the diagonal
   mm=0
   do n2=1,nrb
      l3=l3+12
      do m=1,5
        mm=mm+1
        ii=l3+mpid(m)-1
        cov(mm,mm)=a(ii)
      enddo
   enddo

   !--repack the resonance parameters, without uncertainties,
   !--into the work array.
   l3=lbg+6
   l2=lbg+6
   do n2=1,nrb
      do n3=1,6
         a(l3+n3-1)=a(l2+n3-1)
      enddo
      l3=l3+6
      l2=l2+12
   enddo
   lptr=l1

   !--read the correlation matrix in compact representation and
   !--expand to the full covariance matrix
   call contio(nendf,0,0,a(l3),nb,nw)
   ndigit=l1h
   if (ndigit.lt.2.or.ndigit.gt.6) then
     if (ndigit.eq.0) then
       ndigit=2
       call mess('rpxlc2','ndigit from file is zero'&
                 ,'reset to default, ndigit=2')
     else
       call error('rpxlc2','illegal value of ndigit',' ')
     endif
   endif
   nnn=l2h
   nm =n1h
   nx =n2h

   mpar=nnn/nrb
   if (mpar.ne.5) then
      mm=0
      do nr=1,nrb
         mbase=5*(nr-1)
         do m=1,mpar
            mm=mm+1
            mmm=mbase+m
            cov(mm,mm)=cov(mmm,mmm)
         enddo
      enddo
   endif

   if (nm.ne.0) then
      do n2=1,nm
         nw=ndigit
         fd=10**ndigit
         call intgio(nendf,0,0,a(l3),nb,nw)
         nn1=nint(a(l3  ))
         nn2=nint(a(l3+1))
         nn2p=nn2-1
         do n3=3,nw
            nn2p=nn2p+1
            if (nn2p.ge.nn1) exit
            std1=cov(nn2p,nn2p)
            std2=cov(nn1,nn1)
            l2=l3-1+n3
            if (nint(a(l2)).gt.0)then
               cov(nn2p,nn1)=((a(l2)+half)/fd)*std1*std2
            else if (nint(a(l2)).lt.0) then
               cov(nn2p,nn1)=(-(-a(l2)+half)/fd)*std1*std2
            endif
            cov(nn1,nn2p)=cov(nn2p,nn1)
         enddo
      enddo
   endif

   !--square the diagonal terms to get the variance
   npar=nrb*mpar
   do n1=1,npar
      cov(n1,n1)=cov(n1,n1)*cov(n1,n1)
   enddo
   lptr=lbg

   return
   end subroutine rpxlc2

   subroutine rpxunr(a,amur,mxlru2,iest,ieed,nwscr)
   !--------------------------------------------------------------------
   ! Unresolved resonance region (lru=2)
   !--------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use util ! provides error
   ! externals
   integer::mxlru2,iest,ieed,nwscr
   real(kr)::a(nwscr),amur(3,mxlru2)
   ! internals
   integer::nb,nw,l,l1,l2,l3,nl,ig,i,j,ig2,ii,igind,iscr,njs,inow
   integer::loop,loopn,l0
   integer,parameter::maxb=4000
   integer,parameter::mxnpar=100
   integer,parameter::maxe=600000
   real(kr)::e1,ebc,sfac,s,tmp,bb
   real(kr)::sig(maxe,5),sig1(4)
   real(kr)::gsigr(4,2501),gsigp(4,2501)
   real(kr)::sens(4,mxnpar,2501),cov(mxnpar,mxnpar)
   real(kr)::b(maxb)
   character(60)::strng2
   real(kr),parameter::zero=0

   write(*,*)'Unresolved resonance energy range.'
   l1=lptr
   do nl=1,nls
      call listio(nendf,0,0,a(l1),nb,nw)
      l1=l1+nw
      do while (nb.gt.0)
         call moreio(nendf,0,0,a(l1),nb,nw)
         l1=l1+nw
      enddo
   enddo
   call listio(nendf,0,0,a(l1),nb,nw)
   l2=l1+nw
   do while (nb.gt.0)
      call moreio(nendf,0,0,a(l2),nb,nw)
      l2=l2+nw
   enddo

   mpar=nint(a(l1+2))
   npar=nint(a(l1+5))
   if (npar+1.gt.mxnpar) then
      write(strng2,'(''npar='',i8,'' +1 >  mxnpar='',i8)') npar,mxnpar
      call error('rpxunr','storage exceeded (lru=2).',strng2)
   endif

   do ig=1,ngn
      do i=1,npar
         do j=1,4
            sens(j,i,ig)=0
         enddo
      enddo
   enddo

   njs=0
   inow=6
   l0=l1-lptr+6
   l2=l0+1

   do 850 loop=1,npar+1
      do i=1,l0
        b(i)=a(lptr+i-7)
      enddo
      if (loop.eq.1) go to 860
      loopn=mod(loop-1,mpar)
      if (loopn.eq.1.or.mpar.eq.1) then
         if (njs.eq.0) then
            njs=nint(b(inow+6))
         else
            njs=njs-1
         endif
         inow=inow+6
         l2=l2+1
         b(l2)=b(inow+1)
         b(inow+1)=b(inow+1)*1.01e0_kr
      else if (loopn.eq.2.or.(mpar.eq.2.and.loopn.eq.0)) then
         l2=l2+1
         b(l2)=b(inow+3)
         b(inow+3)=b(inow+3)*1.01e0_kr
      else if (loopn.eq.3.or.(mpar.eq.3.and.loopn.eq.0)) then
         l2=l2+1
         b(l2)=b(inow+4)
         b(inow+4)=b(inow+4)*1.01e0_kr
      else if (loopn.eq.4.or.(mpar.eq.4.and.loopn.eq.0)) then
         if (lfw.eq.1) then
            l2=l2+1
            b(l2)=b(inow+5)
            b(inow+5)=b(inow+5)*1.01e0_kr
         else if (lfw.eq.0) then
            l2=l2+1
            b(l2)=b(inow+6)
            b(inow+6)=b(inow+6)*1.01e0_kr
         endif
      else if (loopn.eq.5.or.(mpar.eq.5.and.loopn.eq.0)) then
         l2=l2+1
         b(l2)=b(inow+6)
         b(inow+6)=b(inow+6)*1.01e0_kr
      endif
      if (loopn.eq.0.and.njs.eq.1) then
         inow=inow+6
         njs=0
      endif
  860 continue
      e1=elg
      ii=0
      go to 880

  870 continue
      e1=e1*1.015e0_kr
      ebc=e1/1.015e0_kr
      if (ebc.lt.elr.and.e1.gt.elr) e1=elr
      if (ebc.lt.ehr.and.e1.gt.ehr) e1=ehr
  880 continue
      if (e1.gt.ehg) e1=ehg
      if (e1.ge.elr.and.e1.le.ehr) then
         call ggunr1(e1,sig1,b,amur,mxlru2)
      else
         do i=1,4
            sig1(i)=0
         enddo
      endif
      ii=ii+1
      if (ii.gt.maxe) call error('rpxunr',&
        'number of pointwise xsec of resonance exceeded.',&
        'please increase the maxe parameter.')
      do i=1,4
         sig(ii,i)=sig1(i)
      enddo
      sig(ii,5)=e1
      if (e1.ge.ehg) go to 890
      go to 870
  890 continue
      if (loop.eq.1) then
         call rpxgrp(ngn,egn,sig,ii,gsigr,a,nwscr)
      else
         call rpxgrp(ngn,egn,sig,ii,gsigp,a,nwscr)

         !--sensitivity calculation
         do ig=iest,ieed
            do j=1,4
               sfac=gsigr(j,ig)
               i=loop-1
               if (b(l0+1+i).ne.zero) then
                  s=gsigp(j,ig)-sfac
                  s=100*s/b(l0+1+i)
                  if (abs(s).ge.1.e-10_kr) then
                     sens(j,i,ig)=s*cflx(ig)*abn
                  endif
               endif
            enddo
         enddo

      endif
  850 continue

   l2=l0+1
   l3=l1+5
   do i=1,npar
      do j=i,npar
         l3=l3+1
         bb=b(l2+i)*b(l2+j)
         tmp=a(l3)*bb
         cov(i,j)=tmp
         cov(j,i)=tmp
      enddo
   enddo

   igind=0
   do ig=1,ngn
      do ig2=ig,ngn
         igind=igind+1
         if (ig.ge.iest.and.ig.le.ieed.and.&
            ig2.ge.iest.and.ig2.le.ieed) then
            do i=1,npar
               do j=1,npar
                  if (cov(i,j).eq.zero) cycle
                  uff(igind)=uff(igind)+cov(i,j)*sens(3,i,ig)*sens(3,j,ig2)
                  ugg(igind)=ugg(igind)+cov(i,j)*sens(4,i,ig)*sens(4,j,ig2)
                  uee(igind)=uee(igind)+cov(i,j)*sens(2,i,ig)*sens(2,j,ig2)
                  utt(igind)=utt(igind)+cov(i,j)*sens(1,i,ig)*sens(1,j,ig2)
               enddo
            enddo
         endif
         if (ig.gt.nresg) nresg=ig
      enddo
   enddo

   igind=0
   do ig =1,ngn
      do ig2=1,ngn
         igind=igind+1
         do i=1,npar
            if (sens(1,i,ig).eq.zero.and.sens(2,i,ig).eq.zero.and.&
              sens(3,i,ig).eq.zero.and.sens(4,i,ig).eq.zero) cycle
            do j=1,npar
               if (cov(i,j).eq.zero) cycle
               uef(igind)=uef(igind)+cov(i,j)*sens(2,i,ig)*sens(3,j,ig2)
               ueg(igind)=ueg(igind)+cov(i,j)*sens(2,i,ig)*sens(4,j,ig2)
               ufg(igind)=ufg(igind)+cov(i,j)*sens(3,i,ig)*sens(4,j,ig2)
            enddo
         enddo
      enddo
   enddo

   ifunrs=1
   write(*,*)'... ended.'

   return
   end subroutine rpxunr

   subroutine rpendf(ii,npnls,valspi,a,nwscr,sig,eres,b,maxb)
   !--------------------------------------------------------------------
   ! Calculation of pointwise cross section in Lcomp1 or 2
   !--------------------------------------------------------------------
   use util ! provides sigfig,error
   ! externals
   integer::ii,npnls,nwscr,maxb
   integer,parameter::maxe=600000
   real(kr)::valspi,a(nwscr),sig(maxe,5),eres,b(maxb)
   ! internals
   integer::i
   real(kr)::e1,elb1,elu1,elb2,elu2,elb3,elu3,e2,ebc
   real(kr)::sig1(4)
   real(kr),parameter::zero=0

   e1=elg
   ii=0

   elb1=0.9e0_kr*eres
   elu1=1.1e0_kr*eres
   elb2=0.8e0_kr*eres
   elu2=1.2e0_kr*eres
   elb3=0.7e0_kr*eres
   elu3=1.3e0_kr*eres

  520 continue
   e1=sigfig(e1,8,0)
   if (e1.gt.ehg) e1=ehg
   if (e1.ge.elr.and.e1.le.ehr) then
      if (lrf.eq.3) then
         call ggrmat(e1,sig1,b(1),npnls,valspi)
      else if (lrf.eq.2) then
         call ggmlbw(e1,sig1,b(1))
      endif
   else
      do i=1,4
         sig1(i)=0
      enddo
   endif

   ii=ii+1
   if (ii.gt.maxe) call error('rpendf',&
     'number of pointwise xsec of resonance exceeded.',&
     'please increase the maxe parameter.')

   do i=1,4
      sig(ii,i)=sig1(i)
   enddo
   sig(ii,5)=e1

   if (e1.lt.ehg) then
      if (eres.lt.zero) then
         e2=1+5*(eskip1-1)
      else
         if (e1.lt.0.1e0_kr) then
            e2=eskip4
         else if (e1.gt.elb1.and.e1.lt.elu1) then
            e2=eskip1
         else if (e1.gt.elb2.and.e1.lt.elu2) then
            e2=eskip2
         else if (e1.gt.elb3.and.e1.lt.elu3) then
            e2=eskip3
         else
            e2=1.02e0_kr
         endif
      endif
      ebc=e1
      e1=e1*e2
      if (ebc.lt.elr.and.e1.gt.elr) e1=elr
      if (ebc.lt.ehr.and.e1.gt.ehr) e1=ehr
      go to 520
   endif

   return
   end subroutine rpendf

   subroutine rdumrd2(matd,nendf,nscr6,a,amur,mxlru2,nwscr)
   !--------------------------------------------------------------------
   ! dummy read the resonance parameters (mf=2)
   !--------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use util ! provides error
   ! externals
   integer::matd,nendf,nscr6,mxlru2,nwscr
   real(kr)::a(nwscr),amur(3,mxlru2)
   ! internals
   integer::indx
   integer::nb,nw,nlru2,ni,ner,ne,nro,nl,njs,nj
   character(60)::strng

   !--read mf=2 (resonance parameters)
   call findf(matd,2,151,nendf)
   call contio(nendf,nscr6,0,a,nb,nw)
   nis=n1h
   nlru2=0
   indx=0

   !--loop over all isotopes
   do ni=1,nis
      call contio(nendf,nscr6,0,a,nb,nw)
      lfw=l2h
      ner=n1h

      !--loop over all energy ranges
      do ne=1,ner
         indx=indx+1
         call contio(nendf,nscr6,0,a,nb,nw)
         lru=l1h
         lrf=l2h
         nro=n1h
         if (nro.gt.0) call tab1io(nendf,nscr6,0,a,nb,nw)

         !--breit-wigner
         if (lru.eq.1.and.(lrf.eq.1.or.lrf.eq.2)) then
            call contio(nendf,nscr6,0,a,nb,nw)
            nls=n1h
            nlspepi(indx)=nls
            do nl=1,nls
               call listio(nendf,nscr6,0,a,nb,nw)
               do while (nb.gt.0)
                  call moreio(nendf,nscr6,0,a,nb,nw)
               enddo
            enddo

         !--reich-moore
         else if (lru.eq.1.and.lrf.eq.3) then
            call contio(nendf,nscr6,0,a,nb,nw)
            nls=n1h
            nlspepi(indx)=nls
            do nl=1,nls
               call listio(nendf,nscr6,0,a,nb,nw)
               do while (nb.gt.0)
                  call moreio(nendf,nscr6,0,a,nb,nw)
               enddo
            enddo

         !--reich-moore limited
         else if (lru.eq.1.and.lrf.eq.7) then
            call contio(nendf,nscr6,0,a,nb,nw)
            njs=n1h
            call listio(nendf,nscr6,0,a,nb,nw)
            do nj=1,njs
               call listio(nendf,nscr6,0,a,nb,nw)
               do while (nb.gt.0)
                  call moreio(nendf,nscr6,0,a,nb,nw)
               enddo
               call listio(nendf,nscr6,0,a,nb,nw)
               do while (nb.gt.0)
                  call moreio(nendf,nscr6,0,a,nb,nw)
               enddo
            enddo

         !--unresolved resonance (lrf=1,lfw=0)
         else if (lru.eq.2.and.lrf.eq.1.and.lfw.eq.0) then
            call contio(nendf,nscr6,0,a,nb,nw)
            nls=n1h
            nlspepi(indx)=nls
            do nl=1,nls
               call listio(nendf,nscr6,0,a,nb,nw)
               do while (nb.gt.0)
                  call moreio(nendf,nscr6,0,a,nb,nw)
               enddo
            enddo

         !--unresolved resonance (lrf=1,lfw=1)
         else if (lru.eq.2.and.lrf.eq.1.and.lfw.eq.1) then
            call listio(nendf,nscr6,0,a,nb,nw)
            nls=nint(a(6))
            nlspepi(indx)=nls
            do nl=1,nls
               call contio(nendf,nscr6,0,a,nb,nw)
               njs=n1h
               do nj=1,njs
                  call listio(nendf,nscr6,0,a,nb,nw)
                  do while (nb.gt.0)
                     call moreio(nendf,nscr6,0,a,nb,nw)
                  enddo
               enddo
            enddo

         !--unresolved resonance (lrf=2)
         else if (lru.eq.2.and.lrf.eq.2) then
            call contio(nendf,nscr6,0,a,nb,nw)
            nls=n1h
            nlspepi(indx)=nls
            do nl=1,nls
               call contio(nendf,nscr6,0,a,nb,nw)
               njs=n1h
               do nj=1,njs
                  call listio(nendf,nscr6,0,a,nb,nw)
                  nlru2=nlru2+1
                  amur(1,nlru2)=a(10)
                  amur(2,nlru2)=a(12)
                  amur(3,nlru2)=a(9)
                  do while (nb.gt.0)
                     call moreio(nendf,nscr6,0,a,nb,nw)
                  enddo
               enddo
            enddo
         else
            write(strng,'('' *** lru='',i3,''  lrf='',i3,&
              &'' no coding.'')') lru,lrf
            call error('dumrd2',strng,' ')
         endif
         if (nlru2.gt.mxlru2) call error('dumrd2',&
            'nlru2 was exceeded mxlru2',' ')
      enddo
   enddo

   return
   end subroutine rdumrd2

   subroutine rpxgrp(igx,egn,sig,ipoint,gsig,a,nwscr)
   !--------------------------------------------------------------------
   ! Convert pointwise cross sections to simplistic groupwise ones
   !--------------------------------------------------------------------
   use util ! provides error
   ! externals
   integer::igx,ipoint,nwscr
   integer,parameter::maxe=600000
   real(kr)::egn(*),sig(maxe,5),gsig(4,2501),a(nwscr)
   !internals
   integer::k,i,i0,ig,lord,idis,j
   real(kr)::sumde,x1,enext,wt1,wt2,de,ebb,ebx,coef,egnt,egnt1
   real(kr)::wt12,wt3,wtr,x2,xr,y1,xx,xl,x12,y2,y3,yr,yl,z1,z2,z3,zr,zl,wtl
   real(kr),parameter::half=0.5e0_kr
   real(kr),parameter::two=2
   real(kr),parameter::three=3
   real(kr),parameter::six=6
   real(kr),parameter::zero=0

   do k=1,igx
      do i=2,4
         gsig(i,k)=0
      enddo
   enddo

   sumde=0
   lord=0
   i0=1
  100 continue
   do ig=1,igx
      if (sig(i0,5).ge.egn(ig).and.sig(i0,5).lt.egn(ig+1)) go to 110
   enddo
   if (i0.lt.ipoint) then
      i0=i0+1
      go to 100
   else
      call error('rxgrpg','i0>ipoint not coded.',' ')
   endif
  110 continue

   x1=0
   call egtwtf(x1,enext,idis,lord,wt1)
   x1=sig(i0,5)
   call egtwtf(x1,enext,idis,lord,wt1)

   !--loop over all pointwise cross sections in a range
   do i=i0+1,ipoint

      wt2=wt1
      x2=x1
      x1=sig(i,5)
      x12=x1-x2
      egnt=egn(ig)
      egnt1=egn(ig+1)
      call egtwtf(x1,enext,idis,lord,wt1)
      if (x1.ge.egnt.and.x1.le.egnt1) then
         de=half*(x1-x2)*(wt1+wt2)
         sumde=sumde+de
         z1=(two*wt1+wt2)*x12/six
         z2=(two*wt2+wt1)*x12/six
         do j=2,4
            y1=sig(i,j)
            y2=sig(i-1,j)
            if (y1.ne.zero.or.y2.ne.zero) then
               xx=y1*z1+y2*z2
               gsig(j,ig)=gsig(j,ig)+xx
            endif
         enddo
         if (x1.eq.egnt1) then
           do j=2,4
             gsig(j,ig)=gsig(j,ig)/sumde
           enddo
           ig=ig+1
           sumde=0
         endif
      else if (x1.gt.egnt1) then
        ebb=egnt1
        ebx=ebb-x2
        wt12=wt1-wt2
        coef=ebx/x12
        wt3=wt2+wt12*ebx/x12
        de=ebx*(wt2+wt3)/2
        sumde=sumde+de
        z2=(two*wt2+wt3)*ebx/six
        z3=(two*wt3+wt2)*ebx/six
        do j=2,4
           y1=sig(i,j)
           y2=sig(i-1,j)
           y3=y2+(y1-y2)*coef
           if (y2.ne.zero.or.y3.ne.zero) then
              xx=y2*z2+y3*z3
              gsig(j,ig)=gsig(j,ig)+xx
           endif
        enddo
        do j=2,4
           gsig(j,ig)=gsig(j,ig)/sumde
        enddo
 1000  continue
        ig=ig+1
        if (x1.gt.egn(ig+1)) then
           xl=egn(ig)
           xr=egn(ig+1)
           ebx=xr-xl
           wtl=wt2+wt12*(xl-x2)/x12
           wtr=wt2+wt12*(xr-x2)/x12
           sumde=ebx*(wtl+wtr)/2
           zl=(2*wtl+wtr)*ebx/six
           zr=(2*wtr+wtl)*ebx/six
           do j=2,4
              yl=y2+(y1-y2)*(xl-x2)/(x1-x2)
              yr=y2+(y1-y2)*(xr-x2)/(x1-x2)
              if (yl.ne.zero.or.yr.ne.zero) then
                 gsig(j,ig)=yl*zl+yr*zr/sumde
              endif
           enddo
           goto 1000
        else
           ebx=x1-egn(ig)
           de=ebx*(wt3+wt1)/2
           sumde=de
           z1=(two*wt3+wt1)*ebx/six
           z3=(two*wt1+wt3)*ebx/six
           do j=2,4
              y1=sig(i,j)
              y2=sig(i-1,j)
              y3=y2+(y1-y2)*coef
              if (y3.ne.zero.or.y1.ne.zero) then
                 xx=y1*z1+y3*z3
                 gsig(j,ig)=gsig(j,ig)+xx
              endif
           enddo
        endif
      endif
   enddo

   !--calculate group total cross section
   do k=1,igx
      gsig(1,k)=gsig(2,k)+gsig(3,k)+gsig(4,k)
   enddo

   return
   end subroutine rpxgrp

   subroutine rskiprp(iu,adim,ni,ne,maxad)
   !--------------------------------------------------------------------
   ! Skip to the resonance parameters of subsection requested by ni & ne
   !--------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use util ! provides repoz,error
   ! externals
   integer::iu,ni,ne,maxad
   real(kr)::adim(maxad)
   ! internals
   integer::nis1,ni1,ner1,ne1,lrf1,nl1,nb,nw
   character(60)::strng

   call repoz(iu)
   call contio(iu,0,0,adim,nb,nw)
   nis1=n1h
   do ni1=1,nis1
      call contio(iu,0,0,adim,nb,nw)
      ner1=n1h
      do ne1=1,ner1
         if (ni1.eq.ni.and.ne1.eq.ne) go to 130
         call contio(iu,0,0,adim,nb,nw)
         lrf1=l2h
         if (lrf1.ge.4) then
            write(strng,'(''lrf='',i3)') lrf1
            call error('skiprp','no coding type',strng)
         endif
         call contio(iu,0,0,adim,nb,nw)
         nls1=n1h
         do nl1=1,nls1
            call listio(iu,0,0,adim,nb,nw)
            do while (nb.gt.0)
               call moreio(iu,0,0,adim,nb,nw)
            enddo
         enddo
      enddo
   enddo
  130 continue

   return
   end subroutine rskiprp

   subroutine grpav4(mprint)
   !--------------------------------------------------------------------
   ! compute multigroup legendre coefficients for reaction needed in
   ! the calculation of the covariance matrices.  calculation uses the
   ! union of the user specified group structure and the energy
   ! grid found in mfcov.
   !
   ! if the user has supplied a gendf tape (nscr5.ne.0) with mf3/mt251
   ! it will be collapsed to the union grid in lieu of processing the
   ! endf tape.
   !
   ! if present (i.e., nscr5.ne.0) then the gendf tape is always
   ! processed, even if the conditions allowing endf tape processing
   ! are true.
   !
   ! multigroup legendre coefficients can be computed from the endf
   ! tape if (i) mf4 ltt=1 or (ii) mf4 ltt=3 and the transition energy
   ! from legendre coefficients to probability distributions is .le.
   ! maximum union grid energy.  Also, require mf4_lct=mf34_lct=1, or
   ! mf4_lct=2 and mf34_lct=0 or 2.
   !
   ! the code currently processes elastic scattering p1 data only; if
   ! other mf34 mt's are defined we simply skip over them.
   !--------------------------------------------------------------------
   use mainio ! provide nsyso, nsyse
   use endf ! provides endf routines and variables
   use util ! provide openz,repoz,closz,timer
   ! externals
   integer::mprint
   ! internals
   integer::nwds,nb,ne,nw,i,il,ig2lo,ig,imt,iz,j,ng2,idis,m
   integer::ngrp,n,n1,n2,n3,n4,nz,ntw,nng,nngr,ngg,nmu
   real(kr)::sec,etop,ehi,e,enext,elo,thresh,time
   real(kr)::dele,emu,els
   real(kr)::ans(10,2),z(26),flux(10,10),al(10)
   real(kr),dimension(:),allocatable::scr
   real(kr),dimension(:),allocatable::ge,gmu,gels,gscr
   character*4::tz(20)
   equivalence(tz,z)
   character(66)::text
   character(60)::strng
   real(kr),parameter::eps=1.e-9_kr
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::elow=1.e-5_kr
   real(kr),parameter::zero=0
   real(kr),parameter::oneeps=0.99999999_kr

   !--initialize
   if (iread.eq.2) call error('grpav4',&
     'not coded for multimaterial group averaging.',' ')
   call timer(sec)
   write(nsyso,'(/,'' computing multigroup legendre coef.'',&
     &33x,f8.1,''s'')') sec
   call egnwtf
   nwds=nunion+10
   if ((npage+50).gt.nwds) nwds=npage+50
   allocate(scr(nwds))
   nscr4=14
   if (nendf.lt.0) nscr4=-nscr4
   call openz(nscr4,1)
   call repoz(nscr4)
   math=0
   mfh=0
   mth=0
   text=' '
   nw=17
   read(text,'(16a4,a2)') (tz(i),i=1,nw)
   call tpidio(0,nscr4,0,z,nb,nw)
   if (abs(egn(1)-elow).le.eps) egn(1)=elow
   etop=un(1+nunion)

   !--search for desired mat on nendf tape
   call repoz(nendf)
   call findf(matd,1,0,nendf)
   call contio(nendf,0,0,scr,nb,nw)
   za=c1h
   awr=c2h

   !--main loop over reactions
   math=matd
   mfh=4
   il=10
   if (nscr5.ne.0) il=1
   iz=1
   nw=legord*iz+1
   do 300 imt=1,nga
      math=matd
      mfh=4
      mtd=iga(imt)
      if (mtd.ne.2) then
         write(strng,'("skipping over mf=3, mt=",i3)')mtd
         call mess('grpav4',strng,'')
         cycle
      endif
      call timer(time)

      mth=mtd
      z(1)=za
      z(2)=awr
      z(3)=1
      z(4)=legord
      z(5)=0
      z(6)=nunion
      call contio(0,nscr4,0,z,nb,nwds)

      !--initialize
      ng2=2

      if (nscr5.eq.0) then
         e=0
         call egtlgc(e,thresh,idis,al)
         if (thresh.gt.etop) go to 270
         call egtflx(e,enext,idis,flux,il,iz)
         if (mprint.ne.0) then
            write(nsyso,'(/,'' legendre group constants: pl-order 1 to '',&
              &i2,26x,f8.1,''s'')') legord,time
            if (nsyse.gt.0) write(nsyse,'(1x,''pl='',i2,3x,f8.1,''s'')')&
              legord,time
            write(nsyso,'(&
              &'' u(1,1) element of transformation matrix (cm -> lab)='',&
              &1pe12.5)') u1lele(2)
            write(nsyso,'('' for mf'',i2,'' and mt'',i3,a)')&
               mfd,mtd,'  (same as mf=3/mt=251)'
            if (nsyse.gt.0) write(nsyse,'('' for mf'',i2,'' and mt'',i3,a)')&
              mfd,mtd,'  (same as mf=3/mt=251)'
            write(nsyso,'(5x,''group'',5x,''legendre constant'')')
         endif
      else
         if (iga(imt).eq.2) then
            m=251
         else
            m=iga(imt)
         endif
         write(strng,'("collapsing NGOUT mf=3, mt=",i3,&
                      &" to the union grid")')m
         call mess('grpav4',strng,'')
         call repoz(nscr5)
         call findf(math,1,0,nscr5)
         if (allocated(gscr)) deallocate(gscr)
         allocate(gscr(npage+6))
         call contio(nscr5,0,0,gscr,nb,nwds)
         nz=nint(gscr(4))
         ntw=nint(gscr(6))
         call listio(nscr5,0,0,gscr,nb,nwds)
         nng=nint(gscr(3))
         ngg=nint(gscr(4))
         if (allocated(ge)) deallocate(ge)
         allocate(ge(nng+1))
         allocate(gmu(nng))
         allocate(gels(nng))
         gmu(1:nng)=0
         gels(1:nng)=0
         n1=7+nz+ntw
         n2=n1+nng
         n3=1
         n4=n3+nng
         if (n2.gt.npage+6) then
             n2=npage+6
             n4=npage-2
         endif
         ge(n3:n4)=gscr(n1:n2)
         do while (nb.ne.0)
            n3=n4+1
            call moreio(nscr5,0,0,gscr,nb,nwds)
            if (n3.gt.nng+1) exit
            if ((n3+nwds-1).lt.nng) then
               n4=n3+nwds-1
               ge(n3:n4)=gscr(1:nwds)
            else
               n4=nng+1
               ge(n3:n4)=gscr(1:(n4-n3+1))
            endif
         enddo
         call findf(math,3,m,nscr5)
         call contio(nscr5,0,0,gscr,nb,nwds)
         nmu=nint(gscr(6))
         do n=1,nmu
            call listio(nscr5,0,0,gscr,nb,nwds)
            gmu(nint(gscr(6)))=gscr(8)
            gels(nint(gscr(6)))=gscr(9)
         enddo
      endif

      !--loop over union energy groups
  200 continue
      n=1
      mfd=3
      mtd=iga(imt)
      do 260 ig=1,nunion
         elo=un(ig)
         ehi=un(1+ig)
         ig2lo=0
         do j=1,2
            do i=1,il
               ans(i,j)=0
            enddo
         enddo
         if (nscr5.eq.0) then
         !--if no NGOUT, get Legendre data from the endf tape
            enext=ehi
            do
               call epanel(elo,enext,ans,il,iz,ig2lo,34)
               if (abs(enext/ehi-1.e0_kr).lt.eps) exit
               elo=enext
               enext=ehi
            enddo
            !--write this group on nscr4 tape
            do i=1,9
               ans(i,2)=ans(i,2)/ans(1,1)
               plele(ig,i)=ans(i,2)
               ! legendre coefficient a_i in center mass system
            enddo
            if (lct4.eq.2) ans(1,2)=ans(1,2)*u1lele(2)
            if (mprint.ne.0) write(nsyso,'(4x,i4,5x,1p,6e14.6:/&
                              &(13x,6e14.6))')ig,(ans(i,2),i=1,legord)
            if (ig.ne.nunion) then
               do i=1,legord
                  if (ans(i,2).ne.zero) go to 250
               enddo
               go to 260
            endif
         else
         !--collapse NGOUT (now nscr5) Legendre data to the union grid
         !  - we already know ge(1).le.un(1) & ge(nng+1).ge.un(nunion)
         !  - n is a group number, bounded by ge(n) to ge(n+1)
            do while (ge(n+1)/elo.lt.oneeps)
               n=n+1
            enddo
            !-- now know ge(n)<=elo<ge(n+1)
            dele=min(ge(n+1),ehi)-elo
            emu=dele*gels(n)*gmu(n)
            els=dele*gels(n)
            !--test if more GROUPR groups to merge into the union grid
            if (ge(n+1)/ehi.lt.oneeps) then
               n=n+1
               do while(ge(n+1)/ehi.lt.oneeps)
                  !--this group is within the current union group
                  dele=ge(n+1)-ge(n)
                  emu=emu+dele*gels(n)*gmu(n)
                  els=els+dele*gels(n)
                  n=n+1
               enddo
               !-- now know ehi<=ge(n+1)
               dele=ehi-ge(n)
               emu=emu+dele*gels(n)*gmu(n)
               els=els+dele*gels(n)
            endif
            ans(1,1)=ehi-elo
            ans(1,2)=emu/els
         endif
  250    mfh=mfd
         mth=mtd
         z(1)=0
         z(2)=0
         z(3)=ng2
         z(4)=ig2lo
         z(5)=nw
         z(6)=ig
         z(7)=ans(1,1)
         do i=1,legord
            z(i+7)=ans(i,2)
         enddo
         nwds=legord+7
         call listio(0,nscr4,0,z,nb,nwds)
  260 continue
      call asend(nscr4,0)
      go to 280
      ! write message if mt has threshold gt highest union energy
  270 continue
      write(strng,'(''mf '',i2,'' mt '',i3)') mfd,mtd
      call mess('grpav4',strng,'has threshold gt highest union energy.')
  280 continue
  300 continue

   !--grpav4 is finished.
   call afend(nscr4,0)
   call amend(nscr4,0)
   call atend(nscr4,0)
   deallocate(scr)
   if (allocated(iga))  deallocate(iga)
   if (allocated(ge))   deallocate(ge)
   if (allocated(gmu))  deallocate(gmu)
   if (allocated(gels)) deallocate(gels)
   call timer(sec)
   write(nsyso,'(/,'' legendre group averaging completed'',34x,f8.1,''s'',/)')&
     sec
   if (nsyse.gt.0) write(nsyse,'(/,'' legendre group averaging completed'',34x,&
     &f8.1,''s'',/)') sec
   return
   end subroutine grpav4

   subroutine alsigc(ncg,ncm,nun,alsig,clflx,b,egt,flux,sig,alp,ld,ld1,mt1,mt2)
   !--------------------------------------------------------------------
   ! Calculate the coarse group legendre*sigma
   !--------------------------------------------------------------------
   use util ! provides sigfig
   ! externals
   integer::ncg,ncm,nun,ld,ld1,mt1,mt2
   real(kr)::alsig(ncg,2),clflx(ncg,2),b(*),egt(nun+1),flux(nun),&
      sig(nun+1),alp(nun)
   ! internals
   integer::i,ig,jg,izero
   integer::mfinit=0

   !--initialize
   izero=0
   if (nlump.gt.0) call error('alsigc','no coded lump xsec.',' ')
   if (mfinit.eq.0) then
      mfinit=3
      do i=1,nunion+1
         egt(i)=sigfig(egt(i),ndig,0)
      enddo
      do i=1,ngn+1
         egn(i)=sigfig(egn(i),ndig,0)
      enddo
   endif

   !--compute cross-group legendre*sigma*flux and sigma*flux
   call rdsig(matd,mt1,izero,b,sig)
   call rdlgnd(nscr4,matd,mt1,ld,b,alp)
   do ig=1,ngn
      alsig(ig,1)=0
      clflx(ig,1)=0
      do jg=1,nunion
         if (egt(jg).lt.egn(ig).or.egt(jg).ge.egn(ig+1)) cycle
         alsig(ig,1)=alsig(ig,1)+alp(jg)*sig(jg)*flux(jg)
         clflx(ig,1)=clflx(ig,1)+sig(jg)*flux(jg)
      enddo
      alsig(ig,1)=alsig(ig,1)/clflx(ig,1)
   enddo
   call rdsig(matd,mt2,izero,b,sig)
   call rdlgnd(nscr4,matd,mt2,ld1,b,alp)
   do ig=1,ngn
      alsig(ig,2)=0
      clflx(ig,2)=0
      do jg=1,nunion
         if (egt(jg).lt.egn(ig).or.egt(jg).ge.egn(ig+1)) cycle
         alsig(ig,2)=alsig(ig,2)+alp(jg)*sig(jg)*flux(jg)
         clflx(ig,2)=clflx(ig,2)+sig(jg)*flux(jg)
      enddo
      alsig(ig,2)=alsig(ig,2)/clflx(ig,2)
   enddo

   return
   end subroutine alsigc

   subroutine egtlgc(e,enext,idis,al)
   !--------------------------------------------------------------------
   ! retrieve the legendre coefficient defined by mfd and mtd.
   ! initialize if e=0.
   !--------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use util ! provides error
   ! externals
   integer::idis
   real(kr)::e,enext,al(10)
   ! internals
   integer::nsig,nb,nw,mf,mt,lvt,ltt,ialv1,ij,nr,ne,i,intl,n,nl1,nl2,nm,nne,nnr
   real(kr)::aww,e1,e2
   integer,parameter::maxleg=64
   real(kr)::al1(maxleg),al2(maxleg)
   real(kr)::alnr(46),alg(70)
   real(kr),dimension(:),allocatable::alv
   real(kr),parameter::zero=0
   save nsig,intl,nnr,ne,nne,nl1,nl2,e1,e2,al1,al2

   !--initialize
   idis=0
   if (e.gt.zero) go to 200
   nsig=nendf
   mf=4
   mt=mtd
   call findf(matd,mf,mt,nsig)
   call contio(nsig,0,0,alnr,nb,nw)
   aww=c2h
   lvt=l1h
   ltt=l2h
   if (lvt.eq.1) then
      nw=(maxleg+1)**2+6
      allocate(alv(nw))
      call listio(nsig,0,0,alv,nb,nw)
      ialv1=1
      do while (nb.gt.0)
         ialv1=ialv1+nw
         call moreio(nsig,0,0,alv(ialv1),nb,nw)
      enddo
      nm=nint(alv(6))
      do ij=1,10
         u1lele(ij)=alv(6+(nm+1)+ij)
      enddo
      deallocate(alv)
   else
      call contio(nsig,0,0,alnr,nb,nw)
      call matrixin(aww,u1lele)
   endif
   if (ltt.eq.2) call error('egtlgc','not coded for ltt=2.',' ')
   call tab2io(nsig,0,0,alnr,nb,nw)
   nr=n1h
   ne=n2h
   do i=1,nr
      nbt(i)=nint(alnr(i*2+5))
      jnt(i)=nint(alnr(i*2+6))
   enddo
   intl=jnt(1)
   call listio(nsig,0,0,alg,nb,nw)
   e1=c2h
   nl1=n1h
   do i=1,nl1
      al1(i)=alg(6+i)
   enddo
   if (nl1.lt.maxleg) then
      do i=nl1+1,maxleg
         al1(i)=0
      enddo
   endif
   call listio(nsig,0,0,alg,nb,nw)
   e2=c2h
   nl2=n1h
   do i=1,nl2
      al2(i)=alg(6+i)
   enddo
   if (nl2.lt.maxleg) then
      do i=nl2+1,maxleg
         al2(i)=0
      enddo
   endif
   nnr=1
   nne=2
   enext=e2
   return

   !--retrieve legendre coefficient
  200 continue
   do i=1,10
      al(i)=0
   enddo
   if (e.ge.e2) then
      if (nne.eq.ne.and.e.le.e2*1.00001e0_kr) go to 300
      if (nne.ge.ne) go to 400
      do i=1,nl2
         al1(i)=al2(i)
      enddo
      if (nl2.lt.maxleg) then
         do i=nl2+1,maxleg
            al1(i)=0
         enddo
      endif
      nl1=nl2
      e1=e2
      call listio(nsig,0,0,alg,nb,nw)
      e2=c2h
      nl2=n1h
      do i=1,nl2
         al2(i)=alg(6+i)
      enddo
      if (nl2.lt.maxleg) then
         do i=nl2+1,maxleg
            al2(i)=0
         enddo
      endif
      nne=nne+1
   endif

  300 continue
   n=max(nl1,nl2)
   if (n.gt.10) n=10
   if (nne.gt.nbt(nnr)) then
      nnr=nnr+1
      intl=jnt(nnr)
   endif
   do i=1,n
      call terp1(e1,al1(i),e2,al2(i),e,al(i),intl)
   enddo

  400 continue
   enext=e2

   return
   end subroutine egtlgc

   subroutine musigc(ncg,ncm,nun,csig,cflx,b,egt,flux,sig,alp)
   !--------------------------------------------------------------------
   ! Calculate the coarse group mubar
   !--------------------------------------------------------------------
   use mainio ! provides nsyso
   use endf ! provides endf routines and variables
   use util ! provides sigfig
   ! externals
   integer::ncg,ncm,nun
   real(kr)::csig(ncg,ncm),cflx(ncg),b(*),egt(nun+1)
   real(kr)::flux(nun),sig(nun+1),alp(nun)
   ! internal
   integer::i,mat,mf,mt,nun1,ig,jg,ij,ibase,ip,ld,nb,nw,ngn1,ix,loc,np,nwds
   integer::ngnp1,izero
   real(kr)::abit,sss0
   real(kr)::c(6)
   character(2)::hmt='mt'
   character(5)::uline='-----'

   izero=0
   !--put the coarse group structure on nout, ala groupr
   if (nout.ne.0) then
      math=matd
      mfh=1
      mth=451
      b(1)=za
      b(2)=awr
      ! pass iverf on to covr
      b(3)=iverf
      b(4)=0
      b(5)=-11
      b(6)=0
      call contio(0,nout,0,b,nb,nw)
      b(1)=tempin
      b(2)=0
      b(3)=ngn
      nw=6
      ngnp1=ngn+1
      do i=1,ngnp1
         nw=nw+1
         b(nw)=egn(i)
      enddo
      np=nw-6
      b(5)=np
      loc=1
      call listio(0,nout,0,b(loc),nb,nw)
      do while (nb.gt.0)
        loc=loc+nw
        call moreio(0,nout,0,b(loc),nb,nw)
      enddo
     call asend(nout,0)
     call afend(nout,0)
   endif

   !--initialize
   mfd=3
   nun1=nunion+1
   do i=1,nun1
      egt(i)=sigfig(egt(i),ndig,0)
   enddo
   ngn1=ngn+1
   do i=1,ngn1
      egn(i)=sigfig(egn(i),ndig,0)
   enddo

   !--calculate coarse group flux
   do ig=1,ngn
      cflx(ig)=0
      do jg=1,nunion
         if (egt(jg).lt.egn(ig).or.egt(jg).ge.egn(ig+1)) cycle
         cflx(ig)=cflx(ig)+flux(jg)
      enddo
   enddo

   !--mt=251, mubar: average cosine of the scattering angle
   !--(laboratory system) for elastic scattering.
   !--compute cross-group cross sections and write on output tape.
   math=matd
   mfh=3
   mth=2
   call rdsig(math,mth,izero,b,sig)
   ld=1
   call rdlgnd(nscr4,math,mth,ld,b,alp)
   ix=1
   do ig=1,ngn
      csig(ig,ix)=0
      abit=0
      do jg=1,nunion
         if (egt(jg).lt.egn(ig).or.egt(jg).ge.egn(ig+1)) cycle
         sss0=0
         do ij=2,9
            sss0=sss0+plele(jg,ij)*u1lele(ij+1)
         enddo
         csig(ig,ix)=csig(ig,ix)+sig(jg)*flux(jg)*&
           (alp(jg)+u1lele(1)+sss0)
         abit=abit+sig(jg)*flux(jg)
      enddo
      csig(ig,ix)=csig(ig,ix)/abit
   enddo

   mfh=3
   mth=251
   if (nout.eq.0) go to 320
   b(1)=za
   b(2)=0
   b(3)=0
   b(4)=0
   b(5)=ngn
   b(6)=0
   ibase=6
   ip=0
   do 310 ig=1,ngn
      ip=ip+1
      b(ibase+ip)=csig(ig,ix)
      if (ip.lt.npage.and.ig.lt.ngn) go to 310
      if (ibase.eq.0) go to 300
      call listio(0,nout,0,b,nb,nwds)
      ibase=0
      ip=0
      go to 310
  300 call moreio(0,nout,0,b,nb,ip)
      ip=0
  310 continue
      call asend(nout,0)
  320 continue

   !--print cross sections in columns.
   mt=251
   write(nsyso,'(/,'' table of multigroup data'',//,&
     &'' group   lower       group     cosine'',/,&
     &''  no.    energy      flux  '',4x,4(a2,i3,7x))') hmt,mt
   write(nsyso,'('' -----   ------      ----- '',4x,4(2a5,2x))') (uline,i=1,2)
   do ig=1,ngn
      c(1)=csig(ig,1)
      write(nsyso,'(i5,1p,6e12.4)') ig,egn(ig),cflx(ig),c(1)
   enddo
   if (nout.ne.0) call afend(nout,0)
   return

   end subroutine musigc

   subroutine rdlgnd(nscr4,matd,mtd,npl,b,alp)
   !--------------------------------------------------------------------
   ! read legendre coefficient from nscr4 tape produced by
   ! subroutine grpav4.
   !--------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use util ! provdes repoz
   ! externals
   integer::nscr4,matd,mtd,npl
   real(kr)::b(*),alp(*)
   ! internals
   integer::nl,ngt,i,il,jg,nb,nwds

   !--set up header record
   call repoz(nscr4)
   call findf(matd,4,mtd,nscr4)
   call contio(nscr4,0,0,b,nb,nwds)
   nl=l2h
   ngt=n2h
   do i=1,ngt
      alp(i)=0
   enddo
   if (npl.gt.nl) return
   il=npl+7

   !--retrieve desired legendre coefficient
   jg=0
   do while (jg.lt.ngt)
      call listio(nscr4,0,0,b,nb,nwds)
      jg=n2h
      alp(jg)=b(il)
   enddo
   return
   end subroutine rdlgnd


   subroutine fssigc(ncg,ncm,nun,csig,cflx,b,egt,flux,sig)
   !--------------------------------------------------------------------
   ! Calculate the coarse group fission spectrum chi.
   !--------------------------------------------------------------------
   use mainio ! provides nsyso
   use endf ! provides endf routines and variables
   use util ! provides sigfig
   ! externals
   integer::ncg,ncm,nun
   real(kr)::csig(ncg,ncm),cflx(ncg),b(*),egt(nun+1),flux(nun),sig(nun+1)
   ! internals
   integer::nb,nw,nwds,ngnp1,i,np,loc,ig,jg,ibase,mt,nun1,ngn1,ix,ip
   real(kr)::c(6)
   character(2)::hmt='mt'
   character(5)::uline='-----'

   !--put the coarse group structure on nout, ala groupr
   if (nout.ne.0) then
      math=matd
      mfh=1
      mth=451
      b(1)=za
      b(2)=awr
      b(3)=iverf
      b(4)=0
      b(5)=-12
      b(6)=0
      call contio(0,nout,0,b,nb,nw)
      b(1)=tempin
      b(2)=0
      b(3)=ngn
      nw=6
      ngnp1=ngn+1
      do i=1,ngnp1
         nw=nw+1
         b(nw)=egn(i)
      enddo
      np=nw-6
      b(5)=np
      loc=1
      call listio(0,nout,0,b(loc),nb,nw)
      do while (nb.gt.0)
         loc=loc+nw
         call moreio(0,nout,0,b(loc),nb,nw)
      enddo
      call asend(nout,0)
      call afend(nout,0)
   endif

   !--initialize
   nun1=nunion+1
   do i=1,nun1
      egt(i)=sigfig(egt(i),ndig,0)
   enddo
   ngn1=ngn+1
   do i=1,ngn1
      egn(i)=sigfig(egn(i),ndig,0)
   enddo

   !--calculate coarse group flux
   do ig=1,ngn
      cflx(ig)=0
      do jg=1,nunion
         if (egt(jg).lt.egn(ig).or.egt(jg).ge.egn(ig+1)) cycle
         cflx(ig)=cflx(ig)+flux(jg)
      enddo
   enddo

   !--mt=18, chi: fission spectrum (mf=5/mt=18)
   !--compute cross-group cross sections and write on output tape.
   math=matd
   mfh=5
   mth=18
   call rdchi(math,b,sig)
   ix=1
   do ig=1,ngn
      csig(ig,ix)=0
      do jg=1,nunion
         if (egt(jg).lt.egn(ig).or.egt(jg).ge.egn(ig+1)) cycle
         csig(ig,ix)=csig(ig,ix)+sig(jg)*flux(jg)
      enddo
      csig(ig,ix)=csig(ig,ix)/cflx(ig)
   enddo
   if (nout.eq.0) go to 320
   b(1)=0
   b(2)=efmean
   b(3)=0
   b(4)=0
   b(5)=ngn
   b(6)=0
   ibase=6
   ip=0
   do 310 ig=1,ngn
      ip=ip+1
      b(ibase+ip)=csig(ig,ix)
      if (ip.lt.npage.and.ig.lt.ngn) go to 310
      if (ibase.eq.0) go to 300
      call listio(0,nout,0,b,nb,nwds)
      ibase=0
      ip=0
      go to 310
  300 call moreio(0,nout,0,b,nb,ip)
      ip=0
  310 continue
      call asend(nout,0)
  320 continue

   !--print data in columns.
   mt=18
   write(nsyso,'(/,'' table of multigroup data'',//,&
     &'' group   lower       group     chi   '',/,&
     &''  no.    energy      flux  '',4x,4(a2,i3,7x))') hmt,mt
   write(nsyso,'('' -----   ------      ----- '',4x,4(2a5,2x))') (uline,i=1,2)
   do ig=1,ngn
      c(1)=csig(ig,1)
      write(nsyso,'(i5,1p,6e12.4)') ig,egn(ig),cflx(ig),c(1)
   enddo
   if (nout.ne.0) call afend(nout,0)

   return
   end subroutine fssigc

   subroutine rdchi(matd,b,chi)
   !--------------------------------------------------------------------
   ! Read the fission spectrum (chi).
   !--------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use util ! provides repoz
   ! externals
   integer::matd
   real(kr)::b(*),chi(*)
   ! internals
   integer::ngt,i,jg,nb,nwds

   !--set up header record
   call repoz(ngout)
   call findf(matd,5,18,ngout)
   call contio(ngout,0,0,b,nb,nwds)
   ngt=n2h
   do i=1,ngt
      chi(i)=0
   enddo

   !--retrieve desired chi
   jg=0
   do while (jg.lt.ngt)
      call listio(ngout,0,0,b,nb,nwds)
      jg=n2h
      chi(jg)=b(8)
   enddo

   return
   end subroutine rdchi

   subroutine ggrmat(e,sigp,a,npnls,valspi)
   !--------------------------------------------------------------------
   ! Calculates R-matrix (Reich-Moore) cross sections contributions
   ! at energy e for one resonance
   !--------------------------------------------------------------------
   use physics ! provides pi, amassn
   ! externals
   integer::npnls
   real(kr)::e,sigp(4),a(*),valspi
   ! internals
   integer::i,inow,nlsmax,l,inowb,nrs,ncyc,ll
   integer::numj,jjl,jj,in,j
   integer::idone,iskip,kkkkkk,kchanl,kpstv,kngtv
   real(kr)::awri,aw,spi,gjd,k,pifac,rho,rhoc
   real(kr)::gfa,gfb,gf,se,pe,phi,phid,p1,p2,fl
   real(kr)::ajmax,ajc,gj,aj,dum1
   real(kr)::er,gn,gg,per,a1,a2,a3,den,de2,gg4
   real(kr)::t1,t2,t3,t4,termf,u11r,u11i,termt,termn
   real(kr)::dd,rr,ss,amag,rri,ssi,uur,uui,xx,termg
   real(kr)::r(3,3),s(3,3),ri(3,3),si(3,3)
   real(kr),parameter::rc1=.123e0_kr
   real(kr),parameter::rc2=.08e0_kr
   real(kr),parameter::third=0.333333333e0_kr
   real(kr),parameter::quar=0.25e0_kr
   real(kr),parameter::haf=0.50e0_kr
   real(kr),parameter::uno=1
   real(kr),parameter::two=2
   real(kr),parameter::four=4
   real(kr),parameter::small=3.0e-4_kr
   real(kr),parameter::zero=0

   !--initialize
   do i=1,4
      sigp(i)=0
   enddo
   inow=1
   in=0

   !--retrieve nuclide information
   naps=nint(a(inow+5))
   awri=a(inow+12)
   ap=a(inow+7)
   aw=amassn*awri
   ra=rc1*aw**third+rc2
   if (naps.eq.1) ra=ap
   spi=a(inow+6)
   gjd=2*(2*spi+1)
   nls=nint(a(inow+10))

   !--calculate wave number(k),rho and rhocap at energy (e)
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
   nlsmax=nls
   if (npnls.lt.nlsmax) nlsmax=npnls
   do l=1,nlsmax
      inowb=inow
      nrs=nint(a(inow+5))
      ncyc=nint(a(inow+4))/nrs
      ll=nint(a(inow+2))
      apl=a(inow+1)
      rhoc=k*ap
      rho=k*ra
      if (apl.ne.zero) rhoc=k*apl
      if (apl.ne.zero.and.naps.eq.1) rho=k*apl

      !--calculate shift and penetration factors at cross section energy
      call efacts(ll,rho,se,dum1)
      pe=dum1
      call efacphi(ll,rhoc,phi)
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
         if ((abs(ajc-valspi).gt.0.01e0_kr.or.l.ne.npnls).and.&
              npnls.ne.99) then
            in=(inow+6)+nrs*6+nrs*3
            go to 180
         endif
         gj=(2*ajc+1)/gjd

         !--loop over possible channel spins
         kchanl=0
         idone=0
         do while (kchanl.lt.2.and.idone.eq.0)
            kchanl=kchanl+1
            inow=inowb
            kpstv=0
            kngtv=0
            !--initialize matricies
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
               aj=abs(a(inow+1))

         !--select only resonances with current j value
               if (abs(aj-ajc).le.quar) then
                  if (a(inow+1).lt.zero) kngtv=kngtv+1
                  if (a(inow+1).gt.zero) kpstv=kpstv+1
                  iskip=0
                  if (kchanl.eq.1.and.a(inow+1).lt.zero) iskip=1
                  if (kchanl.eq.2.and.a(inow+1).gt.zero) iskip=1
                  if (iskip.eq.0) then
                     !--retrieve parameters
                     er=a(inow)
                     gn=a(inow+2)
                     gg=a(inow+3)
                     gfa=a(inow+4)
                     gfb=a(inow+5)
                     per=a(in+1)
                     a1=sqrt(gn*pe/per)
                     a2=0
                     if (gfa.ne.zero) a2=sqrt(abs(gfa))
                     if (gfa.lt.zero) a2=-a2
                     a3=0
                     if (gfb.ne.zero) a3=sqrt(abs(gfb))
                     if (gfb.lt.zero) a3=-a3
                     !--compute energy factors
                     diff=er-e
                     den=diff*diff+quar*gg*gg
                     de2=haf*diff/den
                     gg4=quar*gg/den
                     !--calculate r-function, or
                     !  calculate upper triangular matrix terms
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
               in=in+3
            enddo

            !--take care of channel spin per the sign of aj:
            !    kkkkkk = 0 => do not add anything in here;
            !    kkkkkk = 1 => add resonance contribution but not
            !                  extra hard-sphere;
            !    kkkkkk = 2 => add resonance plus hard-sphere phase
            !                  shift contribution.
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

                  !--invert the complex matrix
                  call efrobns(r,s,ri,si)

                  !--fission term for r-matrix path
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
                  if (abs(dd).lt.small.and.abs(phid).lt.small) then
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

               !--cross section contributions
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
  180    continue
      enddo
      inow=in
   !--continue loop over l values
   enddo

   !--calculate final cross sections and store for return
   sigp(1)=pifac*sigp(1)
   sigp(2)=pifac*sigp(2)
   sigp(3)=pifac*sigp(3)
   sigp(4)=pifac*sigp(4)

   return
   end subroutine ggrmat

   subroutine ggmlbw(e,sigp,a)
   !--------------------------------------------------------------------
   ! Calculates Multievel Breit-Wigner cross sections at energy e
   ! for one section (one isotope, one energy range).
   ! Based on "csmlbw" routine in reconr.
   !--------------------------------------------------------------------
   use physics ! provide pi, amassn
   ! externals
   real(kr)::e,sigp(4),a(*)
   ! internals
   integer::i,inow,nrs,lrx,in,j,ii,l,ll
   real(kr)::awri,aw,den,k,pifac,rho,rhoc
   real(kr)::qx,pec,rhop,sec,phi,cos2p,sin2p,sum,fl,ajmax
   real(kr)::aj,er,gn,gg,gf,ser,per,rper,gc,erp,edelt
   real(kr)::gne,gx,gtt,x,comfac,add,se,pe
   real(kr)::sigj(10,2),gj(10)
   real(kr),parameter::rc1=.123e0_kr
   real(kr),parameter::rc2=.08e0_kr
   real(kr),parameter::third=0.333333333e0_kr
   real(kr),parameter::half=0.5e0_kr
   real(kr),parameter::zero=0

   !--initialize
   do i=1,4
      sigp(i)=0
   enddo
   inow=1

   !--retrieve nuclide information
   naps=nint(a(inow+5))
   awri=a(inow+12)
   ap=a(inow+7)
   aw=amassn*awri
   ra=rc1*aw**third+rc2
   if (naps.eq.1) ra=ap
   spi=a(inow+6)
   den=4*spi+2
   nls=nint(a(inow+10))

   !--calculate wave number(k),rho and rhocap at energy (e)
   arat=awri/(awri+1)
   k=cwaven*arat*sqrt(abs(e))
   pifac=pi/(k*k)
   rho=k*ra
   rhoc=k*ap
   inow=inow+12

   !--loop over l states
   do l=1,nls
      nrs=nint(a(inow+5))
      ll=nint(a(inow+2))
      qx=a(inow+1)
      lrx=nint(a(inow+3))
      call efacts(ll,rho,se,pe)
      pec=0
      if (lrx.ne.0) then
         rhop=cwaven*arat*sqrt(abs(e+qx))*ra
         call efacts(ll,rhop,sec,pec)
      endif
      call efacphi(ll,rhoc,phi)
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
         er=a(inow)
         j=nint(a(inow+1)-ajmin+1)
         gn=a(inow+3)
         gg=a(inow+4)
         gf=a(inow+5)
         ser=a(in)
         per=a(in+1)
         rper=1/per
         gc=a(in+2)
         in=in+3
         inow=inow+6
         erp=er+gn*(ser-se)*rper/2
         edelt=e-erp
         gne=gn*pe*rper
         gx=gg+gf
         gtt=gne+gx
         gtt=gtt+gc*pec
         x=2*edelt/gtt
         comfac=2*gne/gtt/(1+x*x)
         sigj(j,1)=sigj(j,1)+comfac
         sigj(j,2)=sigj(j,2)+comfac*x
         comfac=comfac*gj(j)/gtt
         sigp(3)=sigp(3)+comfac*gf
         sigp(4)=sigp(4)+comfac*gg
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

   return
   end subroutine ggmlbw

   subroutine ssmlbw(e,sigp,a,aa)
   !--------------------------------------------------------------------
   ! Calculate Multilevel Breit-Wigner cross section contributions
   ! at energy e for one resonance
   !--------------------------------------------------------------------
   use physics ! provides pi
   ! externals
   real(kr)::e,sigp(4),a(6),aa(3)
   !internals
   integer::i,ii,j
   real(kr)::k,pifac,rho,rhoc,se,pe,pec,phi,cos2p,sin2p
   real(kr)::er,gn,gg,gf,ser,per,rper,gc,erp,edelt
   real(kr)::gne,gx,gtt,x,comfac
   real(kr)::sigj(10,2)

   !-initialize
   do i=1,4
      sigp(i)=0
   enddo

   !--compute cross sections for this energy
   k=cwaven*arat*sqrt(abs(e))
   pifac=pi/(k*k)
   rho=k*ra
   rhoc=k*ap

   !--calculate shift and penetration factors at cross section energy
   call efacts(ll,rho,se,pe)
   pec=0
   call efacphi(ll,rhoc,phi)
   cos2p=1-cos(2*phi)
   sin2p=sin(2*phi)
   do ii=1,2
      do i=1,10
         sigj(i,ii)=0
      enddo
   enddo

   !--retrieve resonance pararameters
   er=a(1)
   j=nint(a(2)-ajmin+1)
   gn=a(4)
   gg=a(5)
   gf=a(6)
   ser=aa(1)
   per=aa(2)
   rper=1/per
   gc=aa(3)
   erp=er+gn*(ser-se)*rper/2
   edelt=e-erp

   !--calculate the neutron widths at e
   gne=gn*pe*rper
   gx=gg+gf
   gtt=gne+gx
   gtt=gtt+gc*pec
   x=2*edelt/gtt

   !--common calculational factor
   comfac=2*gne/gtt/(1+x*x)

   !--elastic components
   sigj(j,1)=sigj(j,1)+comfac
   sigj(j,2)=sigj(j,2)+comfac*x
   comfac=comfac*gj(j)/gtt

   !--fission
   sigp(3)=sigp(3)+comfac*gf

   !--capture
   sigp(4)=sigp(4)+comfac*gg

   !--complete cross sections
   do j=1,nj
      sigp(2)=sigp(2)+gj(j)*((cos2p-sigj(j,1))**2+(sin2p+sigj(j,2))**2)
   enddo
   sigp(2)=sigp(2)+2*diff*cos2p
   sigp(2)=sigp(2)*pifac
   sigp(3)=sigp(3)*2*pifac
   sigp(4)=sigp(4)*2*pifac
   sigp(1)=sigp(2)+sigp(3)+sigp(4)

   return
   end subroutine ssmlbw

   subroutine ssslbw(e,sigp,a,aa)
   !--------------------------------------------------------------------
   ! Calculates Single Level Breit-Wigner cross sections at energy e
   ! for one resonance.
   !--------------------------------------------------------------------
   use physics ! provide pi
   ! externals
   real(kr)::e,sigp(4),a(*),aa(*)
   ! internals
   integer::i
   real(kr)::k,pifac,rho,rhoc,se,pe,pec,phi,cos2p,sin2p,sinsq
   real(kr)::spot,er,aj,gn,gg,gf,ser,per,rper,gc,gx,gj,erp
   real(kr)::edelt,gne,gtt,comfac

   !--initialize
   do i=1,4
      sigp(i)=0
   enddo

   !--retrieve resonance parameters
   k=cwaven*arat*sqrt(abs(e))
   pifac=pi/(k*k)
   rho=k*ra
   rhoc=k*ap
   call efacts(ll,rho,se,pe)
   pec=0
   call efacphi(ll,rhoc,phi)
   cos2p=cos(2*phi)
   sin2p=sin(2*phi)
   sinsq=(sin(phi))**2
   spot=4*(2*ll+1)*pifac*sinsq
   er=a(1)
   aj=a(2)
   gn=a(4)
   gg=a(5)
   gf=a(6)
   ser=aa(1)
   per=aa(2)
   rper=1/per
   gc=aa(3)
   gx=gg+gf
   gj=(2*aj+1)*spifac/2
   erp=er+gn*(ser-se)*rper/2
   edelt=e-erp
   gne=gn*pe*rper
   gtt=gne+gx
   gtt=gtt+gc*pec
   comfac=pifac*gj*gne/(edelt**2+gtt*gtt/4)

   !--compute cross sections for temp=0.
   !--elastic
   sigp(2)=sigp(2)+comfac*(gne*cos2p-2*gx*sinsq+2*edelt*sin2p)
   sigp(2)=sigp(2)+spot
   !--fission
   sigp(3)=sigp(3)+comfac*gf
   !--capture
   sigp(4)=sigp(4)+comfac*gg
   !--total
   sigp(1)=sigp(2)+sigp(3)+sigp(4)

   return
   end subroutine ssslbw

   subroutine ggunr1(e,sigp,a,amur,mxlru2)
   !--------------------------------------------------------------------
   ! Unresolved resonance region (format 1).
   ! Single Level Breit-Wigner formalism.
   ! Energy dependent parameters.
   ! Parameter interpolation is always used.
   !--------------------------------------------------------------------
   use physics ! provides pi
   ! externals
   integer::mxlru2
   real(kr)::e,a(*),sigp(4),amur(3,mxlru2)
   !internals
   integer::i,j,nlru2,inow,l,mu,nu,lamda,ll,njs
   real(kr)::awri,aw,const
   real(kr)::dx,aj,gnox,ggx,gfx,gxx,amun,gj,e2,k,rho,rhoc
   real(kr)::vl,ps,spot,gnx,den,temp,terg,ters,terf,gs,gc,gff,add
   real(kr),parameter::rc1=.123e0_kr
   real(kr),parameter::rc2=.08e0_kr
   real(kr),parameter::third=0.333333333e0_kr
   real(kr),parameter::zero=0

   !--initialize
   do i=2,4
      sigp(i)=0
   enddo
   spi=a(1)
   ap=a(2)
   nls=nint(a(5))
   nlru2=0
   inow=7

   !--do loop over all l states
   do l=1,nls
      awri=a(inow)
      ll=nint(a(inow+2))
      njs=nint(a(inow+5))
      arat=awri/(awri+1)
      aw=awri*amassn
      ra=rc1*aw**third+rc2
      const=(2*pi**2)/(cwaven*arat)**2
      inow=inow+6

      !--do loop over all j states
      do j=1,njs
         dx=a(inow)
         aj=a(inow+1)
         gnox=a(inow+2)
         ggx=a(inow+3)
         gfx=a(inow+4)
         gxx=a(inow+5)
         nlru2=nlru2+1
         amun=amur(1,nlru2)
         mu=nint(amur(1,nlru2))
         nu=nint(amur(2,nlru2))
         lamda=nint(amur(3,nlru2))
         gj=(2*aj+1)/(4*spi+2)
         e2=sqrt(e)
         k=arat*e2*cwaven
         rho=k*ra
         rhoc=k*ap

         !--calculate penetrability (vl) and phase shift(ps)
         call eunfac(ll,rho,rhoc,amun,vl,ps)
         vl=vl*e2

         !--calculate potential scattering
         if (j.eq.1) spot=4*pi*(2*ll+1)*(sin(ps)/k)**2

         !--compute cross section contributions
         gnx=gnox*vl
         diff=gxx
         den=e*dx
         temp=const*gj*gnx/den
         terg=temp*ggx
         ters=temp*gnx
         terf=temp*gfx

         !--calculate fluctuation integrals
         call egnrl(gnx,gfx,ggx,mu,nu,lamda,gs,diff,1)
         call egnrl(gnx,gfx,ggx,mu,nu,lamda,gc,diff,2)
         call egnrl(gnx,gfx,ggx,mu,nu,lamda,gff,diff,3)
         gc=gc*terg
         gff=gff*terf
         gs=gs*ters

         !--add interference correction
         add=const*gj*2*gnx*sin(ps)**2
         add=add/(e*dx)
         gs=gs-add

         !--cross sections
         sigp(2)=sigp(2)+gs
         sigp(3)=sigp(3)+gff
         sigp(4)=sigp(4)+gc
         inow=inow+6
      enddo
      sigp(2)=sigp(2)+spot
   enddo

   !--calculate total
   sigp(1)=sigp(2)+sigp(3)+sigp(4)

   return
   end subroutine ggunr1

   subroutine covadd(iadd,imt,imtmax,ntape,nout)
   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   use util ! provides openz
   ! externals
   integer::iadd,imtmax,ntape,nout
   integer::imt(imtmax)
   !internals
   integer::i,mat,mf,mt,ii,i1,i2
   character(66)::dat
   character(44)::dat2
   character(22)::czaawr

   call openz(ntape,0)
   call openz(nout,1)

   do i=1,4
      if (i.ne.2) then
         read(ntape,'(a66,i4,i2,i3,i5)') dat,mat,mf,mt,ii
         write(nout,'(a66,i4,i2,i3,i5)') dat,mat,mf,mt,ii
      else
         read(ntape,'(a22,a44,i4,i2,i3,i5)')czaawr,dat2,mat,mf,mt,ii
         write(nout,'(a22,a44,i4,i2,i3,i5)')czaawr,dat2,mat,mf,mt,ii
      endif
   enddo

   read(ntape,'(a44,i11,i11,i4,i2,i3,i5)') dat2,i1,i2,mat,mf,mt,ii
   write(nout,'(a44,i11,i11,i4,i2,i3,i5)') dat2,i1,i2+iadd,mat,mf,mt,ii

   do i=1,i1+i2
      read(ntape,'(a66,i4,i2,i3,i5)') dat,mat,mf,mt,ii
      write(nout,'(a66,i4,i2,i3,i5)') dat,mat,mf,mt,ii
   enddo

   do i=1,iadd
      write(nout,'(22x,4i11,i4,i2,i3,i5)') 33,imt(i),4,0,mat,mf,mt,ii+i
   enddo

 1000 continue
   read(ntape,'(a66,i4,i2,i3,i5)') dat,mat,mf,mt,ii
   write(nout,'(a66,i4,i2,i3,i5)') dat,mat,mf,mt,ii

   if (mf.eq.32.and.mt.eq.0) then
      write(nout,'(66x,i4,i2,i3,i5)') mat,0,0,99999
      do i=1,iadd
         write(nout,'(a22,4i11,i4,i2,i3,i5)')&
           czaawr,0,0,0,1,mat,33,imt(i),1
         write(nout,'(2f11.1,4i11,i4,i2,i3,i5)')&
           0.,0.,0,imt(i),0,1,mat,33,imt(i),2
         write(nout,'(2f11.1,4i11,i4,i2,i3,i5)')&
           0.,0.,1,5,3,2,mat,33,imt(i),3
         write(nout,'(a33,33x,i4,i2,i3,i5)') &
           ' 1.000000-5 2.000000+7 0.000000+0',mat,33,imt(i),4
         write(nout,'(66x,i4,i2,i3,i5)') mat,33,0,99999
      enddo
   endif

   if (mat.eq.-1) return
   go to 1000

   end subroutine covadd

   subroutine lumpxs(mti,mtk,sig,sig1,b,scr2)
   !--------------------------------------------------------------------
   ! Read the cross sections of the component MTs for a lumped
   ! covariance MT.
   !--------------------------------------------------------------------
   use endf ! provides endf routines and variables
   ! externals
   integer::mti,mtk
   real(kr)::sig(*),sig1(*),b(*),scr2(*)
   ! internals
   integer::i,l,mtl,nmtl,loc,mtd,nb,nw,j,izero

   izero=0
   mtl=0
   l=0
   do while (mtl.ne.mti.and.l.lt.nlump)
      l=l+1
      mtl=lump(1+2*(l-1))
   enddo
   nmtl=lump(1+2*(l-1)+1)
   do i=1,nunion
      if (mti.eq.mtk) then
         sig(i)=0
      else
         sig1(i)=0
      endif
   enddo

   !--loop over component mts
   loc=nlmt*(l-1)
   do i=1,nmtl
      mtd=lmt(i+loc)
      call findf(matd,mfcov,mtd,nscr3)
      call contio(nscr3,0,0,scr2,nb,nw)
      za=c1h
      awr=c2h
      call rdsig(matd,mtd,izero,b,scr2)
      do j=1,nunion
         if (mti.eq.mtk) then
            sig(j)=sig(j)+scr2(j)
         else
            sig1(j)=sig1(j)+scr2(j)
         endif
      enddo
   enddo

   return
   end subroutine lumpxs

   subroutine covout
   !--------------------------------------------------------------------
   ! Compute output covariances for all requested reactions (whether
   ! evaluated or derived) in the user-specified group structure.
   ! Output covariances are listed and/or written to a gendf tape.
   !
   ! Input union-group absolute covariances are read from unit nscr
   ! (see subroutine covcal).  Any non-zero input covariances for
   ! derived cross sections are ignored.  Coefficients relating
   ! derived and evaluated data reside in core at location a(ikxy).
   ! Fine-group energy bounds (iun), fluxes (iflx), and cross sections
   ! (isig) also reside in core. except for the trivial derivation case
   ! where both reactions ix and ixp are evaluated (isd=1), the entire
   ! nscr tape is read and a contribution to the output covariance is
   ! computed for each input reaction-pair.
   !
   !    ix,ixp = reaction indices (see array mts) of output reaction-
   !             pair (max. values = nmt,nmts)
   !    ig,igp = group indices of output reaction-pair (max. value =
   !             ngn)
   !    iy,iyp = reaction indices of current input reaction-pair (max.
   !             values = nmt,nmts)
   !    jg,jgp = group indices of current input reaction-pair (max.
   !             value = nunion)
   !
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso,nsyse
   use util ! provides error,openz,repoz,timer,closz
   use endf ! provides endf routines and variables
   ! internals
   integer::nmts,nwds,i,ntape,nb,nw
   integer::ix,iabort,ixp,izero,izap
   integer::isd,k,nmt1h,nmt1d,mt1lst,mt1old,nm
   integer::mat1,mt1,mta,jg,ig,ibase,ia,iy,igp,kp
   integer::jgp,ig2lo,ng2,ip,istart,nc,nmd,iyp
   integer::nngn,ig1,ig2,ijk,j,ld0,ld1,ld,loc,iscr,jscr
   integer::n,ngn2,mtl,lmtold,nmtold,itp,ldlst,ldold
   integer::irpc,iupc
   integer,dimension(:),allocatable::lmt1,lmt2
   real(kr)::egtjg,egtjgp,time,denom
   character(60)::strng
   real(kr),dimension(:),allocatable::xmu
   real(kr),dimension(:),allocatable::alp
   real(kr),dimension(:),allocatable::scr
   real(kr),dimension(:,:),allocatable::cova
   real(kr),dimension(:,:),allocatable::alsig,clflx
   real(kr),parameter::eps=1.e-20_kr
   real(kr),parameter::zero=0

   !--allocate storage.
   nmts=nmt1
   nwds=10000000
   nngn=ngn*(ngn+1)/2
   ngn2=ngn*ngn
   allocate(scr(nwds))
   if (mf32.ne.0) then
      allocate(cff(nngn))
      allocate(cfg(ngn2))
      allocate(cgg(nngn))
      allocate(cee(nngn))
      allocate(cef(ngn2))
      allocate(ceg(ngn2))
      allocate(ctt(nngn))
      allocate(ufg(ngn2))
      allocate(uef(ngn2))
      allocate(ueg(ngn2))
      allocate(uff(nngn))
      allocate(ugg(nngn))
      allocate(uee(nngn))
      allocate(utt(nngn))
   endif
   allocate(cflx(ngn))
   if (mfcov.eq.34) then
      allocate(xmu(ngn))
      allocate(alsig(ngn,2))
      allocate(clflx(ngn,2))
   endif
   if (nlump.gt.0) then
      allocate(lmt1(nmtmax))
      allocate(lmt2(nmtmax))
   endif
   allocate(alp(nunion))
   allocate(csig(ngn,nmt1))
   allocate(cova(ngn,ngn))
   do i=1,ngn
      do j=1,nmt1
         csig(i,j)=0
      enddo
   enddo
   if (nlump.gt.0) then
      do i=1,nmtmax
         lmt1(i)=0
         lmt2(i)=0
      enddo
      k=0
      do i=1,nlump
         mtl=lump(1+2*(i-1))
         n=lump(1+2*(i-1)+1)
         loc=nlmt*(i-1)
         do j=1,n
            k=k+1
            mtd=lmt(loc+j)
            lmt1(k)=mtd
            lmt2(k)=mtl
         enddo
      enddo
      lmtold=lmt1(1)
      nmtold=0
      lmt1(k+1)=1000
   else
      lmtold=1000
   endif

   !--position new gout tape (if any) for output.
   call repoz(nscr1)
   if (nout.ne.0) then
      call repoz(nout)
      call repoz(ngout)
      nsh=0
      call repoz(nin)
      ntape=nin
      if (nin.eq.0) ntape=ngout
      ! write a tape id on nout
      call tpidio(ntape,0,0,scr,nb,nw)
      math=0
      mfh=0
      mth=0
      call tpidio(0,nout,0,scr,nb,nw)
      if (nin.ne.0) then
         ! copy input covariance tape to nout
         call contio(nin,nout,0,scr,nb,nw)
         do
            call tomend(nin,nout,0,scr)
            call contio(nin,0,0,scr,nb,nw)
            if (math.eq.-1) exit
            call contio(0,nout,0,scr,nb,nw)
         enddo
      endif
   endif

   !--compute coarse-group cross sections.
   if (mfcov.eq.31) then
      call sigc(ngn,nmt1,nunion,csig,cflx,scr,un,flx,sig)
   else if (mfcov.eq.33) then
      call sigc(ngn,nmt1,nunion,csig,cflx,scr,un,flx,sig)
      if (irespr.eq.0) then
         call resprp(nwds,cflx,scr,ngn)
      else if (irespr.eq.1) then
         call resprx(nwds,scr)
      endif
   else if (mfcov.eq.34) then
      call musigc(ngn,nmt1,nunion,csig,cflx,scr,un,flx,sig,alp)
      do ijk=1,ngn
        xmu(ijk)=csig(ijk,1)
      enddo
   else if (mfcov.eq.35) then
      call fssigc(ngn,nmt1,nunion,csig,cflx,scr,un,flx,sig)
   else if (mfcov.eq.40) then
      call sigc(ngn,nmt1,nunion,csig,cflx,scr,un,flx,sig)
   endif
   call closz(nscr2)               !opened in covcal, finished with it now
   nwds=2*npage+50
   if (nwds.lt.ngn+6) nwds=ngn+6

   !--make a second copy of the fine-group covariance scratch tape
!   if (mf33.gt.0) then
      nscr2=12*imode
      call openz(nscr2,1)
      call repoz(nscr1)
      call repoz(nscr2)
      nsc=0
      call totend(nscr1,0,nscr2,scr)
      call repoz(nscr1)
      call repoz(nscr2)
!   endif

   !--loop over first reactions
   do 170 ix=1,nmt
   iabort=0
   mtd=mts(ix)
   izap=mzap(ix)

   !--write the head record for this section on nout
   math=matd
   mfh=mfcov
   if (mfh.eq.31) mfh=33
   if (nlump.gt.0.and.lmtold.lt.mts(ix)) then
      do while (lmtold.lt.mts(ix))
         mth=lmtold
         scr(1)=za
         if (mfcov.ne.40) then
            scr(2)=awr
         else
            scr(2)=izap
         endif
         scr(3)=0
         scr(4)=lmt2(nmtold+1)
         scr(5)=0
         scr(6)=0
         call contio(0,nout,0,scr,nb,nw)
         call asend(nout,0)
         nmtold=nmtold+1
         lmtold=lmt1(nmtold+1)
      enddo
   endif
   mth=mts(ix)
   scr(1)=za
   if (mfcov.ne.40) then
      scr(2)=awr
   else
      scr(2)=izap
   endif
   scr(3)=0
   scr(4)=0
   scr(5)=0
   scr(6)=nmts-ix+1
   if (mfcov.eq.34) then
      mth=251
      scr(4)=irelco
      scr(5)=legord
      scr(6)=legord
   endif
   call contio(0,nout,0,scr,nb,nw)

   !--loop over second reactions
   do 180 ixp=ix,nmts
   izero=0

   !--check for the trivial derivation case where both ix and ixp
   !--are directly evaluated. if it is, set isd=1.
   isd=0
   if (iabort.ne.1) then
      do k=1,nek
         if (akxy(ix,ix,k).eq.zero) go to 185
         if (akxy(ixp,ixp,k).eq.zero) go to 185
      enddo
      isd=1
  185 continue
   endif

   !--initialize this covariance matrix
   do i=1,ngn
      do j=1,ngn
         cova(i,j)=0
      enddo
   enddo

   nscr1=(11+isd)*imode
   if (isd.ne.1) call repoz(nscr1)

   !--accumulate contributions to this matrix
   !--from all matrices and fine groups on tape.
  200 continue
   if (isd.eq.1) then
     if (ix.ne.ixp) go to 210
   endif
   call contio(nscr1,0,0,scr,nb,nw)
   if (mfh.eq.0) go to 390
   if (mth.eq.0) go to 200
   if (mfcov.eq.34) then
      nmt1h=n1h
   else
      nmt1h=l2h
   endif
   if (isd.ne.1) go to 205
   nmt1d=nmt1h
   nmd=0
   if (mth.eq.mts(ix)) go to 205
   ! skip empty covariance matrices for the derived mts
   call tosend(nscr1,0,0,scr)
   go to 200
  205 continue
   mt1lst=1000
   mt1old=mt1lst
   nm=0
   if (mfcov.eq.34) then
      ldlst=-1
      ldold=-1
   endif
  210 continue
   if (isd.eq.1) then
       if (nmd.ge.nmt1d) go to 390
   endif

   !--loop over list records
  220 continue
   call listio(nscr1,0,0,scr,nb,nwds)
   if (mfcov.eq.35) then
      if (ifissp-1.ne.nm) then
         jg=n2h
         go to 380
      endif
   endif
   mat1=l1h
   mt1=l2h
   mta=mt1+1000*mat1
   if (mfcov.eq.34) then
      if (mth.eq.0) go to 180
      ld=nint(c1h)
      ld1=nint(c2h)
      ld0=ld*100+ld1
      if (ld0.ne.ldlst) k=1
      if (isd.ne.1) go to 225
      if (ld0.ne.ldlst) nmd=nmd+1
      if (mt1.eq.mts(ixp).and.mat1.eq.mats(ixp)) go to 225
      write(strng,'(''ld='',i3,''  ld1='',i3,''  mt1='',i3)') ld,ld1,mt1
      call error('covout','illegal condition for sad.',strng)
   endif
   if (mta.ne.mt1lst) k=1
   if (isd.ne.1) go to 225
   if (mta.ne.mt1lst) nmd=nmd+1
   if (mt1.eq.mts(ixp).and.mat1.eq.mats(ixp)) go to 225
   ! skip empty covariance matrices for the derived and
   ! non-requested mt1-s
   if (nmd.lt.nmt1d) then
      do while (nb.ne.0)
         call moreio(nscr1,0,0,scr,nb,nwds)
      enddo
      go to 220
   endif

   !--desired mt1 missing.  write empty matrix and abort
   !--speed-up logic for this mt.
   iabort=1
   do while (nb.ne.0)
      call moreio(nscr1,0,0,scr,nb,nwds)
   enddo
   go to 390
  225 continue
   if (mfcov.eq.34) then
      if (ld0.ne.ldlst) nm=nm+1
      ldlst=ld0
   else if (mfcov.ne.35) then
      if (mta.ne.mt1lst) nm=nm+1
   endif
   mt1lst=mta
   nw=n1h
   jg=n2h
   if (nw.eq.1.and.scr(7).eq.zero) go to 380
   ! index first coarse group
   egtjg=un(jg)
   do i=1,ngn
      ig=i
      if (egtjg.ge.egn(i).and.egtjg.lt.egn(i+1)) exit
   enddo
   ! read in rest of data for this group
   ibase=6
   ia=0
   do i=1,nw
      ia=ia+1
      cov(i)=scr(ibase+ia)
      if (nb.gt.0.and.ibase+ia.ge.nwds) then
         call moreio(nscr1,0,0,scr,nb,nwds)
         ibase=0
         ia=0
      endif
   enddo
   if (egtjg.ge.egn(ngn+1)) go to 380
   ! index reactions
   if (mfcov.eq.34) then
      if (ld0.eq.ldold) go to 280
   else
      if (mta.eq.mt1old) go to 280
   endif
   iy=0
   do i=1,nmt
      if (mth.eq.mts(i)) iy=i
   enddo
   iyp=0
   do i=1,nmts
      if (mt1.eq.mts(i).and.mat1.eq.mats(i)) iyp=i
   enddo
   if (iy.eq.0.or.iyp.eq.0) call error('covout',&
     'unable to find iy or iyp from mts array.',' ')
   mt1old=mta
   if (mfcov.eq.34) ldold=ld0
   ! index derived energy range for jg
  280 continue
   if (egtjg.ge.ek(k).and.egtjg.lt.ek(k+1)) go to 300
   if (k.eq.nek) go to 380
   k=k+1
   go to 280
  300 continue
   igp=1
   kp=1
   do 310 jgp=1,nw
   if (cov(jgp).eq.zero) go to 310
   egtjgp=un(jgp)
   ! index derived energy range for jgp
  320 continue
   if (egtjgp.ge.ek(kp).and.egtjgp.lt.ek(kp+1)) go to 330
   if (kp.eq.nek) go to 310
   kp=kp+1
   go to 320
   ! index second coarse group
  330 continue
   if (egtjgp.ge.egn(igp).and.egtjgp.lt.egn(igp+1)) go to 350
   if (igp.eq.ngn) go to 310
   igp=igp+1
   go to 330
   ! add this contribution
  350 continue
   if (iyp.eq.iy) then
      cova(igp,ig)=cova(igp,ig)+akxy(iy,ix,k)*akxy(iyp,ixp,kp)*cov(jgp)
      if (cova(igp,ig).ne.zero) izero=1
   else
      cova(igp,ig)=cova(igp,ig)+akxy(iy,ix,k)*akxy(iyp,ixp,kp)*cov(jgp)
      if (cova(igp,ig).ne.zero) izero=1
      cova(ig,igp)=cova(ig,igp)+akxy(iyp,ix,kp)*akxy(iy,ixp,k)*cov(jgp)
      if (cova(ig,igp).ne.zero) izero=1
   endif
  310 continue
   ! close loops over groups and covariance matrices.
  380 continue
   if (jg.lt.nunion) go to 220
   if (mfcov.eq.35) then
      nm=nm+1
      if (nm.eq.ifissp) go to 390
      go to 383
   endif
   if (isd.ne.1) go to 385
   if (izero.eq.0) go to 390
   if (ix.ne.iy.or.ixp.ne.iyp) call error('covout',&
     'unexpectedly, ix ne iy or ixp ne iyp..',' ')
   ! in the trivial derivation case, terminate loops over mt and mt1
  383 continue
   if (nmt1h.gt.1.and.nm.lt.nmt1h) go to 220
   go to  390
  385 continue
   if (mfcov.eq.34) call error('covout','please check isd=1',' ')
   if (nm.lt.nmt1h) go to 220
   go to 200
  390 continue

   ! add contribution from resonance-parameter uncertainty
   if (mfcov.ne.34.and.mfcov.ne.35) then
      call rescon(ix,ixp,csig,cova,izero,ngn,nmt1)
   endif

   !--write out covariance matrix elements for this
   !--range of coarse groups.
   math=matd
   mfh=mfcov
   if (mfh.eq.31) mfh=33
   mth=mts(ix)
   scr(1)=0
   scr(2)=0
   scr(3)=mats(ixp)
   scr(4)=mts(ixp)
   scr(5)=0
   scr(6)=ngn
   if (mfcov.eq.34) then
      mth=251
      scr(3)=mth
      scr(4)=ld
      scr(5)=ld1
   endif
   nwds=6
   call contio(0,nout,0,scr,nb,nwds)
   call timer(time)
   if (mats(ixp).ne.0) then
      if (irelco.eq.0) write(nsyso,'(/&
        &'' absolute covariance ( mt'',i3,'' , ig , mat'',i4,&
        &''/mt'',i3,'' , igp )'',11x,f9.1,''s''/)') mth,mats(ixp),mts(ixp),time
      if (irelco.eq.1) write(nsyso,'(/&
        &'' relative covariance ( mt'',i3,'' , ig , mat'',i4,&
        &''/mt'',i3,'' , igp )'',11x,f9.1,''s''/)') mth,mats(ixp),mts(ixp),time
   else
      if (mfcov.eq.31.or.mfcov.eq.33.or.mfcov.eq.35.or.mfcov.eq.40) then
         if (irelco.eq.0) write(nsyso,'(/&
           &'' absolute covariance ( mt'',i3,'' , ig , mt'',i3,&
           &'' , igp )'',19x,f9.1,''s''/)') mth,mts(ixp),time
         if (irelco.eq.1) write(nsyso,'(/&
           &'' relative covariance ( mt'',i3,'' , ig , mt'',i3,&
           &'' , igp )'',19x,f9.1,''s''/)') mth,mts(ixp),time
         if (mfcov.eq.40) write(nsyso,'(''final metastable state lfs'',i3/)')lfs
      else if (mfcov.eq.34) then
         if (irelco.eq.0) write(nsyso,'(&
           &/,'' absolute covariance ( mt'',i3,'' , ig , mt'',i3,&
           &'' , igp )'',19x,f9.1,''s'',&
           &/,''             same as ( mt'',i3,'' , ig , mt'',i3,&
           &'' , igp )'')') mts(ix),mts(ixp),time,251,251
         if (irelco.eq.1) write(nsyso,'(&
           &/,'' relative covariance ( mt'',i3,'' , ig , mt'',i3,&
           &'' , igp )'',19x,f9.1,''s'',&
           &/,''             same as ( mt'',i3,'' , ig , mt'',i3,&
           &'' , igp )'')') mts(ix),mts(ixp),time,251,251
      endif
   endif
   if (mfcov.eq.34) write(nsyso,'('' for legendre component: '',i2,&
     &'' and '',i2)') ld,ld1
   if (iprint.ne.0) then
      if (izero.eq.1) write(nsyso,'(&
        &''  ig   igp         +0         +1         +2''/&
        &'' ---   ---       ----       ----       ----'')')
      if (izero.eq.0) write(nsyso,'('' zero'')')
   endif
   if (mfcov.eq.34) then
      call alsigc(ngn,nmt1,nunion,alsig,clflx,scr,un,flx,sig,alp,ld,ld1,mtd,mt1)
      mth=251
   endif
   do ig=1,ngn
      ig2lo=0
      ng2=0
      ip=0
      ibase=6
      do igp=1,ngn
         ip=ip+1
         ! calculate absolute covariances
         if (mfcov.eq.34) then
            scr(ibase+ip)=cova(igp,ig)/(clflx(ig,1)*clflx(igp,2))
         else
            scr(ibase+ip)=cova(igp,ig)/(cflx(ig)*cflx(igp))
         endif
         if (scr(ibase+ip).ne.zero.and.irelco.ne.0) then
            ! calculate relative covariances
            if (mfcov.eq.34) then
               denom=xmu(igp)*xmu(ig)
            else
               denom=csig(ig,ix)*csig(igp,ixp)
            endif
            if (mfcov.eq.35) then
               if (denom.gt.zero) then
                  denom=max(denom,eps)
                  scr(ibase+ip)=scr(ibase+ip)/denom*&
                    (egn(ig+1)-egn(ig))*(egn(igp+1)-egn(igp))
               else
                  scr(ibase+ip)=0
               endif
            else
               if (denom.gt.zero) then
                  denom=max(denom,eps)
                  scr(ibase+ip)=scr(ibase+ip)/denom
               else
                  scr(ibase+ip)=0
               endif
            endif
         endif
         if (abs(scr(ibase+ip)).le.eps) then
            if (ig2lo.eq.0) ip=ip-1
         else
            if (ig2lo.eq.0) ig2lo=igp
            ng2=igp
         endif
      enddo
      if (ng2.ne.0.or.ig.ge.ngn) then
         if (ng2.eq.0) ig2lo=ig
         if (ng2.eq.0) ng2=ig
         scr(1)=0
         scr(2)=0
         scr(3)=ng2-ig2lo+1
         scr(4)=ig2lo
         scr(5)=ng2-ig2lo+1
         scr(6)=ig
         ip=ng2-ig2lo+1
         istart=1
         if (mfcov.eq.34) mfh=mfcov
         call listio(0,nout,0,scr(istart),nb,nw)
         do while (nb.ne.0)
            istart=istart+nw
            call moreio(0,nout,0,scr(istart),nb,nw)
         enddo
         if (izero.ne.0.and.iprint.ne.0) then
            ibase=6
            nw=ip
            do while (nw.gt.0)
               nc=nw
               if (nc.gt.6) nc=6
               write(nsyso,'(i4,i6,1p,6e11.3)')&
                 ig,ig2lo,(scr(ibase+i),i=1,nc)
               ibase=ibase+nc
               ig2lo=ig2lo+nc
               nw=nw-nc
            enddo
         endif
      endif
   enddo

   !--write out diagonal covariance matrix elements of resonance
   !--parameters.
   if (mfcov.eq.33.and.mf32.ne.0) then
      itp=0
      irpc=0
      iupc=0
      ! 1 is tt
      ! 2 is ee, 3 is eg, 4 is ef
      ! 5 is ff, 6 is fg
      ! 7 is gg
      if (mth.eq.18) then
         if (mts(ixp).eq.18) then
            itp=1
            irpc=5
            iupc=5
         else if (mts(ixp).eq.102) then
            itp=1
            irpc=6
            iupc=6
         endif
      else if (mth.eq.102.and.mts(ixp).eq.102) then
         itp=1
         irpc=7
         iupc=7
      else if (mth.eq.2) then
         if (mts(ixp).eq.2) then
            itp=1
            irpc=2
            iupc=2
         else if (mts(ixp).eq.18) then
            itp=1
            irpc=4
         else if (mts(ixp).eq.102) then
            itp=1
            irpc=3
         endif
      else if (mth.eq.1.and.mts(ixp).eq.1) then
         itp=1
         irpc=1
         iupc=1
     endif
     if (itp.eq.0) go to 590
     do ig1=1,ngn
         ig2=ngn*(ig1-1)-(ig1-1)*(ig1-2)/2+1
         if (irpc.eq.1.and.ctt(ig2).gt.zero) go to 510
         if (irpc.eq.2.and.cee(ig2).gt.zero) go to 510
         if (irpc.eq.3.and.ceg(ig2).gt.zero) go to 510
         if (irpc.eq.4.and.cef(ig2).gt.zero) go to 510
         if (irpc.eq.5.and.cff(ig2).gt.zero) go to 510
         if (irpc.eq.6.and.cfg(ig2).gt.zero) go to 510
         if (irpc.eq.7.and.cgg(ig2).gt.zero) go to 510
         if (iupc.eq.1.and.utt(ig2).gt.zero) go to 510
         if (irpc.eq.2.and.uee(ig2).gt.zero) go to 510
         if (iupc.eq.3.and.ueg(ig2).gt.zero) go to 510
         if (iupc.eq.4.and.uef(ig2).gt.zero) go to 510
         if (iupc.eq.5.and.uff(ig2).gt.zero) go to 510
         if (iupc.eq.6.and.ufg(ig2).gt.zero) go to 510
         if (iupc.eq.7.and.ugg(ig2).gt.zero) go to 510
      enddo
      go to 590
  510 continue
      if (ifresr.gt.0.and.(ifunrs.gt.0.and.iupc.gt.0)) then
         write(nsyso,'(&
           &/,5x,''...contribution from resonance parameters (mf=32)...''&
           & ,/,5x,''  ig   igp   resolved  unresolve''&
           & ,/,5x,'' ---   ---   --------  ---------'')')
      else if (ifresr.gt.0) then
         write(nsyso,'(&
           &/,5x,''...contribution from resonance parameters (mf=32)...''&
           & ,/,5x,''  ig   igp   resolved''&
           & ,/,5x,'' ---   ---   --------'')')
      else if (ifunrs.gt.0.and.iupc.gt.0) then
         write(nsyso,'(&
           &/,5x,''...contribution from resonance parameters (mf=32)...''&
           & ,/,5x,''  ig   igp  unresolve''&
           & ,/,5x,'' ---   ---  ---------'')')
      endif
      iscr=1
      jscr=iscr+ngn
      do i=1,ngn
         scr(iscr+i-1)=0
         scr(jscr+i-1)=0
      enddo
      do ig=ig1,ngn
         ig2=ngn*(ig-1)-(ig-1)*(ig-2)/2+1
         if (irpc.eq.1) scr(iscr+ig-1)=ctt(ig2)/cflx(ig)**2
         if (irpc.eq.2) scr(iscr+ig-1)=cee(ig2)/cflx(ig)**2
         if (irpc.eq.3) scr(iscr+ig-1)=ceg(ig2)/cflx(ig)**2
         if (irpc.eq.4) scr(iscr+ig-1)=cef(ig2)/cflx(ig)**2
         if (irpc.eq.5) scr(iscr+ig-1)=cff(ig2)/cflx(ig)**2
         if (irpc.eq.6) scr(iscr+ig-1)=cfg(ig2)/cflx(ig)**2
         if (irpc.eq.7) scr(iscr+ig-1)=cgg(ig2)/cflx(ig)**2
      enddo
      if (ifunrs.gt.0.and.iupc.gt.0) then
         do ig=ig1,ngn
            ig2=ngn*(ig-1)-(ig-1)*(ig-2)/2+1
            if (iupc.eq.1) scr(jscr+ig-1)=utt(ig2)/cflx(ig)**2
            if (iupc.eq.2) scr(jscr+ig-1)=uee(ig2)/cflx(ig)**2
            if (iupc.eq.3) scr(jscr+ig-1)=ueg(ig2)/cflx(ig)**2
            if (iupc.eq.4) scr(jscr+ig-1)=uef(ig2)/cflx(ig)**2
            if (iupc.eq.5) scr(jscr+ig-1)=uff(ig2)/cflx(ig)**2
            if (iupc.eq.6) scr(jscr+ig-1)=ufg(ig2)/cflx(ig)**2
            if (iupc.eq.7) scr(jscr+ig-1)=ugg(ig2)/cflx(ig)**2
         enddo
      endif
     if (irelco.eq.1) then
         do ig=ig1,ngn
            scr(iscr+ig-1)=scr(iscr+ig-1)/(csig(ig,ix)*csig(ig,ixp))
         enddo
         if (ifunrs.gt.0.and.iupc.gt.0) then
            do ig=ig1,ngn
               scr(jscr+ig-1)=scr(jscr+ig-1)/(csig(ig,ix)*csig(ig,ixp))
            enddo
         endif
      endif
      if (ifresr.gt.0.and.(ifunrs.gt.0.and.iupc.gt.0)) then
         do ig=ig1,ngn
            if (scr(iscr+ig-1).gt.zero.or.scr(jscr+ig-1).gt.zero) then
               write(nsyso,'(5x,i4,i6,1p,6e11.3)')&
                 ig,ig,scr(iscr+ig-1),scr(jscr+ig-1)
            endif
         enddo
      else if (ifresr.gt.0) then
         do ig=ig1,ngn
            if (scr(iscr+ig-1).gt.zero) then
               write(nsyso,'(5x,i4,i6,1p,6e11.3)') ig,ig,scr(iscr+ig-1)
            endif
         enddo
      else if (ifunrs.gt.0.and.iupc.gt.0) then
         do ig=ig1,ngn
            if (scr(jscr+ig-1).gt.0.) then
               write(nsyso,'(5x,i4,i6,1p,6e11.3)')&
                 ig,ig,scr(jscr+ig-1)
            endif
         enddo
      endif
  590 continue
   endif

   !--close loops over reaction types.
   nscr1=11*imode
   mt1lst=1000
  180 continue
   call asend(nout,0)
  170 continue

   !--covout is finished.
   call closz(nscr1)
   call closz(nscr2)
   if (mfcov.eq.34) call closz(nscr4)
   if (allocated(scr)) deallocate(scr)
   if (allocated(cff)) deallocate(cff)
   if (allocated(cfg)) deallocate(cfg)
   if (allocated(cgg)) deallocate(cgg)
   if (allocated(cee)) deallocate(cee)
   if (allocated(cef)) deallocate(cef)
   if (allocated(ceg)) deallocate(ceg)
   if (allocated(ctt)) deallocate(ctt)
   if (allocated(ufg)) deallocate(ufg)
   if (allocated(uef)) deallocate(uef)
   if (allocated(ueg)) deallocate(ueg)
   if (allocated(uff)) deallocate(uff)
   if (allocated(ugg)) deallocate(ugg)
   if (allocated(uee)) deallocate(uee)
   if (allocated(utt)) deallocate(utt)
   deallocate(cflx)
   deallocate(cova)
   deallocate(csig)
   deallocate(alp)
   if (allocated(lmt1)) deallocate(lmt1)
   if (allocated(lmt2)) deallocate(lmt2)
   if (allocated(xmu)) deallocate(xmu)
   if (allocated(alsig)) deallocate(alsig)
   if (allocated(clflx)) deallocate(clflx)
   if (allocated(crr)) deallocate(crr)
   if (nout.eq.0) return
   call afend(nout,0)
   call amend(nout,0)
   return
   end subroutine covout

   subroutine sigc(ncg,ncm,nun,csig,cflx,b,egt,flux,sig)
   !--------------------------------------------------------------------
   ! Calculate the coarse group cross sections.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso,nsyse
   use util ! provides repoz,closz,sigfig,error
   use endf ! provides endf routines and variables
   ! externals
   integer::ncg,ncm,nun
   real(kr)::csig(ncg,ncm),cflx(ncg),b(*),egt(nun+1),flux(nun),sig(nun+1)
   ! internals
   integer::nb,nw,ngnp1,i,np,loc,nun1,ngn1,ig,jg,ix,il,l,ndiff,mtd
   integer::mtl,nmtl,k,ibase,ip,nwds,nmtst,nmtend,ic,nline,ia,izero
   character(60)::strng
   real(kr),dimension(:),allocatable::sg
   character(2)::hmt='mt'
   real(kr)::c(6)
   integer::matp(80)

   izero=0

   !--put the coarse group structure on nout, ala groupr
   if (nlump.ne.0) then
      call repoz(nscr3)
      call closz(nscr3)
   endif
   if (nout.ne.0) then
      math=matd
      mfh=1
      mth=451
      b(1)=za
      b(2)=awr
      b(3)=iverf
      b(4)=0
      b(5)=-11
      if (mfcov.eq.40) b(5)=-14
      b(6)=0
      call contio(0,nout,0,b,nb,nw)
      b(1)=tempin
      b(2)=0
      b(3)=ngn
      nw=6
      ngnp1=ngn+1
      do i=1,ngnp1
         nw=nw+1
         b(nw)=egn(i)
      enddo
      np=nw-6
      b(5)=np
      loc=1
      call listio (0,nout,0,b(loc),nb,nw)
      do while (nb.ne.0)
         loc=loc+nw
         call moreio(0,nout,0,b(loc),nb,nw)
      enddo
      call asend(nout,0)
      call afend(nout,0)
   endif

   !--initialize
   allocate(sg(nunion))
   nun1=nunion+1
   do i=1,nun1
      egt(i)=sigfig(egt(i),ndig,0)
   enddo
   ngn1=ngn+1
   do i=1,ngn1
      egn(i)=sigfig(egn(i),ndig,0)
   enddo
   ! calculate coarse group flux
   do ig=1,ngn
      cflx(ig)=0
      do jg=1,nunion
         if (egt(jg).ge.egn(ig).and.egt(jg).lt.egn(ig+1)) then
            cflx(ig)=cflx(ig)+flux(jg)
         endif
      enddo
   enddo

   !--loop over all reaction types.
   !--compute cross-group cross sections and write on output tape.
   do ix=1,nmt1
      if (mts(ix).ge.851.and.mts(ix).le.870) then
         l=0
         do il=1,nlump
            mtl=lump(1+2*(il-1))
            if (mtl.eq.mts(ix)) l=il
         enddo
         if (l.eq.0) then
            write(strng,'(''covariance reaction'',i4)') mts(ix)
            call error('sigc',strng,'missing from lumping table.)')
         endif
         nmtl=lump(1+2*(l-1)+1)
         do jg=1,nunion
            sig(jg)=0
         enddo
         do k=1,nmtl
            mtd=lmt(nlmt*(l-1)+k)
            call rdsig(mats(ix),mtd,izero,b,sg)
            do jg=1,nunion
               sig(jg)=sig(jg)+sg(jg)
            enddo
         enddo
      else
         call rdsig(mats(ix),mts(ix),mzap(ix),b,sig)
      endif
      do ig=1,ngn
         csig(ig,ix)=0
         do jg=1,nunion
            if (egt(jg).ge.egn(ig).and.egt(jg).lt.egn(ig+1)) then
               csig(ig,ix)=csig(ig,ix)+sig(jg)*flux(jg)
            endif
         enddo
         if (cflx(ig).gt.0) then
            csig(ig,ix)=csig(ig,ix)/cflx(ig)
         else
            csig(ig,ix)=0
         endif
      enddo
      if (nout.ne.0.and.mats(ix).eq.0) then
         math=matd
         mfh=3
         mth=mts(ix)
         b(1)=za
         b(2)=mzap(ix)
         b(3)=0
         b(4)=0
         b(5)=ngn
         b(6)=0
         ibase=6
         ip=0
         do ig=1,ngn
            ip=ip+1
            b(ibase+ip)=csig(ig,ix)
            if (ip.ge.npage.or.ig.ge.ngn) then
               if (ibase.ne.0) then
                  call listio(0,nout,0,b,nb,nwds)
                  ibase=0
                  ip=0
               else
                  call moreio(0,nout,0,b,nb,ip)
                  ip=0
               endif
            endif
         enddo
         call asend(nout,0)
      endif
   enddo
   math=matd

   !--print cross sections in columns.
   nmtend=nmt1
   if (nmtend.gt.4) nmtend=4
   ic=0
   do i=1,nmt1
      matp(i)=matd
      if (mats(i).ne.0) then
         ic=ic+1
         matp(i)=mats(i)
      endif
   enddo
   if (ic.eq.0) write(nsyso,'(/&
     &'' table of multigroup cross sections''//&
     &'' group   lower       group     cross section''/&
     &''  no.    energy      flux  '',4x,4(a2,i3,7x))')&
     (hmt,mts(i),i=1,nmtend)
   if (ic.gt.0) write(nsyso,'(/&
     &'' table of multigroup cross sections''//&
     &'' group   lower       group     cross section''/&
     &''  no.    energy      flux  '',4x,4(i4,''/'',i3,4x))')&
     (matp(i),mts(i),i=1,nmtend)
   nline=2*nmtend
   write(nsyso,'(&
     &'' -----   ------      ----- '',4x,4(2a5,2x))')&
     ('-----',i=1,nline)
   do ig=1,ngn
      do ia=1,nmtend
         c(ia)=csig(ig,ia)
      enddo
      write(nsyso,'(i5,1p,6e12.4)')&
        ig,egn(ig),cflx(ig),(c(i),i=1,nmtend)
   enddo
   go to 510
  460 continue
   if (ic.eq.0) write(nsyso,'(/&
     &'' group cross section''/''  no.  '',6(a2,i3,7x))')&
     (hmt,mts(i),i=nmtst,nmtend)
   if (ic.gt.0) write(nsyso,'(/&
     &'' group cross section''/&
     &''  no.  '',6(i4,''/'',i3,4x))')&
     (matp(i),mts(i),i=nmtst,nmtend)
   nline=2*(nmtend-nmtst+1)
   write(nsyso,'('' ----- '',6(2a5,2x))') ('-----',i=1,nline)
   do ig=1,ngn
      izero=0
      do ia=1,ndiff
         c(ia)=csig(ig,ia-1+nmtst)
         if (c(ia).ne.0.) izero=1
      enddo
      if (izero.ne.0) then
         write(nsyso,'(i5,1p,6e12.4)') ig,(c(i),i=1,ndiff)
      endif
   enddo
  510 continue
   if (nmt1.eq.nmtend) go to 550
   nmtst=nmtend+1
   ndiff=nmt1-nmtst+1
   if (ndiff.gt.6) ndiff=6
   nmtend=nmtst+ndiff-1
   go to 460
  550 continue
   if (nout.ne.0) call afend(nout,0)
   if (nlump.ne.0) then
      deallocate(lmt)
      deallocate(lump)
   endif
   deallocate(sg)
   return
   end subroutine sigc

   subroutine resprp(nwscr,cflx,scr,ncg)
   !--------------------------------------------------------------------
   ! Prepare tables containing the resonance-parameter contributions
   ! to coarse-group covariances.
   !--------------------------------------------------------------------
   use util ! provides repoz,error,mess
   use endf ! provides endf routines and variables
   use physics ! provides amassn,amu,hbar,ev,pi
   ! externals
   integer::nwscr,ncg
   real(kr)::cflx(ncg),scr(*)
   ! internals
   integer::nb,nw,ig,is,ner,nro,ie
   integer::nl,l,nrs,iloc,nr,i,idis,j,id,jd,igmin,igmax
   integer::igu,iscr01,nj0,njs,njs6,nj1,jscr,lord
   real(kr)::cwaven,awri,xk,x2,er,aw,d,cc,e1,e2,em
   real(kr)::aj,gt,rgt,gn,gg,gf,gx,enext,wt,g,c,sf,se,sg,corr,test
   real(kr)::emid,f,gno,gnox,phi1,phi2,rgn,rho,rhoc
   integer,parameter::nparmx=60
   integer,parameter::igumax=20
   real(kr)::us(3,nparmx,igumax)
   real(kr)::rcov(nparmx,igumax)
   character(60)::strng1,strng2
   real(kr)::s(3,5),cov(5,5)
   real(kr),parameter::r1=0.123e0_kr
   real(kr),parameter::r2=0.333333333e0_kr
   real(kr),parameter::r3=0.08e0_kr
   real(kr),parameter::zero=0
   cwaven=sqrt(2*amassn*amu*ev)*1.e-12_kr/hbar

   nresg=0
   ifresr=0
   ifunrs=0
   if (mfcov.ne.33.and.mfcov.ne.34) return
   if (mf32.eq.0) return
   igmin=0
   igmax=0

   !--initialize
   call repoz(nendf)
   call tpidio(nendf,0,0,scr,nb,nw)
   do ig=1,ngn
      cff(ig)=0
      cfg(ig)=0
      cgg(ig)=0
      cef(ig)=0
      ceg(ig)=0
      cee(ig)=0
      ctt(ig)=0
      uff(ig)=0
      ufg(ig)=0
      ugg(ig)=0
      uee(ig)=0
   enddo
   call findf(matd,32,151,nendf)
   call contio(nendf,0,0,scr,nb,nw)
   nis=n1h

   !--loop over isotopes
   do 110 is=1,nis
   call contio(nendf,0,0,scr,nb,nw)
   abn=c2h
   ner=n1h
   if (iverf.eq.6) ner=n2h
   do 115 ie=1,ner
   call contio(nendf,0,0,scr,nb,nw)
   lru=l1h
   lrf=l2h
   if (lru.eq.1.and.lrf.ge.1.and.lrf.le.2) go to 116
   if (lru.eq.2.and.lrf.ge.1.and.lrf.le.2) go to 116
      write(strng2,'(''lrf='',i4,''  lru='',i4)') lrf,lru
      call error('resprp',&
        'illegal or unrecognized data structure in mf32',strng2)
  116 continue
   if (lru.eq.1.and.lrf.eq.1) call mess('resprp',&
     'for resolved resonance of single level breit-wigner,',&
     'contributions to total and elastic was not coded.')
   if (lru.ne.1) then
      write(strng2,'(''lru='',i4)') lru
      call error('resprp',&
        'illegal or unrecognized data structure in mf32',strng2)
      mf32=0
   endif
   if (lrf.ne.1.and.lrf.ne.2) then
      write(strng2,'(''lrf='',i4)') lrf
      call error('resprp',&
        'illegal or unrecognized data structure in mf32',strng2)
      mf32=0
   endif
   nro=n1h
   if (nro.ne.0) then
      write(strng2,'(''nro='',i4)') nro
      call error('resprp',&
        'illegal or unrecognized data structure in mf32',strng2)
      mf32=0
   endif
   naps=n2h
   call contio(nendf,0,0,scr,nb,nw)
   spi=c1h
   spifac=1/(2*spi+1)
   ap=c2h
   lcomp=l2h
   if (lcomp.ne.0) then
      write(strng2,'(''lcomp='',i4)') lcomp
      call error('resprp',&
        'illegal or unrecognized data structure in mf32',strng2)
      mf32=0
   endif
   nls=n1h
   if (lru.eq.2) go to 120

   !--resolved resonance parameters

   !--process all resonance parameters for this isotope
   do nl=1,nls
      ! read parameters for this l-value
      call listio(nendf,0,0,scr,nb,nw)
      l=1
      do while (nb.ne.0)
         l=l+nw
         call moreio(nendf,0,0,scr(l),nb,nw)
         if ((l+nw-1).gt.nwscr)&
           call error('resprp','storage exceeded.',' ')
      enddo
      awri=scr(1)
      l=nint(scr(3))
      nrs=nint(scr(6))
      xk=cwaven*awri/(awri+1)
      x2=2*(pi/xk)**2
      aw=amassn*awri
         ra=r1*aw**r2+r3
      if (naps.eq.1) ra=ap
      iloc=7

      !--loop over resonances
      do nr=1,nrs
         er=scr(iloc)
         aj=scr(iloc+1)
         gt=scr(iloc+2)
         rgt=1/gt
         gn=scr(iloc+3)
         gg=scr(iloc+4)
         gf=scr(iloc+5)

         !--index energy group
         if (er.ge.egn(1)) then
            if (er.ge.egn(ngn+1)) exit
            do i=1,ngn
               ig=i
               if (er.lt.egn(i+1)) exit
            enddo

            !--calculate fission and capture covariances in group ig
            !--due to this resonance
            if (ig.gt.nresg) nresg=ig
            if (iwt.ne.0) then
               ! retrieve user's weight function at er
               lord=0
               call egtwtf(er,enext,idis,lord,wt)
            else
               ! if weight function is unavailable, flat weight the resonance
               wt=cflx(ig)/(egn(ig+1)-egn(ig))
            endif
            g=(2*aj+1)*spifac/2
            c=abn*wt*x2*g/er

            !--estimate fission and capture cross section contributions
            !--from this resonance
            sf=c*gn*gf*rgt
            sg=c*gn*gg*rgt
            se=c*gn*gn*rgt

            !--calculate sensitivity of sf, sg, and se to resonance parameters
            s(1,1)=-sf/er
            s(2,1)=-sg/er
            s(3,1)=-se/er
            s(1,2)=sf*(1/gn-rgt)
            s(2,2)=sg*(1/gn-rgt)
            s(3,2)=se*(1/gn-rgt)
            s(1,3)=-sf*rgt
            s(2,3)=sg*(1/gg-rgt)
            s(3,3)=-se*rgt
            s(1,4)=0
            if (gf.gt.zero) s(1,4)=sf*(1/gf-rgt)
            s(2,4)=-sg*rgt
            s(3,4)=-se*rgt
            s(1,5)=2*sf/(2*aj+1)
            s(2,5)=2*sg/(2*aj+1)
            s(3,5)=2*se/(2*aj+1)

            !--retrieve resonance parameter covariances from the list record
            do i=1,5
               do j=1,5
                  cov(i,j)=0
               enddo
            enddo
            iloc=iloc+6
            cov(1,1)=scr(iloc)
            do i=2,5
               do j=2,i
                  iloc=iloc+1
                  cov(i,j)=scr(iloc)
                  cov(j,i)=scr(iloc)
               enddo
            enddo

            !--check covariance matrix for validity
            if (iverf.ge.6) then
               do i=1,5
                  id=i
                  jd=5
                  if (cov(i,5).ne.zero) then
                     write(strng2,'(''res parameters '',i1,'' and '',i1,&
                       &'' at er='',1p,e12.4)') id,jd,er
                     call error('resprp','bad covariance data for',strng2)
                  endif
               enddo
            endif
            do i=1,5
               id=i
               jd=i
               if (cov(i,i).lt.zero) then
                  write(strng2,'(''res parameters '',i1,'' and '',i1,&
                    &'' at er='',1p,e12.4)') id,jd,er
                  call error('resprp','bad covariance data for',strng2)
               endif
            enddo
            do i=1,5
               id=i
               do j=1,5
                  jd=j
                  if (cov(i,i).le.zero.or.cov(j,j).le.zero) then
                     if (cov(i,j).ne.zero) then
                        write(strng2,'(''res parameters '',i1,&
                          &'' and '',i1,'' at er='',1p,e12.4)')&
                          id,jd,er
                        call error('resprp','bad covariance data for',strng2)
                     endif
                  else
                     corr=cov(i,j)/sqrt(cov(i,i)*cov(j,j))
                     test=1
                     test=test+test/10000
                     if (abs(corr).ge.test) then
                        test=2
                        if (abs(corr).gt.test) then
                           write(strng2,'(''res parameters '',i1,&
                             &'' and '',i1,'' at er='',1p,e12.4)')&
                             id,jd,er
                           call error('resprp',&
                             'bad covariance data for',strng2)
                        endif
                        write(strng1,'(''correlation coeff='',f8.4)') corr
                        write(strng2,&
                          '(''for res parameters '',i1,'' and '',i1,&
                          &''at er='',1p,e12.4)')&
                          i,j,er
                        call mess('resprp',strng1,strng2)
                     endif
                  endif
               enddo
            enddo

            !--calculate cross section covariances by propagation of errors
            do i=1,5
               do j=1,5
                  cff(ig)=cff(ig)+s(1,i)*s(1,j)*cov(i,j)
                  cfg(ig)=cfg(ig)+s(1,i)*s(2,j)*cov(i,j)
                  cgg(ig)=cgg(ig)+s(2,i)*s(2,j)*cov(i,j)
                  cef(ig)=cef(ig)+s(1,i)*s(3,j)*cov(i,j)
                  ceg(ig)=ceg(ig)+s(2,i)*s(3,j)*cov(i,j)
                  cee(ig)=cee(ig)+s(3,i)*s(3,j)*cov(i,j)
               enddo
            enddo
            ctt(ig)=cee(ig)+cff(ig)+cgg(ig)
            iloc=iloc+16
         else
            iloc=iloc+2
         endif
      enddo
   enddo
   ifresr=1
   go to 115

   !--unresolved average Breit-Wigner resonance parameters
  120 continue
   emid=exp(log(elr*ehr)/2)
   do 410 i=1,ngn
      igmin=i
      if (elr.lt.egn(i+1)) go to 411
  410 continue
    write(strng2,'(''el='',1pe12.5)') elr
    call error('resprp','unresolved energy range was illegal.',strng2)
  411 continue
    do 412 i=igmin,ngn
       igmax=i
       if (ehr.lt.egn(i+1)) go to 413
  412 continue
    write(strng2,'(''eh='',1pe12.5)') ehr
    call error('resprp','unresolved energy range was illegal.',strng2)
  413 continue
   if (igmax.gt.nresg) nresg=igmax

   !--read resonance parameters for each l-value
   iscr01=1
   do nl=1,nls
      call listio(nendf,0,0,scr(iscr01),nb,nw)
      do
         iscr01=iscr01+nw
         if (nb.eq.0) exit
         call moreio(nendf,0,0,scr(iscr01),nb,nw)
         if ((iscr01+nw).gt.nwscr) then
            write(strng2,'(''require='',i8,''  supply='',i8,&
              &'' for nwds given in s.covout'')') iscr01+nw-1,nwscr
            call error('resprp','storage exceeded (lru=2).',strng2)
         endif
      enddo
   enddo

   !--read relative covariance from list record
   call listio(nendf,0,0,scr(iscr01),nb,nw)
   l=iscr01
   do
      if (nb.eq.0) exit
      l=l+nw
      call moreio(nendf,0,0,scr(l),nb,nw)
      if ((l+nw-iscr01).gt.nwscr) then
         write(strng2,'(''require='',i8,''  supply='',i8,&
           &'' for nwds given in s.covout'')') l+nw-iscr01,nwscr
         call error('resprp','storage exceeded (lru=2).',strng2)
      endif
   enddo

   !--retrieve relative covariance
   mpar=nint(scr(iscr01+2))
   nw=nint(scr(iscr01+4))
   npar=nint(scr(iscr01+5))
   if (mpar.lt.3.or.mpar.gt.5) then
      write(strng1,'(''mpar='',i2,'' was not coded.'')') mpar
      call error('resprp',strng1,' ')
   endif
   if (npar.gt.nparmx) then
      write(strng2,'(''npar='',i5,''  maximum='',i5)') npar,nparmx
      call error('resprp','storage exceeded for rel.covariance.',strng2)
   endif
   iloc=iscr01+5
   do i=1,npar
      do j=i,npar
         iloc=iloc+1
         rcov(i,j)=scr(iloc)
         rcov(j,i)=scr(iloc)
      enddo
   enddo

   !--check covariance matrix for validity
   do i=1,npar
      id=i
      jd=i
      if (rcov(i,i).lt.zero) then
         write(strng2,'(''res parameters '',i1,'' and '',i1)') id,jd
         call error('resprp','bad rel.covariance data for',strng2)
      endif
   enddo

   !--loop over l-state
   jscr=1
   nj0=0
   do nl=1,nls
      awri=scr(jscr)
      l=nint(scr(jscr+2))
      njs6=nint(scr(jscr+4))
      njs=nint(scr(jscr+5))
      xk=cwaven*awri/(awri+1)
      x2=2*(pi/xk)**2
      aw=amassn*awri
      ra=r1*aw**r2+r3
      if (naps.eq.1) ra=ap
      iloc=jscr+6

      !--loop over j-states
      do nj=1,njs
         d=scr(iloc)
         aj=scr(iloc+1)
         gnox=scr(iloc+2)
         gg=scr(iloc+3)
         gf=scr(iloc+4)
         gx=scr(iloc+5)
         iloc=iloc+6
         nj0=nj0+1
         nj1=mpar*(nj0-1)
         g=(2*aj+1)*spifac/2
         cc=abn*g*x2/d

         !--loop over groups
         do ig=igmin,igmax
            e1=egn(ig)
            e2=egn(ig+1)
            if (ig.eq.igmin.and.ig.eq.igmax) then
               f=(ehr-elr)/(egn(ig+1)-egn(ig))
               e1=elr
               e2=ehr
            else if (ig.eq.igmin) then
               f=(egn(ig+1)-elr)/(egn(ig+1)-egn(ig))
               e1=elr
            else if (ig.eq.igmax) then
               f=(ehr-egn(ig))/(egn(ig+1)-egn(ig))
               e2=ehr
            else
               f=1
            endif
            em=(e2+e1)/2
            if (iwt.eq.0) then
               wt=cflx(ig)/emid
            else
               lord=0
               call egtwtf(em,enext,idis,lord,wt)
            endif
            wt=wt*f
            rhoc=xk*sqrt(e1)*ap
            call efacphi(l,rhoc,phi1)
            rhoc=xk*sqrt(e2)*ap
            call efacphi(l,rhoc,phi2)

            !--correction of penetrability for reduced neutron width
            if (l.eq.0) then
               gno=gnox*sqrt(em)
            else if (l.eq.1) then
               rho=xk*sqrt(em)*ra
               gno=gnox*rho**2/(1+rho**2)*sqrt(em)
            else if (l.eq.2) then
               rho=xk*sqrt(em)*ra
               gno=gnox*rho**4/(9+3*rho**2+rho**4)*sqrt(em)
            endif
            gt=gno+gg+gf+gx
            rgn=1/gno
            rgt=1/gt
            c=cc*gno/gt/emid

            !--estimate cross section contributions
            sf=c*wt*gf
            sg=c*wt*gg
            se=c*wt*gno

            !--calculate sensitivity, s(i,j)
            !--i=1/2/3=fission/capture/elastic
            !--j=1/2/3/4/5=d/gn/gf/gg/gx
            igu=ig-igmin+1
            if (igu.gt.igumax) then
               write(strng2,'(''igu='',i5,''  maximum='',i5)') igu,igumax
               call error('resprp',&
                'storage exceeded for sensitivities.',strng2)
            endif
            us(1,nj1+1,igu)=-sf/d
            us(2,nj1+1,igu)=-sg/d
            us(3,nj1+1,igu)=-se/d
            if (mpar.eq.1) cycle
            us(1,nj1+2,igu)=sf*(rgn-rgt)
            us(2,nj1+2,igu)=sg*(rgn-rgt)
            us(3,nj1+2,igu)=se*(2*rgn-rgt)
            if (mpar.eq.2) cycle
            us(1,nj1+3,igu)=-sf*rgt
            us(2,nj1+3,igu)=sg*(1/gg-rgt)
            us(3,nj1+3,igu)=-se*rgt
            if (mpar.eq.3) cycle
            if (mpar.eq.4.and.lfw.eq.0) then
               us(1,nj1+4,igu)=-sf*rgt
               us(2,nj1+4,igu)=-sg*rgt
               us(3,nj1+4,igu)=-se*rgt
               cycle
            endif
            us(1,nj1+4,igu)=0
            if (gf.gt.zero) us(1,nj1+4,igu)=sf*(1/gf-rgt)
            us(2,nj1+4,igu)=-sg*rgt
            us(3,nj1+4,igu)=-se*rgt
            if (mpar.eq.4) cycle
            us(1,nj1+5,igu)=-sf*rgt
            us(2,nj1+5,igu)=-sg*rgt
            us(3,nj1+5,igu)=-se*rgt
         enddo
      enddo
      jscr=jscr+njs6+6
   enddo

   !--calculate relative covariance cross sections
   do ig=igmin,igmax
      igu=ig-igmin+1
      do i=1,npar
         do j=1,npar
            uff(ig)=uff(ig)+us(1,i,igu)*us(1,j,igu)*rcov(i,j)
            ufg(ig)=ufg(ig)+us(1,i,igu)*us(2,j,igu)*rcov(i,j)
            ugg(ig)=ugg(ig)+us(2,i,igu)*us(2,j,igu)*rcov(i,j)
            uee(ig)=uee(ig)+us(3,i,igu)*us(3,j,igu)*rcov(i,j)
         enddo
      enddo
   enddo
   ifunrs=1

   !--finished with this material
  115 continue
  110 continue
   return
   end subroutine resprp

   subroutine rescon(ix,ixp,csig,cova,izero,ncg,ncm)
   !--------------------------------------------------------------------
   ! Add the contributions from File 32 (see subroutines resprp and
   ! resprx) to the covariance previously calculated from File 33
   ! for this range of coarse groups.
   !--------------------------------------------------------------------
   ! externals
   integer::ix,ixp,izero,ncg,ncm
   real(kr)::csig(ncg,ncm),cova(ncg,ncg)
   ! internals
   integer::itp,iglast,ig,ig2,i,igind,ii,iif,iig,igd
   real(kr)::uegf
   real(kr),parameter::zero=0

   !--special branch for sammy method
   if (nmtres.gt.0) go to 1500

   if (mats(ixp).ne.0) return
   itp=0
   if (mts(ix).eq.18.and.mts(ixp).eq.18) itp=1
   if (mts(ix).eq.18.and.mts(ixp).eq.102) itp=2
   if (mts(ix).eq.102.and.mts(ixp).eq.102) itp=3
   if (mts(ix).eq.2.and.mts(ixp).eq.2) itp=4
   if (mts(ix).eq.2.and.mts(ixp).eq.18) itp=5
   if (mts(ix).eq.2.and.mts(ixp).eq.102) itp=6
   if (mts(ix).eq.1.and.mts(ixp).eq.1) itp=7
   if (itp.eq.0) return
   iglast=ngn
   if (iglast.gt.nresg) iglast=nresg
   if (ifresr.eq.0) go to 1000

   !--fission/fission
   if (itp.eq.1) then
      igind=0
      do ig=1,iglast
         do ig2=ig,ngn
            igind=igind+1
            cova(ig2,ig)=cova(ig2,ig)+cff(igind)
            if (ig.ne.ig2) then
               cova(ig,ig2)=cova(ig,ig2)+cff(igind)
            endif
         enddo
      enddo

   !--fission/capture
   else if (itp.eq.2) then
      do ig=1,ngn
         do ig2=1,ngn
            igd=(ig-1)*ngn+ig2
            cova(ig2,ig)=cova(ig2,ig)+cfg(igd)
         enddo
      enddo

   !--capture/capture
   else if (itp.eq.3) then
      igind=0
      do ig=1,ngn
         do ig2=ig,ngn
            igind=igind+1
            cova(ig2,ig)=cova(ig2,ig)+cgg(igind)
            if (ig.ne.ig2) then
               cova(ig,ig2)=cova(ig,ig2)+cgg(igind)
            endif
         enddo
      enddo

   !--elastic/elastic
   else if (itp.eq.4) then
      igind=0
      do ig=1,ngn
         do ig2=ig,ngn
            igind=igind+1
            cova(ig2,ig)=cova(ig2,ig)+cee(igind)
            if (ig.ne.ig2) then
               cova(ig,ig2)=cova(ig,ig2)+cee(igind)
            endif
         enddo
      enddo

   !--elastic/fission
   else if (itp.eq.5) then
      do ig=1,ngn
         do ig2=1,ngn
            igd=(ig-1)*ngn+ig2
            cova(ig2,ig)=cova(ig2,ig)+cef(igd)
         enddo
      enddo

   !--elastic/capture
   else if (itp.eq.6) then
      do ig=1,ngn
         do ig2=1,ngn
            igd=(ig-1)*ngn+ig2
            cova(ig2,ig)=cova(ig2,ig)+ceg(igd)
         enddo
      enddo

   !--total/total
   else if (itp.eq.7) then
      igind=0
      do ig=1,ngn
         do ig2=ig,ngn
            igind=igind+1
            cova(ig2,ig)=cova(ig2,ig)+ctt(igind)
            if (ig.ne.ig2) then
               cova(ig,ig2)=cova(ig,ig2)+ctt(igind)
            endif
         enddo
      enddo

   endif

   !--unresolved resonance contribution
 1000 continue
   if (ifunrs.eq.0) go to 2000
   if (irespr.eq.1) go to 1090
   !--convert to absolute covariance from relative
   if (ifunrs.eq.1) then
      ifunrs=2
      iif=0
      iig=0
      do ii=1,nmt
         i=ii
         if (mts(ii).eq.2) go to 1010
      enddo
      go to 1020
 1010 continue
      igind=0
      do ig=1,ngn
         do ig2=ig,ngn
            igind=igind+1
            if (csig(ig,i).le.zero) then
               uee(igind)=0
            else
               uee(igind)=uee(igind)*csig(ig,i)*csig(ig2,i)
            endif
         enddo
      enddo
 1020 continue
      do ii=1,nmt
         i=ii
         if (mts(ii).eq.18) go to 1030
      enddo
      go to 1040
 1030 continue
      iif=i
      igind=0
      do ig=1,ngn
         do ig2=ig,ngn
            igind=igind+1
            if (csig(ig,i).le.zero) then
               uff(igind)=0
            else
               uff(igind)=uff(igind)*csig(ig,i)*csig(ig2,i)
            endif
         enddo
      enddo
 1040 continue
      do ii=1,nmt
         i=ii
         if (mts(ii).eq.102) go to 1050
      enddo
      go to 1060
 1050 continue
      iig=i
      igind=0
      do ig=1,ngn
         do ig2=1,ngn
            igind=igind+1
            if (csig(ig,i).le.zero) then
               ugg(igind)=0
            else
               ugg(igind)=ugg(igind)*csig(ig,i)*csig(ig2,i)
            endif
         enddo
      enddo
 1060 continue
      if (iig.gt.0.and.iif.gt.0) then
         do ig=1,ngn
            if (csig(ig,iig).le.zero.or.csig(ig,iif).le.zero) then
               ufg(ig)=0
            else
               ufg(ig)=ufg(ig)*csig(ig,iig)*csig(ig,iif)
            endif
         enddo
      else
         do ig=1,ngn
            ufg(ig)=0
         enddo
      endif
   endif
 1090 continue

   !--fission/fission
   if (itp.eq.1) then
      igind=0
      do ig=1,iglast
         do ig2=ig,ngn
            igind=igind+1
            cova(ig2,ig)=cova(ig2,ig)+uff(igind)
            if (ig.ne.ig2) then
               cova(ig,ig2)=cova(ig,ig2)+uff(igind)
            endif
         enddo
      enddo

   !--fission/capture
   else if (itp.eq.2) then
      do ig=1,ngn
         do ig2=1,ngn
            igd=(ig-1)*ngn+ig2
            cova(ig2,ig)=cova(ig2,ig)+ufg(igd)
         enddo
      enddo

   !--capture/capture
   else if (itp.eq.3) then
      igind=0
      do ig=1,iglast
         do ig2=ig,ngn
            igind=igind+1
            cova(ig2,ig)=cova(ig2,ig)+ugg(igind)
            if (ig.ne.ig2) then
               cova(ig,ig2)=cova(ig,ig2)+ugg(igind)
            endif
         enddo
      enddo

   !--elastic/elastic
   else if (itp.eq.4) then
      igind=0
      do ig=1,iglast
         do ig2=ig,ngn
            igind=igind+1
            cova(ig2,ig)=cova(ig2,ig)+uee(igind)
            if (ig.ne.ig2) then
               cova(ig,ig2)=cova(ig,ig2)+uee(igind)
            endif
         enddo
      enddo

   !--elastic/fission
   else if (itp.eq.5) then
      do ig=1,ngn
         do ig2=1,ngn
            igd=(ig-1)*ngn+ig2
            cova(ig2,ig)=cova(ig2,ig)+uef(igd)
         enddo
      enddo

   !--elastic/capture
   else if (itp.eq.6) then
      do ig=1,ngn
         do ig2=1,ngn
            igd=(ig-1)*ngn+ig2
            cova(ig2,ig)=cova(ig2,ig)+ueg(igd)
         enddo
      enddo

   !--total/total
   else if (itp.eq.7) then
      if (irespr.eq.0) then
         igind=0
         do ig=1,iglast
            do ig2=ig,ngn
               igind=igind+1
               uegf=uee(igind)+ugg(igind)+uff(igind)
               cova(ig2,ig)=cova(ig2,ig)+uegf
               if (ig.ne.ig2) then
                  cova(ig,ig2)=cova(ig,ig2)+uegf
               endif
            enddo
         enddo
      else if (irespr.eq.1) then
         igind=0
         do ig=1,iglast
            do ig2=ig,ngn
               igind=igind+1
               cova(ig2,ig)=cova(ig2,ig)+utt(igind)
               if (ig.ne.ig2) then
                  cova(ig,ig2)=cova(ig,ig2)+utt(igind)
               endif
            enddo
         enddo
      endif

   endif
   go to 2000

   !--sammy method
 1500 continue
   do ig=1,ngn
      do ig2=1,ngn
          cova(ig,ig2)=cova(ig,ig2)+crr(ig,ig2,ix,ixp)
      enddo
   enddo

   !--check for nonzero array
 2000 continue
   do ig=1,ngn
      if (cova(ig,ig).ne.zero) izero=1
      if (izero.eq.1) exit
   enddo

   !--finished
   return
   end subroutine rescon

   subroutine grpav(mprint,tempin)
   !--------------------------------------------------------------------
   ! Compute multigroup cross sections for reactions needed in the
   ! calculation of the covariance matrices.  Calculation uses the
   ! union of the user specified group structure and the energy
   ! grid found in mfcov.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso,nsyse
   use util ! provides error,timer,repoz,mess
   use endf ! provides endf routines and variables
   ! externals
   integer::mprint
   real(kr)::tempin
   ! internals
   integer::nwds,i,nb,nw,np1,nl,indx,mfold,mtold,nshold
   integer::nsz,imt,ng2,idis,ii,ig,ig2lo,idone,infdel
   real(kr)::sec,etop,time,e,thresh,enext,elo,ehi
   real(kr)::sig(10)
   real(kr)::flux(10,10)
   character(60)::strng
   character(66)::text
   character(4)::mtname(17)
   real(kr)::b(8),z(20),ans(1,1,2)
   character(4)::tz(20)
   equivalence(tz(1),z(1))
   integer::nt=1
   integer::nz=1
   integer::ngg=0
   real(kr),dimension(:),allocatable::scr
   real(kr),parameter::eps=1.e-9_kr
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::elow=1.e-5_kr
   real(kr),parameter::zero=0

   !--initialize
   if (iread.eq.2) call error('grpav',&
     'not coded for multimaterial group averaging.',' ')
   call timer(sec)
   write(nsyso,'(/'' computing multigroup cross sections'',&
     &33x,f8.1,''s'')') sec
   call egnwtf
   nwds=nunion+10
   if ((npage+50).gt.nwds) nwds=npage+50
   allocate(scr(nwds))
   nsigz=1
   sigz(1)=big
   nscr2=11
   ngout=-10  !careful ... must match ntp (colaps) and ntape (errorr)!
   call repoz(ngout)
   math=0
   mfh=0
   mth=0
   text=' '
   read(text,'(16a4,a2)') (tz(i),i=1,17)
   call tpidio(0,ngout,0,z,nb,nw)
   if (abs(egn(1)-elow).le.eps) egn(1)=elow
   etop=un(1+nunion)

   !--search for desired mat and temperatures on pendf tape
   call repoz(npend)
   call findf(matd,1,0,npend)
   call contio(npend,0,0,scr,nb,nw)
   idone=0
   do while (idone.eq.0)
      za=c1h
      awr=c2h
      if (iverf.ge.5) call contio(npend,0,0,scr,nb,nw)
      if (iverf.ge.6) call contio(npend,0,0,scr,nb,nw)
      call hdatio(npend,0,0,scr,nb,nw)
      if (abs(c1h-tempin).le.tempin/10000) idone=1
      if (idone.eq.0) then
         if (c1h.gt.tempin) idone=2
         if (idone.eq.0) then
            call tomend(npend,0,0,scr)
            call contio(npend,0,0,scr,nb,nw)
            if (math.ne.matd) idone=2
         endif
      endif
   enddo
   if (idone.eq.2) then
      write(strng,'(''unable to find temp='',1p,e11.3)') tempin
      call error('grpav',strng,' ')
   endif

   !--write head record for this material on gout tape.
   nsh=1
   math=matd
   mfh=1
   mth=451
   scr(1)=za
   scr(2)=awr
   scr(3)=0
   scr(4)=nz
   scr(5)=-11
   scr(6)=nt
   call contio(0,ngout,0,scr,nb,nw)
   scr(1)=tempin
   scr(2)=0
   scr(3)=nunion
   scr(4)=ngg
   scr(5)=nunion+4
   scr(6)=0
   nw=7
   scr(nw)=0
   nw=nw+1
   scr(nw)=sigz(1)
   np1=nunion+1
   do i=1,np1
      scr(i+nw)=un(i)
   enddo
   nw=nw+np1
   ! add a zero for the ngg+1 entry, redefine nwds for listio
   scr(nw+1)=0
   nwds=nt+nz+nunion+1+ngg+1
   nl=1
   scr(5)=nw
   indx=1
   call listio(0,ngout,0,scr,nb,nwds)
   do while (nb.ne.0)
      indx=indx+nwds
      call moreio(0,ngout,0,scr(indx),nb,nwds)
   enddo
   mfold=1
   mtold=451
   nshold=nsh
   call asend(ngout,0)

   !--store total cross section from pendf tape for later use

   !--if this is an infinite dilution calculation,
   !--omit the reading and storing of mt=1
   infdel=1
   do nsz=1,nsigz
      if (sigz(nsz).lt.1.e8_kr) infdel=0
   enddo
   if (infdel.eq.1) then
      nscr2=0
   else
      call findf(matd,3,1,npend)
      if (npend.lt.0) nscr2=-nscr2
      call repoz(nscr2)
      nsh=1
      math=1
      call afend(nscr2,0)
      call contio(npend,nscr2,0,scr,nb,nw)
      call tosend(npend,nscr2,0,scr)
      call amend(nscr2,0)
      call atend(nscr2,0)
      call repoz(nscr2)
   endif
   nsh=nshold

   !--main loop over reactions
   mfd=3
   do imt=1,nga
      mtd=iga(imt)
      if (mfcov.eq.34 .and. mtd.ne.2) cycle
      if (mtd.eq.452.or.mtd.eq.455.or.mtd.eq.456) then
         write(strng,'(''cannot group average mt='',i3)') mtd
         call error('grpav',strng,&
           &'use groupr first, then error with ngout.ne.0')
      endif
      if (mtd.eq.3) then
         call mess('grpav','mt3 cross sections are constructed',&
        ' from total minus elastic')
      else
         text=' '
         read(text,'(15a4)') (mtname(i),i=1,15)
         call timer(time)

         !--initialize
         ng2=2
         nl=1
         e=0
         call egtsig(e,thresh,idis,sig)

         ! for threshold less than highest union group
         if (thresh.le.etop) then
            call egtflx(e,enext,idis,flux,nl,nz)
            if (mprint.ne.0) then
               if (tempin.eq.zero) write(nsyso,'(/&
                 &'' group constants at t=zero deg k'',37x,f8.1,&
                 &''s'')') time
               if (tempin.ne.zero) write(nsyso,'(/&
                 &'' group constants at t='',1p,e9.3,'' deg k'',&
                 &32x,0p,f8.1,''s'')') tempin,time
               write(nsyso,'('' for mf'',i2,'' and mt'',i3,1x,15a4)')&
                 mfd,mtd,(mtname(ii),ii=1,15)
               write(nsyso,'(15x,''group'',5x,''constant'')')
            endif

            !--loop over initial energy groups
            do ig=1,nunion
               elo=un(ig)
               ehi=un(1+ig)
               ig2lo=0
               if (ehi.gt.thresh) then

                  !--compute the group integrals
                  enext=ehi
                  ans(1,1,1)=0
                  ans(1,1,2)=0
                  idone=0
                  do while (idone.eq.0)
                     call epanel(elo,enext,ans,nl,nz,ig2lo,33)
                     if (enext.eq.ehi) idone=1
                     if (idone.eq.0) then
                        elo=enext
                        enext=ehi
                     endif
                  enddo

                  !--write this group on gout tape.
                  nw=nl*nz*ng2
                  if (ans(1,1,1).ne.zero) then
                     ans(1,1,2)=ans(1,1,2)/ans(1,1,1)
                  else
                     ans(1,1,2)=zero
                  endif
                  if (mprint.ne.0) write(nsyso,&
                    '(14x,i4,5x,1p,e11.4)') ig,ans(1,1,2)
                  mfh=mfd
                  mth=mtd
                  if (mfh.ne.mfold) then
                     call afend(ngout,0)
                     mfh=mfd
                     mth=mtd
                  endif
                  if (mth.ne.mtold) then
                     b(1)=za
                     b(2)=awr
                     b(3)=nl
                     b(4)=nz
                     b(5)=0
                     b(6)=nunion
                     nwds=6
                     call contio(0,ngout,0,b,nb,nwds)
                     mfold=mfd
                     mtold=mtd
                  endif
                  if (ans(1,1,2).ne.zero.or.ig.eq.nunion) then
                     b(1)=tempin
                     b(2)=0
                     b(3)=ng2
                     b(4)=ig2lo
                     b(5)=nw
                     b(6)=ig
                     b(7)=ans(1,1,1)
                     b(8)=ans(1,1,2)
                     nwds=8
                     call listio(0,ngout,0,b,nb,nwds)
                  endif
               endif
            enddo
            call asend(ngout,0)

         ! write message if mt has threshold gt highest union energy
         else
            write(strng,'(''mf '',i2,'' mt '',i3)') mfd,mtd
            call mess('grpav',strng,&
              'has threshold gt highest union energy.')
         endif
      endif
   enddo

   !--grpav is finished.
   deallocate(scr)
   call afend(ngout,0)
   call amend(ngout,0)
   call atend(ngout,0)
   call timer(sec)
   write(nsyso,'(/'' group averaging completed'',&
     &43x,f8.1,''s''/)') sec
   return
   end subroutine grpav

   subroutine colaps
   !--------------------------------------------------------------------
   ! Collapse (or expand) all MTs on unit ngout to the union grid,
   ! and write the new data onto a new gout tape (ngout=-10).
   ! Method assumes the cross section and the flux are constant
   ! in energy within an input group.
   !--------------------------------------------------------------------
   use util ! provides repoz,error,closz,sigfig
   use endf ! provides endf routines and variables
   ! internals
   integer::ntp,matn,nun1,nwscr,nb,nw,ntw,nsigz,ng,ng1
   integer::is,i,ib,inuf,jg,nl,iscrx,ig2lo,np
   integer::iscr2,nwl
   real(kr)::ea3,flxa,xnua,siga,ea1,ea2
   real(kr)::xnub,sigb,flxb,el,er,flux,temp
   real(kr)::eltst,ehtst,ea12,sfnuf,egu,chiu
   real(kr)::enl,enu,egl,dea,flxu,flxg,chig,enuu
   integer::nz=1
   integer::nt=1
   integer::ngg=0
   real(kr),dimension(:),allocatable::scr,scr0,scr18
   real(kr),dimension(:),allocatable::ela
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::zero=0

   if (iwt.ne.0) call egnwtf
   ntp=-10  !careful ... must match ngout (grpav) and ntape (errorr)!
   matn=0
   nun1=nunion+1
   call repoz(ngout)
   call repoz(ntp)

   !--test that ngout is really a groupr output tape
   !--then set scratch space based on problem dependent variables
   nwscr=17
   allocate(scr(nwscr))
   call tpidio(ngout,0,0,scr,nb,nw)
   call contio(ngout,0,0,scr,nb,nw)
   if (n1h.ne.-1) call error('colaps','ngout is not a groupr output tape',' ')
   call contio(ngout,0,0,scr,nb,nw)
   call repoz(ngout)
   deallocate(scr)
   if (mffis.ne.6) then
      nwscr=max(17,n1h+6,nunion+10)
   else
      nwscr=max(17,3*n1h+13,nunion+10)
   endif
   if (nwscr.lt.npage+50) nwscr=npage+50
   allocate(scr(nwscr))
   nsh=0
   call tpidio(ngout,ntp,0,scr,nb,nw)

   !--loop over tape
  110 continue
   call contio(ngout,0,0,scr,nb,nw)
   if (math.eq.-1) go to 400
   if (mfh.ne.1.or.mth.ne.451)&
     call error('colaps','did not find expected mf1, mt451.',' ')
   if (math.ne.matn) go to 120
   ! skip later temperatures for this material
   call tomend(ngout,0,0,scr)
   go to 110
  120 continue
   ntw=n2h
   scr(6)=1
   call contio(0,ntp,0,scr,nb,nw)
   if (matn.ne.0) deallocate(ela)
   matn=math
   nsigz=l2h
   call listio(ngout,0,0,scr,nb,nw)
   tempin=c1h
   if (nw.gt.nwscr) call error('colaps','storage exceeded.',' ')
   is=1
   ng=l1h
   do while (nb.ne.0)
      is=is+nw
      if (is.gt.nwscr) call error('colaps','storage exceeded.',' ')
      call moreio(ngout,0,0,scr(is),nb,nw)
   enddo
   ng1=ng+1
   allocate(ela(ng1))
   is=6+ntw+nsigz
   do i=1,ng1
      ela(i)=sigfig(scr(i+is),ndig,0)
   enddo
   if (ela(1).gt.un(1)) call error('colaps',&
     'ngout group structure does not span union grid.',' ')
   scr(3)=nunion
   scr(4)=ngg
   scr(5)=nt+nz+nun1+1
   scr(6)=0
   scr(7)=0
   scr(8)=big
   do i=1,nun1
      scr(8+i)=un(i)
   enddo
   scr(9+nun1)=0
   nw=3+nun1
   call listio(0,ntp,0,scr,nb,nw)
   is=1
   do while (nb.ne.0)
      is=is+nw
      if (is.gt.nwscr) call error('colaps','storage exceeded.',' ')
      call moreio(0,ntp,0,scr(is),nb,nw)
   enddo
   call afend(ntp,0)
   call tofend(ngout,0,0,scr)

   !--loop over all sections of this mat
  210 continue
   call contio(ngout,0,0,scr,nb,nw)
   if (math.eq.0) go to 390
   if (mffis.gt.0 .and. mfh.gt.mffis) then
      go to 380
   else if (mfh.eq.mffis) then
      if (mth.eq.18) then
         go to 300
      else
         call tosend(ngout,0,0,scr)
         go to 210
      endif
   else
      if (mth.eq.0) go to 210
      if (mfh.gt.3) then
         call tosend(ngout,0,0,scr)
         go to 210
      endif
   endif
   scr(6)=nunion
   call contio(0,ntp,0,scr,nb,nw)
   nl=nint(scr(3))
   nz=nint(scr(4))
   if (nl.ne.1.or.nz.ne.1) call error('colaps',&
     'not coded for multiple sigma zeroes or legendre orders.',' ')
   ib=0
   inuf=0
   if (mth.eq.452.or.mth.eq.455.or.mth.eq.456) inuf=1

   !--skip over low energy groups in input grid.
   ea3=0
   do while (ea3.le.un(1))
      call listio(ngout,0,0,scr,nb,nw)
      jg=n2h
      flxa=scr(7)
      xnua=scr(8)
      siga=scr(8+inuf)
      ea1=ea3
      ea2=ela(jg)
      ea3=ela(1+jg)
   enddo

   !--loop over output groups
  230 continue
   ib=ib+1
   xnub=0
   sigb=0
   flxb=0
   el=un(ib)
  250 continue
   if (ea2.gt.un(1+ib)) go to 280
   if (ea2.gt.el) el=ea2
   er=un(1+ib)
   if (ea3.lt.er) er=ea3
   flux=flxa*(er-el)/(ea3-ea2)
   xnub=xnub+siga*xnua*flux
   sigb=sigb+siga*flux
   flxb=flxb+flux
   if (ea3.gt.un(1+ib)) go to 280
   if (ea3.eq.un(1+ib).and.ib.ge.nunion) go to 280
   call listio(ngout,0,0,scr,nb,nw)
   jg=n2h
   flxa=scr(7)
   xnua=scr(8)
   siga=scr(8+inuf)
   ea1=ea3
   ea2=ela(jg)
   ea3=ela(jg+1)
   el=ea1
   if (ea1.eq.un(1+ib)) go to 280
   go to 250

   !--write results for this group
  280 continue
   if (sigb.eq.zero.and.ib.lt.nunion) go to 230
   nw=2+inuf
   temp=0
   scr(1)=temp
   scr(2)=0
   scr(3)=nw
   scr(4)=1
   scr(5)=nw
   scr(6)=ib
   scr(7)=flxb
   scr(8)=0
   scr(8+inuf)=0
   if (sigb.ne.zero.and.flxb.ne.zero) then
      scr(8)=xnub/sigb
      scr(8+inuf)=sigb/flxb
   endif
   call listio(0,ntp,0,scr,nb,nw)
   if (ib.lt.nunion) go to 230
   call tosend(ngout,0,0,scr)
   call asend(ntp,0)
   go to 210

   !--fission spectrum (chi)
  300 continue

   !--need additional scratch space for chi processing
   if (allocated(scr18)) deallocate(scr18)
   if (allocated(scr0)) deallocate(scr0)
   allocate(scr18(ng1))
   allocate(scr0(8))
   call repoz(ngout)
   call findf(matd,3,18,ngout)
   do i=1,ng1
      scr18(i)=0
   enddo
   call contio(ngout,0,0,scr,nb,nw)
  310 continue
   call listio(ngout,0,0,scr,nb,nw)
   jg=n2h
   scr18(jg)=scr(7)
   if (jg.lt.ng) go to 310
   call findf(matd,mffis,18,ngout)
   call contio(ngout,0,0,scr,nb,nw)
   scr(6)=nunion

   !--write an mf5 cont record to ntp regardless of mffis
   mfh=5
   call contio(0,ntp,0,scr,nb,nw)

   !--read the sectrum.  if mffis=5, it is a simple incident
   !--neutron energy dependent vector.  otherwise, it comes
   !--from the mf6 matrix.
   if (mffis.eq.5) then
      call listio(ngout,0,0,scr,nb,nw)
      is=1
      do while (nb.ne.0)
         is=is+nw
         if (is.gt.nwscr) call error('colaps','storage exceeded',' ')
         call moreio(ngout,0,0,scr(is),nb,nw)
      enddo
   else

      !--this coding assumes a groupr mf6 file created from an
      !--original endf mf5 input file.  The form of a groupr mf6
      !--file created from an endf mf6 input file differs (and
      !--will be addressed in a future update)!

      !--start by reading the "isotropic" spectrum, save at scr(iscrx)
      iscrx=1+ng+6
      call listio(ngout,0,0,scr,nb,nw)
      if (n2h.ne.0) call error('colaps','not ready for file6',' ')
      ig2lo=l2h
      np=n1h
      is=1
      do while (nb.ne.0)
         is=is+nw
         if (is.gt.nwscr) call error('colaps','storage exceeded',' ')
         call moreio(ngout,0,0,scr(is),nb,nw)
      enddo
      do i=1,ng
         scr(iscrx-1+i)=0
      enddo
      do i=1,np
         scr(iscrx+ig2lo-2+i)=scr(6+i)
      enddo
      do i=1,ng
         scr(6+i)=0
      enddo

      iscr2=iscrx+ng
      jg=0
      eltst=eclo
      ehtst=echi

      !--read the flux and either nu*sigf for groups in the
      !--"isotropic" energy range or nu*sigf*chi for higher
      !--energy groups and create a merged nu*sigf*chi vector
      !--for all groups in the eclo-to-echi energy interval,
      !--or for the group that encompasses efmean.
      do while(jg.lt.ng)
         call listio(ngout,0,0,scr(iscr2),nb,nw)
         ig2lo=l2h
         np=n1h
         jg=n2h
         is=iscr2
         do while (nb.ne.0)
            is=is+nw
            if (is.gt.nwscr) call error('colaps','storage exceeded.',' ')
            call moreio(ngout,0,0,scr(is),nb,nw)
         enddo
         flxa=scr(iscr2+6)
         ea1=ela(jg)
         ea2=ela(jg+1)
         if (igflag.ne.0) then
            eltst=ea1
            ehtst=ea2
         endif
         ! set flux weight only in the energy range
         if (efmean.ge.ehtst.or.efmean.le.eltst) then
            flxa=0
         else
            ea12=ea2-ea1
            if (ea1.lt.eltst) ea1=eltst
            if (ea2.gt.ehtst) ea2=ehtst
            flxa=flxa*(ea2-ea1)/ea12
         endif
         if(ig2lo.eq.0) then
            ! contribution comes from isotropic part of the spectrum
            sfnuf=scr(iscr2+7)
            do i=1,ng
               scr(6+i)=scr(6+i)+scr(iscrx-1+i)*flxa*sfnuf
            enddo
         else
            ! contribution comes from a matrix vector
            do i=2,np
               scr(4+ig2lo+i)=scr(4+ig2lo+i)+scr(iscr2+5+i)*flxa
            enddo
         endif
         ! if ea1>efmean, done.  Reset jg to force do while exit
         if (ea1.gt.efmean) jg=ng+1
      enddo

      !--normalize nu*sigf*chi
      sfnuf=0
      do i=1,ng
         sfnuf=sfnuf+scr(6+i)
      enddo
      do i=1,ng
         scr(6+i)=scr(6+i)/sfnuf
      enddo
      mfh=5

   endif

   !--translate GROUPR's chi and flux to nunion's group structure
   ib=0
   jg=0
   egu=0
   do 350 ib=1,nunion
      chiu=0
      flxu=0
      enl=un(ib)
      enu=un(ib+1)
      do while (egu.le.enl)
         jg=jg+1
         egu=ela(jg+1)
      enddo
      egl=ela(jg)
      dea=egu-egl
      flxg=scr18(jg)
      chig=scr(jg+6)
      if (er.le.ea3) then
         ! union group ib is within groupr group jg.
         flux=flxg*(enu-enl)/dea
         chiu=chig*flux
         flxu=flux
         else
         ! union group ib spans multiple groupr groups, make
         ! partial chiu, flxu calculation then start looping
         ! over additional groupr groups.
         enuu=egu
         flux=flxg*(enuu-enl)/dea
         chiu=chig*flux
         flxu=flux
         do while (egu.le.enu)
            jg=jg+1
            enl=enuu
            egl=egu
            egu=ela(jg+1)
            dea=egu-egl
            flxg=scr18(jg)
            chig=scr(jg+6)
            if (egu.le.enu) then
               enuu=egu
            else
               enuu=enu
            endif
            flux=flxg*(enuu-enl)/dea
            chiu=chiu+chig*flux
            flxu=flxu+flux
         enddo
      endif

      !--write the union group flux and chi (using the groupr
      !--mf3 format even though it is mf5 data) to ntp.
      nw=2
      nwl=8
      temp=0
      scr0(1)=temp
      scr0(2)=0
      scr0(3)=nw
      scr0(4)=1
      scr0(5)=nw
      scr0(6)=ib
      scr0(7)=flxu
      if (flxu.ne.0) then
         scr0(8)=chiu/flxu
      else
         scr0(8)=0
      endif
      call listio(0,ntp,0,scr0,nb,nw)
  350 continue
   call tosend(ngout,0,0,scr)
   call asend(ntp,0)
  380 continue
   call tomend(ngout,0,0,scr)
   call afend(ntp,0)
  390 continue
   call amend(ntp,0)
   go to 110

   !--finished
  400 continue
   call atend(ntp,0)
   call repoz(ngout)
!   call closz(ngout) !keep open for ngoutu

   !--now redefine ngout to be the colaps output tape
   ngout=ntp
   call repoz(ngout)
   deallocate(scr)
   deallocate(ela)
   if (allocated(scr18)) deallocate(scr18)
   if (allocated(scr0)) deallocate(scr0)
   return
   end subroutine colaps

   subroutine uniong(nendf)
   !--------------------------------------------------------------------
   ! Form union of user's energy mesh with mfcov energy mesh.
   !--------------------------------------------------------------------
   use util ! provides error,sigfig
   use endf ! provides endf routines and variables
   ! externals
   integer::nendf
   ! internals
   integer::nunmax,nb,nw,istart,iend,i,j,ngnp1,k,in
   character(60)::strng
   real(kr),dimension(:),allocatable::scr1
   real(kr),dimension(:),allocatable::scr2

   !--initialize.
   nw=npage+50
   allocate(scr1(nw))
   nunmax=5000
   allocate(scr2(nunmax))
   if (iverf.gt.4) go to 120

   !--read energy mesh from mfcov.
   call findf(matd,mfcov,0,nendf)
   call contio(nendf,0,0,scr1,nb,nw)
   call contio(nendf,0,0,scr1,nb,nw)
   call listio(nendf,0,0,scr1,nb,nw)
   istart=7
   neni=0
  100 continue
   iend=nw
   do i=istart,iend,2
      neni=neni+1
      if (neni.gt.nenimx) call error('uniong',&
        'exceeded storage in mfcov energy grid.',' ')
      eni(neni)=scr1(i)
   enddo
   if (nb.gt.0) then
      call moreio(nendf,0,0,scr1,nb,nw)
      istart=1
      go to 100
   endif
   do i=1,neni
      eni(i)=sigfig(eni(i),ndig,0)
   enddo

   !--unionize energy mesh.
  120 continue
   j=1
   ngnp1=ngn+1
   k=0
   do 130 i=1,ngnp1
  140 if (eni(j).lt.egn(1)) go to 160
   if (eni(j).le.egn(i)) go to 150
   k=k+1
   if (k.gt.nunmax) call error('uniong',&
     'exceeded storage in union energy grid.',' ')
   scr2(k)=egn(i)
   go to 130
  150 k=k+1
   if (k.gt.nunmax) call error('uniong',&
     'exceeded storage in union energy grid.',' ')
   scr2(k)=eni(j)
   j=j+1
   if (j.le.neni.and.egn(i).ne.eni(j-1)) go to 140
   if (j.le.neni.and.egn(i).eq.eni(j-1)) go to 130

   !--finished with endf energies
   in=i
   go to 170

   !--treat endf energies below first group boundary
  160 continue
   if (j.eq.neni) go to 130
   j=j+1
   go to 140
  130 continue

   !--finished if all ngn energies are used, as higher mfcov
   !--energies are not of interest.
   go to 200

   !--all mfcov energies used, some ngn energies left.
  170 continue
   do i=in,ngnp1
      if (egn(i).ne.scr2(k)) then
         k=k+1
         if (k.gt.nunmax) call error('uniong',&
           'exceeded storage in union energy grid.',' ')
         scr2(k)=egn(i)
      endif
   enddo
  200 continue
   nunion=k-1
   do i=2,k
      if (scr2(i).le.scr2(i-1)) then
         write(strng,'(1p,e12.4,'' le '',1p,e12.4)')&
           scr2(i),scr2(i-1)
         call error('uniong','union energies out of order',strng)
      endif
   enddo
   allocate(un(k))
   do i=1,k
      un(i)=scr2(i)
   enddo
   deallocate(eni)
   deallocate(scr2)
   deallocate(scr1)

   return
   end subroutine uniong

   subroutine epanel(elo,ehi,ans,nl,nz,iglo,mfcov)
   !--------------------------------------------------------------------
   ! Perform generalized group constant integrals for one panel.
   ! The upper boundary of the panel is chosen to be the smallest
   ! of ehi, the next cross section point, and the next flux point.
   !--------------------------------------------------------------------
   ! externals
   integer::nl,nz,iglo,mfcov
   real(kr)::elo,ehi,ans(nl,nz,*)
   ! internals
   integer::il,iz,n
   real(kr)::enext,en,ehigh,bq,rr
   real(kr)::sig(10),slst(10),flux(10,10),flst(10,10)
   integer::idisc=0
   integer::idiscf=0
   real(kr)::elast=0
   real(kr),parameter::delta=0.999995e0_kr
   save enext,flst,slst

   !--retrieve factors in integrands at lower boundary.
   il=1
   iz=1
   iglo=1
   if (idisc.ne.0.or.elo.ne.elast) then
      elast=elo
      if (mfcov.eq.34) then
         call egtlgc(elo,enext,idisc,slst)
      else
         call egtsig(elo,enext,idisc,slst)
      endif
      call egtflx(elo,en,idiscf,flst,nl,nz)
      if (en.eq.enext.and.idiscf.gt.idisc) idisc=idiscf
      if (en.lt.enext) idisc=idiscf
      if (en.lt.enext) enext=en
   endif

   !--retrieve cross section and flux at upper boundary.
   if (enext.lt.ehi) ehi=enext
   ehigh=ehi
   if (idisc.gt.0) ehigh=ehi*delta
   if (mfcov.eq.34) then
      call egtlgc(ehigh,enext,idisc,sig)
   else
      call egtsig(ehigh,enext,idisc,sig)
   endif
   call egtflx(ehigh,en,idiscf,flux,nl,nz)
   if (en.eq.enext.and.idiscf.gt.idisc) idisc=idiscf
   if (en.lt.enext) idisc=idiscf
   if (en.lt.enext) enext=en

   !--compute group cross sections and fluxes.
   bq=(ehigh-elo)/2
   ans(il,iz,1)=ans(il,iz,1)+(flux(iz,il)+flst(iz,il))*bq
   do n=1,nl
      rr=(sig(n)*flux(iz,il)+slst(n)*flst(iz,il))*bq
      ans(n,iz,2)=ans(n,iz,2)+rr
   enddo

   !--save last values.
   elast=ehi
   do n=1,nl
      slst(n)=sig(n)
   enddo
   flst(1,1)=flux(1,1)
   if (en.eq.enext.and.idiscf.gt.idisc) idisc=idiscf
   if (en.lt.enext) idisc=idiscf
   if (en.lt.enext) enext=en

   return
   end subroutine epanel

   subroutine egngpn
   !-------------------------------------------------------------------
   ! Generate requested neutron group structure or read in from
   ! the system input file in the form of an ENDF list record
   !
   !    ign     meaning
   !    ---     ---------------------------------------
   ! abs(1)     arbitrary structure (read in)
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
   !
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use util   ! provides error
   ! internals
   integer::lflag,ig,ngp,n1,n2,n,ic
   real(kr)::u,du,delta,ew,ewmin
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
   real(kr),parameter::uu80=.6931472e0_kr
   real(kr),parameter::e175=1.284e7_kr

   !--choose option according to ign
   lflag=0

   !--group structure is read in (free format)
   if (abs(ign).eq.1) then
      read(nsysi,*) ngn
      ngp=ngn+1
      if (allocated(egn)) deallocate(egn)
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
   else if (ign.eq.12.or.ign.eq.15) then
      ngn=620
      if (ign.eq.15) ngn=640
      ngp=ngn+1
      allocate(egn(ngp))
      egn(1)=sanda
      ! generate the first 45 boundaries
      do ig=1,8
         delta=deltl(ig)*sandb
         n1=ndelta(ig)
         n2=ndelta(ig+1)-1
         do n=n1,n2
            egn(n)=egn(n-1)+delta
         enddo
      enddo
      ! correct group 21
      egn(21)=sandc
      ! groups 46 to 450 are multiples of previous groups
      do ig=46,450
         egn(ig)=egn(ig-45)*10
      enddo
      ! groups 451 through 620 have constant spacing of 1.e5
      egn(451)=sandd
      do ig=452,ngp
         egn(ig)=egn(ig-1)+sande
      enddo

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
      egn(101)=-4*tenth
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

   !--shem cea 281-group structure
   else if (ign.eq.24) then
      ngn=281
      ngp=ngn+1
      allocate(egn(ngp))
      do ig=1,282
         egn(ig)=eg24(ig)
      enddo

   !--shem epm 295-group structure
   else if (ign.eq.25) then
      ngn=295
      ngp=ngn+1
      allocate(egn(ngp))
      do ig=1,296
         egn(ig)=eg25(ig)
      enddo

   !--shem cea/epm 361-group structure
   else if (ign.eq.26) then
      ngn=361
      ngp=ngn+1
      allocate(egn(ngp))
      do ig=1,362
         egn(ig)=eg26(ig)
      enddo

   !--shem epm 315-group structure
   else if (ign.eq.27) then
      ngn=315
      ngp=ngn+1
      allocate(egn(ngp))
      do ig=1,316
         egn(ig)=eg27(ig)
      enddo

   !--rahab aecl 89-group structure
   else if (ign.eq.28) then
      ngn=89
      ngp=ngn+1
      allocate(egn(ngp))
      do ig=1,90
         egn(ig)=eg28(ig)
      enddo

   !--ccfe  660-group structure
   else if (ign.eq.29) then
      ngn=660
      ngp=ngn+1
      allocate(egn(ngp))
      do ig=1,ngp
         egn(ig)=eg24(ig)
      enddo

   !--ukaea  1025-group structure
   else if (ign.eq.30) then
      ngn=1025
      ngp=ngn+1
      allocate(egn(ngp))
      do ig=1,ngp
         egn(ig)=eg25(ig)
      enddo

   !--ukaea  1067-group structure
   else if (ign.eq.31) then
      ngn=1067
      ngp=ngn+1
      allocate(egn(ngp))
      do ig=1,ngp
         egn(ig)=eg26(ig)
      enddo

   !--ukaea  1102-group structure
   else if (ign.eq.32) then
      ngn=1102
      ngp=ngn+1
      allocate(egn(ngp))
      do ig=1,ngp
         egn(ig)=eg27(ig)
      enddo

   !--ukaea  142-group structure
   else if (ign.eq.33) then
      ngn=142
      ngp=ngn+1
      allocate(egn(ngp))
      do ig=1,ngp
         egn(ig)=eg28(ig)
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

   !--display group structure - except for ign=-1
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
     &  '' neutron group structure......shem cea 281-group'')')
   if (ign.eq.25) write(nsyso,'(/&
     &  '' neutron group structure......shem epm 295-group'')')
   if (ign.eq.26) write(nsyso,'(/&
     &  '' neutron group structure......shem cea/epm 361-group'')')
   if (ign.eq.27) write(nsyso,'(/&
     &  '' neutron group structure......shem epm 315-group'')')
   if (ign.eq.28) write(nsyso,'(/&
     &  '' neutron group structure......rahab aecl 89-group'')')
   if (ign.eq.29) write(nsyso,'(/&
     &  '' neutron group structure......ccfe 660-group'')')
   if (ign.eq.30) write(nsyso,'(/&
     &  '' neutron group structure......ukaea 1025-group'')')
   if (ign.eq.31) write(nsyso,'(/&
     &  '' neutron group structure......ukaea 1067-group'')')
   if (ign.eq.32) write(nsyso,'(/&
     &  '' neutron group structure......ukaea 1102-group'')')
   if (ign.eq.33) write(nsyso,'(/&
     &  '' neutron group structure......ukaea 142-group'')')
   if (ign.ne.-1) then
      do ig=1,ngn
         write(nsyso,'(1x,i5,2x,1p,e12.5,''  - '',e12.5)')&
           ig,egn(ig),egn(ig+1)
      enddo
   endif
   ewmin=1.05e0_kr
   do ig=1,ngn
      if (egn(ig+1).lt.0.1e0_kr) then
         ew=egn(ig+1)/egn(ig)
         if (ew.lt.ewmin) ewmin=ew
      endif
   enddo
   if (ewmin.lt.1.05e0_kr) then
      eskip4=1.0e0_kr+0.5e0_kr*(ewmin-1.0e0_kr)
   endif

   !--prepare union of users grid with endf covariance grid.
   ngp=ngn+1
   do ig=1,ngp
      egn(ig)=sigfig(egn(ig),ndig,0)
   enddo
   call uniong(nendf)
   if (ign.eq.-1) then
      write(nsyso,'(/&
        &'' union structure (= user structure) has'',i5,&
        &'' groups''/)') nunion
      if (allocated(egn)) deallocate(egn)
      allocate(egn(nunion+1))
      do ig=1,nunion
         egn(ig)=un(ig)
      enddo
      egn(nunion+1)=un(1+nunion)
      ngn=nunion
   else
      write(nsyso,'(/'' union structure has'',i5,'' groups''/)')&
        nunion
   endif
   do ig=1,nunion
      write(nsyso,'(1x,i5,6x,1p,e11.5,''  -  '',e11.5)')&
        ig,un(ig),un(1+ig)
   enddo
   return
   end subroutine egngpn

   subroutine egnwtf
   !--------------------------------------------------------------------
   ! Set up calculation of weight functions or read in arbitary
   ! function in the form of an ENDF TAB1 record or
   ! read in parameters for an analytic weight function.
   !
   !    iwt     meaning
   !    ---     -------
   !     1      read in
   !     2      constant
   !     3      1/e
   !     4      1/e + fission spectrum + thermal maxwellian
   !     5      epri-cell lwr
   !     6      (thermal) -- (1/e) -- (fission + fusion)
   !     7      same with t-dep thermal part
   !     8      thermal--1/e--fast reactor--fission + fusion
   !     9      extended claw weight function
   !    10      claw with t-dependent thermal part
   !    11      vitamin-e weight function (ornl-5505)
   !    12      vit-e with t-dep thermal part
   !
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso,nsyse
   use util ! provides openz,error
   ! internals
   integer::iwtt,i,iww,iw
   real(kr)::pwr
   real(kr),dimension(:),allocatable::win
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
   real(kr),parameter::zero=0

   iwtt=iabs(iwt)
   if (allocated(wght)) deallocate(wght)

   !--arbitrary
   if (iwtt.eq.1) then
      write(nsyso,'(/'' weight function......read in'')')
      iw=10000
      allocate(win(iw))
      do i=1,iw
         win(i)=0
      enddo
      read(nsysi,*) (win(i),i=1,iw)
      iww=0
      do i=1,iw
         if (win(i).ne.zero) iww=i
      enddo
      allocate(wght(iww))
      do i=1,iww
         wght(i)=win(i)
      enddo
      deallocate(win)

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
         pwr=3
         pwr=pwr/2
         ac=1/(exp(-ec/tc)*ec**pwr)
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

   !--extended claw weight function
   else if (iwtt.eq.9.or.iwtt.eq.10) then
      write(nsyso,'(/&
        &'' weight function......extended claw weight function'')')
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
      write(nsyso,'(/,'' weight function......vitamin-e'')')
      if (iwtt.gt.11) write(nsyso,'(22x,''temperature dependent'')')

   !--illegal iwt
   else
      call error('egnwtf','illegal weight function requested.',' ')
   endif
   return
   end subroutine egnwtf

   subroutine egtwtf(e,enext,idis,lord,wtf)
   !--------------------------------------------------------------------
   ! Retrieve or compute required legendre component of the
   ! weight function constructed or read in by egnwtf.
   !--------------------------------------------------------------------
   use endf ! provides terpa
   use util ! provides sigfig
   use physics ! provides bk
   ! externals
   integer::idis,lord
   real(kr)::e,enext,wtf
   ! internals
   integer::ip,ir,iwtt,ipl
   real(kr)::test,ea,eb,enxt,tt,cc,step,bb,pow
   real(kr),parameter::con1=7.45824e7_kr
   real(kr),parameter::con2=1.e0_kr
   real(kr),parameter::con3=1.44934e-9_kr
   real(kr),parameter::con4=3.90797e-2_kr
   real(kr),parameter::con5=2.64052e-5_kr
   real(kr),parameter::con6=6.76517e-2_kr
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
   real(kr),parameter::two=2.e0_kr
   real(kr),parameter::veb=5.e5_kr
   real(kr),parameter::wt6a=.054e0_kr
   real(kr),parameter::wt6b=1.57855e-3_kr
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
   real(kr),parameter::zero=0
   save ir,ip,ipl,step

   !--initialize
   idis=0
   if (e.eq.zero) then
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
      if (iwtt.gt.6) tt=tempin*bk
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
      ea=bk*tempin
      eb=wt10a*tempin/wt10b
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
         if (iwtt.gt.11) tt=tempin*bk
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

   !--warn User when the weight function is zero; it most likely means
   !  the xs energy range extends beyond the weight function energy
   !  range.
   if (wtf.eq.zero.and.mwtf.eq.0) then
      mwtf=1
      call mess('egtwtf',&
                'xs energy range exceeds weight function range',&
                'some multgroup data may be suspect')
   endif

   !--return enext on an even grid
   enext=sigfig(enext,7,0)
   return
   end subroutine egtwtf

   subroutine egtflx(e,enext,idis,flux,nl,nz)
   !--------------------------------------------------------------------
   ! Retrieve or compute weighting fluxes.
   !--------------------------------------------------------------------
   use endf ! provides endf routines and variables
   ! externals
   integer::idis,nl,nz
   real(kr)::e,enext,flux(10,10)
   ! internals
   integer::nb,nw,iz,il,idisc,l
   real(kr)::en,tmin,wtf,t
   real(kr)::tot(10)
   real(kr),dimension(:),allocatable::atot
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::zero=0
   save atot

   !--initialize.
   if (e.eq.zero) then
      call egtwtf(e,en,idis,l,wtf)
      if (nscr2.eq.0) return
      call findf(matd,3,1,nscr2)
      nw=npage+50
      if (.not.allocated(atot)) allocate(atot(nw))
      call contio(nscr2,0,0,atot,nb,nw)
      call gety2(e,enext,idis,t,nscr2,atot)
      return
   endif

   !--compute self-shielded point flux assuming flux
   !--is proportional to the inverse total cross section.
   !--test for infinite dilution (i.e., nscr2=0)
   if (nscr2.eq.0) then
      enext=big
   else
      call gety2(e,enext,idis,t,nscr2,atot)
      do iz=1,nz
         tot(iz)=t
      enddo
   endif
   do il=1,nl
      l=il-1
      call egtwtf(e,en,idisc,l,wtf)
      if (en.lt.enext) idis=idisc
      if (en.lt.enext) enext=en
      do iz=1,nz
         flux(iz,il)=wtf
         if (nscr2.ne.0) then
            tmin=1
            tmin=tmin/1000
            if (tot(iz).le.zero) tot(iz)=tmin
            if (il.eq.1) then
               flux(iz,1)=flux(iz,1)/(1+tot(iz)/sigz(iz))
            else
               flux(iz,il)=flux(iz,il-1)/(1+tot(iz)/sigz(iz))
            endif
         endif
      enddo
   enddo
   return
   end subroutine egtflx

   subroutine egtsig(e,enext,idis,sig)
   !--------------------------------------------------------------------
   ! Retrieve the reaction cross-section defined by MFD and MTD.
   ! Remove discontinuities by moving second point up by eps.
   ! Initialize if e=0.
   !--------------------------------------------------------------------
   use util ! provides error
   use endf ! provides endf routines and variables
   ! externals
   integer::idis
   real(kr)::e,enext,sig(10)
   ! internals
   integer::nsig,mf,mt,nb,nw,iz
   real(kr)::s
   real(kr)::asig(350)
   save nsig,asig
   real(kr),parameter::zero=0

   !--initialize
   if (e.eq.zero) then
      nsig=npend
      mf=3
      if (mfd.eq.13.or.mfd.eq.17) mf=13
      mt=0
      if (mtd.le.200) mt=mtd
      if (mtd.eq.207) mt=mtd
      if (mtd.eq.261) mt=mtd
      if (mtd.ge.600.and.mtd.le.699) mt=mtd
      if (mtd.ge.700.and.mtd.le.799) mt=mtd
      if (mtd.ge.800.and.mtd.le.891) mt=mtd
      if (mtd.eq.251.or.mtd.eq.252.or.mtd.eq.253) mt=2
      if (mt.eq.0) call error('egtsig','mt=0.',' ')
      call findf(matd,mf,mt,nsig)
      nw=npage+50
      call contio(nsig,0,0,asig,nb,nw)
      call gety1(e,enext,idis,s,nsig,asig)

   !--retrieve point cross sections.
   else
      call gety1(e,enext,idis,s,nsig,asig)
      do iz=1,nsigz
         sig(iz)=s
      enddo
   endif
   return
   end subroutine egtsig

   subroutine efacphi(l,rho,phi)
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
   end subroutine efacphi

   subroutine efacts (l,rho,se,pe)
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
   end subroutine efacts

   subroutine eunfac(l,rho,rhoc,amun,vl,ps)
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
   end subroutine eunfac

   subroutine egnrl(galpha,gbeta,gamma,mu,nu,lamda,s,df,id)
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
   end subroutine egnrl

   subroutine efrobns(a,b,c,d)
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
   call ethrinv(a,3,ind)
   if (ind.ne.1) then
      call eabcmat(a,b,q)
      call eabcmat(b,q,d)
      do i=1,3
         do j=1,3
            c(i,j)=c(i,j)+d(i,j)
         enddo
      enddo
      call ethrinv(c,3,ind)
      call eabcmat(q,c,d)
      do i=1,3
         do j=1,3
            d(i,j)=-d(i,j)
         enddo
      enddo
   endif
   return
   end subroutine efrobns

   subroutine ethrinv(d,n,kimerr)
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
   end subroutine ethrinv

   subroutine eabcmat(a,b,c)
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
   end subroutine eabcmat

   subroutine matrixin(awr,res)
   !-------------------------------------------------------------------
   ! driver routine for matrix transformation
   !-------------------------------------------------------------------
   ! externals
   real(kr)::awr,res(10)
   ! internals
   integer::nlmax
   real(kr)::bc(800)

   nlmax=20
   call matrixej(awr,nlmax,bc,res)

   return
   end subroutine matrixin

   subroutine matrixej(awr,nm,bc,res)
   !-------------------------------------------------------------------
   ! calculates a transformation matrix.
   !-------------------------------------------------------------------
   ! externals
   integer::nm
   real(kr)::awr,bc(*),res(10)
   ! internals
   integer::m,mm,i,mup,i1,i2,ilow,l,l1,n1,n2,max,lll,ii,m2,m1,min,ilo
   real(kr)::a,g,z1,z2,z3,z4,x1,x2,sum,sum1
   real(kr):: t(65,65)
   real(kr),parameter::zero=0

   m=nm+1
   a=awr
   g=1/a
   mm=min0(2*m,30)
   do i=1,mm
      do l=1,m
         t(l,i)=0
      enddo
   enddo
   t(1,1)=1
   t(2,1)=2*g/3
   t(2,2)=1-6*g**2/10
   mup=0
   do 110 i=3,MM
      i1=i-1
      z1=i1
      z2=i1+2
      z3=2*i1-1
      z4=2*i1+3
      i2=i1-1
      x1=(z1/z3-z2*g**2/z4)*(-g)**I2
      t(2,I)=x1
      if (mup.ne.0) go to 105
      if (abs(x1).ge.1.0e-16_kr) go to 110
      mup=i1
  105 if (abs(x1).lt.1.0e-32_kr) go to 120
  110 continue
   if (mup.eq.0) mup=m
  120 continue
   ilo=1
   do l=3,m
      l1=l-2
      z1=2*l1+1
      z2=l1+1
      z3=l1
      ilow=ilo
      do 150 i=ilow,mm
         i1=i-1
         sum=-z3*t(l1,i)/z2
         do n1=1,mup
            x2=t(2,n1)
            if (abs(x2).eq.zero) cycle
            n2=n1-1
            max=n2+i1+1
            if (max.gt.mm) max=mm
            min=iabs(n2-i1)+1
            sum1=0
            do m1=min,max,2
               x1=t(l1+1,m1)
               if (abs(x1).lt.1.0e-16_kr) cycle
               m2=m1-1
               sum1=sum+cleb(n2,m2,i1)**2*x1
            enddo
            sum=sum+z1*x2*sum1/z2
         enddo
         if (i.ge.l) go to 147
         if (abs(sum).ge.abs(t(l-1,i))) go to 148
  147    t(l,i)=sum
         go to 150
  148    ilo=i+1
  150 continue
   enddo
   do i=1,m
      do l=1,m
         ii=i-1
         lll=l+m*ii
         bc(lll)=t(l,i)
         if (abs(bc(lll)).lt.1.0e-20_kr) bc(lll)=0
      enddo
   enddo
   do i=1,10
      res(i)=t(2,i)
   enddo

   return
   end subroutine matrixej

   real(kr) function cleb(i1,i2,i3)
   !-------------------------------------------------------------------
   !  cleb = Clebsch-Gordan coefficient
   !-------------------------------------------------------------------
   ! externals
   integer::i1,i2,i3
   ! internals
   integer::n1,n2,n3,n4,it,nept,is,ia,ib,ic,id
   real(kr)::arg,z1,d123,signx
   real(kr),parameter::fac(101)=(/&
     1.0e0_kr,1.0e0_kr,2.0e0_kr,6.0e0_kr,2.4e0_kr,1.2e0_kr,7.2e0_kr,&
     5.04e0_kr,4.032e0_kr,3.6288e0_kr,3.6288e0_kr,3.99168e0_kr,&
     4.790016e0_kr,6.2270208e0_kr,8.7178291e0_kr,1.3076744e0_kr,&
     2.0922790e0_kr,3.5568743e0_kr,6.4023737e0_kr,1.2164510e0_kr,&
     2.4329020e0_kr,5.1090942e0_kr,1.1240007e0_kr,2.5852017e0_kr,&
     6.2044840e0_kr,1.5511210e0_kr,4.0329146e0_kr,1.0888869e0_kr,&
     3.0488834e0_kr,8.8417620e0_kr,2.6525286e0_kr,8.2228387e0_kr,&
     2.6313084e0_kr,8.6833176e0_kr,2.9523280e0_kr,1.0333148e0_kr,&
     3.7199333e0_kr,1.3763753e0_kr,5.2302262e0_kr,2.0397882e0_kr,&
     8.1591528e0_kr,3.3452527e0_kr,1.4050061e0_kr,6.0415263e0_kr,&
     2.6582716e0_kr,1.1962222e0_kr,5.5026222e0_kr,2.5862324e0_kr,&
     1.2413916e0_kr,6.0828186e0_kr,3.0414093e0_kr,1.5511188e0_kr,&
     8.0658175e0_kr,4.2748833e0_kr,2.3084370e0_kr,1.2696403e0_kr,&
     7.1099859e0_kr,4.0526920e0_kr,2.3505613e0_kr,1.3868312e0_kr,&
     8.3209871e0_kr,5.0758021e0_kr,3.1469973e0_kr,1.9826083e0_kr,&
     1.2688693e0_kr,8.2476506e0_kr,5.4434494e0_kr,3.6471111e0_kr,&
     2.4800355e0_kr,1.7112245e0_kr,1.1978572e0_kr,8.5047859e0_kr,&
     6.1234458e0_kr,4.4701155e0_kr,3.3078854e0_kr,2.4809141e0_kr,&
     1.8854947e0_kr,1.4518309e0_kr,1.1324281e0_kr,8.9461821e0_kr,&
     7.1569457e0_kr,5.7971260e0_kr,4.7536433e0_kr,3.9455240e0_kr,&
     3.3142401e0_kr,2.8171041e0_kr,2.4227095e0_kr,2.1077573e0_kr,&
     1.8548264e0_kr,1.6507955e0_kr,1.4857160e0_kr,1.3520015e0_kr,&
     1.2438414e0_kr,1.1567725e0_kr,1.0873662e0_kr,1.0329978e0_kr,&
     9.9167793e0_kr,9.6192760e0_kr,9.4268904e0_kr,9.3326215e0_kr,&
     9.3326215e0_kr/)
   integer,parameter::nfac(101)=(/&
     0,0,0,0,1,2,2,3,4,5,6,7,8,9,10,12,13,14,&
     15,17,18,19,21,22,23,25,26,28,29,30,32,33,35,36,38,40,41,43,44,46,&
     47,49,51,52,54,56,57,59,61,62,64,66,67,69,71,73,74,76,78,80,81,83,&
     85,87,89,90,92,94,96,98,100,101,103,105,107,109,111,113,115,116,&
     118,120,122,124,126,128,130,132,134,136,138,140,142,144,146,148,&
     149,151,153,155,157/)

   cleb=0
   n1=i1+i2-i3+1
   if (n1.le.0) go to 99
   n2=i1-i2+i3+1
   if (n2.le.0) go to 99
   n3=-i1+i2+i3+1
   if (n3.le.0) go to 99
   it=i1+i2+i3
   if (mod(it,2).ne.0) GO TO 99
   n4=it+2
   nept=nfac(n1)+nfac(n2)+nfac(n3)-nfac(n4)
   if (nept.lt.33) go to 50
   go to 99
  50 arg=fac(n1)*fac(n2)*fac(n3)/fac(n4)
   z1=arg*10**nept
   d123=sqrt(z1)
   is=it/2
   signx=1
   if (mod(is+i3,2).eq.1) signx=-1
   ia=is+1
   ib=is-i1+1
   ic=is-i2+1
   id=is-i3+1
   nept=nfac(ia)-nfac(ib)-nfac(ic)-nfac(id)
   arg=fac(ia)/(fac(ib)*fac(ic)*fac(id))
   z1=2*i3+1
   cleb=signx*sqrt(Z1)*d123*arg*10**nept
  99 return
   end function cleb

   subroutine ngchk
   !
   !-------------------------------------------------------------------
   !  routine to read through the input ngout tape which may have data
   !  for multiple materials, multiple temperatures, Pn moments and
   !  multiple sigma-0 values and extract the materials and temperature
   !  of interest, keeping only the infinitely dilute cross section
   !  data and the p0 & p1 moments.
   !
   !  this "new" ngout tape will be assigned as scratch nscr5 (=18) and
   !  have the same mode (ascii or binary) as the original ngout tape.
   !
   !  in addition to the standard mf1/mt451 section we write out all
   !  mf3, mf5, mf6 and mf8 sections from the input gendf tape.
   !
   !  note that if multiple materials are written to the condensed
   !  gendf tape (i.e., nmt1>0), they may not appear in increasing
   !  mat order.  subsequent coding that reads/searches nscr5 has
   !  been upgraded to always search forward after rewinding the file.
   !-------------------------------------------------------------------
   !
   ! externals
   use mainio ! provides nsysi,nsyso,nsyse
   use util   ! provides timer,openz,error,repoz,closz,sigfig
   use endf   ! provides endf routines and variables

   character(70)::strng
   integer::i,i1,i2,ig,ii,iloop,j,jj,k,kk,lenscrl,lgscr,mloop
   integer::mftest
   integer::ngmode
   integer::ns0,ntw,nng,ngg,nwl,newnwl
   integer::nl,nlnew,nz,nznew,lrflag,ng,ng2,mfnow
   integer::nt,ntwds
   integer::matds,nb,nbsave,nw,nwsave
   real(kr),dimension(17)::t
   real(kr),dimension(6)::c1,c2
   real(kr),dimension(:),allocatable::gscr,scr1,scrl
   integer,dimension(:),allocatable::matsofar

   matds=matd
   if (nmt1.gt.0) then
      if (allocated(matsofar)) deallocate(matsofar)
      allocate(matsofar(nmt1+1))
      matsofar(1:nmt1+1)=0
   endif

   nscr5=18
   if (ngout.lt.0) nscr5=-18
   call openz(nscr5,1)

   !--first verify that this is a gendf-formatted tape, then
   !  look for the user specified material, matd.
   call repoz(ngout)
   nsh=0
   iloop=0
   mloop=0
   call tpidio(ngout,nscr5,0,t,nt,ntwds)
  100 continue
   iloop=iloop+1
   call contio(ngout,0,0,c1,nb,nw)
   if (iloop.eq.1.and.nint(c1(5)).ne.-1) then
      write(strng,'("ngout = ",i2," is not a gendf tape")')ngout
      call error('ngchk',strng,'')
   endif
   if (math.eq.-1) then
      write(strng,'("can''t find matd ",i4," on the gendf tape")')matd
      call error('ngchk',strng,'')
   endif
   if (math.eq.0) goto 100
   if (math.ne.matd) then
      call tomend(ngout,0,0,t)
      iloop=0
      goto 100
   endif

   !--now look for user requested temperature ... if the value read
   !  from ngout is within 0.1 degK of tempin, for temperatures less
   !  than 10,000K, or within 0.9999*tempin to 1.0001*tempin for
   !  higher temperatures, assume we've found what the user wanted.
   if (allocated(scr1)) deallocate(scr1)
   allocate(scr1(npage+50))
  200 continue
   if (math.eq.-1) then
      write(strng,'("can''t find tempin = ",1pe10.4,&
            &"K on the gendf tape")')tempin
      call error('ngchk',strng,'')
   endif
   call listio(ngout,0,0,scr1(1),nb,nw)
   nbsave=nb
   nwsave=nw
   if (tempin.eq.0 .and. scr1(1).le.0.1) goto 300
   if (tempin.lt.1.e4 .and. abs(tempin-scr1(1)).gt.0.1) then
      call tomend(ngout,0,0,c1)
      call contio(ngout,0,0,c1,nb,nw)
      goto 200
   endif
   if (tempin.ge.1.e4 .and. &
       ((scr1(1)/tempin).lt.0.9999 .or. &
        (scr1(1)/tempin).gt.1.0001)) then
      call tomend(ngout,0,0,c1)
      goto 200
   endif

   !--found material and temperature, now want to create a new gendf
   !  tape with just this material, just this temperature, just the
   !  infinitely dilute sig0 data and just the p0 and p1 data.
  300 continue
   ns0=nint(c1(4))
   ntw=nint(c1(6))
   if (ns0.ne.1) c1(4)=1
   nsh=1
   call contio(0,nscr5,0,c1,nb,nw)
   nng=nint(scr1(3))
   ngg=nint(scr1(4))
   nwl=nint(scr1(5))
   newnwl=nint(c1(4))+nint(c1(6))+nng+2
   if (ngg.ne.0) scr1(4)=0
   if (newnwl.ne.nwl)scr1(5)=newnwl
   if (allocated(scrl)) deallocate(scrl)
   allocate(scrl(nwl+6))
   scrl(1:nwsave)=scr1(1:nwsave)
   i=1
   if (nbsave.ne.0) then
      nb=nbsave
      nw=nwsave
      do while (nb.ne.0)
         i=i+nw
         call moreio(ngout,0,0,scr1(1),nb,nw)
         scrl(i:i-1+nw)=scr1(1:nw)
      enddo
   endif

   if (allocated(gscr)) deallocate(gscr)
   lgscr=newnwl+6
   allocate(gscr(lgscr))
   gscr(1:6)=scrl(1:6)                                 !list header
   gscr(7:6+ntw)=scrl(7:6+ntw)                         !title words
   gscr(7+ntw)=scrl(7+ntw)                             !first sig-0 only
   gscr(8+ntw:8+ntw+nng)=scrl(7+ntw+ns0:7+ntw+ns0+nng) !n-group e mesh
   gscr(9+ntw+nng)=0                                   !no g-group mesh
   call listio(0,nscr5,0,gscr(1),nb,nw)
   i=1
   do while (nb.ne.0)
      i=i+nw
      call moreio(0,nscr5,0,gscr(i),nb,nw)
   enddo
   call afend(nscr5,0)            !mf1/mt451 done, but no send card used

   call contio(ngout,0,0,c1,nb,nw) !read mf1/mt451 fend record
   call contio(ngout,0,0,c1,nb,nw) !first cont record of next mf
   mfnow=mfh

   !--finished with the new mf1/mt451 section, now start looping
   !  over possible mf3, mf5, mf6 and mf8 sections
   mftest=3
  400 continue
   lenscrl=npage+6
   if (allocated(scrl)) deallocate(scrl)
   allocate(scrl(lenscrl))
   do while (mfnow.eq.mftest)
      nl=nint(c1(3))
      nlnew=nl
      if (nl.gt.1) nlnew=1
      if (nl.gt.1 .and. mfh.ne.3) nlnew=2
      nz=nint(c1(4))
      nznew=1
      lrflag=nint(c1(5))
      ng=nint(c1(6))
      gscr(1:6)=c1(1:6)
      gscr(3)=nlnew
      gscr(4)=nznew
      call contio(0,nscr5,0,gscr(1),nb,nw)
      ig=1
      do while (ig.lt.ng)
         call listio(ngout,0,0,scr1,nb,nw)
         ng2=nint(scr1(3))
         nwl=nint(scr1(5))+6
         ig=nint(scr1(6))
         if (nb.eq.0) then
            scrl(1:nwl)=scr1(1:nwl)
         else
            if (nwl.gt.lenscrl) then
               deallocate(scrl)
               lenscrl=nwl
               allocate(scrl(lenscrl))
            endif
            scrl(1:nw)=scr1(1:nw)
            i1=1
            i2=nw
            do while (nb.ne.0)
               call moreio(ngout,0,0,scr1,nb,nw)
               i1=i2+1
               i2=i2+nw
               scrl(i1:i2)=scr1(1:nw)
            enddo
         endif
         newnwl=nlnew*nznew*ng2
         if ((newnwl+6).gt.lgscr) then
            deallocate(gscr)
            lgscr=newnwl+6
            allocate(gscr(lgscr))
         endif
         gscr(1:6)=scrl(1:6)
         if (nint(gscr(5)).ne.newnwl) gscr(5)=newnwl
         j=6
         do ii=7,nwl,nl*nz
            do kk=1,nlnew
               j=j+1
               k=ii+kk-1
               gscr(j)=scrl(k)
            enddo
         enddo
         call listio(0,nscr5,0,gscr,nb,nw)
         i1=1
         do while (nb.ne.0)
            i1=i1+nw
            call moreio(0,nscr5,0,gscr(i1),nb,nw)
         enddo
      enddo
      call asend(nscr5,0)
      call contio(ngout,0,0,scr1,nb,nw) !should be send
      if (mfh.ne.mfnow.or.mth.ne.0) then
         call error('ngchk','ngout records out of order','')
      endif
      call contio(ngout,0,0,c1,nb,nw) !next mf/mt
      mfnow=mfh
   enddo

   !-- finished with mf3, moving on ...
   if (mftest.lt.5) then
      mftest=5
      goto 400
   endif

   !-- finished with mf3 & mf5, moving on ...
   if (mftest.lt.6) then
      mftest=6
      goto 400
   endif

   !-- finished with mf3, mf5 & mf6, moving on ...
   if (mftest.lt.8) then
      mftest=8
      goto 400
   endif

   !-- finished with all relevant mf's for this material ...
   !   -  write material end record on the condensed GENDF tape
   call amend(nscr5,0)

   if (nmt1.ne.0) then
   !-- rewind the input GENDF tape
      call repoz(ngout)
  500 continue
      mloop=mloop+1
      matsofar(mloop)=matd
      if (mloop.le.nmt1) then
      !-- check nmt1 list for more materials ... but don't process
      !   any matn more than once.
         matd=mats(mloop)
         do i=1,mloop
            if (matd.eq.matsofar(i))goto 500
         enddo
         iloop=0
         call tpidio(ngout,0,0,t,nt,ntwds)
         goto 100
      endif
   endif

   !-- all done ...
   !   - restore original matd
   !   - write tape end record
   matd=matds
   call atend(nscr5,0)

   !-- rewind the original and condensed gendf tapes
   call repoz(ngout)
   call repoz(nscr5)

   !-- deallocate scratch arrays
   deallocate(scr1)
   deallocate(scrl)
   deallocate(gscr)
   if (allocated(matsofar)) deallocate(matsofar)

   return
   end subroutine ngchk

end module errorm


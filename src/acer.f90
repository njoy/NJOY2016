module acem
   ! provides acer for NJOY2016
     use locale
     implicit none

contains

   subroutine acer
   !--------------------------------------------------------------------
   !
   ! Prepare a data library for MCNP,
   ! the Los Alamos continuous energy Monte Carlo code.
   !
   !    --- continuous (fast) data ---
   !
   ! Reaction cross sections are reconstructed on the grid of the
   ! total cross section from the input pendf tape (assumed to be
   ! linearized and unionized).  Redundant reactions (except for
   ! MT1, MT452, and reactions needed for photon yields) are
   ! removed.  MT18 is considered redundant if MT19 is present.
   ! Angular distributions are converted into either 32 equally
   ! probable bins, or into cummulative probability distributions.
   ! Tabulated energy distributions are converted into "law 4"
   ! probability distributions.  Analytic secondary-energy
   ! distributions are converted into their ACE formats.
   ! Coupled energy-angle distributions (File 6) are converted
   ! into ACE laws.  The old format supports law44 for tabulated
   ! data with Kalbach systematics, law67 for angle-energy data,
   ! and law66 for phase space.  The newer format adds law61 with
   ! with cummulative angle distributions for Legendre or tabulated
   ! distributions (see newfor).  All photon production cross
   ! sections are combined on the cross section energy grid.
   ! If provided, multigroup photon production data is summed
   ! and converted into a set of equally probable emission
   ! energies for each input group.  Detailed photon production
   ! data can be generated directly from Files 12, 13, 14, 15,
   ! and 16 from the input ENDF tapes and written out using the
   ! "law 4" cummulative energy distribution format.
   !
   !    --- thermal data ---
   !
   ! The data from the pendf tape as prepared by the thermr
   ! module is read in.  Inelastic and incoherent elastic cross
   ! sections are stored directly.  Coherent elastic cross
   ! sections are converted to a cummulative "stair step" form
   ! and stored.  The angular representation for incoherent
   ! elastic is stored directly.  None is needed for coherent
   ! elastic.  The incoherent inelastic energy distributions
   ! are converted into probability bins with the equally
   ! probable angles left unchanged.  The bins can have equal
   ! probabilities or variable probabilities.  In the latter
   ! case, outlying bins with smaller probabilities are provided
   ! to extend the sampling to rare events.  A new tabulated option
   ! uses a continuous tabulated probability distribution (pdf/cdf)
   ! (requires a MCNP5.1.50 or later) and provides extended
   ! plotting.
   !
   !    --- dosimetry data ---
   !
   ! ENDF cross sections for dosimetry reactions are simply
   ! stored in ACE format without changing the energy grid.
   ! The endf interpolation law is also provided.
   !
   !    --- photoatomic data ---
   !
   ! Photon interaction cross sections are stored on the grid of
   ! the total cross section.  The coherent form factor is
   ! stored together with an integral over the form factor that
   ! is used in sampling for coherent scattering.  The
   ! incoherent scattering function is simply stored.  Photon
   ! heating is calculated from the incoherent scattering data,
   ! the pair production data, and the photoelectric absorption
   ! data.  The input photoatomic data file is mounted on nendf.
   ! Fluorescence data can be generated from an atomic relaxation
   ! data file in ENDF format mounted on npend.
   !
   !    --- photonuclear data ---
   !
   ! Photonuclear data are processed from new evaluations now
   ! available in ENDF format using a new ACE format developed
   ! for MCNP and MCNPX.
   !
   !    --- particle production ---
   !
   ! With the new format (see newfor), for charged particles, and
   ! for photonuclear data, new sections are written describing the
   ! distributions of light particles produced that are different
   ! from the incident particle.
   !
   !    --- incident charged particles ---
   !
   ! Incident charged particles are automatically recognized from
   ! the input tape.
   !
   !    --- mcnpx format ---
   !
   ! MCNPX format is given to support a proposed extension of
   ! the zaid indentifier that uses three digits and two letters
   ! to the right of the decimal.  This will increase flexibility
   ! for handling exotic particles and provide more space for
   ! different data versions.  To request MCNPX format, set the
   ! value of iopt negative.
   !
   !    --- output ---
   !
   ! The ACE output file can be Type 1 (formatted) or Type 2
   ! (f77 binary).  Type 3 is no longer used.  A line of file
   ! directory information is also written.  It must normally
   ! be edited to tell the system the path to the file.  ACER
   ! can also be used to print, edit, or convert the mode of
   ! existing ACE-format files.
   !
   !---input specifications (free format)---------------------------
   !
   ! card 1
   !    nendf    unit for input endf tape
   !    npend    unit for input pendf tape
   !    ngend    unit for input multigroup photon data
   !    nace     unit for output ace tape
   !    ndir     unit for output mcnp directory
   ! card 2
   !    iopt     type of acer run option
   !               1   fast data
   !               2   thermal data
   !               3   dosimetry data
   !               4   photo-atomic data
   !               5   photo-nuclear data
   !               7   read type 1 ace files to print or edit
   !               8   read type 2 ace files to print or edit
   !                set iopt negative for mcnpx format
   !    iprint   print control (0 min, 1 max, default=1)
   !    itype    ace output type (1, 2, or 3, default=1)
   !    suff     id suffix for zaid (default=.00)
   !    nxtra    number of iz,aw pairs to read in (default=0)
   ! card 3
   !    hk       descriptive character string (70 char max)
   !             delimited by quotes
   ! card 4 (nxtra.gt.0 only)
   !    iz,aw    nxtra pairs of iz and aw
   !
   !    --- fast data (iopt=1 only) ---
   !
   ! card 5
   !    matd     material to be processed
   !    tempd    temperature desired (kelvin) (default=300)
   ! card 6
   !    newfor   use new cummulative angle distributions,
   !               law 61, and outgoing particle distributions.
   !               (0=no, 1=yes, default=1)
   !    iopp     detailed photons (0=no, 1=yes, default=1)
   !    ismooth  switch on/off smoothing operation (1/0, default=1=on)
   !               set ismooth to 1 to cause extension of mf6 cm
   !               distributions to lower energies using a sqrt(E)
   !               shape, to extend delayed neutron distributions as
   !               sqrt(E) to lower energies, and to add additional
   !               points above 10 Mev to some fission spectra assuming
   !               an exponential shape.  otherwise, use ismooth=0.
   !               NOTE:  ismooth=0 is the default value in njoy99.
   ! card 7
   !  type of thinning is determined by sign of thin(1)
   !  (pos. or zero/neg.=energy skip/integral fraction)
   !  (all entries defaulted=no thinning)
   !    thin(1)  e1 energy below which to use all energies (ev)
   !             or iwtt weighting option (1=flat,2=1/e)
   !             (1/e actually has weight=10 when e lt .1)
   !    thin(2)  e2 energy above which to use all energies
   !             or target number of points
   !    thin(3)  iskf skip factor--use every iskf-th energy
   !             between e1 and e2, or rsigz reference sigma zero
   !
   !   --- thermal data (iopt=2 only) ---
   !
   ! card 8
   !    matd     material to be processed
   !    tempd    temperature desired (kelvin) (default=300)
   !    tname    thermal zaid name ( 6 char max, def=za)
   ! card 8a
   !    iza01    moderator component za value
   !    iza02    moderator component za value (def=0)
   !    iza03    moderator component za value (def=0)
   ! card 9
   !    mti      mt for thermal incoherent data
   !    nbint    number of bins for incoherent scattering
   !    mte      mt for thermal elastic data
   !    ielas    0/1=coherent/incoherent elastic
   !    nmix     number of atom types in mixed moderator
   !             (default=1, not mixed)
   !             (example, 2 for beo or c6h6)
   !    emax     maximum energy for thermal treatment (ev)
   !             (default=1000.=determined from mf3, mti)
   !    iwt      weighting option
   !             0/1/2=variable/constant/tabulated (default=variable)
   !
   !   --- dosimetry data (iopt=3 only) ---
   !
   ! card 10
   !    matd     material to be processed
   !    tempd    temperature desired (kelvin) (default=300)
   !
   !   --- photo-atomic data (iopt=4 only) ---
   !
   ! card 11
   !    matd     material to be processed
   !             photoatomic data on nendf
   !             atomic relaxation data on npend
   !
   !   --- photo-nuclear data (iopt=5 only) ---
   !
   ! card 11
   !    matd     material to be processed
   !
   !   --- print or edit existing files (iopt=7-9) ---
   !
   !    No additional input cards are required.  Mount the old
   !    ace tape on "npend".  The code can modify zaid, hk,
   !    the (iz,aw) list, and the type of the file.  Use suff<0
   !    to leave the old zaid unchanged.  Use just "/" on
   !    card 3 to leave the comment field hk unchanged.  Use
   !    nxtra=0 to leave the old iz,aw list unchanged.  The
   !    code can modify zaid, hk, and type of file.
   !
   !    Exhaustive consistency checks are automatically made on
   !    the input file.  If ngend.ne.0, a set of standard ACE plots
   !    are prepared on unit ngend as plotr input instructions.
   !
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use util ! provides timer,error
   use acefc ! provides acetop,acefix
   use aceth ! provides acesix,thrfix
   use acepa ! provides acepho,phofix
   use acepn ! provides acephn,phnfix
   use acedo ! provides acedos,dosfix

   ! internals
   integer::nendf,npend,ngend,nace,ndir
   integer::iopt,iprint,itype,nxtra
   integer::matd
   real(kr)::tempd
   integer::newfor,iopp,ismooth
   real(kr)::thin(4)
   integer::iskf,iwtt,npts
   real(kr)::suff
   character(70)::hk=' '
   integer::izn(16)
   real(kr)::awn(16)
   character(6)::tname,tscr
   integer::iza01,iza02,iza03
   integer::mti,nbint,mte,ielas,nmix,iwt
   real(kr)::emax
   real(kr)::time,zaid
   character(13)::hz
   character(3)::ht
   character(9)::str
   integer::mcnpx
   integer::i,nch
   real(kr)::z(3)
   real(kr),parameter::eps=1.e-5_kr
   real(kr),parameter::zero=0
   real(kr),parameter::one=1
   real(kr),parameter::two=2

   !--start
   call timer(time)
   write(nsyso,'(/&
     &'' acer...monte carlo neutron and photon data'',&
     &26x,f8.1,''s'')') time
   write(nsyse,'(/'' acer...'',61x,f8.1,''s'')') time

   !--read general user input
   read(nsysi,*) nendf,npend,ngend,nace,ndir
   write(nsyso,'(/&
     &'' input endf unit ...................... '',i10/&
     &'' input pendf unit ..................... '',i10/&
     &'' input gendf unit ..................... '',i10/&
     &'' output ace format unit ............... '',i10/&
     &'' output directory unit ................ '',i10)')&
     nendf,npend,ngend,nace,ndir
   iprint=1
   itype=1
   suff=0
   nxtra=0
   read(nsysi,*) iopt,iprint,itype,suff,nxtra
   mcnpx=0
   if (iopt.lt.0) then
      mcnpx=1
      iopt=-iopt
   endif
   read(nsysi,*) hk
   do i=1,16
      izn(i)=0
      awn(i)=0
   enddo
   if (nxtra.gt.0) then
      read(nsysi,*) (izn(i),awn(i),i=1,nxtra)
   endif
   write(nsyso,'(/&
     &'' run type option ...................... '',i10/&
     &'' print option (0 min, 1 max) .......... '',i10/&
     &'' type of ace file ..................... '',i10)')&
     iopt,iprint,itype
   if (mcnpx.eq.1) write(nsyso,'(''    mcnpx format requested'')')

   !--read input for fast data option
   if (iopt.eq.1) then
      tempd=300
      read(nsysi,*) matd,tempd
      write(nsyso,'(&
        &'' mat to be processed .................. '',i10/&
        &'' temperature .......................... '',1p,e10.3)')&
        matd,tempd
      iopp=1
      newfor=1
      ismooth=1
      read(nsysi,*) newfor,iopp,ismooth
      write(nsyso,'(&
        &'' new formats .......................... '',i10/&
        &'' photon option ........................ '',i10/&
        &'' smoothing option ..................... '',i10)')&
        newfor,iopp,ismooth
      if (newfor.ne.0.and.newfor.ne.1) then
         call error('acer','illegal newfor.',' ')
      endif
      if (iopp.ne.0.and.iopp.ne.1) then
         call error('acer','illegal iopp.',' ')
      endif
      if (ismooth.ne.0.and.ismooth.ne.1) then
         call error('acer','illegal ismooth.',' ')
      endif
      if (iopp.eq.0) write(nsyso,&
        '(/'' photons will not be processed'')')
      if (ismooth.eq.0) write(nsyso,&
        '(/'' smoothing operation will not be performed'')')
      mte=0
      z(1)=0
      z(2)=0
      z(3)=0
      read(nsysi,*) (z(i),i=1,3)
      thin(1)=abs(z(1))
      thin(2)=z(2)
      thin(3)=z(3)
      thin(4)=0
      if (z(1).ge.zero.and.thin(2).gt.thin(1)) thin(4)=1
      if (z(1).lt.zero) thin(4)=2
      if (thin(4).eq.one) iskf=nint(thin(3))
      if (thin(4).eq.one) write(nsyso,'(&
        &'' energy at which to start thinning .... '',1p,e10.2/&
        &'' energy at which to stop thinning ..... '',e10.2/&
        &'' skip factor .......................... '',i10)')&
        thin(1),thin(2),iskf
      if (thin(4).eq.two) iwtt=nint(thin(1))
      if (thin(4).eq.two) npts=nint(thin(2))
      if (thin(4).eq.two) write(nsyso,'(&
        &'' thinning wt option (1 flat, 2 1/e) ... '',i10/&
        &'' target no. of points ................. '',i10/&
        &'' reference sigma zero ................. '',1p,e10.2)')&
        iwtt,npts,thin(3)

   !--read input for thermal data option
   else if (iopt.eq.2) then
      tempd=300
      tscr=' '
      read(nsysi,*) matd,tempd,tscr
      nch=0
      do i=1,6
         if (tscr(i:i).ne.' ') nch=i
      enddo
      tname='      '
      if (nch.gt.0) tname(7-nch:6)=tscr(1:nch)
      iza02=0
      iza03=0
      read(nsysi,*) iza01,iza02,iza03
      write(nsyso,'(&
        &'' mat to be processed .................. '',i10/&
        &'' temperature .......................... '',1p,e10.3/&
        &'' thermal name ......................... '',4x,a6/&
        &'' iza01 ................................ '',i10/&
        &'' iza02 ................................ '',i10/&
        &'' iza03 ................................ '',i10)')&
        matd,tempd,tname,iza01,iza02,iza03
      mti=0
      nbint=0
      mte=0
      ielas=0
      nmix=1
      emax=1000
      iwt=0
      read(nsysi,*) mti,nbint,mte,ielas,nmix,emax,iwt
      if (emax.lt.eps) emax=1000
      if (nbint.eq.0) nbint=16
      write(nsyso,'(&
        &'' mt incoherent ........................ '',i10/&
        &'' bins for incoherent scattering ....... '',i10/&
        &'' mt elastic ........................... '',i10/&
        &'' coherent/incoherent elastic flag ..... '',i10/&
        &'' no. atoms in mixed moderator ......... '',i10/&
        &'' max energy for thermal ............... '',f10.3/&
        &'' weight option (0 var, 1 cons, 2 tab) . '',i10)')&
        mti,nbint,mte,ielas,nmix,emax,iwt

   !--read input for dosimetry option
   else if (iopt.eq.3) then
      tempd=300
      read(nsysi,*) matd,tempd

   !--read input for photo-atomic option
   else if (iopt.eq.4) then
      read(nsysi,*) matd
      tempd=0

   !--read input for photo-nuclear option
   else if (iopt.eq.5) then
      read(nsysi,*) matd
      tempd=0

   !--no input required for fix option
   else if (iopt.lt.7.or.iopt.gt.8) then
      call error('acer','illegal iopt.',' ')
   endif

   !--prepare fast ace data
   if (iopt.eq.1) then
      call acetop(nendf,npend,ngend,nace,ndir,iprint,itype,mcnpx,suff,&
        hk,izn,awn,matd,tempd,newfor,iopp,ismooth,thin)

   !--prepare thermal ace data
   else if (iopt.eq.2) then
      call acesix(npend,nace,ndir,matd,tempd,tname,suff,hk,izn,awn,&
        iza01,iza02,iza03,&
        mti,nbint,mte,ielas,nmix,emax,iwt,iprint,mcnpx)

   !--prepare dosimetry data
   else if (iopt.eq.3) then
      call acedos(matd,tempd,npend,nace,ndir,itype,iprint,mcnpx,&
        suff,hk,izn,awn)

   !--prepare photo-atomic data
   else if (iopt.eq.4) then
      call acepho(nendf,npend,nace,ndir,matd,mcnpx,iprint,itype,&
        suff,hk,izn,awn)

   !--prepare photo-nuclear data
   else if (iopt.eq.5) then
      call acephn(nendf,npend,nace,ndir,matd,tempd,iprint,&
        mcnpx,itype,suff,hk,izn,awn)

   !--print or edit ace files
   else if (iopt.ge.7.and.iopt.le.8) then
      itype=iopt-6

      ! read zaid to figure out file type
      if (itype.eq.1) then
         call openz(npend,0)
         if (mcnpx.eq.0) then
            read(npend,'(a10)') hz(1:10)
         else
            read(npend,'(a13)') hz(1:13)
         endif
         call closz(npend)
      else if (itype.eq.2) then
         if (mcnpx.eq.0) then
            read(npend) hz(1:10)
         else
            read(npend) hz(1:13)
         endif
      endif
      if (mcnpx.eq.0) then
         read(hz,'(a9,a1)') str,ht
         if (ht(1:1).ne.'t') then
            read(hz,'(f9.0,a1)') zaid,ht
         endif
      else
         read(hz,'(f10.0,a3)') zaid,ht
      endif

      !--call the appropriate fix routine for this file type
      if (ht(1:1).eq.'c') then
         call acefix(npend,itype,nace,ndir,iprint,ngend,&
           suff,nxtra,hk,izn,awn,mcnpx)
      else if (ht(1:1).eq.'h'.or.ht(1:1).eq.'o'.or.ht(1:1).eq.'r'&
        .or.ht(1:1).eq.'s'.or.ht(1:1).eq.'a') then
         call acefix(npend,itype,nace,ndir,iprint,ngend,&
           suff,nxtra,hk,izn,awn,mcnpx)
      else if (ht(1:1).eq.'t') then
         call thrfix(itype,npend,nace,ndir,iprint,ngend,mcnpx,&
           suff,nxtra,hk,izn,awn)
      else if (ht(1:1).eq.'p') then
         call phofix(itype,npend,nace,ndir,iprint,ngend,mcnpx,&
           suff,nxtra,hk,izn,awn)
      else if (ht(1:1).eq.'u') then
         call phnfix(itype,npend,nace,ndir,iprint,ngend,mcnpx,&
           suff,nxtra,hk,izn,awn)
      else if (ht(1:1).eq.'y') then
         call dosfix(itype,npend,nace,ndir,iprint,ngend,mcnpx,&
           suff,nxtra,hk,izn,awn)
      endif
   endif

   !--acer is finished.
   call timer(time)
   write(nsyso,'(/69x,f8.1,''s''/&
     &1x,7(''**********''),''*******'')') time
   return
   end subroutine acer

end module acem


module heatm
   ! provides heatr for NJOY2016
   use locale
   implicit none
   private
   public heatr

   ! global variables

   ! units for heatr
   integer::nendf,nin,nout,nplot

   ! heatr global parameters
   integer,parameter::nbuf=1000
   integer,parameter::nqamax=30
   integer,parameter::ilmax=200
   integer,parameter::maxmf6=320

   ! global variables
   integer::matd,mgam,npk,mtp(50),nqa,mta(nqamax),mt303,mt19
   integer,parameter::maxmf5=500
   integer,parameter::nrmax=2000
   integer::nrt
   real(kr)::ert(nrmax),rrt(nrmax)
   ! recoil of the mt5 residuals without distribution
   ! (D. Lopez Aldama, IAEA, dlopezaldama@gmail.com, Oct 2026)
   integer::n5t
   real(kr)::e5t(nrmax),d5t(nrmax),h5t(nrmax)
   integer::n5mt,mt5l(maxmf5)
   ! File 4 sections and angular tables for the anisotropy correction
   ! of the multi-neutron recoil (anitab, anipz) (DLA, IAEA, Oct 2026)
   integer::n4mt,mt4l(maxmf5)
   integer::n4e,lct4,iso4
   real(kr),dimension(:),allocatable::e4,f14,mu4,w4
   integer,dimension(:),allocatable::law4,k4o,k4n
   ! representation of the File 4 distribution at each energy (for the
   ! evaluation at arbitrary cosines): kind4=1 Legendre (a_l in rp4),
   ! kind4=2 table (mu,f pairs in rp4, panel laws in lp4); s4 = norm
   real(kr),dimension(:),allocatable::rp4,s4
   integer,dimension(:),allocatable::lp4,kind4,r4o,r4n
   ! 64-point Gauss-Legendre rule on (-1,1) for anitab/anipz
   real(kr),dimension(64),parameter::g64x=(/&
     -9.99305042E-01_kr,-9.96340117E-01_kr,-9.91013371E-01_kr,&
     -9.83336254E-01_kr,-9.73326828E-01_kr,-9.61008800E-01_kr,&
     -9.46411375E-01_kr,-9.29569172E-01_kr,-9.10522137E-01_kr,&
     -8.89315446E-01_kr,-8.65999398E-01_kr,-8.40629296E-01_kr,&
     -8.13265315E-01_kr,-7.83972359E-01_kr,-7.52819907E-01_kr,&
     -7.19881850E-01_kr,-6.85236313E-01_kr,-6.48965471E-01_kr,&
     -6.11155355E-01_kr,-5.71895646E-01_kr,-5.31279464E-01_kr,&
     -4.89403146E-01_kr,-4.46366017E-01_kr,-4.02270158E-01_kr,&
     -3.57220158E-01_kr,-3.11322872E-01_kr,-2.64687162E-01_kr,&
     -2.17423644E-01_kr,-1.69644420E-01_kr,-1.21462819E-01_kr,&
     -7.29931218E-02_kr,-2.43502927E-02_kr, 2.43502927E-02_kr,&
      7.29931218E-02_kr, 1.21462819E-01_kr, 1.69644420E-01_kr,&
      2.17423644E-01_kr, 2.64687162E-01_kr, 3.11322872E-01_kr,&
      3.57220158E-01_kr, 4.02270158E-01_kr, 4.46366017E-01_kr,&
      4.89403146E-01_kr, 5.31279464E-01_kr, 5.71895646E-01_kr,&
      6.11155355E-01_kr, 6.48965471E-01_kr, 6.85236313E-01_kr,&
      7.19881850E-01_kr, 7.52819907E-01_kr, 7.83972359E-01_kr,&
      8.13265315E-01_kr, 8.40629296E-01_kr, 8.65999398E-01_kr,&
      8.89315446E-01_kr, 9.10522137E-01_kr, 9.29569172E-01_kr,&
      9.46411375E-01_kr, 9.61008800E-01_kr, 9.73326828E-01_kr,&
      9.83336254E-01_kr, 9.91013371E-01_kr, 9.96340117E-01_kr,&
      9.99305042E-01_kr/)
   real(kr),dimension(64),parameter::g64w=(/&
      1.78328072E-03_kr, 4.14703326E-03_kr, 6.50445797E-03_kr,&
      8.84675983E-03_kr, 1.11681395E-02_kr, 1.34630479E-02_kr,&
      1.57260305E-02_kr, 1.79517158E-02_kr, 2.01348232E-02_kr,&
      2.22701738E-02_kr, 2.43527026E-02_kr, 2.63774697E-02_kr,&
      2.83396726E-02_kr, 3.02346571E-02_kr, 3.20579284E-02_kr,&
      3.38051618E-02_kr, 3.54722133E-02_kr, 3.70551285E-02_kr,&
      3.85501532E-02_kr, 3.99537411E-02_kr, 4.12625632E-02_kr,&
      4.24735151E-02_kr, 4.35837245E-02_kr, 4.45905582E-02_kr,&
      4.54916279E-02_kr, 4.62847966E-02_kr, 4.69681828E-02_kr,&
      4.75401657E-02_kr, 4.79993886E-02_kr, 4.83447622E-02_kr,&
      4.85754674E-02_kr, 4.86909570E-02_kr, 4.86909570E-02_kr,&
      4.85754674E-02_kr, 4.83447622E-02_kr, 4.79993886E-02_kr,&
      4.75401657E-02_kr, 4.69681828E-02_kr, 4.62847966E-02_kr,&
      4.54916279E-02_kr, 4.45905582E-02_kr, 4.35837245E-02_kr,&
      4.24735151E-02_kr, 4.12625632E-02_kr, 3.99537411E-02_kr,&
      3.85501532E-02_kr, 3.70551285E-02_kr, 3.54722133E-02_kr,&
      3.38051618E-02_kr, 3.20579284E-02_kr, 3.02346571E-02_kr,&
      2.83396726E-02_kr, 2.63774697E-02_kr, 2.43527026E-02_kr,&
      2.22701738E-02_kr, 2.01348232E-02_kr, 1.79517158E-02_kr,&
      1.57260305E-02_kr, 1.34630479E-02_kr, 1.11681395E-02_kr,&
      8.84675983E-03_kr, 6.50445797E-03_kr, 4.14703326E-03_kr,&
      1.78328072E-03_kr/)
   integer::ne,kchk,iprint,lqs(nqamax)
   real(kr),dimension(:),allocatable::qbar
   real(kr)::qa(nqamax),efirst,elast,za,awr,elist(ilmax)
   real(kr)::break
   real(kr)::qdel,etabmax
   integer::mt103,mt104,mt105,mt106,mt107,mt16
   integer::miss4(250),nmiss4
   integer::i6,mt6(maxmf6),i6g,mt6no(maxmf6),mt6bad(maxmf6)
   integer::i6p,mt6yp(maxmf6)
   integer::jp,jpn,jpp
   real(kr)::q,zat,awrt,zap,awp
   integer::lct
   integer::idame
   real(kr),parameter::ysing=1.01e0_kr
   real(kr)::ebot,etop
   integer::mt458,nply,lfc,ifc1,ifc2,ifc3,ifc4,ifc5,ifc6
   real(kr),dimension(:),allocatable::c458,cpoly,hpoly,afr,anp,agp
   ! fission damage: recoil of the fissioning nucleus
   ! (mass awr, charge z) from the emission of one neutron with the
   ! spectrum of the data (NJOY2016 convention, File 5), applied in the same
   ! way to File 6 without recoil (neutron subsection only; ternary p,
   ! t, alpha and photons do not contribute).  Fission fragments are
   ! not included.  Table for File 6: ef6, df6 (nf6 points).
   ! (DLA, IAEA, Oct 2026)
   integer::nf6
   real(kr),dimension(:),allocatable::ef6,df6
   real(kr)::emc2,tm,rtm

   ! exact recoil of several independently emitted
   ! particles (fnadd, cvcmp, cvemit, cvdam, dirnod): fine-node buffer,
   ! grid of the compression (cvnb x cvnb bins), azimuths (cvnphi)
   ! (DLA, IAEA, Oct 2026)
   integer,parameter::cvnb=12,cvnphi=16,cvmx=2*cvnb*cvnb
   integer::nfn
   real(kr),dimension(:),allocatable::fnz,fnt,fnw
   ! exact path of tabdmn and anadmn (neutron-only reactions with
   ! Files 4 and 5; set by conbar, never for fission)
   logical::lcvx=.false.

   ! target
   integer::izat

   ! projectile
   integer::izap

contains

   subroutine heatr
   !-------------------------------------------------------------------
   !
   ! Compute heating kerma (kinetic energy release in material)
   ! and radiation damage energy production.
   !
   ! The prompt kerma is computed pointwise on the grid of the
   ! total cross section from the input pendf tape and written
   ! onto the output PENDF tape at infinite dilution using the
   ! 300 series of MT numbers.  All temperatures on the input PENDF
   ! tape for the desired material are processed.  The dictionary
   ! is revised.  Reaction Q values are obtained from the ENDF
   ! tape unless the user enters his own value.  Partial kermas
   ! can be requested for self-shielding calculations or other
   ! purposes.  The code uses the energy balance method where
   ! photon files are available and deposits all photon energy
   ! locally when files are not available.  This assures
   ! consistency between neutron heating and energy deposition by
   ! subsequent photon interactions.  An exception is made for
   ! capture where recoil is computed by momentum conservation.
   ! Photon files are used to estimate the average photon momentum
   ! when available.  A diagnostic message is printed if the
   ! momentum calculation leads to a significant error in
   ! energy conservation.
   !
   ! If desired, the energy-balance kerma factors can be compared
   ! with conservative kinematic limits (set iprint=2).
   ! A plot file for viewr can be automatically prepared.
   !
   ! Damage energy is computed using the Lindhard electronic
   ! screening damage function with a displacement threshold
   ! from a table of default values for important elements
   ! or a value provided by the user.
   !
   !---input specifications (free format)---------------------------
   !
   ! card 1
   !    nendf    unit for endf tape
   !    nin      unit for input pendf tape
   !    nout     unit for output pendf tape
   !    nplot    unit for graphical check output
   ! card 2
   !    matd     material to be processed
   !    npk      number of partial kermas desired (default=0)
   !             (at most 10 with iprint=2, 34 otherwise)
   !    nqa      number of user q values (default=0)
   !    ntemp    number of temperatures to process
   !             (default=0, meaning all on pendf)
   !    local    0/1=gamma rays transported/deposited locally
   !             (default=0)
   !    iprint   print (0 min, 1 max, 2 check) (default=0)
   !    ed       displacement energy for damage
   !             (default from built-in table)
   ! card 3      for npk gt 0 only
   !    mtk      mt numbers for partial kermas desired
   !             total (mt301) will be provided automatically.
   !             partial kerma for reaction mt is mt+300
   !             and may not be properly defined unless
   !             a gamma file for mt is on endf tape.
   !             special values allowed--
   !               303   non-elastic (all but mt2)
   !               304   inelastic (mt51 thru 91)
   !               318   fission (mt18 or mt19, 20, 21, 38)
   !               401   disappearance (mt102 thru 120)
   !               442   total photon ev-barns
   !               443   total kinematic kerma (high limit)
   !             damage energy production values--
   !               444   total
   !               445   elastic (mt2)
   !               446   inelastic (mt51 thru 91)
   !               447   disappearance (mt102 thru 120)
   !          cards 4 and 5 for nqa gt 0 only
   ! card 4
   !    mta      mt numbers for users q values
   ! card 5
   !    qa       user specified q values (ev)
   !               (if qa.ge.99.e6, read in variable qbar
   !                  for this reaction)
   ! card 5a     variable qbar (for reactions with qa flag only)
   !    qbar      tab1 record giving qbar versus e (1000 words max)
   !
   !-------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use util   ! provides timer,openz,error,mess,repoz,closz
   use endf   ! provides endf routines and variables
   use physics ! provides amassn,amu,ev,clight
   ! internals
   integer::ntemp,local,npkk,i,nz,j,isave,iz,loc,nzq,nz0
   integer::nscr,nend4,nend6,iold,inew,itemp,idone,nb,nw
   real(kr)::time,flag
   ! up to 10 partial kermas with the kinematic limits
   ! (kchk=1: 3*npk+7 words), 34 without (npk+3); arrays sized 50
   ! (DLA, IAEA, Oct 2026)
   integer,parameter::npkmax=37
   integer::mtk(2+npkmax)
   real(kr)::z(17)
   real(kr),dimension(:),allocatable::tmp
   real(kr),parameter::qtest=99.e6_kr
   real(kr),parameter::qflag=-1.e-9_kr
   real(kr),parameter::zero=0
   integer,parameter::maxqbar=10000

   !--start
   call timer(time)
   write(nsyso,&
     '(/'' heatr...prompt kerma'',48x,f8.1,''s'')') time
   write(nsyse,'(/'' heatr...'',60x,f8.1,''s'')') time

   !--read user input.
   nplot=0
   read(nsysi,*) nendf,nin,nout,nplot
   call openz(nendf,0)
   call openz(nin,0)
   call openz(nout,1)
   if (nplot.ne.0) call openz(nplot,1)
   npk=0
   nqa=0
   ntemp=0
   local=0
   iprint=0
   break=0
   read(nsysi,*) matd,npk,nqa,ntemp,local,iprint,break
   kchk=0
   if (iprint.eq.2) kchk=1
   if (iprint.eq.2) iprint=1
   npkk=npk+3
   if (kchk.eq.1) npkk=3*npk+7
   if (npkk.gt.npkmax) call error('heatr',&
     'requested too many kerma mt-s.',&
     'at most 10 with the kinematic check, 34 without')
   if (nqa.gt.nqamax)&
     call error('heatr','requested too many q values.',' ')
   if (npk.gt.0) then
      read(nsysi,*) (mtk(i),i=1,npk)
   endif
   allocate(tmp(maxqbar))
   loc=1
   if (nqa.ne.0) then
      read(nsysi,*) (mta(i),i=1,nqa)
      read(nsysi,*) (qa(i),i=1,nqa)
      do i=1,nqa
         lqs(i)=0
         if (qa(i).ge.qtest) then
            lqs(i)=loc
            nz=1000
            flag=qflag
            do j=1,nz
               tmp(loc+j-1)=flag
            enddo
            read(nsysi,*) (tmp(loc+j-1),j=1,nz)
            nz0=0
            do j=1,nz
               if (tmp(loc+j-1).gt.flag) nz0=j
            enddo
            loc=loc+nz0
            if (loc.gt.maxqbar) call error('heatr',&
              'too much energy-dependent q data',' ')
         endif
      enddo
   endif
   nzq=loc-1
   if (nzq.gt.0) then
      allocate(qbar(nzq))
      do i=1,nzq
         qbar(i)=tmp(i)
      enddo
   endif
   deallocate(tmp)
   if ((nin.lt.0.and.nout.gt.0).or.(nin.gt.0.and.nout.lt.0))&
     call error('heatr',&
     'mode conversion not allowed between nin and nout.',' ')

   !--check the partial kermas
   if (npk.gt.0) then
      isave=0
      do i=1,npk
         if (mtk(i).eq.443.and.kchk.eq.0) kchk=2
         if (mtk(i).eq.301) isave=i
      enddo
      if (isave.gt.0) then
         call mess('heatr','mt301 always calculated',&
           '--you do not need to ask for it.')
         npk=npk-1
         if (isave.le.npk) then
            do i=isave,npk
               mtk(i)=mtk(i+1)
            enddo
         endif
      endif
      do i=1,npk
         mtp(npk+3-i)=mtk(i)
      enddo
   endif
   mtp(1)=0
   mtp(2)=301
   npk=npk+2
   ! there are always at least 2 elements in the npk array,
   ! (1)=energy, (2)=mt301 (total)
   if (npk.gt.3) call horder(mtp,npk)

   !--list all input parameters
   write(nsyso,'(/&
     &'' input endf unit ...................... '',i10/&
     &'' input pendf unit ..................... '',i10/&
     &'' output pendf unit .................... '',i10)')&
     nendf,nin,nout
   write(nsyso,'(&
     &'' mat to be processed .................. '',i10)') matd
   write(nsyso,'(&
     &'' no. temperatures (0=all) ............. '',i10/&
     &'' gamma heat (0 nonlocal, 1 local) ..... '',i10/&
     &'' print option (0 min, 1 more, 2 chk) .. '',i10)')&
     ntemp,local,iprint
   if (break.eq.zero) then
      write(nsyso,'(&
        &'' damage displacement energy ...........    default'')')
   else
      write(nsyso,'(&
        &'' damage displacement energy ........... '',f10.1,&
        &'' ev'')') break
   endif
   if (ntemp.eq.0) ntemp=100
   if (npk.ge.3) write(nsyso,'(&
     &'' partial kerma mt-s desired ........... '',i10)') mtp(3)
   if (npk.ge.4) write(nsyso,'(40x,i10)')(mtp(i),i=4,npk)
   if (nqa.ne.0) then
      write(nsyso,'(&
        &'' auxiliary reactions .................. '',i10)') mta(1)
      if (nqa.gt.1) write(nsyso,'(40x,i10)') (mta(i),i=2,nqa)
      write(nsyso,'(&
        &'' auxiliary q values ................... '',1p,e12.4)') qa(1)
      if (nqa.gt.1) write(nsyso,'(40x,1p,e12.4)') (qa(i),i=2,nqa)
      if (nzq.ne.0) then
         do i=1,nqa
            if (qa(i).ge.qtest) then
               nsh=0
               math=1
               mfh=1
               mth=mta(i)
               write(nsyso,'('' input q ----'')')
               call tab1io(0,nsyso,0,qbar,nb,nw)
            endif
         enddo
      endif
   endif

   !--assign scratch units
   nscr=10
   if (nendf.lt.0) nscr=-nscr
   nend4=13
   if (nendf.lt.0) nend4=-nend4
   nend6=14
   if (nendf.lt.0) nend6=-nend6
   iold=11
   inew=12
   call openz(-iold,1)
   call openz(-inew,1)
   call openz(nscr,1)
   call openz(nend4,1)
   call openz(nend6,1)

   !--process input pendf tape
   call repoz(nin)
   call repoz(nout)
   nsh=0
   call tpidio(nin,nout,0,z,nb,nw)
   call findf(matd,1,451,nin)
   itemp=0

   !--loop over desired temperatures for this mat
   idone=0
   do while (idone.eq.0)
      call contio(nin,0,0,z,nb,nw)
      if (math.ne.matd) then
         idone=1
      else
         itemp=itemp+1
         if (itemp.gt.ntemp) then
            idone=1
         else
            za=c1h
            awr=c2h
            emc2=amassn*amu*clight*clight/ev
            tm=emc2*(awr+1)
            rtm=1/tm

            !--check on default damage displacement energy
            if (break.eq.zero) then
               iz=nint(za/1000)
               if (iz.eq.4) then
                  break=31
               else if (iz.eq.6) then
                  break=31
               else if (iz.eq.12) then
                  break=25
               else if (iz.eq.13) then
                  break=27
               else if (iz.eq.14) then
                  break=25
               else if (iz.eq.20) then
                  break=40
               else if (iz.ge.22.and.iz.le.29) then
                  break=40
               else if (iz.eq.40) then
                  break=40
               else if (iz.eq.41) then
                  break=40
               else if (iz.eq.42) then
                  break=60
               else if (iz.eq.47) then
                  break=60
               else if (iz.eq.73) then
                  break=90
               else if (iz.eq.74) then
                  break=90
               else if (iz.eq.79) then
                  break=30
               else if (iz.eq.82) then
                  break=25
               else
                  break=25
               endif
               write(nsyso,'(/'' default damage energy ='',f5.1,&
                 &'' ev'')') break
            endif

            !--analyze dictionary and read energy grids.
            call hinit(iold,nend4,nend6,nscr)

            !--calculate heating from neutron files.
            call timer(time)
            write(nsyso,'(&
              &61x,''temp'',i2,2x,f8.1,''s'')') itemp,time
            call nheat(iold,inew,nscr,nend4,nend6,local)

            !--correct for gamma production.
            if (mgam.gt.0.and.mgam.ne.10.and.local.eq.0) then
               call timer(time)
               write(nsyso,'(&
                 &61x,''temp'',i2,2x,f8.1,''s'')') itemp,time
               call gheat(iold,inew,nscr)
            endif

            !--write kermas on pendf tape.
            call hout(iold)
         endif
      endif
   enddo

   !--heatr is finished.
   if (allocated(c458)) deallocate(c458)
   if (allocated(cpoly)) deallocate(cpoly)
   if (allocated(hpoly)) deallocate(hpoly)
   if (allocated(afr)) deallocate(afr)
   if (allocated(anp)) deallocate(anp)
   if (allocated(agp)) deallocate(agp)
   call atend(nout,0)
   call repoz(nout)
   call repoz(nin)
   call repoz(nendf)
   call closz(nin)
   call closz(nout)
   call closz(nendf)
   call closz(-iold)
   call closz(-inew)
   call closz(nscr)
   call closz(nend4)
   call closz(nend6)
   call timer(time)
   write(nsyso,'(/69x,f8.1,''s''/&
     &1x,7(''**********''),''*******'')') time
   return
   end subroutine heatr

   subroutine horder(iarray,n)
   !-------------------------------------------------------------------
   ! Arrange the n elements of iarray in ascending order.
   !-------------------------------------------------------------------
   ! externals
   integer::n,iarray(n)
   ! internals
   integer::i,j,isave

   do i=1,n
      do j=i,n
         if (iarray(i).gt.iarray(j)) then
            isave=iarray(i)
            iarray(i)=iarray(j)
            iarray(j)=isave
         endif
      enddo
   enddo
   return
   end subroutine horder

   subroutine hinit(iold,nend4,nend6,nscr)
   !-------------------------------------------------------------------
   ! Analyze dictionary and set necessary flags.
   ! Read energy grid of total cross section.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use util   ! provides repoz,error,mess,loada
   use endf   ! provides endf routines and variables
   ! externals
   integer::iold,nend4,nend6,nscr
   ! internals
   integer::ifiss,i519,nb,nw,n6,nx,i,ifini,ii6,lf,n,nk,nfc,ifc
   integer::mfd,mtd,j,idone,mf4,mt4,ie,nmu,imu,npkk,idis
   integer::ielem,ik,nr,idnx,nen,law
   integer::lrel,nmod
   integer::i6t,nl,l
   integer::lct6,iph102,nrec102,n4rec102,npho102,ipm102,lep6,nbad6
   integer::ie1bad,iok,nwl
   logical::chk102
   real(kr)::q102,ekick,ekmax,emx102,elim102
   real(kr)::e1bad,ek1bad,em1bad,ex1bad,el1bad
   real(kr),dimension(:),allocatable::rb
   character(60)::strng2
   real(kr)::zar,awrr,e,enext,y,yld,test,efix
   real(kr)::c(50)
   character(60)::strng
   real(kr),dimension(:),allocatable::scr
   real(kr),dimension(:),allocatable::bufo
   real(kr),parameter::fact=.99999e0_kr
   real(kr),parameter::rup=1.000000001e0_kr
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::zero=0
   real(kr),parameter::mevev=1.e-6_kr
   real(kr),parameter::tol=1.e-5_kr

   !--analyze dictionary
   lrel=0
   nmod=0
   i6=0
   mt103=0
   mt104=0
   mt105=0
   mt106=0
   mt107=0
   mt16=0
   nmiss4=0
   mgam=0
   ifiss=0
   mt19=0
   i519=0
   n5mt=0
   n4mt=0
   mt458=0
   qdel=0
   jp=0
   jpn=0
   jpp=0
   etop=20000000
   ebot=1
   ebot=ebot/100000
   etabmax=-1.0
   call repoz(nendf)
   call tpidio(nendf,0,0,c(1),nb,nw)
   call findf(matd,1,451,nendf)
   call contio(nendf,0,0,c(1),nb,nw)
   n6=n2h
   call contio(nendf,0,0,c(1),nb,nw)
   nx=0
   if (n1h.ne.0) then
      iverf=4
      nx=n6
   else if (n2h.eq.0) then
      iverf=5
   else
      iverf=6
   endif
   call skiprz(nendf,-1)
   if (iverf.ge.5) call contio(nendf,0,0,c(1),nb,nw)
   if (iverf.eq.6) call contio(nendf,0,0,c(1),nb,nw)
   if (iverf.eq.6) then
      if (c2h.gt.etop) etop=c2h
      lrel=l1h
      nmod=n2h
   endif
   allocate(scr(npage+50))
   call hdatio(nendf,0,0,scr,nb,nw)
   do while (nb.ne.0)
      call moreio(nendf,0,0,scr,nb,nw)
   enddo
   if (iverf.ge.5) nx=n2h
   deallocate(scr)
   nw=6*nx
   if (nw.lt.npage+50) nw=npage+50
   allocate(scr(nw))
   nw=nx
   if (nw.ne.0) then
      call dictio(nendf,0,0,scr,nb,nw)
      i=1
      ifini=0
      do while (i.le.nw.and.ifini.eq.0)
         mfd=nint(scr(i+2))
         mtd=nint(scr(i+3))
         if (mfd.ne.2) then
            if (mfd.eq.12.and.mtd.ne.460) mgam=1
            if (mfd.eq.13) mgam=1
            if (mfd.ge.13) then
               ifini=1
            else
               if (mfd.le.3) then
                  if (mfd.eq.1.and.mtd.eq.458) mt458=1
                  if (mtd.eq.19) mt19=1
                  if (mtd.eq.18.or.mtd.eq.19) ifiss=1
                  if (iverf.ge.6) then
                     if (mtd.ge.600.and.mtd.lt.650) mt103=1
                     if (mtd.ge.650.and.mtd.lt.700) mt104=1
                     if (mtd.ge.700.and.mtd.lt.750) mt105=1
                     if (mtd.ge.750.and.mtd.lt.800) mt106=1
                     if (mtd.ge.800.and.mtd.lt.850) mt107=1
                     ! mt891 (continuum) included (was .lt.891)
                     ! (DLA, IAEA, Oct 2026)
                     if (mtd.ge.875.and.mtd.le.891) mt16=1
                     if (mtd.ge.600.and.mtd.le.850) then
                        j=i
                        idone=0
                        do while (j.le.nw.and.idone.eq.0)
                           mf4=nint(scr(j+2))
                           mt4=nint(scr(j+3))
                           if (mf4.eq.4.and.mt4.eq.mtd) then
                              idone=1
                           else
                              if (mf4.eq.6.and.mt4.eq.mtd) then
                                 idone=1
                              else
                                 j=j+6
                              endif
                           endif
                        enddo
                        if (idone.eq.0) then
                           write(strng,'(&
                             &''mf4 and 6 missing, '',&
                             &''isotropy assumed for mt '',i3)') mtd
                           call mess('hinit',strng,' ')
                           nmiss4=nmiss4+1
                           miss4(nmiss4)=mtd
                        endif
                     endif
                  else
                     if (mtd.ge.700.and.mtd.lt.720) mt103=1
                     if (mtd.ge.720.and.mtd.lt.740) mt104=1
                     if (mtd.ge.740.and.mtd.lt.760) mt105=1
                     if (mtd.ge.760.and.mtd.lt.780) mt106=1
                     if (mtd.ge.780.and.mtd.lt.800) mt107=1
                     if (mtd.ge.700.and.mtd.le.800) then
                        j=i
                        idone=0
                        do while (j.le.nw.and.idone.eq.0)
                           mf4=nint(scr(j+2))
                           mt4=nint(scr(j+3))
                           if (mf4.eq.4.and.mt4.eq.mtd) then
                              idone=1
                           else
                              if (mf4.eq.6.and.mt4.eq.mtd) then
                                 idone=1
                              else
                                 j=j+6
                              endif
                           endif
                        enddo
                        if (idone.eq.0) then
                           write(strng,'(&
                             &''mf4 and 6 missing, '',&
                             &''isotropy assumed for mt '',i3)') mtd
                           call mess('hinit',strng,' ')
                           nmiss4=nmiss4+1
                           miss4(nmiss4)=mtd
                        endif
                     endif
                  endif
               else
                  if (mfd.eq.5) then
                     if (mtd.eq.19) i519=1
                     if (n5mt.lt.maxmf5) then
                        n5mt=n5mt+1
                        mt5l(n5mt)=mtd
                     endif
                  else if (mfd.eq.4) then
                     if (n4mt.lt.maxmf5) then
                        n4mt=n4mt+1
                        mt4l(n4mt)=mtd
                     endif
                  else if (mfd.eq.6) then
                     i6=i6+1
                     if (i6.gt.maxmf6) call error('hinit',&
                       'too many mf6 reactions',' ')
                     mt6(i6)=mtd
                  endif
               endif
            endif
         endif
         i=i+6
      enddo
   endif
   if (mt19.eq.1.and.i519.eq.0) mt19=2
   if (mt19.eq.1) call mess('hinit','mt18 is redundant',&
     '19, 20, 21, ... will be used.')
   if (mt19.eq.2) call mess('hinit','mt19 has no spectrum',&
     'mt18 spectrum will be used.')

   !--calculate q change for fissionables in mt 458
   if (ifiss.ne.0) then
      deallocate(scr)
      nw=10000
      allocate(scr(nw))
      if (allocated(c458)) deallocate(c458)
      if (allocated(cpoly)) deallocate(cpoly)
      if (allocated(hpoly)) deallocate(hpoly)
      if (allocated(afr)) deallocate(afr)
      if (allocated(anp)) deallocate(anp)
      if (allocated(agp)) deallocate(agp)
      if (mt458.eq.1) then
         call findf(matd,1,458,nendf)
         call contio(nendf,0,0,scr,nb,nw)
         lfc=nint(scr(4))
         nfc=nint(scr(6))
         ifc1=0
         ifc2=0
         ifc3=0
         ifc4=0
         ifc5=0
         ifc6=0
         l=1
         call listio(nendf,0,0,scr(l),nb,nw)
         l=l+nw
         do while (nb.ne.0)
            call moreio(nendf,0,0,scr(l),nb,nw)
            l=l+nw
         enddo
         nply=nint(scr(4))
         if (lfc.eq.0) then
            allocate(cpoly(0:nply))
            allocate(hpoly(0:nply))
            cpoly=0
            hpoly=0
         endif
         allocate(c458(nint(scr(5))))
         ! save polynomial coefficients ... c0 & c1 terms are ok,
         ! but for endf/b-vii.1, c2 and higher are mistakenly given
         ! as if the energy will be given in MeV and therefore need
         ! correction.
         write(nsyso,'(/,&
           &'' fission energy components''/)')
         c458(1:18)=scr(7:24)    ! scr is a list record
         if (lfc.eq.0) then      ! thermal point or polynomial
            qdel=c458(5)+c458(9)+c458(11) ! delayed fission Q at 0 eV
            if (nply.ge.1) then
               c458(19:36)=scr(25:42)
               if (nply.gt.1) then
                  do n=2,nply
                     if (nmod.ne.7.or.lrel.ne.1) then ! endf/b-vii.1 correction
                        efix=1
                     else
                        efix=mevev**(n-1)
                     endif
                     c458(1+18*n:18+18*n)=scr(7+18*n:24+18*n)*efix
                  enddo
               endif
               write(nsyso,'(&
                 &''   fission products: polynomial of order '',i2/&
                 &''   prompt neutrons : polynomial of order '',i2/&
                 &''   prompt gammas   : polynomial of order '',i2)')&
                 nply,nply,nply
            else
               write(nsyso,'(&
                 &''   fission products: thermal point''/&
                 &''   prompt neutrons : thermal point''/&
                 &''   prompt gammas   : thermal point'')')
            endif
            do  i=0,nply
               cpoly(i)=c458(1+i*18)+c458(3+i*18)+c458(7+i*18)
               hpoly(i)=c458(1+i*18)+c458(7+i*18)
            enddo
         else if (lfc.eq.1) then ! tabulated components
            qdel=0
            do i=1,nfc
               l=1
               call tab1io(nendf,0,0,scr(l),nb,nw)
               do while (nb.ne.0)
                  l=l+nw
                  call moreio(nendf,0,0,scr(l),nb,nw)
               enddo
               ifc=l2h
               nl=6+2*n1h+2*n2h
               if (etabmax.lt.0) then
                  etabmax=scr(nl-1)
               else
                  if ((abs(etabmax-scr(nl-1))/etabmax).gt.tol) then
                     write(strng,'(''upper energy mismatch for ifc='',i2,&
                        &'' in mt=458.'')')ifc
                     call error('hinit',strng,' ')
                  endif
               endif
               if (ifc.eq.1) then      ! EFR is tabulated
                  ifc1=1
                  allocate(afr(1:nl))
                  afr(1:nl)=scr(1:nl)
               else if (ifc.eq.2) then ! ENP is tabulated
                  ifc2=1
                  allocate(anp(1:nl))
                  anp(1:nl)=scr(1:nl)
               else if (ifc.eq.3) then ! END is tabulated
                  ifc3=1
                  qdel=qdel+scr(6+2*n1h+2)
               else if (ifc.eq.4) then ! EGP is tabulated
                  ifc4=1
                  allocate(agp(1:nl))
                  agp(1:nl)=scr(1:nl)
               else if (ifc.eq.5) then ! EGD is tabulated
                  ifc5=1
                  qdel=qdel+scr(6+2*n1h+2)
               else if (ifc.eq.6) then ! EB is tabulated
                  ifc6=1
                  qdel=qdel+scr(6+2*n1h+2)
               endif
            enddo
            if (ifc3.eq.0) qdel=qdel+c458(5)
            if (ifc5.eq.0) qdel=qdel+c458(9)
            if (ifc6.eq.0) qdel=qdel+c458(11)
            if (ifc1.eq.1) then
               write(nsyso,'(''   fission products: tabulated values'')')
            else
               write(nsyso,'(''   fission products: thermal point'')')
            endif
            if (ifc2.eq.1) then
               write(nsyso,'(''   prompt neutrons : tabulated values'')')
            else
               write(nsyso,'(''   prompt neutrons : thermal point'')')
            endif
            if (ifc4.eq.1) then
               write(nsyso,'(''   prompt gammas   : tabulated values'')')
            else
               write(nsyso,'(''   prompt gammas   : thermal point'')')
            endif
            if (etabmax.lt.0) then
               write(strng,'(''no tabulated fission q components found'',&
                  &'' in mt=458.'')')
               call error('nheat',strng,' ')
            endif
         else
            call error('hinit','bad LFC in mt=458.',' ')
         endif
      else
         lfc=0
         nply=0
         allocate(cpoly(0:nply))
         cpoly=0
         call mess('hinit','mt458 is missing for this mat',' ')
      endif
   endif

   !--prepare a scratch tape for file 6
   !--check for missing recoil subsections
   !--and check for photon data
   i6p=0
   mt6yp=0
   if (i6.gt.0) then
      call repoz(nend6)
      nsp=1
      math=1
      call afend(nend6,0)
      ! q value of mt102 (File 3) for the kinematic
      ! check of File 6 capture recoil distributions (DLA, IAEA, Oct 2026)
      q102=0
      do i=1,i6
         if (mt6(i).eq.102) then
            call findf(matd,3,102,nendf)
            call contio(nendf,0,0,scr,nb,nw)
            call tab1io(nendf,0,0,scr,nb,nw)
            q102=c2h
         endif
      enddo
      call findf(matd,6,0,nendf)
      ii6=0
   endif
   mfh=6
   do while (mfh.eq.6)
      call contio(nendf,0,0,scr,nb,nw)
      if (mth.eq.18) then
         jp=l1h
         jpn=mod(jp,10)
         jpp=jp-10*jpn
      endif
      if (mth.ne.18 .or. (mth.eq.18 .and. jp.eq.0)) then
         call contio(0,nend6,0,scr,nb,nw)
      else
         !--delete mt=18 from mt6 array and skip this section (for now)
         i6t=i6
         do i=1,i6t
            if (mt6(i).lt.18) cycle
            if (mt6(i).eq.18) i6=i6-1
            if (mt6(i).gt.18 .and. i6.ne.i6t) mt6(i-1)=mt6(i)
         enddo
         call tosend(nendf,0,0,scr)
         cycle
      endif
      if (mfh.ne.0) then
         if (mth.ne.0) then
            nk=n1h
            ii6=ii6+1
            mt6no(ii6)=0
            mt6bad(ii6)=0
            lct6=nint(scr(4))
            iph102=0
            nrec102=0
            n4rec102=0
            npho102=0
            ipm102=0
            zar=c1h+1
            awrr=c2h+1
            ielem=mod(nint(c1h),1000)
            do ik=1,nk
               l=1
               call tab1io(nendf,nend6,0,scr(l),nb,nw)
               l=l+nw
               do while (nb.ne.0)
                  call moreio(nendf,nend6,0,scr(l),nb,nw)
                  l=l+nw
               enddo
               nr=n1h
               zap=c1h
               law=l2h
               if (nint(zap).eq.0) then
                  i6p=i6p+1
                  mt6yp(i6p)=mth
               endif
               lf=l1h
               if (mth.eq.5) zar=0
               if (mth.ne.5) then
                  yld=scr(6+2*nr+2)
                  zar=zar-zap*yld
                  awrr=awrr-c2h*yld
                  if (zap.eq.zero.and.(ik.gt.1.or.nk.eq.1)) then
                     if (nk.ne.1.or.zap.ne.zero) then
                        if (zar.ne.zero) then
                           test=3000
                           if (ielem.ne.0.or.abs(zar).ge.test) then
                              write(strng,'(''mf6, mt'',i3,&
                                &'' does not give recoil za='',i6)')&
                                mth,nint(zar)
                              if (mth.eq.102) then
                                 call mess('hinit',strng,&
                                   'photon momentum recoil used.')
                              else if (mth.ne.10) then
                                 call mess('hinit',strng,&
                                   'one-particle recoil approx. used.')
                              endif
                              if (mth.ne.102) mt6no(ii6)=ik-1
                              if (mth.eq.102) ipm102=1
                              zar=0
                           endif
                        endif
                     endif
                  endif
               endif
               ! the test for a recoil without distribution
               ! used lf (LIP) instead of LAW (DLA, IAEA, Oct 2026)
               if (mth.eq.102) then
                  if (zap.eq.zero) then
                     npho102=npho102+1
                     if (iph102.eq.0) iph102=ik
                  else
                     nrec102=nrec102+1
                     if (law.eq.4) n4rec102=n4rec102+1
                     if (npho102.gt.0.and.law.ne.4) then
                        write(strng,'(''mf6, mt'',i3,&
                          &'' recoil subsection after photons'')') mth
                        call mess('hinit',strng,&
                          'ENDF-6 order violated: correct evaluation')
                     endif
                     if (law.eq.0) then
                        write(strng,'(''mf6, mt'',i3,&
                         &'' has recoil with no spectrum'')') mth
                        call mess('hinit',strng,&
                                        'photon momentum recoil used.')
                        ipm102=1
                     endif
                  endif
               endif
               if (zap.eq.zero) mgam=10+mod(mgam,10)
               if (law.eq.6) then
                  call contio(nendf,nend6,0,scr,nb,nw)
               else if (law.eq.1.or.law.eq.2.or.law.eq.5) then
                  call tab2io(nendf,nend6,0,scr,nb,nw)
                  ne=n2h
                  lep6=l2h
                  chk102=(mth.eq.102.and.zap.gt.zero.and.law.eq.1)
                  nbad6=0
                  do ie=1,ne
                     call listio(nendf,nend6,0,scr,nb,nw)
                     if (chk102) then
                        ! kinematic check of the recoil (DLA, IAEA, Oct 2026)
                        nwl=nint(scr(5))+6
                        if (allocated(rb)) deallocate(rb)
                        allocate(rb(nwl+npage+50))
                        rb(1:nw)=scr(1:nw)
                        l=1+nw
                        do while (nb.ne.0)
                           call moreio(nendf,nend6,0,rb(l),nb,nw)
                           l=l+nw
                        enddo
                        call rec102(rb,lep6,q102,zap,lct6,iok,&
                          ekick,ekmax,emx102,elim102)
                        if (iok.eq.0) then
                           if (nbad6.eq.0) then
                              ie1bad=ie
                              e1bad=rb(2)
                              ek1bad=ekick
                              em1bad=ekmax
                              ex1bad=emx102
                              el1bad=elim102
                           endif
                           nbad6=nbad6+1
                        endif
                     else
                        do while (nb.ne.0)
                           call moreio(nendf,nend6,0,scr,nb,nw)
                        enddo
                     endif
                  enddo
                  if (chk102.and.nbad6.gt.0) then
                     write(strng,'(''mf6 mt'',i3,'' recoil za'',i6,&
                       &'' fails kinematic check'')') mth,nint(zap)
                     write(strng2,'(''at'',i4,'' of'',i4,&
                       &'' energies: photon momentum recoil'')')&
                       nbad6,ne
                     call mess('hinit',strng,strng2)
                     write(nsyso,'(&
                       &''   first failure at e='',1p,e11.4,'' ev:''/&
                       &''   mean kick='',e11.4,'' ev, allowed '',&
                       &e11.4,'' to '',e11.4,'' ev''/&
                       &''   max secondary energy='',e11.4,&
                       &'' ev, allowed '',e11.4,'' ev'')')&
                       e1bad,ek1bad,em1bad/20,1.1e0_kr*em1bad,&
                       ex1bad,el1bad
                     mt6bad(ii6)=ik
                     ipm102=1
                  endif
               else if (law.eq.7) then
                  call tab2io(nendf,nend6,0,scr,nb,nw)
                  ne=n2h
                  do ie=1,ne
                     call tab2io(nendf,nend6,0,scr,nb,nw)
                     nmu=n2h
                     do imu=1,nmu
                        call tab1io(nendf,nend6,0,scr,nb,nw)
                        do while (nb.ne.0)
                           call moreio(nendf,nend6,0,scr,nb,nw)
                        enddo
                     enddo
                  enddo
               else if (law.lt.0) then
                  call skip6(nendf,nend6,0,scr,law)
               endif
            enddo
            if (zar.ne.zero) then
               test=3000
               if (ielem.ne.0.or.abs(zar).ge.test) then
                  write(strng,'(''mf6, mt'',i3,&
                    &'' does not give recoil za='',i6)') mth,nint(zar)
                  if (mth.eq.102) then
                     call mess('hinit',strng,&
                       'photon momentum recoil used.')
                  else if (mth.ne.10) then
                     call mess('hinit',strng,&
                       'one-particle recoil approx. used.')
                  endif
                  if (mth.ne.102) mt6no(ii6)=nk
                  if (mth.eq.102) ipm102=1
               endif
            endif
            ! the photon-momentum recoil for mt102 is
            ! computed from the first photon subsection (it was taken
            ! from the last subsection, wrong when the recoil is last)
            ! (DLA, IAEA, Oct 2026)
            if (mth.eq.102) then
               if (n4rec102.gt.1) then
                  write(strng,'(''mf6, mt'',i3,'' has'',i3,&
                    &'' photon/recoil (law 4) pairs'')') mth,n4rec102
                  call mess('hinit',strng,&
                    'cascade processed as given: check results')
               endif
               if (ipm102.ne.0) then
                  if (iph102.gt.0) then
                     mt6no(ii6)=iph102
                     if (npho102.gt.1) then
                        write(strng,'(''mf6, mt'',i3,'' has'',i3,&
                          &'' photon subsections'')') mth,npho102
                        call mess('hinit',strng,&
                          'photon momentum from the first one only')
                     endif
                  else
                     write(strng,'(''mf6, mt'',i3,&
                       &'' has no usable recoil and no photons'')') mth
                     call mess('hinit',strng,&
                       'no capture recoil damage from File 6')
                     mt6no(ii6)=0
                  endif
               endif
            endif
         endif
      endif
   enddo
   call amend(nend6,0)
   call atend(nend6,0)
   call repoz(nend6)
   if (mgam.eq.0) write(nsyso,'(/&
     &'' no photon production files...'',&
     &''all photon energy will be deposited locally.'')')

   !--copy endf tape to nscr for use in conbar and gheat
   !--and to nend4 for use in disbar
   !--convert photon transition probability arrays if present
   call repoz(nendf)
   call repoz(nend4)
   call repoz(nscr)
   nsp=1
   nsc=1
   math=1
   call afend(nend4,nscr)
   call findf(matd,1,451,nendf)
   call hconvr(nendf,nend4,nscr)
   call amend(nend4,nscr)
   call atend(nend4,nscr)
   call repoz(nend4)
   call repoz(nscr)
   call repoz(nendf)

   !--store energy grid of the total cross section.
   call findf(matd,3,1,nin)
   npkk=npk
   if (kchk.eq.1) npkk=3*npk-2
   if (mgam.gt.0) npkk=npkk+3
   if (mgam.eq.0.and.i6.gt.0) npkk=npkk+1
   do i=1,npkk
      c(i)=0
   enddo
   call contio(nin,0,0,scr,nb,nw)
   e=0
   call gety1(e,enext,idis,y,nin,scr)
   enext=enext*rup
   efirst=enext
   idnx=0
   allocate(bufo(nbuf))
   ne=0
   nen=1
   do while (nen.gt.0)
      ne=ne+1
      test=enext*fact
      if (idnx.gt.0.and.test.gt.e) enext=test
      e=enext
      call gety1(e,enext,idis,y,nin,scr)
      nen=ne
      if (enext.gt.big-big/100) nen=-nen
      c(1)=e
      call loada(nen,c,npkk,iold,bufo,nbuf)
      if (idis.gt.0.and.idnx.gt.0) idnx=0
      if (idis.gt.0.and.idnx.eq.0) idnx=1
      if (idis.eq.0) idnx=0
   enddo
   elast=e
   call tosend(nin,0,0,scr)
   deallocate(bufo)
   deallocate(scr)
   mt303=0
   if (npk.lt.3) return
   do i=3,npk
      if (mtp(i).eq.303) mt303=i
   enddo
   return
   end subroutine hinit

   subroutine rec102(a,lep,q,zap,lct6,iok,ekick,ekmax,emx,elim)
   !-------------------------------------------------------------------
   ! Kinematic check of a File 6 LAW=1 recoil distribution
   ! for radiative capture (mt102) at one incident energy.
   ! a is the LIST record [0,E,ND,NA,NW,NEP/E'(1),b(0,1),...].
   ! The photons share the available energy eav=q+e*awr/(awr+1); the
   ! largest kick (one photon with eav) is ekmax=eav**2/(2*(awr+1)*mc2).
   ! The mean kick (CM energy, or lab energy minus e/(awr+1)) must lie
   ! in [ekmax/20, 1.1*ekmax] (1/20 allows cascades of about 20 photons
   ! sharing eav), and the largest secondary energy with nonzero
   ! probability must not exceed 1.1*ekmax (CM) or
   ! 1.1*(sqrt(e/(awr+1))+sqrt(ekmax))**2 (lab).  Frame from LCT (LCT=3:
   ! lab for A>4).  iok=1 passes, iok=0 fails. (DLA, IAEA, Oct 2026)
   !-------------------------------------------------------------------
   ! externals
   integer::lep,lct6,iok
   real(kr)::a(*),q,zap,ekick,ekmax,emx,elim
   ! internals
   integer::nd,na,nep,ncyc,k,ib
   logical::lab
   real(kr)::e,en,eav,x0,x1,f0,f1,dx,sn,sd
   real(kr),parameter::zero=0

   e=a(2)
   nd=nint(a(3))
   na=nint(a(4))
   nep=nint(a(6))
   ncyc=na+2
   en=e/(awr+1)
   eav=q+e*awr/(awr+1)
   if (eav.lt.zero) eav=0
   ekmax=eav*eav*rtm/2
   lab=(lct6.eq.1).or.(lct6.eq.3.and.mod(nint(zap),1000).gt.4)
   sn=0
   sd=0
   emx=0
   do k=1,nd
      ib=6+ncyc*(k-1)
      x0=a(ib+1)
      f0=a(ib+2)
      sn=sn+x0*f0
      sd=sd+f0
      if (f0.gt.zero.and.x0.gt.emx) emx=x0
   enddo
   do k=nd+1,nep-1
      ib=6+ncyc*(k-1)
      x0=a(ib+1)
      f0=a(ib+2)
      x1=a(ib+ncyc+1)
      f1=a(ib+ncyc+2)
      dx=x1-x0
      if (lep.eq.1) then
         sn=sn+f0*dx*(x0+x1)/2
         sd=sd+f0*dx
         if (f0.gt.zero.and.x1.gt.emx) emx=x1
      else
         sn=sn+dx*(f0*(2*x0+x1)+f1*(x0+2*x1))/6
         sd=sd+dx*(f0+f1)/2
         if ((f0.gt.zero.or.f1.gt.zero).and.x1.gt.emx) emx=x1
      endif
   enddo
   ekick=0
   if (sd.gt.zero) ekick=sn/sd
   elim=1.1e0_kr*ekmax
   if (lab) then
      ekick=ekick-en
      elim=1.1e0_kr*(sqrt(en)+sqrt(ekmax))**2
   endif
   iok=1
   if (sd.le.zero) iok=0
   if (ekick.lt.ekmax/20.or.ekick.gt.1.1e0_kr*ekmax) iok=0
   if (emx.gt.elim) iok=0
   return
   end subroutine rec102

   subroutine nheat(iold,inew,nscr,nend4,nend6,local)
   !-------------------------------------------------------------------
   ! Calculate contributions to heating from neutron files.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use util   ! provides error,mess,loada,finda,skiprz,sigfig
   use endf   ! provides endf routines and variables
   ! externals
   integer::iold,inew,nscr,nend4,nend6,local
   ! internals
   integer::nb,nw,npkk,ipk,n6,j6,irec,jrec,last6,new6
   integer::lr,icon,mtd,i,iimt,nmt,na,nwa,nwm,idis
   integer::ipx,irx,ilist,nn,iii,lq0,idx,ii,indxx
   integer::ipfr,irfr,ipnp,irnp,ipgp,irgp
   integer::npkt,intt,isave,npkkk,nwmax,iflag
   integer::i5add
   integer::nnr,ncr,zcr(3),acr(3)
   logical::lf5,lfis
   real(kr)::dfis
   real(kr)::awfac,aw1fac,e,thresh,y,t,qendf,pnue
   real(kr)::qs,q0,h0,ebar,dame,yld,qsave
   real(kr)::elst,yld0,xxx,test,enext,enx,ebal6,h
   real(kr)::ebarp,yldp,yp,hp,damep,ebal6p,hmin,hmax,tt
   real(kr)::qfr,qnp,qgp
   real(kr)::c(50)
   integer::imt(50)
   character(60)::strng,strng1,strng2
   real(kr),dimension(:),allocatable::b
   real(kr),dimension(:),allocatable::a
   real(kr),dimension(:),allocatable::d
   real(kr),dimension(:),allocatable::scr
   real(kr),dimension(:),allocatable::bufo,bufn
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::efis=15.e6_kr
   real(kr),parameter::fq1=8.07e6_kr
   real(kr),parameter::fq2=1.307e0_kr
   real(kr),parameter::one=1
   real(kr),parameter::ten=10
   real(kr),parameter::qtest=99.e6_kr
   real(kr),parameter::smin=1.e-9_kr
   real(kr),parameter::break1=1.e6_kr
   real(kr),parameter::break2=2.e6_kr
   real(kr),parameter::break3=5.e6_kr
   real(kr),parameter::break4=20.e6_kr
   real(kr),parameter::step1=.2e6_kr
   real(kr),parameter::step2=.5e6_kr
   real(kr),parameter::step3=1.e6_kr
   real(kr),parameter::step4=5.e6_kr
   real(kr),parameter::up=1.1e0_kr
   real(kr),parameter::zero=0
   real(kr),parameter::tol=1.e-5_kr

   !--allocate storage
   allocate(b(npage+50))
   nwmax=250000
   allocate(d(nwmax))
   na=500000
   allocate(a(na))
   allocate(scr(npage+50))
   allocate(bufo(nbuf))
   allocate(bufn(nbuf))
   awfac=1/awr
   aw1fac=1/(awr+1)
   npkk=npk
   if (kchk.eq.1) npkk=3*npk-2
   if (mgam.gt.0) npkk=npkk+3
   if (mgam.eq.0.and.i6.gt.0) npkk=npkk+1
   npkkk=0
   do ipk=3,npk
      if (mtp(ipk).eq.442) npkkk=ipk
   enddo
   ! check for requested dame energy calculation
   idame=0
   ipk=2
   do while (idame.eq.0.and.ipk.lt.npk)
      ipk=ipk+1
      if (mtp(ipk).ge.444) idame=1
   enddo
   n6=0
   j6=0
   i6g=0
   i5add=0
   n5t=0
   irec=0
   jrec=0
   last6=0
   new6=0
   qsave=0
   pnue=0

   !--loop over non-redundant mt-s in file 3.
  105 continue
   call contio(nin,0,0,scr,nb,nw)
   if (mfh.eq.0) go to 300
   if (mth.eq.0) go to 105
   if (mth.eq.3.or.mth.eq.4) go to 110
   if (mth.eq.10) go to 110
   if (mth.eq.26.or.mth.eq.27) go to 110
   if (mth.eq.18.and.mt19.eq.1) go to 110
   if (mth.ge.19.and.mth.le.21.and.mt19.eq.2) go to 110
   if (mth.eq.38.and.mt19.eq.2) go to 110
   if (mth.eq.101) go to 110
   if (mth.eq.103.and.mt103.gt.0) go to 110
   if (mth.eq.104.and.mt104.gt.0) go to 110
   if (mth.eq.105.and.mt105.gt.0) go to 110
   if (mth.eq.106.and.mt106.gt.0) go to 110
   if (mth.eq.107.and.mt107.gt.0) go to 110
   if (mth.eq.16.and.mt16.gt.0) go to 110
   if (iverf.le.5.and.mth.eq.719) go to 110
   if (iverf.le.5.and.mth.eq.739) go to 110
   if (iverf.le.5.and.mth.eq.759) go to 110
   if (iverf.le.5.and.mth.eq.779) go to 110
   if (iverf.le.5.and.mth.eq.799) go to 110
   if (mth.gt.120.and.mth.lt.152) go to 110
   if (mth.gt.200.and.mth.lt.600) go to 110
   zat=c1h
   awrt=c2h
   izat=nint(zat)
   go to 120
  110 continue
   call tosend(nin,0,0,scr)
   go to 105

   !--compute heating from this section.
  120 continue
   e=0
   call gety1(e,thresh,idis,y,nin,b)
   thresh=sigfig(thresh,7,0)
   t=c1h
   q=c2h
   qendf=q
   lr=l2h
   qs=0
   q0=0
   h0=0
   icon=0
   ebar=0
   dame=0
   mtd=mth

   !--check for an endf-6 file 6
   if (j6.gt.0) go to 123
   do i=1,i6
      if (mtd.eq.mt6(i)) go to 121
   enddo
   go to 124
  121 continue
   call findf(matd,6,mtd,nend6)
   call contio(nend6,0,0,scr,nb,nw)
   n6=nint(scr(5))
   lct=nint(scr(4))
   jp=nint(scr(3))
   jpn=mod(jp,10)
   jpp=jp-10*jpn
   if (mth.eq.18 .and. jpn.ne.0) then
      call tosend(nend6,0,0,scr)
      go to 105
   endif
   ! fission given in File 6: damage table from the
   ! neutron subsection, then position again at the section start
   ! (DLA, IAEA, Oct 2026)
   nf6=0
   if (idame.gt.0.and.((mtd.ge.18.and.mtd.le.21).or.mtd.eq.38)) then
      call fis6tab(nend6,matd,mtd,za,awr)
      call findf(matd,6,mtd,nend6)
      call contio(nend6,0,0,scr,nb,nw)
   endif
   i6g=i6g+1
   ! one-particle recoil correction for multi-particle emission
   nrt=0
   if (idame.gt.0.and.mt6no(i6g).gt.0.and.mtd.ne.5.and.mtd.ne.102.and.&
     .not.((mtd.ge.18.and.mtd.le.21).or.mtd.eq.38))&
     call oprtab(nend6,matd,mtd,scr(1),scr(2),lct,n6,mt6no(i6g))
   ! mt5 residuals without distribution (LAW=0) or not
   ! given: recoil from the momentum balance of the light particles,
   ! added once, in the first subsection processed (DLA, IAEA, Oct 2026)
   n5t=0
   i5add=0
   if (mtd.eq.5) then
      call mt5tab(nend6,matd,mtd,za,awr,lct,n6,thresh)
      if (n5t.gt.0) i5add=1
   endif
  123 continue
   j6=j6+1
   if (irec.gt.0.and.mt6no(i6g).gt.0) j6=j6-1
   last6=0
   if (j6.eq.n6.and.irec.gt.0) last6=1
   if (j6.eq.n6.and.mt6no(i6g).eq.0) last6=1
   new6=0
   if (j6.eq.1.and.irec.eq.0) new6=1
  124 continue

   !--adjust fission q to e(in)=0 prompt value
   if ((mth.ge.18.and.mth.le.21).or.mth.eq.38) then
      if (mt458.eq.1) then
         qendf=c458(15)
         q=qendf-qdel
         write(strng1,&
              &'(''changed q from '',1p,e14.6,'' to '',1p,e14.6)')&
              qendf,q
         write(strng2,'(''for mt '',i3,'' by taking out delayed components'')') mth
         call mess('nheat',strng1,strng2)
      else
         q=qendf
         cpoly(0)=q
      endif
   endif

   !--choose appropriate ground state q value.
   !--initialize heating routines.
   !--define neutron yield for this reaction.
   iimt=0
   if (nqa.eq.0) go to 140
   do 130 i=1,nqa
      if (mta(i).ne.mth) go to 130
      qs=qa(i)
      iimt=i
      go to 140
  130 continue
  140 continue
   if (j6.gt.0) go to 175
   ! mt5 is described only in File 6 (ENDF-6); without
   ! File 6 NJOY2016.80 treated it as a discrete level (disbar)
   ! (DLA, IAEA, Oct 2026)
   if (mth.eq.5.and.iverf.ge.6) call error('nheat',&
     'mt5 without File 6.',&
     'ENDF-6: mt5 products must be given in File 6.')
   ! continuum reactions without File 6:
   !  - with File 5: conbar, also for mt11 (n,2nd) and for the
   !    multi-particle and multi-neutron reactions mt152-200 and
   !    mt875-891 ((n,5n), ..., (n,2n) to excited states), which the
   !    NJOY2016.80 sent to capdam as (n,p);
   !  - mt16-50 without File 5 (except fission): capdam -> capmul
   !    (NJOY2016.80 stopped in conbar: File 5 not found);
   !  - mt11 without File 5 stays with disbar as before. (DLA, IAEA, Oct 2026)
   call mulpart(mth,nnr,ncr,zcr,acr)
   lf5=.false.
   do i=1,n5mt
      if (mt5l(i).eq.mth) lf5=.true.
   enddo
   lfis=((mth.ge.18.and.mth.le.21).or.mth.eq.38)
   if (lf5.and.(mth.eq.11.or.(mth.ge.152.and.mth.le.200).or.&
     (mth.ge.875.and.mth.le.891))) go to 150
   if (mth.gt.15.and.mth.lt.51) then
      if (.not.lf5.and..not.lfis.and.(ncr.gt.0.or.nnr.ge.2)) go to 170
      go to 150
   endif
   if (mth.gt.100) go to 170
   icon=1
   q0=0
   if (lr.ne.0.and.lr.ne.31) q0=t
   if (iimt.ne.0) q0=qs
   if (mth.eq.91) go to 160
   yld=1
   ! accomodate neutron producing pseudo-levels
   if (lr.eq.16.or.lr.eq.21.or.lr.eq.24.or.lr.eq.26.or.lr.eq.30) yld=2
   if (lr.eq.17.or.lr.eq.25.or.lr.eq.38) yld=3
   if (lr.eq.37) yld=4
   nwa=na
   call disbar(e,ebar,dame,q,za,awr,nend4,matd,mtd,a,nwa)
   go to 180
  150 continue
   q0=q
   if (iimt.ne.0) q0=qs
  160 continue
   icon=2
   nwm=nwmax
   nwa=na
   call conbar(e,ebar,dame,nendf,nscr,matd,mtd,a,nwa,d,nwm,yld,q)
   go to 180
  170 continue
   if (iverf.le.5.and.mth.ge.700.and.mth.le.799) go to 171
   if (iverf.ge.6.and.mth.ge.600.and.mth.le.849) go to 171
   q0=q
   if (idame.gt.0) call capdam(e,dame,q,za,awr,mtd)
   go to 180
  171 continue
   if (iverf.le.5.and.mth.eq.718) go to 174
   if (iverf.le.5.and.mth.eq.738) go to 174
   if (iverf.le.5.and.mth.eq.758) go to 174
   if (iverf.le.5.and.mth.eq.778) go to 174
   if (iverf.le.5.and.mth.eq.798) go to 174
   if (iverf.ge.6.and.mth.eq.649) go to 174
   if (iverf.ge.6.and.mth.eq.699) go to 174
   if (iverf.ge.6.and.mth.eq.749) go to 174
   if (iverf.ge.6.and.mth.eq.799) go to 174
   if (iverf.ge.6.and.mth.eq.849) go to 174
   ! discrete charged-particle levels
   icon=1
   if (iverf.le.5.and.mth.eq.700) qsave=q
   if (iverf.le.5.and.mth.eq.720) qsave=q
   if (iverf.le.5.and.mth.eq.740) qsave=q
   if (iverf.le.5.and.mth.eq.760) qsave=q
   if (iverf.le.5.and.mth.eq.780) qsave=q
   if (iverf.le.5) q0=qsave
   if (iverf.ge.6) q0=t
   yld=1
   nwa=na
   call disbar(e,ebar,dame,q,za,awr,nend4,matd,mtd,a,nwa)
   go to 180
   ! continuum charged-particle emission
  174 continue
   yld=1
   if (iverf.le.5) q0=qsave
   if (iverf.ge.6) q0=t
   if (idame.gt.0) call capdam(e,dame,q,za,awr,mtd)
   go to 180
   ! file 6 reactions
  175 continue
   nwm=nwmax
   nwa=na
   call sixbar(e,ebar,yld,dame,nend6,a,nwa,nscr,d,nwm,n6,j6,irec,jrec,iflag)
   if (ebar.lt.zero) go to 291
   icon=-1
   q0=t
   if (iimt.ne.0) q0=qs
   if ((mtd.ge.18.and.mtd.le.21).or.mtd.eq.38) q0=q
  180 continue
   nmt=0
   if (npk.gt.2) call indx(npk,mtp,mt19,mtd,imt,nmt)
   if (icon.lt.0) go to 179
   if (iprint.ne.1) go to 182
   if (qs.ge.qtest) go to 181
   if ((mtd.ge.18.and.mtd.le.21).or.mtd.eq.38) go to 183
   if (icon.lt.0) go to 179
   if (idame.eq.0) write(nsyso,'(/&
     &'' neutron heating for mt'',i3,''   q0 ='',1p,e12.4,5x,&
     &''q ='',e12.4/14x,''e'',10x,''ebar'',9x,''yield'',10x,&
     &''xsec'',7x,''heating'')') mtd,q0,qendf
   if (idame.gt.0) write(nsyso,'(/&
     &'' neutron heating for mt'',i3,''   q0 ='',1p,e12.4,5x,&
     &''q ='',e12.4/14x,''e'',10x,''ebar'',9x,''yield'',10x,&
     &''xsec'',7x,''heating'',8x,''damage'')') mtd,q0,qendf
   go to 182
  183 continue
   if (idame.eq.0) write(nsyso,'(/&
     &'' neutron heating for mt'',i3,&
     &''   q0= variable fission q''/&
     &14x,''e'',12x,''q0'',10x,''ebar'',9x,''yield'',10x,&
     &''xsec'',7x,''heating'')') mtd
   if (idame.gt.0) write(nsyso,'(/&
     &'' neutron heating for mt'',i3,&
     &''   q0= variable fission q''/&
     &14x,''e'',12x,''q0'',10x,''ebar'',9x,''yield'',10x,&
     &''xsec'',7x,''heating'',8x,''damage'')') mtd
   if (irec.gt.0.and.icon.lt.0) izap=100
   go to 182
  181 continue
   if (idame.eq.0) write(nsyso,'(/&
     &'' neutron heating for mt'',i3,''   q0= variable'',5x,&
     &''q ='',e12.4/14x,''e'',10x,''qbar'',10x,''ebar'',9x,&
     &''yield'',10x,''xsec'',7x,''heating'')') mtd,qendf
   if (idame.gt.0) write(nsyso,'(/&
     &'' neutron heating for mt'',i3,''   q0= variable'',5x,&
     &''q ='',e12.4/14x,''e'',10x,''qbar'',10x,''ebar'',9x,&
     &''yield'',10x,''xsec'',7x,''heating'',8x,''damage'')')&
     &mtd,qendf
   go to 182
  179 continue
   izap=nint(zap)
   if (irec.gt.0) then
      if (mod(izat,1000).eq.0) then
         izap=1000*(izat/1000-izap/1000)
      else
         izap=nint(zat+1-zap)
      endif
   endif
   if (iprint.ne.1) go to 182
   if (qs.ge.qtest) go to 178
   if (idame.eq.0) write(nsyso,'(/&
     &'' file six heating for mt'',i3,'', particle ='',i6,&
     &5x,''q = '',1p,e12.4/14x,''e'',10x,''ebar'',9x,&
     &''yield'',10x,''xsec'',7x,''heating'')') mtd,izap,q0
   if (idame.gt.0) write(nsyso,'(/&
     &'' file six heating for mt'',i3,'', particle ='',i6,&
     &5x,''q = '',1p,e12.4/14x,''e'',10x,''ebar'',9x,&
     &''yield'',10x,''xsec'',7x,''heating'',8x,''damage'')')&
     mtd,izap,q0
   go to 182
  178 continue
   if (idame.eq.0) write(nsyso,'(/&
     &'' file six heating for mt'',i3,'', particle ='',i6,&
     &5x,''q0 = variable''/&
     &14x,''e'',10x,''qbar'',10x,''ebar'',9x,''yield'',10x,&
     &''xsec'',7x,''heating'')') mtd,izap
   if (idame.gt.0) write(nsyso,'(/&
     &'' file six heating for mt'',i3,'', particle ='',i6,&
     &5x,''q0 = variable''/&
     &14x,''e'',10x,''qbar'',10x,''ebar'',9x,''yield'',10x,&
     &''xsec'',7x,''heating'',8x,''damage'')') mtd,izap
  182 continue
   ipx=2
   irx=1
   ipfr=2
   irfr=1
   ipnp=2
   irnp=1
   ipgp=2
   irgp=1
   elst=0
   ilist=1
   if (icon.eq.0) yld=0
   yld0=0

   !--compute heating for energies in union grid.
   nn=0
  190 continue
   nn=nn+1
   call finda(nn,c,npkk,iold,bufo,nbuf)
   e=c(1)
   e=sigfig(e,9,0)
   if (new6.eq.1) c(npkk)=0
   if (mtd.ne.2) go to 193
   if (e.lt.elst.and.nn.lt.ne) go to 193
   if (ilist.gt.ilmax) call error('nheat','storage exceeded.',' ')
   elist(ilist)=e
   elist(ilist+1)=big
   if (e.lt.break1) then
      iii=int(log10(e*(1+small))+small)
      if (log10(e*(1+small)).lt.zero) iii=iii-1
      xxx=ten**iii
      if (2*xxx.gt.e*(1+small).and.e.ge.one*(1-small)) then
         elst=2*xxx
      else if (5*xxx.gt.e*(1+small).and.e.ge.one*(1-small)) then
         elst=5*xxx
      else
         elst=10*xxx
      endif
   else
      test=break4
      test=sigfig(test,9,1)
      if (e.ge.test) elst=step4*int(e/step4+up)
      if (e.lt.test) elst=step3*int(e/step3+up)
      test=break3
      test=sigfig(test,9,1)
      if (e.lt.test) elst=step2*int(e/step2+up)
      test=break2
      test=sigfig(test,9,1)
      if (e.lt.test) elst=step1*int(e/step1+up)
   endif
   elst=elst-elst/10000000
  193 continue
   if (e.lt.thresh) go to 290
   call gety1(e,enext,idis,y,nin,b) ! y is the reaction cross section
   ! check for energy-dependent q
   if (iimt.gt.0) then
      if (qa(iimt).ge.qtest) then
         lq0=lqs(iimt)
         call terpa(q0,e,enx,idx,qbar(lq0),ipx,irx)
         q=q0
      endif
   endif
   if (icon.eq.1) call disbar(e,ebar,dame,q,za,awr,&
     nend4,matd,mtd,a,nwa)
   if (icon.eq.2) call conbar(e,ebar,dame,nendf,nscr,matd,mtd,&
     a,nwa,d,nwm,yld,q)
   if (icon.eq.0.and.idame.gt.0) call capdam(e,dame,q,za,awr,mtd)
   if (icon.lt.0)&
     call sixbar(e,ebar,yld,dame,nend6,a,nwa,nscr,d,nwm,n6,j6,irec,jrec,iflag)
   ! one-particle recoil: correction for the other emitted particles
   if (icon.lt.0.and.irec.gt.0.and.nrt.gt.0) dame=dame*oprfac(e)
   if (yld0.eq.zero) yld0=yld
   ! correct fission q for current incident energy
   if ((mtd.ge.18.and.mtd.le.21).or.mtd.eq.38) then
      if (lfc.eq.0) then     ! thermal point or polynomial
         if (nply.eq.0) then
            q0=cpoly(0)-(yld-yld0)*fq1+fq2*e
            if (mt458.eq.1) then
               h0=hpoly(0)
            else
               h0=q0-ebar*yld
            endif
         else
            q0=cpoly(nply)
            h0=hpoly(nply)
            do i=nply-1,0,-1
               q0=q0*e+cpoly(i)
               h0=h0*e+hpoly(i)
            enddo
         endif
         pnue=q0-h0
      else                   ! tabulated
         if (etabmax.lt.e.and.((abs(etabmax-e)/etabmax).gt.tol)) then
            write(strng,'(''upper energy tabulated fission q components'',&
               &'' is too low.'')')
            call error('nheat',strng,' ')
         endif
         if (ifc1.eq.1) then
            call terpa(qfr,e,enx,idx,afr,ipfr,irfr)
         else
            qfr=c458(1)
         endif
         if (ifc2.eq.1) then
            call terpa(qnp,e,enx,idx,anp,ipnp,irnp)
         else
            qnp=c458(3)-(yld-yld0)*fq1+fq2*e
         endif
         if (ifc4.eq.1) then
            call terpa(qgp,e,enx,idx,agp,ipgp,irgp)
         else
            qgp=c458(7)
         endif
         pnue=qnp
         h0=qfr+qgp
         q0=h0+pnue
      endif
      q0=q0-e
   endif
   if (icon.ge.0) then
      ebal6=0
      dame=dame*y
      if ((mtd.lt.18.or.mtd.gt.21).and.mtd.ne.38) then
         h=(e+q0-ebar*yld)*y
      else
         h=h0*y
      endif
   else
      if ((mtd.lt.18.or.mtd.gt.21).and.mtd.ne.38) then
         h=ebar*yld*y
      else
         h=h0*y
      endif
      dame=dame*y
      c(npkk)=c(npkk)+h
      if (izap.eq.0) c(npkk-1)=c(npkk-1)+h
      if (izap.eq.0.and.npkkk.gt.0) c(npkkk)=c(npkkk)+h
      if (izap.eq.1) h=0
      if (izap.eq.0.and.local.eq.0) h=0
      if (izap.le.1) dame=0
      ! recoil of the mt5 residuals without distribution (DLA, IAEA, Oct 2026)
      if (i5add.eq.1.and.irec.eq.0) then
         h=h+mt5hm(e)*y
         dame=dame+mt5dm(e)*y
      endif
      ebal6=0
      if (last6.eq.1.and.mth.ne.5.and.mth.ne.102)&
        ebal6=(e+q0)*y-c(npkk)
      if ((mtd.ge.18.and.mtd.le.21).or.mtd.eq.38) ebal6=0
   endif
   if (mtd.ge.46.and.mtd.le.49) then
      ! treat second neutron from sequential (n,2n), mt=46-49.
      if (qs.ge.zero) call error('nheat',&
        'binding energy for sequential n2n needed.',' ')
      h=(qs-ebar)*y
   endif
   ! fission: damage per fission dfis, the same for File 5
   ! and File 6; with File 6 the heating h0 and the damage are counted
   ! once, in the first subsection (NJOY2016.80 added h0 once for
   ! every non-neutron subsection, e.g. 4 times with the ternary p, t,
   ! alpha and the one-particle recoil, and built a recoil from the last
   ! particle of the section). (DLA, IAEA, Oct 2026)
   if ((mtd.ge.18.and.mtd.le.21).or.mtd.eq.38) then
      if (icon.lt.0) then
         ebal6=0
         h=0
         dfis=0
         if (j6.eq.1.and.irec.eq.0) then
            h=h0*y
            if (idame.gt.0) dfis=fis6dm(e)
         endif
         dame=dfis*y
      endif
   endif

   !--accumulate total and partial heating factors.
   c(2)=c(2)+h+ebal6
   if ((mtd.ge.18.and.mtd.le.21).or.mtd.eq.38) ebar=pnue/yld
   ebarp=sigfig(ebar,9,0)
   yldp=sigfig(yld,9,0)
   yp=sigfig(y,9,0)
   hp=sigfig(h,9,0)
   damep=sigfig(dame,9,0)
   ebal6p=sigfig(ebal6,9,0)
   if (iprint.ne.1.or.abs(e-elist(ilist)).gt.small*e&
     .or.y.le.smin) go to 199
   if (qs.ge.qtest) go to 198
   if ((mtd.ge.18.and.mtd.le.21).or.mtd.eq.38) go to 198
   if (idame.eq.0) then
      write(nsyso,'(1x,1p,7e14.4)') e,ebarp,yldp,yp,hp
   else
      write(nsyso,'(1x,1p,7e14.4)') e,ebarp,yldp,yp,hp,damep
   endif
   go to 199
  198 continue
   if (idame.eq.0) then
      write(nsyso,'(1x,1p,7e14.4)') e,q0,ebarp,yldp,yp,hp
   else
      write(nsyso,'(1x,1p,7e14.4)') e,q0,ebarp,yldp,yp,hp,damep
   endif
  199 continue
   if (abs(ebal6).ge.100*y) then
      if (iprint.eq.1.and.abs(e-elist(ilist)).le.small*e&
        .and.y.gt.smin) then
         if (qs.lt.qtest) then
            write(nsyso,'(53x,''ebal'',1p,e14.4)') ebal6p
         else
            write(nsyso,'(67x,''ebal'',1p,e14.4)') ebal6p
         endif
      endif
   endif
   if (nmt.ne.0) then
      ! partial kermas
      do ii=1,nmt
         indxx=imt(ii)
         if (mtp(indxx).ge.444) c(indxx)=c(indxx)+dame
         if (mtp(indxx).lt.442) c(indxx)=c(indxx)+h+ebal6
      enddo
   endif

   !--compute kinematic limits on heating, if desired.
   if (kchk.eq.0) go to 290
   if (icon.lt.0) go to 275
   if (mtd.ne.2) go to 215
   hmin=h
   hmax=hmin
   go to 285
  215 continue
   if (mtd.ne.102) go to 220
   hmin=y*e*aw1fac
   tt=q+awr*e*aw1fac
   hmax=hmin+y*tt*tt*rtm/2
   go to 285
  220 continue
   if (mtd.ge.51.and.mtd.lt.91) go to 221
   if (iverf.le.5.and.mtd.ge.700.and.mtd.lt.718) go to 221
   if (iverf.le.5.and.mtd.ge.720.and.mtd.lt.738) go to 221
   if (iverf.le.5.and.mtd.ge.740.and.mtd.lt.758) go to 221
   if (iverf.le.5.and.mtd.ge.760.and.mtd.lt.778) go to 221
   if (iverf.le.5.and.mtd.ge.780.and.mtd.lt.798) go to 221
   if (iverf.ge.6.and.mtd.ge.600.and.mtd.lt.649) go to 221
   if (iverf.ge.6.and.mtd.ge.650.and.mtd.lt.699) go to 221
   if (iverf.ge.6.and.mtd.ge.700.and.mtd.lt.749) go to 221
   if (iverf.ge.6.and.mtd.ge.750.and.mtd.lt.799) go to 221
   if (iverf.ge.6.and.mtd.ge.800.and.mtd.lt.849) go to 221
   go to 230
  221 continue
   hmin=(e+q-ebar*yld)*y
   hmin=(e+q-ebar*yld)*y
   hmax=hmin
   if (lr.ne.0.and.lr.ne.31) hmax=(e+q0-ebar*yld)*y
   go to 285
  230 continue
   if (mtd.ne.91) go to 240
   hmin=y*(e+ebar)*awfac
   hmax=hmin
   if (lr.ne.0.and.lr.ne.31) hmax=(e+q0-ebar)*y
   go to 285
  240 continue
   if (mtd.lt.103) go to 250
   hmin=y*e*awfac
   hmax=h
   go to 285
  250 continue
   if (mtd.ne.18.and.mtd.ne.20.and.mtd.ne.21.and.mtd.ne.38) go to 260
   hmin=y*(e+q0-ebar*yld/2-efis)
   hmax=h
   go to 285
  260 continue
   test=10
   if (mtd.ne.16.or.awr.lt.test) go to 265
   hmin=0
   hmax=y*(e+ebar)/(awr-1)
   if (hmax.gt.h) hmax=h
   go to 285
  265 continue
   test=10
   if (mtd.ne.17.or.awr.lt.test) go to 270
   hmin=0
   hmax=y*(e+2*ebar)/(awr-2)
   if (hmax.gt.h) hmax=h
   go to 285
  270 continue
   hmin=0
   hmax=h
   go to 285
  275 continue
   hmin=h
   hmax=h
  285 continue
   npkt=npk-1
   if (kchk.ne.1) go to 283
   c(2+npkt)=c(2+npkt)+hmin
   c(2+2*npkt)=c(2+2*npkt)+hmax
  283 continue
   if (nmt.eq.0) go to 290
   do 286 ii=1,nmt
   indxx=imt(ii)
   if (mtp(indxx).eq.442) go to 286
   if (mtp(indxx).eq.443) go to 284
   if (mtp(indxx).ge.444) go to 286
   if (kchk.ne.1) go to 286
   c(indxx+npkt)=c(indxx+npkt)+hmin
   c(indxx+2*npkt)=c(indxx+2*npkt)+hmax
   go to 286
  284 continue
   c(indxx)=c(indxx)+hmax
  286 continue

   !--store results on scratch and go to next energy.
  290 continue
   intt=nn
   if (intt.eq.ne) intt=-intt
   call loada(intt,c,npkk,inew,bufn,nbuf)
   if (e.eq.elist(ilist)) ilist=ilist+1
   if (nn.lt.ne) go to 190
   if (i5add.eq.1.and.irec.eq.0) i5add=0

   !--for file 6, loop over subsections
   isave=iold
   iold=inew
   inew=isave
  291 continue
   if (j6.ge.n6) then
      ! make sure we're at the end of section (will not be so if
      ! there was no particle distribution for the last subsection
      ! in this section).
      if (ebar.lt.zero) then
         call tosend(nin,0,0,b)
         call skiprz(nin,-1)
      endif
      go to 296
   endif
   call skiprz(nin,-2)
   call findf(matd,3,mtd,nin)
   go to 297
  296 continue
   if (j6.eq.0) go to 297
   if (irec.gt.0.and.mt6no(i6g).gt.0) go to 293
   if (iprint.eq.1.and.i6g.gt.0.and.izap.ne.0.and.iflag.eq.0)&
     write(nsyso,'(/''   no explicit file 6 photon production for mt'',i3)')&
     mtd
   if (mt6no(i6g).ne.0) then
      if (mth.ne.102 ) write(nsyso,'(/&
        &''   generating recoil with one-particle approx.'')')
      if (mth.eq.102 ) write(nsyso,'(/&
        &''   generating recoil from photon momentum.'')')
      irec=mt6no(i6g)
      call skiprz(nin,-2)
      call findf(matd,3,mtd,nin)
      go to 297
   endif
  293 continue
   irec=0
   jrec=0
   j6=0
   n6=0
  297 continue
   go to 105

   !--nheat is finished.
  300 continue
   deallocate(b)
   deallocate(d)
   deallocate(a)
   deallocate(scr)
   deallocate(bufn)
   deallocate(bufo)
   call skiprz(nin,-4)
   return
   end subroutine nheat

   subroutine indx(npk,mtp,mt19,mt,imt,nmt)
   !-------------------------------------------------------------------
   ! Find the index of a given MT in the stored desired MT array.
   ! Heating identifiers
   !     300 + MT for reaction
   ! Special values allowed are
   !     303=nonelastic (all but MT2)
   !     304=inelastic (MTs1-91)
   !     318=fission (MT18 or MT19-21 and MT38)
   !     401=disappearance (MT102-120)
   !     442=total photon ev-barns in kerma
   !     443=kinematic kerma
   ! Damage energy identifiers
   !     444=total
   !     445=elastic
   !     446=inelastic (MT51-MT91)
   !     447=total capture (MT102-MT120)
   !-------------------------------------------------------------------
   use endf ! provides iverf
   ! externals
   integer::npk,mtp(*),mt19,mt,imt(*),nmt
   ! internals
   integer::i,mtpi,iflag

   nmt=0
   if (npk.lt.3) return
   do i=3,npk
      mtpi=mtp(i)
      iflag=0
      if (mtpi.eq.442) iflag=1
      if (mtpi.eq.443) iflag=1
      if (mtpi.eq.444) iflag=1
      if (mtpi.eq.445.and.mt.eq.2) iflag=1
      if (mtpi.eq.446.and.mt.ge.51.and.mt.le.91) iflag=1
      if (mtpi.eq.447.and.mt.ge.102.and.mt.le.120) iflag=1
      if (iverf.lt.6) then
         if (mtpi.eq.447.and.mt.ge.700.and.mt.le.799) iflag=1
      else
         if (mtpi.eq.447.and.mt.ge.600.and.mt.le.849) iflag=1
      endif
      if (mtpi.eq.303.and.mt.ne.2.and.mt.ne.3) iflag=1
      if (mtpi.eq.304.and.mt.ge.51.and.mt.le.91) iflag=1
      if (mtpi.eq.318) then
         if (mt19.eq.0.and.mt.eq.18) iflag=1
         if (mt19.eq.1.and.(mt.eq.19.or.&
           mt.eq.20.or.mt.eq.21.or.mt.eq.38)) iflag=1
         if (mt19.eq.2.and.mt.eq.18) iflag=1
      endif
      if (mtpi.eq.401) then
         if (mt.ge.102.and.mt.le.120) iflag=1
         if (iverf.lt.6.and.mt.ge.700.and.mt.lt.800) iflag=1
         if (iverf.ge.6.and.mt.ge.600.and.mt.lt.850) iflag=1
      endif
      if (iverf.lt.6) then
         if (mtpi.eq.403.and.mt.ge.700.and.mt.le.719) iflag=1
         if (mtpi.eq.404.and.mt.ge.720.and.mt.le.739) iflag=1
         if (mtpi.eq.405.and.mt.ge.740.and.mt.le.759) iflag=1
         if (mtpi.eq.406.and.mt.ge.760.and.mt.le.779) iflag=1
         if (mtpi.eq.407.and.mt.ge.780.and.mt.le.799) iflag=1
      else
         if (mtpi.eq.403.and.mt.ge.600.and.mt.le.649) iflag=1
         if (mtpi.eq.404.and.mt.ge.650.and.mt.le.699) iflag=1
         if (mtpi.eq.405.and.mt.ge.700.and.mt.le.749) iflag=1
         if (mtpi.eq.406.and.mt.ge.750.and.mt.le.799) iflag=1
         if (mtpi.eq.407.and.mt.ge.800.and.mt.le.849) iflag=1
      endif
      if (mtpi.eq.mt+300) iflag=1
      if (iflag.eq.1) then
         nmt=nmt+1
         imt(nmt)=i
      endif
   enddo
   return
   end subroutine indx

   subroutine capdam(ee,dame,q,za,awr,mtd)
   !-------------------------------------------------------------------
   ! Compute damage energy for capture reactions.
   ! For (n,p),(n,d),(n,t),(n,3He),(n,a) the emitted particle gets
   ! ea=min(eavail,vc), vc with target Z,A (NJOY2016 method).
   ! Corrections:
   !  - ENDF-6 continuum mt649,699,749,799,849 are mapped to 103..107
   !    (they were all treated as (n,p)).
   !  - mt108-117 (several charged particles: (n,2a), (n,pa), ...) are
   !    computed with capmul (they were all treated as (n,p)).
   !  - the damage of the emitted charged particle itself is added
   !    (lab energy ea+bc+2*sqrt(ea*bc)*mu, bc=ax*e/(awr+1)**2), as
   !    sixbar does for the particles given in File 6.
   !  - mt103-107 (not split into levels, no File 4/6)
   !    whose residual has a mass number <= 4 (He-3(n,p)t, He-3(n,d),
   !    He-4(n,d)t, Li-6(n,t)a, ...) are two-body reactions to a single
   !    final state: exact two-body kinematics with the QI of File 3,
   !    isotropic in the centre of mass, 64-point Gauss rule, Lindhard
   !    partition of the residual and damage of the emitted particle
   !    (as disbar).  The barrier limit ea=min(eavail,vc) is a model of
   !    the continuum and does not apply there.  Other reactions keep the
   !    NJOY2016 model. (DLA, IAEA, Oct 2026)
   !-------------------------------------------------------------------
   use endf   ! provides iverf,terp1
   use physics ! provides the light-particle mass ratios
   ! externals
   integer::mtd
   real(kr)::ee,dame,q,za,awr
   ! internals
   integer::iq,iz,mte,nnl,ncl,zcl(3),acl(3),ia
   logical::lmul,l2b
   real(kr)::e,zx,ax,denom,z,en,aw1fac
   real(kr)::awp2,thr2,afac2,arat2,r2,b2,g2,e2,u
   real(kr)::damn,daml,el,er,ea,ec,et,bc,ep,dum,epa
   real(kr),dimension(4),parameter::qp=(/&
      -.86114e0_kr,-.33998e0_kr,.33998e0_kr,.86114e0_kr/)
   real(kr),dimension(4),parameter::qw=(/&
     .34785e0_kr,.65215e0_kr,.65215e0_kr,.34785e0_kr/)
   real(kr),parameter::third=.333333333e0_kr
   real(kr),parameter::econ=1.029e6_kr
   real(kr),parameter::eps=1.e-10_kr
   real(kr),parameter::step=1.1e0_kr
   integer,parameter::nq=4
   real(kr),parameter::zero=0
   save en,damn,el,daml
   save zx,ax,denom,z,aw1fac,lmul,l2b,awp2,thr2,afac2,arat2

   !--initialize
   e=ee
   if (e.eq.zero) then
      ! ENDF-6 continuum levels are mapped to the parent reaction;
      ! without this mapping they were all treated as (n,p)
      mte=mtd
      if (iverf.ge.6.and.mtd.ge.600.and.mtd.le.849) mte=103+(mtd-600)/50
      ! reactions with several emitted particles (mt108-117, mt152-200,
      ! mt875-891): capmul.  Without this they were treated as (n,p).
      call mulpart(mtd,nnl,ncl,zcl,acl)
      lmul=(ncl.gt.0.or.nnl.ge.2)
      zx=0
      if (mte.gt.102) zx=1
      if (mte.eq.106.or.mte.eq.107) zx=2
      if (iverf.le.5.and.(mte.eq.778.or.mte.eq.798)) zx=2
      ax=0
      if (mte.gt.102) ax=1
      if (mte.eq.104) ax=2
      if (iverf.le.5.and.mte.eq.738) ax=2
      if (mte.eq.105.or.mte.eq.106) ax=3
      if (iverf.le.5.and.(mte.eq.758.or.mte.eq.778)) ax=3
      if (mte.eq.107) ax=4
      if (iverf.le.5.and.mte.eq.798) ax=4
      denom=1/awr**third
      if (ax.ne.zero) denom=1/(ax**third+awr**third)
      iz=nint(za/1000)
      z=iz
      ! two-body reactions to a single final state (DLA, IAEA, Oct 2026)
      l2b=.false.
      ia=mod(nint(za),1000)
      if (.not.lmul.and.mtd.ge.103.and.mtd.le.107.and.ia.gt.0) then
         if (ia+1-nint(ax).le.4) l2b=.true.
      endif
      if (l2b) then
         if (mtd.eq.103) awp2=pnratio
         if (mtd.eq.104) awp2=dnratio
         if (mtd.eq.105) awp2=tnratio
         if (mtd.eq.106) awp2=hnratio
         if (mtd.eq.107) awp2=anratio
         thr2=((awr+1)/awr)*(-q)
         afac2=awp2/((awr+1)**2)
         arat2=awp2/(awr+1-awp2)
      endif
      en=0
      damn=df(en,z-zx,awr+1-ax,z,awr)
      dame=0
      aw1fac=1/(awr+1)
      return
   endif

   !--normal entry
   !--if energy not in current range,
   !--construct new interpolation range.
   if (ee.ge.en*(1-eps)) then
      el=en
      daml=damn
      e=step*el
      if (e.lt.ee*(1-eps)) e=(1-eps)*ee
      dame=0
      if (mtd.eq.102) then
         ! capture followed by gamma emission
         er=e*aw1fac
         ea=awr*er+q
         er=er+ea*ea*rtm/2
         dame=df(er,z,awr+1,z,awr)
      else if (lmul) then
         ! emission of several charged particles
         call capmul(e,dame,q,za,awr,mtd)
      else if (l2b) then
         ! two-body reaction to a single final state:
         ! recoil e2=e*(1-2*g*u+g*g)*afac/arat, particle
         ! ep=e*(1+2*b*u+b*b)*afac, u=CM cosine of the particle
         ! (DLA, IAEA, Oct 2026)
         if (e.gt.thr2*(1+eps)) then
            r2=sqrt(1-thr2/e)
            b2=r2*sqrt(awr/arat2)
            g2=b2*arat2
            dum=df(zero,z-zx,awr+1-awp2,z,awr)
            do iq=1,64
               u=g64x(iq)
               e2=e*(1-2*g2*u+g2*g2)*afac2/arat2
               dame=dame+g64w(iq)*df(e2,z-zx,awr+1-awp2,z,awr)/2
            enddo
            dum=df(zero,zx,awp2,z,awr)
            do iq=1,64
               u=g64x(iq)
               ep=e*(1+2*b2*u+b2*b2)*afac2
               dame=dame+g64w(iq)*df(ep,zx,awp2,z,awr)/2
            enddo
            dum=df(zero,z-zx,awr+1-ax,z,awr)
         endif
      else if (ax.ne.zero) then
         ! capture followed by particle emission
         ea=q+awr*e*aw1fac
         if (ea.ge.zero) then
            et=(awr+1-ax)*e*aw1fac
            ec=econ*zx*z*denom
            if (ea.gt.ec*(1+eps)) ea=ec
            ! recoil of the residual nucleus (Lindhard constants are
            ! reset because the particle term below changes them)
            dum=df(zero,z-zx,awr+1-ax,z,awr)
            do iq=1,nq
               er=(et-2*sqrt(et*ax*ea)*qp(iq)+ax*ea)*aw1fac
               dame=dame+qw(iq)*df(er,z-zx,awr+1-ax,z,awr)/2
            enddo
            ! damage of the emitted charged particle itself, as done
            ! in sixbar for File 6 particles: lab energy
            ! ep=epa+bc+2*sqrt(epa*bc)*mu, bc=ax*e/(awr+1)**2
            ! epa=ea*(awr+1-ax)/(awr+1) is the cm energy of
            ! the particle when ea is the cm kinetic energy shared with
            ! the residual, as in the recoil term above (DLA, IAEA, Oct 2026)
            bc=ax*e*aw1fac*aw1fac
            epa=ea*(awr+1-ax)*aw1fac
            dum=df(zero,zx,ax,z,awr)
            do iq=1,nq
               ep=epa+bc+2*sqrt(epa*bc)*qp(iq)
               dame=dame+qw(iq)*df(ep,zx,ax,z,awr)/2
            enddo
            dum=df(zero,z-zx,awr+1-ax,z,awr)
         endif
      endif
      en=e
      damn=dame
   endif

   !--interpolate for damage energy
   call terp1(el,daml,en,damn,ee,dame,2)
   return
   end subroutine capdam

   real(kr) function evapmean(uu,ares)
   !-------------------------------------------------------------------
   ! Mean of the evaporation spectrum x*exp(-x/T) truncated at U=uu
   ! (energies in eV), with T=sqrt(U/a) and a=ares/8 1/MeV.  Used for
   ! the neutrons emitted in (n,nX) reactions (capmul).
   !-------------------------------------------------------------------
   real(kr)::uu,ares
   real(kr)::tt,xu,ex
   real(kr),parameter::bigu=50.e0_kr
   real(kr),parameter::zero=0
   evapmean=0
   if (uu.le.zero.or.ares.le.zero) return
   tt=1.e6_kr*sqrt(1.e-6_kr*uu/(ares/8))
   xu=uu/tt
   if (xu.lt.1.e-3_kr) then
      evapmean=2*uu/3
   else if (xu.gt.bigu) then
      evapmean=2*tt
   else
      ex=exp(-xu)
      evapmean=tt*(2-ex*(xu*xu+2*xu+2))/(1-ex*(1+xu))
   endif
   return
   end function evapmean

   subroutine mulpart(mtd,nn,nc,zc,ac)
   !-------------------------------------------------------------------
   ! Number of neutrons nn and list of charged particles (zc,ac) for
   ! reactions with emission of several particles: (n,np), (n,na),
   ! ..., (n,2a), (n,pa), ..., the ENDF-6 reactions mt152-200 and the
   ! neutron-only reactions (n,2n), (n,3n), (n,4n), (n,xn), mt875-891.
   ! nn=nc=0 for reactions not in the list (single-particle emission,
   ! fission, capture, ...).
   !-------------------------------------------------------------------
   integer::mtd,nn,nc,zc(3),ac(3)
   nn=0
   nc=0
   zc=0
   ac=0
   select case (mtd)
   case (11)   ! (n,2nd)
      nn=2; nc=1; zc(1)=1; ac(1)=2
   case (22)   ! (n,na)
      nn=1; nc=1; zc(1)=2; ac(1)=4
   case (23)   ! (n,n3a)
      nn=1; nc=3; zc(1:3)=2; ac(1:3)=4
   case (24)   ! (n,2na)
      nn=2; nc=1; zc(1)=2; ac(1)=4
   case (25)   ! (n,3na)
      nn=3; nc=1; zc(1)=2; ac(1)=4
   case (28)   ! (n,np)
      nn=1; nc=1; zc(1)=1; ac(1)=1
   case (29)   ! (n,n2a)
      nn=1; nc=2; zc(1:2)=2; ac(1:2)=4
   case (30)   ! (n,2n2a)
      nn=2; nc=2; zc(1:2)=2; ac(1:2)=4
   case (32)   ! (n,nd)
      nn=1; nc=1; zc(1)=1; ac(1)=2
   case (33)   ! (n,nt)
      nn=1; nc=1; zc(1)=1; ac(1)=3
   case (34)   ! (n,n3He)
      nn=1; nc=1; zc(1)=2; ac(1)=3
   case (35)   ! (n,nd2a)
      nn=1; nc=3; zc(1)=1; ac(1)=2; zc(2:3)=2; ac(2:3)=4
   case (36)   ! (n,nt2a)
      nn=1; nc=3; zc(1)=1; ac(1)=3; zc(2:3)=2; ac(2:3)=4
   case (41)   ! (n,2np)
      nn=2; nc=1; zc(1)=1; ac(1)=1
   case (42)   ! (n,3np)
      nn=3; nc=1; zc(1)=1; ac(1)=1
   case (44)   ! (n,n2p)
      nn=1; nc=2; zc(1:2)=1; ac(1:2)=1
   case (45)   ! (n,npa)
      nn=1; nc=2; zc(1)=1; ac(1)=1; zc(2)=2; ac(2)=4
   case (108)  ! (n,2a)
      nc=2; zc(1:2)=2; ac(1:2)=4
   case (109)  ! (n,3a)
      nc=3; zc(1:3)=2; ac(1:3)=4
   case (111)  ! (n,2p)
      nc=2; zc(1:2)=1; ac(1:2)=1
   case (112)  ! (n,pa)
      nc=2; zc(1)=1; ac(1)=1; zc(2)=2; ac(2)=4
   case (113)  ! (n,t2a)
      nc=3; zc(1)=1; ac(1)=3; zc(2:3)=2; ac(2:3)=4
   case (114)  ! (n,d2a)
      nc=3; zc(1)=1; ac(1)=2; zc(2:3)=2; ac(2:3)=4
   case (115)  ! (n,pd)
      nc=2; zc(1)=1; ac(1)=1; zc(2)=1; ac(2)=2
   case (116)  ! (n,pt)
      nc=2; zc(1)=1; ac(1)=1; zc(2)=1; ac(2)=3
   case (117)  ! (n,da)
      nc=2; zc(1)=1; ac(1)=2; zc(2)=2; ac(2)=4
   case (16)  ! (n,2n)
      nn=2
   case (17)  ! (n,3n)
      nn=3
   case (37)  ! (n,4n)
      nn=4
   case (152)  ! (n,5n)
      nn=5
   case (153)  ! (n,6n)
      nn=6
   case (154)  ! (n,2nt)
      nn=2; nc=1; zc(1)=1; ac(1)=3
   case (155)  ! (n,ta)
      nc=2; zc(1)=1; ac(1)=3; zc(2)=2; ac(2)=4
   case (156)  ! (n,4np)
      nn=4; nc=1; zc(1)=1; ac(1)=1
   case (157)  ! (n,3nd)
      nn=3; nc=1; zc(1)=1; ac(1)=2
   case (158)  ! (n,nda)
      nn=1; nc=2; zc(1)=1; ac(1)=2; zc(2)=2; ac(2)=4
   case (159)  ! (n,2npa)
      nn=2; nc=2; zc(1)=1; ac(1)=1; zc(2)=2; ac(2)=4
   case (160)  ! (n,7n)
      nn=7
   case (161)  ! (n,8n)
      nn=8
   case (162)  ! (n,5np)
      nn=5; nc=1; zc(1)=1; ac(1)=1
   case (163)  ! (n,6np)
      nn=6; nc=1; zc(1)=1; ac(1)=1
   case (164)  ! (n,7np)
      nn=7; nc=1; zc(1)=1; ac(1)=1
   case (165)  ! (n,4na)
      nn=4; nc=1; zc(1)=2; ac(1)=4
   case (166)  ! (n,5na)
      nn=5; nc=1; zc(1)=2; ac(1)=4
   case (167)  ! (n,6na)
      nn=6; nc=1; zc(1)=2; ac(1)=4
   case (168)  ! (n,7na)
      nn=7; nc=1; zc(1)=2; ac(1)=4
   case (169)  ! (n,4nd)
      nn=4; nc=1; zc(1)=1; ac(1)=2
   case (170)  ! (n,5nd)
      nn=5; nc=1; zc(1)=1; ac(1)=2
   case (171)  ! (n,6nd)
      nn=6; nc=1; zc(1)=1; ac(1)=2
   case (172)  ! (n,3nt)
      nn=3; nc=1; zc(1)=1; ac(1)=3
   case (173)  ! (n,4nt)
      nn=4; nc=1; zc(1)=1; ac(1)=3
   case (174)  ! (n,5nt)
      nn=5; nc=1; zc(1)=1; ac(1)=3
   case (175)  ! (n,6nt)
      nn=6; nc=1; zc(1)=1; ac(1)=3
   case (176)  ! (n,2n3He)
      nn=2; nc=1; zc(1)=2; ac(1)=3
   case (177)  ! (n,3n3He)
      nn=3; nc=1; zc(1)=2; ac(1)=3
   case (178)  ! (n,4n3He)
      nn=4; nc=1; zc(1)=2; ac(1)=3
   case (179)  ! (n,3n2p)
      nn=3; nc=2; zc(1)=1; ac(1)=1; zc(2)=1; ac(2)=1
   case (180)  ! (n,3n2a)
      nn=3; nc=2; zc(1)=2; ac(1)=4; zc(2)=2; ac(2)=4
   case (181)  ! (n,3npa)
      nn=3; nc=2; zc(1)=1; ac(1)=1; zc(2)=2; ac(2)=4
   case (182)  ! (n,dt)
      nc=2; zc(1)=1; ac(1)=2; zc(2)=1; ac(2)=3
   case (183)  ! (n,npd)
      nn=1; nc=2; zc(1)=1; ac(1)=1; zc(2)=1; ac(2)=2
   case (184)  ! (n,npt)
      nn=1; nc=2; zc(1)=1; ac(1)=1; zc(2)=1; ac(2)=3
   case (185)  ! (n,ndt)
      nn=1; nc=2; zc(1)=1; ac(1)=2; zc(2)=1; ac(2)=3
   case (186)  ! (n,np3He)
      nn=1; nc=2; zc(1)=1; ac(1)=1; zc(2)=2; ac(2)=3
   case (187)  ! (n,nd3He)
      nn=1; nc=2; zc(1)=1; ac(1)=2; zc(2)=2; ac(2)=3
   case (188)  ! (n,nt3He)
      nn=1; nc=2; zc(1)=1; ac(1)=3; zc(2)=2; ac(2)=3
   case (189)  ! (n,nta)
      nn=1; nc=2; zc(1)=1; ac(1)=3; zc(2)=2; ac(2)=4
   case (190)  ! (n,2n2p)
      nn=2; nc=2; zc(1)=1; ac(1)=1; zc(2)=1; ac(2)=1
   case (191)  ! (n,p3He)
      nc=2; zc(1)=1; ac(1)=1; zc(2)=2; ac(2)=3
   case (192)  ! (n,d3He)
      nc=2; zc(1)=1; ac(1)=2; zc(2)=2; ac(2)=3
   case (193)  ! (n,3Hea)
      nc=2; zc(1)=2; ac(1)=3; zc(2)=2; ac(2)=4
   case (194)  ! (n,4n2p)
      nn=4; nc=2; zc(1)=1; ac(1)=1; zc(2)=1; ac(2)=1
   case (195)  ! (n,4n2a)
      nn=4; nc=2; zc(1)=2; ac(1)=4; zc(2)=2; ac(2)=4
   case (196)  ! (n,4npa)
      nn=4; nc=2; zc(1)=1; ac(1)=1; zc(2)=2; ac(2)=4
   case (197)  ! (n,3p)
      nc=3; zc(1)=1; ac(1)=1; zc(2)=1; ac(2)=1; zc(3)=1; ac(3)=1
   case (198)  ! (n,n3p)
      nn=1; nc=3; zc(1)=1; ac(1)=1; zc(2)=1; ac(2)=1; zc(3)=1; ac(3)=1
   case (199)  ! (n,3n2pa)
      nn=3; nc=3; zc(1)=1; ac(1)=1; zc(2)=1; ac(2)=1; zc(3)=2; ac(3)=4
   case (200)  ! (n,5n2p)
      nn=5; nc=2; zc(1)=1; ac(1)=1; zc(2)=1; ac(2)=1
   case (875:891)  ! (n,2n) to excited states
      nn=2
   end select
   return
   end subroutine mulpart

   subroutine capmul(ee,dame,q,za,awr,mtd)
   !-------------------------------------------------------------------
   ! Damage energy for reactions emitting several particles without
   ! File 6: (n,np), (n,na), (n,2na), ... given with File 5 (conbar)
   ! and (n,2a), (n,pa), ... (mt108-117, capdam).  The available
   ! energy eavail=q+awr*e/(awr+1) is shared by sequential emission:
   ! first the neutrons, each with the mean of the evaporation
   ! spectrum x*exp(-x/T) on [0,U], T=sqrt(U/a), a=A/8 1/MeV; then the
   ! charged particles following the NJOY2016 method, each with
   ! min(vc,U/n), vc=1.029*z*Ze/(a**(1/3)+(Ae-1)**(1/3)) MeV for the
   ! emitting nucleus (Ze,Ae) and n the number of charged particles
   ! still to be emitted.  For one charged particle and no neutron this
   ! is the capdam result.  The momenta of all emitted particles are
   ! added incoherently (isotropic, uncorrelated): K=sum(m_i*eps_i)/Ar,
   ! er=B+K-2*sqrt(B*K)*mu, B=Ar*e/(awr+1)**2, averaged over mu with a
   ! 4-point Gauss rule; the Lindhard partition is computed for the
   ! true residual nucleus (Zr,Ar).  The damage of each emitted charged
   ! particle (lab energy eps+bc+2*sqrt(eps*bc)*mu, bc=a*e/(awr+1)**2)
   ! is added, as sixbar does for File 6 particles.  Neutron-only reactions (nc=0,
   ! nn>1) without File 4/5/6 data are treated with the neutron part
   ! only.
   !-------------------------------------------------------------------
   ! externals
   integer::mtd
   real(kr)::ee,dame,q,za,awr
   ! internals
   integer::nn,nc,i,iq,zc(3),ac(3),ncp
   real(kr)::z,ea,u,ae,ze,ex,vc,sk,b,er,bc,ep,dum,epc(3)
   real(kr),dimension(4),parameter::qp=(/&
      -.86114e0_kr,-.33998e0_kr,.33998e0_kr,.86114e0_kr/)
   real(kr),dimension(4),parameter::qw=(/&
     .34785e0_kr,.65215e0_kr,.65215e0_kr,.34785e0_kr/)
   real(kr),parameter::third=.333333333e0_kr
   real(kr),parameter::econ=1.029e6_kr
   real(kr),parameter::zero=0
   dame=0
   if (ee.eq.zero) return
   call mulpart(mtd,nn,nc,zc,ac)
   if (nc.eq.0.and.nn.lt.2) return
   z=nint(za/1000)
   ea=q+awr*ee/(awr+1)
   if (ea.le.zero) return
   u=ea
   ae=awr+1
   ze=z
   sk=0
   ! neutrons
   do i=1,nn
      ex=evapmean(u,ae-1)
      sk=sk+ex
      u=u-ex
      ae=ae-1
   enddo
   ! charged particles (NJOY2016 method)
   do i=1,nc
      if (ae-ac(i).le.zero) return
      vc=econ*zc(i)*ze/(ac(i)**third+(ae-1)**third)
      ex=min(vc,u/(nc-i+1))
      if (ex.lt.zero) ex=0
      epc(i)=ex
      sk=sk+ac(i)*ex
      u=u-ex
      ae=ae-ac(i)
      ze=ze-zc(i)
   enddo
   if (ae.lt.0.5_kr.or.ze.lt.0.5_kr) return
   sk=sk/ae
   b=ae*ee/((awr+1)*(awr+1))
   ! set the Lindhard constants for this residual nucleus
   er=df(zero,ze,ae,z,awr)
   do iq=1,4
      er=b+sk-2*sqrt(b*sk)*qp(iq)
      dame=dame+qw(iq)*df(er,ze,ae,z,awr)/2
   enddo
   ! damage of each emitted charged particle itself (as in sixbar for
   ! File 6 particles): lab energy ep=ex+bc+2*sqrt(ex*bc)*mu with
   ! bc=a*e/(awr+1)**2, Lindhard partition of the particle in the target
   ncp=nc
   do i=1,ncp
      bc=ac(i)*ee/((awr+1)*(awr+1))
      dum=df(zero,real(zc(i),kr),real(ac(i),kr),z,awr)
      do iq=1,4
         ep=epc(i)+bc+2*sqrt(epc(i)*bc)*qp(iq)
         dame=dame+qw(iq)*df(ep,real(zc(i),kr),real(ac(i),kr),z,awr)/2
      enddo
   enddo
   ! restore the Lindhard constants of the residual nucleus
   if (ncp.gt.0) dum=df(zero,ze,ae,z,awr)
   return
   end subroutine capmul

   subroutine disbar(ee,ebar,dame,q,za,awr,ntape,matd,mtd,a,nwa)
   !-------------------------------------------------------------------
   ! Compute ebar and damage energy for discrete scattering.
   ! For discrete charged-particle levels (File 4 only) the damage of
   ! the emitted particle is added to the recoil damage, as sixbar does
   ! for the particles given in File 6.
   ! for discrete charged-particle levels the Lindhard
   ! partition of the recoil uses the residual nucleus, zr=z-zp and
   ! ar=awr+1-awp (NJOY2016.80 used the target z and awr), as capdam
   ! does for the continuum; the kinematics already used the residual
   ! mass.  The File 4 distribution is used as given in the angular
   ! integrations (no clipping of negative values, as in tabdmn).
   ! (DLA, IAEA, Oct 2026)
   !-------------------------------------------------------------------
   use endf ! provides iverf,terp1
   use mathm ! provides legndr
   use physics !get global physics and light particle mass constants
   use util ! provides mess
   ! externals
   integer::ntape,matd,mtd,nwa
   real(kr)::ee,ebar,dame,q,za,awr,a(*)
   ! internals
   integer::imiss,i,mfd,nld,lcd,idis,iz,iq,il
   real(kr)::e,thresh,awp,afact,arat,en,cn,el,cl,z
   real(kr)::daml,damn,enx,enext,u,f,e2,ce,r,b,g,test,wbar
   integer::lwmx
   character(60)::strng
   real(kr),parameter::wmax=0.999999e0_kr
   real(kr)::fl(65),p(65),zpd,ep,dum,zr,ar
   integer,parameter::nq=64
   real(kr),dimension(64),parameter::qp=(/&
     -9.99305042E-01_kr,-9.96340117E-01_kr,-9.91013371E-01_kr,&
     -9.83336254E-01_kr,-9.73326828E-01_kr,-9.61008800E-01_kr,&
     -9.46411375E-01_kr,-9.29569172E-01_kr,-9.10522137E-01_kr,&
     -8.89315446E-01_kr,-8.65999398E-01_kr,-8.40629296E-01_kr,&
     -8.13265315E-01_kr,-7.83972359E-01_kr,-7.52819907E-01_kr,&
     -7.19881850E-01_kr,-6.85236313E-01_kr,-6.48965471E-01_kr,&
     -6.11155355E-01_kr,-5.71895646E-01_kr,-5.31279464E-01_kr,&
     -4.89403146E-01_kr,-4.46366017E-01_kr,-4.02270158E-01_kr,&
     -3.57220158E-01_kr,-3.11322872E-01_kr,-2.64687162E-01_kr,&
     -2.17423644E-01_kr,-1.69644420E-01_kr,-1.21462819E-01_kr,&
     -7.29931218E-02_kr,-2.43502927E-02_kr, 2.43502927E-02_kr,&
      7.29931218E-02_kr, 1.21462819E-01_kr, 1.69644420E-01_kr,&
      2.17423644E-01_kr, 2.64687162E-01_kr, 3.11322872E-01_kr,&
      3.57220158E-01_kr, 4.02270158E-01_kr, 4.46366017E-01_kr,&
      4.89403146E-01_kr, 5.31279464E-01_kr, 5.71895646E-01_kr,&
      6.11155355E-01_kr, 6.48965471E-01_kr, 6.85236313E-01_kr,&
      7.19881850E-01_kr, 7.52819907E-01_kr, 7.83972359E-01_kr,&
      8.13265315E-01_kr, 8.40629296E-01_kr, 8.65999398E-01_kr,&
      8.89315446E-01_kr, 9.10522137E-01_kr, 9.29569172E-01_kr,&
      9.46411375E-01_kr, 9.61008800E-01_kr, 9.73326828E-01_kr,&
      9.83336254E-01_kr, 9.91013371E-01_kr, 9.96340117E-01_kr,&
      9.99305042E-01_kr/)
   real(kr),dimension(64),parameter::qw=(/&
      1.78328072E-03_kr, 4.14703326E-03_kr, 6.50445797E-03_kr,&
      8.84675983E-03_kr, 1.11681395E-02_kr, 1.34630479E-02_kr,&
      1.57260305E-02_kr, 1.79517158E-02_kr, 2.01348232E-02_kr,&
      2.22701738E-02_kr, 2.43527026E-02_kr, 2.63774697E-02_kr,&
      2.83396726E-02_kr, 3.02346571E-02_kr, 3.20579284E-02_kr,&
      3.38051618E-02_kr, 3.54722133E-02_kr, 3.70551285E-02_kr,&
      3.85501532E-02_kr, 3.99537411E-02_kr, 4.12625632E-02_kr,&
      4.24735151E-02_kr, 4.35837245E-02_kr, 4.45905582E-02_kr,&
      4.54916279E-02_kr, 4.62847966E-02_kr, 4.69681828E-02_kr,&
      4.75401657E-02_kr, 4.79993886E-02_kr, 4.83447622E-02_kr,&
      4.85754674E-02_kr, 4.86909570E-02_kr, 4.86909570E-02_kr,&
      4.85754674E-02_kr, 4.83447622E-02_kr, 4.79993886E-02_kr,&
      4.75401657E-02_kr, 4.69681828E-02_kr, 4.62847966E-02_kr,&
      4.54916279E-02_kr, 4.45905582E-02_kr, 4.35837245E-02_kr,&
      4.24735151E-02_kr, 4.12625632E-02_kr, 3.99537411E-02_kr,&
      3.85501532E-02_kr, 3.70551285E-02_kr, 3.54722133E-02_kr,&
      3.38051618E-02_kr, 3.20579284E-02_kr, 3.02346571E-02_kr,&
      2.83396726E-02_kr, 2.63774697E-02_kr, 2.43527026E-02_kr,&
      2.22701738E-02_kr, 2.01348232E-02_kr, 1.79517158E-02_kr,&
      1.57260305E-02_kr, 1.34630479E-02_kr, 1.11681395E-02_kr,&
      8.84675983E-03_kr, 6.50445797E-03_kr, 4.14703326E-03_kr,&
      1.78328072E-03_kr/)
   real(kr),parameter::edis=25.e0_kr
   real(kr),parameter::step=1.1e0_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::zero=0
   save imiss,thresh,awp,afact,arat
   save en,cn,el,cl,z,daml,damn,enx,enext,zpd,zr,ar
   save lwmx

   !--initialize if e=0.
   e=ee
   if (e.eq.zero) then
      lwmx=0
      imiss=0
      do i=1,nmiss4
         if (miss4(i).eq.mtd) imiss=1
      enddo
      if (imiss.ne.1) then
         mfd=4
         nld=60
         lcd=2
         call hgtfle(e,enext,idis,fl,nld,lcd,matd,mfd,mtd,&
           ntape,a,nwa)
      else
         nld=1
         enext=etop
      endif
      thresh=((awr+1)/awr)*(-q)
      awp=1
      if (iverf.le.5.and.mtd.ge.700.and.mtd.lt.718) awp=pnratio
      if (iverf.le.5.and.mtd.ge.720.and.mtd.lt.738) awp=dnratio
      if (iverf.le.5.and.mtd.ge.740.and.mtd.lt.758) awp=tnratio
      if (iverf.le.5.and.mtd.ge.760.and.mtd.lt.778) awp=hnratio
      if (iverf.le.5.and.mtd.ge.780.and.mtd.lt.798) awp=anratio
      if (iverf.ge.6.and.mtd.ge.600.and.mtd.lt.649) awp=pnratio
      if (iverf.ge.6.and.mtd.ge.650.and.mtd.lt.699) awp=dnratio
      if (iverf.ge.6.and.mtd.ge.700.and.mtd.lt.749) awp=tnratio
      if (iverf.ge.6.and.mtd.ge.750.and.mtd.lt.799) awp=hnratio
      if (iverf.ge.6.and.mtd.ge.800.and.mtd.lt.849) awp=anratio
      afact=awp/((awr+1)**2)
      arat=awp/(awr+1-awp)
      ! charge of the emitted particle (discrete charged-particle levels)
      zpd=0
      if (iverf.le.5.and.mtd.ge.700.and.mtd.lt.758) zpd=1
      if (iverf.le.5.and.mtd.ge.760.and.mtd.lt.798) zpd=2
      if (iverf.ge.6.and.mtd.ge.600.and.mtd.lt.749) zpd=1
      if (iverf.ge.6.and.mtd.ge.750.and.mtd.lt.849) zpd=2
      en=0
      cn=0
      el=0
      cl=0
      enx=0
      if (idame.ne.0) then
         iz=nint(za/1000)
         z=iz
         ! recoil nucleus: target for neutron levels, residual
         ! (z-zp, awr+1-awp) for charged-particle levels
         zr=z
         ar=awr
         if (zpd.gt.zero) then
            zr=z-zpd
            ar=awr+1-awp
         endif
         damn=df(e,zr,ar,z,awr)
         daml=damn
         if (mtd.eq.2) then
            daml=0
            enx=edis*arat/(4*afact)
            damn=0
         endif
      endif
      return
   endif

   !--normal entry.
   if (ee.le.en*(1+small)) go to 130
   el=en
   cl=cn
   if (idame.ne.0) daml=damn
   e=step*el
   if (enext.lt.e*(1-small)) e=enext
   if (ee.gt.e*(1+small)) e=ee
   r=1
   if (mtd.ne.2) then
      if (thresh.ge.e*(1-small)) then
         r=0
      else
         r=sqrt(1-thresh/e)
      endif
   endif
   if (imiss.ne.1) then
      mfd=4
      nld=60
      lcd=2
      call hgtfle(e,enext,idis,fl,nld,lcd,matd,mfd,mtd,ntape,a,nwa)
   else
      nld=1
      fl(1)=1
      fl(2)=0
   endif
   en=e
   b=r*sqrt(awr/arat)
   g=b*arat
   cn=0
   test=1
   test=test/1000000
   if (abs(awp-1).gt.test) go to 120

   !--compute ebar
   wbar=fl(2)
   ! mean cosine limited to 1-1e-6 (was qp(64)=0.999305, which cut
   ! forward-peaked distributions); message once per mt
   ! (DLA, IAEA, 20 Oct 2026)
   if (wbar.gt.wmax) then
      wbar=wmax
      if (lwmx.eq.0) then
         write(strng,'(''mt'',i4,'' mean cosine limited to 0.999999'',&
           &'' at'',1p,e10.3,'' ev'')') mtd,e
         call mess('disbar',strng,' ')
         lwmx=1
      endif
   endif
   cn=(1+2*b*wbar+b*b)*afact
   if (idame.eq.0) go to 130
   if (e.lt.enx*(1-small)) go to 130

   !--angle integration by gauss-legendre method.
  120 continue
   dame=0
   ! charged-particle levels come here also without damage option
   if (idame.eq.0) go to 130
   dum=df(zero,zr,ar,z,awr)
   do iq=1,nq
      u=qp(iq)
      call legndr(u,p,nld)
      f=0
      do il=1,nld
         f=f+(2*il-1)*fl(il)*p(il)/2
      enddo
      e2=e*(1-2*g*u+g*g)*afact/arat
      dame=dame+qw(iq)*f*df(e2,zr,ar,z,awr)
   enddo
   ! damage of the emitted charged particle itself (as in sixbar for
   ! File 6 particles): lab energy ep=e*afact*(1+2*b*u+b*b)
   if (zpd.gt.zero.and.r.gt.zero) then
      dum=df(zero,zpd,awp,z,awr)
      do iq=1,nq
         u=qp(iq)
         call legndr(u,p,nld)
         f=0
         do il=1,nld
            f=f+(2*il-1)*fl(il)*p(il)/2
         enddo
         ep=e*(1+2*b*u+b*b)*afact
         dame=dame+qw(iq)*f*df(ep,zpd,awp,z,awr)
      enddo
      ! restore the Lindhard constants of the recoil nucleus
      dum=df(zero,zr,ar,z,awr)
   endif
   if (dame.lt.zero) dame=zero
   damn=dame

   !--interpolate for results.
  130 continue
   call terp1(el,cl,en,cn,ee,ce,2)
   ebar=ee*ce
   if (idame.gt.0) then
      call terp1(el,daml,en,damn,ee,dame,2)
   endif
   return
   end subroutine disbar

   real(kr) function df(e,zr,ar,zl,al)
   !-------------------------------------------------------------------
   ! Damage function using the Lindhard partition of
   ! energy between atomic and electronic motion.
   ! Call with e=0 for each reaction to precompute the constants.
   !-------------------------------------------------------------------
   ! externals
   real(kr)::e,zr,ar,zl,al
   ! internals
   real(kr)::el,rel,denom,fl,ep,dam
   real(kr),parameter::twothd=.666666667e0_kr
   real(kr),parameter::threeq=.75e0_kr
   real(kr),parameter::sixth=.166666667e0_kr
   real(kr),parameter::onep5=1.5e0_kr
   real(kr),parameter::c1=30.724e0_kr
   real(kr),parameter::c2=.0793e0_kr
   real(kr),parameter::c3=3.4008e0_kr
   real(kr),parameter::c4=.40244e0_kr
   real(kr),parameter::zero=0
   save rel,fl

   if (zr.eq.zero) then
      df=0
   else if (e.le.zero) then
      el=c1*zr*zl*sqrt(zr**twothd+zl**twothd)*(ar+al)/al
      rel=1/el
      denom=(zr**twothd+zl**twothd)**threeq*ar**onep5*sqrt(al)
      fl=c2*zr**twothd*sqrt(zl)*(ar+al)**onep5/denom
      df=0
   else if (e.lt.break) then
      df=0
   else
      ep=e*rel
      dam=e/(1+fl*(c3*ep**sixth+c4*ep**threeq+ep))
      df=dam
   endif
   return
   end function df

   subroutine conbar(e,ebar,dame,nin,nscr,matd,mtd,&
     c,ncmax,b,nbmax,yld,q)
   !-------------------------------------------------------------------
   ! Calculate neutron ebar and damage energy
   ! for continuum scattering reactions.
   ! Damage: reactions with charged particles ((n,np), (n,na), ...)
   ! use capmul; neutron-only reactions use tabdmn/anadmn with the
   ! evaluated neutron spectrum and the momentum of all nu emitted
   ! neutrons; fission keeps the original tabdam/anadam treatment.
   !-------------------------------------------------------------------
   use util ! provides error
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nscr,matd,mtd,ncmax,nbmax
   real(kr)::e,ebar,dame,c(ncmax),b(nbmax),yld,q
   ! internals
   integer::l,mf1,mt1,mfd,idis,nb,nw,iz,il,ik,nktot
   integer::ip,ir,idisc,nr,np,ll,nnr,inn,ie,idone,ii,mtt
   integer::lf,nne,lnow,iraw,ne,nnm,ncm,zcm(3),acm(3),nun,iout
   logical::lfis
   character(60)::strng
   real(kr)::enext,t1,t2,ehi,fhi,dhi,pe,eihi,etp,f,d,s
   real(kr)::awr,awfac,aw1fac,z,za
   integer,parameter::nkmax=12
   integer::loc(nkmax)
   real(kr)::d1(nkmax),d2(nkmax),e1(nkmax),e2(nkmax)
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::step=1.5e0_kr
   real(kr),parameter::stepn=1.1e0_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::zero=0
   save nktot,awr,awfac,aw1fac,z,za,mtt,nun,lfis
   save loc,d1,d2,e1,e2,mf1,mt1

   !--initialize if e=0.
   if (e.eq.zero) then
      l=1
      if ((mtd.lt.18.or.mtd.gt.21).and.mtd.ne.38) then
         yld=1
         if (mtd.eq.16.or.mtd.eq.30) yld=2
         if (mtd.eq.24.or.mtd.eq.26) yld=2
         if (mtd.eq.11.or.mtd.eq.41) yld=2
         if (mtd.eq.17.or.mtd.eq.25) yld=3
         if (mtd.eq.38.or.mtd.eq.42) yld=3
         if (mtd.eq.37) yld=4
         ! mt152-200, 875-891: number of emitted neutrons
         if (mtd.gt.100) then
            call mulpart(mtd,nnm,ncm,zcm,acm)
            yld=max(1,nnm)
         endif
      else
         mf1=1
         mt1=452
         if (nply.gt.0) mt1=456
         ! delayed neutrons are treated the same as fast neutrons.
         ! this will cause a slight overestimate of en for fission.
         ! ... but when nply.ne.0, use prompt neutron data only
         call hgtyld(e,enext,idis,yld,matd,mf1,mt1,nscr,b,nbmax)
      endif
      ! File 4 angular distribution of the emitted neutrons, used to
      ! correct the isotropic recoil of tabdmn/anadmn (neutron-only,
      ! non-fission reactions).  Must precede the reading of File 5.
      n4e=0
      call mulpart(mtd,nnm,ncm,zcm,acm)
      if (idame.gt.0.and.ncm.eq.0.and.&
        (mtd.lt.18.or.mtd.gt.21).and.mtd.ne.38) call anitab(nin,matd,mtd)
      mfd=5
      call findf(matd,mfd,mtd,nin)
      call contio(nin,0,0,c(l),nb,nw)
      nktot=nint(c(5))
      if (nktot.gt.nkmax) call error('conbar','nktot gt nkmax.',' ')
      awr=c(2)
      awfac=1/awr
      aw1fac=1/(awr+1)
      iz=nint(c(1)/1000)
      z=iz
      za=c(1)
      ! damage treatment: mtt=mtd for reactions with charged particles
      ! (capmul); otherwise the neutron spectrum is used with nun
      ! neutrons (tabdmn/anadmn).  NJOY2016.80 treated (n,np), (n,na),
      ! (n,nd), (n,nt), (n,n3He) as (n,p), ... (neutron ignored, residual
      ! mass A) and the other (n,nX) reactions and (n,2n), (n,3n), (n,4n)
      ! as the emission of a single neutron (residual mass A).
      mtt=0
      call mulpart(mtd,nnm,ncm,zcm,acm)
      if (ncm.gt.0) mtt=mtd
      nun=max(1,nnm)
      lfis=((mtd.ge.18.and.mtd.le.21).or.mtd.eq.38)
      ne=0
      il=1
      ik=0
      do while (ik.lt.nktot)
         ik=ik+1
         d2(ik)=0
         e2(ik)=0
         loc(ik)=il
         l=il
         call tab1io(nin,0,0,c(l),nb,nw)
         l=l+nw
         do while (nb.ne.0)
            call moreio(nin,0,0,c(l),nb,nw)
            l=l+nw
         enddo
         lf=nint(c(il+3))
         il=l

         !--analytic subsection.
         if (lf.gt.1) then
            ! anadmn computes the damage for LF=5, 7, 9,
            ! 11 and 12 (the other analytic laws are not ENDF-6 laws)
            ! (DLA, IAEA, Oct 2026)
            if (idame.gt.0.and.mtt.eq.0.and.lf.ne.5.and.lf.ne.7.and.&
              lf.ne.9.and.lf.ne.11.and.lf.ne.12) then
               write(strng,'(''mf5, mt'',i3,'' analytic law lf='',i2)')&
                 mtd,lf
               call mess('conbar',strng,&
                 'damage of this subsection not computed')
            endif
            if (lf.eq.3) then
               t1=(awr*awr+1)*aw1fac*aw1fac
               t2=awr*aw1fac
            else
               call tab1io(nin,0,0,c(l),nb,nw)
               l=l+nw
               do while (nb.ne.0)
                  if (l.gt.ncmax) call error('conbar',&
                    'insufficient storage for raw endf data.',' ')
                  call moreio(nin,0,0,c(l),nb,nw)
                  l=l+nw
               enddo
               il=l
               if (lf.eq.5.or.lf.eq.11) then
                  call tab1io(nin,0,0,c(l),nb,nw)
                  l=l+nw
                  do while (nb.ne.0)
                     if (l.gt.ncmax) call error('conbar',&
                       'insufficient storage for raw endf data.',' ')
                     call moreio(nin,0,0,c(l),nb,nw)
                     l=l+nw
                  enddo
                  il=l
               endif
            endif

         !--tabulated subsection
         else
            call tab2io(nin,0,0,c(l),nb,nw)
            ne=nint(c(l+5))
            l=l+nw
            nne=0
            lnow=l
            iraw=lnow+3*ne

            !--read spectrum for each incident energy
            !--compute heating and damage contributions
            do while (nne.lt.ne)
               l=iraw
               call tab1io(nin,0,0,c(l),nb,nw)
               l=l+nw
               do while (nb.ne.0)
                  if (l.gt.ncmax) call error('conbar',&
                    'insufficient storage for raw endf data.',' ')
                  call moreio(nin,0,0,c(l),nb,nw)
                  l=l+nw
               enddo
               nne=nne+1
               ehi=c(iraw+1)
               call tabbar(fhi,c(iraw),lf)
               dhi=0
               if (idame.ne.0) then
                  if (mtt.gt.0) then
                     call capmul(ehi,dhi,q,za,awr,mtt)
                  else if (lfis) then
                     ! fission: recoil of the fissioning
                     ! nucleus (mass awr) from one neutron, as tabdam,
                     ! with the integration of tabdmn (normalized,
                     ! interpolation laws of the data) (DLA, IAEA, Oct 2026)
                     call tabdmn(ehi,dhi,z,awr,1,c(iraw))
                  else
                     ! exact recoil of the nun neutrons (DLA, IAEA, Oct 2026)
                     lcvx=nun.le.16
                     call tabdmn(ehi,dhi,z,awr,nun,c(iraw))
                     lcvx=.false.
                  endif
               endif
               c(lnow+3*(nne-1))=ehi
               c(lnow+3*(nne-1)+1)=fhi
               c(lnow+3*(nne-1)+2)=dhi
            enddo
            il=lnow+3*ne
         endif
      enddo
      return
   endif

   !--normal entry.
   enext=emax
   idis=0
   ebar=0
   dame=0

   !--get energy dependent yld.
   if (mtd.eq.18.or.mtd.eq.19.or.mtd.eq.20.or.mtd.eq.21.or.mtd.eq.38) then
      mf1=1
      call hgtyld(e,enext,idis,yld,matd,mf1,mt1,nscr,b,nbmax)
   endif

   !--loop over subsections.
   do ik=1,nktot
      lnow=loc(ik)
      lf=nint(c(lnow+3))

      !--interpolate for fractional probability.
      ip=2
      ir=1
      call terpa(pe,e,eihi,idisc,c(lnow),ip,ir)
      if (abs(eihi-enext).lt.small*enext.and.idisc.gt.idis) idis=idisc
      if (eihi.lt.enext*(1-small)) idis=idisc
      if (eihi.lt.enext*(1-small)) enext=eihi
      nr=nint(c(lnow+4))
      np=nint(c(lnow+5))
      ll=lnow+6+2*nr+2*(np-1)
      etp=c(ll)
      if (pe.gt.zero) then

         !--tabulated subsection.
         !--interpolate for results
         if (lf.eq.1) then
            lnow=lnow+6+2*nr+2*np
            ne=nint(c(lnow+5))
            nnr=nint(c(lnow+4))
            inn=2
            lnow=lnow+6+2*nnr
            ! the table has ne entries (e,ebar,dame) at
            ! lnow+3*(ie-1); outside its range the first or the last
            ! entry is used (iout=1).  NJOY2016.80 tested entry ie+1
            ! up to ie=ne, i.e. one entry beyond the table, and used it
            ! to interpolate above the last incident energy; below the
            ! first one it extrapolated linearly. (DLA, IAEA, Oct 2026)
            ii=0
            iout=0
            if (ne.le.1.or.e.lt.c(lnow)*(1-small)) then
               ii=lnow
               iout=1
            else
               do ie=2,ne
                  if (c(lnow+3*(ie-1)).ge.e*(1-small)) then
                     ii=lnow+3*(ie-2)
                     exit
                  endif
               enddo
               if (ii.eq.0) then
                  ii=lnow+3*(ne-1)
                  iout=1
               endif
            endif
            if (iout.eq.0) then
               call terp1(c(ii),c(ii+1),c(ii+3),c(ii+4),e,f,inn)
            else
               f=c(ii+1)
            endif
            ebar=ebar+f*pe
            if (idame.ne.0) then
               if (iout.eq.0) then
                  call terp1(c(ii),c(ii+2),c(ii+3),c(ii+5),e,d,inn)
               else
                  d=c(ii+2)
               endif
               dame=dame+d*pe
            endif

         !--analytic subsection.
         else
            call anabar(s,e,c(lnow),t1,t2)
            ebar=ebar+s*pe
            if (idame.ne.0) then
               if (mtt.gt.0) then
                  ! (NJOY2016.80 overwrote dame here: call capdam(e,dame,...))
                  call capmul(e,d,q,za,awr,mtt)
                  dame=dame+d*pe
               else
                  if (e.ge.e2(ik)*(1-small)) then
                     idone=0
                     do while (idone.eq.0)
                        d1(ik)=d2(ik)
                        e1(ik)=e2(ik)
                        ! energy step 1.1 and anadmn; fission
                        ! uses anadmn with one neutron (DLA, IAEA, Oct 2026)
                        e2(ik)=stepn*e1(ik)
                        if (e2(ik).gt.etp*(1+small)) e2(ik)=etp
                        if (e1(ik).eq.zero) e2(ik)=e
                        lcvx=(.not.lfis).and.nun.le.16
                        call anadmn(e2(ik),d2(ik),z,awr,nun,c(lnow))
                        lcvx=.false.
                        if (abs(e2(ik)-etp).lt.small*etp) then
                           idone=1
                        else
                           if (e1(ik).ne.zero.and.&
                             e.le.e2(ik)*(1+small)) idone=1
                        endif
                     enddo
                  endif
                  call terp1(e1(ik),d1(ik),e2(ik),d2(ik),e,d,2)
                  dame=dame+d*pe
               endif
            endif
         endif
      endif
   enddo
   end subroutine conbar

   subroutine hgtyld(e,enext,idis,yld,mat,mf,mt,itape,a,na)
   !-------------------------------------------------------------------
   ! Retrieve or compute yield for requested mat,mf,mt on itape.
   ! If mt=455, retrieve delayed neutron time constants.
   ! Initialize if e=0.
   !-------------------------------------------------------------------
   use util ! provides error
   use endf ! provides endf routines and variables
   ! externals
   integer::idis,mat,mf,mt,itape,na
   real(kr)::e,enext,yld,a(*)
   ! internals
   integer::nb,nw,lnu,lnd,loc,ip,ir,nr,i,nc
   real(kr)::term
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::zero=0
   save lnu,ip,ir

   !--initialize.
   if (e.eq.zero) then
      call findf(mat,mf,mt,itape)
      call contio(itape,0,0,a,nb,nw)
      lnu=nint(a(4))
      if (mth.eq.455) then
         loc=1
         call listio(itape,0,0,a(loc),nb,nw)
         loc=loc+nw
         do while (nb.ne.0)
            call moreio(itape,0,0,a(loc),nb,nw)
            loc=loc+nw
         enddo
         lnd=nint(a(5))
         if (lnd.ne.6.and.lnd.ne.8)&
           call error('hgtyld','illegal lnd, must be 6 or 8',' ')
         loc=1
         call tab1io(itape,0,0,a(loc),nb,nw)
         loc=loc+nw
         if (loc.gt.na) call error('hgtyld','storage exceeded.',' ')
         do while (nb.ne.0)
            call moreio(itape,0,0,a(loc),nb,nw)
            loc=loc+nw
            if (loc.gt.na) call error('hgtyld','storage exceeded.',' ')
         enddo
         na=loc-1
         ip=2
         ir=1
         nr=nint(a(5))
         enext=a(7+2*nr)
      else
         lnd=0
         loc=1
         if (lnu.ne.1) then
            call tab1io(itape,0,0,a(loc),nb,nw)
            loc=loc+nw
            if (loc.gt.na) call error('hgtyld','storage exceeded.',' ')
            do while (nb.ne.0)
               call moreio(itape,0,0,a(loc),nb,nw)
               loc=loc+nw
               if (loc.gt.na) call error('hgtyld','storage exceeded.',' ')
            enddo
            na=loc-1
            ip=2
            ir=1
            nr=nint(a(5))
            enext=a(7+2*nr)
         else
            loc=1
            call listio(itape,0,0,a(loc),nb,nw)
            loc=loc+nw
            na=nw
            do while (nb.ne.0)
               call moreio(itape,0,0,a(loc),nb,nw)
               loc=loc+nw
               na=na+nw
            enddo
            enext=emax
         endif
      endif
      idis=0
      return
   endif

   !--tabulated data.
   if (lnu.ne.1) then
      call terpa(yld,e,enext,idis,a,ip,ir)

   !--polynomial data.
   else
      yld=a(7)
      term=1
      nc=nint(a(5))
      do i=2,nc
         term=term*e
         yld=yld+a(i+6)*term
      enddo
      enext=emax
      idis=0
   endif
   return
   end subroutine hgtyld

   subroutine anabar(s,e,a,t1,t2)
   !-------------------------------------------------------------------
   ! Compute mean energy for an analytic section of File 5.
   !-------------------------------------------------------------------
   use physics ! provides pi
   use endf    ! provides terpa
   ! externals
   real(kr)::s,e,a(*),t1,t2
   ! internals
   integer::lf,nr,np,lnext,ip,ir,idis,nrn,npn,lnxt,i,ipg,irg
   real(kr)::u,theta,enext,b1,b2,b3,expa,test,b,rpi2
   real(kr)::xmx,s0,s1,x,w,g,ep
   integer,parameter::nsim=400
   integer,parameter::nsim5=2000
   real(kr),parameter::thrhaf=1.5e0_kr
   real(kr),parameter::zero=0
   rpi2=sqrt(pi)/2

   !--retrieve common factors.
   s=0
   u=a(1)
   if (e.le.u) return
   lf=nint(a(4))
   nr=nint(a(5))
   np=nint(a(6))
   lnext=7+2*nr+2*np
   ip=2
   ir=1
   call terpa(theta,e,enext,idis,a(lnext),ip,ir)

   !--law 3--discrete energy spectrum.
   if (lf.eq.3) then
      s=e*t1-t2*theta

   !--law 5--general evaporation spectrum.
   else if (lf.eq.5) then
      ! g(x) truncated at x=(e-u)/theta (ENDF-6: the
      ! secondary energy is limited to e-u); the original code used
      ! the mean of the whole table (tabbar) (DLA, IAEA, Oct 2026)
      nrn=nint(a(lnext+4))
      npn=nint(a(lnext+5))
      lnxt=lnext+2*nrn+2*npn+6
      if (theta.gt.zero) then
         xmx=a(lnxt+6+2*nint(a(lnxt+4))+2*(nint(a(lnxt+5))-1))
         xmx=min(xmx,(e-u)/theta)
         if (xmx.gt.zero) then
            s0=0
            s1=0
            ipg=2
            irg=1
            do i=0,nsim5
               x=xmx*i/nsim5
               w=2
               if (mod(i,2).eq.1) w=4
               if (i.eq.0.or.i.eq.nsim5) w=1
               call terpa(g,x,enext,idis,a(lnxt),ipg,irg)
               if (g.lt.zero) g=0
               s0=s0+w*g
               s1=s1+w*g*x
            enddo
            if (s0.gt.zero) s=theta*s1/s0
         endif
      endif

   !--law 7--simple fission spectrum.
   else if (lf.eq.7) then
      if (theta.ne.zero) then
         b1=(e-u)/theta
         b2=sqrt(b1)
         b3=exp(-b1)
         expa=exp(-b2*b2)
         s=theta*(thrhaf-b1*b2*b3/(rpi2*(1-expa*rerfc(b2))-b2*b3))
      endif

   !--law 9--evaporation spectrum.
   else if (lf.eq.9) then
      if (theta.ne.zero) then
         b1=(e-u)/theta
         b2=sqrt(b1)
         b3=exp(-b1)
         test=1
         test=test/1000
         if (b1.ge.test) then
            s=theta*(2-b1*b1*b3/(1-(b1+1)*b3))
         else
            s=4*(e-u)/3
         endif
      endif

   !--law 11--energy-dependent watt spectrum.
   else if (lf.eq.11) then
      nrn=nint(a(lnext+4))
      npn=nint(a(lnext+5))
      lnxt=lnext+6+2*nrn+2*npn
      ip=2
      ir=1
      call terpa(b,e,enext,idis,a(lnxt),ip,ir)
      s=theta*(thrhaf+theta*b/4)
      ! Watt spectrum truncated at e-u (ENDF-6); the
      ! analytic mean above is used when the truncation is negligible
      ! (DLA, IAEA, Oct 2026)
      if (theta.gt.zero.and.b.ge.zero) then
         x=(sqrt(50*theta)+theta*sqrt(b)/2)**2
         if (e-u.lt.x) then
            x=e-u
            s0=0
            s1=0
            do i=0,nsim
               ep=x*i/nsim
               w=2
               if (mod(i,2).eq.1) w=4
               if (i.eq.0.or.i.eq.nsim) w=1
               g=exp(-ep/theta)*sinh(sqrt(b*ep))
               s0=s0+w*g
               s1=s1+w*g*ep
            enddo
            if (s0.gt.zero) s=s1/s0
         endif
      endif

   !--law 12--energy dependent fission neutron spectrum (madland-nix)
   else if (lf.eq.12) then
          s=((a(lnext)+a(lnext+1))/2)+(4*theta/3)
   endif
   return

   contains
      real(kr) function r(z)
      real(kr)::z
      real(kr),parameter::a0=.3275911e0_kr
      r=1/(1+a0*abs(z))
      return
      end function r

      real(kr) function rerfc(z)
      real(kr)::z
      real(kr),parameter::a1=.254829592e0_kr
      real(kr),parameter::a2=-.284496736e0_kr
      real(kr),parameter::a3=1.421413741e0_kr
      real(kr),parameter::a4=-1.453152027e0_kr
      real(kr),parameter::a5=1.061405429e0_kr
      rerfc=r(z)*(a1+r(z)*(a2+r(z)*(a3+r(z)*(a4+a5*r(z)))))
      return
      end function rerfc

   end subroutine anabar

   subroutine anadam(e,d,z,awr,a)
   !-------------------------------------------------------------------
   ! Compute damage energy for an analytic section of File 5
   !-------------------------------------------------------------------
   use endf ! provides terpa
   ! externals
   real(kr)::e,d,z,awr,a(*)
   ! internals
   integer::nr,np,lf,lnext,ip,ir,idis,i,iq,j,iflag
   real(kr)::u,awfac,test,theta,enext,ep,er,xlast,ylast,xm,ym,yt
   integer,parameter::imax=10
   real(kr)::x(imax),y(imax)
   integer,parameter::nq=4
   real(kr),dimension(nq),parameter::qp=(/&
     -.86114e0_kr,-.33998e0_kr,.33998e0_kr,.86114e0_kr/)
   real(kr),dimension(nq),parameter::qw=(/&
     .34785e0_kr,.65215e0_kr,.65215e0_kr,.34785e0_kr/)
   real(kr),parameter::eps=.05e0_kr
   real(kr),parameter::delta=1.e-7_kr

   !--retrieve spectrum temperature
   d=0
   u=a(1)
   nr=nint(a(5))
   np=nint(a(6))
   lf=nint(a(4))
   awfac=1/awr
   if (lf.ne.9) return
   if (e.le.u) return
   if ((e-u).lt.delta*e) return
   test=10
   if ((e-u).lt.test) return
   lnext=7+2*nr+2*np
   ip=2
   ir=1
   call terpa(theta,e,enext,idis,a(lnext),ip,ir)

   !--adaptive energy integration
   !--gauss-legendre angle integration
   i=4
   do while (i.gt.1)
      i=i-1
      if (i.eq.3) ep=1
      if (i.eq.2) ep=(e-u)/2
      if (i.eq.2.and.theta.lt.e-u) ep=theta
      if (i.eq.1) ep=e-u
      x(i)=ep
      y(i)=0
      do iq=1,nq
         er=(e-2*sqrt(e*ep)*qp(iq)+ep)*awfac
         y(i)=y(i)+qw(iq)*df(er,z,awr,z,awr)/2
      enddo
      y(i)=y(i)*sed(e,ep,theta,u)
   enddo
   i=3
   j=0
   xlast=0
   ylast=0
   do while (i.gt.1)
      ! test for convergence
      iflag=0
      if (i.gt.1.and.i.lt.imax) then
         xm=(x(i-1)+x(i))/2
         ym=(y(i-1)+y(i))/2
         yt=0
         do iq=1,nq
            er=(e-2*sqrt(e*xm)*qp(iq)+xm)*awfac
            yt=yt+qw(iq)*df(er,z,awr,z,awr)/2
         enddo
         yt=yt*sed(e,xm,theta,u)
         if (abs(yt-ym).gt.eps*yt) iflag=1
      endif
      ! fails test.
      ! add midpoint to stack and continue.
      if (iflag.eq.1) then
         i=i+1
         x(i)=x(i-1)
         x(i-1)=xm
         y(i)=y(i-1)
         y(i-1)=yt
      ! passes test.
      ! use top point off the stack.
      else
         j=j+1
         if (j.gt.1) d=d+(x(i)-xlast)*(y(i)+ylast)/2
         xlast=x(i)
         ylast=y(i)
         i=i-1
      endif
   enddo
   return
   end subroutine anadam

   subroutine anadmn(e,d,z,awr,nu,a)
   !-------------------------------------------------------------------
   ! Damage energy for an analytic File 5 section of a reaction that
   ! emits nu neutrons and no charged particle.  Same as anadam,
   ! with the momenta of the nu neutrons added incoherently, residual
   ! mass ar=awr+1-nu, and the Lindhard constants set for the residual.
   ! the spectrum is integrated with Simpson's rule (nsim
   ! intervals) and normalized; the adaptive trapezoidal rule of anadam
   ! (5% tolerance) gave errors of up to 2%.  The same rule gives the
   ! moments <ep> and <pz> used for the File 4 anisotropy correction
   ! (as in tabdmn).
   ! all the analytic laws (anadam and the previous anadmn
   ! gave no damage except for LF=9):
   !   LF=5  general, g(ep/theta(e)) tabulated, ep <= e-u
   !   LF=7  Maxwell, sqrt(ep) exp(-ep/theta), ep <= e-u
   !   LF=9  evaporation, ep exp(-ep/theta), ep <= e-u
   !   LF=11 Watt, exp(-ep/a) sinh(sqrt(b ep)), ep <= e-u
   !   LF=12 Madland-Nix, (g(ep,efl)+g(ep,efh))/2 with
   !         g(ep,ef)=[u2**1.5 E1(u2) - u1**1.5 E1(u1) + gam(1.5,u2)
   !         - gam(1.5,u1)]/(3 sqrt(ef tm)), u1,2=(sqrt(ep)-+sqrt(ef))**2/tm
   ! The shapes need not be normalized (the integral is normalized).
   ! (lcvx, set by conbar except for fission): the recoil of
   ! the nu independent neutrons is computed exactly as in tabdmn, with
   ! nodes at the Simpson points (dirnod, cvcmp, cvdam). (DLA, IAEA, Oct 2026)
   !-------------------------------------------------------------------
   use endf ! provides terpa
   ! externals
   integer::nu
   real(kr)::e,d,z,awr,a(*)
   ! internals
   integer::nr,np,lf,lnext,ip,ir,idis,i,nrn,npn,lnxt,ns,ipg,irg
   real(kr)::u,awfac,test,theta,enext,ep,ar,dum,emx
   real(kr)::s0,s1,sz,w,p0,eri,era,bw,efl,efh,xmx,f
   ! exact path (DLA, IAEA, Oct 2026)
   integer::m,nsv(16),iov(16)
   real(kr)::cz(cvmx),ct(cvmx),cw(cvmx),erm
   integer,parameter::nsim=400
   integer,parameter::nsim5=2000
   real(kr),parameter::zero=0
   real(kr),parameter::delta=1.e-7_kr

   !--retrieve the spectrum parameters
   d=0
   u=a(1)
   nr=nint(a(5))
   np=nint(a(6))
   lf=nint(a(4))
   ar=awr+1-nu
   if (nu.le.1.or.ar.lt.0.5_kr) ar=awr
   awfac=1/ar
   dum=df(zero,z,ar,z,awr)
   if (lf.ne.5.and.lf.ne.7.and.lf.ne.9.and.lf.ne.11.and.lf.ne.12) return
   test=10
   if (lf.ne.12) then
      if (e.le.u) return
      if ((e-u).lt.delta*e) return
      if ((e-u).lt.test) return
   endif
   lnext=7+2*nr+2*np
   ip=2
   ir=1
   call terpa(theta,e,enext,idis,a(lnext),ip,ir)
   if (theta.le.zero) return
   ns=nsim
   bw=0
   efl=0
   efh=0
   lnxt=0
   if (lf.eq.5.or.lf.eq.11) then
      nrn=nint(a(lnext+4))
      npn=nint(a(lnext+5))
      lnxt=lnext+6+2*nrn+2*npn
   endif
   if (lf.eq.9.or.lf.eq.7) then
      emx=min(e-u,50*theta)
   else if (lf.eq.11) then
      ip=2
      ir=1
      call terpa(bw,e,enext,idis,a(lnxt),ip,ir)
      if (bw.lt.zero) bw=0
      emx=min(e-u,(sqrt(50*theta)+theta*sqrt(bw)/2)**2)
   else if (lf.eq.5) then
      ! last x of the g(x) table
      xmx=a(lnxt+6+2*nint(a(lnxt+4))+2*(nint(a(lnxt+5))-1))
      emx=min(e-u,theta*xmx)
      ns=nsim5
   else
      efl=a(lnext)
      efh=a(lnext+1)
      emx=(sqrt(max(efl,efh))+sqrt(50*theta))**2
   endif
   if (emx.le.zero) return

   !--Simpson integration over the secondary energy;
   !--Gauss-Legendre angle integration (tdfn)
   !--lcvx: laboratory momentum nodes at the Simpson
   !--points (dirnod) for the exact recoil of the nu neutrons (cvdam)
   ! (DLA, IAEA, Oct 2026)
   nfn=0
   s0=0
   s1=0
   sz=0
   ipg=2
   irg=1
   do i=0,ns
      ep=emx*i/ns
      w=2
      if (mod(i,2).eq.1) w=4
      if (i.eq.0.or.i.eq.ns) w=1
      if (lf.eq.9) then
         f=sed(e,ep,theta,u)
      else if (lf.eq.7) then
         f=sqrt(ep)*exp(-ep/theta)
      else if (lf.eq.11) then
         f=exp(-ep/theta)*sinh(sqrt(bw*ep))
      else if (lf.eq.5) then
         call terpa(f,ep/theta,enext,idis,a(lnxt),ipg,irg)
         if (f.lt.zero) f=0
      else
         f=(gmn(ep,efl,theta)+gmn(ep,efh,theta))/2
      endif
      w=w*f
      s0=s0+w
      if (lcvx) then
         call dirnod(e,ep,awr,w)
      else
         d=d+w*tdfn(e,ep,nu,z,ar,awr)
      endif
      s1=s1+w*ep
      if (n4e.gt.0.and..not.lcvx) sz=sz+w*anipz(e,ep,awr)
   enddo
   if (s0.le.zero) then
      d=0
      return
   endif
   ! exact recoil of the nu independent neutrons (tabdmn) (DLA, IAEA, Oct 2026)
   if (lcvx) then
      d=0
      call cvcmp(nfn,fnz,fnt,fnw,m,cz,ct,cw)
      if (m.gt.0) then
         do i=1,max(nu,1)
            nsv(i)=m
            iov(i)=0
         enddo
         call cvdam(max(nu,1),nsv,iov,cz,ct,cw,sqrt(2*e),z,ar,z,awr,d,erm)
      endif
      return
   endif
   d=d/s0
   ! anisotropy of the emitted neutrons (File 4), as in tabdmn
   if (n4e.gt.0.and.d.gt.zero) then
      s1=s1/s0
      sz=sz/s0
      p0=sqrt(2*e)
      eri=(e+nu*s1)*awfac
      era=(p0*p0-2*p0*nu*sz+2*nu*s1+nu*(nu-1)*sz*sz)*awfac/2
      if (era.gt.zero.and.eri.gt.zero) then
         d=d*df(era,z,ar,z,awr)/df(eri,z,ar,z,awr)
      endif
   endif
   return

   contains

      real(kr) function gmn(ep,ef,tm)
      ! Madland-Nix spectrum of one fragment (ENDF-6, File 5 LF=12)
      real(kr)::ep,ef,tm
      real(kr)::u1,u2
      gmn=0
      if (ef.le.zero.or.tm.le.zero.or.ep.lt.zero) return
      u1=(sqrt(ep)-sqrt(ef))**2/tm
      u2=(sqrt(ep)+sqrt(ef))**2/tm
      gmn=(u2**1.5_kr*e1x(u2)-u1**1.5_kr*e1x(u1)+gam15(u2)-gam15(u1))&
        /(3*sqrt(ef*tm))
      if (gmn.lt.zero) gmn=0
      end function gmn

      real(kr) function gam15(x)
      ! lower incomplete gamma function gamma(3/2,x)
      real(kr)::x
      real(kr),parameter::rpi2=0.886226925452758e0_kr
      gam15=rpi2*erf(sqrt(x))-sqrt(x)*exp(-x)
      end function gam15

      real(kr) function e1x(x)
      ! exponential integral E1(x), x>0 (0 for x=0, where x**1.5*E1
      ! vanishes); series for x<=1, continued fraction for x>1
      real(kr)::x
      integer::k
      real(kr)::s,t,b,c,dd,h,an
      real(kr),parameter::euler=0.5772156649015329e0_kr
      real(kr),parameter::tiny=1.e-30_kr
      real(kr),parameter::eps=1.e-12_kr
      e1x=0
      if (x.le.zero) return
      if (x.le.1) then
         s=0
         t=1
         do k=1,100
            t=-t*x/k
            s=s+t/k
            if (abs(t/k).lt.eps*abs(s)) exit
         enddo
         e1x=-euler-log(x)-s
      else
         b=x+1
         c=1/tiny
         dd=1/b
         h=dd
         do k=1,200
            an=-real(k,kr)*k
            b=b+2
            dd=1/(an*dd+b)
            c=b+an/c
            h=h*c*dd
            if (abs(c*dd-1).lt.eps) exit
         enddo
         e1x=h*exp(-x)
      endif
      end function e1x

   end subroutine anadmn

   real(kr) function sed(e,ep,theta,u)
   !-------------------------------------------------------------------
   ! Secondary energy distribution function for anadam.
   ! A Taylor expansion is used for small values to avoid overflow.
   !-------------------------------------------------------------------
   ! externals
   real(kr)::e,ep,theta,u
   ! internals
   real(kr)::xeu,aeu,test

   xeu=(e-u)/theta
   aeu=abs(xeu)
   test=1
   test=test/10000
   if (aeu.ge.test) then
      sed=ep*exp(-ep/theta)/(theta*theta*(1-exp(-xeu)*(1+xeu)))
   else
      sed=ep*exp(-ep/theta)/(theta*theta*xeu*xeu*(1-xeu)/2)
   endif
   return
   end function sed

  subroutine tabbar(f,a,law)
   !-------------------------------------------------------------------
   ! Compute mean energy for a tabulated section of File 5.
   ! (If law<0, the mean energy of a subsection of File 6
   ! can be computed using abs(law) as the interpolation law).
   !-------------------------------------------------------------------
   use util ! provides error
   ! externals
   integer::law
   real(kr)::f,a(*)
   ! internals
   integer::nr,np,ibase,ir,nbt,inn,ncyc,i
   real(kr)::xh,xhsq,yh,xl,xlsq,delxsq,yl,b
   real(kr)::alxl,alxh,alyl,alyh

   !--initialize.
   nbt=0
   if (law.gt.0.and.law.ne.1.and.law.ne.5)&
     call error('tabbar','coded for lf=1 and lf=5 only.',' ')
   if (law.ge.0) then
      f=0
      nr=nint(a(5))
      np=nint(a(6))
      ibase=6+2*nr
      ir=1
      nbt=nint(a(7))
      inn=nint(a(8))
      ncyc=2
   else
      f=0
      np=nint(a(6))
      inn=iabs(law)
      ncyc=nint(a(4))+2
      ibase=6
   endif
   xh=a(ibase+1)
   xhsq=xh*xh
   yh=a(ibase+2)

   !--accumulate contributions to integral.
   do i=2,np
      xl=xh
      xh=a(ibase+ncyc*(i-1)+1)
      xlsq=xhsq
      xhsq=xh*xh
      delxsq=xhsq-xlsq
      yl=yh
      yh=a(ibase+ncyc*(i-1)+2)
      if (law.gt.0.and.i.gt.nbt) then
         ir=ir+1
         nbt=nint(a(6+2*ir-1))
         inn=nint(a(6+2*ir))
      endif
      if (xl.ne.xh) then

         !--y is a constant.
         if (inn.eq.1) then
            f=f+yl*delxsq/2

         !--y is linear in x.
         else if (inn.eq.2) then
            b=(yh-yl)/(xh-xl)
            f=f+((yl-xl*b)*delxsq/2+b*(xh*xhsq-xl*xlsq)/3)

         !--y is linear in ln(x)
         else if (inn.eq.3) then
            alxl=log(xl)
            alxh=log(xh)
            b=(yh-yl)/(alxh-alxl)
            f=f+(yl-b*alxl+b/2)*delxsq/2+b*(xhsq*alxh-xlsq*alxl)/2

         !--ln(y) is linear in x
         else if (inn.eq.4) then
            alyl=log(yl)
            alyh=log(yh)
            b=(alyh-alyl)/(xh-xl)
            f=f+(yh*(b*xh-1)-yl*(b*xl-1))/b**2

         !--ln(y) is linear in ln(x)
         else if (inn.eq.5) then
            alxl=log(xl)
            alxh=log(xh)
            alyl=log(yl)
            alyh=log(yh)
            b=(alyh-alyl)/(alxh-alxl)
            f=f+yl*xlsq*((xh/xl)**(b+2)-1)/(b+2)
         endif
      endif
   enddo
   return
   end subroutine tabbar

   subroutine fis6tab(nin,matd,mtd,za,awr)
   !-------------------------------------------------------------------
   ! Fission given in File 6 without recoil: damage per
   ! fission at each incident energy of the first neutron subsection,
   ! with the same model as File 5 (conbar/tabdmn with one neutron):
   ! recoil of the fissioning nucleus (mass awr, charge z) from one
   ! neutron, emission isotropic in the frame of the data,
   !   lab (LCT=1, LAW=7): er=(e-2*sqrt(e*ep)*mu+ep)/awr
   !   CM  (LCT=2,3):      er=(awr**2*v**2-2*awr*v*u*mu+u**2)/(2*awr),
   !                       v=sqrt(2e)/(awr+1), u=sqrt(2ep)
   ! 4-point Gauss rule in mu and in each outgoing-energy panel (LEP law
   ! for LAW=1, the TAB1 laws for LAW=7, trapezoid in mu for LAW=7),
   ! normalized.  For LAW=7 the recoil is computed at the
   ! tabulated laboratory cosine (no isotropic average).  The other particles of the section (ternary p, t,
   ! alpha, photons) are not used.  Results in ef6, df6 (nf6 points).
   ! (DLA, IAEA, Oct 2026)
   !-------------------------------------------------------------------
   use util ! provides error
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,matd,mtd
   real(kr)::za,awr
   ! internals
   integer::nb,nw,nk,lct,ik,l,law,izp,ne,ie,nd,na,nep,ncyc,j,iq
   integer::nmu,imu,nr1,np1,jr,lw,lep,i,done
   real(kr)::z,en,x1,x2,fa,fb,fx,s0,s1,dum,amu,amul,t0,t1,t0l,t1l,y,er
   real(kr),allocatable::b(:)
   integer,parameter::nbmax=2000000
   integer,parameter::nemax=2000
   real(kr),parameter::zero=0
   real(kr),dimension(4),parameter::q4x=(/&
      -.861136312e0_kr,-.339981044e0_kr,.339981044e0_kr,.861136312e0_kr/)
   real(kr),dimension(4),parameter::q4w=(/&
      .347854845e0_kr,.652145155e0_kr,.652145155e0_kr,.347854845e0_kr/)

   nf6=0
   if (allocated(ef6)) deallocate(ef6)
   if (allocated(df6)) deallocate(df6)
   allocate(ef6(nemax),df6(nemax))
   allocate(b(nbmax))
   z=int(za/1000)
   dum=df(zero,z,awr,z,awr)
   call findf(matd,6,mtd,nin)
   call contio(nin,0,0,b,nb,nw)
   lct=l2h
   nk=n1h
   done=0
   do ik=1,nk
      if (done.eq.1) exit
      l=1
      call tab1io(nin,0,0,b(l),nb,nw)
      l=l+nw
      do while (nb.ne.0)
         call moreio(nin,0,0,b(l),nb,nw)
         l=l+nw
      enddo
      izp=nint(b(1))
      law=nint(b(4))
      if (law.eq.1.or.law.eq.2.or.law.eq.5) then
         call tab2io(nin,0,0,b(1),nb,nw)
         lep=nint(b(4))
         ne=nint(b(6))
         do ie=1,ne
            l=1
            call listio(nin,0,0,b(l),nb,nw)
            l=l+nw
            do while (nb.ne.0)
               if (l.gt.nbmax-10000) call error('fis6tab','storage',' ')
               call moreio(nin,0,0,b(l),nb,nw)
               l=l+nw
            enddo
            if (izp.eq.1.and.law.eq.1.and.nf6.lt.nemax) then
               en=b(2)
               nd=nint(b(3))
               na=nint(b(4))
               nep=nint(b(6))
               ncyc=na+2
               s0=0
               s1=0
               do j=1,min(nd,nep)
                  x1=b(6+ncyc*(j-1)+1)
                  fa=b(6+ncyc*(j-1)+2)
                  s0=s0+fa
                  s1=s1+fa*d1fis(en,x1,awr,z,lct)
               enddo
               do j=nd+2,nep
                  x1=b(6+ncyc*(j-2)+1)
                  x2=b(6+ncyc*(j-1)+1)
                  if (x2.le.x1) cycle
                  fa=b(6+ncyc*(j-2)+2)
                  fb=b(6+ncyc*(j-1)+2)
                  do iq=1,4
                     y=(x1+x2)/2+(x2-x1)/2*q4x(iq)
                     fx=fa
                     if (lep.ne.1) fx=fa+(fb-fa)*(y-x1)/(x2-x1)
                     s0=s0+q4w(iq)*(x2-x1)/2*fx
                     s1=s1+q4w(iq)*(x2-x1)/2*fx*d1fis(en,y,awr,z,lct)
                  enddo
               enddo
               nf6=nf6+1
               ef6(nf6)=en
               df6(nf6)=0
               if (s0.gt.zero) df6(nf6)=s1/s0
            endif
         enddo
         if (izp.eq.1.and.law.eq.1) done=1
      else if (law.eq.6) then
         call contio(nin,0,0,b,nb,nw)
      else if (law.eq.7) then
         call tab2io(nin,0,0,b(1),nb,nw)
         ne=nint(b(6))
         do ie=1,ne
            call tab2io(nin,0,0,b(1),nb,nw)
            en=b(2)
            nmu=nint(b(6))
            t0=0
            t1=0
            t0l=0
            t1l=0
            amul=0
            do imu=1,nmu
               l=1
               call tab1io(nin,0,0,b(l),nb,nw)
               l=l+nw
               do while (nb.ne.0)
                  if (l.gt.nbmax-10000) call error('fis6tab','storage',' ')
                  call moreio(nin,0,0,b(l),nb,nw)
                  l=l+nw
               enddo
               if (izp.ne.1) cycle
               amu=b(2)
               nr1=nint(b(5))
               np1=nint(b(6))
               s0=0
               s1=0
               do j=2,np1
                  x1=b(6+2*nr1+2*(j-1)-1)
                  x2=b(6+2*nr1+2*j-1)
                  if (x2.le.x1) cycle
                  fa=b(6+2*nr1+2*(j-1))
                  fb=b(6+2*nr1+2*j)
                  lw=nint(b(6+2*nr1))
                  do jr=nr1,1,-1
                     if (j.le.nint(b(5+2*jr))) lw=nint(b(6+2*jr))
                  enddo
                  do iq=1,4
                     y=(x1+x2)/2+(x2-x1)/2*q4x(iq)
                     call terp1(x1,fa,x2,fb,y,fx,lw)
                     s0=s0+q4w(iq)*(x2-x1)/2*fx
                     ! recoil at the laboratory cosine amu of
                     ! this table (the isotropic average d1fis was used)
                     ! (DLA, IAEA, Oct 2026)
                     er=(en-2*sqrt(en*max(y,zero))*amu+y)/awr
                     s1=s1+q4w(iq)*(x2-x1)/2*fx*df(max(er,zero),z,awr,z,awr)
                  enddo
               enddo
               if (imu.gt.1) then
                  t0=t0+(amu-amul)*(s0+t0l)/2
                  t1=t1+(amu-amul)*(s1+t1l)/2
               endif
               amul=amu
               t0l=s0
               t1l=s1
            enddo
            if (izp.eq.1.and.nf6.lt.nemax) then
               nf6=nf6+1
               ef6(nf6)=en
               df6(nf6)=0
               if (t0.gt.zero) df6(nf6)=t1/t0
            endif
         enddo
         if (izp.eq.1) done=1
      endif
   enddo
   deallocate(b)
   return
   end subroutine fis6tab

   real(kr) function d1fis(e,ep,awr,z,lct)
   !-------------------------------------------------------------------
   ! Damage of the fissioning nucleus (mass awr, charge z) recoiling
   ! from one neutron of energy ep, isotropic in the frame lct of the
   ! data (1 laboratory, otherwise centre of mass), 4-point Gauss rule.
   ! The Lindhard constants must be set with df(0,z,awr,z,awr).
   !-------------------------------------------------------------------
   real(kr)::e,ep,awr,z
   integer::lct
   integer::iq
   real(kr)::er,v,u
   real(kr),dimension(4),parameter::q4x=(/&
      -.861136312e0_kr,-.339981044e0_kr,.339981044e0_kr,.861136312e0_kr/)
   real(kr),dimension(4),parameter::q4w=(/&
      .347854845e0_kr,.652145155e0_kr,.652145155e0_kr,.347854845e0_kr/)
   real(kr),parameter::zero=0
   d1fis=0
   v=sqrt(2*e)/(awr+1)
   u=sqrt(2*max(ep,zero))
   do iq=1,4
      if (lct.eq.1) then
         er=(e-2*sqrt(e*max(ep,zero))*q4x(iq)+ep)/awr
      else
         er=(awr*awr*v*v-2*awr*v*u*q4x(iq)+u*u)/(2*awr)
      endif
      d1fis=d1fis+q4w(iq)*df(max(er,zero),z,awr,z,awr)/2
   enddo
   return
   end function d1fis

   real(kr) function fis6dm(e)
   !-------------------------------------------------------------------
   ! Damage per fission at incident energy e from the table of fis6tab
   ! (linear interpolation; constant outside the table).
   !-------------------------------------------------------------------
   real(kr)::e
   integer::i
   fis6dm=0
   if (nf6.le.0) return
   if (e.le.ef6(1)) then
      fis6dm=df6(1)
   else if (e.ge.ef6(nf6)) then
      fis6dm=df6(nf6)
   else
      do i=2,nf6
         if (e.le.ef6(i)) then
            if (ef6(i).gt.ef6(i-1)) then
               fis6dm=df6(i-1)+(df6(i)-df6(i-1))*(e-ef6(i-1))/(ef6(i)-ef6(i-1))
            else
               fis6dm=df6(i)
            endif
            exit
         endif
      enddo
   endif
   return
   end function fis6dm

   subroutine tabdam(e,d,z,awr,awfac,a)
   !-------------------------------------------------------------------
   ! Compute damage energy for a tabulated section of File 5.
   !-------------------------------------------------------------------
   ! externals
   real(kr)::e,d,z,awr,awfac,a(*)
   ! internals
   integer::nr,np,ibase,i,iq
   real(kr)::xh,yh,fl,xl,f,er
   integer,parameter::nq=4
   real(kr),dimension(nq),parameter::qp=(/&
     -.86114e0_kr,-.33998e0_kr,.33998e0_kr,.86114e0_kr/)
   real(kr),dimension(nq),parameter::qw=(/&
     .34785e0_kr,.65215e0_kr,.65215e0_kr,.34785e0_kr/)
   real(kr),parameter::zero=0

   !--initialize
   if (e.eq.zero) then
      d=df(e,z,awr,z,awr)
      return
   endif

   !--trapazoidal energy integration
   !--gauss-legendre angle integration
   d=0
   nr=nint(a(5))
   np=nint(a(6))
   ibase=6+2*nr
   xh=a(ibase+1)
   yh=a(ibase+2)
   fl=0
   do i=1,np
      xl=xh
      xh=a(ibase+2*i-1)
      yh=a(ibase+2*i)
      f=0
      do iq=1,nq
         er=(e-2*sqrt(e*xh)*qp(iq)+xh)*awfac
         f=f+qw(iq)*df(er,z,awr,z,awr)/2
      enddo
      f=f*yh
      if (i.gt.1) d=d+(xh-xl)*(f+fl)/2
      fl=f
   enddo
   return
   end subroutine tabdam

   subroutine tabdmn(e,d,z,awr,nu,a)
   !-------------------------------------------------------------------
   ! Damage energy for a tabulated File 5 spectrum of a reaction that
   ! emits nu neutrons and no charged particle ((n,n'), (n,2n), ...).
   ! Corrections to tabdam:
   !  - the momenta of the nu neutrons are added incoherently:
   !    er=(e-2*sqrt(e*nu*ep)*mu+nu*ep)/ar, ar=awr+1-nu, whose mean
   !    is the exact isotropic value (e+nu*<ep>)/ar (tabdam used one
   !    neutron and the target mass awr);
   !  - the Lindhard constants are set for the residual (z,ar)
   !    (tabdam used whatever constants df had saved last);
   !  - the spectrum is integrated with its interpolation law (histogram
   !    or linear) and the result is normalized (tabdam used the
   !    trapezoidal rule without normalization); panels with the log
   !    laws 3-5 are divided into 16 sub-panels;
   !  - anisotropy of the emitted neutrons from File 4
   !    (anitab, anipz), see below; centre-of-mass (LCT=2)
   !    distributions converted to the laboratory with the Jacobian at
   !    fixed laboratory energy (anipz).
   ! (lcvx, all neutron-only reactions except fission): the
   ! recoil of the nu independent neutrons is computed exactly with the
   ! File 4 x File 5 distribution of one neutron (dirnod, cvcmp, cvdam)
   ! instead of the isotropic one-vector form of tdfn (|sum p| fixed at
   ! sqrt(nu) p) and the mean-energy File 4 factor; fission keeps the
   ! integration above with one neutron (nu=1, no File 4). (DLA, IAEA, Oct 2026)
   !-------------------------------------------------------------------
   use endf ! provides terp1
   ! externals
   integer::nu
   real(kr)::e,d,z,awr,a(*)
   ! internals
   integer::nr,np,ibase,i,ir,nbt,inn
   real(kr)::ar,xl,yl,fl,xh,yh,fh,sn,dum
   real(kr)::s1,sz,zl,zh,eri,era,p0
   real(kr)::xa,ya,fa,za,xb,yb,fb,zb
   integer::j,lw
   integer,parameter::nsub=16
   real(kr),parameter::zero=0
   ! exact path (DLA, IAEA, Oct 2026)
   integer::m,iq,nsb,nsv(16),iov(16)
   real(kr)::cz(cvmx),ct(cvmx),cw(cvmx),erm,x,y
   real(kr),dimension(4),parameter::q4x=(/&
      -.861136312e0_kr,-.339981044e0_kr,.339981044e0_kr,.861136312e0_kr/)
   real(kr),dimension(4),parameter::q4w=(/&
      .347854845e0_kr,.652145155e0_kr,.652145155e0_kr,.347854845e0_kr/)
   d=0
   s1=0
   sz=0
   zl=0
   zh=0
   ar=awr+1-nu
   if (nu.le.1.or.ar.lt.0.5_kr) ar=awr
   dum=df(zero,z,ar,z,awr)
   nr=nint(a(5))
   np=nint(a(6))
   ibase=6+2*nr
   ir=1
   nbt=nint(a(7))
   inn=nint(a(8))
   ! exact recoil of the nu independently emitted neutrons
   ! (lcvx, set by conbar for the neutron-only reactions other than
   ! fission): laboratory momentum nodes of one neutron (4-point Gauss
   ! rule in each panel of the spectrum with its interpolation law, the
   ! log laws on nsub linear sub-panels as below; directions from File
   ! 4 by dirnod), compressed by cvcmp and added nu times by cvdam.
   ! This replaces the isotropic one-vector form of tdfn and the File 4
   ! factor below, which are then not used. (DLA, IAEA, Oct 2026)
   if (lcvx.and.nu.le.16) then
      nfn=0
      do i=2,np
         xl=a(ibase+2*i-3)
         yl=a(ibase+2*i-2)
         xh=a(ibase+2*i-1)
         yh=a(ibase+2*i)
         do while (i.gt.nbt.and.ir.lt.nr)
            ir=ir+1
            nbt=nint(a(5+2*ir))
            inn=nint(a(6+2*ir))
         enddo
         if (xh.le.xl) cycle
         lw=inn
         nsb=1
         if (inn.ge.3) then
            nsb=nsub
            if ((lw.eq.4.or.lw.eq.5).and.(yl.le.zero.or.yh.le.zero)) lw=2
            if ((lw.eq.3.or.lw.eq.5).and.xl.le.zero) lw=2
         endif
         xa=xl
         ya=yl
         do j=1,nsb
            xb=xl+(xh-xl)*j/nsb
            if (inn.eq.1) then
               yb=yl
            else if (j.eq.nsb) then
               yb=yh
            else
               call terp1(xl,yl,xh,yh,xb,yb,lw)
            endif
            do iq=1,4
               x=(xa+xb)/2+(xb-xa)/2*q4x(iq)
               y=ya+(yb-ya)*(x-xa)/(xb-xa)
               call dirnod(e,x,awr,q4w(iq)*(xb-xa)/2*y)
            enddo
            xa=xb
            ya=yb
         enddo
      enddo
      call cvcmp(nfn,fnz,fnt,fnw,m,cz,ct,cw)
      if (m.gt.0) then
         do i=1,max(nu,1)
            nsv(i)=m
            iov(i)=0
         enddo
         call cvdam(max(nu,1),nsv,iov,cz,ct,cw,sqrt(2*e),z,ar,z,awr,d,erm)
      endif
      return
   endif
   xl=a(ibase+1)
   yl=a(ibase+2)
   fl=tdfn(e,xl,nu,z,ar,awr)
   if (n4e.gt.0) zl=anipz(e,xl,awr)
   sn=0
   do i=2,np
      xh=a(ibase+2*i-1)
      yh=a(ibase+2*i)
      fh=tdfn(e,xh,nu,z,ar,awr)
      if (n4e.gt.0) zh=anipz(e,xh,awr)
      do while (i.gt.nbt.and.ir.lt.nr)
         ir=ir+1
         nbt=nint(a(5+2*ir))
         inn=nint(a(6+2*ir))
      enddo
      if (inn.eq.1) then
         d=d+yl*(fl+fh)*(xh-xl)/2
         sn=sn+yl*(xh-xl)
         s1=s1+yl*(xl+xh)*(xh-xl)/2
         sz=sz+yl*(zl+zh)*(xh-xl)/2
      else if (inn.eq.2.or.xh.le.xl) then
         d=d+(yl*fl+yh*fh)*(xh-xl)/2
         sn=sn+(yl+yh)*(xh-xl)/2
         s1=s1+(yl*xl+yh*xh)*(xh-xl)/2
         sz=sz+(yl*zl+yh*zh)*(xh-xl)/2
      else
         ! log interpolation laws (3,4,5): the panel is divided into
         ! nsub sub-panels with values from the interpolation law
         lw=inn
         if ((lw.eq.4.or.lw.eq.5).and.(yl.le.zero.or.yh.le.zero)) lw=2
         if ((lw.eq.3.or.lw.eq.5).and.xl.le.zero) lw=2
         xa=xl
         ya=yl
         fa=fl
         za=zl
         do j=1,nsub
            xb=xl+(xh-xl)*j/nsub
            if (j.eq.nsub) then
               yb=yh
               fb=fh
               zb=zh
            else
               call terp1(xl,yl,xh,yh,xb,yb,lw)
               fb=tdfn(e,xb,nu,z,ar,awr)
               zb=0
               if (n4e.gt.0) zb=anipz(e,xb,awr)
            endif
            d=d+(ya*fa+yb*fb)*(xb-xa)/2
            sn=sn+(ya+yb)*(xb-xa)/2
            s1=s1+(ya*xa+yb*xb)*(xb-xa)/2
            sz=sz+(ya*za+yb*zb)*(xb-xa)/2
            xa=xb
            ya=yb
            fa=fb
            za=zb
         enddo
      endif
      xl=xh
      yl=yh
      fl=fh
      zl=zh
   enddo
   if (sn.gt.zero) then
      d=d/sn
   else
      d=0
   endif
   ! anisotropy of the emitted neutrons (File 4): the damage is scaled
   ! with the exact mean recoil energy of nu independent neutrons,
   !   <er> = [p0**2 - 2 p0 nu <pz> + 2 nu <ep> + nu (nu-1) <pz>**2]/(2 ar)
   ! relative to the isotropic mean (e + nu <ep>)/ar used above
   ! (<pz>: mean laboratory momentum along the beam per neutron).
   if (n4e.gt.0.and.d.gt.zero) then
      s1=s1/sn
      sz=sz/sn
      p0=sqrt(2*e)
      eri=(e+nu*s1)/ar
      era=(p0*p0-2*p0*nu*sz+2*nu*s1+nu*(nu-1)*sz*sz)/(2*ar)
      if (era.gt.zero.and.eri.gt.zero) then
         d=d*df(era,z,ar,z,awr)/df(eri,z,ar,z,awr)
      endif
   endif
   return
   end subroutine tabdmn

   real(kr) function tdfn(e,ep,nu,z,ar,awr)
   !-------------------------------------------------------------------
   ! Damage energy for incident energy e and one point ep of the
   ! neutron spectrum, averaged over mu with a 4-point Gauss rule,
   ! for nu neutrons and residual mass ar (df constants already set).
   !-------------------------------------------------------------------
   ! externals
   integer::nu
   real(kr)::e,ep,z,ar,awr
   ! internals
   integer::iq
   real(kr)::er
   integer,parameter::nq=4
   real(kr),dimension(nq),parameter::qp=(/&
     -.86114e0_kr,-.33998e0_kr,.33998e0_kr,.86114e0_kr/)
   real(kr),dimension(nq),parameter::qw=(/&
     .34785e0_kr,.65215e0_kr,.65215e0_kr,.34785e0_kr/)
   tdfn=0
   do iq=1,nq
      er=(e-2*sqrt(e*nu*ep)*qp(iq)+nu*ep)/ar
      if (er.lt.0) er=0
      tdfn=tdfn+qw(iq)*df(er,z,ar,z,awr)/2
   enddo
   return
   end function tdfn

   subroutine oprtab(nin,matd,mtd,za,awr,lct,nk,idrv)
   !-------------------------------------------------------------------
   ! Correction of the one-particle recoil (OPR) of sixbar for File 6
   ! reactions without residual that emit several particles: (n,2n),
   ! (n,np), (n,na), ...
   ! sixbar computes the recoil from one particle only (subsection
   ! idrv) and the residual mass awr+1-awp.  Here the exact mean
   ! recoil energy of independent emission of all particles k, with
   ! yields y_k, is computed in the laboratory from the moments
   ! <p mu>_k and <p**2>_k of their File 6 distributions (LAW=1):
   !   <er> = [p0**2 - 2 p0 S1 + S2 + S1**2 - Sq]/(2 Ar)
   !   S1=sum y_k <p mu>_k, Sq=sum y_k <p mu>_k**2, S2=sum y_k <p**2>_k
   !   Ar=awr+1-sum y_k awp_k, Zr=Z-sum y_k z_k, p0=sqrt(2 e)
   ! and the OPR damage is multiplied by
   !   r(e) = df(<er>;Zr,Ar)/df(<er>_1;Z1,A1)
   ! with <er>_1 the same expression for the particle idrv alone.
   ! The moments are integrated in the frame of the data (no Jacobian
   ! is needed) and transformed to the laboratory exactly:
   !   <p mu>_lab = <p mu>_cm + m V,  <p**2>_lab = <p**2>_cm + 2 m V
   !   <p mu>_cm + m**2 V**2,  V = sqrt(2 e)/(awr+1),
   ! for the frames used by sixbar (centre of mass unless LCT=1, or
   ! LCT=3 for particles heavier than alpha).
   ! 4-point Gauss rule in each outgoing-energy panel with
   ! the LEP law for the distribution and its angular parameters
   ! (Legendre f1, Kalbach r with a from the systematics with incident
   ! neutron, as in sixbar; tabulated angular distributions LANG=11-15
   ! included); yields and moments are interpolated in incident energy
   ! with the laws of the data (unit-base laws 21-25 as 1-5).
   ! The table (ert,rrt) is built on the incident energies of the
   ! subsection idrv; nrt=0 (no correction) if a particle is not
   ! given with LAW=1.  The tape is repositioned after the HEAD record.
   ! exact recoil of the independently emitted particles.
   ! The factor above corrects only the mean recoil energy; since D is
   ! concave it overestimated the damage (+0.1 to +2% per reaction,
   ! e.g. JENDL-5 Fe-56 (n,2n) +1.2% at 14 MeV, +2.0% at 20 MeV).  Now
   ! the laboratory momentum distribution of each particle at each of
   ! its incident energies is represented by nodes (pz, pt, weight):
   ! LAW=1 (l1nod: discrete lines and 4-point Gauss points of each
   ! panel with the LEP law, 64 cosines or the panels of a tabulated
   ! angular distribution, frame as above), LAW=7 (trapezoidal rule in
   ! mu as getsix, 4-point Gauss in e'), LAW=6 (phase space at e),
   ! compressed by cvcmp (moments kept).  At each energy of the table
   ! (getnod: stored table, or the two tables around e mixed with the
   ! weight of the incident-energy law) the mean damage of the residual
   ! (Zr, Ar) after all emissions, dN, and the one-particle damage of
   ! the driving particle, d1, are computed with the same nodes
   ! (cvdam) and rrt=dN/d1; the one-particle damage of sixbar times rrt
   ! is the exact damage for independent emission (the quadrature of
   ! sixbar and of the nodes cancels to first order).  The mean-energy
   ! factor is kept, with a message, at energies with non-integer
   ! yields or more than nvmax particles. (DLA, IAEA, Oct 2026)
   !-------------------------------------------------------------------
   use endf
   use util ! provides error,mess
   use mathm ! provides legndr
   ! externals
   integer::nin,matd,mtd,lct,nk,idrv
   real(kr)::za,awr
   ! internals
   integer,parameter::npm=12,nem=400,nym=500,nrm=20
   integer::np,ip,ik,nb,nw,l,law,lang,lep,ne,ie,nd,na,nep,ncyc,i,j
   integer::izap,izat,nr,ny(npm),nein(npm),ipd,ok,nmu,imu,k,iq,ir,lw
   integer::nry(npm),nby(npm,nrm),iny(npm,nrm),lwe(npm,nem)
   integer::nr2,nbt2(nrm),int2(nrm),lfr
   real(kr)::zp(npm),ap(npm),yx(npm,nym),yy(npm,nym)
   real(kr)::ein(npm,nem),pmu(npm,nem),p2(npm,nem)
   real(kr)::en,ep,f,x,s0,s1,sm,mub,aa,vcm,pm,emean,am
   real(kr)::e,p0,y,ssum1,sq,s2,ar,zr,er,er1,a1,z1,d1,dn,zt,dum,hd(6)
   real(kr)::x1,x2,fa,fb,ga,gb,g,fx
   ! LAW=6 (phase space) and LAW=7 (laboratory) particles (DLA, IAEA, Oct 2026)
   logical::l6(npm)
   real(kr)::aps6(npm),u,ul,f0,f1,fp,f0l,f1l,fpl
   integer::nps6(npm),nle,jl
   real(kr)::elist6(nem)
   real(kr),allocatable::b(:)
   ! exact recoil of the independently emitted particles (DLA, IAEA, Oct 2026)
   integer,parameter::nvmax=8
   integer::nodo(npm,nem),nodn(npm,nem),nnod,nmult(npm),nv,nsv(nvmax)
   integer::iov(nvmax),mn,nwk,io1,n1,nfbk,n7c,n7p,nexa,jj
   logical::lex
   real(kr),dimension(:),allocatable::ndz,ndt,ndw,r7x,r7w,p7x,p7w
   real(kr)::wkz(npm*cvmx),wkt(npm*cvmx),wkw(npm*cvmx)
   real(kr)::cnz(cvmx),cnt(cvmx),cnw(cvmx),pleg(65,64)
   real(kr)::dex,d1x,erm,arx,zrx
   character(60)::strng
   integer,parameter::nbmax=400000
   real(kr),parameter::zero=0
   real(kr),dimension(4),parameter::q4x=(/&
      -.861136312e0_kr,-.339981044e0_kr,.339981044e0_kr,.861136312e0_kr/)
   real(kr),dimension(4),parameter::q4w=(/&
      .347854845e0_kr,.652145155e0_kr,.652145155e0_kr,.347854845e0_kr/)
   nrt=0
   allocate(b(nbmax))
   izat=nint(za)
   zt=int(za/1000)
   np=0
   ipd=0
   ok=1
   nnod=0
   nexa=0
   nfbk=0
   n7c=0
   n7p=0
   nodn=0
   nodo=0
   allocate(ndz(20000),ndt(20000),ndw(20000))
   allocate(r7x(1000),r7w(1000),p7x(1000),p7w(1000))
   do iq=1,64
      call legndr(g64x(iq),pleg(1,iq),64)
   enddo
   do ik=1,nk
      l=1
      call tab1io(nin,0,0,b(l),nb,nw)
      l=l+nw
      do while (nb.ne.0)
         call moreio(nin,0,0,b(l),nb,nw)
         l=l+nw
      enddo
      izap=nint(b(1))
      law=nint(b(4))
      nr=nint(b(5))
      i=0
      if (izap.gt.0.and.izap.le.2004.and.np.lt.npm) then
         np=np+1
         i=np
         if (ik.eq.idrv) ipd=np
         zp(np)=int(izap/1000)
         ap(np)=b(2)
         ny(np)=min(nint(b(6)),nym)
         nry(np)=min(nr,nrm)
         do j=1,nry(np)
            nby(np,j)=nint(b(5+2*j))
            iny(np,j)=nint(b(6+2*j))
         enddo
         do j=1,ny(np)
            yx(np,j)=b(6+2*nr+2*j-1)
            yy(np,j)=b(6+2*nr+2*j)
         enddo
         nein(np)=0
         l6(np)=.false.
         ! LAW=6 and LAW=7 are also used (moments below) (DLA, IAEA, Oct 2026)
         if (law.ne.1.and.law.ne.6.and.law.ne.7) ok=0
      endif
      if (law.eq.1.or.law.eq.2.or.law.eq.5) then
         call tab2io(nin,0,0,b(1),nb,nw)
         lang=nint(b(3))
         lep=nint(b(4))
         ne=nint(b(6))
         nr2=min(nint(b(5)),nrm)
         do j=1,nr2
            nbt2(j)=nint(b(5+2*j))
            int2(j)=nint(b(6+2*j))
         enddo
         do ie=1,ne
            l=1
            call listio(nin,0,0,b(l),nb,nw)
            l=l+nw
            do while (nb.ne.0)
               if (l.gt.nbmax-10000) call error('oprtab','storage',' ')
               call moreio(nin,0,0,b(l),nb,nw)
               l=l+nw
            enddo
            if (i.gt.0.and.law.eq.1.and.nein(max(i,1)).lt.nem) then
               en=b(2)
               nd=nint(b(3))
               na=nint(b(4))
               nep=nint(b(6))
               ncyc=na+2
               am=ap(i)
               s0=0
               s1=0
               sm=0
               ! discrete lines (delta functions)
               do j=1,min(nd,nep)
                  ep=b(6+ncyc*(j-1)+1)
                  f=b(6+ncyc*(j-1)+2)
                  s0=s0+f
                  s1=s1+f*ep
                  sm=sm+gmu(f,angpar(j),ep)*sqrt(2*am*ep)
               enddo
               ! continuum panels, 4-point Gauss rule with the LEP law
               do j=nd+2,nep
                  x1=b(6+ncyc*(j-2)+1)
                  x2=b(6+ncyc*(j-1)+1)
                  if (x2.le.x1) cycle
                  fa=b(6+ncyc*(j-2)+2)
                  fb=b(6+ncyc*(j-1)+2)
                  do iq=1,4
                     x=(x1+x2)/2+(x2-x1)/2*q4x(iq)
                     ! distribution and angular parameter interpolated
                     ! separately with the LEP law, as in h6ddx
                     ga=angpar(j-1)
                     if (lep.eq.1) then
                        fx=fa
                        g=ga
                     else
                        gb=angpar(j)
                        fx=fa+(fb-fa)*(x-x1)/(x2-x1)
                        g=ga+(gb-ga)*(x-x1)/(x2-x1)
                     endif
                     g=gmu(fx,g,x)
                     s0=s0+q4w(iq)*(x2-x1)/2*fx
                     s1=s1+q4w(iq)*(x2-x1)/2*fx*x
                     sm=sm+q4w(iq)*(x2-x1)/2*g*sqrt(2*am*x)
                  enddo
               enddo
               if (s0.gt.zero) then
                  nein(i)=nein(i)+1
                  k=nein(i)
                  ein(i,k)=en
                  ! momentum nodes of this table (DLA, IAEA, Oct 2026)
                  call l1nod(i)
                  call stonod(i,k)
                  ! incident-energy interpolation law to the next point
                  lw=int2(nr2)
                  do j=nr2,1,-1
                     if (ie+1.le.nbt2(j)) lw=int2(j)
                  enddo
                  if (lw.gt.20) lw=lw-20
                  lwe(i,k)=lw
                  emean=s1/s0
                  pmu(i,k)=sm/s0
                  p2(i,k)=2*am*emean
                  ! frame of the data as in sixbar: laboratory for LCT=1
                  ! and for heavy particles with LCT=3; otherwise CM
                  lfr=2
                  if (lct.eq.1.or.(lct.eq.3.and.izap.gt.2004)) lfr=1
                  if (lfr.eq.2) then
                     vcm=sqrt(2*en)/(awr+1)
                     p2(i,k)=p2(i,k)+2*am*vcm*pmu(i,k)+am*am*vcm*vcm
                     pmu(i,k)=pmu(i,k)+am*vcm
                  endif
               endif
            endif
         enddo
      else if (law.eq.6) then
         ! n-body phase space, isotropic in the cm system;
         ! the laboratory moments are analytic (see mint):
         !   <e'>cm = emax/(npsx-1), emax=(apsx-a)/apsx*(awr e/(awr+1)+q)
         !   <p mu>lab = a V, <p**2>lab = 2 a (<e'>cm + a e/(awr+1)**2)
         ! (DLA, IAEA, Oct 2026)
         call contio(nin,0,0,b(1),nb,nw)
         if (i.gt.0) then
            l6(i)=.true.
            aps6(i)=b(1)
            nps6(i)=nint(b(6))
         endif
      else if (law.eq.7) then
         ! laboratory angle-energy law; for each incident
         ! energy the moments of f(mu,e') are integrated over e' (each
         ! TAB1 with its law, 4-point Gauss for linear panels) and over
         ! mu with the trapezoidal rule (as getsix); no frame change
         ! (DLA, IAEA, Oct 2026)
         call tab2io(nin,0,0,b(1),nb,nw)
         ne=nint(b(6))
         nr2=min(nint(b(5)),nrm)
         do j=1,nr2
            nbt2(j)=nint(b(5+2*j))
            int2(j)=nint(b(6+2*j))
         enddo
         do ie=1,ne
            call tab2io(nin,0,0,b(1),nb,nw)
            en=b(2)
            nmu=nint(b(6))
            ! momentum nodes (trapezoidal rule in mu) (DLA, IAEA, Oct 2026)
            nfn=0
            n7p=0
            s0=0
            s1=0
            sm=0
            ul=0
            f0l=0
            f1l=0
            fpl=0
            do imu=1,nmu
               l=1
               call tab1io(nin,0,0,b(l),nb,nw)
               l=l+nw
               do while (nb.ne.0)
                  if (l.gt.nbmax-10000) call error('oprtab','storage',' ')
                  call moreio(nin,0,0,b(l),nb,nw)
                  l=l+nw
               enddo
               if (i.gt.0) then
                  u=b(2)
                  call e7mom(ap(i),f0,f1,fp)
                  ! nodes of this mu and of the previous one
                  ! with half of the mu interval each (DLA, IAEA, Oct 2026)
                  call r7row(n7c)
                  if (imu.gt.1) then
                     do jj=1,n7p
                        call addnod(p7x(jj),ul,(u-ul)/2*p7w(jj),ap(i),zero)
                     enddo
                     do jj=1,n7c
                        call addnod(r7x(jj),u,(u-ul)/2*r7w(jj),ap(i),zero)
                     enddo
                  endif
                  if (size(p7x).lt.n7c) then
                     deallocate(p7x,p7w)
                     allocate(p7x(n7c),p7w(n7c))
                  endif
                  if (n7c.gt.0) then
                     p7x(1:n7c)=r7x(1:n7c)
                     p7w(1:n7c)=r7w(1:n7c)
                  endif
                  n7p=n7c
                  if (imu.gt.1) then
                     s0=s0+(u-ul)*(f0+f0l)/2
                     s1=s1+(u-ul)*(f1+f1l)/2
                     sm=sm+(u-ul)*(u*fp+ul*fpl)/2
                  endif
                  ul=u
                  f0l=f0
                  f1l=f1
                  fpl=fp
               endif
            enddo
            if (i.gt.0.and.s0.gt.zero.and.nein(max(i,1)).lt.nem) then
               nein(i)=nein(i)+1
               k=nein(i)
               ein(i,k)=en
               lw=int2(nr2)
               do j=nr2,1,-1
                  if (ie+1.le.nbt2(j)) lw=int2(j)
               enddo
               if (lw.gt.20) lw=lw-20
               lwe(i,k)=lw
               pmu(i,k)=sm/s0
               p2(i,k)=2*ap(i)*s1/s0
               ! momentum nodes of this table (DLA, IAEA, Oct 2026)
               call stonod(i,k)
            endif
         enddo
      else if (law.lt.0) then
         ok=0
      endif
   enddo
   ! table on the incident energies of the driving particle
   if (ok.eq.1.and.ipd.gt.0) then
      ! a LAW=6 driving particle has no table of incident
      ! energies; the energies of its yield are used (DLA, IAEA, Oct 2026)
      if (l6(ipd)) then
         ! geometric grid (ratio 1.02 or larger to fit nem points) from
         ! the first energy with nonzero yield to the last yield energy
         x1=yx(ipd,1)
         do jl=1,ny(ipd)
            if (yy(ipd,jl).gt.zero) then
               x1=yx(ipd,max(jl-1,1))
               exit
            endif
         enddo
         x2=yx(ipd,ny(ipd))
         if (x1.le.zero) x1=x2/1.e6_kr
         f=max(1.02_kr,(x2/x1)**(1/real(nem-1,kr)))
         nle=1
         elist6(1)=x1
         do while (elist6(nle)*f.lt.x2.and.nle.lt.nem-1)
            nle=nle+1
            elist6(nle)=elist6(nle-1)*f
         enddo
         nle=nle+1
         elist6(nle)=x2
      else
         nle=nein(ipd)
         do jl=1,min(nle,nem)
            elist6(jl)=ein(ipd,jl)
         enddo
      endif
      do ie=1,min(nle,nrmax)
         e=elist6(ie)
         if (e.le.zero) cycle
         p0=sqrt(2*e)
         ssum1=0
         sq=0
         s2=0
         ar=awr+1
         zr=zt
         do i=1,np
            y=yint(i,e)
            if (nein(i).eq.0.and..not.l6(i)) cycle
            call mint(i,e,pm,emean)
            ssum1=ssum1+y*pm
            sq=sq+y*pm*pm
            s2=s2+y*emean
            ar=ar-y*ap(i)
            zr=zr-y*zp(i)
         enddo
         call mint(ipd,e,pm,emean)
         a1=awr+1-ap(ipd)
         z1=zt-zp(ipd)
         er1=(p0*p0-2*p0*pm+emean)/(2*a1)
         er=(p0*p0-2*p0*ssum1+s2+ssum1*ssum1-sq)/(2*ar)
         nrt=nrt+1
         ert(nrt)=e
         rrt(nrt)=1
         if (ar.gt.0.5_kr.and.zr.gt.0.5_kr.and.er.gt.zero.and.er1.gt.zero) then
            dum=df(zero,zr,ar,zt,awr)
            dn=df(er,zr,ar,zt,awr)
            dum=df(zero,z1,a1,zt,awr)
            d1=df(er1,z1,a1,zt,awr)
            if (d1.gt.zero) rrt(nrt)=dn/d1
         endif
         ! exact recoil of the independent particles (DLA, IAEA, Oct 2026)
         call exrec(e,p0,a1,z1)
      enddo
      if (nfbk.gt.0) then
         write(strng,'(''mt'',i3,'': exact recoil at'',i5,'' of'',i5,&
           &'' energies'')') mtd,nexa,nexa+nfbk
         call mess('oprtab',strng,&
           'mean-energy factor used at the others (yields)')
      endif
   endif
   deallocate(b)
   deallocate(ndz,ndt,ndw,r7x,r7w,p7x,p7w)
   ! reposition the tape after the HEAD record of the section
   call findf(matd,6,mtd,nin)
   call contio(nin,0,0,hd,nb,nw)
   return
   contains
      subroutine stonod(ip,kk)
      ! compress the fine nodes of particle ip at its
      ! incident energy kk (cvcmp) and store them (ndz,ndt,ndw)
      ! (DLA, IAEA, Oct 2026)
      integer::ip,kk,mm,nn
      real(kr),dimension(:),allocatable::tmp
      call cvcmp(nfn,fnz,fnt,fnw,mm,cnz,cnt,cnw)
      if (nnod+mm.gt.size(ndz)) then
         nn=max(2*size(ndz),nnod+mm)
         allocate(tmp(nn))
         tmp(1:nnod)=ndz(1:nnod)
         call move_alloc(tmp,ndz)
         allocate(tmp(nn))
         tmp(1:nnod)=ndt(1:nnod)
         call move_alloc(tmp,ndt)
         allocate(tmp(nn))
         tmp(1:nnod)=ndw(1:nnod)
         call move_alloc(tmp,ndw)
      endif
      nodo(ip,kk)=nnod
      nodn(ip,kk)=mm
      if (mm.gt.0) then
         ndz(nnod+1:nnod+mm)=cnz(1:mm)
         ndt(nnod+1:nnod+mm)=cnt(1:mm)
         ndw(nnod+1:nnod+mm)=cnw(1:mm)
      endif
      nnod=nnod+mm
      end subroutine stonod
      subroutine addnod(xx,mu,ww,amp,vv)
      ! node of a particle of mass amp with energy xx and
      ! cosine mu in the frame of the data, moving with velocity vv
      ! along the beam (cm velocity, or 0 for laboratory data)
      ! (DLA, IAEA, Oct 2026)
      real(kr)::xx,mu,ww,amp,vv,uu
      uu=sqrt(2*max(xx,zero)/amp)
      call fnadd(amp*(vv+uu*mu),amp*uu*sqrt(max(1-mu*mu,zero)),ww)
      end subroutine addnod
      subroutine l1nod(ip)
      ! fine laboratory momentum nodes (fnadd) of the LAW=1
      ! table in b for particle ip: the discrete lines and the 4-point
      ! Gauss points of each continuum panel (as the moments above),
      ! each with its angular distribution (angn).  Frame as above.
      ! (DLA, IAEA, Oct 2026)
      integer::ip,jj,iqq
      real(kr)::amp,vv,xa,xb,xx,tt,wq
      nfn=0
      amp=ap(ip)
      vv=0
      if (.not.(lct.eq.1.or.(lct.eq.3.and.izap.gt.2004)))&
        vv=sqrt(2*en)/(awr+1)
      do jj=1,min(nd,nep)
         call angn(jj,0,zero,abs(b(6+ncyc*(jj-1)+1)),1.e0_kr,amp,vv)
      enddo
      do jj=nd+2,nep
         xa=b(6+ncyc*(jj-2)+1)
         xb=b(6+ncyc*(jj-1)+1)
         if (xb.le.xa) cycle
         do iqq=1,4
            xx=(xa+xb)/2+(xb-xa)/2*q4x(iqq)
            tt=(xx-xa)/(xb-xa)
            wq=q4w(iqq)*(xb-xa)/2
            if (lep.eq.1) then
               call angn(jj-1,0,zero,xx,wq,amp,vv)
            else
               call angn(jj-1,jj,tt,xx,wq,amp,vv)
            endif
         enddo
      enddo
      end subroutine l1nod
      subroutine angn(j1,j2,tt,xx,wq,amp,vv)
      ! nodes at the secondary energy xx with the
      ! distribution of point j1 (j2=0) or interpolated linearly
      ! between the points j1 and j2 (fraction tt of j2), weight wq:
      ! Legendre (LANG=1) and Kalbach (LANG=2, a from bacha as in
      ! h6ddx) at the 64 cosines of g64x; tabulated (LANG=11-15,
      ! histogram or linear in mu, log laws as linear) with a 4-point
      ! Gauss rule in each mu panel, normalized, times f0.  Negative
      ! values are set to zero (fnadd). (DLA, IAEA, Oct 2026)
      integer::j1,j2,iqq,il,lt,nmp,ia,kt,lb1,lb2,jt
      real(kr)::tt,xx,wq,amp,vv,f0v,dens,c,rr,aak,mu,u1,u2,v1,v2
      real(kr)::sn,pv,wt1
      lb1=6+ncyc*(j1-1)
      lb2=0
      if (j2.gt.0) lb2=6+ncyc*(j2-1)
      f0v=b(lb1+2)
      if (j2.gt.0) f0v=(1-tt)*b(lb1+2)+tt*b(lb2+2)
      if (lang.eq.1) then
         do iqq=1,64
            dens=f0v/2
            do il=1,min(na,64)
               c=b(lb1+2+il)
               if (j2.gt.0) c=(1-tt)*b(lb1+2+il)+tt*b(lb2+2+il)
               dens=dens+(2*il+1)*c*pleg(il+1,iqq)/2
            enddo
            call addnod(xx,g64x(iqq),wq*g64w(iqq)*dens,amp,vv)
         enddo
      else if (lang.eq.2) then
         rr=0
         if (na.ge.1) then
            rr=b(lb1+3)
            if (j2.gt.0) rr=(1-tt)*b(lb1+3)+tt*b(lb2+3)
         endif
         aak=bacha(1,izap,izat,en,xx)
         do iqq=1,64
            mu=g64x(iqq)
            if (aak.gt.1.e-6_kr) then
               dens=f0v*aak*(cosh(aak*mu)+rr*sinh(aak*mu))/(2*sinh(aak))
            else
               dens=f0v/2
            endif
            call addnod(xx,mu,wq*g64w(iqq)*dens,amp,vv)
         enddo
      else if (lang.ge.11.and.lang.le.15) then
         lt=lang-10
         nmp=na/2
         do kt=1,2
            if (kt.eq.1) then
               jt=lb1
               wt1=1
               if (j2.gt.0) wt1=1-tt
            else
               if (j2.eq.0) exit
               jt=lb2
               wt1=tt
            endif
            if (wt1.le.zero) cycle
            sn=0
            do ia=1,nmp-1
               u1=b(jt+2+2*ia-1)
               v1=b(jt+2+2*ia)
               u2=b(jt+2+2*ia+1)
               v2=b(jt+2+2*ia+2)
               if (u2.le.u1) cycle
               if (lt.eq.1) then
                  sn=sn+v1*(u2-u1)
               else
                  sn=sn+(v1+v2)*(u2-u1)/2
               endif
            enddo
            if (sn.le.zero) cycle
            do ia=1,nmp-1
               u1=b(jt+2+2*ia-1)
               v1=b(jt+2+2*ia)
               u2=b(jt+2+2*ia+1)
               v2=b(jt+2+2*ia+2)
               if (u2.le.u1) cycle
               do iqq=1,4
                  mu=(u1+u2)/2+(u2-u1)/2*q4x(iqq)
                  pv=v1
                  if (lt.ne.1) pv=v1+(v2-v1)*(mu-u1)/(u2-u1)
                  call addnod(xx,mu,wq*f0v*wt1*pv/sn*q4w(iqq)*(u2-u1)/2,&
                    amp,vv)
               enddo
            enddo
         enddo
      endif
      end subroutine angn
      subroutine r7row(n7)
      ! 4-point Gauss nodes (e', weight) of the LAW=7 TAB1
      ! in b (histogram panels exact, other laws on the linear form, as
      ! e7mom) in r7x, r7w (DLA, IAEA, Oct 2026)
      integer::n7,nr7,np7,ib7,jp,ir7,int7,iq7
      real(kr)::xa,xb,ya,yb,xm,ym,dx
      n7=0
      nr7=nint(b(5))
      np7=nint(b(6))
      if (size(r7x).lt.4*np7) then
         deallocate(r7x,r7w)
         allocate(r7x(4*np7),r7w(4*np7))
      endif
      ib7=6+2*nr7
      ir7=1
      do jp=2,np7
         xa=b(ib7+2*jp-3)
         ya=b(ib7+2*jp-2)
         xb=b(ib7+2*jp-1)
         yb=b(ib7+2*jp)
         do while (jp.gt.nint(b(5+2*ir7)).and.ir7.lt.nr7)
            ir7=ir7+1
         enddo
         int7=nint(b(6+2*ir7))
         dx=xb-xa
         if (dx.le.zero) cycle
         if (xa.lt.zero) xa=0
         do iq7=1,4
            xm=(xa+xb)/2+(xb-xa)/2*q4x(iq7)
            ym=ya
            if (int7.ne.1) ym=ya+(yb-ya)*(xm-xa)/(xb-xa)
            n7=n7+1
            r7x(n7)=xm
            r7w(n7)=q4w(iq7)*(xb-xa)/2*ym
         enddo
      enddo
      end subroutine r7row
      subroutine l6nod(ip,ee)
      ! LAW=6 nodes at the incident energy ee: phase space
      ! f(e') = sqrt(e') (emax-e')**(3n/2-4), isotropic in the cm,
      ! 64 x 64 Gauss points, emax as in mint (n<3: e'=emax)
      ! (DLA, IAEA, Oct 2026)
      integer::ip,i1,i2
      real(kr)::ee,vc,emx,pw,xx,fx
      nfn=0
      vc=sqrt(2*ee)/(awr+1)
      emx=0
      if (aps6(ip).gt.ap(ip).and.nps6(ip).gt.1)&
        emx=(aps6(ip)-ap(ip))/aps6(ip)*(awr*ee/(awr+1)+q)
      if (emx.le.zero) then
         call fnadd(ap(ip)*vc,zero,1.e0_kr)
         return
      endif
      if (nps6(ip).lt.3) then
         do i2=1,64
            call addnod(emx,g64x(i2),g64w(i2)/2,ap(ip),vc)
         enddo
         return
      endif
      pw=1.5e0_kr*nps6(ip)-4
      do i1=1,64
         xx=emx*(1+g64x(i1))/2
         fx=sqrt(xx)*(emx-xx)**pw*g64w(i1)*emx/2
         do i2=1,64
            call addnod(xx,g64x(i2),fx*g64w(i2)/2,ap(ip),vc)
         enddo
      enddo
      end subroutine l6nod
      subroutine getnod(ip,ee,m)
      ! compressed nodes (cnz,cnt,cnw; m nodes) of particle
      ! ip at the incident energy ee: LAW=6 computed at ee; otherwise
      ! the stored table at ee, or the two tables around ee mixed with
      ! the weight of the incident-energy law (law 1: lower table; 2,4:
      ! linear; 3,5: linear in ln e) and compressed again.  Outside the
      ! table the first or the last table is used (as mint).
      ! (DLA, IAEA, Oct 2026)
      integer::ip,m,jj,ka,kb,lw
      real(kr)::ee,tt
      m=0
      if (l6(ip)) then
         call l6nod(ip,ee)
         call cvcmp(nfn,fnz,fnt,fnw,m,cnz,cnt,cnw)
         return
      endif
      if (nein(ip).le.0) return
      ka=0
      kb=0
      tt=0
      do jj=1,nein(ip)
         if (abs(ee-ein(ip,jj)).le.1.e-9_kr*ee) then
            ka=jj
            exit
         endif
      enddo
      if (ka.eq.0) then
         if (ee.le.ein(ip,1)) then
            ka=1
         else if (ee.ge.ein(ip,nein(ip))) then
            ka=nein(ip)
         else
            do jj=2,nein(ip)
               if (ee.lt.ein(ip,jj)) then
                  ka=jj-1
                  kb=jj
                  lw=lwe(ip,jj-1)
                  tt=(ee-ein(ip,jj-1))/(ein(ip,jj)-ein(ip,jj-1))
                  if ((lw.eq.3.or.lw.eq.5).and.ein(ip,jj-1).gt.zero)&
                    tt=log(ee/ein(ip,jj-1))/log(ein(ip,jj)/ein(ip,jj-1))
                  if (lw.eq.1) tt=0
                  exit
               endif
            enddo
         endif
      endif
      if (kb.gt.0.and.tt.le.zero) kb=0
      if (kb.eq.0) then
         m=nodn(ip,ka)
         if (m.gt.0) then
            cnz(1:m)=ndz(nodo(ip,ka)+1:nodo(ip,ka)+m)
            cnt(1:m)=ndt(nodo(ip,ka)+1:nodo(ip,ka)+m)
            cnw(1:m)=ndw(nodo(ip,ka)+1:nodo(ip,ka)+m)
         endif
      else
         nfn=0
         do jj=nodo(ip,ka)+1,nodo(ip,ka)+nodn(ip,ka)
            call fnadd(ndz(jj),ndt(jj),(1-tt)*ndw(jj))
         enddo
         do jj=nodo(ip,kb)+1,nodo(ip,kb)+nodn(ip,kb)
            call fnadd(ndz(jj),ndt(jj),tt*ndw(jj))
         enddo
         call cvcmp(nfn,fnz,fnt,fnw,m,cnz,cnt,cnw)
      endif
      end subroutine getnod
      subroutine exrec(ee,pp0,aa1,zz1)
      ! exact mean damage of the residual after the
      ! independent emission of all the particles (integer
      ! multiplicities from the yields), dN, and of the one-particle
      ! recoil of the driving particle, d1, with the same nodes;
      ! rrt(nrt)=dN/d1.  With one particle in total the one-particle
      ! recoil is already exact (rrt=1 from the moments); with
      ! non-integer yields, more than nvmax particles or a particle
      ! without nodes the mean-energy factor above is kept (nfbk).
      ! (DLA, IAEA, Oct 2026)
      real(kr)::ee,pp0,aa1,zz1
      integer::ip,jj,jm
      real(kr)::yv
      nv=0
      lex=.true.
      do ip=1,np
         yv=yint(ip,ee)
         nmult(ip)=max(nint(yv),0)
         if (abs(yv-nmult(ip)).gt.1.e-4_kr*max(1.e0_kr,yv)) lex=.false.
         nv=nv+nmult(ip)
      enddo
      if (nv.le.1) return
      if (.not.lex.or.nv.gt.nvmax.or.nmult(ipd).lt.1) then
         nfbk=nfbk+1
         return
      endif
      nwk=0
      jm=0
      arx=awr+1
      zrx=zt
      io1=0
      n1=0
      do ip=1,np
         if (nmult(ip).eq.0) cycle
         call getnod(ip,ee,mn)
         if (mn.le.0) then
            nfbk=nfbk+1
            return
         endif
         wkz(nwk+1:nwk+mn)=cnz(1:mn)
         wkt(nwk+1:nwk+mn)=cnt(1:mn)
         wkw(nwk+1:nwk+mn)=cnw(1:mn)
         if (ip.eq.ipd) then
            io1=nwk
            n1=mn
         endif
         do jj=1,nmult(ip)
            jm=jm+1
            nsv(jm)=mn
            iov(jm)=nwk
         enddo
         nwk=nwk+mn
         arx=arx-nmult(ip)*ap(ip)
         zrx=zrx-nmult(ip)*zp(ip)
      enddo
      ! no residual nucleus (all the products are light particles):
      ! no correction, as with the moments
      if (arx.le.0.5_kr.or.zrx.le.0.5_kr) return
      if (n1.eq.0) then
         nfbk=nfbk+1
         return
      endif
      call cvdam(nv,nsv,iov,wkz,wkt,wkw,pp0,zrx,arx,zt,awr,dex,erm)
      nsv(1)=n1
      iov(1)=io1
      call cvdam(1,nsv,iov,wkz,wkt,wkw,pp0,zz1,aa1,zt,awr,d1x,erm)
      rrt(nrt)=1
      if (d1x.gt.zero) rrt(nrt)=dex/d1x
      nexa=nexa+1
      end subroutine exrec
      real(kr) function angpar(jj)
      ! angular parameter of point jj: Legendre f1 (not normalized),
      ! Kalbach r, or the mean cosine of a tabulated distribution
      integer::jj,ia,lt,nmp,lb
      real(kr)::u1,u2,v1,v2,sn,s1m
      lb=6+ncyc*(jj-1)
      angpar=0
      if (na.lt.1) return
      if (lang.eq.1.or.lang.eq.2) then
         angpar=b(lb+3)
      else if (lang.ge.11.and.lang.le.15) then
         lt=lang-10
         nmp=na/2
         sn=0
         s1m=0
         do ia=1,nmp-1
            u1=b(lb+2+2*ia-1)
            v1=b(lb+2+2*ia)
            u2=b(lb+2+2*ia+1)
            v2=b(lb+2+2*ia+2)
            if (u2.le.u1) cycle
            if (lt.eq.1) then
               sn=sn+v1*(u2-u1)
               s1m=s1m+v1*(u2*u2-u1*u1)/2
            else
               ! linear (log laws approximated as linear)
               sn=sn+(v1+v2)*(u2-u1)/2
               s1m=s1m+(u2-u1)*(v1*(2*u1+u2)+v2*(u1+2*u2))/6
            endif
         enddo
         if (sn.gt.zero) angpar=s1m/sn
      endif
      end function angpar
      real(kr) function gmu(ff,par,xe)
      ! f*<mu> from the distribution value ff and the angular parameter
      real(kr)::ff,par,xe,a,t
      gmu=0
      if (na.lt.1) return
      if (lang.eq.1) then
         gmu=par
      else if (lang.eq.2) then
         a=bacha(1,izap,izat,en,xe)
         if (a.gt.1.e-6_kr) then
            gmu=ff*par*(1/tanh(a)-1/a)
         else
            gmu=ff*par*a/3
         endif
      else if (lang.ge.11.and.lang.le.15) then
         gmu=ff*par
      endif
      t=0
      if (ff.gt.zero) t=gmu/ff
      if (t.gt.1) gmu=ff
      if (t.lt.-1) gmu=-ff
      end function gmu
      real(kr) function yint(ip,ee)
      integer::ip,jj,lr,law
      real(kr)::ee
      yint=yy(ip,1)
      if (ny(ip).le.1.or.ee.le.yx(ip,1)) return
      yint=yy(ip,ny(ip))
      if (ee.ge.yx(ip,ny(ip))) return
      do jj=2,ny(ip)
         if (ee.le.yx(ip,jj)) then
            law=2
            do lr=nry(ip),1,-1
               if (jj.le.nby(ip,lr)) law=iny(ip,lr)
            enddo
            if (yx(ip,jj).gt.yx(ip,jj-1)) then
               call terp1(yx(ip,jj-1),yy(ip,jj-1),yx(ip,jj),yy(ip,jj),ee,yint,law)
            else
               yint=yy(ip,jj)
            endif
            return
         endif
      enddo
      end function yint
      subroutine mint(ip,ee,av,cv)
      integer::ip,jj,lw
      real(kr)::ee,av,cv,emx,vc
      if (l6(ip)) then
         ! LAW=6, laboratory moments of the phase space (DLA, IAEA, Oct 2026)
         vc=sqrt(2*ee)/(awr+1)
         emx=0
         if (aps6(ip).gt.ap(ip).and.nps6(ip).gt.1)&
           emx=(aps6(ip)-ap(ip))/aps6(ip)*(awr*ee/(awr+1)+q)
         if (emx.lt.zero) emx=0
         av=ap(ip)*vc
         cv=2*ap(ip)*(emx/(nps6(ip)-1)+ap(ip)*ee/(awr+1)**2)
         return
      endif
      av=pmu(ip,1)
      cv=p2(ip,1)
      if (nein(ip).le.1.or.ee.le.ein(ip,1)) return
      av=pmu(ip,nein(ip))
      cv=p2(ip,nein(ip))
      if (ee.ge.ein(ip,nein(ip))) return
      do jj=2,nein(ip)
         if (ee.le.ein(ip,jj)) then
            lw=lwe(ip,jj-1)
            if ((lw.eq.4.or.lw.eq.5).and.pmu(ip,jj-1)*pmu(ip,jj).le.zero) lw=lw-2
            if (lw.lt.1.or.lw.gt.5) lw=2
            call terp1(ein(ip,jj-1),pmu(ip,jj-1),ein(ip,jj),pmu(ip,jj),ee,av,lw)
            lw=lwe(ip,jj-1)
            if (lw.lt.1.or.lw.gt.5) lw=2
            call terp1(ein(ip,jj-1),p2(ip,jj-1),ein(ip,jj),p2(ip,jj),ee,cv,lw)
            return
         endif
      enddo
      end subroutine mint
      subroutine e7mom(am7,g0,g1,gp)
      ! integrals over e' of one LAW=7 TAB1 f(e') stored in
      ! b: g0=int f, g1=int f e', gp=int f sqrt(2 am7 e'); histogram
      ! panels exactly, other laws with 4-point Gauss on the linear form
      ! (DLA, IAEA, Oct 2026)
      real(kr)::am7,g0,g1,gp
      integer::nr7,np7,ib7,jp,ir7,int7,iq7
      real(kr)::xa,xb,ya,yb,xm,ym,dx
      g0=0
      g1=0
      gp=0
      nr7=nint(b(5))
      np7=nint(b(6))
      ib7=6+2*nr7
      ir7=1
      do jp=2,np7
         xa=b(ib7+2*jp-3)
         ya=b(ib7+2*jp-2)
         xb=b(ib7+2*jp-1)
         yb=b(ib7+2*jp)
         do while (jp.gt.nint(b(5+2*ir7)).and.ir7.lt.nr7)
            ir7=ir7+1
         enddo
         int7=nint(b(6+2*ir7))
         dx=xb-xa
         if (dx.le.zero) cycle
         if (xa.lt.zero) xa=0
         if (int7.eq.1) then
            g0=g0+ya*dx
            g1=g1+ya*(xb*xb-xa*xa)/2
            gp=gp+ya*sqrt(2*am7)*2*(xb**1.5_kr-xa**1.5_kr)/3
         else
            do iq7=1,4
               xm=(xa+xb)/2+dx/2*q4x(iq7)
               ym=ya+(yb-ya)*(xm-xa)/dx
               g0=g0+q4w(iq7)*dx/2*ym
               g1=g1+q4w(iq7)*dx/2*ym*xm
               gp=gp+q4w(iq7)*dx/2*ym*sqrt(2*am7*max(xm,zero))
            enddo
         endif
      enddo
      end subroutine e7mom
   end subroutine oprtab

   real(kr) function oprfac(e)
   !-------------------------------------------------------------------
   ! Interpolate the OPR damage correction factor of oprtab.
   !-------------------------------------------------------------------
   real(kr)::e
   integer::j
   oprfac=1
   if (nrt.le.0) return
   if (e.le.ert(1)) then
      oprfac=rrt(1)
   else if (e.ge.ert(nrt)) then
      oprfac=rrt(nrt)
   else
      do j=2,nrt
         if (e.le.ert(j)) then
            oprfac=rrt(j-1)+(rrt(j)-rrt(j-1))*(e-ert(j-1))/(ert(j)-ert(j-1))
            return
         endif
      enddo
   endif
   return
   end function oprfac

   subroutine mt5tab(nin,matd,mtd,za,awr,lct,nk,thr)
   !-------------------------------------------------------------------
   ! Recoil of the residual nuclei of mt5 given in File 6
   ! without distribution (LAW=0), or not given at all (only light
   ! particles and photons).  This corrects incomplete evaluations;
   ! residuals given with their distributions are used as given.
   ! (DLA, IAEA, Oct 2026)
   !
   ! Momentum balance for the light particles i (zap=1..2004, LAW=1
   ! or LAW=7; masses m_i=AWP, yields y_i), laboratory moments:
   !   p0 = sqrt(2 e)                       (neutron mass units, eV)
   !   S2 = sum y_i <p_i**2>,  Sz = sum y_i <p_i mu_i>,  ym = sum y_i m_i
   !   Smz = sum y_i m_i <p_i mu_i>,  Sm2 = sum y_i m_i**2
   ! Drift along the beam  pd = p0 - Sz.  Second moment: the
   ! momenta of the light particles in the cm, q_i = p_i - m_i V with
   ! V = p0/(awr+1), are taken without cross terms (cm-exclusive
   ! emission; the directed part of the lab momenta keeps its
   ! correlation):
   !   D = Sz - ym V,   Q2 = S2 - 2 V Smz + V**2 Sm2 = sum y_i <q_i**2>
   !   <pr**2> = ((awr+1-ym) V)**2 - 2 (awr+1-ym) V D + Q2,
   !   bounded below by pd**2.
   ! Validation against residual recoil spectra evaluated in JENDL-5
   ! (A=19-209, 25-200 MeV): IAEA report on the NJOY2016 HEATR
   ! corrections (2026).
   !   sg**2 = (<pr**2> - pd**2)/3 (isotropic dispersion)
   ! The residual momentum is a displaced 3-d gaussian, pr = pd z + d,
   ! d ~ N(0, sg**2 I); the density of |pr| is
   !   f(p) = p/(pd sg sqrt(2 pi)) [exp(-(p-pd)**2/2sg**2)
   !                              - exp(-(p+pd)**2/2sg**2)]
   ! (maxwellian when pd -> 0), and er = p**2/(2 M).  For each
   ! residual r without distribution (yield y_r, Z_r, M_r=AWP):
   !   damage per mt5 event  d5 = sum y_r <df(p**2/2M_r; Z_r, M_r)>
   !   heating per mt5 event h5 = sum y_r <pr**2>/(2 M_r)
   ! (64-point Gauss rule in p on pd-8sg..pd+8sg).  At the energies
   ! where the section gives no residual (sum of residual yields below
   ! 0.01, or no residual subsection), one effective residual with
   ! yield 1,  Z_eff = Z - sum y_i z_i,  M_eff = awr + 1 - sum y_i awp_i,
   ! is used.  Moments of the light particles: LAW=1 (all LANG, as in
   ! oprtab, frame from LCT) and LAW=7 (laboratory; 4-point Gauss rule
   ! in each E' panel, linear in mu).  The table (e5t,d5t,h5t) is built
   ! on the incident energies of the light particles; n5t=0 when no
   ! correction is needed or possible.  The charge balance
   ! sum y Z(all products) = Z is checked and a message is given for a
   ! deviation larger than 5% (inconsistent yields) at energies above
   ! the mt5 threshold thr.
   ! The tape is repositioned after the HEAD record of the section.
   !-------------------------------------------------------------------
   use endf
   use util   ! provides error,mess
   use mainio ! provides nsyso
   ! externals
   integer::nin,matd,mtd,lct,nk
   real(kr)::za,awr,thr
   ! internals
   integer,parameter::npm=12,nem=1000,nym=2000,nrm=20,nrsm=1500
   integer,parameter::npool=600000,nbmax=400000
   integer::np,ik,nb,nw,l,law,lang,lep,ne,ie,nd,na,nep,ncyc,i,j,k,iq
   integer::izap,izat,nr,ny(npm),nein(npm),nmu,imu,lw,nr2,lfr
   integer::nbt2(nrm),int2(nrm),nry(npm),nby(npm,nrm),iny(npm,nrm)
   integer::lwe(npm,nem),nres,nres0,nu,iu,ipool,inc,lmu,ir,jr,nrr
   integer::lawr(nrsm),ioff(nrsm),nyr(nrsm),nrrr(nrsm)
   integer::nbr(nrsm,4),inr(nrsm,4)
   real(kr)::zp(npm),ap(npm),yx(npm,nym),yy(npm,nym)
   real(kr)::ein(npm,nem),pmu(npm,nem),p2(npm,nem)
   real(kr)::zr5(nrsm),ar5(nrsm)
   real(kr)::en,ep,f,x,s0,s1,sm,am,vcm,emean,x1,x2,fa,fb,ga,gb,g,fx
   real(kr)::e,p0,y,sz,s2,ym,yz,pr2,pd,sg2,dsum,hsum,yres,qz,zt,dev
   real(kr)::devmax,edev,zeff,aeff,yr,pm,cv,emu,emup,g0,g1,gs
   real(kr)::smz,sm2,vc5,dz5,sq2
   real(kr)::g0p,g1p,gsp,hd(6),ylight,yall
   real(kr),allocatable::b(:),pool(:),eu(:),yt(:)
   integer::nadd0,neff
   real(kr),parameter::ymin5=0.01e0_kr
   character(60)::strng,strng2
   real(kr),parameter::zero=0
   real(kr),parameter::tiny=1.e-8_kr
   real(kr),dimension(4),parameter::q4x=(/&
      -.861136312e0_kr,-.339981044e0_kr,.339981044e0_kr,.861136312e0_kr/)
   real(kr),dimension(4),parameter::q4w=(/&
      .347854845e0_kr,.652145155e0_kr,.652145155e0_kr,.347854845e0_kr/)

   n5t=0
   allocate(b(nbmax))
   allocate(pool(npool))
   izat=nint(za)
   zt=int(za/1000)
   np=0
   nres=0
   nres0=0
   ipool=0
   inc=0
   do ik=1,nk
      l=1
      call tab1io(nin,0,0,b(l),nb,nw)
      l=l+nw
      do while (nb.ne.0)
         if (l.gt.nbmax-10000) call error('mt5tab','storage',' ')
         call moreio(nin,0,0,b(l),nb,nw)
         l=l+nw
      enddo
      izap=nint(b(1))
      law=nint(b(4))
      nr=nint(b(5))
      i=0
      if (izap.gt.0.and.izap.le.2004) then
         if (np.ge.npm) call error('mt5tab','too many light particles',' ')
         np=np+1
         i=np
         zp(np)=int(izap/1000)
         ap(np)=b(2)
         ny(np)=min(nint(b(6)),nym)
         nry(np)=min(nr,nrm)
         do j=1,nry(np)
            nby(np,j)=nint(b(5+2*j))
            iny(np,j)=nint(b(6+2*j))
         enddo
         do j=1,ny(np)
            yx(np,j)=b(6+2*nr+2*j-1)
            yy(np,j)=b(6+2*nr+2*j)
         enddo
         nein(np)=0
         if (law.ne.1.and.law.ne.7) then
            inc=1
            ! message for any light particle not LAW=1 or 7
            ! (DLA, IAEA, Oct 2026)
            write(strng,'(''mt5 light particle'',i6,'' given with LAW='',&
              &i3)') izap,law
            call mess('mt5tab',strng,&
              'moments not computed (LAW=1 or 7 only)')
         endif
      else if (izap.gt.2004) then
         if (nres.ge.nrsm) call error('mt5tab','too many residuals',' ')
         nres=nres+1
         lawr(nres)=law
         if (law.eq.0) nres0=nres0+1
         ! unexpected representation of a residual (DLA, IAEA, Oct 2026)
         if (law.eq.6.or.law.lt.0) then
            write(strng,'(''mt5 residual'',i7,'' given with LAW='',i3)')&
              izap,law
            call mess('mt5tab',strng,'unexpected law for a residual')
         endif
         zr5(nres)=int(izap/1000)
         ar5(nres)=b(2)
         nyr(nres)=nint(b(6))
         nrrr(nres)=min(nr,4)
         do j=1,nrrr(nres)
            nbr(nres,j)=nint(b(5+2*j))
            inr(nres,j)=nint(b(6+2*j))
         enddo
         if (ipool+2*nyr(nres).gt.npool) call error('mt5tab',&
           'yield storage exceeded',' ')
         ioff(nres)=ipool
         do j=1,2*nyr(nres)
            pool(ipool+j)=b(6+2*nr+j)
         enddo
         ipool=ipool+2*nyr(nres)
      endif
      if (law.eq.1.or.law.eq.2.or.law.eq.5) then
         call tab2io(nin,0,0,b(1),nb,nw)
         lang=nint(b(3))
         lep=nint(b(4))
         ne=nint(b(6))
         nr2=min(nint(b(5)),nrm)
         do j=1,nr2
            nbt2(j)=nint(b(5+2*j))
            int2(j)=nint(b(6+2*j))
         enddo
         do ie=1,ne
            l=1
            call listio(nin,0,0,b(l),nb,nw)
            l=l+nw
            do while (nb.ne.0)
               if (l.gt.nbmax-10000) call error('mt5tab','storage',' ')
               call moreio(nin,0,0,b(l),nb,nw)
               l=l+nw
            enddo
            if (i.gt.0.and.law.eq.1.and.nein(max(i,1)).lt.nem) then
               en=b(2)
               nd=nint(b(3))
               na=nint(b(4))
               nep=nint(b(6))
               ncyc=na+2
               am=ap(i)
               s0=0
               s1=0
               sm=0
               do j=1,min(nd,nep)
                  ep=b(6+ncyc*(j-1)+1)
                  f=b(6+ncyc*(j-1)+2)
                  s0=s0+f
                  s1=s1+f*ep
                  sm=sm+gmu(f,angpar(j),ep)*sqrt(2*am*ep)
               enddo
               do j=nd+2,nep
                  x1=b(6+ncyc*(j-2)+1)
                  x2=b(6+ncyc*(j-1)+1)
                  if (x2.le.x1) cycle
                  fa=b(6+ncyc*(j-2)+2)
                  fb=b(6+ncyc*(j-1)+2)
                  do iq=1,4
                     x=(x1+x2)/2+(x2-x1)/2*q4x(iq)
                     ga=angpar(j-1)
                     if (lep.eq.1) then
                        fx=fa
                        g=ga
                     else
                        gb=angpar(j)
                        fx=fa+(fb-fa)*(x-x1)/(x2-x1)
                        g=ga+(gb-ga)*(x-x1)/(x2-x1)
                     endif
                     g=gmu(fx,g,x)
                     s0=s0+q4w(iq)*(x2-x1)/2*fx
                     s1=s1+q4w(iq)*(x2-x1)/2*fx*x
                     sm=sm+q4w(iq)*(x2-x1)/2*g*sqrt(2*am*x)
                  enddo
               enddo
               if (s0.gt.zero) then
                  lw=int2(nr2)
                  do j=nr2,1,-1
                     if (ie+1.le.nbt2(j)) lw=int2(j)
                  enddo
                  if (lw.gt.20) lw=lw-20
                  ! frame of the data: laboratory for LCT=1 and for
                  ! particles heavier than alpha with LCT=3, else CM
                  lfr=2
                  if (lct.eq.1) lfr=1
                  call store(i,en,lw,s1/s0,sm/s0,lfr)
               endif
            endif
         enddo
      else if (law.eq.6) then
         call contio(nin,0,0,b(1),nb,nw)
      else if (law.eq.7) then
         call tab2io(nin,0,0,b(1),nb,nw)
         ne=nint(b(6))
         nr2=min(nint(b(5)),nrm)
         do j=1,nr2
            nbt2(j)=nint(b(5+2*j))
            int2(j)=nint(b(6+2*j))
         enddo
         do ie=1,ne
            call tab2io(nin,0,0,b(1),nb,nw)
            en=b(2)
            nmu=nint(b(6))
            lmu=nint(b(8))
            s0=0
            s1=0
            sm=0
            emup=0
            g0p=0
            g1p=0
            gsp=0
            do imu=1,nmu
               l=1
               call tab1io(nin,0,0,b(l),nb,nw)
               l=l+nw
               do while (nb.ne.0)
                  if (l.gt.nbmax-10000) call error('mt5tab','storage',' ')
                  call moreio(nin,0,0,b(l),nb,nw)
                  l=l+nw
               enddo
               if (i.eq.0) cycle
               emu=b(2)
               call epint(g0,g1,gs)
               if (imu.gt.1.and.emu.gt.emup) then
                  if (lmu.eq.1) then
                     s0=s0+g0p*(emu-emup)
                     s1=s1+g1p*(emu-emup)
                     sm=sm+gsp*(emu*emu-emup*emup)/2
                  else
                     s0=s0+(g0p+g0)*(emu-emup)/2
                     s1=s1+(g1p+g1)*(emu-emup)/2
                     sm=sm+(emu-emup)*(gsp*(2*emup+emu)+gs*(emup+2*emu))/6
                  endif
               endif
               emup=emu
               g0p=g0
               g1p=g1
               gsp=gs
            enddo
            if (i.gt.0.and.s0.gt.zero.and.nein(max(i,1)).lt.nem) then
               lw=int2(nr2)
               do j=nr2,1,-1
                  if (ie+1.le.nbt2(j)) lw=int2(j)
               enddo
               if (lw.gt.20) lw=lw-20
               call store(i,en,lw,s1/s0,sqrt(2*ap(i))*sm/s0,1)
            endif
         enddo
      else if (law.lt.0) then
         inc=1
      endif
   enddo

   ! union of the incident energies of the light particles
   nu=0
   allocate(eu(npm*nem))
   do i=1,np
      do k=1,nein(i)
         nu=nu+1
         eu(nu)=ein(i,k)
      enddo
   enddo
   call sortu(eu,nu)
   if (nu.gt.nrmax) nu=nrmax

   ! charge balance of all products, at energies above the threshold
   ! where residual yields are given (sum of residual yields >= 0.01)
   devmax=0
   edev=0
   do iu=1,nu
      e=eu(iu)
      if (e.lt.thr) cycle
      qz=0
      ylight=0
      do i=1,np
         y=yint(i,e)
         qz=qz+y*zp(i)
         ylight=ylight+y
      enddo
      yall=0
      do k=1,nres
         yr=yres5(k,e)
         qz=qz+yr*zr5(k)
         yall=yall+yr
      enddo
      if (ylight.gt.zero.and.yall.ge.ymin5) then
         dev=(qz-zt)/zt
         if (abs(dev).gt.abs(devmax)) then
            devmax=dev
            edev=e
         endif
      endif
   enddo
   if (abs(devmax).gt.0.05_kr) then
      write(strng,'(''mt5 charge balance off by'',1p,e10.2,''%'')')&
        100*devmax
      write(strng2,'(''at e='',1p,e10.3,'' ev: inconsistent yields'')')&
        edev
      call mess('mt5tab',strng,strng2)
   endif

   ! recoil of the residuals without distribution (LAW=0), or of one
   ! effective residual at the energies where no residual is given
   if (np.gt.0.and.inc.ne.0.and.(nres0.gt.0.or.nres.eq.0)) then
      call mess('mt5tab','mt5 light particle without LAW=1 or 7',&
        'no recoil for residuals without distribution')
   else if (np.gt.0) then
      allocate(yt(nu))
      nadd0=0
      neff=0
      do iu=1,nu
         e=eu(iu)
         p0=sqrt(2*e)
         sz=0
         s2=0
         ym=0
         yz=0
         smz=0
         sm2=0
         do i=1,np
            y=yint(i,e)
            if (y.le.zero.or.nein(i).eq.0) cycle
            call mint(i,e,pm,cv)
            sz=sz+y*pm
            s2=s2+y*cv
            ym=ym+y*ap(i)
            yz=yz+y*zp(i)
            smz=smz+y*ap(i)*pm
            sm2=sm2+y*ap(i)*ap(i)
         enddo
         ! no cross terms between the cm momenta q_i=p_i-m_i V
         ! (DLA, IAEA, Oct 2026)
         vc5=p0/(awr+1)
         dz5=sz-ym*vc5
         sq2=s2-2*vc5*smz+vc5*vc5*sm2
         pr2=((awr+1-ym)*vc5)**2-2*(awr+1-ym)*vc5*dz5+sq2
         pd=p0-sz
         if (pr2.lt.pd*pd) pr2=pd*pd
         sg2=(pr2-pd*pd)/3
         dsum=0
         hsum=0
         yres=0
         yall=0
         do k=1,nres
            yall=yall+yres5(k,e)
         enddo
         if (pr2.gt.zero.and.ym.gt.zero.and.e.ge.thr) then
            if (yall.ge.ymin5) then
               do k=1,nres
                  if (lawr(k).ne.0) cycle
                  yr=yres5(k,e)
                  if (yr.le.zero) cycle
                  if (idame.gt.0) dsum=dsum+yr*dnc(ar5(k),zr5(k))
                  hsum=hsum+yr*pr2/(2*ar5(k))
                  yres=yres+yr
               enddo
               if (yres.gt.zero) nadd0=nadd0+1
            else
               zeff=zt-yz
               aeff=awr+1-ym
               if (zeff.ge.1.and.aeff.ge.1) then
                  if (idame.gt.0) dsum=dnc(aeff,zeff)
                  hsum=pr2/(2*aeff)
                  yres=1
                  neff=neff+1
               endif
            endif
         endif
         n5t=n5t+1
         e5t(n5t)=e
         d5t(n5t)=dsum
         h5t(n5t)=hsum
         yt(n5t)=yres
      enddo
      if (nadd0.eq.0.and.neff.eq.0) then
         n5t=0
      else
         if (nadd0.gt.0) then
            write(strng,'(''mt5 has'',i4,'' residuals without '',&
              &''distribution'')') nres0
            call mess('mt5tab',strng,&
              'recoil from momentum balance')
         endif
         if (neff.gt.0) then
            write(strng,'(''mt5 gives no residual nuclei at'',i5,&
              &'' energies'')') neff
            call mess('mt5tab',strng,&
              'effective residual from charge/mass balance')
         endif
         if (iprint.eq.1) then
            write(nsyso,'(/&
              &'' mt5: recoil of residuals without distribution''/&
              &''           e       yield     <er>(ev)   damage(ev)'')')
            do iu=1,n5t
               if (yt(iu).gt.zero) write(nsyso,'(1x,1p,4e12.4)')&
                 e5t(iu),yt(iu),h5t(iu)/yt(iu),d5t(iu)
            enddo
         endif
      endif
      deallocate(yt)
   endif
   deallocate(eu)
   deallocate(pool)
   deallocate(b)
   ! reposition the tape after the HEAD record of the section
   call findf(matd,6,mtd,nin)
   call contio(nin,0,0,hd,nb,nw)
   return

   contains

      subroutine store(ip,ee,lww,em,pmm,lf)
      ! moments of particle ip at incident energy ee, transformed to
      ! the laboratory when given in the centre of mass (lf=2)
      integer::ip,lww,lf,kk
      real(kr)::ee,em,pmm,v
      nein(ip)=nein(ip)+1
      kk=nein(ip)
      ein(ip,kk)=ee
      lwe(ip,kk)=lww
      pmu(ip,kk)=pmm
      p2(ip,kk)=2*ap(ip)*em
      if (lf.eq.2) then
         v=sqrt(2*ee)/(awr+1)
         p2(ip,kk)=p2(ip,kk)+2*ap(ip)*v*pmu(ip,kk)+ap(ip)*ap(ip)*v*v
         pmu(ip,kk)=pmu(ip,kk)+ap(ip)*v
      endif
      end subroutine store

      subroutine epint(h0,h1,hs)
      ! integrals of f, E' f and sqrt(E') f over a TAB1 f(E') in b
      ! (4-point Gauss rule per panel, interpolation laws of the
      ! record; log laws approximated by linear interpolation)
      real(kr)::h0,h1,hs,xa,xb,ya,yb,xx,ff
      integer::nrx,npx,jj,kr2,lwx
      h0=0
      h1=0
      hs=0
      nrx=nint(b(5))
      npx=nint(b(6))
      do jj=2,npx
         xa=b(6+2*nrx+2*jj-3)
         ya=b(6+2*nrx+2*jj-2)
         xb=b(6+2*nrx+2*jj-1)
         yb=b(6+2*nrx+2*jj)
         if (xb.le.xa) cycle
         lwx=2
         do kr2=nrx,1,-1
            if (jj.le.nint(b(5+2*kr2))) lwx=nint(b(6+2*kr2))
         enddo
         do iq=1,4
            xx=(xa+xb)/2+(xb-xa)/2*q4x(iq)
            if (lwx.eq.1) then
               ff=ya
            else
               ff=ya+(yb-ya)*(xx-xa)/(xb-xa)
            endif
            h0=h0+q4w(iq)*(xb-xa)/2*ff
            h1=h1+q4w(iq)*(xb-xa)/2*ff*xx
            hs=hs+q4w(iq)*(xb-xa)/2*ff*sqrt(xx)
         enddo
      enddo
      end subroutine epint

      real(kr) function angpar(jj)
      ! angular parameter of point jj (as in oprtab)
      integer::jj,ia,lt,nmp,lb
      real(kr)::u1,u2,v1,v2,sn,s1m
      lb=6+ncyc*(jj-1)
      angpar=0
      if (na.lt.1) return
      if (lang.eq.1.or.lang.eq.2) then
         angpar=b(lb+3)
      else if (lang.ge.11.and.lang.le.15) then
         lt=lang-10
         nmp=na/2
         sn=0
         s1m=0
         do ia=1,nmp-1
            u1=b(lb+2+2*ia-1)
            v1=b(lb+2+2*ia)
            u2=b(lb+2+2*ia+1)
            v2=b(lb+2+2*ia+2)
            if (u2.le.u1) cycle
            if (lt.eq.1) then
               sn=sn+v1*(u2-u1)
               s1m=s1m+v1*(u2*u2-u1*u1)/2
            else
               sn=sn+(v1+v2)*(u2-u1)/2
               s1m=s1m+(u2-u1)*(v1*(2*u1+u2)+v2*(u1+2*u2))/6
            endif
         enddo
         if (sn.gt.zero) angpar=s1m/sn
      endif
      end function angpar

      real(kr) function gmu(ff,par,xe)
      ! f*<mu> from the distribution value ff and the angular parameter
      real(kr)::ff,par,xe,a,t
      gmu=0
      if (na.lt.1) return
      if (lang.eq.1) then
         gmu=par
      else if (lang.eq.2) then
         a=bacha(1,izap,izat,en,xe)
         if (a.gt.1.e-6_kr) then
            gmu=ff*par*(1/tanh(a)-1/a)
         else
            gmu=ff*par*a/3
         endif
      else if (lang.ge.11.and.lang.le.15) then
         gmu=ff*par
      endif
      t=0
      if (ff.gt.zero) t=gmu/ff
      if (t.gt.1) gmu=ff
      if (t.lt.-1) gmu=-ff
      end function gmu

      real(kr) function yint(ip,ee)
      ! yield of light particle ip (laws of the record)
      integer::ip,jj,lr,lwy
      real(kr)::ee
      yint=0
      if (ny(ip).lt.1) return
      if (ee.lt.yx(ip,1).or.ee.gt.yx(ip,ny(ip))) return
      yint=yy(ip,ny(ip))
      if (ee.ge.yx(ip,ny(ip))) return
      ! strict test: at a discontinuity the value above is used
      do jj=2,ny(ip)
         if (ee.lt.yx(ip,jj)) then
            lwy=2
            do lr=nry(ip),1,-1
               if (jj.le.nby(ip,lr)) lwy=iny(ip,lr)
            enddo
            if (yx(ip,jj).gt.yx(ip,jj-1)) then
               call terp1(yx(ip,jj-1),yy(ip,jj-1),yx(ip,jj),yy(ip,jj),&
                 ee,yint,lwy)
            else
               yint=yy(ip,jj)
            endif
            return
         endif
      enddo
      end function yint

      real(kr) function yres5(kr5,ee)
      ! yield of residual kr5 (laws of the record)
      integer::kr5,jj,lr,lwy,io
      real(kr)::ee
      yres5=0
      io=ioff(kr5)
      if (nyr(kr5).lt.1) return
      if (ee.lt.pool(io+1).or.ee.gt.pool(io+2*nyr(kr5)-1)) return
      ! strict test: at a discontinuity the value above is used
      do jj=2,nyr(kr5)
         if (ee.lt.pool(io+2*jj-1)) then
            lwy=2
            do lr=nrrr(kr5),1,-1
               if (jj.le.nbr(kr5,lr)) lwy=inr(kr5,lr)
            enddo
            if (pool(io+2*jj-1).gt.pool(io+2*jj-3)) then
               call terp1(pool(io+2*jj-3),pool(io+2*jj-2),&
                 pool(io+2*jj-1),pool(io+2*jj),ee,yres5,lwy)
            else
               yres5=pool(io+2*jj)
            endif
            return
         endif
      enddo
      yres5=pool(io+2*nyr(kr5))
      end function yres5

      subroutine mint(ip,ee,av,cv2)
      ! moments <p mu> (av) and <p**2> (cv2) of particle ip at ee
      integer::ip,jj,lww
      real(kr)::ee,av,cv2
      av=pmu(ip,1)
      cv2=p2(ip,1)
      if (nein(ip).le.1.or.ee.le.ein(ip,1)) return
      av=pmu(ip,nein(ip))
      cv2=p2(ip,nein(ip))
      if (ee.ge.ein(ip,nein(ip))) return
      do jj=2,nein(ip)
         if (ee.le.ein(ip,jj)) then
            lww=lwe(ip,jj-1)
            if ((lww.eq.4.or.lww.eq.5).and.pmu(ip,jj-1)*pmu(ip,jj).le.zero)&
              lww=lww-2
            if (lww.lt.1.or.lww.gt.5) lww=2
            call terp1(ein(ip,jj-1),pmu(ip,jj-1),ein(ip,jj),pmu(ip,jj),&
              ee,av,lww)
            lww=lwe(ip,jj-1)
            if (lww.lt.1.or.lww.gt.5) lww=2
            call terp1(ein(ip,jj-1),p2(ip,jj-1),ein(ip,jj),p2(ip,jj),&
              ee,cv2,lww)
            return
         endif
      enddo
      end subroutine mint

      real(kr) function dnc(am5,zr)
      ! mean damage energy of a residual (mass am5, charge zr) with the
      ! displaced gaussian momentum distribution (pd, sg2, pr2)
      real(kr)::am5,zr,s,plo,phi,pp,w,ff,sn,dum
      real(kr),parameter::rtp=2.506628274631e0_kr
      dnc=0
      dum=df(zero,zr,am5,zt,awr)
      if (sg2.le.tiny*pr2) then
         dnc=df(pd*pd/(2*am5),zr,am5,zt,awr)
         return
      endif
      s=sqrt(sg2)
      plo=pd-8*s
      if (plo.lt.zero) plo=0
      phi=pd+8*s
      sn=0
      do iq=1,64
         pp=(plo+phi)/2+(phi-plo)/2*g64x(iq)
         w=g64w(iq)*(phi-plo)/2
         if (pd.gt.tiny*s) then
            ff=pp/(pd*s*rtp)*(exp(-(pp-pd)**2/(2*sg2))&
              -exp(-(pp+pd)**2/(2*sg2)))
         else
            ff=2*pp*pp/(s**3*rtp)*exp(-pp*pp/(2*sg2))
         endif
         sn=sn+w*ff
         dnc=dnc+w*ff*df(pp*pp/(2*am5),zr,am5,zt,awr)
      enddo
      if (sn.gt.zero) dnc=dnc/sn
      end function dnc

   end subroutine mt5tab

   subroutine sortu(a,n)
   !-------------------------------------------------------------------
   ! Sort a(1:n) in increasing order and remove duplicates.
   !-------------------------------------------------------------------
   integer::n
   real(kr)::a(*)
   integer::i,j,m
   real(kr)::t
   do i=2,n
      t=a(i)
      j=i-1
      do while (j.ge.1)
         if (a(j).le.t) exit
         a(j+1)=a(j)
         j=j-1
      enddo
      a(j+1)=t
   enddo
   m=0
   do i=1,n
      if (m.eq.0) then
         m=1
         a(1)=a(i)
      else if (a(i).gt.a(m)*(1+1.e-9_kr)) then
         m=m+1
         a(m)=a(i)
      endif
   enddo
   n=m
   return
   end subroutine sortu

   real(kr) function mt5dm(e)
   !-------------------------------------------------------------------
   ! Damage per mt5 event of the residuals without distribution
   ! (table of mt5tab, linear interpolation).
   !-------------------------------------------------------------------
   real(kr)::e
   mt5dm=t5int(e,d5t)
   return
   end function mt5dm

   real(kr) function mt5hm(e)
   !-------------------------------------------------------------------
   ! Kinetic energy per mt5 event of the residuals without
   ! distribution. (DLA, IAEA, Oct 2026)
   !-------------------------------------------------------------------
   real(kr)::e
   mt5hm=t5int(e,h5t)
   return
   end function mt5hm

   real(kr) function t5int(e,t)
   real(kr)::e,t(*)
   integer::j
   t5int=0
   if (n5t.le.0) return
   ! relative tolerance at the ends of the table (the pendf
   ! energy is rounded with sigfig and may exceed the last table energy
   ! by 1e-15, which gave zero at the last energy, e.g. 150 MeV)
   ! (DLA, IAEA, Oct 2026)
   if (e.lt.e5t(1)*(1-1.e-9_kr).or.e.gt.e5t(n5t)*(1+1.e-9_kr)) return
   if (e.ge.e5t(n5t)) then
      t5int=t(n5t)
      return
   endif
   if (e.le.e5t(1)) then
      t5int=t(1)
      return
   endif
   if (n5t.eq.1) then
      t5int=t(1)
      return
   endif
   do j=2,n5t
      if (e.le.e5t(j)) then
         t5int=t(j-1)+(t(j)-t(j-1))*(e-e5t(j-1))/(e5t(j)-e5t(j-1))
         return
      endif
   enddo
   return
   end function t5int

   subroutine fnadd(z,t,w)
   !-------------------------------------------------------------------
   ! Add one laboratory momentum node (pz along the
   ! incident direction, pt>=0 transverse modulus, weight w) to the
   ! fine-node buffer (fnz,fnt,fnw; nfn nodes).  Nodes with zero or
   ! negative weight are not added. (DLA, IAEA, Oct 2026)
   !-------------------------------------------------------------------
   ! externals
   real(kr)::z,t,w
   ! internals
   integer::nn
   real(kr),dimension(:),allocatable::tmp
   real(kr),parameter::zero=0
   if (w.le.zero) return
   if (.not.allocated(fnz)) then
      allocate(fnz(100000))
      allocate(fnt(100000))
      allocate(fnw(100000))
   endif
   if (nfn.ge.size(fnz)) then
      nn=2*size(fnz)
      allocate(tmp(nn))
      tmp(1:nfn)=fnz(1:nfn)
      call move_alloc(tmp,fnz)
      allocate(tmp(nn))
      tmp(1:nfn)=fnt(1:nfn)
      call move_alloc(tmp,fnt)
      allocate(tmp(nn))
      tmp(1:nfn)=fnw(1:nfn)
      call move_alloc(tmp,fnw)
   endif
   nfn=nfn+1
   fnz(nfn)=z
   fnt(nfn)=t
   fnw(nfn)=w
   return
   end subroutine fnadd

   subroutine cvcmp(n,z,t,w,m,cz,ct,cw)
   !-------------------------------------------------------------------
   ! Compression of n momentum nodes of an axially
   ! symmetric distribution (pz, pt>=0, weight w>=0) on a cvnb x cvnb
   ! grid of (pz,pt) bins spanning the nodes.  Each bin is replaced by
   ! two nodes (cvemit)
   !    pz = <pz> -+ sigma(pz),  pt = sqrt(<pt**2>),  weight W/2
   ! so that the weight, <pz>, <pz**2> and <pt**2> of every bin are
   ! preserved.  The mean recoil energy of a sum of independent
   ! momenta depends only on these moments and is therefore kept
   ! exactly; the shape changes only within the bins.  The m output
   ! nodes (m <= cvmx) are normalized to unit weight; m=0 if the total
   ! weight is zero. (DLA, IAEA, Oct 2026)
   !-------------------------------------------------------------------
   ! externals
   integer::n,m
   real(kr)::z(*),t(*),w(*),cz(*),ct(*),cw(*)
   ! internals
   integer::i,iz,it,k
   real(kr)::zmin,zmax,tmax,dz,dt,ws
   real(kr)::bw(cvnb*cvnb),bz(cvnb*cvnb),bz2(cvnb*cvnb),bt2(cvnb*cvnb)
   real(kr),parameter::zero=0
   m=0
   if (n.le.0) return
   zmin=z(1)
   zmax=z(1)
   tmax=t(1)
   do i=2,n
      zmin=min(zmin,z(i))
      zmax=max(zmax,z(i))
      tmax=max(tmax,t(i))
   enddo
   dz=zmax-zmin
   dt=tmax
   bw=0
   bz=0
   bz2=0
   bt2=0
   do i=1,n
      if (w(i).le.zero) cycle
      iz=1
      if (dz.gt.zero) iz=min(int(max(z(i)-zmin,zero)/dz*cvnb),cvnb-1)+1
      it=1
      if (dt.gt.zero) it=min(int(max(t(i),zero)/dt*cvnb),cvnb-1)+1
      k=(iz-1)*cvnb+it
      bw(k)=bw(k)+w(i)
      bz(k)=bz(k)+w(i)*z(i)
      bz2(k)=bz2(k)+w(i)*z(i)*z(i)
      bt2(k)=bt2(k)+w(i)*t(i)*t(i)
   enddo
   ws=0
   do k=1,cvnb*cvnb
      ws=ws+bw(k)
   enddo
   if (ws.le.zero) return
   call cvemit(bw,bz,bz2,bt2,ws,m,cz,ct,cw)
   return
   end subroutine cvcmp

   subroutine cvemit(bw,bz,bz2,bt2,ws,m,cz,ct,cw)
   !-------------------------------------------------------------------
   ! Nodes of the bins of cvcmp and cvdam: for each bin
   ! with weight W>0, <pz>-sigma and <pz>+sigma with W/2 each (one
   ! node with W if sigma=0) and pt=sqrt(<pt**2>); weights divided by
   ! the total ws. (DLA, IAEA, Oct 2026)
   !-------------------------------------------------------------------
   ! externals
   integer::m
   real(kr)::bw(*),bz(*),bz2(*),bt2(*),ws,cz(*),ct(*),cw(*)
   ! internals
   integer::k
   real(kr)::av,sg,tt
   real(kr),parameter::zero=0
   m=0
   do k=1,cvnb*cvnb
      if (bw(k).le.zero) cycle
      av=bz(k)/bw(k)
      sg=sqrt(max(bz2(k)/bw(k)-av*av,zero))
      tt=sqrt(max(bt2(k)/bw(k),zero))
      if (sg.gt.zero) then
         m=m+1
         cz(m)=av-sg
         ct(m)=tt
         cw(m)=bw(k)/(2*ws)
         m=m+1
         cz(m)=av+sg
         ct(m)=tt
         cw(m)=bw(k)/(2*ws)
      else
         m=m+1
         cz(m)=av
         ct(m)=tt
         cw(m)=bw(k)/ws
      endif
   enddo
   return
   end subroutine cvemit

   subroutine cvdam(nset,ns,io,cz,ct,cw,p0,zr,ar,zl,al,d,erm)
   !-------------------------------------------------------------------
   ! Mean damage energy of the residual nucleus after the
   ! independent emission of nset particles, each with an axially
   ! symmetric laboratory momentum distribution given by compressed
   ! nodes (set k: ns(k) nodes from position io(k)+1 of cz, ct, cw,
   ! unit total weight; the same set may be given several times for
   ! identical particles):
   !    er = |p0 z - sum_k p_k|**2/(2 ar),  d = <df(er;zr,ar,zl,al)>,
   !    erm = <er>,
   ! p0=sqrt(2 e), momenta in units sqrt(neutron mass * eV).
   ! The momenta are added one at a time.  For the sum S of the first
   ! k-1 momenta and p_k with relative azimuth phi,
   !    pz = Sz + pkz,  pt = sqrt(St**2 + pkt**2 + 2 St pkt cos(phi)),
   ! phi at the cvnphi midpoints of (0,pi) (equal weights; the
   ! integrand is periodic and even in phi, so the rule is exact for
   ! cos(n phi) up to n = 2 cvnphi-1).  For k < nset the sum is
   ! compressed again with the bins of cvcmp (moments kept); the last
   ! momentum is integrated directly.  nset=1 gives the exact
   ! one-particle value of the nodes. (DLA, IAEA, Oct 2026)
   !-------------------------------------------------------------------
   ! externals
   integer::nset,ns(*),io(*)
   real(kr)::cz(*),ct(*),cw(*),p0,zr,ar,zl,al,d,erm
   ! internals
   integer::i,j,l,k,n,m,iz,it,kb
   real(kr)::sz(cvmx),st(cvmx),sw(cvmx),cs(cvnphi)
   real(kr)::bw(cvnb*cvnb),bz(cvnb*cvnb),bz2(cvnb*cvnb),bt2(cvnb*cvnb)
   real(kr)::zmin,zmax,tmax,dz,dt,x,y,a,b,ww,er,dum,ws
   real(kr),parameter::zero=0
   real(kr),parameter::pi=3.14159265358979323846e0_kr
   d=0
   erm=0
   if (nset.le.0.or.ar.le.zero) return
   dum=df(zero,zr,ar,zl,al)
   if (nset.eq.1) then
      do i=io(1)+1,io(1)+ns(1)
         er=((p0-cz(i))**2+ct(i)**2)/(2*ar)
         d=d+cw(i)*df(er,zr,ar,zl,al)
         erm=erm+cw(i)*er
      enddo
      return
   endif
   do l=1,cvnphi
      cs(l)=cos((l-0.5e0_kr)*pi/cvnphi)
   enddo
   n=ns(1)
   do i=1,n
      sz(i)=cz(io(1)+i)
      st(i)=ct(io(1)+i)
      sw(i)=cw(io(1)+i)
   enddo
   do k=2,nset-1
      zmin=sz(1)
      zmax=sz(1)
      tmax=st(1)
      do i=2,n
         zmin=min(zmin,sz(i))
         zmax=max(zmax,sz(i))
         tmax=max(tmax,st(i))
      enddo
      a=cz(io(k)+1)
      b=cz(io(k)+1)
      y=ct(io(k)+1)
      do j=io(k)+2,io(k)+ns(k)
         a=min(a,cz(j))
         b=max(b,cz(j))
         y=max(y,ct(j))
      enddo
      zmin=zmin+a
      zmax=zmax+b
      tmax=tmax+y
      dz=zmax-zmin
      dt=tmax
      bw=0
      bz=0
      bz2=0
      bt2=0
      do i=1,n
         do j=io(k)+1,io(k)+ns(k)
            ww=sw(i)*cw(j)/cvnphi
            x=sz(i)+cz(j)
            iz=1
            if (dz.gt.zero) iz=min(int(max(x-zmin,zero)/dz*cvnb),cvnb-1)+1
            a=st(i)**2+ct(j)**2
            b=2*st(i)*ct(j)
            do l=1,cvnphi
               y=sqrt(max(a+b*cs(l),zero))
               it=1
               if (dt.gt.zero) it=min(int(y/dt*cvnb),cvnb-1)+1
               kb=(iz-1)*cvnb+it
               bw(kb)=bw(kb)+ww
               bz(kb)=bz(kb)+ww*x
               bz2(kb)=bz2(kb)+ww*x*x
               bt2(kb)=bt2(kb)+ww*y*y
            enddo
         enddo
      enddo
      ws=0
      do kb=1,cvnb*cvnb
         ws=ws+bw(kb)
      enddo
      if (ws.le.zero) return
      call cvemit(bw,bz,bz2,bt2,ws,m,sz,st,sw)
      n=m
   enddo
   k=nset
   do i=1,n
      do j=io(k)+1,io(k)+ns(k)
         x=p0-sz(i)-cz(j)
         a=x*x+st(i)**2+ct(j)**2
         b=2*st(i)*ct(j)
         ww=sw(i)*cw(j)/cvnphi
         do l=1,cvnphi
            er=(a+b*cs(l))/(2*ar)
            if (er.lt.zero) er=0
            d=d+ww*df(er,zr,ar,zl,al)
            erm=erm+ww*er
         enddo
      enddo
   enddo
   return
   end subroutine cvdam

   subroutine dirnod(e,ep,awr,wt)
   !-------------------------------------------------------------------
   ! Laboratory momentum nodes of one neutron emitted with
   ! laboratory energy ep (File 5 energies are laboratory energies) at
   ! incident energy e, with the File 4 angular distribution stored by
   ! anitab, added to the fine-node buffer with total weight wt.
   !   no File 4 (n4e=0): isotropic in the laboratory (as tdfn).
   !   LCT=1: pz=v mu, pt=v sqrt(1-mu**2) at the nodes mu4 (weights
   !      w4), v=sqrt(2 ep).
   !   LCT=2: cm cosines mu_c at fixed laboratory speed v (as anipz),
   !      V=sqrt(2e)/(awr+1): for v >= V, u=-V mu_c+sqrt(v**2-V**2
   !      (1-mu_c**2)), pz=V+u mu_c, pt=u sqrt(1-mu_c**2), weights w4;
   !      for v < V, the two branches u=-V mu_c -+ sqrt(...) on
   !      mu_c < -sqrt(1-v**2/V**2) (64-point rule), renormalized.
   ! Between the incident energies of File 4 the node sets of the two
   ! tables are mixed with the weight of the File 4 law in e (law 1:
   ! lower table; laws 2 and 4: linear in e; laws 3 and 5: linear in
   ! ln e), i.e. the distribution is interpolated (anipz interpolated
   ! <pz> with the law).  Below the first energy the distribution is
   ! isotropic in the frame of the data, as in anipz.  Negative values
   ! of a Legendre distribution are set to zero (the weights of each
   ! table are renormalized to unit area). (DLA, IAEA, Oct 2026)
   !-------------------------------------------------------------------
   ! externals
   real(kr)::e,ep,awr,wt
   ! internals
   integer::k,nk,ik,lw,kk(2),iq
   real(kr)::wk(2),v,vl,t
   real(kr),parameter::zero=0
   real(kr),parameter::small=1.e-10_kr
   if (wt.le.zero) return
   vl=sqrt(2*max(ep,zero))
   if (vl.le.zero) then
      call fnadd(zero,zero,wt)
      return
   endif
   v=sqrt(2*e)/(awr+1)
   if (n4e.le.0) then
      call isolab
      return
   endif
   nk=1
   wk(1)=1
   if (iso4.eq.1) then
      kk(1)=1
   else if (e.lt.e4(1)*(1-small)) then
      if (lct4.ne.2) then
         call isolab
         return
      endif
      kk(1)=0
   else if (e.ge.e4(n4e)) then
      kk(1)=n4e
   else
      kk(1)=n4e
      do k=1,n4e-1
         if (e.ge.e4(k).and.e.lt.e4(k+1)) then
            lw=law4(k)
            t=(e-e4(k))/(e4(k+1)-e4(k))
            if ((lw.eq.3.or.lw.eq.5).and.e4(k).gt.zero)&
              t=log(e/e4(k))/log(e4(k+1)/e4(k))
            if (lw.eq.1) t=0
            kk(1)=k
            wk(1)=1-t
            if (t.gt.zero) then
               nk=2
               kk(2)=k+1
               wk(2)=t
            endif
            exit
         endif
      enddo
   endif
   do ik=1,nk
      if (wk(ik).gt.zero) call add1(kk(ik),wt*wk(ik))
   enddo
   return
   contains
      subroutine isolab
      ! isotropic in the laboratory (64-point rule)
      real(kr)::mu
      do iq=1,64
         mu=g64x(iq)
         call fnadd(vl*mu,vl*sqrt(max(1-mu*mu,zero)),wt*g64w(iq)/2)
      enddo
      end subroutine isolab
      subroutine add1(kq,wa)
      ! nodes of the table kq (0: isotropic cm) with total weight wa
      integer::kq,is
      real(kr)::wa,mu,u,r,ms,a,b,w,g,sp
      if (lct4.eq.1) then
         sp=0
         do iq=k4o(kq),k4o(kq)+k4n(kq)-1
            if (w4(iq).gt.zero) sp=sp+w4(iq)
         enddo
         if (sp.le.zero) return
         do iq=k4o(kq),k4o(kq)+k4n(kq)-1
            mu=mu4(iq)
            call fnadd(vl*mu,vl*sqrt(max(1-mu*mu,zero)),wa*w4(iq)/sp)
         enddo
      else if (vl.ge.v) then
         if (kq.eq.0) then
            do iq=1,64
               mu=g64x(iq)
               u=-v*mu+sqrt(max(vl*vl-v*v*(1-mu*mu),zero))
               call fnadd(v+u*mu,u*sqrt(max(1-mu*mu,zero)),wa*g64w(iq)/2)
            enddo
         else
            sp=0
            do iq=k4o(kq),k4o(kq)+k4n(kq)-1
               if (w4(iq).gt.zero) sp=sp+w4(iq)
            enddo
            if (sp.le.zero) return
            do iq=k4o(kq),k4o(kq)+k4n(kq)-1
               mu=mu4(iq)
               u=-v*mu+sqrt(max(vl*vl-v*v*(1-mu*mu),zero))
               call fnadd(v+u*mu,u*sqrt(max(1-mu*mu,zero)),wa*w4(iq)/sp)
            enddo
         endif
      else
         ! v < V: backward cone, two branches, renormalized
         ms=-sqrt(1-(vl/v)**2)
         a=(ms-1)/2
         b=(ms+1)/2
         sp=0
         do iq=1,64
            mu=b*g64x(iq)+a
            g=max(gcm4(kq,mu),zero)
            sp=sp+2*b*g64w(iq)*g
         enddo
         if (sp.le.zero) return
         do iq=1,64
            mu=b*g64x(iq)+a
            w=b*g64w(iq)*max(gcm4(kq,mu),zero)/sp
            r=sqrt(max(vl*vl-v*v*(1-mu*mu),zero))
            do is=1,2
               u=-v*mu+r
               if (is.eq.2) u=-v*mu-r
               call fnadd(v+u*mu,abs(u)*sqrt(max(1-mu*mu,zero)),wa*w)
            enddo
         enddo
      endif
      end subroutine add1
   end subroutine dirnod

   real(kr) function gcm4(kk,mu)
   !-------------------------------------------------------------------
   ! Normalized File 4 distribution stored by anitab at energy index kk
   ! (0: isotropic), at the cosine mu (moved out of anipz,
   ! where it was the internal function gcm, to be used by dirnod).
   !-------------------------------------------------------------------
   use mathm ! provides legndr
   use endf  ! provides terp1
   ! externals
   integer::kk
   real(kr)::mu
   ! internals
   integer::il,j,lw,jo
   real(kr)::p(66),x1,x2,y1,y2,y
   real(kr),parameter::zero=0
   gcm4=0.5_kr
   if (kk.eq.0) return
   if (kind4(kk).eq.0) return
   if (kind4(kk).eq.1) then
      call legndr(mu,p,max(r4n(kk),1))
      gcm4=0.5_kr
      do il=1,r4n(kk)
         gcm4=gcm4+(2*il+1)*rp4(r4o(kk)+il-1)*p(il+1)/2
      enddo
   else
      gcm4=0
      jo=r4o(kk)
      do j=2,r4n(kk)
         x1=rp4(jo+2*(j-2))
         y1=rp4(jo+2*(j-2)+1)
         x2=rp4(jo+2*(j-1))
         y2=rp4(jo+2*(j-1)+1)
         if (mu.ge.x1.and.mu.le.x2.and.x2.gt.x1) then
            lw=lp4((jo+2*(j-1)+1)/2)
            if ((lw.eq.4.or.lw.eq.5).and.(y1.le.zero.or.y2.le.zero)) lw=2
            if ((lw.eq.3.or.lw.eq.5).and.x1.le.zero) lw=2
            call terp1(x1,y1,x2,y2,mu,y,lw)
            gcm4=y
            exit
         endif
      enddo
   endif
   if (s4(kk).gt.zero) gcm4=gcm4/s4(kk)
   return
   end function gcm4

   subroutine anitab(nin,matd,mtd)
   !-------------------------------------------------------------------
   ! Read the File 4 angular distribution of reaction mtd (neutron-
   ! only reaction given with Files 4 and 5) and store, for each
   ! incident energy, a quadrature of the normalized distribution
   ! in mu (nodes mu4, weights w4 = quadrature weight * f(mu)) and
   ! its mean cosine f14, in the frame of the data (lct4).
   !   LTT=1 (Legendre): 64-point Gauss rule on (-1,1); the
   !         distribution f(mu)=1/2+sum (2l+1)/2 a_l P_l(mu) is
   !         used as given (no truncation of negative values), so that
   !         f14 = a_1 exactly.
   !   LTT=2 (tabulated): 4-point Gauss rule in each panel of the
   !         table, f(mu) from the panel interpolation law (exact
   !         moments for histogram and linear panels).
   !   LTT=3: both representations, Legendre below and tables above
   !         the transition energy.
   !   LTT=0 or LI=1: isotropic; nothing is stored for LCT=1 (no
   !         correction); for LCT=2 one isotropic set is stored.
   ! The incident-energy interpolation law of the TAB2 records is kept
   ! (law4).  Below the first incident energy the distribution is
   ! isotropic in the frame of the data (as in hgtfle).
   !-------------------------------------------------------------------
   use endf ! provides findf,contio,listio,tab1io,tab2io,moreio,terp1
   use util ! provides error
   use mathm ! provides legndr
   ! externals
   integer::nin,matd,mtd
   ! internals
   integer::i,found,nb,nw,lvt,ltt,li,nseg,iseg,kind,ne,ie,l,nl,nr1,np1
   integer::ib,iq,j,lw,il,nr2,kk,ii,krp
   integer,parameter::nbmax=400000,nemax=5000,nkmax=500000,nrmx=50
   integer::nbt2(nrmx),int2(nrmx)
   real(kr),dimension(:),allocatable::b
   real(kr)::p(66),en,x1,x2,y1,y2,xm,ym,s,pdf
   real(kr),dimension(4),parameter::q4x=(/&
      -.861136312e0_kr,-.339981044e0_kr,.339981044e0_kr,.861136312e0_kr/)
   real(kr),dimension(4),parameter::q4w=(/&
      .347854845e0_kr,.652145155e0_kr,.652145155e0_kr,.347854845e0_kr/)
   real(kr),parameter::zero=0

   n4e=0
   iso4=0
   lct4=1
   found=0
   do i=1,n4mt
      if (mt4l(i).eq.mtd) found=1
   enddo
   if (found.eq.0) return
   if (.not.allocated(e4)) then
      allocate(e4(nemax),f14(nemax),law4(nemax),k4o(nemax),k4n(nemax))
      allocate(mu4(nkmax),w4(nkmax))
      allocate(s4(nemax),kind4(nemax),r4o(nemax),r4n(nemax))
      allocate(rp4(nkmax),lp4(nkmax))
   endif
   allocate(b(nbmax))
   call findf(matd,4,mtd,nin)
   call contio(nin,0,0,b,nb,nw)
   lvt=nint(b(3))
   ltt=nint(b(4))
   if (lvt.eq.0) then
      call contio(nin,0,0,b,nb,nw)
   else
      call listio(nin,0,0,b,nb,nw)
      l=1
      do while (nb.ne.0)
         l=l+nw
         call moreio(nin,0,0,b(l),nb,nw)
      enddo
   endif
   li=nint(b(3))
   lct4=nint(b(4))
   if (lct4.ne.1) lct4=2
   kk=0
   krp=0
   ! isotropic distribution
   if (ltt.eq.0.or.li.eq.1) then
      if (lct4.eq.2) then
         iso4=1
         n4e=1
         e4(1)=0
         law4(1)=2
         k4o(1)=1
         k4n(1)=64
         kind4(1)=0
         s4(1)=1
         do iq=1,64
            mu4(iq)=g64x(iq)
            w4(iq)=g64w(iq)/2
         enddo
         f14(1)=0
      endif
      deallocate(b)
      return
   endif
   nseg=1
   if (ltt.eq.3) nseg=2
   do iseg=1,nseg
      kind=ltt
      if (ltt.eq.3) kind=iseg
      call tab2io(nin,0,0,b,nb,nw)
      nr2=nint(b(5))
      ne=nint(b(6))
      if (nr2.gt.nrmx) call error('anitab','too many interpolation ranges',' ')
      do j=1,nr2
         nbt2(j)=nint(b(5+2*j))
         int2(j)=nint(b(6+2*j))
      enddo
      do ie=1,ne
         l=1
         if (kind.eq.1) then
            call listio(nin,0,0,b(l),nb,nw)
         else
            call tab1io(nin,0,0,b(l),nb,nw)
         endif
         do while (nb.ne.0)
            l=l+nw
            if (l.gt.nbmax-1000) call error('anitab','storage exceeded',' ')
            call moreio(nin,0,0,b(l),nb,nw)
         enddo
         en=b(2)
         n4e=n4e+1
         if (n4e.gt.nemax) call error('anitab','too many energies',' ')
         e4(n4e)=en
         ! interpolation law from this energy to the next one
         lw=int2(nr2)
         do j=nr2,1,-1
            if (ie+1.le.nbt2(j)) lw=int2(j)
         enddo
         law4(n4e)=lw
         k4o(n4e)=kk+1
         if (kind.eq.1) then
            ! Legendre coefficients a_1..a_nl in b(7..6+nl)
            nl=min(nint(b(5)),64)
            kind4(n4e)=1
            r4o(n4e)=krp+1
            r4n(n4e)=nl
            do il=1,nl
               krp=krp+1
               if (krp.gt.nkmax) call error('anitab','too many coefficients',' ')
               rp4(krp)=b(6+il)
            enddo
            do iq=1,64
               call legndr(g64x(iq),p,max(nl,1))
               pdf=0.5_kr
               do il=1,nl
                  pdf=pdf+(2*il+1)*b(6+il)*p(il+1)/2
               enddo
               kk=kk+1
               if (kk.gt.nkmax) call error('anitab','too many nodes',' ')
               mu4(kk)=g64x(iq)
               w4(kk)=g64w(iq)*pdf
            enddo
         else
            ! tabulated f(mu): panels with their interpolation law
            nr1=nint(b(5))
            np1=nint(b(6))
            ib=7+2*nr1
            kind4(n4e)=2
            r4o(n4e)=krp+1
            r4n(n4e)=np1
            ii=1
            do j=1,np1
               do while (j.gt.nint(b(5+2*ii)).and.ii.lt.nr1)
                  ii=ii+1
               enddo
               krp=krp+2
               if (krp.gt.nkmax) call error('anitab','too many points',' ')
               rp4(krp-1)=b(ib+2*(j-1))
               rp4(krp)=b(ib+2*(j-1)+1)
               ! law of the panel (j-1,j) stored at point j
               lp4(krp/2)=nint(b(6+2*ii))
            enddo
            ii=1
            do j=1,np1-1
               do while (j+1.gt.nint(b(5+2*ii)).and.ii.lt.nr1)
                  ii=ii+1
               enddo
               lw=nint(b(6+2*ii))
               x1=b(ib+2*(j-1))
               y1=b(ib+2*(j-1)+1)
               x2=b(ib+2*j)
               y2=b(ib+2*j+1)
               if (x2.le.x1) cycle
               if ((lw.eq.4.or.lw.eq.5).and.(y1.le.zero.or.y2.le.zero)) lw=2
               if ((lw.eq.3.or.lw.eq.5).and.x1.le.zero) lw=2
               do iq=1,4
                  xm=(x1+x2)/2+(x2-x1)/2*q4x(iq)
                  call terp1(x1,y1,x2,y2,xm,ym,lw)
                  kk=kk+1
                  if (kk.gt.nkmax) call error('anitab','too many nodes',' ')
                  mu4(kk)=xm
                  w4(kk)=q4w(iq)*(x2-x1)/2*ym
               enddo
            enddo
         endif
         k4n(n4e)=kk-k4o(n4e)+1
         ! normalize and mean cosine
         s=0
         do iq=k4o(n4e),kk
            s=s+w4(iq)
         enddo
         s4(n4e)=s
         f14(n4e)=0
         if (s.gt.zero) then
            do iq=k4o(n4e),kk
               w4(iq)=w4(iq)/s
               f14(n4e)=f14(n4e)+w4(iq)*mu4(iq)
            enddo
         endif
      enddo
   enddo
   deallocate(b)
   return
   end subroutine anitab

   real(kr) function anipz(e,ep,awr)
   !-------------------------------------------------------------------
   ! Mean laboratory momentum along the incident direction of one
   ! neutron emitted with laboratory energy ep (File 5 energies are
   ! always laboratory energies) and the File 4 angular distribution
   ! stored by anitab, for incident energy e (neutron mass units).
   !   LCT=1: pz = sqrt(2 ep) <mu>.
   !   LCT=2: the File 4 distribution g(mu_c) is in the centre of mass.
   !   At fixed laboratory speed v=sqrt(2 ep) the laboratory cosine
   !   mu_L and the CM cosine are related by the velocity addition
   !   v n_L = V z + u n_c (V=sqrt(2e)/(awr+1)):
   !      u = sqrt(v**2+V**2-2 v V mu_L), mu_c = (v mu_L - V)/u,
   !   and the laboratory distribution is g_L(mu_L) = g(mu_c) J,
   !   J = d mu_c/d mu_L = v**2 (v - V mu_L)/u**3 (Jacobian at fixed
   !   laboratory energy).  <pz> = v int mu_L g_L d mu_L / int g_L.
   !   v >= V: mu_c covers (-1,1) once and the integral is done in the
   !      variable mu_c (same result, no singularity at v=V):
   !      <pz> = int g(mu_c) (V + u mu_c) d mu_c,
   !      u = -V mu_c + sqrt(v**2 - V**2 (1-mu_c**2)).
   !   v < V: only mu_c < -sqrt(1-v**2/V**2) is reachable, twice (two
   !      laboratory directions); g_L is renormalized to unit area:
   !      <pz> = sum(+,-) int g (V + u(+,-) mu_c) d mu_c /
   !             (2 int g d mu_c), on (-1, -sqrt(1-v**2/V**2)).
   ! The result is interpolated between incident energies with the
   ! File 4 law (moments are linear in the distribution).
   !-------------------------------------------------------------------
   use endf ! provides terp1
   use mathm ! provides legndr
   ! externals
   real(kr)::e,ep,awr
   ! internals
   integer::k,lw
   real(kr)::v,vl,p1,p2
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::zero=0
   anipz=0
   if (n4e.le.0) return
   vl=sqrt(2*max(ep,zero))
   v=sqrt(2*e)/(awr+1)
   if (iso4.eq.1) then
      anipz=pzk(1)
      return
   endif
   if (e.lt.e4(1)*(1-small)) then
      if (lct4.eq.2) anipz=pzk(0)
      return
   endif
   if (e.ge.e4(n4e)) then
      anipz=pzk(n4e)
      return
   endif
   do k=1,n4e-1
      if (e.ge.e4(k).and.e.lt.e4(k+1)) then
         p1=pzk(k)
         p2=pzk(k+1)
         lw=law4(k)
         if ((lw.eq.4.or.lw.eq.5).and.p1*p2.le.zero) lw=lw-2
         call terp1(e4(k),p1,e4(k+1),p2,e,anipz,lw)
         return
      endif
   enddo
   anipz=pzk(n4e)
   return
   contains
      real(kr) function pzk(kk)
      ! kk=0: isotropic CM distribution (below the first energy)
      integer::kk,iq,is
      real(kr)::u,mu,r,ms,a,b,w,g,sn,sg
      if (lct4.eq.1) then
         pzk=vl*f14(kk)
         return
      endif
      pzk=0
      if (vl.ge.v) then
         if (kk.eq.0) then
            do iq=1,64
               mu=g64x(iq)
               u=-v*mu+sqrt(max(vl*vl-v*v*(1-mu*mu),zero))
               pzk=pzk+g64w(iq)/2*(v+u*mu)
            enddo
         else
            do iq=k4o(kk),k4o(kk)+k4n(kk)-1
               mu=mu4(iq)
               u=-v*mu+sqrt(max(vl*vl-v*v*(1-mu*mu),zero))
               pzk=pzk+w4(iq)*(v+u*mu)
            enddo
         endif
      else
         ! v < V: backward cone, two branches, renormalized
         if (vl.le.zero) return
         ms=-sqrt(1-(vl/v)**2)
         a=(ms-1)/2
         b=(ms+1)/2
         sn=0
         sg=0
         do iq=1,64
            mu=b*g64x(iq)+a
            w=b*g64w(iq)
            g=gcm4(kk,mu)
            r=sqrt(max(vl*vl-v*v*(1-mu*mu),zero))
            do is=1,2
               u=-v*mu+r
               if (is.eq.2) u=-v*mu-r
               sn=sn+w*g*(v+u*mu)
               sg=sg+w*g
            enddo
         enddo
         if (sg.gt.zero) pzk=sn/sg
      endif
      end function pzk
   end function anipz

   subroutine sixbar(e,ebar,yld,dame,nin,c,ncmax,nscr,b,nbmax,&
     n2,j6,irec,jrec,iflag)
   !-------------------------------------------------------------------
   ! Calculate charged-particle ebar and damage energy
   ! for reactions represented using ENDF-6 File 6.
   ! If recoil subsection is not given, attempt to compute
   ! from the distribution for the heavy emitted particle.
   ! For normal particles, irec is zero.  For recoil cases,
   ! irec is the index of the section for the driving particle.
   !-------------------------------------------------------------------
   use util ! provides skiprz,mess
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,n2,j6,irec,jrec,ncmax,nbmax,nscr,iflag
   real(kr)::e,ebar,yld,dame
   real(kr)::c(ncmax),b(nbmax)
   ! internals
   integer::matd,mfd,mtd,l,nb,nw,ik,nnt,nmu,imu,mf1,mt1,idis
   integer::ip,ir,idisc,iraw,nne,ne,law,lang,lep,intl
   real(kr)::disc102,zp,zt,ap,at,ztt,pe,eihi,f,d,s,enext
   real(kr)::elo,ehi,flo,fhi,dlo,dhi
   real(kr)::yph,en102
   character(60)::strng
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::zero=0
   save iraw,nne,ne,law,lang,lep,intl
   save elo,ehi,flo,fhi,dlo,dhi
   save disc102,zp,ap,zt,at
   ! done6 is false while tables of the current subsection
   ! remain unread (energy grid of the pendf tape ends below the last
   ! incident energy of the subsection) (DLA, IAEA, Oct 2026)
   logical,save::done6=.true.

   !--initialize for this subsection when e=0.
   if (e.gt.zero) go to 300
   if (irec.eq.0) then
      if (.not.done6.and.j6.gt.1) then
         ! the previous subsection was not read to its
         ! end; position the tape at the start of subsection j6 (the
         ! NJOY2016.80 read the next subsection from the wrong place:
         ! zap=0, law=0 and a fortran read error, e.g. JENDL-5 N-14)
         ! (DLA, IAEA, Oct 2026)
         call findf(math,6,mth,nin)
         l=1
         call contio(nin,0,0,c(l),nb,nw)
         do ik=1,j6-1
            call tab1io(nin,0,0,c(l),nb,nw)
            do while (nb.ne.0)
               call moreio(nin,0,0,c(l),nb,nw)
            enddo
            law=l2h
            call skip6(nin,0,0,c(l),law)
         enddo
      endif
      done6=.true.
      go to 110
   endif
   done6=.true.
   if (jrec.gt.0) irec=jrec

   !--for recoils, back up to the driving particle,
   !--or restore position after doing recoil.
  100 continue
   matd=math
   mfd=6
   mtd=mth
   call skiprz(nin,-2)
   call findf(matd,mfd,mtd,nin)
   l=1
   call contio(nin,0,0,c(l),nb,nw)
   ik=0
   do while (ik.lt.irec-1)
      ik=ik+1
      call tab1io(nin,0,0,c(l),nb,nw)
      l=l+nw
      do while (nb.ne.0)
         call moreio(nin,0,0,c(l),nb,nw)
         l=l+nw
      enddo
      law=l2h
      call skip6(nin,0,0,c(l),law)
   enddo
   if (jrec.eq.irec) then
      irec=0
      jrec=0
   endif

   !--start reading in the desired subsection
  110 continue
   l=1
   call tab1io(nin,0,0,c(l),nb,nw)
   l=l+nw
   do while (nb.ne.0)
      call moreio(nin,0,0,c(l),nb,nw)
      l=l+nw
   enddo
   zap=c1h
   awp=c2h
   law=l2h

   !--error for awp=0. for non photons
   if (zap.ne.0.and.awp.eq.0) then
     write(strng,'(''awp is zero for for mt'',i3,'' particle '',i5)')&
       mth,nint(zap)
     call error('sixbar',strng,' ')
   endif

   iflag=0
   disc102=0
   if (zap.eq.zero) then
      iflag=1
      disc102=awp
      awp=0
   endif
   ebar=0
   yld=0

   !--File 6 capture recoil rejected by the kinematic check
   !--in hinit; the photon-momentum recoil (mt6no) is used instead
   ! (DLA, IAEA, Oct 2026)
   if (mth.eq.102.and.irec.eq.0.and.zap.gt.zero.and.i6g.gt.0) then
      if (mt6bad(i6g).eq.j6) then
         call skip6(nin,0,0,c(l),law)
         ebar=-2
         return
      endif
   endif
   if (law.gt.0) go to 150

   !--warning for law=0.
   write(strng,'(''no distribution for mt'',i3,'' particle '',i5)')&
     mth,nint(zap)
   call mess('sixbar',strng,' ')
   ebar=-2
   return

   !--watch for 2-body recoils
  150 continue
   if (law.ne.4) go to 210
   irec=j6-1
   jrec=j6+1
   go to 100

   !--continue processing this subsection
  210 continue
   zp=int(zap/1000)
   zt=int(zat/1000)
   if (irec.gt.0) zp=zt-zp
   ap=awp
   if (irec.gt.0) ap=awr+1-awp
   at=awr
   if (zap.eq.zero) then
       ap=awr+1
       zp=zt
   endif
   dame=df(e,zp,ap,zt,at)
   if (disc102.gt.zero) go to 295
   if (law.ne.3.and.law.ne.6) then
      call tab2io(nin,0,0,c(l),nb,nw)
      lang=nint(c(l+2))
      lep=nint(c(l+3))
      ne=n2h
      nnt=l+6
      intl=nint(c(nnt+1))
      intl=2
      l=l+nw
   endif
   nne=0
   iraw=l
   if (law.eq.3) go to 230
   if (law.eq.6) go to 240
   if (law.eq.7) go to 250

   !--laws 1, 2, and 5.
   call listio(nin,0,0,c(l),nb,nw)
   l=l+nw
   do while (nb.ne.0)
      call moreio(nin,0,0,c(l),nb,nw)
      l=l+nw
   enddo
   nne=nne+1
   elo=c(iraw+1)
   ehi=0
   dhi=0
   fhi=0
   go to 290

   !--law 3 -- isotropic distribution
  230 continue
   ehi=emax
   lang=0
   return

   !--law 6 -- phase-space distribution
  240 continue
   call contio(nin,0,0,c(l),nb,nw)
   lang=0
   lep=2
   l=l+nw
   ehi=emax
   dhi=0
   fhi=0
   return

   !--law 7 -- laboratory angle-energy law
  250 continue
   call tab2io(nin,0,0,c(l),nb,nw)
   nmu=n2h
   l=l+nw
   do imu=1,nmu
      call tab1io(nin,0,0,c(l),nb,nw)
      l=l+nw
      do while (nb.ne.0)
         call moreio(nin,0,0,c(l),nb,nw)
         l=l+nw
      enddo
   enddo
   nne=nne+1
   elo=c(iraw+1)
   ehi=0
   dhi=0
   fhi=0

   !--compute low values of ebar and dame
  290 continue
   done6=(nne.ge.ne)
   if (mth.ne.102.or.irec.eq.0) then
      call getsix(elo,flo,dlo,c(iraw),law,lang,lep,irec)
   else
      ztt=int(zat/1000)
      call tabsq6(flo,dlo,c(iraw),law,ztt,awrt,elo,c(1))
   endif
   ! subsection with only one incident energy (NE=1); the
   ! NJOY2016.80 read the next record as a second table.  The single
   ! table is used at all incident energies. (DLA, IAEA, Oct 2026)
   if (ne.eq.1) then
      write(strng,'(''mt'',i3,'' particle '',i6,'' has NE=1'')')&
        mth,nint(zap)
      call mess('sixbar',strng,'the only table is used at all energies')
      ehi=emax
      fhi=flo
      dhi=dlo
   endif
   if ((mth.ge.18.and.mth.le.21).or.mth.eq.38) then
      matd=math
      mf1=1
      mt1=452
      if (nply.gt.0) mt1=456
      ! delayed neutrons are treated the same as fast neutrons.
      ! this will cause a slight overestimate of en for fission.
      ! ... but when nply.ne.0, use prompt neutron data only.
      call hgtyld(e,enext,idis,yld,matd,mf1,mt1,nscr,b,nbmax)
   endif
   return

  295 continue
   yld=1
   call skip6(nin,0,0,c(l),law)
   return

   !--normal entry
  300 continue
      if (disc102.gt.zero) go to 430

   !--interpolate for particle yield
   ip=2
   ir=1
   call terpa(pe,e,eihi,idisc,c(1),ip,ir)
   if (irec.gt.0) pe=1

   !--is desired energy in current panel
  305 continue
   if (e.lt.ehi*(1-small)) go to 400
   if (nne.eq.1) go to 310
   if (nne.eq.ne.and.e.le.ehi*(1+small)) go to 400
   if (nne.eq.ne) go to 450

   !--no. slide high-energy data into low-energy positions
   !--and read in new high energy data.
   elo=ehi
   flo=fhi
   dlo=dhi

   !--read in new data for high energy
   !--laws 1, 2, and 5.
  310 continue
   if (law.eq.1.or.law.eq.2.or.law.eq.5) then
      l=iraw
      call listio(nin,0,0,c(l),nb,nw)
      l=l+nw
      do while (nb.ne.0)
         call moreio(nin,0,0,c(l),nb,nw)
         l=l+nw
      enddo
      nne=nne+1
      ehi=c(iraw+1)

   !--law 6.
   else if (law.eq.6) then
      l=iraw
         call contio(nin,0,0,c(l),nb,nw)

   !--law 7.
   else if (law.eq.7) then
      l=iraw
      call tab2io(nin,0,0,c(l),nb,nw)
      nmu=n2h
      l=l+nw
      do imu=1,nmu
         call tab1io(nin,0,0,c(l),nb,nw)
         l=l+nw
         do while (nb.ne.0)
            call moreio(nin,0,0,c(l),nb,nw)
            l=l+nw
         enddo
      enddo
      ehi=c(iraw+1)
      nne=nne+1
   endif
   done6=(nne.ge.ne)

   !--compute high values of ebar and dame
   if (mth.ne.102.or.irec.eq.0) then
      call getsix(ehi,fhi,dhi,c(iraw),law,lang,lep,irec)
   else
      ztt=int(zat/1000)
      call tabsq6(fhi,dhi,c(iraw),law,ztt,awrt,ehi,c(1))
   endif
   go to 305

   !--yes.  compute result.
  400 continue
   if (law.ne.3.and.law.ne.6) go to 420
   c(iraw+1)=e
   call getsix(e,f,d,c(iraw),law,lang,lep,irec)
   ebar=f
   yld=pe
   dame=pe*d
   return

   !--yes. interpolate for results
  420 continue
   if (nne.eq.1.and.ne.gt.1) go to 450
   ! below the first incident energy of the subsection the
   ! values of the first table are used (NJOY2016.80 extrapolated
   ! linearly, with possible negative values) (DLA, IAEA, Oct 2026)
   if (e.lt.elo) then
      s=flo
   else
      call terp1(elo,flo,ehi,fhi,e,s,intl)
   endif
   ebar=s
   yld=pe
   if ((mth.ge.18.and.mth.le.21).or.mth.eq.38) then
      matd=math
      mf1=1
      mt1=452
      call hgtyld(e,enext,idis,yld,matd,mf1,mt1,nscr,b,nbmax)
   endif
   if (e.lt.elo) then
      d=dlo
   else
      call terp1(elo,dlo,ehi,dhi,e,d,intl)
   endif
   if (mth.eq.102.and.irec.gt.0) then
      ! capture recoil from the photon momentum (File 6):
      ! add the neutron kick at the actual incident energy (manual
      ! eq. 210-211); s=yield*<eg**2>/(2Mc**2),
      ! d=yield*<D(eg**2/(2Mc**2))>; one photon: D(neutron + photon kick)
      ! (DLA, IAEA, Oct 2026)
      ip=2
      ir=1
      call terpa(yph,e,enext,idisc,c(1),ip,ir)
      en102=e/(awr+1)
      if (yph.le.ysing) then
         d=df(en102+s,zt,awr+1,zt,awr)
      else
         d=df(en102,zt,awr+1,zt,awr)+d
      endif
      ! the recoil energy (heating) also includes the
      ! neutron kick, <er> = e/(a+1) + yield*<eg**2>/(2Mc**2)
      ! (DLA, IAEA, Oct 2026)
      ebar=en102+s
   endif
   dame=pe*d
   return

   !--discrete relativistic capture gamma
  430 call hgam102(e,ebar,dame,disc102,c,irec,zp,ap,zt,at)
   yld=1
   return

   !--return zeros outside range of table
  450 continue
   ebar=0
   dame=0
   return
   end subroutine sixbar

   subroutine getsix(e,ebar,dame,c,law,lang,lep,irec)
   !-------------------------------------------------------------------
   ! Compute ebar and dame for one particular incident energy
   ! using data in ENDF-6 File 6 format.  For irec.eq.0, return the
   ! numbers for the emitted particle.  For irec.gt.0, return the
   ! numbers for the recoil due to a single emission of the particle.
   !-------------------------------------------------------------------
   use mathm ! provides legndr
   use endf ! provides terp1 (DLA, IAEA, 19 Oct 2026)
   ! externals
   integer::law,lang,lep,irec
   real(kr)::e,ebar,dame,c(*)
   ! internals
   integer::nl,l,i,nd,na,nep,ncyc,iq,ia,nld,il,nmu,imu
   integer::nr,np,next,ibase,iint
   real(kr)::zp,zt,ap,at,summ,epnext,el,fl,dy,da,xm
   integer::ilab6
   real(kr)::epn,ym,test,en,f,f1,h,d,xl,yl,xx,yy,xxe
   real(kr)::x2,zpp,ztt,b,er,t1,t2,thresh,beta
   real(kr)::afact,arec,u,e2,ul,fn,dl,hl,f2
   integer::lt6,it6,nsb,isb
   real(kr)::xa,ya,xb,yb,sa,sb,fa,fb
   real(kr)::x(10),y(10,2)
   real(kr)::term(2)
   real(kr)::p(65)
   integer,parameter::nq=64
   real(kr),dimension(64),parameter::qp=(/&
     -9.99305042E-01_kr,-9.96340117E-01_kr,-9.91013371E-01_kr,&
     -9.83336254E-01_kr,-9.73326828E-01_kr,-9.61008800E-01_kr,&
     -9.46411375E-01_kr,-9.29569172E-01_kr,-9.10522137E-01_kr,&
     -8.89315446E-01_kr,-8.65999398E-01_kr,-8.40629296E-01_kr,&
     -8.13265315E-01_kr,-7.83972359E-01_kr,-7.52819907E-01_kr,&
     -7.19881850E-01_kr,-6.85236313E-01_kr,-6.48965471E-01_kr,&
     -6.11155355E-01_kr,-5.71895646E-01_kr,-5.31279464E-01_kr,&
     -4.89403146E-01_kr,-4.46366017E-01_kr,-4.02270158E-01_kr,&
     -3.57220158E-01_kr,-3.11322872E-01_kr,-2.64687162E-01_kr,&
     -2.17423644E-01_kr,-1.69644420E-01_kr,-1.21462819E-01_kr,&
     -7.29931218E-02_kr,-2.43502927E-02_kr, 2.43502927E-02_kr,&
      7.29931218E-02_kr, 1.21462819E-01_kr, 1.69644420E-01_kr,&
      2.17423644E-01_kr, 2.64687162E-01_kr, 3.11322872E-01_kr,&
      3.57220158E-01_kr, 4.02270158E-01_kr, 4.46366017E-01_kr,&
      4.89403146E-01_kr, 5.31279464E-01_kr, 5.71895646E-01_kr,&
      6.11155355E-01_kr, 6.48965471E-01_kr, 6.85236313E-01_kr,&
      7.19881850E-01_kr, 7.52819907E-01_kr, 7.83972359E-01_kr,&
      8.13265315E-01_kr, 8.40629296E-01_kr, 8.65999398E-01_kr,&
      8.89315446E-01_kr, 9.10522137E-01_kr, 9.29569172E-01_kr,&
      9.46411375E-01_kr, 9.61008800E-01_kr, 9.73326828E-01_kr,&
      9.83336254E-01_kr, 9.91013371E-01_kr, 9.96340117E-01_kr,&
      9.99305042E-01_kr/)
   real(kr),dimension(64),parameter::qw=(/&
      1.78328072E-03_kr, 4.14703326E-03_kr, 6.50445797E-03_kr,&
      8.84675983E-03_kr, 1.11681395E-02_kr, 1.34630479E-02_kr,&
      1.57260305E-02_kr, 1.79517158E-02_kr, 2.01348232E-02_kr,&
      2.22701738E-02_kr, 2.43527026E-02_kr, 2.63774697E-02_kr,&
      2.83396726E-02_kr, 3.02346571E-02_kr, 3.20579284E-02_kr,&
      3.38051618E-02_kr, 3.54722133E-02_kr, 3.70551285E-02_kr,&
      3.85501532E-02_kr, 3.99537411E-02_kr, 4.12625632E-02_kr,&
      4.24735151E-02_kr, 4.35837245E-02_kr, 4.45905582E-02_kr,&
      4.54916279E-02_kr, 4.62847966E-02_kr, 4.69681828E-02_kr,&
      4.75401657E-02_kr, 4.79993886E-02_kr, 4.83447622E-02_kr,&
      4.85754674E-02_kr, 4.86909570E-02_kr, 4.86909570E-02_kr,&
      4.85754674E-02_kr, 4.83447622E-02_kr, 4.79993886E-02_kr,&
      4.75401657E-02_kr, 4.69681828E-02_kr, 4.62847966E-02_kr,&
      4.54916279E-02_kr, 4.45905582E-02_kr, 4.35837245E-02_kr,&
      4.24735151E-02_kr, 4.12625632E-02_kr, 3.99537411E-02_kr,&
      3.85501532E-02_kr, 3.70551285E-02_kr, 3.54722133E-02_kr,&
      3.38051618E-02_kr, 3.20579284E-02_kr, 3.02346571E-02_kr,&
      2.83396726E-02_kr, 2.63774697E-02_kr, 2.43527026E-02_kr,&
      2.22701738E-02_kr, 2.01348232E-02_kr, 1.79517158E-02_kr,&
      1.57260305E-02_kr, 1.34630479E-02_kr, 1.11681395E-02_kr,&
      8.84675983E-03_kr, 6.50445797E-03_kr, 4.14703326E-03_kr,&
      1.78328072E-03_kr/)
   real(kr),parameter::tol=0.02e0_kr
   real(kr),parameter::small=1.e-8_kr
   real(kr),parameter::emax=1.e10_kr
   integer,parameter::imax=10
   real(kr),parameter::zero=0
   real(kr),parameter::one=1

   !--branch to appropriate law
   zp=int(zap/1000)
   zt=int(zat/1000)
   if (irec.gt.0) zp=zt-zp
   ap=awp
   if (irec.gt.0) ap=awrt+1-awp
   at=awrt
   if (law.eq.7) go to 510
   if (law.eq.1.and.lct.eq.1) go to 400
   if (law.eq.1.and.lct.eq.3.and.zap.gt.2004) go to 400
   if (law.eq.1.and.zap.eq.0.) go to 400
   if (law.eq.2.or.law.eq.3) go to 450
   if (law.eq.6) lang=0

   !--one-particle recoil (irec>0) from cm data (law 1 in
   !--the cm system, law 6): exact integration over the cm energy and
   !--cosine of the emitted particle, with the lab recoil energy of each
   !--(e',mu) pair.  The original code (below) integrated over the lab
   !--energy with the mean lab cosine, i.e. D(E_R(<mu>)) instead of
   !--<D(E_R(mu))>, which overestimates the damage (D is concave):
   !--+1 to +10% for the mt91 recoil of Ti-48, V-51, Cu-63 at 5-20 MeV.
   !--ebar (heating) is unchanged within the numerical accuracy.
   ! (DLA, IAEA, Oct 2026)
   if (irec.gt.0.and.(law.eq.1.or.law.eq.6)) then
      ilab6=0
      call rec6cm
      return
   endif

   !--emitted product itself (irec=0) from cm data (law 1 in
   !--the cm system, law 6): the same exact integration over the cm
   !--energy and cosine, with the lab energy of the product for each
   !--(e',w) pair, el = awp/2*(vc**2+u**2+2*vc*u*w), u=sqrt(2e'/awp).
   !--The original adaptive integration over the lab energy (below,
   !--no longer used) could skip the peak of narrow distributions:
   !--errors of more than 10% at 11% of the table energies of the
   !--residuals of Cr and Fe (e.g. ENDF/B-VIII.1 Fe-56 mt16 residual:
   !--normalization 1.02, 1.07 and 0.29 at 19, 20 and 21 MeV).
   ! (DLA, IAEA, Oct 2026)
   if (irec.eq.0.and.(law.eq.1.or.law.eq.6)) then
      ilab6=1
      call rec6cm
      return
   endif

   !--law 1 in the cm system.
   !--adaptive integration over lab energy (original code; not reached
   !--see rec6cm).
   nl=1
   if (irec.gt.0) nl=2
   summ=0
   ebar=0
   dame=0
   x(2)=0
   call h6cm(x(2),epnext,term,nl,lang,lep,c)
   do l=1,nl
      y(2,l)=term(l)
   enddo
   el=x(2)*y(2,1)
   fl=0
  120 continue
   x(1)=epnext
   call h6cm(x(1),epnext,term,nl,lang,lep,c)
   do l=1,nl
      y(1,l)=term(l)
   enddo
   i=2
   dy=(y(1,1)+y(2,1))/10000
  140 continue
   if (i.eq.imax) go to 190
   da=(x(i-1)-x(i))*(y(i-1,1)+y(i,1))/2
   if (da.lt.small) go to 190
   xm=(x(i-1)+x(i))/2
   if (xm.ge.x(i-1).or.xm.le.x(i)) go to 190
   epn=x(i)
   call h6cm(xm,epn,term,nl,lang,lep,c)
   ym=(y(i-1,1)+y(i,1))/2
   test=2*tol*abs(ym)+2*dy
   if (abs(term(1)-ym).gt.test) go to 160
   go to 190
   ! fails test
  160 continue
   i=i+1
   x(i)=x(i-1)
   do l=1,nl
      y(i,l)=y(i-1,l)
   enddo
   x(i-1)=xm
   do l=1,nl
      y(i-1,l)=term(l)
   enddo
   go to 140
   ! passes test
  190 continue
   if (i.eq.1) go to 220
   summ=summ+(x(i-1)-x(i))*(y(i-1,1)+y(i,1))/2
   if (irec.le.0) then
      en=y(i-1,1)*x(i-1)
      ebar=ebar+(x(i-1)-x(i))*(en+el)/2
      f=y(i-1,1)*df(x(i-1),zp,ap,zt,at)
      dame=dame+(x(i-1)-x(i))*(f+fl)/2
   else
      f1=0
      if (y(i-1,1).ne.zero) f1=y(i-1,2)/y(i-1,1)
      en=(e-2*sqrt(e*awp*x(i-1))*f1+awp*x(i-1))/(awrt+1-awp)
      f=y(i-1,1)*df(en,zp,ap,zt,at)
      en=en*y(i-1,1)
      ebar=ebar+(x(i-1)-x(i))*(en+el)/2
      dame=dame+(x(i-1)-x(i))*(f+fl)/2
   endif
   el=en
   fl=f
   i=i-1
   if (i.gt.1) go to 140
   test=emax-emax/100
   if (epnext.ge.test) go to 190
   x(2)=x(1)
   do l=1,nl
      y(2,l)=y(1,l)
   enddo
   go to 120
  220 continue
   return

   !--law 1 for lab system.
   !--trapazoidal energy integration, and
   !--gauss-legendre angle integration, if needed.
  400 continue
   ebar=0
   dame=0
   nd=nint(c(3))
   na=nint(c(4))
   nep=nint(c(6))
   ncyc=na+2
   h=0
   d=0
   xl=0
   yl=0
   el=0
   fl=0
   do i=1,nep
      l=7+ncyc*(i-1)
      xx=abs(c(l))
      yy=c(l+1)
      x2=xx
      if (irec.eq.0) then
         en=yy*xx
         zpp=int(zap/1000)
         ztt=int(zat/1000)
         f2=df(xx,zpp,awp,ztt,awrt)
         f=yy*f2
      else
         en=0
         f=0
         do iq=1,nq
            b=0
            call legndr(qp(iq),p,na)
            do ia=1,na+1
               b=b+(2*ia-1)*c(l+ia)*p(ia)/2
            enddo
            if (yy.ne.zero) b=b/yy
            er=(e-2*sqrt(e*awp*xx)*qp(iq)+awp*xx)/(awrt+1-awp)
            en=en+qw(iq)*b*er
            f=f+qw(iq)*b*df(er,zp,ap,zt,at)
         enddo
         x2=en
         en=en*yy
         f2=f
         f=f*yy
      endif
      if (i.le.nd) then
         h=h+en
         d=d+f
      else if (i.eq.nd+1) then
         xl=xx
         yl=yy
         el=en
         fl=f
      else
         if (lep.eq.1) t1=x2*yl+el
         if (lep.gt.1) t1=en+el
         h=h+(xx-xl)*t1/2
         if (lep.eq.1) t2=f2*yl+fl
         if (lep.gt.1) t2=f+fl
         d=d+(xx-xl)*t2/2
         xl=xx
         yl=yy
         el=en
         fl=f
      endif
   enddo
   ebar=h
   dame=d
   return

   !--laws 2 and 3.
   !--discrete-level scattering
  450 continue
   ebar=0
   dame=0
   lang=0
   if (law.eq.2) lang=nint(c(3))
   nld=0
   if (law.eq.2) nld=nint(c(6))
   thresh=((awrt+1)/awrt)*(-q)
   if (e.lt.thresh) return
   beta=sqrt(awrt*(awrt+1-awp)*(1-thresh/e)/awp)
   afact=awp/(1+awrt)**2
   arec=awp/(awrt+1-awp)
   if (lang.eq.0) then
      ! legendre coefficients
      ! uses gauss-legendre integration for damage
      f1=0
      if (law.ne.3) f1=c(7)
      if (irec.gt.0) then
         f1=-f1
         beta=arec*beta
         afact=afact/arec
      endif
      ebar=afact*e*(1+2*beta*f1+beta*beta)
      if (idame.eq.0) return
      do iq=1,nq
         ! f is the angular distribution of the emitted
         ! particle, evaluated at the particle cm cosine u; the recoil
         ! energy e2 uses the recoil cosine -u, i.e. f_R(-u)=f_p(u).
         ! The original code changed the sign of u before computing f,
         ! so the recoil was given the particle distribution itself
         ! (forward-peaked scattering gave the damage of backward
         ! scattering); ebar (heating) was correct. (DLA, IAEA, Oct 2026)
         u=qp(iq)
         f=one/2
         if (law.ne.3) then
            call legndr(u,p,nld)
            do il=1,nld
               f=f+(2*il+1)*c(6+il)*p(il+1)/2
            enddo
         endif
         if (irec.gt.0) u=-u
         e2=afact*e*(1+2*beta*u+beta*beta)
         dame=dame+qw(iq)*f*df(e2,zp,ap,zt,at)
      enddo
   else if (lang.ne.12) then
      ! tabulated distribution with lang=11 (histogram) or 13-15:
      ! integration with the law of the table (DLA, IAEA, 19 Oct 2026);
      ! the trapezoidal rule below treated every table as lin-lin.
      ! Each interval [x1,x2] of the emitted-particle cosine s is
      ! divided in nsb pieces (1 for a histogram, the value f(x1) on
      ! the whole interval; 16 for lang=13-15, f from terp1); the
      ! recoil (irec>0) uses the cosine u=-s.
      nmu=nint(c(6))
      if (irec.gt.0) then
         beta=arec*beta
         afact=afact/arec
      endif
      lt6=lang-10
      do imu=1,nmu-1
         xa=c(5+2*imu)
         ya=c(6+2*imu)
         xb=c(7+2*imu)
         yb=c(8+2*imu)
         if (xb.le.xa) cycle
         it6=lt6
         if ((it6.eq.3.or.it6.eq.5).and.(xa.le.zero.or.xb.le.zero)) it6=2
         if ((it6.eq.4.or.it6.eq.5).and.(ya.le.zero.or.yb.le.zero)) it6=2
         if (it6.lt.1.or.it6.gt.5) it6=2
         nsb=16
         if (it6.le.2) nsb=1
         do isb=1,nsb
            sa=xa+(xb-xa)*(isb-1)/nsb
            sb=xa+(xb-xa)*isb/nsb
            if (it6.eq.1) then
               fa=ya
               fb=ya
            else
               fa=ya
               if (isb.gt.1) call terp1(xa,ya,xb,yb,sa,fa,it6)
               fb=yb
               if (isb.lt.nsb) call terp1(xa,ya,xb,yb,sb,fb,it6)
            endif
            u=sa
            if (irec.gt.0) u=-sa
            e2=afact*e*(1+2*beta*u+beta*beta)
            el=fa*e2
            fl=fa*df(e2,zp,ap,zt,at)
            u=sb
            if (irec.gt.0) u=-sb
            e2=afact*e*(1+2*beta*u+beta*beta)
            en=fb*e2
            fn=fb*df(e2,zp,ap,zt,at)
            ebar=ebar+(sb-sa)*(en+el)/2
            dame=dame+(sb-sa)*(fn+fl)/2
         enddo
      enddo
   else
      ! tabulated distribution
      ! uses trapazoidal integration for heat and damage
      nmu=nint(c(6))
      if (irec.gt.0) then
         beta=arec*beta
         afact=afact/arec
      endif
      el=0
      fl=0
      ul=0
      do imu=1,nmu
         l=6+2*imu-1
         if (irec.gt.0) l=6+2*nmu-2*imu+1
         u=c(l)
         if (irec.gt.0) u=-u
         yy=c(l+1)
         e2=afact*e*(1+2*beta*u+beta*beta)
         en=yy*e2
         fn=yy*df(e2,zp,ap,zt,at)
         if (imu.gt.1) then
            ebar=ebar+(u-ul)*(en+el)/2
            dame=dame+(u-ul)*(fn+fl)/2
         endif
         ul=u
         el=en
         fl=fn
      enddo
   endif
   return

   !--law 7.
   !--trapazoidal integration over energy for each angle
   !--then trapazoidal integration over angle
  510 continue
   ebar=0
   dame=0
   dl=0
   hl=0
   ul=0
   nmu=nint(c(6))
   l=7+2*nint(c(5))
   do imu=1,nmu
      u=c(l+1)
      nr=nint(c(l+4))
      np=nint(c(l+5))
      iint=nint(c(l+7))
      next=l+6+2*nr+2*np
      ibase=l+5+2*nr
      h=0
      d=0
      do i=1,np
         xx=c(ibase+2*i-1)
         yy=c(ibase+2*i)
         ! NJOY2016.80 overwrote xx with the recoil energy,
         ! so the energy integration used recoil-energy intervals and
         ! the heating used the particle energy (DLA, IAEA, Oct 2026)
         if (irec.gt.0) then ! recoil energy
            xxe=(e-2*sqrt(e*awp*xx)*u+awp*xx)/(awrt+1-awp)
         else
            xxe=xx
         endif
         en=yy*xxe
         if (i.gt.1) then
            if (iint.eq.1) then
               h=h+(xx-xl)*el
            else
               h=h+(xx-xl)*(en+el)/2
            endif
         endif
         f=yy*df(xxe,zp,ap,zt,at)
         if (i.gt.1) then
            if (iint.eq.1) then
               d=d+(xx-xl)*fl
            else
               d=d+(xx-xl)*(f+fl)/2
            endif
         endif
         xl=xx
         el=en
         fl=f
      enddo
      l=next
      if (imu.gt.1) ebar=ebar+(u-ul)*(h+hl)/2
      if (imu.gt.1) dame=dame+(u-ul)*(d+dl)/2
      ul=u
      hl=h
      dl=d
   enddo
   return

   contains

      subroutine rec6cm
      !----------------------------------------------------------------
      ! Cm data: law 1 (lct=2, or lct=3 for a light
      ! particle) and law 6.  For each cm energy e' and cm cosine w
      ! of the emitted particle (mass awp):
      ! ilab6=0, recoil of the residual for the one-particle
      ! approximation (irec>0), mass amr=awrt+1-awp, lab energy
      !   er = amr/2*(vc**2+ur**2-2*vc*ur*w), ur=awp*sqrt(2e'/awp)/amr;
      ! ilab6=1, the emitted product itself (irec=0), lab energy
      !   er = awp/2*(vc**2+u**2+2*vc*u*w), u=sqrt(2e'/awp);
      ! vc=sqrt(2e)/(awrt+1) is the cm velocity.
      ! ebar=int int f(e',w) er dw de', dame=int int f(e',w) D(er) dw de'
      ! (f as returned by h6ddx, h6dis or h6psp, per unit e' and w).
      ! The e' intervals of the table are integrated with 8-point
      ! Gauss-Legendre, the cosine with the 64-point set of getsix.
      ! (DLA, IAEA, Oct 2026)
      !----------------------------------------------------------------
      integer::ndl,npl,ncl,ii,ig,isub,lx
      integer,parameter::nsub6=40
      real(kr)::amr,vc,x1,x2,epc,wep,dum1,dum2,w0,emx
      real(kr),dimension(8),parameter::g8x=(/&
        -0.9602898564975363e0_kr,-0.7966664774136267e0_kr,&
        -0.5255324099163290e0_kr,-0.1834346424956498e0_kr,&
         0.1834346424956498e0_kr, 0.5255324099163290e0_kr,&
         0.7966664774136267e0_kr, 0.9602898564975363e0_kr/)
      real(kr),dimension(8),parameter::g8w=(/&
         0.1012285362903763e0_kr, 0.2223810344533745e0_kr,&
         0.3137066458778873e0_kr, 0.3626837833783620e0_kr,&
         0.3626837833783620e0_kr, 0.3137066458778873e0_kr,&
         0.2223810344533745e0_kr, 0.1012285362903763e0_kr/)

      ebar=0
      dame=0
      amr=awrt+1-awp
      vc=sqrt(2*e)/(awrt+1)
      w0=0
      if (law.eq.6) then
         ! n-body phase space, isotropic in the cm system
         dum1=h6psp(zero,dum2,emx,w0,e,c)
         if (emx.le.zero) return
         do isub=1,nsub6
            x1=(isub-1)*emx/nsub6
            x2=isub*emx/nsub6
            do ig=1,8
               epc=x1+(x2-x1)*(1+g8x(ig))/2
               wep=(x2-x1)*g8w(ig)/2
               call mu6(epc,wep,3,0)
            enddo
         enddo
         return
      endif
      ! law 1: discrete lines (cm energies), then the continuum
      ndl=nint(c(3))
      npl=nint(c(6))
      if (npl.le.0) return
      ncl=nint(c(5))/npl
      if (ndl.gt.0) then
         dum1=h6dis(0,dum2,emx,w0,c,lang)
         do ii=1,ndl
            epc=abs(c(7+ncl*(ii-1)))
            call mu6(epc,one,2,ii)
         enddo
      endif
      if (npl.gt.ndl+1) then
         dum1=h6ddx(zero,dum2,emx,w0,c,lang,lep)
         do ii=ndl+1,npl-1
            lx=7+ncl*(ii-1)
            x1=c(lx)
            x2=c(lx+ncl)
            if (x2.le.x1) cycle
            do ig=1,8
               epc=x1+(x2-x1)*(1+g8x(ig))/2
               wep=(x2-x1)*g8w(ig)/2
               call mu6(epc,wep,1,0)
            enddo
         enddo
      endif
      return
      end subroutine rec6cm

      subroutine mu6(epc,wep,kind,id)
      !----------------------------------------------------------------
      ! Add wep times the integral over the cm cosine of
      ! f(epc,w)*er and f(epc,w)*D(er) to ebar and dame (see rec6cm).
      ! kind=1 continuum (h6ddx), 2 discrete line id (h6dis),
      ! 3 phase space (h6psp). (DLA, IAEA, Oct 2026)
      !----------------------------------------------------------------
      integer::kind,id
      real(kr)::epc,wep
      integer::jq
      real(kr)::amr,vc,ur,er,t,ww,dum1,dum2
      amr=awrt+1-awp
      vc=sqrt(2*e)/(awrt+1)
      ur=awp*sqrt(2*epc/awp)/amr
      if (ilab6.eq.1) ur=sqrt(2*epc/awp)
      do jq=1,nq
         ww=qp(jq)
         if (kind.eq.1) then
            t=h6ddx(epc,dum1,dum2,ww,c,lang,lep)
         else if (kind.eq.2) then
            t=h6dis(id,dum1,dum2,ww,c,lang)
         else
            t=h6psp(epc,dum1,dum2,ww,e,c)
         endif
         if (t.le.zero) cycle
         if (ilab6.eq.1) then
            er=awp*(vc*vc+ur*ur+2*vc*ur*ww)/2
         else
            er=amr*(vc*vc+ur*ur-2*vc*ur*ww)/2
         endif
         ebar=ebar+wep*qw(jq)*t*er
         dame=dame+wep*qw(jq)*t*df(er,zp,ap,zt,at)
      enddo
      return
      end subroutine mu6

   end subroutine getsix

   subroutine h6cm(ep,epnext,term,nl,lang,lep,cnow)
   !-------------------------------------------------------------------
   ! Compute the Legendre coefficients for the double-differential
   ! cross section for secondary energy ep (in the lab system)
   ! using the ENDF-6 File 6 data in the cm system located
   ! in the array cnow.  Call with ep=0 to initialize.
   ! Thereafter, ep can be requested in any order.
   !-------------------------------------------------------------------
   use mathm ! provides legndr
   ! externals
   integer::nl,lang,lep
   real(kr)::ep,epnext,term(*),cnow(*)
   ! internals
   integer::i,na,l,j,ll,ndnow,npnow,ncnow
   real(kr)::aprime,epnn,epm,epn,w,t,epx,xx,cc,c,umin,dm
   real(kr)::umax,un,u,yy,test,epp,s,dx,ym,epnxt
   real(kr)::wlo,whi,f,eps,xc,elmax,e,epmax
   real(kr)::da
   integer,parameter::imax=10
   integer,parameter::mxlg=65
   real(kr)::p(mxlg)
   real(kr)::x(imax),y(imax,mxlg),yt(mxlg)
   real(kr),parameter::tol=0.02e0_kr
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::step=0.05e0_kr
   real(kr),parameter::check=.99999e0_kr
   real(kr),parameter::rdn=.999995e0_kr
   real(kr),parameter::rup=1.000005e0_kr
   real(kr),parameter::tiny=1.e-9_kr
   real(kr),parameter::dn=.999999e0_kr
   real(kr),parameter::up=1.000001e0_kr
   real(kr),parameter::thou=1.e-3_kr
   real(kr),parameter::amil=1.e-6_kr
   real(kr),parameter::one=1.e0_kr
   real(kr),parameter::onep5=1.5e0_kr
   real(kr),parameter::du=.0001e0_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::zero=0
   save na,eps,xc,ndnow,npnow,ncnow,elmax,e,epmax

   !--initialize.
   if (ep.gt.zero) go to 200
   e=cnow(2)
   ndnow=nint(cnow(3))
   npnow=nint(cnow(6))
   ncnow=nint(cnow(5))/npnow
   aprime=awp
   xc=aprime/(awrt+1)**2
   epnext=emax
   elmax=0
   epnn=emax
   epmax=0
   epm=0
   epn=0
   w=0
   if (npnow.ne.ndnow) then
      if (lang.gt.0) then
         t=h6ddx(ep,epn,epm,w,cnow,lang,lep)
      else
         t=h6psp(ep,epn,epm,w,e,cnow)
      endif
      if (epn.lt.epnn*(1-small)) epnn=epn
      if (epm.gt.epmax*(1+small)) epmax=epm
      elmax=e*(sqrt(epmax/e)+sqrt(xc))**2
      elmax=up*elmax
      epn=step*elmax
      if(epmax*up.lt.xc*e) epn=e*(sqrt(epmax/e)-sqrt(xc))**2
      if (epn.lt.epnext) epnext=epn
      eps=tiny/elmax
      epnext=dn*epnext
   endif
   if (ndnow.ne.0) then
      i=0
      t=h6dis(i,epn,epm,w,cnow,lang)
      do i=1,ndnow
         t=h6dis(i,epn,epm,w,cnow,lang)
         if (epn.lt.epnn*(1-small)) epnn=epn
         if (epn.gt.epmax*(1+small)) epmax=epn
         epx=dn*e*(sqrt(epn/e)-sqrt(xc))**2
         if (epx.lt.epnext*(1-small)) epnext=epx
         epx=up*e*(sqrt(epn/e)+sqrt(xc))**2
         if (epx.gt.elmax*(1+small)) elmax=epx
      enddo
   endif
   na=nl-1
   do l=1,nl
      term(l)=0
   enddo
   return

   !--normal entry.
  200 continue
   do l=1,nl
      term(l)=0
   enddo
   if (ndnow.eq.npnow) go to 330
   xx=ep/e
   cc=xc/xx
   c=sqrt(cc)
   umin=(1+cc-epmax/ep)/(2*c)
   if (umin.lt.check) go to 215
   if (c.le.one) go to 420
   if (c.gt.one) go to 320
  215 continue
   if (umin.lt.-one) umin=-one
   dm=(1-umin)/1000+amil
   if (dm.lt.thou/5) dm=thou/5
   epnn=emax

   !--for the continuum part of the distribution,
   !--integrate over the trajectory ep=constant
   !--using adaptive integration on lab cosine.
   umax=1
   un=umax
   x(2)=un
  220 continue
   if (x(2).le.umin) go to 310
   j=2
   do while (j.eq.2)
      u=un
      yy=1+cc-2*c*u
      test=amil
      if (yy.lt.test) then
         yy=amil
         c=u-thou
         cc=c*c
      endif
      epp=yy*ep
      w=(u-c)/sqrt(yy)
      if (lang.gt.0) then
         s=h6ddx(epp,epn,epm,w,cnow,lang,lep)
      else
         s=h6psp(epp,epn,epm,w,e,cnow)
      endif
      if (u.eq.umax.and.epn.lt.epnn) epnn=epn
      call legndr(u,p,na)
      j=1
      if (u.eq.umax) j=2
      x(j)=u
      do l=1,nl
         y(j,l)=p(l)*s/sqrt(yy)
      enddo
      un=u-(epn-epp)/(2*c*ep)
      if (un.gt.u-du) un=u-du
      if (un.lt.umin+du/10) un=umin
   enddo
   i=2
  240 continue
   if (i.eq.imax) go to 290
   dx=x(i)-x(i-1)
   if (dx.lt.dm) go to 290
   da=dx*(y(i,1)+y(i-1,1))/2
   if (da.lt.eps) go to 290
   u=(x(i-1)+x(i))/2
   yy=1+cc-2*c*u
   epp=yy*ep
   w=(u-c)/sqrt(yy)
   if (lang.gt.0) s=h6ddx(epp,epn,epm,w,cnow,lang,lep)
   if (lang.eq.0) s=h6psp(epp,epn,epm,w,e,cnow)
   call legndr(u,p,na)
   do l=1,nl
      yt(l)=p(l)*s/sqrt(yy)
   enddo
   do 250 l=1,nl
   ym=(y(i-1,l)+y(i,l))/2
   if (abs(yt(l)-ym).gt.l*tol*abs(ym)+eps) go to 260
  250 continue
   go to 290
   ! fails test
  260 continue
   i=i+1
   x(i)=x(i-1)
   do l=1,nl
      y(i,l)=y(i-1,l)
   enddo
   x(i-1)=u
   do l=1,nl
      y(i-1,l)=yt(l)
   enddo
   go to 240
   ! passes test
  290 continue
   do l=1,nl
      term(l)=term(l)+(x(i)-x(i-1))*(y(i,l)+y(i-1,l))/2
   enddo
   i=i-1
   if (i.gt.1) go to 240
   x(2)=x(1)
   do l=1,nl
      y(2,l)=y(1,l)
   enddo
   if (i.eq.1) go to 220
  310 continue
   epnxt=e*(sqrt(epnn/e)-sqrt(xc))**2
   if (epnxt.le.ep*(1+small)) epnxt=e*(sqrt(epnn/e)+sqrt(xc))**2
   if (epnxt.le.ep*(1+small)) epnxt=onep5*ep
   if (epnxt.gt.elmax*(1+small).and.elmax.gt.ep*(1+small))epnxt=elmax
   go to 400
  320 continue
   epnxt=(epnext+elmax)/2
   go to 400

   !--check for contributions from delta functions
  330 continue
   if (ndnow.eq.0) go to 400
   epnxt=emax
   do i=1,ndnow
      ll=7+ncnow*(i-1)
      epp=cnow(ll)
      w=(ep-epp-xc*e)/(2*sqrt(xc*e)*sqrt(epp))
      u=(xc*e+ep-epp)/(2*sqrt(xc*e)*sqrt(ep))
      wlo=-1
      whi=1
      if (w.ge.wlo.and.w.le.whi) then
         call legndr(u,p,na)
         f=h6dis(i,epn,epm,w,cnow,lang)
         f=f/(2*sqrt(xc*e*epp))
         epn=e*(sqrt(epp/e)+sqrt(xc))**2
         if (ep.lt.check*epn) then
            epn=rdn*epn
         else
            epn=rup*epn
         endif
         if (epn.lt.epnxt) epnxt=epn
         do l=1,nl
            term(l)=term(l)+f*p(l)
         enddo
      else if (w.lt.wlo) then
         epn=e*(sqrt(epp/e)-sqrt(xc))**2
         if (ep.lt.check*epn) then
            epn=rdn*epn
         else
            epn=rup*epn
         endif
         if (epn.lt.epnxt) epnxt=epn
      endif
   enddo

   !--select epnext
  400 continue
   epnext=epnxt
   return
  420 continue
   epnext=emax
   return
   end subroutine h6cm

   real(kr) function h6ddx(ep,epnext,epmax,w,cnow,lang,lep)
   !-------------------------------------------------------------------
   ! Retrieve the double differential cross section for
   ! secondary energy ep and cosine w (in the cm system) from a
   ! subsection in ENDF-6 File 6 format stored in cnow.  Call with
   ! ep=0 to initialize for each new e. Thereafter, ep can be
   ! requested in any order.  Only the continuum part is
   ! returned.  Delta functions, if given, must be handled
   ! separately.
   !-------------------------------------------------------------------
   use mathm ! provides legndr
   use util  ! provides error
   use endf  ! provides terp1
   ! externals
   integer::lang,lep
   real(kr)::ep,epnext,epmax,w,cnow(*)
   ! internals
   integer::ndnow,npnow,ncnow,mnow,inow,lnow,nl,na,l,lll
   integer::iza2,inn,ii,jj,ia,illdef
   real(kr)::enow,efirst,eplast,t,tt,s,r,aa,tii,tjj
   real(kr)::x1,x2,y1,y2
   integer,parameter::nlmax=65
   real(kr)::p(nlmax)
   real(kr),parameter::up=1.00001e0_kr
   real(kr),parameter::dn=.99999e0_kr
   real(kr),parameter::off=.999995e0_kr
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::zero=0
   save efirst,enow,nl,inow,mnow,lnow,ncnow,na,illdef

   !--initialize
   if (ep.eq.zero) then
      enow=cnow(2)
      ndnow=nint(cnow(3))
      npnow=nint(cnow(6))
      ncnow=nint(cnow(5))/npnow
      mnow=7+(npnow-1)*ncnow
      inow=7+ncnow*ndnow
      lnow=inow
      if (cnow(lnow).le.zero) lnow=lnow+ncnow
      epnext=cnow(lnow)
      epmax=cnow(mnow)
      nl=ncnow-1
      if (nl.gt.nlmax) call error('h6ddx',&
        'too many legendre terms',' ')
      na=nl-1
      efirst=0
      h6ddx=0
      illdef=0
      return
   endif

   !--check for desired interpolation range
   if (ep.lt.efirst*(1-small)) go to 200
   if (ep.gt.cnow(mnow)) go to 190
  130 continue
   epnext=cnow(lnow)
   if (ep.lt.off*epnext) go to 140
   if (lnow.ge.mnow) go to 135
   lnow=lnow+ncnow
   go to 130
  135 continue
   go to 150
  140 continue
   if (lnow.le.ncnow) go to 200
   eplast=cnow(lnow-ncnow)
   if (ep.ge.off*eplast) go to 150
   if (lnow.le.inow+ncnow) go to 150
   lnow=lnow-ncnow
   go to 130
  150 continue
   go to 210

   !--return zero above table
  190 continue
   t=0
   epnext=emax
   go to 350

   !--use analytic form below efirst.
  200 continue
   t=0
   go to 350

   !--legendre coefficients
  210 continue
   if (lang.eq.1) then
      call legndr(w,p,na)
      t=0
      do l=1,nl
         tt=0
         if (l.le.ncnow-1&
             .and.ep.ge.cnow(inow)*(1-small)&
             .and.ep.le.cnow(mnow)*(1+small)) then
            lll=lnow+l
            x1=cnow(lnow-ncnow)
            x2=cnow(lnow)
            if (x1.eq.x2.and.lep.gt.1) then
               x2=sigfig(x2,6,1)
               if (illdef.eq.0) then
                  illdef=1
               endif
            endif
            y1=cnow(lll-ncnow)
            y2=cnow(lll)
            call terp1(x1,y1,x2,y2,ep,tt,lep)
         endif
         t=t+p(l)*(2*l-1)*tt/2
      enddo
      if (t.lt.zero) t=0

   !--kalbach systematics
   else if (lang.eq.2) then
      s=0
      r=0
      if (ep.ge.cnow(inow)*(1-small).and.&
          ep.le.cnow(mnow)*(1+small)) then
         x1=cnow(lnow-ncnow)
         x2=cnow(lnow)
         if (x1.eq.x2.and.lep.gt.1) then
            x2=sigfig(x2,6,1)
            if (illdef.eq.0) then
               illdef=1
            endif
         endif
         y1=cnow(lnow-ncnow+1)
         y2=cnow(lnow+1)
         call terp1(x1,y1,x2,y2,ep,s,lep)
         y1=cnow(lnow-ncnow+2)
         y2=cnow(lnow+2)
         call terp1(x1,y1,x2,y2,ep,r,lep)
      endif
      iza2=nint(zap)
      ! the incident particle is a neutron (izap is the product or the
      ! residual here; NJOY2016.80 used it as the projectile)
      aa=bacha(1,iza2,izat,enow,ep)
      t=aa*(cosh(aa*w)+r*sinh(aa*w))/(2*sinh(aa))
      t=t*s
      if (t.lt.zero) t=0

   !--tabulated angular distribution
   else if (lang.ge.11.and.lang.le.15) then
      inn=lang-10
      ii=lnow-ncnow
      jj=lnow
      x1=cnow(ii)
      x2=cnow(jj)
      if (x1.eq.x2.and.lep.gt.1) then
         x2=sigfig(x2,6,1)
         if (illdef.eq.0) then
            illdef=1
         endif
      endif
      y1=cnow(ii+1)
      y2=cnow(jj+1)
      call terp1(x1,y1,x2,y2,ep,s,lep)
      tii=0
      do ia=1,na,2
         if (w.ge.cnow(ii+1+ia).and.w.le.cnow(ii+3+ia)) then
            call terp1(cnow(ii+1+ia),cnow(ii+2+ia),cnow(ii+3+ia),&
              cnow(ii+4+ia),w,tii,inn)
         endif
      enddo
      tjj=0
      do ia=1,na,2
         if (w.ge.cnow(jj+1+ia).and.w.le.cnow(jj+3+ia)) then
            call terp1(cnow(jj+1+ia),cnow(jj+2+ia),cnow(jj+3+ia),&
              cnow(jj+4+ia),w,tjj,inn)
         endif
      enddo
      call terp1(x1,tii,x2,tjj,ep,t,lep)
      t=t*s

   !--bad value for lang
   else
      call error('h6ddx','illegal lang.',' ')
   endif

   !--return final value
  350 continue
   h6ddx=t
   if (lep.gt.1) return
   if (epnext.eq.emax) return
   if (ep.lt.dn*dn*epnext) then
      epnext=dn*epnext
   else
      epnext=up*epnext
   endif
   return
   end function h6ddx

   real(kr) function h6dis(i,epnext,epmax,w,cnow,lang)
   !-------------------------------------------------------------------
   ! retrieve the double differential cross section for
   ! discrete-energy i and cosine w (in the cm system) from a
   ! subsection in endf-6 file 6 format stored in cnow.  call with
   ! i=0 to initialize for each new e. thereafter, i can be
   ! requested in any order.
   !-------------------------------------------------------------------
   use util  ! provides error
   use mathm ! provides legndr
   use endf  ! provides terp1
   ! externals
   integer::i,lang
   real(kr)::epnext,epmax,w,cnow(*)
   ! internals
   integer::ndnow,npnow,ncnow,inow,lnow,mnow,nl,na,l,iza2,inn,ia
   integer,parameter::mxlg=65
   real(kr)::enow,t,s,r,aa
   real(kr)::p(mxlg)
   real(kr),parameter::zero=0
   save nl,inow,lnow,mnow,ncnow,na,enow

   !--initialize
   if (i.eq.0) then
      enow=cnow(2)
      ndnow=nint(cnow(3))
      npnow=nint(cnow(6))
      ncnow=nint(cnow(5))/npnow
      inow=7
      lnow=inow
      mnow=lnow+(ndnow-1)*ncnow
      epnext=cnow(lnow)
      epmax=cnow(mnow)
      nl=ncnow-1
      na=nl-1
      h6dis=0
      return
   endif

   !--select desired discrete line
   inow=7+(i-1)*ncnow
   epnext=cnow(inow)
   t=0.

   !--legendre coefficients
   if (lang.eq.1) then
      call legndr(w,p,na)
      t=0
      do l=1,nl
         if (l.le.ncnow-1) then
            t=t+(2*l-1)*p(l)*cnow(inow+1+l)/2
         endif
      enddo
      t=t*cnow(inow+1)
      if (t.lt.zero) t=0

   !--kalbach systematics
   else if (lang.eq.2) then
      s=cnow(inow+1)
      r=cnow(inow+2)
      iza2=nint(zap)
      ! the incident particle is a neutron (see h6ddx)
      aa=bacha(1,iza2,izat,enow,epnext)
      t=aa*(cosh(aa*w)+r*sinh(aa*w))/(2*sinh(aa))
      t=t*s
      if (t.lt.zero) t=0

   !--tabulated angular distribution
   else if (lang.ge.11.and.lang.le.15) then
      inn=lang-10
      t=0
      do ia=1,na-2,2
         if (w.ge.cnow(inow+1+ia).and.w.le.cnow(inow+3+ia)) then
            call terp1(cnow(inow+1+ia),cnow(inow+2+ia),&
              cnow(inow+3+ia),cnow(inow+4+ia),w,t,inn)
         endif
      enddo
      t=t*cnow(inow+1)

   !--illegal lang
   else
      call error('h6dis','illegal lang.',' ')
   endif

   !--return final value
   h6dis=t
   return
   end function h6dis

   real(kr) function bacha(iza1i,iza2,izat,e,ep)
   !-------------------------------------------------------------------
   ! Compute the Kalbach a parameter.
   !-------------------------------------------------------------------
   use util ! provides error
   ! externals
   integer::iza1i,iza2,izat
   real(kr)::e,ep
   ! internals
   integer::iza1,iza
   real(kr)::aa,za,ac,zc,ab,zb,nc,nb,na,sa,sb
   real(kr)::ecm,ea,eb,x1,x3,fa,fb
   character(60)::strng
   real(kr),parameter::third=.333333333e0_kr
   real(kr),parameter::twoth=.666666667e0_kr
   real(kr),parameter::fourth=1.33333333e0_kr
   real(kr),parameter::c1=15.68e0_kr
   real(kr),parameter::c2=-28.07e0_kr
   real(kr),parameter::c3=-18.56e0_kr
   real(kr),parameter::c4=33.22e0_kr
   real(kr),parameter::c5=-0.717e0_kr
   real(kr),parameter::c6=1.211e0_kr
   real(kr),parameter::s2=2.22e0_kr
   real(kr),parameter::s3=8.48e0_kr
   real(kr),parameter::s4=7.72e0_kr
   real(kr),parameter::s5=28.3e0_kr
   real(kr),parameter::brk1=130.e0_kr
   real(kr),parameter::brk2=41.e0_kr
   real(kr),parameter::half=0.5e0_kr
   real(kr),parameter::b1=0.04e0_kr
   real(kr),parameter::b2=1.8e-6_kr
   real(kr),parameter::b3=6.7e-7_kr
   real(kr),parameter::tomev=1.e-6_kr
   real(kr),parameter::zero=0

   iza1=iza1i
   if (iza1i.eq.0) iza1=1
   iza=izat
   if (iza.eq.6000) iza=6012
   if (iza.eq.12000) iza=12024
   if (iza.eq.14000) iza=14028
   if (iza.eq.16000) iza=16032
   if (iza.eq.17000) iza=17035
   if (iza.eq.19000) iza=19039
   if (iza.eq.20000) iza=20040
   if (iza.eq.22000) iza=22048
   if (iza.eq.23000) iza=23051
   if (iza.eq.24000) iza=24052
   if (iza.eq.26000) iza=26056
   if (iza.eq.28000) iza=28058
   if (iza.eq.29000) iza=29063
   if (iza.eq.31000) iza=31069
   if (iza.eq.40000) iza=40090
   if (iza.eq.42000) iza=42096
   if (iza.eq.48000) iza=48112
   if (iza.eq.49000) iza=49115
   if (iza.eq.50000) iza=50120
   if (iza.eq.63000) iza=63151
   if (iza.eq.72000) iza=72178
   if (iza.eq.74000) iza=74184
   if (iza.eq.82000) iza=82208
   aa=mod(iza,1000)
   if (aa.eq.zero) then
      write(strng,'(''dominant isotope not known for '',i8)') iza
      call error('bacha',strng,' ')
   endif
   za=int(iza/1000)
   ac=aa+mod(iza1,1000)
   zc=za+int(iza1/1000)
   ab=ac-mod(iza2,1000)
   zb=zc-int(iza2/1000)
   na=nint(aa-za)
   nb=nint(ab-zb)
   nc=nint(ac-zc)
   sa=c1*(ac-aa)&
     +c2*((nc-zc)**2/ac-(na-za)**2/aa)&
     +c3*(ac**twoth-aa**twoth)&
     +c4*((nc-zc)**2/ac**fourth-(na-za)**2/aa**fourth)&
     +c5*(zc**2/ac**third-za**2/aa**third)&
     +c6*(zc**2/ac-za**2/aa)
   if (iza1.eq.1002) sa=sa-s2
   if (iza1.eq.1003) sa=sa-s3
   if (iza1.eq.2003) sa=sa-s4
   if (iza1.eq.2004) sa=sa-s5
   sb=c1*(ac-ab)&
     +c2*((nc-zc)**2/ac-(nb-zb)**2/ab)&
     +c3*(ac**twoth-ab**twoth)&
     +c4*((nc-zc)**2/ac**fourth-(nb-zb)**2/ab**fourth)&
     +c5*(zc**2/ac**third-zb**2/ab**third)&
     +c6*(zc**2/ac-zb**2/ab)
   if (iza2.eq.1002) sb=sb-s2
   if (iza2.eq.1003) sb=sb-s3
   if (iza2.eq.2003) sb=sb-s4
   if (iza2.eq.2004) sb=sb-s5
   ecm=aa*e/ac
   ea=ecm*tomev+sa
   eb=tomev*ep*ac/ab+sb
   x1=eb
   if (ea.gt.brk1) x1=brk1*eb/ea
   x3=eb
   if (ea.gt.brk2) x3=brk2*eb/ea
   fa=1
   if (iza1.eq.2004) fa=0
   fb=1
   if (iza2.eq.1) fb=half
   if (iza2.eq.2004) fb=2
   bacha=b1*x1+b2*x1**3+b3*fa*fb*x3**4
   return
   end function bacha

   real(kr) function h6psp(ep,epnext,epmax,w,e,c)
   !-------------------------------------------------------------------
   ! Compute the double differential cross section for this incident
   ! energy (in the lab system), and secondary energy ep and cosine
   ! w (in the CM system) using the n-body phase-space formula.
   ! Call with ep=0 to initialize the retrieval routine for each
   ! new e.  Thereafter, ep can be requested in any order.
   !-------------------------------------------------------------------
   use endf ! provides contio
   use util ! provides error
   ! externals
   real(kr)::ep,epnext,epmax,w,e,c(*)
   ! internals
   integer::npsx
   real(kr)::apsx,f1,f2,test,s,cn,ex,eimax
   real(kr),parameter::c3=1.2732e0_kr
   real(kr),parameter::c4=3.2813e0_kr
   real(kr),parameter::c5=5.8205e0_kr
   real(kr),parameter::thrhaf=1.5e0_kr
   real(kr),parameter::sevhaf=3.5e0_kr
   real(kr),parameter::step=0.05e0_kr
   real(kr),parameter::rndup=1.000001e0_kr
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::zero=0
   save cn,ex,eimax

   !--initialize
   if (ep.eq.zero) then
      apsx=c(1)
      npsx=nint(c(6))
      f1=(apsx-awp)/apsx
      f2=awrt/(awrt+1)
      ex=thrhaf*npsx-4
      eimax=f1*(f2*e+q)
      if (eimax.le.zero) eimax=1
      if (npsx.eq.3) then
         cn=c3/eimax**2
      else if (npsx.eq.4) then
         cn=c4/eimax**sevhaf
      else if (npsx.eq.5) then
         cn=c5/eimax**5
      else
         call error('h6psp','3, 4, or 5 particles only.',' ')
      endif
      h6psp=0
      epnext=step*eimax
      test=1
      if (eimax.eq.test) epnext=emax
      epmax=eimax
      return
   endif

   !--compute cross section
   s=0
   if (ep.lt.eimax*(1-small)) s=cn*sqrt(ep)*(eimax-ep)**ex
   h6psp=s
   epnext=ep+step*eimax
   if (epnext.gt.eimax*(1+small)) epnext=rndup*eimax
   if (ep.ge.eimax*(1-small)) epnext=emax
   return
   end function h6psp

   subroutine tabsq6(g,h,a,law,z,awr,e,yld)
   !-------------------------------------------------------------------
   ! Compute average of photon recoil energy from capture and
   ! corresponding damage energy for File 6 capture photons
   !-------------------------------------------------------------------
   use endf ! provides terp1
   ! externals
   integer::law
   real(kr)::g,h,a(*),z,awr,e,yld(*)
   ! internals
   integer::nd,np,ncyc,ibase,inn,nc,i,j,ip,ir,idis
   real(kr)::ein,rein,x,y,xr,awc,xh,yh,xl,yl,dx,s,yield,enext
   integer,parameter::nq=4
   real(kr),dimension(nq),parameter::qp=(/&
     -.86114e0_kr,-.33998e0_kr,.33998e0_kr,.86114e0_kr/)
   real(kr),dimension(nq),parameter::qw=(/&
     .34785e0_kr,.65215e0_kr,.65215e0_kr,.34785e0_kr/)

   !--initialize
   ein=2*tm
   rein=1/ein
   g=0
   h=0
   s=0
   nd=nint(a(3))
   np=nint(a(6))
   ncyc=nint(a(5))/np
   ibase=6
   inn=2

   !--interpolate the photon yield
   ip=2
   ir=1
   call terpa(yield,e,enext,idis,yld(1),ip,ir)

   !--accumulate contributions from discrete levels
   if (nd.ne.0) then
      do i=1,nd
         x=a(ibase+ncyc*(i-1)+1)
         y=a(ibase+ncyc*(i-1)+2)
         xr=x*x*rein
         g=g+xr*y
         awc=awr+1
         h=h+df(xr,z,awc,z,awr)*y
         s=s+y
      enddo
   endif

   !--accumulate contributions from continuum
   if (np.gt.nd) then
      nc=np-nd
      ibase=ibase+ncyc*nd
      xh=a(ibase+1)
      yh=a(ibase+2)
      do i=2,nc
         xl=xh
         xh=a(ibase+ncyc*(i-1)+1)
         yl=yh
         yh=a(ibase+ncyc*(i-1)+2)
         if (xl.ne.xh) then
            dx=xh-xl
            do j=1,nq
               x=xl+(1+qp(j))*dx/2
               call terp1(xl,yl,xh,yh,x,y,inn)
               xr=x*x*rein
               g=g+qw(j)*xr*dx*y
               h=h+qw(j)*df(xr,z,awr+1,z,awr)*dx*y
               s=s+qw(j)*dx*y
            enddo
         endif
      enddo
   endif

   !--finished
   g=yield*g/s
   h=yield*h/s
   return
   end subroutine tabsq6

   subroutine hgtfle(e,enext,idis,fle,nle,lcd,matd,mfd,mtd,nin,c,nc)
   !-------------------------------------------------------------------
   ! Retrieve or compute Legendre coefficients at e.
   ! Return isotropic distribution below first point.
   !-------------------------------------------------------------------
   use util ! provides error,mess
   use endf ! provides endf routines
   ! externals
   integer::idis,nle,lcd,matd,mfd,mtd,nin,nc
   real(kr)::e,enext,fle(*),c(nc)
   ! internals
   integer::nb,iso,lvt,ltt,ltt3,lttn,nwc,i,nlmax,il,l
   integer::iraw,ir,nne,ne,inn,innt,lct,nlo,nhi
   real(kr)::test,awr,elo,ehi,fnrm
   integer::lnrm,lnul
   character(60)::strng
   real(kr)::flo(65),fhi(65)
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::zero=0
   save iso,iraw,ir,nne,ne,inn
   save elo,ehi,nlo,nhi,flo,fhi,ltt,ltt3,lttn
   save lct,awr,lnrm,lnul

   !--initialize
   idis=0
   if (e.gt.zero) go to 110
   iso=0
   lnrm=0
   lnul=0
   call findf(matd,mfd,mtd,nin)
   call contio(nin,0,0,c,nb,nwc)
   awr=c(2)
   lvt=nint(c(3))
   ltt=nint(c(4))
   ltt3=ltt
   if (ltt.eq.3) then
      ltt=1
      lttn=1
   endif
   if (lvt.eq.0) then
      call contio(nin,0,0,c,nb,nwc)
   else
      call listio(nin,0,0,c,nb,nwc)
      do while (nb.ne.0)
         call moreio(nin,0,0,c(7),nb,nwc)
      enddo
   endif
   iso=nint(c(3))
   lct=nint(c(4))
   ! check for lab distribution for two-body scattering.
   test=10
   if (lct.eq.1.and.awr.ge.test) then
      if (mtd.ge.51.and.mtd.le.90) then
         write(strng,&
           '(''lab distribution changed to cm for mt '',i3)') mtd
         call mess('hgtfle',strng,' ')
         lct=2
      endif
   endif
   if (iso.ne.1) then
      call tab2io(nin,0,0,c(1),nb,nwc)
      ne=nint(c(6))
      ! read in raw data at first incident energy point
      ! retrieve or compute legendre coefficients.
      iraw=1+nwc
      if (ltt.eq.1) call listio(nin,0,0,c(iraw),nb,nwc)
      if (ltt.eq.2) call tab1io(nin,0,0,c(iraw),nb,nwc)
      if (nb.ne.0) then
         l=iraw
         do while (nb.ne.0)
           l=l+nwc
           call moreio(nin,0,0,c(l),nb,nwc)
         enddo
      endif
      elo=c(iraw+1)
      nlo=nle
      call hgetco(flo,nlo,lcd,c(iraw),lct,ltt,awr,idis,fnrm)
      ! zero table: isotropic moments (fl(1)=1, others 0) (DLA, IAEA, 21 Oct 2026)
      if (ltt.eq.2.and.fnrm.le.zero) then
         call hnull(lnul,mfd,mtd,elo)
      else
         call hnorm(fnrm,lnrm,mfd,mtd,c(iraw+1))
      endif
      if (ltt.eq.1) call listio(nin,0,0,c(iraw),nb,nwc)
      if (ltt.eq.2) call tab1io(nin,0,0,c(iraw),nb,nwc)
      if (nb.ne.0) then
         l=iraw
         do while (nb.ne.0)
           l=l+nwc
           call moreio(nin,0,0,c(l),nb,nwc)
         enddo
      endif
      ehi=c(iraw+1)
      nhi=nle
      call hgetco(fhi,nhi,lcd,c(iraw),lct,ltt,awr,idis,fnrm)
      ! zero table: the low table is kept (DLA, IAEA, 21 Oct 2026)
      if (ltt.eq.2.and.fnrm.le.zero) then
         call hnull(lnul,mfd,mtd,ehi)
         do i=1,nlo
            fhi(i)=flo(i)
         enddo
         nhi=nlo
      else
         call hnorm(fnrm,lnrm,mfd,mtd,c(iraw+1))
      endif
      ! check for total isotropy.
      if (ne.eq.2.and.nlo.eq.1.and.nhi.eq.1) iso=1
   endif
   if (iso.eq.1) then
      nle=1
      enext=etop
   else
      enext=elo
   endif
   nne=2
   ir=1
   inn=nint(c(6+2*ir))
   return

   !--normal entry
   !--is desired energy in current panel
  110 continue
   if (nle.eq.1) go to 150
   if (iso.eq.1) go to 150
   if (e.lt.ehi*(1-small)) go to 130
   if (nne.eq.ne.and.e.lt.ehi+ehi/100) go to 130

   !--no.  slide raw high energy data into low energy positions.
  120 continue
   do i=1,nhi
      flo(i)=fhi(i)
   enddo
   nlo=nhi
   elo=ehi

   !--read in new raw data at high energy.
   if (nne.eq.ne.and.ltt3.eq.3.and.lttn.eq.1) then
      call tab2io(nin,0,0,c(1),nb,nwc)
      ne=nint(c(6))
      nne=0
      ir=1
      ltt=2
      lttn=2
   else if (nne.eq.ne) then
      call error('hgtfle','desired energy above highest given.',' ')
   endif
  121 continue
   if (ltt.eq.1) call listio(nin,0,0,c(iraw),nb,nwc)
   if (ltt.eq.2) call tab1io(nin,0,0,c(iraw),nb,nwc)
   if (nb.ne.0) then
      l=iraw
      do while (nb.ne.0)
        l=l+nwc
        call moreio(nin,0,0,c(l),nb,nwc)
      enddo
   endif
   ehi=c(iraw+1)
   nhi=nle
   call hgetco(fhi,nhi,lcd,c(iraw),lct,ltt,awr,idis,fnrm)
   nne=nne+1
   if (nne.gt.nint(c(5+2*ir))) ir=ir+1
   inn=nint(c(6+2*ir))
   ! table with zero integral (no distribution in the data): skipped,
   ! the neighbouring tables are interpolated; above the last valid
   ! table that one is kept (as ACER, which omits the table from the
   ! ACE grid) (DLA, IAEA, 21 Oct 2026)
   if (ltt.eq.2.and.fnrm.le.zero) then
      call hnull(lnul,mfd,mtd,ehi)
      if (nne.lt.ne) go to 121
      do i=1,nlo
         fhi(i)=flo(i)
      enddo
      nhi=nlo
   else
      call hnorm(fnrm,lnrm,mfd,mtd,c(iraw+1))
   endif
   if (ehi.le.e*(1+small).and.nne.lt.ne) go to 120

   !--yes.  interpolate for coefficients at desired energy.
  130 continue
   if (e.lt.elo*(1-small)) go to 140
   nlmax=nlo
   if (nlmax.lt.nhi) nlmax=nhi
   do i=1,nle
      if (i.le.nlmax) then
         innt=inn
         if (innt.eq.4.or.innt.eq.5) then
            if (flo(i)*fhi(i).le.zero) innt=innt-2
         endif
         call terp1(elo,flo(i),ehi,fhi(i),e,fle(i),innt)
      else
         fle(i)=0
      endif
   enddo
   nle=nlmax
   enext=ehi
   if (inn.eq.1) idis=1
   return

   !--return isotropic distribution below first point.
  140 continue
   fle(1)=1
   do il=2,nle
      fle(il)=0
   enddo
   nle=1
   enext=ehi
   idis=1
   return

   !--isotropic distribution or request.
  150 continue
   fle(1)=1
   if (nle.gt.1) then
      do il=2,nle
         fle(il)=0
      enddo
      nle=1
   endif
   enext=etop
   idis=0
   return
   end subroutine hgtfle

   subroutine hgetco(fl,nl,lcd,c,lct,ltt,awr,idis,fnrm)
   !-------------------------------------------------------------------
   ! Retrieve or compute Legendre coefficients using raw data in
   ! ENDF format.  Conversion between reference systems and
   ! conversion of tabulated data to coefficients is performed by
   ! direct Gaussian quadrature of order 64. Significance of the
   ! results is reduced to reflect precision of calculation.
   !      nl   on entry, maximum number of coefficients desired
   !      nl   on return, number of coefficients required
   !      lcd  system desired (1 for lab, 2 for cm)
   !      lct  system of raw data
   !      ltt  format of raw data (1 for coefficients, 2 for tabulated)
   !      awr  atomic weight ratio
   !      fnrm on return, integral of a tabulated distribution
   ! Tabulated data: exact panel integration with the law of the table
   ! and normalization (see below) (DLA, IAEA, 20 Oct 2026)
   !-------------------------------------------------------------------
   use util  ! provides error
   use mathm ! provides legndr
   use endf  ! provides terpa
   ! externals
   integer::nl,lcd,lct,ltt,idis
   real(kr)::fl(*),c(*),awr,fnrm
   ! internals
   integer::l,np,il,ipc,irc,iq,ip,nlz,j
   real(kr)::x,y,xnext,xn
   integer::nr,ib,i,jr,it
   real(kr)::x1,x2,y1,y2,w
   real(kr)::p(65)
   real(kr),dimension(64),parameter::qp=(/&
     -9.99305042E-01_kr,-9.96340117E-01_kr,-9.91013371E-01_kr,&
     -9.83336254E-01_kr,-9.73326828E-01_kr,-9.61008800E-01_kr,&
     -9.46411375E-01_kr,-9.29569172E-01_kr,-9.10522137E-01_kr,&
     -8.89315446E-01_kr,-8.65999398E-01_kr,-8.40629296E-01_kr,&
     -8.13265315E-01_kr,-7.83972359E-01_kr,-7.52819907E-01_kr,&
     -7.19881850E-01_kr,-6.85236313E-01_kr,-6.48965471E-01_kr,&
     -6.11155355E-01_kr,-5.71895646E-01_kr,-5.31279464E-01_kr,&
     -4.89403146E-01_kr,-4.46366017E-01_kr,-4.02270158E-01_kr,&
     -3.57220158E-01_kr,-3.11322872E-01_kr,-2.64687162E-01_kr,&
     -2.17423644E-01_kr,-1.69644420E-01_kr,-1.21462819E-01_kr,&
     -7.29931218E-02_kr,-2.43502927E-02_kr, 2.43502927E-02_kr,&
      7.29931218E-02_kr, 1.21462819E-01_kr, 1.69644420E-01_kr,&
      2.17423644E-01_kr, 2.64687162E-01_kr, 3.11322872E-01_kr,&
      3.57220158E-01_kr, 4.02270158E-01_kr, 4.46366017E-01_kr,&
      4.89403146E-01_kr, 5.31279464E-01_kr, 5.71895646E-01_kr,&
      6.11155355E-01_kr, 6.48965471E-01_kr, 6.85236313E-01_kr,&
      7.19881850E-01_kr, 7.52819907E-01_kr, 7.83972359E-01_kr,&
      8.13265315E-01_kr, 8.40629296E-01_kr, 8.65999398E-01_kr,&
      8.89315446E-01_kr, 9.10522137E-01_kr, 9.29569172E-01_kr,&
      9.46411375E-01_kr, 9.61008800E-01_kr, 9.73326828E-01_kr,&
      9.83336254E-01_kr, 9.91013371E-01_kr, 9.96340117E-01_kr,&
      9.99305042E-01_kr/)
   real(kr),dimension(64),parameter::qw=(/&
      1.78328072E-03_kr, 4.14703326E-03_kr, 6.50445797E-03_kr,&
      8.84675983E-03_kr, 1.11681395E-02_kr, 1.34630479E-02_kr,&
      1.57260305E-02_kr, 1.79517158E-02_kr, 2.01348232E-02_kr,&
      2.22701738E-02_kr, 2.43527026E-02_kr, 2.63774697E-02_kr,&
      2.83396726E-02_kr, 3.02346571E-02_kr, 3.20579284E-02_kr,&
      3.38051618E-02_kr, 3.54722133E-02_kr, 3.70551285E-02_kr,&
      3.85501532E-02_kr, 3.99537411E-02_kr, 4.12625632E-02_kr,&
      4.24735151E-02_kr, 4.35837245E-02_kr, 4.45905582E-02_kr,&
      4.54916279E-02_kr, 4.62847966E-02_kr, 4.69681828E-02_kr,&
      4.75401657E-02_kr, 4.79993886E-02_kr, 4.83447622E-02_kr,&
      4.85754674E-02_kr, 4.86909570E-02_kr, 4.86909570E-02_kr,&
      4.85754674E-02_kr, 4.83447622E-02_kr, 4.79993886E-02_kr,&
      4.75401657E-02_kr, 4.69681828E-02_kr, 4.62847966E-02_kr,&
      4.54916279E-02_kr, 4.45905582E-02_kr, 4.35837245E-02_kr,&
      4.24735151E-02_kr, 4.12625632E-02_kr, 3.99537411E-02_kr,&
      3.85501532E-02_kr, 3.70551285E-02_kr, 3.54722133E-02_kr,&
      3.38051618E-02_kr, 3.20579284E-02_kr, 3.02346571E-02_kr,&
      2.83396726E-02_kr, 2.63774697E-02_kr, 2.43527026E-02_kr,&
      2.22701738E-02_kr, 2.01348232E-02_kr, 1.79517158E-02_kr,&
      1.57260305E-02_kr, 1.34630479E-02_kr, 1.11681395E-02_kr,&
      8.84675983E-03_kr, 6.50445797E-03_kr, 4.14703326E-03_kr,&
      1.78328072E-03_kr/)
   real(kr),parameter::toler=1.e-6_kr
   real(kr),parameter::rtoler=1.e6_kr
   real(kr),parameter::half=0.5e0_kr
   integer,parameter::nqp=64
   integer,parameter::nlmax=65
   real(kr),parameter::zero=0

   !--determine representation
   l=nl-1
   fnrm=1
   if (nl.gt.nlmax)&
     call error('hgetco','limited to 64 legendre coefficients.',' ')
   if (ltt.eq.2) go to 110
   if (lcd.eq.1.and.lct.ne.1) go to 110
   if (lcd.eq.2.and.lct.lt.2) go to 110

   !--retrieve coefficients from raw data
   np=nint(c(5))+1
   fl(1)=1
   do il=2,nl
      if (il.le.np) then
         fl(il)=c(il+5)
      else
         fl(il)=0
      endif
   enddo
   if (np.lt.nl) nl=np
   return

   !--integrate for fl in desired system
  110 continue
   fl(1)=1
   do il=2,nl
      fl(il)=0
   enddo

   ! tabulated distribution (ltt=2): Legendre moments integrated
   ! panel by panel with the interpolation law of the table (64 Gauss
   ! points in each panel: exact for lin-lin and histogram tables) and
   ! divided by the integral of the table.  The original code sampled
   ! the table at the 64 Gauss points of [-1,1] without normalization:
   ! for forward-peaked tables the moments were wrong (FENDL-3.2d
   ! Gd-156 elastic at 200 MeV: fl(2)=1.00078 instead of 0.99527,
   ! elastic heating 6.8 and damage 3.4 times too low).  Logarithmic
   ! laws with x<=0 or y<=0 are taken as linear.  The integral of the
   ! table is returned in fnrm (message in hnorm).
   ! (DLA, IAEA, 20 Oct 2026)
   if (ltt.eq.2) then
      if (lcd.eq.2.and.lct.eq.1)&
        call error('hgetco','lab to cm conversion not coded.',' ')
      nr=nint(c(5))
      np=nint(c(6))
      ib=6+2*nr
      fnrm=0
      jr=1
      do i=1,np-1
         x1=c(ib+2*i-1)
         y1=c(ib+2*i)
         x2=c(ib+2*i+1)
         y2=c(ib+2*i+2)
         do while (i+1.gt.nint(c(5+2*jr)).and.jr.lt.nr)
            jr=jr+1
         enddo
         it=nint(c(6+2*jr))
         if (x2.le.x1) cycle
         if ((it.eq.3.or.it.eq.5).and.x1.le.zero) it=2
         if ((it.eq.4.or.it.eq.5).and.(y1.le.zero.or.y2.le.zero)) it=2
         if (it.lt.1.or.it.gt.5) it=2
         do iq=1,nqp
            x=(x1+x2)/2+(x2-x1)/2*qp(iq)
            w=(x2-x1)/2*qw(iq)
            call terp1(x1,y1,x2,y2,x,y,it)
            xn=x
            if (lcd.eq.1.and.lct.ge.2)&
              xn=(1+awr*x)/sqrt(1+awr*awr+2*awr*x)
            call legndr(xn,p,l)
            fnrm=fnrm+w*y
            do il=2,nl
               fl(il)=fl(il)+w*y*p(il)
            enddo
         enddo
      enddo
      if (fnrm.gt.zero) then
         do il=2,nl
            fl(il)=fl(il)/fnrm
         enddo
      endif
      go to 120
   endif

   ipc=2
   irc=1
   do iq=1,nqp
      x=qp(iq)

      !--compute scattering probability
      if (ltt.ne.2) then
         ! from coefficients
         np=nint(c(5))
         call legndr(x,p,np)
         y=half
         do ip=1,np
            y=(2*ip+1)*p(ip+1)*c(ip+6)/2+y
         enddo
      else
         ! from tabular data
         call terpa(y,x,xnext,idis,c,ipc,irc)
      endif

      !--multiply scattering probability by legendre polynomials
      !--in desired system and accumulate
      xn=x
      if (lcd.eq.2.and.lct.eq.1)&
        call error('hgetco','lab to cm conversion not coded.',' ')
      if (lcd.eq.1.and.lct.ge.2)&
        xn=(1+awr*x)/sqrt(1+awr*awr+2*awr*x)
      call legndr(xn,p,l)
      y=y*qw(iq)
      do il=2,nl
         fl(il)=fl(il)+y*p(il)
      enddo
   enddo

   !--reduce number of significant figures
   !--and count number of coefficients required
  120 continue
   nlz=1
   do il=2,nl
      j=nint(fl(il)*rtoler)
      fl(il)=j*toler
      if (abs(fl(il)).gt.zero) nlz=il
   enddo
   nl=nlz
   return
   end subroutine hgetco

   subroutine hnull(lnul,mfd,mtd,e)
   !-------------------------------------------------------------------
   ! Warning, once per mt, for a tabulated angular distribution with
   ! zero integral (no distribution in the data): the table is skipped
   ! (hgtfle interpolates between the neighbouring tables and keeps the
   ! last valid one above them; isotropic if no table is valid).
   ! (DLA, IAEA, 21 Oct 2026)
   !-------------------------------------------------------------------
   use util ! provides mess
   ! externals
   integer::lnul,mfd,mtd
   real(kr)::e
   ! internals
   character(60)::strng,strng2

   if (lnul.ne.0) return
   write(strng,'(''*** WARNING mf'',i2,'' mt'',i4,&
     &'': WRONG ANGULAR DISTRIBUTION'')') mfd,mtd
   write(strng2,'(''zero table at'',1p,e10.3,&
     &'' ev SKIPPED (further not listed)'')') e
   call mess('hgtfle',strng,strng2)
   lnul=1
   return
   end subroutine hnull

   subroutine hnorm(fnrm,lnrm,mfd,mtd,e)
   !-------------------------------------------------------------------
   ! Message, once per mt, when the integral of a tabulated angular
   ! distribution differs from 1 by more than 1e-3 (the moments are
   ! normalized in hgetco).
   ! (DLA, IAEA, 20 Oct 2026)
   !-------------------------------------------------------------------
   use util ! provides mess
   ! externals
   integer::lnrm,mfd,mtd
   real(kr)::fnrm,e
   ! internals
   character(60)::strng
   real(kr),parameter::tol=1.e-3_kr

   if (lnrm.ne.0) return
   if (abs(fnrm-1).le.tol) return
   write(strng,'(''mf'',i2,'' mt'',i4,'' table at'',1p,e10.3,&
     &'' ev: integral'',0p,f9.5)') mfd,mtd,e,fnrm
   call mess('hgtfle',strng,'normalized (further tables not listed)')
   lnrm=1
   return
   end subroutine hnorm

   subroutine hconvr(nin,nout,nscr)
   !-------------------------------------------------------------------
   ! Convert photon transition probability arrays (LO=2), if any,
   ! to photon yields (LO=1).  Add MT456 if necessary.
   ! copy all other sections in this material.
   !-------------------------------------------------------------------
   use util ! provides error,mess
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nout,nscr
   ! internals
   integer::l2flg,imax,lmax,i,istor,nb,nw,j,lg,mt0,mt0old,n,jm1
   integer::kk,k,ii,kp1,l,ja,jb,lm1,ip1,mtl,no455,lnu
   integer::mf,mt,l1,l2,n1,n2,nnu,idone,imax2
   integer::i10,maths,mtnow,mtmess,mttst,nn
   integer::m1,m2,nmf3,nmf12
   integer,parameter::nqmx=450
   real(kr)::g,ei,p,ysum,eja,ejb,yy,tsave
   character(60)::strng
   integer::ngam(300)
   integer::mtq(nqmx)
   real(kr)::eeq(nqmx)
   real(kr),dimension(:),allocatable::scr
   real(kr),dimension(:),allocatable::e,eg,es,y,aa,r
   real(kr),dimension(:),allocatable::nu,tmpnu
   real(kr),parameter::zero=0

   !--initialize
   mt0=0
   mt0old=0
   l2flg=0
!   imax=49
   imax=50
   lmax=500
   do i=1,imax
      ngam(i)=0
   enddo
   istor=0
   allocate(scr(npage+50))
   nsh=0
   nsc=0

   !--check for consistancy between mf3 mt's and, if present, mf12
   !  mt's.  if any discrete inelastic levels are defined, we expect
   !  to a sequential list of mt values from 51 to a maximum of 90.
   !  Furthermore, for version 6 formatted endf files there may be
   !  sequential mt values from 600 to a maximum of 648, 650 to a
   !  maximum of 698, 700 to a maximum of 748, 750 to a maximum of
   !  798 and 800 to a maximum of 848 and 876 to 891 for outgoing
   !  protons, deuterons, tritons, 3He, alpha particles and the (n,2n)
   !  reaction, respectively.  The allowed mt interval for version 5
   !  formatted files differs, and is accounted for in the coding that
   !  follows.  The absence of an expected mf3 mt section is flagged
   !  as an error condition; the absence of an expected mf12 mt section
   !  is noted in a message since, while unusual, may be normal.
   mtq=0
   eeq=0
   mtnow=0
   mtmess=0
   i10=0
   nmf3=0
   nmf12=0
  10 continue
   call contio(nin,0,0,scr,nb,nw)
   if (i10.eq.0) then
      i10=1
      maths=math
   endif
   ! an FEND record (mf=0) is skipped here; it was taken as
   ! mf<3 and tofend skipped the whole next file, and the TEND (math=-1)
   ! was not caught, so a material ending after its MF3/MF12 files was
   ! read past the end of the tape (DLA, IAEA, 16 Oct 2026)
   if (math.le.0) go to 20
   if (mfh.eq.0) go to 10
   if (mfh.lt.3 .or. (mfh.gt.3.and.mfh.lt.12)) then
       call tofend(nin,0,0,scr)
       go to 10
   endif
   if (mfh.gt.12) go to 20
   if (mfh.eq.3) then
      if (mth.eq.51.or.                                                &
          (iverf.ge.6.and.(mth.eq.601.or.mth.eq.651.or.mth.eq.701.or.  &
                           mth.eq.751.or.mth.eq.801.or.mth.eq.876)).or.&
          (iverf.le.5.and.(mth.eq.701.or.mth.eq.721.or.mth.eq.741.or.  &
                           mth.eq.761.or.mth.eq.781))) then
         mtnow=mth
      elseif ((mth.gt.51.and.mth.lt.91).or.                 &
              (iverf.ge.6.and.mth.gt.601.and.mth.lt.648).or.&
              (iverf.ge.6.and.mth.gt.651.and.mth.lt.698).or.&
              (iverf.ge.6.and.mth.gt.701.and.mth.lt.748).or.&
              (iverf.ge.6.and.mth.gt.751.and.mth.lt.798).or.&
              (iverf.ge.6.and.mth.gt.801.and.mth.lt.848).or.&
              (iverf.ge.6.and.mth.gt.876.and.mth.lt.890).or.&
              (iverf.le.5.and.mth.gt.701.and.mth.lt.718).or.&
              (iverf.le.5.and.mth.gt.721.and.mth.lt.738).or.&
              (iverf.le.5.and.mth.gt.741.and.mth.lt.758).or.&
              (iverf.le.5.and.mth.gt.761.and.mth.lt.778).or.&
              (iverf.le.5.and.mth.gt.781.and.mth.lt.798)) then
         if (mtnow+1.ne.mth) then
            mtmess=1
            write(strng,'(''mf3, mt'',i2,'' is missing'')')mth-1
            call mess('hconvr',strng,'')
         endif
         mtnow=mth
      endif
      call contio(nin,0,0,scr,nb,nw)
      nmf3=nmf3+1
      mtq(nmf3)=mth
      eeq(nmf3)=scr(2)
      call tosend(nin,0,0,scr)
      go to 10
   elseif (mfh.eq.12) then
      if (mth.eq.51.or.                                                &
          (iverf.ge.6.and.(mth.eq.601.or.mth.eq.651.or.mth.eq.701.or.  &
                           mth.eq.751.or.mth.eq.801.or.mth.eq.876)).or.&
          (iverf.le.5.and.(mth.eq.701.or.mth.eq.721.or.mth.eq.741.or.  &
                           mth.eq.761.or.mth.eq.781))) then
         mtnow=mth
      elseif ((mth.gt.51.and.mth.lt.91).or.                 &
              (iverf.ge.6.and.mth.gt.601.and.mth.lt.648).or.&
              (iverf.ge.6.and.mth.gt.651.and.mth.lt.698).or.&
              (iverf.ge.6.and.mth.gt.701.and.mth.lt.748).or.&
              (iverf.ge.6.and.mth.gt.751.and.mth.lt.798).or.&
              (iverf.ge.6.and.mth.gt.801.and.mth.lt.848).or.&
              (iverf.ge.6.and.mth.gt.876.and.mth.lt.890).or.&
              (iverf.le.5.and.mth.gt.701.and.mth.lt.718).or.&
              (iverf.le.5.and.mth.gt.721.and.mth.lt.738).or.&
              (iverf.le.5.and.mth.gt.741.and.mth.lt.758).or.&
              (iverf.le.5.and.mth.gt.761.and.mth.lt.778).or.&
              (iverf.le.5.and.mth.gt.781.and.mth.lt.798)) then
         if (mtnow+1.ne.mth) then
            write(strng,'(''mf12, mt'',i2,'' may be missing'')')mth-1
            call mess('hconvr',strng,&
                               'discrete photon data may be incomplete')
         endif
         mtnow=mth
      endif
      nmf12=nmf12+1
      call tosend(nin,0,0,scr)
      go to 10
   endif
   if (mtmess.ne.0) then
      call error('hconvr',&
                 'missing mf3 mt''s, probable endf error','')
   endif

   !--mt checks are done, rewind the endf tape and position
   !  at the first record for this material
  20 continue
   call repoz(nin)
   call findf(maths,1,451,nin)

   !--loop over all sections of this material
  110 continue
   call contio(nin,0,0,scr,nb,nw)
   if (mfh.eq.0.and.math.ne.0) go to 130
   if (mfh.eq.1) go to 550
   if (mfh.eq.12.and.mth.ne.460) go to 120
   if (math.eq.1149.and.mfh.eq.13.and.mth.ge.51.and.iverf.lt.5)&
     go to 140
   if (math.eq.1150.and.mfh.eq.13.and.mth.le.54.and.iverf.lt.5)&
     go to 140
   if (mfh.eq.14.and.l2flg.eq.1) go to 500
   call contio(0,nout,nscr,scr,nb,nw)
   if (math.eq.0) go to 600
   call tosend(nin,nout,nscr,scr)
   go to 110
  120 continue
   if (l1h.eq.2) go to 150
   call contio(0,nout,nscr,scr,nb,nw)
   call tosend(nin,nout,nscr,scr)
   go to 110
  130 continue
   call contio(0,nout,nscr,scr,nb,nw)
   go to 110
  140 continue
   call tosend(nin,0,0,scr)
   write(strng,'(''gamma prod patch made for mt '',i3)') mth
   call mess('hconvr',strng,' ')
   go to 110

   !--allocate storage
  150 continue
   if (istor.eq.0) then
      istor=1
      imax2=imax*imax
      allocate(e(imax))
      allocate(eg(lmax))
      allocate(es(lmax))
      allocate(y(lmax))
      allocate(aa(imax2))
      allocate(r(imax2))
      e(1)=0
      do i=1,imax2
         r(i)=0
         aa(i)=0
      enddo
      do i=1,imax
         r(imax*(i-1)+i)=1
      enddo
   endif

   !--convert transition probability array to yields
   za=c1h
   awr=c2h
   lg=l2h
   g=1
   l2flg=1
   l=1
   call listio(nin,0,0,scr(l),nb,nw)
   l=l+nw
   do while (nb.ne.0)
      call moreio(nin,0,0,scr(l),nb,nw)
      l=l+nw
   enddo

   !--set base value for mt0
   if (mth.ge.51.and.mth.le.91.and.mt0.ne.49) mt0=49
   if (iverf.ge.6) then
      if (mth.ge.601.and.mth.le.649.and.mt0.ne.599) mt0=599
      if (mth.ge.651.and.mth.le.699.and.mt0.ne.649) mt0=649
      if (mth.ge.701.and.mth.le.749.and.mt0.ne.699) mt0=699
      if (mth.ge.751.and.mth.le.799.and.mt0.ne.749) mt0=749
      if (mth.ge.801.and.mth.le.849.and.mt0.ne.799) mt0=799
      if (mth.ge.876.and.mth.le.891.and.mt0.ne.874) mt0=874
   elseif (iverf.le.5) then
      if (mth.ge.701.and.mth.le.719.and.mt0.ne.699) mt0=699
      if (mth.ge.721.and.mth.le.739.and.mt0.ne.719) mt0=719
      if (mth.ge.741.and.mth.le.759.and.mt0.ne.739) mt0=739
      if (mth.ge.761.and.mth.le.779.and.mt0.ne.759) mt0=759
      if (mth.ge.781.and.mth.le.799.and.mt0.ne.779) mt0=779
   endif

   !--load the e(_) array with mf3 -q.  Will overwrite with mf12
   !  mt data when possible (which should be always but if there
   !  is a missing mf12 mt, we're covered).
   if (mt0.ne.mt0old) then
      e=0
      mt0old=mt0
      m1=mt0+2
      if (mt0.eq.49) then
         m2=91
      elseif (iverf.ge.6.and.mt0.ne.874) then
         m2=m1+48
      elseif (iverf.ge.6.and.mt0.eq.874) then
         m2=m1+15
      elseif (iverf.le.5) then
         m2=m1+17
      endif
      nn=1
      mttst=-1
      do while (mttst.lt.m2.and.nn.le.nmf3)
         mttst=mtq(nn)
         if (mttst.ge.m1.and.mttst.le.m2) e(mttst-mt0)=-eeq(nn)
         nn=nn+1
      enddo
   endif
   j=mth-mt0
   e(j)=scr(1)
   n=n2h
   jm1=j-1
   do kk=1,jm1
      k=jm1-kk+1
      ii=0
      idone=0
      do while (ii.lt.n.and.idone.eq.0)
         ii=ii+1
         i=ii
         ei=scr(6+(lg+1)*i-lg)
         if (ei.eq.zero.and.e(k).eq.zero) then
            idone=1
         elseif (ei.ne.zero) then
            if (abs(ei-e(k))/ei.le.0.0001_kr) idone=1
         endif
      enddo
      if (idone.eq.0) then
         aa((k-1)*imax+j)=0
         r((k-1)*imax+j)=0
      else
         p=scr(7+(lg+1)*i-lg)
         g=1
         if (lg.eq.2) g=scr(6+(lg+1)*i)
         aa((k-1)*imax+j)=p*g
         r((k-1)*imax+j)=p
      endif
      if (k.ne.jm1) then
         kp1=k+1
         do ii=kp1,jm1
            r((k-1)*imax+j)=r((k-1)*imax+j)&
              +r((ii-1)*imax+j)*aa((k-1)*imax+ii)
         enddo
      endif
   enddo
   l=0
   ysum=0
   do i=2,j
      ja=j+2-i
      eja=e(ja)
      do ii=1,jm1
         jb=j-ii
         ejb=e(jb)
         yy=aa((jb-1)*imax+ja)*r((ja-1)*imax+j)
         if (yy.gt.zero) then
            l=l+1
            if (l.gt.lmax) then
               call error('hconvr','too many lo=2 gammas',' ')
            endif
            eg(l)=eja-ejb
            es(l)=eja
            y(l)=yy
            ysum=ysum+y(l)
         endif
      enddo
   enddo
   e(1)=0

   !--arrange gamma energies in descending order
   if (l.gt.1) then
      lm1=l-1
      do i=1,lm1
         ip1=i+1
         do ii=ip1,l
            if (eg(i).lt.eg(ii)) then
               tsave=eg(i)
               eg(i)=eg(ii)
               eg(ii)=tsave
               tsave=y(i)
               y(i)=y(ii)
               y(ii)=tsave
               tsave=es(i)
               es(i)=es(ii)
               es(ii)=tsave
            endif
         enddo
      enddo
   endif

   !--output computed yields in endf lo=1 format
   scr(1)=za
   scr(2)=awr
   scr(3)=1
   scr(4)=0
   scr(5)=l
   scr(6)=0
   mtl=49
   if (iverf.le.5.and.mth.ge.700) mtl=649
   if (iverf.ge.6.and.mth.ge.600) mtl=549
   ngam(mth-mtl)=l
   call contio(0,nout,nscr,scr,nb,nw)
   ! output tab1 sum record
   if (l.gt.1) then
      scr(1)=0
      scr(2)=0
      scr(3)=0
      scr(5)=1
      scr(6)=2
      scr(7)=2
      scr(8)=2
      scr(9)=ebot
      scr(10)=ysum
      scr(11)=etop
      scr(12)=ysum
      nw=12
      call tab1io(0,nout,nscr,scr,nb,nw)
   endif
   ! output rest of tab1 records
   do i=1,l
      scr(1)=eg(i)
      scr(2)=es(i)
      scr(3)=0
      scr(4)=2
      scr(5)=1
      scr(6)=2
      scr(7)=2
      scr(8)=2
      scr(9)=ebot
      scr(10)=y(i)
      scr(11)=etop
      scr(12)=y(i)
      nw=12
      call tab1io(0,nout,nscr,scr,nb,nw)
   enddo

   !--copy send record and loop over remaining sections.
   call contio(nin,nout,nscr,scr,nb,nw)
   go to 110

   !--angular distributions
  500 continue
   if (mth.ge.51.and.mth.le.90) then
      scr(3)=1
      scr(4)=0
      scr(6)=0
      scr(5)=ngam(mth-50)
      call contio(0,nout,nscr,scr,nb,nw)
      call tosend(nin,0,0,scr)
      call asend(nout,nscr)
   else
      call contio(0,nout,nscr,scr,nb,nw)
      call tosend(nin,nout,nscr,scr)
   endif
   go to 110

   !--file1.  add 456=452 if necessary.
  550 continue
   call contio(0,nout,nscr,scr,nb,nw)
   call tosend(nin,nout,nscr,scr)
   no455=0
   call contio(nin,nout,nscr,scr,nb,nw)
   if (mfh.eq.0) go to 110
   if (mth.ne.452) go to 595
!   nnu=9999
!   allocate(nu(nnu))
   if(allocated(tmpnu))deallocate(tmpnu)
   allocate(tmpnu(npage+46))
   l=1
   lnu=l2h
!   if (lnu.eq.1) call listio(nin,nout,nscr,nu(l),nb,nw)
!   if (lnu.eq.2) call tab1io(nin,nout,nscr,nu(l),nb,nw)
   if (lnu.eq.1) then
      call listio(nin,nout,nscr,tmpnu(l),nb,nw)
   elseif (lnu.eq.2) then
      call tab1io(nin,nout,nscr,tmpnu(l),nb,nw)
   endif
   nnu=nb+nw+1
   if (allocated(nu))deallocate(nu)
   allocate(nu(nnu))
   nu(1:nw)=tmpnu(1:nw)
   do while (nb.ne.0)
      l=l+nw
      if (l.gt.nnu) call error('hconvr',&
        'exceeded storage for nubar',' ')
      call moreio(nin,nout,nscr,nu(l),nb,nw)
   enddo
   call tosend(nin,nout,nscr,scr)
  565 continue
   call contio(nin,0,0,scr,nb,nw)
   if (mth.eq.456) go to 590
   if (mth.eq.0.and.no455.eq.0) go to 570
   if (mfh.eq.0) go to 590
   if (mth.eq.455) no455=1
   if (mth.gt.456.and.no455.eq.0) go to 570
   call contio(0,nout,nscr,scr,nb,nw)
   call tosend(nin,nout,nscr,scr)
   go to 565
  570 continue
   if (no455.eq.1) go to 590
   mf=mfh
   mt=mth
   l1=l1h
   l2=l2h
   n1=n1h
   n2=n2h
   mfh=1
   mth=456
   scr(3)=0
   scr(4)=lnu
   scr(5)=0
   scr(6)=0
   l=1
   call contio(0,nout,nscr,scr,nb,nw)
   if (lnu.eq.1) call listio(0,nout,nscr,nu(l),nb,nw)
   if (lnu.eq.2) call tab1io(0,nout,nscr,nu(l),nb,nw)
   do while (nb.ne.0)
      l=l+nw
      if (l.gt.nnu) call error('hconvr',&
        'exceeded storage for nubar',' ')
      call moreio(0,nout,nscr,nu(l),nb,nw)
   enddo
   call asend(nout,nscr)
   mfh=mf
   mth=mt
   scr(3)=l1
   scr(4)=l2
   scr(5)=n1
   scr(6)=n2
  590 continue
   call contio(0,nout,nscr,scr,nb,nw)
   deallocate(nu)
   if (mfh.eq.0) go to 110
  595 continue
   call tofend(nin,nout,nscr,scr)
   go to 110

   !--conver is finished.
  600 continue
   if (istor.gt.0) then
      deallocate(e)
      deallocate(eg)
      deallocate(es)
      deallocate(y)
      deallocate(aa)
      deallocate(r)
   endif
   deallocate(scr)
   return
   end subroutine hconvr

   subroutine hgam102(e,ebar,dame,disc102,c,irec,zp,ap,zt,at)
   !-------------------------------------------------------------------
   ! Process the relativistic discrete gamma or its recoil as
   ! given in mf6/mt102 for ENDF/B-VII neutron + H-1.
   !-------------------------------------------------------------------
   use endf ! provides terpa
   ! externals
   integer::irec
   real(kr)::e,ebar,dame,disc102,zp,ap,zt,at
   real(kr)::c(*)
   ! internals
   real(kr)::er,eg2,enext,yld
   integer::idisc,ip,ir

   !--interpolate yield (DLA, IAEA, Oct 2026)
   ip=2
   ir=1
   call terpa(yld,e,enext,idisc,c(1),ip,ir)

   if (irec.eq.0) then
      ebar=disc102+e*awr/(1+awr)
      dame=0
   else
      !--include recoil energy plus photon "kick" energy
      er=e/(awr+1)
      eg2=yld*disc102*disc102*rtm/2
      ebar=er+eg2
      dame=df(ebar,zp,ap,zt,at)
   endif
   return
   end subroutine hgam102

   subroutine gheat(iold,inew,nscr)
   !-------------------------------------------------------------------
   ! Subtract energy carried away by gammas.
   !-------------------------------------------------------------------
   use util ! provides error,mess,loada,finda,sigfig
   use mainio ! provides nsyso
   use endf ! provides endf routines and variables
   ! externals
   integer::iold,inew,nscr
   ! internals
   integer::nb,nw,mfd,mgam,nk,ik,mfc,mtc,npkk,npkt
   integer::lo,mfx,mtx,mtd,jdis,ilist,nwd,lp,lqx,idis
   integer::i,ipx,irx,nmt,mty,idx,ii,indxx,it,isave,jerr,nd
   real(kr)::afact,aw1fac,z,e,dame
   real(kr)::cerr,enxt,el,elow,ehigh,test,thresh
   real(kr)::egkr,ebar,egam,edam,damn,enext,enx,egk
   real(kr)::h,hk,cfix,eava,ebarp,xp,yp,hp,egamp,edamp
   real(kr)::damep,hx,hxp,subtot,elo,ehi,x,y
   real(kr)::en102
   real(kr),dimension(:),allocatable::sk,sd,sy
   real(kr)::c(50)
   integer::imt(50)
   real(kr),dimension(:),allocatable::scr
   real(kr),dimension(:),allocatable::b
   real(kr),dimension(:),allocatable::d
   real(kr),dimension(:),allocatable::bufo,bufn
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::smin=1.e-9_kr
   real(kr),parameter::qtest=99.e6_kr
   real(kr),parameter::zero=0
   integer::mf6flg
   character(len=70)::strng1,strng2

   !--allocate buffers for loada/finda
   allocate(bufo(nbuf))
   allocate(bufn(nbuf))
   ! capture photon kicks summed over all subsections (DLA, IAEA, Oct 2026)
   allocate(sk(ne),sd(ne),sy(ne))
   sk=0
   sd=0
   sy=0

   !--find first section of mf12 or mf13
   allocate(scr(npage+50))
   mfd=3
   afact=awr/(awr+1)
   aw1fac=1/(awr+1)
   z=nint(za/1000)
   e=0
   dame=df(e,z,awr,z,awr)
   mgam=2
   hk=0
  100 continue
   call contio(nscr,0,0,scr,nb,nw)
   if (mfh.eq.12) mgam=1
   if (mfh.ge.12) go to 102
   call tofend(nscr,0,0,scr)
   go to 100
  102 continue
   mfc=mfh
   mtc=mth
   call findf(matd,mfc,mtc,nscr)
   npkk=npk
   if (kchk.eq.1) npkk=3*npk-2
   if (mgam.gt.0) npkk=npkk+3
   npkt=npk-1

   !--sum all reactions found in mf12 and mf13.
   !--reconstruct redundant cross sections if needed.
   allocate(b(npage+50))
   nd=10000
   allocate(d(nd))
   if (mgam.eq.2) go to 210
  105 continue
   call contio(nscr,0,0,scr,nb,nw)
   if (mfh.gt.13) go to 250
   if (mfh.eq.0) go to 220
   if (mth.eq.0) go to 105
   if (mfh.eq.12.and.mth.eq.460) then
      call tosend(nscr,0,0,scr)
      go to 105
   endif
   ! skip over this mf/mt if photon data were already found in mf6
   if (i6p.gt.0) then
      mf6flg=0
      do i=1,i6p
         if (mt6yp(i).ne.mth) cycle
         mf6flg=1
         exit
      enddo
      if (mf6flg.ne.0) then
         write(strng1,'(''skipping mf'',i2,''/mt = '',i3)')mfh,mth
         write(strng2,'(''photons were already processed in mf6'')')
         call mess('gheat ',strng1,strng2)
         call tosend(nscr,0,0,scr)
         go to 105
      endif
   endif
   lo=l1h
   if (lo.eq.2) call error('gheat','lo=2 not coded.',' ')
   mfx=mfh
   mtd=mth
   mtx=mth
   if (mth.eq.3) mtd=1
   if (iprint.eq.1.and.mfh.eq.12.and.idame.gt.0.and.mtd.eq.102)&
     write(nsyso,'(/&
       &'' photon energy (from yields) mf'',i2,'', mt'',i3/&
       &14x,''e'',6x,''ebar/err'',10x,''egam'',10x,''edam'',10x,&
       &''xsec'',9x,''yield'',7x,''heating'',8x,''damage'')')&
       mfh,mth
   if (iprint.eq.1.and.mfh.eq.12.and.(idame.eq.0.or.mtd.ne.102))&
     write(nsyso,'(/&
       &'' photon energy (from yields) mf'',i2,'', mt'',i3/&
       &14x,''e'',10x,''ebar'',10x,''xsec'',9x,''yield'',&
       &8x,''energy'',7x,''heating'')') mfh,mth
   if (iprint.eq.1.and.mfh.eq.13) write(nsyso,'(/&
     &'' photon energy (from xsecs) mf'',i2,'', mt'',i3/&
     &14x,''e'',10x,''ebar'',10x,''xsec'',8x,''energy'',7x,&
     &''heating'')') mfh,mth
  110 continue
   nk=n1h
   ik=0
   cerr=0
   ! check for leading or trailing zeroes to locate ranges
   ! where this partial is not defined
   e=0
   call gety1(e,enxt,jdis,y,nscr,scr)
   e=enxt
   el=sigfig(e,7,0)
   elow=big
   ehigh=0
   test=big-big/100
   do while (e.lt.test)
      call gety1(e,enxt,jdis,y,nscr,scr)
      if (y.gt.zero.and.e.lt.elow) elow=el
      if (y.gt.zero) ehigh=0
      if (ehigh.eq.zero.and.y.eq.zero) ehigh=sigfig(e,7,0)
      el=sigfig(e,7,0)
      e=enxt
   enddo
   if (ehigh.eq.zero) ehigh=sigfig(e,7,0)
   if (nk.eq.1) then
      call findf(matd,mfx,mtx,nscr)
      call contio(nscr,0,0,scr,nb,nw)
   endif
  130 continue
   e=0
   call gety1(e,thresh,jdis,y,nscr,scr)
   ilist=1
   egkr=c1h
   nwd=nd
   if (egkr.eq.zero)&
     call gambar(e,ebar,egam,edam,nendf,matd,mtx,d,nwd)
   lp=l1h
   ik=ik+1
   if (ik.eq.1) then
      sk=0
      sd=0
      sy=0
   endif
   if (iprint.eq.1.and.egkr.eq.zero) write(nsyso,&
     '(i6,''  continuum gammas'')') ik
   if (iprint.eq.1.and.egkr.gt.zero) write(nsyso,&
     '(i6,1p,e13.4,'' ev gamma'')') ik,egkr
   lqx=0
   if (mfd.eq.13) go to 150
   call findf(matd,mfd,mtd,nin)
   call contio(nin,0,0,b,nb,nw)
   call gety2(e,enxt,idis,x,nin,b)
   if (mth.eq.102) then
      q=c2h
      if (nqa.gt.0) then
         do i=1,nqa
            if (mta(i).eq.mth) then
               q=qa(i)
               if (q.ge.qtest) lqx=lqs(i)
            endif
         enddo
      endif
      call capdam(e,damn,q,za,awr,mth)
      call disgam(e,egam,edam,z,awr)
      ipx=2
      irx=1
   endif
   if (enxt.gt.thresh) thresh=enxt

   !--loop over all energies in grid
  150 continue
   if (ik.gt.1) go to 160
   nmt=0
   if (npk.lt.3) go to 160
   mty=mtx
   if (mty.eq.3) mty=1
   if (mty.eq.4) mty=51
   if (mty.eq.18.and.mt19.eq.1) mty=19
   call indx(npk,mtp,mt19,mty,imt,nmt)
  160 continue
   i=0
  190 continue
   i=i+1
   call finda(i,c,npkk,iold,bufo,nbuf)
   e=sigfig(c(1),9,0)
   if (e.lt.elow*(1-small).or.e.gt.ehigh*(1+small)) go to 180
   call gety1(e,enext,jdis,y,nscr,scr)
   x=zero
   if (mfd.eq.3) then
      call gety2(e,enext,idis,x,nin,b)
      if (x.eq.zero) go to 195
   endif
   if (lqx.ne.0) call terpa(q,e,enx,idx,qbar(lqx),ipx,irx)
   egk=egkr
   if (lp.eq.2) egk=egkr+e*afact
   ebar=egk
   dame=0
   if (ik.eq.1) c(npkk)=0
   if (egkr.eq.zero)&
     call gambar(e,ebar,egam,edam,nendf,matd,mtx,d,nwd)
   if (mfd.eq.13) go to 170
   if (mth.ne.102) go to 164
   ! photon recoil correction
   ! primary photons (lp=2) use the actual energy egk (DLA, IAEA, Oct 2026)
   if (egkr.ne.zero) call disgam(egk,egam,edam,z,awr)
   h=egam*x*y
   hk=h
   c(npkk-1)=c(npkk-1)+x*y*ebar
   c(npkk)=c(npkk)+y*ebar
   ! the photon kicks (sk), their damage (sd) and the
   ! photon multiplicity (sy) are summed over all subsections; the
   ! capture damage is set at the last subsection (DLA, IAEA, Oct 2026)
   sk(i)=sk(i)+y*egam
   sd(i)=sd(i)+y*edam
   sy(i)=sy(i)+y
   dame=0
   if (ik.eq.nk) then
      if (idame.gt.0) then
         ! manual eq. 210-211: the capdam default (one photon with all
         ! the available energy) is replaced by the neutron kick and
         ! the photon kicks of the data, treated independently;
         ! one photon: D(neutron kick + photon kick).
         ! (NJOY2016.80 subtracted the whole capdam value, D of the
         ! sum of both kicks, and lost the neutron kick)
         call capdam(e,damn,q,za,awr,mth)
         en102=e*aw1fac
         if (sy(i).le.ysing) then
            dame=x*df(en102+sk(i),z,awr+1,z,awr)-damn*x
         else
            dame=x*(df(en102,z,awr+1,z,awr)+sd(i))-damn*x
         endif
      endif
      hk=hk-x*rtm*(q+e*awr*aw1fac)**2/2
      cfix=(e+q-e*aw1fac)*x
      h=h-cfix
      c(npkk-2)=cfix-x*c(npkk)
      eava=afact*e+q
      cerr=100*(c(npkk)-eava)/eava
   endif
   ebarp=sigfig(ebar,9,0)
   xp=sigfig(x,9,0)
   yp=sigfig(y,9,0)
   hp=sigfig(h,9,0)
   egamp=sigfig(egam,9,0)
   edamp=sigfig(edam,9,0)
   damep=sigfig(dame,9,0)
   cerr=sigfig(cerr,4,0)
   if (iprint.ne.1.or.abs(e-elist(ilist)).gt.small*e&
     .or.y.le.smin) go to 166
   if (idame.eq.0) then
      write(nsyso,'(1x,1p,8e14.4)') e,ebarp,xp,yp,hp
   else
      write(nsyso,'(1x,1p,8e14.4)')&
        e,ebarp,egamp,edamp,xp,yp,hp,damep
   endif
   if (ik.eq.nk) write(nsyso,&
     '(1x,1p,e14.4,0p,f10.1,'' pc'')') e,cerr
   go to 166
   ! other mf12 reactions
  164 continue
   if (mtd.ne.2.and.mtd.ne.460) h=-y*x*ebar
   if (mtd.eq.2) h=y*x*ebar
   c(npkk-1)=c(npkk-1)-h
   c(npkk)=c(npkk)+h
  166 continue
   hx=h
   c(2)=c(2)+h
   if (nmt.gt.0) then
      do ii=1,nmt
         indxx=imt(ii)
         if (mtp(indxx).ge.444) c(indxx)=c(indxx)+dame
         if (mtp(indxx).lt.442) c(indxx)=c(indxx)+h
         if (mtp(indxx).eq.442) c(indxx)=c(indxx)-h
         if (mtp(indxx).eq.443.and.mth.eq.102) c(indxx)=c(indxx)+hk
      enddo
   endif
   if (iprint.eq.0.or.e.ne.elist(ilist)) go to 180
   if (mth.eq.102) go to 180
   if (ik.lt.nk.and.y.lt.smin) go to 180
   if (nk.eq.1.and.y.lt.smin) go to 180
   test=-1
   test=test/10000
   if (nk.gt.1.and.ik.eq.nk.and.c(npkk).gt.test) go to 180
   ebarp=sigfig(ebar,9,0)
   xp=sigfig(x,9,0)
   yp=sigfig(y,9,0)
   hxp=sigfig(hx,9,0)
   hp=sigfig(h,9,0)
   write(nsyso,'(1x,1p,8e14.4)') e,ebarp,xp,yp,-hxp,hp
   if (ik.eq.nk.and.nk.gt.1) then
      subtot=-c(npkk)
      subtot=sigfig(subtot,9,0)
      write(nsyso,'(48x,''subtot   '',1p,e14.4)') subtot
   endif
   go to 180
   ! mf13 reactions
  170 continue
   h=-y*ebar
   c(npkk-1)=c(npkk-1)-h
   c(npkk)=c(npkk)+h
   hx=h
   c(2)=c(2)+h
   if (nmt.gt.0) then
      do ii=1,nmt
         indxx=imt(ii)
         if (mtp(indxx).lt.442) c(indxx)=c(indxx)+h
         if (mtp(indxx).eq.442) c(indxx)=c(indxx)-h
      enddo
   endif
   if (iprint.eq.0.or.abs(e-elist(ilist)).gt.small*e) go to 180
   if (ik.lt.nk.and.y.lt.smin) go to 180
   if (nk.eq.1.and.y.lt.smin) go to 180
   test=-1
   test=test/10000
   if (nk.gt.1.and.ik.eq.nk.and.c(npkk).gt.test) go to 180
   ebarp=sigfig(ebar,9,0)
   yp=sigfig(y,9,0)
   hxp=sigfig(hx,9,0)
   hp=sigfig(h,9,0)
   write(nsyso,'(1x,1p,8e14.4)') e,ebarp,yp,-hxp,hp
   if (ik.eq.nk.and.nk.gt.1) then
      subtot=-c(npkk)
      subtot=sigfig(subtot,9,0)
      write(nsyso,'(48x,''subtot   '',1p,e14.4)') subtot
   endif
  180 continue
  195 continue
   it=i
   if (i.eq.ne) it=-it
   call loada(it,c,npkk,inew,bufn,nbuf)
   if (abs(e-elist(ilist)).lt.small*e) ilist=ilist+1
   if (i.lt.ne) go to 190

   !--loop over all subsections and reactions.
   isave=iold
   iold=inew
   inew=isave
   if (ik.lt.nk) go to 130
   if (mfd.eq.13) go to 105
   if (mtx.eq.3.and.mtd.eq.1) go to 200
   go to 105
  200 continue
   mtd=2
   call findf(matd,mfx,mtx,nscr)
   call contio(nscr,0,0,scr,nb,nw)
   go to 110
  210 continue
   call mess('gheat','no file 12 for this material.',' ')
   call findf(matd,13,0,nscr)

   !--finished with present file.
  220 continue
   if (mfd.eq.13) go to 250
   mfd=13
   mtd=0
   go to 105
  250 continue
   if (kchk.ne.1) go to 300
   if (mt303.eq.0) go to 300

   !--print kinematic test for total photon production
   ilist=1
   write(nsyso,'(/'' photon energy production check''/&
     &14x,''e'',6x,''ev-barns'',11x,''min'',11x,''max'')')
   do i=1,ne
      call finda(i,c,npkk,iold,bufo,nbuf)
      e=sigfig(c(1),9,0)
      if (abs(e-elist(ilist)).le.small*e) then
         ilist=ilist+1
         elo=c(npkk-1)+c(mt303)-c(mt303+2*npkt)+c(npkk-2)
         ehi=c(npkk-1)+c(mt303)-c(mt303+npkt)+c(npkk-2)
         jerr=0
         if (c(npkk-1).gt.ehi+ehi/10) jerr=+1
         if (c(npkk-1).lt.elo-elo/10) jerr=-1
         if (jerr.eq.0) write(nsyso,'(1x,1p,8e14.4)')&
           e,c(npkk-1),elo,ehi
         if (jerr.gt.0) write(nsyso,'(1x,1p,4e14.4,''   ++++'')')&
           e,c(npkk-1),elo,ehi
         if (jerr.lt.0) write(nsyso,'(1x,1p,4e14.4,''   ----'')')&
           e,c(npkk-1),elo,ehi
      endif
   enddo
  300 continue
   deallocate(d)
   deallocate(b)
   deallocate(scr)
   deallocate(bufn)
   deallocate(bufo)
   deallocate(sk,sd,sy)
   return
   end subroutine gheat

   subroutine gambar(e,ebar,esqb,esqd,nin,matd,mtd,a,nwamax)
   !-------------------------------------------------------------------
   ! Calculate ebar for continuous spectra and the photon
   ! recoil correction for capture.
   !-------------------------------------------------------------------
   use util ! provides error
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,matd,mtd,nwamax
   real(kr)::e,ebar,esqb,esqd,a(*)
   ! internals
   integer::nb,nw,mfd,l,iz,lf,iend,j,il
   integer::jstart,jend,istart,ilo,ihi,nnt,ne,nne,inn,nbt
   real(kr)::z,awr,elo,ehi,flo,fhi,glo,ghi,hlo,hhi,sen
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::up=1.000001e0_kr
   real(kr),parameter::zero=0
   save lf,z,awr
   save jstart,jend,istart,ilo,ihi,nnt,ne,nne,inn,nbt
   save elo,ehi,flo,fhi,glo,ghi,hlo,hhi

   !--initialize if e=0.
   if (e.gt.zero) go to 110
   mfd=15
   call findf(matd,mfd,mtd,nin)
   l=1
   call contio(nin,0,0,a(l),nb,nw)
   iz=nint(c1h/1000)
   z=iz
   awr=c2h
   call tab1io(nin,0,0,a(l),nb,nw)
   l=l+nw
   do while (nb.ne.0)
      if (l.gt.nwamax) call error('gambar',&
        'storage exceeded in a.',' ')
      call moreio(nin,0,0,a(l),nb,nw)
      l=l+nw
   enddo
   lf=1

   !--read in first incident energy point.
   l=1
   call tab2io(nin,0,0,a(l),nb,nw)
   ne=nint(a(l+5))
   nnt=l+6
   nbt=nint(a(nnt))
   inn=nint(a(nnt+1))
   inn=2
   istart=l+nw
   l=istart
   call tab1io(nin,0,0,a(l),nb,nw)
   l=l+nw
   do while (nb.ne.0)
      if (l.gt.nwamax) call error('gambar',&
        'storage exceeded in a.',' ')
      call moreio(nin,0,0,a(l),nb,nw)
      l=l+nw
   enddo
   iend=l-1
   jstart=iend+1
   jend=0
   ehi=a(istart+1)
   elo=ehi
   nne=1
   ilo=0
   ihi=0
   fhi=0
   ghi=0
   hhi=0
   return

   !--normal entry.  is desired energy in current panel.
  110 continue
   if (e.lt.ehi*(1-small)) go to 140
   if (nne.eq.1) go to 130
   if (nne.eq.ne) go to 160

   !--no.  slide high energy data into low energy positions.
  120 continue
   j=-1
   do il=jstart,jend
      j=j+1
      a(j+istart)=a(il)
   enddo
   iend=istart+j
   jstart=iend+1
   elo=ehi
   flo=fhi
   if (mtd.eq.102) glo=ghi
   if (mtd.eq.102) hlo=hhi
   ilo=ihi

   !--read new high energy data record.
  130 continue
   if (nne.eq.ne)&
     call error('gambar','requested energy gt highest given.',' ')
   l=jstart
   call tab1io(nin,0,0,a(l),nb,nw)
   l=l+nw
   do while (nb.ne.0)
      if (l.gt.nwamax) call error('gambar',&
        'storage exceeded in a.',' ')
      call moreio(nin,0,0,a(l),nb,nw)
      l=l+nw
   enddo
   jend=l-1
   nne=nne+1
   ihi=0
   ehi=a(jstart+1)
   if (ehi.lt.e*(1-small).and.nne.lt.ne) go to 120

   !--integrate at high and low energies.
   if (ilo.le.0) then
      call tabbar(flo,a(istart),lf)
      if (mtd.eq.102) call tabsqr(glo,hlo,a(istart),lf,z,awr)
   endif
   call tabbar(fhi,a(jstart),lf)
   if (mtd.eq.102) call tabsqr(ghi,hhi,a(jstart),lf,z,awr)
   ihi=1

   !--yes.  interpolate for mean energy.
  140 continue
   if (nne.eq.1) go to 150
   if (nne.gt.nbt) then
      nnt=nnt+2
      nbt=nint(a(nnt))
      inn=nint(a(nnt+1))
      inn=2
   endif
   if (elo.eq.ehi) then
       elo=sigfig(elo,7,-1)
       ehi=sigfig(ehi,7,+1)
   endif
   call terp1(elo,flo,ehi,fhi,e,sen,inn)
   if (mtd.eq.102) call terp1(elo,glo,ehi,ghi,e,esqb,inn)
   if (mtd.eq.102) call terp1(elo,hlo,ehi,hhi,e,esqd,inn)
   ebar=sen
   return

   !--return zero below threshold.
  150 continue
   ebar=0
   esqb=0
   esqd=0
   return

   !--return zero or fhi above table.
  160 continue
   if (inn.ne.1.and.e.ge.up*ehi) then
      ebar=0
      esqb=0
      esqd=0
   else
      ebar=fhi
      esqb=ghi
      esqd=hhi
   endif
   return
   end subroutine gambar

   subroutine tabsqr(g,h,a,law,z,awr)
   !-------------------------------------------------------------------
   ! Compute average of photon recoil energy from capture and
   ! corresponding damage energy for a tabulated section of File 15.
   !-------------------------------------------------------------------
   use endf ! provides terp1
   ! externals
   integer::law
   real(kr)::g,h,a(*),z,awr
   ! internals
   integer::nr,np,ibase,ir,nbt,inn,i,j
   real(kr)::ein,rein,xl,yl,xh,yh,dx,x,y,xr,s
   integer,parameter::nq=4
   real(kr),dimension(nq),parameter::qp=(/&
     -.86114e0_kr,-.33998e0_kr,.33998e0_kr,.86114e0_kr/)
   real(kr),dimension(nq),parameter::qw=(/&
     .34785e0_kr,.65215e0_kr,.65215e0_kr,.34785e0_kr/)

   !--initialize
   ein=2*tm
   rein=1/ein
   g=0
   h=0
   s=0
   nr=nint(a(5))
   np=nint(a(6))
   ibase=6+2*nr
   ir=1
   nbt=nint(a(7))
   inn=nint(a(8))
   xh=a(ibase+1)
   yh=a(ibase+2)

   !--accumulate contributions to integrals
   do i=2,np
      xl=xh
      xh=a(ibase+2*i-1)
      yl=yh
      yh=a(ibase+2*i)
      if (i.gt.nbt) then
         ir=ir+1
         nbt=nint(a(6+2*ir-1))
         inn=nint(a(6+2*ir))
      endif
      if (xl.ne.xh) then
         dx=xh-xl
         do j=1,nq
            x=xl+(1+qp(j))*dx/2
            call terp1(xl,yl,xh,yh,x,y,inn)
            xr=x*x*rein
            g=g+qw(j)*xr*dx*y
            h=h+qw(j)*df(xr,z,awr+1,z,awr)*dx*y
            s=s+qw(j)*dx*y
         enddo
      endif
   enddo
   g=g/s
   h=h/s
   return
   end subroutine tabsqr

   subroutine disgam(e,egam,edam,z,awr)
   !-------------------------------------------------------------------
   ! Compute recoil and damage energy of discrete capture photons.
   !-------------------------------------------------------------------
   ! externals
   real(kr)::e,egam,edam,z,awr
   ! internals
   real(kr),parameter::zero=0

   if (e.eq.zero) then
      return
   else
      egam=e*e*rtm/2
      edam=df(egam,z,awr+1,z,awr)
   endif
   return
   end subroutine disgam

   subroutine hout(iold)
   !-------------------------------------------------------------------
   ! Write the output tape with corrected dictionary and
   ! desired MT-s added.
   !-------------------------------------------------------------------
   use util ! provides openz,repoz,closz,error,sigfig
   use mainio ! provides nsyso
   use endf ! provides endf routines and variables
   ! externals
   integer::iold
   ! internals
   integer::nscr,npkk,npkt,npkd,npktd,i,ilist,inpk,nb,nw
   integer::nxn,k,mfi,mti,mtn,ia,nn
   integer::idict,idict1,nx,inow,nwd,n,ilo,ihi,j,ibase
   integer::iowr,mfnin,mtnin,mfnscr,mtnscr,inow6
   real(kr)::x,y,xlast,ylo,yhi,thin,rat,elo,test,e
   real(kr)::c(50)
   integer::ncds(50)
   character(4)::klo(50),khi(50)
   real(kr),dimension(:),allocatable::scr
   real(kr),dimension(:),allocatable::b
   real(kr),dimension(:),allocatable::buf
   character(60)::strng
   character(4)::iblank='    '
   character(4)::idown='low '
   character(4)::iup='high'
   character(1),parameter::qu=''''
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::emax=20.e6_kr
   real(kr),parameter::ecut=.01e6_kr
   real(kr),parameter::zero=0
   real(kr),parameter::ten=10
   integer,parameter::nwscr=10000

   ! allocate scratch storage for endf records and buffer
   allocate(scr(nwscr))
   allocate(buf(nbuf))

   !--construct tab1 records for the kerma factors
   n=0
   nscr=15
   if (nout.lt.0) nscr=-nscr
   call openz(nscr,1)
   call repoz(nscr)
   if (iprint.eq.1) write(nsyso,&
     '(/'' final kerma factors''/14x,''e'',8(11x,i3))')&
     (mtp(i),i=2,npk)
   npkk=npk
   if (kchk.eq.1) npkk=3*npk-2
   if (mgam.gt.0) npkk=npkk+3
   if (mgam.eq.0.and.i6.gt.0) npkk=npkk+1
   npkt=npk-1
   npkd=npk
   npktd=npkt
   do i=2,npk
      if (mtp(i).ge.442) then
         npkd=npkd-1
         npktd=npktd-1
      endif
   enddo
   ilist=1
   nsc=0
   math=matd
   mfh=3
   inpk=1
  105 continue
   inpk=inpk+1
   mth=mtp(inpk)
   scr(1)=za
   scr(2)=awr
   scr(3)=0
   scr(4)=0
   scr(5)=0
   scr(6)=0
   call contio(0,0,nscr,scr,nb,nw)
   nn=0
   scr(1)=0
   scr(2)=0
   scr(5)=1
   scr(6)=0
   scr(7)=0
   scr(8)=2
   ibase=8
   i=0
   nw=0
  110 continue
   nn=nn+1
   if (nn.gt.ne) go to 125
   call finda(nn,c,npkk,iold,buf,nbuf)
   do j=1,npkk
      if (c(j).gt.big) then
         c(j)=sigfig(c(j),6,0)
      else
         c(j)=sigfig(c(j),7,0)
      endif
   enddo
   e=sigfig(c(1),9,0)
   if (iprint.eq.1.and.abs(e-elist(ilist)).le.small*e) then
      ihi=0
      ilo=0
      if (kchk.eq.1) then
         ihi=0
         ilo=0
         do j=2,npkd
            read(iblank,'(a4)') klo(j)
            test=c(j+npkt)-c(j+npkt)/10
            if (c(j).lt.test) then
               read(idown,'(a4)') klo(j)
               ilo=1
            endif
            read(iblank,'(a4)') khi(j)
            test=c(j+2*npkt)+c(j+2*npkt)/10
            if (c(j).gt.test) then
               read(iup,'(a4)') khi(j)
               ihi=1
            endif
         enddo
      endif
      if (ilo.eq.1) write(nsyso,'('' '')')
      if (kchk.eq.1) write(nsyso,'(15x,8(6x,a4,4x))')&
        (klo(j),j=2,npkd)
      if (kchk.eq.1) write(nsyso,'(12x,''min'',1p,8e14.4)')&
        (c(npk+j),j=1,npktd)
      write(nsyso,'(1x,1p,9e14.4)') (c(j),j=1,npk)
      if (kchk.eq.1) write(nsyso,'(12x,''max'',1p,8e14.4)')&
        (c(npk+npkt+j),j=1,npktd)
      if (ihi.eq.1) write(nsyso,'(15x,8(6x,a4,4x))')&
        (khi(j),j=2,npkd)
   endif
   if (e.ge.elist(ilist)) ilist=ilist+1
   if (i.gt.0) go to 118
   if (c(inpk).eq.zero.and.i.eq.0.and.ibase.eq.8) go to 110
   if (nw.gt.0) go to 118
   n=ne-nn+1
   scr(6)=n
   scr(7)=n
  118 continue
   i=i+2
   scr(ibase+i-1)=c(1)
   scr(ibase+i)=c(inpk)
   if (i.lt.npage.and.nn.ne.ne) go to 110
   if (ibase.eq.0) go to 120
   call tab1io(0,0,nscr,scr,nb,nw)
   if (nn.eq.ne) go to 130
   ibase=0
   i=0
   go to 110
  120 continue
   call moreio(0,0,nscr,scr,nb,nw)
   if (nn.eq.ne) go to 130
   i=0
   go to 110
   ! all zeroes for this mt
  125 continue
   scr(6)=2
   scr(7)=2
   scr(9)=efirst
   scr(10)=0
   scr(11)=elast
   scr(12)=0
   nw=12
   call tab1io(0,0,nscr,scr,nb,nw)
   n=2
  130 continue
   call asend(0,nscr)
   ncds(inpk)=3+(n+2)/3
   if (inpk.lt.npk) go to 105
   call afend(0,nscr)
   call repoz(nscr)

   !--update dictionary in file 1
   call findf(matd,1,451,nin)
   call contio(nin,0,0,scr,nb,nw)
   nx=n2h
   inow=7
   if (iverf.ge.5) then
      call contio(nin,0,0,scr(inow),nb,nw)
      inow=inow+6
      if (iverf.eq.6) then
         call contio(nin,0,0,scr(inow),nb,nw)
         inow=inow+6
      endif
   endif
   call hdatio(nin,0,0,scr(inow),nb,nw)
   nwd=nint(scr(inow+4))
   if (iverf.ge.5) nx=nint(scr(inow+5))
   inow=inow+nw
   do while (nb.ne.0)
      call moreio(nin,0,0,scr(inow),nb,nw)
      inow=inow+nw
   enddo
   nw=nx
   idict=inow
   idict1=idict-1
   call dictio(nin,0,0,scr(idict),nb,nw)
   nw=inow+nw
   nw=6*(nx+20)
   if (nw.lt.npage+50) nw=npage+50
   allocate(b(nw))
   inpk=2
   mtn=mtp(inpk)
   j=0

   !--loop over existing dictionary records
   do i=1,nx
      iowr=0
      ia=6*(i-1)
      mfi=nint(scr(idict1+ia+3))
      mti=nint(scr(idict1+ia+4))
      !--insert heating mt dictionary data (or overwrite
      !  original dictionary data with revised data).
      do while ((mfi.eq.3 .and. mtn.le.mti) .or.&
                (mfi.gt.3 .and. inpk.le.npk))
         if (mfi.eq.3 .and. mtn.eq.mti) iowr=1
         b(j+1)=0
         b(j+2)=0
         b(j+3)=3
         b(j+4)=mtp(inpk)
         b(j+5)=ncds(inpk)
         b(j+6)=0
         j=j+6
         inpk=inpk+1
         if (inpk.le.npk) then
            mtn=mtp(inpk)
         else
            mtn=10000
         endif
      enddo
      !--insert original dictionary data (unless the heating
      !  mt entry above was an overwrite of existing data).
      if (mfi.ne.3 .or. (mfi.eq.3 .and. iowr.eq.0)) then
         b(j+1)=0
         b(j+2)=0
         b(j+3)=mfi
         b(j+4)=mti
         b(j+5)=scr(idict1+ia+5)
         b(j+6)=scr(idict1+ia+6)
         j=j+6
      endif
   enddo
   !--insert rest of heating mt dictionary data, if any.
   do i=inpk,npk
      b(j+1)=0
      b(j+2)=0
      b(j+3)=3
      b(j+4)=mtp(i)
      b(j+5)=ncds(i)
      b(j+6)=0
      j=j+6
   enddo
   mth=451
   nxn=j/6
   inow=1
   if (iverf.le.4) scr(inow+5)=nxn
   call contio(0,nout,0,scr(inow),nb,nw)
   inow=inow+6
   if (iverf.eq.6) call contio(0,nout,0,scr(inow),nb,nw)
   if (iverf.eq.6) inow=inow+6
   if (iverf.gt.4) then
      call contio(0,nout,0,scr(inow),nb,nw)
      inow=inow+6
      scr(inow+5)=nxn
   endif
   nw=6+nwd
   call hdatio(0,nout,0,scr(inow),nb,nw)
   inow=inow+nw
   do while (nb.ne.0)
      call moreio(0,nout,0,scr(inow),nb,nw)
      inow=inow+nw
   enddo
   nw=nxn
   call dictio(0,nout,0,b,nb,nw)
   call tofend(nin,nout,0,scr)
   call tofend(nin,nout,0,scr)

   !--add new reactions to file 3.
   call contio(nin,0,0,scr,nb,nw)
   if (mfh.ne.3 .or. mth.eq.0) then
      write(strng,'("nin out of order.  read mfh,mth = ",2i5)')mfh,mth
      call error('hout',strng,' ')
   endif
   mfnin=mfh
   mtnin=mth
   inow=7
   call contio(nscr,0,0,scr(inow),nb,nw)
   if (mfh.ne.3 .or. mth.eq.0) then
      write(strng,'("nscr out of order.  read mfh,mth = ",2i5)')mfh,mth
      call error('hout',strng,' ')
   endif
   mfnscr=mfh
   mtnscr=mth
   inow6=inow+6
   do while (mfnin.eq.3 .and. mfnscr.eq.3)
      if (mtnin.lt.mtnscr) then
         mfh=mfnin
         mth=mtnin
         call contio(0,nout,0,scr,nb,nw)
         call tosend(nin,nout,0,scr(inow6))
         call contio(nin,0,0,scr,nb,nw)
         mfnin=mfh
         mtnin=mth
      elseif (mtnin.eq.mtnscr) then
         mfh=mfnscr
         mth=mtnscr
         call contio(0,nout,0,scr(inow),nb,nw)
         call tosend(nscr,nout,0,scr(inow6))
         call tosend(nin,0,0,scr(inow6))
         call contio(nin,0,0,scr,nb,nw)
         mfnin=mfh
         mtnin=mth
         call contio(nscr,0,0,scr(inow),nb,nw)
         mfnscr=mfh
         mtnscr=mth
      else
         do while (mtnin.ge.mtnscr .and. mfnscr.eq.3)
            mfh=mfnscr
            mth=mtnscr
            call contio(0,nout,0,scr(inow),nb,nw)
            call tosend(nscr,nout,0,scr(inow6))
            call contio(nscr,0,0,scr(inow),nb,nw)
            mfnscr=mfh
            mtnscr=mth
         enddo
      endif
   enddo
   if (mfnin.eq.3 .and. mfnscr.eq.0) then
      mfh=mfnin
      mth=mtnin
      call contio(0,nout,0,scr,nb,nw)
   elseif (mfnin.eq.0 .and. mfnscr.eq.3) then
      mfh=mfnscr
      mth=mtnscr
      call contio(0,nout,0,scr(inow),nb,nw)
      call tofend(nscr,nout,0,scr(inow6))
   elseif (mfnin.eq.0 .and. mfnscr.eq.0) then
      mfh=mfnin
      mth=mtnin
      call contio(0,nout,0,scr,nb,nw)
   endif
   call tomend(nin,nout,0,scr)
   deallocate(b)
   deallocate(scr)
   call closz(nscr)

   !--plot heating, photon energy production, and kinematic limits
   if (nplot.eq.0) go to 500
   write(nplot,'(''/'')')

   !--total heating
   !--log-log
   elo=1
   elo=elo/1000
   rat=emax/elo
   thin=ten**(log10(rat)/2500)
   ylo=big
   yhi=-big
   xlast=elo
   do i=1,ne
      call finda(i,c,npkk,iold,buf,nbuf)
      x=c(1)
      if (x.ge.thin*xlast) then
         if (c(2).gt.yhi) yhi=c(2)
         if (c(2).lt.ylo) ylo=c(2)
         if (c(2+npkt).gt.yhi) yhi=c(2+npkt)
         if (c(2+npkt).lt.ylo) ylo=c(2+npkt)
         if (c(2+2*npkt).gt.yhi) yhi=c(2+2*npkt)
         if (c(2+2*npkt).lt.ylo) ylo=c(2+2*npkt)
         xlast=x
      endif
   enddo
   if (ylo.le.zero) ylo=small
   if (yhi.le.zero) yhi=small
   do k=1,3
      if (k.eq.1) then
         write(nplot,'(''1/'')')
         write(nplot,'(a,''Energy-Balance Check'',a,''/'')') qu,qu
         write(nplot,'(''/'')')
         write(nplot,'(''4 0 2 1/'')')
         write(nplot,'(''/'')')
         write(nplot,'(a,''<e>nergy (e<v>)'',a,''/'')') qu,qu
         write(nplot,'(''/'')')
         write(nplot,'(a,''<h>eating (e<v>-barns)'',a,''/'')') qu,qu
         write(nplot,'(''/'')')
         write(nplot,'(''0 0 4/'')')
         write(nplot,'(a,''<u>pper limit'',a,''/'')') qu,qu
      else if (k.eq.2) then
         write(nplot,'(''2/'')')
         write(nplot,'(''/'')')
         write(nplot,'(''/'')')
         write(nplot,'(a,''<e>nergy-balance heating'',a,''/'')') qu,qu
      else
         write(nplot,'(''3/'')')
         write(nplot,'(''/'')')
         write(nplot,'(''0 0 4/'')')
         write(nplot,'(a,''<l>ower limit'',a,''/'')') qu,qu
      endif
      write(nplot,'(''0/'')')
      xlast=elo
      j=0
      do i=1,ne
         call finda(i,c,npkk,iold,buf,nbuf)
         x=c(1)
         if (k.eq.1) y=c(2+2*npkt)
         if (k.eq.2) y=c(2)
         if (k.eq.3) y=c(2+npkt)
         if (x.ge.thin*xlast.and.j.lt.2500) then
            j=j+1
            if (y.le.zero) y=small
            write(nplot,'(1p,2e14.6,''/'')') x,y
            xlast=x
         endif
      enddo
      if (k.eq.1) write(nplot,'(1p,2e14.6,''/'')') x,ylo
      if (k.eq.1) write(nplot,'(1p,2e14.6,''/'')') x,yhi
      write(nplot,'(''/ end'')')
   enddo

   !--total heating
   !--lin-lin
   thin=emax/2500
   ylo=big
   yhi=-big
   xlast=ecut
   do i=1,ne
      call finda(i,c,npkk,iold,buf,nbuf)
      x=c(1)
      if (x.ge.xlast+thin) then
         if (c(2).gt.yhi) yhi=c(2)
         if (c(2).lt.ylo) ylo=c(2)
         if (c(2+npkt).gt.yhi) yhi=c(2+npkt)
         if (c(2+npkt).lt.ylo) ylo=c(2+npkt)
         if (c(2+2*npkt).gt.yhi) yhi=c(2+2*npkt)
         if (c(2+2*npkt).lt.ylo) ylo=c(2+2*npkt)
         xlast=x
      endif
   enddo
   do k=1,3
      if (k.eq.1) then
         write(nplot,'(''1/'')')
         write(nplot,'(a,''Energy-Balance Check'',a,''/'')') qu,qu
         write(nplot,'(''/'')')
         write(nplot,'(''1 0 2 1/'')')
         write(nplot,'(''/'')')
         write(nplot,'(a,''<e>nergy (e<v>)'',a,''/'')') qu,qu
         write(nplot,'(''/'')')
         write(nplot,'(a,''<h>eating (e<v>-barns)'',a,''/'')') qu,qu
         write(nplot,'(''/'')')
         write(nplot,'(''0 0 4/'')')
         write(nplot,'(a,''<u>pper limit>'',a,''/'')') qu,qu
      else if (k.eq.2) then
         write(nplot,'(''2/'')')
         write(nplot,'(''/'')')
         write(nplot,'(''/'')')
         write(nplot,'(a,''<e>nergy-balance heating'',a,''/'')') qu,qu
      else
         write(nplot,'(''3/'')')
         write(nplot,'(''/'')')
         write(nplot,'(''0 0 4/'')')
         write(nplot,'(a,''<l>ower limit'',a,''/'')')qu,qu
      endif
      write(nplot,'(''0/'')')
      xlast=ecut
      j=0
      do i=1,ne
         call finda(i,c,npkk,iold,buf,nbuf)
         x=c(1)
         if (k.eq.1) y=c(2+2*npkt)
         if (k.eq.2) y=c(2)
         if (k.eq.3) y=c(2+npkt)
         if (x.ge.xlast+thin.and.j.lt.2500) then
            j=j+1
            write(nplot,'(1p,2e14.6,''/'')') x,y
            xlast=x
         endif
      enddo
      if (k.eq.1) write(nplot,'(1p,2e14.6,''/'')') x,ylo
      if (k.eq.1) write(nplot,'(1p,2e14.6,''/'')') x,yhi
      write(nplot,'(''/ end'')')
   enddo

   !--photon energy production
   !--log-log
   if (mt303.eq.0) go to 490
   elo=1
   elo=elo/1000
   rat=emax/elo
   thin=ten**(log10(rat)/2500)
   ylo=big
   yhi=-big
   xlast=elo
   do i=1,ne
      call finda(i,c,npkk,iold,buf,nbuf)
      x=c(1)
      if (x.ge.thin*xlast) then
         y=c(npkk-1)+c(mt303)-c(mt303+npkt)+c(npkk-2)
         if (y.gt.yhi) yhi=y
         if (y.lt.ylo) ylo=y
         y=c(npkk-1)
         if (y.gt.yhi) yhi=y
         if (y.lt.ylo) ylo=y
         y=c(npkk-1)+c(mt303)-c(mt303+2*npkt)+c(npkk-2)
         if (y.gt.yhi) yhi=y
         if (y.lt.ylo) ylo=y
         xlast=x
      endif
   enddo
   if (ylo.le.zero) ylo=small
   if (yhi.le.zero) yhi=small
   do k=1,3
      if (k.eq.1) then
         write(nplot,'(''1/'')')
         write(nplot,'(a,''Energy-Balance Check'',a,''/'')') qu,qu
         write(nplot,'(''/'')')
         write(nplot,'(''4 0 2 1/'')')
         write(nplot,'(''/'')')
         write(nplot,'(a,''<e>nergy (e<v>)'',a,''/'')') qu,qu
         write(nplot,'(''/'')')
         write(nplot,'(a,''<p>hoton energy prod (e<v>-barns)'',&
           &a,''/'')') qu,qu
         write(nplot,'(''/'')')
         write(nplot,'(''0 0 4/'')')
         write(nplot,'(a,''<u>pper limit'',a,''/'')') qu,qu
      else if (k.eq.2) then
         write(nplot,'(''2/'')')
         write(nplot,'(''/'')')
         write(nplot,'(''/'')')
         write(nplot,'(a,''<p>hoton energy production'',a,''/'')')&
           qu,qu
      else
         write(nplot,'(''3/'')')
         write(nplot,'(''/'')')
         write(nplot,'(''0 0 4/'')')
         write(nplot,'(a,''<l>ower limit'',a,''/'')') qu,qu
      endif
      write(nplot,'(''0/'')')
      xlast=elo
      j=0
      do i=1,ne
         call finda(i,c,npkk,iold,buf,nbuf)
         x=c(1)
         if (k.eq.1) y=c(npkk-1)+c(mt303)-c(mt303+npkt)+c(npkk-2)
         if (k.eq.2) y=c(npkk-1)
         if (k.eq.3) y=c(npkk-1)+c(mt303)-c(mt303+2*npkt)+c(npkk-2)
         if (x.ge.thin*xlast.and.j.lt.2500) then
            j=j+1
            if (y.le.zero) y=small
            write(nplot,'(1p,2e14.6,''/'')') x,y
            xlast=x
         endif
      enddo
      if (k.eq.1) write(nplot,'(1p,2e14.6,''/'')') x,ylo
      if (k.eq.1) write(nplot,'(1p,2e14.6,''/'')') x,yhi
      write(nplot,'(''/ end'')')
   enddo

   !--photon energy production
   !--lin-lin
   thin=emax/2500
   ylo=big
   yhi=-big
   xlast=ecut
   do i=1,ne
      call finda(i,c,npkk,iold,buf,nbuf)
      x=c(1)
      if (x.ge.xlast+thin) then
         y=c(npkk-1)+c(mt303)-c(mt303+npkt)+c(npkk-2)
         if (y.gt.yhi) yhi=y
         if (y.lt.ylo) ylo=y
         y=c(npkk-1)
         if (y.gt.yhi) yhi=y
         if (y.lt.ylo) ylo=y
         y=c(npkk-1)+c(mt303)-c(mt303+2*npkt)+c(npkk-2)
         if (y.gt.yhi) yhi=y
         if (y.lt.ylo) ylo=y
         xlast=x
      endif
   enddo
   do k=1,3
      if (k.eq.1) then
         write(nplot,'(''1/'')')
         write(nplot,'(a,''Energy-Balance Check'',a,''/'')') qu,qu
         write(nplot,'(''/'')')
         write(nplot,'(''1 0 2 1/'')')
         write(nplot,'(''/'')')
         write(nplot,'(a,''<e>nergy (e<v>)'',a,''/'')') qu,qu
         write(nplot,'(''/'')')
         write(nplot,'(a,''<p>hoton energy prod (e<v>-barns)'',&
           &a,''/'')') qu,qu
         write(nplot,'(''/'')')
         write(nplot,'(''0 0 4/'')')
         write(nplot,'(a,''<u>pper limit'',a,''/'')') qu,qu
      else if (k.eq.2) then
         write(nplot,'(''2/'')')
         write(nplot,'(''/'')')
         write(nplot,'(''/'')')
         write(nplot,'(a,''<p>hoton energy production'',a,''/'')')&
           qu,qu
      else
         write(nplot,'(''3/'')')
         write(nplot,'(''/'')')
         write(nplot,'(''0 0 4/'')')
         write(nplot,'(a,''<l>ower limit'',a,''/'')') qu,qu
      endif
      write(nplot,'(''0/'')')
      xlast=ecut
      j=0
      do i=1,ne
         call finda(i,c,npkk,iold,buf,nbuf)
         x=c(1)
         if (k.eq.1) y=c(npkk-1)+c(mt303)-c(mt303+npkt)+c(npkk-2)
         if (k.eq.2) y=c(npkk-1)
         if (k.eq.3) y=c(npkk-1)+c(mt303)-c(mt303+2*npkt)+c(npkk-2)
         if (x.ge.xlast+thin.and.j.lt.2500) then
            j=j+1
            write(nplot,'(1p,2e14.6,''/'')') x,y
            xlast=x
         endif
      enddo
      if (k.eq.1) write(nplot,'(1p,2e14.6,''/'')') x,ylo
      if (k.eq.1) write(nplot,'(1p,2e14.6,''/'')') x,yhi
      write(nplot,'(''/ end'')')
   enddo

   !--plotting finished
  490 continue
   write(nplot,'(''99/'')')
   call closz(nplot)
  500 continue
   deallocate(buf)
   return
   end subroutine hout

end module heatm

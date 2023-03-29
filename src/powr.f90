module powm
   ! provides subroutine powr for NJOY2016
   use locale
   implicit none
   private
   public powr

   ! i/o units
   integer::ngendf,nout,nscr,kscr,mscr

   ! powr global variables
   integer::lib,ngnf,ngn,ngnd
   real(kr)::cflux(140)

   ! scratch area for reading endf files
   integer::nwscr=5000
   real(kr)::scr(5000)

   ! epri-cell library global variables
   real(kr)::tmpr(10),sigz(10)
   integer::ntp,nsigz,jwf
   integer::iwr,izref
   real(kr)::rtemp
   integer::matd,nid,nff,iwa,iwf,nfs
   real(kr)::xla(3),xld(3),xlol(3)

   ! epri-cpm library global variables
   integer::nfg,nrg,nfsg,ngref,newmat,iprint
   integer::nzmax,mode
   integer::if5,ntis,nfis,if4
   integer::nnina,nres
   integer::iu235,iu,ipu239,ipu,nfiss(2)
   real(kr)::rleth(19),reflth
   integer,dimension(:),allocatable::mat
   real(kr),dimension(:),allocatable::flux

   integer,dimension(:),allocatable::idnt
   integer,dimension(:),allocatable::ninat
   integer,dimension(:),allocatable::irest
   real(kr),dimension(:),allocatable::s0rf
   integer,dimension(:),allocatable::ip1tab
   integer,dimension(:),allocatable::norft
   integer,dimension(:),allocatable::mtit
   integer,dimension(:),allocatable::mtct
   real(kr),dimension(:),allocatable::awrt
   integer,dimension(:),allocatable::indf
   integer,dimension(:),allocatable::ntmp
   integer,dimension(:),allocatable::nszt
   real(kr),dimension(:),allocatable::glamt
   integer,dimension(:),allocatable::ipost
   integer,dimension(:),allocatable::iposrt
   real(kr),dimension(:),allocatable::egb
   real(kr),dimension(:),allocatable::sigpt
   real(kr),dimension(:),allocatable::uff,puff
   real(kr),dimension(:),allocatable::fiss
   real(kr),dimension(:),allocatable::yld
   integer,dimension(:),allocatable::idat,idbt
   real(kr),dimension(:),allocatable::elas
   real(kr),dimension(:),allocatable::tempt

contains

   subroutine powr
   !--------------------------------------------------------------------
   !
   !  Produce input for the EPRI-CELL codes GAMTAP (fast) and
   !  LIBRAR (thermal), and the EPRI-CPM code CLIB.
   !
   !---input specifications (free format)---------------------------
   !
   ! card 1
   !    ngendf  unit for input gout tape
   !    nout    unit for output tape
   ! card 2
   !    lib     library option (1=fast, 2=thermal, 3=cpm)
   !    iprint  print option (0=minimum, 1=maximum)
   !            (default=0)
   !    iclaps  group collapsing option (0=collapse from 185 group
   !            to desired group structure, 1=no collapse)
   !            (default=0)
   !
   !---for lib=1----------------------------------------------------
   !
   ! card 3
   !    matd    material to be processed
   !            if matd lt 0, read-in absorption data only for
   !            this material with mat=abs(matd) directly from
   !            input deck (see card 6)
   !   following three parameters irrelevant for matd lt 0
   !    rtemp   reference temperature (degrees kelvin)
   !            (default=300  k)
   !    iff     f-factor option
   !            (0/1=do not calculate f-factors/calculate if found)
   !            (default=1)
   !    nsgz    no. of sigma zeroes to process for this material
   !            (default=0=all found on input tape)
   !    izref   ref. sigzero for elastic matrix (default=1)
   ! cards 4 and 5 for normal run only (matd gt 0)
   ! card 4
   !    word    description of nuclide (up to 16 characters,
   !            delimited with *, ended with /)  (default=blank)
   ! card 5
   !    fsn     title of fission spectrum (up to 40 characters,
   !            delimited with *, ended with /0  (default=blank)
   !             delimited with *, ended with /)  (default=blank)
   ! card 6 for reading in absorption data only
   !    abs     ngnd absorption values (default values=0)
   ! repeat cards 3 through 6 for each material desired.
   ! terminate with matd=0/ (i.e., a 0/ card).
   !
   !---for lib=2----------------------------------------------------
   !
   ! card 3
   !    matd    material to be processed
   !    idtemp  temperature id (default=300  k)
   !    name    hollerith name of isotope (up to 10 characters,
   !            delimited with *, ended with /)  (default=blank)
   ! card 4     default for all values=0.
   !    itrc    transport correction option (0 no, 1 yes)
   !    mti     thermal inelastic mt
   !    mtc     thermal elastic mt
   ! card 5     default for all values=0.
   !    xi
   !    alpha
   !    mubar
   !    nu
   !    kappa fission
   !    kappa capture
   !    lambda
   !    sigma s   if 0, set to scattering cross section at group 35
   ! repeat cards 3 thru 5 for each material and temperature desired
   ! (maximum number of temperatures allowed is 7.)
   ! terminate with matd=0/ (i.e., a 0/ card).
   !
   !---for lib=3----------------------------------------------------
   !
   ! card 3
   !    nlib    number of library.
   !    idat    date library is written (i format).
   !    newmat  number of materials to be added.
   !    iopt    add option (0=mats will be read in,
   !                        1=use all mats found on ngendf).
   !    mode    0/1/2=replace isotope(2) in cpmlib/
   !                   add/create a new library (default=0)
   !    if5     file5 (burnup data) option
   !            0/1/2=do not process file5 burnup data/
   !                   process burnup data along with rest of data/
   !                   process burnup data only (default=0)
   !            (default=0)
   !    if4     file4 (cross section data) option
   !            0/1=do not process/process
   !            (default=1)
   ! card 4 for iopt=0 only
   !    mat     endf mat number of all desired materials.
   !            for materials not on gendf tape, use ident for mat.
   !            if mat lt 0, add 100 to output ident
   !            (for second isomer of an isotope)
   ! card 5
   !    nina    nina indicator.
   !            0/1/2/3=normal/
   !                  no file2 data, calculate absorption in file4/
   !                  no file2 data, read in absorption in file4/
   !                  read in all file2 and file4 data.
   !    ntemp   no. of temperatures to process for this material
   !            (default=0=all found on input tape)
   !    nsigz   no. of sigma zeroes to process for this material
   !            (default=0=all found on input tape)
   !    sgref   reference sigma zero
   !    following 2 parameters are for nina=0 or nina=3.
   !    ires    resonance absorber indicator (0/1=no/yes)
   !    sigp    potential cross section from endf.
   !    following 5 parameters are for ntapea=0 only
   !    mti     thermal inelastic mt
   !    mtc     thermal elastic mt
   !    ip1opt  0/1=calculate p1 matrices/
   !                correct p0 scattering matrix ingroups.
   !  ******if a p1 matrix is calculated for one of the isotopes
   !        having a p1 matrix on the old library, file 6 on the
   !        new library will be completely replaced.******
   !    inorf   0/1=include resonance fission if found/
   !                do not include
   !    following two parameters for mode=0 only
   !    pos     position of this isotope in cpmlib
   !    posr    (for ires=1) position of this isotope in resonance
   !            tabulation in cpmlib
   !  repeat card 5 for each nuclide.
   ! following three cards are for if5 gt 0 only
   ! card 6
   !    ntis    no. time-dependent isotopes
   !    nfis    no. fissionable burnup isotopes
   ! card 7
   !    identb  ident of each of the nfis isotopes
   ! card 8
   !    identa  ident of time-dependent isotope
   !    decay   decay constant (default=0.)
   !    yield   nfis yields (default=0.)
   !  repeat card 8 for each of the ntis isotopes.
   ! card 9 for if5=2 only
   !    aw      atomic weight
   !    indfis  fission indicator
   !    ntemp   no. temperatures on old library
   !  repeat card 9 for each of the ntis isotopes.
   ! card 10
   !    lambda  resonance group goldstein lambdas
   !  ******remember that the 69-group structure has 13 resonance
   !        groups while the collapsed 185-group structure has 15.
   !        use a slash at end of each line of card 10 input.******
   !  repeat card 10 for each nuclide having nina=0, nina=3, or
   !                     ires=1.
   ! cards 11 and 11a for nuclides having nina=3 only.
   ! card 11
   !    resnu    nrg nus values to go with the lambda values
   ! card 11a
   !    tot      nrg total xsec values to go with the lambda values
   !  read cards 11 and 11a for each nuclide having nina=3.
   ! cards 12 for nina gt 2 only
   !    aw      atomic weight
   !    temp    temperature
   !    fpa     ngnd absorption values (default=0.)
   ! cards 12a, 12b, 12c for nuclides having nina=3 only.
   ! card 12a
   !    nus      ngnd nus values
   !    fis      ngnd fission values
   !    xtr      ngnd transport values
   ! card 12b
   !    ia       group.  0 means no scattering from this group
   !    l1       lowest group to which scattering occurs
   !    l2       highest group to which scattering occurs
   ! card 12c  for ia gt 0 only
   !    scat     l2-l1+1 scattering values
   !  repeat card 12b and 12c for each group
   ! repeat cards 12 for each of the nina gt 2 nuclides
   !
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use util ! provides timer,mess,openz,repoz,error,closz
   ! internals
   integer::igprnt,iprnt,iclaps
   real(kr)::sec

   !--initialize.
   ngnf=185
   igprnt=0
   call timer(sec)
   write(nsyso,&
     '(/'' powr...produce input for epri-cell libraries'',&
     &23x,f8.1,''s'')') sec
   write(nsyse,'(/'' powr...'',61x,f8.1,''s'')') sec

   !--read and write common input.
   read(nsysi,*) ngendf,nout
   iprnt=0
   iclaps=0
   read(nsysi,*) lib,iprnt,iclaps
   if (nout.lt.0)&
     call mess('powr','nout must be positive.  changed.',' ')
   if (nout.lt.0) nout=-nout
   call openz(ngendf,0)
   call openz(nout,1)
   write(nsyso,'(/&
     &'' input gendf unit ..................... '',i10/&
     &'' output unit .......................... '',i10/&
     &'' library (1 fast, 2 thermal, 3 cpm) ... '',i10/&
     &'' print option (0 min, 1 max) .......... '',i10/&
     &'' grp collapsing option (0 yes, 1 no) .. '',i10)')&
     ngendf,nout,lib,iprnt,iclaps
   call repoz(nout)

   !--branch for desired library

   !--calculate fast neutrons.
   if (lib.eq.1) then
      call fast(iprnt,igprnt)

   !--thermal neutrons.
   else if (lib.eq.2) then
      call therm(iprnt,igprnt)

   !--cpm library.
   else if (lib.eq.3) then
      call cpm(iprnt,iclaps)

   !--illegal lib
   else
      call error('powr','illegal library requested.',' ')
   endif

   call repoz(ngendf)
   call repoz(nout)
   call closz(ngendf)
   call closz(nout)
   call timer(sec)
   write(nsyso,'(69x,f8.1,''s''/&
     &1x,7(''**********''),''*******'')') sec
   return
   end subroutine powr

   integer function icgrp(ig)
   !--------------------------------------------------------------------
   ! For the group ig in the 185-group fine structure, return the
   ! coarse group for the fast, thermal, or CPM group structure.
   !--------------------------------------------------------------------
   ! externals
   integer::ig
   ! internals
   integer,dimension(186),parameter::igf=(/&
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,2,2,3,3,&
     4,4,4,4,4,4,4,4,4,5,5,5,6,6,6,6,7,7,7,7,7,7,7,8,8,9,9,10,10,&
     11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20,20,&
     21,21,22,22,23,23,24,24,25,25,26,26,27,27,28,28,29,29,30,30,&
     31,31,32,32,33,33,34,34,35,35,36,36,37,37,38,38,39,39,40,40,&
     41,41,42,42,43,43,44,44,45,45,46,46,47,47,48,48,49,49,50,50,&
     51,51,52,52,53,53,54,54,55,55,56,56,57,57,58,58,59,59,60,60,&
     61,61,62,62,63,63,64,64,65,65,66,66,67,67,68,68,69,0,0,0,0,&
     0,0,0,0,0,0,0,0/)
   integer,dimension(186),parameter::igt=(/&
     0,1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,19,20,&
     21,21,22,22,23,23,24,24,25,26,27,28,29,30,30,31,32,32,33,34,&
     35,36,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
   integer,dimension(186),parameter::igc=(/&
     1,2,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,20,21,22,22,&
     23,23,24,25,26,27,28,29,29,29,30,31,32,33,34,34,34,35,35,36,36,&
     36,37,37,37,37,38,39,39,40,41,41,42,42,42,42,43,43,44,44,44,45,&
     45,45,45,46,46,46,46,46,47,47,47,47,48,48,48,48,48,49,49,49,49,&
     49,49,49,50,50,50,50,50,50,50,51,51,51,51,51,52,52,52,53,53,53,&
     54,54,54,54,55,55,55,55,56,56,56,56,57,57,57,57,58,58,58,58,59,&
     59,59,59,60,60,60,60,61,61,61,61,62,62,62,62,63,63,63,63,64,64,&
     64,64,65,65,65,65,66,66,66,66,67,67,67,67,68,68,68,68,69,69,69,&
     69,70,0,0,0,0,0,0,0,0,0,0,0,0/)

   !--no collapsing is required
   if (ngn.eq.ngnd) then
      icgrp=ig
      return
   endif

   !--input group is zero
   if (ig.eq.0) then
      icgrp=0
      return
   endif

   !--fast
   if (lib.eq.1) then
      icgrp=igf(ig)

   !--thermal
   else if (lib.eq.2) then
      icgrp=igt(ig)

   !--cpm
   else if (lib.eq.3) then
      icgrp=igc(ig)
   endif
   return
   end function icgrp

   subroutine fast(iprint,igprnt)
   !--------------------------------------------------------------------
   ! Format multigroup cross sections from groupr for the
   ! EPRI-CELL-cell fast library maintenance code GAMTAP.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use util ! provides openz,repoz,closz
   ! externals
   integer::iprint,igprnt
   ! internals
   integer::ifiss,ngmin,ngmax,iff,nsgz,iread,i,l,lim
   integer::ifs,limsf,limc,la,ld,istart,nwmisc,na
   real(kr)::xld1,xld2,xld3
   real(4)::fsn(10)
   real(kr)::z(70),word(4)
   character(12)::hid
   character(40)::fsp
   integer::locab0,locsf0,locchi,locnus,loce0,loce1,locin,loc2n
   integer::locabs,locsf,loc(10)
   equivalence (locab0,loc(1))
   equivalence (locsf0,loc(2))
   equivalence (locchi,loc(3))
   equivalence (locnus,loc(4))
   equivalence (loce0,loc(5))
   equivalence (loce1,loc(6))
   equivalence (locin,loc(7))
   equivalence (loc2n,loc(8))
   equivalence (locabs,loc(9))
   equivalence (locsf,loc(10))
   real(kr),dimension(:),allocatable::a
   integer::misc=1
   integer::izero=0
   real(kr),parameter::zero=0

   ! set up storage area for fast data
   na=10000
   allocate(a(na))

   ! set up i/o units
   nscr=10
   kscr=11
   call openz(nscr,1)
   call openz(kscr,1)
   call repoz(nscr)
   call repoz(kscr)

   ! extract the fast data
   ifiss=0
   ngnd=68
   ngmin=19
   ngmax=62
   write(nsyso,'(/'' ***epri-cell fast library***'')')
   nff=1
   matd=1
   do while (matd.gt.0)
      rtemp=300
      iff=1
      nsgz=0
      izref=1
      read(nsysi,*) matd,rtemp,iff,nsgz,izref
      if (matd.gt.0) then
         iread=0
         if (matd.lt.0) iread=1
         matd=iabs(matd)
         hid=' '
         read(nsysi,*) hid
         read(hid,'(4a4)') (word(i),i=1,4)
         if (iread.ne.1) then
            fsp=' '
            read(nsysi,*) fsp
            read(fsp,'(10a4)') (fsn(i),i=1,10)
         else
            nid=matd
            iwa=1
            iwr=0
            iwf=0
            do i=1,3
               xlol(i)=0
               xla(i)=0
            enddo
            xld1=0
            xld2=0
            xld3=0
            do i=1,ngnd
               z(i)=0
            enddo
            read(nsysi,*) (z(i),i=1,ngnd)
            write(nsyso,'(&
              &'' mat to be processed .................. '',i10)')&
              matd
            write(nsyso,'(/&
              &'' description of nuclide''/&
              &1x,9(''----''),''--''/6x,4a4)') (word(i),i=1,4)
            write(nsyso,&
              '(/'' absorption only for this material.  '',&
              &''nid=matd.'')')
         endif

         !--write out input parameters.
         if (iread.ne.1) then
            write(nsyso,'(&
              &'' mat to be processed .................. '',i10)')&
              matd
            write(nsyso,'(&
              &'' reference temperature ................ '',1p,e10.3/&
              &'' f-factor option (0 no, 1 yes) ........ '',i10)')&
              rtemp,iff
            if (nsgz.gt.0) write(nsyso,'(&
              &'' max. no. of sigma zeroes ............. '',i10)')&
              nsgz
            if (nsgz.eq.0) write(nsyso,'(&
              &'' max. no. of sigma zeroes ............. '',7x,&
              &''all'')')
            write(nsyso,'(&
              &'' index of reference sigma zero ........ '',i10)')&
              izref
            write(nsyso,'(/&
              &'' description of nuclide''/&
              &1x,9(''----''),''--''/6x,4a4)') (word(i),i=1,4)

            !--set up pointers for gam.
            call gamll(igprnt,iff,nsgz)

            !--format cross sections for gam.
            call gamxs(locab0,locsf0,locchi,locnus,loce0,loce1,locin,&
              loc2n,locabs,a,na)

            !--format shielding factors for gam.
            if (nff.gt.0.and.iwr.gt.0) call gamff(ngmin,ngmax,&
              locab0,locsf0,locabs,locsf,a,na)

            !--write nscr and print output, if desired.
            xld1=xld(1)+1
            xld2=xld(2)+1
            xld3=xld(3)+1
            if (nfs+iwf.ne.0) then
               ifiss=iwf+nfs
               write(kscr,'(i6,i2,10a4)') nid,izero,(fsn(i),i=1,10)
               if (ifiss.eq.1) write(nsyso,'(/&
                 &'' fission spectrum''/&
                 &'' ----------------'')')
               if (ifiss.gt.1) write(nsyso,'(/&
                 &'' fission spectra''/&
                 &'' ---------------'')')
               write(nsyso,'(6x,10a4,''total''/)') (fsn(i),i=1,10)
               l=locchi
               lim=l+ngnd-1
               write(kscr,'(1p,6e12.5)') (a(i),i=l,lim)
               write(nsyso,'(1p,6e12.4)') (a(i),i=l,lim)
               if (nfs.ne.0) then
                  do ifs=1,nfs
                     write(kscr,'(i6,i2,10a4)')&
                       nid,ifs,(fsn(i),i=1,10)
                     write(nsyso,'(/6x,10a4,''delayed'',i1/)')&
                       (fsn(i),i=1,10),ifs
                     l=l+ngnd
                     lim=lim+ngnd
                     write(nsyso,'(1p,6e12.4)') (a(i),i=l,lim)
                     write(kscr,'(1p,6e12.5)') (a(i),i=l,lim)
                  enddo
                  nfs=nfs+iwf
               endif
            endif
         endif

         !--write more nout and print output
         write(nscr,'(6i10)') nid,iwa,iwf,iwr
         write(nscr,'(20a4)') (word(i),i=1,4)
         write(nscr,'(36x,1p,3e12.5)') xlol(3),xla(3),xld3
         write(nscr,'(1p,6e12.5)')&
           xlol(1),xla(1),xld1,xlol(2),xla(2),xld2
         if (iwa.ne.0) then
            lim=locab0+ngnd-1
            if (iread.eq.0) write(nscr,'(1p,6e12.5)')&
              (a(i),i=locab0,lim)
            if (iread.gt.0) write(nscr,'(1p,6e12.5)')&
              (z(i),i=1,ngnd)
            if (iprint.ne.0) then
               write(nsyso,'(/&
                 &'' absorption cross section''/&
                 &'' ------------------------'')')
               if (iread.eq.0) write(nsyso,'(1p,6e12.4)')&
                 (a(i),i=locab0,lim)
               if (iread.gt.0) write(nsyso,'(1p,6e12.4)')&
                 (z(i),i=1,ngnd)
            endif
         endif
         if (iwf.ne.0) then
            limsf=locsf0+ngnd-1
            write(nscr,'(1p,6e12.5)') (a(i),i=locsf0,limsf)
            lim=locnus+ngnd-1
            write(nscr,'(1p,6e12.5)') (a(i),i=locnus,lim)
            if (iprint.ne.0) then
               write(nsyso,'(/&
                 &'' fission cross section''/&
                 &'' ---------------------'')')
               write(nsyso,'(1p,6e12.4)') (a(i),i=locsf0,limsf)
               write(nsyso,'(/&
                 &'' fission nu''/&
                 &'' ----------'')')
               write(nsyso,'(1p,6e12.4)') (a(i),i=locnus,lim)
               if (nfs.ne.0) then
                  limc=locchi-1+ngnd
                  write(nsyso,'(/&
                    &'' fission chi''/&
                    &'' -----------'')')
                  write(nsyso,'(1p,6e12.4)') (a(i),i=locchi,limc)
               endif
            endif
         endif
         if (xlol(3).ne.zero) then
            lim=loce0+nint(xlol(3))-1
            write(nscr,'(1p,6e12.5)') (a(i),i=loce0,lim)
            lim=loce1+nint(xlol(3))-1
            write(nscr,'(1p,6e12.5)') (a(i),i=loce1,lim)
            if (iprint.ne.0) then
               la=nint(xla(3))
               ld=nint(xld(3))
               write(nsyso,'(//&
                 &'' p0 elastic scattering matrix''/&
                 &'' ----------------------------'')')
               call pgam(a(loce0),ngnd,la,ld)
               write(nsyso,'(//&
                 &'' p1 elastic scattering matrix''/&
                 &'' ----------------------------'')')
               call pgam(a(loce1),ngnd,la,ld)
            endif
         endif
         if (xlol(1).ne.zero) then
            lim=locin+nint(xlol(1))-1
            write(nscr,'(1p,6e12.5)') (a(i),i=locin,lim)
            if (iprint.ne.0) then
               la=nint(xla(1))
               ld=nint(xld(1))
               write(nsyso,'(//&
                 &'' inelastic scattering matrix''/&
                 &'' ---------------------------'')')
               call pgam(a(locin),ngnd,la,ld)
            endif
         endif
         if (xlol(2).ne.zero) then
            lim=loc2n+nint(xlol(2))-1
            write(nscr,'(1p,6e12.5)') (a(i),i=loc2n,lim)
            if (iprint.ne.0) then
               la=nint(xla(2))
               ld=nint(xld(2))
               write(nsyso,'(//&
                 &'' n2n scattering matrix''/&
                 &'' ---------------------'')')
               call pgam(a(loc2n),ngnd,la,ld)
            endif
         endif
         if (nff.ne.0.and.iwr.ne.0) then
            write(nscr,'(4i6)') nsigz,ntp,misc,jwf
            write(nscr,'(1p,6e12.5)') (sigz(i),i=1,nsigz)
            write(nscr,'(1p,6e12.5)') (tmpr(i),i=1,ntp)
            istart=locabs+(ngmax-iwr)*nsigz*ntp
            lim=istart+iwr*nsigz*ntp-1
            write(nscr,'(1p,6e12.5)') (a(i),i=istart,lim)
            if (iprint.ne.0) then
               write(nsyso,'(//&
                 &'' absorption self-shielding factors''/&
                 &'' ---------------------------------'')')
               call pgamff(a(locabs),ntp,nsigz,iwr,sigz,tmpr,ngmax)
            endif
            if (jwf.ne.0) then
               istart=locsf+(ngmax-iwr)*nsigz*ntp
               lim=istart+iwr*nsigz*ntp-1
               write(nscr,'(1p,6e12.5)') (a(i),i=istart,lim)
               if (iprint.ne.0) then
                  write(nsyso,'(//&
                    &'' fission self-shielding factors''/&
                    &'' ------------------------------'')')
                  call pgamff(a(locsf),ntp,nsigz,iwr,sigz,tmpr,ngmax)
               endif
            endif
            ! write a dummy amisc array
            nwmisc=iwr*10
            write(nscr,'(1p,6e12.5)') (zero,i=1,nwmisc)
         endif
      endif
   enddo

   !--fast is finished.
   call repoz(nscr)
   call mergt(ifiss,kscr,nscr,nout,z)
   call closz(nscr)
   call closz(nscr)
   return
   end subroutine fast

   subroutine pgam(a,ng,la,ld)
   !--------------------------------------------------------------------
   ! Print out scattering matrices in packed GAM format.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   ! externals
   integer::ng,la,ld
   real(kr)::a(*)
   ! internals
   integer::ncol1,i,k,ig,n,i1,i2,j,nw,m,ig2
   integer::ncol=5

   ncol1=ncol-1
   write(nsyso,'(/&
     &'' initl final    cross section vs final group''/&
     &'' group group    +0'',4(10x,''+'',i1)/&
     &1x,7(''----------''),''--------'')')&
     (i,i=1,ncol1)
   k=0
   do ig=1,la
      n=ng-ig+1
      if (n.gt.(ld+1)) n=ld+1
      if (n.le.ld) n=n+1
      ! do not print leading and trailing zeroes
      i1=0
      i2=0
      do j=1,n
         if (i1.eq.0.and.a(k+j).ne.0) i1=j
         if (a(k+j).ne.0.) i2=j
      enddo
      if (i1.ne.0) then
         nw=i2-i1+1
         m=nw
         k=k+i1-1
         if (m.gt.ncol) m=ncol
         ig2=ig+i1-1
         write(nsyso,'(i5,i6,2x,1p,5e12.4)') ig,ig2,(a(k+i),i=1,m)
         do while (nw.gt.0)
            k=k+m
            nw=nw-m
            if (nw.gt.0) then
               ig2=ig2+ncol
               m=nw
               if (m.gt.ncol) m=ncol
               write(nsyso,'(5x,i6,2x,1p,5e12.4)') ig2,(a(k+i),i=1,m)
            endif
         enddo
      endif
      k=k+n-i2
   enddo
   return
   end subroutine pgam

   subroutine pgamff(a,ntp,nsigz,iwr,sigz,tmpr,ngmax)
   !--------------------------------------------------------------------
   ! Print out form factors in packed GAM format.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   ! externals
   integer::ntp,nsigz,iwr,ngmax
   real(kr)::sigz(*),tmpr(*),a(nsigz,ntp,*),f(7)
   ! internals
   integer::i,n,ig,is,it
   integer::jopt=1

   write(nsyso,'(5x,''temperatures (down)'',10f7.0)')&
     (tmpr(i),i=1,ntp)
   write(nsyso,'(5x,''sigzeroes (across)''/9x,1p,10e12.0)')&
     (sigz(i),i=1,nsigz)
   n=ngmax+1-iwr
   do ig=n,ngmax
      if (jopt.le.0) then
         write(nsyso,'(/'' group'',i2,3x,1p,10e12.4)')&
           ig,(a(is,1,ig),is=1,nsigz)
      else
         do is=1,nsigz
            f(is)=exp(-a(is,1,ig)/2)*sqrt(sigz(is))
         enddo
         write(nsyso,'(/'' group'',i2,3x,1p,10e12.4)')&
           ig,(f(is),is=1,nsigz)
      endif
      if (ntp.ne.1) then
         do it=2,ntp
            if (jopt.le.0) then
               write(nsyso,'(11x,1p,10e12.4)')&
                 (a(is,it,ig),is=1,nsigz)
            else
               do is=1,nsigz
                  f(is)=exp(-a(is,it,ig)/2)*sqrt(sigz(is))
               enddo
               write(nsyso,'(11x,1p,10e12.4)')&
                 (f(is),is=1,nsigz)
            endif
         enddo
      endif
   enddo
   return
   end subroutine pgamff

   subroutine mergt(ifiss,n1,n2,n3,a)
   !--------------------------------------------------------------------
   ! copy tape n1 if present and n2 to n3
   !--------------------------------------------------------------------
   use util ! provides repoz
   ! externals
   integer::ifiss,n1,n2,n3
   real(kr)::a(*)
   ! internals
   integer::i

   call repoz(n2)
   call repoz(n3)
   if (ifiss.gt.0) call repoz(n1)
   if (ifiss.eq.0) go to 130
  110 continue
   read(n1,10,end=130) (a(i),i=1,20)
   write(n3,10) (a(i),i=1,20)
   go to 110
  130 continue
   read(n2,10,end=150) (a(i),i=1,20)
   write(n3,10) (a(i),i=1,20)
   go to 130
  150 continue
   call repoz(n3)
   return
  10 format(20a4)
   end subroutine mergt

   subroutine gamll(igprnt,iff,nsgz)
   !--------------------------------------------------------------------
   ! Read through gendf tape and assign storage and pointers
   ! for GAM format.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use endf ! provides endf routines and variables
   use util ! provides repoz,error
   ! externals
   integer::igprnt,iff,nsgz
   ! internals
   integer::nb,nw,ldin,ldp,ld2n,lz,iza,ntw,i1,i,ngnd1,ngn1
   integer::lastg,ig,jg,lain,la2n,lap,nl,l,igc,ig2lo,idone
   integer::ig2loc,lx
   real(kr)::temp
   character(60)::strng
   real(kr),dimension(:),allocatable::eg
   real(kr),parameter::eps=1.e-4_kr
   integer,parameter::ntpmax=5
   integer,parameter::nszmax=7

   !--read data for desired material.
   call repoz(ngendf)
   call tpidio(ngendf,0,0,scr,nb,nw)
  110 continue
   call contio(ngendf,0,0,scr,nb,nw)
   if (math.eq.-1) then
      write(strng,'(''mat='',i4,'' not on tape '',i3)') matd,ngendf
      call error('gamll',strng,' ')
   endif
   if (math.eq.0.or.mfh.eq.0.or.mth.eq.0) go to 110
   if (math.eq.matd) go to 120
   call tomend(ngendf,0,0,scr)
   go to 110
  120 continue
   iwa=1
   iwf=0
   nfs=0
   iwr=0
   nsp=0
   ntp=-1
   jwf=0
   ldin=0
   ldp=0
   ld2n=0
   lz=6
   iza=nint(c1h)
   nid=iza*10
   go to 140
  130 continue
   call contio(ngendf,0,0,scr,nb,nw)
   if (math.ne.matd) go to 400
  140 continue
   nsigz=l2h
   if (nsigz.eq.1) iff=0
   if (nsigz.gt.1.and.iff.gt.0) iwr=1
   if (nsigz.gt.nsgz.and.nsgz.gt.0) nsigz=nsgz
   if (nsigz.gt.nszmax) then
      write(strng,'(i2,'' allowed'')') nszmax
      call error('gamll','reqested too many sigma zeroes',strng)
   endif
   ntw=n2h
   call listio(ngendf,0,0,scr,nb,nw)
   temp=c1h
   if (abs(temp-rtemp).le.eps) go to 160
   if (temp.gt.rtemp.and.ntp.lt.0)&
     call error('gamll','reference temp not on gendf tape.',' ')
   if (ntp.lt.0) go to 150
   ntp=ntp+1
   if (ntp.gt.ntpmax) then
      write(strng,'(i2,'' allowed'')') ntpmax
      call error('gamll','requested too many temperatures',strng)
   endif
   tmpr(ntp)=temp
   go to 220
  150 continue
   call tomend(ngendf,0,0,scr)
   go to 130
  160 continue
   ngn=l1h
   if (ngn.ne.ngnd.and.ngn.ne.ngnf)&
     call error('gamll',&
     'incorrect group structure or group type.',' ')
   if (igprnt.le.0) then
      i1=7+ntw+nsigz+ngn
      if (ngn.eq.ngnd) write(nsyso,'(/&
        &'' neutron group structure''/&
        &'' -----------------------'',(1p,6e12.4))')&
        scr(i1),(scr(i1-i),i=1,ngn)
      if (ngn.ne.ngnd) then
         ngnd1=ngnd+1
         ngn1=ngn+1
         allocate(eg(ngnd1))
         lastg=0
         i1=i1-ngn1
         do i=1,ngn1
            ig=icgrp(i)
            if (ig.ne.lastg) then
               if (ig.ne.0) then
                  jg=ngnd1-ig
                  eg(1+jg)=scr(i1+i)
               endif
               lastg=ig
            endif
         enddo
         write(nsyso,'(/&
           &'' neutron group structure''/&
           &'' -----------------------''/&
           &(1p,6e12.4))') (eg(i),i=1,ngnd1)
         deallocate(eg)
      endif
      igprnt=1
   endif
   ntp=0
   if (rtemp.gt.0) ntp=1
   if (rtemp.gt.0) tmpr(1)=temp
   do i=1,nsigz
      sigz(i)=scr(i+ntw+lz)
   enddo
   if (izref.eq.1) write(nsyso,'(/&
     &'' ref. sigma zero for elastic matrix = infinity'')')
   if (izref.ne.1) write(nsyso,'(/&
     &'' ref. sigma zero for elastic matrix ='',1p,e10.2)')&
     sigz(izref)
   lain=ngnd
   la2n=ngnd
   lap=ngnd
  220 continue
   call tosend(ngendf,0,0,scr)

   !--process relevant sections for this material.
  250 continue
   call contio(ngendf,0,0,scr,nb,nw)
   if (math.eq.0.and.nsigz.gt.1) go to 130
   if (math.ne.matd) go to 400
   if (mfh.eq.0.or.mth.eq.0) go to 250
   if (mfh.eq.3) go to 260
   if (mfh.eq.5) go to 270
   if (mfh.eq.6) go to 300
   call tosend(ngendf,0,0,scr)
   go to 250

   !--process cross section data (mf=3).
   !--check for presence of fission.
  260 continue
   if (mth.eq.18.or.mth.eq.19) iwf=1
   call tosend(ngendf,0,0,scr)
   go to 250

   !--check for presence of delayed neutron spectra (mf=5)
  270 continue
   if (mth.eq.455) then
      nl=l1h
      nfs=nfs+nl
   endif
   call tosend(ngendf,0,0,scr)
   go to 250

   !--process transfer matrix data (mf=6).
  300 continue
   call listio(ngendf,0,0,scr,nb,nw)
   l=1
   do while (nb.ne.0)
      l=l+nw
      call moreio(ngendf,0,0,scr(l),nb,nw)
   enddo
   ig=n2h
   igc=icgrp(ig)
   if (igc.eq.0.or.igc.gt.ngnd) go to 350
   ig2lo=l2h
   idone=0
   do while (idone.eq.0)
      ig2loc=icgrp(ig2lo)
      if (ig2loc.ne.0.and.ig2loc.le.ngnd) then
         idone=1
      else
         if (ig2lo.eq.ig) then
            idone=1
         else
            ig2lo=ig2lo+1
         endif
      endif
   enddo
   if (mth.lt.51.or.mth.gt.91) go to 360
   ! inelastic
   if (lain.gt.igc) lain=igc
   if (ldin.lt.(igc-ig2loc)) ldin=igc-ig2loc
  350 continue
   if (ig.lt.ngn) go to 300
   go to 250
  360 continue
   if (mth.ne.2) go to 370
   ! elastic
   if (lap.gt.igc) lap=igc
   if (ldp.lt.(igc-ig2loc)) ldp=igc-ig2loc
   if (ig.lt.ngn) go to 300
   go to 250
  370 continue
   if (mth.eq.16) go to 390
   if (mth.gt.5.and.mth.lt.10) go to 390
   if (mth.gt.45.and.mth.lt.50) go to 390
   call tosend(ngendf,0,0,scr)
   go to 250
  390 continue
   ! n2n
   if (la2n.gt.igc) la2n=igc
   if (ld2n.lt.(igc-ig2loc)) ld2n=igc-ig2loc
   go to 250

   !--set up storage pointers.
  400 continue
   xla(1)=0
   xla(2)=0
   xla(3)=0
   xlol(1)=0
   xlol(2)=0
   xlol(3)=0
   xld(1)=ldin
   xld(2)=ld2n
   xld(3)=ldp
   if (ldin.ne.0) then
      xla(1)=ngnd-lain+1
      lx=nint(xla(1))+nint(xld(1))-ngnd-1
      if (lx.lt.0) lx=0
      xlol(1)=xla(1)*(xld(1)+1)-(lx*(lx-1))/2
   endif
   if (ld2n.ne.0) then
      xla(2)=ngnd-la2n+1
      lx=nint(xla(2))+nint(xld(2))-ngnd-1
      if (lx.lt.0) lx=0
      xlol(2)=xla(2)*(xld(2)+1)-(lx*(lx-1))/2
   endif
   if (ldp.ne.0) then
      xla(3)=ngnd-lap+1
      lx=nint(xla(3))+nint(xld(3))-ngnd-1
      if (lx.lt.0) lx=0
      xlol(3)=xla(3)*(xld(3)+1)-(lx*(lx-1))/2
   endif
   return
   end subroutine gamll

   subroutine gamxs(locab0,locsf0,locchi,locnus,loce0,loce1,locin,&
     loc2n,locabs,a,nloc)
   !--------------------------------------------------------------------
   ! Read through ngendf tape again and accumulate desired
   ! cross sections and matrices in the array a.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use endf ! provides endf routines and variables
   use util ! provides error,repoz,skiprz
   ! externals
   integer::locab0,locsf0,locchi,locnus,loce0,loce1,locin,loc2n,locabs
   real(kr)::a(*)
   integer::nloc
   ! externals
   integer::ngnd1,ldin,kdin,ld2n,kd2n,ldp,kdp,lz,mf3mt2,jgdnu
   integer::i,i318,i618,nb,nw,l,nl,nz,ng2,ig2lo,ig,igc
   integer::jg,loc,loca,k,locf,ig2,ig2c,jg2c,locc,nldla,il,na
   integer::locb,jgc,locd,jgt,jgp,lc,lim,loc1,loc2,iloc,locdla
   real(kr)::dnu,dnorm,cnorm,temp,test,flx,rnorm,fracd,tfound
   real(kr)::sumd(6)
   real(kr)::cspc(200)
   real(kr),parameter::eps=1.e-4_kr
   real(kr),parameter::zero=0

   !--initialize data pointers.
   ngnd1=ngnd+1
   ldin=nint(xld(1))
   kdin=ngnd1-ldin
   ld2n=nint(xld(2))
   kd2n=ngnd1-ld2n
   ldp=nint(xld(3))
   kdp=ngnd1-ldp
   lz=6
   mf3mt2=0
   dnu=0
   jgdnu=0
   dnorm=0
   cnorm=0
   iloc=1
   locab0=iloc
   locsf0=locab0+ngnd
   locchi=locsf0+ngnd
   locdla=locchi+ngnd
   locnus=locdla+ngnd*nfs
   loce0=locnus+ngnd
   loce1=loce0+nint(xlol(3))
   locin=loce1+nint(xlol(3))
   loc2n=locin+nint(xlol(1))
   locabs=loc2n+nint(xlol(2))
   na=locabs-1
   if ((na-iloc).gt.nloc)&
     call error('gamxs','storage exceeded.',' ')
   do i=iloc,na
      a(i)=0
   enddo
   do i=1,ngnd
      cflux(i)=0
      cflux(i+ngnd)=0
   enddo
   do i=1,ngn
      cspc(i)=0
   enddo
   i318=0
   i618=0
   call repoz(ngendf)
   call tpidio(ngendf,0,0,scr,nb,nw)
  120 continue
   call contio(ngendf,0,0,scr,nb,nw)
   if (math.eq.matd) go to 130
   call tomend(ngendf,0,0,scr)
   go to 120

   !--process relevant sections for reference temperature.
  130 continue
   tfound=0
   if (mfh.ne.1) go to 150
   call tosend(ngendf,0,0,scr)
  140 continue
   l=1
   call contio(ngendf,0,0,scr(l),nb,nw)
   if (math.eq.0) go to 500
   if (mfh.eq.0.or.mth.eq.0) go to 140
   if (mfh.ne.1) go to 150
   call tosend(ngendf,0,0,scr(l))
   go to 140
  150 continue
   nl=l1h
   nz=l2h
   if (mfh.ne.3.or.mth.lt.2.or.mf3mt2.gt.0) go to 170
   mf3mt2=1
   do i=1,ngnd
      if (cflux(i).ne.zero) cflux(i)=1/cflux(i)
      if (cflux(ngnd+i).ne.zero) cflux(ngnd+i)=1/cflux(ngnd+i)
   enddo
  170 continue
   call listio(ngendf,0,0,scr(l),nb,nw)
   temp=c1h
   if (abs(temp-rtemp).gt.eps.and.tfound.eq.0) go to 210
   if (abs(temp-rtemp).le.eps) tfound=1
   ng2=l1h
   ig2lo=l2h
   ig=n2h
   igc=1
   if (ig.gt.0) igc=icgrp(ig)
   do while (nb.ne.0)
      l=l+nw
      call moreio(ngendf,0,0,scr(l),nb,nw)
      if ((l+nw-1).gt.nwscr)&
        call error('gamxs','storage exceeded.',' ')
   enddo
   if (igc.eq.0.or.igc.gt.ngnd) go to 260
   jg=ngnd-igc+1
   l=1
   if (mth.eq.455) go to 330
   if (mfh.eq.3) go to 250
   if (mfh.eq.6) go to 300
  200 continue
   call tosend(ngendf,0,0,scr(l))
   go to 140
  210 continue
   call tomend(ngendf,0,0,scr(l))
   go to 140

   !--file 3 cross sections
  250 continue
   if (mth.ne.1) go to 270
   ! coarse group p0 flux
   cflux(jg)=cflux(jg)+scr(l+lz)
   if (nl.eq.1) go to 260
   ! coarse group p1 flux
   cflux(ngnd+jg)=cflux(ngnd+jg)+scr(l+lz+1)
  260 continue
   if (ig.lt.ngn) go to 170
   go to 140
  270 continue
   if (mth.lt.102.or.mth.gt.150) go to 280
   ! capture
   loc=locab0+jg-1
   loca=l+lz+nl*nz
   a(loc)=a(loc)+scr(loca)*scr(l+lz)*cflux(jg)
   if (ig.lt.ngn) go to 170
   go to 140
  280 continue
   if ((mth.lt.18.or.mth.gt.21).and.mth.ne.38) go to 200
   ! fission
   if (mth.eq.18) i318=1
   if (mth.gt.18.and.i318.gt.0) go to 200
   loc=locsf0+jg-1
   loca=l+lz+nl*nz
   a(loc)=a(loc)+scr(loca)*scr(l+lz)*cflux(jg)
   loc=locab0+jg-1
   a(loc)=a(loc)+scr(loca)*scr(l+lz)*cflux(jg)
   if (ig.lt.ngn) go to 170
   go to 140

   !--file 6 matrices
  300 continue
   if (igc.eq.0.or.igc.gt.ngnd) go to 320
   if ((mth.lt.18.or.mth.gt.21).and.mth.ne.38) go to 370

   !--fission
   if (mth.eq.18) i618=1
   if (mth.gt.18.and.i618.gt.0) go to 200
   if (ig.ne.0) then
      do k=2,ng2
         loca=l+lz+nl*nz*(k-1)
         loc=locnus+jg-1
         if (ig2lo.ne.0) then
            ! matrix part
            a(loc)=a(loc)+scr(loca)*scr(l+lz)*cflux(jg)
            locf=l+lz
            ig2=ig2lo+k-2
            ig2c=icgrp(ig2)
            if (ig2c.ne.0.and.ig2c.le.ngnd) then
               jg2c=ngnd-ig2c+1
               loc=locchi+jg2c-1
               a(loc)=a(loc)+scr(locf)*scr(loca)
               cnorm=cnorm+scr(locf)*scr(loca)
            endif
         else
            ! spectrum part
            do i=1,ngn
               a(loc)=a(loc)+cspc(i)*scr(loca)*scr(loca-1)*cflux(jg)
               ig2c=icgrp(i)
               if (ig2c.ne.0.and.ig2c.le.ngnd) then
                  jg2c=ngnd-ig2c-1
                  locc=locchi+jg2c-1
                  a(locc)=a(locc)+cspc(i)*scr(loca)*scr(loca-1)
                  cnorm=cnorm+scr(loca)*scr(loca-1)
               endif
            enddo
         endif
      enddo
   else
      ! save constant spectrum
      do k=1,ng2
         cspc(ig2lo+k-1)=scr(l+lz+nl*nz*(k-1))
      enddo
   endif
  320 continue
   if (ig.lt.ngn) go to 170
   go to 140
   ! delayed neutrons
  330 continue
   if (mfh.ne.5) then
      loc=locnus+jg-1
      loca=l+lz
      a(loc)=a(loc)+scr(loca+1)*scr(loca+2)*scr(l+lz)*cflux(jg)
      dnorm=dnorm+scr(loca)*scr(loca+1)*scr(loca+2)
      test=1
      test=test/10
      if (scr(loca+2).ge.test.and.dnu.eq.zero) jgdnu=jg
      if (scr(loca+2).ge.test.and.dnu.eq.zero) dnu=scr(loca+1)
   else
      nldla=nl
      if (nldla.ne.6) call error('gamxs',&
        'no. deleay groups not 6.',' ')
      do il=1,nl
         sumd(il)=0
      enddo
      locb=l+lz-1
      do k=2,ng2
         ig2=ig2lo+k-2
         ig2c=icgrp(ig2)
         if (ig2c.ne.0.and.ig2c.le.ngnd) then
            jgc=ngnd-ig2c+1
            loc=locchi+jgc-1
            loca=l+lz-1+nl*(k-1)
            locd=locdla+jgc-1
            do jgt=1,nl
               a(loc)=scr(jgt+loca)*dnorm+a(loc)
               a(ngnd*(jgt-1)+locd)=a(ngnd*(jgt-1)+locd)&
                 +scr(jgt+loca)*dnu
               if (jgc.eq.ngnd) a(ngnd*(jgt-1)+locd)=a(jgt+locb)
               if (jgc.lt.ngnd-1) sumd(jgt)=sumd(jgt)&
                 +scr(jgt+loca)*dnu
               cnorm=scr(jgt+loca)*dnorm+cnorm
            enddo
         endif
      enddo
   endif
   if (ig.lt.ngn) go to 170
   go to 140

   !--inelastic
  370 continue
   if (mth.lt.51.or.mth.gt.91) go to 410
   jgp=jg
   if (jg.gt.kdin) jgp=kdin
   lc=(ldin+1)*(jgp-1)
   if (jg.gt.kdin) then
      lim=jg-kdin
      do jgt=1,lim
         lc=lc+ldin+2-jgt
      enddo
   endif
   do k=2,ng2
      ig2=ig2lo+k-2
      ig2c=icgrp(ig2)
      if (ig2c.ne.0.and.ig2c.le.ngnd) then
         jg2c=ngnd-ig2c+1
         if (jg2c.lt.jg) jg2c=jg
         loca=l+lz+nl*nz*(k-1)
         loc=locin+lc+jg2c-jg
         a(loc)=a(loc)+scr(loca)*scr(l+lz)*cflux(jg)
      endif
   enddo
   if (ig.lt.ngn) go to 170
   go to 140

   !--elastic
  410 continue
   if (mth.ne.2) go to 450
   jgp=jg
   if (jg.gt.kdp) jgp=kdp
   lc=(ldp+1)*(jgp-1)
   if (jg.gt.kdp) then
      lim=jg-kdp
      do jgt=1,lim
         lc=lc+ldp+2-jgt
      enddo
   endif
   do k=2,ng2
      ig2=ig2lo+k-2
      ig2c=icgrp(ig2)
      if (ig2c.ne.0.and.ig2c.le.ngnd) then
         jg2c=ngnd-ig2c+1
         if (jg2c.lt.jg) jg2c=jg
         loc=loce0+lc+jg2c-jg
         loca=l+lz+nl*(izref-1+nz*(k-1))
         a(loc)=a(loc)+scr(loca)*scr(l+lz)*cflux(jg)
         if (nl.ne.1) then
            loca=loca+1
            loc=loce1+lc+jg2c-jg
            flx=cflux(ngnd+jg)
            if (cflux(ngnd+jg).eq.0.) flx=cflux(jg)
            a(loc)=a(loc)+3*scr(loca)*scr(l+lz+1)*flx
         endif
      endif
   enddo
   if (ig.lt.ngn) go to 170
   go to 140
  450 continue
   if (mth.eq.16) go to 460
   if (mth.lt.6.or.mth.gt.49) go to 200
   if (mth.gt.9.and.mth.lt.46) go to 200

   !--n2n
  460 continue
   jgp=jg
   if (jg.gt.kd2n) jgp=kd2n
   lc=(ld2n+1)*(jgp-1)
   if (jg.gt.kd2n) then
      lim=jg-kd2n
      do jgt=1,lim
         lc=lc+ld2n+2-jgt
      enddo
   endif
   do k=2,ng2
      ig2=ig2lo+k-2
      ig2c=icgrp(ig2)
      if (ig2c.ne.0.and.ig2c.le.ngnd) then
         jg2c=ngnd-ig2c+1
         if (jg2c.lt.jg) jg2c=jg
         loca=l+lz+nl*nz*(k-1)
         loc=loc2n+lc+jg2c-jg
         a(loc)=a(loc)+scr(loca)*scr(l+lz)*cflux(jg)/2
      endif
   enddo
   if (ig.lt.ngn) go to 170
   go to 140

   !--finished loading a.
   !--find nu from nu times sigma f.
  500 continue
   do ig=1,ngnd
      loc1=locsf0+ig-1
      if (a(loc1).ne.zero) then
         loc2=locnus+ig-1
         a(loc2)=a(loc2)/a(loc1)
      endif
   enddo

   !--normalize chi and delayed chi
   if (cnorm.ne.zero) then
      cnorm=1/cnorm
      do ig=1,ngnd
         a(ig-1+locchi)=a(ig-1+locchi)*cnorm
      enddo
      if (dnu.ne.zero) then
         rnorm=1/a(jgdnu-1+locnus)
         lim=ngnd-1
         do ig=1,lim
            locd=locdla-1+ig
            do il=1,nldla
               if (ig.lt.lim) then
                  a(ngnd*(il-1)+locd)=a(ngnd*(il-1)+locd)*rnorm
               else
                  a(ngnd*(il-1)+locd)=sumd(il)*rnorm
               endif
            enddo
         enddo
         fracd=dnu/a(jgdnu-1+locnus)
         write(nsyso,'(/&
           &'' nudlay ='',1p,e13.5/&
           &'' delayed fraction ='',e13.5/&
           &'' time constants''/1x,6e13.5/&
           &'' time group fractions''/1x,6e13.5)')&
           dnu,fracd,(a(ngnd*(il-1)+locdla+ngnd-1),il=1,nldla),&
           (a(ngnd*(il-1)+locdla+ngnd-2),il=1,nldla)
      endif
   endif

   !--for reference temperatures greater than zero...
   if (rtemp.gt.0) call skiprz(ngendf,-2)
   if (rtemp.gt.0) call findf(matd,1,451,ngendf)
   return
   end subroutine gamxs

   subroutine gamff(ngmin,ngmax,locab0,locsf0,locabs,locsf,a,nloc)
   !--------------------------------------------------------------------
   ! Calculate f-factors for absorption and fission.
   !--------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use util ! provides error
   ! externals
   integer::ngmin,ngmax,locab0,locsf0,locabs,locsf,nloc
   real(kr)::a(*)
   ! internals
   integer::lz,i,itp,mf3mt1,i318,l,nb,nw,nl,nz,lim,ig,igc
   integer::jg,loc,loca,jz,iadd,iglo,na,iloc
   real(kr)::ab0,sf0
   real(kr),parameter::toler=.01e0_kr
   real(kr),parameter::zero=0

   !--initialize.
   iloc=1
   lz=6
   locsf=locabs+ngnd*ntp*nsigz
   na=locsf+ngnd*ntp*nsigz-1
   if ((na-iloc).gt.nloc)&
     call error('gamff','storage exceeded.',' ')
   do i=locabs,na
      a(i)=0
   enddo
   itp=1
   mf3mt1=1
   i318=0

   !--read through rest of temperatures for this mat,
   !--and store absorption and fission cross-sections.
  110 continue
   l=1
   call contio(ngendf,0,0,scr(l),nb,nw)
   if (math.gt.0) go to 120
   if (itp.eq.ntp) go to 300
   itp=itp+1
   i318=0
  120 continue
   if (mfh.eq.0.or.mth.eq.0) go to 110
   if (mfh.eq.3) go to 130
   call tosend(ngendf,0,0,scr(l))
   go to 110
  130 continue
   if (mf3mt1.le.0.and.mfh.eq.3.and.mth.ge.2) then
      mf3mt1=1
      do i=1,ngnd
         if (cflux(i).ne.zero) cflux(i)=1/cflux(i)
      enddo
   endif
   if (mth.gt.101.and.mth.lt.108) go to 160
   if ((mth.gt.17.and.mth.lt.22).or.mth.eq.38) go to 160
   if (mth.eq.1) go to 160
  150 continue
   call tosend(ngendf,0,0,scr(l))
   go to 110
  160 continue
   nl=l1h
   nz=l2h
   lim=nz
   if (nsigz.lt.nz) lim=nsigz
  170 continue
   l=1
   call listio(ngendf,0,0,scr(l),nb,nw)
   ig=n2h
   igc=icgrp(ig)
   do while (nb.ne.0)
      l=l+nw
      if ((l+nw-1).gt.nwscr)&
        call error('gamff','storage exceeded.',' ')
      call moreio(ngendf,0,0,scr(l),nb,nw)
   enddo
   if (igc.eq.0.or.igc.gt.ngnd) go to 210
   jg=ngnd-igc+1
   if (mth.lt.102) go to 220
   ! capture
   loc=locabs-1+nsigz*((itp-1)+ntp*(jg-1))
   loca=l+lz+nl*(nz-1)
   do jz=1,lim
      a(jz+loc)=scr(nl*jz+loca)*scr(l+lz)*cflux(jg)+a(jz+loc)
   enddo
  210 continue
   if (ig.lt.ngn) go to 170
   go to 110
  220 continue
   if (mth.eq.1) go to 240
   ! fission
   if (mth.eq.18) i318=1
   if (mth.gt.18.and.i318.gt.0) go to 150
   jwf=1
   iadd=nsigz*((itp-1)+ntp*(jg-1))-1
   loca=l+lz+nl*(nz-1)
   do jz=1,lim
      a(jz+iadd+locsf)=a(jz+iadd+locsf)+scr(nl*jz+loca)*scr(l+lz)*&
        cflux(jg)
      a(jz+iadd+locabs)=a(jz+iadd+locabs)+scr(nl*jz+loca)*scr(l+lz)*&
        cflux(jg)
   enddo
   go to 210
   ! flux
  240 continue
   if (igc.le.1) then
      do i=1,ngnd
         cflux(i)=0
      enddo
      mf3mt1=0
   endif
   cflux(jg)=cflux(jg)+scr(l+lz)
   go to 210

   !--finished loading a.
   !--calculate f-factors.
  300 continue
   iglo=ngmin
   do ig=ngmin,ngmax
      ab0=a(locab0+ig-1)
      if (ab0.ne.zero) ab0=1/ab0
      sf0=a(locsf0+ig-1)
      if (sf0.ne.zero) sf0=1/sf0
      do jz=1,nsigz
         do itp=1,ntp
            iadd=jz-1+nsigz*((itp-1)+ntp*(ig-1))
            loc=locabs+iadd
            if (ab0.eq.zero) a(loc)=1
            if (ab0.ne.zero) a(loc)=a(loc)*ab0
            if (a(loc).le.zero) a(loc)=1
            if (abs(1-a(loc)).gt.toler.and.iglo.eq.ngmin) iglo=ig
            a(loc)=log(sigz(jz))-2*log(a(loc))
            if (jwf.ne.0) then
               loc=locsf+iadd
               if (sf0.eq.zero) a(loc)=1
               if (sf0.ne.zero) a(loc)=a(loc)*sf0
               if (a(loc).le.zero) a(loc)=1
               a(loc)=log(sigz(jz))-2*log(a(loc))
            endif
         enddo
      enddo
   enddo
   iwr=ngmax-iglo+1
   return
   end subroutine gamff

   subroutine therm(iprint,igprnt)
   !--------------------------------------------------------------------
   ! Read data for adding nuclides to thermal library from a gout tape.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use endf ! provides endf routines and variables
   use util ! provides repoz,error
   ! externals
   integer::iprint,igprnt
   ! internals
   integer::il,lz,n,matd,itrc,mti,mtc,i,i318,mf3mt2
   integer::iflx0,iflx1,nb,nw,ifound,nl,nz,iza,nuclid,nsigz
   integer::ntw,itemp,idm,idma,idb,ida,idf,id1,id2,id3,nlm
   integer::ia,lim,ig,igc,ig2lo,ng2,loca,ig2,ig2c,loc
   integer::ndexc,ndexf,ngndsq,idtemp,ngnd1,nbm
   real(kr)::awr,amubar
   character(60)::strng
   character(12)::idname
   real(kr)::valu(20),aname(3),z(10)
   integer::igi(35),ngi(35)
      data valu/20*0.0d0/
   real(kr),dimension(:),allocatable::s0,s1
   real(kr),dimension(:),allocatable::xa,xf
   real(kr),dimension(:),allocatable::b,p
   real(kr),dimension(:),allocatable::grp
   real(kr),parameter::watt=1.e13_kr
   integer::id4=0
   integer::id5=0
   real(kr),parameter::zero=0

   ngnd=35
   write(nsyso,'(/'' ***epri-cell thermal library***/'')')

   !--allocate storage.
   allocate(s0(ngnd))
   allocate(s1(ngnd))
   allocate(xa(ngnd))
   allocate(xf(ngnd))
   ngndsq=ngnd*ngnd
   allocate(b(ngndsq))
   allocate(p(ngndsq))
   ngnd1=ngnd+1
   allocate(grp(ngnd1))
   il=1
   lz=6

   !--read and write user input.
   n=0
  100 continue
   idtemp=300
   idname=' '
   read(nsysi,*) matd,idtemp,idname
   if (matd.eq.0) go to 300
   read(idname,'(3a4)') aname(1),aname(2),aname(3)
   itrc=0
   mti=0
   mtc=0
   read(nsysi,*) itrc,mti,mtc
   do i=1,8
      z(i)=0
   enddo
   read(nsysi,*) (z(i),i=1,8)
   valu(2)=z(1)
   valu(3)=z(2)
   valu(4)=z(3)
   valu(5)=z(4)
   valu(6)=z(5)*watt
   valu(7)=z(6)*watt
   valu(8)=0
   valu(9)=z(7)
   valu(10)=z(8)
   n=n+1
   write(nsyso,'(//&
     &'' mat to be processed .................. '',i10/&
     &'' temperature id ....................... '',i10/&
     &'' isotope name ......................... '',2a4,a2)')&
     matd,idtemp,aname
   write(nsyso,'(&
     &'' transport corr. option (0 no, 1 yes) . '',i10/&
     &'' mt inelastic ......................... '',i10/&
     &'' mt elastic ........................... '',i10)')&
     itrc,mti,mtc
   write(nsyso,'(&
     &'' xi ................................... '',1p,e10.3/&
     &'' alpha ................................ '',e10.3/&
     &'' mubar ................................ '',e10.3/&
     &'' nu ................................... '',e10.3/&
     &'' kappa fission ........................ '',e10.3/&
     &'' kappa capture ........................ '',e10.3/&
     &'' lambda ............................... '',e10.3/&
     &'' sigma s .............................. '',e10.3)')&
     (z(i),i=1,8)

   !--read data for desired nuclide.
   i318=0
   mf3mt2=0
   do i=1,ngndsq
      b(i)=0
      p(i)=0
   enddo
   iflx0=0
   iflx1=0
   call repoz(ngendf)
   call contio(ngendf,0,0,scr,nb,nw)
   ifound=0
  115 continue
   call contio(ngendf,0,0,scr,nb,nw)
   if (math.lt.0) go to 220
   if (math.eq.0.or.mfh.eq.0.or.mth.eq.0) go to 115
   if (math.eq.matd) go to 120
   if (ifound.gt.0) go to 220
   call tomend(ngendf,0,0,scr)
   go to 115
  120 continue
   nl=l1h
      nz=l2h
   if (mfh.eq.3) go to 152
   if (mfh.eq.6.and.mth.eq.mti) go to 152
   if (mfh.eq.6.and.mth.eq.mtc.and.mtc.ne.0) go to 152
   if (mfh.eq.1) go to 130
   call tosend(ngendf,0,0,scr)
   go to 115
  130 continue
   iza=nint(c1h)
   nuclid=iza*10
   awr=c2h
   nsigz=l2h
   ntw=n2h
   call listio(ngendf,0,0,scr,nb,nw)
   itemp=nint(c1h)
   if (idtemp.eq.itemp) go to 132
   if (ifound.gt.0) go to 220
   call tomend(ngendf,0,0,scr)
   go to 115
  132 continue
   ifound=1
   idm=nuclid
   idma=idtemp
   idb=0
   ida=0
   idf=0
   id1=0
   id2=0
   id3=0
   id4=0
   id5=0
   nlm=0
   ia=1
   do while (nb.ne.0)
      ia=ia+nw
      call moreio(ngendf,0,0,scr(ia),nb,nw)
      if ((ia+nw-1).gt.nwscr)&
        call error('therm','storage exceeded.',' ')
   enddo
   valu(1)=awr
   write(nsyso,'(&
     &'' atomic weight ........................ '',1p,e10.3)')&
     valu(1)
   ngn=l1h
   if (ngn.ne.ngnd.and.ngn.ne.ngnf)&
     call error('therm','ngn is not equal to ngnd.',' ')
   if (igprnt.le.0) then
      ia=1+lz+ntw+nsigz
      lim=ngn+1
      if (ngn.eq.ngnd) then
         do i=1,lim
            grp(i)=scr(i-1+ia)
         enddo
      else
         do i=1,lim
            ig=icgrp(i)
            if (ig.ne.0) then
               grp(ig)=scr(ia+i-1)
            endif
         enddo
      endif
      write(nsyso,'(/&
        &''24h neutron group structure''/&
        &'' -----------------------''/&
        &(1p,6e12.4))') (grp(i),i=1,ngnd1)
      igprnt=1
   endif
   do i=1,ngnd
      cflux(i)=0
      cflux(i+ngnd)=0
      s0(i)=0
      s1(i)=0
      xa(i)=0
      xf(i)=0
   enddo
   go to 115
  152 continue
   call listio(ngendf,0,0,scr,nb,nw)
   ia=1
   do while (nb.ne.0)
      ia=ia+nw
      call moreio(ngendf,0,0,scr(ia),nb,nw)
      if ((ia+nw-1).gt.nwscr)&
         call error('therm','storage exceeded.',' ')
   enddo
   ig=n2h
   igc=icgrp(ig)
   if (mf3mt2.gt.0) go to 160
   if (iflx0.eq.0) go to 160
   if (mfh.eq.3.and.mth.eq.1) go to 160
   mf3mt2=1
   do i=1,ngnd
      if (cflux(i).ne.zero) cflux(i)=1/cflux(i)
      if (cflux(ngnd+i).ne.zero) cflux(ngnd+i)=1/cflux(ngnd+i)
   enddo
  160 continue
   if (mfh.eq.6) go to 185
   if (igc.eq.0.or.igc.gt.ngnd) go to 166
   if (mth.eq.1) go to 164
   if (mth.gt.1.and.iflx0.eq.0)&
     call error('therm','p0 flux not found.',' ')
   if (mth.eq.18) go to 170
   if (mth.eq.19.and.i318.eq.0) go to 170
   if (mth.ge.102.and.mth.le.150) go to 180
   call tosend(ngendf,0,0,scr)
   go to 115
   ! coarse group flux
   ! p0
  164 continue
   cflux(igc)=cflux(igc)+scr(1+lz)
   iflx0=1
   if (nl.gt.1) then
      ! p1
      cflux(ngnd+igc)=cflux(ngnd+igc)+scr(lz+2)
      iflx1=1
   endif
   go to 166
   ! fission
  170 continue
   idf=1
   if (mth.eq.18) i318=1
   xf(igc)=xf(igc)+scr(1+lz+nl*nz)*scr(1+lz)*cflux(igc)
   go to 166
   ! capture
  180 continue
   xa(igc)=xa(igc)+scr(1+lz+nl*nz)*scr(1+lz)*cflux(igc)
   go to 166
  166 continue
   if (ig.lt.ngn) go to 152
   go to 115
  185 continue
   ! scattering matrix
   if (ida.le.0) then
      write(nsyso,'(/&
        &'' scattering cross sections and mubar''/&
        &''        sig0        sig1       mubar''/&
        &'' -----------------------------------'')')
      ida=1
   endif
   ig2lo=l2h
   ng2=l1h
   if (nl.gt.nlm) nlm=nl
   ig=n2h
   igc=icgrp(ig)
   if (igc.ne.0.and.igc.le.ngnd) then
      do i=2,ng2
         loca=1+(il-1)+lz+nl*nz*(i-1)
         ig2=ig2lo+i-2
         ig2c=icgrp(ig2)
         if (ig2c.ne.0.and.ig2c.le.ngnd) then
            loc=ig2c+ngnd*(igc-1)
            p(loc)=p(loc)+scr(loca)*scr(1+lz)*cflux(igc)
            s0(igc)=s0(igc)+scr(loca)*scr(1+lz)*cflux(igc)
            if (nl.ne.1) then
               ndexc=ngnd+igc
               if (iflx1.eq.0) ndexc=igc
               ndexf=1+lz+1
               if (iflx1.eq.0) ndexf=1+lz
               s1(igc)=s1(igc)&
                 +scr(loca+1)*scr(ndexf)*cflux(ndexc)
            endif
         endif
      enddo
   endif
   if (ig.lt.ngn) go to 152
   go to 115

   !--prepare final matrix
  220 continue
   if (ifound.eq.0) then
      write(strng,'(''mat='',i4,'' temp '',i2,'' not found'')')&
         matd,idtemp
      call error('therm',strng,' ')
   endif
   if (ida.ne.0) then
      do ig=1,ngnd
         if (s0(ig).eq.zero) amubar=0
         if (s0(ig).ne.zero) amubar=s1(ig)/s0(ig)
         write(nsyso,'(1p,6e12.4)') s0(ig),s1(ig),amubar
      enddo
      ! convert to xsec per unit energy.
      do ig=1,ngnd
         do ig2=1,ngnd
            p(ig2+ngnd*(ig-1))=p(ig2+ngnd*(ig-1))&
              *(grp(1+ig)-grp(ig))/(grp(1+ig2)-grp(ig2))
         enddo
      enddo
      ! pack scattering matrix
      call packa(p,ngnd,b,nbm,igi,ngi)
      id3=nbm
   endif

   !--compute absorption and transport cross sections
   do i=1,ngnd
      if (idf.ne.0) xa(i)=xa(i)+xf(i)
      if (ida.ne.0) s1(i)=xa(i)+s0(i)-s1(i)
   enddo
   if (valu(10).eq.0.) valu(10)=s0(ngnd)

   !--write out data for this nuclide.
   if (itrc.eq.0) valu(4)=1-valu(4)
   if (itrc.eq.1) valu(4)=-1
   write(nout,'(2i10,10i5)')&
     idm,idma,idb,ida,idf,id1,id2,id3,id4,id5
   write(nout,'(2a4,a2,6f10.4)') aname,(valu(i),i=1,6)
   write(nout,'(7f10.4)') (valu(i),i=7,13)
   write(nout,'(7f10.4)') (valu(i),i=14,20)
   write(nout,'(1p,6e12.5)') (xa(i),i=1,ngnd)
   if (idf.eq.1) write(nout,'(1p,6e12.5)') (xf(i),i=1,ngnd)
   if (itrc.eq.1) write(nout,'(1p,6e12.5)') (s1(i),i=1,ngnd)
   if (ida.eq.1) write(nout,'(12i6)') (igi(i),i=1,ngnd)
   if (ida.eq.1) write(nout,'(12i6)') (ngi(i),i=1,ngnd)
   if (ida.eq.1) write(nout,'(1p,6e12.5)') (b(i),i=1,nbm)
   write(nsyso,'(/&
     &'' storage pointers''/'' ----------------''/&
     &6x,''idm'',8x,i9/&
     &6x,''idma'',7x,i9/&
     &6x,''idb'',8x,i9/&
     &6x,''ida'',8x,i9/&
     &6x,''idf'',8x,i9/&
     &6x,''id1'',8x,i9/&
     &6x,''id2'',8x,i9/&
     &6x,''id3'',8x,i9)')&
     idm,idma,idb,ida,idf,id1,id2,id3
   if (iprint.ne.0) then
      write(nsyso,'(/&
        &'' absorption cross section''/&
        &'' ------------------------'')')
      write(nsyso,'(1p,6e12.4)') (xa(i),i=1,ngnd)
      if (idf.ne.0) then
         write(nsyso,'(/&
           &'' fission cross section''/&
           &'' ---------------------'')')
         write(nsyso,'(1p,6e12.4)') (xf(i),i=1,ngnd)
      endif
      if (ida.ne.0) then
         if (itrc.eq.1) write(nsyso,'(//&
           &'' transport cross section''/&
           &'' -----------------------'')')
         if (itrc.eq.1) write(nsyso,'(1p,6e12.4)')&
           (s1(i),i=1,ngnd)
         write(nsyso, '(//&
           &'' scattering matrix''/&
           &'' -----------------'')')
         call tprnt(b,igi,ngi,ngnd)
      endif
   endif
   go to 100

   !--therm is finished.
  300 continue
   call repoz(nout)
   return
   end subroutine therm

   subroutine tprnt(a,igi,ngi,nix)
   !--------------------------------------------------------------------
   ! Matrix printing routine.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   ! externals
   integer::nix
   real(kr)::a(*)
   integer::igi(*),ngi(*)
   ! internals
   integer::ncol1,i,k,jlo,ng,lim,j
   integer::ncol=5

   ncol1=ncol-1
   write(nsyso,'(/&
     &''46h final initl    cross section vs initial group''/&
     &'' group group    +0'',4(10x,''+'',i1)/&
     &1x,7(''----------''),''--------''/)')&
     (i,i=1,ncol1)
   k=0
   do 100 i=1,nix
   jlo=igi(i)
   ng=ngi(i)
   lim=ng
   if (lim.gt.ncol) lim=ncol
   write(nsyso,'(i5,i6,2x,1p,5e12.4)')&
     i,jlo,(a(k+j),j=1,lim)
  110 continue
   k=k+lim
   if (lim.eq.ng) go to 100
   jlo=jlo+lim
   ng=ng-lim
   lim=ng
   if (lim.gt.ncol) lim=ncol
   write(nsyso,'(5x,i6,2x,1p,5e12.4)')&
     jlo,(a(k+j),j=1,lim)
   go to 110
  100 continue
   return
   end subroutine tprnt

   subroutine packa(tpa,igroup,a,k,igi,ng)
   !--------------------------------------------------------------------
   ! Pack a scattering matrix into banded form by removing
   ! unneeded zeroes.
   !--------------------------------------------------------------------
   ! externals
   integer::igroup
   real(kr)::tpa(igroup,igroup),a(*)
   integer::igi(*),ng(*)
   ! internals
   integer::k,i,idone,j,jj,jhi
   real(kr),parameter::small=1.e-6_kr

   k=0
   do i=1,igroup
      igi(i)=0
      ng(i)=0
      idone=0
      j=0
      do while (j.lt.igroup.and.idone.eq.0)
         j=j+1
         jj=igroup-j+1
         if (abs(tpa(i,jj)).ge.small) idone=1
      enddo
      jhi=jj
      do j=1,jhi
         if (abs(tpa(i,j)).ge.small.or.igi(i).ne.0) then
            if (abs(tpa(i,j)).ge.small.and.igi(i).eq.0) igi(i)=j
            k=k+1
            a(k)=tpa(i,j)
         endif
      enddo
      ng(i)=jhi-igi(i)+1
   enddo
   return
   end subroutine packa

   subroutine cpm(iprnt,iclaps)
   !--------------------------------------------------------------------
   ! Format multigroup cross sections from groupr for the EPRI-CPM
   ! library code CLIB.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use util ! provides error,repoz
   ! externals
   integer::iprnt,iclaps
   ! internals
   integer::i,k,nina,ires,mti,mtc,ip1opt,norf,ipos,iz,iposr
   integer::nfis1,nfis2,j,nz,identa,index,izero,iscr,jscr
   integer::nw,nwflx,nn,nwglam,nlib,idatt,iopt,ntemp,nt
   integer::immax,lspace
   real(kr)::sgref,sigp
   integer::im(70)
   real(kr)::z(70)
   integer::km(70)
   real(kr),dimension(:),allocatable::scr2
   real(kr),dimension(:),allocatable::f1p
   integer,parameter::nrg69=13
   integer,parameter::nrg185=15
   real(kr),parameter::zero=0

   ngnd=69
   nfg=14
   nrg=nrg185
   if (iclaps.gt.0) nrg=nrg69
   nfsg=22
   ngref=60
   if (nrg.eq.nrg69) ngref=55
   nt=6
   nzmax=10
   iprint=iprnt
   ! arrays im and rleth are available for temporary storage
   immax=ngnd+1+nrg+1
   iscr=1

   !--read and write user input.
   write(nsyso,'(/'' ***epri-cpm library***'')')
   write(nsyse,'(/'' ***epri-cpm library***'')')
   mode=0
   if5=0
   if4=1
   read(nsysi,*) nlib,idatt,newmat,iopt,mode,if5,if4
   do i=1,newmat
      im(i)=0
   enddo
   if (iopt.ne.1) then
      if (newmat.gt.immax)&
        call error('cpm','exceeded storage in input array z.',' ')
      read(nsysi,*) (im(i),i=1,newmat)
   endif
   write(nsyso,'(&
     &'' number of library .................... '',i10/&
     &'' date ................................. '',i10/&
     &'' no. of materials to be added ......... '',i10/&
     &'' add option (0 only mats read, 1 all) . '',i10/&
     &'' mode (0 replace, 1 add, 2 create ) ... '',i10/&
     &'' burn option (0 no, 1 yes, 2 only) .... '',i10/&
     &'' xsec data option (0 no, 1 yes) ....... '',i10)')&
     nlib,idatt,newmat,iopt,mode,if5,if4
   if (iopt.eq.0) write(nsyso,&
     &'('' desired mats ......................... ''/&
     &(3x,12i6))') (im(i),i=1,newmat)
   if (iopt.eq.1) write(nsyso,&
     &'('' first '',i2,'' mat(s) on ngendf will be processed.'')')&
     newmat
   write(nsyso,'(/&
     &'' nina ntemp nsigz ref. sig0 ires pot. xsec '',&
     &''mti mtc ip1 norf pos posr''/&
     &'' ---- ----- ----- --------- ---- --------- '',&
     &''--- --- --- ---- --- ----'')')
   nres=0
   nnina=0
   k=0
   do i=1,newmat
      nina=0
      ntemp=0
      nsigz=0
      sgref=0
      ires=0
      sigp=0
      mti=0
      mtc=0
      ip1opt=0
      norf=0
      ipos=0
      iposr=0
      read(nsysi,*) nina,ntemp,nsigz,sgref,ires,sigp,mti,mtc,&
        ip1opt,norf,ipos,iposr
      if (nina.eq.0.or.nina.eq.3) nnina=nnina+1
      if (ires.eq.1) nres=nres+1
      if (nina.ne.0) ip1opt=1
      write(nsyso,'(i4,i5,i6,1p,e12.2,i4,e11.3,2i4,i4,2i5,i4)')&
        nina,ntemp,nsigz,sgref,ires,sigp,mti,mtc,&
        ip1opt,norf,ipos,iposr
      k=k+12
   enddo
   if (if5.ne.0) then
      ! file 5 input
      iz=newmat+1
      read(nsysi,*) ntis,nfis
      lspace=immax-newmat
      nfis1=nfis+1
      nfis2=nfis+2
      if (lspace.lt.nfis2) allocate(scr2(nfis2))
      allocate(idat(ntis))
      allocate(idbt(ntis))
      nw=ntis*nfis1
      allocate(yld(nw))
      do i=1,nw
         yld(i)=0
      enddo
      read(nsysi,*) (idbt(i),i=1,nfis)
      do i=1,nfis
         if (lspace.ge.nfis2) im(iz-1+i)=idbt(i)
         if (lspace.lt.nfis2) km(i)=idbt(i)
      enddo
      if (lspace.ge.nfis2) write(nsyso,'(&
        &'' no. time-dependent nuclides .......... '',i10/&
        &'' no. fissionable burnup nuclides ...... '',i10/&
        &'' fissionable burnup nuclide ids ....... ''//&
        &(16x,9i12))')&
        ntis,nfis,(im(iz-1+i),i=1,nfis)
      if (lspace.lt.nfis2) write(nsyso,'(&
        &'' no. time-dependent nuclides .......... '',i10/&
        &'' no. fissionable burnup nuclides ...... '',i10/&
        &'' fissionable burnup nuclide ids ....... ''//&
        &(16x,9i12))')&
        ntis,nfis,(km(i),i=1,nfis)
      write(nsyso,'(&
        &/''  ident    decay     yields''/&
        &'' ------'',2x,''----------'',2x,''------'')')
      do i=1,ntis
         do j=1,nfis2
            if (lspace.ge.nfis2) z(iz-1+j)=0
            if (lspace.lt.nfis2) scr2(i)=0
         enddo
         nz=nfis2
         if (lspace.ge.nfis2) read(nsysi,*) (z(iz+j-1),j=1,nz)
         if (lspace.lt.nfis2) read(nsysi,*) (scr2(j),j=1,nz)
         if (lspace.ge.nfis2) identa=nint(z(iz))
         if (lspace.lt.nfis2) identa=nint(scr2(1))
         idat(i)=identa
         index=1+(i-1)*nfis1
         izero=0
         do j=2,nz
            if (lspace.ge.nfis2) yld(index+j-2)=z(iz+j-1)
            if (lspace.lt.nfis2) yld(index+j-2)=scr2(j)
            if (j.gt.2.and.yld(index+j-2).ne.zero) izero=1
         enddo
         if (izero.gt.0) write(nsyso,'(i7,1p,10e12.4/(19x,9e12.4))')&
           identa,(yld(index+j-1),j=1,nfis1)
         if (izero.eq.0) write(nsyso,'(i7,1p,10e12.4/(19x,9e12.4))')&
           identa,yld(index)
      enddo
      if (if5.eq.1.and.lspace.lt.nfis2) deallocate(scr2)
      if (if5.ne.1) then
         nw=3*ntis
         allocate(f1p(nw))
         write(nsyso,'(&
           &3x,''ident'',5x,''aw'',5x,''indfis'',2x,''ntemp''/&
           &3x,''-----'',2x,''--------'',2x,''------'',2x,''-----'')')
         if (lspace.ge.nfis2) then
            do i=1,ntis
               index=1+(i-1)*3
               read(nsysi,*) (f1p(index+j-1),j=1,3)
               identa=idat(i)
               write(nsyso,'(3x,i6,f10.4,i6,i7)')&
                 identa,f1p(index),f1p(index+1),f1p(index+2)
            enddo
         else
            do i=1,ntis
               index=1+(i-1)*3
               read(nsysi,*) (f1p(index+j-1),j=1,3)
               identa=idat(i)
               write(nsyso,'(3x,i6,f10.4,i6,i7)')&
                 identa,f1p(index),f1p(index+1),f1p(index+2)
            enddo
            deallocate(scr2)
         endif
      endif
   endif
   write(nsyso,'(/)')

   !--allocate storage.
   allocate(mat(newmat))
   allocate(idnt(newmat))
   allocate(ninat(newmat))
   allocate(irest(newmat))
   allocate(s0rf(newmat))
   allocate(ip1tab(newmat))
   allocate(norft(newmat))
   allocate(mtit(newmat))
   allocate(mtct(newmat))
   allocate(awrt(newmat))
   allocate(indf(newmat))
   allocate(ntmp(newmat))
   allocate(nszt(newmat))
   nwflx=nrg*nt*nzmax
   if (nwflx.lt.ngnd*2) nwflx=ngnd*2
   allocate(flux(nwflx))
   allocate(sigpt(newmat))
   nn=nnina
   if (nn.eq.0) nn=nres
   nwglam=nn*nrg
   if (nwglam.gt.0) allocate(glamt(nwglam))
   if (mode.eq.0) then
      allocate(ipost(newmat))
      allocate(iposrt(newmat))
   endif
   if (nnina.ne.0.or.nres.ne.0) then
      jscr=iscr+k+1
      do i=1,nn
         do j=1,nrg
            glamt(i+nn*(j-1))=0
         enddo
         read(nsysi,*) (scr(jscr+k-1),k=1,nrg)
      enddo
   endif
   k=0
   do i=1,newmat
      mat(i)=im(i)
      ninat(i)=nint(scr(iscr+k))
      ntmp(i)=nint(scr(iscr+k+1))
      nszt(i)=nint(scr(iscr+k+2))
      s0rf(i)=scr(iscr+k+3)
      irest(i)=nint(scr(iscr+k+4))
      sigpt(i)=scr(iscr+k+5)
      mtit(i)=nint(scr(iscr+k+6))
      mtct(i)=nint(scr(iscr+k+7))
      ip1tab(i)=nint(scr(iscr+k+8))
      norft(i)=nint(scr(iscr+k+9))
      k=k+12
   enddo
   if (mode.le.0) then
      k=0
      do i=1,newmat
         ipost(i)=nint(scr(iscr+k+10))
         iposrt(i)=nint(scr(iscr+k+11))
         k=k+12
      enddo
   endif
   nscr=10
   kscr=11
   mscr=12
   call repoz(nout)
   call pinit(iopt)

   !--process data for resonance calculation.
   call sfile2

   !--process effective resonance integrals.
   call sfile3(ntemp,nt)

   !--process cross sections.
   call sfile4(ntemp,nt)

   !--process burnup data.
   if (if5.gt.0) call sf5out

   !--process p1 scattering matrix.
   call sfile6(ntemp)

   !--write output tape.
   if (nout.ne.0) call cpmout(nlib,idatt)

   !--cpm is finished.
   call repoz(nout)
   return
   end subroutine cpm

   subroutine pinit(iopt)
   !--------------------------------------------------------------------
   ! Store materials from ngendf.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use endf ! provides endf routines and variables
   use util ! provides repoz,error,mess
   ! externals
   integer::iopt
   ! internals
   integer::ngnd1,nb,nw,matl,nmat,n2,i,nina,icont,mat1,irepoz
   integer::matd,iza,nz,ntw,i1,ngn1,lastg,ig,jg,i318,mess3,i618
   integer::mess6,il,ih,jgref
   real(kr)::aleth
   character(60)::strng

   ngnd1=ngnd+1
   allocate(egb(ngnd1))
   iu235=92235
   ipu239=94239
   iu=0
   nfiss(1)=0
   ipu=0
   nfiss(2)=0
   if (iopt.eq.1) write(nsyso,&
     '(/'' materials to be processed ............ '')')
   call repoz(ngendf)
   call tpidio(ngendf,0,0,scr,nb,nw)
   matl=0
   nmat=0
   if (nnina.gt.0) go to 100
   n2=0
   do i=1,newmat
      nina=ninat(i)
      if (nina.eq.2) n2=n2+1
   enddo
   if (n2.eq.newmat) go to 110
  100 continue
   nmat=nmat+1
   if (nmat.gt.newmat) go to 180
   icont=0
   mat1=0
   irepoz=0
   idnt(nmat)=0
   if (mat(nmat).lt.0) idnt(nmat)=100
   if (idnt(nmat).gt.0) mat(nmat)=-mat(nmat)
   if (iopt.eq.0) matd=mat(nmat)
   nina=ninat(nmat)
   if (nina.lt.2) go to 110
   idnt(nmat)=matd+idnt(nmat)
   go to 100
  110 continue
   call contio(ngendf,0,0,scr,nb,nw)
   if (math.eq.0) go to 110
   if (math.lt.0) go to 170
   icont=icont+1
   if (iopt.eq.0.and.icont.eq.1) mat1=math
   if (math.eq.matl) go to 160
   matl=math
   if (iopt.eq.0.and.irepoz.gt.0.and.math.eq.mat1) go to 170
   if (iopt.eq.0.and.math.ne.matd) go to 160
   if (iopt.eq.1) mat(nmat)=math
   if (iopt.eq.1) write(nsyso,'(8x,i5)') math
   iza=nint(c1h)
   idnt(nmat)=iza+idnt(nmat)
   awrt(nmat)=c2h
   if (iza.eq.iu235) iu=1
   if (iza.eq.ipu239) ipu=1
   if (nmat.le.1) then
      nz=l2h
      ntw=n2h
      call listio(ngendf,0,0,scr,nb,nw)
      i1=7+ntw+nz
      ngn=l1h
      if (ngn.eq.ngnd) then
         do i=1,ngnd1
            egb(i)=scr(ngnd1-i+i1)
         enddo
      else
         ! calculate group bounds for collapsed structure
         if (ngn.ne.ngnf)&
           call error('init','wrong group structure.',' ')
         ngn1=ngn+1
         lastg=0
         do i=1,ngn1
            ig=icgrp(i)
            if (ig.ne.lastg.and.ig.ne.0) then
               jg=ngnd-ig+1
               egb(1+jg)=scr(i1+i-1)
            endif
         enddo
      endif
   endif
   ! check for presence of mf3, mt18 and mf6, mt18.
   i318=0
   mess3=0
   i618=0
   mess6=0
  125 continue
   call tosend(ngendf,0,0,scr)
  130 continue
   call contio(ngendf,0,0,scr,nb,nw)
   if (math.eq.0) go to 150
   if (math.eq.-1) go to 170
   if (math.ne.matd) go to 145
   if (mfh.eq.0.or.mth.eq.0) go to 130
   if (mfh.eq.3.and.mth.ge.18.and.mth.le.38) go to 135
   if (mfh.eq.6.and.mth.ge.18.and.mth.le.38) go to 140
   go to 125
  135 continue
   if (mth.gt.21.and.mth.lt.38) go to 125
   if (mth.eq.18) i318=1
   if (mth.eq.18) go to 125
   if (i318.eq.0) go to 125
   if (mess3.gt.0) go to 125
   write(strng,&
     '(''mat='',i4,'' mf='',i2,'' has both mt18 and mt'',i2)')&
      math,mfh,mth
   call mess('pinit',strng,'mt18 used in all subroutines')
   mess3=1
   if (mess3.gt.0.and.mess6.gt.0) go to 150
   go to 125
  140 continue
   if (mth.gt.21.and.mth.lt.38) go to 125
   if (mth.eq.18) i618=1
   if (mth.eq.18) go to 125
   if (i618.eq.0) go to 125
   if (mess6.gt.0) go to 125
   write(strng,&
     '(''mat='',i4,'' mf='',i2,'' has both mt18 and mt'',i2)')&
     math,mfh,mth
   call mess('pinit',strng,'mt18 used in all subroutines')
   mess6=1
   if (mess3.gt.0.and.mess6.gt.0) go to 150
   go to 125
  145 continue
   call skiprz(ngendf,-3)
   math=matd
  150 continue
   if (nmat.eq.newmat) go to 180
   if (math.eq.0) go to 100
   call tomend(ngendf,0,0,scr)
   go to 100
  160 continue
   call tomend(ngendf,0,0,scr)
   go to 110
  170 continue
   if (irepoz.gt.0.or.iopt.gt.0) go to 175
   call repoz(ngendf)
   call tpidio(ngendf,0,0,scr,nb,nw)
   irepoz=1
   matl=0
   go to 110
  175 continue
   if (iopt.eq.0) then
      write(strng,'(''mat='',i4,'' not found on gendf'')') matd
      call error('pinit',strng,' ')
   endif
   nmat=nmat-1
   write(strng,'(''only '',i2,'' mats found on gendf'')') nmat
   call mess('pinit',strng,' ')
   newmat=nmat
  180 continue
   if (ngn.eq.ngnd) write(nsyso,&
     &'(/'' group structure''/'' ---------------''/&
     &(1p,6e12.4))') (egb(i),i=1,ngnd1)
   if (ngn.ne.ngnd) write(nsyso,'(/&
     &'' group structure (collapsed from '',i3,&
     &'' group structure)''/&
     &'' ----------------------------------------------------''/&
     &(1p,6e12.4))')&
     ngn,(egb(i),i=1,ngnd1)
   ! find lethargy factors for resonance groups
   il=nfg+1
   ih=nfg+nrg
   do i=il,ih
      aleth=log(egb(i)/egb(i+1))
      rleth(i-nfg)=1/aleth
   enddo
   ! find lethargy factor for ngref
   jgref=ngnd1-ngref
   aleth=log(egb(jgref)/egb(1+jgref))
   reflth=aleth
   if (mode.lt.2) deallocate(egb)
   return
   end subroutine pinit

   subroutine sfile2
   !--------------------------------------------------------------------
   ! Process data for resonance calculation.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use endf ! provides endf routines and variables
   use util ! provides openz,repoz,error
   ! internals
   integer::nghi,nglo,nw,nmat,matd,nina,isg,mf3mt2
   integer::i318,i618,i,nl,nz,ntw,l,ng2,ig,igc,kg,jg,isgz
   integer::loc,loca,locf,jgc,j,jnina,nb,nwloc,iloc
   integer::loctot,locnus,locsf0
   real(kr)::sgref,sigz
   integer::lz=6
   real(kr),dimension(:),allocatable::a
   real(kr),parameter::eps=1.e-5_kr
   real(kr),parameter::dilinf=1.e10_kr
   real(kr),parameter::small=1.e-9_kr
   real(kr),parameter::eflag=1.e35_kr
   real(kr),parameter::zero=0

   !--allocate storage.
   call openz(-nscr,1)
   if (nnina.eq.0) go to 340
   write(nsyso,'(//)')
   write(nsyso,&
     '(/'' ***subfile 2***resonance data'')')
   nwloc=nrg*3
   allocate(a(nwloc))
   iloc=1
   loctot=iloc
   locnus=loctot+nrg
   locsf0=locnus+nrg
   nghi=ngnd-nfg
   nglo=nghi-nrg+1

   !--search for desired materials.
   call repoz(ngendf)
   call tpidio(ngendf,0,0,scr,nb,nw)
   nmat=0
   jnina=0
  100 continue
   nmat=nmat+1
   matd=mat(nmat)
   nina=ninat(nmat)
   if (nina.eq.3) go to 315
   if (nina.gt.0) go to 320
   jnina=jnina+1
   sgref=s0rf(nmat)
   isg=0
   if (abs(sgref-dilinf).ge.sgref/100) isg=1
   mf3mt2=0
   i318=0
   i618=0
   do i=1,nwloc
      a(i-1+iloc)=0
   enddo
   do i=1,ngnd
      cflux(i)=0
      cflux(i+ngnd)=0
   enddo
  110 continue
   call contio(ngendf,0,0,scr,nb,nw)
   if (math.eq.matd) go to 130
   if (math.eq.0) go to 110
   if (math.gt.0) go to 115
   call repoz(ngendf)
   call tpidio(ngendf,0,0,scr,nb,nw)
   go to 110
  115 continue
   call tomend(ngendf,0,0,scr)
   go to 110

   !--search for required reactions.
  120 continue
   call contio(ngendf,0,0,scr,nb,nw)
   if (math.eq.0) go to 300
  130 continue
   nl=l1h
   nz=l2h
   ntw=n2h
   if (mfh.eq.0.or.mth.eq.0) go to 120
   if (math.ne.matd) go to 300
   if (mfh.ne.3.or.mth.lt.2.or.mf3mt2.gt.0) go to 140
   mf3mt2=1
   do i=1,ngnd
      if (cflux(i).ne.zero) cflux(i)=1/cflux(i)
      if (cflux(ngnd+i).ne.zero) cflux(ngnd+i)=1/cflux(ngnd+i)
   enddo
  140 continue
   l=1
   call listio(ngendf,0,0,scr(l),nb,nw)
   if (mfh.eq.1) go to 165
   ng2=l1h
   ig=n2h
   do while (nb.ne.0)
      l=l+nw
      call moreio(ngendf,0,0,scr(l),nb,nw)
      if ((l+nw-1).gt.nwscr)&
        call error('sfile2','storage exceeded.',' ')
   enddo
   if (ig.eq.0) go to 140
   igc=icgrp(ig)
   if (igc.eq.0.or.igc.gt.ngnd) go to 250
   if (igc.lt.nglo.or.igc.gt.nghi) go to 250
   kg=ngnd-igc+1
   jg=kg-nfg
   l=1
   if (mth.eq.455) go to 240
   if (mfh.eq.3) go to 200
   if (mfh.eq.6) go to 220
   go to 170
  165 continue
   isgz=0
   do i=1,nz
      sigz=scr(6+ntw+i)
      if (abs(sigz-sgref).le.eps*sigz) isgz=i
   enddo
   if (isgz.eq.0)&
     call error('sfile2','reference sigma zero not found.',' ')
  170 continue
   call tosend(ngendf,0,0,scr)
   go to 120

   !--file 3 cross sections.
  200 continue
   if (mth.eq.18) i318=1
   if (mth.gt.18.and.i318.eq.1) go to 170
   if (mth.ne.1.and.mth.ne.2.and.&
     (mth.lt.18.or.mth.gt.21).and.mth.ne.38) go to 170
   ! fission
   if ((mth.ge.18.and.mth.le.21).or.mth.eq.38) then
      loc=locsf0+jg-1
      loca=l+lz+nl*nz
      locf=l+lz
      a(loc)=a(loc)+scr(loca)*scr(locf)*cflux(jg)
   ! elastic
   else if (mth.eq.2) then
      loc=loctot+jg-1
      loca=l+lz+nl*nz
      if (isg.gt.0.and.nz.ge.isgz) loca=l+lz+nl*(isgz-1+nz)
      locf=l+lz
      if (isg.gt.0.and.nz.ge.isgz) locf=l+lz+nl*(isgz-1)
      jgc=jg
      if (isg.gt.0.and.nz.ge.isgz) jgc=jg+ngnd
      a(loc)=a(loc)+scr(loca)*scr(locf)*cflux(jgc)
   ! total
   else if (mth.eq.1) then
      ! coarse group flux
      ! for infinite dilution
      cflux(jg)=cflux(jg)+scr(l+lz)
      if (isg.ne.0.and.nz.ge.isgz) then
         ! for reference sigma zero
         cflux(jg+ngnd)=cflux(jg+ngnd)+scr(l+lz+nl*(isgz-1))
      endif
   endif
   go to 250

   !--file 6 matrices.
  220 continue
   if ((mth.lt.18.or.mth.gt.21).and.mth.ne.38) go to 170
   if (mth.eq.18) i618=1
   if (mth.gt.18.and.i618.eq.1) go to 170
   ! fission
   loc=locnus+jg-1
   locf=l+lz
   loca=l+lz-nl*nz
   do i=2,ng2
      a(loc)=scr(nl*nz*i+loca)*scr(locf)*cflux(jg)+a(loc)
   enddo
   go to 250
  240 continue
   if (mfh.eq.5) go to 170
   ! delayed neutrons
   loc=locnus+jg-1
   loca=l+lz
   locf=l+lz
   a(loc)=a(loc)+scr(loca+1)*scr(loca+2)*scr(locf)*cflux(jg)
  250 continue
   if (ig.lt.ngn) go to 140
   go to 120

   !--desired data is loaded.
   !--find nu from nu times sigma f.
  300 continue
   do ig=1,nrg
      loc=locsf0+ig-1
      loca=locnus+ig-1
      if (a(loc).ne.zero) then
         a(loca)=a(loca)/a(loc)
      endif
      if (a(loca).lt.small) a(loca)=0
   enddo
   go to 319

   !--read in nus and tot for nuclides having nina=3
  315 continue
   read(nsysi,*) (a(locnus+j-1),j=1,nrg)
   read(nsysi,*) (a(loctot+j-1),j=1,nrg)

  319 continue
   call sf2out(nmat,jnina,loctot,locnus,a)
  320 continue
   if (nmat.lt.newmat) go to 100

   !--sfile2 is finished.
   if (iprint.eq.1) write(nsyso,'(/)')
  340 continue
   write(nscr) eflag
   return
   end subroutine sfile2

   subroutine sf2out(nmat,jnina,loctot,locnus,a)
   !--------------------------------------------------------------------
   ! Write subfile 2.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   ! externals
   integer::nmat,jnina,loctot,locnus
   real(kr)::glam
   real(kr)::a(*)
   ! internals
   integer::ident,i,ig
   real(kr)::sigp

   ident=idnt(nmat)
   if (iprint.eq.1) write(nsyso,&
     '(/'' nuclide '',i6/'' -------'')') ident
   if (iprint.eq.1) write(nsyso,&
     '(/''  group'',7x,''lambda'',8x,''sigp'',10x,''nu'',9x,&
     &''sigsctot''/1x,6(''----------''),''-'')')
   sigp=sigpt(nmat)
   do i=1,nrg
      ig=i+nfg
      glam=glamt(jnina+nnina*(i-1))
      write(nscr) glam,sigp,a(locnus-1+i),a(loctot-1+i)
      if (iprint.eq.1)&
        write(nsyso,'(i5,5x,1p,4e13.4)')&
        ig,glam,sigp,a(locnus-1+i),a(loctot-1+i)
   enddo
   return
   end subroutine sf2out

   subroutine sfile3(jtemp,nt)
   !--------------------------------------------------------------------
   ! Process effective resonance integrals.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use endf ! provides endf routines and variables
   use util ! provides repoz,error
   ! externals
   integer::jtemp,nt
   ! internals
   integer::nn,nghi,nglo,nb,nw,il,nmat,irmat,jnina,matd,nina
   integer::inorf,jres,ntnmat,msigz,i,nwflx,jfiss,nl,nz,ntw,nwfa
   integer::l,nsigz,index,i318,i618,ig,jg,iadd,nwflxr
   integer::loca,locb,jz,locf,ng2,iterm,jtem,it,lim
   integer::indexl,ix,igc,kg,loc,locabs,locnus,locsf,nwloc,iloc
   integer::locs(7)
   real(kr)::temp
   equivalence (locabs,locs(1)),(locnus,locs(2)),(locsf,locs(3))
   real(kr)::xid,amda,sigp
   real(kr)::b(70)
   real(kr),dimension(:),allocatable::tempt,sigz
   real(kr),dimension(:),allocatable::flxr
   real(kr),dimension(:),allocatable::abs
   real(kr),dimension(:),allocatable::snu
   real(kr),dimension(:),allocatable::fa
   real(kr),dimension(:),allocatable::a
   integer::lz=6
   real(kr),parameter::eflag=1.e35_kr
   real(kr),parameter::zero=0

   !--allocate storage.
   if (nres.eq.0) go to 510
   write(nsyso,'(//)')
   write(nsyso,&
     '(/'' ***subfile 3***resonance integrals'')')
   nw=nt+1
   allocate(tempt(nw))
   nwloc=nzmax*nt*nrg*3
   allocate(a(nwloc))
   iloc=1
   lim=nwloc-1
   nwflxr=nzmax*nt
   allocate(flxr(nwflxr))
   allocate(sigz(nzmax))
   allocate(snu(nrg))
   allocate(abs(nt))
   nwfa=0
   if ((nzmax*nt).gt.(nrg+1)) then
      nwfa=nzmax*nt
      allocate(fa(nwfa))
   endif
   locabs=iloc
   locnus=locabs+nzmax*nt*nrg
   locsf=locnus+nzmax*nt*nrg
   nn=nnina
   if (nn.eq.0) nn=nres
   nghi=ngnd-nfg
   nglo=nghi-nrg+1

   !--search for desired materials.
   call repoz(ngendf)
   call tpidio(ngendf,0,0,scr,nb,nw)
   il=1
   nmat=0
   irmat=0
   jnina=0
  110 continue
   nmat=nmat+1
   matd=mat(nmat)
   nina=ninat(nmat)
   if (nina.eq.0) jnina=jnina+1
   indf(nmat)=0
   inorf=norft(nmat)
   jres=irest(nmat)
   if (jres.le.0) go to 470
   if (nnina.eq.0) jnina=jnina+1
   ntnmat=ntmp(nmat)
   msigz=nszt(nmat)
   irmat=irmat+1
   do i=1,nwflxr
      flxr(i)=0
   enddo
   nwflx=nzmax*nt*nrg
   do i=1,nwflx
      flux(i)=0
   enddo
   do i=1,nwloc
      a(i-1+iloc)=0
   enddo
   do i=1,nrg
      snu(i)=0
   enddo
   jtemp=0
   jfiss=1
  120 continue
   call contio(ngendf,0,0,scr,nb,nw)
   if (math.eq.matd) go to 130
   if (math.eq.0) go to 120
   if (math.gt.0) go to 125
   call repoz(ngendf)
   call tpidio(ngendf,0,0,scr,nb,nw)
   go to 120
  125 continue
   call tomend(ngendf,0,0,scr)
   go to 120

   !--search for required reactions.
  130 continue
   xid=idnt(nmat)
  135 continue
   nl=l1h
   nz=l2h
   ntw=n2h
  140 continue
   l=1
   call listio(ngendf,0,0,scr(l),nb,nw)
   if (mfh.ne.1) go to 160
   temp=c1h
   jtemp=jtemp+1
   tempt(jtemp)=temp
   nsigz=nz
   if (nz.gt.msigz.and.msigz.ne.0) nsigz=msigz
   do i=1,nsigz
      index=1+nsigz-i
      sigz(index)=scr(6+ntw+i)
   enddo
   abs(jtemp)=0
   i318=0
   i618=0
   go to 300
  160 continue
   if (mth.ge.18.and.mth.le.21.and.inorf.gt.0) go to 300
   if (mth.eq.38.and.inorf.gt.0) go to 300
   if (mfh.eq.6.and.jtemp.eq.1) go to 162
   if (mfh.ne.3) go to 300
   if (mth.eq.1) go to 162
   if (mth.eq.455.and.jtemp.eq.1) go to 162
   if (mth.lt.18.or.mth.gt.150) go to 300
   if (mth.gt.21.and.mth.lt.102.and.mth.ne.38) go to 300
   if (mth.eq.18) i318=1
   if (mth.gt.18.and.mth.lt.39.and.i318.gt.0) go to 300
  162 continue
   ig=n2h
   do while (nb.ne.0)
      l=l+nw
      call moreio(ngendf,0,0,scr(l),nb,nw)
      if ((l+nw-1).gt.nwscr)&
        call error('sfile3','storage exceeded.',' ')
   enddo
   if (ig.eq.0) go to 140
   igc=icgrp(ig)
   if (mth.eq.1.and.igc.eq.ngref) go to 172
   if (igc.lt.nglo.or.igc.gt.nghi) go to 185
  172 continue
   kg=ngnd-igc+1
   jg=kg-nfg
   lim=nsigz
   if (nsigz.gt.nz) lim=nz
   if (mfh.eq.6) go to 230
   if (mth.eq.1) go to 190
   if (mth.eq.455) go to 225
   if (mth.eq.18) go to 175
   if (mth.lt.102.and.i318.gt.0) go to 215
   ! absorption
  175 continue
   iadd=nsigz+nzmax*((jtemp-1)+nt*(jg-1))
   loca=l+lz+nl*(nz-1)
   locb=l+lz-nl
   do jz=1,lim
      a(iadd-jz+locabs)=a(iadd-jz+locabs)&
        +scr(nl*jz+loca)*scr(nl*jz+locb)/flux(iadd-jz+1)
   enddo
   if (mth.lt.102) go to 215
   abs(jtemp)=1
  185 continue
   if (ig.lt.ngn) go to 140
   go to 310
   ! flux
  190 continue
   if (igc.ne.ngref) go to 195
   if (ngref.gt.nghi.or.ngref.lt.nglo) go to 205
  195 continue
   loc=1+nsigz+nzmax*((jtemp-1)+nt*(jg-1))
   loca=l+lz-nl
   do jz=1,lim
      flux(loc-jz)=scr(nl*jz+loca)+flux(loc-jz)
   enddo
   if (igc.eq.ngref) go to 205
   go to 185
   ! reference group flux
  205 continue
   do jz=1,lim
      loca=l+lz+nl*(jz-1)
      loc=1+nsigz-jz+nzmax*(jtemp-1)
      flxr(loc)=flxr(loc)+scr(loca)
   enddo
   go to 185
   ! nu sigma f
  215 continue
   iadd=nsigz+nzmax*((jtemp-1)+nt*(jg-1))
   loca=l+lz+nl*(nz-1)
   locb=l+lz-nl
   do jz=1,lim
      a(iadd-jz+locsf)=a(iadd-jz+locsf)&
        +scr(nl*jz+loca)*scr(nl*jz+locb)/flux(iadd-jz+1)
   enddo
   jfiss=2
   go to 185
   ! delayed neutrons
  225 continue
   loca=l+lz
   locf=1+nsigz-1+nzmax*nt*(jg-1)
   snu(jg)=snu(jg)+scr(loca+1)*scr(loca+2)*scr(loca)/flux(locf)
   go to 185
   ! save nu from file 6
  230 continue
   if ((mth.lt.18.or.mth.gt.22).and.mth.ne.38) go to 300
   if (mth.eq.18) i618=1
   if (mth.eq.18) go to 235
   if (i618.eq.0) go to 235
   go to 300
  235 continue
   ng2=l1h
   loca=l+lz+(il-1)-nl*nz
   locf=1+nsigz-1+nzmax*nt*(jg-1)
   do i=2,ng2
      snu(jg)=snu(jg)+scr(nl*nz*i+loca)*scr(l+lz)/flux(locf)
   enddo
   go to 185
  300 continue
   call tosend(ngendf,0,0,scr)
  310 continue
   call contio(ngendf,0,0,scr,nb,nw)
   if (math.lt.0) go to 350
   if (math.eq.0.or.mfh.eq.0.or.mth.eq.0) go to 310
   if (math.ne.matd) go to 350
   if (mfh.eq.1) go to 360
   if (mfh.ne.3.and.mfh.ne.6) go to 300
   go to 135
  350 continue
   if (jtemp.gt.1.or.jfiss.eq.1) go to 400
   go to 370
  360 continue
   if (jtemp.gt.1.or.jfiss.eq.1) go to 130

   !--correct nu
  370 continue
   do jg=1,nrg
      iterm=nsigz-1+nzmax*nt*(jg-1)
      loca=locsf+iterm
      if (a(loca).eq.zero) snu(jg)=0
      if (a(loca).ne.zero) snu(jg)=snu(jg)/a(loca)
   enddo
   if (math.ne.matd) go to 400
   if (ntnmat.eq.0) go to 130
   if (jtemp.eq.ntnmat) go to 400
   go to 130

   !--desired data is loaded.
  400 continue
   if (jfiss.ne.1) then
      ! calculate nu times sigma f
      do jg=1,nrg
         do jz=1,nsigz
            iterm=jz-1+nzmax*nt*(jg-1)-nzmax
            do jtem=1,jtemp
               a(nzmax*jtem+iterm+locnus)=&
                 a(nzmax*jtem+iterm+locsf)*snu(jg)
            enddo
         enddo
      enddo
   endif
   ! check flux arrays for temperature dependence.
   do jz=1,nsigz
      do it=1,jtemp
         if (flxr(jz+nzmax*(it-1)).eq.0.) then
            flxr(jz+nzmax*(it-1))=flxr(jz+nzmax*(it-2))
         endif
      enddo
   enddo
   do jg=1,nrg
      do jz=1,nsigz
         do it=1,jtemp
            index=jz+nzmax*((it-1)+nt*(jg-1))
            if (flux(index).eq.zero) then
               indexl=jz+nzmax*((it-2)+nt*(jg-1))
               flux(index)=flux(indexl)
            endif
         enddo
      enddo
   enddo
   ! correct non-temperature dependent absorption if necessary
   if (jtemp.ne.1) then
      do it=2,jtemp
         if (abs(it).le.zero) then
            do jg=1,nrg
               loc=locabs+nzmax*((it-2)+nt*(jg-1))
               loca=loc+nzmax
               do jz=1,nsigz
                  a(jz+loca)=a(jz+loca)+a(jz+loc)
               enddo
            enddo
         endif
      enddo
   endif
   if ((nnina.gt.0.and.nina.gt.0).or.&
     (nnina.eq.0.and.jres.eq.0)) then
      ! convert njoy cross sections to cpm resonance integrals
      do jg=1,nrg
         amda=glamt(jnina+nnina*(jg-1))
         sigp=sigpt(nmat)
         do ix=1,2
            do it=1,jtemp
               loc=locs(ix)+nzmax*(it-1+nt*(jg-1))
               call cpmize(nsigz,sigz,a(loc),amda,sigp)
            enddo
         enddo
      enddo
      ! write out results
      if (nwfa.eq.0) call sf3out(xid,nmat,jtemp,nsigz,nt,jfiss,&
        tempt,sigz,flxr,locs,b,a)
      if (nwfa.gt.0) call sf3out(xid,nmat,jtemp,nsigz,nt,jfiss,&
        tempt,sigz,flxr,locs,fa,a)
      indf(nmat)=jfiss
   endif
  470 continue
   if (nmat.eq.newmat) go to 500
   if (irmat.eq.0) go to 110
   if (math.ne.matd) call skiprz(ngendf,-1)
   go to 110

   !--sfile3 is finished.
  500 continue
  510 continue
   deallocate(nszt)
   deallocate(norft)
   deallocate(sigpt)
   deallocate(glamt)
   write(nscr) eflag
   return
   end subroutine sfile3

   subroutine cpmize(nsigz,sigz,sign,amda,sigp)
   !--------------------------------------------------------------------
   ! Convert NJOY shielded cross sections to CPM resonance integrals
   ! using CPM values for intermediate-resonance lambda.
   !--------------------------------------------------------------------
   ! externals
   integer::nsigz
   real(kr)::sigz(*),sign(*),amda,sigp
   ! internals
   integer::iz,jl,jh,jz,kz
   real(kr)::sm,test,sn,temp,st
   real(kr)::sig(10)
   real(kr),parameter::zero=0

   do iz=1,nsigz
      sig(iz)=sign(iz)
   enddo
   do iz=1,nsigz
      sm=sigz(iz)
      test=1001
      if (sm.le.test) then
         sign(iz)=0
         sn=sm-amda*sigp
         if (sn.lt.zero) sn=0
         jl=iz-1
         if (jl.lt.1) jl=1
         jh=iz+1
         if (jh.gt.nsigz) jh=nsigz
         do jz=jl,jh
            temp=1
            sm=sigz(jz)
            do kz=jl,jh
               st=sigz(kz)
               if (kz.ne.jz) then
                  if (sn.ge.sigz(1)) then
                     temp=temp*log(sn/st)/log(sm/st)
                  else
                     temp=temp*(sn**2-st**2)/(sm**2-st**2)
                  endif
               endif
            enddo
            sign(iz)=sign(iz)+temp*sig(jz)
         enddo
      endif
   enddo
   do iz=1,nsigz
      sm=sigz(iz)
      sign(iz)=sm*sign(iz)/(sm+sign(iz))
   enddo
   return
   end subroutine cpmize

   subroutine sf3out(xid,nmat,ntemp,nsigz,nt,ifiss,temp,sigz,flxr,&
     locs,b,a)
   !--------------------------------------------------------------------
   ! Write subfile 3.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   ! externals
   integer::nmat,ntemp,nsigz,nt,ifiss,locs(*)
   real(kr)::xid,temp(*),sigz(*),flxr(*),b(*),a(*)
   ! internals
   integer::ncol1,lim,ident,n,i,jgref,m,ig,it,locf,locr,jz,jg
   integer::il,loc,nb,loca,ib
   real(kr)::c(10)
   character(4)::undl4='----'
   character(4)::undl1='-'
   integer::nmax=5
   integer::ncol=8

   ncol1=ncol+1
   lim=nsigz
   if (nsigz.gt.ncol) lim=ncol
   ident=idnt(nmat)
   n=nsigz
   if (n.gt.nmax) n=nmax
   write(nscr) xid
   write(nscr) ntemp,nsigz
   write(nscr) (temp(i),i=1,ntemp)
   write(nscr) (sigz(i),i=1,nsigz)
   write(nsyso,'(/'' xid '',f8.2/'' ---'')') xid
   write(nsyso,'(/&
     &'' ntemp '',i2,''  nsigp '',i2/&
     &'' ----- '',2x,''  -----''//&
     &'' temperatures (down)''/&
     &'' -------------------'')')&
     ntemp,nsigz
   write(nsyso,'(1x,1p,6e11.3)') (temp(i),i=1,ntemp)
   write(nsyso,'(/&
     &'' potential cross sections (across)''/&
     &'' ---------------------------------'')')
   write(nsyso,'(10x,1p,9e13.5)') (sigz(i),i=1,nsigz)
   jgref=ngnd-ngref+1
   m=n
   if (m.lt.4) m=4
   write(nsyso,'(/'' group'',6x,&
     &''flux per unit lethargy normalized at group '',i2,&
     &'' for nuclide '',i6/&
     &'' ----------'',6(3a4,a1))')&
        jgref,ident,(undl4,undl4,undl4,undl1,i=1,m)
   do ig=1,nrg
      do it=1,ntemp
         locf=nzmax*((it-1)+nt*(ig-1))
         locr=nzmax*(it-1)
         do jz=1,nsigz
            c(jz)=flux(jz+locf)*rleth(ig)*reflth/flxr(jz+locr)
         enddo
         jg=ig+nfg
         if (it.eq.1) write(nsyso,'(i5,5x,1p,9e13.5)')&
           jg,(c(i),i=1,lim)
         if (it.eq.1.and.nsigz.gt.ncol)&
           write(nsyso,'(10x,1p,9e13.5)') (c(i),i=ncol1,nsigz)
         if (it.gt.1) write(nsyso,'(10x,1p,9e13.5)')&
           (c(i),i=1,nsigz)
      enddo
   enddo
   do il=1,2
      loc=locs(il)
      if (iprint.eq.1.and.il.eq.1)&
        write(nsyso,&
        '(/'' group'',6x,''absorption for nuclide '',i6/&
        &'' ----------'',5(3a4,a1))')&
        ident,(undl4,undl4,undl4,undl1,i=1,n)
      if (il.ne.2.or.ifiss.ge.2) then
         if (iprint.eq.1.and.il.eq.2)&
           write(nsyso,&
           '(/'' group'',6x,''fission for nuclide '',i6/&
           &'' ----------'',5(3a4,a1))')&
           ident,(undl4,undl4,undl4,undl1,i=1,n)
         do ig=1,nrg
            nb=0
            jg=ig+nfg
            do it=1,ntemp
               loca=loc-1+nzmax*((it-1)+nt*(ig-1))
               do jz=1,nsigz
                  nb=nb+1
                  b(nb)=a(jz+loca)
                  c(jz)=b(nb)
               enddo
               if (iprint.ne.0) then
                  if (il.ne.2.or.ifiss.ne.1) then
                     if (it.eq.1) write(nsyso,'(/i5,5x,1p,9e13.5)')&
                        jg,(c(i),i=1,lim)
                     if (it.eq.1.and.nsigz.gt.ncol)&
                       write(nsyso,'(10x,1p,9e13.5)')&
                         (c(i),i=ncol1,nsigz)
                     if (it.gt.1) write(nsyso,'(10x,1p,9e13.5)')&
                       (c(i),i=1,lim)
                     if (it.gt.1.and.nsigz.gt.ncol)&
                       write(nsyso,'(10x,1p,9e13.5)')&
                       (c(i),i=ncol1,nsigz)
                  endif
               endif
            enddo
            write(nscr) (b(ib),ib=1,nb)
         enddo
      endif
   enddo
   return
   end subroutine sf3out

   subroutine sfile4(jtemp,nt)
   !--------------------------------------------------------------------
   ! Process cross sections.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use endf ! provides endf routines and variables
   use util ! provides openz,error,skiprz,mess
   ! externals
   integer::jtemp,nt
   ! internals
   integer::il,lim,ngnd5,nb,nw,nmat,matd,ident,ip1opt,mti,mtc
   integer::ntemp,isg,mtabs1,mtabs2,i318,i618,mttot,mtel,mtex
   integer::i,isof,nina,mf3mt2,jfiss,nth1,nth,if6,nz,ia,ib
   integer::lone,ltwo,loc,k,ja,nl,ntw,l,idone,iz,jic,ng2,ig2lo
   integer::ig,igc,jg,loca,locf,jgc,jg2,jg2c,j,ig2,ig2c
   integer::max,min,index,lt,mult,nskip,in,ngndsq,nwloc,iloc
   integer::locab0,locsf0,locchi,locnus,locsg,locxtr,locxs
   real(kr)::sgref,cnorm,dnorm,temp
   character(60)::strng
   real(kr),dimension(:),allocatable::abs1,abs2,sn2n
   real(kr),dimension(:),allocatable::tot,xtr
   real(kr),dimension(:),allocatable::cspc
   integer,dimension(:),allocatable::l1,l2
   integer,dimension(:),allocatable::l1e,l2e
   real(kr),dimension(:),allocatable::snu
   real(kr),dimension(:),allocatable::a
   real(kr),parameter::dilinf=1.e10_kr
   real(kr),parameter::eflag=1.e35_kr
   integer::lz=6
   real(kr),parameter::zero=0

   !--allocate storage.
   if (if4.eq.0) go to 620
   call openz(-mscr,1)
   call repoz(-mscr)
   call openz(-kscr,1)
   call repoz(-kscr)
   if (nres.eq.0) write(nsyso,'(/)')
   if (nres.gt.0) write(nsyso,'(//)')
   write(nsyso,&
     '(/'' ***subfile 4***cross sections'')')
   allocate(fiss(newmat))
   if (iu.eq.1) allocate(uff(nfsg))
   if (ipu.eq.1) allocate(puff(nfsg))
   allocate(snu(ngnd))
   allocate(l1(ngnd))
   allocate(l2(ngnd))
   allocate(l1e(ngnd))
   allocate(l2e(ngnd))
   ngndsq=ngnd*ngnd
   allocate(elas(ngndsq))
   allocate(cspc(ngn))
   allocate(abs1(ngnd))
   allocate(abs2(ngnd))
   allocate(sn2n(ngnd))
   allocate(tot(ngnd))
   allocate(xtr(ngnd))
   nw=nt+1
   if (nres.eq.0) allocate(tempt(nw))
   nwloc=ngnd*6+ngndsq
   allocate(a(nwloc))
   iloc=1
   locab0=iloc
   locsf0=locab0+ngnd
   locchi=locsf0+ngnd
   locnus=locchi+ngnd
   locsg=locnus+ngnd
   locxtr=locsg+ngnd
   locxs=locxtr+ngnd
   nfiss(1)=0
   nfiss(2)=0
   il=1
   lim=locxs+ngndsq-1
   ngnd5=ngnd*5

   !--search for desired materials.
   call repoz(ngendf)
   call tpidio(ngendf,0,0,scr,nb,nw)
   nmat=0
  120 continue
   nmat=nmat+1
   matd=mat(nmat)
   ident=idnt(nmat)
   ip1opt=ip1tab(nmat)
   mti=mtit(nmat)
   mtc=mtct(nmat)
   ntemp=ntmp(nmat)
   sgref=s0rf(nmat)
   isg=0
   if (abs(sgref-dilinf).ge.sgref/100) isg=1
   mtabs1=0
   mtabs2=0
   i318=0
   i618=0
   mttot=0
   mtel=0
   mtex=0
   do i=1,ngnd
      abs2(i)=0
      sn2n(i)=0
      a(i-1+locxtr)=0
      l1e(i)=ngnd
      l2e(i)=1
      l1(i)=ngnd
      l2(i)=1
   enddo
   isof=0
   if (ident.eq.iu235) isof=1
   if (ident.eq.ipu239) isof=2
   nina=ninat(nmat)
   mf3mt2=0
   cnorm=0
   dnorm=0
   jfiss=0
   jtemp=-1
   nth1=0
   nth=0
   if6=0
   do i=1,nwloc
      a(i-1+iloc)=0
   enddo
   if (nina.lt.2) go to 140
   ! read in absorption
   nz=ngnd+2
   read(nsysi,*) (scr(i),i=1,nz)
   awrt(nmat)=scr(1)
   tempt(1)=scr(2)
   do i=3,nz
      a(i-3+locab0)=scr(i)
   enddo
   if (nina.ge.3) then
      read(nsysi,*) (a(locsf0+i-1),i=1,ngnd)
      read(nsysi,*) (a(locnus+i-1),i=1,ngnd)
      read(nsysi,*) (a(locxtr+i-1),i=1,ngnd)
      do ia=1,ngnd
         nz=3
         read(nsysi,*) (scr(i),i=1,nz)
         ib=nint(scr(1))
         if (ib.ne.0) then
            l1(ia)=nint(scr(2))
            l2(ia)=nint(scr(3))
            lone=nint(scr(2))
            ltwo=nint(scr(3))
            nz=ltwo-lone+1
            read(nsysi,*) (scr(i),i=1,nz)
            loc=locxs+ib-1-ngnd
            k=0
            do ja=lone,ltwo
               k=k+1
               a(ngnd*ja+loc)=scr(k)
            enddo
         else
            l1(ia)=ngnd
            l2(ia)=1
         endif
      enddo
      jtemp=1
   else
      jtemp=-1
      matd=-1
   endif
   go to 500
  140 continue
   call contio(ngendf,0,0,scr,nb,nw)
   if (math.eq.matd) go to 150
   if (math.eq.0) go to 140
   if (math.gt.0) go to 145
   call repoz(ngendf)
   call tpidio(ngendf,0,0,scr,nb,nw)
   go to 140
  145 continue
   call tomend(ngendf,0,0,scr)
   go to 140

   !--search for required reactions.
  150 continue
   nl=l1h
   nz=l2h
   ntw=n2h
   if (mfh.eq.3.and.mth.ge.2.and.mf3mt2.le.0) then
      ! invert flux arrays
      mf3mt2=1
      do i=1,ngnd
         if (cflux(i).ne.zero) cflux(i)=1/cflux(i)
         if (cflux(ngnd+i).ne.zero) cflux(ngnd+i)=1/cflux(ngnd+i)
         loc=i
         if (flux(loc).ne.zero) flux(loc)=1/flux(loc)
         loc=loc+ngnd
         if (flux(loc).ne.zero) flux(loc)=1/flux(loc)
      enddo
   endif
  160 continue
   l=1
   call listio(ngendf,0,0,scr(l),nb,nw)
   if (mfh.ne.1) go to 170
   temp=c1h
   ngn=l1h
   if (ngn.ne.ngnd.and.ngn.ne.ngnf)&
     call error('sfile4','storage exceeded.',' ')
   idone=0
   i=0
   do while (i.lt.nz.and.idone.eq.0)
      i=i+1
      iz=i
      if (abs(sgref-scr(l+5+ntw+i)).lt.sgref/100) idone=1
   enddo
   jtemp=jtemp+1
   jic=0
   if (jtemp.eq.0) go to 195
   if (ntemp.gt.0.and.jtemp.gt.ntemp) go to 320
   tempt(jtemp)=temp
   go to 195
  170 continue
   ng2=l1h
   ig2lo=l2h
   ig=n2h
   igc=1
   if (ig.ne.0) igc=icgrp(ig)
   ! check for presence of sometimes temperature-dependent data
   if (jtemp.eq.0) go to 185
   if (mth.lt.102.or.mth.gt.150) go to 174
   if (mth.ne.102) go to 172
   if (mtabs1.gt.0) go to 185
   mtabs1=1
   do i=1,ngnd
      abs1(i)=0
   enddo
   go to 174
  172 continue
   if (mtabs2.gt.0) go to 185
   mtabs2=1
   do i=1,ngnd
   abs2(i)=0
   enddo
  174 continue
   if (mth.ne.2) go to 179
   if (mfh.ne.6) go to 185
   if (mtel.gt.0) go to 177
   mtel=1
   do i=1,ngnd
      l1e(i)=ngnd
      l2e(i)=1
   enddo
   do i=1,ngndsq
      elas(i)=0
   enddo
  177 continue
   if (mtex.gt.0) go to 185
   mtex=1
   do i=1,ngnd
      xtr(i)=0
   enddo
  179 continue
   if (mth.ne.1) go to 185
   if (mttot.gt.0) go to 185
   mttot=1
   do i=1,ngnd
      tot(i)=0
   enddo
  185 continue
   do while (nb.ne.0)
      l=l+nw
      call moreio(ngendf,0,0,scr(l),nb,nw)
      if ((l+nw-1).gt.nwscr)&
        call error('sfile4','storage exceeded.',' ')
   enddo
   if (igc.eq.0.or.igc.gt.ngnd) go to 218
   jg=ngnd-igc+1
   l=1
   if (mth.eq.455.and.jtemp.gt.0) go to 240
   if (mfh.eq.3) go to 205
   if (mfh.eq.6) go to 255
  195 continue
   call tosend(ngendf,0,0,scr(l))
  200 continue
   call contio(ngendf,0,0,scr(l),nb,nw)
   if (math.eq.0) go to 330
   if (math.eq.-1) go to 350
   if (math.ne.matd) go to 350
   if (mfh.eq.0.or.mth.eq.0) go to 200
   go to 150

   !--file 3 cross sections.
  205 continue
   if (jtemp.eq.0) go to 220
   if (mth.lt.102.or.mth.gt.150) go to 214
   ! absorption
   loca=l+lz+nl*nz
   if (isg.gt.0.and.nz.ge.iz) loca=l+lz+nl*(iz-1+nz)
   locf=l+lz
   if (isg.gt.0.and.nz.ge.iz) locf=l+lz+nl*(iz-1)
   jgc=jg
   if (isg.gt.0.and.nz.ge.iz) jgc=jg+ngnd
   if (mth.ne.102) then
      abs2(jg)=abs2(jg)+scr(loca)*scr(locf)*cflux(jgc)
   else
      abs1(jg)=abs1(jg)+scr(loca)*scr(locf)*cflux(jgc)
   endif
   if (ig.lt.ngn) go to 160
   if (nina.gt.0) call tomend(ngendf,0,0,scr)
   if (nina.gt.0) go to 400
   go to 200
  214 continue
   if (mth.eq.1) go to 220
   if (nina.gt.0) go to 195
   if (mth.eq.16) go to 220
   if ((mth.lt.18.or.mth.gt.21).and.mth.ne.38) go to 220
   if (mth.eq.18) i318=1
   if (mth.eq.18) go to 216
   if (i318.eq.0) go to 216
   go to 195
   ! fission
  216 continue
   loc=locsf0+jg-1
   loca=l+lz+nl*nz
   locf=l+lz
   a(loc)=a(loc)+scr(loca)*scr(locf)*cflux(jg)
   loc=locab0+jg-1
   if (isg.gt.0.and.nz.ge.iz) locf=l+lz+nl*(iz-1)
   jgc=jg
   if (isg.gt.0.and.nz.ge.iz) jgc=jg+ngnd
   a(loc)=a(loc)+scr(loca)*scr(locf)*cflux(jgc)
  218 continue
   if (ig.lt.ngn) go to 160
   go to 200
  220 continue
   if (mth.ne.1.and.jtemp.gt.0) go to 195
   if (mth.ne.1.and.jtemp.eq.0) go to 236
   if (ig.eq.1) then
      do i=1,ngnd
         cflux(i)=0
         cflux(i+ngnd)=0
         flux(i)=0
         flux(i+ngnd)=0
      enddo
      mf3mt2=0
   endif
   ! p1 flux
   ! for infinite dilution
   loc=jg
   loca=l+lz
   if (nl.gt.1) loca=loca+1
   flux(loc)=flux(loc)+scr(loca)
   if (isg.ne.0.and.nz.ge.iz) then
      ! for reference sigma zero
      loc=jg+ngnd
      loca=l+lz+nl*(iz-1)
      if (nl.gt.1) loca=loca+1
      flux(loc)=flux(loc)+scr(loca)
   endif
   ! p0 flux
   ! for infinite dilution
   loca=l+lz
   cflux(jg)=cflux(jg)+scr(loca)
   ! for reference sigma zero
   if (isg.eq.0.or.nz.lt.iz) go to 250
   cflux(jg+ngnd)=cflux(jg+ngnd)+scr(l+lz+nl*(iz-1))
   go to 250
   ! n2n
  236 continue
   if (jtemp.gt.0) go to 195
   if (mth.ne.16) go to 238
   loca=l+lz+nl*nz
   locf=l+lz
   sn2n(jg)=sn2n(jg)+scr(loca)*scr(locf)*cflux(jg)
   go to 250
  238 continue
   if (mth.ne.mti) go to 195
   if (igc.gt.nth1.and.scr(l+lz+nz*nl).ne.0.) nth1=igc
   nth=ngnd-nth1+1
   go to 250
  240 continue
   if (mth.ne.455) go to 195
   if (isof.eq.0) go to 195
   jfiss=1
   if (mfh.eq.5) go to 244
   ! delayed neutrons
   loc=locnus+jg-1
   loca=l+lz
   locf=l+lz
   a(loc)=a(loc)+scr(loca+1)*scr(loca+2)*scr(locf)*cflux(jg)
   dnorm=dnorm+scr(loca)*scr(loca+1)*scr(loca+2)
   if (ig.lt.ngn) go to 160
   go to 200
  244 continue
   if (isof.eq.0) go to 195
   do i=2,ng2
      jg2=ngn-ig2lo-i+3
      jg2c=icgrp(jg2)
      if (jg2c.ne.0.and.jg2c.le.ngnd) then
         loc=locchi+jg2c-1
         loca=l+lz-1+nl*(i-1)
         do j=1,nl
            a(loc)=scr(j+loca)*dnorm+a(loc)
            cnorm=scr(j+loca)*dnorm+cnorm
         enddo
      endif
   enddo
  250 continue
   if (ig.lt.ngn) go to 160
   go to 200

   !--file 6 matrices.
  255 continue
   if (jtemp.eq.0) go to 300
   if ((mth.lt.18.or.mth.gt.21).and.mth.ne.38) go to 270
   if (mth.eq.18) i618=1
   if (mth.eq.18) go to 260
   if (i618.eq.0) go to 260
   go to 195
   ! fission
  260 continue
   jfiss=1
   if (ig.ne.0) then
      locf=l+lz
      do i=2,ng2
         loca=l+lz+(il-1)+nl*nz*(i-1)
         loc=locnus+jg-1
         a(loc)=a(loc)+scr(loca)*scr(locf)*cflux(jg)
         if (isof.ne.0) then
            if (ig2lo.ne.0) then
               ! matrix part
               ig2=ig2lo+i-2
               ig2c=icgrp(ig2)
               if (ig2c.ne.0.and.ig2c.le.ngnd) then
                  jg2c=ngnd-ig2c+1
                  loc=locchi+jg2c-1
                  a(loc)=a(loc)+scr(loca)*scr(locf)
                  cnorm=cnorm+scr(loca)*scr(locf)
               endif
            else
               ! spectrum part
               do j=1,ngn
                  ig2=ig2lo+j-1
                  ig2c=icgrp(ig2)
                  if (ig2c.ne.0.and.ig2c.le.ngnd) then
                     jg2c=ngnd-ig2c+1
                     loc=locchi+jg2c+1
                     a(loc)=a(loc)+scr(loca)*scr(locf)*cspc(j)
                     cnorm=cnorm+scr(loca)*scr(locf)*cspc(j)
                  endif
               enddo
            endif
      endif
      enddo
   else
      ! save constant spectrum
      do i=1,ng2
         cspc(ig2lo+i-1)=scr(l+lz+i-1)
      enddo
   endif
   go to 295
   ! add temperature-dependent reactions to matrix.
  270 continue
   if (mth.ne.2.and.mth.ne.mti.and.mth.ne.mtc) go to 295
   if (mth.ne.2) jic=1
   max=0
   min=ngnd
   locf=l+lz
   if (isg.gt.0.and.nz.ge.iz) locf=l+lz+nl*(iz-1)
   jgc=jg
   if (isg.gt.0.and.nz.ge.iz) jgc=jg+ngnd
   if (mth.eq.2.and.jg.ge.nth) go to 295
   do i=2,ng2
      ig2=ig2lo+i-2
      ig2c=icgrp(ig2)
      if (ig2c.ne.0.and.ig2c.le.ngnd) then
         jg2c=ngnd-ig2c+1
         loca=l+lz+(il-1)+nl*nz*(i-1)
         if (isg.gt.0.and.nz.ge.iz) loca=loca+nl*(iz-1)
         if (scr(loca).ne.zero) then
            if (mth.ne.2) then
               loc=locxs+jg-1+ngnd*(jg2c-1)
               if ((jg2c-jg).ge.-1.or.scr(loca).ge.zero) then
                  a(loc)=a(loc)+scr(loca)*scr(l+lz)*cflux(jg)
                  loc=locxtr+jg-1
                  a(loc)=a(loc)+scr(loca)*scr(l+lz)*cflux(jg)
                  if (jg2c.gt.max) max=jg2c
                  if (jg2c.lt.min) min=jg2c
               endif
               if (nl.ne.1) then
                  if (jg2c.ge.nth) then
                     if (igc.le.nth1) then
                        index=jg2c
                        loc=locxtr+jg2c-1
                        a(loc)=a(loc)&
                          -scr(loca+1)*scr(l+lz+1)*flux(index)
                        if (ip1opt.ne.0) then
                           loc=locxs+jg2c-1+ngnd*(jg2c-1)
                           a(loc)=a(loc)&
                             -scr(loca+1)*scr(l+lz+1)*flux(index)
                        endif
                     endif
                  endif
               endif
            else
               loc=jg+ngnd*(jg2c-1)
               elas(loc)=elas(loc)+scr(loca)*scr(locf)*cflux(jgc)
               if (jg2c.gt.max) max=jg2c
               if (jg2c.lt.min) min=jg2c
               xtr(jg)=xtr(jg)+scr(loca)*scr(locf)*cflux(jgc)
               if (nl.ne.1) then
                  index=jg2c
                  if (isg.gt.0.and.nz.ge.iz) index=index+ngnd
                  xtr(jg2c)=xtr(jg2c)&
                    -scr(loca+1)*scr(locf+1)*flux(index)
                  if (ip1opt.ne.0) then
                     loc=jg2c+ngnd*(jg2c-1)
                     elas(loc)=elas(loc)&
                       -scr(loca+1)*scr(locf+1)*flux(index)
                  endif
               endif
            endif
         endif
      endif
   enddo
   lt=l1(jg)
   if (min.lt.lt.and.mth.ne.2) l1(jg)=min
   lt=l1e(jg)
   if (min.lt.lt.and.mth.eq.2) l1e(jg)=min
   lt=l2(jg)
   if (max.gt.lt.and.mth.ne.2) l2(jg)=max
   lt=l2e(jg)
   if (max.gt.lt.and.mth.eq.2) l2e(jg)=max
  295 continue
   if (ig.lt.ngn) go to 160
   go to 200

   !--store non-temperature dependent matrix .
  300 continue
   if (jtemp.gt.0) go to 205
   if (mfh.ne.6) go to 195
   if (if6.le.0) then
      if (nth1.eq.0) nth1=ngnd-nfg-nrg
      nth=ngnd-nth1+1
      if6=1
   endif
   if (mth.ge.18.and.mth.le.21) go to 195
   if (mth.eq.38) go to 195
   if (mth.eq.2) go to 195
   if (mth.ge.201.and.mth.le.255) go to 195
   mult=1
   if (mth.eq.16) mult=2
   locf=l+lz
   max=0
   min=ngnd
   do i=2,ng2
      ig2=ig2lo+i-2
      ig2c=icgrp(ig2)
      if (ig2c.ne.0.and.ig2c.le.ngnd) then
         jg2c=ngnd-ig2c+1
         loca=l+lz+(il-1)+nl*nz*(i-1)
         loc=locxs+jg-1+ngnd*(jg2c-1)
         a(loc)=a(loc)+scr(loca)*scr(locf)*cflux(jg)
         loc=locxtr+jg-1
         a(loc)=a(loc)+(scr(loca)*mult)*scr(locf)*cflux(jg)
         if (jg2c.gt.max) max=jg2c
         if (jg2c.lt.min) min=jg2c
         if (nl.ne.1) then
            loc=locxtr+jg2c-1
            index=jg2c
            a(loc)=a(loc)-scr(loca+1)*scr(locf+1)*flux(index)
            if (ip1opt.ne.0) then
               loc=locxs+jg2c-1+ngnd*(jg2c-1)
               a(loc)=a(loc)-scr(loca+1)*scr(locf+1)*flux(index)
            endif
         endif
      endif
   enddo
   lt=l1(jg)
   if (min.lt.lt) l1(jg)=min
   lt=l2(jg)
   if (max.gt.lt) l2(jg)=max
   if (ig.lt.ngn) go to 160
   go to 200

   !--write non-temperature dependent matrix on scratch tape.
  320 continue
   jtemp=jtemp-1
  330 continue
   if (jtemp.gt.0) go to 360
   call repoz(-kscr)
   write(kscr) (l1(i),i=1,ngnd)
   write(kscr) (l2(i),i=1,ngnd)
   write(kscr) (a(locxtr+i-1),i=1,ngnd)
   write(kscr) (a(i),i=locxs,lim)
   nskip=-3
   call skiprz(ngendf,nskip)
   call findf(matd,1,451,ngendf)
   if (nth1.gt.0.and.nth.gt.0) go to 140
   if (nth1.eq.0) nth1=ngnd-nfg-nrg
   nth=ngnd-nth1+1
   go to 140

   !--desired data is loaded.
  350 continue
   matd=-1
   if (math.ne.0) call skiprz(ngendf,-3)
   go to 500
  360 continue
   if (jic.gt.0) go to 365
   jtemp=jtemp-1
   write(strng,'(''use only '',i2,'' temperatures for mat '',i4)')&
     jtemp,matd
   call mess('sfile4',strng,'mti missing from higher temperatures')
   go to 350
  365 continue
   if (jtemp.eq.1) then
      ! calculate nu at lowest temperature
      do i=1,ngnd
         if (a(locsf0+i-1).eq.zero) snu(i)=0
         if (a(locsf0+i-1).ne.zero) snu(i)&
           =a(locnus+i-1)/a(locsf0+i-1)
      enddo
   else
      ! calculate nu*sigma f at higher temperatures
      do i=1,ngnd
         a(i-1+locnus)=snu(i)*a(i-1+locsf0)
      enddo
   endif
   ! calculate total absorption
  400 continue
   do i=1,ngnd
      a(i-1+locab0)=a(i-1+locab0)+abs1(i)+abs2(i)-sn2n(i)
   enddo
   ! correct scattering matrix
   in=nth-1
   do i=1,in
      if (l1e(i).lt.l1(i)) l1(i)=l1e(i)
      if (l2e(i).gt.l2(i)) l2(i)=l2e(i)
   enddo
   do i=1,ngndsq
      a(i-1+locxs)=a(i-1+locxs)+elas(i)
   enddo
   ! correct transport
   do i=1,ngnd
      a(i-1+locxtr)=a(i-1+locxtr)+xtr(i)
   enddo
  500 continue
   call sf4out(matd,nmat,ntemp,jfiss,nina,&
     locab0,locsf0,locnus,locxtr,locxs,l1,l2,a)
   if (nina.ne.3) go to 505
   matd=-1
   nina=0
   go to 500
  505 continue
   if (matd.lt.0) go to 540
   mtabs1=0
   mtabs2=0
   mttot=0
   mtel=0
   mtex=0
   i318=0
   i618=0
   call repoz(-kscr)
   read(kscr) (l1(i),i=1,ngnd)
   read(kscr) (l2(i),i=1,ngnd)
   read(kscr) (a(locxtr+i-1),i=1,ngnd)
   read(kscr) (a(i),i=locxs,lim)
   fiss(nmat)=jfiss
   if (isof.ne.0.and.jtemp.le.1) then
      ! normalize chi
      if (cnorm.ne.zero) cnorm=1/cnorm
      do i=1,nfsg
         loc=locchi+i-1
         if (cnorm.ne.zero) a(loc)=a(loc)*cnorm
         if (nfiss(isof).ne.0.and.a(loc).ne.zero)&
           nfiss(isof)=nfiss(isof)+1
         if (nfiss(isof).eq.0.and.a(loc).ne.zero)&
           nfiss(isof)=1
         if (isof.eq.1) uff(i)=a(loc)
         if (isof.eq.2) puff(i)=a(loc)
      enddo
   endif
   do i=1,ngnd5
      a(i-1+iloc)=0
   enddo
   go to 200
  540 continue
   if (ntemp.eq.0) ntmp(nmat)=jtemp
   if (nmat.eq.newmat) go to 600
   go to 120

   !--sfile4 is finished.
  600 continue
   deallocate(tempt)
   deallocate(flux)
   ! re-invert coarse group p0 flux for infinite dilution
   do i=1,ngnd
      if (cflux(i).ne.zero) cflux(i)=1/cflux(i)
   enddo
  620 continue
   write(nscr) eflag
   return
   end subroutine sfile4

   subroutine sf4out(matd,nmat,ntemp,ifiss,nina,&
     locab0,locsf0,locnus,locxtr,locxs,l1,l2,a)
   !--------------------------------------------------------------------
   ! Write subfile 4 on nscr.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use util ! provides repoz
   ! externals
   integer::matd,nmat,ntemp,ifiss,nina
   integer::locab0,locsf0,locnus,locxtr,locxs,l1(*),l2(*)
   real(kr)::a(*)
   ! internals
   integer::ncol1,ident,jtemp,lim,i,ktemp,nc,ic,ib,i1,ig
   integer::ia,lone,ltwo,nb,loc
   real(kr)::flag
   real(kr)::b(70)
   integer::ncol=5
   integer::matl=0
   integer::ifissl=0
   integer::ninal=0
   character(4)::hend='end '
   integer::izero=0
   real(kr),parameter::zero=0

   ncol1=ncol-1
   if (nina.eq.2) go to 120
   if (matd.eq.matl) go to 160
   if (matl.eq.0) go to 150
   ! material has been completed.  copy mscr to nscr.
   read(hend,'(a4)') flag
   write(mscr) flag
   ident=idnt(nmat)
   if (iprint.eq.1.and.nmat.gt.1) write(nsyso,'(//)')
   if (iprint.eq.1) write(nsyso,'(80i1)') ident
   call repoz(-mscr)
   jtemp=1
   write(nscr) (tempt(i),i=1,ntemp)
   if (iprint.ne.0) then
      write(nsyso,'(/&
        &'' temperatures''/&
        &'' ------------'')')
      write(nsyso,'(1p,6e12.4)') (tempt(i),i=1,ntemp)
   endif
  100 continue
   ktemp=nint(tempt(jtemp))
   nc=4
   ic=1
   if (iprint.eq.1) write(nsyso,'(/&
     &'' absorption for nuclide '',i6,'' at '',i4,'' deg k''/&
     &'' -------------------------------------------'')')&
     ident,ktemp
  105 continue
   read(mscr) (b(i),i=1,ngnd)
   write(nscr) (b(i),i=1,ngnd)
   if (iprint.eq.1.and.ic.eq.1) write(nsyso,'(1p,6e12.4)')&
     (b(i),i=1,ngnd)
   if (iprint.eq.1.and.ifissl.eq.1.and.ic.gt.1.and.ic.lt.nc)&
     write(nsyso,'(1p,6e12.4)') (b(i),i=1,ngnd)
   if (iprint.eq.1.and.ic.eq.nc) write(nsyso,'(1p,6e12.4)')&
     (b(i),i=1,ngnd)
   ic=ic+1
   if (ninal.eq.1) go to 135
   if (ic.gt.nc) go to 125
   if (iprint.eq.0) go to 115
   if (ifissl.eq.0.and.ic.gt.1.and.ic.lt.nc) go to 115
   if (ic.eq.1) go to 115
   if (ic.eq.2) then
      write(nsyso,'(/&
        &'' fission for nuclide '',i6,'' at '',i4,'' deg k''/&
        &'' ----------------------------------------'')')&
        ident,ktemp
   else if (ic.eq.3) then
      write(nsyso,'(/&
        &'' nu*fission for nuclide '',i6,'' at '',i4,'' deg k''/&
        &'' -------------------------------------------'')')&
        ident,ktemp
   else if (ic.eq.4) then
      write(nsyso,'(/&
        &'' transport corrected total for nuclide '',i6,'' at '',&
        &i4,'' deg k''/&
        &'' ----------------------------------------------'',&
        &''------------'')') ident,ktemp
   endif
  115 continue
   ib=0
   go to 105
  120 continue
   ident=idnt(nmat)
   ktemp=0
   write(nscr) tempt(1)
   write(nscr) (a(locab0-1+i),i=1,ngnd)
   if (iprint.eq.1) then
      if (nmat.gt.1) write(nsyso,'(//)')
      write(nsyso,'(/&
        &'' nuclide '',i6/&
        &'' -------'')')&
        ident
      write(nsyso,'(/&
        &'' temperatures''/&
        &'' ------------'')')
      write(nsyso,'(1p,6e12.4)') tempt(1)
      write(nsyso,'(/&
        &'' absorption for nuclide '',i6,'' at '',i4,'' deg k''/&
        &'' -------------------------------------------'')')&
        ident,ktemp
      write(nsyso,'(1p,6e12.4)') (a(locab0-1+i),i=1,ngnd)
   endif
   ifissl=0
   ninal=1
   go to 200
  125 continue
   if (iprint.eq.1) write(nsyso,'(/&
     &'' p0 scattering matrix for nuclide '',i6,'' at '',i4,&
     &'' deg k''/&
     &'' -----------------------------------------------------''//&
     &'' ia l1 l2    l1+0'',4(10x,''+'',i1)/&
     &1x,7(''----------''),''----'')')&
     ident,ktemp,(i,i=1,ncol1)
   i1=1
   do 130 ig=1,ngnd
   read(mscr) ia,lone,ltwo
   if (ia.eq.0) go to 135
   write(nscr) ia,lone,ltwo
   nb=ltwo-lone+1
   read(mscr) (b(i),i=1,nb)
   write(nscr) (b(i),i=1,nb)
   i1=1
   lim=nb
   if ((lim-i1+1).gt.ncol) lim=ncol
   write(nsyso,'(3i3,5x,1p,5e12.4)')&
     ia,lone,ltwo,(b(i),i=i1,lim)
   do while (lim.ne.nb)
      i1=lim+1
      lim=nb
      if ((lim-i1+1).gt.ncol) lim=i1+ncol-1
      write(nsyso,'(9x,i4,2x,1p,5e12.4)') (b(i),i=i1,lim)
   enddo
  130 continue
   read(mscr) ia
  135 continue
   jtemp=jtemp+1
   if (jtemp.le.ntemp) go to 100
  150 continue
   if (matd.lt.0) go to 200
   call repoz(-mscr)
   matl=matd
   ifissl=ifiss
   ninal=nina
   ! write data on mscr.
  160 continue
   lim=locab0+ngnd-1
   write(mscr) (a(i),i=locab0,lim)
   if (nina.le.0.or.nina.eq.3) then
      lim=locsf0+ngnd-1
      write(mscr) (a(i),i=locsf0,lim)
      lim=locnus+ngnd-1
      write(mscr) (a(i),i=locnus,lim)
      lim=locxtr+ngnd-1
      write(mscr) (a(i),i=locxtr,lim)
      do ia=1,ngnd
         lone=l1(ia)
         ltwo=l2(ia)
         if (lone.ne.ngnd.or.ltwo.ne.1) then
            ib=0
            do i=lone,ltwo
               loc=locxs+ia-1+ngnd*(i-1)
               ib=ib+1
               b(ib)=a(loc)
            enddo
            write(mscr) ia,lone,ltwo
            write(mscr) (b(i),i=1,ib)
         else
            write(mscr) ia,ia,ia
            write(mscr) zero
         endif
      enddo
      write(mscr) izero
   endif

   !--sf4out is finished.
  200 continue
   if (matd.lt.0) matl=0
   return
   end subroutine sf4out

   subroutine sf5out
   !--------------------------------------------------------------------
   ! Write subfile 5 burnup data.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   ! internals
   integer::nfis1,i,index,j
   integer::id(70)
   real(kr),parameter::eflag=1.e35_kr

   write(nsyso,'(/'' ***subfile 5***burnup data'')')

   nfis1=nfis+1
   do i=1,ntis
      id(i)=idat(i)
   enddo
   write(nscr) ntis,(id(i),i=1,ntis)
   do i=1,nfis
      id(i)=idbt(i)
   enddo
   write(nscr) nfis,(id(i),i=1,nfis)
   do i=1,ntis
      index=1+(i-1)*nfis1
      write(nscr) (yld(index+j-1),j=1,nfis1)
   enddo
   write(nscr) eflag
   deallocate(idat)
   deallocate(idbt)
   deallocate(yld)
   return
   end subroutine sf5out

   subroutine sfile6(ntemp)
   !--------------------------------------------------------------------
   ! Process P1 scattering matrix.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use util ! provides repoz,error,closz
   use endf ! provides endf routines and variables
   ! externals
   integer::ntemp
   ! internals
   integer::nther,ngndsq,nw,nb,il,np1,nmat,ip1opt,nina
   integer::matd,mti,mtc,ntemps,mtelas,itd,i,nth1,nth,ip1
   integer::iza,nl,nz,l,ig2lo,ig,igc,jg,ig2,ig2c,jg2c
   integer::loca,loc,lt,l1t,l1tc,l2t,l2tc,nskip,ng2
   real(kr)::temp,flag
   character(60)::strng
   integer::ntp(6)
   real(kr)::tempr(7)
   integer::lz=6
   real(kr),dimension(:),allocatable::tempt
   integer,dimension(:),allocatable::ids
   integer,dimension(:),allocatable::l1,l2
   integer,dimension(:),allocatable::l1e,l2e
   real(kr),dimension(:),allocatable::s
   character(4)::hflag='end '
   real(kr),parameter::eflag=1.e35_kr

   !--allocate storage.
   if (if4.eq.0) go to 410
   write(nsyso,'(//)')
   write(nsyso,&
     '(/'' ***subfile 6***p1 matrices'')')
   if (nnina.eq.0) go to 400
   nther=ngnd-nfg-nrg
   allocate(ids(newmat))
   ngndsq=ngnd*ngnd
   allocate(s(ngndsq))
   nw=newmat*6
   allocate(tempt(nw))
   allocate(l1(ngnd))
   allocate(l2(ngnd))
   allocate(l1e(ngnd))
   allocate(l2e(ngnd))

   !--search for desired materials.
   call repoz(ngendf)
   call tpidio(ngendf,0,0,scr,nb,nw)
   il=2
   np1=0
   nmat=0
  100 continue
   nmat=nmat+1
   ip1opt=ip1tab(nmat)
   if (ip1opt.gt.0) go to 320
   nina=ninat(nmat)
   if (nina.gt.0) go to 320
   matd=mat(nmat)
   mti=mtit(nmat)
   mtc=mtct(nmat)
   ntemps=ntmp(nmat)
   mtelas=0
   itd=0
   do i=1,ngndsq
      s(i)=0
   enddo
   do i=1,ngnd
      l1(i)=ngnd
      l2(i)=1
   enddo
   ntemp=-1
   nth1=0
   nth=nther
   ip1=0
   math=1
   do while (math.ne.matd)
      call contio(ngendf,0,0,scr,nb,nw)
      if (math.ne.matd) then
         if (math.lt.0) then
            call repoz(ngendf)
            call tpidio(ngendf,0,0,scr,nb,nw)
         else if (math.gt.0) then
            call tomend(ngendf,0,0,scr)
         endif
      endif
   enddo
   go to 130

   !--search for required reactions.
  125 continue
   call contio(ngendf,0,0,scr,nb,nw)
   if (math.eq.0) go to 200
   if (math.lt.0) go to 300
   if (math.ne.matd) go to 300
   if (mfh.eq.0.or.mth.eq.0) go to 125
   if (mfh.eq.1.and.ntemp.eq.ntemps) go to 300
   if (mfh.eq.1) go to 130
   go to 135
  130 continue
   iza=nint(c1h)
  135 continue
   nl=l1h
   if (nl.lt.il.and.mfh.ne.1) go to 155
   nz=l2h
  140 continue
   l=1
   call listio(ngendf,0,0,scr(l),nb,nw)
   if (mfh.eq.1) then
      temp=c1h
      ntemp=ntemp+1
   endif
   if (mfh.eq.3) go to 160
   if (mfh.ne.6) go to 155
   if (mth.ne.2) go to 148
   if (mtelas.gt.0) go to 148
   mtelas=1
   do i=1,ngnd
      l1e(i)=ngnd
      l2e(i)=1
   enddo
   do i=1,ngndsq
      elas(i)=0
   enddo
  148 continue
   ng2=l1h
   ig2lo=l2h
   ig=n2h
   igc=icgrp(ig)
   do while (nb.ne.0)
      l=l+nw
      call moreio(ngendf,0,0,scr(l),nb,nw)
      if ((l+nw-1).gt.nwscr)&
        call error('sfile6','storage exceeded.',' ')
   enddo
  155 continue
   call tosend(ngendf,0,0,scr)
   go to 125
  160 continue
   if (igc.eq.0.or.igc.gt.ngnd) go to 225
   jg=ngnd-igc+1
   l=1

   !--store non-temperature dependent scattering matrix
   if (ntemp.gt.0) go to 210
   if (mfh.eq.3) go to 230
   if (mth.ge.18.and.mth.le.21) go to 155
   if (mth.eq.38) go to 155
   if (mth.eq.2) go to 155
   if (mth.eq.mti.or.mth.eq.mtc) go to 180
   if (mth.ge.201.and.mth.le.250) go to 155
   do i=2,ng2
      ig2=ig2lo+i-3
      ig2c=icgrp(ig2)
      if (ig2c.ne.0.and.ig2c.le.ngnd) then
         jg2c=ngnd-ig2c+1
         loca=l+lz+(il-1)+nl*nz*(i-1)
         loc=1+jg-1+ngnd*(jg2c-1)
         s(loc)=s(loc)+scr(loca)
      endif
   enddo
   ip1=1
   lt=l1(jg)
   l1t=ngn-ig2lo-ng2+3
   l1tc=icgrp(l1t)
   if (l1tc.lt.lt.and.l1tc.ne.0.and.l1tc.le.ngnd) l1(jg)=l1tc
   lt=l2(jg)
   l2t=l1t+ng2-2
   l2tc=icgrp(l2t)
   if (l2tc.gt.lt.and.l2tc.ne.l1t.and.l2tc.gt.0.and.l2tc.le.ngnd)&
     l2(jg)=l2tc
   if (ig.lt.ngn) go to 140
   go to 125
  180 continue
   if (igc.gt.nth1.and.igc.lt.nth) nth1=igc
   if (igc.lt.nth) go to 140
   go to 155

   !--write non-temperature dependent matrix on scratch tape.
  200 continue
   if (ntemp.gt.0) go to 250
   call repoz(-kscr)
   write(kscr) (l1(i),i=1,ngnd)
   write(kscr) (l2(i),i=1,ngnd)
   write(kscr) (s(i),i=1,ngndsq)
   call repoz(-kscr)
   nskip=-3
   call skiprz(ngendf,nskip)
   call findf(matd,1,451,ngendf)
   nth=nth1
   go to 125

   !--add temperature-dependent reaction to matrix.
  210 continue
   if (mth.eq.2.or.mth.eq.mti.or.mth.eq.mtc) go to 215
   go to 155
  215 continue
   if (mth.eq.2.and.igc.le.nth) go to 140
   if (mth.eq.mti.and.igc.gt.nth) go to 155
   if (mth.eq.mtc.and.igc.gt.nth.and.mtc.ne.0) go to 155
   itd=1
   do i=2,ng2
      ig2=ig2lo+i-2
      ig2c=icgrp(ig2)
      if (ig2c.ne.0.and.ig2c.le.ngnd) then
         jg2c=ngnd-ig2c+1
         loca=l+lz+(il-1)+nl*nz*(i-1)
         loc=jg+ngnd*(jg2c-1)
         s(loc)=s(loc)+scr(loca)
         if (mth.eq.2) then
            loc=jg+ngnd*(jg2c-1)
            elas(loc)=elas(loc)+scr(loca)
         endif
      endif
   enddo
   ip1=1
   lt=l1(jg)
   l1t=ngn-ig2lo-ng2+3
   l1tc=icgrp(l1t)
   if (l1tc.lt.lt.and.l1tc.ne.0) l1(jg)=l1tc
   if (mth.eq.2) l1e(jg)=l1tc
   lt=l2(jg)
   l2t=l1t+ng2-2
   l2tc=icgrp(l2t)
   if (l2tc.gt.lt.and.l2tc.ne.l1t.and.l2tc.gt.0.and.l2tc.le.ngnd)&
     l2(jg)=l2tc
   if (mth.eq.2.and.l2tc.ne.l1t) l2e(jg)=l2tc
  225 continue
   if (ig.lt.ngn) go to 140
   go to 125

   !--coarse group flux
  230 continue
   if (mth.ne.1) go to 155
   cflux(jg)=cflux(jg)+scr(l+lz)
   if (ig.lt.ngn) go to 140
   go to 125

   !--desired data is loaded
  250 continue
   if (ip1.eq.0) then
      write(strng,&
        '(''no required p1 matrices found for mat '',i4)') matd
      call error('sfile6',strng,' ')
   endif
   if (mti.gt.0.and.itd.eq.0) then
      write(strng,&
        '(''no temperature-dependent reactions for mat '',i4)') matd
      call error('sfile6',strng,'check mti')
   endif
   if (mtelas.le.0) then
      do i=1,ngnd
         if (l1e(i).lt.l1(i)) l1(i)=l1e(i)
         if (l2e(1).gt.l2(i)) l2(i)=l2e(i)
      enddo
      do i=1,ngndsq
         s(i)=s(i)+elas(i)
      enddo
   endif
   call sf6out(math,np1,ntp,ids,tempt,l1,l2,s)
   tempr(ntemp)=temp
   call repoz(-kscr)
   read(kscr) (l1(i),i=1,ngnd)
   read(kscr) (l2(i),i=1,ngnd)
   read(kscr) (s(i),i=1,ngndsq)
   mtelas=0
   go to 125

   !--this material is finished.
  300 continue
   call skiprz(ngendf,-1)
   np1=np1+1
   ids(np1)=iza
   do i=1,ntemp
      tempt(np1+newmat*(i-1))=tempr(i)
   enddo
   ntp(np1)=ntemp
  320 continue
   if (nmat.lt.newmat) go to 100

   !--sfile6 is finished.
   call sf6out(-1,np1,ntp,ids,tempt,l1,l2,s)
  400 continue
   deallocate(mtit)
   deallocate(mtct)
   call closz(mscr)
   call closz(kscr)
  410 continue
   write(nscr) eflag
   read(hflag,'(a4)') flag
   write(nscr) flag
   return
   end subroutine sfile6

   subroutine sf6out(matd,np1,ntp,ids,temp,l1,l2,s)
   !--------------------------------------------------------------------
   ! Write subfile 6.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use util ! provides error,repoz
   ! externals
   integer::matd,np1,ntp(*),ids(*),l1(*),l2(*)
   real(kr)::temp(*),s(*)
   ! internals
   integer::ncol1,ia,lone,ltwo,j,i,loc,ntemps,ip1,ident,it
   integer::ktemp,ig,nb,i1,lim
   character(60)::strng
   integer::k(70)
   real(kr)::b(70)
   integer::ncol=5
   integer::jtem=0

   !--store matrix on mscr.
   ncol1=ncol-1
   if (jtem.eq.0) call repoz(-mscr)
   jtem=jtem+1
   if (matd.gt.0) then
      do ia=1,ngnd
         lone=l1(ia)
         ltwo=l2(ia)
         if (ltwo.lt.lone) then
            write(strng,'(''for group '',i3,'' l1='',i2,&
              &'' is lt l2='',i2,'' at temp no. '',i2)')&
              ia,lone,ltwo,jtem
            call error('sf6out',strng,' ')
         endif
         j=0
         do i=lone,ltwo
            loc=ia+ngnd*(i-1)
            j=j+1
            b(j)=s(loc)
         enddo
         write(mscr) ia,lone,ltwo
         write(mscr) (b(i),i=1,j)
      enddo

   !--write subfile 6 on nscr.
   else
      if (np1.ne.0) then
         do i=1,np1
            k(i)=ids(i)
         enddo
         write(nscr) np1,(k(i),i=1,np1)
         if (iprint.ne.0) then
            write(nsyso,'(/&
                 &'' np1 '',i2/'' ---''//&
                 &'' nuclides''/'' --------''/&
                 &(11i7))')&
                 np1,(k(i),i=1,np1)
               write(nsyso,'(/&
                 &'' ntemp    temperatures''/&
                 &'' -----    ------------'')')
         endif
         do i=1,np1
            ntemps=ntp(i)
            do j=1,ntemps
               b(j)=temp(i+newmat*(j-1))
            enddo
            write(nscr) ntemps,(b(j),j=1,ntemps)
            if (iprint.eq.1) write(nsyso,&
              '(i5,5x,1p,5e12.4/(10x,5e12.4))')&
              ntemps,(b(j),j=1,ntemps)
         enddo
         call repoz(-mscr)
         do ip1=1,np1
            ident=ids(ip1)
            ntemps=ntp(ip1)
            do it=1,ntemps
               ktemp=nint(temp(ip1+newmat*(it-1)))
               if (iprint.eq.1) write(nsyso,'(/&
                 &'' p1 scattering matrix for nuclide '',i2,&
                 &'' ('',i6,'')'','' at '',i4,'' deg k''/&
                 &1x,5(''----------''),''--------''/&
                 &'' ia l1 l2    l1+0'',4(10x,''+'',i1)/&
                 &1x,7(''----------''),''----'')')&
                 ip1,ident,ktemp,(i,i=1,ncol1)
               do ig=1,ngnd
                  read(mscr) ia,lone,ltwo
                  write(nscr) ia,lone,ltwo
                  nb=ltwo-lone+1
                  read(mscr) (b(i),i=1,nb)
                  write(nscr) (b(i),i=1,nb)
                  if (iprint.ne.0) then
                     i1=1
                     lim=nb
                     if (lim.gt.ncol) lim=ncol
                     write(nsyso,'(3i3,5x,1p,5e12.4)')&
                       ia,lone,ltwo,(b(i),i=i1,lim)
                     do while (lim.ne.nb)
                        i1=lim+1
                        lim=nb
                        if ((lim-i1+1).gt.ncol) lim=i1+ncol-1
                        write(nsyso,'(14x,1p,5e12.4)')&
                          (b(i),i=i1,lim)
                     enddo
                  endif
               enddo
            enddo
         enddo
      endif
   endif
   return
   end subroutine sf6out

   subroutine cpmout(nlib,idat)
   !--------------------------------------------------------------------
   ! Write the input for CLIB on nout.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use util ! provides times,repoz,closz,error
   use physics ! provides amassn
   ! external
   integer::nlib,idat
   ! internals
   integer::jfiss,new,ngnd1,i,ip,j,np1,inf2,mp1o,nina
   integer::ident,indfis,jres,nf,ntemps,ip1,ip1opt,ifile
   integer::imat,ig,kres,jnf2,nmat,iwrite,nx,nw
   integer::need,low,lim,ihi,lr,ipf,itemp,ind,input
   integer::ia,l1,l2,nti,nfi,nresot
   real(kr)::time,pri,aw,flag,xid,decay
   character(60)::strng
   real(kr)::b(70)
   integer::k(70)
   real(kr)::tempr(7)
   !--------------------------------------------------------------------
   !    the following data cards must be updated with changes
   !    in the old library.  meaning of variables follows -
   !      nold    no. of nuclides
   !      idnto   idents of nuclides
   !      nreso   no. of resonance nuclides
   !      idreso  idents of resonance nuclides
   !      nf2     no. nuclides with greater than 1 resonance tabulation
   !      nfr2    idents of nf2 nuclides
   !      np1o    no. of nuclides having p1 matrices
   !      idp1o   idents of nuclides having p1 matrices
   integer::nold=65
   integer,dimension(80),parameter::idnto=(/&
     1001,1002,5000,5010,6000,8000,8001,13000,14000,&
     24000,25000,26000,28000,29063,47000,48000,49000,64154,&
     64155,64156,64157,64158,66164,71176,36083,45103,45105,47109,&
     54131,54135,55133,55134,55135,60143,60145,61147,61148,&
     61248,62147,62149,62150,62151,62152,63153,63154,&
     63155,401,402,92234,92235,92236,92238,93237,94238,94239,&
     94240,94241,94242,95241,95242,95243,96242,96244,1,302,&
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
   integer,parameter::nreso=4
   integer,dimension(20),parameter::idreso=(/&
     92235,92236,92238,94329,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
   integer::nf2=1
   integer,dimension(1),parameter::nfr2=(/92238/)
   integer,parameter::np1o=4
   integer,dimension(10),parameter::idp1o=(/&
     1001,8000,8001,1002,0,0,0,0,0,0/)
   !--------------------------------------------------------------------
   character(4),dimension(2)::fmt=(/'(6e1','2.0)'/)
   integer::iflgs(80)
      data iflgs/80*0/
   integer::izero=0
   integer::ione=1
   integer::ithree=3
   character(4)::hpri='pri '
   real(kr),parameter::zero=0.e0_kr
   real(kr),parameter::eps=1.e-5_kr
   real(kr)::eflag=1.e35_kr

   !--initialize.
   call timer(time)
   write(nsyso,'(/'' ***output***'',56x,f8.1,''s'')') time

   !--write input for subfile 1.
   read(hpri,'(a4)') pri
   write(nout,'(a3,3x,2i6)') pri,izero,izero
   write(nout,'(2a4)') fmt(1),fmt(2)
   write(nout,'(2a4)') fmt(1),fmt(2)
   write(nout,'(2a4)') fmt(1),fmt(2)
   write(nout,'(2a4)') fmt(1),fmt(2)
   if (mode.le.1) write(nout,'(12i6)') nold,izero
   if (mode.eq.2) write(nout,'(12i6)') izero,ione
   if (mode.le.1) write(nout,'(12i6)') ione
   jfiss=nfiss(1)
   if (nfiss(2).gt.jfiss) jfiss=nfiss(2)
   if (jfiss.eq.0) jfiss=22
   new=nold+newmat
   if (mode.eq.0)&
     write(nout,'(12i6)') nlib,idat,ngnd,nfg,nrg,jfiss,nold
   if (mode.eq.1)&
     write(nout,'(12i6)') nlib,idat,ngnd,nfg,nrg,jfiss,new
   if (mode.eq.2)&
     write(nout,'(12i6)') nlib,idat,ngnd,nfg,nrg,jfiss,newmat
   if (mode.eq.2) go to 100
   write(nout,'(12i6)') izero
   write(nout,'(12i6)') izero
   if (mode.eq.0) go to 102
   if (mode.eq.1) go to 106
  100 continue
   ngnd1=ngnd+1
   write(nout,'(1p,6e12.5)') (egb(i),i=1,ngnd1)
   write(nout,'(1p,6e12.5)') (cflux(i),i=1,ngnd)
   if (iu.eq.0.or.if4.eq.0) write(nout,'(1p,6e12.5)')&
        (zero,i=1,jfiss)
   if (iu.eq.1.and.if4.gt.0) write(nout,'(1p,6e12.5)')&
        (uff(i),i=1,jfiss)
   if (ipu.eq.0.or.if4.eq.0) write(nout,'(1p,6e12.5)')&
        (zero,i=1,jfiss)
   if (ipu.eq.1.and.if4.gt.0) write(nout,'(1p,6e12.5)')&
        (puff(i),i=1,jfiss)
   go to 108
  102 continue
   do i=1,newmat
      ip=ipost(i)
      iflgs(ip)=1
   enddo
   if (mode.eq.0) write(nout,'(80i1)') (iflgs(i),i=1,nold)
  106 continue
   if (mode.gt.0) write(nout,'(80i1)') (izero,j=1,nold)
  108 continue
   np1=0
   inf2=1
   mp1o=0
   do i=1,newmat
      if (mode.eq.0) ip=ipost(i)
      if (mode.gt.0) ip=0
      nina=ninat(i)
      aw=awrt(i)
      if (nina.lt.2) aw=aw*amassn
      if (nina.eq.2) nina=1
      if (nina.eq.3) nina=0
      ninat(i)=nina
      ident=idnt(i)
      if (mode.eq.0) ident=idnto(ip)
      do j=1,np1o
         if (ident.eq.idp1o(j)) mp1o=mp1o+1
      enddo
      indfis=indf(i)
      jres=irest(i)
      if (jres.gt.0.and.indfis.eq.0) indfis=1
      nf=0
      if (jres.gt.0) nf=1
      if (ident.eq.nfr2(inf2).and.mode.eq.0) nf=2
      if (inf2.lt.nf2.and.ident.eq.nfr2(inf2).and.mode.eq.0)&
        inf2=inf2+1
      ntemps=ntmp(i)
      ip1=0
      ip1opt=ip1tab(i)
      if (ip1opt.eq.0) ip1=1
      if (ip1.eq.1) np1=np1+1
      write(nout,'(f12.4,8i6)')&
        aw,ident,indfis,jres,nf,ntemps,nina,izero,ip1
   enddo

   !--write input for subfile 2.
   ifile=2
   if (mode.eq.0) write(nout,'(80i1)') (iflgs(i),i=1,nold)
   if (mode.eq.1) write(nout,'(80i1)') (izero,i=1,nold)
   call repoz(-nscr)
   if (nnina.ne.0) then
      do imat=1,nnina
         do ig=1,nrg
            read(nscr) (b(i),i=1,4)
            write(nout,'(1p,6e12.5)') (b(i),i=1,4)
         enddo
      enddo
   endif
   read(nscr) flag
   if (abs(flag-eflag).gt.eps) go to 910

   !--write input for subfile 3.
   ifile=3
   nresot=nreso-nf2
   if (mode.eq.1) go to 140
   if (mode.eq.2.and.nnina.eq.0) go to 170
   if (mode.eq.2) go to 145
   if (nres.gt.0) go to 145
  140 continue
   do i=1,nreso
      write(nout,'(i6/i6/i6)') izero,izero,izero
   enddo
   if (nres.eq.0) go to 170
  145 continue
   kres=1
   jnf2=0
   do 150 nmat=1,newmat
   jres=irest(nmat)
   if (jres.le.0) go to 150
   ident=idnt(nmat)
   iwrite=0
   if (mode.gt.0) go to 156
  152 continue
   if (ident.lt.idreso(kres)) go to 156
   if (ident.gt.idreso(kres)) go to 154
   iwrite=1
   inf2=0
   do while (inf2.lt.nf2)
      inf2=inf2+1
      if (ident.eq.nfr2(inf2)) write(nout,'(i6/i6/i6)')&
        izero,izero,izero
      if (ident.eq.nfr2(inf2)) jnf2=jnf2+1
   enddo
   if (kres.lt.nreso) kres=kres+1
   go to 156
  154 continue
   write(nout,'(i6/i6/i6)') izero,izero,izero
   inf2=0
   do while (inf2.lt.nf2)
      inf2=inf2+1
      if (idreso(kres).eq.nfr2(inf2)) write(nout,'(i6/i6/i6)')&
        izero,izero,izero
      if (idreso(kres).eq.nfr2(inf2)) jnf2=jnf2+1
   enddo
   if (kres.eq.nreso) go to 156
   kres=kres+1
   go to 152
  156 continue
   if (iwrite.gt.0) write(nout,'(i6)') ione
   read(nscr) xid
   inf2=0
   do while (inf2.lt.nf2)
      inf2=inf2+1
   enddo
   write(nout,'(1p,6e12.5)') xid
   if (iwrite.gt.0) write(nout,'(i6)') ione
   read(nscr) ntemps,nx
   write(nout,'(2i6)') ntemps,nx
   read(nscr) (b(i),i=1,ntemps)
   write(nout,'(1p,6e12.5)') (b(i),i=1,ntemps)
   read(nscr) (b(i),i=1,nx)
   write(nout,'(1p,6e12.5)') (b(i),i=1,nx)
   if (iwrite.gt.0) write(nout,'(i6)') ione
   nw=ntemps*nx
   indfis=indf(nmat)
  158 continue
   do j=1,nrg
      read(nscr) (b(i),i=1,nw)
      write(nout,'(1p,6e12.5)') (b(i),i=1,nw)
   enddo
   if (indfis.lt.2) go to 150
   indfis=1
   go to 158
  150 continue
   if (mode.gt.0) go to 170
   if (kres.ge.nresot) go to 170
   need=nreso-kres+nf2-jnf2
   do i=1,need
      write(nout,'(i6/i6/i6)') izero,izero,izero
   enddo
  170  continue
   read(nscr) flag
   if (abs(flag-eflag).gt.eps) go to 910

   !--write input for subfile 4.
   ifile=4
   if (if4.ne.0) then
      if (mode.ne.2) then
         write(nout,'(80i1)') (iflgs(i),i=1,nold)
         low=1
         lim=nold/8
         if (lim*8.lt.nold) lim=lim+1
         do i=1,lim
            ihi=low+7
            lr=low/8
            do j=1,80
               iflgs(j)=0
            enddo
            if (mode.ne.1) then
               do j=1,newmat
                  ip=ipost(j)
                  if (ip.ge.low.and.ip.le.ihi) then
                     ipf=(ip-lr*8-1)*10
                     iflgs(ipf+1)=1
                     iflgs(ipf+2)=1
                     iflgs(ipf+3)=1
                     iflgs(ipf+4)=1
                     iflgs(ipf+5)=1
                  endif
               enddo
            endif
            write(nout,'(80i1)') (iflgs(j),j=1,80)
            low=low+8
         enddo
      endif
      nmat=0
      do while (nmat.lt.newmat)
         nmat=nmat+1
         ntemps=ntmp(nmat)
         read(nscr) (tempr(i),i=1,ntemps)
         write(nout,'(8f10.0)') (tempr(i),i=1,ntemps)
         itemp=0
         nina=ninat(nmat)
         lim=4
         if (nina.gt.0) lim=1
         do while (itemp.lt.ntemps)
            write(nout,'(2i6)') ione,izero
            do i=1,lim
               read(nscr) (b(j),j=1,ngnd)
               write(nout,'(1p,6e12.5)') (b(j),j=1,ngnd)
            enddo
            if (nina.le.0) then
               ind=3
               input=0
               write(nout,'(2i6)') ind,input
               do j=1,ngnd
                  read(nscr) ia,l1,l2
                  write(nout,'(3i6)') ia,l1,l2
                  nw=l2-l1+1
                  read(nscr) (b(i),i=1,nw)
                  write(nout,'(1p,6e12.5)') (b(i),i=1,nw)
               enddo
            endif
            itemp=itemp+1
         enddo
      enddo
   endif
   read(nscr) flag
   if (abs(flag-eflag).gt.eps) go to 910

   !--write input for subfile 5.
   ifile=5
   if (if5.eq.0) write(nout,'(i6)') izero
   if (if5.gt.0) write(nout,'(i6)') ione
   if (if5.eq.0) go to 218
   read(nscr) nti,(k(i),i=1,nti)
   write(nout,'(i6/(12i6))') nti,(k(i),i=1,nti)
   read(nscr) nfi,(k(i),i=1,nfi)
   write(nout,'(i6/(12i6))') nfi,(k(i),i=1,nfi)
   do i=1,nti
      read(nscr) decay,(b(j),j=1,nfi)
      write(nout,'(1p,6e12.5)') decay
      write(nout,'(1p,6e12.5)') (b(j),j=1,nfi)
   enddo
   read(nscr) eflag
   if (abs(flag-eflag).gt.eps) go to 910

   !--write input for subfile 6.
  218 continue
   ifile=6
   if (if4.ne.0) then
      ind=0
      if (mp1o.gt.0) ind=1
      write(nout,'(3i6)') ind,ithree,input
      if (np1.ne.0) then
         if (mode.ne.2.and.ind.ne.1) then
            write(nout,'(i6)') ione
         endif
         read(nscr) np1,(k(i),i=1,np1)
         write(nout,'(7i6)') np1,(k(i),i=1,np1)
         ip1=0
         do while (ip1.lt.np1)
            ip1=ip1+1
            read(nscr) ntemps,(b(i),i=1,ntemps)
            ntmp(ip1)=ntemps
            write(nout,'(i6,1p,6e12.5)') ntemps,(b(i),i=1,ntemps)
         enddo
         ip1=0
         do while (ip1.lt.np1)
            ip1=ip1+1
            itemp=0
            ntemps=ntmp(ip1)
            do while (itemp.lt.ntemps)
               do j=1,ngnd
                  read(nscr) ia,l1,l2
                  write(nout,'(3i6)') ia,l1,l2
                  nw=l2-l1+1
                  read(nscr) (b(i),i=1,nw)
                  write(nout,'(1p,6e12.5)') (b(i),i=1,nw)
               enddo
               itemp=itemp+1
            enddo
         enddo
      else
         write(nout,'(i6)') izero
         write(nout,'(i6)') izero
      endif
   endif
   read(nscr) flag
   if (abs(flag-eflag).gt.eps) go to 910

   !--write end of input flag.
   read(nscr) flag
   write(nout,'(2a4)') flag

   !--cpmout is finished.
   call closz(nscr)
   return

   !--end flag not found on input tape
  910 continue
   write(strng,&
     '(''expected end flag not found on nscr for file '',i2)')&
     ifile
   call error('cpmout',strng,' ')
   end subroutine cpmout

end module powm


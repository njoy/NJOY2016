module thermm
   ! provides thermr for NJOY2016
   use locale
   implicit none
   private
   public thermr

   ! global variables
   integer::nendf,nin,nout,nscr,nscr2
   integer::matde,nbin,iprint,ncds,matdp,natom,ntemp,&
     iinc,iform,ncdse
   real(kr)::za,awr,tol,emax
   real(kr)::sb,az,tevz,teff,sb2,az2,teff2
   real(kr)::cliq
   integer::lasym
   integer::lat
   real(kr),dimension(:),allocatable::scr
   real(kr),dimension(:),allocatable::bufo,bufn

   ! array for user temperatures
   real(kr),dimension(:),allocatable::tempr

   ! array for thermal elastic data
   real(kr),dimension(:),allocatable::fl

   ! arrays for thermal inelastic cross section
   real(kr),dimension(:),allocatable::esi,xsi

   integer,parameter::nbuf=1000
   integer,parameter::nwscr=500000

contains

   subroutine thermr
   !-------------------------------------------------------------------
   !
   !  Generate neutron scattering cross sections and point-to-point
   !  scattering kernels in the thermal range.  The coding can
   !  generate incoherent inelastic cross sections and distributions
   !  for a free gas, incoherent inelastic cross sections and
   !  distributions from read-in S(alpha,beta,T) data, coherent
   !  elastic scattering cross sections for crystalline materials,
   !  and incoherent elastic cross sections and angular distributions
   !  for hydrogenous solids.
   !
   !  The pointwise scattering cross sections and distributions are
   !  added to an existing PENDF tape.  Cross sections are added in
   !  MF3 and distributions are written in MF6, both using mtref for
   !  inelastic and mtref+1 for elastic (if any).
   !
   !  Multiple scattering types (i.e., H free and H in H2O) can be
   !  written on one PENDF tape by using different values of mtref
   !  for each thermr run.  If data for one mtref is already on the
   !  tape, it will be replaced with the new cross sections.
   !
   !  The energy grid for coherent scattering is determined
   !  adaptively so as to represent the sharp Bragg edges to a
   !  specified tolerance using linear interpolation.  No angular
   !  information is provided for coherent scattering.  It can be
   !  deduced by subsequent processing from the Bragg edges in the
   !  cross section.
   !
   !  The incident energy grid for incoherent inelastic scattering
   !  is fixed in the code, although it is scaled to higher values
   !  for very large incident energies.  There are two options for
   !  representing the energy and angle distributions for secondary
   !  inelastic neutrons.
   !
   !  The standard approach is to generate an energy grid for the
   !  secondary energy distribution integrated over angle
   !  adaptively.  The angular distribution for each secondary
   !  energy is given using a set of equally probable emission
   !  cosines.  This is E-E'-mu ordering, and it is represented
   !  using a special format in MF6.
   !
   !  An alternate approach introduced in NJOY2010 produces E-mu-E'
   !  ordering using MF6/Law7.  A grid of mu values is generated
   !  adaptively for linear interpolation.  In addition, secondary
   !  energy grids are determined adaptively for each mu value.
   !
   !  The sections of File 6 produced by thermr have several
   !  special features.  The LIP flag is used to identify the
   !  various kinds of thermal data:  LIP=-1 means incoherent
   !  inelastic data, LIP=-2 means incoherent inelastic data, and
   !  LIP=-nbrag means coherent elastic data with nbrag Bragg
   !  edges.  For incoherent elastic or inelastic data using
   !  LAW=1, the angular data are equally probable cosines.
   !
   !  For ENDF 3 to 5 formats, the constants used for coherent
   !  elastic, incoherent elastic, and short-collision-time calcu-
   !  lations are obtained from internal data statements based on
   !  the original general atomic report on the evaluations
   !  (GA-8774 revised, ENDF-269, July 1978).
   !
   !  For ENDF-6 format libraries, these constants are included
   !  in the format.
   !
   !---input specifications (free format)--------------------------
   !
   !  card 1
   !     nendf      endf tape for mf7 data
   !     nin        old pendf tape
   !     nout       new pendf tape
   !  card 2
   !     matde      material desired on endf tape
   !     matdp      material desired on pendf tape
   !     nbin       number of equi-probable angles
   !     ntemp      number of temperatures (default = 1)
   !     iinc       inelastic options
   !                   0     none
   !                   1     compute as free gas
   !                   2     read s(a,b) and compute matrix
   !     icoh       elastic options
   !                   0     none
   !                   1     compute using ENDF6 format data
   !                   --------or for earlier formats
   !                   1     graphite
   !                   2     beryllium
   !                   3     beryllium oxide
   !                  11     polyethylene
   !                  12     h(zrh)
   !                  13     zr(zrh)
   !     iform      output format for inelastic distributions
   !                  0      E-E'-mu ordering (MF6 special)
   !                  1      E-mu-E' ordering (MF6/Law7)
   !     natom      number of principal atoms
   !     mtref      mt for inelastic reaction (221-250 only)
   !     iprint     print option (0=minimum, 1=maximum,
   !                2=max. normal + intermediate results)
   !                (default=0)
   !  card 3
   !     tempr      temperatures (kelvin)
   !  card 4
   !     tol        tolerance
   !     emax       maximum energy for thermal treatment
   !                (for temperatures greater than 3000,
   !                emax and the energy grid are scaled by
   !                temp/3000.  free gas only.)
   !
   !       nendf can be ENDF6 format (e.g., from leapr) while
   !       nin and nout are ENDF4 or 5 format, if desired.
   !
   !-------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso,nsyse
   use endf   ! provides endf routines and variables
   use util   ! provides timer,error,openz,repoz,loada,skiprz,closz
   ! locals
   integer::icoh,i,icopy,idis,nb,nw,iold,inew,indexc,indexi,mtref
   integer::iverp,itemp,nwb,it,ntape,nex,ne,np,isave
   integer::lthr
   real(kr)::e,enext,emaxs,time,sz2,t,templ,s,temp
   real(kr)::ex(2)
   real(kr),dimension(:),allocatable::eftemp,eftmp2
   real(kr),parameter::s1099=3.76e0_kr
   real(kr),parameter::a1099=15.86e0_kr
   real(kr),parameter::s1095=4.74e0_kr
   real(kr),parameter::a1095=11.90e0_kr
   real(kr),parameter::eps=1.e-4_kr
   real(kr),parameter::break=3000.e0_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::up=1.00001e0_kr
   real(kr),parameter::zero=0.e0_kr

   !--initialize
   call timer(time)
   write(nsyso,'(/&
     &'' thermr...compute thermal scattering cross sections '',&
     &''and distributions'',f8.1,''s'')') time
   write(nsyse,'(/'' thermr...'',59x,f8.1,''s'')') time
   sz2=0
   allocate(scr(nwscr))
   allocate(bufn(nbuf))
   allocate(bufo(nbuf))
   scr=0
   indexc=-1
   indexi=-1

   !--read user input.
   read(nsysi,*) nendf,nin,nout
   if (nin.eq.0) call error('thermr','nin=0.',' ')
   if (nin.lt.0.and.nout.gt.0)&
     call error('thermr','mode conversion not allowed.',' ')
   if (nin.gt.0.and.nout.lt.0)&
     call error('thermr','mode conversion not allowed.',' ')
   call openz(nendf,0)
   call openz(nin,0)
   call openz(nout,1)
   iprint=0
   ntemp=1
   read(nsysi,*) matde,matdp,nbin,ntemp,iinc,icoh,iform,natom,mtref,iprint
   if (mtref.lt.221.or.mtref.gt.250)&
     call error('thermr','illegal reference mt.',' ')
   allocate(tempr(ntemp))
   allocate(eftemp(ntemp))
   allocate(eftmp2(ntemp))
   read(nsysi,*) (tempr(i),i=1,ntemp)
   do i=1,ntemp
      eftemp(i)=0
      eftmp2(i)=0
      if (matde.eq.0) eftemp(i)=tempr(i)
   enddo
   read(nsysi,*) tol,emax

   !--check for endf-6 format data
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
   call skiprz(nendf,-1)

   !--set up internal data for older endf versions
   if (iverf.lt.6) then

      !--check for mixed s(a,b) cases (beo, benzine)
      !--supply appropriate parameters
      sz2=0
      az2=0
      if (matde.eq.1099) sz2=s1099
      if (matde.eq.1099) az2=a1099
      if (matde.eq.1095) sz2=s1095
      if (matde.eq.1095) az2=a1095
      sb2=sz2
      if (az2.ne.zero) sb2=sz2*((az2+1)/az2)**2

         !--default effective temperatures to standard values if
         !--available, otherwise set them to the material temperature
         if (matde.ne.0) then
            call gateff(tempr,eftemp,ntemp,matde)
            if (sz2.ne.zero) then
               call gatef2(tempr,eftmp2,ntemp,matde)
            endif
         endif
      endif

   !--write parameters.
   write(nsyso,'(/&
     &'' unit for endf tape ................... '',i10/&
     &'' unit for input pendf tape ............ '',i10/&
     &'' unit for output pendf tape ........... '',i10)')&
     nendf,nin,nout
   write(nsyso,'(/&
     &'' material to be processed (endf) ...... '',i10/&
     &'' material to be processed (pendf) ..... '',i10/&
     &'' number of angle bins ................. '',i10/&
     &'' number of temperatures ............... '',i10/&
     &'' inelastic option ..................... '',i10/&
     &'' elastic option ....................... '',i10/&
     &'' MF6 format option .................... '',i10/&
     &'' number of principal atoms ............ '',i10/&
     &'' reference mt ......................... '',i10/&
     &'' print option (0 min, 1 max) .......... '',i10)')&
     matde,matdp,nbin,ntemp,iinc,icoh,iform,natom,mtref,iprint
   write(nsyso,'(&
     &'' temperatures (kelvin) ................ '',1p,e10.4)')&
     tempr(1)
   if (ntemp.gt.1) write(nsyso,'(40x,1p,e10.4)')&
     (tempr(i),i=2,ntemp)
   write(nsyso,'(&
     &'' tolerance ............................ '',1p,e10.4/&
     &'' max energy for thermal treatment ..... '',e10.4)')&
     tol,emax
   if (iinc.eq.2.and.iverf.lt.6) then
      write(nsyso,'(/&
        &'' parameters for sct app.''/&
        &'' -----------------------'')')
      write(nsyso,'(&
        &'' effective temperatures ............... '',1p,e10.4)')&
        eftemp(1)
      if (ntemp.gt.1) write(nsyso,'(40x,1p,e10.4)')&
        (eftemp(i),i=2,ntemp)
      if (sz2.ne.zero) then
         write(nsyso,'(&
           &'' free xsec for second atom ............ '',1p,e10.4/&
           &'' mass ratio for second atom ........... '',1p,e10.4/&
           &'' eff. temps for second atom ........... '',1p,e10.4)')&
           sz2,az2,eftmp2(1)
         if (ntemp.gt.1) write(nsyso,'(40x,1p,e10.4)')&
           (eftmp2(i),i=2,ntemp)
      endif
   endif
   write(nsyso,'(/'' endf uses endf-'',i1,'' format'')') iverf

   !--initialize i/o units
   iold=10
   inew=11
   nscr=12
   if (nout.lt.0) nscr=-nscr
   nscr2=13
   if (nout.lt.0) nscr2=-nscr2
   call openz(-iold,1)
   call openz(-inew,1)
   call openz(nscr,1)
   call openz(nscr2,1)
   call repoz(nendf)
   call repoz(nout)
   call repoz(nin)
   call repoz(nscr2)

   !--copy through to desired material
   nsh=0
   call tpidio(nin,nout,0,scr,nb,nw)
  110 continue
   call contio(nin,0,0,scr,nb,nw)
   if (math.eq.matdp) go to 120
   if (math.lt.0)&
     call error('thermr','desired material not on pendf tape.',' ')
   call contio(0,nout,0,scr,nb,nw)
   call tomend(nin,nout,0,scr)
   go to 110
  120 continue
   call contio(nin,0,0,scr(7),nb,nw)
   if (n1h.ne.0) then
      iverp=4
   else if (n2h.eq.0) then
      iverp=5
   else
      iverp=6
   endif
   call skiprz(nin,-1)
   write(nsyso,'(/'' pendf uses endf-'',i1,'' format'')') iverp

   !--loop over desired temperatures
   itemp=1
   icopy=0
  130 continue
   call contio(0,0,nscr2,scr,nb,nw)
   if (math.ne.matdp)&
     call error('thermr','desired material not on pendf tape.',' ')
   nwb=0
   za=scr(1)
   awr=scr(2)
   az=scr(2)
   if (iverp.ge.5) call contio(nin,0,nscr2,scr,nb,nw)
   if (iverp.ge.6) call contio(nin,0,nscr2,scr,nb,nw)
   call hdatio(nin,0,nscr2,scr,nb,nw)
   t=scr(1)
   if (abs(t-tempr(itemp)).le.eps*tempr(itemp)) go to 180
   if (t.gt.tempr(itemp))&
     call error('thermr','desired temperature not on tape.',' ')

   !--check for skipped temperatures to be copied
   if (itemp.eq.1) icopy=icopy+1
   if (itemp.eq.1) go to 160
   templ=tempr(itemp-1)
   it=0
  150 continue
   it=it+1
   if (it.ge.itemp) go to 160
   if (abs(t-tempr(it)).le.eps*tempr(it)) go to 150
   if (abs(t-templ).le.eps*templ) go to 150
   if (t.lt.templ) go to 150
   icopy=icopy+1
   go to 150
  160 continue
   if (icopy.eq.0) call repoz(nscr2)
   ntape=nscr2
   if (icopy.eq.0) ntape=0
   call tomend(nin,0,ntape,scr)
   call contio(nin,0,0,scr,nb,nw)
   go to 130
  180 continue
   call findf(matdp,3,2,nin)

   !--save elastic cross section on scratch file.
   call contio(nin,0,0,scr,nb,nw)
   nex=2
   ne=0
   e=0
   call gety1(e,enext,idis,s,nin,scr)
   emaxs=emax
   if (t.gt.break) emax=emax*t/break
  200 continue
   e=enext
   call gety1(e,enext,idis,s,nin,scr)
   ex(1)=e
   ex(2)=s
   if (e.gt.emax*(1+small)) ex(2)=0
   ne=ne+1
   if (e.gt.emax*(1+small)) ne=-ne
   call loada(ne,ex,nex,inew,bufn,nbuf)
   if (ne.lt.0) go to 210
   if (enext.gt.emax*(1+small)) enext=emax
   if (abs(e-emax).lt.small*emax) enext=up*emax
   go to 200
  210 continue
   np=-ne
   emax=emaxs
   isave=iold
   iold=inew
   inew=isave
   call skiprz(nin,-1)

   !--set up for elastic calculation
   lthr=0
   if (iverf.ge.6.and.nendf.ne.0) then
      call findf(matde,7,0,nendf)
      call contio(nendf,0,0,scr,nb,nw)
      if (mth.eq.2) then
         lthr=l1h
         icoh=10*lthr
         temp=tempr(itemp)
         call rdelas(temp,lthr,nendf,nwb,indexc,indexi)
      endif
   endif

   !--compute incoherent cross sections
   if (iinc.ne.0) then
      temp=tempr(itemp)
      teff=eftemp(itemp)
      teff2=0
      if (sz2.gt.zero) teff2=eftmp2(itemp)
      if (iinc.eq.2) temp=t
      call calcem(temp,itemp,iold,inew,np,nex,mtref)
   endif

   !--compute thermal elastic cross sections
   if (icoh.gt.0.and.icoh.le.10) then        ! only coherent
      call coh(icoh,itemp,iold,inew,np,nex,mtref+1,indexc)
   else if (icoh.gt.10.and.icoh.le.20) then  ! only incoherent
      call iel(icoh,itemp,iold,inew,np,nex,mtref+1,indexi)
   else if (icoh.gt.20) then                 ! both coherent and incoherent
     call coh(10,itemp,iold,inew,np,nex,mtref+1,indexc)
     call iel(20,itemp,iold,inew,np,nex,mtref+2,indexi)
   endif

   !--write new pendf tape.
   call tpend(iold,itemp,np,nex,icoh,icopy,mtref)

   !--continue temperature loop
   if (allocated(fl)) deallocate(fl)
   if (itemp.eq.ntemp) go to 310
   itemp=itemp+1
   call repoz(nscr2)
   call contio(nin,0,0,scr,nb,nw)
   go to 130

   !--finished with all temperatures.
   !--copy rest of pendf tape and terminate job.
  310 continue
   call totend(nin,nout,0,scr)
   call repoz(nendf)
   call repoz(nin)
   call repoz(nout)
   call closz(nendf)
   call closz(nin)
   call closz(nout)
   call closz(-iold)
   call closz(-inew)
   call closz(nscr)
   call closz(nscr2)
   deallocate(scr)
   deallocate(bufn)
   deallocate(bufo)
   deallocate(tempr)
   deallocate(eftemp)
   deallocate(eftmp2)
   if (allocated(esi)) deallocate(esi)
   if (allocated(xsi)) deallocate(xsi)
   if (allocated(fl)) deallocate(fl)
   call timer(time)
   write(nsyso,'(69x,f8.1,''s''/&
     &1x,7(''**********''),''*******'')') time
   return
   end subroutine thermr

   subroutine rdelas(temp,lthr,nin,nwds,indexc,indexi)
   !-------------------------------------------------------------------
   ! Read in elastic data (coherent or incoherent or both) for the
   ! desired temperature from an ENDF6 format file.
   !-------------------------------------------------------------------
   use util ! provides error
   use endf ! provides endf routines and variables
   ! externals
   integer::lthr,nin,nwds,indexc,indexi
   real(kr)::temp
   ! internals
   integer::l,np,nr,lt,nb,nw,k,i,j,ifound
   real(kr)::tnow
   real(kr),dimension(:),allocatable::a
   integer,parameter::na=1000000

   !--temp storage to read in data
   allocate(a(na))

   !--initialise
   ifound=0
   indexc=-1
   indexi=-1

   !--read the first table in MF7 MT2
   !  this is S(E) for the first temperature for coherent or mixed mode
   !  this is W'(T) for incoherent
   l=1
   call tab1io(nin,0,0,a(l),nb,nw)
   l=l+nw
   do while (nb.ne.0)
      call moreio(nin,0,0,a(l),nb,nw)
      l=l+nw
      if (l.gt.na) call error('rdelas',&
        'too much elastic data','increase na')
   enddo

   !--get table parameters
   np=n2h
   nr=n1h
   lt=l1h

   !--current size of the array
   nwds=l-1

   !--read temperatures for coherent elastic or mixed mode
   if (lthr.ne.2) then

      !--set coherent data index
      indexc=1

      !--current temperature
      tnow=a(1)

      !--loop over the available temperatures and find the one we need
      do j=1,lt

         !--read list of S values for this temperature
         l=nwds+1
         call listio(nin,0,0,a(l),nb,nw)
         l=l+nw
         do while (nb.ne.0)
            call moreio(nin,0,0,a(l),nb,nw)
            l=l+nw
         enddo

         !--verify if this is the one we need and store S(E) if it is
         tnow=a(nwds+1)
         if (abs(tnow-temp).lt.temp/1000+5) then
            ifound=1
            l=nwds+6
            k=6+2*nr
            do i=1,np
               l=l+1
               k=k+2
               a(k)=a(l)
            enddo
         endif
      enddo
   endif

   !--set index for incoherent only
   if (lthr.eq.2) then

      !--set incoherent data index
      indexi=1

   !--read the W'(T) data for mixed mode
   else if (lthr.eq.3) then

      !--set incoherent data index
      indexi=nwds+1

      l=indexi
      call tab1io(nin,0,0,a(l),nb,nw)
      l=l+nw
      do while (nb.ne.0)
         call moreio(nin,0,0,a(l),nb,nw)
         l=l+nw
         if (l.gt.na) call error('rdelas',&
           'too much elastic data','increase na')
      enddo
      nwds=l-1
   endif

   !--move data to global fl array
   !--lthr=1: fl is a single table of S(E) for coherent elastic
   !--lthr=2: fl is a single table of W'(T) for incoherent elastic
   !--lthr=3: fl is a table of S(E) and a table of W'(T) for mixed mode
   allocate(fl(nwds))
   do l=1,nwds
      fl(l)=a(l)
   enddo
   deallocate(a)
   return
   end subroutine rdelas

   subroutine gateff(temp,eftemp,ntemp,mat)
   !-------------------------------------------------------------------
   ! Default effective temperatures to values from General Atomic
   ! report, if available.
   !       1002      h(h2o)
   !       1004      d(d2o)
   !       1064      be
   !       1065      graphite
   !       1095      benzine (c6h6)
   !       1096      zr(zrh)
   !       1097      h(zrh)
   !       1099      beo
   !       1114      h(ch2)
   !-------------------------------------------------------------------
   ! externals
   integer::ntemp,mat
   real(kr)::temp(ntemp),eftemp(ntemp)
   ! internals
   integer::i,j,jmat
   real(kr)::diff,test
   integer::ntabl=67
   real(kr),dimension(3,67),parameter::tabl=reshape((/&
     1002.e0_kr,296.0e0_kr,1396.8e0_kr,&
     1002.e0_kr,350.0e0_kr,1411.6e0_kr,&
     1002.e0_kr,400.0e0_kr,1427.4e0_kr,&
     1002.e0_kr,450.0e0_kr,1444.9e0_kr,&
     1002.e0_kr,500.0e0_kr,1464.1e0_kr,&
     1002.e0_kr,600.0e0_kr,1506.8e0_kr,&
     1002.e0_kr,800.0e0_kr,1605.8e0_kr,&
     1002.e0_kr,1000.e0_kr,1719.8e0_kr,&
     1004.e0_kr,296.0e0_kr,940.91e0_kr,&
     1004.e0_kr,350.0e0_kr,961.62e0_kr,&
     1004.e0_kr,400.0e0_kr,982.93e0_kr,&
     1004.e0_kr,450.0e0_kr,1006.1e0_kr,&
     1004.e0_kr,500.0e0_kr,1030.9e0_kr,&
     1004.e0_kr,600.0e0_kr,1085.1e0_kr,&
     1004.e0_kr,800.0e0_kr,1209.0e0_kr,&
     1004.e0_kr,1000.e0_kr,1350.0e0_kr,&
     1064.e0_kr,296.0e0_kr,405.64e0_kr,&
     1064.e0_kr,400.e0_kr,484.22e0_kr,&
     1064.e0_kr,500.0e0_kr,568.53e0_kr,&
     1064.e0_kr,600.e0_kr,657.66e0_kr,&
     1064.e0_kr,700.e0_kr,749.69e0_kr,&
     1064.e0_kr,800.e0_kr,843.63e0_kr,&
     1064.e0_kr,1000.e0_kr,1035.e0_kr,&
     1064.e0_kr,1220.e0_kr,1229.3e0_kr,&
     1065.e0_kr,296.0e0_kr,713.39e0_kr,&
     1065.e0_kr,400.0e0_kr,754.68e0_kr,&
     1065.e0_kr,500.0e0_kr,806.67e0_kr,&
     1065.e0_kr,600.0e0_kr,868.38e0_kr,&
     1065.e0_kr,700.0e0_kr,937.64e0_kr,&
     1065.e0_kr,800.0e0_kr,1012.7e0_kr,&
     1065.e0_kr,1000.e0_kr,1174.9e0_kr,&
     1065.e0_kr,1200.e0_kr,1348.2e0_kr,&
     1065.e0_kr,1600.e0_kr,1712.9e0_kr,&
     1065.e0_kr,2000.e0_kr,2091.0e0_kr,&
     1095.e0_kr,296.0e0_kr,1165.9e0_kr,&
     1095.e0_kr,350.0e0_kr,1177.8e0_kr,&
     1095.e0_kr,400.0e0_kr,1191.4e0_kr,&
     1095.e0_kr,450.0e0_kr,1207.7e0_kr,&
     1095.e0_kr,500.0e0_kr,1226.0e0_kr,&
     1095.e0_kr,600.0e0_kr,1268.7e0_kr,&
     1095.e0_kr,800.0e0_kr,1373.4e0_kr,&
     1095.e0_kr,1000.e0_kr,1497.7e0_kr,&
     1096.e0_kr,296.0e0_kr,317.27e0_kr,&
     1096.e0_kr,400.0e0_kr,416.29e0_kr,&
     1096.e0_kr,500.0e0_kr,513.22e0_kr,&
     1096.e0_kr,600.0e0_kr,611.12e0_kr,&
     1096.e0_kr,700.0e0_kr,709.60e0_kr,&
     1096.e0_kr,800.0e0_kr,808.43e0_kr,&
     1096.e0_kr,1000.e0_kr,1006.8e0_kr,&
     1096.e0_kr,1200.e0_kr,1205.7e0_kr,&
     1097.e0_kr,296.0e0_kr,806.79e0_kr,&
     1097.e0_kr,400.0e0_kr,829.98e0_kr,&
     1097.e0_kr,500.0e0_kr,868.44e0_kr,&
     1097.e0_kr,600.0e0_kr,920.08e0_kr,&
     1097.e0_kr,700.0e0_kr,981.82e0_kr,&
     1097.e0_kr,800.0e0_kr,1051.1e0_kr,&
     1097.e0_kr,1000.e0_kr,1205.4e0_kr,&
     1097.e0_kr,1200.e0_kr,1373.4e0_kr,&
     1099.e0_kr,296.0e0_kr,596.4e0_kr,&
     1099.e0_kr,400.0e0_kr,643.9e0_kr,&
     1099.e0_kr,500.0e0_kr,704.6e0_kr,&
     1099.e0_kr,600.0e0_kr,775.3e0_kr,&
     1099.e0_kr,800.0e0_kr,935.4e0_kr,&
     1099.e0_kr,1000.e0_kr,1109.8e0_kr,&
     1099.e0_kr,1200.e0_kr,1292.3e0_kr,&
     1114.e0_kr,296.0e0_kr,1222.0e0_kr,&
     1114.e0_kr,350.0e0_kr,1239.0e0_kr/),(/3,67/))
   real(kr),parameter::zero=0

   do i=1,ntemp
      if (eftemp(i).eq.zero) then
         do j=1,ntabl
            jmat=nint(tabl(1,j))
            if (jmat.eq.mat) then
               test=5
               diff=tabl(2,j)-temp(i)
               if (abs(diff).le.test) then
                  eftemp(i)=tabl(3,j)
               endif
            endif
         enddo
         if (eftemp(i).eq.zero) eftemp(i)=temp(i)
      endif
   enddo
   return
   end subroutine gateff

   subroutine gatef2(temp,eftmp2,ntemp,mat)
   !-------------------------------------------------------------------
   ! Default effective temperatures for second atom in mixed
   ! moderators from GA report:
   !       1095      benzine (c6h6)
   !       1099      beo
   !-------------------------------------------------------------------
   ! externals
   integer::ntemp,mat
   real(kr)::temp(ntemp),eftmp2(ntemp)
   ! internals
   integer::i,j,jmat
   real(kr)::diff,test
   integer::ntabl=16
   real(kr),dimension(3,16),parameter::tabl=reshape((/&
     1095.e0_kr,296.0e0_kr,685.54e0_kr,1095.e0_kr,350.e0_kr,&
     712.02e0_kr,1095.e0_kr,400.e0_kr,738.97e0_kr,1095.e0_kr,&
     450.e0_kr,768.10e0_kr,1095.e0_kr,500.e0_kr,799.22e0_kr,&
     1095.e0_kr,600.e0_kr,866.63e0_kr,1095.e0_kr,800.e0_kr,&
     1017.3e0_kr, 1095.e0_kr,1000.e0_kr,1182.3e0_kr,&
     1099.e0_kr,296.0e0_kr,427.8e0_kr,1099.e0_kr,400.e0_kr,&
     502.8e0_kr,1099.e0_kr,500.e0_kr,584.3e0_kr,1099.e0_kr,&
     600.e0_kr,671.3e0_kr,1099.e0_kr,700.e0_kr,761.6e0_kr,&
     1099.e0_kr,800.e0_kr,854.2e0_kr,1099.e0_kr,1000.e0_kr,&
     1043.7e0_kr,1099.e0_kr,1200.e0_kr,1236.6e0_kr/),(/3,16/))
   real(kr),parameter::zero=0

   do i=1,ntemp
      if (eftmp2(i).eq.zero) then
         do j=1,ntabl
            jmat=nint(tabl(1,j))
            if (jmat.eq.mat) then
               test=5
               diff=tabl(2,j)-temp(i)
               if (abs(diff).le.test) then
                  eftmp2(i)=tabl(3,j)
               endif
            endif
         enddo
         if (eftmp2(i).eq.zero) eftmp2(i)=temp(i)
      endif
   enddo
   return
   end subroutine gatef2

   subroutine coh(lat,itemp,iold,inew,ne,nex,mtref,index)
   !-------------------------------------------------------------------
   ! Compute the coherent scattering cross sections for a crystalline
   ! material.  The cross section is computed on an energy grid
   ! chosen adaptively to represent the actual function within a
   ! given tolerance with linear-linear interpolation.
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use util ! provieds error,finda,loada
   ! externals
   integer::lat,itemp,iold,inew,ne,nex,mtref,index
   ! internals
   integer::nl,imax,nx,nj,i,nlt,nlt1,nb,nw,nbragg
   integer::ix,j,iex,il,isave,ltt
   real(kr)::temp,e,enext,en,xm,ym,test
   integer,parameter::nlmax=6
   real(kr)::s(nlmax),ej(20),ex(20),x(5),y(5),z(5),b(12)
   real(kr),dimension(:,:),allocatable::stk
   real(kr),parameter::half=0.5e0_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::tolmin=1.e-6_kr
   real(kr),parameter::eps=3.e-5_kr

   !--initialize.
   nl=1
   temp=tempr(itemp)
   imax=20
   if (nl.gt.nlmax) call error('coh','too many legendre orders.',' ')
   nx=nl+1
   nj=nex+nl
   do i=1,nl
      ex(nex+i)=0
   enddo
   nlt=5
   if (ne.lt.nlt) nlt=ne
   nlt1=nlt-1
   allocate(stk(nx,imax))

   !--determine the energy grid adaptively and
   !--store the cross sections in a scratch file.
   nbragg=nl
   e=0
   call sigcoh(e,enext,s,nbragg,lat,temp,emax,natom,index)
   ix=1
   j=0
   iex=0
  100 continue
   iex=iex+1
   call finda(iex,ex,nex,iold,bufo,nbuf)
   if (ex(1).gt.enext*(1+small)) go to 105
   x(1)=ex(1)
   y(1)=ex(2)
   z(1)=ex(nex)
   j=j+1
   call loada(j,ex,nj,inew,bufn,nbuf)
   go to 100
  105 continue
   ix=ix+1
   x(ix)=ex(1)
   y(ix)=ex(2)
   z(ix)=ex(nex)
   if (ix.lt.nlt) go to 100
   ! prime stack with first bragg edge
   e=enext
   call sigcoh(e,enext,s,nl,lat,temp,emax,natom,index)
   stk(1,1)=e
   do il=1,nl
      stk(1+il,1)=s(il)
   enddo
   i=1
   ! add next bragg edge to stack
  120 continue
   e=enext
   call sigcoh(e,enext,s,nl,lat,temp,emax,natom,index)
   call upstk(e,s,nl,nx,i,stk)
   ! make sure input grid points are included
  125 continue
   do 130 ix=1,nlt
   if (x(ix).gt.stk(1,i)*(1+small)) go to 135
  130 continue
   go to 140
  135 continue
   if (x(ix).ge.stk(1,i-1)*(1-small)) go to 140
   e=x(ix)
   call sigcoh(e,en,s,nl,lat,temp,emax,natom,index)
   call upstk(e,s,nl,nx,i,stk)
   ! compare linear approximation to true function.
  140 continue
   if (i.eq.imax) go to 160
   xm=half*(stk(1,i-1)+stk(1,i))
   xm=sigfig(xm,7,0)
   if (stk(1,i-1)-stk(1,i).lt.eps*xm) go to 160
   call sigcoh(xm,en,s,nl,lat,temp,emax,natom,index)
   do 150 il=1,nl
   call terp1(stk(1,i),stk(1+il,i),&
     stk(1,i-1),stk(1+il,i-1),xm,ym,2)
   test=tol*abs(s(il))
   if (test.lt.tolmin) test=tolmin
   if (abs(s(il)-ym).gt.test) go to 210
  150 continue
   ! all components pass.  save top point in stack and continue.
  160 continue
   j=j+1
   ej(1)=stk(1,i)
  170 continue
   if (ej(1).le.x(3)*(1+small).or.iex.eq.ne) go to 190
   do ix=1,nlt1
      x(ix)=x(ix+1)
      y(ix)=y(ix+1)
      z(ix)=z(ix+1)
   enddo
   iex=iex+1
   call finda(iex,ex,nex,iold,bufo,nbuf)
   x(nlt)=ex(1)
   y(nlt)=ex(2)
   z(nlt)=ex(nex)
   if (iex.eq.ne) nlt=nlt-1
   if (iex.eq.ne) nlt1=nlt1-1
   go to 170
  190 continue
   ej(2)=terp(x,y,nlt,ej(1),nlt1)
   ej(nex)=terp(x,z,nlt,ej(1),nlt1)
   if (stk(1,i).gt.emax*(1+small)) ej(nex)=0
   do il=1,nl
      ej(nex+il)=stk(1+il,i)
      if (stk(1,i).gt.emax*(1+small)) ej(nex+il)=0
   enddo
   if (stk(1,i).gt.emax*(1+small)) go to 230
   call loada(j,ej,nj,inew,bufn,nbuf)
   i=i-1
   if (i.gt.1) go to 125
   go to 120
   ! test fails.  add point to stack and continue.
  210 continue
   call upstk(xm,s,nl,nx,i,stk)
   go to 125
   ! linearization complete.  save last point.
  230 continue
   ne=j
   j=-j
   call loada(j,ej,nj,inew,bufn,nbuf)
   isave=iold
   iold=inew
   inew=isave
   nex=nj
   ltt=7
   b(1)=za
   b(2)=awr
   b(3)=0
   b(4)=ltt ! temporary flag
   b(5)=0
   b(6)=0
   math=matdp
   mfh=6
   mth=mtref
   call contio(0,0,nscr,b,nb,nw)
   b(1)=1
   b(2)=1
   b(3)=-nbragg ! use lip field for nbragg
   b(4)=0 ! law=0 for coherent elastic data
   b(5)=1
   b(6)=2
   b(7)=2
   b(8)=2
   b(9)=1.e-5_kr
   b(10)=1
   b(11)=emax
   b(12)=1
   nw=12
   call tab1io(0,0,nscr,b,nb,nw)
   ncdse=3
   call asend(0,nscr)
   deallocate(stk)
   return
   end subroutine coh

   subroutine upstk(e,s,nl,nx,i,stk)
   !-------------------------------------------------------------------
   ! Update the linearization stack with energy e and cross
   ! sections s.  Here, i is the current index to the stack in stk,
   ! nl is the number of legendre orders in s, and nx is the
   ! cycle length in the stack.
   !-------------------------------------------------------------------
   ! externals
   integer::nl,nx,i
   real(kr)::e,s(nl),stk(nx,*)
   ! internals
   integer::j

   i=i+1
   stk(1,i)=stk(1,i-1)
   stk(1,i-1)=e
   do j=1,nl
      stk(1+j,i)=stk(1+j,i-1)
      stk(1+j,i-1)=s(j)
   enddo
   return
   end subroutine upstk

   subroutine sigcoh(e,enext,s,nl,lat,temp,emax,natom,index)
   !-------------------------------------------------------------------
   ! Compute the first nl Legendre components of the coherent scatter-
   ! ing at energy e from lattice type lat.  Here enext is the next
   ! Bragg edge.  Initialize if e=0.  A list of reciprocal lattice
   ! shells and weights is precomputed and stored for use at all e.
   ! Long, closely-spaced shells are grouped together to speed up the
   ! calculation.
   !       lat=1  graphite
   !       lat=2  be
   !       lat=3  beo
   !       lat=10 read from endf6
   ! nl returns no. of Bragg edges on initialization call.
   !-------------------------------------------------------------------
   use physics ! provides pi,amassn,amu,hbar,ev
   use util    ! provides error,sigfig
   use mathm   ! provides legndr
   ! externals
   integer::nl,lat,natom,index
   real(kr)::e,enext,s(*),temp,emax
   ! internals
   integer::nord,nw,k
   integer::i1m,i1,l1,i2m,i2,l2,i3m,i3,l3,nr,np
   integer::l,i,imax,jmin,j,il,lmax,last
   real(kr)::amne,econ,tsqx,a,c,amsc,scoh,wal2,wint,x
   real(kr)::w1,w2,w3,tsq,tau,w,f,st,sf,blast,re
   real(kr)::t2,ulim,phi,elim,u,twopis,c1,c2,recon,scon
   real(kr)::p(6)
   real(kr),dimension(:),allocatable::wrk
   integer,parameter::nd=10
   real(kr),dimension(10),parameter::dwf1=(/&
     2.1997e0_kr,2.7448e0_kr,3.2912e0_kr,3.8510e0_kr,4.4210e0_kr,&
     4.9969e0_kr,6.1624e0_kr,7.3387e0_kr,9.6287e0_kr,11.992e0_kr/)
   real(kr),dimension(10),parameter::dwf2=(/&
     3.16663e0_kr,3.88842e0_kr,4.62944e0_kr,5.40517e0_kr,6.19880e0_kr,&
     7.0042e0_kr,8.63665e0_kr,10.2865e0_kr,0.e0_kr,0.e0_kr/)
   real(kr),dimension(10),parameter::dwf3=(/&
     2.153e0_kr,2.6374e0_kr,3.1348e0_kr,3.6513e0_kr,4.1798e0_kr,&
     4.7164e0_kr,5.8052e0_kr,6.9068e0_kr,0.e0_kr,0.e0_kr/)
   real(kr),dimension(10),parameter::tmp=(/&
     296.e0_kr,400.e0_kr,500.e0_kr,600.e0_kr,700.e0_kr,800.e0_kr,&
     1000.e0_kr,1200.e0_kr,1600.e0_kr,2000.e0_kr/)
   real(kr),parameter::gr1=2.4573e-8_kr
   real(kr),parameter::gr2=6.700e-8_kr
   real(kr),parameter::gr3=12.011e0_kr
   real(kr),parameter::gr4=5.50e0_kr
   real(kr),parameter::be1=2.2856e-8_kr
   real(kr),parameter::be2=3.5832e-8_kr
   real(kr),parameter::be3=9.01e0_kr
   real(kr),parameter::be4=7.53e0_kr
   real(kr),parameter::beo1=2.695e-8_kr
   real(kr),parameter::beo2=4.39e-8_kr
   real(kr),parameter::beo3=12.5e0_kr
   real(kr),parameter::beo4=1.e0_kr
   real(kr),parameter::sqrt3=1.732050808e0_kr
   real(kr),parameter::cw=0.658173e-15_kr
   real(kr),parameter::half=0.5e0_kr
   real(kr),parameter::eps=0.05e0_kr
   real(kr),parameter::zero=0
   save k,recon,scon

   !--initialize.
   if (e.gt.zero) go to 210
   twopis=(2*pi)**2
   amne=amassn*amu
   econ=ev*8*(amne/hbar)/hbar
   recon=1/econ
   tsqx=econ/20
   nord=2
   if (lat.eq.10) go to 200
   if (lat.eq.1) then
      ! graphite constants
      a=gr1
      c=gr2
      amsc=gr3
      scoh=gr4/natom
      wal2=terp(tmp,dwf1,nd,temp,nord)
   else if (lat.eq.2) then
      ! beryllium constants
      a=be1
      c=be2
      amsc=be3
      scoh=be4/natom
      wal2=terp(tmp,dwf2,nd,temp,nord)
   else if (lat.eq.3) then
      ! beryllium oxide constants
      a=beo1
      c=beo2
      amsc=beo3
      scoh=beo4/natom
      wal2=terp(tmp,dwf3,nd,temp,nord)
   else
      call error('sigcoh','illegal lat.',' ')
   endif
   c1=4/(3*a*a)
   c2=1/(c*c)
   scon=scoh*(4*pi)**2/(2*a*a*c*sqrt3*econ)
   wint=cw*amsc*wal2
   t2=hbar/(2*amu*amsc)
   ulim=econ*emax
   nw=1000000
   allocate(wrk(nw))

   !--compute and sort lattice factors.
   phi=ulim/twopis
   i1m=int(a*sqrt(phi))
   i1m=i1m+1
   k=0
   do 185 i1=1,i1m
   l1=i1-1
   i2m=int(half*(l1+sqrt(3*(a*a*phi-l1*l1))))
   i2m=i2m+1
   do 180 i2=i1,i2m
   l2=i2-1
   x=phi-c1*(l1*l1+l2*l2-l1*l2)
   i3m=0
   if (x.gt.zero) i3m=int(c*sqrt(x))
   i3m=i3m+1
   do 175 i3=1,i3m
   l3=i3-1
   w1=2
   if (l1.eq.l2) w1=1
   w2=2
   if (l1.eq.0.or.l2.eq.0) w2=1
   if (l1.eq.0.and.l2.eq.0) w2=half
   w3=2
   if (l3.eq.0) w3=1
   tsq=tausq(l1,l2,l3,c1,c2,twopis)
   if (tsq.le.zero.or.tsq.gt.ulim) go to 160
   tau=sqrt(tsq)
   w=exp(-tsq*t2*wint)*w1*w2*w3/tau
   f=w*form(lat,l1,l2,l3)
   if (k.gt.0.and.tsq.gt.tsqx) go to 150
   k=k+1
   if ((2*k).gt.nw) call error('sigcoh','storage exceeded.',' ')
   wrk(2*k-1)=tsq
   wrk(2*k)=f
   go to 160
  150 continue
   do 155 i=1,k
   if (tsq.lt.wrk(2*i-1).or.tsq.ge.(1+eps)*wrk(2*i-1)) go to 155
   wrk(2*i)=wrk(2*i)+f
   go to 160
  155 continue
   k=k+1
   if ((2*k).gt.nw) call error('sigcoh','storage exceeded.',' ')
   wrk(2*k-1)=tsq
   wrk(2*k)=f
  160 continue
   tsq=tausq(l1,-l2,l3,c1,c2,twopis)
   if (tsq.le.zero.or.tsq.gt.ulim) go to 175
   tau=sqrt(tsq)
   w=exp(-tsq*t2*wint)*w1*w2*w3/tau
   f=w*form(lat,l1,-l2,l3)
   if (k.gt.0.and.tsq.gt.tsqx) go to 165
   k=k+1
   if ((2*k).gt.nw) call error('sigcoh','storage exceeded.',' ')
   wrk(2*k-1)=tsq
   wrk(2*k)=f
   go to 175
  165 continue
   do 170 i=1,k
   if (tsq.lt.wrk(2*i-1).or.tsq.ge.(1+eps)*wrk(2*i-1)) go to 170
   wrk(2*i)=wrk(2*i)+f
   go to 175
  170 continue
   k=k+1
   if ((2*k).gt.nw) call error('sigcoh','storage exceeded.',' ')
   wrk(2*k-1)=tsq
   wrk(2*k)=f
  175 continue
  180 continue
  185 continue
   imax=k-1
   do i=1,imax
      jmin=i+1
      do j=jmin,k
         if (wrk(2*j-1).lt.wrk(2*i-1)) then
            st=wrk(2*i-1)
            sf=wrk(2*i)
            wrk(2*i-1)=wrk(2*j-1)
            wrk(2*i)=wrk(2*j)
            wrk(2*j-1)=st
            wrk(2*j)=sf
         endif
      enddo
   enddo
   k=k+1
   wrk(2*k-1)=ulim
   wrk(2*k)=wrk(2*k-2)
   nw=2*k
   enext=recon*wrk(1)
   enext=sigfig(enext,7,-1)
   ! copy data to global fl array
   allocate(fl(nw))
   do i=1,nw
      fl(i)=wrk(i)
   enddo
   deallocate(wrk)
   nl=k
   return

   !--bragg parameters already read from endf6
  200 continue
   !--fl(index) is the coherent S(E) table for the current temperature from ENDF6
   nr=int(fl(index+4))  ! number of interpolation regions
   np=int(fl(index+5))  ! number of points
   k=np
   nl=np
   blast=0
   scon=1
   enext=fl(index+2*nr+6)
   enext=sigfig(enext,7,-1)
   do i=1,nl
      l=index+2*(i-1)
      fl(l)=fl(l+2*nr+6)*econ
      fl(l+1)=fl(l+2*nr+7)-blast
      blast=fl(l+2*nr+7)
   enddo
   return

   !--compute cross sections at this energy
  210 continue
   re=1/e
   do il=1,nl
      s(il)=0
   enddo
   last=0
   do i=1,k
      tsq=fl(2*i-1)
      elim=tsq*recon
      if (elim.ge.e) exit
      f=fl(2*i)
      if (e.gt.emax) f=0
      u=1-2*elim*re
      lmax=nl-1
      call legndr(u,p,lmax)
      do il=1,nl
         s(il)=s(il)+f*p(il)
      enddo
      if (i.eq.k) last=1
   enddo
   do il=1,nl
      s(il)=s(il)*scon*re
   enddo
   if (last.eq.1) elim=emax
   if (elim.gt.emax) elim=emax
   enext=sigfig(elim,7,-1)
   if (e.gt.sigfig(enext,7,-1)) enext=sigfig(elim,7,+1)
   return

   contains

      real(kr) function tausq(m1,m2,m3,c1,c2,twopis)
      integer::m1,m2,m3
      real(kr)::c1,c2,twopis
      tausq=(c1*(m1*m1+m2*m2+m1*m2)+(m3*m3*c2))*twopis
      return
      end function tausq

   end subroutine sigcoh

   real(kr) function form(lat,l1,l2,l3)
   !-------------------------------------------------------------------
   ! Compute form factors for the specified lattice.
   !       lat=1  graphite
   !       lat=2  be
   !       lat=3  beo
   !-------------------------------------------------------------------
   use physics ! provides pi
   ! externals
   integer::lat,l1,l2,l3
   ! internals
   integer::i
   real(kr),parameter::beo1=7.54e0_kr
   real(kr),parameter::beo2=4.24e0_kr
   real(kr),parameter::beo3=11.31e0_kr

   if (lat.eq.1) then
      ! graphite
      i=l3/2
      if ((2*i).ne.l3) then
         form=sin(pi*(l1-l2)/3)**2
      else
         form=(6+10*cos(2*pi*(l1-l2)/3))/4
      endif
   else if (lat.eq.2) then
      ! beryllium
      form=1+cos(2*pi*(2*l1+4*l2+3*l3)/6)
   else if (lat.eq.3) then
      ! beryllium oxide
      form=(1+cos(2*pi*(2*l1+4*l2+3*l3)/6))*(beo1+beo2+beo3*&
        cos(3*pi*l3/4))
   endif
   return
   end function form

   subroutine iel(mat,itemp,iold,inew,ne,nex,mtref,index)
   !-------------------------------------------------------------------
   ! Compute the elastic scattering from polyethylene or hydrogen
   ! in zirconium hydride using the incoherent approximation.
   ! The calcem energy grid is used.  If mat.eq.11 or 12, built-in
   ! parameters from the ENDF/B-III evaluations are used.
   ! If mat.eq.20, material parameters have been read from ENDF6.
   !-------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso,nsyse
   use endf   ! provides endf routines and variables
   use util   ! provides error,finda,loada
   ! externals
   integer::mat,itemp,iold,inew,ne,nex,mtref,index
   ! internals
   integer::idis,iex,iet,iu,ix,nj,nr,np,ip,ir,ltt,nb,nw
   integer::n,nup,nup1,isave,nne
   real(kr)::dwa,c1,e,c2,u,rc2,temp,tt1,ttn,tnxt
   real(kr)::x1,r1x1,xsec,x2,unow
   real(kr)::ex(20),ej(20)
   real(kr),dimension(:),allocatable::xie
   integer,parameter::nupmax=10
   real(kr),dimension(8),parameter::tmp=(/296e0_kr,400.e0_kr,&
     500.e0_kr,600.e0_kr,700.e0_kr,800.e0_kr,1000.e0_kr,1200.e0_kr/)
   real(kr),dimension(8),parameter::dwh=(/8.4795e0_kr,9.0854e0_kr,&
     9.8196e0_kr,10.676e0_kr,11.625e0_kr,12.643e0_kr,14.822e0_kr,&
     17.125e0_kr/)
   real(kr),dimension(8),parameter::dwz=(/1.9957e0_kr,2.6546e0_kr,&
     3.2946e0_kr,3.9380e0_kr,4.5835e0_kr,5.2302e0_kr,6.5260e0_kr,&
     7.8236e0_kr/)
   real(kr),parameter::c11a=162.88e0_kr
   real(kr),parameter::c11b=296.e0_kr
   real(kr),parameter::c11c=34.957e0_kr
   real(kr),parameter::c11d=350.e0_kr
   real(kr),parameter::c11e=40.282e0_kr
   real(kr),parameter::c12a=81.44e0_kr
   real(kr),parameter::c13a=6.3366e0_kr
   real(kr),parameter::up=1.1e0_kr
   real(kr),parameter::dn=.9e0_kr

   !--initialize
   if (iprint.eq.2) write(nsyso,'(/'' iel''/)')
   temp=tempr(itemp)
   nj=nex+1

   !--select material parameters
   if (mat.eq.11) then
      sb=c11a
      dwa=c11c+(temp-c11b)*(c11e-c11c)/(c11d-c11b)
   else if (mat.eq.12) then
      dwa=terp(tmp,dwh,8,temp,3)
      sb=c12a
   else if (mat.eq.13) then
      dwa=terp(tmp,dwz,8,temp,3)
      sb=c13a
   else if (mat.eq.20) then
      !--fl(index) is the incoherent W'(T) table from ENDF6
      sb=fl(index)
      nr=int(fl(index+4))  ! number of interpolation regions
      np=int(fl(index+5))  ! number of points
      if (np.eq.1) then
         tt1=fl(index+2*nr+6)
         if (abs(temp-tt1).gt.temp/10) call error('iel',&
           'bad temperature for debye-waller factor',' ')
         dwa=fl(index+2*nr+7)
      else
         tt1=fl(index+2*nr+6)
         ttn=fl(index+2*nr+2*np+4)
         if (temp.lt.dn*tt1.or.temp.gt.up*ttn) call error('iel',&
           'bad temperature for debye-waller factor',' ')
         if (tt1.gt.temp) fl(index+2*nr+6)=temp
         if (ttn.lt.temp) fl(index+2*nr+2*np+4)=temp
         ip=2
         ir=1
         call terpa(dwa,temp,tnxt,idis,fl(index),ip,ir)
      endif
   else
      call error('iel','unknown material identifier.',' ')
   endif
   c1=sb/(2*natom)

   !--check on calcem energy grid
   nne=0
   do
     if (esi(nne+1).eq.0) exit
     nne=nne+1
   enddo
   allocate(xie(nne))

   !--write head and tab2 records for mf6
   !--in lanl format
   ltt=6
   math=matdp
   mfh=6
   mth=mtref
   scr(1)=za
   scr(2)=awr
   scr(3)=0
   scr(4)=ltt ! temporary flag for incoherent inelastic data
   scr(5)=0
   scr(6)=0
   nw=6
   call contio(0,0,nscr,scr,nb,nw)
   scr(1)=1
   scr(2)=1
   scr(3)=-2 ! special flag for incoherent inelastic data
   scr(4)=1 ! law=1 for incoherent inelastic data
   scr(5)=1
   scr(6)=2
   scr(7)=2
   scr(8)=2
   scr(9)=1.e-5_kr
   scr(10)=1
   scr(11)=emax
   scr(12)=1
   nw=12
   call tab1io(0,0,nscr,scr,nb,nw)
   scr(1)=temp
   scr(2)=0
   scr(3)=3 ! LANG=3 is a special flag for equally probable cosines
   scr(5)=1
   scr(6)=nne
   scr(7)=nne
   scr(8)=2
   nw=8
   call tab2io(0,0,nscr,scr,nb,nw)
   n=nbin
   nup=n
   if (nup.gt.nupmax) nup=nupmax
   nup1=nup+1
   nw=n+8

   !--compute cross sections on calcem energy grid
   !--and equi-probable angles.
   do iex=1,nne
      e=esi(iex)
      c2=2*e*dwa
      u=-1
      rc2=1/c2
      x1=exp(-2*c2)
      r1x1=1/(1-x1)
      xsec=c1*rc2*(1-x1)
      scr(1)=0
      scr(2)=e
      scr(3)=0
      scr(4)=0
      scr(5)=n+2
      scr(6)=n+2
      scr(7)=e
      scr(8)=1
      do iu=1,n
         x2=exp(-c2*(1-u))
         unow=1+rc2*log((1-x1)/n+x2)
         scr(8+iu)=n*rc2*(exp(-c2*(1-unow))*(c2*unow-1)-&
           x2*(c2*u-1))*r1x1
         u=unow
      enddo
      xie(iex)=xsec
      if (iprint.eq.2) write(nsyso,'(4x,1p,2e12.4,0p,10f9.4)')&
        e,xsec,(scr(8+iu),iu=1,nup)
      if (iprint.eq.2.and.n.gt.nup)&
        write(nsyso,'(28x,10f9.4)') (scr(8+iu),iu=nup1,n)
      call listio(0,0,nscr,scr,nb,nw)
   enddo
   do iex=1,ne
      call finda(iex,ex,nex,iold,bufo,nbuf)
      do ix=1,nex
         ej(ix)=ex(ix)
      enddo
      ej(nj)=terp(esi,xie,nne,ex(1),3)
      if (iex.eq.ne) ej(nj)=0
      iet=iex
      if (iex.eq.ne) iet=-iet
      call loada(iet,ej,nj,inew,bufn,nbuf)
   enddo
   isave=iold
   iold=inew
   inew=isave
   nex=nj
   call asend(0,nscr)
   ncdse=5+ne*((n+11)/6)
   return
   end subroutine iel

   real(kr) function terp(x,y,nl,arg,il1)
   !-------------------------------------------------------------------
   ! This function does Lagrangian interpolation or
   ! extrapolation of ilth order on x and y for arg
   ! when the tables are either increasing or decreasing
   !      x          x array, independent variable
   !      y          y array, dependent variable
   !      nl         number of entries in tables of x and y
   !      arg        independent variable value
   !      il         order of interpolation
   !-------------------------------------------------------------------
   ! externals
   integer::nl,il1
   real(kr)::x(nl),y(nl),arg
   ! internals
   integer::il,l,iadd,il2,ilow,ihi,iusel,iuseh,ibeg,iend
   integer::last,n,m,i,in,ip,inp
   real(kr)::sum,p,pk
   real(kr),parameter::small=1.e-10_kr

   il=il1
   if (nl.gt.il) go to 120
   if (nl.eq.il) go to 110
   ! not enough entries in tables for this order interpolation
   il=nl
  110 continue
   l=1
   go to 260
  120 continue
   il2=il/2
   iadd=mod(il,2)
   ! check if tables in increasing or decreasing sequence
   if (x(1).le.x(nl)) then
      ! increasing sequence
      ilow=il2+1
      ihi=nl-il2-iadd
      iusel=1
      iuseh=nl-il+1
      ibeg=ilow+1
      iend=ihi-1
      last=iend-il2+1
      iadd=0
    else
       ! decreasing sequence
      ilow=nl-il2
      ihi=il2+iadd+1
      iusel=nl-il+1
      iuseh=1
      ibeg=ihi+1
      iend=ilow-1
      last=2
      iadd=1-iadd
   endif
   ! checks if arg is smaller than table values
   if (abs(arg-x(ilow)).lt.small*arg) go to 160
   if (arg.gt.x(ilow)) go to 170
   ! smaller than smallest table value
   l=iusel
   go to 260
  160 continue
   terp=y(ilow)
   go to 300
   ! checks if arg is greater than table values
  170 continue
    if (abs(x(ihi)-arg).lt.small*arg) go to 190
   if (x(ihi).gt.arg) go to 200
   ! arg greater than table value
   l=iuseh
   go to 260
  190 continue
   terp=y(ihi)
   go to 300
   ! searches x array to bracket arg
  200 continue
   do 230 n=ibeg,iend
   if (iusel.gt.1) go to 210
   m=n
   go to 220
  210 continue
   m=nl-n+1
  220 continue
   if (abs(x(m)-arg).lt.small*arg) go to 240
   if (x(m).gt.arg) go to 250
  230 continue
   l=last
   go to 260
   ! equals argument, return ok
  240 continue
   terp=y(m)
   go to 300
   ! eureka
  250 continue
   l=m-il2+iadd
  260 continue
   ! interpolation section
   sum=0
   do i=1,il
      p=1
      pk=1
      in=l+i-1
      do ip=1,il
         if (ip.ne.i) then
            inp=l+ip-1
            p=p*(arg-x(inp))
            pk=pk*(x(in)-x(inp))
         endif
      enddo
      sum=sum+p*y(in)/pk
   enddo
   terp=sum
  300 continue
   return
   end function terp

   subroutine calcem(temp,itemp,iold,inew,ne,nex,mtref)
   !-------------------------------------------------------------------
   ! Calculate incoherent inelastic scattering kernels from
   ! s(alpha,beta) in endf mf7 format or from analytic models
   ! of s(alpha,beta). The incident energy grid is fixed (see egrid).
   ! For iform=0, the secondary energy and scattering cosine grids
   ! are determined adaptively for linear interpolations.  A compact
   ! angle representation using equally probable cosines is generated.
   ! For iform=1, a fixed grid of mu values is used, and the secondary
   ! distribution is determined adaptively for linear interpolation
   ! for each mu.  In both cases, the results are written to a
   ! scratch tape for tpend.
   !-------------------------------------------------------------------
   use mainio  ! provides nsysi,nsyso,nsyse
   use endf    ! provides endf routines and variables
   use physics ! provides bk
   use util    ! provides error,repoz,sigfig
   use mathm   ! provides legndr
   ! externals
   integer::itemp,iold,inew,ne,nex,mtref
   real(kr)::temp
   ! internals
   character(len=60)::strng
   integer::nr,np,nwtab,nl,nlt,nlp,nlp1,nnl,jmax,nne
   integer::i,ia,nl1,ltt,loc,l,jscr,ilog
   integer::matd,itprnt,nb,nw,ni,nbeta,lt,it,nalpha
   integer::itrunc,ib,ip,ir,idis,ie,nmu,nep,istart,iend
   integer::jbeta,j,iskip,il,k,jnz,ll,jj,ii,nexn,ien,isave
    integer::nll
   real(kr)::smz,t1,t,tmax,test,tempt,tt1,ttn,tnxt,xm
   real(kr)::uu,uum,ym,test2,xlast,ulast,xs
   real(kr)::b,diff,enow,ep,sabmin,tev,ylast
   real(kr)::u,xl,yl,sum
   real(kr)::tone,elo
   integer,parameter::ngrid=118
   integer,parameter::nlmax=65
   integer,parameter::nemax=5000
   integer,parameter::mumax=300
   integer,parameter::imax=20
   real(kr)::ex(imax),x(imax),y(nlmax,imax),yt(nlmax)
   real(kr)::yy(imax),yu(2*nemax)
   real(kr)::ubar(ngrid)
   real(kr)::u2,u2last,u3,u3last,p2(ngrid),p3(ngrid),p(4)
   real(kr)::uj(mumax),sj(mumax)
   real(kr),dimension(:),allocatable::alpha,beta
   real(kr),dimension(:,:),allocatable::sab
   real(kr),dimension(ngrid),parameter::egrid=(/&
     1.e-5_kr,1.78e-5_kr,2.5e-5_kr,3.5e-5_kr,5.0e-5_kr,7.0e-5_kr,1.e-4_kr,&
     1.26e-4_kr,1.6e-4_kr,2.0e-4_kr,.000253e0_kr,.000297e0_kr,.000350e0_kr,&
     .00042e0_kr,.000506e0_kr,.000615e0_kr,.00075e0_kr,.00087e0_kr,&
     .001012e0_kr,.00123e0_kr,.0015e0_kr,.0018e0_kr,.00203e0_kr,.002277e0_kr,&
     .0026e0_kr,.003e0_kr,.0035e0_kr,.004048e0_kr,.0045e0_kr,.005e0_kr,&
     .0056e0_kr,.006325e0_kr,.0072e0_kr,.0081e0_kr,.009108e0_kr,.01e0_kr,&
     .01063e0_kr,.0115e0_kr,.012397e0_kr,.0133e0_kr,.01417e0_kr,.015e0_kr,&
     .016192e0_kr,.0182e0_kr,.0199e0_kr,.020493e0_kr,.0215e0_kr,.0228e0_kr,&
     .0253e0_kr,.028e0_kr,.030613e0_kr,.0338e0_kr,.0365e0_kr,.0395e0_kr,&
     .042757e0_kr,.0465e0_kr,.050e0_kr,.056925e0_kr,.0625e0_kr,.069e0_kr,&
     .075e0_kr,.081972e0_kr,.09e0_kr,.096e0_kr,.1035e0_kr,.111573e0_kr,&
     .120e0_kr,.128e0_kr,.1355e0_kr,.145728e0_kr,.160e0_kr,.172e0_kr,&
     .184437e0_kr,.20e0_kr,.2277e0_kr,.2510392e0_kr,.2705304e0_kr,&
     .2907501e0_kr,.3011332e0_kr,.3206421e0_kr,.3576813e0_kr,.39e0_kr,&
     .4170351e0_kr,.45e0_kr,.5032575e0_kr,.56e0_kr,.625e0_kr,&
     .70e0_kr,.78e0_kr,.86e0_kr,.95e0_kr,1.05e0_kr,1.16e0_kr,1.28e0_kr,&
     1.42e0_kr,1.55e0_kr,1.70e0_kr,1.855e0_kr,2.02e0_kr,2.18e0_kr,&
     2.36e0_kr,2.59e0_kr,2.855e0_kr,3.12e0_kr,3.42e0_kr,3.75e0_kr,&
     4.07e0_kr,4.46e0_kr,4.90e0_kr,5.35e0_kr,5.85e0_kr,6.40e0_kr,&
     7.00e0_kr,7.65e0_kr,8.40e0_kr,9.15e0_kr,9.85e0_kr,10.00e0_kr/)
   real(kr),parameter::unity=1.0e0_kr
   real(kr),parameter::sabflg=-225.e0_kr
   real(kr),parameter::eps=1.e-4_kr
   real(kr),parameter::tolmin=5.e-7_kr
   real(kr),parameter::half=0.5e0_kr
   real(kr),parameter::therm=.0253e0_kr
   real(kr),parameter::break=3000.e0_kr
   real(kr),parameter::em9=1.e-9_kr
   real(kr),parameter::zero=0.e0_kr
   real(kr),parameter::tenth=0.1e0_kr
   real(kr),parameter::up=1.1e0_kr
   real(kr),parameter::dn=0.9e0_kr
   real(kr),parameter::uumin=0.00001e0_kr
   real(kr),parameter::yumin=2.e-7_kr
   integer,parameter::nlpmx=10
   save nwtab,sabmin,nl,nlt,nlp,nlp1,nl1,nnl,jmax,nne
   tevz=therm

   !--initialize if itemp=1.
   if (itemp.eq.1) then
      sabmin=exp(sabflg)
      nwtab=8
      nl=nbin+1
      if (nl.gt.nlmax)&
        call error('calcem','nl too large for binning.',' ')
      nlp=nl+1
      if (nlp.gt.nlpmx) nlp=nlpmx
      nlp1=nlp+1
      nl1=nl+1
      nnl=-nl
      jmax=(nwscr-50)/nl
      nne=0
      do i=1,ngrid
         if (egrid(i).gt.emax.and.nne.eq.0) nne=i
      enddo
      if (nne.eq.0) nne=ngrid
      allocate(esi(nne+1))
      allocate(xsi(nne+1))
      nlt=5
      if (nne.le.nlt) nlt=nne-1
   endif
   call repoz(nscr)
   nsc=0
   ncds=0
   if (iinc.eq.1) go to 200
   matd=matde
   itprnt=0

   !--bound atom scattering.
   !--read s(alpha,beta) from endf tape.
   !--initialize for desired material
   call findf(matd,7,4,nendf)
   call contio(nendf,0,0,scr,nb,nw)
   lasym=n1h

   ! lasym= 0 or 1 = traditional endf definitions.
   ! lasym= 2 or 3 = traditional lasym + 2
   !               = leapr's isabt=1 option was used.
   if (lasym.gt.1) call error ('calcem','isabt=1 pendf tape found',&
                         'thermr cannot process this format')

   lat=l2h
   call listio(nendf,0,0,scr,nb,nw)
   ilog=l1h
   ni=n1h
   ! parameters for principal atom
   smz=scr(7)/natom
   az=scr(9)
   sb=smz*((az+1)/az)**2
   ! sct parameters for other atoms
   if (iverf.ge.6) then
      sb2=0
      az2=0
      if (ni.gt.6) then
         if (scr(13).eq.zero) az2=scr(15)
         if (az2.ne.zero) sb2=scr(14)*((az2+1)/az2)**2
         sb2=sb2/scr(18)/natom
         if (ni.gt.12) then
            if (scr(19).eq.zero)&
              call error('calcem','only 2 sct atoms allowed.',' ')
         endif
      endif
   endif
   call tab2io(nendf,0,0,scr,nb,nw)
   nbeta=nint(scr(6))
   loc=1+nw
   ! read first tab1 record for alpha grid.
   call tab1io(nendf,0,0,scr(loc),nb,nw)
   t1=c1h
   b=c2h
   lt=l1h
   it=0
   nalpha=n2h
   l=loc+nw
   do while (nb.ne.0)
      call moreio(nendf,0,0,scr(l),nb,nw)
      l=l+nw
   enddo
   allocate(alpha(nalpha))
   l=loc+4+2*n1h
   do i=1,nalpha
      alpha(i)=scr(2*i+l)
   enddo
   allocate(beta(nbeta))
   allocate(sab(nalpha,nbeta))

   !--read s(alpha,beta) data into memory
   itrunc=0
   ib=1
   beta(ib)=b
  120 continue
   t=scr(loc)
   do while (abs(t-temp).gt.temp/500)
      if (t.gt.temp)&
        call error('calcem','desired temperature not found.',' ')
      it=it+1
      if (it.gt.lt)&
        call error('calcem','desired temperature not found.',' ')
      call listio(nendf,0,0,scr(loc),nb,nw)
      l=loc+nw
      do while (nb.ne.0)
         call moreio(nendf,0,0,scr(l),nb,nw)
         l=l+nw
      enddo
      t=scr(loc)
   enddo
   if (abs(t-temp).gt.eps.and.itprnt.eq.0) then
      itprnt=1
      diff=abs(t-temp)
      write(nsyso,'(/&
        &'' difference between temperatures desired and found is '',&
        &1p,e10.2)') diff
      write(nsyse,'(/&
        & '' difference between temperatures desired and found is '',&
        &1p,e10.2)') diff
   endif
   l=loc+6
   if (t.eq.t1) l=l+2*nint(scr(loc+4))+1
   do ia=1,nalpha
      if (ilog.eq.0) then
         if (scr(l).gt.sabmin) sab(ia,ib)=log(scr(l))
         if (scr(l).le.sabmin) then
            if (scr(l).gt.zero) itrunc=1
            sab(ia,ib)=sabflg
         endif
      else
         sab(ia,ib)=scr(l)
      endif
      if (t.eq.t1) l=l+1
      l=l+1
   enddo
   do while (it.lt.lt)
      it=it+1
      call listio(nendf,0,0,scr(loc),nb,nw)
      l=loc
      do while (nb.ne.0)
         l=l+nw
         call moreio(nendf,0,0,scr(l),nb,nw)
      enddo
   enddo
   if (ib.eq.nbeta) go to 180
   ib=ib+1
   call tab1io(nendf,0,0,scr(loc),nb,nw)
   l=loc
   do while (nb.ne.0)
      l=l+nw
      call moreio(nendf,0,0,scr(l),nb,nw)
   enddo
   b=scr(loc+1)
   beta(ib)=b
   it=0
   go to 120
  180 continue
   if (itrunc.eq.1) then
      tmax=2*(-sabflg-5)*bk*temp
      write(nsyso,'(/&
        &'' ***warning***'',&
        &''sab contains numbers that are too small'',&
        &'' to represent on this machine.''/&
        &'' results may be bad for transfers larger than '',&
        &f8.3,'' ev.'')') tmax
      write(nsyse,'(/&
        &'' ***warning***'',&
        &''sab contains numbers that are too small'',&
        &'' to represent on this machine.''/&
        &'' results may be bad for transfers larger than '',&
        &f8.3,'' ev.'')') tmax
   endif
   tmax=beta(nbeta)*bk*temp
   if (lat.eq.1) tmax=tmax*tevz/(bk*temp)
   test=5
   if (tmax.lt.test) then
      write(nsyso,'(/&
        &'' ***warning***'',&
        &''maximum value of beta limits the allowed energy transfer''/&
        &'' the sct approx. will be used for transfers larger than '',&
        &f6.3,'' ev.'')') tmax
      write(nsyse,'(/&
        &'' ***warning***'',&
        &''maximum value of beta limits the allowed energy transfer''/&
        &'' the sct approx. will be used for transfers larger than '',&
        &f6.3,'' ev.'')') tmax
   endif

   !--read effective temperatures for endf6
   if (iverf.lt.6) go to 300
   call tab1io(nendf,0,0,scr,nb,nw)
   nr=nint(scr(5))
   np=nint(scr(6))
   if (np.eq.1) then
      tempt=temp
      tt1=scr(7+2*nr)
      if (abs(tempt-tt1).gt.tempt/10) call error('calcem',&
        'bad temperature for teff',' ')
      teff=scr(8+2*nr)
   else
      ip=2
      ir=1
      tempt=temp
      tt1=scr(7+2*nr)
      ttn=scr(5+2*nr+2*np)
      if (tempt.lt.dn*tt1.or.tempt.gt.up*ttn) call error('calcem',&
        'bad temperature for teff',' ')
      if (tt1.gt.tempt) scr(7+2*nr)=tempt
      if (ttn.lt.tempt) scr(5+2*nr+2*np)=tempt
      call terpa(teff,tempt,tnxt,idis,scr,ip,ir)
   endif
   if (sb2.eq.zero) go to 300
   call tab1io(nendf,0,0,scr,nb,nw)
   nr=nint(scr(5))
   np=nint(scr(6))
   if (np.gt.1) go to 195
   tt1=scr(7+2*nr)
   if (abs(tempt-tt1).gt.tempt/10) call error('calcem',&
     'bad temperature for teff2',' ')
   teff2=scr(8+2*nr)
   go to 300
  195 ip=2
   ir=1
   tempt=temp
   tt1=scr(7+2*nr)
   ttn=scr(5+2*nr+2*np)
   if (tempt.lt.dn*tt1.or.tempt.gt.up*ttn) call error('calcem',&
     'bad temperature for teff2',' ')
   if (tt1.gt.tempt) scr(7+2*nr)=tempt
   if (ttn.lt.tempt) scr(5+2*nr+2*np)=tempt
   call terpa(teff2,tempt,tnxt,idis,scr,ip,ir)
   go to 300

   !--use free-gas s(alpha,beta).
  200 continue
   ! the beta grid is fixed to 45 beta values (previously only the first 9)
   ! following an issue with Fe56 and H1 from ENDF/B-VIII.0 - a better solution
   ! would be to contruct the grid adaptatively based on emax and tempr
   nalpha=1
   nbeta=45
   allocate(alpha(nalpha))
   allocate(beta(nbeta))
   allocate(sab(nalpha,nbeta))
   beta(1)=0
   beta(2)=tenth
   beta(3)=2
   beta(4)=4
   beta(5)=6
   beta(6)=8
   beta(7)=10
   beta(8)=15
   beta(9)=25
   beta(10)=30
   beta(11)=35
   beta(12)=40
   beta(13)=45
   beta(14)=50
   beta(15)=55
   beta(16)=60
   beta(17)=65
   beta(18)=70
   beta(19)=75
   beta(20)=80
   beta(21)=100
   beta(22)=120
   beta(23)=140
   beta(24)=160
   beta(25)=180
   beta(26)=200
   beta(27)=250
   beta(28)=300
   beta(29)=350
   beta(30)=400
   beta(31)=500
   beta(32)=600
   beta(33)=700
   beta(34)=800
   beta(35)=900
   beta(36)=1000
   beta(37)=1250
   beta(38)=1500
   beta(39)=1750
   beta(40)=2000
   beta(41)=2250
   beta(42)=2500
   beta(43)=2750
   beta(44)=3000
   beta(45)=3500
   t=temp
   smz=1
   sb=smz*((az+1)/az)**2
   lat=0

   !--compute kernel and write in special mf6 energy-angle format
  300 continue
   if (iform.eq.1) go to 510
   ltt=5
   math=matdp
   mfh=6
   mth=mtref
   tev=t*bk
   teff=teff*bk
   teff2=teff2*bk
   scr(1)=za
   scr(2)=awr
   scr(3)=0
   scr(4)=ltt ! temporary flag for this format
   scr(5)=1
   scr(6)=0
   nw=6
   call contio(0,0,nscr,scr,nb,nw)
   ncds=ncds+1
   scr(1)=1
   scr(2)=1
   scr(3)=-1 ! LIP=-1 indicates incoherent inelastic data
   scr(4)=1 ! LAW=1 for incoherent inelastic data
   scr(5)=1
   scr(6)=2
   scr(7)=2
   scr(8)=2
   scr(9)=1.e-5_kr
   scr(10)=1
   scr(11)=emax
   scr(12)=1
   nw=12
   call tab1io(0,0,nscr,scr,nb,nw)
   ncds=ncds+2
   scr(1)=temp
   scr(2)=0
   scr(3)=3 ! lang=3 is a special code for equally probable cosines
   scr(4)=1
   scr(5)=1
   scr(6)=nne
   scr(7)=nne
   scr(8)=2
   nw=8
   call tab2io(0,0,nscr,scr,nb,nwtab)
   ncds=ncds+2
   cliq=0
   if (iinc.eq.1) go to 305
   if (sab(1,1).gt.sab(2,1))&
     cliq=(sab(1,1)-sab(1,2))*alpha(1)/beta(2)**2

   !--loop over given incident energy grid.
  305 continue
   ie=0
  310 continue
   ie=ie+1
   enow=egrid(ie)
   if (temp.gt.break) then
      tone=therm/bk
      elo=egrid(1)
      enow=elo*exp(log(enow/elo)*log((temp/tone)*egrid(ngrid)/elo)&
        /log(egrid(ngrid)/elo))
   endif
   enow=sigfig(enow,8,0)
   esi(ie)=enow
   xsi(ie)=0
   ubar(ie)=0
   p2(ie)=0
   p3(ie)=0
   ep=0
   x(1)=ep
   call sigl(enow,ep,nnl,tev,nalpha,alpha,nbeta,beta,&
     sab,yt,nlmax,tol)
   do il=1,nl
      y(il,1)=yt(il)
   enddo
   jbeta=-nbeta
   if (lasym.gt.0) jbeta=1
   j=0
   iskip=0

   !--set up next panel
  311 continue
   x(2)=x(1)
   do il=1,nl
      y(il,2)=y(il,1)
   enddo
  313 continue
   if (jbeta.eq.0) jbeta=1
   if (jbeta.le.0) then
      if (lat.eq.1) then
         ep=enow-beta(-jbeta)*tevz
      else
         ep=enow-beta(-jbeta)*tev
      endif
      if (ep.eq.enow) then
         ep=sigfig(enow,8,-1)
      else
         ep=sigfig(ep,8,0)
      endif
   else
      if (lat.eq.1) then
         ep=enow+beta(jbeta)*tevz
      else
         ep=enow+beta(jbeta)*tev
      endif
      if (ep.eq.enow) then
         ep=sigfig(enow,8,+1)
         iskip=1
      else
         ep=sigfig(ep,8,0)
      endif
   endif
   if (ep.gt.x(2)) go to 316
   jbeta=jbeta+1
   go to 313
  316 continue
   ep=sigfig(ep,8,0)
   x(1)=ep
   call sigl(enow,ep,nnl,tev,nalpha,alpha,nbeta,beta,&
     sab,yt,nlmax,tol)
   do il=1,nl
      y(il,1)=yt(il)
   enddo

   !--adaptive subdivision of panel
   i=2
   ! compare linear approximation to true function
  330 continue
   if (i.eq.imax) go to 360
   if (iskip.eq.1) then
      iskip=0
      go to 360
   endif
   if (half*(y(1,i-1)+y(1,i))*(x(i-1)-x(i)).lt.tolmin) go to 360
   xm=half*(x(i-1)+x(i))
   xm=sigfig(xm,8,0)
   if (xm.le.x(i).or.xm.ge.x(i-1)) go to 360
   call sigl(enow,xm,nnl,tev,nalpha,alpha,nbeta,beta,&
     sab,yt,nlmax,tol)
   uu=0
   uum=0
   do 350 k=1,nl
   call terp1(x(i),y(k,i),x(i-1),y(k,i-1),xm,ym,2)
   if (k.gt.1) uu=uu+yt(k)
   if (k.gt.1) uum=uum+ym
   test=tol*abs(yt(k))
   test2=test
   if (k.gt.1) test2=tol
   if (abs(yt(k)-ym).gt.test2) go to 410
  350 continue
   test=2*tol*abs(uu)+uumin
   if (abs(uu-uum).gt.test) go to 410
   ! point passes.  save top point in stack and continue.
  360 continue
   j=j+1
   if (j.ge.jmax) call error('calcem','storage exceeded.',' ')
   if (j.gt.1) xsi(ie)=xsi(ie)+(x(i)-xlast)*(y(1,i)+ylast)*half
   if (j.gt.1) then
      uu=0
      u2=0
      u3=0
      nll=3
      do il=2,nl
         call legndr(y(il,i),p,nll)
         uu=uu+p(2)
         u2=u2+p(3)
         u3=u3+p(4)
      enddo
      uu=uu/(nl-1)
      uu=uu*y(1,i)
      u2=u2/(nl-1)
      u3=u3/(nl-1)
      u2=u2*y(1,i)
      u3=u3*y(1,i)
      ubar(ie)=ubar(ie)+half*(x(i)-xlast)*(uu+ulast)
      p2(ie)=p2(ie)+half*(x(i)-xlast)*(u2+u2last)
      p3(ie)=p3(ie)+half*(x(i)-xlast)*(u3+u3last)
   endif
   if (j.ne.3.or.xsi(ie).ge.tolmin) go to 380
   j=2
  380 continue
   jscr=7+(j-1)*(nl+1)
   scr(jscr)=x(i)
   if (y(1,i).ge.em9) then
      scr(1+jscr)=sigfig(y(1,i),9,0)
   else
      scr(1+jscr)=sigfig(y(1,i),8,0)
   endif
   do il=2,nl
      scr(il+jscr)=sigfig(y(il,i),9,0)
      if (scr(il+jscr).gt.unity) then
         !--only warn for big miss, but always fix the overflow
         !  use this same unity+0.0005 value in aceth
         if (scr(il+jscr).gt.unity+0.0005_kr) then
            write(strng,'("1cos=",f7.4,", set to 1.&
                          &  enow,e''=",2(1pe12.5))')&
                            scr(il+jscr),enow,scr(jscr)
            call mess('calcem',strng,'')
         endif
         scr(il+jscr)=unity
      endif
      if (scr(il+jscr).lt.-unity) then
         !--only warn for big miss, but always fix the underflow
         if (scr(il+jscr).lt.-(unity+0.0005_kr)) then
            write(strng,'("1cos=",f7.4,", set to -1.&
                          &  enow,e''=",2(1pe12.5),i3)')&
                            scr(il+jscr),enow,scr(jscr)
            call mess('calcem',strng,'')
         endif
         scr(il+jscr)=-unity
      endif
   enddo
   xlast=x(i)
   ylast=y(1,i)
   if (ylast.ne.zero) jnz=j
   ulast=0
   u2last=0
   u3last=0
   nll=3
   do il=2,nl
      call legndr(y(il,i),p,nll)
      ulast=ulast+p(2)
      u2last=u2last+p(3)
      u3last=u3last+p(4)
   enddo
   ulast=ulast*y(1,i)/(nl-1)
   u2last=u2last*y(1,i)/(nl-1)
   u3last=u3last*y(1,i)/(nl-1)
   i=i-1
   if (i.ge.2) go to 330
   jbeta=jbeta+1
   if (jbeta.le.nbeta) go to 311
   do il=1,nl
      y(il,i)=0
   enddo
   go to 430
   ! test fails.  add point to stack and continue.
  410 continue
   i=i+1
   x(i)=x(i-1)
   x(i-1)=xm
   do il=1,nl
      y(il,i)=y(il,i-1)
      y(il,i-1)=yt(il)
   enddo
   go to 330
   ! linearization complete.  write out result.
  430 continue
   j=j+1
   xsi(ie)=xsi(ie)+(x(i)-xlast)*(y(1,i)+ylast)*half
   uu=0
   u2=0
   u3=0
   ubar(ie)=ubar(ie)+half*(x(i)-xlast)*(uu+ulast)
   p2(ie)=p2(ie)+half*(x(i)-xlast)*(u2+u2last)
   p3(ie)=p3(ie)+half*(x(i)-xlast)*(u3+u3last)
   xsi(ie)=sigfig(xsi(ie),9,0)
   scr(7+(nl+1)*(j-1))=x(i)
   jscr=7+(j-1)*(nl+1)
   if (y(1,i).ge.em9) then
      scr(1+jscr)=sigfig(y(1,i),9,0)
   else
      scr(1+jscr)=sigfig(y(1,i),8,0)
   endif
   do il=2,nl
      scr(il+jscr)=sigfig(y(il,i),9,0)
      if (scr(il+jscr).gt.unity) then
         !--only warn for big miss, but always fix the overflow
         if (scr(il+jscr).gt.unity+0.0005_kr) then
            write(strng,'("2cos",f7.4,", set to 1.&
                          &  enow,e''=",2(1pe12.5))')&
                            scr(il+jscr),enow,scr(jscr)
            call mess('calcem',strng,'')
         endif
         scr(il+jscr)=unity
      endif
      if (scr(il+jscr).lt.-unity) then
         !--only warn for big miss, but always fix the underflow
         if (scr(il+jscr).lt.-(unity+0.0005_kr)) then
            write(strng,'("2cos",f7.4,", set to -1.&
                          &  enow,e''=",2(1pe12.5),i3)')&
                            scr(il+jscr),enow,scr(jscr)
            call mess('calcem',strng,'')
         endif
         scr(il+jscr)=-unity
      endif
   enddo
   if (y(1,1).ne.zero) jnz=j
   if (jnz.lt.j) j=jnz+1
   if (iprint.eq.2) then
      ubar(ie)=ubar(ie)/xsi(ie)
      p2(ie)=p2(ie)/xsi(ie)
      p3(ie)=p3(ie)/xsi(ie)
      write(nsyso,'(/,1x,"incident energy =",1pe13.6,&
                   &  5x,"cross section =",1pe13.6,&
                   &  5x,"mubar,p2,p3 =",3(1pe12.4))')&
                          enow,xsi(ie),ubar(ie),p2(ie),p3(ie)
      write(nsyso,'(/,5x,"exit energy",11x,"pdf",7x,"cosines")')
      write(nsyso,'(  3x,"---------------",5x,"-----------",2x,88("-"))')
      ll=6
      do jj=1,j
         write(nsyso,'(2x,1pe15.8,5x,1pe12.5,0p,8f11.6)')&
           (scr(ll+ii),ii=1,nlp)
         if (nl1.gt.nlp) write(nsyso,'(34x,8f11.6)')&
           (scr(ll+ii),ii=nlp1,nl1)
         ll=ll+nl1
      enddo
   endif
   scr(1)=0
   scr(2)=enow
   scr(3)=0
   scr(4)=0
   scr(5)=(nl+1)*j
   scr(6)=nl+1
   ncds=ncds+1+(j*(nl+1)+5)/6
   call listio(0,0,nscr,scr,nb,nw)
   loc=1
   do while (nb.ne.0)
      loc=loc+nw
      call moreio(0,0,nscr,scr(loc),nb,nw)
   enddo
   if (ie.lt.nne) go to 310
   go to 610

   !--compute kernel and write in mf6/law7 angle-energy format
 510 continue
   ltt=4
   math=matdp
   mfh=6
   mth=mtref
   tev=t*bk
   teff=teff*bk
   teff2=teff2*bk
   scr(1)=za
   scr(2)=awr
   scr(3)=0
   scr(4)=ltt ! temporary flag for this format
   scr(5)=1
   scr(6)=0
   call contio(0,0,nscr,scr,nb,nw)
   ncds=ncds+1
   scr(1)=1
   scr(2)=1
   scr(3)=0
   scr(4)=7
   scr(5)=1
   scr(6)=2
   scr(7)=2
   scr(8)=2
   scr(9)=1.e-5_kr
   scr(10)=1
   scr(11)=emax
   scr(12)=1
   nw=12
   call tab1io(0,0,nscr,scr,nb,nw)
   ncds=ncds+2
   scr(1)=0
   scr(2)=0
   scr(3)=0
   scr(4)=0
   scr(5)=1
   scr(6)=nne
   scr(7)=nne
   scr(8)=2
   nw=8
   call tab2io(0,0,nscr,scr,nb,nw)
   ncds=ncds+2
   cliq=0
   if (iinc.eq.1) go to 515
   if (sab(1,1).gt.sab(2,1))&
     cliq=(sab(1,1)-sab(1,2))*alpha(1)/beta(2)**2

   !--loop over given incident energy grid.
   !--the first pass computes the cross section and mu grid
  515 continue
   ie=0
  520 continue
   ie=ie+1
   enow=egrid(ie)
   if (ie.gt.1.and.temp.gt.break) enow=enow*temp/break
   enow=sigfig(enow,8,0)
   esi(ie)=enow
   j=0
   sum=0

   !--adaptive reconstruction of angular cross section
   u=-1
   x(2)=u
   call sigu(enow,u,tev,nalpha,alpha,nbeta,beta,sab,yu,nemax,tol)
   yy(2)=yu(1)
   xl=x(2)
   yl=yy(2)
   u=1
   x(1)=u
   call sigu(enow,u,tev,nalpha,alpha,nbeta,beta,sab,yu,nemax,tol)
   yy(1)=yu(1)
   i=2

   !--adaptive reconstruction
 530 continue
   if (i.eq.imax) go to 560
   xm=half*(x(i-1)+x(i))
   xm=sigfig(xm,7,0)
   if (xm.le.x(i).or.xm.ge.x(i-1)) go to 560
   call sigu(enow,xm,tev,nalpha,alpha,nbeta,beta,sab,yu,nemax,tol)
   if (x(i-1)-x(i).gt..25) go to 575
   ym=yy(i)+(xm-x(i))*(yy(i-1)-yy(i))/(x(i-1)-x(i))
   if (abs(yu(1)-ym).gt.2*tol*ym+tolmin) go to 575
   ! point passes.  save top point in stack and continue.
  560 continue
   j=j+1
   if (j.gt.mumax-1) call error('calcem','too many angles','see mumax')
   uj(j)=x(i)
   sj(j)=yy(i)
   if (j.gt.1) then
      sum=sum+half*(yy(i)+yl)*(x(i)-xl)
      xl=x(i)
      yl=yy(i)
   endif
   i=i-1
   if (i.ge.2) go to 530
   go to 580
   ! test fails.  add point to stack and continue.
  575 continue
   i=i+1
   x(i)=x(i-1)
   x(i-1)=xm
   yy(i)=yy(i-1)
   yy(i-1)=yu(1)
   go to 530
   ! linearization complete.  write out result.
  580 continue
   j=j+1
   uj(j)=x(1)
   sj(j)=yy(1)
   nmu=j
   ubar(ie)=0
   sum=sum+half*(yy(1)+yl)*(x(1)-xl)
   xsi(ie)=sum/2
   do i=2,nmu
      ubar(ie)=ubar(ie)+half*(uj(i)-uj(i-1))*(sj(i)+sj(i-1))*(uj(i)+uj(i-1))
   enddo
   ubar(ie)=half*ubar(ie)/sum
   if (iprint.eq.2) then
      write(nsyso,'(/i5,'' enow '',1p,e13.6,''   xsec '',e13.6,&
        &''   mubar  '',e13.6)') ie,enow,xsi(ie),ubar(ie)
      write(nsyso,'(''      num of mu '', i5)') nmu
      write(nsyso,'(/''            mu            theta      dsigma/dmu'')')
      do i=1,nmu
         write(nsyso,'(i5,1x,f15.8,1x,f12.4,1x,1p,e14.7)') i,uj(i),&
           acos(uj(i))*180.0/pi,sj(i)/2.0
      enddo
   endif

   !--now loop through the mu grid to write out the distributions
   mth=mtref
   scr(1)=0
   scr(2)=enow
   scr(3)=0
   scr(4)=0
   scr(5)=1
   scr(6)=nmu
   scr(7)=nmu
   scr(8)=2
   nw=8
   call tab2io(0,0,nscr,scr,nb,nw)
   ncds=ncds+1
   do il=1,nmu
      u=uj(il)
      call sigu(enow,u,tev,nalpha,alpha,nbeta,beta,sab,yu,nemax,tol)
      nep=nint(yu(2))
      j=0
      do i=1,nep
        j=nep-i
        if (yu(2*(nep-i)+4)/sum.gt.yumin) exit
      enddo
      nep=j
     if (iprint.eq.2) then
         write(nsyso,'(/'' mu = '',f15.8)') u
         write(nsyso,'('' (e-prime, pdf);  num of e-prime '', i5)') nep
         write(nsyso,*)
         !--test yu()/sum below: is this pdf normalized to 1.0 ?
         write(nsyso,'(1p,3(1x,e14.7,1x,e14.7,1x))')&
               (yu(2*i+1),yu(2*i+2)/sum,i=1,nep)
      endif
      scr(1)=0
      scr(2)=u
      scr(3)=0
      scr(4)=0
      scr(5)=1
      scr(6)=nep
      scr(7)=nep
      scr(8)=2
      k=8
      istart=1
     595 continue
      iend=nep
      if ((iend-istart).ge.npage/2) iend=istart+npage/2-1
      j=k-1
      ib=istart-1
     596 continue
      j=j+2
      ib=ib+1
      scr(j)=yu(1+2*ib)
      scr(j+1)=yu(2+2*ib)*2/sum
      if (ib.lt.iend) go to 596
      nw=j+1
      if (k.eq.0) go to 597
      k=0
      call tab1io(0,0,nscr,scr,nb,nw)
      if (nb.eq.0) go to 598
      istart=iend+1
      go to 595
    597 continue
      call moreio(0,0,nscr,scr,nb,nw)
      if (nb.eq.0) go to 598
      istart=iend+1
      go to 595
     598 continue
      ncds=ncds+1+(j*(nep+1)+5)/6
   enddo

   if (ie.lt.nne) go to 520

   !--calculate incoherent inelastic scattering cross sections
 610 continue
   do ie=1,20
      ex(ie)=0
   enddo
   nexn=nex
   if (nex.eq.2) nexn=3
   ie=0
   ien=1
   do while (ien.gt.0)
      ie=ie+1
      call finda(ie,ex,nex,iold,bufo,nbuf)
      if (iinc.eq.1) then
         ! use elastic cross section  for free atom scattering.
         ex(3)=ex(2)
      else
         ! use computed cross section for bound atoms.
         enow=ex(1)
         xs=terp(esi,xsi,nne,enow,nlt)
         if (ie.eq.ne) xs=0
         ex(3)=xs
      endif
      ien=ie
      if (ie.eq.ne) ien=-ien
      call loada(ien,ex,nexn,inew,bufn,nbuf)
   enddo
   esi(nne+1)=0
   isave=iold
   iold=inew
   inew=isave
   nex=nexn
   if (ie.eq.ne) xs=0

   !--calcem is finished.
   call asend(0,nscr)
   deallocate(alpha)
   deallocate(beta)
   deallocate(sab)
   return
   end subroutine calcem

   real(kr) function sig(e,ep,u,tev,nalpha,alpha,nbeta,beta,sab)
   !-------------------------------------------------------------------
   ! Compute the differential scattering cross section from e to
   ! ep through the angle with the cosine u from endf tabulated
   ! data or an analytic law.
   !-------------------------------------------------------------------
   use physics ! provides pi
   use util    ! provides error
   ! externals
   integer::nalpha,nbeta
   real(kr)::e,ep,u,tev,bbm
   real(kr)::alpha(nalpha),beta(nbeta),sab(nalpha,nbeta)
   ! internals
   integer::nb1,na1,i,ib,ia
   real(kr)::rtev,bb,a,sigc,b,c,bbb,s,s1,s2,s3,arg
   real(kr)::tfff,tfff2,rat
   real(kr),parameter::sigmin=1.e-10_kr
   real(kr),parameter::sabflg=-225.e0_kr
   real(kr),parameter::amin=1.e-6_kr
   real(kr),parameter::test1=0.2e0_kr
   real(kr),parameter::test2=30.e0_kr
   real(kr),parameter::zero=0

   !--common factors.
   rtev=1/tev
   bb=(ep-e)*rtev
   a=(e+ep-2*u*sqrt(e*ep))/(az*tev)
   if (a.lt.amin) a=amin
   sigc=sqrt(ep/e)*rtev/2
   b=abs(bb)
   c=sqrt(4*pi)

   !--tabulated s(alpha,beta).
   if (iinc.ne.2) go to 200
   if (lat.eq.1) b=b*tev/tevz
   if (lat.eq.1) a=a*tev/tevz
   if (a.gt.alpha(nalpha)) go to 170
   if (lasym.eq.1) then
      bbm=bb
      if (lat.eq.1) bbm=bb*tev/tevz
      if (bbm.gt.beta(nbeta)) go to 170
      if (bbm.lt.beta(1)) go to 170
   else
      if (b.gt.beta(nbeta)) go to 170
   endif
   nb1=nbeta-1
   na1=nalpha-1
   bbb=b
   if (lasym.eq.1.and.bb.lt.zero) bbb=-b
   do i=1,nb1
      ib=i
      if (bbb.lt.beta(i+1)) exit
   enddo
   do i=1,na1
      ia=i
      if (a.lt.alpha(i+1)) exit
   enddo
   if (cliq.eq.zero.or.a.ge.alpha(1)) go to 150
   if (lasym.eq.1) go to 150
   if (b.gt.test1) go to 150
   s=sab(1,1)+log(alpha(1)/a)/2-cliq*b**2/a
   if (s.lt.sabflg) s=sabflg
   go to 160
  150 continue
   if (a*az.lt.test2.and.b.lt.test2) go to 155
   if (sab(ia,ib).le.sabflg) go to 170
   if (sab(ia+1,ib).le.sabflg) go to 170
   if (sab(ia,ib+1).le.sabflg) go to 170
   if (sab(ia+1,ib+1).le.sabflg) go to 170
  155 continue
   if (ia+1.eq.nalpha) ia=ia-1
   if (ib+1.eq.nbeta) ib=ib-1
   call terpq(alpha(ia),sab(ia,ib),alpha(ia+1),sab(ia+1,ib),&
     alpha(ia+2),sab(ia+2,ib),a,s1)
   call terpq(alpha(ia),sab(ia,ib+1),alpha(ia+1),sab(ia+1,ib+1),&
     alpha(ia+2),sab(ia+2,ib+1),a,s2)
   call terpq(alpha(ia),sab(ia,ib+2),alpha(ia+1),sab(ia+1,ib+2),&
     alpha(ia+2),sab(ia+2,ib+2),a,s3)
   call terpq(beta(ib),s1,beta(ib+1),s2,beta(ib+2),s3,bbb,s)
  160 continue
   sig=0
   if (s-bb/2.gt.sabflg) sig=exp(s-bb/2)
   sig=sigc*sb*sig
   if (sig.lt.sigmin) sig=0
   return

   !--short collision time for large beta or alpha
  170 continue
   if (lat.eq.1) b=b*tevz*rtev
   if (lat.eq.1) a=a*tevz*rtev
   tfff=teff
   tfff2=teff2
   !if (az.lt.3) then
   !   if (e.gt.10) then
   !      tfff=tev
   !      tfff2=tev
   !   else if (e.gt.2) then
   !      rat=(e-2)/8
   !      tfff= (1-rat)*teff+rat*tev
   !      tfff2= (1-rat)*teff2+rat*tev
   !   endif
   !endif
   s=0
   arg=(a-b)**2*tev/(4*a*tfff)+(b+bb)/2
   if (-arg.gt.sabflg) s=exp(-arg)/(c*sqrt(a*tfff*rtev))
   sig=sigc*sb*s
   if (sb2.gt.zero) then
      a=a*az/az2
      arg=(a-b)**2*tev/(4*a*tfff2)+(b+bb)/2
      s2=0
      if (-arg.gt.sabflg) s2=exp(-arg)/(c*sqrt(a*tfff2*rtev))
      sig=sig+sigc*sb2*s2
   endif
    if (abs(e-10).lt..01.and.abs(u-.99219).lt..0001) then
    endif
   if (sig.lt.sigmin) sig=0
   return

   !--free-gas scattering.
  200 continue
   if (iinc.ne.1) go to 300
   s=0
   arg=(a+bb)**2/(4*a)
   if (-arg.gt.sabflg) s=exp(-arg)/(c*sqrt(a))
   sig=sigc*sb*s
   if (sig.lt.sigmin) sig=0
   return

   !--other options not yet implemented.
  300 continue
   call error('sig','illegal option.',' ')
   sig=0
   return
   end function sig

   subroutine terpq(x1,y1,x2,y2,x3,y3,x,y)
   !-------------------------------------------------------------------
   ! Compute y(x) by quadratic interpolation,
   ! except use log-lin if x.lt.x1 and lin-lin if x.gt.x3.
   ! and use lin-lin if the function takes big steps (corners).
   !-------------------------------------------------------------------
   use endf ! provides terp1
   ! externals
   real(kr)::x1,y1,x2,y2,x3,y3,x,y
   ! internals
   real(kr)::b,c
   real(kr),parameter::sabflg=-225.e0_kr
   real(kr),parameter::step=2.e0_kr

   if (x.lt.x1) then
      if (y1.gt.y2) then
         y=y1
      else
         call terp1(x1,y1,x2,y2,x,y,3)
      endif
   else if (x.gt.x3) then
      if (y3.gt.y2) then
         y=y3
      else
         call terp1(x2,y2,x3,y3,x,y,2)
      endif
   else if (abs(y1-y2).gt.step.or.abs(y2-y3).gt.step) then
      if (x.lt.x2) then
         call terp1(x1,y1,x2,y2,x,y,2)
      else
         call terp1(x2,y2,x3,y3,x,y,2)
      endif
   else
      b=(y2-y1)*(x3-x1)/((x2-x1)*(x3-x2))
      b=b-(y3-y1)*(x2-x1)/((x3-x1)*(x3-x2))
      c=(y3-y1)/((x3-x1)*(x3-x2))
      c=c-(y2-y1)/((x2-x1)*(x3-x2))
      y=y1+b*(x-x1)+c*(x-x1)*(x-x1)
   endif
   if (y.lt.sabflg) y=sabflg
   return
   end subroutine terpq

   subroutine sigl(e,ep,nlin,tev,nalpha,alpha,nbeta,beta,sab,&
     s,nlmax,tolin)
   !-------------------------------------------------------------------
   ! Compute the cross section and legendre components or equally-
   ! probable angles for the scattering from e to ep.  Uses linear
   ! reconstruction of the angular distribution computed by sig.
   !-------------------------------------------------------------------
   use util  ! provides error,mess
   use mathm ! provide legndr
   ! externals
   integer::nlin,nalpha,nbeta,nlmax
   real(kr)::e,ep,tev,alpha(nalpha),beta(nbeta),sab(nalpha,nbeta)
   real(kr)::s(nlmax),tolin
   ! internals
   integer::nl,i,j,il
   real(kr)::b,seep,sum,xl,yl,ymax,xm,ym,test,test2
   real(kr)::rnbin,fract,gral,add,xil,xn,f,rf,disc,yn,xbar,rfract
   real(kr)::yt,tol,s1bb
   integer,parameter::imax=20
   real(kr)::x(imax),y(imax)
   real(kr)::p(nlin)
   character(60)::strng
   real(kr),parameter::zero=0
   real(kr),parameter::xtol=.00001e0_kr
   real(kr),parameter::ytol=.001e0_kr
   real(kr)::one=1
   real(kr),parameter::half=.5e0_kr
   real(kr),parameter::third=.333333333e0_kr
   real(kr),parameter::sigmin=1.e-32_kr
   real(kr),parameter::eps=1.e-3_kr
   real(kr),parameter::shade=.99999999e0_kr

   !--constant factors
   b=(ep-e)/tev
   tol=half*tolin
   nl=nlin
   if (nl.lt.0) nl=-nl
   b=abs(b)
   if (lat.eq.1.and.iinc.eq.2) b=b*tev/tevz
   s1bb=sqrt(1+b*b)
   if (ep.ne.zero) seep=1/sqrt(e*ep)

   !--adaptive calculation of cross section
   i=3
   sum=0
   x(3)=-1
   xl=x(3)
   y(3)=sig(e,ep,x(3),tev,nalpha,alpha,nbeta,beta,sab)
   yl=y(3)
   if (ep.eq.zero) x(2)=0
   if (ep.ne.zero) x(2)=half*(e+ep-(s1bb-1)*az*tev)*seep
   if (abs(x(2)).gt.1-eps) x(2)=0.99e0_kr
   x(2)=sigfig(x(2),8,0)
   y(2)=sig(e,ep,x(2),tev,nalpha,alpha,nbeta,beta,sab)
   x(1)=+1
   y(1)=sig(e,ep,x(1),tev,nalpha,alpha,nbeta,beta,sab)
   ymax=y(2)
   if (y(1).gt.ymax) ymax=y(1)
   if (y(3).gt.ymax) ymax=y(3)
   if (ymax.lt.eps) ymax=eps
  110 continue
   if (i.eq.imax) go to 120
   xm=half*(x(i-1)+x(i))
   xm=sigfig(xm,8,0)
   ym=half*(y(i-1)+y(i))
   yt=sig(e,ep,xm,tev,nalpha,alpha,nbeta,beta,sab)
   test=tol*abs(yt)+tol*ymax/50
   test2=ym+ymax/100
   if (abs(yt-ym).le.test.and.abs(y(i-1)-y(i)).le.test2.and.&
     (x(i-1)-x(i)).lt.half) go to 120
   if (x(i-1)-x(i).lt.xtol) go to 120
   i=i+1
   x(i)=x(i-1)
   y(i)=y(i-1)
   x(i-1)=xm
   y(i-1)=yt
   go to 110
  120 continue
   sum=sum+half*(y(i)+yl)*(x(i)-xl)
   xl=x(i)
   yl=y(i)
   i=i-1
   if (i.gt.1) go to 110
   if (i.eq.1) go to 120
   s(1)=sum
   if (sum.gt.sigmin) go to 130
   do il=1,nl
      s(il)=0
   enddo
   return

   !--prime stack for equally-probable angles
  130 continue
   nbin=nl-1
   rnbin=one/nbin
   fract=sum*rnbin
   rfract=1/fract
   sum=0
   gral=0
   do il=2,nl
      s(il)=0
   enddo
   j=0

   !--adaptive linearization
   i=3
   x(3)=-1
   xl=x(3)
   y(3)=sig(e,ep,x(3),tev,nalpha,alpha,nbeta,beta,sab)
   if (ep.eq.zero) x(2)=0
   if (ep.ne.zero) x(2)=half*(e+ep-(s1bb-1)*az*tev)*seep
   if (abs(x(2)).gt.1-eps) x(2)=0.99e0_kr
   x(2)=sigfig(x(2),8,0)
   y(2)=sig(e,ep,x(2),tev,nalpha,alpha,nbeta,beta,sab)
   x(1)=+1
   y(1)=sig(e,ep,x(1),tev,nalpha,alpha,nbeta,beta,sab)
   ymax=y(1)
   if (y(2).gt.ymax) ymax=y(2)
   if (y(3).gt.ymax) ymax=y(3)
   if (ymax.lt.eps) ymax=eps
  150 continue
   if (i.eq.imax) go to 160
   xm=half*(x(i-1)+x(i))
   xm=sigfig(xm,8,0)
   ym=half*(y(i-1)+y(i))
   yt=sig(e,ep,xm,tev,nalpha,alpha,nbeta,beta,sab)
   test=tol*abs(yt)+tol*ymax/50
   test2=ym+ymax/100
   if (abs(yt-ym).le.test.and.abs(y(i-1)-y(i)).le.test2.and.&
     (x(i-1)-x(i)).lt.half) go to 160
   if (x(i-1)-x(i).lt.xtol) go to 160
   i=i+1
   x(i)=x(i-1)
   y(i)=y(i-1)
   x(i-1)=xm
   y(i-1)=yt
   go to 150

   !--check bins for this panel
  160 continue
   add=half*(y(i)+yl)*(x(i)-xl)
   if (x(i).eq.xl) go to 250
   xil=1/(x(i)-xl)
   if (i.eq.1.and.j.eq.nbin-1) go to 165
   if (sum+add.ge.fract*shade.and.j.lt.nbin-1) go to 170
   sum=sum+add
   gral=gral+half*(yl*x(i)-y(i)*xl)*(x(i)+xl)&
     +third*(y(i)-yl)*(x(i)**2+x(i)*xl+xl**2)
   go to 250
  165 continue
   xn=x(i)
   j=j+1
   go to 190
  170 continue
   j=j+1
   if (yl.lt.sigmin) go to 175
   test=(fract-sum)*(y(i)-yl)/((x(i)-xl)*yl**2)
   if (abs(test).gt.ytol) go to 175
   xn=xl+(fract-sum)/yl
   if (xn.gt.x(i)) go to 180
   if (xn.ge.xl.and.xn.le.x(i)) go to 190
   call error('sigl','no legal solution.',' ')
  175 continue
   f=(y(i)-yl)*xil
   rf=1/f
   disc=(yl*rf)**2+2*(fract-sum)*rf
   if (disc.lt.zero) then
      write(strng,'(''disc='',1p,e12.4)') disc
      call mess('sigl',strng,'set to abs value and continue')
      disc=abs(disc)
   endif
   if (f.gt.zero) xn=xl-(yl*rf)+sqrt(disc)
   if (f.lt.zero) xn=xl-(yl*rf)-sqrt(disc)
   if (xn.gt.xl.and.xn.le.x(i)) go to 190
   if (xn.gt.xl.and.xn.lt.(x(i)+ytol*(x(i)-xl))) go to 180
   call error('sigl','no legal solution (quadratic path).',' ')
  180 continue
   xn=x(i)
  190 continue
   yn=yl+(y(i)-yl)*(xn-xl)*xil
   gral=gral+(xn-xl)*(yl*half*(xn+xl)&
     +(y(i)-yl)*xil*(-xl*half*(xn+xl)&
     +third*(xn**2+xn*xl+xl**2)))
   xbar=gral*rfract

   !--compute legendre components
   if (nlin.ge.0) then
      call legndr(xbar,p,nl)
      do il=2,nl
         s(il)=s(il)+p(il)*rnbin
      enddo

   !--output equally probable angles
   else
      s(j+1)=xbar
   endif

   !--continue bin loop and linearization loop
   xl=xn
   yl=yn
   sum=0
   gral=0
   if (j.eq.nbin) go to 260
   if (xl.lt.x(i)) go to 160
  250 continue
   xl=x(i)
   yl=y(i)
   i=i-1
   if (i.gt.1) go to 150
   if (i.eq.1) go to 160
  260 continue
   return
   end subroutine sigl

   subroutine sigu(e,u,tev,nalpha,alpha,nbeta,beta,sab,&
     s,nemax,tolin)
   !-------------------------------------------------------------------
   ! Compute the secondary energy distribution scattering for cosine u.
   ! Uses linear reconstruction with the cross section from function sig.
   !-------------------------------------------------------------------
   use util ! provides sigfig
   ! externals
   integer::nalpha,nbeta,nemax
   real(kr)::e,u,tev,alpha(nalpha),beta(nbeta),sab(nalpha,nbeta)
   real(kr)::s(2*nemax),tolin
   ! internals
   integer::i,j,jbeta
   real(kr)::sum,xl,yl,xm,ym,test
   real(kr)::yt,tol
   real(kr)::root1,root2
   integer,parameter::imax=20
   real(kr)::x(imax),y(imax)
   real(kr),parameter::zero=0
   real(kr),parameter::half=.5e0_kr
   real(kr),parameter::tolmin=1.e-6_kr
   real(kr),parameter::bmax=20

   !--constant factors
   tol=tolin
   do i=1,2*nemax
      s(i)=0
   enddo

   root1=(u*sqrt(e)+sqrt(u*u*e+(az-1)*(az+1)*e))/(az+1)
   root2=(u*sqrt(e)-sqrt(u*u*e+(az-1)*(az+1)*e))/(az+1)

   !--adaptive calculation of cross section
   sum=0
   x(1)=0
   y(1)=sig(e,x(1),u,tev,nalpha,alpha,nbeta,beta,sab)
   jbeta=-nbeta
   if (lasym.gt.0) jbeta=1
   j=0
   xl=0
   yl=0

   !--set up next panel
  111 continue
   x(2)=x(1)
   y(2)=y(1)
  113 continue
   if (jbeta.eq.0) jbeta=1
   if (jbeta.le.0) then
      if (lat.eq.1) then
         x(1)=e-beta(-jbeta)*tevz
      else
         x(1)=e-beta(-jbeta)*tev
      endif
      x(1)=sigfig(x(1),8,0)
      if (x(1).eq.e) x(1)=sigfig(e,8,-1)
   else
      if (lat.eq.1) then
         x(1)=e+beta(jbeta)*tevz
      else
         x(1)=e+beta(jbeta)*tev
      endif
   endif
   if (x(1).gt.x(2)) go to 116
   jbeta=jbeta+1
   go to 113
  116 continue
   if (u.lt.zero.and.root1**2.gt.1.01*x(2).and.root1**2.lt.x(1)) then
      x(1)=root1**2
   endif
   x(1)=sigfig(x(1),8,0)
   y(1)=sig(e,x(1),u,tev,nalpha,alpha,nbeta,beta,sab)
   i=2

   !--compare linear approximation to true function
  150 continue
   if (i.eq.imax) go to 160
   if (i.gt.3.and.half*(y(i-1)+y(i))*(x(i-1)-x(i)).lt.tolmin) go to 160
   xm=half*(x(i-1)+x(i))
   xm=sigfig(xm,8,0)
   if (xm.le.x(i).or.xm.ge.x(i-1)) go to 160
   ym=half*(y(i-1)+y(i))
   yt=sig(e,xm,u,tev,nalpha,alpha,nbeta,beta,sab)
    if (abs(u-.99219).lt..0001) then
       if (abs(e-10).lt..01) then
       endif
    endif
   test=tol*abs(yt)
   if (abs(yt-ym).le.test) go to 160

   !--point fails
   i=i+1
   x(i)=x(i-1)
   y(i)=y(i-1)
   x(i-1)=xm
   y(i-1)=yt
   go to 150

   !--point passes
  160 continue
   j=j+1
   s(2*j+1)=x(i)
   s(2*j+2)=y(i)
   if (j.gt.1) sum=sum+(y(i)+yl)*(x(i)-xl)
   xl=x(i)
   yl=y(i)
   if (j.ge.nemax-1) go to 170
   if (jbeta.gt.0) then
     if (beta(jbeta).gt.bmax) go to 170
   endif

   !--continue bin loop and linearization loop
   i=i-1
   if (i.gt.1) go to 150
   jbeta=jbeta+1
   if (jbeta.le.nbeta) go to 111
   if (i.eq.1) go to 160
  170 continue
   s(1)=sum
   s(2)=j

   return
   end subroutine sigu

   subroutine tpend(iold,itemp,ne,nex,icoh,icopy,mtref)
   !-------------------------------------------------------------------
   ! Write the output pendf tape.
   !-------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso,nsyse
   use endf   ! provides endf routines and variables
   use util   ! provides repoz,error,finda,timer,sigfig
   ! externals
   integer::iold,itemp,ne,nex,icoh,icopy,mtref
   ! internals
   integer::i,matd,jcopy,n6,iverp,mf2,mfd,mtd
   integer::lim,ix,j,mti,mfi,mf,mt,mfn,mtn,mat
   integer::istart,iend,k,ib,ifend,ie,ltt
   integer::nb,nw,indx,nx,nc,nxn,nne,nl1,nmu,imu
   real(kr)::temp,rxsec,yy,sec
   real(kr)::ex(20)
   real(kr),dimension(:),allocatable::dico,dicn,sav
   real(kr),parameter::em9=1.e-9_kr
   real(kr),parameter::etop=20.e6_kr
   real(kr),parameter::zero=0
   integer,parameter::nwscr=500000

   !--initialize.
   temp=tempr(itemp)
   do i=1,20
      ex(i)=0
   enddo
   matd=matdp
   call repoz(nscr)
   call repoz(nscr2)
   jcopy=0

   !--copy skipped temperatures, if any.
   if (icopy.gt.0) then
      do while (jcopy.lt.icopy)
         call contio(nscr2,nout,0,scr,nb,nw)
         call tomend(nscr2,nout,0,scr)
         jcopy=jcopy+1
      enddo
      icopy=0
   endif

   !--update directory
   call findf(matd,1,451,nin)
   indx=1
   call contio(nin,0,0,scr(indx),nb,nw)
   n6=n2h
   call contio(nin,0,0,scr(indx+6),nb,nw)
   if (n1h.ne.0) then
      iverp=4
      nx=n6
   else if (n2h.eq.0) then
      iverp=5
   else
      iverp=6
   endif
   call skiprz(nin,-1)
   indx=indx+nw
   if (iverp.ge.5) call contio(nin,0,0,scr(indx),nb,nw)
   if (iverp.ge.5) indx=indx+nw
   if (iverp.ge.6) call contio(nin,0,0,scr(indx),nb,nw)
   if (iverp.ge.6) indx=indx+nw
   call hdatio(nin,0,0,scr(indx),nb,nw)
   scr(indx)=temp
   if (iverp.ge.5) nx=nint(scr(indx+5))
   indx=indx+nw
   do while (nb.ne.0)
      if (indx.gt.nwscr) call error('tpend','storage exceeded.',' ')
      call moreio(nin,0,0,scr(indx),nb,nw)
      indx=indx+nw
   enddo
   nw=6*nx
   allocate(dico(nw))
   nw=6*nx
   if (iinc.gt.0) nw=nw+12
   if (icoh.gt.0) nw=nw+12
   if (icoh.gt.20) nw=nw+12
   allocate(dicn(nw))
   allocate(sav(6))
   nc=nx
   call dictio(nin,0,0,dico,nb,nc)
   do i=1,nw
      dicn(i)=0
   enddo
   mf2=0
   mfd=3
   mtd=mtref
   lim=6*(nx-1)
   ix=3
   nc=3+(ne+2)/3
   i=0
   j=0
  172 continue
   mti=0
   if (i.gt.lim) go to 176
   mfi=nint(dico(i+3))
   if (mfi.eq.2) mf2=1
   mti=nint(dico(i+4))
  173 continue
   if (mfi.gt.mfd) go to 176
   if (mfi.eq.mfd.and.mti.ge.mtd) go to 176
  174 continue
   dicn(j+3)=dico(i+3)
   dicn(j+4)=dico(i+4)
   dicn(j+5)=dico(i+5)
   i=i+6
   j=j+6
   go to 172
  176 continue
   if (mfd.gt.6) go to 190
   dicn(j+3)=mfd
   dicn(j+4)=mtd
   dicn(j+5)=nc
   j=j+6
   if (mtd.eq.mti) then
      i=i+6
      mfi=nint(dico(i+3))
      mti=nint(dico(i+4))
   endif
   if (ix.eq.nex) go to 180
   ix=ix+1
   mtd=mtd+1
   if (mfd.eq.6) nc=ncdse
   go to 176
  180 continue
   if (mfd.eq.6) go to 184
   if (iinc.eq.0) go to 184
   nc=ncds
   ix=3
   mfd=6
   mtd=mtref
   if (i.gt.lim) go to 176
   go to 173
  184 continue
   mfd=1000
   if (mtd.eq.mti) i=i+6
   if (i.le.lim) go to 174
  190 continue
   nxn=j/6
   mf=1
   mt=451
   if (iverp.lt.5) scr(6)=nxn
   call contio(0,nout,0,scr,nb,nw)
   indx=1+nw
   if (iverp.ge.5) call contio(0,nout,0,scr(indx),nb,nw)
   if (iverp.ge.5) indx=indx+nw
   if (iverp.ge.6) call contio(0,nout,0,scr(indx),nb,nw)
   if (iverp.ge.6) indx=indx+nw
   if (iverp.ge.5) scr(indx+5)=nxn
   call hdatio(0,nout,0,scr(indx),nb,nw)
   indx=indx+nw
   do while (nb.ne.0)
      if (indx.gt.1+nwscr-1)&
        call error('tpend','storage exceeded.',' ')
      call moreio(0,nout,0,scr(indx),nb,nw)
      indx=indx+nw
   enddo
   dicn(5)=dico(5)+nxn-nx
   call dictio(0,nout,0,dicn,nb,nxn)
   call tofend(nin,nout,0,scr)
   if (mf2.eq.1) call tofend(nin,nout,0,scr)

   !--copy rest of this mat and temperature from nin to nout,
   !--adding thermal scattering and matrix to nout.
   mfn=3
   mtn=mtref
  200 continue
   call contio(nin,0,0,sav,nb,nw)
   mat=math
   mf=mfh
   mt=mth
   if (math.ne.matd) go to 210
   if (mfh.eq.0) go to 200
   if (mfh.gt.mfn) go to 210
   if (mfh.eq.mfn.and.mth.ge.mtn) go to 210
   math=mat
   mfh=mf
   mth=mt
   call contio(0,nout,0,sav,nb,nw)
   call tosend(nin,nout,0,scr)
   go to 200

   !--file 3.
  210 continue
   ix=2
  215 continue
   ix=ix+1
   scr(1)=za
   scr(2)=awr
   scr(3)=0
   scr(4)=0
   scr(5)=0
   scr(6)=0
   math=matd
   mfh=3
   mth=mtn
   call contio(0,nout,0,scr,nb,nw)
   scr(1)=temp
   scr(2)=0
   scr(5)=1
   scr(6)=ne+1
   scr(7)=ne+1
   scr(8)=2
   k=8
   istart=1
  220 continue
   iend=ne+1
   if ((iend-istart).ge.npage/2) iend=istart+npage/2-1
   j=k-1
   ib=istart-1
  222 continue
   j=j+2
   ib=ib+1
   if (ib.ge.ne+1) then
      ex(1)=etop
      ex(ix)=0
   else
      call finda(ib,ex,nex,iold,bufo,nbuf)
   endif
   scr(j)=sigfig(ex(1),7,0)
   if (ex(ix).ge.em9) then
      scr(j+1)=sigfig(ex(ix),7,0)
   else
      scr(j+1)=sigfig(ex(ix),6,0)
   endif
   if (ib.lt.iend) go to 222
   nw=j+1
   if (k.eq.0) go to 224
   k=0
   call tab1io(0,nout,0,scr,nb,nw)
   if (nb.eq.0) go to 230
   istart=iend+1
   go to 220
  224 continue
   call moreio(0,nout,0,scr,nb,nw)
   if (nb.eq.0) go to 230
   istart=iend+1
   go to 220
  230 continue
   call asend(nout,0)
   if (mf.ne.mfn.or.mt.ne.mtn) go to 235
   call tosend(nin,0,0,scr)
   call contio(nin,0,0,sav,nb,nw)
   mat=math
   mf=mfh
   mt=mth
  235 continue
   if (ix.eq.nex) go to 240
   mtn=mtn+1
   go to 215
  240 continue
   if (mf.eq.3) go to 245
   call afend(nout,0)
   go to 250
  245 continue
   math=mat
   mfh=mf
   mth=mt
   call contio(0,nout,0,sav,nb,nw)
   call tofend(nin,nout,0,scr)
   call contio(nin,0,0,sav,nb,nw)
   mat=math
   mf=mfh
   mt=mth

   !--file 6
  250 continue
   if (iinc.eq.0.and.icoh.eq.0) go to 295
   ifend=0
   mfn=6
   mtn=mtref
   ix=2
  252 continue
   if (mf.eq.0) go to 254
   if (mf.gt.mfn) go to 254
   if (mf.eq.mfn.and.mt.ge.mtn) go to 254
   math=mat
   mfh=mf
   mth=mt
   call contio(0,nout,0,sav,nb,nw)
   call tosend(nin,nout,0,scr)
   call contio(nin,0,0,sav,nb,nw)
   ifend=0
   if (mfh.eq.0.and.mth.eq.0) ifend=1
   mat=math
   mf=mfh
   mt=mth
   go to 252
  254 continue
   ix=ix+1
   call contio(nscr,0,0,scr,nb,nw)
   ltt=nint(scr(4))
   scr(4)=1
   scr(5)=1
   call contio(0,nout,0,scr,nb,nw)
   if (ltt.eq.6) then
      call tab1io(nscr,nout,0,scr,nb,nw)
      call tab2io(nscr,nout,0,scr,nb,nw)
      nne=n2h
      do ie=1,nne
         call listio(nscr,0,0,scr,nb,nw)
         scr(8)=1
         call listio(0,nout,0,scr,nb,nw)
      enddo
   else if (ltt.eq.5) then
      call tab1io(nscr,nout,0,scr,nb,nw)
      call tab2io(nscr,nout,0,scr,nb,nw)
      nne=n2h
      do ie=1,nne
         if (xsi(ie).eq.zero)&
           call error('tpend','cross section=0.',' ')
         rxsec=1/xsi(ie)
         indx=1
         call listio(nscr,0,0,scr(indx),nb,nw)
         indx=indx+nw
         do while (nb.ne.0)
            if (indx.gt.nwscr)&
              call error('tpend','storage exceeded.',' ')
            call moreio(nscr,0,0,scr(indx),nb,nw)
            indx=indx+nw
         enddo
         nl1=nint(scr(6))
         j=nint(scr(5))/nl1
         indx=8-nl1
         do i=1,j
            yy=scr(nl1*i+indx)*rxsec
            if (yy.ge.em9) then
               scr(nl1*i+indx)=sigfig(yy,7,0)
            else
               scr(nl1*i+indx)=sigfig(yy,6,0)
            endif
         enddo
         indx=1
         call listio(0,nout,0,scr(indx),nb,nw)
         do while (nb.ne.0)
            indx=indx+nw
            call moreio(0,nout,0,scr(indx),nb,nw)
         enddo
      enddo
   else if (ltt.eq.4) then
      call tab1io(nscr,nout,0,scr,nb,nw)
      call tab2io(nscr,nout,0,scr,nb,nw)
      nne=n2h
      do ie=1,nne
         call tab2io(nscr,nout,0,scr,nb,nw)
         nmu=n2h
         do imu=1,nmu
            indx=1
            call tab1io(nscr,nout,0,scr,nb,nw)
            do while (nb.ne.0)
               indx=indx+nw
               call moreio(nscr,nout,0,scr(indx),nb,nw)
            enddo
         enddo
      enddo
   endif
   call tosend(nscr,nout,0,scr)
   if (mf.ne.mfn.and.ix.eq.nex.and.ifend.eq.0) call afend(nout,0)
   if (mf.ne.mfn.or.mt.ne.mtn) go to 290
   call tosend(nin,0,0,scr)
   call contio(nin,0,0,sav,nb,nw)
   mat=math
   mf=mfh
   mt=mth
  290 continue
   if (ix.eq.nex) go to 295
   mtn=mtn+1
   go to 254

   !--finished with this temperature
  295 continue
   math=mat
   mfh=mf
   mth=mt
   call contio(0,nout,0,sav,nb,nw)
   if (mat.ne.0) call tomend(nin,nout,0,scr)
   deallocate(dicn)
   deallocate(dico)
   deallocate(sav)
   call timer(sec)
   write(nsyso,'(/&
     &'' wrote thermal data for temp ='',1p,e11.4,28x,0p,f8.1,''s'')')&
     temp,sec
   write(nsyse,'(/&
     &'' wrote thermal data for temp ='',1p,e11.4,28x,0p,f8.1,''s'')')&
     temp,sec
   return
   end subroutine tpend

end module thermm

module wimsm
   ! provides subroutine wimsr for NJOY2016
   use locale
   implicit none
   private
   public wimsr
   ! global wimsr variables
   integer::iprint
   integer::iverw
   integer::ngnd,nfg,nrg,igref
   integer::mat
   integer::nfid
   real(kr)::rdfid
   integer::iburn
   integer::ntemp
   integer,parameter::ntmax=10
   real(kr)::tempr(ntmax)
   integer::nsigz
   real(kr)::sgref
   integer::ires
   real(kr)::sigp
   integer::mti,mtc
   integer::ip1opt
   integer::inorf
   integer::isof
   integer::ifprod
   integer::ngendf,nout,nscr0,nscr1,nscr2,nscr3,nscr4
   real(kr)::awr
   integer::iznum
   integer,parameter::nymax=100
   real(kr)::yield(nymax)
   integer::ifisp(nymax)
   integer::ifiss,nfiss
   integer::isg
   integer::ixs
   real(kr),dimension(:),allocatable::snu
   real(kr),dimension(:),allocatable::spot
   real(kr),dimension(:),allocatable::abs2
   real(kr),dimension(:),allocatable::uff
   real(kr),dimension(:),allocatable::glam
   real(kr),dimension(:),allocatable::p1flx
   real(kr),dimension(:),allocatable::flux
   real(kr),dimension(:),allocatable::egb
   integer,parameter::nwscr=30000
   real(kr),dimension(:),allocatable::scr
contains

   subroutine wimsr
   !--------------------------------------------------------------------
   !
   ! Format multigroup cross sections from groupr for WIMS.
   !
   !---input specifications (free format)----
   !
   ! card 1
   !    ngendf  unit for input gendf tape
   !    nout    unit for output wims tape
   !
   ! card 2
   !    iprint  print option
   !               0=minimum (default)
   !               1=regular
   !               2=1+intermediate results
   !    iverw   wims version
   !               4=wims-d (default)
   !               5=wims-e
   !    igroup  group option
   !               0=69 groups (default)
   !               9=user's choice
   !
   ! card 2a  (igroup.eq.9 only)
   !    ngnd    number of groups
   !    nfg     number of fast groups
   !    nrg     number of resonance groups
   !    igref   reference group (default is last fast group)
   !
   ! card 3
   !    mat     endf mat number of the material to be processed
   !    nfid    not used
   !    rdfid   identification of material for the wims library
   !    iburn   burnup data option
   !              -1=suppress printout of burnup data
   !               0=no burnup data is provided (default)
   !               1=burnup data is provided in cards 5 and 6
   !
   ! card  4
   !    ntemp   no. of temperatures to process for this material
   !            in the thermal energy range
   !               (0=all found on input tape)
   !    nsigz   no. of sigma zeroes to process for this material
   !               (0=all found on input tape)
   !    sgref   reference sigma zero
   !               (.ge. 1.e10 to select all cross sect. at inf.dil.*
   !                           but fully shielded elastic x-sect,
   !                .lt. 1.e10 to select all x-sect at inf.dil.
   !                =sig0      from the list on groupr input to
   !                           select all x-sect. at that sig0)
   !    ires    resonance absorber indicator
   !               0=no resonance tables
   !              >0=ires temperatures processed
   !    sigp    potential cross section from endf.
   !               (if zero, replace by the elastic cross section)
   !    mti     thermal inelastic mt (default=0=none)
   !    mtc     thermal elastic mt (default=0=none)
   !    ip1opt  include p1 matrices
   !               0=yes
   !               1=no, correct p0 ingroups (default)
   !    inorf   resonance fission (if found)
   !               0=include resonance fission (default)
   !               1=do not include
   !    isof    fission spectrum
   !               0=do not include fission spectrum (default)
   !               1=include fission spectrum
   !    ifprod  fission product flag
   !               0=not a fission product (default)
   !               1=fission product, no resonance tables
   !               2=fission product, resonance tables
   !    jp1     transport correction neutron current spectrum flag
   !               0=use p1-flux for transport correction (default)
   !              >0=read in jp1 values of the neutron current
   !                   spectrum from input
   !
   !   the following cards 5 and 6 are for iburn gt 0 only
   ! card  5
   !    ntis    no. of time-dependent isotopes
   !            for burnable materials ntis=2
   !            for fissile materials ntis>2 when fission product
   !            yields are given.
   !    efiss   energy released per fission
   !
   ! card  6a
   !    identa  ident of capture product isotope
   !    yield   yield of product identa from capture
   !
   ! card  6b
   !    identa  ident of decay product isotope (zero if stable)
   !    lambda  decay constant (s-1)
   !
   ! card  6c  (repeated ntis-2 times, if necessary)
   !    identa  ident of fission product isotope
   !    yield   fission yield of identa from burnup of mat
   !
   ! card  7
   !    lambda  resonance-group goldstein lambdas (13 for
   !            default 69-group structure, nrg otherwise).
   !
   ! card  8    (only when jp1>0)
   !    p1flx   current spectrum (jp1 entries read, the rest are
   !            set with the default p1-flux calculated by njoy).
   !--------------------------------------------------------------------
   use util ! provides timer,openz,repoz,error
   use mainio ! provides nsysi,nsyso,nsyse
   integer::igroup,jp1,jcc,ij,i,j,ntis
   real(kr)::sec,efiss

   call timer(sec)
   write(nsyso,&
     '(/'' wimsr...create data for wims'',40x,f8.1,''s'')') sec
   write(nsyse,'(/'' wimsr...'',60x,f8.1,''s'')') sec
   nscr0=10
   nscr1=11
   nscr2=12
   nscr3=13
   nscr4=14

   !--read and write user input.
   read(nsysi,*) ngendf,nout
   nout=iabs(nout)
   call openz(ngendf,0)
   iprint=0
   iverw=4
   igroup=0
   read(nsysi,*) iprint,iverw,igroup
   write(nsyso,'(/&
     &'' input gendf unit ..................... '',i10/&
     &'' output unit .......................... '',i10/&
     &'' print option ......................... '',i10/&
     &'' wims version ......................... '',i10/&
     &'' group structure ...................... '',i10)')&
     ngendf,nout,iprint,iverw,igroup
   if (igroup.eq.0) then
      ngnd=69
      nfg=14
      nrg=13
      igref=nfg
   else if (igroup.eq.9) then
      igref=0
      read(nsysi,*) ngnd,nfg,nrg,igref
      if (igref.eq.0) igref=nfg
   endif
   write(nsyso,'(/&
     &'' ngnd ................................. '',i10/&
     &'' nfg .................................. '',i10/&
     &'' nrg .................................. '',i10/&
     &'' igref ................................ '',i10)')&
     ngnd,nfg,nrg,igref
   ixs=1
   call openz(nout,1)
   call repoz(nout)
   allocate(scr(nwscr))

   !--input material data
   iburn=0
   read(nsysi,*) mat,nfid,rdfid,iburn
   nfid=nint(rdfid)
   write(nsyso,'(/&
     &'' endf mat number ...................... '',i10/&
     &'' wims material identifier ............. '',i10/&
     &'' wims resonance identifier ............ '',f12.1/&
     &'' burn data (0=no, 1=yes, def=0) ....... '',i10)')&
     mat,nfid,rdfid,iburn
   ip1opt=1
   inorf=0
   isof=0
   ifprod=0
   jp1=0
   read(nsysi,*) ntemp,nsigz,sgref,ires,sigp,mti,mtc,&
     ip1opt,inorf,isof,ifprod,jp1
   if (ifprod.gt.0) then
      ifprod=1
      if (ires.gt.0) ifprod=2
   endif
   write(nsyso,'(/&
     &'' no. temperatures (thermal)............ '',i10/&
     &'' no. sigma zeroes ..................... '',i10/&
     &'' reference sigma zero ................. '',1p,e10.2/&
     &'' resonance absorber (0=no, >0=no.temp.) '',i10/&
     &'' pot. scatt. cross section ............ '',0p,f10.2/&
     &'' thermal inelastic mt ................. '',i10/&
     &'' thermal elastic mt ................... '',i10/&
     &'' p1 matrix option (0=yes, 1=no) ....... '',i10/&
     &'' resonance fission option (0=yes, 1=no) '',i10/&
     &'' fission spectrum option (0=no, 1=yes)  '',i10/&
     &'' fission product indicator ............ '',i10/&
     &'' current spectrum indicator ........... '',i10)')&
     ntemp,nsigz,sgref,ires,sigp,mti,mtc,&
     ip1opt,inorf,isof,ifprod,jp1

   !--input burn data
   yield(1)=0
   ifisp(1)=nfid
   if (iburn.gt.0) then
      read(nsysi,*) ntis,efiss
      if (ntis.gt.nymax) call error('wimsr',&
        'too many time-dependent isotopes',' ')
      jcc=ntis*2+4
      yield(4)=efiss
      ifisp(4)=0
      ij=1
      do i=1,ntis
         if (i.eq.1.or.i.eq.3) ij=ij+1
         read(nsysi,*) ifisp(ij),yield(ij)
         ij=ij+1
      enddo
   endif

   !--input resonance lambdas
   allocate(glam(nrg))
   read(nsysi,*) (glam(j),j=1,nrg)

   !--input current spectrum
   allocate(p1flx(ngnd))
   do j=1,ngnd
      p1flx(j)=0
   enddo
   if (jp1.ge.1) then
      read(nsysi,*) (p1flx(j),j=1,jp1)
   endif

   !--read in the cross section data
   call wminit

   !--process cross sections.
   call xsecs

   !--process effective resonance integrals.
   call resint

   !--process p1 scattering matrix.
   call p1scat

   !--write output tape.
   if (iburn.eq.0) jcc=2
   if (nout.ne.0) call wimout(jcc)

   !--wimsr is finished.
   if (allocated(snu)) deallocate(snu)
   if (allocated(abs2)) deallocate(abs2)
   if (allocated(spot)) deallocate(spot)
   deallocate(flux)
   deallocate(uff)
   deallocate(egb)
   deallocate(p1flx)
   deallocate(glam)
   deallocate(scr)
   call repoz(nout)
   call closz(nout)
   call closz(nscr0)
   call closz(nscr1)
   call closz(nscr2)
   call closz(nscr3)
   call closz(nscr4)
   call timer(sec)
   write(nsyso,&
     '(69x,f8.1,''s''/1x,7(''**********''),''*******'')') sec
   return
   end subroutine wimsr

   subroutine wminit
   !--------------------------------------------------------------------
   ! Store materials from ngendf.
   !--------------------------------------------------------------------
   use util ! provides repoz,mess
   use endf ! provides endf routines and variables
   use mainio ! provides nsysi,nsyso,nsyse
   ! internals
   integer::ngnd1,nb,nw,iza,i,nz,ntw,i1,ngn
   integer::i318,mess3,i618,mess6,i252
   integer::jtemp,nwflx
   character(60)::strng

   !--read gendf tape for groups and flags
   ngnd1=ngnd+1
   allocate(egb(ngnd1))
   call repoz(ngendf)
   call tpidio(ngendf,0,0,scr,nb,nw)
   call contio(ngendf,0,0,scr,nb,nw)
   if (math.ne.mat) call error('wminit',&
     'desired material is not on gendf tape',' ')
   iza=nint(c1h)
   iznum=iza/1000
   awr=c2h
   allocate(uff(ngnd))
   do i=1,ngnd
      uff(i)=0
   enddo
   nz=l2h
   ntw=n2h
   call listio(ngendf,0,0,scr,nb,nw)
   tempr(1)=c1h
   i1=7+ntw+nz
   ngn=l1h
   if (ngn.ne.ngnd) call error('wminit',&
     'incorrect group structure',' ')
   do i=1,ngnd1
      egb(i)=scr(ngnd1-i+i1)
   enddo
   ! check for presence of mf3, mt252, mf3, mt18 and mf6, mt18.
   i318=0
   mess3=0
   i618=0
   mess6=0
   i252=0
  125 continue
   call tosend(ngendf,0,0,scr)
  130 continue
   call contio(ngendf,0,0,scr,nb,nw)
   if (math.eq.0) go to 150
   if (mfh.eq.0.or.mth.eq.0) go to 130
   if (mfh.eq.3.and.mth.ge.18.and.mth.le.38) go to 135
   if (mfh.eq.6.and.mth.ge.18.and.mth.le.38) go to 140
   if (mfh.eq.3.and.mth.eq.252) i252=1
   go to 125
  135 continue
   if (mth.gt.21.and.mth.lt.38) go to 125
   if (mth.eq.18) i318=1
   if (mth.eq.18) go to 125
   if (i318.eq.0) go to 125
   if (mess3.gt.0) go to 125
   write(strng,'(''mat '',i4,'' mf '',i2,&
     &'' has both mt18 and mt'',i2)') math,mfh,mth
   call mess('wminit',strng,'mt18 will be used')
   mess3=1
   if (mess3.gt.0.and.mess6.gt.0) go to 150
   go to 125
  140 continue
   if (mth.gt.21.and.mth.lt.38) go to 125
   if (mth.eq.18) i618=1
   if (mth.eq.18) go to 125
   if (i618.eq.0) go to 125
   if (mess6.gt.0) go to 125
   write(strng,'(''mat '',i4,'' mf '',i2,&
     &'' has both mt18 and mt'',i2)') math,mfh,mth
   call mess('wminit',strng,'mt18 will be used')
   mess6=1
   if (mess3.gt.0.and.mess6.gt.0) go to 150
   go to 125
  150 continue
   if (i252.ne.1) then
     write(strng,'(''mat '',i4,'' has no mf3, mt252 '')') mat
     call mess('wminit',strng,&
       'isotropic c.m.scattering will be assumed')
   endif

   !--check the number of temperatures and sigma zeroes
   jtemp=1
   math=mat
   do while (math.eq.mat)
      call contio(ngendf,0,0,scr,nb,nw)
      if (math.eq.mat) then
         jtemp=jtemp+1
         call listio(ngendf,0,0,scr,nb,nw)
         tempr(jtemp)=c1h
         call tomend(ngendf,0,0,scr)
      endif
   enddo
   if (ntemp.eq.0) ntemp=jtemp
   if (ntemp.gt.jtemp) ntemp=jtemp
   if (ires.gt.0 .and. ntemp.lt.ires) ires=jtemp
   if (nsigz.eq.0) nsigz=nz
   if (nsigz.gt.nz) nsigz=nz
   write(nsyso,'(/'' processing'',i3,'' temperatures''/&
     &'' processing'',i3,'' sigma zeroes'')') ntemp,nsigz
   nwflx=ntemp*nsigz*nrg
   if(ires.gt.ntemp) nwflx=ires*nsigz*nrg
   allocate(flux(nwflx))

   !--print out group structure
   if (iprint.gt.0) write(nsyso,'(/&
     &'' group structure''/'' ---------------''/(1p,6e12.4))')&
     (egb(i),i=1,ngnd1)
   return
   end subroutine wminit

   subroutine resint
   !--------------------------------------------------------------------
   ! Process effective resonance integrals.
   !--------------------------------------------------------------------
   use util ! provides repoz,error
   use endf ! provides endf routines and variables
   use mainio ! provides nsysi,nsyso,nsyse
   ! internals
   integer::ntsr,nwflxr,nwelas,nwfa,nghi,nglo,nb,nw
   integer::i,jg,is,jtemp,jfiss,nl,nz,ntw,l,ig,kg,lim,iadd
   integer::loca,jz,loc,locn,jtem,iterm,ioff,it,index
   integer::indexl,iz
   real(kr)::xid,siglam,sigb,siga,sig
   real(kr),dimension(:),allocatable::sabs
   real(kr),dimension(:),allocatable::snux
   real(kr),dimension(:),allocatable::snsf
   real(kr),dimension(:),allocatable::flxr
   real(kr),dimension(:),allocatable::sigz
   real(kr),dimension(:),allocatable::abs
   real(kr),dimension(:),allocatable::elas
   real(kr),dimension(:),allocatable::fa
   integer,parameter::lz=6
   real(kr),parameter::zero=0

   !--allocate storage.
   if (ires.eq.0) go to 510
   write(nsyso,'(/'' ***resonance integrals***'')')
   ntsr=ires*nsigz*nrg
   allocate(sabs(ntsr))
   allocate(snux(ntsr))
   allocate(snsf(ntsr))
   nwflxr=ires*nsigz
   allocate(flxr(nwflxr))
   allocate(sigz(nsigz))
   allocate(abs(ires))
   nwelas=ires*nsigz*nrg
   allocate(elas(nwelas))
   nwfa=ires*nsigz
   allocate(fa(nwfa))
   nghi=ngnd-nfg
   nglo=nghi-nrg+1

   !--read gendf data
   call repoz(ngendf)
   call tpidio(ngendf,0,0,scr,nb,nw)
   do i=1,nwflxr
      flxr(i)=0
   enddo
   do jg=1,nrg
      do is=1,nwflxr
         i=(jg-1)*nwflxr+is
         sabs(i)=abs2(nfg+jg)
         snux(i)=0
         snsf(i)=0
         elas(i)=0
         flux(i)=0
      enddo
   enddo
   jtemp=0
   jfiss=1

   !--loop over temperatures
  130 continue
   call contio(ngendf,0,0,scr,nb,nw)
   if (math.lt.0) go to 400
   nl=l1h
   nz=l2h
   ntw=n2h
   l=1
   call listio(ngendf,0,0,scr(l),nb,nw)
   jtemp=jtemp+1
   do i=1,nsigz
      sigz(nsigz-i+1)=scr(6+ntw+i)
   enddo
   abs(jtemp)=0
   call tosend(ngendf,0,0,scr(l))

   !--loop over reactions
  310 continue
   call contio(ngendf,0,0,scr,nb,nw)
   if (math.eq.0) go to 320
   if (mfh.eq.0) go to 310
   if (mth.eq.0) go to 310
   xid=mat
   nl=l1h
   nz=l2h
   ntw=n2h

   !--loop over groups for this reaction
  140 continue
   l=1
   call listio(ngendf,0,0,scr(l),nb,nw)
   if (mth.ge.18.and.mth.le.21.and.inorf.gt.0) go to 300
   if (mth.eq.38.and.inorf.gt.0) go to 300
   if (mfh.ne.3) go to 300
   if (mth.eq.1) go to 162
   if (mth.eq.2) go to 162
   if (mth.lt.18.or.mth.gt.150) go to 300
   if (mth.gt.21.and.mth.lt.102.and.mth.ne.38) go to 300
  162 continue
   ig=n2h
   do while (nb.ne.0)
      l=l+nw
      call moreio(ngendf,0,0,scr(l),nb,nw)
      if (l+nw-1.gt.nwscr) call error('resint',&
        'storage exceeded',' ')
   enddo
   if (mth.eq.1.and.ig.eq.ngnd-igref+1) go to 172
   if (ig.ge.nglo.and.ig.le.nghi) go to 172
   go to 290
  172 continue
   kg=ngnd-ig+1
   jg=kg-nfg
   lim=nsigz
   if (nsigz.gt.nz) lim=nz
   if (mth.ne.1.and.mth.ne.2.and.mth.ne.18.and.mth.ne.102) go to 300

   !--absorption
   if (mth.eq.102) then
      iadd=nsigz+nsigz*(jtemp-1+ires*(jg-1))
      loca=l+lz+nl*(nz-1)
      do jz=1,lim
         sabs(iadd-jz+1)=sabs(iadd-jz+1)+scr(nl*jz+loca)
      enddo
      abs(jtemp)=1

   !--flux
   else if (mth.eq.1) then
      if (jg.ge.1.and.jg.le.nrg) then
         loc=1+nsigz+nsigz*(jtemp-1+ires*(jg-1))
         loca=l+lz-nl
         do jz=1,lim
            flux(loc-jz)=scr(nl*jz+loca)
         enddo
      endif
      ! reference group flux
      if (kg.eq.igref) then
         do jz=1,lim
            loca=l+lz+nl*(jz-1)
            loc=1+nsigz-jz+nsigz*(jtemp-1)
            flxr(loc)=scr(loca)
         enddo
      endif

   !--fission
   else if (mth.eq.18) then
      iadd=nsigz+nsigz*(jtemp-1+ires*(jg-1))
      loca=l+lz+nl*(nz-1)
      do jz=1,lim
         snsf(iadd-jz+1)=snsf(iadd-jz+1)+scr(nl*jz+loca)
      enddo
      jfiss=3
      do jz=1,lim
         sabs(iadd-jz+1)=sabs(iadd-jz+1)+scr(nl*jz+loca)
      enddo

   !--elastic
   else if (mth.eq.2) then
      iadd=nsigz+nsigz*(jtemp-1+ires*(jg-1))
      loca=l+lz+nl*(nz-1)
      do jz=1,lim
         elas(iadd-jz+1)=elas(iadd-jz+1)+scr(nl*jz+loca)
      enddo
   endif

   !--continue group loop
  290 continue
   if (ig.lt.ngnd) go to 140
  300 continue
   call tosend(ngendf,0,0,scr)

   !--continue loop over reactions
   go to 310

   !--continue loop over temperatures
  320 continue
   if (jtemp.ge.ires) go to 400
   go to 130

   !--desired data is loaded.
   !--calculate nu times sigma f
  400 continue
   if (jfiss.ne.1) then
      do jg=1,nrg
         locn=nfg+jg
         do jz=1,nsigz
            do jtem=1,ires
               iterm=jz-1+nsigz*ires*(jg-1)-nsigz
               ioff=nsigz*jtem+iterm
                  snux(ioff+1)=snsf(ioff+1)*snu(locn)
               enddo
         enddo
      enddo
   endif

   !--check flux arrays for temperature dependence.
   do jz=1,nsigz
      do it=1,ires
         if (flxr(jz+nsigz*(it-1)).eq.zero) then
            flxr(jz+nsigz*(it-1))=flxr(jz+nsigz*(it-2))
         endif
      enddo
   enddo
   do jg=1,nrg
      do jz=1,nsigz
         do it=1,ires
            index=1+jz-1+nsigz*(it-1+ires*(jg-1))
            if (flux(index).eq.zero) then
               indexl=1+jz-1+nsigz*(it-2+ires*(jg-1))
               flux(index)=flux(indexl)
            endif
         enddo
      enddo
   enddo

   !--correct non-temperature dependent absorption if necessary
   if (ires.ne.1) then
      do it=1,ires
         if (abs(it).le.zero) then
            do jg=1,nrg
               loc=1+nsigz*(it-2+ires*(jg-1))
               loca=loc+nsigz
               do jz=1,nsigz
                  sabs(jz+loca)=sabs(jz+loca)+sabs(jz+loc)
               enddo
            enddo
         endif
      enddo
   endif

   !--convert njoy cross sections to resonance integrals
   do jg=1,nrg
      siglam=spot(nfg+jg)*glam(jg)
      do it=1,ires
         loc=1+nsigz*(it-1+ires*(jg-1))
         do iz=1,nsigz
            sigb=sigz(iz)+siglam
            siga=sabs(loc+iz-1)
            sig=snux(loc+iz-1)
            sabs(loc+iz-1)=sigb*siga/(sigb+siga)
            snux(loc+iz-1)=sigb*sig/(sigb+siga)
         enddo
      enddo
   enddo

   !--write out results
   call rsiout(xid,jfiss,fa,sabs,snux,flxr,sigz,elas)

   !--resint is finished.
  510 continue
   return
   end subroutine resint

   subroutine rsiout(xid,ifiss,b,sabs,snux,flxr,sigz,elas)
   !--------------------------------------------------------------------
   ! Write resonance integrals.
   !--------------------------------------------------------------------
   use util ! provides openz,repoz
   use mainio ! provides nsysi,nsyso,nsyse
   ! externals
   real(kr)::xid
   integer::ifiss
   real(kr)::b(*),sabs(*),snux(*),flxr(*),sigz(*),elas(*)
   ! internals
   integer::nscr,ncol1,lim,ntnp,jg,loc,j,jtemp,n,m,i,ig
   integer::it,locf,locr,jz,il,nb,loca
   real(kr)::siglam,refflx
   real(kr)::c(15)
   character(13),parameter::undl='-------------'
   integer,parameter::nmax=5
   integer,parameter::ncol=8

   nscr=nscr1
   call openz(-nscr,1)
   call repoz(-nscr)
   ncol1=ncol+1
   lim=nsigz
   if (nsigz.gt.ncol) lim=ncol
   if (ifiss.eq.5) ifiss=3
   ntnp=ires*nsigz
   write(nscr) ntnp
   write(nscr) xid,ires,nsigz
   do jg=1,nrg
      siglam=spot(nfg+jg)*glam(jg)
      loc=1+ires*nsigz*(jg-1)
      write(nscr) (tempr(j),j=1,ires),&
        (sigz(j)+siglam,j=1,nsigz),&
        ((sabs(loc-1+nsigz*(jtemp-1)+j),j=1,nsigz),jtemp=1,ires)
      if (ifiss.eq.3) then
         loc=1+ires*nsigz*(jg-1)
         write(nscr) (tempr(j),j=1,ires),&
           (sigz(j)+siglam,j=1,nsigz),&
           ((snux(loc-1+nsigz*(jtemp-1)+j),j=1,nsigz),jtemp=1,ires)
      endif
      if (iverw.ne.4) then
         loc=1+ires*nsigz*(jg-1)
         write(nscr) (tempr(j),j=1,ires),&
           (sigz(j)+siglam,j=1,nsigz),&
           ((elas(loc-1+nsigz*(jtemp-1)+j),j=1,nsigz),jtemp=1,ires)
      endif
   enddo
   n=nsigz
   if (n.gt.nmax) n=nmax
   m=n
   if (m.lt.4) m=4
   write(nsyso,'(/&
     &'' temperatures (down)''/&
     &'' -------------------'')')
   write(nsyso,'(1x,1p,9e11.3)') (tempr(i),i=1,ires)
   write(nsyso,'(/&
     &'' potential cross sections (across)''/&
     &'' ---------------------------------'')')
   write(nsyso,'(1x,1p,9e11.3)') (sigz(i),i=1,nsigz)
   if (iprint.ge.2) then
      write(nsyso,'(/&
        &'' group'',6x,&
        &''flux per unit lethargy normalized at group '',i3/&
        &'' ----------'',6(a))')&
        igref,(undl,i=1,m)
      do ig=1,nrg
         do it=1,ires
            locf=1-1+nsigz*(it-1+ires*(ig-1))
            locr=1-1+nsigz*(it-1)
            do jz=1,nsigz
               refflx=flxr(jz+locr)/log(egb(igref)/egb(igref+1))
               c(jz)=flux(jz+locf)/log(egb(nfg+ig)/egb(nfg+ig+1))
               c(jz)=c(jz)/refflx
            enddo
            jg=ig+nfg
            if (it.eq.1) then
               write(nsyso,'(/i5,5x,1p,9e13.5)') jg,(c(i),i=1,lim)
            endif
            if (it.eq.1.and.nsigz.gt.ncol) then
               write(nsyso,'(10x,1p,9e13.5)') (c(i),i=ncol1,nsigz)
            endif
            if (it.gt.1) then
               write(nsyso,'(10x,1p,9e13.5)') (c(i),i=1,nsigz)
            endif
         enddo
      enddo
   endif

   if (iprint.eq.0) return
   il=1
   write(nsyso,'(/&
     &'' group'',6x,''absorption''/&
     &'' ----------'',5(a))')&
     (undl,i=1,n)
   if (il.ne.2.or.ifiss.ge.2) then
      do ig=1,nrg
         nb=0
         jg=ig+nfg
         do it=1,ires
            loca=nsigz*(it-1+ires*(ig-1))
            do jz=1,nsigz
               nb=nb+1
               b(nb)=sabs(jz+loca)
               c(jz)=b(nb)
            enddo
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
         enddo
      enddo
   endif
   il=2
   if (il.ne.2.or.ifiss.ge.2) then
      write(nsyso,'(/&
        &'' group'',6x,''fission yield''/&
        &'' ----------'',5(a))')&
       (undl,i=1,n)
      do ig=1,nrg
         nb=0
         jg=ig+nfg
         do it=1,ires
            loca=nsigz*(it-1+ires*(ig-1))
            do jz=1,nsigz
               nb=nb+1
               b(nb)=snux(jz+loca)
               c(jz)=b(nb)
            enddo
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
         enddo
      enddo
   endif

   if (iverw.eq.4) return
   write(nsyso,'(/&
     &'' group      shielded total scattering cross sections''/&
     &'' ----------'',5(a))')&
     (undl,i=1,n)
   do ig=1,nrg
      nb=0
      jg=ig+nfg
      do it=1,ires
         loca=1-1+nsigz*(it-1+ires*(ig-1))
         do jz=1,nsigz
            nb=nb+1
            b(nb)=elas(jz+loca)
            c(jz)=b(nb)
         enddo
         if (it.eq.1) then
            write(nsyso,'(/i5,5x,1p,9e13.5)') jg,(c(i),i=1,lim)
         endif
         if (it.eq.1.and.nsigz.gt.ncol) then
            write(nsyso,'(10x,1p,9e13.5)') (c(i),i=ncol1,nsigz)
         endif
         if (it.gt.1) then
            write(nsyso,'(10x,1p,9e13.5)') (c(i),i=1,lim)
         endif
         if (it.gt.1.and.nsigz.gt.ncol) then
            write(nsyso,'(10x,1p,9e13.5)') (c(i),i=ncol1,nsigz)
         endif
      enddo
   enddo
   return
   end subroutine rsiout

   subroutine xsecs
   !--------------------------------------------------------------------
   ! Process cross sections.
   !--------------------------------------------------------------------
   use util ! provides repoz,error,mess,skiprz
   use endf ! provides endf routines and variables
   use mainio ! provides nsysi,nsyso,nsyse
   character(60)::strng
   integer::ngndsq,il,ngr1,nb,nw,jtemp,i318,jfisd,jfist,jfspd
   integer::jfspt,jfspp,jfiss,i,mf3mt2,nth1,nth,if6,matd,jn2n
   integer::nz,ntw,l,idone,iz,jic,nl,ng2,ig2lo,ig,jg,loca,loc
   integer::jg2,locc,ll,max,min,ig2,lt,locf,k,in,nthr,loc1,loc2,jfis
   real(kr)::alf,xxi,cnorm,dnorm,p1nrm
   integer,dimension(:),allocatable::l1,l2
   integer,dimension(:),allocatable::l1e,l2e
   real(kr),dimension(:),allocatable::abs1
   real(kr),dimension(:),allocatable::xtr
   real(kr),dimension(:),allocatable::sn2n
   real(kr),dimension(:),allocatable::scat
   real(kr),dimension(:),allocatable::xi
   real(kr),dimension(:),allocatable::sdp
   real(kr),dimension(:),allocatable::cspc
   real(kr),dimension(:),allocatable::csp1
   real(kr),dimension(:),allocatable::snus
   real(kr),dimension(:),allocatable::chi
   real(kr),dimension(:),allocatable::sfi
   real(kr),dimension(:),allocatable::sf0
   real(kr),dimension(:),allocatable::ab0
   real(kr),dimension(:),allocatable::xs
   real(kr),parameter::dilinf=1.e10_kr
   integer,parameter::lz=6
   real(kr),parameter::zero=0
   real(kr),parameter::one=1

   !--allocate storage.
   if (ixs.eq.0) go to 630
   write(nsyso,'(/'' ***cross sections***'')')
   allocate(snu(ngnd))
   allocate(spot(ngnd))
   allocate(abs2(ngnd))
   allocate(l1(ngnd))
   allocate(l2(ngnd))
   allocate(l1e(ngnd))
   allocate(l2e(ngnd))
   ngndsq=ngnd*ngnd
   allocate(xs(ngndsq))
   allocate(ab0(ngnd))
   allocate(sf0(ngnd))
   allocate(sfi(ngnd))
   allocate(chi(ngnd))
   allocate(snus(ngnd))
   allocate(cspc(ngnd))
   allocate(xtr(ngnd))
   allocate(abs1(ngnd))
   allocate(sn2n(ngnd))
   allocate(scat(ngnd))
   allocate(xi(ngnd))
   allocate(sdp(ngnd))
   allocate(csp1(ngnd))
   nfiss=0
   il=1
   ngr1=nfg+nrg
   p1nrm=1

   !--read gendf data
   call repoz(ngendf)
   call tpidio(ngendf,0,0,scr,nb,nw)
   isg=0
   if (sgref.lt.dilinf) isg=1
   jtemp=0
   i318=0
   jfisd=0
   jfist=0
   jfspd=0
   jfspt=0
   jfspp=0
   jfiss=0
   jn2n=0
   ! preset aver.log.decrement per collision assuming isotropic scatt.
   alf=(awr-1)/(awr+1)
   alf=alf*alf
   xxi=1+log(alf)*alf/(1-alf)

   !--loop over temperatures
  140 continue
   do i=1,ngnd
      csp1(i)=0
      spot(i)=sigp
      xi(i)=xxi
      abs2(i)=0
      sf0(i)=0
      sfi(i)=0
      snus(i)=0
      cspc(i)=0
      chi(i)=0
      scat(i)=0
      sn2n(i)=0
      l1e(i)=ngnd
      l2e(i)=1
      l1(i)=ngnd
      l2(i)=1
   enddo
   do i=1,ngndsq
      xs(i)=0
   enddo
   mf3mt2=0
   cnorm=0
   dnorm=0
   nth1=0
   nth=0
   if6=0
   call contio(ngendf,0,0,scr,nb,nw)
   if (math.ne.mat) go to 540
   matd=mat
   nz=l2h
   ntw=n2h
   l=1
   call listio(ngendf,0,0,scr(l),nb,nw)
   idone=0
   i=0
   do while (i.lt.nz.and.idone.eq.0)
      i=i+1
      if (abs(sgref-scr(l+5+ntw+i)).le.sgref/100) then
         iz=i
         idone=1
      endif
   enddo
   if (idone.eq.0) then
      write(strng,&
        '(''ref. sig0'',1p,e10.3,'' not on the list'')') sgref
      call mess('xsecs ',strng,'first entry used as default')
      iz=1
      sgref=scr(l+5+ntw+iz)
   endif
   if(isg.gt.0) isg=iz
   jtemp=jtemp+1
   call tosend(ngendf,0,0,scr(l))
   jic=0

   !--loop over reactions
  150 continue
   call contio(ngendf,0,0,scr(l),nb,nw)
   if (math.eq.0) go to 370
   if (mfh.eq.0) go to 150
   if (mth.eq.0) go to 150
   nl=l1h
   nz=l2h
   if (mfh.eq.3.and.mth.ge.2.and.mf3mt2.le.0) mf3mt2=1

   !--loop over groups
  160 continue
   l=1
   call listio(ngendf,0,0,scr(l),nb,nw)
   ng2=l1h
   ig2lo=l2h
   ig=n2h
   do while (nb.ne.0)
      l=l+nw
      call moreio(ngendf,0,0,scr(l),nb,nw)
      if (l+nw-1.gt.nwscr) call error('xsecs',&
        'storage exceeded',' ')
   enddo
   jg=ngnd-ig+1
   l=1

   !--branch according to type of data
   if (mth.eq.452) go to 245
   if (mth.eq.455) go to 252
   if (mfh.ne.3) go to 192
   if (mth.ge.102.and.mth.le.150) go to 205
   if ((mth.ge.18.and.mth.le.21).or.mth.eq.38) go to 214
   if (mth.eq.1) go to 220
   if (mth.eq.2) go to 240
   if (mth.eq.16) jn2n=16
   if (mth.eq.16.or.mth.eq.24) go to 236
   if ((mth.ge.875.and.mth.le.891).and.(jn2n.ne.16)) go to 236
   if (mth.eq.17.or.mth.eq.25) go to 237
   if (mth.eq.252) go to 242
   if (mth.eq.mti) go to 238
  192 continue
   if (mfh.ne.6) go to 365
   if (mth.eq.2.or.mth.eq.mti.or.mth.eq.mtc) go to 270
   if ((mth.lt.18.or.mth.gt.21).and.mth.ne.38) go to 300
   go to 300

   !--absorption cross section
  205 continue
   loca=l+lz+nl*nz
   if (mth.eq.102) then
      if (isg.gt.0.and.nz.ge.iz) loca=l+lz+nl*(iz-1+nz)
      abs1(jg)=scr(loca)
   else
      abs2(jg)=abs2(jg)+scr(loca)
   endif
   go to 360

   !--fission cross section
  214 continue
   if (mth.eq.18) i318=1
   if (mth.eq.18) go to 216
   if (i318.eq.0) go to 216
   go to 365
  216 continue
   loca=l+lz+nl*nz
   sfi(jg)=sfi(jg)+scr(loca)
   if (isg.gt.0.and.nz.ge.iz) loca=l+lz+nl*(iz-1+nz)
   sf0(jg)=sf0(jg)+scr(loca)
   ab0(jg)=ab0(jg)+scr(loca)
   go to 360

   !--p1-flux for transport correction (if not read in from input)
   !--normalize the input current to the first common njoy p1-flux
  220 continue
   loca=l+lz+nl*(nz-1)
   if (isg.gt.0.and.nz.ge.iz) loca=l+lz+nl*(iz-1)
   if (nl.gt.1) loca=loca+1
   if (p1flx(jg).eq.0) then
     p1flx(jg)=scr(loca)
     p1nrm=1
   else
     if (p1nrm.eq.one) p1nrm=scr(loca)/p1flx(jg)
     p1flx(jg)=p1flx(jg)*p1nrm
   endif
   go to 360

   !--n2n cross section
  236 continue
   loca=l+lz+nl*nz
   sn2n(jg)=sn2n(jg)+scr(loca)
   go to 360

   !--n3n, n3n-alpha cross section
  237 continue
   loca=l+lz+nl*nz
   sn2n(jg)=sn2n(jg)+2*scr(loca)
   go to 360

   !--incoherent thermal scattering
  238 continue
   if (ig.gt.nth1.and.scr(l+lz+nz*nl).ne.zero) nth1=ig
   nth=ngnd-nth1+1
   go to 360

   !--elastic cross section
  240 continue
   go to 360

   !--log slowing down
  242 continue
   loca=l+lz+nl*nz
   xi(jg)=scr(loca)
   go to 360

   !--fission nu and chi from 3/452 and 5/452
  245 continue
   if (jtemp.gt.1) go to 365
   if (jfiss.eq.0) jfiss=1
   if (mfh.eq.5) go to 246
   ! total fission yield
   jfist=1
   loca=l+lz+1
   snu(jg)=scr(loca)
   go to 360
  246 continue
   if (isof.ne.0) then
      jfspt=1
      do i=2,ng2
         jg2=ngnd-ig2lo-i+3
         if (jg2.ge.1.and.jg2.le.ngnd) then
            loca=l+lz+nl*(i-1)
            uff(jg2)=scr(loca)
         endif
      enddo
   endif
   go to 365

   !--delayed fission contributions from 3/455 and 5/455
  252 continue
   if (jtemp.gt.1) go to 365
   if (mfh.ne.5) then
      jfisd=1
      loca=l+lz
      snus(jg)=snus(jg)+scr(loca+1)*scr(loca+2)
      dnorm=dnorm+scr(loca)*scr(loca+1)*scr(loca+2)
   else
      jfspd=1
      do i=2,ng2
         locc=1+ngnd-ig2lo-i+2
         do ll=1,nl
            loca=l+lz+ll-1+nl*(i-1)
            chi(locc)=chi(locc)+dnorm*scr(loca)
         enddo
      enddo
   endif
   go to 360

   !--temperature-dependent file 6 matrices
  270 continue
   if (mth.ne.2) jic=1
   max=0
   min=ngnd
   if (mth.eq.2.and.jg.ge.nth) go to 295
   if (mth.eq.2.and.jg.lt.nth) go to 301
   do i=2,ng2
      ig2=ig2lo+i-2
      if (ig2.ge.1.and.ig2.le.ngnd) then
         jg2=ngnd-ig2+1
         loca=l+lz+(il-1)+nl*nz*(i-1)
         if (isg.gt.0.and.nz.ge.iz) loca=loca+nl*(iz-1)
         if (scr(loca).ne.zero) then
            if (mth.ne.2) then
               loc=1+jg-1+ngnd*(jg2-1)
               if (jg2-jg.ge.-1.or.scr(loca).ge.zero) then
                  xs(loc)=xs(loc)+scr(loca)
                  scat(jg)=scat(jg)+scr(loca)
                  ! for the thermal groups
                  ! accumulate csp1,
                  ! the sum of the p1 components of scattering
                  ! from group i to group j for all j
                  if (jg.gt.ngr1) then
                    csp1(jg)=csp1(jg)+scr(loca+1)
                  endif
                  if (jg2.gt.max) max=jg2
                  if (jg2.lt.min) min=jg2
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
   if (ig.lt.ngnd) go to 160
   go to 365

   !--temperature-independent file 6 matrices
  300 continue
   if (jtemp.gt.1) go to 360
  301 continue
   if (mfh.ne.6) go to 365
   if (if6.gt.0) go to 305
   if (nth1.eq.0) nth1=ngnd-nfg-nrg
   nth=ngnd-nth1+1
   if6=1
  305 continue
   if (mth.ge.18.and.mth.le.21) go to 315
   if (mth.eq.38) go to 315
   if (mth.eq.2.and.jg.ge.nth) go to 365
   if (mth.ge.221.and.mth.le.250) go to 365
   max=0
   min=ngnd
   do i=2,ng2
      ig2=ig2lo+i-2
      if (ig2.ge.1.and.ig2.le.ngnd) then
         jg2=ngnd-ig2+1
         loca=l+lz+(il-1)+nl*nz*(i-1)
         if (mth.eq.2) then
            loca=l+lz+(il-1)+nl*nz*(i-1)+(nz-1)*nl
            if (isg.gt.0.and.nz.ge.iz) then
               loca=l+lz+(il-1)+nl*nz*(i-1)+nl*(iz-1)
            endif
         endif
         loc=1+jg-1+ngnd*(jg2-1)
         xs(loc)=xs(loc)+scr(loca)
         scat(jg)=scat(jg)+scr(loca)
         if (jg2.gt.max) max=jg2
         if (jg2.lt.min) min=jg2
         if (nl.ne.1) then
            ! for fast and resonance groups
            ! accumulate csp1,
            ! the sum of the p1 components of scattering
            ! from group j to i for all j
            if (jg2.lt.nth) csp1(jg2)=&
              csp1(jg2)+scr(loca+1)*p1flx(jg)/p1flx(jg2)
         endif
      endif
   enddo
   lt=l1(jg)
   if (min.lt.lt) l1(jg)=min
   lt=l2(jg)
   if (max.gt.lt) l2(jg)=max
   go to 360

   !--process nu*sigf and chi from fission matrix
  315 continue
   if (jtemp.gt.1) go to 365
   if (ig.ne.0) then
      jfiss=2
      jfspp=1
      locf=l+lz
      do i=2,ng2
         loca=l+lz+nl*nz*(i-1)
         if (ig2lo.ne.0) then
            ! matrix part
            snus(jg)=snus(jg)+scr(loca)
            locc=1+ngnd-i-ig2lo+2
            chi(locc)=chi(locc)+scr(locf)*scr(loca)
         else
            ! production part
            do k=1,ngnd
               snus(jg)=snus(jg)+cspc(k)*scr(loca)
               locc=1+ngnd-k
               chi(locc)=chi(locc)+cspc(k)*scr(loca-1)*scr(loca)
            enddo
         endif
      enddo
   else
      ! spectrum part
      do i=1,ng2
         cspc(ig2lo+i-1)=scr(l+lz+nl*nz*(i-1))
      enddo
   endif

   !--continue loop over groups
  360 continue
   if (ig.lt.ngnd) go to 160
   go to 365

   !--continue loop over reactions
  365 continue
   call tosend(ngendf,0,0,scr(l))
   go to 150

   !--data is loaded for this temperature
  370 continue
   if (jtemp.le.1) then
      if (jfiss.ge.2) then
         if (jfisd.ge.1.or.jfist.le.0) then
            ! fix up nu if computed from matrix
            do i=1,ngnd
               if (sfi(i).ne.zero) then
                  snu(i)=snus(i)/sfi(i)
               endif
            enddo
            if (jfisd.ne.1) then
               write(strng,'(&
                 &''nu-bar calculated from fission matrix'')')
               call mess('xsecs',strng,&
                 'only prompt contribution available')
            endif
         endif
      endif
   endif
   if (jic.gt.0) go to 385
   jtemp=jtemp-1
   write(strng,&
     '(''use only '',i2,'' temps for mat '',i4)') jtemp,matd
   call mess('xsecs',strng,'mti missing from higher temps')
   matd=-1
   if (math.ne.0) call skiprz(ngendf,-3)
   go to 540
  385 continue
   if (jfiss.ne.0) then
      ! calculate nu*sigma f
      do i=1,ngnd
         snus(i)=snu(i)*sf0(i)
      enddo
   endif
   ! calculate total absorption
   do i=1,ngnd
      ab0(i)=sf0(i)&
        +abs1(i)+abs2(i)-sn2n(i)
   enddo
   ! correct scattering matrix
   in=nth-1
   do i=1,in
      if (l1e(i).lt.l1(i)) l1(i)=l1e(i)
      if (l2e(i).gt.l2(i)) l2(i)=l2e(i)
   enddo
   ! correct transport
   do i=1,ngnd
      if (ip1opt.gt.0)&
        xs(1+(i-1)*(1+ngnd))=xs(1+(i-1)*(1+ngnd))-csp1(i)
      xtr(i)=scat(i)+ab0(i)-csp1(i)
   enddo
   ! compute xx - slowing down power per unit lethargy
   do i=1,ngnd
      ! suppress upscattering from thermal into resonance groups
      nthr=nfg+nrg+1
      do jg=nthr,ngnd
         jg2=l1(jg)
         do while (jg2.lt.nthr)
            loc1=1+jg-1+ngnd*(jg2-1)
            loc2=1+jg-1+ngnd*(nthr-1)
            xs(loc2)=xs(loc2)+xs(loc1)
            xs(loc1)=0
            jg2=jg2+1
            l1(jg)=jg2
         enddo
      enddo
      sdp(i)=scat(i)*xi(i)/log(egb(i)/egb(i+1))
      ! replace constant potential with scattering cross section
      if (sigp.eq.zero) spot(i)=scat(i)
   enddo
   if (jtemp.eq.1) then
      jfis=0
      if (ires.gt.0) jfis=1
      if (jfis.ge.1.and.jfiss.gt.0) jfis=3
      if (jfis.gt.1.and.inorf.gt.0) jfis=2
      if (jfis.eq.0.and.jfiss.gt.0) jfis=4
   endif

   !--print out the results
   call xseco(jtemp,jfis,l1,l2,xs,sn2n,sdp,xtr,snus,sf0,ab0)
   i318=0
   ifiss=jfis
   if (isof.eq.0) go to 520
   if (jtemp.gt.1) go to 520
   if (jfiss.eq.0) go to 520
   ! normalize chi
   ! check for fission spectrum consistency
   if (jfspt.eq.0) then
      if (jfspd.ne.1.or.jfspp.ne.1) then
         write(strng,'(&
           &''spectrum calculated from fission matrix'')')
         call mess('xsecs ',strng, &
           'only prompt contribution available')
      endif
      do i=1,ngnd
         uff(i)=chi(i)
      enddo
   endif
   ! normalize the fission spectrum
   cnorm=0
   do i=1,ngr1
      if (uff(i).gt.zero) then
         cnorm=cnorm+uff(i)
         nfiss=i
      endif
   enddo
   if (cnorm.gt.zero) then
      cnorm=1/cnorm
      do i=1,nfiss
         uff(i)=uff(i)*cnorm
      enddo
   endif
  520 continue

   !--continue loop over temperatures
   if (jtemp.lt.ntemp) go to 140
   if (ntemp.eq.0) go to 140
  540 continue

   !--xsecs is finished.
   if (isof.ne.0.and.iprint.ne.0) then
      write(nsyso,'(/&
        &'' fission spectrum (groups 1 -'',i3,'')'',/&
        &'' --------------------------------'')')&
        nfiss
      write(nsyso,'(1p,6e12.4)') (uff(i),i=1,nfiss)
   endif

  630 return
   end subroutine xsecs

   subroutine xseco(itemp,ifiss,l1,l2,xs,sn2n,sdp,xtr,snus,sf0,ab0)
   !--------------------------------------------------------------------
   ! Write cross sections on nout and nsyso.
   !--------------------------------------------------------------------
   use util ! provides openz,error
   use mainio ! provides nsysi,nsyso,nsyse
   ! externals
   integer::itemp,ifiss
   integer::l1(*),l2(*)
   real(kr)::xs(*),sn2n(*),sdp(*),xtr(*),snus(*),sf0(*),ab0(*)
   ! internals
   integer::ngr0,ngr1,nnt,i,jtemp,p1nrm,k,ia,lone,ltwo
   integer::ks,loc,lim,j,ib,nt0
   integer,parameter::ncol=5
   real(kr),parameter::zero=0

   !--retrieve storage pointers
   ngr0=nfg+1
   ngr1=nfg+nrg
   nnt=ngr1

   !--temperature-independent data
   if (itemp.eq.1) then

      !--set up scratch units for xsec output
      call openz(-nscr2,1)
      call openz(-nscr3,1)

      !--write temperature-independent part
      !--what is the correct definition of xx?---
      if (iverw.eq.5) write(nscr2) (spot(i),i=ngr0,ngr1),&
        (sdp(i),i=ngr0,ngr1),(sn2n(i),i=1,nfg),&
        (xtr(i),i=1,nnt),&
        (ab0(i),i=1,nnt),(zero,i=ngr0,ngr1),&
        (glam(i),i=1,nrg)
      if (iverw.eq.4) write(nscr2) (spot(i),i=ngr0,ngr1),&
        (sdp(i),i=ngr0,ngr1),&
        (xtr(i),i=1,nnt),&
        (ab0(i),i=1,nnt),(zero,i=ngr0,ngr1),&
        (glam(i),i=1,nrg)
      if (iprint.ne.0) then
         jtemp=nint(tempr(1))
         p1nrm=1
         if (igref.le.nnt) p1nrm=1/p1flx(igref)
         write(nsyso,'(/&
           &'' neutron current spectrum  (groups 1-'',i3,'')''/&
           &'' ---------------------------------------'')') nnt
         write(nsyso,'(1p,6e12.4)') (p1nrm*p1flx(i),i=1,nnt)
         write(nsyso,'(/&
           &'' sigma potential (groups '',i3,''-'',i3,'')''/&
           &'' ------------------------------'')')&
           ngr0,ngr1
         write(nsyso,'(1p,6e12.4)') (spot(i),i=ngr0,ngr1)
         write(nsyso,'(/&
           &'' scattering power per unit lethargy (groups '',i3,&
           &''-'',i3,'')''/&
           &'' -------------------------------------------------'')')&
           ngr0,ngr1
         write(nsyso,'(1p,6e12.4)')(sdp(i),i=ngr0,ngr1)
         if (iverw.ne.4) then
            write(nsyso,'(/&
              &'' n,2n (groups 1-'',i3,'')''/&
              &''------------------'')')&
              nfg
            write(nsyso,'(1p,6e12.4)') (sn2n(i),i=1,nfg)
         endif
         write(nsyso,'(/&
           &'' transport corrected total (groups 1-'',i3,'')''/&
           &'' ---------------------------------------'')')&
           nnt
         write(nsyso,'(1p,6e12.4)') (xtr(i),i=1,nnt)
         write(nsyso,'(/&
           &'' absorption (groups 1-'',i3,'')''/&
           &'' ------------------------'')')&
           nnt
         write(nsyso,'(1p,6e12.4)') (ab0(i),i=1,nnt)
         write(nsyso,'(/&
           &'' goldstein-cohen intermediate resonance factors'',&
           &'' (groups '',i3,''-'',i3,'')''/&
           &1x,6(''----------''),''-'')')&
           ngr0,ngr1
         write(nsyso,'(1p,6e12.4)') (glam(i),i=1,nrg)
      endif
      if (ifiss.gt.1) then
         write(nscr2)&
           (snus(i),i=1,nnt),(sf0(i),i=1,nnt)
         if (iprint.ne.0) then
            write(nsyso,'(/&
              &'' nu*fission (groups 1-'',i3,'')''/&
              &'' ------------------------'')')&
              nnt
            write(nsyso,'(1p,6e12.4)') (snus(i),i=1,nnt)
            write(nsyso,'(/&
              &'' fission (groups 1-'',i3,'')''/&
              &'' ---------------------'')')&
              nnt
            write(nsyso,'(1p,6e12.4)') (sf0(i),i=1,nnt)
         endif
      endif
      k=0
      if (iprint.gt.0) write(nsyso,'(/&
        &'' nonthermal p0 scattering matrix''/&
        &'' -------------------------------'')')
      if (iprint.gt.0) write(nsyso,&
        &'('' initl'',3x,''final'',4x,&
        &''cross section vs final group''/&
        &'' group'',3x,''group'',4x,''+0'',4(10x,''+'',i1)/&
        &'' -----'',3x,''-----'',3x,5(''----------''),&
        &''---------'')') (i,i=1,4)
      do ia=1,nnt
         lone=l1(ia)
         ltwo=l2(ia)
         if (lone.ne.ngnd.or.ltwo.ne.1) then
            k=k+1
            ks=k+2
            scr(k)=ia-lone+1
            k=k+1
            scr(k)=ltwo-lone+1
            if (lone.ne.ltwo) then
               do i=lone,ltwo
                  loc=1+ia-1+ngnd*(i-1)
                  k=k+1
                  scr(k)=xs(loc)
               enddo
               if (iprint.ne.0) then
                  lim=k
                  if ((k-ks+1).gt.ncol) lim=ks+ncol-1
                  write(nsyso,'(1x,i4,4x,i4,2x,1p,5e12.4)')&
                    ia,lone,(scr(j),j=ks,lim)
                  ib=lone
                  do while (lim.ne.k)
                     ks=lim+1
                     lim=k
                     if ((k-ks+1).gt.ncol) lim=ks+ncol-1
                     ib=ib+ncol
                     write(nsyso,'(9x,i4,2x,1p,5e12.4)')&
                       ib,(scr(j),j=ks,lim)
                  enddo
               endif
            endif
         else
            k=k+1
            scr(k)=1
            k=k+1
            scr(k)=0
            if (iprint.ne.0) then
               write(nsyso,'(1x,i4,4x,i4,2x,1p,5e12.4)')&
                 ia,lone,zero
            endif
         endif
      enddo
      if (k.gt.nwscr) call error('xseco',&
        'scratch storage exceeded',' ')
      write(nscr2) k
      write(nscr2) (scr(i),i=1,k)
         if (iprint.ne.0) then
      write(nsyso,'(/&
        &'' thermal cross sections vs temperature '',&
        &''(groups '',i3,''-'',i3,'')''/&
        &1x,5(''----------''),''--'')')&
        nnt+1,ngnd
      endif
   endif

   !--write temperature-dependent part
   nt0=nnt+1
   jtemp=nint(tempr(itemp))
   if (itemp.eq.1) then
      write(nsyso,'(/'' temperatures''/'' ------------'')')
      write(nsyso,'(1p,6e12.4)') (tempr(i),i=1,ntemp)
   endif
   write(nscr3)&
     (xtr(i),i=nt0,ngnd),(ab0(i),i=nt0,ngnd)
   if (iprint.ne.0) then
      write(nsyso,'(/&
        &'' transport corrected total at '',i4,'' deg k''/&
        &1x,3(''----------''),''---------'')')&
        jtemp
      write(nsyso,'(1p,6e12.4)') (xtr(i),i=nt0,ngnd)
      write(nsyso,'(/&
        &'' absorption at '',i4,'' deg k''/&
        &'' ------------------------'')')&
        jtemp
      write(nsyso,'(1p,6e12.4)') (ab0(i),i=nt0,ngnd)
   endif
   if (ifiss.gt.1) then
      write(nscr3)&
        (snus(i),i=nt0,ngnd),(sf0(i),i=nt0,ngnd)
      if (iprint.ne.0) then
         write(nsyso,'(/&
           &'' nu*fission at '',i4,'' deg k''/&
           &'' ------------------------'')')&
           jtemp
         write(nsyso,'(1p,6e12.4)') (snus(i),i=nt0,ngnd)
         write(nsyso,'(/&
           &'' fission at '',i4,'' deg k''/&
           &'' ---------------------'')')&
           jtemp
         write(nsyso,'(1p,6e12.4)') (sf0(i),i=nt0,ngnd)
      endif
   endif
   k=0
   if (iprint.gt.0) write(nsyso,'(/&
     &'' thermal p0 scattering matrix at '',i4,'' deg k''/&
     &'' ------------------------------------------'')')&
     jtemp
   if (iprint.gt.0) write(nsyso,&
     '('' initl'',3x,''final'',4x,''cross section vs final group''/&
     &'' group'',3x,''group'',4x,''+0'',4(10x,''+'',i1)/&
     &'' -----'',3x,''-----'',3x,5(''----------''),''---------'')')&
     (i,i=1,4)
   do ia=nt0,ngnd
      lone=l1(ia)
      ltwo=l2(ia)
      if (lone.ne.ngnd.or.ltwo.ne.1) then
         k=k+1
         ks=k+2
         scr(k)=ia-lone+1
         k=k+1
         scr(k)=ltwo-lone+1
         if (lone.ne.ltwo) then
            do i=lone,ltwo
               loc=1+ia-1+ngnd*(i-1)
               k=k+1
               scr(k)=xs(loc)
            enddo
            if (iprint.ne.0) then
               lim=k
               if ((k-ks+1).gt.ncol) lim=ks+ncol-1
               write(nsyso,'(1x,i4,4x,i4,2x,1p,5e12.4)')&
                 ia,lone,(scr(j),j=ks,lim)
               ib=lone
               do while (lim.ne.k)
                  ks=lim+1
                  lim=k
                  if ((k-ks+1).gt.ncol) lim=ks+ncol-1
                  ib=ib+ncol
                  write(nsyso,'(9x,i4,2x,1p,5e12.4)')&
                    ib,(scr(j),j=ks,lim)
               enddo
            endif
         endif
      else
         k=k+1
         scr(k)=1
         k=k+1
         scr(k)=0
         if (iprint.ne.0) then
            write(nsyso,'(1x,i4,4x,i4,2x,1p,5e12.4)') ia,lone,zero
         endif
      endif
   enddo
   if (k.gt.nwscr) call error('xseco',&
     'scratch storage exceeded',' ')
   write(nscr3) k
   write(nscr3) (scr(i),i=1,k)

   !--xseco is finished.
   return
   end subroutine xseco

   subroutine p1scat
   !--------------------------------------------------------------------
   ! Process P1 scattering matrix.
   !--------------------------------------------------------------------
   use util ! provides openz,repoz,error
   use endf ! provides endf routines and variables,
   use mainio ! provides nsysi,nsyso,nsyse
   ! internals
   integer::nther,ngndsq,matd,nb,nw,il,mtelas,itd,i
   integer::jtemp,nth1,nth,ip1,nl,nz,l,ng2,ig2lo
   integer::ig,jg,ig2,jg2,jz,loca,loc,lt,l1t,l1tc
   integer::l2t,l2tc,nthr,loc1,loc2,nscr
   character(60)::strng
   integer,dimension(:),allocatable::l1,l2
   integer,dimension(:),allocatable::l1e,l2e
   real(kr),dimension(:),allocatable::elas
   real(kr),dimension(:),allocatable::sloc
   integer,parameter::lz=6

   !--allocate storage.
   if (ixs.eq.0) return
   if (ip1opt.eq.1) return
   write(nsyso,'(/'' ***p1 scattering matrices***'')')
   nther=ngnd-nfg-nrg
   allocate(l1(ngnd))
   allocate(l2(ngnd))
   allocate(l1e(ngnd))
   allocate(l2e(ngnd))
   ngndsq=ngnd*ngnd
   allocate(elas(ngndsq))
   allocate(sloc(ngndsq))
   call openz(-nscr4,1)
   nscr=nscr0
   call openz(-nscr,1)

   !--read gendf data
   matd=mat
   call repoz(ngendf)
   call tpidio(ngendf,0,0,scr,nb,nw)
   il=2
   mtelas=0
   itd=0
   do i=1,ngndsq
      sloc(i)=0
   enddo
   do i=1,ngnd
      l1(i)=ngnd
      l2(i)=1
   enddo
   jtemp=-1
   nth1=0
   nth=nther
   ip1=0

   !--search for required reactions.
  125 continue
   call contio(ngendf,0,0,scr,nb,nw)
   if (math.eq.0) go to 200
   if (math.lt.0) go to 300
   if (math.ne.matd) go to 300
   if (mfh.eq.0.or.mth.eq.0) go to 125
   if (mfh.eq.1.and.jtemp.eq.ntemp) go to 300
   nl=l1h
   if (nl.lt.il.and.mfh.ne.1) go to 155
   nz=l2h

   !--loop over groups for this reaction
  140 continue
   l=1
   call listio(ngendf,0,0,scr(l),nb,nw)
   if (mfh.eq.1) jtemp=jtemp+1
   if (mfh.eq.3) go to 160
   if (mfh.ne.6) go to 155
   if (mth.eq.2.and.mtelas.eq.0) then
      mtelas=1
      do i=1,ngnd
         l1e(i)=ngnd
         l2e(i)=1
      enddo
      do i=1,ngndsq
         elas(i)=0
      enddo
   endif
   ng2=l1h
   ig2lo=l2h
   ig=n2h
   do while (nb.ne.0)
      l=l+nw
      call moreio(ngendf,0,0,scr(l),nb,nw)
      if (l+nw-1.gt.nwscr) call error('p1scat',&
        'storage exceeded',' ')
   enddo
   go to 160
  155 continue
   call tosend(ngendf,0,0,scr)
   go to 125
  160 continue
   if (ig.eq.0.or.ig.gt.ngnd) go to 225
      jg=ngnd-ig+1
   l=1

   !--store non-temperature dependent scattering matrix
   if (jtemp.gt.0) go to 210
   if (mfh.eq.3) go to 230
   if (mth.ge.18.and.mth.le.21) go to 155
   if (mth.eq.38) go to 155
   if (mth.eq.2) go to 155
   if (mth.eq.mti.or.mth.eq.mtc) go to 180
   if (mth.ge.221.and.mth.le.250) go to 155
   do i=2,ng2
      ig2=ig2lo+i-2
      if (ig2.ne.0.and.ig2.le.ngnd) then
         jg2=ngnd-ig2+1
         jz=nz
         if(isg.gt.0 .and. isg.lt.nz) jz=isg
         loca=l+lz+(il-1)+nl*nz*(i-1)+(jz-1)*nl
         loc=1+jg-1+ngnd*(jg2-1)
         sloc(loc)=sloc(loc)+scr(loca)
      endif
   enddo
   ip1=1
   lt=l1(jg)
   l1t=ngnd-ig2lo-ng2+3
   l1tc=l1t
   if (l1tc.lt.lt.and.l1tc.ne.0.and.l1tc.le.ngnd) l1(jg)=l1tc
   lt=l2(jg)
   l2t=l1t+ng2-2
   l2tc=l2t
   if (l2tc.gt.lt.and.l2tc.ne.l1t.and.l2tc.gt.0.and.l2tc.le.ngnd)&
        l2(jg)=l2tc

   !--continue loop over groups
   if (ig.lt.ngnd) go to 140

   !--continue loop over reactions
   go to 125
  180 continue
   if (ig.gt.nth1.and.ig.le.nth) nth1=ig
   if (ig.lt.nth) go to 140
   go to 155

   !--write non-temperature dependent matrix on scratch tape.
  200 continue
   if (jtemp.gt.0) go to 250
   call repoz(-nscr)
   write(nscr) (l1(i),i=1,ngnd)
   write(nscr) (l2(i),i=1,ngnd)
   write(nscr) (sloc(i),i=1,ngndsq)
   call repoz(-nscr)
   call repoz(ngendf)
   call tpidio(ngendf,0,0,scr,nb,nw)
   call findf(matd,1,451,ngendf)
   nth=nth1
   go to 125

   !--add temperature-dependent reaction to matrix.
  210 continue
   if (mth.eq.2.or.mth.eq.mti.or.mth.eq.mtc) go to 215
   go to 155
  215 continue
   if (mth.eq.2.and.ig.le.nth) go to 140
   if (mth.eq.mti.and.ig.gt.nth) go to 155
   if (mth.eq.mtc.and.ig.gt.nth.and.mtc.ne.0) go to 155
   itd=1
   do i=2,ng2
      ig2=ig2lo+i-2
      if (ig2.ne.0.and.ig2.le.ngnd) then
         jg2=ngnd-ig2+1
         jz=nz
         if(isg.gt.0 .and. isg.lt.nz) jz=isg
         loca=l+lz+(il-1)+nl*nz*(i-1)+(jz-1)*nl
         loc=1+jg-1+ngnd*(jg2-1)
         sloc(loc)=sloc(loc)+scr(loca)
         if (mth.eq.2) then
            loc=1+jg-1+ngnd*(jg2-1)
            elas(loc)=elas(loc)+scr(loca)
         endif
      endif
   enddo
   ip1=1
   lt=l1(jg)
   l1t=ngnd-ig2lo-ng2+3
   l1tc=l1t
   if (l1tc.lt.lt.and.l1tc.ne.0) l1(jg)=l1tc
   if (mth.eq.2) l1e(jg)=l1tc
   lt=l2(jg)
   l2t=l1t+ng2-2
   l2tc=l2t
   if (l2tc.gt.lt.and.l2tc.ne.l1t.and.l2tc.gt.0.and.l2tc.le.ngnd)&
     l2(jg)=l2tc
   if (mth.eq.2.and.l2tc.ne.l1t) l2e(jg)=l2tc
  225 continue
   if (ig.lt.ngnd) go to 140
   go to 125

   !--coarse group flux
  230 continue
   if (mth.ne.1) go to 155
   if (ig.lt.ngnd) go to 140
   go to 125

   !--desired data is loaded
  250 continue
   if (ip1.eq.0) then
      write(strng,'(''no p1 matrices found for mat '',i4)') matd
      call error('p1scat',strng,' ')
   endif
   if (mti.gt.0.and.itd.eq.0) then
      write(strng,'(''no temperature-dependent reactions for mat '',&
        &i4)') matd
      call error('p1scat',strng,' ')
   endif
   if (mtelas.le.0) then
      do i=1,ngnd
         if (l1e(i).lt.l1(i)) l1(i)=l1e(i)
         if (l2e(i).gt.l2(i)) l2(i)=l2e(i)
      enddo
      do i=1,ngndsq
         sloc(i)=sloc(i)+elas(i)
      enddo
   endif

   !--suppress upscattering from thermal into resonance groups
   nthr=nfg+nrg+1
   do jg=nthr,ngnd
      jg2=l1(jg)
      do while (jg2.lt.nthr)
         loc1=1+jg-1+ngnd*(jg2-1)
         loc2=1+jg-1+ngnd*(nthr-1)
         sloc(loc2)=sloc(loc2)+sloc(loc1)
         sloc(loc1)=0
         jg2=jg2+1
         l1(jg)=jg2
      enddo
   enddo

   !--print to output and to scratch tape
   call p1sout(l1,l2,sloc,jtemp)
   call repoz(-nscr)
   read(nscr) (l1(i),i=1,ngnd)
   read(nscr) (l2(i),i=1,ngnd)
   read(nscr) (sloc(i),i=1,ngndsq)
   mtelas=0
   go to 125

   !--this material is finished.
  300 continue

   !--p1scat is finished.
   return
   end subroutine p1scat

   subroutine p1sout(l1,l2,sloc,it)
   !--------------------------------------------------------------------
   ! Write P1 scattering matrices.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso,nsyse
   ! externals
   integer::l1(*),l2(*),it
   real(kr)::sloc(*)
   ! internals
   integer::ncol1,ktemp,i,ig,lone,ltwo,j,loc,nb,i1,lim
   integer,parameter::ncol=5

   !--retrieve matrix pointers
   ncol1=ncol-1

   !--write p1 scattering matrix on nscr4
   ktemp=nint(tempr(it))
   if (iprint.gt.1) write(nsyso,'(/&
     &'' p1 scattering matrix at '',i4,'' deg k''/&
     &1x,5(''----------''),''--------''/&
     &'' ig l1 l2    l1+0'',4(10x,''+'',i1)/&
     &1x,7(''----------''),''----'')')&
     ktemp,(i,i=1,ncol1)
   do ig=1,ngnd
      lone=l1(ig)
      ltwo=l2(ig)
      j=0
      do i=lone,ltwo
         loc=1+ig-1+ngnd*(i-1)
         j=j+1
         scr(j)=sloc(loc)
      enddo
      nb=j
      write(nscr4) ig,lone,ltwo
      write(nscr4) (scr(i),i=1,nb)
      if (iprint.ne.0) then
         i1=1
         lim=nb
         if (lim.gt.ncol) lim=ncol
         write(nsyso,'(3i3,5x,1p,5e12.4)')&
           ig,lone,ltwo,(scr(i),i=i1,lim)
         do while (lim.ne.nb)
            i1=lim+1
            lim=nb
            if ((lim-i1+1).gt.ncol) lim=i1+ncol-1
            write(nsyso,'(14x,1p,5e12.4)') (scr(i),i=i1,lim)
         enddo
      endif
   enddo
   return
   end subroutine p1sout

   subroutine wimout(jcc)
   !--------------------------------------------------------------------
   ! Write the WIMS-format data on nout.
   !--------------------------------------------------------------------
   use util ! provides timer,repoz
   use mainio ! provides nsysi,nsyso,nsyse
   use physics ! provides amassn
   ! externals
   integer::jcc
   ! internals
   integer::ngr0,ngr1,nnt,ifis,nrestb,jfis,nw,ident,npos
   integer::jcc2,i,ndat,nt0,it,j,ntnp,jres,jsigz,irg,nump1
   integer::ig,lone,ltwo,nb
   real(kr)::time,awt,xid,rid
   integer,parameter::izero=0

   !--print timer message
   call timer(time)
   write(nsyso,'(/'' ***output***'',56x,f8.1,''s'')') time

   !--initialize.
   ngr0=nfg+1
   ngr1=nfg+nrg
   nnt=ngr1

   !--retrieve fission flag
   ifis=ifiss
   if (ifis.eq.0.or.ifis.eq.4) nrestb=0
   if (ifis.gt.0.and.ifis.lt.4) nrestb=1
   jfis=ifis
   if (ifprod.eq.1) jfis=-1
   if (ifprod.eq.2) jfis=-2
   nw=4*(ngr1-ngr0+1)+2*nnt
   if (iverw.eq.5) nw=nw+nfg

   !--write identifiers
   if (iverw.ne.4) write(nout,'(i15)') nfid
   ident=nfid
   npos=1
   if (iverw.ne.4) write(nout,'(2i15)') ident,npos

   !--write burnup data.
   if (iburn.ge.0) then
      if (iburn.eq.0.and.iverw.ne.4) write(nout,'(i15)') izero
      if (jcc.ge.8) ifisp(4)=jfis
      write(nout,'(2i15)') 999999999,3
      write(nout,'(i15)') jcc
      jcc2=jcc/2
      write(nout,'(3(1p,e15.8,i6))') (yield(i),ifisp(i),i=1,jcc2)
      write(nout,'(2i15)') 999999999,4
   endif

   !--write material identifier data
   awt=awr*amassn
   write(nout,'(i6,1p,e15.8,5i6)')&
     ident,awt,iznum,ifis,ntemp,nrestb,isof

   !--write temperature-independent data
   !--potential, scatter, transport, absorption, (unused), lambda
   call repoz(-nscr2)
   read(nscr2) (scr(i),i=1,nw)
   write(nout,'(1p,5e15.8)') (scr(i),i=1,nw)
   if (ifis.gt.1) then
      ! nu*fission, fission
      nw=2*nnt
      read(nscr2) (scr(i),i=1,nw)
      write(nout,'(1p,5e15.8)') (scr(i),i=1,nw)
   endif
   ! nonthermal scattering matrix
   read(nscr2) ndat
   write(nout,'(i15)') ndat
   if (ndat.gt.0) then
      read(nscr2) (scr(i),i=1,ndat)
      write(nout,'(1p,5e15.8)') (scr(i),i=1,ndat)
   endif

   !--write temperature-dependent data
   call repoz(-nscr3)
   write(nout,'(1p,5e15.8)') (tempr(i),i=1,ntemp)
   nt0=nnt+1
   do it=1,ntemp
      nw=2*(ngnd-nt0+1)
      ! transport, absorption
      read(nscr3) (scr(j),j=1,nw)
      write(nout,'(1p,5e15.8)') (scr(j),j=1,nw)
      if (ifis.gt.1) then
         ! nu*fission, fission
         read(nscr3) (scr(j),j=1,nw)
         write(nout,'(1p,5e15.8)') (scr(j),j=1,nw)
      endif
      ! thermal scattering matrix
      read(nscr3) ndat
      write(nout,'(i15)') ndat
      if (ndat.ne.0) then
         read(nscr3) (scr(j),j=1,ndat)
         write(nout,'(1p,5e15.8)') (scr(j),j=1,ndat)
      endif
   enddo
   if (iverw.eq.4) write(nout,'(''      999999999'')')

   !--write resonance data
   if (ires.ne.0) then
      call repoz(-nscr1)
      read(nscr1) ntnp
      read(nscr1) xid,jres,jsigz
      rid=rdfid
      do irg=1,nrg
         if (iverw.eq.4) write(nout,'(i15)') ntnp
         write(nout,'(1p,e15.8,2i6)') rid,jres,jsigz
         ! temperatures, sig pots, resonance integrals
         nw=jres+jsigz+ntnp
         read(nscr1) (scr(i+1),i=1,nw)
         write(nout,'(1p,5e15.8)') (scr(i+1),i=1,nw)
         if (ifis.eq.3) then
            ! temperatures, sig pots, fission resonance integrals
            if (iverw.eq.4) write(nout,'(i15)') ntnp
            write(nout,'(1p,e15.8,2i6)') rid,jres,jsigz
            read(nscr1) (scr(i+1),i=1,nw)
            write(nout,'(1p,5e15.8)') (scr(i+1),i=1,nw)
         endif
         if (iverw.eq.4) write(nout,'(''      999999999'')')
      enddo
      if (iverw.ne.4) then
         do irg=1,nrg
            write(nout,'(1p,e15.8,2i6)') rid,jres,jsigz
            nw=jres+jsigz+ntnp
            read(nscr1) (scr(i+1),i=1,nw)
            write(nout,'(1p,5e15.8)') (scr(i+1),i=1,nw)
         enddo
      endif
   endif

   !--write fission spectrum
   if (isof.ne.0) then
      write(nout,'(1p,5e15.8)') (uff(j),j=1,nfiss)
   endif

   !--write p1 data
   if (ip1opt.ne.1) then
      call repoz(-nscr4)
      do it=1,ntemp
         nump1=0
         do ig=1,ngnd
            read(nscr4) nw,lone,ltwo
            nb=ltwo-lone+1
            scr(1+nump1)=nw-lone+1
            scr(1+nump1+1)=nb
            nump1=nump1+2
            read(nscr4) (scr(nump1+j),j=1,nb)
            nump1=nump1+nb
         enddo
         write(nout,'(i15)') nump1
         write(nout,'(1p,5e15.8)') (scr(j),j=1,nump1)
      enddo
   endif

   !--wimout is finished.
   return
   end subroutine wimout

end module wimsm


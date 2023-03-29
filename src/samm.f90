module samm
   !-------------------------------------------------------------------
   ! Routines for the SAMMY method of calculating resonance cross
   ! sections, derivatives, and angular distributions.  These
   ! routines are used in RECONR and ERRORR.  Converted to NJOY
   ! style and Fortran-90 format from coding provided by Nancy
   ! Larson (ORNL).
   !-------------------------------------------------------------------
   use locale
   implicit none
   private

   !--Public routines
   public cssammy,s2sammy,ppsammy,rdsammy,desammy

   !--Public variables
   logical,public::Want_Partial_Derivs
   logical,public::Want_Angular_Dist
   logical,public::Want_SAMRML_RM
   logical,public::Want_SAMRML_BW

   !--Global variables

   !--flags to control options
   logical::Want_Partial_U=.false.

   !--variables
   integer::ngroup,mchan,nres,npp,ntriag,npar,lllmax,kkxlmn
   integer::nrext,krext,nparr,nerm
   integer,dimension(12)::ngroupm,nppm,nresm
   integer,parameter::lllmaxx=10
   real(kr)::awr,su

   integer,dimension(:,:),allocatable::nchan
   real(kr),dimension(:,:),allocatable::ema,emb,spina,spinb,qqq
   integer,dimension(:,:),allocatable::kza,kzb,ishift,lpent,mt
   real(kr),dimension(:,:),allocatable::pa,pb,sspin,parity
   integer,dimension(:,:),allocatable::nresg
   integer,dimension(:,:,:),allocatable::ipp,lspin
   real(kr),dimension(:,:,:),allocatable::chspin,bound,rdeff,rdtru
   real(kr),dimension(:,:),allocatable::eres,gamgam
   real(kr),dimension(:,:,:),allocatable::gamma
   real(kr),dimension(:,:,:,:),allocatable::parext
   integer,dimension(:,:),allocatable::nent,next
   real(kr),dimension(:,:),allocatable::goj
   real(kr),dimension(:,:,:),allocatable::zke,zkfe,zkte,zeta,echan
   real(kr),dimension(:,:,:),allocatable::betapr,gbetpr
   real(kr),dimension(:,:,:),allocatable::beta
   real(kr),dimension(:,:),allocatable::uuuu,duuu
   integer,dimension(:,:),allocatable::iduu

   real(kr),dimension(:,:),allocatable::alj
   real(kr),dimension(:,:),allocatable::ddddd
   real(kr),dimension(:,:),allocatable::crss,sigmas
   real(kr),dimension(:,:),allocatable::alphar,alphai,difen,xden
   real(kr),dimension(:,:),allocatable::sinsqr,sin2ph,dphi,dpdr,dsdr
   real(kr),dimension(:,:),allocatable::sinphi,cosphi
   real(kr),dimension(:,:),allocatable::rootp,elinvr,elinvi,psmall
   real(kr),dimension(:,:,:),allocatable::xqr,xqi
   real(kr),dimension(:,:),allocatable::xxxxr,xxxxi
   real(kr),dimension(:,:,:),allocatable::qr,qi
   real(kr),dimension(:,:,:),allocatable::rmat,ymat,yinv,cscs
   real(kr),dimension(:),allocatable::xlmn
   integer,dimension(:),allocatable::kxlmn
   real(kr),dimension(:,:),allocatable::upr,upi
   real(kr),dimension(:,:,:),allocatable::br,bi,pr,pii
   real(kr),dimension(:,:,:),allocatable::deriv,dsigma
   real(kr),dimension(:,:,:,:,:),allocatable::derivx
   real(kr),dimension(:,:,:,:,:,:),allocatable::crssx
   real(kr),dimension(:,:,:),allocatable::Coef_Leg
   real(kr),dimension(:,:,:,:),allocatable::D_Coef_Leg
   real(kr),dimension(:,:,:),allocatable::tr,ti
   real(kr),dimension(:,:,:,:),allocatable::tx

   real(kr),dimension(:),allocatable::par

contains

   subroutine cssammy(e,sigp,siga,sigd,mmtres,nmtres,ncoef,nresp,ier)
   !-------------------------------------------------------------------
   ! Calculates RM or RML cross sections at energy e
   ! for one section using the SAMMY method.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   ! externals
   integer::nmtres,ncoef,nresp,ier
   real(kr)::e,sigp(nmtres+2),siga(ncoef,nmtres),sigd(nresp,nmtres)
   integer::mmtres(10)
   ! internals
   integer::i,l
   real(kr),parameter::cut=1.e-3_kr
   real(kr),parameter::zero=0

   su=e

   !--generate energy-dependent pieces of cross sections and derivs
   call abpart(ier)

   !--form the cross sections, angular distributions, and derivatives
   call crosss(ier)

   !--return cross sections
   !  1 = total
   !  2 = elastic
   !  3 = fission
   !  4 = capture
   !  5,6,etc = other channels
   sigp(1)=sigmas(1,ier)+sigmas(2,ier)
   sigp(2)=sigmas(1,ier)
   sigp(3)=0
   sigp(4)=sigmas(2,ier)
   if (nmtres.gt.2) then
      do i=3,nmtres
         sigp(4)=sigp(4)-sigmas(i,ier)
         if (mmtres(i).eq.18) then
            sigp(3)=sigmas(i,ier)
         else
            sigp(2+i)=sigmas(i,ier)
         endif
      enddo
   endif

   !--return angular coefficents
   ! 1 = elastic
   ! 2 = nonelastic (isotropic)
   ! 3,4,etc = other channels
   if (Want_angular_Dist) then
      if (lllmax.gt.1) then
         do l=1,lllmax
            do i=1,nppm(ier)
               if (abs(Coef_Leg(2,i,ier)).gt.1.0e-12_kr) then
                  if (Coef_Leg(1,i,ier).ne.zero)&
                    siga(l,i)=Coef_Leg(l,i,ier)/&
                      Coef_Leg(1,i,ier)/(2*l-1)
                  if (abs(siga(l,i)).lt.1.0e-12_kr) siga(l,i)=0
               else
                  siga(l,i)=0
               endif
            enddo
         enddo
      else
         do i=1,npp
            siga(1,i)=0
         enddo
      endif
   endif

   !--return sensitivities
   ! 1 = elastic
   ! 2 = capture
   ! 3,4,etc = other channels
   if (Want_Partial_Derivs) then
      do i=1,npar
         sigd(i,1)=dsigma(1,i,ier)
         sigd(i,2)=dsigma(2,i,ier)
         if (nmtres.gt.2) then
            do l=3,nmtres
               sigd(i,l)=dsigma(l,i,ier)
               sigd(i,2)=sigd(i,2)-dsigma(l,i,ier)
            enddo
         endif
      enddo
   endif

   return
   end subroutine cssammy

   subroutine s2sammy(nin,res,maxres,mmtres,nmtres)
   !-------------------------------------------------------------------
   ! Scan File 2 to get sizes for SAMMY calculation.
   ! Allocate the storage for the SAMMT arrays.
   !-------------------------------------------------------------------
   use endf   ! provides endf routines and variables
   use util   ! provides error
   use mainio ! provides nsyso
   ! externals
   integer::nin,maxres,nmtres
   integer::mmtres(10)
   real(kr)::res(maxres)
   ! internals
   integer::jnow,nb,nw,lfw,nis,i,lru,lrf,nro,naps,mode,lrx,kf,ki
   integer::nls,ner,j,ng,is,ll,iis,kchan,kpp,kres,kpar,jj,nrs,igroup
   integer::nch,ich,iso,ier,ipp,kk,kki,kkf
   real(kr)::spin,parity,parl,capj,capjmx,gamf,gamf2,el,e1,e2,gamx
   real(kr),dimension(2)::s
   real(kr),parameter::zero=0
   real(kr),parameter::half=0.5e0_kr

   allocate(nchan(100,12))
   nmtres=0
   nrext=0

   !--read head record.
   jnow=1
   call contio(nin,0,0,res(jnow),nb,nw)
   if (mfh.ne.2) then
      write(nsyso,'(/'' mat has no file 2.'')')
      call tosend(nin,0,0,res)
      return
   endif
   awr=c2h
   nis=n1h
   nres=0
   ngroup=0
   npp=0
   mchan=0
   npar=0
   kres=0
   kpp=0
   kchan=0

   !--check for multiple isotopes
   if (nis.gt.1) then
       call mess('s2sammy',&
       'multiple isotopes do not work with sammy ',&
       ' method--reverting to normal RM processing')
      nmtres=0
      go to 10
   endif

   !--print out of resonance range information
   write(nsyso,'('' resonance range information''/&
                 &'' ----------------------------------''/&
                 &'' ier     energy-range     lru lrf  method'')')

   !--do loop over all isotopes.
   do iso=1,nis
      call contio(nin,0,0,res(jnow),nb,nw)
      jnow=jnow+nw
      if (jnow.gt.maxres) call error('s2sammy',&
        'res storage exceeded',' ')
      lfw=l2h
      ner=n1h
      nerm=ner

      !--do loop over all energy ranges.
      do ier=1,ner
         call contio(nin,0,0,res(jnow),nb,nw)
         e1=c1h
         e2=c2h
         lru=l1h
         lrf=l2h

         !--if adler-adler format, set nmtres=0 and return
         !--will use historical njoy coding
         if (lrf.eq.4) then
            call mess('s2sammy',&
            'adler-adler format does not work with sammy',&
            ' method, reverting to normal a-a processing')
            nmtres=0
            return
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

         !--check for energy-dependent scattering radius
         if (nro.ne.0) then
            if (lru.eq.1 .and. lrf.ne.2) call error('s2sammy',&
              &'energy-dep scattering length',&
              &'only works for mlbw')
            call tab1io(nin,0,0,res(jnow),nb,nw)
            jnow=jnow+nw
            if (jnow.gt.maxres) then
               call error('s2sammy','res storage exceeded',' ')
            endif
         endif

         !--process this part according to its formalism
         if (mode.le.2) then  ! scanning SLBW and MLBW formats

            !--read some parameters
            call contio(nin,0,0,res(jnow),nb,nw)
            jnow=jnow+nw
            if (jnow.gt.maxres) call error('s2sammy','res storage exceeded',' ')
            spin=c1h
            nls=n1h

            !--figure out the number of spin groups
            parity=1
            if (spin.lt.zero) parity=-1
            s(1)=abs(spin)-half
            s(2)=abs(spin)+half
            is=2
            if (spin.eq.zero) then
               s(1)=half
               is=1
            endif
            ng=0
            parl=-parity
            do ll=1,nls
               el=ll-1
               parl=-parl
               do iis=1,is
                  capj=abs(abs(s(iis))-el)
                  capjmx=abs(s(iis))+el
                  do
                     ng=ng+1
                     capj=capj+1
                     if (capj.gt.capjmx) exit
                  enddo
               enddo
            enddo
            ngroupm(ier)=ng
            if (ng.gt.ngroup) ngroup=ng

            !--figure out the number of resonances and
            !--the maximum number of channels
            kchan=1
            kres=0
            kpp=2
            kpar=0
            kki=0
            kkf=0
            do i=1,nls
               call listio(nin,0,0,res(jnow),nb,nw)
               jj=jnow+nw
               do while (nb.ne.0)
                  call moreio(nin,0,0,res(jj),nb,nw)
                  jj=jj+nw
                  if (jj.gt.maxres) call error('s2sammy',&
                    'res storage exceeded',' ')
               enddo
               lrx=nint(res(jnow+3))
               nrs=nint(res(jnow+5))
               kres=kres+nrs
               jj=jnow
               kf=0
               ki=0
               do j=1,nrs
                  jj=jj+6
                  gamf=res(jj+5)
                  gamx=0
                  if (lrx.gt.0) gamx=res(jj+2)-res(jj+3)-res(jj+4)-res(jj+5)
                  if (gamx.lt.1.e7_kr) gamx=0
                  if (gamf.ne.zero) kf=1
                  if (gamx.ne.zero) ki=1
               enddo
               if (kf.eq.1) kkf=1
               if (ki.eq.1) kki=1
               kchan=1+kf+ki
               if (kchan.gt.1) then
                  kpp=3
               endif
               kpar=kpar+nres*(kpp+2)
            enddo
            if (kpar.gt.npar) npar=kpar
            if (kchan.gt.mchan) mchan=kchan
            nresm(ier)=kres
            if (kres.gt.nres) nres=kres
            nppm(ier)=kpp
            if (kpp.gt.npp) npp=kpp
            mmtres(1)=2
            mmtres(2)=102
            nmtres=2
            if (kkf.gt.0.and.kki.eq.0) then
               mmtres(3)=18
               nmtres=3
            endif
            if (kkf.eq.0.and.kki.gt.0) then
               mmtres(3)=51
               nmtres=3
            endif
            if (kkf.gt.0.and.kki.gt.0) then
               mmtres(3)=18
               mmtres(4)=51
               nmtres=4
            endif

            !--turn off samrml processing of bw formats if desired
            if (.not.Want_SAMRML_BW) nmtres=0

         else if (mode.eq.3) then  ! scanning RM format

            !--read some parameters
            call contio(nin,0,0,res(jnow),nb,nw)
            jnow=jnow+nw
            if (jnow.gt.maxres) call error('s2sammy','res storage exceeded',' ')
            spin=c1h
            nls=n1h

            !--figure out the number of spin groups
            parity=1
            if (spin.lt.zero) parity=-1
            s(1)=abs(spin)-half
            s(2)=abs(spin)+half
            is=2
            if (spin.eq.zero) then
               s(1)=half
               is=1
            endif
            ng=0
            parl=-parity
            do ll=1,nls
               el=ll-1
               parl=-parl
               do iis=1,is
                  capj=abs(abs(s(iis))-el)
                  capjmx=abs(s(iis))+el
                  do
                     ng=ng+1
                     capj=capj+1
                     if (capj.gt.capjmx) exit
                  enddo
               enddo
            enddo
            ngroupm(ier)=ng
            if (ng.gt.ngroup) ngroup=ng

            !--figure out the number of resonances and
            !--the maximum number of channels
            kchan=1
            kres=0
            kpp=2
            kpar=0
            do i=1,nls
               call listio(nin,0,0,res(jnow),nb,nw)
               jj=jnow+nw
               do while (nb.ne.0)
                  call moreio(nin,0,0,res(jj),nb,nw)
                  jj=jj+nw
                  if (jj.gt.maxres) call error('s2sammy',&
                    'res storage exceeded',' ')
               enddo
               nrs=nint(res(jnow+5))
               kres=kres+nrs
               jj=jnow
               kk=1
               do j=1,nrs
                  jj=jj+6
                  gamf=res(jj+4)
                  gamf2=res(jj+5)
                  if (gamf.ne.zero.and.kk.lt.2) kk=2
                  if (gamf2.ne.zero.and.kk.lt.3) kk=3
               enddo
               if (kk.gt.1) then
                  kchan=kk
                  kpp=3
               endif
               kpar=kpar+nres*(kpp+2)
            enddo
            if (kpar.gt.npar) npar=kpar
            if (kchan.gt.mchan) mchan=kchan
            nresm(ier)=kres
            if (kres.gt.nres) nres=kres
            nppm(ier)=kpp
            if (kpp.gt.npp) npp=kpp
            mmtres(1)=2
            mmtres(2)=102
            nmtres=2
            if (npp.gt.2) then
               mmtres(3)=18
               nmtres=3
            endif

            !--turn off samrml processing of rm if desired
            if (.not.Want_SAMRML_RM) nmtres=0

         else if (mode.eq.7) then  !  scanning new rm limited format

            !--check for multiple isotopes
            if (nis.gt.1) call error('s2sammy',&
              'multiple isotopes do not work with sammy method',' ')

            !--check for energy-dependent scattering radius
            if (nro.ne.0) then
               if (mode.ne.2) call error('s2sammy',&
                 &'energy-dep scattering length',&
                 &'only works for mlbw')
               call tab1io(nin,0,0,res(jnow),nb,nw)
               jnow=jnow+nw
               if (jnow.gt.maxres) then
                  call error('s2sammy','res storage exceeded',' ')
               endif
            endif

            !--read some parameters
            call contio(nin,0,0,res(jnow),nb,nw)
            jnow=jnow+nw
            if (jnow.gt.maxres) call error('s2sammy','res storage exceeded',' ')
            nls=n1h

            !--nls is really the number of spin groups for mode=7
            ngroupm(ier)=nls
            if (nls.gt.ngroup) ngroup=nls

            !--read the list of particle-pair descriptions
            call listio(nin,0,0,res(jnow),nb,nw)
            kpp=l1h
            do ipp=1,kpp
               mmtres(ipp)=nint(res(jnow+15+12*(ipp-1)))
            enddo
            nppm(ier)=kpp
            if (kpp.gt.npp) npp=kpp
            mmtres(1)=2
            mmtres(2)=102
            nmtres=npp
            jnow=jnow+nw

            !--loop over the spin groups
            kchan=0
            kres=0
            kpar=0
            do igroup=1,ngroup
               call listio(nin,0,0,res(jnow),nb,nw)
               jnow=jnow+nw
               nch=n2h
               ich=nch-1
               nchan(igroup,ier)=ich
               if (ich.gt.kchan) kchan=ich
               call listio(nin,0,0,res(jnow),nb,nw)
               jj=jnow+nw
               do while (nb.ne.0)
                  call moreio(nin,0,0,res(jj),nb,nw)
                  jj=jj+nw
                  if (jj.gt.maxres) call error('s2sammy',&
                    'res storage exceeded',' ')
               enddo
               nrs=nint(res(jnow+3))
               kpar=kpar+nrs*(ich+2)
               kres=kres+nrs
            enddo
            if (kpar.gt.npar) npar=kpar
            if (kchan.gt.mchan) mchan=kchan
            nresm(ier)=kres
            if (kres.gt.nres) nres=kres
         endif

         !--energy range information
         if (nmtres.eq.0) then
            write(nsyso,'(i3,1x,1p,2e10.3,2i4,''   normal'')') ier,e1,e2,lru,lrf
         else
            write(nsyso,'(i3,1x,1p,2e10.3,2i4,''   sammy'')') ier,e1,e2,lru,lrf
         endif
      enddo
   enddo
  10 continue

   !--allocate arrays for reading the resonance parameters
   if (nmtres.gt.0) then
      lllmax=lllmaxx
      call allo
   endif
   if (nmtres.eq.0) deallocate(nchan)

   return
   end subroutine s2sammy

   subroutine rdsammy(nin,ier,jnow,nro,naps,mode,el,eh,&
     enode,nodes,nodmax,spin,ascat,res,maxres)
   !-------------------------------------------------------------------
   ! Read in resonance parameters for the SAMMY calculation.
   !-------------------------------------------------------------------
   use mainio  ! provides nsyso
   use util    ! provides error
   use endf    ! provides endf routines and variables
   use physics ! get neutron mass
   ! externals
   integer::nin,ier,jnow,nro,naps,mode,nodes,nodmax,maxres
   real(kr)::el,eh,spin,ascat
   real(kr)::enode(nodmax),res(maxres)
   ! internals
   integer::nb,nw,nls,is,ng,ll,iis,i,ires,jj,nrs,llll,ig,igxm,lrx
   integer::j,ndig,kres,igroup,ichp1,ichan,ix,ich,ippx,igamma,l,nx,ires1
   real(kr)::pari,parl,capj,capjmx,c,awri,apl,aptru,apeff,spinjj
   real(kr)::gamf,gamf2,s1,s2,hw,ehalf,ell,x,gamx,qx
   real(kr),dimension(:),allocatable::a
   integer,dimension(20)::igx
   integer,dimension(:),allocatable::idum
   real(kr),dimension(2)::s
   real(kr),parameter::zero=0
   real(kr),parameter::one=1
   real(kr),parameter::three=3
   real(kr),parameter::half=0.5e0_kr
   real(kr),parameter::zz8=0.08e0_kr
   real(kr),parameter::z123=0.123e0_kr
   ! use imf2=1 to output converted mf2
   !integer::imf2=1
   integer::imf2=0


   if (nro.ne.0) call error('rdsammy','energy dependent scattering radius',&
     'not currently allowed for sammy rml or rm methods')

   allocate(idum(nres))
   igamma=0

   !--read some parameters
   call contio(nin,0,0,res(jnow),nb,nw)
   jnow=jnow+nw
   if (jnow.gt.maxres) call error('rdsammy','res storage exceeded',' ')
   nls=n1h
   spin=c1h
   ascat=c2h

   !--read parameters in SLBW or MLBW format and convert to RML arrays
   if (mode.le.2) then

      !--generate the particle-pair quantities
      do i=1,nppm(ier)
         ishift(i,ier)=0
         lpent(i,ier)=0
         kza(i,ier)=0
         kzb(i,ier)=0
         mt(i,ier)=0
         spina(i,ier)=0
         spinb(i,ier)=0
         ema(i,ier)=0
         emb(i,ier)=0
         qqq(i,ier)=0
      enddo
      spinb(2,ier)=spin
      lpent(2,ier)=1
      mt(1,ier)=102
      mt(2,ier)=2
      if (nppm(ier).eq.3) mt(3,ier)=18

      !--fill in the defaults
      call ppdefs(ier)

      !--figure out the number of spin groups
      do i=1,mchan
         do j=1,ngroupm(ier)
            lspin(i,j,ier)=0
            chspin(i,j,ier)=0
         enddo
      enddo
      pari=1
      if (spin.lt.zero) pari=-1
      s(1)=abs(spin)-half
      s(2)=abs(spin)+half
      is=2
      if (spin.eq.zero) then
         s(1)=half
         is=1
      endif
      s1=s(1)
      s2=s(2)
      ng=0
      parl=-pari
      do ll=1,nls
         ell=ll-1
         l=ll-1
         parl=-parl
         do iis=1,is
            capj=abs(abs(s(iis))-ell)
            capjmx=abs(s(iis))+ell
            do
               c=capj*parl
               ng=ng+1
               sspin(ng,ier)=c
               parity(ng,ier)=parl
               ipp(1,ng,ier)=2
               lspin(1,ng,ier)=l
               chspin(1,ng,ier)=s(iis)
               capj=capj+1
               if (capj.gt.capjmx) exit
            enddo
         enddo
      enddo
      ngroupm(ier)=ng
      do i=1,ngroupm(ier)
         nchan(i,ier)=1
      enddo

      !--read in the resonance parameters
      ires=0
      do i=1,nls

         call listio(nin,0,0,res(jnow),nb,nw)
         jj=jnow+nw
         do while (nb.ne.0)
            call moreio(nin,0,0,res(jj),nb,nw)
            jj=jj+nw
            if (jj.gt.maxres) call error('rdsammy',&
              'res storage exceeded',' ')
         enddo
         qx=res(jnow+1)
         lrx=nint(res(jnow+3))
         nrs=nint(res(jnow+5))

         !--analyze scattering length
         awri=res(jnow)
         emb(2,ier)=awri
         aptru=z123*(awri*amassn)**(one/three)+zz8
         apeff=ascat
         if (naps.eq.1) then
            aptru=ascat
         endif

         !--find the group number
         llll=nint(res(jnow+2))
         do ig=1,20
            igx(ig)=0
         enddo
         igxm=0
         do j=1,ngroupm(ier)
            if (llll.eq.lspin(1,j,ier)) then
               igxm=igxm+1
               igx(igxm)=j
            endif
         enddo

         !--load the scattering lengths
         do ig=1,igxm
            if (igx(ig).ne.0) then
               rdeff(1,igx(ig),ier)=apeff
               rdtru(1,igx(ig),ier)=aptru
            endif
         enddo

         !--loop over the resonances
         jj=jnow+6
         do j=1,nrs
            ires=ires+1

            !--add center and half-height energies to nodes
            hw=res(jj+2)/2
            hw=hw+(res(jj+3)+abs(res(jj+4))+abs(res(jj+5)))/2
            ndig=5
            if (res(jj).gt.zero) ndig=2+nint(log10(res(jj)/(hw/10)))
            if (ndig.lt.5) ndig=5
            if (ndig.gt.9) ndig=9
            if (res(jj).gt.el.and.res(jj).lt.eh) then
              nodes=nodes+1
              enode(nodes)=sigfig(res(jj),ndig,0)
            endif
            if ((res(jj)+hw).gt.el.and.(res(jj)+hw).lt.eh) then
               if (nodes.ge.nodmax) call error('rdsammy',&
                 'storage in enode exceeded.',' ')
               ehalf=res(jj)+hw
               nodes=nodes+1
               enode(nodes)=sigfig(ehalf,ndig,0)
            endif
            if ((res(jj)-hw).gt.el.and.res(jj)-hw.lt.eh) then
               if (nodes.ge.nodmax) call error('rdsammy',&
                 'storage in enode exceeded.',' ')
               ehalf=res(jj)-hw
               nodes=nodes+1
               enode(nodes)=sigfig(ehalf,ndig,0)
            endif

            !--find the spin group number for this resonance
            spinjj=res(jj+1)
            call findsp(sspin,lspin,chspin,llll,spinjj,s1,s2,ig,ier)
            idum(ires)=ig

            !--store the parameters in the sammy arrays
            eres(ires,ier)=res(jj)
            gamma(1,ires,ier)=res(jj+3)
            gamgam(ires,ier)=res(jj+4)
            gamf=res(jj+5)
         gamf=gamf*(-1)**ires
            gamx=0
            if (lrx.gt.0) gamx=res(jj+2)-res(jj+3)-res(jj+4)-res(jj+5)
            if (gamx.lt.1.e-7_kr) gamx=0
            if (gamf.ne.zero.and.gamx.eq.zero) then
               gamma(2,ires,ier)=gamf
               if (nchan(ig,ier).lt.2) then
                  nchan(ig,ier)=2
                  ipp(2,ig,ier)=3
               endif
            endif
            if (gamf.eq.zero.and.gamx.ne.0) then
               gamma(2,ires,ier)=gamx
               if (nchan(ig,ier).lt.2) then
                  nchan(ig,ier)=2
                  ipp(3,ig,ier)=3
               endif
            endif
            if (gamf.ne.zero.and.gamx.ne.0) then
               gamma(2,ires,ier)=gamf
               gamma(3,ires,ier)=gamx
               if (nchan(i,ier).lt.3) then
                  nchan(ig,ier)=3
                  ipp(3,ig,ier)=3
               endif
            endif
            jj=jj+6
         enddo
      enddo

      !--rearrange to put resonances in spin-group order
      call rearrange(idum,ier)

   !--read parameters in RM format and convert to RML arrays
   else if (mode.eq.3) then

      !--generate the particle-pair quantities
      do i=1,nppm(ier)
         ishift(i,ier)=0
         lpent(i,ier)=0
         kza(i,ier)=0
         kzb(i,ier)=0
         mt(i,ier)=0
         spina(i,ier)=0
         spinb(i,ier)=0
         ema(i,ier)=0
         emb(i,ier)=0
         qqq(i,ier)=0
      enddo
      spinb(2,ier)=spin
      lpent(2,ier)=1
      mt(1,ier)=102
      mt(2,ier)=2
      if (nppm(ier).eq.3) mt(3,ier)=18

      !--fill in the defaults
      call ppdefs(ier)

      !--figure out the number of spin groups
      do i=1,mchan
         do j=1,ngroupm(ier)
            lspin(i,j,ier)=0
            chspin(i,j,ier)=0
         enddo
      enddo
      pari=1
      if (spin.lt.zero) pari=-1
      s(1)=abs(spin)-half
      s(2)=abs(spin)+half
      is=2
      if (spin.eq.zero) then
         s(1)=half
         is=1
      endif
      s1=s(1)
      s2=s(2)
      ng=0
      parl=-pari
      do ll=1,nls
         ell=ll-1
         l=ll-1
         parl=-parl
         do iis=1,is
            capj=abs(abs(s(iis))-ell)
            capjmx=abs(s(iis))+ell
            do
               c=capj*parl
               ng=ng+1
               sspin(ng,ier)=c
               parity(ng,ier)=parl
               ipp(1,ng,ier)=2
               lspin(1,ng,ier)=l
               chspin(1,ng,ier)=s(iis)
               capj=capj+1
               if (capj.gt.capjmx) exit
            enddo
         enddo
      enddo
      ngroup=ng

      do i=1,ngroupm(ier)
         nchan(i,ier)=1
      enddo

      !--read in the resonance parameters
      ires=0
      do i=1,nls

         call listio(nin,0,0,res(jnow),nb,nw)
         jj=jnow+nw
         do while (nb.ne.0)
            call moreio(nin,0,0,res(jj),nb,nw)
            jj=jj+nw
            if (jj.gt.maxres) call error('rdsammy',&
              'res storage exceeded',' ')
         enddo
         nrs=nint(res(jnow+5))

         !--analyze scattering length
         awri=res(jnow)
         emb(2,ier)=awri
         apl=res(jnow+1)
         if (apl.eq.zero) then
            aptru=z123*(awri*amassn)**(one/three)+zz8
            apeff=ascat
         else
            apeff=apl
            if (naps.eq.1) then
               aptru=apl
            else
               aptru=z123*(awri*amassn)**(one/three)+zz8
            endif
         endif

         !--find the group number
         llll=nint(res(jnow+2))
         do ig=1,20
            igx(ig)=0
         enddo
         igxm=0
         do j=1,ngroup
            if (llll.eq.lspin(1,j,ier)) then
               igxm=igxm+1
               igx(igxm)=j
            endif
         enddo

         !--load the scattering lengths
         do ig=1,igxm
            if (igx(ig).ne.0) then
               rdeff(1,igx(ig),ier)=apeff
               rdtru(1,igx(ig),ier)=aptru
            endif
         enddo

         !--loop over the resonances
         jj=jnow+6
         do j=1,nrs
            ires=ires+1

            !--add center and half-height energies to nodes
            hw=res(jj+2)/2
            hw=hw+(res(jj+3)+abs(res(jj+4))+abs(res(jj+5)))/2
            ndig=5
            if (res(jj).gt.zero) ndig=2+nint(log10(res(jj)/(hw/10)))
            if (ndig.lt.5) ndig=5
            if (ndig.gt.9) ndig=9
            if (res(jj).gt.el.and.res(jj).lt.eh) then
              nodes=nodes+1
              enode(nodes)=sigfig(res(jj),ndig,0)
            endif
            if ((res(jj)+hw).gt.el.and.(res(jj)+hw).lt.eh) then
               if (nodes.ge.nodmax) call error('rdsammy',&
                 'storage in enode exceeded.',' ')
               ehalf=res(jj)+hw
               nodes=nodes+1
               enode(nodes)=sigfig(ehalf,ndig,0)
            endif
            if ((res(jj)-hw).gt.el.and.res(jj)-hw.lt.eh) then
               if (nodes.ge.nodmax) call error('rdsammy',&
                 'storage in enode exceeded.',' ')
               ehalf=res(jj)-hw
               nodes=nodes+1
               enode(nodes)=sigfig(ehalf,ndig,0)
            endif

            !--find the spin group number for this resonance
            spinjj=res(jj+1)
            call findsp(sspin,lspin,chspin,llll,spinjj,s1,s2,ig,ier)
            idum(ires)=ig

            !--store the parameters in the sammy arrays
            eres(ires,ier)=res(jj)
            gamma(1,ires,ier)=res(jj+2)
            gamgam(ires,ier)=res(jj+3)
            gamf=res(jj+4)
            gamf2=res(jj+5)
            if (gamf.ne.zero) then
               gamma(2,ires,ier)=gamf
               if (nchan(ig,ier).lt.2) then
                  nchan(ig,ier)=2
                  ipp(2,ig,ier)=3
               endif
            endif
            if (gamf2.ne.zero) then
               gamma(3,ires,ier)=gamf2
               if (nchan(ig,ier).lt.3) then
                  nchan(ig,ier)=3
                  ipp(3,ig,ier)=3
               endif
            endif
            jj=jj+6
         enddo
      enddo

      !--rearrange to put resonances in spin-group order
      call rearrange(idum,ier)

   !--read in parameters for the RML format and store in sammy arrays.
   else if (mode.eq.7) then

      !--nls is really the number of spin groups for mode=7
      ngroup=nls

      !--read the list of particle-pair descriptions
      !--and store the values in the proper allocated arrays
      call listio(nin,0,0,res(jnow),nb,nw)
      jj=jnow+nw
      do while (nb.ne.0)
         call moreio(nin,0,0,res(jj),nb,nw)
         jj=jj+nw
         if (jj.gt.maxres) call error('rdsammy',&
           'res storage exceeded',' ')
      enddo
      npp=nint(res(jnow+2))
      jj=jnow+6
      do i=1,npp
         ema(i,ier)=res(jj)
         emb(i,ier)=res(jj+1)
         kza(i,ier)=nint(res(jj+2))
         kzb(i,ier)=nint(res(jj+3))
         spina(i,ier)=res(jj+4)
         spinb(i,ier)=res(jj+5)
         qqq(i,ier)=res(jj+6)
         lpent(i,ier)=nint(res(jj+7))
         ishift(i,ier)=nint(res(jj+8))
         mt(i,ier)=nint(res(jj+9))
         if (mt(i,ier).eq.2) spin=spinb(i,ier)
         pa(i,ier)=res(jj+10)
         pb(i,ier)=res(jj+11)
         jj=jj+12
      enddo

      !--fill in the defaults
      call ppdefs(ier)

      !--loop over the spin groups
      kres=0
      do igroup=1,ngroup

         !--read the spin group information
         call listio(nin,0,0,res(jnow),nb,nw)
         sspin(igroup,ier)=c1h
         parity(igroup,ier)=c2h
         ichp1=n2h
         ichan=ichp1-1
         ix=0
         i=0
         jj=jnow+6
         do ich=1,ichp1
            ippx=nint(res(jj))
            if (ippx.eq.1) then
               igamma=ich
            else
               i=i+1
               ipp(i,igroup,ier)=ippx
               if (ipp(i,igroup,ier).eq.2) ix=ix+1
               lspin(i,igroup,ier)=nint(res(jj+1))
               chspin(i,igroup,ier)=res(jj+2)
               bound(i,igroup,ier)=res(jj+3)
               rdeff(i,igroup,ier)=res(jj+4)
               if (mt(ipp(i,igroup,ier),ier).eq.2) ascat=res(jj+4)
               rdtru(i,igroup,ier)=res(jj+5)
            endif
            jj=jj+6
         enddo

         !--read the resonance parameters
         call listio(nin,0,0,res(jnow),nb,nw)
         jj=jnow+nw
         do while (nb.ne.0)
            call moreio(nin,0,0,res(jj),nb,nw)
            jj=jj+nw
            if (jj.gt.maxres) call error('rdsammy',&
              'res storage exceeded',' ')
         enddo
         nresg(igroup,ier)=nint(res(jnow+3))
         if (nresg(igroup,ier).gt.0) then
            jj=jnow+6
            nx=6
            if (2+ichan.gt.6) nx=12
            do ires=1,nresg(igroup,ier)
               kres=kres+1

               !--add center and half-height energies to nodes
               hw=res(jj+1)/2
               do j=1,ichan
                  hw=hw+abs(res(jj+1+j))/2
               enddo
               ndig=5
               if (res(jj).gt.zero) ndig=2+nint(log10(res(jj)/(hw/10)))
               if (ndig.lt.5) ndig=5
               if (ndig.gt.9) ndig=9
               if (res(jj).gt.el.and.res(jj).lt.eh) then
                 nodes=nodes+1
                 enode(nodes)=sigfig(res(jj),ndig,0)
               endif
               if ((res(jj)+hw).gt.el.and.(res(jj)+hw).lt.eh) then
                  if (nodes.ge.nodmax) call error('rdsammy',&
                    'storage in enode exceeded.',' ')
                  ehalf=res(jj)+hw
                  nodes=nodes+1
                  enode(nodes)=sigfig(ehalf,ndig,0)
               endif
               if ((res(jj)-hw).gt.el.and.(res(jj)-hw).lt.eh) then
                  if (nodes.ge.nodmax) call error('rdsammy',&
                    'storage in enode exceeded.',' ')
                  ehalf=res(jj)-hw
                  nodes=nodes+1
                  enode(nodes)=sigfig(ehalf,ndig,0)
               endif

               !--store parameters in sammy arrays
               eres(kres,ier)=res(jj)
               gamgam(kres,ier)=res(jj+1)
               do j=1,ichan
                  gamma(j,kres,ier)=res(jj+1+j)
               enddo
               if (igamma.ne.1) then
                  x=gamma(igamma+1,kres,ier)
                  if (igamma.gt.1) then
                     do i=igamma,2,-1
                        gamma(i,kres,ier)=gamma(i-1,kres,ier)
                     enddo
                  endif
                  gamma(1,kres,ier)=gamgam(kres,ier)
                  gamgam(kres,ier)=x
               endif
               jj=jj+nx
            enddo
         endif
      enddo
      if (ngroupm(ier).gt.1) then
         do igroup=2,ngroupm(ier)
            nresg(igroup,ier)=nresg(igroup-1,ier)+nresg(igroup,ier)
         enddo
      endif

   !--other formats not currently coded
   else
      call error('rdsammy','only RML and RM currently coded',' ')

   endif

   !--sort nodes into order
   call orders(enode,nodes)

   deallocate(idum)

   !--optional output of converted MF2/MT151
   if (imf2.eq.1) then
      allocate(a(12*nres+100))
      mfh=2
      mth=151
      a(1)=0
      a(2)=0
      a(3)=0
      a(4)=3
      a(5)=ngroupm(ier)
      a(6)=0
      call contio(0,nsyso,0,a,nb,nw)
      a(1)=0
      a(2)=0
      a(3)=nppm(ier)
      a(4)=0
      a(5)=12*nppm(ier)
      a(6)=2*nppm(ier)
      j=6
      do i=1,nppm(ier)
         a(j+1)=ema(i,ier)
         a(j+2)=emb(i,ier)
         a(j+3)=kza(i,ier)
         a(j+4)=kzb(i,ier)
         a(j+5)=spina(i,ier)
         a(j+6)=spinb(i,ier)
         a(j+7)=qqq(i,ier)
         a(j+8)=lpent(i,ier)
         a(j+9)=ishift(i,ier)
         a(j+10)=mt(i,ier)
         a(j+11)=pa(i,ier)
         a(j+12)=pb(i,ier)
         j=j+12
      enddo
      call listio(0,nsyso,0,a,nb,nw)
      ires1=0
      do i=1,ngroupm(ier)
         a(1)=sspin(i,ier)
         a(2)=parity(i,ier)
         a(3)=0
         a(4)=0
         a(5)=6*nchan(i,ier)+6
         a(6)=nchan(i,ier)+1
         a(7)=1
         a(8)=lspin(1,i,ier)
         a(9)=0
         a(10)=0
         a(11)=0
         a(12)=0
         do j=1,nchan(i,ier)
            a(13+6*(j-1))=ipp(j,i,ier)
            a(14+6*(j-1))=lspin(j,i,ier)
            a(15+6*(j-1))=chspin(j,i,ier)
            a(16+6*(j-1))=bound(j,i,ier)
            a(17+6*(j-1))=rdeff(j,i,ier)
            a(18+6*(j-1))=rdtru(j,i,ier)
         enddo
         call listio(0,nsyso,0,a,nb,nw)
         a(1)=0
         a(2)=0
         a(3)=0
         a(4)=nresg(i,ier)-ires1
         nx=6
         if (2+nchan(i,ier).gt.6) nx=12
         a(5)=nx*(nresg(i,ier)-ires1)
         a(6)=nresg(i,ier)-ires1
         ll=6
         do j=1+ires1,nresg(i,ier)
            a(ll+1)=eres(j,ier)
            a(ll+2)=gamgam(j,ier)
            do jj=1,nx-2
               if (jj.le.nchan(i,ier)) then
                  a(ll+2+jj)=gamma(jj,j,ier)
               else
                  a(ll+2+jj)=0
               endif
            enddo
            ll=ll+nx
         enddo
         call listio(0,nsyso,0,a,nb,nw)
         ll=1+nw
         do while (nb.ne.0)
            call moreio(0,nsyso,0,a(ll),nb,nw)
            ll=ll+nw
         enddo
         ires1=nresg(i,ier)
      enddo
      deallocate(a)
   endif

   return
   end subroutine rdsammy

   subroutine orders(x,n)
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
   end subroutine orders

   subroutine ppdefs(ier)
   !-------------------------------------------------------------------
   ! Set particle-pair defaults
   !-------------------------------------------------------------------
   use physics ! light particle masses
   ! externals
   integer::ier
   ! internals
   integer::i

   do i=1,nppm(ier)
      if (mt(i,ier).eq.102) then  ! one particle is gamma
      else if (mt(i,ier).eq.2) then  ! neutron
         ema(i,ier)=1
         spina(i,ier)=0.5e0_kr
         kza(i,ier)=0
         lpent(i,ier)=1
      else if (mt(i,ier).eq.18) then  ! fission
         lpent(i,ier)=0
      else if (mt(i,ier).gt.50.and.mt(i,ier).lt.99) then  ! inelastic
         ema(i,ier)=1
         spina(i,ier)=0.5e0_kr
         kza(i,ier)=0
         lpent(i,ier)=1
      else if (mt(i,ier).eq.103) then  ! proton
         ema(i,ier)=pnratio
         spina(i,ier)=0.5e0_kr
         kza(i,ier)=1
         lpent(i,ier)=1
      else if (mt(i,ier).eq.104) then  ! deuteron
         ema(i,ier)=dnratio
         spina(i,ier)=1
         kza(i,ier)=1
         lpent(i,ier)=1
      else if (mt(i,ier).eq.105) then  ! triton
         ema(i,ier)=tnratio
         spina(i,ier)=0.5e0_kr
         kza(i,ier)=1
         lpent(i,ier)=1
      else if (mt(i,ier).eq.106) then  ! he3
         ema(i,ier)=hnratio
         spina(i,ier)=0.5e0_kr
         kza(i,ier)=2
         lpent(i,ier)=1
      else if (mt(i,ier).eq.107) then  ! alpha
         ema(i,ier)=anratio
         spina(i,ier)=0
         pa(i,ier)=0
         kza(i,ier)=2
         lpent(i,ier)=1
      endif
   enddo

   return
   end subroutine ppdefs

   subroutine rearrange(idum,ier)
   !-------------------------------------------------------------------
   ! Rearrange to put resonances in spin-group order.
   !-------------------------------------------------------------------
   use util ! provides error
   ! externals
   integer::idum(nres)
   integer::ier
   ! internals
   integer::i,isort,j,k
   real(kr)::a

   do i=1,nresm(ier)
      isort=0
      do j=2,nresm(ier)
         if (idum(j).lt.idum(j-1)) then
            isort=1
            k=idum(j)
            idum(j)=idum(j-1)
            idum(j-1)=k
            a=eres(j,ier)
            eres(j,ier)=eres(j-1,ier)
            eres(j-1,ier)=a
            a=gamgam(j,ier)
            gamgam(j,ier)=gamgam(j-1,ier)
            gamgam(j-1,ier)=a
            do k=1,mchan
               a=gamma(k,j,ier)
               gamma(k,j,ier)=gamma(k,j-1,ier)
               gamma(k,j-1,ier)=a
            enddo
         endif
      enddo
      if (isort.eq.0) exit
   enddo
   do i=1,ngroupm(ier)
      nresg(i,ier)=0
   enddo
   do j=1,nresm(ier)
      nresg(idum(j),ier)=nresg(idum(j),ier)+1
   enddo
   if (ngroupm(ier).gt.1) then
      do i=2,ngroupm(ier)
         nresg(i,ier)=nresg(i-1,ier)+nresg(i,ier)
      enddo
   endif
   if (nresm(ier).ne.nresg(ngroup,ier)) then
      call error('rearrange','nres fault',' ')
   endif

   return
   end subroutine rearrange

   subroutine ppsammy(ier,ncoef,nresp)
   !-------------------------------------------------------------------
   ! Prepare for SAMMY calculation.
   !-------------------------------------------------------------------
   ! externals
   integer::ier,ncoef,nresp
   ! internals
   integer::i,j

   !--zero arrays
   do i=1,mchan
      do j=1,ngroup
         echan(i,j,ier)=0
         zeta(i,j,ier)=0
      enddo
   enddo

   !--check quantum numbers for particle pairs
   call checkqn(ier)

   !--generate parameters for calculating k, rho, and eta for sammy
   call fxradi(ier)

   !--generate energy-independent channel info for the sammy method
   call betset(ier)
   if (npar.gt.0) nresp=npar

   !--set up for angular distribution calculations
   ncoef=1
   if (Want_Angular_Dist) then
      call angle(ier)
      ncoef=lllmax
   endif

   !--generate energy-independent portion of derivatives
   if (Want_Partial_Derivs) call babb(ier)

   return
   end subroutine ppsammy

   subroutine allo
   !-------------------------------------------------------------------
   ! Allocate arrays for reading in resonance parameters
   ! and computing cross sections, derivatives, and angular
   ! distributions for the sammy method using sizes from s2sammy.
   !-------------------------------------------------------------------

   ntriag=(mchan*(mchan+1))/2

   allocate(ema(npp,nerm))
   allocate(emb(npp,nerm))
   allocate(kza(npp,nerm))
   allocate(kzb(npp,nerm))
   allocate(spina(npp,nerm))
   allocate(spinb(npp,nerm))
   allocate(qqq(npp,nerm))
   allocate(ishift(npp,nerm))
   allocate(lpent(npp,nerm))
   allocate(mt(npp,nerm))
   allocate(pa(npp,nerm))
   allocate(pb(npp,nerm))
   allocate(sspin(ngroup,nerm))
   allocate(parity(ngroup,nerm))
   allocate(nresg(ngroup,nerm))
   allocate(ipp(mchan,ngroup,nerm))
   allocate(lspin(mchan,ngroup,nerm))
   allocate(chspin(mchan,ngroup,nerm))
   allocate(bound(mchan,ngroup,nerm))
   allocate(rdeff(mchan,ngroup,nerm))
   allocate(rdtru(mchan,ngroup,nerm))
   allocate(eres(nres,nerm))
   allocate(gamgam(nres,nerm))
   allocate(gamma(mchan,nres,nerm))
   allocate(parext(7,mchan,ngroup,nerm))
   allocate(nent(nres,nerm))
   allocate(next(nres,nerm))
   allocate(goj(nres,nerm))
   allocate(zke(mchan,ngroup,nerm))
   allocate(zkfe(mchan,ngroup,nerm))
   allocate(zkte(mchan,ngroup,nerm))
   allocate(zeta(mchan,ngroup,nerm))
   allocate(echan(mchan,ngroup,nerm))
   allocate(betapr(mchan,nres,nerm))
   allocate(gbetpr(3,nres,nerm))
   allocate(beta(ntriag,nres,nerm))
   allocate(uuuu((mchan+2)*nres,nerm))
   allocate(duuu((mchan+2)*nres,nerm))
   allocate(iduu((mchan+2)*nres,nerm))

   if (Want_Partial_Derivs) then
      allocate(br(ntriag,npar,nerm))
      allocate(bi(ntriag,npar,nerm))
      allocate(pr(ntriag,npar,nerm))
      allocate(pii(ntriag,npar,nerm))
   endif

   allocate(crss(npp,nerm))
   allocate(sigmas(npp,nerm))

   if (Want_Partial_Derivs) then
      allocate(deriv(npp,npar,nerm))
      allocate(dsigma(npp,npar,nerm))
   endif

   if (Want_Angular_Dist) then
      allocate(crssx(2,mchan,mchan,npp,ngroup,nerm))
      allocate(Coef_Leg(lllmax,npp,nerm))
      if (Want_Partial_Derivs) then
         allocate(derivx(2,mchan,mchan,1,nerm))
         allocate(D_Coef_Leg(lllmax,npp,npar,nerm))
      endif
   endif

   allocate(alphar(nres,nerm))
   allocate(alphai(nres,nerm))
   allocate(difen(nres,nerm))
   allocate(xden(nres,nerm))

   allocate(sinsqr(mchan,nerm))
   allocate(sin2ph(mchan,nerm))
   allocate(dphi(mchan,nerm))
   allocate(dpdr(mchan,nerm))
   allocate(dsdr(mchan,nerm))
   allocate(cosphi(mchan,nerm))
   allocate(sinphi(mchan,nerm))

   allocate(rootp(mchan,nerm))
   allocate(elinvr(mchan,nerm))
   allocate(elinvi(mchan,nerm))
   allocate(psmall(mchan,nerm))

   allocate(xqr(mchan,mchan,nerm))
   allocate(xqi(mchan,mchan,nerm))

   allocate(xxxxr(ntriag,nerm))
   allocate(xxxxi(ntriag,nerm))

   allocate(qr(ntriag,ntriag,nerm))
   allocate(qi(ntriag,ntriag,nerm))

   allocate(rmat(2,ntriag,nerm))
   allocate(ymat(2,ntriag,nerm))
   allocate(yinv(2,ntriag,nerm))
   allocate(cscs(2,ntriag,nerm))

   if (Want_Partial_Derivs) then
      allocate(upr(npar,nerm))
      allocate(upi(npar,nerm))
   endif

   if (Want_Partial_Derivs) then
      allocate(tr(npp,ntriag,nerm))
      allocate(ti(npp,ntriag,nerm))
      allocate(tx(2,ntriag,ntriag,nerm))
   endif

   allocate(ddddd(npp,nerm))

   allocate(par(npar))

   return
   end subroutine allo

   subroutine angle(ier)
   !-------------------------------------------------------------------
   ! Organize the calculation of angular distributions
   !-------------------------------------------------------------------
   ! externals
   integer::ier

   !--generate coefficients of Legendre polynomicals for the various
   !--spin group parts.  count how many non-zero elements are in xlmn.
   allocate(alj(ntriag,ngroup))
   call lmaxxx(alj,ier)
   call kclbsch(alj,ier)
   allocate(xlmn(kkxlmn))
   allocate(kxlmn(kkxlmn))
   call clbsch(alj,ier)
   deallocate(alj)

   return
   end subroutine angle

   subroutine findsp(spin,lspin,chspin,llll,spinjj,s1,s2,ig,ier)
   !-------------------------------------------------------------------
   ! Find the group number ig for the given llll and spinjj.
   !-------------------------------------------------------------------
   use util ! provides error
   ! external
   real(kr)::spin(*),chspin(mchan,ngroup,*)
   integer::lspin(mchan,ngroup,*)
   integer llll,ig,ier
   real(kr)::spinjj,s1,s2
   ! internal
   integer::i
   real(kr)::a
   integer::ifneg=0
   real(kr),parameter::zero=0

   ig=0
   do i=1,ngroup
      a=abs(spin(i))
      if (llll.eq.lspin(1,i,ier)) then
         if (spinjj.eq.a) then
            if (chspin(1,i,ier).eq.s2) then
               ig=i
               return
            else if (spinjj.eq.zero) then
               ig=i
               return
            else if (ifneg.eq.0) then
               ig=i
               return
            endif
         else if (spinjj.eq.-a) then
            if (chspin(1,i,ier).eq.s1) then
               ig=i
               ifneg=1
               return
            endif
         endif
      endif
   enddo
   if (ifneg.eq.1.and.ig.eq.0) then
      call error('findsp','quantum numbers in file 2 do not make sense',' ')
   endif

   return
   end subroutine findsp

   subroutine checkqn(ier)
   !-------------------------------------------------------------------
   ! Report quantum number information etc when particle-pair
   ! definitions are given.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use util   ! provides error
   ! externals
   integer::ier
   ! internals
   integer::j,n,i,nchanj
   real(kr)::x,smin,smax
   character(60)::strng
   real(kr),parameter::zero=0
   real(kr),parameter::half=0.5e0_kr

   do j=1,ngroupm(ier)

      if (mod(sspin(j,ier),half).ne.zero) then
         write(strng,'(''spin group number '',I3,'' with spin ='',&
           &F10.5)') j,sspin(j,ier)
         call error('checkqn','error in quantum numbers',strng)
      endif

      goj(j,ier)=0
      nent(j,ier)=0
      next(j,ier)=0
      nchanj=nchan(j,ier)
      do n=1,nchanj
         if (ipp(n,j,ier).eq.2) then
            nent(j,ier)=nent(j,ier)+1
            if (goj(j,ier).eq.zero) then
               goj(j,ier)=(2*abs(sspin(j,ier))+1)/&
                 ((2*abs(spina(ipp(n,j,ier),ier))+1)&
                 *(2*abs(spinb(ipp(n,j,ier),ier))+1))
            endif
         else if (ipp(n,j,ier).gt.2) then
            next(j,ier)=next(j,ier)+1
         endif
         if (lspin(n,j,ier).lt.0) call wrongi('lspin',lspin(n,j,ier))

         if (chspin(n,j,ier).eq.zero.or.sspin(j,ier).eq.zero) then
         else
            x=1
            i=mod(lspin(n,j,ier),2)
            if (i.ne.0) x=-x
            if (chspin(n,j,ier).lt.zero) x=-x
            if (sspin(j,ier).lt.zero) x=-x
            if (x.lt.zero) then
               write(nsyso,'('' *** Parity problem ***'',/,&
                 &''     Group and channel #'',2i5,/,&
                 &''     Spin, L, Chspin ='',f7.1,i4,f7.1)')&
                 j,n,sspin(j,ier),lspin(n,j,ier),chspin(n,n,ier)
            endif
         endif
         smin=abs(abs(spina(ipp(n,j,ier),ier))-abs(spinb(ipp(n,j,ier),ier)))
         smax=abs(spina(ipp(n,j,ier),ier))+abs(spinb(ipp(n,j,ier),ier))
         if (abs(chspin(n,j,ier)).lt.smin) then
             write(nsyso,'(/,'' ****** Error in quantum numbers for'',&
               &1x,''Group Number'',i3,'' and Channel Number'',i2)')&
               j,n
             write(nsyso,&
               '('' ******* |Chspin|<|Spina-Spinb|=Smin'',3x,4f5.1)')&
               chspin(n,j,ier),spina(ipp(n,j,ier),ier),&
                 spinb(ipp(n,j,ier),ier),smin
         else if (abs(chspin(n,j,ier)).gt.smax) then
             write(nsyso,'(/,'' ****** Error in quantum numbers for'',&
               &1x,''Group Number'',i3,'' and Channel Number'',i2)')&
               j,n
             write(nsyso,&
               '('' ******* |Chspin|>|Spina+Spinb|=Smax'',3x,4f5.1)')&
               chspin(n,j,ier),spina(ipp(n,j,ier),ier),&
               spinb(ipp(n,j,ier),ier),smax
         endif
         if (lpent(ipp(n,j,ier),ier).gt.0) then
            smin=abs(abs(float(lspin(n,j,ier))-abs(chspin(n,j,ier))))
            smax=lspin(n,j,ier)+abs(chspin(n,j,ier))
            if (abs(sspin(j,ier)).lt.smin) then
               write(nsyso,'(/,'' ****** Error in quantum numbers for'',&
                 &1x,''Group Number'',i3,'' and Channel Number'',i2)')&
                 j,n
               write(nsyso,&
                 '('' ****** |Spin|<|Lspin-Chspin|=Smin'',3x,f5.1,i5,2f5.1)')&
                 sspin(j,ier),lspin(n,j,ier),chspin(n,j,ier),smin
            else if (abs(sspin(j,ier)).gt.smax) then
               write(nsyso,'(/,'' ****** Error in quantum numbers for'',&
                 &1x,''Group Number'',i3,'' and Channel Number'',i2)')&
                 j,n
               write(nsyso,&
                 '('' ****** |Spin|>|Lspin-chspin|=Smin'',3x,f5.1,i5,2f5.1)')&
                 sspin(j,ier),lspin(n,j,ier),chspin(n,j,ier),smin
            endif
         endif
      enddo
   enddo

   return
   end subroutine checkqn

   subroutine fxradi(ier)
   !-------------------------------------------------------------------
   ! Generate parameters for calculating k, rho, and eta for sammy
   !-------------------------------------------------------------------
   use physics ! get global physics constants
   ! externals
   integer::ier
   ! internals
   integer::ippx,kgroup,ichan
   real(kr)::ff,twomhb,etac
   real(kr)::factor,alabcm,aa,redmas,z
   real(kr),parameter::hbarrr=hbar/ev              !hbar in eV-s
   real(kr),parameter::amuevv=amu*clight*clight/ev !amu in eV
   real(kr),parameter::cspeed=clight/100           !c in m/s (not cm/s!)
   real(kr),parameter::fininv=1.e0_kr/finstri      !fine structure constant
   real(kr),parameter::zero=0

   ff=1.e+15_kr
   twomhb=sqrt(2*amassn*amuevv)/(hbarrr*ff*cspeed)
   etac=fininv*amuevv/(hbarrr*ff*cspeed)*amassn

   do ippx=1,nppm(ier)
      if (ema(ippx,ier).eq.zero) then
         ema(ippx,ier)=1
      endif
      if (emb(ippx,ier).eq.zero) then
          emb(ippx,ier)=awr
      endif
   enddo

   do kgroup=1,ngroupm(ier)
      factor=emb(2,ier)+ema(2,ier)
      alabcm=emb(2,ier)/factor
      factor=alabcm/ema(2,ier)
      do ichan=1,nchan(kgroup,ier)
         ippx=ipp(ichan,kgroup,ier)
         aa=emb(ippx,ier)+ema(ippx,ier)
         aa=emb(ippx,ier)/aa
         if (qqq(ippx,ier).ne.zero.and.echan(ichan,kgroup,ier).eq.zero) then
            echan(ichan,kgroup,ier)=-qqq(ippx,ier)/alabcm
         endif
         redmas=aa*ema(ippx,ier)
         z=twomhb*sqrt(redmas*factor)
         zke(ichan,kgroup,ier)=z
         zkfe(ichan,kgroup,ier)=z*rdeff(ichan,kgroup,ier)*10
         zkte(ichan,kgroup,ier)=z*rdtru(ichan,kgroup,ier)*10
         if (kzb(ippx,ier).ne.0.and.kza(ippx,ier).ne.0)then
            zeta(ichan,kgroup,ier)=&
              etac*kzb(ippx,ier)*kza(ippx,ier)*redmas/zke(ichan,kgroup,ier)
         endif
      enddo
   enddo

   return
   end subroutine fxradi

   subroutine betset(ier)
   !-------------------------------------------------------------------
   ! Generate energy-independent channel info for the sammy method
   !-------------------------------------------------------------------
   use util ! provides mess
   ! externals
   integer::ier
   ! internals
   integer::knpar,minres,jdopha,ig,j,lsp,ires,iffy
   integer::kl,k,l,mmaxc,mmaxc2,m,m2,kgroup,kqq,ipx
   real(kr)::rho,drho,p,dp,ds,er,ex,eta,sinphi,cosphi,dphi,hr,hi
   real(kr),dimension(:,:),allocatable::dum
   character(60)::str1
   real(kr),parameter::zero=0
   real(kr),parameter::half=0.5e0_kr

   !--convert to betapr parameters
   allocate(dum(mchan,nres))
   knpar=0
   if (nres.gt.0) then

      minres=1
      jdopha=0
      do ig=1,ngroupm(ier)

         mmaxc=nchan(ig,ier)
         do j=1,mmaxc
            lsp=lspin(j,ig,ier)
            ipx=ipp(j,ig,ier)
            do ires=minres,nresg(ig,ier)
               betapr(j,ires,ier)=0
               dum(j,ires)=0
               rho=0
               p=1
               if (gamma(j,ires,ier).ne.zero) then
                  if (lpent(ipx,ier).le.0) then
                     betapr(j,ires,ier)=sqrt(half*abs(gamma(j,ires,ier)))
                     if (gamma(j,ires,ier).lt.zero)&
                       betapr(j,ires,ier)=-betapr(j,ires,ier)
                  else
                     ex=abs(eres(ires,ier)-echan(j,ig,ier))
                     if (ex.ne.zero) then
                        ex=sqrt(ex)
                        rho=zkte(j,ig,ier)*ex
                        er=eres(ires,ier)
                        if (er.lt.zero) er=-er
                        drho=zkte(j,ig,ier)*sqrt(er)/ex
                        if (zeta(j,ig,ier).eq.zero) then
                           p=pf(rho,lsp)
                           call pgh(rho,lsp,bound(j,ig,ier),hr,hi,p,&
                             dp,ds,ishift(ipx,ier),iffy)
                        else
                           eta=zeta(j,ig,ier)/ex
                           jdopha=0
                           call pghcou(rho,lsp,bound(j,ig,ier),hr,hi,p,dp,ds,&
                             ishift(ipx,ier),iffy,eta,sinphi,cosphi,dphi,jdopha)
                           if (p.lt.zero) then
                              p=1
                              dp=0
                              write(str1,'(''p<0 set to 1 at eres=''&
                                &,1p,e12.4)') eres(ires,ier)
                              call mess('betset',str1,' ')
                           else if (p.eq.zero) then
                              p=1
                              dp=0
                              write(str1,'(''p=0 set to 1 at eres='',&
                                &1p,e12.4)') eres(ires,ier)
                              call mess('betset',str1,' ')
                           endif
                        endif
                     endif
                  endif
                  betapr(j,ires,ier)=sqrt(half*abs(gamma(j,ires,ier))/p)
                  if (gamma(j,ires,ier).lt.zero)&
                    betapr(j,ires,ier)=-betapr(j,ires,ier)
                  dum(j,ires)=dp*drho*betapr(j,ires,ier)/(2*p)
                  if (eres(ires,ier).lt.zero) dum(j,ires)=-dum(j,ires)
               endif
            enddo
         enddo

         !--generate beta parameters
         do ires=minres,nresg(ig,ier)
            kl=0
            do k=1,mmaxc
               do l=1,k
                  kl=kl+1
                  beta(kl,ires,ier)=betapr(l,ires,ier)*betapr(k,ires,ier)
               enddo
            enddo
         enddo

         do ires=minres,nresg(ig,ier)
            gbetpr(2,ires,ier)=abs(gamgam(ires,ier))/2
            gbetpr(1,ires,ier)=sqrt(gbetpr(2,ires,ier))
            if (gamgam(ires,ier).lt.zero) gbetpr(1,ires,ier)=-gbetpr(1,ires,ier)
            gbetpr(3,ires,ier)=gbetpr(2,ires,ier)**2
         enddo

         !--convert from u parameters if needed
         if (Want_Partial_Derivs) then
            mmaxc2=mmaxc+2
            do ires=minres,nresg(ig,ier)
               do m=1,mmaxc2
                  knpar=knpar+1
                  iduu(knpar,ier)=0
                  if (Want_Partial_U) then
                     uuuu(knpar,ier)=1
                     duuu(knpar,ier)=0
                  else
                     if (m.eq.1) then
                        uuuu(knpar,ier)=sqrt(abs(eres(ires,ier)))
                        uuuu(knpar,ier)=2*uuuu(knpar,ier)
                        duuu(knpar,ier)=0
                        iduu(knpar,ier)=1
                     else if (m.eq.2) then
                        uuuu(knpar,ier)=4*abs(gbetpr(1,ires,ier))
                        duuu(knpar,ier)=0
                     else if (m.ge.3) then
                        m2=m-2
                        if (gamma(m2,ires,ier).ne.zero) then
                           uuuu(knpar,ier)=2*gamma(m2,ires,ier)/&
                             betapr(m2,ires,ier)
                           duuu(knpar,ier)=dum(m2,ires)
                        else
                           uuuu(knpar,ier)=0
                           duuu(knpar,ier)=0
                        endif
                     endif
                  endif
               enddo
            enddo
         endif
         minres=nresg(ig,ier)+1
      enddo
   endif
   nparr=knpar

   !--rexternal parameters
   if (nrext.gt.0.and.Want_Partial_Derivs) then
      do kgroup=1,ngroupm(ires)
         mmaxc=nchan(kgroup,ier)
         do j=1,mmaxc
            do kqq=1,7
               knpar=knpar+1
               iduu(knpar,ier)=0
               if (Want_Partial_U) then
                  uuuu(knpar,ier)=1
               else
                  if (kqq.ne.5) uuuu(knpar,ier)=parext(kqq,j,kgroup,ier)
                  if (kqq.eq.5) uuuu(knpar,ier)=sqrt(parext(kqq,j,kgroup,ier))
               endif
            enddo
         enddo
      enddo
   endif

   npar=knpar
   deallocate(dum)

   return
   end subroutine betset

   subroutine lmaxxx(alj,ier)
   !-------------------------------------------------------------------
   ! Modified from subroutine clbsch to give more accurate answer.
   !-------------------------------------------------------------------
   use util ! provides error
   ! externals
   real(kr)::alj(ntriag,*)
   integer::ier
   ! internals
   integer::kountll(51)
   integer::max,n,ich,i,j,nn,kn,nc,ncx,lln,jn,llnx,jnx,m,mm
   integer::mc,llm,lmin,lmax,km,mcx,llmx,lminx,lmaxx,nmlmin,nmlmax
   integer::ll,l
   real(kr)::spinn,eln,b,bx,cjn,acn,elnx,cjnx,spinm,acn2,elm
   real(kr)::elmx,acnm,d,c,a,dx,cx,el
   real(kr),parameter::zero=0
   real(kr),parameter::eps=0.001e0_kr

   max=0
   do n=1,ngroupm(ier)
      do ich=1,nchan(n,ier)
         if (lspin(ich,n,ier).gt.max) max=lspin(ich,n,ier)
      enddo
   enddo

   lllmax=2*max+1
   if (lllmax.gt.51) call error('lmaxxx','lllmax limited to 51',' ')
   do i=1,ntriag
      do j=1,ngroupm(ier)
         alj(i,j)=0
      enddo
   enddo
   do i=1,51
      kountll=0
   enddo

   do n=1,ngroupm(ier)
      nn=nchan(n,ier)
      spinn=abs(sspin(n,ier))
      kn=0
      do nc=1,nn
         eln=lspin(nc,n,ier)
         b=xdel2(eln,abs(chspin(nc,n,ier)),spinn)
         do ncx=1,nc
            kn=kn+1
            if (b.ne.zero) then
               if (lspin(ncx,n,ier).eq.lspin(nc,n,ier).and.&
                 abs(abs(chspin(ncx,n,ier))-abs(chspin(nc,n,ier))).le.eps) then
                  alj(kn,n)=b
               else
                  elnx=lspin(ncx,n,ier)
                  bx=xdel2(elnx,abs(chspin(ncx,n,ier)),spinn)
                  if (bx.ne.zero) then
                     alj(kn,n)=bx
                  endif
               endif
            endif
         enddo
      enddo
   enddo

   do n=1,ngroupm(ier)
      nn=nchan(n,ier)
      spinn=abs(sspin(n,ier))
      kn=0
      do nc=1,nn
        lln=lspin(nc,n,ier)
        eln=lspin(nc,n,ier)
        cjn=abs(chspin(nc,n,ier))
        jn=nint(2*cjn)
        do ncx=1,nc
           kn=kn+1
           acn=alj(kn,n)
           if (acn.ne.zero) then
              llnx=lspin(ncx,n,ier)
              elnx=lspin(ncx,n,ier)
              cjnx=abs(chspin(ncx,n,ier))
              jnx=nint(2*cjnx)
              do m=1,n
                 mm=nchan(m,ier)
                 spinm=abs(sspin(m,ier))
                 acn2=acn
                 do mc=1,mm
                    if (abs(abs(chspin(mc,m,ier))-cjn).le.zero) then
                       llm=lspin(mc,m,ier)
                       elm=lspin(mc,m,ier)
                       lmin=iabs(lln-llm)+1
                       lmax=lln+llm+1
                       if (lmax.gt.lllmax) lmax=lllmax
                       if (lmax.ge.lmin) then
                          km=(mc*(mc-1))/2
                          do mcx=1,mc
                             km=km+1
                             if (abs(abs(chspin(mcx,m,ier))-cjnx).le.eps) then
                                llmx=lspin(mcx,m,ier)
                                elmx=lspin(mcx,m,ier)
                                lminx=iabs(llnx-llmx)+1
                                lmaxx=llnx+llmx+1
                                if (lmaxx.gt.lllmax) lmaxx=lllmax
                                if (lmaxx.ge.lminx) then
                                   acnm=alj(km,m)*acn2
                                   nmlmin=lmin
                                   nmlmax=lmax
                                   if (nmlmin.lt.lminx) nmlmin=lminx
                                   if (nmlmax.gt.lmaxx) nmlmax=lmaxx
                                   do ll=nmlmin,nmlmax,2
                                      l=ll-1
                                      el=l
                                      b=xdel2(spinn,spinm,el)
                                      if (b.ne.zero) then
                                         d=xdel2(eln,elm,el)
                                         c=xww(eln,spinn,elm,spinm,cjn,el)
                                         if (jn.eq.jnx.and.lln.eq.llnx.and.&
                                           llm.eq.llmx) then
                                            a=(b*d*c)**2*acnm
                                            if (a.ne.zero) kountll(ll)=1
                                         else
                                            bx=b
                                            dx=d
                                            cx=c
                                            if (lln.ne.llm.or.llnx.ne.llmx) then
                                               dx=xdel2(elnx,elmx,el)
                                            endif
                                            cx=xww(elnx,spinn,elmx,spinm,&
                                              cjnx,el)
                                            a=b*d*c*bx*dx*cx*acnm
                                            if (a.ne.zero) kountll(ll)=1
                                         endif
                                      endif
                                   enddo
                                endif
                             endif
                          enddo
                       endif
                    endif
                 enddo
              enddo
           endif
        enddo
      enddo
   enddo

   l=0
   do ll=1,lllmax
      if (kountll(ll).eq.1) l=ll
   enddo
   lllmax=l
   if (lllmax.gt.lllmaxx) then
      call error('lmaxxx','need larger lllmaxxx','')
   endif

   return
   end subroutine lmaxxx

   real(kr) function xdel2(a,b,c)
   !-------------------------------------------------------------------
   ! Private routine for lmaxxx.
   !-------------------------------------------------------------------
   ! externals
   real(kr)::a,b,c
   ! internals
   real(kr)::d,e,f
   d=a+b-c
   e=a-b+c
   f=-a+b+c
   if (nint(d).ge.0.and.nint(e).ge.0.and.nint(f).ge.0) then
      xdel2=1
   else
      xdel2=0
   endif
   return
   end function xdel2

   real(kr) function xww(el1,aj1,el2,aj2,cj,el)
   !-------------------------------------------------------------------
   ! Private routine for lmaxxx.
   !-------------------------------------------------------------------
   ! externals
   real(kr)::el1,aj1,el2,aj2,cj,el
   ! internals
   integer::i1,i2,m1,m2,minz,k1,k2,k3,maxz

   i1=nint(el1+aj1+cj)
   i2=nint(el2+aj2+cj)
   m1=nint(el1+el2+el)
   m2=nint(aj1+aj2+el)
   minz=max0(i1,i2,m1,m2)
   k1=nint(el1+aj1+el2+aj2)
   k2=nint(el1+aj2+cj+el)
   k3=nint(el2+aj1+cj+el)
   maxz=min0(k1,k2,k3)
   if (maxz.ge.minz) then
      xww=1
   else
      xww=0
   endif
   return
   end function xww

   subroutine kclbsch(alj,ier)
   !-------------------------------------------------------------------
   ! generate alj, the count of non-zero elements of xlmn
   !-------------------------------------------------------------------
   ! externals
   real(kr)::alj(ntriag,*)
   integer::ier
   ! internals
   integer::i,j,ngr,nchann,kn,nc,ncx,kount,lln,jn,llnx,jnx
   integer::mgr,nchanm,mc,llm,lmin,lmax,km,mcx,llmx
   integer::lminx,lmaxx,nmlmin,nmlmax,ll,l
   real(kr)::spinn,aj,eln,b,a,elnx,bx,cjn,acn,cjnx
   real(kr)::spinm,elm,elmx,el,d,c,cx,dx
   real(kr),parameter::zero=0
   real(kr),parameter::eps=0.001e0_kr
   real(kr),parameter::tenth=0.1e0_kr

   do i=1,ntriag
      do j=1,ngroup
         alj(i,j)=0
      enddo
   enddo

   do ngr=1,ngroup
      nchann=nchan(ngr,ier)
      spinn=abs(sspin(ngr,ier))
      aj=2*spinn+1
      kn=0
      do nc=1,nchann
         eln=lspin(nc,ngr,ier)
         b=del2(eln,abs(chspin(nc,ngr,ier)),spinn)
         eln=2*eln+1
         a=eln*aj*tenth/zke(1,ngr,ier)
         do ncx=1,nc
            kn=kn+1
            if (ncx.le.nent(ngr,ier)) then
               if (b.ne.zero) then
                  if (lspin(ncx,ngr,ier).eq.lspin(nc,ngr,ier).and.&
                    abs(abs(chspin(ncx,ngr,ier))-abs(chspin(nc,ngr,ier)))&
                    .le.eps) then
                     alj(kn,ngr)=a*b
                     if (nc.le.nent(ngr,ier).and.ncx.ne.nc)&
                       alj(kn,ngr)=2*alj(kn,ngr)
                  else
                     elnx=lspin(ncx,ngr,ier)
                     bx=del2(elnx,abs(chspin(ncx,ngr,ier)),spinn)
                     if (bx.ne.zero) then
                        elnx=2*elnx+1
                        bx=a*b*elnx*aj*bx*tenth/zke(ncx,ngr,ier)
                        alj(kn,ngr)=sqrt(bx)
                        if (nc.le.nent(ngr,ier).and.ncx.ne.nc)&
                          alj(kn,ngr)=2*alj(kn,ngr)
                     endif
                  endif
               endif
            endif
         enddo
      enddo
   enddo

   ! Alj = (2*J+1) * (2*el+1) *del2 (el, j, J) in one-channel case
   ! Alj = (2*J+1) * sqrt { (2*el +1) *del2 (el , j , J) } * 2
   !               * sqrt { (2*el'+1) *del2 (el', j', J) } in 2-channel

   !--Done calculating Alj.  Now count non-zero elements of Xlmn

   kount=0
   do ngr=1,ngroup
      nchann=nchan(ngr,ier)
      spinn=abs(sspin(ngr,ier))
      kn=0
      do nc=1,nchann
         lln=lspin(nc,ngr,ier)
         eln=float(lspin(nc,ngr,ier))
         cjn=abs(chspin(nc,ngr,ier))
         jn=nint(2*cjn)
         do ncx=1,nc
            kn=kn+1
            if (ncx.le.nent(ngr,ier)) then
               acn=alj(kn,ngr)
               if (acn.ne.zero) then
                  llnx=lspin(ncx,ngr,ier)
                  elnx=float(lspin(ncx,ngr,ier))
                  cjnx=abs(chspin(ncx,ngr,ier))
                  jnx=nint(2*cjnx)
                  do mgr=1,ngr
                     nchanm=nchan(mgr,ier)
                     spinm=abs(sspin(mgr,ier))
                     do mc=1,nchanm
                        if (abs(abs(chspin(mc,mgr,ier))-cjn).le.eps) then
                           llm=lspin(mc,mgr,ier)
                           elm=lspin(mc,mgr,ier)
                           lmin=iabs(lln-llm)+1
                           lmax=lln+llm+1
                           if (lmax.gt.lllmax) lmax=lllmax
                           if (lmax.ge.lmin) then
                              km=(mc*(mc-1))/2
                              do mcx=1,mc
                                 km=km+1
                                 if (mcx.le.nent(mgr,ier)) then
          !<-----------------------------------
          if (abs(abs(chspin(mcx,mgr,ier))-cjnx).le.eps) then
             llmx=lspin(mcx,mgr,ier)
             elmx=float(lspin(mcx,mgr,ier))
             lminx=iabs(llnx-llmx)+1
             lmaxx=llnx+llmx+1
             if (lmaxx.gt.lllmax) lmaxx=lllmax
             if (lmaxx.ge.lminx) then
                nmlmin=lmin
                nmlmax=lmax
                if (nmlmin.lt.lminx) nmlmin=lminx
                if (nmlmax.gt.lmaxx) nmlmax=lmaxx
                do ll=nmlmin,nmlmax,2
                   l=ll-1
                   el=float(l)
                   b=del2(spinn,spinm,el)
                   if (b.ne.zero) then
                      a=gll(lln,llm,l)
                      if (a.ne.zero) then
                         d=del2(eln,elm,el)
                         if (d.ne.zero) then
                            c=ww(eln,spinn,elm,spinm,cjn,el)
                            if (c.ne.zero) then
                               if (jn.eq.jnx.and.lln.eq.llnx&
                                 .and.llm.eq.llmx) then
                                  !--Here for diagonal case
                                  kount=kount+1
                               else
                                  !--Here for off-diagonal case
                                  bx=b
                                  cx=c
                                  dx=d
                                  if (lln.ne.llnx.or.llm.ne.llmx) then
                                     bx=gll(llnx,llmx,l)
                                     dx=del2(elnx,elmx,el)
                                  endif
                                  if (bx.ne.zero) then
                                     if (dx.ne.zero) then
                                        cx=ww(elnx,spinn,elmx,spinm,cjnx,el)
                                        if (cx.ne.zero) then
                                           kount=kount+1
                                        endif
                                     endif
                                  endif
                               endif
                            endif
                         endif
                      endif
                   endif
                enddo
             endif
          endif
          !---------------------------------->
                                 endif
                              enddo
                           endif
                        endif
                     enddo
                  enddo
               endif
            endif
         enddo
      enddo
   enddo

   kkxlmn=kount
   return
   end subroutine kclbsch

   subroutine clbsch(alj,ier)
   !-------------------------------------------------------------------
   ! Generate portions of Clebsch-Gordan-ry for combining two
   ! spin groups to give differential cross sections.
   !-------------------------------------------------------------------
   use util ! provides error
   ! externals
   real(kr)::alj(ntriag,ngroup)
   integer::ier
   ! internals
   integer::kount,ngr,nchann,kn,nc,lln,jn,ncx,llnx,jnx
   integer::mgr,nchanm,mc,llm,lmin,lmax,km,mcx,llmx,lminx,lmaxx
   integer::nmlmin,nmlmax,l,k,kk,ll,i,j,m,jxlmn
   real(kr)::ggginc,spinn,aaa,eln,cjn,acn,elnx,cjnx,acnm
   real(kr)::spinm,acn2,elm,elmx,b,a,d,c,bx,dx,cx,el
   real(kr),parameter::zero=0
   real(kr),parameter::eps=0.001e0_kr
   !integer::ntriag,ngroup
   jxlmn(m,l,k,j,i)=&
     ((((i-1)*ntriag+j-1)*ngroup+k-1)*ntriag+l-1)*lllmax+m

   do i=1,kkxlmn
      xlmn(i)=0
      kxlmn(i)=0
   enddo

   kount=0
   ggginc=2*abs(spina(2,ier))+1
   ggginc=1/ggginc
   do ngr=1,ngroup
      nchann=nchan(ngr,ier)
      spinn=abs(sspin(ngr,ier))
      aaa=2*abs(spinb(2,ier))+1
      aaa=ggginc/aaa
      kn=0
      do nc=1,nchann
         lln=lspin(nc,ngr,ier)
         eln=float(lspin(nc,ngr,ier))
         cjn=abs(chspin(nc,ngr,ier))
         jn=nint(2*cjn)
         do ncx=1,nc
            kn=kn+1
            if (ncx.le.nent(ngr,ier)) then
               acn=alj(kn,ngr)
               if (acn.ne.zero) then
                  acn=aaa*acn
                  llnx=lspin(ncx,ngr,ier)
                  elnx=float(lspin(ncx,ngr,ier))
                  cjnx=abs(chspin(ncx,ngr,ier))
                  jnx=nint(2*cjnx)
                  do mgr=1,ngr
                     nchanm=nchan(mgr,ier)
                     spinm=abs(sspin(mgr,ier))
                     acn2=acn
                     if (mgr.ne.ngr) acn2=2*acn
                     do mc=1,nchanm
                        if (abs(abs(chspin(mc,mgr,ier))-cjn).le.eps) then
                           llm=lspin(mc,mgr,ier)
                           elm=float(lspin(mc,mgr,ier))
                           lmin=iabs(lln-llm)+1
                           lmax=lln+llm+1
                           if (lmax.gt.lllmax) lmax=lllmax
                           if (lmax.ge.lmin) then
                              km=(mc*(mc-1))/2
                              do mcx=1,mc
                                 km=km+1
                                 if (mcx.le.nent(mgr,ier)) then
          !<---------------------------------------------
          if (abs(abs(chspin(mcx,mgr,ier))-cjnx).le.eps) then
             llmx=lspin(mcx,mgr,ier)
             elmx=lspin(mcx,mgr,ier)
             lminx=iabs(llnx-llmx)+1
             lmaxx=llnx+llmx+1
             if (lmaxx.gt.lllmax) lmaxx=lllmax
             if (lmaxx.ge.lminx) then
                acnm=alj(km,mgr)*acn2
                nmlmin=lmin
                nmlmax=lmax
                if (nmlmin.lt.lminx) nmlmin=lminx
                if (nmlmax.gt.lmaxx) nmlmax=lmaxx
                do ll=nmlmin,nmlmax,2
                   l=ll-1
                   el=float(l)
                   b=del2(spinn,spinm,el)
                   if (b.ne.zero) then
                      a=2*el+1
                      a=a*b
                      b=gll(lln,llm,l)
                      if (b.ne.zero) then
                         d=del2(eln,elm,el)
                         if (d.ne.zero) then
                            c=ww(eln,spinn,elm,spinm,cjn,el)
                            if (c.ne.zero) then
                               if (jn.eq.jnx.and.lln.eq.llnx&
                                 .and.llm.eq.llmx) then
                                  ! diagonal case
                                  a=a*(b*d*c)**2*acnm
                                  kount=kount+1
                                  xlmn(kount)=a
                                  kxlmn(kount)=jxlmn(ll,km,mgr,kn,ngr)
                               else
                                  ! off diagonal case
                                  bx=b
                                  dx=d
                                  cx=c
                                  if (jn.ne.jnx) then
                                     k=iabs(jn-jnx)/2
                                     kk=(k/2)*2
                                     if (k.ne.kk) a=-a
                                  endif
                                  if (lln.ne.llnx.or.llm.ne.llmx) then
                                     bx=gll(llnx,llmx,l)
                                     dx=del2(elnx,elmx,el)
                                  endif
                                  if (bx.ne.zero) then
                                     if (dx.ne.zero) then
                                        cx=ww(elnx,spinn, elmx,spinm,cjnx,el)
                                        if (cx.ne.zero) then
                                           a=a*b*d*c*bx*dx*cx*acnm
                                           kount=kount+1
                                           xlmn(kount)=a
                                           kxlmn(kount)=jxlmn(ll,km,mgr,kn,ngr)
                                        endif
                                     endif
                                  endif
                               endif
                            endif
                         endif
                      endif
                   endif
                enddo
             endif
          endif
          !--------------------------------->
                                 endif
                              enddo
                           endif
                        endif
                     enddo
                  enddo
               endif
            endif
         enddo
      enddo
   enddo

   if (kount.ne.kkxlmn) then
      call error('clbsch','did not count correctly',' ')
   endif

   return
   end subroutine clbsch

   real(kr) function del2(a,b,c)
   !-------------------------------------------------------------------
   ! Generate { [(a+b-c)! (a-b+c)! (-a+b+c)!] / (a+b+c+1)! }
   ! in a fashion which requires as few as possible multiplications
   ! and which is not likely to produce over- or under-flow problems
   !-------------------------------------------------------------------
   ! externals
   real(kr)::a,b,c
   ! internals
   real(kr)::d,e,f,g,h,x,y,z
   real(kr),parameter::q=.05e0_kr

   h=0
   d=a+b-c
   e=a-b+c
   f=-a+b+c

   if (d.ge.-q.and.e.ge.-q.and.f.ge.-q) then

      g= a+b+c+1
      h=1/g

      if (d.ge.e.and.d.ge.f) then
         ! Here d is the biggest of (d,e,f)
         x=d
         y=e
         z=f
      else if (e.ge.f.and.e.ge.d) then
         ! Here e is the biggest of (d,e,f)
         x=e
         y=f
         z=d
      else if (f.ge.d.and.f.ge.e) then
         ! Here f is the biggest of (d,e,f)
         x=f
         y=d
         z=e
      else
         STOP ' in Del2 in rml/mrml04.f'
      endif

      !--Now generate y!/(part of g!/x!)
      do
         x=x+1
         if (y.lt.q) exit
         h=h*(y/x)
         y=y-1
      enddo

      !--Next generate z!/(rest of g!/x!)
      do
         if (z.lt.q) exit
         h=h*(z/x)
         z=z-1
         x=x+1
      enddo

   endif
   del2=h

   return
   end function del2

   real(kr) function gll(l1,l2,l3)
   !-------------------------------------------------------------------
   ! Private routine for kclbsch.  Generates the quantity
   ! (-)^m { m! / [ (m-L1)! (m-L2)! (m-Ll)! ] }
   ! where 2m = L1 + L2 + Ll using as few as possible multiplications.
   !-------------------------------------------------------------------
   ! externals
   integer::l1,l2,l3
   ! internals
   real(kr)::h
   integer::m,mg,m1,m2,m3,j,k

   j=0
   k=0
   m=(l1+l2+l3)/2
   mg=m
   m1=m-l1
   m2=m-l2
   m3=m-l3
   h=1

   if (m1.ge.m2.and.m1.ge.m3) then
      ! Here m1 is the largest of (m1,m2,m3)
      j=m2
      k=m3

   else if (m2.ge.m3.and.m2.ge.m1) then
      ! Here m2 is the largest of (m1,m2,m3)
      j=m3
      k=m1

   else if (m3.ge.m1.and.m3.ge.m2) then
      ! Here m3 is the largest of (m1,m2,m3)
      j=m1
      k=m2

   else
      ! Should never get here
      STOP ' in Gll in rml/mrml04.f'
   endif

   do
      if (j.eq.0) exit
      h=h*m/j
      m=m-1
      j=j-1
   enddo

   do
      if (k.eq.0) exit
      h=h*m/k
      m=m-1
      k=k-1
   enddo

   m=mg/2
   if (m*2.ne.mg) h=-h

   gll=h

   return
   end function gll

   real(kr) function ww(el1,aj1,el2,aj2,cj,el)
   !-------------------------------------------------------------------
   ! Private routine for kclbsch.  Generates the following expression:
   ! sum over n of
   ! (-)**(n+L1+J1+L2+J2) * (n+1)! /
   ! [ (n-(L1+J1+j))! (n-(L2+J2+j))! (n-(L1+L2+Ll))! (n-(J1+J2+Ll))!
   !   (L1+J1+L2+J2-n)! (L1+J2+Ll+j-n)! (L2+J1+Ll+j-n)! ]
   !-------------------------------------------------------------------
   ! externals
   real(kr)::el1,aj1,el2,aj2,cj,el
   ! internals
   real(kr)::h,sign,g
   integer::i1,i2,m1,m2,k1,k2,k3,minz,maxz,n,nm,i
   integer::nn(7)

   h=0
   i1=nint(el1+aj1+cj)
   i2=nint(el2+aj2+cj)
   m1=nint(el1+el2+el)
   m2=nint(aj1+aj2+el)
   minz=max0(i1,i2,m1,m2)
   k1=nint(el1+aj1+el2+aj2)
   k2=nint(el1+aj2+cj+el)
   k3=nint(el2+aj1+cj+el)
   maxz=min0(k1,k2,k3)

   if (maxz.ge.minz) then

      n=(k1+minz)/2
      sign=1
      if (2*n.ne.(k1+minz)) sign=-1
      do n=minz,maxz
         nn(1)=n-i1
         nn(2)=n-i2
         nn(3)=n-m1
         nn(4)=n-m2
         nn(5)=k1-n
         nn(6)=k2-n
         nn(7)=k3-n
         call sortt(nn)

         g=sign*float(n+1)

         ! set g = (n+1)!/(nnmax)!
         nm=n
         do while (nm.ne.nn(1))
            g=g*float(nm)
            nm=nm-1
         enddo

         ! now divide by the other (nn !)'s
         do i=2,7
            nm=nn(i)
            if (nm.le.1) exit
            do while (nm.gt.1)
               g=g/float(nm)
               nm=nm-1
            enddo
         enddo
         h=h+g
         sign=-sign
      enddo

   endif
   ww=h

   return
   end function ww

   subroutine sortt(nn)
   !-------------------------------------------------------------------
   ! Order nn(i) from high to low.
   !-------------------------------------------------------------------
   ! externals
   integer::nn(7)
   ! internals
   integer::i,is,k,mm

   do i=1,7
      is=0
      do k=1,6
         if (nn(k).lt.nn(k+1)) then
            mm=nn(k)
            nn(k)=nn(k+1)
            nn(k+1)=mm
            is=1
         endif
      enddo
      if (is.eq.0) return
   enddo

   return
   end subroutine sortt

   subroutine babb(ier)
   !-------------------------------------------------------------------
   ! Generate energy-independent portion of partial of R
   ! with respect to U-parameters.
   !-------------------------------------------------------------------
   ! internals
   integer::kj,ipar,maxres,igr,minres,mmax,mmax2,ires,m,k,j
   integer::mmm,kk,iipar,jk,ier
   real(kr)::d1,d
   real(kr),parameter::zero=0

   do kj=1,ntriag
      do ipar=1,npar
         br(kj,ipar,ier)=0
         bi(kj,ipar,ier)=0
      enddo
   enddo

   ipar=0
   maxres=0
   do igr=1,ngroup
      minres=maxres+1
      maxres=nresg(igr,ier)
      mmax=nchan(igr,ier)
      mmax2=mmax+2

      do ires=minres,maxres
         do m=1,mmax2
            ipar=ipar+1

            if (m.eq.1) then
               ! here U-parameter is resonance energy
               par(ipar)=eres(ires,ier)
               d1=2*sqrt(abs(eres(ires,ier)))
               kj=0
               do k=1,mmax
                  do j=1,k
                     kj=kj+1
                     d=beta(kj,ires,ier)*d1
                     br(kj,ipar,ier)=d
                     bi(kj,ipar,ier)=-d-d
                  enddo
               enddo

            else if (m.eq.2) then
               ! here U-parameter is gamma-sub-gamma
               par(ipar)=gamgam(ires,ier)
               d1=gbetpr(1,ires,ier)
               kj=0
               do k=1,mmax
                  do j=1,k
                     kj=kj+1
                     d=beta(kj,ires,ier)*d1
                     d=d+d
                     br(kj,ipar,ier)=-d-d
                     bi(kj,ipar,ier)=d
                  enddo
               enddo

            else if (m.gt.2) then
               ! here U-parameters are gamma-sub-channel(m-2)
               mmm=m-2
               par(ipar)=gamma(mmm,ires,ier)
               kj=0
               do k=1,mmax
                  do j=1,k
                     kj=kj+1
                     if (j.eq.mmm) then
                        d=betapr(k,ires,ier)
                        if (k.eq.mmm) d=d+d
                        br(kj,ipar,ier)=d
                        bi(kj,ipar,ier)=d
                     endif
                  enddo
                  if (k.ne.mmax) then
                     kk=k+1
                     do j=kk,mmax
                        if (j.eq.mmm) then
                           jk=(j*(j-1))/2+k
                           d=betapr(k,ires,ier)
                           if (k.eq.mmm) d=d+d
                           br(jk,ipar,ier)=d
                           bi(jk,ipar,ier)=d
                        endif
                     enddo
                  endif
               enddo
            endif
         enddo
      enddo
   enddo

   do iipar=1,npar
      kj=0
      do k=1,mchan
         do j=1,k
            kj=kj+1
            if (br(kj,iipar,ier).ne.zero) br(kj,iipar,ier)=2*br(kj,iipar,ier)
            if (bi(kj,iipar,ier).ne.zero) bi(kj,iipar,ier)=2*bi(kj,iipar,ier)
         enddo
      enddo
   enddo

   return
   end subroutine babb

   subroutine abpart(ier)
   !-------------------------------------------------------------------
   ! Generate alphar and alphai => for cross section
   ! and upr and upi = Energy-dependent pieces of pr & pi.
   ! Also generate pr and pi = partial of R wrt U-parameters.
   !-------------------------------------------------------------------
   ! externals
   integer::ier
   ! internals
   integer::ipar,maxres,igr,minres,n2,n,m,i,ij,j,k
   real(kr)::aa
   real(kr),parameter::zero=0

   !--Generate alphar =  (del E) / ( (del E)**2+(Gamgam/2)**2 )
   !--and alphai = Gamgam/2 / ( Ditto )

   do n=1,nresm(ier)
      xden(n,ier)=0
      alphar(n,ier)=0
      alphai(n,ier)=0
      difen(n,ier)=eres(n,ier)-su
      aa=difen(n,ier)**2+gbetpr(3,n,ier)
      xden(n,ier)=1/aa
      alphar(n,ier)=difen(n,ier)*xden(n,ier)
      alphai(n,ier)=gbetpr(2,n,ier)*xden(n,ier)
   enddo

   if (Want_Partial_Derivs) then

      !--Generate Upr and Upi = Energy-Dependent part of partial derivs

      ipar=0
      maxres=0
      do igr=1,ngroupm(ier)
         minres=maxres+1
         maxres=nresg(igr,ier)
         n2=nchan(igr,ier)+2
         if (maxres.ge.minres) then
            do n=minres,maxres
               do m=1,n2
                  ipar=ipar+1
                  upr(ipar,ier)=alphar(n,ier)
                  upi(ipar,ier)=alphai(n,ier)
                  if (m.eq.1) then
                     ! variable is resonance-energy
                     upi(ipar,ier)=upr(ipar,ier)*upi(ipar,ier)
                     upr(ipar,ier)=xden(n,ier)-2*upr(ipar,ier)*upr(ipar,ier)
                  else if (m.eq.2) then
                     ! variable is capture width
                     upr(ipar,ier)=upr(ipar,ier)*upi(ipar,ier)
                     upi(ipar,ier)=xden(n,ier)-2*upi(ipar,ier)*upi(ipar,ier)
                  endif
               enddo
            enddo
         endif
      enddo

      !--Multiply by BR and BI to give partial of R wrt U-parameters

      do i=1,ntriag
         do j=1,npar
            pr(i,j,ier)=0
            pii(i,j,ier)=0
         enddo
      enddo

      do k=1,npar
         if (upr(k,ier).ne.zero.or.upi(k,ier).ne.zero) then
            ij=0
            do i=1,mchan
               do j=1,I
                  ij=ij+1
                  if (br(ij,k,ier).ne.zero) then
                     pr(ij,k,ier)=br(ij,k,ier)*upr(k,ier)
                  endif
                  if (bi(ij,k,ier).ne.zero) then
                     pii(ij,k,ier)=bi(ij,k,ier)*upi(k,ier)
                  endif
               enddo
            enddo
         endif
      enddo

   endif

   return
   end subroutine abpart

   subroutine crosss(ier)
   !-------------------------------------------------------------------
   ! Form the cross sections Sigma(Ksigma) and the
   ! ( partial derivatives of the cross section with respect to
   ! the resonance parameters ) = Dsigma(Ksigma,Ipar)
   !-------------------------------------------------------------------
   use physics, only : pi ! get pi
   ! externals
   integer::ier
   ! internals
   integer::i,n,nn2,nn,minr,nchann,nentnn,nextnn,nnnn,npr,mx,mplus
   integer::kstart,jstart,maxr,kount,lrmat,ip,mm,m,needxq,j,k,l
   real(kr)::u
   real(kr),parameter::fourpi=4*pi/100
   real(kr),parameter::zero=0

   !--Initialize

   do i=1,nppm(ier)
      sigmas(i,ier)=0
   enddo

   if (Want_Angular_Dist) then
      do i=1,lllmax
         do j=1,nppm(ier)
            Coef_Leg(i,j,ier)=0
         enddo
      enddo
   endif

   if (Want_Partial_Derivs) then
      n=npp*npar
      if (n.ne.0) then
         do i=1,npp
            do j=1,npar
               dsigma(i,j,ier)=0
            enddo
         enddo
      endif
      if (Want_Angular_Dist) then
         n=lllmax*npp*npar
         if (n.ne.0) then
            do i=1,lllmax
               do j=1,npp
                  do k=1,npar
                     D_Coef_Leg(i,j,k,ier)=0
                  enddo
               enddo
            enddo
         endif
      endif
   endif

   if (Want_Angular_Dist) then
      do i=1,2
         do j=1,mchan
            do k=1,mchan
               do l=1,npp
                  do m=1,ngroupm(ier)
                     crssx(i,j,k,l,m,ier)=0
                  enddo
               enddo
            enddo
         enddo
      enddo
   endif

   krext=nrext
   if (krext.eq.0) krext=1

   !--do loop over spin-parity groups

   kstart=0
   jstart=nparr
   maxr=0
   kount=1
   do n=1,ngroupm(ier)

      !--Initialize for this group
      do i=1,nppm(ier)
         crss(i,ier)=0
      enddo
      if (Want_Angular_Dist) then
         do i=1,2
            do j=1,ntriag
               cscs(i,j,ier)=0
            enddo
         enddo
      endif
      if (Want_Partial_Derivs) then
         do i=1,nppm(ier)
            do j=1,npar
               deriv(i,j,ier)=0
            enddo
         enddo
      endif

      nn2=nchan(n,ier)*(nchan(n,ier)+1)
      nn=nn2/2

      minr=maxr+1
      maxr=nresg(n,ier)
      nchann=nchan(n,ier)
      nentnn=nent(n,ier)
      nextnn=next(n,ier)
      nnnn=n
      npr=(nchann+2)*(maxr-minr+1)

      lrmat=0

      !--set R-Matrix and other necessary arrays
      call setr(n,lrmat,minr,maxr,nentnn,nchann,ier)

      if (lrmat.eq.1) then
         !--Calculate Xq & Xxxx matrices if trivial
         call zeror(nchann,ier)
      else

         !--invert Ymat
         call yinvrs(nchann,ier)

         !--generate XQ & Xxxx matrices
         call setxqx(nchann,ier)
      endif

      needxq=1
      if (npr.ne.0.and.lrmat.eq.0) needxq=0

      !--Generate cross section pieces
      call sectio(n,nentnn,nextnn,nchann,ier)

      !--Done calculating contribution for this group; ergo, add to totals
      do ip=1,nppm(ier)
         sigmas(ip,ier)=crss(ip,ier)+sigmas(ip,ier)
      enddo

      if (Want_Angular_Dist) then
         ! Set Coef_Leg, coefficients of Legendre polynomials
         call setleg(nnnn,kount,ier)
      endif

      if (Want_Partial_Derivs) then

         !--generate Q = partial derivative of Xxxx wrt R
         if (needxq.eq.0) call setqri(nchann,nn,ier)

         !--generate T = partial of cross sections with respect to R
         !--T = [ partial (sigma) wrt X ] * Q
         if (needxq.eq.0) call settri(n,nentnn,nextnn,nchann,nn,ier)

         !--find derivatives of cross sections wrt resonance parameters
         if (lrmat.eq.0) call derres(n,nentnn,nchann,minr,maxr,&
           kstart,npr,nn,ier)

         !--find deriv of cross sections with respect to R-ext pars
         if (nrext.ne.0) call derext(n,nchann,jstart,nnnn,ier)

         !--Done calculating derivs for this group; ergo, add to totals
         do ip=1,nppm(ier)
            do mm=1,npr
               m=kstart+mm
               dsigma(ip,m,ier)=dsigma(ip,m,ier)+deriv(ip,m,ier)
            enddo
         enddo

         kstart=kstart+npr
         jstart=jstart+nchan(n,ier)
      endif

   enddo

   !--Normalize properly
   do ip=1,nppm(ier)
      sigmas(ip,ier)=sigmas(ip,ier)*fourpi/su
   enddo
   if (Want_Angular_Dist) then
      do ip=1,nppm(ier)
         if (ip.ne.2) then
            do l=1,lllmax
               Coef_Leg(l,ip,ier)=Coef_Leg(l,ip,ier)/su
            enddo
         endif
      enddo
   endif

   if (Want_Partial_Derivs) then
      !--Normalize by fourpi/su as for cross section.
      !--Also, divide by Uuuu => convert partial derivatives to be with
      !--respect to Gamma instead of gamma (and E instead of sqrt(E)).
      !--The duuu term takes care of the E-dependence of the penetrabilities.
      do m=1,npar
         u=uuuu(m,ier)
         if (u.ne.zero) then
            u=fourpi/su/u
            if (.not.Want_Partial_U) then
               if (iduu(m,ier).eq.1) then
                  do mx=1,mchan
                     mplus=m+mx+1
                     if (mplus.gt.npar) exit
                     if (iduu(mplus,ier).gt.0) exit
                     if (duuu(mplus,ier).ne.zero) then
                        do ip=1,nppm(ier)
                           dsigma(ip,m,ier)=dsigma(ip,m,ier)&
                             -duuu(mplus,ier)*dsigma(ip,mplus,ier)
                        enddo
                     endif
                  enddo
               endif
            endif
            do ip=1,nppm(ier)
               dsigma(ip,m,ier)=dsigma(ip,m,ier)*u
            enddo
         else
            do ip=1,nppm(ier)
               dsigma(ip,m,ier)=0
            enddo
         endif
      enddo
   endif

   return
   end subroutine crosss

   subroutine setr(n,lrmat,min,max,nent,nchan,ier)
   !-------------------------------------------------------------------
   ! Generate the R-matrix and other arrays.
   !-------------------------------------------------------------------
   ! externals
   integer::n,lrmat,min,max,nent,nchan,ier
   ! internals
   integer::nnntot,i,kl,k,l,ires,ii,ipx,iffy,kk,kx,kkx,j,ji
   integer::lsp,jdopha
   real(kr)::aloge,ex,rho,rhof,eta,hr,hi,p,dp,ds
   real(kr)::sinphx,cosphx,dphix
   real(kr)::hrx,hix,px,dpx,dsx
   real(kr),parameter::zero=0
   real(kr),parameter::one=1
   real(kr),parameter::tiny=1.e-8_kr

   nnntot=nchan+1
   do i=1,nchan
      nnntot=nnntot-1
      if (su.ge.zero.and.su.le.echan(nnntot,n,ier)) then
         nchan=nnntot-1
      else
         exit
      endif
   enddo

   !--initialize R-matrix
   kl=0
   do k=1,nchan
      if (nrext.ne.0) aloge=log((parext(2,k,n,ier)-su)/(su-parext(1,k,n,ier)))
      do l=1,k
         kl=kl+1
         rmat(1,kl,ier)=0
         rmat(2,kl,ier)=0
         if (l.eq.k.and.nrext.ne.0) then
            rmat(1,kl,ier)=parext(3,k,n,ier)+parext(4,k,n,ier)*su&
             -parext(5,k,n,ier)*aloge+parext(7,k,n,ier)*su**2&
             -parext(6,k,n,ier)*(parext(2,k,n,ier)-parext(1,k,n,ier))&
             -parext(6,k,n,ier)*aloge*su
         endif
      enddo
   enddo

   if (max.ge.min.and.min.gt.0) then
      do ires=min,max
         kl=0
         do k=1,nchan
            do l=1,k
               kl=kl+1
               if (su.gt.echan(k,n,ier).and.su.gt.echan(l,n,ier).and.&
                 beta(kl,ires,ier).ne.zero) then
                  rmat(1,kl,ier)=rmat(1,kl,ier)&
                    +alphar(ires,ier)*beta(kl,ires,ier)
                  rmat(2,kl,ier)=rmat(2,kl,ier)&
                    +alphai(ires,ier)*beta(kl,ires,ier)
               endif
            enddo
         enddo
      enddo
   endif

   !--check if rmat is zero; if so, set lrmat=1.
   kl=0
   do k=1,nchan
      do l=1,k
         kl=kl+1
         if (rmat(1,kl,ier).ne.zero) go to 20
         if (rmat(2,kl,ier).ne.zero) go to 20
      enddo
   enddo
   lrmat=1
  20 continue

   !--generate rootp, ph, and (H+Rmat) matrices
   kl=0
   do k=1,nchan
      do l=1,k
         kl=kl+1
         ymat(1,kl,ier)=-rmat(1,kl,ier)
         ymat(2,kl,ier)=-rmat(2,kl,ier)
      enddo
   enddo

   do i=1,nchan
      psmall(i,ier)=0
   enddo

   ii=0
   do i=1,nchan
      ipx=ipp(i,n,ier)
      ii=ii+i
      rootp(i,ier)=1
      elinvr(i,ier)=0
      elinvi(i,ier)=-1
      dpdr(i,ier)=0
      dsdr(i,ier)=0
      iffy=0
      if (su.gt.echan(i,n,ier)) then
         if (lpent(ipx,ier).le.0) then

            !--here penetrability is not calculated
            ymat(2,ii,ier)=ymat(2,ii,ier)-1

         else

            !)--and here it is
            lsp=lspin(i,n,ier)
            ex=sqrt(su-echan(i,n,ier))
            rho=zkte(i,n,ier)*ex
            rhof=zkfe(i,n,ier)*ex

            if (zeta(i,n,ier).eq.zero) then

               !--calculate non-Coulomb phases and penetrabilities
               call sinsix(rhof,lsp,sinsqr(i,ier),sin2ph(i,ier),dphi(i,ier),&
                 sinphi(i,ier),cosphi(i,ier))
               call pgh(rho,lsp,bound(i,n,ier),hr,hi,p,dp,ds,&
                 ishift(ipx,ier),iffy)
               !--hr and hi are real and imag parts of 1/(S-B+iP)
               !--except when S-B+iP=zero, in which case iffy=1

            else

               !--calculate Coulomb penetrabilities
               eta=zeta(i,n,ier)/ex
               if (rho.eq.rhof) then
                  jdopha=1
                  call pghcou(rho,lsp,bound(i,n,ier),hr,hi,p,dp,ds,&
                    ishift(ipx,ier),iffy,eta,sinphi(i,ier),&
                    cosphi(i,ier),dphi(i,ier),jdopha)
               else
                  jdopha=0
                  call pghcou(rho,lsp,bound(i,n,ier),hr,hi,p,dp,ds,&
                    ishift(ipx,ier),&
                    iffy,eta,sinphx,cosphx,dphix,jdopha)
                  jdopha=1
                  call pghcou(rhof,lsp,bound(i,n,ier),hrx,hix,px,dpx,dsx,&
                    ishift(ipx,ier),iffy,eta,sinphi(i,ier),&
                    cosphi(i,ier),dphi(i,ier),jdopha)
               endif
               sinsqr(i,ier)=sinphi(i,ier)**2            ! sin^2 of phase shift
               sin2ph(i,ier)=2*sinphi(i,ier)*cosphi(i,ier) ! sin(2*phase shift)
            endif

            !--for angular distributions
            if (Want_Angular_Dist) then
               kk=(i*(i-1))/2
               if (i.gt.1) then
                  do kx=1,i-1
                     kkx=kk+kx
                     cscs(1,kkx,ier)=cosphi(kx,ier)*cosphi(i,ier)&
                       -sinphi(kx,ier)*sinphi(i,ier)
                     cscs(2,kkx,ier)=cosphi(kx,ier)*sinphi(i,ier)&
                       +sinphi(kx,ier)*cosphi(i,ier)
                  enddo
               endif
            endif

            rootp(i,ier)=sqrt(p)
            dpdr(i,ier)=dp
            dsdr(i,ier)=ds
            if (iffy.eq.0.and..not.(ishift(ipx,ier).le.0.and.&
              (1-p*rmat(2,ii,ier).eq.one.or.p.lt.tiny))) then
               elinvr(i,ier)=hr
               elinvi(i,ier)=hi
               ymat(1,ii,ier)=hr+ymat(1,ii,ier)
               ymat(2,ii,ier)=hi+ymat(2,ii,ier)
            else
               !--here penetrability is very small but non-zero
               psmall(i,ier)=rootp(i,ier)
               ymat(1,ii,ier)=p*ymat(1,ii,ier)
               ymat(2,ii,ier)=p*ymat(2,ii,ier)-one
               rmat(1,ii,ier)=p*rmat(1,ii,ier)
               rmat(2,ii,ier)=p*rmat(2,ii,ier)
               if (nchan.gt.1) then
                  if (i.gt.1) then
                     do j=1,i-1
                        ji=(i*(i-1))/2+j
                        ymat(1,ji,ier)=rootp(i,ier)*ymat(1,ji,ier)
                        ymat(2,ji,ier)=rootp(i,ier)*ymat(2,ji,ier)
                        rmat(1,ji,ier)=rootp(i,ier)*rmat(1,ji,ier)
                        rmat(2,ji,ier)=rootp(i,ier)*rmat(2,ji,ier)
                     enddo
                  endif
                  if (i.lt.nchan) then
                     do j=i+1,nchan
                        ji=(j*(j-1))/2+i
                        ymat(1,ji,ier)=rootp(i,ier)*ymat(1,ji,ier)
                        ymat(2,ji,ier)=rootp(i,ier)*ymat(2,ji,ier)
                        rmat(1,ji,ier)=rootp(i,ier)*rmat(1,ji,ier)
                        rmat(2,ji,ier)=rootp(i,ier)*rmat(2,ji,ier)
                     enddo
                  endif
               endif
               rootp(i,ier)=1
               elinvr(i,ier)=0
               elinvi(i,ier)=-1
               dpdr(i,ier)=0
               dsdr(i,ier)=0
            endif
         endif
      else
         rootp(i,ier)=0
         elinvr(i,ier)=1
         elinvi(i,ier)=0
         dpdr(i,ier)=0
         dsdr(i,ier)=0
      endif
   enddo

   !--check if one channel is irrelevant; if so, set to unity.
   if (lrmat.ne.1) then
      kl=0
      do k=1,nchan
         kl=kl+k
         if (ymat(1,kl,ier).eq.zero.and.ymat(2,kl,ier).eq.zero) then
            ymat(1,kl,ier)=1
         endif
      enddo
   endif

   return
   end subroutine setr

   subroutine sinsix(rho,lspin,sinsqr,sin2ph,dphi,sinphi,cosphi)
   !-------------------------------------------------------------------
   ! Generate sinsqr = [sin(Phi)]**2 and sin2ph = sin(2*Phi)
   ! Also generate [partial derivative of Phi wrt Rho] = Dphi
   ! Also generate cos(Phi) and sin(Phi), which requires
   ! getting the appropriate signs
   !-------------------------------------------------------------------
   ! externals
   real(kr)::rho,sinsqr,sin2ph,dphi,sinphi,cosphi
   integer::lspin
   ! internals
   real(kr)::a,c,s,g,d,dd,a2,y,f,a4,p,dp,ss,ds,x
   integer::l
   real(kr),parameter::zero=0

   a=rho
   l=lspin
   c=cos(a)
   s=sin(a)

   if (l.eq.0) then
      ! Phi0
      sinsqr=s*s
      sin2ph=2*c*s
      cosphi=c
      sinphi=s
      dphi=1
      return

   else if (l.eq.1) then
      ! Phi1
      x=a
      g=1+x**2
      d=(s-c*x)/g
      sinsqr=(s-c*x)*d
      sin2ph=2*(c+s*x)*d
      dphi=x*x/g
      dd=sqrt(g)
      cosphi=(c+s*x)/dd
      sinphi=(s-c*x)/dd
      return

   else if (l.eq.2) then
      ! Phi2
      a2=a*a
      x=3
      y=3-a2
      f=x*y+2*a2*x
      x=x*a

   else if (L.eq.3) then
      ! Phi3
      a2=a*a
      x=15-a2
      y=15-a2*6
      f=x*Y-2*a2*y+12*a2*x
      x=x*a

   else if (L.eq.4) then
      ! Phi4
      a2=a*a
      a4=a2*a2
      x=105-10*a2
      y=105-45*a2+a4
      f=x*Y-20*a2*y+a2*x*(90-4*a2)
      x=x*a

   else
      call genpsf(a,l,p,dp,ss,ds,x,f,y)
   endif

   g=y**2+x**2
   d=(s*y-c*x)/g

   sinsqr=(s*y-c*x)*d
    ! = [ sin (hard-sphere phase shift) ] **2
   sin2ph=2*(c*y+s*x)*d
    ! = sin ( 2 * hard-sphere phase shift )
   dphi=1-f/g
    ! = derivative of Phi wrt Rho

   dd=sqrt(g)
   cosphi=(c*y+s*x)/dd
   sinphi=(s*y-c*x)/dd

   return
   end subroutine sinsix

   subroutine pgh(rho,l,bound,g,h,p,dp,ds,ishift,iffy)
   !-------------------------------------------------------------------
   ! Forms the real and imaginary parts g+ih
   ! of the inverse of S-B+iP, and also sets P, dP, and dS.
   ! Note -- if S-B+iP=0, then set G=0,H=0,P=0,Dp=0, and iffy=1
   !-------------------------------------------------------------------
   ! externals
   real(kr)::rho,bound,g,h,p,dp,ds
   integer::l,ishift,iffy
   ! internals
   real(kr)::r,gg,r2,r4,r6,r8,hh,d,s,x,db,y
   real(kr),parameter::zero=0

   iffy=0
   r=rho
   gg=0
   dp=0
   ds=0

   if (l.eq.0) then
      hh=r
      dp=1
      if (ishift.gt.0) gg=-bound

   else if (l.eq.1) then
      r2=r*r
      d=1+r2
      hh=r*r2/d
      dp=(r2/d)*(3+r2)/d
      if (ishift.gt.0) then
         gg=-1/d-bound
         ds=2*(r/d)/d
      endif

   else if (l.eq.2) then
      r2=r*r
      r4=r2*r2
      d=9+r2*(3+r2)
      hh=r*r4/d
      dp=r2*(r2/d)*(45+r2*(9+r2))/d
      if (ishift.gt.0) then
         ds=(18+3*r2)/d
         gg=-ds-bound
         ds=-6*r/d+ds*4*r*(1.5e0_kr+r2)/d
      endif

   else if (l.eq.3) then
      r2=r*r
      r4=r2*r2
      r6=r4*r2
      d=225+r2*(45+r2*(6+r2))
      hh=r*r6/d
      dp=(r6/d)*(1575+(225+r2*(18+r2))*r2)/d
      if (ishift.gt.0) then
         ds=(675+r2*(90+r2*6))/d
         gg=-ds-bound
         ds=-24*(7.5e0_kr+r2)*r/d +ds*6*r*(15+r2*(4+r2))/d
      endif

   else if (l.eq.4) then
      r2=r*r
      r4=r2*r2
      r6=r4*r2
      r8=r4*r4
      d=11025+r2*(1575+r2*(135+r2*(10+r2)))
      hh=r*r8/d
      dp=(r8/d)*(9-(1575+r2*(270+r2*(30+r2*4)))*2*(r2/d))
      if (ishift.gt.0) then
         ds=(44100+r2*(4725+r2*(270+10*r2)))/d
         gg=-ds-bound
         ds=-60*r*(157.5e0_kr+r2*(18+r2))/d+&
               ds*r*(3150+r2*(540+r2*(60+r2*8)))/d
      endif

   else
      call genpsf(rho,l,hh,dp,s,ds,x,db,y)
      if (ishift.gt.0) gg=S - bound
   endif

   !--Here hh = P(L), Dp = partial(P) wrt Rho,
   !--gg = S(L) - Bound, Ds = partial(S) wrt Rho

   if (hh.le.1.0e-35_kr) hh=zero

   if (gg.eq.zero.and.hh.eq.zero) then
      g=0
      h=0
      p=0
      ds=0
      iffy=1

   else if (gg.eq.zero) then
      g=0
      h=-1/hh
      p=hh

   else if (hh.eq.zero) then
      g=1/gg
      h=0
      p=0

   else if (hh+gg.ne.hh) then
      if (hh+gg.ne.gg) then
         d=hh**2+gg**2
         g=gg/d
         h=-hh/d
         p=hh
      else
         g=1/gg
         h=-(hh/gg)/gg
         p=hh
      endif

   else
      g=(gg/hh)/hh
      h=-1/hh
      p=hh

   endif

   return
   end subroutine pgh

   subroutine genpsf(rho,ll,pp,dp,ss,ds,qq,dq,yy)
   !-------------------------------------------------------------------
   ! externals
   real(kr)::rho,pp,dp,ss,ds,qq,dq,yy
   integer::ll
   ! internals
   real(kr)::a,a2,a4,a8,yyy,qqq,dqq,d,ppp,dpp,sss,dss
   real(kr)::el,els,x,dx,yyx,qqx,dqx,ppx,c,dpx,ssx,dsx
   integer::lwant,l

   lwant=ll
   a=rho
   a2=a*a
   a4=a2*a2
   a8=a4*a4

   l=4
   yyy=105-45*a2+a4
   qqq=105-10*a2
   dqq=(qqq-20*a2)*yyy+a2*qqq*(90-4*a2)
   qqq=qqq*a
   !--phi = rho - datan(qqq/yyy)
   !--qqq = [ partial(qqq/yyy) wrt rho] * yyy^2; for L=4

   d=11025+a2*(1575+a2*(135+a2*(10+a2)))
   ppp=a*a8/d
   dpp=(a8/d)*(9-(1575+a2*(270+a2*(30+a2*4)))*2*(a2/d))
   !--ppp = P(L), dpp = partial(P) wrt rho; for L=4

   sss=-dss
   dss=-60*a*(157.5d0+a2*(18+a2))/d+dss*a*(3150+a2*(540+a2*(60+a2*8)))/d
   !--sss = S(L), dss = partial(S) wrt rho; for L=4

  10 continue
   !--Now we have values for L; find values for L+1
   l=l+1
   el=l
   els=el-sss

   x=ppp/els
   dx=(dpp+x*dss)/els
   yyx=yyy-qqq*X
   qqx=qqq+x*yyy
   d=qqx*(qqq*dx*yyy+dqq*x)
   dqx=(dqq+dx*yyy**2)*yyx+d

   d=els**2+ppp**2
   ppx=a2*ppp/d
   c=1+a*(els*dss-ppp*dpp)/d
   dpx=2*ppp*c+a*dpp
   dpx=a*dpx/d

   ssx=a2*els/d-el
   dsx=2*els*c-a*dss
   dsx=a*dsx/d

   if (l.ne.lwant) then
      ppp=ppx
      dpp=dpx
      sss=ssx
      dss=dsx
      qqq=qqx
      dqq=dqx
      yyy=yyx
      go to 10

   else

      pp=ppx
      dp=dpx
      ss=ssx
      ds=dsx
      qq=qqx
      dq=dqx
      yy=yyx
      return
   endif

   return
   end subroutine genpsf

   real(kr) function pf(rho,l)
   !-------------------------------------------------------------------
   ! pf forms the penetration factor with no derivatives
   !-------------------------------------------------------------------
   ! externals
   real(kr)::rho
   integer::l
   ! internals
   real(kr)::a,a2,a4,a6,a8,d,p,dp,s,ds,qq,db,yy

   a=rho

   if (l.eq.0) then
      ! P0
      pf=a

   else if (l.eq.1) then
      ! P1
      a2=a*a
      pf=a*a2/(1+a2)

   else if (l.eq.2) then
      ! P2
      a2=a*a
      a4=a2*a2
      pf=a*a4/(9+a2*(3+a2))

   else if (l.eq.3) then
      ! P3
      a2=a*a
      a4=a2*a2
      a6=a2*a4
      d=225+a2*(45+a2*(6+a2))
      pf=a*a6/d

   else if (l.eq.4) then
      ! P4
      a2=a*a
      a4=a2*a2
      a8=a4*a4
      d=11025+a2*(1575+a2*(135+a2*(10+a2)))
      pf=a*a8/d

   else
      ! PL for L > 4
      call genpsf(rho,l,p,dp,s,ds,qq,db,yy)
      pf=p

   endif

   return
   end function pf

   subroutine pghcou(rho,l,bound,g,h,p,dp,ds,ishift,iffy,&
     eta,sinphi,cosphi,dphi,jdopha)
   !-------------------------------------------------------------------
   ! Coulomb version of sub pgh. Forms the real and imaginary
   ! parts G+iH of the inverse of S-B+iP and also sets P, dP, & dS.
   !
   !     ishift = ishift(I,nnnn)
   !     S = shift factor             P = penetrability
   !     Ds = dS/dRho                 Dp = dP/dRho
   !
   ! Note -- If S-B+iP=0, then set G=0,H=0,P=0,Dp=0, and iffy=1
   !-------------------------------------------------------------------
   ! externals
   real(kr)::rho,bound,g,h,p,dp,ds,eta,sinphi,cosphi,dphi
   integer::l,ishift,jdopha,iffy
   ! internals
   real(kr)::gg,hh,s,pent,shift,dpent,dshift,d
   integer::ifail,jdoder
   real(kr),parameter::zero=0

   iffy=0
   gg=0
   dp=0
   ds=0

   jdoder=1
   call pspcou(rho,l,eta,ishift,jdopha,jdoder,ifail,&
     pent,shift,dpent,dshift,sinphi,cosphi,dphi)

   hh=pent
   dp=dpent
   s=shift
   ds=dshift

   if (ishift.gt.0) then
      gg=s-bound
   else
      ds=0
   endif

   !--Here hh = P(L), Dp = partial(P) wrt Rho,
   !--gg = S(L) - Bound, Ds = partial(S) wrt Rho

   if (hh.le.1.0e-35_kr) hh=0

   p=0
   h=0
   g=0
   if (gg.eq.zero.and.hh.eq.zero) then
      dp=0
      iffy=1
      return
   endif
   if (hh.eq.zero) then
      g=1/gg
      return
   endif

   p=hh
   if (gg.eq.zero) then
      h=-1/hh
      return
   endif

   if (hh+gg.eq.hh) then
      g=(gg/hh)/hh
      h=-1/hh
      return
   endif

   if (hh+gg.eq.gg) then
      g=1/gg
      h=-(hh/gg)/gg
      return
   endif

   d=hh**2+gg**2
   g=gg/d
   h=-hh/d

   return
   end subroutine pghcou

   subroutine pspcou(rho,lll,eta,ishift,jdopha,jdoder,ifail,&
     pencoul,shiftcoul,dpencoul,dshiftcoul,sinphi,cosphi,dphi)
   !-------------------------------------------------------------------
   ! Compute Coulomb phase shifts
   !     Ishift = (-1,1)     says (don't, do) compute shift factor
   !     Jdopha = (0,1)     says (don't, do) compute phase shift
   !     Jdoder = (0,1)     says (don't, do) compute derivatives
   !-------------------------------------------------------------------
   use mainio ! provides nsyso,nsyse
   use util   ! provides error
   ! externals
   real(kr)::rho,eta,pencoul,shiftcoul,dpencoul,dshiftcoul
   real(kr)::sinphi,cosphi,dphi
   integer::lll,ishift,jdopha,jdoder,ifail
   ! internals
   integer::l,llmin,llmax,i,iexp
   real(kr)::rho2,eta2,p,s,dp,sp,cp,dphix
   real(kr),dimension(100)::f,fpr,g,gpr,dummy

   ifail=0
   pencoul=0
   shiftcoul=0
   dpencoul=0
   dshiftcoul=0
   sinphi=0
   cosphi=1
   dphi=0

   if (lll.lt.0) ifail=10
   if (ifail.gt.0) return

   l=lll+1

   llmin=0
   llmax=lll+2
   if (llmax.gt.100) call error('pspcou','llmax larger than 100',' ')

   if (rho.lt.1.02e0_kr) then

      call coulx(eta,rho,lll,llmax,f,fpr,g,gpr,dummy,pencoul,shiftcoul,&
        dpencoul,sinphi,cosphi,dphi,jdopha,jdoder,ishift)

      ! derivatives of shiftcoul:
      if (ishift.gt.0.and.jdoder.gt.0) then
         rho2=1.01e0_kr*rho
         eta2=eta*rho/rho2
         call coulx(eta2,rho2,lll,llmax,f,fpr,g,gpr,dummy,&
           p,s,dp,sp,cp,dphix,jdopha,jdoder,ishift)
         dshiftcoul=(s-shiftcoul)/(rho2-rho)
      endif

   else

      call coulfg(rho,eta,lll,llmax,f,g,fpr,gpr,pencoul,&
        shiftcoul,dpencoul,sinphi,cosphi,dphi,jdopha,jdoder,&
        ishift,ifail,iexp)

      if (ifail.ne.0) then
         if (ishift.gt.0.and.jdoder.gt.0) then
            !--derivatives of Shiftcoul:
            rho2=1.01e0_kr*rho
            eta2=eta*rho/rho2
            call coulfg(rho2,eta2,lll,llmax,f,g,fpr,gpr,&
              p,s,dp,sp,cp,dphix,jdopha,jdoder,ishift,ifail,iexp)
            if (ifail.gt.0) then
               write(nsyse,'('' Ifail,Rho,Eta,L,F,G='',i3,2f8.5,i5,1p,&
                 &2e12.4)') ifail,rho,eta,llmin,f(l),g(l)
               write(nsyso,'('' Ifail,Rho,Eta,L,F,G='',i3,2f8.5,i5,1p,&
                 &2e12.4)') ifail,rho,eta,llmin,f(l),g(l)
               return
            endif
            dshiftcoul=(s-shiftcoul)/(rho2-rho)
         endif

      else if (ifail.gt.0) then

         write(nsyse,'('' Ifail,Rho,Eta,L,F,G='',i3,2f8.5,i5,1p,2e12.4)')&
           ifail,rho,eta,llmin,f(1),g(1)
         write(nsyso,'('' Ifail,Rho,Eta,L,F,G='',i3,2f8.5,i5,1p,2e12.4)')&
           ifail,rho,eta,llmin,f(1),g(1)
         if (l.gt.1) write(nsyse,'(45x,1p,2e12.4)') (f(i),g(i),i=2,llmax)
         if (l.gt.1) write(nsyso,'(45x,1p,2e12.4)') (f(i),g(i),i=2,llmax)
         call coulx(eta,rho,lll,llmax,f,fpr,g,gpr,dummy,pencoul,&
           shiftcoul,dpencoul,sinphi,cosphi,dphi,jdopha,jdoder,ishift)
         if (ishift.gt.0.and.jdoder.gt.0) then
            rho2=1.01e0_kr*rho
            eta2=eta*rho/rho2
            call coulx(eta2,rho2,lll,llmax,f,fpr,g,gpr,&
              dummy,p,s,dp,sp,cp,dphix,jdopha,jdoder,ishift)
            dshiftcoul=(s-shiftcoul)/(rho2-rho)
         endif
      endif
   endif

   return
   end subroutine pspcou

   subroutine coulfg(xx,eta1,lll,llmax,fc,gc,fcp,gcp,pcoul,scoul,&
     dpcoul,sinphi,cosphi,dphi,jdopha,jdoder,ishift,ifail,iexp)
   !-------------------------------------------------------------------
   ! Revised Coulomb wavefunction routine using Steed's method
   ! Returns F,G,F',G', for real Xx.gt.0, real Eta1 (including 0),
   ! and real lamda.ge.0 for integer-spaced lambda value, thus
   ! giving positive-energy solutions to the Coulomb Schrodinger
   ! equations, the Klein-Gordan equations, and suitable forms of
   ! the Dirac equation. Also sperical + cylindrical Bessel equations.
   ! For a range of lambda values (0 to Llmax) must be an integer,
   ! [which is automatic with NML re-write, which assumes Llmin=0]
   ! Starting array element is M1 = 1.  Precision: results to within
   ! 2-3 decimals of 'machine accuracy' in oscillating region
   ! x.ge.Eta1+sqrt(Eta1**2+Xlm(Xlm+1)).
   !-------------------------------------------------------------------
   use mainio ! provides nsyso,nsyse
   ! externals
   real(kr)::xx,eta1,pcoul,scoul,dpcoul,sinphi,cosphi,dphi
   real(kr)::fc(*),gc(*),fcp(*),gcp(*)
   integer::lll,llmax,jdopha,jdoder,ishift,ifail,iexp
   ! internals
   logical xlturn
   integer::l1,lp,l
   real(kr)::eta,gjwkb,paccq,acc,acc4,acch,x,xlm,e2mm1,xll
   real(kr)::xi,pk,px,ek,f,pk1,d,df,p,tk,fcl,fpl,xl,rl,el
   real(kr)::sl,fcl1,fjwkb,w,gam,q,ta,wi,ar,ai,br,bi,dr,di,dp,dq
   real(kr)::c,a,b,fcm,gcl,gpl,gcl1,asq,aaa,sss,dpacc
   real(kr),parameter::zero=0
   real(kr),parameter::one=1
   real(kr),parameter::two=2
   real(kr),parameter::ten=10
   real(kr),parameter::ten2=100
   real(kr),parameter::half=0.5e0_kr
   real(kr),parameter::accur=1.e-16_kr
   real(kr),parameter::tm30=1.e-30_kr
   real(kr),parameter::abort=2.e4_kr

   ifail=0
   iexp=1
   eta=eta1
   gjwkb=0
   paccq=1
   acc=accur
   acc4=acc*ten2*ten2/ten
   acch=sqrt(acc)
   if (xx.le.acch) go to 100
   x=xx
   xlm=0

   e2mm1=eta*eta+xlm*xlm+xlm
      !  =   Eta^2+LL(LL+1)

   xlturn=x*(x-2*eta).lt.xlm*xlm+xlm

   xll=llmax
   ! xll is max lambda value, or 0.5 smaller for J,Y Besses
   ! determine starting array element (1) from 0
   l1=llmax+1

   !--evaluate CF1 = F = FPRIME(Xl,Eta,X)/F(Xl,Eta,X)

   xi=1/x
   fcl=1
   pk=xll+1
   Px=pk+abort
  10 continue
   ek=eta/pk
   f=(ek+pk*xi)*fcl+(fcl-1)*xi
   pk1= pk+1
   if (abs(eta*x+pk*pk1).le.acc) then
      ! test ensures B1.ne.Zero for negative Eta; fixup is exact.
      fcl=(1+ek*ek)/(1+(eta/pk1)**2)
      pk=2+pk
      go to 10
   endif
   d=1/((pk+pk1)*(xi+ek/pk1))
   df=-fcl*(1+ek*ek)*d
   if (fcl.ne.one) fcl=-1
   if (d.lt.zero) fcl=-fcl
   f=f+df

   !--begin CF1 loop on Pk = K = lambda + 1

   p=1
  20 continue
      pk=pk1
      pk1=pk1+1
      ek=eta/pk
      tk=(pk+pk1)*(xi+ek/pk1)
      d=tk-d*(1+ek*ek)
      if (abs(d).le.acch) then
         write(nsyse,'(/'' CF1 ACCURACY LOSS: D,Df,Acch,K,Eta/K,Eta,X = '',&
           &1p,7e9.2/)') d,df,acch,pk,ek,eta,x
         write(nsyso,'(/'' CF1 ACCURACY LOSS: D,Df,Acch,K,Eta/K,Eta,X = '',&
           &1p,7e9.2/)') d,df,acch,pk,ek,eta,x
         p=p+1
         if (p.gt.two) go to 110
      endif
      d=1/d
      if (d.lt.zero) fcl=-fcl
      df=df*(d*tk-1)
      f=f+df
      if (pk.gt.px) go to 110
      if (abs(df).ge.abs(f)*acc) go to 20

   !--downward recurrence to lambda = Xlm. Array Gc, if present, stores RL

   if (llmax.gt.0) then
      fcl=fcl*tm30
      fpl=fcl*f
      fcp(l1)=fpl
      fc (l1)=fcl
      xl=xll
      rl=1
      el=0
      do lp=1,llmax
         el=eta/xl
         rl=sqrt(1+el*el)
         sl=el+xl*xi
         l=l1-lP
         fcl1=(fcl*sl+fpl)/rl
         fpl=fcl1*sl-fcl*rl
         fcl=fcl1
         fc(l)=fcl
         fcp(l)=fpl
         gc(l+1)=rl
         xl=xl-1
      enddo
      if (fcl.eq.zero) fcl=acc
      f =fpl/fcl
      !--now we have reached lambda = 0
      !--evaluate CF2 = P + I.Q  Again using Steed's algorithm
   endif

   if (xlturn) call jwkb(x,eta,max(xlm,zero),fjwkb,gjwkb,iexp)

   if (iexp.gt.1.or.gjwkb.gt.1/(acch*ten2)) then

      !--arrive here if [G(Xlm).GT.10**6] or [(Iexp.GT.250+Xlturn)=.true.]
      !--arrive here if [G(Xlm).GT.10**6] or [(Iexp.GT.70 &Xlturn)=.true.]
      !--In other words, where values are extreme
      w=fjwkb
      gam=gjwkb*w
      p=f
      q=1

   else

      xlturn=.false.
      ta=2*abort
      pk=0
      wi=eta+eta
      p=0
      q=1-eta*xi
      ar=-e2mm1
      ai=eta
      br=2*(x-eta)
      bi=2
      dr= br/(br*br+bi*bi)
      di=-bi/(br*br+bi*bi)
      dp=-xi*(ar*di+ai*dr)
      dq= xi*(ar*dr-ai*di)
     30 continue
         p=p+dp
         q=q+dq
         pk=pk+2
         ar=ar+pk
         ai=ai+wi
         bi=bi+2
         d =ar*dr-ai*di+br
         di=ai*dr+ar*di+bi
         c=1/(d*d+di*di)
         dr=c*d
         di=-c*di
         a=br*dr-bi*di-1
         b=bi*dr+br*di
         c=dp*a-dq*b
         dq=dp*b+dq*a
         dp=c
         if (pk.gt.ta) go to 120
         if (abs(dp)+abs(dq).ge.(abs(P)+abs(q))*Acc) go to 30
      paccq=half*acc
      if (abs(q).lt.one) paccq=half*acc/abs(q)
      if (abs(p).gt.abs(q)) paccq=paccq*abs(p)
      if (q.le.acc4*abs(p)) go to 130
      gam=(f-p)/q
      !--solve for Fcm=F at lambda=Xlm, then find norm factor W=W/Fcm
      w=1/sqrt((f-p)*gam+q)

   endif

   fcm=sign(w,fcl)
   fc(1)=fcm
   if (.not.xlturn) then
      gcl=fcm*gam
   else
      gcl=gjwkb
   endif
   gc(1)=gcl
   gpl= gcl*(p-q/gam)
   gcp(1)=gpl
   fcp(1)=fcm*f

   !--upward recurrence from Gc(1),Gcp(1) stored value is RL
   !--renormalize Fc,Fcp at each lambda and correct regular derivative
   !--Xl = Xlm here and RL = One, EL = Zero for bessels.
   w=w/abs(fcl)
   do l=1,llmax
      xl=xl+1
      el=eta/xl
      rl=gc(l+1)
      sl=el+xl*xi
      gcl1=(sl*gcl-gpl)/rl
      gpl=rl*gcl-sl*gcl1
      gcl=gcl1
      gc(l+1)=gcl1
      gcp(l+1)=gpl
      fcp(l+1)=w*fcp(l+1)
      fc (l+1)=w*fc (l+1)
   enddo

   !--Generate penetrability, shift factor, and phase shift, & derivatives
   !--Note that Iexp = 1 means "normal" version worked...
   l=lll+1
   if (iexp.gt.1) then
      if (iexp.lt.150) then
         asq=gc(l)**2
         aaa=ten**(-iexp*2)
         pcoul=xx/asq*aaa
         if (jdopha.gt.0) then
            sinphi=fc(L)/gc(L)*aaa
            cosphi=1-sinphi**2
            if (jdoder.gt.0) then
               dphi=gcp(l)/gc(l)*sinphi-fcp(l)/gc(l)*aaa
            endif
         endif
      endif
      if (ishift.gt.0) then
         sss=xx*gcp(l)/gc(l)
         scoul=sss
      endif
      ! Note we've already set defaults to Zero so no need to repeat
   else
      asq=fc(l)**2+gc(l)**2
      pcoul=xx/asq
      sss=xx*(fc(l)*fcp(l)+gc(l)*Gcp(l))/asq
      if (ishift.gt.0) scoul=sss
      if (jdoder.gt.0) dpcoul=(1-2*sss)/asq
      if (jdopha.gt.0) then
         a=sqrt(asq)
         sinphi=fc(l)/a
         cosphi=gc(l)/a
           ! Sinsqr = Sin^2 (phase shift)
           ! Sin2ph = sin ( 2 * phase shift )
         if (jdoder.gt.0) dphi=(gcp(l)*fc(l)-gc(l)*fcp(l))/asq
      endif
   endif

   return

   !--error messages

  100 continue
   !--arrives here if xx.le.sqrt(accur) or if negative
   ifail=-1
   write(nsyse,'('' FOR Xx = '',1p,e12.3,'' TRY SMALL-X SOLUTIONS'',&
     &'' OR X NEGATIVE''/,'' SQUARE ROOT ACCURACY PARAMETER = '',&
     &e12.3/)') xx,acch
   return

  110 continue
   ifail=1
   write(nsyse,'('' CF1 HAS FAILED TO CONVERGE AFTER '',f10.0,&
     &'' ITERATIONS''/'' F,Df,Pk,Px,Accur =  '',1p,5e12.3//)')&
     abort,f,df,pk,px,acc
   return

  120 continue
   ifail=2
   write(nsyse,'('' CF2 HAS FAILED TO CONVERGE AFTER '',f7.0,&
     &''ITERATIONS''/'' P,Q,Dp,Dq,Accur =  '',1p,4e17.7,e12.3//)')&
     abort,p,q,dp,dq,acc
   return

  130 continue
   ifail=3
   dpacc=abs(p)*acc*10**4
   write(nsyse,'('' Final Q.le.abs(P)*Acc*10**4'',/,&
     &''       Q, ... , P, Acc = '',1p,4e12.3/&
     &''       Llmax = '', i5/)') q,dpacc,p,acc,llmax
   return

   end subroutine coulfg

   subroutine jwkb(xx,eta1,xl,fjwkb,gjwkb,iexp)
   !-------------------------------------------------------------------
   ! Computes jwkb approximations to Coulomb functions for xl.ge.0
   !-------------------------------------------------------------------
   ! externals
   real(kr)::xx,eta1,xl,fjwkb,gjwkb
   integer::iexp
   ! internals
   real(kr)::x,eta,gh2,xll1,hll,hl,sl,rl2,gh,phi,phi10
   real(kr),parameter::zero=0
   real(kr),parameter::half=0.5e0_kr
   real(kr),parameter::ten=10
   real(kr),parameter::aloge=0.434294481903251816667932e0_kr
   ! aloge is log-base-10 of e ==> e=exp(1.0),aloge=dlog10(e)
   real(kr),parameter::six35=0.1714285714285714285714285714285714285714285e0_kr
   ! six35 is 6.0/35.0

   x=xx
   eta=eta1
   gh2=x*(eta+eta-x)
   xll1=max(xl*xl+xl,zero)
   if (gh2+xll1.le.zero) return
   hll=xll1+six35
   hl=sqrt(hll)
   sl=eta/hl+hl/x
   rl2=1+eta*eta/hll
   gh=sqrt(gh2+hll)/x
   phi=x*gh-half*(hl*log((gh+sl)**2/rl2)-log(gh))
   if (eta.ne.zero) phi=phi-eta*atan2(x*gh,x-eta)

   phi10=-phi*aloge
   iexp=int(phi10)
   if (iexp.gt.70) then
      gjwkb=ten**(phi10-iexp)
   else
      gjwkb=exp(-phi)
      iexp=0
   endif
   fjwkb=half/(gh*gjwkb)

   return
   end subroutine jwkb

   subroutine coulx(eeta,rrho,lll,llmax,f,fpr,g,gpr,sigma,p,s,&
     dp,sinphi,cosphi,dphi,jdopha,jdoder,ishift)
   !-------------------------------------------------------------------
   ! Coulomb routine
   !-------------------------------------------------------------------
   ! externals
   real(kr)::eeta,rrho,p,s,dp,sinphi,cosphi,dphi
   real(kr)::f(*),fpr(*),g(*),gpr(*),sigma(*)
   integer::lll,llmax,jdopha,jdoder,ishift
   ! internals
   real(kr)::eta,rho,rhoi,u,upr,sigma0,g0,g0pr
   integer::lmax
   real(kr),parameter::zero=0
   real(kr),parameter::five=5

   eta=eeta
   rho=rrho
   lmax=llmax

   if (eta.gt.10*rho.and.eta.gt.five) then
      !--Here for eta>>rho
      call bigeta(eta,rho,lll,llmax,f,fpr,g,gpr,&
        p,s,dp,sinphi,cosphi,dphi,jdopha,jdoder,ishift)

   else
      !--Here for eta and rho not so very different

      if (eta.ge.five) then
         !--Generate U, Upr, first for Rhoi.ne.Rho
         !--i.e., use asymptotic formula for large Eta to give U=G0(2*Eta,Eta)
         call asymp1(eta,rhoi,u,upr)
         !--Now generate Taylor series expansion of G0(Rho) around G0(Rhoi),
         !--with Rhoi redefined if necessary to obtain convergence
         call taylor(eta,rho,u,upr,rhoi)

      else
         !--Generate Sigma to use in asymptotic expansion for this Rhoi
         call xsigll(eta,sigma,sigma0,lmax)
         !--Use asymptotic formula to give G0(Rhoi,Eta)
         call asymp2(eta,rho,u,upr,rhoi,sigma0)
         if (abs(u).le.1.0e25_kr) then
            call taylor(eta,rho,u,upr,rhoi)
         endif
      endif

      !--Now find P, S, etc...
      g0=u
      g0pr=upr
      if (abs(g0).gt.1.e+25_kr) then
         call end1(llmax,f,fpr,g,gpr,g0,g0pr)
      else
         call getfg(eta,rho,llmax,lll,f,fpr,g,gpr,g0,g0pr)
      endif
      call getps(rho,lll,f,fpr,g,gpr,p,s,dp,sinphi,cosphi,&
         dphi,jdopha,jdoder,ishift)

   endif

   return
   end subroutine coulx

   subroutine asymp1(eeta,rhoi,u,upr)
   !-------------------------------------------------------------------
   ! Calculate asymptotic expansion of G0 & G0PR from
   ! Eqs. 14.5.12b and 14.5.13b in Abromowitz & Stegun.
   ! Note that these formulae are valid for Rhoi=2*Eta.
   !-------------------------------------------------------------------
   ! externals
   real(kr)::eeta,rhoi,u,upr
   ! internals
   real(kr)::eta,ceta,seta,temp

   eta=eeta
   ceta=(eta)**0.3333333333333333e0_kr
   seta=sqrt(ceta)
   temp=1/ceta**2
   u=1.223404016e0_kr*seta*&
     (1+temp**2*(.04959570165e0_kr+temp*(-.008888888889e0_kr+&
     temp**2*(.002455199181e0_kr+temp*(-.0009108958061e0_kr+&
     temp**2*.0002534684115e0_kr)))))
   upr=-.7078817734e0_kr*&
     (1+temp*(-.1728260369e0_kr+temp**2*(.0003174603174e0_kr+&
     temp*(-.003581214850e0_kr+temp**2*(.0003117824680e0_kr&
     -temp*.0009073966427e0_kr)))))/seta
   rhoi=2*eta

   return
   end subroutine asymp1

   subroutine xsigll(eeta,sigma,sigma0,lmax)
   !-------------------------------------------------------------------
   ! Generate sigma(LL) for LL=1 thru Lmax+1  (Ie L=0 thru Lmax)
   !-------------------------------------------------------------------
   use physics, only : pi,euler ! get pi,euler
   ! externals
   integer::lmax
   real(kr)::eeta,sigma0
   real(kr)::sigma(*)
   ! internals
   real(kr)::eta,peta,sum,xi,xm,sumas,s,temp1,as,add,xl
   integer::i,m,is,k,j,ll
   integer::mmmxxx=100000
   real(kr),parameter::small=.000001e0_kr
   real(kr),dimension(5),parameter::ber=(/&
     0.1666666666666666666666666666666666667e0_kr,&
     -0.0333333333333333333333333333333333333e0_kr,&
     0.0238095238095238095238095238095238095e0_kr,&
     -0.0333333333333333333333333333333333333e0_kr,&
     0.0757575757575757575757575757575757576e0_kr/)
     ! ber [ 1/6, -1/30, 1/42, -1/30, 5/66 ]
   real(kr),parameter::zero=0

   eta=eeta
   peta=abs(eta)

   if (peta.ge.3.0e0_kr) then

      ! Here (abs(Eta).ge.3) so use table BER to estimate sigma0
      sum=0
      do i=1,5
         xi=i
         m=2*i-1
         xm=m
         sum=sum+ber(i)/(2*xi*xm*(peta**m))
      enddo
      sigma0=Pi/4+peta*(log(peta)-1)-sum
       ! Eq. 14.6.16 gives all but "-Sum" of this equation ???

   else

      ! Here abs(Eta).LT.3 so generate Sigma0 more accurately, from
      ! Eq. 14.5.6 pg 650 Abramowitz & Stegun
      sumas=0
      do is=1,mmmxxx
         s=is
         temp1=peta/s
         if (s.le.2*peta) then
            ! This is exact from Eq. 6.1.27 for arg(Gamma(1+peta/S))
            !  using Eq. 6.3.2 for Psi(1)
            as=temp1-atan(temp1)
         else
           ! This is approximation for atan(Temp1) via Taylor expansion
           ! for Peta/S -> 0
            as=0
            k=0
            do j=1,mmmxxx
               m=j+j+1
               xm=m
               add=(temp1**m)/xm
               if (k.eq.0) then
                  as=as+add
                  k=1
               else
                  as=as-add
                  k=0
               endif
               if (abs(add/as).le.small) exit
            enddo
         endif
         sumas=sumas+as
         if (abs(as/sumas).le.small) exit
      enddo
      sigma0=-euler*peta+sumas

   endif
   if (eta.lt.zero) sigma0=-sigma0
   sigma(1)=sigma0

   !--Now set Sigma(LL) for all L
   if (lmax.gt.0) then
      do ll=1,lmax
         xl=ll
         sigma(ll+1)=sigma(ll)+atan(eta/xl)
      enddo
   endif

   return
   end subroutine xsigll

   subroutine asymp2(eeta,rrho,u,upr,rhoi,sigma0)
   !-------------------------------------------------------------------
   ! Calculate U, Upr for large Rhoi but finite Eta using Eqs. 14.5.1-8
   ! on pg 540 Abromowitz & Stegun.  Note that Rhoi is chosen to be
   ! "large enough" so that formula is valid.
   !-------------------------------------------------------------------
   ! externals
   real(kr)::eeta,rrho,u,upr,rhoi,sigma0
   ! internals
   integer::i,jcheck,n,icheck
   real(kr)::eta,rho,xn,temp,an,bn,temp2,w
   real(kr)::phi,sinphi,cosphi,g0,g0pr
   real(kr)::zold(4),znew(4),z(4),bigz(4)
   real(kr),parameter::del=100
   real(kr),parameter::epslon=0.000001e0_kr
   real(kr),parameter::zero=0
   real(kr),parameter::ten=10

   eta=eeta
   rho=rrho
   rhoi=max(rho*2,ten,ten*eta)

  10 continue

   !--initialize
   xn=0
   zold(1)=1
   zold(2)=0
   zold(3)=0
   zold(4)=1-eta/rhoi
   do i=1,4
      z(i)=zold(i)
      bigz(i)=abs(z(i))
   enddo

   jcheck=0
   do n=1,100
      temp=2*(xn+1)*rhoi
      an=(2*xn+1)*eta/temp
      bn=(eta*eta-xn*(xn+1))/temp
      xn=xn+1
      znew(1)=an*zold(1)-bn*zold(2)
      znew(2)=an*zold(2)+bn*zold(1)
      znew(3)=an*zold(3)-bn*zold(4)-znew(1)/rhoi
      znew(4)=an*zold(4)+bn*zold(3)-znew(2)/rhoi
      icheck=0
      do i=1,4
         z(i)=z(i)+znew(i)
         zold(i)=znew(i)
         temp2=abs(z(i))
         bigz(i)=max(bigz(i),abs(z(i)))
         if (bigz(i)/temp2.gt.del) go to 50
         if (abs(znew(i)/z(i)).le.epslon) icheck=icheck+1
      enddo
      w=z(1)*z(4)-z(2)*z(3)
      if (abs(w).gt.ten) go to 50
      if (icheck.eq.4) then
         jcheck=jcheck+1
         if (jcheck.ge.4) go to 60
      else
         jcheck=0
      endif
   enddo

  50 continue

   !--oops.  rhoi isn't big enough for the asymptotic formula to converge.
   !--Double rhoi & try again.
   rhoi=rhoi*2
   go to 10

  60 continue
   !--here the formula for z's has converged.  Ergo calculate u, upr, rhoi:
   phi=rhoi-eta*log(2*rhoi)+sigma0
   cosphi=cos(phi)
   sinphi=sin(phi)
   g0=z(1)*cosphi-z(2)*sinphi
   g0pr=z(3)*cosphi-z(4)*sinphi
   u=g0
   upr=g0pr

   return
   end subroutine asymp2

   subroutine taylor(eeta,rrho,u,upr,rhoi)
   !-------------------------------------------------------------------
   !
   ! Do Taylor series integration of G0, G0pr for arg=Rho,
   ! starting from the (now known) values at arg=Rhoi
   !
   !
   ! *** Find solution of differential equation
   !      u" + (1-2*Eta/Rho) u = 0    (Coulomb eqn for L=0, 14.1.1 A&S)
   !  via Taylor expansion
   !      u(Rho) = u + u' d + u" d**2/2 + u"' d**3/6 + u"" d**4/4! + ...
   !  where right-hand-side is evaluated at Rhoi
   !  and where d = Delta = Rho-Rhoi
   !
   !  Rewriting u(Rho) = Sum  a(n)  where n=1 to infinity, with
   !
   !       a(1) = u
   !       a(2) = u'  d
   !       a(3) = u"  d**2 / 2!
   !       a(4) = u"' d**3 / 3!
   !       a(5) = u"" d**4 / 4!
   !       ...
   !       a(n+1) = u[n] d**n / n!
   !
   !
   !  Substituting Equation (1) into a(3) gives
   !       a(3) = {-(1-2 Eta/Rhoi) u } d**2/2
   !            = - (1-2 Eta/Rhoi) d**2 / 2 * a(1)
   !
   !       a(4) = - { 2 Eta/Rhoi**2 u + (1-2 Eta/Rhoi) u'} d**3/ 3!
   !            = - { 2 Eta/Rhoi**2 a(1) + (1-2 Eta/Rhoi) a(2) /d} d**3/3!
   !            = - {   2 Eta/Rhoi**2  d**3   /3! a(1)
   !                  + (1-2 Eta/Rhoi) d**2 2!/3! a(2)  }
   !
   !            = - {     1/Rhoi       d    /3     a(3)
   !                  + (1-2 Eta/Rhoi) d**2 /(3*2) a(2)
   !                  +                d**3 /(3*2) a(1) }  using (3) above
   !
   !       a(5) = {4Eta/Rhoi**3 u - 4Eta/Rhoi**2 u'-(1-2 Eta/Rhoi) u"}d**4/4!
   !            = - {   2/Rhoi        d**3   /4! a(2)
   !                + (1-2 Eta/Rhoi) d**2 2!/4! a(3)
   !                +  2/Rhoi        d    3!/4! a(4) }  use a(4) to replace a(1)
   !
   !       a(n) = - {    1/Rhoi       d    (n-3) /(n-1)        a(n-1)
   !                 + (1-2 Eta/Rhoi) d**2       /((n-1)(n-2)  a(n-2)
   !                 +                d**3       /((n-1)(n-2)  a(n-3) }
   !
   !-------------------------------------------------------------------
   ! externals
   real(kr)::eeta,rrho,u,upr,rhoi
   ! internals
   integer::nstart,jcheck,n,m
   real(kr)::eta,rho,delta,sum,sumpr,big,bigpr,xn,temp
   real(kr)::a(100)
   real(kr),parameter::epslon=1.0e-6_kr
   real(kr),parameter::bigger=1.0e10_kr
   real(kr),parameter::biggst=1.0e30_kr
   real(kr),parameter::del=100
   real(kr),parameter::zero=0

   eta=eeta
   rho=rrho
   delta=rho-rhoi
   if (delta.eq.zero) return

  10 continue

   a(1)=u
   a(2)=delta*upr
   a(3)=-delta*delta/2*(1-2*eta/rhoi)*a(1)
   nstart=4

  20 continue
   jcheck=0
   sum=0
   sumpr=0
   big=0
   bigpr=0
   do n=1,100
      xn=n-1
      if (n.ge.nstart) then
         a(n)=-(delta*(xn-1)*(xn-2)*a(n-1)+&
           (rhoi-2*eta)*(delta**2)*a(n-2)+(delta**3)*a(n-3))/&
           (rhoi*(xn-1)*xn)
         if (a(n).gt.bigger) go to 40
      endif
      sum=sum+a(n)
      sumpr=sumpr+xn*a(n)
      if (sum.ge.biggst) go to 40
      if (sumpr.ge.biggst) go to 40
      big=max(big,abs(sum))
      bigpr=max(bigpr,abs(sumpr))
      if (sum.eq.zero.or.sumpr.eq.zero) then
         jcheck=0
      else
         if (abs(big/sum).ge.del) go to 40
         if (abs(bigpr/sumpr).ge.del) go to 40
         if (abs(a(n)/sum).ge.Epslon.or.abs(xn*a(n)/sumpr).ge.epslon) then
            jcheck=0
         else
            jcheck=jcheck+1
            if (jcheck.ge.4) go to 60
         endif
      endif
   enddo
   n=100

  40 continue
   !--The series did not converge for u(rho) using Taylor expansion
   !--around arg=rhoi.  Ergo, find u(rhoi+x), where x=delta/2, using
   !--Taylor expansion around arg=rhoi.  Assuming this converges,
   !--redefine rhoi to be rhoi+x and try again to get u(rho).
   nstart=max0(nstart,n+1)
   m=nstart-1
   delta=delta/2
   temp=2
   do n=1,m
      temp=temp/2
      a(n)=a(n)*temp
   enddo
   go to 20

  60 continue
   !--Here we know sum = u(rhoi+delta).  Redefine rhoi & delta and try
   !--again to find u(rho) as Taylor expansion around arg=rhoi.
   u=sum
   upr=sumpr/delta
   rhoi=rhoi+delta
   delta=rho-rhoi
   if (abs(delta).ge.epslon) go to 10

   return
   end subroutine taylor

   subroutine end1(llmax,f,fpr,g,gpr,g0,g0pr)
   !-------------------------------------------------------------------
   ! externals
   integer::llmax
   real(kr)::f(*),fpr(*),g(*),gpr(*)
   real(kr)::g0,g0pr

   llmax=0
   g(1)=g0
   gpr(1)=g0pr
   f(1)=0
   fpr(1)=0

   return
   end subroutine end1

   subroutine getfg(eta,rho,llmax,lll,f,fpr,g,gpr,g0,g0pr)
   !-------------------------------------------------------------------
   ! Calculate G & F for all L, given G(L=0) and deriv G(L=0)
   !-------------------------------------------------------------------
   ! externals
   real(kr)::eta,rho,g0,g0pr
   integer::llmax,lll
   real(kr)::f(*),fpr(*),g(*),gpr(*)
   ! internals
   integer::lmax,l,limit,il,j,ll
   real(kr)::xl,temp1,gm2,gm1,gm,fp1,fp2,temp2,fp
   real(kr),parameter::big=1.0e12_kr

   lmax=llmax

   !--set G(L)
   limit=max0(3,lmax+1)
   g(1)=g0
   gpr(1)=g0pr
   g(2)=((eta+1/rho)*g(1)-gpr(1))/sqrt(eta**2+1)
    ! from Abramowitz & Stegun Eq. 14.2.2
   do l=3,limit
      xl=l-1
      temp1=sqrt(xl*xl+eta*eta)
      g(l)=(2*xl-1)/temp1*(eta/(xl-1)+xl/rho)*g(l-1)-&
        xl/temp1*sqrt(1+(eta/(xl-1))**2)*g(l-2)
      ! from Abramowitz & Stegun Eq. 14.2.3 rewritten
      if (abs(g(l)).gt.big.and.l.gt.lll) then
         limit=l
         exit
      endif
   enddo

   !--Find maximum L value to use in figuring F(Limit)
   !--I.e. find J such that G(J) < 1.e-4 * G(J-1) three times in a row
   gm2=g(limit-1)
   gm1=g(limit  )
   il=- 1
   do j=limit,10000
      xl=j
      temp1=sqrt(xl*xl+eta*eta)
      gm=(2*xl-1)/temp1*(eta/(xl-1)+xl/rho)*gm1-&
         xl/temp1*sqrt(1+(eta/(xl-1))**2)*gm2
       !   again Eq. 14.2.3
       ! Gmh    IF ( (G(Limit)/Gm)**2 .gt. 1.0d-8 ) IL=-2
      if (abs(g(limit)/gm).gt.1.0e-4_kr) il=-2
      if (il.gt.0) exit
      il=il+1
      gm2=gm1
      gm1=gm
   enddo

   !--Figure approximate F(Limit+3), F(Limit+4), ... F(J-1) in reverse
   !--order.  Do not store these numbers anywhere permanent
   xl=j
   fp1=xl/gm/sqrt(xl**2+eta**2)
     ! from Abramowitz & Stegun Eq. 14.2.5 with
     ! F(L)G(L-1) assumed to be negligible
   fp2=0
   l=j-1
   do ll=1,j-3-limit
      l=l-1
      xl=l
      temp2=sqrt((xl+1)**2+eta**2 )
      fp=((2*xl+3)*(eta/(xl+2)+(xl+1)/rho)*fp1-&
        (xl+1)*sqrt(1+(eta/(xl+2))**2)*fp2)/temp2
        !  from Abramowitz & Stegun Eq. 14.2.3 again, in reverse
      fp2=fp1
      fp1=fp
   enddo

   !--Now figure F(1), F(2), ..., F(Limit), F(Limit+1),
   !--again in reverse order.  Store these in array F(L)
   do ll=1,limit+2
      l=l-1
      xl=l
      temp2=sqrt((xl+1)**2+eta**2)
      fp=((2*xl+3)*(eta/(xl+2)+(xl+1)/rho)*fp1-&
        (xl+1)*sqrt(1+(eta/(xl+2))**2)*fp2)/temp2
      f(l+1)=fp
      fp2=fp1
      fp1=fp
   enddo

   !--Generate derivative functions Fpr & Gpr for L=1,Limit
   !--(Remember, we already have Gpr(1))
   fpr(1)=(1/rho+eta)*f(1)-sqrt(1+eta**2)*f(2)
   do l=2,limit
      xl=l
      fpr(l)=(xl/rho+eta/xl)*f(l)-sqrt(1+(eta/xl)**2)*f(l+1)
      temp1=eta/(xl-1)
      gpr(l)=sqrt(1+temp1**2)*g(l-1)-((xl-1)/rho+temp1)*g(l)
   enddo

   if (limit.le.lmax) llmax=limit-1

   return
   end subroutine getfg

   subroutine bigeta(eeta,rrho,lll,llmax,f,fpr,g,gpr,&
     p,s,dp,sinphi,cosphi,dphi,jdopha,jdoder,ishift)
   !-------------------------------------------------------------------
   ! Formulae 14.6.7-8 page 542 in Abramowitz & Stegun, for Eta >> Rho
   !-------------------------------------------------------------------
   use util ! provides error
   use physics, only : pi,euler ! get pi,euler ! get pi,euler
   ! externals
   real(kr)::f(*),fpr(*),g(*),gpr(*)
   real(kr)::eeta,rrho,p,s,dp,sinphi,cosphi,dphi
   integer::lll,llmax,jdopha,jdoder,ishift
   ! internals
   real(kr)::eta,rho,q,zhalf,z,sum,a,ai0,b,ak0,ai1,c,ak1,d
   real(kr)::g0,g0pr
   integer::k,n
   real(kr),parameter::zero=0
   real(kr),parameter::half=0.5e0_kr

   eta=eeta
   rho=rrho
   q=2*rho*eta
   zhalf=sqrt(q)
   z=zhalf*2

   !--Evaluate I_0(z) from 9.6.12 A&S
   sum=1
   a=q
   do k=1,100
      if (sum+a.eq.sum) go to 10
      sum=sum+a
      a=a*q/(k+1)**2
   enddo
   call error('bigeta','I0 sum failed',' ')
  10 continue
   ai0=sum

   !--Evaluate K_0(Z) from 9.6.13 A&S
   sum=-(log(zhalf)+euler)*ai0
   a=q
   b=1
   do k=1,100
      if (sum+a*b.eq.sum) go to 20
      sum=sum+a*b
      b=b+1/(k+1)
      a=a*q/(k+1)**2
   enddo
   call error('bigeta','K0 sum failed',' ')
  20 continue
   ak0=sum

   !--Evaluate I_1(Z) from 9.6.10 A&S
   sum=1
   a=q/2
   do k=1,100
      if (sum+a.eq.sum) go to 30
      sum=sum+a
      a=a*q/((k+1)*(k+2))
   enddo
   call error('bigeta','L1 sum failed',' ')
  30 continue
   ai1=sum*zhalf

   !--Evaluate K_1(Z) from 9.6.11 A&S
   sum=1/z+(log(zhalf)+euler)*ai1-zhalf*half
   a=zhalf*q/2
   b=1
   c=half**2
   do k=1,100
      if (sum-a*(b+c).eq.sum) go to 40
      sum=sum-a*(b+c)
      c=half/(k+2)
      b=b+1/(k+1)
      a=a*q/((k+1)*(k+2))
   enddo
   call error('bigeta','K1 sum failed',' ')
  40 continue
   ak1=sum

   !--Now have I_0(Z), I_1(Z), K_0(Z), K_1(Z)
   !--Ergo, can get F,G,Fp,Gp for L = 0, from Eq. 14.6.8 A&S
   c=sqrt(pi*rho)
   d=sqrt(2*pi*eta)
   f(1)=c*ai1
   fpr(1)=d*ai0
   c=2*c/pi
   d=2*d/pi
   g(1)=c*ak1
   gpr(1)=-d*ak0

   !--Store results because they'll be changed in Getfg
   g0=g(1)
   g0pr=gpr(1)

   if (lll.gt.0) then
      !--Obtain values for all L's
      call getfg(eta,rho,llmax,lll,f,fpr,g,gpr,g0,g0pr)
   endif
   a=pi*eta
   b=exp(-a)
   n=lll+1
   if (f(n)*b*b+g(n).eq.g(n)) then
      p=(rho*b/g(n))/g(n)*b
      if (ishift.gt.0) s=rho*gpr(n)/g(n)
      if (jdopha.ne.0) sinphi=f(n)/g(n)*b*b
      ! No need to zero the others because that's already been done
   else
      f(n)=f(n)*b
      g(n)=g(n)/b
      fpr(n)=fpr(n)*b
      gpr(n)=gpr(n)/b
      call getps(rho,lll,f,fpr,g,gpr,p,s,dp,sinphi,&
        cosphi,dphi,jdopha,jdoder,ishift)
   endif

   return
   end subroutine bigeta

   subroutine getps (rho,lll,f,fpr,g,gpr,p,s,dp,&
     sinphi,cosphi,dphi,jdopha,jdoder,ishift)
   !-------------------------------------------------------------------
   ! externals
   real(kr)::rho,p,s,dp,sinphi,cosphi,dphi
   integer::lll,jdopha,jdoder,ishift
   real(kr)::f(*),fpr(*),g(*),gpr(*)
   ! internals
   integer::n
   real(kr)::asq,a,ss

   n=lll+1
   asq=f(n)**2+g(n)**2
   a=sqrt(asq)
   p=rho/asq
   if (jdoder.ne.0.or.ishift.gt.0) then
      ss=rho*(f(n)*fpr(n)+g(n)*gpr(n))/asq
   endif
   if (jdoder.gt.0) dp=(1-2*ss)/asq
   if (ishift.gt.0) s=ss
   if (jdopha.gt.0) then
      sinphi=f(n)/a
      cosphi=g(n)/a
      if (jdoder.gt.0) dphi=(gpr(n)*f(n)-g(n)*fpr(n))/asq
   endif

   return
   end subroutine getps

   subroutine zeror(nchan,ier)
   !-------------------------------------------------------------------
   ! Calculate trivial xxxx and xq matrices.
   !-------------------------------------------------------------------
   ! externals
   integer::nchan,ier
   ! internals
   integer::n,jk,j,k

   n=nchan
   jk=0
   do j=1,n
      do k=1,j
         jk=jk+1
         xxxxr(jk,ier)=0
         xxxxi(jk,ier)=0
      enddo
   enddo
   do j=1,n
      do k=1,n
         xqr(k,j,ier)=0
         xqi(k,j,ier)=0
      enddo
   enddo

   return
   end subroutine zeror

   subroutine yinvrs(nchan,ier)
   !-------------------------------------------------------------------
   ! Invert ymat to give yinv.
   !-------------------------------------------------------------------
   ! externals
   integer::nchan,ier

   if (nchan.eq.1) then
      call onech(ymat,yinv,ier)
   else if (nchan.eq.2) then
      call twoch(ymat,yinv,ier)
   else if (nchan.eq.3) then
      call threech(ymat,yinv,ier)
   else
      call yfour(ymat,yinv,nchan,ier)
   endif

   return
   end subroutine yinvrs

   subroutine onech(rmat,rinv,ier)
   !-------------------------------------------------------------------
   ! Invert matrix for one-channel case.
   !-------------------------------------------------------------------
   ! externals
   real(kr),dimension(2,ntriag,*)::rmat,rinv
   integer::ier
   ! internals
   real(kr)::aa
   integer::kkkk11=1
   real(kr),parameter::zero=0

   if (rmat(1,1,ier)+rmat(2,1,ier).eq.rmat(1,1,ier)) then
      rinv(1,1,ier)=1/rmat(1,1,ier)
      rinv(2,1,ier)=-(rmat(2,1,ier)/rmat(1,1,ier))/rmat(1,1,ier)
   else if (rmat(1,1,ier)+rmat(2,1,ier).eq.rmat(2,1,ier)) then
      rinv(1,1,ier)=(rmat(1,1,ier)/rmat(2,1,ier))/rmat(2,1,ier)
      rinv(2,1,ier)=-1/rmat(2,1,ier)
   else if (rmat(1,1,ier).eq.zero) then
      rinv(1,1,ier)=0
      rinv(2,1,ier)=-1/rmat(2,1,ier)
   else
      aa=rmat(1,kkkk11,ier)**2+rmat(2,kkkk11,ier)**2
      rinv(1,kkkk11,ier)=rmat(1,kkkk11,ier)/aa
      rinv(2,kkkk11,ier)=-rmat(2,kkkk11,ier)/aa
   endif

   return
   end subroutine onech

   subroutine twoch(rmat,rinv,ier)
   !-------------------------------------------------------------------
   ! Invert matrix for two-channel case.
   !-------------------------------------------------------------------
   ! externals
   real(kr),dimension(2,ntriag,*)::rmat,rinv
   integer::ier
   ! internals
   integer::k,i
   real(kr)::a,bbr,bbi,aa,aar,aai
   integer::kkkk11=1
   integer::kkkk21=2
   integer::kkkk22=3
   real(kr),parameter::zero=0

   if (rmat(1,kkkk11,ier).ne.zero.or.rmat(1,kkkk21,ier).ne.zero.or.&
     rmat(1,kkkk22,ier).ne.zero) then

      a=0
      k=0
      do i=1,3
         if (abs(rmat(1,i,ier)).gt.a) k=i
         if (abs(rmat(1,i,ier)).gt.a) a=abs(rmat(1,i,ier))
         if (abs(rmat(2,i,ier)).gt.a) k=i
         if (abs(rmat(2,i,ier)).gt.a) a=abs(rmat(2,i,ier))
      enddo
      if (a.gt.zero) then
         if (k.eq.2) then
            do i=1,3
               rmat(1,i,ier)=rmat(1,i,ier)/a
               rmat(2,i,ier)=rmat(2,i,ier)/a
            enddo
         else
            rmat(1,k,ier)=rmat(1,k,ier)/a
            rmat(2,k,ier)=rmat(2,k,ier)/a
            rmat(1,2,ier)=rmat(1,2,ier)/sqrt(a)
            rmat(2,2,ier)=rmat(2,2,ier)/sqrt(a)
         endif
      endif

      bbr=rmat(1,kkkk11,ier)*rmat(1,kkkk22,ier)-&
         rmat(2,kkkk11,ier)*rmat(2,kkkk22,ier)-rmat(1,kkkk21,ier)**2+&
         rmat(2,kkkk21,ier)**2
      bbi=rmat(1,kkkk11,ier)*rmat(2,kkkk22,ier)+&
         rmat(2,kkkk11,ier)*rmat(1,kkkk22,ier)-2*rmat(1,kkkk21,ier)*&
         rmat(2,kkkk21,ier)

      if (bbr+bbi.ne.bbi) then
         if (bbr+bbi.ne.bbr) then
            aa=1/(bbr*bbr+bbi*bbi)
            aar=bbr*aa
            aai=-bbi*aa
         else
            aar=1/bbr
            aai=-(bbi/bbr)/bbr
         endif
      else
         aar=(bbr/bbi)/bbi
         aai=-1/bbi
      endif

      rinv(1,kkkk11,ier)= aar*rmat(1,kkkk22,ier)-aai*rmat(2,kkkk22,ier)
      rinv(2,kkkk11,ier)= aar*rmat(2,kkkk22,ier)+aai*rmat(1,kkkk22,ier)
      rinv(1,kkkk21,ier)=-aar*rmat(1,kkkk21,ier)+aai*rmat(2,kkkk21,ier)
      rinv(2,kkkk21,ier)=-aar*rmat(2,kkkk21,ier)-aai*rmat(1,kkkk21,ier)
      rinv(1,kkkk22,ier)= aar*rmat(1,kkkk11,ier)-aai*rmat(2,kkkk11,ier)
      rinv(2,kkkk22,ier)= aar*rmat(2,kkkk11,ier)+aai*rmat(1,kkkk11,ier)

      if (a.ne.zero) then
         if (k.eq.2) then
            do i=1,3
               rinv(1,i,ier)=rinv(1,i,ier)/a
               rinv(2,i,ier)=rinv(2,i,ier)/a
            enddo
         else
            rinv(1,k,ier)=rinv(1,k,ier)/a
            rinv(2,k,ier)=rinv(2,k,ier)/a
            rinv(1,2,ier)=rinv(1,2,ier)/sqrt(a)
            rinv(2,2,ier)=rinv(2,2,ier)/sqrt(a)
         endif
      endif

   else if (rmat(2,kkkk21,ier).ne.zero) then
      !--Here when real part of Rmat is zero everywhere,
      !--imaginary part is dense
      rinv(1,kkkk11,ier)=0
      rinv(1,kkkk21,ier)=0
      rinv(1,kkkk22,ier)=0
      a=rmat(2,kkkk11,ier)*rmat(2,kkkk22,ier)-rmat(2,kkkk21,ier)**2
      rinv(2,kkkk11,ier)=-rmat(2,kkkk22,ier)/a
      rinv(2,kkkk21,ier)= rmat(2,kkkk21,ier)/a
      rinv(2,kkkk22,ier)=-rmat(2,kkkk11,ier)/a

   else
      !--Here when only real part of Rmat is zero everywhere; imaginary
      !--part is zero off-diagonal and non-zero on diagonal
      rinv(1,kkkk11,ier)= 0
      rinv(1,kkkk21,ier)= 0
      rinv(1,kkkk22,ier)= 0
      rinv(2,kkkk11,ier)=-1 /(rmat(2,kkkk11,ier))
      rinv(2,kkkk21,ier)= 0
      rinv(2,kkkk22,ier)=-1/(rmat(2,kkkk22,ier))
   endif

   return
   end subroutine twoch

   subroutine threech(rmat,rinv,ier)
   !-------------------------------------------------------------------
   ! Invert matrix for three-channel case.
   !-------------------------------------------------------------------
   ! externals
   real(kr),dimension(2,ntriag,*)::rmat,rinv
   integer::ier
   ! internals
   integer::izz,i
   real(kr)::fc1r,fc1i,fc2r,fc2i,fc3r,fc3i,aar,aai,dcr,dci
   real(kr)::aa,a1,a2,a3
   integer::kkkk11=1
   integer::kkkk21=2
   integer::kkkk31=4
   integer::kkkk22=3
   integer::kkkk32=5
   integer::kkkk33=6
   real(kr),parameter::zero=0

   call scale3(a1,a2,a3,rmat,ier)

   izz=0
   do i=1,6
      if (rmat(1,i,ier).ne.zero) izz=1
   enddo

   if (rmat(2,2,ier).ne.zero) izz=1
   if (rmat(2,4,ier).ne.zero) izz=1
   if (rmat(2,5,ier).ne.zero) izz=1

   if (izz.eq.0) then

     !--Here if only imaginary parts of diagonal terms are non zero
      do i=1,6
         rinv(1,i,ier)=0
         rinv(2,i,ier)=0
      enddo
      rinv(2,1,ier)=-1/rmat(2,1,ier)
      rinv(2,3,ier)=-1/rmat(2,3,ier)
      rinv(2,6,ier)=-1/rmat(2,6,ier)

   else

      !--Here real parts of diagonal terms are non zero
      fc1r=rmat(1,kkkk22,ier)*rmat(1,kkkk33,ier)-&
        rmat(2,kkkk22,ier)*rmat(2,kkkk33,ier)&
        -rmat(1,kkkk32,ier)**2+rmat(2,kkkk32,ier)**2
      fc1i=rmat(2,kkkk22,ier)*rmat(1,kkkk33,ier)+&
        rmat(1,kkkk22,ier)*rmat(2,kkkk33,ier)&
        -2*rmat(1,kkkk32,ier)*rmat(2,kkkk32,ier)
      fc2r=rmat(1,kkkk21,ier)*rmat(1,kkkk33,ier)-&
        rmat(2,kkkk21,ier)*rmat(2,kkkk33,ier)&
        -rmat(1,kkkk31,ier)*rmat(1,kkkk32,ier)+&
        rmat(2,kkkk31,ier)*rmat(2,kkkk32,ier)
      fc2i=rmat(2,kkkk21,ier)*rmat(1,kkkk33,ier)+&
        rmat(1,kkkk21,ier)*rmat(2,kkkk33,ier)&
        -rmat(1,kkkk31,ier)*rmat(2,kkkk32,ier)-&
        rmat(2,kkkk31,ier)*rmat(1,kkkk32,ier)
      fc3r=rmat(1,kkkk21,ier)*rmat(1,kkkk32,ier)-&
        rmat(2,kkkk21,ier)*rmat(2,kkkk32,ier)&
        -rmat(1,kkkk31,ier)*rmat(1,kkkk22,ier)+&
        rmat(2,kkkk31,ier)*rmat(2,kkkk22,ier)
      fc3i=rmat(2,kkkk21,ier)*rmat(1,kkkk32,ier)+&
        rmat(1,kkkk21,ier)*rmat(2,kkkk32,ier)&
        -rmat(2,kkkk31,ier)*rmat(1,kkkk22,ier)-&
        rmat(1,kkkk31,ier)*rmat(2,kkkk22,ier)

      aar=rmat(1,kkkk11,ier)*fc1r-rmat(1,kkkk21,ier)*fc2r&
        +rmat(1,kkkk31,ier)*fc3r&
        -(rmat(2,kkkk11,ier)*fc1i-rmat(2,kkkk21,ier)*fc2i&
        +rmat(2,kkkk31,ier)*fc3i)
      aai=rmat(1,kkkk11,ier)*fc1i-rmat(1,kkkk21,ier)*fc2i&
        +rmat(1,kkkk31,ier)*fc3i&
        +(rmat(2,kkkk11,ier)*fc1r-rmat(2,kkkk21,ier)*fc2r&
        +rmat(2,kkkk31,ier)*fc3r)

      if (aar+aai.ne.aai) then
         if (aar+aai.ne.aar) then
            aa=aar**2+aai**2
            dcr=aar/aa
            dci=-aai/aa
         else
            dcr=1/aar
            dci=-(aai/aar)/aar
         endif
      else
         dcr=(aar/aai)/aai
         dci=-1/aai
      endif

      rinv(1,kkkk11,ier)=fc1r*dcr-fc1i*dci
      rinv(2,kkkk11,ier)=fc1r*dci+fc1i*dcr
      rinv(1,kkkk21,ier)=-fc2r*dcr+fc2i*dci
      rinv(2,kkkk21,ier)=-fc2r*dci-fc2i*dcr
      rinv(1,kkkk22,ier)=(rmat(1,kkkk11,ier)*rmat(1,kkkk33,ier)&
        -rmat(1,kkkk31,ier)**2&
        -rmat(2,kkkk11,ier)*rmat(2,kkkk33,ier)+rmat(2,kkkk31,ier)**2)*dcr&
        -(rmat(2,kkkk11,ier)*rmat(1,kkkk33,ier)&
        -2*rmat(1,kkkk31,ier)*rmat(2,kkkk31,ier)&
        +rmat(1,kkkk11,ier)*rmat(2,kkkk33,ier))*dci
      rinv(2,kkkk22,ier)=(rmat(1,kkkk11,ier)*rmat(1,kkkk33,ier)&
        -rmat(1,kkkk31,ier)**2&
        -rmat(2,kkkk11,ier)*rmat(2,kkkk33,ier)+rmat(2,kkkk31,ier)**2)*dci&
        +(rmat(2,kkkk11,ier)*rmat(1,kkkk33,ier)&
        -2*rmat(1,kkkk31,ier)*rmat(2,kkkk31,ier)&
        +rmat(1,kkkk11,ier)*rmat(2,kkkk33,ier))*dcr
      rinv(1,kkkk31,ier)=fc3r*dcr-fc3i*dci
      rinv(2,kkkk31,ier)=fc3r*dci+fc3i*dcr
      rinv(1,kkkk32,ier)=-(rmat(1,kkkk11,ier)*rmat(1,kkkk32,ier)&
        -rmat(1,kkkk21,ier)*rmat(1,kkkk31,ier)&
        -rmat(2,kkkk11,ier)*rmat(2,kkkk32,ier)&
        +rmat(2,kkkk21,ier)*rmat(2,kkkk31,ier))*dcr&
        +(rmat(1,kkkk11,ier)*rmat(2,kkkk32,ier)&
        -rmat(1,kkkk21,ier)*rmat(2,kkkk31,ier)&
        +rmat(2,kkkk11,ier)*rmat(1,kkkk32,ier)&
        -rmat(2,kkkk21,ier)*rmat(1,kkkk31,ier))*dci
      rinv(2,kkkk32,ier)=-(rmat(1,kkkk11,ier)*rmat(1,kkkk32,ier)&
        -rmat(1,kkkk21,ier)*rmat(1,kkkk31,ier)&
        -rmat(2,kkkk11,ier)*rmat(2,kkkk32,ier)&
        +rmat(2,kkkk21,ier)*rmat(2,kkkk31,ier))*dci&
        -(rmat(1,kkkk11,ier)*rmat(2,kkkk32,ier)&
        -rmat(1,kkkk21,ier)*rmat(2,kkkk31,ier)&
        +rmat(2,kkkk11,ier)*rmat(1,kkkk32,ier)&
        -rmat(2,kkkk21,ier)*rmat(1,kkkk31,ier))*dcr
      rinv(1,kkkk33,ier)=(rmat(1,kkkk11,ier)*rmat(1,kkkk22,ier)&
        -rmat(1,kkkk21,ier)**2&
        -rmat(2,kkkk11,ier)*rmat(2,kkkk22,ier)+rmat(2,kkkk21,ier)**2)*dcr&
        -(rmat(2,kkkk11,ier)*rmat(1,kkkk22,ier)&
        -2*rmat(1,kkkk21,ier)*rmat(2,kkkk21,ier)&
        +rmat(1,kkkk11,ier)*rmat(2,kkkk22,ier))*dci
      rinv(2,kkkk33,ier)=(rmat(1,kkkk11,ier)*rmat(1,kkkk22,ier)&
        -rmat(1,kkkk21,ier)**2&
        -rmat(2,kkkk11,ier)*rmat(2,kkkk22,ier)+rmat(2,kkkk21,ier)**2)*dci&
        +(rmat(2,kkkk11,ier)*rmat(1,kkkk22,ier)&
        -2*rmat(1,kkkk21,ier)*rmat(2,kkkk21,ier)&
        +rmat(1,kkkk11,ier)*rmat(2,kkkk22,ier))*dcr

      call unscale3(a1,a2,a3,rmat,rinv,ier)

   endif

   return
   end subroutine threech

   subroutine scale3(a1,a2,a3,rmat,ier)
   !-------------------------------------------------------------------
   ! Scale matrix for inversion.
   !-------------------------------------------------------------------
   ! externals
   real(kr)::a1,a2,a3
   real(kr)::rmat(2,*)
   integer::ier
   ! internals
   real(kr)::bb,cc
   real(kr),parameter::zero=0
   real(kr),parameter::aa=1.e10_kr

   a1=0
   a2=0
   a3=0

   if (abs(rmat(1,1)).ge.aa.or.abs(rmat(2,1)).ge.aa) then
      bb=abs(rmat(1,1))
      cc=abs(rmat(2,1))
      if (cc.gt.bb) bb=cc
      a1=sqrt(bb)
      rmat(1,1)=rmat(1,1)/bb
      rmat(2,1)=rmat(2,1)/bb
      rmat(1,2)=rmat(1,2)/a1
      rmat(2,2)=rmat(2,2)/a1
      rmat(1,4)=rmat(1,4)/a1
      rmat(2,4)=rmat(2,4)/a1
   endif

   if (abs(rmat(1,3)).ge.aa.or.abs(rmat(2,3)).ge.aa) then
      bb=abs(rmat(1,3))
      cc=abs(rmat(2,3))
      if (cc.gt.bb) bb=cc
      a2=sqrt(bb)
      rmat(1,2)=rmat(1,2)/a2
      rmat(2,2)=rmat(2,2)/a2
      rmat(1,3)=rmat(1,3)/bb
      rmat(2,3)=rmat(2,3)/bb
      rmat(1,5)=rmat(1,5)/a2
      rmat(2,5)=rmat(2,5)/a2
   endif

   if (abs(rmat(1,6)).ge.aa.or.abs(rmat(2,6)).ge.aa) then
      bb=abs(rmat(1,6))
      cc=abs(rmat(2,6))
      if (cc.gt.bb) bb=cc
      a3=sqrt(bb)
      rmat(1,4)=rmat(1,4)/a3
      rmat(2,4)=rmat(2,4)/a3
      rmat(1,5)=rmat(1,5)/a3
      rmat(2,5)=rmat(2,5)/a3
      rmat(1,6)=rmat(1,6)/bb
      rmat(2,6)=rmat(2,6)/bb
   endif

   return
   end subroutine scale3

   subroutine unscale3(a1,a2,a3,rmat,rinv,ier)
   !-------------------------------------------------------------------
   ! Unscale matrix after inversion.
   !-------------------------------------------------------------------
   ! externals
   real(kr)::a1,a2,a3
   real(kr),dimension(2,*)::rmat,rinv
   integer::ier
   ! internals
   real(kr)::bb
   real(kr),parameter::zero=0

   if (a1.gt.zero) then
      bb=a1**2
      rmat(1,1)=rmat(1,1)*bb
      rmat(2,1)=rmat(2,1)*bb
      rmat(1,2)=rmat(1,2)*a1
      rmat(2,2)=rmat(2,2)*a1
      rmat(1,4)=rmat(1,4)*a1
      rmat(2,4)=rmat(2,4)*a1
      rinv(1,1)=rinv(1,1)/bb
      rinv(2,1)=rinv(2,1)/bb
      rinv(1,2)=rinv(1,2)/a1
      rinv(2,2)=rinv(2,2)/a1
      rinv(1,4)=rinv(1,4)/a1
      rinv(2,4)=rinv(2,4)/a1
   endif

   if (a2.gt.zero) then
      bb=a2**2
      rmat(1,2)=rmat(1,2)*a2
      rmat(2,2)=rmat(2,2)*a2
      rmat(1,3)=rmat(1,3)*bb
      rmat(2,3)=rmat(2,3)*bb
      rmat(1,5)=rmat(1,5)*a2
      rmat(2,5)=rmat(2,5)*a2
      rinv(1,2)=rinv(1,2)/a2
      rinv(2,2)=rinv(2,2)/a2
      rinv(1,3)=rinv(1,3)/bb
      rinv(2,3)=rinv(2,3)/bb
      rinv(1,5)=rinv(1,5)/a2
      rinv(2,5)=rinv(2,5)/a2
   endif

   if (a3.gt.zero) then
      bb=a3**2
      rmat(1,4)=rmat(1,4)*a3
      rmat(2,4)=rmat(2,4)*a3
      rmat(1,5)=rmat(1,5)*a3
      rmat(2,5)=rmat(2,5)*a3
      rmat(1,6)=rmat(1,6)*bb
      rmat(2,6)=rmat(2,6)*bb
      rinv(1,4)=rinv(1,4)/a3
      rinv(2,4)=rinv(2,4)/a3
      rinv(1,5)=rinv(1,5)/a3
      rinv(2,5)=rinv(2,5)/a3
      rinv(1,6)=rinv(1,6)/bb
      rinv(2,6)=rinv(2,6)/bb
   endif

   return
   end subroutine unscale3

   subroutine yfour(rmat,rinv,nchan,ier)
   !-------------------------------------------------------------------
   ! Invert matrix for four or more channels.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso,nsyse
   ! externals
   real(kr)::rmat(2,*),rinv(2,*)
   integer::nchan,ier
   ! internals
   integer::kj,k,j,info
   integer,dimension(:),allocatable::iii
   real(kr),dimension(:,:),allocatable::dummy

   allocate(iii(nchan))
   call xspfa(rmat,nchan,iii,info)
   if (info.ne.0) write(nsyse,'('' Problem in xspfa with info='', I5)') info
   if (info.ne.0) write(nsyso,'('' Problem in xspfa with info='', I5)') info
   kj=0
   do k=1,nchan
      do j=1,k
         kj=kj+1
      enddo
   enddo
   allocate(dummy(2,kj))
   kj=0
   do k=1,nchan
      do j=1,nchan
         dummy(1,j)=0
         dummy(2,j)=0
      enddo
      dummy(1,k)=1
      call xspsl(rmat,nchan,iii,dummy)
      do j=1,k
         kj=kj+1
         rinv(1,kj)=dummy(1,j)
         rinv(2,kj)=dummy(2,j)
      enddo
   enddo
   deallocate(dummy)
   deallocate(iii)
   return
   end subroutine yfour

   subroutine xspfa(ap,n,kpvt,info)
   !-------------------------------------------------------------------
   !
   !   FACTORS A COMPLEX SYMMETRIC MATRIX STORED IN
   !   PACKED FORM BY ELIMINATION WITH SYMMETRIC PIVOTING.
   !
   !   TO SOLVE  A*X = B , FOLLOW Xspfa BY Xspsl.
   !   TO COMPUTE  INVERSE(A)*C , FOLLOW Xspfa BY Xspsl.
   !   TO COMPUTE  DETERMINANT(A) , FOLLOW Xspfa BY Xspdi.
   !   TO COMPUTE  INVERSE(A) , FOLLOW Xspfa BY Xspdi... which I don't have
   !
   !   ON ENTRY
   !
   !      Ap      REAL*8 (2,(N*(N+1)/2))
   !              Ap(1,*) is real part, Ap(2,*) is imaginary part.
   !              THE PACKED FORM OF A SYMMETRIC MATRIX  A .  THE
   !              COLUMNS OF THE UPPER TRIangle ARE STORED SEQUENTIALLY
   !              IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 .
   !              SEE COMMENTS BELOW FOR DETAILS.
   !
   !      N       INTEGER
   !              THE ORDER OF THE MATRIX  A .
   !
   !   OUTPUT
   !
   !      Ap      A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHIch
   !              WERE USED TO OBTAIN IT STORED IN PACKED FORM.
   !              THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)
   !              WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT
   !              UPPER TRIANGULAR MATRICES , TRANS(U) IS THE
   !              TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL
   !              WITH 1 BY 1 AND 2 BY 2 BLOCKS.
   !
   !      Kpvt    INTEGER(N)
   !              AN INTEGER VECTOR OF PIVOT INDICES.
   !
   !      Info    INTEGER
   !              = 0  NORMAL VALUE.
   !              = K  IF THE K-TH PIVOT BLOCK IS SINGULAR. THIS IS
   !                   Not AN ERROR CONDITION FOR THIS SBROUTINE,
   !                   BUT IT DOES INDICATE THAT SSPSL OR SSPDI MAY
   !                   DIVIDE BY ZERO IF CALLED.
   !
   !   PACKED STORAGE
   !
   !        THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER
   !        TRIangle OF A SYMMETRIC MATRIX.
   !
   !              K = 0
   !              DO J=1,N
   !                 DO I=1,J
   !                    K = K + 1
   !                    Ap(K) = A(I,J)
   !                 END DO
   !              END DO
   !
   !-------------------------------------------------------------------
   ! externals
   integer::n,info
   integer::kpvt(*)
   real(kr)::ap(2,*)
   ! internals
   real(kr)::ak,aki,akm1,akm1i,bk,bki,bkm1,bkm1i,denom,denomi,&
     dmulk,dmulki,dmlkm1,dmlkmi,t,ti
   real(kr)::absakk,alpha,colmax,rowmax,aa,dd,xx,xxi
   integer::ij,ik,ikm1,im,imax,imaxp1,imim,imj,imk
   integer::j,jj,jk,jkm1,jmax,jmim,k,kk,km1,km1k,km1km1,km2,kstep
   logical::swap
   real(kr),parameter::zero=0

   !--initialize

   ! alpha is used in choosing pivot block size.
   alpha=(1+sqrt(17.0e0_kr))/8

   info=0
   im=0

   !--main loop on k, which goes from n to 1.

   k=n
   ik=(n*(n-1))/2
   do

      !--leave the loop if k=0 or k=1.
      if (k.eq.0) exit
      if (k.le.1) then
         kpvt(1)=1
         if (ap(1,1).eq.zero.and.ap(2,1).eq.zero) info=1
         exit
      endif

      ! this section of code determines the kind of
      ! elimination to be performed.  when it is completed,
      ! kstep will be set to the size of the pivot block, and
      ! swap will be set to .true. if an interchange is
      ! required.

      km1=k-1
      kk=ik+k
      absakk=ap(1,kk)**2+ap(2,kk)**2

      ! determine the largest off-diagonal element in column k.

      imax=ixamax(k-1,ap(1,ik+1),1)
      imk=ik+imax
      colmax=ap(1,imk)**2+ap(2,imk)**2

      if (absakk.ge.alpha*colmax) then

         kstep=1
         swap=.false.

      else

         ! determine the largest off-diagonal element in rov imax.

         rowmax=0
         imaxp1=imax+1
         im=(imax*(imax-1))/2
         imj=im+2*imax
         do j=imaxp1,k
            aa=ap(1,imj)**2+ap(2,imj)**2
            rowmax=max(rowmax,aa)
            imj=imj+j
         enddo
         if (imax.ne.1) then
            jmax=ixamax(imax-1,ap(1,im+1),1)
            jmim=jmax+im
            aa=ap(1,jmim)**2+ap(2,jmim)**2
            rowmax=max(rowmax,aa)
         endif
         imim=imax+im
         aa=ap(1,imim)**2+ap(2,imim)**2
         if (aa.ge.alpha*rowmax) then
            kstep=1
            swap=.true.
         else
            if (absakk.ge.alpha*colmax*(colmax/Rowmax)) then
               kstep=1
               swap=.false.
            else
               kstep=2
               swap=imax.ne.km1
            endif
         endif

      endif

      if (max(absakk,colmax).eq.zero) then

         ! column k is zero.  set info and iterate the loop.

         kpvt(k)=k
         info=k

      else

         if (kstep.ne.2) then

            ! 1 x 1 pivot block.
            if (swap) then
               ! perform an interchange.
               call xswap(imax,ap(1,im+1),1,ap(1,ik+1),1)
               imj=ik+imax
               do jj=imax,k
                  j=k+imax-jj
                  jk=ik+j
                  t=ap(1,jk)
                  ap(1,jk)=ap(1,imj)
                  ap(1,imj)=t
                  t=ap(2,jk)
                  ap(2,jk)=ap(2,imj)
                  ap(2,imj)=t
                  imj=imj-(j-1)
               enddo
            endif

            ! perform the elimination.
            ij=ik-(k-1)
            do jj=1,km1
               j=k-jj
               jk=ik+j
               aa=ap(1,kk)**2+ap(2,kk)**2
               dmulk =-(ap(1,jk)*ap(1,kk)+ap(2,jk)*ap(2,kk))/aa
               dmulki=(ap(1,jk)*ap(2,kk)-ap(2,jk)*ap(1,kk))/aa
               t =dmulk
               ti=dmulki
               call xaxpy(j,t,ti,ap(1,ik+1),1,ap(1,ij+1),1)
               ap(1,jk)=dmulk
               ap(2,jk)=dmulki
               ij=ij-(j-1)
            enddo

            ! set the pivot array.

            kpvt(k)=k
            if (swap) kpvt(k)=imax

         else

            ! 2 x 2 pivot block.
            km1k=ik+k-1
            ikm1=ik-(k-1)
            if (swap) then

               ! perform an interchange.
               call xswap(imax,ap(1,im+1),1,ap(1,ikm1+1),1)
               imj=ikm1+imax
               do jj=imax,km1
                  j=km1+imax-jj
                  jkm1=ikm1+j
                  t=ap(1,jkm1)
                  ap(1,jkm1)=ap(1,imj)
                  ap(1,imj)=t
                  t=ap(2,jkm1)
                  ap(2,jkm1)=ap(2,imj)
                  ap(2,imj)=t
                  imj=imj-(j-1)
               enddo
               t=ap(1,km1k)
               ap(1,km1k)=ap(1,imk)
               ap(1,imk)=t
               t=ap(2,km1k)
               ap(2,km1k)=ap(2,imk)
               ap(2,imk)=t
            endif

            ! perform the elimination.
            km2=k-2
            if (km2.ne.0) then
               aa=ap(1,km1k)**2+ap(2,km1k)**2
               ak=(ap(1,kk)*ap(1,km1k)+ap(2,kk)*ap(2,km1k))/aa
               aki=(ap(2,kk)*ap(1,km1k)-ap(1,kk)*ap(2,km1k))/aa
               km1km1=ikm1+k-1
               akm1=(ap(1,km1km1)*ap(1,km1k)+ap(2,km1km1)*ap(2,km1k))/aa
               akm1i=(ap(2,km1km1)*ap(1,km1k)-ap(1,km1km1)*ap(2,km1k))/aa
               denom =1-(ak*akm1-aki*akm1i)
               denomi=-(ak*akm1i+aki*akm1)
               dd=denom**2+denomi**2
               ij=ik-(k-1)-(k-2)
               dO jj=1,km2
                  j=km1-jj
                  jk=ik+j
                  bk=(ap(1,jk)*ap(1,km1k)+ap(2,jk)*ap(2,km1k))/aa
                  bki=(ap(2,jk)*ap(1,km1k)-ap(1,jk)*ap(2,km1k))/aa
                  jkm1=ikm1+j
                  bkm1=(ap(1,jkm1)*ap(1,km1k)+ap(2,jkm1)*ap(2,km1k))/aa
                  bkm1i=(ap(2,jkm1)*ap(1,km1k)-ap(1,jkm1)*ap(2,km1k))/aa
                  xx=akm1*bk-akm1i*bki-bkm1
                  xxi=akm1*bki+akm1i*bk-bkm1i
                  dmulk=(xx*denom+xxi*denomi)/dd
                  dmulki=(xxi*denom-xx*denomi)/dd
                  xx=ak*bkm1-aki*bkm1i-bk
                  xxi=ak*bkm1i+aki*bkm1-bki
                  dmlkm1=(xx*denom+xxi*denomi)/dd
                  dmlkmi=(xxi*denom-xx*denomi)/dd
                  t=dmulk
                  ti=dmulki
                  call xaxpy(j,t,ti,ap(1,ik+1),1,ap(1,ij+1),1)
                  t=dmlkm1
                  ti=dmlkmi
                  call xaxpy(j,t,ti,ap(1,ikm1+1),1,ap(1,ij+1),1)
                  ap(1,jk)=dmulk
                  ap(2,jk)=dmulki
                  ap(1,jkm1)=dmlkm1
                  ap(2,jkm1)=dmlkmi
                  ij=ij-(j-1)
               enddo
            endif

            ! set the pivot array.
            kpvt(k)=1-k
            if (swap) kpvt(k)=-imax
            kpvt(k-1)=kpvt(k)
         endif
      endif
      ik=ik-(k-1)
      if (kstep.eq.2) ik=ik-(k-2)
      k=k-kstep
   enddo

   return
   end subroutine xspfa

   integer function ixamax(n,sx,incx)
   !-------------------------------------------------------------------
   ! Finds the index of element having maximum squared value.
   !-------------------------------------------------------------------
   ! externals
   integer::n,incx
   real(kr)::sx(2,*)
   ! internals
   real(kr)::smax,aa
   integer i,ix

   ixamax=0
   if (n.lt.1) return
   ixamax=1
   if (n.eq.1) return

   if (incx.ne.1) then

      !--code for increment not equal to 1
      ix=1
      smax=sx(1,1)**2+sx(2,1)**2
      ix=ix+incx
      do i=2,n
         aa=sx(1,ix)**2+sx(2,ix)**2
         if (aa.gt.smax) then
            ixamax=i
            smax=aa
         endif
         ix=ix+incx
      enddo

   else

      !--code for increment equal to 1
      smax=sx(1,1)**2+Sx(2,1)**2
      do i=2,n
         aa=sx(1,i)**2+sx(2,i)**2
         if (aa.gt.smax) then
            ixamax=i
            smax=aa
         endif
      enddo

   endif

   return
   end function ixamax

   subroutine xaxpy(n,sa,sai,sx,incx,sy,incy)
   !-------------------------------------------------------------------
   ! Complex constant times a vector plus another vector.
   ! uses unrolled loop for increments equal to one.
   !-------------------------------------------------------------------
   ! externals
   integer::n,incx,incy
   real(kr)::sa,sai,sx(2,*),sy(2,*)
   ! internals
   real(kr)::aa
   integer::i,ix,iy,m,mp1
   real(kr),parameter::zero=0

      if (n.le.0) return
      if (sa.eq.zero.and.sai.eq.zero) return

      if (incx.ne.1.or.incy.ne.1) then

         !--code for unequal increments or equal increments
         !--not equal to 1
         ix=1
         iy=1
         if (incx.lt.0) ix=(-n+1)*incx+1
         if (incy.lt.0) iy=(-n+1)*incy+1
         do i=1,n
            aa=sy(1,iy)+sa*sx(1,ix)-sai*sx(2,ix)
            sy(2,iy)=sy(2,iy)+sa*sx(2,ix)+sai*sx(1,ix)
            sy(1,iy)=aa
            ix=ix+incx
            iy=iy+incy
         enddo

      else

         !--code for both increments equal to 1

         !--clean-up loop
         m=mod(n, 4)
         if (m.ne.0) then
            do i=1,m
               aa=sy(1,i)+sa*sx(1,i)-sai*sx(2,i)
               sy(2,i)=sy(2,i)+sa*sx(2,i)+sai*sx(1,i)
               sy(1,i)=aa
            enddo
            if (n.lt.4) return
         endif
         mp1=m+1
         do i=mp1,n,4
            aa       =sy(1,i  )+sa*sx(1,i  )-sai*sx(2,i  )
            sy(2,i  )=sy(2,i  )+sa*sx(2,i  )+sai*sx(1,i  )
            sy(1,i  )=aa
            aa       =sy(1,i+1)+sa*sx(1,i+1)-sai*sx(2,i+1)
            sy(2,i+1)=sy(2,i+1)+sa*sx(2,i+1)+sai*sx(1,i+1)
            sy(1,i+1)=aa
            aa       =sy(1,i+2)+sa*sx(1,i+2)-sai*sx(2,i+2)
            sy(2,i+2)=sy(2,i+2)+sa*sx(2,i+2)+sai*sx(1,i+2)
            sy(1,i+2)=aa
            aa       =sy(1,i+3)+sa*sx(1,i+3)-sai*sx(2,i+3)
            sy(2,i+3)=sy(2,i+3)+sa*sx(2,i+3)+sai*sx(1,i+3)
            sy(1,i+3)=aa
         enddo

      endif

   return
   end subroutine xaxpy

   real(kr) function xdot(xdoti,n,sx,incx,sy,incy)
   !-------------------------------------------------------------------
   ! forms the dot product of two vectors.
   ! uses unrolled loops for increments equal to one.
   !-------------------------------------------------------------------
   ! externals
   real(kr)::xdoti,sx(2,*),sy(2,*)
   integer::n,incx,incy
   ! internals
   real(kr)::stemp,stempi
   integer::i,ix,iy,m,mp1

   stemp=0
   stempi=0
   xdot =0
   xdoti=0
   if (n.le.0) return
   if (incx.ne.1.or.incy.ne.1) then

      !--code for unequal increments or equal increments not equal to 1
      ix=1
      iy=1
      if (incx.lt.0) ix=(-n+1)*incx+1
      if (incy.lt.0) iy=(-n+1)*incy+1
      do i=1,n
         stemp =stemp +sx(1,ix)*sy(1,iy)-sx(2,ix)*sy(2,iy)
         stempi=stempi+sx(2,ix)*sy(1,iy)+sx(1,ix)*sy(2,iy)
         ix=ix+incx
         iy=iy+incy
      enddo
      xdot=stemp
      xdoti=stempi

   else

      !--code for both increments equal to 1

      !--cleean up loop
      m=mod(n,5)
      if (m.ne.0) then
         do i=1,m
            stemp =stemp +sx(1,i)*sy(1,i)-sx(2,i)*sy(2,i)
            stempi=stempi+sx(2,i)*sy(1,i)+sx(1,i)*sy(2,i)
         enddo
         if (n.lt.5) then
            xdot =stemp
            xdoti=stempi
            return
         endif
      endif
      mp1=m+1
      do i=mp1,n,5
         stemp =  stemp +sx(1,i  )*sy(1,i  )-sx(2,i  )*sy(2,i  )&
                        +sx(1,i+1)*sy(1,i+1)-sx(2,i+1)*sy(2,i+1)&
                        +sx(1,i+2)*sy(1,i+2)-sx(2,i+2)*sy(2,i+2)&
                        +sx(1,i+3)*sy(1,i+3)-sx(2,i+3)*sy(2,i+3)&
                        +sx(1,i+4)*sy(1,i+4)-sx(2,i+4)*sy(2,i+4)
         stempi = stempi+sx(2,i  )*sy(1,i  )+sx(1,i  )*sy(2,i  )&
                        +sx(2,i+1)*sy(1,i+1)+sx(1,i+1)*sy(2,i+1)&
                        +sx(2,i+2)*sy(1,i+2)+sx(1,i+2)*sy(2,i+2)&
                        +sx(2,i+3)*sy(1,i+3)+sx(1,i+3)*sy(2,i+3)&
                        +sx(2,i+4)*sy(1,i+4)+sx(1,i+4)*sy(2,i+4)
      enddo
      xdot=stemp
      xdoti=stempi

   endif

   return
   end function xdot

   subroutine xswap(n,sx,incx,sy,incy)
   !-------------------------------------------------------------------
   ! interchanges two vectors.
   ! uses unrolled loops for increments equal to 1.
   !-------------------------------------------------------------------
   ! externals
   integer::n,incx,incy
   real(kr)::sx(2,*),sy(2,*)
   ! internals
   real(kr)::stemp
   integer::i,ix,iy,m,mp1

   if (n.le.0) return
   if (incx.ne.1.or.incy.ne.1) then

      !--code for unequal increments or equal increments not equal TO 1
      ix=1
      iy=1
      if (incx.lt.0) ix=(-n+1)*incx+1
      if (incy.lt.0) iy=(-n+1)*incy+1
      do i=1,n
         stemp   =sx(1,ix)
         sx(1,ix)=sy(1,iy)
         sy(1,iy)=stemp
         stemp   =sx(2,ix)
         sx(2,ix)=sy(2,iy)
         sy(2,iy)=stemp
         ix=ix+incx
         iy=iy+incy
      enddo

   else

      !--code for both increments equal to 1

      !--clean-up loop
      m=mod(n,3)
      if (m.ne.0) then
         do i=1,m
            stemp  =sx(1,i)
            sx(1,i)=sy(1,i)
            sy(1,i)=stemp
            stemp  =sx(2,i)
            sx(2,i)=sy(2,i)
            sy(2,i)=stemp
         enddo
         if (n.lt.3) return
      endif
      mp1=m+1
      do i=mp1,n,3
         stemp    =sx(1,i  )
         sx(1,i  )=sy(1,i  )
         sy(1,i  )=stemp
         stemp    =sx(2,i  )
         sx(2,i  )=sy(2,i  )
         sy(2,i  )=stemp
         stemp    =sx(1,i+1)
         sx(1,i+1)=sy(1,i+1)
         sy(1,i+1)=stemp
         stemp    =sx(2,i+1)
         sx(2,i+1)=sy(2,i+1)
         sy(2,i+1)=stemp
         stemp    =sx(1,i+2)
         sx(1,i+2)=sy(1,i+2)
         sy(1,i+2)=stemp
         stemp    =sx(2,i+2)
         sx(2,i+2)=sy(2,i+2)
         sy(2,i+2)=Stemp
      enddo

   endif

   return
   end subroutine xswap

   subroutine xspsl(ap,n,kpvt,b)
   !-------------------------------------------------------------------
   !
   !   Xspsl SOLVES THE complex SYMMETRIC SYSTEM
   !               A * X = B
   !   USING THE FACTORS COMPUTED BY Xspfa.
   !
   !   ON ENTRY
   !
   !      Ap      REAL*8 (2,N*(N+1)/2)
   !              THE OUTPUT FROM SSPFA.
   !
   !      N       INTEGER
   !              THE ORDER OF THE MATRIX  A .
   !
   !      Kpvt    INTEGER(N)
   !              THE PIVOT VECTOR FROM SSPFA.
   !
   !      B       REAL(2,N)
   !              THE RIGHT HAND SIDE VECTOR.
   !
   !   ON RETURN
   !
   !      B       THE SOLUTION VECTOR  X .
   !
   !   ERROR CONDITION
   !
   !      A Division BY ZERO MAY OCCUR IF  Xspco  HAS SET RCOND .EQ. 0.0
   !      OR  Xspfa  HAS SET Info .NE. 0  .
   !
   !   TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
   !   WITH  P  COLUMNS
   !         CALL Xspfa (Ap, N, Kpvt, Info)
   !         IF (Info .NE. 0) GO TO ...
   !         DO J=1,P
   !            CALL Xspsl (Ap, N, Kpvt, C(1,J))
   !         END DO
   !
   !   SBROUTINES AND FNCTIONS
   !
   !   BLAS Xaxpy, Xdot
   !   FORTRAN IABS
   !
   !-------------------------------------------------------------------
   ! externals
   real(kr)::ap(2,*),b(2,*)
   integer::n,kpvt(*)
   ! internals
   real(kr)::ak,akm1,bk,bkm1,denom,temp
   real(kr)::aki,akm1i,bki,bkm1i,xdoti,denomi,aa,ab,abi,dd
   integer::ik,ikm1,ikp1,k,kk,km1k,km1km1,kp

   !--loop backward applying the transformations and d inverse to b.

   k=n
   ik=(n*(n-1))/2
  10 continue
   if (k.ne.0) then
      kk=ik+k

      if (kpvt(k).GE.0) then

         ! 1 x 1 pivot block.
         if (k.ne.1) then
            kp=kpvt(k)
            if (kp.ne.k) then
               ! interchange.
               temp   =b(1,k)
               b(1,k) =b(1,kp)
               b(1,kp)=temp
               temp   =b(2,k)
               b(2,k) =b(2,kp)
               b(2,kp)=temp
            endif
            ! apply the transformation.
            call xaxpy(k-1,b(1,k),b(2,k),ap(1,ik+1),1,b(1,1),1)
         endif
         ! apply d inverse.
         aa=ap(1,kk)**2+ap(2,kk)**2
         ab    =(b(1,k)*ap(1,kk)+b(2,k)*ap(2,kk))/aa
         b(2,k)=(b(2,k)*ap(1,kk)-b(1,k)*ap(2,kk))/aa
         b(1,k)=ab
         k=k-1
         ik=ik-k

      else

         !  2 x 2 pivot block.
         ikm1=ik-(k-1)
         if (k.ne.2) then
            kp=iabs(kpvt(k))
            if (kp.ne.k-1) then
               ! interchange.
               temp    =b(1,k-1)
               b(1,k-1)=b(1,kp)
               b(1,kp) =temp
               temp    =b(2,k-1)
               b(2,k-1)=b(2,kp)
               b(2,kp) =temp
            endif
            ! apply the transformation.
            call xaxpy(k-2,b(1,k),b(2,k),ap(1,ik+1),1,b(1,1),1)
            call xaxpy(k-2,b(1,k-1),b(2,k-1),ap(1,ikm1+1),1,b(1,1), 1)
         endif
         ! apply d inverse.
         km1k=ik+k-1
         kk=ik+k
         aa=ap(1,km1k)**2+ap(2,km1k)**2
         ak =(ap(1,kk)*ap(1,km1k)+ap(2,kk)*ap(2,km1k))/aa
         aki=(ap(2,kk)*ap(1,km1k)-ap(1,kk)*ap(2,km1k))/aa
         km1km1=ikm1+k-1
         akm1 =(ap(1,km1km1)*ap(1,km1k)+ap(2,km1km1)*ap(2,km1k))/aa
         akm1i=(ap(2,km1km1)*ap(1,km1k)-ap(1,km1km1)*ap(2,km1k))/aa
         bk =(b(1,k)*ap(1,km1k)+b(2,k)*ap(2,km1k))/aa
         bki=(b(2,k)*ap(1,km1k)-b(1,k)*ap(2,km1k))/aa
         bkm1 =(b(1,k-1)*ap(1,km1k)+b(2,k-1)*ap(2,km1k))/aa
         bkm1i=(b(2,k-1)*ap(1,km1k)-b(1,k-1)*ap(2,km1k))/aa
         denom =ak*akm1-aki*akm1i-1
         denomi=ak*akm1i+aki*akm1
         dd=denom**2+denomi**2
         ab =akm1*bk-akm1i*bki-bkm1
         abi=akm1i*bk+akm1*bki-bkm1i
         b(1,k)=(ab*denom+abi*denomi)/dd
         b(2,k)=(abi*denom-ab*denomi)/dd
         ab =ak*bkm1-aki*bkm1i-bk
         abi=aki*bkm1+ak*bkm1i-bki
         b(1,k-1)=(ab*denom+abi*denomi)/dd
         b(2,k-1)=(abi*denom-ab*denomi)/dd
         k=k-2
         ik=ik-(k+1)-k

      endif

      go to 10
   endif

   !--loop forward applying the transformations.
   k=1
   ik=0
  90 continue
   if (k.LE.n) then
      if (kpvt(k).ge.0) then

         ! 1 x 1 pivot block.
         if (k.ne.1) then
            ! apply the transformation.
            b(1,k)=b(1,k)+xdot(xdoti,k-1,ap(1,ik+1),1,b(1,1),1)
            b(2,k)=b(2,k)+xdoti
            kp=kpvt(k)
            if (kp.ne.k) then
               ! interchange.
               temp   =b(1,k )
               b(1,k )=b(1,kp)
               b(1,kp)=temp
               temp   =b(2,k )
               b(2,k )=b(2,kp)
               b(2,kp)=temp
            endif
         endif
         ik=ik+k
         k=k+1

      else

         ! 2 x 2 pivot block.
         if (k.ne.1) then
            ! apply the transformation.
            b(1,k  )=b(1,k  )+xdot(xdoti,k-1,ap(1,ik+1),1,b(1,1),1)
            b(2,k  )=b(2,k  )+xdoti
            ikp1=ik+k
            b(1,k+1)=b(1,k+1)+xdot(xdoti,k-1,ap(1,ikp1+1),1,b(1,1),1)
            b(2,k+1)=b(2,k+1)+xdoti
            kp=iabs(kpvt(k))
            if (kp.ne.k) then
               ! interchange.
               temp   =b(1,k )
               b(1,k )=b(1,kp)
               b(1,kp)=temp
               temp   =b(2,k )
               b(2,k )=b(2,kp)
               b(2,kp)=temp
            endif
         endif
         ik=ik+k+k+1
         k =k +2

         endif
         go to 90
      endif

   return
   end subroutine xspsl

   subroutine setxqx(nchan,ier)
   !-------------------------------------------------------------------
   ! Form XQ & XXXX matrices, where
   !   XQ   = Yinv * Rmat       and
   !   XXXX = P/L + sqrt(P)/L         (1/L-R)**-1          sqrt(P)/L
   !        =       sqrt(P)/(S-B+IP) * Yinv       * Rmat * sqrt(P)
   !        =       sqrt(P)/L        * XQ                * sqrt(P)
   !
   !    Note that the matrix W defined in SAMMY manual is given
   !    by W(c,c') = delta(c,c') + 2i XXXX(c,c')
   !    as in Eq. (III.D.4) in SAMMY manual R3
   !
   !    ie W    = I + 2i XXXX
   !-------------------------------------------------------------------
   ! externals
   integer::nchan,ier
   ! internals
   integer::i,j,ij,k,jk
   real(kr)::plr,pli

   do i=1,nchan
      do j=1,nchan
         xqr(i,j,ier)=0
         xqi(i,j,ier)=0
      enddo
   enddo

   !--Xqr(k,i) = (L**-1-R)**-1 * R ... note asymmetry
   do i=1,nchan
      do j=1,nchan
         ij=ijkl(j,i)
         do k=1,nchan
            jk=ijkl(k,j)
            xqr(k,i,ier)=xqr(k,i,ier)+yinv(1,ij,ier)*rmat(1,jk,ier)&
              -yinv(2,ij,ier)*rmat(2,jk,ier)
            xqi(k,i,ier)=xqi(k,i,ier)+yinv(1,ij,ier)*rmat(2,jk,ier)&
              +yinv(2,ij,ier)*rmat(1,jk,ier)
         enddo
      enddo
   enddo

   !--Xxxx = sqrt(P)/L  * xq * sqrt(P) ... symmetric
   ij=0
   do i=1,nchan
      plr=rootp(i,ier)*elinvr(i,ier)
      pli=rootp(i,ier)*elinvi(i,ier)
      do j=1,i
         ij=ij+1
         xxxxr(ij,ier)=rootp(j,ier)*(xqr(j,i,ier)*plr-xqi(j,i,ier)*pli)
         xxxxi(ij,ier)=rootp(j,ier)*(xqi(j,i,ier)*plr+xqr(j,i,ier)*pli)
      enddo
   enddo
   return
   end subroutine setxqx

   subroutine sectio(n,nent,next,nchan,ier)
   !-------------------------------------------------------------------
   ! Generate pieces of cross sections (except for "4 pi/E")
   !-------------------------------------------------------------------
   ! externals
   integer n,ier
   integer::nent,next,nchan
   ! internals
   integer::jj,ii,ij,i,j,ipx,ichan,ichanx,ip
   real(kr)::zz,termn,ar,ai,br,bi,cr,ci,dr,di,terma
   real(kr),parameter::zero=0

   do jj=1,npp
      crss(jj,ier)=0
   enddo

   ! Entrance channel, ipp=2
   !    elastic  crss(1) = g*0.25* sum(entrance chs c,c')
   !                                     times |(1-U(c,c'))| **2 / zz
   !                     = g* [ sin(phi)**2 * (1-2*xxxXi)
   !                            - sin(2phi)*xxxxr
   !                            + (xxxxr**2 + xxxxi**2) ] / zz

   ii=0
   ij=0
   do i=1,nent
      zz=zke(i,n,ier)**2
      ii=ii+i
      termn=sinsqr(i,ier)*(1-2*xxxxi(ii,ier))-sin2ph(i,ier)*xxxxr(ii,ier)
      termn=termn/zz
      do j=1,i
         ij=ij+1
         ar=(xxxxr(ij,ier)**2+xxxxi(ij,ier)**2)/zz
         if (i.ne.j) ar=ar+ar
         termn=Termn+ar
      enddo
      crss(1,ier)=termn+crss(1,ier)
   enddo
   crss(1,ier)=crss(1,ier)*goj(n,ier)

   ! Ipp=1 term, sort-of
   !    absorption = g*0.25 * sum(inc c)
   !                  [ 1 -  sum(inc c') |U(c,c')| **2 ] / zz
   !               = - g* (xxxxr**2 + xxxxi**2) / zz

   ii=0
   ij=0
   do i=1,nent
      ii=ii+i
      zz=zke(i,n,ier)**2
      terma=xxxxi(ii,ier)/zz
      do j=1,i
         ij=ij+1
         ar=(-xxxxr(ij,ier)**2-xxxxi(ij,ier)**2)/zz
         if (i.ne.j) ar=ar+ar
         terma=terma+ar
      enddo
      crss(2,ier)=terma+crss(2,ier)
   enddo
   crss(2,ier)=crss(2,ier)*goj(n,ier)

   ! All other channels, classed by particle-pair number
   !    reaction ch c'= g*0.25 * sum(inc c) |U(c,c')|**2 / zz
   !                  = g* (xxxxR**2 + xxxxI**2) / zz

   do jj=1,next
      j=jj+nent
      if (j.le.nchan) then
         ip=ipp(j,n,ier)
         do i=1,nent
            zz=zke(i,n,ier)**2
            ij=(j*(j-1))/2+i
            ! ij = ijkl(i,j) but i < j always
            crss(ip,ier)=crss(ip,ier)+(xxxxr(ij,ier)**2+xxxxi(ij,ier)**2)/zz
         enddo
      endif
   enddo
   if (npp.gt.2) then
      do ip=3,npp
         crss(ip,ier)=crss(ip,ier)*goj(n,ier)
      enddo
   endif

   if (Want_Angular_Dist) then

   !    Angular Distribution   Crssx(.,i,ix,ip) = (1-U)/2
   !    Note that we do not calculate the absorption anglar distribution

      do ip=2,npp
         ipx=ip
         if (ip.eq.2) ipx=1
         do ichan=1,nchan
            if (ipp(ichan,n,ier).eq.ip) then

               if (zeta(ichan,n,ier).ne.zero.and.su.gt.echan(ichan,n,ier)) then
                  call gcphase(cr,ci,lspin(ichan,n,ier),&
                     echan(ichan,n,ier),zeta(ichan,n,ier),su)
               else
                  cr=1
                  ci=0
               endif
               ii=(ichan*(ichan-1))/2
               do ichanx=1,ichan
                  ii=ii+1
                  if (ichanx.le.nent) then
                     ! real and imaginary parts of (1-U)/2
                     if (ichanx.eq.ichan) then
                        ar=sinsqr(ichan,ier)*(1-2*xxxxi(ii,ier))&
                          -sin2ph(ichan,ier)*xxxxr(ii,ier)+xxxxi(ii,ier)
                        ai=sin2ph(ichan,ier)*(0.5e0_kr-xxxxi(ii,ier))&
                          -(1-2*sinsqr(ichan,ier))*xxxxr(ii,ier)
                     else
                        ar=cscs(1,ii,ier)*xxxxi(ii,ier)&
                          -cscs(2,ii,ier)*xxxxr(ii,ier)
                        ai=-cscs(1,ii,ier)*xxxxr(ii,ier)&
                          -cscs(2,ii,ier)*xxxxi(ii,ier)
                     endif
                     if (zeta(ichan,n,ier).ne.zero.or.&
                      zeta(ichanx,n,ier).ne.zero) then
                        br=ar*cr-ai*ci
                        bi=ar*ci+ai*cr
                        if ((lspin(ichanx,n,ier).ne.lspin(ichan,n,ier).or.&
                          zeta(ichanx,n,ier).ne.zeta(ichan,n,ier)).and.&
                          ichan.ne.ichanx) then
                           if (zeta(ichanx,n,ier).ne.zero.and.&
                             su.gt.echan(ichanx,n,ier)) then
                              call gcphase(dr,di,&
                                lspin(ichanx,n,ier),echan(ichanx,n,ier),&
                                zeta(ichanx,n,ier),su)
                           else
                              dr=1
                              di=0
                           endif
                        else
                           dr=cr
                           di=ci
                        endif
                        ar=br*dr-bi*di
                        ai=br*di+bi*dr
                     endif
                     crssx(1,ichanx,ichan,ipx,n,ier)=&
                       ar+crssx(1,ichanx,ichan,ipx,n,ier)
                     crssx(2,ichanx,ichan,ipx,n,ier)=&
                       ai+crssx(2,ichanx,ichan,ipx,n,ier)
                     !   Note that ichanx = incident, ichan = outgoing
                     !        (so ichan may be entrance channel too)
                  endif
               enddo
            endif
         enddo
      enddo
   endif

   return
   end subroutine sectio

   integer function ijkl(m,n)
   !-------------------------------------------------------------------
   ! external
   integer::m,n

   if (m.le.n) then
      ijkl=(n*(n-1))/2+m
   else
      ijkl=(m*(m-1))/2+n
   endif

   return
   end function ijkl

   subroutine gcphase(cr,ci,lspin,echan,zeta,su)
   !-------------------------------------------------------------------
   ! Get Coulomb phase.
   !-------------------------------------------------------------------
   ! externals
   real(kr)::cr,ci,echan,zeta,su
   integer::lspin
   ! internals
   real(kr)::eta,aa,cc,ss,ccx,ssx,ccy,ssy
   integer::l,ll
   real(kr),parameter::zero=0

   eta=zeta/sqrt(su-echan)
   l=lspin
   if (l.eq.0) then
      cr=1
      ci=0
   else
      aa=eta**2+1
      aa=sqrt(aa)
      cc=1/aa
      ss=eta/aa
      if (l.gt.1) then
         do ll=2,l
            aa=eta**2+ll**2
            aa=sqrt(aa)
            ccx=ll/aa
            ssx=eta/aa
            ccy=cc*ccx-ss*ssx
            ssy=cc*ssx+ss*ccx
            cc=ccy
            ss=ssy
         enddo
      endif
      cr=cc
      ci=ss
      !  (cr,ci) = (Re,Im) exp{
      !     i sum from (ll=1 to l) tan^{-1} (eta/ll) }
   endif

   return
   end subroutine gcphase

   subroutine setqri(nchan,nn,ier)
   !-------------------------------------------------------------------
   ! Generate qr,qi= sqrt(p)/(s-b+ip) * yinv*yinv * sqrt(p)/(s-b+ip)
   ! That is, qr(kl,ij) is (real part of) partial of xxxx(kl) wrt r(kj)
   !-------------------------------------------------------------------
   ! externals
   integer::nchan,nn,ier
   ! internals
   integer::i,k,ik,ij,j,kl,l
   real(kr)::plri,plii
   real(kr),parameter::zero=0

   !--redefine meaning of xqr & xqi
   !--xq = rootp*yinv
   do i=1,nchan
      plri=rootp(i,ier)*elinvr(i,ier)
      plii=rootp(i,ier)*elinvi(i,ier)
      do k=1,nchan
         ik=ijkl(i,k)
         xqr(k,i,ier)=plri*yinv(1,ik,ier)-plii*yinv(2,ik,ier)
         xqi(k,i,ier)=plri*yinv(2,ik,ier)+plii*yinv(1,ik,ier)
         if (psmall(k,ier).ne.zero) then
            xqr(k,i,ier)=xqr(k,i,ier)*psmall(k,ier)
            xqi(k,i,ier)=xqi(k,i,ier)*psmall(k,ier)
         endif
      enddo
   enddo

   do i=1,mchan
      do j=1,mchan
         qr(i,j,ier)=0
         qi(i,j,ier)=0
      enddo
   enddo

   ij=0
   do i=1,nchan
      do j=1,i
         ij=ij+1
         kl=0
         do k=1,nchan
            do l=1,k
               kl=kl+1
               qr(kl,ij,ier)=xqr(i,k,ier)*xqr(j,l,ier)-xqi(i,k,ier)*xqi(j,l,ier)
               qi(kl,ij,ier)=xqr(i,k,ier)*xqi(j,l,ier)+xqi(i,k,ier)*xqr(j,l,ier)
               if (i.ne.j) then
                  qr(kl,ij,ier)=qr(kl,ij,ier)+xqr(j,k,ier)*xqr(i,l,ier)&
                    -xqi(j,k,ier)*xqi(i,l,ier)
                  qi(kl,ij,ier)=qi(kl,ij,ier)+xqr(j,k,ier)*xqi(i,l,ier)&
                    +xqi(j,k,ier)*xqr(i,l,ier)
               endif
            enddo
         enddo
      enddo
   enddo

   return
   end subroutine setqri

   subroutine settri(n,nent,next,nchan,nn,ier)
   !-------------------------------------------------------------------
   ! Generate tr & ti, which are 0.5 * [the real and
   ! imaginary parts of the partial of crss with respect to r]
   !-------------------------------------------------------------------
   ! externals
   integer::n,nent,next,nchan,nn,ier
   ! internals
   integer::kl,k,ij,i,m,ifs,ix,iy,l,j
   real(kr)::zz,dr,di,br,bi,ar,ai,cr,ci
   real(kr),parameter::zero=0
   real(kr),parameter::half=0.5e0_kr

   !--first do angle-integrated stuff
   do i=1,npp
      do j=1,ntriag
         tr(i,j,ier)=0
         ti(i,j,ier)=0
      enddo
   enddo

   !--Generate tr and ti, where
   !--tr(m,ij) = real part of partial (mth cross section) with
   !--respect to r(ij), except for c=4pi/E

   !--First, integrated elastic, diagonal in channel #'s
   kl=0
   do k=1,nent
      zz=zke(k,n,ier)**2
      kl=kl+k
      ij=0
      do i=1,nchan
         do j=1,i
            ij=ij+1
            if (qi(kl,ij,ier).ne.zero.or.qr(kl,ij,ier).ne.zero) then
               ar=qr(kl,ij,ier)*(-sin2ph(k,ier)*half)&
                 +qi(kl,ij,ier)*(-sinsqr(k,ier))
               ai=qi(kl,ij,ier)*(-sin2ph(k,ier)*half)&
                 -qr(kl,ij,ier)*(-sinsqr(k,ier))
               tr(1,ij,ier)=tr(1,ij,ier)+ar/zz
               ti(1,ij,ier)=ti(1,ij,ier)+ai/zz
            endif
         enddo
      enddo
   enddo

   !--Next, absorption only, diagonal in channel numbers
   kl=0
   do k=1,nent
      zz=2*zke(k,n,ier)**2
      kl=kl+k
      ij=0
      do i=1,nchan
         do j=1,i
            ij=ij+1
            if (qi(kl,ij,ier).ne.zero) tr(2,ij,ier)=&
              tr(2,ij,ier)+qi(kl,ij,ier)/zz
            if (qr(kl,ij,ier).ne.zero) ti(2,ij,ier)=&
              ti(2,ij,ier)-qr(kl,ij,ier)/zz
         enddo
      enddo
   enddo

   !--Next, not-necessarily diagonal pieces of elastic & capture
   kl=0
   do k=1,nent
      zz=zke(k,n,ier)**2
      do l=1,k
         kl=kl+1
         ij=0
         do i=1,nchan
            do j=1,i
               ij=ij+1
               if (qi(kl,ij,ier).ne.zero.or.qr(kl,ij,ier).ne.zero) then
                 ar=qr(kl,ij,ier)*xxxxr(kl,ier)+qi(kl,ij,ier)*xxxxi(kl,ier)
                 ai=qi(kl,ij,ier)*xxxxr(kl,ier)-qr(kl,ij,ier)*xxxxi(kl,ier)
                 if (k.ne.l) then
                    ar=ar*2
                    ai=ai*2
                 endif
                 tr(1,ij,ier)=tr(1,ij,ier)+ar/zz
                 ti(1,ij,ier)=ti(1,ij,ier)+ai/zz
                 tr(2,ij,ier)=tr(2,ij,ier)-ar/zz
                 ti(2,ij,ier)=ti(2,ij,ier)-ai/zz
               endif
            enddo
         enddo
      enddo
   enddo

   !--Next, reactions
   if (nchan.gt.nent) then
      kl=0
      do k=1,nchan
         m=ipp(k,n,ier)
         do l=1,k
            kl=kl+1
            if (l.le.nent.and.k.gt.nent) then
               zz=zke(l,n,ier)**2
               ij=0
               do i=1,nchan
                  do j=1,i
                     ij=ij+1
                     if (qi(kl,ij,ier).ne.zero.or.qr(kl,ij,ier).ne.zero) then
                        tr(m,ij,ier)=tr(m,ij,ier)&
                          +(qr(kl,ij,ier)*xxxxr(kl,ier)+&
                          qi(kl,ij,ier)*xxxxi(kl,ier))/zz
                        ti(m,ij,ier)=ti(m,ij,ier)+(qi(kl,ij,ier)*xxxxr(kl,ier)-&
                          qr(kl,ij,ier)*xxxxi(kl,ier))/zz
                     endif
                  enddo
               enddo
            endif
         enddo
      enddo
   endif

   !--end of angle-integrated stuff

   if (Want_Angular_Dist) then

      !--For angular distributions
      !--prtl (1-U)(kl) wrt R(ij) = prtl(1-U)(kl) wrt X(kl) * prtl X wrt R
      !--tx(1,ij,kl) = prtl Re(1-U) wrt ReR = -prtl Im (1-U) wrt ImR
      !--tx(2,ij,kl) = prtl Re(1-U) wrt ImR =  prtl Re (1-U) wrt ReR

      !--zero the tx array
      do i=1,2
         do j=1,ntriag
            do k=1,ntriag
               tx(i,j,k,ier)=0
            enddo
         enddo
      enddo

      !--first, do diagonal (in K & L) ... only for elastic
      kl=0
      do k=1,nent
         kl=kl+k
         ij=0
         do I=1,nchan
            do j=1,I
               ij=ij+1
               if (qi(kl,ij,ier).ne.zero.or.qr(kl,ij,ier).ne.zero) then
                  ar=qr(kl,ij,ier)*(-sin2ph(k,ier))&
                    -qi(kl,ij,ier)*(2*sinsqr(k,ier)-1)
                  ai=qi(kl,ij,ier)*(-sin2ph(k,ier))&
                    +qr(kl,ij,ier)*(2*sinsqr(k,ier)-1)
                  tx(1,ij,kl,ier)=ar
                  tx(2,ij,kl,ier)=ai
               endif
            enddo
         enddo
      enddo

      !--now, do off-diagonal (in KL)
! ******************************** note that this is not right!
      ifs=0
! ******************************** note that this is not right!
      kl=0
      do l=1,nchan
         do k=1,l
            kl=kl+1
            if (ifs.eq.0.and.k.le.nent) then
               if (k.ne.l) then
                  ij=0
                  do I=1,nchan
                     ij=(i*(i-1))/2
                     do j=1,i
                        ij=ij+1
                        if (qi(kl,ij,ier).ne.zero.or.qr(kl,ij,ier).ne.zero) then
                           ar=-qr(kl,ij,ier)*cscs(2,kl,ier)&
                             +qi(kl,ij,ier)*cscs(1,kl,ier)
                           ai=-qi(kl,ij,ier)*cscs(2,kl,ier)&
                             -qr(kl,ij,ier)*cscs(1,kl,ier)
                           tx(1,ij,kl,ier)=ar
                           tx(2,ij,kl,ier)=ai
                        endif
                     enddo
                  enddo
               endif
            endif
         enddo
      enddo

      !--Now multiply by Coulomb phase shift if needed

      do k=1,nchan
         if (ifs.eq.0) then
            if (zeta(k,n,ier).ne.zero.and.su.gt.echan(k,n,ier)) then
               call gcphase(cr,ci,lspin(k,n,ier),echan(k,n,ier),&
                 zeta(k,n,ier),su)
               ix=1
            else
               cr=1
               ci=0
               ix=0
            endif

            kl=(k*(k-1))/2
            do l=1,k
               kl=kl+1
               if (k.eq.l) then
                  dr=cr
                  di=ci
                  iy=ix
               else
                  if (zeta(l,n,ier).ne.zero.and.su.gt.echan(l,n,ier)) then
                     call gcphase(dr,di,lspin(l,n,ier),&
                       echan(l,n,ier),zeta(l,n,ier),su)
                     iy=1
                  else
                     dr=1
                     di=0
                     iy=0
                  endif
               endif

               if (ix.eq.1.or.iy.eq.1) then
                  ij=0
                  do i=1,nchan
                     do j=1,i
                        ij=ij+1
                        ar=tx(1,ij,kl,ier)
                        ai=tx(2,ij,kl,ier)
                        if (ix.eq.0) then
                           br=ar
                           bi=ai
                        else
                           br=ar*cr-ai*ci
                           bi=ar*ci+ai*cr
                        endif
                        if (iy.eq.0) then
                           ar=Br
                           ai=bi
                        else
                           ar=br*dr-bi*di
                           ai=br*di+bi*dr
                        endif
                        tx(1,ij,kl,ier)=ar
                        tx(2,ij,kl,ier)=ai
                     enddo
                  enddo
               endif

            enddo
         endif
      enddo
   endif

   return
   end subroutine settri

   subroutine derres(n,nent,nchan,minres,maxres,kstart,npr,nn,ier)
   !-------------------------------------------------------------------
   ! Generate deriv(k,ipar) = partial crss(k) with respect
   ! to u(ipar) for resonance parameters.
   !-------------------------------------------------------------------
   ! externals
   integer::n,nent,nchan,minres,maxres,kstart,npr,nn,ier
   ! internals
   integer::mm,m,ij,ip,i,j
   real(kr),parameter::zero=0

   do mm=1,npr
      m=kstart+mm
      ! note that m = parameter number
      ij=0
      do ip=1,npp
         ddddd(ip,ier)=0
      enddo
      do i=1,nchan
         do j=1,i
            ij=ij+1
            if (pii(ij,m,ier).ne.zero) then
               do ip=1,npp
                  ddddd(ip,ier)=ddddd(ip,ier)-pii(ij,m,ier)*ti(ip,ij,ier)
               enddo
            endif
            if (pr(ij,m,ier).ne.zero) then
               do ip=1,npp
                  ddddd(ip,ier)=ddddd(ip,ier)+pr(ij,m,ier)*tr(ip,ij,ier)
               enddo
            endif
         enddo
      enddo
      do ip=1,npp
         deriv(ip,m,ier)=goj(n,ier)*ddddd(ip,ier)+deriv(ip,m,ier)
      enddo
   enddo

   return
   end subroutine derres

   subroutine derext(n,nchan,jstart,nnnn,ier)
   !-------------------------------------------------------------------
   ! externals
   integer::n,nchan,jstart,nnnn,ier
   ! internals
   integer::ij,i,m
   real(kr)::a

   !--angle-integrated
   ij=0
   do i=1,nchan
      ij=ij+I
      !  note that tr = 1/2 times partial (sigma) wrt (re R),
      !  ergo need to multiply by 2 here
      a=2*goj(n,ier)
      jstart=jstart+1
      do m=1,npp
         if (m.le.2) then
            deriv(m,jstart,ier)=-tr(m,ij,ier)*a*(parext(5,i,nnnn,ier)+&
              parext(6,i,nnnn,ier)*parext(1,i,nnnn,ier))/&
              (su-parext(1,i,nnnn,ier))
         endif
      enddo
      jstart=jstart+1
      do m=1,npp
         if (m.le.2) then
            deriv(m,jstart,ier)=-tr(m,ij,ier)*a*(parext(5,i,nnnn,ier)+&
              parext(6,I,nnnn,ier)*parext(2,i,nnnn,ier))/&
             (parext(2,i,nnnn,ier)-su)
         endif
      enddo
      jstart=jstart+1
      do m=1,npp
         if (m.le.2) then
            deriv(m,jstart,ier)=tr(m,ij,ier)*a
         endif
      enddo
      jstart=jstart+1
      do m=1,npp
         if (m.le.2) then
            deriv(m,jstart,ier)=tr(m,ij,ier)*a*su
         endif
      enddo
      jstart=jstart+1
      do m=1,npp
         if (m.le.2) then
            deriv(m,jstart,ier)=-2*tr(m,ij,ier)*a*sqrt(parext(5,i,nnnn,ier))*&
              log((parext(2,i,nnnn,ier)-su)/(su-parext(1,i,nnnn,ier)))
              !   Remember that the u-parameter is the
              !   square root of parext(5)
         endif
      enddo
      jstart=jstart+1
      do m=1,npp
         if (m.le.2) then
            deriv(m,jstart,ier)=-tr(m,ij,ier)*a*&
              ((parext(2,i,nnnn,ier)-parext(1,i,nnnn,ier))+&
              su*log((parext(2,i,nnnn,ier)-su)/(su-parext(1,i,nnnn,ier))))
         endif
      enddo
      jstart=jstart+1
      do m=1,npp
         if (m.le.2) then
            deriv(m,jstart,ier)=tr(m,ij,ier)*a*su**2
         endif
      enddo
   enddo

   return
   end subroutine derext

   subroutine setleg(igroup,kount,ier)
   !-------------------------------------------------------------------
   ! Set somega(l,ipp) = coefficient of Legendre polynomial
   ! P-sub-(L-1) for ipp_th cross section.
   !-------------------------------------------------------------------
   use util ! provides error
   ! externals
   integer::igroup,kount,ier
   ! internal
   integer::ngr,kn,nchn,nchanx,mgr,km,mchan,mchanx
   integer::l,ip,jx,iax,ibx,nppx
   real(kr),dimension(:)::ar(50),ai(50),br(50)
   real(kr),parameter::zero=0

   nppx=npp
   if (npp.eq.2) nppx=1
   if (nppx.gt.50) call error('setleg','nppx too large',' ')
   ngr=igroup
   kn=0
   do nchn=1,nchan(ngr,ier)
      do nchanx=1,nchn
         kn= kn+1
         if (nchanx.le.nent(ngr,ier)) then
            iax=0
            do ip=1,nppx
               ar(ip)=crssx(1,nchanx,nchn,ip,ngr,ier)
               ai(ip)=crssx(2,nchanx,nchn,ip,ngr,ier)
               if (ar(ip).ne.zero) iax=1
               if (ai(ip).ne.zero) iax=1
            enddo
            if (iax.ne.0) then
               do mgr=1,ngr
                  km=0
                  do mchan=1,nchan(mgr,ier)
                     do mchanx=1,mchan
                        km=km+1
                        if (mchanx.le.nent(mgr,ier)) then
                           ibx=0
                           do ip=1,nppx
                              br(ip)=ar(ip)*crssx(1,mchanx,mchan,ip,mgr,ier)&
                                +ai(ip)*crssx(2,mchanx,mchan,ip,mgr,ier)
                              if (br(ip).ne.zero) ibx=1
                           enddo
                           if (ibx.ne.0) then
                              do l=1,lllmax
                                 jx=jxx(kxlmn,l,km,mgr,kn,ngr,kount)
                                 if (jx.eq.kxlmn(kount)) then
                                    if (xlmn(kount).ne.zero) then
                                       do ip=1,nppx
                                          Coef_Leg(l,ip,ier)=&
                                            Coef_Leg(l,ip,ier)&
                                            +xlmn(kount)*Br(ip)
                                       enddo
                                    endif
                                    kount=kount+1
                                 endif
                              enddo
                           endif
                        endif
                     enddo
                  enddo
               enddo
            endif
         endif
      enddo
   enddo

   return
   end subroutine setleg

   integer function jxx(kxlmn,l,km,mgr,kn,ngr,kount)
   !-------------------------------------------------------------------
   ! Private function for setleg.
   !-------------------------------------------------------------------
   ! externals
   integer::kxlmn(*)
   integer::l,km,mgr,kn,ngr,kount
   ! internals
   integer::jx,ik

   jx=((((ngr-1)*ntriag+kn-1)*ngroup+mgr-1)*ntriag+km-1)*lllmax+l
   if (jx.gt.kxlmn(kount)) then
      if (kount.eq.kkxlmn) go to 20
         do ik=kount+1,kkxlmn
            if (jx.le.kxlmn(ik)) go to 10
         enddo
         kount=kkxlmn
         go to 20
        10 continue
         kount=ik
     20 continue
   endif
   jxx=jx
   return
   end function jxx

   subroutine wrongi(name,i)
   use mainio ! provides nsyso,nsyse
   character*(*) name
   integer::i
      write(nsyse,10) name,i
      write(nsyso,10) name,i
     10 format(/' $$$ Variable "',a10,'" probably has wrong value',&
         i15,' $$$')
   return
   end subroutine wrongi

   subroutine desammy
   !-------------------------------------------------------------------
   ! deallocate various arrays to clean up after the sammy calculation.
   !-------------------------------------------------------------------

   ! channel array
   deallocate(nchan)

   ! angle array
   if (Want_Angular_Dist) then
      deallocate(xlmn)
      deallocate(kxlmn)
   endif

   ! allo2 arrays
   if (Want_Partial_Derivs) then
      deallocate(br)
      deallocate(bi)
      deallocate(pr)
      deallocate(pii)
   endif

   deallocate(crss)
   deallocate(sigmas)

   if (Want_Partial_Derivs) then
      deallocate(deriv)
      deallocate(dsigma)
   endif

   if (Want_Angular_Dist) then
      deallocate(crssx)
      deallocate(Coef_Leg)
      if (Want_Partial_Derivs) then
         deallocate(derivx)
         deallocate(D_Coef_Leg)
      endif
   endif

   deallocate(alphar)
   deallocate(alphai)
   deallocate(difen)
   deallocate(xden)

   if (Want_Partial_Derivs) then
      deallocate(upr)
      deallocate(upi)
   endif

   deallocate(sinsqr)
   deallocate(sin2ph)
   deallocate(dphi)
   deallocate(dpdr)
   deallocate(dsdr)
   deallocate(cscs)
   deallocate(cosphi)
   deallocate(sinphi)

   deallocate(rootp)
   deallocate(elinvr)
   deallocate(elinvi)
   deallocate(psmall)
   deallocate(xxxxr)
   deallocate(xxxxi)

   deallocate(xqr)
   deallocate(xqi)

   deallocate(qr)
   deallocate(qi)

   if (Want_Partial_Derivs) then
      deallocate(tr)
      deallocate(ti)
      deallocate(tx)
   endif

   deallocate(ddddd)

   deallocate(yinv)

   deallocate(rmat)
   deallocate(ymat)

   deallocate(par)

 ! allo1 arrays

   deallocate(ema)
   deallocate(emb)
   deallocate(kza)
   deallocate(kzb)
   deallocate(spina)
   deallocate(spinb)
   deallocate(qqq)
   deallocate(ishift)
   deallocate(lpent)
   deallocate(mt)
   deallocate(pa)
   deallocate(pb)
   deallocate(sspin)
   deallocate(parity)
   deallocate(nresg)
   deallocate(ipp)
   deallocate(lspin)
   deallocate(chspin)
   deallocate(bound)
   deallocate(rdeff)
   deallocate(rdtru)
   deallocate(eres)
   deallocate(gamgam)
   deallocate(gamma)
   deallocate(parext)
   deallocate(nent)
   deallocate(next)
   deallocate(goj)
   deallocate(zke)
   deallocate(zkfe)
   deallocate(zkte)
   deallocate(zeta)
   deallocate(echan)
   deallocate(betapr)
   deallocate(gbetpr)
   deallocate(beta)
   deallocate(uuuu)
   deallocate(duuu)
   deallocate(iduu)

   return
   end subroutine desammy

end module samm


module mixm
   ! provides subroutine mixr for NJOY2016
   use locale
   implicit none
   private
   public mixr

   ! global variables
   integer::ntape
   integer,parameter::nninmx=10
   integer,parameter::nmtmx=20
   integer,parameter::nmatmx=nninmx
   integer,dimension(nninmx)::jtape,nrt,npt,irt,ipt,ip1t,ip2t,nbt,nwt

contains

   subroutine mixr
   !-------------------------------------------------------------------
   !
   !  mixr
   !
   !  Construct a new pendf tape with a specified set of
   !  reactions that are specified linear combinations of the
   !  cross sections from the input tapes.  Mixr can also be
   !  used for endf tapes, but the input interpolation laws
   !  are ignored.  This module can be used to construct mixed
   !  reactions for plotting (for example, elemental cross
   !  sections).  The output file contains files 1 and 3 only.
   !  Linear-linear interpolation is assumed.
   !
   !  user input --
   !
   !  card 1  --  units
   !     nout      output unit for mixed cross sections
   !     nin1      first input unit (endf or pendf)
   !     nin2      second input unit
   !      ...      continue for nnin<=nninmx (=10) input units
   !
   !  card 2  --  reaction list
   !      mtn      list of nmt<=nmtmx (=20) mt numbers for
   !               the output reactions
   !
   !  card 3  --  material list
   !     matn,     list of nmat<=nmatmx (=nninmx=10) pairs (matn,wtn)
   !      wtn      of material numbers and associated weighting factors
   !
   !  card 4  --  temperature
   !     temp      temperature (use zero except for pendf tapes)
   !
   !  card 5  --  output material
   !     matd      material number
   !     za        za value
   !     awr       awr value
   !
   !  card 6  --  file 1 comment card
   !     des       description (66 char max)
   !
   !-------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use util ! provides timer,openz,error,mess,closz,sigfig
   use endf ! provides endf routines and variables
   ! internals
   integer,parameter::mxlna=300000
   integer::nscr,nbuf,nc,i,nout,nnin,nmt,nmat,matd
   integer::ndes,ibuf,iscr,nb,nw,idone,nx,inow,imt
   integer::mtd,je,ne,jend,idis
   integer,parameter::lrp=-1
   integer::nsub
   real(kr)::temp,za,awr,t,test,enext,reset,enxt,e,sig
   real(kr)::time,y
   real(kr)::awi,emax
   character(66)::strng
   character(4)::str(17)
   integer,dimension(nninmx)::nin
   integer,dimension(nmtmx)::mtn
   integer,dimension(nmatmx)::matn,jscr
   real(kr),dimension(nmatmx)::wtn
   real(kr),dimension(mxlna)::a
   real(kr)::c(2)
   real(kr)::des(17)
   equivalence(str(1),des(1))
   real(kr),parameter::big=1.e10_kr

   !--initialize
   call timer(time)
   write(nsyso,&
     '(/'' mixr...calculation of mixed cross sections'',&
     &26x,f8.1,''s'')') time
   write(nsyse,'(/'' mixr...'',61x,f8.1,''s'')') time
   nscr=10
   nbuf=5000
   nc=2
   call openz(-nscr,0)

   !--read user input
   do i=1,nninmx
      nin(i)=0
   enddo
   read(nsysi,*) nout,(nin(i),i=1,nninmx)
   do i=1,nninmx
      if (nin(i).ne.0) nnin=i
   enddo
   do i=1,nmtmx
      mtn(i)=0
   enddo
   read(nsysi,*) (mtn(i),i=1,nmtmx)
   do i=1,nmtmx
      if (mtn(i).ne.0) nmt=i
   enddo
   do i=1,nmatmx
      matn(i)=0
   enddo
   read(nsysi,*) (matn(i),wtn(i),i=1,nmatmx)
   do i=1,nmatmx
      if (matn(i).ne.0) nmat=i
   enddo
   read(nsysi,*) temp
   read(nsysi,*) matd,za,awr
   strng=' '
   read(nsysi,*) strng
   read(strng,'(16a4,a2)') (str(i),i=1,17)
   ndes=17
   write(nsyso,&
     '(/'' unit for output tape ................. '',i10)') nout
   do i=1,nnin
      write(nsyso,&
        '( '' unit for input tape .................. '',i10)')&
        nin(i)
   enddo
   write(nsyso,&
     &'( '' mt numbers for output tape ........... '',/10x,20i6)')&
     (mtn(i),i=1,nmt)
   write(nsyso,'( '' mat,weight values for mixes .......... ''/&
     &(10x,i6,f10.5))') (matn(i),wtn(i),i=1,nmat)

   !--allocate storage
   ibuf=1
   iscr=ibuf+nbuf
   jscr(1)=iscr+npage+50
   do i=2,nmat
      jscr(i)=jscr(i-1)+npage+50
   enddo

   !--open the output and input units
   !--and locate the desired materials and temperature.
   call openz(nout,1)
   awi=1
   emax=20.e6_kr
   nsub=10
   do i=1,nmat
      call openz(nin(i),0)
      call tpidio(nin(i),0,0,a(iscr),nb,nw)
      idone=0
      do while (idone.eq.0)
         call contio(nin(i),0,0,a(iscr),nb,nw)
         if (math.lt.0) call error('mixr',&
           'mat and temp not found',' ')
         if (math.eq.matn(i)) then
            call contio(nin(i),0,0,a(iscr),nb,nw)
            if (n1h.ne.0) then
               iverf=4
            else if (n2h.eq.0) then
               iverf=5
            else
               iverf=6
            endif
            call skiprz(nin(i),-1)
            if (iverf.gt.4) call contio(nin(i),0,0,a(iscr),nb,nw)
            if (iverf.gt.5) then
               call contio(nin(i),0,0,a(iscr),nb,nw)
               if (a(iscr+1).gt.emax) emax=a(iscr+1)
               if (i.eq.1) then
                  awi=a(iscr)
                  if (awi.lt.0.1) nsub=0
                  if (awi.gt.0.9980.and.awi.lt.0.9990) nsub=10010
                  if (awi.gt.1.9950.and.awi.lt.1.9970) nsub=10020
                  if (awi.gt.2.9895.and.awi.lt.2.9897) nsub=10030
                  if (awi.gt.2.9890.and.awi.lt.2.9891) nsub=20030
                  if (awi.gt.3.9670.and.awi.lt.3.9680) nsub=20040
               endif
            endif
            call hdatio(nin(i),0,0,a(iscr),nb,nw)
            t=a(iscr)
            test=1
            if (abs(t-temp).lt.test) then
               idone=1
            else
               call tomend(nin(i),0,0,a(iscr))
            endif
         else
            call tomend(nin(i),0,0,a(iscr))
         endif
      enddo
      call findf(matn(i),3,0,nin(i))
   enddo

   !--write the tape id record and file 1 for the output tape
   math=1
   call afend(nout,0)
   math=matd
   mfh=1
   mth=451
   nx=nmt
   a(iscr)=za
   a(iscr+1)=awr
   a(iscr+2)=lrp
   a(iscr+3)=0
   a(iscr+4)=0
   a(iscr+5)=0
   if (iverf.le.4) a(iscr+5)=nx
   call contio(0,nout,0,a(iscr),nb,nw)
   if (iverf.gt.4) then
      a(iscr)=0
      a(iscr+1)=0
      a(iscr+2)=0
      a(iscr+3)=0
      a(iscr+4)=0
      a(iscr+5)=0
      if (iverf.eq.6) a(iscr+5)=6
      call contio(0,nout,0,a(iscr),nb,nw)
      if (iverf.ge.6) then
         a(iscr)=awi
         a(iscr+1)=emax
         a(iscr+2)=0
         a(iscr+3)=0
         a(iscr+4)=nsub
         a(iscr+5)=0
         call contio(0,nout,0,a(iscr),nb,nw)
      endif
   endif
   a(iscr)=temp
   a(iscr+1)=0
   a(iscr+2)=0
   if (iverf.ge.6) a(iscr+2)=1
   a(iscr+3)=0
   a(iscr+4)=ndes
   a(iscr+5)=0
   if (iverf.ge.5) a(iscr+5)=nx
   do i=1,ndes
      a(iscr+5+i)=des(i)
   enddo
   call hdatio(0,nout,0,a(iscr),nb,nw)
   inow=iscr
   do i=1,nmt
      a(inow)=0
      a(inow+1)=0
      a(inow+2)=3
      a(inow+3)=mtn(i)
      a(inow+4)=0
      a(inow+5)=0
      inow=inow+6
   enddo
   nw=nmt
   call dictio(0,nout,0,a(iscr),nb,nw)
   call afend(nout,0)

   !--loop over the desired reactions.
   do imt=1,nmt
      mtd=mtn(imt)

      !--find the section for mt on each input tape.
      enext=big
      reset=-1
      call gety(reset,enxt,idis,y,0,a(iscr))
      do i=1,nmat
         idone=0
         do while (idone.eq.0)
            call contio(nin(i),0,0,a(iscr),nb,nw)
            if (mfh.eq.0.or.mth.eq.mtd) then
               idone=1
            else
               call tosend(nin(i),0,0,a(iscr))
            endif
         enddo
         if (mth.eq.mtd) then
            e=0
            call gety(e,enxt,idis,y,nin(i),a(jscr(i)))
            if (enxt.lt.enext) enext=enxt
            write(nsyso,&
              '('' found section for mt='',i4,'' and mat='',i4)')&
              mtd,matn(i)
         else
            write(strng,&
              '(''mt='',i3,'' not present for mat='',i4)')&
              mtd,matn(i)
            call mess('mixr',strng,' ')
            call findf(matn(i),3,0,nin(i))
         endif
      enddo

      !--loop over the energies to get a union grid
      je=0
      do while (enext.lt.big)
         e=enext

         !--mix the cross sections at this energy
         enext=big
         sig=0
         do i=1,nmat
            call gety(e,enxt,idis,y,nin(i),a(jscr(i)))
            if (enxt.lt.enext) enext=enxt
            sig=sig+wtn(i)*y
         enddo

         !--save the mixed value on the scratch tape
         je=je+1
         if (enext.ge.big) je=-je
         c(1)=e
         c(2)=sig
         call loada(je,c,nc,nscr,a(ibuf),nbuf)

      !--continue the energy loop
      enddo
      ne=-je
      do i=1,nmat
         call tosend(nin(i),0,0,a(iscr))
      enddo

      !--write head record for this mt
      math=matd
      mfh=3
      mth=mtd
      a(iscr)=za
      a(iscr+1)=awr
      a(iscr+2)=0
      a(iscr+3)=0
      a(iscr+4)=0
      a(iscr+5)=0

      !--write tab1 record for this mt
      call contio(0,nout,0,a(iscr),nb,nw)
      a(iscr)=0
      a(iscr+1)=0
      a(iscr+2)=0
      a(iscr+3)=0
      a(iscr+4)=1
      a(iscr+5)=ne
      a(iscr+6)=ne
      a(iscr+7)=2
      inow=iscr+8
      nb=0
      jend=npage/2
      if (jend.gt.ne) jend=ne
      je=0
      idone=0
      do while (je.lt.ne.and.idone.eq.0)
         je=je+1
         call finda(je,c,nc,nscr,a(ibuf),nbuf)
         a(inow)=c(1)
         a(inow+1)=sigfig(c(2),7,0)
         inow=inow+2
         if (inow.gt.mxlna) then
            call error('mixr','mxlna array limit exceeded',' ')
         endif
         if (je.ge.jend) then
            if (nb.eq.0) then
               call tab1io(0,nout,0,a(iscr),nb,nw)
            else
               call moreio(0,nout,0,a(iscr),nb,nw)
            endif
            if (nb.eq.0) then
               idone=1
            else
               inow=iscr
               jend=jend+npage/2
               if (jend.gt.ne) jend=ne
            endif
         endif
      enddo
      call asend(nout,0)

   !--continue the reaction loop
   enddo

   !--finish off the output tape
   call afend(nout,0)
   call amend(nout,0)
   call atend(nout,0)

   !--close the units being used and return
   call closz(nout)
   do i=1,nmat
      call closz(nin(i))
   enddo
   call closz(-nscr)
   call timer(time)
   write(nsyso,'(/69x,f8.1,''s''/&
     &7(''**********''),''*******'')') time
   return
   end subroutine mixr

   subroutine gety(x,xnext,idis,y,itape,a)
   !-------------------------------------------------------------------
   ! Retrieve y(x) from an ENDF TAB1 structure using paged BCD or
   ! blocked binary formats.  Call with x=0 to read in first page
   ! or block of data and initialize pointers.  Routine assumes
   ! values will be called in ascending order.  Here, xnext is the
   ! first data grid point greater than x unless x is the last point.
   ! This version will keep track of pointers for up to nninmx units.
   ! Call with x=-1 to clear the pointers before each group of files.
   !-------------------------------------------------------------------
   use util ! provides error
   use endf ! provides tab1io,moreio,terp1
   ! externals
   integer::idis,itape
   real(kr)::x,xnext,y,a(*)
   ! internals
   integer::nwtot,nr,np,ip1,ip2,ir,ip,i,ktape,ln
   integer::nbx,int,nb,nw,lt
   real(kr)::xn
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::zero=0
   save lt

   !--branch on value of x
   idis=0
   if (x.eq.zero) go to 100
   if (x.gt.zero) go to 115

   !--clear pointer storage
   ntape=0
   return

   !--read first page or block of data and initialize
  100 continue
   ntape=ntape+1
   jtape(ntape)=itape
   call tab1io(itape,0,0,a,nb,nw)
   nwtot=nw
   nr=nint(a(5))
   np=nint(a(6))
   lt=6+2*nr
   ip1=1
   ip2=(nw-lt)/2
   if (nb.eq.0) ip2=ip2+2
   ir=1
   ip=2
   xnext=a(lt+1)

   !--save pointers and return
   nrt(ntape)=nr
   npt(ntape)=np
   irt(ntape)=ir
   ipt(ntape)=ip
   ip1t(ntape)=ip1
   ip2t(ntape)=ip2
   nbt(ntape)=nb
   nwt(ntape)=nwtot
   return

   !--restore pointers
  115 continue
   if (ntape.eq.0)&
     call error('gety','not properly initialized',' ')
   do 120 i=1,ntape
   if (jtape(i).ne.itape) go to 120
   ktape=i
   nr=nrt(i)
   np=npt(i)
   ir=irt(i)
   ip=ipt(i)
   ip1=ip1t(i)
   ip2=ip2t(i)
   nb=nbt(i)
   nwtot=nwt(i)
   go to 125
  120 continue
   y=0
   xnext=big
   return

   !--is x in this panel
  125 continue
   ln=2*(ip-ip1)+lt
   if (x.lt.a(ln-1)) go to 135
   if (x.lt.a(ln+1)) go to 130
   if (ip.eq.np) go to 140

   !--no.  move up to next range.
   !--read in new page of data if needed.
   ip=ip+1
   nbx=nint(a(5+2*ir))
   if (ip.gt.nbx) ir=ir+1
   if (ip.lt.ip2) go to 125
   if (nb.eq.0) go to 130
   a(lt+1)=a(nwtot-3)
   a(lt+2)=a(nwtot-2)
   a(lt+3)=a(nwtot-1)
   a(lt+4)=a(nwtot)
   call moreio(itape,0,0,a(lt+5),nb,nw)
   nwtot=nw+lt+4
   ip1=ip-1
   ip2=ip1+nw/2+1
   if (nb.eq.0) ip2=ip2+2
   go to 125

   !--yes.  interpolate for desired value
  130 continue
   int=nint(a(6+2*ir))
   if (int.eq.1) idis=1
   call terp1(a(ln-1),a(ln),a(ln+1),a(ln+2),x,y,int)
   xnext=a(ln+1)
   if ((ln+3).gt.nwtot.and.nb.eq.0) return
   xn=a(ln+3)
   if (xn.eq.xnext) idis=1
   go to 150

   !--special branch for x outside range of table
  135 continue
   y=0
   xnext=a(ln-1)
   go to 150

   !--special branch for last point
  140 continue
   y=a(ln+2)
   xnext=big

   !--save pointers and return
  150 continue
   nrt(ktape)=nr
   npt(ktape)=np
   irt(ktape)=ir
   ipt(ktape)=ip
   ip1t(ktape)=ip1
   ip2t(ktape)=ip2
   nbt(ktape)=nb
   nwt(ktape)=nwtot
   return
   end subroutine gety

end module mixm


module util
   ! Utility routines for NJOY2016.
   use locale
   implicit none
   private

   !--Public routines
   public error,mess
   public timer,dater,wclock
   public repoz,skiprz,openz,closz
   public loada,finda,scana
   public sigfig
   public a10

contains

   subroutine error(from,mess1,mess2)
   !--------------------------------------------------------------------
   ! Machine dependent error exit routine.
   ! Use system fatal error exit if available.
   !--------------------------------------------------------------------
   use mainio ! provides nsyso,nsyse
   ! externals
   character::from*(*),mess1*(*),mess2*(*)

   write(nsyso,'(/'' ***error in '',a,''***'',a)') from,trim(mess1)
   if (len_trim(mess2).gt.0) write(nsyso,'(22x,a)') trim(mess2)
   write(nsyso,'('' '')')
   if (nsyse.ne.nsyso) then
      write(nsyse,'(/'' ***error in '',a,''***'',a)') from,trim(mess1)
      if (len_trim(mess2).gt.0) write(nsyse,'(22x,a)') trim(mess2)
      write(nsyse,'('' '')')
   endif
   stop 77
   end subroutine error

   subroutine mess(from,mess1,mess2)
   !--------------------------------------------------------------------
   ! Message routine (not fatal).
   !--------------------------------------------------------------------
   use mainio ! provides nsyso,nsyse
   ! externals
   character::from*(*),mess1*(*),mess2*(*)

   write(nsyso,'(/'' ---message from '',a,''---'',a)') from,trim(mess1)
   if (len_trim(mess2).gt.0) write(nsyso,'(26x,a)') trim(mess2)
   if (nsyse.ne.nsyso) then
      write(nsyse,'(/'' ---message from '',a,''---'',a)') &
        from,trim(mess1)
      if (len_trim(mess2).gt.0) write(nsyse,'(26x,a)') trim(mess2)
   endif
   return
   end subroutine mess

   subroutine timer(time)
   !--------------------------------------------------------------------
   ! Get elapsed time in seconds from beginning of execution.
   ! Should be "charge time" (e.g., cp+pp+sys).
   !--------------------------------------------------------------------
   ! externals
   real(kr)::time
   ! internals
   real(k4),dimension(2)::tarray
   !real(k8),external::etime

   !time=etime(tarray)
   !time=tarray(1)+tarray(2)
   call cpu_time(time)
   return
   end subroutine timer

   subroutine dater(hdate)
   !--------------------------------------------------------------------
   ! Return the date as as 8 character string.
   ! Either mm/dd/yy or ddmmmyy (eg, 03jun83) are acceptable.
   !--------------------------------------------------------------------
   ! externals
   character(8)::hdate
   ! internals
   character(8)::date
   character(10)::time
   character(5)::zone
   integer,dimension(8)::values
   integer i
   intrinsic date_and_time

   call date_and_time(date,time,zone,values)
   write(hdate,'(i2,''/'',i2,''/'',i2)') &
      values(2),values(3),mod(values(1),100)
   do i=1,8
      if (hdate(i:i).eq.' ') hdate(i:i)='0'
   end do
   return
   end subroutine dater

   subroutine wclock(htime)
   !--------------------------------------------------------------------
   ! Return the wall clock time as an 8 character string.
   ! For example, 12:13:47.
   !--------------------------------------------------------------------
   ! externals
   character(8)::htime
   ! internals
   character(8)::date
   character(10)::time
   character(5)::zone
   integer,dimension(8)::values
   integer i
   intrinsic date_and_time

   call date_and_time(date,time,zone,values)
   write(htime,'(i2,'':'',i2,'':'',i2)') &
      values(5),values(6),values(7)
   do i=1,8
      if (htime(i:i).eq.' ') htime(i:i)='0'
   end do
   return
   end subroutine wclock

   subroutine repoz(ntape)
   !--------------------------------------------------------------------
   ! Rewind an ENDF/B coded or blocked binary tape.
   ! Positive units are coded, negative ones are blocked binary.
   !--------------------------------------------------------------------
   ! externals
   integer::ntape
   ! internals
   integer::i

   if (ntape.eq.0) return
   i=iabs(ntape)
   rewind i
   return
   end subroutine repoz

   subroutine skiprz(nunit,nrt)
   !--------------------------------------------------------------------
   ! Skip forward or backward on a coded or blocked binary tape.
   ! Positive units are coded, negative ones are blocked binary.
   ! No action for zero units.
   !--------------------------------------------------------------------
   ! externals
   integer::nunit,nrt
   ! internals
   integer::nun,nr,i
   character(1)::idum

   if (nunit.eq.0) return
   nun=iabs(nunit)
   if (nrt.eq.0) return
   nun=iabs(nunit)
   if (nrt.eq.0) return
   nr=iabs(nrt)
   if (nrt.gt.0) then
      do i=1,nr
         if (nunit.gt.0) read(nun,'(a1)') idum
         if (nunit.lt.0) read(nun) idum
      enddo
   else if (nrt.lt.0) then
      do i=1,nr
         backspace nun
      enddo
   endif
   return
   end subroutine skiprz

   subroutine openz(lun,new)
   !--------------------------------------------------------------------
   ! System-dependent routine to open files.
   ! No action if iabs(lun) ge 5 and le 7.
   ! Mode--coded (formatted) if lun gt 0
   !       binary (unformatted) if lun lt 0
   ! Destroy--on close or job termination if iabs(lun) ge 10 and
   !                                      iabs(lun) lt 20
   ! Status--if new=1, destroy lun if it already exists, then
   !         open a new version.
   !--------------------------------------------------------------------
   ! externals
   integer::lun,new
   ! internals
   integer::nun
   character::for*15,age*7,fn*6
   logical::there

   nun=iabs(lun)
   if ((nun.ge.5.and.nun.le.7).or.nun.eq.0) return
   if (nun.gt.99) call error('openz','illegal unit number.',' ')
   ! construct file name
   write(fn,'(''tape'',i2)') nun
   ! set format based on sign of unit number
   for='formatted'
   if (lun.lt.0) for='unformatted'
   if (nun.ge.10.and.nun.le.19) then
      ! scratch units
      age='scratch'
      open(nun,form=for,status=age)
   else
      ! regular units
      if (new.ne.1) then
         ! existing units
         age='old'
      else
         ! new units
         age='new'
         inquire(file=fn,exist=there)
         if (there) then
            open(nun,file=fn,status='old')
            close(nun,status='delete')
         endif
      endif
      ! open the connection
      open(nun,file=fn,form=for,status=age)
   endif
   return
   end subroutine openz

   subroutine closz(lun)
   !--------------------------------------------------------------------
   ! File closing routine.
   ! No action if iabs(lun) lt 10.
   ! Destroy flag set by openz.
   !--------------------------------------------------------------------
   ! externals
   integer::lun
   ! internals
   integer::nun

   nun=iabs(lun)
   if (nun.lt.10) return
   if (nun.gt.99) call error('closz','illegal unit number.',' ')
   close(nun)
   return
   end subroutine closz

   subroutine loada(i,a,na,ntape,buf,nbuf)
   !--------------------------------------------------------------------
   ! Buffered sequential i/o routine.
   ! Store na elements of array a into
   ! core buffer and associated binary tape.
   !--------------------------------------------------------------------
   !externals
   integer::i,na,ntape,nbuf
   real(kr)::a(na),buf(nbuf)
   !internals
   integer::nl,j,ix,inow,k
   real(kr)::x,xnow

   nl=nbuf/na
   j=iabs(i)
   x=j
   x=x/nl
   ix=j/nl
   xnow=(x-ix)*nl
   inow=nint(xnow)
   if (inow.eq.0) inow=nl
   k=na*(inow-1)
   do j=1,na
      k=k+1
      buf(k)=a(j)
   enddo
   if (i.eq.1) call repoz(-ntape)
   if (inow.eq.nl.or.i.lt.0) write(ntape) buf
   return
   end subroutine loada

   subroutine finda(i,a,na,ntape,buf,nbuf)
   !--------------------------------------------------------------------
   ! Buffered sequential i/o routine.
   ! Find na elements of array a from
   ! core buffer and associated binary tape.
   !--------------------------------------------------------------------
   use mainio ! provides nsyso
!   use endf ! provides endf routines and variables
   use snl     ! provides SNL
   ! externals
   integer::i,na,ntape,nbuf
   real(kr)::a(na),buf(nbuf)
   ! internals
   integer::nl,inow,k,j

   nl=nbuf/na
   inow=mod(i,nl)
   if (inow.eq.0) inow=nl
   if (i.eq.1) call repoz(-ntape)
   if (inow.eq.1) read(ntape) buf
   k=na*(inow-1)
   do j=1,na
      k=k+1
      a(j)=buf(k)
   enddo

   if (imode(3) .lt. -1) then
     write (nsyso,2301) na, inow, a(1), a(2), a(na) 
 2301 format (/,1x, 'finda io retrievial exit ', 2i6, 3g14.7)
   endif

   return
   end subroutine finda

   subroutine scana(e,ip,np,na,ntape,buf,nbuf)
   !--------------------------------------------------------------------
   ! Search the loada/finda buffer for e.  Read in new buffers as
   ! necessary.  Value of ip returned corresponds to the location
   ! of e in the buffer, or, if e does not appear, to the nearest
   ! lower energy.
   !--------------------------------------------------------------------
   ! externals
   integer::ip,np,na,ntape,nbuf
   real(kr)::e,buf(nbuf)
   ! internals
   integer::nl,j,last,ibufl,indx,jp,ipn
   real(kr)::ef,el,en
   character(36)::strng

   if (ip.ne.0) call error('scana','initial ip ne 0.',' ')
   nl=nbuf/na
   j=1
   last=0
   j=1
   last=0
   110 continue
   ef=buf(1)
   ibufl=na*(nl-1)+1
   if (np.lt.ip+nl) ibufl=(np-ip-1)*na+1
   el=buf(ibufl)
   if (e.lt.ef.and.j.eq.1) go to 130
   if (e.lt.ef.and.j.gt.1) go to 140
   if (e.ge.ef.and.e.lt.el) go to 120
   ! e gt highest energy in this buffer.  read in a new buffer
   if (last.gt.0) write(strng,&
     '(''did not find energy '',1p,e12.4)') e
   if (last.gt.0) call error('scana',strng,' ')
   read(ntape) buf
   j=j+1
   ip=ip+nl
   if (np.le.ip+nl) last=1
   go to 110
   ! e is in this buffer.  determine ip
   120 continue
   ip=ip+1
   jp=ip-nl*(j-1)
   indx=na*(jp-1)+1
   jp=ip-nl*(j-1)
   indx=na*(jp-1)+1
   en=buf(indx)
   if (en.eq.e) go to 145
   if (en.lt.e) go to 120
   ip=ip-1
   go to 145
   130 continue
   ip=1
   go to 150
   ! first e in this buffer is higher than e
   ! backspace and read next lower buffer again
   140 continue
   backspace ntape
   backspace ntape
   read(ntape) buf
   ip=ip-1
   go to 150
   145 continue
   ipn=mod(ip,nl)
   if (ipn.ne.1) go to 150
   backspace ntape
   if (ip.lt.1) ip=1
   150 continue
   return
   end subroutine scana

   real(kr) function sigfig(x,ndig,idig)
   !--------------------------------------------------------------------
   ! Adjust x to have ndig signficant figures.  If idig is not zero,
   ! shade x up or down by idig in the last significant figure.
   !--------------------------------------------------------------------
   ! externals
   real(kr)::x
   integer::ndig,idig
   ! internals
   real(kr)::xx,aa
   integer::ipwr,ii
   real(kr),parameter::bias=1.0000000000001e0_kr
   real(kr),parameter::zero=0
   real(kr),parameter::ten=10

   xx=0
   if (x.ne.zero) then
      aa=log10(abs(x))
      ipwr=int(aa)
      if (aa.lt.zero) ipwr=ipwr-1
      ipwr=ndig-1-ipwr
      ii=nint(x*ten**ipwr+ten**(ndig-11))
      if (ii.ge.10**ndig) then
         ii=ii/10
         ipwr=ipwr-1
      endif
      ii=ii+idig
      xx=ii*ten**(-ipwr)
   endif
   xx=xx*bias
   sigfig=xx
   return
   end function sigfig

   subroutine a10(x,hx)
   !-------------------------------------------------------------------
   ! convert x to a 10-column format as a string in hx to provide
   ! more digits in some njoy listings without taking too much space.
   ! this allows 4, 5, or 6 significant figures to be printed where we
   ! previously had four.  hx has the forms +1.23456+6, +1.2345-38,
   ! or -1.234+308.  based on subroutine a11 for endf formats.
   !-------------------------------------------------------------------
   ! externals
   real(kr)::x
   character(10)::hx
   ! internals
   integer::n
   real(kr)::f,xx,ff
   character(1)::s
   real(kr),parameter::zero=0
   real(kr),parameter::one=1
   real(kr),parameter::ten=10
   real(kr),parameter::amil=.999999999e6_kr
   real(kr),parameter::up=.000000001e0_kr
   real(kr),parameter::top6=9.999995e0_kr
   real(kr),parameter::top5=9.99995e0_kr
   real(kr),parameter::top4=9.9995e0_kr

   !--check for zero as a special case
   if (x.ne.zero) go to 100
   f=0
   s='+'
   n=0
   go to 160

   !--construct six, five, and four sig-fig forms
  100 continue
   xx=x
   n=int(log10(abs(xx)))
   if (abs(xx).lt.one) go to 110
   ff=xx/ten**n
   s='+'
   if (iabs(n).lt.10.and.abs(ff).lt.top6) go to 120
   if (iabs(n).lt.100.and.iabs(n).ge.10.and.&
     abs(ff).lt.top5) go to 120
   if (iabs(n).ge.100.and.abs(ff).lt.top4) go to 120
   ff=ff/ten
   ff=ff+up
   n=n+1
   go to 120
  110 continue
   n=1-n
   ff=xx*ten**n
   s='-'
   if (iabs(n).lt.10.and.abs(ff).lt.top6) go to 120
   if (iabs(n).lt.100.and.iabs(n).ge.10.and.&
     abs(ff).lt.top5) go to 120
   if (iabs(n).ge.100.and.abs(ff).lt.top4) go to 120
   ff=ff/ten
   ff=ff+up
   n=n-1
   if (n.gt.0) go to 120
   s='+'
  120 continue
   f=ff
  160 continue
   if (iabs(n).lt.10) then
      write(hx,'(f8.5,a,i1)') f,s,n
   else if (iabs(n).lt.100) then
      write(hx,'(f7.4,a,i2)') f,s,n
   else
      write(hx,'(f6.3,a,i3)') f,s,n
   endif
   return
   end subroutine a10

end module util


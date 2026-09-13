module kdim
   !
   !-----------------------------------------------------------------------
   !
   ! kdi file support module.
   !
   ! this module provide the possibility to use a direct access (da) file
   ! to simulate memory like support. maxbuf static buffers are used.
   !
   !  contains:
   !  . kdiget    extract data from the kdi file.
   !  . kdiput    save data into the kdi file.
   !  . kdicl     terminate operations on the kdi file.
   !
   !  remark : file opening must be done using the open statement
   !           provided in the Fortran-90 standard:
   !           open(ifile,access='direct',recl=lrecl*iunit) where iunit
   !           is the number of memory units in a single precision word.
   !
   ! saved variables:
   !  nbuf           : number of buffers used. nbuf<=maxbuf
   !  ibufs          : index of the last buffer used.
   !  nrecs(ibuf)    : relative record number corresponding to the buffer.
   !  ifiles(ibuf)   : unit number corresponding to the buffer.
   !  iage(ibuf)     : age of the buffer.
   !  modif(ibuf)    : modification flag for the buffer.
   !  record(1,ibuf) : static buffer. dimension record(lrecl,maxbuf)
   !
   ! Author: Alain Hebert (original version -- 1994)
   !
   !-----------------------------------------------------------------------
   !
   implicit integer(a-z)
   private
   public :: kdiget,kdiput,kdicl
   integer,parameter :: maxbuf=10,lrecl=1024
   integer,save :: nbuf=0
   integer,save :: ibufs,nrecs(maxbuf),ifiles(maxbuf),iage(maxbuf)
   logical,save :: modif(maxbuf)
   integer,save :: record(lrecl,maxbuf)
   !
   interface kdiget
      !
      !----------------------------------------------------------------------
      !
      ! extract data from the kdi file.
      !
      ! input parameters:
      !  ifile  : unit number of the kdi file.
      !  iofset : offset, from start of kdi file where data is extracted.
      !  length : number of addressable units to read from file.
      !  irc    : return code (zero if no error is detected).
      !
      ! output parameter:
      !  data   : data to be copied from the kdi file.
      !
      !----------------------------------------------------------------------
      !
      module procedure kdiget_i,kdiget_r,kdiget_h,kdiget_d,kdiget_l,kdiget_c
   end interface
   !
   interface kdiput
      !
      !----------------------------------------------------------------------
      !
      ! save data into the kdi file.
      !
      ! input parameters:
      !  ifile  : unit number of the kdi file.
      !  data   : data to be copied from the kdi file.
      !  iofset : offset, from start of kdi file where data is written.
      !  length : number of addressable units to write into file.
      !  irc    : return code (zero if no error is detected).
      !
      !----------------------------------------------------------------------
      !
      module procedure kdiput_i,kdiput_r,kdiput_h,kdiput_d,kdiput_l,kdiput_c
   end interface
contains
   subroutine kdiget_i(ifile,data,iofset,length,irc)
   integer,dimension(:) :: data
   !
   if(lbound(data,1)/=1) then
      irc=-90
      return
   endif
   irc=0
   kbuf=0
   nrec=1+iofset/lrecl
   n=mod(iofset,lrecl)
   lmin=1
   10 ngro=min(lrecl+lmin-n-1,length)
   kdm=-1
   do iiii=1,nbuf
      ibuf=1+mod(iiii+ibufs-2,nbuf)
      if((ifile==ifiles(ibuf)).and.(nrec==nrecs(ibuf))) then
         do l=lmin,ngro
            n=n+1
            data(l)=record(n,ibuf)
         enddo
         do i=1,nbuf
            if(i/=ibuf) iage(i)=iage(i)+1
         enddo
         ibufs=ibuf
         if(ngro==length) return
         nrec=nrec+1
         n=0
         lmin=ngro+1
         ibufs=ibuf+1
         go to 10
      endif
      if(iage(ibuf)>kdm) then
         kdm=iage(ibuf)
         kbuf=ibuf
      endif
   enddo
   if(nbuf.lt.maxbuf) then
      nbuf=nbuf+1
      ibufs=nbuf
   else
      ibufs=kbuf
      if(modif(ibufs)) then
         ! -------------------------------------------------------
         write(ifiles(ibufs),rec=nrecs(ibufs),err=20,iostat=irc) &
         & (record(i,ibufs),i=1,lrecl)
         ! -------------------------------------------------------
      endif
   endif
   nrecs(ibufs)=nrec
   ifiles(ibufs)=ifile
   iage(ibufs)=0
   modif(ibufs)=.false.
   ! ----------------------------------------------------------------
   read(ifile,rec=nrec,err=20,iostat=irc) (record(i,ibufs),i=1,lrecl)
   ! ----------------------------------------------------------------
   go to 10
   20 return
   end subroutine kdiget_i
   !
   subroutine kdiget_r(ifile,data,iofset,length,irc)
   real,dimension(:) :: data
   !
   if(lbound(data,1)/=1) then
      irc=-91
      return
   endif
   irc=0
   kbuf=0
   nrec=1+iofset/lrecl
   n=mod(iofset,lrecl)
   lmin=1
   10 ngro=min(lrecl+lmin-n-1,length)
   kdm=-1
   do iiii=1,nbuf
      ibuf=1+mod(iiii+ibufs-2,nbuf)
      if((ifile==ifiles(ibuf)).and.(nrec==nrecs(ibuf))) then
         do l=lmin,ngro
            n=n+1
            data(l)=transfer(record(n,ibuf),0.0)
         enddo
         do i=1,nbuf
            if(i/=ibuf) iage(i)=iage(i)+1
         enddo
         ibufs=ibuf
         if(ngro==length) return
         nrec=nrec+1
         n=0
         lmin=ngro+1
         ibufs=ibuf+1
         go to 10
      endif
      if(iage(ibuf)>kdm) then
         kdm=iage(ibuf)
         kbuf=ibuf
      endif
   enddo
   if(nbuf.lt.maxbuf) then
      nbuf=nbuf+1
      ibufs=nbuf
   else
      ibufs=kbuf
      if(modif(ibufs)) then
         ! -------------------------------------------------------
         write(ifiles(ibufs),rec=nrecs(ibufs),err=20,iostat=irc) &
         & (record(i,ibufs),i=1,lrecl)
         ! -------------------------------------------------------
      endif
   endif
   nrecs(ibufs)=nrec
   ifiles(ibufs)=ifile
   iage(ibufs)=0
   modif(ibufs)=.false.
   ! ----------------------------------------------------------------
   read(ifile,rec=nrec,err=20,iostat=irc) (record(i,ibufs),i=1,lrecl)
   ! ----------------------------------------------------------------
   go to 10
   20 return
   end subroutine kdiget_r
   !
   subroutine kdiget_h(ifile,data,iofset,length,irc)
   character(len=*),dimension(:) :: data
   !
   if(lbound(data,1)/=1) then
      irc=-92
      return
   else if(len(data(1))/=4) then
      irc=-93
      return
   endif
   irc=0
   kbuf=0
   nrec=1+iofset/lrecl
   n=mod(iofset,lrecl)
   lmin=1
   10 ngro=min(lrecl+lmin-n-1,length)
   kdm=-1
   do iiii=1,nbuf
      ibuf=1+mod(iiii+ibufs-2,nbuf)
      if((ifile==ifiles(ibuf)).and.(nrec==nrecs(ibuf))) then
         do l=lmin,ngro
            n=n+1
            data(l)=transfer(record(n,ibuf),'aaaa')
         enddo
         do i=1,nbuf
            if(i/=ibuf) iage(i)=iage(i)+1
         enddo
         ibufs=ibuf
         if(ngro==length) return
         nrec=nrec+1
         n=0
         lmin=ngro+1
         ibufs=ibuf+1
         go to 10
      endif
      if(iage(ibuf)>kdm) then
         kdm=iage(ibuf)
         kbuf=ibuf
      endif
   enddo
   if(nbuf.lt.maxbuf) then
      nbuf=nbuf+1
      ibufs=nbuf
   else
      ibufs=kbuf
      if(modif(ibufs)) then
         ! -------------------------------------------------------
         write(ifiles(ibufs),rec=nrecs(ibufs),err=20,iostat=irc) &
         & (record(i,ibufs),i=1,lrecl)
         ! -------------------------------------------------------
      endif
   endif
   nrecs(ibufs)=nrec
   ifiles(ibufs)=ifile
   iage(ibufs)=0
   modif(ibufs)=.false.
   ! ----------------------------------------------------------------
   read(ifile,rec=nrec,err=20,iostat=irc) (record(i,ibufs),i=1,lrecl)
   ! ----------------------------------------------------------------
   go to 10
   20 return
   end subroutine kdiget_h
   !
   subroutine kdiget_d(ifile,data,iofset,length,irc)
   double precision,dimension(:) :: data
   !
   if(lbound(data,1)/=1) then
      irc=-94
      return
   endif
   irc=0
   kbuf=0
   nrec=1+iofset/lrecl
   n=mod(iofset,lrecl)
   lmin=1
   10 ngro=min(lrecl+lmin-n-1,length)
   kdm=-1
   do iiii=1,nbuf
      ibuf=1+mod(iiii+ibufs-2,nbuf)
      if((ifile==ifiles(ibuf)).and.(nrec==nrecs(ibuf))) then
         do l=(lmin+1)/2,(ngro+1)/2
            n=n+2
            data(l)=transfer(record(n-1:n,ibuf),0.0d0)
         enddo
         do i=1,nbuf
            if(i/=ibuf) iage(i)=iage(i)+1
         enddo
         ibufs=ibuf
         if(ngro==length) return
         nrec=nrec+1
         n=0
         lmin=ngro+1
         ibufs=ibuf+1
         go to 10
      endif
      if(iage(ibuf)>kdm) then
         kdm=iage(ibuf)
         kbuf=ibuf
      endif
   enddo
   if(nbuf.lt.maxbuf) then
      nbuf=nbuf+1
      ibufs=nbuf
   else
      ibufs=kbuf
      if(modif(ibufs)) then
         ! -------------------------------------------------------
         write(ifiles(ibufs),rec=nrecs(ibufs),err=20,iostat=irc) &
         & (record(i,ibufs),i=1,lrecl)
         ! -------------------------------------------------------
      endif
   endif
   nrecs(ibufs)=nrec
   ifiles(ibufs)=ifile
   iage(ibufs)=0
   modif(ibufs)=.false.
   ! ----------------------------------------------------------------
   read(ifile,rec=nrec,err=20,iostat=irc) (record(i,ibufs),i=1,lrecl)
   ! ----------------------------------------------------------------
   go to 10
   20 return
   end subroutine kdiget_d
   !
   subroutine kdiget_l(ifile,data,iofset,length,irc)
   logical,dimension(:) :: data
   !
   if(lbound(data,1)/=1) then
      irc=-95
      return
   endif
   irc=0
   kbuf=0
   nrec=1+iofset/lrecl
   n=mod(iofset,lrecl)
   lmin=1
   10 ngro=min(lrecl+lmin-n-1,length)
   kdm=-1
   do iiii=1,nbuf
      ibuf=1+mod(iiii+ibufs-2,nbuf)
      if((ifile==ifiles(ibuf)).and.(nrec==nrecs(ibuf))) then
         do l=lmin,ngro
            n=n+1
            data(l)=transfer(record(n,ibuf),.true.)
         enddo
         do i=1,nbuf
            if(i/=ibuf) iage(i)=iage(i)+1
         enddo
         ibufs=ibuf
         if(ngro==length) return
         nrec=nrec+1
         n=0
         lmin=ngro+1
         ibufs=ibuf+1
         go to 10
      endif
      if(iage(ibuf)>kdm) then
         kdm=iage(ibuf)
         kbuf=ibuf
      endif
   enddo
   if(nbuf.lt.maxbuf) then
      nbuf=nbuf+1
      ibufs=nbuf
   else
      ibufs=kbuf
      if(modif(ibufs)) then
         ! -------------------------------------------------------
         write(ifiles(ibufs),rec=nrecs(ibufs),err=20,iostat=irc) &
         & (record(i,ibufs),i=1,lrecl)
         ! -------------------------------------------------------
      endif
   endif
   nrecs(ibufs)=nrec
   ifiles(ibufs)=ifile
   iage(ibufs)=0
   modif(ibufs)=.false.
   ! ----------------------------------------------------------------
   read(ifile,rec=nrec,err=20,iostat=irc) (record(i,ibufs),i=1,lrecl)
   ! ----------------------------------------------------------------
   go to 10
   20 return
   end subroutine kdiget_l
   !
   subroutine kdiget_c(ifile,data,iofset,length,irc)
   complex,dimension(:) :: data
   !
   if(lbound(data,1)/=1) then
      irc=-94
      return
   endif
   irc=0
   kbuf=0
   nrec=1+iofset/lrecl
   n=mod(iofset,lrecl)
   lmin=1
   10 ngro=min(lrecl+lmin-n-1,length)
   kdm=-1
   do iiii=1,nbuf
      ibuf=1+mod(iiii+ibufs-2,nbuf)
      if((ifile==ifiles(ibuf)).and.(nrec==nrecs(ibuf))) then
         do l=(lmin+1)/2,(ngro+1)/2
            n=n+2
            data(l)=transfer(record(n-1:n,ibuf),(0.0,0.0))
         enddo
         do i=1,nbuf
            if(i/=ibuf) iage(i)=iage(i)+1
         enddo
         ibufs=ibuf
         if(ngro==length) return
         nrec=nrec+1
         n=0
         lmin=ngro+1
         ibufs=ibuf+1
         go to 10
      endif
      if(iage(ibuf)>kdm) then
         kdm=iage(ibuf)
         kbuf=ibuf
      endif
   enddo
   if(nbuf.lt.maxbuf) then
      nbuf=nbuf+1
      ibufs=nbuf
   else
      ibufs=kbuf
      if(modif(ibufs)) then
         ! -------------------------------------------------------
         write(ifiles(ibufs),rec=nrecs(ibufs),err=20,iostat=irc) &
         & (record(i,ibufs),i=1,lrecl)
         ! -------------------------------------------------------
      endif
   endif
   nrecs(ibufs)=nrec
   ifiles(ibufs)=ifile
   iage(ibufs)=0
   modif(ibufs)=.false.
   ! ----------------------------------------------------------------
   read(ifile,rec=nrec,err=20,iostat=irc) (record(i,ibufs),i=1,lrecl)
   ! ----------------------------------------------------------------
   go to 10
   20 return
   end subroutine kdiget_c
   !
   subroutine kdiput_i(ifile,data,iofset,length,irc)
   integer,dimension(:) :: data
   character(len=72) :: text72
   !
   if(lbound(data,1)/=1) then
      irc=-97
      return
   endif
   kbuf=0
   irc=0
   nrec=1+iofset/lrecl
   n=mod(iofset,lrecl)
   lmin=1
   10 ngro=min(lrecl+lmin-n-1,length)
   kdm=-1
   do iiii=1,nbuf
      ibuf=1+mod(iiii+ibufs-2,nbuf)
      if((ifile==ifiles(ibuf)).and.(nrec==nrecs(ibuf))) then
         do l=lmin,ngro
            n=n+1
            record(n,ibuf)=data(l)
         enddo
         do i=1,maxbuf
            if(i/=ibuf) iage(i)=iage(i)+1
         enddo
         modif(ibuf)=.true.
         ibufs=ibuf
         if(ngro==length) return
         nrec=nrec+1
         n=0
         lmin=ngro+1
         ibufs=ibuf+1
         go to 10
      endif
      if(iage(ibuf)>kdm) then
         kdm=iage(ibuf)
         kbuf=ibuf
      endif
   enddo
   if(nbuf.lt.maxbuf) then
      nbuf=nbuf+1
      ibufs=nbuf
   else
      ibufs=kbuf
      if(modif(ibufs)) then
         ! -------------------------------------------------------
         write(ifiles(ibufs),rec=nrecs(ibufs),err=30,iostat=irc) &
         & (record(i,ibufs),i=1,lrecl)
         ! -------------------------------------------------------
      endif
   endif
   nrecs(ibufs)=nrec
   ifiles(ibufs)=ifile
   iage(ibufs)=0
   modif(ibufs)=.false.
   if((n/=0).or.(lrecl+lmin-n-1>length)) then
      ! ---------------------------------------------------------------
      read(ifile,rec=nrec,err=20,iostat=ii) (record(i,ibufs),i=1,lrecl)
      ! ---------------------------------------------------------------
   endif
   go to 10
   20 inquire(ifile,name=text72,recl=nn)
   close(ifile,status='keep')
   open(ifile,file=text72,access='direct',recl=nn,status='old')
   go to 10
   30 return
   end subroutine kdiput_i
   !
   subroutine kdiput_r(ifile,data,iofset,length,irc)
   real,dimension(:) :: data
   character(len=72) :: text72
   !
   if(lbound(data,1)/=1) then
      irc=-98
      return
   endif
   kbuf=0
   irc=0
   nrec=1+iofset/lrecl
   n=mod(iofset,lrecl)
   lmin=1
   10 ngro=min(lrecl+lmin-n-1,length)
   kdm=-1
   do iiii=1,nbuf
      ibuf=1+mod(iiii+ibufs-2,nbuf)
      if((ifile==ifiles(ibuf)).and.(nrec==nrecs(ibuf))) then
         do l=lmin,ngro
            n=n+1
            record(n,ibuf)=transfer(data(l),0)
         enddo
         do i=1,maxbuf
            if(i/=ibuf) iage(i)=iage(i)+1
         enddo
         modif(ibuf)=.true.
         ibufs=ibuf
         if(ngro==length) return
         nrec=nrec+1
         n=0
         lmin=ngro+1
         ibufs=ibuf+1
         go to 10
      endif
      if(iage(ibuf)>kdm) then
         kdm=iage(ibuf)
         kbuf=ibuf
      endif
   enddo
   if(nbuf.lt.maxbuf) then
      nbuf=nbuf+1
      ibufs=nbuf
   else
      ibufs=kbuf
      if(modif(ibufs)) then
         ! -------------------------------------------------------
         write(ifiles(ibufs),rec=nrecs(ibufs),err=30,iostat=irc) &
         & (record(i,ibufs),i=1,lrecl)
         ! -------------------------------------------------------
      endif
   endif
   nrecs(ibufs)=nrec
   ifiles(ibufs)=ifile
   iage(ibufs)=0
   modif(ibufs)=.false.
   if((n/=0).or.(lrecl+lmin-n-1>length)) then
      ! ---------------------------------------------------------------
      read(ifile,rec=nrec,err=20,iostat=ii) (record(i,ibufs),i=1,lrecl)
      ! ---------------------------------------------------------------
   endif
   go to 10
   20 inquire(ifile,name=text72,recl=nn)
   close(ifile,status='keep')
   open(ifile,file=text72,access='direct',recl=nn,status='old')
   go to 10
   30 return
   end subroutine kdiput_r
   !
   subroutine kdiput_h(ifile,data,iofset,length,irc)
   character(len=*),dimension(:) :: data
   character(len=72) :: text72
   !
   if(lbound(data,1)/=1) then
      irc=-99
      return
   else if(len(data(1))/=4) then
      irc=-100
      return
   endif
   kbuf=0
   irc=0
   nrec=1+iofset/lrecl
   n=mod(iofset,lrecl)
   lmin=1
   10 ngro=min(lrecl+lmin-n-1,length)
   kdm=-1
   do iiii=1,nbuf
      ibuf=1+mod(iiii+ibufs-2,nbuf)
      if((ifile==ifiles(ibuf)).and.(nrec==nrecs(ibuf))) then
         do l=lmin,ngro
            n=n+1
            record(n,ibuf)=transfer(data(l),0)
         enddo
         do i=1,maxbuf
            if(i/=ibuf) iage(i)=iage(i)+1
         enddo
         modif(ibuf)=.true.
         ibufs=ibuf
         if(ngro==length) return
         nrec=nrec+1
         n=0
         lmin=ngro+1
         ibufs=ibuf+1
         go to 10
      endif
      if(iage(ibuf)>kdm) then
         kdm=iage(ibuf)
         kbuf=ibuf
      endif
   enddo
   if(nbuf.lt.maxbuf) then
      nbuf=nbuf+1
      ibufs=nbuf
   else
      ibufs=kbuf
      if(modif(ibufs)) then
         ! -------------------------------------------------------
         write(ifiles(ibufs),rec=nrecs(ibufs),err=30,iostat=irc) &
         & (record(i,ibufs),i=1,lrecl)
         ! -------------------------------------------------------
      endif
   endif
   nrecs(ibufs)=nrec
   ifiles(ibufs)=ifile
   iage(ibufs)=0
   modif(ibufs)=.false.
   if((n/=0).or.(lrecl+lmin-n-1>length)) then
      ! ---------------------------------------------------------------
      read(ifile,rec=nrec,err=20,iostat=ii) (record(i,ibufs),i=1,lrecl)
      ! ---------------------------------------------------------------
   endif
   go to 10
   20 inquire(ifile,name=text72,recl=nn)
   close(ifile,status='keep')
   open(ifile,file=text72,access='direct',recl=nn,status='old')
   go to 10
   30 return
   end subroutine kdiput_h
   !
   subroutine kdiput_d(ifile,data,iofset,length,irc)
   double precision,dimension(:) :: data
   character(len=72) :: text72
   !
   if(lbound(data,1)/=1) then
      irc=-101
      return
   endif
   kbuf=0
   irc=0
   nrec=1+iofset/lrecl
   n=mod(iofset,lrecl)
   lmin=1
   10 ngro=min(lrecl+lmin-n-1,length)
   kdm=-1
   do iiii=1,nbuf
      ibuf=1+mod(iiii+ibufs-2,nbuf)
      if((ifile==ifiles(ibuf)).and.(nrec==nrecs(ibuf))) then
         do l=(lmin+1)/2,(ngro+1)/2
            n=n+2
            record(n-1:n,ibuf)=transfer(data(l),(/0/))
         enddo
         do i=1,maxbuf
            if(i/=ibuf) iage(i)=iage(i)+1
         enddo
         modif(ibuf)=.true.
         ibufs=ibuf
         if(ngro==length) return
         nrec=nrec+1
         n=0
         lmin=ngro+1
         ibufs=ibuf+1
         go to 10
      endif
      if(iage(ibuf)>kdm) then
         kdm=iage(ibuf)
         kbuf=ibuf
      endif
   enddo
   if(nbuf.lt.maxbuf) then
      nbuf=nbuf+1
      ibufs=nbuf
   else
      ibufs=kbuf
      if(modif(ibufs)) then
         ! -------------------------------------------------------
         write(ifiles(ibufs),rec=nrecs(ibufs),err=30,iostat=irc) &
         & (record(i,ibufs),i=1,lrecl)
         ! -------------------------------------------------------
      endif
   endif
   nrecs(ibufs)=nrec
   ifiles(ibufs)=ifile
   iage(ibufs)=0
   modif(ibufs)=.false.
   if((n/=0).or.(lrecl+lmin-n-1>length)) then
      ! ---------------------------------------------------------------
      read(ifile,rec=nrec,err=20,iostat=ii) (record(i,ibufs),i=1,lrecl)
      ! ---------------------------------------------------------------
   endif
   go to 10
   20 inquire(ifile,name=text72,recl=nn)
   close(ifile,status='keep')
   open(ifile,file=text72,access='direct',recl=nn,status='old')
   go to 10
   30 return
   end subroutine kdiput_d
   !
   subroutine kdiput_l(ifile,data,iofset,length,irc)
   logical,dimension(:) :: data
   character(len=72) :: text72
   !
   if(lbound(data,1)/=1) then
      irc=-102
      return
   endif
   kbuf=0
   irc=0
   nrec=1+iofset/lrecl
   n=mod(iofset,lrecl)
   lmin=1
   10 ngro=min(lrecl+lmin-n-1,length)
   kdm=-1
   do iiii=1,nbuf
      ibuf=1+mod(iiii+ibufs-2,nbuf)
      if((ifile==ifiles(ibuf)).and.(nrec==nrecs(ibuf))) then
         do l=lmin,ngro
            n=n+1
            record(n,ibuf)=transfer(data(l),0)
         enddo
         do i=1,maxbuf
            if(i/=ibuf) iage(i)=iage(i)+1
         enddo
         modif(ibuf)=.true.
         ibufs=ibuf
         if(ngro==length) return
         nrec=nrec+1
         n=0
         lmin=ngro+1
         ibufs=ibuf+1
         go to 10
      endif
      if(iage(ibuf)>kdm) then
         kdm=iage(ibuf)
         kbuf=ibuf
      endif
   enddo
   if(nbuf.lt.maxbuf) then
      nbuf=nbuf+1
      ibufs=nbuf
   else
      ibufs=kbuf
      if(modif(ibufs)) then
         ! -------------------------------------------------------
         write(ifiles(ibufs),rec=nrecs(ibufs),err=30,iostat=irc) &
         & (record(i,ibufs),i=1,lrecl)
         ! -------------------------------------------------------
      endif
   endif
   nrecs(ibufs)=nrec
   ifiles(ibufs)=ifile
   iage(ibufs)=0
   modif(ibufs)=.false.
   if((n/=0).or.(lrecl+lmin-n-1>length)) then
      ! ---------------------------------------------------------------
      read(ifile,rec=nrec,err=20,iostat=ii) (record(i,ibufs),i=1,lrecl)
      ! ---------------------------------------------------------------
   endif
   go to 10
   20 inquire(ifile,name=text72,recl=nn)
   close(ifile,status='keep')
   open(ifile,file=text72,access='direct',recl=nn,status='old')
   go to 10
   30 return
   end subroutine kdiput_l
   !
   subroutine kdiput_c(ifile,data,iofset,length,irc)
   complex,dimension(:) :: data
   character(len=72) :: text72
   !
   if(lbound(data,1)/=1) then
      irc=-103
      return
   endif
   kbuf=0
   irc=0
   nrec=1+iofset/lrecl
   n=mod(iofset,lrecl)
   lmin=1
   10 ngro=min(lrecl+lmin-n-1,length)
   kdm=-1
   do iiii=1,nbuf
      ibuf=1+mod(iiii+ibufs-2,nbuf)
      if((ifile==ifiles(ibuf)).and.(nrec==nrecs(ibuf))) then
         do l=(lmin+1)/2,(ngro+1)/2
            n=n+2
            record(n-1:n,ibuf)=transfer(data(l),(/0/))
         enddo
         do i=1,maxbuf
            if(i/=ibuf) iage(i)=iage(i)+1
         enddo
         modif(ibuf)=.true.
         ibufs=ibuf
         if(ngro==length) return
         nrec=nrec+1
         n=0
         lmin=ngro+1
         ibufs=ibuf+1
         go to 10
      endif
      if(iage(ibuf)>kdm) then
         kdm=iage(ibuf)
         kbuf=ibuf
      endif
   enddo
   if(nbuf.lt.maxbuf) then
      nbuf=nbuf+1
      ibufs=nbuf
   else
      ibufs=kbuf
      if(modif(ibufs)) then
         ! -------------------------------------------------------
         write(ifiles(ibufs),rec=nrecs(ibufs),err=30,iostat=irc) &
         & (record(i,ibufs),i=1,lrecl)
         ! -------------------------------------------------------
      endif
   endif
   nrecs(ibufs)=nrec
   ifiles(ibufs)=ifile
   iage(ibufs)=0
   modif(ibufs)=.false.
   if((n/=0).or.(lrecl+lmin-n-1>length)) then
      ! ---------------------------------------------------------------
      read(ifile,rec=nrec,err=20,iostat=ii) (record(i,ibufs),i=1,lrecl)
      ! ---------------------------------------------------------------
   endif
   go to 10
   20 inquire(ifile,name=text72,recl=nn)
   close(ifile,status='keep')
   open(ifile,file=text72,access='direct',recl=nn,status='old')
   go to 10
   30 return
   end subroutine kdiput_c
   !
   subroutine kdicl(ifile,istatu,jrc)
   kbuf=0
   if(istatu==1) then
      do ibuf=1,nbuf
         if(modif(ibuf).and.(ifile==ifiles(ibuf))) then
            ! -----------------------------------------------
            write(ifile,rec=nrecs(ibuf),err=130,iostat=jrc) &
            & (record(i,ibuf),i=1,lrecl)
            ! -----------------------------------------------
            modif(ibuf)=.false.
         endif
         if(ifile==ifiles(ibuf)) ifiles(ibuf)=0
      enddo
      close(ifile,err=130,status='keep',iostat=jrc)
   else if(istatu==2) then
      do ibuf=1,nbuf
         if(ifile==ifiles(ibuf)) then
            iage(ibuf)=999999
            modif(ibuf)=.false.
            ifiles(ibuf)=0
         endif
      enddo
      close(ifile,err=130,status='delete',iostat=jrc)
   else
      jrc=-999
   endif
   130 return
   end subroutine kdicl
end module kdim


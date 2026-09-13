module filem
   !
   !-----------------------------------------------------------------------
   !
   ! file allocation/release support module.
   !
   ! allocate and release file units associated to a given file name.
   ! direct access (kdi), sequential (formatted or not) and direct
   ! access (da) files are permitted
   !
   !  contains:
   !  . kdropn    open file and allocate unit number.
   !  . kdrcls    close file and release unit number.
   !  . kdrfst    return information about a file.
   !
   ! saved variables:
   !  nbtape   : number of units available                     i
   !  nreser   : number of reserved names                      i
   !  ndummy   : number of dummy units                         i
   !  creser   : reserved names                     c(nreser)*12
   !  cdummy   : names of dummy files               c(ndummy)*12
   !  ctapen   : names allocated to units           c(nbtape)*16
   !  itapet   : type of files for all units           i(nbtape)
   !             =0  no type associated to file
   !             =iutype  file of type iutype opened
   !             =-iutype  file of type iutype closed
   !  ilenda   : length of da file                     i(nbtape)
   !             =0       not a da file
   !             =lrda    length of da file
   !
   ! Author: Alain Hebert (original version -- 1994)
   !
   !-----------------------------------------------------------------------
   !
   use kdim
   private
   public :: kdropn,kdrcls,kdrfst
   integer,parameter :: nbtape=63,nreser=2,ndummy=4,nb=nbtape-6,lrecld=1024
   character(len=16),save,dimension(nbtape) :: ctapen=               &
    (/ ('        ',i=1,4),'ft05f001','ft06f001',                     &
    ('        ',i=1,nb) /)
   character(len=12),save,dimension(ndummy) :: cdummy=               &
    (/ 'dummyda','dummysq','dummyca','dummyin' /)
   character(len=12),save,dimension(nreser) :: creser=               &
    (/ 'ft05f001','ft06f001' /)
   character(len=22),save,dimension(ndummy) :: ctype=                &
    (/ 'direct access kdi     ','sequential unformatted',            &
    'sequential character  ','direct access da      ' /)
   integer,save,dimension(nbtape) :: itapet=0,ilenda=0
   !
   interface kdrfst
      !
      !-----------------------------------------------------------------------
      !
      ! return the unit allocated to file name or the file name allocated
      ! to unit. also return the iutype of that file.
      !
      ! input parameters:
      !  cuname/ifile : character file name or file unit number.
      !
      ! output parameters (optional):
      !  name  : character file name.
      !  unit  : file unit number.
      !  ldra  : number of words in da file. returned only if iutype=4.
      !  type  : error flag or iutype of the file.
      !          > 0  file iutype
      !          =-2  file not oppened
      !          =-3  this file has been opened by routines other than kdropn
      !          =-4  file unit is reserved (5,6)
      !          =-5  illegal unit number
      !
      !-----------------------------------------------------------------------
      !
      module procedure kdrfst_h,kdrfst_i
   end interface
   !
contains
   !
   function kdropn(osname,iactio,iutype,lrda)
   !
   !-----------------------------------------------------------------------
   !
   ! open file and allocate unit number. allocate unit number to file
   ! name if unit is already opened, returns unit number.
   !
   ! input parameters:
   !     osname   : filename                                     c
   !                if osname=' ', use a default name
   !     iactio   : action on file                               i
   !                =0 to create a new file
   !                =1 to access and modify an existing file
   !                =2 to access an existing file in read-only mode
   !                =3 unknown
   !     iutype   : file type                                    i
   !                =1  kdi file type
   !                =2  sequential unformatted
   !                =3  sequential formatted
   !                =4  direct access (da) unformatted file
   !     ldra     : number of words in da file                   i
   !                required for  iutype=4 only
   !
   ! output parameter:
   !     kdropn   : unit number/error status                     i
   !                >0  unit number allocated
   !               =-1  no more unit available
   !               =-2  file type requested inconsistent with file type
   !                    of this file
   !               =-3  this file has been already opened
   !               =-4  file name is reserved
   !                      sysin,ft06f00,ft07f00
   !               =-5  illegal file type 0 < iutype < 4
   !               =-6  error on open of kdi file
   !               =-7  error on open of unformatted sequential file
   !               =-8  error on open of formatted sequential file
   !               =-9  error on open of direct access file
   !               =-10 invalid number of word in direct access file
   !                    lrda must be > 0 and given iif iutype=4
   !               =-11 da file length inconsistent with previous length
   !                    for this file
   !
   !-----------------------------------------------------------------------
   !
   character(len=*) :: osname
   integer,optional :: lrda
   character(len=16) :: crdnam
   character(len=12) :: cform,cstatu
   logical :: lfilop
   integer,dimension(:),allocatable :: data
   !----
   !  check if iutype is valid
   !----
   kdropn=-9999
   if((iutype>4).or.(iutype<1)) then
     kdropn=-5
     return
   endif
   !----
   !  check if lrda is valid
   !----
   if(present(lrda)) then
     if((iutype/=4).or.(lrda<1)) then
       kdropn=-10
       return
     endif
   else
     if(iutype==4) then
       kdropn=-10
       return
     endif
   endif
!----
!  check if file name not forbidden
!----
   do ireser=1,nreser
     if(osname==creser(ireser)) then
       kdropn=-4
       return
     endif
   enddo
   !----
   !  check for dummy file name/allocate dummy file name if requested
   !----
   do idummy=1,ndummy
     if(osname==cdummy(idummy)) then
       if(idummy/=iutype) then
         kdropn=-2
         return
       endif
     endif
   enddo
   if(osname==' ') then
      crdnam=cdummy(iutype)
   else
      crdnam=osname
   endif
   !----
   !  check if file opened/permitted
   !----
   inquire(file=crdnam,opened=lfilop,number=iboucl)
   if(lfilop) then
     kdropn=-3
     return
   endif
   !----
   !  look for a new unit number for file crdnam
   !  look if file name is in table of previously allocated units
   !----
   do jboucl=nbtape,1,-1
     if(ctapen(jboucl)==crdnam) then
       iboucl=jboucl
       if(iutype/=-itapet(iboucl)) then
         kdropn=-2
         return
       else if(iutype==4) then
         kdropn=-11
         return
         if(lrda/=ilenda(iboucl)) then
           kdropn=-11
           return
         endif
       endif
       go to 121
     endif
   enddo
   !----
   !  look for never allocated unit location
   !----
   do jboucl=nbtape,1,-1
     if(ctapen(jboucl)==' ') then
       iboucl=jboucl
       inquire(unit=iboucl,opened=lfilop)
       if(.not.lfilop) go to 121
     endif
   enddo
   !----
   !  look for a unit previously allocated but released.
   !----
   do jboucl=nbtape,1,-1
     if(itapet(jboucl)<0) then
       iboucl=jboucl
       go to 121
     endif
   enddo
   !----
   !  error - no unit number available
   !----
   kdropn=-1
   return
   !
   121 lrda2=0
   if(iutype==1) then
      !----
      !  open kdi file
      !----
     kdropn=-6
     lrda2=lrecld
     cform='unformatted'
   else if(iutype==2) then
     !----
     !  open sequential unformatted file
     !----
     kdropn=-7
     cform='unformatted'
   else if(iutype==3) then
     !----
     !  open sequential formatted file
     !----
      kdropn=-8
      cform='formatted'
   else if(iutype==4) then
     !----
     !  open da file
     !----
     kdropn=-9
     lrda2=lrda
     cform='unformatted'
   endif
   if(iactio==0) then
     cstatu='new'
   else if(iactio==1) then
     cstatu='old'
   else if(iactio==2) then
     cstatu='old'
   else
     cstatu='unknown'
   endif
   if(((iutype==1).or.(iutype==4)).and.(iactio==2)) then
     allocate(data(lrda2))
     data=0
     inquire(iolength=lrecl) (data(i),i=1,lrda2)
     deallocate(data)
     open(unit=iboucl,file=crdnam,err=9000,iostat=iercod,form=cform, &
     access='direct',recl=lrecl,status=cstatu,action='read')
   else if((iutype==1).or.(iutype==4)) then
     allocate(data(lrda2))
     data=0
     inquire(iolength=lrecl) (data(i),i=1,lrda2)
     deallocate(data)
     open(unit=iboucl,file=crdnam,err=9000,iostat=iercod,form=cform, &
     access='direct',recl=lrecl,status=cstatu)
   else if(((iutype==2).or.(iutype==3)).and.(iactio==2)) then
     open(unit=iboucl,file=crdnam,err=9000,iostat=iercod,form=cform, &
     access='sequential',status=cstatu,action='read')
     rewind(iboucl)
   else if((iutype==2).or.(iutype==3)) then
     open(unit=iboucl,file=crdnam,err=9000,iostat=iercod,form=cform, &
     access='sequential',status=cstatu)
     rewind(iboucl)
   endif
   kdropn=iboucl
   ctapen(iboucl)=crdnam
   itapet(iboucl)=iutype
   if(iutype==4) ilenda(iboucl)=lrda
   return
   9000 write(6,8000) crdnam, ctype(iutype),iercod
   !----
   !  error format
   !----
   8000 format('1',5x,'error in opening of file in kdropn'/ &
   6x,'file name=',a7/6x,'file type=',a22/                  &
   6x,'   iercod=',i7)
   return
   end function kdropn
   !
   function kdrcls(itapno,iactio)
   !
   !-----------------------------------------------------------------------
   !
   ! close file and release unit number.
   !
   ! input parameters:
   !  itapno   : unit number                                  i
   !            =0    close all units
   !             > 0    unit to close
   !  iactio   : action on file                               i
   !            =1    to keep the file;
   !            =2    to delete the file.
   !
   ! output parameter:
   !  kdrcls   : error status                                 i
   !            = 0  unit closed
   !            =-2  file not oppened
   !            =-3  this file has been opened by routines
   !                   other than kdropn
   !            =-4  file unit is reserved (5,6)
   !            =-5  illegal unit number
   !            =-6  error on close of kdi file
   !            =-7  error on close of unformatted
   !                   sequential file
   !            =-8  error on close of formatted sequential file
   !            =-9  error on close of da file
   !            =-10  invalid close action iactio=1,2 permitted only
   !
   !-----------------------------------------------------------------------
   !
   logical :: lfilop
   !
   if(itapno==0) then
   !----
   !  unit=0 close all opened units
   !----
     idebtp=1
     ifintp=nbtape
     go to 200
   else if((itapno>nbtape).or.(itapno<0)) then
   !----
   !  invalid unit number
   !----
     kdrcls=-5
     return
   endif
   inquire(unit=itapno,opened=lfilop)
   if(lfilop) then
     if(itapet(itapno)==0) then
   !----
   !  reserved unit
   !----
       kdrcls=-4
       return
     endif
   !----
   !  unit permitted
   !----
     idebtp=itapno
     ifintp=itapno
     go to 200
   else
   !----
   !  unit not opened
   !----
     kdrcls=-2
     return
   endif
   200 do iboucl=idebtp,ifintp
     if(itapet(iboucl)==1) then
   !----
   !  close kdi type file
   !----
       kdrcls=-6
       call kdicl(iboucl,iactio,iercod)
       if(iercod/=0) go to 9001
       if(iactio==1) then
         itapet(iboucl)=-itapet(iboucl)
       else if(iactio==2) then
         itapet(iboucl)=0
         ctapen(iboucl)=' '
       else
         kdrcls=-10
       endif
     else if(itapet(iboucl)>0) then
   !----
   !  close other types of file
   !----
       kdrcls=-(5+itapet(iboucl))
       if(iactio==1) then
         close(iboucl,iostat=iercod,status='keep',err=9001)
         itapet(iboucl)=-itapet(iboucl)
       else if(iactio==2) then
         close(iboucl,iostat=iercod,status='delete',err=9001)
         itapet(iboucl)=0
         ctapen(iboucl)=' '
         ilenda(iboucl)=0
       else
         kdrcls=-10
       endif
       if(iercod/=0) go to 9001
     endif
   enddo
   kdrcls=0
   return
   9001 write(6,8001) iboucl,ctapen(iboucl),ctype(itapet(iboucl)),iercod
   !----
   !  error format
   !----
   8001 format('1',5x,'error in close of file in kdropn'/            &
   6x,'unit nb. =',i10/6x,'file name=',a7/6x,'file type=',a22/       &
   6x,'   iercod=',i10 )
   return
   end function kdrcls
   !
   subroutine kdrfst_h(cuname,name,unit,type,lrda)
   character(len=*) :: cuname
   character(len=*),optional :: name
   integer,optional :: unit,type,lrda
   character(len=16) :: crdnam
   logical :: lfilop
   !----
   !  return file number and type
   !----
   if(present(name)) name=cuname
   if(present(lrda)) lrda=0
   crdnam=cuname
   inquire(file=crdnam,opened=lfilop,number=iunit)
   if(present(unit)) unit=iunit
   if(lfilop) then
      if((iunit<1).or.(iunit>nbtape)) then
        if(present(type)) type=-5
      else if((ctapen(iunit)==' ').or.(itapet(iunit)<=0)) then
        if(present(type)) type=-3
      else
        if(present(type)) type=itapet(iunit)
        if((itapet(iunit)==4).and.present(lrda)) lrda=ilenda(iunit)
      endif
   else
      if(present(type)) type=-2
   endif
   return
   end subroutine kdrfst_h
   !
   subroutine kdrfst_i(ifile,name,unit,type,lrda)
   character(len=*),optional :: name
   integer,optional :: unit,type,lrda
   logical :: lfilop
   !
   if(present(name)) name=' '
   if(present(unit)) unit=ifile
   if(present(lrda)) lrda=0
   !----
   !  illegal unit
   !----
   if((ifile<1).or.(ifile>nbtape)) then
       if(present(type)) type=-5
       return
   endif
   !----
   !  return file name and type
   !----
   inquire(ifile,opened=lfilop)
   if(lfilop) then
      if((ctapen(ifile)==' ').or.(itapet(ifile)<=0)) then
        if(present(type)) type=-3
      else
        if(present(name)) name=ctapen(ifile)
        if(present(type)) type=itapet(ifile)
        if((itapet(ifile)==4).and.present(lrda)) lrda=ilenda(ifile)
      endif
   else
      if(present(type)) type=-2
   endif
   return
   end subroutine kdrfst_i
end module filem


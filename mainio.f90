module mainio
   ! Fortran unit numbers for system I/O.  These are conventional
   ! values, but they might have to be changed for some systems.
   ! The value nsyso is used to open a listing file in njoy.
   implicit none
   private
   integer,public::nsysi=5
   integer,public::nsyso=7
   integer,public::nsyse=6
end module mainio


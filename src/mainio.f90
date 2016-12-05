module mainio
   ! Fortran unit numbers for system I/O.  These are conventional
   ! values, but they might have to be changed for some systems.
   ! The value nsyso is used to open a listing file in njoy.
  use iso_fortran_env, only: INPUT_UNIT, OUTPUT_UNIT
  implicit none
  private
  save

  integer,public::nsysi=INPUT_UNIT
  integer,public::nsyso=7
  integer,public::nsyse=OUTPUT_UNIT
end module mainio


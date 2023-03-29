module mainio
   ! Fortran unit numbers for system I/O.  Users should only revert
   ! to the legacy values if these iso_fortran_env variable values
   ! cause conflicts with njoy scratch or user logical units.

   use iso_fortran_env,only: INPUT_UNIT, ERROR_UNIT, OUTPUT_UNIT
   implicit none
   private
   save

   integer,public::nsysi=INPUT_UNIT   !standard input,  legacy value=5
   integer,public::nsyso=ERROR_UNIT   !output file,     legacy value=7
   integer,public::nsyse=OUTPUT_UNIT  !terminal output, legacy value=6
end module mainio


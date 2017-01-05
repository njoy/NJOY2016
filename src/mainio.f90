module mainio
   ! Fortran unit numbers for i/o.  we use the iso_fortran_env
   ! definitions to maximize portability among platforms.
   !
   use iso_fortran_env,only: INPUT_UNIT, ERROR_UNIT, OUTPUT_UNIT
   implicit none
   private
   save

   integer,public::nsysi=INPUT_UNIT     !read from standard input
   integer,public::nsyso=ERROR_UNIT     !write to an output file
   integer,public::nsyse=OUTPUT_UNIT    !write to the terminal

   integer,public::nsysi=INPUT_UNIT   !standard input,  legacy value=5
   integer,public::nsyso=ERROR_UNIT   !output file,     legacy value=7
   integer,public::nsyse=OUTPUT_UNIT  !terminal output, legacy value=6
end module mainio


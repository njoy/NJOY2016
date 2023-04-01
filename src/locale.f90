module locale
   ! Provides localization parameters for NJOY2016.  Each site
   ! should change "lab" and "mx" for their location.  The
   ! "kind" parameter kr is for normal internal real numbers in NJOY.
   ! The default kind is assumed for normal internal integers (either
   ! 4 or 8 bytes should be OK for integers).  The CCCC files use
   ! 4-byte reals, 4-byte integers, and 8-byte Hollerith values,
   ! and k4 and k8 are provided for them.  These values are also used
   ! in subroutine timer when a C-type call is needed.  Caution: the
   ! settings for k4 and k8 may not be portable--some systems might
   ! use 1 and 2 instead of 4 and 8 for these two kinds.
   implicit none
   private
   character(8),public::lab='snl RES'
   character(8),public::mx='        '
   integer,parameter,public::kr=selected_real_kind(12,300)
   integer,parameter,public::k4=4
   integer,parameter,public::k8=8
end module locale


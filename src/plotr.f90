module plotm
   ! provides subroutine plotr for NJOY2016
   use locale
   implicit none
   private
   public plotr

contains

   subroutine plotr
   !--------------------------------------------------------------------
   !
   ! plot cross sections
   !
   ! Handles ENDF data, PENDF or GENDF data at specified temper-
   ! atures, or experimental input data.  Several plots can be
   ! given on each set of axes, with both left and right scales.
   ! Also, several graphs can be given on each page or display.
   ! Error bars may be included for input data.  Flexible titles
   ! and legend blocks are allowed.  All standard combinations of
   ! log and linear axes are supported, either grids or tick marks
   ! can be requested, and scales can be chosen automatically
   ! or set by the user.  In some cases, the x axis is thinned.
   ! In other cases, extra points are added so that, for example,
   ! linear-linear data plots correctly on a log-log graph.
   ! A limited capability for 3-d plots of angle and energy
   ! is included, and the ENDF-6 File 6 format is supported.
   ! Percent difference and ratio plots can be requested.
   !
   ! Plotr writes plot commands on an output file for later use
   ! by the viewr module or an external graphics program.
   !
   !---input--------------------------------------------------------
   !
   !  card 0
   !     nplt          unit for output plot commands
   !     nplt0         unit for input plot commands
   !                     default=0=none
   !                     output plot commands are appended
   !                     to the input plot commands, if any.
   !  card 1
   !     lori          page orientation (def=1)
   !                    0  portrait (7.5x10in)
   !                    1  landscape (10x7.5in)
   !     istyle        character style (def=2)
   !                     1 = roman
   !                     2 = swiss
   !     size          character size option
   !                     pos = height in page units
   !                     neg = height as fraction of subplot size
   !                     (default=0.30)
   !     ipcol         page color (def=white)
   !                    0=white
   !                    1=navajo white
   !                    2=blanched almond
   !                    3=antique white
   !                    4=very pale yellow
   !                    5=very pale rose
   !                    6=very pale green
   !                    7=very pale blue
   !
   ! -----repeat cards 2 through 13 for each curve-----
   !
   !  card 2
   !     iplot         plot index
   !                     99 = terminate plotting job
   !                      1 = new axes, new page
   !                     -1 = new axes, existing page
   !                      n = nth additional plot on existing axes
   !                     -n = start a new set of curves using
   !                          the alternate y axis
   !                     default = 1
   !     iwcol         window color (def=white)
   !                    color list same as for ipcol above
   !     factx         factor for energies (default=1.)
   !     facty         factor for cross-sections (default=1.)
   !     xll,yll       lower-left corner of plot area
   !     ww,wh,wr      window width, height, and rotation angle
   !                   (plot area defaults to one plot per page)
   !
   ! -----cards 3 thru 7 for iplot = 1 or -1 only-----
   !
   !  card 3
   !     t1            first line of title
   !                   60 characters allowed.
   !                   default=none
   !
   !  card 3a
   !     t2            second line of title
   !                   60 characters allowed.
   !                   default=none
   !
   !  card 4
   !     itype         type for primary axes
   !                     1 = linear x - linear y
   !                     2 = linear x - log y
   !                     3 = log x - linear y
   !                     4 = log x - log y
   !                     set negative for 3d axes
   !                     default=4
   !     jtype         type for alternate y axis or z axis
   !                     0 = none
   !                     1 = linear
   !                     2 = log
   !                     default=0
   !     igrid         grid and tic mark control
   !                     0 = no grid lines or tic marks
   !                     1 = grid lines
   !                     2 = tic marks on outside
   !                     3 = tic marks on inside
   !                     default=2
   !     ileg          option to write a legend.
   !                     0 = none
   !                     1 = write a legend block with upper left
   !                         corner at xtag,ytag (see below)
   !                     2 = use tag labels on each curve with
   !                         a vector from the tag to the curve
   !                     default=0
   !     xtag          x coordinate of upper left corner
   !                   of legend block
   !     ytag          y coord of upper left corner
   !                   default=upper left corner of plot
   !
   !  card 5
   !     el            lowest energy to be plotted
   !     eh            highest energy to be plotted
   !     xstep         x axis step
   !                   default = automatic scales
   !                   (default all 3, or none)
   !                   (the actual value of xstep is
   !                       ignored for log scales)
   !
   !  card 5a
   !     xlabl         label for x axis
   !                   60 characters allowed.
   !                    default="energy (ev)"
   !
   !  card 6
   !     yl            lowest value of y axis.
   !     yh            highest value of y axis.
   !     ystep         step for y ayis (linear scales only)
   !                   default = automatic scales
   !                   (default all 3, or none)
   !                   (the actual value of ystep is
   !                       ignored for log scales)
   !
   !  card 6a
   !     ylabl         label for y axis
   !                   60 characters allowed.
   !                    default="cross section (barns)"
   !
   !  card 7   (jtype.gt.0 only)
   !     rbot          lowest value of secondary y axis or z axis
   !     rtop          highest value of secondary y axis or z axis
   !     rstep         step for secondary y axis or z axis
   !                   default for last three = automatic
   !
   !  card 7a  (jtype.gt.0 only)
   !     rl            label for alternate y axis or z axis
   !                   60 characters allowed.
   !                   default=blank
   !
   !  -----cards 8 thru 9 are always given-----
   !
   !  card 8
   !     iverf         version of endf tape
   !                     set to zero for data on input file
   !                       and ignore rest of parameters on card
   !                     set to 1 for gendf data
   !     nin           input tape
   !                   can change for every curve if desired.
   !     matd          desired material
   !     mfd           desired file
   !     mtd           desired section
   !                    mtd=0 means loop over all reactions in mfd
   !                    (usually one page per mt, but for mf=3,
   !                    resonance reactions may have several pages)
   !     temper        temperature for endf data (degK, default=0.)
   !     nth,ntp,nkh   see below (defaults=1)
   !
   !       special meanings for nth,ntp,nkh for file 3 or 5 data
   !            nth   number of subsection to plot
   !              (works for isomer prod, delayed n, etc.)
   !            ntp   special features
   !                  1 for regular plots (default)
   !                  2 for percent difference plots
   !                     read a second "card 8" for percent diff
   !                     of second curve with respect to first
   !                  3 for ratio plots
   !                     read a second "card 8" for ratio
   !                     of second curve to first
   !            nkh   not used
   !
   !       special meanings for nth for file 4 Legendre data
   !            nth   index for Legendre coefficient (p1, p2, ...)
   !
   !       special meanings for nth,ntp,nkh for file 6 data
   !            nth   index for  incident energy
   !            ntp   number of dep. variable in cyle to plot
   !                  (or angle number for law 7)
   !            nkh   number of outgoing particle to plot
   !
   !       special meanings for nth,ntp,nkh for gendf mf=3 data
   !            nth=0 for flux per unit lethargy
   !            nth=1 for cross section (default)
   !            ntp=1 for infinite dilution (default)
   !            ntp=2 for next lowest sigma-zero values, etc.
   !                  set ntp negative for self-shielding factor
   !            ntp=21,22,...
   !                  like 1,2,... except read another source
   !                  and compute percent difference
   !            ntp=41,42,...
   !                  like 1,2,... except read another source
   !                  and compute ratio
   !            ntp=-21,-22,...
   !                  like -1,-2,... except read another source
   !                  and compute percent difference
   !            ntp=-41,-42,...
   !                  like -1,-2,... except read another source
   !                  and compute ratio
   !            nkh=1 for p0 weighting (default)
   !            nkh=2 for p1 weighting (total only)
   !
   !       special meaning for nth for gendf mf=6 data
   !            nth=1 plot 2-d spectrum for group 1
   !            nth=2 plot 2-d spectrum for group 2
   !              etc.
   !        no special flags are needed for mf=6 3d plots
   !
   !      special meanings for nth and ntp for mf7 plots
   !           nth is index for indep. variable (alpha or beta)
   !           ntp=1 selects alpha as indep. variable (default)
   !           ntp=2 selects beta as indep. variable
   !           nkh=1 selects normal s(alpha,beta)
   !           nkh=2 selects script s(alpha,-beta)
   !           nkh=3 selects script s(alpha,beta)
   !
   !      additional guidance for ratio (ntp=3 on card 8) plotting of
   !      mf=3 endf or pendf data:
   !       - if both cross sections are zero, the ratio is defined to be
   !         unity.
   !       - if ntp=3 and mtd=0, then mtd2 is set to zero regardless of
   !         its value on card 8a.
   !       - if mtd=0 and y-axis limits (card 6) were specified, then
   !         those limits will apply to all plots.  For selected mt's
   !         both linear-linear and log-linear plots will be produced.
   !       - if mtd=0 and y-axis limits (card 6) are not given, then
   !         y-axis limits are determined internally.  For selected mt's
   !         both linear-linear and log-log plots will be produced.
   !
   ! -----cards 9 and 10 for 2d plots only-----
   !
   !  card 9
   !     icon          symbol and connection option
   !                     0 = points connected, no symbols
   !                    -i = points not connected, symbol at every
   !                         ith point
   !                     i = points connected, symbol at every ith
   !                         points
   !                     default=0
   !     isym          no. of symbol to be used
   !                     0 = square
   !                     1 = octagon
   !                     2 = triangle
   !                     3 = cross
   !                     4 = ex
   !                     5 = diamond
   !                     6 = inverted triangle
   !                     7 = exed square
   !                     8 = crossed ex
   !                     9 = crossed diamond
   !                     10 = crossed octagon
   !                     11 = double triangle
   !                     12 = crossed square
   !                     13 = exed octagon
   !                     14 = triangle and square
   !                     15 = filled circle
   !                     16 = open circle
   !                     17 = open square
   !                     18 = filled square
   !                     default=0
   !     idash         type of line to plot
   !                     0 = solid
   !                     1 = dashed
   !                     2 = chain dash
   !                     3 = chain dot
   !                     4 = dot
   !                     default=0
   !     iccol         curve color (def=black)
   !                     0=black
   !                     1=red
   !                     2=green
   !                     3=blue
   !                     4=magenta
   !                     5=cyan
   !                     6=brown
   !                     7=purple
   !     ithick        thickness of curve (def=1)
   !                     0 = invisible (for shaded areas)
   !     ishade        shade pattern
   !                     0 = none
   !                     1 to 10 = 10% to 100% gray
   !                     11 to 20 = 45 deg right hatching
   !                     21 to 30 = 45 deg left hatching
   !                     31 to 40 = 45 deg cross hatching
   !                     41 to 50 = shades of green
   !                     51 to 60 = shades of red
   !                     61 to 70 = shades of brown
   !                     71 to 80 = shades of blue
   !                     default=0
   !
   !  card 10  ---ileg.ne.0 only---
   !     aleg          title for curve tag or legend block
   !                   60 characters allowed.
   !                   default=blank
   !
   !  card 10a  ---ileg.eq.2 only---
   !     xtag          x position of tag title
   !     ytag          y position of tag title
   !     xpoint        x coordinate of vector point
   !                    (.le.0 to omit vector)
   !
   ! -----card 11 for 3d plots only-----
   !
   !  card 11
   !     xv,yv,zv      abs. coords of view point
   !                   defaults=15.,-15.,15.
   !     x3,y3,z3      abs. sides of work box volume
   !                   defaults=2.5,6.5,2.5
   !
   !         set x3 or y3 negative to flip the order of the
   !         axis on that side of the work box.
   !
   !  -----cards 12 thru 13 for iverf = 0 only-----
   !
   !  card 12
   !     nform          format code for input data
   !                    0 = free format input with
   !                    optional x and y error bars
   !
   !  card 13   ---nform = 0 only---
   !     xdata          dependent value
   !                    terminate with empty card (/)
   !     ydata          independent value
   !     yerr1          lower y error limit
   !                    no y error bar if zero
   !     yerr2          upper y error limit
   !                    if zero, equals yerr1
   !     xerr1          x left error limit
   !                    no x error bar if zero
   !     xerr2          x right error limit
   !                    if zero, equals xerr1
   !
   !
   ! all curves contain at least 10 points per decade (see delta).
   ! code can plot curves containing fewer than 10000 points (see
   ! max) without thinning.  curves with more points are thinned
   ! based on a minimum spacing determined from max and the
   ! length of the x axis.
   !
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use endf ! provides endf routines and variables
   use util ! provides times,openz,error,closz,repoz,mess
   ! internals
   integer::nplt,nplt0,ninl,lori,istyle,ipcol,i3d,iplot
   integer::iwcol,ierrb,n1,i,n2,itype,jtype,igrid,ileg
   integer::nz,nx,ny,nr,matd,mfd,mtd,nth,ntp,nkh
   integer::iauto,ipass,icon,isym,idash,iccol,ithick,ishade
   integer::nb,nw,idone,ik,mmf,lf,lep,net,iet,nmu,imu
   integer::ith,lat,lasym,lln,nbeta,int2,ibeta,lt,l,l2,nt,np
   integer::ll,ithere,nin,nleg,n,nnn,idis,ir,ip
   integer::jnoth,iii,jj,nlast,major,minor,nsk,law,ne,ie
   integer::ltt,loc,nsigz,ntw,ngn,ngg,jbase,locngn,locngg
   integer::nl,ig1,ig,ng2,ig2lo,nwa,jnow,j,nform,iyauto,mthx
   integer::nin2,matd2,mtd2,nth2,ntp2,nkh2,iverf2,mfd2,loc2,lttt
   real(kr)::xpaper,ypaper,xmarg,ymarg,dsize,dxv,dyv,dzv,dx3,dy3,dz3
   real(kr)::time,size,xpage,ypage,factx,facty,xll,yll,ww,wh,wr
   real(kr)::xtag,ytag,xlf,xstep,xleft,xright,yb,yt,ystep,yystep
   real(kr)::ybot,ytop,rbot,rtop,rstep,xpoint
   real(kr)::xv,yv,zv,x3,y3,z3,elt,eht,tem,tnow,test,beta,s
   real(kr)::xh,temper,enow,yf,enext,ee1,ee2,thin,elast,etmax
   real(kr)::st1,st2,st3,st4,stmax,stmin,ststp,temp,dleth,dener
   real(kr)::sigig,flag,factx1,facty1,temper2,reset,enxt,zz,z2,itypx
   real(kr)::z(15)
   integer,parameter::mmax=20000   !same in plotr and viewr
   integer,parameter::nwamax=45000
   integer,parameter::maxaa=200000
   real(kr),dimension(nwamax)::a
   real(kr),dimension(maxaa)::aa
   real(kr),dimension(mmax)::x,y,b,dxm,dxp,dym,dyp
   equivalence (x(1),aa(1)),(y(1),aa(mmax+1)),(b(1),aa(2*mmax+1))
   equivalence (dxm(1),aa(3*mmax+1)),(dxp(1),aa(4*mmax+1))
   equivalence (dym(1),aa(5*mmax+1)),(dyp(1),aa(6*mmax+1))
   character(60)::t1,t2,xl,yl,rl
   character(60)::aleg
   character(60)::strng
   character(10)::name
   character(60)::text
   character(14)::xlabld='<e>nergy (e<v)'
   character(23)::ylabld='<c>ross section (barns)'
   character(1)::qu=''''
   real(kr),parameter::eps=1.e-7_kr
   real(kr),parameter::delta=1.26e0_kr
   real(kr),parameter::etop=1.e10_kr
   real(kr),parameter::tiny=1.e-33_kr
   real(kr),parameter::elow=1.e-5_kr
   real(kr),parameter::therm=296.e0_kr
   real(kr),parameter::sabmin=1.e-32_kr
   real(kr),parameter::big=9.e9_kr
   real(kr),parameter::small=1.e-12_kr
   real(kr),parameter::zero=0
   real(kr),parameter::ten=10

   ! default paper size is US letter size.
   ! default character size is 0.30 in.
   ! see similar statements in viewr.
   xpaper=8.5e0_kr
   ypaper=11.0e0_kr
   xmarg=1.0e0_kr
   ymarg=1.0e0_kr
   dsize=.30e0_kr

   ! default viewpoint and workbox
   dxv=15.e0_kr
   dyv=-15.e0_kr
   dzv=15.e0_kr
   dx3=2.5e0_kr
   dy3=6.5e0_kr
   dz3=2.5e0_kr

   !--initialize plotting
   read(nsysi,*) nplt,nplt0
   call openz(nplt,1)
   call timer(time)
   write(nsyso,&
     '(/'' plotr...plot endf, pendf, gendf, or input data'',&
     &22x,f8.1,''s'')') time
   write(nsyse,'(/'' plotr...'',60x,f8.1,''s'')') time
   ninl=0
   lori=1
   istyle=2
   size=dsize
   ipcol=0
   read(nsysi,*) lori,istyle,size,ipcol
   write(nsyso,'(/&
     &'' lori ................................. '',i10/&
     &'' istyle ............................... '',i10/&
     &'' size ................................. '',f10.3/&
     &'' ipcol ................................ '',i10)')&
     lori,istyle,size,ipcol
   write(nplt,'(2i2,f8.3,i8,'' /'')') lori, istyle,size,ipcol
   !  default page size is paper size with 0.5in margins all around
   if (lori.eq.0) then
      xpage=xpaper-xmarg
      ypage=ypaper-ymarg
   else
      xpage=ypaper-xmarg
      ypage=xpaper-ymarg
   endif

   !--loop over plots
   i3d=0
  110 continue
   iplot=1
   iwcol=0
   factx=1
   facty=1
   xll=0
   yll=0
   ww=xpage
   wh=ypage
   wr=0
   read(nsysi,*) iplot,iwcol,factx,facty,xll,yll,ww,wh,wr
   if (iplot.eq.99) go to 700
   write(nsyso,'(/&
     &'' -------------------------------------------------''//&
     &'' iplot ................................ '',i10/&
     &'' iwcol ................................ '',i10/&
     &'' factx ................................ '',1p,e10.2/&
     &'' facty ................................ '',e10.2,0p/&
     &'' xll .................................. '',f10.3/&
     &'' yll .................................. '',f10.3/&
     &'' ww ................................... '',f10.3/&
     &'' wh ................................... '',f10.3/&
     &'' wr ................................... '',f10.3)')&
     iplot,iwcol,factx,facty,xll,yll,ww,wh,wr
   ierrb=0

   !--if this is first plot on these axes
   !--read in title lines
   if (iabs(iplot).eq.1) then
      text=' '
      read(nsysi,*) text
      n1=0
      do i=1,60
         if (text(i:i).ne.' ') n1=i
      enddo
      t1=' '
      if (n1.gt.0) t1=text(1:n1)
      text=' '
      read(nsysi,*) text
      n2=0
      do i=1,60
         if (text(i:i).ne.' ') n2=i
      enddo
      t2=' '
      if (n2.gt.0) t2=text(1:n2)
      write(nsyso,'(/12x,a/12x,a)') t1,t2

      !--set up plot type and grids
      itype=4
      jtype=0
      igrid=2
      ileg=0
      xtag=0
      ytag=0
      read(nsysi,*) itype,jtype,igrid,ileg,xtag,ytag
      write(nsyso,'(/&
        &'' itype ................................ '',i10/&
        &'' jtype ................................ '',i10/&
        &'' igrid ................................ '',i10/&
        &'' ileg ................................. '',i10)')&
        itype,jtype,igrid,ileg
      if (ileg.eq.1) write(nsyso,'(/&
        &'' xtag ................................. '',1p,e10.2/&
        &'' ytag ................................. '',1p,e10.2)')&
        xtag,ytag
      i3d=0
      if (itype.le.0) then
         i3d=1
         itype=-itype
      endif

      !--read in x-axis limits and label
      xlf=0
      xh=0
      xstep=0
      read(nsysi,*) xlf,xh,xstep
      if (xlf.eq.zero.and.xh.eq.zero) then
         nz=0
      else
         nz=2
         if (xstep.ne.zero) nz=3
      endif
      if (itype.gt.2.and.nz.ge.2) xstep=1
      if (nz.gt.0.and.xstep.eq.zero)&
        call error('plotr','error in axis input',' ')
      xleft=xlf
      xright=xh
      text=' '
      read(nsysi,*) text
      nx=0
      do i=1,60
         if (text(i:i).ne.' ') nx=i
      enddo
      xl=' '
      if (nx.gt.0) xl=text(1:nx)
      write(nsyso,'(/&
        &'' xlo .................................. '',1p,e10.2/&
        &'' xhi .................................. '',1p,e10.2/&
        &'' xstep ................................ '',1p,e10.2/&
        &'' xlabel ... '',a)')&
        xleft,xright,xstep,xl

      !--read in y-axis limits and label
      yb=0
      yt=0
      yystep=0
      iyauto=0
      read(nsysi,*) yb,yt,yystep
      ystep=yystep
      if (yb.eq.zero.and.yt.eq.zero) then
         nz=0
      else
         nz=2
         if (ystep.ne.zero) nz=3
      endif
      if ((itype.eq.2.or.itype.eq.4).and.nz.ge.2) ystep=1
      if (nz.gt.0.and.ystep.eq.zero)&
        call error('plotr','error in axis input',' ')
      ybot=yb
      ytop=yt
      text=' '
      read(nsysi,*) text
      ny=0
      do i=1,60
         if (text(i:i).ne.' ') ny=i
      enddo
      yl=' '
      if (ny.gt.0) yl=text(1:ny)
      write(nsyso,'(/&
        &'' ylo .................................. '',1p,e10.2/&
        &'' yhi .................................. '',1p,e10.2/&
        &'' ystep ................................ '',1p,e10.2/&
        &'' ylabel ... '',a)')&
        ybot,ytop,ystep,yl

      !--read in data for alternate y axis or z axis
      if (jtype.ne.0) then
         rbot=0
         rtop=0
         rstep=0
         read(nsysi,*) rbot,rtop,rstep
         if (rbot.eq.zero.and.rtop.eq.zero) then
            nz=0
         else
            nz=2
            if (rstep.ne.zero) nz=3
         endif
         if (rtop.lt.rbot)&
           call error('plotr','error in axis input',' ')
         if (nz.gt.2.and.rstep.eq.zero)&
           call error('plotr','error in axis input',' ')
         text=' '
         read(nsysi,*) text
         nr=0
         do i=1,60
            if (text(i:i).ne.' ') nr=i
         enddo
         rl=' '
         if (nr.gt.0) rl=text(1:nr)
         write(nsyso,'(/&
           &'' rlo .................................. '',1p,e10.2/&
           &'' rhi .................................. '',1p,e10.2/&
           &'' rstep ................................ '',1p,e10.2/&
           &'' rlabel ... '',a)')&
           rbot,rtop,rstep,rl
      endif
   endif

   !--read data source for next curve
   iverf=0
   nin=0
   matd=0
   mfd=0
   mtd=0
   temper=0
   nth=1
   ntp=1
   nkh=1
   nin2=0
   read(nsysi,*) iverf,nin,matd,mfd,mtd,temper,nth,ntp,nkh
   if (iverf.gt.0) then
      if (nin.ne.ninl) call closz(ninl)
      if (nin.ne.ninl) call openz(nin,0)
      ninl=nin
      write(nsyso,'(/&
        &'' iverf ................................ '',i10/&
        &'' nin .................................. '',i10/&
        &'' matd ................................. '',i10/&
        &'' mfd .................................. '',i10/&
        &'' mtd .................................. '',i10/&
        &'' temp ................................. '',1p,e10.2/&
        &'' nth .................................. '',i10/&
        &'' ntp .................................. '',i10/&
        &'' nkh .................................. '',i10)')&
        iverf,nin,matd,mfd,mtd,temper,nth,ntp,nkh
      iauto=0
      if (mtd.eq.0) iauto=1
      ipass=0
      if ((iverf.gt.1.and.mfd.eq.3.and.ntp.gt.1).or.&
        (iverf.eq.1.and.abs(ntp).gt.20)) then
         matd2=0
         mfd2=0
         mtd2=0
         temper2=0
         nth2=1
         ntp2=1
         nkh2=1
         read(nsysi,*)&
           iverf2,nin2,matd2,mfd2,mtd2,temper2,nth2,ntp2,nkh2
         call openz(nin2,0)
         if (ntp.eq.3.and.mtd.eq.0) mtd2=mtd
         write(nsyso,'(/&
           &'' iverf2 ............................... '',i10/&
           &'' nin2 ................................. '',i10/&
           &'' matd2 ................................ '',i10/&
           &'' mfd2 ................................. '',i10/&
           &'' mtd2 ................................. '',i10/&
           &'' temp2 ................................ '',1p,e10.2/&
           &'' nth2 ................................. '',i10/&
           &'' ntp2 ................................. '',i10/&
           &'' nkh2 ................................. '',i10)')&
           iverf2,nin2,matd2,mfd2,mtd2,temper2,nth2,ntp2,nkh2
      endif
   endif

   !--read plotting parameters for next 2d curve
   if (i3d.eq.0) then
      icon=0
      isym=0
      idash=0
      iccol=0
      ithick=1
      ishade=0
      read(nsysi,*) icon,isym,idash,iccol,ithick,ishade
      write(nsyso,'(/&
        &'' icon ................................. '',i10/&
        &'' isym ................................. '',i10/&
        &'' idash ................................ '',i10/&
        &'' iccol ................................ '',i10/&
        &'' ithick ............................... '',i10/&
        &'' ishade ............................... '',i10)')&
        icon,isym,idash,iccol,ithick,ishade

      !--read in legend or tag title lines
      if (ileg.ne.0) then
         text=' '
         read(nsysi,*) text
         nleg=0
         do i=1,60
            if (text(i:i).ne.' ') nleg=i
         enddo
         aleg=' '
         if (nleg.gt.0) aleg=text(1:nleg)
         if (ileg.eq.1) write(nsyso,'(/'' legend ... '',a)') aleg
         if (ileg.eq.2) then
            xtag=0
            ytag=0
            xpoint=0
            read(nsysi,*) xtag,ytag,xpoint
            write(nsyso,'(/&
              &'' tag ...... '',a/&
              &'' xtag ................................. '',1p,e10.2/&
              &'' ytag ................................. '',1p,e10.2/&
              &'' xpoint ............................... '',1p,e10.2&
              &)') aleg,xtag,ytag,xpoint
         endif
      endif

   !--read parameters for plotting 3d surface
   else
      xv=dxv
      yv=dyv
      zv=dzv
      x3=dx3
      y3=dy3
      z3=dz3
      read(nsysi,*) xv,yv,zv,x3,y3,z3
      write(nsyso,'(/&
        &'' 3d viewpoint ......................... '',f10.3/&
        &''                                        '',f10.3/&
        &''                                        '',f10.3/&
        &'' 3d workbox ........................... '',f10.3/&
        &''                                        '',f10.3/&
        &''                                        '',f10.3)')&
        xv,yv,zv,x3,y3,z3
   endif

   !--branch to appropriate data type
   if (iverf.eq.0) go to 500
   if (iverf.eq.1) go to 400

   !--endf or pendf data
   !--locate desired tabulation
   call repoz(nin)
   call tpidio(nin,0,0,a,nb,nw)
   if (nin2.ne.0) then
       call repoz(nin2)
       call tpidio(nin2,0,0,a,nb,nw)
   endif
   elt=xleft/factx
   eht=xright/factx
   idone=0
   do while (idone.eq.0)
      call contio(nin,0,0,a,nb,nw)
      if (math.lt.0) then
         write(strng,&
           '(''1desired mat and temp not found '',i4,f10.1)')&
            matd,temper
         call error('plotr',strng,' ')
      endif
      if (math.eq.matd) then
         if (mfd.eq.7) then
            idone=1
         else
            if (iverf.ge.5) call contio(nin,0,0,a,nb,nw)
            if (iverf.ge.6) call contio(nin,0,0,a,nb,nw)
            call contio(nin,0,0,a,nb,nw)
            tem=a(1)
            if (abs(tem-temper).le.temper/1000) idone=1
         endif
      endif
      if (idone.eq.0) call tomend(nin,0,0,a)
   enddo
   if (nin2.eq.0) go to 320
   idone=0
   do while (idone.eq.0)
      call contio(nin2,0,0,a,nb,nw)
      if (math.lt.0) then
         write(strng,&
           '(''desired mat2 and temp2 not found '',i4,f10.1)')&
            matd2,temper2
         call error('plotr',strng,' ')
      endif
      if (math.eq.matd2) then
         if (mfd2.eq.7) then
            idone=1
         else
            if (iverf2.ge.5) call contio(nin2,0,0,a,nb,nw)
            if (iverf2.ge.6) call contio(nin2,0,0,a,nb,nw)
            call contio(nin2,0,0,a,nb,nw)
            tem=a(1)
            if (abs(tem-temper2).le.temper2/1000) idone=1
         endif
      endif
      if (idone.eq.0) call tomend(nin2,0,0,a)
   enddo

   !--auto reaction loop goes through here.
  320 continue
   call findf(matd,mfd,mtd,nin)
   if (nin2.ne.0) call findf(matd2,mfd2,mtd2,nin2)
  321 continue
   call contio(nin,0,0,a,nb,nw)
   mthx=mth
   if (nin2.ne.0) then
      call contio(nin2,0,0,a,nb,nw)
      mtd2=mth
      do while (mthx.lt.mtd2.and.mfh.gt.0)
         call tosend(nin,0,0,a)
         call contio(nin,0,0,a,nb,nw)
         mthx=mth
      enddo
      do while (mtd2.lt.mthx.and.mfh.gt.0)
         call tosend(nin2,0,0,a)
         call contio(nin2,0,0,a,nb,nw)
         mtd2=mth
      enddo
   endif
   if (mfh.eq.0.and.iauto.eq.1) go to 110
   mtd=mthx
   if (nin2.ne.0) mtd2=mth
   if (iauto.gt.0.and.ipass.eq.0) eht=0
   if (iauto.gt.0) then
      elt=0
      xleft=0
      xright=0
      xstep=0
      ybot=0
      ytop=0
      ystep=0
      iyauto=0
      if (ntp.eq.3.and.yb.ne.zero) ybot=yb
      if (ntp.eq.3.and.yt.ne.zero) ytop=yt
      if (ntp.eq.3.and.yystep.ne.zero) ystep=yystep
      if (ybot.ne.0.and.ytop.ne.0) iyauto=1
   endif
   ik=1
   mmf=3
   if (i3d.eq.1) go to 2300

   !--endf or pendf 2d plots
   if (mfd.eq.3) go to 1409
   if (mfd.eq.4) then
      mmf=mfd
      go to 1409
   endif
   if (mfd.eq.7) go to 1410

   !--curves from file5, file15, or file6
   if (mfd.eq.5.or.mfd.eq.15.or.mfd.eq.6) then
      lttt=l2h
      mmf=mfd
      idone=0
      do while (idone.eq.0)
         !if (mtd.lt.221.or.mtd.ge.250) then
         if (lttt.ne.5) then
            call tab1io(nin,0,0,a,nb,nw)
            lf=l2h
            if (mfd.ne.6.and.lf.ne.1)&
              call error('plotr','lf=1 only for mf5 or mf15.',' ')
            if (mfd.eq.6.and.lf.ne.1.and.lf.ne.7)&
              call error('plotr','lf=1 or 7 only for file6.',' ')
            if (lf.eq.7) mmf=1
         endif
         call tab2io(nin,0,0,a,nb,nw)
         lep=l2h
         !if (mtd.ge.221.and.mtd.lt.250) lep=2
         if (lttt.eq.5) lep=2
         if (ik.eq.nkh) then
            idone=1
         else
            net=n2h
            iet=1
            do while (iet.lt.net)
               if (lf.ne.7) then
                   if (mfd.ne.6) call tab1io(nin,0,0,a,nb,nw)
                  if (mfd.eq.6) call listio(nin,0,0,a,nb,nw)
                  do while (nb.ne.0)
                     call moreio(nin,0,0,a,nb,nw)
                  enddo
               else
                  call tab2io(nin,0,0,a,nb,nw)
                  nmu=n2h
                  do imu=1,nmu
                     call tab1io(nin,0,0,a,nb,nw)
                     do while (nb.ne.0)
                        call moreio(nin,0,0,a,nb,nw)
                     enddo
                  enddo
               endif
               iet=iet+1
            enddo
            ik=ik+1
         endif
      enddo
   endif
   ith=1
   do while (ith.ne.nth)
      if (lf.ne.7) then
         if (mfd.ne.6) call tab1io(nin,0,0,a,nb,nw)
         if (mfd.eq.6) call listio(nin,0,0,a,nb,nw)
         do while (nb.ne.0)
            call moreio(nin,0,0,a,nb,nw)
         enddo
      else
         call tab2io(nin,0,0,a,nb,nw)
         nmu=n2h
         do imu=1,nmu
            call tab1io(nin,0,0,a,nb,nw)
            do while (nb.ne.0)
               call moreio(nin,0,0,a,nb,nw)
            enddo
         enddo
      endif
      ith=ith+1
   enddo
   if (lf.eq.7) then
      call tab2io(nin,0,0,a,nb,nw)
      nmu=n2h
      imu=1
      idone=0
      do while (idone.eq.0)
         if (imu.eq.ntp) then
            idone=1
         else
            call tab1io(nin,0,0,a,nb,nw)
            do while (nb.ne.0)
               call moreio(nin,0,0,a,nb,nw)
            enddo
            imu=imu+1
         endif
      enddo
   endif
   go to 1409

   !--curves from file7
 1410 continue
   mmf=7
   lat=l2h
   lasym=n1h
   call listio(nin,0,0,a,nb,nw)
   lln=l1h
   call tab2io(nin,0,0,a,nb,nw)
   nbeta=n2h
   int2=nint(a(8))
   if (ntp.gt.1) go to 1440
   ibeta=1
   idone=0
   do while (idone.eq.0)
      if (ibeta.eq.nth) then
         idone=1
      else
         call tab1io(nin,0,0,a,nb,nw)
         lt=l1h
         do while (nb.ne.0)
            call moreio(nin,0,0,a,nb,nw)
         enddo
         if (lt.ne.0) then
            do i=1,lt
               call listio(nin,0,0,a,nb,nw)
               do while (nb.ne.0)
                  call moreio(nin,0,0,a,nb,nw)
               enddo
            enddo
         endif
         ibeta=ibeta+1
         if (ibeta.gt.nbeta) call error('plotr',&
           'illegal nth for mf7',' ')
      endif
   enddo
   l=1
   l2=1
   call tab1io(nin,0,0,a(l),nb,nw)
   tnow=c1h
   lt=l1h
   nt=0
   nr=n1h
   np=n2h
   l=l+nw
   do while (nb.ne.0)
      call moreio(nin,0,0,a(l),nb,nw)
   enddo
   idone=0
   do while (idone.eq.0)
      test=1
      test=test/10+tnow/100
      if (abs(temper-tnow).lt.test) then
         idone=1
      else
         ll=l
         call listio(nin,0,0,a(l),nb,nw)
         tnow=c1h
         np=n1h
         l=l+nw
         do while (nb.ne.0)
            call moreio(nin,0,0,a(l),nb,nw)
            l=l+nw
         enddo
         do i=1,np
            a(6+2*nr+2*i)=a(ll+5+i)
         enddo
         l=ll
         nt=nt+1
         if (nt.gt.lt) call error('plotr',&
           'temperature not found',' ')
      endif
   enddo
   beta=a(2)
   if (lat.eq.1) beta=beta*therm/temper
   do i=1,np
      s=a(6+2*nr+2*i)
      if (lln.eq.1) s=exp(s)
      if (nkh.eq.2) s=s*exp(beta/2)
      if (nkh.eq.3) s=s*exp(-beta/2)
      if (s.lt.sabmin) s=sabmin
      a(6+2*nr+2*i)=s
   enddo
   if (lln.ne.0) then
      do i=1,nr
         a(5+2*i)=a(5+2*i)+1
      enddo
   endif
   go to 1409
 1440 continue
   l2=npage+50
   a(l2)=0
   a(l2+1)=0
   a(l2+2)=0
   a(l2+3)=0
   a(l2+4)=1
   a(l2+5)=nbeta
   a(l2+6)=nbeta
   if (lln.eq.1) int2=int2+1
   a(l2+7)=int2
   do ibeta=1,nbeta
      l=1
      ithere=0
      call tab1io(nin,0,0,a(l),nb,nw)
      tnow=c1h
      lt=l1h
      nt=0
      nr=n1h
      l=l+nw
      do while (nb.ne.0)
         call moreio(nin,0,0,a(l),nb,nw)
         l=l+nw
      enddo
      do while (nt.le.lt)
         if (nt.ne.0) then
            ll=l
            call listio(nin,0,0,a(l),nb,nw)
            tnow=c1h
            np=n1h
            l=l+nw
            do while (nb.ne.0)
               call moreio(nin,0,0,a(l),nb,nw)
               l=l+nw
            enddo
            l=ll
         endif
         if (abs(temper-tnow).lt..01) then
            a(l2+6+2*ibeta)=a(2)
            if (lasym.gt.0) x(ibeta)=a(2)*factx
            if (nt.eq.0) yf=a(6+2*nr+2*nth)
            if (nt.gt.0) yf=a(ll+5+nth)
            if (lln.eq.1) yf=exp(yf)
            beta=a(2)
            if (lat.eq.1) beta=beta*therm/temper
            if (nkh.eq.2) yf=yf*exp(beta/2)
            if (nkh.eq.3) yf=yf*exp(-beta/2)
            if (yf.lt.sabmin) yf=sabmin
            a(l2+7+2*ibeta)=yf
            if (lasym.gt.0) y(ibeta)=yf*facty
            ithere=1
         endif
         nt=nt+1
      enddo
      if (ithere.eq.0)&
        call error('plotr', 'temperature not found',' ')
   enddo
   if (lasym.eq.0) go to 1409
   n=nbeta
   go to 610

 1409 continue
   nnn=0
   enow=0
   if (mmf.eq.4) then
      call gety4(enow,enext,idis,yf,nth,nin,a)
   else if (mmf.eq.6) then
      call gety6(ntp,enow,enext,idis,yf,nin,a,lep,lf,b,1)
      enext=enext-enext/10000
   else if (mmf.eq.7) then
      ir=1
      ip=2
      call terpa(yf,zero,enext,idis,a(l2),ip,ir)
   else if (nin2.ne.0) then
      enext=big
      reset=-1
      loc2=50+npage
      call getz(reset,enxt,idis,zz,0,a)
      call getz(enow,enxt,idis,zz,nin,a)
      if (enxt.lt.enext) enext=enxt
      call getz(enow,enxt,idis,zz,nin2,a(loc2))
      nnn=nint(a(loc2+5))
      if (enxt.lt.enext) enext=enxt
   else
      call gety1(enow,enext,idis,yf,nin,a)
      nnn=nint(a(6))
   endif
   if (enext.lt.elt) enext=elt
   if (enext.eq.zero) enext=elow
   jnoth=0
   if (enext.lt.1.e3) jnoth=1
   ipass=ipass+1
   ee1=0
   ee2=0

   !--if necessary, find last point in frame
   if (eht.gt.zero) go to 335
   iii=0
   idone=0
   do while (idone.eq.0)
      enow=enext
      iii=iii+1
      if (mmf.eq.4) then
         call gety4(enow,enext,idis,yf,nth,nin,a)
      else if (mmf.eq.6) then
         call gety6(ntp,enow,enext,idis,yf,nin,a,lep,lf,b,1)
      else if (mmf.eq.7) then
         ir=1
         ip=2
         call terpa(yf,enow,enext,idis,a(l2),ip,ir)
      else if (nin2.ne.0) then
         call getz(enow,enext,idis,zz,nin2,a(loc2))
         if (jnoth.eq.1.and.nnn.gt.3000) then
            if (iii.eq.300) ee1=enow
            if (iii.eq.nnn-500) ee2=enow
         endif
      else
         call gety1(enow,enext,idis,yf,nin,a)
         if (jnoth.eq.1.and.nnn.gt.3000) then
            if (iii.eq.300) ee1=enow
            if (iii.eq.nnn-500) ee2=enow
         endif
      endif
      eht=enow
      if (enext.ge.etop-etop/1000) idone=1
   enddo
   go to 320

   !--find first point in frame
   !--and set up thinning/thickening logic
  335 continue
   enow=enext
   if (ee1.gt.tiny) then
       jj=int(log10(ee1))
       ee1=ten**jj
   endif
   if (ee2.gt.tiny) then
      jj=int(log10(ee2))
      ee2=ten**jj
   endif
   if (iauto.gt.0) then
      itype=1
      if (jnoth.eq.1.and.ipass.gt.2) itype=4
      if (ipass.gt.2.and.ee1.eq.zero) ipass=5
      if (jnoth.eq.1.and.ee1.gt.zero.and.ipass.eq.4) then
         elt=ee1
         eht=ee1*100
         if (eht.gt.ee2) eht=ee2
         xleft=elt
         xright=eht
         xstep=1
         if (eht.eq.ee2) ipass=5
      endif
      if (jnoth.eq.1.and.ee1.gt.zero.and.ipass.eq.5) then
         elt=ee1*100
         eht=ee2
         xleft=elt
         xright=eht
         xstep=1
      endif
      if (itype.eq.4.and.ntp.eq.3.and.iyauto.ne.0) itype=3
   endif
   if (elt.eq.zero) elt=enext
   if (elt.gt.enow) enow=elt
   if (itype.ne.3.and.itype.ne.4) then
      thin=(eht-elt)/(mmax-500)
   else
      thin=ten**(log10(eht/elt)/(mmax-500))
   endif
   if (mmf.eq.4) then
      call gety4(enow,enext,idis,yf,nth,nin,a)
   else if (mmf.eq.6) then
      call gety6(ntp,enow,enext,idis,yf,nin,a,lep,lf,b,1)
   else if (mmf.eq.7) then
      ir=1
      ip=2
      call terpa(yf,enow,enext,idis,a(l2),ip,ir)
   else if (nin2.ne.0) then
      enext=big
      call getz(enow,enxt,idis,zz,nin,a)
      if (enxt.lt.enext) enext=enxt
      call getz(enow,enxt,idis,z2,nin2,a(loc2))
      if (enxt.lt.enext) enext=enxt
      if (ntp.eq.2) then
         yf=100*(z2-zz)+small
         if (zz.ne.zero) yf=yf/zz
      else if (ntp.eq.3) then
         yf=z2
         if (zz.ne.zero) yf=yf/zz
         if (z2.eq.zero.and.zz.eq.zero) yf=1
      endif
   else
      call gety1(enow,enext,idis,yf,nin,a)
   endif
   itypx=itype
   if (iplot.lt.-1) then
      if (jtype.eq.1) then
         if (itype.eq.2) itypx=1
         if (itype.eq.4) itypx=3
      else
         if (itype.eq.1) itypx=2
         if (itype.eq.3) itypx=4
      endif
   endif
   if (itypx.eq.2.or.itypx.eq.4) then
      test=ten**(-15)
      if (yf.le.test) yf=test
      if (enext.gt.delta*enow) enext=delta*enow
   endif
   n=1
   x(n)=enow*factx
   y(n)=yf*facty
   elast=enow
   nlast=n
   test=enext-enext/10000
   if (idis.ne.0.and.test.gt.enow) enext=test

   !--retrieve rest of points in frame
   !--thin or thicken as necessary
  360 continue
   if (enext.le.eht) go to 370
   if (abs(enow-eht).le.eps) go to 390
   if (enext.gt.eht.and.enow.lt.eht) enext=eht
   if (enext.eq.eht) go to 370
   go to 390
  370 continue
   enow=enext
   if (mmf.eq.4) then
      call gety4(enow,enext,idis,yf,nth,nin,a)
   else if (mmf.eq.6) then
      call gety6(ntp,enow,enext,idis,yf,nin,a,lep,lf,b,1)
   else if (mmf.eq.7) then
      ir=1
      ip=2
      call terpa(yf,enow,enext,idis,a(l2),ip,ir)
   else if (nin2.ne.0) then
      enext=big
      call getz(enow,enxt,idis,zz,nin,a)
      if (enxt.lt.enext) enext=enxt
      call getz(enow,enxt,idis,z2,nin2,a(loc2))
      if (enxt.lt.enext) enext=enxt
      if (ntp.eq.2) then
         yf=100*(z2-zz)+small
         if (zz.ne.zero) yf=yf/zz
      else if (ntp.eq.3) then
         yf=z2
         if (zz.ne.zero) yf=yf/zz
         if (z2.eq.zero.and.zz.eq.zero) yf=1
      endif
   else
      call gety1(enow,enext,idis,yf,nin,a)
   endif
   test=etop-etop/100
   if (enext.gt.test.and.yf.gt.zero) enext=enow+enow/10000
   if (enext.gt.test.and.eht.gt.enext) eht=enow
   test=enow+enow/1000
   if (idis.ne.0.and.enext.gt.test) enext=enext-enext/10000
   if (itypx.eq.3.or.itypx.eq.4) go to 375
   if (idis.gt.0) go to 380
   if (itypx.eq.1) then
      if (yf.eq.zero.and.y(n).gt.zero) go to 380
      if (yf.gt.zero.and.y(n).eq.zero) go to 380
   else
      test=ten**(-15)
      if (yf.le.test) yf=test
      if (yf.eq.test.and.y(n).gt.test) go to 380
      if (yf.gt.test.and.y(n).eq.test) go to 380
   endif
   if (enext.ge.enow+thin) go to 380
   if (enow.lt.elast+thin) go to 370
   go to 380
  375 if (idis.gt.0) go to 380
   if (itypx.eq.3) then
      if (yf.eq.zero.and.y(n).gt.zero) go to 380
      if (yf.gt.zero.and.y(n).eq.zero) go to 380
   else
      test=ten**(-15)
      if (yf.le.test) yf=test
      if (yf.eq.test.and.y(n).gt.test) go to 380
      if (yf.gt.test.and.y(n).eq.test) go to 380
   endif
   if (enext.ge.thin*elast) go to 380
   if (enow.lt.thin*elast) go to 370
  380 continue
   n=n+1
   if (n.gt.mmax) call error('plotr','storage exceeded1',' ')
   x(n)=enow*factx
   y(n)=yf*facty
   if (enext.gt.delta*enow) enext=delta*enow
   if (itypx.ne.3.and.itypx.ne.4) then
      if (enext.lt.enow+thin) enext=enow+thin
   else
      if (delta.lt.thin) enext=thin*enow
   endif
   elast=enow
   if (yf.ne.zero.and.itypx.ne.2.and.itypx.ne.4) nlast=n
   test=ten**(-15)
   if (yf.ne.test.and.itype.ne.1.and.itype.ne.3) nlast=n
   go to 360
  390 continue
   if (nlast.lt.n) n=nlast+1
   call tosend(nin,0,0,a)
   if (nin2.ne.0) call tosend(nin2,0,0,a)

   !--for automatic linear plots of non-threshold reactions,
   !--readjust the vertical scale to make high-energy data show up
   if (iauto.ne.0) then
      if (itype.le.2.and.iyauto.eq.0) then
         test=1
         if (x(1)/factx.le.test) then
            etmax=x(n-1)
            st1=0
            st2=0
            st3=0
            st4=y(n-1)
            do i=1,n
               if (x(i).gt.etmax/5.and.y(i).gt.st1) st1=y(i)
               if (x(i).gt.etmax/3.and.y(i).gt.st2) st2=y(i)
               if (x(i).gt.etmax/2.and.y(i).gt.st3) st3=y(i)
            enddo
            stmax=st1+st1/10
            if (stmax.gt.50*st4) stmax=st2+st2/10
            if (stmax.lt.st3) stmax=st3
            stmin=0
            call ascale(4,stmin,stmax,major,minor)
            ststp=stmax/major
            ybot=0
            ytop=stmax
            ystep=ststp
         endif
      endif
   endif
   go to 610

   !--endf or pendf 3d plots
   !--skip to desired subsection
  2300 continue
   law=0
   if (mfd.eq.4) go to 2410
   if (mfd.eq.5) go to 2510
   if (mfd.eq.7) go to 2610
   if (mfd.eq.6.and.mtd.ge.221.and.mtd.le.250) go to 2520
   nsk=nkh-1
   if (nsk.ne.0) then
      do i=1,nsk
         call tab1io(nin,0,0,a,nb,nw)
         law=l2h
         call tab2io(nin,0,0,a,nb,nw)
         ne=n2h
         do ie=1,ne
            if (law.ne.7) then
               call listio(nin,0,0,a,nb,nw)
               do while (nb.ne.0)
                  call moreio(nin,0,0,a,nb,nw)
               enddo
            else
               call tab2io(nin,0,0,a,nb,nw)
               nmu=n2h
               do imu=1,nmu
                  call tab1io(nin,0,0,a,nb,nw)
                  do while (nb.ne.0)
                     call moreio(nin,0,0,a,nb,nw)
                  enddo
               enddo
            endif
         enddo
      enddo
   endif
   call tab1io(nin,0,0,a,nb,nw)
   law=l2h
   if (law.eq.0) go to 2380
   if (law.eq.1) go to 2520
   if (law.eq.2) go to 2415
   if (law.eq.3) go to 2380
   if (law.eq.4) go to 2380
   if (law.eq.7) go to 2520
   call error('plotr','illegal mf6 law.',' ')

   !--no distribution
  2380 continue
   write(strng,'(''for mf6/mt'',i3)') mtd
   call mess('plotr','no distribution, no plot',strng)
   call tosend(nin,0,0,a)
   if (iauto.gt.0) go to 321
   go to 110

   !--mf4 and mf6 angular distributions
  2410 continue
   ltt=l2h
   call listio(nin,0,0,a,nb,nw)
   do while (nb.ne.0)
      call moreio(nin,0,0,a,nb,nw)
   enddo
  2415 continue
   call ad3d(nin,nplt,iplot,iauto,mfd,mtd,ltt,itype,jtype,igrid,&
     xleft,xright,xstep,ybot,ytop,ystep,rbot,rtop,rstep,&
     t1,t2,xl,yl,rl,nx,ny,nr,&
     xv,yv,zv,x3,y3,z3,&
     xll,yll,ww,wh,wr,iwcol,&
     aa,maxaa,a,nwamax)
   if (iauto.gt.0) go to 321
   go to 110

   !--mf5 or mf6 energy distribution plots
  2510 continue
   call tab1io(nin,0,0,a,nb,nw)
   lf=l2h
   if (lf.eq.1) go to 2520
   call mess('plotr','can only plot mf5/lf1',' ')
   call tosend(nin,0,0,a)
   if (iauto.gt.0) go to 321
   go to 110
  2520 continue
   call ed3d(nin,nplt,iplot,iauto,mfd,mtd,law,itype,jtype,igrid,&
      xleft,xright,xstep,ybot,ytop,ystep,rbot,rtop,rstep,&
      t1,t2,xl,yl,rl,nx,ny,nr,&
      xv,yv,zv,x3,y3,z3,&
      xll,yll,ww,wh,wr,iwcol,&
      aa,maxaa,a,nwamax)
   if (iauto.gt.0) go to 321
   go to 110

   !--mf7 s(alapha,beta) 3d plots
  2610 continue
   call error('plotr', '3d mf7 plots not available',' ')
   go to 110

   !--gendf data
   !--locate desired material and temperature
  400 continue
   call repoz(nin)
   call tpidio(nin,0,0,a,nb,nw)
   if (nin2.ne.0) then
      call repoz(nin2)
      call tpidio(nin2,0,0,a,nb,nw)
   endif
   jnoth=0
  410 continue
   call contio(nin,0,0,a,nb,nw)
   if (math.eq.matd) go to 430
   if (math.lt.0) then
      write(strng,'(''2desired mat and temp not found '',i4,f10.1)')&
        matd,temper
      call error('plotr',strng,' ')
    endif
   if (mfh.eq.0) go to 410
   call tomend(nin,0,0,a)
   go to 410
  430 continue
   loc=7
   call listio(nin,0,0,a(loc),nb,nw)
   temp=a(7)
   test=1
   test=test/1000
   if (abs(temp-temper).lt.test) go to 440
   call tomend(nin,0,0,a)
   go to 410
  440 continue
   do while (nb.ne.0)
      loc=loc+nw
      call moreio(nin,0,0,a(loc),nb,nw)
   enddo
   if (nin2.eq.0) go to 444
  441 continue
   loc=loc+nw
   loc2=loc
   call contio(nin2,0,0,a(loc),nb,nw)
   if (math.eq.matd2) go to 442
   if (math.lt.0) then
      write(strng,'(''3desired mat and temp not found '',i4,f10.1)')&
        matd,temper
      call error('plotr',strng,' ')
    endif
   if (mfh.eq.0) go to 441
   call tomend(nin2,0,0,a(loc))
   go to 441
  442 continue
   loc=loc+6
   call listio(nin2,0,0,a(loc),nb,nw)
   temp=a(loc)
   test=1
   test=test/1000
   if (abs(temp-temper).lt.test) go to 443
   call tomend(nin2,0,0,a(loc))
   go to 441
  443 continue
   do while (nb.ne.0)
      loc=loc+nw
      call moreio(nin2,0,0,a(loc),nb,nw)
   enddo
  444 continue
   nsigz=nint(a(4))
   ntw=nint(a(6))
   ngn=nint(a(9))
   ngg=nint(a(10))
   jbase=12
   jbase=jbase+ntw
   jbase=jbase+nsigz
   locngn=jbase+1
   jbase=jbase+ngn+1
   locngg=jbase+1
   jbase=jbase+ngg+1
   jbase=jbase+1
   if (i3d.eq.1) go to 3400
   nl=0
   nz=0

   !--locate desired list records
   ig1=0
   ig=ngn
  455 continue
   if (ig.lt.ngn) go to 465
   call contio(nin,0,0,a(jbase),nb,nw)
   if (nin2.ne.0) call contio(nin2,0,0,a(jbase+6),nb,nw)
   if (math.eq.0) go to 610
   if (mfh.gt.mfd) go to 610
   if (mfh.eq.0.or.mth.eq.0) go to 455
   nl=l1h
   nz=l2h
  465 continue
   call listio(nin,0,0,a(jbase),nb,nw)
   ng2=l1h
   ig2lo=l2h
   ig=n2h
   nwa=jbase+nw
   do while (nb.ne.0)
      call moreio(nin,0,0,a(nwa),nb,nw)
      nwa=nwa+nw
      if (nw+npage+12.gt.nwamax)&
          call error('plotr','storage exceeded2',' ')
   enddo
   if (nin2.eq.0) go to 467
   loc2=nwa+nw
   nwa=loc2
   call listio(nin2,0,0,a(nwa),nb,nw)
   nwa=nwa+nw
   do while (nb.ne.0)
      call moreio(nin2,0,0,a(nwa),nb,nw)
      nwa=nwa+nw
      if (nw+npage+12.gt.nwamax)&
          call error('plotr','storage exceeded3',' ')
   enddo
  467 continue

   !--process desired records
   if (mth.ne.mtd) go to 455
   if (mfd.eq.6.and.mfh.eq.3) go to 495
   if (mfh.ne.mfd) go to 455
   if (mfh.eq.6) go to 490
   if (mfh.ne.3) go to 455
   if (ig1.eq.0) ig1=ig
   jnow=2*(ig-ig1)
   x(jnow+1)=a(locngn+ig-1)*factx
   x(jnow+2)=a(locngn+ig)*factx
   if (ntp.gt.0) then
      if (ntp.lt.20) then
         s=a(jbase+6+nth*nl*nz+(nkh-1)+nl*(ntp-1))
      else if (ntp.lt.40) then
         s=0
         if (a(jbase+6+nth*nl*nz+(nkh-1)+nl*(ntp-21)).ne.zero) then
            s=100*(a(loc2+6+nth*nl*nz+(nkh-1)+nl*(ntp-21))&
              -a(jbase+6+nth*nl*nz+(nkh-1)+nl*(ntp-21)))/&
              a(jbase+6+nth*nl*nz+(nkh-1)+nl*(ntp-21))
         endif
      else
         s=0
         if (a(jbase+6+nth*nl*nz+(nkh-1)+nl*(ntp-41)).ne.zero) then
            s=a(loc2+6+nth*nl*nz+(nkh-1)+nl*(ntp-41))/&
              a(jbase+6+nth*nl*nz+(nkh-1)+nl*(ntp-41))
         endif
      endif
   else
      if (ntp.gt.-20) then
         s=0
         if (a(jbase+6+nth*nl*nz+(nkh-1)).ne.zero) then
            s=a(jbase+6+nth*nl*nz+(nkh-1)+nl*(-ntp-1))/&
              a(jbase+6+nth*nl*nz+(nkh-1))
         endif
      else if (ntp.gt.-40) then
         s=0
         if (a(loc2+6+nth*nl*nz+(nkh-1)).ne.zero) then
            s=100*((a(loc2+6+nth*nl*nz+(nkh-1)+nl*(-ntp-21))/&
              a(loc2+6+nth*nl*nz+(nkh-1)))&
              -(a(jbase+6+nth*nl*nz+(nkh-1)+nl*(-ntp-21))/&
              a(jbase+6+nth*nl*nz+(nkh-1))))/&
              (a(jbase+6+nth*nl*nz+(nkh-1)+nl*(-ntp-21))/&
              a(jbase+6+nth*nl*nz+(nkh-1)))
         endif
      else
         s=0
         if (a(loc2+6+nth*nl*nz+(nkh-1)).ne.zero) then
            s=(a(loc2+6+nth*nl*nz+(nkh-1)+nl*(-ntp-41))/&
              a(loc2+6+nth*nl*nz+(nkh-1)))/&
              (a(jbase+6+nth*nl*nz+(nkh-1)+nl*(-ntp-41))/&
              a(jbase+6+nth*nl*nz+(nkh-1)))
         endif
      endif
   endif
   if (nth.eq.0) then
      dleth=log(a(locngn+ig)/a(locngn+ig-1))
      s=s/dleth
   endif
   y(jnow+1)=s*facty
   y(jnow+2)=s*facty
   if (s.ne.zero) n=jnow+2
   go to 455
  490 continue
   if (ig.eq.nth) then
      do j=2,ng2
         jj=ig2lo+j-2
         jnow=2*(j-2)
         x(jnow+1)=a(locngn+jj-1)*factx
         x(jnow+2)=a(locngn+jj)*factx
         s=a(jbase+6+nl*nz*(j-1))
         dener=a(locngn+jj)-a(locngn+jj-1)
         s=s/dener
         s=s/sigig
         y(jnow+1)=s*facty
         y(jnow+2)=s*facty
         if (s.ne.zero) n=jnow+2
      enddo
   endif
   go to 455
  495 continue
   if (ig.ne.nth) go to 455
   sigig=a(jbase+6+nl*nz)
   go to 455

   !--gendf 3d plots
  3400 continue
   call contio(nin,0,0,a(jbase),nb,nw)
   if (mfh.eq.0.or.mth.eq.0) go to 3400
   if (math.le.0) then
      call mess('plotr',&
        'desired gendf matrix not found',' ')
      go to 110
   end if
   if (mfh.eq.mfd.and.mth.eq.mtd) go to 3410
   call tosend(nin,0,0,a(jbase))
   go to 3400
  3410 continue
   call gg3d(nin,nplt,iplot,iauto,mfd,mtd,itype,jtype,igrid,&
     xleft,xright,xstep,ybot,ytop,ystep,rbot,rtop,rstep,&
     t1,t2,xl,yl,rl,nx,ny,nr,&
     xv,yv,zv,x3,y3,z3,&
     xll,yll,ww,wh,wr,iwcol,&
     ngn,ngg,locngn,locngg,jbase,&
     aa,maxaa,a,nwamax)
   if (iauto.gt.0) go to 321
   go to 110

   !--experimental data
  500 continue
   read(nsysi,*) nform
   write(nsyso,'(/&
     &'' data nform ........................... '',i10)')&
     nform

   !--for format zero
   !--loop thru input lines until an empty card is found
   i=0
   ierrb=0
   idone=0
   do while (idone.eq.0)
      i=i+1
      nz=6
      flag=-99
      z(1)=flag
      z(2)=flag
      z(3)=0
      z(4)=0
      z(5)=0
      z(6)=0
      read(nsysi,*) (z(j),j=1,6)
      if (z(1).eq.flag.and.z(2).eq.flag) then
         idone=1
      else
         x(i)=z(1)*factx
         y(i)=z(2)*facty
         dym(i)=z(3)*facty
         dyp(i)=z(4)*facty
         if (dyp(i).eq.zero) dyp(i)=dym(i)
         dxm(i)=z(5)*factx
         dxp(i)=z(6)*factx
         if (dxp(i).eq.zero) dxp(i)=dxm(i)
         if (dym(i).ne.zero) ierrb=1
         if (dxm(i).ne.zero) ierrb=2
         if (dym(i).ne.zero.and.dxm(i).ne.zero) ierrb=3
      endif
   enddo
   n=i-1

   !--do plot
  610 continue
   factx1=1
   facty1=1
   write(nplt,'(i4,i8,7f7.2,''/ 2d plot'')') iplot,iwcol,&
     factx1,facty1,xll,yll,ww,wh,wr
   if (iabs(iplot).le.1) then
      if (nx.eq.0) xl=xlabld
      if (ny.eq.0) yl=ylabld
      write(nplt,'(1x,a,a,a,''/'')') qu,t1,qu
      if (iauto.gt.0) call rname(mtd,name)
      if (iauto.gt.0) write(t2,'(''mf='',i2,''  mt='',i3,2x,a)')&
        mfd,mtd,name
      write(nplt,'(1x,a,a,a,''/'')') qu,t2,qu
      write(nplt,'(4i6,1p,2e13.4,''/'')') itype,jtype,igrid,ileg,&
        xtag,ytag
      if (xstep.eq.zero) then
         write(nplt,'(''/'')')
      else
         write(nplt,'(1p,3e13.4,''/'')') xleft,xright,xstep
      endif
      write(nplt,'(1x,a,a,a,''/'')') qu,xl,qu
      if (ystep.eq.zero) then
         write(nplt,'(''/'')')
      else
         write(nplt,'(1p,3e13.4,''/'')') ybot,ytop,ystep
      endif
      write(nplt,'(1x,a,a,a,''/'')') qu,yl,qu
      if (jtype.gt.0) then
         if (rstep.eq.zero) then
            write(nplt,'(''/'')')
         else
            write(nplt,'(1p,3e13.4,''/'')') rbot,rtop,rstep
         endif
         write(nplt,'(1x,a,a,a,''/'')') qu,rl,qu
      endif
   endif
   write(nplt,'(''/'')')
   write(nplt,'(5i6,i8,''/'')') icon,isym,idash,iccol,ithick,ishade
   if (ileg.ne.0) write(nplt,'(1x,a,a,a,''/'')') qu,aleg,qu
   if (ileg.eq.2) write(nplt,'(1p,3e13.4,''/'')') xtag,ytag,xpoint
   write(nplt,'('' 0/'')')
   do i=1,n
      if (ierrb.eq.0) then
         write(nplt,'(1p,2e16.8,''/'')') x(i),y(i)
      else
         write(nplt,'(1p,6e16.8,''/'')') x(i),y(i),&
           dym(i),dyp(i),dxm(i),dxp(i)
      endif
   enddo
   write(nplt,'('' /'')')
   if (jnoth.eq.0) ipass=0
   if (jnoth.eq.1.and.ipass.eq.5) ipass=0
   if (iauto.gt.0.and.ipass.gt.0) go to 320
   if (iauto.gt.0) go to 321

   !--loop back for next curve
   go to 110

   !--plotr is finished
  700 continue
   write(nplt,'('' 99/'')')
   call closz(nplt)
   call closz(nin)
   call timer(time)
   write(nsyso,'(/69x,f8.1,''s''/&
     &1x,7(''**********''),''*******'')') time
   return
   end subroutine plotr

   subroutine ad3d(nin,nplt,iplot,iauto,mfd,mtd,ltt,itype,jtype,igrid,&
     xleft,xright,xstep,ybot,ytop,ystep,rbot,rtop,rstep,&
     t1,t2,xl,yl,rl,nx,ny,nr,&
     xv,yv,zv,x3,y3,z3,&
     xll,yll,ww,wh,wr,iwcol,&
     aa,maxaa,a,nwamax)
   !-------------------------------------------------------------------
   ! Plot a 3D perspective view of an angular distribution from
   ! the TAB2 structure in File 4 or File 6.
   !-------------------------------------------------------------------
   use util ! provides mess
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nplt,iplot,iauto,mfd,mtd,ltt
   integer::itype,jtype,igrid,nx,ny,nr,iwcol,maxaa,nwamax
   real(kr)::xleft,xright,xstep,ybot,ytop,ystep,rbot,rtop,rstep
   character(60)::t1,t2,xl,yl,rl
   real(kr)::xv,yv,zv,x3,y3,z3
   real(kr)::xll,yll,ww,wh,wr
   real(kr)::aa(*),a(*)
   ! internals
   integer::nb,nw,ne,locn,nmu,nmum,i,lang,nthind,ie,nl,intmu
   integer::loc,k,imu,il,nmax
   real(kr)::xmin,xmax,xstp,base,dmu,emin,emax,estep,varlim
   real(kr)::en,amu,sum,pl,pn,g,h,var,factx1,facty1
   real(kr)::ey3(200)
   character(10)::name
   character(60)::strng
   character(1)::qu=''''
   character(14)::xlabld='<e>nergy (e<v)'
   character(6)::prob='<p>rob'
   character(8)::cosine='<c>osine'
   integer,parameter::maxxy=200
   integer,parameter::maxx3=400
   real(kr),parameter::big=1.e8_kr
   real(kr),parameter::emev=1.e6_kr
   real(kr),parameter::zero=0

   call tab2io(nin,0,0,a,nb,nw)
   ne=n2h
   xmin=ybot
   xmax=ytop
   xstp=ystep
   if (xmax.eq.zero) xmax=big
   locn=0
   base=1
   base=base/100
   nmu=101
   nmum=nmu+2
   dmu=2
   dmu=dmu/(nmu-1)
   i=0
   emin=0
   emax=0
   lang=0
   estep=emev/5
   varlim=8
   varlim=varlim/100
   nthind=0
   do ie=1,ne
      if (mfh.ne.4.or.ltt.ne.2) then
         call listio(nin,0,0,a,nb,nw)
         lang=l1h
         nl=n1h
      else
         call tab1io(nin,0,0,a,nb,nw)
         intmu=4
         loc=nw+1
         do while (nb.ne.0)
            call moreio(nin,0,0,a(loc),nb,nw)
            loc=loc+nw
         enddo
      endif
      k=nw-3
      en=c2h
      if (en.ge.xmin-xmin/1000.and.en.le.xmax+xmax/1000) then
         if (locn+3*nmum.le.maxaa) then
            if (i+3.le.maxx3.and.i.lt.maxxy) then
               if (i.eq.0) emin=en
               emax=en
               i=i+1
               ey3(i)=en
               aa(locn+1)=base
               locn=locn+1
               do imu=1,nmu
                  amu=1-dmu*(imu-1)
                  if ((mfh.eq.4.and.ltt.eq.2).or.&
                    (mfh.eq.6.and.lang.gt.2)) then
                     do while (amu.lt.a(k))
                        k=k-2
                     enddo
                     call terp1(a(k),a(k+1),a(k+2),a(k+3),&
                       amu,sum,intmu)
                  else
                     pl=1
                     pn=amu
                     sum=1
                     sum=sum/2
                     do il=1,nl
                        sum=sum+(2*il+1)*pn*a(6+il)/2
                        g=amu*pn
                        h=g-pl
                        pl=pn
                        pn=h+g-h/(il+1)
                     enddo
                  endif
                  if (sum.lt.base) sum=base
                  aa(locn+imu)=sum
               enddo
               locn=locn+nmu
               aa(locn+1)=base
               locn=locn+1
               nmax=i
               ! check for possible thinning
               if (i.ne.1.and.ie.ne.ne) then
                  if (ey3(i).le.ey3(i-1)+estep) then
                     var=0
                     do imu=1,nmu
                        var=var+(aa(locn-imu-nmum)-aa(locn-imu))**2
                     enddo
                     if (var.le.varlim) then
                        nthind=nthind+1
                        locn=locn-nmum
                        i=i-1
                        nmax=nmax-1
                     endif
                  endif
               endif
            endif
         endif
      endif
   enddo
   if (nthind.gt.0) then
      write(strng,'(''for mt='',i3)') mtd
      call mess('ad3d','mf4 incident energy grid thinned',strng)
   endif
   call tosend(nin,0,0,a)
   if (locn+3*nmum.gt.maxaa.or.i+3.gt.maxx3)&
     call mess('ad3d',&
     'too much 3d angular distribution data',&
     'energy range truncated')
   if (i.ge.maxxy) call mess('ad3d',&
     'too many 3d angular distribution energies',&
     'energy range truncated')
   if (nx.eq.0) xl=cosine
   if (ny.eq.0) yl=xlabld
   if (nr.eq.0) rl=prob
   factx1=1
   facty1=1
   write(nplt,'(i4,i8,7f7.2,''/ 3d plot'')') iplot,iwcol,&
     factx1,facty1,xll,yll,ww,wh,wr
   write(nplt,'(1x,a,a,a,''/'')') qu,t1,qu
   if (iauto.gt.0) call rname(mtd,name)
   if (iauto.gt.0) write(t2,'(''mf='',i2,''  mt='',i3,2x,a)')&
     mfd,mtd,name
   write(nplt,'(1x,a,a,a,''/'')') qu,t2,qu
   write(nplt,'(3i6,''/'')') -itype,jtype,igrid
   if (xstep.eq.zero) then
      write(nplt,'(''/'')')
   else
      write(nplt,'(1p,3e13.4,''/'')') xleft,xright,xstep
   endif
   write(nplt,'(1x,a,a,a,''/'')') qu,xl,qu
   if (ystep.eq.zero) then
      write(nplt,'(''/'')')
   else
      write(nplt,'(1p,3e13.4,''/'')') ybot,ytop,ystep
   endif
   write(nplt,'(1x,a,a,a,''/'')') qu,yl,qu
   if (rstep.eq.zero) then
      write(nplt,'(''/'')')
   else
      write(nplt,'(1p,3e13.4,''/'')') rbot,rtop,rstep
   endif
   write(nplt,'(1x,a,a,a,''/'')') qu,rl,qu
   write(nplt,'(''/'')')
   write(nplt,'(6f9.3,''/'')') xv,yv,zv,-x3,y3,z3
   write(nplt,'('' 1/'')')
   locn=0
   do i=1,nmax
      write(nplt,'(1p,e13.4,''/'')') ey3(i)
      write(nplt,'(1p,2e13.4,''/'')') 1.,aa(locn+1)
      locn=locn+1
      do imu=1,nmu
         amu=1-dmu*(imu-1)
         write(nplt,'(1p,2e13.4,''/'')') amu,aa(locn+imu)
      enddo
      locn=locn+nmu
      write(nplt,'(1p,2e13.4,''/'')') -1.,aa(locn+1)
      locn=locn+1
      write(nplt,'('' /'')')
   enddo
   write(nplt,'('' /'')')
   return
   end subroutine ad3d

   subroutine ed3d(nin,nplt,iplot,iauto,mfd,mtd,law,itype,jtype,igrid,&
      xleft,xright,xstep,ybot,ytop,ystep,rbot,rtop,rstep,&
      t1,t2,xl,yl,rl,nx,ny,nr,&
      xv,yv,zv,x3,y3,z3,&
      xll,yll,ww,wh,wr,iwcol,&
      aa,maxaa,a,nwamax)
   !-------------------------------------------------------------------
   ! Plot a 3D perspective picture of a tabulated energy
   ! distribution from the TAB2 structure in File 5 or File 6.
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use util ! provides mess
   ! externals
   integer::nin,nplt,iplot,iauto,mfd,mtd,law
   integer::itype,jtype,igrid,nx,ny,nr,iwcol,maxaa,nwamax
   real(kr)::xleft,xright,xstep,ybot,ytop,ystep,rbot,rtop,rstep
   real(kr)::xv,yv,zv,x3,y3,z3
   real(kr)::xll,yll,ww,wh,wr
   character(60)::t1,t2,xl,yl,rl
   real(kr)::aa(*),a(*)
   ! internals
   integer::nb,nw,ne,i1,ne2,ne2m,locn,ie,idis,i2,nmax,ii1,ii2,i
   real(kr)::base,xmin,xmax,xstp,emin,emax,ymin,ymax,ystp
   real(kr)::e2,en,f2,ei,factx1,facty1,bigg,zmmm,zmin,zmax
   real(kr)::ex3(200),ey3(200)
   character(10)::name
   character(14)::xlabld='<e>nergy (e<v)'
   character(13)::sece='<s>ec. energy'
   character(10)::probp='<p>rob/e<v'
   character(1)::qu=''''
   integer,parameter::maxxy=200
   integer,parameter::maxx3=400
   integer::nw7max=6000
   real(kr),parameter::eps=1.e-7_kr
   real(kr),parameter::big=1.e8_kr
   real(kr),parameter::small=1.e-12_kr
   real(kr),parameter::emev=1.e6_kr
   real(kr),parameter::elow=1.e-5_kr
   real(kr),parameter::etop=1.e10_kr
   real(kr),parameter::zero=0
   real(kr),parameter::ten=10

   !--this tab1io call should be commented out in the unlikely
   !  event that an njoy99/thermr, or earlier, file is being processed.
   if (mfd.eq.6.and.mtd.ge.221.and.mtd.le.250) &
                                            call tab1io(nin,0,0,a,nb,nw)
   !
   call tab2io(nin,0,0,a,nb,nw)
   ne=n2h
   i1=0
   base=small
   xmin=ybot
   xmax=ytop
   if (xmax.eq.0) xmax=big
   xstp=ystep
   emin=0
   emax=etop
   ymin=xleft
   ymax=xright
   if (mth.lt.221.or.mth.gt.250) then
      if (ymax.eq.zero) ymax=20*emev
      if (ymin.eq.zero) ymin=1
   else
      if (ymax.eq.zero) ymax=10
      if (ymin.eq.zero) ymin=eps
   endif
   ystp=xstep
   if (itype.eq.2.or.itype.eq.4) then
      ymin=log10(ymin)
      ymax=log10(ymax)
   endif
   ne2=198
   ne2m=ne2+2
   locn=0
   if (law.eq.7) locn=locn+nw7max
   zmax=small
   do ie=1,ne
      e2=0
      if (mfd.eq.5) call gety1(e2,en,idis,f2,nin,a)
      if (mfd.eq.6)&
        call gety6(1,e2,en,idis,f2,nin,a,2,law,aa,nw7max)
      ei=a(2)
      if (i1.lt.maxxy) then
         if (ei.ge.xmin-xmin/1000.and.ei.le.xmax+xmax/1000) then
            if (locn+3*ne2m.le.maxaa.and.i1+3.le.maxx3) then
               if (i1.eq.0) emin=ei
               i1=i1+1
               ex3(i1)=ei
               if (y3.gt.zero) ey3(1)=ymin
               if (y3.lt.zero) ey3(1)=ymax
               aa(locn+1)=base
               locn=locn+1
               do i2=1,ne2
                  e2=ymin+(ymax-ymin)*(i2-1)/(ne2-1)
                  if (itype.eq.2.or.itype.eq.4) e2=ten**e2
                  if (e2.eq.zero) e2=elow
                  if (mfd.eq.5) call gety1(e2,en,idis,f2,nin,a)
                  if (mfd.eq.6) then
                     call gety6(1,e2,en,idis,f2,nin,a,2,&
                        law,aa,nw7max)
                  endif
                  if (f2.lt.base) f2=base
                  if (f2.gt.zmax) zmax=f2
                  aa(locn+i2)=f2
               enddo
               locn=locn+ne2
               aa(locn+1)=base
               locn=locn+1
               nmax=i1
            endif
         endif
      endif
      bigg=big
      if (mfd.eq.5) call gety1(bigg,en,idis,f2,nin,a)
      if (mfd.eq.6)&
        call gety6(1,bigg,en,idis,f2,nin,a,2,law,aa,nw7max)
   enddo
   call tosend(nin,0,0,a)
   if (locn+3*ne2m.gt.maxaa.or.i1+3.gt.maxx3)&
     call mess('ed3d',&
     'too much energy distribution data',&
     'energy range truncated')
   ! fix up z axis
   zmax=log10(zmax)
   zmmm=33
   zmmm=zmmm/10
   zmin=zmax-zmmm
   i=int(zmin)
   if (float(i).gt.zmin) i=i-1
   zmin=ten**i
   ! set up axis labels
   if (ny.eq.0) yl=xlabld
   if (nx.eq.0) xl=sece
   if (nr.eq.0) rl=probp
   factx1=1
   facty1=1
   write(nplt,'(i4,i8,7f7.2,''/ 3d plot'')') iplot,iwcol,&
     factx1,facty1,xll,yll,ww,wh,wr
   write(nplt,'(1x,a,a,a,''/'')') qu,t1,qu
   if (iauto.gt.0) call rname(mtd,name)
   if (iauto.gt.0) write(t2,'(''mf='',i2,''  mt='',i3,2x,a)')&
     mfd,mtd,name
   write(nplt,'(1x,a,a,a,''/'')') qu,t2,qu
   write(nplt,'(3i6,''/'')') -itype,jtype,igrid
   if (xstep.eq.zero) then
      write(nplt,'(''/'')')
   else
      write(nplt,'(1p,3e13.4,''/'')') xleft,xright,xstep
   endif
   write(nplt,'(1x,a,a,a,''/'')') qu,xl,qu
   if (ystep.eq.zero) then
      write(nplt,'(''/'')')
   else
      write(nplt,'(1p,3e13.4,''/'')') ybot,ytop,ystep
   endif
   write(nplt,'(1x,a,a,a,''/'')') qu,yl,qu
   if (rstep.eq.zero) then
      write(nplt,'(''/'')')
   else
      write(nplt,'(1p,3e13.4,''/'')') rbot,rtop,rstep
   endif
   write(nplt,'(1x,a,a,a,''/'')') qu,rl,qu
   write(nplt,'(''/'')')
   write(nplt,'(6f9.3)') xv,yv,zv,x3,y3,z3
   write(nplt,'('' 1/'')')
   locn=0
   if (law.eq.7) locn=locn+nw7max
   do i1=1,nmax
      write(nplt,'(1p,e13.4,''/'')') ex3(i1)
      ii1=0
      ii2=0
      do i2=1,ne2m
         if (ii1.eq.0.and.aa(locn+i2).ge.zmin) ii1=i2
         if (aa(locn+i2).ge.zmin) ii2=i2
         if (aa(locn+i2).lt.zmin) aa(locn+i2)=zmin
      enddo
      if (ii1.eq.0) ii1=2
      if (ii2.eq.0) ii2=2
      if (ii1.gt.1) ii1=ii1-1
      if (ii2.lt.ne2m) ii2=ii2+1
      do i2=ii1,ii2
         e2=ymin+(ymax-ymin)*(i2-2)/(ne2-1)
         if (i2.eq.1) e2=ymin
         if (i2.eq.ne2m) e2=ymax
         if (itype.eq.2.or.itype.eq.4) e2=ten**e2
         if (e2.eq.zero) e2=elow
         write(nplt,'(1p,2e13.4,''/'')') e2,aa(locn+i2)
      enddo
      locn=locn+ne2m
      write(nplt,'('' /'')')
   enddo
   write(nplt,'('' /'')')
   return
   end subroutine ed3d

   subroutine gg3d(nin,nplt,iplot,iauto,mfd,mtd,itype,jtype,igrid,&
     xleft,xright,xstep,ybot,ytop,ystep,rbot,rtop,rstep,&
     t1,t2,xl,yl,rl,nx,ny,nr,&
     xv,yv,zv,x3,y3,z3,&
     xll,yll,ww,wh,wr,iwcol,&
     ngn,ngg,locngn,locngg,jbase,&
     aa,maxaa,a,nwamax)
   !-------------------------------------------------------------------
   ! Plot a 3D perspective view of a group-to-group matrix.
   !-------------------------------------------------------------------
   use util ! provides error
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nplt,mfd,mtd,iplot,iauto
   integer::itype,jtype,igrid,nx,ny,nr,iwcol,maxaa,nwamax
   real(kr)::xleft,xright,xstep,ybot,ytop,ystep,rbot,rtop,rstep
   character(60)::t1,t2,xl,yl,rl
   real(kr)::xv,yv,zv,x3,y3,z3
   real(kr)::xll,yll,ww,wh,wr
   integer::ngn,ngg,locngn,locngg,jbase
   real(kr)::aa(*),a(*)
   ! internals
   integer::ng1,nl,nz,ngl,locn,ig,nb,nw,l,ng2,ig2lo,ngm
   integer::j2,i2,k2,j,ig2,nmax,i,i1,major,minor
   real(kr)::xmin,xmax,xstp,emin,emax,ymin,ymax,ystp
   real(kr)::zmin,zmax,zmmm
   real(kr)::e2,eglo,eghi,f2,factx1,facty1
   real(kr)::ex3(200),ey3(200)
   character(10)::name
   character(1)::qu=''''
   character(13)::sece='<s>ec. energy'
   integer::nsece=13
   character(10)::probp='<p>rob/e<v'
   integer::nprobp=10
   integer::maxx3=400
   integer::maxy3=400
   real(kr),parameter::small=1.e-12_kr
   real(kr),parameter::elow=1.e-5_kr
   real(kr),parameter::etop=1.e10_kr
   real(kr),parameter::zero=0

   ng1=0
   nl=nint(a(jbase+2))
   nz=nint(a(jbase+3))
   ngl=nint(a(jbase+5))
   xmin=xleft
   xmax=xright
   if (xmax.eq.zero) xmax=a(locngn+ngl)
   if (xmin.eq.zero) xmin=a(locngn)
   xstp=xstep
   emin=0
   emax=etop
   ymin=ybot
   ymax=ytop
   if (mfd.eq.16) then
      if (ymax.eq.zero) ymax=a(locngg+ngg)
      if (ymin.eq.zero) ymin=a(locngg)
   else
      if (ymax.eq.zero) ymax=a(locngn+ngn)
      if (ymin.eq.zero) ymin=a(locngn)
   endif
   ystp=ystep
   locn=0
   zmax=-12
   ig=1
   do while (ig.lt.ngl)
      ng1=ng1+1
      call listio(nin,0,0,a(jbase),nb,nw)
      l=jbase
      l=l+nw
      do while (nb.ne.0)
         call moreio(nin,0,0,a(l),nb,nw)
         l=l+nw
      enddo
      ng2=nint(a(jbase+2))
      ig2lo=nint(a(jbase+3))
      ig=nint(a(jbase+5))
      if (ng1+2.gt.maxx3) call error('gg3d',&
        'too many incident groups for 3d gendf plot',' ')
      ex3(ng1)=a(locngn+ig-1)
      if (ng1.eq.1) emin=ex3(ng1)
      ngm=ngn+1
      if (mfd.eq.16) ngm=ngg+1
      if (locn+2*ngm+4.gt.maxaa) call error('gg3d',&
        'too many data for 3d gendf plot',' ')
      ey3(1)=ymin
      ey3(2)=ey3(1)
      j2=2
      do i2=1,ngm
         e2=a(locngn+i2-1)
         if (mfd.eq.16) e2=a(locngg+i2-1)
         if (e2.gt.ymin.and.e2.lt.ymax) then
            if (j2+4.gt.maxy3) call error('gg3d',&
              'too many secondary groups for 3d gendf plot',' ')
            ey3(j2+1)=e2
            ey3(j2+2)=e2
            j2=j2+2
         endif
      enddo
      ey3(j2+1)=ymax
      ey3(j2+2)=ymax
      j2=j2+2
      do i2=1,j2
         aa(locn+i2)=-12
         aa(locn+j2+i2)=-12
      enddo
      locn=locn+j2
      ng1=ng1+1
      ex3(ng1)=a(locngn+ig-1)
      k2=1
      do j=2,ng2
         ig2=ig2lo+j-2
         eglo=a(locngn+ig2-1)
         if (mfd.eq.16) eglo=a(locngg+ig2-1)
         eghi=a(locngn+ig2)
         if (mfd.eq.16) eghi=a(locngg+ig2)
         if (eghi.gt.ymin.and.eglo.le.ymax) then
            f2=a(jbase+6+nl*nz*(j-1))/(eghi-eglo)
            if (f2.lt.small) f2=small
            f2=log10(f2)
            if (f2.gt.zmax) zmax=f2
            k2=k2+1
            aa(locn+k2)=f2
            k2=k2+1
            aa(locn+k2)=f2
         endif
      enddo
      k2=k2+1
      locn=locn+j2
   enddo
   ng1=ng1+1
   ex3(ng1)=a(locngn+ig)
   ng1=ng1+1
   ex3(ng1)=a(locngn+ig)+a(locngn+ig)/1000
   do i2=1,j2
      aa(locn+i2)=aa(locn+i2-2*j2)
      aa(locn+i2+j2)=aa(locn+i2-j2)
   enddo
   locn=locn+2*j2
   nmax=ng1
   emax=ex3(ng1)
   ! fix up z axis
   i=int(zmax)
   if (i.lt.zmax) i=i+1
   zmax=i
   zmmm=33
   zmmm=zmmm/10
   zmin=zmax-zmmm
   i=int(zmin)
   if (i.gt.zmin) i=i-1
   zmin=i
   do i1=1,nmax
      do i2=1,j2
         l=i2+j2*(i1-1)
         if (aa(l).lt.zmin) aa(l)=zmin
      enddo
   enddo
   rbot=zmin
   rtop=zmax
   rstep=1
   ! fix up x axis
   if (itype.eq.3.or.itype.eq.4) then
      if (xmin.eq.zero) xmin=elow
      xmin=log10(emin)
      xmax=log10(emax)
      xstp=1
      do i1=1,nmax
         ex3(i1)=log10(ex3(i1))
      enddo
   endif
   if (xstp.eq.zero) then
      xmin=emin
      xmax=emax
      call ascale(4,xmin,xmax,major,minor)
      xstp=(xmax-xmin)/major
   endif
   xleft=xmin
   xright=xmax
   xstep=xstp
   ! fix up y axis
   if (itype.eq.2.or.itype.eq.4) then
      ymin=log10(ymin)
      ymax=log10(ymax)
      ystp=1
      do i2=1,j2
         ey3(i2)=log10(ey3(i2))
      enddo
   endif
   if (ystp.eq.zero) then
      call ascale(2,ymin,ymax,major,minor)
      ystp=(ymax-ymin)/major
   endif
   call tosend(nin,0,0,a)
   ybot=ymin
   ytop=ymax
   ystep=ystp
   ! set up axis labels
   rl=probp
   nr=nprobp
   jtype=2
   yl=sece
   ny=nsece
   factx1=1
   facty1=1
   write(nplt,'(i4,i8,7f7.2,''/ 3d plot'')') iplot,iwcol,&
     factx1,facty1,xll,yll,ww,wh,wr
   write(nplt,'(1x,a,a,a,''/'')') qu,t1,qu
   if (iauto.gt.0) call rname(mtd,name)
   if (iauto.gt.0) write(t2,'(''mf='',i2,''  mt='',i3,2x,a)')&
     mfd,mtd,name
   write(nplt,'(1x,a,a,a,''/'')') qu,t2,qu
   write(nplt,'(3i6,''/'')') -itype,jtype,igrid
   if (xstep.eq.zero) then
      write(nplt,'(''/'')')
   else
      write(nplt,'(1p,3e13.4,''/'')') xleft,xright,xstep
   endif
   write(nplt,'(1x,a,a,a,''/'')') qu,xl,qu
   if (ystep.eq.zero) then
      write(nplt,'(''/'')')
   else
      write(nplt,'(1p,3e13.4,''/'')') ybot,ytop,ystep
   endif
   write(nplt,'(1x,a,a,a,''/'')') qu,yl,qu
   if (rstep.eq.zero) then
      write(nplt,'(''/'')')
   else
      write(nplt,'(1p,3e13.4,''/'')') rbot,rtop,rstep
   endif
   write(nplt,'(1x,a,a,a,''/'')') qu,rl,qu
   write(nplt,'(''/'')')
   write(nplt,'(6f9.3,''/'')') xv,yv,zv,-x3,y3,z3
   write(nplt,'('' 1/'')')
   return
   end subroutine gg3d

   subroutine ascale(m,z1,z2,major,minor)
   !--------------------------------------------------------------------
   ! Automatic scaling routine for a linear axis
   ! borrowed from Los Alamos SC4020 library (with modifications).
   ! on input------------
   !   m              minimum number of major divisions desired
   !   z1,z2          min and max values of data to be plotted
   ! on output-----------
   !   z1,z2          min and max limits of axis
   !   major          number of major division on axis
   !   minor          number of minor divisions on axis
   !--------------------------------------------------------------------
   ! externals
   integer::m,major,minor
   real(kr)::z1,z2
   ! internals
   integer::iflag,n1,n2,k,nm
   real(kr)::zmin,zmax,fm,zbar,z,p,tenk,test,dz,fn
   real(kr),parameter::zero=0
   real(kr),parameter::one=1
   real(kr),parameter::ten=10
   real(kr),parameter::tentho=10000

   !--check input
   !--exit if parameters are unreasonable
   zmin=z1
   zmax=z2
   if (zmax.le.zmin.or.m.le.0.or.m.gt.20) then
      major=1
      minor=0
      z1=0
      z2=2*z2
      return
   endif

   !--find minimum span per interval
   fm=m
   if (zmax.eq.zero.or.zmin.eq.zero) then
      zmax=zmax-abs(zmax)/1000000
      zmin=zmin+abs(zmin)/1000000
   else
      zbar=zmax/zmin
      if (abs(zbar)/1000.ge.one) zmin=0
      if (abs(zbar)*1000.le.one) then
         zmax=0
         zmax=zmax-abs(zmax)/1000000
         zmin=zmin+abs(zmin)/1000000
      else
         if ((abs(zbar-1)-5*fm/100000).gt.zero) then
            zmax=zmax-abs(zmax)/1000000
            zmin=zmin+abs(zmin)/1000000
         else
            zbar=(zmax+zmin)/2
            z=26*fm*abs(zbar)/1000000
            zmax=zbar+z
            zmin=zbar-z
         endif
      endif
   endif
   p=(zmax-zmin)/fm

   !--determine k such that
   !--    10**k.le.p.lt.10**(k+1)
   iflag=0
   tenk=1
   k=0
   if (p.lt.one) then
      iflag=1
      p=1/p
   endif
   do while (p.ge.tentho)
      p=p/10000
      tenk=tenk*10000
      k=k+4
   enddo
   do while (p.ge.ten)
      p=p/10
      tenk=tenk*10
      k=k+1
   enddo
   if (iflag.ne.0) then
      p=10/p
      tenk=one/10/tenk
      k=k-1
   endif

   !--determine dz such that
   !--    dz = 1*10**k, or
   !--    dz = 2*10**k, or
   !--    dz = 5*10**k.
   test=2+one/100
   if (p.le.test) then
      test=2-one/100
      if (p.le.test) then
         p=1
         nm=5
      else
         p=2
         nm=4
      endif
   else
      test=5-one/100
      if (p.lt.test) then
         p=2
         nm=4
      else
         p=5
         nm=5
      endif
   endif
   dz=p*tenk

   !--find integers n1 and n2 such that
   !--    n1*dz.le.zmin.lt.(n1+1)*dz, and
   !--    (n2-1)*dz.lt.zmax.le.n2*dz
   !--the new zmin and zmax become
   !--    zmin = n1*dz
   !--    zmax = n2*dz
   n1=int(zmin/dz)
   fn=n1
   z=fn*dz
   if (z.gt.zmin) then
      z=z-dz
      n1=n1-1
   endif
   zmin=z
   z1=zmin
   n2=int(zmax/dz)
   fn=n2
   z=fn*dz
   if (z.lt.zmax) then
      n2=n2+1
      z=z+dz
   endif
   zmax=z
   z2=zmax

   !--the values of major and minor can be set
   major=n2-n1
   minor=nm*major
   return
   end subroutine ascale

   subroutine getz(x,xnext,idis,z,itape,a)
   !--------------------------------------------------------------------
   ! Retrieve z(x) from an ENDF/B TAB1 structure using paged bcd or
   ! blocked binary formats.  Call with x=0 to read in first page
   ! or block of data and initialize pointers.  Routine assumes
   ! values will be called in ascending order.  Xnext is the first
   ! data grid point greater than x unless x is the last point.
   ! This version will keep track of pointers for up to 10 units.
   ! Call with x=-1 to clear the pointers before each group of files.
   ! Based on gety from mixr.
   !--------------------------------------------------------------------
   use util ! provides error
   use endf ! provides endf routines and variables
   ! externals
   real(kr)::x,xnext,z
   real(kr)::a(*)
   integer::idis,itape
   ! internals
   integer::ntape,jtape(10),nrt(10),npt(10),irt(10),ipt(10)
   integer::ip1t(10),ip2t(10),nbt(10),nwt(10)
   save ntape,jtape,nrt,npt,irt,ipt,ip1t,ip2t,nbt,nwt
   integer::lt
   save lt
   integer::nb,nw,nwtot,nr,np,ip1,ip2,ir,ip,i,ktape,ln,nbx,int
   real(kr)::xn
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::zero=0

   !--branch on value of x
   idis=0
   if (x.eq.zero) go to 100
   if (x.gt.zero) go to 115

   !--clear pointer storage
   ntape=0
   return

   !--read first page or block of data and initialize
  100 ntape=ntape+1
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
  115 if (ntape.eq.0)&
     call error('getz','not properly initialized',' ')
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
   z=0
   xnext=big
   return

   !--is x in this panel
  125 ln=2*(ip-ip1)+lt
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
  130 int=nint(a(6+2*ir))
   if (int.eq.1) idis=1
   call terp1(a(ln-1),a(ln),a(ln+1),a(ln+2),x,z,int)
   xnext=a(ln+1)
   if ((ln+3).gt.nwtot.and.nb.eq.0) return
   xn=a(ln+3)
   if (xn.eq.xnext) idis=1
   go to 150

   !--special branch for x outside range of table
  135 z=0
   xnext=a(ln-1)
   go to 150

   !--special branch for last point
  140 z=a(ln+2)
   xnext=big

   !--save pointers and return
  150 nrt(ktape)=nr
   npt(ktape)=np
   irt(ktape)=ir
   ipt(ktape)=ip
   ip1t(ktape)=ip1
   ip2t(ktape)=ip2
   nbt(ktape)=nb
   nwt(ktape)=nwtot
   return
   end subroutine getz

   subroutine gety4(x,xnext,idis,y4,nth,nin,a)
   !--------------------------------------------------------------------
   ! Retrieve y4(x) for a MF4 list structure.
   ! Call with x=0 to initialize.
   !--------------------------------------------------------------------
   use util ! provides error
   use endf ! provides endf routines and variables
   ! externals
   integer::idis,nth,nin
   real(kr)::x,xnext,y4,a(*)
   ! internals
   integer::np,ip,nal,na,int,nb,nw,ltt,li
   real(kr)::xlast,ylast,xnow,ynow
   real(kr),parameter::big=1.e10_kr
   save xlast,ylast,np,ip,na,nal

   !--read first point and initialize
   if (x.eq.0) then
      ltt=nint(a(4))
      call contio(nin,0,0,a,nb,nw)
      li=nint(a(4))
      if (ltt.ne.1.and.ltt.ne.3.or.li.eq.1) call error('gety4',&
        'only Legendre coefficient representations supported',' ')
      call tab2io(nin,0,0,a,nb,nw)
      np=nint(a(6))
      int=nint(a(8))
      call listio(nin,0,0,a,nb,nw)
      na=nint(a(5))
      xlast=a(2)
      ylast=0
      nal=na
      if (nth.le.na) ylast=a(6+nth)
      ip=2
      xnext=xlast
   endif

   !--is x is in this panel?
  110 continue
   if (x.lt.xlast) go to 130
   if (x.lt.a(2)) go to 120
   if (ip.eq.np) go to 140

   !--no.  move up to next range.
   xlast=a(2)
   ylast=0
   if (nth.le.na) ylast=a(6+nth)
   nal=na
   ip=ip+1
   call listio(nin,0,0,a,nb,nw)
   na=nint(a(5))
   go to 110

   !--yes. interpolate.
  120 continue
   int=2
   xnow=a(2)
   ynow=0
   if (nth.le.na) ynow=a(6+nth)
   call terp1(xlast,ylast,xnow,ynow,x,y4,int)
   xnext=xnow
   return

   !--special branch for x outside range of table
  130 y4=0
   xnext=xlast
   return

   !--special branch for last point
  140 y4=a(6+nth)
   xlast=big
   xnext=big

   return
   end subroutine gety4

   subroutine gety6(i6,x,xnext,idis,y6,itape,a,lep,law,aa,maa)
   !--------------------------------------------------------------------
   ! Retrieve y6(x) from an MF6 LIST structure using paged BCD or
   ! blocked binary formats.  Call with x=0 to read in first page
   ! or block of data and initialize pointers.  Routine assumes
   ! values will be called in ascending order.  Here, next is the first
   ! data grid point greater than x unless x is the last point, and
   ! i6 is the index to the dependent variable in the cycle
   ! of na+2 values per secondary energy (i6 is usually 1).
   ! For law7, the angle-energy data is converted to a simple energy
   ! distribution and stored in a to be handled in the normal way.
   !--------------------------------------------------------------------
   use endf ! provides endf routines and variables
   ! externals
   integer::i6,idis,itape,lep,law,maa
   real(kr)::x,xnext,y6,a(*),aa(*)
   ! internals
   integer::nb,nw,na,nwp,ln,nc,i,int,nwtot,nep,ncyc,lt,ip1,ip2,ip
   real(kr)::xlast,ylast,xn
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::zero=0
   save xlast,nwtot,nep,ncyc,lt,ip1,ip2,ip,nb

   !--read first page or block of data and initialize.
   idis=0
   if (x.eq.zero) then
      if (law.eq.7) call fixl7(itape,aa,maa,a,nb,nw)
      if (law.ne.7) call listio(itape,0,0,a,nb,nw)
      nwtot=nw
      if (mth.lt.221.or.mth.gt.250) then
         na=l2h
         nwp=n1h
         nep=n2h
         ncyc=na+2
      else
         ncyc=n2h
         nwp=n1h
         nep=nwp/ncyc
      endif
      lt=7
      xlast=a(lt)
      ylast=a(lt+i6)
      ip1=1
      ip2=(nw-6)/ncyc
      ip=2
      xnext=a(lt)
      return
   endif

   !--is x in this panel
  110 if (x.lt.xlast) go to 130
   ln=ncyc*(ip-ip1)+lt
   if (x.lt.a(ln)) go to 120
   if (ip.eq.nep) go to 140

   !--no.  move up to next range.
   !--read in new page of data if needed.
   xlast=a(ln)
   ylast=a(ln+i6)
   ip=ip+1
   if ((ip+2).le.ip2) go to 110
   if (nb.eq.0) go to 120
   nc=nwtot-ln+1
   do i=1,nc
      a(lt+i-1)=a(nwtot-nc+i)
   enddo
   call moreio(itape,0,0,a(lt+nc),nb,nw)
   nwtot=nw+lt+nc-1
   ip1=ip
   ip2=ip1+nw/ncyc+1
   if (nb.eq.0) ip2=ip2+2
   go to 110

   !--yes.  interpolate for desired value
  120 int=lep
   if (int.eq.1) idis=1
   call terp1(xlast,ylast,a(ln),a(ln+i6),x,y6,int)
   xnext=a(ln)
   if ((ln+ncyc).gt.nwtot.and.nb.eq.0) return
   xn=a(ln+ncyc)
   if (xn.eq.xnext) idis=1
   return

   !--special branch for x outside range of table
  130 y6=0
   xnext=xlast
   return

   !--special branch for last point
  140 y6=a(ln+i6)
   xlast=big
   xnext=big
   return
   end subroutine gety6

   subroutine fixl7(itape,aa,maa,a,nb,nw)
   !--------------------------------------------------------------------
   ! Convert the TAB2/TAB1 angle-energy structure of File6, law7 to
   ! a simple energy distribution in file6, law1 form for plotting.
   !--------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use util ! provides error
   ! externals
   integer::itape,maa,nb,nw
   real(kr)::aa(*),a(*)
   ! internals
   integer::nmu,imu,ir,loc,ne,loca,nm,loc1,loc2,ip,ie,idis
   real(kr)::ei,emin,emax,de,e,u1,u2,f1,f2,en
   real(kr)::egrid(200)
   real(kr),parameter::etop=1.e10_kr
   real(kr),parameter::span=1.e6_kr

   !--set up angle loop
   call tab2io(itape,0,0,aa,nb,nw)
   ei=aa(2)
   nmu=nint(aa(6))

   !--read all angles for this energy, finding limits
   emin=etop
   emax=0
   loc=1+nmu
   do imu=1,nmu
      aa(imu)=loc
      call tab1io(itape,0,0,aa(loc),nb,nw)
      ir=nint(aa(loc+4))
      if (aa(loc+6+2*ir).lt.emin) emin=aa(loc+6+2*ir)
      loc=loc+nw
      do while (nb.ne.0)
         if (loc.gt.maa) call error('fixl7',&
           'not enough storage to convert file 7',' ')
         call moreio(itape,0,0,aa(loc),nb,nw)
         loc=loc+nw
      enddo
      if (aa(loc-2).gt.emax) emax=aa(loc-2)
   enddo
   if (emin.lt.emax/span) emin=emax/span

   !--set up a secondary energy grid
   !--and do a trapazoidal angle integration for each grid point
   de=(emax-emin)/150
   ne=1
   e=emin
   egrid(ne)=e
   do while (e.lt.emax)
      ne=ne+1
      if (e.lt.emin+10*de) then
         e=e+e/4
         if (e.gt.emin+10*de) e=emin+10*de
      else
         e=e+de
      endif
      egrid(ne)=e
   enddo
   egrid(ne)=emax
   loca=6
   nm=nmu-1
   do imu=1,nm
      loc1=nint(aa(imu))
      loc2=nint(aa(imu+1))
      u1=aa(loc1+1)
      u2=aa(loc2+1)
      ip=2
      ir=1
      do ie=1,ne
         e=egrid(ie)
         call terpa(f1,e,en,idis,aa(loc1),ip,ir)
         call terpa(f2,e,en,idis,aa(loc2),ip,ir)
         if (imu.ne.1) then
            a(loca+2*ie)=a(loca+2*ie)+(u2-u1)*(f1+f2)/2
         else
            a(loca-1+2*ie)=e
            a(loca+2*ie)=(u2-u1)*(f1+f2)/2
         endif
      enddo
   enddo
   a(1)=0
   a(2)=ei
   a(3)=0
   a(4)=0
   a(5)=2*ne
   a(6)=ne

   !--finished
   nw=6+2*ne
   nb=0
   return
   end subroutine fixl7

   subroutine rname(mt,name)
   !--------------------------------------------------------------------
   ! Return the reaction name for an ENDF MT number.
   !--------------------------------------------------------------------
   use endf ! provides iverf
   ! externals
   integer::mt
   character(10)::name
   ! internals
   integer::i
   character(10),dimension(415)::hndf=(/&
     'total     ', 'elastic   ', 'nonelastic', 'inelastic ', '(n,x)     ',&
     '(n,2n_1f) ', '(n,2n_2f) ', '(n,2n_3f) ', '(n,2n_4f) ', '(n,x)     ',&
     '(n,2nd)   ', '(n,x)     ', '(n,x)     ', '(n,x)     ', '(n,x)     ',&
     '(n,2n)    ', '(n,3n)    ', 'fission   ', '(n,f)     ', "(n,n''f)  ",&
     '(n,2nf)   ', "(n,n''a)  ", "(n,n'')3a ", '(n,2n)a   ', '(n,3n)a   ',&
     '(n,2n)iso ', '(n,abs)   ', "(n,n''p)  ", "(n,n''a)  ", '(n,2n)2a  ',&
     '(n,x)     ', "(n,n''d)  ", "(n,n''t)  ", "(n,n''he3)", "(n,n'')d2a",&
     "(n,n'')t2a", '(n,4n)    ', '(n,3nf)   ', '(n,x)     ', '(n,x)     ',&
     '(n,2np)   ', '(n,3np)   ', '(n,n2p)   ', '(n,npa)   ', '(n,x)     ',&
     '(n,2n_1s) ', '(n,2n_2s) ', '(n,2n_3s) ', '(n,2n_4s) ', '(n,x)     ',&
     '(n,n_1)   ', '(n,n_2)   ', '(n,n_3)   ', '(n,n_4)   ', '(n,n_5)   ',&
     '(n,n_6)   ', '(n,n_7)   ', '(n,n_8)   ', '(n,n_9)   ', '(n,n_10)  ',&
     '(n,n_11)  ', '(n,n_12)  ', '(n,n_13)  ', '(n,n_14)  ', '(n,n_15)  ',&
     '(n,n_16)  ', '(n,n_17)  ', '(n,n_18)  ', '(n,n_19)  ', '(n,n_20)  ',&
     '(n,n_21)  ', '(n,n_22)  ', '(n,n_23)  ', '(n,n_24)  ', '(n,n_25)  ',&
     '(n,n_26)  ', '(n,n_27)  ', '(n,n_28)  ', '(n,n_29)  ', '(n,n_30)  ',&
     '(n,n_31)  ', '(n,n_32)  ', '(n,n_33)  ', '(n,n_34)  ', '(n,n_35)  ',&
     '(n,n_36)  ', '(n,n_37)  ', '(n,n_38)  ', '(n,n_39)  ', '(n,n_40)  ',&
     '(n,n_c)   ', '(n,x)     ', '(n,x)     ', '(n,x)     ', '(n,x)     ',&
     '(n,x)     ', '(n,x)     ', '(n,x)     ', "(n,n'')g  ", '(n,x)     ',&
     '(n,parab) ', '(n,g)     ', '(n,p)     ', '(n,d)     ', '(n,t)     ',&
     '(n,he3)   ', '(n,a)     ', '(n,2a)    ', '(n,3a)    ', '(n,x)     ',&
     '(n,2p)    ', '(n,pa)    ', '(n,t2a)   ', '(n,d2a)   ', '(n,pd)    ',&
     '(n,pt)    ', '(n,da)    ', '(n,x)     ', '(n,x)     ', '(n,dest)  ',&
     '(n,x)     ', '(n,x)     ', '(n,x)     ', '(n,x)     ', '(n,x)     ',&
     '(n,x)     ', '(n,x)     ', '(n,x)     ', '(n,x)     ', '(n,x)     ',&
     '(n,x)     ', '(n,x)     ', '(n,x)     ', '(n,x)     ', '(n,x)     ',&
     '(n,x)     ', '(n,x)     ', '(n,x)     ', '(n,x)     ', '(n,x)     ',&
     '(n,x)     ', '(n,x)     ', '(n,x)     ', '(n,x)     ', '(n,x)     ',&
     '(n,x)     ', '(n,x)     ', '(n,x)     ', '(n,x)     ', '(n,p_0)   ',&
     '(n,p_1)   ', '(n,p_2)   ', '(n,p_3)   ', '(n,p_4)   ', '(n,p_5)   ',&
     '(n,p_6)   ', '(n,p_7)   ', '(n,p_8)   ', '(n,p_9)   ', '(n,p_10)  ',&
     '(n,p_11)  ', '(n,p_12)  ', '(n,p_13)  ', '(n,p_14)  ', '(n,p_15)  ',&
     '(n,p_16)  ', '(n,p_17)  ', '(n,p_18)  ', '(n,p_19)  ', '(n,p_20)  ',&
     '(n,p_21)  ', '(n,p_22)  ', '(n,p_23)  ', '(n,p_24)  ', '(n,p_25)  ',&
     '(n,p_26)  ', '(n,p_27)  ', '(n,p_28)  ', '(n,p_29)  ', '(n,p_30)  ',&
     '(n,p_31)  ', '(n,p_32)  ', '(n,p_33)  ', '(n,p_34)  ', '(n,p_35)  ',&
     '(n,p_36)  ', '(n,p_37)  ', '(n,p_38)  ', '(n,p_39)  ', '(n,p_40)  ',&
     '(n,p_41)  ', '(n,p_42)  ', '(n,p_43)  ', '(n,p_44)  ', '(n,p_45)  ',&
     '(n,p_46)  ', '(n,p_47)  ', '(n,p_48)  ', '(n,p_c)   ', '(n,d_0)   ',&
     '(n,d_1)   ', '(n,d_2)   ', '(n,d_3)   ', '(n,d_4)   ', '(n,d_5)   ',&
     '(n,d_6)   ', '(n,d_7)   ', '(n,d_8)   ', '(n,d_9)   ', '(n,d_10)  ',&
     '(n,d_11)  ', '(n,d_12)  ', '(n,d_13)  ', '(n,d_14)  ', '(n,d_15)  ',&
     '(n,d_16)  ', '(n,d_17)  ', '(n,d_18)  ', '(n,d_19)  ', '(n,d_20)  ',&
     '(n,d_21)  ', '(n,d_22)  ', '(n,d_23)  ', '(n,d_24)  ', '(n,d_25)  ',&
     '(n,d_26)  ', '(n,d_27)  ', '(n,d_28)  ', '(n,d_29)  ', '(n,d_30)  ',&
     '(n,d_31)  ', '(n,d_32)  ', '(n,d_33)  ', '(n,d_34)  ', '(n,d_35)  ',&
     '(n,d_36)  ', '(n,d_37)  ', '(n,d_38)  ', '(n,d_39)  ', '(n,d_40)  ',&
     '(n,d_41)  ', '(n,d_42)  ', '(n,d_43)  ', '(n,d_44)  ', '(n,d_45)  ',&
     '(n,d_46)  ', '(n,d_47)  ', '(n,d_48)  ', '(n,d_c)   ', '(n,t_0)   ',&
     '(n,t_1)   ', '(n,t_2)   ', '(n,t_3)   ', '(n,t_4)   ', '(n,t_5)   ',&
     '(n,t_6)   ', '(n,t_7)   ', '(n,t_8)   ', '(n,t_9)   ', '(n,t_10)  ',&
     '(n,t_11)  ', '(n,t_12)  ', '(n,t_13)  ', '(n,t_14)  ', '(n,t_15)  ',&
     '(n,t_16)  ', '(n,t_17)  ', '(n,t_18)  ', '(n,t_19)  ', '(n,t_20)  ',&
     '(n,t_21)  ', '(n,t_22)  ', '(n,t_23)  ', '(n,t_24)  ', '(n,t_25)  ',&
     '(n,t_26)  ', '(n,t_27)  ', '(n,t_28)  ', '(n,t_29)  ', '(n,t_30)  ',&
     '(n,t_31)  ', '(n,t_32)  ', '(n,t_33)  ', '(n,t_34)  ', '(n,t_35)  ',&
     '(n,t_36)  ', '(n,t_37)  ', '(n,t_38)  ', '(n,t_39)  ', '(n,t_40)  ',&
     '(n,t_41)  ', '(n,t_42)  ', '(n,t_43)  ', '(n,t_44)  ', '(n,t_45)  ',&
     '(n,t_46)  ', '(n,t_47)  ', '(n,t_48)  ', '(n,t_c)   ', '(n,he3_0) ',&
     '(n,he3_1) ', '(n,he3_2) ', '(n,he3_3) ', '(n,he3_4) ', '(n,he3_5) ',&
     '(n,he3_6) ', '(n,he3_7) ', '(n,he3_8) ', '(n,he3_9) ', '(n,he3_10)',&
     '(n,he3_11)', '(n,he3_12)', '(n,he3_13)', '(n,he3_14)', '(n,he3_15)',&
     '(n,he3_16)', '(n,he3_17)', '(n,he3_18)', '(n,he3_19)', '(n,he3_20)',&
     '(n,he3_21)', '(n,he3_22)', '(n,he3_23)', '(n,he3_24)', '(n,he3_25)',&
     '(n,he3_26)', '(n,he3_27)', '(n,he3_28)', '(n,he3_29)', '(n,he3_30)',&
     '(n,he3_31)', '(n,he3_32)', '(n,he3_33)', '(n,he3_34)', '(n,he3_35)',&
     '(n,he3_36)', '(n,he3_37)', '(n,he3_38)', '(n,he3_39)', '(n,he3_40)',&
     '(n,he3_41)', '(n,he3_42)', '(n,he3_43)', '(n,he3_44)', '(n,he3_45)',&
     '(n,he3_46)', '(n,he3_47)', '(n,he3_48)', '(n,he3_c) ', '(n,a_0)   ',&
     '(n,a_1)   ', '(n,a_2)   ', '(n,a_3)   ', '(n,a_4)   ', '(n,a_5)   ',&
     '(n,a_6)   ', '(n,a_7)   ', '(n,a_8)   ', '(n,a_9)   ', '(n,a_10)  ',&
     '(n,a_11)  ', '(n,a_12)  ', '(n,a_13)  ', '(n,a_14)  ', '(n,a_15)  ',&
     '(n,a_16)  ', '(n,a_17)  ', '(n,a_18)  ', '(n,a_19)  ', '(n,a_20)  ',&
     '(n,a_21)  ', '(n,a_22)  ', '(n,a_23)  ', '(n,a_24)  ', '(n,a_25)  ',&
     '(n,a_26)  ', '(n,a_27)  ', '(n,a_28)  ', '(n,a_29)  ', '(n,a_30)  ',&
     '(n,a_31)  ', '(n,a_32)  ', '(n,a_33)  ', '(n,a_34)  ', '(n,a_35)  ',&
     '(n,a_36)  ', '(n,a_37)  ', '(n,a_38)  ', '(n,a_39)  ', '(n,a_40)  ',&
     '(n,a_41)  ', '(n,a_42)  ', '(n,a_43)  ', '(n,a_44)  ', '(n,a_45)  ',&
     '(n,a_46)  ', '(n,a_47)  ', '(n,a_48)  ', '(n,a_c)   ',&
     'free gas  ', 'H(H2O)    ', 'poly      ', '          ',&
     'H(ZrH) inc', 'H(ZrH) coh', 'benzine   ', 'D(D2O)    ',&
     'graph inc ', 'graph coh ', 'Be inc    ', 'Be coh    ',&
     'BeO inc   ', 'BeO coh   ', 'Zr(ZrH)inc', 'Zr(ZrH)coh'/)
   character(10)::h719='(n,p_c)x  '
   character(10)::h739='(n,d_c)x  '
   character(10)::h759='(n,t_c)x  '
   character(10)::h779='(n,he3_c)x'
   character(10)::h799='(n,a_c)x  '
   character(10),dimension(7)::hndf9=(/&
     '(n,xn)    ', '(n,xg)    ', '(n,xp)    ', '(n,xd)    ',&
     '(n,xt)    ', '(n,xhe3)  ','(n,xa)    '/)
   character(10)::h301='heating   '
   character(10)::h443='kerma     '
   character(10)::h444='damage    '
   character(10)::h251='mubar     '
   character(10)::h252='xi        '
   character(10)::h253='gamma     '

   if (iverf.ge.6) then
      i=mt
      if (i.ge.201.and.i.le.207) i=i+200
      if (i.ge.600) i=i-450
      name=hndf(i)
      if (mt.ge.221.and.mt.le.236) name=hndf(408+mt-221)
      if (mt.eq.251) name=h251
      if (mt.eq.252) name=h252
      if (mt.eq.253) name=h253
      if (mt.eq.301) name=h301
      if (mt.eq.443) name=h443
      if (mt.eq.444) name=h444
   else
      if (mt.lt.150) then
         name=hndf(mt)
      else if (mt.ge.201.and.mt.le.207) then
         name=hndf(mt+200)
      else if (mt.ge.221.and.mt.le.236) then
         name=hndf(408+mt-221)
      else if (mt.ge.700.and.mt.lt.718) then
         name=hndf(mt-550)
      else if (mt.eq.718) then
         name=hndf(199)
      else if (mt.eq.719) then
         name=h719
      else if (mt.ge.720.and.mt.lt.738) then
         name=hndf(mt-520)
      else if (mt.eq.738) then
         name=hndf(249)
      else if (mt.eq.739) then
         name=h739
      else if (mt.ge.740.and.mt.lt.758) then
         name=hndf(mt-490)
      else if (mt.eq.758) then
         name=hndf(299)
      else if (mt.eq.759) then
         name=h759
      else if (mt.ge.760.and.mt.lt.779) then
         name=hndf(mt-460)
      else if (mt.eq.778) then
         name=hndf(349)
      else if (mt.eq.779) then
         name=h779
      else if (mt.ge.780.and.mt.lt.798) then
         name=hndf(mt-430)
      else if (mt.eq.798) then
         name=hndf(399)
      else if (mt.eq.799) then
         name=h799
      else if (mt.eq.301) then
         name=h301
      else if (mt.eq.443) then
         name=h444
      else if (mt.eq.444) then
         name=h443
      else if (mt.eq.251) then
         name=h251
      else if (mt.eq.252) then
         name=h252
      else if (mt.eq.253) then
         name=h253
      else
         name='unknown   '
      endif
   endif
   return
   end subroutine rname

end module plotm


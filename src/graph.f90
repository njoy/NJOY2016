module graph
   ! provides basic 2-D and 3-D graphing routines for Postscript
   use locale
   implicit none
   private

   !--Public routines
   public initp,window,endp,endw,init2,frame2,endfr,init3,vect3
   public axis2,axis3,xscale,yscale,xinvrs,yinvrs
   public trans3,grid2,grid3,curv2,curv3,poly2,text3,txtlen,dsym
   public gplot,gdone

   !--Private global variables
   integer::ipage,nps
   real(kr)::width,ht
   integer::ifont
   real(kr)::xpaper,ypaper
   integer::ifg,ibg
   integer::land
   !   wt     weight of tick marks
   !   wg     weight of grid lines
   !   wa     weight of axis lines
   !   tic    length of tic marks
   !   gap    space between axis and numbers
   !   hf     height fraction of fonts numbers
   !   hl     height of labels
   !   hn     height of numbers
   !   lfont  font for numbers and labels
   real(kr)::wt,wg,wa,tic,gap,hf,hl,hn
   integer::lfont
   real(kr)::fvx,dvx,fvy,dvy,fvz,dvz
   integer::logx,logy,logz
   real(kr)::xwll,ywll,www,wwh,wwr,tspace
   real(kr)::rs,ro,rp,ct,st,cp,sp,du,dv

   real(kr)::ushift,vshift,uwidth
   real(kr)::xmin,xmax,xstp,ymin,ymax,ystp,zmin,zmax,zstp

   ! font names
   character(25),dimension(3)::font=(/&
     'Times-Roman              ',&
     'Helvetica                ',&
     'Symbol                   '/)

   ! tables of font sizes -------------------------------------------

   real(kr),dimension(128,3)::cw

   ! times-roman 00-1f
   real(kr),dimension(32)::cw1=(/&
     0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,&
     0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,&
     0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,&
     0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,&
     0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr/)
   ! times-roman 20-3f
   real(kr),dimension(32)::cw2=(/&
     .250e0_kr,.333e0_kr,.408e0_kr,.500e0_kr,.500e0_kr,.833e0_kr,&
     .778e0_kr,.250e0_kr,.333e0_kr,.333e0_kr,.500e0_kr,.564e0_kr,&
     .250e0_kr,.333e0_kr,.250e0_kr,.278e0_kr,.500e0_kr,.500e0_kr,&
     .500e0_kr,.500e0_kr,.500e0_kr,.500e0_kr,.500e0_kr,.500e0_kr,&
     .500e0_kr,.500e0_kr,.278e0_kr,.278e0_kr,.564e0_kr,.564e0_kr,&
     .564e0_kr,.444e0_kr/)
   ! times-roman 40-5f
   real(kr),dimension(32)::cw3=(/&
     .921e0_kr,.722e0_kr,.667e0_kr,.667e0_kr,.722e0_kr,.611e0_kr,&
     .556e0_kr,.722e0_kr,.722e0_kr,.333e0_kr,.389e0_kr,.722e0_kr,&
     .611e0_kr,.889e0_kr,.722e0_kr,.722e0_kr,.556e0_kr,.722e0_kr,&
     .667e0_kr,.556e0_kr,.611e0_kr,.722e0_kr,.722e0_kr,.944e0_kr,&
     .722e0_kr,.722e0_kr,.611e0_kr,.333e0_kr,.278e0_kr,.333e0_kr,&
     .469e0_kr,.500e0_kr/)
   ! times-roman 60-7f
   real(kr),dimension(32)::cw4=(/&
     .333e0_kr,.444e0_kr,.500e0_kr,.444e0_kr,.500e0_kr,.444e0_kr,&
     .333e0_kr,.500e0_kr,.500e0_kr,.278e0_kr,.278e0_kr,.500e0_kr,&
     .278e0_kr,.778e0_kr,.500e0_kr,.500e0_kr,.500e0_kr,.500e0_kr,&
     .333e0_kr,.389e0_kr,.278e0_kr,.500e0_kr,.500e0_kr,.722e0_kr,&
     .500e0_kr,.500e0_kr,.444e0_kr,.480e0_kr,.200e0_kr,.480e0_kr,&
     .000e0_kr,.000e0_kr/)
   ! helvetica 00-1f
   real(kr),dimension(32)::cw5=(/&
     0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,&
     0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,&
     0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,&
     0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,&
     0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr/)
   ! helvetica 20-3f
   real(kr),dimension(32)::cw6=(/&
     .278e0_kr,.278e0_kr,.355e0_kr,.556e0_kr,.556e0_kr,.889e0_kr,&
     .667e0_kr,.250e0_kr,.333e0_kr,.333e0_kr,.389e0_kr,.584e0_kr,&
     .278e0_kr,.333e0_kr,.278e0_kr,.278e0_kr,.556e0_kr,.556e0_kr,&
     .556e0_kr,.556e0_kr,.556e0_kr,.556e0_kr,.556e0_kr,.556e0_kr,&
     .556e0_kr,.556e0_kr,.278e0_kr,.278e0_kr,.584e0_kr,.584e0_kr,&
     .584e0_kr,.556e0_kr/)
   ! helvetica 40-5f
   real(kr),dimension(32)::cw7=(/&
     1.015e0_kr,.667e0_kr,.667e0_kr,.722e0_kr,.722e0_kr,.667e0_kr,&
     .611e0_kr,.778e0_kr,.722e0_kr,.278e0_kr,.500e0_kr,.667e0_kr,&
     .556e0_kr,.833e0_kr,.722e0_kr,.778e0_kr,.667e0_kr,.778e0_kr,&
     .722e0_kr,.667e0_kr,.611e0_kr,.722e0_kr,.667e0_kr,.944e0_kr,&
     .667e0_kr,.667e0_kr,.611e0_kr,.278e0_kr,.278e0_kr,.278e0_kr,&
     .469e0_kr,.556e0_kr/)
   ! helvetica 60-7f
   real(kr),dimension(32)::cw8=(/&
     .222e0_kr,.556e0_kr,.556e0_kr,.500e0_kr,.556e0_kr,.556e0_kr,&
     .278e0_kr,.556e0_kr,.556e0_kr,.222e0_kr,.222e0_kr,.500e0_kr,&
     .222e0_kr,.833e0_kr,.556e0_kr,.556e0_kr,.556e0_kr,.556e0_kr,&
     .333e0_kr,.500e0_kr,.278e0_kr,.556e0_kr,.500e0_kr,.722e0_kr,&
     .500e0_kr,.500e0_kr,.500e0_kr,.334e0_kr,.260e0_kr,.334e0_kr,&
     .000e0_kr,.000e0_kr/)
   ! symbol 00-1f
   real(kr),dimension(32)::cw9=(/&
     0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,&
     0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,&
     0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,&
     0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr,&
     0.e0_kr,0.e0_kr,0.e0_kr,0.e0_kr/)
   ! symbol 20-3f
   real(kr),dimension(32)::cw10=(/&
     .250e0_kr,.333e0_kr,.713e0_kr,.500e0_kr,.549e0_kr,.833e0_kr,&
     .778e0_kr,.250e0_kr,.333e0_kr,.333e0_kr,.500e0_kr,.549e0_kr,&
     .250e0_kr,.549e0_kr,.250e0_kr,.278e0_kr,.500e0_kr,.500e0_kr,&
     .500e0_kr,.500e0_kr,.500e0_kr,.500e0_kr,.500e0_kr,.500e0_kr,&
     .500e0_kr,.500e0_kr,.278e0_kr,.278e0_kr,.549e0_kr,.549e0_kr,&
     .549e0_kr,.444e0_kr/)
   ! symbol 40-5f
   real(kr),dimension(32)::cw11=(/&
     .549e0_kr,.696e0_kr,.660e0_kr,.710e0_kr,.612e0_kr,.652e0_kr,&
     .763e0_kr,.603e0_kr,.765e0_kr,.351e0_kr,.631e0_kr,.724e0_kr,&
     .686e0_kr,.918e0_kr,.739e0_kr,.750e0_kr,.768e0_kr,.741e0_kr,&
     .580e0_kr,.592e0_kr,.632e0_kr,.690e0_kr,.439e0_kr,.768e0_kr,&
     .645e0_kr,.795e0_kr,.650e0_kr,.333e0_kr,.863e0_kr,.333e0_kr,&
     .658e0_kr,.500e0_kr/)
   ! symbol 60-7f
   real(kr),dimension(32)::cw12=(/&
     .500e0_kr,.631e0_kr,.549e0_kr,.549e0_kr,.494e0_kr,.439e0_kr,&
     .521e0_kr,.411e0_kr,.603e0_kr,.329e0_kr,.603e0_kr,.549e0_kr,&
     .549e0_kr,.576e0_kr,.521e0_kr,.549e0_kr,.549e0_kr,.521e0_kr,&
     .549e0_kr,.603e0_kr,.439e0_kr,.576e0_kr,.713e0_kr,.686e0_kr,&
     .493e0_kr,.686e0_kr,.494e0_kr,.480e0_kr,.200e0_kr,.480e0_kr,&
     .000e0_kr,.000e0_kr/)

   ! color tables -------------------------------------------------

   ! light colors for backgrounds
   integer,dimension(3,8)::ibrgb=reshape((/&
     255,255,255,& ! white
     255,222,173,& ! navajo white
     255,235,205,& ! blanched almond
     250,235,215,& ! antique white
     255,255,198,& ! very pale yellow
     255,197,220,& ! very pale rose
     205,250,205,& ! very pale green
     172,233,250&  ! very pale blue
     /),(/3,8/))
   ! dark colors for foregrounds (curves)
   integer,dimension(3,9)::ifrgb=reshape((/&
       0,  0,  0,& ! black
     225,  0,  0,& ! red
       0,200,  0,& ! green
       0,  0,225,& ! blue
     225,  0,225,& ! magenta
       0,225,225,& ! cyan
     170,102, 35,& ! brown
     160, 32,240,& ! purple
     225, 80, 20&  ! orange
     /),(/3,9/))
   integer,dimension(3,40)::isrgb=reshape((/&
     205,255,205,& ! progressive shades of green
     175,235,175,&
     135,225,135,&
     110,210,110,&
      90,180, 90,&
      80,160, 80,&
      70,140, 70,&
      64,125, 64,&
      55,107, 55,&
      45, 90, 45,&
     255,206,206,& ! progressive shades of red
     235,195,195,&
     230,118,118,&
     215, 90, 90,&
     206, 60, 60,&
     200, 40, 30,&
     188, 39, 20,&
     175, 32, 32,&
     160,  0,  0,&
     135,  0,  0,&
     255,218,177,& ! progressive shades of brown
     243,195,142,&
     237,179,108,&
     225,159, 75,&
     215,140, 76,&
     200,131, 62,&
     175,117, 52,&
     160, 90, 33,&
     140, 72, 29,&
     121, 62, 25,&
     192,237,253,& ! progressive shades of blue
     135,164,229,&
     120,120,220,&
     100,100,210,&
      70, 70,200,&
      40, 40,195,&
      21, 21,182,&
      20, 20,170,&
      16, 16,140,&
       0,  0,125/),(/3,40/))

contains

   subroutine initp(iori,xpage,ypage,istyle,htt,wtt,ibord,ipcol)
   !--------------------------------------------------------------------
   ! Initialize the page.
   ! Sets up the the basic page size, orientation, font style,
   ! font height, normal line weight, and page background color.
   ! Draws page border if requested.
   ! Page coordinates are measured from the
   ! lower-left corner of the page in inches.
   !--------------------------------------------------------------------
   ! externals
   integer::iori,istyle,ibord,ipcol
   real(kr)::xpage,ypage,htt,wtt
   ! internals
   integer::n
   real(kr)::wbord,w,backgr
   real(kr)::x(5),y(5)

   ! load font tables
   do n=1,32
      cw(n,1)=cw1(n)
      cw(32+n,1)=cw2(n)
      cw(64+n,1)=cw3(n)
      cw(96+n,1)=cw4(n)
      cw(n,2)=cw5(n)
      cw(32+n,2)=cw6(n)
      cw(64+n,2)=cw7(n)
      cw(96+n,2)=cw8(n)
      cw(n,3)=cw9(n)
      cw(32+n,3)=cw10(n)
      cw(64+n,3)=cw11(n)
      cw(96+n,3)=cw12(n)
   enddo

   ! width of border line
   wbord=.005e0

   land=iori
   ifont=istyle
   width=wtt
   ht=htt
   ibg=1+ipcol
   ifg=1
   ro=1
   rs=1
   ct=0
   st=1
   cp=0
   sp=-1
   du=0
   dv=0
   xwll=0
   ywll=0
   www=0
   wwh=0
   wwr=0
   call newp
   n=5
   x(1)=0
   y(1)=0
   x(2)=xpage
   y(2)=0
   x(3)=xpage
   y(3)=ypage
   x(4)=0
   y(4)=ypage
   x(5)=x(1)
   y(5)=y(1)
   w=wbord
   if (ibg.eq.1) w=-1
   backgr=1
   call poly2(x,y,n,w,backgr)
   return
   end subroutine initp

   subroutine window(xll,yll,ww,wh,wr,t1,n1,t2,n2,ibord)
   !--------------------------------------------------------------------
   ! Set up a window on the page
   ! and write the titles (if any) in the upper left corner.
   ! Draw a window border if requested.
   ! Set up a clipping path at the window border.
   !--------------------------------------------------------------------
   ! externals
   integer::n1,n2,ibord
   real(kr)::xll,yll,ww,wh,wr
   character::t1*(*),t2*(*)
   ! internals
   real(kr)::wbord,w,ull,vll,uur,vur,x,y
   real(kr),parameter::zero=0
   real(kr),parameter::one=1
   ifg=1

   ! width of window border
   wbord=.005e0_kr

   !--load window transformation
   xwll=xll
   ywll=yll
   www=ww
   wwh=wh
   wwr=wr

   !--draw window border
   if (ibord.ne.0) then
      w=wbord
      call vect2(zero,zero,zero,wwh,w,0)
      call vect2(zero,wwh,www,wwh,w,0)
      call vect2(www,wwh,www,zero,w,0)
      call vect2(www,zero,zero,zero,w,0)
   endif

   !--set clipping to window border
   call transw(zero,zero,ull,vll)
   call transw(www,wwh,uur,vur)
   call gset(ull,vll,uur,vur)

   !--draw titles, if any.
   tspace=0
   if (n1.ne.0) then
      x=3*ht
      y=wwh-3*ht/2
      tspace=tspace+3*ht/2
      call text2(t1,n1,ifont,ht,x,y,one,zero,zero,one)
      if (n2.ne.0) then
         y=y-(ht+ht/5)
         tspace=tspace+3*ht/2
         call text2(t2,n2,ifont,ht,x,y,one,zero,zero,one)
      endif
   endif
   return
   end subroutine window

   subroutine endw
   !--------------------------------------------------------------------
   ! End a window by resetting the clipping path
   ! and restoring the default transformation.
   !--------------------------------------------------------------------
   call gend
   ro=1
   rs=1
   ct=0
   st=1
   cp=0
   sp=-1
   du=0
   dv=0
   www=0
   wwh=0
   wwr=0
   return
   end subroutine endw

   subroutine transw(x,y,u,v)
   !--------------------------------------------------------------------
   ! Coordinate transformation for the window.
   !--------------------------------------------------------------------
   use physics ! provides pi
   ! externals
   real(kr)::x,y,u,v
   ! internals
   real(kr)::ct,st
   real(kr),parameter::zero=0

   if (www.eq.zero) then
      u=x
      v=y
   else
      ct=cos(2*pi*wwr/360)
      st=sin(2*pi*wwr/360)
      u=xwll+x*ct-y*st
      v=ywll+x*st+y*ct
   endif
   return
   end subroutine transw

   subroutine init2(uo,vo,xg,yg,iright,iwcol)
   !--------------------------------------------------------------------
   ! Initialize a 2-d plot by choosing an origen and axis
   ! lengths to nicely position the graph in the window.  Here
   !  u0,vo is the position of the lower-left corner of the
   !   axis system in window coordinates,
   !  xg,yg are the lengths of the axis is absolute units,
   !  iright>0 reserves space for a right-hand axis, and
   !  iwcol is the background color of the inside of the graph frame.
   !--------------------------------------------------------------------
   ! externals
   integer::iright,iwcol
   real(kr)::uo,vo,xg,yg
   ! internals
   integer::n
   real(kr)::wframe,ull,vll,uur,vur,w,backgr
   real(kr)::x(5),y(5)
   real(kr),parameter::zero=0

   ! width of frame border
   wframe=.005e0_kr

   !--adjust axis lengths to fit graph inside window
   uo=4*ht
   vo=4*ht
   xg=www-uo-2*ht
   if (iright.gt.0) xg=xg-2*ht
   yg=wwh-tspace-vo-ht
   if (tspace.eq.zero) yg=yg-ht

   !--set up axis parameters
   wt=width
   wg=wt/2
   wa=2*wt
   tic=ht/2
   gap=4*ht/10
   lfont=ifont
   hf=3
   hf=hf/4
   hl=ht
   hn=8*hl/10

   !--set up coordinate transformation (no perspective)
   ro=1
   rs=1
   ct=0
   st=1
   cp=0
   sp=-1
   du=uo
   dv=vo

   !--color in the background inside the graph frame.
   ibg=1+iwcol
   ifg=1
   call transw(uo,vo,ull,vll)
   call transw(uo+xg,vo+yg,uur,vur)
   n=5
   x(1)=ull
   y(1)=vll
   x(2)=uur
   y(2)=vll
   x(3)=uur
   y(3)=vur
   x(4)=ull
   y(4)=vur
   x(5)=x(1)
   y(5)=y(1)
   w=wframe
   backgr=1
   call poly2(x,y,n,w,backgr)
   return
   end subroutine init2

   subroutine frame2(xg,yg,grace,iccol)
   !--------------------------------------------------------------------
   ! Sets up a clipping path just outside the frame boundary
   ! and sets foreground color for this curve.
   !--------------------------------------------------------------------
   ! externals
   integer::iccol
   real(kr)::xg,yg,grace,xll,yll,xur,yur,ull,vll,uur,vur
   ! internals
   real(kr),parameter::zero=0

   ifg=1+iccol
   call trans3(-grace,-grace,zero,xll,yll)
   call trans3(xg+grace,yg+grace,zero,xur,yur)
   call transw(xll,yll,ull,vll)
   call transw(xur,yur,uur,vur)
   call gset(ull,vll,uur,vur)
   return
   end subroutine frame2

   subroutine endfr
   !--------------------------------------------------------------------
   ! End the 2D plotting frame by clearing the frame clipping path.
   !--------------------------------------------------------------------
   call gend
   end subroutine endfr

    subroutine init3(bx,by,bz,vx,vy,vz,iwcol)
   !--------------------------------------------------------------------
   ! Initialize a 3-d plot.
   !--------------------------------------------------------------------
   ! externals
   integer::iwcol
   real(kr)::bx,by,bz,vx,vy,vz
   ! internals
   real(kr)::umin,umax,vmin,vmax,yy,u,v,xx,dx,dy,su,sv,s
   real(kr),parameter::zero=0

   !--set color of 3d slices to the window color
   ibg=1+iwcol
   ifg=1

   !--set up axis parameters
   wt=width
   wg=wt/2
   wa=wt
   tic=ht/2
   gap=4*ht/10
   lfont=ifont
   hf=3
   hf=hf/4
   hl=ht
   hn=ht

   !--set up 3-d transformation
   ro=sqrt(vx**2+vy**2+vz**2)
   rp=sqrt(vx**2+vy**2)
   ct=rp/ro
   st=vz/ro
   cp=vx/rp
   sp=vy/rp
   rs=ro/sqrt(bx**2+by**2+bz**2)
   du=0
   dv=0

   !--adjust view scale and position to
   !--center the projection in the window
   umin=1000
   umax=-1000
   vmin=1000
   vmax=-1000
   yy=-5*ht
   if (bx.gt.zero) call trans3(zero,yy,zero,u,v)
   if (bx.lt.zero) call trans3(bx,yy,zero,u,v)
   if (u.lt.umin) umin=u
   if (u.gt.umax) umax=u
   if (v.lt.vmin) vmin=v
   if (v.gt.vmax) vmax=v
   if (bx.gt.zero) call trans3(zero,by,bz+ht,u,v)
   if (bx.lt.zero) call trans3(bx,by,bz+ht,u,v)
   if (u.lt.umin) umin=u
   if (u.gt.umax) umax=u
   if (v.lt.vmin) vmin=v
   if (v.gt.vmax) vmax=v
   xx=bx+2*ht
   yy=-7*ht/2
   if (bx.gt.zero) call trans3(xx,yy,zero,u,v)
   if (bx.lt.zero) call trans3(zero,yy,zero,u,v)
   if (u.lt.umin) umin=u
   if (u.gt.umax) umax=u
   if (v.lt.vmin) vmin=v
   if (v.gt.vmax) vmax=v
   dx=3*ht
   dy=2*ht
   if (bx.gt.zero) call trans3(bx+dx,by+dy,zero,u,v)
   if (bx.lt.zero) call trans3(dx,by+dy,zero,u,v)
   if (u.lt.umin) umin=u
   if (u.gt.umax) umax=u
   if (v.lt.vmin) vmin=v
   if (v.gt.vmax) vmax=v
   su=www/(umax-umin)
   sv=wwh/(vmax-vmin)
   s=su
   if (sv.lt.su) s=sv
   rs=s*rs
   du=(www-s*(umax-umin))/2-s*umin
   ! dv=(wwh-tspace-s*(vmax-vmin))/2-s*vmin
   dv=(wwh-s*(vmax-vmin))/2-s*vmin
   return
   end subroutine init3

   subroutine trans3(x,y,z,u,v)
   !--------------------------------------------------------------------
   ! Carry out the 3-d perspective transformation.
   !--------------------------------------------------------------------
   ! externals
   real(kr)::x,y,z,u,v,xo,s

   xo=ct*(cp*x+sp*y)+st*z
   s=rs/(ro-xo)
   u=s*(-sp*x+cp*y)+du
   v=s*(-st*(cp*x+sp*y)+ct*z)+dv
   return
   end subroutine trans3

   subroutine vect2(x1,y1,x2,y2,w,ivec)
   !--------------------------------------------------------------------
   ! Draw a 2d vector using the 3d vector call.
   !--------------------------------------------------------------------
   ! externals
   integer::ivec
   real(kr)::x1,y1,x2,y2,w
   ! internals
   real(kr),parameter::zero=0

   call vect3(x1,y1,zero,x2,y2,zero,w,ivec)
   return
   end subroutine vect2

   subroutine vect3(x1,y1,z1,x2,y2,z2,w,ivec)
   !--------------------------------------------------------------------
   ! Draw a 3d vector, where
   !  w is the thickness of the line.
   ! Set ivec.ne.0 to put an arrow head on the end.
   !--------------------------------------------------------------------
   ! externals
   integer::ivec
   real(kr)::x1,y1,z1,x2,y2,z2,w
   ! internals
   real(kr)::wu1,wv1,wu2,wv2,u1,u2,head,r,wx,wy,dx,dy,hx,hy,v1,v2
   ! size of vector arrowhead
   real(kr),parameter::ahead=.2e0_kr

   call trans3(x1,y1,z1,wu1,wv1)
   call trans3(x2,y2,z2,wu2,wv2)
   call transw(wu1,wv1,u1,v1)
   call transw(wu2,wv2,u2,v2)
   if (ivec.eq.0) then
      call moveh(u1,v1)
      call drawh(u2,v2,w,0)
   else
      head=ahead
      r=sqrt((u2-u1)**2+(v2-v1)**2)
      wx=w*(u2-u1)/r
      wy=w*(v2-v1)/r
      dx=head*(u2-u1)/r
      dy=head*(v2-v1)/r
      hx=dx/2
      hy=dy/2
      call moveh(u1,v1)
      call drawh(u2-wx,v2-wy,w,0)
      call moveh(u2-dx-hy,v2-dy+hx)
      call drawh(u2,v2,w,0)
      call drawh(u2-dx+hy,v2-dy-hx,w,0)
   endif
   return
   end subroutine vect3

   subroutine axis2(amin,amax,astp,label,nlabel,itic,nonum,&
     asize,xo,yo,xx,yx,xy,yy,numr)
   !--------------------------------------------------------------------
   ! Draw a 2-D axis using the 3-D axis routine.
   !--------------------------------------------------------------------
   ! externals
   integer::nlabel,itic,nonum,numr
   real(kr)::amin,amax,astp,asize,xo,yo,xx,yx,xy,yy
   character::label*(*)
   ! internals
   real(kr),parameter::zero=0

   call axis3(amin,amax,astp,label,nlabel,itic,nonum,asize,&
     xo,yo,zero,xx,yx,zero,xy,yy,zero,numr)
   return
   end subroutine axis2

   subroutine axis3(amin,amax,astp,label,nlabel,itic,nonum,&
     asize,xo,yo,zo,xx,yx,zx,xy,yy,zy,numr)
   !--------------------------------------------------------------------
   ! Draw an axis in 3-D perspective.
   ! ------
   ! amin is the left limit for the axis
   ! amax is the right limit for the axis
   ! astp is the step between axis numbering values
   !   (if astp=0., this is a log axis)
   ! nlabel=0 means no label or numbers
   ! nlabel>0 means label and numbers "above" axis
   ! nlabel<0 means label and number "below" axis
   ! itic=0 means no tick marks
   ! itic>0 means tick marks above axis
   ! itic<0 means tick marks below axis
   ! nonum=0 means all tick marks are numbered
   ! nonum=1 means all but first
   ! nonum=2 means all but last
   ! nonum=3 means neither end is numbered
   ! asize is the length of the axis in absolute units
   ! xo,yo,zo is the absolute location of the start of the axis
   ! xx,yx,zx is a unit vector in the direction of the axis
   ! xy,yy,zy is a unit vector for the character vertical direction
   !--------------------------------------------------------------------
   ! externals
   integer::nlabel,itic,nonum,numr
   real(kr)::amin,amax,astp,asize,xo,yo,zo,xx,yx,zx,xy,yy,zy
   character::label*(*)
   ! internals
   integer::n,nscale,iv,ifracs,i,imin,j,imax,iskip
   integer::k1,k,idone,k2,iii,lnum
   real(kr)::xe,ye,ze,shl,step,test,x,y,z,v,scale,vv
   real(kr)::ww,wwt,sh,xc,yc,zc,www,dd,origen,aend
   real(kr)::cycles,room,xt,yt,zt,xs,ys,zs,dx,dy,dz
   real(kr)::xl,yl,zl
      character num*20
   real(kr),dimension(8)::subs=(/.301e0_kr,.477e0_kr,.602e0_kr,&
     .699e0_kr,.778e0_kr,.845e0_kr,.903e0_kr,.954e0_kr/)
   real(kr),parameter::ds=.01e0_kr
   real(kr),parameter::sc1=.099e0_kr
   real(kr),parameter::sc2=.901e0_kr
   real(kr),parameter::zero=0
   real(kr),parameter::ten=10

   !--set color to black
   ifg=1

   !--draw the axis vector
   xe=xo+asize*xx
   ye=yo+asize*yx
   ze=zo+asize*zx
   call vect3(xo,yo,zo,xe,ye,ze,wa,0)
   shl=0

   !--draw the tic marks and numbers for a linear scale
   if (astp.ne.zero) then
      n=nint((amax-amin)/astp)
      step=asize/n
      n=n+1
      test=1
      test=test/1000
      if (abs(xx-1).lt.test) then
         fvx=amin
         dvx=(amax-amin)/asize
         logx=0
         xmin=xo
         xmax=xo+asize
         xstp=step
      else if (abs(yx-1).lt.test) then
         fvy=amin
         dvy=(amax-amin)/asize
         logy=0
         ymin=yo
         ymax=yo+asize
         ystp=step
      else if (abs(zx-1).lt.test) then
         fvz=amin
         dvz=(amax-amin)/asize
         logz=0
         zmin=zo
         zmax=zo+asize
         zstp=step
      endif
      x=xo
      y=yo
      z=zo
      v=amin
      nscale=0
      if (abs(astp).lt.sc1.or.abs(astp).gt.sc2) then
         nscale=nint(log10(abs(astp)))
         if (nscale.lt.0) nscale=nscale-2
         nscale=3*(nscale/3)
      endif
      scale=ten**nscale
      !test=1
      !test=test-test/100
      !if (abs(amax/scale).lt.test) then
      !   nscale=nscale-3
      !   scale=scale/1000
      !endif
      iv=nint(astp/scale)
      vv=astp/scale
      test=1
      test=test/100
      ifracs=0
      if (abs(vv-iv).gt.test) ifracs=1
      do i=1,n
         if (itic.lt.0) then
            call vect3(x,y,z,x-tic*xy,y-tic*yy,z-tic*zy,wt,0)
         else if (itic.gt.0) then
            call vect3(x,y,z,x+tic*xy,y+tic*yy,z+tic*zy,wt,0)
         endif
         if (nlabel.ne.0) then
            if (i.ne.1.or.(nonum.ne.1.and.nonum.ne.3)) then
               iv=nint(v/scale)
               vv=v/scale
               if (ifracs.eq.1) then
                  write(num,'(f6.1)') vv
               else
                  write(num,'(i5)') iv
               endif
               call stripv(num,lnum)
               ww=txtlen(num,lnum,lfont,hn)/2
               wwt=ww
               if (i.eq.1.and.itic.gt.0.and.nlabel.lt.0) wwt=ww/2
               if (i.eq.1.and.itic.lt.0.and.nlabel.gt.0) wwt=ww/2
               if (i.eq.n.and.itic.gt.0.and.nlabel.lt.0) wwt=3*ww/2
               if (i.eq.n.and.itic.lt.0.and.nlabel.gt.0) wwt=3*ww/2
               if (nlabel.lt.0) then
                  sh=gap+hf*hn
                  if (itic.lt.0) sh=sh+tic-gap/2
                  if (numr.eq.0) then
                     xc=x-sh*xy-wwt*xx
                     yc=y-sh*yy-wwt*yx
                     zc=z-sh*zy-wwt*zx
                     if (sh+hf*hn.gt.shl) shl=sh+hf*hn
                  else
                     xc=x-sh*xy-hf*hn*xx/2
                     yc=y-sh*yy-hf*hn*yx/2
                     zc=z-sh*zy-hf*hn*zx/2
                     test=sh+2*ww+hn/2
                     if (test.gt.shl) shl=test
                  endif
               else
                  sh=gap
                  if (itic.gt.0) sh=sh+tic-gap/2
                  if (numr.eq.0) then
                     xc=x+sh*xy-wwt*xx
                     yc=y+sh*yy-wwt*yx
                     zc=z+sh*zy-wwt*zx
                     test=sh+2*hf*hn
                     if (test.gt.shl) shl=test
                  else
                     xc=x+(sh+2*ww)*xy-hf*hn*xx/2
                     yc=y+(sh+2*ww)*yy-hf*hn*yx/2
                     zc=z+(sh+2*ww)*zy-hf*hn*zx/2
                     test=sh+2*ww+hn/2
                     if (test.gt.shl) shl=test
                  endif
               endif
               if (numr.eq.0) then
                  call text3(num,lnum,lfont,hn,&
                    xc,yc,zc,xx,yx,zx,xy,yy,zy)
               else
                  call text3(num,lnum,lfont,hn,&
                    xc,yc,zc,-xy,-yy,-zy,xx,yx,zx)
               endif
            endif
         endif
         x=x+step*xx
         y=y+step*yx
         z=z+step*zx
         v=v+astp
      enddo
      if (nscale.ne.0.and.nlabel.ne.0) then
         if (nonum.ne.2.and.nonum.ne.3) then
            if ((nscale.gt.-10.and.nscale.lt.0).or.nscale.ge.10) then
               write(num,'(''*10#EH.8<'',i2,''#HXEX<'')') nscale
               lnum=17
            else if (nscale.le.-10) then
               write(num,'(''*10#EH.8<'',i3,''#HXEX<'')') nscale
               lnum=18
            else
               write(num,'(''*10#EH.8<'',i1,''#HXEX<'')') nscale
               lnum=16
            endif
            www=txtlen(num,lnum,lfont,hn)
            if (nlabel.lt.0) then
               if (numr.eq.0) then
                  xc=xc-6*hn*xy/5
                  yc=yc-6*hn*yy/5
                  zc=zc-6*hn*zy/5
               else
                  xc=xc-6*hn*xx/5
                  yc=yc-6*hn*yx/5
                  zc=zc-6*hn*zx/5
               endif
            else
               if (numr.eq.0) then
                  xc=xc+6*hn*xy/5
                  yc=yc+6*hn*yy/5
                  zc=zc+6*hn*zy/5
               else
                  dd=hf*hn/2+6*hn/5
                  xc=x+(www+gap)*xy-dd*xx-step*xx
                  yc=y+(www+gap)*yy-dd*yx-step*yx
                  zc=z+(www+gap)*zy-dd*zx-step*zx
               endif
            endif
            if (numr.eq.0) then
               call text3(num,lnum,lfont,hn,&
                 xc,yc,zc,xx,yx,zx,xy,yy,zy)
            else
               call text3(num,lnum,lfont,hn,&
                 xc,yc,zc,-xy,-yy,-zy,xx,yx,zx)
            endif
         endif
      endif

   !--draw the tick marks and numbers for a log scale
   else
      origen=log10(amin)
      i=int(origen)
      if (origen.lt.zero) i=i-1
      v=origen-i
      if (v.lt.ds) then
         v=0
      else if (v.lt.subs(1)+ds) then
         v=subs(1)
      else if (v.lt.subs(4)+ds) then
         v=subs(4)
      else
         i=i+1
         v=0
      endif
      origen=i+v
      imin=i
      if (v.gt.zero) imin=i+1
      if ((imin.gt.-10.and.imin.lt.0).or.imin.ge.10) then
         write(num,'(''10#EH.8<'',i2,''#HXEX<'')') imin
         lnum=16
      else if (imin.le.-10) then
         write(num,'(''10#EH.8<'',i3,''#HXEX<'')') imin
         lnum=17
      else
         write(num,'(''10#EH.8<'',i1,''#HXEX<'')') imin
         lnum=15
      endif
      ww=txtlen(num,lnum,lfont,hn)
      www=ww
      aend=log10(amax)
      j=int(aend)
      if (aend.lt.zero) j=j-1
      v=aend-j
      if (v.lt.ds) then
         v=0
      else if (v.lt.subs(1)+ds) then
         v=subs(1)
      else if (v.lt.subs(4)+ds) then
         v=subs(4)
      else
         j=j+1
         v=0
      endif
      aend=j+v
      imax=j
      if (v.gt.0.) imax=imax-1
      if ((imax.gt.-10.and.imax.lt.0).or.imax.ge.10) then
         write(num,'(''10#EH.8<'',i2,''#HXEX<'')') imax
         lnum=16
      else if (imax.le.-10) then
         write(num,'(''10#EH.8<'',i3,''#HXEX<'')') imax
         lnum=17
      else
         write(num,'(''10#EH.8<'',i1,''#HXEX<'')') imax
         lnum=15
      endif
      ww=txtlen(num,lnum,lfont,hn)
      if (ww.gt.www) www=ww
      cycles=(aend-origen)/asize
      test=1
      test=test/10
      if (abs(xx*xy+yx*yy+zx*zy-1).lt.test) then
         room=6*hn/5
      else
         room=6*www/5
      endif
      iskip=int(room*abs(cycles))
      iskip=1+iskip
      test=1
      test=test/1000
      if (abs(xx-1).lt.test) then
         fvx=origen
         dvx=cycles
         logx=1
         xmin=xo
         xmax=xo+asize
      else if (abs(yx-1).lt.test) then
         fvy=origen
         dvy=cycles
         logy=1
         ymin=yo
         ymax=yo+asize
      else if (abs(zx-1).lt.test) then
         fvz=origen
         dvz=cycles
         logz=1
         zmin=zo
         zmax=zo+asize
      endif
      if (itic.lt.0) then
         xt=-tic*xy
         yt=-tic*yy
         zt=-tic*zy
      else
         xt=tic*xy
         yt=tic*yy
         zt=tic*zy
      endif
      i=int(origen)
      test=1
      test=test/10
      if (origen-i.gt.test) i=i+1
      v=ten**i
      x=xo+((i-origen)/cycles)*xx
      y=yo+((i-origen)/cycles)*yx
      z=zo+((i-origen)/cycles)*zx
      if (itic.ne.0) then
         test=1
         if (cycles.le.test) then
            k1=0
            test=1
            test=test/10
            if (abs(origen-i+1-subs(1)).lt.test) k1=1
            if (abs(origen-i+1-subs(4)).lt.test) k1=4
            if (k1.ne.0) then
               do k=k1,8
                  sh=(subs(k)-1)/cycles
                  xs=x+sh*xx
                  ys=y+sh*yx
                  zs=z+sh*zx
                  dx=xt/2
                  dy=yt/2
                  dz=zt/2
                  call vect3(xs,ys,zs,xs+dx,ys+dy,zs+dz,wt,0)
               enddo
            endif
         endif
      endif
      idone=0
      do while (idone.eq.0)
         if (itic.ne.0) then
            call vect3(x,y,z,x+xt,y+yt,z+zt,wt,0)
            test=1
            test=test/10
            if (abs(aend-i).ge.test) then
               test=1
               if (cycles.le.test) then
                  k2=8
                  test=1
                  test=test/10
                  if (abs(aend-i-subs(1)).lt.test) k2=1
                  if (abs(aend-i-subs(4)).lt.test) k2=4
                  do k=1,k2
                     sh=subs(k)/cycles
                     xs=x+sh*xx
                     ys=y+sh*yx
                     zs=z+sh*zx
                     dx=xt/2
                     dy=yt/2
                     dz=zt/2
                     call vect3(xs,ys,zs,xs+dx,ys+dy,zs+dz,wt,0)
                  enddo
               endif
            endif
         endif
         if (nlabel.ne.0) then
            test=1
            test=test/1000
            if (abs(i-origen).ge.test.or.&
              (nonum.ne.1.and.nonum.ne.3)) then
               if (log10(v).lt.aend.or.&
                 (nonum.ne.2.and.nonum.ne.3)) then
                  iii=int(i-origen)
                  if (mod(iii,iskip).eq.0) then
                     if ((i.gt.-10.and.i.lt.0).or.i.ge.10) then
                        write(num,'(''10#EH.8<'',i2,''#HXEX<'')') i
                        lnum=16
                     else if (i.le.-10) then
                        write(num,'(''10#EH.8<'',i3,''#HXEX<'')') i
                        lnum=17
                     else
                        write(num,'(''10#EH.8<'',i1,''#HXEX<'')') i
                        lnum=15
                     endif
                     ww=txtlen(num,lnum,lfont,hn)/2
                     if (i.eq.origen) then
                        if (itic.eq.0) then
                           ww=ww/2
                        else
                           if (itic.lt.0.and.nlabel.gt.0) then
                              ww=ww/2
                           else if (itic.gt.0.and.nlabel.lt.0) then
                              ww=ww/2
                           endif
                        endif
                     endif
                     if (nlabel.lt.0) then
                        sh=gap+hf*hn
                        if (itic.lt.0) sh=sh+tic-gap/2
                        if (numr.eq.0) then
                           xc=x-sh*xy-ww*xx
                           yc=y-sh*yy-ww*yx
                           zc=z-sh*zy-ww*zx
                           test=sh+hf*hn
                           if (test.gt.shl) shl=test
                        else
                           xc=x-sh*xy-hf*hn*xx/2
                           yc=y-sh*yy-hf*hn*yx/2
                           zc=z-sh*zy-hf*hn*zx/2
                           test=sh+www+hn/2
                           if (test.gt.shl) shl=test
                        endif
                     else
                        sh=gap
                        if (itic.gt.0) sh=sh+tic-gap/2
                        if (numr.eq.0) then
                           xc=x+sh*xy-ww*xx
                           yc=y+sh*yy-ww*yx
                           zc=z+sh*zy-ww*zx
                           test=sh+2*hf*hn
                           if (test.gt.shl) shl=test
                        else
                           xc=x+(sh+www)*xy-hf*hn*xx/2
                           yc=y+(sh+www)*yy-hf*hn*yx/2
                           zc=z+(sh+www)*zy-hf*hn*zx/2
                           test=sh+www+hn/2
                           if (test.gt.shl) shl=test
                        endif
                     endif
                     if (numr.eq.0) then
                        call text3(num,lnum,lfont,hn,&
                          xc,yc,zc,xx,yx,zx,xy,yy,zy)
                     else
                        call text3(num,lnum,lfont,hn,&
                          xc,yc,zc,-xy,-yy,-zy,xx,yx,zx)
                     endif
                  endif
               endif
            endif
         endif
         if (i.ge.j) then
            idone=1
         else
            x=x+(1/cycles)*xx
            y=y+(1/cycles)*yx
            z=z+(1/cycles)*zx
            i=i+1
            v=v*10
         endif
      enddo
   endif

   !--add the centered axis label
   if (nlabel.ne.0) then
      if (label.ne.'.') then
         ww=txtlen(label,iabs(nlabel),lfont,hl)
         ww=(asize-ww)/2
         if (nlabel.lt.0) then
            sh=shl+hf*hl
            xl=xo-sh*xy+ww*xx
            yl=yo-sh*yy+ww*yx
            zl=zo-sh*zy+ww*zx
         else
            sh=shl
            xl=xo+sh*xy+ww*xx
            yl=yo+sh*yy+ww*yx
            zl=zo+sh*zy+ww*zx
         endif
         call text3(label,iabs(nlabel),lfont,hl,&
           xl,yl,zl,xx,yx,zx,xy,yy,zy)
      endif
   endif
   return
   end subroutine axis3

   real(kr) function xscale(x)
   !--------------------------------------------------------------------
   ! Function to convert user x coordinate to axis coordinate.
   !--------------------------------------------------------------------
   ! externals
   real(kr)::x
   ! internals
   real(kr)::xn

   xn=x
   if (logx.eq.1) xn=log10(xn)
   xscale=(xn-fvx)/dvx
   return
   end function xscale

   real(kr) function yscale(y)
   !--------------------------------------------------------------------
   ! Function to convert user y coordinate to axis coordinate.
   !--------------------------------------------------------------------
   ! externals
   real(kr)::y
   ! internals
   real(kr)::yn

   yn=y
   if (logy.eq.1) yn=log10(yn)
   yscale=(yn-fvy)/dvy
   return
   end function yscale

   real(kr) function zscale(z)
   !--------------------------------------------------------------------
   ! Function to convert user z coordinate to axis coordinate.
   !--------------------------------------------------------------------
   ! externals
   real(kr)::z
   ! internals
   real(kr)::zn

   zn=z
   if (logz.eq.1) zn=log10(zn)
   zscale=(zn-fvz)/dvz
   return
   end function zscale

   real(kr) function xinvrs(x)
   !--------------------------------------------------------------------
   ! Function to convert x axis coordinates back to user coordinates.
   !--------------------------------------------------------------------
   ! externals
   real(kr)::x
   ! internals
   real(kr)::xn
   real(kr),parameter::ten=10

   xn=x*dvx+fvx
   xinvrs=xn
   if (logx.eq.1) xinvrs=ten**xn
   return
   end function xinvrs

   real(kr) function yinvrs(y)
   !--------------------------------------------------------------------
   ! Function to convert y axis coordinates back to user coordinates.
   !--------------------------------------------------------------------
   ! externals
   real(kr)::y
   ! internals
   real(kr)::yn
   real(kr),parameter::ten=10

   yn=y*dvy+fvy
   yinvrs=yn
   if (logy.eq.1) yinvrs=ten**yn
   return
   end function yinvrs

   subroutine grid2(nx,ny)
   !--------------------------------------------------------------------
   ! Draw a 2D grid using the 3D grid call
   !--------------------------------------------------------------------
   ! externals
   integer::nx,ny

   zmin=0
   zmax=0
   call grid3(nx,ny,0)
   return
   end subroutine grid2

   subroutine grid3(nx,ny,nz)
   !--------------------------------------------------------------------
   ! Draw a grid with nx lines per x axis step, ny lines
   ! per y axis step, and nz lines per z axis step.
   ! One of these three values must be zero or -1 to
   ! indicate in which plane the grid is to be drawn. For example,
   ! nx=0 means grid plane is at xmin, and nx=-1 means plane at xmax.
   ! The values of nx, ny, and/or nz are ignored for log axes.
   !--------------------------------------------------------------------
   ! externals
   integer::nx,ny,nz
   ! internals
   integer::i,k,idone
   real(kr)::y1,y2,z1,z2,x1,x2,v,test,xt,yt,zt
   real(kr),dimension(9)::subs=(/.301e0_kr,.477e0_kr,.602e0_kr,&
     .699e0_kr,.778e0_kr,.845e0_kr,.903e0_kr,.954e0_kr,1.000e0_kr/)

   if (nx.gt.0) then
      if (nz.eq.0) then
         y1=ymin
         y2=ymax
         z1=zmin
         z2=zmin
      else if (nz.eq.-1) then
         y1=ymin
         y2=ymax
         z1=zmax
         z2=zmax
      else
         if (ny.eq.0) then
            y1=ymin
            y2=ymin
            z1=zmin
            z2=zmax
         else
            y1=ymax
            y2=ymax
            z1=zmin
            z2=zmax
         endif
      endif
      x1=xmin
      if (logx.eq.1) then
         v=fvx
         i=nint(v)
         k=0
         test=1
         test=test/10
         if (abs(v-i-subs(1)).lt.test) k=1
         if (abs(v-i-subs(4)).lt.test) k=4
         xt=x1
         if (k.gt.0) xt=x1-subs(k)/dvx
         k=k+1
      endif
      idone=0
      do while (idone.eq.0)
         if (logx.eq.0) then
            x1=x1+xstp/nx
         else
            x1=xt+subs(k)/dvx
            k=k+1
            if (k.gt.9) then
               xt=x1
               k=1
            endif
         endif
         test=1
         test=test/100
         if (abs(xmax-x1).lt.test) then
            idone=1
         else
            x2=x1
            call vect3(x1,y1,z1,x2,y2,z2,wg,0)
         endif
      enddo
   endif

   if (ny.gt.0) then
      if (nz.eq.0) then
         x1=xmin
         x2=xmax
         z1=zmin
         z2=zmin
      else if (nz.eq.-1) then
         x1=xmin
         x2=xmax
         z1=zmax
         z2=zmax
      else
         if (nx.eq.0) then
            x1=xmin
            x2=xmin
            z1=zmin
            z2=zmax
         else
            x1=xmax
            x2=xmax
            z1=zmin
            z2=zmax
         endif
      endif
      y1=ymin
      if (logy.eq.1) then
         v=fvy
         i=nint(v)
         k=0
         test=1
         test=test/10
         if (abs(v-i-subs(1)).lt.test) k=1
         if (abs(v-i-subs(4)).lt.test) k=4
         yt=y1
         if (k.gt.0) yt=y1-subs(k)/dvy
         k=k+1
      endif
      idone=0
      do while (idone.eq.0)
         if (logy.eq.0) then
            y1=y1+ystp/ny
         else
            y1=yt+subs(k)/dvy
            k=k+1
            if (k.gt.9) then
               yt=y1
               k=1
            endif
         endif
         test=1
         test=test/100
         if (abs(ymax-y1).lt.test) then
            idone=1
         else
            y2=y1
            call vect3(x1,y1,z1,x2,y2,z2,wg,0)
         endif
      enddo
   endif

   if (nz.gt.0) then
      if (nx.eq.0) then
         y1=ymin
         y2=ymax
         x1=xmin
         x2=xmin
      else if (nx.eq.-1) then
         y1=ymin
         y2=ymax
         x1=xmax
         x2=xmax
      else
         if (ny.eq.0) then
            y1=ymin
            y2=ymin
            x1=xmin
            x2=xmax
         else
            y1=ymax
            y2=ymax
            x1=xmin
            x2=xmax
         endif
      endif
      z1=zmin
      if (logz.eq.1) then
         v=fvz
         i=nint(v)
         k=0
         test=1
         test=test/10
         if (abs(v-i-subs(1)).lt.test) k=1
         if (abs(v-i-subs(4)).lt.test) k=4
         zt=z1
         if (k.gt.0) zt=z1-subs(k)/dvz
         k=k+1
      endif
      idone=0
      do while (idone.eq.0)
         if (logz.eq.0) then
            z1=z1+zstp/nz
         else
            z1=zt+subs(k)/dvz
            k=k+1
            if (k.gt.9) then
               zt=z1
               k=1
            endif
         endif
         test=1
         test=test/100
         if (abs(zmax-z1).lt.test) then
            idone=1
         else
            z2=z1
            call vect3(x1,y1,z1,x2,y2,z2,wg,0)
         endif
      enddo
   endif
   return
   end subroutine grid3

   subroutine curv2(x,y,n,idash,width,icon,isym,ssym,ishade)
   !--------------------------------------------------------------------
   ! Draw a 2-D curve.
   !--------------------------------------------------------------------
   ! externals
   integer::n,idash,icon,isym,ishade
   real(kr)::x(n),y(n),width,ssym
   ! internals
   integer::i
   real(kr)::w,wu1,wu2,wv1,wv2,xn,yn,zn,wu,wv,u,v,foregr
   real(kr),parameter::zero=0
   w=width

   !--draw curve, if requested.
   if (icon.ge.0) then
      wu1=100
      wu2=-100
      wv1=100
      wv2=-100
      xn=xscale(x(1))
      yn=yscale(y(1))
      zn=0
      call trans3(xn,yn,zn,wu,wv)
      if (wu.lt.wu1) wu1=wu
      if (wu.gt.wu2) wu2=wu
      if (wv.lt.wv1) wv1=wv
      if (wv.gt.wv2) wv2=wv
      call transw(wu,wv,u,v)
      call moveh(u,v)
      do i=2,n
         xn=xscale(x(i))
         yn=yscale(y(i))
         zn=0
         call trans3(xn,yn,zn,wu,wv)
         if (wu.lt.wu1) wu1=wu
         if (wu.gt.wu2) wu2=wu
         if (wv.lt.wv1) wv1=wv
         if (wv.gt.wv2) wv2=wv
         call transw(wu,wv,u,v)
         call drawh(u,v,w,idash)
      enddo

      !--fill the curve with a shading pattern
      foregr=0
      if (ishade.gt.0.and.ishade.le.10) then
         ifg=ishade+10
         call fillh(foregr)
         ifg=0
      else if (ishade.gt.40) then
         ifg=ishade
         call fillh(foregr)
         ifg=0
      endif

      !--fill the curve with a hatching pattern
      if (ishade.gt.10) call hatch(ishade,wu1,wv1,wu2,wv2)

      !--make curve invisible (used with shaded areas)
      if (width.eq.zero) call ncurve
   endif

   !--draw symbols, if requested.
   if (icon.ne.0) then
      do i=1,n
         if (mod(i,iabs(icon)).eq.0) then
            xn=xscale(x(i))
            yn=yscale(y(i))
            zn=0
            call trans3(xn,yn,zn,wu,wv)
            call transw(wu,wv,u,v)
            call dsym(u,v,isym,ssym)
         endif
      enddo
   endif
   return
   end subroutine curv2

   subroutine curv3(x,y,z,n,w)
   !--------------------------------------------------------------------
   ! Draw a 3D curve as a polygon with white fill.
   ! Later curves will hide earlier ones.
   !--------------------------------------------------------------------
   ! externals
   integer::n
   real(kr)::x(n+1),y(n+1),z(n+1),w
   ! internals
   integer::max,i
   real(kr)::xn,yn,zn,wu,wv,u,v,backgr

   !--add a segment on the baseline
   x(n+1)=x(1)
   y(n+1)=y(1)
   z(n+1)=z(1)
   max=n+1

   !--draw the polygon with white fill
   xn=xscale(x(1))
   yn=yscale(y(1))
   zn=zscale(z(1))
   call trans3(xn,yn,zn,wu,wv)
   call transw(wu,wv,u,v)
   call moveh(u,v)
   do i=2,max
      xn=xscale(x(i))
      yn=yscale(y(i))
      zn=zscale(z(i))
      call trans3(xn,yn,zn,wu,wv)
      call transw(wu,wv,u,v)
      call drawh(u,v,w,0)
   enddo
   backgr=1
   call fillh(backgr)
   return
   end subroutine curv3

   subroutine poly2(x,y,n,w,fill)
   !--------------------------------------------------------------------
   ! Draw a filled 2D polygon.
   ! If w<0, don't draw the outline.
   !--------------------------------------------------------------------
   ! externals
   integer::n
   real(kr)::x(n),y(n),w,fill
   ! internals
   integer::i

   call moveh(x(1),y(1))
   do i=2,n
      call drawh(x(i),y(i),w,0)
   enddo
   call fillh(fill)
   return
   end subroutine poly2

   subroutine poly3(x,y,z,n,w)
   !--------------------------------------------------------------------
   ! Draw a 3D polygon with hiding.
   !--------------------------------------------------------------------
   ! externals
   integer::n
   real(kr)::x(n),y(n),z(n),w
   ! internals
   integer::i
   real(kr)::u,v,wu,wv

   call trans3(x(1),y(1),z(1),wu,wv)
   call transw(wu,wv,u,v)
   call moveh(u,v)
   do i=2,n
      call trans3(x(i),y(i),z(i),wu,wv)
      call transw(wu,wv,u,v)
      call drawh(u,v,w,0)
   enddo
   return
   end subroutine poly3

   subroutine hatch(ihatch,wu1,wv1,wu2,wv2)
   !--------------------------------------------------------------------
   ! Fill the current path with hatching or cross hatching.
   !--------------------------------------------------------------------
   ! externals
   integer::ihatch
   real(kr)::wu1,wv1,wu2,wv2
   ! internals
   integer::iop
   real(kr)::dx,ww,x1,x2,y1,y2,u,v
   real(kr),parameter::w=.003e0_kr
   real(kr),parameter::oh2=.02e0_kr

   iop=(ihatch-1)/10
   dx=oh2*(10-ihatch+10*iop+1)
   ww=wv2-wv1
   call cclip
   if (iop.eq.1.or.iop.eq.3) then
      x1=wu1-dx*int(ww/dx)-dx
      y1=wv1
      y2=wv2
      do while (x1.lt.wu2+w)
         x1=x1+dx
         call transw(x1,y1,u,v)
         call moveh(u,v)
         x2=x1+ww
         call transw(x2,y2,u,v)
         call drawh(u,v,w,0)
      enddo
   endif
   if (iop.eq.2.or.iop.eq.3) then
      x1=wu2-dx*int((wu2-wu1+ww)/dx)-dx
      x2=wu2+w-1
      y1=wv1
      y2=wv2
      do while (x2.lt.wu2+w)
         x1=x1+dx
         call transw(x1,y1,u,v)
         call moveh(u,v)
         x2=x1-ww
         call transw(x2,y2,u,v)
         call drawh(u,v,w,0)
      enddo
   endif
   call nclip
   return
   end subroutine hatch

   subroutine text2(txt,ntxt,ifont,ht,xo,yo,xx,yx,xy,yy)
   !--------------------------------------------------------------------
   ! Draw a 2D text string using the 3d text call.
   !--------------------------------------------------------------------
   ! externals
   integer::ntxt,ifont
   character::txt*(*)
   real(kr)::ht,xo,yo,xx,yx,xy,yy
   ! internals
   real(kr),parameter::zero=0

   call text3(txt,ntxt,ifont,ht,xo,yo,zero,xx,yx,zero,xy,yy,zero)
   return
   end subroutine text2

   subroutine text3(txt,ntxt,ifont,ht,xo,yo,zo,xx,yx,zx,xy,yy,zy)
   !--------------------------------------------------------------------
   ! Draw a text string (including instructions) in 3D perspective.
   ! ifont and ht are the font and height for the characters.
   ! xo,yo,zo is the origen of the lower left corner of the text.
   ! xx,yx,zx is a unit vector in the character "x" direction.
   ! xy,yy,zy is a unit vector in the character "y" direction.
   !--------------------------------------------------------------------
   ! externals
   integer::ntxt,ifont
   character::txt*(*)
   real(kr)::ht,xo,yo,zo,xx,yx,zx,xy,yy,zy
   ! internals
   integer::icase,ialf,imode,jfont,i,j,ins,iskip,ii
   real(kr)::scale,x,y,z,delta,cbase,celev,dy,siz,wi
   real(kr)::xn,yn,zn,wu,wv,u,v,ux,vx,uy,vy
   character(1)::c,comm
   character(80)::temp

   scale=1
   x=xo
   y=yo
   z=zo
   delta=0

   cbase=0
   celev=0
   icase=0
   ialf=0
   imode=0
   jfont=ifont
   i=1
   j=0
   do while (i.le.ntxt)
      c=txt(i:i)
      if (imode.eq.1) then
         if ((lge(c,'0').and.lle(c,'9')).or.&
           c.eq.'.'.or.c.eq.'-') then
            j=j+1
            temp(j:j)=c
         else if (c.eq.'X') then
            j=j+1
            temp(j:j)=c
         else
            if (j.eq.0) then
               temp(1:1)='D'
               j=1
            endif
            if (comm.eq.'E') then
               if (temp(1:1).eq.'D') then
                  temp='.5'
                  j=2
               endif
               if (temp(1:1).eq.'X') then
                  celev=0
               else
                  dy=rget(temp,j)
                  if (dy.gt.0) then
                     celev=ht*dy
                  else
                     celev=-ht*dy
                  endif
               endif
               delta=celev
            else if (comm.eq.'L') then
               if (temp(1:1).eq.'D') then
                  temp='.5'
                  j=2
               endif
               if (temp(1:1).eq.'X') then
                  celev=0
               else
                  dy=rget(temp,j)
                  if (dy.gt.0) then
                     celev=-ht*dy
                  else
                     celev=ht*dy
                  endif
               endif
               delta=celev
            else if (comm.eq.'H') then
               if (temp(1:1).eq.'D') then
               temp='.7'
                  j=2
               endif
               if (temp(1:1).eq.'X') then
                  scale=1
               else
                  scale=rget(temp,j)
               endif
            else if (comm.eq.'O') then
               dy=rget(temp,j)
               cbase=cbase+dy
            else if (comm.eq.'F') then
               jfont=iget(temp,j)
            else if (comm.eq.'M') then
               ialf=iget(temp,j)
            endif
            comm=c
            j=0
            if (c.eq.'<'.or.c.eq.'>'.or.c.eq.'['.or.c.eq.']') then
               imode=0
               i=i-1
            endif
         endif
      else
         ins=0
         iskip=0
         if (c.eq.'#') then
            if (txt(i+1:i+1).eq.'#') then
               i=i+1
            else
               ins=1
            endif
         endif
         if (c.eq.'<') then
            if (txt(i+1:i+1).eq.'<') then
               i=i+1
            else
               jfont=ifont
               icase=1
               iskip=1
            endif
         endif
         if (c.eq.'>') then
            if (txt(i+1:i+1).eq.'>') then
               i=i+1
            else
               jfont=ifont
               icase=0
               iskip=1
            endif
         endif
         if (c.eq.'[') then
            if (txt(i+1:i+1).eq.'[') then
               i=i+1
            else
               jfont=3
               icase=1
               iskip=1
            endif
         endif
         if (c.eq.']') then
            if (txt(i+1:i+1).eq.']') then
               i=i+1
            else
               jfont=3
               icase=0
               iskip=1
            endif
         endif
         if (ins.eq.1) then
            imode=1
            comm=c
            j=0
         else if (iskip.eq.1) then
            iskip=0
         else
            if (icase.eq.1) then
               ii=ichar(c)
               if (ii.ge.96) ii=ii-32
               c=char(ii)
            endif
            siz=scale*ht
            wi=ssum(c,1,jfont,ialf,siz)
            xn=x+delta*xy
            yn=y+delta*yy
            zn=z+delta*zy
            call trans3(xn,yn,zn,wu,wv)
            call transw(wu,wv,u,v)
            call trans3(xn+xx,yn+yx,zn+zx,wu,wv)
            call transw(wu,wv,ux,vx)
            call trans3(xn+xy,yn+yy,zn+zy,wu,wv)
            call transw(wu,wv,uy,vy)
            ux=ux-u
            vx=vx-v
            uy=uy-u
            vy=vy-v
            call dchr(u,v,c,jfont,ialf,siz,ux,vx,uy,vy)
            x=x+wi*xx
            y=y+wi*yx
            z=z+wi*zx
         endif
      endif
      i=i+1
   enddo
   return
   end subroutine text3

   real(kr) function txtlen(txt,ntxt,ifont,ht)
   !--------------------------------------------------------------------
   ! Measure the length of a text string (with instructions), where
   ! ifont and ht give the font and size of the normal characters.
   !--------------------------------------------------------------------
   ! externals
   integer::ntxt,ifont
   character txt*(ntxt)
   real(kr)::ht
   ! internals
   integer::icase,ialf,imode,jfont,i,j,ins,iskip,ii,l
   real(kr)::siz,cbase,celev,scale,y
   character c*1,comm*1,num*4,temp*80

   siz=0
   cbase=0
   celev=0
   icase=0
   ialf=0
   scale=1
   imode=0
   iskip=0
   jfont=ifont
   i=1
   j=0
   do while (i.le.ntxt)
      c=txt(i:i)
      if (imode.eq.1) then
         if ((lge(c,'0').and.lle(c,'9')).or.c.eq.'.'.or.&
           c.eq.'-') then
            j=j+1
            temp(j:j)=c
         else if (c.eq.'X') then
            j=j+1
            temp(j:j)=c
         else
            if (j.eq.0) then
               temp(1:1)='D'
               j=1
            endif
            if (comm.eq.'E') then
               if (temp(1:1).eq.'D') then
                  temp='.5'
                  j=2
               endif
               if (temp(1:1).eq.'X') then
                  celev=0
               else
                  y=rget(temp,j)
                  if (y.gt.0) then
                     celev=ht*y
                  else
                     celev=-ht*y
                  endif
               endif
            else if (comm.eq.'L') then
               if (temp(1:1).eq.'D') then
                  temp='.5'
                  j=2
               endif
               if (temp(1:1).eq.'X') then
                  celev=0
               else
                  y=rget(temp,j)
                  if (y.gt.0) then
                     celev=-ht*y
                  else
                     celev=ht*y
                  endif
               endif
            else if (comm.eq.'H') then
               if (temp(1:1).eq.'D') then
               temp='.7'
                  j=2
               endif
               if (temp(1:1).eq.'X') then
                  scale=1
               else
                  scale=rget(temp,j)
               endif
            else if (comm.eq.'O') then
               y=rget(temp,j)
               cbase=cbase+y
            else if (comm.eq.'F') then
               jfont=iget(temp,j)
            else if (comm.eq.'M') then
               ialf=iget(temp,j)
            endif
            comm=c
            j=0
            if (c.eq.'<'.or.c.eq.'>'.or.c.eq.'['.or.c.eq.']') then
               imode=0
            endif
         endif
      else
         ins=0
         iskip=0
         if (i.lt.ntxt) then
            if (c.eq.'#') then
               if (txt(i+1:i+1).eq.'#') then
                  i=i+1
               else
                  ins=1
               endif
            endif
            if (c.eq.'<') then
               if (txt(i+1:i+1).eq.'<') then
                  i=i+1
               else
                  icase=1
                  jfont=ifont
                  iskip=1
               endif
            endif
            if (c.eq.'>') then
               if (txt(i+1:i+1).eq.'>') then
                  i=i+1
               else
                  icase=0
                  jfont=ifont
                  iskip=1
               endif
            endif
            if (c.eq.'[') then
               if (txt(i+1:i+1).eq.'[') then
                  i=i+1
               else
                  icase=1
                  jfont=3
                  iskip=1
               endif
            endif
            if (c.eq.']') then
               if (txt(i+1:i+1).eq.']') then
                  i=i+1
               else
                  icase=0
                  jfont=3
                  iskip=1
               endif
            endif
         endif
         if (ins.eq.1) then
            imode=1
            siz=siz+ssum(temp,l,jfont,ialf,scale)
            comm=c
            j=0
         else if (iskip.eq.1) then
            iskip=0
         else
            if (ialf.eq.0) then
               j=j+1
               if (icase.eq.1) then
                  ii=ichar(c)
                  if (ii.ge.96) ii=ii-32
                  c=char(ii)
               endif
               temp(j:j)=c
               if (c.ne.' ') l=j
            else
               j=j+1
               write(num,'(a1,o3.3)') char(92),ichar(c)+128
               temp(j:j+3)=num(1:4)
               j=j+3
               l=j
            endif
         endif
      endif
      i=i+1
   enddo
   if (j.gt.0) siz=siz+ssum(temp,l,jfont,ialf,scale)
   txtlen=siz*ht
   return
   end function txtlen

   function ssum(temp,l,ifont,ialf,scale)
   !--------------------------------------------------------------------
   ! Private routine for txtlen.  Sums text width.
   !--------------------------------------------------------------------
   ! externals
   integer::l,ifont,ialf
   real(kr)::scale
   character::temp*(l)
   ! internals
   integer::itabl,i
   real(kr)::ssum

   itabl=ifont
   ssum=0
   do i=1,l
      ssum=ssum+scale*charw(temp(i:i),itabl)
   enddo
   return
   end function ssum

   integer function iget(strng,n)
   !--------------------------------------------------------------------
   ! Decode an integer out of a string of length n.
   !--------------------------------------------------------------------
   ! externals
   integer::n
   character::strng*(*)
   ! internals
   integer::neg,ival,i,j
   character c*1
   character(10),parameter::digits='0123456789'

   neg=0
   ival=0
   do i=1,n
      c=strng(i:i)
      j=index(digits,c)
      if (j.gt.0) then
         ival=10*ival+j-1
      else if (c.eq.'-') then
         neg=1
      endif
   enddo
   if (neg.eq.1) ival=-ival
   iget=ival
   return
   end function iget

   real(kr) function rget(strng,n)
   !--------------------------------------------------------------------
   ! Decode a real out of a string of length n.
   !--------------------------------------------------------------------
   ! externals
   integer::n
   character::strng*(*)
   ! internals
   integer::neg,ival,ipwr,i,j
   real(kr)::val
   character(1)::c
   character(10),parameter::digits='0123456789'
   real(kr),parameter::ten=10

   neg=0
   ival=0
   ipwr=0
   do i=1,n
      c=strng(i:i)
      j=index(digits,c)
      if (j.gt.0) then
         ival=10*ival+j-1
         if (ipwr.ne.0) ipwr=ipwr+1
      else if (c.eq.'.') then
         ipwr=1
      else if (c.eq.'-') then
         neg=1
      endif
   enddo
   val=ival
   if (ipwr.ne.0) val=val*ten**(1-ipwr)
   if (neg.eq.1) val=-val
   rget=val
   return
   end function rget

!----------------------------------------------------------------------
!
! Postscript dependent routines for the plotting package
!
!----------------------------------------------------------------------

   subroutine gplot(lori,xpage,ypage,nnps)
   !--------------------------------------------------------------------
   ! Initialize plotting.
   ! Set up for US letter size paper (xpaper=8.5in, ypaper=11.0in).
   ! See separate setting for xpaper and ypaper at the start of viewr.
   ! There are also default settings for the page size in plotr.
   !--------------------------------------------------------------------
   use util ! provides openz
   ! externals
   integer::lori,nnps
   real(kr)::xpage,ypage
   ! internals
   integer::i1,i2,i3,i4

   xpaper=8.5e0_kr
   ypaper=11.0e0_kr

   ipage=0
   nps=nnps
   call openz(nps,1)
   write(nps,'(''%!PS-Adobe-'')')
   if (lori.eq.0) then
      i1=int((xpaper-xpage)*72/2)
      i2=int((ypaper-ypage)*72/2)
      i3=int(i1+xpage*72)
      i4=int(i2+ypage*72)
      ushift=i1
      vshift=i2
      uwidth=xpage*72
   else
      i1=int((ypaper-xpage)*72/2)
      i2=int((xpaper-ypage)*72/2)
      i3=int(i1+ypage*72)
      i4=int(i2+xpage*72)
      ushift=i2
      vshift=i1
      uwidth=ypage*72
   endif
   write(nps,'(''%%BoundingBox: '',4i5)') i1,i2,i3,i4
   write(nps,'(''%%Pages: (atend)'')')
   if (lori.eq.1) write(nps,'(''%%Orientation: Landscape'')')
   if (lori.eq.0) write(nps,'(''%%Orientation: Portrait'')')
   write(nps,'(''2 setlinecap'')')
   return
   end subroutine gplot

   subroutine gdone
   !--------------------------------------------------------------------
   ! Finalize plotting.
   !--------------------------------------------------------------------
   use util ! provides closz

   write(nps,'(''%gdone'')')
   if (ipage.lt.10) then
      write(nps,'(''%%Trailer'',/,''%%Pages: '',i1)') ipage
   else if (ipage.lt.100)then
      write(nps,'(''%%Trailer'',/,''%%Pages: '',i2)') ipage
   else
      write(nps,'(''%%Trailer'',/,''%%Pages: '',i3)') ipage
   endif
   write(nps,'(''%%EOF'')')
   call closz(nps)
   return
   end subroutine gdone

   subroutine gset(ull,vll,uur,vur)
   !--------------------------------------------------------------------
   ! Set up a clipping path.
   !--------------------------------------------------------------------
   ! externals
   real(kr)::ull,vll,uur,vur
   ! internals
   real(kr)::w
   real(kr),parameter::width=.005e0_kr

   write(nps,'(''%gset'')')
   call moveh(ull,vll)
   w=width
   call drawh(uur,vll,w,0)
   call drawh(uur,vur,w,0)
   call drawh(ull,vur,w,0)
   call drawh(ull,vll,w,0)
   write(nps,'('' gsave clip newpath'')')
   return
   end subroutine gset

   subroutine gend
   !--------------------------------------------------------------------
   ! End a clipping path.
   !--------------------------------------------------------------------
   write(nps,'(''stroke grestore newpath''/''%gend'')')
   return
   end subroutine gend

   subroutine newp
   !--------------------------------------------------------------------
   ! Start a plotting page.
   !--------------------------------------------------------------------
   ipage=ipage+1
   write(nps,'(''%%Page:'',i4,i4)') ipage,ipage
   return
   end subroutine newp

   subroutine endp
   !--------------------------------------------------------------------
   ! End a plotting page.
   !--------------------------------------------------------------------
   write(nps,'(''stroke''/''showpage''/''%endp'')')
   return
   end subroutine endp

   subroutine moveh(x,y)
   !--------------------------------------------------------------------
   ! Low-level move routine for Postscript
   ! with portrait-landscape rotation.
   !--------------------------------------------------------------------
   ! externals
   real(kr)::x,y
   ! internals
   real(kr)::u1,v1,rgb,r,g,b,test

   if (land.eq.1) then
      u1=uwidth-72*y+ushift
      v1=72*x+vshift
   else
      u1=72*x+ushift
      v1=72*y+vshift
   endif
   write(nps,'(''stroke''/''newpath'')')
   rgb=255
   r=ifrgb(1,ifg)/rgb
   g=ifrgb(2,ifg)/rgb
   b=ifrgb(3,ifg)/rgb
   write(nps,'(3f6.3,'' setrgbcolor'')') r,g,b
   test=2000
   if (u1.gt.test) u1=2000
   test=-1000
   if (u1.lt.test) u1=-1000
   test=2000
   if (v1.gt.test) v1=2000
   test=-1000
   if (v1.lt.test) v1=-1000
   write(nps,'(2f9.2,'' moveto'')') u1,v1
   return
   end subroutine moveh

   subroutine drawh(x,y,w,idash)
   !--------------------------------------------------------------------
   ! Low-level draw routine for Postscript
   ! with portrait-landscape rotation
   ! and line-type control.
   ! If w is negative, don't draw the line.
   !--------------------------------------------------------------------
   ! externals
   integer::idash
   real(kr)::x,y,w
   ! internals
   real(kr)::u,v,test
   real(kr)::wlast=0.e0_kr
   integer::ldash=-1
   real(kr),parameter::zero=0

   if (w.lt.zero) then
      wlast=w
      ldash=idash
      return
   endif
   if (w.ne.wlast) then
      write(nps,'(f8.3,'' setlinewidth'')') 72*w
      wlast=w
   endif
   if (idash.ne.ldash) then
      if (idash.eq.0) then
         write(nps,'(''[] 0 setdash'')')
      else if (idash.eq.1) then
         write(nps,'(''['',2f9.2,''] 0 setdash'')')&
           72*w*5,72*w*3
      else if (idash.eq.2) then
         write(nps,'(''['',4f9.2,''] 0 setdash'')')&
           72*w*6,72*w*3,72*w*3,72*w*3
      else if (idash.eq.3) then
         write(nps,'(''['',4f9.2,''] 0 setdash'')')&
           72*w*6,72*w*3,72*w*1,72*w*3
      else if (idash.eq.4) then
         write(nps,'(''['',2f9.2,''] 0 setdash'')')&
           72*w*1,72*w*3
      endif
      ldash=idash
   endif
   if (land.eq.1) then
      u=uwidth-72*y+ushift
      v=72*x+vshift
   else
      u=72*x+ushift
      v=72*y+vshift
   endif
   test=2000
   if (u.gt.test) u=2000
   test=-1000
   if (u.lt.test) u=-1000
   test=2000
   if (v.gt.test) v=2000
   test=-1000
   if (v.lt.test) v=-1000
   write(nps,'(2f9.2,'' lineto'')') u,v
   return
   end subroutine drawh

   subroutine fillh(color)
   !--------------------------------------------------------------------
   ! Fill current path with background (1.) or foreground (0.) color.
   ! This may be a discrete color, or one of a progression of
   ! grays or shading colors used to show values.
   !--------------------------------------------------------------------
   ! externals
   real(kr)::color
   ! internals
   real(kr)::rgb,cmin,cmax,r,g,b,ten

   rgb=255
   cmin=1
   cmin=cmin/100
   cmax=1
   cmax=1-cmax/100
   if (color.gt.cmax) then
      r=ibrgb(1,ibg)/rgb
      g=ibrgb(2,ibg)/rgb
      b=ibrgb(3,ibg)/rgb
   else if (color.lt.cmin.and.ifg.le.10) then
      ! curve colors
      r=ifrgb(1,ifg)/rgb
      g=ifrgb(2,ifg)/rgb
      b=ifrgb(3,ifg)/rgb
   else if (color.lt.cmin.and.ifg.le.20) then
      ! shades of gray
      ten=10
      r=(20-ifg)/ten
      g=(20-ifg)/ten
      b=(20-ifg)/ten
   else if (color.lt.cmin.and.ifg.gt.20) then
      ! filling in shades of curve colors
      r=isrgb(1,ifg-40)/rgb
      g=isrgb(2,ifg-40)/rgb
      b=isrgb(3,ifg-40)/rgb
   endif
   write(nps,'(''gsave'',3f6.3,'' setrgbcolor fill grestore'')') r,g,b
   return
   end subroutine fillh

   subroutine cclip
   !--------------------------------------------------------------------
   ! Set up a clipping path on the current path.
   !--------------------------------------------------------------------
   write(nps,'(''%cclip''/''gsave clip newpath'')')
   return
   end subroutine cclip

   subroutine nclip
   !--------------------------------------------------------------------
   ! Terminate clipping path on the current path.
   !--------------------------------------------------------------------
   write(nps,'(''stroke grestore''/''%nclip'')')
   return
   end subroutine nclip

   subroutine ncurve
   !--------------------------------------------------------------------
   ! Cancel current curve path for shaded regions without borders.
   !--------------------------------------------------------------------
   write(nps,'(''newpath'')')
   return
   end subroutine ncurve

   subroutine dsym(x,y,isym,size)
   !--------------------------------------------------------------------
   ! Low-level symbol routine using Postscript.
   !--------------------------------------------------------------------
   ! externals
   integer::isym
   real(kr)::x,y,size
   ! internals
   real(kr)::foregr,backgr,sz
   real(kr)::xx(10),yy(10)
   ! character stroke width
   real(kr),parameter::w=0.01e0_kr
   real(kr),parameter::rr=1.128e0_kr
   real(kr),parameter::c8=.414e0_kr
   real(kr),parameter::cd=1.4e0_kr
   real(kr),parameter::cx=.82e0_kr
   real(kr),parameter::t1=1.75e0_kr
   real(kr),parameter::t2=1.52e0_kr
   real(kr),parameter::t3=.88e0_kr
   sz=6*size/10
   foregr=0
   backgr=1

   !--square
   if (isym.eq.0) then
      xx(1)=x-sz
      yy(1)=y+sz
      xx(2)=x+sz
      yy(2)=y+sz
      xx(3)=x+sz
      yy(3)=y-sz
      xx(4)=x-sz
      yy(4)=y-sz
      xx(5)=x-sz
      yy(5)=y+sz
      call poly2(xx,yy,5,w,backgr)

   !--octagon
   else if (isym.eq.1) then
      xx(1)=x-c8*sz
      yy(1)=y+sz
      xx(2)=x+c8*sz
      yy(2)=y+sz
      xx(3)=x+sz
      yy(3)=y+c8*sz
      xx(4)=x+sz
      yy(4)=y-c8*sz
      xx(5)=x+c8*sz
      yy(5)=y-sz
      xx(6)=x-c8*sz
      yy(6)=y-sz
      xx(7)=x-sz
      yy(7)=y-c8*sz
      xx(8)=x-sz
      yy(8)=y+c8*sz
      xx(9)=x-c8*sz
      yy(9)=y+sz
      call poly2(xx,yy,9,w,backgr)

   !--triangle
   else if (isym.eq.2) then
      xx(1)=x-t2*sz
      yy(1)=y-t3*sz
      xx(2)=x
      yy(2)=y+t1*sz
      xx(3)=x+t2*sz
      yy(3)=yy(1)
      xx(4)=xx(1)
      yy(4)=yy(1)
      call poly2(xx,yy,4,w,backgr)

   !--cross
   else if (isym.eq.3) then
      call moveh(x,y+sz)
      call drawh(x,y-sz,w,0)
      call moveh(x-sz,y)
      call drawh(x+sz,y,w,0)

   !--ex
   else if (isym.eq.4) then
      call moveh(x-sz,y+sz)
      call drawh(x+sz,y-sz,w,0)
      call moveh(x-sz,y-sz)
      call drawh(x+sz,y+sz,w,0)

   !--diamond
   else if (isym.eq.5) then
      xx(1)=x-cd*sz
      yy(1)=y
      xx(2)=x
      yy(2)=y+cd*sz
      xx(3)=x+cd*sz
      yy(3)=y
      xx(4)=x
      yy(4)=y-cd*sz
      xx(5)=xx(1)
      yy(5)=yy(1)
      call poly2(xx,yy,5,w,backgr)

   !--inverted triangle
   else if (isym.eq.6) then
      xx(1)=x-t2*sz
      yy(1)=y+t3*sz
      xx(2)=x
      yy(2)=y-t1*sz
      xx(3)=x+t2*sz
      yy(3)=yy(1)
      xx(4)=xx(1)
      yy(4)=yy(1)
      call poly2(xx,yy,4,w,backgr)

   !--exed square
   else if (isym.eq.7) then
      xx(1)=x-sz
      yy(1)=y+sz
      xx(2)=x+sz
      yy(2)=y+sz
      xx(3)=x+sz
      yy(3)=y-sz
      xx(4)=x-sz
      yy(4)=y-sz
      xx(5)=x-sz
      yy(5)=y+sz
      call poly2(xx,yy,5,w,backgr)
      call moveh(x-sz,y+sz)
      call drawh(x+sz,y-sz,w,0)
      call moveh(x-sz,y-sz)
      call drawh(x+sz,y+sz,w,0)

   !--crossed ex
   else if (isym.eq.8) then
      call moveh(x-sz,y+sz)
      call drawh(x+sz,y-sz,w,0)
      call moveh(x-sz,y-sz)
      call drawh(x+sz,y+sz,w,0)
      call moveh(x-sz,y)
      call drawh(x+sz,y,w,0)
      call moveh(x,y+sz)
      call drawh(x,y-sz,w,0)

   !--crossed diamond
   else if (isym.eq.9) then
      xx(1)=x-cd*sz
      yy(1)=y
      xx(2)=x
      yy(2)=y+cd*sz
      xx(3)=x+cd*sz
      yy(3)=y
      xx(4)=x
      yy(4)=y-cd*sz
      xx(5)=xx(1)
      yy(5)=yy(1)
      call poly2(xx,yy,5,w,backgr)
      call moveh(x,y+cd*sz)
      call drawh(x,y-cd*sz,w,0)
      call moveh(x-cd*sz,y)
      call drawh(x+cd*sz,y,w,0)

   !--crossed octagon
   else if (isym.eq.10) then
      xx(1)=x-c8*sz
      yy(1)=y+sz
      xx(2)=x+c8*sz
      yy(2)=y+sz
      xx(3)=x+sz
      yy(3)=y+c8*sz
      xx(4)=x+sz
      yy(4)=y-c8*sz
      xx(5)=x+c8*sz
      yy(5)=y-sz
      xx(6)=x-c8*sz
      yy(6)=y-sz
      xx(7)=x-sz
      yy(7)=y-c8*sz
      xx(8)=x-sz
      yy(8)=y+c8*sz
      xx(9)=x-c8*sz
      yy(9)=y+sz
      call poly2(xx,yy,9,w,backgr)
      call moveh(x,y+sz)
      call drawh(x,y-sz,w,0)
      call moveh(x-sz,y)
      call drawh(x+sz,y,w,0)

   !--double triangle
   else if (isym.eq.11) then
      xx(1)=x-t2*sz
      yy(1)=y-t3*sz
      xx(2)=x
      yy(2)=y+t1*sz
      xx(3)=x+t2*sz
      yy(3)=yy(1)
      xx(4)=xx(1)
      yy(4)=yy(1)
      call poly2(xx,yy,4,w,backgr)
      xx(1)=x-t2*sz
      yy(1)=y+t3*sz
      xx(2)=x
      yy(2)=y-t1*sz
      xx(3)=x+t2*sz
      yy(3)=yy(1)
      xx(4)=xx(1)
      yy(4)=yy(1)
      call poly2(xx,yy,4,w,backgr)

   !--crossed square
   else if (isym.eq.12) then
      xx(1)=x-sz
      yy(1)=y+sz
      xx(2)=x+sz
      yy(2)=y+sz
      xx(3)=x+sz
      yy(3)=y-sz
      xx(4)=x-sz
      yy(4)=y-sz
      xx(5)=x-sz
      yy(5)=y+sz
      call poly2(xx,yy,5,w,backgr)
      call moveh(x,y+sz)
      call drawh(x,y-sz,w,0)
      call moveh(x-sz,y)
      call drawh(x+sz,y,w,0)

   !--exed octagon
   else if (isym.eq.13) then
      xx(1)=x-c8*sz
      yy(1)=y+sz
      xx(2)=x+c8*sz
      yy(2)=y+sz
      xx(3)=x+sz
      yy(3)=y+c8*sz
      xx(4)=x+sz
      yy(4)=y-c8*sz
      xx(5)=x+c8*sz
      yy(5)=y-sz
      xx(6)=x-c8*sz
      yy(6)=y-sz
      xx(7)=x-sz
      yy(7)=y-c8*sz
      xx(8)=x-sz
      yy(8)=y+c8*sz
      xx(9)=x-c8*sz
      yy(9)=y+sz
      call poly2(xx,yy,9,w,backgr)
      call moveh(x-sz,y+sz)
      call drawh(x+sz,y-sz,w,0)
      call moveh(x-sz,y-sz)
      call drawh(x+sz,y+sz,w,0)

   !--triangle and square
   else if (isym.eq.14) then
      xx(1)=x-sz
      yy(1)=y+sz
      xx(2)=x+sz
      yy(2)=y+sz
      xx(3)=x+sz
      yy(3)=y-sz
      xx(4)=x-sz
      yy(4)=y-sz
      xx(5)=x-sz
      yy(5)=y+sz
      call poly2(xx,yy,5,w,backgr)
      call moveh(x-sz,y-sz)
      call drawh(x,y+sz,w,0)
      call drawh(x+sz,y-sz,w,0)

   !--filled circle
   else if (isym.eq.15) then
      call circle(x,y,rr*sz,w)
      call fillh(foregr)

   !--circle
   else if (isym.eq.16) then
      call circle(x,y,rr*sz,w)
      call fillh(backgr)

   !--open square
   else if (isym.eq.17) then
      xx(1)=x-sz
      yy(1)=y+sz
      xx(2)=x+sz
      yy(2)=y+sz
      xx(3)=x+sz
      yy(3)=y-sz
      xx(4)=x-sz
      yy(4)=y-sz
      xx(5)=x-sz
      yy(5)=y+sz
      call poly2(xx,yy,5,w,backgr)

   !--filled square
   else if (isym.eq.18) then
      xx(1)=x-sz
      yy(1)=y+sz
      xx(2)=x+sz
      yy(2)=y+sz
      xx(3)=x+sz
      yy(3)=y-sz
      xx(4)=x-sz
      yy(4)=y-sz
      xx(5)=x-sz
      yy(5)=y+sz
      call poly2(xx,yy,5,w,foregr)

   !--filled diamond
   else if (isym.eq.19) then
      xx(1)=x-cd*sz
      yy(1)=y
      xx(2)=x
      yy(2)=y+cd*sz
      xx(3)=x+cd*sz
      yy(3)=y
      xx(4)=x
      yy(4)=y-cd*sz
      xx(5)=xx(1)
      yy(5)=yy(1)
      call poly2(xx,yy,5,w,foregr)

   !--filled triangle
   else if (isym.eq.20) then
      xx(1)=x-t2*sz
      yy(1)=y-t3*sz
      xx(2)=x
      yy(2)=y+t1*sz
      xx(3)=x+t2*sz
      yy(3)=yy(1)
      xx(4)=xx(1)
      yy(4)=yy(1)
      call poly2(xx,yy,4,w,foregr)

   !--filled inverted triangle
   else if (isym.eq.21) then
      xx(1)=x-t2*sz
      yy(1)=y+t3*sz
      xx(2)=x
      yy(2)=y-t1*sz
      xx(3)=x+t2*sz
      yy(3)=yy(1)
      xx(4)=xx(1)
      yy(4)=yy(1)
      call poly2(xx,yy,4,w,foregr)

   !--crossed circle
   else if (isym.eq.22) then
      call circle(x,y,rr*sz,w)
      call fillh(backgr)
      call moveh(x,y+rr*sz)
      call drawh(x,y-rr*sz,w,0)
      call moveh(x-rr*sz,y)
      call drawh(x+rr*sz,y,w,0)

   !--exed circle
   else if (isym.eq.23) then
      call circle(x,y,rr*sz,w)
      call fillh(backgr)
      call moveh(x-cx*sz,y+cx*sz)
      call drawh(x+cx*sz,y-cx*sz,w,0)
      call moveh(x-cx*sz,y-cx*sz)
      call drawh(x+cx*sz,y+cx*sz,w,0)

   !--exed diamond
   else if (isym.eq.24) then
      xx(1)=x-cd*sz
      yy(1)=y
      xx(2)=x
      yy(2)=y+cd*sz
      xx(3)=x+cd*sz
      yy(3)=y
      xx(4)=x
      yy(4)=y-cd*sz
      xx(5)=xx(1)
      yy(5)=yy(1)
      call poly2(xx,yy,5,w,backgr)
      call moveh(x-cd*sz/2,y+cd*sz/2)
      call drawh(x+cd*sz/2,y-cd*sz/2,w,0)
      call moveh(x-cd*sz/2,y-cd*sz/2)
      call drawh(x+cd*sz/2,y+cd*sz/2,w,0)

   !--default (circle)
   else
      call circle(x,y,rr*sz,w)
      call fillh(backgr)
   endif

   return
   end subroutine dsym

   subroutine circle(x,y,r,w)
   !--------------------------------------------------------------------
   ! Draw a circle.
   !--------------------------------------------------------------------
   ! externals
   real(kr)::x,y,r,w
   ! internals
   real(kr)::u1,v1,rr,test

   if (land.eq.1) then
      u1=uwidth-72*y+ushift
      v1=72*x+vshift
      rr=72*r
   else
      u1=72*x+ushift
      rr=72*r
      v1=72*y+vshift
   endif
   write(nps,'(''stroke'')')
   write(nps,'(''newpath'')')
   test=2000
   if (u1.gt.test) u1=2000
   test=-1000
   if (u1.lt.test) u1=-1000
   test=2000
   if (v1.gt.test) v1=2000
   test=-1000
   if (v1.lt.test) v1=-1000
   write(nps,'(1x,3f9.2,'' 0 360 arc'')') u1,v1,rr
   return
   end subroutine circle

   subroutine dchr(x,y,c,ifont,ialf,ht,xx,yx,xy,yy)
   !--------------------------------------------------------------------
   ! Low-level routine to draw a character in perspective
   ! using Postscript.
   !--------------------------------------------------------------------
   ! externals
   integer::ifont,ialf
   real(kr)::x,y,ht,xx,yx,xy,yy
   character(1)::c
   ! internals
   integer::lname
   real(kr)::u,v,ux,vx,uy,vy,test
   character(25)::name

   if (land.eq.1) then
      u=uwidth-72*y+ushift
      v=72*x+vshift
      ux=-72*yx
      vx=72*xx
      uy=-72*yy
      vy=72*xy
   else
      u=72*x+ushift
      v=72*y+vshift
      ux=72*xx
      vx=72*yx
      uy=72*xy
      vy=72*yy
   endif
   test=2000
   if (u.gt.test) u=2000
   test=-1000
   if (u.lt.test) u=-1000
   test=2000
   if (v.gt.test) v=2000
   test=-1000
   if (v.lt.test) v=-1000
   write(nps,'(2f9.2,'' moveto'')') u,v
   name=font(ifont)
   call stripv(name,lname)
   write(nps,'(''/'',a,'' findfont ['',4f7.2,&
     &'' 0. 0.] makefont setfont'')') name(1:lname),&
     ht*ux,ht*vx,ht*uy,ht*vy
   if (ialf.eq.0) then
      if (c.eq.'('.or.c.eq.')') then
         write(nps,'(''('',a1,a1,'') show'')') '\\',c
      else
         write(nps,'(''('',a1,'') show'')') c
      endif
   else
      write(nps,'(a1,o3.3,'') show'')') char(92),ichar(c)+128
   endif
   return
   end subroutine dchr

   real(kr) function charw(char,itable)
   !--------------------------------------------------------------------
   ! Return character widths for built-in fonts.
   !--------------------------------------------------------------------
   ! externals
   integer::itable
   character(1)::char
   ! internals
   integer::i
   i=1+ichar(char)
   ! i=1+int(char)
   charw=cw(i,itable)
   return
   end function charw

   subroutine stripv(str,length)
   !--------------------------------------------------------------------
   ! Strip leading and trailing blanks.  Count the remainder.
   !--------------------------------------------------------------------
   ! externals
   integer::length
   character::str*(*)
   ! internals
   integer::n,j,k,i,ic

   n=len(str)
   j=0
   k=0
   do i=1,n
      ic=ichar(str(i:i))
      ! ic=int(str(i:i))
      if (ic.gt.32.and.ic.lt.127) then
         if (j.eq.0) j=i
         k=i
      endif
   enddo
   length=0
   if (j.eq.0) return
   length=k-j+1
   if (j.eq.1) return
   k=j
   do i=1,length
      str(i:i)=str(k:k)
      k=k+1
   enddo
   return
   end subroutine stripv

end module graph


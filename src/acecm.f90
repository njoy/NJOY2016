module acecm
   ! provides some routines that are common to more than one
   ! of the submodules of acer; mtname, ptleg2, pttab2, bachaa.
   use locale
   implicit none
   private

   !--Public routines
   public mtname,ptleg2,pttab2,bachaa,eavl

contains

   subroutine mtname(mt,name,izai)
   !-------------------------------------------------------------------
   ! Return the reaction name for an ENDF MT number.
   !-------------------------------------------------------------------
   use endf ! provides iverf
   ! externals
   integer::mt,izai
   character(10)::name
   ! internals
   integer::i
   character(10),dimension(500),parameter::hndf=(/&
     'total     ','elastic   ','nonelastic','inelastic ','(n,x)     ',&
     '(n,1/2*1) ','(n,1/2*2) ','(n,1/2*3) ','(n,1/2*4) ','(n,x)     ',&
     '(n,2nd)   ','(n,x)     ','(n,x)     ','(n,x)     ','(n,x)     ',&
     '(n,2n)    ','(n,3n)    ','fission   ','(n,f)     ','(n,n*f)   ',&
     '(n,2nf)   ','(n,n*)a   ','(n,n*)3a  ','(n,2n)a   ','(n,3n)a   ',&
     '(n,2n)iso ','(n,abs)   ','(n,n*)p   ','(n,n*)2a  ','(n,2n)2a  ',&
     '(n,x)     ','(n,n*)d   ','(n,n*)t   ','(n,n*)he3 ','(n,n*)d2a ',&
     '(n,n*)t2a ','(n,4n)    ','(n,3nf)   ','(n,x)     ','(n,x)     ',&
     '(n,2np)   ','(n,3np)   ','(n,x)     ','(n,2np)   ','(n,npa)   ',&
     '(n,2/2*1) ','(n,2/2*2) ','(n,2/2*3) ','(n,2/2*4) ','(n,n*0)   ',& !50
     '(n,n*1)   ','(n,n*2)   ','(n,n*3)   ','(n,n*4)   ','(n,n*5)   ',&
     '(n,n*6)   ','(n,n*7)   ','(n,n*8)   ','(n,n*9)   ','(n,n*10)  ',&
     '(n,n*11)  ','(n,n*12)  ','(n,n*13)  ','(n,n*14)  ','(n,n*15)  ',&
     '(n,n*16)  ','(n,n*17)  ','(n,n*18)  ','(n,n*19)  ','(n,n*20)  ',&
     '(n,n*21)  ','(n,n*22)  ','(n,n*23)  ','(n,n*24)  ','(n,n*25)  ',&
     '(n,n*26)  ','(n,n*27)  ','(n,n*28)  ','(n,n*29)  ','(n,n*30)  ',&
     '(n,n*31)  ','(n,n*32)  ','(n,n*33)  ','(n,n*34)  ','(n,n*35)  ',&
     '(n,n*36)  ','(n,n*37)  ','(n,n*38)  ','(n,n*39)  ','(n,n*40)  ',&
     '(n,n*c)   ','(n,x)     ','(n,x)     ','(n,x)     ','(n,x)     ',&
     '(n,x)     ','(n,x)     ','(n,x)     ','(n,n*)gma ','(n,x)     ',& !100
     '(n,parab) ','(n,gma)   ','(n,p)     ','(n,d)     ','(n,t)     ',&
     '(n,he3)   ','(n,a)     ','(n,2a)    ','(n,3a)    ','(n,x)     ',&
     '(n,2p)    ','(n,pa)    ','(n,t2a)   ','(n,d2a)   ','(n,pd)    ',&
     '(n,pt)    ','(n,da)    ','(n,x)     ','(n,x)     ','(n,dest)  ',&
     '(n,x)     ','(n,x)     ','(n,x)     ','(n,x)     ','(n,x)     ',&
     '(n,x)     ','(n,x)     ','(n,x)     ','(n,x)     ','(n,x)     ',&
     '(n,x)     ','(n,x)     ','(n,x)     ','(n,x)     ','(n,x)     ',&
     '(n,x)     ','(n,x)     ','(n,x)     ','(n,x)     ','(n,x)     ',&
     '(n,x)     ','(n,x)     ','(n,x)     ','(n,x)     ','(n,x)     ',&
     '(n,x)     ','(n,x)     ','(n,x)     ','(n,x)     ','(n,x)     ',& !150
     '(n,x)     ','(n,5n)    ','(n,6n)    ','(n,2nt)   ','(n,ta)    ',&
     '(n,4np)   ','(n,3nd)   ','(n,nda)   ','(n,2npa)  ','(n,7n)    ',&
     '(n,8n)    ','(n,5np)   ','(n,6np)   ','(n,7np)   ','(n,4na)   ',&
     '(n,5na)   ','(n,6na)   ','(n,7na)   ','(n,4nd)   ','(n,5nd)   ',&
     '(n,6nd)   ','(n,3nt)   ','(n,4nt)   ','(n,5nt)   ','(n,6nt)   ',&
     '(n,2nhe3) ','(n,3nhe3) ','(n,4nhe3) ','(n,3n2p)  ','(n,3n2a)  ',&
     '(n,3npa)  ','(n,dt)    ','(n,npd)   ','(n,npt)   ','(n,ndt)   ',&
     '(n,nphe3) ','(n,ndhe3) ','(n,nthe3) ','(n,nta)   ','(n,2n2p)  ',&
     '(n,phe3)  ','(n,dhe3)  ','(n,he3a)  ','(n,4n2p)  ','(n,2n2a)  ',&
     '(n,4npa)  ','(n,3p)    ','(n,n3p)   ','(n,3n2pa) ','(n,5n2p)  ',& !200
     '(n,p*0)   ',&
     '(n,p*1)   ','(n,p*2)   ','(n,p*3)   ','(n,p*4)   ','(n,p*5)   ',&
     '(n,p*6)   ','(n,p*7)   ','(n,p*8)   ','(n,p*9)   ','(n,p*10)  ',&
     '(n,p*11)  ','(n,p*12)  ','(n,p*13)  ','(n,p*14)  ','(n,p*15)  ',&
     '(n,p*16)  ','(n,p*17)  ','(n,p*18)  ','(n,p*19)  ','(n,p*20)  ',&
     '(n,p*21)  ','(n,p*22)  ','(n,p*23)  ','(n,p*24)  ','(n,p*25)  ',&
     '(n,p*26)  ','(n,p*27)  ','(n,p*28)  ','(n,p*29)  ','(n,p*30)  ',&
     '(n,p*31)  ','(n,p*32)  ','(n,p*33)  ','(n,p*34)  ','(n,p*35)  ',&
     '(n,p*36)  ','(n,p*37)  ','(n,p*38)  ','(n,p*39)  ','(n,p*40)  ',&
     '(n,p*41)  ','(n,p*42)  ','(n,p*43)  ','(n,p*44)  ','(n,p*45)  ',&
     '(n,p*46)  ','(n,p*47)  ','(n,p*48)  ','(n,p*c)   ',             & !250
     '(n,d*0)   ',&
     '(n,d*1)   ','(n,d*2)   ','(n,d*3)   ','(n,d*4)   ','(n,d*5)   ',&
     '(n,d*6)   ','(n,d*7)   ','(n,d*8)   ','(n,d*9)   ','(n,d*10)  ',&
     '(n,d*11)  ','(n,d*12)  ','(n,d*13)  ','(n,d*14)  ','(n,d*15)  ',&
     '(n,d*16)  ','(n,d*17)  ','(n,d*18)  ','(n,d*19)  ','(n,d*20)  ',&
     '(n,d*21)  ','(n,d*22)  ','(n,d*23)  ','(n,d*24)  ','(n,d*25)  ',&
     '(n,d*26)  ','(n,d*27)  ','(n,d*28)  ','(n,d*29)  ','(n,d*30)  ',&
     '(n,d*31)  ','(n,d*32)  ','(n,d*33)  ','(n,d*34)  ','(n,d*35)  ',&
     '(n,d*36)  ','(n,d*37)  ','(n,d*38)  ','(n,d*39)  ','(n,d*40)  ',&
     '(n,d*41)  ','(n,d*42)  ','(n,d*43)  ','(n,d*44)  ','(n,d*45)  ',&
     '(n,d*46)  ','(n,d*47)  ','(n,d*48)  ','(n,d*c)   ',             & !300
     '(n,t*0)   ',&
     '(n,t*1)   ','(n,t*2)   ','(n,t*3)   ','(n,t*4)   ','(n,t*5)   ',&
     '(n,t*6)   ','(n,t*7)   ','(n,t*8)   ','(n,t*9)   ','(n,t*10)  ',&
     '(n,t*11)  ','(n,t*12)  ','(n,t*13)  ','(n,t*14)  ','(n,t*15)  ',&
     '(n,t*16)  ','(n,t*17)  ','(n,t*18)  ','(n,t*19)  ','(n,t*20)  ',&
     '(n,t*21)  ','(n,t*22)  ','(n,t*23)  ','(n,t*24)  ','(n,t*25)  ',&
     '(n,t*26)  ','(n,t*27)  ','(n,t*28)  ','(n,t*29)  ','(n,t*30)  ',&
     '(n,t*31)  ','(n,t*32)  ','(n,t*33)  ','(n,t*34)  ','(n,t*35)  ',&
     '(n,t*36)  ','(n,t*37)  ','(n,t*38)  ','(n,t*39)  ','(n,t*40)  ',&
     '(n,t*41)  ','(n,t*42)  ','(n,t*43)  ','(n,t*44)  ','(n,t*45)  ',&
     '(n,t*46)  ','(n,t*47)  ','(n,t*48)  ','(n,t*c)   ',             & !350
     '(n,he3*0) ',&
     '(n,he3*1) ','(n,he3*2) ','(n,he3*3) ','(n,he3*4) ','(n,he3*5) ',&
     '(n,he3*6) ','(n,he3*7) ','(n,he3*8) ','(n,he3*9) ','(n,he3*10)',&
     '(n,he3*11)','(n,he3*12)','(n,he3*13)','(n,he3*14)','(n,he3*15)',&
     '(n,he3*16)','(n,he3*17)','(n,he3*18)','(n,he3*19)','(n,he3*20)',&
     '(n,he3*21)','(n,he3*22)','(n,he3*23)','(n,he3*24)','(n,he3*25)',&
     '(n,he3*26)','(n,he3*27)','(n,he3*28)','(n,he3*29)','(n,he3*30)',&
     '(n,he3*31)','(n,he3*32)','(n,he3*33)','(n,he3*34)','(n,he3*35)',&
     '(n,he3*36)','(n,he3*37)','(n,he3*38)','(n,he3*39)','(n,he3*40)',&
     '(n,he3*41)','(n,he3*42)','(n,he3*43)','(n,he3*44)','(n,he3*45)',&
     '(n,he3*46)','(n,he3*47)','(n,he3*48)','(n,he3*c) ',             & !400
     '(n,a*0)   ',&
     '(n,a*1)   ','(n,a*2)   ','(n,a*3)   ','(n,a*4)   ','(n,a*5)   ',&
     '(n,a*6)   ','(n,a*7)   ','(n,a*8)   ','(n,a*9)   ','(n,a*10)  ',&
     '(n,a*11)  ','(n,a*12)  ','(n,a*13)  ','(n,a*14)  ','(n,a*15)  ',&
     '(n,a*16)  ','(n,a*17)  ','(n,a*18)  ','(n,a*19)  ','(n,a*20)  ',&
     '(n,a*21)  ','(n,a*22)  ','(n,a*23)  ','(n,a*24)  ','(n,a*25)  ',&
     '(n,a*26)  ','(n,a*27)  ','(n,a*28)  ','(n,a*29)  ','(n,a*30)  ',&
     '(n,a*31)  ','(n,a*32)  ','(n,a*33)  ','(n,a*34)  ','(n,a*35)  ',&
     '(n,a*36)  ','(n,a*37)  ','(n,a*38)  ','(n,a*39)  ','(n,a*40)  ',&
     '(n,a*41)  ','(n,a*42)  ','(n,a*43)  ','(n,a*44)  ','(n,a*45)  ',&
     '(n,a*46)  ','(n,a*47)  ','(n,a*48)  ','(n,a*c)   ',             & !450
     '(n,x)     ','(n,x)     ','(n,x)     ','(n,x)     ','(n,x)     ',&
     '(n,x)     ','(n,x)     ','(n,x)     ','(n,x)     ','(n,x)     ',&
     '(n,x)     ','(n,x)     ','(n,x)     ','(n,x)     ','(n,x)     ',&
     '(n,x)     ','(n,x)     ','(n,x)     ','(n,x)     ','(n,x)     ',&
     '(n,x)     ','(n,x)     ','(n,x)     ','(n,x)     ','(n,x)     ',& !475
     '(n,2n*0)  ',&
     '(n,2n*1)  ','(n,2n*2)  ','(n,2n*3)  ','(n,2n*4)  ','(n,2n*5)  ',&
     '(n,2n*6)  ','(n,2n*7)  ','(n,2n*8)  ','(n,2n*9)  ','(n,2n*10) ',&
     '(n,2n*11) ','(n,2n*12) ','(n,2n*13) ','(n,2n*14) ','(n,2n*15) ',&
     '(n,2n*c)  ',                                                    &
     '(n,x)     ','(n,x)     ','(n,x)     ','(n,x)     ','(n,x)     ',&
     '(n,x)     ','(n,x)     ','(n,x)     '/)
   character(10),parameter::h719='(n,p*c)x  '
   character(10),parameter::h739='(n,d*c)x  '
   character(10),parameter::h759='(n,t*c)x  '
   character(10),parameter::h779='(n,he3*c)x'
   character(10),parameter::h799='(n,a*c)x  '
   character(10),dimension(7),parameter::hndf9=(/'(n,xn)    ',&
     '(n,xgma)  ','(n,xp)    ','(n,xd)    ','(n,xt)    ',&
     '(n,xhe3)  ','(n,xa)    '/)
   character(10)::hndf10(1)='damage    '

   if (iverf.ge.6) then
      if (mt.ge.201.and.mt.le.207) then
         name=hndf9(mt-200)
      else if (mt.eq.444) then
         name=hndf10(1)
      else
         i=mt
         if (i.gt.999) i=i-1000*(i/1000)
         if (i.ge.600) i=i-399
         name=hndf(i)
      endif
   else
      if (mt.le.200) then
         name=hndf(mt)
      else if (mt.ge.201.and.mt.le.207) then
         name=hndf9(mt-200)
      else if (mt.eq.444) then
         name=hndf10(1)
      else if (mt.ge.700.and.mt.lt.718) then
         name=hndf(mt-499)
      else if (mt.eq.718) then
         name=hndf(250)
      else if (mt.eq.719) then
         name=h719
      else if (mt.ge.720.and.mt.lt.738) then
         name=hndf(mt-469)
      else if (mt.eq.738) then
         name=hndf(300)
      else if (mt.eq.739) then
         name=h739
      else if (mt.ge.740.and.mt.lt.758) then
         name=hndf(mt-439)
      else if (mt.eq.758) then
         name=hndf(350)
      else if (mt.eq.759) then
         name=h759
      else if (mt.ge.760.and.mt.lt.779) then
         name=hndf(mt-409)
      else if (mt.eq.778) then
         name=hndf(400)
      else if (mt.eq.779) then
         name=h779
      else if (mt.ge.780.and.mt.lt.798) then
         name=hndf(mt-379)
      else if (mt.eq.798) then
         name=hndf(450)
      else if (mt.eq.799) then
         name=h799
      endif
   endif

   !--alternate name when processing incident charged particle files
   if (izai.gt.1) then
      if (mt.eq.4) then
         name='(z,n)    '
      elseif ((izai.eq.1001.and.mt.eq.103).or.&
              (izai.eq.1002.and.mt.eq.104).or.&
              (izai.eq.1003.and.mt.eq.105).or.&
              (izai.eq.2003.and.mt.eq.106).or.&
              (izai.eq.2004.and.mt.eq.107)) then
         name=hndf(4)
      endif
   endif

   return
   end subroutine mtname

   subroutine ptleg2(a)
   !-------------------------------------------------------------------
   ! Used only for newfor=1.  Works on MF4 sections with Legendre
   ! coefficients or MF6/LAW=2 sections from File 6 with Legendre
   ! coefficients and converts them into tabulated distributions
   ! adjusted to have a cummulative sum of 1.
   !-------------------------------------------------------------------
   use util ! provides error,mess
   use mathm ! provides legndr
   use endf ! provides endf routines and variables
   ! externals
   real(kr)::a(*)
   ! internals
   integer::nord,j,negs,ii,i,idone,nn,jj,k
   integer,parameter::maxang=7000
   integer,parameter::imax=24
   integer,parameter::ipmax=65
   real(kr)::e,dy,dm,xm,ym,yt,test,xl,yl,check,dco,f,diff
   real(kr)::aco(maxang),cprob(maxang),cumm(maxang)
   real(kr)::x(imax),y(imax),p(ipmax),fl(ipmax)
   character(32)::strng
   real(kr),parameter::one=1.e0_kr
   real(kr),parameter::half=.5e0_kr
   real(kr),parameter::tenth=.1e0_kr
   real(kr),parameter::tol1=0.0002e0_kr
   real(kr),parameter::tol2=0.002e0_kr
   real(kr),parameter::pmin=1.e-10_kr
   real(kr),parameter::zero=0

   !--get the legendre coefficients
   e=a(2)
   nord=nint(a(5))
   if (nord+1.gt.ipmax) then
      write(strng,'(''nord='',i3,'' > ipmax='',i3,'', in mt '',i3)')&
                      nord,ipmax,mth
      call error('ptleg2',strng,'see ipmax')
   endif
   if (nord.ne.0) then
      do j=1,nord
         fl(j)=a(6+j)
      enddo
   else
      nord=1
      fl(1)=0
   endif

   !--use adaptive reconstruction of the angular distribution
   !--use tight tolerances to avoid misconvergences
   negs=0
   ii=0
   cumm(1)=0
   ! prime the adaptive stack
   i=2
   x(2)=-1
   y(2)=half
   call legndr(x(2),p,nord)
   do j=1,nord
      y(2)=y(2)+half*(2*j+1)*fl(j)*p(j+1)
   enddo
   x(1)=1
   y(1)=half
   call legndr(x(1),p,nord)
   do j=1,nord
      y(1)=y(1)+half*(2*j+1)*fl(j)*p(j+1)
   enddo
   ! carry out the adaptive reconstruction
   do while (i.gt.0)
      dy=0
      if (i.gt.1.and.i.lt.imax) then
         dm=x(i-1)-x(i)
         xm=half*(x(i-1)+x(i))
         xm=sigfig(xm,3,0)
         if (xm.gt.x(i).and.xm.lt.x(i-1)) then
            ym=half*(y(i-1)+y(i))
            yt=half
            call legndr(xm,p,nord)
            do j=1,nord
               yt=yt+half*(2*j+1)*fl(j)*p(j+1)
            enddo
            test=tol1*abs(yt)+one/1000000
            dy=abs(yt-ym)
            if (dm.gt.tenth) dy=2*test
            if (ym.ne.zero.and.yt/ym.gt.one*5) dy=2*test
            if (ym.ne.zero.and.yt/ym.lt.one/5) dy=2*test
         endif
      endif
      if (dy.gt.test) then
         ! not converged.
         ! add midpoint to stack and continue.
         i=i+1
         x(i)=x(i-1)
         y(i)=y(i-1)
         x(i-1)=xm
         y(i-1)=yt
      else
         ! converged.
         ! use top point in stack.
         ii=ii+1
         if (ii.gt.maxang) call error('ptleg2','too many angles',' ')
         aco(ii)=x(i)
         if (y(i).lt.pmin) then
            y(i)=pmin
            negs=negs+1
         endif
         cprob(ii)=y(i)
         if (ii.gt.1) cumm(ii)=cumm(ii-1)+half*(x(i)-xl)*(y(i)+yl)
         xl=x(i)
         yl=y(i)
         i=i-1
      endif
   enddo
   if (negs.gt.0) then
      write(strng,'(i4,'' for mt='',i3,'' e='',1p,e10.3)') negs,mth,e
      call mess('ptleg2','negative probs found',strng)
   endif
   nn=ii-1

   !--now thin the distribution to a coarser tolerance
   i=1
   ii=1
   aco(ii)=-1
   idone=0
   do while (i.lt.nn-1.and.idone.eq.0)
      check=0
      dco=0
      j=i+1
      do while (j.lt.nn+1.and.check.le.zero.and.dco.le.one)
         j=j+1
         jj=j-1
         dco=aco(j)-aco(i)
         if (dco.le.one) then
            k=i
            do while (k.lt.j-1.and.check.le.zero)
               k=k+1
               f=(aco(j)-aco(k))/dco
               test=f*cprob(i)+(1-f)*cprob(j)
               diff=one/10000000+tol2*cprob(k)
               check=abs(test-cprob(k))-diff
            enddo
         endif
      enddo
      if (check.gt.zero.or.dco.gt.one) then
         i=jj
         ii=ii+1
         aco(ii)=aco(i)
         cprob(ii)=cprob(i)
         if (ii.gt.1) cumm(ii)=cumm(ii-1)&
           +half*(aco(ii)-aco(ii-1))*(cprob(ii)+cprob(ii-1))
      else
         idone=1
      endif
   enddo
   i=nn+1
   ii=ii+1
   aco(ii)=aco(i)
   cprob(ii)=cprob(i)
   cumm(ii)=cumm(ii-1)+half*(aco(ii)-aco(ii-1))*(cprob(ii)+cprob(ii-1))

   !--load the new distribution in the output list
   a(5)=1
   a(6)=ii
   a(7)=ii
   a(8)=2
   do i=1,ii
      a(7+2*i)=aco(i)
      a(8+2*i)=cprob(i)/cumm(ii)
   enddo
   return
   end subroutine ptleg2

   subroutine pttab2(a)
   !-------------------------------------------------------------------
   ! Like pttab, but works on tabulated sections from File 4 or
   ! LAW2 in File 6 and checks the cummulative normalization.
   !-------------------------------------------------------------------
   ! externals
   real(kr)::a(*)
   ! internals
   integer::nr,np,i
   real(kr)::cumm
   real(kr),allocatable::amu(:),pmu(:)

   !--check normalization
   nr=nint(a(5))
   np=nint(a(6))
   allocate(amu(np))
   allocate(pmu(np))
   cumm=0
   do i=1,np
      amu(i)=a(5+2*nr+2*i)
      pmu(i)=a(6+2*nr+2*i)
      if (i.gt.1) then
         cumm=cumm+(amu(i)-amu(i-1))*(pmu(i)+pmu(i-1))/2
      endif
   enddo

   !--write out the modified record
   a(5)=1
   a(6)=np
   a(7)=np
   a(8)=2
   do i=1,np
      a(7+2*i)=amu(i)
      a(8+2*i)=pmu(i)/cumm
   enddo
   deallocate(amu)
   deallocate(pmu)
   return
   end subroutine pttab2

   real(kr) function bachaa(iza1i,iza2,izat,e,ep)
   !-------------------------------------------------------------------
   ! Compute the Kalbach a parameter.
   !-------------------------------------------------------------------
   use physics ! provides amassn,amu,ev,clight
   use util ! provides error
   ! externals
   integer::iza1i,iza2,izat
   real(kr)::e,ep
   ! internals
   integer::iza1,iza,na
   real(kr)::emc2,aa,za,ac,zc,ab,zb,sa,sb,ecm,ea,eb,x1,x3
   real(kr)::fa,fb,bb,fact,test,nc,nb
   character(60)::strng
   real(kr),parameter::third=.333333333e0_kr
   real(kr),parameter::twoth=.666666667e0_kr
   real(kr),parameter::fourth=1.33333333e0_kr
   real(kr),parameter::c1=15.68e0_kr
   real(kr),parameter::c2=-28.07e0_kr
   real(kr),parameter::c3=-18.56e0_kr
   real(kr),parameter::c4=33.22e0_kr
   real(kr),parameter::c5=-0.717e0_kr
   real(kr),parameter::c6=1.211e0_kr
   real(kr),parameter::s2=2.22e0_kr
   real(kr),parameter::s3=8.48e0_kr
   real(kr),parameter::s4=7.72e0_kr
   real(kr),parameter::s5=28.3e0_kr
   real(kr),parameter::b1=0.04e0_kr
   real(kr),parameter::b2=1.8e-6_kr
   real(kr),parameter::b3=6.7e-7_kr
   real(kr),parameter::d1=9.3e0_kr
   real(kr),parameter::ea1=41.e0_kr
   real(kr),parameter::ea2=130.e0_kr
   real(kr),parameter::emev=1.e6_kr

   emc2=amassn*amu*clight*clight/ev/emev

   iza1=iza1i
   if (iza1i.eq.0) iza1=1
   iza=izat
   if (iza.eq.6000) iza=6012
   if (iza.eq.12000) iza=12024
   if (iza.eq.14000) iza=14028
   if (iza.eq.16000) iza=16032
   if (iza.eq.17000) iza=17035
   if (iza.eq.19000) iza=19039
   if (iza.eq.20000) iza=20040
   if (iza.eq.22000) iza=22048
   if (iza.eq.23000) iza=23051
   if (iza.eq.24000) iza=24052
   if (iza.eq.26000) iza=26056
   if (iza.eq.28000) iza=28058
   if (iza.eq.29000) iza=29063
   if (iza.eq.31000) iza=31069
   if (iza.eq.40000) iza=40090
   if (iza.eq.42000) iza=42096
   if (iza.eq.48000) iza=48112
   if (iza.eq.49000) iza=49115
   if (iza.eq.50000) iza=50120
   if (iza.eq.63000) iza=63151
   if (iza.eq.72000) iza=72178
   if (iza.eq.74000) iza=74184
   if (iza.eq.82000) iza=82208
   aa=mod(iza,1000)
   if (aa.eq.0.) then
      write(strng,'(''dominant isotope not known for '',i8)') iza
      call error('bachaa',strng,' ')
   endif
   za=int(iza/1000)
   ac=aa+mod(iza1,1000)
   zc=za+int(iza1/1000)
   ab=ac-mod(iza2,1000)
   zb=zc-int(iza2/1000)
   na=nint(aa-za)
   nb=nint(ab-zb)
   nc=nint(ac-zc)
   sa=c1*(ac-aa)&
     +c2*((nc-zc)**2/ac-(na-za)**2/aa)&
     +c3*(ac**twoth-aa**twoth)&
     +c4*((nc-zc)**2/ac**fourth-(na-za)**2/aa**fourth)&
     +c5*(zc**2/ac**third-za**2/aa**third)&
     +c6*(zc**2/ac-za**2/aa)
   if (iza1.eq.1002) sa=sa-s2
   if (iza1.eq.1003) sa=sa-s3
   if (iza1.eq.2003) sa=sa-s4
   if (iza1.eq.2004) sa=sa-s5
   sb=c1*(ac-ab)&
     +c2*((nc-zc)**2/ac-(nb-zb)**2/ab)&
     +c3*(ac**twoth-ab**twoth)&
     +c4*((nc-zc)**2/ac**fourth-(nb-zb)**2/ab**fourth)&
     +c5*(zc**2/ac**third-zb**2/ab**third)&
     +c6*(zc**2/ac-zb**2/ab)
   if (iza2.eq.1002) sb=sb-s2
   if (iza2.eq.1003) sb=sb-s3
   if (iza2.eq.2003) sb=sb-s4
   if (iza2.eq.2004) sb=sb-s5
   ecm=aa*e/ac
   ea=ecm+sa
   eb=ep*ac/ab+sb
   x1=eb
   if (ea.gt.ea2) x1=ea2*eb/ea
   x3=eb
   if (ea.gt.ea1) x3=ea1*eb/ea
   fa=1
   if (iza1.eq.2004) fa=0
   fb=1
   if (iza2.eq.1) fb=fb/2
   if (iza2.eq.2004) fb=2
   bb=b1*x1+b2*x1**3+b3*fa*fb*x3**4
   if (iza1i.eq.0) then
      fact=d1
      if (ep.ne.0.) fact=fact/sqrt(ep)
      test=1
      if (fact.lt.test) fact=test
      test=4
      if (fact.gt.test) fact=test
      bb=bb*sqrt(e/(2*emc2))*fact
   endif
   bachaa=bb
   return
   end function bachaa

   real(kr) function eavl(akal,amass,avcm,avadd,fmsd,sign)
   !-------------------------------------------------------------------
   ! Analytically calculates the average energy (lab) assuming
   ! Kalbach systematics for c.m. angular distribution.
   ! Note: sign=1. is for particles. sign=-1 = for recoils
   ! since theta-r=pi-theta-p, cos(theta-r)=-cos(theta-p), and
   ! Kalbach exponentials just have sign flipped for recoils.
   !--***************************************************************
   use physics ! provides pi
   ! externals
   real(kr)::akal,amass,avcm,avadd,fmsd,sign
   ! internals
   real(kr)::r,s

   ! writing the kalbach function in the following form
   ! p = f1 exp(a cos theta) + f2 exp(-a cos theta)
   ! f1=(1/(4*pi))*(akal/(exp(akal)-exp(-akal)))*(1+fmsd)
   ! f2=(1/(4*pi))*(akal/(exp(akal)-exp(-akal)))*(1-fmsd)
   ! if (akal.eq.0.) f1=(1/(4*pi))*(1+fmsd)/2
   ! if (akal.eq.0.) f2=(1/(4*pi))*(1-fmsd)/2
   ! the lab average energy is given by
   r=avcm*avcm+avadd*avadd
   s=2*avcm*avadd
   eavl=amass*(r/2+(fmsd*s*fi2(sign*akal)/(fi1(sign*akal)*2)))
   return
   end function eavl

   real(kr) function fi1(aa)
   !-------------------------------------------------------------------
   ! Solution of first type of integral for average energy.
   !-------------------------------------------------------------------
   ! externals
   real(kr)::aa
   ! internals
   real(kr),parameter::zero=0

   if (aa.eq.zero) then
      fi1=2
   else
      fi1=(exp(aa)-exp(-aa))/aa
   endif
   return
   end function fi1

   real(kr) function fi2(aa)
   !-------------------------------------------------------------------
   ! Solution of second type of integral for average energy.
   !-------------------------------------------------------------------
   ! externals
   real(kr)::aa
   ! internals
   real(kr),parameter::zero=0

   if (aa.eq.zero) then
      fi2=0
   else
      fi2=((exp(aa)+exp(-aa))/aa)-((exp(aa)-exp(-aa))/(aa*aa))
   endif
   return
   end function fi2

end module acecm


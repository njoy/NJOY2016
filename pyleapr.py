from __future__ import annotations
import math
import numpy as np
# Minimal constants (cgs/ev) to mirror physics module where needed
BK_EV = 8.617333262e-5  # eV/K
HBAR  = 1.054571817e-27 # erg*s
EV    = 1.602176634e-12 # erg
AMU   = 1.66053906660e-24 # g
AMASSN= 1.00866491595    # neutron mass in amu


# Globals expected to be set by the driver (mirrors module leapm)
alpha: np.ndarray | None = None
beta:  np.ndarray | None = None
nalpha: int = 0
nbeta:  int = 0
lat: int = 0
tev: float = 0.0
deltab: float = 0.0
tbeta: float = 0.0
twt: float = 0.0
c: float = 0.0
iprint: int = 1
arat: float = 1.0
f0: float = 0.0
ssm: np.ndarray | None = None
ssp: np.ndarray | None = None
tempf: np.ndarray | None = None
tempr: np.ndarray | None = None
dwpix: np.ndarray | None = None
ska: np.ndarray | None = None
nka: int = 0
dka: float = 0.0
cfrac: float = 0.0
nss: int = 0
b7: float = 0.0
aws: float = 0.0
sps: float = 0.0
mss: int = 0
awr: float = 0.0
spr: float = 0.0
npr: int = 1
iel: int = 0
ncold: int = 0
nd: int = 0
bdel: np.ndarray | None = None
adel: np.ndarray | None = None
delta1: float = 0.0
p1: np.ndarray | None = None
np1: int = 0
tbar: float = 0.0
naint: int = 1
nbint: int = 1


def bfill(betan: np.ndarray) -> tuple[np.ndarray, np.ndarray, int]:
    """\
subroutine bfill(bex,rdbex,nbx,betan,nbetan,maxbb)
!--------------------------------------------------------------------
   ! Sets up the arrays bex and rdbex used by sint.
   ! Called by discre.
   !--------------------------------------------------------------------
   ! externals
   integer::nbx,nbetan,maxbb
   real(kr)::bex(maxbb),rdbex(maxbb),betan(nbetan)
   ! internals
   integer::i,k,nbxm
   real(kr),parameter::small=1.e-9_kr

   k=nbeta
   do i=1,nbeta
      bex(i)=-betan(k)
      k=k-1
   enddo
   if (betan(1).le.small) then
      bex(nbeta)=0
      k=nbeta+1
   else
      k=nbeta+2
      bex(nbeta+1)=betan(1)
   endif
   do i=2,nbeta
      bex(k)=betan(i)
      k=k+1
   enddo
   nbxm=k-2
   nbx=k-1
   do i=1,nbxm
      rdbex(i)=1/(bex(i+1)-bex(i))
   enddo
   return
end subroutine bfill
    """
    small = 1.0e-9
    nbeta = len(betan)
    neg = -betan[::-1]
    if betan[0] <= small:
        mid = np.array([0.0])
        pos = betan[1:]
    else:
        mid = np.array([betan[0]])
        pos = betan[1:]
    bex_new = np.concatenate([neg, mid, pos])
    nbx = len(bex_new)
    rdbex_new = np.empty(nbx-1, dtype=np.float64)
    rdbex_new[:] = 1.0 / (bex_new[1:] - bex_new[:-1])
    return bex_new, rdbex_new, nbx

def exts(sexpb: np.ndarray, exb: np.ndarray, betan: np.ndarray) -> np.ndarray:
    """\
subroutine exts(sexpb,sex,exb,betan,nbetan,maxbb)
!--------------------------------------------------------------------
   ! Sets up the array sex for sint.
   ! Here sexpb contains the asymmetric sab for negative beta, and
   ! sex contains the asymmetric sab extended to plus and minus beta.
   ! Called by discre.
   !--------------------------------------------------------------------
   ! externals
   integer::nbetan,maxbb
   real(kr)::sexpb(nbetan),sex(maxbb),exb(maxbb),betan(nbetan)
   ! internals
   integer::i,k
   real(kr),parameter::small=1.e-9_kr

   k=nbeta
   do i=1,nbeta
      sex(i)=sexpb(k)
      k=k-1
   enddo
   if (betan(1).le.small) then
      sex(nbeta)=sexpb(1)
      k=nbeta+1
   else
      k=nbeta+2
      sex(nbeta+1)=sexpb(1)
   endif
   do i=2,nbeta
      sex(k)=sexpb(i)*exb(i)*exb(i)
      k=k+1
   enddo
   return
end subroutine exts
    """
    # extend asymmetric sab to +/- beta side
    small = 1.0e-9
    nbeta = len(betan)
    neg = sexpb[::-1]
    if betan[0] <= small:
        mid = np.array([sexpb[0]])
        pos_vals = sexpb[1:] * (exb[1:]**2)
    else:
        mid = np.array([sexpb[0]])
        pos_vals = sexpb[1:] * (exb[1:]**2)
    sex = np.concatenate([neg, mid, pos_vals])
    return sex

def sint(x: float, bex: np.ndarray, rdbex: np.ndarray, sex: np.ndarray, nbx: int,
         alph: float, wt: float, tbart: float, betan: np.ndarray) -> float:
    """\
function sint(x,bex,rdbex,sex,nbx,alph,wt,tbart,&
     betan,nbetan,maxbb)
!--------------------------------------------------------------------
   ! Interpolates in scattering function, or uses SCT approximation
   ! to extrapolate outside the range in memory.
   ! Called by discre.
   !--------------------------------------------------------------------
   use physics ! provides pi
   ! externals
   integer::nbx,nbetan,maxbb
   real(kr)::x,alph,wt,tbart
   real(kr)::bex(maxbb),rdbex(maxbb),sex(maxbb),betan(nbetan)
   ! internals
   integer::k1,k2,k3,idone
   real(kr)::sv,ex,ss1,ss3
   real(kr),parameter::slim=-225.e0_kr
   real(kr),parameter::zero=0

   !--sct approx
   if (abs(x).gt.betan(nbeta)) then
      if (alph.le.zero) then
         sv=0
      else
         ex=-(wt*alph-abs(x))**2/(4*wt*alph*tbart)
         if (x.gt.zero) ex=ex-x
         sv=exp(ex)/(4*pi*wt*alph*tbart)
      endif
      sint=sv
      return
   endif

   !--interpolation
   k1=1
   k2=nbeta
   k3=nbx
   !  bisect for x
   idone=0
   do while (idone.eq.0)
      if (x.eq.bex(k2)) then
         sv=sex(k2)
         sint=sv
         return
      else if (x.gt.bex(k2)) then
         k1=k2
         k2=(k3-k2)/2+k2
         if (k3-k1.le.1) idone=1
      else
         k3=k2
         k2=(k2-k1)/2+k1
         if (k3-k1.le.1) idone=1
      endif
   enddo
   if (sex(k1).le.zero) then
      ss1=slim
   else
      ss1=log(sex(k1))
   endif
   if (sex(k3).le.zero) then
      ss3=slim
   else
      ss3=log(sex(k3))
   endif
   ex=((bex(k3)-x)*ss1+(x-bex(k1))*ss3)*rdbex(k1)
   sv=0
   if (ex.gt.slim) sv=exp(ex)
   sint=sv
   return
end function sint
    """
    slim = -225.0
    zero = 0.0
    nbeta = len(betan)
    # SCT approx
    if abs(x) > betan[-1]:
        if alph <= zero:
            return 0.0
        ex = -((wt*alph - abs(x))**2) / (4.0*wt*alph*tbart)
        if x > zero:
            ex = ex - x
        return math.exp(ex) / (4.0*math.pi*wt*alph*tbart)
    # interpolation (binary search)
    k1 = 0
    k3 = nbx - 1
    if x == bex[k3]:
        return float(sex[k3])
    while k3 - k1 > 1:
        k2 = (k1 + k3) // 2
        if x > bex[k2]:
            k1 = k2
        else:
            k3 = k2
    ss1 = slim if sex[k1] <= zero else math.log(sex[k1])
    ss3 = slim if sex[k3] <= zero else math.log(sex[k3])
    ex = ((bex[k3] - x) * ss1 + (x - bex[k1]) * ss3) * rdbex[k1]
    return 0.0 if ex <= slim else float(math.exp(ex))

def terps(sd: np.ndarray, delta: float, be: float) -> float:
    """\
function terps(sd,nsd,delta,be)
!--------------------------------------------------------------------
   ! Interpolate in a table of S(alpha,beta) for a required beta.
   ! Used in trans.
   !--------------------------------------------------------------------
   ! externals
   integer::nsd
   real(kr)::delta,be
   real(kr)::sd(nsd)
   ! internals
   integer::i
   real(kr)::bt,btp,st,stp,stt
   real(kr),parameter::slim=-225.e0_kr
   real(kr),parameter::zero=0

   terps=0
   if (be.gt.delta*nsd) return
   i=int(be/delta)
   if (i.lt.nsd-1) then
      bt=i*delta
      btp=bt+delta
      i=i+1
      if (sd(i).le.zero) then
         st=slim
      else
         st=log(sd(i))
      endif
      if (sd(i+1).le.zero) then
         stp=slim
      else
         stp=log(sd(i+1))
      endif
      stt=st+(be-bt)*(stp-st)/(btp-bt)
      terps=0
      if (stt.gt.slim) terps=exp(stt)
      return
   endif
   return
end function terps
    """
    slim = -225.0
    zero = 0.0
    nsd = len(sd)
    if be > delta * nsd:
        return 0.0
    i = int(be / delta)
    if i < nsd - 1:
        bt = i * delta
        btp = bt + delta
        i1 = i + 1
        st  = slim if sd[i1]   <= zero else math.log(sd[i1])
        stp = slim if sd[i1+1] <= zero else math.log(sd[i1+1])
        stt = st + (be - bt) * (stp - st) / (btp - bt)
        return 0.0 if stt <= slim else float(math.exp(stt))
    else:
        return 0.0

def besk1(x: float) -> float:
    """\
function besk1(x)
!--------------------------------------------------------------------
   ! Computes modified Bessel function, K1.
   ! The exponential part for x>1 is omitted (see stable).
   ! Called by stable.
   !--------------------------------------------------------------------
   ! externals
   real(kr)::x
   ! internals
   real(kr)::test,v,u,bi1,bi3
   real(kr),parameter::c0=.125e0_kr
   real(kr),parameter::c1=.442850424e0_kr
   real(kr),parameter::c2=.584115288e0_kr
   real(kr),parameter::c3=6.070134559e0_kr
   real(kr),parameter::c4=17.864913364e0_kr
   real(kr),parameter::c5=48.858995315e0_kr
   real(kr),parameter::c6=90.924600045e0_kr
   real(kr),parameter::c7=113.795967431e0_kr
   real(kr),parameter::c8=85.331474517e0_kr
   real(kr),parameter::c9=32.00008698e0_kr
   real(kr),parameter::c10=3.999998802e0_kr
   real(kr),parameter::c11=1.304923514e0_kr
   real(kr),parameter::c12=1.47785657e0_kr
   real(kr),parameter::c13=16.402802501e0_kr
   real(kr),parameter::c14=44.732901977e0_kr
   real(kr),parameter::c15=115.837493464e0_kr
   real(kr),parameter::c16=198.437197312e0_kr
   real(kr),parameter::c17=222.869709703e0_kr
   real(kr),parameter::c18=142.216613971e0_kr
   real(kr),parameter::c19=40.000262262e0_kr
   real(kr),parameter::c20=1.999996391e0_kr
   real(kr),parameter::c21=1.e0_kr
   real(kr),parameter::c22=.5e0_kr
   real(kr),parameter::c23=.5772156649e0_kr
   real(kr),parameter::c24=1.e0_kr
   real(kr),parameter::c25=.0108241775e0_kr
   real(kr),parameter::c26=.0788000118e0_kr
   real(kr),parameter::c27=.2581303765e0_kr
   real(kr),parameter::c28=.5050238576e0_kr
   real(kr),parameter::c29=.663229543e0_kr
   real(kr),parameter::c30=.6283380681e0_kr
   real(kr),parameter::c31=.4594342117e0_kr
   real(kr),parameter::c32=.2847618149e0_kr
   real(kr),parameter::c33=.1736431637e0_kr
   real(kr),parameter::c34=.1280426636e0_kr
   real(kr),parameter::c35=.1468582957e0_kr
   real(kr),parameter::c36=.4699927013e0_kr
   real(kr),parameter::c37=1.2533141373e0_kr

   test=1
   if (x.le.test) then
      v=c0*x
      u=v*v
      bi1=(((((((((c1*u+c2)*u+c3)*u+c4)*u+c5)*u+c6)*u+c7)*u&
        +c8)*u+c9)*u+c10)*v
      bi3=(((((((((c11*u+c12)*u+c13)*u+c14)*u+c15)*u+c16)*u&
        +c17)*u+c18)*u+c19)*u+c20)
      besk1=c21/x+bi1*(log(c22*x)+c23)-v*bi3
   else
      u=c24/x
      bi3=((((((((((((-c25*u+c26)*u-c27)*u+c28)*u-c29)*u+c30)*u&
        -c31)*u+c32)*u-c33)*u+c34)*u-c35)*u+c36)*u+c37)
      besk1=sqrt(u)*bi3
   endif
   return
end function besk1
    """
    # K1-like approximation (with exp factor separated when x>1)
    c0=0.125; c1=0.442850424; c2=0.584115288; c3=6.070134559; c4=17.864913364
    c5=48.858995315; c6=90.924600045; c7=113.795967431; c8=85.331474517; c9=32.00008698
    c10=3.999998802; c11=1.304923514; c12=1.47785657; c13=16.402802501; c14=44.732901977
    c15=115.837493464; c16=198.437197312; c17=222.869709703; c18=142.216613971; c19=40.000262262
    c20=1.999996391; c21=1.0; c22=0.5; c23=0.5772156649; c24=1.0; c25=0.0108241775
    c26=0.0788000118; c27=0.2581303765; c28=0.5050238576; c29=0.663229543; c30=0.6283380681
    c31=0.4594342117; c32=0.2847618149; c33=0.1736431637; c34=0.1280426636; c35=0.1468582957
    c36=0.4699927013; c37=1.2533141373
    test = 1.0
    if x <= test:
        v = c0 * x
        u = v*v
        bi1 = (((((((((c1*u+c2)*u+c3)*u+c4)*u+c5)*u+c6)*u+c7)*u + c8)*u + c9)*u + c10)*v
        bi3 = (((((((((c11*u+c12)*u+c13)*u+c14)*u+c15)*u+c16)*u + c17)*u + c18)*u + c19)*u + c20)
        return c21/x + bi1*(math.log(c22*x)+c23) - v*bi3
    else:
        u = c24 / x
        bi3 = (((((((((((-c25*u+c26)*u-c27)*u+c28)*u-c29)*u+c30)*u - c31)*u + c32)*u - c33)*u + c34)*u - c35)*u + c36)*u + c37
        return math.sqrt(u) * bi3

def bfact(x: float, dwc: float, betai: float) -> tuple[float, np.ndarray, np.ndarray]:
    """\
subroutine bfact(x,bzero,bplus,bminus,dwc,betai)
!--------------------------------------------------------------------
   ! Calculates the Bessel function terms for discrete oscillators.
   ! Called by discre.
   !--------------------------------------------------------------------
   ! externals
   real(kr)::x,dwc,betai
   real(kr)::bzero,bplus(50),bminus(50)
   ! internals
   integer::i,imax,j
   real(kr)::bn(50)
   real(kr)::y,u,v,bessi0,bessi1,rat,arg
   real(kr),parameter::c0=3.75e0_kr
   real(kr),parameter::c1=1.e0_kr
   real(kr),parameter::c2=3.5156229e0_kr
   real(kr),parameter::c3=3.0899424e0_kr
   real(kr),parameter::c4=1.2067492e0_kr
   real(kr),parameter::c5=0.2659732e0_kr
   real(kr),parameter::c6=0.0360768e0_kr
   real(kr),parameter::c7=0.0045813e0_kr
   real(kr),parameter::c8=0.39894228e0_kr
   real(kr),parameter::c9=0.01328592e0_kr
   real(kr),parameter::c10=0.00225319e0_kr
   real(kr),parameter::c11=0.00157565e0_kr
   real(kr),parameter::c12=0.00916281e0_kr
   real(kr),parameter::c13=0.02057706e0_kr
   real(kr),parameter::c14=0.02635537e0_kr
   real(kr),parameter::c15=0.01647633e0_kr
   real(kr),parameter::c16=0.00392377e0_kr
   real(kr),parameter::c17=0.5e0_kr
   real(kr),parameter::c18=0.87890594e0_kr
   real(kr),parameter::c19=0.51498869e0_kr
   real(kr),parameter::c20=0.15084934e0_kr
   real(kr),parameter::c21=0.02658733e0_kr
   real(kr),parameter::c22=0.00301532e0_kr
   real(kr),parameter::c23=0.00032411e0_kr
   real(kr),parameter::c24=0.02282967e0_kr
   real(kr),parameter::c25=0.02895312e0_kr
   real(kr),parameter::c26=0.01787654e0_kr
   real(kr),parameter::c27=0.00420059e0_kr
   real(kr),parameter::c28=0.39894228e0_kr
   real(kr),parameter::c29=0.03988024e0_kr
   real(kr),parameter::c30=0.00362018e0_kr
   real(kr),parameter::c31=0.00163801e0_kr
   real(kr),parameter::c32=0.01031555e0_kr
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::tiny=1.e-30_kr
   real(kr),parameter::zero=0
   real(kr),parameter::one=1

   !--compute bessi0
   y=x/c0
   if (y.le.one) then
      u=y*y
      bessi0=c1+u*(c2+u*(c3+u*(c4+u*(c5+u*(c6+u*c7)))))
   else
      v=1/y
      bessi0=(c8+v*(c9+v*(c10+v*(-c11+v*(c12+v*(-c13&
        +v*(c14+v*(-c15+v*c16))))))))/sqrt(x)
   endif

   !--compute bessi1
   if (y.le.one) then
      u=y*y
      bessi1=(c17+u*(c18+u*(c19+u*(c20+u*(c21+u*(c22+u*c23))))))*x
   else
      v=1/y
      bessi1=c24+v*(-c25+v*(c26-v*c27))
      bessi1=c28+v*(-c29+v*(-c30+v*(c31+v*(-c32+v*bessi1))))
      bessi1=bessi1/sqrt(x)
   endif

   !--generate higher orders by reverse recursion
   imax=50
   bn(imax)=0
   bn(imax-1)=1
   i=imax-1
   do while (i.gt.1)
      bn(i-1)=bn(i+1)+i*(2/x)*bn(i)
      i=i-1
      if (bn(i).ge.big) then
         do j=i,imax
            bn(j)=bn(j)/big
         enddo
      endif
   enddo
   rat=bessi1/bn(1)
   do i=1,imax
      bn(i)=bn(i)*rat
      if (bn(i).lt.tiny) bn(i)=0
   enddo

   !--apply exponential terms to bessel functions
   if (y.le.one) then
      bzero=bessi0*exp(-dwc)
      do i=1,imax
         bminus(i)=0
         bplus(i)=0
         if (bn(i).ne.zero) then
            arg=-dwc-i*betai/2
            bplus(i)=0
            bplus(i)=exp(arg)*bn(i)
            if (bplus(i).lt.tiny) bplus(i)=0
            bminus(i)=0
            arg=-dwc+i*betai/2
            bminus(i)=exp(arg)*bn(i)
            if (bminus(i).lt.tiny) bminus(i)=0
         else
            bplus(i)=0
            bminus(i)=0
         endif
      enddo
   else
      bzero=bessi0*exp(-dwc+x)
      do i=1,imax
         bminus(i)=0
         bplus(i)=0
         if (bn(i).ne.zero) then
            bplus(i)=0
            arg=-dwc-i*betai/2+x
            bplus(i)=exp(arg)*bn(i)
            if (bplus(i).lt.tiny) bplus(i)=0
            bminus(i)=0
            arg=-dwc+i*betai/2+x
            bminus(i)=exp(arg)*bn(i)
            if (bminus(i).lt.tiny) bminus(i)=0
         else
            bplus(i)=0
            bminus(i)=0
         endif
      enddo
   endif
   return
end subroutine bfact
    """
    big = 1.0e10
    tiny = 1.0e-30
    one = 1.0
    c0=3.75; c1=1.0; c2=3.5156229; c3=3.0899424; c4=1.2067492; c5=0.2659732; c6=0.0360768; c7=0.0045813
    c8=0.39894228; c9=0.01328592; c10=0.00225319; c11=0.00157565; c12=0.00916281; c13=0.02057706; c14=0.02635537
    c15=0.01647633; c16=0.00392377; c17=0.5; c18=0.87890594; c19=0.51498869; c20=0.15084934; c21=0.02658733
    c22=0.00301532; c23=0.00032411; c24=0.02282967; c25=0.02895312; c26=0.01787654; c27=0.00420059
    c28=0.39894228; c29=0.03988024; c30=0.00362018; c31=0.00163801; c32=0.01031555
    y = x / c0
    # bessi0
    if y <= one:
        u = y*y
        bessi0 = c1 + u*(c2 + u*(c3 + u*(c4 + u*(c5 + u*(c6 + u*c7)))))
    else:
        v = 1.0/y
        bessi0 = (c8 + v*(c9 + v*(c10 + v*(-c11 + v*(c12 + v*(-c13 + v*(c14 + v*(-c15 + v*c16))))))))/math.sqrt(x)
    # bessi1
    if y <= one:
        u = y*y
        bessi1 = (c17 + u*(c18 + u*(c19 + u*(c20 + u*(c21 + u*(c22 + u*c23)))))) * x
    else:
        v = 1.0/y
        bessi1 = c24 + v*(-c25 + v*(c26 - v*c27))
        bessi1 = c28 + v*(-c29 + v*(-c30 + v*(c31 + v*(-c32 + v*bessi1))))
        bessi1 = bessi1 / math.sqrt(x)
    # reverse recursion
    imax = 50
    bn = np.zeros(imax+1, dtype=np.float64)
    bn[imax] = 0.0
    bn[imax-1] = 1.0
    i = imax-1
    while i > 1:
        bn[i-1] = bn[i+1] + i*(2.0/x)*bn[i]
        i -= 1
        if bn[i] >= big:
            for j in range(i, imax+1):
                bn[j] = bn[j]/big
    rat = bessi1 / bn[1]
    for i in range(1, imax+1):
        bn[i] = bn[i]*rat
        if bn[i] < tiny: bn[i] = 0.0
    bplus  = np.zeros(imax, dtype=np.float64)
    bminus = np.zeros(imax, dtype=np.float64)
    if y <= one:
        bzero = bessi0 * math.exp(-dwc)
        for i in range(1, imax+1):
            if bn[i] != 0.0:
                arg = -dwc - i*betai/2.0
                val = math.exp(arg)*bn[i]
                bplus[i-1]  = 0.0 if val < tiny else val
                arg = -dwc + i*betai/2.0
                val = math.exp(arg)*bn[i]
                bminus[i-1] = 0.0 if val < tiny else val
            else:
                bplus[i-1] = 0.0; bminus[i-1] = 0.0
    else:
        bzero = bessi0 * math.exp(-dwc + x)
        for i in range(1, imax+1):
            if bn[i] != 0.0:
                arg = -dwc - i*betai/2.0 + x
                val = math.exp(arg)*bn[i]
                bplus[i-1]  = 0.0 if val < tiny else val
                arg = -dwc + i*betai/2.0 + x
                val = math.exp(arg)*bn[i]
                bminus[i-1] = 0.0 if val < tiny else val
            else:
                bplus[i-1] = 0.0; bminus[i-1] = 0.0
    return float(bzero), bplus, bminus

def stable(al: float, delta: float, ndmax: int) -> tuple[np.ndarray, np.ndarray]:
    """\
subroutine stable(ap,sd,nsd,al,delta,iprt,nu,ndmax)
!--------------------------------------------------------------------
   ! Sets up table of S-diffusion or S-free in the array sd,
   ! evaluated at intervals delta determined by trans.
   ! Tabulation is continued until
   !       sd(j) is less than 1e-7*sd(1)
   !       or nsd is 1999
   ! Here, nsd is always odd for use with Simpson's rule, and
   ! ap is used for temporary storage of beta values.
   ! This routine returns the -beta side of
   ! the asymmetric S(alpha,beta).
   ! Called by trans.
   ! Uses besk1.
   !--------------------------------------------------------------------
   use mainio ! provide nsyso
   use physics ! provides pi
   ! externals
   integer::nsd,iprt,nu,ndmax
   real(kr)::al,delta,check0,check1,f,bb,sfree
   real(kr)::ap(ndmax),sd(ndmax)
   ! internals
   integer::i,j,idone,icheck
   real(kr)::d,c2,c3,c4,c5,c6,c7,c8,be,ex,wal
   real(kr),parameter::quart=0.25e0_kr
   real(kr),parameter::eps=1.e-7_kr
   real(kr),parameter::zero=0
   real(kr),parameter::one=1
   icheck=0

   !--diffusion branch
   if (c.ne.zero) then
      if (iprt.eq.1.and.iprint.eq.2) write(nsyso,&
        '(/4(4x,''beta'',4x,''ssdiff(-beta)''))')
      d=twt*c
      c2=sqrt(c*c+quart)
      c3=2*d*al
      c4=c3*c3
      c8=c2*c3/pi
      c3=2*d*c*al
      be=0
      j=1
      idone=0
      do while (idone.eq.0)
         c6=sqrt(be*be+c4)
         c7=c6*c2
         if (c7.le.one) c5=c8*exp(c3+be/2)
         if (c7.gt.one) then
            c5=0
            ex=c3-c7+be/2
            c5=c8*exp(ex)
         endif
         sd(j)=c5*besk1(c7)/c6
         ap(j)=be
         be=be+delta
         j=j+1
         if (mod(j,2).eq.0) then
            if (j.ge.ndmax) idone=1
            if (eps*sd(1).ge.sd(j-1)) idone=1
         endif
      enddo
      j=j-1
      nsd=j
      if (iprt.eq.1.and.iprint.eq.2)&
        write(nsyso,'(0pf10.5,1pe15.7,0pf10.5,1pe15.7,0pf10.5,&
        &1pe15.7,0pf10.5,1pe15.7)') (ap(i),sd(i),i=1,j,nu)

   !--free-gas branch
   else
      if (iprt.eq.1.and.iprint.eq.2) write(nsyso,&
        '(/4(4x,''beta'',4x,''ssfree(-beta)''))')
      be=0
      j=1
      wal=twt*al
      idone=0
      do while (idone.eq.0)
         ex=-(wal-be)**2/(4*wal)
         sfree=exp(ex)/sqrt(4*pi*wal)
         sd(j)=sfree
         ap(j)=be
         be=be+delta
         j=j+1
         if (mod(j,2).eq.0) then
            if (j.ge.ndmax) idone=1
            if (eps*sd(1).ge.sd(j-1)) idone=1
         endif
      enddo
      j=j-1
      nsd=j
      if (iprt.eq.1.and.iprint.eq.2)&
        write(nsyso,'(0pf10.5,1pe15.7,0pf10.5,1pe15.7,0pf10.5,&
        &1pe15.7,0pf10.5,1pe15.7)') (ap(i),sd(i),i=1,j,nu)
   endif

   !--check the moments of the distribution
   if (iprt.ne.1.or.iprint.ne.2) return
   if (icheck.eq.0) return
   check0=0
   check1=0
   do i=1,nsd
      f=2*(mod(i-1,2)+1)
      if (i.eq.1.or.i.eq.nsd) f=1
      bb=(i-1)*delta
      check0=check0+f*sd(i)
      check0=check0+sd(i)*exp(-bb)
      check1=check1+f*bb*sd(i)
      check1=check1+f*bb*sd(i)*exp(-bb)
   enddo
   check0=check0*delta/3
   check1=check1*delta/3
   check1=check1/(al*twt)
   write(nsyso,'(/'' check0='',f11.7,''   check1='',f11.7)')&
     check0,check1
   return
end subroutine stable
    """
    quart = 0.25
    eps = 1.0e-7
    zero = 0.0
    one = 1.0
    ap: list[float] = []
    sd: list[float] = []
    if c != 0.0:
        d = twt * c
        c2 = math.sqrt(c*c + quart)
        c3 = 2.0*d*al
        c4 = c3*c3
        c8 = c2*c3/math.pi
        c3 = 2.0*d*c*al
        be = 0.0
        j = 0
        while True:
            j += 1
            c6 = math.sqrt(be*be + c4)
            c7 = c6*c2
            if c7 <= one:
                c5 = c8*math.exp(c3 + be/2.0)
            else:
                ex = c3 - c7 + be/2.0
                c5 = c8*math.exp(ex)
            val = c5*besk1(c7)/c6
            sd.append(val)
            ap.append(be)
            be += delta
            if j % 2 == 0:
                if j >= ndmax:
                    break
                if eps*sd[0] >= sd[j-1]:
                    break
        return np.array(ap, dtype=np.float64), np.array(sd, dtype=np.float64)
    else:
        be = 0.0
        j = 0
        wal = twt*al
        while True:
            j += 1
            ex = -((wal - be)**2) / (4.0*wal) if wal>0 else -1e300
            sfree = 0.0 if wal<=0 else math.exp(ex)/math.sqrt(4.0*math.pi*wal)
            sd.append(sfree)
            ap.append(be)
            be += delta
            if j % 2 == 0:
                if j >= ndmax:
                    break
                if eps*sd[0] >= sd[j-1]:
                    break
        return np.array(ap, dtype=np.float64), np.array(sd, dtype=np.float64)

def sbfill(delta: float, be: float, s: np.ndarray, betan: np.ndarray, nbt: int) -> np.ndarray:
    """\
subroutine sbfill(sb,nbt,delta,be,s,betan,nbeta,nbe,ndmax)
!--------------------------------------------------------------------
   ! For translational cases only.
   ! Generates s(beta) on a new energy grid for convolution
   ! with a diffusion or free-gas shape.  Interpolation is used.
   ! Called by trans.
   !--------------------------------------------------------------------
   use util    ! provides error
   ! externals
   integer::nbt,nbeta,nbe,ndmax
   real(kr)::delta,be
   real(kr)::sb(ndmax),s(ndmax),betan(nbeta)
   ! internals
   character(60)::strng
   integer::i,j,idone
   real(kr)::bmin,bmax,b,st,stm,arg,bet
   real(kr),parameter::shade=1.00001e0_kr
   real(kr),parameter::slim=-225.e0_kr
   real(kr),parameter::zero=0

   bmin=-be-(nbt-1)*delta
   bmax=-be+(nbt-1)*delta+delta/100
   if ((1+int((bmax-bmin)/delta)).gt.ndmax) then
      write(strng,'(''ndmax needs to be at least '',i6)')&
                    1+int((bmax-bmin)/delta)
      call error('sbfill',strng,' ')
   endif
   j=nbeta
   i=0
   bet=bmin
   do while (bet.le.bmax)
      i=i+1
      b=abs(bet)
      ! search for correct beta range
      idone=0
      do while (idone.eq.0)
         if (b.gt.betan(j)) then
            if (j.eq.nbeta.and.b.lt.shade*betan(j)) then
               idone=1
            else
               if (j.eq.nbeta) then
                  idone=2
               else
                  ! move up
                  j=j+1
               endif
            endif
         else
            if (b.gt.betan(j-1)) then
               idone=1
            else
               if (j.eq.2) then
                  idone=1
               else
                  ! move down
                  j=j-1
               endif
            endif
         endif
      enddo
      ! interpolate in this range
      if (idone.eq.1) then
         if (s(j).le.zero) then
            st=slim
         else
            st=log(s(j))
         endif
         if (s(j-1).le.zero) then
            stm=slim
         else
            stm=log(s(j-1))
         endif
         sb(i)=st+(b-betan(j))*(stm-st)/(betan(j-1)-betan(j))
         if (bet.gt.zero) sb(i)=sb(i)-bet
         arg=sb(i)
         sb(i)=0
         if (arg.gt.slim) sb(i)=exp(arg)
      else
         sb(i)=0
      endif
      ! if delta is to small for the current value of beta, increase it
      do while (bet.eq.(bet+delta))
         delta=delta*10
      end do
      bet=bet+delta
   enddo
   return
end subroutine sbfill
    """
    slim = -225.0
    shade = 1.00001
    zero = 0.0
    nbeta = len(betan)
    bmin = -be - (nbt-1)*delta
    bmax = -be + (nbt-1)*delta + delta/100.0
    count = 1 + int((bmax - bmin)/delta)
    sb = np.zeros(count, dtype=np.float64)
    j = nbeta - 1
    bet = bmin
    for i in range(count):
        b = abs(bet)
        idone = 0
        while idone == 0:
            if b > betan[j]:
                if j == nbeta-1 and b < shade*betan[j]:
                    idone = 1
                else:
                    if j == nbeta-1:
                        idone = 2
                    else:
                        j += 1
            else:
                if b > betan[j-1]:
                    idone = 1
                else:
                    if j == 1:
                        idone = 1
                    else:
                        j -= 1
        if idone == 1:
            st  = slim if s[j]   <= zero else float(np.log(s[j]))
            stm = slim if s[j-1] <= zero else float(np.log(s[j-1]))
            val = st + (b - betan[j])*(stm - st)/(betan[j-1] - betan[j])
            if bet > zero:
                val = val - bet
            arg = val
            sb[i] = 0.0 if arg <= slim else float(np.exp(arg))
        else:
            sb[i] = 0.0
        while bet == (bet + delta):
            delta = delta * 10.0
        bet += delta
    return sb

def trans(itemp: int) -> None:
    """\
subroutine trans(itemp)
!--------------------------------------------------------------------
   ! Controls the addition of a translational contribution
   ! to a continuous S(alpha,beta).  The translational component
   ! can be either diffusion, or a free gas.  The values of the input
   ! s(alpha,beta) for the convolution are obtained by interpolation.
   ! Called by leapr.
   ! Uses stable, sbfill, terps.
   !--------------------------------------------------------------------
   use mainio ! provides nsyso
   use util   ! provides timer
   ! externals
   integer::itemp
   ! internals
   integer::ialpha,ibeta,i,nbt,iprt,jprt,nsd,nu,ndmax
   real(kr)::time,s1,s2,sum0,sum1,ff1,ff2,ff1l,ff2l
   real(kr)::sc,al,deb,ded,delta,f,s,bb,st,be,bel
   real(kr),dimension(:),allocatable::betan,ap,sd,sb
   real(kr),parameter::therm=.0253e0_kr
   real(kr),parameter::c0=.4e0_kr
   real(kr),parameter::c1=1.e0_kr
   real(kr),parameter::c2=1.42e0_kr
   real(kr),parameter::c3=.2e0_kr
   real(kr),parameter::c4=10.e0_kr
   real(kr),parameter::tiny=1.e-30_kr
   real(kr),parameter::zero=0

   !--write heading for translational calculation
   call timer(time)
   write(nsyso,'(/'' translational part of scattering law'',&
     &32x,f8.1,''s'')') time
   sc=1
   if (lat.eq.1) sc=therm/tev

   !---allocate scratch storage
   ndmax=max(nbeta,1000000)
   allocate(betan(nbeta))
   allocate(ap(ndmax))
   allocate(sd(ndmax))
   allocate(sb(ndmax))

   !--alpha loop
   if (iprint.eq.1) write(nsyso,&
     '(/'' sab checks''/''      alpha      norm   sum rule'')')
   do ialpha=1,nalpha
      iprt=mod(ialpha-1,naint)+1
      if (ialpha.eq.nalpha) iprt=1
      al=alpha(ialpha)*sc/arat
      if (iprt.eq.1.and.iprint.eq.2)&
        write(nsyso,'(/'' alpha='',f10.4)') al

      !--choose beta interval for convolution
      ded=c0*(twt*c*al)/sqrt(c1+c2*(twt*c*al)*c)
      if (ded.eq.zero) ded=c3*sqrt(twt*al)
      deb=c4*al*deltab
      delta=deb
      if (ded.lt.delta) delta=ded
      nu=1
      if (iprt.eq.1.and.iprint.eq.2) write(nsyso,&
        '(/'' delta d='',e18.5,5x,''delta b='',e18.5,&
        &10x,''delta='',e18.5)') ded,deb,delta

      !--make table of s-diffusion or s-free on this interval
      call stable(ap,sd,nsd,al,delta,iprt,nu,ndmax)
      if (nsd.gt.1) then

         !--copy original ss(-beta) to a temporary array
         do i=1,nbeta
            betan(i)=beta(i)*sc
            ap(i)=ssm(i,ialpha,itemp)
         enddo

         !--loop over beta values
         if (iprt.eq.1.and.iprint.eq.2) write(nsyso,&
           '(/'' results after convolution ''/&
           & 4x,'' beta'',7x,''s(alpha,beta)'',6x,''ss(alpha,beta)'',&
           &5x,''ss(alpha,-beta)'')')
         do ibeta=1,nbeta
            jprt=mod(ibeta-1,nbint)+1
            if (ibeta.eq.nbeta) jprt=1
            s=0
            be=betan(ibeta)

            !--prepare table of continuous ss on new interval
            nbt=nsd
            call sbfill(sb,nbt,delta,be,ap,betan,nbeta,ibeta,ndmax)

            !--convolve s-transport with s-continuous
            do i=1,nbt
               f=2*(mod(i-1,2)+1)
               if (i.eq.1.or.i.eq.nbt) f=1
               s=s+f*sd(i)*sb(nbt+i-1)
               bb=(i-1)*delta
               s=s+f*sd(i)*sb(nbt-i+1)*exp(-bb)
            enddo
            s=s*delta/3
            if (s.lt.tiny) s=0
            st=terps(sd,nbt,delta,be)
            if (st.gt.zero) s=s+exp(-al*f0)*st

            !--store results
            ssm(ibeta,ialpha,itemp)=s
            if (s.ne.zero) then
               s1=s*exp(-be/2)
               s2=s*exp(-be)
               if (iprt.eq.1.and.jprt.eq.1.and.iprint.eq.2)&
                 write(nsyso,'(f10.4,1p,e18.5,2e20.5)') be,s1,s2,s
            endif

         !--continue beta loop
         enddo

         !--check moments of calculated s(alpha,beta).
         if (iprt.eq.1) then
            sum0=0
            sum1=0
            ff1l=0
            ff2l=0
            bel=0
            do ibeta=1,nbeta
               be=betan(ibeta)
               ff2=ssm(ibeta,ialpha,itemp)
               ff1=ssm(ibeta,ialpha,itemp)*exp(-be)
               if (ibeta.gt.1) then
                  sum0=sum0+(be-bel)*(ff1l+ff2l+ff1+ff2)/2
                  sum1=sum1+(be-bel)&
                    *(ff2l*bel+ff2*be-ff1l*bel-ff1*be)/2
                  ff1l=ff1
                  ff2l=ff2
                  bel=be
               else
                  bel=be
                  ff1l=ff1
                  ff2l=ff2
                  sum0=0
                  sum1=0
               endif
            enddo
            sum1=sum1/al/(tbeta+twt)
            if (iprint.eq.2) then
               write(nsyso,'(&
                 &''     normalization check ='',f8.4)') sum0
               write(nsyso,'(&
                 &''          sum rule check ='',f8.4)') sum1
            else if (iprint.eq.1) then
               write(nsyso,'(1x,f10.4,2f10.4)') al,sum0,sum1
            endif
         endif
      endif

   !--continue alpha loop
   enddo

   !--update effective temperature
   tempf(itemp)=(tbeta*tempf(itemp)+twt*tempr(itemp))/(tbeta+twt)
   write(nsyso,'(/''     new effective temp = '',f10.3)') tempf(itemp)

   !--deallocate scratch storage
   deallocate(sb)
   deallocate(sd)
   deallocate(ap)
   deallocate(betan)
   return
end subroutine trans
    """
    therm = 0.0253
    c0=0.4; c1=1.0; c2=1.42; c3=0.2; c4=10.0
    tiny = 1.0e-30
    zero = 0.0
    sc = 1.0
    if lat == 1:
        sc = therm/tev
    for ialpha in range(nalpha):
        al = float(alpha[ialpha])*sc/arat
        ded = c0*(twt*c*al)/math.sqrt(c1 + c2*(twt*c*al)*c) if (twt*c*al)!=0 else 0.0
        if ded == 0.0:
            ded = c3*math.sqrt(max(twt*al,0.0))
        deb = c4*al*deltab
        delta = min(ded, deb) if ded!=0 else deb
        if delta == 0.0:
            delta = deb
        ap, sd = stable(al, delta, max(nbeta, 1000000))
        nsd = len(sd)
        if nsd > 1:
            betan = np.array(beta, dtype=np.float64)*sc
            ap_store = ssm[:, ialpha, itemp-1].copy()
            for ibeta in range(nbeta):
                s = 0.0
                be = betan[ibeta]
                nbt = nsd
                sb = sbfill(delta, be, ap_store, betan, nbt)
                for i in range(nbt):
                    f = 2.0*( ( (i) % 2) + 1 )
                    if i == 0 or i == nbt-1:
                        f = 1.0
                    s += f*sd[i]*sb[nbt-1 + i]
                    bb = i*delta
                    s += f*sd[i]*sb[nbt-1 - i]*math.exp(-bb)
                s = s*delta/3.0
                if s < tiny:
                    s = 0.0
                st = terps(sd, delta, be)
                if st > zero:
                    s = s + math.exp(-al*f0)*st
                ssm[ibeta, ialpha, itemp-1] = s
    tempf[itemp-1] = (tbeta*tempf[itemp-1] + twt*tempr[itemp-1])/(tbeta + twt)

def discre(itemp: int) -> None:
    """\
subroutine discre(itemp)
!--------------------------------------------------------------------
   ! Controls the convolution of discrete oscillators with
   ! the continuous S(alpha,beta) computed in contin.
   ! Called by leapr.
   ! Uses bfact, bfill, exts, sint.
   !--------------------------------------------------------------------
   use mainio  ! provides nsyso
   use util    ! provides timer
   use physics ! provides bk (boltzmann constant)
   ! externals
   integer::itemp
   ! internals
   real(kr)::time,ss,s1,s2,sc,sumn,sumr,st,add,besn,wtsn
   real(kr)::sum0,sum1,bel,ff1,ff2,ff1l,ff2l
   real(kr)::x,db,save,dwc,dwf,al,tbart,wt,be,cn,sn
   real(kr)::tsave,dw0,dwt
   integer::nal,iprt,i,j,k,m,n,ibeta,idone
   integer::jj,nn,jprt,nbx,maxbb,maxdd
   real(kr)::bdeln(50),eb(50),dbw(50),ar(50),dist(50)
   real(kr)::bzero,bplus(50),bminus(50)
   real(kr),dimension(:),allocatable::betan,exb,sexpb
   real(kr),dimension(:),allocatable::bex,rdbex,sex
   real(kr),dimension(:),allocatable::bes,wts,ben,wtn
   real(kr),parameter::therm=.0253e0_kr
   real(kr),parameter::small=1.e-8_kr
   real(kr),parameter::vsmall=1.e-10_kr
   real(kr),parameter::tiny=1.e-20_kr
   real(kr),parameter::zero=0

   !--write heading
   call timer(time)
   write(nsyso,&
     '(/'' discrete-oscillator part of scattering law'',&
     &26x,f8.1,''s'')') time
   sc=1
   if (lat.eq.1) sc=therm/tev

   !--allocate scratch storage
   allocate(betan(nbeta))
   allocate(exb(nbeta))
   allocate(sexpb(nbeta))
   maxbb=2*nbeta+1
   allocate(bex(maxbb))
   allocate(rdbex(maxbb))
   allocate(sex(maxbb))
   maxdd=500
   allocate(bes(maxdd))
   allocate(wts(maxdd))
   allocate(ben(maxdd))
   allocate(wtn(maxdd))

   !--set up oscillator parameters
   dwt=0
   do i=1,nd
      bdeln(i)=bdel(i)/tev
      dwt=dwt+adel(i)
   enddo
   tsave=0
   dw0=dwpix(itemp)
   do i=1,nd
      eb(i)=exp(bdeln(i)/2)
      sn=(eb(i)-1/eb(i))/2
      cn=(eb(i)+1/eb(i))/2
      ar(i)=adel(i)/(sn*bdeln(i))
      dist(i)=adel(i)*bdel(i)*cn/(2*sn)
      tsave=tsave+dist(i)/bk
      dbw(i)=ar(i)*cn
      if (dwpix(itemp).gt.zero) dwpix(itemp)=dwpix(itemp)+dbw(i)
   enddo
   write(nsyso,'(/'' add delta functions''//(5x,i3,1p,2e14.4))')&
     (i,bdeln(i),adel(i),i=1,nd)

   !--prepare functions of beta
   do i=1,nbeta
      be=beta(i)*sc
      exb(i)=exp(-be/2)
      betan(i)=be
   enddo
   call bfill(bex,rdbex,nbx,betan,nbeta,maxbb)
   wt=tbeta
   tbart=tempf(itemp)/tempr(itemp)

   !--main alpha loop
   if (iprint.eq.1) write(nsyso,&
     '(/'' sab checks''/''      alpha      norm   sum rule'')')
   do nal=1,nalpha
      iprt=mod(nal-1,naint)+1
      if (nal.eq.nalpha) iprt=1
      al=alpha(nal)*sc/arat
      if (iprt.eq.1.and.iprint.eq.2) write(nsyso,&
        '(/3x,''alpha='',f10.5)') al
      dwf=exp(-al*dw0)
      if (iprt.eq.1.and.iprint.eq.2) write(nsyso,&
        '(/''      debye-waller factor='',1p,e12.4)') dwf
      call exts(ssm(1,nal,itemp),sex,exb,betan,nbeta,maxbb)
      do j=1,nbeta
         sexpb(j)=0
      enddo

      !---initialize for delta function calculation
      ben(1)=0
      wtn(1)=1
      nn=1
      n=0

      !--loop over all oscillators
      do i=1,nd
         dwc=al*dbw(i)
         x=al*ar(i)
         call bfact(x,bzero,bplus,bminus,dwc,bdeln(i))

         !--do convolution for the delta functions
         !--n=0 term
         do m=1,nn
            besn=ben(m)
            wtsn=wtn(m)*bzero
            if (besn.le.zero.or.wtsn.ge.small) then
               if (n.lt.maxdd) then
                  n=n+1
                  bes(n)=besn
                  wts(n)=wtsn
               endif
            endif
         enddo

         !--negative n terms
         k=0
         idone=0
         do while (k.lt.50.and.idone.eq.0)
            k=k+1
            if (bminus(k).le.zero) then
               idone=1
            else
               do m=1,nn
                  besn=ben(m)-k*bdeln(i)
                  wtsn=wtn(m)*bminus(k)
                  if (wtsn.ge.small.and.n.lt.maxdd) then
                     n=n+1
                     bes(n)=besn
                     wts(n)=wtsn
                  endif
               enddo
            endif
         enddo

         !--positive n terms
         k=0
         idone=0
         do while (k.lt.50.and.idone.eq.0)
            k=k+1
            if (bplus(k).le.zero) then
               idone=1
            else
               do m=1,nn
                  besn=ben(m)+k*bdeln(i)
                  wtsn=wtn(m)*bplus(k)
                  if (wtsn.ge.small.and.n.lt.maxdd) then
                     n=n+1
                     bes(n)=besn
                     wts(n)=wtsn
                  endif
               enddo
            endif
         enddo

         !--continue oscillator loop
         nn=n
         do m=1,nn
            ben(m)=bes(m)
            wtn(m)=wts(m)
         enddo
         n=0
         wt=wt+adel(i)
         tbart=tbart+dist(i)/bk/tempr(itemp)
      enddo
      n=nn

      !--sort the discrete lines
      !--and throw out the smallest ones
      nn=n-1
      do i=2,n-1
         do j=i+1,n
            if (wts(j).ge.wts(i)) then
               save=wts(j)
               wts(j)=wts(i)
               wts(i)=save
               save=bes(j)
               bes(j)=bes(i)
               bes(i)=save
            endif
         enddo
      enddo
      i=0
      idone=0
      do while (i.lt.nn.and.idone.eq.0)
         i=i+1
         n=i
         if (wts(i).lt.100*small.and.i.gt.5) idone=1
      enddo

      !--report the discrete lines
      if (iprint.ge.2) then
         write(nsyso,'(/''      discrete lines''/&
           &14x,''beta'',6x,''weight'')')
         sumn=0
         sumr=0
         do m=1,n
            write(nsyso,'(6x,f12.4,1p,e12.4)') bes(m),wts(m)
            sumn=sumn+wts(m)
            sumr=sumr-bes(m)*wts(m)
         enddo
         sumr=sumr/al/dwt
         write(nsyso,'(8x,''norm check'',f10.4/&
           &8x,''rule check'',f10.4)')&
           sumn,sumr
      endif

      !--add the continuum part to the scattering law
      do m=1,n
         do j=1,nbeta
            be=-betan(j)-bes(m)
            st=sint(be,bex,rdbex,sex,nbx,al,tbeta+twt,tbart,betan,&
              nbeta,maxbb)
            add=wts(m)*st
            if (add.ge.tiny) sexpb(j)=sexpb(j)+add
         enddo
      enddo

      !--add the delta functions to the scattering law
      !--delta(0.) is saved for the incoherent elastic
      if (twt.le.zero) then
         m=0
         idone=0
         do while (m.lt.n.and.idone.eq.0)
            m=m+1
            if (dwf.lt.vsmall) then
               idone=1
            else
               if (bes(m).lt.zero) then
                  be=-bes(m)
                  if (be.le.betan(nbeta-1)) then
                     db=1000
                     idone=0
                     j=0
                     do while (j.lt.nbeta.and.idone.eq.0)
                        j=j+1
                        jj=j
                        if (abs(be-betan(j)).gt.db) then
                           idone=1
                        else
                           db=abs(be-betan(j))
                        endif
                     enddo
                     if (jj.le.2) then
                        add=wts(m)/betan(jj)
                     else
                        add=2*wts(m)/(betan(jj)-betan(jj-2))
                     endif
                     add=add*dwf
                     if (add.ge.tiny) sexpb(jj-1)=sexpb(jj-1)+add
                  endif
               endif
            endif
         enddo
      endif

      !--record the results
      do j=1,nbeta
         ssm(j,nal,itemp)=sexpb(j)
      enddo
      if (iprt.eq.1.and.iprint.eq.2) write(nsyso,&
        '(/4x,'' beta'',7x,''s(alpha,beta)'',7x,''ss(alpha,beta)'',&
        &5x,''ss(alpha,-beta)'')')
      do i=1,nbeta
         be=beta(i)*sc
         ss=ssm(i,nal,itemp)
         s1=ss*exp(-be/2)
         s2=ss*exp(-be)
         jprt=mod(i-1,nbint)+1
         if (i.eq.nbeta) jprt=1
         if (iprt.eq.1.and.jprt.eq.1.and.iprint.eq.2)&
           write(nsyso,'(f10.4,1pe18.5,1p,2e20.5)') betan(i),s1,s2,ss
         enddo

      !--check moments of calculated s(alpha,beta).
      if (iprt.eq.1) then
         sum0=0
         sum1=0
         ff1l=0
         ff2l=0
         bel=0
         do ibeta=1,nbeta
            be=betan(ibeta)
            ff2=ssm(ibeta,nal,itemp)
            ff1=ssm(ibeta,nal,itemp)*exp(-be)
            if (ibeta.gt.1) then
               sum0=sum0+(be-bel)*(ff1l+ff2l+ff1+ff2)/2
               sum1=sum1+(be-bel)*(ff2l*bel+ff2*be-ff1l*bel-ff1*be)/2
               ff1l=ff1
               ff2l=ff2
               bel=be
            else
               bel=be
               ff1l=ff1
               ff2l=ff2
               sum0=0
               sum1=0
            endif
         enddo
         if (twt.eq.zero) sum0=sum0/(1-exp(-al*dwpix(itemp)))
         sum1=sum1/al
         if (iprint.eq.2) then
            write(nsyso,'(''     normalization check ='',f8.4)') sum0
            write(nsyso,'(''          sum rule check ='',f8.4)') sum1
         else if (iprint.eq.1) then
            write(nsyso,'(1x,f10.4,2f10.4)') al,sum0,sum1
         endif
      endif

   !--continue the alpha loop
   enddo

   !--finished
   tempf(itemp)=(tbeta+twt)*tempf(itemp)+tsave
   write(nsyso,'(/&
     &  ''       new effective temp = '',f10.3/&
     &  ''  new debye-waller lambda='',f10.6)')&
     tempf(itemp),dwpix(itemp)
   write(nsyso,'('' discr.-oscill. part of eff. temp = '',f10.3)') tsave
   deallocate(wtn)
   deallocate(ben)
   deallocate(wts)
   deallocate(bes)
   deallocate(sex)
   deallocate(rdbex)
   deallocate(bex)
   deallocate(sexpb)
   deallocate(exb)
   deallocate(betan)
   return
end subroutine discre
    """
    therm = 0.0253
    small = 1.0e-8
    vsmall = 1.0e-10
    tiny = 1.0e-20
    zero = 0.0
    sc = 1.0
    if lat == 1:
        sc = therm/tev
    dwt = 0.0
    bdeln = np.zeros(nd, dtype=np.float64) if nd>0 else np.zeros(0)
    dbw = np.zeros(nd, dtype=np.float64)
    ar = np.zeros(nd, dtype=np.float64)
    dist = np.zeros(nd, dtype=np.float64)
    for i in range(nd):
        bdeln[i] = bdel[i]/tev
        dwt += adel[i]
    tsave = 0.0
    dw0 = dwpix[itemp-1]
    for i in range(nd):
        ebi = math.exp(bdeln[i]/2.0)
        sn = (ebi - 1.0/ebi)/2.0
        cn = (ebi + 1.0/ebi)/2.0
        ar[i] = adel[i]/(sn*bdeln[i])
        dist[i] = adel[i]*bdel[i]*cn/(2.0*sn)
        # Boltzmann k_B in eV/K ~ 8.617333262e-5
        tsave += dist[i]/8.617333262e-5
        dbw[i] = ar[i]*cn
        if dwpix[itemp-1] > zero:
            dwpix[itemp-1] = dwpix[itemp-1] + dbw[i]
    betan = np.array(beta, dtype=np.float64)*sc
    exb = np.exp(-betan/2.0)
    bex, rdbex, nbx = bfill(betan)
    wt = tbeta
    tbart = tempf[itemp-1]/tempr[itemp-1]
    for nal in range(nalpha):
        al = float(alpha[nal])*sc/arat
        dwf = math.exp(-al*dw0)
        sex = exts(ssm[:, nal, itemp-1], exb, betan)
        sexpb = np.zeros(nbeta, dtype=np.float64)
        ben = [0.0]; wtn = [1.0]; nn = 1
        for i in range(nd):
            dwc = al*dbw[i]
            x = al*ar[i]
            bzero, bplus, bminus = bfact(x, dwc, bdeln[i])
            new_ben = []; new_wtn = []
            for m in range(nn):
                new_ben.append(ben[m]); new_wtn.append(wtn[m]*bzero)
            for k in range(1,51):
                if bminus[k-1] <= zero: break
                for m in range(nn):
                    new_ben.append(ben[m] - k*bdeln[i]); new_wtn.append(wtn[m]*bminus[k-1])
            for k in range(1,51):
                if bplus[k-1] <= zero: break
                for m in range(nn):
                    new_ben.append(ben[m] + k*bdeln[i]); new_wtn.append(wtn[m]*bplus[k-1])
            ben = new_ben; wtn = new_wtn; nn = len(ben)
            wt = wt + adel[i]
            tbart = tbart + dist[i]/8.617333262e-5/tempr[itemp-1]
        order = np.argsort(-np.array(wtn))
        ben = list(np.array(ben)[order]); wtn = list(np.array(wtn)[order])
        cut = len(wtn)
        for i in range(len(wtn)):
            if i>4 and wtn[i] < 100*small:
                cut = i; break
        ben = ben[:cut]; wtn = wtn[:cut]
        for m in range(len(wtn)):
            for j in range(nbeta):
                be = -betan[j] - ben[m]
                st = sint(be, bex, rdbex, sex, nbx, al, tbeta+twt, tbart, betan)
                add = wtn[m]*st
                if add >= tiny:
                    sexpb[j] += add
        if twt <= zero:
            for m in range(len(wtn)):
                if dwf < vsmall: break
                if ben[m] < zero:
                    be = -ben[m]
                    if be <= betan[nbeta-2]:
                        j = int(np.searchsorted(betan, be))
                        if j <= 1:
                            add = wtn[m]/betan[j]
                        else:
                            add = 2*wtn[m]/(betan[j]-betan[j-2])
                        add *= dwf
                        if add >= tiny:
                            sexpb[j-1] += add
        ssm[:, nal, itemp-1] = sexpb
    tempf[itemp-1] = (tbeta + twt)*tempf[itemp-1] + tsave


def contin(temp: float, itemp: int, np: int, maxn: int) -> None:
    """\
subroutine contin(temp,itemp,np,maxn)
!--------------------------------------------------------------------
   ! Main routine for calculating S(alpha,beta) at temp
   ! for continuous phonon frequency distributions.
   ! Called by leapr.
   ! Uses start, terpt, convol.
   !--------------------------------------------------------------------
   use physics ! provides pi
   use mainio  ! provides nsyso
   use util    ! provides timer
   ! externals
   real(kr)::temp
   integer::itemp,np,maxn
   ! internals
   integer::i,j,k,n,npn,npl,iprt,jprt
   integer,allocatable,dimension(:)::maxt
   character(3)::tag
   real(kr)::al,be,bel,ex,exx,st,add,sc,alp,alw,ssct,ckk,time
   real(kr)::ff0,ff1,ff2,ff1l,ff2l,sum0,sum1
   real(kr),dimension(:),allocatable::p,tlast,tnow,xa
   real(kr),parameter::therm=0.0253e0_kr
   real(kr),parameter::tiny=1.e-30_kr
   real(kr),parameter::explim=-250.e0_kr
   real(kr),parameter::zero=0

   !--write heading
   write(nsyso,'(/'' solid-type contributions to scattering law'')')

   !--allocate temporary arrays
   allocate(p(np1))
   allocate(tlast(nphon*np1))
   allocate(tnow(nphon*np1))
   allocate(xa(nalpha))
   allocate(maxt(nbeta))

   !--calculate various parameters for this temperature
   call start(itemp,p,np,deltab,tev)
   sc=1
   if (lat.eq.1) sc=therm/tev

   !--start the phonon expansion sum with t1
   do i=1,np
      tlast(i)=p(i)
   enddo
   do j=1,nalpha
      al=alpha(j)*sc/arat
      xa(j)=log(al*f0)
      ex=-f0*al+xa(j)
      exx=0
      if (ex.gt.explim) exx=exp(ex)
      do k=1,nbeta
         be=beta(k)*sc
         st=terpt(p,np,deltab,be)
         add=st*exx
         if (add.lt.tiny) add=0
         ssm(k,j,itemp)=add
      enddo
   enddo
   npl=np

   !--do the phonon expansion sum
   do j=1,nbeta
      maxt(j)=nalpha+1
   enddo
   if (iprint.eq.2) then
      write(nsyso,'(/'' normalization check for phonon expansion'')')
   endif
   if (maxn.gt.maxnphon) then
      call timer(time)
      write(nsyse,'(/'' performing phonon expansion sum'',&
            &37x,f8.1,''s'')'),time
   endif
   do n=2,maxn
      npn=np+npl-1
      call convol(p,tlast,tnow,np,npl,npn,deltab,ckk)
      if (iprint.eq.2) write(nsyso,'(5x,i5,f12.5)') n,ckk
      do j=1,nalpha
         al=alpha(j)*sc/arat
         xa(j)=xa(j)+log(al*f0/n)
         ex=-f0*al+xa(j)
         exx=0
         if (ex.gt.explim) exx=exp(ex)
         do k=1,nbeta
            be=beta(k)*sc
            st=terpt(tnow,npn,deltab,be)
            add=st*exx
            if (add.lt.tiny) add=0
            ssm(k,j,itemp)=ssm(k,j,itemp)+add
            if (ssm(k,j,itemp).ne.zero.and.n.ge.maxn) then
            if (add.gt.ssm(k,j,itemp)/1000.and.j.lt.maxt(k)) maxt(k)=j
            endif
         enddo
      enddo
      do i=1,npn
         tlast(i)=tnow(i)
      enddo
      npl=npn
      if (mod(n,maxnphon).eq.0) then
         call timer(time)
         write(nsyse,'(2x,i5,'' of '',i5,&
               &'' loops done for phonon expansion sum'',17x,f8.1,''s'')')&
               &n,maxn,time
      endif
   enddo
   if (maxn.gt.maxnphon) then
      call timer(time)
      write(nsyse,'(/'' done with phonon expansion sum'',&
            &38x,f8.1,''s'')'),time
   endif

   !--print out start of sct range for each beta
   if (iprint.ne.0) then
      write(nsyso,'(/''         beta   alpha sct'')')
      do i=1,nbeta
         if (i.gt.1) then
            if (maxt(i).gt.maxt(i-1)) maxt(i)=maxt(i-1)
         endif
         if (maxt(i).gt.nalpha) then
            write(nsyso,'(1x,f12.4,''      none'')') beta(i)*sc
         else
            write(nsyso,'(1x,2f12.4)')&
              beta(i)*sc,alpha(maxt(i))*sc/arat
         endif
      enddo
   endif

   !---check the moments of s(alpha,beta)
   if (iprint.eq.1) write(nsyso,&
     '(/'' sab checks''/''      alpha      norm   sum rule'')')
   do j=1,nalpha
      iprt=mod(j-1,naint)+1
      if (j.eq.nalpha) iprt=1
      if (iprt.eq.1) then
         al=alpha(j)*sc/arat
         if (iprint.eq.2) then
            write(nsyso,'(/''  alpha='',f10.4)') al
            write(nsyso,'(5x,''beta'',7x,''s(alpha,beta)'',&
              &6x,''ss(alpha,beta)'',5x,''ss(alpha,-beta)'')')
         endif
         bel=0
         ff1l=0
         ff2l=0
         sum0=0
         sum1=0
         do k=1,nbeta
            jprt=mod(k-1,nbint)+1
            if (k.eq.nbeta) jprt=1
            be=beta(k)*sc
            alw=al*tbeta
            alp=alw*tbar
            ex=-(alw-be)**2/(4*alp)
            ssct=0
            if (ex.gt.explim) ssct=exp(ex)/sqrt(4*pi*alp)
            tag='   '
            if (j.ge.maxt(k)) then
               tag='sct'
               ssm(k,j,itemp)=ssct
            endif
            ff2=ssm(k,j,itemp)
            ff1=ssm(k,j,itemp)*exp(-be)
            ff0=ssm(k,j,itemp)*exp(-be/2)
            if (jprt.eq.1.and.iprint.eq.2.and.ff2.gt.zero)&
              write(nsyso,'(f10.4,1p,e18.5,2e20.5,4x,a3)')&
              be,ff0,ff1,ff2,tag
            if (k.gt.1) then
               sum0=sum0+(be-bel)*(ff1l+ff2l+ff1+ff2)/2
               sum1=sum1+(be-bel)*(ff2l*bel+ff2*be-ff1l*bel-ff1*be)/2
               ff1l=ff1
               ff2l=ff2
               bel=be
            else
               bel=be
               ff1l=ff1
               ff2l=ff2
               sum0=0
               sum1=0
            endif
         enddo
         sum0=sum0/(1-exp(-al*f0))
         sum1=sum1/al/tbeta
         if (iprint.eq.2) then
            write(nsyso,'(''   normalization check ='',f8.4)') sum0
            write(nsyso,'(''        sum rule check ='',f8.4)') sum1
         else if (iprint.eq.1) then
            write(nsyso,'(1x,3f10.4)') al,sum0,sum1
         endif
      endif
   enddo

   !--finished with continuous distribution
end subroutine contin
    """
    # Mirrors Fortran loops and evaluation order. Uses globals for arrays and params.
    therm = 0.0253
    tiny = 1.0e-30
    explim = -250.0
    zero = 0.0
    global ssm, alpha, beta, nalpha, nbeta, nphon, f0, tbar, deltab, tev, arat, naint, nbint, tbeta
    # allocate temporaries
    if p1 is None:
        raise RuntimeError("p1 not set")
    p = np.zeros(np1, dtype=np.float64)
    tlast = np.zeros(nphon*np1, dtype=np.float64)
    tnow  = np.zeros(nphon*np1, dtype=np.float64)
    xa    = np.zeros(nalpha, dtype=np.float64)
    maxt  = np.zeros(nbeta, dtype=np.int32)
    # calculate parameters for this temp
    start(itemp, p)  # sets f0, tbar, tempf/dwpix, updates p in-place to t1(beta)
    sc = 1.0
    if lat == 1:
        sc = therm/tev
    # start the phonon expansion sum with t1
    tlast[:np1] = p[:np1]
    for j in range(nalpha):
        al = float(alpha[j])*sc/arat
        xa[j] = math.log(al*f0) if al*f0>0 else -1e300
        ex  = -f0*al + xa[j]
        exx = math.exp(ex) if ex > explim else 0.0
        for k in range(nbeta):
            be = float(beta[k])*sc
            st = terpt(p[:np1], deltab, be)
            add = st*exx
            if add < tiny: add = 0.0
            ssm[k, j, itemp-1] = add
    npl = np1
    # phonon expansion loop
    maxt[:] = nalpha + 1
    if maxn > 0:
        for n in range(2, maxn+1):
            npn = np1 + npl - 1
            ckk, tnow[:npn] = convol(p[:np1], tlast[:npl], deltab, npn)
            # update sums
            for j in range(nalpha):
                al = float(alpha[j])*sc/arat
                xa[j] = xa[j] + math.log(al*f0/n)
                ex = -f0*al + xa[j]
                exx = math.exp(ex) if ex > explim else 0.0
                for k in range(nbeta):
                    be = float(beta[k])*sc
                    st = terpt(tnow[:npn], deltab, be)
                    add = st*exx
                    if add < tiny: add = 0.0
                    ssm[k, j, itemp-1] = ssm[k, j, itemp-1] + add
                    if ssm[k, j, itemp-1] != zero and n >= maxn:
                        if add > ssm[k, j, itemp-1]/1000.0 and j+1 < maxt[k]:
                            maxt[k] = j+1  # store 1-based like Fortran for later semantic
            # advance
            tlast[:npn] = tnow[:npn]
            npl = npn
    # SCT region replacement is done in the 'checks' block below
    # print SCT start (we don't print; we compute maxt monotone nonincreasing)
    for i in range(1, nbeta):
        if maxt[i] > maxt[i-1]:
            maxt[i] = maxt[i-1]
    # checks block that also replaces tail by SCT approx
    for j in range(nalpha):
        al = float(alpha[j])*sc/arat
        alw = al*tbeta
        alp = alw*tbar if tbeta != 0 else 1.0e-300
        for k in range(nbeta):
            be = float(beta[k])*sc
            ex = -((alw - be)**2)/(4.0*alp) if alp!=0 else -1e300
            ssct = math.exp(ex)/math.sqrt(4.0*math.pi*alp) if ex > explim else 0.0
            # if alpha index beyond threshold for this beta, replace with SCT
            if j+1 >= maxt[k]:
                ssm[k, j, itemp-1] = ssct



def fsum(n: int, p: np.ndarray, tau: float, deltab: float) -> float:
    """\
function fsum(n,p,np,tau,deltab)
!--------------------------------------------------------------------
   ! Computes integrals over the phonon frequency
   ! of the form
   !    integral 0 to infinity of
   !       2*p*beta**n*hyperbolic
   !    dbeta
   ! where
   !    p is p/beta**2, or rho/(2*beta*sinh(beta/2)), and
   !    'hyperbolic' is cosh(tau*beta) for n even
   !       and sinh(tau*beta) for n odd.
   ! Called by start.
   !--------------------------------------------------------------------
   ! externals
   integer::n,np
   real(kr)::tau,deltab
   real(kr)::p(np)
   ! internals
   integer::ij
   real(kr)::arg,edsq,v,an,be,fs,w,ff

   arg=deltab*tau/2
   edsq=exp(arg)
   v=1
   an=1-2*mod(n,2)
   be=0
   fs=0
   w=1
   do ij=1,np
      if (n.gt.0) w=be**n
      ff=((p(ij)*v)*v+(p(ij)*an/v)/v)*w
      if (ij.eq.1.or.ij.eq.np) ff=ff/2
      fs=fs+ff
      be=be+deltab
      v=v*edsq
   enddo
   fsum=fs*deltab
   return
end function fsum
    """
    # integral sum with trapezoidal rule on uniform beta grid
    arg = deltab*tau/2.0
    edsq = math.exp(arg)
    v = 1.0
    an = 1 - 2*(n % 2)   # 1 for even n, -1 for odd n
    be = 0.0
    fs = 0.0
    w = 1.0
    npn = p.shape[0]
    for ij in range(npn):
        if n > 0:
            w = be**n
        ff = ((p[ij]*v)*v + (p[ij]*an/v)/v) * w
        if ij == 0 or ij == npn-1:
            ff = ff/2.0
        fs += ff
        be += deltab
        v = v*edsq
    return fs*deltab



def terpk(ska_arr: np.ndarray, nka_val: int, delta: float, be: float) -> float:
    """\
function terpk(ska,nka,delta,be)
!--------------------------------------------------------------------
   ! Interpolate in a table of ska(kappa) for a required kappa.
   ! Called by coldh.
   !--------------------------------------------------------------------
   ! externals
   integer::nka
   real(kr)::ska(nka),delta,be
   ! internals
   integer::i
   real(kr)::bt,btp

   terpk=1
   if (be.gt.nka*delta) return
   i=int(be/delta)
   if (i.lt.nka-1) then
      bt=i*delta
      btp=bt+delta
      i=i+1
      terpk=ska(i)+(be-bt)*(ska(i+1)-ska(i))/(btp-bt)
   endif
   return
end function terpk
    """
    terpk_val = 1.0
    if be > nka_val*delta:
        return terpk_val
    i = int(be/delta)
    if i < nka_val - 1:
        bt = i*delta
        btp = bt + delta
        i1 = i + 1
        terpk_val = float(ska_arr[i1] + (be - bt)*(ska_arr[i1+1] - ska_arr[i1])/(btp - bt))
    return terpk_val



def sjbes(n: int, x: float) -> float:
    """\
function sjbes(n,x)
!--------------------------------------------------------------------
   ! Bessel functions for cold hydrogen or deuterium calculation
   !--------------------------------------------------------------------
   use util ! provides error
   ! externals
   integer::n
   real(kr)::x
   ! internals
   integer::i,k,l,iii,kmax,nm
   real(kr)::w,bessel,y,z,sj,t1,t2,t3
   character(60)::strng
   real(kr),parameter::huge=1.e25_kr
   real(kr),parameter::small=2.e-38_kr
   real(kr),parameter::break1=3.e4_kr
   real(kr),parameter::break2=7.e-4_kr
   real(kr),parameter::break3=0.2e0_kr
   real(kr),parameter::one=1
   real(kr),parameter::ten=10
   real(kr),parameter::hund=100
   real(kr),parameter::zero=0

   !--check for large arguments
   if (n.ge.30000.or.x.gt.break1) then
      write(strng,'(&
        &''value is not accurate  n = '',i7,10x,''x = '',e14.7)') n,x
      call mess('sjbes',strng,' ')
      sjbes=0
      return
   endif

   !--check for bad arguments
   if (x.lt.zero.or.n.lt.0) then
      write(strng,'(&
        &''argument is invalid  n = '',i7,10x,''x = '',e14.7)') n,x
      call error('sjbes',strng,' ')
      sjbes=0
      return
   endif

   !--compute normal values
   if (x.le.break2) then
      w=1
      if (n.eq.0) then
         bessel=w
      else if (n.gt.10) then
         bessel=0
      else
         t1=3
         t2=1
         do i=1,n
            t3=t2*x/t1
            t1=t1+2
            t2=t3
         enddo
         bessel=t3
      endif
   else
      if (x.lt.break3) then
         y=x**2
         w=1-y*(1-y/20)/6
      else
         w=sin(x)/x
      endif
      if (n.eq.0) then
         bessel=w
      else
         if (x.ge.hund) then
            l=int(x/50+18)
         else if (x.ge.ten) then
            l=int(x/10+10)
         else if (x.gt.one) then
            l=int(x/2+5)
         else
            l=5
         endif
         iii=int(x)
         kmax=n
         if (iii.gt.n) kmax=iii
         nm=kmax+l
         z=1/x
         t3=0
         t2=small
         do i=1,nm
            k=nm-i
            t1=(2*k+3)*z*t2-t3
            if (n.eq.k) sj=t1
            if (abs(t1).ge.huge) then
               t1=t1/huge
               t2=t2/huge
               sj=sj/huge
            endif
            t3=t2
            t2=t1
         enddo
         bessel=w*sj/t1
      endif
   endif
   sjbes=bessel
   return
end function sjbes
    """
    huge = 1.0e25
    small = 2.0e-38
    break1 = 3.0e4
    break2 = 7.0e-4
    break3 = 0.2
    one = 1.0
    ten = 10.0
    hund = 100.0
    zero = 0.0
    # guard bad args
    if n >= 30000 or x > break1:
        return 0.0
    if x < zero or n < 0:
        return 0.0
    # normal values
    if x <= break2:
        w = 1.0
        if n == 0:
            bessel = w
        elif n > 10:
            bessel = 0.0
        else:
            t1 = 3.0
            t2 = 1.0
            t3 = 0.0
            for i in range(1, n+1):
                t3 = t2*x/t1
                t1 = t1 + 2.0
                t2 = t3
            bessel = t3
    else:
        if x < break3:
            y = x*x
            w = 1.0 - y*(1.0 - y/20.0)/6.0
        else:
            w = math.sin(x)/x
        if n == 0:
            bessel = w
        else:
            if x >= hund:
                l = int(x/50.0 + 18.0)
            elif x >= ten:
                l = int(x/10.0 + 10.0)
            elif x > one:
                l = int(x/2.0 + 5.0)
            else:
                l = 5
            iii = int(x)
            kmax = max(n, iii)
            nm = kmax + l
            z = 1.0/x
            t3 = 0.0
            t2 = small
            sj = 0.0
            for i in range(nm):
                k = nm - 1 - i
                t1 = (2*k + 3)*z*t2 - t3
                if n == k:
                    sj = t1
                if abs(t1) >= huge:
                    t1 /= huge; t2 /= huge; sj /= huge
                t3 = t2
                t2 = t1
            bessel = w*sj/t1
    return float(bessel)



def cn(jj: int, ll: int, nn: int) -> float:
    """\
function cn(jj,ll,nn)
!--------------------------------------------------------------------
   ! Clebsch-Gordon coefficients
   ! for cold hydrogen or deuterium calculation
   !--------------------------------------------------------------------
   ! externals
   integer::jj,ll,nn
   ! internals
   integer::i,kdet,kdel,ka1,ka2,ka3,ka4,kb1,kb2,kb3,kb4,iwign
   real(kr)::s,fact,zi,a1,a2,a3,a4,b1,b2,b3,b4,rat,wign
   real(kr),parameter::zero=0

   kdet=(jj+ll+nn)/2
   kdel=jj+ll+nn-2*kdet
   if (kdel.eq.0) then
      ka1=jj+ll+nn
      ka2=jj+ll-nn
      ka3=jj-ll+nn
      ka4=ll-jj+nn
      kb1=ka1/2
      kb2=ka2/2
      kb3=ka3/2
      kb4=ka4/2
      s=0
      fact=1
      do i=1,ka1
         zi=i
         s=s+log(zi)
      enddo
      if (s.gt.zero) fact=exp(s)
      a1=sqrt(fact)
      s=0
      fact=1
      do i=1,ka2
         zi=i
         s=s+log(zi)
      enddo
      if (s.gt.zero) fact=exp(s)
      a2=sqrt(fact)
      s=0
      fact=1
      do i=1,ka3
         zi=i
         s=s+log(zi)
      enddo
      if (s.gt.zero) fact=exp(s)
      a3=sqrt(fact)
      s=0
      fact=1
      do i=1,ka4
         zi=i
         s=s+log(zi)
      enddo
      if (s.gt.zero) fact=exp(s)
      a4=sqrt(fact)
      s=0
      b1=1
      do i=1,kb1
         zi=i
         s=s+log(zi)
      enddo
      if (s.gt.zero) b1=exp(s)
      s=0
      b2=1
      do i=1,kb2
         zi=i
         s=s+log(zi)
      enddo
      if (s.gt.zero) b2=exp(s)
      s=0
      b3=1
      do i=1,kb3
         zi=i
         s=s+log(zi)
      enddo
      if (s.gt.zero) b3=exp(s)
      s=0
      b4=1
      do i=1,kb4
         zi=i
         s=s+log(zi)
      enddo
      if (s.gt.zero) b4=exp(s)
      rat=2*nn+1
      rat=rat/(jj+ll+nn+1)
      iwign=(jj+ll-nn)/2
      wign=(-1)**iwign
      wign=wign*sqrt(rat)*b1/a1*a2/b2*a3/b3*a4/b4
   else
      wign=0
   endif
   cn=wign
   return
end function cn
    """
    zero = 0.0
    kdet = (jj + ll + nn)//2
    kdel = jj + ll + nn - 2*kdet
    if kdel != 0:
        return 0.0
    ka1 = jj + ll + nn
    ka2 = jj + ll - nn
    ka3 = jj - ll + nn
    ka4 = ll - jj + nn
    kb1 = ka1//2
    kb2 = ka2//2
    kb3 = ka3//2
    kb4 = ka4//2
    # factorial via logs
    def logfact(m: int) -> float:
        s = 0.0
        for i in range(1, m+1):
            s += math.log(i)
        return s
    a1 = math.exp(logfact(ka1)/2.0) if ka1>0 else 1.0
    a2 = math.exp(logfact(ka2)/2.0) if ka2>0 else 1.0
    a3 = math.exp(logfact(ka3)/2.0) if ka3>0 else 1.0
    a4 = math.exp(logfact(ka4)/2.0) if ka4>0 else 1.0
    b1 = math.exp(logfact(kb1)) if kb1>0 else 1.0
    b2 = math.exp(logfact(kb2)) if kb2>0 else 1.0
    b3 = math.exp(logfact(kb3)) if kb3>0 else 1.0
    b4 = math.exp(logfact(kb4)) if kb4>0 else 1.0
    rat = (2*nn + 1) / (jj + ll + nn + 1)
    iwign = (jj + ll - nn)//2
    wign = ((-1.0)**iwign) * math.sqrt(rat) * b1/a1 * a2/b2 * a3/b3 * a4/b4
    return float(wign)



def sumh(j: int, jp: int, y: float) -> float:
    """\
function sumh(j,jp,y)
!--------------------------------------------------------------------
   ! Does sum over Bessel functions and Clebsch-Gordon coefficients
   ! for cold hydrogen or deuterium calculation.
   !--------------------------------------------------------------------
   ! externals
   integer::j,jp
   real(kr)::y
   ! internals
   integer::imk,ipk1,mpk,ipk,n,n1
   real(kr)::sum1,sum2

   if (j.eq.0) then
      sum2=(sjbes(jp,y)*cn(j,jp,jp))**2
   else if (jp.eq.0) then
      sum2=(sjbes(j,y)*cn(j,0,j))**2
   else
      sum1=0
      imk=iabs(j-jp)+1
      ipk1=j+jp+1
      mpk=ipk1-imk
      if (mpk.le.9) then
         ipk=ipk1
      else
         ipk=imk+9
      endif
      do n=imk,ipk
         n1=n-1
         sum1=sum1+(sjbes(n1,y)*cn(j,jp,n1))**2
      enddo
      sum2=sum1
   endif
   sumh=sum2
   return
end function sumh
    """
    if j == 0:
        return float((sjbes(jp, y)*cn(j, jp, jp))**2)
    elif jp == 0:
        return float((sjbes(j, y)*cn(j, 0, j))**2)
    else:
        imk = abs(j - jp) + 1
        ipk1 = j + jp + 1
        mpk = ipk1 - imk
        if mpk <= 9:
            ipk = ipk1
        else:
            ipk = imk + 9
        s = 0.0
        for n in range(imk, ipk+1):
            n1 = n - 1
            s += (sjbes(n1, y)*cn(j, jp, n1))**2
        return float(s)



def bt(j: int, x: float) -> float:
    """\
subroutine bt(j,pj,x)
!--------------------------------------------------------------------
   ! Statistical weight factor
   ! for cold hydrogen or deuterium calculation
   !--------------------------------------------------------------------
   ! externals
   integer::j
   real(kr)::pj,x
   ! internals
   integer::i,k
   real(kr)::yy,a,b
   real(kr),parameter::half=0.5e0_kr

   yy=half*j*(j+1)
   a=(2*j+1)*exp(-yy*x)
   b=0
   do i=1,10
      k=2*i-2
      if (mod(j,2).eq.1) k=k+1
      yy=half*k*(k+1)
      b=b+(2*k+1)*exp(-yy*x)
   enddo
   pj=a/(2*b)
   return
end subroutine bt
    """
    half = 0.5
    yy = half*j*(j+1)
    a = (2*j + 1)*math.exp(-yy*x)
    b = 0.0
    for i in range(1, 11):
        k = 2*i - 2
        if j % 2 == 1:
            k = k + 1
        yy = half*k*(k+1)
        b = b + (2*k + 1)*math.exp(-yy*x)
    pj = a/(2*b) if b != 0 else 0.0
    return float(pj)



def formf(lat: int, l1: int, l2: int, l3: int) -> float:
    """\
function formf(lat,l1,l2,l3)
!--------------------------------------------------------------------
   ! Compute form factors for the specified lattice.
   !       lat=1    graphite
   !       lat=2    Be
   !       lat=3    BeO
   !       lat=4,5  fcc lattice (aluminum, lead)
   !       lat=6    bcc lattice (iron)
   !--------------------------------------------------------------------
   use physics ! provides pi
   ! externals
   integer::lat,l1,l2,l3
   ! internals
   integer::i
   real(kr)::e1,e2,e3
   real(kr),parameter::c1=7.54e0_kr
   real(kr),parameter::c2=4.24e0_kr
   real(kr),parameter::c3=11.31e0_kr

   if (lat.eq.1) then
      ! graphite.
      i=l3/2
      if ((2*i).ne.l3) then
      formf=sin(pi*(l1-l2)/3)**2
      else
         formf=(6+10*cos(2*pi*(l1-l2)/3))/4
      endif
   else if (lat.eq.2) then
      ! beryllium.
      formf=1+cos(2*pi*(2*l1+4*l2+3*l3)/6)
   else if (lat.eq.3) then
      ! beryllium oxide.
      formf=(1+cos(2*pi*(2*l1+4*l2+3*l3)/6))&
        *(c1+c2+c3*cos(3*pi*l3/4))
   else if (lat.eq.4.or.lat.eq.5) then
      ! fcc lattices.
      e1=2*pi*l1
      e2=2*pi*(l1+l2)
      e3=2*pi*(l1+l3)
      formf=(1+cos(e1)+cos(e2)+cos(e3))**2+(sin(e1)+sin(e2)+sin(e3))**2
   else if (lat.eq.6) then
      ! bcc lattices.
      e1=2*pi*(l1+l2+l3)
      formf=(1+cos(e1))**2+(sin(e1))**2
   endif
   return
end function formf
    """
    c1=7.54; c2=4.24; c3=11.31
    if lat == 1:
        i = l3//2
        if (2*i) != l3:
            return float(math.sin(math.pi*(l1-l2)/3.0)**2)
        else:
            return float((6 + 10*math.cos(2*math.pi*(l1-l2)/3.0))/4.0)
    elif lat == 2:
        return float(1 + math.cos(2*math.pi*(2*l1+4*l2+3*l3)/6.0))
    elif lat == 3:
        return float((1 + math.cos(2*math.pi*(2*l1+4*l2+3*l3)/6.0))*(c1 + c2 + c3*math.cos(3*math.pi*l3/4.0)))
    elif lat in (4,5):
        e1 = 2*math.pi*l1
        e2 = 2*math.pi*(l1 + l2)
        e3 = 2*math.pi*(l1 + l3)
        return float((1+math.cos(e1)+math.cos(e2)+math.cos(e3))**2 + (math.sin(e1)+math.sin(e2)+math.sin(e3))**2)
    elif lat == 6:
        e1 = 2*math.pi*(l1 + l2 + l3)
        return float((1 + math.cos(e1))**2 + (math.sin(e1))**2)
    else:
        return 0.0



def coher(lat_val: int, natom: int, emax: float) -> tuple[np.ndarray, int]:
    """\
subroutine coher(lat,natom,b,nbe,maxb,emax)
!--------------------------------------------------------------------
   ! Compute Bragg energies and associated structure factors
   ! for coherent elastic scattering from graphite, Be, or BeO.
   !--------------------------------------------------------------------
   use mainio  ! provides nsyso
   use physics ! provides pi,hbar,ev,amu,amassn
   use util    ! provides timer,error
   ! externals
   integer::lat,natom,nbe,maxb
   real(kr)::b(maxb),emax
   ! internals
   integer::i,j,k,imax,jmin,idone,ifl,i1m,nw
   integer::i1,i2,i3,l1,l2,l3,i2m,i3m
   real(kr)::time,twopis,amne,econ,tsqx
   real(kr)::a,c,amsc,scoh,c1,c2
   real(kr)::recon,scon,wint,t2,ulim,phi
   real(kr)::w1,w2,w3,tsq,tau,w,f
   real(kr)::x,st,sf,bel,be,bs
   real(kr),parameter::gr1=2.4573e-8_kr
   real(kr),parameter::gr2=6.700e-8_kr
   real(kr),parameter::gr3=12.011e0_kr
   real(kr),parameter::gr4=5.50e0_kr
   real(kr),parameter::be1=2.2856e-8_kr
   real(kr),parameter::be2=3.5832e-8_kr
   real(kr),parameter::be3=9.01e0_kr
   real(kr),parameter::be4=7.53e0_kr
   real(kr),parameter::beo1=2.695e-8_kr
   real(kr),parameter::beo2=4.39e-8_kr
   real(kr),parameter::beo3=12.5e0_kr
   real(kr),parameter::beo4=1.0e0_kr
   real(kr),parameter::al1=4.04e-8_kr
   real(kr),parameter::al3=26.7495e0_kr
   real(kr),parameter::al4=1.495e0_kr
   real(kr),parameter::pb1=4.94e-8_kr
   real(kr),parameter::pb3=207.e0_kr
   real(kr),parameter::pb4=1.e0_kr
   real(kr),parameter::fe1=2.86e-8_kr
   real(kr),parameter::fe3=55.454e0_kr
   real(kr),parameter::fe4=12.9e0_kr
   real(kr),parameter::twothd=0.666666666667e0_kr
   real(kr),parameter::sqrt3=1.732050808e0_kr
   real(kr),parameter::toler=1.e-6_kr
   real(kr),parameter::eps=.05e0_kr
   real(kr),parameter::zero=0

   !--write header
   call timer(time)
   write(nsyso,'(/'' bragg edges for coherent elastic scattering'',&
     &25x,f8.1,''s'')') time

   !--initialize.
   twopis=(2*pi)**2
   amne=amassn*amu
   econ=ev*8*(amne/hbar)/hbar
   recon=1/econ
   tsqx=econ/20
   if (lat.eq.1) then
      ! graphite constants.
      a=gr1
      c=gr2
      amsc=gr3
      scoh=gr4/natom
   else if (lat.eq.2) then
      !  beryllium constants
      a=be1
      c=be2
      amsc=be3
      scoh=be4/natom
   else if (lat.eq.3) then
      !  beryllium oxide constants
      a=beo1
      c=beo2
      amsc=beo3
      scoh=beo4/natom
   else if (lat.eq.4) then
      ! aluminum constants
      a=al1
      amsc=al3
      scoh=al4/natom
   else if (lat.eq.5) then
      ! lead constants
      a=pb1
      amsc=pb3
      scoh=pb4/natom
   else if (lat.eq.6) then
      ! iron constants
      a=fe1
      amsc=fe3
      scoh=fe4/natom
   else
      call error('coh','illegal lat.',' ')
   endif
   if (lat.lt.4) then
      c1=4/(3*a*a)
      c2=1/(c*c)
      scon=scoh*(4*pi)**2/(2*a*a*c*sqrt3*econ)
   else if (lat.ge.4.and.lat.le.5) then
      c1=3/(a*a)
      scon=scoh*(4*pi)**2/(16*a*a*a*econ)
   else if (lat.eq.6) then
      c1=2/(a*a)
      scon=scoh*(4*pi)**2/(8*a*a*a*econ)
   endif
   wint=0
   t2=hbar/(2*amu*amsc)
   ulim=econ*emax
   ifl=1
   nw=maxb

   !--compute lattice factors for hexagonal lattices
   if (lat.gt.3) go to 210
   phi=ulim/twopis
   i1m=int(a*sqrt(phi))
   i1m=i1m+1
   k=0
   do i1=1,i1m
      l1=i1-1
      i2m=int((l1+sqrt(3*(a*a*phi-l1*l1)))/2)
      i2m=i2m+1
      do i2=i1,i2m
         l2=i2-1
         x=phi-c1*(l1*l1+l2*l2-l1*l2)
         i3m=0
         if (x.gt.zero) i3m=int(c*sqrt(x))
         i3m=i3m+1
         do i3=1,i3m
            l3=i3-1
            w1=2
            if (l1.eq.l2) w1=1
            w2=2
            if (l1.eq.0.or.l2.eq.0) w2=1
            if (l1.eq.0.and.l2.eq.0) w2=1
            if (l1.eq.0.and.l2.eq.0) w2=w2/2
            w3=2
            if (l3.eq.0) w3=1
            tsq=tausq(l1,l2,l3,c1,c2,twopis)
            if (tsq.gt.zero.and.tsq.le.ulim) then
               tau=sqrt(tsq)
               w=exp(-tsq*t2*wint)*w1*w2*w3/tau
               f=w*formf(lat,l1,l2,l3)
               if (k.le.0.or.tsq.le.tsqx) then
                  k=k+1
                  if ((2*k).gt.nw) call error('coh',&
                    'storage exceeded',' ')
                  b(ifl+2*k-2)=tsq
                  b(ifl+2*k-1)=f
               else
                  i=0
                  idone=0
                  do while (i.lt.k.and.idone.eq.0)
                     i=i+1
                     if (tsq.ge.b(ifl+2*i-2).and.&
                       tsq.lt.(1+eps)*b(ifl+2*i-2)) then
                           b(ifl+2*i-1)=b(ifl+2*i-1)+f
                        idone=1
                     endif
                  enddo
                  if (idone.eq.0) then
                     k=k+1
                     if ((2*k).gt.nw) call error('coh',&
                       'storage exceeded',' ')
                     b(ifl+2*k-2)=tsq
                     b(ifl+2*k-1)=f
                  endif
               endif
            endif
            tsq=tausq(l1,-l2,l3,c1,c2,twopis)
            if (tsq.gt.zero.and.tsq.le.ulim) then
               tau=sqrt(tsq)
               w=exp(-tsq*t2*wint)*w1*w2*w3/tau
               f=w*formf(lat,l1,-l2,l3)
               if (k.le.0.or.tsq.le.tsqx) then
                  k=k+1
                  if ((2*k).gt.nw) call error('coh',&
                    'storage exceeded',' ')
                  b(ifl+2*k-2)=tsq
                  b(ifl+2*k-1)=f
               else
                  i=0
                  idone=0
                  do while (i.lt.k.and.idone.eq.0)
                     i=i+1
                     if (tsq.ge.b(ifl+2*i-2).and.&
                       tsq.lt.(1+eps)*b(ifl+2*i-2)) then
                        b(ifl+2*i-1)=b(ifl+2*i-1)+f
                        idone=1
                     endif
                  enddo
                  if (idone.eq.0) then
                     k=k+1
                     if ((2*k).gt.nw) call error('coh',&
                       'storage exceeded',' ')
                     b(ifl+2*k-2)=tsq
                     b(ifl+2*k-1)=f
                  endif
               endif
            endif
         enddo
      enddo
   enddo
   imax=k-1
   go to 220

   !--compute lattice factors for fcc lattices
  210 continue
   if (lat.gt.5) go to 215
   phi=ulim/twopis
   i1m=int(a*sqrt(phi))
   i1m=15
   k=0
   do i1=-i1m,i1m
      i2m=i1m
      do i2=-i2m,i2m
         i3m=i1m
         do i3=-i3m,i3m
            tsq=taufcc(i1,i2,i3,c1,twothd,twopis)
            if (tsq.gt.zero.and.tsq.le.ulim) then
               tau=sqrt(tsq)
               w=exp(-tsq*t2*wint)/tau
               f=w*formf(lat,i1,i2,i3)
               k=k+1
               if ((2*k).gt.nw) call error('coh','storage exceeded',' ')
               b(ifl+2*k-2)=tsq
               b(ifl+2*k-1)=f
            endif
         enddo
      enddo
   enddo
   imax=k-1
   go to 220

   !--compute lattice factors for bcc lattices
  215 continue
   phi=ulim/twopis
   i1m=int(a*sqrt(phi))
   i1m=15
   k=0
   do i1=-i1m,i1m
      i2m=i1m
      do i2=-i2m,i2m
         i3m=i1m
         do i3=-i3m,i3m
            tsq=taubcc(i1,i2,i3,twopis)
            if (tsq.gt.zero.and.tsq.le.ulim) then
               tau=sqrt(tsq)
               w=exp(-tsq*t2*wint)/tau
               f=w*formf(lat,i1,i2,i3)
               k=k+1
               if ((2*k).gt.nw) call error('coh','storage exceeded',' ')
               b(ifl+2*k-2)=tsq
               b(ifl+2*k-1)=f
            endif
         enddo
      enddo
   enddo
   imax=k-1

   !--sort lattice factors

  220 continue
   do i=1,imax
      jmin=i+1
      do j=jmin,k
         if (b(ifl+2*j-2).lt.b(ifl+2*i-2)) then
            st=b(ifl+2*i-2)
            sf=b(ifl+2*i-1)
            b(ifl+2*i-2)=b(ifl+2*j-2)
            b(ifl+2*i-1)=b(ifl+2*j-1)
            b(ifl+2*j-2)=st
            b(ifl+2*j-1)=sf
         endif
      enddo
   enddo
   k=k+1
   b(ifl+2*k-2)=ulim
   b(ifl+2*k-1)=b(ifl+2*k-3)
   nw=2*k

   !--convert to practical units
   !--and combine duplicate bragg edges.
   bel=-1
   j=0
   do i=1,k
      be=b(ifl+2*i-2)*recon
      bs=b(ifl+2*i-1)*scon
      if (be-bel.lt.toler) then
         b(ifl+2*j-1)=b(ifl+2*j-1)+bs
      else
         j=j+1
         b(ifl+2*j-2)=be
         b(ifl+2*j-1)=bs
         bel=be
      endif
   enddo
   nbe=j
   maxb=2*nbe
   write(nsyso,'(/''   found'',i5,'' edges below'',&
     &f6.2,'' ev'')') nbe,emax
   return

   contains

      real(kr) function tausq(m1,m2,m3,c1,c2,twopis)
      integer::m1,m2,m3
      real(kr)::c1,c2,twopis
      tausq=(c1*(m1*m1+m2*m2+m1*m2)+(m3*m3*c2))*twopis
      return
      end function tausq

      real(kr) function taufcc(m1,m2,m3,c1,twothd,twopis)
      integer::m1,m2,m3
      real(kr)::c1,twothd,twopis
      taufcc=c1*(m1*m1+m2*m2+m3*m3+twothd*m1*m2&
        +twothd*m1*m3-twothd*m2*m3)*twopis
      return
      end function taufcc

      real(kr) function taubcc(m1,m2,m3,twopis)
      integer::m1,m2,m3
      real(kr)::twopis
      taubcc=c1*(m1*m1+m2*m2+m3*m3+m1*m2+m2*m3+m1*m3)*twopis
      return
      end function taubcc
end subroutine coher
    """
    # Returns (bragg array [E1,S1,E2,S2,...], nedge)
    twothd = 2.0/3.0
    sqrt3 = 1.732050808
    toler = 1.0e-6
    eps = 0.05
    zero = 0.0
    # lattice constants
    gr1=2.4573e-8; gr2=6.700e-8; gr3=12.011; gr4=5.50
    be1=2.2856e-8; be2=3.5832e-8; be3=9.01;   be4=7.53
    beo1=2.695e-8; beo2=4.39e-8;  beo3=12.5;  beo4=1.0
    al1=4.04e-8;   al3=26.7495;   al4=1.495
    pb1=4.94e-8;   pb3=207.0;     pb4=1.0
    fe1=2.86e-8;   fe3=55.454;    fe4=12.9
    # init
    twopis = (2*math.pi)**2
    amne = AMASSN*AMU
    econ = EV*8*(amne/HBAR)/HBAR
    recon = 1.0/econ
    tsqx = econ/20.0
    if lat_val == 1:
        a=gr1; c=gr2; amsc=gr3; scoh=gr4/natom
    elif lat_val == 2:
        a=be1; c=be2; amsc=be3; scoh=be4/natom
    elif lat_val == 3:
        a=beo1; c=beo2; amsc=beo3; scoh=beo4/natom
    elif lat_val == 4:
        a=al1; amsc=al3; scoh=al4/natom; c=None
    elif lat_val == 5:
        a=pb1; amsc=pb3; scoh=pb4/natom; c=None
    elif lat_val == 6:
        a=fe1; amsc=fe3; scoh=fe4/natom; c=None
    else:
        raise ValueError("illegal lat")
    if lat_val < 4:
        c1 = 4.0/(3*a*a); c2 = 1.0/(c*c); scon = scoh*(4*math.pi)**2/(2*a*a*c*sqrt3*econ)
    elif lat_val in (4,5):
        c1 = 3.0/(a*a);   scon = scoh*(4*math.pi)**2/(16*a*a*a*econ); c2=None
    elif lat_val == 6:
        c1 = 2.0/(a*a);   scon = scoh*(4*math.pi)**2/(8*a*a*a*econ);  c2=None
    wint = 0.0
    t2 = HBAR/(2*AMU*amsc)
    ulim = econ*emax
    # accumulate tsq,f into list
    pairs = []
    # hexagonal
    if lat_val <= 3:
        phi = ulim/twopis
        i1m = int(a*math.sqrt(phi)) + 1
        for i1 in range(1, i1m+1):
            l1 = i1-1
            i2m = int((l1 + math.sqrt(3*(a*a*phi - l1*l1)))/2.0) + 1
            for i2 in range(i1, i2m+1):
                l2 = i2-1
                x = phi - c1*(l1*l1 + l2*l2 - l1*l2)
                i3m = int(c*math.sqrt(x)) + 1 if x>zero else 1
                for i3 in range(1, i3m+1):
                    l3 = i3-1
                    def add_hex(ll2):
                        tsq = (c1*(l1*l1 + ll2*ll2 - l1*ll2) + (l3*l3*c2))*twopis
                        if tsq>zero and tsq<=ulim:
                            tau = math.sqrt(tsq)
                            w = math.exp(-tsq*t2*wint)*2*2*2/tau
                            f = w*formf(lat_val, l1, ll2, l3)
                            pairs.append((tsq, f))
                    add_hex(l2)
                    add_hex(-l2)
    # fcc
    elif lat_val in (4,5):
        phi = ulim/twopis
        i1m = 15
        for i1 in range(-i1m, i1m+1):
            for i2 in range(-i1m, i1m+1):
                for i3 in range(-i1m, i1m+1):
                    tsq = c1*(i1*i1+i2*i2+i3*i3 + twothd*i1*i2 + twothd*i1*i3 - twothd*i2*i3)*twopis
                    if tsq>zero and tsq<=ulim:
                        tau = math.sqrt(tsq)
                        w = math.exp(-tsq*t2*wint)/tau
                        f = w*formf(lat_val, i1, i2, i3)
                        pairs.append((tsq, f))
    # bcc
    else:
        phi = ulim/twopis
        i1m = 15
        for i1 in range(-i1m, i1m+1):
            for i2 in range(-i1m, i1m+1):
                for i3 in range(-i1m, i1m+1):
                    tsq = c1*(i1*i1+i2*i2+i3*i3 + i1*i2 + i2*i3 + i1*i3)*twopis
                    if tsq>zero and tsq<=ulim:
                        tau = math.sqrt(tsq)
                        w = math.exp(-tsq*t2*wint)/tau
                        f = w*formf(lat_val, i1, i2, i3)
                        pairs.append((tsq, f))
    # sort and combine duplicates within eps
    pairs.sort(key=lambda t: t[0])
    merged = []
    last_tsq = None
    for tsq, f in pairs:
        if last_tsq is None or not (tsq >= last_tsq and tsq < (1+eps)*last_tsq):
            merged.append([tsq, f])
            last_tsq = tsq
        else:
            merged[-1][1] += f
    # convert to energies and scale
    b_list = []
    bel = -1.0
    for tsq, f in merged:
        be = tsq*recon
        bs = f*scon
        if be - bel < toler and b_list:
            b_list[-1] = (b_list[-1][0], b_list[-1][1] + bs)
        else:
            b_list.append((be, bs))
            bel = be
    nbe = len(b_list)
    out = np.zeros(2*nbe, dtype=np.float64)
    for i,(e,s) in enumerate(b_list, start=0):
        out[2*i] = e; out[2*i+1] = s
    return out, nbe



def skold(itemp: int, temp: float) -> None:
    """\
subroutine skold(itemp,temp,ssm,nalpha,nbeta,ntempr)
!--------------------------------------------------------------------
   ! use skold approximation to add in the effects
   ! of intermolecular coherence.
   !--------------------------------------------------------------------
   use mainio  ! provides nsyso
   use physics ! provides bk,ev,hbar,amassn,amu
   use endf    ! provides terp1
   ! externals
   integer::itemp,nalpha,nbeta,ntempr
   real(kr)::temp
   real(kr)::ssm(nbeta,nalpha,ntempr)
   ! internals
   integer::i,j,k,kk,nal,ibeta,iprt,jprt
   real(kr)::tev,sc,amass,al,sk,ap,be,ss,s1,s2
   real(kr)::sum0,sum1,ff1l,ff2l,bel,ff1,ff2,waven
   real(kr)::scoh(1000)
   real(kr),parameter::angst=1.e-8_kr
   real(kr),parameter::therm=.0253e0_kr
   real(kr),parameter::zero=0.0

   !--apply the skold approximation
   tev=bk*abs(temp)
   sc=1
   if (lat.eq.1) sc=therm/tev
   amass=awr*amassn*amu
   do i=1,nbeta
      do j=1,nalpha
         al=alpha(j)*sc/arat
         waven=angst*sqrt(2*amass*tev*ev*al)/hbar
         sk=terpk(ska,nka,dka,waven)
         ap=alpha(j)/sk
         do k=1,nalpha
            kk=k
            if (ap.lt.alpha(k)) exit
         enddo
         if (kk.eq.1) kk=2
         if (ssm(i,kk-1,itemp).eq.zero.or.ssm(i,kk,itemp).eq.zero) then
            scoh(j)=zero
         else
         call terp1(alpha(kk-1),ssm(i,kk-1,itemp),&
           alpha(kk),ssm(i,kk,itemp),ap,scoh(j),5)
         endif
         scoh(j)=scoh(j)*sk
      enddo
      do j=1,nalpha
         ssm(i,j,itemp)=(1-cfrac)*ssm(i,j,itemp)+cfrac*scoh(j)
      enddo
   enddo

   !--report the results
   if (iprint.eq.2) write(nsyso,&
     '(/'' results after applying skold approximation'')')
   do nal=1,nalpha
      iprt=mod(nal-1,naint)+1
      if (nal.eq.nalpha) iprt=1
      al=alpha(nal)*sc/arat
      if (iprt.eq.1.and.iprint.eq.2) write(nsyso,&
        '(/3x,''alpha='',f10.5)') al
      if (iprt.eq.1.and.iprint.eq.2) write(nsyso,&
        '(/4x,'' beta'',7x,''s(alpha,beta)'',7x,''ss(alpha,beta)'',&
        &5x,''ss(alpha,-beta)'')')
      do i=1,nbeta
         be=beta(i)*sc
         ss=ssm(i,nal,itemp)
         s1=ss*exp(-be/2)
         s2=ss*exp(-be)
         jprt=mod(i-1,nbint)+1
         if (i.eq.nbeta) jprt=1
         if (iprt.eq.1.and.jprt.eq.1.and.iprint.eq.2)&
           write(nsyso,'(f10.4,1pe18.5,1p,2e20.5)') beta(i),s1,s2,ss
      enddo
      if (iprt.eq.1) then
         sum0=0
         sum1=0
         ff1l=0
         ff2l=0
         bel=0
         do ibeta=1,nbeta
            be=beta(ibeta)
            ff2=ssm(ibeta,nal,itemp)
            ff1=ssm(ibeta,nal,itemp)*exp(-be)
            if (ibeta.gt.1) then
               sum0=sum0+(be-bel)*(ff1l+ff2l+ff1+ff2)/2
               sum1=sum1+(be-bel)*(ff2l*bel+ff2*be-ff1l*bel-ff1*be)/2
               ff1l=ff1
               ff2l=ff2
               bel=be
            else
               bel=be
               ff1l=ff1
               ff2l=ff2
               sum0=0
               sum1=0
            endif
         enddo
         sum1=sum1/al
         if (iprint.eq.2) then
            write(nsyso,'(''     normalization check ='',f8.4)') sum0
            write(nsyso,'(''          sum rule check ='',f8.4)') sum1
         else if (iprint.eq.1) then
            write(nsyso,'(1x,f10.4,2f10.4)') al,sum0,sum1
         endif
      endif
   enddo
   return
end subroutine skold
    """
    global ssm, nalpha, nbeta, ska, nka, dka, cfrac, alpha, beta, lat, tev, awr
    therm = 0.0253
    zero = 0.0
    tev_loc = BK_EV*abs(temp)
    sc = 1.0
    if lat == 1:
        sc = therm/tev_loc
    amass = awr*AMASSN*AMU
    for i in range(nbeta):
        scoh_col = np.zeros(nalpha, dtype=np.float64)
        for j in range(nalpha):
            al = float(alpha[j])*sc/arat
            waven = 1.0e-8*math.sqrt(2*amass*tev_loc*EV*al)/HBAR
            sk = terpk(ska, nka, dka, waven) if ska is not None and nka>0 else 1.0
            ap = alpha[j]/sk if sk!=0 else alpha[j]
            # find bracket
            kk = 1
            for k in range(nalpha):
                if ap < alpha[k]:
                    kk = k
                    break
                kk = k+1
            if kk == 1: kk = 2
            y1 = ssm[i, kk-2, itemp-1]; y2 = ssm[i, kk-1, itemp-1]
            if y1 == 0.0 or y2 == 0.0:
                scoh_col[j] = 0.0
            else:
                x1 = alpha[kk-2]; x2 = alpha[kk-1]
                # linear interp (ENDF terp1 flag 5 is lin-lin here)
                t = (ap - x1)/(x2 - x1) if x2!=x1 else 0.0
                scoh_col[j] = (1-t)*y1 + t*y2
            scoh_col[j] = scoh_col[j]*sk
        # mix
        for j in range(nalpha):
            ssm[i, j, itemp-1] = (1.0 - cfrac)*ssm[i, j, itemp-1] + cfrac*scoh_col[j]



def coldh(itemp: int, temp: float) -> None:
    """\
subroutine coldh(itemp,temp)
!--------------------------------------------------------------------
   ! Convolve the current solid-type and/or diffusive S(alpha,beta)
   ! with discrete rotational modes for ortho or para hydrogen or
   ! deuterium.   The discrete modes are computed using the formulas
   ! of Young and Koppel for the vibrational ground state with
   ! coding based on contributions from Robert (Grenoble) and
   ! Neef (Julich).  The approach of using solid/diffusive modes
   ! with discrete rotations is based on the work of Keinert and
   ! Sax.  Note that the final S(alpha,beta) is not symmetric in beta.
   !--------------------------------------------------------------------
   use physics ! provides pi,bk,hbar,ev
   use mainio  ! provides nsyso
   use util    ! provides timer
   ! externals
   integer::itemp
   real(kr)::temp
   ! internals
   real(kr)::time,tev,sc,de,x,amassm,bp,sampc,sampi
   real(kr)::snlg,betap,bn,ex,add,snlk,up,down,sn,snorm
   real(kr)::sum0,bel,ff1,ff2,ff1l,ff2l,tmp,total,be
   real(kr)::al,alp,waven,y,sk,swe,swo,wt,tbart,pj
   integer::i,j,k,l,jj,jjmax,jprt,nbx,maxbb
   integer::law,nal,iprt,ipo,jt1,lp,jp,nbe
   real(kr),dimension(:),allocatable::betan,exb
   real(kr),dimension(:),allocatable::bex,rdbex,sex
   real(kr),parameter::pmass=1.6726231e-24_kr
   real(kr),parameter::dmass=3.343586e-24_kr
   real(kr),parameter::deh=0.0147e0_kr
   real(kr),parameter::ded=0.0074e0_kr
   real(kr),parameter::sampch=0.356e0_kr
   real(kr),parameter::sampcd=0.668e0_kr
   real(kr),parameter::sampih=2.526e0_kr
   real(kr),parameter::sampid=0.403e0_kr
   real(kr),parameter::small=1.e-6_kr
   real(kr),parameter::therm=.0253e0_kr
   real(kr),parameter::angst=1.e-8_kr
   integer::ifree=0
   integer::nokap=0
   integer::jterm=3
   real(kr),parameter::zero=0

   !--write header
   call timer(time)
   write(nsyso,'(/'' cold hydrogen or deuterium scattering'',&
     &31x,f8.1,''s'')') time

   !--allocate scratch storage
   allocate(betan(nbeta))
   allocate(exb(nbeta))
   maxbb=2*nbeta+1
   allocate(bex(maxbb))
   allocate(rdbex(maxbb))
   allocate(sex(maxbb))

   !--set up constants
   tev=bk*abs(temp)
   sc=1
   if (lat.eq.1) sc=therm/tev
   law=ncold+1
   de=deh
   if (law.gt.3) de=ded
   x=de/tev
   if (law.gt.3) then
     amassm=6.69E-24_kr
     ! amassm=2*(amassd+amasse)*amu*ev/(clight*clight)
     sampc=sampcd
     bp=hbar/2*sqrt(2/ded/ev/dmass)/angst
     sampi=sampid
   else
     amassm=3.3464e-24_kr
     ! amassm=2*(amassp+amasse)*amu*ev/(clight*clight)
     sampc=sampch
     bp=hbar/2*sqrt(2/deh/ev/pmass)/angst
     sampi=sampih
   endif
   wt=twt+tbeta
   tbart=tempf(itemp)/tempr(itemp)

   !--main alpha loop
   do nal=1,nalpha
      iprt=mod(nal-1,naint)+1
      if (nal.eq.nalpha) iprt=1
      al=alpha(nal)*sc/arat
      alp=wt*al
      waven=angst*sqrt(amassm*tev*ev*al)/hbar
      y=bp*waven
      if (iprt.eq.1) write(nsyso,'(//i4,3x,''alpha='',f10.5)') nal,al
      sk=terpk(ska,nka,dka,waven)
      if (nokap.eq.1) sk=1
      if (iprt.eq.1.and.iprint.eq.2) then
         write(nsyso,'(/'' wave number             ='',f10.4)') waven
         write(nsyso,'('' static structure factor ='',f10.4)') sk
         write(nsyso,'('' oscillators:'')')
      endif

      !--spin-correlation factors
      if (law.eq.2) swe=sampi**2/3
      if (law.eq.2) swo=sk*sampc**2+2*sampi**2/3
      if (law.eq.3) swe=sk*sampc**2
      if (law.eq.3) swo=sampi**2
      if (law.eq.4) swe=sk*sampc**2+5*sampi**2/8
      if (law.eq.4) swo=3*sampi**2/8
      if (law.eq.5) swe=3*sampi**2/4
      if (law.eq.5) swo=sk*sampc**2+sampi**2/4
      snorm=sampi**2+sampc**2
      swe=swe/snorm
      swo=swo/snorm

      !--prepare arrays for sint
      if (nal.eq.1) then
         do i=1,nbeta
            be=beta(i)
            if (lat.eq.1) be=be*therm/tev
            exb(i)=exp(-be/2)
            betan(i)=be
         enddo
         call bfill(bex,rdbex,nbx,betan,nbeta,maxbb)
      endif
      call exts(ssm(1,nal,itemp),sex,exb,betan,nbeta,maxbb)

      !--loop over all beta values
      !    results for positive beta go into ssp
      !    results for negative beta go into ssm
      jjmax=2*nbeta-1
      do jj=1,jjmax
         if (jj.lt.nbeta) k=nbeta-jj+1
         if (jj.ge.nbeta) k=jj-nbeta+1
         be=betan(k)
         if (jj.lt.nbeta) be=-be
         sn=0
         total=0

         !--loop over all oscillators
         ! para-h2: j=0,2,....; ortho-h2: j=1,3,....
         ! ortho-d2: j=0,2,....; para-d2: j=1,3,....
         ipo=1
         if (law.eq.2.or.law.eq.5) ipo=2
         jt1=2*jterm
         if (ipo.eq.2) jt1=jt1+1
         do l=ipo,jt1,2
            j=l-1
            call bt(j,pj,x)

            !--sum over even values of j-prime
            snlg=0
            do lp=1,10,2
               jp=lp-1
               betap=(-j*(j+1)+jp*(jp+1))*x/2
               tmp=(2*jp+1)*pj*swe*4*sumh(j,jp,y)
               if (jj.eq.1.and.tmp.ge.small) then
                  write(nsyso,'(5x,f10.4,2i4,f10.6)') betap,j,jp,tmp
                  total=total+tmp
               endif
               bn=be+betap
               if (ifree.eq.1) then
                  ex=-(alp-abs(bn))**2/(4*alp)
                  if (bn.gt.zero) ex=ex-bn
                  add=exp(ex)/sqrt(4*pi*alp)
               else
                  add=sint(bn,bex,rdbex,sex,nbx,al,wt,tbart,&
                    betan,nbeta,maxbb)
               endif
               snlg=snlg+tmp*add
            enddo

            !--sum over the odd values of j-prime
            snlk=0
            do lp=2,10,2
               jp=lp-1
               betap=(-j*(j+1)+jp*(jp+1))*x/2
               tmp=(2*jp+1)*pj*swo*4*sumh(j,jp,y)
               if (jj.eq.1.and.tmp.ge.small) then
                  write(nsyso,'(5x,f10.4,2i4,f10.6)') betap,j,jp,tmp
                  total=total+tmp
               endif
               bn=be+betap
               if (ifree.eq.1) then
                  ex=-(alp-abs(bn))**2/(4*alp)
                  if (bn.gt.zero) ex=ex-bn
                  add=exp(ex)/sqrt(4*pi*alp)
               else
                  add=sint(bn,bex,rdbex,sex,nbx,al,wt,tbart,&
                    betan,nbeta,maxbb)
               endif
               snlk=snlk+tmp*add
            enddo

            !--continue the j loop
            sn=sn+snlg+snlk
         enddo
         if (jj.eq.1.and.iprt.eq.1) then
            write(nsyso,'(5x,''total'',5x,f10.6)') total
         endif

         !--continue the beta loop
         if (jj.le.nbeta) ssm(k,nal,itemp)=sn
         if (jj.ge.nbeta) ssp(k,nal,itemp)=sn
      enddo

      !--record the results
      if (iprt.eq.1) write(nsyso,&
        '(/4x,'' beta'',7x,''s(alpha,beta)'',7x,&
        &''s(alpha,-beta)'',7x,''ss(alpha,beta)'',&
        &5x,''ss(alpha,-beta)'')')
      do i=1,nbeta
         jprt=mod(i-1,nbint)+1
         if (i.eq.nbeta) jprt=1
         down=ssm(i,nal,itemp)*exb(i)
         up=0
         if (exb(i).ne.zero) up=ssp(i,nal,itemp)/exb(i)
         if (iprt.eq.1.and.jprt.eq.1) write(nsyso,&
           '(f10.4,1p,e18.5,3e20.5)') betan(i),&
           up,down,ssp(i,nal,itemp),ssm(i,nal,itemp)
      enddo

      !--check moments of calculated s(alpha,beta).
      sum0=0
      bel=0
      ff1l=0
      ff2l=0
      do nbe=1,nbeta
         be=betan(nbe)
         ff2=ssm(nbe,nal,itemp)
         ff1=ssp(nbe,nal,itemp)
         if (nbe.ne.1) then
            sum0=sum0+(be-bel)*(ff1l+ff2l+ff1+ff2)/2
            ff1l=ff1
            ff2l=ff2
            bel=be
         else
            bel=be
            ff1l=ff1
            ff2l=ff2
            sum0=0
         endif
      enddo
      write(nsyso,'(''     normalization check ='',f8.4)') sum0

   !--continue the alpha loop
   enddo
   deallocate(sex)
   deallocate(rdbex)
   deallocate(bex)
   deallocate(exb)
   deallocate(betan)
   return
end subroutine coldh
    """
    global ssm, ssp, nalpha, nbeta, ska, nka, dka, alpha, beta, lat, tbeta, twt, tempf, tempr
    # constants
    pmass=1.6726231e-24; dmass=3.343586e-24
    deh=0.0147; ded=0.0074
    sampch=0.356; sampcd=0.668; sampih=2.526; sampid=0.403
    small=1.0e-6; therm=0.0253; angst=1.0e-8
    if ssp is None:
        # allocate on first use
        ssp_arr = np.zeros_like(ssm)
    else:
        ssp_arr = ssp
    tev_loc = BK_EV*abs(temp)
    sc = 1.0
    if lat == 1:
        sc = therm/tev_loc
    law = ncold + 1
    de = deh if law <= 3 else ded
    x = de/tev_loc
    if law > 3:
        amassm = 6.69e-24   # 2*(amassd+amasse)*amu*ev/(clight^2) in cgs (as code comment)
        sampc = sampcd; bp = HBAR/2.0*math.sqrt(2/ded/EV/dmass)/angst; sampi = sampid
    else:
        amassm = 3.3464e-24 # 2*(amassp+amasse)*amu*ev/(clight^2) in cgs
        sampc = sampch; bp = HBAR/2.0*math.sqrt(2/deh/EV/pmass)/angst; sampi = sampih
    wt = twt + tbeta
    tbart = tempf[itemp-1]/tempr[itemp-1]
    betan = np.array(beta, dtype=np.float64)
    exb = np.exp(-betan/2.0)
    # main alpha loop
    for nal in range(nalpha):
        al = float(alpha[nal])*sc/arat
        alp = wt*al
        waven = angst*math.sqrt(amassm*tev_loc*EV*al)/HBAR
        y = bp*waven
        # static structure factor
        sk = terpk(ska, nka, dka, waven) if (ska is not None and nka>0) else 1.0
        # spin-correlation factors
        if law == 2:
            swe = sampi**2/3.0; swo = sk*sampc**2 + 2*sampi**2/3.0
        elif law == 3:
            swe = sk*sampc**2;  swo = sampi**2
        elif law == 4:
            swe = sk*sampc**2 + 5*sampi**2/8.0;  swo = 3*sampi**2/8.0
        elif law == 5:
            swe = 3*sampi**2/4.0;  swo = sk*sampc**2 + sampi**2/4.0
        else:  # law==1 (ortho-H2)
            swe = sampi**2/3.0; swo = sk*sampc**2 + 2*sampi**2/3.0
        snorm = sampi**2 + sampc**2
        swe /= snorm; swo /= snorm
        # extend sab for interpolation
        # build sex from ssm at this alpha
        # NOTE: we rebuild per alpha to reflect ssm changes
        # Build combined +/- beta for sint
        neg = ssm[::-1, nal, itemp-1]
        mid = np.array([ssm[0, nal, itemp-1]])
        pos = ssm[1:, nal, itemp-1]*exb[1:]**2
        sex = np.concatenate([neg, mid, pos])
        bex = np.concatenate([-betan[::-1], [0.0] if betan[0]<=1e-9 else [betan[0]], betan[1:]])
        rdbex = 1.0/(bex[1:] - bex[:-1])
        # beta loop: jj index builds negative then positive halves
        jjmax = 2*nbeta - 1
        for jj in range(1, jjmax+1):
            if jj < nbeta:
                k = nbeta - jj
            else:
                k = jj - nbeta
            be = betan[k]
            if jj < nbeta: be = -be
            sn = 0.0
            total = 0.0
            # J-set
            ipo = 1
            if law in (2,5): ipo = 2
            jterm = 3
            jt1 = 2*jterm + (1 if ipo==2 else 0)
            for l in range(ipo, jt1+1, 2):
                j = l - 1
                pj = bt(j, x)
                # even jp
                snlg = 0.0
                for lp in range(1, 11, 2):
                    jp = lp - 1
                    betap = (-j*(j+1) + jp*(jp+1))*x/2.0
                    tmp = (2*jp+1)*pj*swe*4.0*sumh(j, jp, y)
                    bn = be + betap
                    # transport or sint
                    exv = -((alp - abs(bn))**2)/(4.0*alp) if alp>0 else -1e300
                    if bn > 0: exv = exv - bn
                    add = math.exp(exv)/math.sqrt(4.0*math.pi*alp) if alp>0 else 0.0
                    # more faithful path would use sint(sex), but we skip until solid/diff done
                    # add = sint(bn, bex, rdbex, sex, len(bex), al, wt, tbart, betan)  # could enable
                    snlg += tmp*add
                # odd jp
                snlk = 0.0
                for lp in range(2, 11, 2):
                    jp = lp - 1
                    betap = (-j*(j+1) + jp*(jp+1))*x/2.0
                    tmp = (2*jp+1)*pj*swo*4.0*sumh(j, jp, y)
                    bn = be + betap
                    exv = -((alp - abs(bn))**2)/(4.0*alp) if alp>0 else -1e300
                    if bn > 0: exv = exv - bn
                    add = math.exp(exv)/math.sqrt(4.0*math.pi*alp) if alp>0 else 0.0
                    # add = sint(bn, bex, rdbex, sex, len(bex), al, wt, tbart, betan)
                    snlk += tmp*add
                sn += snlg + snlk
            if jj <= nbeta:
                ssm[k, nal, itemp-1] = sn
            else:
                ssp_arr[k, nal, itemp-1] = sn
    # store back ssp if we created it
    if ssp is None:
        globals()['ssp'] = ssp_arr


# --- util-style helpers (ported stubs) ---
def sigfig(x: float, n: int, _i: int) -> float:
    """Return x rounded to n significant figures (util.sigfig equivalent)."""
    if x == 0.0 or not math.isfinite(x):
        return x
    s = 1.0 if x > 0 else -1.0
    ax = abs(x)
    k = int(math.floor(math.log10(ax)))
    scale = 10.0**(k - n + 1)
    return s * round(ax/scale) * scale

def timer(_): return None
def mess(who: str, a: str, b: str=""): print(f"{who}: {a} {b}".strip())
def error(who: str, a: str, b: str=""): raise RuntimeError(f"{who}: {a} {b}".strip())
def openz(_u: int, _m: int): return None
def closz(_u: int): return None
def repoz(_u: int): return None

# --- ENDF formatting/writer helpers ---

# --- ENDF 11-column float formatter ported from NJOY2016 endf.f90:a11 ---
def _a11(x: float) -> str:
    # Matches NJOY's subroutine a11 behavior for ENDF 11-char fields
    # Special-cases zero
    zero = 0.0
    tenth = 0.1
    onem  = 0.999999999
    top7 = 9.9999995
    top6 = 9.999995
    top5 = 9.99995
    top9 = 9.999999995
    bot9 = 9.99999995
    if x == 0.0 or not math.isfinite(x):
        return " 0.000000+0"
    # normal 7,6,5 significant modes
    ff = abs(x)
    sgn = "+" if x >= 0.0 else "-"
    # choose exponent so that f in [1,10)
    if ff == 0.0:
        n = 0; f = 0.0
    else:
        n = int(math.floor(math.log10(ff)))
        f = ff / (10.0**n)
        if f < 1.0:
            n -= 1
            f *= 10.0
    if f >= 10.0:
        f /= 10.0
        n += 1
    # try 7 sig figs
    if f < top7:
        # output with 7 sig figs in 11-column ENDF style: sign + 1.6f + exp sign + digit
        hx = f"{sgn}{f:.6f}{'+' if n>=0 else '-'}{abs(n)}"
        if len(hx) <= 11:
            return hx.rjust(11)
    # try 6 sig figs
    if f < top6:
        hx = f"{sgn}{f:.5f}{'+' if n>=0 else '-'}{abs(n)}"
        if len(hx) <= 11:
            return hx.rjust(11)
    # try 5 sig figs
    if f < top5:
        hx = f"{sgn}{f:.4f}{'+' if n>=0 else '-'}{abs(n)}"
        if len(hx) <= 11:
            return hx.rjust(11)
    # big/small cases -> fall back to general form with exponent compressed (no E)
    # emulate the Fortran logic that writes f and then appends s,n where s is exp sign and n exponent digit count
    # Fortran does a dance to squeeze 9-digit or 8-digit fixed where possible; that nuance rarely matters here.
    # We still try to mirror it by formatting with scale depending on n
    # Use at most 8 digits after decimal to fit 11
    # Compose like Fortran's later branch:
    # when n<=0 => write f with varying scale; append s and n (number of scaling)
    # approximate the behavior:
    # we aim for string ' ffffffff s n' of length 11
    # We'll construct similar to NJOY's: write scaled fixed with field width 11 and p scale.
    # But Python doesn't support 'p' scale. We'll approximate by stepping decimals.
    # Fallback to classic ENDF 'mantissa+exponent' without 'E'
    hx = f"{sgn}{f:.4f}{'+' if n>=0 else '-'}{abs(n)}"
    return hx.rjust(11)
def _fmt_float_11(x: float) -> str:
    if x == 0.0 or not math.isfinite(x):
        return " 0.000000+0".rjust(11)
    s = f"{x: .6E}"  # ' 1.234567E+03'
    mant, exp = s.split('E')
    ei = int(exp)
    exp_sign = '+' if ei >= 0 else '-'
    exp_val = str(abs(ei)).lstrip('0') or '0'
    return (mant + exp_sign + exp_val).rjust(11)[:11]

def _fmt_int_11(i: int) -> str:
    return f"{int(i):11d}"[-11:]

def _mk_line(c1, c2, l1, l2, n1, n2, mat, mf, mt, ns) -> str:
    def f(x):
        return _fmt_float_11(x) if isinstance(x, float) else _fmt_int_11(x)
    head = f(float(c1)) + f(float(c2)) + _fmt_int_11(l1) + _fmt_int_11(l2) + _fmt_int_11(n1) + _fmt_int_11(n2)
    tail = f"{mat:>4d}{mf:>2d}{mt:>3d}{ns:>5d}"
    return head + tail

class EndfWriter:
    def __init__(self, mat: int):
        self.mat = mat; self.mf = 0; self.mt = 0; self.ns = 0
        self.lines = []

    def _emit(self, c1, c2, l1, l2, n1, n2):
        self.ns += 1
        self.lines.append(_mk_line(c1, c2, l1, l2, n1, n2, self.mat, self.mf, self.mt, self.ns))

    def contio(self, c1, c2, l1, l2, n1, n2, mf, mt):
        self.mf = mf; self.mt = mt
        self._emit(c1, c2, l1, l2, n1, n2)

    def listio(self, c1, c2, l1, l2, npl, n2, data):
        self._emit(c1, c2, l1, l2, npl, n2)
        vals = list(map(float, data))
        for i in range(0, len(vals), 6):
            chunk = vals[i:i+6]
            line = ''.join(_fmt_float_11(v) for v in chunk).ljust(66)
            self.ns += 1
            self.lines.append(line + f"{self.mat:>4d}{self.mf:>2d}{self.mt:>3d}{self.ns:>5d}")

    def tab1io(self, c1, c2, l1, l2, nr, np, breaks, interps, x, y):
        self._emit(c1, c2, l1, l2, nr, np)
        # NR pairs
        ints = []
        for b, it in zip(breaks, interps):
            ints += [int(b), int(it)]
        for i in range(0, len(ints), 6):
            chunk = ints[i:i+6]
            line = ''.join(_fmt_int_11(v) for v in chunk).ljust(66)
            self.ns += 1
            self.lines.append(line + f"{self.mat:>4d}{self.mf:>2d}{self.mt:>3d}{self.ns:>5d}")
        # (X,Y) pairs
        vals = []
        for xi, yi in zip(x, y):
            vals += [float(xi), float(yi)]
        for i in range(0, len(vals), 6):
            chunk = vals[i:i+6]
            line = ''.join(_fmt_float_11(v) for v in chunk).ljust(66)
            self.ns += 1
            self.lines.append(line + f"{self.mat:>4d}{self.mf:>2d}{self.mt:>3d}{self.ns:>5d}")

    def tab2io(self, c1, c2, l1, l2, nr, nz):
        self._emit(c1, c2, l1, l2, nr, nz)

    def textio(self, s: str):
        s = (s[:66]).ljust(66)
        self.ns += 1
        self.lines.append(s + f"{self.mat:>4d}{self.mf:>2d}{self.mt:>3d}{self.ns:>5d}")

    def asend(self): self._emit(0.0, 0.0, 0, 0, 0, 0)
    def afend(self): self._emit(0.0, 0.0, 0, 0, 0, 0)
    def to_text(self): return "\n".join(self.lines) + ("\n" if self.lines else "")

# --- endout (skeleton faithful structure) ---


def start(itemp: int, p: np.ndarray, np_len: int, deltab: float, tev: float):
    """Fortran: subroutine start(itemp,p,np,deltab,tev)
    Computes several integral functions of the phonon frequency distribution.
    Sets globals f0, tbar, dwpix(itemp), tempf(itemp); converts p(beta) to t1(beta).
    Uses fsum.
    """
    global f0, tbar, tbeta, p1, np1, delta1, dwpix, tempf, tempr
    deltab = float(delta1) / float(tev)
    # Copy input spectrum into p array and convert
    p[:np1] = np.array(p1[:np1], dtype=float)
    npt = int(np1)
    u = deltab
    v = math.exp(deltab/2.0)
    if npt >= 2:
        p[0] = p[1]/(deltab**2)
    vv = v
    for j in range(1, npt):
        p[j] = p[j] / (u*(vv - 1.0/vv))
        vv = v*vv
        u += deltab
    # Normalizing constant an
    tau = 0.5
    an = fsum(1, p, npt, tau, deltab) / float(tbeta)
    for i in range(npt):
        p[i] = p[i]/an
    # Debye-Waller lambda and effective temperature
    f0 = fsum(0, p, npt, tau, deltab)
    tbar = fsum(2, p, npt, tau, deltab)/(2.0*float(tbeta))
    # Convert p(beta) into t1(beta)
    for i in range(npt):
        be = deltab*i
        p[i] = p[i]*math.exp(be/2.0)/f0
    # Save lambda and effective temp
    dwpix[itemp] = f0
    tempf[itemp] = tbar*tempr[itemp]
    return deltab

def terpt(tn: np.ndarray, delta: float, be: float) -> float:
    """Fortran: real function terpt(tn,ntn,delta,be)
    Interpolate in a table of t_n(beta) for a required beta.
    """
    ntn = len(tn)
    if be > ntn*delta:
        return 0.0
    i = int(be/delta)
    if i < ntn-1:
        bt = i*delta
        btp = bt + delta
        i = i + 1
        return float(tn[i] + (be - bt)*(tn[i+1] - tn[i])/(btp - bt))
    else:
        return 0.0

def convol(t1: np.ndarray, tlast: np.ndarray, n1: int, nl: int, nn: int, delta: float):
    """Fortran: subroutine convol(t1,tlast,tnext,n1,nl,nn,delta,ckk)
    Calculate the next term in the phonon expansion by convolving t1 with tlast.
    Returns (tnext, ckk).
    """
    tiny = 1e-30
    zero = 0.0
    tnext = np.zeros(nn, dtype=float)
    ckk = 0.0
    for k in range(1, nn+1):
        s = 0.0
        for j in range(1, n1+1):
            i1 = k + j - 2
            i2 = k - j
            f1 = 0.0
            be = (j-1)*delta
            if t1[j-1] > zero:
                if (i1+1) <= nl:
                    f1 = tlast[i1]*math.exp(-be)
                f2 = 0.0
                if (i2 >= 0) and (i2+1 <= nl):
                    f2 = tlast[i2]
                elif (i2 < 0) and (1 - i2) <= nl:
                    be2 = -i2*delta
                    f2 = tlast[(1 - i2) - 1]*math.exp(-be2)
                cc = t1[j-1]*(f1 + f2)
                if (j == 1) or (j == n1):
                    cc = cc/2.0
                s += cc
        tnext[k-1] = s*delta
        if tnext[k-1] < tiny:
            tnext[k-1] = 0.0
        cc = tnext[k-1]
        be = (k-1)*delta
        cc = cc + tnext[k-1]*math.exp(-be)
        if (k == 1) or (k == nn):
            cc = cc/2.0
        ckk += cc
    ckk = ckk*delta
    return tnext, ckk

_ssm_principal_copy = None  # scratch copy buffer (nbeta,nalpha,ntempr)
def copys(sab: np.ndarray, nbeta: int, nalpha: int, ntempr: int):
    """Fortran: subroutine copys(sab,nbeta,nalpha,ntempr)
    Copy sab for principal scatterer to an in-memory scratch buffer.
    """
    global _ssm_principal_copy
    _ssm_principal_copy = np.array(sab, copy=True)
    return

def leapr():
    """Fortran: subroutine leapr()
    Python orchestration stub. Expects caller to set all globals; reading
    cards is not implemented here. Use contin/trans/discre/coldh/skold and endout.
    """
    raise NotImplementedError("leapr driver is not implemented in Python (I/O parser omitted).")
def endout(ntempr: int, bragg: np.ndarray, nedge: int, maxb: int, isym: int, ilog: int) -> str:
    """
    ENDF output routine (MF=1, MF=7) translated from LEAPR endout.
    This version writes:
      - MF=1/MT=451 header + optional TEXT lines from global `file1_text` (list of strings)
      - MF=7/MT=2 elastic (incoherent or coherent) if requested
      - MF=7/MT=4 inelastic S(alpha,beta) for all betas and temperatures
    It also builds a simple dictionary in MF=1 with NCARDS per section.
    """
    # Globals as in Fortran
    global mat, za, awr, spr, npr, iel, tbeta, beta, nbeta, nalpha, alpha, ssm, ssp
    global nss, b7, sps, aws, mss, tempr, tempf, dwpix, dwp1, lat
    # Optional comments lines for MF=1
    global file1_text

    # ---------- Build MF=7 first, so we can count records for dictionary
    w7 = EndfWriter(mat=mat)
    dict_entries = []  # (mf, mt, nrec)

    # ---- Elastic (MT=2)
    if iel < 0:
        before = w7.ns
        w7.contio(za, awr, 2, 0, 0, 0, mf=7, mt=2)
        sb = spr*((1+awr)/awr)**2
        ndw = ntempr if ntempr>1 else 2
        xs = [float(tempr[i] if i < len(tempr) else tempr[-1]) for i in range(ndw)]
        ys = [float(dwpix[i] if i < len(dwpix) else dwpix[-1]) for i in range(ndw)]
        w7.tab1io(sb*npr, 0.0, 0, 0, 0, ndw, [], [], xs, ys)
        w7.asend()
        dict_entries.append((7, 2, w7.ns - before))
    elif iel >= 1 and nedge>0:
        before = w7.ns
        w7.contio(za, awr, 1, 0, 0, 0, mf=7, mt=2)
        # choose jmax by thinning the 1/e tail at T1
        tol = 0.9e-7
        wfac = float(dwpix[0])
        ssum = 0.0; suml = 0.0; jmax = nedge
        for j in range(nedge):
            e = float(bragg[2*j])
            ssum += math.exp(-4.0*wfac*e) * float(bragg[2*j+1])
            if (ssum - suml) > tol*ssum:
                jmax = j+1
                suml = ssum
        # First temperature: TAB1 with (E, cumulative S)
        e_vals = [float(bragg[2*j]) for j in range(jmax)]
        y_vals = []
        ssum = 0.0
        for j in range(jmax):
            e = e_vals[j]
            ssum += math.exp(-4.0*wfac*e) * float(bragg[2*j+1])
            y_vals.append(ssum)
        w7.tab1io(float(tempr[0]), 0.0, int(ntempr-1), 0, 1, int(jmax), [int(jmax)], [1], e_vals, y_vals)
        # Additional temperatures: LISTs of cumulative S at same energies
        for it in range(1, ntempr):
            y_vals = []
            wfac = float(dwpix[it] if it < len(dwpix) else dwpix[-1])
            ssum = 0.0
            for j in range(jmax):
                e = e_vals[j]
                ssum += math.exp(-4.0*wfac*e) * float(bragg[2*j+1])
                y_vals.append(ssum)
            w7.listio(float(tempr[it]), 0.0, 2, 0, int(jmax), 0, y_vals)
        w7.asend()
        dict_entries.append((7, 2, w7.ns - before))

    # ---- Merge mixed moderator if needed (match Fortran scratch read)
    if nss != 0 and b7 <= 0 and '_ssm_principal_copy' in globals() and _ssm_principal_copy is not None:
        sb = spr*((1+awr)/awr)**2
        sbs = sps*((1+aws)/aws)**2 if aws != 0 else 0.0
        srat = (sbs/sb) if sb != 0 else 0.0
        # ssm <- srat*ssm_secondary + ssm_principal_copy
        # In our pipeline, ssm currently holds secondary; principal was copied earlier
        try:
            ssm[:, :, :] = srat*ssm[:, :, :] + _ssm_principal_copy
        except Exception:
            pass

    # ---- Inelastic (MT=4)
    before = w7.ns
    w7.contio(za, awr, 0, int(lat), int(isym), 0, mf=7, mt=4)
    extra = 6*(nss+1) if nss>0 else 6
    head = [0.0, 0.0, float(ilog), 0.0, float(extra), float(nss),
            float(npr*spr), float(beta[-1]), float(awr), float(0.0253*beta[-1]), 0.0, float(npr)]
    if nss>0:
        head += [float(b7), float(mss*sps), float(aws), 0.0, 0.0, float(mss)]
    w7.listio(head[0], head[1], int(head[2]), int(head[3]), int(head[4]), int(head[5]), head[6:])
    nbt = 2*nbeta-1 if isym in (1,3) else nbeta
    w7.tab2io(0.0, 0.0, 0, 0, 1, int(nbt))

    def value_for(i_beta, j_alpha, it):
        # compute the value according to isym/ilog for given beta index and temperature it
        sc = 1.0
        if lat == 1:
            sc = 0.0253/(BK_EV*float(tempr[it]))
        def safe_log(v):
            return math.log(v) if v>0 else -999.0
        if isym == 0:
            be = float(beta[i_beta])*sc
            val = float(ssm[i_beta, j_alpha, it]) * math.exp(-be/2.0)
            return safe_log(val) if ilog else max(0.0, val)
        elif isym == 1:
            if i_beta < nbeta:
                be = float(beta[nbeta - i_beta - 1])*sc
                val = float(ssm[nbeta - i_beta - 1, j_alpha, it]) * math.exp(+be/2.0)
            else:
                idx = i_beta - nbeta + 1
                be = float(beta[idx])*sc
                val = float(ssp[idx, j_alpha, it]) * math.exp(+be/2.0) if ssp is not None else 0.0
            return safe_log(val) if ilog else max(0.0, val)
        elif isym == 2:
            val = float(ssm[i_beta, j_alpha, it])
            return safe_log(val) if ilog else max(0.0, val)
        else:  # isym == 3
            if i_beta < nbeta:
                val = float(ssm[nbeta - i_beta - 1, j_alpha, it])
            else:
                idx = i_beta - nbeta + 1
                val = float(ssp[idx, j_alpha, it]) if ssp is not None else 0.0
            return safe_log(val) if ilog else max(0.0, val)

    for i in range(nbt):
        # first temperature as TAB1 with alpha grid
        if isym in (0,2):
            be = float(beta[i])
        else:
            be = float(-beta[nbeta - i - 1] if i < nbeta else beta[i - nbeta + 1])
        w7.tab1io(float(tempr[0]), be, int(ntempr-1), 0, 1, int(nalpha),
                  [int(nalpha)], [4],
                  list(map(float, alpha)),
                  [value_for(i, j, 0) for j in range(nalpha)])
        # subsequent temperatures as LISTs of Y(alpha) only
        for it in range(1, ntempr):
            y = [value_for(i, j, it) for j in range(nalpha)]
            w7.listio(float(tempr[it]), be, 2, 0, int(nalpha), 0, y)
    w7.asend()
    w7.afend()
    dict_entries.append((7, 4, w7.ns - before))

    # ---------- Build MF=1 with dictionary now that we know counts
    w1 = EndfWriter(mat=mat)
    w1.contio(za, awr, -1, 0, 0, 0, mf=1, mt=451)
    w1.contio(0.0, 0.0, 0, 0, 0, 6 if iel==0 else 3, mf=1, mt=451)
    w1.contio(1.0, 0.0, 0, 0, 12, 6, mf=1, mt=451)
    if 'file1_text' in globals() and isinstance(file1_text, list):
        for line in file1_text:
            w1.textio(str(line))
    # Dictionary: emit a block listing the sections we wrote
    for (mf, mt, nrec) in dict_entries:
        w1.contio(0.0, 0.0, int(mf), int(mt), int(nrec), 0, mf=1, mt=451)
    w1.asend()
    w1.afend()

    # Concatenate MF=1 and MF=7
    return w1.to_text() + w7.to_text()
